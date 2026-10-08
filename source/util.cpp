/*
*****************************************************************************
UTILITY
---------------------


*****************************************************************************
*/

#include "lisflood.h"
#include <sys/stat.h>

//-----------------------------------------------------------------------------------
// CHECK FOR DRYING ELEMENTS
// If dV is going to make the water depth -ve, scale all flows out of cell
// so that the water volume in the cell goes to 0 in
void DryCheck(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
	int i, j;
	NUMERIC_TYPE WDweight, dV;
	NUMERIC_TYPE *hptr, *q1, *q2, *q3, *q4, cv;

	// Check for drying elements
	//#pragma omp parallel for private( i,hptr,q1,q2,q3,q4,cv,WDweight,dV)
	for (j = 0; j < Parptr->ysz; j++)
	{
		hptr = Arrptr->H + j*Parptr->xsz;
		for (i = 0; i<Parptr->xsz; i++, hptr++) if (*hptr>Solverptr->DepthThresh)
		{

			//if (j == 0 && i == 80) {

			//	bool test;
			//	test = 1;
			//}

			q1 = Arrptr->Qx + j*(Parptr->xsz + 1) + i;
			q2 = Arrptr->Qx + j*(Parptr->xsz + 1) + i + 1;
			q3 = Arrptr->Qy + j*(Parptr->xsz + 1) + i;
			q4 = Arrptr->Qy + (j + 1)*(Parptr->xsz + 1) + i;

			dV = Solverptr->Tstep*(*q1 - *q2 + *q3 - *q4);
			cv = *hptr*Parptr->dA;

			if (cv + dV < 0)
			{
				WDweight = -cv / dV; // C(0.5) for improved drying stability
				if (*q1 < 0) *q1 *= WDweight;
				if (*q2 > 0) *q2 *= WDweight;
				if (*q3 < 0) *q3 *= WDweight;
				if (*q4 > 0) *q4 *= WDweight;
			}
		}
	}

	return;
}


//---------------------------------------------------------------------------
// GENERAL UTILITY BOBBINS
NUMERIC_TYPE x_centre(Pars *Parptr, const int i)
{
	return Parptr->blx + C(0.5)*Parptr->dx + Parptr->dx*i;
}

NUMERIC_TYPE y_centre(Pars *Parptr, const int j)
{
	return Parptr->tly - C(0.5)*Parptr->dy - Parptr->dy*j;
}

NUMERIC_TYPE x_vertex(Pars *Parptr, const int i)
{
	return Parptr->blx + Parptr->dx*i;
}

NUMERIC_TYPE y_vertex(Pars *Parptr, const int j)
{
	return Parptr->tly - Parptr->dy*j;
}

//---------------------------------------------------------------------------
// CALCULATE VOLUME OF WATER IN CHANNEL AND FLOODPLAIN
NUMERIC_TYPE DomainVol(States *Statesptr, Pars *Parptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr)
{
	NUMERIC_TYPE vol = C(0.0);
	int i, j, pi, pj, chseg, pH0, p0;
	NUMERIC_TYPE HeightOverBank, por0, dAPor;
	ChannelSegmentType *csp;

#pragma omp parallel for private( i,pH0,por0,dAPor,p0) reduction( + : vol )
	for (j = 0; j < Parptr->ysz; j++)
		for (i = 0; i<Parptr->xsz; i++)
		{
			p0 = i + j*Parptr->xsz;
			if (Statesptr->porosity == ON)
			{
				if (Parptr->Por_Ident == 1 || Parptr->Por_Ident == 3)
				{
					por0 = Arrptr->paerial[p0];
					dAPor = Parptr->dA*por0;
					if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * dAPor;
				}
				else if (Parptr->Por_Ident == 2 || Parptr->Por_Ident == 4)
				{
					pH0 = (int)(Arrptr->H[p0] / Parptr->zlev);
					if (pH0 > (Parptr->maxelev / Parptr->zlev)) pH0 = (int)(Parptr->maxelev / Parptr->zlev);
					por0 = Arrptr->paerial[i + j*Parptr->xsz + pH0*Parptr->xsz*Parptr->ysz];
					dAPor = Parptr->dA*por0;
					if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * dAPor;
				}
			}
			else if (Statesptr->SGC == ON)
			{
				// calculate add on volume
				vol += Arrptr->SGCVol[p0];
			}
			else
			{
				if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * Parptr->dA;
			}
		}

	// Calculate Channel Volume if present
	if (Statesptr->ChannelPresent == ON)
	{
		for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++)
		{
			csp = ChannelSegments + chseg;
			for (i = 0; i < csp->chsz; i++)
			{
				pi = csp->ChanX[i];
				pj = csp->ChanY[i];
				vol += Arrptr->H[pi + pj*Parptr->xsz] * csp->ChanWidth[i] * csp->Chandx[i];

				HeightOverBank = Arrptr->DEM[pi + pj*Parptr->xsz] + Arrptr->H[pi + pj*Parptr->xsz] - csp->BankZ[i];

				if (Statesptr->NCFS && HeightOverBank>0 && csp->ChanWidth[i] < Parptr->dx)
					vol += csp->Chandx[i] * (Parptr->dx - csp->ChanWidth[i])*HeightOverBank;
			}
		}
	}

	return(vol);
}
//---------------------------------------------------------------------------
// REMOVE LOW LYING BANK AREAS
// Searches through floodplain cells neighbouring the channel - if it finds
// one less than htol above bed elevation it is raised to htol above bed
// elevation. Low lying bank pixels can lead to instabilities.
void SmoothBanks(Pars *Parptr, Solver *Solverptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, const int verbose)
{
	int i, pi, pj, di, dj;
	NUMERIC_TYPE zchan;
	int chseg;
	ChannelSegmentType *csp;

	if (verbose == ON) printf("Smoothing bank cells with tolerance:\t%" NUM_FMT"m\t\n", Solverptr->htol);

	for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
	{
		csp = ChannelSegments + chseg;

		for (i = 0; i < csp->chsz - 1; i++) {
			pi = csp->ChanX[i];
			pj = csp->ChanY[i];
			zchan = Arrptr->DEM[pi + pj*Parptr->xsz];

			for (di = -1; di <= 1; di++){
				for (dj = -1; dj <= 1; dj++){
					if (!(di == 0 && dj == 0)) {
						if (pi + di >= 0 && pi + di < Parptr->xsz && pj + dj >= 0 && pj + dj < Parptr->ysz){
							if (Arrptr->ChanMask[pi + di + (pj + dj)*Parptr->xsz] == -1){
								if (Arrptr->DEM[pi + di + (pj + dj)*Parptr->xsz] < (zchan + Solverptr->htol)){
									Arrptr->DEM[pi + di + (pj + dj)*Parptr->xsz] = zchan + Solverptr->htol;
								}
							}
						}
					}
				}
			}
		}
	}

	if (verbose == ON) printf("Smoothing bank cells done.\n\n");

	return;
}
//-----------------------------------------------------------------------------------
// (MT) check if file exists, return 0 if it doesn't and 1 if it does
int fexist(char *filename) {
	struct stat buffer;
	if (stat(filename, &buffer)) return 0;
	return 1;
}
//-----------------------------------------------------------------------------------
// Calculates velocity and hazard
void UpdateV(States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr)
{
	int i, j, p0, pxy0, px1, pyl;
	NUMERIC_TYPE Vc, Xv, Yv, Haz;

#pragma omp parallel for private( i,p0,pxy0,px1,pyl,Xv,Yv,Vc,Haz)
	// calc at Ps
	for (j = 0; j < Parptr->ysz; j++)
	{
		for (i = 0; i<Parptr->xsz; i++)
		{
			p0 = i + j*Parptr->xsz;
			pxy0 = i + j*(Parptr->xsz + 1);
			px1 = i + 1 + j*(Parptr->xsz + 1);
			pyl = i + (j + 1)*(Parptr->xsz + 1);

			Xv = getmax(fabs(Arrptr->Vx[pxy0]), fabs(Arrptr->Vx[px1]));
			Yv = getmax(fabs(Arrptr->Vy[pxy0]), fabs(Arrptr->Vy[pyl]));

			Vc = sqrt(Xv*Xv + Yv*Yv);
			Haz = Arrptr->H[p0] * (Vc + C(1.5)); // Changed to equation from DEFRA 2006 (ALD)
			Arrptr->maxHaz[p0] = getmax(Arrptr->maxHaz[p0], Haz);
			if (Vc > Arrptr->maxVc[p0])
			{
				Arrptr->maxVc[p0] = Vc;
				Arrptr->maxVcH[p0] = Arrptr->H[p0];
			}
		}
	}
	return;
}
// Calculate discharge past a cross section
// this version uses the new 
NUMERIC_TYPE CalcVirtualGauge(const int gauge_i, const int grid_cols_padded, 
	const NUMERIC_TYPE * Qx_grid, const NUMERIC_TYPE * Qy_grid,
	Stage *Locptr)
{
	int j, q_index;
	NUMERIC_TYPE discharge = C(0.0);
	
	int x = Locptr->gauge_grid_x[gauge_i];
	int y = Locptr->gauge_grid_y[gauge_i];

	// work out the discharge from Q's
	if (Locptr->gauge_dir[gauge_i] == North || Locptr->gauge_dir[gauge_i] == South) // North or South flow
	{
		for (j = 0; j < Locptr->gauge_cells[gauge_i]; j++)
		{
			q_index = x + y * grid_cols_padded;
			// North or South flow so use Qy flows
			discharge += Qy_grid[q_index];
			x++;
		}
	}
	else if (Locptr->gauge_dir[gauge_i] == East || Locptr->gauge_dir[gauge_i] == West) // East or West flow
	{
		for (j = 0; j<Locptr->gauge_cells[gauge_i]; j++)
		{
			q_index = x + y * grid_cols_padded;
			// East or West flow so use Qx flows
			discharge += Qx_grid[q_index];
			y++;
		}
	}

	// switch signs for west and north flow
	if ((Locptr->gauge_dir[gauge_i] == North || Locptr->gauge_dir[gauge_i] == West) &&
		FABS(discharge) > C(0.0))
	{
		discharge = -discharge;
	}

	return(discharge);
}
//-----------------------------------------------------------------------------------
// CCS Populate dx, dy and dA arrays.  Calculates true cell dimensions if using a lat-long grid.
void CalcArrayDims(States *Statesptr, Pars *Parptr, Arrays *Arrptr)
{
	int i, j, p0;
	NUMERIC_TYPE pi, earth, deg_top, deg_bottom, deg_left, deg_right, rad_top, rad_bottom, rad_left, rad_right, delta_lat, delta_long, a, c;
	NUMERIC_TYPE dx, dy, min_dx, min_dy;

	if (Statesptr->latlong == OFF) // Regular grid with uniform square cells
	{
		for (j = 0; j < Parptr->ysz; j++)
		{
			for (i = 0; i < Parptr->xsz; i++)
			{
				p0 = i + j*Parptr->xsz; // location of cell
				Arrptr->dx[p0] = Parptr->dx;
				Arrptr->dy[p0] = Parptr->dy;
				Arrptr->dA[p0] = Parptr->dA;
			}
		}
		Parptr->min_dx = Parptr->dx;
		Parptr->min_dy = Parptr->dy;
		Parptr->min_dx_dy = Parptr->dx; // dx=dy for regular grid.
	}

	else // Lat-Long grid where dimensions of cells vary in both x and y:
	{
		//pi:
		pi = C(4.0)*atan(C(1.0));

		//earth's radius in metres (WGS 84 arithmetic mean radius of Earth)
		earth = 6371009;

		/* Loop through domain, calc cell boundary coordinates, calc cell dimensions:*/
		for (j = 0; j < Parptr->ysz; j++)
		{
			for (i = 0; i < Parptr->xsz; i++)
			{
				p0 = i + j*Parptr->xsz; // location of cell

				deg_bottom = Parptr->bly + (Parptr->ysz - j - 1)*Parptr->dx; // cell bottom boundary latitude
				deg_top = deg_bottom + Parptr->dx; // cell top boundary latitude
				deg_left = Parptr->blx + i*Parptr->dx; // cell left boundary longitude
				deg_right = deg_left + Parptr->dx; // cell right boundary longitude

				if (deg_left > 180)
				{
					deg_left = -180 + (deg_left - 180); // correct if left boundary longitude exceeds +180 degrees.
				}
				if (deg_right > 180)
				{
					deg_right = -180 + (deg_right - 180); // correct if right boundary longitude exceeds +180 degrees.
				}

				// convert to radians:
				rad_top = deg_top*(pi / 180);
				rad_bottom = deg_bottom*(pi / 180);
				rad_left = deg_left*(pi / 180);
				rad_right = deg_right*(pi / 180);

				// calc differences:
				delta_lat = rad_top - rad_bottom;
				delta_long = rad_left - rad_right;

				// calc dx distance using Haversine formula:
				a = cos(rad_bottom)*cos(rad_bottom)*pow(sin(delta_long / C(2.0)), C(2.0));
				c = 2 * atan2(sqrt(a), sqrt(C(1.0) - a));
				dx = earth*c;
				Arrptr->dx[p0] = dx;

				// calc dy distance using Haversine formula:
				a = pow(sin(delta_lat / C(2.0)), C(2.0));
				c = 2 * atan2(sqrt(a), sqrt(C(1.0) - a));
				dy = earth*c;
				Arrptr->dy[p0] = dy;

				// calc cell dA:
				Arrptr->dA[p0] = Arrptr->dx[p0] * Arrptr->dy[p0];

				if (j == 0 && i == 0) // if first cell, initialise min values:
				{
					min_dx = dx;
					min_dy = dy;
				}
				if (dx < min_dx) min_dx = dx; // Update min dx value if needed
				if (dy < min_dy) min_dy = dy; // Update min dy value if needed
			}
		}
		// Set Parptr->dx and Parptr->dy values to min values in array; calc min cell dimension:
		Parptr->min_dx = min_dx;
		Parptr->min_dy = min_dy;
		Parptr->min_dx_dy = getmin(min_dx, min_dy);
	}
}
