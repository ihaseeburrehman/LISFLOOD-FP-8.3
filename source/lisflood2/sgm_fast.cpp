/*
* sgm_fast.cpp
*
*  Created on: 14 May 2014
*      Author: td14281
*/

#include "../lisflood.h"
#include "../utility.h"
#include "sgm_fast.h"
#include <math.h>
#include <omp.h>

#include "../sgc.h"
#include "lis2_output.h"
#include "file_tool.h"

#if defined (__INTEL_COMPILER) && _PROFILE_MODE > 0
#include "ittnotify.h"
#endif

/*
OpenMP reduction variables defined as global.
This enables the pragma omp parallel to be defined in a separate
function to the pragma omp for reduction(...)
Generally a bad practice to define variables as global
as they can soon become a mess.
*/
NUMERIC_TYPE reduce_Hmax;
NUMERIC_TYPE reduce_evap_loss;
NUMERIC_TYPE reduce_rain_total;
NUMERIC_TYPE reduce_Qpoint_timestep_pos;
NUMERIC_TYPE reduce_Qpoint_timestep_neg;

NUMERIC_TYPE reduce_flood_area;
NUMERIC_TYPE reduce_domain_volume;

NUMERIC_TYPE evap_deltaH_step;
NUMERIC_TYPE rain_deltaH_step;

int itCount;
double total_write_time = C(0.0);

///
/// Calculate channel flow using inertial wave equation --- in this case q_old and returned q will be in m3s-1
/// 
inline NUMERIC_TYPE CalculateQ(const NUMERIC_TYPE surface_slope,
	NUMERIC_TYPE R, // hydraulic radius
	const NUMERIC_TYPE delta_time,
	const NUMERIC_TYPE g,
	const NUMERIC_TYPE area, //flow area (not cell area)
	const NUMERIC_TYPE g_friction_squared,
	const NUMERIC_TYPE q_old,
	const NUMERIC_TYPE max_Froude)
{

	// calculate flow based on m^3 formula (note power is 4/3, profiling shows performance gain by multiply and cuberoot)
#if _CALCULATE_Q_MODE == 0
	NUMERIC_TYPE pow_tmp1, pow_tmp, abs_q, calc_num, calc_den;

	//R = getmin(R, 10.0); // removed so depth of flow doesn't max out at 10m (JCN)

	pow_tmp1 = R * R * R * R;
	pow_tmp = CBRT(pow_tmp1); // 4 multiplies and 1 cube root profiled faster than POW(R,4/3)

	abs_q = FABS(q_old);

	calc_num = (q_old - g * area * delta_time * surface_slope);
	calc_den = (1 + delta_time * g_friction_squared * abs_q / (pow_tmp * area));
	
	return calc_num / calc_den;
	
#else
#if _CALCULATE_Q_MODE == 1
	NUMERIC_TYPE pow_tmp1, pow_tmp, abs_q, calc_num, calc_den;

	//R = getmin(R, 10.0); // removed so depth of flow doesn't max out at 10m (JCN)

	pow_tmp1 = R * R * R * R;
	pow_tmp = CBRT(pow_tmp1); // 4 multiplies and 1 cube root profiled faster than POW(R,4/3)

	abs_q = FABS(q_old);

	calc_num = (q_old - g * area * delta_time * surface_slope);
	calc_den = (1 + delta_time * g_friction_squared * abs_q / (pow_tmp * area));

	// limit to max Froude
	calc_num = calc_num / calc_den; // calculate Q as calc_num
	calc_den = max_Froude*area*SQRT(g*R); // Calculate max Q for max Froude as calcden

	if (FABS(calc_num) < calc_den) return calc_num; // return calc_num if its less than calc_den
	else return copysign(1.0, calc_num)*calc_den; // else return calc_den but get the sign from the surface slope

#else
#if _CALCULATE_Q_MODE == 2
	NUMERIC_TYPE pow_tmp1, pow_tmp, abs_q, calc_num, calc_den;

	pow_tmp = POW(R, C(4.0)/C(3.0));

	abs_q = FABS(q_old);

	calc_num = (q_old - g * area * delta_time * surface_slope);
	calc_den = (1 + delta_time * g_friction_squared * abs_q / (pow_tmp * area));
	return calc_num / calc_den;
#else

	return (q_old - g * area * delta_time * surface_slope) / ((1 + delta_time * g_friction_squared * FABS(q_old) / (POW(R, C(4.0)/C(3.0)) * area)));
#endif
#endif
#endif
}

inline NUMERIC_TYPE SGC2_CalculateVelocity(const int index, const int index_next,
	const NUMERIC_TYPE * Q_grid,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE width)
{
	if (Q_grid[index_next] != C(0.0) && (h_grid[index] > C(0.0) || h_grid[index_next] > C(0.0)))
	{
		NUMERIC_TYPE h0 = h_grid[index];
		NUMERIC_TYPE h1 = h_grid[index_next];
		NUMERIC_TYPE z0 = dem_grid[index];
		NUMERIC_TYPE z1 = dem_grid[index_next];

		NUMERIC_TYPE surface_elevation0 = z0 + h0;
		NUMERIC_TYPE surface_elevation1 = z1 + h1;
		// Calculating hflow based on floodplain levels
		NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
			- getmax(z0, z1);

		return (Q_grid[index_next] / width) / hflow;
	}
	else
	{
		return C(0.0);
	}
}

//-----------------------------------------------------------------------------------
// This function calculates the hydraulic radius of a sub-grid channel given the channel type
inline NUMERIC_TYPE SGC2_CalcR(int gr, NUMERIC_TYPE h, NUMERIC_TYPE hbf, NUMERIC_TYPE w, NUMERIC_TYPE wbf, NUMERIC_TYPE A, const SGCprams *SGCptr)
{
#if defined _ONLY_RECT && _ONLY_RECT == 1
	NUMERIC_TYPE  R = getmin(h, hbf); // don't exceed bankfull
	return R = A / (w + 2 * R); // calculate hydraulic radius
#else
	// This function calculates the hydraulic redius of the channel given the flow area (A) and depth (h)
	NUMERIC_TYPE R = C(0.0), sl, hp, b1, b2;

	switch (SGCptr->SGCchantype[gr])
	{
	case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		R = getmin(h, hbf); // don't exceed bankfull
		R = A / (w + 2 * R); // calculate hydraulic radius
		break; // Break terminates the switch statement

	case 2: // Exponent channel       
		h = getmin(h, hbf); // don't exceed bankfull
		hp = 2 * h / wbf; // use half width
		// get beta parameters for the channel group number
		if (hp <= SGCptr->SGCbetahmin)
		{
			// beta parameters for flow depths below SGCbetahmin (C(0.05)) depth/bankfull width
			b1 = SGCptr->SGCbeta1[gr]; b2 = SGCptr->SGCbeta2[gr];
			R = b1*hp + b2*hp*hp; // calculate wetted perimiter component
		}
		else
		{
			// beta parameters for flow depths above or equal to SGCbetahmin (C(0.05)) depth/bankfull width
			b1 = SGCptr->SGCbeta3[gr]; b2 = SGCptr->SGCbeta4[gr];
			hp -= SGCptr->SGCbetahmin; // decrease hp to account for for wetted perimiter fraction below C(0.05) depth/bankful depth
			R = SGCptr->SGCbeta5[gr] + b1*hp + b2*hp*hp; // calculate wetted perimiter component
		}
		R = A / (w + R*wbf);  // calculate hydraulic radius
		break;	// Break terminates the switch statement

	case 3: // linear slope - This model has a bed elevation, slope, top width and top elevation. 
		if (h < hbf)	R = A / (h + SQRT(h  *h + w*w));	// within bank flow hydraulic radius
		else			R = A / (hbf + SQRT(hbf*hbf + w*w));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
		break;	// Break terminates the switch statement

	case 4: // triangular channel - This model has the bed elevation, slope, top width and top elevation 
		w = C(0.5)*w;
		if (h < hbf)	R = A / (2 * SQRT(h  *h + w*w));   // within bank flow hydraulic radius
		else			R = A / (2 * SQRT(hbf*hbf + w*w));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
		break;	// Break terminates the switch statement  

	case 5: // parabolic channel       
		h = getmin(h, hbf); // don't exceed bankfull
		w = w / C(2.0); // half width
		R = SQRT(w*w + C(16.0)*h*h);
		R = C(0.5)*R + (w*w) / (C(8.0)*h) * log((C(4.0)*h + R) / w);
		R = C(2.0)*R;
		R = A / R;  // calculate hydraulic radius  */
		break;	// Break terminates the switch statement

	case 6: // Rectangular channel (no banks) - This model has a top width and bed elevation and top 
		R = A / w; // calculate hydraulic radius
		break; // Break terminates the switch statement

	case 7: // trapazoidal channel
		sl = SGCptr->SGCs[gr];
		if (h < hbf)	R = A / (w + 2 * h   * SQRT(1 + sl*sl));   // within bank flow hydraulic radius
		else			R = A / (w + 2 * hbf * SQRT(1 + sl*sl));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
		break;	// Break terminates the switch statement

	default: // its all gone wrong
		printf("should not be here! Something is wrong with the SGC channel model R calculation");
		break;
	}
	return(R);
#endif
}

inline void SGC2_CalcA(int gr, NUMERIC_TYPE hflow, NUMERIC_TYPE bf, NUMERIC_TYPE *A, NUMERIC_TYPE *we, const SGCprams *SGCptr)
{
#if defined _ONLY_RECT && _ONLY_RECT == 1
	(*A) = (*we)*hflow;
#else

	// This function calculates the area of flow (A) for a given flow depth (hflow), in some cases it
	// also returnes the widths of flow (We).
	NUMERIC_TYPE sl;
	// switch depending on the channel type
	switch (SGCptr->SGCchantype[gr])
	{
	case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		(*A) = (*we)*hflow;
		break;	// Break terminates the switch statement

	case 2: // Power shaped channel. h = x^sl channel       
		sl = SGCptr->SGCs[gr];
		if (hflow < bf)
		{
			(*we) = (*we)*pow(hflow / bf, C(1.0) / sl);
			(*A) = hflow*(*we)*(C(1.0) - C(1.0) / (sl + C(1.0)));
		}
		else (*A) = bf*(*we)*(1 - 1 / (sl + 1)) + (*we)*(hflow - bf); // out of bank flow area
		break;	// Break terminates the switch statement

	case 3: // linear slope. 
		if (hflow < bf)
		{
			(*we) = ((*we) / bf)*hflow;
			(*A) = (*we)*hflow*C(0.5);    // within bank flow area
		}
		else  (*A) = (*we)*bf*C(0.5) + (*we)*(hflow - bf); // out of bank flow area
		break;	// Break terminates the switch statement

	case 4: // triangular channel
		if (hflow < bf)
		{
			(*we) = ((*we) / bf)*hflow;
			(*A) = (*we)*hflow*C(0.5);    // within bank flow area
		}
		else  (*A) = (*we)*bf*C(0.5) + (*we)*(hflow - bf); // out of bank flow area
		break;	// Break terminates the switch statement  

	case 5: // parabolic channel       
		if (hflow < bf)
		{
			(*we) = (*we)*SQRT(hflow / bf);
			(*A) = hflow*(*we)*(C(2.0) / C(3.0));
		}
		else    (*A) = bf*   (*we)*(C(2.0) / C(3.0)) + (*we)*(hflow - bf); // out of bank flow area
		break;	// Break terminates the switch statement

	case 6: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		(*A) = (*we)*hflow;
		break;	// Break terminates the switch statement

	case 7: // trapazoidal channel
		sl = SGCptr->SGCs[gr];
		if (hflow < bf) (*A) = ((*we) + sl*hflow)*hflow;    // within bank flow area
		else            (*A) = ((*we) + sl*bf)*bf + ((*we) + sl*bf)*(hflow - bf);
		break;	// Break terminates the switch statement

	default: // its all gone wrong
		printf("should not be here! Something is wrong with the SGC channel model A calculation");
		break;
	}
	return;
#endif
}

/*

case 7: // trapazoidal channel
// calculate hydraulic radius
if (hflow < bf)	R = A / (we + 2* hflow * SQRT(1+sl*sl));   // within bank flow hydraulic radius
else			R = A / (we + 2* bf    * SQRT(1+sl*sl));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
break;	// Break terminates the switch statement

*/


/// cell_width: cross section width
/// 
/// output: Q_FP_corrected flood plain flow minus any region of flow above the sub-grid (Neal 2012 Figure 1 (C) )
/// output: Q_FP_old flood flood plain flow (not corrected, use for next iteration input)
/// output: Q_SG_old flood flood plain flow (use for next iteration input)
inline NUMERIC_TYPE SGC2_CalcPointFREE(const NUMERIC_TYPE hflow, const NUMERIC_TYPE SGC_width, const NUMERIC_TYPE Sf, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE cell_width, const NUMERIC_TYPE g, const NUMERIC_TYPE g_friction_squared_SG, const NUMERIC_TYPE g_friction_squared_FP, const NUMERIC_TYPE SGC_BankFullHeight, const int gr, const int sign,
	NUMERIC_TYPE *Q_FP_old, NUMERIC_TYPE *Q_SG_old,
	const SGCprams *SGCptr, const NUMERIC_TYPE max_Froude)
{
	NUMERIC_TYPE Q_FP_corrected;
	NUMERIC_TYPE A, R, SGC_width_current;
	NUMERIC_TYPE abs_q, SGC_hflow;
	SGC_width_current = SGC_width; // save bank full width
	if (SGC_width > C(0.0)) // channel present
	{
		SGC_hflow = hflow + SGC_BankFullHeight;

		SGC2_CalcA(gr, SGC_hflow, SGC_BankFullHeight, &A, &SGC_width_current, SGCptr); // calculate channel area for SGC
		R = SGC2_CalcR(gr, SGC_hflow, SGC_BankFullHeight, SGC_width_current, SGC_width, A, SGCptr); // calculate hydraulic radius for SGC
		// Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
		abs_q = FABS(*Q_SG_old);

		/// set Sf to be negative to ensure the numerator in CalculateQ will be positive (Sf used in subtraction)
		*Q_SG_old = sign * CalculateQ(-1 * FABS(Sf), R, delta_time, g, A, g_friction_squared_SG, abs_q, max_Froude);
	}
	//	else
	//	{
	//#ifdef _DEBUG
	//		// Toby Dunne removed redundant *Q_SG_old = C(0.0);
	//		if (*Q_SG_old != C(0.0))
	//			printf("SGC2_CalcPointFREE expect non sub grid cell to have no sub grid flow\n");
	//#endif
	//	}
	// multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
	// FABS on Sf and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
	//if(hflow>Solverptr->DepthThresh && SGC_width_current < Parptr->dx) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell //CCS_deletion
	if (hflow > depth_thresh && SGC_width_current < cell_width) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell
	{
		// calculate FP flow
		//*qoldptr = sign*(FABS(*qoldptr) + FABS(g*Tstep*hflow*Sf)) / (1 + Tstep*g_friction_squared_FP*FABS(*qoldptr) / (pow(hflow, (C(7.) / C(3.)))));
		A = cell_width * hflow;
		R = hflow;
		abs_q = FABS(*Q_FP_old);

		NUMERIC_TYPE q;

		/// set Sf to be negative to ensure the numerator in CalculateQ will be positive (Sf used in subtraction)
		q = sign * CalculateQ(-1 * FABS(Sf), R, delta_time, g, A, g_friction_squared_FP, abs_q, max_Froude);
		*Q_FP_old = q;

		// subtract channel flow from flood plain q
		if (SGC_width_current > C(0.0))
		{
			NUMERIC_TYPE channel_ratio = min(SGC_width_current / cell_width, C(1.0));
			q = q - channel_ratio * q;
		}
		Q_FP_corrected = q;
	}
	else
	{
		*Q_FP_old = C(0.0);
		Q_FP_corrected = C(0.0);
	}
	return Q_FP_corrected;
}


//-----------------------------------------------------------------------------------
// This function updates H for a sub-grid channel given the channel type
inline NUMERIC_TYPE SGC2_CalcUpH(const NUMERIC_TYPE V, const NUMERIC_TYPE c, const int channel_group, const SGCprams *SGCptr)
{
#if defined _ONLY_RECT && _ONLY_RECT == 1
	return V / c;
#else
	NUMERIC_TYPE h = C(0.0);

	int SGCchan_type = SGCptr->SGCchantype[channel_group];
	NUMERIC_TYPE sl = SGCptr->SGCs[channel_group];

	// switch to the correct sub-grid channel, default 1 is the rectangular
	switch (SGCchan_type)
	{
	case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		h = V / c;
		break;

	case 2: // y = x^sl channel       
		h = pow(V / c, sl / (sl + C(1.0)));
		break;  // Break terminates the switch statement

	case 3: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		h = V / c;
		break;

	case 4: // Rectangular channel (default) - This model has a top width and bed elevation and top 
		h = V / c;
		break;

	case 5: // parabiolic channel       
		h = pow(V / c, C(2.0) / C(3.0));
		break;  // Break terminates the switch statement

	case 6: // Rectangular channel (no banks) - This model has a top width and bed elevation and top 
		h = V / c;
		break;

	default: // its all gone wrong
		printf("should not be here! Something is wrong with the SGC channel model in SGC_UpdateH");
		break;
	} // end of switch statement
	return(h);
#endif
}

inline NUMERIC_TYPE SGC2_CalcUpV(const NUMERIC_TYPE h, const NUMERIC_TYPE c, const int channel_group, const SGCprams *SGCptr)
{
#if defined _ONLY_RECT && _ONLY_RECT == 1
	return h*c;
#else

	int SGCchan_type = SGCptr->SGCchantype[channel_group];
	NUMERIC_TYPE sl = SGCptr->SGCs[channel_group];

	// This function calculates the volume of a sub-grid channel given a depth
	NUMERIC_TYPE v = C(0.0);
	// switch to the correct sub-grid channel, default 1 is the rectangular
	switch (SGCchan_type)
	{
	case 1:
		v = h*c;
		break;

	case 2:
		v = c*pow(h, C(1.0) / sl + C(1.0));
		break;

	case 3:
		v = h*c;
		break;

	case 4:
		v = h*c;
		break;

	case 5:
		v = c*pow(h, C(3.0) / C(2.0));
		break;

	case 6: // Rectangular channel (no banks)
		v = h*c;
		break;

	default: // its all gone wrong
		printf("should not be here! Something is wrong with the SGC channel model in SGC2_CalcUpV");
		break;
	} // end of switch statement
	return(v);
#endif
}



inline NUMERIC_TYPE CalculateRoutingQ(const NUMERIC_TYPE delta_time,
	const NUMERIC_TYPE h0, NUMERIC_TYPE h1,
	const NUMERIC_TYPE z0, const NUMERIC_TYPE z1,
	const NUMERIC_TYPE route_V_ratio_per_sec, const NUMERIC_TYPE cell_area)
{
	//h1 = Arrptr->H[Arrptr->FlowDir[p0]] - Arrptr->SGCbfH[Arrptr->FlowDir[p0]];
	if (h1 < 0)
	{
		h1 = C(0.0);
	}	// If h1 negative due to below bankful SG channel cell, set h1 to zero 

	//z0 = Arrptr->DEM[p0]; //cell DEM height
	//z1 = Arrptr->DEM[Arrptr->FlowDir[p0]]; //lowest neighbour cell DEM height 

	NUMERIC_TYPE flow = (z0 + h0) - (z1 + h1);/*calculate the maximum possible flow into lowest neigbour cell by
											  comparing water surface elevations:*/
	/*where water surface elevation of neighbour cell is below DEM
	level of current cell, set maxflow to h0*/
	if (flow > h0)
	{
		flow = h0;
	}
	/*where water surface elevation of neighbour cell is above water TFD - can never happen as route only triggered if high surface slope
	surface elevation of current cell, set maxflow to 0*/
	if (flow < 0)
	{
		flow = 0;
	}

	NUMERIC_TYPE flow_fraction_s = route_V_ratio_per_sec; // fraction of cell volume to route in this time step
	if ((flow_fraction_s * delta_time) > 1)
		flow_fraction_s = 1 / delta_time; // more than 100 % of the volume being removed pre time step

	NUMERIC_TYPE flowQ = flow * cell_area * flow_fraction_s; // q is volume per second, when q converted to voume in SGC2_UpdateVol_floodplain_row q is multiplied by delta_time

	return flowQ;
}

//-----------------------------------------------------------------------------------
// CALCULATE WIER FLOW BETWEEN A POINT AND ITS W NEIGHBOUR
inline NUMERIC_TYPE CalcWeirQ(WeirLayout * weirs, const int grid_cols,
	const int grid_index_this, const int grid_index_next,
	const int weir_id,
	const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE delta_time,
	const NUMERIC_TYPE * h_grid, 
	const NUMERIC_TYPE * volume_grid,
	WetDryRowBound* wet_dry_bounds,
	const EDirection dir_positive, const EDirection dir_negative)
{
	NUMERIC_TYPE z0, z1, h0, h1, Q;
	//NUMERIC_TYPE usVel; // MT upstream velocity for energy gradient height calc.
	//NUMERIC_TYPE heg; // MT upstream energy gradient height
	//int p0, p1, pq0, weir_id;

	const int weir_pair_id0 = 2 * weir_id;
	const int weir_pair_id1 = weir_pair_id0 + 1;

	z0 = weirs->cell_pair.sg_cell_dem[weir_pair_id0];
	z1 = weirs->cell_pair.sg_cell_dem[weir_pair_id1];

	h0 = h_grid[grid_index_this];
	h1 = h_grid[grid_index_next];

	NUMERIC_TYPE surfaceElevation0 = h0 + z0;
	NUMERIC_TYPE surfaceElevation1 = h1 + z1;

	h0 += weirs->cell_pair.sg_cell_SGC_BankFullHeight[weir_pair_id0];
	h1 += weirs->cell_pair.sg_cell_SGC_BankFullHeight[weir_pair_id1];

	//int x = weirs->cell_pair.sg_cell_x[weir_pair_id0];
	//int y = weirs->cell_pair.sg_cell_y[weir_pair_id0];

	Q = C(0.0);
	if (h0 > depth_thresh || h1 > depth_thresh)
	{

		if ((surfaceElevation0) > (surfaceElevation1))		// Flow in + direction
		{
			if ((surfaceElevation0) > weirs->Weir_hc[weir_id] && h0 > 0) // check depth is above weir and that the cell is wet
			{
				if (weirs->Weir_Fixdir[weir_id] == DirectionNA || weirs->Weir_Fixdir[weir_id] == dir_positive) // check for one-directional flow (culvert)
				{
					NUMERIC_TYPE hu = surfaceElevation0 - weirs->Weir_hc[weir_id]; // upstream head
					NUMERIC_TYPE hd = surfaceElevation1 - weirs->Weir_hc[weir_id]; // downstream head
					if ((hd / hu) < weirs->Weir_m[weir_id])
						Q = weirs->Weir_Cd[weir_id] * weirs->Weir_w[weir_id] * pow(hu, (C(1.5))); // Free flow
					else
						Q = weirs->Weir_Cd[weir_id] * weirs->Weir_w[weir_id] * hu*(SQRT(hu - hd)) / SQRT(weirs->Weir_m[weir_id]); // Drowned flow
					NUMERIC_TYPE maxQ=((volume_grid[grid_index_this]/delta_time)*0.5);
					Q=getmin(Q,maxQ);				
}
			}
                                 
		}
		else if ((surfaceElevation0) < (surfaceElevation1))		// Flow in - direction
		{
			if ((surfaceElevation1) > weirs->Weir_hc[weir_id] && h1 > 0) // check depth is above weir and that the cell is wet
			{
				if (weirs->Weir_Fixdir[weir_id] == DirectionNA || weirs->Weir_Fixdir[weir_id] == dir_negative) // check for one-directional flow (culvert)
				{
					NUMERIC_TYPE hu = surfaceElevation1 - weirs->Weir_hc[weir_id]; // upstream head
					NUMERIC_TYPE hd = surfaceElevation0 - weirs->Weir_hc[weir_id]; // downstream head
					if ((hd / hu) < weirs->Weir_m[weir_id])
						Q = -weirs->Weir_Cd[weir_id] * weirs->Weir_w[weir_id] * pow(hu, (C(1.5))); // Free flow
					else
						Q = -weirs->Weir_Cd[weir_id] * weirs->Weir_w[weir_id] * hu*(SQRT(hu - hd)) / SQRT(weirs->Weir_m[weir_id]); // Drowned flow
						NUMERIC_TYPE maxQ=((volume_grid[grid_index_this]/delta_time)*0.5);
						Q=getmin(Q,maxQ);				
}
			}
		}
	}
	if (Q != C(0.0))
	{
		int x = weirs->cell_pair.sg_cell_x[weir_pair_id0];
		int y = weirs->cell_pair.sg_cell_y[weir_pair_id0];
		wet_dry_bounds->fp_vol[y].start = min(wet_dry_bounds->fp_vol[y].start, x);
		wet_dry_bounds->fp_vol[y].end = max(wet_dry_bounds->fp_vol[y].end, x + 1);
	}

	

	return(Q);
}


//-----------------------------------------------------------------------------------
// CALCULATE WIER FLOW BETWEEN A POINT AND ITS W NEIGHBOUR
inline NUMERIC_TYPE CalcBridgeQ(WeirLayout * weirs, const int grid_cols,
	const int grid_index_this, const int grid_index_next,
	const int weir_id,
	const NUMERIC_TYPE g, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE cell_length,
	const NUMERIC_TYPE * h_grid,
	WetDryRowBound* wet_dry_bounds,
	const SubGridState * sub_grid_state,
	const NUMERIC_TYPE max_Froude)
{
	NUMERIC_TYPE z0, z1, h0, h1, Q;
	NUMERIC_TYPE usVel; // MT upstream velocity for energy gradient height calc.
	NUMERIC_TYPE heg; // MT upstream energy gradient height
	//int p0, p1, pq0, weir_id;


	const int weir_pair_id0 = 2 * weir_id;
	const int weir_pair_id1 = weir_pair_id0 + 1;

	z0 = weirs->cell_pair.sg_cell_dem[weir_pair_id0];
	z1 = weirs->cell_pair.sg_cell_dem[weir_pair_id1];

	h0 = h_grid[grid_index_this];
	h1 = h_grid[grid_index_next];

	NUMERIC_TYPE surfaceElevation0 = h0 + z0;
	NUMERIC_TYPE surfaceElevation1 = h1 + z1;

	h0 += weirs->cell_pair.sg_cell_SGC_BankFullHeight[weir_pair_id0];
	h1 += weirs->cell_pair.sg_cell_SGC_BankFullHeight[weir_pair_id1];

	//int x = weirs->cell_pair.sg_cell_x[weir_pair_id0];
	//int y = weirs->cell_pair.sg_cell_y[weir_pair_id0];

	Q = C(0.0);
	if (h0 > depth_thresh || h1 > depth_thresh)
	{
		NUMERIC_TYPE Qoc; // open channel flow
		NUMERIC_TYPE Qp; // pressure(orifice) flow
		NUMERIC_TYPE Tz; // transit zone
		NUMERIC_TYPE Cd; // bridge Cd
		NUMERIC_TYPE Soffit; // bridge soffit elevation
		NUMERIC_TYPE Area; // bridge open area - precaclulated in input
		NUMERIC_TYPE Z; //bridge opening
		NUMERIC_TYPE Zratio; // flow depth to bridge opening ratio
		NUMERIC_TYPE Width; // bridge width
		NUMERIC_TYPE dh; // flow head change (for open channel flow calc)
		NUMERIC_TYPE Sf; // friction slope (for open channel flow calc)
		NUMERIC_TYPE hflow; // flow depth
		NUMERIC_TYPE A; // open channel flow area
		NUMERIC_TYPE R; // open channel hydraulic radius
		// get basic bridge params from global arrays;
		Tz = weirs->Weir_m[weir_id];
		Soffit = weirs->Weir_hc[weir_id];
		Cd = weirs->Weir_Cd[weir_id];
		Width = weirs->Weir_w[weir_id];

		// calculate some more bridge parameters from basic ones
		Z = getmin(Soffit - z1, Soffit - z0); // bridge opening (smallest opening)
		Area = Width*Z; // bridge flow area (again smallest opening)

		// get some basic paramters for the open channel flow calc
		NUMERIC_TYPE g_friction_sq = weirs->Weir_g_friction_sq[weir_id];// cn = C(0.5)* (SGCptr->SGCn[SGCgroup_grid[p0]] + SGCptr->SGCn[SGCgroup_grid[p1]]); // mean mannings (note n2)
		dh = surfaceElevation0 - surfaceElevation1; // difference in water level
		//Sf=-dh/Parptr->dx; //CCS_deletion
		Sf = -dh / cell_length; // CCS added for subgrid lat long compatibility.
		hflow = getmax(surfaceElevation0, surfaceElevation1) - getmax(z0, z1); // calculate the max flow depth
		A = Width*hflow; // calc open channel flow area
		R = A / (Width + 2 * hflow); // calc hydraulic radius from open channel flow

		// calculate open channel flow using SGC flow
		Qoc = CalculateQ(Sf, R, delta_time, g, A, g_friction_sq, weirs->Weir_Q_old_SG[weir_id], max_Froude);

		// orifice bridge flow
		if (surfaceElevation0 > surfaceElevation1) { // Positive flow
			Zratio = h0 / Z; // calc current flow depth to bridge opening ratio
			usVel = sub_grid_state->sg_flow_Q[weirs->Weir_pair_stream_flow_index[weir_pair_id0]] / (h0*weirs->cell_pair.sg_cell_SGC_width[weir_pair_id0]); // MT calculate upstream velocity
			heg = (usVel*usVel) / (2 * g); // MT calculate upstream energy gradient height
			Qp = Cd*Area*SQRT(2 * g*(surfaceElevation0 - surfaceElevation1 + heg)); // calc bridge orifice flow
		}
		else { // Negative flow
			Zratio = h1 / Z; // calc current flow depth to bridge opening ratio
			usVel = sub_grid_state->sg_flow_Q[weirs->Weir_pair_stream_flow_index[weir_pair_id1]] / (h1*weirs->cell_pair.sg_cell_SGC_width[weir_pair_id1]); // MT calculate upstream velocity
			heg = (usVel*usVel) / (2 * g); // MT calculate upstream energy gradient height
			Qp = -(Cd*Area*SQRT(2 * g*(surfaceElevation1 - surfaceElevation0 + heg))); // calc bridge orifice flow
		}


		if (surfaceElevation0 < Soffit && surfaceElevation1 < Soffit)
		{
			// flow is below the soffit so use SGC open channel flow 
			Q = Qoc;
		}
		else if (Zratio >= C(1.0) && Zratio <= Tz) // transition flow between open and orifice/pressure
		{
			Q = (Qoc*(Tz - Zratio) / (Tz - C(1.0))) + (Qp*(Zratio - C(1.0)) / (Tz - C(1.0)));
		}
		else if (Zratio > Tz) // pressure flow
		{
			Q = Qp;
		}
		else  // other flow - should not happen, but put here to catch other cases?
		{
			printf("WARNING: Unexpected Bridge flow calc fail at t=%.3" NUM_FMT" , Soffit %" NUM_FMT" m.\n", curr_time, Soffit); //Warning for fail
			Q = Qoc;
		}
	}
	if (Q != C(0.0))
	{
		int x = weirs->cell_pair.sg_cell_x[weir_pair_id0];
		int y = weirs->cell_pair.sg_cell_y[weir_pair_id0];
		wet_dry_bounds->fp_vol[y].start = min(wet_dry_bounds->fp_vol[y].start, x);
		wet_dry_bounds->fp_vol[y].end = max(wet_dry_bounds->fp_vol[y].end, x + 1);
	}
	return(Q);
}

/*
* Qx: next_cell_add = 1
* Qy: next_cell_add = grid_cols_padded
*/
void SGC2_UpdateVelocity_row(const int grid_row_index,
	const int row_start_prev, const int row_end_prev,
	const int row_start, const int row_end,
	const int next_cell_add,
	const NUMERIC_TYPE cell_width,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE * Q_grid,
	NUMERIC_TYPE * V_grid, NUMERIC_TYPE * V_max_grid
	)
{

#pragma ivdep
#pragma simd
	// clear from start of prev bound to start of new bound
	for (int i = row_start_prev; i < row_start; i++)
	{
		int index_next = grid_row_index + i + next_cell_add;
		V_grid[index_next] = C(0.0);
	}
	if (row_end != -1)
	{
#pragma ivdep
#pragma simd
		// clear from end bound to end of prev bound
		for (int i = row_end; i < row_end_prev; i++)
		{
			int index_next = grid_row_index + i + next_cell_add;
			V_grid[index_next] = C(0.0);
		}
	}

#pragma simd
	for (int i = row_start; i < row_end; i++)
	{
		int index = grid_row_index + i;
		//next column
		int index_next = index + next_cell_add;

		NUMERIC_TYPE velocity = SGC2_CalculateVelocity(index, index_next,
			Q_grid,
			h_grid, dem_grid, cell_width);
		V_grid[index_next] = velocity;
		V_max_grid[index_next] = getmax(FABS(velocity), V_max_grid[index_next]);
	}
}

/// update Q for weir flows
/// separate from update Q as the bridges depend on adjacent q values.
void SGC2_UpdateWeirsFlow_row(const int j, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE delta_time,
	const NUMERIC_TYPE * h_grid, 
	const NUMERIC_TYPE * volume_grid,
	NUMERIC_TYPE * Qx_grid, 
	NUMERIC_TYPE * Qy_grid,
	WetDryRowBound * wet_dry_bounds,
	WeirLayout * weirs)
{
	const int weir_row_cols_padded = weirs->row_cols_padded;

	/// weirs are processed after UpdateQ completed
	/// bridges depend on q's upstream/downstream
	/// which if processed in updateQ, will never be updated in E->W, and randomly available in the case of N->S, S>N
	const int row_weir_start = j * weir_row_cols_padded;

	const int weir_Qx_row_count = weirs->weir_Qx_row_count[j];
	for (int weir_row_index = 0; weir_row_index < weir_Qx_row_count; weir_row_index++)
	{
		const int weir_id = weirs->weir_index_qx[weir_row_index + row_weir_start];

		const int grid_index_this = weirs->Weir_grid_index[weir_id];
		const int grid_index_next = grid_index_this + 1;

		NUMERIC_TYPE Q = CalcWeirQ(weirs, grid_cols,
			grid_index_this, grid_index_next, weir_id,
			depth_thresh, delta_time,
			h_grid, volume_grid,
			wet_dry_bounds, East, West);
		// update flow array
		weirs->Weir_Q_old_SG[weir_id] = Q;

		Qx_grid[grid_index_next] = Q;
	}
	const int weir_Qy_row_count = weirs->weir_Qy_row_count[j];
	for (int weir_row_index = 0; weir_row_index < weir_Qy_row_count; weir_row_index++)
	{
		const int weir_id = weirs->weir_index_qy[weir_row_index + row_weir_start];

		const int grid_index_this = weirs->Weir_grid_index[weir_id];
		const int grid_index_next = grid_index_this + grid_cols_padded;

		NUMERIC_TYPE Q = CalcWeirQ(weirs, grid_cols,
			grid_index_this, grid_index_next, weir_id,
			depth_thresh, delta_time,
			h_grid, volume_grid,
			wet_dry_bounds, South, North);
		// update flow array
		weirs->Weir_Q_old_SG[weir_id] = Q;

		Qy_grid[grid_index_next] = Q;
	}
}

/// update Q for weir flows
/// separate from update Q as the bridges depend on adjacent q values.
void SGC2_UpdateBridgesFlow_row(const int j, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE g,
	const NUMERIC_TYPE * dx_col, const NUMERIC_TYPE * dy_col,
	const NUMERIC_TYPE * h_grid,
	NUMERIC_TYPE * Qx_grid, NUMERIC_TYPE * Qy_grid,
	WetDryRowBound * wet_dry_bounds,
	SubGridState * sub_grid_state,
	WeirLayout * weirs,
	const NUMERIC_TYPE max_Froude)
{
	const int weir_row_cols_padded = weirs->row_cols_padded;

	/// weirs are processed after UpdateQ completed
	/// bridges depend on q's upstream/downstream
	/// which if processed in updateQ, will never be updated in E->W, and randomly available in the case of N->S, S>N
	const int row_weir_start = j * weir_row_cols_padded;

	const int weir_Qx_row_count = weirs->weir_Qx_row_count[j];
	for (int weir_row_index = 0; weir_row_index < weir_Qx_row_count; weir_row_index++)
	{
		const int weir_id = weirs->weir_index_qx[weir_row_index + row_weir_start];

		const int grid_index_this = weirs->Weir_grid_index[weir_id];
		const int grid_index_next = grid_index_this + 1;

		NUMERIC_TYPE Q = CalcBridgeQ(weirs, grid_cols,
			grid_index_this, grid_index_next, weir_id,
			g, delta_time, curr_time, depth_thresh,
			dx_col[j], h_grid,
			wet_dry_bounds, sub_grid_state, max_Froude);
		// update flow array
		weirs->Weir_Q_old_SG[weir_id] = Q;

		Qx_grid[grid_index_next] = Q;
	}
	const int weir_Qy_row_count = weirs->weir_Qy_row_count[j];
	for (int weir_row_index = 0; weir_row_index < weir_Qy_row_count; weir_row_index++)
	{
		const int weir_id = weirs->weir_index_qy[weir_row_index + row_weir_start];

		const int grid_index_this = weirs->Weir_grid_index[weir_id];
		const int grid_index_next = grid_index_this + grid_cols_padded;

		NUMERIC_TYPE Q = CalcBridgeQ(weirs, grid_cols,
			grid_index_this, grid_index_next, weir_id,
			g, delta_time, curr_time, depth_thresh,
			dy_col[j], h_grid,
			wet_dry_bounds, sub_grid_state, max_Froude);
		// update flow array
		weirs->Weir_Q_old_SG[weir_id] = Q;

		Qy_grid[grid_index_next] = Q;
	}
}
// Dam Code FEOL 2016
void DamOpQ(const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE * h_grid, NUMERIC_TYPE * Qx_grid, NUMERIC_TYPE * Qy_grid, DamData *Damptr, const NUMERIC_TYPE g, const int ID)
{
	if (Damptr->DamOperationCode[ID] == 1)
	{
		Damptr->DamOperationQ[ID] = Damptr->DamMeanQ[ID];

		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID] / delta_time;
		}
		
	}
	if (Damptr->DamOperationCode[ID] == 2) // Based off matlab code provided by Francesca Pianosi for constant release (regardless of storage)
	{
		NUMERIC_TYPE upper_limit = getmin((Damptr->Volmax[ID]/ delta_time), Damptr->DamMeanQ[ID]);
		NUMERIC_TYPE lower_limit = getmax((Damptr->Volmax[ID] - Damptr->DamVol[ID]), 0) / delta_time;
		Damptr->DamOperationQ[ID] = getmax(lower_limit,upper_limit);
		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID]/delta_time;
		}
	}
	if (Damptr->DamOperationCode[ID] == 3) // Based off matlab code provided by Francesca Pianosi for releases kinearly proportional to of storage
	{
		NUMERIC_TYPE tmp = Damptr->DamMeanQ[ID] * ((Damptr->DamVol[ID] / Damptr->Volmax[ID]) + C(0.5)); // The Qmean flow is associated with 50% storage, tmp varies from 50% to 150% of Qmean
		NUMERIC_TYPE upper_limit = getmin((Damptr->Volmax[ID] / delta_time), tmp);
		NUMERIC_TYPE lower_limit = getmax((Damptr->DamVol[ID] - Damptr->Volmax[ID]), 0) / delta_time;
		Damptr->DamOperationQ[ID] = getmax(lower_limit, upper_limit);
		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID] / delta_time;
		}
	}
	if (Damptr->DamOperationCode[ID] == 4) // Based off paper by Doll et al (2003). Used in WaterGAP
	{
		// Q= kS[S/Smax]^1.5, where k=0.01/d
		Damptr->DamOperationQ[ID] = (C(0.01) * Damptr->DamVol[ID] * pow((Damptr->DamVol[ID] / Damptr->Volmax[ID]), C(1.5))) / delta_time;
		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID] / delta_time;
		}
	}
	if (Damptr->DamOperationCode[ID] == 5) // Based of Wada et al (2014). Used in PCR-GLOBWB
		//Q=(S-Smin)/(Smax-Smin)*Qmean; Smin=10% of Smax
	{
		Damptr->DamOperationQ[ID] = (((Damptr->DamVol[ID] - (Damptr->Volmax[ID] * C(0.2))) / (Damptr->Volmax[ID] * C(0.7))))*Damptr->DamMeanQ[ID]; // /delta_time;
		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID] / delta_time;
		}
	}
	if (Damptr->DamOperationCode[ID] == 6) // Based of Wisser et al (2010). Used in WBMplus
	{
		NUMERIC_TYPE Kappa = C(0.16);
		NUMERIC_TYPE Lambda = C(0.6);

		Damptr->DamOperationQ[ID] = Kappa*(Damptr->DamVin[ID] / delta_time);

		if ((Damptr->DamVin[ID] / delta_time) < Damptr->DamMeanQ[ID])
		{
			Damptr->DamOperationQ[ID] = Lambda*(Damptr->DamVin[ID]/delta_time) + (Damptr->DamMeanQ[ID] - (Damptr->DamVin[ID]/delta_time));
		}

		if (Damptr->DamOperationQ[ID] * delta_time > Damptr->DamVol[ID])
		{
			Damptr->DamOperationQ[ID] = Damptr->DamVol[ID] / delta_time;
		}
	}
	if (Damptr->DamOperationCode[ID] == 7) // To be coded to follow the procedure of Hanasaki et al (2005). Used in HO8
	{
		NUMERIC_TYPE Alpha = C(0.85);
		NUMERIC_TYPE Storage = Damptr->Volmax[ID] / Damptr->DamMeanQ[ID];
		// Yearly Release
		// Do not need to caluclate yearly release only need to calculate Kappa which is the storage at beginning of year divided by alpha times max storage
		if (curr_time <= Damptr->DamYear[ID])
		{
			Damptr->OP7_Kappa[ID] = Damptr->DamVol[ID] / (Alpha*Damptr->Volmax[ID]);
			Damptr->DamYear[ID] += (C(365.0) * C(86400.0)); // corrected from 356 to 365, .0 added to prevent visual studio 2015 error
		}
		// Daily Release
		NUMERIC_TYPE DayRelease = Damptr->DamMeanQ[ID];

		// Switch
		Damptr->DamOperationQ[ID] = Damptr->OP7_Kappa[ID] *DayRelease;
		if (Storage < C(0.5))
		{
			Damptr->DamOperationQ[ID] = pow((Storage/0.5), 2)*Damptr->OP7_Kappa[ID] * DayRelease + (1 - pow((Storage/0.5), 2))*(Damptr->DamVin[ID] / delta_time);
		}

	}
}

void SGC2_UpdateDamFlowVolume(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE * h_grid, NUMERIC_TYPE * volume_grid,
	NUMERIC_TYPE * Qx_grid, NUMERIC_TYPE * Qy_grid, DamData *Damptr, const SGCprams *SGCptr, const NUMERIC_TYPE g, const NUMERIC_TYPE max_Froude)
{
	
	delete[] Damptr->DamVin; // Needed to avoid memory leakage
	Damptr->DamVin = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamMaxH = C(0.0);
	
	for (int i = 0; i < Damptr->TotalEdge; i++)
	{
	//	printf("Value of i: %d\n", i);
		fflush(stdout);
		const int grid_index_this = Damptr->DynamicEdge[2 * i];
		const int reservoir_index = Damptr->DynamicEdge[(2 * i) + 1];
		const int gr = Damptr->DynamicEdge[(2 * i) + 2];
		
		// Floodplain flowing into Dam
				//Calculate Q
				NUMERIC_TYPE h0 = h_grid[grid_index_this];
				
				NUMERIC_TYPE h1 = Damptr->InitialHeight[reservoir_index - 1];  //h1 is absolute height - not relative height(equivilent to h1 + dem[index])
				NUMERIC_TYPE z0 = Damptr->DynamicEdgeData[12 * i + 4];
				NUMERIC_TYPE zb0 = Damptr->DynamicEdgeData[12 * i + 6];
				//	NUMERIC_TYPE z1 = C(0.0); //Dam has infinite Depth NOT NEEDED!!!
				NUMERIC_TYPE friction = Damptr->DynamicEdgeData[12 * i + 3] * Damptr->DynamicEdgeData[12 * i + 3] * g;
				NUMERIC_TYPE friction1= Damptr->DynamicEdgeData[12 * i + 5] * Damptr->DynamicEdgeData[12 * i + 5] * g;
				NUMERIC_TYPE surface_elevation0 = z0 + h0;
				NUMERIC_TYPE surface_elevation1 = h1;
				NUMERIC_TYPE SGC_width = Damptr->DynamicEdgeData[12 * i + 11]; // SGC width
				NUMERIC_TYPE hflow = C(0.0);
				
				
				//
			
				Damptr->DynamicEdgeData[12 * i + 1] = C(0.0);
				NUMERIC_TYPE dh, surface_slope, A, R, deltaV;
					if (SGC_width > C(0.0)) //  check for sub-grid channel
					{
						hflow = getmax(surface_elevation0, surface_elevation1) - zb0;
						if (hflow > depth_thresh)
						{
							dh = surface_elevation1 - surface_elevation0;//, surface_elevation0 - Damptr->DynamicEdgeData[12 * i + 6]);
							surface_slope = dh / Damptr->DynamicEdgeData[12 * i + 8];
							//surface_slope = (-1 * abs(surface_slope));
							SGC2_CalcA(gr, hflow, Damptr->DynamicEdgeData[12 * i + 6], &A, &SGC_width, SGCptr); // calculate channel area for SGC
							R = SGC2_CalcR(gr, hflow, Damptr->DynamicEdgeData[12 * i + 6], SGC_width, SGC_width, A, SGCptr); // calculate hydraulic radius for SGC

							Damptr->DynamicEdgeData[12 * i + 7] = CalculateQ(surface_slope, R, delta_time, g, A, friction1, Damptr->DynamicEdgeData[12 * i + 7], max_Froude);
						}	// DynamicEdgeData{12i +7] = subgrid flow m3/s
						else Damptr->DynamicEdgeData[12 * i + 7] = C(0.0);
					}

					//hflow = getmax(getmax(surface_elevation1 - z0, h0), C(0.0));
					hflow = getmax(surface_elevation0, surface_elevation1) - z0;
					if (hflow > depth_thresh && SGC_width < Damptr->DynamicEdgeData[12 * i + 2])
					{
						dh = surface_elevation1 - surface_elevation0;
						surface_slope = dh / Damptr->DynamicEdgeData[12 * i + 8];
					//	surface_slope = (-1 * abs(surface_slope));
						// Set Cross-Section area equal to the min of cell size or flow cross section
						A = min(Damptr->DynamicEdgeData[12 * i + 2],Damptr->DynamicEdgeData[12*i+8]) * hflow;
						// calculate FP flow
						NUMERIC_TYPE Q;
						Q = CalculateQ(surface_slope, hflow, delta_time, g, A, friction, min(Damptr->DynamicEdgeData[12 * i + 2], Damptr->DynamicEdgeData[12 * i + 8]), max_Froude);
						Damptr->DynamicEdgeData[12 * i] = Q; // FP Flow in m3/s

						if (SGC_width > C(0.0))
						{
							NUMERIC_TYPE channel_ratio = min(SGC_width / min(Damptr->DynamicEdgeData[12 * i + 2], Damptr->DynamicEdgeData[12 * i + 8]), C(1.0));
							Q = Q - channel_ratio * Q;
						}
						
						Damptr->DynamicEdgeData[12*i+1] = Q; // DynamicEdgeData[12i+1] = TotalQ m3/s
					}
					else Damptr->DynamicEdgeData[12 * i] = C(0.0); //FP Q set to 0. 

					Damptr->DynamicEdgeData[12 * i + 1] += Damptr->DynamicEdgeData[12 * i + 7]; // Add Sub-grid Q to Total Q. m3/s

					deltaV = (Damptr->DynamicEdgeData[12 * i + 1] * delta_time);
					
					// Check Statement needed to ensure volume grid >= zero
					NUMERIC_TYPE DVratio = min(volume_grid[grid_index_this]/deltaV, C(1.0));;
					if (deltaV < C(0.0))
					{
						DVratio = C(1.0);
					
					}
					DVratio = C(1.0);
					volume_grid[grid_index_this] -= deltaV*DVratio;
					Damptr->DamVin[reservoir_index - 1] += deltaV*DVratio;
					Damptr->DynamicEdgeData[12 * i] = Damptr->DynamicEdgeData[12 * i] *DVratio;
					Damptr->DynamicEdgeData[12 * 7] = Damptr->DynamicEdgeData[12 * 7] *DVratio;

					Damptr->DamMaxH = getmax(Damptr->DamMaxH, hflow);
			}
	//		
	
	// Calculate Update Volume of Dam
	//NUMERIC_TYPE DamVol
	
	for (int d = 0; d < Damptr->NumDams; d++)
	{
		
		Damptr->DamOperationQ[d] = C(0.0);
		Damptr->DamTotalQ[d] = C(0.0);
		// Operational Outflows will go here;
		
		DamOpQ(delta_time,curr_time, h_grid, Qx_grid, Qy_grid, Damptr, g, d);
			Damptr->DamTotalQ[d] += Damptr->DamOperationQ[d];

		// Spill Way Calcs.
		if (Damptr->InitialHeight[d] <= Damptr->SpillHeight[d])
		{
			Damptr->SpillQ[d] = C(0.0);
		}
		else Damptr->SpillQ[d] = Damptr->Spill_Cd[d] * Damptr->SpillWidth[d] *pow(g,(C(0.5))) *pow((Damptr->InitialHeight[d] - Damptr->SpillHeight[d]), (C(1.5)));

		Damptr->DamTotalQ[d] += Damptr->SpillQ[d];

		// Volume and H update for Dams
		//		Mass balance for Dam is the change in reservoir storage
		Damptr->DamLoss = Damptr->DamLoss + Damptr->DamVol[d];
		// Update Reservoir Volume with the volume change in boundary cells
		Damptr->DamVol[d] += Damptr->DamVin[d];
		// Remove outflow from Reservoir Volume and Update DamLoss calculation
		Damptr->DamVol[d] -= (Damptr->DamTotalQ[d] * delta_time);
		Damptr->DamLoss = Damptr->DamLoss - Damptr->DamVol[d];
		//	Update Height of Reservoir
		Damptr->InitialHeight[d] = (Damptr->DamVol[d] / Damptr->Volmax[d])*Damptr->DamHeight[d] + (Damptr->SpillHeight[d] - Damptr->DamHeight[d]);

			// Output DamTotalQ to output cells
			//convert coordinates to grid cell numbers already done in input.cpp
		// convert to padded_cell_index
			int output = Damptr->OutputCellX[d] + Damptr->OutputCellY[d] * grid_cols_padded;
			volume_grid[output] += (Damptr->DamTotalQ[d] * delta_time);
	}
	
}

// Dam Code FEOL 2016 



/*
 * This function can be used to calculate the diffusive switch.
 * Can be used to handle either the Qx or Qy direction.
 * Qx: next_cell_add = 1
 * Qy: next_cell_add = grid_cols_padded
 * use appropriate Q_grid, Q_old_grid and g_friction_sq
 */
inline NUMERIC_TYPE SGC2_UpdateDiffusiveQ_row(const int grid_row_index, const int row_start, const int row_end,
	const NUMERIC_TYPE row_cell_length, const NUMERIC_TYPE row_cell_width,
	const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE froude_thresh, const NUMERIC_TYPE diffusive_max_hflow,
	const NUMERIC_TYPE g,
	const int next_cell_add, // qx next_cell_add = 1, qy next_cell_add = grid_cols_padded
	const NUMERIC_TYPE * h_grid,
	const NUMERIC_TYPE * dem_grid,
	const NUMERIC_TYPE * friction_grid,
	NUMERIC_TYPE * Q_grid, NUMERIC_TYPE * Q_old_grid)
{
#ifdef __INTEL_COMPILER
	__assume_aligned(h_grid, 64);
	__assume_aligned(dem_grid, 64);
	__assume_aligned(Q_grid, 64);
	__assume_aligned(Q_old_grid, 64);
	__assume_aligned(friction_grid, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

	int route_count = 0;
	for (int i = row_start; i < row_end; i++)
	{
		int index = grid_row_index + i;
		//next column
		int index_next = index + next_cell_add;

		// if diffusive turns out to be better than routing - could use the tmp_row to store hflow instead of surface slope
		NUMERIC_TYPE h0 = h_grid[index];
		NUMERIC_TYPE h1 = h_grid[index_next];
		NUMERIC_TYPE z0 = dem_grid[index];
		NUMERIC_TYPE z1 = dem_grid[index_next];
		NUMERIC_TYPE surface_elevation0 = z0 + h0;
		NUMERIC_TYPE surface_elevation1 = z1 + h1;
		// Calculating hflow based on floodplain levels
		NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
			- getmax(z0, z1);

		if (hflow > depth_thresh)
		{
			NUMERIC_TYPE velocity = (Q_grid[index_next] / row_cell_width) / hflow;
			//NUMERIC_TYPE velocity = Velocity_grid[index_next];
			NUMERIC_TYPE froude = velocity / SQRT(g / hflow);

			/// hflow limiter - diffusive will be unstable if this is not used // TFD
			hflow = getmin(hflow, diffusive_max_hflow);

			if (FABS(froude) > froude_thresh)
			{
				// could be faster if using the regular 
				NUMERIC_TYPE fn = friction_grid[index_next];

				NUMERIC_TYPE Q;
				if (surface_elevation0 > surface_elevation1 && h0 > depth_thresh)
				{
					NUMERIC_TYPE dh = surface_elevation0 - surface_elevation1;
					NUMERIC_TYPE Sf = SQRT(dh / row_cell_length);
					Q = (POW(hflow, (C(5.0) / C(3.0)))*Sf*row_cell_width / fn);
				}
				else if (surface_elevation1 > surface_elevation0 && h1 > depth_thresh)
				{
					NUMERIC_TYPE dh = surface_elevation1 - surface_elevation0;
					NUMERIC_TYPE Sf = SQRT(dh / row_cell_length);
					Q = (-POW(hflow, (C(5.0) / C(3.0)))*Sf*row_cell_width / fn);
				}
				else
				{
					Q = C(0.0);
				}

				//printf("old %" NUM_FMT" new %" NUM_FMT"\n", Q_grid[index_next], Q);
				Q_grid[index_next] = Q;

				// if resetting to zero, the inertial model Q will change dramatically 
				// which in turn changes the froude causing jumping back and forward
				//Q_old_grid[index_next] = C(0.0);

				// update the Q_old with the diffusive Q
				// this will be the input to the inertial model
				Q_old_grid[index_next] = Q;

				// don't update Q_old_grid - the old inertial Q will be used
				// Q_old_grid[index_next] = Q_old_grid[index_next];
			}
			/*	else
				{
				printf("old %d %" NUM_FMT"\n",itCount, Q_grid[index_next]);
				}*/
		}
	}
	//printf("%d --\n", itCount);
	return route_count;
}

/*
// qx: next_cell_add = 1,
// qy: next_cell_add = grid_cols_padded
*/

inline NUMERIC_TYPE SGC2_UpdateRouteQ_row(const int grid_row_index, const int row_start, const int row_end,
	const NUMERIC_TYPE row_cell_area, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE route_slope_thresh,
	const int next_cell_add,
	const NUMERIC_TYPE * surface_slopes_row,
	const NUMERIC_TYPE * h_grid,
	const NUMERIC_TYPE * dem_grid,
	const NUMERIC_TYPE * route_V_ratio_per_sec_grid,
	NUMERIC_TYPE * Q_grid, NUMERIC_TYPE * Q_old_grid,
	int * route_list_i_lookup)
{
#ifdef __INTEL_COMPILER
	__assume_aligned(surface_slopes_row, 64);
	__assume_aligned(h_grid, 64);
	__assume_aligned(dem_grid, 64);
	__assume_aligned(route_V_ratio_per_sec_grid, 64);
	__assume_aligned(Q_grid, 64);
	__assume_aligned(Q_old_grid, 64);
	__assume_aligned(route_list_i_lookup, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

	int route_count = 0;
	for (int i = row_start; i < row_end; i++)
	{
		// always turn off the flood plain flow when slope above route_slope_thresh
		if (FABS(surface_slopes_row[i]) > route_slope_thresh)
		{
			int index = grid_row_index + i;
			//next column
			int index_next = index + next_cell_add;

			if (surface_slopes_row[i] > C(0.0) && route_V_ratio_per_sec_grid[index_next] > C(0.0))
			{
				NUMERIC_TYPE q_route = -CalculateRoutingQ(delta_time,
					h_grid[index_next], h_grid[index],
					dem_grid[index_next], dem_grid[index], route_V_ratio_per_sec_grid[index_next], row_cell_area);
#ifdef _DEBUG
				//printf("[%d,%d] Route Q %" NUM_FMT" -> %" NUM_FMT"\n", i, j, Q_grid[index_next], q_route);
#endif
				Q_grid[index_next] = q_route;

				// 3 options for old Q -- Toby Dunne
				// - leave as inertial Q (could be unstable since switched)
				// - set to zero (sudden change could cause instability)
				// - set to the replaced Q ** chosen as this is the actual previous Q used, should result in more stable inertia
				Q_old_grid[index_next] = q_route;
				//Q_old_grid[index_next] = C(0.0);

				route_list_i_lookup[grid_row_index + route_count] = i;
				route_count++;
			}
			else if (surface_slopes_row[i] < C(0.0) && route_V_ratio_per_sec_grid[index_next] < C(0.0))
			{
				NUMERIC_TYPE q_route = CalculateRoutingQ(delta_time,
					h_grid[index], h_grid[index_next],
					dem_grid[index], dem_grid[index_next], -route_V_ratio_per_sec_grid[index_next], row_cell_area);
#ifdef _DEBUG
				//printf("[%d,%d] Route Q %" NUM_FMT" -> %" NUM_FMT"\n", i, j, Q_grid[index_next], q_route);
#endif
				Q_grid[index_next] = q_route;

				// 3 options for old Q -- Toby Dunne
				// - leave as inertial Q (could be unstable since switched)
				// - set to zero (sudden change could cause instability)
				// - set to the replaced Q ** chosen as this is the actual previous Q used, should result in more stable inertia
				Q_old_grid[index_next] = q_route;
				//Q_old_grid[index_next] = C(0.0);

				route_list_i_lookup[grid_row_index + route_count] = i;
				route_count++;
			}
		}
	}
	return route_count;
}

///
/// checks each routing q
/// adds up all the out flow from the source cell
/// if the total outflow is greater than the cell volume, reduce the routing q
/// 
inline void SGC2_CorrectRouteFlow_row(const int j, const int grid_row_index, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE delta_time,
	const RouteDynamicList * route_dynamic_list,
	const NUMERIC_TYPE * volume_grid,
	NUMERIC_TYPE* Qx_grid, NUMERIC_TYPE* Qy_grid)
{
	int route_count;
	route_count = route_dynamic_list->row_route_qx_count[j];
#pragma ivdep
	for (int route_index = 0; route_index < route_count; route_index++)
	{
		int i = route_dynamic_list->route_list_i_lookup_qx[route_index];
		int q_index = grid_row_index + i + 1;
		NUMERIC_TYPE q_route = Qx_grid[q_index];
		int source_cell_grid_index;
		if (q_route > 0) // route q is on the east face (source cell west)
		{
			source_cell_grid_index = grid_row_index + i;
			int flow_west = source_cell_grid_index; // west flow of source cell

			int flow_north = source_cell_grid_index;
			int flow_south = source_cell_grid_index + grid_cols_padded;

			NUMERIC_TYPE shallow_flow_out = C(0.0);
			shallow_flow_out -= (i > 0 && Qx_grid[flow_west] < C(0.0)) ? Qx_grid[flow_west] : C(0.0);
			shallow_flow_out -= getmin(Qy_grid[flow_north], C(0.0));
			shallow_flow_out += getmax(Qy_grid[flow_south], C(0.0));

			NUMERIC_TYPE total_out_flow;
			total_out_flow = shallow_flow_out + q_route; //+ east

			total_out_flow *= delta_time;
			if (total_out_flow > volume_grid[source_cell_grid_index])
			{
				q_route -= shallow_flow_out;
				q_route = getmax(q_route, C(0.0));
				//printf("routeQ %" NUM_FMT" shallow vol out: %" NUM_FMT" route => %" NUM_FMT" \n", Qx_grid[q_index], shallow_flow_out * delta_time, q_route);
				Qx_grid[q_index] = q_route;
			}
		}
		else //if (q_route < 0) // route q is on the west face (source cell east)
		{
			source_cell_grid_index = grid_row_index + i + 1;
			int flow_east = source_cell_grid_index + 1; // east flow of source cell

			int flow_north = source_cell_grid_index;
			int flow_south = source_cell_grid_index + grid_cols_padded;

			NUMERIC_TYPE shallow_flow_out = C(0.0);

			shallow_flow_out += (i < grid_cols - 2 && Qx_grid[flow_east] > C(0.0)) ? Qx_grid[flow_east] : C(0.0);
			shallow_flow_out -= getmin(Qy_grid[flow_north], C(0.0));
			shallow_flow_out += getmax(Qy_grid[flow_south], C(0.0));

			NUMERIC_TYPE total_out_flow;
			total_out_flow = shallow_flow_out - q_route; //- west

			total_out_flow *= delta_time;
			if (total_out_flow > volume_grid[source_cell_grid_index])
			{
				q_route += shallow_flow_out;
				q_route = getmin(q_route, C(0.0));
				//printf("routeQ %" NUM_FMT" shallow vol out: %" NUM_FMT" route => %" NUM_FMT" \n", Qx_grid[q_index], shallow_flow_out * delta_time, q_route);
				Qx_grid[q_index] = q_route;
				//total_out_flow = (shallow_flow_out - q_route) * delta_time;
				//printf("total out %" NUM_FMT" vol: %" NUM_FMT"\n", total_out_flow, volume_grid[source_cell_grid_index]);
			}
		}
	}

	route_count = route_dynamic_list->row_route_qy_count[j];
#pragma ivdep
	for (int route_index = 0; route_index < route_count; route_index++)
	{
		int i = route_dynamic_list->route_list_i_lookup_qy[route_index];
		int q_index = grid_row_index + i + grid_cols_padded;
		NUMERIC_TYPE q_route = Qy_grid[q_index];
		int source_cell_grid_index;
		if (q_route > 0) // route q is on south face (source cell north)
		{
			source_cell_grid_index = grid_row_index + i;
			int flow_north = source_cell_grid_index; // north face of the source cell

			int flow_west = source_cell_grid_index; // west face of source cell
			int flow_east = source_cell_grid_index + 1; // east face of source cell

			NUMERIC_TYPE shallow_flow_out = C(0.0);
			shallow_flow_out -= (j > 0 && Qy_grid[flow_north] < C(0.0)) ? Qy_grid[flow_north] : C(0.0);
			shallow_flow_out -= getmin(Qx_grid[flow_west], C(0.0));
			shallow_flow_out += getmax(Qx_grid[flow_east], C(0.0));

			NUMERIC_TYPE total_out_flow;
			total_out_flow = shallow_flow_out + q_route; //+ south

			total_out_flow *= delta_time;
			if (total_out_flow > volume_grid[source_cell_grid_index])
			{
				q_route -= shallow_flow_out;
				q_route = getmax(q_route, C(0.0));
				//printf("routeQ %" NUM_FMT" shallow: %" NUM_FMT" => %" NUM_FMT" \n", Qy_grid[q_index], shallow_flow_out, q_route);
				Qy_grid[q_index] = q_route;
			}
		}
		else //if (q_route < 0) //route is on north face (source cell south)
		{
			source_cell_grid_index = grid_row_index + i + grid_cols_padded;
			int flow_south = source_cell_grid_index + grid_cols_padded; // south face of the source cell

			int flow_west = source_cell_grid_index; // west face of source cell
			int flow_east = source_cell_grid_index + 1; // east face of source cell

			NUMERIC_TYPE shallow_flow_out = C(0.0);

			shallow_flow_out += (j < grid_rows - 2 && Qy_grid[flow_south] > C(0.0)) ? Qy_grid[flow_south] : C(0.0);
			shallow_flow_out -= getmin(Qx_grid[flow_west], C(0.0));
			shallow_flow_out += getmax(Qx_grid[flow_east], C(0.0));

			NUMERIC_TYPE total_out_flow;
			total_out_flow = shallow_flow_out - q_route; //- north

			total_out_flow *= delta_time;
			if (total_out_flow > volume_grid[source_cell_grid_index])
			{
				q_route += shallow_flow_out;
				q_route = getmin(q_route, C(0.0));
				//printf("routeQ %" NUM_FMT" shallow: %" NUM_FMT" => %" NUM_FMT" \n", Qy_grid[q_index], shallow_flow_out, q_route);
				Qy_grid[q_index] = q_route;
			}
		}
	}
}

inline void ProcessSubGridQBlock(const int block_index, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE g,
	const SubGridRowList * sub_grid_layout, SubGridState * sub_grid_state, const SGCprams * SGCptr,
	const NUMERIC_TYPE * h_grid,
	WetDryRowBound * wet_dry_bounds,
	NUMERIC_TYPE * Qx_grid, NUMERIC_TYPE * Qy_grid,
	NUMERIC_TYPE * Qx_old_grid, NUMERIC_TYPE * Qy_old_grid,
	const NUMERIC_TYPE max_Froude)
{
	const int * sg_pair_grid_index_lookup = sub_grid_layout->flow_info.flow_pair.sg_cell_grid_index_lookup;
	const NUMERIC_TYPE * sg_pair_SGC_BankFullHeight = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_BankFullHeight;
	const NUMERIC_TYPE * sg_pair_SGC_width = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_width;
	const NUMERIC_TYPE * sg_pair_dem = sub_grid_layout->flow_info.flow_pair.sg_cell_dem;
	const int * sg_pair_SGC_group = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_group;

	//const NUMERIC_TYPE * sg_flow_cell_width = sub_grid_layout->flow_info.sg_flow_cell_width;
	NUMERIC_TYPE * sg_flow_Q = sub_grid_state->sg_flow_Q;
	//NUMERIC_TYPE * sg_flow_ChannelRatio = sub_grid_state->sg_flow_ChannelRatio;
	const NUMERIC_TYPE * sg_flow_effective_distance = sub_grid_layout->flow_info.sg_flow_effective_distance;
	const NUMERIC_TYPE * sg_flow_g_friction_sq = sub_grid_layout->flow_info.sg_flow_g_friction_sq;

	const int sg_row_start_index = block_index * sub_grid_layout->row_cols_padded;
	const int sg_row_pair_start_index = block_index * 2 * sub_grid_layout->row_cols_padded;
	const int flow_end = sub_grid_layout->flow_row_count[block_index];
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(sg_row_start_index % GRID_ALIGN_WIDTH == 0);
	__assume(sg_row_pair_start_index % GRID_ALIGN_WIDTH == 0);
#endif
#ifdef __INTEL_COMPILER
	__assume_aligned(sg_pair_SGC_BankFullHeight, 64);
	__assume_aligned(sg_pair_SGC_width, 64);
	__assume_aligned(sg_pair_dem, 64);
	__assume_aligned(sg_pair_grid_index_lookup, 64);
	__assume_aligned(h_grid, 64);
	__assume_aligned(Qx_grid, 64);
	__assume_aligned(Qy_grid, 64);
	__assume_aligned(Qx_old_grid, 64);
	__assume_aligned(Qy_old_grid, 64);
#endif


#pragma ivdep
#pragma simd
	for (int flow_i = 0; flow_i < flow_end; flow_i++)
	{
		int flow_index = sg_row_start_index + flow_i;
		int flow_pair_index = sg_row_pair_start_index + 2 * flow_i;
		int flow_pair_index_next = flow_pair_index + 1;

		int grid_index0 = sg_pair_grid_index_lookup[flow_pair_index];
		int grid_index1 = sg_pair_grid_index_lookup[flow_pair_index_next]; // also the q index for this flow (in d4)

		NUMERIC_TYPE h0 = h_grid[grid_index0];
		NUMERIC_TYPE h1 = h_grid[grid_index1];

		NUMERIC_TYPE channel_Q = C(0.0);
		//NUMERIC_TYPE channel_ratio = C(0.0);

		NUMERIC_TYPE bankFullHeight0 = sg_pair_SGC_BankFullHeight[flow_pair_index];
		NUMERIC_TYPE bankFullHeight1 = sg_pair_SGC_BankFullHeight[flow_pair_index_next];

		if (h0 + bankFullHeight0 > depth_thresh ||
			h1 + bankFullHeight1 > depth_thresh)
		{
			//NUMERIC_TYPE cell_width = sg_flow_cell_width[flow_index];
			//NUMERIC_TYPE channel_area0, channel_area1;
			NUMERIC_TYPE channel_width0 = sg_pair_SGC_width[flow_pair_index];
			NUMERIC_TYPE channel_width1 = sg_pair_SGC_width[flow_pair_index_next];

			NUMERIC_TYPE z0 = sg_pair_dem[flow_pair_index];
			NUMERIC_TYPE z1 = sg_pair_dem[flow_pair_index_next];

			NUMERIC_TYPE surface_elevation0 = z0 + h0;
			NUMERIC_TYPE surface_elevation1 = z1 + h1;
			// Calculating hflow based on sub-channel elevation levels
			NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
				- getmax(z0 - bankFullHeight0, z1 - bankFullHeight1);
			if (hflow > depth_thresh)
			{
				int channel_group0 = sg_pair_SGC_group[flow_pair_index];
				int channel_group1 = sg_pair_SGC_group[flow_pair_index_next];

				// SGC2_CalcA(channel_group0, hflow, bankFullHeight0, &channel_area0, &channel_width0, SGCptr);
				// SGC2_CalcA(channel_group1, hflow, bankFullHeight1, &channel_area1, &channel_width1, SGCptr);

				NUMERIC_TYPE area;
				NUMERIC_TYPE R;

				// always use smallest flow area
				// PFU, change to using smallest channel width
				if (channel_width0 < channel_width1) // select the appropriate slope and area based on the smallest area
				{
					// calculate hydraulic radius for SGC
					SGC2_CalcA(channel_group0, hflow, bankFullHeight0, &area, &channel_width0, SGCptr);
					R = SGC2_CalcR(channel_group0, hflow, bankFullHeight0, channel_width0, sg_pair_SGC_width[flow_pair_index], area, SGCptr);

				}
				else
				{
					// calculate hydraulic radius for SGC
					SGC2_CalcA(channel_group1, hflow, bankFullHeight1, &area, &channel_width1, SGCptr);
					R = SGC2_CalcR(channel_group1, hflow, bankFullHeight1, channel_width1, sg_pair_SGC_width[flow_pair_index_next], area, SGCptr);
				}

				NUMERIC_TYPE effective_distance = sg_flow_effective_distance[flow_index];
				NUMERIC_TYPE g_friction_squared = sg_flow_g_friction_sq[flow_index];

				NUMERIC_TYPE dh = (surface_elevation0)-(surface_elevation1);
				NUMERIC_TYPE surface_slope = -dh / effective_distance;
				channel_Q = CalculateQ(surface_slope, R, delta_time, g, area, g_friction_squared, sg_flow_Q[flow_index], max_Froude);
			}
		}
		sg_flow_Q[flow_index] = channel_Q;
		//sg_flow_ChannelRatio[flow_index] = getmin(channel_ratio, C(1.0)); //PFU set constant channel ratio in lisflood_processing

		// Update wet-dry 
		if (channel_Q != C(0.0))
		{
			int x = sub_grid_layout->flow_info.flow_pair.sg_cell_x[flow_pair_index];
			int y = sub_grid_layout->flow_info.flow_pair.sg_cell_y[flow_pair_index];
			wet_dry_bounds->fp_vol[y].start = min(wet_dry_bounds->fp_vol[y].start, x);
			wet_dry_bounds->fp_vol[y].end = max(wet_dry_bounds->fp_vol[y].end, x + 1);
		}
	}
}

void SGC2_UpdateVelocitySubGrid_block(const int block_index, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE delta_time,
	const SubGridRowList * sub_grid_layout, SubGridState * sub_grid_state, const SGCprams * SGCptr,
	const NUMERIC_TYPE * h_grid)
{
	const int * sg_pair_grid_index_lookup = sub_grid_layout->flow_info.flow_pair.sg_cell_grid_index_lookup;
	const NUMERIC_TYPE * sg_pair_dem = sub_grid_layout->flow_info.flow_pair.sg_cell_dem;
	const NUMERIC_TYPE * sg_pair_SGC_BankFullHeight = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_BankFullHeight;
	const NUMERIC_TYPE * sg_pair_SGC_width = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_width;
	const int * sg_pair_SGC_group = sub_grid_layout->flow_info.flow_pair.sg_cell_SGC_group;

	NUMERIC_TYPE * sg_flow_Q = sub_grid_state->sg_flow_Q;
	NUMERIC_TYPE * sg_velocity = sub_grid_state->sg_velocity;

	const int sg_row_start_index = block_index * sub_grid_layout->row_cols_padded;
	const int sg_row_pair_start_index = block_index * 2 * sub_grid_layout->row_cols_padded;
	const int flow_end = sub_grid_layout->flow_row_count[block_index];

#pragma ivdep
#pragma simd
	for (int flow_i = 0; flow_i < flow_end; flow_i++)
	{
		int flow_index = sg_row_start_index + flow_i;

		NUMERIC_TYPE q = sg_flow_Q[flow_index];
		if (FABS(q) > depth_thresh)
		{
			int flow_pair_index = sg_row_pair_start_index + 2 * flow_i;
			int flow_pair_index_next = flow_pair_index + 1;

			int grid_index0 = sg_pair_grid_index_lookup[flow_pair_index];
			int grid_index1 = sg_pair_grid_index_lookup[flow_pair_index_next]; // also the q index for this flow (in d4)

			NUMERIC_TYPE h0 = h_grid[grid_index0];
			NUMERIC_TYPE h1 = h_grid[grid_index1];
			NUMERIC_TYPE z0 = sg_pair_dem[flow_pair_index];
			NUMERIC_TYPE z1 = sg_pair_dem[flow_pair_index_next];
			NUMERIC_TYPE surface_elevation0 = z0 + h0;
			NUMERIC_TYPE surface_elevation1 = z1 + h1;
			NUMERIC_TYPE bankFullHeight0 = sg_pair_SGC_BankFullHeight[flow_pair_index];
			NUMERIC_TYPE bankFullHeight1 = sg_pair_SGC_BankFullHeight[flow_pair_index_next];

			// Calculating hflow based on sub-channel elevation levels
			NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
				- getmax(z0 - bankFullHeight0, z1 - bankFullHeight1);
			if (hflow > depth_thresh)
			{
				NUMERIC_TYPE channel_area0, channel_area1;
				NUMERIC_TYPE channel_width0 = sg_pair_SGC_width[flow_pair_index];
				NUMERIC_TYPE channel_width1 = sg_pair_SGC_width[flow_pair_index_next];

				int channel_group0 = sg_pair_SGC_group[flow_pair_index];
				int channel_group1 = sg_pair_SGC_group[flow_pair_index_next];

				SGC2_CalcA(channel_group0, hflow, bankFullHeight0, &channel_area0, &channel_width0, SGCptr);
				SGC2_CalcA(channel_group1, hflow, bankFullHeight1, &channel_area1, &channel_width1, SGCptr);

				NUMERIC_TYPE channel_area = getmin(channel_area0, channel_area1);

				sg_velocity[flow_index] = q / channel_area;
			}
		}
	}
}

inline void SGC2_UpdateQx_row(const int grid_cols,
	const int grid_row_index,
	const int row_start_x_prev, const int row_end_x_prev,
	const int row_start_x, const int row_end_x,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE row_dx, const NUMERIC_TYPE * Fp_ywidth,
	const NUMERIC_TYPE g, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time,
	NUMERIC_TYPE * tmp_row,
	const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE *h_grid,

	const NUMERIC_TYPE * g_friction_sq_x_grid,
	NUMERIC_TYPE *Qx_grid, NUMERIC_TYPE *Qx_old_grid,
	const NUMERIC_TYPE max_Froude)
{
#ifdef _DEBUG

	// checking only
	int check_end = min(grid_cols - 1, row_start_x_prev);
	for (int i = 0; i < check_end; i++)
	{
		int index = grid_row_index + i;
		int index_next = index + 1;
		if (Qx_old_grid[index_next] != C(0.0))
			printf("Error: Qx %" NUM_FMT" @ %d (%d)  \n", Qx_old_grid[index_next], index_next, i);
	}
	// checking only
	if (row_end_x_prev != -1)
		for (int i = row_end_x_prev; i < grid_cols - 1; i++)
		{
			int index = grid_row_index + i;
			int index_next = index + 1;
			if (Qx_old_grid[index_next] != C(0.0))
				printf("Error: Qx %" NUM_FMT" @ %d (%d)  \n", Qx_old_grid[index_next], index_next, i);
		}

#endif

#ifdef __INTEL_COMPILER
	__assume_aligned(h_grid, 64);
	__assume_aligned(dem_grid, 64);
	__assume_aligned(g_friction_sq_x_grid, 64);
	__assume_aligned(Qx_grid, 64);
	__assume_aligned(Qx_old_grid, 64);
	__assume_aligned(Fp_ywidth, 64);
	__assume_aligned(tmp_row, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

#pragma ivdep
#pragma simd
	// clear from start of prev bound to start of new bound
	for (int i = row_start_x_prev; i < row_start_x; i++)
	{
		int index = grid_row_index + i;
		//next row
		int index_next = index + 1;
		Qx_grid[index_next] = C(0.0);
		Qx_old_grid[index_next] = C(0.0);
	}
	if (row_end_x != -1)
	{
#pragma ivdep
#pragma simd
		// clear from end bound to end of prev bound
		for (int i = row_end_x; i < row_end_x_prev; i++)
		{
			int index = grid_row_index + i;
			//next row
			int index_next = index + 1;
			Qx_grid[index_next] = C(0.0);
			Qx_old_grid[index_next] = C(0.0);
		}
	}

	// Calculate Qx (base model)
#pragma ivdep
#pragma simd // note this pragma is here to hint the compiler that this should be vectorized - the compiler will warn if not vectorized
	for (int i = row_start_x; i < row_end_x; i++)
	{
		int index = grid_row_index + i;
		//next column
		int index_next = index + 1;

		// -ve h where sub grid channel
		//NUMERIC_TYPE h0 = getmax(h_grid[index], C(0.));
		//NUMERIC_TYPE h1 = getmax(h_grid[index_next], C(0.));
		NUMERIC_TYPE h0 = h_grid[index];
		NUMERIC_TYPE h1 = h_grid[index_next];
		NUMERIC_TYPE z0 = dem_grid[index];
		NUMERIC_TYPE z1 = dem_grid[index_next];

		NUMERIC_TYPE surface_elevation0 = z0 + h0;
		NUMERIC_TYPE surface_elevation1 = z1 + h1;
		// Calculating hflow based on floodplain levels
		NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
			- getmax(z0, z1);
		NUMERIC_TYPE q_tmp, surface_slope;
		if (hflow > depth_thresh)// && (h0 > depth_thresh || h1 > depth_thresh))
		{
			//NUMERIC_TYPE area = (row_dy)* hflow;
			// PFU use floodplain width corrected for sub grid channel ratio rather than cell width
			NUMERIC_TYPE area = Fp_ywidth[index_next]* hflow;
			NUMERIC_TYPE dh = (surface_elevation0)-(surface_elevation1);
			surface_slope = -dh / row_dx;
			q_tmp = CalculateQ(surface_slope, hflow, delta_time, g, area, g_friction_sq_x_grid[index_next], Qx_old_grid[index_next], max_Froude);
		}
		else
		{
			surface_slope = C(0.0);
			q_tmp = C(0.0);
		}
		Qx_old_grid[index_next] = q_tmp;
		tmp_row[i] = surface_slope;
	}
	int count = row_end_x - row_start_x;
	if (count > 0)
		memcpy(Qx_grid + grid_row_index + row_start_x + 1, Qx_old_grid + grid_row_index + row_start_x + 1, sizeof(NUMERIC_TYPE) * count);
}

inline void SGC2_UpdateQy_row(const int grid_cols,
	const int grid_row_index,
	const int row_start_y_prev, const int row_end_y_prev,
	const int row_start_y, const int row_end_y,
	const int next_cell_add,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE * Fp_xwidth, const NUMERIC_TYPE row_dy,
	const NUMERIC_TYPE g, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time,
	NUMERIC_TYPE * tmp_row,
	const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE *h_grid,

	const NUMERIC_TYPE * g_friction_sq_y_grid,
	NUMERIC_TYPE *Qy_grid, NUMERIC_TYPE *Qy_old_grid,
	const NUMERIC_TYPE max_Froude)
{
#ifdef _DEBUG
	// checking only
	for (int i = 0; i < row_start_y_prev; i++)
	{
		int index = grid_row_index + i;
		int index_next = index + next_cell_add;
		if (Qy_old_grid[index_next] != C(0.0))
			printf("Error: Qy %" NUM_FMT" @ %d (%d)  \n", Qy_old_grid[index_next], index_next, i);
	}
	// checking only
	if (row_end_y_prev != -1)
		for (int i = row_end_y_prev; i < grid_cols; i++)
		{
			int index = grid_row_index + i;
			int index_next = index + next_cell_add;
			if (Qy_old_grid[index_next] != C(0.0))
				printf("Error: Qy %" NUM_FMT" @ %d (%d)  \n", Qy_old_grid[index_next], index_next, i);
		}
#endif

#ifdef __INTEL_COMPILER
	__assume_aligned(h_grid, 64);
	__assume_aligned(dem_grid, 64);
	__assume_aligned(g_friction_sq_y_grid, 64);
	__assume_aligned(Qy_grid, 64);
	__assume_aligned(Qy_old_grid, 64);
	__assume_aligned(Fp_xwidth, 64);
	__assume_aligned(tmp_row, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
	__assume(next_cell_add % GRID_ALIGN_WIDTH == 0);
#endif

#pragma ivdep
#pragma simd
	// clear from start of prev bound to start of new bound
	for (int i = row_start_y_prev; i < row_start_y; i++)
	{
		//next row
		int index_next = grid_row_index + i + next_cell_add;
		Qy_grid[index_next] = C(0.0);
		Qy_old_grid[index_next] = C(0.0);
	}
	if (row_end_y != -1)
	{
#pragma ivdep
#pragma simd
		// clear from end bound to end of prev bound
		for (int i = row_end_y; i < row_end_y_prev; i++)
		{
			//next row
			int index_next = grid_row_index + i + next_cell_add;
			Qy_grid[index_next] = C(0.0);
			Qy_old_grid[index_next] = C(0.0);
		}
	}
	// Calculate Qy (base model)
#pragma ivdep
#pragma simd
	for (int i = row_start_y; i < row_end_y; i++)
	{
		int index = grid_row_index + i;
		//next row
		int index_next = index + next_cell_add;

		//NUMERIC_TYPE h0 = getmax(h_grid[index], C(0.));
		//NUMERIC_TYPE h1 = getmax(h_grid[index_next], C(0.));
		NUMERIC_TYPE h0 = h_grid[index];
		NUMERIC_TYPE h1 = h_grid[index_next];
		NUMERIC_TYPE z0 = dem_grid[index];
		NUMERIC_TYPE z1 = dem_grid[index_next];

		NUMERIC_TYPE surface_elevation0 = z0 + h0;
		NUMERIC_TYPE surface_elevation1 = z1 + h1;
		// Calculating hflow based on floodplain levels
		NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1)
			- getmax(z0, z1);
		NUMERIC_TYPE q_tmp, surface_slope;
		if (hflow > depth_thresh)// && (h0 > depth_thresh || h1 > depth_thresh))
		{
			//NUMERIC_TYPE area = (row_dx)* hflow;
			// PFU use floodplain width corrected for sub grid channel ratio rather than cell width
			NUMERIC_TYPE area = Fp_xwidth[index_next]* hflow;
			NUMERIC_TYPE dh = (surface_elevation0)-(surface_elevation1);
			surface_slope = -dh / row_dy;
			q_tmp = CalculateQ(surface_slope, hflow, delta_time, g, area, g_friction_sq_y_grid[index_next], Qy_old_grid[index_next], max_Froude);
		}
		else
		{
			surface_slope = C(0.0);
			q_tmp = C(0.0);
		}
		Qy_old_grid[index_next] = q_tmp;
		tmp_row[i] = surface_slope;
	}
	int count = row_end_y - row_start_y;
	if (count > 0)
		memcpy(Qy_grid + grid_row_index + row_start_y + next_cell_add, Qy_old_grid + grid_row_index + row_start_y + next_cell_add, sizeof(NUMERIC_TYPE) * count);
}

	//-----------------------------------------------------------------------------
	// FLOODPLAIN DISTRIBUTED INFILTRATION
	// with correction for sub grid channels
	inline NUMERIC_TYPE SGC2_Infil_floodplain_row(
		const int row_start, int row_end,
		const NUMERIC_TYPE depth_thresh,
		const NUMERIC_TYPE row_cell_area,
		const NUMERIC_TYPE evap_deltaH_step,
		const NUMERIC_TYPE * h_row,
		const NUMERIC_TYPE * infil_row,
		NUMERIC_TYPE * volume_row)
	{
#ifdef __INTEL_COMPILER
		__assume_aligned(h_row, 64);
		__assume_aligned(volume_row, 64);
#endif

		NUMERIC_TYPE reduce_evap_loss = C(0.0);
#pragma ivdep
#pragma simd
		for (int i = row_start; i < row_end; i++)
		{
			NUMERIC_TYPE h_new, dV = C(0.0);
			NUMERIC_TYPE h_old = h_row[i];
			NUMERIC_TYPE evap_deltaV_step = infil_row[i] * row_cell_area;
			if (h_old > depth_thresh) // There is water to evaporate on the flood plain
			{
				// update depth by subtracting evap depth
				h_new = h_old - infil_row[i];
				//check for -ve depths
				if (h_new < C(0.0))
				{
					// reduce evap loss to account for dry bed (don't go below 0)
					dV = h_old * row_cell_area;
				}
				else
				{
					dV = evap_deltaV_step;
				}
				volume_row[i] -= dV;
			}
			reduce_evap_loss += dV; //mass-balance for a standard cell
		}
		return reduce_evap_loss;
	}

//-----------------------------------------------------------------------------
// FLOODPLAIN EVAPORATION
// with correction for sub grid channels
inline NUMERIC_TYPE SGC2_Evaporation_floodplain_row(
	const int row_start, int row_end,
	const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE row_cell_area,
	const NUMERIC_TYPE evap_deltaH_step,
	const NUMERIC_TYPE * evap_row,
	const NUMERIC_TYPE * h_row,
	NUMERIC_TYPE * volume_row)
{
#ifdef __INTEL_COMPILER
	__assume_aligned(h_row, 64);
	__assume_aligned(volume_row, 64);
#endif

	NUMERIC_TYPE reduce_evap_loss = C(0.0);
#pragma ivdep
#pragma simd
	for (int i = row_start; i < row_end; i++)
	{
		NUMERIC_TYPE h_new, dV = C(0.0);
		NUMERIC_TYPE h_old = h_row[i];
		NUMERIC_TYPE evap_deltaV_step = evap_row[i] * row_cell_area;
		if (h_old > depth_thresh) // There is water to evaporate on the flood plain
		{
			// update depth by subtracting evap depth
			h_new = h_old - evap_row[i];
			//check for -ve depths
			if (h_new < C(0.0))
			{
				// reduce evap loss to account for dry bed (don't go below 0)
				dV = h_old * row_cell_area;
			}
			else
			{
				dV = evap_deltaV_step;
			}
			volume_row[i] -= dV;
		}
		reduce_evap_loss += dV; //mass-balance for a standard cell
	}
	return reduce_evap_loss;
}

// routine for uniform rainfall
inline NUMERIC_TYPE SGC2_Uniform_Rainfall_row(const int j,
	const NUMERIC_TYPE rain_deltaV_step,
	const NUMERIC_TYPE * dem_row,
	NUMERIC_TYPE * volume_row,
	WetDryRowBound * wet_dry_bounds, const States *Statesptr )
{
	NUMERIC_TYPE loc_rainfall_total = C(0.0);

	const int row_start = wet_dry_bounds->dem_data[j].start;
	const int row_end = wet_dry_bounds->dem_data[j].end;

	// update wet_dry_bounds as all dem cells will now be wet
	wet_dry_bounds->fp_vol[j] = wet_dry_bounds->dem_data[j];

#ifdef __INTEL_COMPILER
	__assume_aligned(dem_row, 64);
	__assume_aligned(volume_row, 64);
#endif

#pragma ivdep
#pragma simd reduction (+:loc_rainfall_total)
	
	for (int i = row_start; i < row_end; i++)
	{
		//NUMERIC_TYPE dV = (dem_row[i] != C(1e10)) ? rain_step_dV : C(0.0);
		NUMERIC_TYPE dV;
		if (dem_row[i] != DEM_NO_DATA)
		{
			dV = rain_deltaV_step;
			volume_row[i] += dV; // add rainfall volume to cell		
			loc_rainfall_total += dV; // mass balance for local cell (cumulative)
		}
	}
	return loc_rainfall_total;
}

// routine for distributed and distributed time varying rainfall
inline NUMERIC_TYPE SGC2_Distrubuted_Rainfall_row(const int j,
	const NUMERIC_TYPE delta_time,
	const NUMERIC_TYPE * dem_row, const NUMERIC_TYPE * rainmask_row,
	NUMERIC_TYPE * volume_row,
	WetDryRowBound * wet_dry_bounds, const States * Statesptr)
{
	NUMERIC_TYPE loc_rainfall_total = C(0.0);

	const int row_start = wet_dry_bounds->dem_data[j].start;
	const int row_end = wet_dry_bounds->dem_data[j].end;

	// update wet_dry_bounds as all dem cells will now be wet
	wet_dry_bounds->fp_vol[j] = wet_dry_bounds->dem_data[j];

#ifdef __INTEL_COMPILER
				__assume_aligned(dem_row, 64);
				__assume_aligned(volume_row, 64);
#endif

#pragma ivdep
#pragma simd reduction (+:loc_rainfall_total)

	for (int i = row_start; i < row_end; i++)
	{
			//NUMERIC_TYPE dV = (dem_row[i] != C(1e10)) ? rain_step_dV : C(0.0);
			NUMERIC_TYPE dV;
			if (dem_row[i] != DEM_NO_DATA)
			{
				dV = rainmask_row[i] * delta_time;
				volume_row[i] += dV; // add rainfall volume to cell		
				loc_rainfall_total += dV; // mass balance for local cell (cumulative)
			}
		}
	return loc_rainfall_total;
}


//-----------------------------------------------------------------------------------
// BOUNDARY CONDITIONS
// Calculate Qx and Qy at edges of the domain in response to boundary
// conditions
void SGC2_BCs(const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE g,
	const NUMERIC_TYPE* dx_col, const NUMERIC_TYPE* dy_col,
	const NUMERIC_TYPE* h_grid,
	NUMERIC_TYPE* Qx_grid, NUMERIC_TYPE* Qy_grid, NUMERIC_TYPE* Qx_old_grid, NUMERIC_TYPE* Qy_old_grid,
	WetDryRowBound* wet_dry_bounds,
	const States *Statesptr, const Pars *Parptr, BoundaryCondition * boundary_cond, const SGCprams *SGCptr, const NUMERIC_TYPE max_Froude)
{
	//int i, j;
	//int BCi, index, sign, dir, edge, q_index, gr;
	//NUMERIC_TYPE h1, z1, hflow, dh, surface_slope, g, SGC_width_current, A, R, cell_length, cell_width;
	//NUMERIC_TYPE *q_FP_old, *q_SG_old, *q_FP_combined;
	NUMERIC_TYPE Q_multiplier = C(1.0);

	// CCS Multiplier for Q boundaries. If using regular grid, Qs are specified as m^2 and need to be multiplied by dx; 
	// if using lat-long Qs are specified in m^3 and therefore multiplier is C(1.) Note intialised above as C(1.0).
	if (Statesptr->latlong == OFF) Q_multiplier = Parptr->dx;

	const int numBCs = boundary_cond->bc_info.count;
	WaterSource ws = boundary_cond->bc_info;

	NUMERIC_TYPE * q_SG_old = ws.Q_SG_old;

	NUMERIC_TYPE q_in = C(0.0);
	NUMERIC_TYPE q_out = C(0.0);

	for (int BCi = 0; BCi < numBCs; BCi++)
	{
		if (ws.Ident[BCi] != NONE0) // if BCi = 0 do nothing
		{
			NUMERIC_TYPE *q_FP_old, *q_FP_combined;
			NUMERIC_TYPE cell_length, cell_width;
			int index, q_index, q_sg_old_index, q_fp_old_index;
			int i, j;
			int sign, gr;

			q_sg_old_index = BCi;
			//q_FP_old = ws.Q_FP_old + BCi;

			// First for each edge number work out where it is on the boundary,
			// the associated edge pixels and whether it's facing in the x or y
			// direction
			if (BCi < grid_cols)
			{
				// N(j=0) edge
				j = 0;
				i = BCi;

				index = i;
				q_index = i;// + j * grid_cols_padded;
				q_fp_old_index = q_index;

				q_FP_combined = Qy_grid;
				q_FP_old = Qy_old_grid;

				//SGC_qptr = SGC_Qy_grid + q_index;
				sign = -1;
				cell_length = dy_col[j];
				cell_width = dx_col[j];
				//printf("North %d, %d\n", i, j);
				//int j2 = (int)floor((double)index / grid_cols_padded);
				//int i2 = index - j*grid_cols_padded;
				//printf("Check %d, %d\n", i2, j2);
			}
			else if (/*BCi >= grid_cols && */ BCi < grid_cols + grid_rows)
			{
				// E edge
				j = BCi - grid_cols;
				i = grid_cols - 1;

				index = i + j * grid_cols_padded;
				q_index = index + 1;
				q_fp_old_index = q_index;

				q_FP_combined = Qx_grid;
				q_FP_old = Qx_old_grid;
				//SGC_qptr = SGC_Qx_grid + q_index;
				sign = 1;
				cell_length = dx_col[j];
				cell_width = dy_col[j];
				//printf("East %d, %d\n", i, j);
				//int j2 = (int)floor((double)index / grid_cols_padded);
				//int i2 = index - j*grid_cols_padded;
				//printf("Check %d, %d\n", i2, j2);
			}
			else if (/*BCi >= grid_cols + grid_rows &&*/ BCi < 2 * grid_cols + grid_rows)
			{
				// S(j=ysz-1) edge
				j = grid_rows - 1;
				i = (2 * grid_cols + grid_rows) - BCi - 1;

				index = i + j * grid_cols_padded;
				q_index = index + grid_cols_padded; // next row
				q_fp_old_index = q_index;

				q_FP_combined = Qy_grid;
				q_FP_old = Qy_old_grid;
				//SGC_qptr = SGC_Qy_grid + q_index;
				sign = 1;
				cell_length = dy_col[j];
				cell_width = dx_col[j];

				//printf("South %d, %d\n", i, j);
				//int j2 = (int)floor((double)index / grid_cols_padded);
				//int i2 = index - j*grid_cols_padded;
				//printf("Check %d, %d\n", i2, j2);
			}
			else
			{
				// W edge
				j = (2 * grid_cols + 2 * grid_rows) - BCi - 1;
				i = 0;

				index = i + j * grid_cols_padded;
				q_index = index;
				q_fp_old_index = q_index;

				q_FP_combined = Qx_grid;
				q_FP_old = Qx_old_grid;
				//SGC_qptr = SGC_Qx_grid + q_index;
				sign = -1;
				cell_length = dx_col[j];
				cell_width = dy_col[j];
				//printf("West  %d, %d\n", i, j);
				//int j2 = (int)floor((double)index / grid_cols_padded);
				//int i2 = index - j*grid_cols_padded;
				//printf("Check %d, %d\n", i2, j2);
			}

			//// CCS Record cell length and cell width relative to direction of channel (for lat-long grids)
			//if (edge == 1 || edge == 3) // N or S boundary; assume flow is N-S or S-N
			//{
			//	cell_length = dy_col[j];
			//	cell_width = dx_col[j];
			//}
			//if (edge == 2 || edge == 4) // E or W boundary; assume flow is E-W or W-E
			//{
			//	cell_length = dx_col[j];
			//	cell_width = dy_col[j];
			//}

			gr = ws.ws_cell.sg_cell_SGC_group[BCi];

			NUMERIC_TYPE g_friction_squared_FP = ws.g_friction_squared_FP[BCi];
			NUMERIC_TYPE g_friction_squared_SG = ws.g_friction_squared_SG[BCi];
			// Now calculate flows
			switch (ws.Ident[BCi])
			{
			case FREE1: // FREE boundary
			{
				if (h_grid[index] + ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi] > depth_thresh)
				{
					NUMERIC_TYPE qcorrected;
					// calcQ
					qcorrected = SGC2_CalcPointFREE(h_grid[index], ws.ws_cell.sg_cell_SGC_width[BCi], ws.Val[BCi],
						depth_thresh, delta_time, cell_width, g, g_friction_squared_SG, g_friction_squared_FP,
						ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi], gr, sign, &q_FP_old[q_fp_old_index], &q_SG_old[q_sg_old_index], SGCptr, max_Froude);

					q_FP_combined[q_index] = qcorrected + q_SG_old[q_sg_old_index];
				}
				else
				{
					q_FP_combined[q_index] = C(0.0);
					q_FP_old[q_fp_old_index] = C(0.0);
					q_SG_old[q_sg_old_index] = C(0.0);
				}
			}
			break;

			case HFIX2:// HFIX & HVAR boundary
			case HVAR3:// HFIX & HVAR boundary
			{
				if (h_grid[index] + ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi] > depth_thresh)
				{
					NUMERIC_TYPE surface_elevation0;
					if (ws.Ident[BCi] == HFIX2)
						surface_elevation0 = ws.Val[BCi];   // boundary depth for HFIX
					else // boundary depth for HVAR 
						surface_elevation0 = InterpolateTimeSeries(ws.timeSeries[BCi], curr_time);
					//note: h0 is absolute height - not relative height (equivilent to h1+dem[index])
					NUMERIC_TYPE h1, z1, SGC_width_current, hflow, dh, surface_slope, R, A;

					h1 = h_grid[index];     // cell depth
					z1 = ws.ws_cell.sg_cell_dem[BCi];   // FP elevation
					SGC_width_current = ws.ws_cell.sg_cell_SGC_width[BCi]; // SGC width

					if (SGC_width_current > C(0.0)) //  check for sub-grid channel
					{
						//surface_elevation0-z1 is depth above flood plain (may be negative) add bankfullheight to get depth above channel bed
						hflow = getmax(surface_elevation0 - z1, h1) + ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi]; // use max of cell depth and boundary depth
						//h0 is a surface elevation
						dh = surface_elevation0 - (h1 + z1);

						surface_slope = dh / (cell_length*SGCptr->SGCm[gr]);
						//if (edge == 1 || edge == 4) surface_slope = -surface_slope;
						surface_slope *= sign;

						SGC2_CalcA(gr, hflow, ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi], &A, &SGC_width_current, SGCptr); // calculate channel area for SGC
						R = SGC2_CalcR(gr, hflow, ws.ws_cell.sg_cell_SGC_BankFullHeight[BCi], SGC_width_current, ws.ws_cell.sg_cell_SGC_width[BCi], A, SGCptr); // calculate hydraulic radius for SGC

						q_SG_old[q_sg_old_index] = CalculateQ(surface_slope, R, delta_time, g, A, g_friction_squared_SG, q_SG_old[q_sg_old_index], max_Froude);
					}

					hflow = getmax(getmax(surface_elevation0 - z1, h1), C(0.0));
					// multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
					// FABS on surface_slope and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
					//if(hflow>depth_thresh && SGC_width_current < Parptr->dx) //CCS_deletion
					if (hflow > depth_thresh && SGC_width_current < cell_width)
					{
						dh = surface_elevation0 - (h1 + z1);

						surface_slope = dh / cell_length;

						surface_slope *= sign;

						A = cell_width * hflow;

						// calculate FP flow
						NUMERIC_TYPE q;
						q = CalculateQ(surface_slope, hflow, delta_time, g, A, g_friction_squared_FP, q_FP_old[q_fp_old_index], max_Froude);
						q_FP_old[q_fp_old_index] = q;

						if (SGC_width_current > C(0.0))
						{
							NUMERIC_TYPE channel_ratio = min(SGC_width_current / cell_width, C(1.0));
							q = q - channel_ratio * q;
						}
						q_FP_combined[q_index] = q;
					}
					else
					{
						q_FP_combined[q_index] = C(0.0);
						q_FP_old[q_fp_old_index] = C(0.0);
					}
					q_FP_combined[q_index] += q_SG_old[q_sg_old_index];
				}
				else
				{
					q_FP_combined[q_index] = C(0.0);
					q_FP_old[q_fp_old_index] = C(0.0);
					q_SG_old[q_sg_old_index] = C(0.0);
				}
			}
			break;

			case QFIX4:// QFIX boundary
			{
				//*qptr=-sign*ws.Val[BCi]*Parptr->dx; //CCS_deletion
				NUMERIC_TYPE q = -sign*ws.Val[BCi] * Q_multiplier;
				q_FP_combined[q_index] = q;
				q_FP_old[q_fp_old_index] = q;
				q_SG_old[q_sg_old_index] = C(0.0);
			}
			break;

			case QVAR5:// QVAR boundary
			{
				//*qptr=-sign*InterpolateTimeSeries(ws.TimeSeries[BCi],curr_time)*Parptr->dx; //CCS_deletion
				NUMERIC_TYPE q = -sign * InterpolateTimeSeries(ws.timeSeries[BCi], curr_time)*Q_multiplier;
				q_FP_combined[q_index] = q;
				q_FP_old[q_fp_old_index] = q;
				q_SG_old[q_sg_old_index] = C(0.0);
			}
			break;
			default:
				break;
			}
			NUMERIC_TYPE q = q_FP_combined[q_index];
			// ensure that any flow or change in volume is processed by subsequent steps
			if (q != C(0.0))
			{
				q *= sign;
				if (q > 0)
					q_out += q;
				else
					q_in -= q;

				wet_dry_bounds->fp_vol[j].start = min(wet_dry_bounds->fp_vol[j].start, i);
				wet_dry_bounds->fp_vol[j].end = max(wet_dry_bounds->fp_vol[j].end, i + 1);
			}
		}
	}
	boundary_cond->Qout = q_out;
	boundary_cond->Qin = q_in;

	return;
}

void SGC2_PointSources_Vol_row(const int y, const int grid_cols,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE g, const NUMERIC_TYPE Q_multiplier,
	const NUMERIC_TYPE *dx_col, const NUMERIC_TYPE *dy_col,
	const NUMERIC_TYPE * h_grid,
	const SGCprams *SGCptr,
	NUMERIC_TYPE * volume_row,
	PointSourceRowList * ps_layout,
	WetDryRowBound* wet_dry_bounds,
	NUMERIC_TYPE * out_Qpoint_timestep_pos, NUMERIC_TYPE * out_Qpoint_timestep_neg,
	const NUMERIC_TYPE max_Froude)
{
	NUMERIC_TYPE Qpoint_timestep_pos = C(0.0);
	NUMERIC_TYPE Qpoint_timestep_neg = C(0.0);
	const int ps_count = ps_layout->ps_row_count[y];
	const int row_cols_padded = ps_layout->row_cols_padded;
	const int row_start = y * row_cols_padded;
	WaterSource ps_info = ps_layout->ps_info;
	for (int i = 0; i < ps_count; i++)
	{
		int ws_index = row_start + i;

		NUMERIC_TYPE h;
		// Set initial dV and himp as zero
		NUMERIC_TYPE dV = C(0.0);
		int grid_index;
		int ps_x = ps_info.ws_cell.sg_cell_x[ws_index];
		int ps_y = ps_info.ws_cell.sg_cell_y[ws_index];
		// location in vector
		grid_index = ps_info.ws_cell.sg_cell_grid_index_lookup[ws_index];
		// different boundary conditions
		switch (ps_info.Ident[ws_index])
		{
			//NOTE HVAR and HFIX applied after update H
		case QVAR5: //QVAR ps.Val already set to the interpolated value
		case QFIX4:
			dV = ps_info.Val[ws_index] * Q_multiplier * delta_time; // QFIX // Calculate change in volume
			break;
		case FREE6:
			h = h_grid[grid_index] + ps_info.ws_cell.sg_cell_SGC_BankFullHeight[ws_index];
			if (h > depth_thresh)
			{
				NUMERIC_TYPE cell_width = getmin(dx_col[ps_y], dy_col[ps_y]);
				NUMERIC_TYPE FP_g_friction_squared = ps_info.g_friction_squared_FP[ws_index];
				NUMERIC_TYPE SGC_g_friction_squared = ps_info.g_friction_squared_SG[ws_index];

				NUMERIC_TYPE q_free_FP_corrected = SGC2_CalcPointFREE(h_grid[grid_index], ps_info.ws_cell.sg_cell_SGC_width[ws_index], ps_info.Val[ws_index],
					depth_thresh, delta_time, cell_width, g,
					SGC_g_friction_squared, FP_g_friction_squared,
					ps_info.ws_cell.sg_cell_SGC_BankFullHeight[ws_index], ps_info.ws_cell.sg_cell_SGC_group[ws_index],
					-1, &ps_info.Q_FP_old[ws_index], &ps_info.Q_SG_old[ws_index], SGCptr, max_Froude);

				NUMERIC_TYPE Qfree = q_free_FP_corrected + ps_info.Q_SG_old[ws_index];
				dV = Qfree * delta_time;
			}
			break;
		}

		if (dV != C(0.0))
		{
			wet_dry_bounds->fp_vol[ps_y].start = min(wet_dry_bounds->fp_vol[ps_y].start, ps_x);
			wet_dry_bounds->fp_vol[ps_y].end = max(wet_dry_bounds->fp_vol[ps_y].end, ps_x + 1);
			// update the cell volume change and in point source Q
			volume_row[ps_x] += dV; // Add volume to SGCdVol for use later by update H e.g wait for main update H before calculating H
			// Update Qpoint
			if (dV > 0)
				Qpoint_timestep_pos += dV;
			else
				Qpoint_timestep_neg += dV;
		}
	}
	(*out_Qpoint_timestep_pos) = Qpoint_timestep_pos;
	(*out_Qpoint_timestep_neg) = Qpoint_timestep_neg;
}

void SGC2_PointSources_H_row(const int y, const int grid_cols,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE *cell_area_col,
	const SGCprams *SGCptr,
	NUMERIC_TYPE * h_grid,
	NUMERIC_TYPE * volume_grid,
	PointSourceRowList * ps_layout,
	WetDryRowBound* wet_dry_bounds,
	NUMERIC_TYPE * out_Qpoint_timestep_pos, NUMERIC_TYPE * out_Qpoint_timestep_neg)
{
	// todo note - could be simpler to just calculate the new total volume based on height, rather than dv adjustment
	// less prone to errors and possibly faster

	NUMERIC_TYPE Qpoint_timestep_pos = C(0.0);
	NUMERIC_TYPE Qpoint_timestep_neg = C(0.0);

	const int ps_count = ps_layout->ps_row_count[y];
	const int row_cols_padded = ps_layout->row_cols_padded;
	const int row_start = y * row_cols_padded;
	WaterSource ps_info = ps_layout->ps_info;
	for (int i = 0; i < ps_count; i++)
	{
		int ws_index = row_start + i;
		// Set initial dV and himp as zero
		NUMERIC_TYPE V;
		NUMERIC_TYPE new_h = C(0.0);
		int gr;

		// location in grid

		const int ps_y = ps_info.ws_cell.sg_cell_y[ws_index];
		// location in vector
		const int grid_index = ps_info.ws_cell.sg_cell_grid_index_lookup[ws_index];
		// different boundary conditions
		if (ps_info.Ident[ws_index] == HFIX2 || ps_info.Ident[ws_index] == HVAR3) // HFIX or HVAR
		{
			// HVAR 'Val' already updated with the interpolated value for current time
			new_h = ps_info.Val[ws_index];
			new_h -= ps_info.ws_cell.sg_cell_dem[ws_index]; // get depth
			
			if (ps_info.ws_cell.sg_cell_SGC_width[ws_index] > C(0.0))
			{
				// sub-grid channel
				// ensure height is not below the bottom of the channel
				new_h = max(-ps_info.ws_cell.sg_cell_SGC_BankFullHeight[ws_index], new_h);

				// Calculate volume after update and subtract from before update
				gr = ps_info.ws_cell.sg_cell_SGC_group[ws_index]; // channel group number
				if (new_h < C(0.0) || ps_info.ws_cell.sg_cell_SGC_is_large[ws_index]) // if below the flood plain or cell is large Calculate channel volume
					V = SGC2_CalcUpV(new_h + ps_info.ws_cell.sg_cell_SGC_BankFullHeight[ws_index], ps_info.ws_cell.sg_cell_SGC_c[ws_index], gr, SGCptr);
				else // out of bank level
					V = ps_info.ws_cell.sg_cell_SGC_BankFullVolume[ws_index] + new_h * cell_area_col[ps_y];
			}
			else
			{
				// floodplain only cell
				new_h = max(C(0.0), new_h);
				V = new_h * cell_area_col[ps_y];
			}

			// ensure this point source is within the wet/dry bound (fp_vol_start,fp_vol_end used in update H to calculate the new fp_h_start)
			if (new_h > depth_thresh)
			{
				int ps_x = ps_info.ws_cell.sg_cell_x[ws_index];
				wet_dry_bounds->fp_vol[ps_y].start = min(wet_dry_bounds->fp_vol[ps_y].start, ps_x);
				wet_dry_bounds->fp_vol[ps_y].end = max(wet_dry_bounds->fp_vol[ps_y].end, ps_x + 1);
			}

			// calculate dV for mass balance
			NUMERIC_TYPE dV = V - volume_grid[grid_index];

			//h_grid[grid_index] = new_h;
			volume_grid[grid_index] = V;

			// Update Qpoint
			if (dV > 0)
				Qpoint_timestep_pos += dV;
			else
				Qpoint_timestep_neg += dV;
		}
	}
	(*out_Qpoint_timestep_pos) = Qpoint_timestep_pos;
	(*out_Qpoint_timestep_neg) = Qpoint_timestep_neg;
}

inline void SGC2_UpdateVol_floodplain_row(const int j, const int grid_row_index, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE row_cell_area,
	const NUMERIC_TYPE * Qx_grid, const NUMERIC_TYPE * Qy_grid,
	NUMERIC_TYPE * volume_grid, NUMERIC_TYPE * delta_volume_row, NUMERIC_TYPE * h_grid,

	const WetDryRowBound * wet_dry_bounds)
{
	bool not_last_row = (j < grid_rows - 1);
	// flow from next cell, previous row, next row
	int row_start = wet_dry_bounds->fp_vol[j].start - 1;
	int row_end = wet_dry_bounds->fp_vol[j].end + 1;
	if (j > 0)
	{
		row_start = min(row_start, wet_dry_bounds->fp_vol[j - 1].start);
		row_end = max(row_end, wet_dry_bounds->fp_vol[j - 1].end);
	}
	if (not_last_row)
	{
		row_start = min(row_start, wet_dry_bounds->fp_vol[j + 1].start);
		row_end = max(row_end, wet_dry_bounds->fp_vol[j + 1].end);
	}

	//ensure row_start, row_end are not out of bounds
	row_start = max(wet_dry_bounds->dem_data[j].start, row_start);
	row_end = min(row_end, wet_dry_bounds->dem_data[j].end);

	// update bounds for subsequent updates
	wet_dry_bounds->fp_vol[j].start = row_start;
	wet_dry_bounds->fp_vol[j].end = row_end;

	//__assume(row_start % GRID_ALIGN_WIDTH == 0);
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_cols_padded % GRID_ALIGN_WIDTH == 0);
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);

#endif
#ifdef __INTEL_COMPILER
	__assume_aligned(delta_volume_row, 64);
	__assume_aligned(Qx_grid, 64);
	__assume_aligned(Qy_grid, 64);
	__assume_aligned(volume_grid, 64);
#endif

#pragma ivdep
#pragma simd
	for (int i = row_start; i < row_end; i++)
	{
		int index = grid_row_index + i;
		int index_right = index + 1;
		int index_below = index + (grid_cols_padded);
		NUMERIC_TYPE dV = delta_time*(Qx_grid[index] - Qx_grid[index_right] + Qy_grid[index] - Qy_grid[index_below]); // compute volume change in cell
		dV += delta_volume_row[i];
		if (dV != C(0.0))
		{
			NUMERIC_TYPE V = volume_grid[index] + dV;
			volume_grid[index] = V;
			//h_grid[index] = V / row_cell_area;
		}
	}
	int count = row_end - row_start;
	if (count > 0)
		memset(delta_volume_row + row_start, 0, sizeof(NUMERIC_TYPE)*count);

#ifdef _DEBUG
	{
		int check_count = 0;
		for (int check = 0; check < grid_cols; check++)
		{
			if (delta_volume_row[check] != C(0.0))
			{
				check_count++;
				printf("dv!=0 %" NUM_FMT" (%d, %d),", delta_volume_row[check], check, j);
			}
		}
		if (check_count > 0)
			printf("\n");
	}
#endif
}

inline NUMERIC_TYPE SGC2_UpdateVol_sub_grid_row(const int sg_row_start, const int cell_count, const int grid_cols, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE evap_deltaH_step, const NUMERIC_TYPE delta_time,
	const int * sg_cell_grid_index_lookup,
	const int * sg_cell_x, const int * sg_cell_y,
	const NUMERIC_TYPE * sg_cell_SGC_BankFullHeight, const NUMERIC_TYPE * sg_cell_SGC_BankFullVolume,
	const NUMERIC_TYPE * sg_cell_cell_area,

	const SubGridFlowLookup * sg_cell_flow_lookup,
	const NUMERIC_TYPE * sg_flow_Q,
	//const NUMERIC_TYPE * sg_flow_ChannelRatio,
	const NUMERIC_TYPE * Qx_grid,
	const NUMERIC_TYPE * Qy_grid,
	const int * sg_cell_SGC_group,
	const NUMERIC_TYPE * sg_cell_SGC_c,
	const int * sg_cell_SGC_is_large, // constant_channel_width > 0.5*(cell_dx + cell_dy)

	NUMERIC_TYPE * volume_grid,
	NUMERIC_TYPE * h_grid,
	NetCDFVariable * evap_grid,
	WetDryRowBound * wet_dry_bounds,
	const SGCprams *SGCptr,
	const int SGCd8flag)
{
	NUMERIC_TYPE row_evap_loss = C(0.0);
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(sg_row_start % GRID_ALIGN_WIDTH == 0);
#endif
#ifdef __INTEL_COMPILER
	__assume_aligned(sg_cell_flow_lookup, 64);
	__assume_aligned(sg_cell_SGC_BankFullHeight, 64);
	__assume_aligned(sg_cell_grid_index_lookup, 64);
	__assume_aligned(sg_cell_SGC_BankFullVolume, 64);
	__assume_aligned(sg_cell_x, 64);
	__assume_aligned(sg_cell_y, 64);
	__assume_aligned(sg_cell_SGC_group, 64);
	__assume_aligned(sg_cell_SGC_c, 64);
	__assume_aligned(sg_cell_cell_area, 64);
	__assume_aligned(sg_cell_SGC_is_large, 64);
	__assume_aligned(sg_flow_Q, 64);

	__assume_aligned(h_grid, 64);
	__assume_aligned(volume_grid, 64);
	__assume_aligned(Qx_grid, 64);
	__assume_aligned(Qy_grid, 64);
#endif

	// Calc sub grid evaporation
	if (evap_deltaH_step > C(0.0))
	{
#pragma ivdep
		for (int cell_i = 0; cell_i < cell_count; cell_i++)
		{
			int cell_index = sg_row_start + cell_i;

			int grid_index = sg_cell_grid_index_lookup[cell_index];

			const NUMERIC_TYPE h_prev = h_grid[grid_index];

			const NUMERIC_TYPE cell_area = sg_cell_cell_area[cell_index];

			NUMERIC_TYPE evap_dV = C(0.0);

			NUMERIC_TYPE SGC_BankFullHeight = sg_cell_SGC_BankFullHeight[cell_index];
			NUMERIC_TYPE h_old = h_prev;
			NUMERIC_TYPE h_new = h_old - evap_grid->data[grid_index];
			// if there is water to evapourate: h_old + SGC_BankFullHeight) > depth_thresh
			// and the new height is below the flood plain h_new < 0 then need to update
			// note - if no water or water is still above the flood plain (h_new >= 0) then the flood plain calculation is correct

			if ((h_old + SGC_BankFullHeight) > depth_thresh)
			{
				NUMERIC_TYPE cell_area = sg_cell_cell_area[cell_index];
				if (h_new < C(0.0))
				{
					// this is a sub-grid channel negative height allowed
					// therefore undo the flood plain calculation (which truncated to zero)
					if (h_old > depth_thresh)
						//if (h_old > C(0.0))
					{
						evap_dV = h_old * cell_area;
					}

					h_old += SGC_BankFullHeight;
					h_new = h_old - evap_grid->data[grid_index];

					//ensure evapouration doesn't remove water that isn't there
					if (h_new < C(0.0))
					{
						h_new = C(0.0);
					}

					int gr = sg_cell_SGC_group[cell_index]; // channel group number
					NUMERIC_TYPE SGC_c = sg_cell_SGC_c[cell_index];
					// sub-grid channel evap or transition evap
					if (h_old < SGC_BankFullHeight || sg_cell_SGC_is_large[cell_index])
					{
						// calculate loss in vol
						evap_dV -= SGC2_CalcUpV(h_old, SGC_c, gr, SGCptr); //Calculate channel volume 
						evap_dV += SGC2_CalcUpV(h_new, SGC_c, gr, SGCptr); //Calculate channel volume 
					}
					// old water level must be above bank height and the channel is smaller than a cell width
					// but the new water level is below bank height, evap mass loss for bank transition
					else
					{
						// mass lost from floodplain
						NUMERIC_TYPE cell_evap = h_old - SGC_BankFullHeight;
						// mass lost from channel
						evap_dV -= SGC2_CalcUpV(SGC_BankFullHeight, SGC_c, gr, SGCptr); //Calculate bankfull area
						evap_dV += SGC2_CalcUpV(h_new, SGC_c, gr, SGCptr); //Calculate channel area
						evap_dV -= (cell_evap * cell_area);
					}
				}
				else if (h_old <= depth_thresh)
				{
					// normal flood plain evap - for the region between depth_thresh and zero.
					// this cell would have been skipped by the flood plain calculation, as the depth is below depth_thresh
					// evepouration needs to be calculated for this cell to allow for evapouration from the channel

					evap_dV -= (evap_grid->data[grid_index] * cell_area);

				} //else water is fully above the sub-grid, no need to update here
			}

			//dV = SGC2_Evaporation_SubGridCell(cell_index, grid_index, evap_deltaH_step, depth_thresh, sub_grid_cell_info, sub_grid_state, SGCptr);
			row_evap_loss -= evap_dV;
			//Cell_V[cell_index] = dV;
			//volume_grid[grid_index] += dV * delta_time;
			//dV += evap_dV;
			volume_grid[grid_index] += evap_dV;
		}
	}

	// PFU separate loop for subgrid flow_dV (needed if evap=0)
	#pragma ivdep
	for (int cell_i = 0; cell_i < cell_count; cell_i++)
	{
		int cell_index = sg_row_start + cell_i;

		int grid_index = sg_cell_grid_index_lookup[cell_index];

		SubGridFlowLookup sg_cell_flow_lookup_item = sg_cell_flow_lookup[cell_index];
		NUMERIC_TYPE flow_dV = C(0.0);

		// PFU use the subgrid flow add and sub variables to work out
		// how the subgrid flow contributes to the volume change
		// Replaces code in ProcessSubGridQBlock where channelQ was added to Q_grid

		// 0 qx (west side flows inward: add)
		int flow_index = sg_cell_flow_lookup_item.flow_add[0];
		//This adds the subgrid flow to the cell volume
		flow_dV += (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0); 
		
		// 0 qx (east side flows outward: subtract)
		flow_index = sg_cell_flow_lookup_item.flow_subtract[0];
		flow_dV -= (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);

		// 1 qy (north side flows inward: add)
		flow_index = sg_cell_flow_lookup_item.flow_add[1];
		flow_dV += (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);

		// 1 qy (south side flows outward: subtract)
		flow_index = sg_cell_flow_lookup_item.flow_subtract[1];
		flow_dV -= (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);

		// PFU: d8 flows
		if (SGCd8flag)
		{
			// 2 north west corner flows inward
			flow_index = sg_cell_flow_lookup_item.flow_add[2];
			flow_dV += (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);
			// 2 south east corner flows outward
			flow_index = sg_cell_flow_lookup_item.flow_subtract[2];
			flow_dV -= (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);

			// 3 north east corner flows inward
			flow_index = sg_cell_flow_lookup_item.flow_add[3];
			flow_dV += (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);
			// 3 south west corner flows outward
			flow_index = sg_cell_flow_lookup_item.flow_subtract[3];
			flow_dV -= (flow_index != -1) ? sg_flow_Q[flow_index] : C(0.0);
		}

		// PFU: removed channel ratio adjustment, now using adjusted floodplain width in Q_grid calculation

		//printf("Flow_dV: %" NUM_FMT,flow_dV);

		flow_dV *= delta_time;
		volume_grid[grid_index] += flow_dV;
	}
	return row_evap_loss;
}

inline NUMERIC_TYPE SGC2_ProcessH_Row(const int j, const int grid_cols, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE row_cell_area,
	const NUMERIC_TYPE sg_row_mem_size, const int sg_cell_row_count,
	NUMERIC_TYPE * h_grid, NUMERIC_TYPE * volume_grid,
	const int * sg_cell_x, const int * sg_cell_y,
	const int * sg_cell_grid_index_lookup,
	const NUMERIC_TYPE * sg_cell_SGC_BankFullHeight,
	const NUMERIC_TYPE * sg_cell_SGC_BankFullVolume,
	const NUMERIC_TYPE * sg_cell_cell_area,
	const int * sg_cell_SGC_group,
	const NUMERIC_TYPE * sg_cell_SGC_c,
	const int * sg_cell_SGC_is_large,
	WetDryRowBound * wet_dry_bounds,
	const SGCprams * SGCptr)
{
	int index, grid_row_index;
	NUMERIC_TYPE row_Hmax = C(0.0);

	int fp_h_start = grid_cols;
	int fp_h_end = -1;

	int row_start = wet_dry_bounds->fp_vol[j].start;
	int row_end = wet_dry_bounds->fp_vol[j].end;

	row_start = max(wet_dry_bounds->dem_data[j].start, row_start);
	row_end = min(row_end, wet_dry_bounds->dem_data[j].end);

	grid_row_index = j * grid_cols_padded;

#ifdef __INTEL_COMPILER
	__assume_aligned(h_grid, 64);
	__assume_aligned(volume_grid, 64);
	__assume_aligned(sg_cell_x, 64);
	__assume_aligned(sg_cell_y, 64);
	__assume_aligned(sg_cell_grid_index_lookup, 64);
	__assume_aligned(sg_cell_SGC_BankFullHeight, 64);
	__assume_aligned(sg_cell_SGC_BankFullVolume, 64);
	__assume_aligned(sg_cell_cell_area, 64);
	__assume_aligned(sg_cell_SGC_group, 64);
	__assume_aligned(sg_cell_SGC_c, 64);
	__assume_aligned(sg_cell_SGC_is_large, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(grid_cols_padded % GRID_ALIGN_WIDTH == 0);
	__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

#pragma ivdep
	for (int i = row_start; i < row_end; i++)
	{
		index = grid_row_index + i;

		NUMERIC_TYPE h;
		if (volume_grid[index] >= C(0.0))
		{
			h = volume_grid[index] / row_cell_area;
		}
		else
		{
#ifdef _DEBUG
			//printf("Volume <0 %" NUM_FMT" (%d,%d) \n", volume_grid[index], i, j);
#endif
			h = C(0.0);
			volume_grid[index] = C(0.0);
		}
		h_grid[index] = h;
	}

	const int sg_row_start = j * sg_row_mem_size;//  sub_grid_layout->row_cols_padded;
	const int cell_end = sg_cell_row_count;//  sub_grid_layout->cell_row_count[j];
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	__assume(sg_row_start % GRID_ALIGN_WIDTH == 0);
#endif
#pragma ivdep
	for (int cell_i = 0; cell_i < cell_end; cell_i++)
	{
		const int cell_index = sg_row_start + cell_i;
		//int x = sg_cell_x[cell_index];
		//int y = sg_cell_y[cell_index];
		const int grid_index = sg_cell_grid_index_lookup[cell_index];// x + y * grid_cols_padded;

		const NUMERIC_TYPE V = volume_grid[grid_index];
		NUMERIC_TYPE h;
		if (V >= sg_cell_SGC_BankFullVolume[cell_index] && !sg_cell_SGC_is_large[cell_index])// there is a sub-grid channel above it bank
		{
			// use V as the volume of water above the flood plain
			// simple h calculation possible
			//h = V / sg_cell_cell_area[cell_index];
			h = (V - sg_cell_SGC_BankFullVolume[cell_index]) / row_cell_area;
			if (h > depth_thresh)
			{
				// note this check and lock is only required when not processing row per thread
				NUMERIC_TYPE old_h = h_grid[grid_index]; // check if the most up to date h is below the depth_threshhold
				if (old_h < depth_thresh)
				{
					//int y = sg_cell_y[cell_index];
					int x = sg_cell_x[cell_index];

					// note this should be rare as only occurrs when a sub grid flow first overflows onto the floodplain
					{
						// flow onto flood plain update flood plain wet/dry
						row_start = min(row_start, x);
						row_end = max(row_end, x + 1);
						//wet_dry_bounds->fp_vol[y].start = min(wet_dry_bounds->fp_vol[y].start, x);
						//wet_dry_bounds->fp_vol[y].end = max(wet_dry_bounds->fp_vol[y].end, x + 1);
					}
				}

			}
		}
		else if (V > C(0.0)) // there is a sub-grid channel and its within bank
		{
			int channel_group = sg_cell_SGC_group[cell_index];
			h = SGC2_CalcUpH(V, sg_cell_SGC_c[cell_index], channel_group, SGCptr);
			h -= sg_cell_SGC_BankFullHeight[cell_index];
		}
		else // V < 0
		{
			//printf("check negative volume (sub grid) %" NUM_FMT"\n", V);
			h = -sg_cell_SGC_BankFullHeight[cell_index];
		}

		h_grid[grid_index] = h;

		h += sg_cell_SGC_BankFullHeight[cell_index];

		if (h > row_Hmax)
			row_Hmax = h;
	}

	row_start = wet_dry_bounds->fp_vol[j].start;
	row_end = wet_dry_bounds->fp_vol[j].end;

	row_start = max(wet_dry_bounds->dem_data[j].start, row_start);
	row_end = min(row_end, wet_dry_bounds->dem_data[j].end);


#pragma ivdep
	for (int i = row_start; i < row_end; i++)
	{
		index = grid_row_index + i;

		NUMERIC_TYPE h = h_grid[index];
		if (h > row_Hmax)
			row_Hmax = h;
		// update wet/dry boundary
		//if (h_grid[index] > C(0.0) /*depth_thresh*/)
		if (h > depth_thresh)
		{
			if (i < fp_h_start)
				fp_h_start = i;
			if (i > fp_h_end)
				fp_h_end = i;
		}
	}

	wet_dry_bounds->fp_h_prev[j].start = wet_dry_bounds->fp_h[j].start;
	wet_dry_bounds->fp_h_prev[j].end = wet_dry_bounds->fp_h[j].end;
	wet_dry_bounds->fp_h[j].start = fp_h_start;
	if (fp_h_end == -1)
	{
		wet_dry_bounds->fp_h[j].end = fp_h_end;
	}
	else
	{
		wet_dry_bounds->fp_h[j].end = min(fp_h_end + 1, wet_dry_bounds->dem_data[j].end);
	}
	return row_Hmax;
}

void SGC2_UpdateHazard_row(const int j, const int grid_row_index,
	const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE * Vx_grid, const NUMERIC_TYPE * Vy_grid,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * dem_grid,
	const NUMERIC_TYPE row_dx, const NUMERIC_TYPE row_dy,
	const NUMERIC_TYPE *SGC_BankFullHeight_grid,
	NUMERIC_TYPE *maxVc_grid, NUMERIC_TYPE *maxVc_height_grid, NUMERIC_TYPE *maxHazard_grid,
	WetDryRowBound *wet_dry_bounds)
{
	bool not_last_row = (j < grid_rows - 1);
	// flow from next cell, previous row, next row
	int row_start = wet_dry_bounds->fp_vol[j].start - 1;
	int row_end = wet_dry_bounds->fp_vol[j].end + 1;
	if (j > 0)
	{
		row_start = min(row_start, wet_dry_bounds->fp_vol[j - 1].start);
		row_end = max(row_end, wet_dry_bounds->fp_vol[j - 1].end);
	}
	if (not_last_row)
	{
		row_start = min(row_start, wet_dry_bounds->fp_vol[j + 1].start);
		row_end = max(row_end, wet_dry_bounds->fp_vol[j + 1].end);
	}

	//ensure row_start, row_end are not out of bounds
	row_start = max(wet_dry_bounds->dem_data[j].start, row_start);
	row_end = min(row_end, wet_dry_bounds->dem_data[j].end);

	for (int i = row_start; i < row_end; i++)
	{
		int index = grid_row_index + i;

		// calculate velocities, ignore boundary velocities
		//NUMERIC_TYPE Vx_west = (i != 0) ? SGC2_CalculateVelocity(index - 1, index, Qx_grid, h_grid, dem_grid, row_dy) : C(0.0);
		//NUMERIC_TYPE Vx_east = (i < grid_cols - 1) ? SGC2_CalculateVelocity(index, index + 1, Qx_grid, h_grid, dem_grid, row_dy) : C(0.0);

		//NUMERIC_TYPE Vy_north = (j != 0) ? SGC2_CalculateVelocity(index - grid_cols_padded, index, Qy_grid, h_grid, dem_grid, row_dx) : C(0.0);
		//NUMERIC_TYPE Vy_south = (not_last_row) ? SGC2_CalculateVelocity(index, index + grid_cols_padded, Qy_grid, h_grid, dem_grid, row_dx) : C(0.0);
		NUMERIC_TYPE Vx_west = Vx_grid[index];
		NUMERIC_TYPE Vx_east = Vx_grid[index + 1];

		NUMERIC_TYPE Vy_north = Vy_grid[index];
		NUMERIC_TYPE Vy_south = Vy_grid[index + grid_cols_padded];

		NUMERIC_TYPE Vx = getmax(FABS(Vx_west), FABS(Vx_east));
		NUMERIC_TYPE Vy = getmax(FABS(Vy_north), FABS(Vy_south));

		NUMERIC_TYPE Vc = SQRT(Vx*Vx + Vy*Vy);
		NUMERIC_TYPE depth = h_grid[index] + SGC_BankFullHeight_grid[index];
		NUMERIC_TYPE hazard = depth * (Vc + C(1.5)); // Changed to equation from DEFRA 2006 (ALD)
		maxHazard_grid[index] = getmax(hazard, maxHazard_grid[index]);

		if (Vc > maxVc_grid[index])
		{
			maxVc_grid[index] = Vc;
			maxVc_height_grid[index] = depth;
		}
	}
}

NUMERIC_TYPE SGC2_UpdateVolumeHeight_block(const int block_index,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE curr_time,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE g,
	const NUMERIC_TYPE evap_deltaH_step, NetCDFVariable * evap_grid, const NUMERIC_TYPE rain_deltaH_step, const NUMERIC_TYPE * rain_grid, const NUMERIC_TYPE * dist_infil_grid,

	WetDryRowBound *wet_dry_bounds, NUMERIC_TYPE *tmp_row,
	const NUMERIC_TYPE* Qx_grid, const NUMERIC_TYPE* Qy_grid,
	const NUMERIC_TYPE * cell_area_col, const NUMERIC_TYPE * dx_col, const NUMERIC_TYPE * dy_col,
	const NUMERIC_TYPE * dem_grid, 

	const SubGridRowList * sub_grid_layout, const SubGridState * sub_grid_state,
	PointSourceRowList * ps_layout,
	const RouteDynamicList * route_dynamic_list,

	NUMERIC_TYPE * h_grid,
	NUMERIC_TYPE * volume_grid,

	const NUMERIC_TYPE * SGC_BankFullHeight_grid,
	const States *Statesptr,
	Pars *Parptr,
	const Solver *Solverptr,
	const SGCprams *SGCptr,
	//DynamicRain<> & dynamic_rain, removed JCN
	VolumeHeightUpdateInfo * update_info)
{
	NUMERIC_TYPE block_evap_loss = C(0.0);
	NUMERIC_TYPE block_rain_total = C(0.0);
	NUMERIC_TYPE block_Qpoint_timestep_pos = C(0.0);
	NUMERIC_TYPE block_Qpoint_timestep_neg = C(0.0);

	NUMERIC_TYPE block_Hmax = C(0.0);

	NUMERIC_TYPE Q_multiplier = C(1.0); // for volume point sources
	// CCS Multiplier for Q inputs. If using regular grid, Qs are specified as m^2 and need to be multiplied by dx; 
	// if using lat-long Qs are specified in m^3 and therefore multiplier is C(1.) Note intialised above as C(1.0).
	if (Statesptr->latlong == OFF) Q_multiplier = Parptr->dx;

	SubGridFlowLookup * sg_cell_flow_lookup = sub_grid_layout->flow_info.sg_cell_flow_lookup;
	const NUMERIC_TYPE * sg_cell_SGC_BankFullHeight = sub_grid_layout->cell_info.sg_cell_SGC_BankFullHeight;
	const int * sg_cell_grid_index_lookup = sub_grid_layout->cell_info.sg_cell_grid_index_lookup;
	const NUMERIC_TYPE * sg_cell_SGC_BankFullVolume = sub_grid_layout->cell_info.sg_cell_SGC_BankFullVolume;
	const int * sg_cell_x = sub_grid_layout->cell_info.sg_cell_x;
	const int * sg_cell_y = sub_grid_layout->cell_info.sg_cell_y;
	const int * sg_cell_SGC_group = sub_grid_layout->cell_info.sg_cell_SGC_group;
	const NUMERIC_TYPE * sg_cell_SGC_c = sub_grid_layout->cell_info.sg_cell_SGC_c;
	const NUMERIC_TYPE * sg_cell_cell_area = sub_grid_layout->cell_info.sg_cell_cell_area;
	const int * sg_cell_SGC_is_large = sub_grid_layout->cell_info.sg_cell_SGC_is_large;

//	const NUMERIC_TYPE * sg_flow_ChannelRatio = sub_grid_state->sg_flow_ChannelRatio;
	const NUMERIC_TYPE * sg_flow_Q = sub_grid_state->sg_flow_Q;

	NUMERIC_TYPE q_pos, q_neg;

	const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
	const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

	for (int j = start_y; j < end_y; j++)
	{
		const int grid_row_index = j * grid_cols_padded;
		const NUMERIC_TYPE row_cell_area = cell_area_col[j];
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
		__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

		NUMERIC_TYPE * delta_volume_row = tmp_row;

#ifdef __INTEL_COMPILER
		__assume_aligned(delta_volume_row, 64);

		__assume_aligned(h_grid, 64);
		__assume_aligned(dist_infil_grid, 64);
		__assume_aligned(volume_grid, 64);
		__assume_aligned(dem_grid, 64);
		__assume_aligned(cell_area_col, 64);
		__assume_aligned(Qx_grid, 64);
		__assume_aligned(Qy_grid, 64);

		__assume_aligned(sg_cell_flow_lookup, 64);
		__assume_aligned(sg_cell_SGC_BankFullHeight, 64);
		__assume_aligned(sg_cell_grid_index_lookup, 64);
		__assume_aligned(sg_cell_SGC_BankFullVolume, 64);
		__assume_aligned(sg_cell_x, 64);
		__assume_aligned(sg_cell_y, 64);
		__assume_aligned(sg_cell_SGC_group, 64);
		__assume_aligned(sg_cell_SGC_c, 64);
		__assume_aligned(sg_cell_cell_area, 64);
		__assume_aligned(sg_cell_SGC_is_large, 64);
//		__assume_aligned(sg_flow_ChannelRatio, 64);
		__assume_aligned(sg_flow_Q, 64);

#endif		

		SGC2_PointSources_Vol_row(j, grid_cols, delta_time, curr_time, depth_thresh, g, Q_multiplier, dx_col, dy_col, h_grid, SGCptr, delta_volume_row, ps_layout, wet_dry_bounds, &q_pos, &q_neg, Parptr->max_Froude);
		block_Qpoint_timestep_pos += q_pos;
		block_Qpoint_timestep_neg += q_neg;

		
		if ((Statesptr->dynamicrainfall == ON || Statesptr->rainfallmask == ON) && wet_dry_bounds->dem_data[j].start > -1)
		{
			block_evap_loss -= SGC2_Distrubuted_Rainfall_row(j, delta_time,
				dem_grid + grid_row_index, rain_grid + grid_row_index, // pointer to start of row
				delta_volume_row, wet_dry_bounds, Statesptr);
		}
		else if (rain_deltaH_step > C(0.0) && wet_dry_bounds->dem_data[j].start > -1)
		{
			const NUMERIC_TYPE rain_step_dV = rain_deltaH_step * row_cell_area; // area constant per row => rainfall constant per row

			block_evap_loss -= SGC2_Uniform_Rainfall_row(j, rain_step_dV,
				dem_grid + grid_row_index, // pointer to start of row
				delta_volume_row, wet_dry_bounds, Statesptr);
		}
		
		// const_cell_evap always zero if evap Evaporation
		if (evap_deltaH_step > C(0.0))
		{
			int evap_row_start = wet_dry_bounds->fp_h[j].start;
			int evap_row_end = min(wet_dry_bounds->fp_h[j].end, wet_dry_bounds->dem_data[j].end);
			block_evap_loss += SGC2_Evaporation_floodplain_row(evap_row_start, evap_row_end, depth_thresh, row_cell_area, evap_deltaH_step,
									   evap_grid->data + grid_row_index,
									   h_grid + grid_row_index, // pointer to start of row
									   delta_volume_row);
		}
		if (Statesptr->calc_distributed_infiltration == ON)
		{
			int evap_row_start = wet_dry_bounds->fp_h[j].start;
			int evap_row_end = min(wet_dry_bounds->fp_h[j].end, wet_dry_bounds->dem_data[j].end);
			block_evap_loss += SGC2_Infil_floodplain_row(evap_row_start, evap_row_end, depth_thresh, row_cell_area, evap_deltaH_step,
				h_grid + grid_row_index, // pointer to start of row
				dist_infil_grid + grid_row_index, // pointer to start of distributed infiltratioon grid
				delta_volume_row);
		}

		SGC2_UpdateVol_floodplain_row(j, grid_row_index, grid_cols, grid_rows, grid_cols_padded, delta_time, row_cell_area, Qx_grid, Qy_grid, volume_grid, delta_volume_row, h_grid, wet_dry_bounds);


#if _DEBUG
		for (int i = 0; i < grid_cols; i++)
		{
			if (delta_volume_row[i] != C(0.0))
				printf("delta_volume_row[%d] not 0\n", i);
		}
#endif

		const int sg_row_start = j * sub_grid_layout->row_cols_padded;
		const int cell_row_count = sub_grid_layout->cell_row_count[j];

		// note could use row_cell_area, currently looks up for each cell (SGC2_UpdateVol_sub_grid_row function can be used for single row or all rows with altered memory structure)
		block_evap_loss += SGC2_UpdateVol_sub_grid_row(sg_row_start, cell_row_count, grid_cols, grid_cols_padded,
			depth_thresh, evap_deltaH_step, delta_time, sg_cell_grid_index_lookup, sg_cell_x, sg_cell_y, sg_cell_SGC_BankFullHeight, sg_cell_SGC_BankFullVolume,
			sg_cell_cell_area, sg_cell_flow_lookup, sg_flow_Q,  Qx_grid, Qy_grid, sg_cell_SGC_group, sg_cell_SGC_c, sg_cell_SGC_is_large,
			volume_grid, h_grid, evap_grid, wet_dry_bounds, SGCptr, Statesptr->SGCd8);

		SGC2_PointSources_H_row(j, grid_cols, delta_time, curr_time, depth_thresh, cell_area_col, SGCptr, h_grid, volume_grid, ps_layout, wet_dry_bounds, &q_pos, &q_neg);
		block_Qpoint_timestep_pos += q_pos;
		block_Qpoint_timestep_neg += q_neg;

		NUMERIC_TYPE row_Hmax = SGC2_ProcessH_Row(j, grid_cols, grid_cols_padded, depth_thresh, row_cell_area,
			sub_grid_layout->row_cols_padded, sub_grid_layout->cell_row_count[j],
			h_grid, volume_grid, sg_cell_x, sg_cell_y, sg_cell_grid_index_lookup,
			sg_cell_SGC_BankFullHeight, sg_cell_SGC_BankFullVolume, sg_cell_cell_area, sg_cell_SGC_group, sg_cell_SGC_c, sg_cell_SGC_is_large, wet_dry_bounds, SGCptr);

		if (row_Hmax > block_Hmax)
			block_Hmax = row_Hmax;
	}

	update_info->evap_loss = block_evap_loss;
	update_info->rain_total = block_rain_total;
	update_info->Qpoint_timestep_pos = block_Qpoint_timestep_pos;
	update_info->Qpoint_timestep_neg = block_Qpoint_timestep_neg;

	return block_Hmax;
}

void SGC2_Inundation_block(const int block_index, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh,
	const NUMERIC_TYPE current_time_hours,
	const NUMERIC_TYPE delta_time_hours,
	const NUMERIC_TYPE * h_grid,
	const SubGridRowList * sub_grid_layout, const SubGridState * sub_grid_state,
	WetDryRowBound* wet_dry_bounds,
	NUMERIC_TYPE * initHtm_grid,
	NUMERIC_TYPE * totalHtm_grid,
	NUMERIC_TYPE * maxH_grid,
	NUMERIC_TYPE * maxHtm_grid)
{
	const int * sg_cell_grid_index_lookup = sub_grid_layout->cell_info.sg_cell_grid_index_lookup;
	const NUMERIC_TYPE * sg_cell_SGC_BankFullHeight = sub_grid_layout->cell_info.sg_cell_SGC_BankFullHeight;
	{
		const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
		const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

		for (int j = start_y; j < end_y; j++)
		{
			int i, index, grid_row_index;
			grid_row_index = j * grid_cols_padded;

			int row_start = wet_dry_bounds->fp_h[j].start;
			int row_end = wet_dry_bounds->fp_h[j].end;

#ifdef __INTEL_COMPILER
			__assume_aligned(h_grid, 64);
			__assume_aligned(initHtm_grid, 64);
			__assume_aligned(totalHtm_grid, 64);
			__assume_aligned(maxH_grid, 64);
			__assume_aligned(maxHtm_grid, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
			__assume(grid_cols_padded % GRID_ALIGN_WIDTH == 0);
			__assume(grid_row_index % GRID_ALIGN_WIDTH == 0);
#endif

#pragma ivdep
#pragma simd
			for (int i = row_start; i < row_end; i++)
			{
				index = grid_row_index + i;
				if (h_grid[index] > depth_thresh)
				{
					if (initHtm_grid[index] == NULLVAL)
						initHtm_grid[index] = current_time_hours;

					totalHtm_grid[index] += delta_time_hours;
					// Update maximum water depths, and time of maximum (in hours)
					maxHtm_grid[index] = (h_grid[index] > maxH_grid[index]) ? current_time_hours : maxHtm_grid[index];
					maxH_grid[index] = (h_grid[index] > maxH_grid[index]) ? h_grid[index] : maxH_grid[index];
				}
			}

#ifdef __INTEL_COMPILER
			__assume_aligned(h_grid, 64);
			__assume_aligned(initHtm_grid, 64);
			__assume_aligned(totalHtm_grid, 64);
			__assume_aligned(maxH_grid, 64);
			__assume_aligned(maxHtm_grid, 64);
			__assume_aligned(sg_cell_grid_index_lookup, 64);
			__assume_aligned(sg_cell_SGC_BankFullHeight, 64);
#endif

			const int sg_row_start = j * sub_grid_layout->row_cols_padded;
			const int cell_end = sub_grid_layout->cell_row_count[j];
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
			__assume(sg_row_start % GRID_ALIGN_WIDTH == 0);
#endif
#pragma ivdep
			for (int cell_i = 0; cell_i < cell_end; cell_i++)
			{
				const int cell_index = sg_row_start + cell_i;
				const int grid_index = sg_cell_grid_index_lookup[cell_index];
				const NUMERIC_TYPE h_tmp = h_grid[grid_index] + sg_cell_SGC_BankFullHeight[cell_index];
				if (h_tmp > depth_thresh) // use depth_thresh as threshhold
					//if (depth > C(0.01)) // used to be hardcoded at 0.01 - maintain consistent output
				{
					if (initHtm_grid[grid_index] == NULLVAL)
						initHtm_grid[grid_index] = current_time_hours;
					//only add to time if it hasn't already been added by the flood plain calculation
					if (h_grid[grid_index] <= depth_thresh)
						totalHtm_grid[grid_index] += delta_time_hours;
					// Update maximum water depths, and time of maximum (in hours)
					if (h_tmp > maxH_grid[grid_index])
					{
						maxH_grid[grid_index] = h_tmp;
						maxHtm_grid[grid_index] = current_time_hours;
					}
				}
			}
		}
	}
}

NUMERIC_TYPE SGC2_InitHBounds(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh,
	const SubGridRowList * sub_grid_layout,
	SubGridState * sub_grid_state,
	const NUMERIC_TYPE * cell_area_col,
	NUMERIC_TYPE * h_grid, NUMERIC_TYPE * volume_grid,
	WetDryRowBound* wet_dry_bounds, const SGCprams * SGCptr)
{
	NUMERIC_TYPE Hmax = C(0.0);

	const int * sg_cell_x = sub_grid_layout->cell_info.sg_cell_x;
	const int * sg_cell_y = sub_grid_layout->cell_info.sg_cell_y;
	const int * sg_cell_grid_index_lookup = sub_grid_layout->cell_info.sg_cell_grid_index_lookup;
	const NUMERIC_TYPE * sg_cell_SGC_BankFullHeight = sub_grid_layout->cell_info.sg_cell_SGC_BankFullHeight;
	const NUMERIC_TYPE * sg_cell_SGC_BankFullVolume = sub_grid_layout->cell_info.sg_cell_SGC_BankFullVolume;
	const NUMERIC_TYPE * sg_cell_cell_area = sub_grid_layout->cell_info.sg_cell_cell_area;

	const NUMERIC_TYPE * sg_cell_SGC_c = sub_grid_layout->cell_info.sg_cell_SGC_c;
	const int * sg_cell_SGC_group = sub_grid_layout->cell_info.sg_cell_SGC_group;
	const int * sg_cell_SGC_is_large = sub_grid_layout->cell_info.sg_cell_SGC_is_large;

	// reduction to find the maximum h (note MS compiler OpenMP version does not support max reduction, use critical to prevent race condition)
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma omp parallel for default(shared) schedule(static)
#else
#pragma omp parallel for default(shared) reduction(max : Hmax) schedule(static)
#endif
	//for (int j = 0; j < grid_rows; j++)
	for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
	{
		const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
		const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

		for (int j = start_y; j < end_y; j++)
		{
			NUMERIC_TYPE row_Hmax = SGC2_ProcessH_Row(j, grid_cols, grid_cols_padded, depth_thresh,
				cell_area_col[j],
				sub_grid_layout->row_cols_padded, sub_grid_layout->cell_row_count[j],
				h_grid, volume_grid, sg_cell_x, sg_cell_y, sg_cell_grid_index_lookup, sg_cell_SGC_BankFullHeight, sg_cell_SGC_BankFullVolume, sg_cell_cell_area, sg_cell_SGC_group, sg_cell_SGC_c, sg_cell_SGC_is_large, wet_dry_bounds, SGCptr);


#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma omp critical
#endif
			{
				if (row_Hmax > Hmax)
					Hmax = row_Hmax;
			}
		}
	}

	for (int j = 0; j < grid_rows; j++)
	{
		wet_dry_bounds->fp_vol[j].start = wet_dry_bounds->fp_h[j].start;
		wet_dry_bounds->fp_vol[j].end = wet_dry_bounds->fp_h[j].end;
	}

	return Hmax;
}


void SGC2_UpdateLoadBalance(const int grid_rows, const int grid_cols_padded,
	const SubGridRowList * sub_grid_layout,
	WetDryRowBound* wet_dry_bounds)
{
	//printf("SGC2_UpdateLoadBalance\n");
	const int block_count = wet_dry_bounds->block_count;

	IndexRange* block_row_bounds = wet_dry_bounds->block_row_bounds;

#ifdef __INTEL_COMPILER
	__assume_aligned(block_row_bounds, 64);
#endif

	//SetArrayValue(block_start_row, -1, THREAD_ROW_BLOCK);
	//SetArrayValue(block_end_row, -1, THREAD_ROW_BLOCK);
	if (grid_rows > block_count)
	{
		int total_work = 0;
		for (int j = 0; j < grid_rows; j++)
		{
			//int grid_row_index = j * grid_cols_padded;
#if _BALANCE_TYPE == 1
			int row_start = wet_dry_bounds->fp_h[j].start;
			int row_end = wet_dry_bounds->fp_h[j].end;
			int work_this_row = ((row_end == -1) ? 0 : (row_end - row_start)) + 10;
#else
			int work_this_row = 1;
#endif

			total_work += work_this_row;
		}


		float target_work_per_block = (float)total_work / block_count;
		int allocated_work = 0;
		int block_index = 0;

		block_row_bounds[0].start = 0;
#ifdef _DEBUG
		//printf("total work: %d block target: %d\n", total_work, target_work_per_block);
#endif

		for (int j = 0; j < grid_rows; j++)
		{
			if (block_index < (block_count - 1) &&
				(allocated_work) >= (block_index + 1) * target_work_per_block)
			{
				block_row_bounds[block_index].end = j;
				block_index++;
				block_row_bounds[block_index].start = j;
			}


			int row_start = wet_dry_bounds->fp_h[j].start;
			int row_end = wet_dry_bounds->fp_h[j].end;
#if _BALANCE_TYPE == 1
			int work_this_row = ((row_end == -1) ? 0 : (row_end - row_start)) + 10;
#else
			int work_this_row = 1;
#endif

			allocated_work += work_this_row;

		}
		block_row_bounds[block_index].end = grid_rows;
#ifdef _DEBUG
		for (int bi = 0; bi < block_count; bi++)
		{
			int rows_in_block = block_row_bounds[bi].end - block_row_bounds[bi].start;
			NUMERIC_TYPE percent_rows_in_block = ((NUMERIC_TYPE)rows_in_block / grid_rows) * 100;

			int work_this_block1 = 0;
			for (int j = block_row_bounds[bi].start; j < block_row_bounds[bi].end; j++)
			{
				int row_start = wet_dry_bounds->fp_h[j].start;
				int row_end = wet_dry_bounds->fp_h[j].end;
#if _BALANCE_TYPE == 1
				int work_this_row = ((row_end == -1) ? 0 : (row_end - row_start)) + 10;
#else
				int work_this_row = 1;
#endif
				work_this_block1 += work_this_row;
			}

			NUMERIC_TYPE percent_work1 = ((NUMERIC_TYPE)work_this_block1 / total_work) * 100;
			//NUMERIC_TYPE percent_work2 = ((NUMERIC_TYPE)work_this_block2 / total_work) * 100;
			printf("block: %d row: %d -> %d [%d](%.2" NUM_FMT"%%) work: [%d : %.2" NUM_FMT"%%]\n",
				bi,
				block_row_bounds[bi].start, block_row_bounds[bi].end,
				rows_in_block, percent_rows_in_block, work_this_block1, percent_work1);

		}
#endif
	}
	else
	{
		//int block_index = 0;
		for (int j = 0; j < grid_rows; j++)
		{
			block_row_bounds[j].start = j;
			block_row_bounds[j].end = j + 1;
		}
	}
#ifdef _DEBUG
	//printf("FloodPlain, %d, %d %lf\n", wetCount, wetBound, (double)wetCount / wetBound);
#endif

}

//---------------------------------------------------------------------------
// CALCULATE VOLUME OF WATER IN CHANNEL AND FLOODPLAIN
// CALCULATE FLOOD AREA
void SGC2_DomainVolumeAndFloodArea_block(const int block_index, const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE depth_thresh,
	const WetDryRowBound* wet_dry_bounds,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * cell_area_col,
	const NUMERIC_TYPE * volume_grid,
	NUMERIC_TYPE *out_flood_area, NUMERIC_TYPE *out_domain_volume)
{
	// Calculate flood area
	const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
	const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

	NUMERIC_TYPE block_flood_area = C(0.0);
	NUMERIC_TYPE block_domain_volume = C(0.0);

	for (int j = start_y; j < end_y; j++)
	{
#ifdef __INTEL_COMPILER
		__assume_aligned(cell_area_col, 64);
		__assume_aligned(h_grid, 64);
		__assume_aligned(volume_grid, 64);
#endif
		NUMERIC_TYPE dA = cell_area_col[j]; // if latlong is on change dA to local cell area

#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
		__assume(grid_cols_padded % GRID_ALIGN_WIDTH == 0);
#endif

#pragma simd reduction(+: block_flood_area)
		for (int i = 0; i < grid_cols; i++)
		{
			int index = i + j*grid_cols_padded;
			// only count height above the flood plain ( If sub-grid used channel depth is negative)
			block_flood_area += (h_grid[index] > depth_thresh) ? dA : C(0.0);
		}
#pragma simd reduction(+: block_domain_volume)
		for (int i = 0; i < grid_cols; i++)
		{
			int index = i + j*grid_cols_padded;
			block_domain_volume += volume_grid[index];
		}
	}
	(*out_flood_area) = block_flood_area;
	(*out_domain_volume) = block_domain_volume;
}


#ifdef RESULT_CHECK
NUMERIC_TYPE Do_Update_old(States *Statesptr, Pars *Parptr, Solver *Solverptr, Arrays * Arrptr, SGCprams * SGCptr, BoundCs *BCptr, ChannelSegmentType *ChannelSegments, vector<ChannelSegmentType> *ChannelSegmentsVecPtr)
{
	// sub grid floodplain models
	SGC_FloodplainQ(Statesptr, Parptr, Solverptr, Arrptr, SGCptr);

	SGC_BCs(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);
	// Infiltration, evaporation and rainfall routines after time step update (TJF)11
	if (Statesptr->calc_evap == ON) SGC_Evaporation(Parptr, Solverptr, Arrptr, SGCptr);
	if (Statesptr->rainfall == ON && Statesptr->routing == OFF) SGC_Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
	if (Statesptr->routing == ON) SGC_Routing(Statesptr, Parptr, Solverptr, Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
	if (Statesptr->hazard == ON) UpdateV(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);
	NUMERIC_TYPE Hmax = SGC_UpdateH(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);

	BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);
	// NOTES: Update Q's handeled within SGC flux equations (SGC_FloodplainQ etc.) time-step calculation intergrated with UpdateH

	return Hmax;
}
#endif

void Do_Update(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE curr_time, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE g,
	NUMERIC_TYPE *h_grid, NUMERIC_TYPE *volume_grid,
	NUMERIC_TYPE *Qx_grid, NUMERIC_TYPE *Qy_grid, NUMERIC_TYPE *Qx_old_grid, NUMERIC_TYPE *Qy_old_grid,
	NUMERIC_TYPE *initHtm_grid, NUMERIC_TYPE *maxHtm_grid, NUMERIC_TYPE *totalHtm_grid, NUMERIC_TYPE *maxH_grid,
	NUMERIC_TYPE * maxVc_grid, NUMERIC_TYPE * maxVc_height_grid, NUMERIC_TYPE * maxHazard_grid,
	NUMERIC_TYPE * Vx_grid, NUMERIC_TYPE * Vy_grid, NUMERIC_TYPE * Vx_max_grid, NUMERIC_TYPE * Vy_max_grid,
	const NUMERIC_TYPE * SGC_BankFullHeight_grid,
	const NUMERIC_TYPE *dem_grid,
	const NUMERIC_TYPE *g_friction_sq_x_grid, const NUMERIC_TYPE *g_friction_sq_y_grid,
	const NUMERIC_TYPE *friction_x_grid, const NUMERIC_TYPE *friction_y_grid,
	const NUMERIC_TYPE *dx_col, const NUMERIC_TYPE *dy_col, const NUMERIC_TYPE *cell_area_col,
	const NUMERIC_TYPE *Fp_xwidth, const NUMERIC_TYPE *Fp_ywidth,

	const SubGridRowList * sub_grid_layout_rows, SubGridState * sub_grid_state_rows,
	const SubGridRowList * sub_grid_layout_blocks, SubGridState * sub_grid_state_blocks,

	TimeSeries * evap_time_series,
	NetCDFVariable * evap_grid,
	TimeSeries * rain_time_series,
	const NUMERIC_TYPE *rain_grid,
	const NUMERIC_TYPE *dist_infil_grid,

	WetDryRowBound* wet_dry_bounds,
	PointSourceRowList * ps_layout, BoundaryCondition * boundary_cond,
	WeirLayout * weir_weirs, WeirLayout * weir_bridges,
	RouteDynamicList * route_dynamic_list,
	const NUMERIC_TYPE *route_V_ratio_per_sec_qx, const NUMERIC_TYPE * route_V_ratio_per_sec_qy,

	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,

	const SGCprams * SGCptr,
	Files *Fptr,
	Stage *Locptr,
	DamData *Damptr,
	SuperGridLinksList *Super_linksptr,
	//DynamicRain<> & dynamic_rain, removed JCN
	NUMERIC_TYPE ** tmp_thread_data,
	const int verbose)
{
	NUMERIC_TYPE * tmp_row = tmp_thread_data[omp_get_thread_num()];

	//each thread clears the tmp_row, before using in SGC2_UpdateQ_block
	memset(tmp_row, 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);

	/// nowait means that some threads may still continue after this function returns.
	/// note the #pragma omp barrier is required before the using or updating q
#pragma omp for schedule(static) nowait
	for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
	{
		const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
		const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

		for (int j = start_y; j < end_y; j++)
		{
			const IndexRange dem_data_bound = wet_dry_bounds->dem_data[j];
			const int grid_row_index = j * grid_cols_padded;
			const NUMERIC_TYPE row_dx = dx_col[j];
			const NUMERIC_TYPE row_dy = dy_col[j];
			const NUMERIC_TYPE row_cell_area = cell_area_col[j];

			// first process Qx
			{
				int row_start_x = max(dem_data_bound.start, wet_dry_bounds->fp_h[j].start - 1); // start one cell to the left (flow from cell on the right)
				int row_end_x = min(wet_dry_bounds->fp_h[j].end, dem_data_bound.end - 1);

				int row_start_x_prev = max(dem_data_bound.start, wet_dry_bounds->fp_h_prev[j].start - 1);
				int row_end_x_prev = min(wet_dry_bounds->fp_h_prev[j].end, dem_data_bound.end - 1);

				SGC2_UpdateQx_row(grid_cols,
					grid_row_index,
					row_start_x_prev, row_end_x_prev,
					row_start_x, row_end_x,
					depth_thresh, row_dx, Fp_ywidth,
					g, delta_time, curr_time,
					tmp_row,
					dem_grid, h_grid,
					g_friction_sq_x_grid, Qx_grid, Qx_old_grid, Parptr->max_Froude);

				if (Statesptr->routing == ON)
				{
					int route_count = SGC2_UpdateRouteQ_row(grid_row_index, row_start_x, row_end_x,
						row_cell_area, delta_time, Parptr->RouteSfThresh,
						1,
						tmp_row, h_grid, dem_grid,
						route_V_ratio_per_sec_qx, Qx_grid, Qx_old_grid,
						route_dynamic_list->route_list_i_lookup_qx);
					route_dynamic_list->row_route_qx_count[j] = route_count;
				}
				else if (Statesptr->diffusive_switch == ON)
				{
					SGC2_UpdateDiffusiveQ_row(grid_row_index, row_start_x, row_end_x,
						row_dy, row_dx,
						depth_thresh, Parptr->DiffusiveFroudeThresh, Solverptr->MaxHflow, g,
						1,
						h_grid, dem_grid, friction_x_grid, Qx_grid, Qx_old_grid);
				}
				if ((curr_time >= Parptr->SaveTotal && Statesptr->voutput == ON) || // calculate v at save interval only
					Statesptr->voutput_max == ON) // v_max_output or hazard require velocity to be calculated each step
				{
					// velocity calculated before sub grid or weir updates
					// only calculates flood-plain value
					SGC2_UpdateVelocity_row(grid_row_index, row_start_x_prev, row_end_x_prev, row_start_x, row_end_x,
						1,
						row_dy, h_grid, dem_grid, Qx_grid, Vx_grid, Vx_max_grid);
				}
			}
			bool not_last_row = (j < grid_rows - 1);
			// PFU: this is processing Qy
			//skip last row
			if (not_last_row)
			{
				int row_start_y = min(wet_dry_bounds->fp_h[j].start, wet_dry_bounds->fp_h[j + 1].start);
				int row_end_y = max(wet_dry_bounds->fp_h[j].end, wet_dry_bounds->fp_h[j + 1].end);

				int row_start_y_prev = min(wet_dry_bounds->fp_h_prev[j].start, wet_dry_bounds->fp_h_prev[j + 1].start);
				int row_end_y_prev = max(wet_dry_bounds->fp_h_prev[j].end, wet_dry_bounds->fp_h_prev[j + 1].end);

				SGC2_UpdateQy_row(grid_cols,
					grid_row_index,
					row_start_y_prev, row_end_y_prev,
					row_start_y, row_end_y,
					grid_cols_padded,
					depth_thresh, Fp_xwidth, row_dy,
					g, delta_time, curr_time,
					tmp_row,
					dem_grid, h_grid,
					g_friction_sq_y_grid, Qy_grid, Qy_old_grid, Parptr->max_Froude);

				if (Statesptr->routing == ON)
				{
					int route_count = SGC2_UpdateRouteQ_row(grid_row_index, row_start_y, row_end_y,
						row_cell_area, delta_time, Parptr->RouteSfThresh,
						grid_cols_padded,
						tmp_row, h_grid, dem_grid,
						route_V_ratio_per_sec_qy, Qy_grid, Qy_old_grid,
						route_dynamic_list->route_list_i_lookup_qy);
					route_dynamic_list->row_route_qy_count[j] = route_count;
				}
				else if (Statesptr->diffusive_switch == ON)
				{
					SGC2_UpdateDiffusiveQ_row(grid_row_index, row_start_y, row_end_y,
						row_dx, row_dy,
						depth_thresh, Parptr->DiffusiveFroudeThresh, Solverptr->MaxHflow, g,
						grid_cols_padded,
						h_grid, dem_grid, friction_y_grid, Qy_grid, Qy_old_grid);
				}
				if ((curr_time >= Parptr->SaveTotal && Statesptr->voutput == ON) || // calculate v at save interval only
					Statesptr->voutput_max == ON) // v_max_output or hazard require velocity to be calculated each step
				{
					// velocity calculated before sub grid or weir updates
					// only calculates flood-plain value
					SGC2_UpdateVelocity_row(grid_row_index, row_start_y_prev, row_end_y_prev, row_start_y, row_end_y,
						grid_cols_padded,
						row_dx, h_grid, dem_grid, Qy_grid, Vy_grid, Vy_max_grid);
				}
			}

#if _SGM_BY_BLOCKS == 0
			ProcessSubGridQBlock(j, grid_cols_padded, depth_thresh, delta_time, g, sub_grid_layout_rows, sub_grid_state_rows, SGCptr, h_grid,
				wet_dry_bounds, Qx_grid, Qy_grid, Qx_old_grid, Qy_old_grid, Parptr->max_Froude);
			if (curr_time >= Parptr->SaveTotal && Statesptr->SGCvoutput == ON)
			{
				SGC2_UpdateVelocitySubGrid_block(j, grid_cols_padded, depth_thresh, delta_time, sub_grid_layout_rows, sub_grid_state_rows, SGCptr, h_grid);
			}
#endif
			if (weir_weirs->row_cols_padded != 0)
			{
				SGC2_UpdateWeirsFlow_row(j, grid_cols, grid_rows, grid_cols_padded,
					depth_thresh, delta_time, h_grid, volume_grid, Qx_grid, Qy_grid,
					wet_dry_bounds, weir_weirs);
			}
		}
		//printf("Q done %d\n", omp_get_thread_num());
	}
#if _SGM_BY_BLOCKS == 1
	if (sub_grid_layout_blocks->row_cols_padded > 0)
	{
#pragma omp barrier
#pragma omp for schedule(static) nowait
		for (int bi = 0; bi < wet_dry_bounds->block_count; bi++)
		{
			ProcessSubGridQBlock(bi, grid_cols_padded, depth_thresh, delta_time, g, sub_grid_layout_blocks, sub_grid_state_blocks, SGCptr, h_grid,
				wet_dry_bounds, Qx_grid, Qy_grid, Qx_old_grid, Qy_old_grid);
			if (curr_time >= Parptr->SaveTotal && Statesptr->SGCvoutput == ON)
			{
				SGC2_UpdateVelocitySubGrid_block(j, grid_cols_padded, depth_thresh, delta_time, sub_grid_layout_blocks, sub_grid_state_blocks, SGCptr, h_grid);
			}
		}
	}	
	}
	{
		
#endif
	
	// this block of code is executed by the first thread that finishes it's updateQ (with nowait clause) 
	// nowait clause added - disable the implicit barrier, since and explicit barrier used (stops both the update q threads as well as the single thread)
#pragma omp single nowait 
	{
		//printf("overlapped single %d\n", omp_get_thread_num());
		reduce_Hmax = C(0.0);
		reduce_evap_loss = C(0.0);
		reduce_rain_total = C(0.0);
		reduce_Qpoint_timestep_pos = C(0.0);
		reduce_Qpoint_timestep_neg = C(0.0);

		reduce_flood_area = C(0.0);
		reduce_domain_volume = C(0.0);

		evap_deltaH_step = C(0.0);
		rain_deltaH_step = C(0.0);

		// pre-calculate the TimeSeries interpolation for all point sources
		// * prevent any threading issues if a TimeSeries is used on multiple points
		// * reduce code complexity in calculating point source
		for (int j = 0; j < grid_rows; j++)
		{
			const int ps_count = ps_layout->ps_row_count[j];
			const int row_cols_padded = ps_layout->row_cols_padded;
			const int row_start = j * row_cols_padded;
			WaterSource ps_info = ps_layout->ps_info;
			for (int i = 0; i < ps_count; i++)
			{
				int ws_index = row_start + i;
				if (ps_info.timeSeries[ws_index] != NULL)
				{
					ps_info.Val[ws_index] = InterpolateTimeSeries(ps_info.timeSeries[ws_index], curr_time);
				}
			}
		}

		// SGC_BCs updates the Q flow at the boundaries
		// Updates Q_x and Q_y on the floodplain and also updates volume for any sub-grid flow
		SGC2_BCs(grid_cols, grid_rows, grid_cols_padded, delta_time, curr_time, depth_thresh, g,
			dx_col, dy_col, h_grid, Qx_grid, Qy_grid, Qx_old_grid, Qy_old_grid,
			wet_dry_bounds,
			Statesptr, Parptr, boundary_cond, SGCptr, Parptr->max_Froude);

		if (Statesptr->calc_evap == ON)
		  {
		    // evap_deltaH_step = InterpolateTimeSeries(evap_time_series, curr_time); //constant rate across whole floodplain
		    for (int j = 0; j < evap_time_series->count; j++) {
		      if (evap_time_series->time[j] > curr_time) {
			evap_deltaH_step = evap_time_series->value[j];
			break;
		      }
		    }
		    evap_deltaH_step *= delta_time;
		    for (int j = 0; j < grid_cols_padded * grid_rows; j++) {
		      evap_grid->data[j] = evap_deltaH_step;
		    }
		  }
		else if (Statesptr->calc_evap == TIME_SPACE)
		  {
		    for (int j =0; j < grid_cols_padded * grid_rows; j++) {
		      evap_grid->data[j] *= delta_time;
		      if (evap_grid->data[j] > 0.0)
			evap_deltaH_step = evap_grid->data[j];
		    }
		  }

		if (Statesptr->rainfall == ON)
		{
			rain_deltaH_step = InterpolateTimeSeries(rain_time_series, curr_time); //constant rate across whole floodplain
			rain_deltaH_step *= delta_time;

		}
		//printf("overlapped done %d\n", omp_get_thread_num());
	}

	//printf("before barrier done %d\n", omp_get_thread_num());

	// Toby suggest DamFlowVolume go here!!! FEOL
#pragma omp barrier // ensure all threads have finished their updateQ (nowait) and single section
	
	if (Statesptr->DamMode == ON)
	{
#pragma omp single
		{
			// Call Dam Function after all Q have being calculated!  FEOL
			SGC2_UpdateDamFlowVolume(grid_cols, grid_rows, grid_cols_padded, depth_thresh, delta_time, curr_time, h_grid, volume_grid, Qx_grid, Qy_grid, Damptr, SGCptr, g, Parptr->max_Froude);
		}
		//End Dam FEOL
	}
	// supergrid channel links... currently in OMP single section
	if (Statesptr->ChanMaskRead == ON)
	{
#pragma omp single
		{
			SGC2_CalcLinksQ(Super_linksptr, volume_grid, h_grid, delta_time, g, depth_thresh, Parptr->max_Froude, wet_dry_bounds);
		}
	}
	//printf("after barrier done %d\n", omp_get_thread_num());
	if (Statesptr->save_stages == ON &&
		Statesptr->voutput_stage == ON &&
		curr_time >= Parptr->MassTotal)
	{
#pragma omp single
		{


			fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", curr_time);
			for (int i = 0; i < Locptr->Nstages; i++)
			{
				if (Locptr->stage_check[i] == 1)
				{
					int y = Locptr->stage_grid_y[i];
					int grid_index = Locptr->stage_grid_x[i] + y * grid_cols_padded;
					NUMERIC_TYPE Vx_west = (i != 0) ? SGC2_CalculateVelocity(grid_index - 1, grid_index, Qx_grid, h_grid, dem_grid, dy_col[y]) : C(0.0);
					NUMERIC_TYPE Vx_east = (i < grid_cols - 1) ? SGC2_CalculateVelocity(grid_index, grid_index + 1, Qx_grid, h_grid, dem_grid, dy_col[y]) : C(0.0);

					NUMERIC_TYPE Vy_north = (y != 0) ? SGC2_CalculateVelocity(grid_index - grid_cols_padded, grid_index, Qy_grid, h_grid, dem_grid, dx_col[y]) : C(0.0);
					NUMERIC_TYPE Vy_south = (y < grid_rows - 1) ? SGC2_CalculateVelocity(grid_index, grid_index + grid_cols_padded, Qy_grid, h_grid, dem_grid, dx_col[y]) : C(0.0);

					// v not calculated at each timestep 
					//NUMERIC_TYPE Vx_west = Vx_grid[grid_index];
					//NUMERIC_TYPE Vx_east = Vx_grid[grid_index + 1];

					//NUMERIC_TYPE Vy_north = Vy_grid[grid_index];
					//NUMERIC_TYPE Vy_south = Vy_grid[grid_index + grid_cols_padded];

					fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", SQRT(POW(getmax(FABS(Vx_east), FABS(Vx_west)), C(2.0)) + POW(getmax(FABS(Vy_north), FABS(Vy_south)), C(2.0))));
				}
				else
					fprintf(Fptr->vel_fp, "-\t");
			}
			fprintf(Fptr->vel_fp, "\n");
			fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		}
	}
	
	if (Statesptr->routing_mass_check == ON || weir_bridges->row_cols_padded > 0 || Statesptr->hazard == ON)
	{
		/// This update must be performed after UpdateQ as it requires all Q's to be updated
		/// Since it also modifies Q, it must be before any reading of Q (in SGC2_UpdateVolumeHeight)
#pragma omp for schedule(static)
		//for (int j = 0; j < grid_rows; j++)
		for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
		{
			const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
			const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

			for (int j = start_y; j < end_y; j++)
			{
				int grid_row_index = j * grid_cols_padded;
				if (Statesptr->routing_mass_check == ON)
				{
					SGC2_CorrectRouteFlow_row(j, grid_row_index, grid_cols, grid_rows, grid_cols_padded, delta_time,
						route_dynamic_list, volume_grid,
						Qx_grid, Qy_grid);
				}
				if (weir_bridges->row_cols_padded > 0)
				{
#if _SGM_BY_BLOCKS == 1
					SGC2_UpdateBridgesFlow_row(j, grid_cols, grid_rows, grid_cols_padded,
						delta_time, curr_time, depth_thresh, g, dx_col, dy_col, h_grid, Qx_grid, Qy_grid,
						wet_dry_bounds, sub_grid_state_blocks, weir_bridges);
#else
					SGC2_UpdateBridgesFlow_row(j, grid_cols, grid_rows, grid_cols_padded,
						delta_time, curr_time, depth_thresh, g, dx_col, dy_col, h_grid, Qx_grid, Qy_grid,
						wet_dry_bounds, sub_grid_state_rows, weir_bridges, Parptr->max_Froude);
#endif
				}
				if (Statesptr->hazard == ON)
				{
					SGC2_UpdateHazard_row(j, grid_row_index, grid_cols, grid_rows, grid_cols_padded,
						Vx_grid, Vy_grid, h_grid, dem_grid,
						dx_col[j], dy_col[j], SGC_BankFullHeight_grid, maxVc_grid, maxVc_height_grid, maxHazard_grid, wet_dry_bounds);
				}
			}
		}
		//#ifdef _DEBUG
		//			int total_routes = 0;
		//			for (int j = 0; j < grid_rows; j++)
		//			{
		//				total_routes += route_dynamic_list->row_route_qx_count[j];
		//				total_routes += route_dynamic_list->row_route_qy_count[j];
		//			}
		//			printf("Routes %d Rain %12.4e %.2" NUM_FMT" %% of cells\n", total_routes, rain_deltaH_step, (NUMERIC_TYPE)total_routes / (grid_cols*grid_rows) * 100);
		//#endif
	}

	//each thread clears the tmp_row, before using in SGC2_UpdateVolumeHeight_block
	memset(tmp_row, 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma omp for reduction (+:reduce_evap_loss, reduce_rain_total, reduce_Qpoint_timestep_pos, reduce_Qpoint_timestep_neg) schedule(static)
#else
#pragma omp for reduction(+:reduce_evap_loss, reduce_rain_total, reduce_Qpoint_timestep_pos, reduce_Qpoint_timestep_neg) reduction( max : reduce_Hmax) schedule(static)
#endif

	for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
	{
		VolumeHeightUpdateInfo update_info;
		NUMERIC_TYPE block_Hmax=C(0.0);

		block_Hmax = SGC2_UpdateVolumeHeight_block(block_index, grid_cols, grid_rows, grid_cols_padded,
			curr_time, delta_time, depth_thresh, g,
			evap_deltaH_step, evap_grid, rain_deltaH_step, rain_grid, dist_infil_grid,
			wet_dry_bounds, tmp_row,
			Qx_grid, Qy_grid,
			cell_area_col, dx_col, dy_col,
			dem_grid,
			sub_grid_layout_rows, sub_grid_state_rows,
			ps_layout, route_dynamic_list,

			h_grid, volume_grid,
			SGC_BankFullHeight_grid,
			Statesptr, Parptr, Solverptr, SGCptr, &update_info);

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma omp critical
#endif
		{
			if (block_Hmax > reduce_Hmax)
				reduce_Hmax = block_Hmax;
		}

		reduce_evap_loss += update_info.evap_loss;
		reduce_rain_total += update_info.rain_total;
		reduce_Qpoint_timestep_pos += update_info.Qpoint_timestep_pos;
		reduce_Qpoint_timestep_neg += update_info.Qpoint_timestep_neg;
	}

	// Update time of initial flood inundation (in hours) and total inundation time (in seconds)
	if ((Statesptr->reset_timeinit == ON) && (curr_time > Parptr->reset_timeinit_time))
	{
		//reset the time of initial inundation if called for in parameter file
		Statesptr->reset_timeinit = OFF;

#pragma omp for schedule(static)
		for (int j = 0; j < grid_rows; j++)
		{
			for (int i = 0; i < grid_cols; i++)
			{
				int index = i + j*grid_cols_padded;
				initHtm_grid[index] = (NULLVAL);
			}
		}
		if (verbose == ON)
		{
#pragma omp single
			printf("\n Time of initial inundation reset \n");
		}
	}

	NUMERIC_TYPE time_next = curr_time + delta_time;
	if (time_next >= Parptr->MassTotal || Statesptr->mint_hk == OFF)
	{
		NUMERIC_TYPE delta_time_hours = delta_time / C(3600.0);
		NUMERIC_TYPE current_time_hours = time_next / C(3600.0);

#pragma omp for schedule(static) nowait
		for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
		{
			SGC2_Inundation_block(block_index, grid_cols, grid_rows, grid_cols_padded,
				depth_thresh,
				current_time_hours,
				delta_time_hours,
				h_grid,
				sub_grid_layout_rows, sub_grid_state_rows,
				wet_dry_bounds,
				initHtm_grid,
				totalHtm_grid,
				maxH_grid,
				maxHtm_grid);
		}
		if (time_next >= Parptr->MassTotal)
		{
#pragma omp for reduction ( + : reduce_flood_area, reduce_domain_volume) schedule(static) nowait
			for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
			{
				NUMERIC_TYPE out_flood_area, out_domain_volume;
				SGC2_DomainVolumeAndFloodArea_block(block_index, grid_cols, grid_rows, grid_cols_padded, depth_thresh,
					wet_dry_bounds,
					h_grid, cell_area_col, volume_grid, &out_flood_area, &out_domain_volume);
				reduce_flood_area += out_flood_area;
				reduce_domain_volume += out_domain_volume;
			}
		}

		// ensure inundation & volume/area threads finished 
#pragma omp barrier 
	}


#if _BALANCE_TYPE == 1
	SGC2_UpdateLoadBalance(grid_rows, grid_cols_padded, sub_grid_layout_rows, wet_dry_bounds);
#endif

	//return reduce_Hmax;
}

//-----------------------------------------------------------------------------
// ITERATE THROUGH TIME STEPS
///
///
///
/*! \fn void Fast_IterateLoop(...)
\brief
\param
\param
\param h_grid matrix grid_cols_padded*grid_rows
\param cell_area_col column 1*grid_rows
*/
void Fast_IterateLoop(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	NUMERIC_TYPE *h_grid, NUMERIC_TYPE *volume_grid, 
	NUMERIC_TYPE *Qx_grid, NUMERIC_TYPE *Qy_grid, NUMERIC_TYPE *Qx_old_grid, NUMERIC_TYPE *Qy_old_grid,
	NUMERIC_TYPE *maxH_grid, NUMERIC_TYPE *maxHtm_grid, NUMERIC_TYPE *initHtm_grid, NUMERIC_TYPE *totalHtm_grid,
	NUMERIC_TYPE *maxVc_grid, NUMERIC_TYPE *maxVc_height_grid, NUMERIC_TYPE *maxHazard_grid,
	NUMERIC_TYPE *Vx_grid, NUMERIC_TYPE *Vy_grid, NUMERIC_TYPE *Vx_max_grid, NUMERIC_TYPE *Vy_max_grid,
	const NUMERIC_TYPE *dem_grid,
	const NUMERIC_TYPE *g_friction_sq_x_grid, const NUMERIC_TYPE *g_friction_sq_y_grid,
	const NUMERIC_TYPE *friction_x_grid, const NUMERIC_TYPE *friction_y_grid,
	const NUMERIC_TYPE *dx_col, const NUMERIC_TYPE *dy_col, const NUMERIC_TYPE *cell_area_col,
	const NUMERIC_TYPE *Fp_xwidth, const NUMERIC_TYPE *Fp_ywidth,

	const SubGridRowList * sub_grid_layout_rows,
	SubGridState * sub_grid_state_rows,
	const SubGridRowList * sub_grid_layout_blocks,
	SubGridState * sub_grid_state_blocks,

	const NUMERIC_TYPE * SGC_BankFullHeight_grid,

	TimeSeries * evap_time_series,
	NetCDFVariable * evap_grid,
	TimeSeries * rain_time_series,
	NUMERIC_TYPE *rain_grid,
	const NUMERIC_TYPE *dist_infil_grid,

	WetDryRowBound* wet_dry_bounds,
	PointSourceRowList * ps_layout, BoundaryCondition * boundary_cond,
	WeirLayout * weirs_weirs, WeirLayout * weirs_bridges,
	RouteDynamicList * route_dynamic_list,
	const NUMERIC_TYPE *route_V_ratio_per_sec_qx, const NUMERIC_TYPE * route_V_ratio_per_sec_qy,

	SuperGridLinksList *Super_linksptr,
	Fnames *Fnameptr,
	Files *Fptr,
	Stage *Locptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	DamData *Damptr,
	SGCprams * SGCptr,
	NUMERIC_TYPE ** tmp_thread_data,
#ifdef RESULT_CHECK
	Arrays * Arrptr, // only for compare results
	BoundCs * BCptr, // only for compare results
	ChannelSegmentType *ChannelSegments,
	vector<ChannelSegmentType> *ChannelSegmentsVecPtr,
#endif

	const int verbose)
{
	int tstep_counter = -1;   // start at -1 so that in first run through we calculate river
	int steadyCount = 0;
	NUMERIC_TYPE tstep_channel = C(0.0); // channel timestep
	NUMERIC_TYPE Previous_t;      // previous time channel was calculated
	NUMERIC_TYPE discharge = C(0.0); // value of discharge for virtual gauge output
	NUMERIC_TYPE loss; //temp variable to keep track of losses since last mass interval
	NUMERIC_TYPE Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;

	const NUMERIC_TYPE depth_thresh = Solverptr->DepthThresh;
	const NUMERIC_TYPE g = Solverptr->g;

	// tmp data for preparing data to write as output
	NUMERIC_TYPE * tmp_grid1 = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	NUMERIC_TYPE * tmp_grid2 = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	NUMERIC_TYPE * tmp_grid3 = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * grid_rows * grid_cols_padded);
	memset(tmp_grid2, 0, sizeof(NUMERIC_TYPE) * grid_rows * grid_cols_padded);
	memset(tmp_grid3, 0, sizeof(NUMERIC_TYPE) * grid_rows * grid_cols_padded);

	{
		NUMERIC_TYPE reduce_flood_area = C(0.0);
		NUMERIC_TYPE reduce_domain_volume = C(0.0);
#pragma omp parallel for default(shared) reduction ( + : reduce_flood_area, reduce_domain_volume) schedule(static)
		for (int block_index = 0; block_index < wet_dry_bounds->block_count; block_index++)
		{
			NUMERIC_TYPE block_flood_area, block_domain_volume;
			SGC2_DomainVolumeAndFloodArea_block(block_index, grid_cols, grid_rows, grid_cols_padded, depth_thresh,
				wet_dry_bounds,
				h_grid, cell_area_col, volume_grid, &block_flood_area, &block_domain_volume);
			reduce_flood_area += block_flood_area;
			reduce_domain_volume += block_domain_volume;
		}
		Solverptr->vol1 = reduce_domain_volume;
	}

	// set previous time to one timestep backwards so that first river calcs uses Solverptr->Tstep for 1st iteration
	// this is because of the way the timestep is calculated as the time difference from the last time the river was run
	Previous_t = Solverptr->t - Solverptr->Tstep;

#if defined (__INTEL_COMPILER) && _PROFILE_MODE > 0
	printf("Intel profiler resume\n");
	__itt_resume();
#endif

	//NUMERIC_TYPE * write_params = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * 1 * 1);

	struct timeval timstr;      /* structure to hold elapsed time */
	double processing_start_time, processing_end_time;
	gettimeofday(&timstr, NULL);
	processing_start_time = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
	

	time_t loop_start;
	time(&loop_start);

	// initialise the h wet_dry_bounds
	SGC2_InitHBounds(grid_cols, grid_rows, grid_cols_padded, depth_thresh, sub_grid_layout_rows,
		sub_grid_state_rows,
		cell_area_col,
		h_grid, volume_grid, wet_dry_bounds, SGCptr);
	SGC2_UpdateLoadBalance(grid_rows, grid_cols_padded, sub_grid_layout_rows, wet_dry_bounds);

	// initalise distributed rainfall
	DynamicRain<> dynamic_rain(Fnameptr->dynamicrainfilename, verbose);

#if defined(RESULT_CHECK)	
	States Statesptr2;
	Pars Parptr2;
	Solver Solverptr2;
	SGCprams SGCptr2;
	NUMERIC_TYPE limit = C(1e-5);
	NUMERIC_TYPE qlimit = C(1e-5);
#endif
#ifdef _MSC_VER
	//Solverptr->t = 54 * 79488;
#endif

	NUMERIC_TYPE curr_time = Solverptr->t;

#if defined(RESULT_CHECK)	
	memcpy(&Statesptr2, Statesptr, sizeof(States));
	memcpy(&Parptr2, Parptr, sizeof(Pars));
	memcpy(&Solverptr2, Solverptr, sizeof(Solver));
	memcpy(&SGCptr2, SGCptr, sizeof(SGCprams));

	for (int j = 0; j < grid_rows; j++)
	{
		for (int i = 0; i < grid_cols; i++)
		{
			int index_old = i + j * grid_cols;
			int index_new = i + j * grid_cols_padded;

			NUMERIC_TYPE h_tmp = h_grid[index_new] + SGC_BankFullHeight_grid[index_new];
			NUMERIC_TYPE delta = h_tmp - Arrptr->H[index_old];
			if (FABS(delta) > limit)
			{
				printf("error h (%d,%d) Solverptr->itCount == %ld \n", i, j, Solverptr->itCount);
				printf("h new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", h_tmp, Arrptr->H[index_old], delta);

				NUMERIC_TYPE tdelta = curr_time - Solverptr2.t;
				printf("t delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", curr_time, Solverptr2.t, tdelta);

				delta = volume_grid[index_new] - Arrptr->SGCVol[index_old];
				printf("vol new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", volume_grid[index_new], Arrptr->SGCVol[index_old], delta);

				//printf("vol new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", volume_grid[index_new], Arrptr->SGCbfH[index_old], delta);

				printf("H bound %d => %d\n", wet_dry_bounds->fp_h[j].start, wet_dry_bounds->fp_h[j].end);
				printf("H prev  %d => %d\n", wet_dry_bounds->fp_h_prev[j].start, wet_dry_bounds->fp_h_prev[j].end);
				printf("vol     %d => %d\n", wet_dry_bounds->fp_vol[j].start, wet_dry_bounds->fp_vol[j].end);

				printf("--\n");
			}
		}
	}

#endif

	int timestep_tripped = 0;
	int error_tripped = 0;

	int stop_loop = OFF;
	itCount = Solverptr->itCount;

#pragma omp parallel default(shared)
	{
		// main iteration loop
		while (curr_time < Solverptr->Sim_Time && stop_loop == OFF)
		{
#if defined(RESULT_CHECK)	
#pragma omp single
			{
				memcpy(&Statesptr2, Statesptr, sizeof(States));
				memcpy(&Parptr2, Parptr, sizeof(Pars));
				memcpy(&Solverptr2, Solverptr, sizeof(Solver));
				memcpy(&SGCptr2, SGCptr, sizeof(SGCprams));
			}
#endif

#if (_NETCDF == 1)			
#pragma omp single
			{
			  // read in evaporation grid
			  if (Statesptr->calc_evap == TIME_SPACE) {
			    read_file_netCDF(evap_grid, curr_time);
			    for (int j = 0; j < grid_cols_padded * grid_rows; j++)
			      evap_grid->data[j] /= (86400. * 1000.);  // convert from mm/day to m/s
			  }

			  // read dynamic rain grid and write to rain grid

			  dynamic_rain.update_rain_grid_SGM(curr_time, rain_grid, dem_grid, wet_dry_bounds, dx_col, dy_col, Parptr->tly, grid_rows, grid_cols_padded, cell_area_col);

			}
#endif

			// sub grid floodplain models
			NUMERIC_TYPE delta_time;
			delta_time = Solverptr->SGCtmpTstep;


			Do_Update(grid_cols, grid_rows, grid_cols_padded, delta_time,
				curr_time, depth_thresh, g, h_grid, volume_grid, Qx_grid, Qy_grid, Qx_old_grid, Qy_old_grid,
				initHtm_grid, maxHtm_grid, totalHtm_grid, maxH_grid,
				maxVc_grid, maxVc_height_grid, maxHazard_grid,
				Vx_grid, Vy_grid, Vx_max_grid, Vy_max_grid,
				SGC_BankFullHeight_grid,
				dem_grid, g_friction_sq_x_grid, g_friction_sq_y_grid,
				friction_x_grid, friction_y_grid,
				dx_col, dy_col, cell_area_col,
				Fp_xwidth,Fp_ywidth,
				sub_grid_layout_rows, sub_grid_state_rows, sub_grid_layout_blocks, sub_grid_state_blocks,
				evap_time_series, evap_grid, rain_time_series, rain_grid, dist_infil_grid, wet_dry_bounds, ps_layout, boundary_cond,
				weirs_weirs, weirs_bridges,
				route_dynamic_list, route_V_ratio_per_sec_qx, route_V_ratio_per_sec_qy,
				Statesptr, Parptr, Solverptr, SGCptr, Fptr, Locptr, Damptr, Super_linksptr,
				tmp_thread_data,
				verbose);


#pragma omp single
			{
				Solverptr->Tstep = delta_time;
#if _DISABLE_WET_DRY == 1
				for (int j = 0; j < grid_rows; j++)
				{
					wet_dry_bounds->fp_h[j].start = 0;
					wet_dry_bounds->fp_h[j].end = grid_cols;
				}
#endif
				for (int j = 0; j < grid_rows; j++)
				{
					wet_dry_bounds->fp_vol[j].start = wet_dry_bounds->fp_h[j].start;
					wet_dry_bounds->fp_vol[j].end = wet_dry_bounds->fp_h[j].end;
				}
				Parptr->EvapTotalLoss += reduce_evap_loss;
				Parptr->RainTotalLoss += reduce_rain_total;
				ps_layout->Qpoint_pos = reduce_Qpoint_timestep_pos / delta_time;
				ps_layout->Qpoint_neg = reduce_Qpoint_timestep_neg / delta_time;

				// Point sources
				boundary_cond->Qin += ps_layout->Qpoint_pos;
				boundary_cond->Qout -= ps_layout->Qpoint_neg;

				// calculate volume in and volume out
				boundary_cond->VolInMT += boundary_cond->Qin*delta_time;
				boundary_cond->VolOutMT += boundary_cond->Qout*delta_time;

				Solverptr->vol2 = reduce_domain_volume;
				Solverptr->FArea = reduce_flood_area;


				/*if (Solverptr->Tstep < 5)
				{
				printf("error TimeStep @ Solverptr->itCount == %ld\n", Solverptr->itCount);
				timestep_tripped++;
				Solverptr->Tstep = 5;
				}
				if(timestep_tripped>0)
				{
				timestep_tripped++;
				write_regular_output(Fnameptr, grid_cols, grid_rows, grid_cols_padded, depth_thresh, tmp_grid1,
				h_grid, dem_grid,
				Qx_grid, Qy_grid,
				SGC_BankFullHeight_grid, maxH_grid
				sub_grid_layout,
				sub_grid_state,

				dx_col, dy_col,
				Solverptr, Statesptr, Parptr, boundary_cond, SGCptr,
				10000 + timestep_tripped,
				ON, ON, ON);

				if(timestep_tripped>10)
				{
				exit(0);
				}
				}*/


#if defined(RESULT_CHECK)	
				Solverptr2.Tstep = Solverptr2.SGCtmpTstep;
				NUMERIC_TYPE Hmax2 = Do_Update_old(&Statesptr2, &Parptr2, &Solverptr2, Arrptr, &SGCptr2, BCptr, ChannelSegments, ChannelSegmentsVecPtr);

				//if (Solverptr->itCount >= 345000)
				{
					int step_error = 0;
					for (int j = 0; j < (grid_rows + 1); j++)
					{
						for (int i = 0; i < grid_cols; i++)
						{
							int index_old_q = i + j * (grid_cols + 1);
							int index_new_q = i + j * grid_cols_padded;

							NUMERIC_TYPE delta = Qy_grid[index_new_q] - Arrptr->Qy[index_old_q];
							if (FABS(delta) > qlimit)
							{
								printf("error Qy (%d,%d) Solverptr->itCount == %ld\n", i, j, Solverptr->itCount);

								NUMERIC_TYPE tdelta = curr_time - Solverptr2.t;
								printf("t delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", curr_time, Solverptr2.t, tdelta);

								NUMERIC_TYPE hmaxdelta = reduce_Hmax - Hmax2;
								printf("hmax delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", reduce_Hmax, Hmax2, hmaxdelta);

								NUMERIC_TYPE q_old_m3 = (Arrptr->Qyold[index_old_q] * dx_col[j]);

								printf("Qy    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_grid[index_new_q], Arrptr->Qy[index_old_q], delta);
								NUMERIC_TYPE old_delta = Qy_old_grid[index_new_q] - q_old_m3;
								printf("Qyold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_old_grid[index_new_q], q_old_m3, old_delta);

								NUMERIC_TYPE new_non_fp = Qy_grid[index_new_q] - Qy_old_grid[index_new_q];
								NUMERIC_TYPE old_non_fp = Arrptr->Qy[index_old_q] - q_old_m3;
								NUMERIC_TYPE non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								printf("--\n");
							}
						}
					}

					for (int j = 0; j < grid_rows; j++)
					{
						for (int i = 0; i < (grid_cols + 1); i++)
						{
							int index_old_q = i + j * (grid_cols + 1);
							int index_new_q = i + j * grid_cols_padded;

							NUMERIC_TYPE delta = Qx_grid[index_new_q] - Arrptr->Qx[index_old_q];
							if (FABS(delta) > qlimit)
							{
								printf("error Qx (%d,%d) Solverptr->itCount == %ld\n", i, j, Solverptr->itCount);

								NUMERIC_TYPE tdelta = curr_time - Solverptr2.t;
								printf("t delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", curr_time, Solverptr2.t, tdelta);

								NUMERIC_TYPE hmaxdelta = reduce_Hmax - Hmax2;
								printf("hmax delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", reduce_Hmax, Hmax2, hmaxdelta);

								NUMERIC_TYPE q_old_m3 = (Arrptr->Qxold[index_old_q] * dy_col[j]);

								printf("Qx    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_grid[index_new_q], Arrptr->Qx[index_old_q], delta);
								NUMERIC_TYPE old_delta = Qx_old_grid[index_new_q] - q_old_m3;
								printf("Qxold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_old_grid[index_new_q], q_old_m3, old_delta);

								NUMERIC_TYPE new_non_fp = Qx_grid[index_new_q] - Qx_old_grid[index_new_q];
								NUMERIC_TYPE old_non_fp = Arrptr->Qx[index_old_q] - q_old_m3;
								NUMERIC_TYPE non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								printf("--\n");
							}
						}
					}

					if (Statesptr->hazard == ON)
					{
						for (int j = 0; j < grid_rows; j++)
						{
							for (int i = 0; i < grid_cols; i++)
							{
								int index_old = i + j * grid_cols;
								int index_new = i + j * grid_cols_padded;

								NUMERIC_TYPE delta = maxVc_grid[index_new] - Arrptr->maxVc[index_old];
								if (FABS(delta) > limit)
								{
									printf("Vc diff\n");
								}

							}
						}
					}

					for (int j = 0; j < grid_rows; j++)
					{
						for (int i = 0; i < grid_cols; i++)
						{
							int index_old = i + j * grid_cols;
							int index_new = i + j * grid_cols_padded;

							NUMERIC_TYPE h_tmp = h_grid[index_new] + SGC_BankFullHeight_grid[index_new];
							NUMERIC_TYPE delta = h_tmp - Arrptr->H[index_old];
							if (FABS(delta) > limit)
							{
								step_error++;
								printf("error h (%d,%d) Solverptr->itCount == %ld \n", i, j, Solverptr->itCount);
								printf("h new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", h_tmp, Arrptr->H[index_old], delta);

								NUMERIC_TYPE tdelta = curr_time - Solverptr2.t;
								printf("t delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", curr_time, Solverptr2.t, tdelta);

								NUMERIC_TYPE hmaxdelta = reduce_Hmax - Hmax2;
								printf("hmax delta %" NUM_FMT" %" NUM_FMT" (%.4e)\n", reduce_Hmax, Hmax2, hmaxdelta);

								delta = volume_grid[index_new] - Arrptr->SGCVol[index_old];
								printf("vol new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", volume_grid[index_new], Arrptr->SGCVol[index_old], delta);

								//printf("vol new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", volume_grid[index_new], Arrptr->SGCbfH[index_old], delta);

								printf("H bound %d => %d\n", wet_dry_bounds->fp_h[j].start, wet_dry_bounds->fp_h[j].end);
								printf("H prev  %d => %d\n", wet_dry_bounds->fp_h_prev[j].start, wet_dry_bounds->fp_h_prev[j].end);
								printf("vol     %d => %d\n", wet_dry_bounds->fp_vol[j].start, wet_dry_bounds->fp_vol[j].end);


								int index_old_q;
								int index_new_q;
								NUMERIC_TYPE delta;
								NUMERIC_TYPE old_delta;
								NUMERIC_TYPE new_non_fp;
								NUMERIC_TYPE old_non_fp;
								NUMERIC_TYPE non_fp_delta;
								NUMERIC_TYPE q_old_m3;

								NUMERIC_TYPE q_total_new = C(0.0);
								NUMERIC_TYPE q_total_old = C(0.0);

								printf("-- Qx left\n");
								index_old_q = i + j * (grid_cols + 1);
								index_new_q = i + j * grid_cols_padded;
								delta = Qx_grid[index_new_q] - Arrptr->Qx[index_old_q];
								q_old_m3 = (Arrptr->Qxold[index_old_q] * dy_col[j]);
								printf("Qx    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_grid[index_new_q], Arrptr->Qx[index_old_q], delta);
								old_delta = Qx_old_grid[index_new_q] - q_old_m3;
								printf("Qxold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_old_grid[index_new_q], q_old_m3, old_delta);
								new_non_fp = Qx_grid[index_new_q] - Qx_old_grid[index_new_q];
								old_non_fp = Arrptr->Qx[index_old_q] - q_old_m3;
								non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								q_total_new += Qx_grid[index_new_q];
								q_total_old += Arrptr->Qx[index_old_q];


								printf("-- Qx right\n");
								index_old_q = i + j * (grid_cols + 1) + 1;
								index_new_q = i + j * grid_cols_padded + 1;
								delta = Qx_grid[index_new_q] - Arrptr->Qx[index_old_q];
								q_old_m3 = (Arrptr->Qxold[index_old_q] * dy_col[j]);
								printf("Qx    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_grid[index_new_q], Arrptr->Qx[index_old_q], delta);
								old_delta = Qx_old_grid[index_new_q] - q_old_m3;
								printf("Qxold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qx_old_grid[index_new_q], q_old_m3, old_delta);

								new_non_fp = Qx_grid[index_new_q] - Qx_old_grid[index_new_q];
								old_non_fp = Arrptr->Qx[index_old_q] - q_old_m3;
								non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								q_total_new -= Qx_grid[index_new_q];
								q_total_old -= Arrptr->Qx[index_old_q];


								printf("-- Qy up\n");
								index_old_q = i + j * (grid_cols + 1);
								index_new_q = i + j * grid_cols_padded;
								delta = Qy_grid[index_new_q] - Arrptr->Qy[index_old_q];
								q_old_m3 = (Arrptr->Qyold[index_old_q] * dx_col[j]);
								printf("Qy    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_grid[index_new_q], Arrptr->Qy[index_old_q], delta);
								old_delta = Qy_old_grid[index_new_q] - q_old_m3;
								printf("Qyold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_old_grid[index_new_q], q_old_m3, old_delta);

								new_non_fp = Qy_grid[index_new_q] - Qy_old_grid[index_new_q];
								old_non_fp = Arrptr->Qy[index_old_q] - q_old_m3;
								non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								q_total_new += Qy_grid[index_new_q];
								q_total_old += Arrptr->Qy[index_old_q];

								printf("-- Qy down\n");
								index_old_q = i + (j + 1) * (grid_cols + 1);
								index_new_q = i + (j + 1) * grid_cols_padded;
								delta = Qy_grid[index_new_q] - Arrptr->Qy[index_old_q];
								q_old_m3 = (Arrptr->Qyold[index_old_q] * dx_col[j]);
								printf("Qy    new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_grid[index_new_q], Arrptr->Qy[index_old_q], delta);
								old_delta = Qy_old_grid[index_new_q] - q_old_m3;
								printf("Qyold new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", Qy_old_grid[index_new_q], q_old_m3, old_delta);
								new_non_fp = Qy_grid[index_new_q] - Qy_old_grid[index_new_q];
								old_non_fp = Arrptr->Qy[index_old_q] - q_old_m3;
								non_fp_delta = new_non_fp - old_non_fp;
								printf("nonfp new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", new_non_fp, old_non_fp, non_fp_delta);
								q_total_new -= Qy_grid[index_new_q];
								q_total_old -= Arrptr->Qy[index_old_q];

								NUMERIC_TYPE dv_new = delta_time * q_total_new;
								NUMERIC_TYPE dv_old = delta_time * q_total_old;

								printf("deltaV from Q new: %" NUM_FMT" original:%" NUM_FMT" (%.4e) \n", dv_new, dv_old, dv_new - dv_old);

								printf("--\n");
							}
						}
					}
					if (step_error > 0)
						error_tripped++;
					if (error_tripped > 10)
					{
						exit(0);
					}
				}

				Solverptr2.SGCtmpTstep = getmin(Solverptr2.cfl*Parptr2.min_dx_dy*Parptr2.SGC_m / (SQRT(g * Hmax2)), Solverptr2.InitTstep);
				// Update t with final Tstep calculated
				if (Solverptr2.t > C(0.0)) Solverptr2.MinTstep = getmin(Solverptr2.MinTstep, Solverptr2.Tstep);
				Solverptr2.t += Solverptr2.Tstep;

#endif
				
				Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->min_dx_dy*Parptr->SGC_m / (SQRT(g * getmax(reduce_Hmax,Damptr->DamMaxH))), Solverptr->InitTstep);
				// Update t with final Tstep calculated
				if (curr_time > C(0.0))
					Solverptr->MinTstep = getmin(Solverptr->MinTstep, delta_time);
				curr_time += delta_time;
				Solverptr->t = curr_time;

				Solverptr->itCount++;
				itCount = Solverptr->itCount;

				// Calculate mass balance error (very rare)
				if (curr_time >= Parptr->MassTotal)
				{
					// calc losses for this mass interval
					loss = (Parptr->InfilTotalLoss - Parptr->InfilLoss) + (Parptr->EvapTotalLoss - Parptr->EvapLoss) - (Parptr->RainTotalLoss - Parptr->RainLoss);

					//Solverptr->Qerror=boundary_cond->Qin-boundary_cond->Qout-(Solverptr->vol2+loss-Solverptr->vol1)/Parptr->MassInt;
					// New version using VolInMT and VolOutMT
					// volume error
					
					Solverptr->Verror = boundary_cond->VolInMT - boundary_cond->VolOutMT - (Solverptr->vol2 + loss - Solverptr->vol1) + Damptr->DamLoss;
					
					// Q error
					Solverptr->Qerror = Solverptr->Verror / Parptr->MassInt;
					// reset to 0.0
					boundary_cond->VolInMT = C(0.0);
					boundary_cond->VolOutMT = C(0.0);
					Damptr->DamLoss = C(0.0);

					// record cumulative loss for next time.
					Parptr->InfilLoss = Parptr->InfilTotalLoss;
					Parptr->EvapLoss = Parptr->EvapTotalLoss;
					Parptr->RainLoss = Parptr->RainTotalLoss;


#ifdef _DEBUGxx
					fprintf(Fptr->mass_fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4" NUM_FMT" %12.4" NUM_FMT"  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT"\n",
						curr_time,
						delta_time,
						Solverptr->MinTstep,
						Solverptr->itCount,
						fix_small_negative(Solverptr->FArea),
						fix_small_negative(Solverptr->vol2),
						fix_small_negative(boundary_cond->Qin),
						fix_small_negative(Solverptr->Hds),
						fix_small_negative(boundary_cond->Qout),
						fix_small_negative(Solverptr->Qerror),
						fix_small_negative(Solverptr->Verror),
						fix_small_negative(Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss)));
#else
					fprintf(Fptr->mass_fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10li %12.4e %12.4e  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4e %12.4e %12.4e\n", curr_time, delta_time, Solverptr->MinTstep, Solverptr->itCount, Solverptr->FArea, Solverptr->vol2, boundary_cond->Qin, Solverptr->Hds, boundary_cond->Qout, Solverptr->Qerror, Solverptr->Verror, Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss));
#endif

					fflush(Fptr->mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
					if (Statesptr->DamMode == ON)
					{
						// Dam Output
						fprintf(Fptr->dam_fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT"", curr_time, delta_time);
						for (int i = 0; i < Damptr->NumDams; i++)
						{
							fprintf(Fptr->dam_fp, "%-10.3" NUM_FMT" %-10.3" NUM_FMT" %-10.3" NUM_FMT" %-10.3" NUM_FMT" %-10.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT"\n", Damptr->DamArea[i], Damptr->DamVol[i], Damptr->DamVin[i], Damptr->InitialHeight[i], Damptr->DamTotalQ[i] * delta_time, Damptr->SpillQ[i], Damptr->DamOperationQ[i], C(0.0));
						}
						fflush(Fptr->dam_fp);
					}
					Solverptr->vol1 = Solverptr->vol2;
					Parptr->MassTotal += Parptr->MassInt;

					//stage output
					if (Statesptr->save_stages == ON)
					{
						fprintf(Fptr->stage_fp, "%12.3" NUM_FMT"", curr_time);
						for (int i = 0; i < Locptr->Nstages; i++)
						{
							if (Locptr->stage_check[i] == 1)
							{
								int grid_index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * grid_cols_padded;
								fprintf(Fptr->stage_fp, "%10.4" NUM_FMT"", h_grid[grid_index] + SGC_BankFullHeight_grid[grid_index]);
							}
							else
								fprintf(Fptr->stage_fp, "-\t");
						}
						fprintf(Fptr->stage_fp, "\n");
						fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
						// added to export scalar velocity
						// velocity moved to do_update - between update q and update h
					}

					//virtual gauge output
					if (Statesptr->gsection == ON)
					{
						fprintf(Fptr->gau_fp, "%12.2" NUM_FMT"", curr_time); // print tiem to file
						for (int i = 0; i < Locptr->Ngauges; i++) // loop through each virtual gauge
						{
							// call discharge calculation function
							discharge = CalcVirtualGauge(i, grid_cols_padded, Qx_grid, Qy_grid, Locptr);
							fprintf(Fptr->gau_fp, " %10.3" NUM_FMT"", discharge); // Print discharge to file
						}
						fprintf(Fptr->gau_fp, "\n"); // print end of line
						fflush(Fptr->gau_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
					}

					// Checkpointing
					//if (Statesptr->checkpoint == ON)
					//{
					//	//iteration time
					//	time(&Solverptr->time_check);
					//	Solverptr->itrn_time_now = Solverptr->itrn_time + difftime(Solverptr->time_check, Solverptr->time_start);
					//	if (Solverptr->itrn_time_now >= Parptr->nextcheck)
					//	{
					//		WriteCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, verbose);
					//		Parptr->nextcheck = Solverptr->itrn_time_now + (Parptr->checkfreq * 3600.0);
					//	}
					//}
				}
				// export maximum depth interval
				if (Statesptr->maxint == ON && curr_time >= Parptr->maxintTotal)
				{
					// update h to be 'depth' above sub-grid-channel
					for (int j = 0; j < grid_rows; j++)
					{
						for (int i = 0; i < grid_cols; i++)
						{
							int index_padded = i + j * grid_cols_padded;
							int index = i + j * grid_cols;
							NUMERIC_TYPE temp = maxH_grid[index_padded] + SGC_BankFullHeight_grid[index_padded];
							// depth in channel or above flood plain = h + BankFullHeight (should not be negative)
							// note previous version just dumped the wd with no depth_thresh truncation
							if (temp <= depth_thresh)
								temp = C(0.0);
							tmp_grid1[index] = temp;
						}
					}
					write_grid(Fnameptr->resrootname, Parptr->maxintcount,
						NETCDF_IGNORE, ".wd_max",
						tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, &Statesptr->output_params);


					// update h to be 'depth' above sub-grid-channel
					for (int j = 0; j < grid_rows; j++)
					{
						for (int i = 0; i < grid_cols; i++)
						{
							int index_padded = i + j * grid_cols_padded;
							int index = i + j * grid_cols;
							maxH_grid[index_padded] = C(0.0); // rests max depth to zero at save interval
						}
					}

					// update interval counter
					Parptr->maxintTotal += Parptr->maxint;
					Parptr->maxintcount += 1;

				}


				// Regular output
				if (curr_time >= Parptr->SaveTotal)
				{
					gettimeofday(&timstr, NULL);
					double write_start_time = timstr.tv_sec + (timstr.tv_usec / 1000000.0);

					time(&Solverptr->time_check);
					Comp_time = (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start) / C(60.0);
					if (Comp_time != 0 && Statesptr->comp_out == ON) // only of t is not zero (can't divide by zero)
					{
						Model_Comp_Ratio = ((curr_time / 60) / Comp_time);
						Model_time_left = (Solverptr->Sim_Time - curr_time) / 60;
						Est_Time_Fin = (Model_time_left / Model_Comp_Ratio);
						Est_Time_Tot = Comp_time + Est_Time_Fin;
						printf("T(mins): M: %.1" NUM_FMT", C: %.1" NUM_FMT", M/C: %.2" NUM_FMT", ETot: %.1" NUM_FMT", EFin: %.1" NUM_FMT"\n", (curr_time / C(60.0)), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
					}
#if _PROFILE_MODE < 2			
					write_regular_output(Fnameptr->resrootname, 
						grid_cols, grid_rows, grid_cols_padded, 
						depth_thresh, curr_time,
						tmp_grid1, tmp_grid2, tmp_grid3,
						h_grid, dem_grid,
						Qx_grid, Qy_grid,
						Qx_old_grid, Qy_old_grid,
						Vx_grid, Vy_grid,
						SGC_BankFullHeight_grid, maxH_grid,
						sub_grid_layout_rows, sub_grid_state_rows,

						dx_col, dy_col,
						Solverptr, Statesptr, Parptr, boundary_cond, SGCptr,
						&Statesptr->output_params,
						Parptr->SaveNo,
						Statesptr->save_depth, Statesptr->save_elev, Statesptr->save_Qs, 
						Statesptr->voutput, Statesptr->SGCvoutput);
#endif
					// if regular output includes max reset max
					if (Statesptr->saveint_max == ON)
					{
						// update h to be 'depth' above sub-grid-channel
						for (int j = 0; j < grid_rows; j++)
						{
							for (int i = 0; i < grid_cols; i++)
							{
								int index_padded = i + j * grid_cols_padded;
								int index = i + j * grid_cols;
								maxH_grid[index_padded] = C(0.0); // rests max depth to zero at save interval
							}
						}
					}
					
					// update interval counter
					Parptr->SaveTotal += Parptr->SaveInt;
					Parptr->SaveNo += 1;

					gettimeofday(&timstr, NULL);
					double write_end_time = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
					total_write_time += (write_end_time - write_start_time);
				}

				// Single overpass
				if (Statesptr->single_op == ON && curr_time >= Parptr->op)
				{
					if (verbose == ON) printf("Writing overpass at %" NUM_FMT" seconds\n", curr_time);
					Statesptr->single_op = OFF;
					// raster depth output
					// update h to be 'depth' above sub-grid-channel
					// note previous version just dumped the wd with no depth_thresh truncation
					// updated version sets depth to zero if less than depth_thresh
					write_sum_grid(Fnameptr->resrootname,
						NETCDF_IGNORE,
						-1, ".op", h_grid, SGC_BankFullHeight_grid, tmp_grid1, depth_thresh, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, &Statesptr->output_params);
					// raster elevation output
					if (Statesptr->save_elev == ON)
					{
						write_sum_grid(Fnameptr->resrootname,
							NETCDF_IGNORE,
							-1, ".opelev", h_grid, dem_grid, tmp_grid1, depth_thresh, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, &Statesptr->output_params);
					}
				}

				// Multiple overpasses
				if (Statesptr->multi_op == ON)
				{
					for (int i = 0; i < Parptr->op_multinum; i++)
					{
						if (curr_time >= Parptr->op_multisteps[i] && Parptr->op_multiswitch[i] == 0)
						{
							Parptr->op_multiswitch[i] = 1;
							if (verbose == ON) printf("Writing overpass %d at %" NUM_FMT" seconds\n", i, Parptr->op_multisteps[i]);

							// raster depth output
							// update h to be 'depth' above sub-grid-channel
							// note previous version just dumped the wd with no depth_thresh truncation
							// updated version sets depth to zero if less than depth_thresh
							write_sum_grid(Fnameptr->resrootname,
								NETCDF_IGNORE,
								i, "-T.op", h_grid, SGC_BankFullHeight_grid, tmp_grid1, depth_thresh, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, &Statesptr->output_params);

							// raster elevation output
							if (Statesptr->save_elev == ON)
							{
								write_sum_grid(Fnameptr->resrootname,
									NETCDF_IGNORE,
									i, "-T.opelev", h_grid, dem_grid, tmp_grid1, depth_thresh, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, &Statesptr->output_params);
							}
						}
					}
				}

				// If requested on command line, check whether we should kill this simulation...
				if (Statesptr->killsim == ON)
				{
					//iteration time
					time(&Solverptr->time_check);
					Solverptr->itrn_time_now = Solverptr->itrn_time + (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start);
					// check if we have reached the kill time
					if (Solverptr->itrn_time_now >= Parptr->killsim_time)
					{
						printf("Simulation kill time reached... ");
						stop_loop = ON;
						//break;
					}
				}

				// If an auto steady-state has been requested:
				if (Statesptr->steadycheck == ON && curr_time >= Parptr->steadyTotal)
				{
					Parptr->steadyQdiff = boundary_cond->Qin - boundary_cond->Qout;
					if (FABS(Parptr->steadyQdiff) < Parptr->steadyQtol)
					{
						steadyCount += 1; // keep going a bit to make sure...
						if (steadyCount == 10)
						{
							printf("Simulation stop on steady state... ");
							stop_loop = ON; //break;

							const int tmp_save_num = 9999;
							write_regular_output(Fnameptr->resrootname,
								grid_cols, grid_rows, grid_cols_padded,
								depth_thresh, curr_time,
								tmp_grid1, tmp_grid2, tmp_grid3,
								h_grid, dem_grid,
								Qx_grid, Qy_grid,
								Qx_old_grid, Qy_old_grid,
								Vx_grid, Vy_grid,
								SGC_BankFullHeight_grid, maxH_grid,
								sub_grid_layout_rows, sub_grid_state_rows,

								dx_col, dy_col,
								Solverptr, Statesptr, Parptr, boundary_cond, SGCptr,
								& Statesptr->output_params,
								tmp_save_num, // set the save number to -1 for the final regular outputs
								Statesptr->save_depth, Statesptr->save_elev, Statesptr->save_Qs,
								Statesptr->voutput, Statesptr->SGCvoutput);
						}
					}
					else
					{
						steadyCount = 0;
					}
					Parptr->steadyTotal += Parptr->steadyInt;
				}
			} //end omp single section
		} // END main while loop ITERATIONS
	} // end of omp parallel


	gettimeofday(&timstr, NULL);
	processing_end_time = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
	double total_processing_time = (processing_end_time - processing_start_time) - total_write_time;

	printf("Elapsed processing time:\t\t\t%.6" NUM_FMT" (s) %ld\n", total_processing_time, Solverptr->itCount);
	printf("Output write time:\t\t\t%.6" NUM_FMT" (s) %d\n", total_write_time, Parptr->SaveNo);

	time_t loop_end;
	time(&loop_end);

	double seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);

#if defined (__INTEL_COMPILER) && _PROFILE_MODE > 0
	__itt_pause();
#endif
#if _PROFILE_MODE > 1
	//force stop 
	exit(0);
#endif


	//output final output
	WriteOutput(Fnameptr, grid_cols, grid_rows, grid_cols_padded,
		depth_thresh, 
		tmp_grid1,
		initHtm_grid, totalHtm_grid, maxH_grid, maxHtm_grid,
		maxVc_grid, maxVc_height_grid, maxHazard_grid,
		Vx_max_grid, Vy_max_grid,
		dem_grid, SGC_BankFullHeight_grid,
		Statesptr, Parptr, &Statesptr->output_params);
}

// version of this function without inline
// used by lis2_output.cpp
void SGC2_CalcA_public(int gr, NUMERIC_TYPE hflow, NUMERIC_TYPE bf, NUMERIC_TYPE *A, NUMERIC_TYPE *we, const SGCprams *SGCptr)
{
	SGC2_CalcA(gr, hflow, bf, A, we, SGCptr);
}

// version of this function without inline
// used by lis2_output.cpp
NUMERIC_TYPE SGC2_CalculateVelocity_public(const int index, const int index_next,
	const NUMERIC_TYPE * Q_grid,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE width)
{
	return SGC2_CalculateVelocity(index, index_next, Q_grid, h_grid, dem_grid, width);
}
void SGC2_CalcLinksQ(SuperGridLinksList * Super_linksptr, NUMERIC_TYPE * volume_grid, const NUMERIC_TYPE * h_grid, 
	const NUMERIC_TYPE delta_time, const NUMERIC_TYPE g, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE max_Froude, WetDryRowBound * wet_dry_bounds)
{
	int i;
	// Loop through links
	for (i = 0; i < Super_linksptr->num_links; i++)
	{
		NUMERIC_TYPE h0 = h_grid[Super_linksptr->link_index_SGC[i]] + Super_linksptr->SGC_bfH[i]; // sub grid channel depth
		NUMERIC_TYPE h1 = h_grid[Super_linksptr->link_index_2D[i]];
		NUMERIC_TYPE z0 = Super_linksptr->SGC_z[i];
		NUMERIC_TYPE z1 = Super_linksptr->DEM_z[i];

		NUMERIC_TYPE surface_elevation0 = z0 + h0;
		NUMERIC_TYPE surface_elevation1 = z1 + h1;
		// Calculating hflow based on floodplain levels
		NUMERIC_TYPE hflow = getmax(surface_elevation0, surface_elevation1) - getmax(z0, z1);
		NUMERIC_TYPE q_tmp, surface_slope;
		if (hflow > depth_thresh)
		{
			NUMERIC_TYPE area = (Super_linksptr->dx[i])* hflow;
			NUMERIC_TYPE dh = (surface_elevation0)-(surface_elevation1);
			surface_slope = -dh / Super_linksptr->dx[i];
			q_tmp = CalculateQ(surface_slope, hflow, delta_time, g, area, Super_linksptr->gn2[i], Super_linksptr->Qold[i], max_Froude);
			NUMERIC_TYPE dv = q_tmp*delta_time;
			volume_grid[Super_linksptr->link_index_SGC[i]] -= dv;
			volume_grid[Super_linksptr->link_index_2D[i]] += dv;
		}
		else
		{
			surface_slope = C(0.0);
			q_tmp = C(0.0);
		}
		Super_linksptr->Qold[i] = q_tmp;
		if (q_tmp != C(0.0))
		{
			int x = Super_linksptr->link_index_2D_i[i];
			int y = Super_linksptr->link_index_2D_j[i];
			wet_dry_bounds->fp_vol[y].start = min(wet_dry_bounds->fp_vol[y].start, x);
			wet_dry_bounds->fp_vol[y].end = max(wet_dry_bounds->fp_vol[y].end, x + 1);
		}
	}

}
