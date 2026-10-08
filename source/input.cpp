/*
*****************************************************************************
FILE INPUT
---------------------

A number of functions used to control the file input to LISFLOOD-FP. A short
definition of each is outlined below.

*****************************************************************************
*/
#include "lisflood.h"
#include "utility.h"
#include "sgc.h"
#include "lisflood2/file_tool.h"
#include "swe/fields.h"
#include "swe/dg2/fields.h"

FILE* fopen_or_die(const char * filename, const char* mode, const char* message, const int verbose)
{
	FILE *fp;
	fp = fopen(filename, mode);

	if (fp == NULL)
	{
		fprintf(stderr, "ERROR: %s. Aborting.\t%s\n", message, filename);
		printf("ERROR: %s. Aborting.\t%s\n", message, filename);
		exit(-1);
	}

	if (verbose == ON) printf("%s", message);

	return fp;
}

/// used to compare 2 floating point numbers
/// using == isn't reliable as float values may vary
inline int AreEqual(NUMERIC_TYPE a, NUMERIC_TYPE b)
{
	return FABS(a - b) < C(1e-5);
}

//-----------------------------------------------------------------------------
// LOAD STAGE DATA
void LoadStages(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Stage *Locptr, const int verbose)
{
	//Added by Matt Wilson, 5 Apr 2004
	//Provides functionality to output regular point measurements of water stage

	int i;
	FILE *fp;

	fp = fopen_or_die(Fnameptr->stagefilename, "r", "Loading stage information ", verbose);

	fscanf(fp, "%d", &Locptr->Nstages);
	fgetc(fp); // Retrieve closing EOL

	Locptr->stage_loc_x = memory_allocate_numeric_legacy(Locptr->Nstages);
	Locptr->stage_loc_y = memory_allocate_numeric_legacy(Locptr->Nstages);
	Locptr->stage_grid_x = new int[Locptr->Nstages];
	Locptr->stage_grid_y = new int[Locptr->Nstages];
	Locptr->stage_check = new int[Locptr->Nstages];

	//scan x,y locations from file
	for (i = 0; i < Locptr->Nstages; i++)
	{
		fscanf(fp, "%" NUM_FMT"", &Locptr->stage_loc_x[i]);
		fscanf(fp, "%" NUM_FMT"", &Locptr->stage_loc_y[i]);
	}
	//convert coordinates to grid cell numbers
	for (i = 0; i < Locptr->Nstages; i++)
	{
		Locptr->stage_grid_x[i] = int(floor((Locptr->stage_loc_x[i] - Parptr->blx) / Parptr->dx));
		Locptr->stage_grid_y[i] = Parptr->ysz - 1 - (int(floor((Locptr->stage_loc_y[i] - Parptr->bly) / Parptr->dy)));
		Locptr->stage_check[i] = 1;
	}
	//check for off-image values
	for (i = 0; i < Locptr->Nstages; i++)
	{
		if (Locptr->stage_grid_x[i] < 0 || Locptr->stage_grid_x[i] >= Parptr->xsz || Locptr->stage_grid_y[i] < 0 || Locptr->stage_grid_y[i] >= Parptr->ysz)
		{
			Locptr->stage_check[i] = 0;
			if (verbose == ON) printf("WARNING: Stage off-image: %d\n", i + 1);
		}
	}
	
	if (verbose == ON) printf("Done.\n\n");

	fclose(fp);
	return;
}
//-----------------------------------------------------------------------------
// LOAD WEIR DATA
void LoadWeir(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	if (strlen(Fnameptr->weirfilename) == 0)
		return;
	//One-directional flow functionality added by Matt Wilson 13 Feb 2504

	FILE *fp;
	int nw, i, xi, yi, p0, p1;
	NUMERIC_TYPE x, y, z0, z1;
	char char_tmp[10], buff[LINE_BUFFER_LEN];
	char tag_w1[] = "W";
	char tag_w2[] = "w";
	char tag_e1[] = "E";
	char tag_e2[] = "e";
	char tag_s1[] = "S";
	char tag_s2[] = "s";
	char tag_n1[] = "N";
	char tag_n2[] = "n";
	//tags for one-directional flow (culverts)
	char tag_wf1[] = "WF";
	char tag_wf2[] = "wf";
	char tag_ef1[] = "EF";
	char tag_ef2[] = "ef";
	char tag_sf1[] = "SF";
	char tag_sf2[] = "sf";
	char tag_nf1[] = "NF";
	char tag_nf2[] = "nf";
	//tags for bridge/culvert
	char tag_wb1[] = "WB";
	char tag_wb2[] = "wb";
	char tag_eb1[] = "EB";
	char tag_eb2[] = "eb";
	char tag_sb1[] = "SB";
	char tag_sb2[] = "sb";
	char tag_nb1[] = "NB";
	char tag_nb2[] = "nb";

	fp = fopen_or_die(Fnameptr->weirfilename, "r", "Loading weir information", verbose);

	Statesptr->weirs = ON;

	fgets(buff, LINE_BUFFER_LEN, fp);
	sscanf(buff, "%i", &nw);

	Arrptr->weir_count = nw;

	Arrptr->Weir_hc = memory_allocate_numeric_legacy(nw);
	Arrptr->Weir_Cd = memory_allocate_numeric_legacy(nw);
	Arrptr->Weir_m = memory_allocate_numeric_legacy(nw);
	Arrptr->Weir_w = memory_allocate_numeric_legacy(nw);
	Arrptr->Weir_Typ = new EWeirType[nw]; // type of structure... weir = 0, bridge = 1;

	Arrptr->Weir_Fixdir = new EDirection[nw];   // Fixed flow directions
	Arrptr->Weir_Identx = new int[(Parptr->xsz + 1)*(Parptr->ysz + 1)];
	Arrptr->Weir_Identy = new int[(Parptr->xsz + 1)*(Parptr->ysz + 1)];

	// Defalut to -1 for no weir link, and 0 for fixed flow direction
	SetArrayValue(Arrptr->Weir_Identx, -1, (Parptr->xsz + 1) * (Parptr->ysz + 1));
	SetArrayValue(Arrptr->Weir_Identy, -1, (Parptr->xsz + 1) * (Parptr->ysz + 1));

	for (i = 0; i < nw; i++)
	{
		if (!fgets(buff, LINE_BUFFER_LEN, fp))
		{
			printf("Weir file expects more lines\n");
			exit(0);
		}

		if (sscanf(buff, "%" NUM_FMT" %" NUM_FMT" %s %" NUM_FMT" %" NUM_FMT" %" NUM_FMT" %" NUM_FMT"",
			&x, &y, char_tmp, Arrptr->Weir_Cd + i, Arrptr->Weir_hc + i, Arrptr->Weir_m + i, Arrptr->Weir_w + i) != 7)
			Arrptr->Weir_w[i] = Parptr->dx;

		xi = (int)((x - Parptr->blx) / Parptr->dx);
		yi = (int)((Parptr->tly - y) / Parptr->dy);
		if (xi < 0 || xi >= Parptr->xsz ||
			yi < 0 || yi >= Parptr->ysz)
		{
			printf("Invalid weir coordinate: %d (%d,%d)\n", i, xi, yi);
			exit(0);
		}

		// Unfixed flow direction = C(0.)
		if (strcmp(char_tmp, tag_w1) == 0 || strcmp(char_tmp, tag_w2) == 0)
		{
			Arrptr->Weir_Identx[xi + 1 + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_e1) == 0 || strcmp(char_tmp, tag_e2) == 0)
		{
			Arrptr->Weir_Identx[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_s1) == 0 || strcmp(char_tmp, tag_s2) == 0)
		{
			Arrptr->Weir_Identy[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_n1) == 0 || strcmp(char_tmp, tag_n2) == 0)
		{
			Arrptr->Weir_Identy[xi + (yi + 1)*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		// Control tags for one-directional flow (culverts)
		// Fixed flow directions: N = 1, E = 2, S = 3, W = 4.
		if (strcmp(char_tmp, tag_wf1) == 0 || strcmp(char_tmp, tag_wf2) == 0)
		{
			Arrptr->Weir_Identx[xi + 1 + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = West;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_ef1) == 0 || strcmp(char_tmp, tag_ef2) == 0)
		{
			Arrptr->Weir_Identx[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = East;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_sf1) == 0 || strcmp(char_tmp, tag_sf2) == 0)
		{
			Arrptr->Weir_Identy[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = South;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		if (strcmp(char_tmp, tag_nf1) == 0 || strcmp(char_tmp, tag_nf2) == 0)
		{
			Arrptr->Weir_Identy[xi + (yi + 1)*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = North;
			Arrptr->Weir_Typ[i] = EWeir_Weir;
		}
		// control tags for bridge
		if (strcmp(char_tmp, tag_wb1) == 0 || strcmp(char_tmp, tag_wb2) == 0)
		{
			Arrptr->Weir_Identx[xi + 1 + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Bridge;
		}
		if (strcmp(char_tmp, tag_eb1) == 0 || strcmp(char_tmp, tag_eb2) == 0)
		{
			Arrptr->Weir_Identx[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Bridge;
		}
		if (strcmp(char_tmp, tag_sb1) == 0 || strcmp(char_tmp, tag_sb2) == 0)
		{
			Arrptr->Weir_Identy[xi + yi*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Bridge;
		}
		if (strcmp(char_tmp, tag_nb1) == 0 || strcmp(char_tmp, tag_nb2) == 0)
		{
			Arrptr->Weir_Identy[xi + (yi + 1)*(Parptr->xsz + 1)] = i;
			Arrptr->Weir_Fixdir[i] = DirectionNA;
			Arrptr->Weir_Typ[i] = EWeir_Bridge;
		}

		// now a check to make sure that Arrptr->Weir_hc is greater than the ground elevation,
		// this is especally important for the SGC model where bed elevations may change.

		// first index the cell of the weir and get z0
		p0 = xi + yi*Parptr->xsz;
		if (Statesptr->SGC == ON && Arrptr->SGCwidth[p0] > C(0.0)) z0 = Arrptr->SGCz[p0];
		else z0 = Arrptr->DEM[p0];
		// Then index the cell it flows from and get the elevation
		if (strcmp(char_tmp, tag_w1) == 0 || strcmp(char_tmp, tag_w2) == 0 || strcmp(char_tmp, tag_wf1) == 0 || strcmp(char_tmp, tag_wf2) == 0 || strcmp(char_tmp, tag_wb1) == 0 || strcmp(char_tmp, tag_wb2) == 0)
		{
			p1 = xi + 1 + yi*Parptr->xsz;
		}
		if (strcmp(char_tmp, tag_e1) == 0 || strcmp(char_tmp, tag_e2) == 0 || strcmp(char_tmp, tag_ef1) == 0 || strcmp(char_tmp, tag_ef2) == 0 || strcmp(char_tmp, tag_eb1) == 0 || strcmp(char_tmp, tag_eb2) == 0)
		{
			p1 = xi - 1 + yi*Parptr->xsz;
		}
		if (strcmp(char_tmp, tag_s1) == 0 || strcmp(char_tmp, tag_s2) == 0 || strcmp(char_tmp, tag_sf1) == 0 || strcmp(char_tmp, tag_sf2) == 0 || strcmp(char_tmp, tag_sb1) == 0 || strcmp(char_tmp, tag_sb2) == 0)
		{
			p1 = xi + (yi - 1)*Parptr->xsz;
		}
		if (strcmp(char_tmp, tag_n1) == 0 || strcmp(char_tmp, tag_n2) == 0 || strcmp(char_tmp, tag_nf1) == 0 || strcmp(char_tmp, tag_nf2) == 0 || strcmp(char_tmp, tag_nb1) == 0 || strcmp(char_tmp, tag_nb2) == 0)
		{
			p1 = xi + (yi + 1)*Parptr->xsz;
		}
		if (Statesptr->SGC == ON && Arrptr->SGCwidth[p1] > C(0.0)) z1 = Arrptr->SGCz[p1];
		else z1 = Arrptr->DEM[p1];

		// now work out if either of the elevations (z0,z1) are above the weir crest hight.
		if (Arrptr->Weir_hc[i] < z0 || Arrptr->Weir_hc[i] < z1)
		{
			if (verbose == ON)
			{
				if (Arrptr->Weir_Typ[i] == EWeir_Weir)
				{
					printf("WARNING: Weir crest height is below DEM\n");
					// for sub-grid model increase the crest height
					//if(Statesptr->SGC==ON)
					//{
					Arrptr->Weir_hc[i] = getmax(z0, z1);
					if (verbose == ON) printf("Weir number %i crest height increased to %.3" NUM_FMT" m\n", i, Arrptr->Weir_hc[i]);
					//}
				}
				//else 
				//{
				//printf("WARNING: Bridge soffit height is below DEM converted to weir!!\n");
				//Arrptr->Weir_Typ[i] = 0;
				//}
			}
		}
		// need check for bridge less than or equal to subgrid width ?.......
	}

	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
//-----------------------------------------------------------------------------
/* LOAD RIVER CHANNEL NETWORK - CCS
Added by Chris Sampson January 2011
Loads one or more river systems into the ChannelSegments vector, using the RiversIndex vector to keep track of where each river is located
within ChannelSegents.  This is done by recording the size of ChannelSegments after each river is loaded as an int in RiversIndex. Example:

River Thames has 5 segments; River Severn has 4 segments.

After first LoadRiver loop:
ChannelSegments[0][1][2][3][4]					RiversIndex[0]
<....Thames...>							   <5>

After second LoadRiver loop:
ChannelSegments[0][1][2][3][4][5][6][7][8]		RiversIndex[0][1]
<....Thames...><..Severn..>				   <5><9>

Note: remember that the size recorded in RiversIndex is always 1 greater than the max index of ChannelSegments because size starts at 1
whereas the index starts at 0!

Just as with the old system, we use a pointer called CSTypePtr throughout the rest of the model code to point at ChannelSegments[0].  We also have
a new pointer called RiversIndexPtr to point RiversIndex[0].  Therefore CSTypePtr+5 points at ChannelSegments[5] etc.

An extra trick with the vectors is that we can have a pointer to the vector itself (rather than a particular element of the vector as above).  These are
called ChannelSegmentsVecPtr and RiversIndexVecPtr, and are used primarily to construct the vectors using member functions such as push.back.  They are
also sometimes employed to determine the size loops need to be using the size() member function.  For example, if we have loaded 5 rivers, then:
RiversIndexVecPtr->size() will equal C(5.) If each of these rivers contained 5 channel segments then ChannelSegmentsVecPtr->size() would equal C(25.)

*/

void LoadRiverNetwork(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, Arrays *Arrptr, vector<QID7_Store> *QID7_Vec_Ptr, vector<int> *RiversIndexVecPtr, const int verbose)
{
	FILE *rfp;  // local file pointer
	int i, n;
	int tmp_size;

	if (Statesptr->multiplerivers == 0)
	{
		LoadRiver(Fnameptr, Statesptr, Parptr, ChannelSegmentsVecPtr, Arrptr, QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver once.
		tmp_size = ChannelSegmentsVecPtr->size();
		RiversIndexVecPtr->push_back(tmp_size); // CCS
	}
	else if (Statesptr->multiplerivers == 1)
	{
		rfp = fopen_or_die(Fnameptr->multiriverfilename, "r", "Loading Rivers", verbose);
		fscanf(rfp, "%i", &n);
		if (verbose == ON) printf("Loading %i Rivers\n\n", n);

		for (i = 0; i < n; i++)
		{
			fscanf(rfp, "%s", Fnameptr->rivername); //Scan the next .river filename from the .rivers file and asign it to Fnameptr->rivername.
			LoadRiver(Fnameptr, Statesptr, Parptr, ChannelSegmentsVecPtr, Arrptr, QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver in loop.
			tmp_size = ChannelSegmentsVecPtr->size();
			RiversIndexVecPtr->push_back(tmp_size); /* Builds the RiversIndex vector so we know where one river stops and the next starts within
			the ChannelSegments vector. */
		}

		if (verbose == ON) printf("%i Rivers Loaded Successfully.\n\n", n);
	}

	return;
}
//-----------------------------------------------------------------------------
// LOAD RIVER DATA FROM FILE
void LoadRiver(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, Arrays *Arrptr, vector<QID7_Store> *QID7_Vec_Ptr, vector<int> *RiversIndexVecPtr, const int verbose)
{
	if (strlen(Fnameptr->rivername) == 0)
		return;

	FILE *fp;  // local file pointer
	int npoints, *xpi, *ypi;
	int *trib; // temp xs array to record any trib connections
	NUMERIC_TYPE *xp, *yp, *wp, *np, *hp, *cp, *rp, *ap, total_length = C(0.0), *qp;
	char buff[800], buff2[800], buff3[800];
	int i, j, pi, pj, ni, nj, oldpi, oldpj, i1, i2;
	NUMERIC_TYPE tmp1, tmp2, tmp3;
	NUMERIC_TYPE grad; // temporary gradient calculation variable eventually stored in csp->Shalf
	char *Q_Name_tmp;
	char buffer[80];
	ESourceType *Q_Ident_tmp;
	int count, chseg, SegOutNo, tmp_int;

	// MSH: csp is a utility pointer, 
	fp = fopen_or_die(Fnameptr->rivername, "r", "Loading Rivers", verbose);

	Statesptr->ChannelPresent = ON;

	fscanf(fp, "%s", buffer);
	if (!STRCMPi(buffer, "Tribs"))
	{
		Statesptr->TribsPresent = ON;

		// MSH: Since we haven't allocated the memory yet, we can't assign the number of channel segments in
		// the first element of the ChannelSegments array - so read into a temp variable
		fscanf(fp, "%i", &tmp_int);
		if (verbose == ON) printf("%i Channel Segments\n", tmp_int);
	}
	else
	{
		rewind(fp);
		tmp_int = 1;
		if (verbose == ON) printf("%i Channel Segment\n", tmp_int);
	}

	for (chseg = 0; chseg < tmp_int; chseg++) // CCS Slight reorganisation of old LoadRiver function but fundamentally unchanged.
	{
		ChannelSegmentType tmp_chan; // CCS
		tmp_chan.Next_Segment = tmp_chan.Next_Segment_Loc = -1; // CCS from old code
		tmp_chan.N_Channel_Segments = tmp_int; // CCS
		ChannelSegmentType *csp = &tmp_chan; // CCS

		fscanf(fp, "%i", &npoints);
		if (verbose == ON) printf("%i points in channel segment %i\n", npoints, chseg);

		//setup local temporary arrays, Note the () at the end ensures all elements are initialised to zero
		xp = memory_allocate_zero_numeric_legacy(npoints);
		yp = memory_allocate_zero_numeric_legacy(npoints);
		wp = memory_allocate_zero_numeric_legacy(npoints);
		np = memory_allocate_zero_numeric_legacy(npoints);
		hp = memory_allocate_zero_numeric_legacy(npoints);
		cp = memory_allocate_zero_numeric_legacy(npoints);
		qp = memory_allocate_zero_numeric_legacy(npoints);
		rp = memory_allocate_zero_numeric_legacy(npoints); // chainage ratio between entered cross sections and cell chainage
		ap = memory_allocate_zero_numeric_legacy(npoints); // actual entered xs chainage 
		trib = new int[npoints]();
		Q_Name_tmp = new char[npoints * 80]();
		Q_Ident_tmp = new ESourceType[npoints]();

		for (i = 0; i < npoints; i++)
		{
			// set default qp and trib value (all the other arrays are pre zeroed with new command and end brackets () )
			qp[i] = -1;
			trib[i] = -1;

			fscanf(fp, "%" NUM_FMT" %" NUM_FMT"", xp + i, yp + i); // Load x,y values.

			// load buffer until EOL
			j = 0;
			do{ buff[j] = fgetc(fp); } while (buff[j++] != '\n');
			buff[j] = '\0';									// Finish off string
			if (sscanf(buff, "%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"", &tmp1, &tmp2, &tmp3) == 3)
			{        										// Only store values if 3 reads successful
				wp[i] = tmp1;
				np[i] = tmp2;
				hp[i] = tmp3;
				if (verbose == ON)
					printf("Xsec %4i\tw=%8.3" NUM_FMT" n=%5.3" NUM_FMT" z=%6.3" NUM_FMT"\n", i, wp[i], np[i], hp[i]);
			}

			if (sscanf(buff, "%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%s%s", &tmp1, &tmp1, &tmp1, buff2, buff3) == 5  // 4+5th item found must be Q_Ident 
				|| sscanf(buff, "%s%s", buff2, buff3) == 2) // OR No channel info - just Q in
			{
				if (!STRCMPi(buff2, "QFIX"))
				{
					Q_Ident_tmp[i] = QFIX4;
					sscanf(buff3, "%" NUM_FMT"", qp + i);
					if (verbose == ON) printf("Xsec %4i\tQFIX at %7.2" NUM_FMT"\n", i, qp[i]);
				}
				if (!STRCMPi(buff2, "QVAR"))
				{
					Q_Ident_tmp[i] = QVAR5;
					strcpy(Q_Name_tmp + i * 80, buff3);
					if (verbose == ON) printf("Xsec %4i\tQVAR from bdy file %s\n", i, Q_Name_tmp + i * 80);
				}
				if (!STRCMPi(buff2, "QOUT"))
				{
					Q_Ident_tmp[i] = FREE6;
					sscanf(buff3, "%i", &SegOutNo);
					if (verbose == ON) printf("Xsec %4i\tQOUT from segment %i discharges into segment %i\n", i, chseg, SegOutNo);
				}
				if (!STRCMPi(buff2, "TRIB"))
				{
					Q_Ident_tmp[i] = TRIB7;
					sscanf(buff3, "%i", &SegOutNo);
					if (verbose == ON) printf("Xsec %4i\tQin from segment %i\n", i, SegOutNo);
					trib[i] = SegOutNo;
				}
				if (!STRCMPi(buff2, "FREE"))
					// normal depth based on slope, if -1 then use last channel segment slope else use slope supplied
					// NOT fully working for Diffusive - stability issues
				{
					Q_Ident_tmp[i] = FREE1;
					sscanf(buff3, "%" NUM_FMT"", qp + i);
					if (qp[i] < C(-0.999)) // ie -1 (done like this as NUMERIC_TYPE)
					{
						if (verbose == ON) printf("Xsec %4i\tFREE using end slope\n", i);
					}
					else
					{
						if (verbose == ON) printf("Xsec %4i\tFREE using slope %7.4" NUM_FMT"\n", i, qp[i]);
					}
				}
				if (!STRCMPi(buff2, "HFIX"))
				{
					Q_Ident_tmp[i] = HFIX2;
					sscanf(buff3, "%" NUM_FMT"", qp + i);
					if (verbose == ON) printf("Xsec %4i\tHFIX at %7.2" NUM_FMT"\n", i, qp[i]);
				}
				if (!STRCMPi(buff2, "HVAR"))
				{
					Q_Ident_tmp[i] = HVAR3;
					strcpy(Q_Name_tmp + i * 80, buff3);
					if (verbose == ON) printf("Xsec %4i\tHVAR from bdy file %s\n", i, Q_Name_tmp + i * 80);
				}
				if (!STRCMPi(buff2, "RATE"))
					// NOT fully working for Diffusive - stability issues
				{
					Q_Ident_tmp[i] = RATE8;
					strcpy(Q_Name_tmp + i * 80, buff3);
					if (verbose == ON) printf("Xsec %4i\tRATE from bdy file %s\n", i, Q_Name_tmp + i * 80);
				}
			}
		}

		if (verbose == ON) printf("Channel data read for segment %i - interpolating values.\n", chseg);



		// Estimate number of channel pixels - to ensure we allocate enough memory for temporary arrays xpi,ypi
		// total length divided by cell size and then NUMERIC_TYPE this value.
		ap[0] = 0;
		for (i = 0; i < npoints - 1; i++)
		{
			// calc straight line chainage between entered cross sections - used for cell independent chainage calcs. Add up chainage
			ap[i + 1] = ap[i] + sqrt((xp[i + 1] - xp[i])*(xp[i + 1] - xp[i]) + (yp[i + 1] - yp[i])*(yp[i + 1] - yp[i]));
		}
		xpi = new int[int(C(2.0)*ap[npoints - 1] / Parptr->dx)]();
		ypi = new int[int(C(2.0)*ap[npoints - 1] / Parptr->dx)]();



		// Insert channel into DEM grid
		oldpi = (int)((xp[0] - Parptr->tlx) / Parptr->dx);
		oldpj = (int)((Parptr->tly - yp[0]) / Parptr->dy);
		total_length = C(0.0);

		count = 0;
		for (i = 1; i < npoints; i++)
		{
			pi = (int)((xp[i] - Parptr->tlx) / Parptr->dx);
			pj = (int)((Parptr->tly - yp[i]) / Parptr->dy);
			for (j = 0; j <= 1000; j++)					// Take very small steps and insert channel
			{															// whenever x,y position changes
				ni = oldpi + ((pi - oldpi)*j / 1000);
				nj = oldpj + ((pj - oldpj)*j / 1000);
				if (ni >= 0 && ni < Parptr->xsz && nj >= 0 && nj < Parptr->ysz) // check it stays within DEM
				{
					if (count == 0)								// Always insert first point
					{
						xpi[count] = ni;
						ypi[count] = nj;
						Arrptr->ChanMask[ni + nj*Parptr->xsz] = 1; // mark mask with value of 1 - will renumber later in order
						count++;
					}
					else if (ni != xpi[count - 1] || nj != ypi[count - 1]) // if grid location changes
					{
						if (Arrptr->ChanMask[ni + nj*Parptr->xsz] == -1)   // channel mask not set
						{
							xpi[count] = ni;
							ypi[count] = nj;
							Arrptr->ChanMask[ni + nj*Parptr->xsz] = 1; // mark mask with value of 1 - will renumber later in order
							total_length += Parptr->dx*sqrt(pow((NUMERIC_TYPE)(ni - xpi[count - 1]), C(2.0)) + pow((NUMERIC_TYPE)(nj - ypi[count - 1]), C(2.0)));
							count++;
						}
						else
						{
							// channel mask set, so likely that it is trib junction point. 
							// NOTE, cannot have crossing channels !!!
							// DO NOT mark mask with value of 1 - as this is the junction 
							// of the trib with main channel so is already marked for main channel
							xpi[count] = ni;
							ypi[count] = nj;
							total_length += Parptr->dx*sqrt(pow((NUMERIC_TYPE)(ni - xpi[count - 1]), 2) + pow((NUMERIC_TYPE)(nj - ypi[count - 1]), 2));
							count++;
						}
					}
				}
			}
			oldpi = pi;
			oldpj = pj;
			cp[i] = total_length;
			rp[i] = (cp[i] - cp[i - 1]) / (ap[i] - ap[i - 1]);
		}
		csp->chsz = count;

		if (count == 0) printf("\nWARNING: no overlap with DEM cells for channel %i.\n", chseg);

		// Set up other channel rasters and fill in values
		csp->ChanX = new int[csp->chsz]();
		csp->ChanY = new int[csp->chsz]();
		csp->Chandx = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->Shalf = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->A = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->NewA = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->Chainage = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->ChanQ = memory_allocate_zero_numeric_legacy(csp->chsz); // only used to record Q values for output in profile - not used in calc
		csp->ChanWidth = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->ChanN = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->Q_Val = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->BankZ = memory_allocate_zero_numeric_legacy(csp->chsz);
		csp->Q_Name = new char[csp->chsz * 80]();
		csp->Q_Ident = new ESourceType[csp->chsz]();

		for (i = 0; i < csp->chsz; i++)
		{
			csp->ChanX[i] = xpi[i];
			csp->ChanY[i] = ypi[i];
		}


		// Find chainage and dx along channel
		for (i = 0; i < csp->chsz - 1; i++)
			csp->Chandx[i] = sqrt(pow(Parptr->dx*(csp->ChanX[i] - csp->ChanX[i + 1]), 2) +
			pow(Parptr->dy*(csp->ChanY[i] - csp->ChanY[i + 1]), 2));

		// assume dx for last cell is same as last segment
		csp->Chandx[csp->chsz - 1] = sqrt(pow(Parptr->dx*(csp->ChanX[csp->chsz - 1] - csp->ChanX[csp->chsz - 2]), 2) +
			pow(Parptr->dy*(csp->ChanY[csp->chsz - 1] - csp->ChanY[csp->chsz - 2]), 2));
		csp->Chainage[0] = 0;

		// add up dx to get chainage
		for (i = 1; i < csp->chsz; i++)csp->Chainage[i] = csp->Chainage[i - 1] + csp->Chandx[i - 1];

		// Fill in channel mask
		for (i = 0; i < csp->chsz; i++)
		{
			pi = csp->ChanX[i];
			pj = csp->ChanY[i];

			// renumber channel mask
			Arrptr->ChanMask[pi + pj*Parptr->xsz] = i;
			// set bank level
			csp->BankZ[i] = Arrptr->DEM[pi + pj*Parptr->xsz];
			// mark segment mask
			Arrptr->SegMask[pi + pj*Parptr->xsz] = chseg;
		}

		// Adjust chainage calcs so that it is independent of cell size. ####
		if (Statesptr->chainagecalc == ON)
		{
			if (verbose == ON) printf("Cell size independent channel chainage calculations are ON.\n");
			// adjust each dx by ratio calculated previously
			for (i = 0; i < npoints; i++)
			{
				for (j = 0; j < csp->chsz; j++)
				{
					if (csp->Chainage[j] >= cp[i] && csp->Chainage[j] < (cp[i + 1]))
					{
						csp->Chandx[j] = csp->Chandx[j] / rp[i + 1]; // adjust by ratio calculated previously
					}
				}
			}
			// calc last seg dx as same as penultimate
			csp->Chandx[csp->chsz - 1] = csp->Chandx[csp->chsz - 2];

			// add up chainage again
			for (i = 1; i < csp->chsz; i++)csp->Chainage[i] = csp->Chainage[i - 1] + csp->Chandx[i - 1];

			// also adjust original xs chainage as this was based on centre cell distance
			for (i = 1; i < npoints; i++) cp[i] = ap[i];
		}
		else if (verbose == ON) printf("Cell size independent channel chainage calculations are OFF.\n");




		// Interpolate width, Mannings n
		i1 = 0; i2 = 0;
		while (i2 < npoints - 1)
		{
			for (i2 = i1 + 1; i2<npoints; i2++)
			{
				if (wp[i2]>C(0.0)) break; // loop until next non zero value found
			}
			for (i = 0; i < csp->chsz; i++)
			{
				if (csp->Chainage[i] >= cp[i1] && csp->Chainage[i] <= (cp[i2] + C(0.1))) // add C(0.1)m to cp[] to get last point
				{
					pi = csp->ChanX[i];
					pj = csp->ChanY[i];
					csp->ChanWidth[i] = wp[i1] + (wp[i2] - wp[i1])*(csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
					csp->ChanN[i] = np[i1] + (np[i2] - np[i1])*(csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
					if (i == csp->chsz - 1 && chseg != 0)
					{
						// special case where trib junction with main channel. Do not overwrite main channel bed elevation but 
						// otherwise record channel info for trib end point in dummy node (allows channel BC link for diffusive 
						// and correct slope calc for kinematic).
						csp->JunctionDEM = hp[i1] + (hp[i2] - hp[i1])*(csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
					}
					else
					{
						Arrptr->DEM[pi + pj*Parptr->xsz] = hp[i1] + (hp[i2] - hp[i1])*(csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
					}
				}
			}
			i1 = i2;
		}


		// Interpolate slope
		for (i = 0; i < csp->chsz; i++)
		{
			// find gradient of segment
			if (chseg == 0) // main channel
			{
				if (i == csp->chsz - 1) // last point
				{
					// special case for last point as we can only know the slope of the segment behind it
					grad = (Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i - 1] + csp->ChanY[i - 1] * Parptr->xsz])
						/ csp->Chandx[i - 1];
				}
				else // all other points
				{
					grad = (Arrptr->DEM[csp->ChanX[i + 1] + csp->ChanY[i + 1] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
						/ csp->Chandx[i];
				}
			}
			else // trib special case at end due to junction
			{
				if (i == csp->chsz - 2) // last but one point
				{
					// tribs use dummy node for last point - ie junction ##
					grad = (csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
						/ csp->Chandx[i];
				}
				else if (i == csp->chsz - 1) // last point
				{
					// tribs use dummy node for last point - ie junction ##
					grad = (csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i - 1] + csp->ChanY[i - 1] * Parptr->xsz])
						/ csp->Chandx[i - 1];
				}
				else // all other points as main channel
				{
					grad = (Arrptr->DEM[csp->ChanX[i + 1] + csp->ChanY[i + 1] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
						/ csp->Chandx[i];
				}
			}


			// check if slope is positive or negative
			if (grad >= 0)
			{
				if (Statesptr->diffusive == ON)
				{
					// only keep uphill slopes for diffusive
					csp->Shalf[i] = -1 * sqrt(fabs(grad));
				}
				else
				{
					// for kinematic we just pretend it is a downhill
					csp->Shalf[i] = sqrt(fabs(grad));
					// also warn user, in case they don't know
					if (verbose == ON) printf("\nWARNING: Kinematic solver BUT uphill slope at point %i for channel %i.\n", i - 1, chseg);
				}
			}
			else
			{
				csp->Shalf[i] = sqrt(fabs(grad));
			}
		}

		// Fill in Q boundary conditions from _tmp arrays
		for (i = 0; i < npoints; i++) // loop through the input cross-sections
		{
			if (Q_Ident_tmp[i] != 0)  // if there is a boundary condition
			{
				for (j = 0; j < csp->chsz; j++) // loop through the cross-sections mapped onto the dem space
				{
					if ((cp[i] + C(0.1)) >= csp->Chainage[j] && (cp[i] + C(0.1)) < (csp->Chainage[j] + csp->Chandx[j]))
						// make sure we only apply the BC to one point
					{
						csp->Q_Ident[j] = Q_Ident_tmp[i]; // copy type across
						csp->Q_Val[j] = qp[i];            // copy value across
						if (Q_Ident_tmp[i] == HVAR3 || Q_Ident_tmp[i] == QVAR5 || Q_Ident_tmp[i] == RATE8)	// check if BC type has name
						{
							strcpy(csp->Q_Name + j * 80, Q_Name_tmp + i * 80); // copy name across
						}
						if (Q_Ident_tmp[i] == TRIB7) // only for tribs
						{
							/*
							record in the trib data the location and segment it links to
							ChannelSegments[trib[i]].Next_Segment_Loc=j; // CCS
							ChannelSegments[trib[i]].Next_Segment=chseg; // CCS

							^^ The above code is left commented out to show why the following QID7_Store vector is needed.
							As these terms need to be written to an instance of ChannelSegmentType not pointed to by csp
							(ChannelSegments[trib[i]]), we store them in the QID7 vector and use the UpdateChannelVector function
							to move the contents to the correct place after LoadRiver has finished. // CCS
							*/

							QID7_Store QID7_tmp; // CCS create temp instance of QID7_Store and populate struct:
							QID7_tmp.chseg = chseg;
							QID7_tmp.Next_Segment_Loc = j;
							QID7_tmp.trib = trib[i];
							if (Statesptr->multiplerivers == ON)
							{
								QID7_tmp.RiverID = RiversIndexVecPtr->size(); /*#CCS# this allows us to keep track of which river we are in when later using QID7.
																			1st river ID will be 0, 2nd will be 1, etc.*/
							}
							QID7_Vec_Ptr->push_back(QID7_tmp); // CCS push_back temp instance into external vector
						}
					}
				}
			}
		}

		// release memory used for temporary variables
		delete[] xp;
		delete[] yp;
		delete[] wp;
		delete[] np;
		delete[] hp;
		delete[] xpi;
		delete[] ypi;
		delete[] cp;
		delete[] Q_Name_tmp;

		if (Statesptr->diffusive == 1 && csp->Q_Ident[csp->chsz - 1] == 0)
		{
			if (verbose == ON) printf("\nWARNING: Channel %i has no d/s BC using FREE with channel slope\n", chseg); // warn user that for diffusive no BC is set so using free
			csp->Q_Ident[csp->chsz - 1] = FREE1;  // copy type across
			csp->Q_Val[csp->chsz - 1] = -1;   // copy value across
		}

		ChannelSegmentsVecPtr->push_back(tmp_chan);

	} // end of channel segment loop


	fclose(fp);
	if (verbose == ON) printf("Done.\n\n");

	return;
}
//-----------------------------------------------------------------------------
/*UPDATE CHANNEL VECTOR FUNCTION // CCS
This function is needed to update the ChannelSegments vector after LoadRiver has run.  It simply copies the  terms held
in QID7_Store to the correct place in ChannelSegments.  It solves the problem of the LoadRiver function attempting to write
to an index of ChannelSegemnts that hasn't yet been created when tagging tributary junctions.*/

void UpdateChannelsVector(States *Statesptr, ChannelSegmentType *CSTypePtr, vector<QID7_Store> *QID7_Vec_Ptr, QID7_Store *QID7Ptr, int *RiversIndexPtr)
{
	int vecsize, i, n;
	if (Statesptr->ChannelPresent == OFF) return;
	vecsize = QID7_Vec_Ptr->size();
	for (i = 0; i < vecsize; i++)
	{
		if (Statesptr->multiplerivers == OFF)
		{
			n = QID7Ptr[i].trib;
			CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg;
			CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
		}
		else if (Statesptr->multiplerivers == ON)
		{
			if (QID7Ptr[i].RiverID == 0)
			{
				n = QID7Ptr[i].trib;
				CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg;
				CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
			}
			else
			{
				n = QID7Ptr[i].trib + RiversIndexPtr[QID7Ptr[i].RiverID - 1]; //Make sure we write to the correct CSTypePtr index when not in first river.
				CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg + RiversIndexPtr[QID7Ptr[i].RiverID - 1]; //Again, when not in first river we need to ensure correct ID of next segment.
				CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
			}
		}
	}

	return;

}
//-----------------------------------------------------------------------------
// LOAD DISTRIBUTED INFILTRATION RATES FROM FILE
void LoadDistInfil(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	if (strlen(Fnameptr->infilfilename) == 0)
		return;

	FILE *fp;
	int i, j;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->infilfilename, "r", "Loading distributed infiltration data", verbose);

	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	Arrptr->dist_infiltration = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%" NUM_FMT"", Arrptr->dist_infiltration + i + j * Parptr->xsz);
		if (AreEqual(Arrptr->dist_infiltration[i + j * Parptr->xsz], no_data_value))
			Arrptr->dist_infiltration[i + j * Parptr->xsz] = 0.0;
		// convert to m/s from mm per hour
		Arrptr->dist_infiltration[i + j * Parptr->xsz] = Arrptr->dist_infiltration[i + j * Parptr->xsz] / 60 / 60 / 1000;
	}
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");
	return;
}

//-----------------------------------------------------------------------------
// LOAD FLOODPLAIN FRICTION FROM FILE
void LoadManningsn(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	if (strlen(Fnameptr->nfilename) == 0)
		return;

	FILE *fp;
	int i, j;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->nfilename, "r", "Loading floodplain Manning's n data", verbose);

	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	Arrptr->Manningsn = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%" NUM_FMT"", Arrptr->Manningsn + i + j*Parptr->xsz);
		if (AreEqual(Arrptr->Manningsn[i + j*Parptr->xsz], no_data_value))
			Arrptr->Manningsn[i + j*Parptr->xsz] = Parptr->FPn;
	}
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
void LoadSGCManningsn(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	if (strlen(Fnameptr->SGCnfilename) == 0)
		return;

	FILE *fp;
	int i, j;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->SGCnfilename, "r", "Loading SGC Manning's n data", verbose);

	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	Arrptr->SGCManningsn = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%" NUM_FMT"", Arrptr->SGCManningsn + i + j*Parptr->xsz);
		if (AreEqual(Arrptr->SGCManningsn[i + j*Parptr->xsz], no_data_value))
			Arrptr->SGCManningsn[i + j*Parptr->xsz] = Parptr->SGC_n;
	}
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
// Load SGC direction from array PFU
void LoadSGCdirn(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	if (strlen(Fnameptr->SGCdirnfilename) == 0)
		return;

	FILE *fp;
	int i, j;
	char dum[800];
	int no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->SGCdirnfilename, "r", "Loading SGC Direction data", verbose);

	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %i", dum, &no_data_value);

	//Arrptr->SGCdirn = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->SGCdirn =  new int[Parptr->xsz*Parptr->ysz];
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%i", Arrptr->SGCdirn + i + j*Parptr->xsz);
		if (AreEqual(Arrptr->SGCdirn[i + j*Parptr->xsz], no_data_value))
			Arrptr->SGCdirn[i + j*Parptr->xsz] = 0;
	}
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
//-----------------------------------------------------------------------------
// LOAD POROSITY FROM FILE, TJF
void LoadPor(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	FILE *fp;
	int i, j, k, m;
	NUMERIC_TYPE incr;
	char dum[800], buff[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen(Fnameptr->porfilename, "r");
	if (fp == NULL)
	{
		if (verbose == ON)
		{
			Statesptr->porosity = OFF;
			printf("Porosity off\n");
		}
		return;
	}
	if (verbose == ON) printf("Loading porosity data:\t%s\n", Fnameptr->porfilename);

	Statesptr->porosity = ON;

	fscanf(fp, "%s", buff);

	// Determining the porosity method
	if (!STRCMPi(buff, "PFIX"))
	{
		Parptr->Por_Ident = 1;
		if (verbose == ON) printf("Fixed Aerial Porosity Method\n");
	}

	if (!STRCMPi(buff, "PVAR"))
	{
		Parptr->Por_Ident = 2;
		if (verbose == ON) printf("Water Height Dependent Aerial Porosity Method\n");
	}

	if (!STRCMPi(buff, "PBOUND"))
	{
		Parptr->Por_Ident = 3;
		if (verbose == ON) printf("Fixed Boundary Porosity Method\n");
	}

	if (!STRCMPi(buff, "PBVAR"))
	{
		Parptr->Por_Ident = 4;
		if (verbose == ON) printf("Water Height Dependent Boundary Porosity Method\n");
	}

	// Loading the porosity values
	// Por_Ident=1
	if (Parptr->Por_Ident == 1)
	{
		Arrptr->paerial = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->paerial + i + j*Parptr->xsz);
			if (AreEqual(Arrptr->paerial[i + j*Parptr->xsz], no_data_value))
				Arrptr->paerial[i + j*Parptr->xsz] = 1;
		}
		fclose(fp);
	}

	// Por_Ident=2
	else if (Parptr->Por_Ident == 2)
	{
		fscanf(fp, "%s %i", dum, &Parptr->zsz);
		fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->maxelev);
		fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->zlev);

		Arrptr->paerial = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz*Parptr->zsz);

		for (k = 0; k < Parptr->zsz; k++)
		{
			fscanf(fp, "%" NUM_FMT"", &incr);
			for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
			{
				fscanf(fp, "%" NUM_FMT"", Arrptr->paerial + i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz);
				if (AreEqual(Arrptr->paerial[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz], no_data_value))
					Arrptr->paerial[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz] = 1;
			}
		}
		fclose(fp);
	}


	// Por_Ident=3
	else if (Parptr->Por_Ident == 3)
	{
		fscanf(fp, "%s", dum);

		Arrptr->paerial = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->pbound = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz * 4);

		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->paerial + i + j*Parptr->xsz);
			if (AreEqual(Arrptr->paerial[i + j*Parptr->xsz], no_data_value))
				Arrptr->paerial[i + j*Parptr->xsz] = 1;
		}

		fscanf(fp, "%s", dum);

		for (k = 0; k < 4; k++)
		{
			for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
			{
				fscanf(fp, "%" NUM_FMT"", Arrptr->pbound + i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz);
				if (AreEqual(Arrptr->pbound[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz], no_data_value))
					Arrptr->pbound[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz] = 1;
			}
		}
		fclose(fp);
	}

	// Por_Ident=4
	else if (Parptr->Por_Ident == 4)
	{
		fscanf(fp, "%s %i", dum, &Parptr->zsz);
		fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->maxelev);
		fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->zlev);
		fscanf(fp, "%s", dum);

		Arrptr->paerial = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz*Parptr->zsz);
		Arrptr->pbound = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz*Parptr->zsz * 4);

		for (k = 0; k < Parptr->zsz; k++)
		{
			fscanf(fp, "%" NUM_FMT"", &incr);
			for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
			{
				fscanf(fp, "%" NUM_FMT"", Arrptr->paerial + i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz);
				if (AreEqual(Arrptr->paerial[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz], no_data_value))
					Arrptr->paerial[i + j*Parptr->xsz + k*Parptr->xsz*Parptr->ysz] = 1;
			}
		}

		fscanf(fp, "%s", dum);

		for (m = 0; m < Parptr->zsz; m++)
		{
			fscanf(fp, "%" NUM_FMT"", &incr);
			for (k = 0; k < 4; k++) for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
			{
				fscanf(fp, "%" NUM_FMT"", Arrptr->pbound + i + j*Parptr->xsz + m*Parptr->xsz*Parptr->ysz + k*Parptr->xsz*Parptr->ysz);
				if (AreEqual(Arrptr->pbound[i + j*Parptr->xsz + m*Parptr->xsz*Parptr->ysz + k*Parptr->xsz*Parptr->ysz], no_data_value))
				{
					Arrptr->pbound[i + j*Parptr->xsz + m*Parptr->xsz*Parptr->ysz + k*Parptr->xsz*Parptr->ysz] = 1;
				}
			}
		}
		fclose(fp);
	}

	if (verbose == ON) printf("Done.\n\n");

	return;
}

void InitStartFile(States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr,
	const NUMERIC_TYPE no_data_value,
	const int verbose)
{
	int gr;
	for (int j = 0; j < Parptr->ysz; j++)
		for (int i = 0; i < Parptr->xsz; i++)
		{
			//fscanf(fp, "%" NUM_FMT"", Arrptr->H + i + j*Parptr->xsz);
			// if no_data set depth to zero
			if (AreEqual(Arrptr->H[i + j*Parptr->xsz], no_data_value)) Arrptr->H[i + j*Parptr->xsz] = C(0.0);
			else if (Statesptr->startelev == ON) // convert water surface elevation to depth is this is being used
			{
				// check to see if SGC is on
				if (Statesptr->SGC == ON)
				{
					gr = Arrptr->SGCgroup[i + j*Parptr->xsz]; // channel group number
					// is SGC is on calculate both a depth from the channel bed and the domain volume
					Arrptr->H[i + j*Parptr->xsz] = getmax(Arrptr->H[i + j*Parptr->xsz] - Arrptr->SGCz[i + j*Parptr->xsz], C(0.0));

					if (Arrptr->H[i + j*Parptr->xsz] <= Arrptr->SGCbfH[i + j*Parptr->xsz])
						Arrptr->SGCVol[i + j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i + j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i + j*Parptr->xsz]);
					else if (Arrptr->SGCwidth[i + j*Parptr->xsz] >= C(0.5)*(Arrptr->dx[i + j*Parptr->xsz] + Arrptr->dy[i + j*Parptr->xsz]))
						Arrptr->SGCVol[i + j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i + j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i + j*Parptr->xsz]);
					else
						Arrptr->SGCVol[i + j*Parptr->xsz] = Arrptr->SGCbfV[i + j*Parptr->xsz] + (Arrptr->H[i + j*Parptr->xsz] 
							- Arrptr->SGCbfH[i + j*Parptr->xsz])*(Arrptr->dx[i + j*Parptr->xsz] * Arrptr->dy[i + j*Parptr->xsz]); // out of bank level * cell area
				}
				else
					Arrptr->H[i + j*Parptr->xsz] = getmax(Arrptr->H[i + j*Parptr->xsz] - Arrptr->DEM[i + j*Parptr->xsz], C(0.0));

			}
			else if (Statesptr->SGC == ON)
			{
				gr = Arrptr->SGCgroup[i + j*Parptr->xsz]; // channel group number
				if (Arrptr->H[i + j*Parptr->xsz] <= Arrptr->SGCbfH[i + j*Parptr->xsz])
					Arrptr->SGCVol[i + j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i + j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i + j*Parptr->xsz]);
				else if (Arrptr->SGCwidth[i + j*Parptr->xsz] >= C(0.5)*(Arrptr->dx[i + j*Parptr->xsz] + Arrptr->dy[i + j*Parptr->xsz]))
					Arrptr->SGCVol[i + j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i + j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i + j*Parptr->xsz]);
				else
					Arrptr->SGCVol[i + j*Parptr->xsz] = Arrptr->SGCbfV[i + j*Parptr->xsz] + (Arrptr->H[i + j*Parptr->xsz] - Arrptr->SGCbfH[i + j*Parptr->xsz])*(Arrptr->dx[i + j*Parptr->xsz] * Arrptr->dy[i + j*Parptr->xsz]); // out of bank level
			}
		}

}

//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM FILE
void LoadStart(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, const int verbose)
{
	FILE *fp;
	int i, j;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->startfilename, "r", "Loading startfile ", verbose);

	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	for (j = 0; j < Parptr->ysz; j++)
		for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->H + i + j*Parptr->xsz);
		}
	fclose(fp);

	InitStartFile(Statesptr, Parptr, Arrptr, SGCptr, no_data_value, verbose);

	if (verbose == ON) printf("Done.\n\n");

	return;
}


//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM BINARY FILE
void LoadBinaryStart(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, const int verbose)
{
	FILE *fp;
	int i, j, tmpi;
	NUMERIC_TYPE tmp, no_data_value = -9999;

	fp = fopen_or_die(Fnameptr->startfilename, "r", "Loading initial binary startfile", verbose);

	// read and dump header information, but check data are compatable
	fread(&tmpi, sizeof(int), 1, fp);
	if (tmpi != Parptr->xsz) printf("\nWARNING: incorrect number of cells in .start file\n");
	fread(&tmpi, sizeof(int), 1, fp);
	if (tmpi != Parptr->ysz) printf("\nWARNING: incorrect number of cells in .start file\n");

	fread(&tmp, sizeof(NUMERIC_TYPE), 1, fp);
	fread(&tmp, sizeof(NUMERIC_TYPE), 1, fp);
	fread(&tmp, sizeof(NUMERIC_TYPE), 1, fp);
	fread(&no_data_value, sizeof(NUMERIC_TYPE), 1, fp);

	for (j = 0; j < Parptr->ysz; j++) //for(i=0;i<Parptr->xsz;i++) no need to loop x with fread.
	{
		fread(Arrptr->H + j*Parptr->xsz, sizeof(NUMERIC_TYPE), Parptr->xsz, fp);
		//// loop through x
		//for (i = 0; i < Parptr->xsz; i++)
		//{
		//	// Set depth to zero if no_data
		//	if ((int)Arrptr->H[i + j*Parptr->xsz] == no_data_value) Arrptr->H[i + j*Parptr->xsz] = C(0.0);
		//	else if (Statesptr->startelev == ON) // if elevation file used for initial depths this needs to be converted to depths
		//	{
		//		// check to see if SGC is on
		//		if (Statesptr->SGC == ON) Arrptr->H[i + j*Parptr->xsz] = getmax(Arrptr->H[i + j*Parptr->xsz] - Arrptr->SGCz[i + j*Parptr->xsz], C(0.0));
		//		else Arrptr->H[i + j*Parptr->xsz] = getmax(Arrptr->H[i + j*Parptr->xsz] - Arrptr->DEM[i + j*Parptr->xsz], C(0.0));
		//	}
		//}
	}
	fclose(fp);

	InitStartFile(Statesptr, Parptr, Arrptr, SGCptr, no_data_value, verbose);

	if (verbose == ON) printf("Done.\n\n");

	return;
}

//----------------------------------------------------------------------------
// LOAD <startfile>.Qx and <startfile>.Qy into HU and HV for FV1 and DG2
// solvers.
void LoadStartQ2D(Fnames* Fnameptr, Pars* Parptr, Arrays* Arrptr,
        const int verbose)
{
	FILE *fp;
	int i, j;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;

    char hu_startfile[800];
    strcpy(hu_startfile, Fnameptr->startfilename);
    strcat(hu_startfile, ".Qx");

    char hv_startfile[800];
    strcpy(hv_startfile, Fnameptr->startfilename);
    strcat(hv_startfile, ".Qy");

	fp = fopen_or_die(hu_startfile, "r", "Loading startfile.Qx ", verbose);
	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
	for (j = 0; j < Parptr->ysz; j++)
    {
		for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->HU + i + j*Parptr->xsz);
		}
    }
	fclose(fp);
	if (verbose == ON) printf("Done.\n\n");

	fp = fopen_or_die(hv_startfile, "r", "Loading startfile.Qy ", verbose);
	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
	for (j = 0; j < Parptr->ysz; j++)
    {
		for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->HV + i + j*Parptr->xsz);
		}
    }
	fclose(fp);
	if (verbose == ON) printf("Done.\n\n");

	return;
}

//----------------------------------------------------------------------------
// LOAD DEM FROM FILE
// Also loads cell size, lower left corner coordinates, and new all rasters
void LoadDEM(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, const int verbose)
{
	int i, j;
	NUMERIC_TYPE no_data_value = -9999;

	FILE *fp = LoadDomainGeometry(Fnameptr->demfilename, Parptr, verbose,
			 no_data_value);

	// allocate memory for arrays, Note the () at the end ensures all elements are initialised to zero
	Arrptr->H = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	// If Roe==on allocate memory (JN/IV) 
	if (Statesptr->Roe == ON)
	{
		Arrptr->HU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->HV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->LSHU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->RSHU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->BSHU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->TSHU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->LSHV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->RSHV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->BSHV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->TSHV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->FHx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->FHUx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->FHVx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->FHy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->FHUy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->FHVy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	}
	else if (Statesptr->fv1 == ON)
	{
		allocate_swe_fields(Parptr, Arrptr);
	}
	else if (Statesptr->dg2 == ON)
	{
		dg2::allocate_fields(Parptr, Arrptr);
	}

	Arrptr->maxH = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->totalHtm = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->Qx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->Qy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	// added to record Hflow (needed for velocity and hazard calculations
	//Arrptr->Hflowx=memory_allocate_zero_numeric_legacy((Parptr->xsz+1)*(Parptr->ysz+1));
	//Arrptr->Hflowy=memory_allocate_zero_numeric_legacy((Parptr->xsz+1)*(Parptr->ysz+1));

	if (Statesptr->voutput == ON)
	{
		Arrptr->Vx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->Vy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->maxVx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->maxVy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	}
	if (Statesptr->hazard == ON)
	{
		Arrptr->maxVc = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->maxVcH = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->maxHaz = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	}
	if (Statesptr->SGC == ON)
	{
		// Geometric variables
		Arrptr->SGCwidth = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->SGCz = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->SGCc = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);

		// Flow variables
		Arrptr->QxSGold = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		Arrptr->QySGold = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
		// Bank full variables
		Arrptr->SGCbfH = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->SGCbfV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		// Volume variables
		Arrptr->SGCVol = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		Arrptr->SGCdVol = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		// Model parameters
		Arrptr->SGCgroup = new int[Parptr->xsz*Parptr->ysz]();
		if (Statesptr->save_Qs == ON)
		{
			Arrptr->SGCFlowWidth = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
		}
	}

	Arrptr->Qxold = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->Qyold = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));

	// allocate memory for velocity arrays U and V
	// currently only used in acceleration version - initialised under all conditions as may want them 
	// for other lisflood versions (TJF)
	Arrptr->U = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	//Arrptr->V=memory_allocate_zero_numeric_legacy((Parptr->xsz+1)*(Parptr->ysz+1));

	// allocate memory for none zero arrays
	Arrptr->maxHtm = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->initHtm = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->TRecx = memory_allocate_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->TRecy = memory_allocate_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->LimQx = memory_allocate_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->LimQy = memory_allocate_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));

	Arrptr->ChanMask = new int[Parptr->xsz*Parptr->ysz];
	Arrptr->SegMask = new int[Parptr->xsz*Parptr->ysz];

	Arrptr->DEM = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);

	// allocate memory for flow direction array for routing very shallow flows from rainfall component CCS 13/03/2012
	if (Statesptr->routing == ON)
	{
		Arrptr->FlowDir = new int[Parptr->xsz*Parptr->ysz]();
		for (i = 0; i < Parptr->xsz*Parptr->ysz; i++) Arrptr->FlowDir[i] = (int)NULLVAL;
	}

	// allocate memory for lat long arrays
	Arrptr->dx = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->dy = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->dA = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);



	// set initial values of elements for some arrays to NULLVAL
	for (i = 0; i < Parptr->xsz*Parptr->ysz; i++) Arrptr->maxHtm[i] = NULLVAL;
	for (i = 0; i < Parptr->xsz*Parptr->ysz; i++) Arrptr->initHtm[i] = NULLVAL;
	for (i = 0; i < (Parptr->xsz + 1)*(Parptr->ysz + 1); i++) Arrptr->TRecx[i] = NULLVAL;
	for (i = 0; i < (Parptr->xsz + 1)*(Parptr->ysz + 1); i++) Arrptr->TRecy[i] = NULLVAL;
	for (i = 0; i < (Parptr->xsz + 1)*(Parptr->ysz + 1); i++) Arrptr->LimQx[i] = NULLVAL;
	for (i = 0; i < (Parptr->xsz + 1)*(Parptr->ysz + 1); i++) Arrptr->LimQy[i] = NULLVAL;

	// set initial values of elements for mask arrays to NULLVAL
	for (i = 0; i < Parptr->xsz*Parptr->ysz; i++) Arrptr->ChanMask[i] = -1;
	for (i = 0; i < Parptr->xsz*Parptr->ysz; i++) Arrptr->SegMask[i] = -1;

	LoadDEMData(Parptr, Arrptr->DEM, fp, no_data_value);
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");

	return;
}

FILE* LoadDomainGeometry
(
	const char* filename,
	Pars *Parptr,
	const int verbose,
	NUMERIC_TYPE& no_data_value
)
{
	char dum[800];

	FILE *fp = fopen_or_die(filename, "rb", "Loading DEM", verbose);

	fscanf(fp, "%s %i", dum, &Parptr->xsz);
	fscanf(fp, "%s %i", dum, &Parptr->ysz);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->blx);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->bly);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->dx);

	Parptr->dx_sqrt = sqrt((NUMERIC_TYPE)Parptr->dx); // sqrt now for later use in flooplain calcs - small speed increase
	Parptr->dy = Parptr->dx; Parptr->dA = Parptr->dx*Parptr->dy;
	Parptr->tlx = Parptr->blx; Parptr->tly = Parptr->bly + Parptr->ysz*Parptr->dy;
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	if (verbose == ON)
	{
		printf("%ix%i\nBL corner\t(%" NUM_FMT",%" NUM_FMT")\nNODATA_value\t%" NUM_FMT"\n",
			Parptr->xsz, Parptr->ysz, Parptr->blx, Parptr->bly, no_data_value);
	}

	return fp;
}

void LoadDEMData
(
	Pars *Parptr,
	NUMERIC_TYPE *DEM,
	FILE *fp,
	NUMERIC_TYPE file_nodata_value
)
{
	for (int j=0; j<Parptr->ysz; j++)
	{
		for (int i=0; i<Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", &DEM[i + j*Parptr->xsz]);
			if (AreEqual(DEM[i + j*Parptr->xsz], file_nodata_value))
			{
				DEM[i + j*Parptr->xsz] = Parptr->nodata_elevation;
			}
		}
	}
}

//-----------------------------------------------------------------------------------
// LOADS FILE GIVING IDENTIFIERS FOR EACH BOUNDARY CELL FROM .bci FILE
// (e.g. HFIX, QVAR etc)
void LoadBCs(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, BoundCs *BCptr, const int verbose)
{
	int numBCs, i, j, BCi1, BCi2, tmpi;
	NUMERIC_TYPE start, finish;
	FILE *fp;
	char buff[800], buff2[800], buff3[800], side;
	NUMERIC_TYPE BC_tmp;
	int pi = -1, maxpi = 20000;  // increase max from 10 to 20K (MT)
	NUMERIC_TYPE px, py;

	int *new_xpi, *new_ypi;
	ESourceType *new_PS_Ident;
	NUMERIC_TYPE *new_PS_Val, *new_PS_Q_FP_old, *new_PS_Q_SG_old;
	char *new_PS_Name;

	// POINT SOURCE STUFF
	BCptr->xpi = new int[maxpi];
	BCptr->ypi = new int[maxpi];
	for (i = 0; i < maxpi; i++)
		BCptr->xpi[i] = BCptr->ypi[i] = -1;
	BCptr->PS_Ident = new ESourceType[maxpi];
	BCptr->PS_Val = memory_allocate_numeric_legacy(maxpi);
	BCptr->PS_Name = new char[maxpi * 80];
	for (i = 0; i < maxpi; i++)
	{
		BCptr->PS_Ident[i] = NONE0;
		BCptr->PS_Val[i] = C(-1.0);
		BCptr->PS_Name[i] = '\0';
	}
	BCptr->numPS = -1;

	// BOUNDARY CONDITION STUFF
	numBCs = 2 * Parptr->xsz + 2 * Parptr->ysz;
	BCptr->BC_Ident = new ESourceType[numBCs];
	BCptr->BC_Val = memory_allocate_numeric_legacy(numBCs);
	BCptr->BC_Name = new char[numBCs * 80];
	BCptr->numBCs = numBCs;

	for (i = 0; i < numBCs; i++)
	{
		BCptr->BC_Ident[i] = NONE0;
		BCptr->BC_Val[i] = C(-1.0);
		BCptr->BC_Name[i*80] = '\0';
	}

	if (strlen(Fnameptr->bcifilename) == 0)
	{
		if (verbose == ON) printf("No bcifile used\n");
		return;
	}
	fp = fopen_or_die(Fnameptr->bcifilename, "rb", "Loading boundary condition IDs\n", verbose);

	while (!feof(fp))
	{
		BCi1 = BCi2 = -1; side = '\0';
		// Read NSEW and location, and determine start/finish of BC-->(BCi1,BCi2)
		fscanf(fp, "%s", buff);
		if (feof(fp)) break;
		if (buff[0] == 'N')
		{
			fscanf(fp, "%" NUM_FMT"%" NUM_FMT"", &start, &finish);

			if (start < Parptr->blx) start = Parptr->blx;
			if (start > Parptr->blx + Parptr->xsz*Parptr->dx) start = Parptr->blx + Parptr->xsz*Parptr->dx;
			if (finish < Parptr->blx) finish = Parptr->blx;
			if (finish > Parptr->blx + Parptr->xsz*Parptr->dx) finish = Parptr->blx + Parptr->xsz*Parptr->dx;

			BCi1 = (int)((start - Parptr->blx) / Parptr->dy);
			BCi2 = (int)((finish - Parptr->blx) / Parptr->dy);
			if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
			BCi2--; side = 'N';
		}

		if (buff[0] == 'W')
		{
			fscanf(fp, "%" NUM_FMT"%" NUM_FMT"", &start, &finish);

			if (start < Parptr->bly) start = Parptr->bly;
			if (start > Parptr->tly) start = Parptr->tly;
			if (finish < Parptr->bly) finish = Parptr->bly;
			if (finish > Parptr->tly) finish = Parptr->tly;

			BCi1 = (int)(2 * Parptr->xsz + 2 * Parptr->ysz - (Parptr->tly - start) / Parptr->dy);
			BCi2 = (int)(2 * Parptr->xsz + 2 * Parptr->ysz - (Parptr->tly - finish) / Parptr->dy);
			if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
			BCi2--; side = 'W';
		}

		if (buff[0] == 'S')
		{
			fscanf(fp, "%" NUM_FMT"%" NUM_FMT"", &start, &finish);

			if (start < Parptr->blx) start = Parptr->blx;
			if (start > Parptr->blx + Parptr->xsz*Parptr->dx) start = Parptr->blx + Parptr->xsz*Parptr->dx;
			if (finish < Parptr->blx) finish = Parptr->blx;
			if (finish > Parptr->blx + Parptr->xsz*Parptr->dx) finish = Parptr->blx + Parptr->xsz*Parptr->dx;

			BCi1 = (int)(2 * Parptr->xsz + Parptr->ysz - (start - Parptr->blx) / Parptr->dy);
			BCi2 = (int)(2 * Parptr->xsz + Parptr->ysz - (finish - Parptr->blx) / Parptr->dy);

			if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
			BCi2--; side = 'S';
		}

		if (buff[0] == 'E')
		{
			fscanf(fp, "%" NUM_FMT"%" NUM_FMT"", &start, &finish);

			if (start < Parptr->bly) start = Parptr->bly;
			if (start > Parptr->tly) start = Parptr->tly;
			if (finish < Parptr->bly) finish = Parptr->bly;
			if (finish > Parptr->tly) finish = Parptr->tly;

			BCi1 = (int)(Parptr->xsz + (Parptr->tly - start) / Parptr->dy);
			BCi2 = (int)(Parptr->xsz + (Parptr->tly - finish) / Parptr->dy);
			if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
			BCi2--; side = 'E';
		}

		// Read locations of point sources or point free
		if (buff[0] == 'P' || buff[0] == 'F')
		{
			pi++;
			fscanf(fp, "%" NUM_FMT"%" NUM_FMT"", &px, &py);
			BCptr->xpi[pi] = (int)((px - Parptr->blx) / Parptr->dx);
			BCptr->ypi[pi] = (int)((Parptr->tly - py) / Parptr->dy);
		}

		// Read free boundary condition locations
		// load buffer until EOL
		j = 0;
		do{ buff2[j] = fgetc(fp); } while (buff2[j++] != '\n' && !feof(fp));
		if (j > 0)
			buff2[j - 1] = '\0';               // Finish off string
		// get buff so you know boundary type
		BC_tmp = -1;
		sscanf(buff2, "%s%" NUM_FMT"", buff, &BC_tmp);

		// If a FREE surface boundary condition
		if (!STRCMPi(buff, "FREE") && BCi1 > -1)
		{
			// if BC_tmp is -1 there is no slope specifed... use local slope from elevation model (origional mehod)
			if (BC_tmp < C(-0.999)) // ie -1 (done like this as NUMERIC_TYPE)
			{
				if (verbose == ON) printf("FREE on %c side start %" NUM_FMT" end %" NUM_FMT"\n", side, start, finish);
			}
			else
			{
				if (verbose == ON) printf("FREE on %c side start %" NUM_FMT" end %" NUM_FMT" using slope %.5" NUM_FMT"\n", side, start, finish, BC_tmp);
			}
			for (i = BCi1; i <= BCi2; i++)
			{
				BCptr->BC_Ident[i] = FREE1;
				// store floodplain slope in BCptr->BC_Val[i]
				if (Statesptr->adaptive_ts == ON || Statesptr->qlim == ON)
				{
					if (BC_tmp < C(-0.999)) BCptr->BC_Val[i] = BC_tmp; // make BC_Val equal to -1 to use local water surface slope (origional lisflood)
					else BCptr->BC_Val[i] = sqrt(BC_tmp); // sqrt of user specified slope for diffusive version (jcn)
				}
				else
				{
					BCptr->BC_Val[i] = BC_tmp; // user specified slope or -1(for local slope) for accelleration or any other version (jcn)
				}
			}
		}
		// Read in fixed values of H for boundary
		else if (!STRCMPi(buff, "HFIX") && BCi1 > -1)
		{
			sscanf(buff2, "%s%" NUM_FMT"", buff, &BC_tmp);
			if (verbose == ON) printf("HFIX at %" NUM_FMT" on %c side start %" NUM_FMT" end %" NUM_FMT"\n", BC_tmp, side, start, finish);

			for (i = BCi1; i <= BCi2; i++)
			{
				BCptr->BC_Val[i] = BC_tmp;
				BCptr->BC_Ident[i] = HFIX2;
			}
		}
		// Read in fixed values if Q for boundary
		else if (!STRCMPi(buff, "QFIX") && BCi1 > -1)
		{
			sscanf(buff2, "%s%" NUM_FMT"", buff, &BC_tmp);
			if (verbose == ON) printf("QFIX at %" NUM_FMT" on %c side start %" NUM_FMT" end %" NUM_FMT"\n", BC_tmp, side, start, finish);

			for (i = BCi1; i <= BCi2; i++)
			{
				BCptr->BC_Val[i] = BC_tmp;
				BCptr->BC_Ident[i] = QFIX4;
			}
		}

		//	Read boundary names for varying values
		else if (!STRCMPi(buff, "QVAR") && BCi1 > -1)
		{
			//fscanf(fp,"%s",buff);
			sscanf(buff2, "%s%s", buff3, buff);
			if (verbose == ON)
				printf("QVAR from bdy file %s on %c side start %" NUM_FMT" end %" NUM_FMT"\n",
				buff, side, start, finish);
			for (i = BCi1; i <= BCi2; i++)
			{
				strcpy(BCptr->BC_Name + i * 80, buff);
				BCptr->BC_Ident[i] = QVAR5;
			}
		}
		else if (!STRCMPi(buff, "HVAR") && BCi1 > -1)
		{
			//fscanf(fp,"%s",buff);
			sscanf(buff2, "%s%s", buff3, buff);
			if (verbose == ON)
				printf("HVAR from bdy file %s on %c side start %" NUM_FMT" end %" NUM_FMT"\n", buff, side, start, finish);
			for (i = BCi1; i <= BCi2; i++)
			{
				strcpy(BCptr->BC_Name + i * 80, buff);
				BCptr->BC_Ident[i] = HVAR3;
			}
		}
		// Fixed/Varying values/names for point sources
		// Note these need to come after the boundary conditions in the code else both will get implemented! 
		// I'm not convinced &&BCptr->xpi[pi]>-1&&pi>-1 is strict enough (JCN) 
		else if (!STRCMPi(buff, "HFIX") && BCptr->xpi[pi] > -1 && pi > -1)
		{
			//fscanf(fp,"%" NUM_FMT"",&BC_tmp);
			sscanf(buff2, "%s%s", buff3, buff);
			if (verbose == ON) printf("HFIX at point [%" NUM_FMT",%" NUM_FMT"] (%d,%d) %" NUM_FMT"\n", px, py, BCptr->xpi[pi], BCptr->ypi[pi], BC_tmp);
			BCptr->PS_Val[pi] = BC_tmp;
			BCptr->PS_Ident[pi] = HFIX2;
		}
		else if (!STRCMPi(buff, "QFIX") && BCptr->xpi[pi] > -1 && pi > -1)
		{
			//fscanf(fp,"%" NUM_FMT"",&BC_tmp);
			if (verbose == ON) printf("QFIX at point [%" NUM_FMT",%" NUM_FMT"] (%d,%d) %" NUM_FMT"\n", px, py, BCptr->xpi[pi], BCptr->ypi[pi], BC_tmp);
			BCptr->PS_Val[pi] = BC_tmp;
			BCptr->PS_Ident[pi] = QFIX4;
		}
		else if (!STRCMPi(buff, "QVAR") && BCptr->xpi[pi] > -1 && pi > -1)
		{
			//fscanf(fp,"%s",buff);
			sscanf(buff2, "%s%s", buff3, buff);
			if (verbose == ON)
				printf("QVAR at point [%" NUM_FMT",%" NUM_FMT"] (%d,%d) %s\n", px, py, BCptr->xpi[pi], BCptr->ypi[pi], buff);
			strcpy(BCptr->PS_Name + pi * 80, buff);
			BCptr->PS_Ident[pi] = QVAR5;
		}
		else if (!STRCMPi(buff, "HVAR") && BCptr->xpi[pi] > -1 && pi > -1)
		{
			//fscanf(fp,"%s",buff);
			sscanf(buff2, "%s%s", buff3, buff);
			if (verbose == ON) printf("HVAR at point [%" NUM_FMT",%" NUM_FMT"] (%d,%d) %s\n", px, py, BCptr->xpi[pi], BCptr->ypi[pi], buff);
			strcpy(BCptr->PS_Name + pi * 80, buff);
			BCptr->PS_Ident[pi] = HVAR3;
		}
		else if (!STRCMPi(buff, "FREE") && BCptr->xpi[pi] > -1 && pi > -1 && Statesptr->SGC == ON)
		{
			// point FREE boundary for internal SGC boundaries
			if (verbose == ON) printf("FREE at point [%" NUM_FMT",%" NUM_FMT"] (%d,%d) with slope %" NUM_FMT"\n", px, py, BCptr->xpi[pi], BCptr->ypi[pi], BC_tmp);
			BCptr->PS_Val[pi] = BC_tmp;
			BCptr->PS_Ident[pi] = FREE6;
		}
		else
		{
			if (verbose == ON) printf("WARNING: Incorrect boundary condition in .bci file\n");
		}

	}

	// allocate new buffer to the size of the data
	if (pi > -1)
	{
		pi++;
		new_xpi = new int[pi];
		new_ypi = new int[pi];
		new_PS_Ident = new ESourceType[pi];
		new_PS_Val = memory_allocate_numeric_legacy(pi);
		new_PS_Q_FP_old = memory_allocate_zero_numeric_legacy(pi);
		new_PS_Q_SG_old = memory_allocate_zero_numeric_legacy(pi);
		new_PS_Name = new char[pi * 80];


		for (i = 0; i < pi; i++)
		{
			new_xpi[i] = BCptr->xpi[i];
			new_ypi[i] = BCptr->ypi[i];
			new_PS_Ident[i] = BCptr->PS_Ident[i];
			new_PS_Val[i] = BCptr->PS_Val[i];
			for (j = 0; j < 80; j++)
				new_PS_Name[i * 80 + j] = BCptr->PS_Name[i * 80 + j];
		}
		delete[] BCptr->xpi;
		delete[] BCptr->ypi;
		delete[] BCptr->PS_Ident;
		memory_free_legacy(&BCptr->PS_Val);
		delete[] BCptr->PS_Name;

		BCptr->xpi = new_xpi;
		BCptr->ypi = new_ypi;
		BCptr->PS_Ident = new_PS_Ident;
		BCptr->PS_Val = new_PS_Val;
		BCptr->PS_Name = new_PS_Name;

		BCptr->PS_Q_FP_old = new_PS_Q_FP_old;
		BCptr->PS_Q_SG_old = new_PS_Q_SG_old;

		BCptr->numPS = pi;
	}

	if (verbose == ON) printf("Done.\n\n");
	//  LoadBCVar();

	fclose(fp);

	return;
}

void FreeTimeSeries(TimeSeries& timeSeries)
{
	NUMERIC_TYPE ** mem;
	mem = &(timeSeries.time);
	memory_free((void**)mem);
	mem = &(timeSeries.value);
	memory_free((void**)mem);
}

void LoadTimeSeries
(
	TimeSeries& timeSeries,
	const char* filename,
	FILE * fp,
	int skipFirstLine
)
{
	int ndata = -1;
	char units[80];
	char line_buffer[LINE_BUFFER_LEN];
	int line_number = 0;
	while (fgets(line_buffer, LINE_BUFFER_LEN, fp))
	{
		line_number++;
		int line_len = strlen(trimwhitespace(line_buffer));
		if (line_len == LINE_BUFFER_LEN - 1)
		{
			printf("Time Series line too long line %s\n", filename);
			exit(-1);
		}
		if ((skipFirstLine == ON && line_number == 1) || line_len == 0 || line_buffer[0] == '#')
		{
			continue;
		}
		sscanf(line_buffer, "%i%s", &ndata, units);
		break;
	}
	if (ndata == -1)
	{
		printf("Invalid Time Series %s\n", filename);
		exit(-1);
	}

	timeSeries.count = ndata;
	timeSeries.prev_index = 0;
	timeSeries.prev_time = C(-1.0);
	timeSeries.time = (NUMERIC_TYPE*)memory_allocate(ndata*sizeof(NUMERIC_TYPE));
	timeSeries.value = (NUMERIC_TYPE*)memory_allocate(ndata*sizeof(NUMERIC_TYPE));

	int i = 0;
	NUMERIC_TYPE prev_time = C(-1.);
	while (i < ndata && fgets(line_buffer, LINE_BUFFER_LEN, fp))
	{
		int line_len = strlen(line_buffer);
		if (line_len == LINE_BUFFER_LEN - 1)
		{
			printf("Time Series line too long line %s\n", filename);
			exit(-1);
		}
		if (line_len == 0 || line_buffer[0] == '#')
		{
			continue;

		}
		sscanf(line_buffer, "%" NUM_FMT"%" NUM_FMT"", &timeSeries.value[i], &timeSeries.time[i]);

		if (timeSeries.time[i] <= prev_time)
		{
			printf("Time Series invalid time values - times should be in increasing order. %s\n", filename);
			exit(-1);
		}
		prev_time = timeSeries.time[i];

		i++;
	}
	if (i < ndata)
	{
		printf("Time Series 'count' is greater than the number of values in the file %s\n", filename);
		exit(-1);
	}

	if (!STRCMPi(units, "minutes") || !STRCMPi(units, "minute"))
		for (i = 0; i < ndata; i++)
			timeSeries.time[i] *= 60;
	else if (!STRCMPi(units, "hours") || !STRCMPi(units, "hour"))
		for (i = 0; i < ndata; i++)
			timeSeries.time[i] *= 3600;
	else if (!STRCMPi(units, "days") || !STRCMPi(units, "day"))
		for (i = 0; i < ndata; i++)
			timeSeries.time[i] *= (3600 * 24);
	else if (!(!STRCMPi(units, "seconds") || !STRCMPi(units, "second"))) // if not seconds, report warning
	{
		printf("Time Series unknown units '%s' expect 'seconds', 'minutes', 'hours', 'days' (use input as seconds)\n", units);
	}
}

TimeSeries * LoadTimeSeries(const char* filename, FILE * fp, int skipFirstLine)
{
	TimeSeries * timeSeries = new TimeSeries();
	LoadTimeSeries(*timeSeries, filename, fp, skipFirstLine);
	return timeSeries;
}

//-----------------------------------------------------------------------------
// LOAD TIME VARYING BOUNDARY CONDITIONS FROM .bdy FILE
void LoadBCVar(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, BoundCs *BCptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, const int verbose)
{
	FILE *fp;
	int i, numBCs, chseg;
	char line_buffer[LINE_BUFFER_LEN];

	char name_buffer[255], units[80];

	numBCs = 2 * Parptr->xsz + 2 * Parptr->ysz;
	BCptr->BC_TimeSeries = new TimeSeries*[numBCs];

	for (i = 0; i < numBCs; i++)
	{
		BCptr->BC_TimeSeries[i] = NULL;
	}

	if (BCptr->numPS > 0)
	{
		BCptr->PS_TimeSeries = new TimeSeries*[BCptr->numPS];
		for (i = 0; i < BCptr->numPS; i++)
		{
			BCptr->PS_TimeSeries[i] = NULL;
		}
	}

	if (Statesptr->ChannelPresent == ON) {
		for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
		{
			ChannelSegments[chseg].Q_TimeSeries = new TimeSeries*[ChannelSegments[chseg].chsz];
			for (i = 0; i < ChannelSegments[chseg].chsz; i++)
				ChannelSegments[chseg].Q_TimeSeries[i] = NULL;
		}
	}

	if (strlen(Fnameptr->bdyfilename) != 0)
	{
		fp = fopen_or_die(Fnameptr->bdyfilename, "r", "Loading time varying boundary conditions\n", verbose);

		int line_number = 0;
		while (fgets(line_buffer, LINE_BUFFER_LEN, fp))
		{
			line_number++;

			int line_len = strlen(trimwhitespace(line_buffer));
			if (line_len == LINE_BUFFER_LEN - 1)
			{
				printf("bdyfile file line too long line %s\n", Fnameptr->bdyfilename);
				exit(-1);
			}
			// skip 1st, any comment line and any empty line
			if (line_number == 1 || line_len == 0 || line_buffer[0] == '#')
				continue;

			sscanf(line_buffer, "%s", name_buffer);

			/* older code (incompatible with GPU solvers)
			TimeSeries* timeSeries = new TimeSeries();
			BCptr->allTimeSeries.push_back(*timeSeries);
			LoadTimeSeries(*timeSeries, Fnameptr->bdyfilename, fp, OFF);
			*/

			// new code
			TimeSeries* timeSeries = new TimeSeries();
			BCptr->allTimeSeries.push_back(*timeSeries);
			LoadTimeSeries(BCptr->allTimeSeries.back(), Fnameptr->bdyfilename, fp, OFF); 
			*timeSeries = BCptr->allTimeSeries.back(); 

			/* James shaw code (incompatible with CPU solvers)
			BCptr->allTimeSeries.push_back(TimeSeries());
			TimeSeries& timeSeries = BCptr->allTimeSeries.back();
			LoadTimeSeries(timeSeries, Fnameptr->bdyfilename, fp, OFF);
			*/

			int timeSeriesUsed = 0;

			// Check through list (2d domain) of boundary names, if found, point to the loaded time series.
			for (i = 0; i < numBCs; i++)
			{
				if (!strcmp(name_buffer, (BCptr->BC_Name + i * 80)))
				{
					BCptr->BC_Val[i] = C(0.0);
					BCptr->BC_TimeSeries[i] = timeSeries; // OLD VERSION
					//BCptr->BC_TimeSeries[i] = &timeSeries; // JAMES SHAW VERSION
					timeSeriesUsed++;
				}
			}

			// Check through list (river channel) of boundary names, if found, point to the loaded time series.
			if (Statesptr->ChannelPresent == ON)
			{
				for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
				{
					for (i = 0; i < ChannelSegments[chseg].chsz; i++)
					{
						if (!strcmp(name_buffer, (ChannelSegments[chseg].Q_Name + i * 80)))
						{
							ChannelSegments[chseg].Q_Val[i] = C(0.0);
							ChannelSegments[chseg].Q_TimeSeries[i] = timeSeries; // OLD VERSION
							//ChannelSegments[chseg].Q_TimeSeries[i] = &timeSeries; // JAMES SHAW VERSION
							timeSeriesUsed++;
						}
					}
				}
			}

			// Check through point source names, if found, point to the loaded time series.
			for (i = 0; i < BCptr->numPS; i++)
			{
				if (!strcmp(name_buffer, (BCptr->PS_Name + i * 80)))
				{
					BCptr->PS_Val[i] = C(-1.);
					BCptr->PS_TimeSeries[i] = timeSeries; // OLD VERSION
					//BCptr->PS_TimeSeries[i] = &timeSeries; // JAMES SHAW VERSION
					timeSeriesUsed++;
				}
			}

			//    for(i=0;i<ndata*2+2;i++){
			//      printf("\nLoadBCVar: i=%d, peterpointer[0]+i*2=%" NUM_FMT", peterpointer[nbdy]+i*2+1=%" NUM_FMT"",i,*(peterpointer[0]+i*2),*(peterpointer[0]+i*2+1));
			//    }

			if (timeSeriesUsed == 0)
			{
				printf("WARNING: bdy %s is unreferenced - data ignored.\n", name_buffer);

				FreeTimeSeries(*timeSeries); // OLD VERSION
				//FreeTimeSeries(timeSeries); // JAMES SHAW VERSION
				BCptr->allTimeSeries.pop_back();

				continue;
			}

			if (verbose == ON) printf("bdy %s read.\n", name_buffer);
		}
		fclose(fp);
	}
	else
	{
		if (verbose == ON) printf("No bdyfile used\n");

		// Check for bdy names not found in bdy file
		if (Statesptr->ChannelPresent == ON)
		{
			for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++)
			{
				for (i = 0; i < ChannelSegments[chseg].chsz; i++)
				{
					if (ChannelSegments[chseg].Q_Ident[i] == QVAR5 && ChannelSegments[chseg].Q_TimeSeries[i] == NULL) // CCS
					{
						printf("WARNING: bdy %s in river file not found in bdy file - disabled.\n", ChannelSegments[chseg].Q_Name + i * 80);
						ChannelSegments[chseg].Q_Ident[i] = NONE0;
					}
				}
			}
		}

		for (i = 0; i < numBCs; i++)
		{
			if ((BCptr->BC_Ident[i] == QVAR5 || BCptr->BC_Ident[i] == HVAR3) && BCptr->BC_TimeSeries[i] == NULL)
			{
				printf("WARNING: bdy %s in bci file not found in bdy file - disabled.\n", BCptr->BC_Name + i * 80);
				BCptr->BC_Ident[i] = NONE0;
			}
		}

		for (i = 0; i < BCptr->numPS; i++)
		{
			if ((BCptr->PS_Ident[i] == QVAR5 || BCptr->PS_Ident[i] == HVAR3) && BCptr->PS_TimeSeries[i] == NULL)
			{
				printf("WARNING: bdy %s in bci file not found in bdy file - disabled.\n", BCptr->BC_Name + i * 80);
				BCptr->PS_Ident[i] = NONE0;
			}
		}
	}


	if (verbose == ON) printf("Done.\n\n");
	return;
}

//-----------------------------------------------------------------------------
// LOAD TIME EVAPORATION FROM .evap FILE
void LoadEvap(Fnames *Fnameptr, Arrays *Arrptr, const int verbose)
{
	FILE *fp;
	int i;
	char buff[255], units[80];

	fp = fopen_or_die(Fnameptr->evapfilename, "r", "Loading time varying evaporation", verbose);

	Arrptr->evap = LoadTimeSeries(Fnameptr->evapfilename, fp, ON);
	fclose(fp);

	// convert evaporation rate from mm/day to m/second
	for (i = 0; i < Arrptr->evap->count; i++)
		Arrptr->evap->value[i] /= (1000 * 24 * 3600);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
//-----------------------------------------------------------------------------
// LOAD TIME VARYING RAINFALL FROM .rain FILE
void LoadRain(Fnames *Fnameptr, Arrays *Arrptr, const int verbose)
{
	FILE *fp;
	int i, j, ndata;
	char buff[255], units[80];

	fp = fopen_or_die(Fnameptr->rainfilename, "r", "Loading time varying rainfall", verbose);

	Arrptr->rain = LoadTimeSeries(Fnameptr->rainfilename, fp, ON);
	fclose(fp);

	// convert rainfall rate from mm/hr to m/second
	for (i = 0; i < Arrptr->rain->count; i++)
		Arrptr->rain->value[i] /= (1000 * 3600);

	if (verbose == ON) printf("Done.\n\n");
	return;
}
//-----------------------------------------------------------------------------
// LOAD RAINFALL DISTRIBUTION MASK
void LoadRainmask(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, States *Statesptr, const int verbose)
{
	FILE *fp;
	char dum[800];
	NUMERIC_TYPE no_data_value = -9999;
	int i, j;

	fp = fopen_or_die(Fnameptr->rainmaskname, "rb", "Loading distributed rain mask", verbose);

	fscanf(fp, "%s %i", dum, &Parptr->xsz);
	fscanf(fp, "%s %i", dum, &Parptr->ysz);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->blx);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->bly);
	fscanf(fp, "%s %" NUM_FMT"", dum, &Parptr->dx);

	Parptr->dx_sqrt = sqrt((NUMERIC_TYPE)Parptr->dx); // sqrt now for later use in flooplain calcs - small speed increase
	Parptr->dy = Parptr->dx; Parptr->dA = Parptr->dx*Parptr->dy;
	Parptr->tlx = Parptr->blx; Parptr->tly = Parptr->bly + Parptr->ysz*Parptr->dy;
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	Arrptr->Rainmask = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);

	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%" NUM_FMT"", Arrptr->Rainmask + i + j*Parptr->xsz);
		if (AreEqual(Arrptr->Rainmask[i + j*Parptr->xsz], no_data_value))
			Arrptr->Rainmask[i + j*Parptr->xsz] = 0;
	}
	fclose(fp);

	if (verbose == ON) printf("Done.\n\n");

	return;
}
//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM FILE
void LoadSGC(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, States *Statesptr, const int verbose)
{
	FILE *fp;
	int i, j;
	char dum[80];
	NUMERIC_TYPE no_data_value = -9999, tmp;

	if (strlen(Fnameptr->SGCwidthfilename) == 0 && strlen(Fnameptr->SGCbankfilename) == 0)
	{
		if (verbose == ON) printf("Creating empty SGCwidth and bed arrays\t");
		memset(Arrptr->SGCwidth, 0, sizeof(NUMERIC_TYPE)*(Parptr->xsz)*(Parptr->ysz));
		memcpy(Arrptr->SGCz, Arrptr->DEM, sizeof(NUMERIC_TYPE)*(Parptr->xsz)*(Parptr->ysz));
	}
	else
	{
		fp = fopen_or_die(Fnameptr->SGCbankfilename, "r", "Loading SGCbank ", verbose);

		for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
		fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->SGCz + i + j*Parptr->xsz);
			if (AreEqual(Arrptr->SGCz[i + j*Parptr->xsz], no_data_value))
				Arrptr->SGCz[i + j*Parptr->xsz] = Arrptr->DEM[i + j*Parptr->xsz]; // In the case of no_data_value set the bank height to the DEM
		}
		fclose(fp);
		if (verbose == ON) printf("Done.\n");

		fp = fopen_or_die(Fnameptr->SGCwidthfilename, "r", "Loading SGCwidth ", verbose);

		for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
		fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", Arrptr->SGCwidth + i + j*Parptr->xsz);
			if (AreEqual( Arrptr->SGCwidth[i + j*Parptr->xsz], no_data_value))
				Arrptr->SGCwidth[i + j*Parptr->xsz] = C(0.0); // In the case of no_data_value set the width to zero
		}
		fclose(fp);

	}
	if (verbose == ON) printf("Done.\n\n");

	// This loads the distributed SGC group information
	if (Statesptr->SGCchangroup == ON)
	{
		fp = fopen_or_die(Fnameptr->SGCchangroupfilename, "r", "Loading SGC channel group", verbose);

		for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
		fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", &tmp); // read a NUMERIC_TYPE in case someone doesn't put intergers in the ascii file.
			Arrptr->SGCgroup[i + j*Parptr->xsz] = (int)tmp;
			if (AreEqual(Arrptr->SGCgroup[i + j*Parptr->xsz], no_data_value) || Arrptr->SGCgroup[i + j*Parptr->xsz] < 0) 
				Arrptr->SGCgroup[i + j*Parptr->xsz] = 0; // In the case of no_data_value set the chan group to zero
		}
		fclose(fp);
		if (verbose == ON) printf("Done.\n\n");
	}
	if (Statesptr->ChanMaskRead == ON) // loads SGC channel mask if present
	{ 
		// creat memory space for channel mask
		memset(Arrptr->ChanMask, 0, sizeof(int)*(Parptr->xsz)*(Parptr->ysz));
		// open file
		fp = fopen_or_die(Fnameptr->ChanMaskfilename, "r", "Creating and Loading SGC Channel Mask ", verbose);
		// read header
		for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
		fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);
		// read data
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", &tmp);
			if (AreEqual(tmp, no_data_value)) Arrptr->ChanMask[i + j*Parptr->xsz] = 0; // In the case of no_data_value set the chan mask to zero
			else if (tmp > C(0.0)) Arrptr->ChanMask[i + j*Parptr->xsz] = 1; // a positive number is a mask cell 1
			else Arrptr->ChanMask[i + j*Parptr->xsz] = 0; // anything else is a no mask cell 0
		}
		fclose(fp);
		if (verbose == ON) printf("Done.\n");
	}
	


	return;
}
//-----------------------------------------------------------------------------
// LOAD GAUGE DATA
void LoadGauges(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Stage *Locptr, const int verbose)
{
	//Added by Jeff Neal, 22 Jul 2011
	//Provides functionality to output regular section measurements of discharge

	int i;
	char dum[10];
	FILE *fp;
	fp = fopen_or_die(Fnameptr->gaugefilename, "r", "Gauge section information", verbose);
	Statesptr->gsection = ON;

	fscanf(fp, "%d", &Locptr->Ngauges);
	fgetc(fp); // Retrieve closing EOL

	Locptr->gauge_loc_x = memory_allocate_zero_numeric_legacy(Locptr->Ngauges);
	Locptr->gauge_loc_y = memory_allocate_zero_numeric_legacy(Locptr->Ngauges);
	Locptr->gauge_grid_x = new int[Locptr->Ngauges]();
	Locptr->gauge_grid_y = new int[Locptr->Ngauges]();
	Locptr->gauge_dir = new EDirection[Locptr->Ngauges]();
	Locptr->gauge_dist = memory_allocate_zero_numeric_legacy(Locptr->Ngauges);
	Locptr->gauge_cells = new int[Locptr->Ngauges]();

	//scan x,y locations from file
	for (i = 0; i < Locptr->Ngauges; i++)
	{
		fscanf(fp, "%" NUM_FMT"", &Locptr->gauge_loc_x[i]);
		fscanf(fp, "%" NUM_FMT"", &Locptr->gauge_loc_y[i]);
		fscanf(fp, "%s", dum);
		if (!STRCMPi(dum, "N")) Locptr->gauge_dir[i] = North;
		if (!STRCMPi(dum, "E")) Locptr->gauge_dir[i] = East;
		if (!STRCMPi(dum, "S")) Locptr->gauge_dir[i] = South;
		if (!STRCMPi(dum, "W")) Locptr->gauge_dir[i] = West;

		fscanf(fp, "%" NUM_FMT"", &Locptr->gauge_dist[i]);
		Locptr->gauge_cells[i] = int(ceil(Locptr->gauge_dist[i] / Parptr->dx)); // work out number of cells
	}
	for (i = 0; i < Locptr->Ngauges; i++)
	{
		// convert coordinates to cells
		Locptr->gauge_grid_x[i] = int(floor((Locptr->gauge_loc_x[i] - Parptr->blx) / Parptr->dx));
		Locptr->gauge_grid_y[i] = Parptr->ysz - 1 - (int(floor((Locptr->gauge_loc_y[i] - Parptr->bly) / Parptr->dy)));

		// check for off-image values and set to domain edge
		if (Locptr->gauge_grid_x[i] < 0) Locptr->gauge_grid_x[i] = 0;
		if (Locptr->gauge_grid_x[i] >= Parptr->xsz) Locptr->gauge_grid_x[i] = Parptr->xsz - 1;
		if (Locptr->gauge_grid_y[i] < 0) Locptr->gauge_grid_y[i] = 0;
		if (Locptr->gauge_grid_y[i] >= Parptr->ysz) Locptr->gauge_grid_y[i] = Parptr->ysz - 1;

		// adjust distances if these will go off the domain - check these (Toby.D Checked and fixed!)
		if ((Locptr->gauge_dir[i] == North || Locptr->gauge_dir[i] == South) && 
			Locptr->gauge_grid_y[i] + Locptr->gauge_cells[i] > Parptr->ysz - 1)
			Locptr->gauge_cells[i] = Parptr->ysz - 1 - Locptr->gauge_grid_y[i];
		else if ((Locptr->gauge_dir[i] == East || Locptr->gauge_dir[i] == West)
			&& Locptr->gauge_grid_x[i] + Locptr->gauge_cells[i] > Parptr->xsz - 1)
			Locptr->gauge_cells[i] = Parptr->xsz - 1 - Locptr->gauge_grid_x[i];

		// increment x and y positions of gauges - to measure the correct cell face
		if (Locptr->gauge_dir[i] == East)
			Locptr->gauge_grid_x[i]++; // start on east face
		if (Locptr->gauge_dir[i] == South)
			Locptr->gauge_grid_y[i]++; // start on south face
	}

	if (verbose == ON) printf("Done.\n\n");

	fclose(fp);
	return;
}

// LOAD SGC PARAMETER DATA
void LoadSGCChanPrams(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, SGCprams *SGCptr, const int verbose)
{
	//Added by Jeff Neal, 06 Aug 2012
	//Provides functionality to import difstrbuted sub-grid channel parameters

	int i, j, tmp, buff_size = 800;
	char buff[800];
	FILE *fp;
	fp = fopen_or_die(Fnameptr->SGCchanpramsfilename, "r", "Loading SGC channel parameter information", verbose);

	//read line
	for (j = 0; j < buff_size; j++)
	{
		buff[j] = fgetc(fp);
		if (buff[j] == '\n' || buff[j] == EOF) break;
	}
	buff[j] = '\0';									// Finish off string
	sscanf(buff, "%i", &SGCptr->NSGCprams);

	// create new variables
	SGCptr->SGCchantype = new int[SGCptr->NSGCprams]();
	SGCptr->SGCp = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);
	SGCptr->SGCr = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);
	SGCptr->SGCs = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);
	SGCptr->SGCn = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);
	SGCptr->SGCm = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);
	SGCptr->SGCa = memory_allocate_zero_numeric_legacy(SGCptr->NSGCprams);

	if (verbose == ON) printf("Num   Type  p     r     sl    n     m     a    \n");
	//scan x,y locations from file
	for (i = 0; i < SGCptr->NSGCprams; i++)
	{
		// initalise with defaults
		SGCptr->SGCchantype[i] = 1;
		SGCptr->SGCp[i] = Parptr->SGC_p;
		SGCptr->SGCr[i] = Parptr->SGC_r;
		SGCptr->SGCs[i] = Parptr->SGC_s;
		SGCptr->SGCn[i] = Parptr->SGC_n;
		SGCptr->SGCm[i] = Parptr->SGC_m;
		SGCptr->SGCa[i] = Parptr->SGC_a;
		// load buffer until EOL
		for (j = 0; j < buff_size; j++)
		{
			buff[j] = fgetc(fp);
			if (buff[j] == '\n' || buff[j] == EOF) break;
		}
		buff[j] = '\0';									// Finish off string
		sscanf(buff, "%i%i%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"", &tmp, &SGCptr->SGCchantype[i], &SGCptr->SGCp[i], &SGCptr->SGCr[i], &SGCptr->SGCs[i], &SGCptr->SGCn[i], &SGCptr->SGCm[i], &SGCptr->SGCa[i]);

		if (verbose == ON && SGCptr->SGCchantype[i] == 2 && SGCptr->SGCs[i] > 20)  printf("Warning channel shape exponent above recomended value");
		if (verbose == ON && SGCptr->SGCchantype[i] == 2 && SGCptr->SGCs[i] < C(1.3)) printf("Warning channel shape exponent below recomended value");
		if (verbose == ON && SGCptr->SGCs[i] < C(0.0)) printf("ERROR SGC meander coefficient is too low!");
		if (verbose == ON) printf("%i     %i     %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT"\n", i, SGCptr->SGCchantype[i], SGCptr->SGCp[i], SGCptr->SGCr[i], SGCptr->SGCs[i], SGCptr->SGCn[i], SGCptr->SGCm[i], SGCptr->SGCa[i]);
	}
	if (verbose == ON) printf("Done.\n\n");

	fclose(fp);
	return;
}

// LOAD DAM PARAMETER DATA //FEOL
void LoadDamPrams(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, DamData *Damptr, const int verbose)
{
	//Added by Fiachra O'Loughlin, 20 July 2016
	//Provides functionality to import Dam parameters

	int i, j, tmp, buff_size = 800;
	char buff[800];
	FILE *fp;
	fp = fopen_or_die(Fnameptr->Damparfilename, "r", "Loading Dam parameter information", verbose);

	//read line
	for (j = 0; j < buff_size; j++)
	{
		buff[j] = fgetc(fp);
		if (buff[j] == '\n' || buff[j] == EOF) break;
	}
	buff[j] = '\0';									// Finish off string
	sscanf(buff, "%i", &Damptr->NumDams);

	// create new variables
	Damptr->Volmax = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamArea = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->InitialHeight = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamHeight = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->SpillWidth = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->Spill_Cd = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->SpillHeight = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamOperationCode = new int[Damptr->NumDams]();
	Damptr->DamMeanQ = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamOperationQ = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamVin = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamVol = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->DamTotalQ = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->SpillQ = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->AnnualRelease = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->OP7_Kappa = memory_allocate_zero_numeric_legacy(Damptr->NumDams);
	Damptr->OutputCellX = memory_allocate_zero_numeric_legacy(Damptr->NumDams);//new int[Damptr->NumDams]();
	Damptr->OutputCellY = memory_allocate_zero_numeric_legacy(Damptr->NumDams);//;
	Damptr->DamMaxH = C(0.0);
	Damptr->DamYear = new int[Damptr->NumDams]();
	
	// Damptr->DynamicEdge = new DamEdge *[Damptr->NumDams];

	

	if (verbose == ON) printf("Num   Vol	Area	Initial_H	Dam_H     Spill_Width     Spill_Cd    Spill_Height     Dam_Op_Q     Output_X	Output_Y    \n");
	//scan x,y locations from file
	for (i = 0; i < Damptr->NumDams; i++)
	{
		// Set Size of DynamicEgde to NULL;
		// Damptr->DynamicEdge[i] = NULL;
		// load buffer until EOL
		for (j = 0; j < buff_size; j++)
		{
			buff[j] = fgetc(fp);
			if (buff[j] == '\n' || buff[j] == EOF) break;
		}
		buff[j] = '\0';									// Finish off string
		sscanf(buff, "%i%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"%i%" NUM_FMT"%" NUM_FMT"%" NUM_FMT"", &tmp, &Damptr->Volmax[i], &Damptr->DamArea[i], &Damptr->InitialHeight[i], &Damptr->DamHeight[i], &Damptr->SpillWidth[i], &Damptr->Spill_Cd[i], &Damptr->SpillHeight[i], &Damptr->DamOperationCode[i], &Damptr->DamMeanQ[i], &Damptr->OutputCellX[i], &Damptr->OutputCellY[i]);
				
		if (verbose == ON) printf("%i     %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT" %.3" NUM_FMT"  %" NUM_FMT"  %" NUM_FMT"\n", i, Damptr->Volmax[i], Damptr->DamArea[i], Damptr->InitialHeight[i], Damptr->DamHeight[i], Damptr->SpillWidth[i], Damptr->Spill_Cd[i], Damptr->SpillHeight[i], Damptr->DamMeanQ[i], Damptr->OutputCellX[i], Damptr->OutputCellY[i]);
		Damptr->OutputCellX[i] = (floor((Damptr->OutputCellX[i] - Parptr->blx) / Parptr->dx));
		Damptr->OutputCellY[i] = Parptr->ysz - 1 - ((floor((Damptr->OutputCellY[i] - Parptr->bly) / Parptr->dy)));
	}
	if (verbose == ON) printf("Done.\n\n");

	fclose(fp);
	
	return;
}
void LoadDamMask(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, DamData *Damptr, const int verbose)  //FEOL
{
	NUMERIC_TYPE no_data_value = -9999;
	int i, j;
	int n; //Check
	int num_cols; int num_rows; NUMERIC_TYPE  xllcorner; NUMERIC_TYPE  yllcorner; NUMERIC_TYPE  cell_size;

		
	read_file(Fnameptr->DamMaskfilename, no_data_value, &num_cols, &num_rows, &Arrptr->DamMask, &xllcorner, &yllcorner, &cell_size);

	// Changes DEM to DEM_NO_DATA where mask is negative (no flow cells). No flow cells are ignored in 2-D solver.
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		if (Arrptr->DamMask[i + j*Parptr->xsz] < C(0.0)) 
		{
			Arrptr->DEM[i + j*Parptr->xsz] = no_data_value;
			if (AreEqual(Arrptr->DEM[i + j*Parptr->xsz], no_data_value))
			Arrptr->DEM[i + j*Parptr->xsz] = DEM_NO_DATA;
		}
	}
	Damptr->Edgenos = new int[Damptr->NumDams]();
	for (j = 0; j <Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++) for (n = 0; n < Damptr->NumDams; n++)
	{

		if (Arrptr->DamMask[i + j*Parptr->xsz]== (n+1))
			Damptr->Edgenos[n]++;		
	}
	return;
}

