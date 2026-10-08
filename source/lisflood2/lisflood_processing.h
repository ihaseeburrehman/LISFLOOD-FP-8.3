#pragma once

#include "../lisflood.h"

void Fast_MainStart(Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, SGCprams *SGCptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, DamData *Damptr, int verbose);