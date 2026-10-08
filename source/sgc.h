#pragma once

#include "lisflood.h"

NUMERIC_TYPE CalcSGC_UpV(int, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE);
void CalcSGC_A(int, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE *, NUMERIC_TYPE *, const SGCprams *);
NUMERIC_TYPE CalcSGC_R(int, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, const SGCprams *);
void SGC_hotstart(States *, Pars *, Solver *, Arrays *);
void CalcSGCz(Fnames *, States *, Pars *, Arrays *, SGCprams *, const int verbose);

#ifdef RESULT_CHECK
// Sub gid channel prototypes
void SGC_FloodplainQ(const States *, const Pars *, const Solver *, Arrays *, const SGCprams *);

NUMERIC_TYPE CalcFPQxSGC(int i, int j, const States *, const Pars *, const Solver *, const Arrays *, const SGCprams *);
NUMERIC_TYPE CalcFPQySGC(int i, int j, const States *, const Pars *, const Solver *, const Arrays *, const SGCprams *);
NUMERIC_TYPE SGC_UpdateH(States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, SGCprams *);
void SGC_BCs(const States *Statesptr, const Pars *Parptr, const Solver *Solverptr, const BoundCs *BCptr, const ChannelSegmentType *ChannelSegments, Arrays *Arrptr, const SGCprams *SGCptr);

void SGC_Evaporation(Pars *, Solver *, Arrays *, SGCprams *);
void SGC_Rainfall(Pars *, Solver *, Arrays *); // CCS May 2013
void SGC_Routing(States *, Pars *, Solver *, Arrays *);// CCS May 2013
NUMERIC_TYPE CalcSGC_UpH(int, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE);

void SGC_wp_prams(SGCprams *);
void CalcSGC_pointFREE(NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, NUMERIC_TYPE, int, int, NUMERIC_TYPE *, NUMERIC_TYPE *, NUMERIC_TYPE *, const SGCprams *);

#endif