#pragma once
namespace lis
{

typedef struct SolverParams
{
	NUMERIC_TYPE nodata_elevation;
	NUMERIC_TYPE cfl;
	NUMERIC_TYPE max_dt;
	NUMERIC_TYPE DepthThresh;
	NUMERIC_TYPE SpeedThresh;
	NUMERIC_TYPE KrivodonovaThresh;
	NUMERIC_TYPE MaxHflow; 
	int routing; 
	NUMERIC_TYPE RouteSfThresh; 
	NUMERIC_TYPE theta; 
	int voutput; 
	NUMERIC_TYPE InitTstep;

	SolverParams() {};

	SolverParams
	(
		::Pars& pars,
		::Solver& solver,
		::States& states 
	)
	:
	nodata_elevation(pars.nodata_elevation),
	cfl(solver.cfl),
	max_dt(solver.InitTstep),
	DepthThresh(solver.DepthThresh),
	SpeedThresh(solver.SpeedThresh),
	KrivodonovaThresh(solver.krivodonova_threshold),
	MaxHflow(solver.MaxHflow), 
		routing(states.routing), 
	RouteSfThresh(pars.RouteSfThresh), 
	theta(solver.theta), 
	voutput(states.voutput), 
	InitTstep(solver.InitTstep)
	{}
} SolverParams;

typedef struct PhysicalParams
{
	NUMERIC_TYPE g;
	NUMERIC_TYPE manning;

	PhysicalParams() {};

	PhysicalParams
	(
		::Pars& pars,
		::Solver& solver
	)
	:
	g(solver.g),
	manning(pars.FPn)
	{}
} PhysicalParams;



}
