//===================================================================
// File:        main.cc
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Program entry point.
// Project:     A Distributed Multidimensional Knapsack Problem Solver
//              using Island Based Particle Swarm Optimization
//              (MKPDIPSO).
//-------------------------------------------------------------------
// Copyright (C) 2013 Drahoslav Zan
//
// This file is part of MKPDIPSO.
//
// MKPDIPSO is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MKPDIPSO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MKPDIPSO. If not, see <http://www.gnu.org/licenses/>.
//===================================================================
// vim: set nowrap sw=2 ts=2


#include "generic.h"
#include "knapsack.h"
#include "data.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>

#include <getopt.h>
#include <unistd.h>
#include <papi.h>

#include <mpi.h>


using namespace std;


static const char * progName;


ostream & operator <<(ostream & os, const MKPSolver::sol_t & solution)
{
	if(!solution.size())
		return os;

	os << (unsigned)solution[0];

	for(size_t i = 1; i < solution.size(); ++i)
	{
		if(!(i % 64)) os << endl;
		else if(!(i % 8)) os << " ";

		os << (unsigned)solution[i];
	}

	return os;
}


inline size_t strton(const char * a, char ** b)
{
	return strtol(a, b, 10);
}

template <typename T, T (*F)(const char *, char **)>
static void split(std::vector<T> &v, char * str)
{
	char * end = str;

	if(!*end) return;

	for(size_t n = 0; n < v.size(); ++n)
	{
 		T i = F(end, &end);

		if(!i) break;

		if(!*end)
		{
			for(; n < v.size(); ++n)
				v[n] = i;
			break;
		}

		v[n] = i;
	}
}

template <typename T>
static std::ostream & operator <<(ostream &os, std::vector<T> v)
{
	if(!v.size())
		return os;

	os << "0:" <<  v[0];

	for(size_t i = 1; i < v.size(); ++i)
		os << ", " << i << ":" << v[i];

	return os;
}


static void usage(FILE * stream)
{
	fprintf(stream,
			"Usage: %s [OPTION]\n"
			"Try `%s --help' for more information.\n"
			, progName, progName);
}

static void help(FILE * stream)
{
	fprintf(stream,
			"Usage: %s [OPTION] <FILE>\n"
			"\n\n"
#ifdef _OPENMP
			"OpenMP:\n"
			"    --threads     <n>   Number of threads\n"
			"\n"
#endif
			"MPI:\n"
			"    --exchange    <n>   Interval for exchange of the global best.\n"
			"                        between processes\n"
			"    --stats       <n>   Interval for statistics exchange\n"
			"\n"
			"Particle Swarm Optimization:\n"
			"    --epochs      <n>   Number of iterations\n"
			"    --stop        <n>   Stopping criterion in number of exchange intervals\n"
			"    --particles   <n>   Number of particles\n"
			"    --inertia     <r>   Inertia weight\n"
			"    --cognitive   <r>   Cognitive coefficient\n"
			"    --social      <r>   Social coefficient\n"
			"\n"
			"Multidimensional Knapsack Problem:\n"
			"    --penalty     <r>   Penalty coefficient\n"
			"    --velocity    <r>   Maximum value for velocity clamping\n"
			"\n"
			"Others:\n"
			"    --verbose           Verbose output\n"
			"    --help              Show this help and exit\n"
			"\n\n"
			"Argument <n> means integer, <r> means real.\n"
			"\n"
			"Multidimensional Knapsack Problem solver using distributed multi-swarm\n"
			"Particle Swarm Optimization.\n"
			, progName);
}


int main(int argc, char **argv)
{
	progName = argv[0];

	int id;
	int mpisz;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisz);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

#ifdef _OPENMP
	size_t oThreads   = OMP_PROCS;
#endif
	size_t oExchange          = PSO_DEFAULT_EXCHANGE_INTERVAL;
	size_t oStats             = PSO_DEFAULT_STATS_INTERVAL;
	size_t oEpochs            = PSO_DEFAULT_EPOCHS;
	size_t oStop              = PSO_DEFAULT_STOP_CRITERION;
	vector<size_t> oParticles(mpisz, PSO_DEFAULT_PARTICLES);
	vector<float>  oInertia(mpisz, PSO_DEFAULT_INERTIA_WEIGHT);
	vector<float>  oCognitive(mpisz, PSO_DEFAULT_COGNITIVE);
	vector<float>  oSocial(mpisz, PSO_DEFAULT_SOCIAL);
	vector<float>  oPenalty(mpisz, MKP_DEFAULT_PENALTY);
	vector<float>  oVelocity(mpisz, MKP_DEFAULT_VELOCITY_CLAMP);
	int    oVerbose           = 0;
	
	for(opterr = 0;;)
	{
		static struct option longOptions[] = {
#ifdef _OPENMP
			{ "threads"   , required_argument , NULL        , 't' },
#endif
			{ "exchange"  , required_argument , NULL        , 'e' },
			{ "stats"     , required_argument , NULL        , 's' },
			{ "epochs"    , required_argument , NULL        , 'm' },
			{ "stop"      , required_argument , NULL        , 'o' },
			{ "particles" , required_argument , NULL        , 'n' },
			{ "inertia"   , required_argument , NULL        , '1' },
			{ "cognitive" , required_argument , NULL        , '2' },
			{ "social"    , required_argument , NULL        , '3' },
			{ "penalty"   , required_argument , NULL        , 'p' },
			{ "velocity"  , required_argument , NULL        , 'v' },
			{ "verbose"   , no_argument       , &oVerbose   ,  1  },
			{ "help"      , no_argument       , NULL        , 'h'	},
			{ NULL        , 0                 , NULL        ,  0  },
		};

		int c = getopt_long(argc, argv, IF_OMP("t:") "e:s:m:o:n:1:2:3:p:v:h",
				longOptions, NULL);

		if (c == -1) break;

		int   iarg;

		switch(c)
		{
			case 0:
				break;

#ifdef _OPENMP
			case 't':
				iarg = atoi(optarg);
				if(iarg >= 0) oThreads = iarg;
				break;
#endif
			case 'e':
				iarg = atoi(optarg);
				if(iarg >= 0) oExchange = iarg;
				break;
			case 's':
				iarg = atoi(optarg);
				if(iarg >= 0) oStats = iarg;
				break;
			case 'm':
				iarg = atoi(optarg);
				if(iarg >= 0) oEpochs = iarg;
				break;
			case 'o':
				iarg = atoi(optarg);
				if(iarg >= 0) oStop = iarg;
				break;
			case 'n':
				split<size_t, strton>(oParticles, optarg);
				break;
			case '1':
				split<float, strtof>(oInertia, optarg);
				break;
			case '2':
				split<float, strtof>(oCognitive, optarg);
				break;
			case '3':
				split<float, strtof>(oSocial, optarg);
				break;
			case 'p':
				split<float, strtof>(oPenalty, optarg);
				break;
			case 'v':
				split<float, strtof>(oVelocity, optarg);
				break;

			case 'h':
				if(id == MPIROOT) help(stdout);
				MPI_Finalize();
				return 0;
			case '?':
				if(id == MPIROOT)
				{
					cerr << "ERROR: getopt: Unrecognized option provided" << endl;
					usage(stderr);
				}
				MPI_Finalize();
				return 0;
			default:
				if(id == MPIROOT) cerr << "ERROR: getopt: " << c << endl;
				MPI_Finalize();
				return 0;
		}
	}

	if(optind != argc - 1)
	{
		cerr << "ERROR: Missing file or more than one provided!" << endl;
		usage(stderr);
		MPI_Finalize();
		return 0;
	}

#ifdef _OPENMP
	omp_set_num_threads(oThreads);
#endif

	for(size_t i = 0; i < oParticles.size(); ++i)
		oParticles[0] = MKPSolver::adjustParticleCount(oParticles[0]);

	try
	{
		MKPData data(argv[optind]);

		MKPSolver mkp(data, oParticles[id]);
		mkp.setInertiaWeight(oInertia[id]);
		mkp.setLearningFactors(oCognitive[id], oSocial[id]);
		mkp.setVelocityMax(oVelocity[id]);
		mkp.setPenalty(oPenalty[id]);

		// PERFORMANCE METRICS

		long_long vtus;

		vtus = PAPI_get_virt_usec();

		const MKPSolver::sol_t & solution =
			mkp.solve(oEpochs, oStop, oExchange, oStats);

		vtus = PAPI_get_virt_usec() - vtus;

		double sec = vtus / 1.e6;

		// VERBOSE

		if(id == MPIROOT)
		{
			if(oVerbose)
			{
				cout << "PARAMS" << endl;
				cout << "--------------------------------------" << endl;
				cout << "Number of processes  :\t" << mpisz << endl;
#ifdef _OPENMP
				cout << "Number of threads    :\t" << oThreads << endl;
#endif
				cout << "Epochs               :\t" << MKPSolver::adjustEpochs(oEpochs) << endl;
				cout << "Stopping criterion   :\t" << oStop << endl;
				cout << "Exchange interval    :\t" << oExchange << endl;
				cout << "Statistics interval  :\t" << oStats << endl;
				cout << "Particles            :\t" << oParticles << endl;
				cout << "Inertia weight       :\t" << oInertia << endl;
				cout << "Cognitive coefficient:\t" << oCognitive << endl;
				cout << "Social coefficient   :\t" << oSocial << endl;
				cout << "Penalty              :\t" << oPenalty << endl;
				cout << "Maximum velocity     :\t" << oVelocity << endl;
				cout << endl;
	
				const MKPSolver::Stats & stats = mkp.statistics();
	
				cout << "STATISTICS" << endl;
				cout << "--------------------------------------" << endl;
				cout << "Iters  :\t" << stats.iterations << endl;
				cout << "Fitness." << endl;
				cout << "  best :\t" << stats.fitness.best << " / " << data.optimum() << endl;
				cout << "  worst:\t" << stats.fitness.worst << endl;
				cout << "  avg  :\t" << stats.fitness.avg << endl;
				cout << "Weight." << endl;
	
				assert(data.knapsacks() == stats.weight.size());
				for(size_t i = 0; i < data.knapsacks(); ++i)
					cout << "  " << stats.weight[i] << " / " << data.capacities()[i] << endl;
	
				cout << endl;
	
				cout << "PERFORMANCE" << endl;
				cout << "--------------------------------------" << endl;
				cout << "Virtual:\t" << sec << "s" << endl;
	
				cout << endl;
				cout << "SOLUTION (" << argv[optind] << ")" << endl;
				cout << "--------------------------------------" << endl;
			}
	
			// PRESENTATION
			//cout << '[' << solution << ']' << endl;
			cout << solution << endl;
		}
	}
	catch(const std::string & error)
	{
		if(id == MPIROOT)
			cerr << "ERROR: " << error << endl;
	}

	MPI_Finalize();

	return 0;
}
