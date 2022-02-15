/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of "fullrecoursekidney" a collection of              */
/*    methods to solve the full recourse robust kidney exchange problem.     */
/*                                                                           */
/*    Copyright (C) 2020-2021  Danny Blom, Christopher Hojny, Bart Smeulders */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.cpp
 * @brief  main file for solving the full recourse kidney exchange problem
 * @author Christopher Hojny
 */

#include "graph.h"
#include "parseOptions.h"
#include "problem_master_kidneyexchange_picef.h"
#include "problem_master_kidneyexchange.h"
#include "problem_kidneyexchange.h"
#include "kidneyParams.h"
#include "kidneyPlugins.h"
#include "read_graph.h"
#include "solveMasterProblem.h"
#include "solveMasterPICEFProblem.h"
#include "find_graphstructures.h"
#include "auxiliaryStructures.h"
#include "typedefs.h"

#include <scip/scip.h>

#include <iostream>
#include <limits>

/**@brief Solve partitioning problem given in file @a filename */
static
SCIP_RETCODE solveKidneyExchangeProblem(
   std::string           filename,           //!< name of instance file
   int                   adversarybound,     //!< bound on adversary attack
   const char*           settings,           //!< Possible name of setting file
   double                timeLimit,          //!< time limit
   double                memLimit,           //!< memory limit
   SCIP_Longint          nodeLimit,          //!< node limit
   int                   displayFreq,        //!< display frequency
   SCIP_Bool             verbose             //!< whether we print SCIP's logs
   )
{  /*lint --e{429}*/

   // initialize SCIP
   SCIP* masterscip;
   int method;

   SCIP_CALL( SCIPcreate(&masterscip) );

   // output SCIP banner
#ifdef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIPprintVersion(masterscip, NULL);
#else
   SCIPprintVersion(NULL);
#endif
   std::cout << "\n" << std::endl;   // (to force flush)

   // load basic plugins
   SCIP_CALL( includeKidneyPlugins(masterscip) );

   // read parameters
   SCIP_CALL( setSCIPParameters(masterscip) );
   SCIP_CALL( addKidneyParameters(masterscip) );

   if ( settings != 0 )
   {
      if ( ! SCIPfileExists(settings) )
      {
         SCIPerrorMessage("Setting file <%s> does not exist.\n\n", settings);
      }
      else
      {
         SCIPinfoMessage(masterscip, 0, "Reading parameters from <%s>.\n\n", settings);
         SCIP_CALL( SCIPreadParams(masterscip, settings) );
      }
   }

   // get the solution method we want to use for solving the problem
   SCIP_CALL( SCIPgetIntParam(masterscip, "kidney/method", &method) );

   // output changed parameters
   SCIPinfoMessage(masterscip, 0, "Changed settings:\n");
   SCIP_CALL( SCIPwriteParams(masterscip, 0, FALSE, TRUE) );
   SCIPinfoMessage(masterscip, 0, "\n");

   // set limits
   if ( timeLimit < 1e20 )
   {
      SCIPinfoMessage(masterscip, 0, "Setting time limit to %g.\n", timeLimit);
      SCIP_CALL( SCIPsetRealParam(masterscip, "limits/time", timeLimit) );
   }
   if ( memLimit < 1e20 )
   {
      SCIPinfoMessage(masterscip, 0, "Setting memory limit to %g.\n", memLimit);
      SCIP_CALL( SCIPsetRealParam(masterscip, "limits/memory", memLimit) );
   }
   if ( nodeLimit < SCIP_LONGINT_MAX )
   {
      SCIPinfoMessage(masterscip, 0, "Setting node limit to %ld.\n", nodeLimit);
      SCIP_CALL( SCIPsetLongintParam(masterscip, "limits/nodes", nodeLimit) );
   }
   if ( displayFreq < INT_MAX )
   {
      SCIPinfoMessage(masterscip, 0, "Setting display frequency to %d.\n", displayFreq);
      SCIP_CALL( SCIPsetIntParam(masterscip, "display/freq", displayFreq) );
   }
   SCIPinfoMessage(masterscip, 0, "\n");

   // create problem
   if ( method != METHOD_BENDERS_PICEF )
   {
      Graph* graph;
      Cycles* cycles;
      Chains* chains;
      Triplets* triplets;

      clock_t begin = clock();

      cout << "Reading the graph" << endl;
      SCIP_CALL( readGraph(masterscip, filename, &graph) );
      cout << "Finding the cycles" << endl;
      SCIP_CALL( findcycles(masterscip, &cycles, graph) );
      cout << "Finding the chains" << endl;
      SCIP_CALL( findchains(masterscip, &chains, graph) );
      cout << "Strengthening the node bounds" << endl;
      SCIP_CALL( strengthenNodeBounds(masterscip, graph, cycles, chains) );

      //does nothing if usetwothirdcliques is set to FALSE (triplets = NULL)
      cout << "Find the triplets" << endl;
      SCIP_CALL( findtriplets(masterscip, cycles, chains, &triplets, graph) );

      SCIP_CALL( generateNode2Cycles(masterscip, graph, cycles) );
      SCIP_CALL( generateNode2Chains(masterscip, graph, chains) );

      clock_t end = clock();
      double readingtime = (double) (end - begin)/CLOCKS_PER_SEC;
      cout << "Time needed to read graph and detecting cycles/chains/triplets:\t" << readingtime << endl;

      cout << "Create the model" << endl;
      begin = clock();
      SCIP_CALL( SCIPcreateMasterModel(masterscip, graph, cycles, chains, adversarybound) );
      end = clock();

      double creatingtime = (double) (end - begin)/CLOCKS_PER_SEC;
      cout << "Time to create master model:\t" << creatingtime << endl;

      // tell SCIP that optimal objective value is always integral
      SCIP_CALL( SCIPsetObjIntegral(masterscip) );

      cout << "Solve the model" << endl;
      SCIP_CALL( solveMasterProblem(masterscip, graph, cycles, chains, triplets, adversarybound,
            timeLimit, readingtime + creatingtime, settings, verbose) );

      // output statistics
      SCIPinfoMessage(masterscip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(masterscip, NULL) );

      // free data
      cout << "Free the model" << endl;
      SCIP_CALL( SCIPfreeMasterModel(masterscip, graph, cycles, chains, triplets) );
   }
   else
   {
      Graph* graph;
      Cycles* cycles;
      PositionedArcs* posarcs;

      clock_t begin = clock();
      cout << "Reading the graph" << endl;
      SCIP_CALL( readGraph(masterscip, filename, &graph) );
      cout << "Finding the cycles" << endl;
      SCIP_CALL( findcycles(masterscip, &cycles, graph) );

      cout << "Finding the position indexed arcs" << endl;
      SCIP_CALL( generatePositionedArcs(masterscip, graph, &posarcs) );

      SCIP_CALL( generateNode2Cycles(masterscip, graph, cycles) );
      graph->node2chains = NULL;
      graph->node2chainsbegin = NULL;
      clock_t end = clock();
      double readingtime = (double) (end - begin)/CLOCKS_PER_SEC;

      cout << "Create the position indexed mastermodel" << endl;
      begin = clock();
      SCIP_CALL( SCIPcreateMasterPICEFModel(masterscip, graph, cycles, posarcs, adversarybound) );
      end = clock();
      double creatingtime = (double) (end - begin)/CLOCKS_PER_SEC;


      // tell SCIP that optimal objective value is always integral
      SCIP_CALL( SCIPsetObjIntegral(masterscip) );

      cout << "Solve the model" << endl;
      SCIP_CALL( solveMasterPICEFProblem(masterscip, graph, cycles, posarcs, adversarybound,
            timeLimit, readingtime + creatingtime, settings, verbose) );

      // output statistics
      SCIPinfoMessage(masterscip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(masterscip, NULL) );

      // free data
      cout << "Free the model" << endl;
      SCIP_CALL( SCIPfreeMasterPICEFModel(masterscip, graph, cycles, posarcs) );
   }
   SCIP_CALL( SCIPfreeTransform(masterscip) );
   SCIP_CALL( SCIPfree(&masterscip) );

   BMScheckEmptyMemory();
   return SCIP_OKAY;
}


/** main function for solving the partitioning problem */
int main(int argc, const char** argv)
{
   // check parameters
   std::vector<std::string> mainArgs;
   std::vector<std::string> otherArgs;

   std::vector<std::string> knownOptions{"s", "t", "m", "n", "d"};

   if ( ! readOptions(argc, argv, 2, 2, mainArgs, otherArgs, knownOptions) )
   {
      std::cerr << "usage: " << argv[0] << " <instance file> <number of attacks> [-s <settings>] [-t <time limit>] [-m <memory limit>] [-n <node limit>] ";
      std::cerr << "[-d <disp. freq>]" << std::endl;
      exit(1);
   }

   // get optional arguments
   const char* settings = 0;
   double timeLimit = 1e20;
   double memLimit = 1e20;
   SCIP_Longint nodeLimit = SCIP_LONGINT_MAX;
   int dispFreq = INT_MAX;
   SCIP_Bool verbose = FALSE;
   for (long unsigned int j = 0; j < otherArgs.size(); ++j)
   {
      if ( otherArgs[j] == "s" )
         settings = otherArgs[++j].c_str();
      else if ( otherArgs[j] == "t" )
         timeLimit = atof(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "m" )
         memLimit = atof(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "n" )
         nodeLimit = atol(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "d" )
         dispFreq = atoi(otherArgs[++j].c_str());
      else
         dispFreq = 1000;
   }

   // check for adversary attack bound
   int adversarybound = atoi(mainArgs[1].c_str());

   // run kidney exchange code
   SCIP_RETCODE retcode;
   retcode = solveKidneyExchangeProblem(mainArgs[0], adversarybound, settings, timeLimit, memLimit, nodeLimit, dispFreq, verbose);

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }
   return 0;
}
