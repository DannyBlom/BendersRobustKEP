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

/**@file   probdata_benders_subproblem.c
 * @brief  Problem data for the subproblem of the Benders type approach (KEP restricted to subgraph)
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * A list of all interface methods can be found in probdata_benders_subproblem.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_benders_subproblem.h"

#include "scip/scip.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the kidney exchange problem, all variables which are created, and all
 * constraints.
 */

struct SCIP_ProbData
{
   Graph*                graph;              /**< pointer to underlying graph */
   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   Chains*               chains;             /**< pointer to chain structure of underlying graph */
   int                   nnodes;             /**< Number of nodes in the graph */

   SCIP_VAR**            cyclevars;          /**< Variable array for all cycles */
   SCIP_VAR**            chainvars;          /**< Variable array for all chains */

   SCIP_CONS**           nodeconss;          /**< Constraint array for each node in the original KEP graph
                                              *   representing a patient-donor pair */
};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE KEPdataCreate(
   SCIP*                 kepscip,            /**< SCIP instance of 3rd stage */
   SCIP_PROBDATA**       kepdata,            /**< problem data of 3rd stage */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   Cycles*               cycles,             /**< cycle structure of graph */
   Chains*               chains,             /**< chain structure of graph */
   SCIP_VAR**            cyclevars,          /**< array of cyclevars: per solution, we keep the cycle
                                              *   if none of its vertices is attacked */
   SCIP_VAR**            chainvars,          /**< array of chainvars: per solution, we keep the chain
                                              *   if none of its vertices is attacked */
   SCIP_CONS**           nodeconss           /**< array of linear inequalities, each corresponding to
                                              *   using at most one cycle/chain per vertex */
   )
{
   int ncycles;
   int nchains;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( nnodes > 0 );
   assert( nnodes == graph->nnodes );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(kepscip, kepdata) );

   (*kepdata)->nnodes = nnodes;
   (*kepdata)->graph = graph;
   (*kepdata)->cycles = cycles;
   (*kepdata)->chains = chains;
   ncycles = cycles->ncycles;
   nchains = chains->nchains;

   if ( cyclevars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->cyclevars, cyclevars, ncycles) );
   }
   else
      (*kepdata)->cyclevars = NULL;

   if ( chainvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->chainvars, chainvars, nchains) );
   }
   else
      (*kepdata)->chainvars = NULL;

   /* allocate memory for constraints */
   if ( nodeconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->nodeconss, nodeconss, nnodes) );
   }
   else
      (*kepdata)->nodeconss = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given KEP problem data */
static
SCIP_RETCODE KEPdataFree(
   SCIP*                 kepscip,            /**< SCIP data structure */
   SCIP_PROBDATA**       kepdata             /**< pointer to KEP data */
   )
{
   int nnodes;
   int ncycles;
   int nchains;
   int i;

   assert( kepscip != NULL );
   assert( kepdata != NULL );

   nnodes = (*kepdata)->nnodes;
   ncycles = (*kepdata)->cycles->ncycles;
   nchains = (*kepdata)->chains->nchains;

   assert( nnodes > 0 );
   assert( ncycles >= 0);
   assert( nchains >= 0);

   /* release all variables */
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(kepscip, &(*kepdata)->cyclevars[i]) );
   }

   for (i = 0; i < nchains; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(kepscip, &(*kepdata)->chainvars[i]) );
   }

   /** release constraints and free constraint memory arrays */
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(kepscip, &(*kepdata)->nodeconss[i]) );
   }

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(kepscip, &(*kepdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArray(kepscip, &(*kepdata)->chainvars, nchains);
   SCIPfreeBlockMemoryArrayNull(kepscip, &(*kepdata)->nodeconss, nnodes);

   /* free kepdata block memory */
   SCIPfreeBlockMemory(kepscip, kepdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialKEPVars(
   SCIP*                 kepscip,            /**< SCIP data structure */
   SCIP_PROBDATA*        kepdata             /**< pointer to KEP data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int ncycles;
   int nchains;
   int nnodes;
   int i;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( kepdata->graph != NULL );
   assert( kepdata->graph->nodelist != NULL );

   ncycles = kepdata->cycles->ncycles;
   nchains = kepdata->chains->nchains;
   nnodes = SCIPKEPdataGetNumNodes(kepdata);

   /* create x_c-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &kepdata->cyclevars, ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &kepdata->chainvars, nchains) );

   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(kepscip, &kepdata->cyclevars[i], name, 0.0, 1.0,
            kepdata->cycles->cycleweights[i] * nnodes + 1,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(kepscip, kepdata->cyclevars[i]) );
   }

   for (i = 0; i < nchains; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainvar_%d", i);

      SCIP_CALL( SCIPcreateVar(kepscip, &kepdata->chainvars[i], name, 0.0, 1.0,
            kepdata->chains->chainweights[i] * nnodes + 1,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(kepscip, kepdata->chainvars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the initial constraints of the problem */
static
SCIP_RETCODE SCIPcreateInitialKEPConstraints(
   SCIP*                 kepscip,            /**< SCIP data structure */
   SCIP_PROBDATA*        kepdata             /**< pointer to KEP data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   int nchains;
   Graph* graph;
   Cycles* cycles;
   Chains* chains;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   int maxnvarsincons;
   int i;
   int c;
   int cnt;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( kepdata->graph != NULL );

   nnodes = kepdata->graph->nnodes;
   graph = kepdata->graph;
   cycles = kepdata->cycles;
   ncycles = cycles->ncycles;
   assert( kepdata->cycles->ncycles == cycles->ncycles );

   chains = kepdata->chains;
   nchains = chains->nchains;
   assert( kepdata->chains->nchains == chains->nchains );
   maxnvarsincons = ncycles + nchains;

   /* Allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(kepscip, &vars, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(kepscip, &vals, maxnvarsincons) );

   for (i = 0; i < maxnvarsincons; ++i)
      vals[i] = 1.0;

   /* create constraints linking x_c and its corresponding u_v's */
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &(kepdata->nodeconss), nnodes) );
   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;
      for (c = graph->node2cyclesbegin[i]; c < graph->node2cyclesbegin[i + 1]; ++c)
         vars[cnt++] = kepdata->cyclevars[graph->node2cycles[c]];

      for (c = graph->node2chainsbegin[i]; c < graph->node2chainsbegin[i + 1]; ++c)
         vars[cnt++] = kepdata->chainvars[graph->node2chains[c]];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecons_%d", i);

      /* sum of incident cyclevars and chainvars does not exceed 1*/
      SCIP_CALL( SCIPcreateConsLinear(kepscip, &kepdata->nodeconss[i], name, cnt, vars, vals, 0.0, 1.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(kepscip, kepdata->nodeconss[i]) );
      /* do not release constraints here, will be done later */
   }

   SCIPfreeBufferArray(kepscip, &vars);
   SCIPfreeBufferArray(kepscip, &vals);

   return SCIP_OKAY;
}

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of Benders' decomposition problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigKEP)
{
   SCIPdebugMsg(scip, "free original KEP data\n");

   SCIP_CALL( KEPdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPKEPdataCreate(
   SCIP*                 kepscip,            /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   )
{
   SCIP_PROBDATA* kepdata;

   assert( kepscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(kepscip, probname) );
   SCIP_CALL( SCIPsetProbDelorig(kepscip, probdelorigKEP) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(kepscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create problem data */
   SCIP_CALL( KEPdataCreate(kepscip, &kepdata, graph, graph->nnodes, cycles, chains,
         NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialKEPVars(kepscip, kepdata) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateInitialKEPConstraints(kepscip, kepdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(kepscip, kepdata) );

   return SCIP_OKAY;
}

/** returns the number of nodes in the KEP instance */
int SCIPKEPdataGetNumNodes(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL );

   return kepdata->graph->nnodes;
}

/** returns array of cycle variables */
SCIP_VAR** SCIPKEPdataGetCyclevars(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->cyclevars;
}

/** returns number of cycles */
int SCIPKEPdataGetNCycles(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->cycles->ncycles;
}

/** returns array of chain variables */
SCIP_VAR** SCIPKEPdataGetChainvars(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL );

   return kepdata->chainvars;
}

/** returns number of chains */
int SCIPKEPdataGetNChains(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->chains->nchains;
}

/** returns array of node constraints */
SCIP_CONS** SCIPKEPdataGetNodeconss(
   SCIP_PROBDATA*       kepdata              /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->nodeconss;
}

/**@} */
