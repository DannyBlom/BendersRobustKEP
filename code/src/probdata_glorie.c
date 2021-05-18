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

/**@file   probdata_glorie.c
 * @brief  Problem data for Glorie et al. approach
 * @author Christopher Hojny
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_glorie.h"

#include "scip/scip.h"
#include "scip/cons_setppc.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the kidney exchange problem, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   Graph*                graph;              /**< pointer to underlying graph */
   int                   nnodes;             /**< number of nodes in graph */

   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   Chains*               chains;             /**< pointer to chain structure of underlying graph */

   SCIP_VAR**            cyclevars;          /**< array of cyclevars: per solution, we keep the cycle
                                              *   if none of its vertices is attacked */
   SCIP_VAR**            chainvars;          /**< array of chainvars: per solution, we keep the chain
                                              *   if none of its vertices is attacked */

   SCIP_CONS**           nodeconss;          /**< node constraints */
};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   Cycles*               cycles,             /**< cycle structure of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   SCIP_VAR**            cyclevars,          /**< array of cyclevars */
   SCIP_VAR**            chainvars,          /**< array of chainvars */
   SCIP_CONS**           nodeconss           /**< node constraints */
   )
{
   int ncycles;
   int nchains;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( nnodes == graph->nnodes );
   assert( nnodes > 0 );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   (*probdata)->graph = graph;
   (*probdata)->nnodes = nnodes;
   (*probdata)->cycles = cycles;
   ncycles = cycles->ncycles;
   (*probdata)->chains = chains;
   nchains = chains->nchains;

   /* possible copy variable arrays */
   if ( cyclevars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cyclevars, cyclevars, ncycles) );
   }
   else
      (*probdata)->cyclevars = NULL;

   if ( chainvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->chainvars, chainvars, nchains) );
   }
   else
      (*probdata)->chainvars = NULL;

   /* possible copy constraint arrays */
   if ( nodeconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->nodeconss, nodeconss, nnodes) );
   }
   else
      (*probdata)->nodeconss = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;
   int nnodes;
   int ncycles;
   int nchains;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->graph != NULL );
   assert( (*probdata)->cycles != NULL );

   nnodes = (*probdata)->nnodes;
   ncycles = (*probdata)->cycles->ncycles;
   nchains = (*probdata)->chains->nchains;
   assert( nnodes > 0 );
   assert( ncycles >= 0 );
   assert( nchains >= 0 );

   /* release all variables */
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->cyclevars[i]) );
   }

   for (i = 0; i < nchains; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->chainvars[i]) );
   }

   /** release constraints and free constraint memory arrays */
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->nodeconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->nodeconss, nnodes);

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->chainvars, nchains);

   /* free bendersdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   Cycles* cycles;
   Chains* chains;
   int ncycles;
   int nchains;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );
   assert( probdata->graph->nodelist != NULL );
   assert( probdata->cycles != NULL );
   assert( probdata->chains != NULL );

   ncycles = probdata->cycles->ncycles;
   nchains = probdata->chains->nchains;
   cycles = probdata->cycles;
   chains = probdata->chains;

   /* create x_c-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->cyclevars, ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->chainvars, nchains) );

   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->cyclevars[i], name, 0.0, 1.0, cycles->cycleweights[i],
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->cyclevars[i]) );
   }

   for (i = 0; i < nchains; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainvar_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->chainvars[i], name, 0.0, 1.0, chains->chainweights[i],
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->chainvars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the constraints of the problem */
static
SCIP_RETCODE SCIPcreateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   Graph* graph;
   Cycles* cycles;
   Chains* chains;
   SCIP_VAR** vars;
   int nnodes;

   int maxnvarsincons;
   int i;
   int c;
   int cnt;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );

   graph = probdata->graph;
   nnodes = probdata->graph->nnodes;
   cycles = probdata->cycles;
   chains = probdata->chains;

   /* Allocate buffer memory arrays to add conss to the model */
   maxnvarsincons = cycles->ncycles + chains->nchains;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvarsincons) );

   /* create node constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->nodeconss), nnodes) );
   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;
      for (c = graph->node2cyclesbegin[i]; c < graph->node2cyclesbegin[i + 1]; ++c)
         vars[cnt++] = probdata->cyclevars[graph->node2cycles[c]];
      for (c = graph->node2chainsbegin[i]; c < graph->node2chainsbegin[i + 1]; ++c)
         vars[cnt++] = probdata->chainvars[graph->node2chains[c]];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecons_%d", i);

      SCIP_CALL( SCIPcreateConsSetpack(scip, &probdata->nodeconss[i], name, cnt, vars,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, probdata->nodeconss[i]) );
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of Glorie model (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigGlorie)
{
   SCIPdebugMsg(scip, "free original Glorie data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed Glorie model (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransGlorie)
{
   int ncycles;
   int nchains;

   /* create transform bendersdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->chains,
         sourcedata->cyclevars, sourcedata->chainvars, sourcedata->nodeconss) );

   ncycles = sourcedata->cycles->ncycles;
   nchains = sourcedata->chains->nchains;

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->graph->nnodes, (*targetdata)->nodeconss, (*targetdata)->nodeconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->cyclevars, (*targetdata)->cyclevars) );
   SCIP_CALL( SCIPtransformVars(scip, nchains, (*targetdata)->chainvars, (*targetdata)->chainvars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransGlorie)
{
   SCIPdebugMsg(scip, "free transformed Glorie data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}
/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPglorieEtAldataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   )
{
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigGlorie) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransGlorie) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransGlorie) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, graph, graph->nnodes, cycles, chains,
         NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateVariables(scip, probdata) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateConstraints(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   return SCIP_OKAY;
}

/** returns cycle variables */
SCIP_VAR** getCyclevarsGlorie(
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   assert( probdata != NULL );

   return probdata->cyclevars;
}

/** returns chain variables */
SCIP_VAR** getChainvarsGlorie(
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   assert( probdata != NULL );

   return probdata->chainvars;
}


/**@} */
