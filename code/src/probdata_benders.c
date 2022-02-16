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

/**@file   probdata_benders.c
 * @brief  Problem data for Benders type approach
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 *
 * @page BENDERS_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_PROBDATA is shown below.
 *
 * A list of all interface methods can be found in probdata_benders.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "graph.h"
#include "probdata_benders.h"
#include "probdata_master_kidney.h"
#include "typedefs.h"

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
   int                   nnodes;             /**< number of nodes in graph */
   int                   adversarybound;     /**< bound on adversary attack */

   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   Chains*               chains;             /**< pointer to chain structure of underlying graph */

   SCIP_VAR*             objvar;             /**< artificial objective variable modeling integer objective values */
   SCIP_VAR**            cyclevars;          /**< array of cyclevars: per solution, we keep the cycle if none of its vertices is attacked */
   SCIP_VAR**            chainvars;          /**< array of chainvars: per solution, we keep the chain if none of its vertices is attacked */
   SCIP_VAR**            uvars;              /**< array of u-vars: adversary variables for attacking nodes */

   int                   nsolconss;          /**< number of generated KEP solutions forming a constraint */
   int                   maxnsolconss;       /**< maximum number of generated KEP solutions that fits into solconss array */

   SCIP_CONS*            attackboundcons;    /**< constraint bounding adversary attacks */
   SCIP_CONS**           solconss;           /**< array of linear inequalities, each corresponding to a kidney exchange solution */
   SCIP_CONS**           cyclelbconss;       /**< cycle constraints (cyclevars[c] = 1 if all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           cycleubconss;       /**< cycle constraints (cyclevars[c] = 0 if any uvars[v] = 0 for v in cycle c) */
   SCIP_CONS**           chainlbconss;       /**< chain constraints (chainvars[c] = 1 if all uvars[v] = 1 for v in chain c)*/
   SCIP_CONS**           chainubconss;       /**< chain constraints (chainvars[c] = 0 if any uvars[v] = 0 for v in chain c) */

};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE bendersdataCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_PROBDATA**       bendersdata,        /**< pointer to Benders data */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   Cycles*               cycles,             /**< cycle structure of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */

   SCIP_VAR**            cyclevars,          /**< array of cyclevars: per solution, we keep the cycle if none of its vertices is attacked */
   SCIP_VAR**            chainvars,          /**< array of chainvars: per solution, we keep the chain if none of its vertices is attacked */
   SCIP_VAR**            uvars,              /**< array of u-vars: adversary variables for attacking nodes */
   SCIP_VAR*             objvar,             /**< artificial objective variable modeling integer objective values */

   int                   nsolconss,          /**< number of generated KEP solutions forming a constraint */
   int                   maxnsolconss,       /**< maximum number of generated KEP solutions that fits into solconss array */
   int                   adversarybound,     /**< bound on adversary attack */

   SCIP_CONS**           solconss,           /**< array of linear inequalities, each corresponding to a kidney exchange solution */
   SCIP_CONS*            attackboundcons,    /**< constraint bounding adversary attacks */
   SCIP_CONS**           cyclelbconss,       /**< cycle constraints (cyclevars[c] = 1 if all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           cycleubconss,       /**< cycle constraints (cyclevars[c] = 0 if any uvars[v] = 0 for v in cycle c) */
   SCIP_CONS**           chainlbconss,       /**< chain constraints (chainvars[c] = 1 if all uvars[v] = 1 for v in chain c)*/
   SCIP_CONS**           chainubconss        /**< chain constraints (chainvars[c] = 0 if any uvars[v] = 0 for v in chain c) */
   )
{
   int ncycles;
   int nchains;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( nsolconss >= 0 );
   assert( adversarybound >= 0 );
   assert( nnodes == graph->nnodes );
   assert( nnodes > 0 );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(bendersscip, bendersdata) );

   (*bendersdata)->graph = graph;
   (*bendersdata)->nnodes = nnodes;
   (*bendersdata)->adversarybound = adversarybound;
   (*bendersdata)->cycles = cycles;
   ncycles = cycles->ncycles;
   (*bendersdata)->chains = chains;
   nchains = chains->nchains;

   /* possible copy variable arrays */
   if( uvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->uvars, uvars, nnodes) );
   }
   else
      (*bendersdata)->uvars = NULL;

   (*bendersdata)->objvar = objvar;

   if( cyclevars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cyclevars, cyclevars, ncycles) );
   }
   else
      (*bendersdata)->cyclevars = NULL;

   if( chainvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->chainvars, chainvars, nchains) );
   }
   else
      (*bendersdata)->chainvars = NULL;

   /* possible copy constraint arrays */
   if( nsolconss > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->solconss, solconss, maxnsolconss) );
      (*bendersdata)->nsolconss = nsolconss;
      (*bendersdata)->maxnsolconss = maxnsolconss;
   }
   else
   {
      (*bendersdata)->solconss = NULL;
      (*bendersdata)->nsolconss = 0;
      (*bendersdata)->maxnsolconss = 0;
   }

   if( attackboundcons != NULL )
      (*bendersdata)->attackboundcons = attackboundcons;
   else
      (*bendersdata)->attackboundcons = NULL;

   if( cyclelbconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cyclelbconss, cyclelbconss, ncycles) );
   }
   else
      (*bendersdata)->cyclelbconss = NULL;

   if( cycleubconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cycleubconss, cycleubconss, cycles->nnodesincycles) );
   }
   else
      (*bendersdata)->cycleubconss = NULL;

   if( chainlbconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->chainlbconss, chainlbconss, nchains) );
   }
   else
      (*bendersdata)->chainlbconss = NULL;

   if( chainubconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->chainubconss, chainubconss, chains->nnodesinchains) );
   }
   else
      (*bendersdata)->chainubconss = NULL;

   return SCIP_OKAY;
}


/** frees the memory of the given problem data */
static
SCIP_RETCODE bendersdataFree(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_PROBDATA**       bendersdata         /**< pointer to Benders data */
   )
{
   int i;
   int nnodes;
   int ncycles;
   int nnodesincycles;
   int nchains;
   int nnodesinchains;
   int nsolconss;
   int policy;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( (*bendersdata)->graph != NULL );
   assert( (*bendersdata)->cycles != NULL );

   nnodes = (*bendersdata)->nnodes;
   ncycles = (*bendersdata)->cycles->ncycles;
   nnodesincycles = (*bendersdata)->cycles->nnodesincycles;
   nchains = (*bendersdata)->chains->nchains;
   nnodesinchains = (*bendersdata)->chains->nnodesinchains;
   nsolconss = (*bendersdata)->nsolconss;

   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );

   assert( nnodes > 0 );
   assert( ncycles >= 0 );
   assert( nchains >= 0 );

   /* release all variables */
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->uvars[i]) );
   }
   SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->objvar) );

   for( i = 0; i < ncycles; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->cyclevars[i]) );
   }

   for( i = 0; i < nchains; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->chainvars[i]) );
   }

   /** release constraints and free constraint memory arrays */
   for( i = 0; i < nsolconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->solconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->solconss, (*bendersdata)->maxnsolconss);

   for( i = 0; i < ncycles; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->cyclelbconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->cyclelbconss, ncycles);

   for( i = 0; i < nchains; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->chainlbconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->chainlbconss, nchains);

   for( i = 0; i < nnodesincycles; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->cycleubconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->cycleubconss, nnodesincycles);

   for( i = 0; i < nnodesinchains; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->chainubconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->chainubconss, nnodesinchains);

   SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->attackboundcons) );

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->uvars, nnodes);
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->chainvars, nchains);

   /* free bendersdata */
   SCIPfreeBlockMemory(bendersscip, bendersdata);

   return SCIP_OKAY;
}

/** returns whether a given cycle intersects a target cycle / chain in initial solution */
SCIP_Bool SCIPcycleIntersectsTarget(
   Cycles*               cycles,
   Chains*               chains,
   int                   cycle_idx,
   int                   target_idx
   )
{
   int i;
   int j;

   assert( 0 <= cycle_idx && cycle_idx < cycles->ncycles );
   assert( 0 <= target_idx && target_idx < cycles->ncycles + chains->nchains );

   for( i = cycles->nodelistsbegin[cycle_idx]; i < cycles->nodelistsbegin[cycle_idx + 1]; ++i )
   {
      if( target_idx < cycles->ncycles )
      {
         for( j = cycles->nodelistsbegin[target_idx]; j < cycles->nodelistsbegin[target_idx + 1]; ++j )
         {
            if( cycles->nodelists[i] == cycles->nodelists[j] )
               return TRUE;
         }
      }
      else
      {
         for( j = chains->nodelistsbegin[target_idx - cycles->ncycles]; j < chains->nodelistsbegin[target_idx - cycles->ncycles + 1]; ++j )
         {
            if( cycles->nodelists[i] == chains->nodelists[j] )
               return TRUE;
         }
      }
   }
   return FALSE;
}

/**< returns whether chain2 is a superchain of chain1 */
SCIP_Bool SCIPchainIsSuperchain(
   Chains*              chains,              /**< chains data structure */
   int                  chain1,              /**< index of chain */
   int                  chain2               /**< index of possible superchain */
)
{
   int current_c;
   int subchainidx;
   int chainlen1;
   int chainlen2;

   chainlen1 = chains->nodelistsbegin[chain1 + 1] - chains->nodelistsbegin[chain1];
   chainlen2 = chains->nodelistsbegin[chain2 + 1] - chains->nodelistsbegin[chain2];

   if( chainlen1 >= chainlen2 )
      return FALSE;

   current_c = chain2;
   subchainidx = chains->subchains[current_c];
   while( subchainidx != -1 )
   {
      if( subchainidx == chain1 )
         return TRUE;
      current_c = subchainidx;
      subchainidx = chains->subchains[current_c];
   }
   return FALSE;
}

/** returns whether a given chain intersects a target cycle / chain in initial solution */
SCIP_Bool SCIPchainIntersectsTarget(
   Cycles*               cycles,
   Chains*               chains,
   int                   chain_idx,
   int                   target_idx
   )
{
   int i;
   int j;

   assert( 0 <= chain_idx && chain_idx < chains->nchains );
   assert( 0 <= target_idx && target_idx < cycles->ncycles + chains->nchains );

   for( i = chains->nodelistsbegin[chain_idx]; i < chains->nodelistsbegin[chain_idx + 1]; ++i )
   {
      if( target_idx < cycles->ncycles )
      {
         for( j = cycles->nodelistsbegin[target_idx]; j < cycles->nodelistsbegin[target_idx + 1]; ++j )
         {
            if( chains->nodelists[i] == cycles->nodelists[j] )
               return TRUE;
         }
      }
      else
      // both chain and target are chain. Only return true if target is NOT a superchain of chain
      {
         if( SCIPchainIsSuperchain(chains, chain_idx, target_idx - cycles->ncycles) )
            return FALSE;
         else
         {
            for( j = chains->nodelistsbegin[target_idx - cycles->ncycles]; j < chains->nodelistsbegin[target_idx - cycles->ncycles + 1]; ++j )
            {
               if( chains->nodelists[i] == chains->nodelists[j] )
                  return TRUE;
            }
         }
      }
   }
   return FALSE;
}


/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialBendersVars(
   SCIP*                 bendersscip,        /**< SCIP instance of Benders model */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   int nchains;
   int i;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );
   assert( bendersdata->graph->nodelist != NULL );

   nnodes = bendersdata->graph->nnodes;
   ncycles = bendersdata->cycles->ncycles;
   nchains = bendersdata->chains->nchains;

   /* create x_c-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->cyclevars, ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->chainvars, nchains) );

   for( i = 0; i < ncycles; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->cyclevars[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->cyclevars[i]) );
   }

   for( i = 0; i < nchains; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainvar_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->chainvars[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->chainvars[i]) );
   }

   /* create u-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->uvars, nnodes) );
   for( i = 0; i < nnodes; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "u_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->uvars[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->uvars[i]) );
   }

   /* create artificial objective variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "objvar");

   SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->objvar, name, 0.0, nnodes, 1.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->objvar) );

   return SCIP_OKAY;
}

/** creates the initial constraints of the benders problem given the full recourse policy */
static
SCIP_RETCODE SCIPcreateInitialBendersConstraintsFR(
   SCIP*                 bendersscip,        /**< SCIP instance of Benders model */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   int nnodesincycles;
   int nchains;
   int nnodesinchains;
   Cycles* cycles;
   Chains* chains;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool initial;
   int policy;

   int maxnvarsincons;
   int i;
   int c;
   int cnt;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );

   nnodes = bendersdata->graph->nnodes;
   maxnvarsincons = nnodes + 1;

   cycles = bendersdata->cycles;
   ncycles = cycles->ncycles;
   nnodesincycles = cycles->nnodesincycles;
   assert( bendersdata->cycles->ncycles == cycles->ncycles );

   chains = bendersdata->chains;
   nchains = chains->nchains;
   nnodesinchains = chains->nnodesinchains;
   assert( bendersdata->chains->nchains == chains->nchains );

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/cyclechainconssinitial", &initial) );
   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );
   assert( policy == POLICY_FULLRECOURSE );

   /* Allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vals, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vars, maxnvarsincons) );


   /* create bound on adversary attack */
   for( i = 0; i < maxnvarsincons; ++i )
      vals[i] = 1.0;

   SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->attackboundcons, "boundattack", nnodes,
         bendersdata->uvars, vals, nnodes - bendersdata->adversarybound, nnodes,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->attackboundcons) );
   /* do not release constraint here, will be done later */

   /* create constraints linking x_c and its corresponding u_v's */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cyclelbconss), ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cycleubconss), nnodesincycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->chainlbconss), nchains) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->chainubconss), nnodesinchains) );

   for( c = 0; c < ncycles; ++c )
   {
      cnt = 0;
      for( i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i )
         vars[cnt++] = bendersdata->uvars[cycles->nodelists[i]];

      vars[cnt] = bendersdata->cyclevars[c];
      vals[cnt] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclelbcons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cyclelbconss[c], name, cnt + 1, vars, vals,
            -1.0, cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cyclelbconss[c]) );

      /* Reset value of vals[cnt] for next cons */
      vals[cnt] = 1.0;

      /* Secondly, add upper bound constraints x_c <= u_v (observe it is initially inactive) */
      vars[0] = bendersdata->cyclevars[c];
      vals[1] = -1.0;
      for( i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c+1]; ++i )
      {
         vars[1] = bendersdata->uvars[cycles->nodelists[i]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cycleubcons_%d", i);

         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cycleubconss[i], name, 2, vars, vals,
            -SCIPinfinity(bendersscip), 0,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cycleubconss[i]) );

      }
      /* Reset negative value */
      vals[1] = 1.0;
   }

   for( c = 0; c < nchains; ++c )
   {
      cnt = 0;
      for( i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c + 1]; ++i )
         vars[cnt++] = bendersdata->uvars[chains->nodelists[i]];

      vars[cnt] = bendersdata->chainvars[c];
      vals[cnt] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainlbcons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->chainlbconss[c], name, cnt + 1, vars, vals,
            -1.0, chains->nodelistsbegin[c + 1] - chains->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->chainlbconss[c]) );

      /* Reset value of vals[cnt] for next cons */
      vals[cnt] = 1.0;

      /* Secondly, add upper bound constraints x_c <= u_v (observe it is initially inactive) */
      vars[0] = bendersdata->chainvars[c];
      vals[1] = -1.0;
      for( i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c+1]; ++i )
      {
         vars[1] = bendersdata->uvars[chains->nodelists[i]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainubcons_%d", i);

         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->chainubconss[i], name, 2, vars, vals,
            -SCIPinfinity(bendersscip), 0,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->chainubconss[i]) );

      }
      /* Reset negative value */
      vals[1] = 1.0;
   }

   SCIPfreeBufferArray(bendersscip, &vars);
   SCIPfreeBufferArray(bendersscip, &vals);

   return SCIP_OKAY;
}

/** creates the initial constraints of the benders problem given the Keep Unaffected Cycle Chain policy */
static
SCIP_RETCODE SCIPcreateInitialBendersConstraintsKUCC(
   SCIP*                 bendersscip,        /**< SCIP instance of Benders model */
   SCIP*                 masterscip,         /**< SCIP instance of master model */
   SCIP_SOL*             sol,                /**< SCIP solution to master model */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   int nchains;
   int nnodesincycles;
   int nnodesinchains;
   Cycles* cycles;
   Chains* chains;
   SCIP_VAR** vars;
   SCIP_VAR** initxvars;
   SCIP_Real* vals;
   int* initcc;
   SCIP_Bool initial;
   int policy;

   int maxnvarsincons;
   int i;
   int c;
   int ninitcc;
   int cnt;
   int cc;
   int start;
   int cc_idx;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );

   nnodes = bendersdata->graph->nnodes;
   maxnvarsincons = nnodes + 1;

   cycles = bendersdata->cycles;
   ncycles = cycles->ncycles;
   nnodesincycles = cycles->nnodesincycles;
   assert( bendersdata->cycles->ncycles == cycles->ncycles );

   chains = bendersdata->chains;
   nchains = chains->nchains;
   nnodesinchains = chains->nnodesinchains;
   assert( bendersdata->chains->nchains == chains->nchains );

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/cyclechainconssinitial", &initial) );
   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );
   assert( policy == POLICY_KEEPUNAFFECTEDCC );

   /* Allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vals, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vars, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &initcc, nnodes / 2) );

   /* create bound on adversary attack */
   for( i = 0; i < maxnvarsincons; ++i )
      vals[i] = 1.0;

   SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->attackboundcons, "boundattack", nnodes,
         bendersdata->uvars, vals, nnodes - bendersdata->adversarybound, nnodes,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->attackboundcons) );

   /* create constraints linking x_c and its corresponding u_v's */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cyclelbconss), ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cycleubconss), nnodesincycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->chainlbconss), nchains) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->chainubconss), nnodesinchains) );

   /* First identify indices of cycles and chains in initial solution and add their corresponding cycle / chain constraints */
   initxvars = masterProblemGetXvarinit(SCIPgetProbData(masterscip));
   ninitcc = 0;

   printf("Initial solution:\n");
   for( c = 0; c < ncycles + nchains; ++c )
   {
      if( SCIPgetSolVal(masterscip, sol, initxvars[c]) > 0.5 )
      {
         initcc[ninitcc++] = c;
         if( c < ncycles)
            printf("Cycle %d\n", c);
         else
            printf("Chain %d\n", c - cycles->ncycles);
      }
   }

   for( c = 0; c < ncycles; ++c )
   {
      /* Start by adding lower bound constraints x_c >= 1 - sum_{v in V(c)}  u_v */
      cnt = 0;
      for( i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i )
         vars[cnt++] = bendersdata->uvars[cycles->nodelists[i]];

      vals[cnt] = -1.0;
      vars[cnt] = bendersdata->cyclevars[c];
      start = cnt;

      if( SCIPgetSolVal(masterscip, sol, initxvars[c]) <= 0.5 )
      {
         // Also take into account initial cycles chains that intersect this cycle, since these have priority in the Keep Unaffected Cycle Chain policy
         for( cc = 0; cc < ninitcc; ++cc )
         {
            cc_idx = initcc[cc];
            if( SCIPcycleIntersectsTarget(cycles, chains, c, cc_idx) )
            {
               cnt++;
               vals[cnt] = -1.0;
               if( cc_idx < cycles->ncycles )
                  vars[cnt] = bendersdata->cyclevars[cc_idx];
               else
                  vars[cnt] = bendersdata->chainvars[cc_idx - cycles->ncycles];
            }
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclelbcons_%d", c);

      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cyclelbconss[c], name, cnt + 1, vars, vals,
         -SCIPinfinity(bendersscip), cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
         initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cyclelbconss[c]) );

      /* Reset negative values of vals[i] to one for next cons */
      for( i = start; i <= cnt; ++i )
         vals[i] = 1.0;


      /* Secondly, add upper bound constraints x_c <= u_v */
      vars[0] = bendersdata->cyclevars[c];
      vals[1] = -1.0;
      for( i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c+1]; ++i )
      {
         vars[1] = bendersdata->uvars[cycles->nodelists[i]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cycleubcons_%d", i);

         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cycleubconss[i], name, 2, vars, vals,
            -SCIPinfinity(bendersscip), 0,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cycleubconss[i]) );

      }
      /* Reset negative value */
      vals[1] = 1.0;

   }

   for( c = 0; c < nchains; ++c )
   {
      /* Start by adding lower bound constraints x_c >= 1 - sum_{v in V(c)} (1 - u_v) */
      cnt = 0;
      for( i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c + 1]; ++i )
         vars[cnt++] = bendersdata->uvars[chains->nodelists[i]];

      vals[cnt] = -1.0;
      vars[cnt] = bendersdata->chainvars[c];
      start = cnt;

      if( SCIPgetSolVal(masterscip, sol, initxvars[ncycles + c]) <= 0.5 )
      {
         // Also take into account initial cycles chains that intersect this cycle, since these have priority in the Keep Unaffected Cycle Chain policy
         for( cc = 0; cc < ninitcc; ++cc )
         {
            cc_idx = initcc[cc];
            if( SCIPchainIntersectsTarget(cycles, chains, c, cc_idx) )
            {
               cnt++;
               vals[cnt] = -1.0;
               if( cc_idx < cycles->ncycles )
                  vars[cnt] = bendersdata->cyclevars[cc_idx];
               else
                  vars[cnt] = bendersdata->chainvars[cc_idx - cycles->ncycles];
            }
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainlbcons_%d", c);

      /* x_c >= 1 - sum u_v */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->chainlbconss[c], name, cnt + 1, vars, vals,
         -SCIPinfinity(bendersscip), chains->nodelistsbegin[c + 1] - chains->nodelistsbegin[c] - 1,
         initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->chainlbconss[c]) );

      /* Reset negative values of vals[i] to one for next cons */
      for( i = start; i <= cnt; ++i )
         vals[i] = 1.0;

      /* Secondly, add upper bound constraints x_c <= u_v */
      vars[0] = bendersdata->chainvars[c];
      vals[1] = -1.0;
      for( i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c+1]; ++i )
      {
         vars[1] = bendersdata->uvars[chains->nodelists[i]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainubcons_%d", i);

         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->chainubconss[i], name, 2, vars, vals,
            -SCIPinfinity(bendersscip), 0,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->chainubconss[i]) );

      }
      /* Reset negative value */
      vals[1] = 1.0;
   }

   SCIPfreeBufferArray(bendersscip, &vars);
   SCIPfreeBufferArray(bendersscip, &vals);
   SCIPfreeBufferArray(bendersscip, &initcc);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of Benders' decomposition problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigBenders)
{
   SCIPdebugMsg(scip, "free original Benders data\n");

   SCIP_CALL( bendersdataFree(scip, probdata) );

   return SCIP_OKAY;
}


/** creates user data of transformed Benders problem by transforming the original Benders problem data
 *  (called after problem was transformed)
 */
static
SCIP_DECL_PROBTRANS(probtransBenders)
{
   int nnodes;
   int ncycles;
   int nnodesincycles;
   int nchains;
   int nnodesinchains;
   int policy;

   /* create transform bendersdata */
   SCIP_CALL( bendersdataCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->chains,
         sourcedata->cyclevars, sourcedata->chainvars, sourcedata->uvars, sourcedata->objvar, sourcedata->nsolconss, sourcedata->maxnsolconss,
         sourcedata->adversarybound, sourcedata->solconss, sourcedata->attackboundcons, sourcedata->cyclelbconss, sourcedata->cycleubconss,
         sourcedata->chainlbconss, sourcedata->chainubconss) );

   nnodes = sourcedata->graph->nnodes;
   ncycles = sourcedata->cycles->ncycles;
   nnodesincycles = sourcedata->cycles->nnodesincycles;
   nchains = sourcedata->chains->nchains;
   nnodesinchains = sourcedata->chains->nnodesinchains;

   SCIP_CALL( SCIPgetIntParam(scip, "kidney/recoursepolicy", &policy) );

   assert( nnodes > 0 );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nsolconss, (*targetdata)->solconss, (*targetdata)->solconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->attackboundcons, &(*targetdata)->attackboundcons) );
   SCIP_CALL( SCIPtransformConss(scip, ncycles, (*targetdata)->cyclelbconss, (*targetdata)->cyclelbconss) );
   SCIP_CALL( SCIPtransformConss(scip, nchains, (*targetdata)->chainlbconss, (*targetdata)->chainlbconss) );
   SCIP_CALL( SCIPtransformConss(scip, nnodesincycles, (*targetdata)->cycleubconss, (*targetdata)->cycleubconss) );
   SCIP_CALL( SCIPtransformConss(scip, nnodesinchains, (*targetdata)->chainubconss, (*targetdata)->chainubconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->cyclevars, (*targetdata)->cyclevars) );
   SCIP_CALL( SCIPtransformVars(scip, nchains, (*targetdata)->chainvars, (*targetdata)->chainvars) );
   SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->uvars, (*targetdata)->uvars) );
   SCIP_CALL( SCIPtransformVar(scip, (*targetdata)->objvar, &(*targetdata)->objvar) );

   return SCIP_OKAY;
}


/** frees user data of transformed Benders problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransBenders)
{
   SCIPdebugMsg(scip, "free transformed Benders' decomposition data\n");

   SCIP_CALL( bendersdataFree(scip, probdata) );

   return SCIP_OKAY;
}
/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPbendersdataCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure of master problem */
   const char*           probname,           /**< problem name */
   SCIP_SOL*             sol,                /**< SCIP solution to masterscip */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* bendersdata;
   int policy;

   assert( bendersscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(bendersscip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(bendersscip, probdelorigBenders) );
   SCIP_CALL( SCIPsetProbTrans(bendersscip, probtransBenders) );
   SCIP_CALL( SCIPsetProbDeltrans(bendersscip, probdeltransBenders) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(bendersscip, SCIP_OBJSENSE_MINIMIZE) );

   /* create problem data */
   SCIP_CALL( bendersdataCreate(bendersscip, &bendersdata, graph, graph->nnodes, cycles, chains,
         NULL, NULL, NULL, NULL, 0, 0, adversarybound, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialBendersVars(bendersscip, bendersdata) );

   /* create initial constraints */
   if( policy == POLICY_FULLRECOURSE )
      SCIP_CALL( SCIPcreateInitialBendersConstraintsFR(bendersscip, bendersdata) );
   else if( policy == POLICY_KEEPUNAFFECTEDCC )
      SCIP_CALL( SCIPcreateInitialBendersConstraintsKUCC(bendersscip, masterscip, sol, bendersdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(bendersscip, bendersdata) );

   return SCIP_OKAY;
}


/** adds given variable to the problem data */
SCIP_RETCODE SCIPbendersdataAddSolCons(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   int                   nvars               /**< number of variables */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_PROBDATA* bendersdata;
   assert( bendersscip != NULL );
   assert( vals != NULL );
   assert( vars != NULL );

   bendersdata = SCIPgetProbData(bendersscip);
   assert( bendersdata != NULL );

   /* check if enough memory is left */
   if( bendersdata->maxnsolconss <= bendersdata->nsolconss )
   {
      int newsize;
      if( bendersdata->maxnsolconss == 0 )
         newsize = 100;
      else
         newsize = bendersdata->maxnsolconss * 2;
      SCIP_CALL( SCIPreallocBlockMemoryArray(bendersscip, &bendersdata->solconss, bendersdata->maxnsolconss, newsize) );
      bendersdata->maxnsolconss = newsize;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "solcons_%d", bendersdata->nsolconss);

   SCIP_CALL( SCIPcreateConsBasicLinear(bendersscip, &bendersdata->solconss[bendersdata->nsolconss], name, nvars,
         vars, vals, -SCIPinfinity(bendersscip), 0) );

   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->solconss[bendersdata->nsolconss]) );
   bendersdata->nsolconss++;

   SCIPdebugMsg(bendersscip, "added solution constraint to bendersdata; nsolconss = %d\n", bendersdata->nsolconss);
   SCIPinfoMessage(bendersscip, NULL, "added solution constraint to bendersdata; nsolconss = %d\n", bendersdata->nsolconss);

   return SCIP_OKAY;
}


/** returns array of cycle lower bound constraints */
SCIP_CONS** SCIPbendersdataGetCyclelbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );
   return bendersdata->cyclelbconss;
}

/** returns array of cycle upper bound constraints */
SCIP_CONS** SCIPbendersdataGetCycleubconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );
   return bendersdata->cycleubconss;
}


/** returns array of cycle variables */
SCIP_VAR** SCIPbendersdataGetCyclevars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );
   return bendersdata->cyclevars;
}


/** returns number of cycles */
int SCIPbendersdataGetNCycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->cycles->ncycles;
}


/** returns array of chain lower bound constraints */
SCIP_CONS** SCIPbendersdataGetChainlbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->chainlbconss;
}

/** returns array of chain upper bound constraints */
SCIP_CONS** SCIPbendersdataGetChainubconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->chainubconss;
}


/** returns array of chain variables */
SCIP_VAR** SCIPbendersdataGetChainvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );
   return bendersdata->chainvars;
}


/** returns number of chains */
int SCIPbendersdataGetNChains(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->chains->nchains;
}


/** returns uvars */
SCIP_VAR** SCIPbendersdataGetUVars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->uvars;
}


/** returns budget of adversary */
int SCIPbendersdataGetAdversaryBound(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert (bendersdata != NULL );

   return bendersdata->adversarybound;
}


/** returns number of nodes in instance */
int SCIPbendersdataGetNumNodes(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert (bendersdata != NULL );

   return bendersdata->graph->nnodes;
}


/** returns objvar */
SCIP_VAR* SCIPbendersdataGetObjvar(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert (bendersdata != NULL );

   return bendersdata->objvar;
}

/**@} */
