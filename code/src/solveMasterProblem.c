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

/**@file   solveMasterProblem.c
 * @brief  Methods to solve the master problem
 * @author Christopher Hojny
 *
 * This file provides routines to solve the master problem.
 */

#include <string.h>
#include <time.h>

#include "auxiliaryStructures.h"
#include "graph.h"
#include "probdata_master_kidney.h"
#include "probdata_benders.h"
#include "probdata_benders_subproblem.h"
#include "probdata_benders_picef.h"
#include "probdata_glorie.h"
#include "problem_kidneyexchange.h"
#include "problem_benders.h"
#include "problem_glorie.h"
#include "solveMasterProblem.h"
#include "kidneyParams.h"
#include "kidneyPlugins.h"
#include "typedefs.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"

/** returns whether the solving loop has to be terminated
 *
 * criteria for termination are:
 * - SCIP instance could not be solved to optimality if applicable
 * - current time exceeds maximum solution time allowed
 */
static
SCIP_Bool solvingloopShallTerminate(
   SCIP*                 scip,               /**< SCIP instance or NULL */
   clock_t               starttime,          /**< start time of loop */
   SCIP_Real             timelimit           /**< maximum solution time */
   )
{
   clock_t curtime;
   SCIP_Real totaltime;

   curtime = clock();
   totaltime = (SCIP_Real) (curtime - starttime) / CLOCKS_PER_SEC;

   /* terminate if the instance could not be solved to optimality */
   if ( scip != NULL && SCIPgetStatus(scip) != SCIP_STATUS_OPTIMAL )
   {
      SCIP_CALL( SCIPprintStatus(scip, NULL) );
      assert( SCIPgetStatus(scip) < SCIP_STATUS_OPTIMAL );

      return TRUE;
   }

   /* terminate if we hit the time limit*/
   if ( totaltime >= timelimit )
      return TRUE;

  return FALSE;
}

/** sets time/memory limit of instance and verbosity parameter */
static
SCIP_RETCODE SCIPsetLimitsAndVerbose(
   SCIP*                 scip,               /**< SCIP instance */
   clock_t               starttime,          /**< start time of loop */
   SCIP_Real             timelimit,          /**< maximum solution time */
   SCIP_Real             memlimit,           /**< memory limit */
   SCIP_Bool             verbose             /**< whether we print log output */
   )
{
   clock_t curtime;
   SCIP_Real totaltime;

   assert( scip != NULL );

   curtime = clock();
   totaltime = (SCIP_Real) (curtime - starttime) / CLOCKS_PER_SEC;

   assert( timelimit >= totaltime );
   assert( timelimit - totaltime <= 1e20 );

   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MAX(0, timelimit - totaltime)) );
   SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", MAX(0, memlimit)) );
   if ( ! verbose )
      SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   return SCIP_OKAY;
}

/** create SCIP instance with kidney parameters and standard parameter settings */
static
SCIP_RETCODE SCIPcreateBasicInstance(
   SCIP*                 parentscip,         /**< parent SCIP instance to copy parameter values from */
   SCIP**                scip,               /**< pointer to SCIP instance */
   SCIP_Bool             verbose             /**< whether we print log output */
   )
{
   assert( parentscip != NULL );
   assert( scip != NULL );

   /* create new SCIP instance for stage 2/3 */
   SCIP_CALL( SCIPcreate(scip) );
   SCIP_CALL( includeKidneyPlugins(*scip) );

   /* read parameters */
   SCIP_CALL( setSCIPParameters(*scip) );
   SCIP_CALL( addKidneyParameters(*scip) );
   SCIP_CALL( SCIPcopyParamSettings(parentscip, *scip) );
   if ( ! verbose )
      SCIP_CALL( SCIPsetIntParam(*scip, "display/verblevel", 0) );

   SCIP_CALL( SCIPsetBoolParam(*scip, "misc/catchctrlc", FALSE) );

   return SCIP_OKAY;
}


/** adapt cycle and chain weights based on solution of master problem */
static
SCIP_RETCODE adaptCycleChainWeights(
   SCIP*                 masterscip,         /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution of masterscip */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< underlying cycles structure */
   Chains*               chains              /**< underlying chains structure */
   )
{
   SCIP_VAR** dummyyvars;
   SCIP_VAR* var;
   SCIP_PROBDATA* probdata;
   int ncycles;
   int nchains;
   int npairs;
   int c;
   int v;

   assert( masterscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   probdata = SCIPgetProbData(masterscip);
   assert( probdata != NULL );

   dummyyvars = masterProblemGetDummyYVars(probdata);
   assert( dummyyvars != NULL );

   ncycles = cycles->ncycles;
   nchains = chains->nchains;
   npairs = graph->npairs;

   /* reset cycle / chain weights */
   for (c = 0; c < ncycles; ++c)
      cycles->cycleweights[c] = 0;

   for (c = 0; c < nchains; ++c)
      chains->chainweights[c] = 0;

   /* For every pair in the graph, check if it is involved in the initial solution
    * If so, increase weights of cycles / chains containing the pair by 1
    */
   for (v = 0; v < npairs; ++v)
   {
      var = dummyyvars[v];

      /* skip node not involved in initial solution */
      if ( SCIPgetSolVal(masterscip, sol, var) < 0.5 )
         continue;

      /* Increment cycle weights */
      for (c = graph->node2cyclesbegin[v]; c < graph->node2cyclesbegin[v + 1]; ++c)
      {
         assert( graph->node2cycles[c] >= 0 );
         assert( graph->node2cycles[c] < ncycles );
         cycles->cycleweights[graph->node2cycles[c]] += 1;
      }

      /* Increment chain weights */
      for (c = graph->node2chainsbegin[v]; c < graph->node2chainsbegin[v + 1]; ++c)
      {
         assert( graph->node2chains[c] >= 0 );
         assert( graph->node2chains[c] < nchains );
         chains->chainweights[graph->node2chains[c]] += 1;
      }
   }

   return SCIP_OKAY;
}

/** updates the right hand sides of the KEP problem based on an attack pattern */
static
SCIP_RETCODE SCIPupdateRhsKEPProblem(
   SCIP*                 kepscip,            /**< SCIP instance of the KEP */
   SCIP_PROBDATA*        kepdata,            /**< problem data */
   int*                  attackpattern,      /**< array encoding attack pattern */
   int                   nattacks            /**< number of attacks in attackpattern */
   )
{
   SCIP_CONS** nodeconss;
   int nnodes;
   int v;

   assert ( kepscip != NULL );
   assert ( kepdata != NULL );
   assert ( attackpattern != NULL );

   nnodes = SCIPKEPdataGetNumNodes(kepdata);
   nodeconss = SCIPKEPdataGetNodeconss(kepdata);

   for (v = 0; v < nnodes; ++v)
   {
      SCIP_CALL( SCIPchgRhsLinear(kepscip, nodeconss[v], 1.0) );
   }

   for (v = 0; v < nattacks; ++v)
   {
      SCIP_CALL ( SCIPchgRhsLinear(kepscip, nodeconss[attackpattern[v]], 0.0) );
   }

   return SCIP_OKAY;
}


/** updates the right hand sides of the KEP problem based on an attack pattern */
static
SCIP_RETCODE SCIPenforceUnaffectedCyclesChains(
   SCIP*                 kepscip,            /**< SCIP instance of the KEP */
   SCIP*                 masterscip,         /**< SCIP instance of master problem */
   SCIP_SOL*             sol,                /**< SCIP solution to master problem */
   SCIP_PROBDATA*        kepdata,            /**< problem data */
   Cycles*               cycles,             /**< cycles data structure */
   Chains*               chains,             /**< chains data structure */
   int*                  attackpattern,      /**< array encoding attack pattern */
   int                   nattacks            /**< number of attacks in attackpattern */
   )
{
   SCIP_VAR** initxvars;
   SCIP_VAR** kepcyclevars;
   SCIP_VAR** kepchainvars;
   // SCIP_Bool attacked;
   int ncycles;
   int nchains;
   int c;
   // int i;
   // int v;

   assert( kepscip != NULL );
   assert( masterscip != NULL );
   assert( kepdata != NULL );
   assert( attackpattern != NULL );

   initxvars = masterProblemGetXvarinit(SCIPgetProbData(masterscip));
   kepcyclevars = SCIPKEPdataGetCyclevars(kepdata);
   kepchainvars = SCIPKEPdataGetChainvars(kepdata);

   ncycles = SCIPKEPdataGetNCycles(kepdata);
   nchains = SCIPKEPdataGetNChains(kepdata);

   for (c = 0; c < ncycles; ++c)
   {
      SCIP_CALL( SCIPchgVarLb(kepscip, kepcyclevars[c], 0.0) );

      if( SCIPgetSolVal(masterscip, sol, initxvars[c]) > 0.5 && !SCIPcycleIsAttacked(attackpattern, nattacks, cycles, c) )
         SCIP_CALL( SCIPchgVarLb(kepscip, kepcyclevars[c], 1.0) );

   }

   for (c = 0; c < nchains; ++c)
   {
      SCIP_CALL( SCIPchgVarLb(kepscip, kepchainvars[c], 0.0) );
      if( SCIPgetSolVal(masterscip, sol, initxvars[ncycles + c]) > 0.5 && !SCIPchainIsAttacked(attackpattern, nattacks, chains, c) )
         SCIP_CALL( SCIPchgVarLb(kepscip, kepchainvars[c], 1.0) );
   }

   return SCIP_OKAY;
}

/** updates the right hand sides of the KEP problem based on an attack pattern */
static
SCIP_RETCODE SCIPupdateObjCoefsKEPProblem(
   SCIP*                 kepscip,            /**< SCIP instance of the KEP */
   SCIP_PROBDATA*        kepdata,            /**< problem data */
   Cycles*               cycles,             /**< cycles in graph */
   Chains*               chains,             /**< chains in graph */
   int*                  attackpattern,      /**< array encoding attack pattern */
   int                   nattacks            /**< number of attacks in attackpattern */
   )
{
   SCIP_VAR** cyclevars;
   SCIP_VAR** chainvars;
   int ncycles;
   int nchains;
   int nnodes;
   int c;
   SCIP_Bool attacked;

   assert ( kepscip != NULL );
   assert ( kepdata != NULL );
   assert ( attackpattern != NULL );

   cyclevars = SCIPKEPdataGetCyclevars(kepdata);
   ncycles = SCIPKEPdataGetNCycles(kepdata);

   chainvars = SCIPKEPdataGetChainvars(kepdata);
   nchains = SCIPKEPdataGetNChains(kepdata);

   nnodes = SCIPKEPdataGetNumNodes(kepdata);

   for (c = 0; c < ncycles; ++c)
   {
      attacked = SCIPcycleIsAttacked(attackpattern, nattacks, cycles, c);
      if ( attacked )
      {
         SCIP_CALL( SCIPchgVarObj(kepscip, cyclevars[c], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj(kepscip, cyclevars[c], cycles->cycleweights[c] * nnodes +1.0) );
      }
   }

   for (c = 0; c < nchains; ++c)
   {
      attacked = SCIPchainIsAttacked(attackpattern, nattacks, chains, c);

      if ( attacked )
      {
         SCIP_CALL( SCIPchgVarObj( kepscip, chainvars[c], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj( kepscip, chainvars[c], chains->chainweights[c] * nnodes + 1.0) );
      }
   }

   return SCIP_OKAY;
}

/** method used to iteratively solve the Benders model and the corresponding KEPs */
static
SCIP_RETCODE SCIPsolveBendersModel(
   SCIP*                 bendersscip,        /**< SCIP data structure for Benders model */
   SCIP*                 kepscip,            /**< SCIP data structure for the subproblem of Benders (a KEP) */
   SCIP*                 masterscip,         /**< SCIP data structure for the masterproblem */
   SCIP_SOL*             mastersol,          /**< solution of the master problem */
   SCIP_Real             masterobj,          /**< objective value of the master problem */
   Cycles*               cycles,             /**< data structure needed for cycle weights when adding Benders cuts */
   Chains*               chains,             /**< data structure needed for chain weights when adding Benders cuts */
   SCIP_Real             timelimit,          /**< time limit for solving the Benders model */
   SCIP_Bool*            didnotfinish,       /**< pointer to store whether we have hit the time limit */
   SCIP_Bool*            optimal,            /**< pointer to store whether the method terminated optimally */
   SCIP_Real*            timestage2,         /**< pointer to store time spent in stage 2 */
   SCIP_Real*            timestage3,         /**< pointer to store time spent in stage 3 */
   SCIP_Real*            kepobj              /**< pointer to store value of final iteration kepscip solution */
   )
{
   int* attackpattern;
   int c;
   int current_c;
   int subchainidx;
   int nattacks;
   int nnodes;
   int adversarybound;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool liftsols;
   int maxnvarsinsolcons;
   int ncycles;
   int nchains;
   SCIP_Real bendersobj;
   // SCIP_Real kepobj;
   SCIP_SOL* sol;
   SCIP_SOL* kepsol;
   SCIP_VAR** cyclevars;
   SCIP_VAR** chainvars;
   SCIP_VAR** benderscyclevars;
   SCIP_VAR** benderschainvars;
   SCIP_PROBDATA* bendersdata;
   SCIP_PROBDATA* kepdata;
   SCIP_Real memlimit;
   SCIP_Real begintime;
   SCIP_Real endtime;
   SCIP_Real smallestub;
   int cnt = 0;
   int policy;

   assert( bendersscip != NULL );
   assert( kepscip != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( didnotfinish != NULL );
   assert( optimal != NULL );
   assert( timestage2 != NULL );
   assert( timestage3 != NULL );

   begintime = clock();
   *didnotfinish = FALSE;
   *optimal = FALSE;
   *timestage3 = 0.0;
   smallestub = SCIPinfinity(bendersscip);

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/liftbenderscuts", &liftsols) );
   bendersdata = SCIPgetProbData(bendersscip);
   kepdata = SCIPgetProbData(kepscip);

   assert( bendersdata != NULL );
   assert( kepdata != NULL );

   adversarybound = SCIPbendersdataGetAdversaryBound(bendersdata);
   ncycles = SCIPbendersdataGetNCycles(bendersdata);
   nchains = SCIPbendersdataGetNChains(bendersdata);
   maxnvarsinsolcons = SCIPKEPdataGetNumNodes(kepdata);
   nnodes = maxnvarsinsolcons;

   SCIP_CALL( SCIPgetRealParam(masterscip, "limits/memory", &memlimit) );
   SCIP_CALL( SCIPgetIntParam(masterscip, "kidney/recoursepolicy", &policy) );

   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &attackpattern, adversarybound) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &vars, maxnvarsinsolcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &vals, maxnvarsinsolcons) );

   /* Stronger bound on objective variable, as it is monotone increasing with respect to future iterations */
   SCIP_CALL( SCIPchgVarUb(bendersscip, SCIPbendersdataGetObjvar(bendersdata), masterobj) );

   /* In this loop we have to make sure that the transformed bendersscip and kepscip are not freed.
    * Otherwise, we cannot access the termination status in solveMasterProblem().
    */
   while ( TRUE )
   {
      int loop;
      SCIP_Real begintimestage3;
      SCIP_Real endtimestage3;

      /* solve Benders model */
      SCIP_CALL( SCIPsetLimitsAndVerbose(bendersscip, begintime, timelimit,
            MAX(1, memlimit - SCIPgetMemUsed(masterscip) / 1048576), TRUE) );
      SCIP_CALL( SCIPsolve(bendersscip) );

      if ( solvingloopShallTerminate(bendersscip, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         break;
      }

      begintimestage3 = clock();

      sol = SCIPgetBestSol(bendersscip);
      bendersobj = SCIPgetSolOrigObj(bendersscip, sol);
      SCIPinfoMessage(bendersscip, NULL, "Benders obj: %f\n", bendersobj);

      if ( SCIPisGE(bendersscip, bendersobj, masterobj - 0.5) )
      {
         SCIPinfoMessage(bendersscip, NULL, "@95 No attack pattern violates the master solution.\n");
         *optimal = TRUE;
         endtimestage3 = clock();
         *timestage3 += (SCIP_Real) (endtimestage3 - begintimestage3) / CLOCKS_PER_SEC;
         break;
      }

      SCIP_CALL( SCIPgetAttackPattern(bendersscip, sol, attackpattern, &nattacks,
            SCIPbendersdataGetNumNodes(bendersdata), METHOD_BENDERS) );
      SCIPinfoMessage(bendersscip, NULL, "Attack pattern: \n");

      for (loop = 0; loop < nattacks; ++loop)
      {
         SCIPinfoMessage(bendersscip, NULL, "%d ", attackpattern[loop]);
      }
      SCIPinfoMessage(bendersscip, NULL, "\n");

      /* make sure that we can modify the KEP SCIP instance */
      SCIP_CALL( SCIPfreeTransform(kepscip) );

      if( policy == POLICY_KEEPUNAFFECTEDCC )
         SCIP_CALL( SCIPenforceUnaffectedCyclesChains(kepscip, masterscip, mastersol, kepdata, cycles, chains, attackpattern, nattacks) );

      if ( ! liftsols )
         SCIP_CALL( SCIPupdateRhsKEPProblem(kepscip, kepdata, attackpattern, nattacks) );
      else
         SCIP_CALL( SCIPupdateObjCoefsKEPProblem(kepscip, kepdata, cycles, chains, attackpattern, nattacks) );

      SCIP_CALL( SCIPsetLimitsAndVerbose(kepscip, begintime, timelimit,
            MAX(1, memlimit - SCIPgetMemUsed(masterscip) / 1048576 - SCIPgetMemUsed(bendersscip) / 1048576), TRUE) );
      SCIP_CALL( SCIPsolve(kepscip) );

      if ( solvingloopShallTerminate(kepscip, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         endtimestage3 = clock();
         *timestage3 += (SCIP_Real) (endtimestage3 - begintimestage3) / CLOCKS_PER_SEC;
         break;
      }

      kepsol = SCIPgetBestSol(kepscip);

      /* Notice that KEPscip picks the solution among the highest weight ones with highest number of cycles / chains.
       * The actual objective value of the maximum weight is therefore different
       */

      /* Find out the solution in kepscip and add a Benders type cut to the benders model */
      benderscyclevars = SCIPbendersdataGetCyclevars(bendersdata);
      benderschainvars = SCIPbendersdataGetChainvars(bendersdata);

      cyclevars = SCIPKEPdataGetCyclevars(kepdata);
      chainvars = SCIPKEPdataGetChainvars(kepdata);

      *kepobj = 0;
      cnt = 0;

      for (c = 0; c < ncycles; ++c)
      {
         /* Only take into account cycles in kepsol that are not attacked by the current attack pattern */
         if ( SCIPgetSolVal(kepscip, kepsol, cyclevars[c]) > 0.5 )
         {
            if ( SCIPvarGetObj(cyclevars[c]) > nnodes )
               *kepobj += cycles->cycleweights[c];

            vars[cnt] = benderscyclevars[c];
            vals[cnt++] = cycles->cycleweights[c];
         }
      }

      for (c = 0; c < nchains; ++c)
      {
         /* Only take into account chains in kepsol that are not attacked by the current attack pattern */
         if ( SCIPgetSolVal(kepscip, kepsol, chainvars[c]) > 0.5 )
         {
            if ( SCIPvarGetObj(chainvars[c]) > nnodes )
               *kepobj += chains->chainweights[c];

            current_c = c;
            subchainidx = chains->subchains[current_c];
            while( TRUE )
            {
               vars[cnt] = benderschainvars[current_c];
               if( subchainidx == -1 )
               {
                  vals[cnt++] = chains->chainweights[current_c];
                  break;
               }
               else
                  vals[cnt++] = chains->chainweights[current_c] - chains->chainweights[subchainidx];
               current_c = subchainidx;
               subchainidx = chains->subchains[current_c];
            }
         }
      }
      SCIPinfoMessage(bendersscip, NULL, "@95 KEP objective: %f\n", *kepobj);

      if ( SCIPisLT(kepscip, *kepobj, masterobj - 0.5) )
      {
         SCIPinfoMessage(kepscip, NULL, "The recourse value %f of the currently considered attack pattern violates the master solution.\n", *kepobj);
         *optimal = TRUE;
         endtimestage3 = clock();
         *timestage3 += (SCIP_Real) (endtimestage3 - begintimestage3) / CLOCKS_PER_SEC;
         break;
      }

      /* Free the Benders model before we can add new cuts */
      SCIP_CALL( SCIPfreeTransform(bendersscip) );

      for (c = 0; c < cnt; ++c)
      {
         SCIP_CALL( SCIPchgVarLb(bendersscip, vars[c], 0.0) );
      }

      /* Also add dummy variable to vars of the cut */
      vals[cnt] = -1.0;
      vars[cnt++] = SCIPbendersdataGetObjvar(bendersdata);

      SCIP_CALL( SCIPbendersdataAddSolCons(bendersscip, vals, vars, cnt) );

      // SCIP_CALL( SCIPwriteOrigProblem(bendersscip, "bendersscip.lp", NULL, FALSE) );

      /* Stronger bound on objective variable, as it is monotone increasing with respect to future iterations */
      SCIP_CALL( SCIPchgVarLb(bendersscip, SCIPbendersdataGetObjvar(bendersdata), bendersobj) );
      if ( *kepobj < smallestub )
      {
         SCIP_CALL( SCIPchgVarUb(bendersscip, SCIPbendersdataGetObjvar(bendersdata), *kepobj) );
         smallestub = *kepobj;
      }

      endtimestage3 = clock();
      *timestage3 += (SCIP_Real) (endtimestage3 - begintimestage3) / CLOCKS_PER_SEC;
   }

   SCIPfreeBlockMemoryArrayNull(bendersscip, &attackpattern, adversarybound);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &vals, maxnvarsinsolcons);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &vars, maxnvarsinsolcons);

   endtime = clock();
   *timestage2 = (SCIP_Real) (endtime - begintime) / CLOCKS_PER_SEC - *timestage3;

   return SCIP_OKAY;
}

/** gets attack pattern from a 2/3 stage solution */
SCIP_RETCODE SCIPgetAttackPattern(
   SCIP*                 scip,               /**< SCIP instance of 2/3 stage */
   SCIP_SOL*             sol,                /**< solution of 2/3 stage */
   int*                  attackpattern,      /**< array to store attack pattern */
   int*                  nattacks,           /**< pointer to store number of attacks */
   int                   nnodes,             /**< number of nodes in underlying graph */
   int                   method              /**< Which method is used for stage 2/3? (1:Benders, 2:B&P) */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VAR** nodevars;
   int v;
   int cnt = 0;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks != NULL );

   *nattacks = 0;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   if ( method == METHOD_BENDERS )
      nodevars = SCIPbendersdataGetUVars(probdata);
   else if ( method == METHOD_BENDERS_PICEF )
   {
      assert( method == METHOD_BENDERS_PICEF );

      nodevars = SCIPbendersdataPICEFGetUVars(probdata);
   }
   assert( nodevars != NULL );

   /* check which nodes are attacked */
   for (v = 0; v < nnodes; ++v)
   {
      if ( SCIPgetSolVal(scip, sol, nodevars[v]) < 0.5 )
         attackpattern[cnt++] = v;
   }
   *nattacks = cnt;

   return SCIP_OKAY;
}

/** solves the model of Glorie et al. */
static
SCIP_RETCODE SCIPsolveGlorieEtAlModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< pointer to underlying cycle structure */
   Chains*               chains,             /**< pointer to underlying chain structure */
   SCIP_Real             timelimit,          /**< time limit to solve problem */
   SCIP_Real             memlimit,           /**< memory limit */
   SCIP_Real             masterobjbound,     /**< bound on the objective value given by master problem */
   int                   attackbound,        /**< upper bound on number of attacks */
   int*                  usedcycleschains,   /**< array containing indices of cycles and chains used in master solution */
   int                   nusedcycleschains,  /**< number of cycles and chains used in master solution */
   int*                  attacksolu,         /**< array to store optimal attack pattern */
   int*                  nattacksolu,        /**< pointer to store number of attacks in optimal solution */
   SCIP_Real*            optobj,             /**< pointer to store objective at termination */
   SCIP_Bool*            didnotfinish,       /**< pointer to store whether problem could not be solved within timelimit */
   SCIP_Bool*            optimal,            /**< pointer to store whether the method terminated optimally */
   SCIP_Real*            timeneeded          /**< pointer to store time needed to solve problem */
   )
{
   BBNode** problemstack;
   SCIP_PROBDATA* probdata;
   int* attackpattern;
   int maxnproblems;
   int i;
   int nproblems = 0;
   SCIP_Bool terminated = FALSE;
   SCIP_Real begintime;
   SCIP_Real endtime;

   assert( scip != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( usedcycleschains != NULL );
   assert( attacksolu != NULL );
   assert( nattacksolu != NULL );
   assert( optobj != NULL );
   assert( didnotfinish != NULL );
   assert( timeneeded != NULL );

   begintime = clock();

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   *optobj = SCIPinfinity(scip);
   *optimal = FALSE;

   maxnproblems = 100;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problemstack, maxnproblems) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &attackpattern, attackbound) );

   /* create initial problem */
   SCIP_CALL( createBBNode(scip, &problemstack[nproblems++], 0, 0, 0, -1.0) );

   /* solve all problems on the problem stack */
   while ( nproblems > 0 )
   {
      SCIP_VAR** cyclevars;
      SCIP_VAR** chainvars;
      BBNode* curproblem;
      SCIP_Real tmpub;
      SCIP_Real lb;
      int nattacks;
      int nattacksused;
      int newattack;
      int c;

      curproblem = problemstack[--nproblems];
      assert( curproblem != NULL );

      newattack = curproblem->branchingcand;

      nattacks = curproblem->nzerofixings;
      for (i = 0; i < nattacks; ++i)
         attackpattern[i] = curproblem->zerofixings[i];
      nattacksused = nattacks;

      /* possibly extend the attack */
      if ( nattacksused < attackbound )
      {
         lb = GlorieExtendAttack(scip, graph, cycles, chains, usedcycleschains, nusedcycleschains,
            curproblem->zerofixings, curproblem->onefixings, curproblem->nzerofixings, curproblem->nonefixings,
            attackbound, attackpattern);

         SCIPinfoMessage(scip, NULL, "Attack Pattern: ");
         for (i = 0; i < attackbound; ++i)
            SCIPinfoMessage(scip, NULL, "%d \t", attackpattern[i]);
         SCIPinfoMessage(scip, NULL, "\n");
         nattacksused = attackbound;
         newattack = attackpattern[curproblem->nzerofixings];
      }

      /* incorporate branching decisions into problem */
      cyclevars = getCyclevarsGlorie(probdata);
      chainvars = getChainvarsGlorie(probdata);

      for (i = 0; i < nattacksused; ++i)
      {
         assert( 0 <= attackpattern[i] && attackpattern[i] < graph->nnodes );

         for (c = graph->node2cyclesbegin[attackpattern[i]]; c < graph->node2cyclesbegin[attackpattern[i] + 1]; ++c)
         {
            SCIP_CALL( SCIPchgVarUb(scip, cyclevars[graph->node2cycles[c]], 0.0) );
         }
         for (c = graph->node2chainsbegin[attackpattern[i]]; c < graph->node2chainsbegin[attackpattern[i] + 1]; ++c)
         {
            SCIP_CALL( SCIPchgVarUb(scip, chainvars[graph->node2chains[c]], 0.0) );
         }
      }

      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         SCIP_CALL( freeBBNode(scip, &curproblem) );
         break;
      }

      /* compute upper bound if it is not known from parent node */
      tmpub = getObjbound(curproblem);
      if ( tmpub < -0.5 )
      {
         SCIP_CALL( SCIPsetLimitsAndVerbose(scip, begintime, timelimit, memlimit, TRUE) );
         SCIP_CALL( SCIPsolve(scip) );
      }

      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         SCIP_CALL( freeBBNode(scip, &curproblem) );
         break;
      }

      if ( tmpub < -0.5 )
         tmpub = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
      assert( tmpub > -0.5 );

      /* possibly update upper bound */
      if ( SCIPisLT(scip, tmpub, *optobj) )
      {
         SCIPinfoMessage(scip, NULL, "update upper bound from %f to %f\n", *optobj, tmpub);
         *optobj = tmpub;

         *nattacksolu = nattacksused;
         for (i = 0; i < nattacksused; ++i)
            attacksolu[i] = attackpattern[i];
      }

      /* optimality check */
      if ( SCIPisLT(scip, *optobj, masterobjbound) )
      {
         SCIPinfoMessage(scip, NULL, "@95 optimal solution found: ub is %f and masterbound is %f\n", *optobj, masterobjbound);

         SCIP_CALL( freeBBNode(scip, &curproblem) );
         terminated = TRUE;

         goto UNDOBRANCHING;
      }

      /* if we have hit the timelimit */
      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         terminated = TRUE;
         *didnotfinish = TRUE;
         SCIP_CALL( freeBBNode(scip, &curproblem) );
         goto UNDOBRANCHING;
      }

      /* pruning */
      if ( nattacks == attackbound || SCIPisGE(scip, lb, *optobj) )
      {
         SCIPinfoMessage(scip, NULL, "@95 prune node: %d / %d or (lb %f >= ub %f)\n", nattacks, attackbound, lb, *optobj);

         SCIP_CALL( freeBBNode(scip, &curproblem) );
         goto UNDOBRANCHING;
      }

      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         SCIP_CALL( freeBBNode(scip, &curproblem) );
         break;
      }

      /* branching */
      if ( nproblems >= maxnproblems - 1 )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, maxnproblems + 1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problemstack, maxnproblems, newsize) );
         maxnproblems = newsize;
      }

      /* if it is possible to attack a further node */
      SCIPinfoMessage(scip, NULL, "@95 branch on node %d\n", newattack);
      if ( newattack < graph->nnodes && curproblem->nzerofixings < attackbound && attackbound + curproblem->nonefixings < graph->nnodes )
      {
         SCIP_CALL( createBBNode(scip, &problemstack[nproblems], curproblem->nzerofixings, curproblem->nonefixings + 1, newattack + 1, tmpub) );
         for (i = 0; i < curproblem->nzerofixings; ++i)
         {
            SCIP_CALL( addZerofixing(scip, problemstack[nproblems], curproblem->zerofixings[i], FALSE) );
         }
         for (i = 0; i < curproblem->nonefixings; ++i)
         {
            SCIP_CALL( addOnefixing(scip, problemstack[nproblems], curproblem->onefixings[i], FALSE) );
         }
         SCIP_CALL( addOnefixing(scip, problemstack[nproblems++], newattack, FALSE) );

         SCIP_CALL( createBBNode(scip, &problemstack[nproblems], curproblem->nzerofixings + 1, curproblem->nonefixings, newattack + 1, -1.0) );
         for (i = 0; i < curproblem->nzerofixings; ++i)
         {
            SCIP_CALL( addZerofixing(scip, problemstack[nproblems], curproblem->zerofixings[i], FALSE) );
         }
         for (i = 0; i < curproblem->nonefixings; ++i)
         {
            SCIP_CALL( addOnefixing(scip, problemstack[nproblems], curproblem->onefixings[i], FALSE) );
         }
         SCIP_CALL( addZerofixing(scip, problemstack[nproblems++], newattack, FALSE) );
      }

      SCIP_CALL( freeBBNode(scip, &curproblem) );

      /* undo 0-branchings */
   UNDOBRANCHING:

      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         break;
      }

      if ( terminated || nproblems == 0 )
      {
         if ( ! *didnotfinish )
            *optimal = TRUE;

         break;
      }

      SCIP_CALL( SCIPfreeTransform(scip) );
      for (i = 0; i < attackbound; ++i)
      {
         assert( 0 <= attackpattern[i] && attackpattern[i] < graph->nnodes );

         for (c = graph->node2cyclesbegin[attackpattern[i]]; c < graph->node2cyclesbegin[attackpattern[i] + 1]; ++c)
         {
            SCIP_CALL( SCIPchgVarUb(scip, cyclevars[graph->node2cycles[c]], 1.0) );
         }
         for (c = graph->node2chainsbegin[attackpattern[i]]; c < graph->node2chainsbegin[attackpattern[i] + 1]; ++c)
         {
            SCIP_CALL( SCIPchgVarUb(scip, chainvars[graph->node2chains[c]], 1.0) );
         }
      }

      if ( solvingloopShallTerminate(NULL, begintime, timelimit) )
      {
         *didnotfinish = TRUE;
         break;
      }
   }

   for (i = nproblems - 1; i >= 0; --i)
   {
      SCIP_CALL( freeBBNode(scip, &problemstack[i]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &attackpattern, attackbound);
   SCIPfreeBlockMemoryArrayNull(scip, &problemstack, maxnproblems);

   endtime = clock();
   *timeneeded = (SCIP_Real) (endtime - begintime) / CLOCKS_PER_SEC;

   return SCIP_OKAY;
}

/** add cut based on initial master solution to bendersscip */
static
SCIP_RETCODE SCIPaddInitialBendersCut(
   SCIP*                 masterscip,         /**< SCIP instance of master problem */
   SCIP*                 bendersscip,        /**< SCIP instance of Benders problem */
   SCIP_SOL*             sol                 /**< solution of master SCIP */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_PROBDATA* bendersdata;
   SCIP_VAR** benderscyclevars;
   SCIP_VAR** benderschainvars;
   SCIP_VAR** initcycleschains;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   Cycles* cycles;
   Chains* chains;
   int nnodes;
   int ncycles;
   int nchains;
   int cnt;
   int c;
   int current_c; // copy chain index
   int subchainidx;

   assert( masterscip != NULL );
   assert( bendersscip != NULL );

   probdata = SCIPgetProbData(masterscip);
   assert( probdata != NULL );

   bendersdata = SCIPgetProbData(bendersscip);
   assert( bendersdata != NULL );

   benderscyclevars = SCIPbendersdataGetCyclevars(bendersdata);
   benderschainvars = SCIPbendersdataGetChainvars(bendersdata);
   initcycleschains = masterProblemGetXvarinit(probdata);

   cycles = masterProblemGetCycles(probdata);
   chains = masterProblemGetChains(probdata);

   nnodes = SCIPbendersdataGetNumNodes(bendersdata);
   ncycles = SCIPbendersdataGetNCycles(bendersdata);
   nchains = SCIPbendersdataGetNChains(bendersdata);

   SCIP_CALL( SCIPallocBufferArray(masterscip, &vars, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &vals, nnodes + 1) );

   cnt = 0;
   for (c = 0; c < ncycles; ++c)
   {
      /* Only take into account cycles in kepsol that are not attacked by the current attack pattern */
      if ( SCIPgetSolVal(masterscip, sol, initcycleschains[c]) > 0.5 )
      {
         vars[cnt] = benderscyclevars[c];
         vals[cnt++] = cycles->cycleweights[c];
      }
   }

   for (c = 0; c < nchains; ++c)
   {
      /* Only take into account cycles in kepsol that are not attacked by the current attack pattern */
      if ( SCIPgetSolVal(masterscip, sol, initcycleschains[cycles->ncycles + c]) > 0.5 )
      {
         current_c = c;
         subchainidx = chains->subchains[current_c];
         while( TRUE )
         {
            vars[cnt] = benderschainvars[current_c];
            if( subchainidx == -1 )
            {
               vals[cnt++] = chains->chainweights[current_c];
               break;
            }
            else
               vals[cnt++] = chains->chainweights[current_c] - chains->chainweights[subchainidx];
            current_c = subchainidx;
            subchainidx = chains->subchains[current_c];
         }
      }
   }

   /* also add dummy variable to vars of the cut */
   vals[cnt] = -1.0;
   vars[cnt++] = SCIPbendersdataGetObjvar(bendersdata);

   /* the vars are being enabled through resetting the lower bound */
   for (c = 0; c < cnt; ++c)
   {
      SCIP_CALL( SCIPchgVarLb(bendersscip, vars[c], 0.0) );
   }

   /* add cut based on initial master solution to bendersscip */
   SCIP_CALL( SCIPbendersdataAddSolCons(bendersscip, vals, vars, cnt) );

   SCIPfreeBufferArray(masterscip, &vals);
   SCIPfreeBufferArray(masterscip, &vars);

   return SCIP_OKAY;
}

/** solves the master problem */
SCIP_RETCODE solveMasterProblem(
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   Triplets*             triplets,           /**< pointer to triple structures of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Real             timelimit,          /**< time limit to solve the problem */
   SCIP_Real             initialtime,        /**< time spent in initialization */
   const char*           settings,           /**< possible name of setting file */
   SCIP_Bool             verbose             /**< whether we print SCIP's logs */
   )
{
   /** General part of master problem solve method */
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   int* attackpattern;
   int* scenarios;
   SCIP_Real masterobj;
   SCIP_Real subobj;
   SCIP_Real memlimit;
   SCIP_Real modtimelimit;
   SCIP_Real newtimelimit;
   SCIP_Real totaltime;
   SCIP_Real kepobj;
   clock_t begintime;
   clock_t endtime;
   int nscenarios;
   int nattacks;
   int i;
   int j;
   int pos;
   int method;
   SCIP_Bool didnotfinish = FALSE;
   int cnt = 0;
   SCIP_Real timesecondstage = 0.0;
   SCIP_Real timethirdstage = 0.0;

   assert( masterscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( triplets != NULL );
   assert( adversarybound >= 0 );

   /* start time measurement */
   begintime = clock();

   SCIP_CALL( SCIPgetIntParam(masterscip, "kidney/method", &method) );
   SCIP_CALL( SCIPgetRealParam(masterscip, "limits/memory", &memlimit) );

   SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &attackpattern, adversarybound) );

   SCIPinfoMessage(masterscip, NULL, "@99 Start master loop of trilevel optimization algorithm.\n");

   totaltime = initialtime;
   modtimelimit = timelimit - initialtime;

   /* the main loop to solve the problem */
   while ( TRUE )
   {
      SCIP* scip;
      SCIP* subscip;
      SCIP_Bool optimal;

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      SCIP_CALL( SCIPsetLimitsAndVerbose(masterscip, begintime, modtimelimit, memlimit, verbose) );
      SCIP_CALL( SCIPsetIntParam(masterscip, "display/verblevel", 3) );
      SCIP_CALL( SCIPsolve(masterscip) );

      if ( solvingloopShallTerminate(masterscip, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      /* obtain best solution and best objective */
      sol = SCIPgetBestSol(masterscip);
      masterobj = SCIPgetSolOrigObj(masterscip, sol);

      SCIPinfoMessage(masterscip, NULL, "@98 masteriteration %d: objective value %f.\n", ++cnt, masterobj);

      /* adapt the weights of cycles and chains based on initial solution of master problem */
      SCIP_CALL( adaptCycleChainWeights(masterscip, sol, graph, cycles, chains) );

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      /* create the SCIP instance based on selected method */
      if ( method == METHOD_BENDERS )
      {
         SCIP_Real beginsubtime;
         SCIP_Real endsubtime;

         beginsubtime = clock();

         SCIP_CALL( SCIPcreateBasicInstance(masterscip, &scip, verbose) );
         SCIP_CALL( SCIPcreateBendersModel(scip, masterscip, sol, graph, cycles, chains, adversarybound) );
         SCIP_CALL( SCIPsetObjIntegral(scip) );

         endsubtime = clock();
         timesecondstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;

         beginsubtime = clock();

         SCIP_CALL( SCIPcreateBasicInstance(masterscip, &subscip, verbose) );
         SCIP_CALL( SCIPcreateBendersSubModel(subscip, graph, cycles, chains) );
         SCIP_CALL( SCIPsetObjIntegral(subscip) );
         SCIPinfoMessage(masterscip, NULL, "@97 Set up Benders model %d.\n", cnt);

         endsubtime = clock();
         timethirdstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;

         /* add cut based on initial master solution to bendersscip */
         SCIP_CALL( SCIPaddInitialBendersCut(masterscip, scip, sol) );
      }
      else
      {
         SCIP_Real beginsubtime;
         SCIP_Real endsubtime;

         assert( method == METHOD_BRANCHANDBOUND );

         beginsubtime = clock();

         SCIP_CALL( SCIPcreateBasicInstance(masterscip, &scip, verbose) );
         SCIP_CALL( SCIPcreateGlorieEtAlModel(scip, graph, cycles, chains) );
         SCIP_CALL( SCIPsetObjIntegral(scip) );
         SCIPinfoMessage(masterscip, NULL, "@97 Set up Glorie model %d.\n", cnt);

         endsubtime = clock();
         timesecondstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;
      }

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREESUBPROBLEM;

      /* solve stage 2/3 */
      if ( method == METHOD_BENDERS )
      {
         clock_t curtime;
         SCIP_Real timestage2;
         SCIP_Real timestage3;

         curtime = clock();
         newtimelimit = modtimelimit - (SCIP_Real) (curtime - begintime) / CLOCKS_PER_SEC;

         SCIP_CALL( SCIPsolveBendersModel(scip, subscip, masterscip, sol, masterobj, cycles, chains,
               newtimelimit, &didnotfinish, &optimal, &timestage2, &timestage3, &kepobj) );

         timesecondstage += timestage2;
         timethirdstage += timestage3;

         if ( optimal )
         {
            subobj = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
            SCIPinfoMessage(masterscip, NULL, "@96 Benders problem iteration %d: objective value %f.\n", cnt, subobj);
         }
      }
      else
      {
         clock_t curtime;
         SCIP_VAR** xvars;
         int* usedcycleschains;
         int nusedcycleschains = 0;
         SCIP_Real timeneeded;

         assert( method == METHOD_BRANCHANDBOUND );

         xvars = masterProblemGetXvarinit(SCIPgetProbData(masterscip));

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &usedcycleschains, graph->nnodes) );
         for (i = 0; i < cycles->ncycles + chains->nchains; ++i)
         {
            if ( SCIPgetSolVal(masterscip, sol, xvars[i]) > 0.5 )
               usedcycleschains[nusedcycleschains++] = i;
         }

         curtime = clock();
         newtimelimit = modtimelimit - (SCIP_Real) (curtime - begintime) / CLOCKS_PER_SEC;

         SCIP_CALL( SCIPsolveGlorieEtAlModel(scip, graph, cycles, chains, newtimelimit,
               MAX(1, memlimit - SCIPgetMemUsed(masterscip) / 1048576), masterobj,  adversarybound,
               usedcycleschains, nusedcycleschains, attackpattern, &nattacks, &subobj, &didnotfinish, &optimal,
               &timeneeded) );

         timesecondstage += timeneeded;

         SCIPfreeBlockMemoryArray(scip, &usedcycleschains, graph->nnodes);
      }

      /* if we have found an optimal solution */
      if ( optimal && SCIPisLE(masterscip, masterobj, subobj) )
      {
         SCIP_CALL( SCIPprintSol(masterscip, SCIPgetBestSol(masterscip), NULL, FALSE) );

         SCIP_CALL( SCIPfreeTransform(scip) );
         SCIP_CALL( SCIPfree(&scip) );
         if ( method == METHOD_BENDERS )
         {
            SCIP_CALL( SCIPfreeTransform(subscip) );
            SCIP_CALL( SCIPfree(&subscip) );
         }
         SCIP_CALL( SCIPfreeTransform(masterscip) );

         break;
      }
      /* if we hit the time limit or could not solve a subSCIP instance */
      else if ( didnotfinish || solvingloopShallTerminate(scip, begintime, modtimelimit)
         || ( method == METHOD_BENDERS && solvingloopShallTerminate(subscip, begintime, modtimelimit) ) )
         goto FREESUBPROBLEM;

      /* extract attack pattern and update master problem */
      SCIP_CALL( SCIPfreeTransform(masterscip) );

      if ( method != METHOD_BRANCHANDBOUND )
      {
         SCIP_Real beginsubtime;
         SCIP_Real endsubtime;

         beginsubtime = clock();

         SCIP_CALL( SCIPgetAttackPattern(scip, SCIPgetBestSol(scip), attackpattern, &nattacks, graph->nnodes, method) );

         SCIP_CALL( SCIPchgVarUb(masterscip, masterProblemGetObjvar(SCIPgetProbData(masterscip)), masterobj) );

         if( SCIPisLT(masterscip, kepobj, masterobj) )
            SCIP_CALL( SCIPchgVarLb(masterscip, masterProblemGetObjvar(SCIPgetProbData(masterscip)), kepobj) );

         endsubtime = clock();
         timesecondstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;
      }
      SCIPinfoMessage(masterscip, NULL, "nattacks: %d\n", nattacks);
      SCIPinfoMessage(masterscip, NULL, "first attack: %d\n", attackpattern[0]);
      SCIP_CALL( SCIPupdateMasterProblem(masterscip, attackpattern, nattacks) );

      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIP_CALL( SCIPfree(&scip) );
      if ( method == METHOD_BENDERS )
      {
         SCIP_CALL( SCIPfreeTransform(subscip) );
         SCIP_CALL( SCIPfree(&subscip) );
      }

      /* gotos if we could not solve a problem or hit the time limit */
      if ( FALSE )
      {
      FREESUBPROBLEM:
         if ( method == METHOD_BENDERS )
         {
            SCIP_CALL( SCIPfreeTransform(subscip) );
            SCIP_CALL( SCIPfree(&subscip) );
         }
         SCIP_CALL( SCIPfreeTransform(scip) );
         SCIP_CALL( SCIPfree(&scip) );

      FREEMASTERPROBLEM:
         SCIP_CALL( SCIPfreeTransform(masterscip) );
         didnotfinish = TRUE;
         break;
      }
   }

   SCIPfreeBlockMemoryArrayNull(masterscip, &attackpattern, adversarybound);

   endtime = clock();
   totaltime += (SCIP_Real) (endtime - begintime) / CLOCKS_PER_SEC;

   /* print statistics */
   SCIPinfoMessage(masterscip, NULL, "@05================================================================================\n");
   SCIPinfoMessage(masterscip, NULL, "   ================================   STATISTICS   ================================\n");
   if ( didnotfinish )
   {
      SCIPinfoMessage(masterscip, NULL, "[suboptimal] terminated after %f seconds (stage 2: %f, stage 3: %f).\n",
         totaltime, timesecondstage, timethirdstage);
   }
   else
   {
      SCIPinfoMessage(masterscip, NULL, "[optimal] terminated after %f seconds (stage 2: %f, stage 3: %f).\n",
      totaltime, timesecondstage, timethirdstage);
   }

   SCIPinfoMessage(masterscip, NULL, "\n[patterns]\n", totaltime);

   /* print attack patterns */
   probdata = SCIPgetProbData(masterscip);
   assert( probdata != NULL );

   scenarios = masterProblemGetScenarios(probdata);
   nscenarios = masterProblemGetNScenarios(probdata);

   SCIPinfoMessage(masterscip, NULL, "Used attack patterns:\n");
   pos = 0;
   for (i = 0; i < nscenarios; ++i)
   {
      SCIPinfoMessage(masterscip, NULL, "[scenario] %4d:\t", i);
      for (j = 0; j < adversarybound; ++j)
      {
         SCIPinfoMessage(masterscip, NULL, "%8d", scenarios[pos++]);
      }
      SCIPinfoMessage(masterscip, NULL, "\n");
   }
   SCIPinfoMessage(masterscip, NULL, "@06================================================================================\n");

   return SCIP_OKAY;
}
