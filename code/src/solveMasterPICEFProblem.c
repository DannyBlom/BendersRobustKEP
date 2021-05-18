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

/**@file   solveMasterPICEFProblem.c
 * @brief  Methods to solve the master problem based on the PICEF formulation for kidney exchange programs
 * @author Danny Blom
 *
 * This file provides routines to solve the master problem with PICEF formulation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <time.h>

#include "auxiliaryStructures.h"
#include "graph.h"
#include "probdata_master_kidney_picef.h"
#include "probdata_benders_picef.h"
#include "probdata_benders_subproblem_picef.h"
#include "problem_kidneyexchange.h"
#include "problem_benders_picef.h"
#include "solveMasterProblem.h"
#include "solveMasterPICEFProblem.h"

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


/* adapt cycle and posarc weights based on solution of master problem */
static
SCIP_RETCODE adaptCyclePosarcWeights(
   SCIP*                 masterscip,         /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution of masterscip */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< underlying cycles structure */
   PositionedArcs*       posarcs             /**< underlying arcs structure */
   )
{
   SCIP_VAR** dummyyvars;
   SCIP_VAR* var;
   SCIP_PROBDATA* probdata;
   int ncycles;
   int nposarcs;
   int npairs;
   int c;
   int v;

   assert( masterscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );

   probdata = SCIPgetProbData(masterscip);
   assert( probdata != NULL );

   dummyyvars = masterPICEFProblemGetDummyYVars(probdata);
   assert( dummyyvars != NULL );

   ncycles = cycles->ncycles;
   nposarcs = posarcs->nposarcs;
   npairs = graph->npairs;

   /* Reset cycle weights */
   for (c = 0; c < ncycles; ++c)
      cycles->cycleweights[c] = 0;

   /* For every pair in the graph, check if it is involved in the initial solution
    * If so, increase weights of cycles / posarcs containing the pair by 1
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
   }

   /* Increment posarc weights */
   for (c = 0; c < nposarcs; ++c)
   {
      v = posarcs->nodelists[2*c+1];
      var = dummyyvars[v];
      if ( SCIPgetSolVal(masterscip, sol, var) > 0.5 )
         posarcs->arcweights[c] = 1;
      else
         posarcs->arcweights[c] = 0;
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

   nnodes = SCIPKEPdataPICEFGetNumNodes(kepdata);
   nodeconss = SCIPKEPdataPICEFGetNodeconss(kepdata);

   for (v = 0; v < nnodes; ++v)
   {
      SCIP_CALL( SCIPchgRhsLinear( kepscip, nodeconss[v], 1.0) );
   }

   for (v = 0; v < nattacks; ++v)
   {
      SCIP_CALL ( SCIPchgRhsLinear( kepscip, nodeconss[attackpattern[v]], 0.0) );
   }
   return SCIP_OKAY;
}

/** updates the right hand sides of the KEP problem based on an attack pattern */
static
SCIP_RETCODE SCIPupdateUbDummyVars(
   SCIP*                 kepscip,            /**< SCIP instance of the KEP */
   SCIP_PROBDATA*        kepdata,            /**< problem data */
   int*                  attackpattern,      /**< array encoding attack pattern */
   int                   nattacks            /**< number of attacks in attackpattern */
   )
{
   SCIP_VAR** dummyarcvars;
   Graph* graph;
   int nnodes;
   int i;
   int j;
   SCIP_Bool attacked;

   assert ( kepscip != NULL );
   assert ( kepdata != NULL );
   assert ( attackpattern != NULL );

   nnodes = SCIPKEPdataPICEFGetNumNodes(kepdata);
   dummyarcvars = SCIPKEPdataPICEFGetDummyArcvars(kepdata);
   graph = SCIPKEPdataPICEFGetGraph(kepdata);

   for (i = 0; i < nnodes; ++i)
   {
      for (j = graph->adjacencylistbegins[i]; j < graph->adjacencylistbegins[i+1]; ++j)
      {
         attacked = SCIPdummyarcIsAttacked(attackpattern, nattacks, graph, i, j);
         if ( attacked )
            SCIP_CALL( SCIPchgVarUb(kepscip, dummyarcvars[j], 0.0) );
         else
            SCIP_CALL( SCIPchgVarUb(kepscip, dummyarcvars[j], 1.0) );
      }
   }
   return SCIP_OKAY;
}

/** returns whether an arc is attacked based on node attack */
SCIP_Bool SCIParcIsAttacked(
   int*                  attackpattern,      /**< Array of attacked vertices */
   int                   nattacks,           /**< Number of attacked vertices */
   PositionedArcs*       posarcs,            /**< Data structure for the position indexed arcs in the graph */
   int                   arcidx              /**< Index of position indexed arc in arcs data structure */
   )
{
   int i;
   int vtx;

   for (i = 0; i < nattacks; ++i)
   {
      vtx = attackpattern[i];
      if ( vtx == posarcs->nodelists[2 * arcidx] || vtx == posarcs->nodelists[ 2 *arcidx + 1] )
         return TRUE;
   }

   return FALSE;
}

/** returns whether a dummy arc is attacked based on node attack */
SCIP_Bool SCIPdummyarcIsAttacked(
   int*                  attackpattern,      /**< Array of attacked vertices */
   int                   nattacks,           /**< Number of attacked vertices */
   Graph*                graph,              /**< Data structure for the position indexed arcs in the graph */
   int                   source,             /**< vertex that the arc is pointing out from */
   int                   j                   /**< Index of arc in graph data structure */
   )
{
   int i;
   int vtx;
   int head;

   head = graph->adjacencylists[j];
   for (i = 0; i < nattacks; ++i)
   {
      vtx = attackpattern[i];
      if ( vtx == head || vtx == source )
         return TRUE;
   }
   return FALSE;
}

/** updates the right hand sides of the KEP problem based on an attack pattern */
static
SCIP_RETCODE SCIPupdateObjCoefsKEPProblem(
   SCIP*                 kepscip,            /**< SCIP instance of the KEP */
   SCIP_PROBDATA*        kepdata,            /**< problem data */
   Cycles*               cycles,             /**< cycles in graph */
   PositionedArcs*       posarcs,            /**< positioned arcs in graph */
   int*                  attackpattern,      /**< array encoding attack pattern */
   int                   nattacks            /**< number of attackes in attackpattern */
   )
{
   int ncycles;
   int nposarcs;
   int nnodes;
   int npairs;
   int c;
   int i;
   int j;
   int v;
   int index;
   SCIP_VAR** cyclevars;
   SCIP_VAR** arcvars;
   SCIP_VAR** dummyarcvars;
   Graph* graph;
   SCIP_Bool attacked;

   assert ( kepscip != NULL );
   assert ( kepdata != NULL );
   assert ( attackpattern != NULL );

   cyclevars = SCIPKEPdataPICEFGetCyclevars(kepdata);
   ncycles = SCIPKEPdataPICEFGetNCycles(kepdata);

   arcvars = SCIPKEPdataPICEFGetArcvars(kepdata);
   nposarcs = SCIPKEPdataPICEFGetNPosarcs(kepdata);

   nnodes = SCIPKEPdataPICEFGetNumNodes(kepdata);

   dummyarcvars = SCIPKEPdataPICEFGetDummyArcvars(kepdata);
   graph = SCIPKEPdataPICEFGetGraph(kepdata);
   npairs = graph->npairs;

   for (c = 0; c < ncycles; ++c)
   {
      attacked = SCIPcycleIsAttacked(attackpattern, nattacks, cycles, c);
      if ( attacked )
      {
         SCIP_CALL( SCIPchgVarObj( kepscip, cyclevars[c], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj( kepscip, cyclevars[c], cycles->cycleweights[c]*nnodes+1.0) );
      }
   }

   for (c = 0; c < nposarcs; ++c)
   {
      SCIP_CALL( SCIPchgVarObj( kepscip, arcvars[c], 0.5 / nnodes) );
   }

   index = 0;
   for (i = 0; i < nnodes; ++i)
   {
      for (j = graph->adjacencylistbegins[i]; j < graph->adjacencylistbegins[i+1]; ++j)
      {
         for (v = 0; v < nposarcs; ++v)
         {
            if ( posarcs->nodelists[2*v] == i && posarcs->nodelists[2*v+1] == graph->adjacencylists[j])
            {
               if (i < npairs)
               {
                  SCIP_CALL( SCIPchgVarObj(kepscip, dummyarcvars[index], posarcs->arcweights[v]*nnodes) );
               }
               else
               {
                  SCIP_CALL( SCIPchgVarObj(kepscip, dummyarcvars[index], posarcs->arcweights[v]*nnodes+1.0) );
               }
               break;
            }
         }
         index++;
      }
   }

   return SCIP_OKAY;
}

/** method used to iteratively solve the Benders model and the corresponding KEPs */
SCIP_RETCODE SCIPsolveBendersPICEFModel(
   SCIP*                 bendersscip,        /**< SCIP data structure for Benders model */
   SCIP*                 kepscip,            /**< SCIP data structure for the subproblem of Benders (a KEP) */
   SCIP*                 masterscip,         /**< SCIP data structure for the master model */
   SCIP_Real             masterobj,          /**< Objective value of the master problem */
   Cycles*               cycles,             /**< Data structure needed for cycle weights when adding Benders cuts */
   PositionedArcs*       posarcs,            /**< Data structure needed for arc weights when adding Benders cuts */
   SCIP_Real             timelimit,          /**< time limit for solving the Benders model */
   SCIP_Bool*            didnotfinish,       /**< pointer to store whether we have hit the time limit */
   SCIP_Bool*            optimal,            /**< pointer to store whether the method terminated optimally */
   SCIP_Real*            timestage2,         /**< pointer to store time spent in stage 2 */
   SCIP_Real*            timestage3          /**< pointer to store time spent in stage 3 */
   )
{
   int* attackpattern;
   int c;
   int cnt = 0;
   int nattacks;
   int nnodes;
   int adversarybound;
   SCIP_VAR** vars;
   SCIP_VAR** varboundvars;
   SCIP_Real* vals;
   SCIP_Real* varboundvals;
   SCIP_Bool liftsols;
   int maxnvarsinsolcons;
   int ncycles;
   int nposarcs;
   int narcs;
   int old_narcvars;
   int new_narcvars;
   int* arcindices;
   SCIP_Real bendersobj;
   SCIP_Real kepobj;
   SCIP_SOL* sol;
   SCIP_SOL* kepsol;
   SCIP_VAR** cyclevars;
   SCIP_VAR** dummyarcvars;
   SCIP_VAR** arcvars;
   SCIP_VAR** benderscyclevars;
   SCIP_VAR** bendersarcvars;
   SCIP_PROBDATA* bendersdata;
   SCIP_PROBDATA* kepdata;
   SCIP_Real memlimit;
   SCIP_Real begintime;
   SCIP_Real endtime;
   SCIP_Real smallestub;
   SCIP_Real totaltime = 0.0;

   assert( bendersscip != NULL );
   assert( kepscip != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
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

   adversarybound = SCIPbendersdataPICEFGetAdversaryBound(bendersdata);
   ncycles = SCIPbendersdataPICEFGetNCycles(bendersdata);
   nposarcs = SCIPKEPdataPICEFGetNPosarcs(kepdata);
   maxnvarsinsolcons = SCIPKEPdataPICEFGetNumNodes(kepdata);
   nnodes = SCIPKEPdataPICEFGetNumNodes(kepdata);
   narcs = SCIPKEPdataPICEFGetNumArcs(kepdata);

   SCIP_CALL( SCIPgetRealParam(masterscip, "limits/memory", &memlimit) );

   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &attackpattern, adversarybound) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &vars, maxnvarsinsolcons) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &vals, maxnvarsinsolcons) );

   /* arrays used to impose constraints on the use of scenario indexed arc variables in the bendersscip */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &varboundvars, 3) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &varboundvals, 3) );

   /* Stronger bound on objective variable, as it is monotone increasing with respect to future iterations */
   SCIP_CALL( SCIPchgVarUb(bendersscip, SCIPbendersdataPICEFGetObjvar(bendersdata), masterobj) );

   /* In this loop we have to make sure that the transformed bendersscip and kepscip are not freed.
    * Otherwise, we cannot access the termination status in solveMasterPICEFProblem
    */
   while ( TRUE )
   {
      int loop;
      SCIP_Real begintimestage3;
      SCIP_Real endtimestage3;

      /* Solve Benders model */
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
            SCIPbendersdataPICEFGetNumNodes(bendersdata), METHOD_BENDERS_PICEF) );
      SCIPinfoMessage(bendersscip, NULL, "Attack pattern: \n");

      for (loop = 0; loop < nattacks; ++loop)
      {
         SCIPinfoMessage(bendersscip, NULL, "%d ", attackpattern[loop]);
      }
      SCIPinfoMessage(bendersscip, NULL, "\n");

      /* Make sure that we can modify the KEP SCIP instance */

      if ( ! liftsols )
      {
         SCIP_CALL( SCIPupdateRhsKEPProblem(kepscip, kepdata, attackpattern, nattacks) );
      }
      else
      {
         SCIP_CALL( SCIPupdateObjCoefsKEPProblem(kepscip, kepdata, cycles, posarcs, attackpattern, nattacks) );
         SCIP_CALL( SCIPupdateUbDummyVars(kepscip, kepdata, attackpattern, nattacks) );
      }

      SCIP_CALL( SCIPsetLimitsAndVerbose(kepscip, begintime, timelimit, MAX(1, memlimit -
            SCIPgetMemUsed(masterscip) / 1048576 - SCIPgetMemUsed(bendersscip) / 1048576), TRUE) );
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
      benderscyclevars = SCIPbendersdataPICEFGetCyclevars(bendersdata);
      cyclevars = SCIPKEPdataPICEFGetCyclevars(kepdata);
      arcvars = SCIPKEPdataPICEFGetArcvars(kepdata);

      kepobj = 0;
      cnt = 0;

      /* Compute the objective of the KEP problem and based on solution, add solution constraint to bendersscip */
      if ( liftsols )
      {
         dummyarcvars = SCIPKEPdataPICEFGetDummyArcvars(kepdata);
         for (c = 0; c < ncycles; ++c)
         {
            /* Only take into account cycles in kepsol that are not attacked by the current attack pattern */
            if ( SCIPgetSolVal(kepscip, kepsol, cyclevars[c]) > 0.5 )
            {
               if ( SCIPvarGetObj(cyclevars[c]) > nnodes )
                  kepobj += cycles->cycleweights[c];
               vars[cnt] = benderscyclevars[c];
               vals[cnt++] = cycles->cycleweights[c];
            }
         }
         for (c = 0; c < narcs; ++c)
         {
            if ( SCIPgetSolVal(kepscip, kepsol, dummyarcvars[c]) > 0.5 )
            {
               /* When the objective coefficient is larger than nnodes, then we know
                * that the arc has weight 1, otherwise it is 0
                */
               if ( SCIPvarGetObj(dummyarcvars[c]) > 1.0)
                  kepobj += 1.0;
            }
         }
      }
      else
      {
         for (c = 0; c < ncycles; ++c)
         {
            /* Only take into account cycles in kepsol that are not attacked */
            if ( SCIPgetSolVal(kepscip, kepsol, cyclevars[c]) > 0.5 )
            {
               kepobj += cycles->cycleweights[c];
               vars[cnt] = benderscyclevars[c];
               vals[cnt++] = cycles->cycleweights[c];
            }
         }

         for (c = 0; c < nposarcs; ++c)
         {
            /* Only take into account posarc in kepsol that are not attacked */
            if ( SCIPgetSolVal(kepscip, kepsol, arcvars[c]) > 0.5 )
               kepobj += posarcs->arcweights[c];
         }
      }

      /* Free the Benders model before we can add new cuts */
      SCIP_CALL( SCIPfreeTransform(bendersscip) );

      old_narcvars = SCIPbendersdataPICEFGetNArcvars(bendersdata);
      SCIPbendersdataPICEFAddArcVars(bendersscip, kepscip, kepsol, varboundvals, varboundvars, FALSE);
      new_narcvars = SCIPbendersdataPICEFGetNArcvars(bendersdata);
      arcindices = SCIPbendersdataPICEFGetArcIndices(bendersdata);

      bendersarcvars = SCIPbendersdataPICEFGetArcvars(bendersdata);
      for (c = old_narcvars; c < new_narcvars; ++c)
      {
         vars[cnt] = bendersarcvars[c];
         vals[cnt++] = posarcs->arcweights[arcindices[c]];
      }

      SCIPinfoMessage(bendersscip, NULL, "@95 KEP objective: %f\n", kepobj);

      if ( SCIPisLT(kepscip, kepobj, masterobj - 0.5) )
      {
         SCIPinfoMessage(kepscip, NULL, "The recourse value %f of the currently considered attack pattern violates the master solution.\n", kepobj);
         *optimal = TRUE;
         break;
      }

      /* the vars are being enabled through resetting the lower bound */
      for (c = 0; c < cnt; ++c)
      {
         SCIP_CALL( SCIPchgVarLb(bendersscip, vars[c], 0.0) );
      }

      /* Also add dummy variable to vars of the cut */
      vals[cnt] = -1.0;
      vars[cnt++] = SCIPbendersdataPICEFGetObjvar(bendersdata);

      SCIP_CALL( SCIPbendersdataPICEFAddSolCons(bendersscip, vals, vars, cnt) );

      /* Stronger bound on objective variable, as it is monotone increasing with respect to future iterations */
      SCIP_CALL( SCIPchgVarLb(bendersscip, SCIPbendersdataPICEFGetObjvar(bendersdata), bendersobj) );
      if ( kepobj < smallestub )
      {
         SCIP_CALL( SCIPchgVarUb(bendersscip, SCIPbendersdataPICEFGetObjvar(bendersdata), kepobj) );
         smallestub = kepobj;
      }
      SCIP_CALL( SCIPfreeTransform(kepscip) );
      endtime = clock();
      totaltime += (SCIP_Real) (endtime - begintime) / CLOCKS_PER_SEC;

      endtimestage3 = clock();
      *timestage3 += (SCIP_Real) (endtimestage3 - begintimestage3) / CLOCKS_PER_SEC;
   }

   SCIPfreeBlockMemoryArrayNull(bendersscip, &attackpattern, adversarybound);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &vals, maxnvarsinsolcons);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &vars, maxnvarsinsolcons);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &varboundvals, 3);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &varboundvars, 3);

   endtime = clock();
   *timestage2 = (SCIP_Real) (endtime - begintime) / CLOCKS_PER_SEC - *timestage3;

   return SCIP_OKAY;
}

/** add cut based on initial master solution to bendersscip */
SCIP_RETCODE SCIPaddInitialPICEFBendersCut(
   SCIP*                masterscip,          /**< SCIP instance of master problem */
   SCIP*                bendersscip,         /**< SCIP instance of benders problem */
   SCIP_SOL*            sol                  /**< solution of master SCIP */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_PROBDATA* bendersdata;
   SCIP_VAR** benderscyclevars;
   SCIP_VAR** bendersarcvars;
   SCIP_VAR** masterinitcycles;
   SCIP_VAR** vars;
   SCIP_VAR** varboundvars;
   SCIP_Real* vals;
   SCIP_Real* varboundvals;
   Cycles* cycles;
   PositionedArcs* posarcs;
   int* arcindices;
   int nnodes;
   int ncycles;
   int new_narcvars, old_narcvars;
   int cnt;
   int c;

   assert( masterscip != NULL );
   assert( bendersscip != NULL );

   probdata = SCIPgetProbData(masterscip);
   assert( probdata != NULL );

   bendersdata = SCIPgetProbData(bendersscip);
   assert( bendersdata != NULL );

   benderscyclevars = SCIPbendersdataPICEFGetCyclevars(bendersdata);
   masterinitcycles = masterPICEFProblemGetXCyclevarinit(probdata);

   cycles = masterPICEFProblemGetCycles(probdata);
   posarcs = masterPICEFProblemGetPosarcs(probdata);

   nnodes = SCIPbendersdataPICEFGetNumNodes(bendersdata);
   ncycles = SCIPbendersdataPICEFGetNCycles(bendersdata);

   SCIP_CALL( SCIPallocBufferArray(masterscip, &vars, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &vals, nnodes + 1) );

   SCIP_CALL( SCIPallocBufferArray(masterscip, &varboundvals, 3) );
   SCIP_CALL( SCIPallocBufferArray(masterscip, &varboundvars, 3) );

   cnt = 0;

   for (c = 0; c < ncycles; ++c)
   {
      /* Only take into account cycles in kepsol that are not attacked by the current attack pattern */
      if ( SCIPgetSolVal(masterscip, sol, masterinitcycles[c]) > 0.5 )
      {
         vars[cnt] = benderscyclevars[c];
         vals[cnt++] = cycles->cycleweights[c];
      }
   }

   old_narcvars = SCIPbendersdataPICEFGetNArcvars(bendersdata);
   SCIPbendersdataPICEFAddArcVars(bendersscip, masterscip, sol, varboundvals, varboundvars, TRUE);
   new_narcvars = SCIPbendersdataPICEFGetNArcvars(bendersdata);
   arcindices = SCIPbendersdataPICEFGetArcIndices(bendersdata);

   bendersarcvars = SCIPbendersdataPICEFGetArcvars(bendersdata);
   for (c = old_narcvars; c < new_narcvars; ++c)
   {
      vars[cnt] = bendersarcvars[c];
      vals[cnt++] = posarcs->arcweights[arcindices[c]];
   }

   /* also add dummy variable to vars of the cut */
   vals[cnt] = -1.0;
   vars[cnt++] = SCIPbendersdataPICEFGetObjvar(bendersdata);

   /* the vars are being enabled through resetting the lower bound */
   for (c = 0; c < cnt; ++c)
   {
      SCIP_CALL( SCIPchgVarLb(bendersscip, vars[c], 0.0) );
   }

   /* add cut based on initial master solution to bendersscip */
   SCIP_CALL( SCIPbendersdataPICEFAddSolCons(bendersscip, vals, vars, cnt) );

   /* Free buffer memory arrays */
   SCIPfreeBufferArray(masterscip, &vals);
   SCIPfreeBufferArray(masterscip, &vars);
   SCIPfreeBufferArray(masterscip, &varboundvals);
   SCIPfreeBufferArray(masterscip, &varboundvars);

   return SCIP_OKAY;
}

/** solves the master problem */
SCIP_RETCODE solveMasterPICEFProblem(
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to arc structures of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Real             timelimit,          /**< time limit to solve the problem */
   SCIP_Real             initialtime,        /**< time spent in initialization */
   const char*           settings,           /**< Possible name of setting file */
   SCIP_Bool             verbose             /**< whether we print SCIP's logs */
   )
{
   /* General part of master problem solve method */
   SCIP_Real masterobj;
   SCIP_Real subobj;
   SCIP_SOL* sol;
   SCIP_PROBDATA* probdata;
   int nattacks;
   int* attackpattern;
   int* scenarios;
   int nscenarios;
   int i;
   int j;
   int pos;
   clock_t begintime;
   clock_t endtime;
   clock_t curtime;
   int method;
   SCIP_Real memlimit;
   SCIP_Real modtimelimit;
   SCIP_Real newtimelimit;
   SCIP_Real totaltime;
   SCIP_Real timesecondstage = 0.0;
   SCIP_Real timethirdstage = 0.0;
   int cnt = 0;
   int itercnt = 0;

   SCIP_Bool didnotfinish = FALSE;

   assert( masterscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
   assert( adversarybound >= 0 );

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
      SCIP_Real beginsubtime;
      SCIP_Real endsubtime;
      SCIP_Real timestage2;
      SCIP_Real timestage3;

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      SCIP_CALL( SCIPsetLimitsAndVerbose(masterscip, begintime, modtimelimit, memlimit, verbose) );
      SCIP_CALL( SCIPsolve(masterscip) );

      if ( solvingloopShallTerminate(masterscip, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      /* Obtain best solution and best objective */
      sol = SCIPgetBestSol(masterscip);
      masterobj = SCIPgetSolOrigObj(masterscip, sol);

      SCIPinfoMessage(masterscip, NULL, "@98 masteriteration %d: objective value %f.\n", ++itercnt, masterobj);

      /* adapt the weights of cycles and chains based on initial solution of master problem */
      SCIP_CALL( adaptCyclePosarcWeights(masterscip, sol, graph, cycles, posarcs) );

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREEMASTERPROBLEM;

      beginsubtime = clock();

      SCIP_CALL( SCIPcreateBasicInstance(masterscip, &scip, verbose) );

      /* create the SCIP instance based on the PICEF benders model */
      SCIP_CALL( SCIPcreateBendersPICEFModel(scip, graph, cycles, posarcs, adversarybound) );
      SCIP_CALL( SCIPsetObjIntegral(scip) );

      endsubtime = clock();
      timesecondstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;

      beginsubtime = clock();

      /* create Benders subproblem model */
      SCIP_CALL( SCIPcreateBasicInstance(masterscip, &subscip, verbose) );

      SCIP_CALL( SCIPcreateBendersPICEFSubModel(subscip, graph, cycles, posarcs) );
      SCIP_CALL( SCIPsetObjIntegral(subscip) );

      endsubtime = clock();
      timethirdstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;

      SCIPinfoMessage(masterscip, NULL, "@97 Set up Benders model %d.\n", cnt);

      /* add cut based on initial master solution to bendersscip */
      SCIP_CALL( SCIPaddInitialPICEFBendersCut(masterscip, scip, sol) );

      if ( solvingloopShallTerminate(NULL, begintime, modtimelimit) )
         goto FREESUBPROBLEM;

      curtime = clock();
      newtimelimit = modtimelimit - (SCIP_Real) (curtime - begintime) / CLOCKS_PER_SEC;

      SCIP_CALL( SCIPsolveBendersPICEFModel(scip, subscip, masterscip, masterobj, cycles, posarcs, newtimelimit,
            &didnotfinish, &optimal, &timestage2, &timestage3) );
      timesecondstage += timestage2;
      timethirdstage += timestage3;

      if ( optimal )
      {
         subobj = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
         SCIPinfoMessage(masterscip, NULL, "@96 Benders problem iteration %d: objective value %f.\n", cnt, subobj);
      }

      /* if we have found an optimal solution */
      if ( optimal && SCIPisLE(masterscip, masterobj, subobj) )
      {
         SCIP_CALL( SCIPprintSol(masterscip, SCIPgetBestSol(masterscip), NULL, FALSE) );

         SCIP_CALL( SCIPfreeTransform(scip) );
         SCIP_CALL( SCIPfree(&scip) );
         SCIP_CALL( SCIPfreeTransform(subscip) );
         SCIP_CALL( SCIPfree(&subscip) );
         SCIP_CALL( SCIPfreeTransform(masterscip) );

         break;
      }
      /* if we hit the time limit or could not solve a subSCIP instance (check for optimal, because otherwise,
       * we transformed version of scip might have already been freed)
       */
      else if ( didnotfinish || ( ! optimal &&
            (solvingloopShallTerminate(scip, begintime, modtimelimit) || solvingloopShallTerminate(subscip, begintime, modtimelimit))) )
         goto FREESUBPROBLEM;

      beginsubtime = clock();

      /* extract attack pattern and update master problem */
      SCIP_CALL( SCIPfreeTransform(masterscip) );

      SCIP_CALL( SCIPgetAttackPattern(scip, SCIPgetBestSol(scip), attackpattern, &nattacks, graph->nnodes, method) );
      SCIPinfoMessage(masterscip, NULL, "nattacks: %d\n", nattacks);
      SCIPinfoMessage(masterscip, NULL, "first attack: %d\n", attackpattern[0]);
      SCIP_CALL( SCIPupdateMasterPICEFProblem(masterscip, attackpattern, nattacks) );

      endsubtime = clock();
      timesecondstage += (SCIP_Real) (endsubtime - beginsubtime) / CLOCKS_PER_SEC;

      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIP_CALL( SCIPfree(&scip) );
      SCIP_CALL( SCIPfreeTransform(subscip) );
      SCIP_CALL( SCIPfree(&subscip) );

      /* gotos if we could not solve a problem or hit the time limit */
      if ( FALSE )
      {
      FREESUBPROBLEM:
         SCIP_CALL( SCIPfreeTransform(subscip) );
         SCIP_CALL( SCIPfree(&subscip) );
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

   scenarios = masterPICEFProblemGetScenarios(probdata);
   nscenarios = masterPICEFProblemGetNScenarios(probdata);

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
