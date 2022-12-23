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

/**@file   probdata_benders_picef.c
 * @brief  Problem data for Benders type approach with PICEF formulation
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
 * A list of all interface methods can be found in probdata_benders_picef.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_benders_picef.h"
#include "probdata_benders_subproblem_picef.h"
#include "probdata_master_kidney_picef.h"
#include "scip/scip.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"

#include "typedefs.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the kidney exchange problem, all variables which are created, and all
 * constraints.
 */

/** Problem data which is accessible in all places */
typedef struct SCIP_ProbData_KeepUnaffectedCC
{
   SCIP*                 masterscip;            /**< SCIP data structure of master problem */
   SCIP_SOL*             mastersol;             /**< SCIP solution of master problem */
   int*                  initcycles;            /**< array of cycle indices selected in master solution */
   int*                  initposarcs;           /**< array of posarc indices selected in master solution */
   int                   ninitcycles;           /**< number of initial cycles */
   int                   ninitposarcs;          /**< number of initial position indexed arcs */

   int                   ncycleubconss;         /**< number of initial cycle upper bound constraints */
   SCIP_VAR**            yvarcovered;           /**< array of vertex variables indicating whether the vertex is covered by a nonattacked cycle / (sub)chain */
   SCIP_VAR**            dummyinitposarcvars;   /**< array of dummy posarcs for the initially selected posarcs indicating whether or not it is affected by an attack */

   SCIP_CONS**           initcycleubconss;      /**< array of conss for upper bounding the variables associated with initial cycles */
   SCIP_CONS**           dummyinitposarclbconss;/**< array of conss for lower bounding the initial posarc dummy variables */
   SCIP_CONS**           dummyinitposarcubconss;/**< array of conss for upper bounding the initial posarc dummy variables */
   SCIP_CONS**           nodecoverconss;        /**< array of conss for linking yvarcovered to dummyinitcyclevars / dummyinitposarcvars */
} SCIP_PROBDATA_KUCC;

struct SCIP_ProbData
{
   Graph*                graph;              /**< pointer to underlying graph */
   int                   nnodes;             /**< number of nodes in graph */
   int                   adversarybound;     /**< bound on adversary attack */

   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   PositionedArcs*       posarcs;            /**< pointer to position indexed arc structure of underlying graph */

   SCIP_VAR*             objvar;             /**< artificial objective variable modeling integer objective values */
   SCIP_VAR**            cyclevars;          /**< array of cyclevars: per solution, we keep the cycle
                                              *   if none of its vertices is attacked */
   SCIP_VAR**            arcvars;            /**< array of arcvars: per solution, we keep the arc if none of its vertices are attacked */
   SCIP_VAR**            uvars;              /**< array of u-vars: adversary variables for attacking nodes */

   int                   nsolconss;          /**< number of generated KEP solutions forming a constraint */
   int                   narcvars;           /**< number of generated solution arc variables */
   int                   maxnsolconss;       /**< maximum number of generated KEP solutions that fits into solconss array */
   int                   maxnarcvars;        /**< maximum number of arc vars based on generated KEP solutions */
   int*                  arcindices;         /**< array that keeps track of arc indices (1 to nposarcs) occurring
                                              *   in scenario indexed arcvars */
   SCIP_CONS*            attackboundcons;    /**< constraint bounding adversary attacks */
   SCIP_CONS**           solconss;           /**< array of linear inequalities, each corresponding
                                              *   to a kidney exchange solution */
   SCIP_CONS**           cyclelbconss;         /**< cycle constraints (cyclevars[c] = 1 if and only if
                                              *   all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           posarcconss;        /**< arc constraints (arcvar[a,k] = 1 if and only if all uvars[v] = 1
                                              *   for v in arc a and there is a preceding arc used)*/

   int                   policytype;         /**< parameter indicating which recourse policy type is considered */
   SCIP_PROBDATA_KUCC*   probdatakucc;       /**< problem data associated with problem if recourse policy POLICY_KEEPUNAFFECTEDCC is considered */
};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE bendersdataPICEFCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_PROBDATA**       bendersdata,        /**< pointer to Benders data */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   Cycles*               cycles,             /**< cycle structure of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure of underlying graph */

   SCIP_VAR**            cyclevars,          /**< array of cyclevars: per solution, we keep the cycle if none of
                                              *   its vertices is attacked */
   SCIP_VAR**            arcvars,            /**< array of arcvars: per solution, we keep the arc if none of
                                              *   its vertices are attacked and preceding arcs are used in the solution */
   SCIP_VAR**            uvars,              /**< array of u-vars: adversary variables for attacking nodes */
   SCIP_VAR*             objvar,             /**< artificial objective variable modeling integer objective values */

   int                   nsolconss,          /**< number of generated KEP solutions forming a constraint */
   int                   narcvars,           /**< number of generated solution arc variables */
   int                   maxnsolconss,       /**< maximum number of generated KEP solutions that fits into solconss array */
   int                   maxnarcvars,        /**< maximum number of arc vars based on generated KEP solutions */
   int*                  arcindices,         /**< array that keeps track of nodes occurring in scenario indexed arcvars */
   int                   adversarybound,     /**< bound on adversary attack */

   SCIP_CONS**           solconss,           /**< array of linear inequalities, each corresponding to a kidney exchange solution */
   SCIP_CONS*            attackboundcons,    /**< constraint bounding adversary attacks */
   SCIP_CONS**           cyclelbconss,       /**< cycle constraints (cyclevars[c] = 1 if and only if
                                              *   all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           posarcconss,        /**< arc varbound constraints (arcvar[a,k] = 1 if and only if all uvars[v] = 1
                                              *   for v in arc a and there is a preceding arc used)*/

   int                   policytype,         /**< parameter indicating which recourse policy type is considered */

   /* Following input parameters correspond to recourse policy POLICY_KEEPUNAFFECTEDCC */

   SCIP*                 masterscip,         /**< SCIP data structure of master problem */
   SCIP_SOL*             mastersol,          /**< SCIP solution of master problem */

   SCIP_VAR**            yvarcovered,           /**< array of vertex variables indicating whether the vertex is covered by a nonattacked cycle / (sub)chain */
   SCIP_VAR**            dummyinitposarcvars,   /**< array of dummy posarcs for the initially selected posarcs indicating whether or not it is affected by an attack */

   SCIP_CONS**           initcycleubconss,      /**< array of conss for upper bounding the variables associated with initial cycles */
   SCIP_CONS**           dummyinitposarclbconss,/**< array of conss for lower bounding the initial posarc dummy variables */
   SCIP_CONS**           dummyinitposarcubconss,/**< array of conss for upper bounding the initial posarc dummy variables */
   SCIP_CONS**           nodecoverconss        /**< array of conss for linking yvarcovered to dummyinitcyclevars / dummyinitposarcvars */

   )
{
   int ncycles;
   int ninitcycles;
   int ninitposarcs;
   int ncycleubconss;
   int i;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
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
   (*bendersdata)->posarcs = posarcs;
   (*bendersdata)->policytype = policytype;

   /* possible copy variable arrays */
   if ( uvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->uvars, uvars, nnodes) );
   }
   else
      (*bendersdata)->uvars = NULL;

   (*bendersdata)->objvar = objvar;

   if ( cyclevars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cyclevars, cyclevars, ncycles) );
   }
   else
      (*bendersdata)->cyclevars = NULL;

   if ( arcvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->arcvars, arcvars, maxnarcvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->arcindices, arcindices, maxnarcvars) );
      (*bendersdata)->narcvars = narcvars;
      (*bendersdata)->maxnarcvars = maxnarcvars;
   }
   else
   {
      (*bendersdata)->arcvars = NULL;
      (*bendersdata)->arcindices = NULL;
      (*bendersdata)->narcvars = 0;
      (*bendersdata)->maxnarcvars = 0;
   }

   /* possible copy constraint arrays */
   if ( nsolconss > 0 )
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

   if ( attackboundcons != NULL )
      (*bendersdata)->attackboundcons = attackboundcons;
   else
      (*bendersdata)->attackboundcons = NULL;

   if ( cyclelbconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cyclelbconss, cyclelbconss, ncycles) );
   }
   else
      (*bendersdata)->cyclelbconss = NULL;

   if ( posarcconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->posarcconss, posarcconss, maxnarcvars) );
   }
   else
      (*bendersdata)->posarcconss = NULL;


   /* The following allocations need only be run whenever the policy type is KeepUnaffectedCyclesChains,
    * since we then need additional constraints to model the second stage of the trilevel problem
    */

   if( policytype == POLICY_KEEPUNAFFECTEDCC )
   {

      SCIP_VAR** mastercyclevars;
      SCIP_VAR** masterarcvars;


      mastercyclevars = masterPICEFProblemGetXCyclevarinit(SCIPgetProbData(masterscip));
      masterarcvars = masterPICEFProblemGetArcvarinit(SCIPgetProbData(masterscip));

      SCIP_CALL( SCIPallocBlockMemory(bendersscip, &(*bendersdata)->probdatakucc) );

      SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->initcycles, nnodes) );
      SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->initposarcs, nnodes) );

      ninitcycles = 0;
      ninitposarcs = 0;
      ncycleubconss = 0;

      for( i = 0; i < cycles->ncycles; ++i )
      {
         if( SCIPgetSolVal(masterscip, mastersol, mastercyclevars[i]) > 0.5 )
         {
            (*bendersdata)->probdatakucc->initcycles[ninitcycles++] = i;
            ncycleubconss += cycles->nodelistsbegin[i+1] - cycles->nodelistsbegin[i];
         }
      }

      for( i = 0; i < posarcs->nposarcs; ++i )
      {
         if( SCIPgetSolVal(masterscip, mastersol, masterarcvars[i]) > 0.5 )
            (*bendersdata)->probdatakucc->initposarcs[ninitposarcs++] = i;
      }

      (*bendersdata)->probdatakucc->masterscip = masterscip;
      (*bendersdata)->probdatakucc->mastersol = mastersol;
      (*bendersdata)->probdatakucc->ninitcycles = ninitcycles;
      (*bendersdata)->probdatakucc->ninitposarcs = ninitposarcs;
      (*bendersdata)->probdatakucc->ncycleubconss = ncycleubconss;

      if ( yvarcovered != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->yvarcovered, yvarcovered, nnodes) );
      }
      else
         (*bendersdata)->probdatakucc->yvarcovered = NULL;

      if ( dummyinitposarcvars != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcvars, dummyinitposarcvars, ninitposarcs) );
      }
      else
         (*bendersdata)->probdatakucc->dummyinitposarcvars = NULL;

      if ( initcycleubconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->initcycleubconss, initcycleubconss, ncycleubconss) );
      }
      else
         (*bendersdata)->probdatakucc->initcycleubconss = NULL;

      if ( dummyinitposarclbconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarclbconss, dummyinitposarclbconss, ninitposarcs) );
      }
      else
         (*bendersdata)->probdatakucc->dummyinitposarclbconss = NULL;

      if ( dummyinitposarcubconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcubconss, dummyinitposarcubconss, 2*ninitposarcs) );
      }
      else
         (*bendersdata)->probdatakucc->dummyinitposarcubconss = NULL;

      if ( nodecoverconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->nodecoverconss, nodecoverconss, nnodes) );
      }
      else
         (*bendersdata)->probdatakucc->nodecoverconss = NULL;

   }

   return SCIP_OKAY;
}


/** frees the memory of the given problem data */
static
SCIP_RETCODE bendersdataPICEFFree(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_PROBDATA**       bendersdata         /**< pointer to Benders data */
   )
{
   int i;
   int nnodes;
   int ncycles;
   int maxnarcvars;
   int narcvars;
   int nsolconss;
   int policytype;


   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( (*bendersdata)->graph != NULL );
   assert( (*bendersdata)->cycles != NULL );

   nnodes = (*bendersdata)->nnodes;
   ncycles = (*bendersdata)->cycles->ncycles;
   narcvars = (*bendersdata)->narcvars;
   maxnarcvars = (*bendersdata)->maxnarcvars;
   nsolconss = (*bendersdata)->nsolconss;
   policytype = (*bendersdata)->policytype;

   assert( nnodes > 0 );
   assert( ncycles >= 0 );
   assert( maxnarcvars >= 0 );

   /* release all variables */
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->uvars[i]) );
   }
   SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->objvar) );

   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->cyclevars[i]) );
   }

   for (i = 0; i < narcvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->arcvars[i]) );
   }

   /** release constraints and free constraint memory arrays */
   for (i = 0; i < nsolconss; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->solconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->solconss, (*bendersdata)->maxnsolconss);

   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->cyclelbconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->cyclelbconss, ncycles);

   for (i = 0; i < narcvars; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->posarcconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->posarcconss, maxnarcvars);

   SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->attackboundcons) );

   if( policytype == POLICY_KEEPUNAFFECTEDCC )
   {
      int ninitposarcs;
      int ncycleubconss;

      ninitposarcs = (*bendersdata)->probdatakucc->ninitposarcs;
      ncycleubconss = (*bendersdata)->probdatakucc->ncycleubconss;

      for (i = 0; i < nnodes; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->probdatakucc->yvarcovered[i]) );
      }

      for (i = 0; i < ninitposarcs; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcvars[i]) );
      }

      /** release constraints and free constraint memory arrays */
      for (i = 0; i < ncycleubconss; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->probdatakucc->initcycleubconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->initcycleubconss, ncycleubconss);

      /** release constraints and free constraint memory arrays */
      for (i = 0; i < ninitposarcs; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarclbconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarclbconss, ninitposarcs);

      for (i = 0; i < 2*ninitposarcs; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcubconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcubconss, 2*ninitposarcs);

      for (i = 0; i < nnodes; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->probdatakucc->nodecoverconss[i]) );
      }

      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->nodecoverconss, nnodes);

      SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->probdatakucc->yvarcovered, nnodes);
      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->dummyinitposarcvars, ninitposarcs);

      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->initcycles, nnodes);
      SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->probdatakucc->initposarcs, nnodes);

      SCIPfreeBlockMemory(bendersscip, &(*bendersdata)->probdatakucc);
   }

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->uvars, nnodes);
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->arcvars, maxnarcvars);
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->arcindices, maxnarcvars);

   /* free bendersdata */
   SCIPfreeBlockMemory(bendersscip, bendersdata);

   return SCIP_OKAY;
}


/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialBendersVars(
   SCIP*                 bendersscip,        /**< SCIP instance */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];

   int nnodes;
   int ncycles;
   int policytype;
   int i;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );
   assert( bendersdata->graph->nodelist != NULL );

   policytype = bendersdata->policytype;

   nnodes = bendersdata->graph->nnodes;
   ncycles = bendersdata->cycles->ncycles;

   /* create x_c-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->cyclevars, ncycles) );

   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->cyclevars[i], name, 1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->cyclevars[i]) );
   }

   /* create u-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->uvars, nnodes) );
   for (i = 0; i < nnodes; ++i)
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

   if( policytype == POLICY_KEEPUNAFFECTEDCC )
   {
      int ninitposarcs = bendersdata->probdatakucc->ninitposarcs;

      /* create node variables indicating if it is covered by a nonattacked initial cycle / (sub)chain */
      SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->probdatakucc->yvarcovered, nnodes) );
      for (i = 0; i < nnodes; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yvarcovered_%d", i);

         SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->probdatakucc->yvarcovered[i], name, 0.0, 1.0, 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->probdatakucc->yvarcovered[i]) );
      }

      /* create dummy initial posarc variables */
      SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &bendersdata->probdatakucc->dummyinitposarcvars, ninitposarcs) );
      for (i = 0; i < ninitposarcs; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyinitposarcvar_%d", i);

         SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->probdatakucc->dummyinitposarcvars[i], name, 0.0, 1.0, 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->probdatakucc->dummyinitposarcvars[i]) );
      }

   }

   return SCIP_OKAY;
}


/** creates the initial constraints of the problem */
static
SCIP_RETCODE SCIPcreateInitialBendersConstraints(
   SCIP*                 bendersscip,        /**< SCIP instance */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   Cycles* cycles;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool initial;

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
   assert( bendersdata->cycles->ncycles == cycles->ncycles );

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/cyclechainconssinitial", &initial) );

   /* allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vals, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vars, maxnvarsincons) );

   /* create bound on adversary attack */
   for (i = 0; i < maxnvarsincons; ++i)
   {
      vals[i] = 1.0;
   }

   SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->attackboundcons, "boundattack", nnodes,
         bendersdata->uvars, vals, nnodes - bendersdata->adversarybound, nnodes - bendersdata->adversarybound,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->attackboundcons) );
   /* do not release constraint here, will be done later */

   /* create constraints linking x_c and its corresponding u_v's */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cyclelbconss), ncycles) );
   for (c = 0; c < ncycles; ++c)
   {
      cnt = 0;
      for (i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i)
      {
         vars[cnt] = bendersdata->uvars[cycles->nodelists[i]];
         ++cnt;
      }
      vars[cnt] = bendersdata->cyclevars[c];
      vals[cnt] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclelbcons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cyclelbconss[c], name, cnt + 1,
            vars, vals, -1.0, cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cyclelbconss[c]) );

      /* Reset value of vals[cnt] for next cons */
      vals[cnt] = 1.0;
   }

   SCIPfreeBufferArray(bendersscip, &vars);
   SCIPfreeBufferArray(bendersscip, &vals);

   return SCIP_OKAY;
}

/** creates the initial constraints of the problem when POLICY_KEEPUNAFFECTEDCC is considered */
static
SCIP_RETCODE SCIPcreateInitialBendersConstraintsKUCC(
   SCIP*                 bendersscip,        /**< SCIP instance */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int* initcycles;
   int* initposarcs;
   int ninitcycles;
   int ninitposarcs;
   int nnodes;
   int npairs;
   int ncycles;
   int ncycleubconss;
   Cycles* cycles;
   PositionedArcs* posarcs;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool* hasnodecovercons;
   SCIP_Bool initial;
   SCIP_Bool hasPredecessor;

   int maxnvarsincons;
   int nextinitcycle;
   int i;
   int j;
   int c;
   int v;
   int cnt;
   int idx;
   int considx;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );

   nnodes = bendersdata->graph->nnodes;
   npairs = bendersdata->graph->npairs;
   maxnvarsincons = nnodes + 1;
   initcycles = bendersdata->probdatakucc->initcycles;
   initposarcs = bendersdata->probdatakucc->initposarcs;
   ninitcycles = bendersdata->probdatakucc->ninitcycles;
   ninitposarcs = bendersdata->probdatakucc->ninitposarcs;
   ncycleubconss = bendersdata->probdatakucc->ncycleubconss;

   cycles = bendersdata->cycles;
   ncycles = cycles->ncycles;

   posarcs = bendersdata->posarcs;

   assert( bendersdata->cycles->ncycles == cycles->ncycles );

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/cyclechainconssinitial", &initial) );

   /* allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vals, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vars, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &hasnodecovercons, nnodes) );

   /* create bound on adversary attack */
   for (i = 0; i < maxnvarsincons; ++i)
      vals[i] = 1.0;

   for( i = 0; i < nnodes; ++i )
      hasnodecovercons[i] = FALSE;

   SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->attackboundcons, "boundattack", nnodes,
         bendersdata->uvars, vals, nnodes - bendersdata->adversarybound, nnodes,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->attackboundcons) );
   /* do not release constraint here, will be done later */


   /* create constraints linking yvarcovered variables to dummy initial cycle / posarc variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->probdatakucc->nodecoverconss), nnodes) );
   for( c = 0; c < ninitcycles; ++c )
   {
      cnt = 0;
      idx = initcycles[c];

      vals[cnt] = 1.0;
      vars[cnt++] = bendersdata->cyclevars[idx];

      for( i = cycles->nodelistsbegin[idx]; i < cycles->nodelistsbegin[idx + 1]; ++i )
      {
         v = cycles->nodelists[i];
         vals[cnt] = -1.0;
         vars[cnt] = bendersdata->probdatakucc->yvarcovered[v];
         // SCIP_CALL( SCIPchgVarUb(bendersscip, bendersdata->probdatakucc->yvarcovered[v], 1.0) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecovercons_%d", v);

         /* x_c - y_v = 0 */
         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->nodecoverconss[v], name, cnt + 1,
               vars, vals, 0.0, 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->nodecoverconss[v]) );
         hasnodecovercons[v] = TRUE;
      }
   }

   for( c = 0; c < ninitposarcs; ++c )
   {
      cnt = 0;

      vals[cnt] = 1.0;
      vars[cnt++] = bendersdata->probdatakucc->dummyinitposarcvars[c];

      if( posarcs->nodelists[2*initposarcs[c]] >= npairs )
      {
         /* Such arcs start from an NDD, hence k = 1. We add an extra constraint to say the NDD is covered as well */
         v = posarcs->nodelists[2*initposarcs[c]];

         vals[cnt] = -1.0;
         vars[cnt] = bendersdata->probdatakucc->yvarcovered[v];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecovercons_%d", v);

         /* \eta_ijl - y_i = 0 */
         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->nodecoverconss[v], name, cnt + 1,
               vars, vals, 0.0, 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->nodecoverconss[v]) );

         hasnodecovercons[v] = TRUE;
      }

      v = posarcs->nodelists[2*initposarcs[c] + 1];

      vals[cnt] = -1.0;
      vars[cnt] = bendersdata->probdatakucc->yvarcovered[v];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecovercons_%d", v);

      /* \eta_ijl - y_j = 0 */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->nodecoverconss[v], name, cnt + 1,
            vars, vals, 0.0, 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->nodecoverconss[v]) );

      hasnodecovercons[v] = TRUE;
   }

   for( i = 0; i < nnodes; ++i )
   {
      if( !hasnodecovercons[i] )
      {
         vars[0] = bendersdata->probdatakucc->yvarcovered[i];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecovercons_%d", i);

         /* \eta_ijl - y_v = 0 */
         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->nodecoverconss[i], name, 1,
               vars, vals, 0.0, 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->nodecoverconss[i]) );
      }
   }


   /* Reset value of vals[cnt] for next cons */
   vals[1] = 1.0;

   /* create constraints that provide lower bounds on the initial cycle variables and initial posarc variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->probdatakucc->initcycleubconss), ncycleubconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->probdatakucc->dummyinitposarcubconss), 2*ninitposarcs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->probdatakucc->dummyinitposarclbconss), ninitposarcs) );

   i = 0;
   for( c = 0; c < ninitcycles; ++c )
   {
      cnt = 0;
      idx = initcycles[c];

      vals[cnt] = 1.0;
      vars[cnt++] = bendersdata->cyclevars[idx];
      for( v = cycles->nodelistsbegin[idx]; v < cycles->nodelistsbegin[idx + 1]; ++v )
      {
         vals[cnt] = -1.0;
         vars[cnt] = bendersdata->uvars[cycles->nodelists[v]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "initcycleubconss_%d_%d", idx, i);

         /* x_c <= u_v */
         SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->initcycleubconss[i], name, cnt + 1,
               vars, vals, -SCIPinfinity(bendersscip), 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->initcycleubconss[i++]) );
      }
   }

   i = 0;
   for( c = 0; c < ninitposarcs; ++c )
   {
      cnt = 0;

      vals[cnt] = 1.0;
      vars[cnt++] = bendersdata->probdatakucc->dummyinitposarcvars[c];

      idx = initposarcs[c];
      hasPredecessor = FALSE;
      for( v = 0; v < c; ++v )
      {
         j = initposarcs[v];
         if( posarcs->nodelists[2*j + 1] == posarcs->nodelists[2*idx])
         {
            hasPredecessor = TRUE;
            break;
         }
      }

      /*  upper bound given by target of arc */
      vals[cnt] = -1.0;
      vars[cnt] = bendersdata->uvars[posarcs->nodelists[2*idx + 1]];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyinitposarcubconss_%d", i);

      /* \eta_ijl <= u_j */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->dummyinitposarcubconss[i], name, cnt + 1,
            vars, vals, -SCIPinfinity(bendersscip), 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->dummyinitposarcubconss[i++]) );

      /* Second upper bound depends on whether we have a posarc with k = 1 or k >= 2 */
      vals[cnt] = -1.0;
      if( hasPredecessor )
         /* \eta_ijl <= \eta_ki(l-1) */
         vars[cnt++] = bendersdata->probdatakucc->dummyinitposarcvars[v];
      else
         /* \eta_ijl <= u_i */
         vars[cnt++] = bendersdata->uvars[posarcs->nodelists[2*idx]];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyinitposarcubconss_%d", i);

      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->dummyinitposarcubconss[i], name, cnt,
            vars, vals, -SCIPinfinity(bendersscip), 0.0, initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->dummyinitposarcubconss[i++]) );


      /* Finally, add target of arc for lower bound constraints */
      vals[cnt] = -1.0;
      vars[cnt++] = bendersdata->uvars[posarcs->nodelists[2*idx + 1]];

      considx = i/2 - 1;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyinitposarclbconss_%d", considx);

      /* x_c <= u_v */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->probdatakucc->dummyinitposarclbconss[considx], name, cnt,
            vars, vals, -1.0, SCIPinfinity(bendersscip), initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->probdatakucc->dummyinitposarclbconss[considx]) );

   }

   /* create constraints that provide lower bounds on the cycle variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cyclelbconss), ncycles) );

   idx = 0;
   nextinitcycle = initcycles[idx];
   for (c = 0; c < ncycles; ++c)
   {
      cnt = 0;
      for (i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i)
      {
         vals[cnt] = 1.0;
         vars[cnt++] = bendersdata->uvars[cycles->nodelists[i]];
      }
      vals[cnt] = -1.0;
      vars[cnt++] = bendersdata->cyclevars[c];

      /* if c is not an initial cycle index, then add yvarcovered variables */
      if( c != nextinitcycle )
      {
         for (i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i)
         {
            vals[cnt] = -1.0;
            vars[cnt++] = bendersdata->probdatakucc->yvarcovered[cycles->nodelists[i]];

         }
      }
      else
      {
         if( idx < ninitcycles - 1 )
            nextinitcycle = initcycles[++idx];
         else
            nextinitcycle = -1;
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclelbcons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cyclelbconss[c], name, cnt,
            vars, vals, -SCIPinfinity(bendersscip), cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cyclelbconss[c]) );
   }

   SCIPfreeBufferArray(bendersscip, &vars);
   SCIPfreeBufferArray(bendersscip, &vals);
   SCIPfreeBufferArray(bendersscip, &hasnodecovercons);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of Benders' decomposition problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigBendersPICEF)
{
   SCIPdebugMsg(scip, "free original Benders data\n");

   SCIP_CALL( bendersdataPICEFFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed Benders problem by transforming the original Benders problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransBendersPICEF)
{
   int nnodes;
   int ncycles;
   int ncycleubconss;
   int ninitposarcs;
   int policy;

   SCIP_CALL( SCIPgetIntParam(scip, "kidney/recoursepolicy", &policy) );

   /* create transform bendersdata */
   if( policy == POLICY_KEEPUNAFFECTEDCC )
   {
      SCIP_CALL( bendersdataPICEFCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->posarcs,
         sourcedata->cyclevars, sourcedata->arcvars, sourcedata->uvars, sourcedata->objvar, sourcedata->nsolconss, sourcedata->narcvars,
         sourcedata->maxnsolconss, sourcedata->maxnarcvars, sourcedata->arcindices, sourcedata->adversarybound, sourcedata->solconss,
         sourcedata->attackboundcons, sourcedata->cyclelbconss, sourcedata->posarcconss, sourcedata->policytype, sourcedata->probdatakucc->masterscip, sourcedata->probdatakucc->mastersol,
         sourcedata->probdatakucc->yvarcovered, sourcedata->probdatakucc->dummyinitposarcvars, sourcedata->probdatakucc->initcycleubconss,
         sourcedata->probdatakucc->dummyinitposarclbconss, sourcedata->probdatakucc->dummyinitposarcubconss, sourcedata->probdatakucc->nodecoverconss) );
   }
   else
   {
      SCIP_CALL( bendersdataPICEFCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->posarcs,
         sourcedata->cyclevars, sourcedata->arcvars, sourcedata->uvars, sourcedata->objvar, sourcedata->nsolconss, sourcedata->narcvars,
         sourcedata->maxnsolconss, sourcedata->maxnarcvars, sourcedata->arcindices, sourcedata->adversarybound, sourcedata->solconss,
         sourcedata->attackboundcons, sourcedata->cyclelbconss, sourcedata->posarcconss, sourcedata->policytype, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   nnodes = sourcedata->graph->nnodes;
   ncycles = sourcedata->cycles->ncycles;
   assert( nnodes > 0 );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nsolconss, (*targetdata)->solconss, (*targetdata)->solconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->attackboundcons, &(*targetdata)->attackboundcons) );
   SCIP_CALL( SCIPtransformConss(scip, ncycles, (*targetdata)->cyclelbconss, (*targetdata)->cyclelbconss) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->narcvars, (*targetdata)->posarcconss, (*targetdata)->posarcconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->cyclevars, (*targetdata)->cyclevars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->narcvars, (*targetdata)->arcvars, (*targetdata)->arcvars) );
   SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->uvars, (*targetdata)->uvars) );
   SCIP_CALL( SCIPtransformVar(scip, (*targetdata)->objvar, &(*targetdata)->objvar) );

   if( policy == POLICY_KEEPUNAFFECTEDCC )
   {
      ncycleubconss = sourcedata->probdatakucc->ncycleubconss;
      ninitposarcs = sourcedata->probdatakucc->ninitposarcs;

      SCIP_CALL( SCIPtransformConss(scip, ncycleubconss, (*targetdata)->probdatakucc->initcycleubconss, (*targetdata)->probdatakucc->initcycleubconss) );
      SCIP_CALL( SCIPtransformConss(scip, ninitposarcs, (*targetdata)->probdatakucc->dummyinitposarclbconss, (*targetdata)->probdatakucc->dummyinitposarclbconss) );
      SCIP_CALL( SCIPtransformConss(scip, 2*ninitposarcs, (*targetdata)->probdatakucc->dummyinitposarcubconss, (*targetdata)->probdatakucc->dummyinitposarcubconss) );
      SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->probdatakucc->nodecoverconss, (*targetdata)->probdatakucc->nodecoverconss) );

      SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->probdatakucc->yvarcovered, (*targetdata)->probdatakucc->yvarcovered) );
      SCIP_CALL( SCIPtransformVars(scip, ninitposarcs, (*targetdata)->probdatakucc->dummyinitposarcvars, (*targetdata)->probdatakucc->dummyinitposarcvars) );
   }

   return SCIP_OKAY;
}

/* frees user data of transformed Benders problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransBendersPICEF)
{
   SCIPdebugMsg(scip, "free transformed Benders' decomposition data\n");

   SCIP_CALL( bendersdataPICEFFree(scip, probdata) );

   return SCIP_OKAY;
}
/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPbendersdataPICEFCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure of master problem */
   SCIP_SOL*             mastersol,          /**< SCIP solution structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* bendersdata;
   int policy;

   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );

   assert( bendersscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(bendersscip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(bendersscip, probdelorigBendersPICEF) );
   SCIP_CALL( SCIPsetProbTrans(bendersscip, probtransBendersPICEF) );
   SCIP_CALL( SCIPsetProbDeltrans(bendersscip, probdeltransBendersPICEF) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(bendersscip, SCIP_OBJSENSE_MINIMIZE) );

   /* create problem data */
   SCIP_CALL( bendersdataPICEFCreate(bendersscip, &bendersdata, graph, graph->nnodes, cycles, posarcs,
         NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, adversarybound, NULL, NULL, NULL, NULL, policy, masterscip, mastersol,
         NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialBendersVars(bendersscip, bendersdata) );

   /* create initial constraints */
   if( policy == POLICY_KEEPUNAFFECTEDCC )
      SCIP_CALL( SCIPcreateInitialBendersConstraintsKUCC(bendersscip, bendersdata) );
   else
      SCIP_CALL( SCIPcreateInitialBendersConstraints(bendersscip, bendersdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(bendersscip, bendersdata) );

   return SCIP_OKAY;
}

/** adds arc variable and its corresponding constraint to the bendersdata */
SCIP_RETCODE SCIPbendersdataPICEFAddArcVars(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 kepscip,            /**< sub SCIP instance for solving 3rd stage */
   SCIP_SOL*             kepsol,             /**< solution of the sub SCIP instance */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Bool             mastersolution      /**< Are we adding a master solution (true) or a kep solution (false)? */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_PROBDATA* bendersdata;
   SCIP_PROBDATA* problemdata; /* Either masterdata or kepdata, depending on whether a master solution is used or not. */
   SCIP_VAR** uvars;
   SCIP_VAR** arcvars;
   SCIP_VAR** yvarcovered;
   int* initposarcs;
   int ninitposarcs;
   int nnodes;
   int npairs;
   int cnt;
   int c;
   int i;
   int j;
   int k;
   int idx;
   int policy;
   int startindex;
   int endindex;
   PositionedArcs* posarcs;
   SCIP_Bool initialchainarc;

   assert( bendersscip != NULL );
   assert( kepscip != NULL );
   assert( vals != NULL );
   assert( vars != NULL );

   bendersdata = SCIPgetProbData(bendersscip);
   problemdata = SCIPgetProbData(kepscip);

   SCIP_CALL( SCIPgetIntParam(bendersscip, "kidney/recoursepolicy", &policy) );

   assert( bendersdata != NULL );
   assert( problemdata != NULL );

   uvars = SCIPbendersdataPICEFGetUVars(bendersdata);
   nnodes = SCIPbendersdataPICEFGetNumNodes(bendersdata);
   npairs = SCIPbendersdataPICEFGetNumPairs(bendersdata);

   if( policy == POLICY_KEEPUNAFFECTEDCC )
   {
      yvarcovered = SCIPbendersdataPICEFGetYvarcovered(bendersdata);
      initposarcs = SCIPbendersdataPICEFGetInitposarcs(bendersdata);
      ninitposarcs = SCIPbendersdataPICEFGetNinitposarcs(bendersdata);
   }

   posarcs = bendersdata->posarcs;

   if ( mastersolution )
      arcvars = masterPICEFProblemGetArcvarinit(problemdata);
   else
      arcvars = SCIPKEPdataPICEFGetArcvars(problemdata);

   /* check if enough memory for arc variables in benders master problem is left */
   if( bendersdata->maxnarcvars <= bendersdata->narcvars + nnodes )
   {
      int newsize;
      if ( bendersdata->maxnarcvars == 0 )
         newsize = nnodes;
      else
         newsize = bendersdata->maxnarcvars * 2;
      SCIP_CALL( SCIPreallocBlockMemoryArray(bendersscip, &bendersdata->arcvars, bendersdata->maxnarcvars, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(bendersscip, &bendersdata->arcindices, bendersdata->maxnarcvars, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(bendersscip, &bendersdata->posarcconss, bendersdata->maxnarcvars, newsize) );
      bendersdata->maxnarcvars = newsize;
   }

   startindex = bendersdata->narcvars;
   endindex = bendersdata->narcvars;

   for (k = 1; k <= posarcs->npositions; ++k)
   {
      for (c = posarcs->positionbegins[k-1]; c < posarcs->positionbegins[k]; c += 2)
      {
         if ( SCIPgetSolVal(kepscip, kepsol, arcvars[c/2]) > 0.5 )
         {
            /* create variable indexed by arc, position and scenario */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arcvar_%d", bendersdata->narcvars);
            SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->arcvars[bendersdata->narcvars], name, 0.0, 1.0, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->arcvars[bendersdata->narcvars]) );
            bendersdata->arcindices[bendersdata->narcvars] = c/2;

            /* create constraint that imposes limitations on being able to use an arc in a chain */
            cnt = 0;

            vals[cnt] = -1.0;
            vars[cnt++] = bendersdata->arcvars[bendersdata->narcvars];

            if (k == 1)
            {
               vals[cnt] = 1.0;
               vars[cnt++] = uvars[posarcs->nodelists[c]];
            }
            else
            {
               for (i = startindex; i < endindex; ++i)
               {
                  if (posarcs->nodelists[2*bendersdata->arcindices[i]+1] == posarcs->nodelists[c])
                  {
                     vals[cnt] = 1.0;
                     vars[cnt++] = bendersdata->arcvars[i];
                     break;
                  }
               }
            }

            vals[cnt] = 1.0;
            vars[cnt++] = uvars[posarcs->nodelists[c+1]];

            /* TODO: check if arc is used in an initial chain with the same subchain leading up to it. If not, we need to add yvarcovered variables in the POLICY_KEEPUNAFFECTEDCC case */
            if( policy == POLICY_KEEPUNAFFECTEDCC )
            {
               /* Can stop whenever c/2 is not an initial position indexed arc at all */
               initialchainarc = FALSE;
               for( j = 0; j < ninitposarcs; ++j )
               {
                  if( initposarcs[j] == c/2 )
                  {
                     initialchainarc = TRUE;
                     idx = j;
                     break;
                  }
               }

               /* Otherwise, check preceding arcs recursively until first vertex is an NDD */
               while( initialchainarc )
               {
                  if( posarcs->nodelists[2*initposarcs[idx]] >= npairs )
                     /* current arc is at position 1, since its source is a NDD */
                     break;

                  for( j = 0; j < idx; ++j )
                  {
                     /* Search j such that initposarcs[j] precedes initposarcs[idx] in initially selected chain */
                     if( posarcs->nodelists[2*initposarcs[idx]] == posarcs->nodelists[2*initposarcs[j] + 1] )
                     {
                        if( SCIPgetSolVal(kepscip, kepsol, arcvars[initposarcs[j]]) > 0.5 )
                           idx = j;
                        else
                           initialchainarc = FALSE;
                        break;
                     }
                  }
               }

               if( !initialchainarc )
               {
                  if( k == 1 )
                  {
                     vals[cnt] = -1.0;
                     vars[cnt++] = yvarcovered[posarcs->nodelists[c]];
                  }

                  vals[cnt] = -1.0;
                  vars[cnt++] = yvarcovered[posarcs->nodelists[c+1]];
               }
            }

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "posarcconss_%d", bendersdata->narcvars);

            SCIP_CALL( SCIPcreateConsBasicLinear(bendersscip, &bendersdata->posarcconss[bendersdata->narcvars], name, cnt,
                  vars, vals, -SCIPinfinity(bendersscip), 1) );

            SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->posarcconss[bendersdata->narcvars]) );
            bendersdata->narcvars++;
         }
      }
      startindex = endindex;
      endindex = bendersdata->narcvars;
   }
   return SCIP_OKAY;
}

/** adds solution constraint to the bendersdata */
SCIP_RETCODE SCIPbendersdataPICEFAddSolCons(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   int                   nvars               /**< number of variables to add */
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
      if (bendersdata->maxnsolconss == 0)
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
   // SCIPinfoMessage(bendersscip, NULL, "added solution constraint to bendersdata; nsolconss = %d\n", bendersdata->nsolconss);

   return SCIP_OKAY;
}

/** returns array of cycle constraints */
SCIP_CONS** SCIPbendersdataPICEFGetCyclelbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->cyclelbconss;
}

/** returns array of cycle variables */
SCIP_VAR** SCIPbendersdataPICEFGetCyclevars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->cyclevars;
}

/** returns number of cycles */
int SCIPbendersdataPICEFGetNCycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->cycles->ncycles;
}

/** returns array of arc constraints */
SCIP_CONS** SCIPbendersdataPICEFGetArcconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->posarcconss;
}

/** returns array of arc variables */
SCIP_VAR** SCIPbendersdataPICEFGetArcvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->arcvars;
}

/** returns number of arc variables */
int SCIPbendersdataPICEFGetNArcvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->narcvars;
}

/** returns array of arc indices */
int* SCIPbendersdataPICEFGetArcIndices(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->arcindices;
}

/** returns number of position indexed arcs */
int SCIPbendersdataPICEFGetNPosarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->posarcs->nposarcs;
}

/** returns uvars */
SCIP_VAR** SCIPbendersdataPICEFGetUVars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->uvars;
}

/** returns yvarcovered variables */
SCIP_VAR** SCIPbendersdataPICEFGetYvarcovered(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->probdatakucc->yvarcovered;
}

/** returns array of initial cycles */
int* SCIPbendersdataPICEFGetInitcycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->probdatakucc->initcycles;
}

/** returns the number of initial cycles */
int SCIPbendersdataPICEFGetNinitcycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->probdatakucc->ninitcycles;
}

/** returns array of initial position indexed arcs */
int* SCIPbendersdataPICEFGetInitposarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->probdatakucc->initposarcs;
}

/** returns number of initial position indexed arcs */
int SCIPbendersdataPICEFGetNinitposarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->probdatakucc->ninitposarcs;
}

/** returns budget of adversary */
int SCIPbendersdataPICEFGetAdversaryBound(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert( bendersdata != NULL );

   return bendersdata->adversarybound;
}

/** returns number of nodes in instance */
int SCIPbendersdataPICEFGetNumNodes(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert( bendersdata != NULL );

   return bendersdata->graph->nnodes;
}

/** returns number of pairs in instance */
int SCIPbendersdataPICEFGetNumPairs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert( bendersdata != NULL );

   return bendersdata->graph->npairs;
}

/** returns objvar */
SCIP_VAR* SCIPbendersdataPICEFGetObjvar(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert( bendersdata != NULL );

   return bendersdata->objvar;
}

/**@} */
