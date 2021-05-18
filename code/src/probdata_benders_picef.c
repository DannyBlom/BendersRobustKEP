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
   SCIP_CONS**           cycleconss;         /**< cycle constraints (cyclevars[c] = 1 if and only if
                                              *   all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           posarcconss;        /**< arc constraints (arcvar[a,k] = 1 if and only if all uvars[v] = 1
                                              *   for v in arc a and there is a preceding arc used)*/
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
   SCIP_CONS**           cycleconss,         /**< cycle constraints (cyclevars[c] = 1 if and only if
                                              *   all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           posarcconss         /**< arc varbound constraints (arcvar[a,k] = 1 if and only if all uvars[v] = 1
                                              *   for v in arc a and there is a preceding arc used)*/
   )
{
   int ncycles;

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

   if ( cycleconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->cycleconss, cycleconss, ncycles) );
   }
   else
      (*bendersdata)->cycleconss = NULL;

   if ( posarcconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->posarcconss, posarcconss, maxnarcvars) );
   }
   else
      (*bendersdata)->posarcconss = NULL;

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

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( (*bendersdata)->graph != NULL );
   assert( (*bendersdata)->cycles != NULL );

   nnodes = (*bendersdata)->nnodes;
   ncycles = (*bendersdata)->cycles->ncycles;
   narcvars = (*bendersdata)->narcvars;
   maxnarcvars = (*bendersdata)->maxnarcvars;
   nsolconss = (*bendersdata)->nsolconss;
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
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->cycleconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->cycleconss, ncycles);

   for (i = 0; i < narcvars; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->posarcconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->posarcconss, maxnarcvars);

   SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->attackboundcons) );

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
   int i;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( bendersdata->graph != NULL );
   assert( bendersdata->graph->nodelist != NULL );

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
         bendersdata->uvars, vals, nnodes - bendersdata->adversarybound, nnodes,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->attackboundcons) );
   /* do not release constraint here, will be done later */

   /* create constraints linking x_c and its corresponding u_v's */
   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->cycleconss), ncycles) );
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

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclecons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cycleconss[c], name, cnt + 1,
            vars, vals, -1.0, cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cycleconss[c]) );

      /* Reset value of vals[cnt] for next cons */
      vals[cnt] = 1.0;
   }

   SCIPfreeBufferArray(bendersscip, &vars);
   SCIPfreeBufferArray(bendersscip, &vals);

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

   /* create transform bendersdata */
   SCIP_CALL( bendersdataPICEFCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->posarcs,
         sourcedata->cyclevars, sourcedata->arcvars, sourcedata->uvars, sourcedata->objvar, sourcedata->nsolconss, sourcedata->narcvars,
         sourcedata->maxnsolconss, sourcedata->maxnarcvars, sourcedata->arcindices, sourcedata->adversarybound, sourcedata->solconss,
         sourcedata->attackboundcons, sourcedata->cycleconss, sourcedata->posarcconss) );

   nnodes = sourcedata->graph->nnodes;
   ncycles = sourcedata->cycles->ncycles;
   assert( nnodes > 0 );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nsolconss, (*targetdata)->solconss, (*targetdata)->solconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->attackboundcons, &(*targetdata)->attackboundcons) );
   SCIP_CALL( SCIPtransformConss(scip, ncycles, (*targetdata)->cycleconss, (*targetdata)->cycleconss) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->narcvars, (*targetdata)->posarcconss, (*targetdata)->posarcconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->cyclevars, (*targetdata)->cyclevars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->narcvars, (*targetdata)->arcvars, (*targetdata)->arcvars) );
   SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->uvars, (*targetdata)->uvars) );
   SCIP_CALL( SCIPtransformVar(scip, (*targetdata)->objvar, &(*targetdata)->objvar) );

   return SCIP_OKAY;
}

// /** frees user data of transformed Benders problem (called when the transformed problem is freed) */
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
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* bendersdata;

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
         NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, adversarybound, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialBendersVars(bendersscip, bendersdata) );

   /* create initial constraints */
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
   int nnodes;
   int c;
   int i;
   int k;
   int startindex;
   int endindex;
   PositionedArcs* posarcs;

   assert( bendersscip != NULL );
   assert( kepscip != NULL );
   assert( vals != NULL );
   assert( vars != NULL );

   bendersdata = SCIPgetProbData(bendersscip);
   problemdata = SCIPgetProbData(kepscip);

   assert( bendersdata != NULL );
   assert( problemdata != NULL );

   uvars = SCIPbendersdataPICEFGetUVars(bendersdata);
   nnodes = SCIPbendersdataPICEFGetNumNodes(bendersdata);
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

   vals[0] = -1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;

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
            vars[0] = bendersdata->arcvars[bendersdata->narcvars];
            if (k == 1)
               vars[1] = uvars[posarcs->nodelists[c]];
            else
            {
               for (i = startindex; i < endindex; ++i)
               {
                  if (posarcs->nodelists[2*bendersdata->arcindices[i]+1] == posarcs->nodelists[c])
                  {
                     vars[1] = bendersdata->arcvars[i];
                     break;
                  }
               }
            }
            vars[2] = uvars[posarcs->nodelists[c+1]];
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "posarcconss_%d", bendersdata->narcvars);

            SCIP_CALL( SCIPcreateConsBasicLinear(bendersscip, &bendersdata->posarcconss[bendersdata->narcvars], name, 3,
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
   SCIPinfoMessage(bendersscip, NULL, "added solution constraint to bendersdata; nsolconss = %d\n", bendersdata->nsolconss);

   return SCIP_OKAY;
}

/** returns array of cycle constraints */
SCIP_CONS** SCIPbendersdataPICEFGetCycleconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->cycleconss;
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

/** returns objvar */
SCIP_VAR* SCIPbendersdataPICEFGetObjvar(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
)
{
   assert( bendersdata != NULL );

   return bendersdata->objvar;
}

/**@} */
