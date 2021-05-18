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
   SCIP_CONS**           cycleconss;         /**< cycle constraints (cyclevars[c] = 1 if and only if all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           chainconss;         /**< chain constraints (chainvars[c] = 1 if and only if all uvars[v] = 1 for v in chain c)*/
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
   SCIP_CONS**           cycleconss,         /**< cycle constraints (cyclevars[c] = 1 if and only if all uvars[v] = 1 for v in cycle c) */
   SCIP_CONS**           chainconss          /**< chain constraints (chainvars[c] = 1 if and only if all uvars[v] = 1 for v in chain c)*/
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

   if ( chainvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->chainvars, chainvars, nchains) );
   }
   else
      (*bendersdata)->chainvars = NULL;

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

   if ( chainconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(bendersscip, &(*bendersdata)->chainconss, chainconss, nchains) );
   }
   else
      (*bendersdata)->chainconss = NULL;

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
   int nchains;
   int nsolconss;

   assert( bendersscip != NULL );
   assert( bendersdata != NULL );
   assert( (*bendersdata)->graph != NULL );
   assert( (*bendersdata)->cycles != NULL );

   nnodes = (*bendersdata)->nnodes;
   ncycles = (*bendersdata)->cycles->ncycles;
   nchains = (*bendersdata)->chains->nchains;
   nsolconss = (*bendersdata)->nsolconss;
   assert( nnodes > 0 );
   assert( ncycles >= 0 );
   assert( nchains >= 0 );

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

   for (i = 0; i < nchains; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(bendersscip, &(*bendersdata)->chainvars[i]) );
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

   for (i = 0; i < nchains; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->chainconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(bendersscip, &(*bendersdata)->chainconss, nchains);

   SCIP_CALL( SCIPreleaseCons(bendersscip, &(*bendersdata)->attackboundcons) );

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->uvars, nnodes);
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArray(bendersscip, &(*bendersdata)->chainvars, nchains);

   /* free bendersdata */
   SCIPfreeBlockMemory(bendersscip, bendersdata);

   return SCIP_OKAY;
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

   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->cyclevars[i], name, 1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->cyclevars[i]) );
   }

   for (i = 0; i < nchains; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chainvar_%d", i);

      SCIP_CALL( SCIPcreateVar(bendersscip, &bendersdata->chainvars[i], name, 1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(bendersscip, bendersdata->chainvars[i]) );
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
   SCIP*                 bendersscip,        /**< SCIP instance of Benders model */
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int ncycles;
   int nchains;
   Cycles* cycles;
   Chains* chains;
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

   chains = bendersdata->chains;
   nchains = chains->nchains;
   assert( bendersdata->chains->nchains == chains->nchains );

   SCIP_CALL( SCIPgetBoolParam(bendersscip, "kidney/cyclechainconssinitial", &initial) );

   /* Allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vals, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(bendersscip, &vars, maxnvarsincons) );

   /* create bound on adversary attack */
   for (i = 0; i < maxnvarsincons; ++i)
      vals[i] = 1.0;

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
         vars[cnt++] = bendersdata->uvars[cycles->nodelists[i]];

      vars[cnt] = bendersdata->cyclevars[c];
      vals[cnt] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclecons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->cycleconss[c], name, cnt + 1, vars, vals,
            -1.0, cycles->nodelistsbegin[c + 1] - cycles->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->cycleconss[c]) );

      /* Reset value of vals[cnt] for next cons */
      vals[cnt] = 1.0;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(bendersscip, &(bendersdata->chainconss), nchains) );
   for (c = 0; c < nchains; ++c)
   {
      cnt = 0;
      for (i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c + 1]; ++i)
         vars[cnt++] = bendersdata->uvars[chains->nodelists[i]];

      vars[cnt] = bendersdata->chainvars[c];
      vals[cnt] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "chaincons_%d", c);

      /* x_c >= 1 - sum (1 - u_v) */
      SCIP_CALL( SCIPcreateConsLinear(bendersscip, &bendersdata->chainconss[c], name, cnt + 1, vars, vals,
            -1.0, chains->nodelistsbegin[c + 1] - chains->nodelistsbegin[c] - 1,
            initial, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(bendersscip, bendersdata->chainconss[c]) );

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
   int nchains;

   /* create transform bendersdata */
   SCIP_CALL( bendersdataCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->cycles, sourcedata->chains,
         sourcedata->cyclevars, sourcedata->chainvars, sourcedata->uvars, sourcedata->objvar, sourcedata->nsolconss, sourcedata->maxnsolconss,
         sourcedata->adversarybound, sourcedata->solconss, sourcedata->attackboundcons, sourcedata->cycleconss,
         sourcedata->chainconss) );

   nnodes = sourcedata->graph->nnodes;
   ncycles = sourcedata->cycles->ncycles;
   nchains = sourcedata->chains->nchains;
   assert( nnodes > 0 );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nsolconss, (*targetdata)->solconss, (*targetdata)->solconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->attackboundcons, &(*targetdata)->attackboundcons) );
   SCIP_CALL( SCIPtransformConss(scip, ncycles, (*targetdata)->cycleconss, (*targetdata)->cycleconss) );
   SCIP_CALL( SCIPtransformConss(scip, nchains, (*targetdata)->chainconss, (*targetdata)->chainconss) );

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
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* bendersdata;

   assert( bendersscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(bendersscip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(bendersscip, probdelorigBenders) );
   SCIP_CALL( SCIPsetProbTrans(bendersscip, probtransBenders) );
   SCIP_CALL( SCIPsetProbDeltrans(bendersscip, probdeltransBenders) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(bendersscip, SCIP_OBJSENSE_MINIMIZE) );

   /* create problem data */
   SCIP_CALL( bendersdataCreate(bendersscip, &bendersdata, graph, graph->nnodes, cycles, chains,
         NULL, NULL, NULL, NULL, 0, 0, adversarybound, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialBendersVars(bendersscip, bendersdata) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateInitialBendersConstraints(bendersscip, bendersdata) );

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
      if ( bendersdata->maxnsolconss == 0 )
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
SCIP_CONS** SCIPbendersdataGetCycleconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );
   return bendersdata->cycleconss;
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


/** returns array of chain constraints */
SCIP_CONS** SCIPbendersdataGetChainconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   )
{
   assert( bendersdata != NULL );

   return bendersdata->chainconss;
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
