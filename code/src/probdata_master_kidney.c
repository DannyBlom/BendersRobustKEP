/**@file   probdata_master_kidney.c
 * @brief  Problem data for master kidney exchange problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 *
 * @page KIDNEY_MASTER_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * A list of all interface methods can be found in probdata_master_kidney.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_master_kidney.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"

#include "typedefs.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of master the kidney exchange problem, all variables which are created,
 * and all constraints.
 */
struct SCIP_ProbData
{
   Graph*                graph;              /**< pointer to underlying graph */
   int                   nnodes;             /**< number of nodes in graph */
   int                   npairs;             /**< number of pairs in graph */
   int                   adversarybound;     /**< bound on adversary attack */

   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   int                   ncycles;            /**< number of cycles in graph */

   Chains*               chains;             /**< pointer to chain structure of underlying graph */
   int                   nchains;            /**< number of chains in graph */

   int*                  scenarios;          /**< array encoding attack scenarios */
   int                   nscenarios;         /**< number of scenarios encoded in scenarios array */
   int                   nmaxscenarios;      /**< maximum number of scenarios that can be encoded in scenarios array */

   SCIP_VAR*             objvar;             /**< variable modeling the objective */
   SCIP_VAR**            xvarinit;           /**< variables for assignments before attack */
   SCIP_VAR***           xvarscenario;       /**< variables for assignments after attack per scenario */
   SCIP_VAR***           yvarscenario;       /**< variables indicating surviving assignments per scenario */
   SCIP_VAR**            dummyyvars;         /**< dummy variables for initial problem */

   SCIP_CONS**           objboundconss;      /**< constraints bounding the objective */
   SCIP_CONS**           boundxinitconss;    /**< constraint bounding sum of initial xvars */
   SCIP_CONS***          boundxscenarioconss; /**< constraint bounding sum of scenario based xvars */
   SCIP_CONS***          boundyinitconss;    /**< constraint bounding scenario based yvars by init xvars */
   SCIP_CONS***          boundyscenarioconss; /**< constraint bounding scenario based yvars */
   SCIP_CONS**           dummyconss;         /**< dummy constraints for initial problem */
   SCIP_CONS*            dummyobjcons;       /**< dummy constraints for initial objective bound problem */
};



/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE masterProbdataCreate(
   SCIP*                 scip,                 /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,             /**< pointer to problem data */
   Graph*                graph,                /**< underlying graph */
   int                   nnodes,               /**< number of nodes in graph */
   int                   npairs,               /**< number of pairs in graph */
   int                   adversarybound,       /**< upper bound on number of attacks */
   Cycles*               cycles,               /**< cycle structure of graph */
   Chains*               chains,               /**< pointer to chain structures of graph */
   int*                  scenarios,            /**< array encoding attack scenarios */
   int                   nscenarios,           /**< number of scenarios encoded in scenarios array */
   int                   nmaxscenarios,        /**< maximum number of scenarios that can be encoded in scenarios array */

   SCIP_VAR*             objvar,               /**< variable modeling the objective */
   SCIP_VAR**            xvarinit,             /**< variables for assignments before attack */
   SCIP_VAR***           xvarscenario,         /**< variables for assignments after attack per scenario */
   SCIP_VAR***           yvarscenario,         /**< variables indicating surviving assignments per scenario */

   SCIP_CONS**           objboundconss,        /**< constraints bounding the objective */
   SCIP_CONS**           boundxinitconss,      /**< constraint bounding sum of initial xvars */
   SCIP_CONS***          boundxscenarioconss,  /**< constraint bounding sum of scenario based xvars */
   SCIP_CONS***          boundyinitconss,      /**< constraint bounding scenario based yvars by init xvars */
   SCIP_CONS***          boundyscenarioconss,  /**< constraint bounding scenario based yvars */

   SCIP_VAR**            dummyyvars,           /**< dummy variables for initial problem */
   SCIP_CONS**           dummyconss,           /**< dummy constraints for initial problem */
   SCIP_CONS*            dummyobjcons          /**< dummy constraints for initial objective bound problem */
   )
{
   int ncycles;
   int nchains;
   int ncycleschains;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( adversarybound >= 0 );
   assert( nnodes == graph->nnodes );
   assert( nnodes > 0 );
   assert( npairs >= 0 );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   (*probdata)->graph = graph;
   (*probdata)->nnodes = nnodes;
   (*probdata)->npairs = npairs;
   (*probdata)->adversarybound = adversarybound;
   (*probdata)->cycles = cycles;
   (*probdata)->ncycles = cycles->ncycles;
   ncycles = cycles->ncycles;

   (*probdata)->chains = chains;
   (*probdata)->nchains = chains->nchains;
   nchains = chains->nchains;
   ncycleschains = ncycles + nchains;

   (*probdata)->nscenarios = nscenarios;
   (*probdata)->nmaxscenarios = nmaxscenarios;

   if ( scenarios != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->scenarios, scenarios,
            adversarybound * nmaxscenarios) );
   }
   else
      (*probdata)->scenarios = NULL;

   /* possible copy variable arrays */
   if ( xvarinit != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xvarinit, xvarinit, ncycleschains) );
   }
   else
      (*probdata)->xvarinit = NULL;

   if ( xvarscenario != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->xvarscenario, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xvarscenario[i],
               xvarscenario[i], ncycleschains) );
      }
   }
   else
      (*probdata)->xvarscenario = NULL;

   if ( yvarscenario != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->yvarscenario, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->yvarscenario[i],
               yvarscenario[i], npairs) );
      }
   }
   else
      (*probdata)->yvarscenario = NULL;

   if ( dummyyvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->dummyyvars,
            dummyyvars, nnodes) );
   }
   else
      (*probdata)->dummyyvars = NULL;

   (*probdata)->objvar = objvar;

   /* duplicate arrays */
   if ( objboundconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->objboundconss,
            objboundconss, nmaxscenarios) );
   }
   else
      (*probdata)->objboundconss = NULL;

   if ( boundxinitconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->boundxinitconss,
            boundxinitconss, nnodes) );
   }
   else
      (*probdata)->boundxinitconss = NULL;

   if ( boundxscenarioconss != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->boundxscenarioconss, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->boundxscenarioconss[i],
               boundxscenarioconss[i], nnodes) );
      }
   }
   else
      (*probdata)->boundxscenarioconss = NULL;

   if ( boundyscenarioconss != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->boundyscenarioconss, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->boundyscenarioconss[i],
               boundyscenarioconss[i], npairs) );
      }
   }
   else
      (*probdata)->boundyscenarioconss = NULL;

   if ( boundyinitconss != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->boundyinitconss, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->boundyinitconss[i],
               boundyinitconss[i], npairs) );
      }
   }
   else
      (*probdata)->boundyinitconss = NULL;

   if ( dummyconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->dummyconss,
            dummyconss, nnodes) );
   }
   else
      (*probdata)->dummyconss = NULL;

   if ( dummyobjcons != NULL )
      (*probdata)->dummyobjcons = dummyobjcons;
   else
      (*probdata)->dummyobjcons = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE masterProbdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;
   int j;
   int nnodes;
   int npairs;
   int ncycles;
   int nchains;
   int ncycleschains;
   int nscenarios;
   int nmaxscenarios;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->graph != NULL );
   assert( (*probdata)->cycles != NULL );

   nnodes = (*probdata)->nnodes;
   npairs = (*probdata)->npairs;
   ncycles = (*probdata)->ncycles;
   nchains = (*probdata)->nchains;
   nscenarios = (*probdata)->nscenarios;
   nmaxscenarios = (*probdata)->nmaxscenarios;

   assert( nnodes > 0 );
   assert( npairs >= 0 );
   assert( ncycles >= 0 );
   assert( nchains >= 0 );
   assert( nscenarios >= 0 );
   assert( nscenarios <= nmaxscenarios );

   ncycleschains = ncycles + nchains;

   /* release all variables */
   for (i = 0; i < ncycleschains; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->xvarinit[i]) );
   }
   for (i = 0; i < nscenarios; ++i)
   {
      for (j = 0; j < ncycleschains; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->xvarscenario[i][j]) );
      }
   }
   for (i = 0; i < nscenarios; ++i)
   {
      for (j = 0; j < npairs; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->yvarscenario[i][j]) );
      }
   }
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->dummyyvars[i]) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->objvar) );

   /* release constraints */
   for (i = 0; i < nscenarios; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->objboundconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->objboundconss, nmaxscenarios);

   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundxinitconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundxinitconss, nnodes);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      for (j = 0; j < nnodes; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundxscenarioconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundxscenarioconss[i], nnodes);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundxscenarioconss, nmaxscenarios);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      for (j = 0; j < npairs; ++j)
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundyscenarioconss[i][j]) );
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyscenarioconss[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyscenarioconss, nmaxscenarios);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      for (j = 0; j < npairs; ++j)
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundyinitconss[i][j]) );
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyinitconss[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyinitconss, nmaxscenarios);

   for (i = nnodes - 1; i >= 0; --i)
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->dummyconss[i]) );

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->dummyconss, nnodes);
   SCIPreleaseCons(scip, &(*probdata)->dummyobjcons);

   /* free memory of arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->xvarinit, ncycleschains);
   for (i = nscenarios - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->xvarscenario[i], ncycleschains);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->xvarscenario, nmaxscenarios);
   for (i = nscenarios - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->yvarscenario[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->yvarscenario, nmaxscenarios);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->dummyyvars, nnodes);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->scenarios, (*probdata)->adversarybound * nmaxscenarios);

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** returns whether a given cycle is attacked */
SCIP_Bool SCIPcycleIsAttacked(
   int*                  attackpattern,      /**< array of attacked vertices */
   int                   nattacks,           /**< number of attacked vertices */
   Cycles*               cycles,             /**< data structure for the cycles in the graph of suitable length */
   int                   cycleidx            /**< index of cycle in cycles data structure */
   )
{
   int i;
   int j;

   for (i = 0; i < nattacks; ++i)
   {
      for (j = cycles->nodelistsbegin[cycleidx]; j < cycles->nodelistsbegin[cycleidx+1]; ++j)
      {
         if (attackpattern[i] == cycles->nodelists[j])
            return TRUE;
      }
   }
   return FALSE;
}

/** returns whether a given chain is attacked */
SCIP_Bool SCIPchainIsAttacked(
   int*                  attackpattern,      /**< array of attacked vertices */
   int                   nattacks,           /**< number of attacked vertices */
   Chains*               chains,             /**< data structure for the chains in the graph of suitable length */
   int                   chainidx,           /**< index of chain in chains data structure */
   int                   vertex              /**< Consider only chain up until vertex (if -1, consider entire chain) */
   )
{
   int i;
   int j;
   int v;

   if (vertex == -1)
      /* if no vertex is specified, consider entire chain */
      vertex = chains->nodelists[chains->nodelistsbegin[chainidx + 1] - 1];

   for (j = chains->nodelistsbegin[chainidx]; j < chains->nodelistsbegin[chainidx + 1]; ++j)
   {
      v = chains->nodelists[j];
      for (i = 0; i < nattacks; ++i)
      {
         if (attackpattern[i] == v)
            return TRUE;
      }
      if (v == vertex)
         return FALSE;
   }   
   return FALSE;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialVars(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nvars;
   SCIP_Real nnodes;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );

   nvars = probdata->chains->nchains + probdata->cycles->ncycles;
   nnodes = (SCIP_Real) probdata->graph->nnodes;

   /* create xvarinit */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xvarinit, nvars) );
   for (i = 0; i < nvars; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "xvarinit_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->xvarinit[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->xvarinit[i]) );
   }

   /* create objective variable */
   SCIP_CALL( SCIPcreateVar(scip, &probdata->objvar, "master_objvar", 0.0, nnodes, 1.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, probdata->objvar) );

   /* create dummy variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->dummyyvars, nnodes) );
   for (i = 0; i < nnodes; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyy_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->dummyyvars[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->dummyyvars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the initial constraints of the problem */
static
SCIP_RETCODE SCIPcreateInitialConstraints(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int npairs;
   int maxnvars;
   int ncycles;
   int* node2cycles;
   int* node2cyclesbegin;
   int* node2chains;
   int* node2chainsbegin;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   int i;
   int j;
   int cnt;
   int nmaxvars = 0;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );

   nnodes = probdata->graph->nnodes;
   npairs = probdata->graph->npairs;
   maxnvars = 2*MAX(nnodes + 1, probdata->ncycles + probdata->nchains + 1);
   ncycles = probdata->ncycles;

   node2cycles = probdata->graph->node2cycles;
   node2cyclesbegin = probdata->graph->node2cyclesbegin;
   node2chains = probdata->graph->node2chains;
   node2chainsbegin = probdata->graph->node2chainsbegin;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundxinitconss), nnodes) );

   /* add bound x init cons for each cycle */
   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;
      for (j = node2cyclesbegin[i]; j < node2cyclesbegin[i + 1]; ++j)
         vars[cnt++] = probdata->xvarinit[node2cycles[j]];
      for (j = node2chainsbegin[i]; j < node2chainsbegin[i + 1]; ++j)
         vars[cnt++] = probdata->xvarinit[node2chains[j] + ncycles];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "xvarinit_cons_%d", i);
      SCIP_CALL( SCIPcreateConsSetpack(scip, &(probdata->boundxinitconss[i]), name, cnt, vars,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundxinitconss[i]) );
   }

   /* create dummy constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->dummyconss), nnodes) );
   vals[0] = -1.0;
   for (i = 0; i < nnodes; ++i)
   {
      vars[0] = probdata->dummyyvars[i];

      cnt = 1;
      for (j = node2cyclesbegin[i]; j < node2cyclesbegin[i + 1]; ++j)
      {
         vars[cnt] = probdata->xvarinit[node2cycles[j]];
         vals[cnt++] = 1.0;
      }
      for (j = node2chainsbegin[i]; j < node2chainsbegin[i + 1]; ++j)
      {
         vars[cnt] = probdata->xvarinit[node2chains[j] + ncycles];
         vals[cnt++] = 1.0;
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummy_cons_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->dummyconss[i]), name, cnt,
            vars, vals, 0.0, 0.0) );
      SCIP_CALL( SCIPaddCons(scip, probdata->dummyconss[i]) );

      if ( nmaxvars < cnt )
         nmaxvars = cnt;
   }

   /* create dummy objective constraint */
   for (i = nmaxvars; i < nnodes + 1; ++i)
      vals[i] = 1.0;

   vars[0] = probdata->objvar;
   for (i = 0; i < npairs; ++i)
      vars[i + 1] = probdata->dummyyvars[i];

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->dummyobjcons), "dummyobjcons",
         npairs + 1, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, probdata->dummyobjcons) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** creates variables for a new attack pattern */
static
SCIP_RETCODE SCIPcreateVariablesAttackpattern(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  attackpattern,      /**< array containing nodes being attacked */
   int                   nattacks            /**< number of attacks encoded in attackpattern */
   )
{
   SCIP_PROBDATA* probdata;
   char name[SCIP_MAXSTRLEN];
   Graph* graph;
   int nscenarios;
   int ncycleschains;
   int npairs;
   int c;
   int i;
   int j;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks >= 0 );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->graph != NULL );
   assert( probdata->cycles != NULL );
   assert( probdata->chains != NULL );

   nscenarios = probdata->nscenarios;
   ncycleschains = probdata->cycles->ncycles + probdata->chains->nchains;
   npairs = probdata->graph->npairs;
   graph = probdata->graph;

   /* add variables for each cycle and chain */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xvarscenario[nscenarios], ncycleschains) );
   for (c = 0; c < ncycleschains; ++c)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "xscenario_%d_%d", nscenarios, c);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->xvarscenario[nscenarios][c], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->xvarscenario[nscenarios][c]) );
   }

   /* add variables for each pair */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->yvarscenario[nscenarios], npairs) );
   for (c = 0; c < npairs; ++c)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yscenario_%d_%d", nscenarios, c);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->yvarscenario[nscenarios][c], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->yvarscenario[nscenarios][c]) );
   }

   /* update variable bounds according to attack pattern*/
   for (i = 0; i < nattacks; ++i)
   {
      /* fix upper bound of y to 0 if a pair node is attacked */
      if ( attackpattern[i] < npairs )
      {
         SCIP_CALL( SCIPchgVarUb(scip, probdata->yvarscenario[nscenarios][attackpattern[i]], 0.0) );
      }

      /* fix upper bound of cycle variables to 0 if the cycle contains an attacked node */
      for (c = graph->node2cyclesbegin[attackpattern[i]]; c < graph->node2cyclesbegin[attackpattern[i] + 1]; ++c)
      {
         j = graph->node2cycles[c];
         SCIP_CALL( SCIPchgVarUb(scip, probdata->xvarscenario[nscenarios][j], 0.0) );
      }

      /* fix upper bound of chain variables to 0 if the cycle contains an attacked node */
      for (c = graph->node2chainsbegin[attackpattern[i]]; c < graph->node2chainsbegin[attackpattern[i] + 1]; ++c)
      {
         j = probdata->cycles->ncycles + graph->node2chains[c];
         SCIP_CALL( SCIPchgVarUb(scip, probdata->xvarscenario[nscenarios][j], 0.0) );
      }
   }

   return SCIP_OKAY;
}


/** creates constraints for a new attack pattern */
static
SCIP_RETCODE SCIPcreateConstraintsAttackpattern(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  attackpattern,      /**< array containing nodes being attacked */
   int                   nattacks            /**< number of attacks encoded in attackpattern */
   )
{
   SCIP_PROBDATA* probdata;
   char name[SCIP_MAXSTRLEN];
   Graph* graph;
   Cycles* cycles;
   Chains* chains;
   int nscenarios;
   int ncycleschains;
   int nnodes;
   int npairs;
   int c;
   int i;
   int j;
   int cnt;
   int successor;
   int nmaxvars;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool attacked_node;
   int policy;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks >= 0 );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->graph != NULL );
   assert( probdata->cycles != NULL );
   assert( probdata->chains != NULL );

   graph = probdata->graph;
   cycles = probdata->cycles;
   chains = probdata->chains;
   nscenarios = probdata->nscenarios;
   ncycleschains = cycles->ncycles + chains->nchains;
   nnodes = probdata->nnodes;
   npairs = graph->npairs;

   SCIP_CALL( SCIPgetIntParam(scip, "kidney/recoursepolicy", &policy) );

   nmaxvars = 2*MAX(nnodes, ncycleschains) + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nmaxvars) );
   vals[0] = -1.0;
   for (c = 1; c < nmaxvars; ++c)
      vals[c] = 1.0;

   /* create objective constraint */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "objscenario_%d", nscenarios);
   vars[0] = probdata->objvar;
   for (c = 0; c < npairs; ++c)
      vars[c + 1] = probdata->yvarscenario[nscenarios][c];
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->objboundconss[nscenarios], name,
         npairs + 1, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, probdata->objboundconss[nscenarios]) );

   /* create linking cons between new yvars and init xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundyinitconss[nscenarios]), npairs) );
   for (c = 0; c < npairs; ++c)
   {
      cnt = 0;
      vars[cnt++] = probdata->yvarscenario[nscenarios][c];

      for (i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i)
      {
         j = graph->node2cycles[i];
         vars[cnt++] = probdata->xvarinit[j];
      }
      for (i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i)
      {
         j = probdata->cycles->ncycles + graph->node2chains[i];
         vars[cnt++] = probdata->xvarinit[j];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_y_init_cons_%d_%d", nscenarios, c);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->boundyinitconss[nscenarios][c], name,
         cnt, vars, vals, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundyinitconss[nscenarios][c]) );
   }

   /* create linking cons between new yvars and new xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundyscenarioconss[nscenarios]), npairs) );
   for (c = 0; c < npairs; ++c)
   {
      cnt = 0;
      vars[cnt++] = probdata->yvarscenario[nscenarios][c];

      if( policy == POLICY_KEEPUNAFFECTEDCC )
      {
         /* Consider the terms for all enforceable initial exchanges where j receives a transplant
          * these consist of all non-interdicted cycles and chains with no interdictions in the subchain up to and including j */
         for( i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i )
         {
            j = graph->node2cycles[i];
            if( !SCIPcycleIsAttacked(attackpattern, nattacks, cycles, j) )
               vars[cnt++] = probdata->xvarinit[j];
         }

         for( i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i )
         {
            j = graph->node2chains[i];
            if( !SCIPchainIsAttacked(attackpattern, nattacks, chains, j, c) )
               vars[cnt++] = probdata->xvarinit[cycles->ncycles + j];
         }
      }

      for (i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i)
      {
         j = graph->node2cycles[i];
         vars[cnt++] = probdata->xvarscenario[nscenarios][j];
      }
      for (i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i)
      {
         j = probdata->cycles->ncycles + graph->node2chains[i];
         vars[cnt++] = probdata->xvarscenario[nscenarios][j];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_y_scenario_cons_%d_%d", nscenarios, c);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->boundyscenarioconss[nscenarios][c], name,
         cnt, vars, vals, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundyscenarioconss[nscenarios][c]) );
   }

   /* create matching constraints for new xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundxscenarioconss[nscenarios]), nnodes) );
   vals[0] = 1.0;

   for (c = 0; c < nnodes; ++c)
   {
      attacked_node = FALSE;
      for( i = 0; i < nattacks; ++i )
      {
         if( attackpattern[i] == c )
         {
            attacked_node = TRUE;
            break;
         }
      }

      cnt = 0;

      if( c < npairs )
      {
         /* Consider the terms for all enforceable initial exchanges where j receives a transplant
          * these consist of all non-interdicted cycles and chains with no interdictions in the subchain up to and including j */
         if( policy == POLICY_KEEPUNAFFECTEDCC && !attacked_node )
         {
            for( i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i )
            {
               j = graph->node2cycles[i];
               if( !SCIPcycleIsAttacked(attackpattern, nattacks, cycles, j) )
                  vars[cnt++] = probdata->xvarinit[j];
            }

            for( i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i )
            {
               j = graph->node2chains[i];
               if( !SCIPchainIsAttacked(attackpattern, nattacks, chains, j, c) )
                  vars[cnt++] = probdata->xvarinit[cycles->ncycles + j];
            }
         }
         
         /* Furthermore, add the variables for each possible scenario cycle and chain
          * If not feasible, these variables are already fixed to zero */
         for( i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i )
         {
            j = graph->node2cycles[i];
            vars[cnt++] = probdata->xvarscenario[nscenarios][j];
         }

         for( i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i )
         {
            j = probdata->cycles->ncycles + graph->node2chains[i];
            vars[cnt++] = probdata->xvarscenario[nscenarios][j];
         }         
      }
      else
      {
         /* In this case, c is a non-directed donor.
          * That means we have to consider
          * (i) initial chains such that c and its direct predecessor are not interdicted
          * (ii) all scenario chain variables involving NDD c (if interdicted, it is fixed to zero already anyway) */
         
         if( policy == POLICY_KEEPUNAFFECTEDCC && !attacked_node)
         {
            for( i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i )
            {
               j = graph->node2chains[i];
               successor = chains->nodelists[chains->nodelistsbegin[j]+1];
               /* We only enforce part of chain if both !attacked_node and its direct successor in the chain (first recipient-donor pair) */

               if( !SCIPchainIsAttacked(attackpattern, nattacks, chains, j, successor) )
                  vars[cnt++] = probdata->xvarinit[cycles->ncycles + j];
            }
         }
         
         for( i = graph->node2chainsbegin[c]; i < graph->node2chainsbegin[c + 1]; ++i )
         {
            j = probdata->cycles->ncycles + graph->node2chains[i];
            vars[cnt++] = probdata->xvarscenario[nscenarios][j];
         }
      }
      

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_x_cons_%d_%d", nscenarios, c);
      
      if( !attacked_node )
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->boundxscenarioconss[nscenarios][c], name, cnt, vars, vals, 0.0, 1.0) );
      else
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->boundxscenarioconss[nscenarios][c], name, cnt, vars, vals, 0.0, 0.0) );
   
      SCIP_CALL( SCIPaddCons(scip, probdata->boundxscenarioconss[nscenarios][c]) );
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigMasterKidney)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( masterProbdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransMasterKidney)
{
   int i;
   int nnodes;
   int npairs;
   int ncycleschains;
   int nscenarios;

   /* create transform probdata */
   SCIP_CALL( masterProbdataCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->npairs, sourcedata->adversarybound,
         sourcedata->cycles, sourcedata->chains, sourcedata->scenarios, sourcedata->nscenarios, sourcedata->nmaxscenarios,
         sourcedata->objvar, sourcedata->xvarinit, sourcedata->xvarscenario, sourcedata->yvarscenario,
         sourcedata->objboundconss, sourcedata->boundxinitconss, sourcedata->boundxscenarioconss,
         sourcedata->boundyinitconss, sourcedata->boundyscenarioconss, sourcedata->dummyyvars, sourcedata->dummyconss, sourcedata->dummyobjcons) );

   nnodes = sourcedata->graph->nnodes;
   npairs = sourcedata->graph->npairs;
   ncycleschains = sourcedata->cycles->ncycles + sourcedata->chains->nchains;
   nscenarios = sourcedata->nscenarios;

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, nscenarios, (*targetdata)->objboundconss, (*targetdata)->objboundconss) );
   SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->boundxinitconss, (*targetdata)->boundxinitconss) );
   for (i = 0; i < nscenarios; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->boundxscenarioconss[i], (*targetdata)->boundxscenarioconss[i]) );
      SCIP_CALL( SCIPtransformConss(scip, npairs, (*targetdata)->boundyinitconss[i], (*targetdata)->boundyinitconss[i]) );
      SCIP_CALL( SCIPtransformConss(scip, npairs, (*targetdata)->boundyscenarioconss[i], (*targetdata)->boundyscenarioconss[i]) );
   }
   SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->dummyconss, (*targetdata)->dummyconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->dummyobjcons, &(*targetdata)->dummyobjcons) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVar(scip, (*targetdata)->objvar, &(*targetdata)->objvar) );
   SCIP_CALL( SCIPtransformVars(scip, ncycleschains, (*targetdata)->xvarinit, (*targetdata)->xvarinit) );
   for (i = 0; i < nscenarios; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, ncycleschains, (*targetdata)->xvarscenario[i], (*targetdata)->xvarscenario[i]) );
      SCIP_CALL( SCIPtransformVars(scip, npairs, (*targetdata)->yvarscenario[i], (*targetdata)->yvarscenario[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->dummyyvars, (*targetdata)->dummyyvars) );


   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransMasterKidney)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( masterProbdataFree(scip, probdata) );

   return SCIP_OKAY;
}


/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPmasterProbdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigMasterKidney) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransMasterKidney) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransMasterKidney) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create problem data */
   SCIP_CALL( masterProbdataCreate(scip, &probdata, graph, graph->nnodes, graph->npairs, adversarybound, cycles, chains,
         NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialVars(scip, probdata) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateInitialConstraints(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   return SCIP_OKAY;
}


/** updates master problem by a new attack pattern */
SCIP_RETCODE SCIPupdateMasterProblem(
   SCIP*                 scip,               /**< SCIP data structure  */
   int*                  attackpattern,      /**< array containing nodes being attacked */
   int                   nattacks            /**< number of attacks encoded in attackpattern */
   )
{
   SCIP_PROBDATA* probdata;
   int adversarybound;
   int i;
   int pos;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks >= 0 );

   /* ensure that new pattern fits into data structures */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   adversarybound = probdata->adversarybound;

   if ( probdata->nscenarios == probdata->nmaxscenarios )
   {
      if ( probdata->nmaxscenarios == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->scenarios, adversarybound * probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xvarscenario, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->yvarscenario, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundxscenarioconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundyscenarioconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundyinitconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->objboundconss, probdata->nnodes) );


         probdata->nmaxscenarios = probdata->nnodes;
      }
      else
      {
         int newsize = 2 * probdata->nmaxscenarios;

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->scenarios,
               adversarybound * probdata->nmaxscenarios, adversarybound * newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->xvarscenario, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->yvarscenario, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundxscenarioconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundyscenarioconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundyinitconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->objboundconss, probdata->nmaxscenarios, newsize) );

         probdata->nmaxscenarios = newsize;
      }
   }
   assert( probdata->nscenarios < probdata->nmaxscenarios );

   /* store attack pattern */
   pos = adversarybound * probdata->nscenarios;
   for (i = 0; i < adversarybound; ++i, ++pos)
   {
      if ( i < nattacks )
         probdata->scenarios[pos] = attackpattern[i];
      else
         probdata->scenarios[pos] = -1;
   }

   /* create new variables */
   SCIP_CALL( SCIPcreateVariablesAttackpattern(scip, attackpattern, nattacks) );

   /* create new constraints */
   SCIP_CALL( SCIPcreateConstraintsAttackpattern(scip, attackpattern, nattacks) );

   probdata->nscenarios += 1;

   return SCIP_OKAY;
}

/** get objvar */
SCIP_VAR* masterProblemGetObjvar(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert ( probdata != NULL );
   return probdata->objvar;
}

/** get xvarinit array */
SCIP_VAR** masterProblemGetXvarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->xvarinit;
}

/** get dummyyvars array */
SCIP_VAR** masterProblemGetDummyYVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->dummyyvars;
}

/** get attack scenarios array */
int* masterProblemGetScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->scenarios;
}

/** get number of attack scenarios */
int masterProblemGetNScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->nscenarios;
}

/** get cycles of instance */
Cycles* masterProblemGetCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->cycles;
}

/** get chains of instance */
Chains* masterProblemGetChains(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->chains;
}

/**@} */
