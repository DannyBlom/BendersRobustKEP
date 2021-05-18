/**@file   probdata_master_kidney_picef.c
 * @brief  Problem data for master kidney exchange problem with position indexed cycle edge formulation
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 *
 * @page KIDNEY_MASTER_PICEF_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * A list of all interface methods can be found in probdata_master_kidney_picef.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_master_kidney_picef.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"

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
   int                   adversarybound;     /**< bound on adversaty attack */

   Cycles*               cycles;             /**< pointer to cycle structure of underlying graph */
   int                   ncycles;            /**< number of cycles in graph */

   PositionedArcs*       posarcs;            /**< pointer to positioned arc structure of underlying graph */
   int                   nposarcs;           /**< number of positioned arcs in underlying graph */
   int                   npositions;         /**< Number of positions in which arcs can be in a chain (i.e. maxlenchain) */

   int*                  scenarios;          /**< array encoding attack scenarios */
   int                   nscenarios;         /**< number of scenarios encoded in scenarios array */
   int                   nmaxscenarios;      /**< maximum number of scenarios that can be encoded in scenarios array */

   SCIP_VAR*             objvar;             /**< variable modeling the objective */
   SCIP_VAR**            xcyclevarinit;      /**< variables for assignments before attack */
   SCIP_VAR***           xcyclevarscenario;  /**< variables for assignments after attack per scenario */
   SCIP_VAR**            arcvarinit;         /**< variables for arc assignment before attack */
   SCIP_VAR***           arcvarscenario;     /**< variables for arc assignment after attack per scenario */

   SCIP_VAR***           yvarscenario;       /**< variables indicating surviving assignments per scenario */
   SCIP_VAR**            dummyyvars;         /**< dummy variables for initial problem */

   SCIP_CONS**           objboundconss;      /**< constraints bounding the objective */
   SCIP_CONS**           boundxinitconss;    /**< constraint bounding sum of initial xvars */
   SCIP_CONS***          boundxscenarioconss;/**< constraint bounding sum of scenario based xvars */
   SCIP_CONS***          boundyinitconss;    /**< constraint bounding scenario based yvars by init xvars */
   SCIP_CONS***          boundyscenarioconss;/**< constraint bounding scenario based yvars */
   SCIP_CONS**           piefinitconss;      /**< constraint imposing the use of arcs only with existing preceding arcs */
   SCIP_CONS***          piefscenarioconss;  /**< constraint imposing the use of arcs only with existing preceding arcs */
   SCIP_CONS**           dummyconss;         /**< dummy constraints for initial problem */
   SCIP_CONS*            dummyobjcons;       /**< dummy constraints for initial objective bound problem */
};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE masterPICEFProbdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   int                   npairs,             /**< number of pairs in graph */
   int                   adversarybound,     /**< upper bound on number of attacks */

   Cycles*               cycles,             /**< cycle structure of graph */
   PositionedArcs*       posarcs,            /**< positioned arc structure of graph */
   int                   nposarcs,           /**< Number of position indexed arcs */
   int                   npositions,         /**< Number of positions in which an arc can be in a chain */

   int*                  scenarios,          /**< array encoding attack scenarios */
   int                   nscenarios,         /**< number of scenarios encoded in scenarios array */
   int                   nmaxscenarios,      /**< maximum number of scenarios that can be encoded in scenarios array */

   SCIP_VAR*             objvar,             /**< variable modeling the objective */
   SCIP_VAR**            xvarinit,           /**< variables for cycle assignments before attack */
   SCIP_VAR***           xvarscenario,       /**< variables for cycle assignments after attack per scenario */
   SCIP_VAR**            zarcposvarinit,     /**< variables for arc assignment before attack */
   SCIP_VAR***           zarcposvarscenario, /**< variables for arc assignment after attack per scenario */
   SCIP_VAR***           yvarscenario,       /**< variables indicating surviving assignments per scenario */
   SCIP_VAR**            dummyyvars,         /**< dummy variables for initial problem */

   SCIP_CONS**           objboundconss,      /**< constraints bounding the objective */
   SCIP_CONS**           boundxinitconss,    /**< constraint bounding sum of initial xvars */
   SCIP_CONS***          boundxscenarioconss, /**< constraint bounding sum of scenario based xvars */
   SCIP_CONS***          boundyinitconss,    /**< constraint bounding scenario based yvars by init xvars */
   SCIP_CONS***          boundyscenarioconss, /**< constraint bounding scenario based yvars */
   SCIP_CONS**           piefinitconss,      /**< PICEF constraints imposing use of preceding arc
                                              *   in order to be eligible for use */
   SCIP_CONS***          piefscenarioconss,  /**< PICEF constraints imposing use of preceding arc
                                              *   in order to be eligible for use in each scenario */
   SCIP_CONS**           dummyconss,         /**< dummy constraints for initial problem */
   SCIP_CONS*            dummyobjcons        /**< dummy constraints for initial objective bound problem */
   )
{
   int ncycles;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
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

   (*probdata)->posarcs = posarcs;
   (*probdata)->nposarcs = nposarcs;
   (*probdata)->npositions = npositions;
   (*probdata)->nscenarios = nscenarios;
   (*probdata)->nmaxscenarios = nmaxscenarios;

   if ( scenarios != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->scenarios, scenarios, adversarybound * nmaxscenarios) );
   }
   else
      (*probdata)->scenarios = NULL;

   /* possible copy variable array for cycle variables, both initially and in each scenario */
   if ( xvarinit != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xcyclevarinit, xvarinit, ncycles) );
   }
   else
      (*probdata)->xcyclevarinit = NULL;

   if ( xvarscenario != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->xcyclevarscenario, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xcyclevarscenario[i],
               xvarscenario[i], ncycles) );
      }
   }
   else
      (*probdata)->xcyclevarscenario = NULL;

   /* possible copy variable array for positioned arc variables, both initially and in each scenario */
   if ( zarcposvarinit != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->arcvarinit, zarcposvarinit, nposarcs) );
   }
   else
      (*probdata)->arcvarinit = NULL;

   if ( zarcposvarscenario != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->arcvarscenario, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->arcvarscenario[i],
               zarcposvarscenario[i], nposarcs) );
      }
   }
   else
      (*probdata)->arcvarscenario = NULL;

   /* possible copy variable array for dummy variables, both initial and for each scenario */
   if ( dummyyvars  != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->dummyyvars,
            dummyyvars, nnodes) );
   }
   else
      (*probdata)->dummyyvars = NULL;

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

   /* Set objective variable */
   (*probdata)->objvar = objvar;

   /* duplicate constraint arrays */
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

   if ( piefinitconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->piefinitconss,
            piefinitconss, graph->npairs*(posarcs->npositions-1)) );
   }
   else
      (*probdata)->piefinitconss = NULL;

   if ( piefscenarioconss != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->piefscenarioconss, nmaxscenarios) );
      for (i = 0; i < nscenarios; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->piefscenarioconss[i],
               piefscenarioconss[i], npairs*(npositions-1)) );
      }
   }
   else
      (*probdata)->piefscenarioconss = NULL;

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
SCIP_RETCODE masterPICEFProbdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;
   int j;
   int nnodes;
   int npairs;
   int ncycles;
   int nposarcs;
   int npositions;
   int npiefconss;
   int nscenarios;
   int nmaxscenarios;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->graph != NULL );
   assert( (*probdata)->cycles != NULL );
   assert( (*probdata)->posarcs != NULL );

   nnodes = (*probdata)->nnodes;
   npairs = (*probdata)->npairs;
   ncycles = (*probdata)->ncycles;
   nscenarios = (*probdata)->nscenarios;
   nmaxscenarios = (*probdata)->nmaxscenarios;
   nposarcs = (*probdata)->nposarcs;
   npositions = (*probdata)->npositions;
   npiefconss = npairs*(npositions-1);

   assert( nnodes > 0 );
   assert( npairs >= 0 );
   assert( ncycles >= 0 );
   assert( nscenarios >= 0 );
   assert( nscenarios <= nmaxscenarios );
   assert( nposarcs >= 0);
   assert( npositions >= 0);

   /* release all variables */
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->xcyclevarinit[i]) );
   }
   for (i = 0; i < nscenarios; ++i)
   {
      for (j = 0; j < ncycles; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->xcyclevarscenario[i][j]) );
      }
   }
   for (i = 0; i < nposarcs; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->arcvarinit[i]) );
   }
   for (i = 0; i < nscenarios; ++i)
   {
      for (j = 0; j < nposarcs; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->arcvarscenario[i][j]) );
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
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundyscenarioconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyscenarioconss[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyscenarioconss, nmaxscenarios);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      for (j = 0; j < npairs; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->boundyinitconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyinitconss[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->boundyinitconss, nmaxscenarios);

   for (j = npiefconss-1; j >= 0; --j)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->piefinitconss[j]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->piefinitconss, npiefconss);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      for (j = 0; j < npiefconss; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->piefscenarioconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->piefscenarioconss[i], npiefconss );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->piefscenarioconss, nmaxscenarios);

   for (i = nnodes - 1; i >= 0; --i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->dummyconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->dummyconss, nnodes);
   SCIPreleaseCons(scip, &(*probdata)->dummyobjcons);

   /* free memory of variable arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->xcyclevarinit, ncycles);
   for (i = nscenarios - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->xcyclevarscenario[i], ncycles);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->xcyclevarscenario, nmaxscenarios);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->arcvarinit, nposarcs);
   for (i = nscenarios - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->arcvarscenario[i], nposarcs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->arcvarscenario, nmaxscenarios);

   for (i = nscenarios - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->yvarscenario[i], npairs);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->yvarscenario, nmaxscenarios);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->dummyyvars, nnodes);

   /* free memory of scenario arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->scenarios, (*probdata)->adversarybound * nmaxscenarios);

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialVars(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nposarcs;
   int ncycles;
   SCIP_Real nnodes;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );

   nposarcs = probdata->nposarcs;
   ncycles = probdata->ncycles;
   nnodes = (SCIP_Real) probdata->graph->nnodes;

   /* create initial variables for each cycle */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xcyclevarinit, ncycles) );
   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "xcyclevarinit_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->xcyclevarinit[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->xcyclevarinit[i]) );
   }

   /* create initial variable for each position indexed arcs */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->arcvarinit, nposarcs) );
   for (i = 0; i < nposarcs; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arcvarinit_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->arcvarinit[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->arcvarinit[i]) );
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
   PositionedArcs* posarcs;
   int nnodes;
   int npairs;
   int maxnvars;
   int nposarcs;
   int npositions;
   int npiefconss;
   int* node2cycles;
   int* node2cyclesbegin;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   int i;
   int j;
   int cnt;
   int index;
   int nmaxvars = 0;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->graph != NULL );

   nnodes = probdata->graph->nnodes;
   npairs = probdata->graph->npairs;
   maxnvars = MAX(nnodes + 1, probdata->ncycles + probdata->nposarcs + 1);
   nposarcs = probdata->nposarcs;

   node2cycles = probdata->graph->node2cycles;
   node2cyclesbegin = probdata->graph->node2cyclesbegin;
   posarcs = probdata->posarcs;
   npositions = posarcs->npositions;
   npiefconss = npairs*(npositions-1);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundxinitconss), nnodes) );

   /* add boundxinitcons for each node based on cycle and chain edge vars */
   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;

      if ( i < npairs )
      {
         for (j = node2cyclesbegin[i]; j < node2cyclesbegin[i + 1]; ++j)
            vars[cnt++] = probdata->xcyclevarinit[node2cycles[j]];

         for (j = 0; j < nposarcs; ++j)
         {
            if ( posarcs->nodelists[2*j+1] == i )
               /* then the variable with index j corresponds to an arc going into i */
               vars[cnt++] = probdata->arcvarinit[j];
         }
      }
      else
      {
         for (j = 0; j < nposarcs; ++j)
         {
            if ( posarcs->nodelists[2*j] == i )
               /* then the variable with index j corresponds to an arc going into i */
               vars[cnt++] = probdata->arcvarinit[j];
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundxinitconss_%d", i);
      SCIP_CALL( SCIPcreateConsSetpack(scip, &(probdata->boundxinitconss[i]), name, cnt, vars,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundxinitconss[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->piefinitconss), npiefconss) );

   /* add initial picef constraints allowing the use of arc (u,v) at pos k only
    * with preceding arcs at pos k-1 going into u.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, maxnvars) );
   index = 0;
   for (i = 0; i < npairs; ++i)
   {
      int k;

      for (k = 2; k <= posarcs->npositions; ++k)
      {
         int p;
         cnt = 0;

         /* Loop to get incoming arcs with pos k-1 */
         for (p = posarcs->positionbegins[k-2]+1; p < posarcs->positionbegins[k-1]; p += 2)
         {
            if ( posarcs->nodelists[p] == i)
            {
               vars[cnt] = probdata->arcvarinit[p/2];
               vals[cnt++] = 1.0;
            }
         }
         /* Loop to get outgoing arcs with pos k */
         for (p = posarcs->positionbegins[k-1]; p < posarcs->positionbegins[k]; p += 2)
         {
            if ( posarcs->nodelists[p] == i)
            {
               vars[cnt] = probdata->arcvarinit[p/2];
               vals[cnt++] = -1.0;
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "piefinitcons_%d", index);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->piefinitconss[index]), name, cnt,
               vars, vals, 0.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, probdata->piefinitconss[index++]) );
      }
   }

   /* create dummy constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->dummyconss), nnodes) );
   vals[0] = -1.0;
   for (i = 0; i < nnodes; ++i)
   {
      vars[0] = probdata->dummyyvars[i];
      cnt = 1;
      for (j = node2cyclesbegin[i]; j < node2cyclesbegin[i + 1]; ++j)
      {
         vars[cnt] = probdata->xcyclevarinit[node2cycles[j]];
         vals[cnt++] = 1.0;
      }

      for (j = 0; j < nposarcs; ++j)
      {
         if ( posarcs->nodelists[2 * j + 1] == i )
         {
            /* then the variable with index j corresponds to an arc going into i */
            vars[cnt] = probdata->arcvarinit[j];
            vals[cnt++] = 1.0;
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummy_cons_%d", i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &(probdata->dummyconss[i]), name, cnt, vars, vals, 0.0, 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
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

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->dummyobjcons), "dummyobjcons", npairs + 1,
         vars, vals, 0.0, SCIPinfinity(scip)) );
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
   PositionedArcs* posarcs;
   char name[SCIP_MAXSTRLEN];
   Graph* graph;
   int nscenarios;
   int nposarcs;
   int ncycles;
   int npairs;
   int c;
   int i;
   int j;
   int v;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks >= 0 );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->graph != NULL );
   assert( probdata->cycles != NULL );

   nscenarios = probdata->nscenarios;
   nposarcs = probdata->nposarcs;
   ncycles = probdata->ncycles;
   npairs = probdata->graph->npairs;
   graph = probdata->graph;
   posarcs = probdata->posarcs;

   /* add variables for each cycle */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xcyclevarscenario[nscenarios], ncycles) );
   for (c = 0; c < ncycles; ++c)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "xscenario_%d_%d", nscenarios, c);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->xcyclevarscenario[nscenarios][c], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->xcyclevarscenario[nscenarios][c]) );
   }

   /* add variables for each position indexed arc var */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->arcvarscenario[nscenarios], nposarcs) );
   for (c = 0; c < nposarcs; ++c)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arcvarscenario_%d_%d", nscenarios, c);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->arcvarscenario[nscenarios][c], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->arcvarscenario[nscenarios][c]) );
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
      v = attackpattern[i];

      /* fix upper bound of y to 0 if a pair node is attacked */
      if ( v < npairs )
      {
         SCIP_CALL( SCIPchgVarUb(scip, probdata->yvarscenario[nscenarios][v], 0.0) );
      }

      /* fix upper bound of cycle variables to 0 if the cycle contains an attacked node */
      for (c = graph->node2cyclesbegin[v]; c < graph->node2cyclesbegin[v + 1]; ++c)
      {
         j = graph->node2cycles[c];
         SCIP_CALL( SCIPchgVarUb(scip, probdata->xcyclevarscenario[nscenarios][j], 0.0) );
      }

      /* fix upper bound of cycle variables to 0 if the cycle contains an attacked node */
      for (c = 0; c < nposarcs; ++c)
      {
         if ( posarcs->nodelists[2*c] == v || posarcs->nodelists[2*c+1] == v )
         {
            SCIP_CALL( SCIPchgVarUb(scip, probdata->arcvarscenario[nscenarios][c], 0.0) );
         }
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
   PositionedArcs* posarcs;
   char name[SCIP_MAXSTRLEN];
   Graph* graph;
   int nscenarios;
   int ncycles;
   int nposarcs;
   int nnodes;
   int npairs;
   int c;
   int i;
   int j;
   int cnt;
   int index;
   int nmaxvars;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   assert( scip != NULL );
   assert( attackpattern != NULL );
   assert( nattacks >= 0 );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->graph != NULL );
   assert( probdata->cycles != NULL );

   graph = probdata->graph;
   posarcs = probdata->posarcs;
   nposarcs = probdata->nposarcs;
   nscenarios = probdata->nscenarios;
   ncycles = probdata->cycles->ncycles;
   nnodes = probdata->nnodes;
   npairs = graph->npairs;

   nmaxvars = MAX(nnodes, ncycles+nposarcs) + 1;
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
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->objboundconss[nscenarios]), name,
         npairs + 1, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, probdata->objboundconss[nscenarios]) );

   /* create linking cons between new yvars and init xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundyinitconss[nscenarios]), npairs) );
   for (c = 0; c < npairs; ++c)
   {
      vars[0] = probdata->yvarscenario[nscenarios][c];
      cnt = 1;

      for (i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i)
      {
         j = graph->node2cycles[i];
         vars[cnt++] = probdata->xcyclevarinit[j];
      }
      for (i = 0; i < nposarcs; ++i)
      {
         if ( posarcs->nodelists[2*i+1] == c )
            vars[cnt++] = probdata->arcvarinit[i];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_y_init_cons_%d_%d", nscenarios, c);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->boundyinitconss[nscenarios][c]), name,
         cnt, vars, vals, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundyinitconss[nscenarios][c]) );
   }

   /* create linking cons between new yvars and new xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundyscenarioconss[nscenarios]), npairs) );
   for (c = 0; c < npairs; ++c)
   {
      vars[0] = probdata->yvarscenario[nscenarios][c];
      cnt = 1;

      for (i = graph->node2cyclesbegin[c]; i < graph->node2cyclesbegin[c + 1]; ++i)
      {
         j = graph->node2cycles[i];
         vars[cnt++] = probdata->xcyclevarscenario[nscenarios][j];
      }
      for (i = 0; i < nposarcs; ++i)
      {
         if ( posarcs->nodelists[2*i+1] == c )
            vars[cnt++] = probdata->arcvarscenario[nscenarios][i];
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_y_scenario_cons_%d_%d", nscenarios, c);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->boundyscenarioconss[nscenarios][c]), name,
         cnt, vars, vals, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundyscenarioconss[nscenarios][c]) );
   }

   /* create matching constraints for new xvars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->boundxscenarioconss[nscenarios]), nnodes) );
   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;

      if ( i < npairs )
      {
         for (j = graph->node2cyclesbegin[i]; j < graph->node2cyclesbegin[i + 1]; ++j)
            vars[cnt++] = probdata->xcyclevarscenario[nscenarios][graph->node2cycles[j]];

         for (j = 0; j < nposarcs; ++j)
         {
            /* if the variable with index j corresponds to an arc going into i */
            if ( posarcs->nodelists[2*j+1] == i )
               vars[cnt++] = probdata->arcvarscenario[nscenarios][j];
         }
      }
      else
      {
         for (j = 0; j < nposarcs; ++j)
         {
            /* if the variable with index j corresponds to an arc going into i */
            if ( posarcs->nodelists[2*j] == i )
               vars[cnt++] = probdata->arcvarscenario[nscenarios][j];
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_scenario_x_cons_%d_%d", nscenarios, i);
      SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &(probdata->boundxscenarioconss[nscenarios][i]), name, cnt, vars) );
      SCIP_CALL( SCIPaddCons(scip, probdata->boundxscenarioconss[nscenarios][i]) );
   }

   index = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->piefscenarioconss[nscenarios]), npairs*(posarcs->npositions-1)) );
   for (i = 0; i < npairs; ++i)
   {
      int k;

      for (k = 2; k <= posarcs->npositions; ++k)
      {
         int p;
         cnt = 0;

         /* loop to get incoming arcs with pos k-1 */
         for (p = posarcs->positionbegins[k-2] + 1; p < posarcs->positionbegins[k-1]; p += 2)
         {
            if ( posarcs->nodelists[p] == i)
            {
               vars[cnt] = probdata->arcvarscenario[nscenarios][p/2];
               vals[cnt++] = 1.0;
            }
         }

         /* loop to get outgoing arcs with pos k */
         for (p = posarcs->positionbegins[k-1]; p < posarcs->positionbegins[k]; p+=2)
         {
            if ( posarcs->nodelists[p] == i)
            {
               vars[cnt] = probdata->arcvarscenario[nscenarios][p/2];
               vals[cnt++] = -1.0;
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "piefscenarioconss__%d_%d", nscenarios, index);
         SCIP_CALL( SCIPcreateConsLinear(scip, &(probdata->piefscenarioconss[nscenarios][index]), name, cnt, vars, vals, 0.0, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, probdata->piefscenarioconss[nscenarios][index++]) );
      }
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
SCIP_DECL_PROBDELORIG(probdelorigMasterPICEFKidney)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( masterPICEFProbdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransMasterPICEFKidney)
{
   int i;
   int nnodes;
   int npairs;
   int nposarcs;
   int npositions;
   int ncycles;
   int nscenarios;

   /* create transform probdata */
   SCIP_CALL( masterPICEFProbdataCreate(scip, targetdata, sourcedata->graph, sourcedata->nnodes, sourcedata->npairs,
         sourcedata->adversarybound, sourcedata->cycles, sourcedata->posarcs, sourcedata->nposarcs,
         sourcedata->npositions, sourcedata->scenarios, sourcedata->nscenarios, sourcedata->nmaxscenarios,
         sourcedata->objvar, sourcedata->xcyclevarinit, sourcedata->xcyclevarscenario, sourcedata->arcvarinit,
         sourcedata->arcvarscenario, sourcedata->yvarscenario, sourcedata->dummyyvars, sourcedata->objboundconss,
         sourcedata->boundxinitconss, sourcedata->boundxscenarioconss, sourcedata->boundyinitconss,
         sourcedata->boundyscenarioconss, sourcedata->piefinitconss, sourcedata->piefscenarioconss,
         sourcedata->dummyconss, sourcedata->dummyobjcons) );

   nnodes = sourcedata->graph->nnodes;
   npairs = sourcedata->graph->npairs;
   nposarcs = sourcedata->nposarcs;
   npositions = sourcedata->npositions;
   ncycles = sourcedata->cycles->ncycles;
   nscenarios = sourcedata->nscenarios;

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, nscenarios, (*targetdata)->objboundconss, (*targetdata)->objboundconss) );
   SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->boundxinitconss, (*targetdata)->boundxinitconss) );
   SCIP_CALL( SCIPtransformConss(scip, npairs*(npositions-1), (*targetdata)->piefinitconss, (*targetdata)->piefinitconss) );
   for (i = 0; i < nscenarios; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->boundxscenarioconss[i], (*targetdata)->boundxscenarioconss[i]) );
      SCIP_CALL( SCIPtransformConss(scip, npairs, (*targetdata)->boundyinitconss[i], (*targetdata)->boundyinitconss[i]) );
      SCIP_CALL( SCIPtransformConss(scip, npairs, (*targetdata)->boundyscenarioconss[i], (*targetdata)->boundyscenarioconss[i]) );
      SCIP_CALL( SCIPtransformConss(scip, npairs*(npositions-1), (*targetdata)->piefscenarioconss[i], (*targetdata)->piefscenarioconss[i]) );
   }
   SCIP_CALL( SCIPtransformConss(scip, nnodes, (*targetdata)->dummyconss, (*targetdata)->dummyconss) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->dummyobjcons, &(*targetdata)->dummyobjcons) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVar(scip, (*targetdata)->objvar, &(*targetdata)->objvar) );
   SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->xcyclevarinit, (*targetdata)->xcyclevarinit) );
   SCIP_CALL( SCIPtransformVars(scip, nposarcs, (*targetdata)->arcvarinit, (*targetdata)->arcvarinit) );

   for (i = 0; i < nscenarios; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, ncycles, (*targetdata)->xcyclevarscenario[i], (*targetdata)->xcyclevarscenario[i]) );
      SCIP_CALL( SCIPtransformVars(scip, nposarcs, (*targetdata)->arcvarscenario[i], (*targetdata)->arcvarscenario[i]) );
      SCIP_CALL( SCIPtransformVars(scip, npairs, (*targetdata)->yvarscenario[i], (*targetdata)->yvarscenario[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, nnodes, (*targetdata)->dummyyvars, (*targetdata)->dummyyvars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransMasterPICEFKidney)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( masterPICEFProbdataFree(scip, probdata) );

   return SCIP_OKAY;
}


/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPmasterPICEFProbdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigMasterPICEFKidney) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransMasterPICEFKidney) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransMasterPICEFKidney) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create problem data */
   SCIP_CALL( masterPICEFProbdataCreate(scip, &probdata, graph, graph->nnodes, graph->npairs, adversarybound,
         cycles, posarcs, posarcs->nposarcs, posarcs->npositions,
         NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialVars(scip, probdata) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateInitialConstraints(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   return SCIP_OKAY;
}


/** updates master problem by a new attack pattern */
SCIP_RETCODE SCIPupdateMasterPICEFProblem(
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
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xcyclevarscenario, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->arcvarscenario, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->yvarscenario, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundxscenarioconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundyscenarioconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->boundyinitconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->piefscenarioconss, probdata->nnodes) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->objboundconss, probdata->nnodes) );

         probdata->nmaxscenarios = probdata->nnodes;
      }
      else
      {
         int newsize;
         newsize = 2 * probdata->nmaxscenarios;

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->scenarios,
               adversarybound * probdata->nmaxscenarios, adversarybound * newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->xcyclevarscenario, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->arcvarscenario, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->yvarscenario, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundxscenarioconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundyscenarioconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->boundyinitconss, probdata->nmaxscenarios, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->piefscenarioconss, probdata->nmaxscenarios, newsize) );
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
SCIP_VAR* masterPICEFproblemGetObjvar(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert ( probdata != NULL );
   return probdata->objvar;
}

/** get cycles */
Cycles* masterPICEFProblemGetCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert ( probdata != NULL );
   return probdata->cycles;
}

/** get xvarinit array */
SCIP_VAR** masterPICEFProblemGetXCyclevarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->xcyclevarinit;
}

/** get posarcs */
PositionedArcs* masterPICEFProblemGetPosarcs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->posarcs;
}

/** get arcvars array */
SCIP_VAR** masterPICEFProblemGetArcvarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->arcvarinit;
}

/** get dummyyvars array */
SCIP_VAR** masterPICEFProblemGetDummyYVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->dummyyvars;
}

/** get attack scenarios array */
int* masterPICEFProblemGetScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->scenarios;
}

/** get number of attack scenarios */
int masterPICEFProblemGetNScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->nscenarios;
}

/**@} */
