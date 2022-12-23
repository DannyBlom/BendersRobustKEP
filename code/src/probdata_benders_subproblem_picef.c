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

/**@file   probdata_benders_subproblem_picef.c
 * @brief  Problem data for the subproblem of the Benders type approach with the PICEF formulation (KEP restricted to subgraph)
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * A list of all interface methods can be found in probdata_benders_subproblem_picef.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "graph.h"
#include "probdata_benders_subproblem_picef.h"

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
   PositionedArcs*       posarcs;            /**< pointer to positioned arc structure of underlying graph */
   int                   nnodes;             /**< number of nodes in the graph */
   int                   narcs;              /**< number of arcs in the graph */
   int                   nprecconss;         /**< number of precedence constraints */

   SCIP_VAR**            cyclevars;          /**< variable array for all cycles */
   SCIP_VAR**            arcvars;            /**< variable array for all solution based arcs */
   SCIP_VAR**            dummyarcvars;       /**< variable array for all arcs, describing whether the arc part
                                              *   of a subsolution on G[V\u] */

   SCIP_CONS**           nodeconss;          /**< constraint array for each node in the original KEP graph
                                              *   representing a patient-donor pair */
   SCIP_CONS**           arcprecedenceconss; /**< Constraint array of linear inequalities, each corresponding to using
                                              *   an arc from i to j at pos k only if ingoing arc at i pos k-1 is used */
   SCIP_CONS**           dummyarcprecedenceconss; /**< Idem, but for dummy arc variables */
   SCIP_CONS**           dummyarcusedconss;  /**< constraint array: can only choose dummy arc
                                              *   if at any position k its respective posarc is chosen */
};

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE KEPdataPICEFCreate(
   SCIP*                 kepscip,            /**< SCIP instance of 3rd stage */
   SCIP_PROBDATA**       kepdata,            /**< problem data of 3rd stage */
   Graph*                graph,              /**< underlying graph */
   int                   nnodes,             /**< number of nodes in graph */
   int                   narcs,              /**< number of arcs in graph */
   int                   nprecconss,         /**< number of predecessor constraints */
   Cycles*               cycles,             /**< cycle structure of graph */
   PositionedArcs*       posarcs,            /**< positioned arcs in graph */
   SCIP_VAR**            cyclevars,          /**< array of cyclevars: per solution, we keep the cycle
                                              *   if none of its vertices is attacked */
   SCIP_VAR**            arcvars,            /**< array of arcvars: per solution, we keep the arc
                                              *   if none of its vertices is attacked and preceding arcs exist */
   SCIP_VAR**            dummyarcvars,       /**< array of dummy arcvars: per solution, we keep the dummy arc
                                              *   if the arc is chosen and part of the subsolution on G[V\u] */
   SCIP_CONS**           nodeconss,          /**< array of linear inequalities, each corresponding to using
                                              *   at most one cycle/chain per vertex */
   SCIP_CONS**           arcprecedenceconss, /**< array of linear inequalities, each corresponding to using an arc
                                              *   from i to j at pos k only if ingoing arc at i pos k-1 is used */
   SCIP_CONS**           dummyarcprecedenceconss,  /**< array of linear inequalities, each corresponding to using
                                              *   an arc from i to j at pos k only if ingoing arc at i pos k-1 is used */
   SCIP_CONS**           dummyarcusedconss   /**< array of linear inequalities, each corresponding to using a dummy
                                              *   arc only if the arc is chosen at some position */
   )
{
   int npairs;
   int ncycles;
   int nposarcs;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
   assert( nnodes > 0 );
   assert( nnodes == graph->nnodes );
   assert( narcs > 0 );
   assert( narcs == graph->narcs );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(kepscip, kepdata) );

   (*kepdata)->graph = graph;
   (*kepdata)->cycles = cycles;
   (*kepdata)->posarcs = posarcs;
   (*kepdata)->nnodes = nnodes;
   (*kepdata)->narcs = narcs;
   (*kepdata)->nprecconss = graph->npairs*(posarcs->npositions-1);

   ncycles = cycles->ncycles;
   nposarcs = posarcs->nposarcs;
   npairs = graph->npairs;

   if ( cyclevars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->cyclevars, cyclevars, ncycles) );
   }
   else
      (*kepdata)->cyclevars = NULL;

   if ( arcvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->arcvars, arcvars, nposarcs) );
   }
   else
      (*kepdata)->arcvars = NULL;

   if ( dummyarcvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->dummyarcvars, dummyarcvars, narcs) );
   }
   else
      (*kepdata)->dummyarcvars = NULL;

   /* allocate memory for constraints */
   if ( nodeconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->nodeconss, nodeconss, nnodes) );
   }
   else
      (*kepdata)->nodeconss = NULL;

   if ( arcprecedenceconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->arcprecedenceconss, arcprecedenceconss, nprecconss) );
   }
   else
      (*kepdata)->arcprecedenceconss = NULL;

   if ( dummyarcprecedenceconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->dummyarcprecedenceconss, dummyarcprecedenceconss, npairs) );
   }
   else
      (*kepdata)->dummyarcprecedenceconss = NULL;

   if ( dummyarcusedconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(kepscip, &(*kepdata)->dummyarcusedconss, dummyarcusedconss, narcs) );
   }
   else
      (*kepdata)->dummyarcusedconss = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given KEP problem data */
static
SCIP_RETCODE KEPdataPICEFFree(
   SCIP*                 kepscip,            /**< SCIP data structure */
   SCIP_PROBDATA**       kepdata             /**< pointer to problem data */
   )
{
   int nnodes;
   int narcs;
   int ncycles;
   int nposarcs;
   int nprecconss;
   int npairs;
   int i;
   SCIP_Bool lifting;


   assert( kepscip != NULL );
   assert( kepdata != NULL );
   SCIP_CALL( SCIPgetBoolParam(kepscip, "kidney/liftbenderscuts", &lifting) );

   nnodes = (*kepdata)->nnodes;
   narcs = (*kepdata)->narcs;
   npairs = (*kepdata)->graph->npairs;
   nprecconss = (*kepdata)->nprecconss;
   ncycles = (*kepdata)->cycles->ncycles;
   nposarcs = (*kepdata)->posarcs->nposarcs;

   assert( nnodes > 0 );
   assert( narcs >= 0 );
   assert( ncycles >= 0);
   assert( nprecconss >= 0);
   assert( nposarcs >= 0);

   /* release all variables */
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(kepscip, &(*kepdata)->cyclevars[i]) );
   }

   for (i = 0; i < nposarcs; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(kepscip, &(*kepdata)->arcvars[i]) );
   }

   if ( lifting )
   {
      for (i = 0; i < narcs; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(kepscip, &(*kepdata)->dummyarcvars[i]) );
      }
   }

   /** release constraints and free constraint memory arrays */
   for (i = 0; i < nnodes; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(kepscip, &(*kepdata)->nodeconss[i]) );
   }

   for (i = 0; i < nprecconss; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(kepscip, &(*kepdata)->arcprecedenceconss[i]) );
   }

   if ( lifting )
   {
      for (i = 0; i < npairs; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(kepscip, &(*kepdata)->dummyarcprecedenceconss[i]) );
      }

      for (i = 0; i < narcs; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(kepscip, &(*kepdata)->dummyarcusedconss[i]) );
      }
   }

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(kepscip, &(*kepdata)->cyclevars, ncycles);
   SCIPfreeBlockMemoryArray(kepscip, &(*kepdata)->arcvars, nposarcs);
   SCIPfreeBlockMemoryArrayNull(kepscip, &(*kepdata)->nodeconss, nnodes);
   SCIPfreeBlockMemoryArrayNull(kepscip, &(*kepdata)->arcprecedenceconss, nprecconss);


   if ( lifting )
   {
      SCIPfreeBlockMemoryArray(kepscip, &(*kepdata)->dummyarcvars, narcs);
      SCIPfreeBlockMemoryArrayNull(kepscip, &(*kepdata)->dummyarcprecedenceconss, narcs);
      SCIPfreeBlockMemoryArrayNull(kepscip, &(*kepdata)->dummyarcusedconss, narcs);
   }

   /* free kepdata block memory */
   SCIPfreeBlockMemory(kepscip, kepdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialKEPVars(
   SCIP*                 kepscip,            /**< SCIP instance of 3rd stage */
   SCIP_PROBDATA*        kepdata,            /**< problem data of 3rd stage */
   SCIP_Bool             lifting             /**< do we want to lift solutions to entire graph (yes/no)? */
   )
{
   char name[SCIP_MAXSTRLEN];
   int ncycles;
   int narcs;
   int nnodes;
   int nposarcs;
   int npositions;
   int i, j, k;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( kepdata->graph != NULL );
   assert( kepdata->graph->nodelist != NULL );

   ncycles = kepdata->cycles->ncycles;
   nposarcs = kepdata->posarcs->nposarcs;
   npositions = kepdata->posarcs->npositions;
   narcs = kepdata->narcs;
   nnodes = kepdata->nnodes;

   /* create x_c-variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &kepdata->cyclevars, ncycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &kepdata->arcvars, nposarcs) );
   if ( lifting )
      SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &kepdata->dummyarcvars, narcs) );

   for (i = 0; i < ncycles; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cyclevar_%d", i);

      SCIP_CALL( SCIPcreateVar(kepscip, &kepdata->cyclevars[i], name, 0.0, 1.0, kepdata->cycles->cycleweights[i],
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(kepscip, kepdata->cyclevars[i]) );
   }

   for( k = 0; k < npositions; ++k )
   {
      for( i = kepdata->posarcs->positionbegins[k]; i < kepdata->posarcs->positionbegins[k+1]; i += 2 )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arcvar_(%d,%d,%d)", kepdata->posarcs->nodelists[i], kepdata->posarcs->nodelists[i+1], k + 1);

         SCIP_CALL( SCIPcreateVar(kepscip, &kepdata->arcvars[i/2], name, 0.0, 1.0, kepdata->posarcs->arcweights[i/2],
               SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(kepscip, kepdata->arcvars[i/2]) );
      }
   }

   if ( lifting )
   {
      for( i = 0; i < nnodes; ++i )
      {
         for( j = kepdata->graph->adjacencylistbegins[i]; j < kepdata->graph->adjacencylistbegins[i + 1]; ++j )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyarcvar_(%d, %d)", i, kepdata->graph->adjacencylists[j]);
            SCIP_CALL( SCIPcreateVar(kepscip, &kepdata->dummyarcvars[j], name, 0.0, 1.0, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(kepscip, kepdata->dummyarcvars[j]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** creates the initial constraints of the problem */
static
SCIP_RETCODE SCIPcreateInitialKEPConstraints(
   SCIP*                 kepscip,            /**< SCIP instance of 3rd stage */
   SCIP_PROBDATA*        kepdata,            /**< problem data of 3rd stage */
   SCIP_Bool             lifting             /**< do we want to lift solutions to entire graph (yes/no)? */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nnodes;
   int npairs;
   int narcs;
   int ncycles;
   int nposarcs;
   Graph* graph;
   Cycles* cycles;
   PositionedArcs* posarcs;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   int maxnvarsincons;
   int i;
   int c;
   int j;
   int k;
   int n;
   int cnt;
   int index;

   assert( kepscip != NULL );
   assert( kepdata != NULL );
   assert( kepdata->graph != NULL );

   nnodes = kepdata->graph->nnodes;
   npairs = kepdata->graph->npairs;
   narcs = kepdata->graph->narcs;
   graph = kepdata->graph;
   assert ( graph->adjacencylistbegins[nnodes] == narcs);
   cycles = kepdata->cycles;
   ncycles = cycles->ncycles;
   assert( kepdata->cycles->ncycles == cycles->ncycles );

   posarcs = kepdata->posarcs;
   nposarcs = posarcs->nposarcs;
   assert( kepdata->posarcs->nposarcs == posarcs->nposarcs );
   maxnvarsincons = ncycles + MAX(narcs, nposarcs);

   /* allocate buffer memory arrays to add conss to the Benders type model */
   SCIP_CALL( SCIPallocBufferArray(kepscip, &vars, maxnvarsincons) );
   SCIP_CALL( SCIPallocBufferArray(kepscip, &vals, maxnvarsincons) );

   /* create node constraints (only one cycle or chain per node) */
   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &(kepdata->nodeconss), nnodes) );

   for (i = 0; i < maxnvarsincons; ++i)
      vals[i] = 1.0;

   for (i = 0; i < nnodes; ++i)
   {
      cnt = 0;
      if ( i < npairs )
      {
         for (c = graph->node2cyclesbegin[i]; c < graph->node2cyclesbegin[i + 1]; ++c)
            vars[cnt++] = kepdata->cyclevars[graph->node2cycles[c]];

         for (c = 0; c < nposarcs; ++c)
         {
            /* Consider incoming arcs for patient donor pairs */
            if ( posarcs->nodelists[2*c+1] == i )
               vars[cnt++] = kepdata->arcvars[c];
         }
      }
      else
      {
         for (c = 0; c < nposarcs; ++c)
         {
            /* Consider outgoing arcs for NDDs */
            if ( posarcs->nodelists[2*c] == i )
               vars[cnt++] = kepdata->arcvars[c];
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nodecons_%d", i);

      /* sum of incident cyclevars and arcvars does not exceed 1*/
      SCIP_CALL( SCIPcreateConsLinear(kepscip, &kepdata->nodeconss[i], name, cnt, vars, vals, 0.0, 1.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(kepscip, kepdata->nodeconss[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &(kepdata->arcprecedenceconss), kepdata->nprecconss) );
   index = 0;
   for (k = 2; k <= posarcs->npositions; ++k)
   {
      for (i = 0; i < npairs; ++i)
      {
         cnt = 0;
         for (c = posarcs->positionbegins[k-2]+1; c < posarcs->positionbegins[k-1]; c += 2)
         {

            if (posarcs->nodelists[c] == i)
            {
               vals[cnt] = 1.0;
               vars[cnt++] = kepdata->arcvars[c/2];
            }
         }

         for (c = posarcs->positionbegins[k-1]; c < posarcs->positionbegins[k]; c += 2)
         {
            if (posarcs->nodelists[c] == i)
            {
               vals[cnt] = -1.0;
               vars[cnt++] = kepdata->arcvars[c/2];
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arcprecedencecons_(%d, %d)", i, k);

         /* sum of incident cyclevars and arcvars does not exceed 1*/
         SCIP_CALL( SCIPcreateConsBasicLinear(kepscip, &kepdata->arcprecedenceconss[index], name, cnt,
               vars, vals, 0.0, SCIPinfinity(kepscip)) );

         SCIP_CALL( SCIPaddCons(kepscip, kepdata->arcprecedenceconss[index++]) );
      }
   }

   /** Dummy arc constraints are only necessary whenever we lift solutions to entire compatibility graph */
   if ( lifting )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &(kepdata->dummyarcprecedenceconss), narcs) );
      for (i = 0; i < npairs; ++i)
      {
         cnt = 0;
         for (c = 0; c < i; ++c)
         {
            for (j = graph->adjacencylistbegins[c]; j < graph->adjacencylistbegins[c+1]; ++j)
            {
               if ( graph->adjacencylists[j] == i )
               {
                  /* These are incoming arc in i (from vertex with smaller label than i), so add them with coef 1.0 to constraint */
                  vals[cnt] = 1.0;
                  vars[cnt++] = kepdata->dummyarcvars[j];
               }
            }
         }

         /* These are outgoing arcs from i, so add them with coef -1.0 to constraint */
         for (j = graph->adjacencylistbegins[i]; j < graph->adjacencylistbegins[i+1]; ++j)
         {
            vals[cnt] = -1.0;
            vars[cnt++] = kepdata->dummyarcvars[j];
         }

         for (c = i + 1; c < nnodes; ++c)
         {
            for (j = graph->adjacencylistbegins[c]; j < graph->adjacencylistbegins[c+1]; ++j)
            {
               if ( graph->adjacencylists[j] == i )
               {
                  /* These are incoming arc in i (from vertex with larger label than i), so add them with coef 1.0 to constraint */
                  vals[cnt] = 1.0;
                  vars[cnt++] = kepdata->dummyarcvars[j];
               }
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyarcprecedencecons_%d", i);

         /* sum of incident cyclevars and arcvars does not exceed 1*/
         SCIP_CALL( SCIPcreateConsBasicLinear(kepscip, &kepdata->dummyarcprecedenceconss[i], name,
               cnt, vars, vals, 0.0, SCIPinfinity(kepscip)) );

         SCIP_CALL( SCIPaddCons(kepscip, kepdata->dummyarcprecedenceconss[i]) );
      }


      SCIP_CALL( SCIPallocBlockMemoryArray(kepscip, &(kepdata->dummyarcusedconss), narcs) );
      for (i = 0; i < nnodes; ++i)
      {
         for (j = graph->adjacencylistbegins[i]; j < graph->adjacencylistbegins[i+1]; ++j)
         {
            cnt = 0;
            vals[cnt] = -1.0;
            vars[cnt++] = kepdata->dummyarcvars[j];

            c = graph->adjacencylists[j];

            for (n = 0; n < nposarcs; ++n)
            {
               if ( posarcs->nodelists[2*n] == i && posarcs->nodelists[2*n+1] == c )
               {
                  vals[cnt] = 1.0;
                  vars[cnt++] = kepdata->arcvars[n];
               }
            }

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dummyarcusedcons_%d", j);

            /* dummy arc is used only whenever any of the corresponding arcs is used */
            SCIP_CALL( SCIPcreateConsBasicLinear(kepscip, &kepdata->dummyarcusedconss[j], name, cnt,
                  vars, vals, 0.0, SCIPinfinity(kepscip)) );

            SCIP_CALL( SCIPaddCons(kepscip, kepdata->dummyarcusedconss[j]) );
         }
      }
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
SCIP_DECL_PROBDELORIG(probdelorigKEPPICEF)
{
   SCIPdebugMsg(scip, "free original KEP data\n");

   SCIP_CALL( KEPdataPICEFFree(scip, probdata) );

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPKEPdataPICEFCreate(
   SCIP*                 kepscip,            /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs             /**< pointer to positioned arc structures of graph */
   )
{
   SCIP_PROBDATA* kepdata;
   SCIP_Bool lifting;

   assert( kepscip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(kepscip, probname) );
   SCIP_CALL( SCIPsetProbDelorig(kepscip, probdelorigKEPPICEF) );
   SCIP_CALL( SCIPgetBoolParam(kepscip, "kidney/liftbenderscuts", &lifting) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(kepscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create problem data */
   SCIP_CALL( KEPdataPICEFCreate(kepscip, &kepdata, graph, graph->nnodes, graph->narcs,
         graph->npairs*(posarcs->npositions-1), cycles, posarcs,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create initial variables */
   SCIP_CALL( SCIPcreateInitialKEPVars(kepscip, kepdata, lifting) );

   /* create initial constraints */
   SCIP_CALL( SCIPcreateInitialKEPConstraints(kepscip, kepdata, lifting) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(kepscip, kepdata) );

   return SCIP_OKAY;
}

/** returns the number of nodes in the KEP instance */
int SCIPKEPdataPICEFGetNumNodes(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->graph->nnodes;
}

/** returns the number of pairs in the KEP instance */
int SCIPKEPdataPICEFGetNumPairs(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->graph->npairs;
}

/** returns the number of arcs in the KEP instance */
int SCIPKEPdataPICEFGetNumArcs(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->graph->narcs;
}

/** returns the graph in the KEP instance */
Graph* SCIPKEPdataPICEFGetGraph(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->graph;
}

/** returns array of cycle variables */
SCIP_VAR** SCIPKEPdataPICEFGetCyclevars(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL);

   return kepdata->cyclevars;
}

/** returns number of cycles */
int SCIPKEPdataPICEFGetNCycles(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->cycles->ncycles;
}

/** returns array of arc variables */
SCIP_VAR** SCIPKEPdataPICEFGetArcvars(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL );

   return kepdata->arcvars;
}

/** returns number of posarcvars */
int SCIPKEPdataPICEFGetNPosarcs(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->posarcs->nposarcs;
}

/** returns array of dummy arc variables */
SCIP_VAR** SCIPKEPdataPICEFGetDummyArcvars(
   SCIP_PROBDATA*        kepdata             /**< problem data */
   )
{
   assert ( kepdata != NULL );

   return kepdata->dummyarcvars;
}

/** returns array of node constraints */
SCIP_CONS** SCIPKEPdataPICEFGetNodeconss(
   SCIP_PROBDATA*       kepdata              /**< problem data */
   )
{
   assert( kepdata != NULL );

   return kepdata->nodeconss;
}

/**@} */
