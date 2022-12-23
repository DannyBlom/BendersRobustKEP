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

/**@file   problem_master_kidneyexchange_picef.h
 * @brief  Problem data for the master kidney exchange problem with position indexed cycle edge formulation
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 */

#include <scip/scip.h>
#include "probdata_master_kidney_picef.h"
#include "problem_master_kidneyexchange_picef.h"


/** creates initial model of the master kidney exchange problem */
SCIP_RETCODE SCIPcreateMasterPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to arc structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   SCIP_Bool printinfo;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( posarcs != NULL );
   assert( adversarybound >= 0 );

   SCIP_CALL( SCIPgetBoolParam(scip,"kidney/printgraphinfo", &printinfo) );

   /* get information on graph */
   if ( printinfo )
   {
      int nnodes;
      int ndonors;
      int npairs;
      int narcs;
      int u;
      int v;
      Node** nodelist;
      int* adjacencylists;
      int* adjacencylistbegins;

      nnodes = graph->nnodes;
      npairs = graph->npairs;
      ndonors = graph->ndonors;
      narcs = graph->narcs;
      adjacencylists = graph->adjacencylists;
      adjacencylistbegins = graph->adjacencylistbegins;
      nodelist = graph->nodelist;

      SCIPinfoMessage(scip, NULL, "Create model for graph with %d nodes (%d pairs, %d donors) and %d arcs.\n\n",
         nnodes, npairs, ndonors, narcs);
      SCIPinfoMessage(scip, NULL, "Arc list:\n");

      for (u = 0; u < nnodes; ++u)
      {
         SCIPinfoMessage(scip, NULL, "node %d: ", u);
         for (v = adjacencylistbegins[u]; v < adjacencylistbegins[u + 1]; ++v)
            SCIPinfoMessage(scip, NULL, "%d, ", adjacencylists[v]);
         SCIPinfoMessage(scip, NULL, "\n");
      }

      SCIPinfoMessage(scip, NULL, "\nNode information:\n");
      for (u = 0; u < nnodes; ++u)
         SCIPinfoMessage(scip, NULL, "node %5d: \tID %5d,\t fail probability %1.3f,\t is donor pair %d\n",
            u, nodelist[u]->id, nodelist[u]->failprob, nodelist[u]->ispair);
      SCIPinfoMessage(scip, NULL, "\n");
   }
   SCIP_CALL( SCIPmasterPICEFProbdataCreate(scip, "name", graph, cycles, posarcs, adversarybound) );

   return SCIP_OKAY;
}

/** free master kidney exchange problem data */
SCIP_RETCODE SCIPfreeMasterPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   PositionedArcs*       posarcs             /**< arc structures of graph */
   )
{
   int i;
   int nnodes;

   assert( scip != NULL );
   assert( graph != NULL );

   nnodes = graph->nnodes;

   for (i = nnodes - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemory(scip, &graph->nodelist[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &graph->nodelist, nnodes);

   SCIPfreeBlockMemoryArrayNull(scip, &graph->adjacencylistbegins, nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->adjacencylists, graph->narcs);

   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cycles, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cyclesbegin, nnodes + 1);

   SCIPfreeBlockMemoryArrayNull(scip, &cycles->nodelists, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &cycles->nodelistsbegin, cycles->ncycles + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &cycles->cycleweights, cycles->ncycles);

   SCIPfreeBlockMemoryArrayNull(scip, &posarcs->nodelists, 2*posarcs->nposarcs);
   SCIPfreeBlockMemoryArrayNull(scip, &posarcs->positionbegins, posarcs->npositions + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &posarcs->arcweights, posarcs->nposarcs);

   SCIPfreeBlockMemory(scip, &graph);
   SCIPfreeBlockMemory(scip, &cycles);
   SCIPfreeBlockMemory(scip, &posarcs);

   return SCIP_OKAY;
}
