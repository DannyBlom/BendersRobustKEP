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

#include <scip/scip.h>
#include "pricer_stableSet.h"
#include "probdata_kidney.h"
#include "problem_kidneyexchange.h"


/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   Triplets*             triplets,           /**< pointer to triplet structure of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Bool             userankseparation,  /**< whether rank separation shall be used */
   int                   maxrankrhs          /**< maximum rhs of rank inequalities */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** cycleconss;
   SCIP_CONS** chainconss;
   SCIP_CONS* objlinkingcons;
   int nnodes;
   int ndonors;
   int npairs;
   int narcs;
   int u;
   int v;
   Node** nodelist;
   int* adjacencylists;
   int* adjacencylistbegins;
   SCIP_Bool printinfo;
   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );
   assert( adversarybound >= 0 );
   SCIP_CALL( SCIPgetBoolParam(scip,"kidney/printgraphinfo", &printinfo) );

   /* get information on graph */
   nnodes = graph->nnodes;
   npairs = graph->npairs;
   ndonors = graph->ndonors;
   narcs = graph->narcs;
   adjacencylists = graph->adjacencylists;
   adjacencylistbegins = graph->adjacencylistbegins;
   nodelist = graph->nodelist;
   if ( printinfo )
   {
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
   SCIP_CALL( SCIPprobdataCreate(scip, "name", graph, cycles, chains, triplets, adversarybound) );

   /* activate pricer */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   cycleconss = SCIPprobdataGetCycleconss(probdata);
   assert( cycleconss != NULL );
   assert( SCIPprobdataGetNCycles(probdata) == cycles->ncycles );

   chainconss = SCIPprobdataGetChainconss(probdata);
   assert( chainconss != NULL );
   assert( SCIPprobdataGetNChains(probdata) == chains->nchains );

   objlinkingcons = SCIPprobdataGetObjlinkingcons(probdata);
   assert( objlinkingcons != NULL );

   SCIP_CALL( SCIPpricerStableSetActivate(scip, cycleconss, cycles, chainconss, chains, triplets,
         objlinkingcons, nnodes, userankseparation, maxrankrhs) );

   return SCIP_OKAY;
}

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains,             /**< chain structures of graph */
   Triplets*             triplets            /**< triplet structures of grpah */
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

   SCIPfreeBlockMemoryArrayNull(scip, &graph->adjacencylistbegins, graph->narcs + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->adjacencylists, graph->narcs);

   SCIPfreeBlockMemoryArrayNull(scip, &cycles->nodelists, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &cycles->nodelistsbegin, cycles->ncycles + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &cycles->cycleweights, cycles->ncycles);

   SCIPfreeBlockMemoryArrayNull(scip, &chains->nodelists, chains->nnodesinchains);
   SCIPfreeBlockMemoryArrayNull(scip, &chains->nodelistsbegin, chains->nchains + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &chains->chainweights, chains->nchains);

   SCIPfreeBlockMemoryArrayNull(scip, &triplets->nodelists, 3*triplets->ntriplets);
   SCIPfreeBlockMemoryArrayNull(scip, &triplets->cycle2triplets, triplets->nindicescycle2triplet);
   SCIPfreeBlockMemoryArrayNull(scip, &triplets->cycle2tripletsbegin, cycles->ncycles+1);
   SCIPfreeBlockMemoryArrayNull(scip, &triplets->chain2triplets, triplets->nindiceschain2triplet);
   SCIPfreeBlockMemoryArrayNull(scip, &triplets->chain2tripletsbegin, chains->nchains+1);

   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cyclesbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cycles, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chainsbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chains, chains->nnodesinchains);

   SCIPfreeBlockMemory(scip, &graph);
   SCIPfreeBlockMemory(scip, &cycles);
   SCIPfreeBlockMemory(scip, &chains);
   SCIPfreeBlockMemory(scip, &triplets);

   return SCIP_OKAY;
}
