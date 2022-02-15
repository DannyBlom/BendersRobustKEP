/**@file   problem_benders.c
 * @brief  Problem data for Benders type approach to the robust kidney exchange problem
 * @author Danny Blom
 */

#include <scip/scip.h>
#include "problem_benders.h"
#include "probdata_benders.h"
#include "probdata_benders_subproblem.h"
#include "problem_kidneyexchange.h"


/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateBendersModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   SCIP_SOL*             sol,                /**< SCIP solution to master problem */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   )
{
   Node** nodelist;
   int* adjacencylists;
   int* adjacencylistbegins;
   int nnodes;
   int ndonors;
   int npairs;
   int narcs;
   int u;
   int v;
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
   SCIP_CALL( SCIPbendersdataCreate(scip, masterscip, "name", sol, graph, cycles, chains, adversarybound) );

   return SCIP_OKAY;
}

/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateBendersSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   )
{
   Node** nodelist;
   int* adjacencylists;
   int* adjacencylistbegins;
   int nnodes;
   int ndonors;
   int npairs;
   int narcs;
   int u;
   int v;
   SCIP_Bool printinfo;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( cycles != NULL );
   assert( chains != NULL );

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
   SCIP_CALL( SCIPKEPdataCreate(scip, "name", graph, cycles, chains) );

   return SCIP_OKAY;
}

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains              /**< chain structures of graph */
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

   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cyclesbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cycles, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chainsbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chains, chains->nnodesinchains);

   SCIPfreeBlockMemory(scip, &graph);
   SCIPfreeBlockMemory(scip, &cycles);
   SCIPfreeBlockMemory(scip, &chains);

   return SCIP_OKAY;
}

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains              /**< chain structures of graph */
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

   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cyclesbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2cycles, cycles->nnodesincycles);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chainsbegin, graph->nnodes + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graph->node2chains, chains->nnodesinchains);

   SCIPfreeBlockMemory(scip, &graph);
   SCIPfreeBlockMemory(scip, &cycles);
   SCIPfreeBlockMemory(scip, &chains);

   return SCIP_OKAY;
}
