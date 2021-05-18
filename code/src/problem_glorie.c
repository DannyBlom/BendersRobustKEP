/**@file   problem_glorie.h
 * @brief  Problem data for Glorie et al. approach to the robust kidney exchange problem
 * @author Christopher Hojny
 * @author Bart Smeulders
 */

#include <scip/scip.h>
#include "problem_glorie.h"
#include "probdata_glorie.h"
#include "problem_kidneyexchange.h"


/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateGlorieEtAlModel(
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
   SCIP_CALL( SCIPglorieEtAldataCreate(scip, "name", graph, cycles, chains) );

   return SCIP_OKAY;
}


/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeGlorieEtAlModel(
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

/** extends the Glorie attack greedily, and also computes a lower bound based on this attack.*/
SCIP_Real GlorieExtendAttack(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains,             /**< chain structures of graph */
   int*                  usedcycleschains,   /**< array containing indices of cycles and chains used in master solution */
   int                   nusedcycleschains,  /**< number of cycles and chains used in master solution */
   int*                  zerofixings,        /**< array of 0-fixings (identified by variable indices) */
   int*                  onefixings,         /**< array of 1-fixings (identified by variable indices) */
   int                   nzerofixings,       /**< number of 0-fixings */
   int                   nonefixings,        /**< number of 1-fixings */
   int                   attackbound,        /**< upper bound on number of attacks */
   int*                  attackpattern
)
{
   /* We first get the weights of all cycles and chains used in the first stage problem.
    * We also figure out which cycles or chains were already attacked.
    */
   int* usedcyclechainweights = NULL;
   int* cyclechainattacked;
   int i;
   int j;
   int k;
   int l;
   int nodeunfixed = 1;
   int maxunattackedweight = 0;
   int newattackflag = 0;
   int newattackcyclechain = 0;
   int newattacknode = 0;
   SCIP_Real lb = 0.0;

   SCIP_CALL( SCIPallocBufferArray(scip, &usedcyclechainweights, nusedcycleschains) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &cyclechainattacked, nusedcycleschains) );

   for (i = 0; i < nusedcycleschains; ++i)
   {
      if ( usedcycleschains[i] < cycles->ncycles )
      {
         usedcyclechainweights[i] = cycles->cycleweights[usedcycleschains[i]];
         cyclechainattacked[i] = 0;
         for (j = cycles->nodelistsbegin[usedcycleschains[i]]; j < cycles->nodelistsbegin[usedcycleschains[i]+1]; ++j)
         {
            for (k = 0; k < nzerofixings; ++k)
            {
               if ( zerofixings[k] == cycles->nodelists[j] )
               {
                  cyclechainattacked[i] = 1;
                  break;
               }
            }
         }
      }
      else
      {
         usedcyclechainweights[i] = chains->chainweights[usedcycleschains[i] - cycles->ncycles];
         cyclechainattacked[i] = 0;
         for (j = chains->nodelistsbegin[usedcycleschains[i]-cycles->ncycles]; j < chains->nodelistsbegin[usedcycleschains[i]+1-cycles->ncycles]; ++j)
         {
            for (k = 0; k < nzerofixings; ++k)
            {
               if ( zerofixings[k] == chains->nodelists[j] )
               {
                  cyclechainattacked[i] = 1;
                  break;
               }
            }
         }
      }
   }

   /* Next, we attack the highest weight cycles or chains - which have not yet been attacked and for which there
    * is at least one unfixed node. This is limited by either the number of attacks, or the maximum numer
    * of variables not fixed to 1.
    */
   for (i = nzerofixings; i < attackbound && i < graph->nnodes - nonefixings; )
   {
      maxunattackedweight = 0;
      newattackflag = 0;
      newattackcyclechain = 0;
      newattacknode = 0;

      /* Identify the best greedy attack */
      for (j = 0; j < nusedcycleschains; j++)
      {
         if ( usedcyclechainweights[j] > maxunattackedweight && cyclechainattacked[j] == 0 )
         {
            /* If this cycle or chain is interesting to attack, we still need to check whether this is possible
             * (it could be the case that every vertex is fixed to 1 - not attacked.)
             */
            if ( usedcycleschains[j] < cycles->ncycles )
            {
               for (k = cycles->nodelistsbegin[usedcycleschains[j]]; k < cycles->nodelistsbegin[usedcycleschains[j]+1]; ++k)
               {
                  nodeunfixed = 1;

                  /* Note that we only need to check the 1-fixings, if some nodes are fixed to zero,
                   * the cycle or chain woul already have been attacked.
                   */
                  for (l = 0; l < nonefixings && nodeunfixed == 1; ++l)
                  {
                     if ( cycles->nodelists[k] == onefixings[l] )
                     {
                        nodeunfixed = 0;
                        break;
                     }
                  }
                  if ( nodeunfixed == 1 )
                  {
                     newattacknode = cycles->nodelists[k];
                     newattackflag = 1;
                     maxunattackedweight = usedcyclechainweights[j];
                     newattackcyclechain = j;
                     break;
                  }
               }
            }
            else /* Chain */
            {
               for (k = chains->nodelistsbegin[usedcycleschains[j]-cycles->ncycles]; k < chains->nodelistsbegin[usedcycleschains[j]+1-cycles->ncycles]; ++k)
               {
                  nodeunfixed = 1;

                  /* Note that we only need to check the 1-fixings, if some nodes are fixed to zero,
                   * the cycle or chain woul already have been attacked.
                   */
                  for (l = 0; l < nonefixings && nodeunfixed == 1; ++l)
                  {
                     if ( chains->nodelists[k] == onefixings[l] )
                     {
                        nodeunfixed = 0;
                        break;
                     }
                  }
                  if ( nodeunfixed == 1 )
                  {
                     newattacknode = chains->nodelists[k];
                     newattackflag = 1;
                     maxunattackedweight = usedcyclechainweights[j];
                     newattackcyclechain = j;
                     break;
                  }
               }
            }
         }
      }

      /* Add the attack. */
      if( newattackflag == 1 )
      {
           attackpattern[i++] = newattacknode;
           cyclechainattacked[newattackcyclechain] = 1;
      }
      /* None of the used cycles and chains can be attacked. We add lowest index nodes
       * that have not yet been fixed (or are not yet in the attack pattern).
       */
      else
      {
         for (j = 0; i < attackbound; ++j)
         {
            int notfixed = 1;

            for (k = 0; k < nzerofixings; ++k)
            {
               if( zerofixings[k] == j )
                  notfixed = 0;
            }
            for( k = 0; k < nonefixings; ++k)
            {
               if ( onefixings[k] == j )
                  notfixed = 0;
            }

            if ( notfixed == 1 )
               attackpattern[i++] = j;
         }
      }
   }

   /* Compute the lower bound */
   for (i = 0; i < nusedcycleschains; ++i)
   {
      if ( cyclechainattacked[i] == 0 )
         lb += usedcyclechainweights[i];
   }

   SCIPfreeBufferArray(scip, &cyclechainattacked);
   SCIPfreeBufferArray(scip, &usedcyclechainweights);

   return lb;
}
