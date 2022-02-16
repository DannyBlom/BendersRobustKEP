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

/**@file   find_graphstructures.cpp
 * @brief  auxiliary methods to find substructures of a graph
 * @author Bart Smeulders
 * @author Danny Blom
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "graph.h"
#include "find_graphstructures.h"

/** check whether a node is contained in a path */
static
bool Is_already_in_path(
   const vector<int>&    path,                //!< vector encoding a path
   int                   v                    //!< node to be checked
   )
{
   for (int i = 0; i < int(path.size()); ++i)
   {
      if ( path[i] == v )
         return true;
   }

   return false;
}

/** computes all cycles up to a specific length in a graph */
SCIP_RETCODE findcycles(
    SCIP*                scip,               //!< SCIP data structure
    Cycles**             cycles,             //!< pointer to store cycles
    Graph*               G                   //!< pointer to graph
    )
{
   assert( scip != NULL );
   assert( G != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, cycles) );
   cout << "Enumerating all cycles" << endl;

   int maxcyclelen = 0;
   int cyclelengthlimit;
   int* adjacencylists;
   adjacencylists = G->adjacencylists;

   SCIP_CALL( SCIPgetIntParam(scip, "kidney/maxcyclelength", &cyclelengthlimit) );

   // This function identifies all cycles in the graph.
   // We first compute shortest paths from all nodes to nodes with a lower index.
   // Next we iteratively build paths from each node, only through nodes with higher index.
   // Paths are pruned if the shortest distance calculated previously shows the path can not
   // be completed to a cycle of low enough length.

   // Distance [i][j] is distance to i from j using only nodes with index k > i.
   // If the length is proven to be more than cyclelengthlimit, we end.
   std::vector<vector<int>> distances(G->npairs - 1);

   // We calculate distances for each graph copy. The distance is the distance back
   // to the vertex = graph copy. Exclude last pair.../
   for (int i = 0; i < G->npairs - 1; i++)
   {
      distances[i].resize(G->npairs, cyclelengthlimit);
      distances[i][i] = 0;

      // First calculate the distance "1" vertices.
      // To this end, iterate over all nodes with a higher label and check whether there is an arc back to i
      for (int j = i + 1; j < G->npairs; ++j)
      {
         for (int k = G->adjacencylistbegins[j]; k < G->adjacencylistbegins[j + 1]; ++k)
         {
            // stop: there shouldn't be arcs arriving at an undirected donor
            if ( adjacencylists[k] >= G->npairs )
            {
               SCIPerrorMessage("there shouldn't be arcs from pairs to undirected donors, but found %d -> %d",
                  j, adjacencylists[k]);
               SCIPfreeBlockMemory(scip, &cycles);
               return SCIP_ERROR;
            }

            // If there is an arc straight to i, distance is 1..
            if (G->adjacencylists[k] == i)
               distances[i][j] = 1;
         }
      }

      // We need to consider at most a number of steps equal to the cyclelength.
      // If the path is longer, the exact distance does not matter.
      for (int j = 2; j <= cyclelengthlimit; j++)
      {
         for (int k = i + 1; k < G->npairs; ++k)
         {
            for (int l = G->adjacencylistbegins[k]; l < G->adjacencylistbegins[k + 1]; ++l)
            {
               // stop: there shouldn't be arcs arriving at an undirected donor
               if ( adjacencylists[l] >= G->npairs )
               {
                  SCIPerrorMessage("there shouldn't be arcs from pairs to undirected donors, but found %d -> %d",
                     k, adjacencylists[l]);
                  SCIPfreeBlockMemory(scip, &cycles);
                  return SCIP_ERROR;
               }

               if (distances[i][adjacencylists[l]] == j - 1 && distances[i][k] > j)
                  distances[i][k] = j;
            }
         }
      }
   }

   int cntcycle = 0;
   int cntnodesincycle = 0;
   vector<int> nodelists(0);
   vector<int> nodelistsbegin(0);

   // We loop through all vertices as possible starting points of the cycles.
   for (int i = 0; i < G->npairs - 1; i++)
   {
      queue<vector<int>> paths; // A queue (to allow deleting elements at the front) including all current paths.
      // We intialize the queue by adding all arcs from the initial vertex (if a short enough return path exists).
      for (int j = G->adjacencylistbegins[i]; j < G->adjacencylistbegins[i + 1]; j++)
      {
         // stop: there shouldn't be arcs arriving at an undirected donor
         if ( adjacencylists[j] >= G->npairs )
         {
            SCIPerrorMessage("there shouldn't be arcs from pairs to undirected donors, but found %d -> %d",
               i, adjacencylists[j]);
            SCIPfreeBlockMemory(scip, &cycles);
            return SCIP_ERROR;
         }

         // Note that due to the distance calc, no arcs to lower numbered vertices will
         // be used in this initialization (since return distance = cyclelength + 1);
         if (distances[i][G->adjacencylists[j]] < cyclelengthlimit)
         {
            vector<int> init_vector;
            init_vector.push_back(i);
            init_vector.push_back(G->adjacencylists[j]);
            paths.push(init_vector);
         }
      }
      while ( ! paths.empty() )
      {
         vector<int> path_vector = paths.front(); // Get the first path out and delete it from the queue.
         paths.pop();

         int endvertex = path_vector[path_vector.size() - 1];

         // Go through the arcs starting from the current last vertex on the path.
         for (int j = G->adjacencylistbegins[endvertex]; j < G->adjacencylistbegins[endvertex + 1]; j++)
         {
            // stop: there shouldn't be arcs arriving at an undirected donor
            if ( adjacencylists[j] >= G->npairs )
            {
               SCIPerrorMessage("there shouldn't be arcs from pairs to undirected donors, but found %d -> %d",
                  endvertex, adjacencylists[j]);
               SCIPfreeBlockMemory(scip, &cycles);
               return SCIP_ERROR;
            }

            // If the path returns to the intial vertex, save the cycle.
            if (G->adjacencylists[j] == i)
            {
               cntcycle = cntcycle + 1;
               cntnodesincycle = cntnodesincycle + int(path_vector.size());
               nodelistsbegin.push_back(nodelists.size());
               for(int k = 0; k < int(path_vector.size()); k++)
               {
                  nodelists.push_back(path_vector[k]);
               }

            }
            else if ( ! Is_already_in_path(path_vector, G->adjacencylists[j])
               && distances[i][G->adjacencylists[j]] + int(path_vector.size()) - 1 < cyclelengthlimit )
            {
               vector<int> new_path = path_vector;
               new_path.push_back(G->adjacencylists[j]);
               paths.push(new_path);
            }
         }
      }
   }
   nodelistsbegin.push_back(nodelists.size());

   // collect cycles in arrays
   (*cycles)->ncycles = cntcycle;
   cout << "Number of cycles: " << cntcycle << endl;

   (*cycles)->nnodesincycles = cntnodesincycle;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*cycles)->nodelists), cntnodesincycle) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*cycles)->nodelistsbegin), cntcycle + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*cycles)->cycleweights), cntcycle) );

   for(int i = 0; i < cntnodesincycle; i++)
      (*cycles)->nodelists[i] = nodelists[i];
   cout << endl;
   for(int i = 0; i < cntcycle; i++)
   {
      (*cycles)->nodelistsbegin[i] = nodelistsbegin[i];
      (*cycles)->cycleweights[i] = nodelistsbegin[i+1] - nodelistsbegin[i];

      /* keep track of maximum cycle length */
      if ( (*cycles)->cycleweights[i] > maxcyclelen )
         maxcyclelen = (*cycles)->cycleweights[i];
   }
   (*cycles)->nodelistsbegin[cntcycle] = int(cntnodesincycle);
   (*cycles)->maxcyclelen = maxcyclelen;

   // for(int i=0; i < cntcycle; ++i)
   // {
   //    for(int j=nodelistsbegin[i]; j < nodelistsbegin[i+1]; ++j)
   //    {
   //       printf("%d ", (*cycles)->nodelists[j]);
   //    }
   //    printf("\n");
   // }

   return SCIP_OKAY;
}

/** computes all chains up to a specific length in a graph */
SCIP_RETCODE findchains(
    SCIP*                scip,               //!< SCIP data structure
    Chains**             chains,             //!< pointer to store chains
    Graph*               G                   //!< pointer to graph
    )
{
   assert( scip != NULL );
   assert( G != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, chains) );
   cout << "Enumerating all chains" << endl;

   SCIP_Bool issubchain;
   int maxchainlen = 0;
   int chainlen;
   int chainlengthlimit;

   SCIP_CALL( SCIPgetIntParam(scip, "kidney/maxchainlength", &chainlengthlimit) );
   cout << "Chain length maximum: " << chainlengthlimit << endl;

   // This function identifies all chains in the graph.
   // Next we iteratively build paths from each NDD.

   int cntchain = 0;
   int cntnodesinchain = 0;
   vector<int> nodelists(0);
   vector<int> nodelistsbegin(0);

   // We loop through all vertices as possible starting points of the chains.
   for (int i = G->npairs; i < G->nnodes; i++)
   {
      queue<vector<int>> paths; // A queue (to allow deleting elements at the front) including all current paths.

      // We intialize the queue by adding all arcs from the initial vertex (if a short enough return path exists).
      for (int j = G->adjacencylistbegins[i]; j < G->adjacencylistbegins[i + 1]; j++)
      {
         // stop: there shouldn't be arcs arriving at an undirected donor
         if ( G->adjacencylists[j] >= G->npairs )
         {
            SCIPerrorMessage("there shouldn't be arcs arriving at an undirected donor, but found %d -> %d",
               i, G->adjacencylists[j]);
            SCIPfreeBlockMemory(scip, &chains);
            return SCIP_ERROR;
         }

         vector<int> init_vector;
         init_vector.push_back(i);
         init_vector.push_back(G->adjacencylists[j]);

         // Add the chain to the list.
         cntchain = cntchain + 1;
         cntnodesinchain = cntnodesinchain + int(init_vector.size());
         nodelistsbegin.push_back(nodelists.size());
         for(int k = 0; k < int(init_vector.size()); k++)
            nodelists.push_back(init_vector[k]);
         paths.push(init_vector);
      }

      while ( ! paths.empty() )
      {
         vector<int> path_vector = paths.front(); // Get the first path out and delete it from the queue.
         paths.pop();
         int endvertex = path_vector[path_vector.size() - 1];

         // Go through the arcs starting from the current last vertex on the path.
         for (int j = G->adjacencylistbegins[endvertex]; j < G->adjacencylistbegins[endvertex + 1]; j++)
         {
            // stop: there shouldn't be arcs arriving at an undirected donor
            if ( G->adjacencylists[j] >= G->npairs )
            {
               SCIPerrorMessage("there shouldn't be arcs arriving at an undirected donor, but found %d -> %d",
                  endvertex, G->adjacencylists[j]);
               SCIPfreeBlockMemory(scip, &chains);
               return SCIP_ERROR;
            }

            // Only consider the next arc if it does no return to pair found earlier in the path.
            if (!Is_already_in_path(path_vector, G->adjacencylists[j]))
            {
               vector<int> new_path = path_vector;
               new_path.push_back(G->adjacencylists[j]);

               // Add the chain to the list.
               cntchain = cntchain + 1;
               cntnodesinchain = cntnodesinchain + int(new_path.size());
               nodelistsbegin.push_back(nodelists.size());
               for(int k = 0; k < int(new_path.size()); k++)
                  nodelists.push_back(new_path[k]);

               // If the maximum chainlength is not yet reached, add path back to the queue for possible later additions.
               if( int(new_path.size()) <= chainlengthlimit )
                  paths.push(new_path);
            }
         }
      }
   }
   nodelistsbegin.push_back(nodelists.size());

   // collect chains in arrays
   (*chains)->nchains = cntchain;
   cout << "Number of chains: " << cntchain << endl;

   (*chains)->nnodesinchains = cntnodesinchain;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*chains)->nodelists), cntnodesinchain) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*chains)->nodelistsbegin), cntchain + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*chains)->chainweights), cntchain) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*chains)->subchains), cntchain) );

   for (int i = 0; i < cntnodesinchain; i++)
      (*chains)->nodelists[i] = nodelists[i];
   cout << endl;

   for(int i = 0; i < cntchain; i++)
   {
      (*chains)->nodelistsbegin[i] = nodelistsbegin[i];
      (*chains)->chainweights[i] = nodelistsbegin[i+1] - nodelistsbegin[i] - 1; // -1 as we must count the number of transplants (NDD does not receive a transplant)

      /* keep track of maximum chain length */
      if ( (*chains)->chainweights[i] > maxchainlen )
         maxchainlen = (*chains)->chainweights[i];
   }
   (*chains)->nodelistsbegin[cntchain] = int(cntnodesinchain);
   (*chains)->maxchainlen = maxchainlen;

   /**< We make use of the fact that smaller chains have smaller indices based on our chain finding algorithm */
   for( int i = 0; i < cntchain; i++ )
   {
      (*chains)->subchains[i] = -1;
      chainlen = nodelistsbegin[i+1] - nodelistsbegin[i];

      // This chain has no proper subchain
      if( chainlen == 2 )
         continue;

      for( int j = 0; j < i; ++j )
      {
         if( chainlen != nodelistsbegin[j+1] - nodelistsbegin[j] + 1 )
            continue;

         issubchain = TRUE;
         for( int k = 0; k < chainlen - 1; ++k )
         {
            if( (*chains)->nodelists[nodelistsbegin[i]+k] != (*chains)->nodelists[nodelistsbegin[j]+k] )
            {
               issubchain = FALSE;
               break;
            }
         }
         if( issubchain )
         {
            (*chains)->subchains[i] = j;
            break;
         }
      }
   }
   return SCIP_OKAY;
}

/** computes all triplets in a graph */
SCIP_RETCODE findtriplets(
    SCIP*                scip,               //!< SCIP data structure
    Cycles*              cycles,             //!< array storing cycles
    Chains*              chains,             //!< array storing chains
    Triplets**           triplets,           //!< pointer to store triplets
    Graph*               G                   //!< pointer to graph
    )
{
   assert( scip != NULL );
   assert( cycles != NULL);
   assert( chains != NULL);
   assert( G != NULL );

   SCIP_Bool usetwothirdcliques;
   SCIP_Longint* adjmat;
   int lb_arcs;
   SCIP_CALL( SCIPallocBlockMemory(scip, triplets) );
   SCIP_CALL( SCIPgetBoolParam(scip, "kidney/usetwothirdcliques", &usetwothirdcliques));
   SCIP_CALL( SCIPgetIntParam(scip, "kidney/minarccounttriplet", &lb_arcs));

   if ( usetwothirdcliques )
   {
      cout << "Enumerating all triplets with at least " << lb_arcs << " arcs." << endl;
      int* adjacencylists = G->adjacencylists;
      int* adjacencylistbegins  =G->adjacencylistbegins;
      int nnodes = G->nnodes;
      int ncycles = cycles->ncycles;
      int nchains = chains->nchains;

      // This function identifies all triplets with a sufficient number of arcs in the graph.

      // We first build an adjacency matrix of the directed graph (V x V matrix M with M[u,v] = 1 if (u,v) in arcset)
      SCIP_CALL( SCIPallocClearBufferArray(scip, &adjmat, nnodes*nnodes) );
      int u, v, w, i, t, c;
      int narcs;
      int tripletcount = 0;

      for (u = 0; u < nnodes; u++)
      {
         for (v = adjacencylistbegins[u]; v < adjacencylistbegins[u+1]; ++v)
            adjmat[u*nnodes+adjacencylists[v]] = 1;
      }

      // Next we iteratively search for each triple the number of arcs in their respective subgraphs.
      std::vector<int> tripletlist(0); // triplets is an array with all the nodes occurring in a triplet
      int npairs = G->npairs;
      for (u = 0; u < npairs; ++u)
      {
         for (v = u + 1; v < npairs; ++v)
         {
            for (w = v + 1; w < nnodes; ++w)
            {
               narcs = adjmat[u*nnodes+v] + adjmat[u*nnodes+w] + adjmat[v*nnodes+u] + adjmat[v*nnodes+w] + adjmat[w*nnodes+u] + adjmat[w*nnodes+v];
               if (narcs >= lb_arcs )
               {
                  tripletlist.push_back(u);
                  tripletlist.push_back(v);
                  tripletlist.push_back(w);
                  ++tripletcount;
               }
            }
         }
      }

      SCIPfreeBufferArray(scip, &adjmat);

      // collect triplets in array
      (*triplets)->ntriplets = tripletcount;
      if (tripletcount == 0)
      {
         (*triplets)->nodelists = NULL;
         (*triplets)->cycle2tripletsbegin = NULL;
         (*triplets)->cycle2triplets = NULL;
         (*triplets)->chain2tripletsbegin = NULL;
         (*triplets)->chain2triplets = NULL;
         (*triplets)->nindicescycle2triplet = 0;
         (*triplets)->nindiceschain2triplet = 0;
         return SCIP_OKAY;
      }

      // Whenever the tripletcount is larger than zero, we want to consider allocation of memory for the nodelists
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*triplets)->nodelists), 3*tripletcount) );
      for(i = 0; i < 3*tripletcount; i++)
         (*triplets)->nodelists[i] = tripletlist[i];

      // Determine for each cycle which triplets have at least two nodes in common
      std::vector<int> cycle2triplets(0);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*triplets)->cycle2tripletsbegin), ncycles + 1) );
      int cyclebegin = 0;
      int intersection; // Keeps track of intersection of cycle and triplet
      int node;
      for (c = 0; c < ncycles; ++c)
      {
         (*triplets)->cycle2tripletsbegin[c] = cyclebegin;
         for( t = 0; t < tripletcount; ++t)
         {
            intersection = 0;
            for (i = cycles->nodelistsbegin[c]; i < cycles->nodelistsbegin[c + 1]; ++i)
            {
               node = cycles->nodelists[i];
               if (node == (*triplets)->nodelists[3*t] || node == (*triplets)->nodelists[3*t+1]
                  || node == (*triplets)->nodelists[3*t+2] )
                  ++intersection;
            }
            if ( intersection >= 2)
            {
               cycle2triplets.push_back(t);
               ++cyclebegin;
            }
         }
      }
      (*triplets)->cycle2tripletsbegin[ncycles] = cyclebegin;

      if( cyclebegin > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*triplets)->cycle2triplets), cyclebegin) );
         for(i = 0; i < cyclebegin; ++i)
            (*triplets)->cycle2triplets[i] = cycle2triplets[i];
      }
      else
      {
         (*triplets)->cycle2triplets = NULL;
      }
      (*triplets)->nindicescycle2triplet = cyclebegin;

      // Determine for each chain which triplets have at least two nodes in common
      std::vector<int> chain2triplets(0);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*triplets)->chain2tripletsbegin), nchains + 1) );
      int chainbegin = 0;
      for (c = 0; c < nchains; ++c)
      {
         (*triplets)->chain2tripletsbegin[c] = chainbegin;
         for( t = 0; t < tripletcount; ++t)
         {
            intersection = 0;
            for (i = chains->nodelistsbegin[c]; i < chains->nodelistsbegin[c + 1]; ++i)
            {
               node = chains->nodelists[i];
               if (node == (*triplets)->nodelists[3*t] || node == (*triplets)->nodelists[3*t+1]
                  || node == (*triplets)->nodelists[3*t+2] )
                  ++intersection;
            }
            if ( intersection >= 2)
            {
               chain2triplets.push_back(t);
               chainbegin++;
            }
         }
      }
      (*triplets)->chain2tripletsbegin[nchains] = chainbegin;

      if( chainbegin > 0)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*triplets)->chain2triplets), chainbegin) );
         for(i = 0; i < chainbegin; ++i)
            (*triplets)->chain2triplets[i] = chain2triplets[i];
      }
      else
         (*triplets)->chain2triplets = NULL;
      (*triplets)->nindiceschain2triplet = chainbegin;
   }
   else
   {
      // Parameter to use "two out of three" cliques is set to false, so no triplets need to be considered
      (*triplets)->ntriplets = 0;
      (*triplets)->nodelists = NULL;
      (*triplets)->cycle2triplets = NULL;
      (*triplets)->cycle2tripletsbegin = NULL;
      (*triplets)->nindicescycle2triplet = 0;
      (*triplets)->chain2triplets = NULL;
      (*triplets)->chain2tripletsbegin = NULL;
      (*triplets)->nindiceschain2triplet = 0;
   }
   return SCIP_OKAY;
}

/** generate positioned arcs in graph */
SCIP_RETCODE generatePositionedArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                G,                  /**< reference to graph */
   PositionedArcs**      posarcs             /**< Pointer to store positioned arcs */
   )
{
   assert( scip != NULL );
   assert( G != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, posarcs) );

   int cntposarcs = 0;
   int chainlengthlimit;
   int i;
   int k;
   int c;
   int s;
   int v;
   int curr_dist = -1;
   vector<int> nodelists(0);
   vector<int> nodelistsbegin(0);

   /* To reduce the size of the problem, we only generate positioned arcs that CAN exist,
    * e.g. a pair that has minimum distance l from all NDDs can only have outgoing arcs with position l+1, ...
    * We use a breadth-first search approach to find out these distances.
    */
   SCIP_CALL( SCIPgetIntParam(scip, "kidney/maxchainlength", &chainlengthlimit) );

   // Mark all the vertices as not visited
   int* distances;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &distances, G->nnodes) );
   for(i = 0; i < G->nnodes; i++)
       distances[i] = chainlengthlimit;

   // Create a queue for BFS
   queue<int> queue;

   // The non-directed donors are at the start of a chain, so set distance from NDD to 0
   for (v = G->npairs; v < G->nnodes; ++v)
   {
      queue.push(v);
      distances[v] = 0;
   }

   while( ! queue.empty() )
   {
      // Dequeue a vertex from queue and print it
      s = queue.front();
      queue.pop();
      curr_dist = distances[s];

      // Get all adjacent vertices of the dequeued
      // vertex s. If a adjacent has not been visited,
      // then mark it visited and enqueue it
      for (i = G->adjacencylistbegins[s]; i < G->adjacencylistbegins[s + 1]; ++i)
      {
         v = G->adjacencylists[i];

         // If not visited
         if ( curr_dist < chainlengthlimit - 1 && distances[v] == chainlengthlimit )
         {
            distances[v] = curr_dist+1;
            queue.push(v);
         }
      }
   }

   for (k = 1; k <= chainlengthlimit; ++k)
   {
      if ( k > 1 )
      {
         nodelistsbegin.push_back(2*cntposarcs);

         // Only non-NDD arcs have value 2 to L
         for (i = 0; i < G->npairs; ++i)
         {
            // Only include arc if BFS distances allow for it
            if (distances[i] <= k - 1)
            {
               for (c = G->adjacencylistbegins[i]; c < G->adjacencylistbegins[i+1]; ++c)
               {
                  cntposarcs++;
                  nodelists.push_back(i);
                  nodelists.push_back(G->adjacencylists[c]);
               }
            }
         }
      }
      else
      {
         nodelistsbegin.push_back(0);

         // Outgoing arcs from NDDs have index 1
         for (i = G->npairs; i < G->nnodes; ++i)
         {
            for (c = G->adjacencylistbegins[i]; c < G->adjacencylistbegins[i+1]; ++c)
            {
               cntposarcs++;
               nodelists.push_back(i);
               nodelists.push_back(G->adjacencylists[c]);
            }
         }
      }
   }
   nodelistsbegin.push_back(2*cntposarcs);

   (*posarcs)->nposarcs = cntposarcs;

   SCIPfreeBlockMemoryArrayNull(scip, &distances, G->nnodes);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*posarcs)->nodelists), 2*cntposarcs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*posarcs)->positionbegins), chainlengthlimit+1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*posarcs)->arcweights), cntposarcs) );

   /* collect position indexed arcs in array*/
   for(i = 0; i < 2*cntposarcs; i++)
      (*posarcs)->nodelists[i] = nodelists[i];

   for(i = 0; i <= chainlengthlimit; i++)
      (*posarcs)->positionbegins[i] = nodelistsbegin[i];

   for (i = 0; i < cntposarcs; i++)
      (*posarcs)->arcweights[i] = 1.0;

   (*posarcs)->npositions = chainlengthlimit;

   return SCIP_OKAY;
}

/** generate map from nodes to the cycles in a graph */
SCIP_RETCODE generateNode2Cycles(
   SCIP*                 scip,               /**< SCIP instance */
   Graph*                G,                  /**< underlying graph */
   Cycles*               cycles              /**< underlying cycles structure */
   )
{
   assert( scip != NULL );
   assert( G != NULL );
   assert( cycles != NULL );

   int nnodes = G->nnodes;
   int ncycles = cycles->ncycles;
   int* nodelists = cycles->nodelists;
   int* nodelistsbegin = cycles->nodelistsbegin;

   /* for each node, store the indices of cycles containing this node */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(G->node2cycles), cycles->nnodesincycles) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(G->node2cyclesbegin), nnodes + 1) );

   int newbegin = 0;
   for (int i = 0; i < nnodes; ++i)
   {
      G->node2cyclesbegin[i] = newbegin;
      for (int c = 0; c < ncycles; ++c)
      {
         for (int j = nodelistsbegin[c]; j < nodelistsbegin[c + 1]; ++j)
         {
            if ( nodelists[j] == i )
            {
               G->node2cycles[newbegin++] = c;
               break;
            }
         }
      }
   }
   assert( cycles->nnodesincycles == newbegin );
   G->node2cyclesbegin[nnodes] = newbegin;

   return SCIP_OKAY;
}

/** generate map from nodes to the chains in a graph */
SCIP_RETCODE generateNode2Chains(
   SCIP*                 scip,               /**< SCIP instance */
   Graph*                G,                  /**< underlying graph */
   Chains*               chains              /**< underlying chains structure */
   )
{
   assert( scip != NULL );
   assert( G != NULL );
   assert( chains != NULL );

   int nnodes = G->nnodes;
   int nchains = chains->nchains;
   int* nodelists = chains->nodelists;
   int* nodelistsbegin = chains->nodelistsbegin;

   /* for each node, store the indices of chains containing this node */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(G->node2chains), chains->nnodesinchains) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(G->node2chainsbegin), nnodes + 1) );

   int newbegin = 0;
   for (int i = 0; i < nnodes; ++i)
   {
      G->node2chainsbegin[i] = newbegin;
      for (int c = 0; c < nchains; ++c)
      {
         for (int j = nodelistsbegin[c]; j < nodelistsbegin[c + 1]; ++j)
         {
            if ( nodelists[j] == i )
            {
               G->node2chains[newbegin++] = c;
               break;
            }
         }
      }
   }
   assert( chains->nnodesinchains == newbegin );
   G->node2chainsbegin[nnodes] = newbegin;

   return SCIP_OKAY;
}
