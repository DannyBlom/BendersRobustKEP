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

/**@file   read_graph.cpp
 * @brief  methods for reading the graph of our problem
 * @author Christopher Hojny
 */

#include <iostream>
#include <fstream>
#include <sstream>

#include "graph.h"
#include "read_graph.h"

using namespace std;

/** reads a directed graph from a file */
SCIP_RETCODE readGraph(
   SCIP*                 scip,               //<! SCIP data structure
   string                filename,           //<! name of file encoding graph
   Graph**               graph               //<! pointer to store graph
   )
{
   assert( scip != NULL );
   assert( graph != NULL );

   cout << "Reading graph from file " << filename << endl;
   ifstream input(filename);
   if (!input)
   {
      cout << "No input file found." << endl;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, graph) );

   // read number of nodes in graph
   input.ignore(50, '=');
   input >> (*graph)->npairs;
   input.ignore(50, '=');
   input >> (*graph)->ndonors;

   (*graph)->nnodes = (*graph)->npairs + (*graph)->ndonors;

   // read number of arc in graph
   input.ignore(50, '=');
   input >> (*graph)->narcs;

   // store node list and collect information on nodes
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*graph)->nodelist), (*graph)->nnodes) );
   int cnt = 0;
   for (cnt = 0; cnt < (*graph)->npairs; ++cnt)
   {
      Node* node;
      SCIP_CALL( SCIPallocBlockMemory(scip, &node) );

      node->id = cnt;
      node->ispair = TRUE;

      // call it twice, because we ignore the labeling given in the file
      input >> node->failprob;
      input >> node->failprob;

      node->ub = 0.0;

      ((*graph)->nodelist)[cnt] = node;
   }
   for (int i = 0; i < (*graph)->ndonors; i++)
   {
      Node* node;
      SCIP_CALL( SCIPallocBlockMemory(scip, &node) );

      node->id = cnt;
      node->ispair = FALSE;

      // call it twice, because we ignore the labeling given in the file
      input >> node->failprob;
      input >> node->failprob;

      node->ub = 0.0;

      ((*graph)->nodelist)[cnt++] = node;
   }

   // collect arcs
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*graph)->adjacencylists), (*graph)->narcs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*graph)->adjacencylistbegins), (*graph)->nnodes + 1) );

   cnt = 0;
   string str0;
   getline(input, str0);

   int u;
   int v;
   int oldu = -1;

   while ( TRUE )
   {
      string str;
      stringstream stream;
      getline(input, str);

      if (str[0] == '(')
      {
         stream.str(str);
         stream.ignore(1, '(');
         stream >> u;
         stream.ignore(1, ',');
         stream >> v;

         // we have found a new start node, update begins position also for intermediate nodes
         // without leaving arcs
         if ( oldu != u )
         {
            for (int i = oldu + 1; i <= u; ++i)
               (*graph)->adjacencylistbegins[i] = cnt;
         }

         (*graph)->adjacencylists[cnt++] = v;
         oldu = u;
      }
      else
         break;
   }
   assert( cnt == (*graph)->narcs );

   for (int i = oldu + 1; i < (*graph)->nnodes; ++i)
      (*graph)->adjacencylistbegins[i] = cnt;
   (*graph)->adjacencylistbegins[(*graph)->nnodes] = (*graph)->narcs;

   return SCIP_OKAY;
}

/** strengthens bounds on node variables to length of maximum cycle/chain */
SCIP_RETCODE strengthenNodeBounds(
   SCIP*                 scip,               //<! SCIP instance
   Graph*                graph,              //<! underlying graph
   Cycles*               cycles,             //<! underlying cycles structure
   Chains*               chains              //<! underlying chain structure
   )
{
   Node** graphnodelist;
   int* cnodelists;
   int* cnodelistsbegin;
   SCIP_Real ub;
   int c;
   int v;
   int ncycles;
   int nchains;
   int length;

   assert( scip != NULL );
   assert( graph != NULL );
   assert( graph->nodelist != NULL );
   assert( cycles != NULL );
   assert( cycles->nodelists != NULL );
   assert( cycles->nodelistsbegin != NULL );
   assert( chains != NULL );
   assert( chains->nodelists != NULL );
   assert( chains->nodelistsbegin != NULL );

   ncycles = cycles->ncycles;
   nchains = chains->nchains;

   graphnodelist = graph->nodelist;

   cnodelists = cycles->nodelists;
   cnodelistsbegin = cycles->nodelistsbegin;

   for (c = 0; c < ncycles; ++c)
   {
      length = cnodelistsbegin[c + 1] - cnodelistsbegin[c];
      for (v = cnodelistsbegin[c]; v < cnodelistsbegin[c + 1]; ++v)
      {
         ub = graphnodelist[cnodelists[v]]->ub;
         if ( SCIPisGT(scip, length, ub) )
            graphnodelist[cnodelists[v]]->ub = length;
      }
   }

   cnodelists = chains->nodelists;
   cnodelistsbegin = chains->nodelistsbegin;

   for (c = 0; c < nchains; ++c)
   {
      length = cnodelistsbegin[c + 1] - cnodelistsbegin[c];
      for (v = cnodelistsbegin[c]; v < cnodelistsbegin[c + 1]; ++v)
      {
         ub = graphnodelist[cnodelists[v]]->ub;
         if ( SCIPisGT(scip, length, ub) )
            graphnodelist[cnodelists[v]]->ub = length;
      }
   }

   return SCIP_OKAY;
}
