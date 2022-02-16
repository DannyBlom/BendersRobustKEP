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

/**@file   graph.h
 * @brief Declaration of graph data structures
 * @author Christopher Hojny
 * @author Danny Blom
 *
 * This file contains the definitions of graph related data structures
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef GRAPH_H
#define GRAPH_H

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of a node */
typedef struct Node
{
   int                   id;                 /**< node ID */
   SCIP_Real             failprob;           /**< failure probability */
   SCIP_Bool             ispair;             /**< whether the node encodes a pair */
   SCIP_Real             ub;                 /**< upper bound for T- and W-variables */
} Node;

/** data of directed graph */
typedef struct Graph
{
   int                   nnodes;             /**< number of nodes in graph */
   int                   npairs;             /**< number of patient/donor pairs encoded in graph */
   int                   ndonors;            /**< number of donors encoded in graph */
   Node**                nodelist;           /**< array containing node information */

   int                   narcs;              /**< number of arcs in graph */
   int*                  adjacencylists;     /**< array containing adjacency (outgoing arcs) list of all nodes */
   int*                  adjacencylistbegins;/**< array containing start positions of adjacencylists */

   int*                  node2cycles;        /**< array containing for each node the indices of cycles containing the node */
   int*                  node2cyclesbegin;   /**< array to encode where the sublist of cycles for a node begins in node2cycles */
   int*                  node2chains;        /**< array containing for each node the indices of chains containing the node */
   int*                  node2chainsbegin;   /**< array to encode where the sublist of chains for a node begins in node2chains */

} Graph;

typedef struct Cycles
{
   int                  ncycles;            /**< Number of cycles in the problem */
   int                  nnodesincycles;     /**< Sum of node appearences over all cycles. */
   int*                 nodelists;          /**< Array containing the list of nodes in all cycles.*/
   int*                 nodelistsbegin;     /**< Array containing start positions of the nodelists per cycle.*/
   int*                 cycleweights;       /**< Array containing the weights of each cycle.*/
   int                  maxcyclelen;        /**< maximum length of a cycle */
} Cycles;

typedef struct Chains
{
   int                  nchains;            /**< Number of chains in the problem */
   int                  nnodesinchains;     /**< Sum of node appearences over all chains. */
   int*                 nodelists;          /**< Array containing the list of nodes in all chains.*/
   int*                 nodelistsbegin;     /**< Array containing start positions of the nodelists per chain.*/
   int*                 chainweights;       /**< Array containing the weights of each chain.*/
   int*                 subchains;          /**< Array of length nchains indicating index of largest subchain (-1 if chain has 1 arc) */
   int                  maxchainlen;        /**< maximum length of a chain */
} Chains;

typedef struct Triplets
{
   int                  ntriplets;          /**< Number of triplets (satisfying some condition of number of arcs in the induced subgraph) */
   int*                 nodelists;          /**< Array containing the list of nodes in all triplets.*/
   int*                 cycle2triplets;     /**< Array containing for each cycle the indices of triplets with at least 2 intersecting nodes */
   int*                 cycle2tripletsbegin;/**< array to encode where the sublist of triplets for a node begins in cycle2triplets */
   int                  nindicescycle2triplet;/**< Number of triplet indices in cycle2triplets */
   int*                 chain2triplets;     /**< Array containing for each chain the indices of triplets with at least 2 intersecting nodes */
   int*                 chain2tripletsbegin;/**< array to encode where the sublist of triplets for a node begins in chain2triplets */
   int                  nindiceschain2triplet;/**< Number of triplet indices in chain2triplets */
} Triplets;

typedef struct PositionedArcs
{
   int                  nposarcs;            /**< Number of position arcs in the graph */
   int*                 nodelists;           /**< List of nodes of arcs that are position indices */
   int*                 positionbegins;      /**< First index in nodelists where an arc with position k starts, for each k */
   int                  npositions;
   int*                 arcweights;          /**< Array containing the weights of each position indexed arc */
} PositionedArcs;

#ifdef __cplusplus
}
#endif

#endif
