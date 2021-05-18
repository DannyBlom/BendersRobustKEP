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

/**@file   find_graphstructures.h
 * @brief  auxiliary methods to find substructures of a graph
 * @author Bart Smeulders
 * @author Danny Blom
 */

#ifndef FIND_GRAPHSTRUCTURES_H
#define FIND_GRAPHSTRUCTURES_H

#include "graph.h"

using namespace std;

/** computes all cycles up to a specific length in a graph */
SCIP_RETCODE findcycles(
   SCIP*                 scip,               /**< SCIP data structure */
   Cycles**              cycles,             /**< pointer to store cycles */
   Graph*                G                   /**< Reference to graph (constant, so graph can not be changed) */
   );

/** computes all chains up to a specific length in a graph */
SCIP_RETCODE findchains(
   SCIP*                 scip,               /**< SCIP data structure */
   Chains**              cycles,             /**< pointer to store cycles */
   Graph*                G                   /**< Reference to graph (constant, so graph can not be changed) */
   );

/** generate positioned arcs in graph */
SCIP_RETCODE generatePositionedArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                G,                  /**< reference to graph */
   PositionedArcs**      posarcs             /**< Pointer to store positioned arcs */
   );

/** computes all triplets in a graph */
SCIP_RETCODE findtriplets(
   SCIP*                 scip,               /**< SCIP data structure */
   Cycles*               cycles,             /**< array storing cycles */
   Chains*               chains,             /**< array storing chains */
   Triplets**            triplets,           /**< pointer to store relevant triplets */
   Graph*                G                   /**< Reference to graph (constant, so graph can not be changed) */
   );

/** generate map from nodes to the cycles in a graph */
SCIP_RETCODE generateNode2Cycles(
   SCIP*                 scip,               /**< SCIP instance */
   Graph*                G,                  /**< underlying graph */
   Cycles*               cycles              /**< underlying cycles structure */
   );

/** generate map from nodes to the chains in a graph */
SCIP_RETCODE generateNode2Chains(
   SCIP*                 scip,               /**< SCIP instance */
   Graph*                G,                  /**< underlying graph */
   Chains*               chains              /**< underlying chains structure */
   );

#endif
