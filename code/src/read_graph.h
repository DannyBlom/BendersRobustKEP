// -*- C++ -*-
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


#ifndef READ_GRAPH_H
#define READ_GRAPH_H

#include <string>
#include "graph.h"

using namespace std;

/** reads a directed graph from a file */
extern
SCIP_RETCODE readGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   string                filename,           /**< name of file encoding graph */
   Graph**               graph               /**< pointer to store graph */
   );

/** strengthens bounds on node variables to length of maximum cycle/chain */
extern
SCIP_RETCODE strengthenNodeBounds(
   SCIP*                 scip,               /**< SCIP instance */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< underlying cycles structure */
   Chains*               chains              /**< underlying chain structure */
   );

#endif
