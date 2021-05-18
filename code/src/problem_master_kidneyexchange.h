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

/**@file   problem_master_kidneyexchange.h
 * @brief  methods to create the master problem
 * @author Christopher Hojny
 *
 */

#ifndef PROBLEM_MASTER_KIDNEYEXCHANGE_H
#define PROBLEM_MASTER_KIDNEYEXCHANGE_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/* creates initial model of the master kidney exchange problem */
extern
SCIP_RETCODE SCIPcreateMasterModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/* free master kidney exchange problem data */
extern
SCIP_RETCODE SCIPfreeMasterModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains,             /**< chain structures of graph */
   Triplets*             triplets            /**< triplet structures of graph */
   );

#ifdef __cplusplus
}
#endif

#endif
