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

/**@file   solveMasterProblem.h
 * @brief  routines to solve the master problem
 * @author Christopher Hojny
 *
 * This file handles the master problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SOLVEMASTERPROBLEM_H
#define SOLVEMASTERPROBLEM_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** solves the master problem */
SCIP_RETCODE solveMasterProblem(
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   Triplets*             triplets,           /**< pointer to triple structures of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Real             timelimit,          /**< time limit to solve the problem */
   SCIP_Real             initialtime,        /**< time spent in initialization */
   const char*           settings,           /**< possible name of setting file */
   SCIP_Bool             verbose             /**< whether we print SCIP's logs */
   );

SCIP_RETCODE SCIPgetAttackPattern(
   SCIP*                 scip,               /**< SCIP instance of 2/3 stage */
   SCIP_SOL*             sol,                /**< solution of 2/3 stage */
   int*                  attackpattern,      /**< array to store attack pattern */
   int*                  nattacks,           /**< pointer to store number of attacks */
   int                   nnodes,             /**< number of nodes in underlying graph */
   int                   method              /**< Which method is used for stage 2/3? (1:Benders, 2:B&P) */
   );

#ifdef __cplusplus
}
#endif

#endif
