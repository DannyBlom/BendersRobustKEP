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

/**@file   problem_kidneyexchange.h
 * @brief  Problem data for kidney exchange problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef PROBLEM_KIDNEYEXCHANGE_H
#define PROBLEM_KIDNEYEXCHANGE_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/* creates initial model of the kidney exchange problem */
extern
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   Triplets*             triplets,           /**< pointer to triplet structures of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Bool             userankseparation,  /**< whether rank separation shall be used */
   int                   maxrankrhs          /**< maximum rhs of rank inequalities */
   );

/* free kidney exchange problem data */
extern
SCIP_RETCODE SCIPfreeModel(
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
