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

/**@file   probdata_glorie.h
 * @brief  Problem data for Glorie et al. approach
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_GLORIE_H
#define PROBDATA_GLORIE_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPglorieEtAldataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   );

/** returns cycle variables */
SCIP_VAR** getCyclevarsGlorie(
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   );

/** returns chain variables */
SCIP_VAR** getChainvarsGlorie(
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
