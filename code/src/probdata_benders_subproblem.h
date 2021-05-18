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

/**@file   probdata_kidney.h
 * @brief  Problem data for kidney exchange problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_BENDERS_SUBPROBLEM_H
#define PROBDATA_BENDERS_SUBPROBLEM_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPKEPdataCreate(
   SCIP*                 kepscip,            /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
  );

/** returns the number of nodes in the KEP instance */
int SCIPKEPdataGetNumNodes(
  SCIP_PROBDATA*        kepdata              /**< problem data */
  );

/** returns array of cycle variables */
SCIP_VAR** SCIPKEPdataGetCyclevars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of cycles */
int SCIPKEPdataGetNCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of chain variables */
SCIP_VAR** SCIPKEPdataGetChainvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of chains */
int SCIPKEPdataGetNChains(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of node constraints */
SCIP_CONS** SCIPKEPdataGetNodeconss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


#ifdef __cplusplus
}
#endif

#endif
