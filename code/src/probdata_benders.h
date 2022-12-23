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

/**@file   probdata_benders.h
 * @brief  Problem data for Benders type approach
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_BENDERS_H
#define PROBDATA_BENDERS_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPbendersdataCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   const char*           probname,           /**< problem name */
   SCIP_SOL*             sol,                /**< SCIP solution to masterscip */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** returns whether a given cycle intersects a target cycle / chain in initial solution */
SCIP_Bool SCIPcycleIntersectsTarget(
   Cycles*               cycles,
   Chains*               chains,
   int                   cycle_idx,
   int                   target_idx
   );

/**< returns whether chain2 is a superchain of chain1 */
SCIP_Bool SCIPchainIsSuperchain(
   Chains*              chains,              /**< chains data structure */
   int                  chain1,              /**< index of chain */
   int                  chain2               /**< index of possible superchain */
   );

/** returns whether a given chain intersects a target cycle / chain in initial solution */
SCIP_Bool SCIPchainIntersectsTarget(
   Cycles*               cycles,
   Chains*               chains,
   int                   chain_idx,
   int                   target_idx
   );

/** adds given variable to the problem data */
SCIP_RETCODE SCIPbendersdataAddSolCons(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   int                   nvars               /**< number of variables */
   );

/** returns array of cycle lower bound constraints */
SCIP_CONS** SCIPbendersdataGetCyclelbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of cycle upper bound constraints */
SCIP_CONS** SCIPbendersdataGetCycleubconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of cycle variables */
SCIP_VAR** SCIPbendersdataGetCyclevars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of cycles */
int SCIPbendersdataGetNCycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of chain lower bound constraints */
SCIP_CONS** SCIPbendersdataGetChainlbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of chain upper bound constraints */
SCIP_CONS** SCIPbendersdataGetChainubconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of chain variables */
SCIP_VAR** SCIPbendersdataGetChainvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of chains */
int SCIPbendersdataGetNChains(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns uvars */
SCIP_VAR** SCIPbendersdataGetUVars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns the budget for the adversary */
int SCIPbendersdataGetAdversaryBound(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns the number of nodes in the instance */
int SCIPbendersdataGetNumNodes(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns objvar */
SCIP_VAR* SCIPbendersdataGetObjvar(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
