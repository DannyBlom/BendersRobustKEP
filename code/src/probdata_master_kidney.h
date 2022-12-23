/**@file   probdata_master_kidney.h
 * @brief  Problem data for master kidney exchange problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_MASTER_KIDNEY_H
#define PROBDATA_MASTER_KIDNEY_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPmasterProbdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** updates master problem by a new attack pattern */
SCIP_RETCODE SCIPupdateMasterProblem(
   SCIP*                 scip,               /**< SCIP data structure  */
   int*                  attackpattern,      /**< array containing nodes being attacked */
   int                   nattacks            /**< number of attacks encoded in attackpattern */
   );

/** returns whether a given cycle is attacked */
SCIP_Bool SCIPcycleIsAttacked(
   int*                  attackpattern,      /**< array of attacked vertices */
   int                   nattacks,           /**< number of attacked vertices */
   Cycles*               cycles,             /**< data structure for the cycles in the graph of suitable length */
   int                   cycleidx            /**< index of cycle in cycles data structure */
   );

/** returns whether a given chain is attacked */
SCIP_Bool SCIPchainIsAttacked(
   int*                  attackpattern,      /**< array of attacked vertices */
   int                   nattacks,           /**< number of attacked vertices */
   Chains*               chains,             /**< data structure for the chains in the graph of suitable length */
   int                   chainidx,           /**< index of chain in chains data structure */
   int                   vertex              /**< Consider only chain up until vertex (if NULL, consider entire chain) */
   );

/** get objvar */
SCIP_VAR* masterProblemGetObjvar(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get xvarinit array */
SCIP_VAR** masterProblemGetXvarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get dummyyvars array */
SCIP_VAR** masterProblemGetDummyYVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get attack scenarios array */
int* masterProblemGetScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get number of attack scenarios */
int masterProblemGetNScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get cycles of instance */
Cycles* masterProblemGetCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get chains of instance */
Chains* masterProblemGetChains(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
