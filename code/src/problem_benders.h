/**@file   problem_benders.h
 * @brief  Problem data for Benders type approach to the robust kidney exchange problem
 * @author Danny Blom
 */

#ifndef PROBLEM_BENDERS_H
#define PROBLEM_BENDERS_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
   extern "C" {
#endif

/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateBendersModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   SCIP_SOL*             sol,                /**< SCIP solution to master problem */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains,             /**< pointer to chain structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
);

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains              /**< chain structures of graph */
   );

/* create kidney exchange problem data */
SCIP_RETCODE SCIPcreateBendersSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   );

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains              /**< chain structures of graph */
   );


#ifdef __cplusplus
}
#endif

#endif
