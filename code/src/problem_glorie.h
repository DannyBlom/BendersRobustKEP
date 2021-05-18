/**@file   problem_glorie.h
 * @brief  Problem data for Glorie et al. approach to the robust kidney exchange problem
 * @author Christopher Hojny
 * @author Bart Smeulders
 */

#ifndef PROBLEM_GLORIE_H
#define PROBLEM_GLORIE_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateGlorieEtAlModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   Chains*               chains              /**< pointer to chain structures of graph */
   );

/** free kidney exchange problem data */
SCIP_RETCODE SCIPfreeGlorieEtAlModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains              /**< chain structures of graph */
   );

/** extends the Glorie attack greedily, and also computes a lower bound based on this attack. */
SCIP_Real GlorieExtendAttack(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   Chains*               chains,             /**< chain structures of graph */
   int*                  usedcycleschains,   /**< array containing indices of cycles and chains used in master solution */
   int                   nusedcycleschains,  /**< number of cycles and chains used in master solution */
   int*                  zerofixings,        /**< array of 0-fixings (identified by variable indices) */
   int*                  onefixings,         /**< array of 1-fixings (identified by variable indices) */
   int                   nzerofixings,       /**< number of 0-fixings */
   int                   nonefixings,        /**< number of 1-fixings */
   int                   attackbound,        /**< upper bound on number of attacks */
   int*                  attackpattern       /**< the current attack pattern */
);

#ifdef __cplusplus
}
#endif

#endif
