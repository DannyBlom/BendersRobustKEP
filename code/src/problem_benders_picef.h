/**@file   problem_benders_picef.h
 * @brief  Problem data for Benders type approach to the robust kidney exchange problem based on PICEF formulation
 * @author Danny Blom
 */

#ifndef PROBLEM_BENDERS_PICEF_H
#define PROBLEM_BENDERS_PICEF_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateBendersPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure of master problem */
   SCIP_SOL*             mastersol,          /**< SCIP solution structure of master solution */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to arc structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   PositionedArcs*       posarcs             /**< arc structures of graph */
   );

/** create kidney exchange problem data */
SCIP_RETCODE SCIPcreateBendersPICEFSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs             /**< pointer to arc structures of graph */
);

/** free kidney exchange problem data */
SCIP_RETCODE SCIPfreeBendersPICEFSubModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   PositionedArcs*       posarcs             /**< arc structures of graph */
   );

#ifdef __cplusplus
}
#endif

#endif
