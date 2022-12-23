/**@file   probdata_master_kidney_picef.h
 * @brief  Problem data for master kidney exchange problem with position indexed cycle edge formulation
 * @author Danny Blom
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_MASTER_KIDNEY_PICEF_H
#define PROBDATA_MASTER_KIDNEY_PICEF_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPmasterPICEFProbdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** updates master problem by a new attack pattern */
SCIP_RETCODE SCIPupdateMasterPICEFProblem(
   SCIP*                 scip,               /**< SCIP data structure  */
   int*                  attackpattern,      /**< array containing nodes being attacked */
   int                   nattacks            /**< number of attacks encoded in attackpattern */
   );

/** get objvar */
SCIP_VAR* masterPICEFproblemGetObjvar(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get cycles */
Cycles* masterPICEFProblemGetCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get xvarinit array */
SCIP_VAR** masterPICEFProblemGetXCyclevarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get posarcs */
PositionedArcs* masterPICEFProblemGetPosarcs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get arcvars array */
SCIP_VAR** masterPICEFProblemGetArcvarinit(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get dummyyvars array */
SCIP_VAR** masterPICEFProblemGetDummyYVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get attack scenarios array */
int* masterPICEFProblemGetScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** get number of attack scenarios */
int masterPICEFProblemGetNScenarios(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
