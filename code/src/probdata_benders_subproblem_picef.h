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

/**@file  probdata_benders_subproblem_picef.h
 * @brief  Problem data for kidney exchange problem (based on the PICEF formulation for kidney exchange)
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_BENDERS_SUBPROBLEM_PICEF_H
#define PROBDATA_BENDERS_SUBPROBLEM_PICEF_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPKEPdataPICEFCreate(
   SCIP*                 kepscip,            /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs             /**< pointer to arc structures of graph */
  );

/** returns the number of nodes in the KEP instance */
int SCIPKEPdataPICEFGetNumNodes(
   SCIP_PROBDATA*        kepdata            /**< problem data */
   );

/** returns the number of pairs in the KEP instance */
int SCIPKEPdataPICEFGetNumPairs(
   SCIP_PROBDATA*        kepdata            /**< problem data */
   );

/** returns the number of arcs in the KEP instance */
int SCIPKEPdataPICEFGetNumArcs(
   SCIP_PROBDATA*        kepdata            /**< problem data */
   );

/** returns the graph in the KEP instance */
Graph* SCIPKEPdataPICEFGetGraph(
   SCIP_PROBDATA*        kepdata            /**< problem data */
   );

/** returns array of cycle variables */
SCIP_VAR** SCIPKEPdataPICEFGetCyclevars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of cycles */
int SCIPKEPdataPICEFGetNCycles(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of arc variables */
SCIP_VAR** SCIPKEPdataPICEFGetArcvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of positioned arcs */
int SCIPKEPdataPICEFGetNPosarcs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of dummy arc variables */
SCIP_VAR** SCIPKEPdataPICEFGetDummyArcvars(
   SCIP_PROBDATA*        kepdata            /**< problem data */
   );

/** returns array of node constraints */
SCIP_CONS** SCIPKEPdataPICEFGetNodeconss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


#ifdef __cplusplus
}
#endif

#endif
