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

/**@file   probdata_benders_picef.h
 * @brief  Problem data for Benders type approach for PICEF formulation
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_BENDERS_PICEF_H
#define PROBDATA_BENDERS_PICEF_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPbendersdataPICEFCreate(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 masterscip,         /**< SCIP data structure of master problem */
   SCIP_SOL*             mastersol,          /**< SCIP solution structure */
   const char*           probname,           /**< problem name */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to position indexed arc structure of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** adds arc variable and its corresponding constraint to the bendersdata */
SCIP_RETCODE SCIPbendersdataPICEFAddArcVars(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP*                 kepscip,            /**< sub SCIP instance for solving 3rd stage */
   SCIP_SOL*             kepsol,             /**< solution of the sub SCIP instance */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Bool             mastersolution      /**< Are we adding a master solution (true) or a kep solution (false)? */
   );

/** adds solution constraint to the bendersdata */
SCIP_RETCODE SCIPbendersdataPICEFAddSolCons(
   SCIP*                 bendersscip,        /**< SCIP data structure */
   SCIP_Real*            vals,               /**< coefficients of variables to add */
   SCIP_VAR**            vars,               /**< variables to add */
   int                   nvars               /**< number of variables to add */
   );

/** returns array of cycle constraints */
SCIP_CONS** SCIPbendersdataPICEFGetCyclelbconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of cycle variables */
SCIP_VAR** SCIPbendersdataPICEFGetCyclevars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of cycles */
int SCIPbendersdataPICEFGetNCycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of arc constraints */
SCIP_CONS** SCIPbendersdataPICEFGetArcconss(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of arc variables */
SCIP_VAR** SCIPbendersdataPICEFGetArcvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of arc variables */
int SCIPbendersdataPICEFGetNArcvars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of arc indices */
int* SCIPbendersdataPICEFGetArcIndices(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of position indexed arcs */
int SCIPbendersdataPICEFGetNPosarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns uvars */
SCIP_VAR** SCIPbendersdataPICEFGetUVars(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns yvarcovered variables */
SCIP_VAR** SCIPbendersdataPICEFGetYvarcovered(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of initial cycles */
int* SCIPbendersdataPICEFGetInitcycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns the number of initial cycles */
int SCIPbendersdataPICEFGetNinitcycles(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns array of initial position indexed arcs */
int* SCIPbendersdataPICEFGetInitposarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of initial position indexed arcs */
int SCIPbendersdataPICEFGetNinitposarcs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns the budget for the adversary */
int SCIPbendersdataPICEFGetAdversaryBound(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns the number of nodes in the instance */
int SCIPbendersdataPICEFGetNumNodes(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns number of pairs in instance */
int SCIPbendersdataPICEFGetNumPairs(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );

/** returns objvar */
SCIP_VAR* SCIPbendersdataPICEFGetObjvar(
   SCIP_PROBDATA*        bendersdata         /**< problem data */
   );



#ifdef __cplusplus
}
#endif

#endif
