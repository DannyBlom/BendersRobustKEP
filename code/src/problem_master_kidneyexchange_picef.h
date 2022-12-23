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

/**@file   problem_master_kidneyexchange_picef.h
 * @brief  Problem data for the master kidney exchange problem with position indexed cycle edge formulation
 * @author Danny Blom
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef PROBLEM_MASTER_KIDNEYEXCHANGE_PICEF_H
#define PROBLEM_MASTER_KIDNEYEXCHANGE_PICEF_H

#include "graph.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates initial model of the master kidney exchange problem */
extern
SCIP_RETCODE SCIPcreateMasterPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< pointer to kidney exchange graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to arc structures of graph */
   int                   adversarybound      /**< bound on adversary attack */
   );

/** free master kidney exchange problem data */
extern
SCIP_RETCODE SCIPfreeMasterPICEFModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Graph*                graph,              /**< underlying graph */
   Cycles*               cycles,             /**< cycle structures of graph */
   PositionedArcs*       posarcs             /**< arc structures of graph */
   );

#ifdef __cplusplus
}
#endif

#endif
