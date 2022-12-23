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

/**@file   solveMasterPICEFProblem.h
 * @brief  routines to solve the master problem based on the PICEF formulation
 * @author Danny Blom
 *
 * This file provides routines to solve the master problem with PICEF formulation.
 */

#ifndef SOLVEMASTERPICEFPROBLEM_H
#define SOLVEMASTERPICEFPROBLEM_H

#include "graph.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** solves the master problem */
SCIP_RETCODE solveMasterPICEFProblem(
   SCIP*                 masterscip,         /**< SCIP data structure for master problem */
   Graph*                graph,              /**< pointer to underlying graph */
   Cycles*               cycles,             /**< pointer to cycle structures of graph */
   PositionedArcs*       posarcs,            /**< pointer to arc structures of graph */
   int                   adversarybound,     /**< bound on adversary attack */
   SCIP_Real             timelimit,          /**< time limit to solve the problem */
   SCIP_Real             initialtime,        /**< time spent in initialization */
   const char*           settings,           /**< Possible name of setting file */
   SCIP_Bool             verbose             /**< whether we print SCIP's logs */
   );

/** add cut based on initial master solution to bendersscip */
SCIP_RETCODE SCIPaddInitialPICEFBendersCut(
   SCIP*                masterscip,          /**< SCIP instance of master problem */
   SCIP*                bendersscip,         /**< SCIP instance of benders problem */
   SCIP_SOL*            sol                  /**< solution of master SCIP */
   );

SCIP_RETCODE SCIPsolveBendersPICEFModel(
   SCIP*                 bendersscip,        /**< SCIP data structure for Benders model */
   SCIP*                 kepscip,            /**< SCIP data structure for the subproblem of Benders (a KEP) */
   SCIP*                 masterscip,         /**< SCIP data structure for the master model */
   SCIP_SOL*             mastersol,          /**< SCIP solution to master problem */
   SCIP_Real             masterobj,          /**< Objective value of the master problem */
   Cycles*               cycles,             /**< Data structure needed for cycle weights when adding Benders cuts */
   PositionedArcs*       posarcs,            /**< Data structure needed for arc weights when adding Benders cuts */
   SCIP_Real             timelimit,          /**< time limit for solving the Benders model */
   SCIP_Bool*            didnotfinish,       /**< pointer to store whether we have hit the time limit */
   SCIP_STATUS*          bendersstatus,      /**< pointer to store the status of the bendersmodel */
   SCIP_Real*            timestage2,         /**< pointer to store time spent in stage 2 */
   SCIP_Real*            timestage3,         /**< pointer to store time spent in stage 3 */
   SCIP_Real*            kepobj              /**< pointer to store value of final iteration kepscip solution */

   );

/** returns whether an arc is attacked based on node attack */
SCIP_Bool SCIParcIsAttacked(
   int*                  attackpattern,      /**< Array of attacked vertices */
   int                   nattacks,           /**< Number of attacked vertices */
   PositionedArcs*       posarcs,            /**< Data structure for the arcs in the graph of suitable length */
   int                   arcidx              /**< Index of arc in arcs data structure */
   );

/** returns whether a dummy arc is attacked based on node attack */
SCIP_Bool SCIPdummyarcIsAttacked(
   int*                  attackpattern,      /**< Array of attacked vertices */
   int                   nattacks,           /**< Number of attacked vertices */
   Graph*                graph,              /**< Data structure for the position indexed arcs in the graph */
   int                   source,             /**< vertex that the arc is pointing out from */
   int                   j                   /**< Index of arc in graph data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
