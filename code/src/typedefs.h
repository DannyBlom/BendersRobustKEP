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

#ifndef __TYPEDEFS_H_
#define __TYPEDEFS_H_

#ifdef __cplusplus
extern "C" {
#endif

/* methods to solve the second stage problem */
#define METHOD_BENDERS                 1      /**< Benders model in cycle formulation */
#define METHOD_BENDERS_PICEF           2      /**< Benders model in PICEF formulation */
#define METHOD_BRANCHANDBOUND          3      /**< branch-and-bound method used by Glorie et al. */

/* recourse policies in the third stage problem */
#define POLICY_FULLRECOURSE            1      /**< Full recourse policy (no extra restrictions on the recourse solution) */
#define POLICY_KEEPUNAFFECTEDCC        2      /**< Policy in which recourse solution preserves cycles / chains from initial sol that are not affected by attack */
#define POLICY_GUARANTEEUNAFFECTEDCC   3      /**< Policy in which recourse solution guarantees that patients of unaffected cycles/chains receive a transplant */

/* attack selection options (Benders only) */
#define ATTACK_SELECTION_FIRSTVIOLATED 1      /**< Select first scenario obtained that violates the master solution */
#define ATTACK_SELECTION_MOSTVIOLATED  2      /**< Selects the scenario that has the largest violation of the master solution */


#ifdef __cplusplus
}
#endif

#endif
