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

/**@file   kidneyParams.h
 * @brief  set up parameters for the kidney exchange problem
 * @author Christopher Hojny
 */

#ifndef KIDNEYPARAMS_H
#define KIDNEYPARAMS_H

/* SCIP include */
#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Set basic SCIP parameters that are relevant for the kidnet exchange problem problem */
extern
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Introduce parameters that are relevant for the kidney exchange problem problem */
extern
SCIP_RETCODE addKidneyParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
