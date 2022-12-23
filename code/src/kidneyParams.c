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

/**@file   kidneyParams.c
 * @brief  set up parameters for the kidney exchange problem
 * @author Christopher Hojny
 * @author Danny Blom
 */

#include "kidneyParams.h"

/* default settings for parameters - for documentation see below */

/* information regarding problem */
#define DEFAULT_PRINTGRAPHINFO               FALSE
#define DEFAULT_MAXCYCLELENGTH               3
#define DEFAULT_MAXCHAINLENGTH               2
#define DEFAULT_CYCHAINCONSSINITIAL          TRUE
#define DEFAULT_USETWOTHIRDCLIQUES           FALSE
#define DEFAULT_MINARCCOUNTTRIPLET           4
#define DEFAULT_METHOD                       1           /**< Method options (1: Benders, 2: Benders + PICEF, 3: Glorie et al.) */
#define DEFAULT_LIFTINGBENDERSCUTS           TRUE
#define DEFAULT_RECOURSEPOLICY               2           /**< Method options (1: full recourse, 2: keep unaffected CC) */
#define DEFAULT_COMPLEXCCCUTS                FALSE       
#define DEFAULT_ATTACKSELECTION              1           /**< Method options (1: first violated attack, 2: most violated attack) */

/** Set basic SCIP parameters that are relevant for the kidney exchange  problem */
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 0) );

   return SCIP_OKAY;
}


/** Introduce parameters that are relevant for the kidney exhange problem */
SCIP_RETCODE addKidneyParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   /* information regarding problem */
   SCIP_CALL( SCIPaddBoolParam(scip, "kidney/printgraphinfo", "Shall information on the underlying graph be printed?",
         NULL, TRUE, DEFAULT_PRINTGRAPHINFO, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/maxcyclelength", "Maximum length of cycles to be considered in model.",
         NULL, TRUE, DEFAULT_MAXCYCLELENGTH, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/maxchainlength", "Maximum length of chains to be considered in model. Counted in number of transplants",
         NULL, TRUE, DEFAULT_MAXCHAINLENGTH, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "kidney/cyclechainconssinitial", "Whether all cycle/chain constraints are initial cuts.",
         NULL, TRUE, DEFAULT_CYCHAINCONSSINITIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "kidney/usetwothirdcliques", "Shall we use clique inequalities based on cycles / chains including at least two out of third KEP nodes?",
         NULL, TRUE, DEFAULT_USETWOTHIRDCLIQUES, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/minarccounttriplet", "Minimum number of arcs in induced subgraph of a triplet.",
         NULL, TRUE, DEFAULT_MINARCCOUNTTRIPLET, 0, 6, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/method", "Which method do we use to solve the master problem?",
         NULL, TRUE, DEFAULT_METHOD, 1, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "kidney/liftbenderscuts", "Shall we lift KEP solutions on a subgraph to a KEP solution on the entire graph?",
         NULL, TRUE, DEFAULT_LIFTINGBENDERSCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/recoursepolicy", "Which recourse policy do we use?",
         NULL, TRUE, DEFAULT_RECOURSEPOLICY, 1, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "kidney/complexcccuts", "Do we want to use the complex version of CC cuts? (imitating PICEF, only effective for both Benders CC)",
         NULL, TRUE, DEFAULT_COMPLEXCCCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "kidney/attackselection", "Do we want to select first violated scenario (1) or most violated scenario (2)? (only relevant for Benders CC and PICEF)",
         NULL, TRUE, DEFAULT_ATTACKSELECTION, 1, 2, NULL, NULL) );

   return SCIP_OKAY;
}
