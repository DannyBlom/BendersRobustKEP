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

/**@file   auxiliaryStructures.h
 * @brief Declaration of auxiliary data structures and methods to handle them
 * @author Christopher Hojny
 *
 * This file contains the definitions of auxiliary data structures and methods to handle them
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef AUXILIARYSTRUCTURES_H
#define AUXILIARYSTRUCTURES_H

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of a node */
typedef struct BBNode
{
   int*                  zerofixings;        /**< array of 0-fixings (identified by variable indices) */
   int*                  onefixings;         /**< array of 1-fixings (identified by variable indices) */
   int                   nzerofixings;       /**< number of 0-fixings */
   int                   nonefixings;        /**< number of 1-fixings */
   int                   nmaxzerofixings;    /**< maximum number of 0-fixings fitting into zerofixings */
   int                   nmaxonefixings;     /**< maximum number of 1-fixings fitting into onefixings */
   int                   branchingcand;      /**< index of branching candidate */
   SCIP_Real             objbound;           /**< upper bound on objective value (or -1.0 if unknown) */
} BBNode;

/** creates branch-and-bound node */
SCIP_RETCODE createBBNode(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode**              node,               /**< pointer to branch-and-bound node */
   int                   nzerofixings,       /**< number of 0-fixings corresponding to node */
   int                   nonefixings,        /**< number of 1-fixings corresponding to node */
   int                   branchingcand,      /**< index of branching candidate */
   SCIP_Real             objbound            /**< upper bound on objective value (or -1.0 if unknown) */
   );

/** frees branch-and-bound node */
SCIP_RETCODE freeBBNode(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode**              node                /**< pointer to branch-and-bound node */
   );

/** adds a 0-fixing to a branch-and-bound node */
SCIP_RETCODE addZerofixing(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode*               node,               /**< branch-and-bound node */
   int                   zerofixing,         /**< 0-fixing to be added to the branch-and-bound node */
   SCIP_Bool             checkfixing         /**< whether we check if the fixing is already present */
   );

/** adds a 1-fixing to a branch-and-bound node */
SCIP_RETCODE addOnefixing(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode*               node,               /**< branch-and-bound node */
   int                   onefixing,          /**< 1-fixing to be added to the branch-and-bound node */
   SCIP_Bool             checkfixing         /**< whether we check if the fixing is already present */
   );

/** returns the array of 0-fixings of a branch-and-bound node */
int* getZerofixings(
   BBNode*               node                /**< branch-and-bound node */
   );

/** returns the array of 1-fixings of a branch-and-bound node */
int* getOnefixings(
   BBNode*               node                /**< branch-and-bound node */
   );

/** returns number of 0-fixings of a branch-and-bound node */
int getNZerofixings(
   BBNode*               node                /**< branch-and-bound node */
   );

/** returns number of 1-fixings of a branch-and-bound node */
int getNOnefixings(
   BBNode*               node                /**< branch-and-bound node */
   );

/** returns upper bound on the objective value */
SCIP_Real getObjbound(
   BBNode*               node                /**< branch-and-bound node */
   );

#ifdef __cplusplus
}
#endif

#endif
