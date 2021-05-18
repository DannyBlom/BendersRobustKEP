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

/**@file   auxiliaryStructures.c
 * @brief  Declaration of auxiliary data structures and methods to handle them
 * @author Christopher Hojny
 *
 * This file contains the definitions of auxiliary data structures and methods to handle them
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "auxiliaryStructures.h"

#include "scip/scip.h"

/** creates branch-and-bound node */
SCIP_RETCODE createBBNode(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode**              node,               /**< pointer to branch-and-bound node */
   int                   nzerofixings,       /**< number of 0-fixings corresponding to node */
   int                   nonefixings,        /**< number of 1-fixings corresponding to node */
   int                   branchingcand,      /**< index of branching candidate */
   SCIP_Real             objbound            /**< upper bound on objective value (or -1.0 if unknown) */
   )
{
   assert( scip != NULL );
   assert( node != NULL );
   assert( nzerofixings >= 0 );
   assert( nonefixings >= 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, node) );

   (*node)->zerofixings = NULL;
   (*node)->onefixings = NULL;
   (*node)->nzerofixings = 0;
   (*node)->nonefixings = 0;
   (*node)->nmaxzerofixings = 0;
   (*node)->nmaxonefixings = 0;
   (*node)->branchingcand = branchingcand;
   (*node)->objbound = objbound;

   /* allocate memory for initial 0-fixings */
   if ( nzerofixings > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*node)->zerofixings, nzerofixings) );
      (*node)->nmaxzerofixings = nzerofixings;
   }

   /* allocate memory for initial 1-fixings */
   if ( nonefixings > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*node)->onefixings, nonefixings) );
      (*node)->nmaxonefixings = nonefixings;
   }

   return SCIP_OKAY;
}

/** frees branch-and-bound node */
SCIP_RETCODE freeBBNode(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode**              node                /**< pointer to branch-and-bound node */
   )
{
   assert( scip != NULL );
   assert( node != NULL );

   if ( (*node)->nmaxzerofixings > 0 )
   {
      assert( (*node)->zerofixings != NULL );

      SCIPfreeBlockMemoryArray(scip, &(*node)->zerofixings, (*node)->nmaxzerofixings);
   }
   if ( (*node)->nmaxonefixings > 0 )
   {
      assert( (*node)->onefixings != NULL );

      SCIPfreeBlockMemoryArray(scip, &(*node)->onefixings, (*node)->nmaxonefixings);
   }

   SCIPfreeBlockMemory(scip, node);

   return SCIP_OKAY;
}

/** adds a 0-fixing to a branch-and-bound node */
SCIP_RETCODE addZerofixing(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode*               node,               /**< branch-and-bound node */
   int                   zerofixing,         /**< 0-fixing to be added to the branch-and-bound node */
   SCIP_Bool             checkfixing         /**< whether we check if the fixing is already present */
   )
{
   assert( scip != NULL );
   assert( node != NULL );

   if ( checkfixing )
   {
      int i;

      for (i = 0; i < node->nzerofixings; ++i)
      {
         if ( node->zerofixings[i] == zerofixing )
            return SCIP_OKAY;
      }
   }

   /* check whether enough memory is available to store fixing */
   if ( node->nzerofixings == node->nmaxzerofixings )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, node->nzerofixings + 1);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &node->zerofixings, node->nmaxzerofixings, newsize) );
      node->nmaxzerofixings = newsize;
   }

   node->zerofixings[(node->nzerofixings)++] = zerofixing;

   return SCIP_OKAY;
}

/** adds a 1-fixing to a branch-and-bound node */
SCIP_RETCODE addOnefixing(
   SCIP*                 scip,               /**< SCIP instance */
   BBNode*               node,               /**< branch-and-bound node */
   int                   onefixing,          /**< 1-fixing to be added to the branch-and-bound node */
   SCIP_Bool             checkfixing         /**< whether we check if the fixing is already present */
   )
{
   assert( scip != NULL );
   assert( node != NULL );

   if ( checkfixing )
   {
      int i;

      for (i = 0; i < node->nonefixings; ++i)
      {
         if ( node->onefixings[i] == onefixing )
            return SCIP_OKAY;
      }
   }

   /* check whether enough memory is available to store fixing */
   if ( node->nonefixings == node->nmaxonefixings )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, node->nonefixings + 1);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &node->onefixings, node->nmaxonefixings, newsize) );
      node->nmaxonefixings = newsize;
   }

   node->onefixings[(node->nonefixings)++] = onefixing;

   return SCIP_OKAY;
}

/** returns upper bound on the objective value */
SCIP_Real getObjbound(
   BBNode*               node                /**< branch-and-bound node */
   )
{
   assert( node != NULL );

   return node->objbound;
}
