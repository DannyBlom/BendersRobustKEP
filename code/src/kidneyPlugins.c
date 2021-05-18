/**@file   kidneyPlugins.c
 * @brief  load SCIP plugins for the kidney exchange problem
 * @author Christopher Hojny
 */

#include "kidneyPlugins.h"

#include "scip/scipdefplugins.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"


/** Include basic plugins needed for the kidney exchange problem */
SCIP_RETCODE includeKidneyPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   return SCIP_OKAY;
}
