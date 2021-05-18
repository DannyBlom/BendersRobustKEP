/**@file   kidneyPlugins.h
 * @brief  load SCIP plugins for the kidney exchange problem
 * @author Christopher Hojny
 */

#ifndef KIDNEYPLUGINS_H
#define KIDNEYPLUGINS_H

/* SCIP include */
#include <scip/scip.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Include basic plugins needed for the kidney exchange problem */
extern
SCIP_RETCODE includeKidneyPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
