/** @file flag_types.h
 *  @brief This file defines flag types.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FLAGS_HEADER__
#define __FLAGS_HEADER__

// Flag helpers
/** \defgroup flag_helpers Flag Helpers
 *  @{
 */
#define ALL_VARS 0xFF    /**< @brief All variables. */
/** @} */

// Flag types
/** \defgroup flag_types Flag Types
 *  @{
 */
#define FLAG_VARS 0x01    /**< @brief Flag type: variables */
#define FLAG_FIXED 0x02   /**< @brief Flag type: fixed */
#define FLAG_BOUNDED 0x04 /**< @brief Flag type: bounded */
#define FLAG_SPARSE 0x08  /**< @brief Flag type: sparse */
/** @} */

#endif
