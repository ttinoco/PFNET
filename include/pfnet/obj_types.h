/** @file obj_types.h
 *  @brief This file defines object types.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __OBJECTS_HEADER__
#define __OBJECTS_HEADER__

// Object types
/** \defgroup object_types Object Types
 *  @{
 */
#define OBJ_BUS 0      /**< @brief Object type: Bus */
#define OBJ_GEN 1      /**< @brief Object type: Generator */
#define OBJ_BRANCH 2   /**< @brief Object type: Branch */
#define OBJ_SHUNT 3    /**< @brief Object type: Shunt */
#define OBJ_LOAD 4     /**< @bried Object type: Load */
#define OBJ_VARGEN 5   /**< @bried Object type: Variable generator */
#define OBJ_BAT 6      /**< @bried Object type: Battery */
#define OBJ_UNKNOWN 7  /**< @brief Object type: Unknown */
/** @} */

#endif
