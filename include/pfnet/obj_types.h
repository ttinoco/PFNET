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
#define OBJ_ALL 0       /**< @brief Object type: All. */
#define OBJ_BUS 1       /**< @brief Object type: Bus. */
#define OBJ_GEN 2       /**< @brief Object type: Generator. */
#define OBJ_BRANCH 3    /**< @brief Object type: Branch. */
#define OBJ_SHUNT 4     /**< @brief Object type: Shunt. */
#define OBJ_LOAD 5      /**< @brief Object type: Load. */
#define OBJ_VARGEN 6    /**< @brief Object type: Variable generator. */
#define OBJ_BAT 7       /**< @brief Object type: Battery. */
#define OBJ_CONVCSC 8   /**< @brief Object type: CSC converter. */
#define OBJ_CONVVSC 9   /**< @brief Object type: VSC converter. */
#define OBJ_BUSDC 10    /**< @brief Object type: DC bus. */
#define OBJ_BRANCHDC 11 /**< @brief Object type: DC branch. */
#define OBJ_FACTS 12    /**< @brief Object typeL FACTS device */
#define OBJ_UNKNOWN 13  /**< @brief Object type: Unknown. */
/** @} */

#endif
