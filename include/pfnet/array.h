/** @file array.h
 *  @brief This file defines MACROS for manipulating arrays.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __ARRAY_HEADER__
#define __ARRAY_HEADER__

#define ARRAY_alloc(ar, type, num) {	      \
    (ar) = (type*)malloc((num)*sizeof(type)); \
}

#define ARRAY_zalloc(ar, type, num) {	      \
    (ar) = (type*)calloc((num),sizeof(type)); \
}

#define ARRAY_clear(ar, type, num) {   \
    memset((ar),0,(num)*sizeof(type)); \
}

#endif
