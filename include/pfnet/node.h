/** @file node.h
 *  @brief This file lists the constants and routines associated with the Node data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2019, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __NODE_HEADER__
#define __NODE_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"


// Node
typedef struct Node Node;

// Node functions
void* NODE_get_item(Node* node);
Node* NODE_get_next(Node* node);
Node* NODE_list_add(Node* node_list, Node* node);
int NODE_list_len(Node* node_list);
void NODE_list_del(Node* node_list);
Node* NODE_new(void);
void NODE_set_item(Node* node, void* item);

#endif
