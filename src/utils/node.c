/** @file node.c
 *  @brief This file defines the Node data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2019, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/node.h>

struct Node {
  void* item;
  Node* next;
};

void* NODE_get_item(Node* node) {
  if (!node)
    return NULL;
  else
    return node->item;
}

Node* NODE_get_next(Node* node) {
  if (!node)
    return NULL;
  else
    return node->next;
}

Node* NODE_list_add(Node* node_list, Node* node) {
  LIST_add(Node,node_list,node,next);
  return node_list;
}

int NODE_list_len(Node* node_list) {
  int len;
  LIST_len(Node,node_list,next,len);
  return len;
}

void NODE_list_del(Node* node_list) {
  LIST_map(Node,node_list,node,next,{free(node);});
}

Node* NODE_new(void) {
  Node* node = (Node*)malloc(sizeof(Node));
  node->item = NULL;
  node->next = NULL;
  return node;
}

void NODE_set_item(Node* node, void* item) {
  if (node)
    node->item = item;
}
