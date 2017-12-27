/** @file list.h
 *  @brief This file defines MACROS for manipulating linked lists.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __LIST_HEADER__
#define __LIST_HEADER__

#define LIST_push(head, item, next) {		\
  if ((item)) { 				\
    (item)->next = (head);			\
    (head) = (item);				\
  }						\
}

#define LIST_add(type, head, item, next) {\
  type* _h_;				  \
  type* _prev_ = NULL;			  \
  char flag = 0;			  \
  if ((item))				  \
    (item)->next = NULL;		  \
  if (!(head))				  \
    (head) = (item);			  \
  else {								\
    for (_h_ = (head); _h_ != NULL; _prev_ = _h_, _h_ = (_h_)->next) {	\
      if (_h_ == (item)) {						\
	flag = 1;							\
	break;								\
      }									\
    }									\
    if (!flag)								\
      (_prev_)->next = (item);						\
  }									\
}

#define LIST_del(type, head, item, next) {\
  type* _h_;				  \
  type* _prev_;				  \
  if (!(item) || !(head))		  \
    ;					  \
  else if ((head) == (item))		  \
    (head) = (item)->next;		  \
  else {				  \
    _h_ = (head)->next;			  \
    _prev_ = (head);			  \
    while (_h_) {			  \
      if (_h_ == (item))		  \
	(_prev_)->next = (_h_)->next;	  \
      _prev_ = _h_;			  \
      _h_ = (_h_)->next;		  \
    }					  \
  }					  \
}

#define LIST_len(type, head, next, len) {\
  type* _h_;				 \
  (len) = 0;						\
  for (_h_ = (head); _h_ != NULL; _h_ = (_h_)->next) {	\
    (len)++;						\
  }							\
}

#define LIST_map(type, head, var, next, func) {\
  type* _h_;				       \
  type* var;				       \
  (var) = (head);			       \
  while ((var) != NULL) {		       \
    _h_ = (var)->next;			       \
    {func;};				       \
    (var) = _h_;			       \
  }					       \
}

#endif
