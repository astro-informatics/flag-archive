#ifndef FLAG_ERROR
#define FLAG_ERROR

#include <stdio.h>

#define FLAG_ERROR_GENERIC(comment) 					\
  printf("ERROR: %s.\n", comment);					\
  printf("ERROR: %s <%s> %s %s %s %d.\n",				\
	 "Occurred in function",					\
	   __PRETTY_FUNCTION__,						\
	   "of file", __FILE__,						\
	   "on line", __LINE__);					\
  exit(1);

#define FLAG_ERROR_MEM_ALLOC_CHECK(pointer)				\
  if(pointer == NULL) {							\
    FLAG_ERROR_GENERIC("Memory allocation failed")			\
  }


#endif