#ifndef MDUNIT_H
#define MDUNIT_H

//#include <omp.h>

#define SAFE_FILEOPEN(fp,name,type){			\
  if(NULL==(fp=fopen(name,type))){			\
    fprintf(stderr,"error: fopen %s failed\n",name);	\
    exit(0);						\
  }							\
}

#define SAFE_MALLOC(ptr,size,type){		\
  if(NULL==(ptr=(type*)malloc(size))){		\
    fprintf(stderr,"error malloc failed\n");	\
    exit(0);					\
  }						\
}

#endif
