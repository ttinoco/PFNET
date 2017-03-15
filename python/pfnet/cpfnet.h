#ifndef __PYX_HAVE__pfnet__cpfnet
#define __PYX_HAVE__pfnet__cpfnet


#ifndef __PYX_HAVE_API__pfnet__cpfnet

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(PyObject) *new_Network(Net *);

#endif /* !__PYX_HAVE_API__pfnet__cpfnet */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initcpfnet(void);
#else
PyMODINIT_FUNC PyInit_cpfnet(void);
#endif

#endif /* !__PYX_HAVE__pfnet__cpfnet */
