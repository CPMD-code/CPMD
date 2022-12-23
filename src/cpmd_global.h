!vw ------------------------------------------------------
!vw Compiler

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_CONTIGUOUS)
#define __CONTIGUOUS
#else
#define __CONTIGUOUS , CONTIGUOUS
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_ASYNCHRONOUS)
#define __ASYNCHRONOUS
#else
#define __ASYNCHRONOUS , ASYNCHRONOUS
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F03_ISO_FORTRAN_ENV)
#define _HASNT_F03_ISO_FORTRAN_ENV
#define _HASNT_F08_ISO_FORTRAN_ENV
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_F08_ISO_FORTRAN_ENV)
#define _HASNT_F08_ISO_FORTRAN_ENV
#endif

!vw whether the compiler supports command line
#if defined(__HASNT_F03_EXECUTE_COMMAND_LINE)
#define _HASNT_F03_EXECUTE_COMMAND_LINE
#endif

#if defined(__HASNT_F08_POINTER_REMAPPING)
#define _HASNT_F08_POINTER_REMAPPING
#endif

!vw whether the mpi library support standard 3.0
#if defined(__HASNT_MPI_30)
#define _HASNT_MPI_30
#endif

!vw whether the compiler supports the attribute
#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#define __COLLAPSE3
#else
#define __COLLAPSE2 collapse(2)
#define __COLLAPSE3 collapse(3)
#endif

!vw whether the compiler supports openmp 4.5
#if defined(__HASNT_OMP_45)
#define _HASNT_OMP_45
#endif

!vw whether the compiler supports omp_set/get_nested
#if defined(__HASNT_OMP_SET_NESTED)
#define _HASNT_OMP_SET_NESTED
#endif

!vw whether the compiler supports stream io
#if defined(__HASNT_BF_STREAM_IO)
#define _HASNT_BF_STREAM_IO
#endif

!vw bypass multithreading if not supported
#if defined(__HASNT_MULTITHREAD_MPI_SUPPORT)
#define _HASNT_MULTITHREAD_MPI_SUPPORT
#endif

!vw whether to use the external c math library
#if defined(__HAS_EXTERNAL_C_ERF)
#define _HAS_EXTERNAL_C_ERF
#endif

!vw ------------------------------------------------------
!vw EXTERNAL LIBS
#if defined(__HAS_LIBXC)
#define _HAS_LIBXC
#endif


!vw ------------------------------------------------------
!vw CUDA 

!vw whether CUDA is supported
#if defined(__HAS_CUDA)
#define _HAS_CUDA
#endif

!vw whether timer are used
#if defined(__HAS_NVTX_TIMER)
#define __NVTX_TIMER_START( name )  CALL nvtx_range_push_a ( name )
#define __NVTX_TIMER_STOP  CALL nvtx_range_pop ( )
#else
#define __NVTX_TIMER_START( name )
#define __NVTX_TIMER_STOP
#endif
