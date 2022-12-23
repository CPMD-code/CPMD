#include "cpmd_global.h"

MODULE rotate_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE nvtx_utils
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotate
  !!public :: rotate_c
  !!public :: rotate_da
  PUBLIC :: rottr

CONTAINS

  ! ==================================================================
  SUBROUTINE rottr(a,c1,gam,transa,nstate,n,tlsd,na,nb)
    ! ==--------------------------------------------------------------==
    ! ==         C1 <= A*C1*GAM                                       ==
    ! ==   SPECIAL CASE FOR GAM UPPER TRIAGONAL                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a
    COMPLEX(real_8)                          :: c1(:,:)
    REAL(real_8)                             :: gam(:,:)
    CHARACTER(len=*)                         :: transa
    INTEGER                                  :: nstate, n
    LOGICAL                                  :: tlsd
    INTEGER                                  :: na, nb

    INTEGER                                  :: isub, naa

!(nstate,nstate)
! Variables
! ==--------------------------------------------------------------==

    IF (n.EQ.0 .OR. nstate.EQ.0) RETURN
    CALL tiset('     ROTTR',isub)
    IF (tlsd) THEN
       naa=na+1
       CALL dtrmm('R','U',transa,'N',2*n,na,a,gam,nstate,c1,2*n)
       CALL dtrmm('R','U',transa,'N',2*n,nb,a,gam(naa,naa),nstate,&
            c1(1,naa),2*n)
    ELSE
       CALL dtrmm('R','U',transa,'N',2*n,nstate,a,gam,nstate,c1,2*n)
    ENDIF
    CALL tihalt('     ROTTR',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE rotate(a,c1,b,c2,gam,nstate,n,tlsd,na,nb)
    ! ==--------------------------------------------------------------==
    ! ==         C2 <= B*C2 + A*C1*GAM                                ==
    ! ==   ALSO FOR LSD                                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a
    COMPLEX(real_8)                          :: c1(:,:)
    REAL(real_8)                             :: b
    COMPLEX(real_8)                          :: c2(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gam(nstate,*)
    INTEGER                                  :: n
    LOGICAL                                  :: tlsd
    INTEGER                                  :: na, nb

    CHARACTER(*), PARAMETER                  :: procedureN = 'rotate'

    INTEGER                                  :: isub, naa

! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    IF (n>0) THEN
       IF (tlsd) THEN
          naa=na+1
          CALL dgemm('N','N',n,na,na,a,c1(1,1),n,gam(1,1),nstate,b,c2(1,1),n)
          CALL dgemm('N','N',n,nb,nb,a,c1(1,naa),n,gam(naa,naa),nstate,b,c2(1,naa),n)
       ELSE
          CALL dgemm('N','N',n,nstate,nstate,a,c1,n,gam,nstate,b,c2,n)
       ENDIF
    ENDIF

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rotate



END MODULE rotate_utils

! ==================================================================
SUBROUTINE rotate_da(a,c1,b,c2,gam,ld,n,nstate,ndd1,ndd2,nddx,mepos,pgroup,nproc,grp,tlsd,na,nb)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes, cuda_dealloc, &
       cuda_memcpy_host_to_device, cuda_memcpy_device_to_host, cuda_d_points_to, cuda_mem_zero_bytes
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_create, cublas_destroy, cublas_dgemm, cublas_dger
  USE sizeof_kinds,                    ONLY: sizeof_real_8
  USE cp_cuwfn_types, ONLY: cp_cuwfn_device_get_ptrs, cp_cuwfn, cp_cuwfn_get
  USE nvtx_utils
  IMPLICIT NONE
  REAL(real_8)                               :: a, b
  INTEGER                                    :: ld, n, nstate
  REAL(real_8)                               :: c1(ld,nstate), c2(ld,nstate), &
                                                gam(nstate,nstate)
  INTEGER                                    :: ndd1(0:*), ndd2(0:*), nddx, &
                                                mepos, pgroup(0:*), nproc, grp
  LOGICAL                                    :: tlsd
  INTEGER                                    :: na, nb

  CHARACTER(*), PARAMETER                    :: procedureN = 'rotate_da'
  LOGICAL, PARAMETER                         :: use_gpu = .FALSE.

  INTEGER                                    :: i, i1, ierr, ii, ip, isub, j, &
                                                n1, nmax, nmin
  INTEGER(int_8)                             :: n_bytes
  REAL(real_8), ALLOCATABLE                  :: aux(:)
  REAL(real_8), ALLOCATABLE, DIMENSION(:, :) :: c2_tmp
  TYPE(cublas_handle_t)                      :: blas_handle
  TYPE(cuda_memory_t)                        :: c1_d, c2_d, gam_d

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  __NVTX_TIMER_START ( procedureN )

  IF (tlsd) THEN
     IF (ndd1(mepos).LT.na) THEN
        nmax=MIN(na,ndd2(mepos))
        DO i=ndd1(mepos),nmax
           ii=i-ndd1(mepos)+1
           DO j=na+1,nstate
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
     IF (ndd2(mepos).GT.na) THEN
        nmin=MAX(na+1,ndd1(mepos))
        DO i=nmin,ndd2(mepos)
           ii=i-ndd1(mepos)+1
           DO j=1,na
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  CALL dscal(nstate*n,b,c2,1)
  ALLOCATE(aux(nstate*nddx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  CALL dcopy(nstate*nddx,gam,1,aux,1)

  IF( use_gpu ) THEN
     WRITE(*,*) 'rotate_utils.mod.F90: use_gpu',use_gpu
     CALL cublas_create ( blas_handle, 0 )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c1_d, n_bytes, 0 )
     CALL cuda_memcpy_host_to_device ( c1, c1_d )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c2_d, n_bytes, 0 )
     CALL cuda_mem_zero_bytes ( c2_d, n_bytes )

     n_bytes = SIZE( gam, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( gam_d, n_bytes, 0 )
  ENDIF

  DO ip=0,nproc-1
     i1=ndd1(ip)
     n1=ndd2(ip)-ndd1(ip)+1
     IF (n1.GT.0) THEN
        CALL dcopy(nstate*nddx,aux,1,gam,1)
        CALL mp_bcast(gam,nstate*nddx,pgroup(ip+1),grp)
        IF (ld>0) THEN
           IF( use_gpu ) THEN
              CALL cuda_memcpy_host_to_device ( gam, gam_d )
              CALL cublas_dgemm ( blas_handle, 'N', 'N', n, n1, nstate, &
                   & a, c1_d, ld, &
                   & gam_d, nstate, &
                   & 1.0_real_8, cuda_d_points_to(c2_d,(i1-1)*ld+1), ld )
           ELSE
              CALL dgemm('N','N',n,n1,nstate,a,c1(1,1),ld,&
                   gam(1,1),nstate,1._real_8,c2(1,i1),ld)
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  IF( use_gpu ) THEN
     ALLOCATE(c2_tmp(SIZE(c2,1),SIZE(c2,2)),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     CALL cuda_memcpy_device_to_host ( c2_d, c2_tmp )
     CALL daxpy ( SIZE(c2), 1.0_real_8, c2_tmp, 1, c2, 1 )
     DEALLOCATE(c2_tmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
     CALL cuda_dealloc ( gam_d )
     CALL cuda_dealloc ( c2_d )
     CALL cuda_dealloc ( c1_d )
     CALL cublas_destroy ( blas_handle )
  ENDIF

  CALL dcopy(nstate*nddx,aux,1,gam,1)
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  __NVTX_TIMER_STOP
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_da
! ==================================================================
SUBROUTINE rotate_c(a,c1,b,c2,gam,nstate)
  ! ==--------------------------------------------------------------==
  ! ==         C2 <= B*C2 + A*C1*GAM                                ==
  ! ==   ALSO FOR LSD                                               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE spin , ONLY:spin_mod
  IMPLICIT NONE
  COMPLEX(real_8)                            :: a, c1(nkpt%ngwk,*), b, &
                                                c2(nkpt%ngwk,*)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: gam(nstate,*)

  INTEGER                                    :: isub, naa

  IF (nkpt%ngwk.EQ.0 .OR. nstate.EQ.0) RETURN
  CALL tiset('  ROTATE_C',isub)
  IF (cntl%tlsd) THEN
     naa=spin_mod%nsup+1
     IF (nkpt%ngwk>0) THEN
        CALL zgemm('N','N',nkpt%ngwk,spin_mod%nsup,spin_mod%nsup,a,c1(1,1),nkpt%ngwk,gam(1,1),&
             nstate,b,c2(1,1),nkpt%ngwk)
        CALL zgemm('N','N',nkpt%ngwk,spin_mod%nsdown,spin_mod%nsdown,a,c1(1,naa),nkpt%ngwk,&
             gam(naa,naa),nstate,b,c2(1,naa),nkpt%ngwk)
     ENDIF
  ELSE
     IF (nkpt%ngwk>0) THEN
        CALL zgemm('N','N',nkpt%ngwk,nstate,nstate,a,c1(1,1),nkpt%ngwk,&
             gam(1,1),nstate,b,c2(1,1),nkpt%ngwk)
     ENDIF
  ENDIF
  CALL tihalt('  ROTATE_C',isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_c
