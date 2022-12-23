#include "cpmd_global.h"

MODULE ovlap_utils
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn,&
                                             cp_cuwfn_device_get_ptrs,&
                                             cp_cuwfn_get
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_redist
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_create,&
                                             cublas_destroy,&
                                             cublas_dgemm,&
                                             cublas_dger
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_dealloc,&
                                             cuda_memcpy_device_to_host,&
                                             cuda_memcpy_host_to_device,&
                                             cuda_z_points_to
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE nvtx_utils
  USE parac,                           ONLY: parai,&
                                             paral
  USE sizeof_kinds,                    ONLY: sizeof_complex_8,&
                                             sizeof_real_8
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ovlap
  !!public :: ovlap2
  PUBLIC :: dmatc
  !!public :: ovlap_c
  !!public :: ovlap_h
  !!public :: ovlap2_c
  PUBLIC :: ovlap_add
  PUBLIC :: ovlap2_dist

CONTAINS

  ! ==================================================================
  SUBROUTINE ovlap(nstate,a,c1,c2)
    ! ==--------------------------------------------------------------==
    ! ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a(nstate,nstate)
    COMPLEX(real_8), TARGET                  :: c1(:,:), c2(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ovlap'

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C1_local, C2_local
    COMPLEX(real_8), POINTER                 :: pc1, pc2
    INTEGER                                  :: ibeg_c0, ierr, isub, isub2, &
                                                isub3, NGW_local
    LOGICAL                                  :: GEQ0_local, symmetric

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    ! ..Check for symmetric overlap
    symmetric=.FALSE.
    pc1 => c1(1,1)
    pc2 => c2(1,1)
    symmetric = ASSOCIATED(pc1, TARGET=pc2)
    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_a',isub2)
    CALL cp_grp_get_sizes(ngw_l=NGW_local,geq0_l=GEQ0_local,&
         first_g=ibeg_c0)
    ALLOCATE(C1_local(NGW_local,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL cp_grp_copy_wfn_to_local(c1,ncpw%ngw,C1_local,NGW_local,&
         ibeg_c0,NGW_local,nstate)
    IF (.NOT.symmetric) THEN
       ALLOCATE(C2_local(NGW_local,nstate),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL cp_grp_copy_wfn_to_local(c2,ncpw%ngw,C2_local,NGW_local,&
            ibeg_c0,NGW_local,nstate)
    ENDIF
    CALL tihalt(procedureN//'_grps_a',isub2)
    ! <<<<<<<

    IF (NGW_local.GT.0.AND.nstate.GT.0) THEN
       IF (cntl%tlsd) THEN
          CALL zeroing(a)!,nstate*nstate)
          ! ==--------------------------------------------------------------==
          ! ..Alpha spin
          IF (symmetric) THEN
             CALL dsyrk('U','T',spin_mod%nsup,2*NGW_local,2._real_8,C1_local,2*ngw_local,&
                  0._real_8,a,nstate)
             CALL dmatc('U',spin_mod%nsup,a,nstate)
          ELSE
             CALL dgemm('T','N',spin_mod%nsup,spin_mod%nsup,2*NGW_local,2.0_real_8,C1_local(1,1),&
                  2*NGW_local,C2_local(1,1),2*ngw_local,0.0_real_8,a(1,1),&
                  nstate)
          ENDIF
          IF (GEQ0_local) THEN
             IF (symmetric) THEN
                CALL dger(spin_mod%nsup,spin_mod%nsup,-1.0_real_8,C1_local(1,1),2*NGW_local,&
                     C1_local(1,1),2*NGW_local,a(1,1),nstate)
             ELSE
                CALL dger(spin_mod%nsup,spin_mod%nsup,-1.0_real_8,C1_local(1,1),2*NGW_local,&
                     C2_local(1,1),2*NGW_local,a(1,1),nstate)
             ENDIF
          ENDIF
          ! ==--------------------------------------------------------------==
          ! ..Beta spin
          IF (symmetric) THEN
             CALL dsyrk('U','T',spin_mod%nsdown,2*NGW_local,2._real_8,C1_local(1,spin_mod%nsup+1),&
                  2*NGW_local,0._real_8,a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
             CALL dmatc('U',spin_mod%nsdown,a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
          ELSE
             CALL dgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,2*NGW_local,2.0_real_8,&
                  C1_local(1,spin_mod%nsup+1),2*NGW_local,C2_local(1,spin_mod%nsup+1),&
                  2*NGW_local,0.0_real_8,a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
          ENDIF
          IF (GEQ0_local) THEN
             IF (symmetric) THEN
                CALL dger(spin_mod%nsdown,spin_mod%nsdown,-1.0_real_8,C1_local(1,spin_mod%nsup+1),&
                     2*NGW_local,C1_local(1,spin_mod%nsup+1),2*ngw_local,&
                     a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
             ELSE
                CALL dger(spin_mod%nsdown,spin_mod%nsdown,-1.0_real_8,C1_local(1,spin_mod%nsup+1),&
                     2*NGW_local,C2_local(1,spin_mod%nsup+1),2*ngw_local,&
                     a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
             ENDIF
          ENDIF
       ELSE
          ! ==--------------------------------------------------------------==
          ! ..LDA
          IF (symmetric) THEN
             CALL dsyrk('U','T',nstate,2*NGW_local,2._real_8,C1_local(1,1),&
                  2*NGW_local,0._real_8,a(1,1),nstate)
             CALL dmatc('U',nstate,a,nstate)
          ELSE
             CALL dgemm('T','N',nstate,nstate,2*NGW_local,2.0_real_8,&
                  C1_local(1,1),2*NGW_local,C2_local(1,1),2*ngw_local,&
                  0.0_real_8,a(1,1),nstate)
          ENDIF
          IF (GEQ0_local) THEN
             IF (symmetric) THEN
                CALL dger(nstate,nstate,-1.0_real_8,C1_local(1,1),2*NGW_local,&
                     C1_local(1,1),2*NGW_local,a(1,1),nstate)
             ELSE
                CALL dger(nstate,nstate,-1.0_real_8,C1_local(1,1),2*NGW_local,&
                     C2_local(1,1),2*NGW_local,a(1,1),nstate)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       CALL zeroing(a)!,nstate**2)
    ENDIF

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_b',isub3)
    ! we reduce so that we get back the cp_grp distribution of the matrix
    CALL cp_grp_redist(a,nstate,nstate)
    DEALLOCATE(C1_local,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)
    IF (.NOT.symmetric) THEN
       DEALLOCATE(C2_local,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN//'_grps_b',isub3)
    ! <<<<<<<

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ovlap
  ! ==================================================================
  SUBROUTINE dmatc(uplo,n,a,lda)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: uplo
    INTEGER                                  :: n, lda
    REAL(real_8)                             :: a(lda,n)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    IF (uplo.EQ.'U'.OR.uplo.EQ.'u') THEN
       !$omp parallel do private(I,J)  !schedule(static,1)
       DO i=1,n
          DO j=i+1,n
             a(j,i)=a(i,j)
          ENDDO
       ENDDO
       !$omp end parallel do
    ELSEIF (uplo.EQ.'L'.OR.uplo.EQ.'l') THEN
       !$omp parallel do private(I,J)  !schedule(static,1)
       DO i=1,n
          DO j=i+1,n
             a(i,j)=a(j,i)
          ENDDO
       ENDDO
       !$omp end parallel do
    ELSE
       CALL stopgm('DMATC',' ERROR ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE dmatc
  ! ==================================================================
  SUBROUTINE ovlap_add(nstate,a,c1,c2)
    ! ==--------------------------------------------------------------==
    ! ==         SUMS UP THE OVERLAP MATRIX A = A + < C1 | C2 >       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a(nstate,*)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,*), c2(ncpw%ngw,*)

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    IF (nkpt%ngwk.LE.0) RETURN
    CALL tiset('     OVLAP',isub)
    IF (cntl%tlsd) THEN
       ! ..Alpha spin
       CALL dgemm('T','N',spin_mod%nsup,spin_mod%nsup,2*ncpw%ngw,2.0_real_8,c1(1,1),2*ncpw%ngw,&
            c2(1,1),2*ncpw%ngw,1.0_real_8,a(1,1),nstate)
       IF (Geq0)&
            CALL dger(spin_mod%nsup,spin_mod%nsup,-1.0_real_8,c1(1,1),2*ncpw%ngw,&
            c2(1,1),2*ncpw%ngw,a(1,1),nstate)
       ! ..Beta spin
       CALL dgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,2*ncpw%ngw,2.0_real_8,c1(1,spin_mod%nsup+1),&
            2*ncpw%ngw,c2(1,spin_mod%nsup+1),2*ncpw%ngw,1.0_real_8,&
            a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
       IF (Geq0)&
            CALL dger(spin_mod%nsdown,spin_mod%nsdown,-1.0_real_8,c1(1,spin_mod%nsup+1),2*ncpw%ngw,&
            c2(1,spin_mod%nsup+1),2*ncpw%ngw,a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)

    ELSE
       CALL dgemm('T','N',nstate,nstate,2*ncpw%ngw,2.0_real_8,c1(1,1),2*ncpw%ngw,&
            c2(1,1),2*ncpw%ngw,1.0_real_8,a(1,1),nstate)
       IF (Geq0)&
            CALL dger(nstate,nstate,-1.0_real_8,c1(1,1),2*ncpw%ngw,c2(1,1),&
            2*ncpw%ngw,a(1,1),nstate)
    ENDIF
    CALL tihalt('     OVLAP',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ovlap_add
  ! ==================================================================

  SUBROUTINE ovlap2_dist ( nstate, a, c1, c2 )
    ! computes the distributed overlap matrix A = < C1 | C2 >
    USE system,                          ONLY: cnti, paraw, parap
    USE cp_grp_utils,                    ONLY: cp_grp_get_cp_rank
    USE jrotation_utils,                 ONLY: set_orbdist
    USE mp_interface,                    ONLY: mp_sum
    INTEGER, INTENT(IN)                      :: nstate
    REAL(real_8), INTENT(INOUT)              :: a(:,:)
    COMPLEX(real_8), INTENT(INOUT)           :: c1(:,:), c2(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ovlap2_dist'
    LOGICAL, PARAMETER                       :: use_gpu = .FALSE.

    INTEGER                                  :: ierr, ip, ipp, isub, n, nstx, &
                                                NSTX_grp, nx
    INTEGER(int_8)                           :: n_bytes
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: NWA12_grp
    LOGICAL                                  :: ready, use_c0_on_gpu
    REAL(real_8)                             :: chksum
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a2mat, a3mat
    TYPE(cublas_handle_t)                    :: blas_handle
    TYPE(cuda_memory_t)                      :: a2mat_d, c1_d
    TYPE(cuda_memory_t), POINTER             :: c2_d_p
    TYPE(cuda_memory_t), TARGET              :: c2_d_t

!LOGICAL                                    :: use_cp_grp

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )


    ALLOCATE(NWA12_grp(0:parai%cp_nproc-1,2),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL set_orbdist( nstate,cnti%nstblk,parai%cp_nproc,NSTX_grp)
    NWA12_grp(0:parai%cp_nproc-1,1:2)=paraw%nwa12(0:parai%cp_nproc-1,1:2)
    CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)!vw needed?

    ALLOCATE(a2mat(nstate,nstx),a3mat(nstate,nstx),stat=ierr)!vw nstx too big needs only NSTX_grp???
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(a)

    IF( use_gpu ) THEN
       WRITE(*,*) 'ovlap_utils.mod.F90: use_gpu', use_gpu
       CALL cublas_create ( blas_handle, 0 )

       n_bytes = SIZE( a2mat, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
       CALL cuda_alloc_bytes ( a2mat_d, n_bytes, 0 )

       n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_complex_8, int_8 )
       CALL cuda_alloc_bytes ( c1_d, n_bytes, 0 )
       CALL cuda_memcpy_host_to_device ( c1, c1_d )

       !vw use c0 from GPU if possible.
       !vw to avoid problem, c0_d should pass as argument to the procedure.
       CALL cp_cuwfn_get ( cp_cuwfn, ready=ready, chksum=chksum )
       use_c0_on_gpu = ready .AND. ABS( SUM( ABS( c2(:,1) ) ) - chksum ) < 1.0e-12_real_8
       WRITE(*,*) 'ovlap_utils.mod.F90: use_c0_on_gpu=',use_c0_on_gpu
       IF( use_c0_on_gpu ) THEN
          CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx=1, c0_d=c2_d_p )
       ELSE
          n_bytes = SIZE( c2, KIND=int_8 ) * INT( sizeof_complex_8, int_8 )
          CALL cuda_alloc_bytes ( c2_d_t, n_bytes, 0 )
          c2_d_p => c2_d_t
          CALL cuda_memcpy_host_to_device ( c2, c2_d_p )
       ENDIF

    ENDIF

    DO ip=0,parai%nproc-1
       ipp = cp_grp_get_cp_rank(ip,parai%cp_inter_me)
       nx = NWA12_grp(ipp,2)-nwa12_grp(ipp,1)+1
       CALL zeroing(a2mat)
       IF (nx.GT.0) THEN
          IF( use_gpu ) THEN
             CALL cublas_dgemm ( blas_handle, 'T', 'N', nstate, nx, 2*ncpw%ngw, &
                  & 2.0_real_8, c1_d, 2*ncpw%ngw, &
                  & cuda_z_points_to(c2_d_p,(NWA12_grp(ipp,1)-1)*ncpw%ngw+1), 2*ncpw%ngw, &
                  & 0.0_real_8, a2mat_d, nstate )
             IF (Geq0) THEN
                CALL cublas_dger ( blas_handle, nstate, nx, -1.0_real_8, &
                     & c1_d, 2*ncpw%ngw, &
                     & cuda_z_points_to(c2_d_p,(NWA12_grp(ipp,1)-1)*ncpw%ngw+1), 2*ncpw%ngw, &
                     & a2mat_d, nstate )
             ENDIF
             CALL cuda_memcpy_device_to_host ( a2mat_d, a2mat )
          ELSE
             CALL ovlap2(ncpw%ngw,nstate,nx,a2mat,c1,c2(1,NWA12_grp(ipp,1)),.FALSE.)
          ENDIF
          CALL mp_sum(a2mat,a3mat,nstate*NSTX_grp,parap%pgroup(ip+1),parai%allgrp)
          IF (parai%me.EQ.parap%pgroup(ip+1)) CALL dcopy(nstate*NSTX_grp,a3mat,1,a,1)
       ENDIF
    ENDDO


    IF( use_gpu ) THEN

       CALL cuda_dealloc ( a2mat_d )
       CALL cuda_dealloc ( c1_d )
       CALL cublas_destroy ( blas_handle )
       IF( .NOT. use_c0_on_gpu ) CALL cuda_dealloc ( c2_d_t )

    ENDIF

    DEALLOCATE(a2mat,a3mat,NWA12_grp,stat=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)

  END SUBROUTINE ovlap2_dist

END MODULE ovlap_utils


! ==================================================================
SUBROUTINE ovlap2_c(ngwk,n1,n2,a,c1,c2)
  ! ==--------------------------------------------------------------==
  ! ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
  ! ==         WITHOUT INVERSION SYMMETRY                           ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: ngwk, n1, n2
  COMPLEX(real_8)                            :: a(n1,n2), c1(ngwk,*), &
                                                c2(ngwk,*)

  COMPLEX(real_8), PARAMETER                 :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

  INTEGER                                    :: isub

  IF (ngwk.LE.0) THEN
     CALL zeroing(a)!,n1*n2)
     RETURN
  ENDIF
  CALL tiset ('  OVLAP2_C',isub)
  IF (n2.GT.1) THEN
     CALL zgemm('C','N',n1,n2,ngwk,zone,c1(1,1),ngwk,c2(1,1),&
          ngwk,zzero,a(1,1),n1)
  ELSE
     CALL zgemv('C',ngwk,n1,zone,c1(1,1),ngwk,c2(1,1),1,&
          zzero,a(1,1),1)
  ENDIF
  CALL tihalt('  OVLAP2_C',isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE ovlap2_c
! ==================================================================
SUBROUTINE ovlap_h(nstate,a,c1,c2)
  ! ==--------------------------------------------------------------==
  ! ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
  ! ==         WITH INVERSION SYMMETRY (HERMITIAN)                  ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:nkpt
  USE parac, ONLY : paral,parai
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: a(nstate,nstate), &
                                                c1(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)

  COMPLEX(real_8)                            :: tmpv
  COMPLEX(real_8), EXTERNAL                  :: zdotc
  INTEGER                                    :: i, isub, j

  IF (nkpt%ngwk.LE.0) THEN
     CALL zeroing(a)!,nstate*nstate)
     RETURN
  ENDIF
  CALL tiset ('   OVLAP_H',isub)
  !$omp parallel do private(J,I,TMPV)
  DO j=1,nstate
     a(j,j)=zdotc(nkpt%ngwk,c1(1,j),1,c2(1,j),1)
     DO i=j+1,nstate
        tmpv=zdotc(nkpt%ngwk,c1(1,i),1,c2(1,j),1)
        a(i,j)=tmpv
        a(j,i)=CONJG(tmpv)
     ENDDO
  ENDDO
  CALL tihalt('   OVLAP_H',isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE ovlap_h
! ==================================================================
SUBROUTINE ovlap2(ngw,n1,n2,a,c1,c2,use_cp_grp)
  ! ==--------------------------------------------------------------==
  ! ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE geq0mod , ONLY:geq0
  USE cp_grp_utils, ONLY : cp_grp_redist
  USE cp_grp_utils, ONLY : cp_grp_get_sizes
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: ngw, n1, n2
  REAL(real_8)                               :: a(n1,n2)
  COMPLEX(real_8)                            :: c1(ngw,*), c2(ngw,*)
  LOGICAL                                    :: use_cp_grp

  CHARACTER(*), PARAMETER                    :: procedureN = 'ovlap2'

  INTEGER                                    :: ibeg_c0, isub, isub3, n, &
                                                NGW_local
  LOGICAL                                    :: GEQ0_local

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)

  ! >>>>>>> cp_grp trick
  IF (use_cp_grp) THEN
     CALL cp_grp_get_sizes(ngw_l=NGW_local,geq0_l=GEQ0_local,&
          first_g=ibeg_c0)
  ELSE
     NGW_local = ngw
     GEQ0_local = Geq0
     ibeg_c0 = 1
  ENDIF
  ! <<<<<<<

  IF (NGW_local.GT.0) THEN
#if defined (__ES)
     ! mb   This is not necessary if your BLAS libraries are optimized for ES 
     CALL dgemm2('T','N',n1,n2,2*NGW_local,2.0_real_8,c1(ibeg_c0,1),&
          2*ngw,c2(ibeg_c0,1),2*ngw,0.0_real_8,a(1,1),n1)
#else
     CALL dgemm('T','N',n1,n2,2*NGW_local,2.0_real_8,c1(ibeg_c0,1),&
          2*ngw,c2(ibeg_c0,1),2*ngw,0.0_real_8,a(1,1),n1)
#endif
     IF (GEQ0_local)&
          CALL dger(n1,n2,-1.0_real_8,c1(ibeg_c0,1),2*ngw,&
          c2(ibeg_c0,1),2*ngw,a(1,1),n1)
  ELSE
     CALL zeroing(a)!,n1*n2)
  ENDIF
  ! >>>>>>> cp_grp trick
  IF (use_cp_grp) THEN
     CALL tiset(procedureN//'_grps_b',isub3)
     ! we reduce so that we get back the cp_grp distribution of the matrix
     CALL cp_grp_redist(a,n1,n2)
     CALL tihalt(procedureN//'_grps_b',isub3)
  ENDIF
  ! <<<<<<<

  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE ovlap2
! ==================================================================
SUBROUTINE ovlap_c(nstate,a,c1,c2)
  ! ==--------------------------------------------------------------==
  ! ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
  ! ==         WITHOUT INVERSION SYMMETRY                           ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE spin , ONLY:spin_mod
  USE cp_grp_utils, ONLY : cp_grp_redist
  USE cp_grp_utils, ONLY : cp_grp_get_sizes
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: a(nstate,nstate), &
                                                c1(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)

  CHARACTER(*), PARAMETER                    :: procedureN = 'ovlap_c'
  COMPLEX(real_8), PARAMETER                 :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

  COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C1_local, C2_local
  INTEGER                                    :: ibeg_c0, ierr, isub, isub2, &
                                                isub3, NGWK_local

! ==--------------------------------------------------------------==

  CALL tiset (procedureN,isub)

  ! >>>>>>> cp_grp trick
  CALL tiset(procedureN//'_grps_a',isub2)
  CALL cp_grp_get_sizes(NGWK_l=NGWK_local,FIRSTK_g=ibeg_c0)
  ! may avoid 1 copy if C1=C2
  ALLOCATE(C1_local(NGWK_local,nstate),&
       C2_local(NGWK_local,nstate),stat=ierr)
  IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
       __LINE__,__FILE__)
  CALL cp_grp_copy_wfn_to_local(c1,nkpt%ngwk,C1_local,NGWK_local,&
       ibeg_c0,NGWK_local,nstate)
  CALL cp_grp_copy_wfn_to_local(c2,nkpt%ngwk,C2_local,NGWK_local,&
       ibeg_c0,NGWK_local,nstate)
  CALL tihalt(procedureN//'_grps_a',isub2)
  ! <<<<<<<

  IF (NGWK_local.GT.0) THEN
     IF (cntl%tlsd) THEN
        CALL zeroing(a)!,nstate*nstate)
        ! ..Alpha spin
        CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,NGWK_local,zone,C1_local(1,1),&
             NGWK_local,C2_local(1,1),ngwk_local,zzero,a(1,1),&
             nstate)
        ! ..Beta spin
        CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,NGWK_local,zone,&
             C1_local(1,spin_mod%nsup+1),NGWK_local,C2_local(1,spin_mod%nsup+1),&
             NGWK_local,zzero,a(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
     ELSE
        CALL zgemm('C','N',nstate,nstate,NGWK_local,zone,&
             C1_local(1,1),NGWK_local,C2_local(1,1),ngwk_local,&
             zzero,a(1,1),nstate)
     ENDIF
  ENDIF

  ! >>>>>>> cp_grp trick
  CALL tiset(procedureN//'_grps_b',isub3)
  ! we reduce so that we get back the cp_grp distribution of the matrix
  CALL cp_grp_redist(a,nstate,nstate)
  DEALLOCATE(C1_local,C2_local,stat=ierr)
  IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
       __LINE__,__FILE__)
  CALL tihalt(procedureN//'_grps_b',isub3)
  ! <<<<<<<

  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE ovlap_c
