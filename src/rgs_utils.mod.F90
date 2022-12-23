MODULE rgs_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_redist
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE state_utils,                     ONLY: zero_wfn
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgs
  PUBLIC :: uinv
  PUBLIC :: rgs_c
  PUBLIC :: uinvc

CONTAINS

  ! ==================================================================
  SUBROUTINE rgs(cp,nstate,smat)
    ! ==--------------------------------------------------------------==
    ! ==  GRAM-SCHMIDT ORTHOGONALIZATION                              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cp(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: smat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rgs'

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: CP_local
    INTEGER                                  :: ibeg_c0, ierr, isub, isub2, &
                                                isub3, NGW_local
    LOGICAL                                  :: GEQ0_local

    IF (nstate.LE.0) RETURN
    CALL tiset(procedureN,isub)

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_a',isub2)
    CALL cp_grp_get_sizes(ngw_l=NGW_local,geq0_l=GEQ0_local,&
         first_g=ibeg_c0)
    ALLOCATE(CP_local(NGW_local,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL cp_grp_copy_wfn_to_local(cp,ncpw%ngw,CP_local,NGW_local,&
         ibeg_c0,NGW_local,nstate)
    CALL tihalt(procedureN//'_grps_a',isub2)
    ! <<<<<<<

    CALL zeroing(smat)!,nstate*nstate)
    IF (NGW_local.GT.0)&
         CALL dsyrk('U','T',nstate,2*NGW_local,2._real_8,CP_local,&
         2*NGW_local,0._real_8,smat,nstate)
    IF (GEQ0_local)&
         CALL dger(nstate,nstate,-1.0_real_8,CP_local,2*NGW_local,&
         CP_local,2*NGW_local,smat,nstate)
    CALL mp_sum(smat,nstate**2,parai%cp_grp)
    CALL uinv('U',smat,nstate,nstate)
    IF (NGW_local.GT.0)&
         CALL dtrmm('R','U','N','N',2*NGW_local,nstate,1._real_8,&
         smat,nstate,CP_local,2*NGW_local)
    IF (GEQ0_local) CALL zclean(CP_local,nstate,NGW_local)

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_b',isub3)
    ! we need to zero the C0 that we can do the reduce
    IF (parai%cp_nogrp.GT.1) CALL zero_wfn(ncpw%ngw,nstate,cp,ncpw%ngw)
    CALL cp_grp_copy_local_to_wfn(CP_local,NGW_local,cp,ncpw%ngw,&
         ibeg_c0,NGW_local,nstate)
    ! we reduce so that we get back the cp_grp distribution of C0
    CALL cp_grp_redist(cp,ncpw%ngw,nstate)
    DEALLOCATE(CP_local,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)
    CALL tihalt(procedureN//'_grps_b',isub3)
    ! <<<<<<<

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgs
  ! ==================================================================
  SUBROUTINE uinv(uplo,smat,lda,n)
    ! ==--------------------------------------------------------------==
    ! == Inversion of a positive definite matrix                      =
    ! == Use for orthogonalisation                                    ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: uplo
    INTEGER                                  :: lda
    REAL(real_8)                             :: smat(lda,*)
    INTEGER                                  :: n

    CHARACTER(*), PARAMETER                  :: procedureN = 'uinv'
    LOGICAL, PARAMETER                       :: get_condition_number = .FALSE.

    INTEGER                                  :: i, info, isub, j
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    REAL(real_8)                             :: anorm, rcond, tmp
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: work

    CALL tiset(procedureN,isub)

    ! debug statements
    IF (get_condition_number) THEN
       anorm=0.0_real_8
       DO i = 1, n
          tmp=0.0_real_8
          DO j = 1, n
             tmp=tmp+ABS(smat(j,i))
          END DO
          anorm=MAX(anorm, tmp)
       END DO
    END IF

    CALL dpotrf(uplo,n,smat,lda,info)

    ! debug statements
    IF (get_condition_number) THEN
       ! Get the condition number
       ALLOCATE(work(3*n), iwork(n))
       CALL dpocon(uplo,n,smat,lda,anorm,rcond,work,iwork,info)
       IF (paral%io_parent) WRITE(*,*)"Condition number in UINV: ",rcond," ANORM: ", anorm
       DEALLOCATE(work, iwork)
    END IF

    IF (info.NE.0) THEN
       IF (info.LT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A)')&
               ' UINV| THE I-TH (',-info,&
               ') ARGUMENT HAD AN ILLEGAL VALUE'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A)')&
               ' UINV| THE LEADING MINOR OF ORDER',info,&
               ' IS NOT POSITIVE DEFINITE, '
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' UINV| AND THE FACTORIZATION COULD NOT BE COMPLETED.'
       ENDIF
       CALL stopgm('UINV','ILLEGAL RESULTS DPOTRF',& 
            __LINE__,__FILE__)
    ENDIF
    CALL dtrtri(uplo,'N',n,smat,lda,info)
    IF (info.NE.0) CALL stopgm('UINV','ILLEGAL RESULTS DTRTRI',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE uinv
  ! ==================================================================
  SUBROUTINE rgs_c(cp,nstate,smat)
    ! ==--------------------------------------------------------------==
    ! ==  GRAM-SCHMIDT ORTHOGONALIZATION                              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cp(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: smat(nstate,nstate)

    INTEGER                                  :: isub

    IF (nstate.LE.0) RETURN
    CALL tiset('     RGS_C',isub)
    CALL zeroing(smat)!,nstate*nstate)
    IF (nkpt%ngwk.GT.0)&
         CALL zherk('U','C',nstate,nkpt%ngwk,1._real_8,cp,nkpt%ngwk,0._real_8,smat,nstate)
    CALL mp_sum(smat,nstate*nstate,parai%allgrp)
    CALL uinvc('U',smat,nstate,nstate)
    IF (ncpw%ngw.GT.0)&
         CALL ztrmm('R','U','N','N',2*ncpw%ngw,nstate,CMPLX(1._real_8,0._real_8,kind=real_8),&
         smat,nstate,cp,2*ncpw%ngw)
    IF (geq0) CALL zclean_k(cp,nstate,ncpw%ngw)
    CALL tihalt('     RGS_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgs_c
  ! ==================================================================
  SUBROUTINE uinvc(uplo,smat,lda,n)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: uplo
    INTEGER                                  :: lda, n
    COMPLEX(real_8)                          :: smat(lda,n)

    INTEGER                                  :: info

    CALL zpotrf(uplo,n,smat,lda,info)
    IF (info.NE.0) THEN
       IF (info.LT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A)')&
               ' UINVC| THE I-TH (',-info,&
               ') ARGUMENT HAD AN ILLEGAL VALUE'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A)')&
               ' UINVC| THE LEADING MINOR OF ORDER',info,&
               ' IS NOT POSITIVE DEFINITE, '
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' UINVC| AND THE FACTORIZATION COULD NOT BE COMPLETED.'
       ENDIF
       CALL stopgm('UINVC','ILLEGAL RESULTS DPOTRF',& 
            __LINE__,__FILE__)
    ENDIF
    CALL ztrtri(uplo,'N',n,smat,lda,info)
    IF (info.NE.0) CALL stopgm('UINVC','ILLEGAL RESULTS ZTRTRI',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE uinvc
  ! ==================================================================


END MODULE rgs_utils
