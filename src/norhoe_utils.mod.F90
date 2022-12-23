MODULE norhoe_utils
  USE csmat_utils,                     ONLY: csmat
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nlps,                            ONLY: imagp
  USE noforce_utils,                   ONLY: give_scr_noforce
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rotate_utils,                    ONLY: rotate
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: norhoe

CONTAINS

  ! ==================================================================
  SUBROUTINE norhoe(c0,sc0,eigv,rhoe,psi,xmat1,xmat2,nstate)
    ! ==--------------------------------------------------------------==
    ! This routine estimate the electronic density for a set of 
    ! non orthogonal orbitals
    ! NNN (Sep.23,2005)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eigv(*), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: xmat2(nstate,*), &
                                                xmat1(nstate,*)
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'norhoe'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, ierr, ij, il_auxc, &
                                                il_ddia, il_gam, il_smat, &
                                                isub, j, lnorho
    INTEGER, SAVE                            :: iconv = 0
    LOGICAL                                  :: debug
    REAL(real_8), ALLOCATABLE                :: auxc(:), gam(:,:), smat(:)

!(nnr1,clsd%nlsd)
!(maxfftn)
! Variables
! 

    CALL tiset('    NORHOE',isub)
    debug=.FALSE.
    IF (lspin2%tlse) CALL stopgm('NORHOE','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2) CALL stopgm('NORHOE','K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! 
    CALL give_scr_noforce(lnorho,il_gam,il_auxc,il_smat,il_ddia,tag,&
         nstate,.FALSE.)
    ! 
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ! 
    ! Compute the FNL array in non orthogonal basis
    IF (pslo_com%tivan) CALL rnlsm(c0,nstate,1,1,.FALSE.)
    ! 
    ALLOCATE(gam(nstate, il_gam/nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(auxc(il_auxc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(smat(il_smat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    ! Compute the overlap matrix
    CALL csmat(gam,c0,fnl,nstate,1)
    ! 
    ! Compute the S**(-1/2) array for doing orthogonalization of C0
    IF (paral%parent) THEN
       CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,gam(1,1),nstate,&
            xmat1(1,1),nstate,0.0_real_8,xmat2(1,1),nstate)
       CALL dgemm('T','N',nstate,nstate,nstate,1.0_real_8,xmat1(1,1),nstate,&
            xmat2(1,1),nstate,0.0_real_8,gam(1,1),nstate)
       CALL jacobi(nstate,nstate,gam,eigv,xmat2,ierr)
       CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,xmat1(1,1),nstate,&
            xmat2(1,1),nstate,0.0_real_8,smat(1),nstate)
       CALL dcopy(nstate*nstate,smat(1),1,xmat1(1,1),1)
       iconv=iconv+1
       DO i=1,nstate
          eigv(i)=1.0_real_8/SQRT(eigv(i))
       ENDDO
       ij=0
       DO i=1,nstate
          DO j=1,nstate
             ij=ij+1
             auxc(ij)=eigv(i)*smat(ij)
          ENDDO
       ENDDO
       CALL dgemm('N','T',nstate,nstate,nstate,1.0_real_8,auxc(1),nstate,&
            smat(1),nstate,0.0_real_8,gam(1,1),nstate)
       CALL dcopy(nstate*nstate,gam(1,1),1,smat(1),1)
       CALL dcopy(nstate*nstate,smat(1),1,xmat2(1,1),1)
    ENDIF
    CALL mp_bcast(smat,SIZE(smat),parai%source,parai%allgrp)
    ! 
    ! Do symmetric orthogonalization. SC0 is having the orthogonal basis
    CALL rotate(1.0_real_8,c0,0.0_real_8,sc0,smat,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ! 
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! 
    ! Compute the FNL array for the orthogonal basis
    CALL rnlsm(sc0,nstate,1,1,.FALSE.)
    ! 
    ! Compute the electronic density for the orthogonal basis
    CALL rhoofr(sc0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    NORHOE',isub)
    RETURN
  END SUBROUTINE norhoe
  ! ==================================================================

END MODULE norhoe_utils
