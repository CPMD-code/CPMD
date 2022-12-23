MODULE rortog_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE reigs_utils,                     ONLY: reigs
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             parap,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE utils,                           ONLY: symma
!!use ovlap_utils, only : ovlap2
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rortog
  !public :: rortho
  !public :: rortho_da
  PUBLIC :: give_scr_rortog

CONTAINS

  ! ==================================================================
  SUBROUTINE rortog(c0,cm,gam,tnon,nstate,prteig)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(*), cm(*)
    REAL(real_8)                             :: gam(*)
    LOGICAL                                  :: tnon
    INTEGER                                  :: nstate
    LOGICAL                                  :: prteig

    INTEGER                                  :: isub, nstx
    REAL(real_8)                             :: pfac

! Variables
! ==--------------------------------------------------------------==
! ==  Compute the eigenvalues needed by the RORTHO orthogonali-   ==
! ==  zation routine                                              ==
! ==                   ***** WARNING *****                        ==
! ==  If running the code with the rescaling of electronic        ==
! ==  velocities option (cntl%tc), and orthogonalizing with rortho,    ==
! ==  then the eigenvalues printed in routine EIGS have no        ==
! ==  meaning unless convergence is very close (small ekinc).     ==
! ==--------------------------------------------------------------==

    CALL tiset('    RORTOG',isub)
    IF (cntl%tdmal) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       CALL rortho_da(c0,cm,gam,tnon,nstate)
       IF (prteig.AND.paral%parent) THEN
          WRITE(6,*) " WARNING: NO EIGENVALUE PRINTING ",&
               "WITH DISTRIBUTED LINALG"
       ENDIF
    ELSE
       CALL rortho(c0,cm,gam,tnon,nstate)
       IF (prteig.AND.paral%parent) THEN
          pfac=2.0_real_8/dt2bye
          CALL dscal(nstate*nstate,pfac,gam(1),1)
          CALL reigs(nstate,gam,crge%f)
          pfac=1.0_real_8/pfac
          CALL dscal(nstate*nstate,pfac,gam(1),1)
       ENDIF
    ENDIF
    CALL tihalt('    RORTOG',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rortog
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE rortho_da(c0,cm,gam,tnon,nstate)
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==       C0 (ORTHONORMAL)                                       ==
    ! ==       CM (NON-ORTHONORMAL)                                   ==
    ! == OUTPUT:                                                      ==
    ! ==       GAM                                                    ==
    ! ==       CM (ORTHONORMAL)                                       ==
    ! == METHOD USED:                                                 ==
    ! ==                                                              ==
    ! == ITERATIVE ALGORITHM PARRINELLO-CAR LES HOUCHES 1988          ==
    ! == EQS. 18-22. IMPLEMENTATION ZURICH 1989, A.NOBILE             ==
    ! ==                                                              ==
    ! == CODE FOR HARMONIC REFERENCE SYSTEM ADDED                     ==
    ! ==                                                              ==
    ! == SPECIAL VERSION FOR DISTRIBUTED MATRICES                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: gam(*)
    LOGICAL                                  :: tnon
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rortho_da'

    INTEGER                                  :: i, i1, ierr, ii, ip, ipp, &
                                                iter, msglen, nstx, nx
    INTEGER, EXTERNAL                        :: idamax
    REAL(real_8)                             :: difgam
    REAL(real_8), ALLOCATABLE                :: del(:,:), gamn(:), rho(:,:), &
                                                scr(:), sig(:,:), sxx(:)

    CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
    ALLOCATE(gamn(nstate*nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO ALIGN FOR BG
    ALLOCATE(scr(nstate*nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(sig(nstate,nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rho(nstate,nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(del(nstate,nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(sxx(nstate*nstx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ..INITIALIZE RHO, SIG, AND DEL
    CALL zeroing(rho)!,nstate*nstx)
    CALL zeroing(sig)!,nstate*nstx)
    CALL zeroing(del)!,nstate*nstx)
    msglen=nstate*nstx * 8
    DO ip=0,parai%nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       CALL zeroing(gam(1:nstate*nstx))!,nstate*nstx)
       IF (nx.GT.0) THEN
          CALL ovlap2(ncpw%ngw,nstate,nx,gam,c0,cm(1,paraw%nwa12(ip,1)),.TRUE.)
          CALL mp_sum(gam,scr,nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
          IF (parai%mepos.EQ.parap%pgroup(ip+1)) CALL dcopy(nstate*nstx,scr,1,rho,1)
       ENDIF
    ENDDO
    ! 
    DO ip=0,parai%nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       CALL zeroing(gam(1:nstate*nstx))!,nstate*nstx)
       IF (nx.GT.0) THEN
          CALL ovlap2(ncpw%ngw,nstate,nx,gam,cm,cm(1,paraw%nwa12(ip,1)),.TRUE.)
          CALL mp_sum(gam,scr,nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
          IF (parai%mepos.EQ.parap%pgroup(ip+1)) CALL dcopy(nstate*nstx,scr,1,sig,1)
       ENDIF
    ENDDO
    IF (tnon) THEN
       DO ip=0,parai%nproc-1
          nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
          CALL zeroing(gam(1:nstate*nstx))!,nstate*nstx)
          IF (nx.GT.0) THEN
             CALL ovlap2(ncpw%ngw,nstate,nx,gam,c0,c0(1,paraw%nwa12(ip,1)),.TRUE.)
             CALL mp_sum(gam,scr,nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
             IF (parai%mepos.EQ.parap%pgroup(ip+1))&
                  CALL dcopy(nstate*nstx,scr,1,del,1)
          ENDIF
       ENDDO
    ENDIF
    CALL dscal(nstate*nstx,-1.0_real_8,rho,1)
    CALL dscal(nstate*nstx,-1.0_real_8,sig,1)
    DO i=paraw%nwa12(parai%mepos,1),paraw%nwa12(parai%mepos,2)
       ii=i-paraw%nwa12(parai%mepos,1)+1
       rho(i,ii)=1.0_real_8+rho(i,ii)
       sig(i,ii)=1.0_real_8+sig(i,ii)
    ENDDO
    IF (tnon) THEN
       DO i=paraw%nwa12(parai%mepos,1),paraw%nwa12(parai%mepos,2)
          ii=i-paraw%nwa12(parai%mepos,1)+1
          del(i,ii)=del(i,ii)-1.0_real_8
       ENDDO
       CALL dscal(nstate*nstx,-0.5_real_8,del,1)
    ENDIF
    ! ..INITIALIZE GAM
    CALL dcopy(nstate*nstx,sig,1,gam,1)
    CALL dscal(nstate*nstx,0.5_real_8,gam,1)
    ! ..INITIALIZE FIXED TERM OF RECURRENCE
    CALL dcopy(nstate*nstx,rho,1,scr,1)
    CALL dcopy(nstate*nstx,rho,1,sxx,1)
    DO ip=0,parai%nproc-1
       i1=paraw%nwa12(ip,1)
       crge%n=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       IF (crge%n.GT.0) THEN
          CALL dcopy(nstate*nstx,scr,1,rho,1)
          CALL mp_bcast(rho,nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
          CALL dgemm('T','N',crge%n,nstx,nstate,0.5_real_8,rho,nstate,&
               sxx,nstate,0.5_real_8,sig(i1,1),nstate)
       ENDIF
    ENDDO
    iter=0
    CALL dcopy(nstate*nstx,scr,1,rho,1)
100 CONTINUE
    iter=iter+1
    CALL dcopy(nstate*nstx,sig,1,gamn,1)
    DO i=1,nstate*nstx
       scr(i)=rho(i,1)-gam(i)
    ENDDO
    CALL dcopy(nstate*nstx,scr,1,sxx,1)
    CALL dcopy(nstate*nstx,scr,1,del,1)
    DO ip=0,parai%nproc-1
       i1=paraw%nwa12(ip,1)
       crge%n=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       IF (crge%n.GT.0) THEN
          CALL dcopy(nstate*nstx,sxx,1,scr,1)
          CALL mp_bcast(scr,nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
          CALL dgemm('T','N',crge%n,nstx,nstate,-0.5_real_8,scr,nstate,&
               del,nstate,1._real_8,gamn(i1),nstate)
       ENDIF
    ENDDO
    IF (tnon) THEN
       CALL zeroing(scr)!,nstate*nstx)
       DO ip=0,parai%nproc-1
          CALL my_shift(del,sxx,msglen,parai%mepos,-ip,parai%allgrp)
          ipp=MOD(parai%mepos+ip,parai%nproc)
          crge%n=paraw%nwa12(ipp,2)-paraw%nwa12(ipp,1)+1
          i1=paraw%nwa12(ipp,1)
          CALL dgemm('N','N',nstate,nstx,crge%n,1._real_8,sxx,nstate,&
               gam(i1),nstate,1._real_8,scr,nstate)
       ENDDO
       DO ip=0,parai%nproc-1
          CALL my_shift(gam,sxx,msglen,parai%mepos,-ip,parai%allgrp)
          ipp=MOD(parai%mepos+ip,parai%nproc)
          crge%n=paraw%nwa12(ipp,2)-paraw%nwa12(ipp,1)+1
          i1=paraw%nwa12(ipp,1)
          CALL dgemm('T','N',crge%n,nstx,nstate,1._real_8,sxx,nstate,&
               scr,nstate,1._real_8,gamn(i1),nstate)
       ENDDO
    ENDIF
    DO i=1,nstate*nstx
       scr(i)=gamn(i)-gam(i)
    ENDDO
    CALL dcopy(nstate*nstx,gamn,1,gam,1)
    difgam=ABS(scr(idamax(nstate*nstx,scr,1)))
    CALL mp_max(difgam,parai%allgrp)
    IF (difgam.GT.cntr%epsog*0.5_real_8.AND.iter.LE.cnti%maxit)GOTO 100
    ! 
    IF (iter.GT.cnti%maxit.AND.paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'RORTHO:DIFGAM=',difgam,' ITER = ',iter,' > ',cnti%maxit
       IF (paral%io_parent)&
            WRITE(6,*) '       MAXIMUM NUMBER OF ITERATIONS EXCEEDED'
       CALL stopgm('RORTOG','NOT CONVERGED',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..UPDATE WAVEFUNCTIONS
    CALL rotate_da(1._real_8,c0,1._real_8,cm,gam,2*ncpw%ngw,2*ncpw%ngw,nstate,&
         paraw%nwa12(0,1),paraw%nwa12(0,2),nstx,parai%mepos,parap%pgroup,parai%nproc,parai%allgrp,&
         cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(gamn,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(sig,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(del,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(sxx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE rortho_da
  ! ==================================================================
  SUBROUTINE give_scr_rortog(lrortog,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrortog
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: nstx

! ==--------------------------------------------------------------==

    IF (cntl%tdmal) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       lrortog=MAX(6*nstate*nstx,nstate*(4+nstx))
       tag   ='MAX(6*NSTATE*NSTX,...)'
    ELSE
       lrortog=MAX(5*nstate*nstate,nstate*(4+nstate))
       tag   ='MAX(5*NSTATE*NSTATE,...)'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rortog
  ! ==================================================================

END MODULE rortog_utils



SUBROUTINE rortho(c0,cm,gam,tnon,nstate)
  ! ==--------------------------------------------------------------==
  ! == INPUT:                                                       ==
  ! ==       C0 (ORTHONORMAL)                                       ==
  ! ==       CM (NON-ORTHONORMAL)                                   ==
  ! == OUTPUT:                                                      ==
  ! ==       GAM                                                    ==
  ! ==       CM (ORTHONORMAL)                                       ==
  ! == METHOD USED:                                                 ==
  ! ==                                                              ==
  ! == ITERATIVE ALGORITHM PARRINELLO-CAR LES HOUCHES 1988          ==
  ! == EQS. 18-22. IMPLEMENTATION ZURICH 1989, A.NOBILE             ==
  ! ==                                                              ==
  ! == WITH RESPECT TO PREVIOUS VERSIONS OF THIS ROUTINE, THIS ONE: ==
  ! == -)FOLLOWS EXACTLY THE REFERENCE;                             ==
  ! == 1)SAVES ONE EVALUATION OF SIGMA BY NOTICING THAT             ==
  ! == THE VALUE OF SIGMA FROM THE UPDATED WAVEFUNCTION IS          ==
  ! == EQUAL TO TWICE THE NEXT VALUE OF GAMN - GAM,                 ==
  ! == WHICH IS EXPECTED TO BE LESS THAN THE CURRENT ONE            ==
  ! == THEREFORE, CONVERGENCE OF GAM GUARANTEES CONVERGENCE OF SIG  ==
  ! == 2)SAVES ONE MATRIX PRODUCT PER ITERATION BECAUSE             ==
  ! == SIG+TRANS(RHO)*G+G*RHO-G*G = SIG + TRANS(RHO)*RHO -          ==
  ! ==                                  - TRANS(RHO-G)*(RHO-G)      ==
  ! ==                                                              ==
  ! == TO FOLLOW CLOSELY THE REFERENCE, THE INPUT VALUE OF GAM IS   ==
  ! == IGNORED, AND AS AN INITIAL VALUE WE USE 0.5*SIG              ==
  ! == THE ALGORITHM, LIKE ITS PREDECESSORS, IGNORES THE REAL       ==
  ! == PARTS OF EVERYTHING (EXCEPT OF THE WAVEFUNCTIONS)            ==
  ! ==                                                              ==
  ! == CODE FOR HARMONIC REFERENCE SYSTEM ADDED                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE system , ONLY:cnti,cntl,cntr,ncpw
  USE parac, ONLY : paral,parai
  USE spin , ONLY:spin_mod
  USE ovlap_utils, ONLY : ovlap
  USE summat_utils, ONLY : summat_parent
  USE utils, ONLY : symma
  USE rotate_utils, ONLY: rotate
  IMPLICIT NONE
  REAL(real_8)                               :: gam(*)
  LOGICAL                                    :: tnon
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'rortho'

  INTEGER                                    :: i, idamax, ierr, ii, iter, j
  REAL(real_8)                               :: difgam
  REAL(real_8), ALLOCATABLE                  :: del(:,:), gamn(:), rho(:,:), &
                                                scr(:), sig(:,:)

  ALLOCATE(gamn(nstate*nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(scr(nstate*nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(sig(nstate, nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(rho(nstate, nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(del(nstate, nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  ! == INITIALIZE RHO                                               ==
  ! ==--------------------------------------------------------------==
  CALL ovlap(nstate,rho,c0,cm)
  CALL symma(rho,nstate)
  CALL summat_parent(rho,nstate)
  CALL dscal(nstate*nstate,-1.0_real_8,rho(1,1),1)
  !$omp parallel do private(I)
  DO i=1,nstate
     rho(i,i)=1.0_real_8+rho(i,i)
  ENDDO
  ! ==--------------------------------------------------------------==
  ! == INITIALIZE SIG                                               ==
  ! ==--------------------------------------------------------------==
  CALL ovlap(nstate,sig,cm,cm)
  CALL summat_parent(sig,nstate)
  CALL dscal(nstate*nstate,-1.0_real_8,sig(1,1),1)
  !$omp parallel do private(I)
  DO i=1,nstate
     sig(i,i)=1.0_real_8+sig(i,i)
  ENDDO
  ! ==--------------------------------------------------------------==
  ! == INITIALIZE DEL                                               ==
  ! ==--------------------------------------------------------------==
  IF (tnon) THEN
     CALL ovlap(nstate,del,c0,c0)
     CALL summat_parent(del,nstate)
     !$omp parallel do private(I)
     DO i=1,nstate
        del(i,i)=del(i,i)-1.0_real_8
     ENDDO
     CALL dscal(nstate*nstate,-0.5_real_8,del(1,1),1)
  ENDIF
  IF (paral%parent) THEN
     ! ==--------------------------------------------------------------==
     ! == INITIALIZE GAM                                               ==
     ! ==--------------------------------------------------------------==
     !$omp parallel do private(I,J,II)
     DO j=1,nstate
        DO i=1,nstate
           ii=i+(j-1)*nstate
           gam(ii)=sig(i,j)*0.5_real_8
        ENDDO
     ENDDO
     ! ==--------------------------------------------------------------==
     ! == INITIALIZE FIXED TERM OF RECURRENCE                          ==
     ! ==--------------------------------------------------------------==
     CALL dgemm('T','N',nstate,nstate,nstate,0.5_real_8,rho(1,1),nstate,&
          rho(1,1),nstate,0.5_real_8,sig(1,1),nstate)
     iter=0
     ! ==--------------------------------------------------------------==
     ! == ITERATIVE LOOP                                               ==
     ! ==--------------------------------------------------------------==
100  CONTINUE
     iter=iter+1
     CALL dcopy(nstate*nstate,sig(1,1),1,gamn(1),1)
     !$omp parallel do private(I,J,II)
     DO j=1,nstate
        DO i=1,nstate
           ii=i+(j-1)*nstate
           scr(ii)=rho(i,j)-gam(ii)
        ENDDO
     ENDDO
     CALL dgemm('T','N',nstate,nstate,nstate,-0.5_real_8,&
          scr,nstate,scr,nstate,1.0_real_8,gamn(1),nstate)
     IF (tnon) THEN
        CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,&
             del(1,1),nstate,gam(1),nstate,0.0_real_8,scr(1),nstate)
        CALL dgemm('T','N',nstate,nstate,nstate,1.0_real_8,&
             gam(1),nstate,scr(1),nstate,1.0_real_8,gamn(1),nstate)
     ENDIF
     !$omp parallel do private(I)
     DO i=1,nstate*nstate
        scr(i)=gamn(i)-gam(i)
     ENDDO
     CALL dcopy(nstate*nstate,gamn(1),1,gam(1),1)
     difgam=ABS(scr(idamax(nstate*nstate,scr,1)))
     IF (difgam.GT.cntr%epsog*0.5_real_8.AND.iter.LE.cnti%maxit)GOTO 100
     IF (iter.GT.cnti%maxit)THEN
        IF (paral%io_parent)&
             WRITE(6,*) 'RORTHO:DIFGAM=',difgam,' ITER = ',iter,' > ',cnti%maxit
        IF (paral%io_parent)&
             WRITE(6,*) '       MAXIMUM NUMBER OF ITERATIONS EXCEEDED'
        CALL stopgm('RORTOG',' ',& 
             __LINE__,__FILE__)
     ENDIF
     CALL symma(gam,nstate)
  ENDIF
  CALL mp_bcast(gam,nstate*nstate,parai%source,parai%allgrp)
  ! ==--------------------------------------------------------------==
  ! == UPDATE WAVEFUNCTIONS                                         ==
  ! ==--------------------------------------------------------------==
  CALL rotate(1.0_real_8,c0,1.0_real_8,cm,gam,nstate,2*ncpw%ngw,cntl%tlsd,&
       spin_mod%nsup,spin_mod%nsdown)
  ! ==--------------------------------------------------------------==
  DEALLOCATE(gamn,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(scr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(sig,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(rho,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  DEALLOCATE(del,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  RETURN
END SUBROUTINE rortho
