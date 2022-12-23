MODULE kdp_diag_utils
  USE cnst,                            ONLY: uimag
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: iteropt
  USE system,                          ONLY: cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kdp_diag

CONTAINS

  SUBROUTINE kdp_diag(nstate,pkdp,xlambda0,xkdp,nkdp,ekdp,akdp,&
       noccup,auxdiag,rauxdiag)
    ! 

    INTEGER                                  :: nstate
    REAL(real_8)                             :: pkdp(3,nstate,nstate)
    COMPLEX(real_8)                          :: xlambda0(nstate,nstate)
    INTEGER                                  :: nkdp
    REAL(real_8)                             :: xkdp(3,nkdp), &
                                                ekdp(nstate,nkdp)
    COMPLEX(real_8)                          :: akdp(nstate,nstate,nkdp)
    INTEGER                                  :: noccup
    COMPLEX(real_8)                          :: auxdiag(2*nstate)
    REAL(real_8)                             :: rauxdiag(3*nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'kdp_diag'

    COMPLEX(real_8), ALLOCATABLE             :: akdp2(:,:), hkdp(:), &
                                                xlambda(:,:)
    INTEGER                                  :: icoor, ierr, ikdp, info, &
                                                istat, isub, jstat, kstat
    REAL(real_8)                             :: x3kdp

! 
! 
! ==================================================================
! == Select diagonalization method, according to the relative     ==
! == dimension of the occupied subspce with respect to the full   ==
! == space.                                                       ==
! ==--------------------------------------------------------------==
! temp      if(nint(nel/2._real_8)+5.lt.nint(0.3_real_8*nstate)) then
! temp        noccup=nint(nel/2._real_8)+5
! temp      else
! temp        noccup=nstate
! temp      endif
! ==--------------------------------------------------------------==
! Allocation of memory.

    ALLOCATE(hkdp(nstate*(nstate+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xlambda(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(akdp2(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    noccup=nstate
    ! 
    IF (iteropt%nfi.EQ.cnti%nomore) THEN
       IF (paral%io_parent)&
            WRITE(44,*)'Eigenvalues'
       IF (paral%io_parent)&
            WRITE(44,*)
       IF (paral%io_parent)&
            WRITE(47,*)'Occupations'
       IF (paral%io_parent)&
            WRITE(47,*)
    ENDIF
    ! 
    CALL tiset('  KDP_DIAG',isub)
    ! ==================================================================
    ! == Do-loop over k-points.                                       == 
    ! ==--------------------------------------------------------------==
    DO ikdp=1,nkdp
       ! 
       ! ==================================================================
       ! == Construct k^2.                                               ==
       ! ==--------------------------------------------------------------==
       x3kdp=0._real_8
       DO icoor=1,3
          x3kdp=x3kdp+0.5_real_8*xkdp(icoor,ikdp)**2
       ENDDO
       ! 
       ! 
       ! ==================================================================
       ! == Build k.p matrix at k-point ikdp, and store in array hkdp.   ==
       ! ==--------------------------------------------------------------==
       DO istat=1,nstate
          DO jstat=1,nstate
             xlambda(istat,jstat)=(0._real_8,0._real_8)
             ! 
             ! Nonlocal part (Not tested yet)
             ! isa=0
             ! do is=1,nspnl
             ! do ia=1,na(is)
             ! isa=isa+1
             ! do igh=1,ngh
             ! xlambda(istat,jstat)= xlambda(istat,jstat) + 
             ! &                                  wsg(is,igh) *
             ! &            CONJG(fnlkdp(isa,istat,igh,ikdp))*
             ! &                   fnlkdp(isa,jstat,igh,ikdp)
             ! enddo
             ! enddo
             ! enddo
             ! 
             ! local + kinetic parts
             ! 
             xlambda(istat,jstat) = xlambda(istat,jstat) + xlambda0(&
                  istat,jstat)
          ENDDO
          xlambda(istat,istat) = xlambda(istat,istat) + x3kdp
       ENDDO
       ! 
       ! k.p term
       ! 
       DO jstat=1,nstate
          DO istat=1,nstate
             DO icoor=1,3
                xlambda(istat,jstat)=xlambda(istat,jstat)+uimag*xkdp(&
                     icoor,ikdp)*pkdp(icoor,istat,jstat)
             ENDDO
          ENDDO
       ENDDO
       ! 
       ! deb        print *
       ! deb        write(6,10)xkdp(1,ikdp),xkdp(2,ikdp),xkdp(3,ikdp)
       ! deb        do jstat=1,nstate
       ! deb          do istat=1,jstat
       ! deb            xp=0._real_8
       ! deb            do icoor=1,3
       ! deb              xp=xp+xkdp(icoor,ikdp)*pkdp(icoor,istat,jstat)
       ! deb            enddo
       ! deb            write(6,20)istat,jstat,pkdp(1,istat,jstat),
       ! deb     .                 pkdp(2,istat,jstat),pkdp(3,istat,jstat),xp
       ! deb          enddo
       ! deb        enddo

       ! 
       ! Upper diagonal form of the symmetric matrix xlambda, for fast diag.
       ! 
       kstat=0
       DO jstat=1,nstate
          DO istat=jstat,nstate
             kstat=kstat+1
             hkdp(kstat)=xlambda(istat,jstat)
          ENDDO
       ENDDO
       ! 
       ! ==================================================================
       ! == Diagonalize k.p matrix.                                      ==
       ! ==--------------------------------------------------------------==
       ! ss+jk: Use zhpsv if the number of occupied states is small compared
       ! to the total number of states (9/3/98).
       ! 
       IF (noccup.EQ.nstate) THEN
          CALL zhpev('V','L',nstate,hkdp,ekdp(1,ikdp),akdp2,nstate,&
               auxdiag,rauxdiag,info)
          IF (info.GT.0)CALL stopgm('KDP_DIAG','DIAG FAILED TO CONVERGE'&
               ,& 
               __LINE__,__FILE__)
          IF (info.LT.0)CALL stopgm('KDP_DIAG','WRONG ARG TO ZHPEV',& 
               __LINE__,__FILE__)
          ! ibm          call ZHPEV(1,hkdp,ekdp(1,ikdp),akdp2,nstate,nstate,
          ! ibm     c               auxdiag,11*nstate)
       ELSE
          CALL stopgm('KDP_DIAG','selected eigenvalues not programmed',& 
               __LINE__,__FILE__)
          ! call ZHPEVx('V','I','L',nstate,hkdp,1.,1.,1,noccup,0.,moccup,
          ! c                 ekdp(1,ikdp),akdp2,nstate,auxdiag,3*nstate,rwork
          ! c                 iwork,ifail,info)
          ! ibm          call zhpsv(1,hkdp,ekdp(1,ikdp),akdp2,nstate,nstate,noccup,
          ! ibm     c               auxdiag,11*nstate)
       ENDIF
       ! ==================================================================
       ! == Invert the transformation matrices: akdp.                    ==
       ! ==--------------------------------------------------------------==
       DO istat=1,nstate
          DO jstat=1,nstate
             akdp(istat,jstat,ikdp)=akdp2(jstat,istat)
          ENDDO
       ENDDO
       ! 
       ! deb        WRITE(6,*)ikdp
       ! deb        do istat=1,nstate
       ! deb          WRITE(6,*)istat,ekdp(istat,ikdp)*2._real_8*ry
       ! deb        enddo
       ! deb        print *
       ! 
       ! ==================================================================
       ! == End loop over k-points.                                      == 
       ! ==--------------------------------------------------------------==
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Deallocation of memory.
    DEALLOCATE(hkdp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xlambda,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(akdp2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  KDP_DIAG',isub)
    ! 
    RETURN
  END SUBROUTINE kdp_diag

END MODULE kdp_diag_utils
