MODULE gsortho_utils
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE rgs_utils,                       ONLY: rgs,&
                                             rgs_c
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gs_ortho
  PUBLIC :: gs_ortho_c
  PUBLIC :: gsortho
  PUBLIC :: gsortho_c

CONTAINS

  ! ==================================================================
  SUBROUTINE gs_ortho(c0,n0,cp,np)
    ! ==--------------------------------------------------------------==
    ! == GRAM-SCHMIDT ORTHOGONALISATION                               ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   N0         Number of wavefunctions                         ==
    ! ==   C0(NGW,N0) are already orthogonalized.                     ==
    ! ==   NP         Number of wavefunctions to be orthogonalized    ==
    ! == OUTPUT:                                                      ==
    ! ==   CP(NGW,NP) Wavefunctions to be orthogonalized              ==
    ! ==   SMAT(NP,*) !MAX(N0,NP) Scratch array                       ==
    ! ==   ISMA(NP)   Unused                                          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n0
    COMPLEX(real_8)                          :: c0(ncpw%ngw,n0)
    INTEGER                                  :: np
    COMPLEX(real_8)                          :: cp(ncpw%ngw,np)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gs_ortho'

    INTEGER                                  :: ierr, isub
    REAL(real_8), ALLOCATABLE                :: smat(:,:)

! Variables
! MAX(N0,NP)
! ==--------------------------------------------------------------==

    IF (np.LT.1) RETURN
    ALLOCATE(smat(MAX(np,n0), MAX(np,n0)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(smat)
    IF (2*spar%ngws-1.LT.(n0+np)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I6,/,A,I6)')&
            ' GS_ORTHO! The number of components (2*NGWS-1)',&
            2*spar%ngws-1,&
            ' GS_ORTHO! is not greater than the number of states',n0+nP
       CALL stopgm('GS_ORTHO','TOO MANY STATES',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset('   GSORTHO',isub)
    IF (n0.GT.0) THEN
       CALL ovlap2(ncpw%ngw,np,n0,smat,cp,c0,.TRUE.)
       CALL mp_sum(smat,np*n0,parai%allgrp)
       CALL dgemm('N','T',2*ncpw%ngw,np,n0,-1._real_8,c0,2*ncpw%ngw,smat,np,&
            1._real_8,cp,2*ncpw%ngw)
    ENDIF
    CALL rgs(cp,np,smat)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('   GSORTHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gs_ortho
  ! ==================================================================
  SUBROUTINE gs_ortho_c(c0,n0,cp,np,smat)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*)
    INTEGER                                  :: n0
    COMPLEX(real_8)                          :: cp(nkpt%ngwk,*)
    INTEGER                                  :: np
    COMPLEX(real_8)                          :: smat(np,*)

    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) 

    INTEGER                                  :: isub

    IF (np.LT.1) RETURN
    IF (spar%ngwks.LT.(n0+np)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I6,/,A,I6)')&
            ' GS_ORTHO_C! The number of components NGWKS',spar%ngwks,&
            ' GS_ORTHO_C! is not greater than the number of states',&
            n0+np
       CALL stopgm('GS_ORTHO_C','TOO MANY STATES',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset(' GSORTHO_C',isub)
    IF (n0.GT.0) THEN
       CALL ovlap2_c(nkpt%ngwk,np,n0,smat,cp,c0)
       CALL mp_sum(smat,2*np*n0,parai%allgrp)
       CALL zgemm('N','C',nkpt%ngwk,np,n0,-zone,c0,nkpt%ngwk,smat,np,&
            zone,cp,nkpt%ngwk)
    ENDIF
    CALL rgs_c(cp,np,smat)
    CALL tihalt(' GSORTHO_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gs_ortho_c
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE gsortho(cp,n,nfirst,nlast)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: cp(n,*)
    INTEGER                                  :: nfirst, nlast

    INTEGER                                  :: i, isub, j
    REAL(real_8)                             :: anorm, s

    CALL tiset('   GSORTHO',isub)
    DO i=nfirst,nlast
       DO j=1,i-1
          s=-dotp(n,cp(:,j),cp(:,i))
          CALL mp_sum(s,parai%allgrp)
          CALL daxpy(2*n,s,cp(1,j),1,cp(1,i),1)
       ENDDO
       ! Because c0(-g)^*=c0(g) (Important).
       IF (geq0) cp(1,i)=CMPLX(REAL(cp(1,i)),0._real_8,kind=real_8)
       anorm=dotp(n,cp(:,i),cp(:,i))
       CALL mp_sum(anorm,parai%allgrp)
       anorm=1.0_real_8/SQRT(anorm)
       CALL dscal(2*n,anorm,cp(1,i),1)
    ENDDO
    CALL tihalt('   GSORTHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gsortho
  ! ==================================================================
  SUBROUTINE gsortho_c(cp,n,nfirst,nlast)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: cp(n,*)
    INTEGER                                  :: nfirst, nlast

    COMPLEX(real_8)                          :: zs
    COMPLEX(real_8), EXTERNAL                :: zdotc
    INTEGER                                  :: i, isub, j, n2
    REAL(real_8)                             :: anorm
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==

    CALL tiset(' GSORTHO_C',isub)
    n2=2*n
    DO i=nfirst,nlast
       DO j=1,i-1
          zs=-zdotc(n,cp(1,j),1,cp(1,i),1)
          CALL mp_sum(zs,parai%allgrp)
          CALL zaxpy(n,zs,cp(1,j),1,cp(1,i),1)
       ENDDO
       ! ..Convention
       IF (geq0) cp(n/2+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
       anorm=ddot(n2,cp(1,i),1,cp(1,i),1)
       CALL mp_sum(anorm,parai%allgrp)
       anorm=1.0_real_8/SQRT(anorm)
       CALL dscal(n2,anorm,cp(1,i),1)
    ENDDO
    CALL tihalt(' GSORTHO_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gsortho_c
  ! ==================================================================

END MODULE gsortho_utils
