MODULE rhodiis_utils
  USE andr,                            ONLY: andr3
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE odiis_utils,                     ONLY: solve
  USE parac,                           ONLY: parai
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhodiis
  PUBLIC :: change_nrdiis

CONTAINS

  ! ==================================================================
  SUBROUTINE rhodiis(np,rnew,rin,rout,riter,driter,&
       amix,dmat,dm,vb,nrdiismax,nrdiis,ilsd,thl)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE NEW DENSITY USING cntl%diis SCHEME                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: np
    REAL(real_8)                             :: rnew(np), rin(np), rout(np), &
                                                amix
    INTEGER                                  :: nrdiismax
    REAL(real_8)                             :: vb(nrdiismax+1), &
                                                dm(nrdiismax+1,nrdiismax+1), &
                                                dmat(nrdiismax,nrdiismax)
    INTEGER                                  :: nrdiis
    REAL(real_8)                             :: driter(np,nrdiis), &
                                                riter(np,nrdiis)
    INTEGER                                  :: ilsd
    REAL(real_8)                             :: thl

    INTEGER                                  :: i, j, niter, nsize
    INTEGER, DIMENSION(2), SAVE              :: ncall = (/0,0/)
    REAL(real_8)                             :: fac
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
! We can reinitialize (for cntl%md).

    IF (nrdiis.EQ.0.OR.np.EQ.0) THEN
       ncall(ilsd) = 0
       RETURN
    ENDIF
    ncall(ilsd)=ncall(ilsd)+1
    IF (ncall(ilsd).LT.nrdiis) THEN
       nsize=ncall(ilsd)
    ELSE
       nsize=nrdiis
    ENDIF
    niter=MOD(ncall(ilsd)-1,nrdiis)+1
    ! Update buffer
    CALL dcopy(np, rin(1),1, riter(1,niter),1)
    CALL dcopy(np,rout(1),1,driter(1,niter),1)
    CALL daxpy(np,-1._real_8,rin(1),1,driter(1,niter),1)
    ! Update cntl%diis Matrix
    DO i=1,nsize
       dmat(i,niter)=ddot(np,driter(1,i),1,driter(1,niter),1)
    ENDDO
    CALL mp_sum(dmat(:,niter),nsize,parai%allgrp)
    DO i=1,nsize
       IF (i.NE.niter) dmat(niter,i)=dmat(i,niter)
    ENDDO
    ! Setup linear system
    CALL zeroing(dm)!,(nrdiismax+1)*(nrdiismax+1))
    CALL zeroing(vb)!,nrdiismax+1)
    DO i=1,nsize
       DO j=1,nsize
          dm(i,j)=dmat(i,j)
       ENDDO
       dm(i,nsize+1)=-1._real_8
       dm(nsize+1,i)=-1._real_8
    ENDDO
    dm(nsize+1,nsize+1)=0._real_8
    vb(nsize+1)=-1._real_8
    ! Solve linear system
    CALL solve(dm,nrdiismax+1,nsize+1,vb)
    ! Estimate new vectors
    CALL zeroing(rnew)!,np)
    fac=1.0_real_8
    IF (nsize.EQ.1) fac=amix
    DO i=1,nsize
       CALL daxpy(np,    vb(i), riter(1,i),1,rnew(1),1)
       CALL daxpy(np,fac*vb(i),driter(1,i),1,rnew(1),1)
    ENDDO
    thl=dmat(niter,niter)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhodiis
  ! ==================================================================
  SUBROUTINE change_nrdiis(drhomax)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: drhomax

    INTEGER                                  :: i
    INTEGER, SAVE                            :: ifirst = 0

! ==--------------------------------------------------------------==

    IF (andr3%ntabnrdiis.LE.1) RETURN
    IF (ifirst.EQ.0) THEN
       andr3%inrdiis=1
       ifirst=1
    ENDIF
    andr3%nrdiis=andr3%tabnrdiis(andr3%inrdiis)
    DO i=andr3%inrdiis,andr3%ntabnrdiis
       IF (drhomax.LT.andr3%densnrdiis(i)) THEN
          andr3%nrdiis=andr3%tabnrdiis(i)
          andr3%inrdiis=i
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE change_nrdiis
  ! ==================================================================

END MODULE rhodiis_utils
