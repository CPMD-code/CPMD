MODULE dipo_utils
  USE bessm_utils,                     ONLY: bessl
  USE cnst,                            ONLY: pi,&
                                             uimag
  USE cppt,                            ONLY: gk,&
                                             gl,&
                                             isptr,&
                                             nzh
  USE dipomod,                         ONLY: moment
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dipo
  PUBLIC :: rsdipo

CONTAINS

  ! ==================================================================
  SUBROUTINE dipo(tau0,eirop,v)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE DIPOLE MOMENT                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), v(*)

    INTEGER                                  :: i, ia, ig, ig1, is, ish, isub
    REAL(real_8)                             :: b2, dfac, dx, dy, dz, rg, &
                                                rho, rmax, rrx, rry, rrz, &
                                                xkr, zz

    CALL tiset('      DIPO',isub)
    ! Define the center of integration
    rrx=0.0_real_8
    rry=0.0_real_8
    rrz=0.0_real_8
    zz=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          rrx=rrx+tau0(1,ia,is)*ions0%zv(is)
          rry=rry+tau0(2,ia,is)*ions0%zv(is)
          rrz=rrz+tau0(3,ia,is)*ions0%zv(is)
          zz=zz+ions0%zv(is)
       ENDDO
    ENDDO
    rrx=rrx/zz
    rry=rry/zz
    rrz=rrz/zz
    ! EHR[
    IF (cntl%tmdeh) THEN
       ! compute dipole of the electron density alone"
       rrx=0._real_8
       rry=0._real_8
       rrz=0._real_8
    ENDIF
    IF (paral%io_parent.AND.(.NOT.cntl%cmplx_wf))&
         WRITE(6,'(/,A,3F9.5,/)')&
         ' CENTER OF INTEGRATION (CORE CHARGE): ', RRX, RRY, RRZ
    ! EHR]
    ! Calculate the radius of the integrated sphere
    IF (parm%ibrav.EQ.1) THEN
       rmax=0.5_real_8*parm%alat
    ELSE IF (parm%ibrav.EQ.2) THEN
       rmax=parm%alat*0.5_real_8/SQRT(2._real_8)
    ELSE IF (parm%ibrav.EQ.8) THEN
       rmax=0.5_real_8*parm%omega**(1.0_real_8/3.0_real_8)
    ELSE IF (parm%ibrav.GT.2) THEN
       CALL stopgm('DIPO',' FIND THE INTEGRATION RADIUS YOURSELF',& 
            __LINE__,__FILE__)
    ENDIF
    DO i=1,3
       moment%dmom(i)=0.0_real_8
       moment%dnuc(i)=0.0_real_8
    ENDDO
    IF (geq0) THEN
       ig1=2
    ELSE
       ig1=1
    ENDIF
    DO ish=ig1,ncpw%nhgl
       dx=0.0_real_8
       dy=0.0_real_8
       dz=0.0_real_8
       DO ig=isptr(ish),isptr(ish+1)-1
          rg=(gk(1,ig)*rrx+gk(2,ig)*rry+gk(3,ig)*rrz)*parm%tpiba
          rho=AIMAG((v(nzh(ig))+eirop(ig))*EXP(uimag*rg))
          dx=dx+rho*gk(1,ig)
          dy=dy+rho*gk(2,ig)
          dz=dz+rho*gk(3,ig)
       ENDDO
       xkr=rmax*parm%tpiba*SQRT(gl(ish))
       b2=bessl(2,xkr)
       moment%dmom(1)=moment%dmom(1)+b2*dx/gl(ish)
       moment%dmom(2)=moment%dmom(2)+b2*dy/gl(ish)
       moment%dmom(3)=moment%dmom(3)+b2*dz/gl(ish)
    ENDDO
    dfac=8.0_real_8*pi/parm%tpiba*rmax**3
    moment%dmom(1)=dfac*moment%dmom(1)
    moment%dmom(2)=dfac*moment%dmom(2)
    moment%dmom(3)=dfac*moment%dmom(3)
    CALL mp_sum(moment%dmom,3,parai%allgrp)
    CALL tihalt('      DIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dipo
  ! ==================================================================
  SUBROUTINE rsdipo(tau0,eirop,psi,rhor)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE DIPOLE MOMENT                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), psi(maxfft)
    REAL(real_8) :: rhor(fpar%kr1,fpar%kr2,fpar%kr3)

    INTEGER                                  :: i, ia, is, isub, ix, iy, iz, &
                                                nx0
    REAL(real_8)                             :: dx(3,3), rho, rr(3), scl, &
                                                x(3), zz

    CALL tiset('    RSDIPO',isub)
    CALL setfftn(0)
    ! Define the center of integration
    DO i=1,3
       rr(i)=0.0_real_8
    ENDDO
    zz=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO i=1,3
             rr(i)=rr(i)+tau0(i,ia,is)*ions0%zv(is)
          ENDDO
          zz=zz+ions0%zv(is)
       ENDDO
    ENDDO
    DO i=1,3
       rr(i)=rr(i)/zz
    ENDDO
    IF (paral%io_parent)&
         WRITE (6,'(/,A,3F9.5,/)')&
         ' CENTER OF INTEGRATION (CORE CHARGE): ', (RR(I), I=1,3)
    CALL  invfftn(psi,.FALSE.,parai%allgrp)
    ! Calculate the dipole moment in real space
    scl=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    DO i=1,3
       moment%dmom(i)=0.0_real_8
       moment%dnuc(i)=0.0_real_8
    ENDDO
    DO i=1,3
       dx(i,1)=parm%a1(i)/spar%nr1s
       dx(i,2)=parm%a2(i)/spar%nr2s
       dx(i,3)=parm%a3(i)/spar%nr3s
    ENDDO
    nx0=parap%nrxpl(parai%mepos,1)
    DO iz=1,parm%nr3
       DO iy=1,parm%nr2
          DO i=1,3
             x(i)=(iz-1)*dx(i,3) + (iy-1)*dx(i,2) + (nx0-1)*dx(i,1)
          ENDDO
          DO ix=1,parm%nr1
             rho= -rhor(ix,iy,iz)
             DO i=1,3
                moment%dmom(i)=moment%dmom(i) + (x(i)-rr(i))*rho
                x(i)=x(i)+dx(i,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL dscal(3,scl,moment%dmom,1)
    CALL mp_sum(moment%dmom,3,parai%allgrp)
    CALL tihalt('    RSDIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rsdipo
  ! ==================================================================

END MODULE dipo_utils
