MODULE espchg_utils
  USE adat,                            ONLY: covrad,&
                                             elem
  USE cnst,                            ONLY: fbohr,&
                                             fpi
  USE cppt,                            ONLY: hg,&
                                             hipz,&
                                             indz,&
                                             inyh,&
                                             nzh,&
                                             scg
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE hip_utils,                       ONLY: hip
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: selectp
  PUBLIC :: printcg
  PUBLIC :: atfield
  PUBLIC :: espsolv

CONTAINS

  ! ==================================================================
  SUBROUTINE selectp(isel,tau0,ippc)
    ! ==--------------------------------------------------------------==
    ! ==            SELECT POINTS IN REAL SPACE FOR FIT               ==
    ! ==--------------------------------------------------------------==
    INTEGER :: isel(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: ippc

    INTEGER                                  :: i, i1, i2, i3, ia, ii, ippcs, &
                                                is, jj, ntest
    REAL(real_8)                             :: da(3,3), rdp, rr, rvdw1, &
                                                rvdw2, x, x1, x2, xpp, y, y1, &
                                                y2, z, z1, z2

    DO i=1,3
       da(i,1)=parm%a1(i)/REAL(spar%nr1s,kind=real_8)
       da(i,2)=parm%a2(i)/REAL(spar%nr2s,kind=real_8)
       da(i,3)=parm%a3(i)/REAL(spar%nr3s,kind=real_8)
    ENDDO
    CALL zeroing(isel)!,nnr1)
    DO i3=1,spar%nr3s
       DO i2=1,spar%nr2s
          DO i1=1,parm%nr1
             isel(i1,i2,i3)=-1
          ENDDO
       ENDDO
    ENDDO
    DO is=1,ions1%nsp
       rvdw1=(2._real_8*covrad(ions0%iatyp(is))*fbohr)**2
       rvdw2=(4._real_8*covrad(ions0%iatyp(is))*fbohr)**2
       DO ia=1,ions0%na(is)
          x1=tau0(1,ia,is)
          y1=tau0(2,ia,is)
          z1=tau0(3,ia,is)
          !$omp parallel do private(I1,I2,I3,II,X,Y,Z,X2,Y2,Z2,RR) shared(ISEL)
          DO i3=1,spar%nr3s
             DO i2=1,spar%nr2s
                DO ii=1,parm%nr1
                   i1=parap%nrxpl(parai%mepos,1)+ii-1
                   x=i1*da(1,1)+i2*da(1,2)+i3*da(1,3)-x1
                   y=i1*da(2,1)+i2*da(2,2)+i3*da(2,3)-y1
                   z=i1*da(3,1)+i2*da(3,2)+i3*da(3,3)-z1
                   CALL pbc(x,y,z,x2,y2,z2,1,parm%apbc,parm%ibrav)
                   rr=x2*x2+y2*y2+z2*z2
                   IF (rr.LE.rvdw1) isel(ii,i2,i3)=0
                   IF (rr.GT.rvdw1.AND.rr.LT.rvdw2 .AND.&
                        isel(ii,i2,i3).NE.0) isel(ii,i2,i3)=1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ..density of points
    rdp=REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)/parm%omega
    IF (rdp.GT.50._real_8) THEN
       !$omp parallel do private(I1,I2,I3,II,NTEST)
       DO i3=1,spar%nr3s
          DO i2=1,spar%nr2s
             DO ii=1,parm%nr1
                i1=parap%nrxpl(parai%mepos,1)+ii-1
                ntest=MOD(i1,2)*MOD(i2,2)*MOD(i3,2)
                IF (ntest.EQ.0) isel(ii,i2,i3)=0
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    xpp=0._real_8
    ippc=0
    !$omp parallel do private(I1,I2,I3) reduction(+:IPPC)
    DO i3=1,spar%nr3s
       DO i2=1,spar%nr2s
          DO i1=1,parm%nr1
             IF (isel(i1,i2,i3).GT.0) THEN
                ippc=ippc+1
             ELSE
                isel(i1,i2,i3)=0
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    xpp=REAL(ippc,kind=real_8)
    CALL mp_sum(xpp,parai%allgrp)
    ippcs=ANINT(xpp)
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,I10,F16.4)')&
         ' ESP CHARGES| NUMBER OF FITTING POINTS ',ippcs, xpp
    IF (ippcs.LT.4*ions1%nat.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'ESP CHARGES| WARNING! NOT ENOUGH FITTING POINTS.'
       IF (paral%io_parent)&
            WRITE(6,*) 'ESP CHARGES| WARNING! BULK SYSTEMS NOT SUPPORTED.'
       IF (paral%io_parent)&
            WRITE(6,*) 'ESP CHARGES| ESP CHARGES DISABLED.'
       CALL zeroing(isel)!,nnr1)
       RETURN
    ENDIF
    ii=0
    jj=0
    DO i3=1,fpar%kr3s
       DO i2=1,fpar%kr2s
          DO i1=1,fpar%kr1
             jj=jj+1
             IF (isel(i1,i2,i3).GT.0) THEN
                ii=ii+1
                CALL settoint(isel,ii,jj)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE selectp
  ! ==================================================================
  SUBROUTINE settoint(a,i,b)
    INTEGER                                  :: a(*), i, b

    a(i)=b
  END SUBROUTINE settoint
  ! ==================================================================
  SUBROUTINE printcg(tau0,achrg,echrg)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), &
                                                achrg(ions1%nat), &
                                                echrg(ions1%nat)

    INTEGER                                  :: ia, iat, is, k
    REAL(real_8)                             :: chksum

    IF (paral%io_parent)&
         WRITE(6,'(/,1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(A,10X,A,/,15X,A)')&
         '    ATOM  ',' COORDINATES                   CHARGES'&
         ,' X         Y         Z             INT         ESP'
    iat=0
    chksum=0.0_real_8
    IF (paral%io_parent) OPEN(unit=49,file='CHOUT',status='UNKNOWN')
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          IF (paral%io_parent)&
               WRITE(6,'(1X,I4,1X,A2,1X,3F10.4,1X,2F12.6)')&
               iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3),achrg(iat),echrg(iat)
          IF (paral%io_parent) WRITE(49,'(F10.6)') echrg(iat)
          chksum=chksum+ABS(achrg(iat))
       ENDDO
    ENDDO
    IF (paral%io_parent) THEN
       CLOSE(49)
       WRITE(6,'(1X,A,E12.5)') 'ChkSum(CHARGES) =',chksum
       WRITE(6,'(/,1X,64("*"))')
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printcg
  ! ==================================================================
  SUBROUTINE atfield(efield,v,vtemp,rhocc,qphi,isel,ippc)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: v(:), vtemp(:)
    REAL(real_8)                             :: rhocc(:)
    COMPLEX(real_8)                          :: qphi(:)
    INTEGER                                  :: isel(:), ippc
    REAL(real_8)                             :: efield(ippc,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'atfield'

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: i, ia, ig, ig1, is, isa, isub
    REAL(real_8)                             :: ei, er, fpibg, qmax, r2max, &
                                                tpit, vol

    CALL tiset(procedureN,isub)

    CALL setfftn(0)
    isa=0
    vol=1._real_8/parm%omega
    DO is=1,ions1%nsp
       r2max=(covrad(ions0%iatyp(is))*fbohr)**2
       DO ig=1,ncpw%nhg
          qmax=0.25_real_8*r2max*hg(ig)*parm%tpiba2
          rhocc(ig)=-vol*EXP(-qmax)
       ENDDO
       DO ia=1,ions0%na(is)
          isa=isa+1
          CALL zeroing(vtemp)!,nhg)
          IF (cntl%bigmem) THEN
             !$omp parallel do private(IG)
             DO ig=1,ncpw%nhg
                vtemp(ig)=rhocc(ig)*eigrb(ig,isa)
             ENDDO
          ELSE
             !$omp parallel do private(IG,EI123,ER,EI)
             DO ig=1,ncpw%nhg
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                vtemp(ig)=CMPLX(er*rhocc(ig),ei*rhocc(ig),kind=real_8)
             ENDDO
          ENDIF
          IF (isos1%tclust) THEN
             IF (isos3%ps_type.EQ.1) THEN
                CALL zeroing(v)!,maxfft)
                !$omp parallel do private(IG)
                !CDIR NODEP
                DO ig=1,ncpw%nhg
                   v(indz(ig)) = CONJG(hipz(ig)*vtemp(ig))
                   v(nzh(ig))  = hipz(ig)*vtemp(ig)
                ENDDO
                IF (geq0) v(nzh(1)) = vtemp(1)
                CALL  invfftn(v,.FALSE.,parai%allgrp)
                !$omp parallel do private(I)
                DO i=1,ippc
                   efield(i,isa+1)=REAL(v(isel(i)))
                ENDDO
                CALL zeroing(v)!,maxfft)
                !$omp parallel do private(IG)
                !CDIR NODEP
                DO ig=1,ncpw%nhg
                   v(indz(ig)) = CONJG(vtemp(ig))
                   v(nzh(ig))  = vtemp(ig)
                ENDDO
                IF (geq0) v(nzh(1)) = vtemp(1)
                CALL  invfftn(v,.FALSE.,parai%allgrp)
                CALL hip(v,qphi)
                !$omp parallel do private(I)
                DO i=1,ippc
                   efield(i,isa+1)=efield(i,isa+1)+REAL(v(isel(i)))
                ENDDO
             ELSEIF (isos3%ps_type.EQ.2.OR.isos3%ps_type.EQ.3) THEN
                IF (geq0) THEN
                   vtemp(1)=CMPLX(0._real_8,0._real_8,kind=real_8)
                   ig1=2
                ELSE
                   ig1=1
                ENDIF
                !$omp parallel do private(IG)
                DO ig=ig1,ncpw%nhg
                   vtemp(ig)=vtemp(ig)*scg(ig)
                ENDDO
                CALL zeroing(v)!,maxfft)
                !$omp parallel do private(IG)
                !CDIR NODEP
                DO ig=1,ncpw%nhg
                   v(indz(ig)) = CONJG(vtemp(ig))
                   v(nzh(ig))  = vtemp(ig)
                ENDDO
                IF (geq0) v(nzh(1)) = vtemp(1)
                CALL  invfftn(v,.FALSE.,parai%allgrp)
                !$omp parallel do private(I)
                DO i=1,ippc
                   efield(i,isa+1)=REAL(v(isel(i)))
                ENDDO
             ELSE
                CALL stopgm('ESPCHG','PS_TYPE not implemented',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSE
             IF (geq0) THEN
                vtemp(1)=CMPLX(0._real_8,0._real_8,kind=real_8)
                ig1=2
             ELSE
                ig1=1
             ENDIF
             tpit=fpi/parm%tpiba2
             !$omp parallel do private(IG,FPIBG)
             DO ig=ig1,ncpw%nhg
                fpibg=tpit/hg(ig)
                vtemp(ig)=CMPLX(fpibg,kind=real_8)*vtemp(ig)
             ENDDO
             CALL zeroing(v)!,maxfft)
             !$omp parallel do private(IG,FPIBG)
             !CDIR NODEP
             DO ig=1,ncpw%nhg
                v(indz(ig)) = CONJG(vtemp(ig))
                v(nzh(ig))  = vtemp(ig)
             ENDDO
             IF (geq0) v(nzh(1)) = vtemp(1)
             CALL  invfftn(v,.FALSE.,parai%allgrp)
             !$omp parallel do private(I)
             DO i=1,ippc
                efield(i,isa+1)=REAL(v(isel(i)))
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE atfield
  ! ==================================================================
  SUBROUTINE espsolv(efield,echrg,ippc)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: echrg(:)
    INTEGER                                  :: ippc
    REAL(real_8)                             :: efield(ippc,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'espsolv'

    INTEGER                                  :: i, ierr, ij, info, ippcs, kk, &
                                                len
    INTEGER, ALLOCATABLE                     :: ipiv(:)
    LOGICAL                                  :: fexist
    REAL(real_8)                             :: den, diff, rrms, rrmso, xpp
    REAL(real_8), ALLOCATABLE                :: amat(:,:), chold(:), &
                                                czero(:), qpot(:)

!ions1%nat+1)

    rrms = HUGE(0.0_real_8)

    ALLOCATE(czero(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(chold(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qpot(ippc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    len=(ions1%nat+1)**2
    ALLOCATE(amat(ions1%nat+1,ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipiv(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(chold)!,ions1%nat+1)
    CALL zeroing(czero)!,ions1%nat+1)
    CALL zeroing(amat)!,len)
    IF (ippc.GT.0) THEN
       CALL dgemv('T',ippc,ions1%nat,1._real_8,efield(1,2),ippc,efield,1,0._real_8,&
            echrg,1)
    ELSE
       CALL zeroing(echrg)!,ions1%nat+1)
    ENDIF
    CALL mp_sum(echrg,ions1%nat,parai%allgrp)
    IF (cnti%lfit.EQ.0) THEN
       echrg(ions1%nat+1)=crge%charge
       len=(ions1%nat+1)**2
       IF(ALLOCATED(amat)) THEN
          DEALLOCATE(amat,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)   
       ENDIF
       ALLOCATE(amat(ions1%nat+1,ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF(ALLOCATED(ipiv)) THEN
          DEALLOCATE(ipiv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(ipiv(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(amat)!,len)
       IF (ippc.GT.0)&
            CALL dgemm('T','N',ions1%nat,ions1%nat,ippc,1._real_8,efield(1,2),ippc,&
            efield(1,2),ippc,0._real_8,amat,ions1%nat+1)
       CALL mp_sum(amat,len,parai%allgrp)
       DO i=1,ions1%nat
          amat(ions1%nat+1,i)=1._real_8
          amat(i,ions1%nat+1)=1._real_8
       ENDDO
       CALL dgesv(ions1%nat+1,1,amat,ions1%nat+1,ipiv,echrg,ions1%nat+1,info)
       DEALLOCATE(amat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ipiv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(czero,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(chold,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qpot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (cnti%lfit.EQ.1) THEN
       IF (paral%parent) THEN
          INQUIRE(file='CZERO',exist=fexist)
          IF (fexist) THEN
             OPEN(unit=49,file='CZERO',status='UNKNOWN')
             DO i=1,ions1%nat
                READ(49,*) czero(i)
             ENDDO
             CLOSE(49)
             WRITE(6,'(A)') ' REFERENCE CHARGES FROM CZERO'
          ELSE
             DO i=1,ions1%nat
                czero(i)=0._real_8
             ENDDO
             WRITE(6,'(A)') ' REFERENCE CHARGES EQUAL ZERO'
          ENDIF
       ENDIF
       xpp=REAL(ippc,kind=real_8)
       CALL mp_sum(xpp,parai%allgrp)
       ippcs=NINT(xpp)
       CALL add_rest1_p(echrg,chold,czero,ions1%nat,ippcs)
       echrg(ions1%nat+1)=crge%charge
       CALL zeroing(amat)!,len)
       IF (ippc.GT.0)&
            CALL dgemm('T','N',ions1%nat,ions1%nat,ippc,1._real_8,efield(1,2),ippc,&
            efield(1,2),ippc,0._real_8,amat,ions1%nat+1)
       CALL mp_sum(amat,len,parai%allgrp)
       CALL add_rest2_p(amat,chold,czero,ions1%nat,ippcs)
       DO i=1,ions1%nat
          amat(ions1%nat+1,i)=1._real_8
          amat(i,ions1%nat+1)=1._real_8
       ENDDO
       CALL dgesv(ions1%nat+1,1,amat,ions1%nat+1,ipiv,echrg,ions1%nat+1,info)
       CALL dcopy(ions1%nat+1,echrg,1,chold,1)
       IF (ippc.GT.0)&
            CALL dgemv('N',ippc,ions1%nat,1._real_8,efield(1,2),ippc,echrg,1,0._real_8,&
            qpot,1)
       rrmso=rrms
       diff=0._real_8
       den=0._real_8
       DO ij=1,ippc
          diff=diff+(efield(ij,1)-qpot(ij))**2
          den=den+(efield(ij,1))**2
       ENDDO
       CALL mp_sum(diff,parai%allgrp)
       CALL mp_sum(den,parai%allgrp)
       rrms=diff/den
       rrms=SQRT(rrms)
       IF (paral%parent) THEN
          WRITE(6,'(A,10X,1F12.8)')&
               '      RRMS=    ',rrmS
       ENDIF
       DEALLOCATE(amat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ipiv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(czero,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qpot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(chold,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (cnti%lfit.EQ.2) THEN
       xpp=REAL(ippc,kind=real_8)
       CALL mp_sum(xpp,parai%allgrp)
       ippcs=NINT(xpp)
       kk=0
111    CONTINUE
       kk=kk+1
       CALL add_rest1_h(echrg,chold,czero,ions1%nat,ippcs)
       echrg(ions1%nat+1)=crge%charge
       CALL zeroing(amat)!,len)
       IF (ippc.GT.0)&
            CALL dgemm('T','N',ions1%nat,ions1%nat,ippc,1._real_8,efield(1,2),ippc,&
            efield(1,2),ippc,0._real_8,amat,ions1%nat+1)
       CALL mp_sum(amat,len,parai%allgrp)
       CALL add_rest2_h(amat,chold,czero,ions1%nat,ippcs)
       DO i=1,ions1%nat
          amat(ions1%nat+1,i)=1._real_8
          amat(i,ions1%nat+1)=1._real_8
       ENDDO
       CALL dgesv(ions1%nat+1,1,amat,ions1%nat+1,ipiv,echrg,ions1%nat+1,info)
       CALL dcopy(ions1%nat+1,echrg,1,chold,1)
       IF (ippc.GT.0)&
            CALL dgemv('N',ippc,ions1%nat,1._real_8,efield(1,2),ippc,echrg,1,0._real_8,&
            qpot,1)
       rrmso=rrms
       diff=0._real_8
       den=0._real_8
       DO ij=1,ippc
          diff=diff+(efield(ij,1)-qpot(ij))**2
          den=den+(efield(ij,1))**2
       ENDDO
       CALL mp_sum(diff,parai%allgrp)
       CALL mp_sum(den,parai%allgrp)
       rrms=diff/den
       rrms=SQRT(rrms)
       IF (paral%parent) THEN
          WRITE(6,'(A,10X,1F12.8)')&
               '      RRMS=    ',rrmS
       ENDIF
       IF (ABS(rrms-rrmso).GT.1.e-4_real_8.AND.kk.LT.5)  GOTO 111
       DEALLOCATE(amat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ipiv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(czero,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qpot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(chold,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE espsolv
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE  add_rest1_p(echrg,chold,czero,nat,ippcs)
    REAL(real_8), DIMENSION(*)               :: echrg, chold, czero
    INTEGER                                  :: nat, ippcs

    INTEGER                                  :: i
    REAL(real_8)                             :: res

    res=0.0000058_real_8*ippcs
    DO i=1,nat
       IF (czero(i).LT.100) THEN
          echrg(i)=echrg(i)+res*czero(i)
       ENDIF
    ENDDO
  END SUBROUTINE add_rest1_p
  ! ==================================================================
  SUBROUTINE  add_rest2_p(amat,chold,czero,nat,ippcs)
    INTEGER                                  :: nat
    REAL(real_8), DIMENSION(*)               :: czero, chold
    REAL(real_8), DIMENSION(nat+1, nat+1)    :: amat
    INTEGER                                  :: ippcs

    INTEGER                                  :: i
    REAL(real_8)                             :: res

    res=0.0000058_real_8*ippcs
    DO i=1,nat
       IF (czero(i).LT.100) THEN
          amat(i,i)=amat(i,i)+res
       ELSE
          amat(i,i)=amat(i,i)+1000._real_8
       ENDIF
    ENDDO
  END SUBROUTINE add_rest2_p
  ! ==================================================================
  SUBROUTINE  add_rest1_h(echrg,chold,czero,nat,ippcs)
    REAL(real_8), DIMENSION(*)               :: echrg, chold, czero
    INTEGER                                  :: nat, ippcs

    INTEGER                                  :: i
    REAL(real_8)                             :: fact, res1, res2

    res1=0.0000058_real_8*ippcs
    res2=0.1
    DO i=1,nat
       fact=0.5_real_8/SQRT(chold(i)**2+res2**2)
       echrg(i)=echrg(i)-res1*chold(i)*fact
    ENDDO
  END SUBROUTINE add_rest1_h
  ! ==================================================================
  SUBROUTINE  add_rest2_h(amat,chold,czero,nat,ippcs)
    INTEGER                                  :: nat
    REAL(real_8), DIMENSION(*)               :: czero, chold
    REAL(real_8), DIMENSION(nat+1, nat+1)    :: amat
    INTEGER                                  :: ippcs

    INTEGER                                  :: i
    REAL(real_8)                             :: fact, res1, res2

    res1=0.0000058*ippcs
    res2=0.1
    DO i=1,nat
       fact=1._real_8/SQRT(chold(i)**2+res2**2)
       amat(i,i)=amat(i,i)+res1*chold(i)*fact
    ENDDO
  END SUBROUTINE add_rest2_h
  ! ==================================================================

END MODULE espchg_utils
