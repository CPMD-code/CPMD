MODULE td_prop_utils
  USE atomc_utils,                     ONLY: atomc
  USE cppt,                            ONLY: indz,&
                                             nr1h,&
                                             nr2h,&
                                             nr3pl,&
                                             nzh
  USE eicalc_utils,                    ONLY: eicalc
  USE elstpo_utils,                    ONLY: elstpo
  USE error_handling,                  ONLY: stopgm
  USE espchg_utils,                    ONLY: atfield,&
                                             espsolv,&
                                             printcg,&
                                             selectp
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tdcharge
  PUBLIC :: tddipo

CONTAINS

  ! ==================================================================
  SUBROUTINE tdcharge(rhoe,tau0,psi)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:), tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'tdcharge'

    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), qphi(:), &
                                                vtemp(:)
    INTEGER                                  :: i, ierr, ig, ippc, isub, len, &
                                                lqphi
    INTEGER, ALLOCATABLE                     :: isel(:)
    REAL(real_8), ALLOCATABLE                :: achrg(:), echrg(:), &
                                                efield(:), reirop(:), rscr(:)

! ==--------------------------------------------------------------==
! ==  CALCULATE ATOMIC CHARGES                                    ==
! ==--------------------------------------------------------------==

    CALL tiset('  TDCHARGE',isub)
    CALL setfftn(0)
    ALLOCATE(achrg(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(echrg(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(reirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtemp(MAX(maxfftn,ncpw%nhg)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (isos1%tclust .AND. (isos3%ps_type.EQ.1)) THEN
       lqphi = (nr1h+1)*nr2h*nr3pl
    ELSE
       lqphi = 1
    ENDIF
    ALLOCATE(qphi(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rscr(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 

    CALL atomc(rhoe,rscr,tau0,achrg)

    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtemp(ig) = psi(nzh(ig))
    ENDDO
    CALL eicalc(eivps,eirop)
    CALL elstpo(vtemp,eirop,eivps,psi,qphi)
    CALL dcopy(2*ncpw%nhg,psi,1,vtemp,1)
    CALL zeroing(psi)!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       psi(indz(ig)) = CONJG(vtemp(ig))
       psi(nzh(ig))  = vtemp(ig)
    ENDDO
    IF (geq0.AND.isos1%tclust) THEN
       psi(nzh(1)) = vtemp(1)
    ELSEIF (geq0) THEN
       psi(nzh(1)) = CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDIF
    CALL  invfftn(psi,.FALSE.,parai%allgrp)
    DO i=1,fpar%nnr1
       rscr(i)=REAL(psi(i))
    ENDDO
    len=fpar%nnr1
    ALLOCATE(isel(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL selectp(isel,tau0,ippc)
    len=(ions1%nat+1)*ippc
    ALLOCATE(efield(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO i=1,ippc
       efield(i)=rscr(isel(i))
    ENDDO
    CALL atfield(efield,psi,vtemp,reirop,qphi,isel,ippc)
    CALL espsolv(efield,echrg,ippc)
    IF (paral%parent) CALL printcg(tau0,achrg,echrg)
    ! 
    DEALLOCATE(achrg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(echrg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(isel,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qphi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(efield,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(reirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! 
    CALL tihalt('  TDCHARGE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tdcharge

  ! ==================================================================
  SUBROUTINE tddipo(rdmom,rhoe,tau0,psi)
    ! ==--------------------------------------------------------------==
    USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
    USE error_handling, ONLY: stopgm
    USE timer, ONLY: tiset, tihalt
    USE mp_interface, ONLY: mp_sum
    USE fft_maxfft,                      ONLY: maxfftn
    USE system , ONLY:maxsys,ncpw,fpar,parm
    USE parac, ONLY : paral,parai
    USE cppt , ONLY:gk,nzh
    USE ions , ONLY:ions0,ions1
    USE isos , ONLY:isos1,isos3
    USE cnst , ONLY:au_deb,uimag
    USE geq0mod , ONLY:geq0
    USE fftmain_utils, ONLY : fwfftn,invfftn
    REAL(real_8), INTENT(OUT)                :: rdmom(:)
    REAL(real_8)                             :: rhoe(:), tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(:)

    COMPLEX(real_8), PARAMETER               :: cone = (1._real_8,0._real_8) 

    COMPLEX(real_8)                          :: cg, dmom(3), dr(3), emax(3), &
                                                emin(3), gmax(3), gmin(3), &
                                                rdr(3)
    INTEGER                                  :: i, ia, ig, ig1, ipm, is, &
                                                isub, k
    REAL(real_8)                             :: d1, d2, d3, dd, g(3), ogk, &
                                                rcc(3), rmax(3), rmin(3), &
                                                x0(3), zvtot

!(fpar%nnr1), &
!(maxfftn)

    CALL tiset('    TDDIPO',isub)
    ! 
    IF (parm%ibrav.NE.1 .AND. parm%ibrav.NE.8) THEN
       CALL stopgm("TDDIPO","only works with cubic cells",& 
            __LINE__,__FILE__)
    ENDIF
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1)&
         CALL stopgm('TDDIPO','HOCKNEY NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! infact we only need to set the origin to the center of box and
    ! also calculate the ionic contribution
    ! 
    x0(1)=0._real_8
    x0(2)=0._real_8
    x0(3)=0._real_8
    zvtot=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          x0(1)=x0(1)+ions0%zv(is)*tau0(1,ia,is)
          x0(2)=x0(2)+ions0%zv(is)*tau0(2,ia,is)
          x0(3)=x0(3)+ions0%zv(is)*tau0(3,ia,is)
       ENDDO
       zvtot=zvtot+ions0%na(is)*ions0%zv(is)
    ENDDO
    x0(1)=x0(1)/zvtot
    x0(2)=x0(2)/zvtot
    x0(3)=x0(3)/zvtot
    ! 
    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    ! 
    rmin(1)=x0(1)-0.5_real_8*parm%a1(1)
    rmax(1)=x0(1)+0.5_real_8*parm%a1(1)
    rmin(2)=x0(2)-0.5_real_8*parm%a2(2)
    rmax(2)=x0(2)+0.5_real_8*parm%a2(2)
    rmin(3)=x0(3)-0.5_real_8*parm%a3(3)
    rmax(3)=x0(3)+0.5_real_8*parm%a3(3)
    ! RCC: coordinates of the center of charge
    DO k=1,3
       rcc(k)=(rmin(k)+rmax(k))/2
       dmom(k)=0._real_8
    ENDDO
    ! calculate integration range relative to center of charge
    ! sum over all components of the wave function
    ig1=1
    IF (geq0) THEN
       ig1=2
       cg=psi(nzh(1))
       DO k=1,3
          rdr(k)=0.5_real_8*(rmax(k)*rmax(k)-rmin(k)*rmin(k))
          dr(k)=rmax(k)-rmin(k)
       ENDDO
       dmom(1)=dmom(1)-cg*dr(2)*dr(3)*(rdr(1)-dr(1)*rcc(1))
       dmom(2)=dmom(2)-cg*dr(1)*dr(3)*(rdr(2)-dr(2)*rcc(2))
       dmom(3)=dmom(3)-cg*dr(1)*dr(2)*(rdr(3)-dr(3)*rcc(3))
    ENDIF
    DO ig=ig1,ncpw%nhg
       ! sum twice (for positive and negative g values)
       DO ipm=-1,1,2
          DO k=1,3
             g(k)=gk(k,ig)*parm%tpiba*ipm
             gmin(k)=uimag*g(k)*rmin(k)
             gmax(k)=uimag*g(k)*rmax(k)
             emin(k)=EXP(gmin(k))
             emax(k)=EXP(gmax(k))
             IF (ABS(gk(k,ig)).GT.1.e-8_real_8) THEN
                ogk=1._real_8/g(k)
                rdr(k)=-((gmax(k)-cone)*emax(k)-&
                     (gmin(k)-cone)*emin(k))*ogk*ogk
                dr(k)=-uimag*(emax(k)-emin(k))*ogk
             ELSE
                rdr(k)=0.5_real_8*(rmax(k)*rmax(k)-rmin(k)*rmin(k))
                dr(k)=rmax(k)-rmin(k)
             ENDIF
          ENDDO
          cg=CMPLX(REAL(psi(nzh(ig))),AIMAG(psi(nzh(ig)))*ipm,kind=real_8)
          dmom(1)=dmom(1)-cg*dr(2)*dr(3)*(rdr(1)-dr(1)*rcc(1))
          dmom(2)=dmom(2)-cg*dr(1)*dr(3)*(rdr(2)-dr(2)*rcc(2))
          dmom(3)=dmom(3)-cg*dr(1)*dr(2)*(rdr(3)-dr(3)*rcc(3))
       ENDDO
    ENDDO
    CALL mp_sum(dmom,3,parai%allgrp)
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' DIPOLE MOMENT '
       d1=REAL(dmom(1))
       d2=REAL(dmom(2))
       d3=REAL(dmom(3))
       dd=SQRT(d1*d1+d2*d2+d3*d3)
       IF (paral%io_parent)&
            WRITE(6,'(3(10X,A),8X,A)') ' X',' Y',' Z',' TOTAL'
       IF (paral%io_parent)&
            WRITE(6,'(3F12.5,F15.5,A)') d1,d2,d3,dd,'  atomic units'
       IF (paral%io_parent)&
            WRITE(6,'(3F12.5,F15.5,A)') au_deb*d1,au_deb*d2,au_deb*d3,&
            au_deb*dd,'  Debye'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    rdmom=REAL(dmom)
    ! 
    CALL tihalt('    TDDIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tddipo

END MODULE td_prop_utils

