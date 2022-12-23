MODULE atrho_utils
  USE atwf,                            ONLY: &
       atchg, atrg, atwf_mod, atwfr, atwp, atwr, catom, loadc_foc_array_size, &
       natsave
  USE bessm_utils,                     ONLY: bessov
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             inyh,&
                                             nzh
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE qspl,                            ONLY: ggnh,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo
  USE reshaper,                        ONLY: reshape_inplace
  USE setbasis_utils,                  ONLY: loadc
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atrho
  PUBLIC :: atdens
  PUBLIC :: give_scr_atrho

CONTAINS

  ! ==================================================================
  SUBROUTINE atrho(rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'atrho'

    CHARACTER(len=8)                         :: form
    COMPLEX(real_8)                          :: tsfac
    COMPLEX(real_8), ALLOCATABLE             :: rhog(:)
    COMPLEX(real_8), POINTER                 :: pcatom(:,:)
    INTEGER                                  :: ia, iaorb, iat, ierr, iform, &
                                                ig, ir, is, isa, isa0, isub, &
                                                ixx, natst
    LOGICAL                                  :: debug
    REAL(real_8)                             :: ar, &
                                                foc(loadc_foc_array_size), &
                                                rdum, rsum, rsum1, sfc, vol
    REAL(real_8), ALLOCATABLE                :: arho(:), datom(:,:), work(:)

    CALL tiset('     ATRHO',isub)
    ! ==--------------------------------------------------------------==
    ! CATOM was allocated as (NGWK, NATTOT), NGWK=2*NGW
    ! but here it is used as (NGW, NATTOT)
    ! also impossible to use f90 pointers as well      
    ! TODO verify with Zurich team: what really happens here with CATOM? 
    CALL reshape_inplace(catom, (/ncpw%ngw, atwp%nattot/), pcatom)
    CALL setfftn(0)
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    ALLOCATE(rhog(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(datom(nsplpo, 2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  INITIALIZE                                                  ==
    ! ==--------------------------------------------------------------==
    vol=1._real_8/parm%omega
    CALL zeroing(rhog)!,nhg)
    isa0=0
    natsave=0
    DO is=1,ions1%nsp
       CALL atdens(is,datom)
#if defined(__VECTOR)
       !$omp parallel do private(IG,IA,ISA,TSFAC,AR) shared(NSPLPO)
#else
       !$omp parallel do private(IG,IA,ISA,TSFAC,AR) shared(NSPLPO)&
       !$omp  schedule(static)
#endif
#ifdef __SR8000
       !poption parallel, tlocal(IG,IA,ISA,TSFAC,AR)
#endif 

       DO ig=1,ncpw%nhg
          tsfac=CMPLX(0._real_8,0._real_8,kind=real_8)
          DO ia=1,ions0%na(is)
             isa=isa0+ia
             tsfac=tsfac+ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
          ENDDO
          ar=curv2(hg(ig),nsplpo,ggnh(1),datom(1,1),datom(1,2),0._real_8)
          rhog(ig)=rhog(ig) + ar*vol*tsfac
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO
    CALL zeroing(psi)!,maxfft)
    !CDIR NODEP
#if defined(__VECTOR)
    !$omp parallel do private(IG)
#else
    !$omp parallel do private(IG) schedule(static)
#endif
#ifdef __SR8000
    !poption parallel, tlocal(IG)
#endif 
    DO ig=1,ncpw%nhg
       psi(indz(ig)) = CONJG(rhog(ig))
       psi(nzh(ig))  = rhog(ig)
    ENDDO
    IF (geq0) psi(nzh(1)) = rhog(1)
    rsum = 0.0_real_8
    IF (geq0) rsum        = REAL(rhog(1))
    CALL invfftn(psi,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! COMPUTE THE INTEGRAL OF THE CHARGE DENSITY IN REAL SPACE
    rsum1=0._real_8
#if defined(__VECTOR)
    !$omp parallel do private(IR,RDUM) reduction(+:RSUM1)
#else
    !$omp parallel do private(IR,RDUM) reduction(+:RSUM1) schedule(static)
#endif
    DO ir=1,fpar%nnr1
       rdum = REAL(psi(ir))
       rhoe(ir)= rdum
       rsum1   = rsum1 + rdum
    ENDDO
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ! ==--------------------------------------------------------------==
    iaorb=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL loadc(pcatom(:,iaorb),foc,ncpw%ngw,ncpw%ngw,atwp%nattot-iaorb+1,&
               SIZE(foc),is,iat,natst)
          DO ixx=iaorb,iaorb+natst-1
             sfc=dotp(ncpw%ngw,pcatom(:,ixx),pcatom(:,ixx))
             CALL mp_sum(sfc,parai%allgrp)
             IF (sfc.EQ.0._real_8) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A,I5,A)') ' ATRHO|',&
                     ' THE NORM OF ATOMIC ORBITAL (',ixx,') IS NULL'
                CALL stopgm('ATRHO','WRONG ATOMIC ORBITAL',& 
                     __LINE__,__FILE__)
             ELSE
                sfc=1._real_8/SQRT(sfc)
             ENDIF
             CALL dscal(2*ncpw%ngw,sfc,pcatom(1,ixx),1)
          ENDDO
          iaorb=iaorb+natst
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! COMPUTE THE INTEGRAL OF THE CHARGE DENSITY IN G-SPACE and R-SPACE
    ! ACROSS ALL PROCESSORS
    rsum = 0.0_real_8
    IF (geq0) rsum = REAL(rhog(1)*parm%omega)
    CALL mp_sum(rsum,parai%allgrp)
    CALL mp_sum(rsum1,parai%allgrp)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       IF ((atwp%nattot*clsd%nlsd).LT.nstate) THEN
          iform=INT(LOG10(REAL(MAX(atwp%nattot*clsd%nlsd,1),kind=real_8)))+1
          WRITE(form,'("(A,I",I1,",A)")') iforM
          WRITE(6,form)&
               ' ATRHO| WARNING! '//&
               ' THE NUMBER OF GENERATED ATOMIC STATE (',atwp%nattot*clsd%nlsd,')'
          iform=INT(LOG10(REAL(MAX(nstate,1),kind=real_8)))+1
          WRITE(form,'("(A,I",I1,",A)")') iforM
          WRITE(6,form)&
               ' ATRHO| WARNING! IS LESS THAN THE NUMBER OF STATE (',&
               nstate,')'
          WRITE(6,'(A,I5,A)') '  ATRHO| WARNING! ',NSTATE-atwp%nattot*clsd%nlsd,&
               ' STATES ARE EQUAL TO ZERO'
          WRITE(6,'(A,I5,A)')&
               ' ATRHO| USE SECTION BASIS TO SPECIFY ATOMIC BASIS'
       ENDIF
       WRITE(6,'(A,F12.6,A,F12.6)')' ATRHO| CHARGE(R-SPACE):',rsum1,&
            ' (G-SPACE):',rsuM
    ENDIF
    ! 
    DEALLOCATE(rhog,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(datom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(arho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt('     ATRHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE atrho
  ! ==================================================================
  SUBROUTINE atdens(is,datom)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: datom(nsplpo,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'atdens'

    INTEGER                                  :: ierr, il, ir, ish, mmax, mused
    REAL(real_8)                             :: currcharge, fac, occu, tmp, xg
    REAL(real_8), ALLOCATABLE                :: arho(:), work(:)

! Variables
! dimension ATRG(maxsys%mmaxx,*)
! dimension ATWFR(maxsys%mmaxx,M1SHLX,*)
! dimension GGNH(*)
! dimension ATCHG(*)
! maxsys%mmaxx
! NSPLPO
! ==--------------------------------------------------------------==

    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(arho)!,maxsys%mmaxx)
    mmax=atwr%meshat(is)
    ! Sum density
    currcharge=atchg(is)
    DO ish=atwf_mod%nshell(is),1,-1
       IF (currcharge.GT.atwf_mod%oc(ish,is)) THEN
          occu=0
          currcharge=currcharge-atwf_mod%oc(ish,is)
       ELSE
          occu=atwf_mod%oc(ish,is)-currcharge
          currcharge=0
       ENDIF
#if defined(__VECTOR)
       !$omp parallel do private(IR)
#else
       !$omp parallel do private(IR) schedule(static)
#endif
#ifdef __SR8000
       !poption parallel, tlocal(IR)
#endif 
       DO ir=1,mmax
          arho(ir)=arho(ir)+occu * (atwfr(ir,ish,is)/atrg(ir,is))**2
       ENDDO
    ENDDO
    xg=0._real_8
    DO ir=mmax,1,-1
       xg=xg+ABS(arho(ir))
       IF (xg.GT.1.e-8_real_8) THEN
          mused = ir
          GOTO 100
       ENDIF
    ENDDO
    mused=mmax
100 CONTINUE
    ! FFT
    CALL zeroing(datom(:,1))!,nsplpo)
    !$omp parallel do private(IL,XG,TMP) shared(MUSED)
#ifdef __SR8000
    !poption parallel, tlocal(IL,XG,TMP)
#endif 
    DO il=nsplpa,nsplpe
       xg=SQRT(ggnh(il))*parm%tpiba
       CALL bessov(atrg(1,is),atwr%clogat(is),mused,arho,0,xg,atrg(mused,is)&
            ,tmp)
       datom(il,1)=tmp
    ENDDO
    CALL mp_sum(datom,nsplpo,parai%allgrp)
    CALL bessov(atrg(1,is),atwr%clogat(is),mused,arho,0,0.0_real_8,atrg(mused,&
         is),tmp)
    IF (tmp.LT.1.e-3_real_8) THEN
       fac=0._real_8
    ELSE
       fac=(ions0%zv(is)-atchg(is))/tmp
    ENDIF
    CALL dscal(nsplpo,fac,datom(1,1),1)
    ! Spline
    ALLOCATE(work(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL curv1(nsplpo,ggnh,datom(1,1),0._real_8,0._real_8,3,datom(1,2),&
         work,0._real_8,ierr)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(arho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE atdens
  ! ==================================================================
  SUBROUTINE give_scr_atrho(latrho,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: latrho
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    latrho=2*ncpw%nhg+2*nsplpo+maxsys%mmaxx+nsplpo+100 ! RHOG DATOM ARHO WORK
    tag  ='2*NHG+3*NSPLPO+MMAXX'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_atrho
  ! ==================================================================

END MODULE atrho_utils
