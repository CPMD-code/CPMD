MODULE copot_utils
  USE corec_utils,                     ONLY: corec
  USE cppt,                            ONLY: inyh,&
                                             nzh
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE gcener_utils,                    ONLY: gcener,&
                                             gclsd
  USE graden_utils,                    ONLY: graden
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corer,&
                                             drhoc,&
                                             vnlcc
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE spin,                            ONLY: clsd
  USE strs,                            ONLY: decc
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             iatpt,&
                                             ncpw,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: copot
  PUBLIC :: give_scr_copot
  PUBLIC :: drhocore

CONTAINS

  ! ==================================================================
  SUBROUTINE copot(rhoe,v,tstress)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE POTENTIAL AND ENERGIES FOR EXCHANGE AND        ==
    ! ==  CORRELATION FOR THE CORE CHARGES                            ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)
    LOGICAL                                  :: tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'copot'

    COMPLEX(real_8), ALLOCATABLE             :: c_vtmp(:), vtemp(:,:)
    INTEGER                                  :: ierr, il_grad, il_vtemp, ir, &
                                                isub, kk
    REAL(real_8)                             :: sgcc, sgcx, sxc, vxc
    REAL(real_8), ALLOCATABLE                :: grad(:,:), vtmp(:,:)

!(nnr1,clsd%nlsd)

    CALL tiset('     COPOT',isub)
    ! ==--------------------------------------------------------------==
    IF (cntl%ttau) THEN
       CALL stopgm('COPOT','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    il_vtemp=2*ncpw%nhg*clsd%nlsd
    il_grad =fpar%nnr1*clsd%nlsd*4      ! GRADEN
    ALLOCATE(vtemp(ncpw%nhg, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(grad(fpar%nnr1, clsd%nlsd*4),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vtmp(fpar%nnr1,clsd%nlsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(c_vtmp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)


    CALL zeroing(vtmp)!,il_vtmp)
    ! ==--------------------------------------------------------------==
    ! == ADD CORE CHARGE TO RHOE                                      ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(rhoe)!,nnr1*clsd%nlsd)
    CALL corec(rhoe,c_vtmp,v)

    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v(:,1))!,nnr1)
    CALL xcener(sxc,vxc,rhoe,rhoe,v)
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    IF (cntl%tgc) THEN
       ! ==--------------------------------------------------------------==
       ! == CALCULATE THE GRADIENT OF THE DENSITY                        ==
       ! ==--------------------------------------------------------------==
       IF (cntl%tlsd) THEN
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,1),vtemp(1,1),nzh)
          CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(1,2),nzh)
          DO ir=1,fpar%nnr1
             rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          ENDDO
          CALL graden(rhoe(:,1),v,grad(1,1),c_vtmp)
          CALL graden(rhoe(:,2),v,grad(1,5),c_vtmp)
       ELSE
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v,vtemp,nzh)
          CALL graden(rhoe(:,1),v,grad,c_vtmp)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! == GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)            ==
       ! ==--------------------------------------------------------------==
       IF (cntl%tlsd) THEN
          CALL gclsd(sgcx,sgcc,rhoe,v,vtemp,vtmp,grad,.FALSE.)
       ELSE
          CALL gcener(sgcx,sgcc,rhoe,v,vtemp,vtmp(:,1),grad,.FALSE.)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == TRANSF. CORE  POTENTIAL IN G-SPACE                           ==
    ! == PUT THE CORE  POTENTIAL IN G-SPACE INTO VNLCC                ==
    ! ==--------------------------------------------------------------==
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    CALL zgthr(ncpw%nhg,v,vnlcc,nzh)
    IF (cntl%tlsd) THEN
       CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,2),vnlcc(1,2),nzh)
    ENDIF
    ! ==--------------------------------------------------------------==
    corer%excco=sxc
    corer%egcxco=sgcx
    corer%egccco=sgcc
    IF (tstress) THEN
       DO kk=1,6
          CALL drhocore(vtemp,drhoc(1,1,kk),eigrb,ei1,ei2,ei3,inyh)
          decc(kk)=-dotp(ncpw%nhg,vtemp(:,1),vnlcc(:,1))
       ENDDO
    ENDIF
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(c_vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt('     COPOT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE copot
  ! ==================================================================
  SUBROUTINE give_scr_copot(lcopot,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lcopot
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lcopot=2*ncpw%nhg*clsd%nlsd+fpar%nnr1*clsd%nlsd*4+MAX(2*ncpw%nhg,fpar%nnr1*clsd%nlsx)+20
    tag='2*NHG*NLSD+NNR1*NLSD*4+...'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_copot
  ! ==================================================================
  SUBROUTINE drhocore(dco,drhoc,eigrb,ei1,ei2,ei3,inyh)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: dco(*)
    REAL(real_8)                             :: drhoc(ncpw%nhg,ions1%nsp)
    COMPLEX(real_8) :: eigrb(ncpw%nhg,ions1%nat), &
      ei1(ions1%nat,(2*spar%nr1s-1)), ei2(ions1%nat,(2*spar%nr2s-1)), &
      ei3(ions1%nat,(2*spar%nr3s-1))
    INTEGER                                  :: inyh(3,ncpw%nhg)

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: ia, ig, is, isa
    REAL(real_8)                             :: ei, er

    CALL zeroing(dco(1:ncpw%nhg))!,nhg)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(ISA,IA,IS,IG,ER,EI)
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa)
             is=iatpt(2,isa)
             er=REAL(eigrb(ig,isa))
             ei=AIMAG(eigrb(ig,isa))
             dco(ig)=dco(ig)+CMPLX(er*drhoc(ig,is),ei*drhoc(ig,is),kind=real_8)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,ISA,IA,IS,ER,EI,EI123)
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa)
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             er=REAL(ei123)
             ei=AIMAG(ei123)
             dco(ig)=dco(ig)+CMPLX(er*drhoc(ig,is),ei*drhoc(ig,is),kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE drhocore
  ! ==================================================================


END MODULE copot_utils
