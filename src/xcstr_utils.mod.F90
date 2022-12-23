MODULE xcstr_utils
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

!!$public :: xcstr

!!$contains

END MODULE xcstr_utils

! ==================================================================
SUBROUTINE xcstr(drhovg,rhoe,v,vgc)
  ! ==--------------------------------------------------------------==
  ! ==                        COMPUTES                              ==
  ! ==           EXCHANGE AND CORRELATION ENERGY CONTRIBUTION       ==
  ! ==           OF VANDERBILT CHARGE TO STRESS TENSOR              ==
  ! ==--------------------------------------------------------------==
  ! 
  ! WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
  ! FIXME: this code appears to not have been fully revised for use
  ! with NEWCODE. Some fixes were applied to make it run without
  ! crashing. The results are not totally wrong, but the code needs
  ! to be thouroughly tested (and best checked by someone, who knows
  ! what he is doing).
  ! 
  ! <axel.kohlmeyer@theochem.ruhr-uni-bochum.de> 12/2003.
  ! 
  ! WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
  ! 
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,ncpw,fpar,parm,spar
  USE parac, ONLY : paral,parai
  USE cnst , ONLY:uimag
  USE cppt , ONLY:indz,nzh
  USE strs , ONLY:dexc
  USE spin , ONLY:clsd,spin_mod
  USE geq0mod , ONLY:geq0
  USE reshaper , ONLY:reshape_inplace
  USE fftnew_utils, ONLY : setfftn
  USE fftmain_utils, ONLY : invfftn
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  COMPLEX(real_8)                            :: drhovg(ncpw%nhg,6*clsd%nlsd)
  REAL(real_8)                               :: rhoe(fpar%nnr1,clsd%nlsd)
  COMPLEX(real_8), TARGET                    :: v(fpar%nnr1,clsd%nlsd)
  REAL(real_8)                               :: vgc(fpar%nnr1,clsd%nlsd)

  CHARACTER(*), PARAMETER                    :: procedureN = 'xcstr'

  INTEGER                                    :: ig, ir, isub, kk, ierr
  REAL(real_8)                               :: droe1, droe2, &
                                                vxc1, vxca, vxcb, sxc, vxc
  REAL(real_8), EXTERNAL                     :: ddot
  REAL(real_8), POINTER                      :: vr(:,:,:)
  COMPLEX(real_8), ALLOCATABLE               :: vtmp(:,:)

  CALL tiset(procedureN,isub)
  ALLOCATE(vtmp(fpar%nnr1,clsd%nlsd),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  CALL zeroing(vtmp)
  IF (cntl%tgc) THEN
     DO ir=1,fpar%nnr1
        vgc(ir,1)=REAL(v(ir,1))
        IF (cntl%tlsd) THEN
           vgc(ir,2)=REAL(v(ir,2))
           rhoe(ir,1)=rhoe(ir,1)+rhoe(ir,2)
        ENDIF
     ENDDO
  ENDIF
  CALL setfftn(0)
  CALL reshape_inplace(v, (/2,fpar%nnr1,clsd%nlsd/), vr)
  CALL xcener(sxc,vxc,rhoe,rhoe,vtmp)
  DO kk=1,6,2
     CALL zeroing(v)!,nnr1*nlsd)
     !CDIR NODEP
#ifdef __SR8000
     !poption parallel, tlocal(IG)
#endif 
     DO ig=1,ncpw%nhg
        v(nzh(ig),1)  = drhovg(ig,kk) + uimag*drhovg(ig,kk+1)
        v(indz(ig),1) = CONJG(drhovg(ig,kk)) + uimag*CONJG(drhovg(ig,kk+1))
     ENDDO
     IF (geq0) v(nzh(1),1)  = drhovg(1,kk) + uimag*drhovg(1,kk+1)
     CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
     IF (cntl%tlsd) THEN
        DO ig=1,ncpw%nhg
           v(nzh(ig),2)  = drhovg(ig,6+kk) + uimag*drhovg(ig,6+kk+1)
           v(indz(ig),2) = CONJG(drhovg(ig,6+kk)) + uimag*CONJG(drhovg(ig,6+kk+1))
        ENDDO
        IF (geq0) v(nzh(1),2)  = drhovg(1,6+kk) + uimag*drhovg(1,6+kk+1)
        CALL  invfftn(v(:,2),.FALSE.,parai%allgrp)
     ENDIF
     vxca=0._real_8
     vxcb=0._real_8
     DO ir=1,fpar%nnr1
        droe1=REAL(v(ir,1))
        droe2=AIMAG(v(ir,1))
        vxc1=REAL(vtmp(ir,1))
        vxca=vxca+vxc1*droe1
        vxcb=vxcb+vxc1*droe2
        IF (cntl%tlsd) THEN
           droe1=REAL(v(ir,2))
           droe2=AIMAG(v(ir,2))
           vxc1=REAL(vtmp(ir,2))
           vxca=vxca+vxc1*droe1
           vxcb=vxcb+vxc1*droe2
        ENDIF
        IF (cntl%tgc) THEN
           vxca=vxca+vgc(ir,1)*vr(1,ir,1)
           vxcb=vxcb+vgc(ir,1)*vr(2,ir,1)
           IF (cntl%tlsd) THEN
              vxca=vxca+vgc(ir,2)*vr(1,ir,2)
              vxcb=vxcb+vgc(ir,2)*vr(2,ir,2)
           ENDIF
        ENDIF
     ENDDO
     dexc(kk)   = vxca * parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
     dexc(kk+1) = vxcb * parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
  ENDDO
  IF (cntl%tgc.AND.cntl%tlsd) THEN
     DO ir=1,fpar%nnr1
        rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
     ENDDO
  ENDIF
  DEALLOCATE(vtmp,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE xcstr
! ==================================================================

