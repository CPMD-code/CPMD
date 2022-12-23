MODULE vhk_utils
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE gcener_utils,                    ONLY: gcener,&
                                             gclsd
  USE graden_utils,                    ONLY: graden
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vhk
  PUBLIC :: give_scr_vhk

CONTAINS

  ! ==================================================================
  SUBROUTINE vhk(rhoe,fhk)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  FOR THE THOMAS-FERMI FUNCTIONAL OF THE KINETIC ENERGY       ==
    ! ==  AND THE EXCHANGE CORRELATION FUNCTIONAL                     ==
    ! ==--------------------------------------------------------------==
    ! == RHOE:  in electronic density in real space                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                fhk(fpar%nnr1,clsd%nlsd)

    CHARACTER(*), PARAMETER                  :: procedureN = 'vhk'

    COMPLEX(real_8), ALLOCATABLE             :: v(:,:), vtemp(:,:), vtmp(:,:)
    INTEGER                                  :: i, ierr, il_grad, il_vtmp, &
                                                ir, isub
    REAL(real_8)                             :: roe, sgcc, sgcx, sxc, vxc
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rvtmp(:,:)

    CALL tiset('       VHK',isub)
    ! ==--------------------------------------------------------------==
    IF (cntl%ttau) THEN
       CALL stopgm('VHK','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    sxc=0._real_8
    sgcx=0._real_8
    sgcc=0._real_8
    vxc=0._real_8
    il_vtmp = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx
    IF (cntl%tgc) THEN
       il_grad=fpar%nnr1*clsd%nlsd*4
    ELSE
       il_grad=0
    ENDIF
    ALLOCATE(v(maxfft, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vtemp(ncpw%nhg , clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vtmp(ncpw%nhg, il_vtmp/ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rvtmp(fpar%nnr1,clsd%nlsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tgc) THEN
       ALLOCATE(grad(fpar%nnr1, clsd%nlsd*4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) CALL stopgm("VHK","NLCC NOT SUPPORTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft*clsd%nlsd)
    CALL xcener(sxc,vxc,rhoe,rhoe,v)
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          roe=MAX(rhoe(ir,1),0.0_real_8)
          vxc = 4.785_real_8*roe**(0.666666666667_real_8)
          v(ir,1)=v(ir,1)+CMPLX(vxc,0.0_real_8,kind=real_8)
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          roe=MAX(rhoe(ir,1)-rhoe(ir,2),0.0_real_8)
          vxc = 7.596_real_8*roe**(0.666666666667_real_8)
          v(ir,1)=v(ir,1)+CMPLX(vxc,0.0_real_8,kind=real_8)
          roe=MAX(rhoe(ir,2),0.0_real_8)
          vxc = 7.596_real_8*roe**(0.666666666667_real_8)
          v(ir,2)=v(ir,2)+CMPLX(vxc,0.0_real_8,kind=real_8)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tgc) THEN
       IF (cntl%tlsd) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          ENDDO
          CALL graden(rhoe(1,1),v,grad(1,1),vtmp)
          CALL graden(rhoe(1,2),v,grad(1,5),vtmp)
       ELSE
          CALL graden(rhoe,v,grad,vtmp)
       ENDIF
       IF (cntl%tlsd) THEN
          CALL gclsd(sgcx,sgcc,rhoe,v,vtemp,rvtmp,grad,.FALSE.)
       ELSE
          CALL gcener(sgcx,sgcc,rhoe,v,vtemp,rvtmp(:,1),grad,.FALSE.)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(I,IR)
    DO i=1,clsd%nlsd
       DO ir=1,fpar%nnr1
          fhk(ir,i) = REAL(v(ir,i))
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(v,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rvtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tgc) THEN
       DEALLOCATE(grad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('       VHK',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vhk
  ! ==================================================================
  SUBROUTINE give_scr_vhk(lvhk,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvhk
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ltgc, ltinlc

    ltinlc=0
    IF (cntl%tgc) THEN
       ltgc=fpar%nnr1*clsd%nlsd*4      ! GRADEN
    ELSE
       ltgc=0
    ENDIF
    lvhk = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx+ltgc+ltinlc ! GCENER
    lvhk = lvhk + 2*maxfft*clsd%nlsd + 2*ncpw%nhg
    lvhk = lvhk + 100  ! For boundary checks in SCRPTR
    tag='MAX(NNR1,NHG*2)*NLSX+NNR1*...'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vhk
  ! ==================================================================

END MODULE vhk_utils
