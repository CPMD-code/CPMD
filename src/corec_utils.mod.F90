MODULE corec_utils
  USE cppt,                            ONLY: indz,&
                                             inyh,&
                                             nzh
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel,&
                                             rhoc
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: corec

CONTAINS

  ! ==================================================================
  SUBROUTINE corec(rhoe,vtmp,v)
    ! ==--------------------------------------------------------------==
    ! == Calculates Core charges in G-space and total one in R space  ==
    ! == RHOE:  in valence charge densities in real space             ==
    ! ==       out total charge densities in real space               ==
    ! == VTMP: out core charges in G space                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: vtmp(ncpw%nhg), v(maxfftn)

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: ia, ig, ir, is, isa, isa0, &
                                                isub
    REAL(real_8)                             :: vvre

    CALL tiset('     COREC',isub)
    CALL setfftn(0)
    ! SUM CORE CHARGES IN G-SPACE
    CALL zeroing(vtmp)!,nhg)
    IF (cntl%bigmem) THEN
       isa0=0
       DO is=1,ions1%nsp
          IF (corel%tnlcc(is)) THEN
             DO ia=1,ions0%na(is)
                isa=isa0+ia
#ifdef __SR8000
                !poption parallel, tlocal(IG)
#endif 
                DO ig=1,ncpw%nhg
                   vtmp(ig)=vtmp(ig)+rhoc(ig,is)*eigrb(ig,isa)
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
    ELSE
       isa0=0
       DO is=1,ions1%nsp
          IF (corel%tnlcc(is)) THEN
             DO ia=1,ions0%na(is)
                isa=isa0+ia
#ifdef __SR8000
                !poption parallel, tlocal(IG,EI123)
#endif 
                DO ig=1,ncpw%nhg
                   ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                        ei3(isa,inyh(3,ig))
                   vtmp(ig)=vtmp(ig)+rhoc(ig,is)*ei123
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDIF
    ! FFT TO R-SPACE
    CALL zeroing(v)!,maxfft)
#ifdef __SR8000
    !poption parallel, tlocal(IG)
#endif 
    DO ig=1,ncpw%nhg
       v(nzh(ig))=vtmp(ig)
       v(indz(ig))=CONJG(vtmp(ig))
    ENDDO
    IF (geq0) v(nzh(1))=vtmp(1)
    CALL invfftn(v,.FALSE.,parai%allgrp)
    ! ADD UP THE TOTAL CHARGE IN R-SPACE
    IF (cntl%tlsd) THEN
       ! Remember: RHOE(*,1) is the total charge,
       ! RHOE(*,2) the beta charge
       !$omp  parallel do private(IR,VVRE)
       DO ir=1,fpar%nnr1
          vvre=REAL(v(ir))
          rhoe(ir,1)=rhoe(ir,1)+vvre
          rhoe(ir,2)=rhoe(ir,2)+0.5_real_8*vvre
       ENDDO
    ELSE
       !$omp  parallel do private(IR)
       DO ir=1,fpar%nnr1
          rhoe(ir,1)=rhoe(ir,1)+REAL(v(ir))
       ENDDO
    ENDIF
    CALL tihalt('     COREC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE corec
  ! ==================================================================

END MODULE corec_utils
