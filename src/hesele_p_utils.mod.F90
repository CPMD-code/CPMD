MODULE hesele_p_utils
  USE cppt,                            ONLY: hg,&
                                             twnl
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE pslo,                            ONLY: pslo_com
  USE response_pmod,                   ONLY: response1
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hesele_p

CONTAINS

  ! ==================================================================
  SUBROUTINE hesele_p (svar2,z11,nstate,vpp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: svar2
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate), &
                                                vpp(ncpw%ngw)

    INTEGER                                  :: ia, ig, is, isa0, iv, j, jv, &
                                                ki, kj, l, l2, li, lj
    REAL(real_8)                             :: dd, ftpi

! Variables
! ==--------------------------------------------------------------==
! ==  Calculate   H0 : The diagonal approximation to the          ==
! ==                   2nd derivative matrix                      ==
! ==--------------------------------------------------------------==

    IF (response1%prec_p) THEN
       ! ==--------------------------------------------------------------==
       ! ==  Nonlocal Potential                                          ==
       ! ==--------------------------------------------------------------==
       CALL zeroing(vpp)!,ngw)
       DO ig=1,ncpw%ngw
          vpp(ig)=0.0_real_8
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      dd=dvan(iv,jv,is)
                      DO ia=1,ions0%na(is)
                         dd=dd+deeq(isa0+ia,iv,jv,1)
                      ENDDO
                      IF (cntl%tlsd) THEN
                         DO ia=1,ions0%na(is)
                            dd=dd+deeq(isa0+ia,iv,jv,2)
                         ENDDO
                         dd=0.5_real_8*dd
                      ENDIF
                      vpp(ig)=vpp(ig)+dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                   ENDDO
                ENDDO
             ELSEIF (sgpp1%tsgp(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         dd=sgpp2%hlsg(ki,kj,l,is)*ions0%na(is)
                         vpp(ig)=vpp(ig)+dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   vpp(ig)=vpp(ig)+ions0%na(is)*wsg(is,iv)&
                        *twnl(ig,iv,is,1)*twnl(ig,iv,is,1)
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  Kinetic Energy and consolidation                            ==
       ! ==--------------------------------------------------------------==
       ftpi=0.5_real_8*parm%tpiba2
       !$omp parallel do private(ig)
       DO ig=1,ncpw%ngw
          vpp(ig)=vpp(ig)+ftpi*hg(ig)
       ENDDO

       ! ==--------------------------------------------------------------==
    ELSE                      ! IF NO PRECONDITIONING then here:
       ! ==--------------------------------------------------------------==
       !$omp parallel do private(J)
       DO j=1,ncpw%ngw
          vpp(j)=svar2
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hesele_p
  ! ==================================================================

END MODULE hesele_p_utils
