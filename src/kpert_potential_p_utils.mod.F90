MODULE kpert_potential_p_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: paral
  USE pslo,                            ONLY: pslo_com
  USE response_pmod,                   ONLY: ddfnl_ddk,&
                                             dfnl_dk,&
                                             dtwnl_dk,&
                                             e2_nl_c0,&
                                             fnl00
  USE sfac,                            ONLY: eigr
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fnl_fkp_p
  PUBLIC :: kin_pert

CONTAINS

  ! ==================================================================
  SUBROUTINE fnl_fkp_p(c0,h1nl,f,nstate,k_,fkpfirst)
    ! ==--------------------------------------------------------------==
    ! == Computes H1NL (from FNL00), the perturbation needed for      ==
    ! == the KPERT perturbation                                       ==
    ! == the derivative of the nonlocal pseudopotential with respect  ==
    ! == to the K_ component of the K vector, applied to c0.          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: h1nl(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: k_
    LOGICAL                                  :: fkpfirst

    COMPLEX(real_8)                          :: ctm, i_to_l
    INTEGER                                  :: ia, iat, ig, isp_, istate, &
                                                isub, iv, jv, k2_, ki, kj, l, &
                                                l2, li, lj
    REAL(real_8)                             :: dd0, dd1

    CALL tiset(' FNL_FKP_P',isub)

    CALL zeroing(h1nl)!,SIZE(h1nl))
    iat = 0
    ! ==--------------------------------------------------------------==
    DO istate = 1,nstate
       DO isp_ = 1,ions1%nsp
          ! ==--------------------------------------------------------------==
          ! ==  VANDERBILT                                                  ==
          ! ==--------------------------------------------------------------==
          IF (pslo_com%tvan(isp_)) THEN
             CALL stopgm('K_PERT_PREP',&
                  'Vanderbilt PP not implemented',& 
                  __LINE__,__FILE__)
          ELSEIF (sgpp1%tsgp(isp_)) THEN
             ! ==--------------------------------------------------------------==
             ! ==  STEFAN GOEDECKER                                            ==
             ! ==--------------------------------------------------------------==
             IF (paral%parent .AND. istate .EQ. 1) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5)')&
                     "STEFAN GOEDECKER PP for species",isp_
             ENDIF
             IF (cntl%tfdist) THEN
                CALL stopgm('K_PERT_PREP',&
                     'distributed DFNLs not supported yet.',& 
                     __LINE__,__FILE__)
             ENDIF
             DO ia=1,ions0%na(isp_)
                iat=iat+1
                DO iv=1,nlps_com%ngh(isp_)
                   l=nghtol(iv,isp_)+1
                   ki=sgpp2%lfval(iv,isp_)
                   li=sgpp2%lpval(iv,isp_)
                   dd0=0.0_real_8
                   dd1 = 0.0_real_8
                   DO jv=1,nlps_com%ngh(isp_)
                      l2=nghtol(jv,isp_)+1
                      lj=sgpp2%lpval(jv,isp_)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,isp_)
                         dd0=dd0+fnl00(1,iat,jv,istate,1)&
                              *sgpp2%hlsg(ki,kj,l,isp_)
                         dd1=dd1+&
                              dfnl_dk(2,iat,jv,istate,k_)&
                              *sgpp2%hlsg(ki,kj,l,isp_)
                      ENDIF
                   ENDDO
                   ctm=-f(istate)*(0.0_real_8,-1.0_real_8)**nghtol(iv,isp_)
                   DO ig=1,ncpw%ngw
                      h1nl(ig,istate)=h1nl(ig,istate)+&
                           ctm * eigr(ig,iat,1)*&
                           (twnl(ig,iv,isp_,1) * dd1 +&
                           (0.0_real_8,-1.0_real_8)*dtwnl_dk(ig,iv,isp_,k_)*&
                           dd0)
                   ENDDO

                   ! E2_NL_C0 quadratic terms depending only on the C0 wavefunctions
                   IF (fkpfirst) THEN
                      DO k2_ = 1,3
                         e2_nl_c0(k_,k2_)=e2_nl_c0(k_,k2_)+&
                              f(istate)*&
                              (dd0*ddfnl_ddk(1,iat,iv,istate,k_,k2_)+&
                              dd1*dfnl_dk(2,iat,iv,istate,k2_))
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE

             ! ==--------------------------------------------------------------==
             ! ==  BACHELET HAMANN SCHLUTER                                    ==
             ! ==--------------------------------------------------------------==
             DO ia = 1,ions0%na(isp_)
                iat = iat + 1
                DO iv = 1,nlps_com%ngh(isp_)
                   i_to_l = (0.0_real_8,-1.0_real_8)**nghtol(iv,isp_)
                   ctm    = - i_to_l * f(istate)*wsg(isp_,iv)
                   DO ig=1,ncpw%ngw
                      h1nl(ig,istate)=h1nl(ig,istate)+&
                           ctm * eigr(ig,iat,1)*&
                           (twnl(ig,iv,isp_,1) *&
                           dfnl_dk(2,iat,iv,istate,k_)+&
                           (0.0_real_8,-1.0_real_8)*dtwnl_dk(ig,iv,isp_,k_)*&
                           fnl00(1,iat,iv,istate,1))
                   ENDDO

                   ! E2_NL_C0 quadratic terms depending only on the C0 wavefunctions

                   IF (fkpfirst) THEN
                      DO k2_ = 1,3
                         e2_nl_c0(k_,k2_)=e2_nl_c0(k_,k2_)+&
                              f(istate)*wsg(isp_,iv)*&
                              (fnl00(1,iat,iv,istate,1)&
                              *ddfnl_ddk(1,iat,iv,istate,k_,k2_)+&
                              dfnl_dk(2,iat,iv,istate,k_)&
                              *dfnl_dk(2,iat,iv,istate,k2_))
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       iat = 0
    ENDDO

    CALL tihalt(' FNL_FKP_P',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnl_fkp_p
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE kin_pert(c0,h1nl,f,nstate,k_)
    ! ==--------------------------------------------------------------==
    ! == Computes H1NL, the perturbation needed for                   ==
    ! == the k-point perturbation.                                    ==
    ! == the gradient is applied to C0                                ==
    ! == it is done only for one direction k (=xyz)                   ==
    ! == i * |HINL> = i * (-) (d |u0> / d rk)
    ! == H1NL(IG,ISTATE) = (-) F(ISTATE) * (-i * Gk(IG)) * C0(IG,ISTATE)
    ! == the (-) is due to the optimization algorithm 
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: h1nl(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: k_

    COMPLEX(real_8)                          :: g
    INTEGER                                  :: ig, istate

! Variables
! ==--------------------------------------------------------------==
! CALL TISET(' KIN_PERT',ISUB)
! ==--------------------------------------------------------------==
! ==  kinetic part of the perturbation                            ==
! ==--------------------------------------------------------------==
! forces are computed with the negative sign because of 
! the optimization algorithm

    DO istate=1,nstate
       DO ig=1,nkpt%ngwk
          g=CMPLX(0._real_8,-parm%tpiba*Gk(k_,ig),kind=real_8)
          h1nl(ig,istate)=h1nl(ig,istate)-&
               f(istate)*g*c0(ig,istate)
       ENDDO
    ENDDO
    ! CALL TIHALT(' KIN_PERT',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kin_pert
  ! ==================================================================














END MODULE kpert_potential_p_utils
