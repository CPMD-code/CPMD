MODULE prpcmove_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com,&
                                             veps
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prpcmove_iso
  PUBLIC :: prpcmove

CONTAINS

  ! ==================================================================
  SUBROUTINE prpcmove_iso(velp,step,pmx)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)

    INTEGER                                  :: ia, is, k
    REAL(real_8)                             :: aa, at3, geps_p, geps_v

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT VELOCITY                                 ==
    ! ==--------------------------------------------------------------==
    geps_p = metr_com%htfp(1,1) + metr_com%htfp(2,2) + metr_com%htfp(3,3)
    geps_v = 0._real_8
    DO k=1,3
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             geps_v = geps_v + (1._real_8 + 3._real_8/at3)*&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    veps = veps + 0.25_real_8*step*(geps_p + geps_v)/cntr%cmass
    ! ==--------------------------------------------------------------==
    ! ==  SCALE VELOCITIES                                            ==
    ! ==--------------------------------------------------------------==
    aa = EXP(-0.5_real_8*step*(1._real_8+3._real_8/at3)*veps)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          velp(1,ia,is) = velp(1,ia,is)*aa
          velp(2,ia,is) = velp(2,ia,is)*aa
          velp(3,ia,is) = velp(3,ia,is)*aa
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT VELOCITY                                 ==
    ! ==--------------------------------------------------------------==
    geps_p = metr_com%htfp(1,1) + metr_com%htfp(2,2) + metr_com%htfp(3,3)
    geps_v = 0._real_8
    DO k=1,3
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             geps_v = geps_v + (1._real_8 + 3._real_8/at3)*&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    veps = veps + 0.25_real_8*step*(geps_p + geps_v)/cntr%cmass
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpcmove_iso
  ! ==================================================================
  SUBROUTINE prpcmove(velp,step,pmx)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)

    INTEGER                                  :: ia, ierr, is, k, l
    REAL(real_8)                             :: aa1, aa2, aa3, at3, &
                                                htfv(3,3), scmat(3,3), &
                                                scmd(3), scmev(3,3), trhv, &
                                                vtmp1(3)

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! ==  GET VELOCITY-DEPENDENT PART OF CELL FORCE                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(htfv)!,9)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,3
                htfv(k,l) = htfv(k,l) +&
                     pmx(is)*velp(k,ia,is)*velp(l,ia,is)
             ENDDO
             htfv(k,k) = htfv(k,k) +&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)/at3
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT (CELL VELOCITY)                          ==
    ! ==--------------------------------------------------------------==
    DO k=1,3
       DO l=1,3
          metr_com%htvel(k,l) = metr_com%htvel(k,l) +&
               0.25_real_8*step*(htfv(k,l) + metr_com%htfp(k,l))/cntr%cmass
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  FOR VELOCITY SCALING, FORM MATRIX v_h + Tr(v_h)/N + v_zeta1 ==
    ! ==--------------------------------------------------------------==
    trhv = metr_com%htvel(1,1) + metr_com%htvel(2,2) + metr_com%htvel(3,3)
    DO k=1,3
       DO l=1,3
          scmat(k,l) = metr_com%htvel(k,l)
       ENDDO
       scmat(k,k) = scmat(k,k) + trhv/at3
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == DIAGONALIZE THIS MATRIX AND SAVE EIGENVALUES AND EIGENVECTORS==
    ! ==--------------------------------------------------------------==
    CALL jacobi(3,3,scmat,scmd,scmev,ierr)
    ! ==--------------------------------------------------------------==
    ! ==  ROTATE, SCALE, AND ROTATE BACK FOR VELOCITIES               ==
    ! ==--------------------------------------------------------------==
    aa1 = EXP(-0.5_real_8*step*scmd(1))
    aa2 = EXP(-0.5_real_8*step*scmd(2))
    aa3 = EXP(-0.5_real_8*step*scmd(3))
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          ! 
          vtmp1(1) = scmev(1,1)*velp(1,ia,is) + scmev(2,1)*velp(2,ia,is)&
               + scmev(3,1)*velp(3,ia,is)
          vtmp1(2) = scmev(1,2)*velp(1,ia,is) + scmev(2,2)*velp(2,ia,is)&
               + scmev(3,2)*velp(3,ia,is)
          vtmp1(3) = scmev(1,3)*velp(1,ia,is) + scmev(2,3)*velp(2,ia,is)&
               + scmev(3,3)*velp(3,ia,is)
          ! 
          vtmp1(1) = vtmp1(1)*aa1
          vtmp1(2) = vtmp1(2)*aa2
          vtmp1(3) = vtmp1(3)*aa3
          ! 
          velp(1,ia,is) = scmev(1,1)*vtmp1(1) + scmev(1,2)*vtmp1(2)&
               + scmev(1,3)*vtmp1(3)
          velp(2,ia,is) = scmev(2,1)*vtmp1(1) + scmev(2,2)*vtmp1(2)&
               + scmev(2,3)*vtmp1(3)
          velp(3,ia,is) = scmev(3,1)*vtmp1(1) + scmev(3,2)*vtmp1(2)&
               + scmev(3,3)*vtmp1(3)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  GET VELOCITY-DEPENDENT PART OF CELL FORCE                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(htfv)!,9)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,3
                htfv(k,l) = htfv(k,l) +&
                     pmx(is)*velp(k,ia,is)*velp(l,ia,is)
             ENDDO
             htfv(k,k) = htfv(k,k) +&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)/at3
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT (CELL VELOCITY)                          ==
    ! ==--------------------------------------------------------------==
    DO k=1,3
       DO l=1,3
          metr_com%htvel(k,l) = metr_com%htvel(k,l) +&
               0.25_real_8*step*(htfv(k,l) + metr_com%htfp(k,l))/cntr%cmass
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpcmove
  ! ==================================================================

END MODULE prpcmove_utils
