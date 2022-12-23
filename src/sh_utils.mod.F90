MODULE sh_utils
  USE cnst,                            ONLY: factem
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       c_sh, dt_sh, shlct, shsigma, shsigmaf, step_diff, td01, tempp_sh, tshf, tshl, &
       v_sh, xfmqc
  USE machine,                         ONLY: m_cputime
  USE mm_dimmod,                       ONLY: mmdim,&
                                             naq
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE soc_types,                       ONLY: do_soc_because_sh
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: propagate_c
  !public :: deriv_sh
  !public :: rk4
  PUBLIC :: do_tullytest
  !public :: check_fermi_golden_rule
  !public :: sh_photon_ex

CONTAINS

  ! ==================================================================
  SUBROUTINE propagate_c(nstate,neltra,shprob)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, neltra
    REAL(real_8)                             :: shprob(neltra,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'propagate_c'

    COMPLEX(real_8), ALLOCATABLE             :: caux(:), cdot(:), cnew(:)
    INTEGER                                  :: ierr, it, iv
    REAL(real_8)                             :: t1, t2
    REAL(real_8), ALLOCATABLE                :: faux(:,:), rdaux(:,:), vaux(:)

    ALLOCATE(rdaux(neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(faux(neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vaux(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(caux(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cnew(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cdot(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(rdaux)!,neltra*neltra)
    CALL zeroing(faux)!,neltra*neltra)
    CALL zeroing(vaux)!,neltra)
    CALL zeroing(caux)!,neltra)
    CALL zeroing(cnew)!,neltra)
    CALL zeroing(cdot)!,neltra)
    ! propagate coefficients
    IF (paral%io_parent)&
         WRITE(6,'(/," PROPAGATION OF COEFFICIENTS")')
    DO it=1,neltra
       caux(it)=c_sh(2,it)
       vaux(it)=v_sh(2,it)
       DO iv=1,neltra
          rdaux(it,iv)=shsigma(2,it,iv)
          faux(it,iv)=shsigmaf(2,it,iv)
       ENDDO
    ENDDO
    CALL deriv_sh(caux,cdot,rdaux,vaux,faux,neltra)
    ![ exact factorization
    IF (tshl%txfmqc) CALL deriv_sh_xf(caux,cdot,neltra)
    !] exact factorization
    t1=m_cputime()
    CALL rk4(caux,neltra,cdot,cnew,shprob)
    t2=m_cputime()
    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F9.3)')&
            " TIME FOR RUNGE-KUTTA INTEGRATION [s]",(t2-t1)*0.001_real_8
    ENDIF
    DO it=1,neltra
       c_sh(1,it)=cnew(it)
    ENDDO
    !
    DEALLOCATE(faux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rdaux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vaux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(caux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cdot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cnew,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE propagate_c
  ! ==================================================================
  SUBROUTINE deriv_sh(caux,cdot,rdaux,vaux,faux,neltra)
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: neltra
    REAL(real_8)                             :: vaux(neltra), &
                                                rdaux(neltra,neltra), &
                                                faux(neltra,neltra)
    COMPLEX(real_8)                          :: cdot(neltra), caux(neltra)

    INTEGER                                  :: it, iv

    CALL zeroing(cdot)!,neltra)
    ! 
    DO it=1,neltra
       cdot(it)=(0.0_real_8,-1.0_real_8)*vaux(it)*caux(it)
    ENDDO
    DO it=1,neltra
       DO iv=1,neltra
          cdot(it)=cdot(it)-caux(iv)*rdaux(it,iv)
       ENDDO
    ENDDO
    DO it=1,neltra
      DO iv=1,neltra
         cdot(it)=cdot(it)+(0.0_real_8,1.0_real_8)*caux(iv)*faux(it,iv)
      ENDDO
    ENDDO
    ! 
    RETURN
  END SUBROUTINE deriv_sh
  ! ==================================================================
  SUBROUTINE deriv_sh_xf(caux,cdot,neltra)
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: neltra
    COMPLEX(real_8)                          :: cdot(neltra), caux(neltra)

    INTEGER                                  :: l

! eq. 19

    DO l=1,neltra
       cdot(l)=cdot(l) + xfmqc%cf(l)
    ENDDO
    ! 
    RETURN
  END SUBROUTINE deriv_sh_xf
  ! ==================================================================
  SUBROUTINE rk4(cold,neltra,dydx,cnew,shprob)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: neltra
    COMPLEX(real_8)                          :: cold(neltra), dydx(neltra), &
                                                cnew(neltra)
    REAL(real_8)                             :: shprob(neltra,neltra)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rk4'

    COMPLEX(real_8), ALLOCATABLE             :: a_matrix(:,:), dym(:), dyt(:)
    INTEGER                                  :: ia, ierr, it, iv, ix, nspline
    REAL(real_8)                             :: delta_fd(neltra,neltra), &
                                                delta_rd(neltra,neltra), &
                                                delta_v(neltra), h6, hh, &
                                                rdaux(neltra,neltra)
    REAL(real_8), ALLOCATABLE                :: faux(:,:), fdata(:,:,:), ffdata(:,:,:), vaux(:), &
                                                vdata(:,:)

! Variables
! interpolation parameters----------------
! 

    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'FOURTH ORDER RUNGE-KUTTA INTEGRATION'
       IF (paral%io_parent)&
            WRITE(6,'(A,E13.5)') ' integration step:', step_diff
    ENDIF
    nspline=NINT(1.0/step_diff)*2
    hh=step_diff*cntr%delt_ions/2._real_8
    h6=step_diff*cntr%delt_ions/6._real_8
    ! 
    ALLOCATE(faux(neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dyt(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dym(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vaux(neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fdata(NINT(1.0/step_diff)*2+1,neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ffdata(NINT(1.0/step_diff)*2+1,neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vdata(NINT(1.0/step_diff)*2+1,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dyt)!,neltra)
    CALL zeroing(dym)!,neltra)
    CALL zeroing(vaux)!,neltra)
    CALL zeroing(vdata)!,(nspline+1)*neltra)
    CALL zeroing(fdata)!,(nspline+1)*neltra*neltra)
    CALL zeroing(ffdata)!,(nspline+1)*neltra*neltra)
    ALLOCATE(a_matrix(neltra,neltra),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(a_matrix)!,neltra*neltra)
    CALL zeroing(shprob)!,neltra*neltra)
    ! 
    ! begin linear interpolation of shsigma and v_sh
    DO it=1,neltra
       DO iv=1,neltra
          delta_rd(it,iv)=0.5*step_diff*(shsigma(1,it,iv)-&
               shsigma(2,it,iv))
          delta_fd(it,iv)=0.5_real_8*step_diff*(shsigmaf(1,it,iv)-&
               shsigmaf(2,it,iv))
       ENDDO
       delta_v(it)=0.5*step_diff*(v_sh(1,it)-v_sh(2,it))
    ENDDO
    DO ia = 1,nspline+1
       DO it=1,neltra
          DO iv=1,neltra
             fdata(ia,it,iv)=shsigma(2,it,iv)+(ia-1)*delta_rd(it,iv)
             ffdata(ia,it,iv)=shsigmaf(2,it,iv)+(ia-1)*delta_fd(it,iv)
          ENDDO
          vdata(ia,it)=v_sh(2,it)+(ia-1)*delta_v(it)
       ENDDO
    ENDDO
    ! end interpolation
    ! 
    ! Comment on the interpolation:
    ! the time evolution occurs between odd interpolation points:
    ! 1 to 3 to 5 to ... each of length h.
    ! the values of "shsigma" and "v_sh" at the even interpolation points
    ! 2, 4, 6, ... are used to compute the derivatives (deriv_sh) at 
    ! the mid-step (time hh)
    ! 
    DO ix=1,nspline-1,2
       ! 
       ! The RK4 STEPS -------------------------------------------------
       ! update coefficients and parameters to the mid-step 
       ! (even spline index)
       ! propagation time hh
       DO it=1,neltra
          cnew(it)=cold(it)+hh*dydx(it)
          vaux(it)=vdata(ix+1,it)
          DO iv=1,neltra
             rdaux(it,iv)=fdata(ix+1,it,iv)
             faux(it,iv)=ffdata(ix+1,it,iv)
          ENDDO
       ENDDO
       ! ... and compute the new derivative
       CALL deriv_sh(cnew,dyt,rdaux,vaux,faux,neltra)
       ! 
       ! Re-update the coefficients with the new derivative
       ! propagation time hh
       DO it=1,neltra
          cnew(it)=cold(it)+hh*dyt(it)
       ENDDO
       ! ... and compute the new derivative
       CALL deriv_sh(cnew,dym,rdaux,vaux,faux,neltra)
       ! 
       ! Update coefficients and parameters to the full-step: h
       DO it=1,neltra
          cnew(it)=cold(it)+2._real_8*hh*dym(it)
          vaux(it)=vdata(ix+2,it)
          DO iv=1,neltra
             rdaux(it,iv)=fdata(ix+2,it,iv)
             faux(it,iv)=ffdata(ix+2,it,iv)
          ENDDO
          ! compute the sum of the derivatives
          dym(it)=dyt(it)+dym(it)
       ENDDO
       ! Compute the new derivatives at the full step
       CALL deriv_sh(cnew,dyt,rdaux,vaux,faux,neltra)
       ! The funal coefficients are
       DO it=1,neltra
          cnew(it)=cold(it)+h6*(dydx(it)+dyt(it)+2._real_8*dym(it))
       ENDDO
       ! The RK4 STEPS -------------------------------------------------
       ! 
       ! Prepare for the next itaration 
       CALL deriv_sh(cnew,dydx,rdaux,vaux,faux,neltra)
       DO it=1,neltra
          cold(it)=cnew(it)
       ENDDO
       ! 
       ! collect quantities for the transition probability
       DO it=1,neltra
          a_matrix(it,td01%fstate+1)=cnew(it)*CONJG(cnew(td01%fstate+1))
       ENDDO
       DO it=1,neltra
          shprob(it,td01%fstate+1)=shprob(it,td01%fstate+1) &
          -2._real_8*(REAL(a_matrix(it,td01%fstate+1),KIND=real_8)*&
                          rdaux(it,td01%fstate+1))&
          /REAL(a_matrix(td01%fstate+1,td01%fstate+1),KIND=real_8) &
         +2._real_8*(DIMAG(a_matrix(it,td01%fstate+1))*&
          faux(it,td01%fstate+1))&
          /REAL(a_matrix(td01%fstate+1,td01%fstate+1),KIND=real_8)
       ENDDO
    ENDDO
    ! 
    ! compute the transition probability
    DO it=1,neltra
       shprob(it,td01%fstate+1)=&
            -2._real_8 * shprob(it,td01%fstate+1) * (cntr%delt_ions * step_diff)
       ! 
    ENDDO
    !
    DEALLOCATE(a_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dyt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dym,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(faux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vaux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fdata,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ffdata,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vdata,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==     
    RETURN
  END SUBROUTINE rk4
  ! ==================================================================
  SUBROUTINE do_tullytest(neltra,shprob)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: neltra
    REAL(real_8)                             :: shprob(neltra,*)

    INTEGER                                  :: cqmd, ia, is, it
    LOGICAL                                  :: debug
    REAL(real_8)                             :: dene, proba_max, proba_min, &
                                                qmdoff, test

    debug=.TRUE.
    proba_max=0._real_8
    proba_min=0._real_8
    test=repprngu()

    IF (paral%parent)THEN
       DO it=1,neltra
          IF (it.NE.td01%fstate+1)THEN
             IF (.NOT.lqmmm%qmmm)THEN
                dt_sh=(v_sh(1,td01%fstate+1)-v_sh(1,it))*factem*2._real_8/glib
             ELSE
                cqmd=0
                DO is=1,mmdim%nspq
                   DO ia=1,NAq(is)
                      cqmd=cqmd+1
                   ENDDO
                ENDDO
                qmdoff = 3*REAL(cqmd,kind=real_8)
                dt_sh=(v_sh(1,td01%fstate+1)-v_sh(1,it))*factem*2._real_8/qmdoff
             ENDIF
             proba_min=proba_max
             IF (shprob(it,td01%fstate+1).LT.0)THEN
                ! write(6,*) 
                ! &        'neg. probability',it,shprob(it,fstate+1),'set to 0'
                shprob(it,td01%fstate+1)=0._real_8
             ENDIF
             test=repprngu()
             proba_max=proba_max+shprob(it,td01%fstate+1)
             IF (debug) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(/,A,I3,A,I3,A)')&
                     ' PROBABILTY OF A SURFACE HOP FROM STATE',&
                     td01%fstate,' TO',it-1,':'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F10.6,A,F10.6,A)')&
                     '   Interval      : [',proba_min,',',proba_max,']'
                IF (paral%io_parent)&
                     WRITE(6,'(A,F10.6)')&
                     '   Diff          :', proba_max-proba_min
                IF (paral%io_parent)&
                     WRITE(6,'(A,F10.6,/)') '   Random number :',test
             ENDIF
             ! if random number in the right interval
             IF (tshl%tdextpot.OR.shlct%sh_lcontrol) THEN
                IF (test.LT.proba_max.AND.test.GT.proba_min)THEN
                   ! the external field must contain a quanta of the correct frequency
                   dene=v_sh(1,it)-v_sh(1,td01%fstate+1)
                   ! call check_fermi_golden_rule(dene,tfg)
                   CALL sh_photon_ex(dene)
                   ! if (tfg) then
                   IF (paral%io_parent)&
                        WRITE(6,'(/,A,I4,A,I4,/)')&
                        ' SH TRANSITION FROM STATE ',td01%fstate,'TO STATE ',it-1
                   td01%fstate=it-1
                   tshl%tully_sh=.TRUE.
                   IF (it.NE.1)THEN
                      tshl%s0_sh=.FALSE.
                   ELSE
                      tshl%s0_sh=.TRUE.
                   ENDIF
                   ! no rescaling of the velocities in this case
                   dt_sh=0._real_8
                   GOTO 100
                   ! endif
                ENDIF
             ELSE
                IF (test.LT.proba_max.AND.test.GT.proba_min)THEN
                   IF (dt_sh.GT.tempp_sh) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(/,A,I4,A,I4,A,/)')&
                           ' SURFACE HOP FROM ',td01%fstate,' TO ',it-1,' IS FRUSTRATED'
                   ELSEIF ((it-1).GT.td01%fstate) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(/,A,I4,A,I4,A,/)')&
                           ' SURFACE HOP FROM ',td01%fstate,' TO ',it-1,' IS FRUSTRATED'
                      ! NACVs are needed for the computation of the frustrated hops.
                      ! (S. Hammes-Schiffer, J.C. Tully, J.Chem.Phys. 101, 4657 (1994))
                      ! NACVs not yet implemented in the TSH dynamics.
                   ELSE
                      IF (paral%io_parent)&
                           WRITE(6,'(/,A,I4,A,I4,/)')&
                           ' SH TRANSITION FROM STATE ',td01%fstate,'TO STATE ',it-1
                      td01%fstate=it-1
                      tshl%tully_sh=.TRUE.
                      IF (it.NE.1)THEN
                         tshl%s0_sh=.FALSE.
                      ELSE
                         tshl%s0_sh=.TRUE.
                      ENDIF
                      do_soc_because_sh=.TRUE.
                      GOTO 100
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDIF
100 CONTINUE
    CALL mp_bcast(td01%fstate,parai%source,parai%allgrp)
    CALL mp_bcast(tshl%tully_sh,parai%source,parai%allgrp)
    CALL mp_bcast(tshl%s0_sh,parai%source,parai%allgrp)
    CALL mp_bcast(do_soc_because_sh,parai%source,parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,'(/,A,I5,/)') "RUN ON SURFACE ", td01%fstate
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE do_tullytest
  ! ==================================================================
  SUBROUTINE check_fermi_golden_rule(dene,tfg)
    REAL(real_8)                             :: dene
    LOGICAL                                  :: tfg

    tfg=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE check_fermi_golden_rule
  ! ==================================================================
  SUBROUTINE sh_photon_ex(dene)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: dene

    IF ((dene.GT.(tshf%afreq-0.1*tshf%afreq)).AND.&
         (dene.LE.(tshf%afreq+0.1*tshf%afreq))) THEN
       tshl%sh_phex=.TRUE.
    ELSE
       tshl%sh_phex=.FALSE.
    ENDIF
    ! 
    IF (paral%io_parent)&
         WRITE(6,*) 'WARNING'
    IF (paral%io_parent)&
         WRITE(6,*) 'sh_phex SETTED to TRUE in sh_util.F line 394'
    tshl%sh_phex=.TRUE.
    ! ==--------------------------------------------------------------==      
    RETURN
  END SUBROUTINE sh_photon_ex
  ! ==================================================================


  ! ==================================================================

END MODULE sh_utils
