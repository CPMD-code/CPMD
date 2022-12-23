#include "cpmd_global.h"

MODULE pcgrad_driver
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_c,&
                                             ener_com,&
                                             ener_d
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE func,                            ONLY: func1
  USE hubbardu,                        ONLY: hubbu
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: lqmmm
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE norm,                            ONLY: cnorm
  USE nvtx_utils
  USE ortho_utils,                     ONLY: ortho,&
                                             preortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rscpot_utils,                    ONLY: rscpot
  USE spin,                            ONLY: lspin2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE vdwcmod,                         ONLY: vdwr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pcgrad
  PUBLIC :: precon
  PUBLIC :: eback

CONTAINS

  ! ==================================================================
  SUBROUTINE pcgrad(c0,c2,sc0,vpp,hnm1,rhoe,psi,tau0,nstate,dinit)
    ! ==--------------------------------------------------------------==
    ! ==  PRECONDITIONED CONJUGATE GRADIENT OPTIMIZATION              ==
    ! ==--------------------------------------------------------------==
    ! == IMPROVED VERSIONS (October 2003)                             ==
    ! ==--------------------------------------------------------------==
    ! == In logscale, we should have a straight line if we plot CNORM ==
    ! == GEMAX versus the number of iterations                        ==
    ! == When CNORM is less than 1.e-5 (or 1.e-6), the error in the   ==
    ! == determination of minimum from energy and from forces are not ==
    ! == negligible.                                                  ==
    ! == As we want decrease CNORM, instead of evaluation             ==
    ! == twice the energy (E1, E2), we evaluate once the force and    ==
    ! == energy of a point (F1,E1) and determnie the minimum along    ==
    ! == the line only using forces (F0 and F1).                      ==
    ! == So we decrease as we want CNORM.                             ==
    ! == An evaluation of forces is three times more cpu-time more    ==
    ! == consuming than an evaluation of energy, so we switch to this ==
    ! == method only when CNORM is small.                             ==
    ! ==--------------------------------------------------------------==
    ! == Some tricks:                                                 ==
    ! == The scalar product <C2.P> is evaluated only for the occupied ==
    ! == states which is coherent with energy.                        ==
    ! ==                                                              ==
    ! == We add a FORGET parameter to forget a little bit.            ==
    ! == This stabilizes the algorithm but it is not crucial.         ==
    ! == At the beginning (ISTEP < 10) we want to forget a lot        ==
    ! == because the steepest descent works fine.                     ==
    ! ==                                                              ==
    ! == We tested to refine the line minimization search if ALAM is  ==
    ! == too big but there is no real improvement.                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: vpp(:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: hnm1(ncpw%ngw,nstate), &
                                                sc0(ncpw%ngw,nstate)
    LOGICAL                                  :: dinit

    CHARACTER(*), PARAMETER                  :: procedureN = 'pcgrad'
    INTEGER, PARAMETER                       :: maxhist = 10 

    INTEGER                                  :: i, ierr, ig, ihist, isub, nocc
    INTEGER, SAVE                            :: ilsr = 0, istep = 0
    LOGICAL                                  :: oldstatus, statusdummy
    REAL(real_8)                             :: a, alam, de, detot, forget, &
                                                gamma, ggnorm
    REAL(real_8), ALLOCATABLE                :: tscr(:,:,:)
    REAL(real_8), SAVE                       :: ehist(maxhist), &
                                                fhist(maxhist), ghist(maxhist)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    ! ==--------------------------------------------------------------==
    ! we need to 'go mm' here so that the TSCR offset is correct.
    IF (lqmmm%qmmm) CALL mm_dim(mm_go_mm,oldstatus)
    ! ==--------------------------------------------------------------==
    ! TODO align for BG
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm) CALL mm_dim(mm_go_qm,statusdummy)
    ! ==--------------------------------------------------------------==
    nocc=0
    DO i=1,nstate
       IF (crge%f(i,1).GT.1.e-5_real_8) THEN
          nocc=nocc+1
       ENDIF
    ENDDO
    IF (dinit.OR.ilsr.LT.0) THEN
       ! We start a new direction.
       istep=0
    ENDIF
    istep=istep+1
    IF (cntl%prec.AND.ilsr.NE.-2) THEN
       CALL precon(c2,vpp,crge%f,ncpw%ngw,nstate,0)
    ENDIF
    ggnorm=0.0_real_8
#if ! defined(__VECTOR)
    !$omp parallel do private(I) reduction(+:GGNORM) schedule(static)
#endif
    DO i=1,nocc
       ggnorm=ggnorm+dotp(ncpw%ngw,c2(:,i),c2(:,i))
    ENDDO
    CALL mp_sum(ggnorm,parai%allgrp)
    CALL mp_bcast(ggnorm,parai%io_source,parai%cp_grp)
    IF (istep.LE.1) THEN
       ! ==--------------------------------------------------------------==
       ! ==  INITIAL CALL                                                ==
       ! ==--------------------------------------------------------------==
       ehist(1)=ener_com%etot
       fhist(1)=dt2bye
       ghist(1)=ggnorm
       IF (cntl%pcgmin) THEN
          CALL dcopy(2*ncpw%ngw*nstate,c2,1,hnm1,1)
          alam=1.0_real_8
          IF (cntl%prec.AND.ilsr.NE.-2) THEN
             ilsr=2
          ELSE
             ilsr=1
          ENDIF
          CALL linesr(nstate,c0,hnm1,c2,sc0,vpp,tau0,tscr,&
               rhoe,psi,&
               fhist(1),de,a,alam,ilsr)
          IF (paral%io_parent) THEN
             WRITE(6,'(A,G8.3,A,T50,G20.13)')&
                  ' LINE SEARCH : LAMBDA=',alam*fhist(1),&
                  ' PREDICTED ENERGY =',dE
          ENDIF
       ELSE
          CALL dcopy(2*ncpw%ngw*nstate,c2,1,hnm1,1)
          CALL daxpy(2*ncpw%ngw*nstate,fhist(1),c2,1,c0,1)
       ENDIF
    ELSE
       ! ==--------------------------------------------------------------==
       ! ==  CONJUGATE GRADIENT                                          ==
       ! ==--------------------------------------------------------------==
       IF (istep.LE.10) THEN
          ! Almost steepest descent
          forget=0.60_real_8
       ELSEIF (istep.LE.20) THEN
          forget=0.91_real_8
       ELSE
          forget=0.95_real_8
       ENDIF
       IF (istep.GT.maxhist) THEN
          ihist=maxhist
          DO i=1,maxhist-1
             ehist(i)=ehist(i+1)
             fhist(i)=fhist(i+1)
             ghist(i)=ghist(i+1)
          ENDDO
       ELSE
          ihist=istep
       ENDIF
       ehist(ihist)=ener_com%etot
       ghist(ihist)=ggnorm
       detot=ehist(ihist)-ehist(ihist-1)
       fhist(ihist)=fhist(ihist-1)
       gamma=ghist(ihist)/ghist(ihist-1)

       ! Multiply by a forget factor
       gamma=forget*gamma

       IF (cntl%pcgmin) THEN
          !$omp parallel do private(I,IG) __COLLAPSE2
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                hnm1(ig,i)=c2(ig,i)+hnm1(ig,i)*gamma
             ENDDO
          ENDDO
          IF (cntl%prec) THEN
             ilsr=4
          ELSE
             ilsr=3
          ENDIF
          CALL linesr(nstate,c0,hnm1,c2,sc0,vpp,tau0,tscr,&
               rhoe,psi,&
               fhist(ihist),de,a,alam,ilsr)
          IF (paral%io_parent) THEN
             WRITE(6,'(A,G8.3,A,T50,G20.13)')&
                  ' LINE SEARCH : LAMBDA=',alam*fhist(ihist),&
                  ' PREDICTED ENERGY =',dE
          ENDIF
       ELSE
          CALL daxpy(2*ncpw%ngw*nstate,gamma,hnm1(1,1),1,c2(1,1),1)
          CALL dcopy(2*ncpw%ngw*nstate,c2(1,1),1,hnm1(1,1),1)
          CALL daxpy(2*ncpw%ngw*nstate,fhist(ihist),hnm1(1,1),1,c0(1,1),1)
       ENDIF
    ENDIF
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,oldstatus)
    DEALLOCATE(tscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pcgrad
  ! ==================================================================
  SUBROUTINE precon(c2,vpp,f,ngw,nstate,mode)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngw
    REAL(real_8)                             :: vpp(ngw)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ngw,nstate)
    REAL(real_8)                             :: f(nstate)
    INTEGER                                  :: mode

    INTEGER                                  :: i, j
    REAL(real_8)                             :: ff, fi(nstate)

! Variables
! ==--------------------------------------------------------------==

#if defined(__VECTOR)
    IF (mode.EQ.0) THEN
       !$omp parallel do private(I,J,FF)
       DO i=1,nstate
          IF (f(i).GT.0.01_real_8) THEN
             ff=1.0_real_8/f(i)
          ELSE
             ff=1.0_real_8
          ENDIF
          DO j=1,ngw
             c2(j,i)=ff*vpp(j)*c2(j,i)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(I,J,FF)
       DO i=1,nstate
          IF (f(i).GT.0.01_real_8) THEN
             ff=1.0_real_8/f(i)
          ELSE
             ff=1.0_real_8
          ENDIF
          DO j=1,ngw
             c2(j,i)=c2(j,i)/(ff*vpp(j))
          ENDDO
       ENDDO
    ENDIF
#else
    !$omp parallel do private(I,J) schedule(static)
    DO i=1,nstate
       IF (f(i).GT.0.01_real_8) THEN
          fi(i)=1.0_real_8/f(i)
       ELSE
          fi(i)=1.0_real_8
       ENDIF
    ENDDO

    IF (mode.EQ.0) THEN
       !$omp parallel do private(I,J) schedule(static) __COLLAPSE2
       DO j=1,ngw
          DO i=1,nstate
             c2(j,i)=fi(i)*vpp(j)*c2(j,i)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(I,J) schedule(static) __COLLAPSE2
       DO j=1,ngw
          DO i=1,nstate
             c2(j,i)=c2(j,i)/(fi(i)*vpp(j))
          ENDDO
       ENDDO
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE precon
  ! ==================================================================
  SUBROUTINE linesr(nstate,x,p,c2,sc0,vpp,tau0,tscr,&
       rhoe,psi,fac,de,a,alam,ilsr)
    ! ==--------------------------------------------------------------==
    ! == ILSR (input/output)                                          ==
    ! ==  1   P is only simple forces without preconditioning         ==
    ! ==  2   P is the force with preconditioning                     ==
    ! ==  3   P is the conjugate gradient scheme without precond.     ==
    ! ==  4   P is the conjugate gradient scheme with preconditioning ==
    ! ==--------------------------------------------------------------==
    ! ==  0   Everything is correct                                   ==
    ! == -1   Reset the direction                                     ==
    ! == -2   Reset the direction without preconditioning             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: x(ncpw%ngw,nstate), &
                                                p(ncpw%ngw,nstate), &
                                                c2(ncpw%ngw,nstate), &
                                                sc0(ncpw%ngw,nstate)
    REAL(real_8)                             :: vpp(*), tau0(:,:,:), &
                                                tscr(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: fac, de, a, alam
    INTEGER                                  :: ilsr

    INTEGER                                  :: i, nfunc
    LOGICAL                                  :: oldstatus, statusdummy
    REAL(real_8)                             :: alam0, b, c, dd, dx, e0, e1, &
                                                e2, ee, eigv(1), force0, &
                                                force1, x0

! ==--------------------------------------------------------------==

    IF (lqmmm%qmmm) CALL mm_dim(mm_go_mm,oldstatus)
    IF (lqmmm%qmmm) CALL mm_dim(mm_go_qm,statusdummy)
    ! ==--------------------------------------------------------------==
    CALL eback(0)
    nfunc=0
    e0=ener_com%etot
    CALL dscal(2*ncpw%ngw*nstate,fac,p,1)
    ! AK  XETOT does not account for QM/MM, so we have to always do it expensively.
    ! VW  ... for HFX.
    IF (cnorm.LT.1.0e-5_real_8.OR.cntl%tinter.OR.lqmmm%qmmm.OR.func1%mhfx.NE.0) THEN
       ! The norm is small.
       ! To avoid the error between forces and energies
       ! We use forces for each point to interpolate the new one.
       alam=1.0_real_8
       alam0=0._real_8
       ! Calculates the force in the direction HNM1
       force0=0.0_real_8
       DO i=1,nstate
          IF (crge%f(i,1).GT.1.e-5_real_8) THEN
             force0=force0+dotp(ncpw%ngw,p(:,i),c2(:,i))
          ENDIF
       ENDDO
       CALL mp_sum(force0,parai%allgrp)
       CALL mp_bcast(force0,parai%io_source,parai%cp_grp)
       ! First point
       CALL daxpy(2*ncpw%ngw*nstate,alam-alam0,p,1,x,1)
       alam0=alam

       ! Calculate the new point
       IF (lqmmm%qmmm) THEN
          CALL mm_qmmm_forcedr(x,c2,sc0,rhoe,psi,tau0,tscr,&
               eigv,nstate,1,.TRUE.,.FALSE.,.FALSE.)
       ELSE
          CALL forcedr(x,c2,sc0,rhoe,psi,tau0,tscr,eigv,&
               nstate,1,.TRUE.,(cntl%tinter.EQV..TRUE.))
       ENDIF

       IF (cntl%prec.AND.ilsr.NE.-2) THEN
          CALL precon(c2,vpp,crge%f,ncpw%ngw,nstate,0)
       ENDIF
       ! Calculates the force in the direction P
       force1=0.0_real_8
       DO i=1,nstate
          ! Only for the occupied states
          IF (crge%f(i,1).GT.1.e-5_real_8) THEN
             force1=force1+dotp(ncpw%ngw,p(:,i),c2(:,i))
          ENDIF
       ENDDO
       CALL mp_sum(force1,parai%allgrp)
       CALL mp_bcast(force1,parai%io_source,parai%cp_grp)
       e1=ener_com%etot
       ! Determine X0 using E0, F0 and F1 (A.X0**2+B.X0+C)
       a = (force1-force0)/(2._real_8*alam)
       b = force0
       c = e0
       IF (ABS(a).LT.1.0e-15_real_8) a=SIGN(1.0e-15_real_8,a)
       x0=-b/(2.0_real_8*a)
       de = 0.5_real_8*(force1-force0)*x0*x0+force0*x0+e0
       CALL mp_bcast(x0,parai%io_source,parai%cp_grp)
       CALL mp_bcast(alam,parai%io_source,parai%cp_grp)
       CALL mp_bcast(alam0,parai%io_source,parai%cp_grp)
       IF (x0.LT.0._real_8) THEN
          IF (ilsr.EQ.1) THEN
             ! The direction should be correct
             ! Decrease the size step
             alam=0.05_real_8
          ELSE
             ! Reset the direction
             IF (e1.LT.e0) THEN
                alam=1._real_8
             ELSE
                alam=0.05_real_8
             ENDIF
          ENDIF
          ilsr=-2
          ! New point
          CALL daxpy(2*ncpw%ngw*nstate,alam-alam0,p,1,x,1)
       ELSE
          ilsr=0
          alam=x0
          IF (alam.GT.3._real_8) THEN
             alam=3.0_real_8
          ENDIF
          ! New point
          CALL daxpy(2*ncpw%ngw*nstate,alam-alam0,p,1,x,1)
       ENDIF

    ELSE

       alam=1.0_real_8
       ! Standard scheme (low time-consuming)
       dd=0.0_real_8
       DO i=1,nstate
          dd=dd+dotp(ncpw%ngw,p(:,i),p(:,i))
       ENDDO
       CALL mp_sum(dd,parai%allgrp)
       CALL mp_bcast(dd,parai%io_source,parai%cp_grp)
       ! First point
       CALL dcopy(2*ncpw%ngw*nstate,x,1,c2(1,1),1)
       CALL daxpy(2*ncpw%ngw*nstate,alam,p,1,x,1)
       alam0=alam
       ! Calculate the new point
       CALL xetot(nstate,x,sc0,tau0,tscr,rhoe,psi)
       IF (ener_com%etot.GT.e0) THEN
          e2=ener_com%etot
          alam=0.5_real_8
       ELSE
          e1=ener_com%etot
          alam=2.0_real_8
       ENDIF
       ! Second point
       CALL dcopy(2*ncpw%ngw*nstate,c2(1,1),1,x,1)
       CALL daxpy(2*ncpw%ngw*nstate,alam,p,1,x,1)
       alam0=alam
       ! Calculate the new point
       CALL xetot(nstate,x,sc0,tau0,tscr,rhoe,psi)
       IF (alam.GT.1._real_8) THEN
          e2=ener_com%etot
          dx=1.0_real_8
       ELSE
          e1=ener_com%etot
          dx=0.5_real_8
       ENDIF
       ! Extrapolation
       IF ((2.0_real_8*e0-4.0_real_8*e1+2.0_real_8*e2).EQ.0._real_8) THEN
          ! Very rare (T.D.)
          x0=dx
       ELSE
          x0=dx*(4*e1-e2-3*e0)/(2*e0-4*e1+2*e2)
       ENDIF

       IF ((2.0_real_8*x0+dx).EQ.0._real_8) THEN
          a=0._real_8
       ELSE
          a=(e1-e0)/dx/(2.0_real_8*x0+dx)
       ENDIF
       ee=e0-a*x0*x0
       alam=-x0
       IF (alam.GT.0._real_8) THEN
          ! Limit the move
          ilsr=0
          IF (alam.GT.3._real_8) THEN
             ! New minimisation
             alam=3.0_real_8
          ENDIF
       ELSE
          ! ALAM is negative, take the lowest energy between E0, E1 and E2
          IF (e2.LT.e0) THEN
             alam=2.0_real_8*dx
             ilsr=0
          ELSEIF (e1.LT.e0) THEN
             alam=dx
             ilsr=0
          ELSE
             ! E0 < E1 and E2
             alam=0.05_real_8
             ! Reset the direction
             ilsr=-2
          ENDIF
       ENDIF
       de=ee+a*(x0+alam)**2
       a=a/dd
       ! New point
       CALL dcopy(2*ncpw%ngw*nstate,c2(1,1),1,x,1)
       CALL daxpy(2*ncpw%ngw*nstate,alam,p,1,x,1)
    ENDIF
    CALL dscal(2*ncpw%ngw*nstate,1._real_8/fac,p,1)
    IF (alam.LT.0.1_real_8) THEN
       fac=fac*0.25_real_8
    ELSEIF (alam.GT.1.0_real_8) THEN
       fac=fac*1.50_real_8
    ELSEIF (alam.LT.0.2_real_8) THEN
       fac=fac*0.75_real_8
    ENDIF
    CALL eback(1)
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,oldstatus)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE linesr
  ! ==================================================================
  SUBROUTINE xetot(nstate,x,sc0,tau0,tscr,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE THE TOTAL ENERGY FOR A GIVEN WAVEFUNCTION X        ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: x(ncpw%ngw,nstate), &
                                                sc0(ncpw%ngw,nstate)
    REAL(real_8)                             :: tau0(:,:,:), tscr(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    IF (cntl%nonort.AND.(.NOT.cntl%quenchb)) THEN
       CALL stopgm('XETOT','NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSE
       IF (pslo_com%tivan) THEN
          CALL rnlsm(x,nstate,1,1,.FALSE.)
       ENDIF
       CALL preortho(x,nstate)
       CALL ortho(nstate,x,sc0)
       CALL rnlsm(x,nstate,1,1,.FALSE.)
       CALL rscpot(x,tau0,tscr,rhoe,psi,.FALSE.,.FALSE.,nstate,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xetot
  ! ==================================================================
  SUBROUTINE eback(mode)
    ! ==--------------------------------------------------------------==
    ! == Save or restore values of energies                           ==
    ! == Mode 0 -> SAVE                                               ==
    ! ==      1 -> RESTORE                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mode

    REAL(real_8), SAVE :: bbc(4), bcas(21), bcsumg, bcsumr, becnstr, begc, &
      behee, behep, behii, beht, bekin, benl, bepseu, berestr, beself, besr, &
      betot, bevdw, bexc, bvxc, behub

    IF (mode.EQ.0) THEN
       betot=ener_com%etot
       bekin=ener_com%ekin
       bepseu=ener_com%epseu
       benl=ener_com%enl
       beht=ener_com%eht
       behep=ener_com%ehep
       behee=ener_com%ehee
       behii=ener_com%ehii
       bexc=ener_com%exc
       bvxc=ener_com%vxc
       beself=ener_com%eself
       besr=ener_com%esr
       begc=ener_com%egc
       becnstr=ener_com%ecnstr
       berestr=ener_com%erestr
       bevdw=vdwr%evdw
       bcsumg=chrg%csumg
       bcsumr=chrg%csumr
       behub=hubbu%ehub
       IF (lspin2%tlse .AND. lspin2%tcas22) THEN
          CALL dcopy(21,ener_c%etot_a,1,bcas,1)
          CALL dcopy(4,ener_d%etot_b,1,bbc,1)
       ENDIF
    ELSE
       ener_com%etot=betot
       ener_com%ekin=bekin
       ener_com%epseu=bepseu
       ener_com%enl=benl
       ener_com%eht=beht
       ener_com%ehep=behep
       ener_com%ehee=behee
       ener_com%ehii=behii
       ener_com%exc=bexc
       ener_com%vxc=bvxc
       ener_com%eself=beself
       ener_com%esr=besr
       ener_com%egc=begc
       ener_com%ecnstr=becnstr
       ener_com%erestr=berestr
       vdwr%evdw=bevdw
       chrg%csumg=bcsumg
       chrg%csumr=bcsumr
       hubbu%ehub=behub
       IF (lspin2%tlse .AND. lspin2%tcas22) THEN
          CALL dcopy(21,bcas,1,ener_c%etot_a,1)
          CALL dcopy(4,bbc,1,ener_d%etot_b,1)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eback

END MODULE pcgrad_driver
