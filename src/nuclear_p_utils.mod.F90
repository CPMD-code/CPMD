#include "cpmd_global.h"

MODULE nuclear_p_utils
  USE cnst,                            ONLY: pi
  USE coor,                            ONLY: fion,&
                                             tau0
  USE densrd_utils,                    ONLY: densrd
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE printp_utils,                    ONLY: printp
  USE recpnew_utils,                   ONLY: recpnew
  USE response_pmod,                   ONLY: response1,&
                                             response2,&
                                             rho0
  USE rhoofr_p_utils,                  ONLY: rhoofr_p
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rhopri_utils,                    ONLY: rhopri
  USE rinforce_nuc_utils,              ONLY: rinforce_nuc
  USE ropt,                            ONLY: ropt_mod
  USE rwfopt_nuc_utils,                ONLY: rwfopt_nuc
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sgpp,                            ONLY: sgpp2
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cnti,&
                                             fpar,&
                                             maxsp,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE testex_utils,                    ONLY: testex
  USE write_pp_utils,                  ONLY: write_pp
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nuclear_p
  !public :: oacp4forces
  !public :: forcepenalty
  !public :: densitypenalty
  !public :: dummy_optcalc
  !public :: steepestdescent
  !public :: pert4nuc
  !public :: dummy_refcalc
  !public :: distance

CONTAINS

  ! =======================================================================
  SUBROUTINE nuclear_p(c0,c1,psi,rhoe,drhoe,&
       eirop,eivps,z11,nstate)
    ! =======================================================================
    COMPLEX(real_8)                          :: psi(maxfftn,clsd%nlsd)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                drhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), eivps(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,*)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'nuclear_p'

    COMPLEX(real_8)                          :: eirop1(ncpw%nhg)
    COMPLEX(real_8), ALLOCATABLE             :: gde(:,:), pme(:,:), sc0(:), &
                                                vtemp(:)
    INTEGER                                  :: ierr, step
    LOGICAL                                  :: status
    REAL(real_8)                             :: crit_tmp, tp, tp1
    REAL(real_8), ALLOCATABLE                :: eigv(:), rho_ref(:), &
                                                rho_tmp(:,:), vpp(:)

! variables for updwf
! variables
! species, species, optimization step
! ==--------------------------------------------------------------==

    ALLOCATE(pme(nkpt%ngwk,cnti%mdiis*(nkpt%ngwk*nstate+8)/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(nkpt%ngwk,cnti%mdiis*(nkpt%ngwk*nstate+8)/(4*nkpt%ngwk)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv((crge%n*nkpt%nkpts)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nkpt%ngwk*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rho_ref(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rho_tmp(fpar%nnr1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,16("OACP"))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"OACP",56X,"OACP")')
       IF (paral%io_parent)&
            WRITE(6,'(" OACP",5X,A,11X,"OACP")')&
            'OPTIMIZATION of ATOM CENTERED POTENTIALS'
       IF (paral%io_parent)&
            WRITE(6,'(" OACP",7X,A,17X,"OACP")')&
            'J. Chem. Phys. 122, 14113 (2005)'
       IF (paral%io_parent)&
            WRITE(6,'(1X,"OACP",56X,"OACP")')
       IF ((response1%tdummyatom_ref).AND.paral%io_parent)&
            WRITE(6,'(" OACP",5X,A,19X,"OACP")')&
            'Calculation of reference density'
       IF ((response1%tdummyatom).AND.paral%io_parent)&
            WRITE(6,'(" OACP",5X,A,27X,"OACP")')&
            'Density Optimization Run'
       IF ((response1%toacp4f).AND.paral%io_parent)&
            WRITE(6,'(" OACP",5X,A,29X,"OACP")')&
            'Force Optimization Run'
       IF (paral%io_parent)&
            WRITE(6,'(1X,"OACP",56X,"OACP")')
       IF (paral%io_parent)&
            WRITE(6,'(1X,16("OACP"),/)')
    ENDIF

    crit_tmp = 100.0_real_8 ! temporary convergence criterium
    step=0

    ! ==--------------------------------------------------------------==
    IF (response1%tdummyatom_ref)&
         CALL dummy_refcalc(c0,psi,nstate,rhoe)
    ! ==--------------------------------------------------------------==
    IF (response1%tdummyatom_tec) THEN
       ! change to QMMM-dimensions
       CALL mm_dim(mm_go_mm,status)
       ropt_mod%rprint = .TRUE.
       ! write out TRAJECTORY
       IF (paral%parent) CALL printp(tau0,tau0,tau0)
       ! change back to QM-dimensions
       CALL mm_dim(mm_go_qm,status)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (response1%toacp4f) THEN

       CALL oacp4forces(c0,c1,sc0,pme,gde,vpp,eigv,&
            rhoe,psi,nstate)

       ! converged
       IF (paral%parent)THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,16("OACP"))')
          IF (paral%io_parent)&
               WRITE(6,'(" OACP",5X,A,11X,"OACP")')&
               'OPTIMIZATION of ATOM CENTERED POTENTIALS'
          IF (paral%io_parent)&
               WRITE(6,*) 'OACP    - is converged'
          IF (paral%io_parent)&
               WRITE(6,'(1X,16("OACP"),/)')
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (response1%tdummyatom) THEN
       step = 1

       ! read in reference density into rho_ref 
       CALL densrd(vtemp,'rhoofr_ref.dat')
       ! CALL FFTTOR(vtemp,rho_ref,scr,nhg,.true.)
       CALL ffttor(vtemp,rho_ref,psi,ncpw%nhg,.TRUE.)! TODO check scr->psi

       ! get the electrondensity from c0 in rho_tmp
       CALL rhoofr(c0,rho_tmp,psi(:,1),nstate)

       ! compute the penalty value by integration over volume, 
       ! tP = \int dr^3 |rho_ref - rho_tmp|^2
       CALL densitypenalty(rho_ref,rho_tmp,tp1,.TRUE.)
       CALL densitypenalty(rho_ref,rho_tmp,tp,.FALSE.)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*) 'OACP density penalty for reference vol. :',tp1
          IF (paral%io_parent)&
               WRITE(6,*) 'OACP density penalty for full vol.      :',tp
       ENDIF

       ! loop to minimize density penalty tP
       DO WHILE ((response2%dumcrit.LT.crit_tmp).AND..NOT.soft_com%exsoft)
          step=step+1

          ! compute all penalty derivatives, change all parameters 
          ! accordingly, optimize new wavefunction
          CALL dummy_optcalc(c0,c1,psi,nstate,rhoe,drhoe,&
               eirop,eirop1,eivps,z11,rho_ref,rho_tmp,vtemp,crit_tmp,&
               pme,gde,sc0,vpp,eigv)

          ! print out old penalty
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,16("OACP"))')
             IF (paral%io_parent)&
                  WRITE(6,*)'OAPC intermediate results for step:',step
             IF (paral%io_parent)&
                  WRITE(6,*)'OACP temporary convergence criterium:',crit_tmp
             IF (paral%io_parent)&
                  WRITE(6,*)'OACP old density penalty for full vol.     :',tp
             IF (paral%io_parent)&
                  WRITE(6,*)'OACP new density penalty for reference vol.:',tp1
          ENDIF

          ! get the new electrondensity from c0 into rho_tmp
          CALL rhoofr(c0,rho_tmp,psi(:,1),nstate)

          ! compute the penalty value by integration over volume, 
          ! tP = \int dr^3 |rho_ref - rho_tmp|^2
          CALL densitypenalty(rho_ref,rho_tmp,tp1,.TRUE.)
          CALL densitypenalty(rho_ref,rho_tmp,tp,.FALSE.)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'OACP new density penalty for full vol.     :',tp
             IF (paral%io_parent)&
                  WRITE(6,*)'OACP new density penalty for reference vol.:',tp1
             IF (paral%io_parent)&
                  WRITE(6,'(1X,16("OACP"),/)')
          ENDIF

          IF (step.GE.cnti%nomore) soft_com%exsoft=.TRUE.
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)

       ENDDO                ! convergence loop

       ! converged
       IF (paral%parent)THEN
          IF (soft_com%exsoft) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,16("OACP"))')
             IF (paral%io_parent)&
                  WRITE(6,'(" OACP",5X,A,11X,"OACP")')&
                  'OPTIMIZATION of ATOM CENTERED POTENTIALS'
             IF (paral%io_parent)&
                  WRITE(6,*) 'OACP   WARNING! WARNING!WARNING!WARNING!WARNING!'
             IF (paral%io_parent)&
                  WRITE(6,*) 'OACP    - is NOT converged'
             IF (paral%io_parent)&
                  WRITE(6,*) 'OACP   WARNING! WARNING!WARNING!WARNING!WARNING!'
             IF (paral%io_parent)&
                  WRITE(6,'(1X,16("OACP"),/)')
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,16("OACP"))')
             IF (paral%io_parent)&
                  WRITE(6,'(" OACP",5X,A,11X,"OACP")')&
                  'OPTIMIZATION of ATOM CENTERED POTENTIALS'
             IF (paral%io_parent)&
                  WRITE(6,*) 'OACP    - is converged'
             IF (paral%io_parent)&
                  WRITE(6,'(1X,16("OACP"),/)')
          ENDIF
       ENDIF
    ENDIF                     ! that TDUMMYATOM

    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(rho_ref,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho_tmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nuclear_p
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE oacp4forces(c0,c1,sc0,pme,gde,vpp,eigv,&
       rhoe,psi,nstate)
    ! ==================================================================
    COMPLEX(real_8)                          :: sc0(*), pme(nkpt%ngwk,*), &
                                                gde(nkpt%ngwk,*)
    REAL(real_8)                             :: vpp(:), eigv(*), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(maxfftn,clsd%nlsd)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(len=2)                         :: numb
    CHARACTER(len=20)                        :: ecpname
    INTEGER                                  :: isp, k, sig
    REAL(real_8) :: crittie, deriv(2,response2%tdumnum), initd, &
      parm_(2,response2%tdumnum), parm_old(2,response2%tdumnum), &
      tp(2,response2%tdumnum), tPold(2,response2%tdumnum)

! variables for updwf
! variables
! ==--------------------------------------------------------------==
! Initializiation

    initd = 0.01_real_8
    k  = 1
    CALL rwfopt_nuc(c0,c1,sc0,pme,gde,vpp,eigv,nstate,psi,rhoe)
    crittie = 0.0_real_8
    DO isp = 2,response2%tdumnum
       parm_(1,isp)=sgpp2%rcnl(response2%ang_mom,isp)
       parm_(2,isp)=sgpp2%hlsg(1,1,response2%ang_mom,isp)
       DO sig = 1,2
          CALL forcepenalty(tp(sig,isp))
          IF (tp(sig,isp).GT.crittie) crittie = tp(sig,isp)
       ENDDO
    ENDDO

    DO WHILE (response2%dumcrit.LT.crittie)
       crittie = 0.0_real_8
       DO isp = 1,response2%tdumnum
          DO sig = 1,2

             ! get finite difference derivative
             IF (k.EQ.1) THEN
                deriv(sig,isp) = initd
             ELSE
                deriv(sig,isp) =&
                     (tPold(sig,isp)-tp(sig,isp))/(parm_old(sig,isp)-parm_(sig,isp))
             ENDIF
             tPold(sig,isp) = tp(sig,isp)
             parm_old(sig,isp) = parm_(sig,isp)

             ! change parameters
             IF (paral%parent)THEN
                IF (sig.EQ.1)THEN
                   CALL steepestdescent(sgpp2%rcnl(response2%ang_mom,isp),&
                        response2%dumstep4,deriv(sig,isp))
                   parm_(sig,isp) = sgpp2%rcnl(response2%ang_mom,isp)
                ELSEIF (sig.EQ.2) THEN
                   CALL steepestdescent(sgpp2%hlsg(1,1,response2%ang_mom,isp),&
                        response2%dumstep5,deriv(sig,isp))
                   parm_(sig,isp) = sgpp2%hlsg(1,1,response2%ang_mom,isp)
                ENDIF
                IF (paral%io_parent)&
                     WRITE (numb,'(I2.2)') isp
                ecpname = 'DUMMY_'//numb
                CALL write_pp(isp,ecpname)
                IF (paral%io_parent)&
                     WRITE(6,*) ' OACP psedopotential written to: ',ecpname
                CALL recpnew(isp,ecpname)
                IF (paral%io_parent)&
                     WRITE(6,'(/,1X,16("OACP"))')
                IF (paral%io_parent)&
                     WRITE(6,'(" OACP",16X,A,18X,"OACP")') '    new parameters    '
                IF (paral%io_parent)&
                     WRITE(6,*)'param',isp,sig,parm_(sig,isp)
                IF (paral%io_parent)&
                     WRITE(6,'(1X,16("OACP")/)')
             ENDIF
             CALL rinforce_nuc(response2%tdumnum)
             CALL rwfopt_nuc(c0,c1,sc0,pme,gde,vpp,eigv,nstate,psi,rhoe)
             CALL forcepenalty(tp(sig,isp))
             IF (tp(sig,isp).GT.crittie) crittie = tp(sig,isp)
             IF (paral%parent)THEN
                IF (paral%io_parent)&
                     WRITE(6,'(/,1X,16("OACP"))')
                IF (paral%io_parent)&
                     WRITE(6,'(" OACP",16X,A,18X,"OACP")') '  FORCE OPTIMIZATION  '
                IF ((k.EQ.1).AND.paral%io_parent)&
                     WRITE(6,*) k,'st OACP-results for ACP',isp
                IF ((k.EQ.2).AND.paral%io_parent)&
                     WRITE(6,*) k,'nd OACP-results for ACP',isp
                IF ((k.EQ.3).AND.paral%io_parent)&
                     WRITE(6,*) k,'rd OACP-results for ACP',isp
                IF ((k.GT.3).AND.paral%io_parent)&
                     WRITE(6,*) k,'th OACP-results for ACP',isp
                IF (paral%io_parent)&
                     WRITE(6,*) '                   parameter                  pen'
                IF (paral%io_parent)&
                     WRITE(6,*)'penalty',sig,parm_(sig,isp),tp(sig,isp)
                IF (paral%io_parent)&
                     WRITE(6,'(1X,16("OACP")/)')
             ENDIF

          ENDDO
       ENDDO
       k = k + 1
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE oacp4forces
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE forcepenalty(tPtmp)
    ! ==================================================================
    REAL(real_8)                             :: tPtmp

    INTEGER                                  :: i, iat, isp

! variables
! ==--------------------------------------------------------------==

    tPtmp = 0.0_real_8
    DO isp = 1,maxsp
       DO iat = 1,ions0%na(isp)
          IF (response2%tfdimen.EQ.0)THEN! choose dimensionality
             DO i = 1,3
                tPtmp = tptmp + ABS(fion(i,iat,isp))
             ENDDO
          ELSEIF (response2%tfdimen.EQ.1)THEN
             tPtmp = tptmp + ABS(fion(1,iat,isp))
          ELSEIF (response2%tfdimen.EQ.2)THEN
             tPtmp = tptmp + ABS(fion(2,iat,isp))
          ELSEIF (response2%tfdimen.EQ.3)THEN
             tPtmp = tptmp + ABS(fion(3,iat,isp))
          ENDIF
       ENDDO
    ENDDO
    CALL mp_sum(tPtmp,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE forcepenalty
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE densitypenalty(rho_ref,rho_tmp,tp,refvol)
    ! ==================================================================
    REAL(real_8)                             :: rho_ref(*), rho_tmp(*), tp
    LOGICAL                                  :: refvol

    INTEGER                                  :: i1, i2, i3, iat, iii, spec
    LOGICAL                                  :: vol
    REAL(real_8)                             :: dist, omevol

! variables
! ==--------------------------------------------------------------==

    tp = 0.0_real_8
    omevol=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    IF (refvol) THEN
       DO i3 = 1,fpar%kr3s
          DO i2 = 1,fpar%kr2s
             iii = 0
             DO i1 = parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
                iii = iii + 1
                vol = .FALSE.
                DO spec = 1,response2%tdumnum! loop1 to check if u r outside RAD_DUM
                   CALL distance(1,spec,i1,i2,i3,dist)
                   IF (dist.LT.response2%rad_dum) vol = .TRUE.
                ENDDO
                DO spec = response2%tdumnum+1,maxsp! loop2 to check if u r outside RAD_NORM
                   DO iat = 1,ions0%na(spec)
                      CALL distance(iat,spec,i1,i2,i3,dist)
                      IF (dist.LT.response2%rad_norm) vol = .TRUE.
                   ENDDO
                ENDDO
                IF (vol) tp = tp&
                     + ((rho_ref(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s)&
                     - rho_tmp(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s))**2) *omevol
             ENDDO
          ENDDO
       ENDDO
       CALL mp_sum(tp,parai%allgrp)
    ELSE
       DO i1=1,fpar%nnr1
          tp = tp + ((rho_ref(i1)-rho_tmp(i1))**2)*omevol
       ENDDO
       CALL mp_sum(tp,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE densitypenalty
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE dummy_optcalc(c0,c1,psi,nstate,rhoe,drhoe,&
       eirop,eirop1,eivps,z11,rho_ref,rho_tmp,vtemp,crit_tmp,&
       pme,gde,sc0,vpp,eigv)
    ! ==================================================================
    COMPLEX(real_8)                          :: psi(maxfftn,clsd%nlsd)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: rhoe(fpar%nnr1), &
                                                drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eirop1(ncpw%nhg), eivps(*)
    REAL(real_8)                             :: z11(nstate,nstate), &
                                                rho_ref(*), rho_tmp(*)
    COMPLEX(real_8)                          :: vtemp(*)
    REAL(real_8)                             :: crit_tmp
    COMPLEX(real_8)                          :: pme(nkpt%ngwk,*), &
                                                gde(nkpt%ngwk,*), sc0(*)
    REAL(real_8)                             :: vpp(:), eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dummy_optcalc'

    CHARACTER(len=2)                         :: numb
    CHARACTER(len=20)                        :: ecpname
    COMPLEX(real_8), ALLOCATABLE             :: pertXphi(:,:,:,:), temp_p(:,:)
    INTEGER                                  :: i1, i2, i3, iat, ierr, iii, &
                                                isp, numofs, sig, spec
    LOGICAL                                  :: vol
    REAL(real_8)                             :: dist, &
                                                dP(5,response2%tdumnum), &
                                                omevol
    REAL(real_8), ALLOCATABLE                :: rho1(:,:,:)

! reference and temporary rho
! temporary convergence criterium value
! variables for rwfopt_nuc
! variables
! species, atom,species
! total number of parameters, parameter
! 5 is the number of sigma
! ==--------------------------------------------------------------==

    ALLOCATE(temp_p(ncpw%ngw,ncpw%ngw*nstate/ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rho1(fpar%nnr1,5,fpar%nnr1*response2%tdumnum/fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pertXphi(ncpw%ngw,nstate,5,ncpw%ngw*nstate*response2%tdumnum/(ncpw%ngw*&
         nstate)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    iat = 1

    ! get perturbation Hamiltonian into pertXphi for x, state is, parameter sig, species isp 
    numofs = 5 ! = number of sig's of dummy isp
    DO isp = 1,response2%tdumnum
       CALL pert4nuc(isp,iat,c0,rhoe,nstate,numofs,pertXphi)
    ENDDO

    ! call linear response to get perturbed density rho1 from pertXphi for sig and isp
    DO isp = 1,response2%tdumnum
       DO sig = 1,numofs
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,16("OACP"))')
             IF (paral%io_parent)&
                  WRITE(6,'(" OACP",16X,A,17X,"OACP")') ' RESPONSE CALCULATION  '
             IF (paral%io_parent)&
                  WRITE(6,'(" OACP",14X,A,I4,A,I4,14X,"OACP")')&
                  'for parameter',sig,' of ACP',isp
             IF (paral%io_parent)&
                  WRITE(6,'(1X,16("OACP")/)')
          ENDIF
          CALL dcopy(2*ncpw%ngw*nstate,pertXphi(1,1,sig,isp),1,temp_p(1,1),1)
          CALL rwfopt_p(c0(1,1),c1(1,1),psi,rhoe,drhoe,eirop,&
               eivps,vtemp(1),temp_p(1,1),z11,nstate,eirop1)
          CALL rhoofr_p(c0,c1,drhoe,psi,nstate)

          DO i3 = 1,fpar%kr3s
             DO i2 = 1,fpar%kr2s
                iii = 0
                DO i1 = parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
                   iii = iii + 1
                   rho1(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s,sig,isp) = 0.0_real_8
                   vol = .FALSE.
                   DO spec = 1,response2%tdumnum! loop1 to check if u r outside RAD_DUM
                      CALL distance(1,spec,i1,i2,i3,dist)
                      IF (dist.LT.response2%rad_dum) vol = .TRUE.
                   ENDDO
                   DO spec = response2%tdumnum+1,maxsp! loop2 to check if u r outside RAD_NORM
                      DO iat = 1,ions0%na(spec)
                         CALL distance(iat,spec,i1,i2,i3,dist)
                         IF (dist.LT.response2%rad_norm) vol = .TRUE.
                      ENDDO
                   ENDDO
                   IF (vol) rho1(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s,sig,isp)&
                        = - 2 * drhoe(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s) ! missing factor from Daniel's routine
                ENDDO
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    ! compute the penalty derivative wrt 5 parameters sig for all species isp,
    ! **** dP = 2*dSIG* \int dr^3 rho1 * (rho_tmp-rho_ref)  
    omevol=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    DO isp = 1,response2%tdumnum
       DO sig = 1,numofs
          dP(sig,isp) = 0.0_real_8
          DO iii = 1,fpar%nnr1
             dP(sig,isp) = dp(sig,isp) + rho1(iii,sig,isp) * 2&
                  * (rho_tmp(iii)-rho_ref(iii)) * omevol
          ENDDO
          CALL mp_sum(dP(sig,isp),parai%allgrp)
       ENDDO
    ENDDO

    ! update all sig for all species and PP-file
    DO isp = 1, response2%tdumnum
       IF (paral%parent)THEN
          CALL steepestdescent(sgpp2%rcsg(isp),response2%dumstep1,dP(1,isp))
          CALL steepestdescent(sgpp2%clsg(1,isp),response2%dumstep2,dP(2,isp))
          CALL steepestdescent(sgpp2%clsg(2,isp),response2%dumstep3,dP(3,isp))
          CALL steepestdescent(sgpp2%rcnl(1,isp),response2%dumstep4,dP(4,isp))
          CALL steepestdescent(sgpp2%hlsg(1,1,1,isp),response2%dumstep5,dP(5,isp))
          IF (paral%io_parent)&
               WRITE (numb,'(I2.2)') isp
          ecpname = 'DUMMY_'//numb
          CALL write_pp(isp,ecpname)
          IF (paral%io_parent)&
               WRITE (6,*)' OACP psedopotential written to: ',ecpname
          CALL recpnew(isp,ecpname)
       ENDIF
    ENDDO

    ! update the wavefunction
    response1%response_running=.FALSE.
    CALL rinforce_nuc(response2%tdumnum)
    CALL rwfopt_nuc(c0,c1,sc0,pme,gde,vpp,eigv,nstate,psi,rhoe)
    response1%response_running=.TRUE.

    ! print out penalty derivatives
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",19X,A,20X,"*")') ' RESPONSE RESULTS      '
       IF (paral%io_parent)&
            WRITE(6,*) '*   ACP  PARAM  PENALTY GRADIENT:'
    ENDIF
    crit_tmp = 0.0_real_8
    DO isp = 1,response2%tdumnum
       DO sig = 1,numofs
          IF (paral%io_parent)&
               WRITE(6,'(" *",I5,2X,I5,F16.8)')isp,sig,dP(sig,isp)
          IF (crit_tmp.LT.ABS(dP(sig,isp))) crit_tmp=ABS(dp(sig,isp))
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"),/)')

    ! ==--------------------------------------------------------------==
    DEALLOCATE(temp_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pertXphi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dummy_optcalc
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE steepestdescent(sig,delta,dP)
    REAL(real_8)                             :: sig, delta, dP

! ==--------------------------------------------------------------==

    sig = sig - delta * ABS(sig) * dP
    RETURN
  END SUBROUTINE steepestdescent
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE pert4nuc(isp,iat,c0,rhoe,nstate,numofs,pertXphi)
    ! ==================================================================
    ! computes the real space perturbational Hamiltonian for atom (isp,iat),
    ! i.e. the derivative wrt all parameters, of the 1996 SG PP in
    ! real space for l = s (hence, there are only 5 sigma parameters).
    ! The perturbing potential needs to be applied to every orbital and
    ! is fftransformed and copied into pertXphi
    ! J. Chem. Phys. 122, 14113 (2005).
    ! ==================================================================
    INTEGER                                  :: isp, iat
    REAL(real_8)                             :: rhoe(fpar%nnr1)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: numofs
    COMPLEX(real_8) :: pertXphi(ncpw%ngw,nstate,numofs,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pert4nuc'
    REAL(real_8), PARAMETER                  :: cutof = 1.e-10_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: aux(:), temp(:,:)
    INTEGER                                  :: i1, i2, i3, ierr, iii, is, sig
    REAL(real_8)                             :: dist, omevol, term(4)
    REAL(real_8), ALLOCATABLE                :: deltaV(:,:,:,:), &
                                                rhoxyz(:,:,:), totxyz(:,:,:,:)

! ==--------------------------------------------------------------==

    ALLOCATE(deltaV(fpar%kr1,fpar%kr2s,fpar%kr3s,numofs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoxyz(fpar%kr1,fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(totxyz(fpar%kr1,fpar%kr2s,fpar%kr3s,numofs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(temp(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! Compute local parts of h'
    DO i3 = 1,fpar%kr3s
       DO i2 = 1,fpar%kr2s
          iii = 0
          DO i1 = parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
             iii  = iii + 1
             CALL distance(iat,isp,i1,i2,i3,dist)
             deltaV(iii,i2,i3,1) = EXP(-0.5_real_8*(dist**2)/(sgpp2%rcsg(isp)**2))&
                  * ( (sgpp2%clsg(1,isp)-2*sgpp2%clsg(2,isp))*(dist**2)/(sgpp2%rcsg(isp)**3)&
                  +                sgpp2%clsg(2,isp)*(dist**4)/(sgpp2%rcsg(isp)**5)&
                  +       ions0%zv(isp)*(2**0.5_real_8)/((pi**(0.5_real_8))*(sgpp2%rcsg(isp)**2)))

             deltaV(iii,i2,i3,2) = EXP(-0.5_real_8*(dist**2)/(sgpp2%rcsg(isp)**2))
             deltaV(iii,i2,i3,3) = deltav(iii,i2,i3,2)*(dist/sgpp2%rcsg(isp))**2
             IF (ABS(deltaV(iii,i2,i3,1)).LT.cutof) THEN
                deltaV(iii,i2,i3,1) = 0.0_real_8
             ENDIF
             IF (ABS(deltaV(iii,i2,i3,2)).LT.cutof) THEN
                deltaV(iii,i2,i3,2) = 0.0_real_8
             ENDIF
             IF (ABS(deltaV(iii,i2,i3,3)).LT.cutof) THEN
                deltaV(iii,i2,i3,3) = 0.0_real_8
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! Compute nonlocal parts of h' and apply to orbital rhoe: h'|phi> 
    ALLOCATE(aux(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    omevol=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    DO is = 1,nstate
       CALL ffttor(c0(1,is),rhoe,aux,ncpw%ngw,.TRUE.)

       ! Compute the first separable part
       CALL zeroing(term)!,4)
       DO i3 = 1,fpar%kr3s
          DO i2 = 1,fpar%kr2s
             iii = 0
             DO i1 = parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
                iii = iii + 1
                CALL distance(iat,isp,i1,i2,i3,dist)
                term(1) = rhoe(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s) *&
                     EXP(-(dist**2)/(2*sgpp2%rcnl(1,isp)**2))*omevol
                term(2) = term(2) + term(1)
                term(3) = term(3) + term(1) * dist**2
                IF (ABS(term(1)).LT.cutof) term(1) = 0.0_real_8
                IF (ABS(term(2)).LT.cutof) term(2) = 0.0_real_8
                IF (ABS(term(3)).LT.cutof) term(3) = 0.0_real_8
             ENDDO
          ENDDO
       ENDDO
       CALL mp_sum(term(2),parai%allgrp)
       CALL mp_sum(term(3),parai%allgrp)

       ! Compute the other part and apply also the local part to rhoe
       CALL zeroing(rhoxyz)!,kr1*kr2s*kr3s)
       DO i3 = 1,fpar%kr3s
          DO i2 = 1,fpar%kr2s
             iii = 0
             DO i1 = parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
                iii = iii + 1
                CALL distance(iat,isp,i1,i2,i3,dist)
                term(4) = EXP(-(dist**2)/(2*sgpp2%rcnl(1,isp)**2))
                IF (ABS(term(4)).LT.cutof) term(4) = 0.0_real_8
                deltaV(iii,i2,i3,4) = term(4) * ( term(2)&
                     * ((dist**2)*sgpp2%hlsg(1,1,1,isp)/((pi**1.5_real_8)*sgpp2%rcnl(1,isp)**6)-&
                     3*sgpp2%hlsg(1,1,1,isp)/((pi**1.5_real_8)*sgpp2%rcnl(1,isp)**4))&
                     + term(3)*sgpp2%hlsg(1,1,1,isp)/((pi**1.5_real_8)*sgpp2%rcnl(1,isp)**6))
                deltaV(iii,i2,i3,5) = term(2)*term(4)/&
                     ((pi**1.5_real_8)*(sgpp2%rcnl(1,isp)**3))
                IF (ABS(deltaV(iii,i2,i3,4)).LT.cutof)&
                     deltaV(iii,i2,i3,4) = 0.0_real_8
                IF (ABS(deltaV(iii,i2,i3,5)).LT.cutof)&
                     deltaV(iii,i2,i3,5) = 0.0_real_8
                DO sig = 1,numofs
                   totxyz(iii,i2,i3,sig) = deltaV(iii,i2,i3,sig)*&
                        rhoe(iii+(i2-1)*fpar%kr1+(i3-1)*fpar%kr1*fpar%kr2s)
                   IF (sig.GE.4) totxyz(iii,i2,i3,sig) = deltaV(iii,i2,i3,sig)! the nonlocal parts are already multiplied with rhoe
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       ! fft of h'|phi> for every state (is), sigma (xx), and species (isp)
       DO sig = 1,numofs
          CALL ffttog(totxyz(1,1,1,sig),temp(1,is),aux,ncpw%ngw,.TRUE.)
          CALL dcopy(2*ncpw%ngw,temp(1,is),1,pertXphi(1,is,sig,isp),1)
       ENDDO

    ENDDO ! over nstate
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    DEALLOCATE(deltaV,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoxyz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(totxyz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pert4nuc
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE dummy_refcalc(c0,psi,nstate,rhoe)
    ! ==================================================================
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: rhoe(fpar%nnr1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dummy_refcalc'

    CHARACTER(len=20)                        :: filename
    COMPLEX(real_8), POINTER                 :: c0_3d_p(:,:,:)
    INTEGER                                  :: ierr
    REAL(real_8), ALLOCATABLE                :: aux(:,:)

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,*)'DUMMY_REFCALC'
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==

    IF (response1%tdummyatom_ref) THEN
       ! get the electrondensity from c0 in rho0
       CALL rhoofr(c0,rho0,psi(:,1),nstate)

       ! write out rhoofr_ref.dat for reference density
#if defined( _HASNT_F08_POINTER_REMAPPING )
       CALL stopgm(procedureN,'compiler needs to support pointer remapping!',&
            __LINE__,__FILE__)
#else
       c0_3d_p(1:SIZE(c0,1),1:SIZE(c0,2),1:1) => c0
#endif
       CALL rhopri(c0_3d_p,tau0,rho0,psi(:,1),nstate,1)
       filename = 'rhoofr_ref.cube'
       ALLOCATE(aux(fpar%kr2s,fpar%kr3s),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL cubefile(filename,rho0,tau0,aux,.FALSE.)
       DEALLOCATE(aux,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE dummy_refcalc
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE distance(iat,isp,i1,i2,i3,dist) ! get distance to an atom
    ! ==================================================================
    INTEGER                                  :: iat, isp, i1, i2, i3
    REAL(real_8)                             :: dist

! ==--------------------------------------------------------------==

    dist =  SQRT((tau0(1,iat,isp)-(i1-1)*parm%a1(1)/spar%nr1s)**2&
         +(tau0(2,iat,isp)-(i2-1)*parm%a2(2)/spar%nr2s)**2&
         +(tau0(3,iat,isp)-(i3-1)*parm%a3(3)/spar%nr3s)**2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE distance
  ! ==================================================================

END MODULE nuclear_p_utils
