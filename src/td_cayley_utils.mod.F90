MODULE td_cayley_utils
  USE cnst,                            ONLY: uimag
  USE coor,                            ONLY: fion,&
                                             tau0
  USE dipo_utils,                      ONLY: dipo
  USE dipomod,                         ONLY: moment
  USE eicalc_utils,                    ONLY: eicalc
  USE ener,                            ONLY: chrg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn
  USE geq0mod,                         ONLY: geq0
  USE hpsi_utils,                      ONLY: hpsi
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw,&
                                             nkpt
  USE td_input,                        ONLY: &
       gaugefield, gndir, gpot, gpotv, itermax, masklog, niterations, tcurr, &
       td_prop
  USE td_utils,                        ONLY: &
       applymask, currentoper, dext_ch, dext_pipot, dext_pot, dipolox, &
       dyn_analysis, eh_initialize_data, gaugepot_laser, &
       gaugepot_laser_circ_pol, kick, pos_oper, tmpwr_prop
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean_k
  USE vofrho_utils,                    ONLY: vofrho
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_cayley

CONTAINS

  ! ==================================================================
  ! Time Dependent Propagation using cntl%cayley interpolation
  ! author: I. Tavernelli
  ! 
  ! ==================================================================
  SUBROUTINE td_cayley(c0,c2,rhoe,psi,sc0,effpot,nstate,eigv,first,do_dipole)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:), sc0(*)
    REAL(real_8)                             :: effpot(fpar%nnr1,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    REAL(real_8)                             :: eigv(nstate)
    LOGICAL                                  :: first, do_dipole

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_cayley'

    CHARACTER(len=30)                        :: f_dipole, f_el, f_gauge_i, &
                                                f_gauge_l, f_time, tag
    COMPLEX(real_8), ALLOCATABLE             :: cor(:,:), dummy(:), eirop(:), &
                                                hphi(:,:), phi(:,:), &
                                                respo(:), rhog(:), v(:,:)
    INTEGER                                  :: ciclo, ierr, ir, irec(100), &
                                                is, ispin, isub, k, lst
    INTEGER, SAVE                            :: icall = 0
    LOGICAL                                  :: debug, second
    LOGICAL, SAVE                            :: write_dipole_diff = .FALSE.
    REAL(real_8) :: dipox, dmom0(3), elfield, fgpot_i0, fgpot_i1, &
      gpot_i0 = 0.0_real_8, gpot_i1 = 0.0_real_8, norm_l, rhoemax, &
      Time_array(10), work_time
    REAL(real_8), ALLOCATABLE                :: hs(:), neigv(:), norms(:), &
                                                pos(:), rhoe0(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: rhoer(:)
    REAL(real_8), SAVE                       :: gpot_l0 = 0.0_real_8, &
                                                gpot_l1 = 0.0_real_8, &
                                                last_time = 0.0_real_8

    CALL tiset(procedureN,isub)

    ! 
    ! Meaning of the time variables:
    ! work_time if icall.eq.0 then 
    ! gets its value from the file
    ! else
    ! gets its value from last_time
    ! endif
    ! In case of cntl%ruku it takes the intermidiate value
    ! work_time+tt2 after half of the step (first=.true.)
    ! last_time if icall.eq.0 then
    ! gets its value from work_time
    ! else
    ! gets its value from the one save after the previous execution
    ! endif
    ! In case of cntl%ruku it only get updated if second=.true.
    ! it never takes the intermidiat value as work_time
    ! 
    ! ttime     is a global time variable linked to the file abstime in ehrenfest.F
    ! it gets update in td_cayley
    ! At present only used with td_extpot
    ! 
    ! tintrvll  is a time intervall (not an absolute time)
    ! if (cntl%ruku and first) tintrvll = tt2
    ! if (cntl%ruku and second) tintrvll = tt
    ! 
    second=.NOT.first
    f_dipole="dipole.dat"
    f_gauge_i="gaugefield_ind.dat"
    f_gauge_l="gaugefield_laser.dat"
    f_el="electricfield.dat"
    f_time="Time.dat"
    ! 
    ispin=clsd%nlsd
    ! 
    ALLOCATE(phi(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hphi(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hs(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(v(maxfft,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoe0(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pos(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cor(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(norms(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(neigv(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (td_prop%tpointch) THEN
       ALLOCATE(respo(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rhog(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (icall.EQ.0) ALLOCATE(rhoer(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! avoid gfortran runtime error 'not allocated'
       ALLOCATE(respo(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rhog(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (icall.EQ.0) ALLOCATE(rhoer(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! 
    debug=.FALSE.
    ! ==--------------------------------------------------------------==
    ! ==-  open perturbation field file                              -==
    ! ==--------------------------------------------------------------==
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,1X,F15.10)') 'Time intervall    :',td_prop%tintrvll
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,1X,F15.10)') 'accuracy          :',td_prop%epsil
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,1X,F15.10)') 'pert. amplitude   :',td_prop%ampl
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,1X,I15)')    '# cycles          :',td_prop%n_cycles
       IF (cntl%tpspec) THEN
          IF (td_prop%pert_type.EQ.1) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,A)')&
                  'EXTERNAL PULSE (DIRAC) AT t=0'
          ELSEIF (td_prop%pert_type.EQ.2)  THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,A)')&
                  'EXTERNAL CONSTANT (HEAVISIDE) FIELD STARTING FROM t=0'
          ENDIF
          IF (td_prop%pertdir.EQ.1) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'Direction perturbing field: x'
          ELSEIF (td_prop%pertdir.EQ.2) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'Direction perturbing field: y'
          ELSEIF (td_prop%pertdir.EQ.3) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'Direction perturbing field: z'
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'Isotropic perturbation'
          ENDIF
       ENDIF
    ENDIF
    ! 
    ! ==--------------------------------------------------------------==
    ! ==- Initialize TD run                                          -==
    ! ==--------------------------------------------------------------==
    ! 
    ! Initialize input orbitals PHI
    ! -----------------------------
    ! 
    DO lst=1,nstate
       CALL dcopy(2*nkpt%ngwk,c0(1,lst),1,phi(1,lst),1)
       IF (cntl%tpspec) CALL dcopy(2*nkpt%ngwk,c0(1,lst),1,hphi(1,lst),1)
    ENDDO
    ! 
    CALL pos_oper(pos,td_prop%pertdir)
    CALL rhoofr_c(phi,rhoe,psi,nstate)
    norm_l=chrg%csumg
    IF ((paral%parent.AND.debug) .AND.paral%io_parent)&
         WRITE(6,'(1X,A,1X,F10.4)') 'norm of the density ',norm_l
    ! 
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoe0,1)
    ! 
    ! ==--------------------------------------------------------------==
    ! ==- Initialize variables and files                             -==
    ! ==--------------------------------------------------------------==
    CALL eh_initialize_data(rhoe,rhoer,rhog,v(:,:),psi,eirop,elfield,&
         nstate,norms,do_dipole,dmom0,work_time,last_time,&
         gpot_i0,gpot_i1,gpot_l0,gpot_l1,ispin,icall,&
         f_dipole,f_time,f_gauge_i,f_gauge_l,f_el)
    ! 
    IF (paral%parent) THEN
       IF ((do_dipole).AND.paral%io_parent)&
            OPEN(unit=50,file=f_dipole,status='unknown',&
            position='append')
       IF ((td_prop%tpointch).AND.paral%io_parent)&
            OPEN(unit=53,file=f_time,status='unknown',&
            position='append')
       IF ((cntl%tgaugep).AND.paral%io_parent)&
            OPEN(unit=72,file=f_gauge_i,status='unknown',&
            position='append')
       IF ((cntl%tgaugef).AND.paral%io_parent)&
            OPEN(unit=73,file=f_gauge_l,status='unknown',&
            position='append')
       IF ((td_prop%td_extpot).AND.paral%io_parent)&
            OPEN(unit=74,file=f_el,status='unknown',&
            position='append')
    ENDIF
    ! 
    ! ==--------------------------------------------------------------==
    ! ==- Perturbation of the KS orbitals                            -==
    ! ==--------------------------------------------------------------==
    ! 
    IF ((cntl%tpspec.AND.td_prop%stextpot.AND.(icall.EQ.0)).OR.&
         (cntl%tmdeh.AND.(td_prop%stextpot.AND.(td_prop%pert_type.EQ.1)).AND.(icall.EQ.0))&
                                ! &    (cntl%tmdeh.and.do_dipole.and.icall.eq.0)
         ) THEN
       ! The reference value of the dipole at time 0 is needed
       ! as reference for the calculation of the delta dipoles.
       ! We need to store this value to be used in case of restart!
       write_dipole_diff=.TRUE.
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(50,'(4D25.15)') 0.0_real_8,moment%dmom(1),moment%dmom(2),moment%dmom(3)
       ENDIF
       IF (td_prop%pert_type.EQ.1) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Applying perturbation ... '
          CALL zeroing(phi)!,SIZE(phi))
          CALL kick(hphi,nstate,norms,pos,phi,psi)
          CALL rhoofr_c(phi,rhoe,psi,nstate)
          norm_l=chrg%csumg
          CALL zeroing(v)!,SIZE(v))
          IF ((paral%parent.AND.debug).AND.paral%io_parent)&
               WRITE(6,*)  'Norm of the pert. density ',norm_l
          !$omp parallel do private(ir)
          DO ir=1,fpar%nnr1
             v(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL dipo(tau0,eirop,v)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1x,a,3f12.6)') 'Dipole after perturbation',&
                  moment%dmom(1),moment%dmom(2),moment%dmom(3)
          ENDIF
       ELSEIF (td_prop%pert_type.EQ.2) THEN
       ENDIF
    ENDIF
    ! 
    icall=icall+1
    ! 
    ! ==--------------------------------------------------------------==
    ! ==- Initialize cycle                                           -==
    ! ==--------------------------------------------------------------==
    ! 
    IF (cntl%ruku) THEN
       IF (first) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1x,A)') 'Computes effective potential...'
          CALL rhoofr_c(phi,rhoe,psi,nstate)
          norm_l=chrg%csumg
          IF (td_prop%td_extpot) THEN
             IF (gaugefield%pi_pulse) THEN
                CALL dext_pipot(td_prop%ttime,elfield)
             ELSE
                CALL dext_pot(td_prop%ttime,elfield,0)
             ENDIF
          ENDIF
          IF (td_prop%tpointch) CALL dext_ch(work_time)
          CALL vofrho(tau0,fion,rhoe,v,.FALSE.,.FALSE.)
          DO is=1,ispin
             CALL dcopy(fpar%nnr1,rhoe(1,is),1,effpot(1,is),1)
          ENDDO
          IF (cntl%tgaugep) THEN
             fgpot_i0=gpot_i0
             fgpot_i1=gpot_i1
          ENDIF
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)') 'Uses previous effective potential...'
       ENDIF
    ENDIF
    ! 
    work_time=last_time
    ! ----------------------------------------------------------------==
    DO ciclo=1,td_prop%n_cycles
       ! ----------------------------------------------------------------==
       ! 
       IF (.NOT.cntl%ruku) THEN
          CALL rhoofr_c(phi,rhoe,psi,nstate)
          norm_l=chrg%csumg
          CALL vofrho(tau0,fion,rhoe,v,.FALSE.,.FALSE.)
       ENDIF
       ! 
       IF (cntl%tpspec.AND.(td_prop%pert_type.EQ.2))&
            CALL daxpy(fpar%nnr1,td_prop%ampl,pos,1,rhoe,1)
       ! 
       IF (cntl%tgaugep) THEN
          DO lst=1,nstate
             CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,c0(1,lst),1)
          ENDDO
       ENDIF
       ! 
       IF (paral%io_parent)&
            WRITE(6,'(/,A,/)')&
            ' PROPAGATION USING CAYLEY ALGORITHM'
       IF (cntl%ruku) THEN
          CALL cayley_step(phi,hphi,psi,nstate,effpot,sc0,neigv,first)
       ELSE
          CALL cayley_step(phi,hphi,psi,nstate,rhoe,sc0,neigv,first)
       ENDIF
       ! 
       ! ==--------------------------------------------------------------==
       ! ==-  update the absolute time variables                        -==
       ! ==--------------------------------------------------------------==
       ! 
       IF (cntl%ruku) THEN
          work_time=work_time+td_prop%tintrvll
          ! if (first) do not update last_time
          IF (second) last_time=work_time
       ELSE
          work_time=work_time+td_prop%tintrvll
          last_time=work_time
       ENDIF
       ! 
       ! ==--------------------------------------------------------------==
       ! ==-  calculation of the dipol                                  -==
       ! ==--------------------------------------------------------------==
       ! 
       IF (do_dipole.AND.second) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)') ' Calculation of the dipole'
          ! 
          ! -- calculate the density 
          CALL zeroing(rhoe)!,nnr1*clsd%nlsd)
          CALL rhoofr_c(phi,rhoe,psi,nstate)
          norm_l=chrg%csumg
          CALL dipolox(pos,rhoe(:,1),dipox)
          IF (paral%io_parent)&
               WRITE(49,'(2D25.15)') work_time,dipox
          CALL zeroing(v)!,SIZE(v))
          CALL zeroing(eirop)!,SIZE(eirop))
          ALLOCATE(dummy(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL eicalc(dummy,eirop)
          DEALLOCATE(dummy,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          rhoemax=0._real_8
          !$omp parallel do private(ir)
          DO ir=1,fpar%nnr1
             v(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
             IF (ABS(rhoe(ir,1)-rhoe0(ir,1)).GT.rhoemax)&
                  rhoemax=ABS(rhoe(ir,1)-rhoe0(ir,1))
          ENDDO
          CALL mp_max(rhoemax,parai%allgrp)
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL dipo(tau0,eirop,v)
          IF (paral%io_parent) THEN
             IF (write_dipole_diff) THEN
                WRITE(50,'(4D25.15)') work_time,&
                     (moment%dmom(1)-dmom0(1)),&
                     (moment%dmom(2)-dmom0(2)),&
                     (moment%dmom(3)-dmom0(3))
             ELSE
                WRITE(50,'(4D25.15)') work_time,&
                     moment%dmom(1),moment%dmom(2),moment%dmom(3)
             ENDIF
          ENDIF
          ! 
          IF ((paral%parent.AND.debug).AND.paral%io_parent)&
               WRITE(6,'(4D14.6)')&
               work_time,moment%dmom(1),moment%dmom(2),moment%dmom(3)
       ENDIF
       ! 
       ! ==--------------------------------------------------------------==
       ! ==-  update external potentials
       ! ==--------------------------------------------------------------==
       ! 
       IF (.NOT.cntl%ruku) THEN

          td_prop%ttime=td_prop%ttime+td_prop%tintrvll
          IF (td_prop%td_extpot) THEN
             IF (gaugefield%pi_pulse) THEN
                CALL dext_pipot(td_prop%ttime,elfield)
             ELSE
                CALL dext_pot(td_prop%ttime,elfield,1)
             ENDIF
          ENDIF
          IF (td_prop%tpointch) CALL dext_ch(work_time)
          IF (td_prop%tpointch) THEN
             CALL dyn_analysis(first,phi,rhoe,rhoer,rhog,psi,&
                  nstate,respo,work_time)
          ENDIF
          IF (tcurr)&
               CALL currentoper(phi,rhoe,nstate,tau0,psi,1)
          IF (cntl%tgaugep) THEN
             DO lst=1,nstate
                CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,c2(1,lst),1)
             ENDDO
             Time_array(1)=td_prop%tintrvll
             ! call gaugepot_ind(gpot_i0,gpot_i1,1,time,c0,c2)
          ENDIF
          IF (cntl%tgaugef.AND.(gndir.EQ.0)) THEN
             CALL gaugepot_laser(gpot_l0,gpot_l1,work_time)
          ENDIF
          IF (cntl%tgaugef.AND.(gndir.GT.0)) THEN
             CALL gaugepot_laser_circ_pol(work_time)
          ENDIF
          IF ((paral%parent.AND.cntl%tgaugep).AND.paral%io_parent)&
               WRITE(72,'(3f20.12)') work_time,gpot_i0,gpot_i1
          IF (paral%io_parent.AND.cntl%tgaugef.AND.(gndir.EQ.0))&
               WRITE(73,'(3f20.12)') work_time,gpot_l0,gpot_l1
          IF (paral%io_parent.AND.cntl%tgaugef.AND.(gndir.GT.0))&
               WRITE(73,'(4f20.12)') work_time,(gpotv(k),k=1,3)
          IF ((paral%parent.AND.td_prop%td_extpot) .AND.paral%io_parent)&
               WRITE(74,'(2f20.12)') td_prop%ttime,elfield
          IF ((paral%parent.AND.td_prop%tpointch) .AND.paral%io_parent)&
               WRITE(53,'(3f20.12)') work_time

       ELSE

          IF (cntl%tgaugep) THEN
             DO lst=1,nstate
                CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,c2(1,lst),1)
             ENDDO
             Time_array(1)=td_prop%tintrvll
             IF (second) THEN
                gpot_i0=fgpot_i0
                gpot_i1=fgpot_i1
             ENDIF
             ! call gaugepot_ind(gpot_i0,gpot_i1,1,time,c0,c2)
          ENDIF
          IF (cntl%tgaugef.AND.(gndir.EQ.0)) THEN
             CALL gaugepot_laser(gpot_l0,gpot_l1,work_time)
          ENDIF
          IF (cntl%tgaugef.AND.(gndir.GT.0)) THEN
             CALL gaugepot_laser_circ_pol(work_time)
          ENDIF

          IF (second) THEN
             IF ((paral%parent.AND.cntl%tgaugep) .AND.paral%io_parent)&
                  WRITE(72,'(3f20.12)') work_time,gpot_i0,gpot_i1
             IF (paral%io_parent.AND.cntl%tgaugef.AND.(gndir.EQ.0))&
                  WRITE(73,'(3f20.12)') work_time,gpot_l0,gpot_l1
             IF (paral%io_parent.AND.cntl%tgaugef.AND.(gndir.GT.0))&
                  WRITE(73,'(4f20.12)') work_time,(gpotv(k),k=1,3)
             IF ((paral%parent.AND.td_prop%tpointch).AND.paral%io_parent)&
                  WRITE(53,'(3f20.12)') work_time
          ENDIF

       ENDIF

       ! if (cntl%tgaugep) gpot=gpot+gpot_i0
       IF (cntl%tgaugep) gpot=gpot_i0
       IF (cntl%tgaugef) THEN
          IF (gndir.EQ.0) THEN
             gpot=gpot+gpot_l0
          ENDIF
       ENDIF
       ! 
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) GOTO 111

       ! ----------------------------------------------------------------==
    ENDDO   ! n_cycles
    ! ----------------------------------------------------------------==
111 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==-  'THAT s ALL FOLKS                                         -==
    ! ==--------------------------------------------------------------==
    ! 
    DO lst=1,nstate
       CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,c0(1,lst),1)
    ENDDO

    IF (masklog%tmask) THEN
       IF ((.NOT.cntl%ruku).OR.(cntl%ruku.AND.second)) THEN
          CALL applymask(c0,psi,nstate,icall,1)
       ENDIF
    ENDIF
    ! 
    IF (cntl%ruku.AND.first) THEN
       CALL rhoofr_c(c0,rhoe,psi,nstate)
       norm_l=chrg%csumg
       td_prop%ttime=td_prop%ttime+td_prop%tt2
       IF (td_prop%td_extpot) THEN
          IF (gaugefield%pi_pulse) THEN
             CALL dext_pipot(td_prop%ttime,elfield)
          ELSE
             CALL dext_pot(td_prop%ttime,elfield,1)
          ENDIF
       ENDIF
       IF (td_prop%tpointch) CALL dext_ch(work_time)
       td_prop%ttime=td_prop%ttime+td_prop%tt2
       CALL vofrho(tau0,fion,rhoe,v,.FALSE.,.FALSE.)
       DO is=1,ispin
          CALL dcopy(fpar%nnr1,rhoe(1,is),1,effpot(1,is),1)
       ENDDO
    ELSEIF (cntl%ruku.AND.second) THEN
       IF (td_prop%tpointch) THEN
          CALL dyn_analysis(first,c0,rhoe,rhoer,rhog,psi,&
               nstate,respo,work_time)
       ENDIF
       IF (tcurr)&
            CALL currentoper(c0,rhoe,nstate,tau0,psi,1)
       IF (paral%parent.AND.(td_prop%td_extpot)) THEN
          IF (paral%io_parent)&
               WRITE(74,'(2f20.12)') td_prop%ttime,elfield
       ENDIF
    ENDIF
    ! 
    ! CALL TCORREL(C0,PHI,NSTATE,COR)
    ! WRITE(6,*)real(COR(1,1)),aimag(COR(1,2))
    ! 
    ! ==--------------------------------------------------------------==
    ! ==-  write phi to RESTART.CWF.1 file                           -==
    ! ==--------------------------------------------------------------== 
    IF ((.NOT.cntl%ruku).OR.(cntl%ruku.AND.second)) THEN
       irec(9)=1
       ! CALL ZHWWF(2,IREC,C0, CM, NSTATE,EIGV,TAU0,VELP,TAUI,NFI)
       tag="wavefunctions"
       CALL tmpwr_prop(c0,nstate,1,tag)
    ENDIF
    CALL mp_sync(parai%allgrp)

    IF (paral%parent) THEN
       IF ((do_dipole).AND.paral%io_parent)&
            CLOSE(50)
       IF ((td_prop%tpointch).AND.paral%io_parent)&
            CLOSE(53)
       IF ((cntl%tgaugep).AND.paral%io_parent)&
            CLOSE(72)
       IF ((cntl%tgaugef).AND.paral%io_parent)&
            CLOSE(73)
       IF ((td_prop%td_extpot).AND.paral%io_parent)&
            CLOSE(74)
    ENDIF
    ! ==--------------------------------------------------------------== 
    IF (td_prop%tpointch) THEN
       DEALLOCATE(rhog,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(respo,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(neigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(norms,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cor,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hphi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(phi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE td_cayley
  ! ==================================================================
  SUBROUTINE cayley_step(phi,hphi,psi,nstate,rhoe,sc0,norms,first)
    ! ==--------------------------------------------------------------== 
    ! 
    COMPLEX(real_8)                          :: phi(:,:), psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: hphi(nkpt%ngwk,nstate)
    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: sc0(*)
    REAL(real_8)                             :: norms(*)
    LOGICAL                                  :: first

    CHARACTER(*), PARAMETER                  :: procedureN = 'cayley_step'

    COMPLEX(real_8)                          :: phiold
    COMPLEX(real_8), ALLOCATABLE             :: hphi0(:,:), phi0(:,:)
    COMPLEX(real_8), EXTERNAL                :: zdotc
    INTEGER                                  :: ialpha, ibeta, ierr, ig, &
                                                iter, lst
    REAL(real_8)                             :: delta, deltamax, deltaold, &
                                                tinc2

! 
! Arguments:
! Variables:
! (ngwk,*)
! (ngwk,*)
! 

    ALLOCATE(phi0(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hphi0(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    itermax=1
    ! 
    tinc2=td_prop%tintrvll/2.0_real_8
    !$omp parallel do private(lst,ig)
    DO lst=1,nstate
       DO ig=1,nkpt%ngwk
          phi0(ig,lst)=phi(ig,lst)
       ENDDO
    ENDDO
    ! 
    CALL hpsi(phi0(:,1:nstate),hphi0,sc0,rhoe,psi,nstate,1,clsd%nlsd)
    CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,hphi0(1,1),1)
    ! 
    iter=0
    deltamax=1._real_8
    niterations=0
    ! ---------------------------------------------------------------
    IF (cntl%ruku.AND.first) THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A50)')&
            '#iter           dt/4       max.grad      precision'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(2X,A50)')&
            '#iter           dt/2       max.grad      precision'
    ENDIF
700 CONTINUE
    niterations=niterations+1
    deltaold=deltamax
    CALL hpsi(phi(:,1:nstate),hphi,sc0,rhoe,psi,nstate,1,clsd%nlsd)
    CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,hphi(1,1),1)

    deltamax=0._real_8
    DO lst=1,nstate
       DO ig=1,nkpt%ngwk
          phiold=phi(ig,lst)

          phi(ig,lst)=phi0(ig,lst)&
               - uimag * tinc2 * hphi0(ig,lst)/cntr%emass&
               - uimag * tinc2 * hphi(ig,lst)/cntr%emass

          delta=ABS(phiold-phi(ig,lst))
          IF (delta.GT.deltamax) deltamax=delta
       ENDDO
    ENDDO
    CALL mp_max(deltamax,parai%allgrp)
    iter=iter+1
    IF (paral%io_parent)&
         WRITE(6,'(2X,I5,3(E15.6))') iter,tinc2,deltamax,td_prop%epsil
    ! if ((deltamax.gt.epsil).and.(deltamax.le.deltaold)) goto 700
    IF (deltamax.GT.td_prop%epsil) GOTO 700
    ! ---------------------------------------------------------------

    IF (geq0) CALL zclean_k(phi,nstate,ncpw%ngw)
    DO lst=1,nstate
       norms(lst)=REAL(zdotc(nkpt%ngwk,phi(1,lst),1,phi(1,lst),1))
       CALL mp_sum(norms(lst),parai%allgrp)
    ENDDO
    ! 5 next lines to check the norms during the propagation (optional)
    ! IF (paral%io_parent) THEN
    !  do lst=1,nstate
    !    write(*,*) 'norms in td_cayley',norms(lst)
    !  enddo
    ! endif
    !
    IF (td_prop%do_not_normalize) THEN
       ! do nothing
    ELSE
       ! renormalize
       DO lst=1,nstate
          CALL dscal(2*nkpt%ngwk,1._real_8/SQRT(norms(lst)),phi(1,lst),1)
       ENDDO
    ENDIF

    IF (iter.GT.itermax) itermax=iter
    ! 
    IF (paral%io_parent) THEN
       WRITE(6,*) '  -- CAYLEY MAX ITERATION STEPS ',itermax
    ENDIF

    IF (td_prop%td_fix_spin_dens) THEN
       ! contrain the spin density assuming a doublet
       ialpha=1
       DO ibeta=1,spin_mod%nsdown
          IF (ialpha.EQ.td_prop%ionized_state) THEN
             ! do nothing
          ELSE
             DO ig=1,nkpt%ngwk
                phi(ig,ialpha)=0.5_real_8 * (phi(ig,ialpha)+phi(ig,ibeta+spin_mod%nsup))
                phi(ig,ibeta+spin_mod%nsup)=phi(ig,ialpha)
             ENDDO
             ialpha=ialpha+1
          ENDIF
       ENDDO

    ENDIF


    DEALLOCATE(hphi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(phi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cayley_step
  ! ==================================================================

END MODULE td_cayley_utils
