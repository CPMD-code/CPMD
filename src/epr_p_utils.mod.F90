MODULE epr_p_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             inyh,&
                                             nzh
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE epr_current_p_utils,             ONLY: epr_current,&
                                             give_scr_epr_current
  USE epr_dv0_utils,                   ONLY: dv0
  USE epr_hyp_utils,                   ONLY: epr_hyp
  USE epr_util_p_utils,                ONLY: epr_do_full
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: localize
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE nmr_full_p_utils,                ONLY: make_listofstates,&
                                             optimize_llc_indep
  USE nmr_position_p_utils,            ONLY: calc_lower_left_new,&
                                             calc_lower_left_optimize,&
                                             force_xyz,&
                                             print_llc,&
                                             set_origin
  USE nmr_util_p_utils,                ONLY: make_123,&
                                             print_wavefunctions,&
                                             printtime
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phfac_utils,                     ONLY: phfac
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE qspl,                            ONLY: ggnh,&
                                             nsplpo
  USE ragg,                            ONLY: raggio
  USE response_pmod,                   ONLY: &
       epr_options, firstshelldistance, iepr_csgt, iepr_default, iepr_iglo, &
       iepr_llc, iepr_novirt, iepr_wc, lower_left, lower_left_value, &
       nmr_options, nris, ownpotvalue, response1, response2, response_read, &
       response_write, timetag, vofrho0, wanniercenters
  USE restart_p_utils,                 ONLY: restart_epr,&
                                             restart_p
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE soft,                            ONLY: soft_com
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun,&
                                             simpsn
  USE wann,                            ONLY: wannl,&
                                             wannr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_p
  PUBLIC :: give_scr_epr
  !public :: trmm
  !public :: epr_effpot

CONTAINS

  ! ==================================================================
  SUBROUTINE epr_p(c0,h1psi0,psi,rhoe,&
       eirop,eivps,z11,nstate)
    ! ==================================================================
    ! PARALLEL
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(fpar%nnr1,2)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: h1psi0(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'epr_p'

    CHARACTER(len=1)                         :: xyz(3) = (/'x','y','z'/)
    CHARACTER(len=10)                        :: tag
    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: c1(:,:,:), c1b(:,:,:), &
                                                c1t(:,:,:), cs(:), sc0(:), &
                                                vt(:), vt2(:)
    INTEGER :: dummy2, i1, i2, i3, i_state, ib, ier, ierr, ii, iii, ik, &
      iL_end, il_start, ip_end, ip_start, is, isub, isub2, l_c1, nc
    INTEGER, ALLOCATABLE                     :: neighborlist(:,:)
    LOGICAL                                  :: full_optimization, &
                                                simple_done(9)
    LOGICAL, ALLOCATABLE                     :: doepr(:), dowfopt(:,:)
    REAL(real_8) :: dg_rmc, dg_sootemp, dg_sotemp, dummy, epr_chi(3,3), &
      geigval(3), rdummy(1), starttime, w_eps_save
    REAL(real_8), ALLOCATABLE :: aux(:,:), BIND(:,:), dg_so(:,:), &
      dg_soo(:,:), diffv0a(:,:), diffv0b(:,:), geigvec(:,:), gtensor(:,:), &
      j_a(:,:), j_b(:,:), rhoetmp(:,:), rootrho(:)

! ==--------------------------------------------------------------==
! Functions

89  FORMAT ("*  FULL calculation: field in O",a1,32x,"*")
96  FORMAT ("*  PERTURBATION: ",a1,"_",i1," |PSI>",38x,"*")
99  FORMAT ("*  Calculation of responses DONE.",31x,"*")
97  FORMAT ("*  TOTAL TIME: ",f7.1,42x,"*")
98  FORMAT (52("*"),"EPR*RESPONSE*")
    ! ==--------------------------------------------------------------==
    CALL tiset('       EPR',isub)

    ALLOCATE(aux(2*maxfft,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    timetag=' '
    IF (paral%parent) CALL printtime
    starttime=m_walltime()/1000._real_8

    nris(1)   = spar%nr1s
    nris(2)   = spar%nr2s
    nris(3)   = spar%nr3s
    response1%pr_energy = epr_options%teverbose     ! these values are not sooo interesting
    full_optimization=.FALSE.
    ! full_optimization=.TRUE.
    ! IF (parent) print *,
    ! &     'DEBUG DEBUG: >> full_optimization: ',full_optimization

    ! ==--------------------------------------------------------------==
    ! Some initialization 
    ! ==--------------------------------------------------------------==
    l_c1 = 3*2*ncpw%ngw*nstate     ! 6: 2 x pert.param.for.3B-direct.
    IF (epr_options%tepr_full) l_c1 = 3*2*ncpw%ngw*nstate
    ip_start=1
    ip_end  =3
    il_start=1
    iL_end  =3

    ALLOCATE(c1(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1)!,SIZE(c1))
    ALLOCATE(c1b(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1b)!,SIZE(c1b))
    ALLOCATE(c1t(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1t)!,SIZE(c1t))
    ALLOCATE(j_a(fpar%nnr1,9),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(j_a)!,nnr1*9)
    ALLOCATE(j_b(fpar%nnr1,9),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(j_b)!,nnr1*9)
    ALLOCATE(BIND(fpar%nnr1,9),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(BIND)!,nnr1*9)
    ALLOCATE(dg_so(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dg_so)!,3*3)
    ALLOCATE(dg_soo(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dg_soo)!,3*3)
    ALLOCATE(gtensor(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gtensor)!,3*3)
    ALLOCATE(geigvec(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(geigvec)!,3*3)
    ALLOCATE(lower_left(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(lower_left)!,SIZE(lower_left))
    ALLOCATE(lower_left_value(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(lower_left_value)!,SIZE(lower_left_value))
    ALLOCATE(diffv0a(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(diffv0a)!,nnr1*3)
    ALLOCATE(diffv0b(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(diffv0b)!,nnr1*3)
    ALLOCATE(rhoetmp(fpar%nnr1,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rhoetmp)!,2*nnr1)
    ALLOCATE(vt(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vt)!,maxfft)
    ALLOCATE(vt2(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vt2)!,nhg)
    ALLOCATE(rootrho(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rootrho)!,nnr1)
    ALLOCATE(doepr(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    DO ib=1,3
       simple_done(ib)=.FALSE.
       simple_done(ib+3)=.FALSE.
       simple_done(ib+6)=.FALSE.
    ENDDO
    IF (paral%parent) THEN
       ALLOCATE(wanniercenters(4,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(wanniercenters)!,4*nstate)
    ENDIF

    IF (paral%parent) CALL prmem('EPR        ')
    firstshelldistance = 1._real_8 ! in a.u.

    IF (paral%io_parent) WRITE(6,*)

    ! ==--------------------------------------------------------------==
    IF (epr_options%iepr_method.EQ.iepr_default) THEN
       epr_options%iepr_method = iepr_csgt
       IF (isos1%tclust) THEN
          epr_options%iepr_virtual = iepr_novirt
       ENDIF
    ENDIF
    IF (epr_options%iepr_virtual.EQ.iepr_default) THEN
       epr_options%iepr_virtual = iepr_llc
       IF (epr_options%iepr_method.EQ.iepr_iglo) epr_options%iepr_virtual = iepr_wc
    ENDIF
    IF (epr_options%iepr_virtual.EQ.iepr_novirt) epr_options%telocalize = .FALSE.
    IF (response1%t_restart_c1 .AND. cntl%wfopt) THEN
       CALL stopgm('EPR/RESTART',&
            'Use of NOOPT is compulsory for EPR RESTART',& 
            __LINE__,__FILE__)
    ENDIF
    ! Synchronize NMR/EPR options
    nmr_options%inmr_virtual = epr_options%iepr_virtual
    nmr_options%inmr_method = epr_options%iepr_method

    ! ==--------------------------------------------------------------==
    ! Set aside NablaV 
    ! 2 options:
    ! TEPROWNPOT = .FALSE. : use local potential from cntl%wfopt
    ! Works fine with GO potentials
    ! Works not so fine with TM potentials
    ! TEPROWNPOT = .TRUE. : creates effective potential
    ! V = -Z_val / r * erf( r / EPR_OWNPOT_VALUE )
    ! + V_Hartree + V_xc
    ! 
    ! V differs for the 2 spin channels
    ! 
    ! NOT THOROUGHLY TESTED
    ! ==--------------------------------------------------------------==

    IF (epr_options%tepr_ownpot) THEN
       IF (paral%io_parent) THEN
          WRITE(6,'(A,T54,E12.3)')&
               '* EPR|OWNPOT Creating eff. pot. with smooth. value : ',&
               ownpotvalue%epr_ownpot_value
          WRITE(6,*)
       ENDIF
       CALL epr_effpot(rhoetmp(1,1),rhoe,vt,vt2)
       ! Calculate gradient
       CALL dv0(rhoetmp(1,1),vt,diffv0a(1,1),vt2)
       CALL dv0(rhoetmp(1,2),vt,diffv0b(1,1),vt2)
    ELSE
       ! Calculate gradient
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') '* EPR|Using local potential from wavefunction optimization'
       ENDIF
       CALL dv0(vofrho0(1,1),vt,diffv0a(1,1),vt2)
       CALL dv0(vofrho0(1,2),vt,diffv0b(1,1),vt2)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate root of norm of spin density
    ! ==--------------------------------------------------------------==

    DO i1=1,fpar%nnr1
       rhoetmp(i1,1) = rhoe(i1,1) - 2._real_8*rhoe(i1,2)
       rootrho(i1) = SQRT(ABS(rhoetmp(i1,1)))
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Calculate Hyperfine Coupling Tensors
    ! ==--------------------------------------------------------------==

    IF (epr_options%tepr_hyp) CALL epr_hyp(c0,rhoetmp,psi(:,1),nstate)

    ! ==--------------------------------------------------------------==
    ! Calculate deltaG-RMC contribution 
    ! ==--------------------------------------------------------------==
    dg_rmc = 0._real_8
    DO i1=1,spin_mod%nsup
       dg_rmc = dg_rmc + kine_one(c0(1,i1))
    ENDDO
    DO i1=spin_mod%nsup+1,nstate
       dg_rmc = dg_rmc - kine_one(c0(1,i1))
    ENDDO
    ! zeeman + hyperfine^2 + units correction
    dg_rmc = 2.0023192778_real_8 * ( 1._real_8 / 137.03602_real_8 )**2 *&
         dg_rmc

    ! ==--------------------------------------------------------------==
    ! Localize C0 
    ! ==--------------------------------------------------------------==

    CALL phfac(tau0)
    CALL eicalc(eivps,eirop)  ! calculate and store EIVPS and EIROP
    IF (epr_options%telocalize) THEN
       wannl%twann=.TRUE.
       prop1%locl=.TRUE.
       prop2%numorb=MAX(prop2%numorb,nstate)
       lenbk=nxxfun(crge%n)
       nc=MAX(2*lenbk*parai%nproc,2*nkpt%ngwk*prop2%numorb)
       ! Create memory associated with ip_sc0, ip_cs
       ALLOCATE(sc0(nc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cs(nc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Save/Set the convergence criteria
       w_eps_save = wannr%w_eps
       wannr%w_eps      = 1.e-4_real_8
       ! localize
       CALL localize(tau0,c0,cs,sc0,nstate)
       wannr%w_eps      = w_eps_save
       store1%swf        = .TRUE.
       ! Free memory associated with ip_sc0, ip_cs
       DEALLOCATE(sc0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    timetag='localization.'
    IF (paral%parent) CALL printtime

    ! ==--------------------------------------------------------------==
    ! If FULL and SMART, leave out some calculations
    ! ==--------------------------------------------------------------==

    DO i1=1,nstate
       doepr(i1)=.TRUE.
    ENDDO

    IF (epr_options%tepr_smart .AND. epr_options%tepr_full) THEN
       dummy2=0
       DO i1=1,nstate
          CALL ffttor(c0(1,i1),rhoe(1,1),psi,ncpw%ngw,.FALSE.)
          dummy=0._real_8
          DO i2=1,fpar%nnr1
             dummy=ABS(rhoe(i2,1))*rootrho(i2)+dummy
          ENDDO
          CALL mp_sum(dummy,parai%allgrp)
          IF (dummy*SQRT(parm%omega)/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8) .LT.&
               ownpotvalue%epr_smart_value) THEN
             doepr(i1)=.FALSE.
             dummy2=dummy2+1
          ENDIF
       ENDDO
       dummy=0._real_8
       IF (paral%io_parent) THEN
          WRITE(6,'(A,T54,E12.3)')&
               '* EPR|SMART Threshold value : ',&
               ownpotvalue%epr_smart_value
          WRITE(6,'(A,T62,I4)')&
               '* EPR|SMART Orbitals left out of full calculation : ',&
               dummy2
          WRITE(6,*)
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate C1 
    ! ==--------------------------------------------------------------==

    CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)
    timetag='lagrange multipliers.'
    IF (paral%parent) CALL printtime
    ! ==--------------------------------------------------------------==
    IF (nmr_options%tprintwfns)&
         CALL print_wavefunctions(c0,nstate,psi(:,1))
    IF (nmr_options%tprintrho)&
         CALL print_orbital_densities(c0,nstate,psi,aux)
    ! ==--------------------------------------------------------------==
    IF (epr_options%iepr_virtual.EQ.iepr_llc) THEN
       IF (epr_options%tepr_overlaps) THEN
          CALL calc_lower_left_optimize(c0,c1,nstate)
       ELSE
          CALL calc_lower_left_new(c0,nstate)
          IF (full_optimization) THEN
             ! TESTING: LLC optimization...
             ALLOCATE(neighborlist(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(dowfopt(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)

             firstshelldistance= MIN(1.5_real_8, parm%alat/20._real_8)
             IF (paral%io_parent) THEN
                WRITE (6,*) 'PROCEDURE:'
                WRITE (6,*) ' new LLC opti / 1 x ',&
                     firstshelldistance,' a.u.'
             ENDIF

             CALL optimize_llc_indep(neighborlist,dowfopt,nstate)

             firstshelldistance=0.1_real_8
             CALL make_listofstates(neighborlist,dowfopt,nstate)

             DEALLOCATE(neighborlist,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             DEALLOCATE(dowfopt,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ENDIF             ! full_optimization
       ENDIF
    ELSEIF (epr_options%iepr_virtual.EQ.iepr_novirt) THEN
       DO is=1,nstate        ! -> All wfns have same virtual box
          DO ik=1,3
             lower_left_value(ik,is) = -nris(ik)/2
             lower_left(ik,is) = 1
          ENDDO
       ENDDO
    ENDIF                     ! virtual orbitals
    CALL force_xyz(nstate)    ! if the corresponding options FORCE(x,...)
    ! are switched on,
    ! the X/Y/Z - LLCs are averaged.
    IF (paral%parent) CALL print_llc(nstate)

    timetag='LLC initialization.'
    IF (paral%parent) CALL printtime
    ! ==--------------------------------------------------------------==
    IF (response1%t_restart_c1) CALL restart_epr(simple_done,&
         nstate,response_read)

    ! ==--------------------------------------------------------------==
    DO ib=il_start,iL_end
       IF (paral%io_parent) WRITE (6,98)
       IF (paral%io_parent) WRITE (6,96) 'L',ib

       tag = 'L_'//xyz(ib)
       IF (response1%t_restart_c1 .AND. simple_done(ib))&
            CALL restart_p(c1b(1,1,ib),tag,nstate,response_read)
       IF (.NOT. simple_done(ib)) THEN
          CALL make_123(ib,ii,iii)
          DO is=1,nstate    ! apply r_ii p_iii - r_iii p_ii:
             CALL set_origin (is)
             CALL apply_op_p(c0(1,is), aux(1,1), ii, ncpw%ngw)
             CALL apply_op_p(c0(1,is), aux(1,2), iii, ncpw%ngw)
             CALL fft2tor(aux(1,1), aux(1,1),&
                  aux(1,2),aux(1,2), psi, ncpw%ngw,.FALSE.)
             CALL apply_op_rx(aux(1,1), -1._real_8, rhoe,  0._real_8, iii)
             CALL apply_op_rx(aux(1,2),  1._real_8, rhoe,  1._real_8,  ii)
             CALL ffttog(rhoe, h1psi0(1,is), psi, ncpw%ngw,.FALSE.)
          ENDDO

          timetag='preparation.'
          IF (paral%parent) CALL printtime
          IF (paral%io_parent) WRITE (6,98)
          CALL tiset('PSI1-WF-OPT',isub2)
          CALL rwfopt_p(c0,c1b(1,1,ib),psi,rhoe,rdummy,&
               eirop,eivps,cdummy,h1psi0,&
               z11,nstate,cdummy)
          CALL tihalt('PSI1-WF-OPT',isub2)
          timetag='calculation of psi^1.'
          IF (paral%parent) CALL printtime
          IF (soft_com%exsoft) GOTO 441
          IF (paral%io_parent) WRITE (6,98)

          tag = 'L_'//xyz(ib)
          CALL restart_p(c1b(1,1,ib),tag,nstate,response_write)
          simple_done(ib) = .TRUE.
          CALL restart_epr(simple_done,&
               nstate,response_write)
       ENDIF                 ! if simple calculation not already done
    ENDDO                     ! iB
    ! ==--------------------------------------------------------------==
    DO ib=ip_start,ip_end
       IF (paral%io_parent) WRITE (6,98)
       IF (paral%io_parent) WRITE (6,96) 'p',ib

       tag = 'p_'//xyz(ib)
       IF (response1%t_restart_c1 .AND. simple_done(ib+3))&
            CALL restart_p(c1(1,1,ib),tag,nstate,response_read)
       IF (.NOT. simple_done(ib+3)) THEN
          DO is=1,nstate
             CALL apply_op_p(c0(1,is), h1psi0(1,is), ib, ncpw%ngw)
          ENDDO
          timetag='preparation.'
          IF (paral%parent) CALL printtime
          IF (paral%io_parent) WRITE (6,98)
          CALL tiset('PSI1-WF-OPT',isub2)
          CALL rwfopt_p(c0,c1(1,1,ib),psi,rhoe,rdummy,&
               eirop,eivps,cdummy,h1psi0,&
               z11,nstate,cdummy)
          CALL tihalt('PSI1-WF-OPT',isub2)
          IF (paral%io_parent) WRITE (6,98)
          timetag='calculation of psi^1.'
          IF (paral%parent) CALL printtime
          IF (soft_com%exsoft) GOTO 441
          tag = 'p_'//xyz(ib)
          CALL restart_p(c1(1,1,ib),tag,nstate,response_write)
          simple_done(ib+3) = .TRUE.
          CALL restart_epr(simple_done,&
               nstate,response_write)
       ENDIF                 ! if simple calculation not already done
    ENDDO                     ! iB
441 CONTINUE
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE (6,99)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (soft_com%exsoft) GOTO 444
    IF (epr_options%tepr_full) THEN
       IF (paral%io_parent) WRITE(6,*) '*** STARTING FULL Calculation'
       ALLOCATE(neighborlist(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dowfopt(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL make_listofstates(neighborlist,dowfopt,nstate)
       ! neighborlist now contains the states in a reasonable 
       ! proximity-determined order, and dowfopt indicates whether or not
       ! rwfopt_p should be called for that state.

       IF (.NOT. full_optimization) THEN
          DO is=1,nstate
             DO ib=1,3
                dowfopt(ib,is) = .TRUE.
             ENDDO
          ENDDO
          IF (paral%io_parent) THEN
             WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
             WRITE(6,*)'*** *** UNCONDITIONAL  RWFOPT_P:    *** *** ***'
             WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
          ENDIF
       ELSE
          IF (paral%parent) THEN
             WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
             WRITE(6,*)'*** *** ONLY SELECTIVE RWFOPT_P:    *** *** ***'
             WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
          ENDIF
       ENDIF
       DO is=1,nstate
          DO ib=1,3
             dowfopt(ib,is)=dowfopt(ib,is) .AND. doepr(is)
          ENDDO
       ENDDO

       response2%tolog_p = response2%tolog_p * 1.e2_real_8
       ! 100 is a large number, but it seems to work...
       response1%cg_analytic = 0

       DO ib=1,3
          IF (.NOT. simple_done(ib+6)) THEN
             DO is=1,nstate      ! "is" is only an INDEX!!!
                ! The state is i_state !
                i_state = neighborlist(ib,is)
                IF (doepr(i_state)) THEN
                   CALL epr_do_full(i_state,ib,&
                        dowfopt(ib,i_state),&
                        h1psi0,c0,c1,c1b,c1t,psi(:,1),rhoe,&
                        eirop,eivps,z11,nstate)
                ELSE
                   IF (paral%io_parent) WRITE(6,'(A,T29,I4,32X,"*")') &
                        '*  EPR|SM ARTOmitting state ',i_state
                ENDIF
                IF (soft_com%exsoft) GOTO 442
             ENDDO         ! i_state = 1,...,nstate.
             tag = 'L_'//xyz(ib)
             CALL restart_p(c1b(1,1,ib),tag,nstate,response_write)
             simple_done(ib+6) = .TRUE.
             CALL restart_epr(simple_done,&
                  nstate,response_write)
          ELSE
             IF (paral%io_parent) THEN
                WRITE (6,98)
                WRITE (6,89) xyz(ib)
                WRITE (6,'(A)') 'RESPONSE wavefunction read from RESTART file'
                WRITE (6,98)
             ENDIF
          ENDIF
       ENDDO                 ! iB = 1,2,3

       ! 442  CALL freem(ip_neighborlist)
442    DEALLOCATE(neighborlist,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! CALL freem(ip_dowfopt)
       DEALLOCATE(dowfopt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent) WRITE(6,*)

    ENDIF                     ! full calculation
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! Calc induced Current and Magnetic Field
    ! ==--------------------------------------------------------------==

    CALL epr_current(c0,c1,c1b,c1t,nstate,BIND(1,1),j_a(1,1),j_b(1,1),&
         psi(:,1),epr_chi(1,1))

    ! ==--------------------------------------------------------------==
    ! Printing output from earlier calculation
    ! ==--------------------------------------------------------------==

    CALL mp_sum(dg_rmc,parai%allgrp)

    IF (paral%io_parent) THEN
       WRITE(6,'(/,A,A)') ' ********************** DEBUG INFORMATION',&
            ' ***********************'

       WRITE(6,*) ' EPR deltaG-MATRIX RMC CORRECTION'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_rmc,0._real_8,0._real_8,' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            0._real_8,dg_rmc,0._real_8,' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            0._real_8,0._real_8,dg_rmc,' ]'
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate deltaG-SOO contribution          
    ! ==--------------------------------------------------------------==

    ! G eq 0. Assuming a spheric sample -> kappa = 2/3 (hard-coded)

    DO i1=1,3
       DO i2=1,3
          dg_sootemp = 0._real_8
          DO i3=1,fpar%nnr1
             dg_sootemp = dg_sootemp - 2._real_8/3._real_8 * epr_chi(i1,i2)&
                  *rhoetmp(i3,1)
          ENDDO
          ! Perform the following:
          ! * OMEGA/real(nr1s*nr2s*nr3s,kind=real_8) ! from int(r x j(r))
          ! * FPI * alpha**2 / OMEGA
          dg_sootemp = dg_sootemp / REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)* fpi *&
               (1._real_8/137.03602_real_8)**2 / 2._real_8
          ! And now the SOO term, for G=0
          dg_soo(i1,i2) = dg_soo(i1,i2) +&
               2._real_8*dg_sootemp*(1._real_8/137.03602_real_8)/&
               REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       ENDDO                 ! i2
    ENDDO                     ! i1

    ! G neq 0

    DO i1=1,3
       DO i2=1,3
          dg_sootemp = 0._real_8
          DO i3=1,fpar%nnr1
             dg_sootemp = dg_sootemp -&
                  BIND(i3, (i1-1)*3 + i2) * rhoetmp(i3,1)
          ENDDO
          dg_soo(i1,i2) = dg_soo(i1,i2) +&
               2._real_8*dg_sootemp*(1._real_8/137.03602_real_8)/&
               REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       ENDDO                 ! i2
    ENDDO                     ! i1

    CALL mp_sum(dg_soo,3*3,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,*) ' EPR deltaG-MATRIX SOO CORRECTION'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_soo(1,1),dg_soo(1,2),dg_soo(1,3),' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_soo(2,1),dg_soo(2,2),dg_soo(2,3),' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_soo(3,1),dg_soo(3,2),dg_soo(3,3),' ]'
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate deltaG-SpinOrbit contribution 
    ! ==--------------------------------------------------------------==

    DO i1=1,3
       dg_sotemp = 0._real_8
       DO i3=1,fpar%nnr1
          dg_sotemp = dg_sotemp-(j_a(i3,(i1-1)*3+2)*&
               diffv0a(i3,3)-j_a(i3,(i1-1)*3+3)*diffv0a(i3,2) -&
               j_b(i3,(i1-1)*3+2)*diffv0b(i3,3) +&
               j_b(i3,(i1-1)*3+3)*diffv0b(i3,2))
       ENDDO
       dg_so(i1,1) = (1._real_8/137.03602_real_8)**2*2.0046385556_real_8/2._real_8*&
            dg_sotemp/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)

       dg_sotemp = 0._real_8
       DO i3=1,fpar%nnr1
          dg_sotemp = dg_sotemp-(j_a(i3,(i1-1)*3+3)*&
               diffv0a(i3,1)-j_a(i3,(i1-1)*3+1)*diffv0a(i3,3) -&
               j_b(i3,(i1-1)*3+3)*diffv0b(i3,1) +&
               j_b(i3,(i1-1)*3+1)*diffv0b(i3,3))
       ENDDO
       dg_so(i1,2) = (1._real_8/137.03602_real_8)**2*2.0046385556_real_8/2._real_8*&
            dg_sotemp/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)

       dg_sotemp = 0._real_8
       DO i3=1,fpar%nnr1
          dg_sotemp = dg_sotemp-(j_a(i3,(i1-1)*3+1)*&
               diffv0a(i3,2)-j_a(i3,(i1-1)*3+2)*diffv0a(i3,1) -&
               j_b(i3,(i1-1)*3+1)*diffv0b(i3,2) +&
               j_b(i3,(i1-1)*3+2)*diffv0b(i3,1))
       ENDDO
       dg_so(i1,3) = (1._real_8/137.03602_real_8)**2*2.0046385556_real_8/2._real_8*&
            dg_sotemp/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ENDDO                     ! i1

    CALL mp_sum(dg_so,3*3,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,*) ' EPR deltaG-MATRIX SO CORRECTION'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_so(1,1),dg_so(1,2),dg_so(1,3),' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_so(2,1),dg_so(2,2),dg_so(2,3),' ]'
       WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
            dg_so(3,1),dg_so(3,2),dg_so(3,3),' ]'
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Print G-tensor 
    ! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,64("*"))')
          WRITE(6,*)
          WRITE(6,'(/,A,A)') ' ************************ FINAL RESULTS',&
               ' *************************'
       ENDIF

       gtensor(1,1) = dg_rmc+dg_so(1,1)+dg_soo(1,1)+2.0023192778_real_8
       gtensor(1,2) = dg_so(1,2)+dg_soo(1,2)
       gtensor(1,3) = dg_so(1,3)+dg_soo(1,3)
       gtensor(2,1) = dg_so(2,1)+dg_soo(2,1)
       gtensor(2,2) = dg_rmc+dg_so(2,2)+dg_soo(2,2)+2.0023192778_real_8
       gtensor(2,3) = dg_so(2,3)+dg_soo(2,3)
       gtensor(3,1) = dg_so(3,1)+dg_soo(3,1)
       gtensor(3,2) = dg_so(3,2)+dg_soo(3,2)
       gtensor(3,3) = dg_so(3,3)+dg_soo(3,3)+dg_rmc+2.0023192778_real_8

       IF (paral%io_parent) THEN
          WRITE(6,*) ' EPR deltaG-MATRIX'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(1,1)-2.0023192778_real_8,gtensor(1,2),gtensor(1,3),' ]'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(2,1),gtensor(2,2)-2.0023192778_real_8,gtensor(2,3),' ]'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(3,1),gtensor(3,2),gtensor(3,3)-2.0023192778_real_8,' ]'

          WRITE(6,*) ' EPR G-MATRIX'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(1,1),gtensor(1,2),gtensor(1,3),' ]'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(2,1),gtensor(2,2),gtensor(2,3),' ]'
          WRITE(6,'(T15,A,3(1X,F13.11),T64,A)') '[   ',&
               gtensor(3,1),gtensor(3,2),gtensor(3,3),' ]'
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate eigenvalues (in ppm)
    ! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent) WRITE(6,'(A,T32,A)') ' PRINCIPAL VALUES (ppm)',&
            'PRINCIPAL DIRECTIONS (degrees)'
       CALL trmm(gtensor,3)
       CALL dsyev('V','U', 3, gtensor, 3, geigval(1), dg_so, 9, ier)
       IF (ier .NE. 0) THEN
          IF (paral%io_parent) WRITE(6,*) 'ERROR IN DSYEV, code ', ier
       ELSE
          DO ii=1,3
             geigval(ii)=(SQRT(geigval(ii))-2.0023192778_real_8)*1.e6_real_8
             IF (dasin(geigvec(2,ii)).GE.0) THEN
                dg_sotemp=dacos(gtensor(1,ii)/SQRT(1-gtensor(3,ii)**2))&
                     /pi*180._real_8 ! abuse of dg_sotemp
             ELSE
                dg_sotemp=360._real_8-dacos(gtensor(1,ii)/SQRT(1-&
                     gtensor(3,ii)**2))/pi*180._real_8 ! abuse of dg_sotemp
             ENDIF
             IF (paral%io_parent)&
                  WRITE(6,'(T5,F8.1,T36,A,F6.2,A,F6.2,A)') geigval(ii),&
                  'th = ' ,dacos(gtensor(3,ii))/pi*180._real_8, '  phi = ',&
                  dg_sotemp, ' '
          ENDDO
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Some finalization 
    ! ==--------------------------------------------------------------==

    ! 444  IF (parent) CALL freem(ip_wanniercenters)
444 IF (paral%parent) DEALLOCATE(wanniercenters,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(j_a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(j_b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(BIND,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dg_so,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dg_so,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gtensor,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(geigvec,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lower_left,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lower_left,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(diffv0a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(diffv0b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoetmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rootrho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (paral%io_parent) WRITE (6,'("@ TIME: ",F8.1," SEC FOR ",A40)')&
         m_walltime()/1000._real_8 - starttime,&
         'THE WHOLE CALCULATION.      '

    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('       EPR',isub)
    RETURN
  END SUBROUTINE epr_p

  ! ==================================================================
  REAL(real_8) FUNCTION kine_one(src)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments:
    COMPLEX(real_8) :: src(ncpw%ngw)
    ! Local variables:
    INTEGER :: ig,isub
    REAL(real_8) :: kine_one_temp
    kine_one_temp = 0._real_8
    CALL tiset('KINE_ONE  ',isub)
    DO ig=1,ncpw%ngw
       kine_one_temp = kine_one_temp +&
            REAL(CONJG(src(ig))*src(ig)*hg(ig))
    ENDDO
    kine_one = kine_one_temp*parm%tpiba2
    CALL tihalt('KINE_ONE  ',isub)
    RETURN
  END FUNCTION kine_one

  ! ==================================================================
  SUBROUTINE give_scr_epr(lscr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lscr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: l_current

    lscr = 1
    IF (.NOT. response1%tepr) RETURN
    CALL give_scr_epr_current(l_current,tag)
    IF (nmr_options%tlocalize) THEN
       wannl%twann = .TRUE.
       prop1%locl  = .TRUE.
       CALL give_scr_ddipo(lscr,tag)
    ENDIF

    lscr=MAX(lscr,2*ncpw%nhg,&
         2*ncpw%ngw*crge%n,             & ! for pcgrad_nmr_p
         12*fpar%nnr1,             & ! for calc_p and calc_L
         l_current,           & ! for NMR_CURRENT
         2*maxfft)            ! for use as fft scratch

    tag  = 'epr-response          '
    RETURN
  END SUBROUTINE give_scr_epr

  ! ==================================================================
  SUBROUTINE trmm(a,n)
    ! ==--------------------------------------------------------------==
    ! Computes Symmetrical Matrix: Transpose(A)*A
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n)

    INTEGER                                  :: i, j, k
    REAL(real_8)                             :: b(n,n)

    DO i=1,n
       DO j=1,n
          b(i,j)=a(i,j)
          a(i,j)=0._real_8
       ENDDO
    ENDDO
    DO i=1,n
       DO j=1,n
          DO k=1,n
             a(i,j)=a(i,j)+b(k,i)*b(k,j)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE trmm


  ! ==================================================================
  SUBROUTINE epr_effpot(effpot,rhoe,vt,vt2)
    ! ==================================================================
    ! PARALLEL
    ! ==--------------------------------------------------------------==
    ! Arguments & common variables
    REAL(real_8)                             :: effpot(fpar%nnr1,2), &
                                                rhoe(fpar%nnr1,2)
    COMPLEX(real_8)                          :: vt(maxfft), vt2(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'epr_effpot'

    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:,:)
    INTEGER                                  :: gridmesh(ions1%nsp), i1, &
                                                ierr, ig, ig1, is
    REAL(real_8)                             :: dummy1, dummy2, vol
    REAL(real_8), ALLOCATABLE :: gridf(:,:), gridr(:,:), gridstep(:), &
      nucpot(:,:,:), potscr1(:), potscr2(:), totclbpot(:,:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL setfftn(0)
    ALLOCATE(nucpot(nsplpo,2,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(nucpot)!,2*nsplpo*maxsys%nsx)
    ALLOCATE(gridr(maxsys%mmaxx,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gridr)!,maxsys%mmaxx*ions1%nsp)
    ALLOCATE(gridf(maxsys%mmaxx,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gridf)!,maxsys%mmaxx*ions1%nsp)
    ALLOCATE(potscr1(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(potscr1)!,maxsys%mmaxx)
    ALLOCATE(potscr2(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(potscr2)!,maxsys%mmaxx)
    ALLOCATE(totclbpot(ions1%nsp,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(totclbpot)!,ions1%nsp*nhg)
    ALLOCATE(vtmp(maxfft,4),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vtmp)!,SIZE(vtmp))
    ALLOCATE(gridstep(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gridstep)!,ions1%nsp)
    ! ==--------------------------------------------------------------==
    ! Calculation of the Coulomb potential
    DO is=1,ions1%nsp
       ! Real grid (radial)
       CALL epr_pot_rg(is,gridr(1,1),gridf(1,1),gridmesh(1),&
            gridstep(1))
    ENDDO
    DO is=1,ions1%nsp
       ! Reciprocal grid (radial)
       CALL epr_pot_gg(is,nucpot(1,1,1),&
            gridr(1,1),gridf(1,1),gridmesh(1),gridstep(1),&
            potscr1(1),potscr2(1))
    ENDDO
    ! Reciprocal grid (normal)
    vol=1._real_8/parm%omega
    DO is=1,ions1%nsp
       DO ig=1,ncpw%nhg
          totclbpot(is,ig)=vol*curv2(hg(ig),nsplpo,ggnh(1),&
               nucpot(1,1,is),nucpot(1,2,is),0.0_real_8)
       ENDDO
    ENDDO
    ! Sum over all atoms
    CALL totclbpotcalc(vt2(1),totclbpot(1,1))

    ig1=1
    !$omp parallel do private(IG)
    DO ig=ig1,ncpw%nhg
       vt(nzh(ig))=vt2(ig)
       vt(indz(ig))=CONJG(vt2(ig))
    ENDDO

    CALL  invfftn(vt,.FALSE.,parai%allgrp)

    DO i1=1,fpar%nnr1
       effpot(i1,1) = REAL(vt(i1))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Calculation of the Hartree potential
    CALL zeroing(vt)!, maxfft)
    CALL zeroing(vt2)!, nhg)

    !$omp parallel do private(I1)
    DO i1=1,fpar%nnr1
       vt(i1) = CMPLX(rhoe(i1,1),0.0_real_8,kind=real_8)
    ENDDO

    CALL  fwfftn(vt,.FALSE.,parai%allgrp)

    ig1=1
    IF (geq0) THEN
       vt2(1) = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
       ig1=2
    ENDIF

    !$omp parallel do private(IG)
    DO ig=ig1,ncpw%nhg
       vt2(ig) = vt(nzh(ig)) * CMPLX(fpi/(parm%tpiba2*hg(ig)),0.0_real_8,kind=real_8)
    ENDDO

    CALL zeroing(vt)!, maxfft)

    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       vt(indz(ig)) = CONJG(vt2(ig))
       vt(nzh(ig))  = vt2(ig)
    ENDDO

    CALL  invfftn(vt,.FALSE.,parai%allgrp)

    DO i1=1,fpar%nnr1
       effpot(i1,1) = effpot(i1,1) + REAL(vt(i1))
    ENDDO

    CALL zeroing(vt)!, maxfft)
    CALL zeroing(vt2)!, nhg)

    ! ==--------------------------------------------------------------==
    ! Calculation of the XC potential
    !$omp parallel do private(I1)
    DO i1=1,fpar%nnr1
       effpot(i1,2) = effpot(i1,1)
    ENDDO

    CALL xcener(dummy1,dummy2,rhoe,rhoe,vtmp)

    !$omp parallel do private(I1)
    DO i1=1,fpar%nnr1
       effpot(i1,1) = effpot(i1,1) + REAL(vtmp(i1,1))
       effpot(i1,2) = effpot(i1,2) + REAL(vtmp(i1,2))
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(nucpot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gridr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gridf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(potscr1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(potscr2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(totclbpot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE epr_effpot

  ! ==================================================================
  SUBROUTINE epr_pot_rg(isp,pot_rg_rv,pot_rg_vr,&
       pot_rg_meshv,pot_rg_step)
    ! ==--------------------------------------------------------------==
    ! == Create one-electron all-electron potentials in real space    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    REAL(real_8) :: pot_rg_rv(maxsys%mmaxx,ions1%nsp), &
      pot_rg_vr(maxsys%mmaxx,ions1%nsp)
    INTEGER                                  :: pot_rg_meshv(ions1%nsp)
    REAL(real_8)                             :: pot_rg_step(ions1%nsp)

    INTEGER                                  :: ir, meshv
    REAL(real_8)                             :: rr1, rr2, rx, step

    step=LOG(1.025_real_8)
    meshv=LOG(7200.0_real_8*REAL(ions0%iatyp(isp),kind=real_8))/step
    IF (meshv.GT.maxsys%mmaxx) THEN
       step=REAL(meshv,kind=real_8)/REAL(maxsys%mmaxx,kind=real_8)*step
       meshv=maxsys%mmaxx
    ENDIF
    rr1=0.00625_real_8/REAL(ions0%iatyp(isp),kind=real_8)
    rx=LOG(rr1)+meshv*step
    rr2=EXP(rx)
    CALL epr_log_grid(rr1,rr2,pot_rg_rv(1,isp),meshv,step)
    DO ir=1,meshv
       rx=pot_rg_rv(ir,isp)
       pot_rg_vr(ir,isp)=-ions0%zv(isp)/rx*&
            cp_erf(rx/ownpotvalue%epr_ownpot_value)
       ! &    DERF(RX/(1.0_real_8/5.0_real_8*COVRAD(IATYP(ISP))*FBOHR))
    ENDDO
    pot_rg_meshv(isp)=meshv
    pot_rg_step(isp)=step
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE epr_pot_rg
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE epr_pot_gg(is,nucpot,pot_rg_rv,pot_rg_vr,pot_rg_meshv,&
       pot_rg_step,dfint,fint)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE COULOMB PSEUDOPOTENTIAL IN G-SPACE             ==
    ! ==  (NOT PARALLEL RIGHT NOW)                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8) :: nucpot(nsplpo,2,ions1%nsp), &
      pot_rg_rv(maxsys%mmaxx,ions1%nsp), pot_rg_vr(maxsys%mmaxx,ions1%nsp)
    INTEGER                                  :: pot_rg_meshv(ions1%nsp)
    REAL(real_8)                             :: pot_rg_step(ions1%nsp), &
                                                dfint(maxsys%mmaxx), &
                                                fint(maxsys%mmaxx)

    INTEGER                                  :: ierr, il, ir, isub
    LOGICAL                                  :: zer
    REAL(real_8)                             :: check, flx, xg, xrg

    CALL tiset('EPR_POT_GG',isub)
    CALL zeroing(fint)!,maxsys%mmaxx)
    DO ir=1,pot_rg_meshv(is)
       fint(ir)=(pot_rg_rv(ir,is)*pot_rg_rv(ir,is)*pot_rg_rv(ir,is)&
            *pot_rg_vr(ir,is)+&
            REAL(ions0%iatyp(is),kind=real_8)*cp_erf(pot_rg_rv(ir,is)/raggio(is))&
            *pot_rg_rv(ir,is)*pot_rg_rv(ir,is))
       zer = ABS(pot_rg_rv(ir,is)).LT.1.e-8_real_8
       IF (.NOT.zer) THEN
          check = pot_rg_vr(ir,is)+REAL(ions0%iatyp(is),kind=real_8)*&
               cp_erf(pot_rg_rv(ir,is)/raggio(is))/pot_rg_rv(ir,is)
       ELSE
          check = pot_rg_vr(ir,is)+2._real_8*&
               REAL(ions0%iatyp(is),kind=real_8)/(SQRT(pi)*raggio(is))
       ENDIF
       IF (ABS(check).LT.1.e-8_real_8) fint(ir)=0._real_8
    ENDDO
    DO il=1,nsplpo
       xg=SQRT(ggnh(il))*parm%tpiba
       IF (xg.GT.1.e-6_real_8) THEN
          DO ir=1,pot_rg_meshv(is)
             xrg = pot_rg_rv(ir,is)*xg
             dfint (ir) = fint(ir)* SIN(xrg)/xrg
          ENDDO
       ELSE
          CALL dcopy(pot_rg_meshv(is),fint(1),1,dfint(1),1)
       ENDIF
       CALL simpsn (pot_rg_meshv(is), dfint, flx)
       nucpot(il,1,is) = pot_rg_step(is)*fpi* flx
    ENDDO
    CALL curv1(nsplpo,ggnh,nucpot(1,1,is),0.0_real_8,0.0_real_8,3,&
         nucpot(1,2,is),fint,0.0_real_8,ierr)
    CALL tihalt('EPR_POT_GG',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE epr_pot_gg
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE totclbpotcalc(eivps,totclbpot)
    ! ==--------------------------------------------------------------==
    ! ==  derived from eicalc.F                                       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eivps(ncpw%nhg)
    REAL(real_8)                             :: totclbpot(ions1%nsp,ncpw%nhg)

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: ia, ig, is, isa, isub
    REAL(real_8)                             :: ei, er

! Variables

#ifdef __NEC
    INTEGER :: isa0
    INTEGER, ALLOCATABLE :: isv(:)

#endif
    ! ==--------------------------------------------------------------==
    CALL tiset('TOTCLBCALC',isub)
#ifdef __VECTOR
#if defined __NEC
    CALL zeroing(eivps)!,nhg)
    isa0=0
    DO is=1,ions1%nsp
       isa0 = isa0 + ions0%na(is)
    ENDDO
    ALLOCATE(isv(isa0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    isa=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          isa=isa+1
          isv(isa) = is
       ENDDO
    ENDDO
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,ISA,IS)
       DO ig=1,nhg
          DO isa = 1, isa0
             is=isv(isa)
             eivps(ig)=eivps(ig)+totclbpot(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,ISA,IS,EI123,ER,EI)
       DO ig=1,nhg
          DO isa = 1, isa0
             is=isv(isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             er=REAL(ei123)
             ei=AIMAG(ei123)
             eivps(ig)=eivps(ig)&
                  +CMPLX(er*totclbpot(is,ig),ei*totclbpot(is,ig),kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(isv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
#elif defined(_vpp_)
    CALL zeroing(eivps)!,nhg)
    IF (cntl%bigmem) THEN
       DO ig=1,nhg
          DO isa=1,ions1%nat
             is=iatpt(2,isa)
             eivps(ig)=eivps(ig)+totclbpot(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       DO ig=1,nhg
          DO isa=1,ions1%nat
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             eivps(ig)=eivps(ig)+totclbpot(is,ig)*ei123
          ENDDO
       ENDDO
    ENDIF
#else
    CALL zeroing(eivps)!,nhg)
    IF (cntl%bigmem) THEN
#ifdef __SR8000
       !poption parallel
       !voption indep(RHOPS,EIGRB,TOTCLBPOT,EIVPS)
#endif
       DO ig=1,ncpw%nhg
          isa=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                isa=isa+1
                eivps(ig)=eivps(ig)+totclbpot(is,ig)*eigrb(ig,isa)
             ENDDO
          ENDDO
       ENDDO
    ELSE
#ifdef __SR8000
       !poption parallel
       !voption indep(EI3,EI2,EI1,INYH,RHOPS,TOTCLBPOT,EIVPS)
#endif
       DO ig=1,ncpw%nhg
          isa=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                isa=isa+1
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                eivps(ig)=eivps(ig)&
                     +CMPLX(er*totclbpot(is,ig),ei*totclbpot(is,ig),kind=real_8)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
#endif
#else
    CALL zeroing(eivps)!,nhg)
    IF (cntl%bigmem) THEN
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             !$omp parallel do private(IG) shared(ISA,EIVPS)
             DO ig=1,ncpw%nhg
                eivps(ig)=eivps(ig)+totclbpot(is,ig)*eigrb(ig,isa)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             !$omp parallel do private(IG,ER,EI,EI123) shared(ISA,EIVPS)
             DO ig=1,ncpw%nhg
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                er=REAL(ei123)
                ei=AIMAG(ei123)
                eivps(ig)=eivps(ig)&
                     +CMPLX(er*totclbpot(is,ig),ei*totclbpot(is,ig),kind=real_8)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
#endif
    CALL tihalt('TOTCLBCALC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE totclbpotcalc
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE epr_log_grid(rp1,rpl,rw,mesh,clog)
    ! ==--------------------------------------------------------------==
    ! == Construct logarithmic grid from RP1 to RPL with MESH points  ==
    ! == in RW(1:MESH).                                               ==
    ! == CLOG is the step                                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rp1, rpl
    INTEGER                                  :: mesh
    REAL(real_8)                             :: rw(mesh), clog

    INTEGER                                  :: ir
    REAL(real_8)                             :: cf

    clog=LOG(rpl/rp1)/REAL(mesh-1,kind=real_8)
    cf=EXP(clog)
    rw(1)=rp1
    DO ir=2,mesh
       rw(ir)=rw(ir-1)*cf
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE epr_log_grid

END MODULE epr_p_utils
