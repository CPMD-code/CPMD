MODULE nmr_p_utils
  USE cnst,                            ONLY: fbohr
  USE coor,                            ONLY: tau0,&
                                             taup,&
                                             velp
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: localize
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sync
  USE nmr_chi_p_utils,                 ONLY: calc_chi_l,&
                                             calc_chi_p,&
                                             print_chi_tensor
  USE nmr_current_p_utils,             ONLY: calc_current,&
                                             give_scr_nmr_current
  USE nmr_full_p_utils,                ONLY: make_listofstates,&
                                             optimize_llc_indep
  USE nmr_para_p_utils,                ONLY: rewrite_output
  USE nmr_position_p_utils,            ONLY: calc_lower_left_new,&
                                             calc_lower_left_optimize,&
                                             force_xyz,&
                                             print_llc,&
                                             set_origin
  USE nmr_shift_p_utils,               ONLY: calc_l,&
                                             calc_p,&
                                             print_shift_matrix
  USE nmr_util_p_utils,                ONLY: make_123,&
                                             print_configuration,&
                                             print_wavefunctions,&
                                             printtime
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phfac_utils,                     ONLY: phfac
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE response_pmod,                   ONLY: &
       chi, chi_factor, chi_si_to_ppmcgs, chi_si_to_shift_ppm, &
       firstshelldistance, inmr_csgt, inmr_default, inmr_iglo, inmr_llc, &
       inmr_novirt, inmr_wc, lower_left, lower_left_value, nmr_options, &
       nmr_para, nris, response1, response2, response_read, response_write, &
       shift_factor, shift_matrix, timetag, wanniercenters
  USE restart_p_utils,                 ONLY: restart_nmr,&
                                             restart_p
  USE ropt,                            ONLY: iteropt
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE setirec_utils,                   ONLY: isetone,&
                                             write_irec
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: irec_pwv,&
                                             store1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE wann,                            ONLY: wannl,&
                                             wannr
  USE wv30_utils,                      ONLY: zhwwf
!!use nmr_full_p_utils, only : do_full
!!use nmr_position_p_utils, only : apply_op_rx
!!use nmr_util_p_utils, only : fft2tor
!!use nmr_util_p_utils, only : apply_op_p
!!use nmr_util_p_utils, only : print_orbital_densities
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nmr_p
  PUBLIC :: give_scr_nmr

CONTAINS

  ! ==================================================================
  SUBROUTINE nmr_p(c0,h1psi0,psi,rhoe,eirop,eivps,z11,nstate)
    ! ==================================================================
    ! PARALLEL
    ! ==--------------------------------------------------------------==
    ! Arguments & common variables
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: h1psi0(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'nmr_p'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: xyz = (/'x','y','z'/)

    CHARACTER(len=10)                        :: tag
    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: c1(:,:,:), cs(:), &
                                                gradpsi(:,:), sc0(:), scr(:,:)
    INTEGER :: i_c1, i_state, i_thread, iB, ierr, ii, iii, ik, iL_end, &
      iL_start, ip_end, ip_start, irec(100), is, isub, isub2, l_c1, nc
    INTEGER, ALLOCATABLE                     :: neighborlist(:,:)
    LOGICAL                                  :: full_optimization, &
                                                simple_done(6)
    LOGICAL, ALLOCATABLE                     :: dowfopt(:,:), full_done(:,:)
    REAL(real_8)                             :: rdummy(1), starttime, &
                                                usedtime(0:6), w_eps_save
    REAL(real_8), ALLOCATABLE :: chi_saved(:,:), chi_sumup(:,:), eigv(:,:), &
      shift_saved(:,:,:,:), shift_sumup(:,:,:,:)

81  FORMAT ("*  TOTAL CORRECTION :     ",38x,"*")
82  FORMAT ("*  RESULT OF FULL CALCULATION:   ",31x,"*")
96  FORMAT ("*  PERTURBATION: ",a1,"_",i1," |PSI>",38x,"*")
99  FORMAT ("*  Calculation of responses DONE.",31x,"*")
97  FORMAT ("*  TOTAL TIME: ",f7.1,42x,"*")
98  FORMAT (52("*"),"NMR*RESPONSE*")
    ! ==--------------------------------------------------------------==
    CALL tiset('       NMR',isub)
    ALLOCATE(scr(maxfft,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    timetag=' '
    IF (paral%parent) CALL printtime
    starttime=m_walltime()/1000._real_8

    nris(1)   = spar%nr1s
    nris(2)   = spar%nr2s
    nris(3)   = spar%nr3s
    response1%pr_energy = nmr_options%tverbose      ! these values are not sooo interesting
    full_optimization=.FALSE.
    ! full_optimization=.true.
    ! if (parent) WRITE(6,*)
    ! &     'DEBUG DEBUG: >> full_optimization: ',full_optimization
    ! ==--------------------------------------------------------------==
    IF (nmr_para%nmr_threads .EQ. 1) THEN
       ip_start=1
       ip_end  =3
       iL_start=1
       iL_end  =3
    ELSEIF (nmr_para%nmr_threads .EQ. 3) THEN
       ip_start=nmr_para%nmr_mygroup+1
       ip_end  =nmr_para%nmr_mygroup+1
       iL_start=nmr_para%nmr_mygroup+1
       iL_end  =nmr_para%nmr_mygroup+1
    ELSEIF (nmr_para%nmr_threads .EQ. 6) THEN
       ip_start=nmr_para%nmr_mygroup+1
       ip_end  =MOD(nmr_para%nmr_mygroup,3)+1
       iL_start=MOD(nmr_para%nmr_mygroup,3)+1
       iL_end  =nmr_para%nmr_mygroup-2
    ELSE
       CALL stopgm('NMR_P','Bad number of threads! Must be 3 or 6',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    l_c1 = 2*ncpw%ngw*nstate
    IF (nmr_options%tcurrent)    l_c1 = 6 * 2*ncpw%ngw*nstate
    IF (nmr_options%tnmr_full)  l_c1 = 6 * 2*ncpw%ngw*nstate
    IF (nmr_options%tcurrent .AND. nmr_options%tnmr_full) l_c1 = 9 * 2*ncpw%ngw*nstate
    ! call memory(ip_shift_matrix,maxsys%nax*maxsys%nsx*3*3,'shifts')
    ALLOCATE(shift_matrix(3,3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   shift_matrix)!,maxsys%nax*maxsys%nsx*3*3)
    ! call memory(ip_chi,3*3,'suscept')
    ALLOCATE(chi(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   chi)!,3*3)
    ! call memory(ip_c1,       l_c1,'C1')
    ALLOCATE(c1(ncpw%ngw,nstate,l_c1/(ncpw%ngw*nstate)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1)!,SIZE(c1))
    ! call memory(ip_lower_left,3*nstate,'lower left')
    ALLOCATE(lower_left(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   lower_left)!,3*nstate)
    ! call memory(ip_lower_left_value,3*nstate,'llc-value')
    ALLOCATE(lower_left_value(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   lower_left_value)!,3*nstate)
    ! call memory(ip_shift_sumup,3*3*maxsys%nax*maxsys%nsx,'SHIFT/sumup')
    ALLOCATE(shift_sumup(3,3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   shift_sumup)!,3*3*maxsys%nax*maxsys%nsx)
    ! call memory(ip_shift_saved,3*3*maxsys%nax*maxsys%nsx,'SHIFT/saved')
    ALLOCATE(shift_saved(3,3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   shift_saved)!,3*3*maxsys%nax*maxsys%nsx)
    ! call memory(ip_chi_sumup,3*3,'SUSC/sumup')
    ALLOCATE(chi_sumup(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   chi_sumup)!,3*3)
    ! call memory(ip_chi_saved,3*3,'SUSC/saved')
    ALLOCATE(chi_saved(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   chi_saved)!,3*3)
    ! call memory(ip_full_done,3*nstate,'done-matrix')
    ALLOCATE(full_done(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gradpsi(ncpw%ngw,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    DO iB=1,3
       simple_done(iB)=.FALSE.
       simple_done(iB+3)=.FALSE.
       DO is=1,nstate
          full_done(iB,is) = .FALSE.
       ENDDO
    ENDDO
    IF (paral%parent) THEN
       ! call memory(ip_wanniercenters,4*nstate,'WANNIERCENT')
       ALLOCATE(wanniercenters(4,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   wanniercenters)!,4*nstate)
    ENDIF

    IF (paral%parent) CALL prmem('NMR        ')
    ! ==--------------------------------------------------------------==
    ! !    shift_factor= - 2._real_8*( 1._real_8 / 137.03602_real_8 )**2   * 1.e6_real_8 / omega
    ! factor for the CHEMICAL SHIFTS: alpha^2 *  ppm.
    ! 2: double-occupation of the orbitals c0  
    shift_factor=- (1.0_real_8 / 137.03602_real_8 )**2  * 1.0e6_real_8 / parm%omega
    ! The orbital occupations are explicitly included

    ! !    chi_factor = 2._real_8 * 1.9727566e-29_real_8 / 1.e-30_real_8 ! -> displayed in 10^-30 J/T^2
    ! = (2:double occ) * 1/4 * e^2/m * a_0 ^2
    chi_factor = 1.9727566e-29_real_8 / 1.e-30_real_8
    ! The orbital occupations are explicitly included

    chi_SI_to_ppmcgs   = 6.022045_real_8 / 1.e2_real_8
    ! converts 10^-30 J/T^2  into  ppm cgs = ppm cm^3/mol
    ! = 10^-30  *  mu_0/4pi  *  N_A  *  10^6 * 10^6
    ! [one 10^6 for ppm, one for m^3 -> cm^3]

    chi_SI_to_shift_ppm = 1.e-30_real_8 * 8.37758041e-7_real_8&
         / (parm%omega / (fbohr*fbohr*fbohr) *1.0e-30_real_8) * 1.0e6_real_8
    ! 10^-30  *  2/3  mu_0 / Omega  * 1/ppm
    firstshelldistance = 1._real_8 ! in a.u.
    ! ==--------------------------------------------------------------==
    IF (nmr_options%inmr_method.EQ.inmr_default) THEN
       nmr_options%inmr_method = inmr_csgt
       IF (isos1%tclust) THEN
          nmr_options%inmr_virtual = inmr_novirt
       ENDIF
    ENDIF
    IF (nmr_options%inmr_virtual.EQ.inmr_default) THEN
       nmr_options%inmr_virtual = inmr_llc
       IF (nmr_options%inmr_method.EQ.inmr_iglo) nmr_options%inmr_virtual = inmr_wc
    ENDIF
    IF (nmr_options%inmr_method.EQ.inmr_iglo) THEN
       CALL stopgm('NMR','IGLO no longer supported',& 
            __LINE__,__FILE__)
    ENDIF
    IF (nmr_options%inmr_virtual.EQ.inmr_novirt) nmr_options%tlocalize = .FALSE.
    IF (response1%t_restart_c1 .AND. cntl%wfopt) THEN
       CALL stopgm('NMR/RESTART',&
            'Use of NOOPT is compulsory for NMR RESTART',& 
            __LINE__,__FILE__)
    ENDIF
    IF (nmr_para%nmr_superparent) CALL print_configuration
    ! ==--------------------------------------------------------------==
    ! ==  LOCALIZATION/C0  -------------------------------------------==
    CALL phfac(tau0)
    CALL eicalc(eivps,eirop)
    IF (nmr_options%tlocalize) THEN
       wannl%twann=.TRUE.
       prop1%locl=.TRUE.
       prop2%numorb=MAX(prop2%numorb,nstate)
       nc=2*nkpt%ngwk*prop2%numorb
       lenbk=nxxfun(prop2%numorb)
       nc=MAX(2*lenbk*parai%nproc,nc)
       IF (.NOT. response1%t_restart_c1) THEN
          ALLOCATE(sc0(nc),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cs(nc),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          w_eps_save = wannr%w_eps
          wannr%w_eps      = 1.e-4_real_8
          CALL localize(tau0,c0,cs,sc0,nstate)
          wannr%w_eps      = w_eps_save
          store1%swf        = .TRUE.
          CALL write_irec(irec)
          irec(irec_pwv) = isetone(.FALSE.)
          ALLOCATE(eigv(nstate,nkpt%nkpts),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          IF (.NOT.cntl%tmdfile) THEN
             ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
          ENDIF
          CALL zhwwf(49,irec,c0,cs,nstate,eigv,tau0,velp,taup,iteropt%nfi)
          IF (.NOT.cntl%tmdfile) THEN
             DEALLOCATE(taup,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
          DEALLOCATE(eigv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(sc0,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cs,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    timetag='localization.'
    IF (paral%parent) CALL printtime

    ! ==--------------------------------------------------------------==
    CALL lag_mult(c0,h1psi0,psi,rhoe,z11,nstate)
    timetag='lagrange multipliers.'
    IF (paral%parent) CALL printtime
    ! ==--------------------------------------------------------------==
    IF (nmr_options%tprintwfns)&
         CALL print_wavefunctions(c0,nstate,psi(:,1))
    IF (nmr_options%tprintrho)&
         CALL print_orbital_densities(c0,nstate,psi,scr)
    ! ==--------------------------------------------------------------==
    IF (nmr_options%inmr_virtual.EQ.inmr_llc) THEN
       IF (nmr_options%tnmr_overlaps) THEN
          CALL calc_lower_left_optimize(c0,c1,nstate)
       ELSE
          CALL calc_lower_left_new(c0,nstate)
          IF (full_optimization) THEN
             ! TESTING: LLC optimization...
             ! call memory(ip_neighborlist,3*nstate,'NMR:NLIST')
             ALLOCATE(neighborlist(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ! call memory(ip_dowfopt,3*nstate,'NMR:CLIST')
             ALLOCATE(dowfopt(3,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)

             firstshelldistance= MIN (1.5_real_8, parm%alat/20._real_8)
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(6,*)'PROCEDURE:'
                IF (paral%io_parent)&
                     WRITE(6,*)' new LLC opti / 1 x ',&
                     firstshelldistance,' a.u.'
             ENDIF

             CALL optimize_llc_indep(neighborlist,dowfopt,nstate)

             ! firstshelldistance=1.0_real_8
             ! call make_listofstates(neighborlist,dowfopt,nstate)
             ! call optimize_llc(neighborlist,dowfopt,nstate)
             ! call make_listofstates(neighborlist,dowfopt,nstate)
             ! call optimize_llc(neighborlist,dowfopt,nstate)
             ! call make_listofstates(neighborlist,dowfopt,nstate)
             ! call optimize_llc(neighborlist,dowfopt,nstate)

             firstshelldistance=0.1_real_8
             CALL make_listofstates(neighborlist,dowfopt,nstate)

             ! call freem(ip_neighborlist)
             DEALLOCATE(neighborlist,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             ! call freem(ip_dowfopt)
             DEALLOCATE(dowfopt,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ENDIF         ! full_optimization
       ENDIF
    ELSEIF (nmr_options%inmr_virtual.EQ.inmr_novirt) THEN
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
    IF (response1%t_restart_c1) CALL restart_nmr(simple_done,full_done,&
         shift_sumup,chi_sumup,&
         shift_matrix,chi,nstate,response_read)
    ! ==--------------------------------------------------------------==
    i_c1 = 1                  ! Index for the C1-array. Constant unless
    DO iB=iL_start,iL_end     ! calculation of the CURRENT is desired.
       IF (.NOT. simple_done(iB)) THEN
          IF (nmr_options%tcurrent) i_c1 = iB+3
          IF (nmr_options%tcurrent .AND. nmr_options%tnmr_full) i_c1 = iB+6
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               WRITE (6,96) 'L',iB
          ! Preparation of L_i |psi_0>:
          CALL make_123(iB,ii,iii)
          DO is=1,nstate! apply r_ii p_iii - r_iii p_ii:
             CALL set_origin (is)
             CALL apply_op_p (c0(:,is), scr(:,1), ii, ncpw%ngw)
             CALL apply_op_p (c0(:,is), scr(:,2), iii, ncpw%ngw)
             CALL fft2tor    (scr(:,1), scr(:,1),&
                  scr(:,2),scr(:,2), psi, ncpw%ngw,.FALSE.)
             CALL apply_op_rx(scr(:,1), -1._real_8, rhoe,  0._real_8, iii)
             CALL apply_op_rx(scr(:,2),  1._real_8, rhoe,  1._real_8,  ii)
             CALL ffttog     (rhoe, h1psi0(1,is), psi, ncpw%ngw,.FALSE.)
             CALL dscal(2*ncpw%ngw,crge%f(is,1),h1psi0(1,is),1)
          ENDDO


          timetag='preparation.'
          IF (paral%parent) CALL printtime
          IF (paral%io_parent)&
               WRITE (6,98)
          CALL tiset('PSI1-WF-OPT',isub2)
          CALL rwfopt_p(c0,c1(1,1,i_c1),psi,rhoe,rdummy,&
               eirop,eivps,cdummy,h1psi0,&
               z11,nstate,cdummy)
          CALL tihalt('PSI1-WF-OPT',isub2)
          timetag='calculation of psi^1.'
          IF (paral%parent) CALL printtime
          IF (soft_com%exsoft) GOTO 441
          IF (paral%io_parent)&
               WRITE (6,98)
          CALL calc_l(c0,c1(1,1,i_c1),nstate,iB,psi(:,1),rhoe)
          timetag='shift contribution.'
          IF (paral%parent) CALL printtime
          IF (nmr_options%tverbose) CALL print_shift_matrix(nstate,.FALSE.)
          CALL calc_chi_l(c0,c1(1,1,i_c1),nstate,iB,psi(:,1))
          timetag='susceptibility contribution.'
          IF (paral%parent) CALL printtime

          simple_done(iB) = .TRUE.
          CALL restart_nmr(simple_done,full_done,&
               shift_sumup,chi_sumup,&
               shift_matrix,chi,nstate,response_write)
       ENDIF               ! if simple calculation not already done
    ENDDO                     ! iB
    ! ==--------------------------------------------------------------==
    DO iB=ip_start,ip_end
       IF (nmr_options%tcurrent) i_c1 = iB
       IF (nmr_options%tnmr_full) i_c1 = iB
       IF (paral%io_parent)&
            WRITE (6,98)
       IF (paral%io_parent)&
            WRITE (6,96) 'p',iB

       tag = 'p_'//xyz(iB)
       IF (response1%t_restart_c1 .AND. nmr_options%tnmr_full .AND. simple_done(iB+3))&
            CALL restart_p(c1(1,1,i_c1),tag,nstate,response_read)
       IF (.NOT. simple_done(iB+3)) THEN
          DO is=1,nstate
             CALL apply_op_p(c0(1,is), h1psi0(1,is), iB, ncpw%ngw)
             CALL dscal(2*ncpw%ngw,crge%f(is,1),h1psi0(1,is),1)
          ENDDO
          timetag='preparation.'
          IF (paral%parent) CALL printtime
          IF (paral%io_parent)&
               WRITE (6,98)
          CALL tiset('PSI1-WF-OPT',isub2)
          CALL rwfopt_p(c0,c1(1,1,i_c1),psi,rhoe,rdummy,&
               eirop,eivps,cdummy,h1psi0,&
               z11,nstate,cdummy)
          CALL tihalt('PSI1-WF-OPT',isub2)
          IF (paral%io_parent)&
               WRITE (6,98)
          timetag='calculation of psi^1.'
          IF (paral%parent) CALL printtime
          IF (soft_com%exsoft) GOTO 441
          CALL calc_p(c0,c1(1,1,i_c1),nstate,iB,psi(:,1))
          timetag='shift contribution.'
          CALL printtime
          CALL calc_chi_p(c0,c1(1,1,i_c1),nstate,iB,psi(:,1),gradpsi)!vw remplace --rhoe-- with --gradpsi--
          timetag='susceptibility contribution.'
          CALL printtime

          tag = 'p_'//xyz(iB)
          CALL restart_p(c1(1,1,i_c1),tag,nstate,response_write)
          simple_done(iB+3) = .TRUE.
          CALL restart_nmr(simple_done,full_done,&
               shift_sumup,chi_sumup,&
               shift_matrix,chi,nstate,response_write)
       ENDIF               ! if simple calculation not already done
       IF (nmr_options%tverbose) CALL print_shift_matrix(nstate,.FALSE.)
    ENDDO                     ! iB
441 CONTINUE

    CALL zeroing(usedtime)!(0),6)
    usedtime(nmr_para%nmr_mygroup) =m_walltime()/1000._real_8 - starttime
    ! ==--------------------------------------------------------------==
    ! SUPERPARALLEL stuff.
    IF (nmr_para%nmr_threads .NE. 1) THEN
       CALL rewrite_output   ! closes, prints & deletes the files out_?

       IF ((nmr_para%nmr_superparent).AND.paral%io_parent)&
            WRITE(6,98)
       DO i_thread=1,nmr_para%nmr_threads
          CALL mp_sync(nmr_para%nmr_supergroup)
          IF (parai%me.EQ.nmr_para%parents(i_thread)) THEN
             IF (paral%io_parent)&
                  WRITE (6,'("@ TIME: ",F8.1,'//&
                  '" sec times ",I3," procs (grp ",I1,") = ",'//&
                  'I7," CPU-SEC.")')&
                  usedtime(i_thread-1),parai%nproc,i_thread-1,&
                  NINT(usedtime(i_thread-1)*parai%nproc)
          ENDIF
       ENDDO

       parai%allgrp = nmr_para%nmr_supergroup
       paral%parent = nmr_para%nmr_superparent
       parai%nproc  = nmr_para%nmr_total_nproc
       parai%source = nmr_para%nmr_supersource
    ENDIF
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,98)
       IF (paral%io_parent)&
            WRITE (6,98)
       IF (paral%io_parent)&
            WRITE (6,99)
       IF (paral%io_parent)&
            WRITE (6,97) usedtime(nmr_para%nmr_mygroup)
       IF (paral%io_parent)&
            WRITE (6,98)
       IF (paral%io_parent)&
            WRITE (6,98)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL print_chi_tensor     ! calculates also the suscept correction,
    ! chi_SI_iso
    CALL print_shift_matrix(nstate,.TRUE.)
    ! ==--------------------------------------------------------------==
    IF (soft_com%exsoft) GOTO 444
    IF (nmr_options%tnmr_full) THEN
       IF (nmr_options%inmr_method .NE. inmr_csgt) CALL stopgm('NMR_P',&
            'Cannot use FULL option in IGLO calculation.',& 
            __LINE__,__FILE__)
       IF (nmr_para%nmr_threads .NE. 1) CALL stopgm('NMR_P',&
            'Cannot use superparallel option in FULL calculation.',& 
            __LINE__,__FILE__)
       CALL dcopy(3*3*maxsys%nax*maxsys%nsx,shift_matrix,1,shift_saved,1)
       CALL dcopy(3*3,chi,1,chi_saved,1)
       CALL zeroing(   shift_matrix)!,3*3*maxsys%nax*maxsys%nsx)
       CALL zeroing(   chi)!,          3*3)

       ! call memory(ip_neighborlist,3*nstate,'NMR:NLIST')
       ALLOCATE(neighborlist(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! call memory(ip_dowfopt,3*nstate,'NMR:CLIST')
       ALLOCATE(dowfopt(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL make_listofstates(neighborlist,dowfopt,nstate)
       ! neighborlist now contains the states in a reasonable 
       ! proximity-determined order, and dowfopt indicates whether or not
       ! rwfopt_p should be called for that state.

       IF (.NOT. full_optimization) THEN
          DO is=1,nstate
             DO iB=1,3
                dowfopt(iB,is) = .TRUE.
             ENDDO
          ENDDO
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** UNCONDITIONAL  RWFOPT_P:    *** *** ***'
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
          ENDIF
       ELSE
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** ONLY SELECTIVE RWFOPT_P:    *** *** ***'
             IF (paral%io_parent)&
                  WRITE(6,*)'*** *** *** *** *** *** *** *** *** *** *** ***'
          ENDIF
       ENDIF

       response2%tolog_p = response2%tolog_p * 1.e2_real_8
       ! 100 is a large number, but it seems to work...
       response1%cg_analytic = 0

       DO iB=1,3
          DO is=1,NSTATE      ! "is" is only an INDEX!!!
             ! The state is i_state !
             i_state = neighborlist(iB,is)
             IF (.NOT. full_done(iB,i_state)) THEN
                CALL do_full(i_state,iB,&
                     dowfopt(iB,i_state),&
                     h1psi0,c0,c1,scr,psi,rhoe,&
                     eirop,eivps,z11,nstate)
                IF (soft_com%exsoft) GOTO 442

                CALL daxpy(3*3,1._real_8,chi,1,chi_sumup,1)
                CALL zeroing(chi)!,3*3)
                CALL daxpy(3*3*maxsys%nax*maxsys%nsx, 1._real_8,&
                     shift_matrix,1,shift_sumup,1)
                CALL zeroing(shift_matrix)!,3*3*maxsys%nax*maxsys%nsx)

                full_done(iB,i_state) = .TRUE.
                CALL restart_nmr(simple_done,full_done,&
                     shift_sumup,chi_sumup,&
                     shift_saved,chi_saved,nstate,response_write)
             ENDIF
          ENDDO         ! i_state = 1,...,nstate.
       ENDDO               ! iB = 1,2,3
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               WRITE (6,99)
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               WRITE (6,*)
       ENDIF

442    CALL dcopy(3*3,chi_sumup,1,chi,1)
       ! call print_chi_tensor   ! calculates also chi_SI_iso
       CALL dcopy(3*3*maxsys%nax*maxsys%nsx,shift_sumup,1,shift_matrix,1)
       ! call print_shift_matrix(nstate,.false.)

       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,*)
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               WRITE (6,82)
          IF (paral%io_parent)&
               WRITE (6,98)
       ENDIF
       CALL dcopy(3*3,    chi_sumup,1,chi,1)
       CALL daxpy(3*3,1._real_8,chi_saved,1,chi,1)
       CALL dcopy(3*3*maxsys%nax*maxsys%nsx,    shift_sumup,1,shift_matrix,1)
       CALL daxpy(3*3*maxsys%nax*maxsys%nsx,1._real_8,shift_saved,1,shift_matrix,1)
       CALL print_chi_tensor ! calculates also chi_SI_iso
       CALL print_shift_matrix(nstate,.TRUE.)


       ! call freem(ip_neighborlist)
       DEALLOCATE(neighborlist,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! call freem(ip_dowfopt)
       DEALLOCATE(dowfopt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! call freem(ip_shift_sumup)
       DEALLOCATE(shift_sumup,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! call freem(ip_shift_saved)
       DEALLOCATE(shift_saved,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! call freem(ip_chi_sumup)
       DEALLOCATE(chi_sumup,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! call freem(ip_chi_saved)
       DEALLOCATE(chi_saved,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

    ENDIF                     ! full calculation
    ! ==--------------------------------------------------------------==
    IF (nmr_options%tcurrent) THEN
       IF (nmr_para%nmr_threads.NE.1)&
            CALL stopgm('NMR_P:','CURRENT NOT SUPERPARALLEL',& 
            __LINE__,__FILE__)
       CALL calc_current(c0,c1,nstate,psi(:,1))
       timetag='current density calculation. '
       IF (paral%parent) CALL printtime
    ENDIF


444 DEALLOCATE(shift_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! call freem(ip_chi)
    DEALLOCATE(chi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! call freem(ip_lower_left_value)
    DEALLOCATE(lower_left_value,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! call freem(ip_lower_left)
    DEALLOCATE(lower_left,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! call freem(ip_c1)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) DEALLOCATE(wanniercenters,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gradpsi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (paral%io_parent)&
         WRITE (6,'("@ TIME: ",F8.1," sec for ",A40)')&
         m_walltime()/1000._real_8 - starttime,&
         'the whole calculation.      '

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('       NMR',isub)
    RETURN
  END SUBROUTINE nmr_p
  ! ==================================================================
  SUBROUTINE give_scr_nmr(lscr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lscr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: l_current

    lscr = 1
    IF (.NOT. response1%tnmr) RETURN
    CALL give_scr_nmr_current(l_current,tag)
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

    tag  = 'nmr-response          '
    RETURN
  END SUBROUTINE give_scr_nmr
  ! ==================================================================


END MODULE nmr_p_utils
