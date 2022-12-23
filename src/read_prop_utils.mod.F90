MODULE read_prop_utils
  USE efld,                            ONLY: textfld
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr
  USE td_input,                        ONLY: &
       bohm, curr_step, gaugefield, masklog, maskreal, pointcharge, ppcorr, &
       real_8, tberry, tcurr, td_prop, ttransg
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_prop

CONTAINS

  ! ==================================================================
  SUBROUTINE read_prop
    ! ==--------------------------------------------------------------==
    ! 
    ! Variables
    CHARACTER(len=100)                       :: line
    INTEGER                                  :: ierr, iunit

! 

    IF (paral%io_parent)&
         WRITE(6,'(/,a,a)')&
         ' ********************** PROPAGATION TD-',&
         'DFT ********************** '
    ! ==--------------------------------------------------------------==
    ! ==-  read input parameters for cntl%md run                          -==
    ! ==--------------------------------------------------------------==
    IF (.NOT.paral%io_parent) GOTO 9999

    iunit = 5
    ierr=inscan(iunit,'&PTDDFT')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' read_prop input section &PTDDFT not found'
       CALL stopgm('read_prop',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (cntl%tmdeh) THEN
       td_prop%nzeros    = 6
       td_prop%tintrvll  = cntr%delt_elec
       td_prop%n_cycles  = 1
       td_prop%epsil     = 1.0e-10_real_8
       td_prop%ampl      = 1.0e-15_real_8
       td_prop%pertdir   = 1
       td_prop%pert_type = 1
       td_prop%read_wf   = 0
       td_prop%tdfreq    = 0.0
       td_prop%stextpot  = .FALSE.
       td_prop%td_extpot  = .FALSE.
       td_prop%tpointch  = .FALSE.
       masklog%tmask     =.FALSE.
       pointcharge%pcht0     = 0._real_8
       pointcharge%pchdt     = 0._real_8
       pointcharge%pchint    = 0._real_8
       pointcharge%pchfx     = 0._real_8
       pointcharge%pchfy     = 0._real_8
       pointcharge%pchfz     = 0._real_8
       pointcharge%pchwr     = 0
       gaugefield%gswtt     = 250
       gaugefield%gdfreq    = 0.00025
       gaugefield%gimax     = 0
       gaugefield%nperiod   = 5
       gaugefield%pi_pulse  = .FALSE.
       gaugefield%muij      = 1.0
       gaugefield%tkgauge   = .FALSE.
       gaugefield%kgauge    = 0
       tcurr     = .FALSE.
       curr_step = 100
       cntl%start_real=.FALSE.
       bohm%tgaussian_wp=.FALSE.
       bohm%gw_par1   =1.0
       bohm%gw_par2   =1.0
       bohm%tplanewave_wp=.FALSE.
       bohm%nf_bohmtraj=1
       bohm%tbohmtraj=.FALSE.
       bohm%restart_bohm=.FALSE.
       bohm%eh_restart_n=0
       bohm%eh_restart_dens=.FALSE.
       td_prop%ionize_from_state=.FALSE.
       td_prop%td_fix_spin_dens=.FALSE.
       td_prop%ionize_from_rdm=.FALSE.
       td_prop%do_not_normalize=.FALSE.
    ELSEIF (cntl%tpspec.OR.cntl%tpdist) THEN
       td_prop%nzeros    = 6
       td_prop%tintrvll  = 0.05
       td_prop%n_cycles  = cnti%nomore
       td_prop%epsil     = 1.0e-10_real_8
       td_prop%ampl      = 1.0e-15_real_8
       td_prop%pertdir   = 1
       td_prop%pert_type = 1
       td_prop%read_wf   = 0
       td_prop%tdfreq    = 0.0
       td_prop%stextpot  = .FALSE.
       td_prop%td_extpot  = .FALSE.
       td_prop%tpointch  = .FALSE.
       pointcharge%pcht0     = 0._real_8
       pointcharge%pchdt     = 0._real_8
       pointcharge%pchint    = 0._real_8
       pointcharge%pchfx     = 0._real_8
       pointcharge%pchfy     = 0._real_8
       pointcharge%pchfz     = 0._real_8
       pointcharge%pchwr     = 0
       gaugefield%gswtt     = 250
       gaugefield%gdfreq    = 0.00025
       gaugefield%gimax     = 0
       gaugefield%nperiod   = 5
       gaugefield%pi_pulse  = .FALSE.
       gaugefield%muij      = 0.0
       gaugefield%tkgauge   = .FALSE.
       gaugefield%kgauge    = 0
       tcurr     = .FALSE.
       curr_step = 100
       cntl%start_real=.FALSE.
       bohm%tgaussian_wp=.FALSE.
       td_prop%ionize_from_state=.FALSE.
       td_prop%td_fix_spin_dens=.FALSE.
       td_prop%ionize_from_rdm=.FALSE.
       td_prop%do_not_normalize=.FALSE.
    ENDIF
    ppcorr=.FALSE.
    tberry=.FALSE.
    ttransg   =.FALSE.

113 CONTINUE
    IF (paral%io_parent)&
         READ(iunit,*,END=144) line
    IF (INDEX(line,'&END').NE.0) GOTO 144

    IF (INDEX(line,"NUMBER_ZEROS").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%nzeros
    ENDIF
    IF (INDEX(line,"ACCURACY").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%epsil
    ENDIF
    IF (INDEX(line,"PROP_TSTEP").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%tintrvll
    ENDIF
    IF (INDEX(line,"N_CYCLES").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%n_cycles
    ENDIF
    IF (INDEX(line,"EXT_POTENTIAL").NE.0)THEN
       td_prop%stextpot = .TRUE.
    ENDIF
    IF (INDEX(line,"EXT_PULSE").NE.0)THEN
       td_prop%stextpot = .TRUE.
       td_prop%pert_type=1
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%ampl
    ENDIF
    IF (INDEX(line,"TD_POTENTIAL").NE.0)THEN
       td_prop%td_extpot = .TRUE.
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%tdfreq
    ENDIF
    IF (INDEX(line,"POINT_CHARGE").NE.0)THEN
       td_prop%tpointch = .TRUE.
       IF (paral%io_parent)&
            READ(iunit,*) pointcharge%pcht0,pointcharge%pchdt,pointcharge%pchint,pointcharge%pchfx, &
            & pointcharge%pchfy,pointcharge%pchfz,pointcharge%pchwr
    ENDIF
    IF (INDEX(line,"PERT_TYPE").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%pert_type
    ENDIF
    IF (INDEX(line,"PERT_AMPLI").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%ampl
    ENDIF
    IF (INDEX(line,"PERT_DIRECTION").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%pertdir
    ENDIF
    IF (INDEX(line,"RESTART").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) td_prop%read_wf
    ENDIF
    IF (INDEX(line,"PIPULSE").NE.0)THEN
       gaugefield%pi_pulse = .TRUE.
       IF (paral%io_parent)&
            READ(iunit,*) gaugefield%muij
    ENDIF
    IF (INDEX(line,"GAUGEFIELD PARA").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) gaugefield%gswtt,gaugefield%gdfreq,gaugefield%gimax
    ENDIF
    IF (INDEX(line,"GAUGEPULSE PARA").NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,*) gaugefield%nperiod
    ENDIF
    IF (INDEX(line,"GAUGE_VECTOR").NE.0)THEN
       gaugefield%tkgauge=.TRUE.
       IF (paral%io_parent)&
            READ(iunit,*) gaugefield%kgauge
    ENDIF
    IF (INDEX(line,"PPCORR").NE.0)THEN
       ppcorr=.TRUE.
    ENDIF
    IF (INDEX(line,"BERRY").NE.0)THEN
       tberry=.TRUE.
    ENDIF
    IF (INDEX(line,"CAL_CURRENT").NE.0)THEN
       tcurr=.TRUE.
       IF (paral%io_parent)&
            READ(iunit,*) curr_step
    ENDIF
    IF(INDEX(line,"ADD_MASK").NE.0)THEN                      
       masklog%tmask=.TRUE.                                            
       IF(paral%io_parent) READ(iunit,*) maskreal%maskpar1,maskreal%maskpar2,maskreal%maskpar3  
    ENDIF
    IF(INDEX(line,"GAUSSIAN_WP").NE.0)THEN                   
       bohm%tgaussian_wp=.TRUE.                                     
       IF(paral%io_parent) READ(iunit,*) bohm%gw_par1,bohm%gw_par2             
       IF(paral%io_parent) READ(iunit,*) bohm%gw_r1,bohm%gw_r2,bohm%gw_r3           
       IF(paral%io_parent) READ(iunit,*) bohm%gw_p1,bohm%gw_p2,bohm%gw_p3           
    ENDIF
    IF(INDEX(line,"PLANE_WAVE").NE.0)THEN                    
       bohm%tplanewave_wp=.TRUE.                                    
       IF(paral%io_parent) READ(iunit,*) bohm%gw_r1,bohm%gw_r2,bohm%gw_r3           
       IF(paral%io_parent) READ(iunit,*) bohm%gw_p1,bohm%gw_p2,bohm%gw_p3           
    ENDIF
    IF(INDEX(line,"BOHMIANTRAJ").NE.0)THEN                   
       bohm%tbohmtraj=.TRUE.                                        
       IF(paral%io_parent) READ(iunit,*) bohm%nf_bohmtraj                 
    ENDIF
    IF(INDEX(line,"BOHM_R").NE.0)THEN                        
       bohm%restart_bohm=.TRUE.                                     
    ENDIF
    IF(INDEX(line,"DENS_R").NE.0)THEN                        
       bohm%eh_restart_dens=.TRUE.                                  
       IF(paral%io_parent) READ(iunit,*) bohm%eh_restart_n                
    ENDIF
    IF(INDEX(line,"IONIZE_FROM_STATE").NE.0)THEN
       td_prop%ionize_from_state=.TRUE.
       cntl%start_real=.TRUE.
       IF(paral%io_parent) READ(iunit,*) td_prop%ionized_state
    ENDIF
    IF (INDEX(line,"IONIZE_FROM_RDM").NE.0)THEN
       td_prop%ionize_from_rdm=.TRUE.
       cntl%start_real=.TRUE.
    ENDIF
    IF(INDEX(line,"FIX_SPIN_DENS").NE.0)THEN
       td_prop%td_fix_spin_dens=.TRUE.
       cntl%start_real=.TRUE.
    ENDIF
    IF(INDEX(line,"DO_NOT_NORMALIZE").NE.0)THEN
       td_prop%do_not_normalize=.TRUE.
    ENDIF

    GOTO 113

144 CONTINUE

    ! ==--------------------------------------------------------------==
    ! ==-  read input parameters for cntl%md run                          -==
    ! ==--------------------------------------------------------------==
    IF (cntl%tmdeh) THEN
       ! 
       IF (cntr%delt_elec.NE.0._real_8) THEN
          IF (td_prop%tintrvll.NE.cntr%delt_elec) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,F5.2)')&
                  ' WARNING: time step setted to ', cntr%delt_elec,' au'
             td_prop%tintrvll=cntr%delt_elec
             cntr%delt_ions=cntr%delt_elec
          ENDIF
       ELSE
          cntr%delt_elec=td_prop%tintrvll
          cntr%delt_ions=td_prop%tintrvll
       ENDIF
       ! 
       IF (cntl%cheby) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,2("="),62("-"))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,("number of zeros           :"),I5)')  td_prop%nzeros
          IF (paral%io_parent)&
               WRITE(6,'(1X,("epsilon for max           :"),F15.9)')td_prop%epsil
          IF (paral%io_parent)&
               WRITE(6,'(1X,("time intervall            :"),F15.9)')td_prop%tintrvll
          IF (paral%io_parent)&
               WRITE(6,'(1X,("number of cycles          :"),I5)') td_prop%n_cycles
          IF ((td_prop%read_wf.EQ.1) .AND.paral%io_parent)&
               WRITE(6,'(1X,("read excit. wfs           :"))')
          IF ((td_prop%read_wf.EQ.2).AND.paral%io_parent)&
               WRITE(6,'(1X,("read wfs from RESTART.1   :"))')
          IF (td_prop%stextpot) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation amplitude    :"),F15.9)')td_prop%ampl
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation direction    :"),I5)') td_prop%pertdir
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("type of perturbation      :"),I5)') td_prop%pert_type
             IF ((td_prop%pert_type.EQ.1).AND.paral%io_parent)&
                  WRITE(6,*) 'PULSE (DIRAC) FIELD'
             IF ((td_prop%pert_type.EQ.2).AND.paral%io_parent)&
                  WRITE(6,*) 'CONSTATN (HEAVISIDE) FIELD'
             IF (paral%io_parent)&
                  WRITE(6,'(1X,2("="),62("-"))')
          ENDIF
          IF (cntl%tgaugef) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation amplitude    :"),F15.9)')td_prop%ampl
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation direction    :"),I5)') td_prop%pertdir
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation frequency    :"),F15.9)') cntr%gfreq
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(1X,2("="),62("-"))')
       ELSEIF (cntl%cayley) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,2("="),62("-"))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,("number of subintervalls   :"),I5)')  td_prop%nzeros
          IF (paral%io_parent)&
               WRITE(6,'(1X,("epsilon for max           :"),F15.9)')td_prop%epsil
          IF (paral%io_parent)&
               WRITE(6,'(1X,("time intervall            :"),F15.9)')td_prop%tintrvll
          IF (paral%io_parent)&
               WRITE(6,'(1X,("number of cycles          :"),I5)') td_prop%n_cycles
          IF ((td_prop%read_wf.EQ.1) .AND.paral%io_parent)&
               WRITE(6,'(1X,("read excit. wfs           :"))')
          IF ((td_prop%read_wf.EQ.2).AND.paral%io_parent)&
               WRITE(6,'(1X,("read wfs from RESTART.1   :"))')
          IF (td_prop%stextpot) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation amplitude    :"),F15.9)')td_prop%ampl
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation direction    :"),I5)') td_prop%pertdir
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("type of perturbation      :"),I5)') td_prop%pert_type
             IF ((td_prop%pert_type.EQ.1).AND.paral%io_parent)&
                  WRITE(6,*) 'PULSE (DIRAC) FIELD'
             IF ((td_prop%pert_type.EQ.2).AND.paral%io_parent)&
                  WRITE(6,*) 'CONSTATN (HEAVISIDE) FIELD'
             IF (paral%io_parent)&
                  WRITE(6,'(1X,2("="),62("-"))')
             IF (td_prop%td_extpot) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) 'TIME DEP. EXTERNAL POTENTIAL FOR CAYLEY'
                IF (paral%io_parent)&
                     WRITE(6,*) 'not implemented'
                STOP
             ENDIF
          ENDIF
          IF (cntl%tgaugef) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation amplitude    :"),F15.9)')td_prop%ampl
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation direction    :"),I5)') td_prop%pertdir
             IF (paral%io_parent)&
                  WRITE(6,'(1X,("perturbation frequency    :"),F15.9)') cntr%gfreq
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(1X,2("="),62("-"))')
       ENDIF
       ! 
    ELSEIF (cntl%tpspec.OR.cntl%tpdist) THEN

       IF (paral%io_parent)&
            WRITE(6,'(1X,2("="),62("-"))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,("number of zeros           :"),I5)')  td_prop%nzeros
       IF (paral%io_parent)&
            WRITE(6,'(1X,("epsilon for max           :"),F15.9)')td_prop%epsil
       IF (paral%io_parent)&
            WRITE(6,'(1X,("time intervall            :"),F15.9)')td_prop%tintrvll
       IF (paral%io_parent)&
            WRITE(6,'(1X,("number of cycles          :"),I5)') td_prop%n_cycles
       IF ((td_prop%read_wf.EQ.1).AND.paral%io_parent)&
            WRITE(6,'(1X,("read excit. wfs           :"))')
       IF (td_prop%stextpot) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,("perturbation amplitude    :"),F15.9)')td_prop%ampl
          IF (paral%io_parent)&
               WRITE(6,'(1X,("perturbation direction    :"),I5)') td_prop%pertdir
          IF (paral%io_parent)&
               WRITE(6,'(1X,("type of perturbation      :"),I5)') td_prop%pert_type
          IF ((td_prop%pert_type.EQ.1).AND.paral%io_parent)&
               WRITE(6,*) 'PULSE (DIRAC) FIELD'
          IF ((td_prop%pert_type.EQ.2).AND.paral%io_parent)&
               WRITE(6,*) 'CONSTATN (HEAVISIDE) FIELD'
          IF (paral%io_parent)&
               WRITE(6,'(1X,2("="),62("-"))')
       ENDIF

    ENDIF

    IF (ppcorr) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') " GAUGEING OF THE NON-LOCAL PSEUDOPOTENTIALS"
    ENDIF

9999 CONTINUE
    ! 
    ! SOME TESTS AND SETTINGS
    ! -----------------------
    IF (td_prop%read_wf.NE.0) cntl%start_real=.TRUE.
    IF (paral%io_parent)&
         WRITE(6,'(a)')&
         ' Start from real wavefunctions in RESTART.1'
    IF ((td_prop%stextpot.AND.(td_prop%pert_type.GT.1)) .OR.&
         td_prop%td_extpot .OR.&
         td_prop%tpointch) THEN
       ! external field add in subroutine 'vofrhob' to the total potential
       ! but the interaction energy with the field is NOT add to the 
       ! total energy.
       textfld=.TRUE.
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(a)') ' External field is computed and added'
       ENDIF
    ENDIF

    IF (td_prop%ionize_from_rdm.OR.masklog%tmask) THEN
       td_prop%do_not_normalize=.TRUE.
    ENDIF

    IF (cntl%tmdeh.AND.bohm%tplanewave_wp.AND.bohm%tbohmtraj) THEN
       IF ( (NINT(SQRT(1.0_real_8*bohm%nf_bohmtraj))**2).NE.bohm%nf_bohmtraj) THEN
          IF(paral%io_parent) THEN
             WRITE(6,'(a)') 'number of bohm trj must be a perfect square'
             STOP
          ENDIF
       ENDIF
    ENDIF
    IF ((cntl%tgaugep.OR.cntl%tgaugef).AND.(td_prop%pertdir.GT.3)) THEN
       td_prop%pertdir=3
    ENDIF
    IF (cntl%tpdist.AND.(.NOT.td_prop%tpointch)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a)') ' PROP DISTURBANCE WITH POINTCHARGE PERTURBATION'
       STOP
    ELSE
       ! cntl%ruku=.false.
    ENDIF
    IF (cntl%tmdeh.AND.cntl%tpdist.AND.(.NOT.cntl%ruku)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a)') ' EHMD with POINT PERTURBATION NOT TESTED'
       STOP
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(a,a,/)')&
         ' **************************************',&
         '************************** '
    !
    ! BROADCAST INPUT PARAMETERS
    ! --------------------------
    CALL mp_bcast_byte(td_prop,size_in_bytes_of(td_prop),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(pointcharge,size_in_bytes_of(pointcharge),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(gaugefield,size_in_bytes_of(gaugefield),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(tcurr,size_in_bytes_of(tcurr),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(curr_step,size_in_bytes_of(curr_step),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(ppcorr,size_in_bytes_of(ppcorr),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(tberry,size_in_bytes_of(tberry),parai%io_source,parai%cp_grp)

    CALL mp_bcast(cntr%delt_ions,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntr%delt_elec,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntl%start_real,parai%io_source,parai%cp_grp)
    CALL mp_bcast(td_prop%ionize_from_state,parai%io_source,parai%cp_grp)
    CALL mp_bcast(td_prop%ionize_from_rdm,parai%io_source,parai%cp_grp)
    CALL mp_bcast(td_prop%do_not_normalize,parai%io_source,parai%cp_grp)
    CALL mp_bcast(td_prop%td_fix_spin_dens,parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(masklog,size_in_bytes_of(masklog),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(maskreal,size_in_bytes_of(maskreal),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(bohm,size_in_bytes_of(bohm),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE read_prop
  ! ==================================================================

END MODULE read_prop_utils
