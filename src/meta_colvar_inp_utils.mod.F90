MODULE meta_colvar_inp_utils
  USE chain_dr_utils,                  ONLY: chain_dr
  USE cnst,                            ONLY: au_kcm
  USE cnst_dyn,                        ONLY: &
       atcvar, bbeta, cscl_fac, cv_dtemp, cv_dyn_0, cv_ist, cv_langamma, &
       cv_langevintemp, cv_mass, cv_temp, cv_temp0, cvpar, det_celvar, &
       det_colvar, dof_incv, fhills, file_str_ab, fmtdres, hhm, hlhm, hthm, &
       hvm0, hwm, iangcv, iatcnga, iatcngb, iatdlmn, ibound, icv_cell, &
       icv_rmsd_ab, icv_spin, imeta, initial_value, inter_hill, &
       inter_hill_max, iqatdlmn, kharm, lcvtc, lkfix, lmeta, ltcglobal, lupd, &
       max_natcngb, max_nnvar, maxrmsdat, mdcellr, natcngb, ncolvar, ncvsys, &
       nlong, nrmsd_ab, nshort, nsubsys, optdist, rcc, rcc0, rccnga, rch, &
       rcw, rmeta, specindex, tad_scf, tcvlangevin, tcvscale, toll_avcv, &
       trmsd_ab, tvolbound, tycvar, vbound
  USE cnstfc_utils,                    ONLY: bndswitch,&
                                             coorn_rf,&
                                             coornum,&
                                             coornumgrp,&
                                             coornumsp
  USE coninp_utils,                    ONLY: raddeg
  USE constr_utils,                    ONLY: &
       diffd, diffdd, diffo, diffr, difft, funcd, funcdd, funco, funcr, funct
  USE cotr,                            ONLY: cotc0,&
                                             dtm,&
                                             duat,&
                                             lskptr
  USE dum2_utils,                      ONLY: dum2
  USE error_handling,                  ONLY: stopgm
  USE fillc_utils,                     ONLY: fillc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE meta_cv_qmmm_utils,              ONLY: coorn_rf_q,&
                                             coornseq_rf,&
                                             coorntot_rf_q,&
                                             rmsd_seq_a
  USE meta_cv_utils,                   ONLY: &
       abs_difcntot_rf, angle_cv, aplane, cdipolephi, cdipolerho, &
       cdipoletheta, coorntot_rf, ctot_chain, cut_hyd_ion_distance, &
       difcntot_rf, disp_lnm, funcl, funcp_mia, hyd_cutoff, hyd_ion_distance, &
       hyd_presence, rmsd_ab, side_cv, vcoors, volume_cv
  USE meta_ex_mul_util_utils,          ONLY: rmtdresm
  USE meta_exlagr_utils,               ONLY: rmtdres
  USE meta_multiple_walkers_utils,     ONLY: mw_filename,&
                                             rmtdres_mw
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert,&
                                             mmdim,&
                                             nat_cpmd,&
                                             nat_grm
  USE mm_input,                        ONLY: lqmmm
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE odiis_utils,                     ONLY: solve
  USE parac,                           ONLY: paral
  USE prcp,                            ONLY: prcp_com
  USE puttau_utils,                    ONLY: gettau
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: meta_colvar_inp
  !public :: settological
  PUBLIC :: colvar_structure
  !public :: get_dof
  PUBLIC :: colvarofr
  PUBLIC :: colvarpr
  !public :: param_meta

  PUBLIC :: cv_forces3


CONTAINS

#ifdef __SR11000
  !option MP(P(0)), LANGLVL(SAVE(0))
#endif
  ! ==================================================================
  SUBROUTINE meta_colvar_inp(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Collective Variables Input for Metadynamics           ==
    ! ==--------------------------------------------------------------==
    ! 
    ! NCOLVAR   :  number of Collective Variables
    ! TYCVAR    :  type of Collective Variables
    ! ATCVAR    :  atoms which define the Collective Variables
    ! CVPAR     :  other parameters that are required for the C.V.
    ! I_META_MAX:  max # of steps in the MetaDynamics
    ! I_META_RES: # of meta-steps performed in a previous meta-run
    ! used in case of RESTART from previous run 
    ! INTER_HILL: # of ions dynamic steps between two meta steps
    ! INTER_AV_CV: # of IONS configurations requied for the calculation 
    ! of the collective variables as average  
    ! HLLH      : hills altitude 
    ! HLLW      : hills amplitude 
    ! CSCL_FAC  : scale factors required for an homogeneous dynamic in the
    ! NCOLVAR directions in the space of C.V.
    ! they can be read from input or determined from the
    ! second derivatives of the potential well
    ! a tuning of this factors along the dynamics is also
    ! to be considered
    ! IBOUND    : 1 if a wall-potential constrains the C.V. 
    ! within a certain region of space 
    ! otherwise 0 (integer, dimension NCOLVAR)
    ! VBOUND    : parameters for the wall-potential definition  see sub. SETWALL
    ! LAGRANGIAN  The dynamics of the Collective Variable is driven by 
    ! an extended Lagrangean method
    ! With this option several parameter are to be defined
    ! CV_TEMP0:  Initial Temperature of Collective Variables 
    ! Maxwell-Boltzman velocity distribution
    ! LCVTC   :  Temperature Control flag
    ! parameters (Velocity Rescaling) CV_TEMP CV_DTEMP
    ! CV_MASS :  Masses of collective Varibles
    ! KHARM   :  Harmonic constants for the Harmonic potential Vharm 
    ! ==--------------------------------------------------------------==
    ! how C.V. definition is read
    ! 
    ! DEFINE VARIABLE
    ! ==   nfix                                                       ==
    ! ==  DIST      IA1  IA2             "SCA/SCF" sc min max  "KCV" k "MCV" m  ==
    ! ==  STRETCH   IA1  IA2             "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  BEND      IA1  IA2  IA3        "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  TORSION   IA1  IA2  IA3  IA4   "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  OUTP      IA1  IA2  IA3  IA4   "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  COORD     IA1   K  RC          "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  DIFFER    IA1  IA2  IA3        "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  COORSP    IA1  ns1  K RC       "SCA/SCF" sc min max  "KCV" k "MCV" m  ==    
    ! ==  COOR_RF   IA1  ns1  exp1  exp2  RC "SCA/SCF" sc min max "KCV" k "MCV" m  ==   
    ! ==  COORGROUP NA   K               "SCA/SCF" sc min max  "KCV" k "MCV" m  ==
    ! ==            in the following line:                                      ==
    ! ==            IA1  NB1                                                    ==
    ! ==                   IB1...IB_NB1  RC_IA1                                 ==
    ! ==  COOR_RF  "INDAT" IA1  numatoms  exp1  exp2  RC                        ==
    ! "SCA/SCF" sc min max  "KCV" k "MCV" m  == 
    ! in the following line write the list of atoms
    ! ==  BNSWT     IA1  IA2  exp1  exp2  RC "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! ==  TOT_COOR  ns1  ns2  exp1  exp2  RC "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! ==  RMSD_AB   numsp ns1 . ns_numsp     "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! ==  DISPL1    numspA numspB spA_1..spA_numspA spB_1..spB_numspB           ==
    ! "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! ==  DISPL1  "INDAT" numatA numatB spA_1  list atA   list atB                ==
    ! "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! if the lists are too long, write in line "LONGLIST" and
    ! in following lines: indexes A (groups of 10 per line) (last line what remains)
    ! in following lines: indexes B (groups of 10 per line)       ''
    ! ==  PLNANG    IA1 IA2 IA3  IB1 IB2 IB3  "SCA/SCF" sc min max  "KCV" k "MCV" m  ==
    ! ==  HBONDCH   NRES1  NATTR NRES2 NH  exp1_H exp2_H r0_h exp1_A esp2_A r0_A     ==
    ! "SCA/SCF" sc min max "KCV" k  "MCV" m ==
    ! in following line : indexes RES1
    ! in following lines: indexes hbond attractors (groups of 10 per line)
    ! in following line : indexes RES2
    ! in following lines: indexes h attractors (groups of 10 per line)
    ! HBONDPAR  "RCW" rcw "RCH" rch  "LENGTH"  nshort nlong "PAR" bbetha  const 'LUPD' lupd
    ! ==  DIFCOOR   IA1 IA2 ns1 exp1  exp2  RC   "SCA/SCF" sc min max "KCV" k "MCV" m ==
    ! ==  DIFCOOR "INDAT" IA1 IA2 numatoms  exp1  exp2  RC                            ==
    ! "SCA/SCF" sc min max  "KCV" k "MCV" m ==
    ! ==  COOR_CHAIN ISP1 ISP2 ISP3 n m rc1 rc2 "SCA/SCF" sc min max  "KCV" k "MCV" m ==
    ! ==  HYDRONIUM ISP1 ISP2     "SCA/SCF" sc min max "KCV" k "MCV" m    "="         ==
    ! in following line : "H:" n m rc "O:" n m rc "FCUT:" n m rc "LAMBDA" L
    ! ==  HYD_CUT ISP1 ISP2 NUML  "SCA/SCF" sc min max "KCV" k "MCV" m    "="         ==
    ! in following line : "H:" nH mH rcH "O:" nO mO rcO "M:" nM mM rcM1 rcM2 "FCUT:" nL mL rcL "LAMBDA" L
    ! in following line : list of numL atomic indexes
    ! ==  DIS_HYD ISP1 ISP2 IAT   "SCA/SCF" sc min max "KCV" k "MCV" m    "="         ==
    ! in following line : "H:" n m rc "O:" n m rc "FCUT:" n m rc "LAMBDA" L
    ! ==  SPIN  IAT                            "SCA/SCF" sc min max "KCV" k "MCV" m   ==
    ! ==  VOLVAR                               "SCA/SCF" sc min max "KCV" k "MCV" m   ==
    ! ==  CELLSIDE 1/2/3                       "SCA/SCF" sc min max "KCV" k "MCV" m   ==
    ! ==  CELLANGLE 1/2/3                      "SCA/SCF" sc min max "KCV" k "MCV" m   ==
    ! ==  ...  
    ! END DEFINITION


    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_colvar_inp'
    CHARACTER(len=15), DIMENSION(36), PARAMETER :: styp = (/'        STRETCH',&
      '           BEND','        TORSION','       DISTANCE','           OUTP',&
      '        COORDIN','      DIFFEREN.','         COORSP','        COOR_RF',&
      '      BN.SWITCH','       TOT C_RF','        RMSD_AB','         DISPL1',&
      '        COORDIS','         PLNANG','        HBONDCH','        DIFCOOR',&
      '      COO_CHAIN','      HYDRONIUM','        DIS_HYD','           SPIN',&
      '         VOLVAR','       CELLSIDE','      CELLANGLE','       COOR_SEQ',&
      '       RMSD_SEQ','        HYD_CUT','      D_HYD_CUT','      COORGROUP',&
      '        DISAXIS','      DIPOLERHO','      DIPOLETHA','      DIPOLEPHI',&
      '         VCOORS','  DIFF TOT C_RF','|DIFF TOT C_RF|'/)

    CHARACTER(len=10)                        :: chnum, chnum2
    CHARACTER(len=100)                       :: lineform
    CHARACTER(len=80)                        :: line, line2
    INTEGER :: i, i1, i2, ia, iaa, iatma, iatmb, ic, idummy, ie, iee, iend, &
      ierr, ii, ij, ind, inum, iout, ipgrpa, ipgrpa1, ipgrpb, ipgrpb1, &
      ipgrpb2, ipgrpb3, ityp, iwalk1, k, kk, my_nat, nloop1, nloop2, nnvar, &
      numspec, numspec2, numtot
    INTEGER, SAVE                            :: firstcorgrp = 0, &
                                                firstdlmn = 0, secondtdlmn = 0
    LOGICAL                                  :: erread, status
    REAL(real_8)                             :: cpos_0, fmax, fmin, r0_shift, &
                                                r_wall, temp, vbar
    REAL(real_8), ALLOCATABLE                :: dummy(:)

    lmeta%lcolvardyn = .TRUE.
    ! get total number of atoms and set indexing.
    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
    ELSE
       my_nat=mmdim%natm
    ENDIF
    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=20,END=20,fmt='(A)') line

    IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'METADYN').NE.0) THEN
       GOTO 30
    ENDIF

    IF (INDEX(line,'DEF').NE.0.AND. INDEX(line,'VARIABLE').NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) nnvar

       ! Allocate Memory
       ALLOCATE(tycvar(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(tycvar)!,nnvar)
       ALLOCATE(atcvar(15,nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(atcvar)!,15*nnvar)
       ALLOCATE(iangcv(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(iangcv)!,nnvar)
       ALLOCATE(cvpar(10,nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvpar)!,10*nnvar)
       ALLOCATE(ibound(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ibound)!,nnvar)
       ALLOCATE(vbound(4,nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(vbound)!,4*nnvar)
       ALLOCATE(cscl_fac(3,nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tad_scf(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cv_mass(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_mass)!,nnvar)
       ALLOCATE(kharm(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kharm)!,nnvar)
       ALLOCATE(specindex(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(specindex)!,nnvar)
       ALLOCATE(cv_dyn_0(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_dyn_0)!,nnvar)
       ALLOCATE(initial_value(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL settological(initial_value,nnvar,.FALSE.)

       IF (lmeta%tmulti) THEN
          ALLOCATE(hwm(nsubsys),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(hhm(nsubsys),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(hlhm(nsubsys),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(hthm(nsubsys),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(hvm0(nsubsys),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i = 1,nsubsys
             hwm(i) = rmeta%hllw
             hhm(i) = rmeta%hllh
             hlhm(i) = rmeta%hlow
             hthm(i) = rmeta%htop
          ENDDO
       ENDIF

       ! Set Default Values
#ifdef __SR11000
       !poption parallel, tlocal(IC)
       !voption indep(CSCL_FAC,TAD_SCF)
#endif
       DO ic = 1,nnvar
          cscl_fac(1,ic) = 1.0_real_8
          cscl_fac(2,ic) = 0.2_real_8
          cscl_fac(3,ic) = 2._real_8
          tad_scf(ic)  = .FALSE.
       ENDDO
       GOTO 12
       ! Read Type and Atoms for the Collective Variables to be set 
11     READ(iunit,err=20,END=20,fmt='(A)') line
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'DEF').NE.0)&
            GOTO 10
12     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       ! write(6,*) LINE
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'DEF').NE.0) THEN
          ! Problem the number of coll. variables NCOLVAR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'M_COLVAR_INP! # OF GIVEN C.V. =',NNVAR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'M_COLVAR_INP! # OF FOUND C.V.=',NCOLVAR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')&
               'M_COLVAR_INP! END C.V. IS REACHED.'
          CALL stopgm('M_COLVAR_INP',&
               'ERROR WHILE READING FIX STRUCTURES',& 
               __LINE__,__FILE__)
          ! ==------------------- DISTANCE --------------------------------==
       ELSEIF (INDEX(line,'DIST').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 4
          i1=INDEX(line,'DIST')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
       ELSEIF (INDEX(line,'DISAXIS').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 30
          i1=INDEX(line,'DISAXIS')+7
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=idummy
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--------------- DIFFERENCE AMONG DISTANCES ------------------==
       ELSEIF (INDEX(line,'DIFFER').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 7
          i1=INDEX(line,'DIFFER')+6
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------------ STRETCH ----------------------------------==
       ELSEIF (INDEX(line,'STRETCH').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 1
          i1=INDEX(line,'STRETCH')+7
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------------  BEND ------------------------------------==
       ELSEIF (INDEX(line,'BEND').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 2
          i1=INDEX(line,'BEND')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------------ TORSION ANGLE  ----------------------------==
       ELSEIF (INDEX(line,'TORSION').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 3
          iangcv(ncolvar) = 1
          i1=INDEX(line,'TORSION')+7
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------------ OUT OF LPANE ANGLE ------------------------==
       ELSEIF (INDEX(line,'OUTP').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 5
          iangcv(ncolvar) = 1
          i1=INDEX(line,'OUTP')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------------ COORDINATION NUMBER ------------------------==
       ELSEIF (INDEX(line,'COORD').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 6
          i1=INDEX(line,'COORD')+5
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER  -----------==
          ! ==                 f = 1/(1+exp(k*(R-R_0)))                      ==
       ELSEIF (INDEX(line,'COORSP').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 8
          i1=INDEX(line,'COORSP')+6
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          ! mb For QMMM use the qmmm CV
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--- SPECIES DEPENDENT COORDINATION NUMBER OF POINT IN SPACE ---=
          ! ==                 f = 1/(1+exp(k*(R-R_0)))                      =
       ELSEIF (INDEX(line,'VCOORS').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 34
          i1=INDEX(line,'VCOORS')+6
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
          IF (paral%io_parent)&
               READ(iunit,*) cvpar(3,ncolvar),cvpar(4,ncolvar),&
               cvpar(5,ncolvar)


          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==-------------- COORDINATION NUMBER FOR SPECIES IN GROUP-----------------==
          ! ==                 f = sum_i^NA sum_j^NB(i) 1/(1+exp(k*(R_ij-R_0(i)) ) )  ==

       ELSEIF (INDEX(line,'COORGROUP').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 29
          i1=INDEX(line,'COORGROUP')+9
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)! number of a atoms
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)! parameter k
          i1=iout
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          IF (firstcorgrp.EQ.0) THEN
             ALLOCATE(iatcnga(atcvar(1,ncolvar)*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)     ! stores index of atom A
             ALLOCATE(iatcngb(atcvar(1,ncolvar)*max_natcngb*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)     ! stores index of atom B
             ALLOCATE(natcngb(atcvar(1,ncolvar)*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)     ! stores number of atoms of B for each A
             ALLOCATE(rccnga(atcvar(1,ncolvar)*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)      ! stores the cut-off distances for all A
             CALL zeroing(iatcnga)!,atcvar(1,ncolvar)*nnvar)
             CALL zeroing(iatcngb)!,atcvar(1,ncolvar)*max_natcngb*nnvar)
             CALL zeroing(rccnga)!,atcvar(1,ncolvar)*nnvar)
             firstcorgrp = 1
             ipgrpa=0
             ipgrpb=0
             ipgrpb2=0
          ENDIF
          DO iatma=1,atcvar(1,ncolvar)
             ipgrpa=ipgrpa+1
             IF (paral%io_parent)&
                  READ(iunit,*) iatcnga(ipgrpa), rccnga(ipgrpa),&
                  natcngb(ipgrpa)
             iatcnga(ipgrpa)=NAT_cpmd(iatcnga(ipgrpa))! A's in CPMD order
             IF (natcngb(ipgrpa).GT.max_natcngb)THEN
                IF (paral%io_parent)&
                     WRITE(6,'(/,A,I8)')&
                     'Maximum number of type B atoms for COORGROUP =',&
                     max_natcngb
                IF (paral%io_parent)&
                     WRITE(6,'(A,I8,A,I8)')&
                     'For coll. variable ',ncolvar,' number of B atoms =',&
                     natcngb(ipgrpa)
                CALL stopgm('M_COLVAR_INP',&
                     'MAX_NATCNGB < NATCNGB(IPGRPA)',& 
                     __LINE__,__FILE__)
             ENDIF
             IF (paral%io_parent)&
                  READ(iunit,*)&
                  (iatcngb(ipgrpb+iatmb),iatmb=1,natcngb(ipgrpa)) ! read in atoms B for a give A atom
             DO iatmb=1,natcngb(ipgrpa)
                iatcngb(ipgrpb+iatmb)=NAT_cpmd(iatcngb(ipgrpb+iatmb))! B's in CPMD order
             ENDDO
             ipgrpb=ipgrpb+natcngb(ipgrpa)
          ENDDO
          ipgrpb2=ipgrpb2+max_natcngb*atcvar(1,ncolvar)
          ipgrpb=ipgrpb2
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER -------------==
          ! ==            f =sum_i (1-(R/R_0)^n)/(1-(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'COOR_RF').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 9
          i1=INDEX(line,'COOR_RF')+7
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') +5
          ELSEIF(INDEX(line,'ELEMENT') .NE. 0 .AND.&
               INDEX(line,'SOLUTE').NE.0) THEN
             specindex(ncolvar) = -1
             i1 = INDEX(line,'SOLUTE') +6
          ELSEIF(INDEX(line,'ELEMENT') .NE. 0 .AND.&
               INDEX(line,'SOLVENT').NE.0) THEN
             specindex(ncolvar) = -2
             i1 = INDEX(line,'SOLVENT') +7
          ELSEIF (INDEX(line,'ELEMENT') .NE. 0 ) THEN
             specindex(ncolvar) = -3
             i1 = INDEX(line,'ELEMENT') +7
          ELSEIF(INDEX(line,'SEQUENCE') .NE. 0 .AND.&
               INDEX(line,'mM').NE.0) THEN
             specindex(ncolvar) = -4
             i1 = INDEX(line,'mM') +2
          ELSEIF (INDEX(line,'SEQUENCE') .NE. 0) THEN
             specindex(ncolvar) = -5
             i1 = INDEX(line,'SEQUENCE') +8
          ELSEIF(INDEX(line,'SEQ_ELE') .NE. 0 .AND.&
               INDEX(line,'mM').NE.0) THEN
             specindex(ncolvar) = -6
             i1 = INDEX(line,'mM') +2
          ELSEIF (INDEX(line,'SEQ_ELE') .NE. 0) THEN
             specindex(ncolvar) = -7
             i1 = INDEX(line,'SEQ_ELE') +7
          ELSE
             specindex(ncolvar) = 1
          ENDIF

          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF
             IF (lqmmm%qmmm)  CALL stopgm('M_COLVAR_INP',&
                  'INDAT not implemented for QMMM',& 
                  __LINE__,__FILE__)
             IF (paral%io_parent)&
                  READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
                  ii=1,atcvar(2,ncolvar))

          ELSEIF (specindex(ncolvar) .LT. -3) THEN
             IF (paral%io_parent)&
                  READ(iunit, *) atcvar(5,ncolvar), atcvar(6,ncolvar)
          ENDIF

          ii=INDEX(line,'2SHELL')
          IF (ii.NE.0) THEN
             i1=ii+6! missing line
             CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
             ia=atcvar(3,ncolvar)
             ie=atcvar(4,ncolvar)
             IF (MOD(ia,2).NE.0 .OR. MOD(ie,2).NE.0 ) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'for COOR_RF with 2nd shell use even exponents',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==--------------- BOND SWITCH WITH RATIONAL F -------------------==
          ! ==            f = (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'BNSWT').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 10
          i1=INDEX(line,'BNSWT')+5
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--------- TOTAL SPECIES DEPENDENT COORDINATION NUMBER ---------==
          ! ==      f =sum_j[sum_i (1+(Rij/R_0)^n)/(1+(Rij/R_0)^(n+m))]/NA   ==
       ELSEIF (INDEX(line,'TOT_COOR').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 11
          i1=INDEX(line,'TOT_COOR') +8
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') +5
          ELSEIF(INDEX(line,'ELEMENT') .NE. 0 .AND.&
               INDEX(line,'SOLUTE').NE.0) THEN
             specindex(ncolvar) = -1
             i1 = INDEX(line,'SOLUTE') +6
          ELSEIF(INDEX(line,'ELEMENT') .NE. 0 .AND.&
               INDEX(line,'SOLVENT').NE.0) THEN
             specindex(ncolvar) = -2
             i1 = INDEX(line,'SOLVENT') +7
          ELSEIF (INDEX(line,'ELEMENT') .NE. 0 ) THEN
             specindex(ncolvar) = -3
             i1 = INDEX(line,'ELEMENT') +7
          ELSEIF(INDEX(line,'SEQUENCE') .NE. 0 .AND.&
               INDEX(line,'mM').NE.0) THEN
             specindex(ncolvar) = -4
             i1 = INDEX(line,'mM') +2
          ELSEIF (INDEX(line,'SEQUENCE') .NE. 0) THEN
             specindex(ncolvar) = -5
             i1 = INDEX(line,'SEQUENCE') +8
          ELSEIF(INDEX(line,'SEQ_ELE') .NE. 0 .AND.&
               INDEX(line,'mM').NE.0) THEN
             specindex(ncolvar) = -6
             i1 = INDEX(line,'mM') +2
          ELSEIF (INDEX(line,'SEQ_ELE') .NE. 0) THEN
             specindex(ncolvar) = -7
             i1 = INDEX(line,'SEQ_ELE') +7

          ELSE
             specindex(ncolvar) = 1
          ENDIF
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF

             IF (lqmmm%qmmm)  CALL stopgm('M_COLVAR_INP',&
                  'INDAT not implemented for QMMM',& 
                  __LINE__,__FILE__)

             IF (paral%io_parent)&
                  READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
                  ii=1,atcvar(2,ncolvar))

          ELSEIF (specindex(ncolvar) .LT. -3) THEN
             IF (paral%io_parent)&
                  READ(iunit, *) atcvar(5,ncolvar), atcvar(6,ncolvar)
          ENDIF

          ii=INDEX(line,'2SHELL')
          IF (ii.NE.0) THEN
             i1=ii+6! missing line
             CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
             ia=atcvar(3,ncolvar)
             ie=atcvar(4,ncolvar)
             IF (MOD(ia,2).NE.0 .OR. MOD(ie,2).NE.0 ) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'for TOT_COOR with 2nd shell use even exponents',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--------- TOTAL SPECIES DEPENDENT COORDINATION NUMBER OF DIFFERENCE ---------==
          ! ==   f =sum_j[sum_i (1+(Rij/R_0)^n)/(1+(Rij/R_0)^(n+m))]/NA-sum_j[sum_k[]]/NA ==
       ELSEIF (INDEX(line,'TOT_DIFCOOR').NE.0) THEN
          ncolvar = ncolvar + 1
          ! Tag the absolute value of the difference 
          IF (INDEX(line,'|TOT_DIFCOOR|').NE.0) THEN
             tycvar(ncolvar) = 36
             i1=INDEX(line,'|TOT_DIFCOOR|') + 13
          ELSE
             tycvar(ncolvar) = 35
             i1=INDEX(line,'TOT_DIFCOOR') + 11
          ENDIF
          specindex(ncolvar) = 1
          ! three species
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          ! n, m for rational function.
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)

          ! scaling factor, k, m, low and upper bounds, turning on lagrange formulations     
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--------- COORDINATION NUMBER BETWEEN SEQUENCES (QMMM) ---------==
          ! ==      f =sum_j[sum_i (1+(Rij/R_0)^n)/(1+(Rij/R_0)^(n+m))]/NA   ==
       ELSEIF (INDEX(line,'COOR_SEQ').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 25
          i1=INDEX(line,'COOR_SEQ') +8

          IF (INDEX(line,'SOLUTE_1') .NE. 0 .AND.&
               INDEX(line,'SOLVENT').NE.0) THEN
             specindex(ncolvar) = 1
             i1 = INDEX(line,'SOLVENT') +7
          ELSEIF(INDEX(line,'SOLUTE_1') .NE. 0 .AND.&
               INDEX(line,'SOLUTE_2').NE.0) THEN
             specindex(ncolvar) = 2
             i1 = INDEX(line,'SOLUTE_2') +8
          ELSEIF(INDEX(line,'SEQ_1') .NE. 0 .AND.&
               INDEX(line,'SOLVENT') .NE. 0) THEN
             specindex(ncolvar) = 3
             i1 = INDEX(line,'SOLVENT') +7
          ELSEIF(INDEX(line,'SEQ_1') .NE. 0 .AND.&
               INDEX(line,'SEQ_2') .NE. 0 ) THEN
             specindex(ncolvar) = 4
             i1 = INDEX(line,'SEQ_2') +5
          ELSEIF(INDEX(line,'SEQ_1') .NE. 0 .AND.&
               INDEX(line,'SEQ_MM2') .NE. 0) THEN
             specindex(ncolvar) = 5
             i1 = INDEX(line,'SEQ_MM2') +7
          ELSEIF(INDEX(line,'SEQ_QM1') .NE. 0 .AND.&
               INDEX(line,'SEQ_MM2') .NE. 0) THEN
             specindex(ncolvar) = 6
             i1 = INDEX(line,'SEQ_MM2') +7
          ELSEIF(INDEX(line,'SEQ_MM1') .NE. 0 .AND.&
               INDEX(line,'SEQ_MM2') .NE. 0 ) THEN
             specindex(ncolvar) = 7
             i1 = INDEX(line,'SEQ_MM2') +7
          ELSEIF(INDEX(line,'SEQ_QM1') .NE. 0 .AND.&
               INDEX(line,'SEQ_2') .NE. 0 ) THEN
             specindex(ncolvar) = 8
             i1 = INDEX(line,'SEQ_2') +5
          ENDIF

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)

          IF (specindex(ncolvar) .EQ. 3) THEN
             IF (paral%io_parent)&
                  READ(iunit, *) atcvar(5,ncolvar), atcvar(6,ncolvar)
          ELSEIF (specindex(ncolvar) .GT. 3) THEN
             IF (paral%io_parent)&
                  READ(iunit, *) atcvar(5,ncolvar), atcvar(6,ncolvar),&
                  atcvar(7,ncolvar), atcvar(8,ncolvar)
          ENDIF

          ii=INDEX(line,'2SHELL')
          IF (ii.NE.0) THEN
             i1=ii+6! missing line
             CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
             ia=atcvar(3,ncolvar)
             ie=atcvar(4,ncolvar)
             IF (MOD(ia,2).NE.0 .OR. MOD(ie,2).NE.0 ) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'for COOR_SEQ with 2nd shell use even exponents',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==------- RMSD wrt 2 configurations A and B (read from file)-----==
          ! ==                f = (RMSD_A-RMSD_B)/(RMSD_A+RMSD_B)            ==
       ELSEIF (INDEX(line,'RMSD_AB').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 12
          IF (.NOT. trmsd_ab) THEN
             IF (max_nnvar.LT.nnvar) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) "NNVAR:",nnvar,&
                     "  MAX_NNVAR:",max_nnvar
                CALL stopgm("FILE_STR_AB",&
                     "PARAMETER MAX_NNVAR too small",& 
                     __LINE__,__FILE__)
             ENDIF

             ALLOCATE(icv_rmsd_ab(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_rmsd_ab)!,nnvar)
             trmsd_ab=.TRUE.
             nrmsd_ab = 0
          ENDIF
          icv_rmsd_ab(ncolvar)=ncolvar
          nrmsd_ab  = nrmsd_ab + 1

          i1=INDEX(line,'RMSD_AB') +7
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          IF (atcvar(1,ncolvar) .GT. 13) CALL stopgm('M_COLVAR_INP',&
               'RMSD_AB! too many species (>13)',& 
               __LINE__,__FILE__)
          numspec = atcvar(1,ncolvar)
          DO i = 1,numspec
             CALL readsi(line,i1,iout,atcvar(i+1,ncolvar),erread)
             i1=iout
          ENDDO
          ii=INDEX(line,'FILEAB')
          IF (ii.NE.0) THEN
             ia = ii+6
             line2  = line(ia:ia+20)
             CALL xstring(line2,ia,ie)
             file_str_ab(ncolvar) = line2(ia:ie)
          ELSE
             file_str_ab(ncolvar) = 'STRUCTURE_AB'
          ENDIF
          atcvar(numspec+2,ncolvar) = nrmsd_ab

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

       ELSEIF (INDEX(line,'RMSD_SEQ').NE. 0) THEN
          ! RMSD for QMMM, only absolute value wrt initial state for the moment
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 26
          i1=INDEX(line,'RMSD_SEQ') +8

          IF (INDEX(line,'SOLUTE') .NE. 0 ) THEN
             specindex(ncolvar) = 1
             i1 = INDEX(line,'SOLUTE') +6
          ELSEIF (INDEX(line,'SOLVENT') .NE. 0) THEN
             specindex(ncolvar) = 2
             i1 = INDEX(line,'SOLVENT') +7
          ELSEIF (INDEX(line,'SEQUENCE') .NE. 0) THEN
             specindex(ncolvar) = 3
             i1 = INDEX(line,'SEQUENCE') +8
          ENDIF

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          ! Determine max number of atoms to allocate the least space
          IF (atcvar(2,ncolvar) .GT. maxrmsdat) THEN
             maxrmsdat =    atcvar(2,ncolvar)
          ENDIF

          IF (atcvar(1,ncolvar) .GT. 10) THEN
             CALL stopgm('M_COLVAR_INP',&
                  'RMSD_SEQ! too many species (>10)',& 
                  __LINE__,__FILE__)
          ENDIF
          numspec = atcvar(1,ncolvar)
          DO i = 1,numspec
             CALL readsi(line,i1,iout,atcvar(i+4,ncolvar),erread)
             i1=iout
          ENDDO


          IF (specindex(ncolvar) .EQ. 3) THEN
             IF (paral%io_parent)&
                  READ(iunit, *) atcvar(3,ncolvar), atcvar(4,ncolvar)
          ENDIF
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! 
          ! ==--- (lmn) Displacement of one species wrt configuration A -----==
          ! FIXME: AK 2005/05/27 convert to use Gromos ordering of atom indices
       ELSEIF (INDEX(line,'DISPL1').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 13
          i1=INDEX(line,'DISPL1') +6
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') + 5
          ELSE
             specindex(ncolvar) = 1
          ENDIF
          IF (specindex(ncolvar) .EQ. 1) THEN
             CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
             i1=iout
             numspec = atcvar(1,ncolvar)
             CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
             i1=iout
             numspec2 = atcvar(2,ncolvar)
             IF (numspec+numspec2 .GT. 10) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'DISPL! too many species (>10)',& 
                     __LINE__,__FILE__)
             ENDIF
             DO i = 1,numspec
                CALL readsi(line,i1,iout,atcvar(i+2,ncolvar),erread)
                i1=iout
             ENDDO
             IF (numspec2 .NE. 0) THEN
                i2=numspec+2
                DO i = 1,numspec2
                   CALL readsi(line,i1,iout,atcvar(i+i2,ncolvar),erread)
                   i1=iout
                ENDDO
             ENDIF
             i2=numspec+numspec2+2
          ELSE
             CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
             i1=iout
             IF (atcvar(1,ncolvar) .GT. my_nat/2) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'DISPL! too many atoms  (>NAT/2)',& 
                     __LINE__,__FILE__)
             ENDIF
             numspec = atcvar(1,ncolvar)
             CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
             i1=iout
             numspec2 = atcvar(2,ncolvar)
             IF (numspec+numspec2 .GT. my_nat) THEN
                CALL stopgm('M_COLVAR_INP','DISPL1: NAT1+NAT2 > NAT',& 
                     __LINE__,__FILE__)
             ENDIF

             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF
             IF (INDEX(line,'LONGLIST') .NE. 0) THEN
                i1 = INDEX(line,'LONGLIST') + 8
                nloop1 = numspec/10+1
                IF (MOD(numspec,10) .EQ. 0.0_real_8) nloop1 = numspec/10

                nloop2 = numspec2/10+1
                IF (MOD(numspec2,10) .EQ. 0.0_real_8) nloop2 = numspec2/10

                DO i = 1,nloop1
                   IF (paral%io_parent)&
                        READ(iunit,err=20,END=20,fmt='(A)') line2
                   iend = 10
                   IF (i .EQ. nloop1 .AND.&
                        numspec .GT. INT(numspec/10)*10) THEN
                      iend = numspec - INT(numspec/10)*10
                   ENDIF

                   i1=1
                   DO ii = 1,iend
                      CALL readsi(line2,i1,iout,&
                           iatdlmn(ii+(i-1)*10+my_nat*(ncolvar-1)),erread)
                      i1=iout
                   ENDDO
                ENDDO

                IF (numspec2 .NE. 0) THEN
                   DO i = 1,nloop2
                      IF (paral%io_parent)&
                           READ(iunit,err=20,END=20,fmt='(A)') line2
                      iend = 10
                      IF (i .EQ. nloop2 .AND.&
                           numspec2 .GT. INT(numspec2/10)*10) THEN
                         iend = numspec2 - INT(numspec2/10)*10
                      ENDIF

                      i1 = 1
                      DO ii = 1,iend
                         CALL readsi(line2,i1,iout,&
                              iatdlmn(ii+(i-1)*10+numspec+my_nat*(ncolvar-1))&
                              ,erread)
                         i1=iout
                      ENDDO
                   ENDDO

                ENDIF

             ELSE
                DO i = 1,numspec
                   CALL readsi(line,i1,iout,&
                        iatdlmn(i+my_nat*(ncolvar-1)),erread)
                   i1=iout
                ENDDO
                IF (numspec2 .NE. 0) THEN
                   DO i = 1,numspec2
                      CALL readsi(line,i1,iout,&
                           iatdlmn(i+numspec+my_nat*(ncolvar-1)),erread)
                      i1=iout
                   ENDDO
                ENDIF
                i2=2
             ENDIF
          ENDIF
          ii=INDEX(line,'MIL')
          IF (ii.NE.0) THEN
             i1 = ii+3
             DO i = 1,3
                CALL readsi(line,i1,iout,atcvar(i+i2,ncolvar),erread)
                i1=iout
             ENDDO
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,A)') 'WARNING: Miller indeces not found,',&
                  '  default used 100'
             atcvar(1+i2,ncolvar)=1
             atcvar(2+i2,ncolvar)=0
             atcvar(3+i2,ncolvar)=0
          ENDIF
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==---- Angle between 2 planes, each defined by 3 given pnts  -----==
       ELSEIF (INDEX(line,'PLNANG').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 15
          iangcv(ncolvar) = 2
          i1=INDEX(line,'PLNANG')+6
          iangcv(ncolvar) = 1
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(5,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(6,ncolvar)=NAT_cpmd(idummy)
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==--------------------- Hydrogen Bonds Chain  --------------------==
       ELSEIF (INDEX(line,'HBONDCH').NE.0) THEN
          ! FIXME: AK 2005/05/27 convert to use Gromos ordering of atom indices
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 16
          i1=INDEX(line,'HBONDCH')+7
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)! NRES1
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)! NATTR
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)! NRES2
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)! NH
          i1=iout

          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)! EXP1_H
          i1=iout
          CALL readsi(line,i1,iout,atcvar(6,ncolvar),erread)! EXP2_H
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)! R0_H
          i1=iout
          CALL readsi(line,i1,iout,atcvar(7,ncolvar),erread)! EXP1_ATTR
          i1=iout
          CALL readsi(line,i1,iout,atcvar(8,ncolvar),erread)! EXP2_ATTR
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)! R0_ATTR

          IF (erread) CALL stopgm('M_COLVAR_INP',&
               'ERROR in reading HBONDCH line',& 
               __LINE__,__FILE__)

          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF
          iend = 0
          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+iend+my_nat*(ncolvar-1)),&
               ii=1,atcvar(1,ncolvar))
          ! write(6,*) (IATDLMN(II+IEND+MY_NAT*(NCOLVAR-1)),
          ! &                           II=1,ATCVAR(1,NCOLVAR))
          iend = iend + atcvar(1,ncolvar)
          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+iend+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))
          ! write(6,*) (IATDLMN(II+IEND+MY_NAT*(NCOLVAR-1)),
          ! &                           II=1,ATCVAR(2,NCOLVAR))
          iend = iend + atcvar(2,ncolvar)

          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+iend+my_nat*(ncolvar-1)),&
               ii=1,atcvar(3,ncolvar))
          ! write(6,*) (IATDLMN(II+IEND+MY_NAT*(NCOLVAR-1)),
          ! &                           II=1,ATCVAR(3,NCOLVAR))

          iend = iend + atcvar(3,ncolvar)
          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+iend+my_nat*(ncolvar-1)),&
               ii=1,atcvar(4,ncolvar))

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          IF (INDEX(line,'PARA').NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line2
             ii=INDEX(line2,'RCW')
             IF (ii.NE.0) THEN
                ia = ii+3
                CALL readsr(line2,ia,ie,rcw,   erread)
             ENDIF
             ii=INDEX(line2,'RCH')
             IF (ii.NE.0) THEN
                ia = ii+3
                CALL readsr(line2,ia,ie,rch,   erread)
             ENDIF
             ii=INDEX(line2,'NSHORT')
             IF (ii.NE.0) THEN
                ia = ii+6
                CALL readsi(line2,ia,ie,nshort,   erread)
             ENDIF
             ii=INDEX(line2,'NLONG')
             IF (ii.NE.0) THEN
                ia = ii+5
                CALL readsi(line2,ia,ie,nlong,   erread)
             ENDIF
             ii=INDEX(line2,'BBETA')
             IF (ii.NE.0) THEN
                ia = ii+6
                CALL readsr(line2,ia,ie,bbeta,   erread)
             ENDIF
             ii=INDEX(line2,'OPTDIST')
             IF (ii.NE.0) THEN
                ia = ii+7
                CALL readsr(line2,ia,ie,optdist,   erread)
             ENDIF
             ii=INDEX(line2,'LUPD')
             IF (ii.NE.0) THEN
                ia = ii+4
                CALL readsi(line2,ia,ie,lupd,   erread)
             ENDIF
             IF (erread) CALL stopgm('M_COLVAR_INP',&
                  'ERROR in reading H-CHAIN paraneters',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ==------------ DIFFERENCE OF COORDINATION NUMBERS -------------==
          ! ==  f =sum_i (1-(Rbi/R_0)^n)/(1-(Rbi/R_0)^(n+m))-
          ! ==     sum_i (1-(Rai/R_0)^n)/(1-(Rai/R_0)^(n+m))
       ELSEIF (INDEX(line,'DIFCOOR').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 17
          i1=INDEX(line,'DIFCOOR')+7
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') + 5
          ELSE
             specindex(ncolvar) = 1
          ENDIF

          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)

          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF

             IF (paral%io_parent)&
                  READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
                  ii=1,atcvar(3,ncolvar))
          ENDIF
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------ COORDINATION OF SECOND NEIGHBORS -------------==
          ! ==    f =sum_jik [(1+(Rij/R_0)^n)/(1+(Rij/R_0)^(n+m))*       ==
          ! ==                1+(Rjk/R_0)^n)/(1+(Rjk/R_0)^(n+m))]/NA/NB  == 

       ELSEIF (INDEX(line,'COOR_CHAIN').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 18
          i1=INDEX(line,'COOR_CHAIN')+10
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==------------ PRESENCE OF THE HYDRONIUM COMPLEX -----------==
          ! ==     F =(1/lambda)*log(sum_k exp{lambda*f(n_H)*n_O})      ==
          ! ==     f(n_h) = 1-(1-(n_h/rc_3)**n_3)/(1-(n_h/rc_3)**m_3)   ==
          ! ==     n_h =  sum_j(1-(r_kj/rc_1)**n_1)/(1-(r_kj/rc_1)**m_1)==
          ! ==     n_o =  sum_i(1-(r_ki/rc_2)**n_2)/(1-(r_ki/rc_2)**m_2)==

       ELSEIF (INDEX(line,'HYDRONIUM').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 19
          i1=INDEX(line,'HYDRONIUM')+9
          ! Defaults
          atcvar(3,ncolvar) = 6
          atcvar(4,ncolvar) = 8
          cvpar(1,ncolvar)  = 2.4_real_8
          atcvar(5,ncolvar) = 6
          atcvar(6,ncolvar) = 6
          cvpar(2,ncolvar)  = 5.5_real_8
          atcvar(7,ncolvar) = 10
          atcvar(8,ncolvar) = 15
          cvpar(3,ncolvar)  = 2.40_real_8
          cvpar(4,ncolvar)  = 10.0_real_8

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)

          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line2
             ii=INDEX(line2,'H:')
             IF (ii.NE.0) THEN
                i1 = ii+2
                CALL readsi(line2,i1,iout,atcvar(3,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(4,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(1,ncolvar),erread)
             ENDIF
             ii=INDEX(line2,'O:')
             IF (ii.NE.0) THEN
                i1 = ii+2
                CALL readsi(line2,i1,iout,atcvar(5,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(6,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(2,ncolvar),erread)
                i1=iout
             ENDIF
             ii=INDEX(line2,'FCUT:')
             IF (ii.NE.0) THEN
                i1 = ii+5
                CALL readsi(line2,i1,iout,atcvar(7,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(8,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(3,ncolvar),erread)
             ENDIF
             ii=INDEX(line2,'LAMBDA:')
             IF (ii.NE.0) THEN
                i1 = ii+7
                CALL readsr(line2,i1,iout,cvpar(4,ncolvar),erread)
             ENDIF
          ENDIF

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==------------ PRESENCE OF THE HYDRONIUM COMPLEX B -----------==
          ! ==     F =fc(nL)(1/lambda)*log(sum_k exp{lambda*M(n_H)*M(n_O)})      ==
          ! ==     M(n_h) = 1-(1-(n_h/rc_3)**n_M)/(1-(n_h/rc_3)**m_M)     ==
          ! ==     M(n_o) = 1-(1-(n_o/rc_4)**n_M)/(1-(n_h/rc_4)**m_M)     ==
          ! ==     n_h =  sum_j(1-(r_kj/rc_1)**n_1)/(1-(r_kj/rc_1)**m_1)  ==
          ! ==     n_o =  sum_i(1-(r_ki/rc_2)**n_2)/(1-(r_ki/rc_2)**m_2)  ==
          ! ==     nL  = sum_l(n_hl) 
          ! ==     n_hl = sum_j(1-(r_lj/rc_1)**n_1)/(1-(r_lj/rc_1)**m_1)  ==
          ! ==     fc(nL) = (1-(nL/rc_5)**n_L)/ (1-(nL/rc_5)**m_L)

       ELSEIF (INDEX(line,'HYD_CUT').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 27
          i1=INDEX(line,'HYD_CUT')+7
          ! Defaults

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(11,ncolvar),erread)

          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line2
          ii=INDEX(line2,'H:')
          IF (ii.NE.0) THEN
             i1 = ii+2
             CALL readsi(line2,i1,iout,atcvar(3,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(4,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(1,ncolvar),erread)
          ENDIF
          ii=INDEX(line2,'O:')
          IF (ii.NE.0) THEN
             i1 = ii+2
             CALL readsi(line2,i1,iout,atcvar(5,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(6,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(2,ncolvar),erread)
             i1=iout
          ENDIF
          ii=INDEX(line2,'M:')
          IF (ii.NE.0) THEN
             i1 = ii+2
             CALL readsi(line2,i1,iout,atcvar(7,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(8,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(3,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(4,ncolvar),erread)
             i1=iout
          ENDIF
          ii=INDEX(line2,'FCUT:')
          IF (ii.NE.0) THEN
             i1 = ii+5
             CALL readsi(line2,i1,iout,atcvar(9,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(10,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(5,ncolvar),erread)
          ENDIF
          ii=INDEX(line2,'LAMBDA:')
          IF (ii.NE.0) THEN
             i1 = ii+7
             CALL readsr(line2,i1,iout,cvpar(6,ncolvar),erread)
          ENDIF

          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF
          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(11,ncolvar))


          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)


          ! ==------------     DIPOLERHO          -------------==
          ! ==Norm of d = [1/Q sum_iq_i(R_i-r_0)]

       ELSEIF (INDEX(line,'DIPOLERHO').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 31
          i1=INDEX(line,'DIPOLERHO')+9


          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF

          IF (secondtdlmn .EQ. 0) THEN
             ALLOCATE(iqatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iqatdlmn)!,SIZE(iqatdlmn))
             secondtdlmn = 1
          ENDIF


          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)= NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)= idummy

          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))



          IF (paral%io_parent)&
               READ(iunit,*)  (iqatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))




          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)


          ! ==------------     DIPOLETHETA          -------------==
          ! ==Theta of d = [1/Q sum_iq_i(R_i-r_0)]

       ELSEIF (INDEX(line,'DIPOLETHA').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 32
          iangcv(ncolvar) = 1
          i1=INDEX(line,'DIPOLETHA')+9


          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF

          IF (secondtdlmn .EQ. 0) THEN
             ALLOCATE(iqatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iqatdlmn)!,SIZE(iqatdlmn))
             secondtdlmn = 1
          ENDIF

          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)= NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)= idummy

          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))

          IF (paral%io_parent)&
               READ(iunit,*)  (iqatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))


          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)


          ! ==------------     DIPOLEPHI          -------------==
          ! ==Phi of d = [1/Q sum_iq_i(R_i-r_0)]

       ELSEIF (INDEX(line,'DIPOLEPHI').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 33
          iangcv(ncolvar) = 2
          i1=INDEX(line,'DIPOLEPHI')+9


          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF

          IF (secondtdlmn .EQ. 0) THEN
             ALLOCATE(iqatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iqatdlmn)!,SIZE(iqatdlmn))
             secondtdlmn = 1
          ENDIF

          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)= NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)= idummy

          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))

          IF (paral%io_parent)&
               READ(iunit,*)  (iqatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(2,ncolvar))


          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

          ! ==------------ DISTANCE BETWEEN ION AND HYDRONIUM ----------==
          ! ==     F =(1/lambda)*log(sum_k exp{lambda*f(n_H)*n_O})      ==
          ! ==     f(n_h) = 1-(1-(n_h/rc_3)**n_3)/(1-(n_h/rc_3)**m_3)   ==
          ! ==     n_h =  sum_j(1-(r_kj/rc_1)**n_1)/(1-(r_kj/rc_1)**m_1)==
          ! ==     n_o =  sum_i(1-(r_ki/rc_2)**n_2)/(1-(r_ki/rc_2)**m_2)==

       ELSEIF (INDEX(line,'DIS_HYD').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 20
          i1=INDEX(line,'DIS_HYD')+7
          ! Defaults
          atcvar(4,ncolvar) = 6
          atcvar(5,ncolvar) = 8
          cvpar(1,ncolvar)  = 2.4_real_8
          atcvar(6,ncolvar) = 6
          atcvar(7,ncolvar) = 6
          cvpar(2,ncolvar)  = 5.5_real_8
          atcvar(8,ncolvar) = 10
          atcvar(9,ncolvar) = 15
          cvpar(3,ncolvar)  = 2.40_real_8
          cvpar(4,ncolvar)  = 100.0_real_8

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)

          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line2
             ii=INDEX(line2,'H:')
             IF (ii.NE.0) THEN
                i1 = ii+2
                CALL readsi(line2,i1,iout,atcvar(4,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(5,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(1,ncolvar),erread)
             ENDIF
             ii=INDEX(line2,'O:')
             IF (ii.NE.0) THEN
                i1 = ii+2
                CALL readsi(line2,i1,iout,atcvar(6,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(7,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(2,ncolvar),erread)
                i1=iout
             ENDIF
             ii=INDEX(line2,'FCUT:')
             IF (ii.NE.0) THEN
                i1 = ii+5
                CALL readsi(line2,i1,iout,atcvar(8,ncolvar),erread)
                i1=iout
                CALL readsi(line2,i1,iout,atcvar(9,ncolvar),erread)
                i1=iout
                CALL readsr(line2,i1,iout,cvpar(3,ncolvar),erread)
             ENDIF
             ii=INDEX(line2,'LAMBDA:')
             IF (ii.NE.0) THEN
                i1 = ii+7
                CALL readsr(line2,i1,iout,cvpar(4,ncolvar),erread)
             ENDIF
          ENDIF! cmb - this endif must stay here as in version 3.11

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)


          ! ==------------ DISTANCE BETWEEN ION AND HYDRONIUM COMPLEX B -----------==
          ! ==     F =fc(nL)[sum_k R_kI e^{Lambda*nHk*nOk}]/[sum_k e^{Lambda*nHk*nOk}]  ==
          ! ==     nHk =  sum_j(1-(r_kj/rc_1)**n_1)/(1-(r_kj/rc_1)**m_1)==
          ! ==     nOk =  sum_i(1-(r_ki/rc_2)**n_2)/(1-(r_ki/rc_2)**m_2)==
          ! ==     nL  = sum_l(n_hl) 
          ! ==     n_hl = sum_j(1-(r_lj/rc_1)**n_1)/(1-(r_lj/rc_1)**m_1)  ==
          ! ==     fc(nL) = (1-(nL/rc_5)**n_L)/ (1-(nL/rc_5)**m_L)

       ELSEIF (INDEX(line,'HYD_D_CUT').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 28
          i1=INDEX(line,'HYD_D_CUT')+9
          ! Defaults

          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(10,ncolvar),erread)

          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line2
          ii=INDEX(line2,'H:')
          IF (ii.NE.0) THEN
             i1 = ii+2
             CALL readsi(line2,i1,iout,atcvar(4,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(5,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(1,ncolvar),erread)
          ENDIF
          ii=INDEX(line2,'O:')
          IF (ii.NE.0) THEN
             i1 = ii+2
             CALL readsi(line2,i1,iout,atcvar(6,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(7,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(2,ncolvar),erread)
             i1=iout
          ENDIF
          ii=INDEX(line2,'FCUT:')
          IF (ii.NE.0) THEN
             i1 = ii+5
             CALL readsi(line2,i1,iout,atcvar(8,ncolvar),erread)
             i1=iout
             CALL readsi(line2,i1,iout,atcvar(9,ncolvar),erread)
             i1=iout
             CALL readsr(line2,i1,iout,cvpar(3,ncolvar),erread)
          ENDIF
          ii=INDEX(line2,'LAMBDA:')
          IF (ii.NE.0) THEN
             i1 = ii+7
             CALL readsr(line2,i1,iout,cvpar(4,ncolvar),erread)
          ENDIF

          IF (firstdlmn .EQ. 0) THEN
             ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(iatdlmn)!,my_nat*nnvar)
             firstdlmn = 1
          ENDIF
          IF (paral%io_parent)&
               READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
               ii=1,atcvar(10,ncolvar))


          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
          ! ==------------ LOCALIZATION OF THE SPIN DENSITY  -----------==
       ELSEIF (INDEX(line,'SPIN').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 21
          IF (.NOT. cntl%tlsd)  CALL stopgm('M_COLVAR_INP',&
               ' ERROR: LSD is required for SPIN CV',& 
               __LINE__,__FILE__)
          IF (.NOT. lmeta%tlocalizespin) THEN
             ALLOCATE(rcc0(3,ncolvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(rcc0)!,3*nnvar)
             ALLOCATE(icv_spin(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_spin)!,nnvar)
             lmeta%tlocalizespin=.TRUE.
          ENDIF
          icv_spin(ncolvar)=ncolvar

          i1=INDEX(line,'SPIN')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          IF (erread) GOTO 22
          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
       ELSEIF (INDEX(line,'VOLVAR').NE.0) THEN
          IF (.NOT. cntl%tprcp) THEN
             CALL stopgm('M_COLVAR_INP',&
                  ' ERROR: Parr.-Rahm. required for variable cell',& 
                  __LINE__,__FILE__)
          ENDIF
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 22
          IF (.NOT. lmeta%tcvcell) THEN
             ALLOCATE(icv_cell(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_cell)!,nnvar)
             lmeta%tcvcell=.TRUE.
          ENDIF

          icv_cell(ncolvar)=ncolvar

          i1=INDEX(line,'VOLVAR')+6
          ! CALL READSI(LINE,I1,IOUT,ATCVAR(1,NCOLVAR),ERREAD)
          ! IF(ERREAD) GOTO 22

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
       ELSEIF (INDEX(line,'CELLSIDE').NE.0) THEN
          IF (.NOT. cntl%tprcp) THEN
             CALL stopgm('M_COLVAR_INP',&
                  ' ERROR: Parrnello-Rahman required for variable cell',& 
                  __LINE__,__FILE__)
          ENDIF
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 23
          IF (.NOT. lmeta%tcvcell) THEN
             ALLOCATE(icv_cell(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_cell)!,nnvar)
             lmeta%tcvcell=.TRUE.
          ENDIF

          icv_cell(ncolvar)=ncolvar

          i1=INDEX(line,'CELLSIDE')+8
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          IF (erread) GOTO 22

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)
       ELSEIF (INDEX(line,'CELLANGLE').NE.0) THEN
          IF (.NOT. cntl%tprcp) THEN
             CALL stopgm('M_COLVAR_INP',&
                  ' ERROR: Parrinello-Rahman required for variable cell',& 
                  __LINE__,__FILE__)
          ENDIF
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 24
          IF (.NOT. lmeta%tcvcell) THEN
             ALLOCATE(icv_cell(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_cell)!,nnvar)
             lmeta%tcvcell=.TRUE.
          ENDIF

          icv_cell(ncolvar)=ncolvar

          i1=INDEX(line,'CELLANGLE')+9
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          IF (erread) GOTO 22

          CALL param_meta(line,tad_scf,vbound,ibound,&
               cv_dyn_0,initial_value,&
               ncolvar,cscl_fac,kharm,cv_mass,lmeta%lextlagrange)

       ENDIF

       IF (erread) GOTO 20

       IF (ncolvar.GT.nnvar) THEN
          ! Problem the number of coll.variables NCOLVAR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'META_COLVAR_INP! # OF GIVEN COL. VAR= ',NNVAR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'META_COLVAR_INP! # OF FOUND COL. VAR=',NCOLVAR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')&
               'THE NUMBER OF GIVEN COL. VAR. IS TOO SMALL.'
          CALL stopgm('M_COLVAR_INP',&
               'ERROR WHILE READING FIX STRUCTURES',& 
               __LINE__,__FILE__)
       ELSEIF (ncolvar.LT.nnvar) THEN
          GOTO 12
       ELSE
          GOTO 11
       ENDIF
       ! ==------------ MAX NUMBER OF METADYNAMICS STEPS ----------------==
    ELSEIF (INDEX(line,'METASTEPNUM').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%i_meta_max
       GOTO 10
       ! ==-------- RESTART META DYNAMICS FROM METASTEPI_META_RES -------==
    ELSEIF (INDEX(line,'META_RESTART').NE.0) THEN
       lmeta%meta_restart = .TRUE.
       IF (INDEX(line,'RFILE').NE.0) THEN
          lmeta%tresfile = .TRUE.
          ALLOCATE(dummy(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (lmeta%tmulti) THEN
             CALL rmtdresm(50,dummy,.FALSE.,1)
          ELSEIF (tmw)THEN
             DO iwalk1=1,mwi%nwalk
                CALL mw_filename('MTD_RESTART_',fmtdres,iwalk1)
                CALL rmtdres_mw(50,dummy,.FALSE.,1,iwalk1)
             ENDDO
          ELSE
             CALL rmtdres(50,dummy,.FALSE.,1)
          ENDIF
          DEALLOCATE(dummy,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSE
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt=*) imeta%i_meta_res
       ENDIF
       GOTO 10
       ! ==---------------- MORE THAN ONE SET OF CV ---------------------==
    ELSEIF(INDEX(line,'MULTI').NE.0 .AND.&
         INDEX(line,'NUM') .NE.0) THEN
       IF (.NOT. lmeta%tmulti) CALL stopgm('M_COLVAR_INP',&
            'Multy Meta Dynamics not active. ',& 
            __LINE__,__FILE__)
       ALLOCATE(ncvsys(nsubsys),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            READ(iunit,*) (ncvsys(i),i = 1,nsubsys)
       GOTO 10
       ! ==------------- DIMENSIONS OF HILLS IN CV SPACE-----------------==
    ELSEIF (INDEX(line,'HILLS').NE.0) THEN
       IF (INDEX(line,'OFF').NE.0) THEN
          lmeta%lhills=.FALSE.
       ELSE
          lmeta%lhills = .TRUE.
          IF (INDEX(line,'LOREN').NE.0) THEN
             lmeta%hlore = .TRUE.
          ELSE IF (INDEX(line,'SPHERE').NE.0) THEN
             lmeta%sphere= .TRUE.
          ELSEIF (INDEX(line,'RATIO').NE.0) THEN
             lmeta%hratio = .TRUE.
             lmeta%hshift = .FALSE.
             ij=INDEX(line,'POW')
             IF (ij.NE.0) THEN
                ia = ij+3
                CALL readsr(line,ia,ie,rmeta%expup,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%expdown,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%fboost,   erread)
                IF (erread) GOTO 20
             ENDIF
          ELSEIF (INDEX(line,'SHIFT').NE.0) THEN
             lmeta%hshift = .TRUE.
             lmeta%hratio = .FALSE.
             ij=INDEX(line,'RCUT')
             IF (ij.NE.0) THEN
                ia = ij+4
                CALL readsr(line,ia,ie,r0_shift,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%fboost,   erread)
                IF (erread) GOTO 20
             ENDIF

          ENDIF
          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             IF (lmeta%tmulti) THEN
                DO i = 1, nsubsys
                   IF (paral%io_parent)&
                        READ(iunit,err=20,END=20,fmt=*)&
                        hwm(i),hhm(i),hlhm(i),hthm(i)
                   hvm0(i) = hhm(i) * rmeta%hllw**REAL(ncvsys(i),kind=real_8)
                ENDDO
             ELSE
                ia = ii+1
                CALL readsr(line,ia,ie,rmeta%hllw,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%hllh,   erread)
                IF (erread) GOTO 20
                rmeta%hvol0 = rmeta%hllh * (rmeta%hllw)**REAL(ncolvar,kind=real_8)
             ENDIF
          ENDIF
       ENDIF
       GOTO 10
       ! ==---------------------- TUNING PARAMETERS ---------------------==
    ELSEIF (INDEX(line,'TUNING').NE.0) THEN
       IF (INDEX(line,'OFF').NE.0) THEN
          lmeta%ltune_hh     = .FALSE.
          lmeta%ltune_cscl   = .FALSE.
       ENDIF
       IF (INDEX(line,'HHEIGHT').NE.0) THEN
          lmeta%ltune_hh     = .TRUE.
          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             ia = ii+1
             CALL readsr(line,ia,ie,rmeta%hlow,   erread)
             IF (erread) GOTO 20
             ia = ie
             CALL readsr(line,ia,ie,rmeta%htop,   erread)
             IF (erread) GOTO 20
             IF (rmeta%htop.LT.rmeta%hlow) THEN
                temp = rmeta%hlow
                rmeta%hlow = rmeta%htop
                rmeta%htop = temp
             ENDIF
          ENDIF
       ENDIF
       GOTO 10
       ! ==---------MAX # OF DYN. STEPS BETWEEN 2 METASTEPS--------------==
    ELSEIF(INDEX(line,'MAXSTEPNUM').NE.0 .AND.&
         INDEX(line,'INTERMETA') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inter_hill_max
       GOTO 10
       ! ==---------MIN # OF DYN. STEPS BETWEEN 2 METASTEPS--------------==
    ELSEIF(INDEX(line,'MINSTEPNUM').NE.0 .AND.&
         INDEX(line,'INTERMETA') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inter_hill
       GOTO 10
       ! ==------------ INTERVAL BETWEEN DISPLACEMENT CHECKS-------------==
    ELSEIF(INDEX(line,'CHECK') .NE.0 .AND.&
         INDEX(line,'DELAY').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%icheck
       GOTO 10
       ! ==------------------------TOLL_FCNST----------------------------==
    ELSEIF(INDEX(line,'MOVEMENT') .NE.0 .AND.&
         INDEX(line,'CHECK').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) toll_avcv
       GOTO 10
       ! ==----------- STORE RESTART AND TRAJECTORIES -------------------==
    ELSEIF (INDEX(line,'METASTORE') .NE.0) THEN
       ! don t use TR_FREQ if NO TRAJECTORY is given.
       IF ((INDEX(line,'TRAJ').NE.0).AND.(INDEX(line,'NO').NE.0))THEN
          lmeta%ttrajovrr=.FALSE.
       ENDIF
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%st_freq,imeta%tr_freq,imeta%qw_freq
       GOTO 10
       ! ==----------- DEF. OF BOUNDARIES IN CV. SPACE  ---------------==
    ELSEIF(INDEX(line,'CVSPACE').NE.0 .AND.&
         INDEX(line,'BOUND') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inum
       IF (erread) GOTO 20
       DO i = 1,inum
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line2
          IF (paral%io_parent)&
               READ(line2,err=20,END=20,fmt=*) ind,vbar,r_wall,cpos_0
          ii=INDEX(line2,'EPS')
          IF (ii.NE.0) THEN
             ia = ii+3
             CALL readsr(line2,ia,ie,vbound(4,ind),erread)
          ELSE
             vbound(4,ind) = 0.01_real_8
          ENDIF
          IF (ibound(ind) .NE. 0 .OR. ind .GT. ncolvar)&
               CALL stopgm('M_COLVAR_INP',&
               'ERROR IN  BOUNDARIES OF COLLECTIVE VAR. ',& 
               __LINE__,__FILE__)
          ibound(ind) = 2
          vbound(1,ind) = vbar
          vbound(2,ind) = r_wall
          vbound(3,ind) = cpos_0
       ENDDO
       GOTO 10
       ! ==----------- COLLECTIVE VARIABLE OUTPUT  -------------------==
    ELSEIF (INDEX(line,'MONITOR') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%wcv_freq
       lmeta%tcvmonitor = .TRUE.
       GOTO 10
       ! ==----- MAX ELECT. KIN. ENERGY (when above QUENCH BO) ----------==
    ELSEIF (INDEX(line,'MAXKINEN').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       CALL readsr(line,1,iout,rmeta%tolkin,erread)
       GOTO 10
       ! ==-------------- EXTENDED LAGRANGIAN INFO ------ ---------------==
    ELSEIF (INDEX(line,'LAGRANGE').NE.0) THEN
       lmeta%lextlagrange=.TRUE.
       IF (INDEX(line,'TEMPERATURE').NE.0) THEN
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line
          CALL readsr(line,1,iout,cv_temp0,erread)
          IF (erread) CALL stopgm('M_COLVAR_INP','ERR READ TEMP0 ',& 
               __LINE__,__FILE__)! GOTO 20
          GOTO 10
       ELSEIF (INDEX(line,'TEMPCONT').NE.0) THEN
          IF (INDEX(line,'GLOBAL').NE.0) THEN
             ltcglobal = .TRUE.
          ELSE
             lcvtc = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line
             CALL readsr(line,1,iout,cv_temp,erread)
             IF (erread)&
                  CALL stopgm('M_COLVAR_INP','ERR READ TEMP ',& 
                  __LINE__,__FILE__) ! GOTO 20
             i1=iout+1
             CALL readsr(line,i1,iout,cv_dtemp,erread)
             IF (erread)&
                  CALL stopgm('M_COLVAR_INP','ERR READ DTEMP ',& 
                  __LINE__,__FILE__) ! GOTO 20
          ENDIF
          GOTO 10
       ELSEIF (INDEX(line,'TSCALE').NE.0) THEN
          tcvscale = .TRUE.
          GOTO 10
       ELSEIF (INDEX(line,'LANGEVIN').NE.0) THEN
          tcvlangevin = .TRUE.
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line
          CALL readsr(line,1,iout,cv_langevintemp,erread)
          IF (erread)&
               CALL stopgm('M_COLVAR_INP','ERR READ CV_LANGEVINTEMP',& 
               __LINE__,__FILE__) ! GOTO 20
          i1=iout+1
          CALL readsr(line,i1,iout,cv_langamma,erread)
          IF (erread)&
               CALL stopgm('M_COLVAR_INP','ERR READ CV_LANGAMMA',& 
               __LINE__,__FILE__) ! GOTO 20
          GOTO 10
       ELSEIF (INDEX(line,'KADAPT').NE.0) THEN
          lkfix = .FALSE.
          GOTO 10
       ELSE
          GOTO 10
       ENDIF
       ! ==-------------------- CV  ANALYSIS -------------- ---------------==
    ELSEIF (INDEX(line,'ANALYS').NE.0) THEN
       lmeta%tcvanalysis  = .TRUE.
       GOTO 10
       ! ==-------------------- Volume of hills kept constant --------------==
    ELSEIF (INDEX(line,'HILLVOLUME').NE.0) THEN
       lmeta%thillvol = .TRUE.
       GOTO 10
    ELSEIF (INDEX(line,'SKIPHILL').NE.0)THEN
       lmeta%skiphill_mw=.TRUE.
       GOTO 10
    ELSEIF (INDEX(line,'RANDWALK').NE.0)THEN
       lmeta%randwalk=.TRUE.
       GOTO 10
       ! =------ RESTRAIN THE CELL VOLUME IN THE RANGE VOLMIN VOLMAX ----==
    ELSEIF (INDEX(line,'RESTR').NE.0.AND.INDEX(line,'VOLU').NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) fmin,fmax,mdcellr%vbarrier
       mdcellr%volmin = prcp_com%omega0*fmin
       mdcellr%volmax = prcp_com%omega0*fmax
       tvolbound = .TRUE.
       GOTO 10
    ELSE
       GOTO 10
    ENDIF  ! file beginning

20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING METADYNAMICS VARIABLES '
    IF (paral%io_parent)&
         WRITE(6,*) '       last line read: '
    IF (paral%io_parent)&
         WRITE(6,*)  line
    CALL stopgm('M_COLVAR_INP',' ',& 
         __LINE__,__FILE__)   
22  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING SCAL. FAC. or MASS or KH.'
    CALL stopgm('M_COLVAR_INP',' ',& 
         __LINE__,__FILE__)

30  CONTINUE
    ! ==--------------------------------------------------------------==
    ! TESTS

    IF (ncolvar .EQ. 0) CALL stopgm('M_COLVAR_INP',&
         '# OF COLLECTIVE VARIABLES IS ZERO',& 
         __LINE__,__FILE__)
    DO ic = 1,ncolvar
       IF (tad_scf(ic)) THEN
          lmeta%ltune_cscl = .TRUE.
          GOTO 24
       ENDIF
    ENDDO
24  CONTINUE
    IF (lmeta%tcvanalysis) THEN
       lmeta%lextlagrange =.TRUE.
       lmeta%ltune_hh     = .FALSE.
       lmeta%ltune_cscl   = .FALSE.
       lmeta%lhills       = .FALSE.
       lcvtc        = .FALSE.
       CALL zeroing(ibound)!,nnvar)
       CALL dcopy(nnvar,1._real_8,0,cv_mass,1)
       CALL zeroing(kharm)!,nnvar)
       toll_avcv  = 0.0_real_8
       cv_temp0 = 0.0_real_8
       inter_hill = 10
    ENDIF
    IF (rmeta%hllw .EQ. 0.0_real_8) THEN
       CALL stopgm('M_COLVAR_INP',&
            'Hill-width = 0 is too dangerous',& 
            __LINE__,__FILE__)
    ENDIF
    IF (imeta%qw_freq.LE.0) THEN
       CALL stopgm('M_COLVAR_INP',&
            'Illegal QUENCH BO frequency',& 
            __LINE__,__FILE__)
    ENDIF
    IF (lmeta%thillvol .AND. .NOT.lmeta%ltune_hh ) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' WARNING: the hill tuning is off,',&
            ' height and wdth are constant'
       lmeta%thillvol = .FALSE.
    ENDIF
    IF (tvolbound .AND. .NOT. cntl%tprcp) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' WARNING: if not Parr.-Rahm. ',&
            ' the volume of the cell is constant'
       tvolbound = .FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Print initialization

    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(72("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(1x,A,A)')&
         '***             META DYNAMICS OF COLLECTIVE VARIABLES',&
         '             ***'
    IF (paral%io_parent)&
         WRITE(6,'(4x,A,I5)') '- Number of Collective Variables =',&
         ncolvar
    IF (lmeta%tmulti) THEN
       IF (paral%io_parent)&
            WRITE(6,'(4x,A,I5)') '- Multiple meta dynamics: # of CV sets =',&
            nsubsys
       IF (paral%io_parent)&
            WRITE(6,'(12x,A,10(I4))') ' # of CV per set =',&
            (ncvsys(i),i=1,nsubsys)
    ENDIF
    IF ((mwi%nwalk.GT.1).AND.paral%io_parent)&
         WRITE(6,'(4x,A,I5)')&
         '- Multiple walker metadnamics: # of walkers=',mwi%nwalk
    IF (paral%io_parent)&
         WRITE(6,'(4x,A)') '- Names of Output Files :'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A)')&
         'colvar_mtd   : CV values and scaling factors'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A/6x,A)')&
         'parvar_mtd  : norm of the displacement along the traj.',&
         '              transversal hill width (fixed),',&
         '              hills heights HLLH'
    IF (lmeta%lextlagrange) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)')&
            'istvar_mtd: istantaneous CV and difference wrt dynamic CV'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/6x,A/6x,A)')&
            'disvar_mtd  : CV displacements and ',&
            '               diffusivities of the CV ',&
            '               force constants for the coupling potential,'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/6x,A/6x,A)')&
            'forfac_mtd  : forces due to coupling potential,',&
            '               forces due to accumulation of hills',&
            '               forces due to boundary potential.'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)') 'velvar_mtd  : CV velocities.'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/10x,A/10x,A)')&
            'enevar_mtd  :  ions Temp., el. Ekin, CV Ekin (EKCV), ',&
            '         coupling pot. (VH), gaussians pot. (VG),',&
            '         i+el. pot. energy (EKS), H0+EKCV+VH,H0+EKCV+VH+VG'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)')&
            'scavar_mtd: scaled CV and scaled diffusivities'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)')&
            'disvar_mtd  :  displacements and diffusivities of the CV.'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/6x,A)')&
            'forfac_mtd  : forces due to  accumulation of hills,',&
            '               forces due to boundary potential.'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/6x,A/9x,A)')&
            'enevar_mtd  :  ions Temp., el. Ekin, ',&
            '               gaussians pot. (VG),',&
            '            ions+el. pot. energy(EKS), H0+VG '

    ENDIF

    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- Collective variables written every ', imeta%wcv_freq, ' MD STEPS'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- RESTART File Saved every ', imeta%st_freq, ' Metasteps'
    IF ((lmeta%ttrajovrr).AND.paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- TRAJECTORY File Appended every', imeta%tr_freq,' Metasteps'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- QUENCH BO performed every', imeta%qw_freq, ' Metasteps'
    IF (cntl%tc .AND. rmeta%tolkin .LT. cntr%ekinw+cntr%toll) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,F12.6)')&
            '- QUENCH BO performed when EL. Ekin is > ', rmeta%tolkin
    ELSEIF (rmeta%tolkin.NE.-1.0_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,F12.6)')&
            '- QUENCH BO performed when EL. Ekin is > ', rmeta%tolkin
    ENDIF
    IF ((lmeta%meta_restart).AND.paral%io_parent)&
         WRITE(6,'(6x,A,1X,I8)')&
         '- Restart from a previous Metadynamics at step',imeta%i_meta_res
    IF (lmeta%lhills) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)')&
            '- Time Dependent Potential Active: '

       IF (lmeta%tmulti) THEN
          DO i = 1,nsubsys
             IF (paral%io_parent)&
                  WRITE(6,'(12x,A,I4,A,F9.5,A,F9.5)')&
                  ' set ', i, ' Transversal Width = ',&
                  hwm(i),'    Initial Height = ',hhm(i)
          ENDDO
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(15x,A,F9.5,A,F9.5)') ' Transversal Width = ',&
               rmeta%hllw,'    Initial Height = ',rmeta%hllh
       ENDIF
       IF (lmeta%hshift) THEN
          IF (paral%io_parent)&
               WRITE(6,'(6x,A/15x,A,F8.3,A,f8.3)')&
               '- Hills Shape : Truncated Gaussians',&
               ' Rcut = ',rmeta%rshift,'     Width Factor = ',rmeta%fboost
       ELSEIF (lmeta%hratio) THEN
          IF (paral%io_parent)&
               WRITE(6,'(6x,A/15x,A,F8.3,F8.3,A,f8.3)')&
               '- Hills Shape : Rational Function',&
               ' Exponents : ',rmeta%expup,rmeta%expdown,&
               '     Width Factor = ',rmeta%fboost
       ENDIF
       IF (lmeta%ltune_hh) THEN
          IF (paral%io_parent)&
               WRITE(6,'(6x,A)')&
               '- Hills Height tuned on the Underlying potential'
          IF (lmeta%tmulti) THEN
             DO i = 1,nsubsys
                IF (paral%io_parent)&
                     WRITE(6,'(12x,A,I4,A,F10.6,A,f10.6,A)') ' set ',i,&
                     'min. Height ',hlhm(i)*au_kcm,' Kcal,    max Height ',&
                     hthm(i)*au_kcm,' Kcal.'
             ENDDO
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(15x,A,F10.6,A,f10.6,A/)')&
                  'min. Height ',rmeta%hlow*au_kcm,' Kcal,    max Height ',&
                  rmeta%htop*au_kcm,' Kcal.'
          ENDIF
          IF (lmeta%thillvol) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(12x,A)')&
                  ' Hill volume Height*(Width)^(NCOLVAR) = constant'
          ENDIF
       ENDIF
       IF ((lmeta%ltune_cscl).AND.paral%io_parent)&
            WRITE(6,'(6x,A)')&
            '- Some scaling factors might be adjusted along the run'
    ENDIF

    IF (lmeta%lhills) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,F12.6)')&
            '- Minimum. CV Displacement Required to Place New a Hill',&
            toll_avcv
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,I5,A)')&
            '- Update of Gaussian Pot. every ',&
            inter_hill,' MD STEPS (WHEN CV MOVED ENOUGH)'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,I10)')&
            '- MAX. NUMBER OF MD STEPS BEFORE NEXT UPDATE',inter_hill_max
    ENDIF
    IF (lmeta%tcvanalysis) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,A)') '- MD RUN FOR THE ANALYSIS OF CV BEHAVIOR ',&
            'without additional potentials'
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,I5,A)')&
            '- Output written every ',&
            inter_hill,' MD STEPS'
    ENDIF
    IF (lmeta%lextlagrange) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/6x,A,A)')&
            '***         CV DYNAMICS VIA EXTENDED LAGRANGIAN',&
            '         ***'
       IF (paral%io_parent)&
            WRITE(6,'(10x,A,f9.3,A)')&
            '- Initial CV Temperature  ',cv_temp0,' K'
       IF ((lcvtc).AND.paral%io_parent)&
            WRITE(6,'(10x,A/6x,A,f9.3,A,A,f9.3,A)')&
            '- CV Temperature controlled by velocity rescaling',&
            '        mean T = ',cv_temp,' K,',&
            '      tolerance= ',cv_dtemp,' K'
       IF ((tcvlangevin).AND.paral%io_parent)&
            WRITE(6,'(10x,A/6x,A,f9.3,A,A,f9.6,A)')&
            '- CV Temperature controlled by Langevin Dynamics',&
            '        mean T = ',cv_langevintemp,' K,',&
            '      Gamma = ',cv_langamma,' a. u.'
       IF ((.NOT.lkfix).AND.paral%io_parent)&
            WRITE(6,'(10x,A,A/)')&
            '- Force Constants for Coupling Pot.',&
            ' Tuned on Underlying Pot. Surf.'
    ENDIF

    IF (tvolbound) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A/A,f12.6,A/A,f12.6,A)')&
            '- PR MD WITH CONTROL ON VOLUME: ',&
            '                   min Volume = ', mdcellr%volmin,' a.u.^3',&
            '                   max Volume = ', mdcellr%volmax, 'a.u.^3'

    ENDIF

    IF (paral%io_parent)&
         WRITE(6,'(/15x,A)') '   >>>>>>> COLLECTIVE VARIABLE INFO <<<<<<<'
    IF (paral%io_parent)&
         WRITE(6,'(2x,A,5x,A,5x,A,20x,A,11x,A)')&
         'TYPE','ATOMS AND SPEC.',&
         'PARAMETERS',&
         'SCALING FAC.','BOX POT '

    ipgrpa1=0
    ipgrpb1=0
    DO i=1,ncolvar
       ityp=tycvar(i)
       ! COORD
       IF (ityp.EQ.6) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4X,I5,15X,2f8.4,4x,3F8.4,5x,4F8.4)') styp(ityp),&
               nat_grm(atcvar(1,i)),cvpar(1,i),cvpar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
       ELSEIF (ityp.EQ.30) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4X,3I5,15X,2f8.4,4x,3F8.4,5x,4F8.4)') styp(ityp),&
               nat_grm(atcvar(1,i)),nat_grm(atcvar(2,i)),atcvar(2,i),&
               cvpar(1,i),cvpar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ! COORSP
       ELSEIF (ityp.EQ.8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,2I5,15X,2f8.4,5x,3F6.2,5x,4F8.4)') styp(ityp),&
               nat_grm(atcvar(1,i)),atcvar(2,i),cvpar(1,i),cvpar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
       ELSEIF (ityp.EQ.29) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,3F16.2,4F8.4)') styp(ityp),atcvar(1,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ipgrpb3=0
          DO k=1,atcvar(1,i)
             ipgrpa1=ipgrpa1+1
             IF (paral%io_parent)&
                  WRITE(6,'(I8,F16.5,I8)')iatcnga(ipgrpa1),rccnga(ipgrpa1),&
                  natcngb(ipgrpa1)
             DO kk=1,natcngb(ipgrpa1)
                ipgrpb3=ipgrpb3+1
                IF (paral%io_parent)&
                     WRITE(6,'(I8)')iatcngb(ipgrpb1+ipgrpb3)
             ENDDO
          ENDDO
          ipgrpb1=ipgrpb1+atcvar(1,i)*max_natcngb
          ! DIPOLERHO 
       ELSEIF (ityp.EQ.31) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,1I5,18X,1I5, 18x,3F6.2,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ! DIPOLETHETA
       ELSEIF (ityp.EQ.32) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,1I5,18X,1I5,18x,3F6.2,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ! DIPOLEPHI
       ELSEIF (ityp.EQ.33) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,1I5,18X,1I5,18x,3F6.2,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ! VCOORS 
       ELSEIF (ityp.EQ.34) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2X,1I5,20X,2f8.4,5x,3F6.2,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),cvpar(1,i),cvpar(2,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
          ! RMSD_AB, DISPL1, COORDIS
       ELSEIF (ityp .GE.9 .AND. ityp.LE.11) THEN
          IF (cvpar(2,i) .NE. 0.0_real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,2I5,5x,A,2I4,2f7.3,2x,3F8.4,3x,4F7.3)')&
                  styp(ityp),nat_grm(atcvar(1,i)),atcvar(2,i),&
                  '2 SHELL',atcvar(3,i),atcvar(4,i),&
                  cvpar(1,i),cvpar(2,i),&
                  cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
                  (vbound(k,i),k=1,4)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,2I5,13x,2I4,f7.3,7x,3F8.4,5x,4F8.4)')&
                  styp(ityp),nat_grm(atcvar(1,i)),atcvar(2,i),&
                  atcvar(3,i),atcvar(4,i),&
                  cvpar(1,i),cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
                  (vbound(k,i),k=1,4)
          ENDIF
       ELSEIF (ityp .EQ. 17) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I4,11x,2I3,f7.3,7x,3F8.4,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cvpar(1,i),cscl_fac(1,i),cscl_fac(2,i),&
               cscl_fac(3,i),(vbound(k,i),k=1,4)
       ELSEIF (ityp .EQ. 18) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I4,8x,2I3,2f7.3,3x,3F8.4,5x,4F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cvpar(1,i),cscl_fac(1,i),cscl_fac(2,i),&
               cscl_fac(3,i),(vbound(k,i),k=1,4)
       ELSEIF (ityp .EQ. 12) THEN
          numspec = atcvar(1,i)
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') numspec
          CALL xstring(chnum,ia,ie)
          lineform ='(A,2X,I4,A,'//chnum(ia:ie)//&
               'I4,10x,3f8.3,5x,4F8.4)'
          IF (paral%io_parent)&
               WRITE(6,lineform) styp(ityp),atcvar(1,i),&
               ' Species ',(atcvar(k,i),k=1,numspec),cscl_fac(1,i),&
               cscl_fac(2,i),cscl_fac(3,i),(vbound(k,i),k=1,4)
       ELSEIF (ityp.EQ.13 .AND. specindex(i) .EQ. 1) THEN
          numspec = atcvar(1,i)
          numspec2 = atcvar(2,i)
          numtot  = numspec+numspec2
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') numspec
          CALL xstring(chnum,ia,ie)
          IF (numspec2 .NE. 0) THEN
             IF (paral%io_parent)&
                  WRITE(chnum2,'(I5)') numspec2
             CALL xstring(chnum2,iaa,iee)
             lineform ='(A,2X,I4,A9,'//chnum(ia:ie)//'I4,A4,I4,A9,'&
                  //chnum2(ia:ie)//'I4,A15,3I3,2X,3f8.4)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),atcvar(1,i),&
                  ' Species: ',(atcvar(k+2,i),k=1,numspec),' vs ',atcvar(2,i),&
                  ' Species: ',(atcvar(k+2+numspec,i),k=1,numspec2),&
                  ' Miller Index: ',(atcvar(k+2+numtot,i),k=1,3),&
                  cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
          ELSE
             lineform ='(A,2X,I4,A9,'//chnum(ia:ie)//&
                  'I4,A7,A15,3I3,2X,3f8.4)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),atcvar(1,i),&
                  ' Species: ',(atcvar(k+2,i),k=1,numspec),' vs all',&
                  ' Miller Index: ',(atcvar(k+2+numtot,i),k=1,3),&
                  cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
          ENDIF
       ELSEIF (ityp.EQ.13 .AND. specindex(i) .NE. 1) THEN
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') atcvar(1,i)
          CALL xstring(chnum,ia,ie)
          numspec2 = atcvar(2,i)
          IF (numspec2 .NE. 0) THEN
             IF (paral%io_parent)&
                  WRITE(chnum2,'(I5)') atcvar(2,i)
             CALL xstring(chnum2,iaa,iee)
             lineform ='(A9,2X,A5,'//chnum(ia:ie)//'I4,A8,'&
                  //chnum2(iaa:iee)//'I4,A10,3I3,2X,A4,3f8.4)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),' at.',&
                  (iatdlmn(k+my_nat*(i-1)),k=1,atcvar(1,i)),&
                  ' vs at.',(iatdlmn(atcvar(1,i)+k+my_nat*(i-1)),&
                  k=1,atcvar(2,i)),&
                  ' Mil. In.',(atcvar(k+2,i),k=1,3),' SCF',&
                  cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
          ELSE
             lineform ='(A,2X,A5,'//chnum(ia:ie)//&
                  'I4,A8,A10,3I3,2X,A4,3f8.4)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),' at.',&
                  (iatdlmn(k+my_nat*(i-1)),k=1,atcvar(1,i)),&
                  ' vs all ',' Mil. In.',(atcvar(k+2,i),k=1,3),&
                  ' SCF', cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
          ENDIF

       ELSEIF (ityp.EQ.15) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,6I5,13x,3F8.4,5x,4F8.4)') styp(ityp),&
               (atcvar(k,i),k=1,6),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
       ELSEIF (ityp.EQ.16) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4I5,23x,3F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
          IF (paral%io_parent)&
               WRITE(6,'(15x,A,2I3,A,F8.4,A,2I3,A,F8.4)')&
               'expH ',atcvar(5,i),atcvar(6,i),' R0_H ',cvpar(1,i),&
               '   expA ',atcvar(7,i),atcvar(8,i),' R0_A ',cvpar(2,i)
          IF (paral%io_parent)&
               WRITE(6,'(15x,A,2I3,A,2F8.4,A,F8.4,A,I4)')&
               'max & min L',nshort,nlong,'neigh. dist.',&
               rcw,rch,'opt. H-A dist. ',optdist,'UPD. NEIGH. ',lupd
       ELSEIF (ityp.EQ.19) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3x,2I3,8x,A,f6.2,1x,A,f5.2,1x,A,f5.2,2x,A,4F6.2)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),&
               'H:',cvpar(1,i),&
               'O:',cvpar(2,i),&
               'F_c: ',cvpar(3,i),&
               'L:',cvpar(4,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
       ELSEIF (ityp.EQ.27) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3x,3I3,2x,6(A,f6.2,1x),2x,3F6.2)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),atcvar(11,i),&
               'H:',cvpar(1,i),&
               'O:',cvpar(2,i),&
               'nH0: ',cvpar(3,i),&
               'nO0: ',cvpar(4,i),&
               'nL0: ',cvpar(5,i),&
               'L:',cvpar(6,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
       ELSEIF (ityp.EQ.20) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3x,3I3,5x,A,f6.2,1x,A,f5.2,1x,A,f5.2,2x,A,4F6.2)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),atcvar(3,i),&
               'H:',cvpar(1,i),&
               'O:',cvpar(2,i),&
               'F_c: ',cvpar(3,i),&
               'L:',cvpar(4,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)
       ELSEIF (ityp.EQ.28) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3x,4I3,5x,A,f6.2,1x,A,f5.2,1x,A,f5.2,2x,A,4F6.2)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),&
               atcvar(3,i),atcvar(10,i),&
               'H:',cvpar(1,i),&
               'O:',cvpar(2,i),&
               'nL0: ',cvpar(3,i),&
               'L:',cvpar(4,i),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i)

       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,4I5,23x,3F8.4,5x,4F8.4)') styp(ityp),&
               nat_grm(atcvar(1,i)),nat_grm(atcvar(2,i)),&
               nat_grm(atcvar(3,i)),nat_grm(atcvar(4,i)),&
               cscl_fac(1,i),cscl_fac(2,i),cscl_fac(3,i),&
               (vbound(k,i),k=1,4)
       ENDIF
       IF (ityp.EQ.2 .OR. ityp.EQ.3 .OR. ityp.EQ.5) THEN
          CALL raddeg(vbound(1,i),-1)
          CALL raddeg(vbound(3,i),-1)
       ENDIF
    ENDDO

    ! stop 'input'
    IF (paral%io_parent)&
         WRITE(6,'(80("*"))')
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,status)
    RETURN
  END SUBROUTINE meta_colvar_inp
  ! ==================================================================
  SUBROUTINE settological( a, n, l )
    LOGICAL                                  :: a( * )
    INTEGER                                  :: n
    LOGICAL                                  :: l

    INTEGER                                  :: i

    DO i=1,n
       a(i)=l
    ENDDO
  END SUBROUTINE settological
  ! ==================================================================
  SUBROUTINE colvar_structure(tau0,tscr)
    ! ==--------------------------------------------------------------==
    ! == Determination of Collective Variable Structures              ==
    ! ==--------------------------------------------------------------==
    ! 
    ! DOF_INCV   : assign 1 to each dof which partecipate to the CV
    ! dim (NODIM,NCOLVAR,12)
    ! because if the dof is used for a certain CV ICV then 
    ! it has an index among the ions used for the same CV
    ! CV_IST     : at each ions dynamic step the CV value is calculated
    ! dim (NCOLVAR)
    ! DET_COLVAR : at each dynmic step the derivatives of the CV 
    ! wrt to each DOF are calculated. 
    ! dim (NODIM,NCOLVAR)
    ! DET_CELVAR : if TCVCELL some CV are functions of the CELL rather than R
    ! at each dynmic step the derivatives of the CV 
    ! wrt to the cell parameters are calculated
    ! dim (3,3,NCOLVAR)

    ! CV_STORE   : storage of the last INTER_AV_CV istantaneous CV values
    ! in order to be able do calculate the average
    ! CV_PATH    : storage of the average values which give the path in
    ! the CV space 
    ! 
    ! ==--------------------------------------------------------------==    


    REAL(real_8)                             :: tau0(:,:,:), &
                                                tscr(3*maxsys%nax*maxsys%nsx)

    CHARACTER(*), PARAMETER :: procedureN = 'colvar_structure'

    INTEGER                                  :: ia, ib, ic, icv, id, ie, &
                                                ierr, ig, ityp, my_nat
    LOGICAL                                  :: status

    CALL dum2(tau0,tscr)
    IF (ncolvar .EQ. 0) RETURN

    ! get total number of atoms and set indexing.
    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
    ELSE
       my_nat=mmdim%natm
    ENDIF

    ALLOCATE(dof_incv(cotc0%nodim,ncolvar,18),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dof_incv)!,18*ncolvar*cotc0%nodim)
    ALLOCATE(cv_ist(ncolvar*mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cv_ist)!,ncolvar*mwi%nwalk)
    ALLOCATE(det_colvar(cotc0%nodim,ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(det_colvar)!,ncolvar*cotc0%nodim)

    IF (lmeta%tcvcell) THEN
       ALLOCATE(det_celvar(3,3,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(det_celvar)!,ncolvar*9)
    ENDIF

    ! Assign the DPF that take part to the defined CV
    DO icv = 1,ncolvar
       ityp = tycvar(icv)
       IF (ityp .NE. 6 .AND. ityp .LT. 8 ) THEN
          ia = atcvar(1,icv)
          ib = atcvar(2,icv)
          ic = atcvar(3,icv)
          id = atcvar(4,icv)

          CALL get_dof(ia,icv,1,my_nat)
          CALL get_dof(ib,icv,2,my_nat)
          CALL get_dof(ic,icv,3,my_nat)
          CALL get_dof(id,icv,4,my_nat)
       ELSEIF (ityp .EQ. 15) THEN
          ia = atcvar(1,icv)
          ib = atcvar(2,icv)
          ic = atcvar(3,icv)
          id = atcvar(4,icv)
          ie = atcvar(5,icv)
          ig = atcvar(6,icv)

          CALL get_dof(ia,icv,1,my_nat)
          CALL get_dof(ib,icv,2,my_nat)
          CALL get_dof(ic,icv,3,my_nat)
          CALL get_dof(id,icv,4,my_nat)
          CALL get_dof(ie,icv,5,my_nat)
          CALL get_dof(ig,icv,6,my_nat)

       ENDIF
    ENDDO


    ! Calculate starting values of CV and derivatives wrt ions coordinates

    CALL colvarofr(tau0,tscr)
    CALL colvarpr

    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,status)
    RETURN
  END SUBROUTINE colvar_structure

  ! ==================================================================
  SUBROUTINE get_dof(ia,icv,n,my_nat)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, icv, n, my_nat

    INTEGER                                  :: ib, ityp, j, jj, jx, jy, jz, &
                                                kx, naa, nn
    REAL(real_8)                             :: ff

    IF (ia.EQ.0) RETURN
    nn=(n-1)*3
    IF (ia.LT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'GET_DOF! IA=',ia,&
            ' COL.VAR. # =',icv,' N=',n
       CALL stopgm('GET_DOF','IA IS LESS THAN ZERO',& 
            __LINE__,__FILE__)
    ELSEIF (ia.LE.my_nat) THEN
       jx=lskptr(1,ia)
       jy=lskptr(2,ia)
       jz=lskptr(3,ia)
       IF (jx.NE.0) dof_incv(jx,icv,nn+1)=1.0_real_8
       IF (jy.NE.0) dof_incv(jy,icv,nn+2)=1.0_real_8
       IF (jz.NE.0) dof_incv(jz,icv,nn+3)=1.0_real_8
    ELSEIF (ia.LE.my_nat+duat%ndat) THEN
       ! WRITE(6,'(1X,A,I5,A)') 'GET_DOF! IA=',IA,
       ! &                'COL. VAR. WITH DUMMY ATOMS NOT IMPLEMENTED YET'
       naa=ia-my_nat
       ityp=duat%listda(naa,1)
       IF (ityp.EQ.1) THEN
          RETURN
       ELSEIF (ityp.EQ.2) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd2(1,kx)
          ff=1.0_real_8/REAL(jj,kind=real_8)
          DO j=1,jj
             ib=duat%listd2(j+1,kx)
             jx=lskptr(1,ib)
             jy=lskptr(2,ib)
             jz=lskptr(3,ib)
             IF (jx.NE.0) dof_incv(jx,icv,nn+1)=ff
             IF (jy.NE.0) dof_incv(jy,icv,nn+2)=ff
             IF (jz.NE.0) dof_incv(jz,icv,nn+3)=ff
          ENDDO
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'GET_DOF! IA=',ia,&
            ' COL.VAR. # =',icv,' N=',n
       CALL stopgm('GET_DOF','INCORRECT VALUE FOR IA',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_dof
  ! ==================================================================
  SUBROUTINE colvarofr(tau0,tscr)
    ! ==--------------------------------------------------------------==
    ! == Calculates istantaneous values of the collective variables   ==
    ! == and the derivatives wrt the ions coordinates                 ==
    ! ==--------------------------------------------------------------==


    REAL(real_8)                             :: tau0(:,:,:), &
                                                tscr(3,maxsys%nax*maxsys%nsx)

    INTEGER                                  :: iat(15), icv, idof, ipatma, &
                                                ipatmb, isub, ityp, k, kmax, &
                                                mexp(3), nexp(3), numH, numl, &
                                                numO
    REAL(real_8) :: aa, bb, c_km, c_rc1, c_rc2, c_rcm, diff, dx(18), lambda, &
      r0_shift, rc_hy(3), sign, x1(3), x2(3), x3(3), x4(3), x5(3), x6(3)

    CALL tiset(' COLVAROFR',isub)

    CALL dum2(tau0,tscr)
    CALL zeroing(det_colvar)!,ncolvar*cotc0%nodim)
    IF (lmeta%tcvcell) THEN
       CALL zeroing(det_celvar)!,ncolvar*9)
    ENDIF

    ipatmb=0
    ipatma=0
    DO icv=1,ncolvar
       diff = 0.0_real_8
       ityp=tycvar(icv)
       iat(1)=atcvar(1,icv)
       iat(2)=atcvar(2,icv)
       iat(3)=atcvar(3,icv)
       iat(4)=atcvar(4,icv)
       iat(5)=atcvar(5,icv)
       iat(6)=atcvar(6,icv)
       iat(7)=atcvar(7,icv)
       iat(8)=atcvar(8,icv)
       iat(9)=atcvar(9,icv)
       iat(10)=atcvar(10,icv)
       iat(11)=atcvar(11,icv)
       IF (ityp.EQ.1) THEN
          ! ..stretch
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL funcr(diff,cv_ist(icv),0.0_real_8,x1,x2)
          CALL diffr(dx,x1,x2)
          kmax=6
          ! ..angle
       ELSEIF (ityp.EQ.2) THEN
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL fillc(iat(3),tscr,x3)
          CALL funct(diff,cv_ist(icv),0.0_real_8,x1,x2,x3)
          CALL difft(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.3) THEN
          ! ..dihedral angle
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL fillc(iat(3),tscr,x3)
          CALL fillc(iat(4),tscr,x4)
          CALL funco(diff,cv_ist(icv),0.0_real_8,x1,x2,x3,x4,sign)
          CALL diffo(dx,x1,x2,x3,x4,sign)
          kmax=12
       ELSEIF (ityp.EQ.4) THEN
          ! ..distance
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL funcd(diff,cv_ist(icv),0.0_real_8,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
          ! ..displacement in a coordinate
       ELSEIF (ityp.EQ.30) THEN
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          DO k=1,3
             IF (k.NE.iat(3))THEN
                x1(k)=0.0_real_8
                x2(k)=0.0_real_8
             ENDIF
          ENDDO
          CALL funcd(diff,cv_ist(icv),0.0_real_8,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
       ELSEIF (ityp.EQ.5) THEN
          ! ..out of plane
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL fillc(iat(3),tscr,x3)
          CALL fillc(iat(4),tscr,x4)
          CALL funcp_mia(diff,cv_ist(icv),0.0_real_8,x1,x2,x3,x4,dx)
          ! CALL DIFFP(DX,X1,X2,X3,X4)
          kmax=12
       ELSEIF (ityp.EQ.6) THEN
          ! ..coordination number
          c_km     = cvpar(1,icv)
          c_rcm    = cvpar(2,icv)
          CALL coornum(iat(1),c_km,c_rcm,tscr,ions1%nat,det_colvar(1,icv),&
               lskptr,diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.7) THEN
          ! ..difference between distances- 3 atoms
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL fillc(iat(3),tscr,x3)
          CALL funcdd(diff,cv_ist(icv),0.0_real_8,x1,x2,x3)
          CALL diffdd(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.8) THEN
          ! ..species dependent coordination number with Fermi function
          c_km     = cvpar(1,icv)
          c_rcm    = cvpar(2,icv)
          CALL coornumsp(iat(1),iat(2),c_km,c_rcm,tscr,&
               det_colvar(1,icv),lskptr,diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.29) THEN
          ! ..  Total Coordination between a specified group
          c_rcm=cvpar(1,icv)! kappa
          CALL coornumgrp(atcvar(1,icv),max_natcngb,&
               natcngb(ipatma+1),iatcnga(ipatma+1),&
               iatcngb(ipatmb+1),&
               c_rcm,rccnga(ipatma+1),tscr,det_colvar(1,icv),lskptr,&
               diff,cv_ist(icv),0.0_real_8)
          ipatma=ipatma+atcvar(1,icv)
          ipatmb=ipatmb+atcvar(1,icv)*max_natcngb
          kmax=0
       ELSEIF (ityp.EQ.9) THEN
          ! ..species dependent coordination number with rational function
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          IF (specindex(icv).EQ.1) THEN
             CALL coorn_rf(iat(1),iat(2),iat(3),iat(4),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,1,iatdlmn)
          ELSE IF (specindex(icv) .LT. 0) THEN
             CALL  coorn_rf_q(iat,c_rcm,r0_shift,tscr,&
                  det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,specindex(icv))
          ELSE
             CALL coorn_rf(iat(1),iat(2),iat(3),iat(4),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,0,iatdlmn((icv-1)*ions1%nat+1))
          ENDIF
          kmax=0
       ELSEIF (ityp.EQ.10) THEN
          ! ..bond switch with rational function
          c_rcm    = cvpar(1,icv)
          CALL bndswitch(iat(1),iat(2),iat(3),iat(4),c_rcm,tscr,&
               det_colvar(1,icv),diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.11) THEN
          ! ..total species dependent coordination number with rational function
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          IF (specindex(icv).EQ.1) THEN
             CALL coorntot_rf(iat(1),iat(2),iat(3),iat(4),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,1,iatdlmn)
          ELSE IF (specindex(icv) .LT. 0) THEN
             CALL  coorntot_rf_q(iat,c_rcm,r0_shift,tscr,&
                  det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,specindex(icv))
          ELSE
             CALL coorntot_rf(iat(1),iat(2),iat(3),iat(4),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,cv_ist(icv),&
                  0.0_real_8,0,iatdlmn((icv-1)*ions1%nat+1))
          ENDIF
          kmax=0
       ELSEIF (ityp.EQ.35) THEN
          ! ..total species dependent difference of coordination numbers with rational function
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          CALL difcntot_rf(iat(1),iat(2),iat(3),iat(4),iat(5),c_rcm,&
               r0_shift,tscr,det_colvar(1,icv),diff,cv_ist(icv),&
               0.0_real_8,1,iatdlmn)
          kmax=0
       ELSEIF (ityp.EQ.36) THEN
          ! ..total species dependent absolute value of difference of coordination numbers with rational function
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          CALL abs_difcntot_rf(iat(1),iat(2),iat(3),iat(4),iat(5),c_rcm,&
               r0_shift,tscr,det_colvar(1,icv),diff,cv_ist(icv),&
               0.0_real_8,1,iatdlmn)
          kmax=0
       ELSEIF (ityp.EQ.25) THEN
          ! ...coordination numbers calculated between sequences of atoms
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          CALL coornseq_rf(iat,c_rcm,r0_shift,tscr,&
               det_colvar(1,icv),&
               diff,cv_ist(icv),0.0_real_8,specindex(icv))
          kmax = 0
       ELSEIF (ityp.EQ.12) THEN
          ! ..RMSD wrt 2 configurations A and B
          CALL rmsd_ab(iat(1),atcvar(2,icv),atcvar(iat(1)+2,icv),&
               tscr,ions1%nat,det_colvar(1,icv),&
               lskptr,diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.26) THEN
          ! ..RMSD wrt initial configuration
          CALL rmsd_seq_a(iat,maxrmsdat,tscr,det_colvar(1,icv),&
               diff,cv_ist(icv),0.0_real_8,specindex(icv))
          kmax=0
       ELSEIF (ityp.EQ.13) THEN
          ! ..Displacement in the direction (lmn) of a species wrt configuration A 

          IF (specindex(icv).EQ.1) THEN
             CALL disp_lnm(iat(1),iat(2),atcvar(3,icv),&
                  atcvar(3+iat(1)+iat(2),icv),&
                  iatdlmn,tscr,ions1%nat,det_colvar(1,icv),lskptr,&
                  cv_ist(icv),specindex(icv))
             kmax=0
          ELSE

             CALL disp_lnm(iat(1),iat(2),atcvar,atcvar(3,icv),&
                  iatdlmn((icv-1)*ions1%nat+1),tscr,ions1%nat,det_colvar(1,icv),&
                  lskptr, cv_ist(icv),specindex(icv))
             kmax=0
          ENDIF
       ELSEIF (ityp.EQ.15) THEN
          ! ..Angle between 2 planes, each defined by 3 pnts  
          CALL fillc(iat(1),tscr,x1)
          CALL fillc(iat(2),tscr,x2)
          CALL fillc(iat(3),tscr,x3)
          CALL fillc(iat(4),tscr,x4)
          CALL fillc(iat(5),tscr,x5)
          CALL fillc(iat(6),tscr,x6)
          CALL aplane(diff,cv_ist(icv),0.0_real_8,x1,x2,x3,x4,x5,x6,dx)
          kmax=18
       ELSEIF (ityp.EQ.16) THEN
          ! ..Hydrogen Bond Chain

          CALL chain_dr(iat(1),iat(2),iat(3),iat(4),&
               iat(5),iat(6),cvpar(1,icv),&
               iat(7),iat(8),cvpar(2,icv),&
               iatdlmn((icv-1)*ions1%nat+1),tscr,det_colvar(1,icv),&
               lskptr,cv_ist(icv))
       ELSEIF (ityp.EQ.17) THEN
          ! ..Difference in Coordination numbers (rational function)
          c_rcm    = cvpar(1,icv)
          r0_shift   = cvpar(2,icv)
          IF (specindex(icv).EQ.1) THEN
             CALL coorn_rf(iat(2),iat(3),iat(4),iat(5),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,bb,0.0_real_8,1,iatdlmn)
             ! write(6,*) DET_COLVAR
             CALL dscal(cotc0%nodim,-1.0_real_8,det_colvar(1,icv),1)
             CALL coorn_rf(iat(1),iat(3),iat(4),iat(5),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,aa,0.0_real_8,1,iatdlmn)
             ! write(6,*) DET_COLVAR
             cv_ist(icv) = aa - bb
          ELSE
             CALL coorn_rf(iat(2),iat(3),iat(4),iat(5),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,bb,&
                  0.0_real_8,0,iatdlmn((icv-1)*ions1%nat+1))
             CALL dscal(cotc0%nodim,-1.0_real_8,det_colvar(1,icv),1)
             CALL coorn_rf(iat(1),iat(3),iat(4),iat(5),c_rcm,r0_shift,&
                  tscr,det_colvar(1,icv),diff,aa,&
                  0.0_real_8,0,iatdlmn((icv-1)*ions1%nat+1))
             cv_ist(icv) = aa - bb
          ENDIF
          kmax=0
       ELSEIF (ityp.EQ.18) THEN
          ! .. Coordination numbers of second neighbors (rational function)
          c_rc1    = cvpar(1,icv)
          c_rc2    = cvpar(2,icv)
          CALL ctot_chain(iat(1),iat(2),iat(3),iat(4),iat(5),&
               c_rc1,c_rc2,tscr,det_colvar(1,icv),lskptr,&
               diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.31) THEN
          ! ...DIPOLERHO
          CALL cdipolerho(iat(1),tscr,diff,cv_ist(icv),&
               det_colvar(1,icv),iat(2),iatdlmn((icv-1)*ions1%nat+1),&
               iqatdlmn((icv-1)*ions1%nat+1),0.0_real_8)
          kmax = 0
       ELSEIF (ityp.EQ.32) THEN
          ! ...  DIPOLETHETA
          CALL cdipoletheta(iat(1),tscr,diff,cv_ist(icv),&
               det_colvar(1,icv),iat(2),iatdlmn((icv-1)*ions1%nat+1),&
               iqatdlmn((icv-1)*ions1%nat+1),0.0_real_8)
          kmax = 0
          ! _DIPOLEPHI
       ELSEIF (ityp.EQ.33) THEN
          CALL cdipolephi(iat(1),tscr,diff,cv_ist(icv),&
               det_colvar(1,icv),iat(2),iatdlmn((icv-1)*ions1%nat+1),&
               iqatdlmn((icv-1)*ions1%nat+1),0.0_real_8)
          kmax = 0

       ELSEIF (ityp.EQ.34) THEN
          ! _VCOORS
          c_km  = cvpar(1,icv)
          c_rcm = cvpar(2,icv)
          x1(1) = cvpar(3,icv)
          x1(2) = cvpar(4,icv)
          x1(3) = cvpar(5,icv)
          CALL vcoors(iat(1),x1,c_km,c_rcm,tscr,&
               det_colvar(1,icv),lskptr,diff,cv_ist(icv),0.0_real_8)
          kmax=0
       ELSEIF (ityp.EQ.19) THEN
          ! .. Presence of Hydronium

          nexp(1) = iat(3)
          mexp(1) = iat(4)
          rc_hy(1)   = cvpar(1,icv)
          nexp(2) = iat(5)
          mexp(2) = iat(6)
          rc_hy(2)   = cvpar(2,icv)
          nexp(3) = iat(7)
          mexp(3) = iat(8)
          rc_hy(3)   = cvpar(3,icv)
          lambda  = cvpar(4,icv)
          CALL hyd_presence(iat(1),iat(2),nexp,mexp,rc_hy,lambda,&
               tscr,det_colvar(1,icv),lskptr,cv_ist(icv))
          kmax=0
       ELSEIF (ityp.EQ.27) THEN
          ! .. Presence of Hydronium B (with cutoff)

          numO = ions0%na(iat(1))
          numH = ions0%na(iat(2))
          numl = iat(11)
          CALL hyd_cutoff(iat,cvpar(1,icv),numO,numH,numl,&
               iatdlmn((icv-1)*ions1%nat+1),&
               tscr,det_colvar(1,icv),lskptr,cv_ist(icv))
          kmax=0

       ELSEIF (ityp.EQ.20) THEN
          ! .. Presence of Hydronium, distance between ion and hydronium
          nexp(1) = iat(4)
          mexp(1) = iat(5)
          rc_hy(1)   = cvpar(1,icv)
          nexp(2) = iat(6)
          mexp(2) = iat(7)
          rc_hy(2)   = cvpar(2,icv)
          nexp(3) = iat(8)
          mexp(3) = iat(9)
          rc_hy(3)   = cvpar(3,icv)
          lambda  = cvpar(4,icv)
          CALL hyd_ion_distance(iat(1),iat(2),iat(3),nexp,mexp,rc_hy,&
               lambda,tscr,det_colvar(1,icv),lskptr,cv_ist(icv))
          kmax=0

       ELSEIF (ityp.EQ.28) THEN
          ! .. Distance Hydronium - ion with cutoff

          numO = ions0%na(iat(1))
          numH = ions0%na(iat(2))
          numl = iat(10)
          CALL cut_hyd_ion_distance(iat,cvpar(1,icv),numO,numH,numl,&
               iatdlmn((icv-1)*ions1%nat+1),&
               tscr,det_colvar(1,icv),lskptr,cv_ist(icv))
          kmax=0

       ELSEIF (ityp.EQ.21) THEN
          ! .. Localization of the spin density 
          ! ..  independent on the ionic positions

          CALL fillc(iat(1),tscr,x1)
          ! RCC0 is the position of the reference atom 
          rcc0(1,icv)=x1(1)
          rcc0(2,icv)=x1(2)
          rcc0(3,icv)=x1(3)
          ! RCC  is the position of the center of the polarized density
          x2(1)= rcc(1)
          x2(2)= rcc(2)
          x2(3)= rcc(3)

          CALL funcl(diff,cv_ist(icv),0.0_real_8,x1,x2,dx)

          kmax=6
       ELSEIF (ityp.EQ.22) THEN
          ! ..  Volume of the Cell (only if cntl%tprcp)
          CALL volume_cv(cv_ist(icv),det_celvar(1,1,icv))
          kmax=0

       ELSEIF (ityp.EQ.23) THEN
          ! ..  Side of the Cell (only if cntl%tprcp)
          CALL side_cv(cv_ist(icv),iat(1),det_celvar(1,1,icv))
          kmax=0

       ELSEIF (ityp.EQ.24) THEN
          ! ..  ANGLE of the Cell (only if cntl%tprcp)
          CALL angle_cv(cv_ist(icv),iat(1),det_celvar(1,1,icv))
          kmax=0

       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) ' COLVAROFR! I=',icv,' TYPE=',ityP
          CALL stopgm('COLVAROFR','UNKNOWN TYPE OF COLLECTIVE VARIABLE',& 
               __LINE__,__FILE__)
       ENDIF

       DO k=1,kmax
          DO idof=1,cotc0%nodim
             det_colvar(idof,icv)=det_colvar(idof,icv)+&
                  dx(k)*dof_incv(idof,icv,k)
          ENDDO
       ENDDO

    ENDDO  ! NCOLVAR
    CALL tihalt(' COLVAROFR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE colvarofr
  ! ==================================================================
  SUBROUTINE colvarpr
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(len=15), DIMENSION(36), PARAMETER :: styp = (/'        STRETCH',&
      '           BEND','        TORSION','       DISTANCE','           OUTP',&
      '        COORDIN','      DIFFEREN.','         COORSP','        COOR_RF',&
      '      BN.SWITCH','       TOT C_RF','        RMSD_AB','         DISPL1',&
      '        COORDIS','         PLNANG','        HBONDCH','        DIFCOOR',&
      '      COO_CHAIN','      HYDRONIUM','        DIS_HYD','           SPIN',&
      '         VOLVAR','       CELLSIDE','      CELLANGLE','       COOR_SEQ',&
      '       RMSD_SEQ','        HYD_CUT','      D_HYD_CUT','      COORGROUP',&
      '        DISAXIS','      DIPOLERHO','      DIPOLETHA','      DIPOLEPHI',&
      '         VCOORS','  DIFF TOT C_RF','|DIFF TOT C_RF|'/)

    CHARACTER(len=10)                        :: chnum, chnum2
    CHARACTER(len=100)                       :: lineform
    INTEGER                                  :: i, ia, iaa, ie, iee, ityp, k, &
                                                numspec, numspec2, numtot
    REAL(real_8)                             :: cval

! ==--------------------------------------------------------------==
! ==  PRINT SOME INFO                                             ==
! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(A)')&
         ' <<<<<<<<<   COLLECTIVE VARIABLE INFO   >>>>>>>>>>>'
    IF (paral%io_parent)&
         WRITE(6,'(4x,A,6x,A,5x,A,20x,A)') 'TYPE','ATOMS AND SPECIES',&
         'PARAMETERS','VALUE'
    DO i=1,ncolvar
       ityp=tycvar(i)
       cval = cv_ist(i)
       IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) THEN
          CALL raddeg(cval,1)
       ENDIF
       IF (ityp.EQ.6) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,20x,2F8.4,14x,F8.4)') styp(ityp),atcvar(1,i),&
               cvpar(1,i),cvpar(2,i),&
               cval
       ELSEIF (ityp.EQ.30) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I5,20x,2F8.4,14x,F8.4)') styp(ityp),atcvar(1,i),&
               atcvar(2,i), atcvar(3,i),cvpar(1,i),cvpar(2,i),cval
       ELSEIF (ityp.EQ.8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2I5,15x,2F8.4,16x,F8.4)') styp(ityp),atcvar(1,i),&
               atcvar(2,i),cvpar(1,i),cvpar(2,i),&
               cval
       ELSEIF (ityp.EQ.29) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,2F8.4)') styp(ityp),atcvar(1,i),cvpar(1,i),cval
       ELSEIF (ityp.GE.9 .AND. ityp.LE. 11) THEN
          IF (cvpar(2,i) .NE. 0.0_real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,4x,2I5,2x,A,2I4,2f6.3,17x,F8.4)')&
                  styp(ityp),atcvar(1,i),atcvar(2,i),&
                  '2 SHELL',atcvar(3,i),atcvar(4,i),&
                  cvpar(1,i),cvpar(2,i),cval
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,4x,2I5,9x,2I4,f6.3,17x,F8.4)') styp(ityp),&
                  atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
                  cvpar(1,i),cval
          ENDIF
       ELSEIF (ityp.EQ.17) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4x,3I4,8x,2I3,f6.3,17x,F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cvpar(1,i),cval
       ELSEIF (ityp.EQ.18) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4x,3I4,4x,2I3,2f6.3,17x,F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cvpar(1,i),cval
       ELSEIF (ityp.EQ.12) THEN
          numspec  = atcvar(1,i)

          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') numspec
          CALL xstring(chnum,ia,ie)

          lineform ='(A,2X,I4,A9,'//chnum(ia:ie)//'I4,18x,F8.4)'
          IF (paral%io_parent)&
               WRITE(6,lineform) styp(ityp),atcvar(1,i),' Species ',&
               (atcvar(k+1,i),k=1,numspec),cval
       ELSEIF (ityp.EQ.13 .AND. specindex(i) .EQ. 1) THEN
          numspec = atcvar(1,i)
          numspec2 = atcvar(2,i)
          numtot  = numspec + numspec2
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') numspec
          CALL xstring(chnum,ia,ie)
          IF (numspec2 .NE. 0) THEN
             IF (paral%io_parent)&
                  WRITE(chnum2,'(I5)') numspec2
             CALL xstring(chnum2,iaa,iee)
             lineform ='(A,1X,A11,'//chnum(ia:ie)//'I4,A11,'&
                  //chnum2(iaa:iee)//'I4,A11,3I3,f12.6)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),&
                  ' Species 1 ',(atcvar(k+2,i),k=1,numspec),&
                  ' Species 2 ',(atcvar(k+2+numspec,i),k=1,numspec2),&
                  ' MillerInd:',(atcvar(k+2+numtot,i),k=1,3),&
                  cval
          ELSE
             lineform ='(A,1X,A,'//chnum(ia:ie)//'I4,A10,3I3,f12.6)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),&
                  ' Species 1 ',(atcvar(k+2,i),k=1,numspec),&
                  ' Mil.Ind.',(atcvar(k+2+numtot,i),k=1,3),&
                  cval
          ENDIF
       ELSEIF (ityp.EQ.13 .AND. specindex(i) .NE. 1) THEN
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') atcvar(1,i)
          CALL xstring(chnum,ia,ie)
          numspec2 = atcvar(2,i)
          IF (numspec2 .NE. 0) THEN
             IF (paral%io_parent)&
                  WRITE(chnum2,'(I5)') atcvar(2,i)
             CALL xstring(chnum2,iaa,iee)
             lineform ='(A9,2X,A5,'//chnum(ia:ie)//'I4,A8,'&
                  //chnum2(iaa:iee)//'I4,A10,3I3,16X,f12.6)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),' at.',&
                  (iatdlmn(k+ions1%nat*(i-1)),k=1,atcvar(1,i)),&
                  ' vs at.',(iatdlmn(atcvar(1,i)+k+ions1%nat*(i-1)),&
                  k=1,atcvar(2,i)),&
                  ' Mil. In.',(atcvar(k+2,i),k=1,3), cval
          ELSE
             lineform ='(A,2X,A5,'//chnum(ia:ie)//&
                  'I4,A8,A10,3I3,16X,f12.6)'
             IF (paral%io_parent)&
                  WRITE(6,lineform) styp(ityp),' at.',&
                  (iatdlmn(k+ions1%nat*(i-1)),k=1,atcvar(1,i)),&
                  ' vs all ',' Mil. In.',(atcvar(k+2,i),k=1,3),&
                  cval
          ENDIF

       ELSEIF (ityp.EQ.15) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,6I5,20x,F10.5)') styp(ityp),&
               (atcvar(k,i),k=1,6),&
               cval


       ELSEIF (ityp.EQ.16) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,4I5,30x,F8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               cval
       ELSEIF (ityp.EQ.19) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,2I4,3x,2I3,f6.2,2I3,f6.2,2I3,2f6.2,F8.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),&
               atcvar(3,i),atcvar(4,i),cvpar(1,i),&
               atcvar(5,i),atcvar(6,i),cvpar(2,i),&
               atcvar(7,i),atcvar(8,i),cvpar(3,i),&
               cvpar(4,i),cval
       ELSEIF (ityp.EQ.27) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I4,3x,2I3,f6.2,2I3,f6.2,2I3,2f6.2,'//&
               '2I3,2f6.2,F8.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),atcvar(11,i),&
               atcvar(3,i),atcvar(4,i),cvpar(1,i),&
               atcvar(5,i),atcvar(6,i),cvpar(2,i),&
               atcvar(7,i),atcvar(8,i),cvpar(3,i),&
               cvpar(4,i),atcvar(9,i),atcvar(10,i),&
               cvpar(5,i),cvpar(6,i),cval
       ELSEIF (ityp.EQ.28) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I4,3x,2I3,f6.2,2I3,f6.2,2I3,f6.2,'//&
               'I3,f6.2,F8.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),atcvar(3,i),&
               atcvar(4,i),atcvar(5,i),cvpar(1,i),&
               atcvar(6,i),atcvar(7,i),cvpar(2,i),&
               atcvar(8,i),atcvar(9,i),cvpar(3,i),&
               atcvar(10,i),cvpar(4,i),cval
       ELSEIF (ityp.EQ.20) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,3I4,1x,2I3,f6.2,2I3,f6.2,2I3,2f6.2,F8.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),atcvar(3,i),&
               atcvar(4,i),atcvar(5,i),cvpar(1,i),&
               atcvar(6,i),atcvar(7,i),cvpar(2,i),&
               atcvar(8,i),atcvar(9,i),cvpar(3,i),&
               cvpar(4,i),cval
       ELSEIF (ityp.EQ.28) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,25x,I4,25x,f6.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),cval


       ELSEIF (ityp.EQ.29) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,25x,I4,25x,f6.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),cval

       ELSEIF (ityp.EQ.30) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,25x,I4,25x,f6.4)')&
               styp(ityp),atcvar(1,i),atcvar(2,i),cval

       ELSEIF (ityp.EQ.31) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,22x,f6.4,2x,f6.4,18x,f6.4)')&
               styp(ityp),atcvar(1,i),cvpar(1,i),cvpar(2,i),&
               cval
       ELSEIF (ityp.EQ.35) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,5I6,25x,F10.5)') styp(ityp),atcvar(1,i),&
               atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cval
       ELSEIF (ityp.EQ.36) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,5I6,25x,F10.5)') styp(ityp),atcvar(1,i),&
               atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cval
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,4I5,25x,F10.5)') styp(ityp),atcvar(1,i),&
               atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               cval
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A,E14.8)') ' ChkSum(COLVAR) = ',&
         SUM(ABS(cv_ist(1:ncolvar)))
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE colvarpr
  ! ==================================================================
  SUBROUTINE param_meta(line,tad_scf,vbound,ibound,&
       cv_dyn_0,initial_value,&
       ncolvar,cscl_fac,kharm,cv_mass,lextlagrange)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=80)                        :: line
    LOGICAL                                  :: tad_scf(*)
    REAL(real_8)                             :: vbound(4,*)
    INTEGER                                  :: ibound(*)
    REAL(real_8)                             :: cv_dyn_0(*)
    LOGICAL                                  :: initial_value(*)
    INTEGER                                  :: ncolvar
    REAL(real_8)                             :: cscl_fac(3,*), kharm(*), &
                                                cv_mass(*)
    LOGICAL                                  :: lextlagrange

    INTEGER                                  :: ia, ie, ii, ij
    LOGICAL                                  :: erread

! 
! 

    ii=INDEX(line,'SCA')
    ij=INDEX(line,'SCF')
    IF (ii.NE.0) THEN
       ia = ii+3
       tad_scf(ncolvar) = .TRUE.
       CALL readsr(line,ia,ie,cscl_fac(1,ncolvar),   erread)
       ia = ie
       CALL readsr(line,ia,ie,cscl_fac(2,ncolvar),   erread)
       ia = ie
       CALL readsr(line,ia,ie,cscl_fac(3,ncolvar),   erread)
       ia = ie
       IF (erread) GOTO 22
    ELSEIF (ij.NE.0) THEN
       ia = ij+3
       tad_scf(ncolvar) = .FALSE.
       CALL readsr(line,ia,ie,cscl_fac(1,ncolvar),   erread)
       IF (erread) GOTO 22
    ENDIF
    ii=INDEX(line,'KCV')
    IF (ii .NE. 0) THEN
       lextlagrange=.TRUE.
       ia = ii+3
       CALL readsr(line,ia,ie,kharm(ncolvar),   erread)
       IF (erread) GOTO 22
    ENDIF
    ij=INDEX(line,'MCV')
    IF (ij .NE. 0) THEN
       lextlagrange=.TRUE.
       ia = ij+3
       CALL readsr(line,ia,ie,cv_mass(ncolvar),   erread)
       IF (erread) GOTO 22
    ENDIF
    ! ==----------- DEF. OF BOUNDARIES IN CV. SPACE  ---------------==
    ii=INDEX(line,'WALL+')
    IF (ii .NE. 0) THEN
       ibound(ncolvar)=1
       ia = ii+5
       ! position of the wall
       CALL readsr(line,ia,ie,vbound(1,ncolvar),   erread)
       ia = ie
       ! height of the wall (hartrees)
       CALL readsr(line,ia,ie,vbound(2,ncolvar),   erread)
    ENDIF
    ii=INDEX(line,'WALL-')
    IF (ii .NE. 0) THEN
       ibound(ncolvar)=1
       ia = ii+5
       ! position of the wall
       CALL readsr(line,ia,ie,vbound(3,ncolvar),   erread)
       ia = ie
       ! height of the wall (hartrees)
       CALL readsr(line,ia,ie,vbound(4,ncolvar),   erread)
    ENDIF
    ! mb-ale - we want to be free to set the initial value !
    ii=INDEX(line,'INITIAL_VALUE')
    IF (ii .NE. 0) THEN
       initial_value(ncolvar)=.TRUE.
       ia=ii+13
       CALL readsr(line,ia,ie,cv_dyn_0(ncolvar),erread)
    ENDIF

    RETURN
22  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING SCAL. FAC. or MASS or KH.'
    IF (paral%io_parent)&
         WRITE(6,'(A80)') line
    CALL stopgm('M_COLVAR_INP',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE param_meta
  ! ==================================================================


  SUBROUTINE cv_forces3(taup,tscr,velp,cv_f_new)

    REAL(real_8)                             :: taup(:,:,:), tscr(:,:,:), &
                                                velp(:,:,:), cv_f_new(ncolvar)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cv_forces3'
    INTEGER, PARAMETER                       :: maxrat = 10, mrdiis = 5 
    REAL(real_8), PARAMETER                  :: tol_cv_f = 1.e-4_real_8, &
                                                tol_cv_x = 1.e-5_real_8 

    EXTERNAL                                 :: dasum, ddot
    INTEGER                                  :: i, ia, icount, icv, idiis, &
                                                idof, ierr, info, is, iter, &
                                                j, jcount, length, m, nactive
    INTEGER, ALLOCATABLE                     :: iflag(:), ipos(:), ipvt(:)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: dasum, ddot, &
                                                diism(mrdiis+1,mrdiis+1), &
                                                errf, errx, fact, sum, &
                                                v(mrdiis+1)
    REAL(real_8), ALLOCATABLE :: asl(:,:), cv_diff(:), cv_new(:), cv_old(:), &
      cv_str(:), dcv_dr_new(:,:), dcv_dr_str(:,:), dx(:), err(:,:), &
      r_new(:,:,:), r_nnew(:,:,:), r_old(:,:,:), xlo(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_f_old(:)

! ==--------------------------------------------------------------==
! Memory Allocation at IFIRST=0

    IF (ifirst .EQ. 0) THEN
       ALLOCATE(cv_f_old(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_f_old)!,ncolvar)

       ALLOCATE(dx(3*cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(asl(ncolvar,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipvt(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xlo(ncolvar,mrdiis),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(err(ncolvar,mrdiis),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(iflag(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipos(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst = 1
    ENDIF
    CALL zeroing(asl)!,ncolvar*ncolvar)
    CALL zeroing(xlo)!,ncolvar*mrdiis)
    CALL zeroing(err)!,ncolvar*mrdiis)
    CALL zeroing(iflag)!,ncolvar)
    CALL zeroing(ipos)!,ncolvar)
    ! ==--------------------------------------------------------------==

    length = 4*ncolvar + 3*3*maxsys%nax*maxsys%nsx + 2*cotc0%nodim*ncolvar

    ! TODO align for BG
    ALLOCATE(cv_str(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(cv_new(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(cv_old(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r_new(3, maxsys%nax, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r_old(3, maxsys%nax, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dcv_dr_str(cotc0%nodim, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dcv_dr_new(cotc0%nodim, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r_nnew(3, maxsys%nax, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(cv_diff(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ! Store Values
    CALL dcopy(ncolvar,cv_ist,1,cv_str,1)
    CALL dcopy(ncolvar*cotc0%nodim,det_colvar,1,dcv_dr_str,1)

    ! Old Ions Positions aCV and DET_CV
    !$omp parallel do private(IS,IA)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          r_old(1,ia,is) = taup(1,ia,is) - velp(1,ia,is)*dt_ions
          r_old(2,ia,is) = taup(2,ia,is) - velp(2,ia,is)*dt_ions
          r_old(3,ia,is) = taup(3,ia,is) - velp(3,ia,is)*dt_ions
       ENDDO
    ENDDO
    CALL colvarofr(r_old,tscr)
    CALL dcopy(ncolvar, cv_ist,1,cv_old,1)

    ! New Ions Positions aCV and DET_CV
    !$omp parallel do private(IS,IA,FACT)
    DO is=1,ions1%nsp
       fact=-dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          r_new(1,ia,is) =  r_old(1,ia,is)+&
               (velp(1,ia,is)*dt_ions&
               +fact*fhills(1,ia,is))
          r_new(2,ia,is) =  r_old(2,ia,is)+&
               (velp(2,ia,is)*dt_ions&
               +fact*fhills(2,ia,is))
          r_new(3,ia,is) =  r_old(3,ia,is)+&
               (velp(3,ia,is)*dt_ions&
               +fact*fhills(3,ia,is))
       ENDDO
    ENDDO
    CALL colvarofr(r_new,tscr)

    CALL dcopy(ncolvar, cv_ist,1,cv_new,1)
    CALL dcopy(ncolvar*cotc0%nodim,det_colvar,1,dcv_dr_new,1)

    ! First Guess for the forces
    CALL zeroing(dx)!,cotc0%nodim)
    DO icv=1,ncolvar
       !$omp parallel do private(IDOF)
       DO idof=1,cotc0%nodim
          dx(idof)=dx(idof)+cv_f_old(icv)*dcv_dr_new(idof,icv)
       ENDDO
    ENDDO

    ! Update the positions
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,dx)
    !$omp parallel do private(IS,IA,FACT)
    DO is=1,ions1%nsp
       fact=-dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          r_nnew(1,ia,is)=r_new(1,ia,is)+fact*tscr(1,ia,is)
          r_nnew(2,ia,is)=r_new(2,ia,is)+fact*tscr(2,ia,is)
          r_nnew(3,ia,is)=r_new(3,ia,is)+fact*tscr(3,ia,is)
       ENDDO
    ENDDO


    CALL dcopy(ncolvar,cv_f_old,1,cv_f_new,1)


    icount = 0

    DO icv = 1,ncolvar
       sum = 0.0_real_8
       DO idof = 1,cotc0%nodim
          sum = sum + det_colvar(idof,icv)*det_colvar(idof,icv)
       ENDDO
       IF (sum .LT. 1.e-2_real_8) THEN
          iflag(icv) = 0
          cv_f_new(icv) = 0.0_real_8
       ELSE
          iflag(icv) = 1
          icount = icount + 1
          ipos(icount) = icv
       ENDIF
    ENDDO
    nactive = icount

    ! Iterativ calculation of lambda
    DO iter=1,maxrat
       ! Calculate Collective Variables values and differences
       ! for the current value of lambda
       CALL zeroing(dx)!,cotc0%nodim)
       CALL colvarofr(r_nnew,tscr)

       DO icount = 1,nactive
          cv_diff(icount) = -(cv_ist(ipos(icount))-cv_old(ipos(icount)))
       ENDDO
       errf=dasum(nactive,cv_diff(1),1)

       IF (errf .LT. tol_cv_f) GOTO 100

       ! Derivatives of sigma wrt lambda
       DO icount=1,nactive
          DO jcount=1,nactive
             asl(icount,jcount)=0.0_real_8
             DO idof=1,cotc0%nodim
                asl(icount,jcount)=asl(icount,jcount)-&
                     dt_ions*dtm(idof)*det_colvar(idof,ipos(icount))&
                     *dcv_dr_new(idof,ipos(jcount))
             ENDDO
          ENDDO
       ENDDO

       ! Solve for asl*cg=fc
       ! LAPACK matrix solver
       CALL dgesv(nactive,1,asl,ncolvar,ipvt,cv_diff,ncolvar,info)
       IF (info.NE.0) CALL stopgm(' CV_FORCES|','ERROR IN DGESV',& 
            __LINE__,__FILE__)
       errx=dasum(nactive,cv_diff(1),1)
       ! cntl%diis!
       idiis=MOD(iter-1,mrdiis)+1
#ifndef __SR11000
       !$omp parallel do private(ICOUNT)
#endif
       DO icount = 1,nactive
          xlo(icount,idiis) = cv_f_new(ipos(icount))
       ENDDO
       ! CALL DCOPY(NACTIVE,CV_F_NEW(1),1,XLO(1,IDIIS),1)
       CALL dcopy(nactive,cv_diff(1),1,err(1,idiis),1)
       IF (iter.GT.mrdiis) THEN
          m=mrdiis+1
          CALL zeroing(diism)!,m*m)
          CALL zeroing(v)!,m)
          DO i=1,mrdiis
             DO j=1,mrdiis
                diism(i,j)=ddot(nactive,err(1,i),1,err(1,j),1)
             ENDDO
             diism(m,i)=1.0_real_8
             diism(i,m)=1.0_real_8
          ENDDO
          v(m)=1.0_real_8
          CALL solve(diism,m,m,v)
          CALL zeroing(cv_diff)!,nactive)
          CALL zeroing(cv_f_new)!,ncolvar)
          DO i=1,mrdiis
             DO jcount=1,nactive
                cv_diff(jcount)=cv_diff(jcount)+v(i)*err(jcount,i)
                cv_f_new(ipos(jcount))=cv_f_new(ipos(jcount))+&
                     v(i)*xlo(jcount,i)
             ENDDO
          ENDDO
       ENDIF
       DO  icount = 1,nactive
          cv_f_new(ipos(icount)) = cv_f_new(ipos(icount))+&
               cv_diff(icount)
       ENDDO

       IF (errx .LT. tol_cv_x) GOTO 100

       ! Update forces
       CALL zeroing(dx)!,cotc0%nodim)
       DO icount=1,nactive
          fact=cv_f_new(ipos(icount))
          DO idof=1,cotc0%nodim
             ! !es        DX(IDOF)=DX(IDOF)+CV_F_NEW(IPOS(ICOUNT))*
             ! !es &                        DCV_DR_NEW(IDOF,IPOS(ICOUNT))
             dx(idof)=dx(idof)+fact*dcv_dr_new(idof,ipos(icount))
          ENDDO
       ENDDO
       ! Update the positions
       CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
       CALL gettau(tscr,dx)
#ifndef __SR11000
       !$omp parallel do private(IS,IA,FACT)
#endif
       DO is=1,ions1%nsp
          fact=-dt_ions*dtb2mi(is)
          DO ia=1,ions0%na(is)
             r_nnew(1,ia,is)=r_new(1,ia,is)+fact*tscr(1,ia,is)
             r_nnew(2,ia,is)=r_new(2,ia,is)+fact*tscr(2,ia,is)
             r_nnew(3,ia,is)=r_new(3,ia,is)+fact*tscr(3,ia,is)
          ENDDO
       ENDDO

    ENDDO

    IF (paral%io_parent)&
         WRITE(6,'(/,A,I7,A)') ' CV_FORCES|  did not converge '
    DO icount=1,nactive
       IF (ABS(cv_diff(icount)) .GT. tol_cv_x ) THEN
          cv_f_new(ipos(icount)) = cv_f_old(ipos(icount))/2.0_real_8
       ENDIF
    ENDDO
100 CONTINUE

    CALL dcopy(ncolvar,cv_f_new,1,cv_f_old,1)

    ! Recall CV and CV_DET Values
    CALL dcopy(ncolvar,cv_str,1,cv_ist,1)
    CALL dcopy(ncolvar*cotc0%nodim,dcv_dr_str,1,det_colvar,1)

    ! ==--------------------------------------------------------------==
    DEALLOCATE(cv_str,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(cv_new,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(cv_old,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r_new,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r_old,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dcv_dr_str,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dcv_dr_new,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r_nnew,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(cv_diff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cv_forces3

END MODULE meta_colvar_inp_utils
