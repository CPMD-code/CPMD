c     ==================================================================
      SUBROUTINE mm_read_qmmm_input
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 29 January 2008 after J. Blumberger's request
      use mp_interface, only: mp_bcast
      use system
      use parac
      use mm_input
      use mm_dimmod
      use efld
      use forcematch
      use readsr_utils, only : readsr, xstring
      use inscan_utils, only : inscan
      use coordsz

      IMPLICIT NONE

      include 'Gromos/vector.h'
      !include 'Gromos/coordsz.h'
      include 'Gromos/forcesz.h'
      include 'Gromos/toposz.h'

      integer IUNIT,IERR
      logical ft,fi,fc
      integer ic,iw,i,MSGLEN,ia,ie
      CHARACTER LINE*80
      
C     ==--------------------------------------------------------------==
      lqmmm%qmmm=.false.
      lqmmm%qmmm_time=.false.
      lqmmm%qmmm_verbose=.false. 
      lqmmm%qmmm_reflex=.false.  ! cmb-reflecting walls
      capp%cap_h = .false.
      IF(paral%io_parent) THEN  
        IUNIT = 5
        IERR=INSCAN(IUNIT,'&QMMM')
        if(IERR.eq.0)then 
          lqmmm%qmmm=.true.
          lqmmm%qmmm_dcd=.false.
          mm_stat=.true.
          TEXTFLD=.false.
          wp_l%write_potential=.false.
          wp_l%write_density=.false.
          wp_i%nfi_wp=1
          wp_i%nfi_wd=1
          wp_i%nfi_lt=-1
          wp_c%cubename_dens='DENSITY.cube'
          wp_c%cubename_pot='POTENTIAL.cube'
          iqmmm%nt_nlist=0
          rqmmm%rcut_el=10.D0
          rqmmm%rcut_small=10.D0
          rqmmm%rcut_esp=10.D0
          rqmmm%box_toll=7.0D0 ! assumed size (3.7 Angstrom) to be added to
                         ! x/y/z min/max to get the dimension of the 
                         ! charge distribution (based on methane "diameter").
                         ! changed: AK 2008/03/03. see also mm_center.F 
          rqmmm%box_wall=7.0D0 ! cmb-default reflection at 3.7 Angstrom.
          mne%maxnat_el=6000
          iqmmm%coupl_model=0
          sppph%mm_split=.false. ! no QM/MM processors splitting by default
          ft=.false.
          fi=.false.
          fc=.false.
          n_sm%n_smear=4
C     write 'interacting' group of files by default only for MD runs.
          IF(CNTL%MD) THEN
            iqmmm%nt_sample=10
          ELSE
            iqmmm%nt_sample=0
          ENDIF
          r_esp%esp_weight=0.1D0
          n_sml%no_hirshfeld=.false.
          cgrest_i%n_cg=0
          cgrest_r%lambda=0.d0
          cgrest_r%q_rest=0.d0
          solqmmi%nwq=0
          solqmmm%all_water=.false.
          solqmmm%tflexsol=.false.
          agr%qmmm_amber=.false.
          ragr%scale_14_amber=1.0d0
          wp_p%n_print_wf=0
          wp_l%write_orbs=.false.
          do ic=1,20
            solqmmi%solv_bond(ic)=-3000
          enddo
          wp_i%n_stride=1
          rtr_l%restart_traj=.false.
          rtr_i%iframe_restart=1
          rtr_r%rev_vel=1.d0
          rtr_c%file_traj_restart='TRAJECTORY'
          excl_comm%nce=0
          excl_comm_l%gromos_excl=.false.
          excl_comm_l%excl_mech=.false.
C_FM[
C        set default values for forcematching
          if (CNTL%FMATCH) then
            fm_fmatch = .true.
            fm_capping = .false.
            fm_ncap = 0
            fm_compute_traj = .false.
            fm_compute_sp =   .true.
            fm_read_forces = .false.
            fm_read_covforces = .false.
            fm_equivalence = .false.
            fm_fit_charges = .true.
            fm_charges_only = .false.
            fm_fit_fc_only = .false.
            fm_fit_bonds=.true.
            fm_fit_angles=.true.
            fm_fit_diheds=.true.
            fm_fit_improps=.true.
            fm_rinit_wf = .true.
            fm_sp_restart = .false.
            fm_compute_rms = .true.
            fm_ref_traj_filename = 'TRAJECTORY_REF'
            fm_chj_ref_filename = 'FM_REF_CHJ'
            fm_pip_ref_filename = 'FM_REF_PIP'
            fm_forces_ref_filename = 'FM_REF_FORCES'
            fm_covforces_ref_filename = 'FM_REF_COVFORCES'
            fm_ref_traj_stride=1
            fm_outtopo_filename = 'FMATCH.top'
            fm_optimize_weights = .true.
            fm_wv = 1.d0
            fm_wf = 0.d0
            fm_wtot = 1.0d06
            fm_max_kfit_iterations = 500
            fm_wq(:) = 0.1d0
C        MANUEL fix charges
C        and define weights for the charge fitting
            fm_wq_unique(:) = 0.1d0
            fm_fix_q=.false.
            fm_n_qw=0
            fm_wq_general=0.1d0
            fm_nexclude = 0
          else
            fm_fmatch = .false.
          endif
C TODO: Check if we have to change some QMMM defaults
C_FM]
   22     continue
          READ(IUNIT,END=24,ERR=24,FMT='(A)') LINE
          IF(INDEX(LINE,'&END').NE.0) GOTO 23
          IF(INDEX(LINE,'ARRAYSIZ').NE.0) THEN
            CALL SETARRAYSZ(IUNIT)
C_FM[
          ELSE IF(INDEX(LINE,'FORCEMATCH').NE.0) THEN
            CALL READ_FMATCH_INPUT(IUNIT)
C_FM]
          ELSE IF(INDEX(LINE,'TOPOLOGY').NE.0) THEN
            read(iunit,'(A)')cqmmm%filetopo
            ft=.true.
          ELSE IF(INDEX(LINE,'COORDIN').NE.0) THEN
            read(iunit,'(A)')cqmmm%fileconf
            fc=.true.
          ELSE IF(INDEX(LINE,'INPUT').NE.0) THEN
            read(iunit,'(A)')cqmmm%filegrin
            fi=.true.
          ELSE IF(INDEX(LINE,'TIMING').NE.0) THEN
            lqmmm%qmmm_time=.true.
          ELSE IF(INDEX(LINE,'VERBOSE').NE.0) THEN
            lqmmm%qmmm_verbose=.true.
          ELSE IF((INDEX(LINE,'BOX').NE.0).AND.
     &            (INDEX(LINE,'TOLER').NE.0)) THEN
            read(iunit,*) rqmmm%box_toll
          ELSE IF((INDEX(LINE,'BOX').NE.0).AND.
     &            (INDEX(LINE,'WALLS').NE.0)) THEN
            read(iunit,*) rqmmm%box_wall
            lqmmm%qmmm_reflex=.true.  ! cmb-reflecting walls
            WRITE(6,'(1x,''REFLECTING WALLS APPLIED TO QM BOX:'')')
            WRITE(6,'(1x,''       box_wall ='',f14.6, ''(au)'')')
     .            rqmmm%box_wall
            WRITE(6,'(1x,''WARNING: CENTERING OF QM SYSTEM DISABLED'')')
          ELSE IF(INDEX(LINE,'ELECTROSTATIC').NE.0) THEN
            iqmmm%coupl_model=1
            iqmmm%nt_nlist=100
            TEXTFLD=.true.
            IF(INDEX(LINE,'LONG').NE.0) THEN
              iqmmm%coupl_model=2
            ENDIF 
          ELSE IF(INDEX(LINE,'RCUT_NN').NE.0) THEN
            read(iunit,*)rqmmm%rcut_small
          ELSE IF(INDEX(LINE,'RCUT_MIX').NE.0) THEN
            read(iunit,*)rqmmm%rcut_el
          ELSE IF(INDEX(LINE,'RCUT_ESP').NE.0) THEN
            read(iunit,*)rqmmm%rcut_esp
          ELSE IF(INDEX(LINE,'LIST').NE.0
     &           .and. INDEX(LINE,'UPD').NE.0) THEN
            read(iunit,*)iqmmm%nt_nlist
          ELSE IF(INDEX(LINE,'ESPWEIGHT').NE.0) THEN
            read(iunit,*) r_esp%esp_weight
          ELSE IF(INDEX(LINE,'HIRSHFELD').NE.0
     &           .and. INDEX(LINE,'OFF').NE.0) THEN
            n_sml%no_hirshfeld=.true.
          ELSE IF(INDEX(LINE,'AMBER').NE.0) THEN
            agr%qmmm_amber=.true.
            ragr%scale_14_amber=1.0d0/1.2d0
          ELSE IF(INDEX(LINE,'HARMO').NE.0) THEN
            agr%qmmm_amber=.true.
          ELSE IF(INDEX(LINE,'SCALE14').NE.0) THEN
            read(iunit,*) ragr%scale_14_amber
          ELSE IF(INDEX(LINE,'GROMOS').NE.0
     &            .and. INDEX(LINE,'EXCLUSION').EQ.0) THEN
            agr%qmmm_amber=.false.
            ragr%scale_14_amber=1.0d0
          ELSE IF(INDEX(LINE,'NOSPLIT').NE.0) THEN
            sppph%mm_split=.false.
          ELSE IF(INDEX(LINE,'WRITE').NE.0) THEN
            IF(INDEX(LINE,'POTENT').NE.0) THEN
              wp_l%write_potential=.true.
              IF(INDEX(LINE,'STRIDE').NE.0) then 
                do ic=1,70
                  if(line(ic:ic+5).eq."STRIDE") then
                    read(line(ic+6:80),*)wp_i%n_stride
                    goto 441
                  endif
                enddo
 441            continue
              ENDIF
              IF(INDEX(LINE,'STEP').NE.0) then 
                do ic=1,70
                  if(line(ic:ic+3).eq."STEP") then
                    read(line(ic+4:80),*)wp_i%nfi_wp
                    goto 442
                  endif
                enddo
 442            continue
              ENDIF
            ELSE IF(INDEX(LINE,'DENS').NE.0)THEN
              wp_l%write_density=.true.
              IF(INDEX(LINE,'STRIDE').NE.0) then 
                do ic=1,70
                  if(line(ic:ic+5).eq."STRIDE") then
                    read(line(ic+6:80),*)wp_i%n_stride
                    goto 443
                  endif
                enddo
 443            continue
              ENDIF
              IF(INDEX(LINE,'STEP').NE.0) then 
                do ic=1,70
                  if(line(ic:ic+3).eq."STEP") then
                    read(line(ic+4:80),*)wp_i%nfi_wd
                    goto 444
                  endif
                enddo
 444            continue
              ENDIF
            ELSE IF(INDEX(LINE,'LOCAL').NE.0)THEN
              wp_i%nfi_lt=1
              IF(INDEX(LINE,'STEP').NE.0) then
                do ic=1,70
                  if(line(ic:ic+3).eq."STEP") then
                    read(line(ic+4:80),*)wp_i%nfi_lt
                    if (wp_i%nfi_lt.LE.0) wp_i%nfi_lt=1
                    goto 4441
                  endif
                enddo
 4441           continue
              ENDIF
            ELSE IF(INDEX(LINE,'ORB').NE.0) THEN
              read(iunit,*)wp_p%n_print_wf
              read(iunit,*)(wp_p%i_wf_pr(i),i=1,wp_p%n_print_wf)
              wp_l%write_orbs=.true.
              wp_i%n_stride=2
            ELSE
              WRITE(6,*) 'MM_READ_QMMM_INPUT| UNKNOWN WRITE KEY: ', LINE
            ENDIF
          ELSE IF(INDEX(LINE,'RESTART').NE.0
     &           .and.INDEX(LINE,'TRAJ').NE.0)THEN
            rtr_l%restart_traj=.true.
            IF(INDEX(LINE,'FRAME').NE.0) then
              do ic=1,70
                if(line(ic:ic+4).eq."FRAME") then
                  read(line(ic+5:80),*)rtr_i%iframe_restart
                  goto 445
                endif
              enddo
445           continue
            ENDIF
            IF(INDEX(LINE,'FILE').NE.0) then
              do ic=1,70
                if(line(ic:ic+3).eq."FILE") then
                  read(line(ic+4:80),*)rtr_c%file_traj_restart
                  goto 446
                endif
              enddo
446           continue
            ENDIF
            if(INDEX(LINE,'REVERSE').NE.0)rtr_r%rev_vel=-1.d0
          ELSE IF((INDEX(LINE,'MAXNN').NE.0)  ! MAXNAT is kept for backward compatibility
     &           .OR.(INDEX(LINE,'MAXNAT').NE.0)) THEN
            read(iunit,*)mne%maxnat_el             ! max number of interacting atoms
          ELSE IF((INDEX(LINE,'SAMPLE').NE.0)
     &           .AND.(INDEX(LINE,'INTERACT').NE.0)) THEN
            IF(INDEX(LINE,'DCD').NE.0) THEN
              lqmmm%qmmm_dcd=.true.
            ENDIF
            IF(INDEX(LINE,'OFF').NE.0) THEN
              iqmmm%nt_sample=0
            ELSE
              read(iunit,*)iqmmm%nt_sample
            ENDIF
          ELSE IF((INDEX(LINE,'CAPPING').NE.0)
     &           .OR.(INDEX(LINE,'CAP_HYDROGEN').NE.0)) THEN
            capp%cap_h = .true.
          ELSE IF(INDEX(LINE,'ADD_HYDROGEN').NE.0) THEN
            READ(IUNIT,*)addh%n_carbon
            WRITE(6,*) addh%n_carbon
            READ(IUNIT,*)(addh%ig_carbon(i),i=1,addh%n_carbon)
            WRITE(6,'(1x,''H-capped C atoms:'')')
            WRITE(6,*)(addh%ig_carbon(i),i=1,addh%n_carbon)
C FIXME:  AK 2005/05/06: undocumented (second check for 'CH' was pointless, too)
          ELSE IF(INDEX(LINE,'CHARGE').NE.0) THEN
            read(iunit,*)cgrest_i%n_cg
            read(iunit,*)(cgrest_i%atom_qm_cg(i),i=1,cgrest_i%n_cg)
            read(iunit,*)cgrest_r%q_rest,cgrest_r%lambda
          ELSE IF(INDEX(LINE,'FLEXIBLE').NE.0) THEN
            solqmmm%tflexsol=.true.
            IF(INDEX(LINE,'WATER').NE.0)THEN
              solqmmi%solv_bond(1)=35
              solqmmi%solv_bond(2)=35
              solqmmi%solv_bond(3)=41
              IF(INDEX(LINE,'BONDTYPE').NE.0)THEN
                READ(iunit,*) (solqmmi%solv_bond(i),i=1,3)
              ENDIF
              IF(INDEX(LINE,'ALL').NE.0)THEN
                solqmmm%all_water=.true.
              ELSE
                solqmmm%all_water=.false.
                read(iunit,*)solqmmi%nwq
                read(iunit,*)(solqmmi%ig_wq(iw),iw=1,solqmmi%nwq)
              ENDIF         
            ELSE
              call stopgm('MM_READ_QMMM_INPUT',
     &             'PLEASE PARAMETRIZE YOUR SOLVENT',
     & __LINE__,"__file__")
            ENDIF
          ELSE IF(INDEX(LINE,'EXCLUSION').NE.0) THEN
            IF(INDEX(LINE,'LIST').NE.0) THEN
              excl_comm_l%gromos_excl=.false.
              read(iunit,*)excl_comm%nce
              do i=1,excl_comm%nce
                read(iunit,*)excl_comm%atom_qm_excl(i),
     &                excl_comm%atom_mm_excl(i)
              enddo 
              IF(excl_comm%nce.gt.NCe_max) THEN
                call stopgm('MM_READ_QMMM_INPUT','too many exclusions',
     & __LINE__,"__file__")
              ENDIF
              IF(INDEX(LINE,'NORESP').NE.0) THEN
                 excl_comm_l%excl_mech=.true.
              ENDIF
            ELSE IF(INDEX(LINE,'GROMOS').NE.0) THEN
              excl_comm_l%gromos_excl=.true.
            ELSE
C           leave room for future methods of specifying exclusions.
              call stopgm('mm_read_qmmm_input',
     &           'use EXCLUSION LIST or EXCLUSION GROMOS',
     & __LINE__,"__file__")
            ENDIF
          ELSE IF(INDEX(LINE,'SPLIT').NE.0) THEN
            sppph%mm_split=.true.
          ELSE
C           print warning about unrecognized, non-empty lines 
            CALL XSTRING(LINE,IA,IE)
            IF(IE.GE.IA) THEN
              write(6,*) 'WARNING!! UNRECOGNIZED LINE IN &QMMM SECTION'
              write(6,*) LINE
            ENDIF
          ENDIF
          goto 22
24        call stopgm('MM_READ_QMMM_INPUT','ERROR READING &QMMM BLOCK',
     & __LINE__,"__file__")
23        continue

C     check what we got.
          if(iqmmm%coupl_model.eq.0)mne%maxnat_el=1
C     required keywords
          if(.not.ft) call stopgm('MM_READ_QMMM_INPUT',
     &         'NO TOPOLOGY FILE NAME GIVEN',
     & __LINE__,"__file__")
          if(.not.fc) call stopgm('MM_READ_QMMM_INPUT',
     &         'NO COORDINATES FILE NAME GIVEN',
     & __LINE__,"__file__")
          if(.not.fi) call stopgm('MM_READ_QMMM_INPUT',
     &         'NO GROMOS INPUT FILE',
     & __LINE__,"__file__")
        endif
        

        if((excl_comm%nce.ne.0 .or.  excl_comm_l%gromos_excl)  
     &       .and. cgrest_i%n_cg.ne.0)
     c         call stopgm('MM_READ_QMMM_INPUT',
     c         'exclusion and charge restrain not compatible',
     & __LINE__,"__file__")
C_FM[
        IF(CNTL%FMATCH) THEN
          WRITE(6,*)
     c' fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
          WRITE(6,*)' fm  Force-Matching'
          WRITE(6,*)' fm '
          IF (fm_compute_traj) THEN
            WRITE(6,*)' fm  Will compute reference trajectory '
          ELSE IF (fm_compute_sp) THEN
            WRITE(6,*)' fm  Will recompute reference forces'
            WRITE(6,*)' fm  read ref frames from ',
     &             trim(fm_ref_traj_filename)
            WRITE(6,*)' fm  stride: ',fm_ref_traj_stride
            IF (fm_sp_restart) THEN
                WRITE(6,*)' fm  restarting SP series '
                WRITE(6,*)' fm  appending FM_REF* files '
            END IF
          ELSE IF (fm_read_forces) THEN
            WRITE(6,*)' fm  Will read reference forces from ',
     &                 trim(fm_forces_ref_filename)
          ENDIF
          IF (fm_fit_charges) THEN
            IF (fm_charges_only) THEN
              WRITE(6,*)' fm  Will fit ONLY charges'
            ELSE
              WRITE(6,*)' fm  Will fit charges'
            ENDIF
            IF (fm_fix_q) THEN
               WRITE(6,*)' fm  fixing charges for gromos atoms : ',
     &           fm_fixed_q_grm_indx(1:fm_n_fix_q)
C_DEBUG[
               WRITE(6,*) ' fm  ',fm_fixed_q_trgt(1:fm_n_fix_q)
C_DEBUG]
            ENDIF
            IF (fm_nexclude.gt.0) THEN
              WRITE(6,*)' fm  Nr. of atoms that are excluded from '
     c                  // 'charge fitting:', fm_nexclude
              WRITE(6,*)' fm  They are: ', fm_exclude(1:fm_nexclude)
            ENDIF
          ELSE
            WRITE(6,*)' fm  Will not fit charges'
          ENDIF
          IF (fm_optimize_weights) THEN
            WRITE(6,*)' fm  Will determine optimal weights for charge'
     c             // ' fitting'
          ELSE
            WRITE(6,*)' fm  Weights used for charge fitting:'
            WRITE(6,331) fm_wv
            WRITE(6,332) fm_wf
            WRITE(6,333)
c            WRITE(6,333) fm_wq
            WRITE(6,334) fm_wtot
          ENDIF
          IF (fm_fit_bonds) THEN
             WRITE(6,*)' fm  Will fit bond parameters'
          ELSE
             WRITE(6,*)' fm  Will not fit bond parameters'
          END IF
          IF (fm_fit_angles) THEN
             WRITE(6,*)' fm  Will fit angle parameters'
          ELSE
             WRITE(6,*)' fm  Will not fit angle parameters'
          END IF
          IF (fm_fit_diheds) THEN
             WRITE(6,*)' fm  Will fit dihedral parameters'
          ELSE
             WRITE(6,*)' fm  Will not fit dihedral parameters'
          END IF
          IF (fm_fit_improps) THEN
             WRITE(6,*)' fm  Will fit improper dihedral parameters'
          ELSE
             WRITE(6,*)' fm  Will not fit improper dihedral parameters'
          END IF
          IF (fm_fit_fc_only) THEN
             WRITE(6,*)' fm  Will fit force constants only'
          ElSE
             WRITE(6,*)' fm  Will optimize force constants and ',
     &                       'reference values'
          END IF
          WRITE(6,*)' fm '
          WRITE(6,*)' fm  Update NN LIST in every step'
          iqmmm%nt_nlist = 1
          WRITE(6,*)
     c' fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
          WRITE(6,*)
        ENDIF

 331     format("  fm  Potential:           WV = ", e9.3)
 332     format("  fm  Field:               WF = ", e9.3)
 333     format("  fm  Charge restraint:    WQ given below ")
C Manuel, weights for fixed charges 
c 333     format("  fm  Charge restraint:    WQ = ", e9.3)
 334     format("  fm  Total Charge:      WTOT = ", e9.3)
C_FM]       
! prevent splitting off the mm node if not standard cpmd
! this is checked later, but we do it here also, to have
! a consistent output.
        if(.not.CNTL%MD .or. CNTL%TPRCP) sppph%mm_split=.false.

        rqmmm%rcut_el=max(rqmmm%rcut_el,rqmmm%rcut_small)
        rqmmm%rcut_esp=max(rqmmm%rcut_esp,rqmmm%rcut_el)

          write(6,*)
     c' qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm '
          write(6,*)' qmmm by Alessandro Laio,'
          write(6,*)' qmmm    Joost VandeVondele, and'
          write(6,*)' qmmm    Ursula Rothlisberger'
          write(6,*)' qmmm '
          WRITE(6,*)' qmmm Revised dynamical Gromos arrays version:'
          WRITE(6,*)' qmmm Tsukuba/Philadelphia/Bochum: 25 Mar 2008'
          WRITE(6,*)' qmmm Mauro Boero, Axel Kohlmeyer, Nisanth Nair'
          WRITE(6,*)' qmmm '
          WRITE(6,*)' qmmm (Rewritten and vectorized for ES & Co.)'
          WRITE(6,*)' qmmm (Large topology support)'
          WRITE(6,*)' qmmm '
          write(6,*)
     c' qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm '
          write(6,*)' qmmm '
          call xstring(cqmmm%filetopo,ia,ie)
          write(6,*)' qmmm GROMOS TOPOLOGY FILE    : ',
     &         cqmmm%filetopo(ia:ie)
          call xstring(cqmmm%fileconf,ia,ie)
          write(6,*)' qmmm GROMOS COORDINATES FILE : ',
     &         cqmmm%fileconf(ia:ie)
          call xstring(cqmmm%filegrin,ia,ie)
          write(6,*)' qmmm GROMOS INPUT FILE       : ',
     &         cqmmm%filegrin(ia:ie)
          write(6,*)' qmmm '
          if(agr%qmmm_amber) then 
            write(6,*) ' qmmm USING AMBER MM FORCE FIELD'
          else
            write(6,*) ' qmmm USING GROMOS MM FORCE FIELD'
          endif
          write(6,*)' qmmm '
          IF (parai%CP_NPROC.GT.1 .and. sppph%mm_split) then
            write(6,*)' qmmm SPLITTING OFF ONE NODE FOR MM TASK'
          ELSE
            write(6,*)
     .       ' qmmm QM AND MM TASKS PERFORMED IN THE SAME GROUP'
          ENDIF
          write(6,*)' qmmm '
          if(iqmmm%coupl_model.eq.0) then
            write(6,*)' qmmm MECHANICAL COUPLING WITH THE MM ATOMS '
          else if(iqmmm%coupl_model.eq.1) then
            write(6,*)' qmmm ELECTROSTATIC COUPLING WITH THE NN ATOMS '
            write(6,*)' qmmm MECHANICAL COUPLING WITH THE FAR MM ATOMS'
          else if(iqmmm%coupl_model.eq.2) then
            write(6,*)' qmmm ELECTROSTATIC COUPLING WITH THE NN ATOMS '
            write(6,*)' qmmm MULTIPOLE COUPLING WITH THE FAR MM ATOMS '
          endif
          if(iqmmm%coupl_model.gt.0)then
            write(6,*)' qmmm '
            write(6,*)' qmmm CUTOFFS FOR THE ELECTROSTATIC ',
     &           'COUPLING [A.U.]:'
            write(6,53)rqmmm%rcut_small
            write(6,54)rqmmm%rcut_el
            write(6,55)rqmmm%rcut_esp
            write(6,*)' qmmm '
            write(6,51)iqmmm%nt_nlist
            write(6,
     c     '("  qmmm {MAXNN} MAXIMUM NUMBER OF NN ATOMS ",i6)')
     c                mne%maxnat_el
            write(6,
     c     '("  qmmm {ESPWEIGHT} WEIGHT FOR THE ESP FIT = ",f10.4)')
     c                r_esp%esp_weight
            if(n_sml%no_hirshfeld)  write(6,
     c           '("  qmmm RESTRAINT TO HIRSHFELD CHARGES OFF")')
            if(excl_comm_l%excl_mech)  write(6,
     c           '("  qmmm EXCLUSION LIST COUPLED VIA MM-CHARGES")')
            if(excl_comm_l%excl_mech)  write(6,
     c           '("  qmmm EXCLUSION LIST COUPLED VIA RESP CHARGES")')
            if(iqmmm%nt_sample.gt.0) then
              write(6,50)iqmmm%nt_sample
              if(lqmmm%qmmm_dcd) then
                write(6,'("  qmmm WRITE INTERACTING IN {DCD} FORMAT")')
              endif
            else
              write(6,*)' qmmm {SAMPLE INTERACTING} OFF'
            endif
          endif
          write(6,*)' qmmm '
          IF(solqmmm%tflexsol) THEN
            write(6,*)' qmmm CONVERTING RIGID SOLVENT MOLECULES:'
            write(6,*)' qmmm CREATING {FLEXIBLE WATER} WITH BOND TYPES'
            write(6,'("  qmmm",3I5)') (solqmmi%solv_bond(i),i=1,3)
            IF(solqmmm%all_water) THEN
              write(6,*)' qmmm {ALL} WATER MOLECULES ARE MADE FLEXIBLE.'
            ELSE
              write(6,*)' qmmm NUMBER OF CONVERTED WATER MOLECULES:',
     &              solqmmi%nwq
              write(6,*)' qmmm OXYGEN ATOM INDEX NUMBERS:'
              write(6,'("  qmmm",8I7)') (solqmmi%ig_wq(iw),
     &             iw=1,solqmmi%nwq)
            ENDIF
          ENDIF
          write(6,*)
     c' qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm qmmm '
          write(6,*)

      ENDIF
50    format("  qmmm {SAMPLE INTERACTING} DATA SAVED EVERY ",
     c          i6," STEPS")
51    format("  qmmm {UPDATE LIST} NN LIST UPDATED EVERY ",
     c          i6," STEPS")
53    format("  qmmm {RCUT_NN}  RCUT_NN  = ",f7.2)
54    format("  qmmm {RCUT_MIX} RCUT_MIX = ",f7.2)
55    format("  qmmm {RCUT_ESP} RCUT_ESP = ",f7.2)
!
      if(clc%classical) then
        if(iqmmm%coupl_model.gt.0)write(6,*)
     c    'ONLY MECHANICAL COUPLING IF FULLY CLASSICAL MD !!!'
        iqmmm%coupl_model=0
        CNTL%QUENCHB=.false.
        cnti%INWFUN=1
        cntr%ECUT=2.D0
      endif
C     ==---------------- qmmm broadcast of data ----------------------==
      CALL Mp_BCAST(lqmmm%qmmm,parai%io_source,parai%cp_grp)
      if(lqmmm%qmmm)then

C     dynamic arraysizes
        CALL Mp_BCAST_byte(MAXATT,NRTOPO*8/2,parai%io_source,
     &        parai%cp_grp)
        CALL Mp_BCAST_byte(MAXP,NRVECT*8/2,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST_byte(coordsz_t%MAXNAT,NRCOORD*8/2,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST_byte(MAXINB,NRFORCE*8/2,parai%io_source,
     &       parai%cp_grp)
C
        cqmmm%fileconf=ADJUSTL(cqmmm%fileconf)
        cqmmm%filetopo=ADJUSTL(cqmmm%filetopo)
        cqmmm%filegrin=ADJUSTL(cqmmm%filegrin)
        CALL mp_bcast(cqmmm%fileconf,parai%io_source,parai%cp_grp)
        CALL mp_bcast(cqmmm%filetopo,parai%io_source,parai%cp_grp)
        CALL mp_bcast(cqmmm%filegrin,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(TEXTFLD,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(lqmmm%qmmm_time,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(lqmmm%qmmm_verbose,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(lqmmm%qmmm_reflex,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(lqmmm%qmmm_dcd,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_l%write_potential,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(n_sml%no_hirshfeld,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_l%write_density,parai%io_source,
     &       parai%cp_grp)
        CALL mp_bcast(wp_c%cubename_dens,parai%io_source,
     &       parai%cp_grp)
        CALL mp_bcast(wp_c%cubename_pot,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_i%nfi_wd,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(wp_i%nfi_wp,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(wp_i%nfi_lt,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(iqmmm%nt_sample,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_i%n_stride,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(r_esp%esp_weight,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rqmmm%rcut_small,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rqmmm%rcut_el,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rqmmm%rcut_esp,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rqmmm%box_toll,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rqmmm%box_wall,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(iqmmm%nt_nlist,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(n_sm%n_smear,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(solqmmm%tflexsol,parai%io_source,
     &       parai%cp_grp) 
        CALL Mp_BCAST(solqmmm%all_water,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(solqmmi%nwq,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(solqmmi%ig_wq,solqmmi%nwq,parai%io_source,
     &       parai%cp_grp) 
        CALL Mp_BCAST(solqmmi%solv_bond,20,parai%io_source,
     &       parai%cp_grp) 
        CALL Mp_BCAST(iqmmm%coupl_model,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(sppph%mm_split,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(agr%qmmm_amber,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(capp%cap_h,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(ragr%scale_14_amber,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(addh%n_carbon,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(addh%ig_carbon,addh%n_carbon,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(mne%maxnat_el,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(excl_comm_l%gromos_excl,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(excl_comm_l%excl_mech,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(excl_comm%nce,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(excl_comm%atom_qm_excl,excl_comm%nce,
     &       parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(excl_comm%atom_mm_excl,excl_comm%nce,
     &       parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(cgrest_i%n_cg,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(cgrest_i%atom_qm_cg,cgrest_i%n_cg,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(cgrest_r%q_rest,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(cgrest_r%lambda,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rtr_i%iframe_restart,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rtr_l%restart_traj,parai%io_source,
     &       parai%cp_grp)
        CALL mp_bcast(rtr_c%file_traj_restart,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(rtr_r%rev_vel,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(wp_p%n_print_wf,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_l%write_orbs,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(wp_p%i_wf_pr,wp_p%n_print_wf,
     &       parai%io_source,parai%cp_grp)
C_FM[
        CALL Mp_BCAST(fm_fmatch,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(fm_compute_sp,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(fm_compute_traj,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(fm_read_forces,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(fm_rinit_wf,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(fm_sp_restart,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(fm_compute_rms,parai%io_source,
     &       parai%cp_grp)
        CALL Mp_BCAST(fm_nexclude,parai%io_source,parai%cp_grp)
        CALL Mp_BCAST(fm_exclude,fm_max_nexclude,parai%io_source,
     &       parai%cp_grp)
C        CALL SENDCHAR(fm_pip_ref_filename,80,parai%io_source,parai%cp_grp)
C        CALL SENDCHAR(fm_forces_ref_filename,80,parai%io_source,parai%cp_grp)
C        CALL SENDCHAR(fm_chj_ref_filename,80,parai%io_source,parai%cp_grp)
C_FM]
      endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SETARRAYSZ(IUNIT)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
C     read the required array sizes for a topology from 
C     input and set the rest.
C     ==--------------------------------------------------------------==
      use system
      use readsr_utils, only : readsr, xstring, readsi
      use coordsz
      implicit none

      !include 'Gromos/coordsz.h'
      include 'Gromos/forcesz.h'
      include 'Gromos/box.h'
      include 'Gromos/vector.h'
      include 'Gromos/toposz.h'

      INTEGER IUNIT, VALUE, IA, IE
      LOGICAL IERR
      CHARACTER*80 LINE
      CHARACTER*6 VARNAME

C     some default values.
C     FIXME: AK 2005/05/25
C     these defaults need to be tuned, so that only
C     those, that have a large impact on memory use
C     are set to minimal values (i.e. need to be set in the input).
C     and the rest are good for most cases.

C     from vector.h
      MAXP   = 20001
      MFORCE = 20001
      MAXANG = 4001
      JRCNJA = 12
      IVECT  = 64
      MXGGRP = 24
      MAXGP  = 90001

c     from coordsz.h
      coordsz_t%MAXNAT = 500000 ! cmb - enlarged
      coordsz_t%MAXNDR = 300
      coordsz_t%MAXNDL = 1
      coordsz_t%MAXCON = 1
      coordsz_t%MAXJ   = 9

c     from forcesz.h
      MAXINB = 320000
      MAXJNB = 18000000

C     from toposz.h
      MAXATT = 40
      MAXAA2 = 30000
      MAXNRP = 400000
      MAXNBT = 120
      MAXBNH = 400000
      MAXBON = 40000
      MAXTTY = 250
      MXQHEH = 20000
      MAXTHE = 20000
      MAXQTY = 10
      MAXHIH = 10
      MAXPTY = 90
      MXPHIH = 40000
      MAXPHI = 40000
      MAXCAG = 40000
      MAXAEX = 160000
      MXEX14 = 45000
      MAXNRS = 3
      MXCONS = 9
      
C     loop over file 
 100  CONTINUE

      READ(UNIT=IUNIT, END=115, ERR=115, FMT='(A80)') LINE
C     end of block
      IF(INDEX(LINE,'END ARRAYSI').NE.0) GOTO 200

C     split into string and integer
      CALL XSTRING(LINE,IA,IE)

C     empty line. read next
      IF(IA.GT.IE) GOTO 100

      VARNAME=LINE(IA:IE)
      IA=IE+1
      CALL READSI(LINE,IA,IE,VALUE,IERR)
      IF(IERR) GOTO 110

      IF(VALUE.LE.0) THEN
        WRITE(6,'(A,I8,3A)') 'WARNING!! ILLEGAL VALUE ', VALUE,
     &       'FOR ARRAYSIZE ', VARNAME, ' IGNORED.'
      ENDIF
      
      IF(INDEX(VARNAME,'MAXNAT').NE.0) THEN
        coordsz_t%MAXNAT=VALUE
      ELSE IF(INDEX(VARNAME,'MAXNDR').NE.0) THEN
        coordsz_t%MAXNDR=VALUE
      ELSE IF(INDEX(VARNAME,'MAXCON').NE.0) THEN
        coordsz_t%MAXCON=VALUE
      ELSE IF(INDEX(VARNAME,'MFORCE').NE.0) THEN
        MFORCE=VALUE
      ELSE IF(INDEX(VARNAME,'MAXANG').NE.0) THEN
        MAXANG=VALUE
      ELSE IF(INDEX(VARNAME,'JRCNJA').NE.0) THEN
        JRCNJA=VALUE
      ELSE IF(INDEX(VARNAME,'IVECT').NE.0) THEN
        IVECT=VALUE
      ELSE IF(INDEX(VARNAME,'MXGGRP').NE.0) THEN
        MXGGRP=VALUE
      ELSE IF(INDEX(VARNAME,'MAXGP').NE.0) THEN
        MAXGP=VALUE
      ELSE IF((INDEX(VARNAME,'MAXATT').NE.0)
     &       .OR. (INDEX(VARNAME,'NRATT').NE.0)) THEN
        MAXATT=VALUE
      ELSE IF(INDEX(VARNAME,'MAXAA2').NE.0) THEN
        MAXAA2=VALUE
      ELSE IF(INDEX(VARNAME,'MAXNRP').NE.0) THEN
        MAXNRP=VALUE
      ELSE IF(INDEX(VARNAME,'MAXNBT').NE.0) THEN
        MAXNBT=VALUE
      ELSE IF(INDEX(VARNAME,'MAXBNH').NE.0) THEN
        MAXBNH=VALUE
      ELSE IF(INDEX(VARNAME,'MAXBON').NE.0) THEN
        MAXBON=VALUE
      ELSE IF(INDEX(VARNAME,'MAXTTY').NE.0) THEN
        MAXTTY=VALUE
      ELSE IF(INDEX(VARNAME,'MXQHEH').NE.0) THEN
        MXQHEH=VALUE
      ELSE IF(INDEX(VARNAME,'MAXTH').NE.0) THEN
        MAXTHE=VALUE
      ELSE IF(INDEX(VARNAME,'MAXQTY').NE.0) THEN
        MAXQTY=VALUE
      ELSE IF(INDEX(VARNAME,'MAXHIH').NE.0) THEN
        MAXHIH=VALUE
      ELSE IF(INDEX(VARNAME,'MAXQHI').NE.0) THEN
        MAXQHI=VALUE
      ELSE IF(INDEX(VARNAME,'MAXPTY').NE.0) THEN
        MAXPTY=VALUE
      ELSE IF(INDEX(VARNAME,'MXPHIH').NE.0) THEN
        MXPHIH=VALUE
      ELSE IF(INDEX(VARNAME,'MAXPHI').NE.0) THEN
        MAXPHI=VALUE
      ELSE IF(INDEX(VARNAME,'MAXCAG').NE.0) THEN
        MAXCAG=VALUE
      ELSE IF(INDEX(VARNAME,'MAXAEX').NE.0) THEN
        MAXAEX=VALUE
      ELSE IF(INDEX(VARNAME,'MXEX14').NE.0) THEN
        MXEX14=VALUE
      ELSE IF(INDEX(VARNAME,'MAXNRS').NE.0) THEN
        MAXNRS=VALUE
      ELSE IF(INDEX(VARNAME,'MXCONS').NE.0) THEN
        MXCONS=VALUE
      ELSE IF(INDEX(VARNAME,'MAXINB').NE.0) THEN
        MAXINB=VALUE
      ELSE IF(INDEX(VARNAME,'MAXJNB').NE.0) THEN
        MAXJNB=VALUE
C     MAXJ and MAXP have to come last, so they do not
C     match on MAXPHI or MAXJNB. AK.
      ELSE IF(INDEX(VARNAME,'MAXJ').NE.0) THEN
        coordsz_t%MAXJ=VALUE
      ELSE IF(INDEX(VARNAME,'MAXP').NE.0) THEN
        MAXP=VALUE
      ELSE
        WRITE(6,'(3A)') ' WARNING!! UNRECOGNIZED ARRAYSIZE ',
     &       VARNAME, ' IGNORED.'
      ENDIF
      GOTO 100
      
C     parse error. warn but continue
 110  CONTINUE
      WRITE(6,'(2A)') ' WARNING!!! COULD NOT PARSE ARRAYSIZE: ', LINE
      GOTO 100
      
C     error stop.
 115  CONTINUE
      call stopgm('SETARRAYS','ERROR READING ARRAYSIZES SECTION',
     & __LINE__,"__file__")
      
C     regular exit
 200  CONTINUE
      
C     set derived values
C     from vector.h
      IEPAIR = MXGGRP * (MXGGRP + 1) / 2
      JRCNPC = 5*MAXP

C     from coordsz.h
      coordsz_t%MAXXCO = NDIM*coordsz_t%MAXNAT
      coordsz_t%MAX4DA = coordsz_t%MAXNAT
      coordsz_t%MAXXC  = coordsz_t%MAXNAT
      coordsz_t%MAXTDR = coordsz_t%MAXNDR

C     from toposz.h
      MXATT2 = MAXATT*(MAXATT+1)/2
      MAXNP2 = 2*MAXNRP
C     
      VALUE  = MAX(MAXQHI,MAXHIH)
      MAXQHI = VALUE 
      MAXHIH = VALUE

      WRITE(6,'(/," ",6("-ARRAYSIZES-"))')
      WRITE(6,333) 'IEPAIR', IEPAIR
      WRITE(6,333) 'IVECT',  IVECT
      WRITE(6,333) 'JRCNJA', JRCNJA
      WRITE(6,333) 'JRCNPC', JRCNPC
      WRITE(6,333) 'MAX4DA', coordsz_t%MAX4DA
      WRITE(6,333) 'MAXAA2', MAXAA2
      WRITE(6,333) 'MAXAEX', MAXAEX
      WRITE(6,333) 'MAXANG', MAXANG
      WRITE(6,333) 'MAXATT', MAXATT
      WRITE(6,333) 'MAXBNH', MAXBNH
      WRITE(6,333) 'MAXBON', MAXBON
      WRITE(6,333) 'MAXCAG', MAXCAG
      WRITE(6,333) 'MAXCON', coordsz_t%MAXCON
      WRITE(6,333) 'MAXGP',  MAXGP
      WRITE(6,333) 'MAXHIH', MAXHIH
      WRITE(6,333) 'MAXJ',   coordsz_t%MAXJ
      WRITE(6,333) 'MAXNAT', coordsz_t%MAXNAT
      WRITE(6,333) 'MAXNBT', MAXNBT
      WRITE(6,333) 'MAXNDL', coordsz_t%MAXNDL
      WRITE(6,333) 'MAXNDR', coordsz_t%MAXNDR
      WRITE(6,333) 'MAXNP2', MAXNP2
      WRITE(6,333) 'MAXNRP', MAXNRP
      WRITE(6,333) 'MAXNRS', MAXNRS
      WRITE(6,333) 'MAXP',   MAXP
      WRITE(6,333) 'MAXPHI', MAXPHI
      WRITE(6,333) 'MAXPTY', MAXPTY
      WRITE(6,333) 'MAXQTY', MAXQTY
      WRITE(6,333) 'MAXTDR', coordsz_t%MAXTDR
      WRITE(6,333) 'MAXTHE', MAXTHE
      WRITE(6,333) 'MAXTTY', MAXTTY
      WRITE(6,333) 'MAXXC',  coordsz_t%MAXXC
      WRITE(6,333) 'MAXXCO', coordsz_t%MAXXCO
      WRITE(6,333) 'MFORCE', MFORCE
      WRITE(6,333) 'MXATT2', MXATT2
      WRITE(6,333) 'MXCONS', MXCONS
      WRITE(6,333) 'MXEX14', MXEX14
      WRITE(6,333) 'MXGGRP', MXGGRP
      WRITE(6,333) 'MXPHIH', MXPHIH
      WRITE(6,333) 'MXQHEH', MXQHEH
      WRITE(6,333) 'MAXINB', MAXINB
      WRITE(6,333) 'MAXJNB', MAXJNB

 333  FORMAT (4X,A8,2X,I8)
      WRITE(6,'(" ",6("-ARRAYSIZES-"),/)')

C     ==--------------------------------------------------------------==
      RETURN
      END
C_FM[
C     ==================================================================
      SUBROUTINE READ_FMATCH_INPUT(IUNIT)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
C     read forcematch block
C     ==--------------------------------------------------------------==
C_DEBUG[
      use mm_input
C_DEBUG]
      use forcematch
      IMPLICIT NONE
C
      INTEGER IUNIT, VALUE, IA, IE
      INTEGER IEQ, I, IAT, NREQ,IQFIX,iqweight
      LOGICAL IERR
      CHARACTER*80 LINE
C
C     loop over input block
 300  CONTINUE
C
      READ(UNIT=IUNIT, END=315, ERR=315, FMT='(A80)') LINE
      IF(INDEX(LINE,'END FORCEMATCH').NE.0) GOTO 399
C
      IF(INDEX(LINE,'READ REF FORCES').NE.0) THEN
        fm_compute_traj = .FALSE.
        fm_compute_sp = .FALSE.
        fm_read_forces = .TRUE.
        IF(INDEX(LINE,'FILE').NE.0) THEN
          READ(IUNIT,'(A)') fm_forces_ref_filename
        ENDIF
        IF(INDEX(LINE,'COVALENT').NE.0)THEN
          fm_read_covforces = .TRUE.
        ENDIF
      ELSE IF(INDEX(LINE,'READ REF TRAJ').NE.0) THEN
        fm_compute_sp = .TRUE.
        fm_compute_traj = .FALSE.
        fm_read_forces = .FALSE.
C        write(6,*) 'default traj ref name ',fm_ref_traj_filename
        IF(INDEX(LINE,'FILE').NE.0) THEN
          READ(IUNIT,'(A)') fm_ref_traj_filename
C          write(6,*) 'read traj ref name ',fm_ref_traj_filename
        ENDIF
      ELSE IF(INDEX(LINE,'COMPUTE REF TRAJ').NE.0) THEN
        fm_compute_sp = .FALSE.
        fm_compute_traj = .TRUE.
        fm_read_forces = .FALSE.
      ELSE IF(INDEX(LINE,'RESTART SP').NE. 0) THEN
        fm_sp_restart=.TRUE.
        fm_compute_sp = .TRUE.
        fm_compute_traj = .FALSE.
        fm_read_forces = .FALSE.
      ELSE IF(INDEX(LINE,'READ REF STRIDE').NE.0) THEN
         READ(IUNIT,*) fm_ref_traj_stride
      ELSE IF(INDEX(LINE,'TOPOL').NE.0 
     $        .AND. INDEX(LINE,'OUT').NE.0) THEN
        READ(IUNIT,'(A)') fm_outtopo_filename
      ELSE IF(INDEX(LINE,'INITWF').NE.0) THEN
        IF(INDEX(LINE,'OFF').NE.0) THEN
          fm_rinit_wf = .FALSE.
        ELSE
          fm_rinit_wf = .TRUE.
        ENDIF
      ELSE IF (INDEX(LINE,'CHARGES').NE.0) THEN
        IF(INDEX(LINE,'ONLY').NE.0) THEN
          fm_charges_only = .true.
        ENDIF
        IF (INDEX(LINE,'NO').NE.0) THEN
          fm_fit_charges = .false.
        ELSE 
          fm_fit_charges = .true.
        ENDIF
        IF (INDEX(LINE,'FIX').NE.0) THEN
           fm_fit_charges = .true.
           fm_fix_q =.true.
           DO IQFIX=1,fm_max_q_fixers
             fm_fixed_q_indx(IQFIX)=0
             fm_fixed_q_trgt(IQFIX)=0.0d0
           END DO
           READ(IUNIT, *) fm_n_fix_q
           IF(fm_n_fix_q .GT. fm_max_q_fixers) THEN 
            call stopgm('read_fmatch_input',
     &             'fm_n_fix_q > fm_max_q_fixers!',
     & __LINE__,"__file__")
           END IF    
           DO IQFIX=1,fm_n_fix_q
            READ(IUNIT, *) fm_fixed_q_grm_indx(IQFIX),
     &          fm_fixed_q_trgt(IQFIX)
           END DO  
         END IF  
      ELSE IF(INDEX(LINE,'OPT').NE.0 .AND. INDEX(LINE,'QFIT').NE.0 .AND.
     $        INDEX(LINE,'WEIGHT').NE.0) THEN
        fm_optimize_weights = .TRUE.
      ELSE IF(INDEX(LINE,'WV').NE.0) THEN
        fm_optimize_weights = .FALSE.
        READ(IUNIT, *) fm_wv
      ELSE IF(INDEX(LINE,'WF').NE.0) THEN
        fm_optimize_weights = .FALSE.
        READ(IUNIT, *) fm_wf
      ELSE IF(INDEX(LINE,'WQ INDIVIDUAL').NE.0) THEN
        fm_optimize_weights = .FALSE.
        READ(IUNIT, *) fm_n_qw
        do iqweight=1,fm_n_qw
         READ(IUNIT, *) fm_wq_grm_indx(iqweight),
     &                 fm_wq(fm_wq_grm_indx(iqweight))
        end do  
      ELSE IF(INDEX(LINE,'WQ GENERAL').NE.0) THEN
        fm_optimize_weights = .FALSE.
        READ(IUNIT, *) fm_wq_general        
      ELSE IF(INDEX(LINE,'WTOT').NE.0) THEN
        READ(IUNIT, *) fm_wtot
      ELSE IF(INDEX(LINE,'EQUIV').NE.0) THEN
        fm_equivalence = .TRUE.
        DO IEQ = 1, fm_maxeq
          fm_equiv(1,IEQ) = IEQ
          fm_equiv(2,IEQ) = IEQ
        ENDDO
        READ(IUNIT, *) fm_nreq
        IF(fm_NREQ.GT.fm_maxeq) call stopgm('read_fmatch_input', 
     $                                   'fm_nreq > fm_maxeq!',
     & __LINE__,"__file__")
        DO IEQ = 1, fm_nreq
          READ(IUNIT, *) fm_equiv(1,IEQ), fm_equiv(2,IEQ)
        ENDDO
      ELSE IF (INDEX(LINE,'OPT').NE. 0 .AND. INDEX(LINE,'FC').NE.0 
     $        .AND. INDEX(LINE,'ONLY') .NE. 0) THEN
        fm_fit_fc_only=.true.
      ELSE IF (INDEX(LINE,'NO BONDS') .NE. 0) THEN
         fm_fit_bonds=.false.
      ELSE IF (INDEX(LINE,'NO ANGLES') .NE. 0) THEN
         fm_fit_angles=.false.
      ELSE IF (INDEX(LINE,'NO DIHEDRALS') .NE. 0) THEN
         fm_fit_diheds=.false.
      ELSE IF (INDEX(LINE,'NO IMPROPERS') .NE. 0) THEN
         fm_fit_improps=.false.
      ELSE IF(INDEX(LINE,'FMCAP').NE.0) THEN
        fm_capping = .true.
        READ(IUNIT, *) fm_ncap
        DO i = 1, fm_ncap
          READ(IUNIT, *) fm_cap(1, i), fm_cap(2, i)
        ENDDO
      ELSE IF(INDEX(LINE,'MAXITER').NE.0) THEN
        IF(INDEX(LINE,'LM').NE.0) THEN
          READ(IUNIT,*) fm_max_kfit_iterations
        ELSE
          READ(IUNIT, *)
        ENDIF
      ELSE IF(INDEX(LINE,'COMPUTE RMS').NE.0) THEN
        IF (INDEX(LINE,'NO').NE.0) THEN
          fm_compute_rms = .false.
        ELSE
          fm_compute_rms = .true.
        ENDIF
      ELSE IF(INDEX(LINE,'FM_EXCLUDE').NE.0) THEN
        READ(IUNIT,*) fm_nexclude
        IF(fm_nexclude.GT.fm_max_nexclude) THEN
          call stopgm('READ_FMATCH_INPUT',
     $                'fm_nexclude > fm_max_nexclude! Change parameter '
     $             // 'fm_max_nexclude in forcematch.inc and recompile',
     & __LINE__,"__file__")
        ELSE
          do i = 1, fm_nexclude
            READ(IUNIT,*) fm_exclude(i)
          enddo
        ENDIF
C This is for testing/debugging
      ELSE IF(INDEX(LINE,'CLASSICAL').NE.0) THEN
        clc%classical = .true.
      ENDIF
      GOTO 300
C     error stop
 315  CONTINUE
      call stopgm('READ_FMATCH_INPUT',
     $            'ERROR READING FORCEMATCH SECTION',
     & __LINE__,"__file__")
C     regular exit
 399  CONTINUE
      RETURN
      END
C     ==================================================================
C_FM]

