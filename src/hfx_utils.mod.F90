#define PRINT_GROUP_INFOS .FALSE.

MODULE hfx_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: scgx
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             inzs,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             nzff,&
                                             nzfs
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE func,                            ONLY: func1,&
                                             func3
  USE geq0mod,                         ONLY: geq0
  USE hfxmod,                          ONLY: hfxc3,&
                                             hfxc4,&
                                             ipoolhfx,&
                                             wcentx
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_cputime,&
                                             m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_elem_to_proc,&
                                             part_1d_get_elem,&
                                             part_1d_nbr_elems,&
                                             part_1d_symm_holds_pair,&
                                             proc_to_grid2d
  USE pbc_utils,                       ONLY: pbc
  USE pslo,                            ONLY: pslo_com
  USE ropt,                            ONLY: infi,&
                                             infw,&
                                             iteropt
  USE rswfmod,                         ONLY: maxstatesx,&
                                             rswfx
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: &
       add_wfn, copy_im_to_im, copy_im_to_re, copy_re_to_im, copy_re_to_re, &
       set_psi_1_state_g, set_psi_2_states_g, set_psi_im_state_r, &
       set_psi_re_state_r, zero_wfn
  USE system,                          ONLY: cntl,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: hfx_old
  !public :: check_occupation
  PUBLIC :: get_wannier_separation
  !public :: hfxab
  !public :: hfxab2
  !public :: hfxaa
  PUBLIC :: hfxpsi_old
  !public :: hfxpb
  PUBLIC :: hfxrpa_old
  !public :: hfxrpav

CONTAINS

  SUBROUTINE hfx_old(c0,c2,f,psia,nstate,ehfx,vhfx)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE                                        ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ehfx, vhfx

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx_old'
    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8), &
      zzero = (0.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE             :: psic(:), vpotg(:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psib
    INTEGER :: geo_iter, i, ia, ib1, ib2, ibin, ibuf, id_states(2,2), ierr, &
      ig, int_i, int_ij, int_j, ispin, isub, isub2, isub3, isub4, isub5, &
      isub6, iwf, j, jb1, jb2, my_pcol, my_prow, n_int, nbr_states(2), &
      npcols, nprows, nspins, nst, nstates(2), st_offst
    INTEGER(int_8) :: nbr_int_skiped_1, nbr_int_skiped_2, nbr_integrals, &
      nbr_rwfn_precomputed, nbr_rwfn_recomputed
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_vals, new_vals, &
                                                psi_in_core
    INTEGER, SAVE                            :: int_count = 0, &
                                                prev_geo_id = -HUGE(0)
    LOGICAL                                  :: init_ints, no_more_states
    REAL(real_8) :: bin_max, bin_range, dab, dt, dt_tot, EHFX_loc, &
      EHFX_loc_1, EHFX_loc_2, GHFX_loc, GHFX_loc_1, GHFX_loc_2, max_dab, &
      old_DWFC, pfl, pfx, pfx1, pfx2, t1, t1_tot, t2, t2_tot
    REAL(real_8), ALLOCATABLE                :: vpotr(:)
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: acc_vals, max_vals
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: int_vals

    ehfx=0._real_8
    vhfx=0._real_8
    IF (func1%mhfx.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    ! decide if reinit intergrals
    ! need to use the *right* geo iter counter (in the cpmd spirit)!
    init_ints=.FALSE.
    geo_iter=0
    IF(.NOT.cntl%wfopt) geo_iter=infi
    IF (prev_geo_id/=geo_iter) THEN
       prev_geo_id=geo_iter
       init_ints=.TRUE.
    ENDIF


    IF (paral%io_parent.AND.hfxc3%twscr.AND..FALSE.) THEN
       WRITE(6,*) '>>>>>>>>>>     WFOPT=',cntl%wfopt
       WRITE(6,*) '>>>>>>>>>>       NFI=',iteropt%nfi
       WRITE(6,*) '>>>>>>>>>>      INFI=',infi!<
       WRITE(6,*) '>>>>>>>>>>     IINFI=',iteropt%iinfi
       WRITE(6,*) '>>>>>>>>>>      INFW=',infw
       WRITE(6,*) '>>>>>>>>>> init_ints=',init_ints
    ENDIF

    IF (paral%io_parent.AND.init_ints.AND.hfxc3%twscr.AND..FALSE.) THEN
       WRITE(6,*) 'TWSCR=',hfxc3%twscr
       WRITE(6,*) 'TWFT=',hfxc3%twft,&
            'DWF_INTEGRAL_THRESH=',hfxc4%dwf_integral_thresh
       WRITE(6,*) 'TWFC=',hfxc3%twfc,' DWFC=',hfxc4%dwfc
       WRITE(6,*) 'TDIAW=',hfxc3%tdiaw
    ENDIF

    t1_tot=m_cputime()

    IF (cntl%cdft.AND.parai%cp_nogrp>1) CALL stopgm(procedureN,&
         'CDFT AND CP_NOGRP>1 NOT YET IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    ! Prepare the indeces for parallelization over the groups
    ! 
    CALL proc_to_grid2d(parai%cp_nogrp,parai%cp_inter_me,nprows,npcols,&
         my_prow,my_pcol)

    ! 
    ! Need to allocate a buffer to accumulate the x part of C2
    ! 
    ALLOCATE(C2_hfx(ncpw%ngw,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)

    CALL zero_wfn(jgw,nstate,C2_hfx,ncpw%ngw)
    ! 
    ! For the moment we synchronize the C0 and C2 if the cp groups are used
    ! That should be done else where
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_a',isub5)
       CALL mp_bcast(c0,ncpw%ngw*nstate,0,parai%cp_inter_grp)
       CALL mp_bcast(c2,ncpw%ngw*nstate,0,parai%cp_inter_grp)
       CALL tihalt(procedureN//'_a',isub5)
    ENDIF

    ! 
    ! randomize the states to improve load balance
    ! 
    ! TODO

    ! 
    t1=m_cputime()

    ! CALL check_occupation(F,NSTATE,cntl%tlsd)

    ! 
    ! set the spins stuff
    ! 
    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
       pfl = 0.5_real_8
    ELSE
       nspins = 1
       nstates(1) = nstate
       pfl = 0.25_real_8
    ENDIF
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(2*jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(2*llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi_in_core(0:nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    psi_in_core(:) = 0

    IF (hfxc3%twscr) THEN
       IF (.NOT.ALLOCATED(int_vals)) THEN
          ALLOCATE(int_vals( nstate*(nstate+1)/2, 2 ),stat=ierr)
          IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
               __LINE__,__FILE__)
          CALL zeroing(int_vals)!,SIZE(int_vals))
       ENDIF
       ALLOCATE(acc_vals( nstate*(nstate+1)/2, 2 ),&
            new_vals( nstate*(nstate+1)/2 ), stat=ierr )
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL zeroing(acc_vals)!,SIZE(acc_vals))
       CALL zeroing(new_vals)!,SIZE(new_vals))

       int_count=int_count+1
    ENDIF
    IF (func1%mhfx.EQ.1) THEN
       ! HARTREE-FOCK
       ALLOCATE(psib(maxfftn,2),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL zeroing(psib)!,2*maxfftn)

       ! ==--------------------------------------------------------------==
       CALL tiset(procedureN//'_outer',isub2)
       ! 
       ! loop over the spin
       ! 
       id_states(:,:) = 0
       nbr_states(:) = 0
       nbr_int_skiped_1 = 0; nbr_int_skiped_2 = 0;
       nbr_integrals = 0
       nbr_rwfn_precomputed = 0
       nbr_rwfn_recomputed = 0
       DO ispin = 1,nspins

          nst = nstates(ispin)
          st_offst = 0
          IF (ispin.EQ.2) st_offst = nstates(1)



          ! ==--------------------------------------------------------------==
          ! ==--------------------------------------------------------------==
          ! 
          ! if possible precompute some real space wavefunctions
          ! 
          IF (cntl%krwfn) THEN
             CALL tiset(procedureN//'_precomp',isub4)

             IF (.FALSE.) THEN

                ! >>>>>>>>>>>>>> OLD
                iwf = 1
                psi_in_core(:) = 0
                DO ia=st_offst+1,st_offst+nst,2
                   IF (iwf.GT.maxstatesx/2) EXIT
                   CALL zeroing(psia)!,maxfftn)
                   IF (ia.EQ.nstate) THEN
                      CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                      psi_in_core(ia) = iwf
                   ELSE
                      CALL set_psi_2_states_g(c0(:,ia),c0(:,ia+1),psia)
                      psi_in_core(ia)   = iwf
                      psi_in_core(ia+1) = -iwf
                   ENDIF
                   nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                   CALL invfftn(psia,.TRUE.,parai%allgrp)
                   CALL dcopy(2*llr1,psia(1),1,rswfx(1,iwf),1)
                   iwf = iwf + 1
                ENDDO         ! IA

             ELSE

                ! >>>>>>>>>>>>>> NEW
                iwf = 1
                psi_in_core(:) = 0
                i_loop: DO i = 1,part_1d_nbr_elems(nst,my_prow,nprows)
                   ia = part_1d_get_elem(i,my_prow,nprows)+st_offst
                   IF (f(ia).LT.1.e-6_real_8) CYCLE
                   IF (psi_in_core(ia).EQ.0) THEN
                      ! here compute ia i f needed
                      ! 
                      CALL zeroing(psia)!,maxfftn)
                      ! psi_in_core(ia) ->  iwf ie psi(ia) stored in real at position iwf
                      ! psi_in_core(ia) -> -iwf ie psi(ia) stored in imag at position iwf
                      ! psi_in_core(ia) ->    0 ie psi(ia) not stored
                      CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                      nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                      CALL invfftn(psia,.TRUE.,parai%allgrp)
                      IF (ABS(iwf).LE.maxstatesx/2) THEN
                         psi_in_core(ia) = iwf
                         IF (iwf.LT.0) THEN
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                         ELSE
                            CALL copy_re_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ENDIF
                      ELSE
                         ! We can exit the loops
                         EXIT i_loop
                      ENDIF
                   ENDIF

                   j = 1
                   DO
                      no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,npcols)
                      IF (no_more_states) EXIT

                      ! find out the correct ib's
                      ib1 = 0; ib2 = 0;
                      DO
                         ! set id1
                         IF (ib1.EQ.0) THEN
                            ib1 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                            j = j + 1
                         ENDIF
                         IF (ia.EQ.ib1) ib1 = 0
                         IF (.NOT.part_1d_symm_holds_pair(ia,ib1)) ib1 = 0
                         IF (ib1.NE.0) THEN
                            IF (f(ib1).LT.1.e-6_real_8) ib1 = 0
                         ENDIF
                         ! First: check that the separation distance is smaller that the max distance allowed
                         IF (ib1.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                            CALL get_wannier_separation(ib1,ia,dab)
                            IF (dab.GT.hfxc4%dwfc) THEN
                               int_i=MIN(ia,ib1)
                               int_j=MAX(ia,ib1)
                               int_ij=int_j*(int_j-1)/2+int_i
                               ! Set the int_vals to something that if the 
                               ! pair comes within the radius any time later
                               ! we properly recompute the integral.
                               int_vals(int_ij,:)=1.0_real_8
                               nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                               ib1=0
                            ENDIF
                         ENDIF
                         ! Second: check that the integral is larger than the integral threshold
                         IF (.NOT.init_ints) THEN
                            IF (hfxc3%twscr.AND.ib1.NE.0) THEN
                               int_i=MIN(ia,ib1)
                               int_j=MAX(ia,ib1)
                               int_ij=int_j*(int_j-1)/2+int_i
                               IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) ib1 = 0
                            ENDIF
                         ENDIF
                         IF (psi_in_core(ib1).NE.0) ib1 = 0
                         no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                              npcols)
                         IF (no_more_states) EXIT
                         ! set ib2
                         IF (ib2.EQ.0) THEN
                            ib2 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                            j = j + 1
                         ENDIF
                         IF (ia.EQ.ib2) ib2 = 0
                         IF (.NOT.part_1d_symm_holds_pair(ia,ib2)) ib2 = 0
                         IF (ib2.NE.0) THEN
                            IF (f(ib2).LT.1.e-6_real_8) ib2 = 0
                         ENDIF
                         ! First: check that the separation distance is smaller that the max distance allowed
                         IF (ib2.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                            CALL get_wannier_separation(ib2,ia,dab)
                            IF (dab.GT.hfxc4%dwfc) THEN
                               int_i=MIN(ia,ib2)
                               int_j=MAX(ia,ib2)
                               int_ij=int_j*(int_j-1)/2+int_i
                               ! Set the int_vals to something that if the 
                               ! pair comes within the radius any time later
                               ! we properly recompute the integral.
                               int_vals(int_ij,:)=1.0_real_8
                               nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                               ib2=0
                            ENDIF
                         ENDIF
                         ! Second: check that the integral is larger than the integral threshold
                         IF (.NOT.init_ints) THEN
                            IF (hfxc3%twscr.AND.ib2.NE.0) THEN
                               int_i=MIN(ia,ib2)
                               int_j=MAX(ia,ib2)
                               int_ij=int_j*(int_j-1)/2+int_i
                               IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) ib2 = 0
                            ENDIF
                         ENDIF
                         IF (psi_in_core(ib2).NE.0) ib2 = 0
                         no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                              npcols)
                         IF (no_more_states) EXIT
                         ! are we set?
                         IF (ib1.NE.0.AND.ib2.NE.0) EXIT
                      ENDDO
                      IF (ib1.EQ.0.AND.ib2.EQ.0) CYCLE
                      IF (ib1.EQ.0) THEN
                         ib1 = ib2; ib2 = 0
                      ENDIF
                      ! compute ib1 and ib2 if needed
                      ! 
                      CALL zeroing(psia)!,maxfftn)
                      ! psi_in_core(ia) ->  iwf ie psi(ia) stored in real at position iwf
                      ! psi_in_core(ia) -> -iwf ie psi(ia) stored in imag at position iwf
                      ! psi_in_core(ia) ->    0 ie psi(ia) not stored
                      IF (ib2.EQ.0) THEN
                         CALL set_psi_1_state_g(zone,c0(:,ib1),psia)
                      ELSE
                         CALL set_psi_2_states_g(c0(:,ib1),c0(:,ib2),psia)
                      ENDIF
                      nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                      CALL invfftn(psia,.TRUE.,parai%allgrp)
                      IF (ib2.EQ.0) THEN
                         psi_in_core(ib1) = iwf
                         IF (iwf.LT.0) THEN
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                         ELSE
                            CALL copy_re_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ENDIF
                      ELSE
                         IF (iwf.LT.0) THEN
                            psi_in_core(ib1) =  iwf
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                            IF (ABS(iwf).GT.maxstatesx/2) EXIT i_loop
                            psi_in_core(ib2) =  iwf
                            CALL copy_im_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ELSE
                            psi_in_core(ib1) =  iwf
                            psi_in_core(ib2) = -iwf
                            CALL dcopy(2*llr1,psia(1),1,rswfx(1,iwf),1)
                            iwf = iwf + 1
                         ENDIF
                      ENDIF
                      IF (ABS(iwf).GT.maxstatesx/2) EXIT i_loop
                   ENDDO
                ENDDO i_loop ! i
             ENDIF
             CALL tihalt(procedureN//'_precomp',isub4)
          ENDIF
          ! 
          ! ==--------------------------------------------------------------==
          ! ==--------------------------------------------------------------==


          DO i = 1,part_1d_nbr_elems(nst,my_prow,nprows)
             ia = part_1d_get_elem(i,my_prow,nprows)+st_offst
             IF (f(ia).LT.1.e-6_real_8) CYCLE

             IF (psi_in_core(ia).NE.0) THEN
                iwf = psi_in_core(ia)
                IF (iwf.LT.0) THEN
                   iwf=-iwf
                   CALL set_psi_im_state_r(zone,rswfx(:,iwf),&
                        zzero,psia)
                ELSE
                   CALL set_psi_re_state_r(zone,rswfx(:,iwf),&
                        zzero,psia)
                ENDIF
             ELSE
                CALL zeroing(psia)!,maxfftn)
                CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                CALL invfftn(psia,.TRUE.,parai%allgrp)
                nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
             ENDIF
             ! Handle the diagonal terms
             ! Only when needed....
             IF (part_1d_elem_to_proc(ia-st_offst,nprows).EQ.my_prow.AND.&
                  part_1d_elem_to_proc(ia-st_offst,npcols).EQ.my_pcol) THEN
                IF (parai%me.EQ.0)nbr_integrals=nbr_integrals+1
                pfx=pfl*f(ia)*f(ia)
                CALL hfxaa(EHFX_loc,GHFX_loc,pfx,psia,vpotg,vpotr,psic,&
                     C2_hfx(1,ia))
                ehfx=ehfx+EHFX_loc

                IF (hfxc3%twscr) THEN
                   int_i=MIN(ia,ia)
                   int_j=MAX(ia,ia)
                   int_ij=int_j*(int_j-1)/2+int_i
                   acc_vals( int_ij, 1 ) = EHFX_loc
                   acc_vals( int_ij, 2 ) = GHFX_loc
                   new_vals( int_ij ) = 1
                ENDIF
                IF (hfxc3%twscr.AND.hfxc3%tdiaw) GOTO 101
             ENDIF
             ! Now for the other terms IA-IB
             ! ==--------------------------------------------------------------==
             CALL tiset(procedureN//'_inner',isub3)
             j = 1
             DO
                no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                     npcols)
                IF (no_more_states) EXIT

                ! 
                ! find out the correct ib's
                ib1 = 0; ib2 = 0;
                DO
                   ! set id1
                   IF (ib1.EQ.0) THEN
                      ib1 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                      j = j + 1
                   ENDIF
                   IF (ia.EQ.ib1) ib1 = 0
                   IF (.NOT.part_1d_symm_holds_pair(ia,ib1)) ib1 = 0
                   IF (ib1.NE.0) THEN
                      IF (f(ib1).LT.1.e-6_real_8) ib1 = 0
                   ENDIF
                   ! First: check that the separation distance is smaller that the max distance allowed
                   IF (ib1.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                      CALL get_wannier_separation(ib1,ia,dab)
                      IF (dab.GT.hfxc4%dwfc) THEN
                         int_i=MIN(ia,ib1)
                         int_j=MAX(ia,ib1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         ! Set the int_vals to something that if the 
                         ! pair comes within the radius any time later
                         ! we properly recompute the integral.
                         int_vals(int_ij,:)=1.0_real_8
                         nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                         ib1=0
                      ENDIF
                   ENDIF
                   ! Second: check that the integral is larger than the integral threshold
                   IF (.NOT.init_ints) THEN
                      IF (hfxc3%twscr.AND.ib1.NE.0) THEN
                         int_i=MIN(ia,ib1)
                         int_j=MAX(ia,ib1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) THEN
                            nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                            ib1 = 0
                         ENDIF
                      ENDIF
                   ENDIF
                   no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                        npcols)
                   IF (no_more_states) EXIT
                   ! set ib2
                   IF (ib2.EQ.0) THEN
                      ib2 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                      j = j + 1
                   ENDIF
                   IF (ia.EQ.ib2) ib2 = 0
                   IF (.NOT.part_1d_symm_holds_pair(ia,ib2)) ib2 = 0
                   IF (ib2.NE.0) THEN
                      IF (f(ib2).LT.1.e-6_real_8) ib2 = 0
                   ENDIF
                   ! First: check that the separation distance is smaller that the max distance allowed
                   IF (ib2.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                      CALL get_wannier_separation(ib2,ia,dab)
                      IF (dab.GT.hfxc4%dwfc) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,ib2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         ! Set the int_vals to something that if the 
                         ! pair comes within the radius any time later
                         ! we properly recompute the integral.
                         int_vals(int_ij,:)=1.0_real_8
                         nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                         ib2=0
                      ENDIF
                   ENDIF
                   ! Second: check that the integral is larger than the integral threshold
                   IF (.NOT.init_ints) THEN
                      IF (hfxc3%twscr.AND.ib2.NE.0) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,ib2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) THEN
                            nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                            ib2 = 0
                         ENDIF
                      ENDIF
                   ENDIF
                   no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                        npcols)
                   IF (no_more_states) EXIT
                   ! are we set?
                   IF (ib1.NE.0.AND.ib2.NE.0) EXIT
                ENDDO

                ! if(me.eq.0) write(6,*)cp_me,'ia,ib1,ib2',ia,ib1,ib2,no_more_states

                IF (.NOT.(ib1.EQ.0.AND.ib2.EQ.0)) THEN
                   IF ((psi_in_core(ib1).NE.0.AND.psi_in_core(ib2).NE.0)&
                        .AND.(ib1.NE.0.OR.ib2.NE.0)) THEN
                      CALL zeroing(psib(:,1))!,maxfftn)
                      IF (ib1.NE.0) THEN
                         iwf = psi_in_core(ib1)
                         IF (iwf.LT.0) THEN
                            iwf=-iwf
                            CALL set_psi_im_state_r(zone,rswfx(:,iwf),&
                                 zone,psib(:,1))
                         ELSE
                            CALL set_psi_re_state_r(zone,rswfx(:,iwf),&
                                 zone,psib(:,1))
                         ENDIF
                      ENDIF
                      IF (ib2.NE.0) THEN
                         iwf = psi_in_core(ib2)
                         IF (iwf.LT.0) THEN
                            iwf=-iwf
                            CALL set_psi_im_state_r(uimag,&
                                 rswfx(:,iwf),zone,psib(:,1))
                         ELSE
                            CALL set_psi_re_state_r(uimag,&
                                 rswfx(:,iwf),zone,psib(:,1))
                         ENDIF
                      ENDIF
                   ELSE
                      IF (ib1.EQ.0) THEN
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_1_state_g(uimag,c0(:,ib2),psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
                      ELSEIF (ib2.EQ.0) THEN
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_1_state_g(zone,c0(:,ib1),psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
                      ELSE
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_2_states_g(c0(:,ib1),c0(:,ib2),&
                              psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 2
                      ENDIF
                      ! Transform the wavefunction to real space
                      CALL invfftn(psib(:,1),.TRUE.,parai%allgrp)
                   ENDIF
                ENDIF
                ! 
                ! ==--------------------------------------------------------------==
                id_states(1,1) = ib1
                id_states(2,1) = ib2
                nbr_states(1) = 0
                IF (ib1.NE.0) nbr_states(1) = nbr_states(1) + 1
                IF (ib2.NE.0) nbr_states(1) = nbr_states(1) + 1

                IF (nbr_states(1).EQ.1) THEN
                   IF (id_states(1,1).NE.0) THEN
                      IF (id_states(1,2).EQ.0) THEN
                         id_states(1,2) = id_states(1,1)
                         CALL copy_re_to_re(llr1,psib(:,1),psib(:,2))
                      ELSE
                         id_states(2,2) = id_states(1,1)
                         CALL copy_re_to_im(llr1,psib(:,1),psib(:,2))
                      ENDIF
                      id_states(1,1) = 0
                   ELSE
                      IF (id_states(1,2).EQ.0) THEN
                         id_states(1,2) = id_states(2,1)
                         CALL copy_im_to_re(llr1,psib(:,1),psib(:,2))
                      ELSE
                         id_states(2,2) = id_states(2,1)
                         CALL copy_im_to_im(llr1,psib(:,1),psib(:,2))
                      ENDIF
                      id_states(2,1) = 0
                   ENDIF
                   nbr_states(1) = nbr_states(1) - 1
                   nbr_states(2) = nbr_states(2) + 1
                ENDIF

                ! ==--------------------------------------------------------------==

                DO ibuf = 1,2
                   IF (.NOT.no_more_states) THEN
                      IF (nbr_states(ibuf).NE.2) CYCLE
                   ENDIF
                   jb1 = id_states(1,ibuf)
                   jb2 = id_states(2,ibuf)
                   ! 
                   EHFX_loc = 0.0_real_8
                   IF (jb1.NE.0 .AND. jb2.NE.0) THEN
                      nbr_integrals=nbr_integrals+2
                      pfx1=pfl*f(ia)*f(jb1)
                      pfx2=pfl*f(ia)*f(jb2)
                      CALL hfxab2(EHFX_loc_1,EHFX_loc_2,&
                           GHFX_loc_1,GHFX_loc_2,&
                           pfx1,pfx2,psia,&
                           psib(1,ibuf),&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb1),c2_hfx(1,jb2))
                      EHFX_loc = EHFX_loc_1 + EHFX_loc_2
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,jb1)
                         int_j=MAX(ia,jb1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc_1
                         acc_vals( int_ij, 2 ) = GHFX_loc_1
                         new_vals( int_ij ) = 1
                         int_i=MIN(ia,jb2)
                         int_j=MAX(ia,jb2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc_2
                         acc_vals( int_ij, 2 ) = GHFX_loc_2
                         new_vals( int_ij ) = 1
                      ENDIF
                   ELSEIF (jb1.NE.0 .AND. jb2.EQ.0) THEN
                      ! Terms IA*IB1
                      nbr_integrals=nbr_integrals+1
                      pfx=pfl*f(ia)*f(jb1)
                      CALL hfxab(EHFX_loc,GHFX_loc,&
                           pfx,psia,psib(1,ibuf),1,&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb1))
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,jb1)
                         int_j=MAX(ia,jb1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc
                         acc_vals( int_ij, 2 ) = GHFX_loc
                         new_vals( int_ij ) = 1
                      ENDIF
                   ELSEIF (jb1.EQ.0 .AND. jb2.NE.0) THEN
                      ! Terms IA*IB2
                      nbr_integrals=nbr_integrals+1
                      pfx=pfl*f(ia)*f(jb2)
                      CALL hfxab(EHFX_loc,GHFX_loc,pfx,&
                           psia,psib(1,ibuf),2,&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb2))
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,jb2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc
                         acc_vals( int_ij, 2 ) = GHFX_loc
                         new_vals( int_ij ) = 1
                      ENDIF
                   ENDIF
                   ehfx=ehfx+EHFX_loc

                   ! reset the buffer
                   id_states(:,ibuf) = 0
                   nbr_states(ibuf) = 0
                ENDDO     ! ibuf
             ENDDO           ! IB
             CALL tihalt(procedureN//'_inner',isub3)
             ! ==--------------------------------------------------------------==
101          CONTINUE
          ENDDO                 ! IA
       ENDDO                   ! ispin
       CALL tihalt(procedureN//'_outer',isub2)
       ! ==--------------------------------------------------------------==
    ELSE
       ! HARTREE
       DO ia=1,nstate
          IF (f(ia).GT.1.e-6_real_8) THEN
             CALL zeroing(psia)!,maxfftn)
             !ocl novrec
             DO ig=1,jgw
                psia(nzfs(ig))=c0(ig,ia)
                psia(inzs(ig))=CONJG(c0(ig,ia))
             ENDDO
             IF (geq0) psia(nzfs(1))=c0(1,ia)
             ! Transform the wavefunction to real space
             CALL invfftn(psia,.TRUE.,parai%allgrp)
             pfx=pfl*f(ia)*f(ia)
             CALL hfxaa(EHFX_loc,GHFX_loc,pfx,psia,vpotg,vpotr,&
                  psic,C2_hfx(1,ia))
             ehfx=ehfx+EHFX_loc
          ENDIF
       ENDDO
    ENDIF
    ! free some memory
    ! 
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (func1%mhfx.EQ.1) THEN
       DEALLOCATE(psib,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF

    t2=m_cputime()
    dt = t2 - t1

    ! 
    ! redistribute EHFX and C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       ! call mp_sync(CP_GRP)
       CALL tiset(procedureN//'_b',isub6)
       CALL mp_sum(ehfx,parai%cp_inter_grp)
       ! call cp_grp_redist_z(C2_hfx,NGW,NSTATE)
       CALL mp_sum(C2_hfx,ncpw%ngw*nstate,parai%cp_inter_grp)
       CALL tihalt(procedureN//'_b',isub6)
    ENDIF
    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)

    IF (hfxc3%twscr) THEN
       CALL mp_sum(acc_vals,SIZE(acc_vals),parai%cp_grp)
       CALL mp_sum(new_vals,SIZE(new_vals),parai%cp_inter_grp)

       bin_range=0.5_real_8
       bin_max=200.0_real_8
       ALLOCATE( max_vals( INT( bin_max / bin_range ) + 1, 2 ),&
            bin_vals( INT( bin_max / bin_range ) + 1 ), stat=ierr )
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)

       max_vals(:,:) = 0.0_real_8
       bin_vals(:) = 0
       n_int=0
       max_dab=0.0_real_8

       ! The following compiler has a problem with the reduction in the omp section
       ! IBM XL Fortran Advanced Edition for Blue Gene/P, V11.1
       ! Version: 11.01.0000.0011
       !$omp parallel do &
       !$omp             default(none) &
       !$omp             private(i,ibin,int_i,int_j,dab) &
       !$omp             shared(new_vals,acc_vals,bin_range) &
       !$omp             shared(int_vals) &
       !$omp             reduction(+:n_int,bin_vals) &
       !$omp             reduction(max:max_vals,max_dab)
       DO i=1,SIZE(int_vals,1)
          IF (new_vals(i)>0) THEN
             int_vals(i,1) = ABS(acc_vals(i,1))
             int_vals(i,2) = acc_vals(i,2)
             n_int=n_int+1
          ENDIF
          int_j=CEILING((-1.0_real_8+SQRT(1.0_real_8+8.0_real_8*REAL(i,kind=real_8)))/2.0_real_8)
          int_i=i-int_j*(int_j-1)/2
          CALL get_wannier_separation(int_i,int_j,dab)
          max_dab=MAX(max_dab,dab)
          ibin=INT(dab/bin_range)+1
          max_vals(ibin,1)=MAX(max_vals(ibin,1),int_vals(i,1))
          max_vals(ibin,2)=MAX(max_vals(ibin,2),int_vals(i,2))
          bin_vals(ibin)=bin_vals(ibin)+1
       ENDDO
       !$omp end parallel do

       ! find new radius      
       old_DWFC=hfxc4%dwfc
       IF (hfxc3%twfc.AND.hfxc3%twft) THEN
          ibin=INT(hfxc4%dwfc/bin_range)+1
          IF (max_vals(ibin  ,1)>hfxc4%dwf_integral_thresh)&
               hfxc4%dwfc=hfxc4%dwfc+bin_range
          IF (max_vals(ibin-1,1)<hfxc4%dwf_integral_thresh.AND.&
               max_vals(ibin-2,1)<hfxc4%dwf_integral_thresh)&
               hfxc4%dwfc=hfxc4%dwfc-bin_range
       ENDIF
       IF (paral%io_parent.AND..FALSE.) THEN
          WRITE(6,*) 'we have ',n_int,' intergrals (re)set'
          WRITE(6,*) 'new DWFC=',hfxc4%dwfc,' old DWFC=',old_DWFC
          DO i=1,INT(max_dab/bin_range)+1
             WRITE(6,'(2F6.1,2E8.2,I8)') (i-1)*bin_range,&
                  i*bin_range,&
                  max_vals(i,:),bin_vals(i)
          ENDDO
       ENDIF

       IF (paral%io_parent.AND..FALSE.) THEN
          DO i=1,SIZE(int_vals,1)
             int_j=CEILING((-1.0_real_8+SQRT(1.0_real_8+8.0_real_8*REAL(i,kind=real_8)))/2.0_real_8)
             int_i=i-int_j*(int_j-1)/2
             CALL get_wannier_separation(int_i,int_j,dab)
             WRITE(6,'(A,I0,5E12.4)') 'INT',int_count,dab,&
                  wcentx(4,int_i),wcentx(4,int_j),&
                  int_vals(i,:)
          ENDDO
       ENDIF

       DEALLOCATE( max_vals,bin_vals,new_vals,acc_vals,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)

    ENDIF

    ! ENERGY
    ehfx=ehfx*parm%omega
    CALL mp_sum(ehfx,parai%allgrp)
    ! POTENTIAL
    DO ia=1,nstate
       vhfx=vhfx+dotp(ncpw%ngw,c0(:,ia),c2(:,ia))
    ENDDO
    CALL mp_sum(vhfx,parai%allgrp)
    DEALLOCATE(psi_in_core,C2_hfx,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)

    ! 
    ! some infos
    ! 
    t2_tot=m_cputime()
    dt_tot = t2_tot - t1_tot

    IF (PRINT_GROUP_INFOS) THEN
       IF (parai%me.EQ.0) THEN
          WRITE(6,'(1X,6(A,I0),A,F0.2,A,F0.2)')&
               procedureN//'| group ',parai%cp_inter_me,&
               ' computed ',nbr_integrals,&
               ' integrals, precomputed ',nbr_rwfn_precomputed,&
               ', recomputed ',nbr_rwfn_recomputed,&
               ' and skipped ',&
               nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
               ', t_loc ',dt,&
               ', t_per_1k_ints ',1.0e3_real_8*dt/(REAL(nbr_integrals,kind=real_8)+1.0e-6_real_8) ! to avoid NANs
          CALL m_flush(6)
       ENDIF

       IF (parai%cp_nogrp.GT.1) THEN
          CALL mp_sum(nbr_integrals,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_1,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_2,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_precomputed,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_recomputed,parai%cp_inter_grp)
          CALL mp_max(dt_tot,parai%cp_inter_grp)
          IF (parai%cp_me.EQ.0) THEN
             WRITE(6,'(1X,5(A,I0),A,F0.2)')&
                  procedureN//'| all the groups computed ',&
                  nbr_integrals,' integrals, precomputed ',&
                  nbr_rwfn_precomputed,&
                  ', recomputed ',nbr_rwfn_recomputed,&
                  ' and skipped ',&
                  nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
                  ', t_tot ',dt_tot
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfx_old
  ! ==================================================================
  SUBROUTINE check_occupation(f,nstate,tlsd)
    ! ==================================================================
    ! == Check the occupation number                                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    LOGICAL                                  :: tlsd

    REAL(real_8), PARAMETER                  :: one = 1.0_real_8, &
                                                two = 2.0_real_8 , &
                                                zero = 0.0_real_8

    INTEGER                                  :: i
    REAL(real_8)                             :: ff

! ==--------------------------------------------------------------==

    IF (tlsd) THEN
       !$omp parallel do default(none) private(I,FF) shared(F,NSTATE)
       DO i=1,nstate
          ff=ABS(f(i)-one)
          IF (ff.GT.zero.AND.ABS(f(i)).GT.zero) THEN
             CALL stopgm('HFX','ONLY OCCUPATIONS (0,1) ALLOWED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
       !$omp end parallel do
    ELSE
       !$omp parallel do default(none) private(I,FF) shared(F,NSTATE)
       DO i=1,nstate
          ff=ABS(f(i)-two)
          IF (ff.GT.zero .AND. ABS(f(i)).GT.zero) THEN
             CALL stopgm('HFX','ONLY OCCUPATIONS (0,2) ALLOWED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
       !$omp end parallel do
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE check_occupation

  ! ==================================================================
  SUBROUTINE get_wannier_separation(icenter,jcenter,dist)
    ! ==================================================================
    ! == Get the distance between two Wannier centers                 ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: icenter, jcenter
    REAL(real_8)                             :: dist

    REAL(real_8)                             :: coordab(3),coordab_(3)

! we should protect here
! IF(ICENTER.GT.NSTATE.OR.ICENTER.LT.1.OR.
! &   JCENTER.GT.NSTATE.OR.JCENTER.LT.1)
! &   CALL STOPGM('GET_WANNIER_SEPARATION',' something wrong here...')

    coordab_(1)=wcentx(1,icenter)-wcentx(1,jcenter)
    coordab_(2)=wcentx(2,icenter)-wcentx(2,jcenter)
    coordab_(3)=wcentx(3,icenter)-wcentx(3,jcenter)
    CALL pbc(coordab_(1),coordab_(2),coordab_(3),coordab(1),coordab(2),coordab(3),1,parm%apbc,parm%ibrav)
    dist=SQRT(coordab(1)**2+coordab(2)**2+coordab(3)**2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_wannier_separation

  ! ==================================================================
  SUBROUTINE hfxab(ehfx,ghfx,pf,psia,psib,iran,vpotg,vpotr,psic,&
       c2a,c2b)
    ! ==================================================================
    REAL(real_8)                             :: ehfx, ghfx, pf
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(llr1)
    COMPLEX(real_8)                          :: psic(:), c2a(jgw), c2b(jgw)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    ghfx=0.0_real_8
    IF (iran.EQ.1) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG) &
    !$omp  reduction(+:EHFX)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
       ehfx=ehfx+REAL(4._real_8*vpotg(ig)*CONJG(psic(nzff(ig))))
    ENDDO
    IF (geq0) ehfx=ehfx-REAL(2._real_8*vpotg(1)*CONJG(psic(nzff(1))))
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    IF (iran.EQ.1) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*(REAL(psia(ir))+uimag*REAL(psib(ir)))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*(REAL(psia(ir))+uimag*AIMAG(psib(ir)))
       ENDDO
    ENDIF
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx=ghfx+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b(ig)=c2b(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx=ghfx/2.0_real_8
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxab
  ! ==================================================================
  SUBROUTINE hfxab2(ehfx_1,ehfx_2,ghfx_1,ghfx_2,pf1,pf2,psia,psib,&
       vpotg,vpotr,psic,c2a,c2b1,c2b2)
    ! ==================================================================
    REAL(real_8)                             :: ehfx_1, ehfx_2, ghfx_1, &
                                                ghfx_2, pf1, pf2
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1), &
                                                vpotg(jhg,2)
    REAL(real_8)                             :: vpotr(llr1,2)
    COMPLEX(real_8)                          :: psic(:), c2a(jgw), c2b1(jgw), &
                                                c2b2(jgw)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab2'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    CALL tiset(procedureN,isub)
    ehfx_1=0.0_real_8
    ehfx_2=0.0_real_8
    ghfx_1=0.0_real_8
    ghfx_2=0.0_real_8
    DO ir=1,llr1
       psic(ir)=REAL(psia(ir))*psib(ir)
    ENDDO
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) &
    !$omp  reduction(+:EHFX_1,EHFX_2)
    DO ig=1,jhg
       fp=psic(nzff(ig))+psic(inzf(ig))
       fm=psic(nzff(ig))-psic(inzf(ig))
       vpotg(ig,1)=-pf1*scgx(ig)*0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       vpotg(ig,2)=-pf2*scgx(ig)*0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
       ehfx_1=ehfx_1+REAL(2._real_8*vpotg(ig,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_2=ehfx_2+REAL(2._real_8*vpotg(ig,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDDO
    IF (geq0) THEN
       fp=psic(nzff(1))+psic(inzf(1))
       fm=psic(nzff(1))-psic(inzf(1))
       ehfx_1=ehfx_1-REAL(vpotg(1,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_2=ehfx_2-REAL(vpotg(1,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDIF
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig,1)+uimag*vpotg(ig,2)
       psic(inzf(ig))=CONJG(vpotg(ig,1))+uimag*CONJG(vpotg(ig,2))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1,1)+uimag*vpotg(1,2)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir,1)=REAL(psic(ir))
       vpotr(ir,2)=AIMAG(psic(ir))
    ENDDO
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir,1)*(REAL(psia(ir))+uimag*REAL(psib(ir)))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX_1)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx_1=ghfx_1+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b1(ig)=c2b1(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx_1=ghfx_1/2.0_real_8
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir,2)*(REAL(psia(ir))+uimag*AIMAG(psib(ir)))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX_2)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx_2=ghfx_2+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b2(ig)=c2b2(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx_2=ghfx_2/2.0_real_8
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxab2
  ! ==================================================================
  SUBROUTINE hfxaa(ehfx,ghfx,pf,psia,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: ehfx, ghfx, pf
    COMPLEX(real_8)                          :: psia(:), vpotg(*)
    REAL(real_8)                             :: vpotr(*)
    COMPLEX(real_8)                          :: psic(:), c2a(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxaa'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    ghfx=0.0_real_8
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=REAL(psia(ir))*REAL(psia(ir))
    ENDDO
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG,FP) &
    !$omp  reduction(+:EHFX)
    DO ig=1,jhg
       fp=psic(nzff(ig))
       vpotg(ig)=-pf*scgx(ig)*fp
       ehfx=ehfx+REAL(2._real_8*vpotg(ig)*CONJG(fp))
    ENDDO
    IF (geq0) ehfx=ehfx-REAL(vpotg(1)*CONJG(psic(nzff(1))))
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir)*REAL(psia(ir))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx=ghfx+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
       c2a(ig)=c2a(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxaa
  ! ==================================================================
  SUBROUTINE hfxpsi_old(c0,cpsi,c2,f,sign,psia,nstate,norb)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE                                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8)                             :: f(:), sign
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate, norb
    COMPLEX(real_8)                          :: c2(ncpw%ngw,norb), &
                                                cpsi(ncpw%ngw,norb)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxpsi_old'

    COMPLEX(real_8), ALLOCATABLE             :: psib(:), psic(:), vpotg(:)
    INTEGER                                  :: i, ia, ib, ib1, ib2, ierr, &
                                                ig, isub
    REAL(real_8)                             :: pfl, pfx
    REAL(real_8), ALLOCATABLE                :: vpotr(:)

    IF (func1%mhfx.EQ.0) RETURN
    IF (func1%mhfx.EQ.2) CALL stopgm('HFXPSI_OLD',&
         'HARTREE METHOD NOT POSSIBLE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tiset('    HFXPSI_OLD',isub)
    IF (tkpts%tkpnt) CALL stopgm('HFXPSI_OLD','K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm('HFXPSI_OLD','NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm('HFXPSI_OLD','NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm('HFXPSI_OLD',&
         'TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO i=1,nstate
       IF (f(1).NE.f(i).AND. f(i).GT.1.e-6_real_8)&
            CALL stopgm('HFXPSI_OLD','OCCUPATION NUMBERS CHANGE',& 
            __LINE__,__FILE__)
    ENDDO
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    pfl=0.5_real_8
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    pfl=pfl*sign
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! HARTREE-FOCK
    ALLOCATE(psib(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ia=1,norb
       CALL zeroing(psia)!,maxfftn)
       !ocl novrec
       DO ig=1,jgw
          psia(nzfs(ig))=cpsi(ig,ia)
          psia(inzs(ig))=CONJG(cpsi(ig,ia))
       ENDDO
       IF (geq0) psia(nzfs(1))=cpsi(1,ia)
       ! Transform the wavefunction to real space
       CALL invfftn(psia,.TRUE.,parai%allgrp)
       DO ib=1,nstate,2
          ib1=ib
          ib2=ib+1
          IF (f(ib1).LT.1.e-6_real_8) ib1=0
          IF (ib2.GT.nstate) THEN
             ib2=0
          ELSEIF (f(ib2).LT.1.e-6_real_8) THEN
             ib2=0
          ENDIF
          IF (ib1.NE.0 .OR. ib2.NE.0) THEN
             IF (ib1.EQ.0) THEN
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=uimag*c0(ig,ib2)
                   psib(inzs(ig))=uimag*CONJG(c0(ig,ib2))
                ENDDO
                IF (geq0) psib(nzfs(1))=uimag*c0(1,ib2)
             ELSEIF (ib2.EQ.0) THEN
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=c0(ig,ib1)
                   psib(inzs(ig))=CONJG(c0(ig,ib1))
                ENDDO
                IF (geq0) psib(nzfs(1))=c0(1,ib1)
             ELSE
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=c0(ig,ib1)+uimag*c0(ig,ib2)
                   psib(inzs(ig))=CONJG(c0(ig,ib1))+&
                        uimag*CONJG(c0(ig,ib2))
                ENDDO
                IF (geq0) psib(nzfs(1))=c0(1,ib1)+uimag*c0(1,ib2)
             ENDIF
             ! Transform the wavefunction to real space
             CALL invfftn(psib,.TRUE.,parai%allgrp)
             ! Terms IA*IB1
             IF (ib1.NE.0) THEN
                pfx=pfl*f(ib1)
                CALL hfxpb(pfx,psia,psib,1,vpotg,vpotr,psic,c2(1,ia))
             ENDIF
             ! Terms IA*IB2
             IF (ib2.NE.0) THEN
                pfx=pfl*f(ib2)
                CALL hfxpb(pfx,psia,psib,2,vpotg,vpotr,psic,c2(1,ia))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(psib,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    HFXPSI_OLD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxpsi_old
  ! ==================================================================
  SUBROUTINE hfxpb(pf,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(llr1)
    COMPLEX(real_8)                          :: psic(llr1), c2a(jgw)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir

    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
    ENDDO
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)+0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxpb
  ! ==================================================================
  SUBROUTINE hfxrpa_old(c0,c1,c2,psia,nstate,tcis)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE CONTRIBUTION TO RPA                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    LOGICAL                                  :: tcis

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxrpa_old'

    COMPLEX(real_8), ALLOCATABLE             :: psib(:), psic(:), vpotg(:)
    INTEGER                                  :: ia, ib, ib1, ib2, ierr, ig, &
                                                isub, nlower, nupper
    REAL(real_8)                             :: pfl
    REAL(real_8), ALLOCATABLE                :: vpotr(:)

! Variables
! (NHG)
! ==--------------------------------------------------------------==

    IF (func1%mhfx.EQ.0) RETURN
    IF (func1%mhfx.EQ.2) CALL stopgm('HFXRPA','HARTREE METHOD NOT POSSIBLE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tiset('    HFXRPA',isub)
    IF (tkpts%tkpnt) CALL stopgm('HFXRPA','K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm('HFXRPA','NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm('HFXRPA','NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm('HFXRPA','TASK GROUPS NOT IMPLEMENTED'&
         ,& 
         __LINE__,__FILE__)
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    pfl=0.5_real_8
    IF (cntl%tlsd) pfl=1.0_real_8
    pfl=pfl*2.0_real_8
    IF (.NOT.tcis) pfl=2._real_8*pfl
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! HARTREE-FOCK
    ALLOCATE(psib(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ia=1,nstate
       CALL zeroing(psia)!,maxfftn)
       !ocl novrec
       DO ig=1,jgw
          psia(nzfs(ig))=c1(ig,ia)+uimag*c0(ig,ia)
          psia(inzs(ig))=CONJG(c1(ig,ia))+&
               uimag*CONJG(c0(ig,ia))
       ENDDO
       IF (geq0) psia(nzfs(1))=c1(1,ia)+uimag*c0(1,ia)
       ! Transform the wavefunctions to real space
       CALL invfftn(psia,.TRUE.,parai%allgrp)
       IF (cntl%tlsd) THEN
          IF (ia.LE.spin_mod%nsup) THEN
             nupper=spin_mod%nsup
             nlower=1
          ELSE
             nupper=nstate
             nlower=spin_mod%nsup+1
          ENDIF
       ELSE
          nupper=nstate
          nlower=1
       ENDIF
       DO ib=nlower,nupper,2
          CALL zeroing(psib)!,maxfftn)
          ib1=ib
          ib2=ib+1
          IF (ib2.GT.nupper) ib2=0
          IF (ib2.EQ.0) THEN
             !ocl novrec
             DO ig=1,jgw
                psib(nzfs(ig))=c0(ig,ib1)
                psib(inzs(ig))=CONJG(c0(ig,ib1))
             ENDDO
             IF (geq0) psib(nzfs(1))=c0(1,ib1)
          ELSE
             !ocl novrec
             DO ig=1,jgw
                psib(nzfs(ig))=c0(ig,ib1)+uimag*c0(ig,ib2)
                psib(inzs(ig))=CONJG(c0(ig,ib1))+&
                     uimag*CONJG(c0(ig,ib2))
             ENDDO
             IF (geq0) psib(nzfs(1))=c0(1,ib1)+uimag*c0(1,ib2)
          ENDIF
          ! Transform the wavefunction to real space
          CALL invfftn(psib,.TRUE.,parai%allgrp)
          ! 
          CALL hfxrpav(pfl,psia,psib,1,vpotg,vpotr,psic,c2(1,ib1))
          ! 
          IF (ib2.NE.0) THEN
             CALL hfxrpav(pfl,psia,psib,2,vpotg,vpotr,psic,c2(1,ib2))
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(psib,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    HFXRPA',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxrpa_old
  ! ==================================================================
  SUBROUTINE hfxrpav(pf,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(llr1)
    COMPLEX(real_8)                          :: psic(llr1), c2a(jgw)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir

    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
    ENDDO
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    DO ir=1,llr1
       psic(ir)=vpotr(ir)*REAL(psia(ir))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)-0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxrpav
  ! ==================================================================
  INTEGER FUNCTION dminloc(n,x,incx)
    IMPLICIT NONE
    INTEGER :: n,incx
    REAL(real_8) :: x(*)
    INTEGER :: i,ix
    REAL(real_8) :: dmin
    dminloc = 0
    IF (n.LT.1 .OR. incx.LE.0) RETURN
    ix = 1
    dmin = ABS(x(ix))
    dminloc = ix
    ix = ix + incx
    DO i = 2,n
       IF (ABS(x(ix)).LT.dmin) THEN
          dminloc = i
          dmin = ABS(x(ix))
       ENDIF
       ix = ix + incx
    ENDDO
  END FUNCTION dminloc

  INTEGER FUNCTION dmaxloc(n,x,incx)
    IMPLICIT NONE
    INTEGER :: n,incx
    REAL(real_8) :: x(*)
    INTEGER :: i,ix
    REAL(real_8) :: dmax
    dmaxloc = 0
    IF (n.LT.1 .OR. incx.LE.0) RETURN
    ix = 1
    dmax = ABS(x(ix))
    dmaxloc = ix
    ix = ix + incx
    DO i = 2,n
       IF (ABS(x(ix)).GT.dmax) THEN
          dmaxloc = i
          dmax = ABS(x(ix))
       ENDIF
       ix = ix + incx
    ENDDO
  END FUNCTION dmaxloc

END MODULE hfx_utils
