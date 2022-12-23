C    ==================================================================
      SUBROUTINE mm_elstat(TAU,mm_FION,mm_rhoe,
     .                     elstat_eint,MAXNAT_ell,update_files)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ----------------------------------------------------------------== 
C     This routine adds an external potential due to the MM charges
C     on the real space grid of the QM box. It also calculates the forces
c     acting on MM atoms due to the (smoothed) electron distribution
cmb - Revised on 28 February 2014
C     ==--------------------------------------------------------------==
      use system
      use parac
      use machine, only: m_walltime
      use mp_interface, only: mp_bcast
      use prmem_utils, only: prmem
      use epot_types
      use cnst
      use cell
      use adat
      use ions
      use ropt
      use store_types
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
!
      use mm_dimmod
      use mm_input
      use mm_parallel
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen
C_FM[
      use forcematch
C_FM]
      use efld
      IMPLICIT NONE

      INTEGER MAXNAT_ell
      REAL*8  elstat_eint,t1,t2
      REAL*8  TAU(3,maxsys%nax,*),mm_FION(3,maxsys%nax,*)
      REAL*8  mm_rhoe(fpar%kr1,fpar%kr2s,*)
      REAL*8, ALLOCATABLE ::  EXTF_l(:,:,:)
cmb -> temporary solution to be eliminated at the time
cmb -> of SAMPLE keyword and timestep printout rewrite
      LOGICAL update_files ! cmb - temporary solution
!
!  save
!
      integer nn_sr,nn_esp,nn_lr,npr, nn_excl_mm
      integer, pointer :: nn_atom(:)
      integer, pointer :: nn_sp(:)
      integer, pointer :: gr_n_print(:)
      real*8, pointer ::  q_mm(:)
      SAVE nn_atom,nn_sp,nn_sr,nn_esp,nn_lr,nn_excl_mm
      SAVE npr,gr_n_print
      SAVE q_mm
!
!   Locals
!
      real*8 , allocatable ::  q_partition(:)
      real*8 , allocatable ::  q_esp(:)
      real*8 , allocatable ::  F_ch(:,:)
      real*8 , allocatable ::  F_atdens_qm(:,:)
      real*8 , allocatable ::  F_el(:,:)
      real*8 , allocatable ::  F_atdens_mm(:,:)
      real*8 , allocatable ::  TAU_nn(:,:)
      real*8 , allocatable ::  v_mm(:)
      real*8 , allocatable ::  q_mm_eff(:)
      real*8 , allocatable ::  q_mm_eff_h(:)
      real*8 , allocatable ::  R_smear(:)
      real*8 , allocatable ::  r_t(:)
      real*8 , allocatable ::  r_n(:)
      real*8 , allocatable ::  v_mm_ion(:)
      real*8 , allocatable ::  v_mm_atdens(:)

C     If gaussian_ions=.true. the ions are represented as
C     gaussians of radius Rg= rg=mm_RAGGIO(IS) otherwise like point charges
      logical gaussian_ions
      parameter(gaussian_ions=.false.)

C     If atdens_fit=.true. the ESP charges are fitted the the fields due to 
C     rho_electron-rho_atom, where rho_atom is the free atom density superposition.
C     Otherwise, the ESPs are fitted to rho_electron-rho_ion
      logical atdens_fit,ferror
      parameter(atdens_fit=.false.)

      logical  long_range_esp
      parameter(long_range_esp=.true.)
!  energies
      real*8 eint_ch,eint_dp,eint_qp,e_coulomb_lr,elstat_lr
      real*8 e_excl,e_charge_rest
      real*8 elstat_sr,e_coulomb_sr,e_ion_sr,e_ele_sr
      real*8 e_lr_esp,e_lr_esp_mm
! long range
      real*8 ch_mm,dp_mm(3),qp_mm(3,3)
      real*8 ch_ai,dp_ai(3),qp_ai(3,3)
!
      integer iaad,iad,isd,ixd
      real*8 e_1,e_2,DDX,f_1,f_2
      real*8 , allocatable :: TAU_mi(:,:,:)
      real*8 xm(3),tt1,tt2,t_sr,t_lc
      real*8 smear_fac,rsq,xtm,ytm,ztm 
      REAL*8 dum1,dum2,dum3,dum4
      integer ix,id,iq
      integer ia,is,iaa,im,i,fo_first
      CHARACTER*2 chemd
      character*3 cd
      INTEGER it
      LOGICAL first_time,do_force
      SAVE first_time,it,fo_first
      DATA first_time /.true./
      DATA it /0/
      DATA fo_first /FO_VMARK/
      CHARACTER DCDTYPE*4,DCDTITLE*80
      DATA DCDTYPE,DCDTITLE /'CORD',
     & 'CHARMM format DCD INTERACTING TRAJECTORY.     '/
      REAL*4  DCDDELTA
      REAL*8  DCDCELL(6)
      REAL*4, ALLOCATABLE :: DCDX(:),DCDY(:),DCDZ(:)
      INTEGER :: err,ierr
      character(*),parameter :: procedureN='mm_elstat'
      tt1=m_walltime()
!
      if (first_time) then
        first_time = .false.
        allocate(nn_atom(mmdim%NATm))
        allocate(nn_sp(mmdim%NATm))
        allocate(gr_n_print(mne%maxnat_el))
        allocate(q_mm(mmdim%NATm))
!$OMP parallel do private(i)
#ifdef _vpp_
!OCL NOALIAS
#endif
        DO i=1,mmdim%NATm
          nn_atom(i)=0
          nn_sp(i)=0
          q_mm(i)=0.d0
        ENDDO 
        IF(paral%io_parent) CALL PRMEM(' MM_ELSTAT')
      end if

      ALLOCATE(q_partition(mmdim%NATq),q_esp(mmdim%NATq),
     &     F_ch(3,mmdim%NATq),
     &  F_atdens_qm(3,mmdim%NATq),F_el(3,MAXNAT_ell),
     &  F_atdens_mm(3,MAXNAT_ell),TAU_nn(3,MAXNAT_ell),v_mm(MAXNAT_ell),
     &  q_mm_eff(MAXNAT_ell),q_mm_eff_h(MAXNAT_ell),R_smear(MAXNAT_ell),
     &  r_t(MAXNAT_ell),r_n(MAXNAT_ell),v_mm_ion(MAXNAT_ell),
     &  v_mm_atdens(MAXNAT_ell),TAU_mi(3,maxsys%nax,maxsys%nsx))

      IF(paral%IO_PARENT.AND.lqmmm%qmmm_time) CALL FILEOPEN(113,
     &     'TIMING',FO_APP,FERROR)

! electrostatic coupling is computed with positions that are the minimal 
! with respect to the center of mass of the quantum system.
      XM(1)=0.5D0*CELL_COM%CELLDM(1)
      XM(2)=XM(1)*CELL_COM%CELLDM(2)
      XM(3)=XM(1)*CELL_COM%CELLDM(3)
      call mm_min_im(TAU,TAU_mi,maxsys%nax,maxsys%nsx,XM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(mod(it,iqmmm%nt_nlist).eq.0) THEN
        CALL mm_nlist_el(TAU_mi,nn_sr,nn_esp,nn_lr,nn_atom,nn_sp,
     .                   nn_excl_mm, q_mm)
        IF(paral%io_parent) THEN
 200      FORMAT (' MM_ELSTAT: ',A,/,1X,A,I8,',  ',A,I8,/,
     .                               1X,A,I8,',  ',A,I8,/)
          write(6,200)'ELECTROSTATIC INTERACTION NN LIST UPDATED',
     .        'NN ATOMS = ',nn_sr,'ESP COUPLED ATOMS = ',nn_esp-nn_sr,
     .        'LR ATOMS = ',nn_lr-nn_esp,'EMM COUPLED ATOMS = ',
     .                                                nn_excl_mm-nn_lr
          CALL PRMEM(' MM_ELSTAT')
C_FM[
          if(cntl%fmatch .and. nn_esp.gt.nn_sr) then
            write(6,'(A)') '  fm======================================'
            write(6,'(A)') '  fm  WARNING !!!'
            write(6,'(A)') '  fm  There are ESP coupled atoms while'
            write(6,'(A)') '  fm  doing force-matching!'
            write(6,'(A)') '  fm  Be sure you really want that!'
            write(6,'(A)') '  fm======================================'
          endif
C_FM]
        ENDIF
      endif
!!!! DEBUG
!          ddx=0.0000001
!          do iaad = 24,28
!            iad=nn_atom(iaad)
!            isd=nn_sp(iaad)
!            do ixd=1,1
!              do id=1,2
!                CALL AZZERO(mm_FION,3*maxsys%nax*maxsys%nsx) 
!                if(id.eq.1)TAU_mi(ixd,iad,isd)=TAU_mi(ixd,iad,isd)-ddx
!                if(id.eq.2)TAU_mi(ixd,iad,isd)=TAU_mi(ixd,iad,isd)+ddx
!!!! END DEBUG

      xtm=0.D0
      ytm=0.D0
      ztm=0.D0
!$OMP parallel do private(is,ia) reduction(+:xtm,ytm,ztm)
#ifdef __SR11000
*poption parallel, tlocal(is,ia), psum(xtm,ytm,ztm)
#endif 
      DO is=1,mmdim%NSPq
        DO ia=1,NAq(is)
          xtm=xtm+TAU_mi(1,ia,is)
          ytm=ytm+TAU_mi(2,ia,is)
          ztm=ztm+TAU_mi(3,ia,is)
        ENDDO
      ENDDO
      xm(1)=xtm/DBLE(mmdim%NATq)
      xm(2)=ytm/DBLE(mmdim%NATq)
      xm(3)=ztm/DBLE(mmdim%NATq)

      IF(paral%qmnode) THEN
        t1=m_walltime()
        CALL  mm_multipole(ch_ai,dp_ai,qp_ai,xm,TAU_mi,mm_rhoe) 
        t2=m_walltime()
        IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
          WRITE(113,*) "QMNODE: mm_multipole ",(t2-t1)/1000.d0
        ENDIF
      ENDIF
      CALL Mp_BCAST(ch_ai,parai%qmmmsource,parai%qmmmgrp) 
      CALL Mp_BCAST(dp_ai,3,parai%qmmmsource,parai%qmmmgrp) 
      CALL Mp_BCAST(qp_ai,9,parai%qmmmsource,parai%qmmmgrp) 

      if(paral%qmnode)then

! smearing function
        smear_fac=1.0d0 ! change this to affect R_smear
        DO iaa =1,nn_sr
          ia=nn_atom(iaa)
          is=nn_sp(iaa)
          R_smear(iaa)=smear_fac*COVRAD(ions0%IATYP(IS))*FBOHR
          IF(R_smear(iaa).le.0.D0) THEN
            call stopgm('MM_ELSTAT','wrong COVRAD',
     & __LINE__,"__file__")
          ENDIF
          r_n(iaa)=R_smear(iaa)**n_sm%n_smear
          r_t(iaa)=r_n(iaa)*R_smear(iaa)
        ENDDO
      endif

      IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
        WRITE(113,*) "SYNCHRONIZE"
      ENDIF
      CALL Mp_BCAST(ions0%zv,mmdim%NSPq,parai%qmmmsource,
     &     parai%qmmmgrp) 
      CALL Mp_BCAST(r_t,nn_sr,parai%qmmmsource,parai%qmmmgrp) 
      CALL Mp_BCAST(r_n,nn_sr,parai%qmmmsource,parai%qmmmgrp) 

      allocate(EXTF_l(fpar%kr1,fpar%kr2s,fpar%kr3s), stat=err)
      if (err /= 0) call stopgm('MM_ELSTAT', 'Error allocating EXTF_l',
     & __LINE__,"__file__")

      if(paral%qmnode)then
! calculate the q_mm_eff
        do_force=.false.
        t1=m_walltime()
        call azzero(q_partition,mmdim%NATq)
        call azzero(q_mm_eff,MAXNAT_ell)
        if (.not.excl_comm_l%excl_mech) then
        call mm_excl(nn_sr,nn_esp,do_force,F_el,F_atdens_mm,
     c            F_atdens_qm,F_ch,
     c            q_mm,q_mm_eff,q_mm_eff_h,
     c            nn_atom,nn_sp,long_range_esp,e_lr_esp,e_lr_esp_mm,
     c            TAU_mi,mm_FION,v_mm,e_charge_rest,
     c            r_t,r_n,q_partition,q_esp,e_excl,EXTF,EXTF_l)
        t2=m_walltime()
        IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
          WRITE(113,*) "QMNODE: mm_excl ",(t2-t1)/1000.d0
        ENDIF
        endif

!$OMP parallel do private(iaa,ia,is)
#ifdef __SR11000
*poption parallel, tlocal(iaa,ia,is)
#endif 
        DO iaa =1,nn_sr
          ia=nn_atom(iaa)
          is=nn_sp(iaa)
          TAU_nn(1,iaa)=TAU_mi(1,ia,is)
          TAU_nn(2,iaa)=TAU_mi(2,ia,is)
          TAU_nn(3,iaa)=TAU_mi(3,ia,is)
        ENDDO

! calculate potential on the grid due to q_mm_eff and the forces on the MM atoms
! in the nn_list (F_el) due to the electronic charge

        CALL AZZERO(EXTF,fpar%kr1*fpar%kr2s*fpar%kr3s)
        CALL AZZERO(EXTF_l,fpar%kr1*fpar%kr2s*fpar%kr3s)

        t1=m_walltime()
        call mm_elstat_sr(TAU_nn,F_el,v_mm,mm_rhoe, 
     c                     q_mm_eff,q_mm,nn_sp,mmdim%NSPq,
     c                     e_ele_sr,nn_sr,r_t,r_n,EXTF,EXTF_l)  
        t2=m_walltime()
        if(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
          WRITE(113,*) "QMNODE: mm_elstat_sr ",(t2-t1)/1000.d0
        ENDIF
        if (.not.excl_comm_l%excl_mech) then
        t1=m_walltime()
        call mm_charges(mm_rhoe,TAU_mi,q_partition,q_mm_eff_h,F_ch,
     c          nn_sr,nn_atom,nn_sp,EXTF_l)
        t2=m_walltime()
        IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
          WRITE(113,*) "QMNODE: mm_charges ",(t2-t1)/1000.d0
        ENDIF
        endif
        t1=m_walltime()
        call mm_elstat_atdens(TAU_nn,TAU_mi,v_mm_atdens,nn_sr,
     c         R_smear,nn_atom,nn_sp,F_atdens_mm,F_atdens_qm,q_mm_eff,
     c         r_t,r_n,atdens_fit,q_mm)
        t2=m_walltime()
        IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
          WRITE(113,*)"QMNODE: mm_elstat_atdens ",(t2-t1)/1000.d0
        ENDIF
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
        DO i=1,nn_sr
          v_mm(i)=v_mm(i)+v_mm_atdens(i)
        ENDDO

! Forces on mm atoms due to electronic charge distribution are added to the
! forces on the MM atoms in the nn list

        IF(paral%io_parent) THEN
!$OMP parallel do private(iaa,ia,is)
#ifdef __SR11000
*poption parallel, tlocal(iaa,ia,is)
#endif 
          DO iaa =1,nn_sr
            ia=nn_atom(iaa)
            is=nn_sp(iaa)
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+q_mm(iaa)*F_el(1,iaa)
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+q_mm(iaa)*F_el(2,iaa)
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+q_mm(iaa)*F_el(3,iaa)
          ENDDO
        ENDIF
      endif

      if(gparal%mmparent)then
        t1=m_walltime()
        call mm_short_range_classic(TAU_mi,q_mm,
     c             nn_sr,nn_atom,nn_sp,mm_FION,
     c             e_coulomb_sr,e_ion_sr,v_mm_ion,r_t,r_n,
     c             gaussian_ions)
        t2=m_walltime()
        t_sr=(t2-t1)/1000.d0

        t1=m_walltime()
        if(iqmmm%coupl_model.eq.2)then
          call mm_long_range_classic(ch_ai,dp_ai,qp_ai,xm,TAU_mi,
     c             nn_esp,nn_lr,nn_atom,nn_sp,mm_FION,
     c             ch_mm,dp_mm,qp_mm,
     c             eint_ch,eint_qp,eint_dp,e_coulomb_lr) 
        endif
        t2=m_walltime()
        t_lc=(t2-t1)/1000.d0
      endif

      call Mp_BCAST(t_sr,gparai%mmsource,parai%qmmmgrp) 
      call Mp_BCAST(t_lc,gparai%mmsource,parai%qmmmgrp) 
      call Mp_BCAST(v_mm_ion,nn_sr,gparai%mmsource,parai%qmmmgrp) 
      call Mp_BCAST(e_coulomb_sr,gparai%mmsource,parai%qmmmgrp) 
      call Mp_BCAST(e_ion_sr,gparai%mmsource,parai%qmmmgrp) 
      if(paral%qmnode)then
        elstat_sr=e_ele_sr+e_ion_sr
        elstat_eint=elstat_sr-e_coulomb_sr
      endif
      IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
        WRITE(113,*)"gparal%mmnode: mm_short_classic",t_sr
        WRITE(113,*)"gparal%mmnode: mm_long_classic",t_lc
        WRITE(113,*) "SYNCHRONIZE"
      ENDIF

      if(iqmmm%coupl_model.eq.2)then
        call Mp_BCAST(eint_ch,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(eint_qp,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(eint_dp,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(e_coulomb_lr,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(ch_mm,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(dp_mm,3,gparai%mmsource,parai%qmmmgrp) 
        call Mp_BCAST(qp_mm,9,gparai%mmsource,parai%qmmmgrp) 
        if(paral%qmnode)then
          elstat_lr=eint_ch+eint_dp+eint_qp
          elstat_eint=elstat_eint+elstat_lr-e_coulomb_lr
        endif
! sum the contributions to mm_fion from gparal%mmnode and one of the qmnodes (parent)
!
        IF(paral%qmnode)THEN
          t1=m_walltime()
          CALL mm_long_range_potential(ch_mm,dp_mm,qp_mm,xm,EXTF)
          t2=m_walltime()
          IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
            WRITE(113,*)"QMNODE: mm_long_range_potential"
     .                  ,(t2-t1)/1000.d0
          ENDIF
        ENDIF
      endif

! The ESP charges are calculated and the excluded coupples contribution to the forces is added.
!

      if(paral%qmnode)then
        do_force=.true.
        t1=m_walltime()
        if (excl_comm_l%excl_mech) then
           call mm_excl_mech(nn_sr,nn_esp,do_force,
     .                           nn_lr,nn_excl_mm,
     .               nn_atom,nn_sp,long_range_esp,
     .               e_lr_esp,e_lr_esp_mm,
     .               TAU_mi,mm_FION
     .               )
        else
        call mm_excl(nn_sr,nn_esp,do_force,F_el,F_atdens_mm,
     c            F_atdens_qm,F_ch,
     c            q_mm,q_mm_eff,q_mm_eff_h,
     c            nn_atom,nn_sp,long_range_esp,e_lr_esp,e_lr_esp_mm,
     c            TAU_mi,mm_FION,v_mm,e_charge_rest,
     c            r_t,r_n,q_partition,q_esp,e_excl,EXTF,EXTF_l)
        endif
        t2=m_walltime()
        IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN ! cmb - bugfix
          WRITE(113,*) "QMNODE: mm_excl_2 ",(t2-t1)/1000.0d0
        ENDIF
        IF(paral%io_parent) THEN
          IF(excl_comm%NCe.gt.0) elstat_eint=elstat_eint+e_excl
          IF(cgrest_i%n_cg.gt.0) elstat_eint=elstat_eint+e_charge_rest
          IF(long_range_esp)elstat_eint=elstat_eint+e_lr_esp-e_lr_esp_mm
        ENDIF
      endif

      deallocate(EXTF_l)
!!!!

!!!! DEBUG
!                if(id.eq.1)TAU_mi(ixd,iad,isd)=TAU_mi(ixd,iad,isd)+ddx
!                if(id.eq.2)TAU_mi(ixd,iad,isd)=TAU_mi(ixd,iad,isd)-ddx
!                if(id.eq.1)e_1=elstat_eint
!                if(id.eq.2)e_2=elstat_eint
!                if(id.eq.1)f_1=mm_FION(ixd,iad,isd)
!                if(id.eq.2)f_2=mm_FION(ixd,iad,isd)
!              enddo
!            write(6,*)iaad,-(e_2-e_1)/(2.*ddx),(f_1+f_2)*0.5
!            write(16,*)iaad,-(e_2-e_1)/(2.*ddx),(f_1+f_2)*0.5
!            enddo
!          enddo
!          stop
!!!! END DEBUG
! write out
      if(paral%qmnode)then
        if(iteropt%NFI.eq.wp_i%nfi_wp .and. wp_l%write_potential) then
            call mm_write_potential(wp_i%n_stride,EXTF,
     &          wp_c%cubename_pot)  
        endif
        if(iteropt%NFI.eq.wp_i%nfi_wd .and. wp_l%write_density) then
            call mm_write_potential(wp_i%n_stride,mm_rhoe,
     &          wp_c%cubename_dens)  
        endif
      endif
 
      tt2=m_walltime()
        
      IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
        WRITE(113,*) "TOTAL TIME: ",(tt2-tt1)/1000
      ENDIF
!     tt1=m_walltime()
      if(paral%io_parent)then
        if(iqmmm%nt_sample.ne.0) then
          CALL FILEOPEN(60,'EL_ENERGY',FO_APP,FERROR)
          CALL FILEOPEN(12,'CHJ',FO_APP,FERROR)
          CALL FILEOPEN(61,'ESP',FO_APP,FERROR)
          CALL FILEOPEN(50,'MULTIPOLE',FO_APP,FERROR)
          if(iteropt%NFI.gt.0.and.update_files)then
            if(iqmmm%coupl_model.eq.2)then
              write(60,160) iteropt%NFI,elstat_eint+
     c              eint_ch+eint_dp+eint_qp,
     c             elstat_sr,   ! interaction with the NN atom
     c             e_lr_esp,    ! with the ESP-coupled
     c             eint_ch+eint_dp+eint_qp, ! with the multipole
     c             e_coulomb_sr+e_coulomb_lr+e_lr_esp_mm !force field counterpart
            else
              write(60,160) iteropt%NFI,elstat_sr,e_coulomb_sr
            endif
            write(50,150) iteropt%NFI,(dp_ai(ix),ix=1,3),qp_ai(1,1),
     c           qp_ai(2,2),
     c           qp_ai(1,2),qp_ai(1,3),qp_ai(2,3)
            write(12,161) iteropt%NFI,(q_partition(ia),ia=1,mmdim%NATq)
            write(61,161) iteropt%NFI,(q_esp(i),i=1,mmdim%NATq)
          endif
          CALL FILECLOSE(61) 
          CALL FILECLOSE(12)
          CALL FILECLOSE(60)
          CALL FILECLOSE(50)
        endif    ! nt_sample.gt.0
C_FM[
C        if (FMATCH.and.CONVWF.ne.0) then
C          CALL FILEOPEN(12,fm_chj_ref_filename,FO_APP,FERROR)
C          if (FERROR) then
C            call stopgm('mm_elstat','Could not open file '
C     &                 // fm_chj_ref_filename,
C     & __LINE__,"__file__")
C          endif
C          write(12,161) fm_iframe,(q_partition(ia),ia=1,NATq)
C          CALL FILECLOSE(12)
C        endif
C_FM]
        if((iqmmm%nt_sample.ne.0).and.(mod(it,iqmmm%nt_NLIST).eq.0))then
          CALL FILEOPEN(71,'INTERACTING.pdb',FO_OLD,FERROR)
          IF(FERROR) call stopgm('MM_ELSTAT','NO INTERACTING.pdb',
     & __LINE__,"__file__")
          do i=1,mne%maxnat_el
cmb            read(71,*,end=111)cd,id,cd,cd,gr_n_print(i)
cmb -> a proper PDB file MUST have the following format or
cmb -> it is unreadable with the format above
            READ(71,171,end=111)id,chemd,cd,id,dum1,dum2,dum3,dum4,
     .                          gr_n_print(i)
          enddo
111       npr=i-1
          CALL FILECLOSE(71)
          CALL FILEOPEN(60,'TRAJECTORY_INTERACTING',
     &         FO_APP+FO_FIRST,FERROR)
          CALL FILECLOSE(60)
          IF(rout1%MOUT) THEN
            CALL FILEOPEN(70,'MOVIE_INTERACTING',
     &           FO_APP+FO_FIRST,FERROR)
            CALL FILECLOSE(70)
          ENDIF
          FO_FIRST=0
        endif
        if(iqmmm%nt_sample.gt.0) then
          if (mod(iteropt%NFI,iqmmm%nt_sample).eq.0 
     &          .and. (iteropt%NFI.ne.0) .and. update_files) then
            CALL FILEOPEN(60,'TRAJECTORY_INTERACTING',FO_APP,FERROR)
            do i=1,npr
              ia=cpat(gr_n_print(i))
              is=cpsp(gr_n_print(i))
              write(60,106) iteropt%NFI,(TAU(ix,ia,is),ix=1,3)
            enddo
            CALL FILECLOSE(60)
            IF(rout1%MOUT.AND.(mod(iteropt%NFI,cnti%IMOVIE).EQ.0)) THEN
              CALL FILEOPEN(70,'MOVIE_INTERACTING',FO_APP,FERROR)
              do i=1,npr
                ia=cpat(gr_n_print(i))
                is=cpsp(gr_n_print(i))
                write(70,107) (TAU(ix,ia,is),ix=1,3),ions0%IATYP(IS)
              enddo
              CALL FILECLOSE(70)
            ENDIF
            CALL FILEOPEN(60,'TRAJECTORY_PIP',FO_APP,FERROR)
            write(60,*)nn_sr,iteropt%NFI
            do i=1,nn_sr
              cd='MM '
              if(nn_sp(i).le.mmdim%NSPq)cd='QM '
              write(60,108) iteropt%NFI,(TAU_nn(ix,i),ix=1,3),cd,
     &             q_mm(i),v_mm(i),R_smear(i),v_mm_ion(i),
     &             v_mm_atdens(i),
     &             (q_mm(i)*(F_el(ix,i)-F_atdens_mm(ix,i)),ix=1,3)
            enddo
            CALL FILECLOSE(60) 
          endif
        endif                   ! nt_sample.gt.0

C_FM[
C
C       this is still in parent!
        if (cntl%fmatch .and. ropt_mod%CONVWF .and. fm_fit_charges) then
          call fileopen(62,fm_pip_ref_filename,FO_APP,FERROR)
C       write only numbers, makes parsing easier when reading
C          write(62,*)'no NN atoms: ',nn_sr,'  FM frame index: ',
          write(62,*) nn_sr, fm_frame_indx
          do i=1,nn_sr
            cd='MM '
            if(nn_sp(i).le.mmdim%NSPq)cd='QM '
            write(62,181) (TAU_nn(ix,i),ix=1,3),cd,
     &             q_mm(i),v_mm(i),R_smear(i),v_mm_ion(i),
     &             v_mm_atdens(i),
C write out Field, instead of Force on NN atoms. Otherwise, if q_mm(i)
C is zero, we lose the information about the field
C     &             (q_mm(i)*(F_el(ix,i)-F_atdens_mm(ix,i)),ix=1,3)
     &             ((F_el(ix,i)-F_atdens_mm(ix,i)),ix=1,3)
          enddo
          call fileclose(62)
C write out the Hirschfeld charges on the QM atoms               
          CALL FILEOPEN(64,fm_chj_ref_filename,FO_APP,FERROR)
C       also here, write only numbers for easier parsing
C          write(64,*) 'FM frame index: ',fm_frame_indx
          write(64,*) fm_frame_indx
          write(64,182) (q_partition(ia),ia=1,mmdim%NATq)
               call fileclose(64)
        end if

181     FORMAT(1X,3f16.10,1x,a3,1x,8f16.10)
182     FORMAT(1x,130(f10.5,1x))
C_FM]

C  detect whether we need to update the .dcd
        ferror=.true.
        if(iqmmm%nt_sample.lt.0) ferror=(mod(iteropt%NFI-1,
     &       -iqmmm%nt_sample)
     &       .ne.0)
        if(iqmmm%nt_sample.gt.0.and.ferror) 
     &       ferror=(mod(iteropt%NFI-1,iqmmm%nt_sample).ne.0)
        if(lqmmm%qmmm_dcd.and.(.not.ferror).and.update_files) then
C     CHECK IF WE NEED TO WRITE A DCD HEADER
          FERROR=.TRUE.
          INQUIRE(FILE='TRAJ_INT.dcd',EXIST=FERROR)
          IF(.NOT.FERROR) THEN
            CALL FILEOPEN(35,'TRAJ_INT.dcd',FO_APP+FO_UFO,FERROR)
            DCDDELTA=cntr%DELT_IONS
            WRITE(35) DCDTYPE,0,iteropt%NFI,cnti%NTRAJ,0,0,0,0,0,0,
     &           DCDDELTA,1,
     &                0,0,0,0,0,0,0,0,24
            WRITE(35) 2,DCDTITLE,DCDTITLE
            WRITE(35) npr
            CALL FILECLOSE(35)
          ENDIF
C..Write dcd file output
          CALL FILEOPEN(35,'TRAJ_INT.dcd',FO_APP+FO_UFO,FERROR)
C..FIXME: we don't update the number of frames field in the header (yet).
          allocate(DCDX((npr+1)/2),STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
          allocate(DCDY((npr+1)/2),STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
          allocate(DCDZ((npr+1)/2),STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
!$OMP parallel do private(i,ia,is)
          DO i=1,npr
              ia=cpat(gr_n_print(i))
              is=cpsp(gr_n_print(i))
              DCDX(I)=(TAU(1,ia,is)-clsaabox%MM_C_TRANS(1))/FBOHR
              DCDY(I)=(TAU(2,ia,is)-clsaabox%MM_C_TRANS(2))/FBOHR
              DCDZ(I)=(TAU(3,ia,is)-clsaabox%MM_C_TRANS(3))/FBOHR
          ENDDO
C       store unitcell. use classical box for QM/MM
          DCDCELL(1)=clsaabox%BOX_au(1)/FBOHR
          DCDCELL(3)=clsaabox%BOX_au(2)/FBOHR
          DCDCELL(6)=clsaabox%BOX_au(3)/FBOHR
          DCDCELL(2)=0.0d0
          DCDCELL(4)=0.0d0
          DCDCELL(5)=0.0d0
          WRITE(35) (DCDCELL(I),I=1,6)
          WRITE(35) (DCDX(I),I=1,npr)
          WRITE(35) (DCDY(I),I=1,npr)
          WRITE(35) (DCDZ(I),I=1,npr)
          deallocate(DCDX,STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')
          deallocate(DCDY,STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')
          deallocate(DCDZ,STAT=ierr)
          if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')
          CALL FILECLOSE(35)
        endif
      endif

      it=it+1
!     tt2=m_walltime()
        
!     if (parent) write(113,*) "WRITING OUT",(tt2-tt1)/1000

      IF(paral%IO_PARENT.and.lqmmm%qmmm_time) THEN
        CALL FILECLOSE(113)
      ENDIF

150   FORMAT(i6,3(1x,f13.8),3x,5(f15.8,1x))
160   FORMAT(i6,6(1x,f15.8))
106   FORMAT(i7,3f16.6)
108   FORMAT(i7,3f16.10,1x,a3,1x,8f16.10) 
107   FORMAT(3f14.4,i4,'   0')
161   FORMAT(i6,1x,130(f10.5,1x))
cmb - Readable PDB files (Amber, gOpenMol, VMD, ...) must have this format:        
171   FORMAT(5x,I6,1x,A2,3x,A3,I6,4x,3F8.3,F6.2,I6)
!171   FORMAT('ATOM ',I6,1x,A2,3x,A3,I6,4x,3F8.3,F6.2,I6)
!cmb - deallocate memory in the proper last-in/first-out order
!cmb - to avoid corruption
      DEALLOCATE(q_partition)
      DEALLOCATE(TAU_mi)
      DEALLOCATE(v_mm_atdens)
      DEALLOCATE(v_mm_ion)
      DEALLOCATE(r_n)
      DEALLOCATE(r_t)
      DEALLOCATE(R_smear)
      DEALLOCATE(q_mm_eff_h)
      DEALLOCATE(q_mm_eff)
      DEALLOCATE(v_mm)
      DEALLOCATE(TAU_nn)
      DEALLOCATE(F_atdens_mm)
      DEALLOCATE(F_el)
      DEALLOCATE(F_atdens_qm)
      DEALLOCATE(F_ch)
      DEALLOCATE(q_esp)

      RETURN
      END
