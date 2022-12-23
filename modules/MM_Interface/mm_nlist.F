C     ==================================================================
      SUBROUTINE mm_nlist_el(TAU,nn_sr,nn_esp,nn_lr,nn_atom,nn_sp,
     c                    nn_excl_mm,q_mm)
C     ==--------------------------------------------------------------==
cmb - Revised on 04 December 2012
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use timer, only: tiset, tihalt
      use system
      use parac
      use mm_dimmod
      use ions
      use cell
      use adat
      use mm_input
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen
      use forcematch
      implicit none
!
      real*8  TAU(3,maxsys%nax,maxsys%nsx)
      real*8  q_mm(*)
      integer nn_sr,nn_esp,nn_lr, nn_excl_mm
      integer nn_atom(*)
      integer nn_sp(*)
!  locals
      integer ig,n_ch
      integer ig0,igf,icg,ia,i,iaa,iss,is1
      integer iaq,isq,ipp
      logical ferror
      integer nn,is,ig1,ix,iw
      real*8 rc2,Rijsq,q_tot,delta_q,r2_min

      integer lab(mmdim%NATm)
      integer lab_DM, lab_QM, lab_NN, lab_ES, lab_MP, lab_EXCL_MM
      parameter (lab_DM=0, lab_QM=1, lab_NN=2, lab_ES=3, lab_MP=4,
     .           lab_EXCL_MM=5)
      integer inc_group,ISUB
      integer inc_FAR, inc_ESP, inc_BIG, inc_ALL, inc_QUA, inc_CUT
      parameter (inc_FAR=0, inc_ESP=1, inc_BIG=2, inc_ALL=3,
     &     inc_QUA=4, inc_CUT=5)
C      character*2 lab(NATm)
      character*2  aa2
      character*3  strp
C      character*3 inc_group
      real*8    q_prt,ch_tot,rc2_l
      real*8    rc2_small,rc2_esp

      CALL TISET('MM_NLIST_EL',ISUB)
      rc2=rqmmm%rcut_el**2
      rc2_small=rqmmm%rcut_small**2
      rc2_esp=rqmmm%rcut_esp**2

      if(paral%io_parent)then
        CALL FILEOPEN(71,'INTERACTING.pdb',FO_NEW,FERROR)
        if(FERROR)then
          CALL FILEOPEN(71,'INTERACTING_NEW.pdb',FO_DEF,FERROR)
        endif
      endif

!$OMP parallel do private(ig)
#ifdef __SR11000
*poption parallel, tlocal(ig)
#endif
      DO ig=1,mmdim%NATm
        lab(ig)=lab_DM
      ENDDO 

! loop over the charge groups
      do icg=1,NCAG_l
        ig0=1
        if(icg.ge.2) ig0=INC_l(icg-1)+1
        igf=INC_l(icg)
! double cutoff for non zero charge charge groups
        ch_tot=0.D0
        do ig=ig0,igf
          is=cpsp(ig)
          ia=cpat(ig)
          ch_tot=ch_tot+mm_charge(ia,is)
        enddo
        rc2_l=rc2
!       if(abs(ch_tot).gt.0.3) rc2_l=4.d0*rc2
! finds where is the group with respect to the quantum system       
        R2_min=1.0D+06
        DO ig=ig0,igf
          is=cpsp(ig)
          ia=cpat(ig)
          IF(is.le.mmdim%NSPq) lab(ig)=lab_QM
          DO iss=1,mmdim%NSPq
            DO iaa=1,NAq(iss)
              Rijsq=(TAU(1,ia,is)-TAU(1,iaa,iss))**2
     +             +(TAU(2,ia,is)-TAU(2,iaa,iss))**2
     +             +(TAU(3,ia,is)-TAU(3,iaa,iss))**2
              R2_min=min(R2_min,Rijsq)
            ENDDO
          ENDDO
        ENDDO
        inc_group=inc_FAR
        if(R2_min.lt.rc2_esp  ) inc_group=inc_ESP
        if(R2_min.lt.rc2_l    ) inc_group=inc_BIG
        if(R2_min.lt.rc2_small) inc_group=inc_ALL
        do ig=ig0,igf
          is=cpsp(ig)
          if(is.le.mmdim%NSPq) inc_group=inc_QUA
        enddo
        if(inc_group.eq.inc_QUA)then
          do ig=ig0,igf
            is=cpsp(ig)
            do ig1=ig0,igf
              is1=cpsp(ig1)
              if((is.gt.mmdim%NSPq .and. is1.le.mmdim%NSPq) .or. 
     c           (is1.gt.mmdim%NSPq .and. is.le.mmdim%NSPq))
     &             inc_group=inc_CUT
            enddo
          enddo
          if(inc_group.eq.inc_CUT .and. paral%io_parent) then 
            write(6,'(2A,I8,A,2I8)')' WARNING! ',
     &           'CUTTING THROUGH CHARGE GROUP',icg,
     &           '  ATOMS:',ig0,igf
          endif
        endif
!   sets the proper label for the atom, according to its group label and its charge
        do ig=ig0,igf
          is=cpsp(ig)
          ia=cpat(ig)
          q_prt=DABS(mm_charge(ia,is))
          if(inc_group.eq.inc_QUA                     )lab(ig)=lab_QM
          if(inc_group.eq.inc_CUT .and. q_prt.ge.1.D-8)lab(ig)=lab_NN
          if(inc_group.eq.inc_ALL .and. q_prt.ge.1.D-8)lab(ig)=lab_NN
          if(inc_group.eq.inc_BIG .and. q_prt.ge.1.D-8)lab(ig)=lab_ES
Cmb  The threshold below is set to be consistent with AMBER
          if(inc_group.eq.inc_BIG .and. q_prt.ge.0.089D0)lab(ig)=lab_NN
          if(inc_group.eq.inc_ESP .and. q_prt.ge.1.D-8)lab(ig)=lab_ES
          if(inc_group.eq.inc_FAR .and. q_prt.ge.1.D-8)lab(ig)=lab_MP
C_FM[
C       For force-matching, we want also atoms with charge = 0 in the NN
C       group (for the charge-fitting)
          if(CNTL%FMATCH) then
            if(inc_group.eq.inc_CUT) lab(ig) = lab_NN
            if(inc_group.eq.inc_ALL) lab(ig) = lab_NN
            if(inc_group.eq.inc_BIG) lab(ig) = lab_NN
          endif
C_FM]
        enddo
! the MM atoms involved in the exclusions are coupled either
! via the ESP COUPLING scheme or via a mechanical coupling (MM-CHARGES)
!cmb        do ig=1,NCe ! do not break loops with IF
          if (excl_comm_l%excl_mech) then
!$OMP parallel do private(ig)
             DO ig=1,excl_comm%NCe
               lab(excl_comm%atom_mm_excl(ig))=lab_EXCL_MM
             ENDDO
          else
!$OMP parallel do private(ig)
             DO ig=1,excl_comm%NCe
               lab(excl_comm%atom_mm_excl(ig))=lab_ES
             ENDDO
          endif
!cmb        enddo ! do not break loops with IF

C_FM[
C       For force-matching, some atoms are not included (those close to
C       capping hydrogens, for example)
        if (CNTL%FMATCH) then
!$OMP parallel do private(ig)
          do ig = 1, fm_nexclude
            lab(fm_exclude(ig)) = lab_DM
          enddo
        endif
C_FM]

! write INTERACTING.pdb
! FIXME: make sure, that this is a valid .pdb file or change the name
! AK 2005/04/09
        if(paral%io_parent)then
          if(inc_group.eq.inc_QUA .or. inc_group.eq.inc_CUT .or.
     c          inc_group.eq.inc_ALL .or. inc_group.eq.inc_BIG)then
            do ig=ig0,igf
              is=cpsp(ig)
              ia=cpat(ig)
              aa2=elem%el(ions0%iatyp(is))
              q_prt=mm_charge(ia,is)
              if(inc_group.eq.inc_QUA) q_prt=0.D0
              strp='PRT'
              if(is.le.mmdim%NSPq)strp='QUA'
              iw=1
              if(is.le.mmdim%NSPq)iw=2
              WRITE(71,171) ! cmb - True PDB format
     c        ig,aa2,strp,iw,(TAU(ix,ia,is)*0.529177248D0,ix=1,3)
     c        ,q_prt,NAT_cpmd(ig)
            enddo
          endif
        endif
      enddo

      if(paral%io_parent) call fileclose(71)

! costructs the nn list
      nn = 0
      do ig=1,mmdim%NATm
        if(lab(ig).eq.lab_NN .or. lab(ig).eq.lab_QM)then
          nn=nn+1
          is=cpsp(ig)
          ia=cpat(ig)
          nn_atom(nn)=ia
          nn_sp(nn)=is
          q_mm(nn)=mm_charge(ia,is)
C          if (parent) print*,' ig=',ig,' nn=',nn,' q=',mm_charge(ia,is)
          if(lab(ig).eq.lab_QM) q_mm(nn)=0.D0
          IF (nn.gt.mne%maxnat_el) THEN
            WRITE(6,*) 'WARNING: NN=',nn,' > MAXNN=',mne%maxnat_el
            call stopgm('MM_NLIST','MAXNN too small ',
     & __LINE__,"__file__")
          ENDIF
        endif
      enddo
      nn_sr=nn

      do ig=1,mmdim%NATm
        if(lab(ig).eq.lab_ES)then
          nn=nn+1
          is=cpsp(ig)
          ia=cpat(ig)
          nn_atom(nn)=ia
          nn_sp(nn)=is
          q_mm(nn)=mm_charge(ia,is)
        endif
      enddo
      nn_esp=nn

      do ig=1,mmdim%NATm
        if(lab(ig).eq.lab_MP)then
          nn=nn+1
          is=cpsp(ig)
          ia=cpat(ig)
          nn_atom(nn)=ia
          nn_sp(nn)=is
          q_mm(nn)=mm_charge(ia,is)
        endif
      enddo
      nn_lr=nn

!mcb exclusion list atoms via MM-charges
      do ig=1,mmdim%NATm
        if(lab(ig).eq.lab_EXCL_MM)then
          nn=nn+1
          is=cpsp(ig)
          ia=cpat(ig)
          nn_atom(nn)=ia
          nn_sp(nn)=is
          q_mm(nn)=mm_charge(ia,is)
        endif
      enddo
      nn_excl_mm=nn

! computes the total charge of the interacting system. If this charge is not integer, prints out
! a warning and distributes it on all the atoms

      q_tot=0.d0
      n_ch=0
      do i=1,nn_sr
        ia=nn_atom(i)
        is=nn_sp(i)
        if(DABS(q_mm(i)).lt.1.d-8)then
          q_prt=mm_charge(ia,is)
          q_tot=q_tot+q_prt
        else
          n_ch=n_ch+1
        endif
      enddo 

      delta_q=q_tot-anint(q_tot)

      if(DABS(delta_q).gt.1.d-6)then
        IF(paral%io_parent) THEN
          write(6,'(/,4X,A,/,4X,A,/,4X,A,F10.6,A,/,4X,A,/,4X,A,/)')
     &         ' !!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!',
     &         ' THE QM SYSTEM DOES NOT HAVE AN INTEGER CHARGE.',
     &         ' A COMPENSATING CHARGE OF ',delta_q,' HAS BEEN',
     &         ' DISTRIBUTED OVER THE NN ATOMS.',
     &         ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        ENDIF
        do i=1,nn_sr
          if(DABS(q_mm(i)).gt.1.d-8)then
            q_mm(i)=q_mm(i)-delta_q/DBLE(n_ch)
          endif
        enddo
      endif

#if 0
      if(excl_all)then
        write(6,*)'!!!!! EXCL_ALL IS SET TO .TRUE. !!!!!'
        excl_comm%NCe=0
        do nn=1,nn_sr
          ia=nn_atom(nn)
          is=nn_sp(nn)
          do isq=1,mmdim%NSPq
            do iaq=1,NAq(isq)
              excl_comm%NCe=excl_comm%NCe+1
              ipp=maxsys%nax*(is-1)+ia
              ig=gratom(ipp)
              excl_comm%atom_mm_excl(excl_comm%NCe)=ig
              ipp=maxsys%nax*(isq-1)+iaq
              ig=gratom(ipp)
              excl_comm%atom_qm_excl(excl_comm%NCe)=ig
            enddo
          enddo
        enddo
      endif
#endif

cmb - Readable PDB files (Amber, gOpenMol, VMD, ...) must have this format:        
171   FORMAT('ATOM ',I6,1x,A2,3x,A3,I6,4x,3F8.3,F6.2,I6)
      CALL TIHALT('MM_NLIST_EL',ISUB)
      return
      end
