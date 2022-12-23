#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#endif
!     ==================================================================
      SUBROUTINE mm_quantum_topo(NSPq,NRPT,NSOLV)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
!     ----------------------------------------------------------------==
!cmb - Revised on 29 January 2008 after J. Blumberger's request
!cmb - Revised on 25 October 2016
! GROMOSDIP
      use prmem_utils, only: prmem
      use system
      use parac
      use mm_input
      use forcematch
      use mm_parallel
      use forcematch_utils, only : fm_setup_grm_equiv
      use forcematch_kfit_utils, only : fm_setup_topology

      implicit none
      include 'Gromos/toposz.h'
      include 'Gromos/topoar.h'
      integer NRPT,NSOLV,NSPq
!
      logical, allocatable ::  mm_q(:)
      logical, allocatable ::  mm_q_nb(:)
      integer mm_bond_0
      integer mm_bat_0
      integer mm_died_0
      integer mm_imp_0
      logical writeout

      integer m1,m2,m3,m4,i,is
      integer ia,ja,ic,icc
      integer MAXINE
!cmb      parameter (MAXINE=100) ! too small
      parameter (MAXINE=2000)
      integer, allocatable :: saved(:)
      integer, allocatable :: ij_excl(:,:)
      integer, allocatable :: ij_excl14(:,:)
      integer kne_l,kne14_l,n14

      ALLOCATE(saved(MAXNRP))
      ALLOCATE(ij_excl(MAXNRP,MAXINE))
      ALLOCATE(ij_excl14(MAXNRP,MAXINE))
      ALLOCATE(mm_q(NRPT+NSOLV))
      ALLOCATE(mm_q_nb(NRPT+NSOLV))
      IF(paral%io_parent) CALL PRMEM('MM_QM_TOPO')

      writeout=paral%parent

      DO i=1,NRPT+NSOLV
        mm_q(i)=.false.
        mm_q_nb(i)=.false.
        do is=1,NSPq
          do ia=1,gqmmm%gr_nasp(is)
            if(gqmmm%gr_atom(is,ia).eq.i)then
              if(i.gt.NRPT) call stopgm('mm_quantum_topo',
     &         'solvent moleules can not be quantum',
     & __LINE__,"__file__")
              mm_q(i)=.true.
              mm_q_nb(i)=.true.
              if(gqmmm%is_dummy(is))then
                mm_q_nb(i)=.false.
                write(6,*)'DUMMY FOUND',i
                if (CNTL%FMATCH) call stopgm('mm_quantum_topo',
     &         'Force-Matching not working with dummy atoms',
     & __LINE__,"__file__")
              endif
            endif
          enddo
        enddo
       ENDDO
!
C_FM[
C we need something more sophisticated for forcematching
      if (CNTL%FMATCH) then
        call fm_setup_grm_equiv
        call fm_setup_topology(NRPT, NSOLV, mm_q_nb) 
        goto 100
      endif
C_FM]
      NBTY=NBTY+1
      NTTY=NTTY+1
      NQTY=NQTY+1
      NPTY=NPTY+1
      if(NBTY.gt.MAXNBT) call stopgm('mm_quantum_topo','NBTY.gt.MAXNBT',
     & __LINE__,"__file__")
      if(NTTY.gt.MAXTTY) call stopgm('mm_quantum_topo','NTTY.gt.MAXTTY',
     & __LINE__,"__file__")
      if(NQTY.gt.MAXQTY) call stopgm('mm_quantum_topo','NQTY.gt.MAXQTY',
     & __LINE__,"__file__")
      if(NPTY.gt.MAXPTY) call stopgm('mm_quantum_topo','NPTY.gt.MAXPTY',
     & __LINE__,"__file__")

      mm_bond_0=NBTY
      cb(NBTY)=0.d0
      b0(NBTY)=1.d0

      mm_bat_0=NTTY
      CT(NTTY)=0.d0
      T0(NTTY)=1.d0

      mm_imp_0=NQTY
      CQ(NQTY)=0.d0
      Q0(NQTY)=1.d0

      mm_died_0=NPTY
      CP(NPTY)=0.d0
      PD(NPTY)=1.d0
      NP(NPTY)=2

      if(writeout)write(6,*) 
     c'THE FOLLOWING INTERACTIONS ARE EXCLUDED FROM THE MM HAMILTONIAN:'
!
! bonds involving hydrogen
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'BONDS INVOLVING HYDROGEN:'
      do i=1,NBONH
        m1=IBH(i)
        m2=JBH(i)
        if(mm_q(M1) .and. mm_q(M2)) then
          ICBH(i)=mm_bond_0
          if(writeout)write(6,62)m1,m2
        endif
      enddo

!
! other bonds 
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'BONDS NOT INVOLVING HYDROGEN:'
      do i=1,NBON
        m1=IB(i)
        m2=JB(i)
        if(mm_q(M1) .and. mm_q(M2)) then
          ICB(i)=mm_bond_0
          if(writeout)write(6,62)m1,m2
        endif
      enddo

!
! angles involving hydrogen
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'ANGLES INVOLVING HYDROGEN:'
      do i=1,NTHEH
        m1 = ITH(i)
        m2 = JTH(i)
        m3 = KTH(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3)) then
          ICTH(i)=mm_bat_0
          if(writeout)write(6,63)m1,m2,m3
        endif
      enddo

!
! other angles 
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'ANGLES NOT INVOLVING HYDROGEN:'
      do i=1,NTHE
        m1 = IT(i)
        m2 = JT(i)
        m3 = KT(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3)) then
          ICT(i)=mm_bat_0
          if(writeout)write(6,63)m1,m2,m3
        endif
      enddo

!
! proper dihedrals invoving hydrogen
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'PROPER DIHEDRALS INVOLVING HYDROGEN:'
      do i=1,NPHIH
        m1 = IPH(i)
        m2 = JPH(i)
        m3 = KPH(i)
        m4 = LPH(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3) .and. mm_q(M4)) then
          ICPH(i)=mm_died_0
          if(writeout)write(6,64)m1,m2,m3,m4
        endif
      enddo

!
! proper dihedrals
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'PROPER DIHEDRALS NOT INVOLVING HYDROGEN:'
      do i=1,NPHI
        m1 = IP(i)
        m2 = JP(i)
        m3 = KP(i)
        m4 = LP(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3) .and. mm_q(M4)) then
          ICP(i)=mm_died_0
          if(writeout)write(6,64)m1,m2,m3,m4
        endif
      enddo

!
! improper dihedrals invoving hydrogen
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'IMPROPER DIHEDRALS INVOLVING HYDROGEN:'
      do i=1,NQHIH
        m1 = IQH(i)
        m2 = JQH(i)
        m3 = KQH(i)
        m4 = LQH(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3) .and. mm_q(M4)) then
          ICQH(i)=mm_imp_0
          if(writeout)write(6,64)m1,m2,m3,m4
        endif
      enddo

!
! improper dihedrals
!
      if(writeout)write(6,*)
      if(writeout)write(6,*)'IMPROPER DIHEDRALS NOT INVOLVING HYDROGEN:'
      do i=1,NQHI
        m1 = IQ(i)
        m2 = JQ(i)
        m3 = KQ(i)
        m4 = LQ(i)
        if(mm_q(M1) .and. mm_q(M2) .and. mm_q(M3) .and. mm_q(M4)) then
          ICQ(i)=mm_imp_0
          if(writeout)write(6,64)m1,m2,m3,m4
        endif
      enddo

C_FM[
 100  CONTINUE
      if (CNTL%FMATCH) then
C       save exclusions for later
        allocate(fm_save_INE(MAXNRP))
        allocate(fm_save_KNE(MAXNRP))
        allocate(fm_save_JSNE(MAXAEX))
        allocate(fm_save_INE14(MAXNRP))
        allocate(fm_save_KNE14(MAXNRP))
        allocate(fm_save_JSNE14(MXEX14))
C
        fm_save_INE = INE
        fm_save_KNE = KNE
        fm_save_JSNE = JSNE
        fm_save_INE14 = INE14
        fm_save_KNE14 = KNE14
        fm_save_JSNE14 = JSNE14
      endif
C_FM]

! new exclusion list
!$OMP parallel do private(ia,ic) __COLLAPSE2
      DO ia=1,MAXNRP
        DO ic=1,MAXINE
          ij_excl(ia,ic)=MAXNRP+1
          ij_excl14(ia,ic)=MAXNRP+1
        ENDDO
      ENDDO

      kne_l=0
      kne14_l=0
      IF(writeout) then
        write(6,*)
        write(6,*)
     &       'NON BONDED INTERACTIONS EXCLUDED FROM THE MM HAMITONIAN:'
      ENDIF
      DO ia=1,NRP
        IF(ine(ia).gt.MAXINE)
     &       call stopgm('mm_quantum_topo','MAXINE too small',
     & __LINE__,"__file__")
        DO ja=1,ine(ia)
          ij_excl(ia,ja)=jsne(kne_l+ja)
        ENDDO
        DO ja=1,ine14(ia)
          ij_excl14(ia,ja)=jsne14(kne14_l+ja)
        ENDDO
        kne_l=kne_l+ine(ia)
        kne14_l=kne14_l+ine14(ia)
        IF(mm_q_nb(ia))THEN 
          DO ja=ia+1,NRP
            IF(mm_q_nb(ja))THEN !both atoms are quantum

              IF(ine(ia)+1.gt.MAXINE)
     &             call stopgm('mm_excl', 'MAXINE too small',
     & __LINE__,"__file__")

              IF(ine(ia).eq.0)THEN
                ij_excl(ia,1)=ja
                ine(ia)=1 
              ELSEIF(ja.lt.ij_excl(ia,1))THEN
                DO ic=ine(ia),1,-1
                  ij_excl(ia,ic+1)=ij_excl(ia,ic) 
                ENDDO
                ij_excl(ia,1)=ja
                ine(ia)=ine(ia)+1 
              ELSEIF (ja.gt.ij_excl(ia,ine(ia))) THEN
                ij_excl(ia,ine(ia)+1)=ja
                ine(ia)=ine(ia)+1 
              ELSE
                DO ic=1,ine(ia)
                  IF(ij_excl(ia,ic).lt.ja
     &                 .and.ij_excl(ia,ic+1).gt.ja)THEN
!$OMP parallel do private(icc) shared(ia)
#ifdef __SR11000
*poption parallel, tlocal(icc)
#endif 
                    DO icc=ic+1,ine(ia)
                      saved(icc)=ij_excl(ia,icc)
                    ENDDO
                    ij_excl(ia,ic+1)=ja
!$OMP parallel do private(icc) shared(ia)
#ifdef __SR11000
*poption parallel, tlocal(icc)
#endif
                    DO icc=ic+1,ine(ia)
                      ij_excl(ia,icc+1)=saved(icc)
                    ENDDO
                    ine(ia)=ine(ia)+1 
                    GOTO 111
                  ENDIF
                ENDDO
              ENDIF
 111          CONTINUE
c     check if ia and ja are 14 excluded, if so delete the 14 exclusion in the 
c     14 exclusion list, and reorder list
              DO ic=1,ine14(ia)
                IF (ij_excl14(ia,ic).eq.ja) THEN
                  IF (ine14(ia).eq.1) THEN
                    ij_excl14(ia,ic)=MAXNRP+1
                    ine14(ia)=0
                  ELSEIF (ic.eq.ine14(ia)) THEN
                    ij_excl14(ia,ic)=MAXNRP+1
                    ine14(ia)=ine14(ia)-1
                  ELSE
                    DO icc=ic+1,ine14(ia)
                      ij_excl14(ia,icc-1)=ij_excl14(ia,icc) 
                    ENDDO
                    ij_excl14(ia,icc)=MAXNRP+1
                    ine14(ia)=ine14(ia)-1
                    GOTO 112 
                  ENDIF
                ENDIF
              ENDDO     
 112          CONTINUE
            ENDIF
          ENDDO
        ENDIF
      ENDDO

c Redefine jsne and jscne14

      kne_l=0
      kne14_l=0
      do ia=1,NRP
        kne(ia)=kne_l
        kne14(ia)=kne14_l         
        do ja=1,ine(ia)
          jsne(kne_l+ja)=ij_excl(ia,ja)
        enddo
        do ja=1,ine14(ia)
          jsne14(kne14_l+ja)=ij_excl14(ia,ja)
        enddo
        kne_l=kne_l+ine(ia)
        kne14_l=kne14_l+ine14(ia)
      enddo
      if(kne_l.gt.MAXAEX)then
        write(6,*)"MAXAEX should be at least ",kne_l
        call stopgm('mm_quantum_topo',' MAXAEX too small',
     & __LINE__,"__file__")
      endif

      if(gparal%mmparent) THEN
         write(6,*) 'NEW EXCLUSION LIST FOR QM ATOMS'
         kne_l=0
         do ia=1,NRP
           if(mm_q_nb(ia))
     &         write(6,*)ia,ine(ia),(jsne(kne_l+ja),ja=1,ine(ia))
           kne_l=kne_l+ine(ia)
         enddo
         write(6,*) 'NEW 14 EXCLUSION LIST FOR QM ATOMS'
         kne_l=0
         do ia=1,NRP
           if(mm_q_nb(ia))
     &         write(6,*)ia,ine14(ia),(jsne14(kne_l+ja),ja=1,ine14(ia))
           kne_l=kne_l+ine14(ia)
         enddo
      ENDIF

! QM-MM couples excluded from the electrostatic interaction
      if(iqmmm%coupl_model.gt.0)then
! gromos exclusions
        if(excl_comm_l%gromos_excl)then
          excl_comm%nce=0
          do M1=1,NRP
            do i=KNE(M1)+1,KNE(M1)+INE(M1)
              M2=JSNE(i)
              if(mm_q(M1) .and. .not. mm_q(M2)) then
                if(DABS(CG(M2)).gt.1.d-6)then
                   excl_comm%nce=excl_comm%nce+1
                   excl_comm%atom_qm_excl(excl_comm%nce)=M1
                   excl_comm%atom_mm_excl(excl_comm%nce)=M2
                endif
              else if(.not. mm_q(M1) .and. mm_q(M2)) then
                if(DABS(CG(M1)).gt.1.d-6)then
                   excl_comm%nce=excl_comm%nce+1
                   excl_comm%atom_qm_excl(excl_comm%nce)=M2
                   excl_comm%atom_mm_excl(excl_comm%nce)=M1
                endif
              endif
            enddo
          enddo
          if(excl_comm%nce.gt.NCe_max)
     &         call stopgm('MM_QUANTUM_TOPO','too many exclusions',
     & __LINE__,"__file__")
        endif
        if(writeout)then
          if(excl_comm%nce.gt.0)then
            write(6,*)
            write(6,*)
     c'QM MM COUPLES EXCLUDED FROM THE ELECTROSTATIC INTERACTION'
            write(6,*)' QM ATOM   MM ATOM '
            do i=1,excl_comm%nce
              write(6,62)excl_comm%atom_qm_excl(i),
     &              excl_comm%atom_mm_excl(i)
            enddo
          else
            write(6,*)
            write(6,*)
     c'NO ELECTROSTATIC EXCLUSION BETWEEN QM AND MM ATOMS'
          endif
        endif
      else
        excl_comm%nce=0
      endif

      DEALLOCATE(ij_excl14,ij_excl,saved,mm_q_nb,mm_q)
62    format(2(i8,2x))
63    format(3(i8,2x))
64    format(4(i8,2x))
      RETURN
      END
