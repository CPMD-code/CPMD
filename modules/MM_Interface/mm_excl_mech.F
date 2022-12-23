C     ==================================================================
      SUBROUTINE mm_excl_mech(nn_sr,nn_esp,do_force,
     .                           nn_lr,nn_excl_mm,
     .               nn_atom,nn_sp,long_range_esp,
     .               e_lr_esp,e_lr_esp_mm,
     .               TAU_mi,mm_FION
     .               )
C     ==--------------------------------------------------------------==
C     ==  The MM atoms involved in exclusions are coupled to the QM   ==
C     ==  system via the MM CHARGES (mechanical coupling)
C     ==--------------------------------------------------------------==
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use timer, only: tiset, tihalt
      use system
      use parac
      use mm_dimmod
      use mm_input
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      implicit none


C Arguments.
      integer nn_sr,nn_esp
      integer nn_lr,nn_excl_mm
      logical do_force
      integer nn_atom(*),nn_sp(*)

      real*8  e_lr_esp,e_lr_esp_mm
      logical long_range_esp
      real*8  TAU_mi(3,maxsys%nax,*),mm_FION(3,maxsys%nax,*)

      logical, allocatable :: excl_mm(:,:)
      integer, allocatable :: atom_qm_excl_l(:),atom_mm_excl_l(:)
      real*8, allocatable :: F_QM(:,:)
C local variables.
      real*8 xij(3),rmm(3),rij_l,qj,qj_t,qi
      real*8 Rijsq,ft,ft_x,ft_y,ft_z
      real*8 qj_mm
      integer atom_qm_cg_l(mmdim%NATq)
      integer NCe_l
      integer iaa,ia,is,ia1,is1,iq,im,ix,ie,i,iim
      integer iqqe,iqe,ime,j,k
      logical found
      integer ISUB 

      real*8 xij1,xij2,xij3

      if( (excl_comm%nce.gt.0 .or. (nn_lr+1.le.nn_excl_mm)) 
     &     .and. cgrest_i%n_cg.gt.0)
     .      call stopgm('MM_EXCL_MECH','incompatible options',
     & __LINE__,"__file__")

      CALL TISET('   MM_EXCL_MECH',ISUB)

      allocate(excl_mm(mmdim%NATq,nn_excl_mm-nn_lr))
      allocate(atom_qm_excl_l(NCe_max))
      allocate(atom_mm_excl_l(NCe_max))
      allocate(F_QM(3,mmdim%NATq))

!mcb excluded QM-MM pair information in new ordering ...
      iim=0
      DO j=nn_lr+1,nn_excl_mm  !mcb loop MM of EXCL_MM list
        iim=iim+1
        DO i=1,mmdim%NATq   !mcb loop QM
          excl_mm(i,iim)=.false.
        ENDDO
      ENDDO
      e_lr_esp=0.d0
      e_lr_esp_mm=0.d0

      NCe_l=0

!mcb identify QM-MM pair from read-in exclusion list
!mcb and store in new ordering
      DO ie=1,excl_comm%nce !mcb loop exclusion list from cpmd.in 
!      QM atom 
        iqe=excl_comm%atom_qm_excl(ie)
        if(cpsp(iqe).gt.mmdim%NSPq) then !mcb test if QM
           call stopgm('MM_EXCL_MECH','MM ATOM EXCLUDED AS QM ATOM',
     & __LINE__,"__file__")
        endif
!      MM atom 
        ime=excl_comm%atom_mm_excl(ie) !mcb MM atom
        ia1=cpat(ime)
        is1=cpsp(ime)
        if(is1.le.mmdim%NSPq) then ! test if MM
           call stopgm('MM_EXCL_MECH','QM ATOM EXCLUDED AS MM ATOM',
     & __LINE__,"__file__")
        endif 

!mcb   determine linear index for QM atom
        iq=0
        DO is=1,cpsp(iqe)-1
          iq=iq+NAq(is)
        ENDDO
        iqqe=iq+cpat(iqe)

!mcb   determine linear index for MM atom from EMM list
        iim=0
        DO iaa =nn_lr+1,nn_excl_mm
          iim=iim+1
          ia=nn_atom(iaa)
          is=nn_sp(iaa)
          if(ia.eq.ia1 .and. is.eq.is1)then
            found=.true.
!mcb       store in new ordering ...
            excl_mm(iqqe,iim)=.true.
          endif
        ENDDO
        IF(.not.found)THEN
         call stopgm('mm_excl_mech','MAYBE TOO SMALL ES CUTOFF',
     & __LINE__,"__file__")
        ENDIF
!mcb---------------------------------------
      ENDDO
          
!
      call azzero(F_QM,3*mmdim%NATq)

! long range
      if(long_range_esp .and. paral%io_parent)then
        iq=0
        e_lr_esp=0.d0
        e_lr_esp_mm=0.d0

!mcb   loop QM atoms ...
        do is1=1,mmdim%NSPq
          do ia1=1,NAq(is1)
            iq=iq+1
            rmm(1)=TAU_mi(1,ia1,is1)
            rmm(2)=TAU_mi(2,ia1,is1)
            rmm(3)=TAU_mi(3,ia1,is1)
            qj=mm_charge(ia1,is1)
!mcb       loop MM atoms of EMM list ...
            iim=0
            do im=nn_lr+1,nn_excl_mm
              iim=iim+1

!mcb         skip pairs from exclusion list ...
              if( .not.excl_mm(iq,iim) )  then

                 ia=nn_atom(im)
                 is=nn_sp(im)
                 qi=mm_charge(ia,is)
                 xij(1) = TAU_mi(1,ia,is) - rmm(1)
                 xij(2) = TAU_mi(2,ia,is) - rmm(2)
                 xij(3) = TAU_mi(3,ia,is) - rmm(3)
                 Rijsq=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
                 Rij_l=DSQRT(Rijsq)
                 e_lr_esp= e_lr_esp + qi*qj / Rij_l
                 ft=qi*qj/(Rijsq*Rij_l)
                 ft_x=ft*xij(1)
                 ft_y=ft*xij(2)
                 ft_z=ft*xij(3)
                 F_QM(1,iq)=F_QM(1,iq)-ft_x
                 F_QM(2,iq)=F_QM(2,iq)-ft_y
                 F_QM(3,iq)=F_QM(3,iq)-ft_z
                 mm_FION(1,ia,is)= mm_FION(1,ia,is)+ft_x
                 mm_FION(2,ia,is)= mm_FION(2,ia,is)+ft_y
                 mm_FION(3,ia,is)= mm_FION(3,ia,is)+ft_z
               endif 
            enddo ! loop MM
          enddo ! loop QM ia
        enddo ! loop QM is
!       if(parent)write(6,*)"LONG RANGE ESP ENERGY ",
!    c                        e_lr_esp_mm,e_lr_esp
      endif
      if(paral%io_parent) then
!mcb   store forces for QM atoms ...
        iaa=0
        DO is=1,mmdim%NSPq
          DO ia=1,NAq(is)
            iaa=iaa+1
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+F_QM(1,iaa)
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+F_QM(2,iaa)
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+F_QM(3,iaa)
          ENDDO
        ENDDO 
      endif
cMcB

      deallocate(excl_mm,atom_qm_excl_l,atom_mm_excl_l,F_QM)

      CALL TIHALT('   MM_EXCL_MECH',ISUB)
      END
