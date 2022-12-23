#if defined(__BG)
#define __VECTOR
#endif
#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#endif
!     ==================================================================
      SUBROUTINE mm_excl(nn_sr,nn_esp,do_force,
     c               F_el,F_atdens_mm,F_atdens_qm,F_q_partition,
     c               q_mm,q_mm_eff,q_mm_eff_h,nn_atom,nn_sp,
     c               long_range_esp,e_lr_esp,e_lr_esp_mm,
     c               TAU_mi,mm_FION,v_mm,e_charge_rest,
     c               r_t,r_n,q_partition,q_esp,e_excl,EXTF,EXTF_l)
!     ==--------------------------------------------------------------==
!     ==  The MM atoms involved in exclusions are coupled to the QM   ==
!     ==  system via the ESP COUPLING!!!!!                            ==
!     ==--------------------------------------------------------------==
!cmb - Revised on 25 October 2013
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use timer, only: tiset, tihalt
      use system
      use parac
      use mm_dimmod
      use mm_input
      use cell
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      implicit none

C Arguments.
      integer nn_sr,nn_esp
      logical do_force
      real*8  F_el(3,*),F_atdens_mm(3,*),F_atdens_qm(3,*)
      real*8  F_q_partition(3,*),q_mm(*),q_mm_eff(*),q_mm_eff_h(*)
      integer nn_atom(*),nn_sp(*)
      real*8  e_lr_esp,e_lr_esp_mm
      logical long_range_esp
      real*8  TAU_mi(3,maxsys%nax,*),mm_FION(3,maxsys%nax,*),v_mm(*),
     &        e_charge_rest
      real*8  r_t(*),r_n(*),q_partition(*),q_esp(*),e_excl
      real*8  EXTF(fpar%kr1,fpar%kr2s,*),EXTF_l(fpar%kr1,fpar%kr2s,*)

C local dyanmical arrays.
C AK: 2005/07/22. these have to be allocatables, or else
C     the intel 8.x and compaq compiler mess up our stack.
      real*8, allocatable :: rij(:,:), vv(:,:), ff(:,:,:)
      real*8, allocatable :: P(:,:),F_mm(:,:),F_qm(:,:)
      real*8, allocatable :: q_partition_l(:,:)
      real*8, allocatable :: AA(:,:),AAm1(:,:),ff_A(:,:,:)
      real*8, allocatable :: res(:),target(:)
      real*8, allocatable :: q_mm_eff_loc(:),phi(:),at_sr(:)
      logical, allocatable :: excl_esp(:,:)
      integer, allocatable :: atom_qm_excl_l(:),atom_mm_excl_l(:)
C local variables.
      real*8 Ptmp,phitmp,AAm1tmp,q_esp_tmp,res_tmp
      real*8 xij(3),rmm(3),rij_l,qj,qj_t,qi,vv_l
      real*8 Rijsq,ft,ft_x,ft_y,ft_z
      real*8 Rijm,Rijn,Rijp,den,nom,fij,Rijl,denm2
      real*8 k_cg,deltaq,ww(2*mmdim%NATq),wq,qj_mm
      integer atom_qm_cg_l(mmdim%NATq)
      integer n_target,n_smearp,NCe_l
      integer iaa,ia,is,ia1,is1,iq,im,ix,ie,i
      integer iqqe,iqe,ime,j,k,iq1
      integer it,it_mm,it_qm,ISUB
      logical found
!
! temporary
      real*8 dx,dy,dz
      integer iy,iz

#ifdef __VECTOR 
      real*8 xij1,xij2,xij3
#endif 

      if( (excl_comm%NCe.gt.0 .or. (nn_sr+1.le.nn_esp)) 
     &     .and. cgrest_i%n_cg.gt.0)
     c      call stopgm('MM_EXCL','incompatible options',
     & __LINE__,"__file__")

      CALL TISET('   MM_EXCL',ISUB)
      ALLOCATE(excl_esp(mmdim%NATq,nn_esp))
      ALLOCATE(atom_qm_excl_l(NCe_max))
      ALLOCATE(atom_mm_excl_l(NCe_max))
      ALLOCATE(q_partition_l(maxsys%nax,mmdim%NSPq))

!$OMP parallel do private(i,j) __COLLAPSE2
#ifdef __SR11000
*poption parallel, tlocal(i,j)
#endif
      DO j=1,nn_esp
        DO i=1,mmdim%NATq
          excl_esp(i,j)=.false.
        ENDDO
      ENDDO
      wq=r_esp%esp_weight
      n_smearp=n_sm%n_smear+1
      e_excl=0.d0
      e_charge_rest=0.d0
      e_lr_esp=0.d0

      NCe_l=0
! exclusion
      DO ie=1,excl_comm%NCe
        iqe=excl_comm%atom_qm_excl(ie)
        iq=0
! test for MM atom as QM exclusion.
        if(cpsp(iqe).gt.mmdim%NSPq) call stopgm('MM_EXCL',
     &       'MM ATOM EXCLUDED AS QM ATOM',
     & __LINE__,"__file__")
        DO is=1,cpsp(iqe)-1
          iq=iq+NAq(is)
        ENDDO
        iqqe=iq+cpat(iqe)
        ime=excl_comm%atom_mm_excl(ie)
        ia1=cpat(ime)
        is1=cpsp(ime)
! test for QM atom as MM exclusion.
        if(is1.le.mmdim%NSPq) call stopgm('MM_EXCL',
     &       'QM ATOM EXCLUDED AS MM ATOM',
     & __LINE__,"__file__")
        found=.false.
        NCe_l=0
        DO iaa =1,nn_sr
          ia=nn_atom(iaa)
          is=nn_sp(iaa)
          if(ia.eq.ia1 .and. is.eq.is1)then
            found=.true.
            NCe_l=NCe_l+1
            atom_qm_excl_l(NCe_l)=iqqe
            atom_mm_excl_l(NCe_l)=iaa
          endif
        ENDDO
!$OMP parallel do private(iaa,ia,is) shared(iqqe,is1,ia1)
        DO iaa =nn_sr+1,nn_esp
          ia=nn_atom(iaa)
          is=nn_sp(iaa)
          if(ia.eq.ia1 .and. is.eq.is1)then
            found=.true.
            excl_esp(iqqe,iaa)=.true.
          endif
        ENDDO
        IF(.not.found)THEN
         call stopgm('mm_excl','MAYBE TOO SMALL ES CUTOFF, MAYBE A BUG',
     & __LINE__,"__file__")
        ENDIF
      ENDDO

!   for every couple of MM atom i in the nn list and QM atom j compute 
!   rij(i,j) 

      iq=0
      DO is1=1,mmdim%NSPq
        DO ia1=1,NAq(is1)
          iq=iq+1
          q_partition_l(ia1,is1)=q_partition(iq)
        ENDDO
      ENDDO
          
      it=0
      DO iaa =1,nn_sr
        ia=nn_atom(iaa)
        is=nn_sp(iaa)
        if(is.le.mmdim%NSPq)it=it+1
        if(DABS(q_mm(iaa)).gt.1.d-6)it=it+1
      ENDDO
      n_target=it

      ALLOCATE(AA(n_target,mmdim%NATq))
      ALLOCATE(AAm1(n_target,mmdim%NATq))
      ALLOCATE(res(n_target))
      ALLOCATE(target(n_target))
      ALLOCATE(ff_A(3,n_target,mmdim%NATq))
      ALLOCATE(at_sr(n_target))
      ALLOCATE(q_mm_eff_loc(n_target))
      ALLOCATE(phi(mmdim%NATq))
      ALLOCATE(rij(nn_sr,mmdim%NATq))
      ALLOCATE(vv(nn_sr,mmdim%Natq))
      ALLOCATE(ff(3,nn_sr,mmdim%NATq))
      ALLOCATE(P(mmdim%NATq,mmdim%NATq))
      ALLOCATE(F_mm(3,nn_sr))
      ALLOCATE(F_qm(3,mmdim%NATq))

      it=0
      DO iaa =1,nn_sr
        ia=nn_atom(iaa)
        is=nn_sp(iaa)
        if(is.le.mmdim%NSPq)then
          it=it+1
          it_qm=it
          at_sr(it)=iaa
          target(it)=wq*q_partition_l(ia,is) 
        endif
        if(DABS(q_mm(iaa)).gt.1.d-6)then
          it=it+1
          it_mm=it
          target(it)=v_mm(iaa)
          at_sr(it)=iaa
        endif
        iq=0
        do is1=1,mmdim%NSPq
          do ia1=1,NAq(is1)
            iq=iq+1
#ifdef __VECTOR 
            xij1 = TAU_mi(1,ia1,is1) - TAU_mi(1,ia,is)
            xij2 = TAU_mi(2,ia1,is1) - TAU_mi(2,ia,is)
            xij3 = TAU_mi(3,ia1,is1) - TAU_mi(3,ia,is)
            Rijsq=xij1*xij1+xij2*xij2+xij3*xij3
#else 
            Rijsq=0.d0
            DO ix=1,3
              xij(ix) = TAU_mi(ix,ia1,is1) - TAU_mi(ix,ia,is) 
              Rijsq=Rijsq+xij(ix)**2
            ENDDO
#endif 
            rij(iaa,iq) = max(1.d-8,DSQRT(Rijsq))
            Rijl=Rij(iaa,iq)
            Rijm=Rijl**(n_sm%n_smear-1)
            Rijn=Rijl*Rijm
            Rijp=Rijl*Rijn
            nom=(r_n(iaa)-Rijn)
            den=(r_t(iaa)-Rijp)
            denm2=1.0d0/(den*den)
            ft=nom*den*denm2
            fij=+(n_sm%n_smear*Rijm*den-n_smearp*Rijn*nom)*denm2
#ifdef __VECTOR
            ff(1,iaa,iq)=fij*xij1/Rij(iaa,iq)
            ff(2,iaa,iq)=fij*xij2/Rij(iaa,iq)
            ff(3,iaa,iq)=fij*xij3/Rij(iaa,iq)
#else 
            DO ix=1,3
              ff(ix,iaa,iq)=fij*xij(ix)/Rij(iaa,iq)
            ENDDO
#endif 
            vv(iaa,iq)=ft
            if(DABS(q_mm(iaa)).gt.1.d-6)then
              AA(it_mm,iq) = ft 
#ifdef __VECTOR
              ff_A(1,it_mm,iq)= ff(1,iaa,iq) 
              ff_A(2,it_mm,iq)= ff(2,iaa,iq) 
              ff_A(3,it_mm,iq)= ff(3,iaa,iq) 
#else
              DO ix=1,3
                ff_A(ix,it_mm,iq)= ff(ix,iaa,iq)
              ENDDO
#endif
            endif
            if(is.le.mmdim%NSPq)then
              if(is1.eq.is .and. ia1.eq.ia) then
                AA(it_qm,iq)=wq
              else
                AA(it_qm,iq)=0.d0
              endif
#ifdef __VECTOR
              ff_A(1,it_qm,iq)= 0.d0
              ff_A(2,it_qm,iq)= 0.d0
              ff_A(3,it_qm,iq)= 0.d0
#else
              DO ix=1,3
                ff_A(ix,it_qm,iq)= 0.d0
              ENDDO
#endif
            endif
          enddo
        enddo
      ENDDO

!
! calculate the esp charges and the pseudoinverse (in the least square
! sense) of 1/rij

      DO i=1,mmdim%NATq
        DO j=1,mmdim%NATq
cmb          P(i,j)=0.D0
          Ptmp=0.D0
          DO k=1,n_target
cmb            P(i,j)=P(i,j)+AA(k,j)*AA(k,i)
            Ptmp=Ptmp+AA(k,j)*AA(k,i)
          ENDDO
          P(i,j)=Ptmp  ! cmb - access only once P(i,j)
        ENDDO
      ENDDO

      call mm_invert(P,mmdim%NATq,mmdim%NATq,WW)  


#define __MM_USE_BLAS_LAPACK 1
cmb - Compute AAm1(i,j) - set to 0 if you have
cmb - problems with BLAS/LAPACK, but it will slow down
cmb - your code...
#if __MM_USE_BLAS_LAPACK
      CALL dgemm('N','N',n_target,mmdim%NATq,mmdim%NATq,1.0D0,
     &           AA,n_target,P,mmdim%NATq,0.0D0,AAm1,n_target)
#else
cmb - Compute AAm1(i,j) - non-BLAS/LAPACK version
      DO i=1,n_target 
        DO j=1,mmdim%Natq
          AAm1tmp=0.d0
          DO k=1,mmdim%NATq 
            AAm1tmp=AAm1tmp+P(k,j)*AA(i,k)
          ENDDO
          AAm1(i,j)=AAm1tmp
        ENDDO
      ENDDO
#endif

! calculate the efficacious charges localized on the MM atoms
cmb - set q_mm_eff_loc to zero at once
      CALL azzero(q_mm_eff_loc,n_target)
      CALL azzero(phi,mmdim%Natq)
  
!  Exclusion
      if(NCe_l.gt.0)then
        do im=1,n_target
          do ie=1,NCe_l
            iqe=atom_qm_excl_l(ie)
            ime=atom_mm_excl_l(ie)
            q_mm_eff_loc(im)=q_mm_eff_loc(im)+
     c          q_mm(ime)*vv(ime,iqe)*AAm1(im,iqe)
          enddo
        enddo

        do iq=1,mmdim%Natq
          q_esp_tmp=0.d0
          do ie=1,NCe_l
            iqe=atom_qm_Excl_l(ie)
            ime=atom_mm_excl_l(ie)
            q_esp_tmp=q_esp_tmp+q_mm(ime)*vv(ime,iqe)*P(iqe,iq)
          enddo
          phi(iq)=phi(iq)+q_esp_tmp
        enddo
      endif

! long range
      if(long_range_esp)then
        iq=0
        do is1=1,mmdim%NSPq
          do ia1=1,NAq(is1)
            iq=iq+1
            rmm(1)=TAU_mi(1,ia1,is1)
            rmm(2)=TAU_mi(2,ia1,is1)
            rmm(3)=TAU_mi(3,ia1,is1)
            vv_l=0.d0
            DO im=nn_sr+1,nn_esp
              if(excl_esp(iq,im))cycle
              ia=nn_atom(im)
              is=nn_sp(im)
cmb              qj=mm_charge(ia,is)
              xij(1) = TAU_mi(1,ia,is) - rmm(1)
              xij(2) = TAU_mi(2,ia,is) - rmm(2)
              xij(3) = TAU_mi(3,ia,is) - rmm(3)
              Rijsq=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
cmb              Rij_l=DSQRT(Rijsq)
cmb              vv_l=vv_l+qj/Rij_l
              vv_l=vv_l+mm_charge(ia,is)/DSQRT(Rijsq)
            ENDDO
!$OMP parallel do private(it)
#ifdef __SR11000
*poption parallel, tlocal(it)
#endif
            DO it=1,n_target
              q_mm_eff_loc(it)=q_mm_eff_loc(it)-vv_l*AAm1(it,iq)
            ENDDO
!$OMP parallel do private(iq1)
#ifdef __SR11000
*poption parallel, tlocal(iq1)
#endif
            DO iq1=1,mmdim%Natq
              phi(iq1)=phi(iq1)-vv_l*P(iq,iq1)
            ENDDO
          enddo
        enddo
      endif

! charge restrain
      if(cgrest_i%N_CG.gt.0)then
        do ie=1,cgrest_i%N_CG
          iqe=cgrest_i%atom_qm_cg(ie)
          iq=0
          do is=1,cpsp(iqe)-1
            iq=iq+NAq(is)
          enddo
          atom_qm_cg_l(ie)=iq+cpat(iqe)
        enddo   
!$OMP parallel do private(im,ie,iqe) __COLLAPSE2
        do im=1,n_target
          do ie=1,cgrest_i%N_CG
            iqe=atom_qm_cg_l(ie)
            q_mm_eff_loc(im)=q_mm_eff_loc(im)+AAm1(im,iqe)
          enddo
        enddo
      endif
!
      it=0
      do iaa =1,nn_sr
        q_mm_eff_h(iaa)=0.d0
        q_mm_eff(iaa)=0.d0
        is=nn_sp(iaa)
        if(is.le.mmdim%NSPq) then
          it=it+1
          q_mm_eff_h(iaa)=q_mm_eff_loc(it)
        endif
        if(DABS(q_mm(iaa)).gt.1.d-6) then
          it=it+1
          q_mm_eff(iaa)=q_mm_eff_loc(it)
        endif
      enddo

      if(.not.do_force) then
        DEALLOCATE(F_qm,F_mm,P,ff,vv,rij,phi)
        DEALLOCATE(q_mm_eff_loc,at_sr,ff_A,target,res,AAm1,AA)
        DEALLOCATE(q_partition_l,atom_mm_excl_l,atom_qm_excl_l,
     .             excl_esp)
        CALL TIHALT('   MM_EXCL',ISUB)
        return
      endif

! calculate the esp
      DO i=1,mmdim%Natq
        q_esp_tmp=0.D0
        DO j=1,n_target
          q_esp_tmp=q_esp_tmp+target(j)*AAm1(j,i)
        ENDDO
        q_esp(i)=q_esp_tmp
      ENDDO

      DO j=1,n_target
        res_tmp=0.D0
        DO ia=1,mmdim%NATq
          res_tmp= res_tmp + AA(j,ia)*q_esp(ia)
        ENDDO
        res(j)=res_tmp-target(j)
      ENDDO

      if((.not. long_range_esp).and.(NCe_l.eq.0).and.
     &     (cgrest_i%N_CG.eq.0))then
        DEALLOCATE(F_qm,F_mm,P,ff,vv,rij,phi)
        DEALLOCATE(q_mm_eff_loc,at_sr,ff_A,target,res,AAm1,AA)
        DEALLOCATE(q_partition_l,atom_mm_excl_l,atom_qm_excl_l,
     .             excl_esp)
        CALL TIHALT('   MM_EXCL',ISUB)
        return
      endif
!
      call azzero(F_QM,3*mmdim%NATq)
      call azzero(F_MM,3*nn_sr)

! long range
      if(long_range_esp .and. paral%io_parent)then
        e_lr_esp=0.d0
!$OMP parallel do private(im) reduction(+:e_lr_esp)
#ifdef __SR11000
*poption parallel, tlocal(im), psum(e_lr_esp)
#endif
        DO im=1,n_target
          e_lr_esp=e_lr_esp-target(im)*q_mm_eff_loc(im)
        ENDDO
!       write(6,*)"LONG RANGE ESP ENERGY ",e_lr_esp
        iq=0
        e_lr_esp=0.d0
        e_lr_esp_mm=0.d0
        do is1=1,mmdim%NSPq
          do ia1=1,NAq(is1)
            iq=iq+1
            rmm(1)=TAU_mi(1,ia1,is1)
            rmm(2)=TAU_mi(2,ia1,is1)
            rmm(3)=TAU_mi(3,ia1,is1)
            qj=q_esp(iq)
            qj_mm=mm_charge(ia1,is1)
            qj_t=qj-qj_mm
            do im=nn_sr+1,nn_esp
              if(excl_esp(iq,im)) cycle
              ia=nn_atom(im)
              is=nn_sp(im)
              qi=mm_charge(ia,is)
              xij(1) = TAU_mi(1,ia,is) - rmm(1)
              xij(2) = TAU_mi(2,ia,is) - rmm(2)
              xij(3) = TAU_mi(3,ia,is) - rmm(3)
              Rijsq=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
              Rij_l=DSQRT(Rijsq)
              e_lr_esp= e_lr_esp + qi*qj / Rij_l
              e_lr_esp_mm= e_lr_esp_mm + qi*qj_mm / Rij_l
              ft=qi*qj_t/(Rijsq*Rij_l)
              F_QM(1,iq)=F_QM(1,iq)-ft*xij(1)
              mm_FION(1,ia,is)= mm_FION(1,ia,is)+ft*xij(1)
              F_QM(2,iq)=F_QM(2,iq)-ft*xij(2)
              mm_FION(2,ia,is)= mm_FION(2,ia,is)+ft*xij(2)
              F_QM(3,iq)=F_QM(3,iq)-ft*xij(3)
              mm_FION(3,ia,is)= mm_FION(3,ia,is)+ft*xij(3)
            enddo
          enddo
        enddo
!       if(parent)write(6,*)"LONG RANGE ESP ENERGY ",
!    c                        e_lr_esp_mm,e_lr_esp
      endif

! exclusion
      if(NCe_l.gt.0)then
!$OMP parallel do private(ie,iqe,ime,ft_x,ft_y,ft_z)
#ifdef __SR11000
*poption parallel, tlocal(ie,iqe,ime,ft_x,ft_y,ft_z)
#endif
        DO ie=1,NCe_l
          iqe=atom_qm_Excl_l(ie)
          ime=atom_mm_excl_l(ie)
          ft_x=q_mm(ime)*q_esp(iqe)*ff(1,ime,iqe)
          ft_y=q_mm(ime)*q_esp(iqe)*ff(2,ime,iqe)
          ft_z=q_mm(ime)*q_esp(iqe)*ff(3,ime,iqe)
          F_QM(1,iqe)=F_QM(1,iqe)-ft_x
          F_QM(2,iqe)=F_QM(2,iqe)-ft_y
          F_QM(3,iqe)=F_QM(3,iqe)-ft_z
          F_MM(1,ime)=F_MM(1,ime)+ft_x
          F_MM(2,ime)=F_MM(2,ime)+ft_y
          F_MM(3,ime)=F_MM(3,ime)+ft_z
        ENDDO

        e_excl=0.d0
!$OMP parallel do private(im) reduction(+:e_excl)
#ifdef __SR11000
*poption parallel, tlocal(im), psum(e_excl)
#endif
        DO im=1,n_target
          e_excl=e_excl-target(im)*q_mm_eff_loc(im)
        ENDDO
        write(6,*)'EXCLUSION ENERGY',e_excl

        e_excl=0.d0
        DO ie=1,NCe_l
          iqe=atom_qm_Excl_l(ie)
          ime=atom_mm_excl_l(ie)
          e_excl=e_excl-q_mm(ime)*q_esp(iqe)*vv(ime,iqe)
        ENDDO
        write(6,*)'EXCLUSION ENERGY',e_excl   
      endif

      if(NCe_l.gt.0 .or. (nn_sr+1.le.nn_esp))then
        do iq=1,mmdim%NATq
          do im=1,n_target
            it=at_sr(im)
            ft=res(im)*phi(iq)+q_esp(iq)*q_mm_eff_loc(im)
            ft_x=ft*ff_A(1,im,iq)
            F_QM(1,iq)=F_QM(1,iq)+ft_x
            F_MM(1,it)=F_MM(1,it)-ft_x
            ft_y=ft*ff_A(2,im,iq)
            F_QM(2,iq)=F_QM(2,iq)+ft_y
            F_MM(2,it)=F_MM(2,it)-ft_y
            ft_z=ft*ff_A(3,im,iq)
            F_QM(3,iq)=F_QM(3,iq)+ft_z
            F_MM(3,it)=F_MM(3,it)-ft_z
          enddo
        enddo
        if(paral%io_parent) then
          iaa=0
          DO is=1,mmdim%NSPq
            DO ia=1,NAq(is)
              iaa=iaa+1
              mm_FION(1,ia,is)=mm_FION(1,ia,is)+
     c                          F_QM(1,iaa)+
     c                          F_atdens_qm(1,iaa)+
     c                          wq*F_q_partition(1,iaa)
              mm_FION(2,ia,is)=mm_FION(2,ia,is)+
     c                          F_QM(2,iaa)+
     c                          F_atdens_qm(2,iaa)+
     c                          wq*F_q_partition(2,iaa)
              mm_FION(3,ia,is)=mm_FION(3,ia,is)+
     c                          F_QM(3,iaa)+
     c                          F_atdens_qm(3,iaa)+
     c                          wq*F_q_partition(3,iaa)
            ENDDO
          ENDDO
!$OMP parallel do private(iaa,ia,is)
#ifdef __SR11000
*poption parallel, tlocal(iaa,ia,is)
#endif
          DO iaa =1,nn_sr
            ia=nn_atom(iaa)
            is=nn_sp(iaa)
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+F_MM(1,iaa)
     .                      +q_mm_eff(iaa)*(F_atdens_mm(1,iaa)
     .                      -F_el(1,iaa)) 
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+F_MM(2,iaa)
     .                      +q_mm_eff(iaa)*(F_atdens_mm(2,iaa)
     .                      -F_el(2,iaa)) 
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+F_MM(3,iaa)
     .                      +q_mm_eff(iaa)*(F_atdens_mm(3,iaa)
     .                      -F_el(3,iaa)) 
          ENDDO
        endif
        CALL DAXPY(fpar%kr1*fpar%kr2s*fpar%kr3s,
     &       1.0d0,EXTF_l(1,1,1),1,EXTF(1,1,1),1)
      endif

! charge restrain
      if(cgrest_i%N_CG.gt.0)then
        deltaq=0.d0
        do ie=1,cgrest_i%N_CG
          iqe=atom_qm_cg_l(ie)
          deltaq=deltaq+q_esp(iqe)
        enddo
        deltaq=deltaq-cgrest_r%q_rest
        k_cg=2.d0*deltaq*cgrest_r%lambda
 
        DO iq=1,mmdim%Natq
          phitmp=0.d0
          DO ie=1,cgrest_i%N_CG
            iqe=atom_qm_cg_l(ie)
            phitmp=phitmp+P(iqe,iq)
          ENDDO
          phi(iq)=phitmp*k_cg
        ENDDO
        CALL DSCAL(nn_sr,k_cg,q_mm_eff(1),1)

        CALL azzero(F_QM,3*mmdim%NATq)
        CALL azzero(F_MM,3*nn_sr)
!$OMP parallel do private(iq,im,ft,ft_x,ft_y,ft_z) __COLLAPSE2
#ifdef __SR11000
*poption parallel, tlocal(iq,im,ft,ft_x,ft_y,ft_z)
#endif
        DO iq=1,mmdim%NATq
          DO im=1,nn_sr
            ft=res(im)*phi(iq)+q_esp(iq)*q_mm_eff(im)
            ft_x=ft*ff(1,im,iq)
            ft_y=ft*ff(2,im,iq)
            ft_z=ft*ff(3,im,iq)
            F_QM(1,iq)=F_QM(1,iq)-ft_x
            F_QM(2,iq)=F_QM(2,iq)-ft_y
            F_QM(3,iq)=F_QM(3,iq)-ft_z
            F_MM(1,im)=F_MM(1,im)+ft_x
            F_MM(2,im)=F_MM(2,im)+ft_y
            F_MM(3,im)=F_MM(3,im)+ft_z
          ENDDO
        ENDDO
 
        e_charge_rest=cgrest_r%lambda*deltaq**2
 
        if(paral%io_parent) then
          write(6,*)'CHARGE RESTRAIN ENERGY',e_charge_rest,deltaq
          iaa=0
          DO is=1,mmdim%NSPq
            DO ia=1,NAq(is)
              iaa=iaa+1
              mm_FION(1,ia,is)=mm_FION(1,ia,is)+F_QM(1,iaa)
     .                        -k_cg*(F_atdens_qm(1,iaa)
     .                        +wq*F_q_partition(1,iaa))
              mm_FION(2,ia,is)=mm_FION(2,ia,is)+F_QM(2,iaa)
     .                        -k_cg*(F_atdens_qm(2,iaa)
     .                        +wq*F_q_partition(2,iaa))
              mm_FION(3,ia,is)=mm_FION(3,ia,is)+F_QM(3,iaa)
     .                        -k_cg*(F_atdens_qm(3,iaa)
     .                        +wq*F_q_partition(3,iaa))
            ENDDO
          ENDDO 
!$OMP parallel do private(iaa,ia,is)
#ifdef __SR11000
*poption parallel, tlocal(iaa,ia,is)
#endif 
          DO iaa =1,nn_sr
            ia=nn_atom(iaa)
            is=nn_sp(iaa)
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+F_MM(1,iaa)
     .                      +q_mm_eff(iaa)*(F_el(1,iaa)
     .                      -F_atdens_mm(1,iaa))
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+F_MM(2,iaa)
     .                      +q_mm_eff(iaa)*(F_el(2,iaa)
     .                      -F_atdens_mm(2,iaa))
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+F_MM(3,iaa)
     .                      +q_mm_eff(iaa)*(F_el(3,iaa)
     .                      -F_atdens_mm(3,iaa))
          ENDDO
        endif
        CALL DAXPY(fpar%kr1*fpar%kr2s*fpar%kr3s,
     &       -k_cg,EXTF_l(1,1,1),1,EXTF(1,1,1),1)
      endif 
      
      DEALLOCATE(F_qm,F_mm,P,ff,vv,rij,phi)
      DEALLOCATE(q_mm_eff_loc,at_sr,ff_A,target,res,AAm1,AA)
      DEALLOCATE(q_partition_l,atom_mm_excl_l,atom_qm_excl_l,
     .             excl_esp)
    
!     do ia=1,NATq
!       phi(ia)=0.
!     enddo
!     do ie=1,NCe_l
!       iqe=atom_qm_Excl_l(ie)
!       ime=atom_mm_excl_l(ie)
!       phi(iqe)= phi(iqe)+ q_mm(ime)*vv(ime,iqe)
!     enddo
!     DX=CELLDM(1)/DFLOAT(NR1S)
!     DY=CELLDM(2)/DFLOAT(NR2S)*CELLDM(1)
!     DZ=CELLDM(3)/DFLOAT(NR3S)*CELLDM(1) 
!     ia=0
!     do is1=1,NSPq
!       do ia1=1,NAq(is1)
!          ia=ia+1
!          ix = TAU_mi(1,ia1,is1) /dx +1
!          iy = TAU_mi(2,ia1,is1) /dy +1
!          iz = TAU_mi(3,ia1,is1) /dz +1
!          write(6,*)ia,EXTF_l(ix,iy,iz),phi(ia)
!       enddo
!     enddo
!     stop

      CALL TIHALT('   MM_EXCL',ISUB)
      RETURN
      END
