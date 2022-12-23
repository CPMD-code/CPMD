#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#endif
!     ==================================================================
      SUBROUTINE mm_charge_rest(nn_sr,do_force,
     .              F_el,F_atdens_mm,F_atdens_qm,F_q_partition,
     .              q_mm,q_mm_eff,nn_atom,nn_sp,TAU_mi,mm_FION,v_mm,
     .              r_t,r_n,q_partition,q_esp,e_excl,EXTF,EXTF_l)
!     ==--------------------------------------------------------------==
!cmb - Revised on 25 October 2013
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use system
      use parac
      use mm_dimmod
      use mm_input
      implicit none

      integer nn_sr
      integer nn_atom(nn_sr)
      integer nn_sp(nn_sr)
      logical do_force
      real*8  F_el(3,nn_sr)
      real*8  F_atdens_mm(3,nn_sr)
      real*8  F_atdens_qm(3,mmdim%NATq)
      real*8  F_q_partition(3,mmdim%NATq)
      real*8  q_partition(mmdim%NATq)
      real*8  q_mm(nn_sr)
      real*8  q_mm_eff(nn_sr)
      real*8  r_t(nn_sr)
      real*8  r_n(nn_sr)
      real*8  v_mm(nn_sr)
      real*8  TAU_mi(3,maxsys%nax,*)
      real*8  mm_FION(3,maxsys%nax,*)
      real*8  EXTF(fpar%kr1,fpar%kr2s,*)
      real*8  EXTF_l(fpar%kr1,fpar%kr2s,*)
!
C FIXME: this should have some allocatables here.
      real*8  F_mm(3,nn_sr)
      real*8  F_qm(3,mmdim%NATq)
      real*8 rij(nn_sr,Mmdim%Natq)
      real*8 AA(nn_sr,Mmdim%Natq)
      real*8 AAm1(nn_sr,Mmdim%Natq)
      real*8 q_esp(Mmdim%Natq)
      real*8 xij(3,nn_sr,Mmdim%Natq)
      real*8 ff(3,nn_sr,Mmdim%Natq)
      real*8 weight(nn_sr)
      real*8 P(Mmdim%Natq,mmdim%natq)
      integer iaa,ia,is,ia1,is1,iq,im,ie,i
      integer iqe,j,k
      real*8 Rijsq,ft,ft1,ft2,ft3,xtmp,ytmp,ztmp
      integer atom_qm_cg_l(Mmdim%Natq)
      integer n_smearp
      real*8 res(nn_sr),phi(Mmdim%Natq)
      real*8 k_cg,deltaq
      real*8 ww(2*Mmdim%Natq),wq
      Real*8 Rijm,Rijn,Rijp,den,nom,fij,Rijl
      real*8 q_partition_l(maxsys%nax,Mmdim%Nspq)
      real*8 e_excl
      REAL*8 res_tmp,phi_tmp,Pij,AAm1ij,q_espit

      wq=r_esp%esp_weight
      n_smearp=n_sm%n_smear+1
      e_excl=0.d0


!   for every couple of MM atom i in the nn list and QM atom j compute 
!   rij(i,j) and xij(ix,i,j)
   

      iq=0
      do is1=1,Mmdim%Nspq
        do ia1=1,NAq(is1)
          iq=iq+1
          q_partition_l(ia1,is1)=q_partition(iq)
        enddo
      enddo
          
      DO iaa =1,nn_sr
        ia=nn_atom(iaa)
        is=nn_sp(iaa)
        iq=0
        weight(iaa)=1.D0
        if(is.le.Mmdim%Nspq)then
          weight(iaa)=wq
          v_mm(iaa)=q_partition_l(ia,is)
        endif
        IF(is.gt.Mmdim%Nspq) THEN
          do is1=1,Mmdim%Nspq
            do ia1=1,NAq(is1)
              iq=iq+1
              xtmp=TAU_mi(1,ia1,is1)-TAU_mi(1,ia,is)
              ytmp=TAU_mi(2,ia1,is1)-TAU_mi(2,ia,is)
              ztmp=TAU_mi(3,ia1,is1)-TAU_mi(3,ia,is)
              Rijsq=xtmp*xtmp+ytmp*ytmp+ztmp*ztmp
              rij(iaa,iq) = max(1.d-8,DSQRT(Rijsq))
              Rijl=Rij(iaa,iq)
              Rijm=Rijl**(n_sm%n_smear-1)
              Rijn=Rijl*Rijm
              Rijp=Rijl*Rijn
              nom=(r_n(iaa)-Rijn)
              den=(r_t(iaa)-Rijp)
              ft=nom/den
              fij=(n_sm%n_smear*Rijm*den-n_smearp*Rijn*nom)/(den*den)
              AA(iaa,iq) = ft * weight(iaa)
              xij(1,iaa,iq)=xtmp/Rij(iaa,iq) ! access only once
              xij(2,iaa,iq)=ytmp/Rij(iaa,iq) ! the memory location
              xij(3,iaa,iq)=ztmp/Rij(iaa,iq) ! xij(i,iaa,iq)
              ff(1,iaa,iq)= xij(1,iaa,iq)*fij
              ff(2,iaa,iq)= xij(2,iaa,iq)*fij
              ff(3,iaa,iq)= xij(3,iaa,iq)*fij
            enddo
          enddo
        ELSE
          do is1=1,Mmdim%Nspq
            do ia1=1,NAq(is1)
              iq=iq+1
              xij(1,iaa,iq)=0.0d0
              xij(2,iaa,iq)=0.0d0
              xij(3,iaa,iq)=0.0d0
              ff(1,iaa,iq) =0.0d0
              ff(2,iaa,iq) =0.0d0
              ff(3,iaa,iq) =0.0d0
              if(is1.eq.is .and. ia1.eq.ia) then
                AA(iaa,iq)=weight(iaa)
              else
                AA(iaa,iq)=0.0d0
              endif
            enddo
          enddo
        ENDIF
      ENDDO

!
! calculate the esp charges and the pseudoinverse (in the least square
! sense) of 1/rij

      DO i=1,Mmdim%Natq
        DO j=1,Mmdim%Natq
          Pij=0.D0
          DO k=1,nn_sr
            Pij=Pij+AA(k,j)*AA(k,i)
          ENDDO
          P(i,j)=Pij ! cmb - access only once array
        ENDDO
      ENDDO
      CALL mm_invert(P,Mmdim%Natq,Mmdim%Natq,WW)  

      DO i=1,nn_sr 
        DO j=1,Mmdim%Natq
          AAm1ij=0.d0
          DO k=1,Mmdim%Natq 
            AAm1ij=AAm1ij+P(k,j)*AA(i,k)
          ENDDO
          AAm1(i,j)=AAm1ij ! cmb - access only once array
        ENDDO
      ENDDO

! calculate the efficacious charges localized on the MM atoms

      do ie=1,cgrest_i%n_cg
        iqe=cgrest_i%atom_qm_cg(ie)
        iq=0
!$OMP parallel do private(is) reduction(+:iq)
        DO is=1,cpsp(iqe)-1
          iq=iq+NAq(is)
        ENDDO
        atom_qm_cg_l(ie)=iq+cpat(iqe)
      enddo

      do im=1,nn_sr
        q_mm_eff(im)=0.d0
!$OMP parallel do private(ie,iqe)
        DO ie=1,cgrest_i%n_cg
          iqe=atom_qm_cg_l(ie)
          q_mm_eff(im)=q_mm_eff(im)+AAm1(im,iqe)
        ENDDO
      enddo

      if(.not.do_force) return

      do i=1,Mmdim%Natq
        q_espit=0.D0
        do j=1,nn_sr
          q_espit=q_espit+weight(j)*v_mm(j)*AAm1(j,i)
        enddo
        q_esp(i)=q_espit ! cmb - access only once vector
      enddo

      DO j=1,nn_sr
        res_tmp=0.D0
        DO ia=1,Mmdim%Natq
          res_tmp=res_tmp+AA(j,ia)*q_esp(ia)
        ENDDO
        res(j)=res_tmp-weight(j)*v_mm(j) ! cmb - access only once
!       write(6,*)j,res_tmp,res(j),v_mm(j)
      ENDDO

      deltaq=0.d0
!$OMP parallel do private(ie,iqe) reduction(+:deltaq)
#ifdef __SR11000
*poption parallel, tlocal(ie,iqe)
*poption psum(deltaq)
#endif
      DO ie=1,cgrest_i%n_cg
        iqe=atom_qm_cg_l(ie)
        deltaq=deltaq+q_esp(iqe)
      ENDDO
      deltaq=deltaq-cgrest_r%q_rest
      k_cg=2.d0*deltaq*cgrest_r%lambda

      DO iq=1,Mmdim%Natq
        phi_tmp=0.d0
        DO ie=1,cgrest_i%n_cg
          iqe=atom_qm_cg_l(ie)
          phi_tmp=phi_tmp+P(iqe,iq)
        ENDDO
        phi(iq)=phi_tmp*k_cg ! cmb- access only once mem loc
      ENDDO
!$OMP parallel do private(im)
#ifdef __SR11000
*poption parallel, tlocal(im)
#endif
      do im=1,nn_sr
        q_mm_eff(im)=q_mm_eff(im)*k_cg
      enddo
!
      CALL azzero(F_QM,3*Mmdim%Natq)
      CALL azzero(F_MM,3*nn_sr)
!$OMP parallel do private(ft,ft1,ft2,ft3) __COLLAPSE2
      DO iq=1,Mmdim%Natq
        DO im=1,nn_sr
          ft=res(im)*phi(iq)+q_esp(iq)*q_mm_eff(im)
          ft1=ft*ff(1,im,iq)
          ft2=ft*ff(2,im,iq)
          ft3=ft*ff(3,im,iq)
          F_QM(1,iq)=F_QM(1,iq)-ft1
          F_QM(2,iq)=F_QM(2,iq)-ft2
          F_QM(3,iq)=F_QM(3,iq)-ft3
          F_MM(1,im)=F_MM(1,im)+ft1
          F_MM(2,im)=F_MM(2,im)+ft2
          F_MM(3,im)=F_MM(3,im)+ft3
        ENDDO
      ENDDO

      e_excl=cgrest_r%lambda*deltaq**2

      if(cgrest_i%n_cg.gt.0)then
        if(paral%io_parent) then
          write(6,*)'CHARGE RESTRAIN ENERGY',e_excl,deltaq
          iaa=0
          DO is=1,Mmdim%Nspq
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
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
        do i=1,fpar%kr1*fpar%kr2s*fpar%kr3s
          EXTF(i,1,1)=EXTF(i,1,1)-k_cg*EXTF_l(i,1,1)
        enddo
      endif
      RETURN
      END
