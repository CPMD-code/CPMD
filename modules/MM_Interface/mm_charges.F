#if defined(__BG)
#define __VECTOR
#endif
#if defined(__SR11000)
*option OPT(O(ss))
#endif
c     =================================================================
      SUBROUTINE mm_charges(mm_rhoe,TAU_mi,ch,q_eff,F_ch,
     .          nn_sr,nn_atom,nn_sp,EXTF)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use mp_interface, only: mp_sum
C     ----------------------------------------------------------------==
cmb - Revised on 25 May 2007
      use prmem_utils, only: prmem
      use system
      use parac
      use cell
      use mm_input
      use atwf
      use mm_dimmod
      use ions
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      IMPLICIT NONE
!
      real*8,pointer :: ARHO_save(:,:)
      common/c_atrho/ARHO_save
      real*8 dr_RHO
      common/c_atrho_r/dr_RHO
      integer NGR_RHO
      common/c_atrho_i/NGR_RHO
!
      integer nn_sr,nn_atom(nn_sr),nn_sp(nn_sr)
      real*8  mm_rhoe(fpar%kr1,fpar%kr2s,*)
      real*8  ch(mmdim%NATq)
      real*8  TAU_mi(3,maxsys%nax,*)
      real*8  F_ch(3,*)
      real*8  q_eff(*)
      real*8  EXTF(fpar%kr1,fpar%kr2s,*)
!
c   Locals
      integer ix,iy,iz,ia,is
      integer ix2,i,iaa
      real*8  dx,dy,dz
      real*8  ddx,ddy,ddz,wwt
      real*8  rx,ry,rz,rsq
      real*8  c_ch,rho_gp,fv,wwti
      real*8  ch_l(mmdim%NATq),ww(mmdim%NATq)
      real*8  TAU_q(3,mmdim%NATq),q_eff_l(mmdim%NATq)
      real*8  fionx(mmdim%NATq),fiony(mmdim%NATq),fionz(mmdim%NATq)
      real*8  F_ch_l(3,mmdim%NATq)
      real*8  f1,f2,r1,r2
      integer SP_q(mmdim%NATq)
!
      real*8  R_max,r,rr,dr_RHOm1
      real*8 rhoqq,qqt,qq,efld,dfdr
!
      real*8 energy,e_1,e_2,f_1,f_2,dddx
      integer id,iad,ixd,ia0
!
      integer NM,ISH,MMAX
      integer ir
      real*8 intf,aa,bb,fourpi
      real*8 r0,rN,alpha,c1,c2,occu
      real*8,allocatable :: ARHO_tmp(:)
      logical first
      save    first
      data    first/.true./
#ifdef __VECTOR
      INTEGER ijk,IXS,IXE,nr1ve,nr2ve,nr3ve
      REAL*8 wwia
      IXS=parap%NRXPL(parai%mepos,1) 
      IXE=parap%NRXPL(parai%mepos,2)
      nr1ve=1+IXE-IXS
      nr2ve=nr1ve*spar%NR2S
      nr3ve=nr2ve*spar%NR3S
#endif
!
      IF(n_sml%no_hirshfeld)THEN
        CALL azzero(ch,mmdim%NATq)
        CALL azzero(F_ch,3*mmdim%NATq)
        RETURN
      ENDIF
C     Sum density
      if(first)then
        NGR_RHO=100 ! these defaults
        R_max=8.0d0 ! might be changed
        dr_RHO=R_max/DBLE(NGR_RHO-1)
        fourpi=8.d0*DACOS(0.d0)
        first=.false.
        allocate(ARHO_save(NGR_RHO,mmdim%NSPq))
        allocate(ARHO_tmp(maxsys%mmaxx))
        do IS=1,mmdim%NSPq
          MMAX=atwr%MESHAT(IS)
          CALL AZZERO(ARHO_tmp,maxsys%mmaxx)
          DO ISH=1,atwf_mod%NSHELL(IS)
            IF(ISH.EQ.atwf_mod%NSHELL(IS)) THEN
              OCCU=(atwf_mod%OC(ISH,IS)-ATCHG(IS))/fourpi
            ELSE
              OCCU=atwf_mod%OC(ISH,IS)/fourpi
            ENDIF
!$OMP parallel do private(IR) shared(OCCU)
#ifdef __SR11000
*poption parallel, tlocal(IR)
#endif 
            DO IR=1,MMAX
              ARHO_tmp(IR)=ARHO_tmp(IR)+OCCU*
     .             (ATWFR(IR,ISH,IS)/ATRG(IR,IS))**2
            ENDDO
          ENDDO
!          DO IR=1,MMAX
!            ARHO_tmp(IR)=ARHO_tmp(IR)/fourpi
!          ENDDO
          NM=atwr%MESHAT(is)
          r0=ATRG(1,is)
          rN=ATRG(NM,is)
          alpha=DBLE(NM-1)/DLOG(rN/r0)
          c1=1.d0-alpha*DLOG(r0)
          c2=0.5d0*alpha
!$OMP parallel do private(i,r,rsq,ir,r1,r2,f1,f2,fv) 
!$OMP+ shared(ARHO_save)
#ifdef __SR11000
*poption parallel, tlocal(i,r,rsq,ir,r1,r2,f1,f2,fv)
#endif 
          DO i=1,NGR_RHO-2
            r=DBLE(i-1)*dr_RHO
            rsq=max(r*r,1.d-8)
            ir=MAX(1,INT(c1+c2*DLOG(rsq)))
            r1=ATRG(IR,IS)
            r2=ATRG(IR+1,IS)
            f1=ARHO_tmp(ir)
            f2=ARHO_tmp(ir+1)
            fv=f1+(r-r1)*(f2-f1)/(r2-r1)
            ARHO_save(I,IS)=fv
!           ARHO_save(I,IS)=exp(-10.d0*rsq) 
          ENDDO
          ARHO_save(NGR_RHO-1,IS)=1.d-200
          ARHO_save(NGR_RHO,IS)=1.d-200 
!  normalization
          intf=0.d0
!$OMP parallel do private(ir,f1,f2,r1,r2,bb,aa) reduction(+:intf)
#ifdef __SR11000
*poption parallel, tlocal(ir,f1,f2,r1,r2,bb,aa)
*poption psum(intf)
#endif 
          DO IR=1,NGR_RHO-1
            f1=ARHO_save(IR,IS)
            f2=ARHO_save(IR+1,IS)
            r2=DBLE(ir)*dr_RHO
            r1=r2-dr_RHO
            bb=(f2-f1)/(r2-r1)
            aa=f1-bb*r1
            intf=intf+aa/3.d0*(r2**3-r1**3)+bb/4.d0*(r2**4-r1**4)
          ENDDO
          intf=ions0%ZV(IS)/(fourpi*intf)
!$OMP parallel do private(IR) shared(intf)
#ifdef __SR11000
*poption parallel, tlocal(IR)
#endif 
          DO IR=1,NGR_RHO
            ARHO_save(IR,IS)=ARHO_save(IR,IS)*intf
          ENDDO
        enddo
        deallocate(ARHO_tmp)
        IF(paral%io_parent) CALL PRMEM('MM_CHARGES')
      endif
!
      DX=CELL_COM%CELLDM(1)/DBLE(spar%nr1s)
      DY=CELL_COM%CELLDM(1)*CELL_COM%CELLDM(2)/DBLE(spar%nr2s)
      DZ=CELL_COM%CELLDM(1)*CELL_COM%CELLDM(3)/DBLE(spar%nr3s)
      iaa=0
      do is=1,mmdim%NSPq
        do ia=1,NAq(is)
          iaa=iaa+1
          TAU_q(1,iaa)=TAU_MI(1,ia,is)
          TAU_q(2,iaa)=TAU_MI(2,ia,is)
          TAU_q(3,iaa)=TAU_MI(3,ia,is)
          SP_q(iaa)=is
          do i=1,nn_sr
            if(nn_atom(i).eq.ia .and. nn_sp(i).eq.is)goto 111 
          enddo
          write(6,'(1X,A,I4,A,I5)') 'WARNING!! SPECIES:',is,  'ATOM:',ia
          call stopgm('MM_CHARGES','QM ATOM not found NN-list',
     & __LINE__,"__file__")
111       continue
          q_eff_l(iaa)=-q_eff(i)
        enddo
      enddo
!!!! DEBUG
!          dddx=0.0000001
!          do iad =1,1
!            do ixd=1,3
!              do id=1,2
!                if(id.eq.1)TAU_q(ixd,iad)=TAU_q(ixd,iad)-dddx
!                if(id.eq.2)TAU_q(ixd,iad)=TAU_q(ixd,iad)+dddx
!!!! END DEBUG
cmb-> too many calls slow down several machines...
cmb      CALL azzero(ch,NATq) 
cmb      CALL azzero(ch_l,NATq) 
cmb      CALL azzero(F_ch,3*NATq)
cmb      CALL azzero(F_ch_l,3*NATq)
cmb-> ...do zeroing locally and in one shot
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO i=1,mmdim%NATq
        ch(i)=0.d0
        ch_l(i)=0.d0
        F_ch(1,i)=0.d0
        F_ch(2,i)=0.d0
        F_ch(3,i)=0.d0
        F_ch_l(1,i)=0.d0
        F_ch_l(2,i)=0.d0
        F_ch_l(3,i)=0.d0
      ENDDO
      c_ch=parm%omega/DBLE(spar%nr1s*spar%nr2s*spar%nr3s)
      dr_RHOm1=1.d0/dr_RHO
#ifdef __VECTOR
      DO ijk=1,nr3ve
        IX=ixs+MOD(ijk-1,nr1ve)
        IY=1+INT(MOD(ijk-1,nr2ve)/nr1ve)
        IZ=1+INT((ijk-1)/nr2ve)
        rz=DBLE(IZ-1)*DZ
        ry=DBLE(IY-1)*DY
        rx=DBLE(IX-1)*DX
        IX2 = IX-IXS+1
        wwt=0.d0
        qqt=0.d0
!$OMP parallel do private(ia,rsq,rr,ir,r1,r2,f1,f2,dfdr,fv) 
!$OMP+ private(qq,ddx,ddy,ddz,is,wwti)
!$OMP+ shared(ARHO_save,fionx,fiony,fionz,ww)
!$OMP+ reduction(+:qqt,wwt)
        DO ia=1,mmdim%NATq
          qq=q_eff_l(ia)
          ddx=rx-TAU_q(1,ia)
          ddy=ry-TAU_q(2,ia)
          ddz=rz-TAU_q(3,ia)
          is=SP_q(ia)
          rsq=ddx*ddx+ddy*ddy+ddz*ddz
          rr=MAX(DSQRT(rsq),1.d-30)
          ir=MIN(NGR_RHO-1,INT(rr*dr_RHOm1+1))
          r2=DBLE(ir)*dr_RHO
          r1=r2-dr_RHO
          f1=ARHO_save(ir,is)
          f2=ARHO_save(ir+1,is)
          dfdr=(f2-f1)*dr_RHOm1
          fv=f1+dfdr*(rr-r1)
          wwt=wwt+fv
          qqt=qqt+qq*fv
          wwti=dfdr/rr
          fionx(ia)=wwti*ddx
          fiony(ia)=wwti*ddy
          fionz(ia)=wwti*ddz
          ww(ia)=fv
        ENDDO
        wwti=1.d0/wwt
        qqt=qqt*wwti
        rho_gp=mm_rhoe(IX2,IY,IZ)*wwti*c_ch
        efld=0.d0
!$OMP parallel do private(ia,qq,wwia,rhoqq) reduction(+:efld)
!$OMP+ shared(ch_l,F_ch_l,fionx,fiony,fionz)
#ifdef __SR11000
*poption parallel, tlocal(ia,qq,wwia,rhoqq) 
*poption psum(efld)
#endif 
        DO ia=1,mmdim%NATq
           qq=q_eff_l(ia)
           wwia=ww(ia)
           ch_l(ia)=ch_l(ia)+wwia*rho_gp
           efld=efld+qq*wwia*wwti
           rhoqq=rho_gp*(qq-qqt)
           F_ch_l(1,ia)=F_ch_l(1,ia)+fionx(ia)*rhoqq
           F_ch_l(2,ia)=F_ch_l(2,ia)+fiony(ia)*rhoqq
           F_ch_l(3,ia)=F_ch_l(3,ia)+fionz(ia)*rhoqq
        ENDDO
        EXTF(IX2,IY,IZ)=EXTF(IX2,IY,IZ)-r_esp%esp_weight*efld
      ENDDO
#else
      DO IZ=1,spar%nr3s
        rz=DBLE(IZ-1)*DZ
        DO IY=1,spar%nr2s
          ry=DBLE(IY-1)*DY
          DO IX=parap%NRXPL(parai%mepos,1),parap%NRXPL(parai%mepos,2)
            rx=DBLE(IX-1)*DX
            IX2 = IX-parap%NRXPL(parai%mepos,1)+1
            wwt=0.d0
            qqt=0.d0
!$OMP parallel do private(ia,rsq,rr,ir,r1,r2,f1,f2,dfdr,fv) 
!$OMP+ private(qq,ddx,ddy,ddz,is,wwti)
!$OMP+ shared(ARHO_save,fionx,fiony,fionz,ww)
!$OMP+ reduction(+:qqt,wwt)
#ifdef __SR11000
*poption parallel 
*poption tlocal(ia,rsq,rr,ir,r1,r2,f1,f2,dfdr,fv)
*poption tlocal(qq,ddx,ddy,ddz,is,wwti)
*poption psum(qqt,wwt)
#endif
            DO ia=1,mmdim%NATq
              qq=q_eff_l(ia)
              ddx=rx-TAU_q(1,ia)
              ddy=ry-TAU_q(2,ia)
              ddz=rz-TAU_q(3,ia)
              is=SP_q(ia)
              rsq=ddx*ddx+ddy*ddy+ddz*ddz
              rr=MAX(DSQRT(rsq),1.d-30)
              ir=MIN(NGR_RHO-1,INT(rr*dr_RHOm1+1))
              r2=DBLE(ir)*dr_RHO
              r1=r2-dr_RHO
              f1=ARHO_save(ir,is)
              f2=ARHO_save(ir+1,is)
              dfdr=(f2-f1)*dr_RHOm1
              fv=f1+dfdr*(rr-r1)
              wwt=wwt+fv
              qqt=qqt+qq*fv
              wwti=dfdr/rr
              fionx(ia)=wwti*ddx
              fiony(ia)=wwti*ddy
              fionz(ia)=wwti*ddz
              ww(ia)=fv
            ENDDO
            qqt=qqt/wwt
            rho_gp=mm_rhoe(IX2,IY,IZ)*c_ch/wwt
            efld=0.d0
            DO ia=1,mmdim%NATq
              qq=q_eff_l(ia)
              ch_l(ia)=ch_l(ia)+ww(ia)*rho_gp
              efld=efld+qq*ww(ia)/wwt
              rhoqq=rho_gp*(qq-qqt)
              F_ch_l(1,ia)=F_ch_l(1,ia)+fionx(ia)*rhoqq
              F_ch_l(2,ia)=F_ch_l(2,ia)+fiony(ia)*rhoqq
              F_ch_l(3,ia)=F_ch_l(3,ia)+fionz(ia)*rhoqq
            ENDDO
            EXTF(IX2,IY,IZ)=EXTF(IX2,IY,IZ)-r_esp%esp_weight*efld
          ENDDO
        ENDDO
      ENDDO
#endif

      call mp_sum(ch_l,ch,mmdim%NATq,parai%CP_GRP)
      call mp_sum(F_ch_l,F_ch,3*mmdim%NATq,parai%CP_GRP)
      iaa=0
      energy=0.D0
      DO IS=1,mmdim%NSPq
        DO IA=1,NAq(IS)
           iaa=iaa+1
           ch(iaa)=ch(iaa)+ions0%zv(IS)
           energy=energy-q_eff_l(iaa)*ch(iaa)
        ENDDO
      ENDDO
!!!! DEBUG
!                if(id.eq.1)TAU_q(ixd,iad)=TAU_q(ixd,iad)+dddx
!                if(id.eq.2)TAU_q(ixd,iad)=TAU_q(ixd,iad)-dddx
!                if(id.eq.1)e_1=energy
!                if(id.eq.2)e_2=energy
!                if(id.eq.1)f_1=F_ch(ixd,iad)
!                if(id.eq.2)f_2=F_ch(ixd,iad)
!              enddo
!            write(6,*)iad,-(e_2-e_1)/(2.*dddx),(f_1+f_2)*0.5
!            write(16,*)iad,-(e_2-e_1)/(2.*dddx),(f_1+f_2)*0.5
!            enddo
!          enddo
!          stop   
!!!! END DEBUG   
!         do ia0=1,NATq
!           rx=TAU_q(1,ia0)
!           ry=TAU_q(2,ia0)
!           rz=TAU_q(3,ia0)
!           wwt=0.d0
!           qqt=0.d0
!           do ia=1,NATq
!             qq=q_eff_l(ia)
!             ddx=rx-TAU_q(1,ia)
!             ddy=ry-TAU_q(2,ia)
!             ddz=rz-TAU_q(3,ia)
!             is=SP_q(ia)
!             rsq=ddx*ddx+ddy*ddy+ddz*ddz
!             rr=sqrt(rsq)
!             ir=rr*dr_RHOm1+1
!             ir=min(NGR_RHO-1,ir)
!             r1=dble(ir-1)*dr_RHO
!             r2=r1+dr_RHO
!             f1=ARHO_save(ir,is)
!             f2=ARHO_save(ir+1,is)
!             dfdr=(f2-f1)*dr_RHOm1
!             fv=f1+dfdr*(rr-r1)
!             ww(ia)=fv
!             wwt=wwt+fv
!           enddo
!           wwti=1.d0/wwt
!           do ia=1,NATq
!             write(6,*)ia0,ia,ww(ia)*wwti
!           enddo
!         ENDDO
!         stop
      RETURN
      END
