C     ==================================================================
      SUBROUTINE  mm_elstat_atdens(TAU_nn,TAU_mi,v_mm_atdens,nn_sr,
     c         R_smear,nn_atom,nn_sp,F_atdens_mm,F_atdens_qm,q_mm_eff,
     c         r_t,r_n,atdens_fit,q_mm)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------== 
cmb - Revised on 28 Apr. 2005
      use mp_interface, only: mp_sum
      use system
      use parac
      use mm_dimmod
      use ions
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      implicit none
!
      real*8, pointer :: ARHO_save(:,:)
      common/c_atrho/ARHO_save
      real*8 dr_RHO
      common/c_atrho_r/dr_RHO
      integer NGR_RHO
      common/c_atrho_i/NGR_RHO
!
      integer nn_sr
      real*8  TAU_nn(3,nn_sr)
      real*8  F_atdens_mm(3,nn_sr)
      real*8  F_atdens_qm(3,mmdim%NATq)
      real*8  v_mm_atdens(nn_sr)
      real*8  F_atdens_mm_sc(3,nn_sr)
      real*8  F_atdens_qm_sc(3,mmdim%NATq)
      real*8  v_mm_atdens_sc(nn_sr)
      real*8  q_e
      real*8  TAU_mi(3,maxsys%nax,*)
      real*8  q_mm_eff(nn_sr)
      real*8  q_mm(nn_sr)
      real*8  R_smear(nn_sr),df
      integer nn_atom(nn_sr)
      integer nn_sp(nn_sr)
      logical atdens_fit
      real*8 r_t(nn_sr),r_n(nn_sr)
!  locals
      real*8 r3,r4,r5,r_nl,r_tl
      integer iaa,ia,is,imax,iqq
      real*8 xm(3),rx,ry,rz,v_mm_l,rc_loc ,rij


      call azzero(F_atdens_mm,3*nn_sr)
      call azzero(F_atdens_qm,3*mmdim%NATq)
      call azzero(v_mm_atdens,nn_sr)
      call azzero(F_atdens_mm_sc,3*nn_sr)
      call azzero(F_atdens_qm_sc,3*mmdim%NATq)
      call azzero(v_mm_atdens_sc,nn_sr)

      IF(atdens_fit) THEN  ! call atdens_pot
      do iaa=1,nn_sr
!       if(nn_sp(iaa).le.NSPq) cycle
        if(DABS(q_mm(iaa)).le.1.d-06) cycle
        xm(1)=TAU_nn(1,iaa)
        xm(2)=TAU_nn(2,iaa)
        xm(3)=TAU_nn(3,iaa)
        rc_loc=R_smear(iaa)
        iqq=0
        r_tl=r_t(iaa)
        r_nl=r_n(iaa) 
        q_e=q_mm_eff(iaa)
        DO is=1,mmdim%NSPq
          DO ia=1,NAq(is)
            iqq=iqq+1
            if(mod((iaa-1)*mmdim%natq+iqq,parai%CP_NPROC).ne. 
     &         parai%CP_ME) cycle
            rx=TAU_mi(1,ia,is)-xm(1)
            ry=TAU_mi(2,ia,is)-xm(2)
            rz=TAU_mi(3,ia,is)-xm(3)
            rij=MAX(1.d-8,DSQRT(rx*rx+ry*ry+rz*rz))
            CALL atdens_pot(rij,is,rc_loc,v_mm_l,df,imax)
            rx=df*rx/rij
            ry=df*ry/rij
            rz=df*rz/rij
            v_mm_atdens_sc(iaa)=v_mm_atdens_sc(iaa)+v_mm_l
            F_atdens_mm_sc(1,iaa)=F_atdens_mm_sc(1,iaa)-rx
            F_atdens_mm_sc(2,iaa)=F_atdens_mm_sc(2,iaa)-ry
            F_atdens_mm_sc(3,iaa)=F_atdens_mm_sc(3,iaa)-rz
            F_atdens_qm_sc(1,iqq)=F_atdens_qm_sc(1,iqq)+q_e*rx
            F_atdens_qm_sc(2,iqq)=F_atdens_qm_sc(2,iqq)+q_e*ry
            F_atdens_qm_sc(3,iqq)=F_atdens_qm_sc(3,iqq)+q_e*rz
          ENDDO
        ENDDO
      enddo
      ELSE  ! do the calculation without calling atdens_pot
      do iaa=1,nn_sr
!       if(nn_sp(iaa).le.NSPq) cycle
        if(DABS(q_mm(iaa)).le.1.d-06) cycle
        xm(1)=TAU_nn(1,iaa)
        xm(2)=TAU_nn(2,iaa)
        xm(3)=TAU_nn(3,iaa)
        rc_loc=R_smear(iaa)
        iqq=0
        r_tl=r_t(iaa)
        r_nl=r_n(iaa) 
        q_e=q_mm_eff(iaa)
        DO is=1,mmdim%NSPq
          DO ia=1,NAq(is)
            iqq=iqq+1
            if(mod((iaa-1)*mmdim%natq+iqq,parai%CP_NPROC).ne. 
     &         parai%cp_me) cycle
            rx=TAU_mi(1,ia,is)-xm(1)
            ry=TAU_mi(2,ia,is)-xm(2)
            rz=TAU_mi(3,ia,is)-xm(3)
            rij=MAX(1.d-8,DSQRT(rx*rx+ry*ry+rz*rz))
            r3=rij**3
            r4=r3*rij
            r5=r4*rij
            v_mm_l=(r_nl-r4)/(r_tl-r5)
            df=ions0%zv(is)*(-r3*4.d0+5.d0*r4*v_mm_l)/(r_tl-r5) 
            rx=df*rx/rij
            ry=df*ry/rij
            rz=df*rz/rij
            v_mm_atdens_sc(iaa)=v_mm_atdens_sc(iaa)+v_mm_l*ions0%zv(is)
            F_atdens_mm_sc(1,iaa)=F_atdens_mm_sc(1,iaa)-rx
            F_atdens_mm_sc(2,iaa)=F_atdens_mm_sc(2,iaa)-ry
            F_atdens_mm_sc(3,iaa)=F_atdens_mm_sc(3,iaa)-rz
            F_atdens_qm_sc(1,iqq)=F_atdens_qm_sc(1,iqq)+q_e*rx
            F_atdens_qm_sc(2,iqq)=F_atdens_qm_sc(2,iqq)+q_e*ry
            F_atdens_qm_sc(3,iqq)=F_atdens_qm_sc(3,iqq)+q_e*rz
          ENDDO
        ENDDO
      enddo
      ENDIF

      CALL mp_sum(F_atdens_qm_sc,F_atdens_qm,
     c     3*mmdim%Natq,parai%cp_grp)
      CALL mp_sum(F_atdens_mm_sc,F_atdens_mm,
     c     3*nn_sr,parai%cp_grp)
      CALL mp_sum(v_mm_atdens_sc,v_mm_atdens,
     c     nn_sr,parai%cp_grp)

      return
      end
C     ==================================================================
      SUBROUTINE atdens_pot(rij,is,rc,fv,dfv,imax)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==

      use system
      use atwf
      use mm_dimmod
      implicit none
!
      real*8, pointer :: ARHO_save(:,:)
      common/c_atrho/ARHO_save
      real*8 dr_RHO
      common/c_atrho_r/dr_RHO
      integer NGR_RHO
      common/c_atrho_i/NGR_RHO
!
      real*8 rij,fv,fvr,r,pi2,muint,rc,dfv,dmuint
      integer i,imax,is
      integer i1,NGR_RHO_2
      real*8 r_max,r1,r2,fv1,fv2,dr_loc

      pi2=4.0d0*DACOS(0.0d0)     

      fv=0.0d0
      dfv=0.0d0
      r_max=DBLE(NGR_RHO-1)*dr_RHO
      imax=NGR_RHO*4
      NGR_RHO_2=NGR_RHO-2
      dr_loc=r_max/DBLE(imax)
      DO i=1,imax
        r=(0.5d0+DBLE(i-1))*dr_loc
        i1=INT(r/dr_RHO)
        i1=MIN(NGR_RHO_2,i1-1)+1
        r1=DBLE(i1)*dr_RHO
        r2=r1+dr_RHO
        fv1=ARHO_save(i1+1,is)
        fv2=ARHO_save(i1+2,is)
        fvr=fv1+((fv2-fv1)/(r2-r1))*(r-r1)
        CALL do_muint(r,rij,rc,muint,dmuint) 
        fvr=fvr*r*r*dr_loc
        fv=fv+fvr*muint
        dfv=dfv+fvr*dmuint
      ENDDO
      fv=fv*pi2
      dfv=dfv*pi2

      RETURN
      END
C     ==================================================================
      SUBROUTINE do_muint(r,rij,rc,muint,dmuint)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
      implicit none
      real*8 r,rij,rc,muint,dmuint
      real*8 low,high
      real*8 ilow,ihigh,dilow,dihigh,dhigh,dlow
      if (r.lt.rij) then
        low =-(r-rij)/rc
        dlow=1.0d0/rc
      else
        low =+(r-rij)/rc
        dlow=-1.0d0/rc
      endif
      high=(r+rij)/rc
      dhigh=1.0d0/rc
      call do_int_v(low,ilow,dilow)  
      call do_int_v(high,ihigh,dihigh)  
      muint=(rc/(rij*r))*(ihigh-ilow)
      dmuint=-muint/rij+(rc/(rij*r))*(dihigh*dhigh - dilow*dlow)
! to normalise the density
!     muint=2
!     dmuint=0
      RETURN
      END
C     ==================================================================
      SUBROUTINE do_int_v(x,fv,dfv)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
      implicit none
      real*8 x,x2,fv,dfv
      real*8 sqrt5,sqrt55,sqrtk1,sqrtk2,atan1,atan2

      sqrt5=DSQRT(5.0d0)
      sqrt55=sqrt5/5.0d0 ! cmb
      sqrtk1=DSQRT(10.0d0+2.0d0*sqrt5)
      sqrtk2=DSQRT(10.0d0-2.0d0*sqrt5)
      atan1=DATAN((4.0d0*x+1.0d0-sqrt5)/sqrtk1)
      atan2=DATAN((4.0d0*x+1.0d0+sqrt5)/sqrtk2)
      x2=2.0d0*(x*x+1.d0)+x ! cmb

      fv=x+(sqrt55*0.5d0)*DLOG((x2-sqrt5*x)/(x2+sqrt5*x)) ! cmb
     .    +(sqrt55-1.0d0)*atan1/sqrtk1
     .    -(sqrt55+1.0d0)*atan2/sqrtk2
cmb      fv=x
cmb      fv=fv+(sqrt5/10.0d0)*DLOG(2.0d0*x*x+(1.0d0-sqrt5)*x+2.0d0)
cmb      fv=fv-(sqrt5/10.0d0)*DLOG(2.0d0*x*x+(1.0d0+sqrt5)*x+2.0d0)
cmb      fv=fv+(sqrt5/5.0d0-1)*atan1/sqrtk1
cmb      fv=fv-(sqrt5/5.0d0+1)*atan2/sqrtk2
      dfv=x*(1.0d0-x**4)/(1.0d0-x**5)

      RETURN
      END
