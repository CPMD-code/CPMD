#if defined(__BG)
#define __VECTOR
#endif
C     ==================================================================
      SUBROUTINE mm_long_range_classic(ch_ai,dp_ai,qp_ai,xm,TAU_mi,
     c             nn_esp,nn_lr,nn_atom,nn_sp,mm_FION,
     c             ch_mm,dp_mm,qp_mm,
     c             eint_ch,eint_qp,eint_dp,e_coulomb_lr)
C     ==--------------------------------------------------------------==
cmb - Revised on 25 Januaty 2007
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use system
      use parac
      use mm_dimmod
      use ions
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      IMPLICIT NONE

      integer nn_esp,nn_lr,nn_atom(mmdim%NATm),nn_sp(mmdim%NATm)
      real*8  TAU_mi(3,maxsys%nax,maxsys%nsx)
      real*8  mm_FION(3,maxsys%nax,maxsys%nsx)
      real*8  xm(3)
      real*8  ch_ai,dp_ai(3),qp_ai(3,3)
      real*8  ch_mm,dp_mm(3),qp_mm(3,3)

c   Locals
      real*8 ft,rij2,rij,rij3,rij4,qj,qi,ft_x,ft_y,ft_z
      real*8 tr_d,tr_q,ct,e_coulomb_lr,ctq,qqjj
      integer i,in,ia,is,iaq,isq
      real*8 eint_ch,eint_qp,eint_dp
#ifdef __VECTOR
      real*8 sumfx,sumfy,sumfz
#endif 

      REAL*8 xij1,xij2,xij3,ftq1,ftq2,ftq3,ftqn1,ftqn2,ftqn3
      REAL*8 rmm1,rmm2,rmm3,v_q1,v_q2,v_q3

      ch_mm=0.D0
      dp_mm(1)=0.D0
      dp_mm(2)=0.D0
      dp_mm(3)=0.D0
      call azzero(qp_mm(1,1),9)

      ftq1=0.D0
      ftq2=0.D0
      ftq3=0.D0
      e_coulomb_lr=0.D0
      DO in=nn_esp+1,nn_lr
        ia=nn_atom(in)
        is=nn_sp(in)
        qj=mm_charge(ia,is)
        rmm1=TAU_mi(1,ia,is)
        rmm2=TAU_mi(2,ia,is)
        rmm3=TAU_mi(3,ia,is)
        xij1=rmm1-xm(1)
        xij2=rmm2-xm(2)
        xij3=rmm3-xm(3)
        rij2=xij1*xij1+xij2*xij2+xij3*xij3
        rij=DSQRT(rij2) 
        rij3=rij2*rij
        rij4=rij3*rij
        xij1=xij1/rij
        xij2=xij2/rij
        xij3=xij3/rij
        qqjj=0.5D0*qj/rij3

        ch_mm=ch_mm+qj/rij
        dp_mm(1)=dp_mm(1)+qj/rij2*xij1
        dp_mm(2)=dp_mm(2)+qj/rij2*xij2
        dp_mm(3)=dp_mm(3)+qj/rij2*xij3
        qp_mm(1,1)=qp_mm(1,1)+qqjj*xij1*xij1
        qp_mm(2,2)=qp_mm(2,2)+qqjj*xij2*xij2
        qp_mm(3,3)=qp_mm(3,3)+qqjj*xij3*xij3
        qp_mm(1,2)=qp_mm(1,2)+qqjj*xij1*xij2
        qp_mm(1,3)=qp_mm(1,3)+qqjj*xij1*xij3
        qp_mm(2,3)=qp_mm(2,3)+qqjj*xij2*xij3

        tr_d=dp_ai(1)*xij1
     .      +dp_ai(2)*xij2
     .      +dp_ai(3)*xij3
        tr_q=(qp_ai(1,1)*xij1
     .      + qp_ai(2,1)*xij2
     .      + qp_ai(3,1)*xij3)*xij1
     .      +(qp_ai(1,2)*xij1
     .      + qp_ai(2,2)*xij2
     .      + qp_ai(3,2)*xij3)*xij2
     .      +(qp_ai(1,3)*xij1
     .      + qp_ai(2,3)*xij2
     .      + qp_ai(3,3)*xij3)*xij3
        v_q1=qp_ai(1,1)*xij1
     .      +qp_ai(1,2)*xij2
     .      +qp_ai(1,3)*xij3
        v_q2=qp_ai(2,1)*xij1
     .      +qp_ai(2,2)*xij2
     .      +qp_ai(2,3)*xij3
        v_q3=qp_ai(3,1)*xij1
     .      +qp_ai(3,2)*xij2
     .      +qp_ai(3,3)*xij3

        ctq=-2.5D0*tr_q/rij4
        ct=-ch_ai/rij2-3.D0*tr_d/rij3
        ft_x=-qj*(ctq*xij1+v_q1/rij4)
        ft_y=-qj*(ctq*xij2+v_q2/rij4)
        ft_z=-qj*(ctq*xij3+v_q3/rij4)
        mm_FION(1,ia,is)=mm_FION(1,ia,is)+ft_x
     .                  -qj*(ct*xij1+dp_ai(1)/rij3)
        mm_FION(2,ia,is)=mm_FION(2,ia,is)+ft_y
     .                  -qj*(ct*xij2+dp_ai(2)/rij3) 
        mm_FION(3,ia,is)=mm_FION(3,ia,is)+ft_z
     .                  -qj*(ct*xij3+dp_ai(3)/rij3)
        ftq1=ftq1-ft_x
        ftq2=ftq2-ft_y
        ftq3=ftq3-ft_z
cmb
! subtract the contribution from the classical coulombic interaction between MM-QM couples.
#ifdef __VECTOR 
        sumfx=0.0D0
        sumfy=0.0D0
        sumfz=0.0D0
#endif 
        DO isq=1,mmdim%NSPq
          DO iaq=1,NAq(isq)
            qi=mm_charge(iaq,isq)
            xij1=rmm1-TAU_mi(1,iaq,isq)
            xij2=rmm2-TAU_mi(2,iaq,isq)
            xij3=rmm3-TAU_mi(3,iaq,isq)
            rij2=xij1*xij1+xij2*xij2+xij3*xij3
            rij=qi*qj/DSQRT(rij2) 
            e_coulomb_lr=e_coulomb_lr+rij
            ft=rij/rij2
#ifdef __VECTOR 
            sumfx=sumfx-ft*xij1
            sumfy=sumfy-ft*xij2
            sumfz=sumfz-ft*xij3
#else 
            mm_FION(1,ia,is)=mm_FION(1,ia,is)-ft*xij1
            mm_FION(2,ia,is)=mm_FION(2,ia,is)-ft*xij2
            mm_FION(3,ia,is)=mm_FION(3,ia,is)-ft*xij3
#endif 
            mm_FION(1,iaq,isq)=mm_FION(1,iaq,isq)+ft*xij1
            mm_FION(2,iaq,isq)=mm_FION(2,iaq,isq)+ft*xij2
            mm_FION(3,iaq,isq)=mm_FION(3,iaq,isq)+ft*xij3
          ENDDO
        ENDDO
#ifdef __VECTOR 
        mm_FION(1,ia,is)=mm_FION(1,ia,is)+sumfx
        mm_FION(2,ia,is)=mm_FION(2,ia,is)+sumfy
        mm_FION(3,ia,is)=mm_FION(3,ia,is)+sumfz
#endif 
      ENDDO
      qp_mm(2,1)=qp_mm(1,2)
      qp_mm(3,1)=qp_mm(1,3)
      qp_mm(3,2)=qp_mm(2,3)
      ftqn1=ftq1/DBLE(mmdim%NATq)
      ftqn2=ftq2/DBLE(mmdim%NATq)
      ftqn3=ftq3/DBLE(mmdim%NATq)
!$OMP parallel do private(is,ia) shared(ftqn1,ftqn2,ftqn3)
      DO is=1,mmdim%NSPq
        DO ia=1,NAq(is)
          mm_FION(1,ia,is)=mm_FION(1,ia,is)+ftqn1
          mm_FION(2,ia,is)=mm_FION(2,ia,is)+ftqn2
          mm_FION(3,ia,is)=mm_FION(3,ia,is)+ftqn3
        ENDDO
      ENDDO
! long range contribution on the forces of the qm ions due to the
! explicit dependence of the ionic density on the position of the ion.
      tr_q=qp_mm(1,1)+qp_mm(2,2)+qp_mm(3,3)
!$OMP parallel do private(isq,iaq,xij1,xij2,xij3,ftq1,ftq2,ftq3)
      DO isq=1,mmdim%NSPq
        DO iaq=1,NAq(isq)
          xij1=TAU_mi(1,iaq,isq)-xm(1)
          xij2=TAU_mi(2,iaq,isq)-xm(2)
          xij3=TAU_mi(3,iaq,isq)-xm(3)
          ftq1= qp_mm(1,1)*xij1
     .         +qp_mm(1,2)*xij2
     .         +qp_mm(1,3)*xij3
          ftq2= qp_mm(2,1)*xij1
     .         +qp_mm(2,2)*xij2
     .         +qp_mm(2,3)*xij3
          ftq3= qp_mm(3,1)*xij1
     .         +qp_mm(3,2)*xij2
     .         +qp_mm(3,3)*xij3
          mm_FION(1,iaq,isq)=mm_FION(1,iaq,isq)-ions0%zv(isq)*(dp_mm(1)
     +              -2.d0* xij1*tr_q+6.d0*ftq1)
          mm_FION(2,iaq,isq)=mm_FION(2,iaq,isq)-ions0%zv(isq)*(dp_mm(2)
     +              -2.d0* xij2*tr_q+6.d0*ftq2)
          mm_FION(3,iaq,isq)=mm_FION(3,iaq,isq)-ions0%zv(isq)*(dp_mm(3)
     +              -2.d0* xij3*tr_q+6.d0*ftq3)
        ENDDO
      ENDDO

! total energy of long range interaction
      eint_ch=ch_ai*ch_mm
      eint_dp = dp_ai(1)*dp_mm(1)
     +         +dp_ai(2)*dp_mm(2)
     +         +dp_ai(3)*dp_mm(3)
      eint_qp = qp_ai(1,1)*qp_mm(1,1)
     +         +qp_ai(1,2)*qp_mm(1,2)
     +         +qp_ai(1,3)*qp_mm(1,3)
     +         +qp_ai(2,1)*qp_mm(2,1)
     +         +qp_ai(2,2)*qp_mm(2,2)
     +         +qp_ai(2,3)*qp_mm(2,3)
     +         +qp_ai(3,1)*qp_mm(3,1)
     +         +qp_ai(3,2)*qp_mm(3,2)
     +         +qp_ai(3,3)*qp_mm(3,3)

      RETURN
      END
