#if defined(__BG)
#define __VECTOR
#endif
C     ==================================================================
      SUBROUTINE mm_multipole(ch_ai,dp_ai,qp_ai,xm,TAU_mi,mm_rhoe)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 12 June 2006

      use mp_interface, only: mp_sum
      use system
      use parac
      use cell
      use mm_dimmod
      use ions
      IMPLICIT NONE
      REAL*8  mm_rhoe(fpar%kr1,fpar%kr2s,fpar%kr3s)
      REAL*8  ch_ai,dp_ai(3),qp_ai(3,3)
      REAL*8  TAU_mi(3,maxsys%nax,maxsys%nsx)

c   Locals
      REAL*8  ch_ais,dp_ais(3),qp_ais(3,3)
      REAL*8  xm(3)
      REAL*8  rhoe_gr
      INTEGER ix,iy,iz
      INTEGER ix2,i,isq,iaq
      INTEGER MSGLEN
      REAL*8  dx,dy,dz
      REAL*8  rx,ry,rz,rsq
      REAL*8  c_ch
#ifdef __VECTOR
      INTEGER ixs,ixe,nr1ve,nr2ve,nr3ve,ijk
      REAL*8 dp_ai1,dp_ai2,dp_ai3,qp_a11,qp_a22,qp_a33,
     .       qp_a12,qp_a23,qp_a13
#else
      REAL*8 rx2,ry2,rz2
#endif

!     -------------------------------------------------
!     define the grid spacings for X,Y and Z directions
!     -------------------------------------------------
      DX=CELL_COM%CELLDM(1)/DBLE(spar%nr1s)
      DY=CELL_COM%CELLDM(2)/DBLE(spar%nr2s)*CELL_COM%CELLDM(1)
      DZ=CELL_COM%CELLDM(3)/DBLE(spar%nr3s)*CELL_COM%CELLDM(1)
      
      ch_ais=0.D0
      ch_ai=0.D0
      DO i=1,3
         dp_ais(i)=0.D0
         dp_ai(i)=0.D0
         qp_ais(i,1)=0.D0
         qp_ais(i,2)=0.D0
         qp_ais(i,3)=0.D0
         qp_ai(i,1)=0.D0
         qp_ai(i,2)=0.D0
         qp_ai(i,3)=0.D0
      ENDDO

      c_ch=parm%omega/DBLE(spar%nr1s*spar%nr2s*spar%nr3s)
#ifdef __VECTOR
      ixs=parap%NRXPL(parai%mepos,1)
      ixe=parap%NRXPL(parai%mepos,2)
      nr1ve=1+ixe-ixs
      nr2ve=nr1ve*spar%NR2S
      nr3ve=nr2ve*spar%NR3S
      dp_ai1=0.d0
      dp_ai2=0.d0
      dp_ai3=0.d0
      qp_a11=0.d0
      qp_a22=0.d0
      qp_a33=0.d0
      qp_a12=0.d0
      qp_a23=0.d0
      qp_a13=0.d0
!$OMP parallel do private(ijk,ix,iy,iz,rx,ry,rz,rsq,ix2)
!$OMP+ private(RHOE_gr) 
!$OMP+ reduction(+:ch_ais,dp_ai1,dp_ai2,dp_ai3)
!$OMP+ reduction(+:qp_a11,qp_a22,qp_a33,qp_a12,qp_a23,qp_a13)
#ifdef __SR11000
*poption parallel
*poption tlocal(ijk,ix,iy,iz,rx,ry,rz,rsq,ix2,RHOE_gr)
*poption psum(ch_ais,dp_ai1,dp_ai2,dp_ai3)
*poption psum(qp_a11,qp_a22,qp_a33,qp_a12,qp_a23,qp_a13)
#endif 
      DO ijk=1,nr3ve
        ix=ixs+MOD(ijk-1,nr1ve)
        iy=1+INT(MOD(ijk-1,nr2ve)/nr1ve)
        iz=1+INT((ijk-1)/nr2ve)
        rz=DBLE(IZ-1)*DZ-xm(3)
        ry=DBLE(IY-1)*DY-xm(2)
        rx=DBLE(IX-1)*DX-xm(1)
        rsq=rx*rx+ry*ry+rz*rz
        ix2=ix-ixs+1
        RHOE_gr=mm_rhoe(ix2,iy,iz)*c_ch
c       TOTAL CHARGE
        ch_ais=ch_ais+RHOE_gr
c       DIPOLE
        dp_ai1=dp_ai1+RHOE_gr*rx
        dp_ai2=dp_ai2+RHOE_gr*ry
        dp_ai3=dp_ai3+RHOE_gr*rz
c       QUADRUPOLE
        qp_a11=qp_a11+RHOE_gr*(3.D0*rx*rx-rsq)
        qp_a22=qp_a22+RHOE_gr*(3.D0*ry*ry-rsq)
        qp_a33=qp_a33+RHOE_gr*(3.D0*rz*rz-rsq)
        qp_a12=qp_a12+RHOE_gr*3.D0*rx*ry
        qp_a23=qp_a23+RHOE_gr*3.D0*ry*rz
        qp_a13=qp_a13+RHOE_gr*3.D0*rx*rz
      ENDDO
      dp_ais(1)=dp_ai1
      dp_ais(2)=dp_ai2
      dp_ais(3)=dp_ai3
      qp_ais(1,1)=qp_a11 
      qp_ais(2,2)=qp_a22 
      qp_ais(3,3)=qp_a33 
      qp_ais(1,2)=qp_a12 
      qp_ais(2,3)=qp_a23 
      qp_ais(1,3)=qp_a13 
#else
      DO IZ=1,spar%nr3s
        rz=DBLE(IZ-1)*DZ-xm(3)
        rz2=rz*rz
        DO IY=1,spar%nr2s
          ry=DBLE(IY-1)*DY-xm(2)
          ry2=ry*ry
          DO IX=parap%NRXPL(parai%mepos,1),parap%NRXPL(parai%mepos,2)
            rx=DBLE(IX-1)*DX-xm(1)
            rx2=rx*rx
            rsq=rx2+ry2+rz2
            IX2=IX-parap%NRXPL(parai%mepos,1)+1
            RHOE_gr=mm_rhoe(IX2,IY,IZ)*c_ch
c           TOTAL CHARGE
            ch_ais=ch_ais+RHOE_gr 
c           DIPOLE
            dp_ais(1)=dp_ais(1)+RHOE_gr*rx
            dp_ais(2)=dp_ais(2)+RHOE_gr*ry
            dp_ais(3)=dp_ais(3)+RHOE_gr*rz
c           QUADRUPOLE
            qp_ais(1,1)=qp_ais(1,1)+RHOE_gr*(3.D0*rx2-rsq)
            qp_ais(2,2)=qp_ais(2,2)+RHOE_gr*(3.D0*ry2-rsq)
            qp_ais(3,3)=qp_ais(3,3)+RHOE_gr*(3.D0*rz2-rsq)
            qp_ais(1,2)=qp_ais(1,2)+RHOE_gr*3.D0*rx*ry
            qp_ais(2,3)=qp_ais(2,3)+RHOE_gr*3.D0*ry*rz
            qp_ais(1,3)=qp_ais(1,3)+RHOE_gr*3.D0*rx*rz
          ENDDO
        ENDDO
      ENDDO
#endif 
      qp_ais(2,1)=qp_ais(1,2)
      qp_ais(3,2)=qp_ais(2,3)
      qp_ais(3,1)=qp_ais(1,3)
      CALL mp_sum(ch_ais,ch_ai,parai%cp_grp)
      CALL mp_sum(dp_ais,dp_ai,3,parai%cp_grp)
      CALL mp_sum(qp_ais,qp_ai,9,parai%cp_grp)

      DO isq=1,mmdim%NSPq
        RHOE_gr=ions0%ZV(isq)
        DO iaq=1,NAq(isq)
          rx=TAU_mi(1,iaq,isq)-xm(1)
          ry=TAU_mi(2,iaq,isq)-xm(2)
          rz=TAU_mi(3,iaq,isq)-xm(3)
          rsq=rx*rx+ry*ry+rz*rz
c         TOTAL CHARGE
          ch_ai=ch_ai + RHOE_gr
c         DIPOLE
          dp_ai(1)=dp_ai(1)+RHOE_gr*rx
          dp_ai(2)=dp_ai(2)+RHOE_gr*ry
          dp_ai(3)=dp_ai(3)+RHOE_gr*rz
c         QUADRUPOLE
          qp_ai(1,1)=qp_ai(1,1)+RHOE_gr*(3.D0*rx*rx-rsq)
          qp_ai(2,2)=qp_ai(2,2)+RHOE_gr*(3.D0*ry*ry-rsq)
          qp_ai(3,3)=qp_ai(3,3)+RHOE_gr*(3.D0*rz*rz-rsq)
          qp_ai(1,2)=qp_ai(1,2)+RHOE_gr*3.D0*rx*ry
          qp_ai(2,3)=qp_ai(2,3)+RHOE_gr*3.D0*ry*rz
          qp_ai(1,3)=qp_ai(1,3)+RHOE_gr*3.D0*rx*rz
        ENDDO
      ENDDO
      qp_ai(2,1)=qp_ai(1,2)
      qp_ai(3,2)=qp_ai(2,3)
      qp_ai(3,1)=qp_ai(1,3)

      RETURN
      END
