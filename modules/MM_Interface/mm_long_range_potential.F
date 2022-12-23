C     ==================================================================
      SUBROUTINE mm_long_range_potential(ch_mm,dp_mm,qp_mm,xm,EXTF)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 25 January 2007

      use system
      use parac
      use cell
      IMPLICIT NONE
      REAL*8  ch_mm,dp_mm(*),qp_mm(3,*),xm(*)

      REAL*8  EXTF(fpar%kr1,fpar%kr2s,*)

c   Locals
      INTEGER ix,iy,iz,ix2
      REAL*8  dx,dy,dz
      REAL*8  rx,ry,rz,rsq,rx2,ry2,rz2

!     -------------------------------------------------
!     define the grid spacings for X,Y and Z directions
!     -------------------------------------------------
      DX=CELL_COM%CELLDM(1)/DBLE(spar%nr1s)
      DY=CELL_COM%CELLDM(2)/DBLE(spar%nr2s)*CELL_COM%CELLDM(1)
      DZ=CELL_COM%CELLDM(3)/DBLE(spar%nr3s)*CELL_COM%CELLDM(1)

!$OMP parallel do private(IX,IY,IZ,rx,ry,rz,rx2,ry2,rz2,rsq,IX2) 
!$OMP+shared(parai)
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
            IX2 = IX-parap%NRXPL(parai%mepos,1)+1
c              TOTAL CHARGE
            EXTF(IX2,IY,IZ)=EXTF(IX2,IY,IZ)-ch_mm
     .          -dp_mm(1)*rx-dp_mm(2)*ry-dp_mm(3)*rz  ! DIPOLE
     .          -qp_mm(1,1)*(3.D0*rx2-rsq)  !
     .          -qp_mm(2,2)*(3.D0*ry2-rsq)    !
     .          -qp_mm(3,3)*(3.D0*rz2-rsq)    ! QUADRUPOLE
     .          -qp_mm(1,2)*6.D0*rx*ry        ! 
     .          -qp_mm(1,3)*6.D0*rx*rz        !
     .          -qp_mm(3,2)*6.D0*rz*ry        !
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
