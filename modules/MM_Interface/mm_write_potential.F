C     ==================================================================
      SUBROUTINE mm_write_potential(n_stride,EXTF,cubename)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use mp_interface, only: mp_sum
C     ==--------------------------------------------------------------==
cmb - Revised on 24 March 2006

      use system
      use parac
      use cell
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen
      IMPLICIT NONE

      integer n_stride
      real*8  EXTF(fpar%kr1,fpar%kr2s,fpar%kr3s)

      real*4 , allocatable :: EXTF_pr(:,:,:) 
      real*4 , allocatable :: EXTF_pr_scr(:,:,:) 

      real*8 dx,dy,dz
      integer i,ix,iy,iz,ix2
      integer K1W,K2W,K3W
      character*30 cubename
      logical ferror

C     FIXME: check whether we can use the cubefile routine in util_p instead.
C     AK 2005/04/09

! electrostatic coupling is computed with positions that are the minimal 
! with respect to the center of mass of the quantum system.

      DX=CELL_COM%CELLDM(1)/DBLE(spar%nr1s)
      DY=CELL_COM%CELLDM(2)/DBLE(spar%nr2s)*CELL_COM%CELLDM(1)
      DZ=CELL_COM%CELLDM(3)/DBLE(spar%nr3s)*CELL_COM%CELLDM(1)

!     ----------------------------------------------------------
!     ZM, YM, XM are the Cartesian coordinates of the grid point
!     in atomic units
!     ----------------------------------------------------------
      allocate (EXTF_pr_scr(fpar%kr1s,fpar%kr2s,fpar%kr3s))
      allocate (EXTF_pr(fpar%kr1s,fpar%kr2s,fpar%kr3s))
!$OMP parallel do private(i)
      DO i=1,fpar%kr1s*fpar%kr2s*fpar%kr3s
        EXTF_pr_scr(i,1,1)=0.0d0
        EXTF_pr(i,1,1)=0.0d0
      ENDDO
!$OMP parallel do private(IZ,IY,IX,IX2)      
      DO IZ=1,spar%nr3s
        DO IY=1,spar%nr2s
          DO IX=parap%NRXPL(parai%mepos,1),parap%NRXPL(parai%mepos,2)
            IX2 = IX-parap%NRXPL(parai%mepos,1)+1   
            EXTF_pr_scr(IX,IY,IZ)=EXTF(IX2,IY,IZ)
          ENDDO
        ENDDO
      ENDDO
      call mp_sum(EXTF_pr_scr,EXTF_pr,fpar%kr1s*fpar%kr2s*fpar%kr3s,
     $     parai%cp_grp)
      K1W=(fpar%kr1s-1)/n_stride+1
      K2W=(fpar%kr2s-1)/n_stride+1
      K3W=(fpar%kr3s-1)/n_stride+1
      IF(paral%io_parent) THEN
        write(6,*)'writing cube file ',cubename
        CALL FILEOPEN(39,cubename,FO_DEF,FERROR)
        write(39,189) 'TEST'
        write(39,189) 'TEST'
        write(39,191) 1,0.0,0.0,0.0
        write(39,191) K1W, DX*n_stride,0.0,0.0
        write(39,191) K2W, 0.0,DY*n_stride,0.0
        write(39,191) K3W, 0.0,0.0,DZ*n_stride 
        write(39,192) 1, 0.0,0.0,0.0,0.0
        DO IX=1,fpar%kr1s,n_stride
          DO IY=1,fpar%kr2s,n_stride
            WRITE(39,190)(EXTF_pr(IX,IY,IZ),IZ=1,fpar%kr3s,n_stride)
          ENDDO
        ENDDO
        CALL FILECLOSE(39)
        write(6,*)'done!!!!'
      ENDIF
189   format(a4)
190   format(6E13.5)
191   format(i5,3f12.6) 
192   format(i5,4f12.6) 
      deallocate(EXTF_pr,EXTF_pr_scr)

      RETURN
      END
