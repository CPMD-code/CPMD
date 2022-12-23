C     ==================================================================
      SUBROUTINE mm_min_im(TAU,TAU_mi,nax,nsx,XM)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ----------------------------------------------------------------==
cmb - Revised on 23 Feb. 2005
      use mm_dimmod
      use cell
      implicit none
!
      integer nax,nsx
      real*8  TAU(3,nax,nsx)
      real*8  TAU_mi(3,nax,nsx)
      real*8  XM(3),dr(3)
!  locals
      real*8 stmp,BOXINV_au(3),Rij(3),sij(3)
      integer is,i,im
      integer ig0,igf,icg,ig
!
!  for all the particles, take the minimal image with respect to XM.
      BOXINV_au(1)= 1.D0/clsaabox%BOX_au(1)
      BOXINV_au(2)= 1.D0/clsaabox%BOX_au(2)
      BOXINV_au(3)= 1.D0/clsaabox%BOX_au(3)
! solute: loop over the charge groups
      DO icg=1,NCAG_l
        ig0=1
        if(icg.ge.2) ig0=INC_l(icg-1)+1
        igf=INC_l(icg)
        is=cpsp(ig0)
        i=cpat(ig0)
        DO im=1,3
          Rij(im)= TAU(im,i,is)-XM(im)
          stmp= Rij(im)*BOXINV_au(im)
          sij(im)= stmp -anint(stmp)
          TAU_mi(im,i,is) =sij(im)*clsaabox%BOX_au(im)+XM(im)
          dr(im)=TAU_mi(im,i,is)-TAU(im,i,is)
        ENDDO          
        DO ig=ig0+1,igf
          is=cpsp(ig)
          i=cpat(ig)
          TAU_mi(1,i,is) =TAU(1,i,is)+dr(1)
          TAU_mi(2,i,is) =TAU(2,i,is)+dr(2)
          TAU_mi(3,i,is) =TAU(3,i,is)+dr(3)
        ENDDO
      ENDDO

      RETURN
      END
