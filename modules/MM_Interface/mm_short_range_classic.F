C     ==================================================================
      SUBROUTINE mm_short_range_classic(TAU_mi,q_mm,
     c             nn_sr,nn_atom,nn_sp,mm_FION,
     c             e_coulomb_sr,e_ion_sr,v_mm_ion,r_t,r_n,
     c             gaussian_ions)
C     ==--------------------------------------------------------------==
cmb - Revised on 13 May 2005
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use system
      use mm_dimmod
      use ions
      use mm_ion_dens
      IMPLICIT NONE

      integer nn_sr,nn_atom(nn_sr),nn_sp(nn_sr)
      real*8 r_t(nn_sr),r_n(nn_sr)
      real*8  TAU_mi(3,maxsys%nax,maxsys%nsx)
      real*8  mm_FION(3,maxsys%nax,maxsys%nsx)
      real*8  q_mm(nn_sr)
      real*8  v_mm_ion(nn_sr),e_coulomb_sr,e_ion_sr
      logical gaussian_ions

c   Locals
      real*8 xij(3),rij2,rij,rij3,rij4,qj,qi,rmm(3)
      real*8 eg,deg,rg,r_tl,r_nl
      real*8 ft1,ft2,r3,r4,r5
      integer i,in,ia,is,ix,jx,iaq,isq
      integer ia1,is1,ig1,ig

      e_coulomb_sr=0.d0
      e_ion_sr=0.d0
      IF(gaussian_ions) THEN
      do in=1,nn_sr
        v_mm_ion(in)=0.d0
        ia=nn_atom(in)
        is=nn_sp(in)
        if(DABS(q_mm(in)).lt.1.d-6)cycle
        qj=q_mm(in)
        r_tl=r_t(in)
        r_nl=r_n(in)
        rmm(1)=TAU_mi(1,ia,is)
        rmm(2)=TAU_mi(2,ia,is)
        rmm(3)=TAU_mi(3,ia,is)
        rg=mm_RAGGIO(is)
        do isq=1,mmdim%NSPq
          DO iaq=1,NAq(isq)
            xij(1)=rmm(1)-TAU_mi(1,iaq,isq)
            xij(2)=rmm(2)-TAU_mi(2,iaq,isq)
            xij(3)=rmm(3)-TAU_mi(3,iaq,isq)
            rij2=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
            rij=DSQRT(rij2) 
            if(rij.lt.1.d-8)cycle
!  subtract coulomb contribution from classical force field within all the couples
! of QM and NN atoms (also the ones that are excluded in the classical force field)
            qi=mm_charge(iaq,isq)*qj/rij
            e_coulomb_sr=e_coulomb_sr+qi
            ft1=qi/rij2
!  deg is the derivative of eg (NOT the force)
            CALL gauss_sr(r_tl,r_nl,rg,rij,eg,deg)
            e_ion_sr=e_ion_sr+ qj*ions0%zv(isq)*eg
            ft2= - qj*ions0%zv(isq)*deg/rij - ft1
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+ft2*xij(1)
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+ft2*xij(2)
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+ft2*xij(3)
            mm_FION(1,iaq,isq)=mm_FION(1,iaq,isq)-ft2*xij(1)
            mm_FION(2,iaq,isq)=mm_FION(2,iaq,isq)-ft2*xij(2)
            mm_FION(3,iaq,isq)=mm_FION(3,iaq,isq)-ft2*xij(3)
            v_mm_ion(in)=v_mm_ion(in)+ions0%zv(isq)*eg
          ENDDO
        enddo
      enddo
      ELSE
      do in=1,nn_sr
        v_mm_ion(in)=0.d0
        ia=nn_atom(in)
        is=nn_sp(in)
        if(DABS(q_mm(in)).lt.1.d-6)cycle
        qj=q_mm(in)
        r_tl=r_t(in)
        r_nl=r_n(in)
        rmm(1)=TAU_mi(1,ia,is)
        rmm(2)=TAU_mi(2,ia,is)
        rmm(3)=TAU_mi(3,ia,is)
        rg=mm_RAGGIO(is)
        do isq=1,mmdim%NSPq
          DO iaq=1,NAq(isq)
            xij(1)=rmm(1)-TAU_mi(1,iaq,isq)
            xij(2)=rmm(2)-TAU_mi(2,iaq,isq)
            xij(3)=rmm(3)-TAU_mi(3,iaq,isq)
            rij2=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
            rij=DSQRT(rij2) 
            if(rij.lt.1.d-8)cycle
!  subtract coulomb contribution from classical force field within all the couples
! of QM and NN atoms (also the ones that are excluded in the classical force field)
            qi=mm_charge(iaq,isq)*qj/rij
            e_coulomb_sr=e_coulomb_sr+qi
            ft1=qi/rij2
!  deg is the derivative of eg (NOT the force)
            r3=rij2*rij
            r4=r3*rij
            r5=r4*rij
            eg=(r_nl-r4)/(r_tl-r5)
            deg=(-r3*4.d0+5.d0*r4*eg)/(r_tl-r5)   
            e_ion_sr=e_ion_sr+qj*ions0%zv(isq)*eg
            ft2= -qj*ions0%zv(isq)*deg/rij - ft1
            mm_FION(1,ia,is)=mm_FION(1,ia,is)+ft2*xij(1)
            mm_FION(2,ia,is)=mm_FION(2,ia,is)+ft2*xij(2)
            mm_FION(3,ia,is)=mm_FION(3,ia,is)+ft2*xij(3)
            mm_FION(1,iaq,isq)=mm_FION(1,iaq,isq)-ft2*xij(1)
            mm_FION(2,iaq,isq)=mm_FION(2,iaq,isq)-ft2*xij(2)
            mm_FION(3,iaq,isq)=mm_FION(3,iaq,isq)-ft2*xij(3)
            v_mm_ion(in)=v_mm_ion(in)+ions0%zv(isq)*eg
          ENDDO
        enddo
      enddo
      ENDIF

! uncomment to include the mechanical coupling interaction 
! between the QM and MM excluded couples
! 
!     do in=1,nn_sr
!       qj=q_mm(in)
!       if(abs(qj).lt.1.e-6)cycle
!       ia=nn_atom(in)
!       is=nn_sp(in)
!       rmm(1)=TAU_mi(1,ia,is)
!       rmm(2)=TAU_mi(2,ia,is)
!       rmm(3)=TAU_mi(3,ia,is)
!       ig=gratom((is-1)*maxsys%nax+ia)
!       do i=KNE(ig)+1,KNE(ig)+INE(ig)
!         ig1=JSNE(i)
!         ia1=cpat(ig1)
!         is1=cpsp(ig1)
!         if((is1.le.NSPq .and. is.gt.NSPq)  .or. 
!    c       (is1.gt.NSPq .and. is.le.NSPq))then
!           xij(1)=rmm(1)-TAU_mi(1,ia1,is1)
!           xij(2)=rmm(2)-TAU_mi(2,ia1,is1)
!           xij(3)=rmm(3)-TAU_mi(3,ia1,is1)
!           rij2=xij(1)*xij(1)+xij(2)*xij(2)+xij(3)*xij(3)
!           rij=sqrt(rij2) 
!           if(abs(rij).lt.1.d-8)cycle
!           qi=mm_charge(ia1,is1)
!           e_coulomb_sr=e_coulomb_sr-qi*qj/rij
!           ft1=qi*qj/(rij*rij2)
!           mm_FION(1,ia,is)=mm_FION(1,ia,is)+ft1*xij(1)
!           mm_FION(2,ia,is)=mm_FION(2,ia,is)+ft1*xij(2)
!           mm_FION(3,ia,is)=mm_FION(3,ia,is)+ft1*xij(3)
!           mm_FION(1,ia1,is1)=mm_FION(1,ia1,is1)-ft1*xij(1)
!           mm_FION(2,ia1,is1)=mm_FION(2,ia1,is1)-ft1*xij(2)
!           mm_FION(3,ia1,is1)=mm_FION(3,ia1,is1)-ft1*xij(3)
!         endif
!       enddo
!     enddo

      RETURN
      END
C     ==================================================================
       SUBROUTINE gauss_sr(r_t,r_n,rg,rij,fint,dfint)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
       implicit none
       real*8 r_t,r_n,rg,rij,fint,dfint
!
       real*8 cn,r0,rf,r,dr
       real*8 rgm1,r4,fp,rijl,r3
       integer ir

       integer ngr_r
       parameter(ngr_r=1000)
       real*8 pi
       parameter(pi=3.141592653589793d0)

! INTERACTION WITH A POINT CHARGE (FOR THE IONIC PART)
!      r3=rij**3
!      r4=r3*rij
!      r5=r4*rij
!      fint=(r_n-r4)/(r_t-r5)
!      dfint=(-r3*4.d0+5.d0*r4*fint)/(r_t-r5)
!      return

! INTERACTION WITH A GAUSSIAN (FOR THE IONIC PART)
       rijl=MAX(1.d-8,rij)
cmb       cn=1.0d0/pi**1.5d0/rg**3
       cn=1.0d0/DSQRT(pi)/rg**3
       
       rgm1=1.d0/rg**2
       r0=-MAX(0.d0,-rijl-5.d0*rg)
       rf=-rijl+5.0d0*rg
       dr=(rf-r0)/DBLE(ngr_r)

       r=r0
       r4=r**4
       fint=
     .  0.5d0*r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r+rijl)**2*rgm1)*dr
       dfint=0.d0
       r=r+dr
       DO ir=2,ngr_r-1
         r4=r**4
         fp=r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r+rijl)**2*rgm1)*dr
         fint=fint+fp
         dfint=dfint-fp*(r+rijl)
         r=r+dr
       ENDDO
       r4=r**4
       fint=fint+
     .  0.5d0*r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r+rijl)**2*rgm1)*dr

       r0=-MAX(0.d0,rijl-5.d0*rg)
       rf=rijl+5.0d0*rg
       dr=(rf-r0)/DBLE(ngr_r)

       r=r0
       r4=r**4
       fint=fint-
     .  0.5d0*r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r-rijl)**2*rgm1)*dr
       r=r+dr
       DO ir=2,ngr_r-1
         r4=r**4
         fp=r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r-rijl)**2*rgm1)*dr
         fint=fint-fp
         dfint=dfint+fp*(rijl-r)
         r=r+dr
       ENDDO
       r4=r**4
       fint=fint-
     .  0.5d0*r*(r_n-r4)/(r_t-r4*r)*DEXP(-(r-rijl)**2*rgm1)*dr
 
cmb       fint=-rg**2/rijl*fint*pi*cn
cmb       dfint=-rg**2/rijl*dfint*pi*cn*2.d0*rgm1-fint/rijl
       fint=-rg**2/rijl*fint*cn
       dfint=-rg**2/rijl*dfint*cn*2.d0*rgm1-fint/rijl

       IF(rijl.eq.1.d-8) dfint=0.d0

       RETURN
       END
