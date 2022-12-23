#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE3
#else
#define __COLLAPSE3 collapse(3)
#endif
#if defined(__SR11000)
*option OPT(O(ss))
#endif
#if defined(__ALTIX)
#define __VECTOR
#endif
#if defined(__BG)
#define __VECTOR
#endif
!     ==================================================================
      SUBROUTINE mm_elstat_sr(TAU_nn,F_el,v_mm,mm_rhoe,
     c                     q_mm_eff,q_mm,nn_sp,NSPq,
     c                     e_ele_sr,nn_elstat,r_t,r_n,EXTF,EXTF_l)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
!     ----------------------------------------------------------------== 
!     This routine adds an external potential due to the MM charges
!     on the real space grid of the QM box. It also calculates the forces
!     acting on MM atoms due to the (smoothed) electron distribution
!cmb - Revised on 28 February 2014
!     ==--------------------------------------------------------------==
      use timer, only: tiset, tihalt
      use mp_interface, only: mp_sum
      use system
      use parac
      use epot_types
      use cell
      use mm_input
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      IMPLICIT NONE

C arguments
      integer nn_elstat
      integer nn_sp(*),NSPq
      real*8  TAU_nn(3,*)
      real*8  F_el(3,*)
      real*8  e_ele_sr
      real*8  q_mm_eff(*)
      real*8  q_mm(*)
      real*8  v_mm(*)
      real*8  mm_rhoe(fpar%kr1,fpar%kr2s,*)
      real*8  r_t(*)
      real*8  r_n(nn_elstat)

      real*8  EXTF(fpar%kr1,fpar%kr2s,*)
      real*8  EXTF_l(fpar%kr1,fpar%kr2s,*)

C  local variables

      real*8 , allocatable :: F_el_sc(:,:), v_mm_sc(:)

      real*8  dx,dy,dz

      real*8  xm2a,xm2b,xm3a,xm3b
#ifdef __VECTOR
      real*8  xij1,xij2a,xij2b,xij3a,xij3b
#else
      real*8  xij1,xij12,xij2a,xij2b,xij3a,xij3b
      real*8  xij2sqa,xij2sqb,xij3sqa,xij3sqb
#endif
      real*8  xij32sqaa,xij32sqab,xij32sqba,xij32sqbb

      real*8  eint_l
      integer ix,iy,iz,ix2
      integer iaa

      real*8  f1,f2,f3,v
      integer MSGLEN
      real*8  c_ch
      real*8  rhofij3a,rhofij3b
      real*8  rhofij2ab, rhofij2aa
      real*8  rhofij2bb, rhofij2ba
      real*8  rijsqaa,rijsqab
      real*8  rijsqba,rijsqbb
      real*8  rijnaa,rijnab
      real*8  rijnba,rijnbb
      real*8  ftaa,ftab
      real*8  ftba,ftbb
      real*8  rijaa,rijab
      real*8  rijba,rijbb
      real*8  nomaa,nomab,denaa,denab,denm2aa,denm2ab
      real*8  nomba,nombb,denba,denbb,denm2ba,denm2bb
      real*8  rhoe_graa,rhoe_grab
      real*8  rhoe_grba,rhoe_grbb
      real*8  rhofijaa,rhofijab
      real*8  rhofijba,rhofijbb
      INTEGER IXS,IXE,ISUB
      REAL*8  xtmp,ytmp,ztmp ! cmb
      REAL*8  dum(4) ! cmb
      real*8  temp
#ifdef __VECTOR
      INTEGER ijk,nr1ve,nr2ve,nr3ve
      REAL*8, allocatable ::  pp(:),pp1(:),pp2(:)
      integer, allocatable :: list_ix(:),list_iy(:),list_iz(:) 
#endif
      CALL TISET('MM_ELSTAT_SR',ISUB)
!
!     -------------------------------------------------
!     define the grid spacings for X,Y and Z directions
!     -------------------------------------------------
      DX=CELL_COM%CELLDM(1)/DBLE(spar%nr1s)
      DY=CELL_COM%CELLDM(1)*CELL_COM%CELLDM(2)/DBLE(spar%nr2s)
      DZ=CELL_COM%CELLDM(1)*CELL_COM%CELLDM(3)/DBLE(spar%nr3s)

!     ----------------------------------------------------------
!     ZM, YM, XM are the Cartesian coordinates of the grid point
!     in atomic units
!     ----------------------------------------------------------

      ALLOCATE(F_el_sc(3,nn_elstat),v_mm_sc(nn_elstat))
cmb-> zeroing in one shot is faster than 4 calls
cmb      CALL azzero(F_el,nn_elstat*3)
cmb      CALL azzero(F_el_sc,nn_elstat*3)
cmb      CALL azzero(v_mm,nn_elstat)
cmb      CALL azzero(v_mm_sc,nn_elstat)
!$OMP parallel do private (iaa)
#ifdef __SR11000
*poption parallel, tlocal(iaa)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO iaa=1,nn_elstat
        F_el(1,iaa)=0.0d0
        F_el(2,iaa)=0.0d0
        F_el(3,iaa)=0.0d0
        F_el_sc(1,iaa)=0.0d0
        F_el_sc(2,iaa)=0.0d0
        F_el_sc(3,iaa)=0.0d0
        v_mm(iaa)=0.0d0
        v_mm_sc(iaa)=0.0d0
      ENDDO
      e_ele_sr=0.D0
      eint_l=0.D0

! nsmear test
      if (n_sm%n_smear.ne.4) then
        call stopgm('mm_elstat_sr','nsmear should be 4',
     & __LINE__,"__file__")
      endif

! also the real space grid should be fine (the same conditions are in gcxc...F)
      if (MOD(parm%nr2,2).ne.0) then
        call stopgm('mm_elstat_sr','NR2 not programmed',
     & __LINE__,"__file__")
      endif
      if (MOD(parm%nr3,2).ne.0) then
        call stopgm('mm_elstat_sr','NR3 not programmed',
     & __LINE__,"__file__")
      endif

! the real loop over all grid points and NN atoms

      c_ch=parm%omega/DBLE(spar%nr1s*spar%nr2s*spar%nr3s)
      CALL AZZERO(EXTF,fpar%kr1*fpar%kr2s*fpar%kr3s)
      IXS=parap%NRXPL(parai%mepos,1)
      IXE=parap%NRXPL(parai%mepos,2)

#if defined(__VECTOR)
      nr1ve=1+IXE-IXS
      nr2ve=nr1ve*parm%NR2
      nr3ve=nr2ve*parm%NR3
      ALLOCATE(pp(fpar%kr1*fpar%kr2s*fpar%kr3s),
     .     pp1(fpar%kr1*fpar%kr2s*fpar%kr3s)
     .        ,pp2(fpar%kr1*fpar%kr2s*fpar%kr3s))
      CALL azzero(pp,fpar%kr1*fpar%kr2s*fpar%kr3s)
      CALL azzero(pp1,fpar%kr1*fpar%kr2s*fpar%kr3s)
      ALLOCATE(list_ix(nr3ve),list_iy(nr3ve),list_iz(nr3ve))
      
#if defined(_vpp_) || defined(__ES)
!$OMP parallel do private(ijk)
      DO ijk=1,nr3ve 
        list_ix(ijk)=ixs+MOD(ijk-1,nr1ve)
        list_iy(ijk)=1+INT(MOD(ijk-1,nr2ve)/nr1ve)
        list_iz(ijk)=1+INT((ijk-1)/nr2ve)
      ENDDO 
#else 
!$OMP parallel do private(IZ,IY,IX,ijk) __COLLAPSE3
      DO IZ=1,parm%NR3
        DO IY=1,parm%NR2
          DO IX=IXS,IXE
            ijk=1+(IX-IXS)+nr1ve*(IY-1)+nr2ve*(IZ-1)
            list_ix(ijk)=IX 
            list_iy(ijk)=IY
            list_iz(ijk)=IZ 
          ENDDO 
        ENDDO 
      ENDDO 
#endif

!$OMP parallel do private(ijk,IX,IY,IZ,IX2)
!$OMP+ shared(pp2,mm_rhoe)
      DO ijk=1,nr3ve
        IX=list_ix(ijk)
        IY=list_iy(ijk)
        IZ=list_iz(ijk)
        IX2=IX-IXS+1
        pp2(ijk)=mm_rhoe(IX2,IY,IZ)
      ENDDO
#endif

#ifdef __SR11000
*poption parallel
*poption tlocal(iaa,ijk,f1,f2,f3,v)
*poption tlocal(xtmp,ytmp,ztmp,dum)
*poption tlocal(IX,IY,IZ,xij3a,xij2a,xij1)
*poption tlocal(Rijsqaa,Rijnaa,Rijaa,nomaa,denaa,temp,denm2aa)
*poption tlocal(ftaa,RHOE_graa,rhofijaa)
*poption psum(eint_l,pp,pp1)
*poption cyclic 
#endif
      DO iaa =1,nn_elstat
        IF(DABS(q_mm(iaa)).lt.1.d-6) cycle
        f1=0D0
        f2=0D0
        f3=0D0
        v=0D0
        xtmp=tau_nn(1,iaa) ! cmb
        ytmp=tau_nn(2,iaa) ! cmb
        ztmp=tau_nn(3,iaa) ! cmb
        dum(1)=r_n(iaa)
        dum(2)=r_t(iaa)
        dum(3)=q_mm(iaa)
        dum(4)=-q_mm_eff(iaa)
#ifdef __VECTOR
!$OMP parallel do private(ijk,IX,IY,IZ)
!$OMP+ private(xm3a,xij3a,xm2a,xij2a,xij1)
!$OMP+ private(Rijsqaa,Rijnaa,Rijaa,nomaa,denaa,denm2aa,ftaa)
!$OMP+ private(RHOE_graa,temp) shared(pp,pp1,pp2,xtmp,ytmp,ztmp)
!$OMP+ reduction(+:v,f1,f2,f3)
        DO ijk=1,nr3ve
          IX=list_ix(ijk)
          IY=list_iy(ijk)
          IZ=list_iz(ijk)
          xij3a=ztmp+DBLE(1-IZ)*DZ
          xij2a=ytmp+DBLE(1-IY)*DY
          xij1 =xtmp+DBLE(1-IX)*DX
          Rijsqaa   =  xij3a*xij3a+xij2a*xij2a+xij1*xij1
          Rijnaa    =  Rijsqaa*Rijsqaa
          Rijaa     =  DSQRT(Rijsqaa)
          nomaa     =  dum(1)-Rijnaa
          denaa     =  dum(2)-Rijaa*Rijnaa
          temp      =  1.0D0/denaa
          denm2aa   =  temp*temp 
          ftaa      =  nomaa*temp 
cmb        IX2 = IX-IXS+1
          pp(ijk) =pp(ijk) -dum(3)*ftaa
          pp1(ijk)=pp1(ijk)-dum(4)*ftaa
          RHOE_graa=pp2(ijk)
          v        =v+rhoe_graa*ftaa
          rhofijaa =rhoe_graa*denm2aa*
     .        Rijsqaa*(4.0D0*denaa-5.0D0*Rijaa*nomaa)
          f1=f1 + rhofijaa*xij1
          f2=f2 + rhofijaa*xij2a
          f3=f3 + rhofijaa*xij3a
        ENDDO
#else
        DO IZ=1,parm%nr3,2
          xm3b=DBLE(IZ)*DZ
          xm3a=xm3b-DZ
cmb          xij3a=TAU_NN(3,iaa)-xm3a
cmb          xij3b=TAU_NN(3,iaa)-xm3b
          xij3a=ztmp-xm3a
          xij3b=ztmp-xm3b
          xij3sqa=xij3a*xij3a
          xij3sqb=xij3b*xij3b
          rhofij3a=0.0D0
          rhofij3b=0.0D0

          DO IY=1,parm%nr2,2
            xm2b=DBLE(IY)*DY
            xm2a=xm2b-DY
cmb            xij2a=TAU_NN(2,iaa)-xm2a
cmb            xij2b=TAU_NN(2,iaa)-xm2b
            xij2a=ytmp-xm2a
            xij2b=ytmp-xm2b
            xij2sqa=xij2a*xij2a
            xij2sqb=xij2b*xij2b
            xij32sqaa=xij3sqa+xij2sqa
            xij32sqab=xij3sqa+xij2sqb
            xij32sqba=xij3sqb+xij2sqa
            xij32sqbb=xij3sqb+xij2sqb
            rhofij2aa=0.0D0
            rhofij2ab=0.0D0
            rhofij2ba=0.0D0
            rhofij2bb=0.0D0

!!! !$OMP parallel do private(IX,xij1,xij12)
!!! !$OMP+ private(Rijsqaa,Rijsqab,Rijsqba,Rijsqbb)
!!! !$OMP+ private(Rijnaa,Rijnab,Rijnba,Rijnbb)
!!! !$OMP+ private(Rijaa,Rijab,Rijba,Rijbb)
!!! !$OMP+ private(nomaa,nomab,nomba,nombb)
!!! !$OMP+ private(denaa,denab,denba,denbb)
!!! !$OMP+ private(denm2aa,denm2ab,denm2ba,denm2bb)
!!! !$OMP+ private(ftaa,ftab,ftba,ftbb,IX2)
!!! !$OMP+ private(RHOE_graa,RHOE_grab,RHOE_grba,RHOE_grbb)
!!! !$OMP+ private(rhofijaa,rhofijab,rhofijba,rhofijbb)
!!! !$OMP+ shared(EXTF,EXTF_l,mm_rhoe)
!!! !$OMP+ reduction(+:v,f1)
!!! !$OMP+ reduction(+:rhofij2aa,rhofij2ab,rhofij2ba,rhofij2bb)
            DO IX=IXS,IXE
              xij1    = xtmp-DBLE(IX-1)*DX
              xij12   = xij1*xij1
!
              Rijsqaa = xij32sqaa+xij12
              Rijsqab = xij32sqab+xij12
              Rijsqba = xij32sqba+xij12
              Rijsqbb = xij32sqbb+xij12
              Rijnaa  = Rijsqaa*Rijsqaa
              Rijnab  = Rijsqab*Rijsqab
              Rijnba  = Rijsqba*Rijsqba
              Rijnbb  = Rijsqbb*Rijsqbb
              Rijaa   = DSQRT(Rijsqaa)
              Rijab   = DSQRT(Rijsqab)
              Rijba   = DSQRT(Rijsqba)
              Rijbb   = DSQRT(Rijsqbb)
!
              nomaa   = dum(1)-Rijnaa
              nomab   = dum(1)-Rijnab
              nomba   = dum(1)-Rijnba
              nombb   = dum(1)-Rijnbb
              denaa   = dum(2)-Rijaa*Rijnaa
              denab   = dum(2)-Rijab*Rijnab
              denba   = dum(2)-Rijba*Rijnba
              denbb   = dum(2)-Rijbb*Rijnbb
              denm2aa = 1.0D0/(denaa*denaa)
              denm2ab = 1.0D0/(denab*denab)
              denm2ba = 1.0D0/(denba*denba)
              denm2bb = 1.0D0/(denbb*denbb)
              ftaa    = nomaa*denaa*denm2aa
              ftab    = nomab*denab*denm2ab
              ftba    = nomba*denba*denm2ba
              ftbb    = nombb*denbb*denm2bb

              IX2 = IX-IXS+1

              EXTF(IX2,IY,IZ)    =EXTF(IX2,IY,IZ)-dum(3)*ftaa
              EXTF(IX2,IY+1,IZ)  =EXTF(IX2,IY+1,IZ)-dum(3)*ftab
              EXTF(IX2,IY,IZ+1)  =EXTF(IX2,IY,IZ+1)-dum(3)*ftba
              EXTF(IX2,IY+1,IZ+1)=EXTF(IX2,IY+1,IZ+1)-dum(3)*ftbb

              EXTF_l(IX2,IY,IZ)    =EXTF_l(IX2,IY,IZ)-dum(4)*ftaa
              EXTF_l(IX2,IY+1,IZ)  =EXTF_l(IX2,IY+1,IZ)-dum(4)*ftab
              EXTF_l(IX2,IY,IZ+1)  =EXTF_l(IX2,IY,IZ+1)-dum(4)*ftba
              EXTF_l(IX2,IY+1,IZ+1)=EXTF_l(IX2,IY+1,IZ+1)-dum(4)*ftbb

              RHOE_graa = mm_rhoe(IX2,IY,IZ) 
              RHOE_grab = mm_rhoe(IX2,IY+1,IZ) 
              RHOE_grba = mm_rhoe(IX2,IY,IZ+1) 
              RHOE_grbb = mm_rhoe(IX2,IY+1,IZ+1) 

              rhofijaa=RHOE_graa*denm2aa*
     +             Rijsqaa*(4.0D0*denaa - 5.0D0*Rijaa*nomaa)
              rhofijab=RHOE_grab*denm2ab*
     +             Rijsqab*(4.0D0*denab - 5.0D0*Rijab*nomab)
              rhofijba=RHOE_grba*denm2ba*
     +             Rijsqba*(4.0D0*denba - 5.0D0*Rijba*nomba)
              rhofijbb=RHOE_grbb*denm2bb*
     +             Rijsqbb*(4.0D0*denbb - 5.0D0*Rijbb*nombb)

              v=v+RHOE_graa*ftaa+RHOE_grab*ftab
     &             +RHOE_grba*ftba+RHOE_grbb*ftbb

              f1=f1+(rhofijaa+rhofijab
     .              +rhofijba+rhofijbb)*xij1

              rhofij2aa=rhofij2aa+rhofijaa
              rhofij2ab=rhofij2ab+rhofijab
              rhofij2ba=rhofij2ba+rhofijba
              rhofij2bb=rhofij2bb+rhofijbb
            ENDDO

            f2=f2+(rhofij2aa+rhofij2ba)*xij2a
     .           +(rhofij2ab+rhofij2bb)*xij2b
            rhofij3a=rhofij3a+rhofij2aa
            rhofij3a=rhofij3a+rhofij2ab
            rhofij3b=rhofij3b+rhofij2ba
            rhofij3b=rhofij3b+rhofij2bb
          ENDDO

          f3=f3+rhofij3a*xij3a+rhofij3b*xij3b   
        ENDDO
#endif
        v=v*c_ch
        F_el_sc(1,iaa)=f1*c_ch
        F_el_sc(2,iaa)=f2*c_ch
        F_el_sc(3,iaa)=f3*c_ch
        eint_l=eint_l+v*q_mm(iaa)
        v_mm_sc(iaa)  =v
      ENDDO

#ifdef __VECTOR
!$OMP parallel do private(ijk,IX,IY,IZ,IX2)
!$OMP+ shared(EXTF,EXTF_l,pp,pp1)
      DO ijk=1,nr3ve
        IX=list_ix(ijk)
        IY=list_iy(ijk)
        IZ=list_iz(ijk)
        IX2=IX-IXS+1
        EXTF(IX2,IY,IZ)  =pp(ijk)
        EXTF_l(IX2,IY,IZ)=pp1(ijk)
      ENDDO
      DEALLOCATE(list_iz,list_iy,list_ix)
      DEALLOCATE(pp2,pp1,pp)

#endif
! combine the info of the different nodes

      call mp_sum(F_el_sc,F_el,3*nn_elstat,parai%cp_grp)
      call mp_sum(v_mm_sc,v_mm,nn_elstat,parai%cp_grp)
      call mp_sum(eint_l,e_ele_sr,parai%cp_grp)
!     if(parent)write(6,*)'ELECTROSTATIC INTERACTION:',e_ele_sr

      DEALLOCATE(v_mm_sc,F_el_sc)
      CALL TIHALT('MM_ELSTAT_SR',ISUB)
      RETURN
      END
