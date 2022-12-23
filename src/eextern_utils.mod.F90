MODULE eextern_utils
  USE atomc_utils,                     ONLY: positx
  USE cnst,                            ONLY: pi
  USE efld,                            ONLY: extf
  USE epot_types,                      ONLY: epot1,&
                                             epot2,&
                                             epot3,&
                                             epot4
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE ragg,                            ONLY: raggio
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eextern
  !public :: vextnu
  !public :: vextnub

CONTAINS

  ! ==================================================================
  SUBROUTINE eextern(rhoe,v,eirop,tau0,fion,eext,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  ENERGY OF INTERACTION OF AN EXTERNAL FIELD (EXTF) WITH      ==
    ! ==  THE QM CHARGE DISTIBUTION (RHOE) AND THE FORCES ON THE      ==
    ! ==  IONS DUE TO THIS INTERACTION                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: v(*), eirop(*)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), eext
    LOGICAL                                  :: tfor

    INTEGER                                  :: ia, is, isub
    REAL(real_8)                             :: eext2, eexti, ffi(3), gfac, &
                                                xi, yi, zi
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('      EEXT',isub)
    ! ..ENERGY VEXT*RHOE
    eext=0._real_8
    IF (tfor.OR.cntl%texadd) THEN
       gfac=parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
       eext=ddot(fpar%nnr1,rhoe,1,extf,1)*gfac
       ! ..ENERGY VEXT*RHOI AND GRADIENT
       eext2=0._real_8
       eexti=0._real_8
       IF (.NOT.lqmmm%qmmm)THEN
          ! EGO-CPMD interface
          IF (cnti%iftype.EQ.1) THEN
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   CALL vextnu(extf,tau0(:,ia,is),ions0%zv(is),&
                        raggio(is),eext2,ffi,tfor)
                   fion(1,ia,is)=fion(1,ia,is)+ffi(1)*gfac
                   fion(2,ia,is)=fion(2,ia,is)+ffi(2)*gfac
                   fion(3,ia,is)=fion(3,ia,is)+ffi(3)*gfac
                ENDDO
             ENDDO
             eext2=eext2*gfac
             ! Gromacs-CPMD interface. biswas/2005
          ELSE IF (cnti%iftype.EQ.2.AND.paral%parent) THEN
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   xi=tau0(1,ia,is)
                   yi=tau0(2,ia,is)
                   zi=tau0(3,ia,is)
                   CALL vextnub(is,ia,xi,yi,zi,ions0%zv(is),eexti,ffi)
                   fion(1,ia,is)=fion(1,ia,is)+ffi(1)
                   fion(2,ia,is)=fion(2,ia,is)+ffi(2)
                   fion(3,ia,is)=fion(3,ia,is)+ffi(3)
                   eext2 = eext2 + eexti
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       eext=eext+eext2
    ENDIF
    CALL tihalt('      EEXT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eextern
  ! ==================================================================
  ! external potential for the EGO-CPMD interface
  ! ==================================================================
  SUBROUTINE vextnu(extf,tau,zval,raggio,eext,forc,tfor)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: extf(fpar%kr1,fpar%kr2,*), &
                                                tau(*), zval, raggio, eext, &
                                                forc(3)
    LOGICAL                                  :: tfor

    REAL(real_8), PARAMETER                  :: reps = 1.e-8_real_8 

    INTEGER                                  :: i, i1, ii, iorgx, iorgy, &
                                                iorgz, ix, j, j1, jj, k, k1, &
                                                kk, mr(3), nptau(3)
    REAL(real_8)                             :: da1(3,3), r2, rcut, rcut2, &
                                                sc1(3), sc2(3), x, x1, &
                                                xptau(3), xval, xx1, xx2, y, &
                                                y1, yval, z, z1, zr1, zr2

! ==--------------------------------------------------------------==

    zr1=-zval/raggio**3*pi**(-1.5_real_8)
    zr2=1._real_8/raggio**2
    forc(1)=0._real_8
    forc(2)=0._real_8
    forc(3)=0._real_8
    ! ..Atom position
    x1=tau(1)
    y1=tau(2)
    z1=tau(3)
    CALL pbc(x1,y1,z1,sc1(1),sc1(2),sc1(3),1,parm%apbc,parm%ibrav)
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    IF (sc2(1).LT.0._real_8) sc2(1)=sc2(1)+1.0_real_8
    IF (sc2(2).LT.0._real_8) sc2(2)=sc2(2)+1.0_real_8
    IF (sc2(3).LT.0._real_8) sc2(3)=sc2(3)+1.0_real_8
    nptau(1)=NINT(sc2(1)*spar%nr1s)+1
    nptau(2)=NINT(sc2(2)*spar%nr2s)+1
    nptau(3)=NINT(sc2(3)*spar%nr3s)+1
    sc1(1)=sc2(1)-REAL(nptau(1)-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
    sc1(2)=sc2(2)-REAL(nptau(2)-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
    sc1(3)=sc2(3)-REAL(nptau(3)-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)
    IF (nptau(1).GT.spar%nr1s) nptau(1)=1
    IF (nptau(2).GT.spar%nr2s) nptau(2)=1
    IF (nptau(3).GT.spar%nr3s) nptau(3)=1
    CALL dgemv('T',3,3,1._real_8,metr_com%ht,3,sc1,1,0._real_8,xptau,1)
    DO i=1,3
       da1(i,1)=parm%a1(i)/REAL(spar%nr1s,kind=real_8)
       da1(i,2)=parm%a2(i)/REAL(spar%nr2s,kind=real_8)
       da1(i,3)=parm%a3(i)/REAL(spar%nr3s,kind=real_8)
    ENDDO
    ! ..Interaction region
    xx1=reps*raggio**3/zval*pi**1.5_real_8
    xx2=LOG(xx1)
    rcut=-raggio**2*xx2
    rcut2=rcut*rcut
    sc1(1)=rcut
    sc1(2)=rcut
    sc1(3)=rcut
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    mr(1)=NINT(sc2(1)*spar%nr1s)+1
    mr(2)=NINT(sc2(2)*spar%nr2s)+1
    mr(3)=NINT(sc2(3)*spar%nr3s)+1
    mr(1)=MIN(mr(1),spar%nr1s/2-1)
    mr(2)=MIN(mr(2),spar%nr2s/2-1)
    mr(3)=MIN(mr(3),spar%nr3s/2-1)
    ! ..Origin of density grid
    iorgx=nptau(1)-mr(1)
    iorgy=nptau(2)-mr(2)
    iorgz=nptau(3)-mr(3)
    IF (iorgx.LT.1) iorgx=iorgx+spar%nr1s
    IF (iorgy.LT.1) iorgy=iorgy+spar%nr2s
    IF (iorgz.LT.1) iorgz=iorgz+spar%nr3s
    ! ..loop over all grid points
    DO i=0,2*mr(1)
       ix=iorgx+i
       IF (ix.GT.spar%nr1s) ix=ix-spar%nr1s
       IF (ix.GE.parap%nrxpl(parai%mepos,1).AND.ix.LE.parap%nrxpl(parai%mepos,2)) THEN
          ii=ix-parap%nrxpl(parai%mepos,1)+1
          i1=mr(1)-i
          DO j=0,2*mr(2)
             jj=iorgy+j
             IF (jj.GT.spar%nr2s) jj=jj-spar%nr2s
             j1=mr(2)-j
             DO k=0,2*mr(3)
                kk=iorgz+k
                IF (kk.GT.spar%nr3s) kk=kk-spar%nr3s
                k1=mr(3)-k
                CALL positx(x,y,z,i1,j1,k1,da1,xptau(1))
                r2=x*x+y*y+z*z
                IF (r2.LT.rcut2) THEN
                   xval=zr1*EXP(-r2*zr2)*extf(ii,jj,kk)
                   eext=eext+xval
                   IF (tfor) THEN
                      yval=-2._real_8*xval*zr2
                      forc(1)=forc(1)-x*yval
                      forc(2)=forc(2)-y*yval
                      forc(3)=forc(3)-z*yval
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vextnu
  ! ==================================================================
  ! ..this subroutine is added for gromacs-cpmd interface: biswas
  ! ..it adopts a localized expansion for the MM atom point charges and a
  ! ..multilayer (inner, intermediate, outer) scheme with periodic updating 
  ! ..of qmmm energy & forecs for the intermediate and outer layer MM-atoms
  ! ..while evaluating the energy and forces for the inner layer at each step.
  ! ==================================================================
  SUBROUTINE vextnub(is,ia,xi,yi,zi,zvp,eext,forc)
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: is, ia
    REAL(real_8)                             :: xi, yi, zi, zvp, eext, forc(3)

    INTEGER                                  :: idnr, inr
    REAL(real_8) :: alam, alfa, beta, efintm(4,100,100), efoutr(4,100,100), &
      eintm, eoutr, fintm_1, fintm_2, fintm_3, foutr_1, foutr_2, foutr_3, &
      rij, rijq, rijs, xj, yj, zfac, zj, zvac

    COMMON /mm_store/ efintm,efoutr
    DATA ALAM/1.3_real_8/ ! !parameter to control smearing of ESP charge.
    ! ==--------------------------------------------------------------==
    forc(1) = 0._real_8
    forc(2) = 0._real_8
    forc(3) = 0._real_8
    eext = 0._real_8
    ! 
    DO inr=1,epot1%nnear   ! NNEAR= NO OF NEAR MM ATOMS TO THE QM SYSTEM
       idnr=epot1%idnear(inr)
       xj=epot2%konear(1,inr)
       yj=epot2%konear(2,inr)
       zj=epot2%konear(3,inr)
       ! 
       ALFA=ALAM*epot2%datnear(3,INR) ! !DATNEAR(3,IT) is column-6 of .run
       beta=2._real_8*alfa
       ! 
       rijs = (xj-xi)**2+(yj-yi)**2+(zj-zi)**2
       rij = SQRT(rijs)
       rijq = rij*rijs
       ! 
       zfac=1._real_8/rijq&
            - EXP(-beta*rij)*(1._real_8/rijq+beta/rijs+alfa*beta/rij)
       zvac=1._real_8/rij - EXP(-beta*rij)*(1._real_8/rij + alfa)
       ! ..   
       forc(1)=forc(1)+epot2%datnear(1,inr)*zvp*(xi-xj)*zfac
       ! ..       
       forc(2)=forc(2)+epot2%datnear(1,inr)*zvp*(yi-yj)*zfac
       ! ..   
       forc(3)=forc(3)+epot2%datnear(1,inr)*zvp*(zi-zj)*zfac
       ! ..
       eext=eext+epot2%datnear(1,inr)*zvp*zvac
       ! ..
    ENDDO
    ! ==--------------------------------------------------------------==      
    ! ..
    IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
       eintm   = 0.0_real_8
       fintm_1  = 0.0_real_8
       fintm_2  = 0.0_real_8
       fintm_3  = 0.0_real_8
       DO inr=1,epot3%nintm        ! NINTM= NO OF INTERMEDIATE MM ATOMS
          idnr=epot3%idintm(inr)
          xj=epot4%kointm(1,inr)
          yj=epot4%kointm(2,inr)
          zj=epot4%kointm(3,inr)
          ! 
          alfa = alam*epot4%qintm(3,inr)
          beta = 2._real_8*alfa
          ! 
          rijs = (xj-xi)**2+(yj-yi)**2+(zj-zi)**2
          rij = SQRT(rijs)
          rijq = rij*rijs
          ! 
          zfac=1._real_8/rijq&
               - EXP(-beta*rij)*(1._real_8/rijq+beta/rijs+alfa*beta/rij)
          zvac=1._real_8/rij - EXP(-beta*rij)*(1._real_8/rij + alfa)
          ! ..   
          fintm_1 = fintm_1 + epot4%qintm(1,inr)*zvp*(xi-xj)*zfac
          ! ..   
          fintm_2 = fintm_2 + epot4%qintm(1,inr)*zvp*(yi-yj)*zfac
          ! ..   
          fintm_3 = fintm_3 + epot4%qintm(1,inr)*zvp*(zi-zj)*zfac
          ! ..   
          eintm = eintm + epot4%qintm(1,inr)*zvp*zvac
          ! ..
       ENDDO
       ! ..
       efintm(1,is,ia) = fintm_1
       efintm(2,is,ia) = fintm_2
       efintm(3,is,ia) = fintm_3
       efintm(4,is,ia) = eintm
       ! .. 
    ENDIF
    ! ..      
    ! ==--------------------------------------------------------------==      
    ! ..     
    IF (MOD(epot3%stpno,epot3%outmf).EQ.0.AND.epot3%noutr.GT.0) THEN
       eoutr   = 0.0_real_8
       foutr_1  = 0.0_real_8
       foutr_2  = 0.0_real_8
       foutr_3  = 0.0_real_8
       DO inr=1,epot3%noutr        ! NOUTR= NO OF OUTERMOST MM ATOMS
          idnr=epot3%idoutr(inr)
          xj=epot4%kooutr(1,inr)
          yj=epot4%kooutr(2,inr)
          zj=epot4%kooutr(3,inr)
          ! 
          alfa = alam*epot4%qoutr(3,inr)
          beta = 2._real_8*alfa
          ! 
          rijs = (xj-xi)**2+(yj-yi)**2+(zj-zi)**2
          rij = SQRT(rijs)
          rijq = rij*rijs
          ! 
          zfac=1._real_8/rijq&
               - EXP(-beta*rij)*(1._real_8/rijq+beta/rijs+alfa*beta/rij)
          zvac=1._real_8/rij - EXP(-beta*rij)*(1._real_8/rij + alfa)
          ! ..   
          foutr_1 = foutr_1 + epot4%qoutr(1,inr)*zvp*(xi-xj)*zfac
          ! ..   
          foutr_2 = foutr_2 + epot4%qoutr(1,inr)*zvp*(yi-yj)*zfac
          ! ..   
          foutr_3 = foutr_3 + epot4%qoutr(1,inr)*zvp*(zi-zj)*zfac
          ! ..
          eoutr = eoutr + epot4%qoutr(1,inr)*zvp*zvac
          ! ..
       ENDDO
       ! ..
       efoutr(1,is,ia) = foutr_1
       efoutr(2,is,ia) = foutr_2
       efoutr(3,is,ia) = foutr_3
       efoutr(4,is,ia) = eoutr
       ! .. 
    ENDIF
    ! ..      
    forc(1)=forc(1)+efintm(1,is,ia)+efoutr(1,is,ia)
    forc(2)=forc(2)+efintm(2,is,ia)+efoutr(2,is,ia)
    forc(3)=forc(3)+efintm(3,is,ia)+efoutr(3,is,ia)
    eext=eext+efintm(4,is,ia)+efoutr(4,is,ia)
    ! ..
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vextnub


END MODULE eextern_utils
