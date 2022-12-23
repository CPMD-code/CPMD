MODULE potmed_utils
  USE adat,                            ONLY: elem
  USE atomc_utils,                     ONLY: positx
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: potmed
  !public :: getpot
  PUBLIC :: printpot

CONTAINS

  ! ==================================================================
  SUBROUTINE potmed(rhoe,rcut,tau0,achrg,ichrg)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES ATOMIC AVERAGED ELECTROSTATIC POTENTIAL INSIDE RC   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: rhoe(fpar%nnr1), rcut, tau0(3,maxsys%nax,maxsys%nsx), &
      achrg(ions1%nat), ichrg(ions1%nat)

    INTEGER                                  :: ia, iat, igrid, is, isub
    REAL(real_8)                             :: chrg

    CALL tiset('    POTMED',isub)
    ! ==--------------------------------------------------------------==
    ! ..initialization
    ! ==--------------------------------------------------------------==
    ! ..charge integration
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL getpot(rhoe,tau0(1,ia,is),rcut,chrg,igrid)
          achrg(iat)=chrg
          ichrg(iat)=REAL(igrid,kind=real_8)
       ENDDO
    ENDDO
    CALL mp_sum(achrg,ions1%nat,parai%allgrp)
    CALL mp_sum(ichrg,ions1%nat,parai%allgrp)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          chrg=achrg(iat)/ichrg(iat)
          achrg(iat)=chrg
       ENDDO
    ENDDO
    CALL tihalt('    POTMED',isub)
    ! ==--------------------------------------------------------------==
    ! Deallocation of local variables
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE potmed
  ! ==================================================================
  SUBROUTINE getpot(rhoe,tau,rcut,chrg,igrid)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%kr1,fpar%kr2,*), &
                                                tau(3), rcut, chrg
    INTEGER                                  :: igrid

    INTEGER                                  :: i, i1, ii, iorgx, iorgy, &
                                                iorgz, ix, j, j1, jj, k, k1, &
                                                kk, mr(3), nptau(3)
    REAL(real_8)                             :: da1(3,3), r, r2, sc1(3), &
                                                sc2(3), x, x1, xptau(3), y, &
                                                y1, z, z1

! Variables
! ==--------------------------------------------------------------==
! ..Atom position

    chrg=0
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
    chrg=0.0_real_8
    igrid=0
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
                r=SQRT(r2)
                IF (r.LT.rcut) THEN
                   igrid=igrid+1
                   chrg=chrg+rhoe(ii,jj,kk)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getpot
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE printpot(tau0,achrg)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), achrg(ions1%nat)

    INTEGER                                  :: ia, iat, is, k

    IF (paral%io_parent)&
         WRITE(6,'(/,1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(A,10X,A,/,15X,A)')&
         '    ATOM  ','     COORDINATES                   AVERAGED POT.'&
         ,'     X         Y         Z           POT'
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          IF (paral%io_parent)&
               WRITE(6,'(1X,I4,1X,A2,5X,3F10.4,2X,1F10.5)')&
               iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3),achrg(iat)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"),/)')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printpot
  ! ==================================================================


END MODULE potmed_utils
