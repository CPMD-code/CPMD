MODULE chain_dr_utils
  USE cnst_dyn,                        ONLY: bbeta,&
                                             lupd,&
                                             nlong,&
                                             nshort,&
                                             optdist,&
                                             rch,&
                                             rcw
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mm_dimmod,                       ONLY: clsaabox
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: paral
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: chain_dr
  !public :: chain_int
  !public :: inv3
  !public :: r_to_s
  !public :: s_to_r

CONTAINS

  ! ==================================================================
  SUBROUTINE chain_dr(nres1,naa,nres2,nh,ex1_h,ex2_h,r0_h,&
       ex1_a,ex2_a,r0_a,iatom,tscr,an,&
       lsk,ptot)
    ! ==--------------------------------------------------------------==
    ! == Calculates the order parameters that estimates the amount of ==
    ! == hydrogen bond chains in the system (PTOT) and the derivatives==
    ! == with respect to the positions of the ions involved (AN)      ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nres1, naa, nres2, nh, ex1_h, &
                                                ex2_h
    REAL(real_8)                             :: r0_h
    INTEGER                                  :: ex1_a, ex2_a
    REAL(real_8)                             :: r0_a
    INTEGER                                  :: iatom(*)
    REAL(real_8)                             :: tscr(3,*), an(*)
    INTEGER                                  :: lsk(3,*)
    REAL(real_8)                             :: ptot

    CHARACTER(*), PARAMETER                  :: procedureN = 'chain_dr'

    INTEGER                                  :: i, ia, ierr, iia, isub, ix, &
                                                l1, l2, l3, natom, nox
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: cell_ch(9)
    REAL(real_8), ALLOCATABLE                :: force(:,:), pos(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('  CHAIN_DR',isub)

    ! Allocation of memory
    IF (ifirst .EQ. 0) THEN
       ALLOCATE(pos(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(pos)!,3*ions1%nat)
       ALLOCATE(force(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(force)!,3*ions1%nat)
       ifirst = 1
    ENDIF

    ! Cell parameters for PBC
    IF (lqmmm%qmmm)THEN
       i=0
       DO ix=1,3
          cell_ch(ix+i)=clsaabox%box_au(ix)
          i = i+3
       ENDDO
    ELSE
       DO i=1,9
          cell_ch(i)=metr_com%ht(i,1)
       ENDDO
    ENDIF

    ! Initialization of the positions array
    DO ia =1,nres1+naa+nres2+nh
       iia=iatom(ia)
       DO ix=1,3
          pos(ix,ia)=tscr(ix,iia)
       ENDDO
    ENDDO


    nox  = nres1+naa
    natom=nox+nres2+nh

    ! Calculation of H-bond chains (length from 1 to 5)
    CALL chain_int(natom,nox,nres1,nres2,ex1_h,ex2_h,r0_h,&
         ex1_a,ex2_a,r0_a,pos,force,cell_ch,ptot)

    ! Derivatives wrt to the ions involved 
    DO ia = 1,natom
       iia = iatom(ia)
       l1 = lsk(1,iia)
       l2 = lsk(2,iia)
       l3 = lsk(3,iia)
       IF (l1 .NE. 0) an(l1) = force(1,ia)
       IF (l2 .NE. 0) an(l2) = force(2,ia)
       IF (l3 .NE. 0) an(l3) = force(3,ia)

    ENDDO

    CALL tihalt('  CHAIN_DR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chain_dr
  ! ==================================================================
  SUBROUTINE chain_int(natom,nox,nres1,nres2,nh,mh,rh,no,mo,ro,&
       pos,force,hpass,ptot)
    ! ==--------------------------------------------------------------==


    INTEGER                                  :: natom, nox, nres1, nres2, nh, &
                                                mh
    REAL(real_8)                             :: rh
    INTEGER                                  :: no, mo
    REAL(real_8)                             :: ro, pos(3,natom), &
                                                force(3,natom), hpass(3,3), &
                                                ptot

    CHARACTER(*), PARAMETER                  :: procedureN = 'chain_int'
    INTEGER, PARAMETER                       :: nc = 5, nhmax = 5, nomax = 15
    REAL(KIND=real_8), PARAMETER             :: small = 1.0e-34_real_8

    INTEGER                                  :: a0, a4, a5, aa1, aa2, aa3, i, &
                                                ierr, ii, is, ix, j, jj, k, &
                                                kf, kk, l, no2, p1, p2, p3, &
                                                p4, p5
    INTEGER, ALLOCATABLE                     :: hlist(:,:), nhlist(:,:,:), &
                                                nwlist(:,:), wlist(:)
    INTEGER, SAVE                            :: it = 0
    REAL(real_8) :: a, b, cfact, d, d0, den, dfunc, dfunc1, dij(3), dm, dn, &
      eve, func, func1, func2, h, hinv, MaxP(nc), num, OrderParameter(nc), &
      pt1, pt2, pt3, pt4, pt5, Ptmp1, Ptmp2, Ptmp3, Ptmp4, Ptmp5, &
      r(4,natom,natom), Rij, rrij(3), sij(3), tmp, tmp1, tmp2, tmp3, tmp4, &
      tmp5, tmpf
    REAL(real_8), ALLOCATABLE                :: dsw(:,:,:,:), dswf(:,:,:,:), &
                                                f(:,:,:), sw(:,:), swf(:,:)

    COMMON/metr_my/h(3,3),hinv(3,3)

    ! ==--------------------------------------------------------------==
    it=it+1
    IF (it.EQ.1) THEN
       no2 = nox+nres2
       ALLOCATE(wlist(no2+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hlist(nox+nres2,nox+nres2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(nwlist(nhmax,nox+nres2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(nhlist(nhmax,nox+nres2,nox+nres2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(sw(nomax,nox),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(swf(nomax,nox),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dsw(9,nhmax,nomax,nox),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dswf(9,nhmax,nomax,nox),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(f(3,natom,nc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    DO i=1,3
       h(1,i)=hpass(1,i)
       h(2,i)=hpass(2,i)
       h(3,i)=hpass(3,i)
    ENDDO
    CALL inv3(h,hinv)

    ! ==--------------------------------------------------------------==
    ! Neighbours list for the center of mass of each Ox pair

    tmp=0
    IF (MOD(it,lupd).EQ.1)THEN
       CALL zeroing(wlist)!,no2)
       CALL zeroing(hlist)!,no2*no2)
       CALL zeroing(nwlist)!,nomax+no2)
       CALL zeroing(nhlist)!,nhmax+no2+no2)

       is=0

       DO i=1,nox+nres2
          k=0
          DO j=1,nox+nres2
             IF ( i .NE. j) THEN
                DO ix=1,3
                   dij(ix)=pos(ix,i)-pos(ix,j)
                ENDDO
                CALL r_to_s(dij,sij)
                DO ix=1,3
                   sij(ix)=sij(ix)-ANINT(sij(ix))
                ENDDO
                CALL s_to_r(sij,rrij)
                Rij=SQRT(rrij(1)*rrij(1)+rrij(2)*rrij(2)+rrij(3)*rrij(3))
                r(4,i,j)=Rij
                DO ix=1,3
                   r(ix,j,i) = rrij(ix)
                ENDDO
                IF (Rij.LT.Rcw) THEN
                   k=k+1
                   wlist(i)=wlist(i)+1
                   nwlist(k,i)=j
                ENDIF
             ENDIF
          ENDDO
          is=MAX(is,wlist(i))
       ENDDO
       IF (is.GT.nomax) THEN
          CALL stopgm('CHAIN_INT','too many O neighbors',& 
               __LINE__,__FILE__)
       ENDIF
       is=0
       DO i=1,nox+nres2
          DO j=nox+nres2+1,natom
             DO ix=1,3
                dij(ix)=pos(ix,i)-pos(ix,j)
             ENDDO
             CALL r_to_s(dij,sij)
             DO ix=1,3
                sij(ix)=sij(ix)-ANINT(sij(ix))
             ENDDO
             CALL s_to_r(sij,rrij)
             Rij=SQRT(rrij(1)**2+rrij(2)**2+rrij(3)**2)
             r(4,j,i)=Rij
             r(4,i,j)=Rij
             DO ix=1,3
                r(ix,j,i) = rrij(ix)
                r(ix,i,j) = -rrij(ix)
             ENDDO
          ENDDO
       ENDDO

       ! ==------------------------------------------------------------==
       ! List of H atoms w.r.t. center of mass of Ox
       DO i=1,nox+nres2
          DO ii=1,wlist(i)
             j=nwlist(ii,i)
             l=0
             DO k=nox+nres2+1,natom
                DO ix=1,3
                   dij(ix)=pos(ix,i)-pos(ix,j)
                ENDDO
                CALL r_to_s(dij,sij)
                DO ix=1,3
                   sij(ix)=sij(ix)-ANINT(sij(ix))
                ENDDO
                CALL s_to_r(sij,rrij)
                DO ix=1,3
                   dij(ix)=pos(ix,k)-pos(ix,j)-rrij(ix)*0.5_real_8
                ENDDO
                CALL r_to_s(dij,sij)
                DO ix=1,3
                   sij(ix)=sij(ix)-ANINT(sij(ix))
                ENDDO
                CALL s_to_r(sij,rrij)
                Rij=SQRT(rrij(1)**2+rrij(2)**2+rrij(3)**2)
                IF (Rij.LT.Rch) THEN
                   l=l+1
                   hlist(j,i)=l
                   nhlist(l,j,i)=k
                   hlist(i,j)=l
                   nhlist(l,i,j)=k
                ENDIF
             ENDDO
             is=MAX(is,hlist(i,j))
          ENDDO
       ENDDO
       IF (is.GT.nhmax) THEN
          CALL stopgm('CHAIN_INT','too many H neighbors',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       DO i =1,nox+nres2
          DO jj=1,wlist(i)
             j=nwlist(jj,i)
             DO ix=1,3
                dij(ix)=pos(ix,i)-pos(ix,j)
             ENDDO
             CALL r_to_s(dij,sij)
             DO ix=1,3
                sij(ix)=sij(ix)-ANINT(sij(ix))
             ENDDO
             CALL s_to_r(sij,rrij)
             Rij=SQRT(rrij(1)**2+rrij(2)**2+rrij(3)**2)
             r(4,j,i)=Rij
             r(4,i,j)=Rij
             DO ix=1,3
                r(ix,j,i) = rrij(ix)
                r(ix,i,j) = -rrij(ix)
             ENDDO
             DO kk=1,hlist(j,i)
                k=nhlist(kk,j,i)
                DO ix=1,3
                   dij(ix)=pos(ix,i)-pos(ix,k)
                ENDDO
                CALL r_to_s(dij,sij)
                DO ix=1,3
                   sij(ix)=sij(ix)-ANINT(sij(ix))
                ENDDO
                CALL s_to_r(sij,rrij)
                Rij=SQRT(rrij(1)**2+rrij(2)**2+rrij(3)**2)
                r(4,k,i)=Rij
                r(4,i,k)=Rij
                DO ix=1,3
                   r(ix,k,i) = rrij(ix)
                   r(ix,i,k) = -rrij(ix)
                   dij(ix)=pos(ix,k)-pos(ix,j)
                ENDDO
                CALL r_to_s(dij,sij)
                DO ix=1,3
                   sij(ix)=sij(ix)-ANINT(sij(ix))
                ENDDO
                CALL s_to_r(sij,rrij)
                Rij=SQRT(rrij(1)**2+rrij(2)**2+rrij(3)**2)
                r(4,j,k)=Rij
                r(4,k,j)=Rij
                DO ix=1,3
                   r(ix,j,k) = rrij(ix)
                   r(ix,k,j) = -rrij(ix)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Compute the switching function on every H2O and active residues
    CALL zeroing(sw)!,nomax+nox)
    CALL zeroing(swf)!,nox+nres2+nres2)
    CALL zeroing(dsw)!,9+nhmax+nomax+nox)
    CALL zeroing(dswf)!,9+nhmax+nox+nres2+nres2)
    CALL zeroing(f)!,3+natom+nc)
    CALL zeroing(force)!,3*natom)
    tmp=0._real_8

    DO 100 i=1,nox
       DO j=1,wlist(i)
          l=nwlist(j,i)
          d0=r(4,l,i)-optdist
          d=d0/ro
          dn=d**no
          dm=d**mo
          num=1._real_8-dn
          den=1._real_8-dm
          func1=num/den

          dfunc1=(-dn*no/den +num/den/den *dm*mo)/d/ro/r(4,l,i)

          func2=0

          DO k=1,hlist(i,l)
             ii=nhlist(k,i,l)
             d0=r(4,ii,i)+r(4,ii,l)-r(4,l,i)
             d=d0/rh
             dn=d**nh
             dm=d**mh
             num=1._real_8-dn
             den=1._real_8-dm
             func=num/den

             func2=func2+func

             dfunc=(-dn*nh/den +num/den/den *dm*mh)/d/rh

             dsw(1,k,j,i)=r(1,l,i)*dfunc1*func  +&
                  r(1,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(1,l ,i)*func1*dfunc/r(4,l,i)
             dsw(2,k,j,i)=r(2,l,i)*dfunc1*func  +&
                  r(2,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(2,l ,i)*func1*dfunc/r(4,l,i)
             dsw(3,k,j,i)=r(3,l,i)*dfunc1*func  +&
                  r(3,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(3,l ,i)*func1*dfunc/r(4,l,i)

             dsw(4,k,j,i)=r(1,i,l)*dfunc1*func  +&
                  r(1,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(1,i ,l)*func1*dfunc/r(4,i,l)
             dsw(5,k,j,i)=r(2,i,l)*dfunc1*func  +&
                  r(2,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(2,i ,l)*func1*dfunc/r(4,i,l)
             dsw(6,k,j,i)=r(3,i,l)*dfunc1*func  +&
                  r(3,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(3,i ,l)*func1*dfunc/r(4,i,l)

             dsw(7,k,j,i)=r(1,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(1,l,ii)*func1*dfunc/r(4,l,ii)
             dsw(8,k,j,i)=r(2,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(2,l,ii)*func1*dfunc/r(4,l,ii)
             dsw(9,k,j,i)=r(3,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(3,l,ii)*func1*dfunc/r(4,l,ii)

          ENDDO             ! vectorize do k=1,hlist(i,l)

          ! write(6,*) 'func ',i,j,func1,func2,dfunc
          a=func1*func2
          sw(j,i)=a
          b=a**20
          ! write(6,*) a,b,a-b,1-b
          sw(j,i)=(a-b)/(1-b)
          cfact=(a-20*b+20*(a-b)*b/(1-b))/a/(1-b)
          DO k=1,hlist(i,l)
             dsw(1,k,j,i)=dsw(1,k,j,i)*cfact
             dsw(2,k,j,i)=dsw(2,k,j,i)*cfact
             dsw(3,k,j,i)=dsw(3,k,j,i)*cfact
             dsw(4,k,j,i)=dsw(4,k,j,i)*cfact
             dsw(5,k,j,i)=dsw(5,k,j,i)*cfact
             dsw(6,k,j,i)=dsw(6,k,j,i)*cfact
             dsw(7,k,j,i)=dsw(7,k,j,i)*cfact
             dsw(8,k,j,i)=dsw(8,k,j,i)*cfact
             dsw(9,k,j,i)=dsw(9,k,j,i)*cfact
          ENDDO

          tmp=tmp+sw(j,i)

       ENDDO                 ! do j=1,wlist(i)
100 CONTINUE

    ! run over all the last terms
    DO jj=1,nres2
       i=nox+jj
       DO l=1,nox
          d0=r(4,l,i)-optdist
          d=d0/ro
          dn=d**no
          dm=d**mo
          num=1._real_8-dn
          den=1._real_8-dm
          func1=num/den
          dfunc1=(-dn*no/den +num/den/den *dm*mo)/d/ro/r(4,l,i)

          func2=0

          DO k=1,hlist(i,l)
             ii=nhlist(k,i,l)
             d0=r(4,ii,i)+r(4,ii,l)-r(4,l,i)
             d=d0/rh
             dn=d**nh
             dm=d**mh
             num=1._real_8-dn
             den=1._real_8-dm
             func=num/den
             func2=func2+func

             dfunc=(-dn*nh/den +num/den/den *dm*mh)/d/rh

             dswf(1,k,l,jj)=r(1,l,i)*dfunc1*func  +&
                  r(1,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(1,l ,i)*func1*dfunc/r(4,l,i)
             dswf(2,k,l,jj)=r(2,l,i)*dfunc1*func  +&
                  r(2,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(2,l ,i)*func1*dfunc/r(4,l,i)
             dswf(3,k,l,jj)=r(3,l,i)*dfunc1*func  +&
                  r(3,ii,i)*func1*dfunc/r(4,ii,i)-&
                  r(3,l ,i)*func1*dfunc/r(4,l,i)

             dswf(4,k,l,jj)=r(1,i,l)*dfunc1*func  +&
                  r(1,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(1,i ,l)*func1*dfunc/r(4,i,l)
             dswf(5,k,l,jj)=r(2,i,l)*dfunc1*func  +&
                  r(2,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(2,i ,l)*func1*dfunc/r(4,i,l)
             dswf(6,k,l,jj)=r(3,i,l)*dfunc1*func  +&
                  r(3,ii,l)*func1*dfunc/r(4,ii,l)-&
                  r(3,i ,l)*func1*dfunc/r(4,i,l)

             dswf(7,k,l,jj)=r(1,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(1,l,ii)*func1*dfunc/r(4,l,ii)
             dswf(8,k,l,jj)=r(2,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(2,l,ii)*func1*dfunc/r(4,l,ii)
             dswf(9,k,l,jj)=r(3,i,ii)*func1*dfunc/r(4,i,ii)+&
                  r(3,l,ii)*func1*dfunc/r(4,l,ii)

          ENDDO

          ! SWF(l,jj)=func2*func1
          a=func1*func2
          swf(l,jj)=a
          b=a**20
          swf(l,jj)=(a-b)/(1-b)
          cfact=(a-20*b+20*(a-b)*b/(1-b))/a/(1-b)
          DO k=1,hlist(i,l)
             dswf(1,k,l,jj)=dswf(1,k,l,jj)*cfact
             dswf(2,k,l,jj)=dswf(2,k,l,jj)*cfact
             dswf(3,k,l,jj)=dswf(3,k,l,jj)*cfact
             dswf(4,k,l,jj)=dswf(4,k,l,jj)*cfact
             dswf(5,k,l,jj)=dswf(5,k,l,jj)*cfact
             dswf(6,k,l,jj)=dswf(6,k,l,jj)*cfact
             dswf(7,k,l,jj)=dswf(7,k,l,jj)*cfact
             dswf(8,k,l,jj)=dswf(8,k,l,jj)*cfact
             dswf(9,k,l,jj)=dswf(9,k,l,jj)*cfact
          ENDDO
          ! tmp=tmp+swf(l,jj)
       ENDDO
    ENDDO                     ! do jj=1,nres2

    ! ==--------------------------------------------------------------==
    ! compute the order parameter d for the chain (P)
    ! start from residue 1 and go to residue 2
    ! WARNING: the chain must be symmetric
    pt1=0._real_8
    pt2=0._real_8
    pt3=0._real_8
    pt4=0._real_8
    pt5=0._real_8
    ! PT6=0._real_8
    CALL zeroing(MaxP)!,nc)
    CALL zeroing(f)!,3+natom+nc)

    DO 5000 l=1,nres1
       a0=l
       DO 4000 p1=1,wlist(a0)
          aa1=nwlist(p1,a0)
          IF (ABS(sw(p1,a0)) .LT. small) THEN
             GOTO 5000
          ENDIF
          Ptmp1=sw(p1,a0)
          DO k=1,nres2
             ii=nox+k
             tmp=swf(aa1,k)
             IF (aa1 .EQ. ii) THEN
                CALL stopgm('CHAIN_INT','atoms too close: no chain',& 
                     __LINE__,__FILE__)
             ENDIF
             eve=Ptmp1*tmp
             pt1 = pt1 + EXP(bbeta*eve)
             eve=eve*EXP(bbeta*eve)
             ! ompute forces
             tmp1=sw(p1,a0)
             tmpf=eve/tmp1
             DO kf=1,hlist(a0,aa1)
                kk=nhlist(kf,a0,aa1)
                f(1,a0,1) =f(1,a0,1) +dsw(1,kf,p1,a0)*tmpf
                f(2,a0,1) =f(2,a0,1) +dsw(2,kf,p1,a0)*tmpf
                f(3,a0,1) =f(3,a0,1) +dsw(3,kf,p1,a0)*tmpf
                f(1,aa1,1)=f(1,aa1,1)+dsw(4,kf,p1,a0)*tmpf
                f(2,aa1,1)=f(2,aa1,1)+dsw(5,kf,p1,a0)*tmpf
                f(3,aa1,1)=f(3,aa1,1)+dsw(6,kf,p1,a0)*tmpf
                f(1,kk,1) =f(1,kk,1) +dsw(7,kf,p1,a0)*tmpf
                f(2,kk,1) =f(2,kk,1) +dsw(8,kf,p1,a0)*tmpf
                f(3,kk,1) =f(3,kk,1) +dsw(9,kf,p1,a0)*tmpf
             ENDDO
             tmpf=eve/swf(aa1,k)
             DO kf=1,hlist(aa1,ii)
                kk=nhlist(kf,aa1,ii)
                f(1,ii,1) =f(1,ii,1) +dswf(1,kf,aa1,k)*tmpf
                f(2,ii,1) =f(2,ii,1) +dswf(2,kf,aa1,k)*tmpf
                f(3,ii,1) =f(3,ii,1) +dswf(3,kf,aa1,k)*tmpf
                f(1,aa1,1)=f(1,aa1,1)+dswf(4,kf,aa1,k)*tmpf
                f(2,aa1,1)=f(2,aa1,1)+dswf(5,kf,aa1,k)*tmpf
                f(3,aa1,1)=f(3,aa1,1)+dswf(6,kf,aa1,k)*tmpf
                f(1,kk,1) =f(1,kk,1) +dswf(7,kf,aa1,k)*tmpf
                f(2,kk,1) =f(2,kk,1) +dswf(8,kf,aa1,k)*tmpf
                f(3,kk,1) =f(3,kk,1) +dswf(9,kf,aa1,k)*tmpf
             ENDDO

             IF (Ptmp1*tmp.GT.MaxP(1)) THEN
                MaxP(1)=Ptmp1*tmp
                ! AK 2005/03/26: not used.
                ! ch(1,1)=a0
                ! ch(2,1)=aa1
                ! ch(3,1)=ii  
             ENDIF
          ENDDO
          IF (nlong.EQ.1) THEN
             GOTO 4000
          ENDIF

          DO 101 p2=1,wlist(aa1)
             aa2=nwlist(p2,aa1)! take a neighbour Ox
             IF (aa2 .EQ. a0) THEN
                GOTO 101
             ENDIF
             IF (ABS(sw(p2,aa1)) .LT. small) THEN
                GOTO 101
             ENDIF
             Ptmp2=Ptmp1*sw(p2,aa1)
             DO k=1,nres2
                ii=nox+k
                tmp=swf(aa2,k)
                IF (aa1 .EQ. ii) Ptmp2=0._real_8
                eve=Ptmp2*tmp
                pt2 = pt2 + EXP(bbeta*eve)
                eve=eve*EXP(bbeta*eve)
                ! ompute forces 
                tmpf=eve/tmp1
                DO kf=1,hlist(a0,aa1)
                   kk=nhlist(kf,a0,aa1)
                   f(1,a0,2) =f(1,a0,2) +dsw(1,kf,p1,a0)*tmpf
                   f(2,a0,2) =f(2,a0,2) +dsw(2,kf,p1,a0)*tmpf
                   f(3,a0,2) =f(3,a0,2) +dsw(3,kf,p1,a0)*tmpf
                   f(1,aa1,2)=f(1,aa1,2)+dsw(4,kf,p1,a0)*tmpf
                   f(2,aa1,2)=f(2,aa1,2)+dsw(5,kf,p1,a0)*tmpf
                   f(3,aa1,2)=f(3,aa1,2)+dsw(6,kf,p1,a0)*tmpf
                   f(1,kk,2) =f(1,kk,2) +dsw(7,kf,p1,a0)*tmpf
                   f(2,kk,2) =f(2,kk,2) +dsw(8,kf,p1,a0)*tmpf
                   f(3,kk,2) =f(3,kk,2) +dsw(9,kf,p1,a0)*tmpf
                ENDDO
                tmp2=sw(p2,aa1)
                tmpf=eve/tmp2
                DO kf=1,hlist(aa1,aa2)
                   kk=nhlist(kf,aa1,aa2)
                   f(1,aa1,2)=f(1,aa1,2)+dsw(1,kf,p2,aa1)*tmpf
                   f(2,aa1,2)=f(2,aa1,2)+dsw(2,kf,p2,aa1)*tmpf
                   f(3,aa1,2)=f(3,aa1,2)+dsw(3,kf,p2,aa1)*tmpf
                   f(1,aa2,2)=f(1,aa2,2)+dsw(4,kf,p2,aa1)*tmpf
                   f(2,aa2,2)=f(2,aa2,2)+dsw(5,kf,p2,aa1)*tmpf
                   f(3,aa2,2)=f(3,aa2,2)+dsw(6,kf,p2,aa1)*tmpf
                   f(1,kk,2) =f(1,kk,2) +dsw(7,kf,p2,aa1)*tmpf
                   f(2,kk,2) =f(2,kk,2) +dsw(8,kf,p2,aa1)*tmpf
                   f(3,kk,2) =f(3,kk,2) +dsw(9,kf,p2,aa1)*tmpf
                ENDDO
                tmpf=eve/swf(aa2,k)
                DO kf=1,hlist(aa2,ii)
                   kk=nhlist(kf,aa2,ii)
                   f(1,ii,2) =f(1,ii,2) +dswf(1,kf,aa2,k)*tmpf
                   f(2,ii,2) =f(2,ii,2) +dswf(2,kf,aa2,k)*tmpf
                   f(3,ii,2) =f(3,ii,2) +dswf(3,kf,aa2,k)*tmpf
                   f(1,aa2,2)=f(1,aa2,2)+dswf(4,kf,aa2,k)*tmpf
                   f(2,aa2,2)=f(2,aa2,2)+dswf(5,kf,aa2,k)*tmpf
                   f(3,aa2,2)=f(3,aa2,2)+dswf(6,kf,aa2,k)*tmpf
                   f(1,kk,2) =f(1,kk,2) +dswf(7,kf,aa2,k)*tmpf
                   f(2,kk,2) =f(2,kk,2) +dswf(8,kf,aa2,k)*tmpf
                   f(3,kk,2) =f(3,kk,2) +dswf(9,kf,aa2,k)*tmpf
                ENDDO

                IF (Ptmp2*tmp.GT.MaxP(2)) THEN
                   MaxP(2)=Ptmp2*tmp
                   ! AK 2005/03/26: not used.
                   ! ch(1,2)=a0
                   ! ch(2,2)=aa1
                   ! ch(3,2)=aa2
                   ! ch(4,2)=ii
                ENDIF
             ENDDO
             IF (nlong.EQ.2) THEN
                GOTO 101
             ENDIF

             DO 102 p3=1,wlist(aa2)
                aa3=nwlist(p3,aa2)! prendo un Ox suo vicino
                IF (aa3 .EQ. a0) THEN
                   GOTO 102
                ENDIF
                IF (aa3 .EQ. aa1) THEN
                   GOTO 102
                ENDIF
                IF (ABS(sw(p3,aa2)) .LT. small) THEN
                   GOTO 102
                ENDIF
                Ptmp3=Ptmp2*sw(p3,aa2)
                DO k=1,nres2
                   ii=nox+k
                   tmp=swf(aa3,k)
                   IF (aa2 .EQ. ii) THEN
                      Ptmp3=0._real_8
                   ENDIF
                   eve=Ptmp3*tmp
                   pt3 = pt3 + EXP(bbeta*eve)
                   eve=eve*EXP(bbeta*eve)
                   ! ompute forces
                   tmpf=eve/tmp1
                   DO kf=1,hlist(a0,aa1)
                      kk=nhlist(kf,a0,aa1)
                      f(1,a0,3) =f(1,a0,3) +dsw(1,kf,p1,a0)*tmpf
                      f(2,a0,3) =f(2,a0,3) +dsw(2,kf,p1,a0)*tmpf
                      f(3,a0,3) =f(3,a0,3) +dsw(3,kf,p1,a0)*tmpf
                      f(1,aa1,3)=f(1,aa1,3)+dsw(4,kf,p1,a0)*tmpf
                      f(2,aa1,3)=f(2,aa1,3)+dsw(5,kf,p1,a0)*tmpf
                      f(3,aa1,3)=f(3,aa1,3)+dsw(6,kf,p1,a0)*tmpf
                      f(1,kk,3) =f(1,kk,3) +dsw(7,kf,p1,a0)*tmpf
                      f(2,kk,3) =f(2,kk,3) +dsw(8,kf,p1,a0)*tmpf
                      f(3,kk,3) =f(3,kk,3) +dsw(9,kf,p1,a0)*tmpf
                   ENDDO
                   tmpf=eve/tmp2
                   DO kf=1,hlist(aa1,aa2)
                      kk=nhlist(kf,aa1,aa2)
                      f(1,aa1,3)=f(1,aa1,3)+dsw(1,kf,p2,aa1)*tmpf
                      f(2,aa1,3)=f(2,aa1,3)+dsw(2,kf,p2,aa1)*tmpf
                      f(3,aa1,3)=f(3,aa1,3)+dsw(3,kf,p2,aa1)*tmpf
                      f(1,aa2,3)=f(1,aa2,3)+dsw(4,kf,p2,aa1)*tmpf
                      f(2,aa2,3)=f(2,aa2,3)+dsw(5,kf,p2,aa1)*tmpf
                      f(3,aa2,3)=f(3,aa2,3)+dsw(6,kf,p2,aa1)*tmpf
                      f(1,kk,3) =f(1,kk,3) +dsw(7,kf,p2,aa1)*tmpf
                      f(2,kk,3) =f(2,kk,3) +dsw(8,kf,p2,aa1)*tmpf
                      f(3,kk,3) =f(3,kk,3) +dsw(9,kf,p2,aa1)*tmpf
                   ENDDO
                   tmp3=sw(p3,aa2)
                   tmpf=eve/tmp3
                   DO kf=1,hlist(aa2,aa3)
                      kk=nhlist(kf,aa2,aa3)
                      f(1,aa2,3)=f(1,aa2,3)+dsw(1,kf,p3,aa2)*tmpf
                      f(2,aa2,3)=f(2,aa2,3)+dsw(2,kf,p3,aa2)*tmpf
                      f(3,aa2,3)=f(3,aa2,3)+dsw(3,kf,p3,aa2)*tmpf
                      f(1,aa3,3)=f(1,aa3,3)+dsw(4,kf,p3,aa2)*tmpf
                      f(2,aa3,3)=f(2,aa3,3)+dsw(5,kf,p3,aa2)*tmpf
                      f(3,aa3,3)=f(3,aa3,3)+dsw(6,kf,p3,aa2)*tmpf
                      f(1,kk,3) =f(1,kk,3) +dsw(7,kf,p3,aa2)*tmpf
                      f(2,kk,3) =f(2,kk,3) +dsw(8,kf,p3,aa2)*tmpf
                      f(3,kk,3) =f(3,kk,3) +dsw(9,kf,p3,aa2)*tmpf
                   ENDDO
                   tmpf=eve/swf(aa3,k)
                   DO kf=1,hlist(aa3,ii)
                      kk=nhlist(kf,aa3,ii)
                      f(1,ii,3) =f(1,ii,3) +dswf(1,kf,aa3,k)*tmpf
                      f(2,ii,3) =f(2,ii,3) +dswf(2,kf,aa3,k)*tmpf
                      f(3,ii,3) =f(3,ii,3) +dswf(3,kf,aa3,k)*tmpf
                      f(1,aa3,3)=f(1,aa3,3)+dswf(4,kf,aa3,k)*tmpf
                      f(2,aa3,3)=f(2,aa3,3)+dswf(5,kf,aa3,k)*tmpf
                      f(3,aa3,3)=f(3,aa3,3)+dswf(6,kf,aa3,k)*tmpf
                      f(1,kk,3) =f(1,kk,3) +dswf(7,kf,aa3,k)*tmpf
                      f(2,kk,3) =f(2,kk,3) +dswf(8,kf,aa3,k)*tmpf
                      f(3,kk,3) =f(3,kk,3) +dswf(9,kf,aa3,k)*tmpf
                   ENDDO

                   IF (Ptmp3*tmp.GT.MaxP(3)) THEN
                      MaxP(3)=Ptmp3*tmp
                      ! AK 2005/03/26: not used.
                      ! ch(1,3)=a0
                      ! ch(2,3)=aa1
                      ! ch(3,3)=aa2
                      ! ch(4,3)=aa3
                      ! ch(5,3)=ii
                   ENDIF
                ENDDO
                IF (nlong.EQ.3) THEN
                   GOTO 102
                ENDIF
                DO 103 p4=1,wlist(aa3)
                   a4=nwlist(p4,aa3)! take an Ox next neighbour
                   IF (a4 .EQ. a0) THEN
                      GOTO 103
                   ENDIF
                   IF (a4 .EQ. aa1) THEN
                      GOTO 103
                   ENDIF
                   IF (a4 .EQ. aa2) THEN
                      GOTO 103
                   ENDIF
                   IF (ABS(sw(p4,aa3)) .LT. small) THEN
                      GOTO 103
                   ENDIF
                   Ptmp4=Ptmp3*sw(p4,aa3)
                   DO k=1,nres2
                      ii=nox+k
                      tmp=swf(a4,k)
                      IF (aa3 .EQ. ii) THEN
                         Ptmp4=0._real_8
                      ENDIF
                      eve=Ptmp4*tmp
                      pt4 = pt4 + EXP(bbeta*eve)
                      eve=eve*EXP(bbeta*eve)
                      ! ompute forces
                      tmpf=eve/tmp1
                      DO kf=1,hlist(a0,aa1)
                         kk=nhlist(kf,a0,aa1)
                         f(1,a0,4) =f(1,a0,4) +dsw(1,kf,p1,a0)*tmpf
                         f(2,a0,4) =f(2,a0,4) +dsw(2,kf,p1,a0)*tmpf
                         f(3,a0,4) =f(3,a0,4) +dsw(3,kf,p1,a0)*tmpf
                         f(1,aa1,4)=f(1,aa1,4)+dsw(4,kf,p1,a0)*tmpf
                         f(2,aa1,4)=f(2,aa1,4)+dsw(5,kf,p1,a0)*tmpf
                         f(3,aa1,4)=f(3,aa1,4)+dsw(6,kf,p1,a0)*tmpf
                         f(1,kk,4) =f(1,kk,4) +dsw(7,kf,p1,a0)*tmpf
                         f(2,kk,4) =f(2,kk,4) +dsw(8,kf,p1,a0)*tmpf
                         f(3,kk,4) =f(3,kk,4) +dsw(9,kf,p1,a0)*tmpf
                      ENDDO
                      tmpf=eve/tmp2
                      DO kf=1,hlist(aa1,aa2)
                         kk=nhlist(kf,aa1,aa2)
                         f(1,aa1,4)=f(1,aa1,4)+dsw(1,kf,p2,aa1)*tmpf
                         f(2,aa1,4)=f(2,aa1,4)+dsw(2,kf,p2,aa1)*tmpf
                         f(3,aa1,4)=f(3,aa1,4)+dsw(3,kf,p2,aa1)*tmpf
                         f(1,aa2,4)=f(1,aa2,4)+dsw(4,kf,p2,aa1)*tmpf
                         f(2,aa2,4)=f(2,aa2,4)+dsw(5,kf,p2,aa1)*tmpf
                         f(3,aa2,4)=f(3,aa2,4)+dsw(6,kf,p2,aa1)*tmpf
                         f(1,kk,4) =f(1,kk,4) +dsw(7,kf,p2,aa1)*tmpf
                         f(2,kk,4) =f(2,kk,4) +dsw(8,kf,p2,aa1)*tmpf
                         f(3,kk,4) =f(3,kk,4) +dsw(9,kf,p2,aa1)*tmpf
                      ENDDO
                      tmpf=eve/tmp3
                      DO kf=1,hlist(aa2,aa3)
                         kk=nhlist(kf,aa2,aa3)
                         f(1,aa2,4)=f(1,aa2,4)+dsw(1,kf,p3,aa2)*tmpf
                         f(2,aa2,4)=f(2,aa2,4)+dsw(2,kf,p3,aa2)*tmpf
                         f(3,aa2,4)=f(3,aa2,4)+dsw(3,kf,p3,aa2)*tmpf
                         f(1,aa3,4)=f(1,aa3,4)+dsw(4,kf,p3,aa2)*tmpf
                         f(2,aa3,4)=f(2,aa3,4)+dsw(5,kf,p3,aa2)*tmpf
                         f(3,aa3,4)=f(3,aa3,4)+dsw(6,kf,p3,aa2)*tmpf
                         f(1,kk,4) =f(1,kk,4) +dsw(7,kf,p3,aa2)*tmpf
                         f(2,kk,4) =f(2,kk,4) +dsw(8,kf,p3,aa2)*tmpf
                         f(3,kk,4) =f(3,kk,4) +dsw(9,kf,p3,aa2)*tmpf
                      ENDDO
                      tmp4=sw(p4,aa3)
                      tmpf=eve/tmp4
                      DO kf=1,hlist(aa3,a4)
                         kk=nhlist(kf,aa3,a4)
                         f(1,aa3,4)=f(1,aa3,4)+dsw(1,kf,p4,aa3)*tmpf
                         f(2,aa3,4)=f(2,aa3,4)+dsw(2,kf,p4,aa3)*tmpf
                         f(3,aa3,4)=f(3,aa3,4)+dsw(3,kf,p4,aa3)*tmpf
                         f(1,a4,4) =f(1,a4,4) +dsw(4,kf,p4,aa3)*tmpf
                         f(2,a4,4) =f(2,a4,4) +dsw(5,kf,p4,aa3)*tmpf
                         f(3,a4,4) =f(3,a4,4) +dsw(6,kf,p4,aa3)*tmpf
                         f(1,kk,4) =f(1,kk,4) +dsw(7,kf,p4,aa3)*tmpf
                         f(2,kk,4) =f(2,kk,4) +dsw(8,kf,p4,aa3)*tmpf
                         f(3,kk,4) =f(3,kk,4) +dsw(9,kf,p4,aa3)*tmpf
                      ENDDO
                      tmpf=eve/swf(a4,k)
                      DO kf=1,hlist(a4,ii)
                         kk=nhlist(kf,a4,ii)
                         f(1,ii,4)=f(1,ii,4)+dswf(1,kf,a4,k)*tmpf
                         f(2,ii,4)=f(2,ii,4)+dswf(2,kf,a4,k)*tmpf
                         f(3,ii,4)=f(3,ii,4)+dswf(3,kf,a4,k)*tmpf
                         f(1,a4,4)=f(1,a4,4)+dswf(4,kf,a4,k)*tmpf
                         f(2,a4,4)=f(2,a4,4)+dswf(5,kf,a4,k)*tmpf
                         f(3,a4,4)=f(3,a4,4)+dswf(6,kf,a4,k)*tmpf
                         f(1,kk,4)=f(1,kk,4)+dswf(7,kf,a4,k)*tmpf
                         f(2,kk,4)=f(2,kk,4)+dswf(8,kf,a4,k)*tmpf
                         f(3,kk,4)=f(3,kk,4)+dswf(9,kf,a4,k)*tmpf
                      ENDDO

                      IF (Ptmp4*tmp.GT.MaxP(4)) THEN
                         MaxP(4)=Ptmp4*tmp
                         ! AK 2005/03/26: not used.
                         ! ch(1,4)=a0
                         ! ch(2,4)=aa1
                         ! ch(3,4)=aa2
                         ! ch(4,4)=aa3
                         ! ch(5,4)=a4
                         ! ch(6,4)=ii
                      ENDIF
                   ENDDO
                   IF (nlong.EQ.4) THEN
                      GOTO 103
                   ENDIF

                   DO 104 p5=1,wlist(a4)
                      a5=nwlist(p5,a4)! prendo un Ox suo vicino
                      IF (a5 .EQ. a0) THEN
                         GOTO 104
                      ENDIF
                      IF (a5 .EQ. aa1) THEN
                         GOTO 104
                      ENDIF
                      IF (a5 .EQ. aa2) THEN
                         GOTO 104
                      ENDIF
                      IF (a5 .EQ. aa3) THEN
                         GOTO 104
                      ENDIF
                      IF (ABS(sw(p5,a4)) .LT. small) THEN
                         GOTO 104
                      ENDIF
                      Ptmp5=Ptmp4*sw(p5,a4)
                      DO k=1,nres2
                         ii=nox+k
                         tmp=swf(a5,k)
                         IF (a4 .EQ. ii) Ptmp5=0._real_8
                         eve=Ptmp5*tmp
                         pt5 = pt5 + EXP(bbeta*eve)
                         eve=eve*EXP(bbeta*eve)

                         ! ompute forces
                         tmpf=eve/tmp1
                         DO kf=1,hlist(a0,aa1)
                            kk=nhlist(kf,a0,aa1)
                            f(1,a0,5) =f(1,a0,5) +dsw(1,kf,p1,a0)*tmpf
                            f(2,a0,5) =f(2,a0,5) +dsw(2,kf,p1,a0)*tmpf
                            f(3,a0,5) =f(3,a0,5) +dsw(3,kf,p1,a0)*tmpf
                            f(1,aa1,5)=f(1,aa1,5)+dsw(4,kf,p1,a0)*tmpf
                            f(2,aa1,5)=f(2,aa1,5)+dsw(5,kf,p1,a0)*tmpf
                            f(3,aa1,5)=f(3,aa1,5)+dsw(6,kf,p1,a0)*tmpf
                            f(1,kk,5) =f(1,kk,5) +dsw(7,kf,p1,a0)*tmpf
                            f(2,kk,5) =f(2,kk,5) +dsw(8,kf,p1,a0)*tmpf
                            f(3,kk,5) =f(3,kk,5) +dsw(9,kf,p1,a0)*tmpf
                         ENDDO
                         tmpf=eve/tmp2
                         DO kf=1,hlist(aa1,aa2)
                            kk=nhlist(kf,aa1,aa2)
                            f(1,aa1,5)=f(1,aa1,5)+dsw(1,kf,p2,aa1)*tmpf
                            f(2,aa1,5)=f(2,aa1,5)+dsw(2,kf,p2,aa1)*tmpf
                            f(3,aa1,5)=f(3,aa1,5)+dsw(3,kf,p2,aa1)*tmpf
                            f(1,aa2,5)=f(1,aa2,5)+dsw(4,kf,p2,aa1)*tmpf
                            f(2,aa2,5)=f(2,aa2,5)+dsw(5,kf,p2,aa1)*tmpf
                            f(3,aa2,5)=f(3,aa2,5)+dsw(6,kf,p2,aa1)*tmpf
                            f(1,kk,5) =f(1,kk,5) +dsw(7,kf,p2,aa1)*tmpf
                            f(2,kk,5) =f(2,kk,5) +dsw(8,kf,p2,aa1)*tmpf
                            f(3,kk,5) =f(3,kk,5) +dsw(9,kf,p2,aa1)*tmpf
                         ENDDO
                         tmpf=eve/tmp3
                         DO kf=1,hlist(aa2,aa3)
                            kk=nhlist(kf,aa2,aa3)
                            f(1,aa2,5)=f(1,aa2,5)+dsw(1,kf,p3,aa2)*tmpf
                            f(2,aa2,5)=f(2,aa2,5)+dsw(2,kf,p3,aa2)*tmpf
                            f(3,aa2,5)=f(3,aa2,5)+dsw(3,kf,p3,aa2)*tmpf
                            f(1,aa3,5)=f(1,aa3,5)+dsw(4,kf,p3,aa2)*tmpf
                            f(2,aa3,5)=f(2,aa3,5)+dsw(5,kf,p3,aa2)*tmpf
                            f(3,aa3,5)=f(3,aa3,5)+dsw(6,kf,p3,aa2)*tmpf
                            f(1,kk,5) =f(1,kk,5) +dsw(7,kf,p3,aa2)*tmpf
                            f(2,kk,5) =f(2,kk,5) +dsw(8,kf,p3,aa2)*tmpf
                            f(3,kk,5) =f(3,kk,5) +dsw(9,kf,p3,aa2)*tmpf
                         ENDDO
                         tmpf=eve/tmp4
                         DO kf=1,hlist(aa3,a4)
                            kk=nhlist(kf,aa3,a4)
                            f(1,aa3,5)=f(1,aa3,5)+dsw(1,kf,p4,aa3)*tmpf
                            f(2,aa3,5)=f(2,aa3,5)+dsw(2,kf,p4,aa3)*tmpf
                            f(3,aa3,5)=f(3,aa3,5)+dsw(3,kf,p4,aa3)*tmpf
                            f(1,a4,5) =f(1,a4,5) +dsw(4,kf,p4,aa3)*tmpf
                            f(2,a4,5) =f(2,a4,5) +dsw(5,kf,p4,aa3)*tmpf
                            f(3,a4,5) =f(3,a4,5) +dsw(6,kf,p4,aa3)*tmpf
                            f(1,kk,5) =f(1,kk,5) +dsw(7,kf,p4,aa3)*tmpf
                            f(2,kk,5) =f(2,kk,5) +dsw(8,kf,p4,aa3)*tmpf
                            f(3,kk,5) =f(3,kk,5) +dsw(9,kf,p4,aa3)*tmpf
                         ENDDO
                         tmp5=sw(p5,a4)
                         tmpf=eve/tmp5
                         DO kf=1,hlist(a4,a5)
                            kk=nhlist(kf,a4,a5)
                            f(1,a4,5)=f(1,a4,5)+dsw(1,kf,p5,a4)*tmpf
                            f(2,a4,5)=f(2,a4,5)+dsw(2,kf,p5,a4)*tmpf
                            f(3,a4,5)=f(3,a4,5)+dsw(3,kf,p5,a4)*tmpf
                            f(1,a5,5)=f(1,a5,5)+dsw(4,kf,p5,a4)*tmpf
                            f(2,a5,5)=f(2,a5,5)+dsw(5,kf,p5,a4)*tmpf
                            f(3,a5,5)=f(3,a5,5)+dsw(6,kf,p5,a4)*tmpf
                            f(1,kk,5)=f(1,kk,5)+dsw(7,kf,p5,a4)*tmpf
                            f(2,kk,5)=f(2,kk,5)+dsw(8,kf,p5,a4)*tmpf
                            f(3,kk,5)=f(3,kk,5)+dsw(9,kf,p5,a4)*tmpf
                         ENDDO
                         tmpf=eve/swf(a5,k)
                         DO kf=1,hlist(a5,ii)
                            kk=nhlist(kf,a5,ii)
                            f(1,ii,5)=f(1,ii,5)+dswf(1,kf,a5,k)*tmpf
                            f(2,ii,5)=f(2,ii,5)+dswf(2,kf,a5,k)*tmpf
                            f(3,ii,5)=f(3,ii,5)+dswf(3,kf,a5,k)*tmpf
                            f(1,a5,5)=f(1,a5,5)+dswf(4,kf,a5,k)*tmpf
                            f(2,a5,5)=f(2,a5,5)+dswf(5,kf,a5,k)*tmpf
                            f(3,a5,5)=f(3,a5,5)+dswf(6,kf,a5,k)*tmpf
                            f(1,kk,5)=f(1,kk,5)+dswf(7,kf,a5,k)*tmpf
                            f(2,kk,5)=f(2,kk,5)+dswf(8,kf,a5,k)*tmpf
                            f(3,kk,5)=f(3,kk,5)+dswf(9,kf,a5,k)*tmpf
                         ENDDO

                         IF (Ptmp5*tmp.GT.MaxP(5)) THEN
                            MaxP(5)=Ptmp5*tmp
                            ! AK 2005/03/26: not used.
                            ! ch(1,5)=a0
                            ! ch(2,5)=aa1
                            ! ch(3,5)=aa2
                            ! ch(4,5)=aa3
                            ! ch(5,5)=a4
                            ! ch(6,5)=a5
                            ! ch(7,5)=ii
                         ENDIF
                      ENDDO
                      ! !$                        if(NLONG.eq.5)cycle                                 
                      ! !$
                      ! !$                        do 105 p6=1,wlist(a5)
                      ! !$                           a6=nwlist(p6,a5) !prendo un Ox suo vicino
                      ! !$                           if (a6 .eq. a0) cycle 
                      ! !$                           if (a6 .eq. aa1) cycle 
                      ! !$                           if (a6 .eq. aa2) cycle 
                      ! !$                           if (a6 .eq. aa3) cycle 
                      ! !$                           if (a6 .eq. a4) cycle 
                      ! !$                           if (sw(p6,a5) .eq. 0.) cycle 
                      ! !$                           Ptmp6=Ptmp5*sw(p6,a5)
                      ! !$                           do k=1,nres2
                      ! !$                              ii=NOx+k
                      ! !$                              tmp=swf(a6,k)
                      ! !$                              if (a5 .eq. ii) Ptmp6=0.
                      ! !$                              eve=Ptmp6*tmp
                      ! !$                              PT6 = PT6 + exp(bbeta*eve)
                      ! !$                              eve=eve*exp(bbeta*eve)
                      ! !$!     ****************************************
                      ! !$!     calcolo la forza
                      ! !$
                      ! !$                              tmpf=eve/tmp1
                      ! !$                              do kf=1,hlist(a0,aa1)
                      ! !$                                 kk=nhlist(kf,a0,aa1)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,a0,6)=F(ix,a0,6)+DSW(ix  ,kf,p1,a0)*tmpf
                      ! !$                                    F(ix,aa1,6)=F(ix,aa1,6)+DSW(ix+3,kf,p1,a0)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p1,a0)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmpf=eve/tmp2
                      ! !$                              do kf=1,hlist(aa1,aa2)
                      ! !$                                 kk=nhlist(kf,aa1,aa2)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,aa1,6)=F(ix,aa1,6)+DSW(ix  ,kf,p2,aa1)*tmpf
                      ! !$                                    F(ix,aa2,6)=F(ix,aa2,6)+DSW(ix+3,kf,p2,aa1)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p2,aa1)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmpf=eve/tmp3
                      ! !$                              do kf=1,hlist(aa2,aa3)
                      ! !$                                 kk=nhlist(kf,aa2,aa3)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,aa2,6)=F(ix,aa2,6)+DSW(ix  ,kf,p3,aa2)*tmpf
                      ! !$                                    F(ix,aa3,6)=F(ix,aa3,6)+DSW(ix+3,kf,p3,aa2)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p3,aa2)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmpf=eve/tmp4
                      ! !$                              do kf=1,hlist(aa3,a4)
                      ! !$                                 kk=nhlist(kf,aa3,a4)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,aa3,6)=F(ix,aa3,6)+DSW(ix  ,kf,p4,aa3)*tmpf
                      ! !$                                    F(ix,a4,6)=F(ix,a4,6)+DSW(ix+3,kf,p4,aa3)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p4,aa3)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmp5=SW(p5,a4)
                      ! !$                              tmpf=eve/tmp5
                      ! !$                              do kf=1,hlist(a4,a5)
                      ! !$                                 kk=nhlist(kf,a4,a5)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,a4,6)=F(ix,a4,6)+DSW(ix  ,kf,p5,a4)*tmpf
                      ! !$                                    F(ix,a5,6)=F(ix,a5,6)+DSW(ix+3,kf,p5,a4)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p5,a4)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmp6=SW(p6,a5)
                      ! !$                              tmpf=eve/tmp6
                      ! !$                              do kf=1,hlist(a5,a6)
                      ! !$                                 kk=nhlist(kf,a5,a6)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,a4,6)=F(ix,a4,6)+DSW(ix  ,kf,p6,a5)*tmpf
                      ! !$                                    F(ix,a5,6)=F(ix,a5,6)+DSW(ix+3,kf,p6,a5)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSW(ix+6,kf,p6,a5)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$                              tmpf=eve/SWF(a6,k)
                      ! !$                              do kf=1,hlist(a6,ii)
                      ! !$                                 kk=nhlist(kf,a6,ii)
                      ! !$                                 do ix=1,3
                      ! !$                                    F(ix,ii,6)=F(ix,ii,6)+DSWF(ix  ,kf,a6,k)*tmpf
                      ! !$                                    F(ix,a5,6)=F(ix,a5,6)+DSWF(ix+3,kf,a6,k)*tmpf
                      ! !$                                    F(ix,kk,6)=F(ix,kk,6)+DSWF(ix+6,kf,a6,k)*tmpf
                      ! !$                                 enddo
                      ! !$                              enddo
                      ! !$!     ****************************************
                      ! !$                              if (Ptmp6*tmp.gt.MaxP(6)) then
                      ! !$                                 MaxP(6)=Ptmp6*tmp
                      ! !$                                 ch(1,6)=a0
                      ! !$                                 ch(2,6)=aa1
                      ! !$                                 ch(3,6)=aa2
                      ! !$                                 ch(4,6)=aa3
                      ! !$                                 ch(5,6)=a4
                      ! !$                                 ch(6,6)=a5
                      ! !$                                 ch(7,6)=a6
                      ! !$                                 ch(8,6)=ii
                      ! !$                              endif
                      ! !$                           enddo
                      ! !$ 105              continue
104                CONTINUE
103             CONTINUE
102          CONTINUE
101       CONTINUE
4000   CONTINUE
5000 CONTINUE



    OrderParameter(1)=pt1
    OrderParameter(1)=pt2
    OrderParameter(3)=pt3
    OrderParameter(4)=pt4
    OrderParameter(5)=pt5
    ! OrderParameter(6)=PT6

    pt1=LOG(pt1)/bbeta
    pt2=LOG(pt2)/bbeta
    pt3=LOG(pt3)/bbeta
    pt4=LOG(pt4)/bbeta
    pt5=LOG(pt5)/bbeta
    ! PT6=log(PT6)/bbeta

    ptot=0._real_8

    DO i=nshort,nlong
       ptot=ptot+OrderParameter(i)
    ENDDO

    DO j=nshort,nlong
       DO i=1,natom
          force(1,i)=force(1,i)+f(1,i,j)/ptot
          force(2,i)=force(2,i)+f(2,i,j)/ptot
          force(3,i)=force(3,i)+f(3,i,j)/ptot
       ENDDO
    ENDDO

    ptot = LOG(ptot)/bbeta

    DEALLOCATE(sw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(swf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dsw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dswf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(f,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chain_int
  ! ==================================================================
  SUBROUTINE inv3(a,c)
    ! ==--------------------------------------------------------------==
    ! ==            COMPUTES THE INVERSE OF A MATRIX 3*3              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a(3,3), c(3,3)

    INTEGER                                  :: i, iperm, ir, j, k, l
    REAL(real_8)                             :: b(3,3), den, s

! ==--------------------------------------------------------------==

    den = 0._real_8
    s   = 1._real_8
    i   = 1
    j   = 2
    k   = 3
1   CONTINUE
    DO iperm=1,3
       den = den + s*a(1,i)*a(2,j)*a(3,k)
       l   = i
       i   = j
       j   = k
       k   = l
    ENDDO

    i = 2
    j = 1
    k = 3
    s = - s
    IF (s.LT.0._real_8) THEN
       go to 1
    ENDIF

    IF (ABS(den).LT.1.e-20_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,*)'singular matrix'
    ENDIF

    i = 1
    j = 2
    k = 3

    DO ir=1,3
       b(ir,1) = (a(2,j)*a(3,k) - a(2,k)*a(3,j)) / den
       b(ir,2) = (a(3,j)*a(1,k) - a(3,k)*a(1,j)) / den
       b(ir,3) = (a(1,j)*a(2,k) - a(1,k)*a(2,j)) / den
       l = i
       i = j
       j = k
       k = l
    ENDDO

    DO l=1,3
       DO k=1,3
          c(k,l) = b(k,l)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE inv3
  ! ==================================================================
  SUBROUTINE r_to_s (r,s)
    ! ==--------------------------------------------------------------==
    ! == TRANSFORMS REAL VECTOR R IN SCALED VECTOR S                  ==
    ! == (THROUGH MATRIX HINV)                                        ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: r(*), s(*)

    REAL(real_8)                             :: h, hinv

    COMMON/metr_my/ h(3,3),hinv(3,3)

    INTEGER :: i,j
    ! ==--------------------------------------------------------------==
    DO i=1,3
       s(i) = 0._real_8
       DO j=1,3
          s(i) = s(i) + r(j)*hinv(j,i)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE r_to_s
  ! ==================================================================
  SUBROUTINE s_to_r (s,r)
    ! ==--------------------------------------------------------------==
    ! ==  TRANSFORMS SCALED VECTOR S IN REAL VECTOR R                 ==
    ! ==  (THROUGH MATRIX H)                                          ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: s(*), r(*)

    REAL(real_8)                             :: h, hinv

    COMMON/metr_my/ h(3,3),hinv(3,3)

    INTEGER :: i,j
    ! ==--------------------------------------------------------------==
    DO i=1,3
       r(i) = 0._real_8
       DO j=1,3
          r(i) = r(i) + s(j)*h(j,i)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s_to_r
  ! ==================================================================
  ! !$      subroutine  print(pos,hlist,nhlist,ch,n1,n2,MaxP)
  ! !$      
  ! !$      implicit none
  ! !$      integer n1,n2,ch(6+2,*)
  ! !$      real(8) :: pos(3,*)
  ! !$      integer hlist(n2,*),nhlist(n1,n2,*)
  ! !$      integer i,j,l,jj,kkk,kk,icont
  ! !$      real(8) :: conv,MaxP(*)
  ! !$!      real(8) :: avg_corr_time(2,*)
  ! !$
  ! !$      conv =.5291772
  ! !$
  ! !$      do i=1,6
  ! !$         icont =0
  ! !$         write(200+i,*)'REMARK ',MaxP(i)
  ! !$         do j=1,i+2
  ! !$            if (ch(j,i).eq.0)cycle
  ! !$            l=ch(j,i)
  ! !$            icont=icont+1
  ! !$            write(200+i,1)'ATOM  ',icont,'O   ','WAT',1,pos(1,l)*conv,pos(2,l)*conv,pos(3,l)*conv
  ! !$         enddo
  ! !$         do j=1,i+1
  ! !$            if (ch(j,i).eq.0)cycle
  ! !$            jj=ch(j+1,i)
  ! !$            l=ch(j,i)
  ! !$            do kk=1,hlist(jj,l)
  ! !$               kkk=nhlist(kk,jj,l)
  ! !$               icont=icont+1               
  ! !$               write(200+i,1)'ATOM  ',icont,'H   ','WAT',1, 
  ! !$     &              pos(1,kkk)*conv,pos(2,kkk)*conv,pos(3,kkk)*conv
  ! !$            enddo
  ! !$         enddo
  ! !$         write(200+i,*)'TER'
  ! !$      enddo
  ! !$
  ! !$      return
  ! !$1     FORMAT(a6,i5,1x,a5,aa3,2x,i4,4x,3f8.3)
  ! !$      end

END MODULE chain_dr_utils
