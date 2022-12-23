MODULE proylm_utils
  USE bessm_utils,                     ONLY: bessl
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: gk,&
                                             hg
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_ufo
  USE gvec,                            ONLY: gvec_com
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop7
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: proylm
  PUBLIC :: loadylm

CONTAINS

  SUBROUTINE proylm(c0,center,ny,nst,rmax,nmax)
    ! ==--------------------------------------------------------------==
    ! == ----  PROJECT WAVEFUNCTIONS ON REAL SPHERICAL HARMONICS -----==      
    ! ==--------------------------------------------------------------==
    ! == -----           PROGRAM CHANGES COEFFICIENTS             ---===
    REAL(real_8)                             :: center(3)
    INTEGER                                  :: ny, nst
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nst)
    REAL(real_8)                             :: rmax
    INTEGER                                  :: nmax

    CHARACTER(*), PARAMETER                  :: procedureN = 'proylm'
    REAL(real_8), PARAMETER :: sqofpi = 0.2820947917738781434_real_8 

    CHARACTER(len=40)                        :: dg
    COMPLEX(real_8), ALLOCATABLE             :: cc(:,:), cy(:,:), cyt(:,:), &
                                                egr(:)
    INTEGER                                  :: i, ierr, isub, j, j1, jmax, &
                                                jmax2, k, l, lmax, m, s
    LOGICAL                                  :: eqm, ferror, logm
    REAL(real_8)                             :: dp, dr, fpidsom, norm, sfc, &
                                                twpi, zero
    REAL(real_8), ALLOCATABLE                :: bsl(:,:), devl(:), gk2(:,:), &
                                                hg2(:), mesh(:), rlm(:), &
                                                SPREAD(:), SUM(:,:), w(:,:), &
                                                xtmat(:,:)

    CALL tiset('   PROWYLM',isub)     
    ALLOCATE(cy(ncpw%ngw,ny),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cyt(ny,ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! TRANSLATION VECTOR
    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)') 'ALL INPUT PARAMETERS :'
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,3(F10.5))') 'CENTER OF CLUSTER :',&
            center(1), center(2),center(3)
    ENDIF
    twpi=ATAN(1.0_real_8)*8.0_real_8
    ! LOAD ALL SPHERICAL HARMONICS 
    CALL loadylm(cyt,ny,ncpw%ngw)
    ! TRANSPOSE CYT
    DO i=1,ny
       DO j=1,ncpw%ngw
          cy(j,i)=cyt(i,j)
       ENDDO
    ENDDO

    DEALLOCATE(cyt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! OTHER ALLOCATION 
    ALLOCATE(egr(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hg2(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gk2(3,ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! INITIALIZE CARTESIAN COORDINATES OF RECIPROCAL LATTICE VECTORS
    DO i=1,ncpw%ngw
       gk2(1,i)=(gk(1,i)*gvec_com%b1(1)+gk(2,i)*gvec_com%b2(1)+gk(3,i)*gvec_com%b3(1))*twpi/parm%alat
       gk2(2,i)=(gk(1,i)*gvec_com%b1(2)+gk(2,i)*gvec_com%b2(2)+gk(3,i)*gvec_com%b3(2))*twpi/parm%alat
       gk2(3,i)=(gk(1,i)*gvec_com%b1(3)+gk(2,i)*gvec_com%b2(3)+gk(3,i)*gvec_com%b3(3))*twpi/parm%alat
       hg2(i)=SQRT(gk2(1,i)*gk2(1,i)+gk2(2,i)*gk2(2,i)+&
            gk2(3,i)*gk2(3,i))
    ENDDO
    ! FIND EXP(I*<G,R>)  
    DO i=1,ncpw%ngw
       dp=center(1)*gk2(1,i)+center(2)*gk2(2,i)+center(3)*gk2(3,i)
       egr(i)=CMPLX(COS(dp),SIN(dp),kind=real_8)
    ENDDO
    DEALLOCATE(gk2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xtmat(ny,nst),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cc(ncpw%ngw,ny),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! FIND OVERLAP MATRIX
    CALL dcopy(2*ncpw%ngw*ny,cy,1,cc,1)

    DO i=1,ncpw%ngw
       DO j=1,ny
          ! CC(I,J)=CC(I,J)*CMPLX(real(EGR(I)),-aimag(EGR(I)))
       ENDDO
    ENDDO

    DO i=1,nst
       DO j=1,ncpw%ngw
          c0(j,i)=c0(j,i)*egr(j)
       ENDDO
    ENDDO

    DO i=1,ny
       sfc=dotp(ncpw%ngw,cc(:,i),cc(:,i))
       CALL mp_sum(sfc,parai%allgrp)
       sfc=1._real_8/SQRT(sfc)
       CALL zdscal(ncpw%ngw,sfc,cc(1,i),1)
    ENDDO

    CALL ovlap2(ncpw%ngw,ny,nst,xtmat,cc,c0,.TRUE.)
    CALL mp_sum(xtmat,ny*nst,parai%allgrp)

    ! WRITE OVERLAP MATRIX
    IF (paral%parent)THEN
       CALL ints40(dg,ny,s)
       IF (paral%io_parent)&
            CALL fileopen(41,'YLMCOEF',fo_def,ferror)
       IF (paral%io_parent)&
            WRITE(41,'('//dg(41-s:40)//'(F18.14))')&
            ((xtmat(l,k),l=1,ny),k=1,nst)
       IF (paral%io_parent)&
            CALL fileclose(41)
    ENDIF
    DEALLOCATE(cc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (prop7%radial)THEN
       ALLOCATE(SPREAD(nst),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (prop7%spread_ham)THEN
          ! READ WANNIER INPUT         
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  CALL fileopen(60,'SPREAD',fo_def,ferror)
             DO i=1,nst
                IF (paral%io_parent)&
                     READ(60,*) j, dr, norm, zero, SPREAD(i)
             ENDDO
             IF (paral%io_parent)&
                  CALL fileclose(60)
             IF (paral%io_parent)&
                  WRITE(6,"(1X,A,F10.5,A,F10.5)")  &
                  'PROJECT STATES ONLY WITH SPREAD FROM ', &
                  prop7%spr_min,' TO ',prop7%spr_max
          ENDIF
          ! SPREAD
          CALL mp_bcast(spread,SIZE(spread),parai%source,parai%allgrp)
       ELSE
          prop7%spr_min=-1.0_real_8
          prop7%spr_max=10.0_real_8
          DO i=1,nst
             SPREAD(i)=0.0_real_8
          ENDDO
       ENDIF
       IF (paral%parent)THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' RECORD OF ALL RADIAL PARTS'
          IF (paral%io_parent)&
               CALL fileopen(41,'RADIAL',fo_def+fo_ufo,ferror)
       ENDIF
       ! FIND RADIAL PART OF L-TH COMPONENT  AND WEIGHTS.
       lmax=INT(SQRT(ABS(1.0_real_8*ny-1.0)))
       ! NUMBER OF SEGMENTS SHOULD BE EVEN
       IF (MOD(nmax,2).EQ.1) nmax=nmax+1
       ALLOCATE(rlm((nmax+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(mesh(nmax+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(SUM(lmax+1,nst),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(devl(lmax+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(w(ny,nst),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(bsl(ncpw%ngw,nmax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(devl)!,lmax+1)
       fpidsom=fpi/SQRT(parm%omega)
       logm=.FALSE.
       eqm=.TRUE.
       mesh(0)=0.0_real_8
       ! LOGARITHMIC 
       IF (logm)THEN
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'LOGARITHMIC GRID IS USED'
          ENDIF
          mesh(1)=0.01
          sfc=LOG(rmax/mesh(1))/REAL(nmax-1,kind=real_8)
          sfc=EXP(sfc)
          DO i=2,nmax
             mesh(i)=mesh(i-1)*sfc
          ENDDO
       ENDIF
       ! EQUIDISTANT
       IF (eqm)THEN
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'EQUIDISTANT GRID IS USED'
          ENDIF
          dr=rmax/REAL(nmax,kind=real_8)
          DO i=1,nmax
             mesh(i)=dr*i
          ENDDO
       ENDIF

       IF (paral%parent)THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I8)') 'NUMBER OF GRIDS =',nmax
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,F10.5)') 'RADIUS OF INTEGRATION =',rmax
       ENDIF

       CALL zeroing(sum)!,(lmax+1)*nst)
       IF (paral%parent)THEN
          ! to radial file all parameters
          IF (paral%io_parent)&
               WRITE(41) logm,eqm,nmax,rmax,&
               prop7%spread_ham, prop7%spr_min, prop7%spr_max,&
               lmax,nst
       ENDIF
       DO l=0,lmax
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'CURRENT L =', l
          ENDIF
          CALL loadbsl(ncpw%ngw,nmax,l,hg2,mesh,bsl)
          DO m=1,2*l+1
             DO i=1,nst
                IF ((SPREAD(i).GT.prop7%spr_min).AND.&
                     (SPREAD(i).LT.prop7%spr_max)) THEN
                   j1=l*l+m
                   ! POINT FROM DR, RMAX
                   IF (l.EQ.0)THEN
                      DO j=1,nmax
                         rlm(j)=0.0_real_8
                         DO k=1,ncpw%ngw
                            IF (hg(k).LT.1.e-6_real_8)THEN
                               rlm(j)=rlm(j)+&
                                    bsl(k,j)*&
                                    REAL(c0(k,i))*&
                                    sqofpi
                            ELSE
                               rlm(j)=rlm(j)+&
                                    2*bsl(k,j)*&
                                    REAL(c0(k,i))*&
                                    sqofpi
                            ENDIF
                         ENDDO
                         rlm(j)=rlm(j)*fpidsom
                      ENDDO
                   ELSE
                      IF (MOD(l,2).EQ.0)THEN
                         DO j=1,nmax
                            rlm(j)=0.0_real_8
                            DO k=1,ncpw%ngw
                               IF (hg(k).GT.1.e-6_real_8)THEN
                                  rlm(j)=rlm(j)+&
                                       bsl(k,j)*&
                                       REAL(c0(k,i))*REAL(cy(k,j1))
                               ENDIF
                            ENDDO
                            rlm(j)=2*rlm(j)*fpidsom
                         ENDDO
                      ELSE
                         DO j=1,nmax
                            rlm(j)=0.0_real_8
                            DO k=1,ncpw%ngw
                               IF (hg(k).GT.1.e-6_real_8)THEN
                                  rlm(j)=rlm(j)+&
                                       bsl(k,j)*&
                                       AIMAG(c0(k,i))*&
                                       AIMAG(cy(k,j1))
                               ENDIF
                            ENDDO
                            rlm(j)=2*rlm(j)*fpidsom
                         ENDDO
                      ENDIF
                   ENDIF
                   CALL mp_sum(rlm(:),nmax+1,parai%allgrp)
                   IF (paral%parent) THEN
                      IF (paral%io_parent)&
                           WRITE(41)(rlm(j),j=0,nmax)
                   ENDIF
                   IF (eqm)THEN
                      CALL intsimp(nmax,dr,rlm(0),w(j1,i))
                   ELSE
                      CALL wint(nmax,mesh,rlm(0),w(j1,i))
                   ENDIF
                   SUM(l,i)=SUM(l,i)+w(j1,i)
                   ! I=1,NST
                ENDIF
             ENDDO
             ! M=1,2*L+1
          ENDDO
          ! L=0,L
       ENDDO

       IF ((paral%parent).AND.paral%io_parent)&
            CALL fileclose(41)
       ! WRITE ALL WEIGHTS
       IF (paral%parent)THEN
          IF (ny.EQ.1)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'ORBITAL AND MAIN HARMONIC(L)'
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A)') 'ORBITAL AND 2 FIRST MAIN HARMONICS(L)'
          ENDIF
          DO i=1,nst
             IF ((SPREAD(i).GT.prop7%spr_min).AND.&
                  (SPREAD(i).LT.prop7%spr_max)) THEN
                CALL maxarr(lmax+1,jmax,jmax2,SUM(0,i))
                IF (ny.EQ.1)THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,I8,A)') i,'-TH ORBITAL'
                ELSE
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,I8,A,2(/,1X,T10,I6,F10.5))')&
                        i,'-TH ORBITAL',&
                        jmax-1, SUM(jmax-1,i), jmax2-1, SUM(jmax2-1,i)
                ENDIF
                DO l=0,lmax
                   DO m=1,2*l+1
                      j1=l*l+m
                      IF (paral%io_parent)&
                           WRITE(6,*) 'W(',j1,')=',w(j1,i)
                   ENDDO
                   IF (paral%io_parent)&
                        WRITE(6,*) 'SUM(',l,')=',SUM(l,i)
                ENDDO
             ENDIF
          ENDDO

          DO i=1,nst
             IF ((SPREAD(i).GT.prop7%spr_min).AND.&
                  (SPREAD(i).LT.prop7%spr_max)) THEN
                CALL maxarr(lmax+1,jmax,jmax2,SUM(0,i))
                IF (ny.EQ.1)THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,I8,A)') i,'-TH ORBITAL'
                ELSE
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,I8,A,2(/,1X,T10,I6,F10.5))')&
                        i,'-TH ORBITAL',&
                        jmax-1, SUM(jmax-1,i), jmax2-1, SUM(jmax2-1,i)
                ENDIF
             ENDIF
          ENDDO
       ENDIF

       IF (paral%parent)THEN
          ! CALCULATE AND WRITE DEVIATIONS 
          DO l=0,lmax
             devl(l)=0.0_real_8
          ENDDO
          jmax=0
          DO i=1,nst
             IF ((SPREAD(i).GT.prop7%spr_min).AND.&
                  (SPREAD(i).LT.prop7%spr_max)) THEN
                jmax=jmax+1
                DO l=0,lmax
                   DO j=l,lmax
                      devl(j)=devl(j)+SUM(l,i)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

          DO j=0,lmax
             devl(j)=1.0_real_8-devl(j)/REAL(jmax,kind=real_8)
          ENDDO

          DO j=0,lmax
             IF (paral%io_parent)&
                  WRITE(6,*)  'E(',j,') = ',devl(j)
          ENDDO

          IF (paral%io_parent)&
               CALL fileopen(41,'WEIGHTS',fo_def+fo_ufo,ferror)
          IF (paral%io_parent)&
               WRITE(41) ((w(k,i),i=1,nst),k=1,ny)
          IF (paral%io_parent)&
               CALL fileclose(41)

          IF (paral%io_parent)&
               CALL fileopen(41,'SUMWEIGHTS',fo_def+fo_ufo,ferror)
          IF (paral%io_parent)&
               WRITE(41) ((SUM(k,i),i=1,nst),k=0,lmax)
          IF (paral%io_parent)&
               CALL fileclose(41)
       ENDIF
    ENDIF
    CALL tihalt('   PROWYLM',isub)
  END SUBROUTINE proylm

  SUBROUTINE loadylm(cyt,ny,nngw)
    INTEGER                                  :: ny, nngw
    COMPLEX(real_8)                          :: cyt(ny,nngw)

    REAL(real_8), PARAMETER :: sqofpi = 0.2820947917738781434_real_8 

    INTEGER                                  :: i, isub, lmax
    REAL(real_8)                             :: norm, x, y, z

    CALL tiset('   LOADYLM',isub)  
    lmax=INT(SQRT(ABS(1.0_real_8*ny-1.0)))
    DO i=1,nngw
       IF (hg(i).LT.1.e-6_real_8)THEN
          CALL zeroing(cyt(:,i))!,ny)
          cyt(1,i)=CMPLX(sqofpi,0.0_real_8,kind=real_8)
       ELSE
          norm=MAX(SQRT(hg(i)),1.0e-6_real_8)
          x=gk(1,i)/norm
          y=gk(2,i)/norm
          z=gk(3,i)/norm
          CALL loadall(x,y,z,cyt(1,i),ny,lmax)
       ENDIF
    ENDDO
    CALL tihalt('   LOADYLM',isub)
  END SUBROUTINE loadylm

  SUBROUTINE loadall(x,y,z,ylm,n,lmax)
    REAL(real_8)                             :: x, y, z
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: ylm(n)
    INTEGER                                  :: lmax

    CHARACTER(*), PARAMETER                  :: procedureN = 'loadall'
    REAL(real_8), PARAMETER :: sq2 = 1.41421356237309504880_real_8 , &
      sqofpi = 0.2820947917738781434_real_8 

    COMPLEX(real_8)                          :: c, il
    INTEGER                                  :: i, ierr, j, k, l, lmm, m
    REAL(real_8)                             :: iks, r, sq, twl1
    REAL(real_8), ALLOCATABLE                :: coph(:), lay1(:), lay2(:), &
                                                lay3(:), siph(:)

    ALLOCATE(coph(lmax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(siph(lmax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lay1(lmax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lay2(lmax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lay3(lmax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    coph(0)=1.0_real_8
    siph(0)=0.0_real_8
    c=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    DO i=1,lmax
       c=c*CMPLX(x,y,kind=real_8)
       coph(i)=REAL(c)
       siph(i)=AIMAG(c)
    ENDDO

    lay1(0)=1.0_real_8
    lay2(0)=z
    lay2(1)=1.0_real_8

    sq=sqofpi*sq2
    DO l=0,lmax
       j=l*l+l+1
       twl1=SQRT(2.0_real_8*REAL(l,kind=real_8)+1.0_real_8)
       ylm(j)=sqofpi*twl1
       DO m=1,l
          j=j+1
          ylm(j)=sq*twl1
       ENDDO
    ENDDO

    ! LEGENDRE POLYNOMIALS
    ylm(1)=ylm(1)*lay1(0)
    IF (lmax.GT.0)THEN
       ylm(2)=ylm(2)*lay2(1)
       ylm(3)=ylm(3)*lay2(0)
       ylm(4)=ylm(4)*lay2(1)
    ENDIF

    j=2*2
    DO l=2,lmax
       twl1=2.0_real_8*REAL(l,kind=real_8)-1.0_real_8
       DO i=0,l-2
          lay3(i)=(twl1*z*lay2(i)-REAL(l+i-1,kind=real_8)*lay1(i))/REAL(l-i,kind=real_8)
       ENDDO
       lay3(l)=1.0_real_8
       lay3(l-1)=z*twl1

       DO i=0,l
          lay1(i)=lay2(i)
          lay2(i)=lay3(i)
       ENDDO

       j=j+l
       DO i=0,l
          j=j+1
          ylm(j)=lay3(i)*ylm(j)
       ENDDO
    ENDDO


    IF (lmax.GT.0)THEN
       ylm(2)=ylm(2)/sq2
       ylm(4)=ylm(4)/sq2
    ENDIF

    DO l=2,lmax
       j=l*l+l
       DO m=0,l
          j=j+1
          r=1.0_real_8
          lmm=l-m
          DO i=1,m
             iks=(2.0_real_8*REAL(i,kind=real_8)-1.0_real_8)
             r=r*iks*iks/((iks+REAL(lmm,kind=real_8)+1.0_real_8)*(iks+REAL(lmm,kind=real_8)))
          ENDDO
          ylm(j)=ylm(j)*SQRT(r)
       ENDDO
    ENDDO

    il=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    DO l=1,lmax
       il=il*CMPLX(0.0_real_8,1.0_real_8,kind=real_8)
       j=l*l
       k=2*l+1+j
       DO m=l,1,-1
          j=j+1
          ylm(j)=ylm(k)*coph(m)*il
          k=k-1
       ENDDO
       j=j+1
       ylm(j)=ylm(j)*il
       DO m=1,l
          j=j+1
          ylm(j)=ylm(j)*siph(m)*il
       ENDDO
    ENDDO
    DEALLOCATE(coph,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(siph,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lay1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lay2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lay3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE loadall

  SUBROUTINE loadbsl(ngw,nmax,l,hg2,mesh,bsl)
    INTEGER                                  :: ngw, nmax, l
    REAL(real_8)                             :: hg2(ngw), mesh(0:nmax), &
                                                bsl(ngw,1:nmax)

    INTEGER                                  :: i, isub, j

! VARIABLES

    CALL tiset('   LOADBSL',isub)

    DO i=1,nmax
       DO j=1,ngw
          bsl(j,i)=bessl(l,hg2(j)*mesh(i))
       ENDDO
    ENDDO

    CALL tihalt('   LOADBSL',isub)
  END SUBROUTINE loadbsl

  SUBROUTINE ints40(digit,n,s)
    CHARACTER(len=40)                        :: digit
    INTEGER                                  :: n, s

    INTEGER                                  :: i, r, y

! VARIABLES

    y=n
    IF (y.LE.0)THEN
       digit(40:40)='0'
       s=1
    ELSE
       i=41
       DO WHILE(y.GT.0)
          i=i-1
          r=MOD(y,10)
          IF (r.EQ.0)THEN
             digit(i:i)='0'
          ELSE IF (r.EQ.1)THEN
             digit(i:i)='1'
          ELSE IF (r.EQ.2)THEN
             digit(i:i)='2'
          ELSE IF (r.EQ.3)THEN
             digit(i:i)='3'
          ELSE IF (r.EQ.4)THEN
             digit(i:i)='4'
          ELSE IF (r.EQ.5)THEN
             digit(i:i)='5'
          ELSE IF (r.EQ.6)THEN
             digit(i:i)='6'
          ELSE IF (r.EQ.7)THEN
             digit(i:i)='7'
          ELSE IF (r.EQ.8)THEN
             digit(i:i)='8'
          ELSE IF (r.EQ.9)THEN
             digit(i:i)='9'
          ENDIF
          y=y/10
       ENDDO
       s=41-i
    ENDIF
  END SUBROUTINE ints40

  SUBROUTINE  wint(n,mesh,rad,w)
    INTEGER                                  :: n
    REAL(real_8)                             :: mesh(n), rad(0:n), w

    INTEGER                                  :: i
    REAL(real_8)                             :: ri, ri1, sum

! RESULT IN W
! VARIABLES

    sum=0.0_real_8
    ri=mesh(1)
    sum=sum+ri*ri*rad(1)*rad(1)*0.5
    ! ZERO POINT
    DO i=1,n-1
       ri=mesh(i)
       ri1=mesh(i+1)
       sum=sum+(rad(i)*rad(i)*ri*ri+rad(i+1)*rad(i+1)*ri1*ri1)*0.5*&
            (ri1-ri)
    ENDDO
    w=sum
  END SUBROUTINE  wint

  SUBROUTINE intsimp(n,dr,rad,w)
    INTEGER                                  :: n
    REAL(real_8)                             :: dr, rad(0:n), w

    INTEGER                                  :: i
    REAL(real_8)                             :: c23, c43, r1, r2, sum

! RESULT IN W
! VARIABLES     

    c43=4.0_real_8/3.0_real_8
    c23=2.0_real_8/3.0_real_8
    r1=dr
    r2=r1+dr
    sum=0.0_real_8
    DO i=1,n-3,2
       sum=sum+c43*r1*r1*rad(i)*rad(i)+c23*r2*r2*rad(i+1)*rad(i+1)
       r1=r2+dr
       r2=r1+dr
    ENDDO

    sum=sum+c43*r1*r1*rad(n-1)*rad(n-1)+&
         1.0_real_8/3.0_real_8*r2*r2*rad(n)*rad(n)

    sum=sum*dr
    w=sum
  END SUBROUTINE intsimp

  SUBROUTINE maxarr(ny,jmax,jmax2,w)
    INTEGER                                  :: ny, jmax, jmax2
    REAL(real_8)                             :: w(ny)

    INTEGER                                  :: i
    REAL(real_8)                             :: mx, mx2

! VARIABLES 

    IF (ny.NE.1)THEN
       IF (w(1)>w(2))THEN
          mx=w(1)
          jmax=1
          mx2=w(2)
          jmax2=2
       ELSE
          mx=w(2)
          jmax=2
          mx2=w(1)
          jmax2=1
       ENDIF
    ELSE
       mx=w(1)
       jmax=1
    ENDIF

    DO i=3,ny
       IF (w(i)>mx)THEN
          mx2=mx
          jmax2=jmax
          mx=w(i)
          jmax=i
       ELSE IF (w(i)>mx2)THEN
          jmax2=i
          mx2=w(i)
       ENDIF
    ENDDO

  END  SUBROUTINE maxarr

END MODULE proylm_utils
