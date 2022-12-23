MODULE symtrz_utils
  USE cnst,                            ONLY: pi
  USE cppt,                            ONLY: gl,&
                                             hg,&
                                             igl,&
                                             indz,&
                                             inyh,&
                                             isptr,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: epsg
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE symm,                            ONLY: &
       igextinc, indshel, irt, irtvec, isymu, isyshel, numshel, symmi, symmr, &
       symmt, tvec
  USE symm4,                           ONLY: isishel,&
                                             ninshel
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: symvec
  PUBLIC :: symmat
  PUBLIC :: symstress
  PUBLIC :: symrho
  PUBLIC :: give_scr_symvec
  PUBLIC :: give_scr_symmat
  PUBLIC :: give_scr_symrho

CONTAINS

  ! ==================================================================
  SUBROUTINE symvec(force)
    ! ==--------------------------------------------------------------==
    ! == Symmetrization of ionic forces                               ==
    ! == Do not depend on fractional translation                      ==
    ! == Depend on translation vectors (if not primitive cell)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: force(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'symvec'
    REAL(real_8), PARAMETER                  :: epsilon = 1.e-15_real_8

    INTEGER                                  :: ia, iat, ib, ierr, ir, is, &
                                                iv, l
    REAL(real_8)                             :: gg
    REAL(real_8), ALLOCATABLE                :: xau(:,:), yau(:,:)

    IF ((symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1.AND.symmi%ntvec.LE.1)) RETURN
    ALLOCATE(xau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(yau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! Transform to crystal axis basis
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL dgemv('T',3,3,parm%alat,metr_com%htm1(1,1),3,force(1,ia,is),1,0.0_real_8,&
               xau(1,iat),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL zeroing(yau)!,3*ions1%nat)
    ! Symmetry from rotations
    IF (symmi%ntvec.EQ.0) THEN
       DO ia=1,ions1%nat
          DO ir=1,symmi%nrot
             ib=irt(ir,ia)
             DO l=1,3
                yau(l,ia)=yau(l,ia)+symmr%xtable(l,1,ir)*xau(1,ib)+&
                     symmr%xtable(l,2,ir)*xau(2,ib)+&
                     symmr%xtable(l,3,ir)*xau(3,ib)
             ENDDO
          ENDDO
       ENDDO
       CALL dscal(3*ions1%nat,1._real_8/REAL(symmi%nrot,kind=real_8),yau(1,1),1)
    ELSE
       DO ia=1,ions1%nat
          DO ir=1,symmi%nrot
             DO iv=1,symmi%ntvec
                ib=irtvec(irt(ir,ia),iv)
                DO l=1,3
                   yau(l,ib)=yau(l,ib)+symmr%xtable(l,1,ir)*xau(1,ia)+&
                        symmr%xtable(l,2,ir)*xau(2,ia)+&
                        symmr%xtable(l,3,ir)*xau(3,ia)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       CALL dscal(3*ions1%nat,1._real_8/REAL(symmi%nrot*symmi%ntvec,kind=real_8),yau(1,1),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Transform back to cartesian coordinates
    iat=0
    gg=1._real_8/parm%alat
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL dgemv('T',3,3,gg,metr_com%ht(1,1),3,yau(1,iat),1,0._real_8,&
               force(1,ia,is),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Set to 0 if < EPSILON
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO l=1,3
             IF (ABS(force(l,ia,is)).LT.epsilon) force(l,ia,is)=0._real_8
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(xau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(yau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE symvec
  ! ==================================================================
  SUBROUTINE symmat(tns,isu)
    ! ==--------------------------------------------------------------==
    ! == Symmetrization of matrix A(3,3) from each couple of atoms    ==
    ! == Do not depend on fractional translations                     ==
    ! == Depend on translation vectors (if not primitive cell)        ==
    ! == TNS(3,NAT,3,NAT) sets of matrices                            ==
    ! == ISU Give factor                                              ==
    ! ==     ISU=0 FAC=1                                              ==
    ! ==     ISU=1 FAC=ISYMU(IA) (n. of equiv. atoms for the line IA) ==
    ! ==     ISU=2 FAC=ISYMU(IB) (n. of equiv. atoms for the row IB)  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tns(3*ions1%nat,3*ions1%nat)
    INTEGER                                  :: isu

    CHARACTER(*), PARAMETER                  :: procedureN = 'symmat'

    INTEGER                                  :: i, ia, ia1, ian, ib, ib1, &
                                                ibn, ierr, ir, iv, j, k, l
    REAL(real_8)                             :: fac, gg, xmat(3,3)
    REAL(real_8), ALLOCATABLE                :: xau(:,:), yau(:,:)

! Variables
! real(8) :: XAU(3,*),YAU(3,*)
! ==--------------------------------------------------------------==

    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1.AND.symmi%ntvec.LE.1) ) RETURN
    ALLOCATE(xau(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(yau(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! Transform to crystal axis basis
    DO ia=1,ions1%nat
       ia1=3*(ia-1)+1
       DO ib=1,ions1%nat
          ib1=3*(ib-1)+1
          CALL dgemm('T','N',3,3,3,parm%alat,metr_com%htm1(1,1),3,tns(ia1,ib1),3*ions1%nat,&
               0._real_8,xmat(1,1),3)
          CALL dgemm('N','N',3,3,3,parm%alat,xmat(1,1),3,metr_com%htm1(1,1),3,0._real_8,&
               xau(ia1,ib1),3*ions1%nat)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL zeroing(yau)!,3*ions1%nat*3*ions1%nat)
    ! Symmetry from rotations
    DO ia=1,ions1%nat
       ia1=3*(ia-1)
       DO ib=1,ions1%nat
          ib1=3*(ib-1)
          IF (isu.EQ.0) THEN
             fac=1._real_8
          ELSEIF (isu.EQ.1) THEN
             fac=REAL(isymu(ia),kind=real_8)
          ELSEIF (isu.EQ.2) THEN
             fac=REAL(isymu(ib),kind=real_8)
          ELSE
             CALL stopgm('SYMMAT',' ISU ',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (fac.NE.0._real_8) THEN
             DO ir=1,symmi%nrot
                DO iv=1,symmi%ntvec
                   ian=3*( irtvec(irt(ir,ia),iv) - 1 )
                   ibn=3*( irtvec(irt(ir,ib),iv) - 1 )
                   CALL zeroing(xmat)!,9)
                   ! Matrix multiplication FAC XTABLE.XTABLE.XAU
                   DO i=1,3
                      DO j=1,3
                         DO k=1,3
                            DO l=1,3
                               xmat(i,j)=xmat(i,j)+&
                                    symmr%xtable(i,k,ir)*symmr%xtable(j,l,ir)*&
                                    xau(ia1+k,ib1+l)*fac
                            ENDDO! L
                         ENDDO! K
                      ENDDO! J
                   ENDDO ! I
                   ! Put in YAU
                   DO i=1,3
                      DO j=1,3
                         yau(ian+i,ibn+j)=yau(ian+i,ibn+j)+xmat(i,j)
                      ENDDO! J
                   ENDDO ! I
                ENDDO     ! IV
             ENDDO         ! IR
          ENDIF             ! IF(FAC.NE.0._real_8)
       ENDDO                 ! IB
    ENDDO                  ! IA
    CALL dscal(3*ions1%nat*3*ions1%nat,1._real_8/REAL(symmi%nrot*symmi%ntvec,kind=real_8),yau(1,1),1)
    ! ==--------------------------------------------------------------==
    ! Transform back to cartesian coordinates
    gg=1._real_8/parm%alat
    DO ia=1,ions1%nat
       ia1=3*(ia-1)+1
       DO ib=1,ions1%nat
          ib1=3*(ib-1)+1
          CALL dgemm('T','N',3,3,3,gg,metr_com%ht(1,1),3,&
               yau(ia1,ib1),3*ions1%nat,0._real_8,xmat(1,1),3)
          CALL dgemm('N','N',3,3,3,gg,xmat(1,1),3,&
               metr_com%ht(1,1),3,0._real_8,tns(ia1,ib1),3*ions1%nat)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(xau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(yau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE symmat
  ! ==================================================================
  SUBROUTINE symstress(stress)
    ! ==--------------------------------------------------------------==
    ! == Symmetrization of stress tensor                              ==
    ! == Do not depend on fractional translations                     ==
    ! == Do not depend on translation vectors                         ==
    ! ==--------------------------------------------------------------==
    ! == STRESS(1:6) is given with Voigt notation                     ==
    ! == ALPHA(1:6)  first index for Voigt notation                   ==
    ! == BETA(1:6)   second index                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: stress(3,3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-15_real_8 

    INTEGER                                  :: i, ir, j
    REAL(real_8)                             :: astr(3,3), bstr(3,3), sum

! ==--------------------------------------------------------------==

    IF ((symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1)) RETURN
    ! Transform to crystal axis basis
    CALL dgemm('T','N',3,3,3,1._real_8,metr_com%htm1,3,stress,3,0._real_8,astr,3)
    CALL dgemm('N','T',3,3,3,1._real_8,astr,3,metr_com%ht,3,0._real_8,stress,3)
    CALL zeroing(bstr)!,3*3)
    DO ir=1,symmi%nrot
       CALL dgemm('N','N',3,3,3,1._real_8,symmr%xtable(1,1,ir),3,stress,3,&
            0._real_8,astr,3)
       CALL dgemm('N','N',3,3,3,1._real_8,astr,3,symmr%xtable(1,1,symmi%inve(ir)),3,&
            1._real_8,bstr,3)
    ENDDO
    CALL dscal(3*3,1._real_8/REAL(symmi%nrot,kind=real_8),bstr,1)
    ! Transform back to cartesian coordinates      
    CALL dgemm('T','N',3,3,3,1._real_8,metr_com%ht,3,bstr,3,0._real_8,astr,3)
    CALL dgemm('N','T',3,3,3,1._real_8,astr,3,metr_com%htm1,3,0._real_8,stress,3)
    ! Symmetric matrix and set to 0 if < EPSILON
    DO i=1,3
       DO j=i,3
          IF (ABS(stress(j,i)).LT.epsilon) THEN
             stress(j,i)=0._real_8
             stress(i,j)=0._real_8
          ELSE
             sum=0.5_real_8*(stress(j,i)+stress(i,j))
             stress(j,i)=sum
             stress(i,j)=sum
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE symstress
  ! ==================================================================
  SUBROUTINE symrho(rho,psi)
    ! ==--------------------------------------------------------------==
    ! == SYMMETRIZE THE CHARGE DENSITY IN FOURIER SPACE               ==
    ! == ONLY THE NHG COMPONENTS GIVEN BY NZH AND INDZ                ==
    ! ==--------------------------------------------------------------==
    ! == XTABLE is in crystal basis and FTABLE is in reciprocal basis ==
    ! == For each operation IR, R^-1 gives the symmetry.              ==
    ! == As R^-1 is in the group, we use XTABLE(IR) and FTAU(INVE(IR))==
    ! ==--------------------------------------------------------------==
    ! == More accurate than symmetrization in real space              ==
    ! == Faster (20-60 times faster)                                  ==
    ! ==--------------------------------------------------------------==
    ! == IF DUAL IS LESS THAN 3, CALL SYMRHO_ALL WHICH SYMMETRIZES    ==
    ! == ALL THE DENSITY ARRAY (MORE EXPENSIVE)                       ==
    ! == ONLY IN SERIAL VERSION                                       ==
    ! ==--------------------------------------------------------------==
    ! ==  BTEST(I,POS) .TRUE. if bit POS is 1 in I integer            ==
    ! ==  Chech if BTEST works with INTEGER*8                         ==
    ! ==  otherwise use -D__NOINT8                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho(:)
    COMPLEX(real_8)                          :: psi(:)

    COMPLEX(real_8)                          :: crhom, crhop
    COMPLEX(real_8), ALLOCATABLE             :: array(:,:), psip(:)
    INTEGER, ALLOCATABLE                     :: irotp(:,:)
    REAL(real_8), ALLOCATABLE                :: scr(:)

! Variables

#if defined(__NOINT8)
    LOGICAL :: btest4
    ! INTEGER :: NINSHEL(2,NHG),ISISHEL(2,NHG) ! defined in symm4.mod
    INTEGER :: imask(2),isign(2),ir
#else
    ! INTEGER*8 :: NINSHEL(NHG),ISISHEL(NHG) ! defined in symm4.mod
    INTEGER(int_8) :: imask,isign,ir
#endif
    INTEGER :: isub,ig,ish,indig,ir2,nh1,nh2,nh3,&
         irftau,npack,ipack,jnpack,nzh0,indz0,&
         ishl,ishig,ishend,lsymrho,ierr
    REAL(real_8) :: vnrot,arg,tmpi,tpi,rgl1,rgl2,rgl3
    LOGICAL :: tconjig
    CHARACTER(len=30) :: tag
    CHARACTER(*),PARAMETER::procedureN='SYMRHO'
    ! ==--------------------------------------------------------------==
    IF ((symmi%indpg.EQ.0.OR.symmi%indpg.EQ.100) .OR.&
         (symmi%nrot.EQ.1.AND.symmi%ntvec.LE.1)) RETURN
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    CALL tiset('    SYMRHO',isub)
    CALL setfftn(0)
    ! We need to initialize some arrays.
    IF (numshel.EQ.0) THEN
       CALL buildshell
    ENDIF
    CALL give_scr_symrho(lsymrho,tag)
    ! ==--------------------------------------------------------------==
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    tpi=2._real_8*pi
    tmpi=-2._real_8*pi
    IF (symmt%tsymmorphic) THEN
       IF (symmi%inversion.EQ.0) THEN
          ALLOCATE(psip(lsymrho/2),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ELSE
       ipack = lsymrho/(3*symmi%nrot)
       ipack = -ipack
       CALL mp_max(ipack,parai%allgrp)
       ipack = -ipack
       ALLOCATE(irotp(symmi%nrot, lsymrho/(3*symmi%nrot)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(array(symmi%nrot, lsymrho/(3*symmi%nrot)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(psi)!,maxfft)
    CALL dcopy(fpar%nnr1,rho(1),1,psi(1),2)
    ! Transform to reciprocal space
    CALL fwfftn(psi,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! Symmetry from translation vectors -> extinction
    IF (igextinc.NE.0) THEN
       ! There is IGEXTINC extinction (RHO(IG)=0)
       !CDIR NODEP
       DO ig=1,ncpw%nhg
          IF (indshel(ig).EQ.0) THEN
             psi(nzh (ig))=CMPLX(0._real_8,0._real_8,kind=real_8)
             psi(indz(ig))=CMPLX(0._real_8,0._real_8,kind=real_8)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! NROT = 1 (no rotation symmetry except the identical one).
    IF ((symmi%nrot.GT.1).AND.symmt%tsymmorphic) THEN
       ! Translate density in order to have 0 as  centre of symmetry
       IF (.NOT.symmt%torigin) THEN
          !CDIR NODEP
          DO ig=1,ncpw%nhg
             IF (indshel(ig).NE.0) THEN
                arg=tpi*(&
                     REAL(inyh(1,ig)-nh1,kind=real_8)*symmr%origin(1)+&
                     REAL(inyh(2,ig)-nh2,kind=real_8)*symmr%origin(2)+&
                     REAL(inyh(3,ig)-nh3,kind=real_8)*symmr%origin(3) )
                crhop=psi(nzh(ig))*CMPLX(COS(arg),SIN(arg),kind=real_8)
                psi(nzh (ig))=crhop
                psi(indz(ig))=CONJG(crhop)
             ENDIF
          ENDDO
       ENDIF
       IF (symmi%inversion.EQ.0) THEN
          vnrot=1._real_8/REAL(symmi%nrot,kind=real_8)
          CALL zeroing(psip(1:numshel))!,numshel)
          ! Loop over G-vectors. For each equivalent shell,
          ! the electronic density is the same.
          DO ig=1,ncpw%nhg
             indig=indshel(ig)
             IF (indig.NE.0) THEN
                psip(indig)=psip(indig)&
                     + isyshel(1,ig)*psi(nzh (ig))&
                     + isyshel(2,ig)*psi(indz(ig))
             ENDIF
          ENDDO
          CALL mp_sum(psip,numshel,parai%allgrp)
          !CDIR NODEP
          DO ig=1,ncpw%nhg
             indig=indshel(ig)
             IF (indig.NE.0) THEN
                IF (isyshel(1,ig).GT.0) THEN
                   IF (isyshel(2,ig).EQ.0) THEN
                      crhop=vnrot*psip(indig)
                      psi(nzh (ig))=crhop
                      psi(indz(ig))=CONJG(crhop)
                   ELSE
                      arg=vnrot*REAL(psip(indig))
                      psi(nzh (ig))=CMPLX(arg,0._real_8,kind=real_8)
                      psi(indz(ig))=CMPLX(arg,0._real_8,kind=real_8)
                   ENDIF
                ELSE
                   crhop=vnrot*psip(indig)
                   psi(nzh (ig))=CONJG(crhop)
                   psi(indz(ig))=crhop
                ENDIF
             ENDIF
          ENDDO
       ELSE
          ! Special case for inversion
          vnrot=1._real_8/REAL(symmi%nrot/2,kind=real_8)
          ALLOCATE(scr(numshel),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL zeroing(scr)!,numshel)
          ! Loop over G-vectors. For each equivalent shell,
          ! the electronic density is the same.
          DO ig=1,ncpw%nhg
             indig=indshel(ig)
             IF (indig.NE.0) THEN
                scr(indig)=scr(indig)+isyshel(1,ig)*REAL(psi(nzh(ig)))
             ENDIF
          ENDDO
          CALL mp_sum(scr,numshel,parai%allgrp)
          !CDIR NODEP
          DO ig=1,ncpw%nhg
             indig=indshel(ig)
             IF (indig.NE.0) THEN
                arg=vnrot*scr(indig)
                psi(nzh (ig))=CMPLX(arg,0._real_8,kind=real_8)
                psi(indz(ig))=CMPLX(arg,0._real_8,kind=real_8)
             ENDIF
          ENDDO
          ! All PSI array is real
          CALL imagset0(psi,maxfft)
          DEALLOCATE(scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
       ! Translate density to -ORIGIN
       IF (.NOT.symmt%torigin) THEN
          !CDIR NODEP
          DO ig=1,ncpw%nhg
             IF (indshel(ig).NE.0) THEN
                arg=tmpi*(&
                     REAL(inyh(1,ig)-nh1,kind=real_8)*symmr%origin(1)+&
                     REAL(inyh(2,ig)-nh2,kind=real_8)*symmr%origin(2)+&
                     REAL(inyh(3,ig)-nh3,kind=real_8)*symmr%origin(3) )
                crhop=psi(nzh(ig))*CMPLX(COS(arg),SIN(arg),kind=real_8)
                psi(nzh (ig))=crhop
                psi(indz(ig))=CONJG(crhop)
             ENDIF
          ENDDO
       ENDIF
    ELSEIF (symmi%nrot.GT.1) THEN
       ! ==------------------------------------------------------------==
       ! == Non-symmorphic case.                                       ==
       ! == Maybe you could have some trouble with the DUAL option     ==
       ! == the only solution is to increase the DUAL parameter
       ! ==------------------------------------------------------------==
       vnrot=1._real_8/REAL(symmi%nrot,kind=real_8)
       ! Loop over equivalent shell.
       ! For each rotation, we collect the corresponding components.
       ! Then we can calculate for each G-vector.
       DO ishl=1,numshel,ipack
          npack=MIN((numshel-ishl+1),ipack)
          CALL zeroing(irotp)!,symmi%nrot*npack)
          CALL zeroing(array)!,symmi%nrot*npack)
          ishend=ishl+(npack-1)
          DO ig=1,ncpw%nhg
             ishig=indshel(ig)
             IF (  (ishig.GE.ishl).AND.(ishig.LE.ishend) ) THEN
                crhop=psi(nzh (ig))
                crhom=psi(indz(ig))
#if defined(__NOINT8)
                imask(1)=ninshel(1,ig)
                imask(2)=ninshel(2,ig)
                isign(1)=isishel(1,ig)
                isign(2)=isishel(2,ig)
                jnpack=(ishig-ishl)+1
                IF ((jnpack.LT.1.OR.jnpack.GT.npack).AND.paral%io_parent)&
                     WRITE(6,*) jnpack
                DO ir=1,symmi%nrot
                   IF (btest4(imask(1),ir)) THEN
                      IF (btest4(isign(1),ir)) THEN
                         array(ir,jnpack)=crhop
                      ELSE
                         array(ir,jnpack)=crhom
                      ENDIF
                      irotp(ir,jnpack)=ig
                   ENDIF
                ENDDO
#else
                imask=ninshel(ig)
                isign=isishel(ig)
                jnpack=(ishig-ishl)+1
                IF ((jnpack.LT.1.OR.jnpack.GT.npack).AND.paral%io_parent)&
                     WRITE(6,*) jnpack
                DO ir=1,symmi%nrot
                   IF (BTEST(imask,ir)) THEN
                      IF (BTEST(isign,ir)) THEN
                         array(ir,jnpack)=crhop
                      ELSE
                         array(ir,jnpack)=crhom
                      ENDIF
                      irotp(ir,jnpack)=ig
                   ENDIF
                ENDDO
#endif
             ENDIF
          ENDDO
100       CONTINUE
          ! Glosum for the NROT components for NPACK equivalent shells.
          CALL mp_sum(array,symmi%nrot*npack,parai%allgrp)
          DO ish=1,npack
             DO ir=1,symmi%nrot
                IF (irotp(ir,ish).NE.0) THEN
                   ig=irotp(ir,ish)
                   crhop=array(ir,ish)
#if defined(__NOINT8)
                   tconjig=btest4(isishel(1,ig),ir)
#else
                   tconjig=BTEST(isishel(ig),ir)
#endif
                   IF (tconjig) THEN
                      rgl1=inyh(1,ig)-nh1
                      rgl2=inyh(2,ig)-nh2
                      rgl3=inyh(3,ig)-nh3
                   ELSE
                      rgl1=nh1-inyh(1,ig)
                      rgl2=nh2-inyh(2,ig)
                      rgl3=nh3-inyh(3,ig)
                   ENDIF
                   ! Give the number of the inverse of the rotation
                   ! which gives RGL1,RGL2,RGL3 from IG
                   DO ir2=1,symmi%nrot
                      IF (ir2.NE.ir) THEN
                         IF (ir.EQ.1) THEN
                            irftau=symmi%inve(ir2)
                         ELSE
                            irftau=symmi%multab(symmi%inve(ir2),ir)
                         ENDIF
                         IF (symmt%tsymmor(irftau)) THEN
                            crhop=crhop+array(ir2,ish)
                         ELSE
                            arg = tmpi *&
                                 ( rgl1*symmr%ftau(1,irftau)&
                                 + rgl2*symmr%ftau(2,irftau)&
                                 + rgl3*symmr%ftau(3,irftau) )
                            crhop=crhop+array(ir2,ish)&
                                 *CMPLX(COS(arg),SIN(arg),kind=real_8)
                         ENDIF
                         IF (ig.EQ.irotp(ir2,ish)) irotp(ir2,ish)=0
                      ENDIF
                   ENDDO
                   nzh0 =nzh (ig)
                   indz0=indz(ig)
                   crhop=vnrot*crhop
                   IF (tconjig) THEN
                      psi(nzh0 )=crhop
                      psi(indz0)=CONJG(crhop)
                   ELSE
                      psi(indz0)=crhop
                      psi(nzh0 )=CONJG(crhop)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO                 ! DO ISHL=1,NUMSHEL,IPACK
       ! Special case for inversion. We do work in ALL ARRAY.
       IF (symmi%inversion.NE.0) THEN
          IF (symmt%tsymmor(symmi%inversion)) CALL imagset0(psi,maxfft)
       ENDIF
    ENDIF                     ! IF(TSYMMORPHIC)
    ! ==--------------------------------------------------------------==
    ! Transform back to real space
    CALL  invfftn(psi,.FALSE.,parai%allgrp)
    CALL dcopy(fpar%nnr1,psi(1),2,rho(1),1)
    CALL tihalt('    SYMRHO',isub)
    ! ==--------------------------------------------------------------==
    IF (symmt%tsymmorphic) THEN
       IF (symmi%inversion == 0) THEN
          DEALLOCATE(psip,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ELSE
       DEALLOCATE(irotp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(array,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE symrho
  ! ==================================================================
  SUBROUTINE buildshell
    ! ==--------------------------------------------------------------==
    ! == Use bit manipulations functions                              ==
    ! ==  IBSET(I,POS) sets the bit POS to 1 in I integer             ==
    ! ==  BTEST(I,POS) .TRUE. if bit POS is 1 in I integer            ==
    ! ==  Chech if BTEST and IBTEST work with INTEGER*8               ==
    ! ==  otherwise use -D__NOINT8                                    ==
    ! ==--------------------------------------------------------------==
    ! == Symetrisation due to translation vectors (some extinction)   ==
    ! == RHO(R+TVEC) =RHO(R)                                          ==
    ! == i.e. RHO(G) = RHO(G) Exp(i G.TVEC)                           ==
    ! == if Exp(i G.TVEC) = 1 INDSHEL gives equivalent shell          ==
    ! == if Exp(i G.TVEC)\= 1 -> RHO(G) = 0                           ==
    ! == By symMetry true for all G-vectors of the equivalent shell   ==
    ! ==                                                              ==
    ! == There are also other extinction associated with rotations    ==
    ! == and glides.                                                  ==
    ! == If there is some extinction IGEXTINC is set to their number  ==
    ! == and INDSHEL(IG)=0                                            ==
    ! ==--------------------------------------------------------------==
    ! == The symetry are determined with a precision of DELTASYM      ==
    ! == in the direct lattice vector basis                           ==
    ! == Furthermore the ordering has a precision of EPSG*NHGS        ==
    ! == So for translation COS(X) = 1._real_8 with precision of DELTASYM  ==
    ! == and for the equivalent shell determination                   ==
    ! == DELTAG=MAX(DELTASYM*SQRT(HG(IG)),EPSG*NHGS)                  ==
    ! ==--------------------------------------------------------------==
    ! == IFIRST=0 first call to BUILDSHELL                            ==
    ! == IFIRST=1 not the first call and before TSYMMORPHIC=.TRUE.    ==
    ! == IFIRST=2 not the first call and before TSYMMORPHIC=.FALSE.   ==
    ! ==--------------------------------------------------------------==
    ! == TO DO:                                                       ==
    ! ==   Parallel version. Instead of broadcasting one G vector per ==
    ! ==   one, broadcast many G vectors to all processors.           ==
    ! ==   In this case BUILDSHELL is more parallelizable.            ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER, PARAMETER                       :: indnoshel = 2**(8*4-2)

    REAL(real_8)                             :: arg, deltag, r(4), tpi

! NINSHEL is INTEGER*8

#if defined(__NOINT8)
    LOGICAL :: btest4
    ! INTEGER :: NINSHEL(2,NHG),ISISHEL(2,NHG) ! defined in symm4.mod
    INTEGER :: isum(2),ir
#else
    ! INTEGER*8 :: NINSHEL(NHG),ISISHEL(NHG) ! defined in symm4.mod
    INTEGER(int_8) :: isum,ir
#endif
    INTEGER :: nh1,nh2,nh3,msglen,ip,igg,ig,ifac,igbeg,igend,invir,&
         igl1,igl2,igl3,ig1,ig2,ig3,isub,ish,nshl,kshell,&
         lbuildshell,ierr
    INTEGER, SAVE :: ifirst = 0

    REAL(real_8), ALLOCATABLE :: scr(:)
    CHARACTER(*),PARAMETER::procedureN='buildshell'
    ! ==--------------------------------------------------------------==
    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1.AND.symmi%ntvec.LE.1) ) RETURN
    ! 63 because we use integer*8
    IF (symmi%nrot.GT.63) CALL stopgm('BUILDSHELL',&
         'THE NUMBER OF ROTATIONS IS TOO BIG',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tpi=2._real_8*pi
    CALL tiset('BUILDSHELL',isub)
    IF (ifirst.EQ.0) THEN
       ! First call.
       ALLOCATE(indshel(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (symmt%tsymmorphic) THEN
          IF (symmi%nrot.GT.1) THEN
             ALLOCATE(isyshel(2,ncpw%nhg),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          ifirst=1
       ELSE
          IF (symmi%nrot.GT.1) THEN
             ! NINSHEL and ISISHEL are INTEGER*8
#if defined(__NOINT8)
             ALLOCATE(ninshel(2,nhg),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(isishel(2,nhg),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
#else
             ALLOCATE(ninshel(2*ncpw%nhg),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(isishel(2*ncpw%nhg),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
#endif
          ENDIF
          ifirst=2
       ENDIF
    ELSEIF (ifirst.EQ.1.AND.(.NOT.symmt%tsymmorphic)) THEN
       ! Before TSYMMORPHIC, now .NOT.TSYMMORPHIC
       IF (symmi%nrot.GT.1) THEN
          ! Deallocation of ISYSHEL
          DEALLOCATE(isyshel,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ! Allocations of NINSHEL and ISISHEL
#if defined(__NOINT8)
          ALLOCATE(ninshel(2,nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(isishel(2,nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
#else
          ALLOCATE(ninshel(2*ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(isishel(2*ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
#endif
       ENDIF
       ifirst=2
    ELSEIF (ifirst.EQ.2.AND.symmt%tsymmorphic) THEN
       ! Before .NOT.TSYMMORPHIC, now TSYMMORPHIC
       IF (symmi%nrot.GT.1) THEN
          ! Deallocations of NINSHEL and ISISHEL
          DEALLOCATE(ninshel,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(isishel,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ! Allocation of ISYSHEL
          ALLOCATE(isyshel(2,ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ifirst=1
    ENDIF
    ! ==--------------------------------------------------------------==
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    igextinc=0
    ! Special case for extinction only (no rotations)
    IF (symmi%nrot.EQ.1) THEN
       ! NTVEC/=0 -> only extinction
       DO ig=1,ncpw%nhg
          r(1)=REAL(inyh(1,ig)-nh1,kind=real_8)
          r(2)=REAL(inyh(2,ig)-nh2,kind=real_8)
          r(3)=REAL(inyh(3,ig)-nh3,kind=real_8)
          ! Check if there is extinction or not
          DO ir=1,symmi%ntvec
             arg=COS( tpi*(r(1)*tvec(1,ir)+&
                  r(2)*tvec(2,ir)+r(3)*tvec(3,ir)))
             IF (ABS(arg-1._real_8).GT.symmr%deltasym) THEN
                ! Extinction
                igextinc=igextinc+1
                indshel(ig)=0
                GOTO 50
             ENDIF
          ENDDO
          ! No extinction
          indshel(ig)=1
50        CONTINUE
       ENDDO
       ! Work is done
       numshel=1
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    DO ig=1,ncpw%nhg
       indshel(ig)=indnoshel
    ENDDO
    IF (symmt%tsymmorphic) THEN
       CALL zeroing(isyshel)!,2*nhg)
    ELSE
       CALL zeroing(ninshel)!,nhg*2)
       CALL zeroing(isishel)!,nhg*2)
    ENDIF
    numshel=0
    msglen= 4 * 8
    ! Each processor sends
    DO ip=1,parai%nproc
       IF (parap%pgroup(ip).EQ.parai%me) THEN
          DO ig=1,ncpw%nhg
             IF (indshel(ig).EQ.indnoshel) THEN
                ! Another equivalent shell
                ig1=inyh(1,ig)-nh1
                ig2=inyh(2,ig)-nh2
                ig3=inyh(3,ig)-nh3
                r(1)=REAL(ig1,kind=real_8)
                r(2)=REAL(ig2,kind=real_8)
                r(3)=REAL(ig3,kind=real_8)
                r(4)=hg(ig)
                IF (symmi%ntvec.NE.0) THEN
                   ! Check if there is extinction or not
                   DO ir=1,symmi%ntvec
                      arg=COS( tpi*(r(1)*tvec(1,ir)+&
                           r(2)*tvec(2,ir)+r(3)*tvec(3,ir)))
                      IF (ABS(arg-1._real_8).GT.symmr%deltasym) THEN
                         ! Extinction
                         igextinc=igextinc+1
                         indshel(ig)=0
                         GOTO 200
                      ENDIF
                   ENDDO
                ENDIF
                IF (.NOT.symmt%tsymmorphic) THEN
                   ! For not symmorphic group, there are
                   ! rotation associated with glide:
                   ! Some extinction could be exist
                ENDIF
                CALL mp_bcast(r,SIZE(r),parap%pgroup(ip),parai%allgrp)
                numshel=numshel+1
                indshel(ig)=numshel
                IF (symmt%tsymmorphic) THEN
                   isyshel(1,ig)=isyshel(1,ig)+1
                ELSE
                   ir=1
#if defined(__NOINT8)
                   CALL ibset4(ninshel(1,ig),ir)
                   CALL ibset4(isishel(1,ig),ir)
#else
                   ninshel(ig)=IBSET(ninshel(ig),ir)
                   isishel(ig)=IBSET(isishel(ig),ir)
#endif
                ENDIF
                igbeg=ig
                nshl=igl(ig)
                ! direct DELTASYM*ALAT
                ! reciprocal 2PI/ALAT*HI(IG)*DELTASYM
                deltag=MAX(r(4)*symmr%deltasym,epsg*REAL(spar%nhgs,kind=real_8))
                DO ish=nshl,ncpw%nhgl
                   IF (ABS(gl(ish)-r(4)).LE.deltag) THEN
                      igend=isptr(ish+1)-1
                   ELSE
                      igend=isptr(ish+1)-1
                      GOTO 500
                   ENDIF
                ENDDO
500             CONTINUE
                ! This processor finds equivalent G-vectors.
                DO ir=2,symmi%nrot
                   ! Rotations conserves norms
                   igl1 = NINT(&
                        symmr%ftable(1,1,ir) * r(1)&
                        + symmr%ftable(1,2,ir) * r(2)&
                        + symmr%ftable(1,3,ir) * r(3) )
                   igl2 = NINT(&
                        symmr%ftable(2,1,ir) * r(1)&
                        + symmr%ftable(2,2,ir) * r(2)&
                        + symmr%ftable(2,3,ir) * r(3) )
                   igl3 = NINT(&
                        symmr%ftable(3,1,ir) * r(1)&
                        + symmr%ftable(3,2,ir) * r(2)&
                        + symmr%ftable(3,3,ir) * r(3) )
                   IF (igl1.LT.0.OR.&
                        (igl1.EQ.0.AND.igl2.LT.0).OR.&
                        (igl1.EQ.0.AND.igl2.EQ.0.AND.igl3.LT.0)) THEN
                      ifac=2
                      igl1=-igl1+nh1
                      igl2=-igl2+nh2
                      igl3=-igl3+nh3
                   ELSE
                      ifac=1
                      igl1=igl1+nh1
                      igl2=igl2+nh2
                      igl3=igl3+nh3
                   ENDIF
                   IF (igl1.EQ.ig1.AND.igl2.EQ.ig2.AND.igl3.EQ.ig3) THEN
                      ! FTABLE.G=G vector
                      IF (symmt%tsymmorphic) THEN
                         isyshel(ifac,ig)=isyshel(ifac,ig)+1
                      ELSE
#if defined(__NOINT8)
                         CALL ibset4(ninshel(1,ig),ir)
                         IF (ifac.EQ.1) THEN
                            CALL ibset4(isishel(1,ig),ir)
                         ENDIF
#else
                         ninshel(ig)=IBSET(ninshel(ig),ir)
                         IF (ifac.EQ.1) THEN
                            isishel(ig)=IBSET(isishel(ig),ir)
                         ENDIF
#endif
                      ENDIF
                   ELSE
                      ! Looking for the image of IG.
                      DO igg=igbeg,igend
                         IF (  igl1.EQ.inyh(1,igg).AND.&
                              igl2.EQ.inyh(2,igg).AND.&
                              igl3.EQ.inyh(3,igg)) THEN
                            ! Find an equivalent G-vector, continue
                            indshel(igg)=numshel
                            IF (symmt%tsymmorphic) THEN
                               isyshel(ifac,igg)=isyshel(ifac,igg)+1
                            ELSE
#if defined(__NOINT8)   
                               CALL ibset4(ninshel(1,igg),ir)
                               IF (ifac.EQ.1) THEN
                                  CALL ibset4(isishel(1,igg),ir)
                               ENDIF
#else                   
                               ninshel(igg)=IBSET(ninshel(igg),ir)
                               IF (ifac.EQ.1)&
                                    isishel(igg)=IBSET(isishel(igg),ir)
#endif                  
                            ENDIF
                            GOTO 1000
                         ENDIF
                      ENDDO
                      IF (parai%nproc.EQ.1) THEN
                         IF (paral%io_parent)&
                              WRITE(6,*) 'IG=',ig,' IGL=',igl1,igl2,igl3
                         CALL stopgm('BUILDSHELL','G-VECTOR NOT FOUND',& 
                              __LINE__,__FILE__)
                      ENDIF
                      ! Another rotation.
                   ENDIF
1000               CONTINUE
                ENDDO
                ! No more equivalent G-vectors, now the next.
             ENDIF
200          CONTINUE
          ENDDO
          ! It is over for this processor, information for the other ones
          r(1)=-1._real_8
          CALL mp_bcast(r,SIZE(r),parap%pgroup(ip),parai%allgrp)
       ELSE
          ! The other processors listen to IP processor
          r(1)=1._real_8
          nshl=1
          DO WHILE(r(1).GE.0._real_8)
             CALL mp_bcast(r,SIZE(r),parap%pgroup(ip),parai%allgrp)
             IF (r(1).LT.0._real_8) GOTO 3000
             ! Another equivalent shell
             numshel=numshel+1
             ! Find first IG possible.
             kshell=0
             ! Due to reordering we cannot trust the shells (see rggen).
             ! direct DELTASYM*ALAT
             ! reciprocal 2PI/ALAT*HG(IG)*DELTASYM
             deltag=MAX(r(4)*symmr%deltasym,epsg*REAL(spar%nhgs,kind=real_8))
             DO ish=1,ncpw%nhgl
                IF (ABS(gl(ish)-r(4)).LE.deltag) THEN
                   IF (kshell.EQ.0) THEN
                      nshl=ish
                      igbeg=isptr(ish)
                      igend=isptr(ish+1)-1
                      kshell=1
                   ELSE
                      igend=isptr(ish+1)-1
                   ENDIF
                ELSEIF (kshell.EQ.1) THEN
                   igend=isptr(ish+1)-1
                ELSEIF (gl(ish).GT.(r(4)+deltag)) THEN
                   GOTO 1500
                ENDIF
             ENDDO
1500         CONTINUE
             ! No equivalent shell for this processor.
             IF (kshell.EQ.0) GOTO 3000
             ig1=NINT(r(1))
             ig2=NINT(r(2))
             ig3=NINT(r(3))
             DO ir=2,symmi%nrot
                ! Rotations conserves norms
                igl1 = NINT(&
                     symmr%ftable(1,1,ir) * r(1)&
                     + symmr%ftable(1,2,ir) * r(2)&
                     + symmr%ftable(1,3,ir) * r(3) )
                igl2 = NINT(&
                     symmr%ftable(2,1,ir) * r(1)&
                     + symmr%ftable(2,2,ir) * r(2)&
                     + symmr%ftable(2,3,ir) * r(3) )
                igl3 = NINT(&
                     symmr%ftable(3,1,ir) * r(1)&
                     + symmr%ftable(3,2,ir) * r(2)&
                     + symmr%ftable(3,3,ir) * r(3) )
                IF (igl1.LT.0.OR.&
                     (igl1.EQ.0.AND.igl2.LT.0).OR.&
                     (igl1.EQ.0.AND.igl2.EQ.0.AND.igl3.LT.0)) THEN
                   ifac=2
                   igl1=-igl1+nh1
                   igl2=-igl2+nh2
                   igl3=-igl3+nh3
                ELSE
                   ifac=1
                   igl1=igl1+nh1
                   igl2=igl2+nh2
                   igl3=igl3+nh3
                ENDIF
                IF (igl1.EQ.ig1.AND.igl2.EQ.ig2.AND.igl3.EQ.ig3) THEN
                   ! IGG=IG
                ELSE
                   DO igg=igbeg,igend
                      IF (  igl1.EQ.inyh(1,igg).AND.&
                           igl2.EQ.inyh(2,igg).AND.&
                           igl3.EQ.inyh(3,igg)) THEN
                         ! Find one.
                         indshel(igg)=numshel
                         IF (symmt%tsymmorphic) THEN
                            isyshel(ifac,igg)=isyshel(ifac,igg)+1
                         ELSE
#if defined(__NOINT8)
                            CALL ibset4(ninshel(1,igg),ir)
                            IF (ifac.EQ.1) THEN
                               CALL ibset4(isishel(1,igg),ir)
                            ENDIF
#else
                            ninshel(igg)=IBSET(ninshel(igg),ir)
                            IF (ifac.EQ.1)&
                                 isishel(igg)=IBSET(isishel(igg),ir)
#endif
                         ENDIF
                         GOTO 2000
                      ENDIF
                   ENDDO
                ENDIF
2000            CONTINUE
             ENDDO         ! NROT
3000         CONTINUE
          ENDDO             ! end of while
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Check if all G-vector has a shell
    DO ig=1,ncpw%nhg
       IF (indshel(ig).EQ.indnoshel) THEN
          IF ((parai%nproc.NE.1).AND.paral%io_parent)&
               WRITE(6,'("[PROC=",I4,"]")',advance="no") parai%me
          IF (paral%io_parent)&
               WRITE(6,*) 'G-VECTOR=',ig,' HAS NO SHELL NUMBER'
          CALL stopgm('BUILDSHELL','NO SHELL NUMBER',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    ALLOCATE(scr(numshel),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(scr)!,numshel)
    IF (symmt%tsymmorphic) THEN
       IF (symmi%inversion.EQ.0) THEN
          DO ig=1,ncpw%nhg
             ish=indshel(ig)
             IF (ish.NE.0) THEN
                scr(ish)=scr(ish)+REAL(isyshel(1,ig)+isyshel(2,ig),kind=real_8)
             ENDIF
          ENDDO
       ELSE
          IF (geq0) isyshel(1,1)=isyshel(1,1)/2
          DO ig=1,ncpw%nhg
             ish=indshel(ig)
             IF (ish.NE.0) THEN
                scr(ish)=scr(ish)+REAL(2*isyshel(1,ig),kind=real_8)
             ENDIF
          ENDDO
       ENDIF
    ELSE
#if defined(__NOINT8)
       DO ig=1,nhg
          ish=indshel(ig)
          IF (ish.NE.0) THEN
             isum(1)=ninshel(1,ig)
             isum(2)=ninshel(2,ig)
             ifac=0
             DO ir=1,symmi%nrot
                IF (btest4(isum(1),ir)) ifac=ifac+1
             ENDDO
             scr(ish)=scr(ish)+REAL(ifac,kind=real_8)
          ENDIF
       ENDDO
#else
       DO ig=1,ncpw%nhg
          ish=indshel(ig)
          IF (ish.NE.0) THEN
             isum=ninshel(ig)
             ifac=0
             DO ir=1,symmi%nrot
                IF (BTEST(isum,ir)) ifac=ifac+1
             ENDDO
             scr(ish)=scr(ish)+REAL(ifac,kind=real_8)
          ENDIF
       ENDDO
#endif
    ENDIF
    CALL mp_sum(scr,numshel,parai%allgrp)
    IF (paral%io_parent) THEN
       DO ish=1,numshel
          IF (ABS(scr(ish)-symmi%nrot).GT.0.0001_real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ish,scr(ish),symmi%nrot
             CALL stopgm('BUILDSHELL','BAD NUMBER OF ROTATIONS',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDIF
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('BUILDSHELL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE buildshell
  ! ==================================================================
  SUBROUTINE give_scr_symvec(lsymvec,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsymvec
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1) ) THEN
       lsymvec=0
    ELSE
       lsymvec=6*ions1%nat+10
       tag='6*NAT'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_symvec
  ! ==================================================================
  SUBROUTINE give_scr_symmat(lsymmat,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsymmat
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1) ) THEN
       lsymmat=0
    ELSE
       lsymmat=2*9*ions1%nat*ions1%nat+10
       tag='2*9*NAT*NAT'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_symmat
  ! ==================================================================
  SUBROUTINE give_scr_buildshell(lbuildshell,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lbuildshell
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1) ) THEN
       lbuildshell=0
    ELSE
       IF (numshel.EQ.0) THEN
          lbuildshell=3*spar%nhgs/symmi%nrot+10
          tag='3*NHGS/NROT'
       ELSE
          lbuildshell=numshel+10
          tag='NUMSHEL'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_buildshell
  ! ==================================================================
  SUBROUTINE imagset0(psi,n)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: n

    INTEGER                                  :: i

    DO i=1,n
       psi(i)=CMPLX(REAL(psi(i),KIND=real_8),0._real_8,KIND=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE imagset0
  ! ==================================================================
  SUBROUTINE give_scr_symrho(lsymrho,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsymrho
    CHARACTER(len=30)                        :: tag

    INTEGER, PARAMETER                       :: maxpack = 20000 

    INTEGER                                  :: lbuildshell, lpack

! ==--------------------------------------------------------------==

    IF ( (symmi%indpg.EQ.0).OR.(symmi%nrot.EQ.1) ) THEN
       lsymrho=0
    ELSE
       IF (symmt%tsymmorphic) THEN
          IF (numshel.EQ.0) THEN
             lsymrho=spar%nhgs/MAX(1,MIN(symmi%nrot/2,symmi%nrot/parai%nproc))
          ELSE
             lsymrho=numshel
          ENDIF
          IF (symmi%inversion.EQ.0) lsymrho=2*lsymrho! Use complex
       ELSE                  ! PSIP, IROTP (if not symmorphic)
          IF (numshel.EQ.0) THEN
             lpack=maxpack
          ELSE
             lpack=MIN(maxpack,numshel*(symmi%nrot*2+symmi%nrot))
          ENDIF
          lsymrho=MAX(lpack,symmi%nrot*2+symmi%nrot)
       ENDIF
       CALL give_scr_buildshell(lbuildshell,tag)
       lsymrho=MAX(lbuildshell,lsymrho)+10
       tag='2*NUMSHEL'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_symrho
#if defined(__NOINT8)
  FUNCTION btest4(i,pos)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i(2), pos
    LOGICAL                                  :: btest4

! ==--------------------------------------------------------------==

    IF (pos.GT.31) THEN
       btest4=BTEST(i(2),pos-32)
    ELSE
       btest4=BTEST(i(1),pos)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION btest4
  ! ==================================================================
  SUBROUTINE ibset4(i,pos)
    INTEGER                                  :: i(2), pos

! ==--------------------------------------------------------------==

    IF (pos.GT.31) THEN
       i(2)=IBSET(i(2),pos-32)
    ELSE
       i(1)=IBSET(i(1),pos)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ibset4
  ! ==================================================================
#endif

END MODULE symtrz_utils
