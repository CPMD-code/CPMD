MODULE efield_utils
  USE ddip,                            ONLY: afield,&
                                             aifield,&
                                             lenbk,&
                                             ngwmax
  USE ddipo_utils,                     ONLY: setdip
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE opeigr_p_utils,                  ONLY: opeigr_p
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: extfield

CONTAINS

  ! ==================================================================
  SUBROUTINE extfield(tau0,fion,c0,c2,nstate)
    ! ==--------------------------------------------------------------==
    ! == APPLIES FINITE HOMOGENEOUS ELECTRIC FIELD USING BERRY PHASE  ==
    ! ==                                                              ==
    ! == 20/03/2009 Rewritten for general use (JAEA-SPring8/IPCMS)    ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    COMPLEX(real_8)                          :: c0(*), c2(*)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'extfield'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8), &
                                                zzero = (0._real_8,0._real_8) 
    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8

    COMPLEX(real_8)                          :: dd, ephase(3)
    COMPLEX(real_8), ALLOCATABLE             :: ddmat_e(:,:), work(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cw(:,:), cwork(:,:,:), &
                                                ffield(:,:), sc0(:)
    INTEGER                                  :: i, ia, ierr, ig, info, is, &
                                                is1, isub, k, n1, nmcol, &
                                                nop1, nop2, nop3, nxx
    INTEGER, ALLOCATABLE                     :: ipiv(:)
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: a, b, c, d, fac, p(3), &
                                                phase1, phase2, phase3, pin, &
                                                ufg, zvtot

! ==--------------------------------------------------------------==

    CALL tiset('    EFIELD',isub)
    ! ..initialization on first call to subroutine
    IF (ifirst.EQ.0) THEN
       n1=0
       DO i=0,parai%nproc-1
          n1=MAX(n1,parap%sparm(3,i))
       ENDDO
       ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ngwmax=n1
       nmcol=parai%nproc*ngwmax
       ALLOCATE(mapcol(nmcol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setdip(mapful,mapcol)
       ALLOCATE(ffield(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cwork(ncpw%ngw,nstate,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       lenbk=nxxfun(nstate)
       nxx=MAX(2*lenbk*parai%nproc,2*ncpw%ngw*nstate+8)
       ALLOCATE(sc0(nxx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cw(ncpw%ngw,nxx/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ..electric field in lattice frame
       IF (cntl%tlsd) THEN
          fac=1._real_8
       ELSE
          fac=2._real_8
       ENDIF
       pin=fac*0.5_real_8/ACOS(-1._real_8)
       afield(1)=pin*(aifield(1)*parm%a1(1)+aifield(2)*parm%a1(2)&
            +AIFIELD(3)*parm%a1(3))
       afield(2)=pin*(aifield(1)*parm%a2(1)+aifield(2)*parm%a2(2)&
            +AIFIELD(3)*parm%a2(3))
       afield(3)=pin*(aifield(1)*parm%a3(1)+aifield(2)*parm%a3(2)&
            +AIFIELD(3)*parm%a3(3))
       ifirst=1
    ENDIF
    ALLOCATE(ddmat_e(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipiv(20*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(nstate*20),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(p)!,3)
    zvtot=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          p(1)=p(1)+ions0%zv(is)*tau0(1,ia,is)
          p(2)=p(2)+ions0%zv(is)*tau0(2,ia,is)
          p(3)=p(3)+ions0%zv(is)*tau0(3,ia,is)
          ! ..corrections to ionic forces
          IF (paral%parent) THEN
             fion(1,ia,is)=fion(1,ia,is)+ions0%zv(is)*aifield(1)
             fion(2,ia,is)=fion(2,ia,is)+ions0%zv(is)*aifield(2)
             fion(3,ia,is)=fion(3,ia,is)+ions0%zv(is)*aifield(3)
          ENDIF
       ENDDO
       zvtot=zvtot+ions0%zv(is)*ions0%na(is)
    ENDDO
    DO i=1,3
       p(i)=p(i)/zvtot
    ENDDO
    ! ..phase factors due to center of charge motion
    ufg=parm%tpiba*nstate
    phase1=ufg*(p(1)*gvec_com%b1(1)+p(2)*gvec_com%b1(2)+p(3)*gvec_com%b1(3))
    phase2=ufg*(p(1)*gvec_com%b2(1)+p(2)*gvec_com%b2(2)+p(3)*gvec_com%b2(3))
    phase3=ufg*(p(1)*gvec_com%b3(1)+p(2)*gvec_com%b3(2)+p(3)*gvec_com%b3(3))
    ephase(1)=CMPLX(COS(phase1),SIN(phase1),kind=real_8)
    ephase(2)=CMPLX(COS(phase2),SIN(phase2),kind=real_8)
    ephase(3)=CMPLX(COS(phase3),SIN(phase3),kind=real_8)
    ! ..electronic contribution
    ener_com%eefield=0._real_8
    DO k=1,3
       IF (ABS(afield(k)).GT.epsilon) THEN
          nop1=k
          nop2=0
          nop3=0
          CALL opeigr_p(c0,cw,sc0,nstate,mapful,mapcol,ddmat_e,&
               nop1,nop2,nop3,dd,cwork)
          CALL zgetrf(nstate,nstate,ddmat_e,nstate,ipiv,info)
          IF (info.EQ.0) THEN
             dd=zone
             !$omp parallel do private(I) reduction(*:DD)
             DO i=1,nstate
                dd=dd*ddmat_e(i,i)
             ENDDO
             CALL zgetri(nstate,ddmat_e,nstate,ipiv,work,nstate,info)
          ELSE
             CALL stopgm('EFIELD','ERROR IN MATRIX INVERSION',& 
                  __LINE__,__FILE__)
          ENDIF
          CALL zgemm('N','N',ncpw%ngw,nstate,nstate,zone,cwork(1,1,1)&
               ,ncpw%ngw,ddmat_e(1,1),nstate,zzero,cw(1,1),ncpw%ngw)
          CALL zcopy(ncpw%ngw*nstate,cw(1,1),1,cwork(1,1,1),1)
          CALL zgemm('N','N',ncpw%ngw,nstate,nstate,zone,cwork(1,1,2)&
               ,ncpw%ngw,ddmat_e(1,1),nstate,zzero,cw(1,1),ncpw%ngw)
          CALL zcopy(ncpw%ngw*nstate,cw(1,1),1,cwork(1,1,2),1)
          !$omp parallel private(I,IG,A,B,C,D) default(shared)
          IF (geq0) THEN
             is1=2
             !$omp do
             DO i=1,nstate
                b=AIMAG(cwork(1,i,1))
                ffield(1,i)=CMPLX(b,0._real_8,kind=real_8)
             ENDDO
             !$omp enddo
          ELSE
             is1=1
          ENDIF
          !$omp do
          DO i=1,nstate
             DO ig=is1,ncpw%ngw
                a=REAL(cwork(ig,i,1))
                b=AIMAG(cwork(ig,i,1))
                c=REAL(cwork(ig,i,2))
                d=AIMAG(cwork(ig,i,2))
                ffield(ig,i)=0.5_real_8*CMPLX(b+d,c-a,kind=real_8)
             ENDDO
          ENDDO
          !$omp enddo
          !$omp end parallel
          dd=dd*ephase(k)
          p(k)=ATAN(AIMAG(dd)/REAL(dd))
          ! ..electric enthalpy
          ener_com%eefield=ener_com%eefield-p(k)*afield(k)
          ! ..corrections to electric forces
          CALL daxpy(2*ncpw%ngw*nstate,afield(k),ffield(1,1),1,c2,1)
       ENDIF
    ENDDO
    ! ..corrections to total energy 
    ener_com%etot=ener_com%etot+ener_com%eefield
    DEALLOCATE(ddmat_e,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ipiv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    EFIELD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE extfield
  ! ==================================================================

END MODULE efield_utils
