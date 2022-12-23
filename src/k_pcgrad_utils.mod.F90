MODULE k_pcgrad_utils
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE ortho_utils,                     ONLY: ortho,&
                                             preortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE pcgrad_driver,                   ONLY: eback,&
                                             precon
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: rnlsm
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k_pcgrad
  PUBLIC :: k_linesr

CONTAINS

  ! ==================================================================
  SUBROUTINE k_pcgrad(c0,c2,sc0,vpp,hnm1,rhoe,psi,tau0,nstate,dinit,ik)
    ! ==--------------------------------------------------------------==
    ! ==  PRECONDITIONED CONJUGATE GRADIENT OPTIMIZATION              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*), &
                                                sc0(nkpt%ngwk,*)
    REAL(real_8)                             :: vpp(*), rhoe(*)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: tau0(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: hnm1(nkpt%ngwk,nstate,*)
    LOGICAL                                  :: dinit
    INTEGER                                  :: ik

    CHARACTER(*), PARAMETER                  :: procedureN = 'k_pcgrad'

    EXTERNAL                                 :: ddot
    INTEGER                                  :: i, ierr, ihist, isub, nocc
    INTEGER, ALLOCATABLE, SAVE               :: istep(:)
    INTEGER, SAVE                            :: ifirst = 0 
    REAL(real_8)                             :: ddot, detot, gamma, ggnorm
    REAL(real_8), ALLOCATABLE, SAVE          :: ehist(:,:), fhist(:,:), &
                                                ghist(:,:)

    CALL tiset('  K_PCGRAD',isub)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Memory allocation
    IF (ifirst .EQ. 0) THEN
       ALLOCATE(ehist(10,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ghist(10,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fhist(10,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(istep(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(istep)!,nkpt%nkpts)
       ifirst = 1
    ENDIF

    nocc=0
    DO i=1,nstate
       IF (crge%f(i,ik).GT.1.e-5_real_8) nocc=nocc+1
    ENDDO

    IF (dinit) THEN
       CALL zeroing(istep)!,nkpt%nkpts)
       dinit=.FALSE.
    ENDIF

    istep(ik)=istep(ik)+1
    IF (cntl%prec) CALL precon(c2,vpp,crge%f,nkpt%ngwk,nstate,0)
    ggnorm=0.0_real_8
    DO i=1,nocc
       ggnorm=ggnorm+ddot(2*nkpt%ngwk,c2(1,i),1,c2(1,i),1)
    ENDDO
    CALL mp_sum(ggnorm,parai%allgrp)
    IF (istep(ik).LE.1) THEN
       ! ==--------------------------------------------------------------==
       ! ==  INITIAL CALL                                                ==
       ! ==--------------------------------------------------------------==
       ehist(1,ik)=ener_com%etot
       fhist(1,ik)=dt2bye
       ghist(1,ik)=ggnorm
       IF (cntl%pcgmin) CALL stopgm('k_pcgrad',&
            'Line Search not implemented yet',& 
            __LINE__,__FILE__)
       ! IF(cntl%pcgmin) THEN
       ! CALL DCOPY(2*NGW*NSTATE,C2(1,1),1,HNM1(1,1),1)
       ! CALL DSCAL(2*NGW*NSTATE,FHIST(1),HNM1(1,1),1)
       ! ALAM=1.0_real_8
       ! IF(PARENT) CALL EBACK(0)
       ! CALL LINESR(NSTATE,C0,HNM1,C2,SC0,TAU0,TSCR,
       ! &                RHOE,PSI,SCR,LSCR2,
       ! *                FHIST(1),DE,A,ALAM)
       ! IF(PARENT) THEN
       ! CALL EBACK(1)
       ! WRITE(6,'(A,G8.3,A,T55,G15.8)') 
       ! *        ' LINE SEARCH : LAMBDA=',ALAM*FHIST(1),
       ! *        '       ENERGY =',DE
       ! ENDIF
       ! ELSE
       CALL dcopy(2*nkpt%ngwk*nstate,c2(1,1),1,hnm1(1,1,ik),1)
       CALL daxpy(2*nkpt%ngwk*nstate,fhist(1,ik),c2(1,1),1,c0(1,1),1)
       ! ENDIF
    ELSE
       ! ==--------------------------------------------------------------==
       ! ==  CONJUGATE GRADIENT                                          ==
       ! ==--------------------------------------------------------------==
       IF (istep(ik).GT.10) THEN
          ihist=10
          DO i=1,9
             ehist(i,ik)=ehist(i+1,ik)
             fhist(i,ik)=fhist(i+1,ik)
             ghist(i,ik)=ghist(i+1,ik)
          ENDDO
       ELSE
          ihist=istep(ik)
       ENDIF
       ehist(ihist,ik)=ener_com%etot
       ghist(ihist,ik)=ggnorm
       detot=ehist(ihist,ik)-ehist(ihist-1,ik)
       IF (detot.GT.0._real_8) THEN
          fhist(ihist,ik)=0.75_real_8*fhist(ihist-1,ik)
       ELSE
          fhist(ihist,ik)=fhist(ihist-1,ik)
       ENDIF
       gamma=ghist(ihist,ik)/ghist(ihist-1,ik)
       ! IF(cntl%pcgmin) THEN
       ! DO I=1,NSTATE
       ! DO IG=1,NGW
       ! HNM1(IG,I)=C2(IG,I)+HNM1(IG,I)*GAMMA
       ! ENDDO
       ! ENDDO
       ! CALL LINESR(NSTATE,C0,HNM1,C2,SC0,TAU0,TSCR,
       ! &                RHOE,PSI,SCR,LSCR2,
       ! *                FHIST(IHIST),DE,A,ALAM)
       ! IF(PARENT) THEN
       ! WRITE(6,'(A,G8.3,A,T55,G15.8)') 
       ! *        ' LINE SEARCH : LAMBDA=',ALAM*FHIST(IHIST),
       ! *        '       ENERGY =',DE
       ! ENDIF
       ! ELSE
       CALL daxpy(2*nkpt%ngwk*nstate,gamma,hnm1(1,1,ik),1,c2(1,1),1)
       CALL dcopy(2*nkpt%ngwk*nstate,c2(1,1),1,hnm1(1,1,ik),1)
       CALL daxpy(2*nkpt%ngwk*nstate,fhist(ihist,ik),&
            hnm1(1,1,ik),1,c0(1,1),1)
       ! ENDIF
    ENDIF
    CALL tihalt('  K_PCGRAD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_pcgrad
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE k_linesr(nstate,x,p,cs,sc0,tau0,tscr,&
       rhoe,psi,fac,de,a,alam)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: x(ncpw%ngw,nstate), &
                                                p(ncpw%ngw,nstate), &
                                                cs(ncpw%ngw*nstate), &
                                                sc0(ncpw%ngw,nstate)
    REAL(real_8)                             :: tau0(3*maxsys%nax*maxsys%nsx),&
                                                tscr(3*maxsys%nax*maxsys%nsx),&
                                                rhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    REAL(real_8)                             :: fac, de, a, alam

    INTEGER                                  :: i
    REAL(real_8)                             :: dd, dx, e0, e1, e2, ee, x0

! ==--------------------------------------------------------------==

    IF (paral%parent) CALL eback(0)
    ! Here one should calculate the E0 relative to one single k_point
    ! with the wavefunction in the array X
    e0=ener_com%etot
    alam=1.0_real_8
    CALL dscal(2*ncpw%ngw*nstate,fac,p(1,1),1)
    dd=0.0_real_8
    DO i=1,nstate
       dd=dd+dotp(ncpw%ngw,p(:,i),p(:,i))
    ENDDO
    CALL mp_sum(dd,parai%allgrp)
    ! First point
    CALL dcopy(2*ncpw%ngw*nstate,x(1,1),1,cs(1),1)
    CALL daxpy(2*ncpw%ngw*nstate,alam,p(1,1),1,x(1,1),1)
    IF (cntl%nonort) THEN
       CALL stopgm('LINESR','NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSE
       IF (pslo_com%tivan) CALL rnlsm(x,nstate,1,1,.FALSE.)
       CALL preortho(x,nstate)
       CALL ortho(nstate,x,sc0)
       CALL rnlsm(x,nstate,1,1,.FALSE.)
       ! Here one should calculate the energy relative to one single k_point
       ! with the updated wavefunction copied in the array X
       ! CALL RSCPOT(X,TAU0,TSCR,RHOE,PSI,
       ! &                .FALSE.,.FALSE.,NSTATE,1)
    ENDIF
    IF (ener_com%etot.GT.e0) THEN
       e2=ener_com%etot
       alam=0.5_real_8
    ELSE
       e1=ener_com%etot
       alam=2.0_real_8
    ENDIF
    ! Second point
    CALL dcopy(2*ncpw%ngw*nstate,cs(1),1,x(1,1),1)
    CALL daxpy(2*ncpw%ngw*nstate,alam,p(1,1),1,x(1,1),1)
    IF (cntl%nonort) THEN
       CALL stopgm('LINESR','NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSE
       IF (pslo_com%tivan) CALL rnlsm(x,nstate,1,1,.FALSE.)
       CALL preortho(x,nstate)
       CALL ortho(nstate,x,sc0)
       CALL rnlsm(x,nstate,1,1,.FALSE.)
       ! Here one should calculate the energy relative to one single k_point
       ! with the updated wavefunction copied in the array X
       ! CALL RSCPOT(X,TAU0,TSCR,RHOE,PSI,
       ! &                .FALSE.,.FALSE.,NSTATE,1)
    ENDIF
    IF (alam.GT.1._real_8) THEN
       e2=ener_com%etot
       dx=1.0_real_8
    ELSE
       e1=ener_com%etot
       dx=0.5_real_8
    ENDIF
    ! Extrapolation
    IF ((2*e0-4*e1+2*e2).EQ.0._real_8) THEN
       ! Very rare (T.D.)
       x0=dx
    ELSE
       x0=dx*(4*e1-e2-3*e0)/(2*e0-4*e1+2*e2)
    ENDIF
    a=(e1-e0)/dx/(2*x0+dx)
    ee=e0-a*x0*x0
    alam=-x0
    IF (alam.GT.0._real_8) THEN
       IF (alam.GT.3._real_8) alam=3.0_real_8
    ELSE
       IF (e2.LT.e0) THEN
          alam=2*dx
       ELSEIF (e1.LT.e0) THEN
          alam=dx
       ELSE
          alam=0.05_real_8
       ENDIF
    ENDIF
    de=ee+a*(x0+alam)**2
    a=a/dd
    CALL dcopy(2*ncpw%ngw*nstate,cs(1),1,x(1,1),1)
    CALL daxpy(2*ncpw%ngw*nstate,alam,p(1,1),1,x(1,1),1)
    CALL dscal(2*ncpw%ngw*nstate,1._real_8/fac,p(1,1),1)
    IF (alam.LT.0.1_real_8) THEN
       fac=fac*0.25_real_8
    ELSEIF (alam.GT.1.00_real_8) THEN
       fac=fac*1.50_real_8
    ELSEIF (alam.LT.0.20_real_8) THEN
       fac=fac*0.75_real_8
    ENDIF
    IF (paral%parent) CALL eback(1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_linesr
  ! ==================================================================
  ! ==================================================================

END MODULE k_pcgrad_utils
