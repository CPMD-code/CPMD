#if defined(__BG)
#define __VECTOR
#endif
C     ==================================================================
      SUBROUTINE mm_solv_const(str,dt,taup,velp,tau0)
C     ==--------------------------------------------------------------== 
cmb - Revised on 13 June 2006
C includes
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use prmem_utils, only: prmem
      use system
      use rmas
      use mm_dimmod
      use cnst
      use azzero_utils, only:azzero=>azzero_dont_use_me_anymore
      implicit none
c
cmb      integer maxsys%nax ! already declared
      real*8 taup(3,maxsys%nax,*)
      real*8 tau0(3,maxsys%nax,*)
      real*8 velp(3,maxsys%nax,*)
      real*8 dt
      character str*(*)
c
      integer, save ::  NCONST
      integer, allocatable, save :: ABl(:)
      integer, allocatable, save :: BBl(:)
      real*8, allocatable, save :: DAB(:)
      real*8, allocatable, save :: INVM(:)
      real*8, allocatable, save :: xs0(:,:),xsp(:,:),vsp(:,:) 
c
c locals
c Must be SCMASS !!!
C      REAL*8     mass_prot
C      parameter(mass_prot=1822.89D0)
      real*8, allocatable :: wmass(:)
      REAL*8     tol
      parameter(tol=1.D-12)
      integer NSM
      real*8 rij
      integer ia,ja,ic,ii,i,sn,NRAM,ierr
      logical, save :: first = .true.
      character(*),parameter :: procedureN='mm_solv_const'
c
      NRAM=solvvv%NRAM_gr
      NSM=solsolv%NSOLV/NRAM

C     nothing to do without solvent molecules
      if(NSM.eq.0)return

      allocate(WMASS(NRAM))
      call azzero(WMASS,NRAM)
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif 
      DO i=1,NRAM
        WMASS(i)=rmass%PMA0(cpsp(solsolv%NRPT+i))
      ENDDO

      if(first)then
        first=.false.
        NCONST=0
        do ia=1,NRAM
          do ja=ia+1,min(ia+4,NRAM)
            NCONST=NCONST+1
          enddo
        enddo
        allocate(ABl(NCONST),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(BBl(NCONST),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(DAB(NCONST),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(INVM(NRAM),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(xs0(3,NRAM),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(xsp(3,NRAM),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        allocate(vsp(3,NRAM),STAT=ierr)
        if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
        call prmem('MM_SOLV_CO')
cmb        ii=NRPT
!$OMP parallel do private(ia,ii)
#ifdef __SR11000
*poption parallel, tlocal(ia,ii)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
        DO ia=1,NRAM
cmb          ii=ii+1
          ii=solsolv%NRPT+ia
          xs0(1,ia)= tau0(1,cpat(ii),cpsp(ii))
          xs0(2,ia)= tau0(2,cpat(ii),cpsp(ii))
          xs0(3,ia)= tau0(3,cpat(ii),cpsp(ii))
        ENDDO
        ic=0
        DO ia=1,NRAM
          INVM(ia)=1.D0/WMASS(ia)/SCMASS
          DO ja=ia+1,min(ia+4,NRAM)
            ic=ic+1
            rij=(xs0(1,ia)-xs0(1,ja))**2
     .         +(xs0(2,ia)-xs0(2,ja))**2
     .         +(xs0(3,ia)-xs0(3,ja))**2
            Abl(ic)=ia
            Bbl(ic)=ja
            DAB(ic)=rij
          ENDDO
        ENDDO
      endif
c
      if(str(1:1).eq.'S')then
#ifdef __VECTOR
        DO i=1,NSM
!$OMP parallel do private(ia,ii) 
#ifdef __NEC
!CDIR NODEP
#endif
#ifdef __SR11000
*poption parallel, tlocal(ia,ii)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
          DO ia=1,NRAM
            ii=solsolv%NRPT+ia+(i-1)*NRAM
            xs0(1,ia)= tau0(1,cpat(ii),cpsp(ii))
            xs0(2,ia)= tau0(2,cpat(ii),cpsp(ii))
            xs0(3,ia)= tau0(3,cpat(ii),cpsp(ii))
            xsp(1,ia)= taup(1,cpat(ii),cpsp(ii))
            xsp(2,ia)= taup(2,cpat(ii),cpsp(ii))
            xsp(3,ia)= taup(3,cpat(ii),cpsp(ii))
            vsp(1,ia)= velp(1,cpat(ii),cpsp(ii))
            vsp(2,ia)= velp(2,cpat(ii),cpsp(ii))
            vsp(3,ia)= velp(3,cpat(ii),cpsp(ii))
          ENDDO
          call shake_sm (NRAM,xsp,vsp,xs0,
     c         NCONST,Abl,Bbl,DAB,INVM,TOL,dt)
!$OMP parallel do private(ia,ii)
#ifdef __SR11000
*poption parallel, tlocal(ia,ii)
#endif 
          DO ia=1,NRAM
            ii=solsolv%NRPT+ia+(i-1)*NRAM
            taup(1,cpat(ii),cpsp(ii)) = xsp(1,ia)
            velp(1,cpat(ii),cpsp(ii)) = vsp(1,ia)
            taup(2,cpat(ii),cpsp(ii)) = xsp(2,ia)
            velp(2,cpat(ii),cpsp(ii)) = vsp(2,ia)
            taup(3,cpat(ii),cpsp(ii)) = xsp(3,ia)
            velp(3,cpat(ii),cpsp(ii)) = vsp(3,ia)
          ENDDO
        ENDDO
#else
        ii=solsolv%NRPT
        DO i=1,NSM
          DO ia=1,NRAM
            ii=ii+1
            xs0(1,ia)= tau0(1,cpat(ii),cpsp(ii))
            xs0(2,ia)= tau0(2,cpat(ii),cpsp(ii))
            xs0(3,ia)= tau0(3,cpat(ii),cpsp(ii))
            xsp(1,ia)= taup(1,cpat(ii),cpsp(ii))
            xsp(2,ia)= taup(2,cpat(ii),cpsp(ii))
            xsp(3,ia)= taup(3,cpat(ii),cpsp(ii))
            vsp(1,ia)= velp(1,cpat(ii),cpsp(ii))
            vsp(2,ia)= velp(2,cpat(ii),cpsp(ii))
            vsp(3,ia)= velp(3,cpat(ii),cpsp(ii))
          ENDDO
          call shake_sm (NRAM,xsp,vsp,xs0,
     c                 NCONST,Abl,Bbl,DAB,
     c                 INVM,TOL,dt)
          ii=ii-NRAM
          DO ia=1,NRAM
            ii=ii+1
            taup(1,cpat(ii),cpsp(ii)) = xsp(1,ia)
            taup(2,cpat(ii),cpsp(ii)) = xsp(2,ia)
            taup(3,cpat(ii),cpsp(ii)) = xsp(3,ia)
            velp(1,cpat(ii),cpsp(ii)) = vsp(1,ia)
            velp(2,cpat(ii),cpsp(ii)) = vsp(2,ia)
            velp(3,cpat(ii),cpsp(ii)) = vsp(3,ia)
          ENDDO
        ENDDO
#endif
      else if(str(1:1).eq.'R')then
#ifdef __VECTOR
        DO i=1,NSM
!$OMP parallel do private(ia,ii)
#ifdef __NEC
!CDIR NODEP
#endif
#ifdef __SR11000
*poption parallel, tlocal(ia,ii)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
          DO ia=1,NRAM
            ii=solsolv%NRPT+ia+(i-1)*NRAM
            xsp(1,ia)= taup(1,cpat(ii),cpsp(ii))
            vsp(1,ia)= velp(1,cpat(ii),cpsp(ii))
            xsp(2,ia)= taup(2,cpat(ii),cpsp(ii))
            vsp(2,ia)= velp(2,cpat(ii),cpsp(ii))
            xsp(3,ia)= taup(3,cpat(ii),cpsp(ii))
            vsp(3,ia)= velp(3,cpat(ii),cpsp(ii))
          ENDDO
!         write(6,*)'RATTLE',i,ii
          call rattle_sm (NRAM,xsp,vsp,NCONST,Abl,Bbl,DAB,INVM,TOL)
!$OMP parallel do private(ia,ii)
#ifdef __SR11000
*poption parallel, tlocal(ia,ii)
#endif 
          DO ia=1,NRAM
            ii=solsolv%NRPT+ia+(i-1)*NRAM
            velp(1,cpat(ii),cpsp(ii)) = vsp(1,ia)
            velp(2,cpat(ii),cpsp(ii)) = vsp(2,ia)
            velp(3,cpat(ii),cpsp(ii)) = vsp(3,ia)
          ENDDO
        ENDDO
#else
        ii=solsolv%NRPT
        DO i=1,NSM
          DO ia=1,NRAM
            ii=ii+1
            xsp(1,ia)= taup(1,cpat(ii),cpsp(ii))
            xsp(2,ia)= taup(2,cpat(ii),cpsp(ii))
            xsp(3,ia)= taup(3,cpat(ii),cpsp(ii))
            vsp(1,ia)= velp(1,cpat(ii),cpsp(ii))
            vsp(2,ia)= velp(2,cpat(ii),cpsp(ii))
            vsp(3,ia)= velp(3,cpat(ii),cpsp(ii))
!            write(6,*)1,ia,xsp(1,ia),vsp(1,ia)
!            write(6,*)2,ia,xsp(2,ia),vsp(2,ia)
!            write(6,*)3,ia,xsp(3,ia),vsp(3,ia)
          ENDDO
!         write(6,*)'RATTLE',i,ii
          call rattle_sm (NRAM,xsp,vsp,NCONST,Abl,Bbl,DAB,INVM,TOL)
          ii=ii-NRAM
          DO ia=1,NRAM
            ii=ii+1
            velp(1,cpat(ii),cpsp(ii)) = vsp(1,ia)
            velp(2,cpat(ii),cpsp(ii)) = vsp(2,ia)
            velp(3,cpat(ii),cpsp(ii)) = vsp(3,ia)
          ENDDO
        ENDDO
#endif
      else
        call stopgm('MM_SOLV_CONST','wrong call parameter',
     & __LINE__,"__file__")
      endif

      DEALLOCATE(wmass)
      RETURN
      END
C     ==================================================================
      SUBROUTINE shake_sm (NA,TAUP,VELP,TAU0,NCONST,Abl,Bbl,DAB,
     .           INVM,TOL,dt)
C     ==--------------------------------------------------------------== 
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      implicit none
C
      integer NA           !number of atoms involved in the constraints
      real*8  TAUP(3,*),VELP(3,*),TAU0(3,*)
      real*8  INVM(*)     !their inverse mass
      integer NCONST       !number of distance constraints
      integer Abl(*)  !first atom involved in the constraint
      integer Bbl(*)  !second atom involved in the constraint
      real*8  DAB(*)  !their distance  (squared)
      real*8  TOL
      real*8  dt           !time step
c
      integer MAXIT
      PARAMETER (MAXIT=10000)
      real*8  RPTOL
      PARAMETER ( RPTOL = 1.0D-6 )
c
      INTEGER i,ix,IT,ibl,A,B
      LOGICAL DONE,MOVING(NA),MOVED(NA)
      REAL*8  R0AB(3,NCONST)
      REAL*8  RMA, RMB,GAB,DX(3)
      REAL*8  DIFFSQ,RPAB,TOL2,RABSQ,RABSQ0,r_AB(3)
C    *******************************************************************
      TOL2 = 2.0D0*TOL
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
      DO i = 1, NA
        MOVING(i) = .FALSE.
        MOVED(i)  = .TRUE.
      ENDDO
C
      IT = 0
      DONE = .FALSE.
!$OMP parallel do private(ibl,A,B)
      DO ibl = 1, NCONST
        A=Abl(ibl)
        B=Bbl(ibl)
        R0AB(1,ibl) = TAU0(1,A) - TAU0(1,B)
        R0AB(2,ibl) = TAU0(2,A) - TAU0(2,B)
        R0AB(3,ibl) = TAU0(3,A) - TAU0(3,B)
      ENDDO
C
1000  IF ((.not.done) .and. (it.le.maxit) ) THEN
        done = .true.
        DO ibl = 1, NCONST
          A=Abl(ibl)
          B=Bbl(ibl)
          RABSQ0  = DAB(ibl)
          IF ( MOVED(A) .OR. MOVED(B) ) THEN
            r_AB(1)=TAUP(1,A) - TAUP(1,B)
            r_AB(2)=TAUP(2,A) - TAUP(2,B)
            r_AB(3)=TAUP(3,A) - TAUP(3,B)
            RABSQ=r_AB(1)*r_AB(1) + 
     .            r_AB(2)*r_AB(2) +
     .            r_AB(3)*r_AB(3)
            DIFFSQ = RABSQ0 - RABSQ
            IF ( DABS(DIFFSQ) .GT. (RABSQ*TOL2) ) THEN
              RPAB = R0AB(1,ibl)*r_AB(1)+R0AB(2,ibl)*r_AB(2)
     &              +R0AB(3,ibl)*r_AB(3)
              IF( RPAB .LT. ( RABSQ0*RPTOL ))THEN
                write(6,*)RPAB, (R0AB(ix,ibl),r_AB(ix),ix=1,3)
                call stopgm('SHAKE_SM','CONSTRAINT FAILURE',
     & __LINE__,"__file__")
              ENDIF
              RMA = INVM(A)
              RMB = INVM(B)
              GAB = DIFFSQ/(2.0D0*(RMA+RMB)*RPAB)
#ifdef __VECTOR
              DX(1) = R0AB(1,ibl) * GAB
              TAUP(1,A) = TAUP(1,A) + RMA * DX(1)
              TAUP(1,B) = TAUP(1,B) - RMB * DX(1)
              VELP(1,A) = VELP(1,A) + RMA * DX(1)/dt
              VELP(1,B) = VELP(1,B) - RMB * DX(1)/dt
              DX(2) = R0AB(2,ibl) * GAB
              TAUP(2,A) = TAUP(2,A) + RMA * DX(2)
              TAUP(2,B) = TAUP(2,B) - RMB * DX(2)
              VELP(2,A) = VELP(2,A) + RMA * DX(2)/dt
              VELP(2,B) = VELP(2,B) - RMB * DX(2)/dt
              DX(3) = R0AB(3,ibl) * GAB
              TAUP(3,A) = TAUP(3,A) + RMA * DX(3)
              TAUP(3,B) = TAUP(3,B) - RMB * DX(3)
              VELP(3,A) = VELP(3,A) + RMA * DX(3)/dt
              VELP(3,B) = VELP(3,B) - RMB * DX(3)/dt
#else
              DO ix=1,3
                DX(ix) = R0AB(ix,ibl)*GAB
                TAUP(ix,A) = TAUP(ix,A) + RMA*DX(ix)
                TAUP(ix,B) = TAUP(ix,B) - RMB*DX(ix)
                VELP(ix,A) = VELP(ix,A) + RMA*DX(ix)/dt
                VELP(ix,B) = VELP(ix,B) - RMB*DX(ix)/dt
              ENDDO
#endif
              MOVING(A) = .TRUE.
              MOVING(B) = .TRUE.
              DONE = .FALSE.
            ENDIF
          ENDIF
        ENDDO
C
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
        DO i = 1, NA
          MOVED(i)  = MOVING(i)
          MOVING(i) = .FALSE.
        ENDDO
        IT = IT + 1
        GOTO 1000
      ENDIF
C
      IF (.NOT. DONE) THEN
        WRITE(6,*) "MOLECULE",I,it,DABS(DIFFSQ),RABSQ*TOL2
        call stopgm('SHAKE_SM','TOO MANY CONSTRAINT ITERATIONS',
     & __LINE__,"__file__")
      ENDIF
      RETURN
      END
C     ==================================================================
      SUBROUTINE RATTLE_SM(NA,TAUP,VELP,NCONST,Abl,Bbl,DAB,INVM,TOL)
C     ==--------------------------------------------------------------== 
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      implicit none
C
      integer NA           !number of atoms involved in the costraints
      real*8  TAUP(3,*),VELP(3,*)
      real*8  INVM(*)     !their inverse mass
      integer NCONST       !number of distance constraints
      integer Abl(*)  !first atom involved in the constraint
      integer Bbl(*)  !second atom involved in the constraint
      real*8  DAB(*)  !their distance  (squared)
      real*8  TOL
c
      integer MAXIT
      PARAMETER (MAXIT=1000)
      integer i,ix,IT,ibl,A,B
      logical DONE,MOVING(NA),MOVED(NA)
      real*8  RMA, RMB,GAB,DX(3)
      real*8  RVAB,r_AB(3,NCONST),v_AB(3)
C     *******************************************************************
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
      DO i = 1, NA
        MOVING(i) = .FALSE.
        MOVED(i)  = .TRUE.
      ENDDO
!$OMP parallel do private(ibl,A,B)
      DO ibl = 1, NCONST
        A=Abl(ibl)
        B=Bbl(ibl)
        r_AB(1,ibl)=TAUP(1,A) - TAUP(1,B)
        r_AB(2,ibl)=TAUP(2,A) - TAUP(2,B)
        r_AB(3,ibl)=TAUP(3,A) - TAUP(3,B)
      ENDDO
C
      IT = 0
      DONE = .FALSE.
C
1000  IF ((.not.done) .and. (it.le.maxit) ) THEN
        done = .true.
        DO ibl = 1, NCONST
          A=Abl(ibl)
          B=Bbl(ibl)
          IF ( MOVED(A) .OR. MOVED(B) ) THEN
            RVAB=0.D0
            v_AB(1)=VELP(1,A) - VELP(1,B)
            v_AB(2)=VELP(2,A) - VELP(2,B)
            v_AB(3)=VELP(3,A) - VELP(3,B)
            RVAB=r_AB(1,ibl)*v_AB(1) +
     .           r_AB(2,ibl)*v_AB(2) +
     .           r_AB(3,ibl)*v_AB(3)
            RMA = INVM(A)
            RMB = INVM(B)
            GAB  = -RVAB/((RMA+RMB)*DAB(ibl))
            IF ( DABS(GAB) .GT. TOL ) THEN
#ifdef __VECTOR
              DX(1)  = r_AB(1,ibl) * GAB
              DX(2)  = r_AB(2,ibl) * GAB
              DX(3)  = r_AB(3,ibl) * GAB
              VELP(1,A)= VELP(1,A) + RMA*DX(1)
              VELP(2,A)= VELP(2,A) + RMA*DX(2)
              VELP(3,A)= VELP(3,A) + RMA*DX(3)
              VELP(1,B)= VELP(1,B) - RMB*DX(1)
              VELP(2,B)= VELP(2,B) - RMB*DX(2)
              VELP(3,B)= VELP(3,B) - RMB*DX(3)
#else
              DO ix=1,3
                DX(ix)  = r_AB(ix,ibl) * GAB
                VELP(ix,A)= VELP(ix,A) + RMA*DX(ix)
                VELP(ix,B)= VELP(ix,B) - RMB*DX(ix)
              ENDDO
#endif
              MOVING(A) = .TRUE.
              MOVING(B) = .TRUE.
              DONE = .FALSE.
            ENDIF
          ENDIF
        ENDDO
C
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif
        DO i = 1, NA
          MOVED(i)  = MOVING(i)
          MOVING(i) = .FALSE.
        ENDDO
        IT = IT + 1
        goto 1000
      ENDIF
C
      IF (.NOT. DONE) THEN
        WRITE(6,*) "MOLECULE",I,it,DABS(GAB),TOL
        call stopgm('RATTLE_SM','TOO MANY CONSTRAINT ITERATIONS',
     & __LINE__,"__file__")
      ENDIF
      RETURN
      END
