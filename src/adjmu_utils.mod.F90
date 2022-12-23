MODULE adjmu_utils
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: adjmu
  PUBLIC :: occup
  PUBLIC :: enert
  PUBLIC :: convfrie

CONTAINS

  ! ==================================================================
  SUBROUTINE adjmu(nstate,nkpoint,nel,betael,we,wk,tspin,amu,tinfo)
    ! ==--------------------------------------------------------------==
    ! == ADJUST THE CHEMICAL POTENTIAL (AMU)                          ==
    ! == WITH THE BISECTION METHOD                                    ==
    ! == WE DOESN''T NEED BE SORTED                                   ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==    NSTATE number of states                                   ==
    ! ==    NKPOINT number of k points                                ==
    ! ==    NEL number of electrons                                   ==
    ! ==    BETAEL 1/kT                                               ==
    ! ==    WE(NSTATE,NKPOINT) eigenvalues (energies for each state)  ==
    ! ==    WK(NKPOINT) weight of k points                            ==
    ! ==    TSPIN .TRUE. if  local spin polarization                  ==
    ! == OUTPUT:                                                      ==
    ! ==    AMU chemical potential                                    ==
    ! == INPUT/OUTPUT:                                                ==
    ! ==    TINFO if .TRUE. return .TRUE. if ADJMU finds AMU          ==
    ! ==          if .FALSE. stops if ADJMU fails                     ==
    ! ==--------------------------------------------------------------==
    ! == myprecision GIVEN DELTA = 1.e-6_real_8 EV  (0.1Kelvin)       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: nel, betael, &
                                                we(nstate,nkpoint), &
                                                wk(nkpoint)
    LOGICAL                                  :: tspin
    REAL(real_8)                             :: amu
    LOGICAL                                  :: tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'adjmu'
    INTEGER, PARAMETER                       :: nit_max = 200
    REAL(KIND=real_8), PARAMETER :: precpar = 1.e-10_real_8, &
      ry = 13.605804_real_8, delta = 1.e-6_real_8/(2._real_8*ry)

    INTEGER                                  :: i, ik, isub, nit
    LOGICAL                                  :: not_converged
    REAL(real_8)                             :: amu1, amu2, damu, &
                                                myprecision, rhint, rhint1, &
                                                rhint2

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! == STARTING CHEMICAL POTENTIAL AS THE AVERAGE BETWEEN           ==
    ! == LAST OCCUPIED STATE AND THE FIRST EMPTY STATE                ==
    ! ==--------------------------------------------------------------==
    ! AMU1(lowest)...AMU...AMU2(biggest)...
    myprecision=precpar*LOG(nel+1._real_8)
    amu1=we(1,1)
    amu2=we(1,1)
    DO ik=1,nkpoint
       DO i=1,nstate
          IF (amu1.GT.we(i,ik)) THEN
             amu1=we(i,ik)
          ENDIF
          IF (amu2.LT.we(i,ik)) THEN
             amu2=we(i,ik)
          ENDIF
       ENDDO
    ENDDO
    rhint1=1._real_8
    DO WHILE(rhint1.GT.0._real_8)
       rhint1=rhoint(we,wk,amu1,betael,nstate,nkpoint,nel,tspin)
       IF (rhint1.GT.0._real_8) amu1=-2._real_8*ABS(amu1)
    ENDDO
    rhint2=rhoint(we,wk,amu2,betael,nstate,nkpoint,nel,tspin)
    IF (rhint2.LT.0._real_8) THEN
       damu=amu2-amu1
       rhint=rhint1
       not_converged=.FALSE.
       GOTO 1000
    ENDIF
    damu=amu2-amu1
    amu=amu1
    ! ==--------------------------------------------------------------==
    DO nit=1,nit_max
       amu=amu1+0.5_real_8*damu
       rhint=rhoint(we,wk,amu,betael,nstate,nkpoint,nel,tspin)
       IF (damu.LT.delta) THEN
          IF (ABS(rhint).LT.myprecision) THEN
             CALL tihalt(procedureN,isub)
             RETURN
          ELSEIF (ABS(rhint1).LT.myprecision) THEN
             amu=amu1
             CALL tihalt(procedureN,isub)
             RETURN
          ELSEIF (ABS(rhint2).LT.myprecision) THEN
             amu=amu2
             CALL tihalt(procedureN,isub)
             RETURN
          ENDIF
       ENDIF
       ! LT in order to have the lowest possible AMU (if gap).
       IF (rhint.LT.0._real_8) THEN
          amu1=amu
          rhint1=rhint
       ELSE
          amu2=amu
          rhint2=rhint
       ENDIF
       damu=amu2-amu1
    ENDDO
    not_converged=.TRUE.
1000 CONTINUE
    IF (tinfo) THEN
       tinfo=.FALSE.
       CALL tihalt(procedureN,isub)
       RETURN
    ELSE
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' ADJMU|    IK     I    EIGENVALUES:'
          DO ik=1,nkpoint
             DO i=1,nstate
                IF (paral%io_parent)&
                     WRITE(6,'(7x,2i6,2x,f14.8)') ik,i,we(i,ik)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,*) '        AMU1=',amu1,'       RHINT1=',rhint1
          IF (paral%io_parent)&
               WRITE(6,*) '        AMU2=',amu2,'       RHINT2=',rhint2
          IF (paral%io_parent)&
               WRITE(6,*) '        DAMU=',damu,'       RHINT =',rhint
          IF (paral%io_parent)&
               WRITE(6,*) '        RHINT-NEL=',rhint
          IF (not_converged) THEN
             CALL stopgm('ADJMU','BISECTION DID NOT CONVERGE',& 
                  __LINE__,__FILE__)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,*) ' ADJMU! THE NUMBER OF STATES [',&
                  nstate,'] IS TOO SMALL'
             CALL stopgm('ADJMU','BISECTION COULD NOT CONVERGE',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE adjmu
  ! ==================================================================
  FUNCTION rhoint(we,wk,amu,betael,nstate,nkpoint,nel,tlsd)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE NUMBER OF ELECTRONS FOR THE GIVEN POTENTIAL   ==
    ! == WORKS WITH T=0K (BETAEL almost infinite) and with the last   ==
    ! == state degenerated                                            ==
    ! ==--------------------------------------------------------------==
    ! == WE(NSTATE,NKPOINT) eigenvalues                               ==
    ! == WK(NKPOINT) weight of k points                               ==
    ! == AMU the given chemical potential                             ==
    ! == BETAEL 1/kT                                                  ==
    ! == NSTATE number of states                                      ==
    ! == NKPOINT number of kpoints                                    ==
    ! == NEL number of electrons                                      ==
    ! == tlsd if false 2 electrons per state                          ==
    ! ==--------------------------------------------------------------==
    ! == myprecision GIVEN DELTA = 1.e-6_real_8 EV  (0.1Kelvin)       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu, betael
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                wk(nkpoint), nel
    LOGICAL                                  :: tlsd
    REAL(real_8)                             :: rhoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhoint'
    REAL(KIND=real_8), PARAMETER :: ry = 13.606_real_8, &
      delta = 1.e-6_real_8/(2._real_8*ry), tollpar = 1.e-13_real_8

    INTEGER                                  :: i, ik, isub
    LOGICAL                                  :: tamu
    REAL(real_8)                             :: arg, argmax, argmin, dramu, &
                                                drwe, remain, wbigwe, xlm

    CALL tiset(procedureN,isub)
    tamu = .FALSE.
    argmin=LOG(tollpar)
    argmax=-argmin
    dramu=dround(amu,delta)
    ! If last occupied states are degenerated, there should be some pb.
    DO ik=1,nkpoint
       DO i=1,nstate
          drwe=dround(we(i,ik),delta)
          IF (ABS(drwe-dramu).LT.delta) THEN
             tamu=.TRUE.
          ENDIF
       ENDDO
    ENDDO
    ! Find the degeneracy of last occupied state.
    wbigwe=0._real_8
    rhoint=0._real_8
    DO ik=1,nkpoint
       DO i=1,nstate
          arg=-betael*(we(i,ik)-amu)
          IF (arg.GT.argmax) THEN
             xlm=1._real_8
          ELSE IF (arg.LT.argmin) THEN
             xlm=0._real_8
          ELSE
             xlm=1._real_8/(EXP(-arg)+1._real_8)
          ENDIF
          drwe=dround(we(i,ik),delta)
          IF (tamu.AND.ABS(drwe-dramu).LT.delta) THEN
             wbigwe=wbigwe+wk(ik)
          ELSEIF (xlm.GT.0._real_8) THEN
             rhoint=rhoint+wk(ik)*xlm
          ENDIF
       ENDDO
    ENDDO
    ! Solve some small problem.
    IF (.NOT.tlsd) THEN
       rhoint=2._real_8*rhoint
    ENDIF
    remain=nel-rhoint
    IF (remain.LE.0._real_8) THEN
       IF (tlsd) THEN
          remain=1._real_8
       ELSE
          remain=2._real_8
       ENDIF
    ENDIF
    IF (wbigwe.GT.0._real_8) THEN
       remain=remain/wbigwe
       IF (tlsd) THEN
          remain=MIN(remain,1._real_8)
       ELSE
          remain=MIN(remain,2._real_8)
       ENDIF
       DO ik=1,nkpoint
          DO i=1,nstate
             drwe=dround(we(i,ik),delta)
             IF (tamu.AND.ABS(drwe-dramu).LT.delta) THEN
                rhoint=rhoint+wk(ik)*remain
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    rhoint=rhoint-nel
    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    RETURN
  END FUNCTION rhoint
  ! ==================================================================
  SUBROUTINE occup(we,wk,amu,betael,focc,nstate,nkpoint,nel,tlsd)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                                                              ==
    ! ==  THE OCCUPATION NUMBERS: F = 2/(1+LAMBDA^P*EXP(BETA MU))     ==
    ! ==  WORKS WITH T=0K AND THE LAST STATE DEGENERATED              ==
    ! ==--------------------------------------------------------------==
    ! == myprecision GIVEN DELTA = 1.e-6_real_8 EV  (0.1Kelvin)       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu, betael
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                focc(nstate,nkpoint), &
                                                wk(nkpoint), nel
    LOGICAL                                  :: tlsd

    REAL(KIND=real_8), PARAMETER :: ry = 13.606_real_8, &
      delta = 1.e-6_real_8/(2._real_8*ry), tollpar = 1.e-13_real_8

    INTEGER                                  :: i, ik
    LOGICAL                                  :: tamu
    REAL(real_8)                             :: arg, argmax, argmin, dramu, &
                                                drwe, remain, rhoint_, &
                                                wbigwe, xlm

    tamu = .FALSE.
    argmin=LOG(tollpar)
    argmax=-argmin
    dramu=dround(amu,delta)
    ! If last occupied states are degenerated, there should be some pb.
    DO ik=1,nkpoint
       DO i=1,nstate
          drwe=dround(we(i,ik),delta)
          IF (ABS(drwe-dramu).LT.delta) THEN
             tamu=.TRUE.
          ENDIF
       ENDDO
    ENDDO
    ! Find the degeneracy of last occupied state.
    wbigwe=0._real_8
    rhoint_=0._real_8
    IF (tlsd) THEN ! LSD
       DO ik=1,nkpoint
          DO i=1,nstate
             arg=-betael*(we(i,ik)-amu)
             IF (arg.GT.argmax) THEN
                xlm=1._real_8
             ELSE IF (arg.LT.argmin) THEN
                xlm=0._real_8
             ELSE
                xlm=1._real_8/(EXP(-arg)+1._real_8)
             ENDIF
             drwe=dround(we(i,ik),delta)
             IF (tamu.AND.ABS(drwe-dramu).LT.delta) THEN
                wbigwe=wbigwe+wk(ik)
             ELSEIF (xlm.GT.0._real_8) THEN
                rhoint_=rhoint_+wk(ik)*xlm
             ENDIF
             focc(i,ik)=xlm
          ENDDO
       ENDDO
    ELSE ! LDA
       DO ik=1,nkpoint
          DO i=1,nstate
             arg=-betael*(we(i,ik)-amu)
             IF (arg.GT.argmax) THEN
                xlm=1._real_8
             ELSE IF (arg.LT.argmin) THEN
                xlm=0._real_8
             ELSE
                xlm=1._real_8/(EXP(-arg)+1._real_8)
             ENDIF
             drwe=dround(we(i,ik),delta)
             IF (tamu.AND.ABS(drwe-dramu).LT.delta) THEN
                wbigwe=wbigwe+wk(ik)
             ELSEIF (xlm.GT.0._real_8) THEN
                rhoint_=rhoint_+wk(ik)*xlm
             ENDIF
             focc(i,ik)=2._real_8*xlm
          ENDDO
       ENDDO
    ENDIF ! LSD/LDA select

    ! Solve some small problem.
    IF (.NOT.tlsd) THEN
       rhoint_=2._real_8*rhoint_
    ENDIF
    remain=nel-rhoint_
    IF (remain.LE.0._real_8) THEN
       IF (tlsd) THEN
          remain=1._real_8
       ELSE
          remain=2._real_8
       ENDIF
    ENDIF
    IF (wbigwe.GT.0._real_8) THEN
       remain=remain/wbigwe
       IF (tlsd) THEN
          remain=MIN(remain,1._real_8)
       ELSE
          remain=MIN(remain,2._real_8)
       ENDIF
       DO ik=1,nkpoint
          DO i=1,nstate
             drwe=dround(we(i,ik),delta)
             IF ( tamu.AND.ABS(drwe-dramu).LT.delta ) THEN
                focc(i,ik)=remain
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE occup
  ! ==================================================================
  SUBROUTINE enert(we,wk,f,amu,betael,ns,nstate,nkpoint,nel,tlsd,&
       free_energy,band_energy,mkt_entropy)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES:                                                    ==
    ! == - FREE_ENERGY - THE ELECTRONIC FREE ENERGY                   ==
    ! == (see A. Alavi, J. Kohanoff, M. Parrinello, D. Frenkel,       ==
    ! ==  PRL 732599 (1994)),                                        ==
    ! == - BAND_ENERGY - THE BAND ENERGY (sum of eigenvalues x occ. n)==
    ! == - MKT_ENTROPY = - T*S                                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu, betael
    INTEGER                                  :: ns, nstate, nkpoint
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                f(nstate,nkpoint), &
                                                wk(nkpoint), nel
    LOGICAL                                  :: tlsd
    REAL(real_8)                             :: free_energy, band_energy, &
                                                mkt_entropy

    REAL(KIND=real_8), PARAMETER             :: betaelmax = 1.e33_real_8, &
                                                tollpar = 1.e-13_real_8

    INTEGER                                  :: i, ik
    REAL(real_8)                             :: arg, argmax, argmin, gibbs, &
                                                gpband, occup, spind, &
                                                spindinv, temp_en

    argmin=LOG(tollpar)
    argmax=-argmin
    gpband = 0._real_8
    band_energy = 0._real_8
    mkt_entropy = 0._real_8
    IF (tlsd) THEN
       spind = 1._real_8
    ELSE
       spind = 2._real_8
    ENDIF
    spindinv=1._real_8/spind
    DO ik=1,nkpoint
       temp_en = 0._real_8
       DO i=1,ns
          arg=-betael*(we(i,ik)-amu)
          IF (arg.GT.argmax) THEN
             gpband=gpband+wk(ik)*spind*(we(i,ik)-amu)
          ELSEIF (arg.LT.argmin) THEN
             ! Do nothing (avoid underflow).
          ELSE
             gpband=gpband-wk(ik)*spind*LOG(EXP(arg)+1._real_8)/betael
          ENDIF
          band_energy = band_energy + wk(ik)*f(i,ik)*we(i,ik)
          occup=f(i,ik)*spindinv
          IF (occup.GT.0._real_8) THEN
             temp_en = temp_en + occup*LOG(occup)
          ENDIF
          IF (occup.LT.1._real_8) THEN
             temp_en = temp_en + (1._real_8-occup)*LOG(1._real_8-occup)
          ENDIF
       ENDDO
       mkt_entropy = mkt_entropy + wk(ik)*temp_en
    ENDDO
    gibbs=amu*nel
    free_energy = gpband + gibbs
    IF (betael.GE.betaelmax) THEN
       mkt_entropy = 0._real_8
    ELSE
       mkt_entropy = 1._real_8/betael * spind * mkt_entropy
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE enert
  ! ==================================================================
  SUBROUTINE convfrie(nstate,nconv,nel,we,tislsd,amu,tconv)
    ! ==--------------------------------------------------------------==
    ! == TCONV=.TRUE. IF THE NUMBER OF STATES NSTATE IS SUFFICIENT    ==
    ! == I.E. THE LAST STATE COULD BE OCCUPANCY NUMBER < TOLLPAR      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nconv
    REAL(real_8)                             :: nel, we(nstate)
    LOGICAL                                  :: tislsd
    REAL(real_8)                             :: amu
    LOGICAL                                  :: tconv

    CHARACTER(*), PARAMETER                  :: procedureN = 'convfrie'
    REAL(KIND=real_8), PARAMETER             :: tollpar = 1.e-6_real_8

    INTEGER                                  :: i, ierr, ntotal
    REAL(real_8)                             :: wk(1)
    REAL(real_8), ALLOCATABLE                :: focc(:), weint(:)

    ALLOCATE(focc(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(weint(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (fint1%ttrot) THEN
       !$omp parallel do private(I)
       DO i=1,nstate
          weint(i)=-LOG(we(i))/fint1%betap
       ENDDO
    ELSE
       CALL dcopy(nstate,we(1),1,weint(1),1)
       CALL dscal(nstate,-1._real_8,weint(1),1)
    ENDIF
    wk(1)=1._real_8
    CALL adjmu(nstate,1,nel,fint1%betael,weint,wk,tislsd,amu,tconv)
    IF (.NOT.tconv) THEN
       RETURN
    ENDIF
    ! Check if the last state has focc<toll.
    CALL occup(weint,wk,amu,fint1%betael,focc,nstate,1,nel,tislsd)
    ntotal=0
    DO i=1,nstate
       IF (focc(i).GT.tollpar) THEN
          ntotal=ntotal+1
       ENDIF
    ENDDO
    tconv=(ntotal.LT.nconv)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(focc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(weint,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE convfrie
  ! ==================================================================
  FUNCTION dround(a,delta)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a, delta, dround

    dround=a-MOD(a,delta)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dround
  ! ==================================================================

END MODULE adjmu_utils
