MODULE bs_forces_diag_utils
  USE bsym,                            ONLY: bsclcs,&
                                             bsfac,&
                                             cnstwgt
  USE bsympnt,                         ONLY: fnbs
  USE cnstfc_utils,                    ONLY: restfc
  USE cotr,                            ONLY: cotr007,&
                                             gsrate,&
                                             resval
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE lsforce_utils,                   ONLY: lsforce
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE setbsstate_utils,                ONLY: setbsstate
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cnti,&
                                             fpar,&
                                             maxsys,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE tpar,                            ONLY: dt_ions
  USE updwf_utils,                     ONLY: updwf
  USE wrccfl_utils,                    ONLY: wrccfl
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bs_forces_diag

CONTAINS

  ! ==================================================================
  SUBROUTINE bs_forces_diag(nstate,c0,c2,cr,sc0,cscr,vpp,eigv,rhoe,&
       psi,tau0,velp,taui,fion,ifcalc,irec,tfor,tinfo)
    ! ==--------------------------------------------------------------==
    ! == ENTER WITH THE IONIC POSITIONS IN TAU0 AND AN INPUT GUESS    ==
    ! == FOR THE DENSITY IN RIN0, THIS ROUTINE RETURNS THE CONVERGED  ==
    ! == DENSITY, WAVEFUNCTIONS AND IONIC FORCES.                     ==
    ! ==--------------------------------------------------------------==
    ! == TAU0: ATOMIC POSITION                                        ==
    ! == VELP and TAUI are necessary for RESTART FILE                 ==
    ! == FION: IONIC FORCES                                           ==
    ! == IFCALC: total number of iterations                           ==
    ! ==--------------------------------------------------------------==
    ! == This subroutine is called in BS calculations instead of      ==
    ! == FORCES_DIAG. Much of that code is repeated here.             ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,bsfac), &
      c2(nkpt%ngwk,crge%n,bsfac), cr(*), sc0(nkpt%ngwk,crge%n,bsfac), cscr(*)
    REAL(real_8)                             :: vpp(*), eigv(crge%n,bsfac), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), &
                                                taui(:,:,:), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'bs_forces_diag'

    INTEGER                                  :: i, ierr, infr
    LOGICAL                                  :: cvbswf, cvhswf
    REAL(real_8)                             :: detot, etot0, etotbs, etoths, &
                                                tcpu, thl(2), time1, time2

    ALLOCATE(fnbs(3*maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    ! ==--------------------------------------------------------------==
    ener_com%ecnstr=0.0_real_8
    ener_com%erestr=0.0_real_8
    ropt_mod%convwf=.FALSE.
    cvbswf=.FALSE.
    cvhswf=.FALSE.
    etot0=0.0_real_8
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*) 'BSYMM: HS STATE'
    IF (.NOT.tinfo.AND.paral%io_parent)&
         WRITE(6,'(T9,A)')&
         'INFR         GEMAX    ETOT         DETOT     TCPU'
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    ! ==================================================================
    ! ==                        MAIN  LOOP                            ==
    ! ==================================================================
    DO infr=1,cnti%nomore_iter
       time1=m_walltime()
       ifcalc=ifcalc+1
       ! Diagonalization
       IF (.NOT.cvhswf) THEN
          ! HS state
          IF (infr.EQ.1) THEN
             bsclcs=2
             CALL setbsstate
          ENDIF
          CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,fion,cr,cscr,&
               vpp,eigv(1,2),rhoe,psi,nstate,.FALSE.,.TRUE.)
          cvhswf=ropt_mod%convwf
       ELSEIF (.NOT.cvbswf) THEN
          ! BS state
          IF (bsclcs.EQ.2) THEN
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) 'BSYMM: HS STATE CONVERGED. ',&
                     'SWITCHING TO BS STATE.'
                etot0=0.0_real_8
                IF (infi.EQ.1)THEN
                   CALL dcopy(2*nkpt%ngwk*crge%n,c0(1,1,2),1,c0,1)
                   CALL dcopy(crge%n,eigv(1,2),1,eigv,1)
                   CALL dcopy(2*nkpt%ngwk*crge%n,c2(1,1,2),1,c2,1)
                ENDIF
             ENDIF
             bsclcs=1
             CALL setbsstate
             CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,cr,cscr,vpp,eigv,rhoe,psi,&
                  nstate,.FALSE.,.TRUE.)
             IF (tinfo.AND.paral%parent) THEN
                CALL wrener
                IF (paral%io_parent)&
                     WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
                     '           ETOT        DETOT      TCPU'
             ENDIF
          ELSE
             CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,cr,cscr,vpp,eigv,rhoe,psi,&
                  nstate,.FALSE.,.TRUE.)
          ENDIF
          cvbswf=ropt_mod%convwf
       ENDIF
       ropt_mod%convwf=cvhswf.AND.cvbswf
       ! Printout the evolution of the iterative optimization
       IF (paral%parent) THEN
          detot=ener_com%etot-etot0
          IF (infr.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF (tinfo) THEN
             CALL wrprint_wfopt(eigv(1,bsclcs),crge%f,ener_com%amu,nstate,ener_com%etot,etot0,&
                  tcpu,gemax,cnorm,thl,.FALSE.,ifcalc,infr)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,&
                  '(T9,I4,T15,1PE12.3,T31,0PF12.6,T45,1PE12.3,T59,0PF7.2) ')&
                  infr,gemax,ener_com%etot,detot,tcpu
          ENDIF
          etot0=ener_com%etot
       ENDIF
       IF (MOD(ifcalc,store1%istore).EQ.0)CALL zhwwf(2,irec,c0,c2,nstate,eigv,&
            tau0,velp,taui,iteropt%nfi)
       ! Check to break the loop.
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) GOTO 200   ! Soft Exit (no force calculation)
       IF (ropt_mod%convwf) GOTO 100   ! Convergence of wavefunctions
    ENDDO
    ! ==================================================================
    ! ==                      END OF MAIN  LOOP                       ==
    ! ==================================================================
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
100 CONTINUE
    ! BS energy calculation 
    ! BS state
    bsclcs=1
    CALL setbsstate
    CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,cr,cscr,vpp,eigv,rhoe,psi,&
         nstate,tfor,.TRUE.)
    IF (paral%parent) THEN
       IF (infi.EQ.cnti%nomore) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)
             WRITE(6,*) 'BSYMM: BS STATE ENERGIES:'
          END IF
          CALL wrener
       ENDIF
       etotbs=ener_com%etot
       CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
    ENDIF
    ! HS state      
    bsclcs=2
    CALL setbsstate
    CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,fion,cr,cscr,vpp,&
         eigv(1,2),rhoe,psi,nstate,tfor,.TRUE.)
    IF (paral%parent) THEN
       IF (infi.EQ.cnti%nomore) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)
             WRITE(6,*) 'BSYMM: HS STATE ENERGIES:'
          END IF
          CALL wrener
       ENDIF
       etoths=ener_com%etot
       ! LS state      
       ener_com%etot = (1 + cnstwgt) * etotbs - cnstwgt * etoths
       CALL wrccfl(ener_com%etot,etotbs,etoths)
    ENDIF
    IF (tfor) THEN
       ! Forces from geometrical restraints
       IF (paral%parent) THEN
          CALL lsforce(fnbs,fion)
          DO i=1,cotr007%mrestr
             resval(i)=resval(i)+gsrate(i)*dt_ions
          ENDDO
          CALL restfc(tau0,fion)
       ENDIF
       CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
    ENDIF
200 CONTINUE
    bsclcs=1
    CALL setbsstate
    DEALLOCATE(fnbs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bs_forces_diag


END MODULE bs_forces_diag_utils
