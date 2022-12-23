MODULE k_diis_rhofix_utils
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE k_forces_driver,                 ONLY: k_forces
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: kdrho,&
                                             kdrhon
  USE machine,                         ONLY: m_walltime
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: ropt_mod
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: cprint
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE wrener_utils,                    ONLY: wrprint_wfopt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k_diis_rhofix

CONTAINS

  ! ==================================================================
  SUBROUTINE k_diis_rhofix(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
       rhoe,psi,nstate,tfor,tinfo,ddrho)
    ! ==--------------------------------------------------------------==
    ! ==               UPDATES THE WAVEFUNCTIONS                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    COMPLEX(real_8)                          :: pme(*), gde(*)
    REAL(real_8)                             :: vpp(:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,nkpt%nkpts)
    COMPLEX(real_8) :: sc0(nkpt%ngwk,nstate), c2(nkpt%ngwk,nstate), &
      c0(nkpt%ngwk,nstate,nkpt%nkpnt)
    LOGICAL                                  :: tfor, tinfo
    REAL(real_8)                             :: ddrho

    INTEGER                                  :: i, jnfi, numldiis
    REAL(real_8)                             :: etot0, &
                                                f_save(nstate*nkpt%nkpnt), &
                                                ggtol_rho, tcpu, thl(2), &
                                                time1, time2

!nnr1)
! Variables
! ==--------------------------------------------------------------==

    IF (cntl%nonort .OR. cntl%tpath) THEN
       CALL stopgm('K_RHOFIX','CNTL%NONORT,TPATH NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==  UPDATE THE WAVEFUNCTIONS                                    ==
    ! ==--------------------------------------------------------------==
    ! copy occupation numbers
    DO i=1,nstate*nkpt%nkpnt
       f_save(i)=crge%f(i,1)
       crge%f(i,1)=1._real_8
    ENDDO
    numldiis = cnti%maxldiis
    etot0 = 0.0_real_8
    cntl%tdiag = .FALSE.
    ggtol_rho=cntr%tolrhofix
    IF (ddrho.GT.10.0_real_8*cntr%tolog) THEN
       ggtol_rho=5.e+2_real_8*cntr%tolrhofix
    ELSEIF (ddrho.GT.5.0_real_8*cntr%tolog) THEN
       ggtol_rho=5.e+1_real_8*cntr%tolrhofix
    ELSEIF (ddrho.GT.2.0_real_8*cntr%tolog) THEN
       ggtol_rho=5._real_8*cntr%tolrhofix
    ENDIF
    IF (paral%io_parent.AND.tinfo) WRITE(6,'(A,1PE14.6,A,1PE14.6)')&
         'DDRHO ',ddrho, '  TOL ',ggtol_rho

    DO jnfi=1,numldiis
       time1=m_walltime()
       CALL k_forces(c0,c2,sc0,tau0,pme,gde,vpp,eigv,fion,rhoe,psi,&
            f_save,nstate,nkpt%nkpnt,&
            .TRUE.,tfor,kdrho,kdrhon, .FALSE.,jnfi)
       IF (.NOT.tfor) CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (gemax.LT.ggtol_rho.OR. soft_com%exsoft) ropt_mod%convwf=.TRUE.
       IF (jnfi.LT.cnti%minldiis)ropt_mod%convwf=.FALSE.
       IF (paral%parent) THEN
          ropt_mod%engpri=MOD(jnfi-1,cprint%iprint_step).EQ.0
       ELSE
          ropt_mod%engpri=.FALSE.
       ENDIF
       ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%parent) THEN
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
               gemax,cnorm,thl,ropt_mod%convwf,jnfi,jnfi)
       ENDIF
       IF (ropt_mod%convwf) GOTO 100
    ENDDO
    IF (paral%io_parent.AND.tinfo) THEN
       WRITE(6,'(1X,64("!"))')
       WRITE(6,'(" !!",A,T64,"!!")')&
            ' RWFOPT| THE MAXIMUM NUMBER OF STEP IS REACHED'
       WRITE(6,'(1X,64("!"))')
    ENDIF
    ropt_mod%convwf=.TRUE.
100 CONTINUE

    CALL k_forces(c0,c2,sc0,tau0,pme,gde,vpp,eigv,fion,rhoe,psi,&
         f_save,nstate,nkpt%nkpnt,&
         .TRUE.,tfor,kdrho,kdrhon,.FALSE.,jnfi)
    !$omp parallel do private(i)
    DO i=1,nstate*nkpt%nkpnt
       crge%f(i,1)=f_save(i)
    ENDDO

    cntl%tdiag = .TRUE.
    ropt_mod%convwf=.FALSE.
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_diis_rhofix

END MODULE k_diis_rhofix_utils
