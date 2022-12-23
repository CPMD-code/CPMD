MODULE finalp_utils
  USE cdft_utils,                      ONLY: cdft_finalize
  USE cnst,                            ONLY: au_deb,&
                                             au_kb
  USE cnstpr_utils,                    ONLY: cnstpr
  USE dipomod,                         ONLY: moment
  USE dum2_utils,                      ONLY: dumpr
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE metr,                            ONLY: metr_com
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE parac,                           ONLY: paral
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: infi
  USE store_types,                     ONLY: cprint,&
                                             iprint_force,&
                                             rout1
  USE strs,                            ONLY: paiu
  USE struc_utils,                     ONLY: pstruc
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wrener_utils,                    ONLY: wreigen,&
                                             wrener,&
                                             wrstress
  USE wrgeo_utils,                     ONLY: wrgeof,&
                                             wrgeox

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finalp

CONTAINS

  ! ==================================================================
  SUBROUTINE finalp(tau0,fion,velp,eigv)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                velp(:,:,:), &
                                                eigv(crge%n,nkpt%nkpts)

    CHARACTER(*), PARAMETER                  :: procedureN = 'finalp'

    INTEGER                                  :: i, ia, is, isub, j, k, l, nxst
    REAL(real_8)                             :: d1, d2, d3, dd, fact, out(3,3)

    CALL tiset(procedureN,isub)
    IF (paral%io_parent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,'(" *",62X,"*")')
       WRITE(6,'(" *",23X,A,24X,"*")') ' FINAL RESULTS '
       WRITE(6,'(" *",62X,"*")')
       WRITE(6,'(1X,64("*"))')
    ENDIF
    CALL wrgeof(tau0,fion)
    IF (rout1%xgout) THEN
       IF (MOD(infi,cnti%ngxyz).NE.0) THEN
          nxst=cnti%ngxyz
          cnti%ngxyz=1
          CALL wrgeox(tau0)
          cnti%ngxyz=nxst
       ENDIF
    ENDIF
    CALL dumpr
    CALL cnstpr
    IF (cntl%tdiag) CALL wreigen(eigv,crge%f,ener_com%amu,crge%n)
    IF (paral%io_parent) WRITE(6,*)
    ! ==--------------------------------------------------------------==
    IF (cntl%tprcp) THEN
       IF (paral%io_parent) THEN
          WRITE(6,*) '  CELL PARAMETERS:'
          DO i=1,3
             WRITE(6,'(T21,3(F15.8))') (metr_com%ht(i,j),j=1,3)
          ENDDO
          WRITE(6,'(/,1X,64("*"))')
       ENDIF
       CALL dcopy(9,paiu,1,out,1)
       DO is=1,ions1%nsp
          fact=rmass%pma(is)
          DO ia=1,ions0%na(is)
             DO k=1,3
                DO l=1,3
                   out(k,l)=out(k,l)+fact*velp(k,ia,is)*velp(l,ia,is)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF (paral%io_parent) THEN
          WRITE(6,*) '  TOTAL STRESS (kB)'
          DO i=1,3
             WRITE(6,'(7X,3(F18.8))') ((out(i,j)/parm%omega)*au_kb,j=1,3)
          ENDDO
       ENDIF
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("*"),/)')
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND..NOT.tkpts%tknoswap) THEN
       WRITE(6,'(A)') ' ELECTRONIC GRADIENT:'
       WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
            gemax,'         NORM =',cnorM
       IF (cprint%iprint(iprint_force).GT.0) THEN
          WRITE(6,'(A)') ' NUCLEAR GRADIENT:'
          WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
               gnmax,'         NORM =',gnorM
          WRITE(6,*)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (tkpts%tknoswap) THEN
       ! Do not display energies.
    ELSE
       ! Write energies
       CALL wrener
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND.cntl%caldip) THEN
       WRITE(6,'(A)') ' DIPOLE MOMENT'
       d1=moment%dmom(1)
       d2=moment%dmom(2)
       d3=moment%dmom(3)
       dd=SQRT(d1*d1+d2*d2+d3*d3)
       WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'   
       WRITE(6,'(4F12.5,A)') d1,d2,d3,dd,'   atomic units'
       WRITE(6,'(4F12.5,A)') au_deb*d1, au_deb*d2, au_deb*d3,&
            au_deb*dd,'   Debye'
       WRITE(6,'(A,E16.8)') ' ChkSum(DIP) = ',SUM(ABS(moment%dmom(1:3)))
       WRITE(6,*)
    ENDIF
    IF (cntl%tpres) CALL wrstress
    CALL pstruc(tau0)
    IF (paral%io_parent) WRITE(6,'(1X,64("*"),/)')
    IF (cntl%cdft)CALL cdft_finalize
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE finalp
  ! ==================================================================

END MODULE finalp_utils
