MODULE rwfopt_nuc_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf,&
                                             textfld
  USE elct,                            ONLY: crge
  USE enbandpri_utils,                 ONLY: enbandpri
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE finalp_utils,                    ONLY: finalp
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: lqmmm
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: rhoo
  USE prmem_utils,                     ONLY: prmem
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rscpot_utils,                    ONLY: rscpot
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE store_types,                     ONLY: cprint,&
                                             restart1,&
                                             rout1,&
                                             store1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             fpar,&
                                             maxsys,&
                                             nacc,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE updrho_utils,                    ONLY: give_scr_updrho,&
                                             updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rwfopt_nuc
  PUBLIC :: give_scr_rwfopt_dum
  !public :: lseprof_dum

CONTAINS

  ! ==================================================================
  SUBROUTINE rwfopt_nuc(c0,c2,sc0,pme,gde,vpp,eigv,nstate,psi,&
       rhoe)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpnt), c2(nkpt%ngwk,crge%n), &
      sc0(nkpt%ngwk,crge%n), pme(*), gde(*)
    REAL(real_8)                             :: vpp(:), &
                                                eigv(crge%n,nkpt%nkpts)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rwfopt_nuc'

    INTEGER                                  :: i, ierr, irec(100), nhpsi, nnx
    LOGICAL                                  :: status, statusdummy, &
                                                update_pot
    REAL(real_8)                             :: ekin1, ekin2, ekincp, ekinh1, &
                                                ekinh2, etot0, tcpu, temp1, &
                                                temp2, thl(2), time1, time2

! Variables
! real(8) :: smat(nstate,nstate)
! pointer (ip_smat,smat)
! QMMM
! ==--------------------------------------------------------------==

    time1 =m_walltime()
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm)THEN
       IF (textfld)THEN
          ! CALL MEMORY(IP_EXTF,KR1*KR2S*KR3S,'EXTF')
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
    ENDIF
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    CALL mm_dim(mm_go_qm,statusdummy)
    nacc = 7
    IF (cntl%tdiag) THEN
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rout0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       rhoo => rin0

    ENDIF
    iteropt%nfi = 0
    ener_com%ecnstr=0.0_real_8
    ener_com%erestr=0.0_real_8
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%convwf=.FALSE.
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    ! CALL INITRUN(IREC,C0,C2,SC0,RHOE,PSI,EIGV)
    ! CALL WRITE_IREC(IREC)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    ! ETOT0=0._real_8
    ! EBOGO=0._real_8
    ! IF(PARENT) THEN
    ! TIME2 = TIMEF()
    ! TCPU = (TIME2 - TIME1)*0.001
    ! WRITE(6,'(A,T50,F8.2,A8)')
    ! &     ' TIME FOR WAVEFUNCTION INITIALIZATION:',TCPU,' SECONDS'
    ! CALL PRMEM('    RWFOPT')
    ! ENDIF
    ! ==--------------------------------------------------------------==
    ! ==      THE BASIC LOOP FOR WAVEFUNCTION OPTIMIZATION            ==
    ! ==--------------------------------------------------------------==
    ! ME... IF statement introduced by M. Eichinger for interface mode
    IF (cntl%tinter.AND.restart1%restart.AND.restart1%rwf) THEN
       ! do nothing
    ELSE
       update_pot=.TRUE.
       DO infi=1,cnti%nomore
          time1=m_walltime()
          iteropt%nfi=iteropt%nfi+1
          IF (infi.GT.1) update_pot=.FALSE.
          IF (cntl%tdiag) THEN
             CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                  rhoe,psi,&
                  crge%n,.FALSE.,.TRUE.,.FALSE.,infi,thl,nhpsi)
          ELSE
             ! UPDATE THE WAVEFUNCTIONS

             ! CALL DSYRK('U','T',NSTATE,2*NGW,2._real_8,C0,2*NGW,0._real_8,SMAT,NSTATE)

             CALL updwf(c0(:,:,1),c2(:,:),sc0,tau0,fion,pme,gde,vpp,eigv,&
                  rhoe,psi,crge%n,.FALSE.,update_pot)

             IF (infi.EQ.cnti%nomore-2) THEN
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(6,*) 'something is weird ...'
                   IF (paral%io_parent)&
                        WRITE(6,*) 'infi =',infi,'   and nomore-2=',cnti%nomore-2
                   CALL stopgm('RWFOPT_nuc','PROGRAM STOP',& 
                        __LINE__,__FILE__)
                ENDIF
             ENDIF

             ! CALL DSYRK('U','T',NSTATE,2*NGW,2._real_8,C0,2*NGW,0._real_8,SMAT,NSTATE)

             DO i = 1,nstate
                ! if (parent) WRITE(6,*) i,i, smat(i,i)
                CALL m_flush(6)
             ENDDO

          ENDIF
          IF (paral%parent) THEN
             ropt_mod%engpri=MOD(infi-1,cprint%iprint_step).EQ.0
          ELSE
             ropt_mod%engpri=.FALSE.
          ENDIF
          ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
          IF (paral%parent) THEN
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
                  gemax,cnorm,thl,ropt_mod%convwf,iteropt%nfi,infi)
          ENDIF
          ! PERIODICALLY WRITE THE RESTART FILE
          IF (MOD(infi,store1%istore).EQ.0.OR.infi.EQ.cnti%nomore.OR.ropt_mod%convwf)&
               CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
          IF (ropt_mod%convwf) GOTO 100
       ENDDO
    ENDIF
    IF ((.NOT.ropt_mod%convwf).AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Probably your stepsize is too small...'
       CALL stopgm('RWFOPT_nuc','WAVEFUNCTIONS NOT CONVERGED',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==     END OF MAIN LOOP                                         ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.gemax.GT.cntr%tolog.AND.&
         .NOT.cntl%ksener.AND..NOT.tkpts%tknoswap) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,T64,"!!")')&
            ' RWFOPT| THE MAXIMUM NUMBER OF STEP IS REACHED'
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,F10.6,A,T64,"!!")')&
            '         BUT NO CONVERGENCE (DRHOMAX=',gemax,')'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
    IF (tkpts%tknoswap) THEN
       ! No calculation of ionic forces.
       ropt_mod%convwf=.FALSE.
    ENDIF
100 CONTINUE
    ! CALL SUBROUTINE TO FIND OPTIMALDUMMY

    IF (rout1%rhoout) CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
    IF (rout1%teband) CALL enbandpri(eigv,crge%f,ener_com%amu,crge%n,nkpt%nkpts)
    IF (tkpts%tonlydiag) GOTO 150
    ! Calculate ionic forces
    IF (ropt_mod%convwf) THEN
       IF (cntl%tdiag) THEN
          CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               crge%n,.TRUE.,.FALSE.,cntl%tpres,infi,thl,nhpsi)
       ELSE
          ropt_mod%calste=cntl%tpres
          CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
               rhoe,psi,crge%n,.TRUE.,update_pot)
       ENDIF
       ! Calculate norm of nuclear gradient
       CALL gsize(fion,gnmax,gnorm)
       IF (cntl%tpres) CALL totstr
    ELSE
       gnmax=0.0_real_8
       gnorm=0.0_real_8
    ENDIF
150 CONTINUE
    ! CALL SUBROUTINE TO FIND OPTIMALDUMMY
    ! call BEST_SGDUMMY(C0,RHOE,PSI,SCR,LSCR,NSTATE,TAU0,FION)
    ! 
    IF (lspin2%teprof.AND.lspin2%tlse) THEN
       CALL lseprof_dum(c0,rhoe,psi,tau0,fion,crge%n)
    ENDIF
    ! 
    IF (paral%parent) THEN
       CALL prmem('    RWFOPT')
       CALL finalp(tau0,fion,tau0,eigv)
       CALL geofile(tau0,fion,'WRITE')
    ENDIF
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tdiag) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,status)
    RETURN
  END SUBROUTINE rwfopt_nuc
  ! ==================================================================
  SUBROUTINE give_scr_rwfopt_dum(lrwfopt,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrwfopt
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: linitrun, lprof, lrhopri, &
                                                lupd, nstate

    CALL give_scr_initrun(linitrun,tag)
    nstate=crge%n
    lupd=0
    lrhopri=0
    IF (cntl%tdiag) THEN
       CALL give_scr_updrho(lupd,tag,nstate,.TRUE.,cntl%tpres)
    ELSE
       CALL give_scr_updwf(lupd,tag,nstate,.TRUE.)
    ENDIF
    lprof=0
    IF (lspin2%tlse.AND.lspin2%teprof) lprof=4*ncpw%ngw+100
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    lrwfopt=MAX(linitrun,lupd,lrhopri)+lprof
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rwfopt_dum

END MODULE rwfopt_nuc_utils

! ==================================================================
SUBROUTINE lseprof_dum(c0,rhoe,psi,tau0,fion,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE fft_maxfft,                      ONLY: maxfftn
  USE system , ONLY:maxsys,ncpw,fpar
  USE parac, ONLY : paral,parai
  USE spin , ONLY:clsd,lspin2,lspin3
  USE ener , ONLY:ener_com
  USE cnst , ONLY:fpi
  USE rscpot_utils, ONLY : rscpot
  IMPLICIT NONE
  REAL(real_8)                               :: rhoe(fpar%nnr1,clsd%nlsd)
  COMPLEX(real_8)                            :: psi(maxfftn,clsd%nlsd)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx),&
                                                fion(3,maxsys%nax,maxsys%nsx)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'lseprof_dum'

  COMPLEX(real_8), ALLOCATABLE               :: ca(:), cb(:), sca(:)
  INTEGER                                    :: ia, ierr
  LOGICAL                                    :: tccc
  REAL(real_8)                               :: ang, angd, cosa, deee, ezero, &
                                                sina, tpi

  IF (lspin2%tcas22) THEN
     lspin2%troks=.TRUE.
     lspin2%tcas22=.FALSE.
     tccc=.TRUE.
  ELSE
     tccc=.FALSE.
  ENDIF
  tpi=0.5_real_8*fpi

  ALLOCATE(ca(ncpw%ngw),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(cb(ncpw%ngw),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)

  CALL dcopy(2*ncpw%ngw,c0(1,clsd%ialpha),1,ca,1)
  CALL dcopy(2*ncpw%ngw,c0(1,clsd%ibeta),1,cb,1)
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(/,1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,*)
     IF (paral%io_parent)&
          WRITE(6,'(T20,A,T60,A)') '  Rotation Angle ','Energy'
  ENDIF
  lspin3%rotab=lspin3%rotab/tpi*360._real_8
  DO ia=0,9
     ang=ia*5._real_8
     angd=ang/360._real_8*tpi
     cosa=COS(angd)
     sina=SIN(angd)
     CALL dcopy(2*ncpw%ngw,ca,1,c0(1,clsd%ialpha),1)
     CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ialpha),1)
     CALL daxpy(2*ncpw%ngw,sina,cb,1,c0(1,clsd%ialpha),1)
     CALL dcopy(2*ncpw%ngw,cb,1,c0(1,clsd%ibeta),1)
     CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ibeta),1)
     CALL daxpy(2*ncpw%ngw,-sina,ca,1,c0(1,clsd%ibeta),1)

     CALL rscpot(c0,tau0,fion,rhoe,psi,&
          .FALSE.,.FALSE.,nstate,1)
     IF (ia.EQ.0) ezero=ener_com%etot
     deee=ener_com%etot-ezero
     IF (paral%io_parent)&
          WRITE(6,'(T28,F8.0,T50,F16.8)') lspin3%rotab+ang,deee
  ENDDO
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-B ",lspin3%hablse
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-O ",lspin3%haolse
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : B-O ",lspin3%hbolse
     IF (paral%io_parent)&
          WRITE(6,'(1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,*)
  ENDIF

  DEALLOCATE(ca,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(cb,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(sca,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  IF (tccc) THEN
     lspin2%tcas22=.TRUE.
     lspin2%troks=.FALSE.
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE lseprof_dum
! ==================================================================
