MODULE hubbardu_utils

USE atom,                           ONLY: ecpfiles
USE atwf,                           ONLY: atwf_mod,&
                                          atwp,&
                                          m1shlx
USE cnst,                           ONLY: fpi
USE cppt,                           ONLY: gk,&
                                          twnl,&
                                          hg,inyh
USE cvan,                           ONLy: qq
USE hubbardu
USE domdr_utils,                    ONLY: domdr
USE elct,                           ONLY: crge
USE error_handling,                 ONLY: stopgm
USE fileopenmod,                    ONLY: fo_app
USE fileopen_utils,                 ONLY: fileclose,&
                                          fileopen
USE fitpack_utils,                  ONLY: curv1,curv2
USE fnlalloc_utils,                 ONLY: fnlalloc,&
                                          fnldealloc,&
                                          fnl_set
USE geq0mod,                        ONLY: geq0
USE gvec,                           ONLY: gvec_com
USE inscan_utils,                   ONLY: inscan
USE ions,                           ONLY: ions0,&
                                          ions1
USE kinds,                          ONLY: real_8
USE kpts,                           ONLY: tkpts
USE lsfbtr_utils,                   ONLY: lsfbtr
USE mm_dimmod,                      ONLY: nat_grm,&
                                          cpsp
USE mp_interface,                   ONLY: mp_bcast,&
                                          mp_sum
USE nlps,                           ONLY: nghtol,&
                                          nlps_com
USE ortho_utils,                    ONLY: give_scr_ortho,&
                                          ortho
USE ovlap_utils,                    ONLY: ovlap
USE  parac,                         ONLY: paral,&
                                          parai
USE phfac_utils,                    ONLY: phfac
USE pslo,                           ONLY: pslo_com
USE qspl,                           ONLY: ggng,&
                                          ggnh,&
                                          nsplpo
USE readsr_utils,                   ONLY: xstring
USE recpnew_utils,                  ONLY: ckgrid,tgrid
USE rnlsm_utils,                    ONLY: give_scr_rnlsm,rnlsm
USE ropt,                           ONLY: iteropt
USE sfac,                           ONLY: dfnl,&
                                          eigr,&
                                          ei1,ei2,ei3,eigrb,&
                                          fnl
USE spin,                           ONLY: spin_mod
USE spsi_utils,                     ONLY: spsi,&
                                          give_scr_spsi
USE sphe,                           ONLY: gcutka
USE system,                         ONLY: cntl,&
                                          iatpt,&
                                          maxsys,&
                                          ncpw,&
                                          parm,parap
USE timer,                          ONLY: tihalt,&
                                          tiset

USE vdbp,                           ONLY: ncpr1,&
                                          vdb_r,&
                                          vdb_pawf 
USE ylmr_utils,                     ONLY: ylmr
USE zeroing_utils,                  ONLY: zeroing

    !     Written by nisanth.nair@theochem.rub.de May/2008
    !     DONE NATTOT -> hubbu%nuproj
    !     TODO NSTATE -> PARALLELIZE OVER PROCESSORS
    !     TODO HUBINTS loop over ions1%nat -> loop over IUATM
    !     TODO C2U summing to C2 done with DAXPY and a subroutine inside this module
    !     TODO VHUB -> reformulate 
    !     TODO better organization of HUBINTS_USPP2
    !
    !     Optimized for large number of U atoms and transfered to CPMD v4 
    !       by niklas.siemer@theochem.rub.de 
    !
    !
    IMPLICIT NONE
   
    LOGICAL, allocatable, save       ::  is_u_sp(:),is_u_at(:)
    REAL(real_8), allocatable,save   ::  C0XCA(:,:)
    PRIVATE                          ::  dftuinit,occmat 
    PUBLIC                           ::  hubbardUcorrection,&
                                         give_scr_hubbardu,&
                                         add_hubbardu

    CONTAINS
    !     ====================================================================
          SUBROUTINE hubbardUcorrection(c0,c2u,tau0,fion,nstate,psi,tfor,ispin0)
    !     == -------------------------------------------------------------- ==
    !     == This routine calcualtes contributions from DFT+U part          ==
    !     == -------------------------------------------------------------- ==
          IMPLICIT NONE
    !
          INTEGER, INTENT(IN)            ::  nstate,ispin0
          COMPLEX(real_8), INTENT(IN)    ::  c0(:,:),psi(:),c2u(:,:)

          REAL(real_8) , INTENT(INOUT)   ::  tau0(:,:,:),fion(:,:,:)
          LOGICAL ,INTENT(IN)            ::  tfor

          CHARACTER(*), PARAMETER        ::  procedureN = 'hubbardUcorrection'

          LOGICAL                        ::  tivan_bak
    !
          INTEGER                        ::  isub,ierr
    !
          IF (hubbu%debug) THEN 
            IF (paral%io_parent) WRITE(6,*) procedureN ,"| Start Hubbard U correction"
          ENDIF
          CALL tiset(procedureN,isub)
    ! Initialize
          IF(hubbu%firstcall.EQ.0)CALL dftuinit(nstate)
          hubbu%firstcall=hubbu%firstcall+1
    !
          IF(hubbu%nuproj.EQ.0)CALL STOPGM(procedureN,&
                'hubbu%nuproj CANNOT BE ZERO', __LINE__,__FILE__)
    !
    ! Memory allocations
          allocate(mycatom(ncpw%ngw,hubbu%nuproj),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
          allocate(c0xca(nstate,hubbu%nuproj),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
          IF(pslo_com%tivan) allocate(mycatom0(ncpw%ngw,hubbu%nuproj)) 
    !
    ! Testing the ISPIN variable 
          IF((ISPIN0.NE.0).AND.TFOR) CALL STOPGM(procedureN,&
                'ISPIN.NE.0 .AND. TFOR', __LINE__,__FILE__)
          IF((ISPIN0.GT.2).OR.(ISPIN0.LT.0))CALL STOPGM(procedureN,&
                'ISPIN0 HAS WRONG VALUE', __LINE__,__FILE__)
    !
    !numerical derivative
    !      CALL DFTPLUSU_NUMDER(C0,TAU0,FION,NSTATE,PSI,TFOR)
    !      STOP
    !
    ! Calculate the occupation matrix OM
          CALL OCCMAT(C0,TAU0,NSTATE,PSI,TFOR,ISPIN0)
    !
    ! Estimate energy contribution from Habbard term
          CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,TFOR,ISPIN0)
    !
          IF(pslo_com%tivan)CALL rnlsm(c0,nstate,1,1,tfor)
    !
!         CALL MEMORY_CHECK
    !
          deallocate(mycatom,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
          deallocate(c0xca,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
          IF(pslo_com%tivan)deallocate(mycatom0)
    !
          CALL TIHALT(procedureN,ISUB)
          RETURN      
          END SUBROUTINE hubbardUcorrection
    !     ====================================================================
          SUBROUTINE dftuinit(nstate)
    !     == -------------------------------------------------------------- ==
    !     == This routine initializes every thing for DFT+U                 ==
    !     == -------------------------------------------------------------- ==
    !
    !     Identify states corresponds to m values of Habbard atoms 
    !
          INTEGER, INTENT(IN)       ::  NSTATE
          CHARACTER(*), PARAMETER   ::  procedureN = 'dftuinit'
    !
          INTEGER                   ::  iorb_min,iorb_max,is,iat,jat,ish,l,&
                                        lmin,lmax,iuatm,ihubl
          LOGICAL                   ::  found_ls,read_proj_log(maxsys%nsx)
          INTEGER                   ::  isub, ierr
    !
          CALL tiset(procedureN,isub)
    ! 
    ! Allocate OM array
          allocate(om(2,hubbu%nuatm,7,7),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(fomi(2,hubbu%nuatm,7,7,3),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(is_u_sp(ions1%nsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(is_u_at(ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(fion_om(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(atwfr_u(maxsys%mmaxx,m1shlx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
          allocate(atrg_u(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                 __LINE__,__FILE__)
    !
          JAT=0
          IORB_MIN=0
          IORB_MAX=0
          hubbu%nuproj=0
          FOUND_LS=.FALSE.
          DO IS=1,ions1%nsp
            IS_U_SP(IS)=.FALSE.
            DO IAT=1,ions0%na(IS)
              JAT=JAT+1
              IS_U_AT(JAT)=.FALSE.
              DO ISH=1,atwf_mod%nshell(IS)
                 L=atwf_mod%lshell(ISH,IS)
                 LMIN=1
                 LMAX=2*L+1
                 IORB_MIN=IORB_MAX+LMIN
                 IORB_MAX=IORB_MAX+LMAX
                 DO IUATM=1,hubbu%nuatm
                    IF(hubbu%uatm(IUATM).EQ.JAT)THEN
                       IS_U_AT(JAT)=.TRUE.
                       IS_U_SP(IS)=.TRUE.
                       DO IHUBL=1,hubbu%nl(IUATM)
                          IF(ISH.EQ.hubbu%s(IUATM,IHUBL).AND.&
                             L.EQ.hubbu%l(IUATM,IHUBL))THEN
                             FOUND_LS=.TRUE.
                             hubbu%muatm(1,1,IUATM)=IORB_MIN
                             hubbu%muatm(2,1,IUATM)=IORB_MIN
                             hubbu%muatm(1,2,IUATM)=IORB_MAX
                             hubbu%muatm(2,2,IUATM)=IORB_MAX
                             hubbu%nuproj=hubbu%nuproj+IORB_MAX-IORB_MIN+1
                          END IF
                       END DO
                    END IF
                 END DO
              END DO
            END DO
          END DO
          IF(.NOT.FOUND_LS) &
            CALL STOPGM(procedureN,'ORBITALS FOR PROJECTION NOT FOUND!', __LINE__,&
             __FILE__)
    !dbg
    !
    !  Load the corresponding the atomic orbitals (projectors)
          DO IUATM=1,hubbu%nuatm
            CALL LOADUPROJ(hubbu%uatm(IUATM),READ_PROJ_LOG)
          END DO
          IF(paral%io_parent)THEN
            WRITE(6,*)'== DFT+U PROJECTORS =='
            WRITE(6,'(A,I8)')' TOTAL NUMBER OF PROJECTORS =',hubbu%nuproj
            DO IUATM=1,hubbu%nuatm
              WRITE(6,'(3I10)')hubbu%uatm(IUATM),hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
            END DO
            WRITE(6,*)'======================'
          ENDIF
          CALL TIHALT(procedureN,ISUB)

          END SUBROUTINE dftuinit
!     ====================================================================
      SUBROUTINE occmat(c0,tau0,nstate,psi,tfor,ispin0)

!     == -------------------------------------------------------------- ==
!     == This routine calcualtes the occupation matrix needed for DFT+U ==
!     == Output occupation matrix OM() is in dftu.inc                   == 
!     == Orthogonalized atomic wavefunctions are used for projections   == 
!     == Nuclear forces due to occupation matrix is also calculated     ==
!     == Ouput occupation matrix forces due to ions in FOMI() (dftu.inc)==
!     == -------------------------------------------------------------- ==
      IMPLICIT NONE
      INTEGER, INTENT(IN)           ::   NSTATE,ISPIN0
      COMPLEX(real_8)               ::   C0(:,:),PSI(:)
      REAL(real_8)                  ::   TAU0(:,:,:)
      LOGICAL                       ::   TFOR
      CHARACTER(len=30)             ::   TAG
      CHARACTER(*), PARAMETER       ::   procedureN = 'occmat'
    
      INTEGER    ISUB,L_RNLSM,L_ORTHO,NOMAX,L_OCCMAT,ISPIN,ISTATE,&
                IUATM,M1,M2,M10,M20,M0,MOFF,MM1,MM2
!
      REAL(real_8)                  ::   FFI,FQQ(2)
      INTEGER                       ::   FIRSTCALL,ierr
      DATA       FIRSTCALL / 0 /
      SAVE       FIRSTCALL
!
      COMPLEX(real_8), ALLOCATABLE  :: myc0(:,:)
!
      LOGICAL                       ::  TLSD_BAK
      INTEGER                       ::  ISPIN_MIN,ISPIN_MAX
!
      CALL TISET(procedureN,ISUB)
      NOMAX=MAX(NSTATE,hubbu%nuproj)
!
!  New phase factors
      CALL PHFAC(TAU0)
!
!  Scratch for orthogonalization of C0
      L_ORTHO=1
      L_RNLSM=1
      CALL GIVE_SCR_ORTHO(L_ORTHO,TAG,NOMAX)
      CALL GIVE_SCR_RNLSM(L_RNLSM,TAG,NOMAX,.FALSE.)
      L_OCCMAT=MAX(L_ORTHO,L_RNLSM)
!
      allocate(myc0(ncpw%ngw,nstate),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

!  C0 will be changed in this routine; backing up C0
      CALL DCOPY(2*ncpw%ngw*NSTATE,C0,1,MYC0,1)
!
!  Projections on beta-functions are stored in FNL
      IF(pslo_com%tivan)THEN
        CALL rnlsm(myc0,nstate,1,1,.false.)
!       |S|psi>
        CALL spsi(nstate,myc0)
      END IF
!
! Load atomic orbitals and store it for further use
! and Norm/Orthogoanlize if necessary
      CALL ORTHOCATOM(TAU0,NSTATE,PSI,TFOR)
!
!     <phi|S|psi>
      CALL OVLAP2(ncpw%ngw,NSTATE,hubbu%nuproj,C0XCA,MYC0,MYCATOM,.True.)
!      CALL GLOSUM(NSTATE*hubbu%nuproj,C0XCA)
      CALL mp_sum(c0xca,nstate*hubbu%nuproj,parai%allgrp)
!
! Occupation number is calculated 
!     n(m1,m2,s) = sum_i(s) f_i *<phi(m1)|S|psi(i,s)>*<phi(m2)|S|psi(i,s)>
!     CALL zeroing(OM,2*hubbu%nuatm*7*7)
      CALL zeroing(om)
      MOFF=0
      DO IUATM=1,hubbu%nuatm
         M0=hubbu%muatm(1,1,IUATM)-1
         DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM) 
            M10=M1-M0
            MM1=MOFF+M10
            DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM) 
               M20=M2-M0
               MM2=MOFF+M20
               CALL zeroing(FQQ)
               DO ISTATE=1,NSTATE
                  FFI=crge%f(istate,1)
                  IF(ISPIN0.GT.0)THEN
!                    FFI=1.D0
                    ISPIN=ISPIN0
                  ELSE
!                    FFI=1.D0
                    ISPIN=1
                    IF(cntl%tlsd.AND.(ISTATE.GT.spin_mod%nsup))ISPIN=2
                  END IF
                  FQQ(ISPIN)=FQQ(ISPIN)+FFI*C0XCA(ISTATE,MM1)*C0XCA(ISTATE,MM2)
               END DO
               OM(1,IUATM,M10,M20)=FQQ(1)
               OM(2,IUATM,M10,M20)=FQQ(2)
            END DO
         END DO
         MOFF=MOFF+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
      END DO
      ISPIN_MIN=1
      ISPIN_MAX=1
      IF(cntl%tlsd)ISPIN_MAX=2
      IF(ISPIN0.GT.0)THEN
        ISPIN_MIN=ISPIN0
        ISPIN_MAX=ISPIN0
      END IF
      DO ISPIN=ISPIN_MIN,ISPIN_MAX
        DO IUATM=1,hubbu%nuatm
          M0=hubbu%muatm(1,1,IUATM)-1
          DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM) 
             M10=M1-M0
             DO M2=M1+1,hubbu%muatm(1,2,IUATM) 
                M20=M2-M0
!                IF(DABS(OM(ISPIN,IUATM,M20,M10)).LT.1.D-5)
!     &                  OM(ISPIN,IUATM,M20,M10)=0.D0
                OM(ISPIN,IUATM,M10,M20)=OM(ISPIN,IUATM,M20,M10)
             END DO
          END DO
        END DO
      END DO
!
! Calculate forces on ions due to occupation matrix
      IF(TFOR)CALL DOMDR(C0,MYCATOM,c0xca,TAU0,NSTATE,PSI,.TRUE.)
!
! Print the occupation matrix
      IF(paral%io_parent)CALL PRINT_OCCMAT(hubbu%tpom,cntl%tlsd,ISPIN0)
!
      FIRSTCALL=FIRSTCALL+1
!
      deallocate(myc0,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
!
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE occmat
!
! ====================================================================
      SUBROUTINE ORTHOCATOM(TAU0,NSTATE,PSI,TFOR)
!     == -------------------------------------------------------------- ==
!     == This routine load and orthogonalizes atomic orbtials CATOM     ==
!     == Output is MYCATOM (dftu.inc)                                   ==
!     == -------------------------------------------------------------- ==
     IMPLICIT NONE
!
      INTEGER, INTENT(IN)           ::  NSTATE
      REAL(real_8), INTENT(IN)      ::  TAU0(:,:,:)
      COMPLEX(real_8), INTENT(IN)   ::  PSI(:)
      LOGICAL                       ::  TFOR
      CHARACTER(*), PARAMETER       ::  procedureN = 'orthocatom'
!
      REAL(real_8)                  ::  FOC(M1SHLX),CANORM,FAC,DOTP
!
      COMPLEX(real_8) ,ALLOCATABLE  ::  myscr(:,:)!, myscr2(:)
      LOGICAL                       ::  TLSD_BAK,TIVAN_BAK
      INTEGER                       ::  ISUB,NOMAX,IAORB,IAT,IS,IA,NATST,I,J,&
                                        K,IZERO,LSCR_SUMMAT,LSCR_DGEMM,IUATM,&
                                        M1,IHUBL,L,ISH
!
      INTEGER                       ::  FIRSTCALL,ierr
      DATA       FIRSTCALL /0/
      SAVE       FIRSTCALL
!
      CALL TISET(procedureN,ISUB)
!
      CALL PHFAC(TAU0)
      CALL zeroing(mycatom)!,ncpw%ngw*hubbu%nuproj)
!
! r -> G space conversion of atomic orbitals. Usually not required as is in setbasis.crge%f 
!      CALL RTOG(1)
!
! Loading only the projector atomic orbitals
      IAORB=1
      DO IUATM=1,hubbu%nuatm
        IAT=hubbu%uatm(IUATM)
        IS=IATPT(2,IAT)
        DO ISH=1,atwf_mod%nshell(IS)
          L=atwf_mod%lshell(ISH,IS)
          DO IHUBL=1,hubbu%nl(IUATM)
            IF(ISH.EQ.hubbu%s(IUATM,IHUBL).AND.L.EQ.hubbu%l(IUATM,IHUBL))THEN
               CALL LOADC_U(MYCATOM(1,IAORB),ncpw%ngw,ncpw%ngw,hubbu%nuproj-IAORB+1,&
                          IS,IAT,ISH,L)
               IAORB=IAORB+2*L+1
            END IF
          END DO
        END DO
      END DO
!
! Normalize the loaded atomic orbitals
      IF(hubbu%pnorm.AND.(.NOT.pslo_com%tivan))THEN
        DO I=1,hubbu%nuproj
          CANORM=DOTP(ncpw%ngw,MYCATOM(1,I),MYCATOM(1,I))
!          CALL GLOSUM(1,CANORM)
          CALL mp_sum(canorm,parai%allgrp)
          CANORM=1.D0/DSQRT(CANORM)
          CALL ZDSCAL(ncpw%ngw,CANORM,MYCATOM(1,I),1)
        END DO
! Normalization factor <CA_m|S|CA_m>
      ELSEIF(hubbu%pnorm.AND.pslo_com%tivan)THEN
        CALL DCOPY(2*ncpw%ngw*hubbu%nuproj,MYCATOM,1,MYCATOM0,1)
        allocate(myscr(ncpw%ngw,hubbu%nuproj),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
        CALL DCOPY(2*ncpw%ngw*hubbu%nuproj,MYCATOM,1,MYSCR,1)
        TLSD_BAK=cntl%tlsd
        cntl%tlsd=.FALSE.
        CALL FNL_SET('SAVE')
        CALL FNLALLOC(hubbu%nuproj,.FALSE.,.FALSE.)
        CALL rnlsm(myscr,hubbu%nuproj,1,1,.false.)
        CALL spsi(hubbu%nuproj,myscr)
        CALL FNLDEALLOC(.FALSE.,.FALSE.)
        CALL FNL_SET('RECV')
        cntl%tlsd=TLSD_BAK
        DO I=1,hubbu%nuproj
          CANORM=DOTP(ncpw%ngw,MYCATOM(1,I),MYSCR(1,I))
!          CALL GLOSUM(1,CANORM)
          CALL mp_sum(canorm,parai%allgrp)
          CANORM=1.D0/DSQRT(CANORM)
          CALL ZDSCAL(ncpw%ngw,CANORM,MYCATOM(1,I),1)
        END DO
        deallocate(myscr,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
        
      END IF
!
! Orthogonalizae |CA_m>
      IF(hubbu%portho)THEN
        IF(TFOR) CALL STOPGM(&
                procedureN,'FORCES DUE TO ORTHOGONALIZATION OF '//&
                'PROJECTORS NOT IMPLEMENTED!', __LINE__,__FILE__)
        TLSD_BAK=cntl%tlsd
        cntl%tlsd=.FALSE.
        IF(pslo_com%tivan)THEN
          CALL FNL_SET('SAVE')
          CALL FNLALLOC(hubbu%nuproj,.FALSE.,.FALSE.)
          CALL rnlsm(mycatom,hubbu%nuproj,1,1,.false.)
        ENDIF
        allocate(myscr(ncpw%ngw,MAX(nstate,hubbu%nuproj)),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
        CALL ORTHO(hubbu%nuproj,MYCATOM,myscr)
        deallocate(myscr,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
        IF(pslo_com%tivan)THEN
          CALL FNLDEALLOC(.FALSE.,.FALSE.)
          CALL FNL_SET('RECV')
        ENDIF 
        cntl%tlsd=TLSD_BAK
      ENDIF
      FIRSTCALL=FIRSTCALL+1
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE ORTHOCATOM

!     ====================================================================
      SUBROUTINE PRINT_OCCMAT(TPRIOC,TLSD1,ISPIN0)
!     ====================================================================
      IMPLICIT NONE
!
      LOGICAL, INTENT(IN)     ::   TPRIOC,TLSD1
      INTEGER, INTENT(IN)     ::   ISPIN0
!
      INTEGER                 ::  IUATM,ISPIN,M1,M2,M10,M20,M0,ISPIN_MIN,&
                                  ISPIN_MAX,M1_MAX,M1_MIN,ISUB
      CHARACTER(len=25)       ::  MYFMT
      CHARACTER(*), PARAMETER ::   procedureN = 'print_occmat'
      LOGICAL                 ::   FERROR 
      REAL(real_8)            ::   OMT(hubbu%nuatm)
!
      INTEGER                 ::  IUNIT
!
      INTEGER                 ::  ICALL
      DATA       ICALL /0/
      SAVE       ICALL
!
      IF(.NOT.TPRIOC)RETURN
!
      CALL TISET(procedureN,ISUB)
!
      IF(hubbu%uverb)THEN
        IUNIT=6
      ELSE
        IUNIT=116
        CALL FILEOPEN(IUNIT,'OCCMAT',FO_APP,FERROR)
      END IF
!
      ICALL=ICALL+1
!
      CALL zeroing(OMT)
      ISPIN_MIN=1
      ISPIN_MAX=1
      IF(TLSD1)ISPIN_MAX=2
      IF(ISPIN0.GT.0)THEN
        ISPIN_MIN=ISPIN0
        ISPIN_MAX=ISPIN0
      END IF
      WRITE(IUNIT,'(/,A,I0)')' OCCUPATION MATRIX (M1 X M2):  STEP: ', iteropt%nfi
      DO IUATM=1,hubbu%nuatm
        M0=hubbu%muatm(1,1,IUATM)-1
        DO ISPIN=ISPIN_MIN,ISPIN_MAX
!          WRITE(116,'(/,A,I8)')' ICALL =',ICALL  
          WRITE(IUNIT,'(A,I8,A,I8)')' IUATM=',IUATM,'   ISPIN=',ISPIN
          DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
            M20=M2-M0
            M1_MIN=hubbu%muatm(1,1,IUATM)-M0
            M1_MAX=hubbu%muatm(1,2,IUATM)-M0
            MYFMT=' '
            WRITE(MYFMT,'(A1,I2,A6)')'(',M1_MAX-M1_MIN+1,'F12.5)'
            WRITE(IUNIT,MYFMT)&
                 (OM(ISPIN,IUATM,M10,M20),M10=M1_MIN,M1_MAX)
            OMT(IUATM)=OMT(IUATM)+OM(ISPIN,IUATM,M20,M20)
          END DO
        END DO
        WRITE(IUNIT,'(A,I5,A,F12.6)')&
                    ' Tr(oc.mat.)_Atm:',IUATM,' =',OMT(IUATM)
      END DO
      WRITE(IUNIT,*)' '
      IF(.NOT.hubbu%uverb)CALL FILECLOSE(IUNIT)
!
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE PRINT_OCCMAT
!
!     ====================================================================
      SUBROUTINE HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,TFOR,ISPIN0)
!     == -------------------------------------------------------------- ==
!     == This routine estimates the energy due to Habbard terms         ==
!     == Forces on ions due to Hubbard functional form is added to FION ==
!     == -------------------------------------------------------------- ==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER   ::  procedureN = 'hubbe'
      INTEGER                   ::  NSTATE,ISPIN0
      LOGICAL                   ::  TFOR
      REAL(real_8)              ::  TAU0(:,:,:),FION(:,:,:)
! NS TODO:  Was REAL variable! check this
      COMPLEX(real_8)           ::  PSI(:)
      COMPLEX(real_8)           ::  C0(:,:),C2U(:,:)!ncpw%ngw,:)
!
      REAL(real_8)              ::  OMSUMM2,OM1M2,OM2M1,OM1M1,OMSUM,FUI(3),&
                                    FOM1M1(3),FOM1M2(3),FOM2M1(3),FOMSUMM2(3),&
                                    EHUBA,OMSUMI
      INTEGER                   ::  IUATM,ISPIN,N_SPIN,M10,M20,M1,M2,M0,IS,IA,&
                                    IAT,K,IG,ISTATE,ISPIN_MIN,ISPIN_MAX,ISUB
!
      INTEGER                   ::  FIRSTCALL
      DATA       FIRSTCALL / 0 /
      SAVE       FIRSTCALL
!
      CALL TISET(procedureN,ISUB)
!
      N_SPIN=1
      IF(cntl%tlsd)N_SPIN=2
!
      IF(ISPIN0.EQ.0)THEN
       ISPIN_MIN=1
       ISPIN_MAX=N_SPIN
      ELSE
       ISPIN_MIN=ISPIN0
       ISPIN_MAX=ISPIN0
      END IF 
!
      hubbu%ehub=0.0
      IF(paral%io_parent)THEN
        DO IUATM=1,hubbu%nuatm
          DO ISPIN=ISPIN_MIN,ISPIN_MAX
            DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
              M10=M1-hubbu%muatm(1,1,IUATM)+1
              hubbu%ehub=hubbu%ehub+0.5D0*hubbu%u(IUATM)*OM(ISPIN,IUATM,M10,M10)+&
                       hubbu%a(IUATM)*OM(ISPIN,IUATM,M10,M10)
              DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
                M20=M2-hubbu%muatm(1,1,IUATM)+1
                hubbu%ehub=hubbu%ehub-0.5d0*hubbu%u(IUATM)*OM(ISPIN,IUATM,M10,M20)*&
                                      OM(ISPIN,IUATM,M20,M10)
              END DO  
            END DO
          END DO
        END DO
      END IF
!
!     Forces
      IF(TFOR)THEN
        CALL zeroing(FION_OM)
        IUATM=0
        IAT=0
        DO IS=1,ions1%nsp
          DO IA=1,ions0%na(IS)
            IAT=IAT+1
            IF(IS_U_AT(IAT))THEN
              IUATM=IUATM+1
              M0=hubbu%muatm(1,1,IUATM)-1
              CALL zeroing(FUI)!,3)
              DO ISPIN=1,N_SPIN
                DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
                  M10=M1-M0
                  DO K=1,3
                    FOM1M1(K)=FOMI(ISPIN,IUATM,M10,M10,K)
                  END DO
                  CALL zeroing(FOMSUMM2)!,3)
                  DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
                     M20=M2-M0
                     OM1M2=OM(ISPIN,IUATM,M10,M20)
                     OM2M1=OM(ISPIN,IUATM,M20,M10)
                     DO K=1,3
                       FOM1M2(K)=FOMI(ISPIN,IUATM,M10,M20,K)
                       FOM2M1(K)=FOMI(ISPIN,IUATM,M20,M10,K)
                       FOMSUMM2(K)=FOMSUMM2(K)+&
                                   FOM1M2(K)*OM2M1+OM1M2*FOM2M1(K)
                     END DO
                  END DO
                  DO K=1,3
                    FUI(K)=FUI(K)+(FOM1M1(K)-FOMSUMM2(K))&
                                 *0.5d0*hubbu%u(IUATM)+&
                                  FOM1M1(K)*hubbu%a(IUATM)
                  END DO
                END DO
              END DO
              DO K=1,3
                FION_OM(K,IA,IS)=FUI(K)
                IF(paral%io_parent)FION(K,IA,IS)=FION(K,IA,IS)-FION_OM(K,IA,IS)
              END DO
            END IF
          END DO
        END DO 
      END IF
!
! H^U|psi> ...the gradients on the wavefunction 
      CALL DEUDPSI2(C0,MYCATOM,C2U,TAU0,NSTATE,PSI,ISPIN0)
!
!dbg
!      IF(paral%io_parent)THEN
!        WRITE(6,'(/,A,F16.8,A)') ' TOTAL HUBBARD ENERGY    =', 
!     &                                        hubbu%ehub,' (A.U.)'
!        WRITE(6,'(/,A,F16.8,A)') ' RESPONSE HUBBARD ENERGY =', 
!     &                                        EHUBA,' (A.U.)'
!      END IF
!      IF(TFOR.AND.paral%io_parent)THEN
!        WRITE(6,'(/,A)') ' HUBBARD GRADIENTS:'
!        CALL WRGEOF(TAU0,FION_OM) 
!        WRITE(6,*)
!      END IF
!      STOP
!
      FIRSTCALL=FIRSTCALL+1
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE HUBBE
!     ====================================================================
      SUBROUTINE DEUDPSI2(C0,CA,C2U,TAU0,NSTATE,PSI,ISPIN0)
!     == -------------------------------------------------------------- ==
!     == This routine calcualtes the derivative of occupation matrix    ==
!     == with respect to the wavefunction;                              ==
!     == S operator is taken care for vanderbilt PPs.                   == 
!     == Ouput is in C2U                                                == 
!     == -------------------------------------------------------------- ==
      IMPLICIT NONE
!
      COMPLEX(real_8)             ::    C0(:,:),CA(:,:),C2U(:,:)
      REAL(real_8)                ::    TAU0(:,:,:)
      CHARACTER(*), PARAMETER     ::    procedureN = 'deudpsi2'
      
      INTEGER                     ::    NSTATE,ISPIN0,ISUB
!
      COMPLEX(real_8),ALLOCATABLE ::    SCA(:,:)
      COMPLEX(real_8)             ::    PSI(:)
!
      REAL(real_8)                ::    KD,PAR,FAC,TMP,FFI
      COMPLEX(real_8)             ::    UUU,UAA,FC,PAC,IM1M2,IM2M1,IM1M1,IFAC
      INTEGER                     ::    M1,M2,M10,M20,MM1,MM2,MOFF,IG,ISTATE,&
                                        ISPIN,IUATM,ISPIN_MIN,ISPIN_MAX
      LOGICAL                     ::    TLSD_BAK
!
      INTEGER                     ::    FIRSTCALL,ierr
      DATA       FIRSTCALL / 0 /
      SAVE       FIRSTCALL
!
      REAL(real_8)                ::     VHUB(2,hubbu%nuatm,7,7)
!
      CALL TISET(procedureN,ISUB)
!
      allocate(SCA(ncpw%ngw,hubbu%nuproj),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      CALL DCOPY(2*ncpw%ngw*hubbu%nuproj,CA,1,SCA,1)

!
       ISPIN_MIN=1
       ISPIN_MAX=1
       IF(cntl%tlsd)ISPIN_MAX=2
       IF(ISPIN0.GT.0)THEN
        ISPIN_MIN=ISPIN0
        ISPIN_MAX=ISPIN0
       END IF
!
! S|phi> --> SCA
      IF(pslo_com%tivan)THEN
        TLSD_BAK=cntl%tlsd
        cntl%tlsd=.FALSE.
        CALL FNL_SET('SAVE')
        CALL FNLALLOC(hubbu%nuproj,.FALSE.,.FALSE.)
!GM FIXME fix the call
        CALL rnlsm(sca,hubbu%nuproj,1,1,.false.)
        CALL spsi(hubbu%nuproj,sca)
        CALL FNLDEALLOC(.FALSE.,.FALSE.)
        CALL FNL_SET('RECV')
        cntl%tlsd=TLSD_BAK
      ENDIF      
!
      CALL zeroing(C2U)!,2*ncpw%ngw*NSTATE)
      CALL zeroing(VHUB)!,2*hubbu%nuatm*7*7)
      DO IUATM=1,hubbu%nuatm
       DO ISPIN=ISPIN_MIN,ISPIN_MAX
         DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
           M10=M1-hubbu%muatm(1,1,IUATM)+1
           VHUB(ISPIN,IUATM,M10,M10)=VHUB(ISPIN,IUATM,M10,M10)+&
            0.5d0*hubbu%u(IUATM)+hubbu%a(IUATM)
           DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
             M20=M2-hubbu%muatm(1,1,IUATM)+1
             VHUB(ISPIN,IUATM,M10,M20)=VHUB(ISPIN,IUATM,M10,M20)-&
               hubbu%u(IUATM)*OM(ISPIN,IUATM,M20,M10)
           END DO
         END DO
       END DO
      END DO
!
!$OMP parallel do private(ISPIN,MOFF,IUATM,M1,M10,MM1,TMP,M2,M20,MM2,IG)
      DO ISTATE=1,NSTATE
        ISPIN=1
        IF(cntl%tlsd.AND.(ISTATE.GT.spin_mod%nsup))ISPIN=2
        IF(ISPIN0.GT.0)ISPIN=ISPIN0
        MOFF=0
        DO IUATM=1,hubbu%nuatm
          DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
            M10=M1-hubbu%muatm(1,1,IUATM)+1
            MM1=MOFF+M10
            TMP=0.D0
            DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
               M20=M2-hubbu%muatm(1,1,IUATM)+1
               MM2=MOFF+M20
               TMP=TMP+VHUB(ISPIN,IUATM,M10,M20)*C0XCA(ISTATE,MM2)
            END DO
            IF((ISPIN0.EQ.0).AND.(crge%f(ISTATE,1).GT.1.D-3))TMP=TMP*crge%f(ISTATE,1)
            DO IG=1,ncpw%ngw
              C2U(IG,ISTATE)=C2U(IG,ISTATE)+&
                              DCMPLX(TMP,0.D0)*SCA(IG,MM1)
            END DO
          END DO
          MOFF=MOFF+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
        END DO
      END DO
!
!      IF(GEQ0)CALL ZCLEAN(C2U,NSTATE,ncpw%ngw)
!
      FIRSTCALL=FIRSTCALL+1
      deallocate(sca,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE DEUDPSI2
!
!     ==================================================================
      SUBROUTINE LOADUPROJ(IA,READ_PROJ_LOG)
!     Read U projectors from pseudo potential files; 
!     convert them from r to G space and store in CAT_U
      IMPLICIT NONE
!
      INTEGER                      ::  IA
      LOGICAL                      ::  READ_PROJ_LOG(maxsys%nsx)
      INTEGER                      ::  I,IS,IR,IL,NPOINT,ISTR1,ISTR2,IERR,&
                                       FIRSTCALL,MMAX,NWORK,N2,N22,ISHELL,&
                                       L,MSGLEN
      INTEGER, PARAMETER           ::  IFNUM=21, IMESHAT=256
      CHARACTER(len=200)           ::  FNAMES
      CHARACTER(*), PARAMETER      ::  procedureN = 'LOADUPROJ'
!
      REAL(real_8)                 ::  RP1,RP2,CLOGAT1,DISC,XMAX,RMIN,GMIN
      REAL(real_8), allocatable    ::  RP(:),TEMP(:,:),FINT(:),BSINT(:),&
                                       GG(:)
      COMPLEX(real_8), allocatable ::  WORK(:,:)
      LOGICAL                      ::  SAVED
      DATA FIRSTCALL /0/
      SAVE FIRSTCALL
!
      IS=cpsp(NAT_grm(IA))
!
      IF(FIRSTCALL.EQ.0)THEN
        FIRSTCALL=FIRSTCALL+1
        DO I=1,ions1%nsp
         READ_PROJ_LOG(I)=.FALSE.
        END DO
      END IF
      IF(READ_PROJ_LOG(IS))RETURN
      READ_PROJ_LOG(IS)=.TRUE.
!
      allocate(rp(maxsys%mmaxx),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      allocate(temp(maxsys%mmaxx,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
!     Vanderbilt PPs
      IF(paral%io_parent)THEN
        IF(pslo_com%tvan(IS))THEN
          NPOINT=ncpr1%meshva(IS)-1
          DO IR=1,NPOINT
            ATRG_U(IR,IS)=VDB_R(IR+1,IS)
            DO IL=1,atwf_mod%nshell(IS)
              ATWFR_U(IR,IL,IS)=VDB_PAWF(IS,IR+1,IL)
            END DO
          END DO
!     Norm Conserving
        ELSE
          FNAMES=ECPFILES(IS)
          CALL XSTRING(FNAMES,ISTR1,ISTR2)
          OPEN(UNIT=IFNUM,FILE=FNAMES(ISTR1:ISTR2),STATUS='UNKNOWN')
            IERR=INSCAN(IFNUM,'&WAVEFUNCTION')
          IF(IERR.NE.0)&
            CALL STOPGM(procedureN,&
              'READING PROJECTORS FROM PP FILE', __LINE__,&
              __FILE__)
          READ(IFNUM,*) NPOINT
          DO IR=1,NPOINT
             READ(IFNUM,ERR=99,FMT=*) ATRG_U(IR,IS),&
                 (ATWFR_U(IR,IL,IS),IL=1,atwf_mod%nshell(IS))
          ENDDO
        END IF
      END IF
!
      CALL mp_bcast(ATWFR_U,maxsys%mmaxx*maxsys%nsx*M1SHLX ,parai%source,parai%allgrp)
      CALL mp_bcast(ATRG_U,maxsys%mmaxx*maxsys%nsx,parai%source,parai%allgrp)
      CALL mp_bcast(NPOINT,parai%source,parai%allgrp)
!
      CALL DCOPY(NPOINT,ATRG_U(1,IS),1,RP(1),1)
      RP1=RP(1)
      RP2=RP(NPOINT)
      CALL CKGRID(RP1,RP2,ATRG_U(1,IS),IMESHAT,CLOGAT1)
      DO IL=1,atwf_mod%nshell(IS)
        CALL TGRID(RP,NPOINT,ATRG_U(1,IS),IMESHAT,&
                   ATWFR_U(1,IL,IS),&
                   maxsys%mmaxx,TEMP(1,1),TEMP(1,2),TEMP(1,3))
      ENDDO
!
      NWORK=2*maxsys%mmaxx
      allocate(work(NWORK,5),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      allocate(fint(nwork),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      allocate(bsint(nwork),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      allocate(gg(nwork),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      IF(.NOT.allocated(CAT_U)) &
        allocate(CAT_U(NSPLPO,2,atwp%m1shl,maxsys%nsx),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__)
      MMAX=IMESHAT
      XMAX=MMAX
      N2=NINT(LOG(XMAX)/LOG(2.D0)+0.499999D0)
      N22=2**N2
      RMIN=LOG(ATRG_U(1,IS))
      GMIN=LOG(SQRT(gvec_com%gcutw+GCUTKA)*parm%tpiba)-(MMAX-1)*CLOGAT1
      DO IL=1,MMAX
        GG(IL)=(EXP(GMIN+(IL-1)*CLOGAT1)/parm%tpiba)**2
      ENDDO
!
      DO ISHELL=1,atwf_mod%nshell(IS)
        SAVED=.FALSE.
        L=atwf_mod%lshell(ISHELL,IS)
        CALL zeroing(FINT)!,N22)
        DO IR=1,MMAX
          FINT(IR)=FPI*ATWFR_U(IR,ISHELL,IS)/ATRG_U(IR,IS)
        ENDDO
!       Fourier transformation
        CALL LSFBTR(FINT,BSINT,L,RMIN,GMIN,CLOGAT1,N2,SAVED,&
                  WORK(1,1),WORK(1,2),WORK(1,3),WORK(1,5),NWORK,DISC)
        CALL TGRID(GG,MMAX,GGNG,NSPLPO,BSINT,&
                     maxsys%mmaxx,TEMP(1,1),TEMP(1,2),TEMP(1,3))
        CALL DCOPY(NSPLPO,BSINT(1),1,CAT_U(1,1,ISHELL,IS),1)
        IF(L.GT.0.AND.GGNG(1).LT.1.D-12) CAT_U(1,1,ISHELL,IS)=0.0D0
        CALL CURV1(NSPLPO,GGNG,CAT_U(1,1,ISHELL,IS),0.0D0,0.0D0,3,&
                  CAT_U(1,2,ISHELL,IS),TEMP,0.0D0,IERR)
      ENDDO
!
      deallocate(WORK,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      deallocate(FINT,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      deallocate(BSINT,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      deallocate(GG,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
!
      deallocate(RP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      deallocate(TEMP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__)
      RETURN
99    CONTINUE
      CALL STOPGM(procedureN,'ERROR READING PP FILE', __LINE__,__FILE__)
      END SUBROUTINE LOADUPROJ
!
!     ==================================================================
      SUBROUTINE LOADC_U(C0,LDC,NNGW,NDIM,IS,IAT,ISH,L)
!     ==--------------------------------------------------------------==
!     == Calculate atomic orbital in plane wave basis.                ==
!     == For a given shell number ISH and angular momentum L,         ==
!     == 2L+1 atomic orbitals will be loaded in C0                    ==
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
!     Arguments
      INTEGER                   ::  NNGW,NDIM,IS,IAT,LDC,L,ISH
      COMPLEX(real_8)           ::  C0(LDC,NDIM)  
      CHARACTER(*), PARAMETER   ::  procedureN = 'LOADC_U'
!     Variables
      COMPLEX(real_8)           ::  CI,EI123
      REAL(real_8)              ::  VOL,CC
      INTEGER                   ::  ISUB,LXX,LL,IV,LY,IG,&
                                    LPP(5)
      DATA       LPP /0,1,4,9,16/
!     ==--------------------------------------------------------------==
      CALL TISET(procedureN,ISUB)
      VOL=1.D0/SQRT(parm%omega)
      LXX=0
      CI=(0.0D0,1.0D0)**L
      LL=2*L+1
      !FIXME check if zeroing is correct
      CALL zeroing(c0)
      DO IV=1,LL
        LXX=LXX+1
! FIXME :  old ZAZZERO CALL for full array?
!        CALL ZAZZERO(C0(1,LXX),LDC)
        LY=LPP(L+1)+IV
        IF(atwf_mod%nbcut.EQ.1) THEN
          IF(cntl%bigmem) THEN
!$OMP     parallel do private(IG,CC) shared(LY,NSPLPO)
           DO IG=1,NNGW
              CC=CURV2(HG(IG),NSPLPO,GGNH(1),CAT_U(1,1,ISH,IS),&
                 CAT_U(1,2,ISH,IS),0.0D0)*VOL
              C0(IG,LXX)=CI*YLMR(LY,IG,GK(1,1))*CC*EIGRB(IG,IAT)
           ENDDO
        ELSE
!$OMP     parallel do private(IG,CC,EI123) shared(LY,NSPLPO)
              DO IG=1,NNGW
                CC=CURV2(HG(IG),NSPLPO,GGNH(1),CAT_U(1,1,ISH,IS),&
                   CAT_U(1,2,ISH,IS),0.0D0)*VOL
                EI123=EI1(IAT,INYH(1,IG))*EI2(IAT,INYH(2,IG))*&
                      EI3(IAT,INYH(3,IG))
                C0(IG,LXX)=CI*YLMR(LY,IG,GK(1,1))*CC*EI123
              ENDDO
            ENDIF
        ELSEIF(atwf_mod%nbcut.EQ.0) THEN
            IF(NNGW.GT.ncpw%ngw) CALL STOPGM(procedureN,'ncpw%ngw', __LINE__,&
              __FILE__)
!$OMP     parallel do private(IG,CC) shared(LY,NSPLPO)
            DO IG=1,NNGW
              CC=CURV2(HG(IG),NSPLPO,GGNG(1),CAT_U(1,1,ISH,IS),&
                 CAT_U(1,2,ISH,IS),0.0D0)*VOL
              C0(IG,LXX)=CI*YLMR(LY,IG,GK(1,1))*CC*EIGR(IG,IAT,1)
            ENDDO
        ELSE
            CALL STOPGM(procedureN,'atwf_mod%nbcut', __LINE__,__FILE__)
        ENDIF
      ENDDO
      CALL TIHALT(procedureN,ISUB)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LOADC_U
!     ====================================================================
      SUBROUTINE DFTPLUSU_NUMDER(C0,C2U,TAU0,FION,NSTATE,PSI,TFOR)
!     == -------------------------------------------------------------- ==
!     == This routine calcualtes contributions from DFT+U part          ==
!     == -------------------------------------------------------------- ==
      IMPLICIT NONE
!
      INTEGER                       ::  NSTATE
      CHARACTER(*), PARAMETER       ::  procedureN = 'DFTPLUSU_NUMDER'
      COMPLEX(real_8)               ::  C0(:,:),PSI(:),C2U(:,:)
      REAL(real_8)                  ::  TAU0(3,maxsys%nax,maxsys%nsx),&
                                        FION(3,maxsys%nax,maxsys%nsx)     
      LOGICAL                       ::  TFOR,TIVAN_BAK
!
      REAL(real_8)                  ::  TAU0_BAK(3*maxsys%nax*maxsys%nsx),&
                                        EHUB_1,EHUB_2,FX,FY,FZ
!
      INTEGER                       ::  ISUB
!
!
      IF(paral%io_parent)WRITE(*,*)'===NUMERICAL DERIVATIVES==='
!
      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR GIVEN GEOMETRY: '
!
      CALL DCOPY(3*maxsys%nax*maxsys%nsx,TAU0,1,TAU0_BAK,1)

! Calculate the occupation matrix OM
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
!
! Estimate energy contribution from Habbard term
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
!
      FX=FION_OM(1,1,2)
      FY=FION_OM(2,1,2)
      FZ=FION_OM(3,1,2)
!
      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR Z(2)+DZ GEOMETRY: '
      TAU0(3,1,2)=TAU0(3,1,2)+0.001
! Calculate the occupation matrix OM
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_1=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR Z(2)-DZ GEOMETRY: '
      TAU0(3,1,2)=TAU0(3,1,2)-2*0.001
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_2=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)'---------> NUMERICAL DERIVATIVES Z:'
      IF(paral%io_parent)WRITE(*,*)(EHUB_1-EHUB_2)/(2.0*0.001),FZ,&
                          FZ-(EHUB_1-EHUB_2)/(2.0*0.001) 
!
      CALL DCOPY(3*maxsys%nax*maxsys%nsx,TAU0_BAK,1,TAU0,1)



      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR Y(2)+DY GEOMETRY: '
      TAU0(2,1,2)=TAU0(2,1,2)+0.001
! Calculate the occupation matrix OM
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_1=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR Y(2)-DY GEOMETRY: '
      TAU0(2,1,2)=TAU0(2,1,2)-2*0.001
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_2=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)'---------> NUMERICAL DERIVATIVES Y:'
      IF(paral%io_parent)WRITE(*,*)(EHUB_1-EHUB_2)/(2.0*0.001),FY,&
                FY-(EHUB_1-EHUB_2)/(2.0*0.001) 
!
      CALL DCOPY(3*maxsys%nax*maxsys%nsx,TAU0_BAK,1,TAU0,1)



      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR X(2)+DX GEOMETRY: '
      TAU0(1,1,2)=TAU0(1,1,2)+0.001
! Calculate the occupation matrix OM
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_1=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)' '
      IF(paral%io_parent)WRITE(*,*)'---------> FOR X(2)-DX GEOMETRY: '
      TAU0(1,1,2)=TAU0(1,1,2)-2*0.001
      CALL OCCMAT(C0,TAU0,NSTATE,PSI,.true.,0)
      CALL HUBBE(TAU0,FION,C0,C2U,PSI,NSTATE,.true.,0)
      EHUB_2=hubbu%ehub
!
      IF(paral%io_parent)WRITE(*,*)'---------> NUMERICAL DERIVATIVES X:'
      IF(paral%io_parent)WRITE(*,*)(EHUB_1-EHUB_2)/(2.0*0.001),FX, &
                          FX-(EHUB_1-EHUB_2)/(2.0*0.001)
!
      CALL DCOPY(3*maxsys%nax*maxsys%nsx,TAU0_BAK,1,TAU0,1)
!
      END SUBROUTINE DFTPLUSU_NUMDER
!
!  ==================================================================
  SUBROUTINE give_scr_hubbardu(NSTATE,L_DFTU,TAG)
!  ==--------------------------------------------------------------==
  IMPLICIT NONE
!  Arguments
    INTEGER                     ::  L_DFTU,NSTATE
    CHARACTER(*), PARAMETER     ::  procedureN = 'give_scr_hubbardu'
    CHARACTER(len=30)           ::  TAG

    INTEGER                     ::  L_ORTHO,L_RNLSM,L_SPSI,NOMAX
!  ==--------------------------------------------------------------==
  IF(cntl%thubb)THEN 
    NOMAX=MAX(hubbu%nuproj,NSTATE)
    CALL GIVE_SCR_ORTHO(L_ORTHO,TAG,NOMAX)
    CALL GIVE_SCR_RNLSM(L_RNLSM,TAG,NOMAX,.FALSE.)
    CALL GIVE_SCR_SPSI(L_SPSI,TAG)
    L_DFTU=MAX(1,L_ORTHO)
    L_DFTU=MAX(L_DFTU,L_RNLSM)
    L_DFTU=MAX(L_DFTU,L_SPSI)
  ELSE
    L_DFTU=0
  END IF
!  ==--------------------------------------------------------------==
  RETURN
  END SUBROUTINE give_scr_hubbardu
!!  ==================================================================
!

  SUBROUTINE add_hubbardu(c2l,c2ul,nstate)
     COMPLEX(real_8), INTENT(INOUT) ::  c2l(:,:)
     COMPLEX(real_8), INTENT(IN)    ::  c2ul(:,:)
     INTEGER , INTENT(IN)           ::  nstate
     INTEGER                        ::  istate,IG
     CHARACTER(*), PARAMETER        ::  procedureN = 'add_hubbardu'
!   !$OMP parallel do 
      DO istate=1,nstate
        DO ig=1,ncpw%ngw
          c2l(ig,istate)=c2l(ig,istate)-c2ul(ig,istate)
        END DO
      END DO
  END SUBROUTINE add_hubbardu
END MODULE hubbardu_utils

