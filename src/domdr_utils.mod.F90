MODULE domdr_utils

USE cppt,                           ONLY: gk,&
                                          twnl,&
                                          hg,inyh
USE cvan,                           ONLy: qq
USE hubbardu
USE elct,                           ONLY: crge
USE error_handling,                 ONLY: stopgm
USE fnlalloc_utils,                 ONLY: fnlalloc,&
                                          fnldealloc,&
                                          fnl_set
USE geq0mod,                        ONLY: geq0
USE ions,                           ONLY: ions0,&
                                          ions1
USE kinds,                          ONLY: real_8
USE kpts,                           ONLY: tkpts
USE mp_interface,                   ONLY: mp_bcast,&
                                          mp_sum
USE nlps,                           ONLY: nghtol,&
                                          nlps_com
USE  parac,                         ONLY: paral,&
                                          parai
USE pslo,                           ONLY: pslo_com
USE rnlsm_utils,                    ONLY: give_scr_rnlsm,rnlsm
USE sfac,                           ONLY: dfnl,&
                                          eigr,&
                                          ei1,ei2,ei3,eigrb,&
                                          fnl
USE spin,                           ONLY: spin_mod
USE spsi_utils,                     ONLY: spsi,&
                                          give_scr_spsi
USE system,                         ONLY: cntl,&
                                          iatpt,&
                                          maxsys,&
                                          ncpw,&
                                          parm,parap
USE timer,                          ONLY: tihalt,&
                                          tiset
USE zeroing_utils,                  ONLY: zeroing


PUBLIC :: domdr

CONTAINS
!     ====================================================================
      SUBROUTINE DOMDR(C0,CA,c0xca,TAU0,NSTATE,PSI,TFOR)
!     == -------------------------------------------------------------- ==
!     == This routine calcualtes the derivative of occupation matrix    ==
!     == with respect to the nuclear coordinates.                       ==
!     == S operator is taken care for vanderbilt PPs.                   == 
!     == Ouput is in FOMI (dftu.inc)                                    == 
!     == -------------------------------------------------------------- ==
      IMPLICIT NONE
!!
      INTEGER                       ::  NSTATE
      COMPLEX(real_8)               ::  C0(:,:),CA(:,:),PSI(:)
      REAL(real_8),INTENT(IN)       ::  C0XCA(:,:)
      REAL(real_8)                  ::  TAU0(:,:,:)
      LOGICAL                       ::  TFOR
      CHARACTER(*), PARAMETER       ::  procedureN = 'DOMDR'
!
      CHARACTER(len=30)             ::  TAG
      INTEGER                       ::  ISUB,ISPIN,ISTATE,K,IA,IS,IUATM,M1,M2,&
                                        IATTOT,IV,JV,M10,M20,M0,ISA0,ISA,ISP,&
                                        IAT,MOFF,MM1,MM2
      INTEGER   , allocatable       ::  IVJV(:,:,:)
      INTEGER                       ::  NIJV(ions1%nsp),I,OFFSET(hubbu%nuatm)
!
      COMPLEX(real_8), allocatable  ::  dca(:,:,:)

      REAL(real_8), allocatable     ::  c0xdca(:,:,:),&
                                        c0xb(:,:,:),c0xdb(:,:,:,:),caxb(:,:,:),&
                                        caxdb(:,:,:,:), dcaxb(:,:,:,:)
!
      REAL(real_8), allocatable     ::  dc0xca(:,:,:,:)     
!
      REAL(real_8)                  ::  DUM,FF
!
!     REAL(real_8)                  ::  tt, tte, tt1,tt1e, tt2, tt2e, tt3, tt3e
!     REAL(real_8)                  ::  tt4, tt4e, tt5,tt5e, tt6,tt6e,tt7,tt7e
!     REAL(real_8)                  ::  TIMEF
!
      INTEGER                       ::  FIRSTCALL, ierr
      DATA       FIRSTCALL / 0 /
      SAVE       FIRSTCALL
!
      CALL TISET(procedureN,ISUB)
!     tt=TIMEF()
!
      allocate(dca(ncpw%ngw,hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(c0xdca(nstate,hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(dc0xca(nstate,hubbu%nuproj,hubbu%nuatm,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      IF(pslo_com%tivan)THEN
        allocate(c0xb(maxsys%nhxs,nstate,hubbu%nuatm),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
        allocate(caxb(maxsys%nhxs,hubbu%nuproj,hubbu%nuatm),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
        allocate(c0xdb(maxsys%nhxs,nstate,3,hubbu%nuatm),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
        allocate(caxdb(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
        allocate(dcaxb(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm),STAT=ierr)
      END IF
!
!  Need to consider the overlap operator of USPP: S = 1 + SUM_a,b,I q_a,b |beta> <beta| 
!     <C0|DCA>+ sum_I->a,b q_a,b * {<C0|B_a>  * <DCA|B_b> +
!                                   <C0|B_a>  * <CA|DB_b> +
!                                   <C0|DB_a> * <CA|B_b>  } 
!
!     tt1=TIMEF()
      OFFSET(1)=0
      DO IUATM=1,hubbu%nuatm-1
        OFFSET(IUATM+1)=OFFSET(IUATM)+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
      ENDDO
      IF(pslo_com%tivan)THEN
        IF(hubbu%portho) THEN
           CALL STOPGM(procedureN,&
                'FORCES FOR REQUESTED PROJECTOR NOT IMPLEMNTED',&
                  __LINE__,__FILE__)
        ELSE
          IF(hubbu%pnorm)THEN
             CALL HUBINTS_USPP2(NSTATE,C0,CA,C0XDCA,C0XB,&
                             CAXB,C0XDB,CAXDB,DCAXB,DCA,PSI)
          ELSE
             CALL HUBINTS_USPP(NSTATE,C0,CA,C0XDCA,C0XB,&
                             CAXB,C0XDB,CAXDB,DCAXB,DCA)
          ENDIF
        END IF
!        tt1e=TIMEF()
!         CALL HUBINTS_CHECK(NSTATE,C0,CA,C0XDCA,C0XB,
!     &                     CAXB,C0XDB,CAXDB,DCAXB,DCA,PSI)


! d<C0|S|CA>
 

!       tt2=TIMEF()
        allocate(IVJV(2,maxsys%nhxs*maxsys%nhxs,ions1%nsp),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
        call zeroing(NIJV)
        DO IS=1,ions1%nsp
! nsiemer  lookup matrxix (IV,JV,IS) for  pairs with QQ(IV,JV,IS) > 1.D-5 
!          so no IF in the ISTATE loop. 
           DO IV=1,nlps_com%ngh(IS)
             DO JV=1,nlps_com%ngh(IS)
               IF(DABS(QQ(IV,JV,IS)).GT.1.D-5)THEN
                 NIJV(IS)=NIJV(IS)+1
                 IVJV(1,NIJV(IS),IS)=IV
                 IVJV(2,NIJV(IS),IS)=JV
               ENDIF
             ENDDO
           ENDDO
        ENDDO
        MOFF=0
!$OMP parallel do private(IAT,IS,IA,K,M1,M10,MM1,ISTATE,DUM,I,IV,JV)
        DO IUATM=1,hubbu%nuatm
          IAT=hubbu%uatm(IUATM)
          IS=IATPT(2,IAT)
          IA=IATPT(1,IAT)
          DO K=1,3
            DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
              M10=M1-hubbu%muatm(1,1,IUATM)+1
!              MM1=MOFF+M10
              MM1=OFFSET(IUATM)+M10
              DO ISTATE=1,NSTATE
                DUM=0.D0
                DO I=1,NIJV(IS)
                  IV=IVJV(1,I,IS)
                  JV=IVJV(2,I,IS)
                  DUM=DUM+QQ(IV,JV,IS)*&
                      (C0XDB(IV,ISTATE,K,IUATM)*CAXB(JV,MM1,IUATM)   +&
                       C0XB(IV,ISTATE,IUATM)   *CAXDB(JV,MM1,K,IUATM)+&
                       C0XB(IV,ISTATE,IUATM)   *DCAXB(JV,MM1,K,IUATM) )
                END DO 
                DC0XCA(ISTATE,M10,IUATM,K)=C0XDCA(ISTATE,MM1,K)+DUM
              END DO
            END DO
          END DO
!          MOFF=MOFF+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
         END DO
!$OMP end parallel do
      ELSE
         CALL HUBINTS(NSTATE,C0,CA,C0XDCA,DCA)
         DO K=1,3
           MOFF=0
           DO IUATM=1,hubbu%nuatm
             DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM) 
               M10=M1-hubbu%muatm(1,1,IUATM)+1
               MM1=MOFF+M10
               CALL DCOPY(NSTATE,C0XDCA(1,MM1,K),1,&
                                 DC0XCA(1,M10,IUATM,K),1)
             END DO
             MOFF=MOFF+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
           END DO
         END DO
      ENDIF
!     tt2e=TIMEF()
!
!
!    dn(m1,m2)/dR_I= sum_i,s d<C0|S|CA>*<C0|S|CA>+<C0|S|CA>*d<C0|S|CA>
!
      CALL zeroing(FOMI)!,2*hubbu%nuatm*7*7*3)
!     tt3=TIMEF()
      MOFF=0
!$OMP parallel do private(K,M1,M10,MM1,M2,M20,MM2,ISTATE,ISPIN)
      DO IUATM=1,hubbu%nuatm
        DO K=1,3
          DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
             M10=M1-hubbu%muatm(1,1,IUATM)+1
             MM1=OFFSET(IUATM)+M10
             DO M2=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
                M20=M2-hubbu%muatm(1,1,IUATM)+1
                MM2=OFFSET(IUATM)+M20
                DO ISTATE=1,NSTATE
                   ISPIN=1
                   IF(ISTATE.GT.spin_mod%nsup)ISPIN=2
                      FOMI(ISPIN,IUATM,M10,M20,K)=&
                      FOMI(ISPIN,IUATM,M10,M20,K)+&
                      crge%f(ISTATE,1)*(DC0XCA(ISTATE,M10,IUATM,K)*&
                          C0XCA(ISTATE,MM2)+&
                          C0XCA(ISTATE,MM1)*DC0XCA(ISTATE,M20,IUATM,K))
                END DO 
             END DO
          END DO
!          MOFF=MOFF+hubbu%muatm(1,2,IUATM)-hubbu%muatm(1,1,IUATM)+1
        END DO
      END DO
!     tt3e=TIMEF()
!
!     tt4=tt3e
      DO K=1,3
        DO ISPIN=1,2
          DO IUATM=1,hubbu%nuatm
            M0=hubbu%muatm(1,1,IUATM)-1
            DO M1=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM) 
               M10=M1-M0
               DO M2=M1+1,hubbu%muatm(1,2,IUATM) 
                 M20=M2-M0
                 FOMI(ISPIN,IUATM,M10,M20,K)=FOMI(ISPIN,IUATM,M20,M10,K)
               END DO
             END DO
          END DO
        END DO
      END DO
!     tt4e=TIMEF()
!
      FIRSTCALL=FIRSTCALL+1
!
      deallocate(DCA,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
      deallocate(C0XDCA,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
      deallocate(DC0XCA,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
      IF(pslo_com%tivan)THEN
        deallocate(C0XB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
        deallocate(CAXB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
        deallocate(C0XDB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
        deallocate(CAXDB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
        deallocate(DCAXB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
      END IF
!
!      tte=TIMEF()
!      if (parent) print *, 'nsiemer Time DOMDR =',tte-tt
!      if (parent) print *, 'nsiemer: HUBINTS Call:', tt1e-tt1
!      if (parent) print *, 'nsiemer: loops:', tt2e-tt2, tt3e-tt3,tt4e-tt4
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE DOMDR
!
!     ====================================================================
      SUBROUTINE HUBINTS_USPP(NSTATE,C0,CA,C0XDCA,&
                         C0XB,CAXB,C0XDB,CAXDB,DCAXB,DCA)
!     ====================================================================
!     == Integrals needed for gradients of occupation matrix; USPP      ==
!     ====================================================================
      IMPLICIT NONE
!
      INTEGER                     :: NSTATE
!
      COMPLEX(real_8)             :: C0(:,:),&
                                     DCA(:,:,:),CA(:,:)
!
      REAL(real_8)                :: C0XDCA(NSTATE,hubbu%nuproj,3),&
                                   C0XB(maxsys%nhxs,NSTATE,hubbu%nuatm),&
                                   C0XDB(maxsys%nhxs,NSTATE,3,hubbu%nuatm),& 
                                   CAXB(maxsys%nhxs,hubbu%nuproj,hubbu%nuatm),&
                                   CAXDB(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm),& 
                                   DCAXB(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm)
      CHARACTER(*), PARAMETER     :: procedureN = 'HUBINTS_USPP'
!
      COMPLEX(real_8)             :: CIL,IBE,IIGBE,IGV
! GM  MYBETA(ncpw%ngw,IV,hubbu%nuatm) MYDBETA(ncpw%ngw,3,IV,hubbu%nuatm) allocatable and alloc
      COMPLEX(real_8),allocatable :: MYBETA(:,:,:), MYDBETA(:,:,:,:) 
      
!     COMPLEX(real_8) CI,MYBETA(ncpw%ngw),MYDBETA(ncpw%ngw,3)
      COMPLEX(real_8), parameter  :: CI=(0.D0,-1.D0)
!
      INTEGER                     :: IA,IS,IAT,IV,K,IATTOT,ISTATE,IG,ISUB,IUATM
      REAL(real_8)                :: BF,GV,DDOT
      REAL(real_8)                :: GEQFAC
      INTEGER                     :: ierr
!     REAL(real_8)      tt, tte, tt1,tt1e, tt2, tt2e, tt3, tt3e
!     REAL(real_8)      tt4, tt4e, tt5,tt5e, tt6, tt6e, tt7, tt7e
!     REAL(real_8)     TIMEF
 
      CALL TISET(procedureN,ISUB)
!    tt=TIMEF()
      GEQFAC=1.D0
      IF(GEQ0)GEQFAC=0.5D0
 
      allocate(MYBETA(ncpw%ngw,maxsys%nhxs,hubbu%nuatm),&
            MYDBETA(ncpw%ngw,maxsys%nhxs,3,hubbu%nuatm),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
!
      IF(tkpts%tkpnt)CALL STOPGM(procedureN,'K POINTS & DFT+U NOT IMPLEMENTED',&
             __LINE__,__FILE__)
!
      CALL zeroing(C0XDCA)!,NSTATE*hubbu%nuproj*3)
!    tt1=TIMEF()
!$OMP parallel do private(K,IATTOT,IGV,IG) collapse(3) schedule(static)
      DO K=1,3
        DO IATTOT=1,hubbu%nuproj
           DO IG=1,ncpw%ngw
             IGV=DCMPLX(0.D0,-1.D0)*DCMPLX(GK(K,IG)*parm%tpiba,0.D0)
             DCA(IG,IATTOT,K)=CA(IG,IATTOT)*IGV
           END DO
        END DO
      END DO
!     tt1e=TIMEF()
!
!     tt2=TIMEF()
      DO K=1,3
        CALL OVLAP2(ncpw%ngw,NSTATE,hubbu%nuproj,C0XDCA(1,1,K),C0,DCA(1,1,K),.True.)
      END DO
!     tt2e=TIMEF()
!
!      CALL GLOSUM(NSTATE*hubbu%nuproj*3,C0XDCA)
      CALL mp_sum(c0xdca,nstate*hubbu%nuproj*3,parai%allgrp)
!
      
!     tt3=TIMEF()
!$OMP parallel do &
!$OMP& private(IUATM,IAT,IS,IA,IV,BF,IBE,GV,IIGBE)
      DO IUATM=1,hubbu%nuatm
        IAT=hubbu%uatm(IUATM)
        IS=IATPT(2,IAT)
        IA=IATPT(1,IAT)
        DO IV=1,nlps_com%ngh(IS)
          CIL=DCMPLX(0.D0,-1.D0)**NGHTOL(IV,IS) 
          DO IG=1,ncpw%ngw
! beta function |beta>
            BF=TWNL(IG,IV,IS,1)     
            IBE=CIL*EIGR(IG,IAT,1)*BF  
            MYBETA(IG,IV,IUATM)=IBE           
! derivative of beta function |d_beta>
            DO K=1,3
              GV=GK(K,IG)*parm%tpiba
              IIGBE=DCMPLX(0.D0,-1.D0)*GV*IBE 
              MYDBETA(IG,IV,K,IUATM)=IIGBE 
            END DO
          END DO 
          MYBETA(1,IV,IUATM)=MYBETA(1,IV,IUATM)*GEQFAC 
          DO K=1,3
            MYDBETA(1,IV,K,IUATM)=MYDBETA(1,IV,K,IUATM)*GEQFAC
          END DO
        ENDDO
! <psi|beta>
        CALL DGEMM('t','n',nlps_com%ngh(IS),NSTATE,2*ncpw%ngw,2.0d0,&
                MYBETA(1,1,IUATM),2*ncpw%ngw,C0(1,1),2*ncpw%ngw,0.0d0,&
                C0XB(1,1,IUATM),maxsys%nhxs)
! <psi|d_beta>
        DO K=1,3
            CALL DGEMM('t','n',nlps_com%ngh(IS),NSTATE,2*ncpw%ngw,2.0d0,&
                MYDBETA(1,1,K,IUATM),2*ncpw%ngw,C0(1,1),2*ncpw%ngw,0.0d0,&
                C0XDB(1,1,K,IUATM),maxsys%nhxs)
        ENDDO
! <phi|beta>
        CALL DGEMM('t','n',nlps_com%ngh(IS),hubbu%nuproj,2*ncpw%ngw,2.0d0,&
                MYBETA(1,1,IUATM),2*ncpw%ngw,CA(1,1),2*ncpw%ngw,0.0d0,&
                CAXB(1,1,IUATM),maxsys%nhxs)
!  <phi|d_beta>
        DO K=1,3
          CALL DGEMM('t','n',nlps_com%ngh(IS),hubbu%nuproj,2*ncpw%ngw,2.0d0,&
               MYDBETA(1,1,K,IUATM),2*ncpw%ngw,CA(1,1),2*ncpw%ngw,0.0d0,&
                CAXDB(1,1,K,IUATM),maxsys%nhxs)
          CALL DGEMM('t','n',nlps_com%ngh(IS),hubbu%nuproj,2*ncpw%ngw,2.0d0,&
                MYBETA(1,1,IUATM),2*ncpw%ngw,DCA(1,1,K),2*ncpw%ngw,0.0d0,&
                DCAXB(1,1,K,IUATM),maxsys%nhxs)
        ENDDO 
      END DO
!$OMP end parallel do
      deallocate(MYBETA,MYDBETA,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
!     tt3e=TIMEF()
      CALL mp_sum(C0XB ,hubbu%nuatm*maxsys%nhxs*NSTATE  ,parai%allgrp)
      CALL mp_sum(C0XDB,hubbu%nuatm*maxsys%nhxs*NSTATE*3,parai%allgrp)
      CALL mp_sum(CAXB ,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj  ,parai%allgrp)
      CALL mp_sum(CAXDB,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)
      CALL mp_sum(DCAXB,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)
!
!     tte=TIMEF()
!     if (parent) print *, 'nsiemer Time HUBINTS_USPP =',tte-tt
!     if (parent) print *, ' loops:', tt1e-tt1, tt2e-tt2, tt3e-tt3 
!     if (parent) print *, ' glosum', tte -tt3e
      CALL TIHALT(procedureN,ISUB)
      END SUBROUTINE HUBINTS_USPP
!     ====================================================================


!     ====================================================================
      SUBROUTINE HUBINTS_USPP2(NSTATE,C0,CA,C0XDCA,&
                        C0XB,CAXB,C0XDB,CAXDB,DCAXB,DCA,PSI)
!     ====================================================================
!     == Integrals needed for gradients of occupation matrix; USPP      ==
!     ====================================================================
      IMPLICIT NONE
!
      INTEGER                   :: NSTATE
!
      COMPLEX(real_8)           :: C0(:,:),&
                                   DCA(:,:,:),&
                                   CA(:,:),PSI(:)
!
      REAL(real_8)              :: C0XDCA(NSTATE,hubbu%nuproj,3),&
                                 C0XB(maxsys%nhxs,NSTATE,hubbu%nuatm),&
                                 C0XDB(maxsys%nhxs,NSTATE,3,hubbu%nuatm),&
                                 CAXB(maxsys%nhxs,hubbu%nuproj,hubbu%nuatm),&
                                 CAXDB(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm),& 
                                 DCAXB(maxsys%nhxs,hubbu%nuproj,3,hubbu%nuatm)
      CHARACTER(*),PARAMETER    :: procedureN = 'HUBINTS_USPP2'
!
      COMPLEX(real_8)           :: CIL,IBE,IIGBE,IGV
      COMPLEX(real_8)           :: MYBETA(ncpw%ngw),MYDBETA(ncpw%ngw,3)
      COMPLEX(real_8),PARAMETER :: CI=(0.D0,-1.D0)
!
      INTEGER                   :: IA,IS,IV,JV,K,IATTOT,ISTATE,IG,IAT
      REAL(real_8)              :: BF,GV,DDOT,DOTP
!
      INTEGER                   :: FIRSTCALL
      DATA    FIRSTCALL /0/
      SAVE    FIRSTCALL
 
      COMPLEX(real_8),allocatable :: SCA(:,:),DCA0(:,:,:)
      REAL(real_8)              :: DUM
      REAL(real_8), allocatable :: DCA0SCA0(:,:,:),CA0SCA0(:),&
                                   DCA0XCA0(:,:),CA0B(:,:,:),&
                                   DCA0B(:,:,:,:),CA0DB(:,:,:,:)
!
      INTEGER                   :: IUATM,ISUB,ierr
      LOGICAL                   :: TLSD_BAK
!     REAL(real_8)      tt, tte, tt1,tt1e, tt2, tt2e, tt3, tt3e
!     REAL(real_8)      tt4, tt4e, tt5,tt5e, tt6, tt6e, tt7, tt7e
!     REAL(real_8)     TIMEF
!     EXTERNAL   TIMEF
!
      CALL TISET(procedureN,ISUB)
      IF(tkpts%tkpnt)CALL STOPGM(procedureN,&
        'K POINTS & DFT+U NOT IMPLEMENTED', __LINE__,__FILE__)
!     tt=TIMEF()
!
      allocate(SCA(ncpw%ngw,hubbu%nuproj),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(DCA0(ncpw%ngw,hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(DCA0SCA0(hubbu%nuproj,hubbu%nuatm,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(CA0SCA0(hubbu%nuproj),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(DCA0XCA0(hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(CA0B(ions1%nat,maxsys%nhxs,hubbu%nuproj),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(DCA0B(ions1%nat,maxsys%nhxs,hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
      allocate(CA0DB(ions1%nat,maxsys%nhxs,hubbu%nuproj,3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
             __LINE__,__FILE__) 
!
      CALL zeroing(C0XDCA)!,NSTATE*hubbu%nuproj*3)
!     tt1=TIMEF()
      DO K=1,3
        DO IATTOT=1,hubbu%nuproj
           DO IG=1,ncpw%ngw
             IGV=DCMPLX(0.D0,-1.D0)*DCMPLX(GK(K,IG)*parm%tpiba,0.D0)
             DCA(IG,IATTOT,K)=CA(IG,IATTOT)*IGV
             DCA0(IG,IATTOT,K)=MYCATOM0(IG,IATTOT)*IGV
           END DO
        END DO
      END DO
!     tt1e=TIMEF()
!
! ---------------------------------------------------------------------
! <phi|S|phi> and <dphi|S|phi>
! ---------------------------------------------------------------------
! 
      CALL DCOPY(2*ncpw%ngw*hubbu%nuproj,MYCATOM0,1,SCA,1)
!
      TLSD_BAK=cntl%tlsd
      cntl%tlsd=.FALSE.
      CALL FNL_SET('SAVE')
      CALL FNLALLOC(hubbu%nuproj,.FALSE.,.FALSE.)
      CALL rnlsm(sca,hubbu%nuproj,1,1,.false.)
!
! S|phi>
      CALL spsi(hubbu%nuproj,sca) 
!
      CALL FNLDEALLOC(.FALSE.,.FALSE.)
      CALL FNL_SET('RECV')
      cntl%tlsd=TLSD_BAK
!
! <phi_m|S|phi_m>
      DO IATTOT=1,hubbu%nuproj
        CA0SCA0(IATTOT)=DOTP(ncpw%ngw,MYCATOM0(1,IATTOT),SCA(1,IATTOT))
      END DO
     ! CALL GLOSUM(hubbu%nuproj,CA0SCA0)
      CALL mp_sum(ca0sca0,hubbu%nuproj,parai%allgrp)
! <dphi_m,k|phi_m>
      DO K=1,3
        DO IATTOT=1,hubbu%nuproj
          DCA0XCA0(IATTOT,K)=DOTP(ncpw%ngw,DCA0(1,IATTOT,K),MYCATOM0(1,IATTOT))
        END DO
      END DO
      !CALL GLOSUM(hubbu%nuproj*3,DCA0SCA0)
      CALL mp_sum(DCA0SCA0,hubbu%nuproj*3,parai%allgrp)
!
      CALL zeroing(CA0B)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj)
      CALL zeroing(DCA0B)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3)
      CALL zeroing(CA0DB)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3)
! <dphi|beta> & <phi|beta> 
!     tt2=TIMEF()
      DO IUATM=1,hubbu%nuatm
        IAT=hubbu%uatm(IUATM)
        IS=IATPT(2,IAT)
        IA=IATPT(1,IAT)
        DO IV=1,nlps_com%ngh(IS)
          CIL=DCMPLX(0.D0,-1.D0)**NGHTOL(IV,IS) ! -i^l
          DO IG=1,ncpw%ngw
             BF=TWNL(IG,IV,IS,1)     
             IBE=CIL*EIGR(IG,IAT,1)*DCMPLX(BF,0.D0) 
             MYBETA(IG)=IBE           
!
             DO K=1,3
               GV=GK(K,IG)*parm%tpiba
               IIGBE=DCMPLX(0.D0,-1.D0)*DCMPLX(GV,0.D0)*IBE 
               MYDBETA(IG,K)=IIGBE 
             END DO
          END DO 
          IF(GEQ0)THEN
            MYBETA(1)=MYBETA(1)*DCMPLX(0.5D0,0.d0)
            DO K=1,3
              MYDBETA(1,K)=MYDBETA(1,K)*DCMPLX(0.5D0,0.D0)
            END DO
          END IF
!
          DO IATTOT=1,hubbu%nuproj
            CA0B(IUATM,IV,IATTOT)=CA0B(IUATM,IV,IATTOT)+&
                 2.D0*DDOT(2*ncpw%ngw,MYBETA,1,MYCATOM0(1,IATTOT),1)
            DO K=1,3
              DCA0B(IUATM,IV,IATTOT,K)=DCA0B(IUATM,IV,IATTOT,K)+&
                  2.D0*DDOT(2*ncpw%ngw,MYBETA,1,DCA0(1,IATTOT,K),1)
              CA0DB(IUATM,IV,IATTOT,K)=CA0DB(IUATM,IV,IATTOT,K)+&
                  2.D0*DDOT(2*ncpw%ngw,MYDBETA(1,K),1,MYCATOM0(1,IATTOT),1)
            END DO
          END DO
        END DO
      END DO
!     tt2e=TIMEF()
      CALL mp_sum(CA0B ,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj  ,parai%allgrp)
      CALL mp_sum(DCA0B,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)
      CALL mp_sum(CA0DB,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)
!
!     tt3=TIMEF()
      DO K=1,3
        DO IUATM=1,hubbu%nuatm
         IAT=hubbu%uatm(IUATM)
         IS=IATPT(2,IAT)
         IA=IATPT(1,IAT)
         DO IATTOT=hubbu%muatm(1,1,IUATM),hubbu%muatm(1,2,IUATM)
           DUM=0.D0
           DO IV=1,nlps_com%ngh(IS)
             DO JV=1,nlps_com%ngh(IS)
               DUM=DUM+QQ(IV,JV,IS)*&
                  (DCA0B(IAT,IV,IATTOT,K)*CA0B(IAT,JV,IATTOT)+&
                   CA0DB(IAT,IV,IATTOT,K)*CA0B(IAT,JV,IATTOT)+&
                   CA0B(IAT,IV,IATTOT)*CA0DB(IAT,JV,IATTOT,K)+&
                   CA0B(IAT,IV,IATTOT)*DCA0B(IAT,JV,IATTOT,K))
             END DO
           END DO
           DCA0SCA0(IATTOT,IUATM,K)=2.d0*DCA0XCA0(IATTOT,K)+DUM
         END DO
        END DO
      END DO
!     tt3e=TIMEF()
!      CALL MEMORY_CHECK
!
!  -<phi^I|/(2*<phi^I|S|phi^I>^3/2) * d/dR_I( <phi^I|S|phi^I>)  
!     tt4=TIMEF()
      DO K=1,3
        DO IUATM=1,hubbu%nuatm
          DO IATTOT=1,hubbu%nuproj
            DO IG=1,ncpw%ngw
              DCA(IG,IATTOT,K)=DCA(IG,IATTOT,K)-0.5D0*CA(IG,IATTOT)*&
                     DCA0SCA0(IATTOT,IUATM,K)/CA0SCA0(IATTOT)
            END DO
          END DO
        END DO
      END DO
!     tt4e=TIMEF()
! ---------------------------------------------------------------------
!
!  
!     tt5=TIMEF()
      DO K=1,3
        CALL OVLAP2(ncpw%ngw,NSTATE,hubbu%nuproj,C0XDCA(1,1,K),C0,DCA(1,1,K),.True.)
      END DO
!     tt5e=TIMEF()
!
      CALL mp_sum(C0XDCA,NSTATE*hubbu%nuproj*3,parai%allgrp)
!
      CALL zeroing(C0XB)!,hubbu%nuatm*maxsys%nhxs*NSTATE)
      CALL zeroing(CAXB)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj)
      CALL zeroing(C0XDB)!,hubbu%nuatm*maxsys%nhxs*NSTATE*3)
      CALL zeroing(CAXDB)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3)
      CALL zeroing(DCAXB)!,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3)

!     tt6=TIMEF()
      DO IUATM=1,hubbu%nuatm
        IAT=hubbu%uatm(IUATM)
        IS=IATPT(2,IAT)
        IA=IATPT(1,IAT)
        DO IV=1,nlps_com%ngh(IS)
          CIL=DCMPLX(0.D0,-1.D0)**NGHTOL(IV,IS) 
          DO IG=1,ncpw%ngw
            BF=TWNL(IG,IV,IS,1)   
            IBE=CIL*EIGR(IG,IAT,1)*DCMPLX(BF,0.D0) 
            MYBETA(IG)=IBE        
!
            DO K=1,3
              GV=GK(K,IG)*parm%tpiba
              IIGBE=DCMPLX(0.D0,-1.D0)*DCMPLX(GV,0.D0)*IBE 
              MYDBETA(IG,K)=IIGBE 
            END DO
          END DO 
!
          IF(GEQ0)THEN
            MYBETA(1)=MYBETA(1)*DCMPLX(0.5D0,0.d0)
            DO K=1,3
              MYDBETA(1,K)=MYDBETA(1,K)*DCMPLX(0.5D0,0.D0)
            END DO
          END IF
!
          DO ISTATE=1,NSTATE
            C0XB(IV,ISTATE,IUATM)=C0XB(IV,ISTATE,IUATM)+&
                  2.D0*DDOT(2*ncpw%ngw,MYBETA,1,C0(1,ISTATE),1)
            DO K=1,3
               C0XDB(IV,ISTATE,K,IUATM)=C0XDB(IV,ISTATE,K,IUATM)+&
                  2.D0*DDOT(2*ncpw%ngw,MYDBETA(1,K),1,C0(1,ISTATE),1)
            END DO
          END DO
!
          DO IATTOT=1,hubbu%nuproj 
            CAXB(IV,IATTOT,IUATM)=CAXB(IV,IATTOT,IUATM)+&
                2.D0*DDOT(2*ncpw%ngw,MYBETA,1,CA(1,IATTOT),1)
            DO K=1,3
              CAXDB(IV,IATTOT,K,IUATM)=CAXDB(IV,IATTOT,K,IUATM)+&
                2.D0*DDOT(2*ncpw%ngw,MYDBETA(1,K),1,CA(1,IATTOT),1)
              DCAXB(IV,IATTOT,K,IUATM)=DCAXB(IV,IATTOT,K,IUATM)+&
                2.D0*DDOT(2*ncpw%ngw,MYBETA,1,DCA(1,IATTOT,K),1)
            END DO
          END DO
!
        END DO
      END DO 
!     tt6e=TIMEF()
!
      CALL mp_sum(C0XB ,hubbu%nuatm*maxsys%nhxs*NSTATE  ,parai%allgrp)
      CALL mp_sum(C0XDB,hubbu%nuatm*maxsys%nhxs*NSTATE*3,parai%allgrp)
      CALL mp_sum(CAXB ,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj  ,parai%allgrp)
      CALL mp_sum(CAXDB,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)
      CALL mp_sum(DCAXB,hubbu%nuatm*maxsys%nhxs*hubbu%nuproj*3,parai%allgrp)

     deallocate(SCA,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(DCA0,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(DCA0SCA0,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(CA0SCA0,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(DCA0XCA0,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(CA0B,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(DCA0B,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 
     deallocate(CA0DB,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
             __LINE__,__FILE__) 

!
      FIRSTCALL=FIRSTCALL+1
!
      CALL TIHALT(procedureN,ISUB)
!     if (parent) print *, 'nsiemer Time HUBINTS_USPP2 =',tte-tt
!     if (parent) print *, ' loops:', tt1e-tt1, tt2e-tt2, tt3e-tt3, &
!                tt4e-tt4 , tt5e-tt5, tt6e,tt6
      END SUBROUTINE HUBINTS_USPP2
!
!     ====================================================================
      SUBROUTINE HUBINTS(NSTATE,C0,CA,C0XDCA,DCA)
!     ====================================================================
!     == Integrals needed for gradients of occupation matrix; non USPP  ==
!     ====================================================================
      IMPLICIT NONE
!
      INTEGER                       ::  NSTATE
!
      COMPLEX(real_8)               ::  C0(:,:),&
                                        DCA(ncpw%ngw,hubbu%nuproj,3),&
                                        CA(:,:)
!
      REAL(real_8)                  ::  C0XDCA(NSTATE,hubbu%nuproj,3)
!
      COMPLEX(real_8)               ::  IGV
      COMPLEX(real_8)               ::  ZZERO,ZONE
      COMPLEX(real_8),PARAMETER     ::  CI=(0.D0,-1.D0)
      CHARACTER(*), PARAMETER       ::  procedureN = 'HUBINTS'
!
      INTEGER                       ::  ISA,ISA0,IA,IS,IV,K,IATTOT,&
                                        ISTATE,IG,ISUB
      REAL(real_8)                  ::  BF,GV,DDOT
!
      CALL TISET(procedureN,ISUB)
!
      IF(tkpts%tkpnt)CALL STOPGM(procedureN,'DFT+U NOT IMPLEMENTED', __LINE__,&
        __FILE__)
!
! compute dphi = -i. G. phi
      CALL zeroing(dca)!,ncpw%ngw*hubbu%nuproj*3)
      CALL zeroing(C0XDCA)!,NSTATE*hubbu%nuproj*3)
      DO K=1,3
        DO IATTOT=1,hubbu%nuproj
           DO IG=1,ncpw%ngw
             IGV=DCMPLX(0.D0,-1.D0)*DCMPLX(GK(K,IG)*parm%tpiba,0.D0)
             DCA(IG,IATTOT,K)=CA(IG,IATTOT)*IGV
           END DO
        END DO
      END DO
!
      DO K=1,3
        CALL OVLAP2(ncpw%ngw,NSTATE,hubbu%nuproj,C0XDCA(1,1,K),C0,DCA(1,1,K),.True.)
      END DO
!
      CALL mp_sum(C0XDCA,NSTATE*hubbu%nuproj*3,parai%allgrp)
!
      CALL TIHALT(procedureN ,ISUB)
      END SUBROUTINE HUBINTS
END MODULE domdr_utils
