MODULE noforce_utils
  USE csmat_utils,                     ONLY: csmat
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hnlmat_utils,                    ONLY: hnlmat
  USE hubbardu,                        ONLY: c2u0
  USE hubbardu_utils,                  ONLY: add_hubbardu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE nlforce_utils,                   ONLY: give_scr_nlforce,&
                                             nlforce
  USE nlps,                            ONLY: imagp
  USE nlsl_utils,                      ONLY: nlsl
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE reigs_utils,                     ONLY: reigs
  USE rgs_utils,                       ONLY: uinv
  USE rnlfl_utils,                     ONLY: rnlfl
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rotate,&
                                             rottr
  USE rscpot_utils,                    ONLY: give_scr_rscpot,&
                                             rscpot
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE summat_utils,                    ONLY: summat
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noforce
  PUBLIC :: give_scr_noforce

CONTAINS

  ! ==================================================================
  SUBROUTINE noforce(c0,c2,sc0,tau0,fion,eigv,rhoe,psi,xmat1,xmat2,&
       nstate,tfor)
    ! ==--------------------------------------------------------------==
    ! ==          COMPUTES FOR  A SET OF NONORTHOGONAL ORBITALS       ==
    ! ==                     THE TOTAL ENERGY                         ==
    ! ==                  THE ELECTRONIC FORCES                       ==
    ! ==                  THE FORCES ON THE IONS                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                eigv(*), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: xmat2(nstate,*), &
                                                xmat1(nstate,*)
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'noforce'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: pab
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:)
    INTEGER                                  :: i, iabsl, ierr, il_auxc, &
                                                il_ddia, il_gam, il_smat, &
                                                isub, lnoforce
    INTEGER, EXTERNAL                        :: izamax
    LOGICAL                                  :: debug
    REAL(real_8), ALLOCATABLE                :: ddia(:), gam(:), smat(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('   NOFORCE',isub)
    debug=.FALSE.
    IF (lspin2%tlse) CALL stopgm('NOFORCE','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    pab = 0._real_8
    IF (imagp.EQ.2) CALL stopgm('NOFORCE','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL give_scr_noforce(lnoforce,il_gam,il_auxc,il_smat,il_ddia,tag,&
         nstate,tfor)
    ! ==--------------------------------------------------------------==
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    gnmax=0.0_real_8
    gnorm=0.0_real_8
    IF (pslo_com%tivan) CALL rnlsm(c0,nstate,1,1,.FALSE.)
    ! ==--------------------------------------------------------------==
    ALLOCATE(gam(il_gam),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(auxc(il_auxc/2+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ddia(il_ddia),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(smat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==   OVERLAP MATRIX                                             ==
    ! ==--------------------------------------------------------------==
    CALL csmat(gam,c0,fnl,nstate,1)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE INVERSE CHOLESKY DECOMPOSITION                   ==
    ! ==   SCHMIDT ORTHOGONALIZATION    SC0 = U**(-1)*C0              ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       CALL dcopy(nstate*nstate,gam,1,smat(1,1),1)
       IF (cntl%tlsd) THEN
          CALL uinv('U',smat(1,1),nstate,spin_mod%nsup)
          CALL uinv('U',smat(spin_mod%nsup+1,spin_mod%nsup+1),nstate,spin_mod%nsdown)
       ELSE
          CALL uinv('U',smat,nstate,nstate)
       ENDIF
       CALL dcopy(nstate*nstate,smat(1,1),1,xmat2(1,1),1)
    ENDIF
    CALL mp_bcast(smat,SIZE(smat),parai%source,parai%allgrp)
    CALL dcopy(2*ncpw%ngw*nstate,c0,1,sc0,1)
    CALL rottr(1._real_8,sc0,smat,"N",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    IF (geq0) CALL zclean(sc0,nstate,ncpw%ngw)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE POTENTIAL AND THE FORCE ON THE IONS          ==
    ! ==--------------------------------------------------------------==
    CALL rnlsm(sc0,nstate,1,1,tfor)
    CALL rscpot(sc0,tau0,fion,rhoe,psi,tfor,ropt_mod%calste,nstate,1)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE ELECTRONIC FORCE                             ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(c2)!,ngw*nstate)
    ! ==--------------------------------------------------------------==
    ! == Compute the force on the electronic degrees of freedom due   ==
    ! == to the local potential (stored in RHOE)                      ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(gam(il_gam),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(auxc(il_auxc/2+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ddia(il_ddia),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(smat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    CALL vpsi(sc0,c2,crge%f(:,1),rhoe,psi(:,1),nstate,1,clsd%nlsd,.TRUE.)
    ! c2u0 is calculated in uprho or rscpot
    IF (cntl%thubb) CALL add_hubbardu(c2,c2u0,nstate)
    IF (pslo_com%tivan) THEN
       CALL ovlap(nstate,gam,c2,sc0)
       CALL hnlmat(gam,crge%f,nstate)
       CALL summat(gam,nstate)
       CALL rotate(-1.0_real_8,sc0,1.0_real_8,c2,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
            spin_mod%nsdown)
       IF (tfor) CALL rnlfl(fion,gam,nstate,1)
       IF (ropt_mod%calste) CALL nlsl(gam,nstate)
       CALL nlforce(c2,crge%f,gam,auxc,ddia,nstate)
    ELSE
       ! ==--------------------------------------------------------------==
       ! == Compute the force on the electronic degrees of freedom due   ==
       ! == to the non-local part of the potential, and add it to the    ==
       ! == other piece, coming from the local contribution.             ==
       ! ==--------------------------------------------------------------==
       CALL fnonloc(c2,crge%f,nstate,1,clsd%nlsd,.TRUE.)
       CALL ovlap(nstate,gam,c2,sc0)
       CALL summat(gam,nstate)
       ! ==--------------------------------------------------------------==
       ! ==   C2(I) = C2(I) - SUM(J) <SC(I) | H | SC(J)> SC(J)           ==
       ! ==--------------------------------------------------------------==
       CALL rotate(-1.0_real_8,sc0,1.0_real_8,c2,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
            spin_mod%nsdown)
    ENDIF
    IF (ropt_mod%prteig.AND.paral%parent) THEN
       CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
       CALL reigs(nstate,gam,crge%f)
       CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==   ROTATE ELECTRONIC FORCE BACK INTO NONORTHOGONAL BASIS      ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL dcopy(nstate*nstate,xmat2(1,1),1,smat(1,1),1)
    CALL mp_bcast(smat,SIZE(smat),parai%source,parai%allgrp)
    CALL rottr(1._real_8,c2,smat,"T",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    IF (geq0) CALL zclean(c2,nstate,ncpw%ngw)
    ! ==--------------------------------------------------------------==
    gemax=0.0_real_8
    cnorm=0.0_real_8
    DO i=1,nstate
       iabsl=izamax(ncpw%ngw,c2(1,i),1)
       gemax=MAX(crge%f(i,1)*ABS(c2(iabsl,i)),gemax)
       cnorm=cnorm+crge%f(i,1)*dotp(ncpw%ngw,c2(:,i),c2(:,i))
    ENDDO
    CALL mp_sum(cnorm,parai%allgrp)
    CALL mp_max(gemax,parai%allgrp)
    cnorm=SQRT(cnorm/(nstate*spar%ngws))
    ! ==--------------------------------------------------------------==
    IF (tfor) THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
       IF (paral%parent) THEN
          CALL symvec(fion)
          CALL taucl(fion)
          CALL gsize(fion,gnmax,gnorm)
       ENDIF
    ENDIF
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('   NOFORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noforce
  ! ==================================================================
  SUBROUTINE give_scr_noforce(lnoforce,il_gam,il_auxc,il_smat,&
       il_ddia,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnoforce, il_gam, il_auxc, &
                                                il_smat, il_ddia
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lrnlsm, lrscpot

    il_smat=0
    IF (pslo_com%tivan) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       il_gam=imagp*nstate*nstate
    ENDIF
    IF (ropt_mod%prteig.AND.paral%parent) THEN
       il_auxc=MAX(il_auxc,nstate*(nstate+1)/2+3*nstate)! REIGS
    ENDIF
    il_auxc=MAX(il_auxc,nstate*nstate)
    il_smat=MAX(il_smat,nstate*nstate)
    il_ddia=MAX(il_ddia,2*maxsys%nax*nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    lnoforce=il_gam+il_auxc+il_smat+il_ddia+100
    lnoforce=MAX(lnoforce,lrnlsm,lrscpot)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_noforce
  ! ==================================================================

END MODULE noforce_utils
