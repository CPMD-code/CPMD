MODULE forces_prop_utils
  USE cnstfc_utils,                    ONLY: restfc
  USE cotr,                            ONLY: cotr007,&
                                             gsrate,&
                                             resval,&
                                             resval_dest
  USE efld,                            ONLY: textfld
  USE ehpsi_utils,                     ONLY: set_b2l
  USE ehrenfest_utils,                 ONLY: ehrenfest
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE k_updwf_utils,                   ONLY: k_updwf
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_cputime
  USE mp_interface,                    ONLY: mp_bcast
  USE nlforce_utils,                   ONLY: give_scr_nlforce
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: infi,&
                                             ropt_mod
  USE summat_utils,                    ONLY: give_scr_summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE td_input,                        ONLY: td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE updrho_utils,                    ONLY: updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE utils,                           ONLY: zclean_k
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forces_prop
  PUBLIC :: give_scr_forces_prop

CONTAINS

  ! ==================================================================
  SUBROUTINE forces_prop(nstate,c0,ct,c2,cr,sc0,cscr,vpp,eigv,&
       rhoe,psi,&
       tau0,velp,taui,fion,ifcalc,&
       irec,tfor,tinfo)
    ! ==--------------------------------------------------------------==
    ! EHR
    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c0(:,:,:), ct(nkpt%ngwk,*), c2(nkpt%ngwk,nstate), &
      cr(*), sc0(nkpt%ngwk,nstate), cscr(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), &
                                                eigv(nstate,1), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(*), &
                                                taui(*), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'forces_prop'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:), gde(:)
    INTEGER                                  :: freq_diag, i, ierr, ik, &
                                                il_auxc, il_ddia, il_gam, &
                                                infr, isub, lsummat, nhpsi, &
                                                update_pot
    INTEGER, SAVE                            :: icall = -1
    LOGICAL                                  :: calcrho
    REAL(real_8)                             :: detot, etot0, tcpu, thl(2), &
                                                time1
    REAL(real_8), ALLOCATABLE                :: ddia(:)

    CALL tiset(procedureN,isub)

    calcrho=.TRUE.
    freq_diag=1
    ik=1
    ! ==--------------------------------------------------------------==
    time1 = m_cputime()
    CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    CALL give_scr_summat(lsummat,tag,nstate)
    il_auxc=MAX(il_auxc,lsummat)
    il_auxc=MAX(il_auxc,nstate**2)
    ALLOCATE(auxc(il_auxc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddia(2*il_ddia),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! TODO check why IL_DDIA*2 ?
    ! =`=--------------------------------------------------------------==
    ener_com%ecnstr=0.0_real_8
    ropt_mod%convwf=.FALSE.
    etot0=0.0_real_8
    IF (cntl%tdiag.AND.cntl%tlanc) THEN
       CALL set_b2l()
    ENDIF
    IF (cntl%tmdeh) THEN
       ALLOCATE(gde(((nkpt%ngwk*nstate+8)*cnti%mdiis)/4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (.NOT.tinfo.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(T9,A)')&
            'INFR         GEMAX            ETOT         DETOT     TCPU'
    ENDIF
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    ! represtinate external field
    IF (td_prop%stextpot) textfld=.TRUE.
    ! ==================================================================
    ! ==                     EHRENFEST DYNAMICS                       ==
    ! ==================================================================
    ! 
    IF (cntl%tmdeh) THEN
       IF (infi.GT.0) THEN
          IF (paral%parent.AND.icall.LT.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,A)') 'OCCUPATION FOR PROP'
             IF (paral%io_parent)&
                  WRITE(6,'(13F5.1,/)') (crge%f(i,ik),i=1,nstate)
          ENDIF
          CALL ehrenfest(c0,c2,rhoe,psi(:,1),sc0,eigv)
       ENDIF
       icall=icall+1
    ENDIF
    ! 
    ! ==================================================================
    ! ==                      END OF MAIN  LOOP                       ==
    ! ==================================================================
    IF (paral%io_parent.AND.icall.EQ.0) THEN
       WRITE(6,'(/,1X,A)') 'OCCUPATION FOR FORCES COMPUTATION'
       WRITE(6,'(13F5.1,/)') (crge%f(i,ik),i=1,nstate)
    ENDIF
    ! 
    IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
    ! 
    IF (tfor) THEN
       IF (cntl%tdiag) THEN
          ! Diagonalization
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
          CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               nstate,.TRUE.,tinfo,ropt_mod%calste,infr,thl,nhpsi)
       ELSE
          IF (.NOT.tkpts%tkpnt) THEN
             CALL updwf(c0(:,:,1),c2(:,:),sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,.TRUE.)
          ELSEIF (tkpts%tkpnt.AND.(.NOT.cntl%tmdeh)) THEN
             update_pot=1
             CALL k_updwf(c0,c2,sc0,tau0,fion,cr,gde,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,update_pot)
          ELSE
             CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,.TRUE.)
          ENDIF
          ! Forces from geometrical restraints
          IF (paral%parent) THEN
             DO i=1,cotr007%mrestr
                resval(i)=resval(i)+gsrate(i)*dt_ions
                IF (resval_dest(i).NE.-999._real_8) THEN
                   IF (gsrate(i).GT.0._real_8.AND.resval(i).GT.resval_dest(i))&
                        resval(i)=resval_dest(i) ! increase
                   IF (gsrate(i).LT.0._real_8.AND.resval(i).LT.resval_dest(i))&
                        resval(i)=resval_dest(i) ! decrease
                ENDIF
             ENDDO
             CALL restfc(tau0,fion)
          ENDIF
          CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
       ENDIF
    ENDIF
200 CONTINUE
    IF (cntl%tmdeh) THEN
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE forces_prop
  ! ==================================================================
  SUBROUTINE give_scr_forces_prop(lforces_prop,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces_prop
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: ltscr

! ==--------------------------------------------------------------==

    lforces_prop=1
    CALL give_scr_updwf(lforces_prop,tag,nstate,tfor)
    IF (tfor) THEN
       ltscr=3*maxsys%nax*maxsys%nsx
    ELSE
       ltscr=0
    ENDIF
    lforces_prop=MAX(lforces_prop,ltscr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_forces_prop
  ! ==================================================================

END MODULE forces_prop_utils
