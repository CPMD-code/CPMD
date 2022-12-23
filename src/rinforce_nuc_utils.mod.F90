MODULE rinforce_nuc_utils
  USE aainit_utils,                    ONLY: aainit
  USE aavan,                           ONLY: indv,&
                                             lpl
  USE atom,                            ONLY: gnl,&
                                             rps,&
                                             rw,&
                                             vr
  USE cnst,                            ONLY: pi
  USE cppt,                            ONLY: qrad,&
                                             twnl,&
                                             ylmb
  USE cvan,                            ONLY: deeq,&
                                             dvan,&
                                             nelev,&
                                             qq
  USE dpot,                            ONLY: dpot_mod
  USE eam,                             ONLY: tieam
  USE eam_pot_utils,                   ONLY: eamin
  USE ehpsi_utils,                     ONLY: eg2
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE formf_utils,                     ONLY: formfn,&
                                             formfv,&
                                             formsg
  USE ions,                            ONLY: al,&
                                             bl,&
                                             ions0,&
                                             ions1,&
                                             maxgau,&
                                             rcl
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mm_input,                        ONLY: lqmmm
  USE mm_ion_dens,                     ONLY: mm_rhops
  USE nlcc,                            ONLY: corecg,&
                                             corel,&
                                             drhoc,&
                                             rcgrid,&
                                             rhoc
  USE nlccset_utils,                   ONLY: nlccset
  USE nlps,                            ONLY: nghcom,&
                                             nlps_com
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE qrada_s_utils,                   ONLY: qrada_s
  USE qspl,                            ONLY: ggng,&
                                             ggnh,&
                                             nqdim,&
                                             nsplpo,&
                                             qspl1
  USE ragg,                            ONLY: raggio
  USE rinforce_utils,                  ONLY: give_scr_putwnl,&
                                             putps,&
                                             putwnl,&
                                             setspline
  USE rnlin_utils,                     ONLY: nlin
  USE rnlset_utils,                    ONLY: rnlset
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE str2,                            ONLY: qrada
  USE system,                          ONLY: &
       cntl, lmaxx, lx, maxsys, mx, nbrx, ncpw, nhx, nkpt, nlx
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdbinit_utils,                   ONLY: qinit,&
                                             vdbinit
  USE vdbp,                            ONLY: &
       betar, dion, ncpr1, qfunc, qqq, qrl, r, rab, rsatom, rscore, ru, rucore
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rinforce_nuc

CONTAINS

  ! ==================================================================
  SUBROUTINE rinforce_nuc(ndummies)
    ! ==--------------------------------------------------------------==

    ! Variables

    INTEGER                                  :: ndummies

    CHARACTER(*), PARAMETER                  :: procedureN = 'rinforce_nuc'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ierr, im, is, istep, isub, &
                                                iv, jv, lmaxv, lp, lscr, &
                                                lval, mrscr, nylmb
    REAL(real_8), ALLOCATABLE                :: rs1(:), rs2(:), rs3(:), scr(:)

    CALL tiset(procedureN,isub)

    ALLOCATE(rcgrid(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(corecg(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(al(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bl(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rcl(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Real space function and grids for Kleinman-Bylander + Potentials
    ALLOCATE(gnl(maxsys%mmaxx,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rps(maxsys%mmaxx,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rw(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vr(maxsys%mmaxx,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! CALL AZZERO(RCGRID,maxsys%mmaxx*maxsys%nsx)
    ! CALL AZZERO(CORECG,maxsys%mmaxx*maxsys%nsx)
    ! call azzero(al,maxsys%nsx*LMAXX*MAXGAU)
    ! call azzero(bl,maxsys%nsx*LMAXX*MAXGAU)
    ! call azzero(rcl,maxsys%nsx*LMAXX*MAXGAU)
    ! call azzero(gnl,maxsys%mmaxx*maxsys%nsx*LMAXX)
    ! call azzero(rps,maxsys%mmaxx*maxsys%nsx*LMAXX)
    ! call azzero(rw,maxsys%mmaxx*maxsys%nsx)
    ! call azzero(vr,maxsys%mmaxx*maxsys%nsx*LMAXX)

    ! free memory for arrays in setspline routine
    DEALLOCATE(ggnh,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ggng,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    ! ==  SELF-ENERGY (NOTHING TO DO WITH MANY-BODY SELF-ENERGY)      ==
    ! ==--------------------------------------------------------------==
    ener_com%eself=0._real_8
    DO is=1,ndummies ! NSP
       ener_com%eself=ener_com%eself+REAL(ions0%na(is),kind=real_8)*ions0%zv(is)*ions0%zv(is)/raggio(is)
    ENDDO
    ener_com%eself=ener_com%eself/SQRT(2._real_8*pi)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOCATION OF MEMORY                                        ==
    ! ==--------------------------------------------------------------==
    CALL setspline
    ! CALL MEMORY(IP_VPS,maxsys%nsx*NHG,'VPS')
    ! CALL MEMORY(IP_RHOPS,maxsys%nsx*NHG,'RHOPS')
    IF (lqmmm%qmmm)  THEN
       ALLOCATE(mm_RHOPS(maxsys%nsx,ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (pslo_com%tivan) THEN
       ALLOCATE(indv(nhx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nqdim    = ncpw%nhgl
       IF (qspl1%qspline) nqdim=2*nsplpo
       ALLOCATE(qrad(nqdim,nbrx,nbrx,lx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! TODO what to do with INDV and QRAD ?? 
    ENDIF
    IF (corel%tinlc) THEN
       ALLOCATE(rhoc(ncpw%nhg,ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (cntl%tpres.OR.cntl%tprcp)  THEN
          ALLOCATE(drhoc(ncpw%nhg,maxsys%nsx,6),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! SCRATCH SPACE ALLOCATION
    mrscr   = MAX(8*maxsys%mmaxx,nsplpo,nhx*nhx)
    ALLOCATE(rs1(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rs2(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rs3(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! SPACE ALLOCATION FOR SPLINE

    ! CALL MEMORY(IP_VOO,2*NSPLPO*maxsys%nsx,'VOO')
    ! CALL MEMORY(IP_TWNS,2*NSPLPO*maxsys%nsx*maxsys%nhxs,'TWNS')

    ! CALL AZZERO(VOO,2*NSPLPO*maxsys%nsx)
    ! CALL AZZERO(TWNS,2*NSPLPO*maxsys%nsx*maxsys%nhxs)

    ! INITIALIZATION FOR VANDERBILT PSEUDOPOTENTIALS
    IF (pslo_com%tivan) THEN
       lmaxv=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             lval=ncpr1%nvales(is)
             IF (lval.GT.3) CALL stopgm('RINFORCE',' L>3, L=',& 
                  __LINE__,__FILE__)
             IF (lval.GT.lmaxv) lmaxv=lval
          ENDIF
       ENDDO
       ! INITIALIZE CLEBSCH-GORDON COEFFICIENTS
       CALL aainit(lmaxv)
       ALLOCATE(nelev(maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(qq(maxsys%nhxs,maxsys%nhxs,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dvan(maxsys%nhxs,maxsys%nhxs,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(deeq(ions1%nat,maxsys%nhxs,maxsys%nhxs,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(deeq)
    ENDIF
    ! ==--------------------------------------------------------------==
    DO is=1,ndummies   ! NSP
       IF (pslo_com%tvan(is)) THEN
          CALL formfv(is,rs1,rs2)
       ELSEIF (sgpp1%tsgp(is)) THEN
          CALL formsg(is,rs1,rs2)
       ELSE
          IF (ions0%igau(is).GE.0) CALL formfn(is,rs1,rs2)
       ENDIF
    ENDDO
    CALL putps
    DO is=1,ndummies ! NSP
       IF (pslo_com%tvan(is)) THEN
          CALL vdbinit(is,rs1,rs2)
       ELSE
          CALL rnlset(is,rs1,rs2,rs3)
          CALL nlin(is,rs1,rs2)
       ENDIF
    ENDDO
    IF (pslo_com%tivan) CALL qinit(rs1,rs2)
    ! ==--------------------------------------------------------------==
    ! Calculation of Maximal value of L
    maxsys%lpmax=0
    IF (pslo_com%tivan) THEN
       DO im=1,mx
          DO iv=1,nlx
             DO jv=iv,nlx
                lp=lpl(iv,jv,im)
                maxsys%lpmax=MAX(maxsys%lpmax,lp)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DO is=1, ndummies  ! NSP
       IF (nlps_com%ngh(is).GT.0) THEN
          IF (dpot_mod%tkb(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,nghcom(iv,is))
          ELSEIF (sgpp1%tsgp(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,sgpp2%lpval(iv,is))
          ELSEIF (.NOT.pslo_com%tvan(is)) THEN
             istep=NINT(nlps_com%rmaxn(is))
             DO iv=1,nlps_com%ngh(is)
                lp=(iv-1)/istep+1
                maxsys%lpmax=MAX(maxsys%lpmax,lp)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    IF (pslo_com%tivan) THEN
       nylmb=nkpt%nhgk*maxsys%lpmax*nkpt%nkpnt
       ALLOCATE(ylmb(nkpt%nhgk,maxsys%lpmax,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ylmb)!,nylmb)
    ELSE
       nylmb=0
    ENDIF

    ! CALL MEMORY(IP_TWNL,maxsys%nsx*maxsys%nhxs*NGWK*NKPTALL,'TWNL')

    CALL zeroing(twnl)!,maxsys%nsx*maxsys%nhxs*nkpt%ngwk*kpts_com%nkptall)
    CALL give_scr_putwnl(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL putwnl
    ! ==--------------------------------------------------------------==
    ! ==  NONLINEAR CORE CORRECTION                                   ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) CALL nlccset
    ! ==--------------------------------------------------------------==
    ! ==  FINITE TEMPERATURE ELECTRONS                                ==
    ! ==--------------------------------------------------------------==
    IF (fint1%ttrot) CALL eg2
    ! ==--------------------------------------------------------------==
    ! ==  EAM POTENTIALS                                              ==
    ! ==--------------------------------------------------------------==
    tieam=.FALSE.
    DO is=1,ndummies  ! NSP
       tieam=tieam.OR.dpot_mod%team(is)
    ENDDO
    IF (tieam) CALL eamin
    ! ==--------------------------------------------------------------==
    ! ==  INITIALIZE DERIVATIVE ARRAYS FOR STRESS CALCULATION         ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tprcp.OR.cntl%tpres) THEN
       IF (pslo_com%tivan) THEN
          nqdim    = ncpw%nhgl
          IF (qspl1%qspline) nqdim=2*nsplpo
          ALLOCATE(qrada(nqdim,nbrx,nbrx,lx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! TODO what to do with qrada?? 
       ENDIF
       IF (pslo_com%tivan) CALL qrada_s(rs1,rs2)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DEALLOCATE ALL MEMORY NOT LONGER NEEDED                     ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('  RINFORCE')
    DEALLOCATE(rs1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rs2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rs3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(rcgrid,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(corecg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(al,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rcl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! call freem(ip_rv)

    ! call freem(ip_vps)
    ! call freem(ip_rhops)
    ! call freem(ip_voo)
    ! call freem(ip_twns)
    ! call freem(ip_twnl)

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (tkpts%tkpnt) DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       DEALLOCATE(rscore,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(dion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(betar,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qqq,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qfunc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qrl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rucore,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ru,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rsatom,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(r,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinforce_nuc

END MODULE rinforce_nuc_utils
