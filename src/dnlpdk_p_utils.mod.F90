MODULE dnlpdk_p_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE fitpack_utils,                   ONLY: curv2
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggng,&
                                             nsplpo,&
                                             twns
  USE response_pmod,                   ONLY: ddtwnl_ddk,&
                                             dtwnl_dk,&
                                             twnl_m,&
                                             twnl_p
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE ylmr2_utils,                     ONLY: ylmr2
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dnlps_dk
  PUBLIC :: give_scr_putwnl_kpert

CONTAINS

  ! ==================================================================
  SUBROUTINE dnlps_dk
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'dnlps_dk'

    INTEGER                                  :: ierr, is, isub, iv
    REAL(real_8)                             :: deltak

! ==--------------------------------------------------------------==
! ==  ALLOCATION OF MEMORY                                        ==
! ==--------------------------------------------------------------==

    CALL tiset('  DNLPS_DK',isub)

    IF (pslo_com%tivan) THEN
       CALL stopgm('DNLPS_DK','Vanderbilt not implemented',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculation of Maximal value of L
    maxsys%lpmax=0
    DO is=1,ions1%nsp
       IF (nlps_com%ngh(is).GT.0) THEN
          IF (dpot_mod%tkb(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,nghcom(iv,is))
          ELSEIF (sgpp1%tsgp(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,sgpp2%lpval(iv,is))
          ENDIF
       ENDIF
    ENDDO

    ALLOCATE(twnl_p(2*ncpw%ngw,maxsys%nhxs,maxsys%nsx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(twnl_p)!,3*maxsys%nsx*maxsys%nhxs*2*ngw)

    ALLOCATE(twnl_m(2*ncpw%ngw,maxsys%nhxs,maxsys%nsx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(twnl_m)!,3*maxsys%nsx*maxsys%nhxs*2*ngw)

    ! Computation of the form factors in G+k and G-k
    ! Saved in TWNL_P and TWNL_M
    ! Computed separatedly for Kx Ky and Kz for each G vector from 1:2*NGW
    ! DELTAK is the increment given in (2*PI/A) units

    deltak = 0.0050_real_8
    CALL putwnl_kpert(deltak)

    ! ==--------------------------------------------------------------==
    ! ==  NONLINEAR CORE CORRECTION                                   ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) CALL stopgm('DNLPS_DK','nonlinear core not implemented'&
         ,& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  FINITE TEMPERATURE ELECTRONS                                ==
    ! ==--------------------------------------------------------------==
    IF (fint1%ttrot) CALL stopgm('DNLPS_DK',&
         'finite electron temperature not implemented',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  INITIALIZE DERIVATIVE ARRAYS FOR STRESS CALCULATION         ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tprcp.OR.cntl%tpres) THEN
       CALL stopgm('DNLPS_DK','stress calculation not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DEALLOCATE ALL MEMORY NOT LONGER NEEDED                     ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('  DNLPS_DK')

    DEALLOCATE(twnl_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(twnl_m,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('  DNLPS_DK',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dnlps_dk
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE putwnl_kpert(deltak)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE TWNL_P(1:2*NGW,1:NGH(IS),1:NSP,3) [response_p.inc] ==
    ! == CALCULATE TWNL_M(1:2*NGW,1:NGH(IS),1:NSP,3) [response_p.inc] ==
    ! ==      Non-Local projectors array                              ==
    ! ==      for each G-components +/- K(I)(Kleinman-Bylander form)  ==
    ! == FOR TIVAN is not implemented                                 ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: deltak

    CHARACTER(*), PARAMETER                  :: procedureN = 'putwnl_kpert'

    INTEGER                                  :: i, idir, idir2, ierr, ig, is, &
                                                isub, iv, lp
    REAL(real_8)                             :: kplus(3), tw_m, tw_mm, tw_mp, &
                                                tw_p, tw_pm, tw_pp, u2dk, &
                                                udk2, vol, xx
    REAL(real_8), ALLOCATABLE :: gkrk(:,:), hgkm_m(:), hgkm_mm(:), hgkm_p(:), &
      hgkm_pm(:,:), hgkm_pp(:), hgkp_m(:), hgkp_mm(:), hgkp_p(:), &
      hgkp_pm(:,:), hgkp_pp(:), ylmb_m(:,:), ylmb_mm(:,:), ylmb_p(:,:), &
      ylmb_pm(:,:,:), ylmb_pp(:,:)
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset(' PUTWNL_KP',isub)
    ! TODO check stat
    ! TODO align for BG
    ALLOCATE(ylmb_p(2*ncpw%nhg, maxsys%lpmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ylmb_m(2*ncpw%nhg, maxsys%lpmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gkrk(3, ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(hgkp_p(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm_p(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkp_m(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm_m(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(hgkp_p)!,nhg)
    CALL zeroing(hgkm_p)!,nhg)
    CALL zeroing(hgkp_m)!,nhg)
    CALL zeroing(hgkm_m)!,nhg)

    CALL zeroing(twnl_p)!,3*2*ngw*maxsys%nhxs*maxsys%nsx)
    CALL zeroing(twnl_m)!,3*2*ngw*maxsys%nhxs*maxsys%nsx)

    ! ==--------------------------------------------------------------==
    vol=1._real_8/SQRT(parm%omega)
    u2dk = 1.0_real_8/(2.0_real_8*deltak*parm%tpiba)
    udk2 = 1.0_real_8/(deltak*deltak*parm%tpiba2)

    DO  idir = 1,3
       CALL zeroing(kplus)!,3)
       kplus(idir)=deltak

       CALL zeroing(ylmb_p)!,2*nhg*maxsys%lpmax)
       CALL zeroing(ylmb_m)!,2*nhg*maxsys%lpmax)

       CALL zeroing(gkrk)!,3*nhg)

       DO lp=1,maxsys%lpmax
          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)+kplus(i)
             ENDDO
             hgkp_p(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          CALL ylmr2(lp,ncpw%nhg,hgkp_p,gkrk,ylmb_p(1,lp))

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)+kplus(i)
             ENDDO
             hgkm_p(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          CALL ylmr2(lp,ncpw%nhg,hgkm_p,gkrk,ylmb_p(1+ncpw%nhg,lp))

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)-kplus(i)
             ENDDO
             hgkp_m(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          CALL ylmr2(lp,ncpw%nhg,hgkp_m,gkrk,ylmb_m(1,lp))

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)-kplus(i)
             ENDDO
             hgkm_m(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          CALL ylmr2(lp,ncpw%nhg,hgkm_m,gkrk,ylmb_m(1+ncpw%nhg,lp))
       ENDDO

       ! Form factors calculated in G+K and in G-K
       ! The K increment is considered in the three cartesian 
       ! directions separatedly
       ! for each species the product "radial part" times " angular part"
       ! is performed 

       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             IF (dpot_mod%tkb(is)) THEN
                ! AK 2005/03/26: unused
                ! L=NGHTOL(IV,IS)+1
                lp=nghcom(iv,is)
             ELSEIF (sgpp1%tsgp(is)) THEN
                lp=sgpp2%lpval(iv,is)
             ENDIF
             xx=dasum(nsplpo,twns(1,1,iv,is),1)
             IF (xx.GT.1.e-12_real_8) THEN
                DO ig=1,ncpw%ngw
                   tw_p=curv2(hgkp_p(ig),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                   tw_m=curv2(hgkp_m(ig),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)

                   twnl_p(ig,iv,is,idir)=ylmb_p(ig,lp)*tw_p*vol
                   twnl_m(ig,iv,is,idir)=ylmb_m(ig,lp)*tw_m*vol
                   ! FIRST  DERIVATIVE
                   dtwnl_dk(ig,iv,is,idir)= u2dk *&
                        (twnl_p(ig,iv,is,idir)-twnl_m(ig,iv,is,idir))
                   ! SECOND  DERIVATIVE

                   ddtwnl_ddk(ig,iv,is,idir,idir)= udk2 *&
                        (twnl_p(ig,iv,is,idir)+twnl_m(ig,iv,is,idir)&
                        -2.0_real_8*twnl(ig,iv,is,1))
                ENDDO

                DO ig=1,ncpw%ngw
                   tw_p=curv2(hgkm_p(ig),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                   tw_m=curv2(hgkm_m(ig),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)

                   twnl_p(ig+ncpw%ngw,iv,is,idir)=ylmb_p(ig+ncpw%nhg,lp)*&
                        tw_p*vol
                   twnl_m(ig+ncpw%ngw,iv,is,idir)=ylmb_m(ig+ncpw%nhg,lp)*&
                        tw_m*vol
                   ! FIRST  DERIVATIVE
                   dtwnl_dk(ig+ncpw%ngw,iv,is,idir)= u2dk *&
                        (twnl_p(ig+ncpw%ngw,iv,is,idir)-&
                        twnl_m(ig+ncpw%ngw,iv,is,idir))
                   ! SECOND DERIVATIVE
                   ddtwnl_ddk(ig+ncpw%ngw,iv,is,idir,idir)= udk2 *&
                        (twnl_p(ig+ncpw%ngw,iv,is,idir)+&
                        twnl_m(ig+ncpw%ngw,iv,is,idir)&
                        -(-1)**(nghtol(iv,is))&
                        *2.0_real_8*twnl(ig,iv,is,1))
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO                     ! END OF 1,IDIR

    ! ==--------------------------------------------------------------==
    ! ==  DEALLOCATE ALL MEMORY NOT LONGER NEEDED                     ==
    ! ==--------------------------------------------------------------==

    DEALLOCATE(hgkp_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkm_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkp_m,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkm_m,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! ==================================================================
    ! ===    MIXED SECOND DERIVATIVES                                 ==
    ! ==================================================================

    ALLOCATE(ylmb_pp(2*ncpw%nhg,maxsys%lpmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ylmb_mm(2*ncpw%nhg,maxsys%lpmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ylmb_pm(2*ncpw%nhg,maxsys%lpmax,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ALLOCATE(hgkp_pp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm_pp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkp_mm(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm_mm(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkp_pm(ncpw%nhg,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm_pm(ncpw%nhg,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ! Calculation of the following twrms
    ! PP : Gx+Delta,Gy+Delta,Gz - Gx+Delta,Gy,Gz+Delta - Gx,Gy+Delta,Gz+Delta 
    ! MM : Gx-Delta,Gy-Delta,Gz - Gx-Delta,Gy,Gz-Delta - Gx,Gy-Delta,Gz-Delta 
    ! PM : Gx+Delta,Gy-Delta,Gz - Gx+Delta,Gy,Gz-Delta - Gx,Gy+Delta,Gz-Delta 
    ! Gx-Delta,Gy+Delta,Gz - Gx-Delta,Gy,Gz+Delta - Gx,Gy-Delta,Gz+Delta 


    DO idir = 1,3
       DO idir2 = idir+1,3
          ! 
          CALL zeroing(ylmb_pp)!,2*nhg*maxsys%lpmax)
          CALL zeroing(ylmb_mm)!,2*nhg*maxsys%lpmax)
          CALL zeroing(ylmb_pm)!,2*nhg*maxsys%lpmax*2)

          CALL zeroing(hgkp_pp)!,nhg)
          CALL zeroing(hgkm_pp)!,nhg)
          CALL zeroing(hgkp_mm)!,nhg)
          CALL zeroing(hgkm_mm)!,nhg)
          CALL zeroing(hgkp_pm)!,2*nhg)
          CALL zeroing(hgkm_pm)!,2*nhg)

          CALL zeroing(gkrk)!,3*nhg)

          ! PP part
          CALL zeroing(kplus)!,3)
          kplus(idir)=deltak
          kplus(idir2)=deltak

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)+kplus(i)
             ENDDO
             hgkp_pp(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkp_pp,gkrk,ylmb_pp(1,lp))
          ENDDO

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)+kplus(i)
             ENDDO
             hgkm_pp(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkm_pp,gkrk,ylmb_pp(1+ncpw%nhg,lp))
          ENDDO

          ! MM part
          CALL zeroing(kplus)!,3)
          kplus(idir)=-deltak
          kplus(idir2)=-deltak

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)+kplus(i)
             ENDDO
             hgkp_mm(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkp_mm,gkrk,ylmb_mm(1,lp))
          ENDDO

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)+kplus(i)
             ENDDO
             hgkm_mm(ig) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkm_mm,gkrk,ylmb_mm(1+ncpw%nhg,lp))
          ENDDO

          ! PM part
          ! part1
          CALL zeroing(kplus)!,3)
          kplus(idir)=deltak
          kplus(idir2)=-deltak

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)+kplus(i)
             ENDDO
             hgkp_pm(ig,1) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkp_pm(1,1),gkrk,ylmb_pm(1,lp,1))
          ENDDO

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)+kplus(i)
             ENDDO
             hgkm_pm(ig,1) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkm_pm(1,1),gkrk,ylmb_pm(1+ncpw%nhg,lp,1))
          ENDDO

          ! part2
          CALL zeroing(kplus)!,3)
          kplus(idir)=-deltak
          kplus(idir2)=deltak

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=gk(i,ig)+kplus(i)
             ENDDO
             hgkp_pm(ig,2) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkp_pm(1,2),gkrk,ylmb_pm(1,lp,2))
          ENDDO

          DO ig = 1 , ncpw%nhg
             DO i=1,3
                gkrk(i,ig)=-gk(i,ig)+kplus(i)
             ENDDO
             hgkm_pm(ig,2) = gkrk(1,ig)*gkrk(1,ig)+&
                  gkrk(2,ig)*gkrk(2,ig)+&
                  gkrk(3,ig)*gkrk(3,ig)
          ENDDO
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hgkm_pm(1,2),gkrk,ylmb_pm(1+ncpw%nhg,lp,2))
          ENDDO

          ! for each species the product "radial part" times " angular part"
          ! is performed 

          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                IF (dpot_mod%tkb(is)) THEN
                   ! AK 2005/03/26: unused
                   ! L=NGHTOL(IV,IS)+1
                   lp=nghcom(iv,is)
                ELSEIF (sgpp1%tsgp(is)) THEN
                   lp=sgpp2%lpval(iv,is)
                ENDIF
                xx=dasum(nsplpo,twns(1,1,iv,is),1)
                IF (xx.GT.1.e-12_real_8) THEN
                   DO ig=1,ncpw%ngw
                      tw_pp=curv2(hgkp_pp(ig),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_mm=curv2(hgkp_mm(ig),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_pm=curv2(hgkp_pm(ig,1),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_mp=curv2(hgkp_pm(ig,2),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)

                      tw_pp=ylmb_pp(ig,lp)*tw_pp*vol
                      tw_mm=ylmb_mm(ig,lp)*tw_mm*vol
                      tw_pm=ylmb_pm(ig,lp,1)*tw_pm*vol
                      tw_mp=ylmb_pm(ig,lp,2)*tw_mp*vol

                      ! SECOND  DERIVATIVE
                      ddtwnl_ddk(ig,iv,is,idir,idir2)= udk2 *&
                           (tw_pp-tw_pm-tw_mp+tw_mm)
                      ddtwnl_ddk(ig,iv,is,idir2,idir)&
                           =ddtwnl_ddk(ig,iv,is,idir,idir2)
                   ENDDO

                   DO ig=1,ncpw%ngw
                      tw_pp=curv2(hgkm_pp(ig),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_mm=curv2(hgkm_mm(ig),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_pm=curv2(hgkm_pm(ig,1),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)
                      tw_mp=curv2(hgkm_pm(ig,2),nsplpo,ggng(1),&
                           twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)

                      tw_pp=ylmb_pp(ig+ncpw%nhg,lp)*tw_pp*vol
                      tw_mm=ylmb_mm(ig+ncpw%nhg,lp)*tw_mm*vol
                      tw_pm=ylmb_pm(ig+ncpw%nhg,lp,1)*tw_pm*vol
                      tw_mp=ylmb_pm(ig+ncpw%nhg,lp,2)*tw_mp*vol

                      ! SECOND DERIVATIVE
                      ddtwnl_ddk(ig+ncpw%ngw,iv,is,idir,idir2)= udk2 *&
                           (tw_pp-tw_pm-tw_mp+tw_mm)
                      ddtwnl_ddk(ig+ncpw%ngw,iv,is,idir2,idir)&
                           =ddtwnl_ddk(ig+ncpw%ngw,iv,is,idir,idir2)
                   ENDDO
                ENDIF
             ENDDO         ! IV
          ENDDO             ! IS
       ENDDO                 ! IDIR2
    ENDDO                     ! IDIR1

    DEALLOCATE(ylmb_pp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ylmb_mm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ylmb_pm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkp_pp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkm_pp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkp_mm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkm_mm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkp_pm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hgkm_pm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(ylmb_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ylmb_m,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gkrk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt(' PUTWNL_KP',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putwnl_kpert

  ! ====================================================================
  ! ====================================================================
  SUBROUTINE give_scr_putwnl_kpert(lputwnl,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lputwnl
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lylmb

    lputwnl= 3*ncpw%nhg
    lylmb=2*ncpw%nhg*maxsys%lpmax
    lputwnl=MAX(1,lputwnl+2*lylmb)
    tag='3*NHG+4*NHG*maxsys%lpmax'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_putwnl_kpert
  ! ==================================================================
  ! ==================================================================


END MODULE dnlpdk_p_utils
