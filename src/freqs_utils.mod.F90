MODULE freqs_utils
  USE cppt,                            ONLY: hg,&
                                             inyh,&
                                             twnl
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE harm,                            ONLY: dcgt,&
                                             dsgtw,&
                                             dtan2w,&
                                             freq,&
                                             ngharm,&
                                             wdsgt,&
                                             xmu
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rggen_utils,                     ONLY: recips
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE tpar,                            ONLY: dtb2me

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: freqs

CONTAINS

  ! ==================================================================
  SUBROUTINE freqs (nstate,tfullh)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfullh

    CHARACTER(*), PARAMETER                  :: procedureN = 'freqs'

    INTEGER                                  :: ia, ierr, ig, ir, iri1, iri2, &
                                                iri3, is, isa0, iv, jv, ki, &
                                                kj, l, l2, li, lj, nh1, nh2, &
                                                nh3
    REAL(real_8)                             :: ar1(3), ar2(3), ar3(3), &
                                                br1(3), br2(3), br3(3), dd, &
                                                gkr1, gkr2, gkr3, hgr
    REAL(real_8), ALLOCATABLE                :: hamd(:)

! Variables
! TODO make sure that 
! NKPNT==1 => TWNL(NGW,maxsys%nhxs,*,1), F(NSTATE,1)
! 
! ==--------------------------------------------------------------==
! ==  Calculate frequencies for harmonic reference system         ==
! ==  Calculate scaled electron masses                            ==
! ==--------------------------------------------------------------==

    IF (.NOT.(cntl%tharm.OR.cntl%tmass)) RETURN
    ! ==--------------------------------------------------------------==
    IF (cntl%tharm) THEN
       ALLOCATE(freq(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dcgt(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dsgtw(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wdsgt(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dtan2w(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tmass) THEN
       ALLOCATE(xmu(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(vploc,parai%igeq0,parai%allgrp)
    ALLOCATE(hamd(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ngharm=0
    IF (tfullh) THEN
       DO ig=1,ncpw%ngw
          hamd(ig) = 0.5_real_8*parm%tpiba2*hg(ig)+vploc
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      dd=dvan(iv,jv,is)
                      DO ia=1,ions0%na(is)
                         dd=dd+deeq(isa0+ia,iv,jv,1)
                      ENDDO
                      IF (cntl%tlsd) THEN
                         DO ia=1,ions0%na(is)
                            dd=dd+deeq(isa0+ia,iv,jv,2)
                         ENDDO
                         dd=0.5_real_8*dd
                      ENDIF
                      hamd(ig)=hamd(ig)+dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                   ENDDO
                ENDDO
             ELSEIF (sgpp1%tsgp(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         dd=sgpp2%hlsg(ki,kj,l,is)*ions0%na(is)
                         hamd(ig)=hamd(ig)+&
                              dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   hamd(ig)=hamd(ig)+&
                        ions0%na(is)*wsg(is,iv)*twnl(ig,iv,is,1)*twnl(ig,iv,is,1)
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
          IF (hamd(ig).LT.cntr%hthrs) ngharm=ngharm+1
       ENDDO
    ELSE
       DO ir=1,3
          ar1(ir) = metr_com%ht0(1,ir)
          ar2(ir) = metr_com%ht0(2,ir)
          ar3(ir) = metr_com%ht0(3,ir)
       ENDDO
       CALL recips(parm%alat,ar1,ar2,ar3,br1,br2,br3)
       nh1=spar%nr1s/2+1
       nh2=spar%nr2s/2+1
       nh3=spar%nr3s/2+1
       DO ig=1,ncpw%ngw
          iri1=inyh(1,ig)-nh1
          iri2=inyh(2,ig)-nh2
          iri3=inyh(3,ig)-nh3
          gkr1=iri1*br1(1)+iri2*br2(1)+iri3*br3(1)
          gkr2=iri1*br1(2)+iri2*br2(2)+iri3*br3(2)
          gkr3=iri1*br1(3)+iri2*br2(3)+iri3*br3(3)
          hgr=gkr1**2+gkr2**2+gkr3**2
          hamd(ig) = 0.5_real_8*parm%tpiba2*hgr
          IF (hamd(ig).LT.cntr%hthrs) ngharm=ngharm+1
       ENDDO
    ENDIF
    DO ig=1,ngharm
       hamd(ig)=cntr%hthrs
    ENDDO
    IF (cntl%tmass) THEN
       DO ig=1,ncpw%ngw
          xmu(ig)=hamd(ig)*cntr%emass/cntr%hthrs
       ENDDO
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') ' FREQS| SHELL DEPENDENT ELECTRON MASSES'
          WRITE(6,'(A,F12.3)') ' FREQS| SHELL    1  ELECTRON MASS  ',&
               xmu(1)
          WRITE(6,'(A,I5,A,F12.3)') ' FREQS| SHELL',ncpw%ngw,&
               '  ELECTRON MASS  ',xmu(ncpw%ngw)
       ENDIF
    ENDIF
    IF (cntl%tharm) THEN
       DO ig=1,ngharm
          freq(ig)    = 0.0_real_8
          dcgt(ig)    = 1.0_real_8
          dsgtw(ig)   = cntr%delt_elec
          wdsgt(ig)   = 0.0_real_8
          dtan2w(ig)  = dtb2me
       ENDDO
       IF (cntl%tmass) THEN
          DO ig=ngharm+1,ncpw%ngw
             freq(ig) = SQRT(crge%f(1,1)*hamd(ig)/xmu(ig))
          ENDDO
          DO ig=ngharm+1,ncpw%ngw
             dcgt(ig)   = COS(freq(ig)*cntr%delt_elec)
             dsgtw(ig)  = SIN(freq(ig)*cntr%delt_elec)/freq(ig)
             wdsgt(ig)  = SIN(freq(ig)*cntr%delt_elec)*freq(ig)
             dtan2w(ig) = dtan(0.5_real_8*freq(ig)*cntr%delt_elec)/&
                  (freq(ig)*xmu(ig))
          ENDDO
       ELSE
          DO ig=ngharm+1,ncpw%ngw
             freq(ig) = SQRT(crge%f(1,1)*hamd(ig)/cntr%emass)
          ENDDO
          DO ig=ngharm+1,ncpw%ngw
             dcgt(ig)   = COS(freq(ig)*cntr%delt_elec)
             dsgtw(ig)  = SIN(freq(ig)*cntr%delt_elec)/freq(ig)
             wdsgt(ig)  = SIN(freq(ig)*cntr%delt_elec)*freq(ig)
             dtan2w(ig) = dtan(0.5_real_8*freq(ig)*cntr%delt_elec)/(freq(ig)*cntr%emass)
          ENDDO
       ENDIF
    ENDIF
    IF (paral%io_parent)WRITE(6,'(A,I5)') ' FREQS| SIZE OF HARMONIC REFERENC ESYSTEM:',ncpw%ngw-ngharm
    DEALLOCATE(hamd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE freqs
  ! ==================================================================

END MODULE freqs_utils
