MODULE getgyr_utils
  USE cnst,                            ONLY: factem
  USE ekinpp_utils,                    ONLY: ekinpp
  USE glemod,                          ONLY: glepar
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE pbc_utils,                       ONLY: pbc
  USE pimd,                            ONLY: &
       avcor, avcora, avgyr, avgyra, avsus, avsusa, fionks, fpcor, fpgyr, &
       fpsus, ipp1, pimd3, pmar, rcor, rcora, rgyr, rgyra, rsus, rsusa
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getgyration

CONTAINS

  ! ==================================================================
  SUBROUTINE getgyration(taup,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taup(:,:,:,:), velp(:,:,:,:)

    INTEGER                                  :: ia, ip, is, itau, k, nphalf
    REAL(real_8)                             :: ekinp1, eviri, rclass, &
                                                rmean(3), rnorm, rnp, rtot, &
                                                tempp1, var, vecpro, xlm, &
                                                ylm, zlm, xlm_,ylm_,zlm_

! Variables
! ==--------------------------------------------------------------==
! ==  RADII OF GYRATION OF DIFFERENT IONS                         ==
! ==    GILLAN, IN COMPUTER MODELLING OF FLUIDS, SOLIDS, AND      ==
! ==            POLYMERS, 1990                                    ==
! ==--------------------------------------------------------------==

    rnp=REAL(pimd3%np_total,kind=real_8)
    CALL zeroing(avgyra)!,maxsp*2)
    CALL zeroing(rgyra)!,maxsys%nsx*maxsys%nax*2)

    ! If we're doing PI with NVE we need to use the instantaneous
    ! temperature here. For details, see:
    ! http://cpmd.org/mailman/htdig/cpmd-list/2008-April/004350.html
    IF (cntl%tnosep.OR.glepar%gle_mode>0) THEN
       tempp1=cntr%tempw
    ELSE
       CALL ekinpp(ekinp1,velp(:,:,:,1))
       tempp1=ekinp1*factem*2.0_real_8/glib
    ENDIF

    DO is=1,ions1%nsp
       ! free particle reference
       fpgyr(is)=SQRT(factem/(4.0_real_8*pmar(is)*tempp1))
       avgyr(is)=0.0_real_8
       DO ia=1,ions0%na(is)
          DO k=1,3
             rmean(k)=0.0_real_8
             DO ip=1,pimd3%np_total
                rmean(k)=rmean(k)+taup(k,ia,is,ip)
             ENDDO
             rmean(k)=rmean(k)/rnp
          ENDDO
          var=0.0_real_8
          DO ip=1,pimd3%np_total
             xlm_=taup(1,ia,is,ip)-rmean(1)
             ylm_=taup(2,ia,is,ip)-rmean(2)
             zlm_=taup(3,ia,is,ip)-rmean(3)
             CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             var=var+xlm**2+ylm**2+zlm**2
          ENDDO
          ! per atom 
          IF (pimd3%np_total.EQ.1) THEN
             rgyr(ia,is)=0._real_8
          ELSE
             rgyr(ia,is)=SQRT(var/rnp)
          ENDIF
          ! accumulate first and second moment 
          rgyra(ia,is,1)=rgyra(ia,is,1)+rgyr(ia,is)
          rgyra(ia,is,2)=rgyra(ia,is,2)+rgyr(ia,is)**2
          ! per species
          avgyr(is)=avgyr(is)+rgyr(ia,is)
       ENDDO
       avgyr(is)=avgyr(is)/REAL(ions0%na(is),kind=real_8)
       ! accumulate first and second moment 
       avgyra(is,1)=avgyra(is,1)+avgyr(is)
       avgyra(is,2)=avgyra(is,2)+avgyr(is)**2
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  SUSCEPTIBILITY RADII OF DIFFERENT IONS                      ==
    ! ==    PARRINELLO,RAHMAN JCP 80, 860 (1984)                      ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(avsusa)!,maxsp*2)
    CALL zeroing(rsusa)!,maxsys%nsx*maxsys%nax*2)
    DO is=1,ions1%nsp
       ! free particle reference
       fpsus(is)=SQRT(factem/(2.0_real_8*pmar(is)*tempp1))
       avsus(is)=0.0_real_8
       DO ia=1,ions0%na(is)
          ! vector product
          rmean(1)=0.0_real_8
          rmean(2)=0.0_real_8
          rmean(3)=0.0_real_8
          DO ip=1,pimd3%np_total
             rmean(1)=rmean(1)+( taup(2,ia,is,ip)*taup(3,ia,is,ipp1(ip))&
                  -taup(3,ia,is,ip)*taup(2,ia,is,ipp1(ip)))
             rmean(2)=rmean(2)+( taup(1,ia,is,ip)*taup(3,ia,is,ipp1(ip))&
                  -taup(3,ia,is,ip)*taup(1,ia,is,ipp1(ip)))
             rmean(3)=rmean(3)+( taup(1,ia,is,ip)*taup(2,ia,is,ipp1(ip))&
                  -taup(2,ia,is,ip)*taup(1,ia,is,ipp1(ip)))
          ENDDO
          vecpro=rmean(1)**2+rmean(2)**2+rmean(3)**2
          vecpro=vecpro*(tempp1*pmar(is))**2/(factem*factem)
          ! spring kinetic contribution (virial representation)
          eviri=0.0_real_8
          DO k=1,3
             rmean(k)=0.0_real_8
             DO ip=1,pimd3%np_total
                rmean(k)=rmean(k)+taup(k,ia,is,ip)
             ENDDO
             rmean(k)=rmean(k)/rnp
          ENDDO
          DO ip=1,pimd3%np_total
             xlm_=taup(1,ia,is,ip)-rmean(1)
             ylm_=taup(2,ia,is,ip)-rmean(2)
             zlm_=taup(3,ia,is,ip)-rmean(3)
             CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             eviri=eviri-xlm*fionks(1,ia,is,ip)&
                  -ylm*fionks(2,ia,is,ip)&
                  -zlm*fionks(3,ia,is,ip)
          ENDDO
          eviri=eviri*2.0_real_8*factem/(3.0_real_8*rnp*rnp*tempp1)
          ! classical term
          rclass=(3.0_real_8*rnp-2.0_real_8)/(rnp**2)
          ! total estimator 
          rtot=rclass-eviri+vecpro
          ! be careful with crazy start configurations
          IF (rtot.LE.0.0_real_8) rtot=0.0_real_8
          ! per atom
          rsus(ia,is)=fpsus(is)*SQRT(rtot)
          ! accumulate first and second moment 
          rsusa(ia,is,1)=rsusa(ia,is,1)+rsus(ia,is)
          rsusa(ia,is,2)=rsusa(ia,is,2)+rsus(ia,is)**2
          ! per species
          avsus(is)=avsus(is)+rsus(ia,is)
       ENDDO
       avsus(is)=avsus(is)/REAL(ions0%na(is),kind=real_8)
       ! accumulate first and second moment  
       avsusa(is,1)=avsusa(is,1)+avsus(is)
       avsusa(is,2)=avsusa(is,2)+avsus(is)**2
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  RADII FROM IMAGINARY TIME CORRELATION FUNCTION              ==
    ! ==    CHANDLER, LES HOUCHES 1989                                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(avcora)!,maxsp*2)
    CALL zeroing(rcora)!,maxsys%nsx*maxsys%nax*2)
    nphalf=pimd3%np_total/2
    DO is=1,ions1%nsp
       ! free particle reference
       fpcor(is)=SQRT(3.0_real_8*factem/(4.0_real_8*pmar(is)*tempp1))
       avcor(is)=0.0_real_8
       DO ia=1,ions0%na(is)
          var=0.0_real_8
          rnorm=0
          DO ip=1,nphalf
             itau=ip+nphalf
             rnorm=rnorm+1
             xlm_=taup(1,ia,is,ip)-taup(1,ia,is,itau)
             ylm_=taup(2,ia,is,ip)-taup(2,ia,is,itau)
             zlm_=taup(3,ia,is,ip)-taup(3,ia,is,itau)
             CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             var=var+xlm**2+ylm**2+zlm**2
          ENDDO
          ! per atom
          IF (pimd3%np_total.EQ.1) THEN
             rcor(ia,is)=0._real_8
          ELSE
             rcor(ia,is)=SQRT(var/rnorm)
          ENDIF
          ! c          RCOR(IA,IS)=SQRT(VAR/RNORM)
          ! accumulate first and second moment
          rcora(ia,is,1)=rcora(ia,is,1)+rcor(ia,is)
          rcora(ia,is,2)=rcora(ia,is,2)+rcor(ia,is)**2
          ! per species
          avcor(is)=avcor(is)+rcor(ia,is)
       ENDDO
       avcor(is)=avcor(is)/REAL(ions0%na(is),kind=real_8)
       ! accumulate first and second moment
       avcora(is,1)=avcora(is,1)+avcor(is)
       avcora(is,2)=avcora(is,2)+avcor(is)**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getgyration
  ! ==================================================================

END MODULE getgyr_utils
