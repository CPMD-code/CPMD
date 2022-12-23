MODULE pi_init_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: scmass,&
                                             factem
  USE error_handling,                  ONLY: stopgm
  USE glemod,                          ONLY: glepar,&
                                             icm2au,&
                                             tglepc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_multiple_comm_init,           ONLY: multiple_comm_init
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: &
       fcorr, flnm, grandparent, ifcorr, ircorr, maxnp, pcg_pos, pimd1, &
       pimd2, pimd3, pma0s, pmar, pmar0, pmars, rcor, rcora, rcorr, rgyr, &
       rgyra, rsus, rsusa, tnm, tnmi, pi_omega
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntr,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
  USE wann,                            ONLY: wan05,&
                                             wannl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_init

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_init
    ! ==--------------------------------------------------------------==
    ! ==    INITIALIZE PATH INTEGRAL RUNS                             ==
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_init'

    INTEGER                                  :: ierr, iopt, ip, is, isub, jp, &
                                                kp, mcnstr, npx, nspac1, &
                                                nspac2
    REAL(real_8)                             :: aux(3*maxnp), rnp
    REAL(real_8), ALLOCATABLE                :: ap(:), aprim(:,:), &
                                                tnmtemp(:,:)

    CALL tiset(procedureN,isub)
    CALL multiple_comm_init("TROTTER")

    ! ..initalize masses
    ! DM1
    ! ..note:  the mass variables PMA(NSP), PMATOT, PMAT0 were already
    ! ..       initialized in 'setsys.F', but without the WMASS scaling
    ! ..       factor, this is now also done: 
    ! DM2
    ! 
    ! KEY TO MASS VARIABLES:
    ! PMAR0   real mass for primitive prop in chemical units
    ! PMAR    real mass for primitive prop in au
    ! PMARS   real mass for mode      prop in au
    ! PMA0    fict mass for primitive prop in chemical units
    ! DM1
    ! PMAT0   total "  "
    ! DM2
    ! PMA     fict mass for primitive prop in au
    ! DM1
    ! PMATOT  total "  "
    ! DM2
    ! PMA0S   fict mass for mode      prop in au
    ! 
    ! DM1
    rmass%pmat0=0._real_8
    rmass%pmatot=0._real_8
    ! DM2
    DO is=1,ions1%nsp
       ! DM1 delete following line
       ! STORE REAL IONIC MASS (IN CHEMICAL UNITS) IN PMAR0
       ! DM2
       ! STORE REAL IONIC MASS (IN CHEMICAL UNITS) IN PMAR0
       pmar0(is)=rmass%pma0(is)
       ! STORE REAL IONIC MASS (IN ATOMIC UNITS) IN PMAR
       pmar(is)=pmar0(is)*scmass
       ! STORE FICTITIOUS IONIC MASS (IN CHEMICAL UNITS) IN PMA0
       rmass%pma0(is)=pmar0(is)*pimd2%wmass
       ! DM1
       ! STORE TOTAL FICTITIOUS MASS (IN CHEMICAL UNITS) IN PMAT0
       rmass%pmat0=rmass%pmat0+rmass%pma0(is)*ions0%na(is)
       ! STORE FICTITIOUS IONIC MASS (IN ATOMIC UNITS) IN PMA
       rmass%pma(is)=rmass%pma0(is)*scmass
       ! STORE TOTAL FICTITIOUS MASS (IN ATOMIC UNITS) IN PMATOT
       rmass%pmatot=rmass%pmatot+rmass%pma(is)*ions0%na(is)
       ! DM2       
    ENDDO
    ! NM
    ! NM:  IF NORMAL MODES, SET UP TRANSFORMATION MATRICES AND
    ! NM   CALCULATE NORMAL MODE FREQUENCIES AND DETERMINE MASSES
    ! NM   FLNM(NPX): contains normal mode frequencies
    ! NM   TNM(NPX,NPX): contains normal mode transformation matrix
    ! NM                 PRIMITIVE->NM
    ! NM   TNMI(NPX,NPX): NM->PRIMITIVE transformation
    ! NM1
    IF (pimd1%tpinm.OR.pimd1%tread_cent) THEN
       IF (pimd3%np_total.EQ.1)CALL stopgm('PI_INIT',&
            'NO MODE TRANSFORMATION WITH NP=1',& 
            __LINE__,__FILE__)
       npx=pimd3%np_total
       rnp=REAL(npx,kind=real_8)
       ALLOCATE(tnm(pimd3%np_total,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tnmi(pimd3%np_total,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(aprim(pimd3%np_total,npx**2/pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ap(((npx+1)*npx)/2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tnmtemp(pimd3%np_total,npx**2/pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(aprim)!,npx*npx)
       DO ip=2,npx-1
          aprim(ip,ip-1) = -1._real_8
          aprim(ip,ip) = 2._real_8
          aprim(ip,ip+1) = -1._real_8
       ENDDO
       aprim(1,1) = 2._real_8
       aprim(1,2) = -1._real_8
       aprim(npx,npx) = 2._real_8
       aprim(npx,npx-1) = -1._real_8
       aprim(1,npx) = -1._real_8
       aprim(npx,1) = -1._real_8
       ! 
       kp=0
       DO ip=1,npx
          DO jp=ip,npx
             kp = kp + 1
             ap(kp) = aprim(ip,jp)
          ENDDO
       ENDDO
       ! 
       iopt=1
       CALL dspevy(iopt,ap,flnm,tnmtemp,npx,npx,aux,3*npx)
       ! 
       DO ip=1,npx
          flnm(ip) = flnm(ip)*rnp
       ENDDO
       ! 
       DO ip=1,npx
          DO jp=1,npx
             tnm(ip,jp) = -tnmtemp(jp,ip)/SQRT(rnp)
             tnmi(ip,jp) = -tnmtemp(ip,jp)*SQRT(rnp)
          ENDDO
       ENDDO
       ! SOMETIMES THE DIAGONALIZER "DSPEVY" GIVES THE WRONG
       ! PHASE OF THE DIAGONAL MATRIX, FIX THE SIGN IN THIS CASE:
       IF (tnmi(1,1) .LE. 0.0_real_8) THEN
          DO ip=1,npx
             DO jp=1,npx
                tnm(ip,jp) = -tnm(ip,jp)
                tnmi(ip,jp) = -tnmi(ip,jp)
             ENDDO
          ENDDO
       ENDIF
       ! MASS SETTINGS FOR NORMAL MODE REPRESENTATION
       ! NORMAL MODE TRAFO FOR REAL PHYSICAL MASSES (FOR POTENTIAL ENERGY)
       DO is=1,ions1%nsp
          pmars(is,1)=pmar(is)
          ! DM          PMARS(IS,1)=PMAR(IS)*FACSTAGE
          DO ip=2,npx
             pmars(is,ip)=pmar(is)*flnm(ip)
          ENDDO
          ! DM1
          ! delete following 3 lines
          ! FICTITIOUS NORMAL MODE MASSES (FOR FICTICIOUS KINETIC ENERGY)
          ! ARE WMASS-TIMES HEAVIER THAN PHYSICAL MASSES AND ARE
          ! SCALED BY THE NORMAL MODE FREQUENCIES.
          ! MASS DISPARITY FOR CENTROID DYNAMICS:
          ! THE NON-CENTROID BEADS ARE FACSTAGE TIMES LIGHTER THAN
          ! CENTROID BEADS
          ! NOTE: THIS CREATES BOLTZMANN DISTRIBUTION
          ! DM2
          ! DM1
          pma0s(is,1)=pmars(is,1)*pimd2%wmass
          DO ip=2,npx
             ! DM           DO IP=1,NPX
             ! DM            PMA0S(IS,IP)=PMARS(IS,IP)*WMASS
             pma0s(is,ip)=pmars(is,ip)*pimd2%wmass/pimd2%facstage
             ! DM2
          ENDDO
       ENDDO
       DEALLOCATE(aprim,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ap,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tnmtemp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (pimd1%tstage) THEN
       ! MASS SETTINGS FOR STAGING REPRESENTATION
       ! STAGING TRAFO FOR REAL PHYSICAL MASSES (FOR POTENTIAL ENERGY)
       npx=pimd3%np_total
       rnp=REAL(npx,kind=real_8)
       DO is=1,ions1%nsp
          ! DM1
          ! DM          PMARS(IS,1)=PMAR(IS)*FACSTAGE
          pmars(is,1)=pmar(is)
          ! DM2
          DO ip=2,npx
             pmars(is,ip)=pmar(is)*ip/(ip-1)
          ENDDO
          ! DM2
          ! delete following 2 lines
          ! FICTITIOUS STAGING MASSES (FOR FICTICIOUS KINETIC ENERGY)
          ! ARE WMASS-TIMES HEAVIER THAN PHYSICAL MASSES
          ! MASS DISPARITY FOR STAGING CENTROID DYNAMICS:
          ! THE NON-CENTROID BEADS ARE FACSTAGE TIMES LIGHTER THAN
          ! CENTROID BEADS
          ! NOTE: THIS CREATES WIGNER DISTRIBUTION
          ! DM2
          ! DM1
          pma0s(is,1)=pmars(is,1)*pimd2%wmass
          DO ip=2,npx
             ! DM          DO IP=1,NPX
             ! DM            PMA0S(IS,IP)=PMARS(IS,IP)*WMASS
             pma0s(is,ip)=pmars(is,ip)*pimd2%wmass/pimd2%facstage
          ENDDO
          ! DM2
       ENDDO
    ENDIF
    ! MASS SETTINGS FOR RING-POLYMER DYNAMICS
    IF (pimd1%tringp) THEN
       npx=pimd3%np_total
       rnp=REAL(npx,kind=real_8)

       rmass%pmat0=0._real_8
       rmass%pmatot=0._real_8
       DO is=1,ions1%nsp
          ! STORE REAL IONIC MASS (IN CHEMICAL UNITS) IN PMAR0
          pmar0(is)=rmass%pma0(is)
          ! STORE REAL IONIC MASS (IN ATOMIC UNITS) IN PMAR
          pmar(is)=pmar0(is)*scmass
          ! STORE FICTITIOUS IONIC MASS (IN CHEMICAL UNITS) IN PMA0
          rmass%pma0(is)=pmar0(is)/rnp
          ! STORE TOTAL FICTITIOUS MASS (IN CHEMICAL UNITS) IN PMAT0
          rmass%pmat0=rmass%pmat0+rmass%pma0(is)*ions0%na(is)
          ! STORE FICTITIOUS IONIC MASS (IN ATOMIC UNITS) IN PMA
          rmass%pma(is)=rmass%pma0(is)*scmass
          ! STORE TOTAL FICTITIOUS MASS (IN ATOMIC UNITS) IN PMATOT
          rmass%pmatot=rmass%pmatot+rmass%pma(is)*ions0%na(is)
       ENDDO
       IF (pimd1%tpinm) THEN
          ! MASS SETTINGS FOR NORMAL MODE REPRESENTATION
          ! NORMAL MODE TRAFO FOR REAL PHYSICAL MASSES (FOR POTENTIAL ENERGY)
          DO is=1,ions1%nsp
             pmars(is,1)=pmar(is)
             DO ip=2,npx
                pmars(is,ip)=pmar(is)*flnm(ip)
             ENDDO
             ! FICTITIOUS NORMAL MODE MASSES (FOR FICTICIOUS KINETIC ENERGY)
             ! ARE PHYSICAL MASSES FOR ALL NORMAL MODES
             DO ip=1,npx
                pma0s(is,ip)=pmar(is)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! DM1
    ! ..print detailed information about final mass settings in PI case
    IF (grandparent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' ****************** FINAL MASS ',&
            'SETTINGS IN AMU ******************'
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' PRIMITIVE MASSES '
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  TYPE       REAL MASS        FICTITIOUS MASS'
       DO is=1,ions1%nsp
          IF (paral%io_parent)&
               WRITE(6,'(2X,A2,6X,F10.4,7X,F10.4)')&
               elem%el(ions0%iatyp(is)),pmar0(is),rmass%pma0(is)
       ENDDO
       IF (pimd1%tstage.OR.pimd1%tpinm) THEN
          DO is=1,ions1%nsp
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' MODE TRANSFORMED MASSES '
             IF (paral%io_parent)&
                  WRITE(6,'(A)')&
                  '  TYPE       REAL MASS        FICTITIOUS MASS FOR IP'
             DO ip=1,pimd3%np_total
                IF (paral%io_parent)&
                     WRITE(6,'(2X,A2,6X,F10.4,7X,F10.4,11X,I4)')&
                     elem%el(ions0%iatyp(is)),&
                     pmars(is,ip)/scmass,pma0s(is,ip)/scmass,ip
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.4)') ' FICTITIOUS ELECTRON MASS: ',&
            cntr%emass/scmass
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' *********************************',&
            '*******************************'
    ENDIF
    ! ..Constraints
    IF (pcg_pos.NE.1) mcnstr=0
    ! ..Allocate space for radius of gyration
    ALLOCATE(rgyr(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rgyra(maxsys%nax,maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rsus(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rsusa(maxsys%nax,maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rcor(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rcora(maxsys%nax,maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..Allocate space for correlation function
    nspac1=maxsys%nax*maxsys%nsx*pimd3%np_total/2 + 1
    nspac2=maxsys%nax*maxsys%nsx*pimd3%np_total
    ALLOCATE(ircorr(maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ifcorr(maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rcorr(maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fcorr(maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..Distribute information to all parents
    ! ..FPATH
    !if (paral%parent) then
    !   call mp_bcast(fo_info%fpath,supersource,parentgroup)
    !endif
    ! ==--------------------------------------------------------------==

    IF (wannl%twann) THEN
       IF (wan05%loc_npgrp>parai%cp_nproc) wan05%loc_npgrp=parai%cp_nproc
    ENDIF

    ! Friction coefficient for GLE thermostat
    IF (glepar%gle_mode>0) THEN
       npx=pimd3%np_total
       rnp=REAL(npx,kind=real_8)
       DO ip=1,npx
          pi_omega(ip)=2._real_8*SQRT(rnp)*cntr%tempw/factem*pimd2%gle_lambda/icm2au
          IF ((pimd1%tpinm.OR.pimd1%tstage).AND.ip.EQ.1) THEN
             IF ((pimd1%tcentro.OR.pimd1%tringp).AND..NOT.tglepc) THEN
                pi_omega(ip) = 0._real_8
             ELSE
                pi_omega(ip) = 2._real_8*glepar%gle_omega*pimd2%gle_lambda
             ENDIF
          ENDIF
          IF (pimd1%tringp.AND.pimd1%tpinm.AND.ip.GT.1) THEN
             pi_omega(ip) = pi_omega(ip)*SQRT(flnm(ip))
          ENDIF
          IF (pimd1%tcentro.AND.(pimd1%tpinm.OR.pimd1%tstage).AND.ip.GT.1) THEN
             pi_omega(ip) = pi_omega(ip)*SQRT(pimd2%facstage)
          ENDIF
       ENDDO
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') '    IP     POT_OMEGA(CM^-1)     GLE_OMEGA(CM^-1)'
          DO ip=1,npx
             WRITE(6,'(1X,I5,5X,F16.4,5X,F16.4)') &
             ip,pi_omega(ip)/pimd2%gle_lambda,pi_omega(ip)
          ENDDO
       ENDIF
    ENDIF
 
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pi_init
  ! ==================================================================

END MODULE pi_init_utils
