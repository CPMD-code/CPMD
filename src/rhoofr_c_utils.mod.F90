MODULE rhoofr_c_utils
  USE augchg_utils,                    ONLY: augchg
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE density_utils,                   ONLY: build_density_sum
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE moverho_utils,                   ONLY: moverho
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE prcp,                            ONLY: prcp_com
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace,&
                                             type_cast
  USE rhov_utils,                      ONLY: rhov
  USE ropt,                            ONLY: ropt_mod
  USE rswfmod,                         ONLY: maxstates,&
                                             rsactive,&
                                             rswf
  USE sfac,                            ONLY: fnl
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE state_utils,                     ONLY: set_psi_1_state_g_kpts
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: symrho
  USE system,                          ONLY: &
       cntl, dual00, fpar, group, kpbeg, ncpw, nkpbl, nkpt, parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhoofr_c

CONTAINS

  ! ==================================================================
  SUBROUTINE rhoofr_c(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==  THE KINETIC ENERGY EKIN. IT IS DONE IN RECIPROCAL SPACE     ==
    ! ==  WHERE THE ASSOCIATED OPERATORS ARE DIAGONAL.                ==
    ! ==  RHOE IS OBTAINED FOURIER TRANSFORMING THE WFN TO REAL       ==
    ! ==  SPACE (PSI).                                                ==
    ! ==--------------------------------------------------------------==
    ! ==  WARNING: IF YOU USE SPECIAL K-POINTS FOR A SPECIAL STRUCTURE==
    ! ==           YOU NEED TO SYMMETRIZE CHARGE DENSITY              ==
    ! ==           FOR THAT -> SPECIFY THE POINT GROUP                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhoofr_c'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)
    REAL(real_8), PARAMETER                  :: delta = 1.e-6_real_8, &
                                                deltakin = 1.e-10_real_8

    INTEGER                                  :: i, ia, iat, ierr, ig, ikind, &
                                                ikk, ikpt, is, is1, ispin1, &
                                                isub, isub3, iwf, kbeg, kend, &
                                                kinc, nkpoint
    LOGICAL                                  :: tfcal
    REAL(real_8)                             :: argm, argp, coef3, g2m, g2p, &
                                                rsum, rsum1, rsum1abs, rsumv, &
                                                sk1, xkin, xskin
    REAL(real_8), ALLOCATABLE                :: qa(:)
    REAL(real_8), EXTERNAL                   :: dasum, ddot
    REAL(real_8), POINTER                    :: psix(:), rhoeg(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    CALL setfftn(0)

    ! ==--------------------------------------------------------------== 

    CALL TYPE_CAST(PSI, SIZE(PSI), PSIX) ! ! PSI IS C(8), PSIX IS REAL(8)
    CALL RESHAPE_INPLACE(RHOE, (/fpar%nnr1 * clsd%nlsd/), RHOEG) ! ! RHOEG(:) => RHOE(:,:)
    ! Initialize
    CALL zeroing(rhoe)!,clsd%nlsd*nnr1)
    ! Accumulate the charge and kinetic energy
    rsum=0._real_8
    xkin=0._real_8
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,'HGKP HGKM C0')
       DO ikind = 1, nkpoint
          ikk=kpbeg(ikpt)+ikind
          DO i=1,nstate
             IF (crge%f(i,ikk).NE.0._real_8) THEN
                rsum=rsum+wk(ikk)*crge%f(i,ikk)*&
                     ddot(nkpt%ngwk*2,c0(1,i,ikind),1,c0(1,i,ikind),1)
                sk1=0.0_real_8
                IF (prcp_com%akin.GT.deltakin) THEN
                   xskin=1._real_8/prcp_com%gskin
                   DO ig=1,ncpw%ngw
                      argp=(hgkp(ig,ikind)-prcp_com%gckin)*xskin
                      argm=(hgkm(ig,ikind)-prcp_com%gckin)*xskin
                      g2p=hgkp(ig,ikind)+prcp_com%gakin*(1._real_8+cp_erf(argp))
                      g2m=hgkm(ig,ikind)+prcp_com%gakin*(1._real_8+cp_erf(argm))
                      sk1=sk1+g2p*ABS(c0(ig,i,ikind))**2+&
                           g2m*ABS(c0(ig+ncpw%ngw,i,ikind))**2
                   ENDDO
                ELSE
                   DO ig=1,ncpw%ngw
                      sk1=sk1+hgkp(ig,ikind)*ABS(c0(ig,i,ikind))**2+&
                           hgkm(ig,ikind)*ABS(c0(ig+ncpw%ngw,i,ikind))**2
                   ENDDO
                ENDIF
                xkin=xkin+0.5_real_8*wk(ikk)*crge%f(i,ikk)*sk1
             ENDIF
          ENDDO

          ! Loop over the electronic states
          DO i = 1,part_1d_nbr_el_in_blk(nstate,parai%cp_inter_me,parai%cp_nogrp)
             is1 = part_1d_get_el_in_blk(i,nstate,parai%cp_inter_me,parai%cp_nogrp)
             tfcal=rsactive.OR.(crge%f(is1,ikk).NE.0._real_8)

             IF (tfcal) THEN
                CALL zeroing(psi)!,maxfft)

                CALL set_psi_1_state_g_kpts(zone,c0(:,is1,ikind),psi)

                ! Fourier transform the wave functions to real space.

                CALL invfftn(psi,.TRUE.,parai%allgrp)

                ! Store real space wavefunctions
                IF (rsactive) THEN
                   iwf=is1+(ikind-1)*nstate
                   IF (iwf.LE.maxstates) THEN
                      CALL dcopy(2*fpar%nnr1,psi(1),1,rswf(1,iwf),1)
                   ENDIF
                ENDIF
                ! Compute the charge density from the wave functions
                ! in real space 
                coef3=wk(ikk)*crge%f(is1,ikk)/parm%omega
                IF (coef3.NE.0._real_8) THEN
                   IF (cntl%tlsd) THEN
                      ispin1=1
                      IF (is1.GT.spin_mod%nsup) ispin1=2
                      CALL build_density_sum(coef3,coef3,psi,&
                           rhoe(:,ispin1),fpar%nnr1)
                   ELSE
                      CALL build_density_sum(coef3,coef3,psi,&
                           rhoe(:,1),fpar%nnr1)
                   ENDIF
                ENDIF
             ENDIF         ! Endif for TFCAL
          ENDDO             ! End loop over electronic states
       ENDDO                 ! End loop over k points (IKIND)
    ENDDO                     ! End loop over IKPT
    ! ==--------------------------------------------------------------==
    ener_com%ekin=xkin*parm%tpiba2
    ! ==--------------------------------------------------------------==

    ! 
    ! redistribute RHOE over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_grps_b',isub3)
       CALL cp_grp_redist(rhoe,fpar%nnr1,clsd%nlsd)
       CALL tihalt(procedureN//'_grps_b',isub3)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! MOVE DENSITY ACCORDING TO MOVEMENT OF ATOMS
    IF (ropt_mod%modens) CALL moverho(rhoe,psi)
    ! CONTRIBUTION OF THE VANDERBILT PP TO RHOE
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          ! ALPHA SPIN
          CALL rhov(nstate,1,spin_mod%nsup,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
          ! BETA SPIN
          CALL rhov(nstate,spin_mod%nsup+1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,2)=rhoe(i,2)+REAL(psi(i))
          ENDDO
       ELSE
          CALL rhov(nstate,1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
       ENDIF
       ! Vanderbilt Charges
       IF (paral%parent) THEN
          ALLOCATE(qa(ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL zeroing(qa)!,ions1%nat)
          CALL augchg(fnl,crge%f,qa,nstate)
          iat=0
          DO is=1,ions1%nsp
             chrg%vdbchg(is)=0.0_real_8
             DO ia=1,ions0%na(is)
                iat=iat+1
                chrg%vdbchg(is)=chrg%vdbchg(is)+qa(iat)
             ENDDO
             chrg%vdbchg(is)=chrg%vdbchg(is)/ions0%na(is)
          ENDDO
          DEALLOCATE(qa,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2)
    IF (cntl%tlsd) THEN
       !$omp parallel do private(I)
       DO i=1,fpar%nnr1
          rhoe(i,1) = rhoe(i,1) + rhoe(i,2)
       ENDDO
    ENDIF
    ! SYMMETRIZE DENSITY IF POINT GROUP SPECIFIED 
    ! (NEED FOR SPECIAL K-POINTS).
    IF (cntl%tsymrho) THEN
       IF (cntl%tlsd) THEN
          CALL symrho(rhoe(:,1),psi)
          CALL symrho(rhoe(:,2),psi)
       ELSE
          CALL symrho(rhoe(:,1),psi)
       ENDIF
    ENDIF
    ! HERE TO CHECK THE INTEGRAL OF THE CHARGE DENSITY 
    rsum1=dasum(fpar%nnr1,rhoe(1,1),1)
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    chrg%csumg=rsum
    chrg%csumr=rsum1
    IF (cntl%tlsd) THEN
       rsum1=0._real_8
       rsum1abs=0._real_8
       DO i=1,fpar%nnr1
          rsum1 = rsum1 + rhoe(i,1)-2._real_8*rhoe(i,2)
          rsum1abs = rsum1abs + ABS(rhoe(i,1)-2._real_8*rhoe(i,2))
       ENDDO
       chrg%csums=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       chrg%csumsabs=rsum1abs*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSE
       chrg%csums=0._real_8
       chrg%csumsabs=0._real_8
    ENDIF
    CALL mp_sum(chrg%csumg,parai%allgrp)
    CALL mp_sum(chrg%csumr,parai%allgrp)
    CALL mp_sum(chrg%csums,parai%allgrp)
    CALL mp_sum(chrg%csumsabs,parai%allgrp)
    IF (paral%io_parent.AND.ABS(chrg%csumr-chrg%csumg).GT.delta) THEN
       WRITE(6,'(A,T45,F20.12)') 'IN FOURIER SPACE:', chrg%csumg
       WRITE(6,'(A,T45,F20.12)') 'IN REAL SPACE:', chrg%csumr
       IF (symmi%indpg.NE.0.AND.dual00%cdual.LT.4._real_8) WRITE(6,*) 'YOUR DUAL NUMBER '&
            ,dual00%cdual,&
            ' COULD BE TOO SMALL WITH DENSITY SYMMETRISATION'
       CALL stopgm(procedureN,'TOTAL DENSITY SUMS ARE NOT EQUAL',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhoofr_c
  ! ==================================================================

END MODULE rhoofr_c_utils
