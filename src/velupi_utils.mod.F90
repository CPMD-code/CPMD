MODULE velupi_utils
  USE cnst,                            ONLY: au_fs
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com,&
                                             veps
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: pma0s
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE puttau_utils,                    ONLY: taucl
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntr,&
                                             iatpt,&
                                             maxsys
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE utils,                           ONLY: invmat
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: velupi
  PUBLIC :: velupidamp
  PUBLIC :: velupishock
  PUBLIC :: velupc1
  PUBLIC :: velupc2
  PUBLIC :: s_velupi
  PUBLIC :: velupif
  PUBLIC :: velupc
  PUBLIC :: s_to_c
  PUBLIC :: c_to_s

CONTAINS

  ! ==================================================================
  SUBROUTINE velupi(velp,fion,nnlst)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), fion(:,:,:)
    INTEGER                                  :: nnlst

    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: fact

!ocl NOALIAS

    !$omp parallel do private(I,IS,IA,FACT) schedule(static)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       fact=REAL(nnlst,kind=real_8)*dtb2mi(is)
       velp(1,ia,is)=velp(1,ia,is)+fact*fion(1,ia,is)
       velp(2,ia,is)=velp(2,ia,is)+fact*fion(2,ia,is)
       velp(3,ia,is)=velp(3,ia,is)+fact*fion(3,ia,is)
    ENDDO
    CALL taucl(velp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupi
  ! ==================================================================
  SUBROUTINE velupidamp(velp,nnlst)
    REAL(real_8)                             :: velp(:,:,:)
    INTEGER                                  :: nnlst

    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: alpha, gamma, ikin, kin, scale

    kin = 0._real_8
    !ocl NOALIAS
    !$omp parallel do private(I,IS,IA,SCALE) schedule(static) &
    !$omp  reduction(+:KIN)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       scale=rmass%pma(is)
       kin = kin + scale*velp(1,ia,is)*velp(1,ia,is)
       kin = kin + scale*velp(2,ia,is)*velp(2,ia,is)
       kin = kin + scale*velp(3,ia,is)*velp(3,ia,is)
    ENDDO
    ! CJM DBG HARDWIRED UNTIL PUT IN INPUT
    gamma = 0.05_real_8 * au_fs
    ! CJM DBG HARDWIRED UNTIL PUT IN INPUT
    ikin = 1._real_8/kin
    alpha = 2._real_8*cntr%cmass*veps*veps*gamma*ikin
    scale = SQRT ( 1._real_8 + alpha * 0.5_real_8 * dt_ions )
    !ocl NOALIAS
    !$omp parallel do private(I,IS,IA) schedule(static)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       velp(1,ia,is)=velp(1,ia,is)*scale
       velp(2,ia,is)=velp(2,ia,is)*scale
       velp(3,ia,is)=velp(3,ia,is)*scale
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupidamp
  ! ==================================================================
  SUBROUTINE velupishock(velp,fion,nnlst)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), fion(:,:,:)
    INTEGER                                  :: nnlst

    REAL(real_8), PARAMETER :: e2 = 1._real_8/6._real_8, e4 = e2/20._real_8, &
      e6 = e4/42._real_8, e8 = e6/72._real_8

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: aa, aax, arg_v, arg_vx, at3, &
                                                fact, poly_v, poly_vx, &
                                                scale_v, scale_vx

! ==--------------------------------------------------------------==

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    aax = veps * (1.0_real_8+1.0_real_8/at3)
    arg_vx = ( 0.25_real_8 * dt_ions * aax ) * ( 0.25_real_8 * dt_ions * aax )
    poly_vx = 1._real_8+e2*arg_vx+e4*arg_vx*arg_vx+e6*arg_vx**3+&
         e8*arg_vx**4

    aa = veps * 1.0_real_8/at3
    arg_v = ( 0.25_real_8 * dt_ions * aa ) * ( 0.25_real_8 * dt_ions * aa )
    poly_v = 1._real_8+e2*arg_v+e4*arg_v*arg_v+e6*arg_v**3+&
         e8*arg_v**4
    ! AT3 = 3._real_8*real(NAT,kind=real_8)
    ! aax = VEPS * (1.0_real_8+1.0_real_8/AT3)
    ! arg_vx = ( 0.25_real_8 * DT_IONS * aax ) * ( 0.25_real_8 * DT_IONS * aax )
    ! poly_vx = 1._real_8+e2*arg_vx+e4*arg_vx*arg_vx+e6*arg_vx**3+
    ! *         e8*arg_vx**4
    ! AT3 = 3._real_8*real(NAT,kind=real_8)
    ! aax = HTVEL ( 1, 1 ) * (1.0_real_8+1.0_real_8/AT3)
    ! arg_vx = ( 0.25_real_8 * DT_IONS * aax ) * ( 0.25_real_8 * DT_IONS * aax )
    ! poly_vx = 1._real_8+e2*arg_vx+e4*arg_vx*arg_vx+e6*arg_vx**3+
    ! *          e8*arg_vx**4
    ! 
    ! aa = HTVEL ( 1, 1 ) * 1.0_real_8/AT3
    ! arg_v = ( 0.25_real_8 * DT_IONS * aa ) * ( 0.25_real_8 * DT_IONS * aa )
    ! poly_v = 1._real_8+e2*arg_v+e4*arg_v*arg_v+e6*arg_v**3+
    ! *          e8*arg_v**4

    scale_vx = EXP( -0.25_real_8 * dt_ions * aax )
    scale_v = EXP( -0.25_real_8 * dt_ions * aa )

    !ocl NOALIAS
    !$omp parallel do private(IS,IA,FACT)
    DO is=1,ions1%nsp
       fact=REAL(nnlst,kind=real_8)*dtb2mi(is)
       DO ia=1,ions0%na(is)
          velp(1,ia,is)=scale_vx*scale_vx*velp(1,ia,is)+&
               scale_vx*poly_vx*fact*fion(1,ia,is)
          velp(2,ia,is)=scale_v*scale_v*velp(2,ia,is)+&
               scale_v*poly_v*fact*fion(2,ia,is)
          velp(3,ia,is)=scale_v*scale_v*velp(3,ia,is)+&
               scale_v*poly_v*fact*fion(3,ia,is)
       ENDDO
    ENDDO
    ! DTH = 0.25_real_8*DT_IONS
    ! for y and z (doesn't couple to shocks)
    ! S1 = EXP(-0.5_real_8*DT_IONS*HTVEL ( 1, 1 )*1.0_real_8/AT3)
    ! ARG = DTH*HTVEL ( 1, 1 )*1.0_real_8/AT3
    ! SE = EXP(-ARG)
    ! S2 = SE*SINHXBX(ARG)
    ! !for x-direction (couples to shock)
    ! S1X = EXP(-0.5_real_8*DT_IONS*HTVEL ( 1, 1 )*(1.0_real_8+1.0_real_8/AT3))
    ! ARGX = DTH*HTVEL ( 1, 1 )*(1.0_real_8+1.0_real_8/AT3)
    ! SEX = EXP(-ARGX)
    ! S2X = SEX*SINHXBX(ARG)
    ! DO IS=1,NSP
    ! FACT=real(NNLST,kind=real_8)*DTB2MI(IS)
    ! DO IA=1,NA(IS)
    ! VELP(1,IA,IS)=S1X*VELP(1,IA,IS)+S2X*FACT*FION(1,IA,IS)
    ! VELP(2,IA,IS)=S1*VELP(2,IA,IS)+S2*FACT*FION(2,IA,IS)
    ! VELP(3,IA,IS)=S1*VELP(3,IA,IS)+S2*FACT*FION(3,IA,IS)
    ! ENDDO
    ! ENDDO
    CALL taucl(velp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupishock
  ! ==================================================================
  SUBROUTINE velupc1(svelp,sfion,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: svelp(:,:,:), sfion(:,:,:), &
                                                velp(:,:,:)

    INTEGER                                  :: ia, info, is, k, l
    REAL(real_8)                             :: aux(1000), dtb, fact, fomega, &
                                                fsd(3), gdot(3,3), gmat(3,3), &
                                                hforce(3,3), hscr(3,3)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL dgemm('N','T',3,3,3,1._real_8,metr_com%htvel(1,1),3,metr_com%ht(1,1)   ,3,0._real_8,&
         gdot(1,1),3)
    CALL dgemm('N','T',3,3,3,1._real_8,metr_com%ht(1,1)   ,3,metr_com%htvel(1,1),3,1._real_8,&
         gdot(1,1),3)
    CALL dgemm('N','T',3,3,3,1._real_8,metr_com%ht(1,1)   ,3,metr_com%ht(1,1)   ,3,0._real_8,&
         gmat(1,1),3)
    CALL invmat(3,gmat,aux,info)
    CALL dgemm('N','N',3,3,3,1._real_8,gmat(1,1) ,3,gdot(1,1) ,3,0._real_8,&
         hscr(1,1),3)
    CALL dcopy(9,hscr(1,1),1,gmat(1,1),1)
    CALL dcopy(9,metr_com%htfor(1,1),1,hforce(1,1),1)
    CALL zeroing(hscr)!,9)
    !ocl NOALIAS
    DO k=1,ions1%nat
       ia=iatpt(1,k)
       is=iatpt(2,k)
       fact=rmass%pma(is)
       hscr(1,1)=hscr(1,1)+fact*velp(1,ia,is)*velp(1,ia,is)
       hscr(2,1)=hscr(2,1)+fact*velp(2,ia,is)*velp(1,ia,is)
       hscr(3,1)=hscr(3,1)+fact*velp(3,ia,is)*velp(1,ia,is)
       hscr(1,2)=hscr(1,2)+fact*velp(1,ia,is)*velp(2,ia,is)
       hscr(2,2)=hscr(2,2)+fact*velp(2,ia,is)*velp(2,ia,is)
       hscr(3,2)=hscr(3,2)+fact*velp(3,ia,is)*velp(2,ia,is)
       hscr(1,3)=hscr(1,3)+fact*velp(1,ia,is)*velp(3,ia,is)
       hscr(2,3)=hscr(2,3)+fact*velp(2,ia,is)*velp(3,ia,is)
       hscr(3,3)=hscr(3,3)+fact*velp(3,ia,is)*velp(3,ia,is)
    ENDDO
    CALL dgemm('N','N',3,3,3,1._real_8,hscr(1,1),3,metr_com%htm1(1,1),3,1._real_8,&
         hforce(1,1),3)

    fact=dt_ions/(2._real_8*cntr%cmass)
    IF (prcpl%tzflex) THEN
       ! Z-scaling of cell only
       fomega=hforce(1,3)*prcp_com%hunit(1,3)+hforce(2,3)*prcp_com%hunit(2,3)&
            +hforce(3,3)*prcp_com%hunit(3,3)
       DO l=1,3
          metr_com%htvel(l,3)=metr_com%htvel(l,3)+fact*fomega*prcp_com%hunit(l,3)
       ENDDO
    ELSE IF (prcpl%tisot) THEN
       ! Isotropic cell
       fomega=ddot(9,hforce,1,prcp_com%hunit,1)/3.0_real_8
       DO k=1,3
          DO l=1,3
             metr_com%htvel(l,k)=metr_com%htvel(l,k)+fact*fomega*prcp_com%hunit(l,k)
          ENDDO
       ENDDO
    ELSE
       ! General
       DO k=1,3
          DO l=1,3
             metr_com%htvel(l,k)=metr_com%htvel(l,k)+fact*hforce(l,k)
          ENDDO
       ENDDO
    ENDIF
    dtb=dt_ions*0.5_real_8
    DO is=1,ions1%nsp
       fact=dtb2mi(is)
       DO ia=1,ions0%na(is)
          CALL dgemv('N',3,3,1._real_8,gmat(1,1),3,svelp(1,ia,is),1,0._real_8,&
               fsd(1),1)
          svelp(1,ia,is)=svelp(1,ia,is)+fact*sfion(1,ia,is)-fsd(1)*dtb
          svelp(2,ia,is)=svelp(2,ia,is)+fact*sfion(2,ia,is)-fsd(2)*dtb
          svelp(3,ia,is)=svelp(3,ia,is)+fact*sfion(3,ia,is)-fsd(3)*dtb
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupc1
  ! ==================================================================
  SUBROUTINE velupc2(svelp,sfion,velp,tscr)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: svelp(:,:,:), sfion(:,:,:), &
                                                velp(:,:,:), tscr(:,:,:)

    INTEGER, PARAMETER                       :: maxits = 100 
    REAL(real_8), PARAMETER                  :: epss = 1.e-10_real_8 

    INTEGER                                  :: ia, iloop, info, is, k, l
    REAL(real_8) :: aux(1000), dtb, ei, fact, fomega, fsd(3), gdot(3,3), &
      gmat(3,3), hc(3,3), hcr(3,3), hforce(3,3), hscr(3,3)
    REAL(real_8), EXTERNAL                   :: dasum, ddot

    CALL dcopy(9,        metr_com%htvel(1,1)  ,1,hc(1,1)    ,1)
    CALL dcopy(9,        metr_com%htvel(1,1)  ,1,hcr(1,1)   ,1)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,svelp(1,1,1),1,tscr(1,1,1),1)
    DO iloop=1,maxits
       CALL s_to_c(svelp,velp)
       CALL dgemm('N','T',3,3,3,1._real_8,metr_com%htvel(1,1),3,metr_com%ht(1,1)   ,3,0._real_8,&
            gdot(1,1),3)
       CALL dgemm('N','T',3,3,3,1._real_8,metr_com%ht(1,1)   ,3,metr_com%htvel(1,1),3,1._real_8,&
            gdot(1,1),3)
       CALL dgemm('N','T',3,3,3,1._real_8,metr_com%ht(1,1)   ,3,metr_com%ht(1,1)   ,3,0._real_8,&
            gmat(1,1),3)
       CALL invmat(3,gmat,aux,info)
       CALL dgemm('N','N',3,3,3,1._real_8,gmat(1,1) ,3,gdot(1,1) ,3,0._real_8,&
            hscr(1,1),3)
       CALL dcopy(9,hscr(1,1) ,1,gmat(1,1)  ,1)
       CALL dcopy(9,metr_com%htfor(1,1),1,hforce(1,1),1)
       CALL dcopy(9,hcr(1,1)  ,1,metr_com%htvel(1,1) ,1)
       CALL zeroing(hscr)!,9)
       !ocl NOALIAS
       DO k=1,ions1%nat
          ia=iatpt(1,k)
          is=iatpt(2,k)
          fact=rmass%pma(is)
          hscr(1,1)=hscr(1,1)+fact*velp(1,ia,is)*velp(1,ia,is)
          hscr(2,1)=hscr(2,1)+fact*velp(2,ia,is)*velp(1,ia,is)
          hscr(3,1)=hscr(3,1)+fact*velp(3,ia,is)*velp(1,ia,is)
          hscr(1,2)=hscr(1,2)+fact*velp(1,ia,is)*velp(2,ia,is)
          hscr(2,2)=hscr(2,2)+fact*velp(2,ia,is)*velp(2,ia,is)
          hscr(3,2)=hscr(3,2)+fact*velp(3,ia,is)*velp(2,ia,is)
          hscr(1,3)=hscr(1,3)+fact*velp(1,ia,is)*velp(3,ia,is)
          hscr(2,3)=hscr(2,3)+fact*velp(2,ia,is)*velp(3,ia,is)
          hscr(3,3)=hscr(3,3)+fact*velp(3,ia,is)*velp(3,ia,is)
       ENDDO
       CALL dgemm('N','N',3,3,3,1._real_8,hscr(1,1),3,metr_com%htm1(1,1),3,1._real_8,&
            hforce(1,1),3)

       fact=dt_ions/(2._real_8*cntr%cmass)
       IF (prcpl%tzflex) THEN
          ! Z-scaling of cell only
          fomega=hforce(1,3)*prcp_com%hunit(1,3)+hforce(2,3)*prcp_com%hunit(2,3)&
               +hforce(3,3)*prcp_com%hunit(3,3)
          DO l=1,3
             metr_com%htvel(l,3)=metr_com%htvel(l,3)+fact*fomega*prcp_com%hunit(l,3)
          ENDDO
       ELSE IF (prcpl%tisot) THEN
          ! Isotropic cell
          fomega=ddot(9,hforce,1,prcp_com%hunit,1)/3.0_real_8
          DO k=1,3
             DO l=1,3
                metr_com%htvel(l,k)=metr_com%htvel(l,k)+fact*fomega*prcp_com%hunit(l,k)
             ENDDO
          ENDDO
       ELSE
          DO k=1,3
             DO l=1,3
                metr_com%htvel(l,k)=metr_com%htvel(l,k)+fact*hforce(l,k)
             ENDDO
          ENDDO
       ENDIF
       dtb=0.5_real_8*dt_ions
       DO is=1,ions1%nsp
          fact=dtb2mi(is)
          DO ia=1,ions0%na(is)
             CALL dgemv('N',3,3,1._real_8,gmat(1,1),3,svelp(1,ia,is),1,0._real_8,&
                  fsd(1),1)
             svelp(1,ia,is)=tscr(1,ia,is)+fact*sfion(1,ia,is)-fsd(1)*dtb
             svelp(2,ia,is)=tscr(2,ia,is)+fact*sfion(2,ia,is)-fsd(2)*dtb
             svelp(3,ia,is)=tscr(3,ia,is)+fact*sfion(3,ia,is)-fsd(3)*dtb
          ENDDO
       ENDDO
       CALL daxpy(9,-1._real_8,metr_com%htvel(1,1),1,hc(1,1),1)
       ei=dasum(9,hc(1,1),1)
       IF (ei.LT.epss) GOTO 100
       CALL dcopy(9,metr_com%htvel(1,1),1,hc(1,1),1)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' VELUPC2| ITERATION DID NOT CONVERGE '
    IF (paral%io_parent)&
         WRITE(6,*) ' VELUPC2| NUMBER OF STEPS            ',maxits
    IF (paral%io_parent)&
         WRITE(6,*) ' VELUPC2| LAST CHANGE IN HTVEL WAS   ',ei
    IF (paral%io_parent)&
         WRITE(6,*) ' VELUPC2| CONVERGENCE CRITERIA IS    ',epss
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupc2
  ! ==================================================================
  SUBROUTINE s_velupi(velp,fion,nnlst,ip)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), fion(:,:,:)
    INTEGER                                  :: nnlst, ip

    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: fact

! Variables
! ==--------------------------------------------------------------==
!ocl NOALIAS

    !$omp parallel do private(I,IS,IA,FACT) schedule(static)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       fact=REAL(nnlst,kind=real_8)*dtb2mi(is)*rmass%pma(is)/pma0s(is,ip)
       velp(1,ia,is)=velp(1,ia,is)+fact*fion(1,ia,is)
       velp(2,ia,is)=velp(2,ia,is)+fact*fion(2,ia,is)
       velp(3,ia,is)=velp(3,ia,is)+fact*fion(3,ia,is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s_velupi
  ! ==================================================================
  SUBROUTINE velupif(velp,fion,nnlst)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), fion(:,:,:)
    INTEGER                                  :: nnlst

    CALL taucl(velp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupif
  ! ==================================================================
  SUBROUTINE velupc(nnlst)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nnlst

    INTEGER                                  :: i, j
    REAL(real_8)                             :: fact, fomega
    REAL(real_8), EXTERNAL                   :: ddot

    fact=REAL(nnlst,kind=real_8)*dt_ions/(2._real_8*cntr%cmass)
    IF (prcpl%tzflex) THEN
       ! Z-scaling of cell only
       fomega=metr_com%htfor(1,3)*prcp_com%hunit(1,3)+metr_com%htfor(2,3)*prcp_com%hunit(2,3)&
            +metr_com%htfor(3,3)*prcp_com%hunit(3,3)
       DO i=1,3
          metr_com%htvel(i,3)=metr_com%htvel(i,3)+fact*fomega*prcp_com%hunit(i,3)
       ENDDO
       veps = metr_com%htvel(3,3)
    ELSE IF (prcpl%tisot) THEN
       ! Isotropic cell
       fomega=ddot(9,metr_com%htfor,1,prcp_com%hunit,1)/3.0_real_8
       DO i=1,3
          DO j=1,3
             metr_com%htvel(j,i)=metr_com%htvel(j,i)+fact*fomega*prcp_com%hunit(j,i)
          ENDDO
       ENDDO
       veps = (metr_com%htvel(1,1)+metr_com%htvel(2,2)+metr_com%htvel(3,3))/3._real_8
    ELSE
       ! General
       DO i=1,3
          DO j=1,3
             metr_com%htvel(j,i)=metr_com%htvel(j,i)+fact*metr_com%htfor(j,i)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupc
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE s_to_c(scoor,coor)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: scoor(:,:,:), coor(:,:,:)

    INTEGER                                  :: ia, is

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,scoor(1,ia,is),1,&
               0.0_real_8,coor(1,ia,is),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s_to_c
  ! ==================================================================
  SUBROUTINE c_to_s(coor,scoor)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: coor(:,:,:), scoor(:,:,:)

    INTEGER                                  :: ia, is

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,coor(1,ia,is),1,&
               0.0_real_8,scoor(1,ia,is),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE c_to_s

END MODULE velupi_utils
