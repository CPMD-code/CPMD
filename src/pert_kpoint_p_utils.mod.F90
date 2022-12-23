MODULE pert_kpoint_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE dnlpdk_p_utils,                  ONLY: dnlps_dk,&
                                             give_scr_putwnl_kpert
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forces_p_utils,                  ONLY: give_scr_forces_p
  USE gvec,                            ONLY: gvec_com
  USE h0psi1_p_utils,                  ONLY: h0psi1_p
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpert_potential_p_utils,         ONLY: fnl_fkp_p,&
                                             kin_pert
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE ks_ener_p_utils,                 ONLY: give_scr_ks_ener_p,&
                                             ks_ener_p
  USE matrix_p_utils,                  ONLY: give_scr_matrix_p,&
                                             matrix_p
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: energy
  USE response_pmod,                   ONLY: &
       ddfnl_ddk, ddtwnl_ddk, dfnl_dk, dtwnl_dk, e2_nl_c0, fnl00, kpertpar, &
       response1, response_read, response_write, rrk, wwk
  USE restart_p_utils,                 ONLY: restart_p
  USE rhoofr_p_utils,                  ONLY: give_scr_rhoofr_p
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rnl_dk_p_utils,                  ONLY: give_scr_rnl_dk_p,&
                                             rnl_dk_p
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE up3_p_utils,                     ONLY: give_scr_up3_p,&
                                             up3_p
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pert_kpoint_p
  PUBLIC :: give_scr_pert_kpoint_p

CONTAINS

  ! ==================================================================
  SUBROUTINE pert_kpoint_p(c0,v1_nonloc,psi,rhoe,&
       eirop,eivps,z11,nstate)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: v1_nonloc(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pert_kpoint_p'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/'x','y','z'/)

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: cu1(:,:,:), force1(:,:)
    INTEGER                                  :: i, ierr, ikk, il_v1_nonloc, &
                                                ir, ir1, ir2, istate, isub, &
                                                k, k2, lcu1, lfnl
    LOGICAL                                  :: fkpfirst
    REAL(real_8)                             :: epertkp, rdummy(1), &
                                                rk_in(3,nkpt%nkpts), &
                                                sder(3,3), sumkk
    REAL(real_8), ALLOCATABLE                :: sumkikj(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('    pert_kpoint_p',isub)

    IF ((clsd%nlsd .NE. 1))&
         CALL stopgm('PERT_KPOINT','Ispin not supported.',& 
         __LINE__,__FILE__)

    IF (tkpts%tkpnt) CALL stopgm&
         ('PERT_KPOINT','K-points are complicated.',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==


    ! Response  wavefunction, computed applying the Hpert
    ! the exact correction is : i * K * cu1
    lcu1 = 2*ncpw%ngw*nstate*3
    ALLOCATE(cu1(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    ! Potentials and related:
    il_v1_nonloc = 2*ncpw%ngw*nstate
    CALL zeroing(v1_nonloc)!,SIZE(v1_nonloc))
    ! NB: In v1_nonloc, we store \op{v1} |psi0> 

    ALLOCATE(dtwnl_dk(2*ncpw%ngw,maxsys%nhxs,maxsys%nsx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dtwnl_dk)!,3*maxsys%nsx*maxsys%nhxs*2*ngw)
    ALLOCATE(ddtwnl_ddk(2*ncpw%ngw,maxsys%nhxs,maxsys%nsx,3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ddtwnl_ddk)!,3*3*maxsys%nsx*maxsys%nhxs*2*ngw)
    ALLOCATE(dfnl_dk(2,ions1%nat,maxsys%nhxs,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(dfnl_dk)!,2*ions1%nat*maxsys%nhxs*nstate*3)
    ALLOCATE(ddfnl_ddk(2,ions1%nat,maxsys%nhxs,nstate,3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ddfnl_ddk)!,2*ions1%nat*maxsys%nhxs*nstate*3*3)

    IF (nkpt%nkpts .EQ. 0) CALL stopgm('PERT_KPOINT_P','NKPTS = 0',& 
         __LINE__,__FILE__)

    CALL dcopy(3*nkpt%nkpts,rrk,1,rk,1)
    CALL dcopy(nkpt%nkpts,wwk,1,wk,1)
    ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
    ! non local potential
    ! ==--------------------------------------------------------------==
    ! FNL  are real because they do not depend on k

    lfnl = nstate*ions1%nat*maxsys%nhxs
    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rnlsm(c0,nstate,1,1,.FALSE.)
    CALL dcopy(lfnl,   fnl,1,  fnl00,1)

    CALL zeroing(sder)!,9)

    ALLOCATE(force1(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(force1)!,SIZE(force1))
    ALLOCATE(sumkikj(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sumkikj)!,9)

    fkpfirst = .TRUE.

    ALLOCATE(e2_nl_c0(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(e2_nl_c0)!,9)

    ! Finite difference 1^ and 2^ of TWNL
    CALL dnlps_dk
    ! Form factors and their 1^ and 2^ derivatives       
    CALL rnl_dk_p(c0,nstate)

    ! Reading RESTA_P: response wfn only
    IF (kpertpar%tkread_c1) THEN
       ! scaling of k-point coordinates wrt BZ primitive vectors
       IF (tkpts%tkscale) THEN
          CALL zeroing(rk_in)!,3*nkpt%nkpts)
          DO ikk=1,nkpt%nkpts
             DO ir=1,3
                rk_in(ir,ikk)=rk(ir,ikk)
             ENDDO
             DO ir=1,3
                rk(ir,ikk)=rk_in(1,ikk)*gvec_com%b1(ir)+rk_in(2,ikk)*gvec_com%b2(ir)&
                     +rk_in(3,ikk)*gvec_com%b3(ir)
             ENDDO
          ENDDO
       ENDIF

       ! if(parent) write(6,'("ENTER ZHRWF_P")')
       ! CALL    ZHRWF_P(7,CU1,NSTATE,TAU0)    
       ! if(parent) write(6,'(" RESTA_P READ")')

       DO k = 1,3
          tag = 'p_'//cdir(k)
          IF (paral%io_parent)&
               WRITE(6,'("READ RESTART FILE NUMBER",I4)') k
          CALL restart_p(cu1(1,1,k),tag,nstate,response_read)
       ENDDO



       GOTO 1000
    ENDIF

    ! 
    ! ==--------------------------------------------------------------==
    ! ==         the basic loops for perturbations                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",22("*"),a,22("*"),/)')&
            '   perturbations    '
    ENDIF
    ! The Response wfn are computed in the 3 cartesian directions
    DO k=1,3
       IF (paral%io_parent)&
            WRITE(6,'(" **** gradient dir:",a4,4x)')&
            cdir(k)
       ! ..   calculate perturbation wrt. coordinate K

       ! ...  calculate the non-local part:
       ! |v1_nonloc> = (d/dK [PP-Projectors]) |0>
       CALL fnl_fkp_p(c0,v1_nonloc,crge%f,nstate,k,fkpfirst)
       ! calculation of |v1_nonloc> = d / d(r_k) |0> + (nl_PP part)
       CALL kin_pert(c0,v1_nonloc,crge%f,nstate,k)

       ! ==--------------------------------------------------------------==
       CALL rwfopt_p(c0,cu1(1,1,k),psi,rhoe,rdummy,&
            eirop,eivps,cdummy,v1_nonloc,&
            z11,nstate,cdummy)
       ! ==--------------------------------------------------------------==

       ! copy the optimized cu1 coefficients into the cu1 array
       ! remind that the exact first order wavefunctions are i * cu1

       sder(k,k) = ener_com%etot

       ! off diagonal elements of the tensor of  second order energy 
       ! 
       IF (k .NE. 1) THEN

          CALL zeroing(force1)!,SIZE(force1))
          CALL h0psi1_p(cu1(1,1,k),c0,rhoe,crge%f,tau0,fion,psi,&
               nstate,force1)
          DO istate=1,nstate
             CALL daxpy(2*ncpw%ngw, z11(istate,istate), cu1(1,istate,k),1,&
                  force1(1,istate),1)
          ENDDO

          DO k2 = 1,k-1

             sder(k2,k) = -energy(cu1(1,1,k2),force1,nstate)-&
                  2._real_8* energy(cu1(1,1,k2),v1_nonloc,nstate)

             CALL mp_sum(sder(k2,k),parai%allgrp)
          ENDDO
       ENDIF
    ENDDO                     ! k

    fkpfirst = .FALSE.
    ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
    ! TEST on symmetry of the tensor
    ! 
    CALL zeroing(force1)!,SIZE(force1))
    CALL h0psi1_p(cu1(1,1,1),c0,rhoe,crge%f,tau0,fion,psi,&
         nstate,force1)
    DO istate=1,nstate
       CALL daxpy(2*ncpw%ngw, z11(istate,istate), cu1(1,istate,1),1,&
            force1(1,istate),1)
    ENDDO

    CALL fnl_fkp_p(c0,v1_nonloc,crge%f,nstate,1,fkpfirst)
    CALL kin_pert(c0,v1_nonloc,crge%f,nstate,1)

    ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
    ! sder(2,1)  should be equal to sder(1,2)

    sder(2,1) = -energy(cu1(1,1,2),force1,nstate)-&
         2._real_8* energy(cu1(1,1,2),v1_nonloc,nstate)

    CALL mp_sum(sder(2,1),parai%allgrp)

    ! sder(3,1)  should be equal to sder(1,3)
    sder(3,1) = -energy(cu1(1,1,3),force1,nstate)-&
         2._real_8* energy(cu1(1,1,3),v1_nonloc,nstate)

    CALL mp_sum(sder(3,1),parai%allgrp)

    CALL zeroing(force1)!,SIZE(force1))

    CALL h0psi1_p(cu1(1,1,2),c0,rhoe,crge%f,tau0,fion,psi,&
         nstate,force1)
    DO istate=1,nstate
       CALL daxpy(2*ncpw%ngw, z11(istate,istate), cu1(1,istate,2),1,&
            force1(1,istate),1)
    ENDDO

    CALL fnl_fkp_p(c0,v1_nonloc,crge%f,nstate,2,fkpfirst)
    CALL kin_pert(c0,v1_nonloc,crge%f,nstate,2)

    ! sder(3,2)  should be equal to sder(2,3)      
    sder(3,2) = -energy(cu1(1,1,3),force1,nstate)-&
         2._real_8* energy(cu1(1,1,3),v1_nonloc,nstate)

    CALL mp_sum(sder(3,2),parai%allgrp)

    ! ================================================================
    ! FINAL RESULTS
    ! ================================================================

    e2_nl_c0(2,1)=e2_nl_c0(1,2)
    e2_nl_c0(3,1)=e2_nl_c0(1,3)
    e2_nl_c0(3,2)=e2_nl_c0(2,3)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a50)')&
            ' energy terms matrix (<cu1(j)|c12(i)>):'
       DO i=1,3

          IF (paral%io_parent)&
               WRITE(6,'("DELTA TENSOR",1x,i2,3(2x,f14.8))')&
               i, sder(i,1),sder(i,2),sder(i,3)

          IF (paral%io_parent)&
               WRITE(6,*) ' '
       ENDDO
       DO i = 1,3
          IF (paral%io_parent)&
               WRITE(6,'("E2_NL_C0 - DIR=",1x,i2,3(2x,f14.8))') i,&
               e2_nl_c0(i,1),e2_nl_c0(i,2),e2_nl_c0(i,3)
       ENDDO
    ENDIF

    ! ==================================================================
    ! calculation of the total perturbaton energy in terms of k_points

    sumkk = 0.0_real_8
    epertkp = 0.0_real_8
    IF (.NOT. tkpts%tmonkp) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'("reciprocal vectors before tkscale")')
          IF (paral%io_parent)&
               WRITE(6,'("b1=",3f14.5)') gvec_com%b1(1),gvec_com%b1(2),gvec_com%b1(3)
          IF (paral%io_parent)&
               WRITE(6,'("b2=",3f14.5)') gvec_com%b2(1),gvec_com%b2(2),gvec_com%b2(3)
          IF (paral%io_parent)&
               WRITE(6,'("b3=",3f14.5)') gvec_com%b3(1),gvec_com%b3(2),gvec_com%b3(3)
       ENDIF
       IF (tkpts%tkscale) THEN
          CALL zeroing(rk_in)!,3*nkpt%nkpts)
          DO ikk=1,nkpt%nkpts
             DO ir=1,3
                rk_in(ir,ikk)=rk(ir,ikk)
             ENDDO
             DO ir=1,3
                rk(ir,ikk)=rk_in(1,ikk)*gvec_com%b1(ir)+rk_in(2,ikk)*gvec_com%b2(ir)&
                     +rk_in(3,ikk)*gvec_com%b3(ir)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF (paral%parent) THEN
       IF (tkpts%tkscale) THEN
          IF (paral%io_parent)&
               WRITE(6,'(a)')&
               "k vector after Tkscale"
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(10x,a,3x,a,2x,a,2x,a,a,a)')&
            "**** num k ",&
            "    kx    ","    ky    ","    kz    ",&
            "      wk    ","****"
       DO ikk=1,nkpt%nkpts
          IF (paral%io_parent)&
               WRITE(6,'(12x,I8,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)')&
               ikk,rk(1,ikk),rk(2,ikk),rk(3,ikk),wk(ikk)
       ENDDO
    ENDIF

    DO ikk=1,nkpt%nkpts
       DO ir1=1,3
          DO ir2=1,3
             sumkikj(ir2,ir1) = sumkikj(ir2,ir1) +&
                  wk(ikk)*rk(ir1,ikk)*rk(ir2,ikk)
          ENDDO
          sumkk = sumkk + wk(ikk)*rk(ir1,ikk)*rk(ir1,ikk)
       ENDDO
    ENDDO
    sumkk = sumkk*parm%tpiba2
    DO istate = 1,nstate
       epertkp = epertkp+0.5_real_8 *  sumkk * crge%f(istate,1)
    ENDDO

    DO ir1=1,3
       DO ir2=1,3
          epertkp = epertkp +&
               parm%tpiba2 * sumkikj(ir2,ir1) *&
               (sder(ir2,ir1) + e2_nl_c0(ir2,ir1))
       ENDDO
    ENDDO

    ! ==-----------------------------------------------------------==
    ! ==  Etot-Egamma = sum_k sum_ab ka*kb*Dab*0.5                 ==
    ! ==-----------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a50)')&
            '========== total perturbation energy:'
       IF (paral%io_parent)&
            WRITE(6,'(1x,a,1x,f16.8)') "EPERTKP=",epertkp*0.5_real_8
       IF (paral%io_parent)&
            WRITE(6,'(1x,a,1x,f16.10)') "sumkk=",sumkk
       IF (paral%io_parent)&
            WRITE(6,'(1x,"sumkikj=")')
       IF (paral%io_parent)&
            WRITE(6,'(8x,3f16.10)')&
            sumkikj(1,1),sumkikj(1,2), sumkikj(1,3)
       IF (paral%io_parent)&
            WRITE(6,'(8x,3f16.10)')&
            sumkikj(2,1),sumkikj(2,2), sumkikj(2,3)
       IF (paral%io_parent)&
            WRITE(6,'(8x,3f16.10)')&
            sumkikj(3,1),sumkikj(3,2), sumkikj(3,3)
    ENDIF

1000 CONTINUE

    DEALLOCATE(force1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sumkikj,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (kpertpar%tpert_hamilt) THEN
       CALL matrix_p(c0,cu1,nstate,psi)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'("        >>>>>>>> K CYCLE <<<<<<<<")')
          IF (paral%io_parent)&
               WRITE(6,'("computation and diagonalization of H(K) ")')
          IF (paral%io_parent)&
               WRITE(6,'("for the following set of k-point        ")')
          IF (paral%io_parent)&
               WRITE(6,'(a,a,a,a,a)')&
               " NK     ","   Kx   ","   Kx   ","   Kx   ",&
               "   Wk   "
          DO ikk = 1,nkpt%nkpts
             IF (paral%io_parent)&
                  WRITE(6,'(i6,2x,4f8.5)')&
                  ikk,rk(1,ikk),rk(2,ikk),rk(3,ikk),wk(ikk)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")')
       ENDIF
       CALL ks_ener_p(cu1,z11,nstate)
    ENDIF

    IF (cntl%tmemchk) THEN
    ENDIF

    DEALLOCATE(dtwnl_dk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddtwnl_ddk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl_dk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddfnl_ddk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! =================================================================
    ! cu1 are now used to increase the Hilbert space

    IF (.NOT. kpertpar%tpert_hamilt) THEN
       IF (paral%io_parent)&
            WRITE(6,'("CHECK IF RESTART")')
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (kpertpar%tk_prep_diag) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>> in PERT_KPOINT:")')
             IF (paral%io_parent)&
                  WRITE(6,'("BUILD_C00 IS TRUE")')
             IF (paral%io_parent)&
                  WRITE(6,'("NUMBER OF STATES OF C00 = ",I7)')&
                  2*nstate
          ENDIF

          CALL up3_p(c0,cu1,tau0,psi,&
               nstate,.FALSE.,.FALSE.)

       ELSEIF (kpertpar%restart_cu1) THEN
          ! IF(PARENT) THEN
          ! WRITE(6,'(">>>> in PERT_KPOINT: CALL ZHWWF_P")')
          ! WRITE(6,'(">>>> THE CU1 ARE TO BE WRITTEN ON RESTA_P")')
          ! ENDIF
          ! CALL ZHWWF_P(5,CU1,NSTATE,TAU0)
          DO k = 1,3
             tag = 'p_'//cdir(k)
             IF (paral%io_parent)&
                  WRITE(6,'(">>>> WRITE RESTRT FILE NUMBER",I4)') k
             CALL restart_p(cu1(1,1,k),tag,nstate,response_write)
          ENDDO
       ELSE
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>> in PERT_KPOINT: ")')
             IF (paral%io_parent)&
                  WRITE(6,'("NO KIND OF RESTART FILE IS WRITTEN")')
             IF (paral%io_parent)&
                  WRITE(6,'("TK_PREP_DIAG and RESTART_CU1 are false")')
          ENDIF
       ENDIF
    ENDIF

    DEALLOCATE(cu1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('    pert_kpoint_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pert_kpoint_p
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE give_scr_pert_kpoint_p(lpert_kpoint_p,tag,nstate)
    INTEGER                                  :: lpert_kpoint_p
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate

    INTEGER :: lforce1, lks_ener_p, lmatrix_p, lnrl_dk_p, lputwnl, lrho, &
      lrhoofr, lrnlsm, lrscpot, lupdrho_p, n2

! ==--------------------------------------------------------------==

    lpert_kpoint_p = 1
    IF (.NOT. response1%tkpert) RETURN

    lpert_kpoint_p=nstate
    ropt_mod%calste=.FALSE.
    n2 = nstate*nstate
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    CALL give_scr_rhoofr(lrho,tag)
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    CALL give_scr_rhoofr_p(lrhoofr,tag)
    CALL give_scr_forces_p(lforce1,tag,nstate)
    CALL give_scr_putwnl_kpert(lputwnl,tag)
    CALL give_scr_rnl_dk_p(lnrl_dk_p,tag,nstate)
    CALL give_scr_up3_p(lupdrho_p,tag,nstate)
    CALL give_scr_matrix_p(lmatrix_p,tag,nstate)
    CALL give_scr_ks_ener_p(lks_ener_p,n2,nstate)
    lrho=lrho+2*ncpw%nhg
    lrhoofr=lrhoofr+2*ncpw%nhg
    lforce1=2*ncpw%ngw*nstate
    lpert_kpoint_p=MAX(lrho,lrscpot,lrnlsm,lrhoofr,&
         lpert_kpoint_p,lforce1,lnrl_dk_p,lupdrho_p,&
         lputwnl,lmatrix_p,lks_ener_p)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_pert_kpoint_p
  ! ==================================================================

END MODULE pert_kpoint_p_utils
