MODULE rwfopt_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE dotp_utils,                      ONLY: dotp
  USE dynit_utils,                     ONLY: dynit
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE hesele_p_utils,                  ONLY: hesele_p
  USE hpsi_utils,                      ONLY: hpsi
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: ortho_p
  USE response_pmod,                   ONLY: dmbi,&
                                             response1,&
                                             response2,&
                                             rho0,&
                                             s_star,&
                                             statemax,&
                                             statemin,&
                                             voa_data
  USE rho1ofr_utils,                   ONLY: rhosofr
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE rpiiint_utils,                   ONLY: rpiiint
  USE simple_model_p_utils,            ONLY: simple_ortho_p
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             group,&
                                             maxsys,&
                                             nacc,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE tpar,                            ONLY: dt2bye
  USE updwf_p_utils,                   ONLY: updwf_p
  USE utils,                           ONLY: unitmx
  USE vofrho_utils,                    ONLY: vofrho
  USE vpsi_p_utils,                    ONLY: allocate_c0real_p,&
                                             release_c0real_p
  USE wrener_utils,                    ONLY: wrener,&
                                             wrprint_wfopt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rwfopt_p


  PUBLIC :: simple_model_p
  PUBLIC :: print_matrix


CONTAINS

  ! ==================================================================
  SUBROUTINE rwfopt_p(c0,c1,psi,rhoe,drhoe,&
       eirop,eivps,v1_loc,v1_nonloc,&
       z11,nstate,eirop1)
    ! ==--------------------------------------------------------------==
    ! COMPUTES the perturbed wavefunctions C1 which are 
    ! the electronic linear response to the perturbation hamiltonian
    ! stored in v1_loc and v1_nonloc.
    ! RETURNS also the first order density (also spin-split) in drhoe, 
    ! which is stored in real space.
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                c1(nkpt%ngwk,*), &
                                                psi(fpar%nnr1,clsd%nlsd)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                drhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg), &
                                                v1_loc(ncpw%nhg), &
                                                v1_nonloc(nkpt%ngwk,*)
    REAL(real_8)                             :: z11(crge%n,crge%n)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: eirop1(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rwfopt_p'

    COMPLEX(real_8), ALLOCATABLE             :: c12(:,:)
    INTEGER                                  :: cg_analytic_original, i, &
                                                ierr, il_c12
    INTEGER, SAVE                            :: reset_counter, &
                                                reset_supercounter
    LOGICAL                                  :: firstcall, reset_cg
    REAL(real_8) :: cg_factor_original, detot, dummy(1), ekin1, ekin2, &
      ekincp, ekinh1, ekinh2, etot0, max_norm_input, norm_input, tcpu, temp1, &
      temp2, thl(2), time1, time2
    REAL(real_8), ALLOCATABLE                :: gde(:,:), pme(:,:), vpp(:)

!rhoe(*),

70  FORMAT ("*  scaling ",a6," with ",f10.3,30(" ")," *")
60  FORMAT ("*  input has zero norm (",e10.2,&
         "< 1.e-8_real_8). Returning zero response.")

    ! ==--------------------------------------------------------------==
    ! == initialization                                               ==
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    nacc = 7
    cg_factor_original=response2%cg_factor
    cg_analytic_original=response1%cg_analytic
    reset_counter=0
    reset_supercounter=0

    response1%trho1 = .TRUE.
    IF ((response1%tnmr .OR. response1%tepr) .OR. response1%tkpert.OR. &
         (response1%tvoa.AND.voa_data%timagpert)) response1%trho1 = .FALSE.
    ! Is there a non-zero first order response density?
    ! For purely imaginary perturbations: No.

    response1%tnonlocal = .TRUE.
    IF (response1%thardness .OR. response1%teigensystem.OR.response1%tfukui) response1%tnonlocal = .FALSE.
    ! Is there a nonlocal contribution (i.e. a non-local perturbation
    ! potential)?

    ener_com%etot = 0._real_8
    etot0 = 0._real_8
    ener_com%ecnstr=0.0_real_8
    ener_com%erestr=0.0_real_8
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%convwf=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    reset_cg=.TRUE.
    firstcall=.NOT.response1%t_initialguess_c1
    ener_com%ebogo=0._real_8
    ! ==--------------------------------------------------------------==


    ! ==-- ORTHOGONALIZATION of the constant part.   -----------------==
    ! DMB
    IF (response1%tnonlocal) THEN
       IF (dmbi%tsimple_model) THEN
          CALL simple_ortho_p(nstate,c0,v1_nonloc,S_star)
       ELSE
          CALL ortho_p(nstate,c0,v1_nonloc)
       ENDIF
    ENDIF
    ! DMB


    ! ==--  SCALING of the input:  -----------------------------------==
    ! dmb: removed scaling if we perform BO dynamics
    ! !    if (TDUMMYATOM.or.(TNMR .OR. TEPR).or. 
    IF (response1%tdummyatom.OR.response1%tepr.OR.&
         (response1%tvoa.AND.voa_data%timagpert).OR.&
         (response1%tinteraction .AND..NOT.cntl%tmdbo)) THEN
       max_norm_input = 0._real_8
       DO i=1,nstate
          norm_input = dotp(ncpw%ngw,v1_nonloc(:,i),v1_nonloc(:,i))
          CALL mp_sum(norm_input,parai%allgrp)
          max_norm_input = MAX(max_norm_input,norm_input)
       ENDDO
       IF ((.NOT.response1%tdummyatom).AND.(ABS(max_norm_input) .LT. 1.e-8_real_8)) THEN
          IF (paral%io_parent)&
               WRITE (6,60) max_norm_input
          CALL zeroing(c1(:,1:nstate))!,ngw*nstate)
          ! DMB       set noresponse flag for SCF-PT calculation
          dmbi%tnoresponse=.TRUE.

          RETURN
       ENDIF
       max_norm_input =  1._real_8 / SQRT(max_norm_input)
       CALL dscal(2*ncpw%ngw*nstate, max_norm_input, v1_nonloc, 1)
       IF (response1%t_initialguess_c1)&
            CALL dscal(2*ncpw%ngw*nstate, max_norm_input, c1, 1)
       IF (paral%io_parent)&
            WRITE (6,70) 'input',max_norm_input
    ENDIF                     ! END scaling NMR

    IF (response1%thardness .OR. response1%teigensystem) THEN
       norm_input = dotp(ncpw%nhg,v1_loc,v1_loc)
       CALL mp_sum(norm_input,parai%allgrp)
       norm_input =  1._real_8 / SQRT(norm_input)
       CALL dscal(2*ncpw%nhg, norm_input, v1_loc, 1)
       IF (response1%t_initialguess_c1)&
            CALL dscal(2*ncpw%ngw*nstate, norm_input, c1,1)
       IF (paral%io_parent)&
            WRITE (6,70)'input',norm_input
    ENDIF
    ! ==--  end of SCALING.  -----------------------------------------==


    ! ==--------------------------------------------------------------==
    ! == initialization                                               ==
    ! ==--------------------------------------------------------------==
    il_c12 = 2*nkpt%ngwk*nstate+8
    ALLOCATE(c12(nkpt%ngwk,il_c12/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c12)!,SIZE(c12))

    ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vpp)!,nkpt%ngwk)
    CALL hesele_p(dt2bye,z11,nstate,vpp)

    IF (response1%diis_p) THEN
       cnti%mdiis = response1%mdiis_p
       ALLOCATE(pme(nkpt%ngwk,cnti%mdiis*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gde(nkpt%ngwk,cnti%mdiis*nstate/4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(pme)!,nkpt%ngwk*nstate*cnti%mdiis)
       CALL zeroing(gde)!,nkpt%ngwk*nstate*cnti%mdiis/4)
    ELSE
       ! avoid fortran 'not allocated' runtime error
       ALLOCATE(pme(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(gde(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF


    ! ==--  REAL SPACE ground state wavefunctions  ---------------==
    IF (response1%tkeeprealspacewfn) THEN
       CALL allocate_c0real_p(c0,psi(:,1),nstate)! inside vpsi_p.F
    ENDIF


    IF (.NOT. response1%t_initialguess_c1) CALL zeroing(c1(:,1:nstate))!,ngw*nstate)
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu=(time2-time1)*1.e-3_real_8
       IF (paral%io_parent)&
            WRITE(6,'(a,f8.2,a8)') ' TIME FOR INITIALIZATION:',&
            tcpu,' SECONDS'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == end initialization                                           ==
    ! ==--------------------------------------------------------------==


    ! ==--------------------------------------------------------------==
    ! ==      the basic loop for wavefunction optimization            ==
    ! ==--------------------------------------------------------------==
    iteropt%nfi=0
    ropt_mod%sdiis  = .TRUE.
    ropt_mod%convwf = .FALSE.
    ener_com%etot   = 0._real_8
    etot0  = 0._real_8
    DO WHILE ((.NOT. ropt_mod%convwf) .AND. (iteropt%nfi .LE. cnti%nomore))
       time1=m_walltime()
       CALL updwf_p(c0,c1,c12,v1_nonloc,z11,&
            v1_loc,eirop1,tau0,fion,pme,gde,vpp,&
            rhoe,drhoe,psi,nstate,firstcall,reset_cg)
       time2=m_walltime()
       tcpu=(time2-time1) * 1.e-3_real_8
       detot = ener_com%etot - etot0
       IF (iteropt%nfi.GE.1) THEN
          CALL wrprint_wfopt(dummy,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
               gemax,cnorm,thl,ropt_mod%convwf,iteropt%nfi,iteropt%nfi)
       ENDIF
       etot0 = ener_com%etot

       ! ERASE the conjugate-direction-memory if the cg-step
       ! goes too much (more than 0.01 a.u.) in the WRONG DIRECTION:
       IF (detot .GT.  1.e-2_real_8 .AND. iteropt%nfi.GT.1 .AND. response1%pcg_p) THEN
          IF (paral%io_parent) WRITE (6,'(65("~"))')
          reset_cg=.TRUE.
          reset_counter = reset_counter+1
          ! a problem:
          IF (reset_counter .GT. 2)&
               CALL dscal(2*ncpw%ngw*nstate,0.3_real_8,c1,1)
          ! a big problem:
          IF (reset_counter .GT. 4) THEN
             IF (response2%cg_factor .GT. 1._real_8) response2%cg_factor = 1._real_8
             response2%cg_factor = 0.8_real_8*response2%cg_factor
             IF (response1%cg_analytic .LT. 3) response1%cg_analytic = 3
             reset_counter = 0
             reset_supercounter = reset_supercounter + 1
             IF (reset_supercounter .GT. 5) THEN
                IF (paral%io_parent) THEN
                   WRITE (6,*) '*** TOO MANY RESETS. EXITING.'
                ENDIF
                soft_com%exsoft = .TRUE.
             ENDIF
          ENDIF
       ENDIF

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf .OR. soft_com%exsoft) GOTO 100

       iteropt%nfi=iteropt%nfi+1
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==     end of main loop                                         ==
    ! ==--------------------------------------------------------------==


    IF (paral%io_parent.AND.gemax.GT.response2%tolog_p) THEN
       WRITE(6,'(1x,64("!"))')
       WRITE(6,'(" !!",a,t64,"!!")')&
            ' RWF1OPT| THE MAXIMUM NUMBER OF STEP IS REACHED'
       WRITE(6,'(" !!",a,f10.6,a,t64,"!!")')&
            '         BUT NO CONVERGENCE (d E / d phi_1 =',gemax,')'
       WRITE(6,'(1x,64("!"))')
    ENDIF
100 CONTINUE

    IF (paral%parent) CALL wrener

    ! dmb scaling removed for BO dynamics
    ! !    if (TDUMMYATOM.or.(TNMR.OR.TEPR).or.
    IF (response1%tdummyatom.OR.response1%tepr.OR.&
         (response1%tvoa.AND.voa_data%timagpert).OR.&
         (response1%tinteraction.AND..NOT.cntl%tmdbo)) THEN
       CALL dscal(2*ncpw%ngw*nstate, 1._real_8/max_norm_input, c1,1)
       ! if (parent) write (6,70) 'output',1._real_8/max_norm_input
       ! DMB: printing the correction energy
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,70) 'output',1._real_8/max_norm_input
          ener_com%etot=ener_com%etot*(1._real_8/max_norm_input)**2
          IF ((response1%tinteraction) .AND.paral%io_parent)&
               WRITE(6,*) 'TOTAL ENERGY CORRECTION =',ener_com%etot
       ENDIF
    ENDIF
    IF (response1%thardness .OR. response1%teigensystem) THEN
       norm_input = 1._real_8/norm_input
       CALL dscal(2*ncpw%ngw*nstate, norm_input, c1,1)
       CALL dscal(2*ncpw%ngw*nstate, norm_input, v1_loc,1)
       IF (paral%io_parent)&
            WRITE (6,70) 'output',norm_input
    ENDIF

    ! ==--  RELEASE MEMORY...   ----------------------------------==
    IF (response1%diis_p) THEN           ! last allocated = first de-allocated
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(pme,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(pme,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c12,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (response1%tkeeprealspacewfn) THEN
       CALL release_c0real_p ! routine is located inside vpsi_p.f
    ENDIF
    ! ==--------------------------------------------------------------==
    response2%cg_factor=cg_factor_original
    response1%cg_analytic=cg_analytic_original
    RETURN
  END SUBROUTINE rwfopt_p
  ! ==--------------------------------------------------------------==


  ! ==================================================================
  SUBROUTINE simple_model_p(c0,c1,psi,rhoe,nstate,eirop,&
       eivps)
    ! ==================================================================
    ! This subroutine calculates the interaction energy between molecules
    ! using non-orthogonal orbitals and a simplified model for perturbation
    ! theory. So far, this requires that each molecule in the system was
    ! previously subjected to a gas phase calculation where both its ground
    ! state wave function (W0) and its density-dependent energy terms 
    ! (EHT+EPSEU+EXC) are saved (ref-wannier-xx.dat and Erhomol files).
    ! DMB.

    ! PARALLEL
    ! ==--------------------------------------------------------------==
    ! Arguments & common variables
    COMPLEX(real_8)                          :: psi(maxfft*group%nogrp,1)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate), &
                                                eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'simple_model_p'

    COMPLEX(real_8)                          :: cdummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), h0_c0(:,:), &
                                                h0_c1(:,:), h0_Eij_c0(:,:), &
                                                total_force(:,:), v1_c0(:,:)
    INTEGER                                  :: i, ierr, j, nfto
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: E0tot, E_rho_diff, E_SUM_h0, &
                                                E_SUM_h0_diag, &
                                                E_SUM_h0_off_diag, &
                                                E_SUM_rho_mol, Efunct
    REAL(real_8), ALLOCATABLE :: E_rho_mol(:), H_matrix(:,:), &
      my_unit_matrix(:,:), one_m_S2_matrix(:,:), one_m_S_matrix(:,:), &
      rho1(:), rho_mol(:), S_matrix(:,:), vrho0(:,:)

!rhoe(*)

    IF (cntl%ttau) CALL stopgm("simple_model_p","No Tau functionals",& 
         __LINE__,__FILE__)
    ! FIRST INITIALISATIONS

    ! Zeroing all energy terms to test the calculation
    ener_com%etot=0
    ener_com%ekin=0
    ener_com%eht=0
    ener_com%epseu=0
    ener_com%enl=0
    ener_com%exc=0

    ! initialising the arrays
    ALLOCATE(E_rho_mol(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(S_matrix(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(h0_c0(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(h0_Eij_c0(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(h0_c1(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(v1_c0(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(total_force(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(H_matrix(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(one_m_S_matrix(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(one_m_S2_matrix(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(my_unit_matrix(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rho_mol(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rho1(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vrho0(fpar%nnr1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! S_star is defined in response_p.inc...
    ALLOCATE(S_star(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! NOW WE CAN CALCULATE....
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*) '|   SIMPLE INTERACTION MODEL    |'
       IF (paral%io_parent)&
            WRITE(6,*) '| USING NON-ORTHOGONAL ORBITALS |'
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '***  CALCULATING E0  ***'
    ENDIF
    ! trying the manual way of calculating the energy...(taken from noforce)
    ! FIRST DENSITY DEPENDENT TERMS ONLY (using non-orthogonal density..)
    ! ion-ion interaction
    nfto=3*maxsys%nax*maxsys%nsx
    CALL zeroing(fion)!,nfto)
    CALL rpiiint(ener_com%esr,tau0,fion,iteropt%iesr,.FALSE.)

    ! calculate the density and the wrong KE (not orthogonal..)
    ! building the density state by state
    CALL zeroing(rho0)!,nnr1)

    E_SUM_rho_mol=0._real_8

    ! opening erhomol file
    IF (paral%io_parent)&
         CALL fileopen(500,'Erhomol',fo_old,ferror)
    ! Did the file open OK?
    IF (ferror) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Error opening the file Erhomol...'
          CALL stopgm('simple_model_p','file open error',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF

    ! constructing total density from molecular fragments
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'Isolated molecular densities'
    ENDIF

    DO j=1,dmbi%nmol
       ! first zeroing the molecular density
       CALL zeroing(rho_mol)!,nnr1)
       ! loops over states of one molecule
       DO i=STATEmin(j),STATEmax(j)
          ! calculating density for molecule j state by state
          CALL rhosofr(c0(:,i),rhoe(:,1),psi(:,1))
          ! adding state density rhoe times occupation number to molecular density rho_mol
          CALL daxpy(fpar%nnr1,crge%f(i,1),rhoe,1,rho_mol,1)
          ! adding state density rhoe times occupation number to total density rho0 
          CALL daxpy(fpar%nnr1,crge%f(i,1),rhoe,1,rho0,1)
       ENDDO

       ! calculating molecular E[rho_mol] terms

       ! for now, we just read them from a file...
       ! Reading energy of molecule j
       IF (paral%io_parent)&
            READ(500,*) E_rho_mol(j)

       ! calculate EHT, EPSEU, EXC
       ! CALL VOFRHO(TAU0,FION,rho_mol,PSI,.false.,.false.)
       ! storing the corresponding energy terms..
       ! E_rho_mol(j)=EHT+EPSEU+EXC
       ! WRITE(6,*)'EHT   (TEST   W0)=        ',EHT,' A.U.'
       ! WRITE(6,*)'EPSEU (TEST   W0)=        ',EPSEU,' A.U.'
       ! WRITE(6,*)'EXC   (TEST   W0)=        ',EXC,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*) 'molecule',j,'E[rho_mol] = ',E_rho_mol(j)

       ! then adding the energies together 
       E_SUM_rho_mol=e_sum_rho_mol+E_rho_mol(j)
    ENDDO

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'Sum of E[rho_mol] =',E_SUM_rho_mol
       IF (paral%io_parent)&
            CALL fileclose(500)
    ENDIF
    ! NOW USING TOTAL DENSITY
    ! first copy rho0 zero to preserve it for later...
    CALL dcopy(fpar%nnr1,rho0,1,vrho0,1)

    ! calculate the density dependent parts.. (EHT,EPSEU,EXC)
    CALL vofrho(tau0,fion,vrho0,psi,.FALSE.,.FALSE.)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'TOTAL DENSITY DEPENDENT TERMS'
       IF (paral%io_parent)&
            WRITE(6,*) '***** ENERGY DECOMPOSITION *****'
       IF (paral%io_parent)&
            WRITE(6,*) 'EHT   (       W0)= ',ener_com%eht,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*) 'EPSEU (       W0)= ',ener_com%epseu,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*) 'EXC   (       W0)= ',ener_com%exc,' A.U.'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'E[rho0] =',ener_com%eht+ener_com%epseu+ener_com%exc
    ENDIF
    ! END  OF DENSITY ONLY TERMS

    ! new non-orthogonal correction thing...
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'WAVEFUNCTION DEPENDENT TERMS'
    ENDIF

    ! Calculating H[rho0]|c0i>  CS is a tmp array
    CALL hpsi(c0,h0_c0,cs,vrho0,psi(:,1),nstate,1,1)

    ! multiplying h0_c0 by the occupation number
    DO i=1,nstate
       CALL dscal(2*ncpw%ngw,crge%f(i,1),h0_c0(1,i),1)
    ENDDO

    ! then H_matrix=<c0j|h0_c0>
    CALL ovlap(nstate,H_matrix,c0,h0_c0)

    ! changing sign of H matrix (force->hamiltonian)
    CALL dscal(nstate*nstate,-1._real_8,H_matrix,1)

    ! now printing...
    CALL print_matrix(H_matrix,nstate,nstate,' H matrix ')

    ! Calculating <c0i|c0j> for the non-orthogonal states
    CALL ovlap(nstate,S_matrix,c0,c0)

    ! calculating (1-S) matrix
    ! creating a unit matrix
    CALL unitmx(my_unit_matrix,nstate)
    ! copying unit matrix to one_m_S_matrix
    CALL dcopy(nstate*nstate,my_unit_matrix,1,one_m_S_matrix,1)
    ! then adding -S to unit matrix
    CALL daxpy(nstate*nstate,-1._real_8,S_matrix,1,one_m_S_matrix,1)
    ! printing
    CALL print_matrix(one_m_S_matrix,nstate,nstate,'1-S matrix')

    ! calculating (1-S)(1-S) matrix
    CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,one_m_S_matrix,&
         NSTATE,one_m_S_matrix,NSTATE,0._real_8,one_m_S2_matrix,NSTATE)
    ! printing
    CALL print_matrix(one_m_S2_matrix,nstate,nstate,'(1-S)2 mat')

    ! Building S_star matrix which has (1-S)**2 on its diagonal
    ! and (1-S) on the off-diagonal elements
    ! copying one_m_S_matrix to S_star
    CALL dcopy(nstate*nstate,one_m_S_matrix(1,1),1,S_star(1,1),1)

    ! now affecting the (1-S)**2 elements
    DO i=1,nstate
       S_star(i,i)=one_m_S2_matrix(i,i)
    ENDDO
    ! printing
    CALL print_matrix(S_star,nstate,nstate,'S_star    ')

    ! Now summing the energy terms
    E_SUM_h0_diag=0._real_8
    E_SUM_h0_off_diag=0._real_8
    DO i=1,nstate
       DO j=1,nstate
          ! off diagonal term (1-S)ij <psi0j|H[rho(0)]|psi0i>
          IF (i .NE. j) THEN
             E_SUM_h0_off_diag=e_sum_h0_off_diag+&
                  one_m_S_matrix(i,j)*H_matrix(j,i)
          ENDIF
       ENDDO
       ! diagonal term [(1-S)(1-S)]ii <psi0i|H[rho(0)]|psi0i>
       E_SUM_h0_diag=e_sum_h0_diag+&
            one_m_S2_matrix(i,i)*H_matrix(i,i)
    ENDDO

    ! wavefunction terms
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'E[H0]ii (diagonal)=',E_SUM_h0_diag
       IF (paral%io_parent)&
            WRITE(6,*) 'E[H0]ij (off diagonal)=',E_SUM_h0_off_diag
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    ! sum of wavefunction correction terms
    E_SUM_h0=E_SUM_h0_diag+E_SUM_h0_off_diag
    IF (paral%io_parent)&
         WRITE(6,*) 'E[H0] (total)=',E_SUM_h0

    ! Density energy difference
    E_rho_diff=ener_com%eht+ener_com%epseu+ener_com%exc-E_SUM_rho_mol
    IF (paral%io_parent)&
         WRITE(6,*) 'E[rho0]-Sum(E[rho_mol])=',E_rho_diff

    ! TOTAL interaction energy E0
    E0tot=E_SUM_h0+E_rho_diff
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'E0tot = E[H0]+E[rho0]-Sum(E[rho_mol]) =',e0tot
       ! the same in mH
       IF (paral%io_parent)&
            WRITE(6,*) 'Interaction energy (mH)=',1000*E0tot
    ENDIF

    ! Done with the first part...
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '*** FINISHED WITH E0 ***'
    ENDIF

    ! OPTIMISATION FOR C1
    ! it s easier to let the usual PT to finish the job....
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '*** CALCULATING DFPT CORRECTIONS ***'
    ENDIF

    ! Use Daniel s PT to optimise C1
    CALL rwfopt_p(c0,c1,psi,rho0,rho1,eirop,eivps,cdummy,&
         h0_c0,H_matrix,nstate,cdummy)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '*** FINISHED WITH DFPT CORRECTIONS ***'
    ENDIF
    ! ETOT now contains E1+E2
    ! END OF OPTIM FOR C1

    ! TOTAL interaction energy E0+(E1+E2)
    Efunct=E0tot+ener_com%etot
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'Efunct =',efunct
       ! the same in mH
       IF (paral%io_parent)&
            WRITE(6,*) 'EFI energy (mH)=',1000*Efunct
    ENDIF


    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*) '|   EXITING...                  |'
       IF (paral%io_parent)&
            WRITE(6,*) '|   SIMPLE INTERACTION MODEL    |'
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------+'
    ENDIF
    ! simple model finishes here...
    ! freeing the memory
    DEALLOCATE(E_rho_mol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(S_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(h0_c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(h0_Eij_c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(h0_c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(H_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(one_m_S_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(one_m_S2_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(S_star,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(my_unit_matrix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho_mol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rho1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vrho0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE simple_model_p
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE optim_diis_print(step,grad_max,c_norm,etot,detot,tcpu)
    ! ==--------------------------------------------------------------==
    ! prints out step-by-step details of the convergence of cntl%diis minimisation
    INTEGER                                  :: step
    REAL(real_8)                             :: grad_max, c_norm, etot, &
                                                detot, tcpu

    IF ((step.EQ.1).AND.paral%io_parent)&
         WRITE(6,'(A,A)')&
         '   NFI      GEMAX     CNORM',&
         '         ETOT(2)        DETOT     TCPU'
    IF (paral%io_parent)&
         WRITE(6,'(I6,F11.6,F10.6,F16.8,G13.4,F9.3)')&
         step,grad_max,c_norm,etot,detot,tcpu
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE optim_diis_print
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE print_matrix(a,d1,d2,title)
    ! ==--------------------------------------------------------------==
    ! Prints out the matrix A(d1,d2) and its name (max 10 character!!)
    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: d1, d2
    REAL(real_8)                             :: a(d1,d2)
    CHARACTER(len=10)                        :: title

    INTEGER                                  :: i, j

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) title
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,11(1X,I12))') 'ISTATE',(j,j=1,d2)
       DO i=1,d1
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,11(1X,F12.6))') i,(a(i,j),j=1,d2)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    RETURN
  END SUBROUTINE print_matrix
  ! ==================================================================

END MODULE rwfopt_p_utils
