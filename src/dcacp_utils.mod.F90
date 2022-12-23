! ==================================================================
! Provides: Routines to add dispersion-corrected atom-centred
!           potentials; re-using the F-channel of GTH or MT
!           pseudo-potentials (quick and dirty).
!
! dcacp_init: Adds DCACP to f-channels, does broadcasting
!             and reporting
!
!                             11.06.2017 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE dcacp_utils
  USE atom,                            ONLY: atom_common,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE cp_xc_utils,                     ONLY: cp_xc_functional_env,&
                                             cp_xc_env_t
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  ! Remove once NEWCODE and OLDCODE are removed for good
  USE func,                            ONLY: func1,&
                                             mfxcx_is_slaterx,&
                                             mfxcc_is_pade,&
                                             mgcx_is_becke88,&
                                             mgcc_is_lyp,&
                                             mgcc_is_perdew86
  ! End remove
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: default_string_length,&
                                             real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             nlps_com,&
                                             rgh,&
                                             wgh
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE recpnew_utils,                   ONLY: ckgrid,&
                                             tgrid
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsp,&
                                             maxsys
  USE timer,                           ONLY: tiset,&
                                             tihalt
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  !
  ! Easy access to the library
  !
  INTEGER, PARAMETER, PRIVATE :: H = 1, He = 2, &
     Li = 3, Be = 4, B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10, &
     Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18, &
     K = 19, Ca = 20, &
     Sc = 21, Ti = 22, V  = 23, Cr = 24, Mn = 25, Fe = 26, Co = 27, Ni = 28, Cu = 29, Zn = 30, &
     Ga = 31, Ge = 32, As = 33, Se = 34, Br = 35, Kr = 36, &
     Rb = 37, Sr = 38, &
     Y  = 39, Zr = 40, Nb = 41, Mo = 42, Tc = 43, Ru = 44, Rh = 45, Pd = 46, Ag = 47, Cd = 48, &
     In = 49, Sn = 50, Sb = 51, Te = 52, I  = 53, Xe = 54, &
     Cs = 55, Ba = 56, &
     La = 57, Ce = 58, Pr = 59, Nd = 60, Pm = 61, Sm = 62, Eu = 63, &
     Gd = 64, Tb = 65, Dy = 66, Ho = 67, Er = 68, Tm = 69, Yb = 70, Lu = 71, &
     Hf = 72, Ta = 73, W  = 74, Re = 75, Os = 76, Ir = 77, Pt = 78, Au = 79, Hg = 80, &
     Tl = 81, Pb = 82, Bi = 83, Po = 84, At = 85, Rn = 86, &
     Fr = 87

  INTEGER, PARAMETER, PRIVATE         :: l_dcacp     =  4
  INTEGER, PARAMETER, PRIVATE         :: max_dcacp   = Rn
  INTEGER, PARAMETER, PRIVATE         :: output_unit =  6
  CHARACTER(len=*), PARAMETER, PUBLIC :: dcacp_custom_tag = 'Values from custom input'

  !
  ! Library (parameters)
  !
  TYPE, PRIVATE :: sigma_t
     REAL(real_8)                         :: sigma_1  = 0.0_real_8
     REAL(real_8)                         :: sigma_2  = 0.0_real_8
     CHARACTER(len=80)                    :: citation = ' ' 
  END TYPE sigma_t

  !
  ! DCACP input
  !
  TYPE, PRIVATE :: dcacp_env_t
     LOGICAL, DIMENSION(maxsp)            :: no_contribution  = .false.
     LOGICAL                              :: has_custom_sigma = .false.
     LOGICAL                              :: include_metals   = .false.
     TYPE(sigma_t), DIMENSION(max_dcacp)  :: data_from_input
  END TYPE dcacp_env_t

  !
  ! DCACP settings
  !
  TYPE, PRIVATE :: dcacp_t
     INTEGER                              :: l_dcacp          = l_dcacp
     LOGICAL, DIMENSION(maxsp)            :: no_contribution  = .false.
     LOGICAL, DIMENSION(4,maxsp)          :: skip             = .false.
     LOGICAL                              :: active           = .false.
     LOGICAL                              :: is_compatible    = .false.
     CHARACTER(len=default_string_length) :: xc_functional    = ' '
     TYPE(sigma_t), DIMENSION(max_dcacp)  :: element
  END TYPE dcacp_t

  !
  TYPE(dcacp_t), PUBLIC              :: dcacp
  TYPE(dcacp_env_t), PUBLIC          :: dcacp_env

  PUBLIC  :: dcacp_init

  PRIVATE :: dcacp_report
  PRIVATE :: bcast_dcacp_pre_init
  PRIVATE :: bcast_dcacp_post_init
  PRIVATE :: get_dcacp_potential
  PRIVATE :: get_dcacp_wavefunction
  PRIVATE :: set_dcacp_library
  PRIVATE :: set_dcacp_projectors
  PRIVATE :: no_contribution

CONTAINS

  ! ==================================================================
  SUBROUTINE dcacp_init()
    ! ==--------------------------------------------------------------==
    ! Adds DCACP, re-using pseudo-potential arrays and algorithms
    ! for ease of implementation

    CHARACTER(len=*), PARAMETER :: procedureN = 'dcacp_init'

    INTEGER                     :: titag

    CALL tiset(procedureN,titag)

    CALL bcast_dcacp_pre_init()
    CALL set_dcacp_library(dcacp,dcacp_env,cp_xc_functional_env)
    CALL set_dcacp_projectors(dcacp)
    CALL dcacp_report(dcacp,dcacp_env)
    CALL bcast_dcacp_post_init()

    CALL tihalt(procedureN,titag)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE dcacp_init
  ! ==================================================================

  ! ==================================================================
  ! Internal routines
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE set_dcacp_library( dcacp,dcacp_env,cp_xc_env )
    ! ==--------------------------------------------------------------==
    ! Selects a library depending on the xc functional used and/or
    ! customisation options

    CHARACTER(len=*), PARAMETER     :: procedureP = 'set_dcacp_library'
    CHARACTER(len=*), PARAMETER     :: procedureN = procedureP

    TYPE(dcacp_t), INTENT(inout)    :: dcacp
    TYPE(dcacp_env_t), INTENT(in)   :: dcacp_env
    TYPE(cp_xc_env_t), INTENT(in)   :: cp_xc_env

    INTEGER                         :: isp, z, titag


    CALL tiset(procedureN,titag)

    dcacp%l_dcacp = l_dcacp

    IF (cntl%use_xc_driver) THEN
       IF (cp_xc_env%n_funcs == 1) THEN
          IF (cp_xc_env%funcs(1) == 'GGA_XC_BLYP') THEN
             dcacp%xc_functional = cp_xc_env%funcs(1)
          ELSEIF (cp_xc_env%funcs(1) == 'GGA_XC_BP86' .OR. &
                  cp_xc_env%funcs(1) == 'GGA_XC_BP') THEN
             dcacp%xc_functional = cp_xc_env%funcs(1)
          ENDIF
       ELSEIF (cp_xc_env%n_funcs == 2) THEN
          IF ( any(cp_xc_env%funcs(:) == 'GGA_X_B88') .AND. &
               any(cp_xc_env%funcs(:) == 'GGA_C_LYP') ) THEN
             dcacp%xc_functional = 'GGA_XC_BLYP'
          ELSEIF ( any(cp_xc_env%funcs(:) == 'GGA_X_B88') .AND. &
               any(cp_xc_env%funcs(:) == 'GGA_C_P86') ) THEN
             dcacp%xc_functional = 'GGA_XC_BP86'
          ENDIF
       ENDIF
    !
    ! Clean up once OLDCODE and NEWCODE are removed
    !
    ELSE
       IF ( func1%mfxcx == mfxcx_is_slaterx .AND. &
            func1%mgcx  == mgcx_is_becke88 .AND. &
            func1%mgcc  == mgcc_is_lyp .AND. &
            func1%mhfx  == 0 ) THEN
          dcacp%xc_functional = 'GGA_XC_BLYP'
       ELSEIF ( (func1%mfxcx == mfxcx_is_slaterx .OR. &
                 func1%mfxcc == mfxcc_is_pade) .AND. &
            func1%mgcx  == mgcx_is_becke88 .AND. &
            func1%mgcc  == mgcc_is_perdew86 .AND. &
            func1%mhfx  == 0 ) THEN
          dcacp%xc_functional = 'GGA_XC_BP86'
       ENDIF
    ENDIF
    !
    ! End of future clean up
    !
    CALL construct_library( dcacp,dcacp_env ) 
    !
    IF (dcacp_env%has_custom_sigma) THEN
       !
       ! Copy values from input if needed
       !
       DO z=1,max_dcacp
          IF ( trim(adjustl(dcacp_env%data_from_input(z)%citation)) == dcacp_custom_tag ) THEN
             dcacp%element(z)%sigma_1  = dcacp_env%data_from_input(z)%sigma_1
             dcacp%element(z)%sigma_2  = dcacp_env%data_from_input(z)%sigma_2
             dcacp%element(z)%citation = dcacp_env%data_from_input(z)%citation
          ENDIF
       ENDDO
    ENDIF
    
    !
    ! Remove contributions
    !
    FORALL (isp=1:size(dcacp%no_contribution)) &
       dcacp%no_contribution(isp) = no_contribution(dcacp_env,isp)

    CALL tihalt(procedureN,titag)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_dcacp_library
  ! ==================================================================
  SUBROUTINE set_dcacp_projectors( dcacp )
    ! ==--------------------------------------------------------------==
    ! Wrapper that sets up the DCACP projectors

    CHARACTER(len=*), PARAMETER     :: procedureP = 'set_dcacp_projectors'
    CHARACTER(len=*), PARAMETER     :: procedureN = procedureP

    INTEGER                         :: isp

    TYPE(dcacp_t), INTENT(inout)    :: dcacp

    !
    ! As in previous routines, many of the parameters are still held by io_parent
    ! only (broadcasting is done in setsys)
    !
    IF (paral%io_parent) THEN
       !
       ! Add DCACP to every species & change values of nghtol, nghcom
       !
       DO isp=1,ions1%nsp
          IF (dcacp%no_contribution(isp)) CYCLE
          !
          CALL reset_compatibility( dcacp )
          CALL check_compatibility( dcacp,isp )
          CALL check_availability( dcacp,isp )
          !
          IF (pslo_com%tnum(isp)) THEN
             CALL set_numerical_dcacp( dcacp,isp )
          ELSEIF (sgpp1%tsgp(isp)) THEN
             CALL set_analytical_dcacp( dcacp,isp )
          ELSE
             CALL stopgm(procedureN,'DCACP are only implemented in combination with '//&
                                    'Goedecker or numerical pseudo-potentials.', &
                                     __LINE__, __FILE__)
          ENDIF
       ENDDO
    ENDIF

    !
    ! This parameter is used in the PP projector routines - if we get through here,
    ! everything should be set up and be fine
    !
    dcacp%active = .TRUE.

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE set_numerical_dcacp( dcacp,isp )
    
       CHARACTER(len=*), PARAMETER  :: procedureN = procedureP//'/set_numerical_dcacp'

       INTEGER, INTENT(in)          :: isp
       INTEGER                      :: il, iv, l, m

       TYPE(dcacp_t), INTENT(inout) :: dcacp

       IF (.NOT. dcacp%is_compatible) CALL stopgm(procedureN,'Tried calling routine with incompatible '//&
                                                  'DCACP and pseudo-potential parameters',&
                                                  __LINE__, __FILE__)

       IF (.NOT. dpot_mod%tkb(isp)) CALL stopgm(procedureN,'DCACP only supports Kleinman-Bylander '//&
                                                'integration when using numerical pseudo-potentials.', &
                                                __LINE__, __FILE__)
       !
       ! Kleinman-Bylander:
       ! In case there is no wavefunction (e.g. for H)
       !
       IF (dpot_mod%lmax(isp) == 1) THEN
           !
           ! Copy potential mesh into wavefunction mesh
           ! (already logarithmic)
           !
           rw(:,isp)              = rv(:,isp)
           atom_common%clogw(isp) = atom_common%clogvp(isp)
           atom_common%meshw(isp) = atom_common%meshvp(isp)
       ENDIF
       CALL get_dcacp_potential( dcacp, isp, vr(1:atom_common%meshvp(isp),:,:) )
       CALL get_dcacp_wavefunction( dcacp, isp, rw(1:atom_common%meshw(isp),isp),&
                                    rps(1:atom_common%meshw(isp),:,:) )
       !
       DO il=(dpot_mod%lmax(isp)+1), 3
          dcacp%skip(il,isp) = .TRUE.
       ENDDO
       !
       ! Decreasing the nbr. of non-local projectors as in recpnew
       !
       ! Reset ngh first
       !
       nlps_com%ngh(isp) = 0
       iv = 0
       DO il=1,dcacp%l_dcacp
          IF ( il == dpot_mod%lloc(isp)  .OR. &
               il == dpot_mod%lskip(isp) .OR. &
               dcacp%skip(il,isp) ) CYCLE
          l = il-1
          DO m=1,2*l+1
             iv = iv+1
             nghtol(iv,isp) = l
             nghcom(iv,isp) = l*l + m
          ENDDO
          nlps_com%ngh(isp) = nlps_com%ngh(isp) + (2*l+1)
       ENDDO

    END SUBROUTINE set_numerical_dcacp
    ! ==--------------------------------------------------------------==
    SUBROUTINE set_analytical_dcacp( dcacp,isp )
   
       CHARACTER(len=*), PARAMETER :: procedureN = procedureP//'/set_analytical_dcacp'

       INTEGER, INTENT(in)          :: isp
       INTEGER                      :: iv, l, lp, m, k

       TYPE(dcacp_t), INTENT(in)    :: dcacp
 
       IF (.NOT. dcacp%is_compatible) CALL stopgm(procedureN,'Tried calling routine with incompatible '//&
                                                  'DCACP and pseudo-potential parameters',&
                                                  __LINE__, __FILE__)

       iv=0
       lp=0
       !
       ! Get iv for the projectors that are already there
       !
       DO l=1,dpot_mod%lmax(isp)-1
          DO m=1,2*l-1
             lp = lp+1
             DO k=1,sgpp2%npro(l,isp)
                iv = iv+1
             ENDDO
          ENDDO
       ENDDO
       !
       ! Fill empty projectors with zeros
       !
       DO l=dpot_mod%lmax(isp),dcacp%l_dcacp-1
          sgpp2%rcnl(l,isp) = 0.0_real_8
          sgpp2%npro(l,isp) = 0
          DO m=1,2*l-1
             lp = lp+1
          ENDDO
       ENDDO
       !
       ! Add DCACP to f-channel
       !
       sgpp2%rcnl(dcacp%l_dcacp,isp)     = dcacp%element(ions0%iatyp(isp))%sigma_1
       sgpp2%npro(dcacp%l_dcacp,isp)     = 1
       sgpp2%hlsg(1,1,dcacp%l_dcacp,isp) = dcacp%element(ions0%iatyp(isp))%sigma_2
       !
       DO m=1,(2*dcacp%l_dcacp-1)
          lp = lp+1
          iv = iv+1
          nghtol(iv,isp)      = dcacp%l_dcacp-1
          sgpp2%lpval(iv,isp) = lp
          sgpp2%lfval(iv,isp) = 1
       ENDDO
       !
       ! Update (in contrast to numerical PPs, these need to be hacked here already)
       !
       dpot_mod%lmax(isp) = dcacp%l_dcacp + 1 
       dpot_mod%lloc(isp) = dcacp%l_dcacp + 1 
       dpot_mod%lskip(isp)= dcacp%l_dcacp + 6 ! just to be sure 
       nlps_com%ngh(isp)  = iv

    END SUBROUTINE set_analytical_dcacp
    ! ==--------------------------------------------------------------==
    SUBROUTINE check_availability( dcacp,isp )

       CHARACTER(len=*), PARAMETER  :: procedureN = procedureP//'/check_availability'

       INTEGER, INTENT(in)          :: isp

       TYPE(dcacp_t), INTENT(inout) :: dcacp

       CHARACTER(len=3)             :: z
       CHARACTER(len=80)            :: or_zero = ' '
    
       LOGICAL                      :: is_compatible

       is_compatible = .FALSE.

       IF (dcacp%element(ions0%iatyp(isp))%sigma_1 == 0.0_real_8 .OR. &
           dcacp%element(ions0%iatyp(isp))%sigma_2 == 0.0_real_8 ) THEN
          WRITE(z,'(I3)') ions0%iatyp(isp)
          IF (dcacp%element(ions0%iatyp(isp))%citation == dcacp_custom_tag) THEN
             or_zero = 'or input parameters are zero'
          ELSE
             or_zero = ''
          ENDIF
          CALL stopgm(procedureN,'Element with Z='//trim(adjustl(z))//' unavailable '//&
                       'in internal DCACP library '//trim(adjustl(or_zero)),&
                       __LINE__, __FILE__)
       ELSE
          is_compatible = .TRUE.
       ENDIF

       dcacp%is_compatible = (is_compatible .AND. dcacp%is_compatible)

    END SUBROUTINE check_availability
    ! ==--------------------------------------------------------------==
    SUBROUTINE check_compatibility( dcacp, isp )

       CHARACTER(len=*), PARAMETER  :: procedureN = procedureP//'/check_compatibility'

       TYPE(dcacp_t), INTENT(inout) :: dcacp
       INTEGER, INTENT(in)          :: isp

       IF (dcacp%is_compatible) CALL stopgm(procedureN,'Wrong DCACP compatibility status '//&
                                            'upon entering the routine',&
 __LINE__, __FILE__)

       IF (dpot_mod%lmax(isp) >= dcacp%l_dcacp .and. pslo_com%tnum(isp)) THEN 
          CALL stopgm(procedureN, 'DCACP are not implemented for LMAX > 3 due to '//&
                                  'memory management simplifications.', &
 __LINE__, __FILE__)
       ELSE IF (dpot_mod%lmax(isp) > dcacp%l_dcacp .and. sgpp1%tsgp(isp)) THEN
          CALL stopgm(procedureN, 'DCACP are not implemented for LMAX > 3 due to '//&
                                  'memory management simplifications.', &
  __LINE__, __FILE__)
       ENDIF

       dcacp%is_compatible = .TRUE.

    END SUBROUTINE check_compatibility
    ! ==--------------------------------------------------------------==
    SUBROUTINE reset_compatibility( dcacp )

       TYPE(dcacp_t), INTENT(inout) :: dcacp

       dcacp%is_compatible = .FALSE.

    END SUBROUTINE reset_compatibility
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_dcacp_projectors
  ! ==================================================================
  SUBROUTINE dcacp_report( dcacp, dcacp_env )
    ! ==--------------------------------------------------------------==
    ! Writes some more or less nice output

    CHARACTER(len=80)            :: index_list

    INTEGER                      :: z, isp
    LOGICAL                      :: no_contribution_at_all

    TYPE(dcacp_t), INTENT(in)    :: dcacp
    TYPE(dcacp_env_t), &
                   INTENT(in)    :: dcacp_env


    IF (paral%io_parent) THEN
       WRITE(output_unit,'(A)')
       WRITE(output_unit,'(1X,A)') 'VAN DER WAALS-CORRECTION:'
       WRITE(output_unit,'(A)')
       WRITE(output_unit,'(4X,A,1X,A)') 'DCACP:', 'DISPERSION-CORRECTED ATOM-CENTRED POTENTIALS'

       IF ( any(dcacp%no_contribution(:)) ) THEN
          index_list = list_of_indices( dcacp )
          IF (.NOT. dcacp_env%include_metals) &
             WRITE(output_unit,'(11X,A)') 'METAL CENTRES WILL NOT CONTRIBUTE BY DEFAULT'
          WRITE(output_unit,'(11X,A,1X,A)') 'NO CONTRIBUTION FROM SPECIES:', &
                                             trim(adjustl(index_list))
       ENDIF
       WRITE(output_unit,'(A)')
       IF (dcacp_env%has_custom_sigma) THEN
          WRITE(output_unit,'(11X,A)') 'USING A MODIFIED SET OF DCACP PROJECTOR PARAMETERS:'
       ELSE
          WRITE(output_unit,'(11X,A)') 'DCACP PROJECTOR PARAMETERS:'
       ENDIF
       WRITE(output_unit,'(11X,A5,2X,A,7X,A)') 'Z', 'RADIUS (SIGMA_1)', 'DEPTH (SIGMA_2)'
       DO z=1,max_dcacp 
          IF ( any(ions0%iatyp(:) == z) ) THEN
             no_contribution_at_all = .TRUE.
             DO isp=1,ions1%nsp
                IF (ions0%iatyp(isp) == z) THEN
                   IF (.not. dcacp%no_contribution(isp)) no_contribution_at_all = .FALSE.
                ENDIF
             ENDDO
             IF (no_contribution_at_all) THEN
                WRITE(output_unit,'(11X,I5,2X,A)') z, 'NO CONTRIBUTION'
             ELSE
                WRITE(output_unit,'(11X,I5,2X,ES15.6,8X,ES15.6)') z, dcacp%element(z)%sigma_1, &
                                                                    dcacp%element(z)%sigma_2
                WRITE(output_unit,'(21X,A)') 'Ref: '//trim(adjustl(dcacp%element(z)%citation))
             ENDIF
          ENDIF
       ENDDO
       WRITE(output_unit,'(A)')
    ENDIF

    CONTAINS

    ! ==--------------------------------------------------------------==
    ELEMENTAL FUNCTION list_of_indices( dcacp ) &
    RESULT (indices)

       TYPE(dcacp_t), INTENT(in) :: dcacp
       INTEGER                   :: isp

       !
       ! (I3) + (,) + (space) = space for 16 indices
       !
       CHARACTER(len=80)         :: indices
       CHARACTER(len=3)          :: atom_index
       LOGICAL                   :: first_in_list

       indices       = ' ' 
       first_in_list = .true.

       DO isp=1,ions1%nsp
          atom_index = ' '
          IF (dcacp%no_contribution(isp)) THEN
             WRITE(atom_index,'(I3)') isp
             IF (first_in_list) THEN
                indices       = atom_index
                first_in_list = .false.
             ELSE
                indices = trim(adjustl(indices))//&
                          ', '//trim(adjustl(atom_index))
             ENDIF
          ENDIF
       ENDDO

    END FUNCTION list_of_indices
    ! ==--------------------------------------------------------------==
  END SUBROUTINE dcacp_report
  ! ==================================================================
  PURE SUBROUTINE get_dcacp_potential( dcacp, isp, pot_array ) 
    ! ==--------------------------------------------------------------==
    ! Adds numerical DCACP for MT-type PPs

    REAL(real_8), DIMENSION(:,:,:), &
                  INTENT(inout)     :: pot_array
    INTEGER, INTENT(in)             :: isp
    TYPE(dcacp_t), INTENT(in)       :: dcacp

    INTEGER                         :: i
    REAL(real_8)                    :: sigma_2

    sigma_2 = dcacp%element(ions0%iatyp(isp))%sigma_2

    pot_array(:,isp,dcacp%l_dcacp) = 0.0_real_8

    forall (i=1:size(pot_array,1)) pot_array(i,isp,dcacp%l_dcacp) = &
                         pot_array(i,isp,dpot_mod%lloc(isp)) + sigma_2

  END SUBROUTINE get_dcacp_potential
  ! ==================================================================
  SUBROUTINE get_dcacp_wavefunction( dcacp, isp, mesh, wf_array ) 
    ! ==--------------------------------------------------------------==
    ! Adds numerical DCACP for MT-type PPs

    REAL(real_8), DIMENSION(:), &
                  INTENT(in)        :: mesh
    REAL(real_8), DIMENSION(:,:,:), &
                  INTENT(inout)     :: wf_array
    INTEGER, INTENT(in)             :: isp
    TYPE(dcacp_t), INTENT(in)       :: dcacp

    INTEGER                         :: i
    REAL(real_8)                    :: sigma_1, prefactor


    sigma_1   = dcacp%element(ions0%iatyp(isp))%sigma_1
    prefactor = prefactor_for(sigma_1,dcacp%l_dcacp)

    wf_array(:,isp,l_dcacp) = 0.0_real_8

    FORALL (i=1:size(mesh,1))  wf_array(i,isp,dcacp%l_dcacp) = prefactor &
                                                               * wavefunction_of(mesh(i),dcacp%l_dcacp,sigma_1)

    CONTAINS

    ! ==--------------------------------------------------------------==
    ELEMENTAL FUNCTION prefactor_for( sigma_1,l_dcacp ) &
    RESULT (dcacp_prefactor)

       REAL(real_8), PARAMETER     :: sqrt_gamma = sqrt(11.6317283966_real_8)
       INTEGER, INTENT(in)         :: l_dcacp
       REAL(real_8), INTENT(in)    :: sigma_1
       REAL(real_8)                :: dcacp_prefactor

       dcacp_prefactor = sqrt(2.0_real_8) &
                         /( sqrt(sigma_1)*(sigma_1**l_dcacp)*sqrt_gamma )

    END FUNCTION prefactor_for  
    ! ==--------------------------------------------------------------==
    ELEMENTAL FUNCTION wavefunction_of( meshpoint,l_dcacp,sigma_1 ) &
    RESULT (dcacp_wf)

       INTEGER, INTENT(in)         :: l_dcacp
       REAL(real_8), INTENT(in)    :: meshpoint
       REAL(real_8), INTENT(in)    :: sigma_1
       REAL(real_8)                :: dcacp_wf

       dcacp_wf = exp( -1.0_real_8*(meshpoint*meshpoint)/(2.0_real_8*sigma_1*sigma_1) )
       dcacp_wf = dcacp_wf * meshpoint**l_dcacp
                  
    END FUNCTION wavefunction_of
    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_dcacp_wavefunction
  ! ==================================================================

  ! ==================================================================
  ! The actual library and utilities
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE construct_library( dcacp,dcacp_env )
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), PARAMETER  :: procedureN = 'construct_library'

    TYPE(dcacp_t), INTENT(inout) :: dcacp
    TYPE(dcacp_env_t), &
                   INTENT(in)    :: dcacp_env

    SELECT CASE(trim(adjustl(dcacp%xc_functional)))
      CASE('GGA_XC_BLYP')
         CALL get_blyp_library( dcacp )
      CASE('GGA_XC_BP','GGA_XC_BP86')
         CALL get_bp86_library( dcacp )
    CASE DEFAULT
      IF (.not. dcacp_env%has_custom_sigma) &
         CALL stopgm(procedureN,'No DCACP parameters available for the selected xc-functional', &
                     __LINE__,__FILE__)
    END SELECT

    CONTAINS

    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE get_blyp_library( dcacp )

       TYPE(dcacp_t), INTENT(inout) :: dcacp

       dcacp%element(  H  )%sigma_1  =  2.706646438534561_real_8
       dcacp%element(  H  )%sigma_2  = -4.0558298753945934e-04_real_8  
       dcacp%element(  H  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element( He  )%sigma_1  =  2.395256856273089_real_8
       dcacp%element( He  )%sigma_2  = -3.9158931229090230e-04_real_8 
       dcacp%element( He  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(Li:B )%sigma_1  =  0.0_real_8
       dcacp%element(Li:B )%sigma_2  =  0.0_real_8 
       dcacp%element(Li:B )%citation = 'No parameters in library'
       !
       dcacp%element(  C  )%sigma_1  =  3.107044432680790_real_8
       dcacp%element(  C  )%sigma_2  = -5.49106479060575e-04_real_8 
       dcacp%element(  C  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  N  )%sigma_1  =  2.910081224358620_real_8
       dcacp%element(  N  )%sigma_2  = -6.05377446947131e-04_real_8 
       dcacp%element(  N  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  O  )%sigma_1  =  2.57314438378083_real_8
       dcacp%element(  O  )%sigma_2  = -7.92030218954537e-04_real_8    
       dcacp%element(  O  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  F  )%sigma_1  =  2.4_real_8
       dcacp%element(  F  )%sigma_2  = -7.030e-04_real_8  
       dcacp%element(  F  )%citation = 'J. Chem. Theory Comput. 9, 955, 2013'
       !
       dcacp%element( Ne  )%sigma_1  =  2.481030426963918_real_8
       dcacp%element( Ne  )%sigma_2  = -6.4123468198010114e-04_real_8
       dcacp%element( Ne  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(Na:Si)%sigma_1  =  0.0_real_8
       dcacp%element(Na:Si)%sigma_2  =  0.0_real_8 
       dcacp%element(Na:Si)%citation = 'No parameters in library' 
       !
       dcacp%element(  P  )%sigma_1  =  3.53870284787246_real_8
       dcacp%element(  P  )%sigma_2  = -1.159905781812608e-03_real_8 
       dcacp%element(  P  )%citation = 'J. Chem. Theory Comput. 5, 2930, 2009'
       !
       dcacp%element(  S  )%sigma_1  =  3.27299135931733_real_8
       dcacp%element(  S  )%sigma_2  = -1.39292417456171e-03_real_8 
       dcacp%element(  S  )%citation = 'J. Chem. Theory Comput. 5, 23, 2008'
       !
       dcacp%element( Cl  )%sigma_1  =  2.67085120345258_real_8
       dcacp%element( Cl  )%sigma_2  = -1.599661826277423e-03_real_8 
       dcacp%element( Cl  )%citation = 'J. Chem. Theory Comput. 9, 955, 2013'
       !
       dcacp%element( Ar  )%sigma_1  =  2.774807649499496_real_8
       dcacp%element( Ar  )%sigma_2  = -1.2959773245069981e-03_real_8 
       dcacp%element( Ar  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element( K:Se)%sigma_1  =  0.0_real_8
       dcacp%element( K:Se)%sigma_2  =  0.0_real_8 
       dcacp%element( K:Se)%citation = 'No parameters in library' 
       !
       dcacp%element( Br  )%sigma_1  =  2.99969633979284467_real_8
       dcacp%element( Br  )%sigma_2  = -1.40751129235295192e-03_real_8  
       dcacp%element( Br  )%citation = 'J. Chem. Theory Comput. 9, 955, 2013'
       !
       dcacp%element( Kr  )%sigma_1  =  3.184808609280901_real_8
       dcacp%element( Kr  )%sigma_2  = -1.2946645689165006e-03_real_8 
       dcacp%element( Kr  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(Rb:Te)%sigma_1  =  0.0_real_8
       dcacp%element(Rb:Te)%sigma_2  =  0.0_real_8 
       dcacp%element(Rb:Te)%citation = 'No parameters in library' 
       !
       dcacp%element(  I  )%sigma_1  =  3.04633350043523_real_8
       dcacp%element(  I  )%sigma_2  = -2.22196291022035e-03_real_8 
       dcacp%element(  I  )%citation = 'J. Chem. Theory Comput. 9, 955, 2013'
       !
       dcacp%element(Xe:Rn)%sigma_1  =  0.0_real_8
       dcacp%element(Xe:Rn)%sigma_2  =  0.0_real_8 
       dcacp%element(Xe:Rn)%citation = 'No parameters in library' 
                          
    END SUBROUTINE get_blyp_library
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE get_bp86_library( dcacp )

       TYPE(dcacp_t), INTENT(inout) :: dcacp

       dcacp%element(  H  )%sigma_1  =  2.658039321465702_real_8
       dcacp%element(  H  )%sigma_2  = -5.5474848476271985e-04_real_8  
       dcacp%element(  H  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element( He  )%sigma_1  =  2.164544749555909_real_8
       dcacp%element( He  )%sigma_2  = -9.9130209654357913e-04_real_8
       dcacp%element( He  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(Li:B )%sigma_1  =  0.0_real_8
       dcacp%element(Li:B )%sigma_2  =  0.0_real_8
       dcacp%element(Li:B )%citation = 'No parameters in library'
       !
       dcacp%element(  C  )%sigma_1  =  3.503483143425686_real_8
       dcacp%element(  C  )%sigma_2  = -3.7065662609148662e-04_real_8
       dcacp%element(  C  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  N  )%sigma_1  =  2.818695064632713_real_8
       dcacp%element(  N  )%sigma_2  = -8.0630621824677711e-04_real_8
       dcacp%element(  N  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  O  )%sigma_1  =  2.636467588394379_real_8
       dcacp%element(  O  )%sigma_2  = -1.0649113779701964e-03_real_8
       dcacp%element(  O  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(  F  )%sigma_1  =  0.0_real_8
       dcacp%element(  F  )%sigma_2  =  0.0_real_8
       dcacp%element(  F  )%citation = 'No parameters in library'
       !
       dcacp%element( Ne  )%sigma_1  =  2.424284266287178_real_8
       dcacp%element( Ne  )%sigma_2  = -1.2919877413127942e-03_real_8
       dcacp%element( Ne  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element(Na:Si)%sigma_1  =  0.0_real_8
       dcacp%element(Na:Si)%sigma_2  =  0.0_real_8
       dcacp%element(Na:Si)%citation = 'No parameters in library'
       !
       dcacp%element(  P  )%sigma_1  =  3.73764682044944_real_8
       dcacp%element(  P  )%sigma_2  = -9.279539236382878e-04_real_8
       dcacp%element(  P  )%citation = 'J. Chem. Theory Comput. 5, 2930, 2009'
       !
       dcacp%element(  S  )%sigma_1  =  3.359413435501356_real_8
       dcacp%element(  S  )%sigma_2  = -1.3273137477856443e-03_real_8
       dcacp%element(  S  )%citation = 'J. Chem. Theory Comput. 5, 23, 2008'
       !
       dcacp%element( Cl  )%sigma_1  =  0.0_real_8
       dcacp%element( Cl  )%sigma_2  =  0.0_real_8
       dcacp%element( Cl  )%citation = 'No parameters in library'
       !
       dcacp%element( Ar  )%sigma_1  =  2.785775938005546_real_8
       dcacp%element( Ar  )%sigma_2  = -1.6422346768878398e-03_real_8
       dcacp%element( Ar  )%citation = 'Phys. Rev. B. 75, 205131, 2007'
       !
       dcacp%element( K:Br)%sigma_1  =  0.0_real_8
       dcacp%element( K:Br)%sigma_2  =  0.0_real_8
       dcacp%element( K:Br)%citation = 'No parameters in library'
       !
       dcacp%element( Kr  )%sigma_1  =  3.218226418090794_real_8
       dcacp%element( Kr  )%sigma_2  = -1.4682105527952909e-03_real_8
       dcacp%element( Kr  )%citation = 'Phys. Rev. B. 75, 205131, 2007' 
       !
       dcacp%element(Rb:Rn)%sigma_1  =  0.0_real_8
       dcacp%element(Rb:Rn)%sigma_2  =  0.0_real_8
       dcacp%element(Rb:Rn)%citation = 'No parameters in library'

    END SUBROUTINE get_bp86_library 
    ! ==--------------------------------------------------------------==
  END SUBROUTINE construct_library
  ! ==================================================================
  ELEMENTAL FUNCTION no_contribution( dcacp_env,species ) &
  RESULT (does_not_contribute)
    ! ==--------------------------------------------------------------==
    ! Check whether DCACP parameters should be requested or not
    ! Separates metals from non-metals
    ! Modifies dcacp%no_contribution for proper printing

    TYPE(dcacp_env_t), &
             INTENT(in)         :: dcacp_env
    INTEGER, INTENT(in)         :: species

    LOGICAL                     :: does_not_contribute

    does_not_contribute = .FALSE.

    IF (.NOT. dcacp_env%include_metals) THEN
        IF ( (ions0%iatyp(species) >= Li .AND. &
              ions0%iatyp(species) <= Be ) .OR.&
             (ions0%iatyp(species) >= Na .AND. &
              ions0%iatyp(species) <= Al ) .OR.&
             (ions0%iatyp(species) >= K  .AND. &
              ions0%iatyp(species) <= Ga ) .OR.&
             (ions0%iatyp(species) >= Rb .AND. &
              ions0%iatyp(species) <= Sn ) .OR.&
             (ions0%iatyp(species) >= Cs .AND. &
              ions0%iatyp(species) <= Po ) .OR.&
              ions0%iatyp(species) >= Fr ) THEN 
          does_not_contribute = .TRUE.
       ENDIF
    ENDIF

    does_not_contribute = ( does_not_contribute .OR. &
                            dcacp_env%no_contribution(species) )

    ! ==--------------------------------------------------------------==
  END FUNCTION no_contribution
  ! ==================================================================

  ! ==================================================================
  ! Broadcasting
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE bcast_dcacp_pre_init()
    ! ==--------------------------------------------------------------==
    ! Broadcasting of variables read in vdwin
 
    CALL mp_bcast_byte(dcacp,size_in_bytes_of(dcacp),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(dcacp_env,size_in_bytes_of(dcacp_env),parai%io_source,parai%cp_grp)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE bcast_dcacp_pre_init
  ! ==================================================================
  SUBROUTINE bcast_dcacp_post_init()
    ! ==--------------------------------------------------------------==
    ! Re-broadcasting of variables that have already been broadcast in
    ! ratom and that were modified
 
    CALL mp_bcast_byte(dcacp,size_in_bytes_of(dcacp),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(atom_common,size_in_bytes_of(atom_common),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(nlps_com,size_in_bytes_of(nlps_com),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nghtol,SIZE(nghtol),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nghcom,SIZE(nghcom),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rgh,SIZE(rgh),parai%io_source,parai%cp_grp)
    CALL mp_bcast(wgh,SIZE(wgh),parai%io_source,parai%cp_grp)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE bcast_dcacp_post_init
  ! ==================================================================
END MODULE dcacp_utils
