MODULE gndstate_p_utils
  USE cppt,                            ONLY: gk
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE gvec,                            ONLY: gvec_com
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE nmr_position_p_utils,            ONLY: bound
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: llc_ind,&
                                             llc_val,&
                                             nmr_options
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: second_moment
  !public :: generate_gaussian_in_g
  PUBLIC :: gauss_gen_read
  !!public :: save_my_wannier_gaussian

CONTAINS

  ! ==================================================================
  SUBROUTINE second_moment(c0,coc,matrix,psi)
    ! NOW PARALLEL!
    ! ==================================================================
    ! CALCULATES the matrix of second moments for ONE STATE C0:
    ! matrix(i,j) = Int[  phi(r) *  r_i r_j * phi(r) ]
    ! \                -  Int[  phi(r) *  r_i * phi(r) ]
    ! \                  *Int[  phi(r) *  r_j * phi(r) ]
    ! \           = second_moment_(i,j)
    ! 
    ! Here, r_i means (r-r0)_i, where r-r0 is the position operator
    ! relative to the center of charge r0. This r0 is returned in COC.
    ! =----------------------------------------------------------------=
    ! *** INPUT: ***
    ! C0:       one orbital state in G space
    ! psi,scr:  scratch arrays
    ! lscr:     dimension (in real(8) :: ) of scr
    ! *** OUTPUT: ***
    ! coc:      center of charge, 3-vector which will hold 
    ! \         the direct space coordinates (rx,ry,rz) in
    ! \         atomic units of the Center Of Charge of the orbital.
    ! matrix:   real(8) :: 3x3 matrix (symmetric) containing 
    ! \         Int[  phi(r) *  r_i r_j * phi(r) ].
    ! coc and matrix can both directly be re-inserted into the
    ! routine generate_gaussian_g below.
    ! =----------------------------------------------------------------=
    ! =----------------------------------------------------------------=
    ! Arguments:
    COMPLEX(real_8)                          :: c0(ncpw%ngw)
    REAL(real_8) :: coc(3), matrix(3,3), psi(fpar%kr1,fpar%kr2s,fpar%kr3s)

    CHARACTER(*), PARAMETER                  :: procedureN = 'second_moment'

    INTEGER                                  :: i1, i2, ierr, imax(3), ir, &
                                                isub, ix, iy, iz, PE_bcast
    LOGICAL                                  :: debug
    REAL(real_8)                             :: delta_omega, global_rhoM, &
                                                PE_bcast_d, projection, &
                                                r_exp(3), rhoM
    REAL(real_8), ALLOCATABLE                :: scr(:,:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

! =----------------------------------------------------------------=
! variables:
! =----------------------------------------------------------------=

    CALL tiset('     r_psi',isub)
    ALLOCATE(scr(fpar%kr1,fpar%kr2s,fpar%kr3s,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    debug = .TRUE.
    IF ((debug).AND.paral%io_parent)&
         WRITE(6,*)'calculating second moments...'

    delta_omega = 1._real_8/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8) ! wfn normalization for
    ! real space integration


    CALL ffttor(c0,psi,scr,ncpw%ngw,.TRUE.)
    !$omp parallel do private(ir)
    DO ir=1,fpar%nnr1
       scr(ir,1,1,1) = psi(ir,1,1)*psi(ir,1,1)
    ENDDO

    rhoM = 0._real_8
    DO ix=1,parm%nr1
       DO iy=1,parm%nr2
          DO iz=1,parm%nr3
             IF (scr(ix,iy,iz,1) .GT. rhoM) THEN
                rhoM = scr(ix,iy,iz,1)
                imax(1)=ix + parap%nrxpl(parai%mepos,1)
                imax(2)=iy
                imax(3)=iz
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF ((debug).AND.paral%io_parent)&
         WRITE(6,*)'PE ',parai%mepos,&
         'MAX point found at ',imax(1),imax(2),imax(3)

    ! find the DIRECT SPACE LOCATION of the center of
    ! charge of the orbital:
    ! this does not make sense here, r_exp is not initialized
    ! and coc will be overwritten in the next statement anyway.
    ! call DCOPY(3,r_exp,1,coc,1)
    coc(1) =&
         (imax(1)-1)*parm%a1(1)/REAL(spar%nr1s,kind=real_8)&
         + (imax(2)-1)*parm%a2(1)/REAL(spar%nr2s,kind=real_8)&
         + (imax(3)-1)*parm%a3(1)/REAL(spar%nr3s,kind=real_8)
    coc(2) =&
         (imax(1)-1)*parm%a1(2)/REAL(spar%nr1s,kind=real_8)&
         + (imax(2)-1)*parm%a2(2)/REAL(spar%nr2s,kind=real_8)&
         + (imax(3)-1)*parm%a3(2)/REAL(spar%nr3s,kind=real_8)
    coc(3) =&
         (imax(1)-1)*parm%a1(3)/REAL(spar%nr1s,kind=real_8)&
         + (imax(2)-1)*parm%a2(3)/REAL(spar%nr2s,kind=real_8)&
         + (imax(3)-1)*parm%a3(3)/REAL(spar%nr3s,kind=real_8)

    IF ((debug).AND.paral%io_parent)&
         WRITE(6,'(A,3F10.4)')&
         'corresponding to point: ',coc


    global_rhoM = rhoM        ! now, determine the maximum
    ! over all processors:
    CALL mp_max(global_rhoM,parai%allgrp)
    IF ( ABS(global_rhoM-rhoM) .LT. 1.e-10_real_8 ) THEN
       PE_bcast_d = REAL(parai%mepos,kind=real_8)
    ELSE
       PE_bcast_d = 0._real_8
    ENDIF
    CALL mp_max(PE_bcast_d,parai%allgrp) ! makes the broadcasting PE unique (!)
    PE_bcast = NINT(PE_bcast_d)
    CALL mp_bcast(imax,SIZE(imax),PE_bcast,parai%allgrp)
    ! dirty but simple.



    ! Initialization necessary for the routine that 
    ! applies the position operator:
    nmr_options%tsmoothing = .FALSE.
    llc_val(1) = - spar%nr1s/2
    llc_val(2) = - spar%nr2s/2
    llc_val(3) = - spar%nr3s/2
    llc_ind(1) = bound(imax(1) + spar%nr1s/2 ,spar%nr1s)
    llc_ind(2) = bound(imax(2) + spar%nr2s/2 ,spar%nr2s)
    llc_ind(3) = bound(imax(3) + spar%nr3s/2 ,spar%nr3s)



    DO i1=1,3
       CALL apply_op_rx(psi,1._real_8,scr(1,1,1,1),0._real_8,i1)
       r_exp(i1) = ddot(fpar%nnr1,psi,1,scr(1,1,1,1),1)*delta_omega
       IF ((debug .AND. paral%parent).AND.paral%io_parent)&
            WRITE(6,'(A,I2,A,F14.8)') 'r expectation ',i1,' = ',r_exp(i1)
    ENDDO


    ! find the DIRECT SPACE LOCATION of the center of
    ! charge of the orbital:
    CALL daxpy(3,1._real_8,r_exp,1,coc,1)
    ! coc(1) = coc(1) + 
    ! &    (imax(1)-1)*a1(1)/real(nr1s,kind=real_8)
    ! &    + (imax(2)-1)*a2(1)/real(nr2s,kind=real_8)
    ! &    + (imax(3)-1)*a3(1)/real(nr3s,kind=real_8)
    ! coc(2) = coc(2) + 
    ! &    (imax(1)-1)*a1(2)/real(nr1s,kind=real_8)
    ! &    + (imax(2)-1)*a2(2)/real(nr2s,kind=real_8)
    ! &    + (imax(3)-1)*a3(2)/real(nr3s,kind=real_8)
    ! coc(3) = coc(3) + 
    ! &    (imax(1)-1)*a1(3)/real(nr1s,kind=real_8)
    ! &    + (imax(2)-1)*a2(3)/real(nr2s,kind=real_8)
    ! &    + (imax(3)-1)*a3(3)/real(nr3s,kind=real_8)

    IF ((debug).AND.paral%io_parent)&
         WRITE(6,'(A,3F10.4)')&
         'NEW CENTER OF CHARGE at ',coc



    ! find the REAL SPACE MESH indexes which define the center of
    ! charge of the orbital, i.e. the desired center for our gaussian:
    ! N_1  =  <r_exp> . b1    and so on.

    ! FIRST DIMENSION (nr1s)
    projection = r_exp(1) * gvec_com%b1(1)  +   r_exp(2) * gvec_com%b1(2)&
         + r_exp(3) * gvec_com%b1(3)
    projection = projection / parm%alat * REAL(spar%nr1s,kind=real_8)
    imax(1) = imax(1) + NINT(projection)
    llc_ind(1) = bound(imax(1) - spar%nr1s/2,spar%nr1s)
    ! SECOND DIMENSION (nr2s)
    projection = r_exp(1) * gvec_com%b2(1)  +   r_exp(2) * gvec_com%b2(2)&
         + r_exp(3) * gvec_com%b2(3)
    projection = projection / parm%alat * REAL(spar%nr2s,kind=real_8)
    imax(2) = imax(2) + NINT(projection)
    llc_ind(2) = bound(imax(2) - spar%nr2s/2,spar%nr2s)
    ! THIRD DIMENSION (nr3s)
    projection = r_exp(1) * gvec_com%b3(1)  +   r_exp(2) * gvec_com%b3(2)&
         + r_exp(3) * gvec_com%b3(3)
    projection = projection / parm%alat * REAL(spar%nr3s,kind=real_8)
    imax(3) = imax(3) + NINT(projection)
    llc_ind(3) = bound(imax(3) - spar%nr3s/2,spar%nr3s)

    IF ((debug).AND.paral%io_parent)&
         WRITE(6,'(A,3I4)')&
         'NEW CENTER POINT at ',imax



    ! Start calculating the second moments...
    DO i1=1,3
       CALL apply_op_rx(psi,1._real_8,scr(1,1,1,1),0._real_8,i1)
       r_exp(i1) = ddot(fpar%nnr1,psi,1,scr(1,1,1,1),1)*delta_omega
       IF ((debug .AND. paral%parent).AND.paral%io_parent)&
            WRITE(6,'(A,I2,A,F14.8)') 'r expectation ',i1,' = ',r_exp(i1)
       DO i2=i1,3
          CALL apply_op_rx(scr(1,1,1,1),1._real_8,scr(1,1,1,2),0._real_8,i2)
          matrix(i1,i2) = ddot(fpar%nnr1,psi,1,scr(1,1,1,2),1)*delta_omega
          matrix(i2,i1) = matrix(i1,i2)
       ENDDO                 ! i2
    ENDDO                     ! i1

    CALL mp_sum(matrix,3*3,parai%allgrp)
    CALL mp_sum(r_exp,3,parai%allgrp)

    DO i1=1,3
       DO i2=1,3
          matrix(i1,i2) = matrix(i1,i2) - r_exp(i1)*r_exp(i2)
       ENDDO
    ENDDO

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('     r_psi',isub)
    RETURN
  END SUBROUTINE second_moment
  ! ==================================================================




  ! ==================================================================
  SUBROUTINE generate_gaussian_in_g(c0,coc,m)
    ! PARALLEL!
    ! ==================================================================
    ! GENERATES a test gaussian in c0 (returned in g space) with
    ! CENTER coc  (direct space position in atomic units)
    ! SPREAD defined by MATRIX (a spread matrix as computed 
    ! with the routine before)

    COMPLEX(real_8)                          :: c0(ncpw%ngw)
    REAL(real_8)                             :: coc(3), m(3,3)

    INTEGER                                  :: ig
    LOGICAL                                  :: debug
    REAL(real_8)                             :: arg, g1, g2, g3, gg, scaling

    debug = .TRUE.
    IF ((debug).AND.paral%io_parent)&
         WRITE(6,*)'constructing gaussian orbital...'

    CALL zeroing(c0)!,ngw)
    !$omp parallel do private(ig,g1,g2,g3,arg,gg)
    DO ig=1,ncpw%ngw
       g1  = gk(1,ig)*parm%tpiba
       g2  = gk(2,ig)*parm%tpiba
       g3  = gk(3,ig)*parm%tpiba

       arg = g1*coc(1) + g2*coc(2) + g3*coc(3)

       ! gg  = g1*g1*m(1,1) + g2*g2*m(2,2) + g3*g3*m(3,3)
       ! &          + 2*(g1*g2*m(1,2) + g1*g3*m(1,3) + g2*g3*m(2,3))

       ! DMB: physically we only need one time xy, xz, and yz
       ! (otherwise the gaussian is twice as large in these directions)
       gg  = g1*g1*m(1,1) + g2*g2*m(2,2) + g3*g3*m(3,3)&
            + (g1*g2*m(1,2) + g1*g3*m(1,3) + g2*g3*m(2,3))

       ! c0(ig) = CMPLX(cos(arg),-sin(arg))*exp(-gg)
       ! DMB: changed factor...
       c0(ig) = CMPLX(COS(arg),-SIN(arg),kind=real_8)*EXP(-0.5_real_8*gg)
    ENDDO

    scaling = SQRT(dotp(ncpw%ngw,c0,c0))
    CALL dscal(2*ncpw%ngw,1._real_8/scaling,c0,1)

    RETURN
  END SUBROUTINE generate_gaussian_in_g
  ! ==================================================================


  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! This subroutine generates the wannier functions from        C
  ! a gaussian fit saved in the file "gauss-wan.dat",           C
  ! and build up the corresponding states for each molecules.   C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE gauss_gen_read(nstate,c0,nmol,STATEmin,STATEmax)

    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: nmol, STATEmin(nmol), &
                                                STATEmax(nmol)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gauss_gen_read'

    COMPLEX(real_8), ALLOCATABLE             :: Ct(:,:)
    INTEGER                                  :: Gindx, i, ierr, istate, j, &
                                                kstate, mol, r_nstate, &
                                                states_to_read
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: etot, r_r0(3), &
                                                R_sigma_mat(3,3), Sxx, Sxy, &
                                                Sxz, Syy, Syz, Szz

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! opening the file
    IF (paral%io_parent)&
         CALL fileopen(500,'gauss-wan.dat',fo_old,ferror)
    ! Did the file open OK?
    IF (ferror) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Error:the file gauss-wan.dat is not there! '
       CALL stopgm('GAUSS_GEN_READ','FILE NOT FOUND',& 
            __LINE__,__FILE__)
    ENDIF

    IF (paral%io_parent)&
         WRITE(6,*) 'Reading data from file...'

    ! Reading index data from file..
    IF (paral%io_parent)&
         READ(500,*) etot
    IF (paral%io_parent)&
         READ(500,*) r_nstate
    IF (paral%io_parent)&
         WRITE(6,*) 'The file contains data for ',r_nstate,' states'
    IF (paral%io_parent)&
         WRITE(6,*) 'Total energy of the saved state Etot = ',etot

    ! Calculating the number of states to read
    states_to_read=nstate/nmol

    IF (paral%io_parent)&
         WRITE(6,*) 'Number of states per molecule ',states_to_read

    ! checking the number of states
    IF (states_to_read.NE.r_nstate) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Need more states for each molecules...exiting'
       RETURN
    ENDIF

    ! Everything is OK...allocating the Ct array (2*dims for CPLX*16)

    ALLOCATE(Ct(ncpw%ngw,r_nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! reading the wannier functions
    IF (paral%io_parent)&
         WRITE(6,*) 'Proceeding... loading the data for ',nmol,&
         ' molecule(s)'

    ! Main loop on the read states
    DO istate=1,r_nstate
       ! CHECK IF SIGMA or SIGMA2 IS SAVED!!!
       ! Reading gaussian parameters, Rx,ry,rz and S matrix
       IF (paral%io_parent)&
            READ(500,*) r_r0(1),r_r0(2),r_r0(3),Sxx,Sxy,Sxz,&
            Syy,Syz,Szz
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,3(1X,F8.5))') 'POS:',r_r0(1),r_r0(2),r_r0(3)
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,6(1X,F8.5))') 'SIG:',Sxx,Sxy,Sxz,Syy,Syz,Szz

       ! Generating the symetric S matrix
       ! (scaled to fit to Wannier parameters,but S is already sigma^2...)
       R_sigma_mat(1,1) = Sxx
       R_sigma_mat(1,2) = Sxy
       R_sigma_mat(2,1) = Sxy
       R_sigma_mat(1,3) = Sxz
       R_sigma_mat(3,1) = Sxz
       R_sigma_mat(2,2) = Syy
       R_sigma_mat(2,3) = Syz
       R_sigma_mat(3,2) = Syz
       R_sigma_mat(3,3) = Szz

       IF (paral%io_parent)&
            WRITE(6,*) 'SIGMA matrix'
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(6,fmt='(3(1X,F8.5))')&
               (R_sigma_mat(i,j),j=1,3)
       ENDDO
       ! Fully anisotropic generation... S is a 3x3 symetric matrix of sigma^2
       CALL generate_gaussian_in_g(Ct(1,istate),r_r0,R_sigma_mat)

       ! Building the Fourrier representation of the Gaussian
       ! using Cg = [E^(.5 G^2 Sig^2) E^(-i G R0)]
       ! DO Gindx=1,NGW
       ! First: calculating the product E^(-i G.r0)
       ! calculating G.r0
       ! GbyR0=GK(1,Gindx)*R_R0(1)+GK(2,Gindx)*R_R0(2)
       !            +GK(3,Gindx)*R_R0(3)
       ! Adjusting the units
       ! GbyR0=-GbyR0*TPIBA
       ! Calculating E^i(G.r) using cos(x)+i*sin(x)
       ! EXPiGR0=CMPLX(cos(GbyR0),sin(GbyR0))

       ! Second: calculating E^(.5 G^2 Sigma^2)
       ! Isotropic case
       ! EGSIGMA = EXP(-0.5_real_8*HG(Gindx)*TPIBA2*R_sigma(1)**2)
       ! Anisotropic case Sx,Sy,Sz
       ! Gbysigma=(GK(1,Gindx)*R_sigma(1))**2 
       !              +(GK(2,Gindx)*R_sigma(2))**2
       !              +(GK(3,Gindx)*R_sigma(3))**2
       ! Adjusting the units
       ! Gbysigma=Gbysigma*TPIBA2
       ! EGSIGMA = EXP(-0.5_real_8*Gbysigma)

       ! Now calculating the G representation 
       ! Ct(Gindx,istate)= EGSIGMA*EXPiGR0
       ! ENDDO
    ENDDO

    ! Closing the file...
    IF (paral%io_parent)&
         CALL fileclose(500)
    IF (paral%io_parent)&
         WRITE(6,*) 'Building',nstate,' Wannier states'
    ! Now filling up the C0 array with the required ammount of data
    ! initialising the global state counter
    kstate=1
    ! Loop on number of molecules
    DO mol=1,nmol
       STATEmin(mol)=kstate
       ! Loop on number of read states
       DO istate=1,r_nstate
          ! Loop on Gvecs
          !$omp parallel do private(Gindx)
          DO Gindx=1,ncpw%ngw
             c0(Gindx,kstate)=Ct(gindx,istate)
          ENDDO
          ! updating the global state counter
          kstate=kstate+1
       ENDDO
       STATEmax(mol)=kstate-1
    ENDDO

    DO mol=1,nmol
       IF (paral%io_parent)&
            WRITE(6,*) 'molecule: ',mol,' STATEmin = ',&
            STATEmin(mol),' STATEmax = ',STATEmax(mol)
    ENDDO

    IF (paral%io_parent)&
         WRITE(6,*) 'Done! '
    ! Freeing memory
    DEALLOCATE(Ct,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

  END SUBROUTINE gauss_gen_read


END MODULE gndstate_p_utils


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! This subroutine saves the wannier functions in the file     C
! "gauss-wan.dat" in gaussian format.                         C 
! Upon completion the file will contain:                      C
! ETOT                                                        C
! NSTATE                                                      C
! Rx Ry Rz Sxx Sxy Sxz Syy Syz Szz                            C
! -------- -----------------------                            C
! position     spread matrix                                  C
! C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE save_my_wannier_gaussian(nstate,c0,etot,psi)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:ncpw
  USE parac, ONLY : paral,parai
  USE fileopenmod , ONLY:fo_def
  USE gndstate_p_utils, ONLY : second_moment
  USE fileopen_utils, ONLY : fileclose
  USE fileopen_utils, ONLY : fileopen
  USE fft_maxfft,                      ONLY: maxfftn
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate)
  REAL(real_8)                               :: etot, psi(maxfftn)

  INTEGER                                    :: i, istate, j
  LOGICAL                                    :: ferror
  REAL(real_8)                               :: r0(3), s_matrix(3,3), Sxx, &
                                                Sxy, Sxz, Syy, Syz, Szz

! Including the "essential" variables...
! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
! of G vectors, HG(NGW) is the array containing the square of the norm of
! each G vector. (according to Daniel...)

  CALL stopgm('ROTATE_..._P','not parallel.',& 
       __LINE__,__FILE__)

  ! opening the file
  IF (paral%io_parent)&
       CALL fileopen(500,'gauss-wan.dat',fo_def,ferror)
  IF (paral%io_parent)&
       WRITE(6,*) 'Saving ',nstate,' states...'
  IF (paral%io_parent)&
       WRITE(6,*) 'Total energy for this set Etot = ',etot

  ! Saving data to file..
  IF (paral%io_parent)&
       WRITE(500,*) etot
  IF (paral%io_parent)&
       WRITE(500,*) nstate

  ! Main loop on the states
  DO istate=1,nstate
     ! Determining the centre+spread matrix
     CALL second_moment(c0(1,istate),r0,s_matrix(1,1),psi)

     ! printing out the data...
     IF (paral%io_parent)&
          WRITE(6,*) '-----------------------------------------------'
     IF (paral%io_parent)&
          WRITE(6,*) 'STATE: ',istate
     IF (paral%io_parent)&
          WRITE(6,fmt='(1X,A,3(1X,F8.5))') 'POS:',r0(1),r0(2),r0(3)
     IF (paral%io_parent)&
          WRITE(6,*) 'SIGMA matrix'
     DO i=1,3
        IF (paral%io_parent)&
             WRITE(6,fmt='(3(1X,F8.5))')&
             (S_matrix(i,j),j=1,3)
     ENDDO

     ! Picking the spread parameters 
     ! (equivalent to sigma^2 with respect to Gaussian parameters)
     Sxx = s_matrix(1,1)
     Syy = s_matrix(2,2)
     Szz = s_matrix(3,3)

     Sxy = s_matrix(1,2)
     Sxz = s_matrix(1,3)

     Syz = s_matrix(2,3)

     IF (paral%io_parent)&
          WRITE(6,fmt='(1X,A,6(F8.5))') 'SIG:',Sxx,Sxy,Sxz,Syy,Syz,Szz
     IF (paral%io_parent)&
          WRITE(6,*) '-----------------------------------------------'

     ! Saving the data
     IF (paral%io_parent)&
          WRITE(500,*) r0(1),r0(2),r0(3),Sxx,Sxy,Sxz,Syy,Syz,Szz

  ENDDO

  ! Closing the file...
  IF (paral%io_parent)&
       CALL fileclose(500)

END SUBROUTINE save_my_wannier_gaussian
