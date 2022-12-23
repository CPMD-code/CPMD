MODULE nmr_position_p_utils
  USE coor,                            ONLY: tau0
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nmr_util_p_utils,                ONLY: make_123
  USE parac,                           ONLY: parai,&
                                             paral
  USE reshaper,                        ONLY: reshape_inplace
  USE response_pmod,                   ONLY: &
       deltarvalues, llc_ind, llc_val, lower_left, lower_left_value, &
       nmr_options, nris, response1, voa_data, wanniercenters
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE util_p_utils,                    ONLY: dot_ga
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgrid
  PUBLIC :: putongrid
  PUBLIC :: set_origin
  PUBLIC :: force_xyz
  PUBLIC :: calc_lower_left_optimize
  PUBLIC :: calc_lower_left_new
  !!public :: apply_op_rx
  PUBLIC :: calc_min_dist
  PUBLIC :: calc_min_dist_alat
  PUBLIC :: calc_min_dist_alat_smooth
  PUBLIC :: print_llc
  !public :: print_dist_r_wc
  PUBLIC :: bound
  PUBLIC :: give_min_dist

CONTAINS

  ! ==================================================================
  ! THIS FILE contains all functions and routines for handling
  ! problems related to the position operator.
  ! ==================================================================



  INTEGER FUNCTION rgrid(position,direction)
    ! RETURNS the integer grid index ix/iy/iz of an arbitrary position
    ! in space (also outside the box). Convention: 1 <= rgrid <= NR#S
    ! where    # = (direction).
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    REAL(real_8) :: projector(3)

    REAL(real_8) :: position(3),projection
    INTEGER :: direction,nris,iprojection
    ! ==--------------------------------------------------------------==
    IF (direction .EQ. 1) THEN
       projector(:) = gvec_com%b1(:)
       nris=spar%nr1s
    ELSEIF (direction .EQ. 2) THEN
       projector(:) = gvec_com%b2(:)
       nris=spar%nr2s
    ELSEIF (direction .EQ. 3) THEN
       projector(:) = gvec_com%b3(:)
       nris=spar%nr3s
    ELSE
       CALL stopgm('RGRID','INVALID DIRECTION',& 
            __LINE__,__FILE__)
    ENDIF

    projection = (position(1)*projector(1) +position(2)*projector(2) +&
         position(3)*projector(3)) /parm%alat

    iprojection = NINT( (projection - NINT(projection - 0.5_real_8))*&
         REAL(nris,kind=real_8)) + 1
    IF (iprojection .GT. nris) iprojection=1

    rgrid = iprojection

    RETURN
  END FUNCTION rgrid
  ! ==================================================================
  SUBROUTINE putongrid(position)
    ! CHANGES the position such that it is placed on the nearest 
    ! real space gridpoint. No unit-cell-wrapping is done.
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: position(3)

    INTEGER                                  :: i1, i2, i3

! ==--------------------------------------------------------------==

    i1 = NINT((position(1)*gvec_com%b1(1)+position(2)*gvec_com%b1(2)+position(3)*gvec_com%b1(3))*&
         spar%nr1s/parm%alat)
    i2 = NINT((position(1)*gvec_com%b2(1)+position(2)*gvec_com%b2(2)+position(3)*gvec_com%b2(3))*&
         spar%nr2s/parm%alat)
    i3 = NINT((position(1)*gvec_com%b3(1)+position(2)*gvec_com%b3(2)+position(3)*gvec_com%b3(3))*&
         spar%nr3s/parm%alat)

    position(1) = i1*parm%a1(1)/REAL(spar%nr1s,kind=real_8)+ i2*parm%a2(1)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(1)/&
         REAL(spar%nr3s,kind=real_8)
    position(2) = i1*parm%a1(2)/REAL(spar%nr1s,kind=real_8)+ i2*parm%a2(2)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(2)/&
         REAL(spar%nr3s,kind=real_8)
    position(3) = i1*parm%a1(3)/REAL(spar%nr1s,kind=real_8) + i2*parm%a2(3)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(3)&
         /REAL(spar%nr3s,kind=real_8)

    RETURN
  END SUBROUTINE putongrid
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE set_origin(istate)
    INTEGER                                  :: istate

    INTEGER                                  :: ik

! =----------------------------------------------------------------=

    DO ik=1,3
       llc_ind(ik) = lower_left(ik,istate)
       llc_val(ik) = lower_left_value(ik,istate)
    ENDDO
    RETURN
  END SUBROUTINE set_origin
  ! ==================================================================
  SUBROUTINE force_xyz(nstate)
    ! =----------------------------------------------------------------=
    INTEGER                                  :: nstate

    INTEGER                                  :: ik, is
    REAL(real_8)                             :: avg

    DO ik=1,3
       IF (nmr_options%tforce_xyz(ik)) THEN
          avg = REAL(lower_left(ik,1),kind=real_8)
          DO is=2,nstate
             avg = avg + REAL(bound( lower_left(ik,is)-lower_left(ik,1)+&
                  nris(ik)/2,nris(ik)) - nris(ik)/2,kind=real_8) / REAL(nstate,kind=real_8)
          ENDDO
          DO is=1,nstate
             lower_left(ik,is) = bound(NINT(avg),nris(ik))
          ENDDO
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE force_xyz
  ! ==================================================================
  SUBROUTINE calc_lower_left_optimize(c0,d0,nstate)
    ! =----------------------------------------------------------------=
    ! PARALLEL
    ! EXPECTS as input the gnd-state wavefunction c0.
    ! OUTPUT  is the integer array lower_left which will
    ! \       be set to the indices of the lower left corner of the
    ! \       box which is best for the wavefunction c0(*),
    ! \       WHERE ITS OVERLAP WITH SURROUNDING ORBITALS IS
    ! \       ALSO TAKEN INTO CONSIDERATION.
    ! ***     NEW VERSION. Uses the G-space variant to calculate
    ! ***     the penalty functional:
    ! ***     Int[  |psi_i|^2   exp  -(x-x0)^2 / 2 sigma^2 ]
    ! ***     Evalutation in G-space.
    ! ***     Final formula for Penalty[x0]:
    ! ***     Sum[g: g.B2=0 AND g.B3=0 AND |g| < sqrt(20)/sigma ]
    ! ***         exp - 0.5 (2 pi sigma g.B1)^2 
    ! ***         Re{ rho(g)  exp - 2 pi i  x0  g.B1 }
    ! 
    ! ***     the sqrt(20) corresponds to an exp(-(..)^2) of 10^-10.
    ! =----------------------------------------------------------------=
    ! sigma[response_p.inc] is the spread of the damping 
    ! Gaussian in units of the 
    ! whole cell. Thus 0.1 means that the Gaussian spread is 0.1*alat.
    ! =----------------------------------------------------------------=
    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: d0(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_lower_left_optimize'

    COMPLEX(real_8), ALLOCATABLE             :: psi(:)
    COMPLEX(real_8), POINTER                 :: rho(:)
    INTEGER                                  :: ggj, ggk, ggl, ierr, ig, ig1, &
                                                ig_indexarray(3,100), ii, ij, &
                                                ik, il, ir, is, js, ks, &
                                                max_jkl
    LOGICAL                                  :: found_new_state
    LOGICAL, ALLOCATABLE                     :: done(:), tag(:)
    REAL(real_8)                             :: alpha, min_penalty
    REAL(real_8), ALLOCATABLE                :: density_overlap(:,:), &
                                                penalty(:,:), psi2(:)

! Variables
! This array contains the numbers IG of the G-vectors G which
! satisfy  G.B2 = 0  and  G.B3 = 0  (for ig_indexarray(1,*))
! resp.    G.B1 = 0  and  G.B3 = 0  (for ig_indexarray(2,*)) etc.
! OVERLAP quantities:
! =----------------------------------------------------------------=

    ALLOCATE(density_overlap(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(density_overlap)!,nstate*nstate)
    ALLOCATE(tag(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    tag(:)=.FALSE.
    !call iazzero(tag,nstate)
    ALLOCATE(done(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    done(:)=.FALSE.
    !call iazzero(done,nstate)
    ALLOCATE(penalty(3,MAX(spar%nr1s,spar%nr2s,spar%nr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   penalty)!, 3*MAX(spar%nr1s,spar%nr2s,spar%nr3s))

    max_jkl=NINT(2._real_8/deltarvalues%sigma)
    IF (max_jkl .GE. 100)  CALL stopgm('NMR/LLCOPT',&
         'INTERNAL ERROR, sigma too small (<0.02).',& 
         __LINE__,__FILE__)


    ! The ig - Index array is now built up. A G-vector (ig) is "good"
    ! if two of its reciprocal Miller indices, G.B1, G.B2 and G.B3 are zero.
    ! Then, the corresponding IG-value is stored in the array
    ! "ig_indexarray".
    ij=0
    il=0
    ik=0
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%ngw
       ggj = NINT(dot_ga(ig,1))
       ggk = NINT(dot_ga(ig,2))
       ggl = NINT(dot_ga(ig,3))
       IF (ggk.EQ.0 .AND. ggl.EQ.0 .AND. ggj .LE. 6.8_real_8/deltarvalues%sigma) THEN
          ij = ij + 1
          ig_indexarray(1,ij) = ig
       ELSEIF (ggj.EQ.0 .AND. ggl.EQ.0 .AND. ggk .LE. 6.8_real_8/deltarvalues%sigma)&
            THEN
          ik = ik + 1
          ig_indexarray(2,ik) = ig
       ELSEIF (ggj.EQ.0 .AND. ggk.EQ.0 .AND. ggl .LE. 6.8_real_8/deltarvalues%sigma)&
            THEN
          il = il + 1
          ig_indexarray(3,il) = ig
       ENDIF
       IF (ij.GT.max_jkl .OR. ik.GT.max_jkl .OR. il.GT.max_jkl) THEN
          IF (paral%io_parent)&
               WRITE (6,'(A,I3,A,I6,4(A,I3))') 'PROC ',parai%mepos,': IG=',ig,&
               '; jkl=(',ggj,',',ggk,',',ggl,'), max=',max_jkl
          CALL stopgm('NMR/LLC',&
               'INTERNAL ERROR. Contact the STUPIDO PROGRAMMER',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    ! ij, ik, and il now contain the number of g-vectors found.
    ! =----------------------------------------------------------------=
    ! The orbital density needs to be computed. The orbital densities are
    ! stored in G space in the D0 array and the overlap between
    ! these densities is done:
    ALLOCATE(psi(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi2(2*maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(psi)!,nnr1)
    CALL zeroing(psi2)!,2*nnr1)
    DO is=1,nstate
       CALL ffttor(c0(1,is),psi2,psi,ncpw%ngw,.FALSE.)
       DO ir=1,fpar%nnr1
          psi2(ir) = psi2(ir)*psi2(ir)
       ENDDO
       CALL ffttog(psi2,d0(1,is),psi,ncpw%ngw,.TRUE.)
    ENDDO
    ! create density overlap matrix, stored in density_overlap:
    DO is=2,nstate
       DO js=1,is-1
          density_overlap(is,js) = dotp(ncpw%ngw,d0(:,is),d0(:,js))
          density_overlap(js,is) = density_overlap(is,js)
       ENDDO
    ENDDO
    CALL mp_sum(density_overlap,nstate*nstate,parai%allgrp)
    DO is=1,nstate
       DO js=1,nstate
          IF (paral%io_parent.AND.density_overlap(is,js).LT.-1._real_8) THEN
             WRITE(6,'(A,I3,A,2I5,A,F20.10)')&
                  'PROC ',parai%me,': STATES ',is,js,&
                  ' ==> ',density_overlap(is,js)
             CALL stopgm('LLC/OPT','negative density overlap',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDDO
    ! =----------------------------------------------------------------=
    ! Now, all those densities are summed together where the overlap is
    ! larger than the predefined value OVERLAP_THREASHOLD
    DO is=1,nstate
       done(is) = .FALSE.    ! initialization
    ENDDO
    ALLOCATE(rho(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO is=1,nstate
       IF (.NOT. done(is)) THEN
          DO js=1,nstate
             tag(js) = .FALSE.! sub-initialization
          ENDDO
          tag(is) = .TRUE.
          found_new_state = .TRUE.
          DO WHILE (found_new_state)
             found_new_state = .FALSE.
             DO js=1,nstate
                IF (tag(js)) THEN
                   DO ks=1,nstate
                      IF (density_overlap(js,ks).GT.deltarvalues%overlap_threashold)&
                           THEN
                         IF (.NOT. tag(ks)) THEN
                            tag(ks) = .TRUE.
                            found_new_state = .TRUE.
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO             ! "find new states" -- loop.
          DO js=1,nstate
             IF (tag(js)) done(js) = .TRUE.
          ENDDO
          ! Now, in "tag()", we have all those states which we want to consider
          ! together. We have to add up their densities:
          CALL zeroing(rho)!,ngw)
          DO js=1,nstate
             IF (tag(js)) THEN
                CALL daxpy(2*ncpw%ngw,1._real_8,d0(1,js),1,rho,1)
                IF (paral%io_parent)&
                     WRITE (6,'(A,I6,A,I6)') 'State group ',is,&
                     ': Joining state ',js
             ENDIF
          ENDDO


          ! =----------------------------------------------------------------=
          ! Next, the penalty functional is calculated
          ! for the accumulated orbital densities
          ! \  foreach point on the three axis
          ! \     using each one of the previously determined G-vectors.
          ! =----------------------------------------------------------------=
          ! Calculate the actual penalty values in the "first" direction.
          ! Formula: sum_{g: k=0, l=0} exp(-0.5 (2pi j sigma)^2)
          ! \                               exp(2pi i  ix/Nx j)
          CALL zeroing(penalty)!,3*MAX(spar%nr1s,spar%nr2s,spar%nr3s))
          ! do ir=NRXPL(MEPOS,1),NRXPL(MEPOS,2)
          DO ir=1,spar%nr1s          ! NB: NOT from NRXPL(1) to (2)!!!
             alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
             IF (geq0) penalty(1,ir) = 0.5_real_8*REAL(rho(1))
             DO ii=1,ij    ! sum over all the g vectors
                ! in the "first" direction
                ig=ig_indexarray(1,ii)
                penalty(1,ir) = penalty(1,ir) +&
                     EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,1))**2) *&
                     REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,1)*alpha,kind=real_8)) * rho(ig))
             ENDDO         ! the g-vectors
          ENDDO             ! the points ir=1..NR1S

          ! The "second" direction: Same formula except that g.A2 != 0
          ! and g.A1 = g.A3 = 0
          DO ir=1,spar%nr2s
             alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
             IF (geq0) penalty(2,ir) = 0.5_real_8*REAL(rho(1))
             DO ii=1,ik    ! sum over all the g vectors
                ! in the "first" direction
                ig=ig_indexarray(2,ii)
                penalty(2,ir) = penalty(2,ir) +&
                     EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,2))**2) *&
                     REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,2)*alpha,kind=real_8)) * rho(ig))
             ENDDO         ! the g-vectors
          ENDDO             ! the points ir=1..NR2S

          ! The "third" direction: Same formula except that g.A3 != 0
          ! and g.A1 = g.A2 = 0
          DO ir=1,spar%nr3s
             alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)
             IF (geq0) penalty(3,ir) = 0.5_real_8*REAL(rho(1))
             DO ii=1,il    ! sum over all the g vectors
                ! in the "first" direction
                ig=ig_indexarray(3,ii)
                penalty(3,ir) = penalty(3,ir) +&
                     EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,3))**2) *&
                     REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,3)*alpha,kind=real_8)) * rho(ig))
             ENDDO         ! the g-vectors
          ENDDO             ! the points ir=1..NR3S


          CALL mp_sum(penalty,3*MAX(spar%nr1s,spar%nr2s,spar%nr3s),parai%allgrp)


          ! The minimum is searched, and its coordinates are stored
          ! in the lower_left - array:
          min_penalty=1.e3_real_8
          DO ir=1,spar%nr1s
             IF (penalty(1,ir) .LT. min_penalty) THEN
                min_penalty = penalty(1,ir)
                lower_left(1,is) = ir
             ENDIF
          ENDDO
          min_penalty=1.e3_real_8
          DO ir=1,spar%nr2s
             IF (penalty(2,ir) .LT. min_penalty) THEN
                min_penalty = penalty(2,ir)
                lower_left(2,is) = ir
             ENDIF
          ENDDO
          min_penalty=1.e3_real_8
          DO ir=1,spar%nr3s
             IF (penalty(3,ir) .LT. min_penalty) THEN
                min_penalty = penalty(3,ir)
                lower_left(3,is) = ir
             ENDIF
          ENDDO


          DO js=1,nstate
             IF (tag(js)) THEN
                lower_left(1,js) = lower_left(1,is)
                lower_left(2,js) = lower_left(2,is)
                lower_left(3,js) = lower_left(3,is)

                lower_left_value(1,js) = - spar%nr1s/2
                lower_left_value(2,js) = - spar%nr2s/2
                lower_left_value(3,js) = - spar%nr3s/2
             ENDIF
          ENDDO
       ENDIF                 ! if not done.
    ENDDO                     ! is

    DEALLOCATE(rho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(density_overlap,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tag,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(done,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(penalty,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE calc_lower_left_optimize
  ! ==================================================================
  SUBROUTINE calc_lower_left_new(c0,nstate)
    ! =----------------------------------------------------------------=
    ! PARALLEL
    ! EXPECTS as input the gnd-state wavefunction c0.
    ! OUTPUT  is the integer array lower_left which will
    ! \       be set to the indices of the lower left corner of the
    ! \       box which is best for the wavefunction c0(*).
    ! ***     NEW VERSION. Uses the G-space variant to calculate
    ! ***     the penalty functional:
    ! ***     Int[  |psi_i|^2   exp  -(x-x0)^2 / 2 sigma^2 ]
    ! ***     Evalutation in G-space.
    ! ***     Final formula for Penalty[x0]:
    ! ***     Sum[g: g.B2=0 AND g.B3=0 AND |g| < sqrt(20)/sigma ]
    ! ***         exp - 0.5 (2 pi sigma g.B1)^2 
    ! ***         Re{ rho(g)  exp - 2 pi i  x0  g.B1 }
    ! 
    ! ***     the sqrt(20) corresponds to an exp(-(..)^2) of 10^-10.
    ! =----------------------------------------------------------------=
    ! sigma[response_p.inc] is the spread of the damping 
    ! Gaussian in units of the 
    ! whole cell. Thus 0.1 means that the Gaussian spread is 0.1*alat.
    ! =----------------------------------------------------------------=
    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_lower_left_new'

    COMPLEX(real_8), ALLOCATABLE             :: psi(:)
    COMPLEX(real_8), POINTER                 :: rho(:)
    INTEGER                                  :: ggj, ggk, ggl, ierr, ig, ig1, &
                                                ig_indexarray(3,100), ii, ij, &
                                                ik, il, ir, is, max_jkl
    REAL(real_8)                             :: alpha, min_penalty
    REAL(real_8), ALLOCATABLE                :: penalty(:,:)
    REAL(real_8), ALLOCATABLE, TARGET        :: psi2(:)

! Variables
! This array contains the numbers IG of the G-vectors G which
! satisfy  G.B2 = 0  and  G.B3 = 0  (for ig_indexarray(1,*))
! resp.    G.B1 = 0  and  G.B3 = 0  (for ig_indexarray(2,*)) etc.

    max_jkl=NINT(2._real_8/deltarvalues%sigma)
    IF (max_jkl .GE. 100)  CALL stopgm('NMR/LLC',&
         'INTERNAL ERROR, sigma too small (<0.02).',& 
         __LINE__,__FILE__)

    ALLOCATE(penalty(3,MAX(spar%nr1s,spar%nr2s,spar%nr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! calculate the damping exponentials
    ! do ij=1,max_jkl
    ! damping(ij) = exp( -  0.5 * (2*3.141592654_real_8*sigma*real(ij,kind=real_8))**2)
    ! enddo


    ! The ig - Index array is now built up. A G-vector (ig) is "good"
    ! if two of its reciprocal Miller indices, G.B1, G.B2 and G.B3 are zero.
    ! Then, the corresponding IG-value is stored in the array
    ! "ig_indexarray".
    ij=0
    il=0
    ik=0
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%ngw
       ggj = NINT(dot_ga(ig,1))
       ggk = NINT(dot_ga(ig,2))
       ggl = NINT(dot_ga(ig,3))
       IF (ggk.EQ.0 .AND. ggl.EQ.0 .AND. ggj .LE. 6.8_real_8/deltarvalues%sigma) THEN
          ij = ij + 1
          ig_indexarray(1,ij) = ig
       ELSEIF (ggj.EQ.0 .AND. ggl.EQ.0 .AND. ggk .LE. 6.8_real_8/deltarvalues%sigma)&
            THEN
          ik = ik + 1
          ig_indexarray(2,ik) = ig
       ELSEIF (ggj.EQ.0 .AND. ggk.EQ.0 .AND. ggl .LE. 6.8_real_8/deltarvalues%sigma)&
            THEN
          il = il + 1
          ig_indexarray(3,il) = ig
       ENDIF
       IF (ij.GT.max_jkl .OR. ik.GT.max_jkl .OR. il.GT.max_jkl) THEN
          IF (paral%io_parent)&
               WRITE (6,'(A,I3,A,I6,4(A,I3))') 'PROC ',parai%mepos,': IG=',ig,&
               '; jkl=(',ggj,',',ggk,',',ggl,'), max=',max_jkl
          CALL stopgm('NMR/LLC',&
               'INTERNAL ERROR. Contact the STRONZO PROGRAMMER',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    ! ij, ik, and il now contain the number of g-vectors found.



    ! Here, the penalty functional is calculated
    ! foreach state
    ! \  foreach point on the three axis
    ! \     using each one of the previously determined G-vectors.
    ALLOCATE(psi(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi2(2*maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(psi2, (/maxfft/), rho)
    DO is=1,nstate

       ! The orbital density needs to be computed:
       CALL zeroing(psi)!,nnr1)
       CALL zeroing(psi2)!,2*nnr1)
       CALL ffttor(c0(1,is),psi2,psi,ncpw%ngw,.FALSE.)
       DO ir=1,fpar%nnr1
          psi2(ir) = psi2(ir)*psi2(ir)
       ENDDO
       CALL ffttog(psi2,psi2,psi,ncpw%nhg,.TRUE.)


       ! Calculate the actual penalty values in the "first" direction.
       ! Formula: sum_{g: k=0, l=0} exp(-0.5 (2pi j sigma)^2)
       ! \                               exp(2pi i  ix/Nx j)
       CALL zeroing(penalty)!,3*MAX(spar%nr1s,spar%nr2s,spar%nr3s))
       ! do ir=NRXPL(MEPOS,1),NRXPL(MEPOS,2)
       DO ir=1,spar%nr1s            ! NB: NOT from NRXPL(1) to (2)!!!
          alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
          IF (geq0) penalty(1,ir) = 0.5_real_8*REAL(rho(1))! the G=0 contribution.
          DO ii=1,ij        ! sum over all the g vectors
             ! in the "first" direction
             ig=ig_indexarray(1,ii)
             penalty(1,ir) = penalty(1,ir)+ EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,&
                  1))**2) * REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,1)*alpha,kind=real_8)) *&
                  rho(ig))
          ENDDO             ! the g-vectors
       ENDDO                 ! the points ir=1..NR1S

       ! The "second" direction: Same formula except that g.A2 != 0
       ! and g.A1 = g.A3 = 0
       DO ir=1,spar%nr2s
          alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
          IF (geq0) penalty(2,ir) = 0.5_real_8*REAL(rho(1))
          DO ii=1,ik
             ig=ig_indexarray(2,ii)
             penalty(2,ir) = penalty(2,ir)+ EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,&
                  2))**2) * REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,2)*alpha,kind=real_8)) *&
                  rho(ig))
          ENDDO             ! the g-vectors
       ENDDO                 ! the points ir=1..NR2S

       ! The "third" direction: Same formula except that g.A3 != 0
       ! and g.A1 = g.A2 = 0
       DO ir=1,spar%nr3s
          alpha = REAL(ir-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)
          IF (geq0) penalty(3,ir) = 0.5_real_8*REAL(rho(1))
          DO ii=1,il
             ig=ig_indexarray(3,ii)
             penalty(3,ir) = penalty(3,ir)+ EXP(-0.5_real_8*(deltarvalues%sigma*dot_ga(ig,&
                  3))**2) * REAL(EXP(CMPLX(0._real_8,+dot_ga(ig,3)*alpha,kind=real_8)) *&
                  rho(ig))
          ENDDO             ! the g-vectors
       ENDDO                 ! the points ir=1..NR3S


       CALL mp_sum(penalty,3*MAX(spar%nr1s,spar%nr2s,spar%nr3s),parai%allgrp)


       ! The minimum is searched, and its coordinates are stored
       ! in the lower_left - array:
       IF (paral%parent) THEN
          min_penalty=1._real_8
          DO ir=1,spar%nr1s
             IF (penalty(1,ir) .LT. min_penalty) THEN
                min_penalty = penalty(1,ir)
                lower_left(1,is) = ir
             ENDIF
          ENDDO
          min_penalty=1._real_8
          DO ir=1,spar%nr2s
             IF (penalty(2,ir) .LT. min_penalty) THEN
                min_penalty = penalty(2,ir)
                lower_left(2,is) = ir
             ENDIF
          ENDDO
          min_penalty=1._real_8
          DO ir=1,spar%nr3s
             IF (penalty(3,ir) .LT. min_penalty) THEN
                min_penalty = penalty(3,ir)
                lower_left(3,is) = ir
             ENDIF
          ENDDO
       ENDIF                 ! parent.


       lower_left_value(1,is) = - spar%nr1s/2
       lower_left_value(2,is) = - spar%nr2s/2
       lower_left_value(3,is) = - spar%nr3s/2

    ENDDO                     ! the states is=1,nstate
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL mp_bcast(lower_left,SIZE(lower_left),parai%source,parai%allgrp)

    DEALLOCATE(penalty,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE calc_lower_left_new
  ! ==================================================================
  INTEGER FUNCTION bound(i,imax)
    ! returns i mod imax
    IMPLICIT NONE
    INTEGER :: i,imax
    IF (i .LT. 1) THEN
       bound = i + imax*(1+INT(-i/imax))
    ELSEIF (i .GT. imax) THEN
       bound = i - imax*INT((i-1)/imax)
    ELSE
       bound = i
    ENDIF
    RETURN
  END FUNCTION bound
  ! ==================================================================
  REAL(real_8) FUNCTION trivial_dist(pos1,pos2)
    IMPLICIT NONE
    REAL(real_8) :: pos1(3),pos2(3),d
    d = (pos1(1)-pos2(1)) **2 + (pos1(2)-pos2(2)) **2 + (pos1(3)-pos2(&
         3)) **2
    trivial_dist = SQRT(d)
    RETURN
  END FUNCTION trivial_dist
  ! ==================================================================
  REAL(real_8) FUNCTION pbc_dist(pos1,pos2)
    ! Calculates and returns the distance (in a.u.) of two given
    ! positions in space adopting the minimum image convention.
    ! The positions must be given as a 3D-position vector in a.u.
    IMPLICIT NONE
    REAL(real_8) :: pos1(3),pos2(3),pos3(3),d,min_dist_
    INTEGER :: n1,n2,n3,ik

    min_dist_ = 100000._real_8
    DO n1=-2,2
       DO n2=-2,2
          DO n3=-2,2
             DO ik=1,3
                pos3(ik) = pos2(ik) + n1*parm%a1(ik) + n2*parm%a2(ik) + n3*parm%a3(ik)
                d        = trivial_dist(pos1,pos3)
                IF ( d .LT. min_dist_ ) min_dist_ = d
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    pbc_dist = min_dist_
    RETURN
  END FUNCTION pbc_dist
  ! ==================================================================
  INTEGER FUNCTION min_dist(pos1,pos2,n_pbc)
    ! computes pos1 - pos2 in the minimum image convention, using
    ! the periodicity defined through N_PBC.
    ! The actual calculation is:
    ! min_dist(r1,r2,N) = min [ r1 - r2,  r1 - r2 + N, r1 - r2 - N ]
    ! The minimum procedure is understood as the minimum absolute value,
    ! whereas the value returned from this function will contain the 
    ! correct sign.
    IMPLICIT NONE
    INTEGER :: pos1,pos2,n_pbc,r12(5),i

    DO i=-2,2
       r12(i+3) = pos1-pos2 + i*n_pbc
       ! d12(i+3) = abs(r12(i+3))
    ENDDO

    min_dist = r12(1)
    DO i=2,5
       IF (ABS(r12(i)) .LT. ABS(min_dist)) THEN
          ! identical to the previous routine would be ".le" instead of .lt.
          min_dist = r12(i)
       ENDIF
    ENDDO


    RETURN
  END FUNCTION min_dist
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_min_dist(is1,is2,dij)
    ! Returns in dij the minimum image convention vector that links the 
    ! two virtual cells of states is1 and is2: dij = r_is2 - r_is1,
    ! i.e. the connection vector FROM is1 TO is2.
    ! NB: This routine assumes that all lower_left_value being identical!!
    ! VARIANT: COMPLETE DIRECT SPACE DISTANCE CALCULATION
    ! Arguments:
    INTEGER                                  :: is1, is2
    REAL(real_8)                             :: dij(3)

    INTEGER                                  :: ir, n1, n2, n3
    REAL(real_8)                             :: dij0(3), dij_test(3), Ndij, &
                                                Ndij0, Ndij_test

! Variables:

    Ndij = 0._real_8
    Ndij0=0.0_real_8
    DO ir=1,3
       dij0(ir) =&
            (lower_left(1,is2)-lower_left(1,is1)) *parm%a1(ir)/spar%nr1s +&
            (lower_left(2,is2)-lower_left(2,is1)) *parm%a2(ir)/spar%nr2s +&
            (lower_left(3,is2)-lower_left(3,is1)) *parm%a3(ir)/spar%nr3s
       Ndij0 = ndij0 + dij0(ir)*dij0(ir)
       IF (.NOT.&
            (lower_left_value(ir,is1).EQ.lower_left_value(ir,is2)))&
            THEN
          IF (paral%io_parent)&
               WRITE (6,*) is1,is2,ir,&
               lower_left_value(ir,is1),lower_left_value(ir,is2)
          CALL stopgm('CALC_MIN_DIST','LLC-values are not equal! !',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    Ndij = Ndij0
    CALL dcopy(3,dij0,1,dij,1)

    DO n1=-1,1
       DO n2=-1,1
          DO n3=-1,1
             Ndij_test = 0._real_8
             DO ir=1,3
                dij_test(ir) = dij0(ir)&
                     + n1*parm%a1(ir) + n2*parm%a2(ir) + n3*parm%a3(ir)
                Ndij_test = ndij_test + dij_test(ir)*dij_test(ir)
             ENDDO
             IF (Ndij_test .LT. Ndij) THEN
                CALL dcopy(3,dij_test,1,dij,1)
                Ndij = Ndij_test
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE calc_min_dist
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_min_dist_alat(is1,is2,dij)
    ! Returns in dij the minimum image convention vector that links the 
    ! two virtual cells of states is1 and is2: dij = r_is2 - r_is1,
    ! i.e. the connection vector FROM is1 TO is2.
    ! NB: This routine assumes that all lower_left_value being identical!!
    ! VARIANT: DISTANCE CALCULATION in units of the lattice vectors.
    ! only identical to the previous routine if the cell is 
    ! cubic/orthorhombic!
    ! Arguments:
    INTEGER                                  :: is1, is2
    REAL(real_8)                             :: dij(3)

    INTEGER                                  :: i, ijk(3)
    REAL(real_8)                             :: ibyN(3)

! Variables:

    DO i=1,3
       ijk(i)  = min_dist(lower_left(i,is2),lower_left(i,is1),nris(i))
       ibyN(i) = ijk(i) / REAL(nris(i),kind=real_8)
    ENDDO

    ! Direction A1:
    DO i=1,3
       dij(i)  =  ibyN(1)*parm%a1(i) + ibyn(2)*parm%a2(i) + ibyn(3)*parm%a3(i)
    ENDDO


    RETURN
  END SUBROUTINE calc_min_dist_alat
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_min_dist_alat_smooth(is1,is2,dij)
    ! Returns in dij the minimum image convention vector that links the 
    ! two virtual cells of states is1 and is2: dij = r_is2 - r_is1,
    ! i.e. the connection vector FROM is1 TO is2.
    ! NB: This routine assumes that all lower_left_value being identical!!
    ! VARIANT: DISTANCE CALCULATION in units of the lattice vectors.
    ! only identical to the previous routine if the cell is 
    ! cubic/orthorhombic!
    ! ADDITIONAL SMOOTHING for the cases where the distance between the
    ! states is roughly half a lattice vector, at least for one of 
    ! the lattice directions.
    ! Arguments:
    INTEGER                                  :: is1, is2
    REAL(real_8)                             :: dij(3)

    REAL(real_8), PARAMETER                  :: spread = 0.1_real_8 

    INTEGER                                  :: i, ijk(3)
    REAL(real_8)                             :: arg, ibyN(3), smoothing(3)

    DO i=1,3
       ijk(i)  = min_dist(lower_left(i,is2),lower_left(i,is1),nris(i))
       ibyN(i) = REAL(ijk(i),kind=real_8) / REAL(nris(i),kind=real_8)
       ! min_dist returns integers (from -N/2 to +N/2)
       ! the smoothing is then based on these integers.
       arg     =  ABS(ibyN(i)) - 0.5_real_8
       smoothing(i) = 1._real_8-EXP( -0.5_real_8* (arg/spread)**2 )
    ENDDO

    ! Direction A1:
    DO i=1,3
       dij(i)  =  smoothing(1)*ibyN(1)*parm%a1(i) + smoothing(2)*ibyn(2)*&
            parm%a2(i)+ smoothing(3)*ibyN(3)*parm%a3(i)
    ENDDO


    RETURN
  END SUBROUTINE calc_min_dist_alat_smooth
  ! ==================================================================
  FUNCTION give_min_dist(is1,is2,iB) RESULT(res)
    ! Returns the minimum image convention DISTANCE (in a.u.) that links the 
    ! two virtual cells of states is1 and is2: dij = r_is2 - r_is1,
    ! i.e. the connection vector FROM is1 TO is2.
    ! The states are FIRST PROJECTED onto the "iB" plane, i.e. the distance
    ! is actually calculated as  || O(iB) x  ( r_is2 - r_is1 ) ||
    ! NB: This routine assumes that all lower_left_value being identical!!
    ! Arguments:
    INTEGER                                  :: is1, is2, iB
    REAL(real_8)                             :: res

    INTEGER                                  :: iBB, iBBB
    REAL(real_8)                             :: dij(3)

! Variables:

    CALL calc_min_dist(is1,is2,dij)
    CALL make_123(iB,iBB,iBBB)
    res = SQRT(dij(iBB)*dij(ibb) + dij(iBBB)*dij(ibbb))
    RETURN
  END FUNCTION give_min_dist
  ! ==================================================================
  SUBROUTINE print_llc(nstate)
    INTEGER                                  :: nstate

    INTEGER                                  :: ik, is
    REAL(real_8)                             :: r_c(3), r_l(3)

    IF (.NOT. paral%parent) RETURN

    IF (paral%io_parent)&
         WRITE(6,89)
    IF (paral%io_parent)&
         WRITE(6,90)
    DO is=1,nstate
       DO ik=1,3
          r_l(ik) = 0._real_8
          r_c(ik) = 0._real_8
          r_l(ik) = r_l(ik) + (lower_left(1,is)-1)*parm%a1(ik)/REAL(nris(1),kind=real_8)
          r_c(ik) = r_c(ik) + (bound(lower_left(1,is)-lower_left_value(&
               1,is),nris(1))-1)* parm%a1(ik)/REAL(nris(1),kind=real_8)
          r_l(ik) = r_l(ik) + (lower_left(2,is)-1)*parm%a2(ik)/REAL(nris(2),kind=real_8)
          r_c(ik) = r_c(ik) + (bound(lower_left(2,is)-lower_left_value(&
               2,is),nris(2))-1)* parm%a2(ik)/REAL(nris(2),kind=real_8)
          r_l(ik) = r_l(ik) + (lower_left(3,is)-1)*parm%a3(ik)/REAL(nris(3),kind=real_8)
          r_c(ik) = r_c(ik) + (bound(lower_left(3,is)-lower_left_value(&
               3,is),nris(3))-1)* parm%a3(ik)/REAL(nris(3),kind=real_8)
          IF(response1%tvoa) voa_data%R_WC(ik,is) = r_c(ik)
       ENDDO

       IF (paral%io_parent)&
            WRITE(6,81) is,&
            (r_l(ik),ik=1,3),&
            (r_c(ik),ik=1,3)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,91)

    RETURN


81  FORMAT ("* [",i4,"]: ",3f7.2,5x,3f7.2,7x,"*")
89  FORMAT ("** VIRTUAL CELLS ",30("*"),"*NMR/EPR*RESPONSE*")
90  FORMAT ("** STATE ** LLC_x  LLC_y  LLC_z  ",&
         "***  CTR_x  CTR_y  CTR_z *******")
91  FORMAT (47("*"),"*NMR/EPR*RESPONSE*")
  END SUBROUTINE print_llc
  ! ==================================================================
  SUBROUTINE print_dist_r_wc(nstate)
    INTEGER                                  :: nstate

    INTEGER                                  :: iat, iatsp, ik, is, isp
    REAL(real_8)                             :: r_c(3)

70  FORMAT ("* DIST |   ATOM/TOTAL   STATE   DISTANCE [au]")
71  FORMAT ("* DIST |   ",i4,"/",i4,4x,i5,3x,f12.6)

    IF (.NOT. paral%parent) RETURN

    IF (paral%io_parent)&
         WRITE (6,*) 'FIRST METHOD: LLC (FROM LOCALIZE)'
    DO is=1,nstate
       IF (paral%io_parent)&
            WRITE (6,70)
       iatsp = 0
       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             iatsp=iatsp+1
             DO ik=1,3
                r_c(ik) = 0._real_8
                r_c(ik) = r_c(ik) + (bound(lower_left(1,is)-&
                     lower_left_value(1,is),nris(1))-1)* parm%a1(ik)/REAL(nris(1),kind=real_8)
                r_c(ik) = r_c(ik) + (bound(lower_left(2,is)-&
                     lower_left_value(2,is),nris(2))-1)* parm%a2(ik)/REAL(nris(2),kind=real_8)
                r_c(ik) = r_c(ik) + (bound(lower_left(3,is)-&
                     lower_left_value(3,is),nris(3))-1)* parm%a3(ik)/REAL(nris(3),kind=real_8)
             ENDDO
             IF (paral%io_parent)&
                  WRITE (6,71) iat,iatsp,is,pbc_dist(tau0(1,iat,isp),r_c)
          ENDDO
       ENDDO
    ENDDO

    IF (paral%io_parent)&
         WRITE (6,*) 'SECOND METHOD: WANNIERCENTERS (FROM LOCALIZE)'
    DO is=1,nstate
       IF (paral%io_parent)&
            WRITE (6,70)
       iatsp = 0
       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             iatsp=iatsp+1
             IF (paral%io_parent)&
                  WRITE (6,71) iat,iatsp,is,&
                  pbc_dist(tau0(1,iat,isp),wanniercenters(1,is))
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE print_dist_r_wc

END MODULE nmr_position_p_utils


! ==================================================================
SUBROUTINE apply_op_rx(inp,beta,outp,alpha,coord)
  ! ==================================================================
  ! PARALLEL
  ! The most tricky routine in the whole NMR calculation.
  ! Read with care. ;-)
  ! =----------------------------------------------------------------=
  ! CALCULATES:
  ! outp = alpha * outp  +  beta * r_coord * inp
  ! the two arrays may be identical (!)
  ! 
  ! applies an (1 - exp-(x^2)/l^2) where x is the distance to the
  ! nearest LLC point. l is nr_s/20.
  ! =----------------------------------------------------------------=
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:fpar,parm,spar,parap
  USE parac, ONLY : paral,parai
  USE response_pmod , ONLY:llc_ind,llc_val,nmr_options
  USE nmr_position_p_utils, ONLY : bound
  IMPLICIT NONE
  REAL(real_8) :: inp(fpar%kr1,fpar%kr2,fpar%kr3), beta, &
      outp(fpar%kr1,fpar%kr2,fpar%kr3), alpha
  INTEGER                                    :: coord

  CHARACTER(*), PARAMETER                    :: procedureN = 'apply_op_rx'
  INTEGER, PARAMETER                         :: divisor = 20

  INTEGER                                    :: i1, i2, i3, ierr, &
                                                ind_ll1_rel, isub
  INTEGER, SAVE                              :: maxnris = 0
  REAL(real_8)                               :: delta2, fsmooth1, fsmooth2, &
                                                fsmooth3, in, lambda2(3), &
                                                Ndelta2(3), out, Rop, Rop0
  REAL(real_8), ALLOCATABLE, SAVE            :: smoothing(:,:)
  REAL(real_8), SAVE                         :: aN1, aN2, aN3

! SMOOTHING: indices go from 0 to maxnris, 0 means that we are
! at the LLC (smoothing-factor has value 0).
! =----------------------------------------------------------------=

  CALL tiset('     r_psi',isub)

  IF (coord .EQ. 0 .AND. maxnris .NE. 0) THEN
     DEALLOCATE(smoothing,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     maxnris = 0
     RETURN
  ENDIF

  IF (maxnris .EQ. 0) THEN  ! smoothing not yet allocated...
     maxnris = MAX(spar%nr1s,spar%nr2s,spar%nr3s)
     ALLOCATE(smoothing(3,0:maxnris+1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     lambda2(1) = REAL(spar%nr1s*spar%nr1s,kind=real_8)/REAL(divisor*divisor,kind=real_8)
     lambda2(2) = REAL(spar%nr2s*spar%nr2s,kind=real_8)/REAL(divisor*divisor,kind=real_8)
     lambda2(3) = REAL(spar%nr3s*spar%nr3s,kind=real_8)/REAL(divisor*divisor,kind=real_8)
     IF (nmr_options%tsmoothing) THEN
        DO i1=0,maxnris+1
           delta2            = REAL(i1*i1,kind=real_8)
           Ndelta2(1)        = REAL((i1-spar%nr1s) * (i1-spar%nr1s),kind=real_8)
           Ndelta2(2)        = REAL((i1-spar%nr2s) * (i1-spar%nr2s),kind=real_8)
           Ndelta2(3)        = REAL((i1-spar%nr3s) * (i1-spar%nr3s),kind=real_8)
           smoothing(1,i1) = (1 - EXP(- delta2   /lambda2(1)))*&
                (1 - EXP(-Ndelta2(1)/lambda2(1)))
           smoothing(2,i1) = (1 - EXP(- delta2   /lambda2(2)))*&
                (1 - EXP(-Ndelta2(2)/lambda2(2)))
           smoothing(3,i1) = (1 - EXP(- delta2   /lambda2(3)))*&
                (1 - EXP(-Ndelta2(3)/lambda2(3)))
        ENDDO
     ELSE                  ! no smoothing:
        DO i1=0,maxnris+1
           smoothing(1,i1) = 1._real_8
           smoothing(2,i1) = 1._real_8
           smoothing(3,i1) = 1._real_8
        ENDDO
     ENDIF
  ENDIF


  IF (coord .LT. 1 .OR. coord .GT. 3) THEN
     IF (paral%io_parent)&
          WRITE (6,*) '*** apply_op_rx: WRONG COORDINATE. Nothing done.'
     RETURN
  ENDIF


  aN1   = parm%a1(coord) / REAL(spar%nr1s,kind=real_8)
  aN2   = parm%a2(coord) / REAL(spar%nr2s,kind=real_8)
  aN3   = parm%a3(coord) / REAL(spar%nr3s,kind=real_8)
  ind_ll1_rel    =  llc_ind(1) - (parap%nrxpl(parai%mepos,1)-1)
  Rop0  = (llc_val(3) +&
       1 + spar%nr3s-llc_ind(3)) * aN3+&
       (llc_val(2) +&
       1 + spar%nr2s-llc_ind(2)) * aN2+&
       (llc_val(1)&
       + bound(2 + spar%nr1s-ind_ll1_rel,spar%nr1s) -1) * aN1


  ! NB: After one X-loop in parallel, the Rop-variable does NOT
  ! have the value it had before the loop: We add nr1 times a1/nr1S, but
  ! we substract a1. The result is NOT zero. (!!)

  IF (ind_ll1_rel .LE. 1) THEN
     DO i3          =  1,parm%nr3
        fsmooth3     =  beta     * smoothing(3,ABS(i3-llc_ind(3)))
        IF (i3 .EQ. llc_ind(3))Rop0     =  rop0 - parm%a3(coord)
        DO i2        =  1,llc_ind(2)-1
           fsmooth2   =  fsmooth3 * smoothing(2,llc_ind(2)-i2)
           Rop        =  Rop0
           DO i1      =  1,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,i1-ind_ll1_rel)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, first loop
        Rop0         =  rop0 - parm%a2(coord)
        DO i2        =  llc_ind(2),parm%nr2
           fsmooth2   =  fsmooth3 * smoothing(2,i2-llc_ind(2))
           Rop        =  Rop0
           DO i1      =  1,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,i1-ind_ll1_rel)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, second loop
        Rop0         =  rop0 + aN3
     ENDDO
  ELSEIF (ind_ll1_rel .GT. parm%nr1) THEN
     ! THEN the (relative) LLC is 
     ! (right) outside of MY cell.
     DO i3          =  1,parm%nr3
        fsmooth3     =  beta     * smoothing(3,ABS(i3-llc_ind(3)))
        IF (i3 .EQ. llc_ind(3))Rop0     =  rop0 - parm%a3(coord)
        DO i2        =  1,llc_ind(2)-1
           fsmooth2   =  fsmooth3 * smoothing(2,llc_ind(2)-i2)
           Rop        =  Rop0
           DO i1      =  1,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,ind_ll1_rel-i1)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, first loop
        Rop0         =  rop0 - parm%a2(coord)
        DO i2        =  llc_ind(2),parm%nr2
           fsmooth2   =  fsmooth3 * smoothing(2,i2-llc_ind(2))
           Rop        =  Rop0
           DO i1      =  1,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,ind_ll1_rel-i1)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, second loop
        Rop0         =  rop0 + aN3
     ENDDO
  ELSE                      ! THEN, 1 < ind_ll1_rel < nr1+1,
     ! i.e. the llc is within MY cell
     DO i3          =  1,parm%nr3
        fsmooth3     =  beta     * smoothing(3,ABS(i3-llc_ind(3)))
        IF (i3 .EQ. llc_ind(3))Rop0     =  rop0 - parm%a3(coord)
        DO i2        =  1,llc_ind(2)-1
           fsmooth2   =  fsmooth3 * smoothing(2,llc_ind(2)-i2)
           Rop        =  Rop0
           DO i1      =  1,ind_ll1_rel-1
              fsmooth1 =  fsmooth2 * smoothing(1,ind_ll1_rel-i1)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1, first loop
           Rop        =  rop - parm%a1(coord)! here is THE wrap.
           DO i1      =  ind_ll1_rel,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,i1-ind_ll1_rel)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1, second loop
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, first loop
        Rop0         =  rop0 - parm%a2(coord)
        DO i2        =  llc_ind(2),parm%nr2
           fsmooth2   =  fsmooth3 * smoothing(2,i2-llc_ind(2))
           Rop        =  Rop0
           DO i1      =  1,ind_ll1_rel-1
              fsmooth1 =  fsmooth2 * smoothing(1,ind_ll1_rel-i1)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1, first loop
           Rop        =  rop - parm%a1(coord)! here is THE wrap.
           DO i1      =  ind_ll1_rel,parm%nr1
              fsmooth1 =  fsmooth2 * smoothing(1,i1-ind_ll1_rel)
              out      =  alpha * outp(i1,i2,i3)
              in       =  fsmooth1 * inp(i1,i2,i3)
              outp(i1,i2,i3) = out + Rop * in
              Rop      =  rop + aN1
           ENDDO         ! i1, second loop
           Rop0       =  rop0 + aN2
        ENDDO             ! i2, second loop
        Rop0         =  rop0 + aN3
     ENDDO
  ENDIF

  CALL tihalt('     r_psi',isub)
  RETURN
END SUBROUTINE apply_op_rx
