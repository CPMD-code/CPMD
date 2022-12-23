MODULE nmr_full_p_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE nmr_position_p_utils,            ONLY: give_min_dist
  USE nmr_util_p_utils,                ONLY: make_123,&
                                             printtime
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: firstshelldistance,&
                                             lower_left,&
                                             nris
  USE system,                          ONLY: parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: make_listofstates
  PUBLIC :: optimize_llc
  PUBLIC :: optimize_llc_indep

CONTAINS

  ! ==================================================================
  ! THIS FILE contains all routines that are related to the "full"
  ! calculation of the NMR chemical shift tensor.
  ! ==================================================================
  SUBROUTINE make_listofstates(neighborlist,dowfopt,nstate)
    ! This routine looks at the lowerleft positions for all states (and
    ! all B-field directions) and makes a neighborlist such that 
    ! states that follow each other in this list have a very close 
    ! lower-left-corner. Note that this proximity search is done 
    ! _projecting_ the distance on the plane orthogonal to the B direction.

    INTEGER                                  :: nstate
    LOGICAL                                  :: dowfopt(3,nstate)
    INTEGER                                  :: neighborlist(3,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'make_listofstates'

    INTEGER                                  :: iB, ierr, initial_state, ir, &
                                                is, next_initial_state, &
                                                present_index
    LOGICAL, ALLOCATABLE                     :: already_done(:)
    REAL(real_8)                             :: distance, min_distance

    ALLOCATE(already_done(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    DO iB=1,3
       initial_state=1
       present_index=2
       neighborlist(iB,1) = 1! INITIAL STATE.
       already_done(1) = .TRUE.
       dowfopt(iB,1)     = .TRUE.
       DO is=2,nstate
          already_done(is) = .FALSE.
          dowfopt(iB,is)     = .TRUE.
       ENDDO

       DO WHILE (present_index .LE. nstate)
          min_distance = parm%alat! dummy: any large value.

          ! Search for states which are within a range of firstshelldistance:
          DO is=1,nstate
             IF (.NOT. already_done(is)) THEN
                distance = give_min_dist(is,initial_state,iB)
                IF (distance .LE. firstshelldistance) THEN
                   neighborlist(iB,present_index) = is
                   present_index = present_index +1
                   already_done(is) = .TRUE.
                   dowfopt(iB,is) = .FALSE.

                   ! Now: search for the NEXT state:
                   ! Take that one whose distance to the present one
                   ! is closest to 2*firstshelldistance.
                ELSEIF (ABS(distance-2*firstshelldistance).LT.&
                     min_distance) THEN
                   min_distance = distance
                   next_initial_state = is
                ENDIF! distance vs. firstshelldistance
             ENDIF   ! ...if this state not already done
          ENDDO         ! istate; find all states within 2.5A.
          IF (present_index .LE. nstate) THEN
             ! ...then there are more states to do
             initial_state = next_initial_state
             neighborlist(iB,present_index) = initial_state
             present_index = present_index +1
             already_done(initial_state) = .TRUE.
             dowfopt(iB,initial_state) = .TRUE.
          ENDIF
       ENDDO               ! do while not finished with all states
    ENDDO                     ! iB


    ! CHECK IT:
    DO iB=1,3
       DO is=1,nstate
          already_done(is) = .FALSE.
       ENDDO
       DO is=1,nstate
          already_done(neighborlist(iB,is)) = .TRUE.
          IF (paral%io_parent)&
               WRITE(6,'(A,I1,A,I4,A,3I4,L1)') 'field in ',iB,&
               ',state ',neighborlist(iB,is),'  --> ',(lower_left(ir,&
               neighborlist(iB,is)),ir=1,3),dowfopt(ib,&
               neighborlist(iB,is))
       ENDDO

       DO is=1,nstate
          IF (.NOT. already_done(is)) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'ERROR ERRORIS: state ',is,' not in list! !!'
          ENDIF
       ENDDO
    ENDDO
    ! checking completed.

    DEALLOCATE(already_done,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE make_listofstates
  ! ==================================================================

  ! ==================================================================
  ! ==================================================================
  SUBROUTINE optimize_llc(neighborlist,dowfopt,nstate)
    INTEGER                                  :: nstate
    LOGICAL                                  :: dowfopt(3,nstate)
    INTEGER                                  :: neighborlist(3,nstate)

    INTEGER :: averagedeviation(3), count, deviation, iB, ibb, ibbb, iprev, &
      ir, is, isi, newllc(3), presentllc(3)

    IF (paral%io_parent) THEN
       DO is=1,nstate
          WRITE (6,'("State ",I3,": LLC before = ",3I5)')is,&
               (lower_left(ir,is),ir=1,3)
       ENDDO
    ENDIF
    DO iB=1,3
       CALL make_123(iB,ibb,ibbb)
       DO isi=1,nstate
          is = neighborlist(iB,isi)
          IF (dowfopt(iB,is)) THEN
             ! THEN... start a new averaging process.
             DO ir=1,3
                presentllc(ir) = lower_left(ir,is)
                averagedeviation(ir) = 0
             ENDDO
             count = 1
          ELSE
             count = count + 1
             DO ir=1,3
                deviation = lower_left(ir,is)-presentllc(ir)
                deviation = MOD(deviation + (NRiS(ir)*3)/2,nris(ir))-&
                     NRiS(ir)/2
                averagedeviation(ir) = averagedeviation(ir)+ deviation
                newllc(ir) = presentllc(ir) +&
                     NINT(averagedeviation(ir)/REAL(count,kind=real_8))
                newllc(ir) = MOD( newllc(ir) + NRiS(ir),nris(ir))
             ENDDO
             ! Set the averages for the previous states:
             DO iprev = 1, count
                is = neighborlist(iB,isi-iprev+1)
                lower_left(ibb,is)  = newllc(ibb)
                lower_left(ibbb,is) = newllc(ibbb)
             ENDDO
          ENDIF
       ENDDO


       IF (paral%parent) THEN
          DO isi=1,nstate
             is = neighborlist(iB,isi)
             IF (paral%io_parent)&
                  WRITE (6,&
                  '("State ",I3,": LLC after (iB=",I1,") : ",3I5)')&
                  is,iB,(lower_left(ir,is),ir=1,3)
          ENDDO
       ENDIF
    ENDDO


    RETURN
  END SUBROUTINE optimize_llc

  ! ==================================================================
  ! ==================================================================
  SUBROUTINE optimize_llc_indep(neighborlist,dowfopt,nstate)
    INTEGER                                  :: nstate
    LOGICAL                                  :: dowfopt(3,nstate)
    INTEGER                                  :: neighborlist(3,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'optimize_llc_indep'

    INTEGER                                  :: count, dist, distsum, &
                                                disttonext, ierr, inext, &
                                                inow, ir, is, iss, istot, &
                                                maxdist(3), newllc
    INTEGER, ALLOCATABLE                     :: old_llc(:,:)
    LOGICAL, ALLOCATABLE                     :: done(:)
    REAL(real_8)                             :: aa(3)

! integer presentllc(3),deviation,averagedeviation(3),count

    ALLOCATE(old_llc(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(done(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    aa(1) = parm%a1(1)*parm%a1(1) + parm%a1(2)*parm%a1(2) + parm%a1(3)*parm%a1(3)
    aa(2) = parm%a2(1)*parm%a2(1) + parm%a2(2)*parm%a2(2) + parm%a2(3)*parm%a2(3)
    aa(3) = parm%a3(1)*parm%a3(1) + parm%a3(2)*parm%a3(2) + parm%a3(3)*parm%a3(3)

    DO ir=1,3
       aa(ir) = SQRT(aa(ir))
       maxdist(ir) =  NINT(firstshelldistance / aa(ir) * NRiS(ir))
       IF (paral%io_parent)&
            WRITE(6,*)&
            ' ***************  maxdist points in ',ir,': ',maxdist(ir)
       DO is=1,nstate
          old_llc(ir,is) = lower_left(ir,is)
          neighborlist(ir,is) = 0
       ENDDO
    ENDDO



    DO ir=1,3
       DO is=1,nstate
          done(is) = .FALSE.
       ENDDO
       inext = 1
       istot = 1
       DO WHILE (istot .LE. nstate)
          inow = inext
          done(inow) = .TRUE.
          neighborlist(ir,istot) = inow
          istot = istot + 1
          count = 1
          disttonext = NRiS(ir)! any large dummy value
          distsum = 0
          ! if (parent) WRITE(6,*)'*** NEW START: ir=',ir,'; is=',inow,
          ! &           '; llc=',lower_left(ir,inow)
          DO is=1,nstate
             IF (.NOT. done(is)) THEN
                dist = lower_left(ir,is) - lower_left(ir,inow) + 3*&
                     NRiS(ir)/2
                dist = MOD(dist, NRiS(ir)) - nris(ir)/2
                ! dist is now the relative grid point displacement of state is
                ! with respect to state inow.
                IF (ABS(dist) .LE. maxdist(ir)) THEN
                   neighborlist(ir,istot) = is
                   distsum = distsum + dist
                   count = count + 1
                   istot = istot + 1
                   done(is) = .TRUE.
                   ! if (parent)
                   ! &                    WRITE(6,*)'*** FOUND CLOSE: ir=',ir,'; is=',is,
                   ! &                    '; llc=',lower_left(ir,is)
                ELSEIF (ABS(ABS(dist) - 2*maxdist(ir)) .LT.&
                     disttonext) THEN
                   inext = is
                   disttonext = ABS(ABS(dist) - 2*maxdist(ir))
                ENDIF
             ENDIF
          ENDDO
          ! NOW the loop is finished. Bring together all those states that
          ! have been found to be close enough:
          newllc = MOD(lower_left(ir,inow) + NINT(distsum/REAL(count,kind=real_8))&
               + 2*NRiS(ir)-1, nris(ir))+1
          ! if (parent) WRITE(6,*)'*** NEW AVG: ir=',ir,
          ! &           '; llc=',newllc
          DO is=istot-count,istot-1
             iss = neighborlist(ir,is)
             lower_left(ir,iss) = newllc
          ENDDO
       ENDDO
    ENDDO


    IF (paral%io_parent) THEN
       DO ir=1,3
          WRITE(6,*)'*** ir=',ir
          DO is=1,nstate
             WRITE(6,*)' old/new (is=',is,&
                  '):  ',old_llc(ir,is), lower_left(ir,is)
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE(old_llc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE optimize_llc_indep

  ! ==================================================================

END MODULE nmr_full_p_utils


! ==================================================================
#define mindistroutine calc_min_dist_alat_smooth
! define mindistroutine calc_min_dist_alat
! ==================================================================
SUBROUTINE do_full(i_state,iB,dowfopt,h1psi0,c0,c1,scr,psi,rhoe,&
     eirop,eivps,z11,nstate)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw
  USE parac, ONLY : paral,parai
  USE soft , ONLY:soft_com
  USE response_pmod , ONLY:nmr_options,response1,timetag
  USE nmr_util_p_utils, ONLY : make_123
  USE nmr_shift_p_utils, ONLY : calc_l_one
  USE nmr_chi_p_utils, ONLY : calc_chi_l_one
  USE nmr_util_p_utils, ONLY : printtime
  USE nmr_position_p_utils, ONLY : calc_min_dist_alat_smooth
  USE perturbation_p_utils, ONLY : ortho_p
  USE rwfopt_p_utils, ONLY : rwfopt_p
  USE zeroing_utils,                   ONLY: zeroing
  USE fft_maxfft,                      ONLY: maxfftn
  IMPLICIT NONE
  INTEGER                                    :: i_state, iB
  LOGICAL                                    :: dowfopt
  COMPLEX(real_8)                            :: scr(*), psi(maxfftn)
  REAL(real_8)                               :: rhoe(*)
  COMPLEX(real_8)                            :: eirop(1), eivps(1)
  REAL(real_8)                               :: z11(*)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c1(ncpw%ngw,nstate,9), &
                                                c0(ncpw%ngw,nstate), &
                                                h1psi0(ncpw%ngw,nstate)

  CHARACTER(len=1), DIMENSION(3), PARAMETER  :: xyz = (/'x','y','z'/)

  COMPLEX(real_8)                            :: cdummy(1)
  INTEGER                                    :: i_c1, ii, iii, is, isub
  INTEGER, DIMENSION(3), SAVE                :: i_state_previous = (/0,0,0/)
  REAL(real_8)                               :: dij(3), rdummy(1)

89 FORMAT ("*  FULL calculation: state ",i4,&
       ", field in O",a1,19x,"*")
98 FORMAT (52("*"),"NMR*RESPONSE*")

  i_c1 = iB+3

  ! Initialization / very first call:
  IF (i_state_previous(iB) .EQ. 0) THEN
     CALL zeroing(c1(:,:,i_c1))!,ngw*nstate)
     response1%t_initialguess_c1 = .FALSE.
     i_state_previous(iB) = i_state
  ELSE
     response1%t_initialguess_c1 = .TRUE.
  ENDIF


  ! Preparation of B_iB . (d_k - d_fix) x p |psi_0_k>:
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE (6,98)
     IF (paral%io_parent)&
          WRITE (6,89) i_state,xyz(iB)
     IF (paral%io_parent)&
          WRITE (6,98)
  ENDIF
  CALL make_123(iB,ii,iii)
  IF (dowfopt) THEN
     CALL zeroing(h1psi0)!,ngw*nstate)
     DO is=1,nstate      ! apply dd_ii p_iii - dd_iii p_ii:
        ! with dd = d_is - d_istate
        CALL mindistroutine&
             (is,i_state,dij)
        CALL apply_op_p (c0(:,is), scr, ii, ncpw%ngw)
        CALL daxpy(2*ncpw%ngw, dij(iii),scr,1,h1psi0(1,is),1)
        CALL apply_op_p (c0(1,is), scr, iii, ncpw%ngw)
        CALL daxpy(2*ncpw%ngw,-dij(ii),scr,1,h1psi0(1,is),1)
     ENDDO
     timetag='preparation.'
     IF (paral%parent) CALL printtime
     IF (paral%io_parent)&
          WRITE (6,98)
     IF (soft_com%exsoft) RETURN
  ENDIF


  ! Prepare a good guess for C1 from the previous optimization:
  ! this should be:
  ! B  .  {d_k - d_i} x |psi1^p>
  ! where k is the new and i is the old correction state.
  ! B = B * O_beta
  ! For an unknown reason, the convergence is better WITHOUT this
  ! adjustment. Bill Gates alone may have an idea why.
  IF (i_state .NE. i_state_previous(iB)) THEN
     CALL mindistroutine&
          (i_state,i_state_previous(iB),dij)
     CALL daxpy(2*nstate*ncpw%ngw, dij(ii),c1(1,1,iii),1,c1(1,1,i_c1),1)
     CALL daxpy(2*nstate*ncpw%ngw,-dij(iii),c1(1,1,ii),1,c1(1,1,i_c1),1)
     CALL ortho_p(nstate,c0,c1(1,1,i_c1))
  ENDIF
  IF (dowfopt) THEN
     CALL tiset('PSI1-WF-OPT',isub)
     CALL rwfopt_p(c0,c1(1,1,i_c1),psi,rhoe,rdummy,eirop,&
          eivps,cdummy,h1psi0,z11,nstate,cdummy)
     CALL tihalt('PSI1-WF-OPT',isub)
     timetag='calculation of psi^1.'
     IF (paral%parent) CALL printtime
  ENDIF

  IF (nmr_options%tcurrent) THEN
     CALL daxpy(2*ncpw%ngw, 1._real_8, c1(1,i_state,i_c1),1,&
          c1(1,i_state,i_c1+3),1)
  ENDIF

  IF (paral%io_parent)&
       WRITE (6,98)
  CALL calc_l_one(c0,c1(1,1,i_c1),nstate,iB,psi,rhoe,i_state)
  timetag='shift contribution.'
  IF (paral%parent) CALL printtime

  CALL calc_chi_l_one(c0,c1(1,1,i_c1),nstate,iB,psi,i_state)
  timetag='susceptibility contribution.'
  IF (paral%parent) CALL printtime


  i_state_previous(iB) = i_state
  RETURN
END SUBROUTINE do_full
! ==================================================================
