MODULE ks_ener_p_utils
  USE cnst,                            ONLY: ry
  USE error_handling,                  ONLY: stopgm
  USE frsblk_c_utils,                  ONLY: reorder_c
  USE kinds,                           ONLY: real_8
  USE kpert_util_p_utils,              ONLY: hamofk_p,&
                                             s12k_p
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE parac,                           ONLY: paral
  USE response_pmod,                   ONLY: ekpert,&
                                             fkpert,&
                                             hamilk
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ks_ener_p
  PUBLIC :: give_scr_ks_ener_p

CONTAINS

  SUBROUTINE  ks_ener_p(cu1,z11,nstate)
    ! ==-------------------------------------------------------------==


    ! ARGUMRNTS
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: cu1(ncpw%ngw,nstate,3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ks_ener_p'

    COMPLEX(real_8), ALLOCATABLE             :: c1k(:,:), work(:)
    INTEGER                                  :: i, ierr, ig, ikk, il_hh, &
                                                istate, isub, j, lwork, n2, &
                                                nleft
    REAL(real_8)                             :: amu, bottom, ksquare
    REAL(real_8), ALLOCATABLE                :: a_temp(:), aux(:), hh(:), &
                                                s12(:,:), w(:), z(:,:)

! VARIABLE
! Diagonalization

    CALL tiset('   KS_ENER_P',isub)

    n2 = nstate*nstate
    amu = 2000.0_real_8
    bottom = 2000.0_real_8

    lwork = MAX(1,2*2*nstate-1)


    ! ==-------------------------------------------------------------==
    ! ==   Allocate Memory for the Hamiltonian and the KS energy     ==
    ! ==   the KS H is complex and 2*NSTATE x 2*NSTATE --> HAMILK    ==
    ! ==   The KS energy are 2*NSTATE*NKPTS           --> EKPERT     ==
    ! ==   The occupation numbers                     --> FKPERT
    ! ==-------------------------------------------------------------==

    ALLOCATE(hamilk(2*nstate,2*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ekpert(2*nstate,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ekpert)!,2*nstate*nkpt%nkpts)
    ALLOCATE(fkpert(2*nstate,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(fkpert)!,2*nstate*nkpt%nkpts)

    ! == ovlap matrix to -1/2 (depending on k vector)
    ALLOCATE(s12(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! == gregarian matrix
    il_hh = 2*nstate*nstate
    ALLOCATE(hh(il_hh),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(hh)!,il_hh)

    ALLOCATE(c1k(ncpw%ngw, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(z(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(a_temp(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*2*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO ikk=1,nkpt%nkpts
       IF (paral%io_parent)&
            WRITE(6,'("K VECTOR ",i10)') ikk

       ! HAMILK = (.0,.0) S12 = (.0)
       CALL zeroing(hamilk)!,2*nstate*2*nstate)
       CALL zeroing(s12)!,nstate*nstate)

       ! ==-------------------------------------------------------------==
       ! ==     Wavefunction not orthogonalized yet                     ==
       ! ==-------------------------------------------------------------==
       DO istate = 1,nstate
          DO ig = 1,ncpw%ngw
             c1k(ig,istate)= (cu1(ig,istate,1)*rk(1,ikk) +&
                  cu1(ig,istate,2)*rk(2,ikk)  +&
                  cu1(ig,istate,3)*rk(3,ikk) )*parm%tpiba
          ENDDO
       ENDDO
       ! ==-------------------------------------------------------------==
       ! ==     Compute the matrix S^(-1/2)[K]                          ==
       ! ==-------------------------------------------------------------==

       CALL s12k_p(c1k,s12,a_temp,w,z,aux,nstate)

       ! ==-------------------------------------------------------------==
       ! ==     Compute     K^2                                         ==
       ! ==-------------------------------------------------------------==

       ksquare = rk(1,ikk)*rk(1,ikk)+rk(2,ikk)*rk(2,ikk)+&
            rk(3,ikk)*rk(3,ikk)

       ! ==-------------------------------------------------------------==
       ! ==        FINAL HAMILTONIAN MATRIX OF K                        ==
       ! ==-------------------------------------------------------------==

       CALL hamofk_p(z11,s12,nstate,ikk,ksquare,z)

       ! ==                     END LOOPS                                   ==
       ! ==-----------------------------------------------------------------==
       ! ==-----------------------------------------------------------------==
       ! HAMILT IS HERMITIAN .....


       ! DIAGONALIZATION
       CALL zheev('V','U',2*nstate,hamilk,2*nstate,ekpert(1,ikk),work,&
            lwork,aux,ierr)
       IF (lwork.LT.INT(aux(1))) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' WARNING| LWORK SUBOPTIMAL FOR ZHEEV'
          IF (paral%io_parent)&
               WRITE(6,*) '   LWORK:',lwork
          IF (paral%io_parent)&
               WRITE(6,*) ' OPTIMAL:',INT(aux(1))
       ENDIF
       CALL reorder_c(2*nstate,2*nstate,ekpert(1,ikk),hamilk)

       ! Order Roots with respect to energy
       CALL dscal(2*nstate,-1._real_8,ekpert(1,ikk),1)
       DO i = 1 ,nstate
          fkpert(i,ikk) = 2.0_real_8
          fkpert(2*nstate-i+1,ikk) = 0.0_real_8
       ENDDO
       DO i = nstate+1,2*nstate
          amu = MIN(amu,ekpert(i,ikk))
       ENDDO
       DO i = 1,nstate
          bottom = MIN(bottom,ekpert(i,ikk))
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'("BOTTOM LEVEL",F14.7)') bottom*2*ry
    IF (paral%io_parent)&
         WRITE(6,'("MINIMUM CB",F14.7)') amu*2*ry

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') ' EIGENVALUES(EV) AND OCCUPATION:'
       DO ikk=1,nkpt%nkpts
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I8,4F12.6)')&
               'K POINT:',ikk,(rk(i,ikk),i=1,3),wk(ikk)
          DO i=1,2*nstate,2
             nleft=MIN(2*nstate-i,1)
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                  (i+j,ekpert(i+j,ikk)*(2*ry),fkpert(i+j,ikk),j=0,nleft)
          ENDDO
       ENDDO
    ENDIF

    ! ==-----------------------------------------------------------------==
    DEALLOCATE(c1k,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(z,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(a_temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==-----------------------------------------------------------------==
    DEALLOCATE(hamilk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ekpert,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hh,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(s12,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==-----------------------------------------------------------------==

    RETURN
  END SUBROUTINE ks_ener_p
  ! ======================================================================
  ! ======================================================================
  SUBROUTINE give_scr_ks_ener_p(lks_ener_p,n2,nstate)

    INTEGER                                  :: lks_ener_p, n2, nstate

    INTEGER                                  :: il_a_temp, il_aux, il_c1k, &
                                                il_w, il_z, lwork

    il_c1k  = 2*ncpw%ngw*nstate
    il_z = n2
    il_w = nstate
    il_a_temp = n2
    il_aux = 3*2*nstate
    lwork = MAX(1,2*2*nstate-1)

    lks_ener_p = il_z+il_w+il_a_temp+il_aux+il_c1k+2*lwork

    RETURN
  END SUBROUTINE give_scr_ks_ener_p
  ! ======================================================================
  ! ======================================================================

END MODULE ks_ener_p_utils
