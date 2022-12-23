#define split_suscept_contribs 0

MODULE nmr_chi_p_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nmr_position_p_utils,            ONLY: set_origin
  USE nmr_util_p_utils,                ONLY: epsi,&
                                             make_123
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: chi,&
                                             chi_factor,&
                                             chi_si_iso,&
                                             chi_si_to_ppmcgs
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_chi_l
  PUBLIC :: calc_chi_l_one
  PUBLIC :: calc_chi_p
  PUBLIC :: print_chi_tensor

CONTAINS

  ! ==================================================================
  SUBROUTINE calc_chi_l(c0,c1,nstate,b,psi)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c1(ncpw%ngw,*)
    INTEGER                                  :: nstate, b
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_chi_l'

    INTEGER                                  :: a, ierr, is, isub, k, l
    REAL(real_8)                             :: contrib
    REAL(real_8), ALLOCATABLE                :: scr(:,:)
    REAL(real_8), EXTERNAL                   :: dotp

    CALL tiset('chi-L     ',isub)
    ALLOCATE(scr(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! =----------------------------------------------------------------=
    ! ==  x_ab = Int[ r x (phi1_b grad phi0 - phi0 grad phi1_b) ]_a   ==
    ! ==  where phi1_a is the response to a perturbation L_a.         ==
    ! ==  The two contributions are equal.                            ==
    ! =----------------------------------------------------------------=
    ! scr(*,k) = r_k psi1_b   (b=fix)
    CALL zeroing(scr)!,3*nnr1)
    DO is=1,nstate
       CALL set_origin(is)
       CALL ffttor(c1(1,is),scr,psi,ncpw%ngw,.FALSE.)
       CALL apply_op_rx(scr,1._real_8,scr(1,3),0._real_8,3)
       CALL apply_op_rx(scr,1._real_8,scr(1,2),0._real_8,2)
       CALL apply_op_rx(scr,1._real_8,scr(1,1),0._real_8,1)
       CALL fft2tog(scr,scr,scr(1,2),scr(1,2),psi,ncpw%ngw,.FALSE.)
       CALL ffttog  (scr(1,3),scr(1,3),psi,ncpw%ngw,.FALSE.)
       DO l=1,3
          DO k=1,3
             CALL apply_op_p(c0(1,is),psi,l,ncpw%ngw)
             contrib = crge%f(is,1)*dotp(ncpw%ngw,psi,scr(1,k))
             DO a=1,3
                chi(a,b) = chi(a,b) +&
                     2.0_real_8*epsi(a,k,l)*chi_factor*contrib ! F(occ) included
                ! !   &            2*2* epsi(a,k,l)*chi_factor*contrib
                ! Factor 2: the two terms (see above) are equal.
                ! Factor 2: double-occupation of orbitals.
             ENDDO         ! a
          ENDDO             ! k
       ENDDO                 ! l
    ENDDO                     ! is
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('chi-L     ',isub)
    RETURN
  END SUBROUTINE calc_chi_l
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_chi_l_one(c0,c1,nstate,b,psi,istate)
    ! IDENTICAL to calc_chi_L, except that the stuff is only done for istate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c1(ncpw%ngw,*)
    INTEGER                                  :: nstate, b
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: istate

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_chi_l_one'

    INTEGER                                  :: a, ierr, isub, k, l
    REAL(real_8)                             :: contrib
    REAL(real_8), ALLOCATABLE                :: scr(:,:)
    REAL(real_8), EXTERNAL                   :: dotp

    CALL tiset('chi-L     ',isub)
    ALLOCATE(scr(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! =----------------------------------------------------------------=
    ! ==  x_ab = Int[ r x (phi1_b grad phi0 - phi0 grad phi1_b) ]_a   ==
    ! ==  where phi1_a is the response to a perturbation L_a.         ==
    ! ==  The two contributions are equal.                            ==
    ! =----------------------------------------------------------------=
    ! scr(*,k) = r_k psi1_b   (b=fix)
    CALL zeroing(scr)!,3*nnr1)
    CALL set_origin(istate)
    CALL ffttor(c1(1,istate),scr,psi,ncpw%ngw,.FALSE.)
    CALL apply_op_rx(scr,1._real_8,scr(1,3),0._real_8,3)
    CALL apply_op_rx(scr,1._real_8,scr(1,2),0._real_8,2)
    CALL apply_op_rx(scr,1._real_8,scr(1,1),0._real_8,1)
    CALL fft2tog(scr,scr,scr(1,2),scr(1,2),psi,ncpw%ngw,.FALSE.)
    CALL ffttog  (scr(1,3),scr(1,3),psi,ncpw%ngw,.FALSE.)
    DO l=1,3
       DO k=1,3
          CALL apply_op_p(c0(1,istate),psi,l,ncpw%ngw)
          contrib = crge%f(istate,1)*dotp(ncpw%ngw,psi,scr(1,k))
          DO a=1,3
             chi(a,b) = chi(a,b) +&
                  2.0_real_8*epsi(a,k,l)*chi_factor*contrib  ! F(occ) included
             ! !   &          2*2* epsi(a,k,l) * chi_factor * contrib
             ! Factor 2: the two terms (see above) are equal.
             ! Factor 2: double-occupation of orbitals.
          ENDDO             ! a
       ENDDO                 ! k
    ENDDO                     ! l
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('chi-L     ',isub)
    RETURN
  END SUBROUTINE calc_chi_l_one
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE calc_chi_p(c0,c1,nstate,j,psi,gradpsi)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c1(ncpw%ngw,*)
    INTEGER                                  :: nstate, j
    COMPLEX(real_8)                          :: psi(:), gradpsi(ncpw%ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_chi_p'

    INTEGER                                  :: a, b, i1, i2, ierr, is, isub, &
                                                k, l1, l2
    REAL(real_8)                             :: chi_diff(3,3), chi_tmp(3,3), &
                                                contrib
    REAL(real_8), ALLOCATABLE                :: scr(:,:)
    REAL(real_8), EXTERNAL                   :: dotp

    CALL tiset('chi-p     ',isub)
    ALLOCATE(scr(fpar%nnr1,4),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(scr)
    ! =----------------------------------------------------------------=
    ! ==  x_ab = Int[ r x (phi1_b grad phi0 - phi0 grad phi1_b) ]_a   ==
    ! ==  where phi1_a = [ r x phi1_p ]_a                             ==
    ! ==  and phi1_p is the vector containing the three responses     ==
    ! ==  to perturbations p_a.                                       ==
    ! ==  The two contributions are NOT equal.                        ==
    ! =----------------------------------------------------------------=
    ! scr(1)     = c0 in R space
    ! scr(2)     = r_i1 r_k c0 in R space
    ! scr(3)     = r_i2 r_k c0 in R space
    ! gradpsi(l) = grad_l psi1_p_j
    ! (-> first term)
    ! scr(4)     = c1 in R space
    ! scr(2)     = r_i1 r_k c1 in R space
    ! scr(3)     = r_i2 r_k c1 in R space
    ! gradpsi(l) = grad_l psi0_p_j
    ! (-> second term)
    ! =----------------------------------------------------------------=
    CALL make_123(j,i1,i2)
    CALL zeroing(chi_tmp)!,3*3)
    DO is=1,nstate
       CALL set_origin(is)
       CALL fft2tor(c0(1,is),scr,c1(1,is),scr(1,4),psi,ncpw%ngw,.FALSE.)
       DO k=1,3
          CALL make_123(k,l1,l2)
          ! ==  first term
          CALL apply_op_p(c1(1,is),gradpsi(1,1),l1,ncpw%ngw)
          CALL apply_op_p(c1(1,is),gradpsi(1,2),l2,ncpw%ngw)
          CALL apply_op_rx(scr(1,1),1._real_8,scr(1,2),0._real_8, k)
          CALL apply_op_rx(scr(1,2),1._real_8,scr(1,3),0._real_8,i2)
          CALL apply_op_rx(scr(1,2),1._real_8,scr(1,2),0._real_8,i1)
          CALL fft2tog(scr(1,2),scr(1,2),scr(1,3),scr(1,3),&
               psi,ncpw%ngw,.FALSE.)
          DO a=1,3
             DO b=1,3
                ! now, i1 means scr(2),     i2 means scr(3),
                ! and  l1 means gradpsi(1), l2 means gradpsi(2):
                chi_tmp(a,b) = chi_tmp(a,b) +&
                                ! !   &            2 * epsi(a,k,l1) * epsi(b,i1,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l1) * epsi(b,i1,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,2),gradpsi(1,1))
                chi_tmp(a,b) = chi_tmp(a,b) +&
                                ! !   &            2 * epsi(a,k,l1) * epsi(b,i2,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l1) * epsi(b,i2,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,3),gradpsi(1,1))
                chi_tmp(a,b) = chi_tmp(a,b) +&
                                ! !   &            2 * epsi(a,k,l2) * epsi(b,i1,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l2) * epsi(b,i1,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,2),gradpsi(1,2))
                chi_tmp(a,b) = chi_tmp(a,b) +&
                                ! !   &            2 * epsi(a,k,l2) * epsi(b,i2,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l2) * epsi(b,i2,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,3),gradpsi(1,2))
             ENDDO         ! b
          ENDDO             ! a
       ENDDO                 ! k

#if split_suscept_contribs
88     FORMAT (3f15.7)
       IF (paral%io_parent)&
            WRITE (6,'(A,I4)') 'CHI contrib 1 state ',is
       DO a=1,3
          IF (paral%io_parent)&
               WRITE (6,88) (chi_tmp(a,b),b=1,3)
       ENDDO
#endif
       CALL daxpy(3*3,1._real_8,chi_tmp,1,chi,1)
       CALL zeroing(chi_tmp)!,3*3)

       DO k=1,3
          CALL make_123(k,l1,l2)

          ! ==  second term
          CALL apply_op_p(c0(1,is),gradpsi(1,1),l1,ncpw%ngw)
          CALL apply_op_p(c0(1,is),gradpsi(1,2),l2,ncpw%ngw)
          CALL apply_op_rx(scr(1,4),1._real_8,scr(1,2),0._real_8,k)
          CALL apply_op_rx(scr(1,2),1._real_8,scr(1,3),0._real_8,i2)
          CALL apply_op_rx(scr(1,2),1._real_8,scr(1,2),0._real_8,i1)
          CALL fft2tog(scr(1,2),scr(1,2),scr(1,3),scr(1,3),&
               psi,ncpw%ngw,.FALSE.)
          DO a=1,3
             DO b=1,3
                ! now, i1 means scr(2),     i2 means scr(3),
                ! and  l1 means gradpsi(1), l2 means gradpsi(2):

                chi_tmp(a,b) = chi_tmp(a,b) -&
                                ! !   &            2 * epsi(a,k,l1) * epsi(b,i1,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l1) * epsi(b,i1,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,2),gradpsi(1,1))
                chi_tmp(a,b) = chi_tmp(a,b) -&
                                ! !   &            2 * epsi(a,k,l1) * epsi(b,i2,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l1) * epsi(b,i2,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,3),gradpsi(1,1))
                chi_tmp(a,b) = chi_tmp(a,b) -&
                                ! !   &            2 * epsi(a,k,l2) * epsi(b,i1,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l2) * epsi(b,i1,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,2),gradpsi(1,2))
                chi_tmp(a,b) = chi_tmp(a,b) -&
                                ! !   &            2 * epsi(a,k,l2) * epsi(b,i2,j) * chi_factor
                     crge%f(is,1)* epsi(a,k,l2) * epsi(b,i2,j) * chi_factor&
                     * dotp(ncpw%ngw,scr(1,3),gradpsi(1,2))
             ENDDO         ! a
          ENDDO             ! b
       ENDDO                 ! k


#if split_suscept_contribs
       IF (paral%io_parent)&
            WRITE (6,'(A,I4)') 'CHI contrib 2 state ',is
       DO a=1,3
          IF (paral%io_parent)&
               WRITE (6,88) (chi_tmp(a,b),b=1,3)
       ENDDO
#endif

       CALL daxpy(3*3,1._real_8,chi_tmp,1,chi,1)
       CALL zeroing(chi_tmp)!,3*3)

#if split_suscept_contribs
       ! NOW THE THEORETICAL DIFFERENCE BETWEEN THE FIRST TWO TERMS:
       CALL zeroing(chi_diff)!,3*3)
       DO k=1,3
          CALL apply_op_rx(scr(1,1),1._real_8,scr(1,2),0._real_8,k)
          CALL ffttog(scr(1,2),scr(1,2),psi,ngw,.FALSE.)
          contrib=dotp(ngw,c1(1,is),scr(1,2)) * chi_factor

          DO a=1,3
             DO b=1,3
                DO l1=1,3 ! corresponds to "beta"
                   chi_diff(a,b)=chi_diff(a,b)+&
                                ! !   &              2* epsi(l1,b,k)*epsi(l1,j,a)*contrib
                        f(is,1)*epsi(l1,b,k)*epsi(l1,j,a)*contrib
                ENDDO     ! l1
             ENDDO         ! a
          ENDDO             ! b
       ENDDO                 ! k


       IF (paral%io_parent)&
            WRITE (6,'(A,I4)') 'CHI DIFF      state ',is
       DO a=1,3
          IF (paral%io_parent)&
               WRITE (6,88) (chi_diff(a,b),b=1,3)
       ENDDO
       ! END of differences...
#endif

    ENDDO                     ! is
    ! Factor 2: double-occupation of orbitals.
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('chi-p     ',isub)
    RETURN
  END SUBROUTINE calc_chi_p
  ! ==================================================================
  SUBROUTINE print_chi_tensor
    INTEGER                                  :: a, b, info
    REAL(real_8)                             :: aniso, chi_matrix(3,3), &
                                                eigenvalues(3), scr(3,3), &
                                                work(3*3)

    CALL dcopy(3*3,chi,1,chi_matrix,1)
    CALL mp_sum(chi_matrix,3*3,parai%allgrp)

    IF (.NOT. paral%parent) RETURN
    IF (paral%io_parent)&
         WRITE (6,'(15("=")," MAGNETIC SUSCEPTIBILITY TENSOR ",18("="))')
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  XX=',chi_matrix(1,1),&
         '  XY=',chi_matrix(1,2),'  XZ=',chi_matrix(1,3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  YX=',chi_matrix(2,1),&
         '  YY=',chi_matrix(2,2),'  YZ=',chi_matrix(2,3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  ZX=',chi_matrix(3,1),&
         '  ZY=',chi_matrix(3,2),'  ZZ=',chi_matrix(3,3)

    DO a=1,3
       DO b=1,3
          scr(a,b) = (chi_matrix(a,b) + chi_matrix(b,a))/2._real_8
       ENDDO                 ! b
    ENDDO                     ! a
    CALL dcopy(3*3,scr,1,chi_matrix,1)

    CALL dsyev('N','U', 3, chi_matrix, 3, eigenvalues,&
         work, 9, info)
    IF ((info .NE. 0).AND.paral%io_parent)&
         WRITE(6, *)'ERROR IN DSYEV, code ',info
    chi_SI_iso = (eigenvalues(1)+ eigenvalues(2)+eigenvalues(3))/3._real_8
    aniso      =  eigenvalues(3)-(eigenvalues(2)+eigenvalues(1))/2._real_8
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')' PV1=',eigenvalues(1),&
         ' PV2=',eigenvalues(2),' PV3=',eigenvalues(3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')' iso=',chi_SI_iso,&
         '                   aniso=',aniso
    IF (paral%io_parent)&
         WRITE (6, '(15("=")," (units = 10^-30 J/T^2) ",26("="))')

    CALL dscal(3*3,chi_SI_to_ppmcgs,scr,1)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  XX=',scr(1,1),&
         '  XY=',scr(1,2),'  XZ=',scr(1,3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  YX=',scr(2,1),&
         '  YY=',scr(2,2),'  YZ=',scr(2,3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')'  ZX=',scr(3,1),&
         '  ZY=',scr(3,2),'  ZZ=',scr(3,3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')&
         ' PV1=',chi_SI_to_ppmcgs*eigenvalues(1),&
         ' PV2=',chi_SI_to_ppmcgs*eigenvalues(2),&
         ' PV3=',chi_SI_to_ppmcgs*eigenvalues(3)
    IF (paral%io_parent)&
         WRITE (6, '(3(A,F15.5))')' iso=',chi_SI_to_ppmcgs*chi_SI_iso,&
         '                   aniso=',chi_SI_to_ppmcgs*aniso
    IF (paral%io_parent)&
         WRITE (6, '(15("=")," (units = ppm-cgs)      ",26("="))')


    RETURN
  END SUBROUTINE print_chi_tensor
  ! ==================================================================

END MODULE nmr_chi_p_utils
