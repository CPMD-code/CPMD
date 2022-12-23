MODULE kpert_util_p_utils
  USE csmat_utils,                     ONLY: csmat
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk
  USE response_pmod,                   ONLY: h0_11,&
                                             h1_00,&
                                             h1_10,&
                                             h2_00,&
                                             hamilk
  USE sfac,                            ONLY: fnl
  USE system,                          ONLY: ncpw,&
                                             parm
  USE utils,                           ONLY: dspevy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: s12k_p
  PUBLIC :: hamofk_p

CONTAINS

  SUBROUTINE  s12k_p(c1k,s12,a_temp,w,z,aux,nstate)
    ! ==-------------------------------------------------------------==
    ! ==                     COMPUTES:                               ==
    ! ==   The inverse of the overlap matrix for the set of wfn      ==
    ! ==   relative to the K-vector IKIND                            ==
    ! ==-------------------------------------------------------------==

    ! ==-------------------------------------------------------------==
    ! ==     Compute the matrix S^(-1/2)[K]                          ==
    ! ==-------------------------------------------------------------==


    INTEGER                                  :: nstate
    REAL(real_8) :: aux(3*2*nstate), z(nstate,nstate), w(nstate), &
      a_temp(nstate*nstate), s12(nstate,nstate)
    COMPLEX(real_8)                          :: c1k(ncpw%ngw,nstate)

    INTEGER                                  :: i, j, k, l, l0, m
    REAL(real_8)                             :: fac

! VARIABLES
! real(8), allocatable :: FNL(:)
! real(8), allocatable :: DFNL(:)

    CALL csmat(s12,c1k,fnl,nstate,1)

    k = 0
    DO j = 1,nstate
       DO i = j,nstate
          k=k+1
          a_temp(k) = s12(i,j)
       ENDDO
    ENDDO

    CALL dspevy(1,a_temp,w,z,nstate,nstate,aux,3*nstate)

    CALL zeroing(a_temp)!,nstate*nstate)
    DO m = 1,nstate
       w(m) = 1._real_8/SQRT(w(m))
       l = 1
       DO j = 1,nstate
          fac = w(m)*z(j,m)
          CALL daxpy(nstate,fac,z(1,m),1,a_temp(l),1)
          l = l+nstate
       ENDDO
    ENDDO
    l = 1
    DO j = 1,nstate
       DO i = 1,nstate
          s12(i,j) = a_temp(l)
          l = l+1
       ENDDO
    ENDDO

    l0=0
    DO j = 1,nstate
       DO i = j,nstate
          l=l0+i
          s12(i,j) = a_temp(l)
          s12(j,i) = s12(i,j)
       ENDDO
       l0=l0+nstate
    ENDDO

    RETURN
  END SUBROUTINE s12k_p
  ! ==---------------------------------------------------------------==


  ! ==---------------------------------------------------------------==
  SUBROUTINE  hamofk_p(z11,s12,nstate,ikk,ksquare,z)
    ! ==---------------------------------------------------------------==
    ! ==                       COMPUTES:                               ==
    ! ==  The       FINAL HAMILTONIAN MATRIX OF K                      ==
    ! ==---------------------------------------------------------------==


    ! ARGUMRNTS
    INTEGER                                  :: nstate
    REAL(real_8)                             :: s12(nstate,nstate), &
                                                z11(nstate,nstate)
    INTEGER                                  :: ikk
    REAL(real_8)                             :: ksquare, z(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hamofk_p'

    COMPLEX(real_8)                          :: ck2
    INTEGER                                  :: i, ierr, ik1, ik2, istate, &
                                                jstate
    REAL(real_8)                             :: cfac(3,3), fac1(3,3), SumIm, &
                                                SumRe
    REAL(real_8), ALLOCATABLE                :: b_temp(:,:), c_temp(:,:), &
                                                hh(:,:,:)

! VARIABLES
! dimension wk(NKPTS)
! dimension RK(3,NKPTS)
! dimension H1_00(NSTATE,NSTATE,3)
! dimension H2_00(NSTATE,NSTATE,3,3)
! dimension H1_10(NSTATE,NSTATE,3,3)
! dimension H0_11(NSTATE,NSTATE,3,3)

    ALLOCATE(hh(2, nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(b_temp(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(c_temp(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(hh)!,2*nstate*nstate)


    ! ==-------------------------------------------------------------==
    ! ==     First block <c0|H|c0> = (z11 - (1/2)K^2) d_ij +         ==
    ! ==           kx    H1_00(x) +Ky    H1_00(y) +Kz    H1_00(z)    ==
    ! ==           Kx Kx H2_00(xx)+Kx Ky H2_00(xy)+Kx Kz H2_00(xz)   ==
    ! ==           Ky Kx H2_00(yx)+Ky Ky H2_00(yy)+Ky Kz H2_00(yz)   ==
    ! ==           Kz Kx H2_00(zx)+Kz Ky H2_00(zy)+Kz Kz H2_00(zz)   ==
    ! ==-------------------------------------------------------------==

    DO istate = 1,nstate
       hamilk(istate,istate) = CMPLX(-z11(istate,istate)*0.5_real_8-0.5_real_8*&
            parm%tpiba2*ksquare,0.0_real_8,kind=real_8)
       DO jstate = istate+1,nstate
          DO ik1 = 1,3
             hamilk(jstate,istate)= hamilk(jstate,istate)+parm%tpiba*rk(&
                  ik1,ikk)*h1_00(jstate,istate,ik1)
             DO ik2 = 1,3
                hamilk(jstate,istate)= hamilk(jstate,istate)+parm%tpiba2*&
                     rk(ik1,ikk)*rk(ik2,ikk)*h2_00(jstate,istate,ik1,ik2)
             ENDDO
          ENDDO
          hamilk(istate,jstate)= CONJG(hamilk(jstate,istate))
       ENDDO
    ENDDO


    ! ==-------------------------------------------------------------==
    ! ==     Second block <c1|H|c0>{1j} = -i *( Sum_ab [             ==
    ! ==           Ka Kb Sum_l [S^(-1/2){il} H1_10(ab){lj}]])        ==
    ! ==            matrix product                                   ==
    ! ==-------------------------------------------------------------==


    CALL zeroing(hh)!,2*nstate*nstate)
    DO ik1 = 1,3
       DO ik2 = 1,3
          ! K-dependent factor
          ! CFAC(IK2,IK1) = -uimag*RK(IK1,IKK)*RK(IK2,IKK)*TPIBA2
          cfac(ik2,ik1) = -rk(ik1,ikk)*rk(ik2,ikk)*parm%tpiba2
       ENDDO
    ENDDO

    DO istate = 1,nstate
       DO jstate = 1,nstate
          DO ik1 = 1,3
             DO ik2 =  1,3
                SumRe = REAL(h1_10(jstate,istate,ik2,ik1))
                SumIm = AIMAG(h1_10(jstate,istate,ik2,ik1))

                hh(1,jstate,istate)=hh(1,jstate,istate)-&
                     SumIm*cfac(ik2,ik1)
                hh(2,jstate,istate)=hh(2,jstate,istate)+&
                     SumRe*cfac(ik2,ik1)
                ! HH(JSTATE,ISTATE)=HH(JSTATE,ISTATE)+
                ! &                    CFAC(IK2,IK1)*H1_10(JSTATE,ISTATE,IK2,IK1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO istate = 1,nstate
       DO jstate = 1,nstate

          SumRe = 0.0_real_8
          SumIm = 0.0_real_8

          DO i = 1,nstate
             SumRe = sumre + s12(i,jstate)*hh(1,i,istate)
             SumIm = sumim + s12(i,jstate)*hh(2,i,istate)
          ENDDO
          hamilk(jstate+nstate,istate) = hamilk(jstate+nstate,istate)+&
               CMPLX(SumRe,SumIm,kind=real_8)
          hamilk(jstate,istate+nstate) = CONJG(hamilk(jstate+nstate,&
               istate))

       ENDDO
    ENDDO

    ! 

    ! ==-------------------------------------------------------------==
    ! ==     Third  block <c1|H|c1>{1j} = ( Sum_ab [                 ==
    ! ==   Ka Kb Sum_lm [S^(-1/2){il} H1_10(ab){lm} S^(-1/2){mj}]])  ==
    ! ==         double   matrix product                             ==
    ! ==-------------------------------------------------------------==

    ck2 = CMPLX(parm%tpiba2*ksquare*0.5_real_8,0.0_real_8,kind=real_8)

    CALL zeroing(z)!,nstate*nstate)
    DO ik1 = 1,3
       DO ik2 = 1,3
          fac1(ik2,ik1) = rk(ik1,ikk)*rk(ik2,ikk)*parm%tpiba2
       ENDDO
    ENDDO

    DO istate = 1,nstate
       DO jstate = 1,nstate
          DO ik1 = 1,3
             DO ik2 =  1,3
                z(jstate,istate)=z(jstate,istate)+&
                     fac1(ik2,ik1)*h0_11(jstate,istate,ik2,ik1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO



    ! ==--------------------------------------------------------------==
    ! ==         N^4 loop parted in 3 loops (sebast)                  ==


    CALL zeroing(c_temp)!,nstate*nstate)
    CALL zeroing(b_temp)!,nstate*nstate)

    CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,z,nstate,s12,nstate,&
         0.0_real_8,c_temp,nstate)

    CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,s12,nstate,c_temp,&
         nstate,0.0_real_8,b_temp,nstate)


    DO istate = 1,nstate
       hamilk(istate+nstate,istate+nstate) =hamilk(istate+nstate,&
            istate+nstate) - ck2

       DO jstate = 1,nstate
          hamilk(jstate+nstate,istate+nstate) = hamilk(jstate+nstate,&
               istate+nstate)+CMPLX(b_temp(jstate,istate),0.0_real_8,kind=real_8)
       ENDDO
    ENDDO
    DEALLOCATE(hh,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(b_temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(c_temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE hamofk_p



END MODULE kpert_util_p_utils
