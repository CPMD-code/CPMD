MODULE hess_eta_p_utils
  USE adat,                            ONLY: elem
  USE d_mat_p_utils,                   ONLY: d_mat_nonloc
  USE eicalc_utils,                    ONLY: eicalc1
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn
  USE fnonloc_p_utils,                 ONLY: fnonloc_p
  USE implhv,                          ONLY: rs_v,&
                                             sd0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE nlps,                            ONLY: imagp
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00
  USE rhoofr_p_utils,                  ONLY: rhoofr_p
  USE rnlsm_p_utils,                   ONLY: rnlsm3
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE spin,                            ONLY: clsd
  USE symm,                            ONLY: symmi
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hess_eta_p

CONTAINS

  ! ==================================================================
  SUBROUTINE hess_eta_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,&
       z11,nstate,icol,iat,dd_fnl,dtmp)
    ! ==--------------------------------------------------------------==
    ! ccc
    ! ccc

    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), drhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: icol, iat
    REAL(real_8) :: &
      dd_fnl(imagp,ions1%nat,maxsys%nhxs,3,3,nstate,nkpt%nkpnt), dtmp(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hess_eta_p'

    COMPLEX(real_8)                          :: ztmp
    COMPLEX(real_8), ALLOCATABLE             :: eirop1(:), eirop2(:), v1(:), &
                                                v1_l(:), v1_loc(:), v1_nl(:), &
                                                v1_nonloc(:)
    INTEGER                                  :: ia2, iat2, ierr, il_eirop1, &
                                                il_v1_loc, il_v1_nonloc, ir, &
                                                is2, isub, ix, k, k2

! variables
! ==--------------------------------------------------------------==

    CALL tiset('    hess_eta_p',isub)
    il_v1_loc=ncpw%nhg*clsd%nlsd
    il_v1_nonloc=nkpt%ngwk*nstate*nkpt%nkpnt
    il_eirop1=ncpw%nhg

    ALLOCATE(v1_loc(il_v1_loc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_loc)!,SIZE(v1_loc))
    ALLOCATE(v1_nonloc(il_v1_nonloc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_nonloc)!,SIZE(v1_nonloc))
    ! NB: In v1_nonloc, we store \op{v1} |psi0> whereas in
    ! v1_loc, v1(g) is stored. Therefore the different dimensions.
    ALLOCATE(v1_l(il_v1_loc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_l)!,SIZE(v1_l))
    ALLOCATE(v1_nl(il_v1_nonloc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_nl)!,SIZE(v1_nl))


    ALLOCATE(eirop1(il_eirop1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop1)!,SIZE(eirop1))
    ! NB: This array is only needed for phonons! Otherwise, it should not
    ! be touched (filled with 0)
    ALLOCATE(eirop2(il_eirop1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop2)!,SIZE(eirop2))

    ALLOCATE(v1(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ! 
    ! ==--------------------------------------------------------------==
    ! ==         the basic loops for perturbations                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/," ",22("*"),a,22("*"),/)')&
         '   perturbations    '
    ! no symmetrization of gradients in this part!
    symmi%nrot=1
    iat2=0
    IF (paral%parent) THEN
       DO is2=1,ions1%nsp
          DO ia2=1,ions0%na(is2)
             iat2=iat2+1
             IF (paral%io_parent)&
                  WRITE(6,'(" **** atom=",i5,3x,a3,3x,3(f10.8,3x,2x))')&
                  iat2,elem%el(ions0%iatyp(is2)),(dtmp(3*(iat2-1)+k2),k2=1,3)
          ENDDO
       ENDDO
    ENDIF

    iat2=0
    ix=0

    DO is2=1,ions1%nsp
       DO ia2=1,ions0%na(is2)
          iat2=iat2+1
          DO k2=1,3
             ix=3*(iat2-1)+k2

             ztmp=CMPLX(dtmp(ix),0._real_8,kind=real_8)

             ! ccc  The local and nonlocal parts of the perturbation are calculated
             ! ccc  coordinate per coordinate, and then accumulated in the arrays
             ! ccc  v1_loc and v1_nonloc, multiplied by the displacements.

             ! .... calculate local part
             CALL zeroing(v1_l)!,SIZE(v1_l))
             CALL zeroing(eirop1)!,SIZE(eirop1))
             CALL eicalc1(k2,is2,iat2,v1_l,eirop1)
             CALL zaxpy(ncpw%nhg,ztmp,v1_l,1,v1_loc,1)
             IF (cntl%tlsd) CALL dcopy(2*ncpw%nhg,v1_loc(1),1,v1_loc(ncpw%nhg+1),1)
             CALL zaxpy(ncpw%nhg,ztmp,eirop1,1,eirop2,1)

             ! ...  calculate the non-local part
             ! ccc here we do not need to "azzero" v1_nl, because it is made in 
             ! ccc the routine pert_nonloc
             CALL fnonloc_p(c0,psi,&
                  v1_nl,crge%f,nstate,k2,is2,iat2,1)
             CALL daxpy(il_v1_nonloc,dtmp(ix),v1_nl,1,&
                  v1_nonloc,1)

          ENDDO
       ENDDO
    ENDDO

    ! ccc ^^^^^^^^^
    ! CCC ATTENTION; THIS WAY OF IMPLEMENTING THE PERTURBATION AS A COLLECTIVE
    ! CCC MOTION IMPLIES THE ALLOCATION FOR TWICE THE V1_(NON)LOC VECTORS.
    ! CCC IT MAY BE TOO MEMORY DEMANDING!


    ! ...  optimization for c1
    ! ==--------------------------------------------------------------==
    ! ccc  The calling to rwfopt_p is done with eirop2 instead of eirop1. 
    ! ccc In this way all routines called by rwfopt_p will have eirop2, where
    ! ccc they have eirop1. In this way there is no need to modify the code.

    CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
         eirop,eivps,v1_loc,v1_nonloc,&
         z11,nstate,eirop2)
    ! ==--------------------------------------------------------------==
    ! update hessian
    ! calculate the linear perturbation in the density 
    ! \   drhoe=<c0|r><r|c1>+cc
    CALL rhoofr_p(c0,c1,drhoe,psi,nstate)
    ! transform drhoe in g space -> v1
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          v1(ir)=CMPLX(drhoe(ir)+drhoe(ir+fpar%nnr1),0._real_8,kind=real_8)
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          v1(ir)=CMPLX(drhoe(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    CALL  fwfftn(v1,.FALSE.,parai%allgrp)

    CALL rnlsm(c1,nstate,1,1,.TRUE.)
    CALL rnlsm3(c0,nstate,dd_fnl)

    CALL dgemv('N',3*ions1%nat,3*ions1%nat,1._real_8,sd0,3*ions1%nat,&
         dtmp,1,0._real_8,rs_v(1,icol+1),1)

    iat2=0
    DO is2=1,ions1%nsp
       DO ia2=1,ions0%na(is2)
          iat2=iat2+1

          DO k2=1,3
             ix=3*(iat2-1)+k2

             ! ...  local part of the hessian matrix       
             CALL d_mat_locps(rs_v(ix,icol+1),v1,iat2,k2,is2)
             CALL d_mat_nonloc(rs_v(ix,icol+1),fnl,dfnl,fnl00,dfnl00,&
                  crge%f,wk,is2,k,k2,iat2,nstate,1)

          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(v1_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_nonloc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_l,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_nl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('    hess_eta',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hess_eta_p
  ! ==================================================================








END MODULE hess_eta_p_utils
