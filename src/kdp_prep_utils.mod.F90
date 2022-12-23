MODULE kdp_prep_utils
  USE calc_pij_utils,                  ONLY: calc_pij
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             nkpt,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kdp_prep

CONTAINS

  ! ==================================================================
  SUBROUTINE kdp_prep(nkdp,c0,nstate,pkdp,xlambda0,ckdp,vpot,psi,&
       betap,eigv)
    ! ==--------------------------------------------------------------==
    ! Compute < c_i | p | c_j > and <c_i| g^2 |c_j> matrix elements for
    ! the k.p equations.
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nkdp, nstate
    COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)
    REAL(real_8)                             :: pkdp(3,nstate,nstate)
    COMPLEX(real_8)                          :: xlambda0(nstate,nstate), &
                                                ckdp(nkpt%ngwk,nstate)
    REAL(real_8)                             :: vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    REAL(real_8)                             :: betap, eigv(nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'kdp_prep'

    COMPLEX(real_8)                          :: cpx, cpy, cpz
    INTEGER                                  :: ierr, istat, jstat
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8), ALLOCATABLE, SAVE          :: xlam(:,:)

    IF (ifirst.EQ.0) THEN
       ALLOCATE(xlam(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    ! 
    ! ==================================================================
    ! == Compute matrix elements < c_i | p | c_j > and store in pkdp. ==
    ! ==--------------------------------------------------------------==
    DO istat=1,nstate
       DO jstat=1,nstate
          CALL calc_pij(c0(1,istat,1),c0(1,jstat,1),cpx,cpy,cpz,1)
          CALL mp_sum(cpx,parai%allgrp)
          CALL mp_sum(cpy,parai%allgrp)
          CALL mp_sum(cpz,parai%allgrp)
          pkdp(1,jstat,istat)=REAL(cpx)*parm%tpiba
          pkdp(2,jstat,istat)=REAL(cpy)*parm%tpiba
          pkdp(3,jstat,istat)=REAL(cpz)*parm%tpiba
       ENDDO
    ENDDO
    ! ==================================================================
    ! == Compute the hamiltonian matrix elements at Gamma.            ==
    ! ==================================================================
    ! 
    ! test      call ehpsi(c0,ckdp,vpot,psi,nstate)
    ! test      call ovlap(nstate,xlam,c0,ckdp)
    CALL zeroing(xlambda0)!,SIZE(xlambda0))
    DO istat=1,nstate
       xlambda0(istat,istat)=CMPLX(-LOG(eigv(istat))/betap,0._real_8,kind=real_8)
       ! test        xlambda0(istat,istat)=CMPLX(-log(xlam(istat,istat))/betap,
       ! test     .                               0._real_8)
       ! deb        WRITE(6,*)istat,real(xlambda0(istat,istat))*2._real_8*ry
    ENDDO
    ! 
    RETURN
  END SUBROUTINE kdp_prep

END MODULE kdp_prep_utils
