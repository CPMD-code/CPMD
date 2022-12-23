MODULE csmat_utils
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  USE nort,                            ONLY: nort_com
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csmat

CONTAINS

  ! ==================================================================
  ! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- CSMAT_C WILL BE WRITTEN   ==
  ! ==================================================================
  SUBROUTINE csmat(a,c0,fnl,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==         COMPUTES THE OVERLAP MATRIX A = < C0 |S| C0 >        ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8) :: fnl(ions1%nat,maxsys%nhxs,nstate,*)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: a(nstate,*)
    INTEGER                                  :: ikind

    INTEGER                                  :: i, ia, is, isa, isa0, isub, &
                                                iv, j, jmax, jv, nstx
    REAL(real_8)                             :: aij, selem

    CALL tiset('     CSMAT',isub)
    IF (cntl%tfdist) CALL stopgm('CSMAT','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    CALL ovlap(nstate,a,c0,c0)
    CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
    IF (pslo_com%tivan) THEN
       DO i=paraw%nwa12(parai%mepos,1),paraw%nwa12(parai%mepos,2)
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          DO j=i,jmax
             isa0=0
             aij=0.0_real_8
             DO is=1,ions1%nsp
                IF (pslo_com%tvan(is)) THEN
                   DO iv=1,nlps_com%ngh(is)
                      DO jv=1,nlps_com%ngh(is)
                         IF (ABS(qq(iv,jv,is)).GT.1.e-5_real_8) THEN
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               aij=aij+qq(iv,jv,is)*&
                                    fnl(isa,iv,i,ikind)*fnl(isa,jv,j,ikind)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
                isa0=isa0+ions0%na(is)
             ENDDO
             a(i,j)=a(i,j)+aij
             IF (i.NE.j) a(j,i)=a(j,i)+aij
          ENDDO
       ENDDO
    ENDIF
    CALL mp_sum(a,nstate*nstate,parai%allgrp)
    nort_com%scond=0.0_real_8
    DO i=1,nstate
       DO j=1,nstate
          selem=ABS(a(i,j))
          IF (i.EQ.j) selem=ABS(a(i,j)-1.0_real_8)
          IF (nort_com%scond.LT.selem) nort_com%scond=selem
       ENDDO
    ENDDO
    CALL tihalt('     CSMAT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csmat
  ! ==================================================================

END MODULE csmat_utils
