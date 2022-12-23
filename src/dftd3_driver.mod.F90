MODULE dftd3_driver
  USE dftd3_api,                       ONLY: dftd3_dispersion,&
                                             dftd3_pbc_dispersion
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE strs,                            ONLY: alpha,&
                                             beta
  USE system,                          ONLY: maxsp,&
                                             parm
  USE vdwcmod,                         ONLY: dftd3c,&
                                             tdftd3noc
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dftd3

CONTAINS

  SUBROUTINE dftd3(tau0,nx,ny,nz,evdw,fion,devdw)

    REAL(real_8), INTENT(in)                 :: tau0(:,:,:)
    INTEGER, INTENT(in)                      :: nx,ny,nz
    REAL(real_8), INTENT(inout)              :: evdw
    REAL(real_8), INTENT(inout)              :: fion(:,:,:)
    REAL(real_8), INTENT(inout)              :: devdw(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dftd3'

    INTEGER                                  :: ierr,isa,ia,is,i,natx,&
                                                n_rep(3)
    INTEGER, ALLOCATABLE                     :: atnum(:),&
                                                iatptx(:,:)
    REAL(real_8), ALLOCATABLE                :: coords(:,:),&
                                                grads(:,:)
    REAL(real_8)                             :: edisp,&
                                                latvecs(3,3),&
                                                stress(3,3),&
                                                x,y,z,x_,y_,z_

    evdw=0.0_real_8
    CALL zeroing(devdw)

    ! Build list vector iatptx
    ALLOCATE(iatptx(2,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
    natx=0
    DO is=1,ions1%nsp
       IF (tdftd3noc(is)) CYCLE
       DO ia=1,ions0%na(is)
          natx=natx+1
          iatptx(1,natx)=ia
          iatptx(2,natx)=is
       ENDDO
    ENDDO 

    ALLOCATE(coords(3,natx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
    ALLOCATE(atnum(natx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
    ALLOCATE(grads(3,natx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)

    DO isa=1,natx
       ia=iatptx(1,isa)
       is=iatptx(2,isa)
       coords(:,isa)=tau0(:,ia,is)
       atnum(isa)=ions0%iatyp(is)
    ENDDO
    latvecs=TRANSPOSE(metr_com%ht)

    IF (parm%ibrav==0) THEN
       ! non-periodic case
       CALL dftd3_dispersion(dftd3c, coords, atnum, edisp, grads)
    ELSE
       DO isa=1,natx
          x_=coords(1,isa); y_=coords(2,isa); z_=coords(3,isa)
          CALL pbc(x_,y_,z_,x,y,z,1,parm%apbc,parm%ibrav)
          coords(1,isa)=x; coords(2,isa)=y; coords(3,isa)=z
       ENDDO
       n_rep(1)=nx; n_rep(2)=ny; n_rep(3)=nz
       ! periodic case
       CALL dftd3_pbc_dispersion(dftd3c, coords, atnum, latvecs, n_rep, edisp, grads, stress)
    ENDIF
    ! Energy
    evdw=edisp
    ! Ionic forces
    DO isa=1,natx
       ia=iatptx(1,isa)
       is=iatptx(2,isa)
       fion(:,ia,is)=fion(:,ia,is)-grads(:,isa)
    ENDDO
    ! Stress tensor
    IF (.NOT.(ANY(alpha==0).OR.ANY(beta==0))) THEN
       DO i=1,6
          devdw(i)=stress(alpha(i),beta(i))
       ENDDO
    ENDIF

    DEALLOCATE(iatptx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(coords,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(atnum,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(grads,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

  END SUBROUTINE dftd3

END MODULE dftd3_driver
