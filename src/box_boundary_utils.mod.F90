MODULE box_boundary_utils
  USE cell,                            ONLY: cell_com
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: box_boundary

CONTAINS

  ! ==================================================================
  SUBROUTINE box_boundary(taup,velp)
    ! ==--------------------------------------------------------------==
    ! == Imposing reflecting walls to a system free from periodic     ==
    ! == boundary conditions to keep atoms inside the box             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: alay, alaz, vcmio(4)

    IF (.NOT.isos1%tclust) RETURN
    ! Modify velocities of the ions near to the boundary 
    IF (isos1%ttwod) THEN ! two-dimensional PBC
       ! We assume orthorhobmic box with surface in xy plane
       alaz=cell_com%celldm(3)*cell_com%celldm(1)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (taup(3,ia,is).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) =  ABS(VELP(3,IA,IS))
             IF ((alaz-taup(3,ia,is)).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) = -ABS(VELP(3,IA,IS))
          ENDDO
       ENDDO
    ELSE IF (isos1%toned) THEN  ! one-dimensional PBC
       ! We assume periodicity along x, orthorhobmic box in y and z
       alaz=cell_com%celldm(3)*cell_com%celldm(1)
       alay=cell_com%celldm(2)*cell_com%celldm(1)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (taup(2,ia,is).LT.isos2%wall_skin)&
                  VELP(2,IA,IS) =  ABS(VELP(2,IA,IS))
             IF ((alay-taup(2,ia,is)).LT.isos2%wall_skin)&
                  VELP(2,IA,IS) = -ABS(VELP(2,IA,IS))
             IF (taup(3,ia,is).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) =  ABS(VELP(3,IA,IS))
             IF ((alaz-taup(3,ia,is)).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) = -ABS(VELP(3,IA,IS))
          ENDDO
       ENDDO
    ELSE ! no PBC along any direction
       alaz=cell_com%celldm(3)*cell_com%celldm(1)
       alay=cell_com%celldm(2)*cell_com%celldm(1)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (taup(1,ia,is).LT.isos2%wall_skin)&
                  VELP(1,IA,IS) =  ABS(VELP(1,IA,IS))
             IF ((cell_com%celldm(1)-taup(1,ia,is)).LT.isos2%wall_skin)&
                  VELP(1,IA,IS) = -ABS(VELP(1,IA,IS))
             IF (taup(2,ia,is).LT.isos2%wall_skin)&
                  VELP(2,IA,IS) =  ABS(VELP(2,IA,IS))
             IF ((alay-taup(2,ia,is)).LT.isos2%wall_skin)&
                  VELP(2,IA,IS) = -ABS(VELP(2,IA,IS))
             IF (taup(3,ia,is).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) =  ABS(VELP(3,IA,IS))
             IF ((alaz-taup(3,ia,is)).LT.isos2%wall_skin)&
                  VELP(3,IA,IS) = -ABS(VELP(3,IA,IS))
          ENDDO
       ENDDO
    ENDIF
    ! Subtract center of mass velocity of the whole system 
    ! after each walls reflection
    IF (comvl%tsubcom) THEN
       CALL comvel(velp,vcmio,.TRUE.)
       CALL comvel(velp,vcmio,.FALSE.)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE box_boundary
  ! ==================================================================

END MODULE box_boundary_utils
