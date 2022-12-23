C     ==================================================================
      SUBROUTINE mm_qm_boundary(TAUP,VELP) 
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
C     == Imposing (soft) reflecting walls to the QM system to confine ==
C     == the QM atoms always inside the QM box without re-centering   ==
C     ==--------------------------------------------------------------==
      use system
      use ions
      use mm_dimmod
      use isos
      use cell
      use mm_input
      use comvel_utils, only : comvel
      IMPLICIT NONE 
C     Arguments 
      REAL*8  TAUP(3,maxsys%nax,ions1%nsp),VELP(3,maxsys%nax,ions1%nsp)
      REAL*8  VCMIO(4)
C     Variables 
      INTEGER IS,IA
C     ==--------------------------------------------------------------==
      IF(.NOT.isos1%TCLUST) RETURN 
C     Modify velocities of QM ions near to the QM boundary 
      DO IS=1,mmdim%NSPq
        DO IA=1,NAq(IS)
          IF(TAUP(1,IA,IS).LT.rqmmm%box_wall) 
     .                        VELP(1,IA,IS) =  DABS(VELP(1,IA,IS))
          IF((CELL_COM%CELLDM(1)-TAUP(1,IA,IS)).LT.rqmmm%box_wall) 
     .                        VELP(1,IA,IS) = -DABS(VELP(1,IA,IS))
          IF(TAUP(2,IA,IS).LT.rqmmm%box_wall) 
     .                        VELP(2,IA,IS) =  DABS(VELP(2,IA,IS))
          IF((CELL_COM%CELLDM(2)*CELL_COM%CELLDM(1)-TAUP(2,IA,IS))
     &         .LT.rqmmm%box_wall)
     .                        VELP(2,IA,IS) = -DABS(VELP(2,IA,IS))
          IF(TAUP(3,IA,IS).LT.rqmmm%box_wall) 
     .                        VELP(3,IA,IS) =  DABS(VELP(3,IA,IS))
          IF((CELL_COM%CELLDM(3)*CELL_COM%CELLDM(1)-TAUP(3,IA,IS))
     &         .LT.rqmmm%box_wall) 
     .                        VELP(3,IA,IS) = -DABS(VELP(3,IA,IS))
        ENDDO 
      ENDDO 
C     Subtract center of mass velocity of the whole system 
C     after each walls reflection
CMcB  according to the logic of comvel(), this should do the job
      CALL COMVEL(VELP,VCMIO,.TRUE.)
      CALL COMVEL(VELP,VCMIO,.FALSE.)
C     ==--------------------------------------------------------------==
      RETURN 
      END 
C     ==================================================================
