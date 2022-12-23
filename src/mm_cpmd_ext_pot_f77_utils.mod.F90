MODULE mm_cpmd_ext_pot_f77_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_cpmd_ext_pot_f77

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_cpmd_ext_pot_f77
    ! ----------------------------------------------------------------== 
    ! This routine adds an external potential due to the MM charges
    ! on the real space grid of the QM box.
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: icoupling_model, ix, ix2, iy, &
                                                iz
    LOGICAL, SAVE                            :: first_time = .TRUE.
    REAL(real_8)                             :: dx, dy, dz, r_pot, xm, ym, zm
    REAL(real_8), POINTER, SAVE              :: extf_3d(:,:,:)

#if defined (__QMECHCOUPL)

    CALL mm_cpmd_elstat_get_setting(icoupling_model)

    IF (icoupling_model .LT.2) RETURN

    IF (first_time) THEN
       first_time = .FALSE.
       ALLOCATE(extf(kr1*kr2s*kr3s),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL reshape_inplace(extf, (/kr1, kr2s, kr3s/), extf_3d)
    ENDIF
    CALL zeroing(extf)!,kr1*kr2s*kr3s)

    ! -------------------------------------------------
    ! define the grid spacings for X,Y and Z directions
    ! -------------------------------------------------
    dx=cell_com%celldm(1)/REAL(nr1s,kind=real_8)
    dy=cell_com%celldm(2)/REAL(nr2s,kind=real_8)*cell_com%celldm(1)
    dz=cell_com%celldm(3)/REAL(nr3s,kind=real_8)*cell_com%celldm(1)

    ! ----------------------------------------------------------
    ! ZM, YM, XM are the Cartesian coordinates of the grid point
    ! in atomic units
    ! ----------------------------------------------------------

    DO iz=1,kr3s
       zm=REAL(iz-1,kind=real_8)*dz
       DO iy=1,kr2s
          ym=REAL(iy-1,kind=real_8)*dy
          DO ix=nrxpl(parai%mepos,1),nrxpl(parai%mepos,2)
             xm=REAL(ix-1,kind=real_8)*dx
             ix2 = ix-nrxpl(parai%mepos,1)+1
             CALL mm_cpmd_elstat_get_pot(xm,ym,zm,r_pot)

             extf_3d(ix2,iy,iz)=r_pot

          ENDDO
       ENDDO
    ENDDO

#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mm_cpmd_ext_pot_f77
  ! ==================================================================
  ! 

END MODULE mm_cpmd_ext_pot_f77_utils
