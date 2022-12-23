MODULE mm_cpmd_esp_charges_f77_utils
  USE cppt,                            ONLY: indz,&
                                             nr1h,&
                                             nr2h,&
                                             nr3pl,&
                                             nzh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_cpmd_esp_charges_f77

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_cpmd_esp_charges_f77 (tau0,c0,acharge,rhoe,v)
    ! This routine calculates the ESP charges from the current wavefunction
    ! 
    ! This is a butchered copy of INTERFACE_WRITE from EGO_INTERFACE
    ! 
    ! TAU0 - current atomic coordinates passed from CPMD
    ! C0 - current wavefunction expansion coefficients pass from CPMD
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: c0(*)
    REAL(real_8)                             :: acharge(ions1%nat), rhoe(*)
    COMPLEX(real_8)                          :: v(maxfftn)

    CHARACTER(*), PARAMETER :: procedureN = 'mm_cpmd_esp_charges_f77'

    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), qphi(:), &
                                                vtemp(:)
    INTEGER                                  :: i, ierr, ig, ippc, mlen
    INTEGER, ALLOCATABLE                     :: isel(:)
    REAL(real_8), ALLOCATABLE                :: efield(:)

    CALL setfftn(0)
    CALL zeroing(acharge)!,ions1%nat)
    ! 
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qphi((nr1h+1)*nr2h*nr3pl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    IF (pslo_com%tivan) CALL rnlsm(c0,crge%n,v,1,1,.FALSE.)
    CALL rhoofr(c0,rhoe,v,crge%n)
    CALL eicalc(eivps,eirop)
    !$omp parallel do private(I)
    DO i=1,fpar%nnr1
       v(i)=CMPLX(rhoe(i),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.)
    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       vtemp(ig) = v(nzh(ig))
    ENDDO
    CALL elstpo(vtemp,eirop,eivps,v,qphi,tau0)
    CALL dcopy(2*ncpw%nhg,v,1,vtemp,1)
    CALL zeroing(v)!,2*maxfft)
    DO ig=1,ncpw%nhg
       v(indz(ig)) = CONJG(vtemp(ig))
       v(nzh(ig))  = vtemp(ig)
    ENDDO
    IF (geq0.AND.isos1%tclust) THEN
       v(nzh(1)) = vtemp(1)
    ELSEIF (geq0) THEN
       v(nzh(1)) = CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDIF
    CALL  invfftn(v,.FALSE.)
    !$omp parallel do private(I)
    DO i=1,fpar%nnr1
       rhoe(i)=REAL(v(i))
    ENDDO
    mlen=fpar%nnr1/2 + 1
    ippc=0
    ALLOCATE(isel(mlen),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL selectp(isel,tau0,ippc)
    mlen=(ions1%nat+1)*ippc
    ALLOCATE(efield(mlen),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !$omp parallel do private(I)
    DO i=1,ippc
       efield(i)=rhoe(isel(i))
    ENDDO
    CALL atfield(efield,v,vtemp,eirop,qphi,isel,ippc)
    CALL espsolv(efield,acharge,ippc)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(isel,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(efield,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(qphi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mm_cpmd_esp_charges_f77
  ! ==================================================================

END MODULE mm_cpmd_esp_charges_f77_utils
