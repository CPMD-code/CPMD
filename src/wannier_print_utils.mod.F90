MODULE wannier_print_utils
  USE cppt,                            ONLY: indzs,&
                                             nzh,&
                                             nzhs
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: cnti,&
                                             fpar,&
                                             ncpw
  USE wann,                            ONLY: sw_list,&
                                             wanni,&
                                             wannl,&
                                             wannr
!!use densto_utils, only : densto
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wannier_print

CONTAINS

  ! ==================================================================
  SUBROUTINE wannier_print(infi,c0,tau0,nstate,psi,center)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: infi
    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: center(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wannier_print'

    CHARACTER(len=25)                        :: filbod, filen, filtr
    CHARACTER(len=8)                         :: numb
    COMPLEX(real_8), ALLOCATABLE             :: rg(:)
    INTEGER                                  :: i, i1, i2, ia, ie, ierr, ig, &
                                                ir, no
    INTEGER, ALLOCATABLE, SAVE               :: orbpri(:)
    INTEGER, SAVE                            :: ifirst = 0

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    IF (wannl%twpri) THEN
       IF (ifirst.EQ.0) THEN
          no=nstate
          ALLOCATE(orbpri(no),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(orbpri)!,nstate)
          IF (wanni%sw_all.NE.0) THEN
             DO i=1,nstate
                orbpri(i)=1
             ENDDO
          ELSEIF (wanni%sw_first.GT.0) THEN
             DO i=wanni%sw_first,wanni%sw_last
                orbpri(i)=1
             ENDDO
          ELSEIF (wanni%sw_orb.GT.0) THEN
             DO i=1,wanni%sw_orb
                orbpri(sw_list(i))=1
             ENDDO
          ENDIF
          ifirst=1
       ENDIF
       IF (ABS(wannr%sw_spread).GT.0.001_real_8) THEN
          CALL zeroing(orbpri)!,nstate)
          DO i=1,nstate
             IF ((wannr%sw_spread.GT.0.0_real_8).AND.(center(4,i).GT.wannr%sw_spread))&
                  orbpri(i)=1
             IF ((wannr%sw_spread.LT.0.0_real_8).AND.(center(4,i).LT.wannr%sw_spread))&
                  orbpri(i)=1
          ENDDO
          CALL mp_bcast(orbpri,SIZE(orbpri),parai%source,parai%allgrp)
       ENDIF
       IF (MOD(infi-1,wanni%ns_wann*cnti%npdip).EQ.0) THEN
          IF (paral%parent) THEN
             filtr='WANNIER_'
             IF (paral%io_parent)&
                  WRITE(numb,'(I8)') infi
             CALL xstring(numb,ia,ie)
             filbod=filtr(1:8)//numb(ia:ie)
          ENDIF

          ALLOCATE(rg(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          DO i=1,nstate
             IF (orbpri(i).NE.0) THEN
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(numb,'(I8)') i
                   CALL xstring(numb,ia,ie)
                   CALL xstring(filbod,i1,i2)
                   filen=filbod(i1:i2)//'.'//numb(ia:ie)
                ENDIF
                IF (wannl%tsdens) THEN
                   CALL zeroing(psi)!,maxfft)
                   !$omp parallel do private(IG)
                   !CDIR NODEP
                   DO ig=1,ncpw%ngw
                      psi(nzhs(ig))=c0(ig,i)
                      psi(indzs(ig))=CONJG(c0(ig,i))
                   ENDDO
                   CALL  invfftn(psi,.TRUE.,parai%allgrp)
                   !$omp parallel do private(IR)
                   DO ir=1,fpar%nnr1
                      psi(ir)=psi(ir)*psi(ir)
                   ENDDO
                   CALL  fwfftn(psi,.FALSE.,parai%allgrp)
                   !$omp parallel do private(IG)
                   DO ig=1,ncpw%nhg
                      rg(ig)=psi(nzh(ig))
                   ENDDO
                ELSE
                   CALL zeroing(rg)!,nhg)
                   CALL dcopy(2*ncpw%ngw,c0(1,i),1,rg,1)
                ENDIF
                CALL densto(rg,tau0,filen)
             ENDIF
          ENDDO
          DEALLOCATE(rg,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wannier_print
  ! ==================================================================

END MODULE wannier_print_utils
