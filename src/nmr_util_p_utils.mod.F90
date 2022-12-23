#define warning_level 5.e-3_real_8

MODULE nmr_util_p_utils
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftutil_utils,                   ONLY: phase
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE machine,                         ONLY: m_walltime
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: deltarvalues,&
                                             inmr_iglo,&
                                             inmr_llc,&
                                             inmr_novirt,&
                                             inmr_wc,&
                                             nmr_options,&
                                             response1,&
                                             timetag
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: print_wavefunctions
  PUBLIC :: make_123
  PUBLIC :: printtime
  PUBLIC :: print_configuration
  PUBLIC :: epsi

CONTAINS


  ! ==================================================================
  SUBROUTINE print_wavefunctions(c0,nstate,psi)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), psi(:)

    CHARACTER(*), PARAMETER :: procedureN = 'print_wavefunctions'

    CHARACTER(len=80)                        :: filename
    INTEGER                                  :: iat, ierr, is, isp, nnat, &
                                                numstates
    REAL(real_8)                             :: center(3)
    REAL(real_8), ALLOCATABLE                :: psiofr(:)

    ALLOCATE(psiofr(fpar%nnr1+16),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    numstates = nstate
    IF (nstate .GT. 50) THEN
       IF (paral%io_parent) WRITE (6,*)&
            ' PRINT C0: TOO MANY STATES. PRINTING ONLY 1..50.'
       numstates = 50
    ENDIF
    center(1) = 0._real_8
    center(2) = 0._real_8
    center(3) = 0._real_8
    nnat=0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          center(1) = center(1) + tau0(1,iat,isp)
          center(2) = center(2) + tau0(2,iat,isp)
          center(3) = center(3) + tau0(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    center(1) =  center(1) / REAL(nnat,kind=real_8)
    center(2) =  center(2) / REAL(nnat,kind=real_8)
    center(3) =  center(3) / REAL(nnat,kind=real_8)

    DO is=1,numstates
       CALL ffttor(c0(1,is),psiofr,psi,ncpw%ngw,.TRUE.)
       IF (paral%io_parent) WRITE (filename,'(A,I2.2,A)') 'psi0-',is,'.cube '
       CALL cubefile(filename,psiofr,center,psi,.FALSE.)
    ENDDO

    DEALLOCATE(psiofr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE print_wavefunctions
  ! ==================================================================
  REAL(real_8) FUNCTION epsi(a,b,c)
    ! returns 0,+-1; the totally antisymmetric tensor epsilon
    IMPLICIT NONE
    INTEGER :: a,b,c
    epsi = 0._real_8
    IF ((b.EQ.a+1 .OR. b.EQ.a-2) .AND. (c.EQ.b+1 .OR. c.EQ.b-2))&
         epsi = 1._real_8
    IF ((b.EQ.a-1 .OR. b.EQ.a+2) .AND. (c.EQ.b-1 .OR. c.EQ.b+2))&
         epsi = -1._real_8
    RETURN
  END FUNCTION epsi
  ! ==================================================================
  SUBROUTINE make_123(first,second,third)
    ! calculates second = first + 1,
    ! \          third  = first + 2, each one modulo 3.
    INTEGER                                  :: first, second, third

    second = first + 1
    third  = first + 2
    IF (second .GT. 3) second = second - 3
    IF (third  .GT. 3) third  = third  - 3
    RETURN
  END SUBROUTINE make_123
  ! ==================================================================
  INTEGER FUNCTION complete_123(first,second)
    ! calculates third: third <> first, third <> second.
    IMPLICIT NONE
    INTEGER :: first, second, third, f2,f3
    IF (first .EQ. second)&
         CALL stopgm('complete_123','invalid operands',& 
         __LINE__,__FILE__)
    CALL make_123(first,f2,f3)
    third = 0
    IF (f2 .EQ. second) third = f3
    IF (f3 .EQ. second) third = f2
    IF (third .EQ. 0) CALL stopgm('complete_123','internal error',& 
         __LINE__,__FILE__)
    complete_123 = third
    RETURN
  END FUNCTION complete_123
  ! ==================================================================
  SUBROUTINE printtime
    REAL(real_8)                             :: delta_t
    REAL(real_8), SAVE                       :: lasttime = 0.0_real_8

    IF (.NOT. paral%parent) RETURN

    IF (lasttime .LT.1.e-10_real_8) THEN
       lasttime =m_walltime()
       RETURN
    ENDIF

    delta_t = (m_walltime()-lasttime)/1000._real_8
    lasttime =m_walltime()
    IF (paral%io_parent)&
         WRITE (6, '("*  TCPU = ",F7.1," sec: ",A40," *")')&
         delta_t,timetag

    RETURN
  END SUBROUTINE printtime
  ! ==================================================================
  SUBROUTINE print_configuration
    ! ==--------------------------------------------------------------==
    CHARACTER(len=32)                        :: text1, text2

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE (6,91)
    text1 = 'RESPONSE TYPE:'
    IF (nmr_options%tnmr_full) THEN
       text2 = 'NMR / FULL calculation'
       IF (paral%io_parent)&
            WRITE (6,92) text1,text2
    ELSE
       text2 = 'NMR chemical shifts'
       IF (paral%io_parent)&
            WRITE (6,92) text1,text2
    ENDIF
    text1 = ' '
    text2 = 'magnetic susceptibility tensors'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2

    IF (response1%t_restart_c1) THEN
       text1 = 'RESTART FULL CALCULATION'
       text2 = 'FROM PREVIOUS RUN'
       IF (paral%io_parent)&
            WRITE (6,92) text1,text2
    ENDIF
    ! ==-------------------------------------------==
    text1 = 'CALCULATION METHOD:'
    text2 = 'CSGT'
    IF (nmr_options%inmr_method .EQ. inmr_iglo) text2 = 'IGLO'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'OUTPUT:'
    text2 = 'standard'
    IF (nmr_options%tverbose) text2 = 'VERBOSE'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'ELECTRONIC CURRENT:'
    text2 = 'not computed'
    IF (nmr_options%tcurrent) text2 = 'COMPUTED -> jB_.dens)'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'ORBITAL PLOTS:'
    text2 = 'not done'
    IF (nmr_options%tprintwfns) text2 = 'WAVEFUNCTIONS'
    IF (nmr_options%tprintrho)  text2 = 'DENSITIES'
    IF (nmr_options%tprintrho.AND.nmr_options%tprintwfns)&
         text2 = 'WAVEFUNCTIONS and DENSITIES'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'FFT[R->G] - LOSSES:'
    text2 = 'warnings issued'
    IF (nmr_options%tfast) text2 = 'IGNORED'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'LOCALIZATION:'
    text2 = 'wannier functions'
    IF (.NOT. nmr_options%tlocalize) text2 = 'CANONICAL ORBITALS'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==
    text1 = 'VIRTUAL CELL MODE:'
    text2 = 'borders at minimum density'
    IF (nmr_options%inmr_virtual .EQ. inmr_wc) text2 = 'USING WANNIERCENTERS'
    IF (nmr_options%inmr_virtual .EQ. inmr_novirt) text2 = 'NO VIRTUAL CELLS'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    IF (nmr_options%inmr_virtual .EQ. inmr_llc) THEN
       text1 = 'GAUSSIAN SPREAD:'
       IF (paral%io_parent)&
            WRITE (text2,'(F10.3)') deltarvalues%sigma
       IF (paral%io_parent)&
            WRITE (6,92) text1,text2
       text1 = ' '
       IF (nmr_options%tnmr_overlaps) THEN
          text2 = 'OVERLAP OPTIMIZATION'
          IF (paral%io_parent)&
               WRITE (6,92) text1,text2
          text1 = 'OVL THREASHOLD:'
          IF (paral%io_parent)&
               WRITE (text2,'(F10.3)') deltarvalues%overlap_threashold
       ELSE
          text2 = 'each state individually'
       ENDIF
       IF (paral%io_parent)&
            WRITE (6,92) text1,text2
    ENDIF
    ! ==-------------------------------------------==
    text1 = 'SMOOTHING MODE:'
    text2 = 'at +-10% of virtual cells'
    IF (.NOT. nmr_options%tsmoothing) text2 = 'NO SMOOTHING'
    IF (paral%io_parent)&
         WRITE (6,92) text1,text2
    ! ==-------------------------------------------==

    IF (paral%io_parent)&
         WRITE (6,91)
    RETURN
    ! ==--------------------------------------------------------------==
91  FORMAT (52("*"),"NMR*RESPONSE*")
92  FORMAT ("* ",a24,1x,a30,"       *")

  END SUBROUTINE print_configuration
  ! ==================================================================


END MODULE nmr_util_p_utils

! ==================================================================
SUBROUTINE ffttog(src,dest,psi,dim,do_phase)
  ! NB: src and dest may be the same!
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:ncpw,fpar,spar
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:nzh,nzhs
  USE response_pmod , ONLY:nmr_options,response1
  USE fftutil_utils, ONLY : phase
  USE fftmain_utils, ONLY : fwfftn
  USE fft_maxfft,                      ONLY: maxfftn
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  REAL(real_8)                               :: src(*)
  COMPLEX(real_8)                            :: dest(*), psi(maxfftn)
  INTEGER                                    :: dim
  LOGICAL                                    :: do_phase

  INTEGER                                    :: ig, ir, isub
  REAL(real_8)                               :: norm1, norm2
  REAL(real_8), EXTERNAL                     :: ddot, dotp

  CALL tiset('  FFT -> G',isub)
  IF (.NOT. nmr_options%tfast) THEN
     norm1 = ddot(fpar%nnr1,src,1,src,1)/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
     CALL mp_sum(norm1,parai%allgrp)
     norm1 = SQRT(norm1)
  ENDIF

  CALL zeroing(psi)!,maxfft)
  !$omp parallel do private(ir)
  DO ir=1,fpar%nnr1
     psi(ir)  = CMPLX(src(ir),0._real_8,kind=real_8)
  ENDDO
  IF (dim .EQ. ncpw%ngw) THEN
     IF (do_phase) CALL phase(psi)
     CALL  fwfftn(psi,.TRUE.,parai%allgrp)
     !$omp parallel do private(ig)
     DO ig=1,dim
        dest(ig)  = psi(nzhs(ig))
     ENDDO
  ELSE
     CALL  fwfftn(psi,.FALSE.,parai%allgrp)
     !$omp parallel do private(ig)
     DO ig=1,dim
        dest(ig)  = psi(nzh(ig))
     ENDDO
  ENDIF

  IF (.NOT. nmr_options%tfast) THEN
     norm2 = dotp(dim, dest, dest)
     CALL mp_sum(norm2,parai%allgrp)
     norm2 = SQRT(norm2)
     IF (ABS(norm1-norm2)/norm2 .GT. warning_level) THEN
        IF ((paral%parent.AND.(.NOT.response1%tdummyatom)).AND.paral%io_parent)&
             WRITE(6, 80)&
             REAL(100._real_8,kind=real_8)*(norm1-norm2)/norm2
     ENDIF
  ENDIF

  CALL tihalt('  FFT -> G',isub)
  RETURN
80 FORMAT ("*  FFT INFORMATION LOSS WARNING: ",f10.3,&
       "% difference in norm *")
END SUBROUTINE ffttog
! ==================================================================
SUBROUTINE fft2tog(src1,dest1,src2,dest2,psi,dim,do_phase)
  ! NB: src and dest may be the same!
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:ncpw,fpar,spar
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:indz,indzs,nzh,nzhs
  USE response_pmod , ONLY:nmr_options
  USE fftutil_utils, ONLY : phase
  USE fftmain_utils, ONLY : fwfftn
  USE zeroing_utils,                   ONLY: zeroing
  USE fft_maxfft,                      ONLY: maxfftn
  IMPLICIT NONE
  REAL(real_8)                               :: src1(*), dest1(2,*), src2(*), &
                                                dest2(2,*)
  COMPLEX(real_8)                            :: psi(maxfftn)
  INTEGER                                    :: dim
  LOGICAL                                    :: do_phase

  INTEGER                                    :: ig, ir, isub
  REAL(real_8)                               :: im1, im2, norm1a, norm1b, &
                                                norm2a, norm2b, re1, re2
  REAL(real_8), EXTERNAL                     :: ddot, dotp

  CALL tiset('FFT^2 -> G',isub)
  IF (.NOT. nmr_options%tfast) THEN
     norm1a = ddot(fpar%nnr1,src1,1,src1,1)/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
     norm2a = ddot(fpar%nnr1,src2,1,src2,1)/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
     CALL mp_sum(norm1a,parai%allgrp)
     CALL mp_sum(norm2a,parai%allgrp)
     norm1a = SQRT(norm1a)
     norm2a = SQRT(norm2a)
  ENDIF


  CALL zeroing(psi)!,maxfft)
  !$omp parallel do private(ir)
  DO ir=1,fpar%nnr1-1,2
     psi(ir) = CMPLX(src1(ir),src2(ir),kind=real_8)
     psi(ir+1) = CMPLX(src1(ir+1),src2(ir+1),kind=real_8)
  ENDDO
  IF (dim .EQ. ncpw%ngw) THEN
     IF (do_phase) CALL phase(psi)
     CALL  fwfftn(psi,.TRUE.,parai%allgrp)
     !$omp parallel do private(ig,re1,re2,im1,im2)
     DO ig=1,dim
        re1 = 0.5_real_8*REAL(psi(nzhs(ig)))
        re2 = 0.5_real_8*REAL(psi(indzs(ig)))
        im1 = 0.5_real_8*AIMAG(psi(nzhs(ig)))
        im2 = 0.5_real_8*AIMAG(psi(indzs(ig)))

        dest1(1,ig) =   re1 + re2
        dest1(2,ig) =   im1 - im2
        dest2(1,ig) =   im1 + im2
        dest2(2,ig) = - re1 + re2
     ENDDO
  ELSE
     CALL  fwfftn(psi,.FALSE.,parai%allgrp)
     !$omp parallel do private(ig,re1,re2,im1,im2)
     DO ig=1,dim
        re1 = 0.5_real_8*REAL(psi(nzh(ig)))
        re2 = 0.5_real_8*REAL(psi(indz(ig)))
        im1 = 0.5_real_8*AIMAG(psi(nzh(ig)))
        im2 = 0.5_real_8*AIMAG(psi(indz(ig)))

        dest1(1,ig) =   re1 + re2
        dest1(2,ig) =   im1 - im2
        dest2(1,ig) =   im1 + im2
        dest2(2,ig) = - re1 + re2
     ENDDO
  ENDIF

  IF (.NOT. nmr_options%tfast) THEN
     norm1b = dotp(dim, dest1, dest1)
     norm2b = dotp(dim, dest2, dest2)
     CALL mp_sum(norm1b,parai%allgrp)
     CALL mp_sum(norm2b,parai%allgrp)
     norm1b = SQRT(norm1b)
     norm2b = SQRT(norm2b)
     IF (ABS(norm1a-norm1b)/norm1b .GT. warning_level) THEN
        IF (paral%io_parent)&
             WRITE (6, 80)&
             REAL(100._real_8,kind=real_8)*(norm1a-norm1b)/norm1b
     ENDIF
     IF (ABS(norm2a-norm2b)/norm2b .GT. warning_level) THEN
        IF (paral%io_parent)&
             WRITE (6, 80)&
             REAL(100._real_8,kind=real_8)*(norm2a-norm2b)/norm2b
     ENDIF
  ENDIF

  CALL tihalt('FFT^2 -> G',isub)
  RETURN

80 FORMAT ("*  FFT INFORMATION LOSS WARNING: ",f10.3,&
       "% difference in norm *")
END SUBROUTINE fft2tog

SUBROUTINE print_orbital_densities(c0,nstate,psi,aux)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:ncpw,fpar
  USE parac, ONLY : paral,parai
  USE coor , ONLY:tau0
  USE ions , ONLY:ions0,ions1
  USE response_pmod , ONLY:vofrho0
  USE fft_maxfft,                      ONLY: maxfftn
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate)
  REAL(real_8)                               :: psi(maxfftn)
  COMPLEX(real_8)                            :: aux(*)

  CHARACTER(*), PARAMETER :: procedureN = 'print_orbital_densities'

  CHARACTER(len=80)                          :: filename
  INTEGER                                    :: iat, ierr, ir, is, isp, nnat, &
                                                numstates
  REAL(real_8)                               :: center(3)
  REAL(real_8), ALLOCATABLE                  :: rhotot(:)

  ALLOCATE(rhotot(fpar%nnr1),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(rhotot)!,nnr1)

  numstates = nstate
  IF (nstate .GT. 50) THEN
     IF (paral%io_parent) THEN
        WRITE (6,*)&
             ' PRINT RHO0: TOO MANY STATES. PRINTING ONLY 1..50.'
     ENDIF
     numstates = 50
  ENDIF

  center(1) = 0._real_8
  center(2) = 0._real_8
  center(3) = 0._real_8
  nnat=0
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        center(1) = center(1) + tau0(1,iat,isp)
        center(2) = center(2) + tau0(2,iat,isp)
        center(3) = center(3) + tau0(3,iat,isp)
        nnat = nnat + 1
     ENDDO
  ENDDO
  center(1) =  center(1) / REAL(nnat,kind=real_8)
  center(2) =  center(2) / REAL(nnat,kind=real_8)
  center(3) =  center(3) / REAL(nnat,kind=real_8)



  DO is=1,numstates
     CALL ffttor(c0(1,is),psi,aux,ncpw%ngw,.FALSE.)
     !$omp parallel do private(ir)
     DO ir=1,fpar%nnr1
        psi(ir)  =  psi(ir)*psi(ir)
     ENDDO
     CALL daxpy(fpar%nnr1,1._real_8,psi,1,rhotot,1)
     IF (paral%io_parent) THEN
        IF (paral%io_parent)&
             WRITE (filename,'(A,I2.2,A)') 'rho0-',is,'.cube '
     ENDIF
     CALL cubefile(filename,psi,center,aux,.TRUE.)
  ENDDO
  IF (paral%io_parent) WRITE (filename,'(A)') 'rho0-all.cube '
  CALL cubefile(filename,rhotot,center,aux,.TRUE.)
  IF (paral%io_parent) WRITE (filename,'(A)') 'vofrho0.cube '
  CALL cubefile(filename,vofrho0,center,aux,.TRUE.)
  DEALLOCATE(rhotot,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  RETURN
END SUBROUTINE print_orbital_densities
! ==================================================================

! ==================================================================
SUBROUTINE fft2tor(src1,dest1,src2,dest2,psi,dim,do_phase)
  ! dim(src) = ngw (cmplx), dim(dest) = nnr1 (real)
  ! NB: src and dest may be the same!
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:ncpw,fpar
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:indz,indzs,nzh,nzhs
  USE fftutil_utils, ONLY : phase
  USE fftmain_utils, ONLY : invfftn
  USE zeroing_utils,                   ONLY: zeroing
  USE fft_maxfft,                      ONLY: maxfftn
  IMPLICIT NONE
  REAL(real_8)                               :: src1(2,*), dest1(*), &
                                                src2(2,*), dest2(*)
  COMPLEX(real_8)                            :: psi(maxfftn)
  INTEGER                                    :: dim
  LOGICAL                                    :: do_phase

  INTEGER                                    :: ig, ir, isub
  REAL(real_8)                               :: im1, im2, re1, re2

  CALL tiset('FFT^2 -> R',isub)
  CALL zeroing(psi)!,maxfft)
  IF (dim .EQ. ncpw%ngw) THEN
     !CDIR NODEP
#ifdef _vpp_
     !OCL NOVREC
#endif 
     DO ig=1,dim
        re1               =   src1(1,ig)
        re2               =   src2(1,ig)
        im1               =   src1(2,ig)
        im2               =   src2(2,ig)
        psi(nzhs(ig))   =   CMPLX(re1 - im2, im1 + re2,kind=real_8)
        psi(indzs(ig))  =   CMPLX(re1 + im2, - im1 + re2,kind=real_8)
     ENDDO
     CALL  invfftn(psi,.TRUE.,parai%allgrp)
     IF (do_phase) CALL phase(psi)
  ELSEIF (dim .EQ. ncpw%nhg) THEN
     !CDIR NODEP
#ifdef _vpp_
     !OCL NOVREC
#endif 
     DO ig=1,dim
        re1              =   src1(1,ig)
        re2              =   src2(1,ig)
        im1              =   src1(2,ig)
        im2              =   src2(2,ig)
        psi(nzh(ig))   =   CMPLX(re1 - im2, im1 + re2,kind=real_8)
        psi(indz(ig))  =   CMPLX(re1 + im2, - im1 + re2,kind=real_8)
     ENDDO
     CALL invfftn(psi,.FALSE.,parai%allgrp)
  ELSE
     CALL stopgm('fft2TOr    ','wrong dimension.',& 
          __LINE__,__FILE__)
  ENDIF

  !$omp parallel do private(ir)
  DO ir=1,fpar%nnr1-1,2
     dest1(ir)  = REAL(psi(ir))
     dest2(ir)  = AIMAG(psi(ir))
     dest1(ir+1)  = REAL(psi(ir+1))
     dest2(ir+1)  = AIMAG(psi(ir+1))
  ENDDO

  CALL tihalt('FFT^2 -> R',isub)
  RETURN
END SUBROUTINE fft2tor

! ==================================================================
SUBROUTINE apply_op_p(src, dest, coordinate, dim)
  ! PARALLEL
  ! =----------------------------------------------------------------=
  ! INPUT:  src:   wavefunctions in G-space
  ! \       coord: whether p_x (coord=1) or p_y/z (coord=2/3) should
  ! \              be applied.
  ! OUTPUT: dest:  the resulting wavefunction in G-space
  ! \              = d/dx or d/dy or d/dz   src
  ! \              previous content of dest is overwritten.
  ! \       src and dest MAY BE IDENTICAL.
  ! =----------------------------------------------------------------=
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:parm
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:gk
  IMPLICIT NONE
  COMPLEX(real_8)                            :: src(*), dest(*)
  INTEGER                                    :: coordinate, dim

  INTEGER                                    :: ig, isub
  REAL(real_8)                               :: coeff, src_IM, src_RE

! =----------------------------------------------------------------=
! Arguments:
! Local variables:
! =----------------------------------------------------------------=

  CALL tiset('     p_psi',isub)
  !$omp parallel do private(ig,src_RE,src_IM,coeff)
  DO ig=1,dim
     src_RE     =   REAL(src(ig),KIND=real_8)
     src_IM     =   AIMAG(src(ig))
     coeff      =   parm%tpiba  * gk(coordinate,ig)
     dest(ig)   = CMPLX(- coeff  * src_IM, coeff  * src_RE, KIND=real_8 )
  ENDDO

  CALL tihalt('     p_psi',isub)
  RETURN
END SUBROUTINE apply_op_p

! ==================================================================
SUBROUTINE consolidate(src1, src2, dest)
  ! PARALLEL
  ! =----------------------------------------------------------------=
  ! INPUT:  src1,2:   wavefunctions in R-space
  ! OUTPUT: dest(r) = src1(r)  * src2(r)
  ! src1,src2,dest MAY BE IDENTICAL.
  ! =----------------------------------------------------------------=
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: fpar
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  REAL(real_8)                               :: src1(fpar%nnr1), &
                                                src2(fpar%nnr1), &
                                                dest(fpar%nnr1)

  INTEGER                                    :: ir
  REAL(real_8)                               :: x1, x2

! =----------------------------------------------------------------=
! Arguments:
! Local variables:
! =----------------------------------------------------------------=

  !$omp parallel do private(ir,x1,x2)
  DO ir=1,fpar%nnr1-1,2
     x1 =  src1(ir)   * src2(ir)
     x2 =  src1(ir+1) * src2(ir+1)
     dest(ir)   = x1
     dest(ir+1) = x2
  ENDDO
  RETURN
END SUBROUTINE consolidate
! ==================================================================
SUBROUTINE ffttor(src,dest,psi,dim,do_phase)
  ! dim(src) = ngw or nhg (cmplx), dim(dest) = nnr1 (real)
  ! NB: src and dest may be the same!
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:ncpw,fpar
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:indz,indzs,nzh,nzhs
  USE fftutil_utils, ONLY : phase
  USE fftmain_utils, ONLY : invfftn
  USE zeroing_utils,                   ONLY: zeroing
  USE fft_maxfft,                      ONLY: maxfftn
  IMPLICIT NONE
  COMPLEX(real_8)                            :: src(*)
  REAL(real_8)                               :: dest(*)
  COMPLEX(real_8)                            :: psi(maxfftn)
  INTEGER                                    :: dim
  LOGICAL                                    :: do_phase

  INTEGER                                    :: ig, ir, isub

! =----------------------------------------------------------------=

  CALL tiset('  FFT -> R',isub)
  CALL zeroing(psi)!,maxfft)
  IF (dim .EQ. ncpw%ngw) THEN
     !CDIR NODEP
#ifdef _vpp_
     !OCL NOVREC
#endif 
     DO ig=1,dim
        psi(nzhs(ig))  = src(ig)
        psi(indzs(ig)) = CONJG(src(ig))
     ENDDO
     CALL  invfftn(psi,.TRUE.,parai%allgrp)
     IF (do_phase) CALL phase(psi)
  ELSEIF (dim .EQ. ncpw%nhg) THEN
     !CDIR NODEP
#ifdef _vpp_
     !OCL NOVREC
#endif 
     DO ig=1,dim
        psi(nzh(ig))  = src(ig)
        psi(indz(ig)) = CONJG(src(ig))
     ENDDO
     CALL invfftn(psi,.FALSE.,parai%allgrp)
  ELSE
     CALL stopgm('fftTOr   ','wrong dimension.',& 
          __LINE__,__FILE__)
  ENDIF
  !$omp parallel do private(ir)
  DO ir=1,fpar%nnr1
     dest(ir)  = REAL(psi(ir))
  ENDDO
  CALL tihalt('  FFT -> R',isub)
  RETURN
END SUBROUTINE ffttor
