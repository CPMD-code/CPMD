MODULE hip_utils
  USE cppt,                            ONLY: hgpot,&
                                             nr1h,&
                                             nr2h,&
                                             nr3h,&
                                             nr3pl
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: nr1m,&
                                             nr3m,&
                                             xf,&
                                             yf
  USE fftmain_utils,                   ONLY: mltfft
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE system,                          ONLY: fpar,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hip
  PUBLIC :: give_qphi

CONTAINS

  ! ==================================================================
  SUBROUTINE hip(v_1d,qphi)
    ! ==--------------------------------------------------------------==
    ! == USED IN VOFRHOA IF TCLUST=.TRUE.                             ==
    ! == HOCKNEY INTERACTION POTENTIAL                                ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: v_1d(:)
    COMPLEX(real_8)                          :: qphi(:)

    CHARACTER(len=*), PARAMETER              :: procedureN = 'hip'

    COMPLEX(real_8), ALLOCATABLE             :: qphia(:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: va_1d
    COMPLEX(real_8), POINTER                 :: v(:,:,:), va(:,:,:)
    INTEGER                                  :: blklen, i, ierr, iqph, isign, &
                                                isub, ix, ixyz, iy, iz, izp, &
                                                kx, ky, kz, nr1hh, nr2hh, &
                                                nr3hh, nx, nxy
    REAL(real_8)                             :: dvol, scale

!v(kr1,kr2s,*)
! Variables
! cmb
! ==--------------------------------------------------------------==

    CALL tiset('       HIP',isub)
    v(1:fpar%kr1,1:fpar%kr2s,1:SIZE(v_1d)/(fpar%kr1*fpar%kr2s)) => v_1d
    !CALL reshape_inplace(v_, (/fpar%kr1, fpar%kr2s, SIZE(v_)/(fpar%kr1*fpar%kr2s) /), v)!vw hacking the hack
    ALLOCATE(va_1d(fpar%kr1*fpar%kr2s*nr3h),STAT=ierr)
    !ALLOCATE(va(fpar%kr1,fpar%kr2s,nr3h),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    va(1:fpar%kr1,1:fpar%kr2s,1:nr3h) => va_1d
    ALLOCATE(qphia((nr1h+1)*nr2h*nr3pl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(va)
    CALL zeroing(qphia)

    nr1hh=nr1h/2
    nr2hh=nr2h/2
    nr3hh=nr3h/2
    dvol = REAL((nr1h/spar%nr1s) * (nr2h/spar%nr2s) * (nr3h/spar%nr3s),kind=real_8) * parm%omega
    CALL zeroing(v(:,:,spar%nr3s+1:spar%nr3s+fpar%nnr1/(fpar%kr1*fpar%kr2s)))!,nnr1)
    CALL zeroing(v(:,:,1))!,kr1*kr2s)
    IF (isos1%ttwod) THEN
    ELSEIF (isos1%toned) THEN
       !$omp parallel do private(IZ,IX) shared(V)
       DO iz=1,spar%nr3s
          DO ix=1,fpar%kr1
             v(ix,1,iz)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IZ,IX,IY) shared(V)
       DO iz=1,spar%nr3s
          DO ix=1,fpar%kr1
             v(ix,1,iz)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
          IF (parap%nrxpl(parai%mepos,1).EQ.1) THEN
             DO iy=1,fpar%kr2s
                v(1,iy,iz)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    CALL dcopy(2*fpar%kr1*fpar%kr2s*nr3h,v(1,1,1),1,va(1,1,1),1)
    ! ==--------------------------------------------------------------==
    ! ==  SOLVE THE POISSON EQUATION                                  ==
    ! ==--------------------------------------------------------------==
    ! FFT ALONG Z
    isign=1
    scale=1._real_8/REAL(nr1h*nr2h*nr3h,kind=real_8)
    CALL mltfft('T','T',va_1d,fpar%kr1*fpar%kr2s,nr3h,v_1d,fpar%kr1*fpar%kr2s,nr3h,&
         nr3h,fpar%kr1*fpar%kr2s,isign,scale)
    ! Pack data for transpose
    CALL zeroing(qphi)!,(nr1h+1)*nr2h*nr3pl)
    nx=parm%nr1
    nxy=nx*spar%nr2s
    !$omp parallel do private(I,IZ,IY,IX,IXYZ)
    DO i=0,parai%nproc-1
#ifdef _vpp_
       !OCL NOALIAS
#endif
       DO iz=parap%nrzpl(i,1),parap%nrzpl(i,2)
          DO iy=1,spar%nr2s
             DO ix=1,parm%nr1
                ixyz=ix+(iy-1)*nx+(iz-parap%nrzpl(i,1))*nxy
                xf(i*nr1m*spar%nr2s*nr3m+ixyz,1)=v(ix,iy,iz)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! All to all communication
    blklen=16*nr1m*spar%nr2s*nr3m
    CALL my_trans(xf,yf,blklen,1)
    ! Unpack data
    !$omp parallel do private(I,IZ,IY,IX,IXYZ,IQPH,nx,nxy)
    DO i=0,parai%nproc-1
       nx=1+parap%nrxpl(i,2)-parap%nrxpl(i,1)
       nxy=nx*spar%nr2s
       DO iz=1,nr3pl
          DO iy=1,spar%nr2s
             DO ix=parap%nrxpl(i,1),parap%nrxpl(i,2)
                ixyz=1+ix-parap%nrxpl(i,1)+(iy-1)*nx+(iz-1)*nxy
                iqph=ix+(iy-1)*(nr1h+1)+(iz-1)*(nr1h+1)*nr2h
                qphi(iqph)=yf(i*nr1m*spar%nr2s*nr3m+ixyz,1)! YF(IXYZ,I+1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL dcopy(2*(nr1h+1)*nr2h*nr3pl,qphi(1),1,qphia(1),1)
    DO iz=1,nr3pl
       izp=(iz-1)*(nr1h+1)*nr2h + 1
       ! FFT ALONG Y
       isign=1
       scale=1._real_8
       CALL mltfft('T','T',qphia(izp:),nr1h+1,nr2h,qphi(izp:),&
            nr1h+1,nr2h,nr2h,spar%nr1s,isign,scale)
       CALL dcopy(2*(nr1h+1)*nr2h,qphi(izp),1,qphia(izp),1)
       ! FFT ALONG X
       isign=1
       scale=1._real_8
       CALL mltfft('N','N',qphia(izp:),nr1h+1,nr2h,qphi(izp:),&
            nr1h+1,nr2h,nr1h,nr2h,isign,scale)
       ! MULTIPLY BY INFLUENCE FUNCTION
       kz=iz-1
       IF (kz.GT.nr3hh) kz=nr3h-kz
       !$omp parallel do private(IY,KY,IX,IXYZ,KX)
       DO iy=1,nr2h
          ky=iy-1
          IF (ky.GT.spar%nr2s) ky=nr2h-ky! important: not NR2HH
          DO ix=1,nr1h
             ixyz=izp+(iy-1)*(nr1h+1)+ix-1
             kx=ix-1
             IF (kx.GT.spar%nr1s) kx=nr1h-kx! important: not NR1HH
             qphi(ixyz)=qphi(ixyz)*hgpot(kx+1,ky+1,kz+1)*dvol
          ENDDO
       ENDDO
       CALL dcopy(2*(nr1h+1)*nr2h,qphi(izp),1,qphia(izp),1)
       ! INVFFT ALONG X
       isign=-1
       scale=1._real_8
       CALL mltfft('N','N',qphia(izp:),nr1h+1,nr2h,qphi(izp:),&
            nr1h+1,nr2h,nr1h,nr2h,isign,scale)
       ! INVFFT ALONG Y
       isign=-1
       scale=1._real_8
       CALL mltfft('T','T',qphi(izp:),nr1h+1,nr2h,qphia(izp:),&
            nr1h+1,nr2h,nr2h,spar%nr1s,isign,scale)
    ENDDO
    CALL dcopy(2*(nr1h+1)*nr2h*nr3pl,qphia(1),1,qphi(1),1)
    ! Pack data for transpose
    !$omp parallel do private(I,IZ,IY,IX,IXYZ,IQPH,NX,NXY)
    DO i=0,parai%nproc-1
       nx=1+parap%nrxpl(i,2)-parap%nrxpl(i,1)
       nxy=nx*spar%nr2s
#ifdef _vpp_
       !OCL NOALIAS
#endif
       DO iz=1,nr3pl
          DO iy=1,spar%nr2s
             DO ix=parap%nrxpl(i,1),parap%nrxpl(i,2)
                ixyz=1+ix-parap%nrxpl(i,1)+(iy-1)*nx+(iz-1)*nxy
                iqph=ix+(iy-1)*(nr1h+1)+(iz-1)*(nr1h+1)*nr2h
                yf(i*nr1m*spar%nr2s*nr3m+ixyz,1)=qphi(iqph)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! All to all communication
    blklen=16*nr1m*spar%nr2s*nr3m
    CALL my_trans(yf,xf,blklen,1)
    ! Unpack data
    CALL zeroing(v)!,kr1*kr2s*nr3h)
    CALL zeroing(va)!,kr1*kr2s*nr3h)
    nx=parm%nr1
    nxy=nx*spar%nr2s
    !$omp parallel do private(I,IZ,IY,IX,IXYZ)
    DO i=0,parai%nproc-1
       DO iz=parap%nrzpl(i,1),parap%nrzpl(i,2)
          DO iy=1,spar%nr2s
             DO ix=1,parm%nr1
                ixyz=ix+(iy-1)*nx+(iz-parap%nrzpl(i,1))*nxy
                va(ix,iy,iz)=xf(i*nr1m*spar%nr2s*nr3m+ixyz,1)! XF(IXYZ,I+1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! INVFFT ALONG Z
    isign=-1
    scale=1._real_8
    CALL mltfft('T','T',va_1d,fpar%kr1*fpar%kr2s,nr3h,v_1d,fpar%kr1*fpar%kr2s,nr3h,nr3h,&
         fpar%kr1*fpar%kr2s,isign,scale)
    DEALLOCATE(va_1d,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qphia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('       HIP',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hip
  ! ==================================================================
  SUBROUTINE give_qphi(il_qphi)
    ! ==--------------------------------------------------------------==
    ! == GIVE THE SIZE OF QPHI                                        ==
    ! == ARRAY used in HIP                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: il_qphi

    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          il_qphi   = 2*(2*spar%nr1s+1)*2*fpar%kr2s*nr3pl
       ELSEIF (isos3%ps_type.EQ.2) THEN
          il_qphi   = 2
       ELSEIF (isos3%ps_type.EQ.3) THEN
          il_qphi   = 2
       ENDIF
    ELSE
       il_qphi = 0
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_qphi
  ! ==================================================================

END MODULE hip_utils
