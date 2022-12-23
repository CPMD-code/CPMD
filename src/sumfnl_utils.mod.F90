
#if defined(__SUMFNL_ISNT_TRHEADED)
#else
#define IS_THREADED 1
#endif

MODULE sumfnl_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlm,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: fnl,&
                                             fnl2,&
                                             tfnl2
  USE system,                          ONLY: ipept,&
                                             maxsys,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sumfnl
  PUBLIC :: sumfnld
  PUBLIC :: give_scr_sumfnl

CONTAINS

#if defined (_vpp_) || defined (__SR8000) || defined (__ES)
  ! ==================================================================
  SUBROUTINE sumfnl(fnl0,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8) :: fnl0(imagp,ions1%nat,maxsys%nhxs,nstate)

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, ia, iat, iat0, ierr, ii, &
                                                im, imax, is, isa, isub, iv, &
                                                length, lsumfnl
    INTEGER, ALLOCATABLE                     :: iat1(:,:,:), iii(:,:,:,:)
    REAL(real_8), ALLOCATABLE                :: scr(:,:,:,:)

    IF (parai%nproc.LE.1) RETURN
    imax=imagp*nlm*nstate
    IF (imax.EQ.0) RETURN
    CALL give_scr_sumfnl(lsumfnl,tag,nstate)

    ALLOCATE(scr(lsumfnl,1,1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL tiset('    SUMFNL',isub)
    length=imagp*ions1%nat*maxsys%nhxs*nstate
    ALLOCATE(iii(imagp,ions1%nat,maxsys%nhxs,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    length=imagp*maxsys%nax*ions1%nsp*2
    ALLOCATE(iat1(imagp*maxsys%nax,maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Fill index array used to shift the atom index
    ! I -> {1,...,NA(1),NA(1)+1,...,NA(1)+NA(2),...}
    ii=0
    DO ia=1,ions0%na(1)
       DO im=1,imagp
          ii=ii+1
          iat1(ii,1,1)=ia
          iat1(ii,1,2)=im
       ENDDO
    ENDDO
    iat0=0
    DO is=2,ions1%nsp
       ii=0
       iat0=iat0+ions0%na(is-1)
       DO ia=1,ions0%na(is)
          DO im=1,imagp
             ii=ii+1
             iat1(ii,is,1)=iat0+ia
             iat1(ii,is,2)=im
          ENDDO
       ENDDO
    ENDDO
    ! Prepare the III array for pseudo-vectorization
    ii=0
    DO i=1,nstate
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO isa=1,ions0%na(is)*imagp
                iat=iat1(isa,is,1)
                im=iat1(isa,is,2)
                ii=ii+1
                iii(im,iat,iv,i)=ii
             ENDDO
          ENDDO
       ENDDO
    ENDDO
#ifdef IS_THREADED
#ifdef __SR8000
    !poption parallel, tlocal(IAT,IM)
    !voption indep(IAT1,III,SCR)
#endif
    !$omp parallel do private(I,IS,IV,ISA,IAT,IM)
#endif
    DO i=1,nstate
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO isa=1,ions0%na(is)*imagp
                iat=iat1(isa,is,1)
                im=iat1(isa,is,2)
                scr(iii(im,iat,iv,i),1,1,1)=fnl0(im,iat,iv,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL zeroing(fnl0)!,imax)
    CALL mp_sum(scr,fnl0,imax,parai%allgrp)
    CALL dcopy(imax,fnl0(1,1,1,1),1,scr,1)
#ifdef IS_THREADED
#ifdef __SR8000
    !poption parallel, tlocal(IAT,IM)
    !voption indep(IAT1,III,SCR)
#endif
    !$omp parallel do private(I,IS,IV,ISA,IAT,IM)
#endif
    DO i=1,nstate
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO isa=1,ions0%na(is)*imagp
                iat=iat1(isa,is,1)
                im=iat1(isa,is,2)
                fnl0(im,iat,iv,i)=scr(iii(im,iat,iv,i),1,1,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(iii,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(iat1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt('    SUMFNL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sumfnl
  ! ==================================================================
#else
  ! ==================================================================
  SUBROUTINE sumfnl(fnl0,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8) :: fnl0(imagp,ions1%nat,maxsys%nhxs,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'sumfnl'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, ia, iat, iat0, ierr, ii, &
                                                imax, is, isub, iv, lsumfnl
    REAL(real_8), ALLOCATABLE                :: scr(:,:,:,:)

    IF (parai%nproc.LE.1) RETURN
    ! ==--------------------------------------------------------------==
    imax=imagp*nlm*nstate
    IF (imax.EQ.0) RETURN
    CALL give_scr_sumfnl(lsumfnl,tag,nstate)
    ALLOCATE(scr(lsumfnl,1,1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL tiset(procedureN,isub)
    ii = 0
    iat0=0
    DO is=1,ions1%nsp
       DO i=1,nstate
          DO iv=1,nlps_com%ngh(is)
             IF (imagp.EQ.2) THEN
#ifdef IS_THREADED
                !$omp parallel do private(ia,iat)
#endif
                DO ia=1,ions0%na(is)
                   iat=iat0+ia
                   scr(ii+2*ia-1,1,1,1)=fnl0(1,iat,iv,i)
                   scr(ii+2*ia,1,1,1)=fnl0(2,iat,iv,i)
                ENDDO
                ii = ii + 2*ions0%na(is)
             ELSE
#ifdef IS_THREADED
                !$omp parallel do private(ia,iat)
#endif
                DO ia=1,ions0%na(is)
                   iat=iat0+ia
                   scr(ii+ia,1,1,1)=fnl0(1,iat,iv,i)
                ENDDO
                ii = ii + ions0%na(is)
             ENDIF
          ENDDO
       ENDDO
       iat0=iat0+ions0%na(is)
    ENDDO
    CALL zeroing(fnl0)!,imax)
    CALL mp_sum(scr,fnl0,imax,parai%allgrp)
    CALL dcopy(imax,fnl0,1,scr,1)
    ii = 0
    iat0=0
    DO is=1,ions1%nsp
       DO i=1,nstate
          DO iv=1,nlps_com%ngh(is)
             IF (imagp.EQ.2) THEN
#ifdef IS_THREADED
                !$omp parallel do private(ia,iat)
#endif
                DO ia=1,ions0%na(is)
                   iat=iat0+ia
                   fnl0(1,iat,iv,i)=scr(ii+2*ia-1,1,1,1)
                   fnl0(2,iat,iv,i)=scr(ii+2*ia,1,1,1)
                ENDDO
                ii = ii + 2*ions0%na(is)
             ELSE
#ifdef IS_THREADED
                !$omp parallel do private(ia,iat)
#endif
                DO ia=1,ions0%na(is)
                   iat=iat0+ia
                   fnl0(1,iat,iv,i)=scr(ii+ia,1,1,1)
                ENDDO
                ii = ii + ions0%na(is)
             ENDIF
          ENDDO
       ENDDO
       iat0=iat0+ions0%na(is)
    ENDDO

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sumfnl
  ! ==================================================================
#endif
  ! ==================================================================
  SUBROUTINE give_scr_sumfnl(lsumfnl,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsumfnl
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    IF (parai%nproc.LE.1) THEN
       lsumfnl=0
    ELSE
       lsumfnl=imagp*nlm*nstate
       tag   ='IMAGP*NLM*NSTATE'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_sumfnl
  ! ==================================================================
  SUBROUTINE sumfnld(fnlx,isa0,nai,iv,i,ns,nstate,ikk,im,grp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isa0, nai, iv, i, ns, nstate, &
                                                ikk, im
    REAL(real_8)                             :: fnlx(im,nai,ns)
    INTEGER                                  :: grp

    CHARACTER(*), PARAMETER                  :: procedureN = 'sumfnld'

    INTEGER                                  :: ia, ii, isa, ist, isub, ix, &
                                                jj, nmax, nmin

    CALL tiset(procedureN,isub)
    CALL mp_sum(fnlx,im*nai*ns,parai%allgrp)
    nmin=MAX(isa0+1,ipept(1,parai%mepos))
    nmax=MIN(isa0+nai,ipept(2,parai%mepos))
#ifdef IS_THREADED
    !$omp parallel do private(ist,ii,isa,ia,ix)
#endif
    DO ist=i,i+ns-1
       ii=ist-i+1
       DO isa=nmin,nmax
          ia=isa-isa0
          ix=isa-ipept(1,parai%mepos)+1
          IF (im.EQ.1) THEN
             fnl(1,ix,iv,ist,ikk)=fnlx(1,ia,ii)
          ELSE
             fnl(1,ix,iv,ist,ikk)=fnlx(1,ia,ii)
             fnl(2,ix,iv,ist,ikk)=fnlx(2,ia,ii)
          ENDIF
       ENDDO
    ENDDO
    IF (tfnl2) THEN
       nmin=MAX(i,parap%nst12(parai%mepos,1))
       nmax=MIN(i+ns-1,parap%nst12(parai%mepos,2))
#ifdef IS_THREADED
       !$omp parallel do private(ist,ii,jj,ia,isa)
#endif
       DO ist=nmin,nmax
          ii=ist-parap%nst12(parai%mepos,1)+1
          jj=ist-i+1
          DO ia=1,nai
             isa=isa0+ia
             IF (im.EQ.1) THEN
                fnl2(1,isa,iv,ii,ikk)=fnlx(1,ia,jj)
             ELSE
                fnl2(1,isa,iv,ii,ikk)=fnlx(1,ia,jj)
                fnl2(2,isa,iv,ii,ikk)=fnlx(2,ia,jj)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sumfnld
  ! ==================================================================

END MODULE sumfnl_utils
