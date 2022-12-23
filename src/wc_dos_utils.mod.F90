MODULE wc_dos_utils
  USE cnst,                            ONLY: ry
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop1,&
                                             prop7
  USE proylm_utils,                    ONLY: proylm
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wann,                            ONLY: hmat_spread,&
                                             minspr
!!use nmr_util_p_utils, only : ffttor
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wc_dos
  !public :: sort_m
  !public :: plot_states

CONTAINS

  ! ==================================================================
  SUBROUTINE wc_dos(c0,c2,nstate,centers)
    ! ==--------------------------------------------------------------==
    ! PRINT ORBITALS
    ! SPHERICAL ANALYSIS

    ! INPUT 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: centers(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wc_dos'

    CHARACTER(len=10)                        :: name
    CHARACTER(len=30)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: c3(:,:), c4(:,:), psi(:)
    INTEGER :: i, iat, ierr, il_psi, il_rhoe, im, im_a, im_b, imin, imin_a, &
      imin_b, ind_c4, info, isp, isub, j, k, laux, nlscr, nnat, nx
    INTEGER, ALLOCATABLE                     :: ind(:)
    INTEGER, SAVE                            :: if1 = fo_verb
    LOGICAL                                  :: ferror, harm_wan
    REAL(real_8)                             :: center(3), s
    REAL(real_8), ALLOCATABLE                :: aux(:), hmat(:,:), &
                                                hmat2(:,:), newscr(:), &
                                                SPREAD(:), w(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    harm_wan = .FALSE. !vw this shall be initialized before use

    laux=12*nstate
    IF (hmat_spread)THEN
       IF (minspr.GT.0.0_real_8)THEN
          harm_wan=.FALSE.
       ELSE
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,*)' PREPARATION FOR HARMONIC ANALYSIS'
          ENDIF
          harm_wan=.TRUE.
          minspr=-minspr
       ENDIF
       IF (paral%parent)THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'ORDER HAMILTONIAN WITH SPREAD'
          IF (paral%io_parent)&
               WRITE(6,'(A,F10.5)') ' MINIMAL SPREAD =', minspr
       ENDIF
       ALLOCATE(SPREAD(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,nstate
          ! TODO check SPREAD
          SPREAD(i)=centers(4,i)
       ENDDO
       CALL mp_bcast(spread,SIZE(spread),parai%source,parai%allgrp)
    ENDIF

    ALLOCATE(hmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(12*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL ovlap(nstate,hmat,c0,c2)
    CALL mp_sum(hmat,nstate*nstate,parai%allgrp)
    IF (paral%parent) THEN
       filen='WANNIER_HAM'
       IF (tmw) THEN
          CALL mw_filename('WANNIER_HAM_',filen,mwi%walker_id)
       ENDIF
       IF (paral%io_parent)&
            CALL fileopen(40,filen,fo_app+if1,ferror)
       DO i=1,nstate
          IF (paral%io_parent)&
               WRITE(40,*)(hmat(i,j)*(-1.0_real_8*ry),j=1,nstate)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(40)
    ENDIF

    IF (hmat_spread)THEN
       ALLOCATE(ind(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hmat2(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,nstate
          ind(i)=i
       ENDDO
       IF (cntl%tlsd)THEN
          CALL sort_m(SPREAD(1),spin_mod%nsup,ind(1))
          CALL sort_m(SPREAD(1),spin_mod%nsdown,ind(spin_mod%nsup+1))
       ELSE
          CALL sort_m(SPREAD(1),nstate,ind(1))
       ENDIF
       ! ordered hamiltonian   
       DO i=1,nstate
          DO j=1,nstate
             hmat2(i,j)=hmat(ind(i),ind(j))
          ENDDO
       ENDDO
       ! check energy
       IF (paral%parent)THEN
          s=0.0_real_8
          DO i=1,nstate
             s=s+hmat(i,i)
          ENDDO
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,F14.6)')' ENERGY BEFORE',s
          ENDIF
       ENDIF

       IF (paral%parent)THEN
          filen='WANNIER_HAM_SPREAD_ORDERED'
          IF (tmw) THEN
             CALL mw_filename('WANNIER_HAM_SPREAD_ORDERED_',filen,mwi%walker_id)
          ENDIF
          IF (paral%io_parent)&
               CALL fileopen(30,filen,fo_app+if1,ferror)
          DO i=1,nstate
             IF (paral%io_parent)&
                  WRITE(30,*)(hmat2(i,j),j=1,nstate)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(30)
       ENDIF

       IF (cntl%tlsd)THEN
          imin=find_min(SPREAD(1),ind(1),spin_mod%nsup,minspr)
          im=spin_mod%nsup-imin
          imin=imin+1
          imin_a=imin
          im_a=im
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I6)')&
                  'NUMBER OF ALPHA WAVEFUNCTIONS WITH HIGH SPREAD',im_a
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I6,A,I6)') 'FROM ',imin_a,' TO ',spin_mod%nsup
          ENDIF
          IF (im.LE.0) STOP 'NO FUNCTIONS WITH SUCH SPREAD'
          ! find eigenvectors of w.f. with not large spread
          CALL dsyev('V','U',spin_mod%nsup-im,hmat2(1,1),nstate,w(1),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF ALPHA STATES WITH SMALL SPREAD FAILED ',& 
               __LINE__,__FILE__)
          ! find eigenvectors of w.f. with large spread 
          CALL dsyev('V','U',im,hmat2(imin,imin),nstate,w(imin),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF ALPHA STATES WITH LARGE SPREAD FAILED ',& 
               __LINE__,__FILE__)
          imin=find_min(SPREAD(1),ind(spin_mod%nsup+1),spin_mod%nsdown,minspr)
          im=spin_mod%nsdown-imin
          imin=imin+1+spin_mod%nsup
          imin_b=imin
          im_b=im
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I6)')&
                  'NUMBER OF BETA  WAVEFUNCTIONS WITH HIGH SPREAD',im_b
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I6,A,I6)') 'FROM ',imin_b,' TO ',nstate
          ENDIF
          IF (im.LE.0) STOP 'NO FUNCTIONS WITH SUCH SPREAD'
          ! find eigenvectors of w.f. with not large spread
          i=spin_mod%nsup+1
          CALL dsyev('V','U',spin_mod%nsdown-im,hmat2(i,i),nstate,w(i),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF BETA STATES WITH SMALL SPREAD FAILED ',& 
               __LINE__,__FILE__)
          ! find eigenvectors of w.f. with large spread 
          CALL dsyev('V','U',im,hmat2(imin,imin),nstate,w(imin),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF BETA STATES WITH LARGE SPREAD FAILED ',& 
               __LINE__,__FILE__)
       ELSE
          imin=find_min(SPREAD(1),ind(1),nstate,minspr)
          im=nstate-imin
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I6)') 'NUMBER OF WAVEFUNCTIONS  WITH HIGH'&
                  // ' SPREAD', im
          ENDIF
          imin=imin+1
          IF (im.LE.0) STOP 'NO FUNCTIONS WITH SUCH SPREAD'
          ! find eigenvectors of w.f. with not large spread
          CALL dsyev('V','U',nstate-im,hmat2(1,1),nstate,w(1),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF STATES WITH SMALL SPREAD FAILED ',& 
               __LINE__,__FILE__)
          ! find eigenvectors of w.f. with large spread 
          CALL dsyev('V','U',im,hmat2(imin,imin),nstate,w(imin),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF STATES WITH LARGE SPREAD FAILED ',& 
               __LINE__,__FILE__)
       ENDIF

       ! check energy
       IF (paral%parent)THEN
          s=0.0_real_8
          DO i=1,nstate
             s=s+w(i)
          ENDDO
          IF (paral%parent)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,F14.6)')' ENERGY AFTER',s
          ENDIF
       ENDIF

       IF (paral%parent)THEN
          filen='WANNIER_DOS_SPREAD'
          IF (tmw) THEN
             CALL mw_filename('WANNIER_DOS_SPREAD_',filen,mwi%walker_id)
          ENDIF
          IF (paral%io_parent)&
               CALL fileopen(30,filen,fo_app+if1,ferror)
          IF (paral%io_parent)&
               WRITE(30,*)(w(j),j=1,nstate)
          DO i=1,nstate
             IF (paral%io_parent)&
                  WRITE(30,*)(hmat2(i,j),j=1,nstate)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(30)
       ENDIF

       IF (paral%parent)THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')  'PLOT WAVEFUNCTIONS'
          IF (cntl%tlsd)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6)') ' NUMBER OF ALPHA STATES', im_a
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6)') ' NUMBER OF BETA  STATES',  im_b
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6)') ' NUMBER OF STATES',im
          ENDIF
       ENDIF

       ! PLOT WAVEFUNCTIONS WITH HIGH SPREAD 
       center(1) = 0._real_8
       center(2) = 0._real_8
       center(3) = 0._real_8
       nnat = 0
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

       ALLOCATE(c3(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c4(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL rhoe_psi_size(il_rhoe,il_psi)
       nlscr=2*maxfftn
       ALLOCATE(newscr(nlscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(newscr)!,nlscr)
       ALLOCATE(psi(il_psi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(psi)!,SIZE(psi))
       ! PRINT ALL STATES WITH HIGH SPREAD
       IF (cntl%tlsd)THEN
          imin=imin_a
          im=im_a
          j=imin
          DO i=1,im
             CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
             j=j+1
          ENDDO
          CALL dgemm('N','N',2*ncpw%ngw,im,im,1.0_real_8,c3,2*ncpw%ngw,&
               hmat2(imin,imin),nstate,0.0_real_8,c2,2*ncpw%ngw)
          name='WANALPHIGH'
          IF (.NOT.harm_wan)THEN
             CALL plot_states(c2,center,im,name,psi,newscr)
          ELSE
             IF (paral%parent)THEN
                IF (paral%io_parent)&
                     WRITE(6,*) 'ALPHA STATES'
             ENDIF
             ind_c4=1
             CALL dcopy(2*im*ncpw%ngw,c2(1,1),1,c4(1,ind_c4),1)
             ind_c4=ind_c4+im
          ENDIF
          imin=imin_b
          im=im_b
          j=imin
          DO i=1,im
             CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
             j=j+1
          ENDDO
          CALL dgemm('N','N',2*ncpw%ngw,im,im,1.0_real_8,c3,2*ncpw%ngw,&
               hmat2(imin,imin),nstate,0.0_real_8,c2,2*ncpw%ngw)
          name='WNBETAHIGH'
          IF (.NOT.harm_wan)THEN
             CALL plot_states(c2,center,im,name,psi,newscr)
          ELSE
             IF (paral%parent)THEN
                IF (paral%io_parent)&
                     WRITE(6,*) 'BETA STATES'
             ENDIF
             CALL dcopy(2*im*ncpw%ngw,c2(1,1),1,c4(1,ind_c4),1)
             ind_c4=ind_c4+im
          ENDIF
       ELSE
          j=imin
          DO i=1,im
             CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
             j=j+1
          ENDDO
          CALL dgemm('N','N',2*ncpw%ngw,im,im,1.0_real_8,c3,2*ncpw%ngw,&
               hmat2(imin,imin),nstate,0.0_real_8,c2,2*ncpw%ngw)
          name='WANNHIGHPS'
          IF (.NOT.harm_wan)THEN
             CALL plot_states(c2,center,im,name,psi,newscr)
          ELSE
             ind_c4=1
             CALL dcopy(2*im*ncpw%ngw,c2(1,1),1,c4(1,ind_c4),1)
             ind_c4=ind_c4+im
          ENDIF
       ENDIF
       IF (.NOT.harm_wan)THEN
          ! PRINT ALL STATES WITH LOW SPREAD
          IF (cntl%tlsd)THEN
             imin=imin_a
             im=im_a
             j=1
             k=imin-1
             DO i=1,k
                CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
                j=j+1
             ENDDO
             j=1
             CALL dgemm('N','N',2*ncpw%ngw,k,k,1.0_real_8,c3,2*ncpw%ngw,&
                  hmat2(j,j),nstate,0.0_real_8,c2,2*ncpw%ngw)
             name='WANALPHLOW'
             CALL plot_states(c2,center,k,name,psi,newscr)
             imin=imin_b
             im=im_b
             j=spin_mod%nsup+1
             k=imin-1-spin_mod%nsup
             DO i=1,k
                CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
                j=j+1
             ENDDO
             j=spin_mod%nsup+1
             CALL dgemm('N','N',2*ncpw%ngw,k,k,1.0_real_8,c3,2*ncpw%ngw,&
                  hmat2(j,j),nstate,0.0_real_8,c2,2*ncpw%ngw)
             name='WANBETALOW'
             CALL plot_states(c2,center,k,name,psi,newscr)
          ELSE
             j=1
             k=imin-1
             DO i=1,k
                CALL dcopy(2*ncpw%ngw,c0(1,ind(j)),1,c3(1,i),1)
                j=j+1
             ENDDO
             j=1
             CALL dgemm('N','N',2*ncpw%ngw,k,k,1.0_real_8,c3,2*ncpw%ngw,&
                  hmat2(j,j),nstate,0.0_real_8,c2,2*ncpw%ngw)
             name='WANNLOWPSI'
             CALL plot_states(c2,center,k,name,psi,newscr)
          ENDIF
       ENDIF
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(newscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c3,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hmat2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ind,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(spread,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' CALCULATE WANNIER FUNCTION PROJECTED DOS'
       IF (cntl%tlsd) THEN
          nx=spin_mod%nsup+1
          CALL dsyev('V','U',spin_mod%nsup,hmat(1,1),nstate,w(1),aux,laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF UP STATES FAILED ',& 
               __LINE__,__FILE__)
          CALL dsyev('V','U',spin_mod%nsdown,hmat(nx,nx),nstate,w(nx),aux,&
               laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION OF DOWN STATES FAILED ',& 
               __LINE__,__FILE__)
       ELSE
          CALL dsyev('V','U',nstate,hmat,nstate,w,aux,laux,info)
          IF (info.NE.0) CALL stopgm(' WANNIER_DOS',&
               ' DIAGONALIZATION FAILED ',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            filen='WANNIER_DOS'
       IF (tmw) THEN
          CALL mw_filename('WANNIER_DOS_',filen,mwi%walker_id)
       ENDIF
       CALL fileopen(30,filen,fo_app+if1,ferror)
       IF (paral%io_parent)&
            WRITE(30,*)(w(j),j=1,nstate)
       DO i=1,nstate
          IF (paral%io_parent)&
               WRITE(30,*)(hmat(i,j),j=1,nstate)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(30)
    ENDIF
    IF (hmat_spread.AND.harm_wan)THEN
       CALL proylm(c4,prop7%centylm,prop7%numylm,ind_c4-1,prop7%rylmax,prop7%nylmax)
       prop1%pylm=.FALSE.
       DEALLOCATE(c3,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c4,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    if1=0
    ! ==--------------------------------------------------------------==
    DEALLOCATE(hmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE wc_dos
  ! ==--------------------------------------------------------------==

  SUBROUTINE sort_m(spread,n,ind)
    INTEGER                                  :: n
    REAL(real_8)                             :: SPREAD(n)
    INTEGER                                  :: ind(n)

    INTEGER                                  :: i, im, j, temp
    REAL(real_8)                             :: minim

    DO i=1,n-1
       minim=SPREAD(ind(i))
       im=i
       DO j=i,n
          IF (SPREAD(ind(j))<minim)THEN
             minim=SPREAD(ind(j))
             im=j
          ENDIF
       ENDDO
       temp=ind(i)
       ind(i)=ind(im)
       ind(im)=temp
    ENDDO

    RETURN
  END SUBROUTINE sort_m

  SUBROUTINE plot_states(c2,center,im,tag,psi,newscr)
    USE fft_maxfft,                      ONLY: maxfftn
    COMPLEX(real_8)                          :: c2(ncpw%ngw,*)
    REAL(real_8)                             :: center(3)
    INTEGER                                  :: im
    CHARACTER(len=10)                        :: tag
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: newscr(*)

    CHARACTER(len=128)                       :: filename
    INTEGER                                  :: i

    DO i=1,im
       CALL ffttor(c2(1,i),newscr,psi,ncpw%ngw,.TRUE.)
       IF (paral%io_parent)&
            WRITE(filename,'(A10,I5.5,A5)')&
            tag,i,'.cube'
       CALL cubefile(filename,newscr,center,psi,.TRUE.)
    ENDDO

    RETURN
  END SUBROUTINE plot_states

  INTEGER FUNCTION find_min(spread,ind,n,minspr)
    IMPLICIT NONE
    ! Arguments
    INTEGER :: n, ind(n)
    REAL(real_8) :: minspr,SPREAD(n)
    ! Variables
    INTEGER :: i,imin
    imin=0
    DO i=1,n
       IF (SPREAD(ind(i))<minspr)THEN
          imin=i
       ENDIF
    ENDDO
    find_min=imin
    RETURN
  END FUNCTION find_min

END MODULE wc_dos_utils
