MODULE kdp_rho_utils
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE kinds,                           ONLY: real_8
  USE rhoofr_kdp_utils,                ONLY: rhoofr_kdp
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             group,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kdp_rho

CONTAINS

  SUBROUTINE kdp_rho(nkdp,c0,ckdp,rhoe,psi,nstate,akdp,&
       noccup,bmix,fkdp,wkdp,bkdp)
    ! 
    INTEGER                                  :: nkdp
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(maxfft*group%nogrp)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: ckdp(nkpt%ngwk,nstate), &
      c0(nkpt%ngwk,nstate,nkpt%nkpnt), akdp(nstate,nstate,nkdp)
    INTEGER                                  :: noccup
    REAL(real_8)                             :: bmix, fkdp(nstate,nkdp), &
                                                wkdp(nkdp)
    COMPLEX(real_8)                          :: bkdp(nstate,nstate,nkdp)

    CHARACTER(*), PARAMETER                  :: procedureN = 'kdp_rho'

    COMPLEX(real_8), ALLOCATABLE             :: brkdp(:,:)
    INTEGER                                  :: ierr, ikdp, istat, isub, &
                                                jstat, lstat

! Variables
! 

    CALL tiset('   KDP_RHO',isub)
    ! 
    ! ==================================================================
    ! == Construct matrix bkdp (b_jl(k) in the paper):                ==
    ! ==                                                              ==
    ! ==   bkdp(i,j) = sum (fkdp(l,k)* akdp^*(l,i,k)*akdp(l,j,k))     ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    ! Allocation of brkdp
    ALLOCATE(brkdp(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(bkdp)!,SIZE(bkdp))
    DO ikdp=1,nkdp
       DO istat=1,nstate
          DO jstat=1,nstate
             DO lstat=1,noccup
                bkdp(istat,jstat,ikdp)=bkdp(istat,jstat,ikdp)+fkdp(lstat,&
                     ikdp)*CONJG(akdp(lstat,istat,ikdp))*akdp(lstat,&
                     jstat,ikdp)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! 
    ! ==================================================================
    ! == Construct the coefficients brkdp of the hybrid wavefunction  ==
    ! == to be used in the construction of the density in real space  ==
    ! == (beta_jl in the paper).                                      ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(brkdp)!,SIZE(brkdp))
    DO istat=1,nstate
       DO jstat=1,nstate
          DO ikdp=1,nkdp
             brkdp(istat,jstat)=brkdp(istat,jstat)+0.5_real_8*wkdp(ikdp)*&
                  CMPLX(REAL(bkdp(istat,jstat,ikdp)+bkdp(jstat,istat,&
                  ikdp)), 0._real_8,kind=real_8)
          ENDDO
       ENDDO
       ! deb        do jstat=1,nstate
       ! deb          print *,istat,jstat,brkdp(istat,jstat)
       ! deb        enddo
    ENDDO
    ! deb      print *
    ! 
    ! ==================================================================
    ! == construct the hybrid wavefunction (stored in ckdp)           ==
    ! == ckdp = brkdp . c0 ( = Sum beta_jl . phi_l in the paper)      ==
    ! ==--------------------------------------------------------------==
    ! ibm        call zgemul(c0,ngwk,'N',brkdp,nstate,'T',ckdp,ngwk,ngw,nstate,
    ! ibm     .              nstate)
    CALL zgemm('N','T',ncpw%ngw,nstate,nstate,CMPLX(1._real_8,0._real_8,kind=real_8),c0,nkpt%ngwk,&
         brkdp,nstate,CMPLX(0._real_8,0._real_8,kind=real_8),ckdp,nkpt%ngwk)
    ! 
    CALL tihalt('   KDP_RHO',isub)
    ! 
    ! ==================================================================
    ! == compute the k-point weighted electronic density              ==
    ! == rho(r) = Sum c0^*_j(r) . ckdp_j(r)                           ==
    ! ==--------------------------------------------------------------==
    CALL rhoofr_kdp(c0,ckdp,rhoe,psi,nstate)
    ! 
    ! ==--------------------------------------------------------------==
    ! Deallocation of brkdp
    DEALLOCATE(brkdp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kdp_rho
  ! ==================================================================

END MODULE kdp_rho_utils
