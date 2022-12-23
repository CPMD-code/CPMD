MODULE up3_p_utils
  USE cnst,                            ONLY: uimag
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gs_ortho_c
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE rkpnt_utils,                     ONLY: rkpnt
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE setirec_utils,                   ONLY: write_irec
  USE sphe,                            ONLY: tsphere
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE utils,                           ONLY: zclean_k
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: up3_p
  PUBLIC :: give_scr_up3_p

CONTAINS

  ! ==================================================================
  SUBROUTINE up3_p(c0,cu1,tau0,psi,nstate,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==          Computes:                                           ==
    ! ==          New Set of basis wafeunctions, depending on k       ==
    ! ==          linear combinations of the C0 and the C1            ==
    ! ==          where the latter are obtained by the perturbation   ==
    ! ==          due to the k point contribution to                  ==
    ! ==             the standard Hamiltonian, considered in Gamma.   ==
    ! ==          the perturbation contribution are imaginary in the  ==
    ! ==          real space and depend on k: i * k * C1              ==
    ! ==          even if in the perturbation algorithm the C1 are    ==
    ! ==          considered as coefficient of real functions, and    ==
    ! ==          they depend on the direction x, y, or z             ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cu1(ncpw%ngw,nstate,3), &
                                                c0(ncpw%ngw,nstate)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'up3_p'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8), ALLOCATABLE             :: adav(:,:), c00(:,:,:), &
                                                c11(:,:)
    INTEGER                                  :: i, ierr, ig, ikind, &
                                                irec(100), istat, nnx, noccup
    REAL(real_8)                             :: rsum
    REAL(real_8), ALLOCATABLE                :: ev(:,:)
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
! New Hilbert space, spanned by the 2*NSTATE wavefunction
! This set is constructed separately per each K point
! 

    ALLOCATE(c00(ncpw%ngw,nstate,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c00)!,2*ngw*2*nstate*nkpt%nkpts)
    ALLOCATE(c11(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c11)!,2*ngw*2*nstate)
    ALLOCATE(ev(2*nstate,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ev)!,2*nstate*nkpt%nkpts)
    ALLOCATE(adav(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(adav)!,4*nstate*nstate)
    ! ==--------------------------------------------------------------==
    ! == ALLOCATE ARRAYS  for the new wavefunction                    ==
    ! ==--------------------------------------------------------------==

    nnx=fpar%nnr1
    tkpts%tkpnt = .TRUE.
    tkpts%tkscale = .FALSE.
    tsphere = .FALSE.

    crge%n = 2*nstate
    IF (paral%io_parent)&
         WRITE(6,'("NGWKS 1 = ",i10)') spar%ngwks
    nkpt%ngwk = 2*ncpw%ngw
    spar%ngwks = 2*spar%ngws
    IF (paral%io_parent)&
         WRITE(6,'("NGWKS = ",i10)') spar%ngwks
    IF (paral%io_parent)&
         WRITE(6,'("NHGS = ",i10)') spar%nhgs
    IF (paral%io_parent)&
         WRITE(6,'("NHGLS = ",i10)') spar%nhgls


    ! ==--------------------------------------------------------------==
    ! == Calculate the Stress Tensor                                  ==
    ! ==--------------------------------------------------------------==

    noccup = nstate

    ! ==-------------------------------------------------------------==
    ! Preparation of the projector for the non local potential
    ! which depend on K point

    CALL rkpnt


    DO ikind = 1, nkpt%nkpts

       DO i =1,nstate
          DO ig = 1,ncpw%ngw
             c00(ig,i,ikind) = c0(ig,i)
             c00(ig+ncpw%ngw,i,ikind)=CONJG(c0(ig,i))
          ENDDO
       ENDDO

       ! For each rk(1:3,1:nkpts) a new HAMILK is created
       ! It is a block matrix, along the diagonal 

       !>>>
       CALL stopgm(procedureN,'OUT OF BOUND, please fix me',&
            __LINE__,__FILE__)
       !vw>>>> out of bound ! CALL  zazzero(c00(1,nstate+1,ikind),2*ngw*nstate)
       !<<<

       DO istat = 1, nstate
          DO ig = 1, ncpw%ngw
             c00(ig,istat+nstate,ikind) = uimag*( cu1(ig,istat,1)*rk(1,&
                  ikind) +cu1(ig,istat,2)*rk(2,ikind)  +cu1(ig,istat,3)*&
                  rk(3, ikind) )*parm%tpiba
             c00(ig+ncpw%ngw,istat+nstate,ikind) =-CONJG(c00(ig,istat+&
                  nstate,ikind))
          ENDDO
       ENDDO

       IF (geq0) CALL zclean_k(c00(1,1,ikind),2*nstate,ncpw%ngw)


       ! DO ISTAT = 1,2*NSTATE
       ! sum =0.0_real_8
       ! DO IG = 1,2*NGW
       ! sum = sum + 
       ! &              real(c00(IG,ISTAT,IKIND))*
       ! &              real(c00(IG,ISTAT,IKIND))+
       ! &              aimag(c00(IG,ISTAT,IKIND))*
       ! &              aimag(c00(IG,ISTAT,IKIND))
       ! ENDDO
       ! write(6,'("mod stat ",i5," = ",f14.8)') ISTAT, sum
       ! ENDDO




       ! write(6,'("costruiti c00 per k =  ",i4)') ikind 


       CALL gs_ortho_c(c00(1,1,ikind),nstate,c00(1,nstate+1,ikind),&
            nstate,adav)



       ! DO ISTAT = 1,2*NSTATE
       ! sum =0.0_real_8
       ! DO IG = 1,2*NGW
       ! sum = sum + 
       ! c     &              real(c00(IG,ISTAT,IKIND))*
       ! &              real(c00(IG,ISTAT,IKIND))+
       ! &              aimag(c00(IG,ISTAT,IKIND))*
       ! &              aimag(c00(IG,ISTAT,IKIND))
       ! ENDDO
       ! write(6,'("mod stat ",i5," = ",f14.8)') ISTAT, sum
       ! ENDDO


    ENDDO                     ! LOOP OVER IKIND


    ! ==--------------------------------------------------------------==
    ! == End loop over kpoints                                        ==
    ! ==--------------------------------------------------------------==
    rsum =0._real_8
    DO ikind = 1,nkpt%nkpts
       IF (paral%io_parent) THEN
          WRITE(6,'(a,i4,a,f14.6)') "WK(",ikind,")= ",wk(ikind)
       ENDIF
       DO i=1,nstate
          IF (crge%f(i,1).NE.0._real_8) THEN
             rsum=rsum+wk(ikind)*crge%f(i,1)*ddot(nkpt%ngwk*2,c00(1,i,ikind),1,&
                  c00(1,i,ikind),1)
          ENDIF
       ENDDO
    ENDDO

    CALL mp_sum(rsum,parai%allgrp)

    IF (paral%io_parent) THEN
       WRITE(6,'(A,F13.6)')"TOTAL INTEGRATED DENSITY IN FOURIER SPACE", rsuM
    ENDIF

    ! ==--------------------------------------------------------------==
    ! == Write Restart File, with new wavefunctions                   ==
    ! ==--------------------------------------------------------------==

    CALL write_irec(irec)

    irec(10) = 0

    ! write(6,'("IREC : ")')
    ! DO i =1,30
    ! write(6,'("IREC ",i4," = ",i4)') i,IREC(i)
    ! ENDDO

    ! CALL ZHWWF(2,IREC,C0,C11,NSTATE,EV,TAU0,TAU0,TAU0,1)
    CALL zhwwf(2,irec,c00,c11,2*nstate,ev,tau0,tau0,tau0,1)

    DEALLOCATE(c00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c11,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE up3_p
  ! ==================================================================
  SUBROUTINE give_scr_up3_p(lupdrho,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lupdrho
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lhpsi, lortho, lrnlsm

    CALL give_scr_rnlsm(lrnlsm,tag,2*nstate,.FALSE.)
    CALL give_scr_ortho(lortho,tag,2*nstate)
    lhpsi =2*nstate+2*maxsys%nax+                   MAX(2*2*ncpw%ngw,lrnlsm)+&
         100

    lupdrho= MAX(lhpsi,lortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_up3_p
  ! ==================================================================

END MODULE up3_p_utils
