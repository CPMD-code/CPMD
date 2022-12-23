#include "cpmd_global.h"

MODULE syscomb_utils
  USE cdftmod,                         ONLY: cdfthda,&
                                             cm_dir,&
                                             cm_dr,&
                                             sccomm
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: fbohr,&
                                             pi
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             velp
  USE cppt,                            ONLY: gk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzs,&
                                             nzfs
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE filnmod,                         ONLY: filn
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE hpsi_utils,                      ONLY: hpsi
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ortho_utils,                     ONLY: give_scr_ortho,&
                                             ortho
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE prden,                           ONLY: mwfn
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rhopri_utils,                    ONLY: rhopri
  USE rmas,                            ONLY: rmass
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: write_irec
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE store_types,                     ONLY: irec_wf,&
                                             restart1,&
                                             rout1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE utils,                           ONLY: invmat
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: syscomb
  !public :: write_sab
  PUBLIC :: write_ksham
  !public :: getcom
  !public :: mirror_0
  !public :: mirror_1
  !public :: translate_x
  !public :: qm_centre

CONTAINS

  ! ==================================================================
  ! Main source file of the wavefunction combination method
  ! also contains the KS Hamiltonian matrix output functions (FODFT)
  ! H. Oberhofer (ho246@cam.ac.uk) 2009
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE syscomb(c0,cm,nl_s0,nl_s1,nl_s0up,nl_s1up)
    ! ==--------------------------------------------------------------==
    ! == Combine two reference systems from RESTART.R* files          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), &
                                                cm(ncpw%ngw,*)
    INTEGER                                  :: nl_s0, nl_s1, nl_s0up, nl_s1up

    CHARACTER(*), PARAMETER                  :: procedureN = 'syscomb'

    CHARACTER(len=20)                        :: filn_back
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cbd(:,:), cbu(:,:), scr(:), &
                                                scr2(:)
    INTEGER                                  :: i, ierr, irec(100), lortho, &
                                                nfi, nstate
    REAL(real_8)                             :: com(3), eigv(crge%n), &
                                                taui(3,maxsys%nax,maxsys%nsx)

    filn_back=filn
    IF (paral%io_parent)&
         WRITE(6,*) 'ATTEMPTING TO COMBINE SYSTEMS'
    IF (.NOT.(restart1%rco.AND.restart1%rwf))THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RESTART OPTION NOT ACTIVATED'
       CALL stopgm('SYSCOMB','        ',& 
            __LINE__,__FILE__)
    ENDIF
    CALL zeroing(irec)!,100)
    irec(irec_wf)=1
    filn = "RESTART.R1"
    crge%n=nl_s0
    IF (cntl%tlsd)THEN
       spin_mod%nsup=nl_s0up
       spin_mod%nsdown=crge%n-spin_mod%nsup
    ENDIF
    CALL zhrwf(1,irec,c0,cm,nl_s0,eigv,tau0,velp,taui,nfi)
    filn = "RESTART.R2"
    crge%n=nl_s1
    IF (cntl%tlsd)THEN
       spin_mod%nsup=nl_s1up
       spin_mod%nsdown=crge%n-spin_mod%nsup
    ENDIF
    CALL zhrwf(1,irec,c0(:,nl_s0+1:,:),cm(1,nl_s0+1),&
         nl_s1,eigv,tau0,velp,taui,nfi)
    IF (restart1%rgeo.AND.paral%parent) CALL geofile(tau0,velp,'READ')
    CALL write_irec(irec)
    nstate = nl_s0 + nl_s1
    crge%n=nstate
    IF (cntl%tlsd)THEN
       spin_mod%nsup=nl_s0up+nl_s1up
       spin_mod%nsdown=crge%n-spin_mod%nsup
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*)"Found ",nstate,"States"

    IF (cntl%tscombm)THEN
       IF (cm_dir.EQ.0)THEN
          IF (paral%io_parent)WRITE(6,*)"MIRRORING WF 2 AROUND CENTRE OF MASS"
          ! CALL GETCOM(COM)
          DO i=nl_s0+1,nstate
             CALL mirror_0(c0(1,i,1),com)
          ENDDO
       ELSEIF (cm_dir.EQ.1)THEN
          IF (paral%io_parent)WRITE(6,*)"MIRRORING WF 2 ALONG X AXIS"
          CALL setfftn(0)
          ALLOCATE(scr(maxfft),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(scr2(maxfft),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i=sccomm%n_s0+1,nstate
             CALL mirror_1(c0(1,i,1),scr,scr2)
          ENDDO
          DEALLOCATE(scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(scr2,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSEIF (cm_dir.EQ.4)THEN
          IF (paral%io_parent)WRITE(6,*)&
               "MIRRORING AND TRANSLATING WF 2 ALONG X AXIS"
          IF (.NOT.cntl%bohr)cm_dr=cm_dr*fbohr
          CALL setfftn(0)
          ALLOCATE(scr(2*maxfft),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(scr2(2*maxfft),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i=nl_s0+1,nstate
             CALL mirror_1(c0(1,i,1),scr,scr2)
             CALL translate_x(c0(1,i,1),cm_dr)
          ENDDO
          DEALLOCATE(scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(scr2,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          IF (paral%io_parent)WRITE(6,*)&
               "TRANSLATED WF BY ",cm_dr," a.u."
       ENDIF
    ENDIF

    IF (cntl%tlsd)THEN
       IF (paral%io_parent)THEN
          WRITE(6,*)"SORTING UP AND DOWN SPINS."
          WRITE(6,*)nl_s1up,nl_s0up,nl_s1,nl_s0
       ENDIF
       ALLOCATE(cbu(ncpw%ngw,nl_s1up),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(2*ncpw%ngw*nl_s1up,c0(1,sccomm%n_s0+1,1),1,cbu,1)
       ALLOCATE(cbd(ncpw%ngw,nl_s0-nl_s0up),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(2*ncpw%ngw*(nl_s0-nl_s0up),c0(1,sccomm%n_s0up+1,1),1,cbd,1)
       CALL dcopy(2*ncpw%ngw*nl_s1up,cbu,1,c0(1,sccomm%n_s0up+1,1),1)
       CALL dcopy(2*ncpw%ngw*(nl_s0-nl_s0up),cbd,1,&
            c0(1,nl_s0up+nl_s1up+1,1),1)
    ENDIF

    IF (sccomm%tsysl.GE.0)CALL write_sab(c0,sccomm%tsysl,sccomm%tsysk)


    IF (cntl%tscortho)THEN
       IF (paral%io_parent)&
            WRITE(6,*)"ORTHOGONALISING COMBINED WF"
       CALL give_scr_ortho(lortho,tag,crge%n)
       CALL ortho(crge%n,c0(:,:,1),cm)
    ENDIF

    IF (paral%io_parent)THEN
       WRITE(6,*)"COMBINED WAVEFUNCTIONS FROM RESTART.R1 AND RESTART.R2"
       WRITE(6,*)"TO RESTART.1."
       WRITE(6,*)"EXITING NOW"
    ENDIF


    CALL zhwwf(2,irec,c0,cm,nstate,eigv,tau0,tau0,taui,nfi)
    filn = filn_back
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE syscomb
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE write_sab(c0,i,j)
    ! ==--------------------------------------------------------------==
    ! == Calculate and write out the KS Hamiltonian Matrix            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    INTEGER                                  :: i, j

    INTEGER                                  :: gl
    REAL(real_8)                             :: res

    IF (geq0)THEN
       res=REAL(c0(1,i))*REAL(c0(1,j))
    ELSE
       res=2.0_real_8*(REAL(c0(1,i))*REAL(c0(1,j))+&
            AIMAG(c0(1,i))*AIMAG(c0(1,j)))
    ENDIF
    !$omp parallel do private(GL) reduction(+:RES)
    DO gl=2,ncpw%ngw
       res=res+2.0_real_8*(REAL(c0(gl,i))*&
            REAL(c0(gl,j))+AIMAG(c0(gl,i))*AIMAG(c0(gl,j)))
    ENDDO
    CALL mp_sum(res,parai%allgrp)
    IF (paral%io_parent)THEN
       WRITE(6,*)
       WRITE(6,'(2X,"SAB  (",I3,",",I3,") :  ",F14.8)')i,j,reS
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_sab
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE write_ksham(c0,cm,sc0,rhoe,psi,eigv)
    ! ==--------------------------------------------------------------==
    ! == Calculate and write out the KS Hamiltonian Matrix            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: c0(ncpw%ngw,crge%n)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,crge%n), sc0(*)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'write_ksham'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), POINTER                 :: C0_3d_p(:,:,:)
    INTEGER                                  :: i, ierr, info, j, k, l, &
                                                lforces
    REAL(real_8)                             :: fsave, haba, habat, habb, &
                                                habbt, hbaa, hbaat, hbab, &
                                                hbabt, hmat(crge%n,crge%n), &
                                                sab, sba
    REAL(real_8), ALLOCATABLE                :: a(:,:), ai(:,:), hmatt(:,:), &
                                                scr(:)

!(nnr1,clsd%nlsd)
! Variables
! HABAT&HABT are the KS Matrices in the non-orthogonal bases
! (N,*)
! (N,*)

    CALL give_scr_forcedr(lforces,tag,crge%n,.FALSE.,.TRUE.)
    ALLOCATE(scr(lforces),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    l=sccomm%n_s0up
    k=sccomm%n_s1up
    IF (cntl%tlowdmat)THEN
       ALLOCATE(a(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ai(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hmatt(crge%n,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)THEN
          OPEN(41,file="LOWDIN_A",status='OLD')
          REWIND(unit=41)
          DO i=1,crge%n
             READ(41,*)(a(i,j),j=1,crge%n)
          ENDDO
          CLOSE(41)
       ENDIF
       CALL mp_bcast(a,SIZE(a),parai%io_source,parai%cp_grp)
       CALL dcopy(crge%n*crge%n,a,1,ai,1)
       CALL invmat(crge%n,ai,hmatt,info)
       CALL dgemm('n','n',crge%n,crge%n,crge%n,1.0_real_8,ai,crge%n,ai,crge%n,0.0_real_8,hmatt,crge%n)
       sab=hmatt(k,l)
       sba=hmatt(l,k)
       IF(paral%io_parent)THEN
          OPEN(41,file="OVERLAP",status='UNKNOWN')
          REWIND(unit=41)
          DO i=1,crge%n
             WRITE(41,*)(hmatt(i,j),j=1,crge%n)
          ENDDO
          CLOSE(41)
       ENDIF
    ENDIF
    IF(paral%io_parent)WRITE(6,*)"CALCULATING KS-HAMILTONIAN, NO RESET OF OCCUPATIONS"
    CALL rhoofr(c0(:,1:crge%n),rhoe,psi(:,1),crge%n)
    CALL forcedr(c0,cm,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.FALSE.,.TRUE.)
    CALL hpsi (c0, cm, sc0, rhoe, psi(:,1), crge%n, 1, clsd%nlsd)
    CALL ovlap(crge%n,hmat,c0,cm)
    CALL mp_sum(hmat,crge%n*crge%n,parai%allgrp)
    IF(cntl%tlsd)THEN
       haba=-hmat(k,l)
       hbaa=-hmat(l,k)
    ELSE
       haba=-hmat(k,l)/2._real_8 !division by two
       hbaa=-hmat(l,k)/2._real_8
    ENDIF
    IF(paral%io_parent)THEN
       OPEN(41,file="KS_HAM",status='UNKNOWN')
       REWIND(unit=41)
       DO i=1,crge%n
          WRITE(41,*)(-hmat(i,j),j=1,crge%n)
       ENDDO
       CLOSE(41)
       WRITE(6,'(1x,"HAB,HBA IN ORTHOG. BASIS    : ",2(F13.8,1X),"a.u.")') haba,hbaa
    ENDIF
    IF(cntl%tlowdmat)THEN
       CALL dgemm('n','n',crge%n,crge%n,crge%n,1.0_real_8,hmat,crge%n,a,crge%n,0.0_real_8,hmatt,crge%n)
       CALL dcopy(crge%n*crge%n,hmatt,1,hmat,1)
       CALL dgemm('n','n',crge%n,crge%n,crge%n,1.0_real_8,ai,crge%n,hmat,crge%n,0.0_real_8,hmatt,crge%n)
       IF(cntl%tlsd)THEN
          habat=-hmatt(k,l)
          hbaat=-hmatt(l,k)
       ELSE
          habat=-hmatt(k,l)/2._real_8 !division by two
          hbaat=-hmatt(l,k)/2._real_8
       ENDIF
       IF(paral%io_parent) THEN
          WRITE(6,'(1x,"HAB,HBA IN NON-ORTHOG. BASIS: ",2(F13.8,1X),"a.u.",/)')habat,hbaat
       ENDIF
    ENDIF

    IF (cntl%tlsd) THEN
       IF(paral%io_parent)WRITE(6,*)"CALCULATING KS-HAMILTONIAN A"
       IF(k.GT.spin_mod%nsup)THEN
          IF(paral%io_parent)&
               WRITE(6,'(1x,"RESETTING OCCUPATION NUMBER OF BETA STATE",I3," TO ZERO")') k-spin_mod%nsup
       ELSE
          IF(paral%io_parent)&
               WRITE(6,'(1x,"RESETTING OCCUPATION NUMBER OF ALPHA STATE ",I3," TO ZERO")') k
       ENDIF
       fsave=crge%f(k,1)
       crge%f(k,1)=0.0_real_8
       CALL rhoofr(c0(:,1:crge%n),rhoe,psi(:,1),crge%n)
       CALL forcedr(c0,cm,sc0,rhoe,psi,tau0,fion,eigv,&
            crge%n,1,.FALSE.,.TRUE.)
       CALL hpsi (c0, cm, sc0, rhoe, psi(:,1), crge%n, 1, clsd%nlsd)
       CALL ovlap(crge%n,hmat,c0,cm)
       CALL mp_sum(hmat,crge%n*crge%n,parai%allgrp)
       haba=-hmat(k,l)
       hbaa=-hmat(l,k)

       IF (paral%io_parent)THEN
          OPEN(41,file="KS_HAM_A",status='UNKNOWN')
          REWIND(unit=41)
          DO i=1,crge%n
             WRITE(41,*)(-hmat(i,j),j=1,crge%n)
          ENDDO
          CLOSE(41)
          WRITE(6,'(1X,"HAB,HBA IN ORTHOG. BASIS    : ",2(F13.8,1x),"a.u.")')haba,hbaa
       ENDIF

       IF (cntl%tlowdmat)THEN
          CALL dgemm('N','N',crge%n,crge%n,crge%n,1.0_real_8,hmat,crge%n,a,crge%n,0.0_real_8,hmatt,crge%n)
          CALL dcopy(crge%n*crge%n,hmatt,1,hmat,1)
          CALL dgemm('N','N',crge%n,crge%n,crge%n,1.0_real_8,ai,crge%n,hmat,crge%n,0.0_real_8,hmatt,crge%n)
          habat=-hmatt(k,l)
          hbaat=-hmatt(l,k)
          IF (paral%io_parent)WRITE&
               (6,'(1X,"HAB,HBA IN NON-ORTHOG. BASIS: ",2(F13.8,1x),"a.u.",/)')habat,hbaat
       ENDIF

       IF (paral%io_parent)THEN
          IF(paral%io_parent)WRITE(6,*)"CALCULATING KS-HAMILTONIAN B"
          IF (l.GT.spin_mod%nsup)THEN
             WRITE(6,'(1X,"RESETTING OCCUPATION NUMBER OF BETA STATE ", I3," TO ZERO")') l-spin_mod%nsup
          ELSE
             WRITE(6,'(1X,"RESETTING OCCUPATION NUMBER OF ALPHA STATE ", I3, " TO ZERO")') L
          ENDIF
       ENDIF
       crge%f(k,1)=fsave
       crge%f(l,1)=0.0_real_8
       CALL rhoofr(c0(:,1:crge%n),rhoe,psi(:,1),crge%n)
       CALL forcedr(c0,cm,sc0,rhoe,psi,tau0,fion,eigv,&
            crge%n,1,.FALSE.,.TRUE.)
       CALL hpsi (c0, cm, sc0, rhoe, psi(:,1), crge%n, 1, clsd%nlsd)
       CALL ovlap(crge%n,hmat,c0,cm)
       CALL mp_sum(hmat,crge%n*crge%n,parai%allgrp)

       habb=-hmat(l,k)
       hbab=-hmat(k,l)

       IF (paral%io_parent)THEN
          OPEN(41,file="KS_HAM_B",status='UNKNOWN')
          REWIND(unit=41)
          DO i=1,crge%n
             WRITE(41,*)(-hmat(i,j),j=1,crge%n)
          ENDDO
          CLOSE(41)
          WRITE(6,'(1X,"HAB,HBA IN ORTHOG. BASIS    : ",2(F13.8,1x),"a.u.")')habb,hbab
       ENDIF
       IF (cntl%tlowdmat)THEN
          CALL dgemm('N','N',crge%n,crge%n,crge%n,1.0_real_8,hmat,crge%n,a,crge%n,0.0_real_8,hmatt,crge%n)
          CALL dcopy(crge%n*crge%n,hmatt,1,hmat,1)
          CALL dgemm('N','N',crge%n,crge%n,crge%n,1.0_real_8,ai,crge%n,hmat,crge%n,0.0_real_8,hmatt,crge%n)
          habbt=-hmatt(k,l)
          hbabt=-hmatt(l,k)
          IF (paral%io_parent)&
               WRITE(6,'(1X,"HAB,HBA IN NON-ORTHOG. BASIS: ",2(F13.8,1x),"a.u.",/)')habbt,hbabt
       ENDIF

       IF (paral%io_parent)THEN
          WRITE(6,'(1X,"AVERAGE HAB ORTHOG BASIS    : ",F13.8,1x,"a.u.")')(haba+habb)/2.0_real_8
          IF (cntl%tlowdmat)&
               WRITE(6,'(1X,"AVERAGE HAB NON-ORTHOG BASIS: ",F13.8,1x,"a.u.")')&
               (habat+hbaat+habbt+hbabt)/4.0_real_8
       ENDIF
    ENDIF !LSD
    IF(cntl%tlowdmat.AND.paral%io_parent) WRITE(6,'(1X,"OVERLAP MATRIX EL. SAB,SBA  : ",2(F13.8,1x),/)')&
         sab,sba
    IF(rout1%rhoout) THEN
       cdfthda%hdafirst=.TRUE.
       mwfn(1)=-k
       mwfn(2)=-l
       IF(k.GT.spin_mod%nsup)mwfn(1)=mwfn(1)+spin_mod%nsup
       IF(l.GT.spin_mod%nsup)mwfn(2)=mwfn(2)+spin_mod%nsup
#if defined( _HASNT_F08_POINTER_REMAPPING )
       CALL stopgm(procedureN,'compiler needs to support pointer remapping!',&
            __LINE__,__FILE__)
#else
       c0_3d_p(1:SIZE(c0,1),1:SIZE(c0,2),1:1) => c0
#endif
       CALL rhopri(c0_3d_p,tau0,rhoe,psi(:,1),crge%n,1)
    ENDIF
    IF(paral%io_parent) WRITE(6,*)"EXITING NOW"
    ! ==--------------------------------------------------------------==
    ! Deallocation of local variables
    IF (cntl%tlowdmat) DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_ksham
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE getcom(com)
    ! ==--------------------------------------------------------------==
    ! == Calculate the position of the Centre of Mass                 ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: com(3)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: co(3), pma00

    pma00=0.0_real_8
    co(:) = 0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          co(1)=co(1)+tau0(1,ia,is)*rmass%pma0(is)
          co(2)=co(2)+tau0(2,ia,is)*rmass%pma0(is)
          co(3)=co(3)+tau0(3,ia,is)*rmass%pma0(is)
          pma00=pma00+rmass%pma0(is)
       ENDDO
    ENDDO
    com(1)=co(1)/pma00
    com(2)=co(2)/pma00
    com(3)=co(3)/pma00
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getcom
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE mirror_0(wave,com)
    ! ==--------------------------------------------------------------==
    ! == Mirror wavefunctions (Inversion around centre)               ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: wave(ncpw%ngw)
    REAL(real_8)                             :: com(3)

    INTEGER                                  :: ig

    DO ig=1,ncpw%ngw
       wave(ig)=CONJG(wave(ig))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mirror_0
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE mirror_1(wave,scr,scr2)
    ! ==--------------------------------------------------------------==
    ! == Mirror wavefunctions (around ZY-plane)                       ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: wave(ncpw%ngw), &
                                                scr(fpar%nnr1), scr2(:)

    INTEGER                                  :: g, gx, gy, gy2, gz, gz2

    CALL zeroing(scr)!,nnr1)
    DO g=1,ncpw%ngw
       scr(nzfs(g))=CONJG(wave(g))
       scr(inzs(g))=wave(g)
    ENDDO
    CALL invfftn(scr,.FALSE.,parai%allgrp)
    DO gx=1,fpar%kr1
       DO gy=1,fpar%kr2s
          gy2=fpar%kr2s-gy
          DO gz=1,fpar%kr3s
             gz2=fpar%kr3s-gz
             scr2(gx+gy2*fpar%kr1+gz2*fpar%kr1*fpar%kr2s)&
                  =scr(gx+(gy-1)*fpar%kr1+(gz-1)*fpar%kr1*fpar%kr2s)
          ENDDO
       ENDDO
    ENDDO
    CALL fwfftn(scr2,.FALSE.,parai%allgrp)
    DO g=1,ncpw%ngw
       wave(g)=scr2(nzfs(g))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mirror_1
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE translate_x(wave,r)
    ! ==--------------------------------------------------------------==
    ! == Translate wavefunctions in G space along X axis              ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: wave(ncpw%ngw)
    REAL(real_8)                             :: r

    INTEGER                                  :: ii
    REAL(real_8)                             :: f, m

! F=B1(1)*2.0_real_8*PI/ALAT

    f=2.0_real_8*pi/parm%alat
    !$omp parallel do private(II)
    DO ii=1,ncpw%ngw
       m=-gk(1,ii)*r*f
       wave(ii)=wave(ii)*CMPLX(COS(m),SIN(m),kind=real_8)
    ENDDO

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE translate_x
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE qm_centre()
    ! ==--------------------------------------------------------------==
    ! == reposition Centre of Mass                                    ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: cx, cy, cz, dcx, dcy, dcz

    cx=0.0_real_8
    cy=0.0_real_8
    cz=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          cx=cx+tau0(1,ia,is)*rmass%pma0(is)/rmass%pmat0
          cy=cy+tau0(2,ia,is)*rmass%pma0(is)/rmass%pmat0
          cz=cz+tau0(3,ia,is)*rmass%pma0(is)/rmass%pmat0
       ENDDO
    ENDDO
    dcx=0.5_real_8*cell_com%celldm(1)-cx
    dcy=0.5_real_8*cell_com%celldm(2)*cell_com%celldm(1)-cy
    dcz=0.5_real_8*cell_com%celldm(3)*cell_com%celldm(1)-cz
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          tau0(1,ia,is)=tau0(1,ia,is)+dcx
          tau0(2,ia,is)=tau0(2,ia,is)+dcy
          tau0(3,ia,is)=tau0(3,ia,is)+dcz
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(/,A,A,/)')&
         ' >>>>>>>> CENTER OF MASS HAS BEEN MOVED',&
         ' TO CENTER OF BOX <<<<<<<<'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE qm_centre
  ! ==================================================================

END MODULE syscomb_utils
