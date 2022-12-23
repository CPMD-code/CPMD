MODULE phonons_p_utils
  USE adat,                            ONLY: elem
  USE coor,                            ONLY: tau0
  USE cotr,                            ONLY: cotc0,&
                                             hess
  USE d_mat_p_utils,                   ONLY: d_mat_diag_nonloc,&
                                             d_mat_diag_real,&
                                             d_mat_loc0,&
                                             d_mat_nonloc,&
                                             d_mat_real
  USE eicalc_utils,                    ONLY: eicalc1
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnonloc_p_utils,                 ONLY: fnonloc_p
  USE forces_p_utils,                  ONLY: give_scr_forces_p
  USE hessout_utils,                   ONLY: hessout
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             ndfnl
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00,&
                                             response1,&
                                             rho0
  USE rhoofr_p_utils,                  ONLY: give_scr_rhoofr_p
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_p_utils,                   ONLY: rnlsm3
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE secder_utils,                    ONLY: molvib,&
                                             purged,&
                                             vibeig
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE spin,                            ONLY: clsd
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: symmat
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: symma,&
                                             unitmx
!!use d_mat_p_utils, only : d_mat_locps
!!use d_mat_p_utils, only : d_mat_diag_loc
!!use nmr_util_p_utils, only : ffttog
!!use nmr_util_p_utils, only : fft2tog
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phonons_p
  !public :: rt_out
  PUBLIC :: setrot
  PUBLIC :: settras
  PUBLIC :: give_scr_phonon

CONTAINS

  ! ==================================================================
  SUBROUTINE phonons_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), &
                                                drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'phonons_p'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/'x','y','z'/)

    COMPLEX(real_8), ALLOCATABLE             :: eirop1(:), rho0G(:), &
                                                v1_loc(:,:), v1_nonloc(:,:)
    INTEGER :: dimtot, i, ia, ia2, iat, iat2, icol, ierr, il_eirop1, &
      il_v1_loc, il_v1_nonloc, info, is, is2, isub, ix, j, k, k2, ldfnl, &
      lfnl, lwork, nsym
    REAL(real_8)                             :: xmass
    REAL(real_8), ALLOCATABLE :: ddfnl00(:,:,:,:,:,:,:), eigen(:,:), &
      hesssave(:), rot(:,:), sder(:,:), sdersave(:), tras(:,:), vibe(:), &
      work(:), xma(:)

! ==--------------------------------------------------------------==

    CALL tiset('    phonon',isub)
    ! ==--------------------------------------------------------------==
    ! Potentials and related:
    il_v1_loc=2*ncpw%nhg*clsd%nlsd
    il_v1_nonloc=2*nkpt%ngwk*nstate*nkpt%nkpnt
    il_eirop1=2*ncpw%nhg
    ALLOCATE(v1_loc(ncpw%nhg,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_loc)!,SIZE(v1_loc))
    ALLOCATE(v1_nonloc(nkpt%ngwk,il_v1_nonloc/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_nonloc)!,SIZE(v1_nonloc))
    ! NB: In v1_nonloc, we store \op{v1} |psi0> whereas in
    ! v1_loc, v1(g) is stored. Therefore the different dimensions.
    ALLOCATE(eirop1(il_eirop1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop1)!,SIZE(eirop1))


    ! ==--------------------------------------------------------------==
    ! the gnd-state density in G-space:
    ALLOCATE(rho0G(ncpw%nhg*clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rho0G)!,nhg*clsd%nlsd)
    IF (cntl%tlsd) THEN
       CALL fft2tog(rho0(:,1),rho0G(1:ncpw%nhg),rho0(:,2),rho0g(ncpw%nhg+1:),&
            psi,ncpw%nhg,.TRUE.)
    ELSE
       CALL ffttog(rho0,rho0G,psi,ncpw%nhg,.TRUE.)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! FNL,DFNL,DDFNL...
    ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    lfnl = imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt

    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dfnl00(imagp,ions1%nat,maxsys%nhxs,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddfnl00(imagp,ions1%nat,maxsys%nhxs,3,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rnlsm(c0,nstate,1,1,.TRUE.)
    CALL rnlsm3(c0,nstate,ddfnl00)
    CALL dcopy(lfnl,   fnl,1,  fnl00,1)
    CALL dcopy(ldfnl, dfnl,1, dfnl00,1)

    ! ==--------------------------------------------------------------==
    ! Specific ionic matrices (hessian, rotation, transl matrix)
    ALLOCATE(tras(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rot(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(tras)!,9*ions1%nat)
    CALL zeroing(rot)!,9*ions1%nat)

    cotc0%nodim=3*ions1%nat
    ALLOCATE(sder(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sder)!,9*ions1%nat*ions1%nat)
    ALLOCATE(hesssave(9*ions1%nat*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(hesssave)!,9*ions1%nat*ions1%nat)
    ALLOCATE(sdersave(9*ions1%nat*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sdersave)!,9*ions1%nat*ions1%nat)
    IF (paral%parent) THEN
       ALLOCATE(eigen(3*ions1%nat,3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xma(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vibe(3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(hess(cotc0%nodim,cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cotc0%nodim.NE.3*ions1%nat) THEN
       CALL stopgm('PHONONS_P','Dimension must be 3*Nat',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==         the basic loops for perturbations                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/," ",22("*"),a,22("*"),/)')&
         '   perturbations    '
    ! no symmetrization of gradients in this part!
    nsym=symmi%nrot
    symmi%nrot=1
    iat=0
    ! ..   calculate perturbation wrt. to atom IAT and coordinate K
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             IF (paral%io_parent)&
                  WRITE(6,'(" **** atom=",i7,4x,a4,4x,a4,4x)')&
                  iat,elem%el(ions0%iatyp(is)),cdir(k)

             ! ...  local part of the perturbation:
             ! v1_loc(ig)  =  d/dR PP_local(r-R)  ==  PP_local(G) * G
             CALL eicalc1(k,is,iat,v1_loc,eirop1)
             IF (cntl%tlsd)&
                  CALL dcopy(2*ncpw%nhg,v1_loc(1,1),1,v1_loc(1,2),1)
             ! This potential is identical for up and down spins

             ! ...  non-local part:
             ! |v1_nonloc> = (d/dR [PP-Projectors]) |0>
             CALL fnonloc_p(c0,psi,&
                  v1_nonloc,crge%f,nstate,k,is,iat,1)

             ! ==--------------------------------------------------------------==
             CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
                  eirop,eivps,v1_loc,v1_nonloc,&
                  z11,nstate,eirop1)
             ! ==--------------------------------------------------------------==

             ! rwfopt_p returns c1 and the first-order perturbation density in
             ! real space representation in the array DRHOE.
             CALL ffttog(drhoe,drhoe,psi,ncpw%nhg,.TRUE.)

             ! second derivatives matrix:
             icol=3*(iat-1)+k
             iat2=0
             DO is2=1,ions1%nsp
                DO ia2=1,ions0%na(is2)
                   iat2=iat2+1
                   DO k2=1,3
                      ix=3*(iat2-1)+k2
                      IF (iat.EQ.iat2) THEN
                         ! local part, i=j (i,j=atoms)
                         CALL d_mat_diag_loc(sder(ix,icol),rho0G,&
                              eirop,eirop1,v1_loc,k2)
                         ! non-local part, i=j (i,j=atoms)
                         CALL d_mat_diag_nonloc(sder(ix,icol),&
                              ddfnl00,fnl00,&
                              dfnl00,crge%f,wk,is,k,k2,iat,nstate,1)
                         ! real space term, i=j
                         CALL d_mat_diag_real(sder(ix,icol),&
                              is,ia,k,k2,tau0)
                         ! real space term, i <> j
                      ELSE
                         CALL d_mat_real(sder(ix,icol),is,ia,k,k2,&
                              is2,ia2,tau0)
                      ENDIF
                      ! local part      
                      CALL d_mat_loc0(sder(ix,icol),&
                           eirop1,iat2,k2,is2)
                      CALL d_mat_locps(sder(ix,icol),&
                           drhoe,iat2,k2,is2)
                      ! non local part
                      CALL rnlsm(c1,nstate,1,1,.TRUE.)
                      CALL d_mat_nonloc(sder(ix,icol),&
                           fnl,dfnl,fnl00,dfnl00,&
                           crge%f,wk,is2,k,k2,iat2,nstate,1)
                   ENDDO! k2
                ENDDO! ia2
             ENDDO   ! is2
             ! end of one C1-optimization-cycle.
          ENDDO         ! k
       ENDDO               ! ia
    ENDDO                     ! is


    dimtot=9*ions1%nat*ions1%nat
    CALL mp_sum(sder,dimtot,parai%allgrp)

    symmi%nrot=nsym
    IF (paral%io_parent) THEN
       IF (symmi%indpg.NE.0 .AND. (symmi%nrot.NE.1 .OR. symmi%ntvec .GT. 1)) THEN
          WRITE (6,*) 'Symmetrizing second derivatives by SYMMAT'
          ! symmetrize hessian
          CALL symmat(sder,2)
       ENDIF
       CALL symma(sder,3*ions1%nat)
       ! write hessian to file hessian
       CALL dcopy(dimtot,sder,1,hess,1)
       CALL dcopy(dimtot,hess,1,hesssave,1)
       CALL hessout
       ! write output file for molvib program
       CALL molvib(tau0,sder)
       ! mass weighted force constants
       ix=0
       DO is=1,ions1%nsp
          xmass=SQRT(1._real_8/rmass%pma0(is))
          DO ia=1,ions0%na(is)
             DO k=1,3
                ix=ix+1
                xma(ix)=xmass
             ENDDO
          ENDDO
       ENDDO
       !$omp parallel do private(i,j)
       DO i=1,3*ions1%nat
          DO j=1,3*ions1%nat
             sder(i,j)=sder(i,j)*xma(i)*xma(j)
          ENDDO
       ENDDO
       CALL dcopy(dimtot,sder,1,sdersave,1)

       ! ---  Purification, version 1 (no projection)  ------------
       WRITE (6,*)
       WRITE (6,'(66("*"))')
       ! diagonalization
       lwork=dimtot
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,work,lwork,info)
       IF (info.NE.0) WRITE(6,*) 'PHONONS_P: DSYEV failed (v1).'
       !$omp parallel do private(i)
       DO i=1,3*ions1%nat
          vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
       ENDDO
       WRITE (6,*) 'Raw second derivative matrix: No'&
            //'trans/rot elimination'
       WRITE (6,*) ' Harmonic frequencies in cm^-1:'
       WRITE (6,'(4(f12.1))') (vibe(i),i=1,3*ions1%nat)
       WRITE (6,'(A,e11.5)') ' ChkSum(PHONON) = ',SUM(ABS(vibe))
       CALL vibeig(vibe,eigen,3*ions1%nat,.TRUE.)
       ! ---  Purification, version 2 (original phonons_p) --------
       IF (response1%projout) THEN
          CALL dcopy(dimtot,sdersave,1,sder,1)
          ! purification from translation
          CALL settras(tras)
          WRITE (6,*) 'Translations eliminated from'//&
               ' second derivative matrix'
          IF (response1%rotout) THEN
             WRITE (6,*) 'Rotations eliminated from'//&
                  ' second derivative matrix'
             CALL setrot(rot,tras)
          ENDIF
          CALL rt_out(tras,rot,sder,ions1%nat,response1%rotout)
          ! diagonalization
          CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,work,lwork,info)
          IF (info.NE.0) WRITE(6,*) 'PHONONS_P: DSYEV failed (v2).'
          !$omp parallel do private(i)
          DO i=1,3*ions1%nat
             vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
          ENDDO
          WRITE (6,*) ' Harmonic frequencies in cm^-1:'
          WRITE (6,'(4(f12.1))') (vibe(i),i=1,3*ions1%nat)
          CALL vibeig(vibe,eigen,3*ions1%nat,.FALSE.)
       ENDIF
       ! ---  Purification, version 3 (lr_* variant) --------------
       CALL DCOPY(dimtot,hesssave,1,sder,1) ! not yet mass-weighted!
       CALL dcopy(dimtot,hesssave,1,hess,1)
       ! tu   --- purification only AFTER mass weighting of the Hessian
       ! tu      call purged(sder,hess,tau0)
       !$omp parallel do private(i,j)
       DO i=1,3*ions1%nat
          DO j=1,3*ions1%nat
             sder(i,j)=sder(i,j)*xma(i)*xma(j)
          ENDDO
       ENDDO
       CALL purged(sder,hess,tau0)
       CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,work,lwork,info)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (info.NE.0) WRITE(6,*) 'PHONONS_P: DSYEV failed (v3).'
       !$omp parallel do private(i)
       DO i=1,3*ions1%nat
          vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
       ENDDO
       WRITE (6,*) 'Projection (by PURGED routine):'
       WRITE (6,*) ' Harmonic frequencies in cm^-1:'
       WRITE (6,'(4(f12.1))') (vibe(i),i=1,3*ions1%nat)
       CALL vibeig(vibe,eigen,3*ions1%nat,.FALSE.)
    ENDIF
    ! 
    DEALLOCATE(rho0G,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_nonloc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       DEALLOCATE(eigen,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vibe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xma,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hess,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(sder,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sdersave,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tras,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('    phonon',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE phonons_p
  ! ==================================================================
  SUBROUTINE rt_out(tras,rot,sder,nat,rotout)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nat
    REAL(real_8)                             :: sder(3*nat,3*nat), &
                                                rot(3*nat,*), tras(3*nat,*)
    LOGICAL                                  :: rotout

    CHARACTER(*), PARAMETER                  :: procedureN = 'rt_out'

    INTEGER                                  :: ierr
    REAL(real_8), ALLOCATABLE                :: proj(:,:), sder1(:,:)

    ALLOCATE(sder1(3*nat,3*nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sder1)!,9*nat*nat)
    ALLOCATE(proj(3*nat,3*nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(proj)!,9*nat*nat)
    ! c ... P <-- 1-P      
    ! do i=1,3*nat
    ! proj(i,i)=proj(i,i)+1._real_8
    ! enddo


    CALL unitmx(proj,3*nat)
    ! ... P <-- P - T.T^t
    CALL dgemm ('N','T',3*nat,3*nat,3,-1._real_8,tras,3*nat,&
         tras,3*nat,1._real_8,proj,3*nat)

    ! ... P <-- P - R.R^t
    IF (rotout) THEN
       CALL dgemm ('N','T',3*nat,3*nat,3,-1._real_8,rot,3*nat,&
            rot,3*nat,1._real_8,proj,3*nat)
    ENDIF


    ! ... (1-P).H
    CALL dgemm('N','N',3*nat,3*nat,3*nat,1._real_8,proj,3*nat,&
         sder,3*nat,0._real_8,sder1,3*nat)

    ! ... H.(1-P)
    CALL zeroing(sder)!,9*nat*nat)
    CALL dgemm('N','N',3*nat,3*nat,3*nat,1._real_8,sder1,3*nat,&
         proj,3*nat,0._real_8,sder,3*nat)

    DEALLOCATE(sder1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(proj,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE rt_out
  ! ==================================================================
  SUBROUTINE setrot(rot,tras)
    ! ==--------------------------------------------------------------==


    REAL(real_8)                             :: rot(3*ions1%nat,3), &
                                                tras(3*ions1%nat,3)

    EXTERNAL                                 :: ddot
    INTEGER                                  :: i, ia, iat, indx(3), is, ix, &
                                                k, nrot, rr
    REAL(real_8)                             :: ccm(3), cx, cy, cz, ddot, &
                                                norm_l, qt(3)

    rmass%pmat0=0._real_8
    DO is=1,ions1%nsp
       rmass%pmat0=rmass%pmat0+rmass%pma0(is)*ions0%na(is)
    ENDDO

    ! ... definition of the COM
    cx=0._real_8
    cy=0._real_8
    cz=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          cx=cx+tau0(1,ia,is)*rmass%pma0(is)/rmass%pmat0
          cy=cy+tau0(2,ia,is)*rmass%pma0(is)/rmass%pmat0
          cz=cz+tau0(3,ia,is)*rmass%pma0(is)/rmass%pmat0
       ENDDO
    ENDDO

    ! ... actual calculation of the rotational vectors
    nrot=3
    DO rr=1,nrot

       indx(1)=rr
       indx(2)=rr+1
       IF (indx(2).GT.3) indx(2)=indx(2)-3
       indx(3)=rr+2
       IF (indx(3).GT.3) indx(3)=indx(3)-3

       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ccm(1)=tau0(1,ia,is)-cx
             ccm(2)=tau0(2,ia,is)-cy
             ccm(3)=tau0(3,ia,is)-cz

             iat=iat+1
             qt(indx(1))=0._real_8
             qt(indx(2))=-ccm(indx(3))
             qt(indx(3))=ccm(indx(2))
             norm_l = SQRT(ddot(3,qt,1,qt,1))

             IF (rr.EQ.2) THEN
                DO k=1,3
                   qt(k)=-qt(k)
                ENDDO
             ENDIF
             DO k=1,3
                ix=3*(iat-1)+k
                rot(ix,rr)=qt(k)
             ENDDO

          ENDDO
       ENDDO

       ! ...  ortohgonalization wrt translation vectors
       DO i=1,3
          norm_l = - ddot(3*ions1%nat,rot(1,rr),1,tras(1,i),1)
          CALL daxpy(3*ions1%nat, norm_l, tras(1,i),1,rot(1,rr),1)
       ENDDO

       ! .... orthogonalization wrt the previous rotation vectors
       IF (rr.GE.2) THEN
          DO i=1,rr-1
             norm_l = - ddot(3*ions1%nat,rot(1,rr),1,rot(1,i),1)
             CALL daxpy(3*ions1%nat, norm_l, rot(1,i),1,rot(1,rr),1)
          ENDDO
       ENDIF

       ! ... normalization to 1
       norm_l = SQRT(ddot(3*ions1%nat,rot(1,rr),1,rot(1,rr),1))
       CALL dscal(3*ions1%nat, 1._real_8/norm_l,rot(1,rr),1)
    ENDDO

    RETURN
  END SUBROUTINE setrot
  ! ==================================================================
  SUBROUTINE settras(tras)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tras(3*ions1%nat,3)

    EXTERNAL                                 :: ddot
    INTEGER                                  :: is, j, l, m
    REAL(real_8)                             :: ddot, norm_l

! ... definition of the vectors describing the COM translation, 
! ... vectors orthogonal by definition, tras(3*nat,3)
! ... stopping COM

    IF (cntl%tinr) THEN
       DO m=1,3
          l=0
          DO is=1,ions1%nsp
             DO j=1,3*ions0%na(is)
                l=l+1
                IF (MOD((l+m-1),3).EQ.0) THEN
                   tras(l,m)=1._real_8*rmass%pma0(is)
                ELSE
                   tras(l,m)=0._real_8
                ENDIF
             ENDDO
          ENDDO
          norm_l = SQRT(ddot(3*ions1%nat,tras(1,m),1,tras(1,m),1))
          CALL dscal(3*ions1%nat, 1._real_8/norm_l,tras(1,m),1)
       ENDDO
    ELSE
       ! ... mass weighted coordinates
       DO m=1,3
          l=0
          DO is=1,ions1%nsp
             DO j=1,3*ions0%na(is)
                l=l+1
                IF (MOD((l+m-1),3).EQ.0) THEN
                   tras(l,m)=1._real_8*SQRT(rmass%pma0(is))
                ELSE
                   tras(l,m)=0._real_8
                ENDIF
             ENDDO
          ENDDO
          norm_l = SQRT(ddot(3*ions1%nat,tras(1,m),1,tras(1,m),1))
          CALL dscal(3*ions1%nat, 1._real_8/norm_l,tras(1,m),1)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE settras
  ! ==================================================================
  SUBROUTINE give_scr_phonon(lphonon,tag,nstate)
    INTEGER                                  :: lphonon
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lforce1, lrho, lrhoofr, &
                                                lrnlsm, lrscpot

! ==--------------------------------------------------------------==

    lphonon=nstate
    ropt_mod%calste=.FALSE.
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    CALL give_scr_rhoofr(lrho,tag)
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    CALL give_scr_rhoofr_p(lrhoofr,tag)
    CALL give_scr_forces_p(lforce1,tag,nstate)
    lrho=lrho+2*ncpw%nhg*clsd%nlsd
    lrhoofr=lrhoofr+2*ncpw%nhg*clsd%nlsd
    lforce1=lforce1
    lphonon=MAX(lrho,lrscpot,lrnlsm,lrhoofr,lphonon,lforce1)
    lphonon=lphonon+nstate*ncpw%ngw*2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_phonon
  ! ==================================================================

END MODULE phonons_p_utils
