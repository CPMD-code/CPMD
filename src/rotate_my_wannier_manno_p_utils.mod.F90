MODULE rotate_my_wannier_manno_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cppt,                            ONLY: gk,&
                                             hg
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE fnlalloc_utils,                  ONLY: fnlalloc
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_driver,                   ONLY: forces
  USE forces_utils,                    ONLY: give_scr_forces
  USE gndstate_p_utils,                ONLY: gauss_gen_read
  USE kinds,                           ONLY: real_8
  USE legendre_p_utils,                ONLY: rsphgen
  USE metr,                            ONLY: metr_com
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: paral
  USE phfac_utils,                     ONLY: phfac
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: &
       bdxvecx, bdxvecy, bdxvecz, btheta, bvecx, bvecy, bvecz, dmbi, &
       statemax, statemin
  USE symm,                            ONLY: sxscale,&
                                             syscale,&
                                             szscale
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
!!use gndstate_p_utils, only : save_my_wannier_gaussian
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotatemywannier_manno
  PUBLIC :: save_my_wannier
  PUBLIC :: read_my_wannier
  PUBLIC :: read_multi_wannier
  !public :: shift_my_wannier_manno
  PUBLIC :: ortho_my_wannier
  PUBLIC :: build_bylm_sphgen
  PUBLIC :: build_rmat
  PUBLIC :: rotatecoeff
  PUBLIC :: buildqmat
  PUBLIC :: buildqmatsgen
  !public :: calc_wannier_energy_manno
  PUBLIC :: IDelta

CONTAINS

  ! ==================================================================
  ! Subroutine providing the entry + necessary initialisation calls for
  ! the rotation of wannier functions expressed as a LC of momentum   
  ! vectors (g vectors)                                              
  ! ==================================================================
  SUBROUTINE rotatemywannier_manno(nstate,c0,psi)
    ! ==--------------------------------------------------------------==
    ! 
    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), psi(:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'rotatemywannier_manno'

    COMPLEX(real_8), ALLOCATABLE             :: BCrot(:,:)
    INTEGER                                  :: combi, i, ierr, imol, indx, &
                                                istate, isub, jstate, Lindx, &
                                                m1, m2, maxcombi
    REAL(real_8) :: bovl2(nstate,nstate), cart_coord(3), DXVnorm, Ebefore, &
      q(4), qmat(3,3), Qnorm, scaled_coord(3), sum, Vnorm
    REAL(real_8), ALLOCATABLE                :: bylm(:,:), rlmat(:,:)

! scratch variables
! cntl%timing variable
! Initialising the cntl%timing for the routine...

    CALL tiset('  main_rot',isub)
    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! Calculating the maximum number of combined LM angular momentum elements
    maxcombi=(dmbi%blmax+1)**2

    ! Allocating the memory for the arrays...
    IF (paral%io_parent)&
         WRITE(6,*) 'Allocating memory'
    ! write(6,*) 'Alloc BCrot'
    ALLOCATE(BCrot(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! write(6,*) 'OK'
    ! write(6,*) 'Alloc BYLM'
    ALLOCATE(bylm(ncpw%ngw,(dmbi%blmax+1)**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! write(6,*) 'OK'
    ! write(6,*) 'Alloc RLMAT'
    ALLOCATE(rlmat((dmbi%blmax+1)**2,(dmbi%blmax+1)**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! write(6,*) 'OK'
    ! write(6,*) 'Alloc STATEmin'
    ALLOCATE(STATEmin(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! write(6,*) 'OK'
    ! write(6,*) 'Alloc STATEmax'
    ALLOCATE(STATEmax(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,*) 'OK'

    ! write(6,*) 'TPIBA = ',tpiba
    ! write(6,*) 'OMEGA = ',OMEGA
    IF (paral%io_parent)&
         WRITE(6,*) 'Lmax = ',dmbi%blmax
    IF (paral%io_parent)&
         WRITE(6,*) 'maxcombi = ',maxcombi
    IF (paral%io_parent)&
         WRITE(6,*) 'NMOL = ',dmbi%nmol
    IF (paral%io_parent)&
         WRITE(6,*)
    ! Checking initialisations

    ! THETA?
    IF (dmbi%in_btheta.EQV..FALSE.) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'No value for theta... assuming 0._real_8'
       ! allocating memory
       ALLOCATE(btheta(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! initialising to zero
       DO imol=1,dmbi%nmol
          btheta(imol)=0._real_8
       ENDDO
    ENDIF

    ! VEC?
    IF (dmbi%in_bvec.EQV..FALSE.) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'No rotation vector... assuming 0._real_8'
       ! allocating memory
       ALLOCATE(BVECx(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(BVECy(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(BVECz(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! initialising to zero
       DO imol=1,dmbi%nmol
          BVECx(imol)=0._real_8
          BVECy(imol)=0._real_8
          BVECz(imol)=0._real_8
       ENDDO
    ENDIF

    ! DXVEC?
    IF (dmbi%in_bdxvec.EQV..FALSE.) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'No displacement vector... assuming 0._real_8'
       ! allocating memory
       ALLOCATE(BDXVECx(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(BDXVECy(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(BDXVECz(dmbi%nmol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! initialising to zero
       !$omp parallel do private(imol)
       DO imol=1,dmbi%nmol
          BDXVECx(imol)=0._real_8
          BDXVECy(imol)=0._real_8
          BDXVECz(imol)=0._real_8
       ENDDO
       !$omp end parallel do
    ENDIF

    ! Scaling the displacement vectors... -------------- TEST
    IF (cntl%tscale) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Scaling DX vectors...'
       CALL dscal(dmbi%nmol,1._real_8/sxscale,BDXVECx,1)
       CALL dscal(dmbi%nmol,1._real_8/syscale,BDXVECy,1)
       CALL dscal(dmbi%nmol,1._real_8/szscale,BDXVECz,1)
       ! now transform the scaled coordinates to cartesian coords 
       ! using the HT matrix (taken from S_to_C ...)
       DO imol=1,dmbi%nmol
          scaled_coord(1)=BDXVECx(imol)
          scaled_coord(2)=BDXVECy(imol)
          scaled_coord(3)=BDXVECz(imol)

          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,scaled_coord,1,&
               0.0_real_8,cart_coord,1)

          BDXVECx(imol)=cart_coord(1)
          BDXVECy(imol)=cart_coord(2)
          BDXVECz(imol)=cart_coord(3)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) 'Scaling done.'
    ENDIF
    ! -------------------------------------------------- TEST

    ! Are we reading a saved original Wannier function?
    IF (dmbi%wann_load) THEN
       ! do we read a distinct reference wannier function for each molecule?
       IF (dmbi%wann_multi) THEN
          ! YES...
          CALL read_multi_wannier(nstate,c0,dmbi%nmol,STATEmin,STATEmax)
       ELSE
          ! NO, one reference for all molecules
          ! GAUSSIAN representation?
          IF (dmbi%wann_gauss) THEN
             ! YES...
             CALL gauss_gen_read(nstate,c0,dmbi%nmol,STATEmin,STATEmax)
          ELSE
             ! No, usual C0 reading...
             CALL read_my_wannier(nstate,c0,dmbi%nmol,STATEmin,STATEmax)
          ENDIF
       ENDIF
    ELSE
       STATEmin(1)=1
       STATEmax(1)=nstate
    ENDIF

    ! Are we saving the original Wannier functions?
    IF (dmbi%wann_save) THEN
       ! Calculating the energy of the starting Wannier function
       IF (paral%io_parent)&
            WRITE(6,*) 'Energy of the starting wannier function'
       CALL calc_wannier_energy_manno(nstate,c0)

       ! GAUSSIAN representation?
       IF (dmbi%wann_gauss) THEN
          ! YES...
          CALL save_my_wannier_gaussian(nstate,c0,ener_com%etot,psi)
       ELSE
          ! No, usual G-space saving
          CALL save_my_wannier(nstate,c0,ener_com%etot)
       ENDIF
    ENDIF

    ! Storing the value of the total energy before rotation
    Ebefore=ener_com%etot

    ! Building the YLM matrix
    IF (paral%io_parent)&
         WRITE(6,*) 'Building the YLM matrix up to Lmax=',dmbi%blmax
    CALL build_bylm_sphgen(dmbi%blmax,bylm)
    IF (paral%io_parent)&
         WRITE(6,*) 'done'

    ! Printing out the real sperical harmonics value for the first 25 G vectors
    ! write(6,*) 'BYLM (for L=2 for 25 1st G vectors):'
    ! Lindx=2
    ! DO M1=1,25 
    ! write(6,FMT='(26(1X,F8.5))')
    !       (BYLM(M1,Lindx**2+Lindx+M2+1),M2=-Lindx,Lindx)
    ! ENDDO

    ! doing the same with the RSPHGEN function
    ! write(6,*) 'Legendre'
    ! Building the YLM matrix
    ! write(6,*) 'Building the YLM matrix up to Lmax=',BLmax
    ! CALL Build_BYLM_SPHGEN(BLmax,BYLM)
    ! write(6,*) 'done'

    ! Printing out the real sperical harmonics value for the first 25 G vectors
    ! write(6,*) 'BYLM (for L=2 for 25 1st G vectors):'
    ! Lindx=2
    ! DO M1=1,25 
    ! write(6,FMT='(26(1X,F8.5))')
    !       (BYLM(M1,Lindx**2+Lindx+M2+1),M2=-Lindx,Lindx)
    ! ENDDO

    ! calculating the overlap between the spherical harmonics
    ! write(6,*) 'OVERLAP S(lm)[YLM(a)*YLM(b)] up to L = ',BLMAX
    ! DO I=1,13
    ! DO J=1,13
    ! sum=0._real_8
    ! DO Lindx=1,maxcombi
    ! sum=sum+BYLM(I,Lindx)*BYLM(J,Lindx)
    ! write(6,*) Lindx,I,BYLM(I,Lindx),J,BYLM(J,Lindx),sum
    ! ENDDO
    ! YLMover(I,J)=sum
    ! ENDDO
    ! write(6,FMT='(13(1X,F8.5))') (YLMover(I,M1),M1=1,13)
    ! ENDDO

    ! Main loop over the molecules
    DO imol=1,dmbi%nmol
       IF (paral%io_parent)&
            WRITE(6,*) '#########################'
       IF (paral%io_parent)&
            WRITE(6,'(1x,a,i6,a)') '# DOING MOLECULE ',imol,' #'
       IF (paral%io_parent)&
            WRITE(6,*) '#########################'
       IF (paral%io_parent)&
            WRITE(6,*)
       ! Constructing the rotation quaternion
       Vnorm=SQRT(BVECx(imol)**2+BVECy(imol)**2+BVECz(imol)**2)
       IF (Vnorm.GT.0.0_real_8) THEN
          q(1)=COS(btheta(imol)/2)
          q(2)=SIN(btheta(imol)/2)*BVECx(imol)/Vnorm
          q(3)=SIN(btheta(imol)/2)*BVECy(imol)/Vnorm
          q(4)=SIN(btheta(imol)/2)*BVECz(imol)/Vnorm
          ! if direction vector is no good then no rotation
       ELSE
          q(1)=1
          q(2)=0
          q(3)=0
          q(4)=0
       ENDIF

       ! Printing out rotation parameters
       IF (paral%io_parent)&
            WRITE(6,*) 'Rotation parameters:'
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,F8.5,A,3(F8.5,A))') 'Theta= ',&
            BTHETA(imol),' VEC: (',BVECx(imol),',',BVECy(imol),&
            ',',BVECz(imol),' )'
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,4(F8.5,A))') 'q = (',q(1),', ',q(2),&
            ', ',q(3),', ',q(4),' )'
       Qnorm=SQRT(q(1)**2+q(2)**2+q(3)**2+q(4)**2)
       IF (paral%io_parent)&
            WRITE(6,*) 'Qnorm = ',qnorm
       IF (paral%io_parent)&
            WRITE(6,*)

       ! Displacement vector
       IF (paral%io_parent)&
            WRITE(6,*) 'Displacement parameters:'
       IF (paral%io_parent)&
            WRITE(6,fmt='(1X,A,3(F8.5,A))') 'DXVEC: (',BDXVECx(imol),&
            ',',BDXVECy(imol),',',BDXVECz(imol),' )'

       ! Calculating norm of displacement vector (if ==0, no shift..)
       DXVnorm=SQRT(BDXVECx(imol)**2+BDXVECy(imol)**2&
            +BDXVECz(imol)**2)
       IF (paral%io_parent)&
            WRITE(6,*) 'DXVnorm = ',dxvnorm
       IF (paral%io_parent)&
            WRITE(6,*)

       ! If all the rotations are performed (even rot by 0 deg)
       IF (dmbi%wann_allrot.OR.btheta(imol).NE.0._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Performing rotation for molecule',imol
          ! just testing the rotation matrices...
          ! CALL buildqmat(qmat,q)
          CALL buildqmatsgen(qmat,q)
          IF (paral%io_parent)&
               WRITE(6,*) 'Matrix representation of the quaternion'
          DO i=1,3
             IF (paral%io_parent)&
                  WRITE(6,fmt='(3(1X,F8.5))') qmat(i,1),qmat(i,2),qmat(i,3)
          ENDDO

          ! initialising the rotation matrx           
          CALL build_rmat(dmbi%blmax,rlmat,q)

          ! printing out all the rotation matrices by increasing L
          IF (paral%io_parent)&
               WRITE(6,*) 'RLMAT:'
          DO Lindx=0,dmbi%blmax
             IF (paral%io_parent)&
                  WRITE(6,*) 'L = ',Lindx
             DO m1=-Lindx,lindx
                combi=Lindx**2+lindx+m1+1
                IF (paral%io_parent)&
                     WRITE(6,fmt='(13(1X,F8.5))') (rlmat(combi,&
                     Lindx**2+Lindx+M2+1),M2=-Lindx,Lindx)
             ENDDO
          ENDDO

          ! Now performing the rotation of the wannier functions
          IF (paral%io_parent)&
               WRITE(6,*) 'Performing the rotation'
          CALL rotatecoeff(STATEmin(imol),STATEmax(imol),nstate,dmbi%blmax,&
               RLMAT,BYLM,C0,BCrot)
          IF (paral%io_parent)&
               WRITE(6,*) 'FINISHED'
       ELSE
          ! replacing the rotation by zero by a simple copy of C0 to BCrot
          IF (paral%io_parent)&
               WRITE(6,*) 'No rotations for molecule',imol,&
               ', copying the coefficients..'
          DO istate=STATEmin(imol),STATEmax(imol)
             DO indx=1,ncpw%ngw
                BCrot(indx,istate)=c0(indx,istate)
             ENDDO
          ENDDO
       ENDIF

       ! Now we can shift the Wannier functions by the DX vector
       ! if its norm is greater than zero
       IF (DXVnorm.GT.0._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Shifting the Wannier function'
          CALL shift_my_wannier_manno(STATEmin(imol),STATEmax(imol),&
               NSTATE,BCrot,BDXVECx(imol),BDXVECy(imol),BDXVECz(imol))
          IF (paral%io_parent)&
               WRITE(6,*) 'FINISHED'
       ENDIF
    ENDDO

    ! NOW: all the wannier functions are rotated and moved according to the 
    ! data of the input file...

    ! if requested orthogonalise the ROTATED states+renormalising...
    IF (dmbi%torthog_wannier.AND..NOT.dmbi%tsimple_model) THEN
       ! IF (TORTHOG_WANNIER.AND.TSIMPLE_MODEL) THEN
       ! if we use simple_model the the orthogonalisation is done later..
       CALL ortho_my_wannier(nstate,BCrot)
    ENDIF

    ! Overlap matrix for the rotated states
    IF (paral%io_parent)&
         WRITE(6,*) 'Building overlap matrix for the ROTATED states'
    IF ((dmbi%torthog_wannier).AND.paral%io_parent)&
         WRITE(6,*) 'after orthogonalisation'
    CALL ovlap(nstate,bovl2,BCrot,bcrot)
    IF (paral%io_parent)&
         WRITE(6,*) 'Checking overlap matrix...'
    DO istate=1,nstate
       sum=0._real_8
       DO jstate=1,nstate
          ! checking matrix symmetry...
          IF (bovl2(istate,jstate).NE.bovl2(jstate,istate)) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Error in matrix symm. for states :',istate,&
                  'and ',jstate
          ENDIF
          ! checking sum of off-diagonal elements
          sum=sum+bovl2(istate,jstate)
       ENDDO

       IF (ABS(bovl2(istate,istate)- 1._real_8).GT.1.e-10_real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Error of normalisation for state ',istate
       ENDIF

       IF (ABS(sum-bovl2(istate,istate)).GT.1.e-10_real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Detected overlap for state ',istate
          IF (paral%io_parent)&
               WRITE(6,*) 'sum of off-diag elements = ',sum
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) 'Finished'

    IF (nstate.LT.10) THEN
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,11(1X,I12))') 'ISTATE',(jstate,jstate=1,nstate)
       DO istate=1,nstate
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,11(1X,F12.6))') istate,(bovl2(istate,jstate)&
               ,jstate=1,NSTATE)
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*)
    ! Copying the rotated Wannier functions to the original C0
    IF (paral%io_parent)&
         WRITE(6,*) 'Now replacing original by rotated C0'
    CALL dcopy(2*ncpw%ngw*nstate,BCrot,1,c0,1)

    ! Now calculating the energy of the rotated wannier funct.
    ! write(6,*) 'Energy of the ROTATED Wannier system'
    ! CALL Calc_Wannier_Energy_manno(NSTATE,C0)

    ! Storing ETOT after rotation
    ! Eafter=ETOT
    ! Ediff=Ebefore-Eafter
    ! write(6,*) 'Energy difference =',Ediff

    ! Stop the clock!!
    CALL tihalt('  main_rot',isub)

    ! Freeing the memory used by the routine
    DEALLOCATE(rlmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bylm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(BCrot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! CALL FREEM(IP_STATEmin)
    ! CALL FREEM(IP_STATEmax)
  END SUBROUTINE rotatemywannier_manno

  ! ==================================================================
  SUBROUTINE save_my_wannier(nstate,c0,etot)
    ! ==--------------------------------------------------------------==
    ! This subroutine saves the wannier functions in the file   
    ! "ref-wannier.dat". Upon completion the file will contain:
    ! ETOT                                                     
    ! NSTATE                                                   
    ! NGW                                                      
    ! C0(1..NGW,1..NSTATE)                                     
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: etot

    INTEGER                                  :: Gindx, istate
    LOGICAL                                  :: ferror

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! opening the file
    IF (paral%io_parent)&
         CALL fileopen(500,'ref-wannier.dat',fo_def,ferror)
    ! Did the file open OK?
    IF (ferror) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Error creating the file ref-wannier.dat...'
       IF (paral%io_parent)&
            WRITE(6,*) 'The file probably exists already! '
       CALL stopgm('Save_My_Wannier','file open error',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*) 'Saving ',nstate,' states...'
    IF (paral%io_parent)&
         WRITE(6,*) 'Total energy for this set Etot = ',etot

    ! Saving data to file..
    IF (paral%io_parent)&
         WRITE(500,*) etot
    IF (paral%io_parent)&
         WRITE(500,*) nstate
    IF (paral%io_parent)&
         WRITE(500,*) ncpw%ngw

    ! Main loop on the states
    DO istate=1,nstate
       ! Loop on the G vectors
       DO Gindx=1,ncpw%ngw
          IF (paral%io_parent)&
               WRITE(500,*) c0(Gindx,istate)
       ENDDO
    ENDDO

    ! Closing the file...
    IF (paral%io_parent)&
         CALL fileclose(500)

  END SUBROUTINE save_my_wannier

  ! ==================================================================
  SUBROUTINE read_my_wannier(nstate,c0,nmol,STATEmin,STATEmax)
    ! ==--------------------------------------------------------------==
    ! This subroutine reads the wannier functions saved in    
    ! "ref-wannier.dat", and build up the corresponding states
    ! for each molecules...                                   
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: nmol, STATEmin(nmol), &
                                                STATEmax(nmol)

    CHARACTER(*), PARAMETER                  :: procedureN = 'read_my_wannier'

    COMPLEX(real_8), ALLOCATABLE             :: Ct(:,:)
    INTEGER                                  :: Gindx, ierr, istate, kstate, &
                                                mol, r_ngw, r_nstate, &
                                                states_to_read
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: etot

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! opening the file
    IF (paral%io_parent)&
         CALL fileopen(500,'ref-wannier.dat',fo_old,ferror)
    ! Did the file open OK?
    IF (ferror) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Error:the file ref-wannier.dat is not there! '
       CALL stopgm('Read_My_Wannier','file open error',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*) 'Reading data from file...'

    ! Reading index data from file..
    IF (paral%io_parent)&
         READ(500,*) etot
    IF (paral%io_parent)&
         READ(500,*) r_nstate
    IF (paral%io_parent)&
         READ(500,*) r_ngw
    IF (paral%io_parent)&
         WRITE(6,*) 'The file contains ',r_nstate,' states'
    IF (paral%io_parent)&
         WRITE(6,*) 'Each state has ',r_ngw,' G vectors'
    IF (paral%io_parent)&
         WRITE(6,*) 'Total energy of the saved state Etot = ',etot
    ! Checking the number of G vectors...
    IF (r_ngw.NE.ncpw%ngw) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'The number of G vectors found in the file (',&
            R_NGW,') do not match NGW in the simulation (',&
            ncpw%ngw,')'
       IF (paral%io_parent)&
            WRITE(6,*) 'At present this is not allowed...aborting load'
       RETURN
    ENDIF

    ! Calculating the number of states to read
    states_to_read=nstate/nmol
    IF (paral%io_parent)&
         WRITE(6,*) 'Number of states per molecule ',states_to_read

    ! checking the number of states
    IF (states_to_read.NE.r_nstate) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Need more states for each molecules...exiting'
       RETURN
    ENDIF

    ! Everything is OK...allocating the Ct array (2*dims for CPLX*16)

    ALLOCATE(Ct(ncpw%ngw,r_ngw*r_nstate/ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! reading the wannier functions
    IF (paral%io_parent)&
         WRITE(6,*) 'Proceeding... loading the data for ',nmol,&
         ' molecule(s)'

    ! Main loop on the read states
    DO istate=1,r_nstate
       ! Loop on the read G vectors
       DO Gindx=1,r_ngw
          IF (paral%io_parent)&
               READ(500,*) Ct(Gindx,istate)
       ENDDO
    ENDDO

    ! Closing the file...
    IF (paral%io_parent)&
         CALL fileclose(500)
    IF (paral%io_parent)&
         WRITE(6,*) 'Building',nstate,' Wannier states'
    ! Now filling up the C0 array with the required ammount of data
    ! initialising the global state counter
    kstate=1
    ! Loop on number of molecules
    DO mol=1,nmol
       STATEmin(mol)=kstate
       ! Loop on number of read states
       DO istate=1,r_nstate
          ! Loop on Gvecs
          DO Gindx=1,r_ngw
             c0(Gindx,kstate)=Ct(gindx,istate)
          ENDDO
          ! updating the global state counter
          kstate=kstate+1
       ENDDO
       STATEmax(mol)=kstate-1
    ENDDO

    DO mol=1,nmol
       IF (paral%io_parent)&
            WRITE(6,*) 'molecule: ',mol,' STATEmin = ',&
            STATEmin(mol),' STATEmax = ',STATEmax(mol)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) 'Done! '
    ! Freeing memory
    DEALLOCATE(Ct,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

  END SUBROUTINE read_my_wannier

  ! ==================================================================
  SUBROUTINE read_multi_wannier(nstate,c0,nmol,STATEmin,STATEmax)
    ! ==--------------------------------------------------------------==
    ! This subroutine reads a wannier functions for each molecule C
    ! saved in "ref-wannier-1.dat", "ref-wannier-2.dat", etc..    C
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: nmol, STATEmin(nmol), &
                                                STATEmax(nmol)

    CHARACTER(len=20)                        :: name
    CHARACTER(len=8)                         :: number
    INTEGER                                  :: Gindx, ia, ie, istate, &
                                                kstate, mol, r_ngw, r_nstate
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: etot

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,*) 'Multi wannier reference file reading subroutine'
    IF (paral%io_parent)&
         WRITE(6,*) 'this will read ',nmol,' .dat files'

    ! initialising the global state counter
    kstate=1
    ! main loop on the molecules...
    DO mol=1,nmol

       ! First creating the file name for the molecule mol
       IF (paral%io_parent)&
            WRITE(number,'(I8)')mol
       ! removing the spaces before and after the number (CPMD routine)
       CALL xstring(number,ia,ie)

       name='ref-wannier-'//number(ia:ie)//'.dat'
       IF (paral%io_parent)&
            WRITE(6,*) 'molecule ',mol,'opening file ',name

       ! opening the file
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(500,name,fo_old,ferror)
       ! Did the file open OK?
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Error:the file',name,'is not there! '
          CALL stopgm('READ_MULT_WANN','NO FILE! ',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*) 'Reading data from file...'

       ! Reading index data from file..
       IF (paral%io_parent)&
            READ(500,*) etot
       IF (paral%io_parent)&
            READ(500,*) r_nstate
       IF (paral%io_parent)&
            READ(500,*) r_ngw
       IF (paral%io_parent)&
            WRITE(6,*) 'The file contains ',r_nstate,' states'
       IF (paral%io_parent)&
            WRITE(6,*) 'Each state has ',r_ngw,' G vectors'
       IF (paral%io_parent)&
            WRITE(6,*) 'Total energy of the saved state Etot = ',etot

       ! Checking the number of G vectors...
       IF (r_ngw.NE.ncpw%ngw) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'The number of G vectors found in the file (',&
               R_NGW,') do not match NGW in the simulation (',&
               ncpw%ngw,')'
          IF (paral%io_parent)&
               WRITE(6,*) 'At present this is not allowed...aborting load'
          RETURN
       ENDIF

       ! first state occupied by molecule mol
       STATEmin(mol)=kstate
       ! Main loop on the read states
       DO istate=1,r_nstate
          ! Loop on the read G vectors
          DO Gindx=1,r_ngw
             IF (paral%io_parent)&
                  READ(500,*) c0(Gindx,STATEmin(mol)+istate-1)
          ENDDO
          ! updating the global state counter
          kstate=kstate+1
       ENDDO

       ! final state occupied by molecule mol
       STATEmax(mol)=kstate-1
       ! Closing the file...
       IF (paral%io_parent)&
            CALL fileclose(500)
    ENDDO

    DO mol=1,nmol
       IF (paral%io_parent)&
            WRITE(6,*) 'molecule: ',mol,' STATEmin = ',&
            STATEmin(mol),' STATEmax = ',STATEmax(mol)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) 'Done! '

  END SUBROUTINE read_multi_wannier

  ! ==================================================================
  SUBROUTINE shift_my_wannier_manno(st_state,fi_state,nstate,&
       C0,Vx,Vy,Vz)
    ! ==--------------------------------------------------------------==
    ! This subroutine shifts the wannier function by a vector
    ! defined by Vx,Vy,Vz.                                   
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: st_state, fi_state, nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: Vx, Vy, Vz

    COMPLEX(real_8)                          :: EXPiGR
    INTEGER                                  :: Gindx, istate
    REAL(real_8)                             :: GbyR

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,*) 'DOING STATE ',st_state,' TO ',fi_state

    ! Main loop on the states
    DO istate=st_state,fi_state
       ! Loop on the G vectors
       DO Gindx=1,ncpw%ngw
          ! Calculating the product G.r
          GbyR=gk(1,Gindx)*Vx+gk(2,gindx)*Vy+gk(3,gindx)*Vz
          ! write(6,*) 'GbyR = ',GbyR

          ! Adjusting the units
          GbyR=-gbyr*parm%tpiba
          ! write(6,*) 'GbyR = ',GbyR

          ! Calculating E^i(G.r) using cos(x)+i*sin(x)
          EXPiGR=CMPLX(COS(GbyR),SIN(gbyr),kind=real_8)
          ! write(6,*) 'EXPiGR = ',EXPiGR

          ! multiplying the Wannier coeffs
          c0(Gindx,istate)=c0(gindx,istate)*EXPiGR
       ENDDO
    ENDDO
  END SUBROUTINE shift_my_wannier_manno

  ! ==================================================================
  SUBROUTINE ortho_my_wannier(nstate,Ct)
    ! ==--------------------------------------------------------------==
    ! This subroutine orthogonalises the Wannier functions using  C
    ! the lowdin method (PSI)i=sum(j)[Sij^-(1/2) * (PSI)j]        C
    ! also called symmetric orthogonalisation method (Szabo p142) C
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: Ct(ncpw%ngw,nstate)

    COMPLEX(real_8)                          :: tmp
    INTEGER                                  :: Gindx, istate, jstate, k
    COMPLEX(real_8)                          :: tmpA(nstate)
    REAL(real_8) :: eigenval(nstate), eigenvecs(nstate,nstate), &
      EVMatrix(nstate,nstate), EVMtmp(nstate,nstate), s(nstate,nstate), &
      work(3*nstate), z(nstate*(nstate+1)/2)

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,*) 'ORTHOGONALISATION, NSTATE = ',nstate

    ! initialising counter for z matrix
    k=1
    ! building up the S matrix(overlap)
    DO istate=1,nstate
       DO jstate=istate,nstate
          tmp=(0._real_8,0._real_8)
          ! Calculating the overlap of the positive g vect coeffs
          DO Gindx=2,ncpw%ngw
             tmp=tmp+CONJG(Ct(Gindx,istate))*ct(gindx,jstate)
          ENDDO
          ! Special case for G=0 and adding the other G vecs
          s(istate,jstate)&
               =REAL(CONJG(Ct(1,istate)*Ct(1,jstate)+2._real_8*tmp),kind=real_8)
          ! Symmetrising the matrix
          s(jstate,istate)=s(istate,jstate)
          ! Filling the special column Z matrix for dspevy
          z(k)=s(istate,jstate)
          k=k+1
       ENDDO
    ENDDO

    ! if the number of states is no too big
    IF (nstate.LT.10) THEN
       ! printing overlap matrix
       IF (paral%io_parent)&
            WRITE(6,*) 'orthog: Starting S matrix'
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,11(1X,I12))') 'ISTATE',&
            (istate,istate=1,NSTATE)
       DO istate=1,nstate
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,11(1X,F12.6))') istate,(s(istate,jstate)&
               ,jstate=1,NSTATE)
       ENDDO
    ENDIF

    ! Diagonalising the S matrix through the z column matrix
    CALL dspevy(1,z,eigenval,eigenvecs,nstate,nstate,work,3*nstate)

    ! if the number of states is no too big
    IF (nstate.LT.10) THEN
       ! printing eigenval and eigenvecs
       IF (paral%io_parent)&
            WRITE(6,*) 'orthog: eigenvalues + eigenvecs'
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,A10,10(1X,I12))') 'ISTATE',' EIGENVAL',&
            (jstate,jstate=1,NSTATE)
       DO istate=1,nstate
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,F10.6,10(1X,F12.6))') istate,&
               eigenval(istate),&
               (eigenvecs(istate,jstate),jstate=1,NSTATE)
       ENDDO
    ENDIF

    ! Preparing eigenvalue matrix for S^-1/2 
    DO istate=1,nstate
       ! checking if all the eigenvalues are greater than 0
       IF (eigenval(istate).LE. 0._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'linear dependency between orbitals detected at ',&
               istate
       ENDIF
       DO jstate=1,nstate
          IF (istate.EQ.jstate) THEN
             EVMatrix(istate,istate)=1._real_8/SQRT(eigenval(istate))
          ELSE
             EVMatrix(istate,jstate)=0._real_8
          ENDIF
       ENDDO
    ENDDO

    ! building the S^-1/2 matrix (evec.1/sqrt(eval).Tevec)
    ! First EVMtmp=1/sqrt(eval).Tevec
    CALL dgemm('N','T',nstate,nstate,nstate,1._real_8,EVMatrix,&
         NSTATE,eigenvecs,NSTATE,0._real_8,EVMtmp,NSTATE)

    ! Now replacing S by evec.EVMtmp
    CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,eigenvecs,&
         NSTATE,EVMtmp,NSTATE,0._real_8,S,NSTATE)

    ! if the number of states is no too big
    IF (nstate.LT.10) THEN
       ! printing final transformation matrix
       IF (paral%io_parent)&
            WRITE(6,*) 'orthog: FINAL S^-1/2 matrix'
       IF (paral%io_parent)&
            WRITE(6,fmt='(A7,11(1X,I12))') 'ISTATE',&
            (istate,istate=1,NSTATE)
       DO istate=1,nstate
          IF (paral%io_parent)&
               WRITE(6,fmt='(I7,11(1X,F12.6))') istate,(s(istate,jstate)&
               ,jstate=1,NSTATE)
       ENDDO
    ENDIF
    ! Finally, performing the orthogonalisation transformation:
    ! (PSI)i=sum(j)[Sij^-(1/2) * (PSI)j]        

    DO Gindx=1,ncpw%ngw
       DO istate=1,nstate
          tmpA(istate)=(0._real_8,0._real_8)
          DO jstate=1,nstate
             tmpA(istate)=tmpa(istate)+&
                  (S(istate,jstate)*Ct(Gindx,jstate))
          ENDDO
       ENDDO
       ! copying the orthogonalised coefficients
       DO istate=1,nstate
          Ct(Gindx,istate)=tmpA(istate)
       ENDDO
    ENDDO

  END SUBROUTINE ortho_my_wannier

  ! ==================================================================
  SUBROUTINE build_bylm_sphgen(Lmax,bylm)
    ! ==--------------------------------------------------------------==
    ! This subroutine calculates the BYLM matrix using RSPHGEN!
    ! BYLM is a (NGW)x(Lmax+1)^2 matrix containing the value of the 
    ! real sperical harmonics for each G vector.                            
    ! ==--------------------------------------------------------------==

    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)


    INTEGER                                  :: Lmax
    REAL(real_8)                             :: bylm(ncpw%ngw,(Lmax+1)**2)

    INTEGER                                  :: combi, gwindx, isub, Lbase, &
                                                Lindx, Mindx

! cntl%timing Variable

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! Initialising the cntl%timing for the routine...
    CALL tiset('B_SGENBYLM',isub)

    ! Starting the loading loops
    DO Lindx=0,Lmax
       ! Lbase is a time saving device to calculate the LM combined index
       ! given by LM=L^2+L+M+1
       Lbase=Lindx**2+lindx+1
       DO Mindx=-Lindx,lindx
          combi=Lbase+Mindx
          ! Loading up the values of the real sperical harmonics for a
          ! particular COMBINED index for each G vector
          DO gwindx=1,ncpw%ngw
             bylm(gwindx,combi)=rsphgen(Lindx,Mindx,&
                  GK(1,gwindx),GK(2,gwindx),GK(3,gwindx))

             bylm(gwindx,combi)&
                  =BYLM(gwindx,combi)/SQRT(REAL(Lmax+1,kind=real_8))
          ENDDO
       ENDDO
    ENDDO

    ! Stop the clock!
    CALL tihalt('B_SGENBYLM',isub)

  END SUBROUTINE build_bylm_sphgen

  ! ==================================================================
  SUBROUTINE build_rmat(Lmax,rlmat,q)
    ! ==--------------------------------------------------------------==
    ! This subroutine calculates the RLMAT matrix
    ! RLMAT is a (Lmax+1)^2x(Lmax+1)^2 matrix for the rotation of
    ! real spherical harmonics by a given quaternion q.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: Lmax
    REAL(real_8)                             :: rlmat((Lmax+1)**2,(lmax+1)**2)&
                                                , q(4)

    INTEGER                                  :: combi, combiT, isub, Lbase, &
                                                Lindx, Mindx, Mtindx
    REAL(real_8)                             :: qmat(3,3)

! cntl%timing Variable

    CALL stopgm('ROTATE_..._P','NOT PARALLEL.',& 
         __LINE__,__FILE__)

    ! Initialising the cntl%timing for the routine...
    CALL tiset('Build_RMAT',isub)

    ! building the initial 3x3 rotation matrix from the quaternion q
    ! CALL buildqmat(QMAT,q)
    CALL buildqmatsgen(qmat,q)

    ! Starting the loading loops
    DO Lindx=0,Lmax
       ! Lbase is a time saving device to calculate the LM combined index
       ! given by LM=L^2+L+M+1
       Lbase=Lindx**2+lindx+1
       DO Mindx=-Lindx,lindx
          combi=Lbase+Mindx
          ! Building the rotation matrix for each COMBINED index
          DO Mtindx=-Lindx,lindx
             combiT=Lbase+Mtindx
             rlmat(combiT,combi)&
                  =Rsph(Lindx,Mtindx,Mindx,QMAT,RLmat,Lmax)
          ENDDO
       ENDDO
    ENDDO

    ! Stop the clock!
    CALL tihalt('Build_RMAT',isub)

  END SUBROUTINE build_rmat

  ! ==================================================================
  SUBROUTINE rotatecoeff(st_state,fi_state,nstate,Lmax,rlmat,bylm,&
       C0,BCrot)
    ! ==--------------------------------------------------------------==
    ! This subroutine rotates a given set of coefficients (C0) of G
    ! vectors using a pre-defined rotation matrix (RLMAT) and pre-defined
    ! G vector sperical harmonic values matrix (BYLM), the resulting
    ! coefficients are found in BCrot
    ! ==--------------------------------------------------------------==
    ! Including the "essential" variables...

    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: st_state, fi_state, nstate, &
                                                Lmax
    REAL(real_8)                             :: rlmat((Lmax+1)**2,(lmax+1)**2)&
                                                , bylm(ncpw%ngw,(Lmax+1)**2)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), &
                                                BCrot(ncpw%ngw,nstate)

    COMPLEX(real_8)                          :: coeff
    INTEGER                                  :: combi, combiT, GFindx, GFmax, &
                                                GFmin, GKindx, istate, isub, &
                                                Lbase, Lindx, Mindx, Mtindx
    REAL(real_8)                             :: sumlk, Uneg, Upos, Uvar

! cntl%timing variable
! real(8) :: prefac
! PARAMETER(prefac=8.0_real_8*3.141592653589793e+0_real_8**3)

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! Initialising the cntl%timing for the routine...
    CALL tiset(' rot_the_C0',isub)

    ! write(6,*) 'prefac =',prefac

    ! Main loop
    DO istate=st_state,fi_state
       ! Initialising GFmin and GFmax...
       GFmin=0
       GFmax=0
       IF (paral%io_parent)&
            WRITE(6,*) 'doing state ',istate
       ! Special case for G(1)
       BCrot(1,istate)=c0(1,istate)
       ! Then proceed with G from 2 to NGW
       DO GFindx=2,ncpw%ngw
          coeff=(0._real_8,0._real_8)
          ! GFmax has to be set to the last G vec having a norm equal to GFmin
          IF (GFindx.GT.GFmax) THEN
             GFmin=GFindx
             GFmax=GFindx
             ! do while((HG(GFmax+1).le.(HG(GFmin)+1)).and.(GFmax.lt.NGW))
             DO WHILE((hg(GFmax+1).EQ.hg(GFmin)).AND.(gfmax.LT.ncpw%ngw))
                GFmax=gfmax+1
             ENDDO
          ENDIF
          ! Now GFmin is the index of the first Gvector and GFmax 
          ! is the index of the last G vector having the same norm
          DO GKindx=GFmin,GFmax
             ! Initialising the sums over Uvar for positive and negative G vectors
             Upos=0._real_8
             Uneg=0._real_8
             DO Lindx=0,Lmax
                Uvar=0._real_8
                ! Lbase is a time saving device to calculate the LM combined index
                ! given by LM=L^2+L+M+1
                Lbase=Lindx**2+lindx+1
                DO Mtindx=-Lindx,lindx
                   combiT=Lbase+Mtindx
                   sumlk=0._real_8
                   DO Mindx=-Lindx,lindx
                      combi=Lbase+Mindx
                      ! Calculating the sum over Mindx of the spherical 
                      ! harmonic for Gk by rotation matrix element
                      sumlk=sumlk+bylm(GKindx,combi)&
                           *RLMAT(combiT,combi)
                   ENDDO
                   ! Now multiplying the result by the spherical
                   ! harmonic of Gf
                   Uvar=uvar+sumlk*bylm(GFindx,combiT)
                ENDDO
                ! Adding the positive and negative contributions
                Upos=upos+Uvar
                Uneg=uneg+(-1)**Lindx*Uvar
             ENDDO
             ! Summing over the positive and negative Gvector coefficients
             coeff=coeff+c0(GKindx,istate)*Upos&
                  +CONJG(C0(GKindx,istate))*Uneg
          ENDDO
          ! Finally multiplicating the result by required pre-factors
          BCrot(GFindx,istate)=coeff
          ! BCrot(GFindx,istate)=(prefac)*coeff
          ! BCrot(GFindx,istate)=(prefac/HG(GFindx))*coeff
       ENDDO
    ENDDO

    ! Stop the clock...
    CALL tihalt(' rot_the_C0',isub)

  END SUBROUTINE rotatecoeff

  ! ==================================================================
  INTEGER FUNCTION IDelta(a)
    ! ==--------------------------------------------------------------==
    ! Integer Dirac delta function with integer argument
    ! ==--------------------------------------------------------------==
    INTEGER :: a

    IDelta=1
    IF (a.NE.0) THEN
       IDelta=0
    ENDIF
  END FUNCTION IDelta
  ! ==================================================================
  SUBROUTINE buildqmat(qmat,q)
    ! ==--------------------------------------------------------------==
    ! Subroutine builds the modified rotation matrix corresponding to a
    ! given unit quaternion q
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: qmat(3,3), q(4)

    qmat(1,1)=q(1)**2+q(2)**2-q(3)**2-q(4)**2
    qmat(2,1)=-2._real_8 *(q(1)*q(3)+q(2)*q(4))
    qmat(3,1)=-2._real_8 *(q(2)*q(3)-q(1)*q(4))
    qmat(1,2)=2._real_8 *(q(1)*q(3)-q(2)*q(4))
    qmat(2,2)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
    qmat(3,2)=2._real_8 *(q(1)*q(2)+q(3)*q(4))
    qmat(1,3)=- 2._real_8 *(q(2)*q(3)+q(1)*q(4))
    qmat(2,3)=-2._real_8 *(q(1)*q(2)-q(3)*q(4))
    qmat(3,3)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
  END SUBROUTINE buildqmat
  ! ==================================================================
  SUBROUTINE buildqmatsgen(qmat,q)
    ! ==--------------------------------------------------------------==
    ! Subroutine builds the modified rotation matrix corresponding to a
    ! given unit quaternion q to use with RSPHGEN spherical harmonics
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: qmat(3,3), q(4)

    qmat(1,1)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
    qmat(2,1)=2._real_8*(-(q(1)*q(2))-q(3)*q(4))
    qmat(3,1)=2._real_8*(-(q(2)*q(3))+q(1)*q(4))
    qmat(1,2)=2._real_8*(q(1)*q(2)-q(3)*q(4))
    qmat(2,2)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
    qmat(3,2)=2._real_8*(q(1)*q(3)+q(2)*q(4))
    qmat(1,3)=2._real_8*(-(q(2)*q(3))-q(1)*q(4))
    qmat(2,3)=2._real_8*(-(q(1)*q(3))+q(2)*q(4))
    qmat(3,3)=q(1)**2+q(2)**2-q(3)**2-q(4)**2

  END SUBROUTINE buildqmatsgen
  ! ==================================================================
  REAL(real_8) FUNCTION Rsph(l,m1,m2,qmat,rlmat,Lmax)
    ! ==--------------------------------------------------------------==
    ! General rotation matrix for real spherical harmonics 
    ! J. Phys. Chem. 100, 6342 (1996) / J. Phys. Chem. A 102, 9099 (1998)
    ! 24/07/00 dmb: removed the recursion relationship to speed-up the code
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2,Lmax
    REAL(real_8) :: qmat(3,3)
    REAL(real_8) :: rlmat((Lmax+1)**2,(lmax+1)**2)
    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    IF ((l.GT.1) .AND. (ABS(m1).LE.l) .AND. (ABS(m2).LE.l)) THEN
       Rsph=smallu(l,m1,m2)*Pfunc(0,l,m1,m2,qmat,rlmat,Lmax) +&
            smallv(L,M1,M2)*bigV(L,M1,M2,qmat,RLmat,Lmax) +&
            smallw(L,M1,M2)*bigW(L,M1,M2,qmat,RLmat,Lmax)
       ! note that the bigU function is equivalent to the P function
    ELSE
       Rsph=0
    ENDIF
    IF (l.EQ.1) THEN
       IF (ABS(m1).LT.2 .AND. ABS(m2).LT.2) THEN
          Rsph=qmat(m1+2,m2+2)
       ELSE
          Rsph=0
       ENDIF
    ENDIF
    IF (l.EQ.0) THEN
       Rsph=1
    ENDIF
  END FUNCTION Rsph
  ! ==================================================================
  REAL(real_8) FUNCTION smallu(l,m1,m2)
    ! ==--------------------------------------------------------------==
    ! Small u function 
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2

    IF (ABS(m2).LT.l) THEN
       smallu=SQRT(REAL((l+m1)*(l-m1),kind=real_8)/REAL((l+m2)*(l-m2),kind=real_8))
    ENDIF
    IF (ABS(m2).EQ.l) THEN
       smallu=SQRT(REAL((l+m1)*(l-m1),kind=real_8)/REAL(2*l*(2*l-1),kind=real_8))
    ENDIF
  END FUNCTION smallu
  ! ==================================================================
  REAL(real_8) FUNCTION smallv(l,m1,m2)
    ! ==--------------------------------------------------------------==
    ! Small v function
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2

    ! write(6,*) 'received L =',L,' M1 =',M1,' M2 =',M2
    IF (ABS(m2).LT.l) THEN
       smallv=0.5_real_8*SQRT(REAL((1+IDelta(m1))*(l+ABS(m1)-1)*&
            (L+ABS(M1)),kind=real_8)/REAL((L+M2)*(L-M2),kind=real_8))&
            *REAL(1-2*IDelta(M1),kind=real_8)
    ENDIF
    IF (ABS(m2).EQ.l) THEN
       smallv=0.5_real_8*SQRT(REAL((1+IDelta(m1))*(l+ABS(m1)-1)*&
            (L+ABS(M1)),kind=real_8)/REAL(2*L*(2*L-1),kind=real_8))&
            *REAL(1-2*IDelta(M1),kind=real_8)
    ENDIF
  END FUNCTION smallv
  ! ==================================================================
  REAL(real_8) FUNCTION smallw(l,m1,m2)
    ! ==--------------------------------------------------------------==
    ! Small w function
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2

    IF (ABS(m2).LT.l) THEN
       smallw=-0.5_real_8*SQRT(REAL((l-ABS(m1)-1)*(l-ABS(m1)),kind=real_8)/&
            REAL((L+M2)*(L-M2),kind=real_8))*REAL(1-IDelta(M1),kind=real_8)
    ENDIF
    IF (ABS(m2).EQ.l) THEN
       smallw=-0.5_real_8*SQRT(REAL((l-ABS(m1)-1)*(l-ABS(m1)),kind=real_8)/&
            REAL(2*L*(2*L-1),kind=real_8))*REAL(1-IDelta(M1),kind=real_8)
    ENDIF
  END FUNCTION smallw
  ! ==================================================================
  REAL(real_8) FUNCTION Pfunc(k,l,Mu,m2,qmat,rlmat,Lmax)
    ! ==--------------------------------------------------------------==
    ! P function (no conversions)   C
    ! ==--------------------------------------------------------------==
    INTEGER :: k,l,Mu,m2,Lmax
    INTEGER :: combiLMu,combi2,combi3
    REAL(real_8) :: qmat(3,3)
    REAL(real_8) :: rlmat((Lmax+1)**2,(lmax+1)**2)

    ! Calculating the common combined index for the L-1,Mu element of RLmat
    combiLMu=(l-1)**2+(l-1)+Mu+1

    ! write(6,*) 'Pfunc received k =',k,' L =',L,' Mu =',Mu,' M2 =',M2
    IF (ABS(m2).LT.l) THEN
       ! calculating the combined indices for RLmat        
       combi2=(l-1)**2+(l-1)+m2+1
       ! now the value of the Pfunction
       Pfunc=qmat(k+2,2)*rlmat(combiLMu,combi2)
    ENDIF

    IF (m2.EQ.l) THEN
       ! calculating the combined indices for RLmat        
       combi2=(l-1)**2+(l-1)+(l-1)+1
       combi3=(l-1)**2+(l-1)+(1-l)+1
       ! now the value of the Pfunction
       Pfunc=qmat(k+2,3)*rlmat(combiLMu,combi2)-&
            qmat(k+2,1)*RLmat(combiLMu,combi3)
    ENDIF

    IF (m2.EQ.-l) THEN
       ! calculating the combined indices for RLmat        
       combi2=(l-1)**2+(l-1)+(1-l)+1
       combi3=(l-1)**2+(l-1)+(l-1)+1
       ! now the value of the Pfunction
       Pfunc=qmat(k+2,3)*rlmat(combiLMu,combi2)+&
            qmat(k+2,1)*RLmat(combiLMu,combi3)
    ENDIF
  END FUNCTION Pfunc
  ! ==================================================================
  REAL(real_8) FUNCTION bigV(l,m1,m2,qmat,rlmat,Lmax)
    ! ==--------------------------------------------------------------==
    ! bigV function
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2,Lmax
    REAL(real_8) :: qmat(3,3)
    REAL(real_8) :: rlmat((Lmax+1)**2,(lmax+1)**2)

    IF (m1.EQ.0) THEN
       bigV=Pfunc(1,l,1,m2,qmat,rlmat,Lmax)+&
            Pfunc(-1,L,-1,M2,qmat,RLmat,Lmax)
    ENDIF
    IF (m1.GT.0) THEN
       bigV=Pfunc(1,l,m1-1,m2,qmat,rlmat,Lmax)*&
            SQRT(1._real_8+REAL(IDelta(M1-1),kind=real_8))-&
            Pfunc(-1,L,1-M1,M2,qmat,RLmat,Lmax)*&
            (1._real_8-REAL(IDelta(M1-1),kind=real_8))
    ENDIF
    IF (m1.LT.0) THEN
       bigV=Pfunc(1,l,m1+1,m2,qmat,rlmat,Lmax)*&
            (1._real_8-REAL(IDelta(M1+1),kind=real_8)) +&
            Pfunc(-1,L,-1-M1,M2,qmat,RLmat,Lmax)*&
            SQRT(1._real_8+REAL(IDelta(M1+1),kind=real_8))
    ENDIF
  END FUNCTION bigV
  ! ==================================================================
  REAL(real_8) FUNCTION bigW(l,m1,m2,qmat,rlmat,Lmax)
    ! ==--------------------------------------------------------------==
    ! bigW function
    ! ==--------------------------------------------------------------==
    INTEGER :: l,m1,m2,Lmax
    REAL(real_8) :: qmat(3,3)
    REAL(real_8) :: rlmat((Lmax+1)**2,(lmax+1)**2)

    IF (m1.EQ.0) THEN
       bigW=0
    ENDIF
    IF (m1.GT.0) THEN
       bigW=Pfunc(1,l,m1+1,m2,qmat,rlmat,Lmax)+&
            Pfunc(-1,L,-1-M1,M2,qmat,RLmat,Lmax)
    ENDIF
    IF (m1.LT.0) THEN
       bigW=Pfunc(1,l,m1-1,m2,qmat,rlmat,Lmax)-&
            Pfunc(-1,L,1-M1,M2,qmat,RLmat,Lmax)
    ENDIF
  END FUNCTION bigW
  ! ==================================================================
  SUBROUTINE calc_wannier_energy_manno(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! This subroutine calculates the energy of a system
    ! based on the value of the coefficients (C0) of  
    ! the G vectors. (Thanks to Aldo)                
    ! ==--------------------------------------------------------------==
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)


    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate,1)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_wannier_energy_manno'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), psi(:,:)
    INTEGER                                  :: ierr, il_psi_1d, il_psi_2d, &
                                                il_rhoe_1d, il_rhoe_2d, &
                                                lforces
    REAL(real_8), ALLOCATABLE                :: dens(:,:), scr(:)

! Assumed defined in includes files or above:
! ip_fion, fion, tau0, NNR1, LFORCES
! NECESSARY DEFs
! end of NECESSARY DEFs

    CALL stopgm('ROTATE_..._P','not parallel.',& 
         __LINE__,__FILE__)

    ! Allocating the array Cs (scratch array?)
    ALLOCATE(cs(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! Allocating scratch space !
    CALL give_scr_forces(lforces,tag,nstate,.FALSE.,.FALSE.)
    IF (paral%io_parent)&
         WRITE(6,*) 'Scratch dimensions = ',lforces
    ALLOCATE(scr(lforces),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! Calculation of the energy of the wannier function (C) Aldo
    IF (paral%io_parent)&
         WRITE(6,*) 'Entering Energy calculation'

    CALL fnlalloc(nstate,.FALSE.,.FALSE.)

    ! Calculating the size of the wf and density arrays
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    IF (paral%io_parent)&
         WRITE(6,*) 'Rho size (',il_rhoe_1d,il_rhoe_2d,&
         ') psi size (',il_psi_1d,il_psi_2d,')'

    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dens(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(fion)!,SIZE(fion))

    ! Calculate energy!
    ! Prefactor
    CALL phfac(tau0)
    ! The forces routine calculates the energy
    CALL forces(c0,cs,tau0,fion,dens,psi,&
         NSTATE,1,.FALSE.,.FALSE.)
    IF (paral%io_parent)&
         WRITE(6,*) 'WANNIER ENERGY'
    IF (paral%io_parent)&
         WRITE(6,*) 'ETOTAL = ',ener_com%etot,' A. U.'
    IF (paral%io_parent)&
         WRITE(6,*) 'EKIN = ',ener_com%ekin,' A. U.'
    IF (paral%io_parent)&
         WRITE(6,*) 'EHT = ',ener_com%eht,' A. U.'
    IF (paral%io_parent)&
         WRITE(6,*) 'E[rho] = ',ener_com%eht+ener_com%epseu+ener_com%exc,' A. U.'
    IF (paral%io_parent)&
         WRITE(6,*)
    ! Freeing the used memory
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dens,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE calc_wannier_energy_manno

END MODULE rotate_my_wannier_manno_p_utils
