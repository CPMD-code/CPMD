MODULE rotate_my_wannier_para_p_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cppt,                            ONLY: gk
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old,&
                                             fo_ufo
  USE fnlalloc_utils,                  ONLY: fnlalloc
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_driver,                   ONLY: forces
  USE forces_utils,                    ONLY: give_scr_forces
  USE geq0mod,                         ONLY: geq0
  USE gndstate_p_utils,                ONLY: gauss_gen_read
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE lowdin_utils,                    ONLY: lowdin
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: bdxvecx,&
                                             bdxvecy,&
                                             bdxvecz,&
                                             dmbi,&
                                             emono,&
                                             statemax,&
                                             statemin
  USE symm,                            ONLY: sxscale,&
                                             syscale,&
                                             szscale
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotatemywannier_para
  !public :: save_wann_restart_p
  !public :: read_wann_restart_p
  !public :: read_multi_wann_restart_p
  PUBLIC :: shift_my_wannier
  PUBLIC :: calc_wannier_energy

CONTAINS

  ! 1         2         3         4         5         6         7    
  ! 23456789012345678901234567890123456789012345678901234567890123456789012
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! Subroutine providing the entry + necessary initialisation calls for   C
  ! the interaction of localised orbitals (PT interaction)                C
  ! C
  ! ------------------------ PARALLEL VERSION -------------------------   C
  ! C
  ! WARNING: this version DOES NOT perform any wavefunction rotation      C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE rotatemywannier_para(nstate,c0,psi)
    ! 
    ! Including the "essential" variables...
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate), psi(:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'rotatemywannier_para'

    COMPLEX(real_8), ALLOCATABLE             :: BCrot(:,:), cscr(:,:)
    INTEGER                                  :: ierr, imol, istate, isub, &
                                                jstate
    REAL(real_8)                             :: bovl2(nstate,nstate), &
                                                cart_coord(3), DXVnorm, &
                                                Ebefore, scaled_coord(3), sum

! scratch variables
! 
! cntl%timing variable
!cmb-bugfix    EXTERNAL rsph,rsphgen  ! not existent here and not used !
! Initialising the cntl%timing for the routine...

    CALL tiset('  main_rot',isub)

    ! Allocating the memory for the arrays...
    IF (paral%io_parent)&
         WRITE(6,*) 'Allocating memory'
    ALLOCATE(BCrot(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cscr(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(STATEmin(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(STATEmax(dmbi%nmol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,*) 'Done with memory allocation'
    IF (paral%io_parent)&
         WRITE(6,*)
    ! Checking initialisations

    ! DXVEC?
    IF (dmbi%in_bdxvec.EQV..FALSE.) THEN
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
       DO imol=1,dmbi%nmol
          BDXVECx(imol)=0._real_8
          BDXVECy(imol)=0._real_8
          BDXVECz(imol)=0._real_8
       ENDDO
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
          CALL read_multi_wann_restart_p(nstate,c0)
       ELSE
          ! NO, one reference for all molecules
          ! GAUSSIAN representation?
          IF (dmbi%wann_gauss) THEN
             ! YES...
             CALL gauss_gen_read(nstate,c0,dmbi%nmol,STATEmin,STATEmax)
          ELSE
             ! No, usual C0 reading...
             CALL read_wann_restart_p(nstate,c0)
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
       CALL calc_wannier_energy(nstate,c0)

       ! GAUSSIAN representation?
       IF (dmbi%wann_gauss) THEN
          ! YES...
          CALL save_my_wannier_gaussian(nstate,c0,ener_com%etot,psi)
       ELSE
          ! No, usual G-space saving
          CALL save_wann_restart_p(nstate,c0,ener_com%etot)
       ENDIF
    ENDIF

    ! Storing the value of the total energy before rotation
    Ebefore=ener_com%etot

    ! There are no rotations in this implementation so we copy C0 to BCrot
    CALL dcopy(2*ncpw%ngw*nstate,c0,1,BCrot,1)

    ! Main loop over the molecules
    DO imol=1,dmbi%nmol
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) '#########################'
          IF (paral%io_parent)&
               WRITE(6,'(1x,a,i6,a)') '# DISPL MOLECULE ',imol,' #'
          IF (paral%io_parent)&
               WRITE(6,*) '#########################'
       ENDIF

       ! Displacement vector
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) '# Displacement parameters:'
          IF (paral%io_parent)&
               WRITE(6,fmt='(1X,A,3(F8.5,A))') '# DXVEC: (',BDXVECx(imol),&
               ',',BDXVECy(imol),',',BDXVECz(imol),' )'
       ENDIF

       ! Calculating norm of displacement vector (if ==0, no shift..)
       DXVnorm=SQRT(BDXVECx(imol)**2+BDXVECy(imol)**2&
            +BDXVECz(imol)**2)
       IF (paral%io_parent)&
            WRITE(6,*) '# DXVnorm = ',dxvnorm
       IF (paral%io_parent)&
            WRITE(6,*) '# '

       ! Now we can shift the Wannier functions by the DX vector
       ! if its norm is greater than zero
       IF (DXVnorm.GT.0._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) '# Shifting the Wannier function'
          CALL shift_my_wannier(STATEmin(imol),STATEmax(imol),&
               NSTATE,BCrot,BDXVECx(imol),BDXVECy(imol),BDXVECz(imol))
          IF (paral%io_parent)&
               WRITE(6,*) '# FINISHED'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) '# NO TRANSLATION for molecule',imol
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*) '#########################'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDDO

    ! NOW: all the wannier functions are moved according to the 
    ! data of the input file...

    ! if requested orthogonalise the ROTATED states+renormalising... (CSCR is a scratch)
    IF (dmbi%torthog_wannier.AND..NOT.dmbi%tsimple_model) THEN
       CALL lowdin(BCrot,cscr,nstate)
       ! bit of cleaning just to be nice...
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ENDIF

    ! Overlap matrix for the rotated states
    IF (paral%io_parent)&
         WRITE(6,*) 'Building overlap matrix for the WANNIER states'

    IF (((dmbi%torthog_wannier).AND.(paral%parent)) .AND.paral%io_parent)&
         WRITE(6,*) 'Wave function fragments are now orthogonalised'
    ! this is already parallel...
    CALL ovlap(nstate,bovl2,BCrot,bcrot)
    CALL mp_sum(bovl2,nstate*nstate,parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,*) 'Checking wave function overlap matrix...'
    DO istate=1,nstate
       sum=0._real_8
       DO jstate=1,nstate
          ! checking matrix symmetry...
          IF (bovl2(istate,jstate).NE.bovl2(jstate,istate)) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Error in matrix symm. for states :',&
                  istate,'and ',jstate
          ENDIF
          ! checking sum of off-diagonal elements
          sum=sum+bovl2(istate,jstate)
       ENDDO
       IF (ABS(bovl2(istate,istate)- 1._real_8).GT.1.e-10_real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Error of normalisation for state ',&
               istate
       ENDIF

       IF (ABS(sum-bovl2(istate,istate)).GT.1.e-10_real_8) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'Detected overlap for state ',istate
             IF (paral%io_parent)&
                  WRITE(6,*) 'sum of off-diag elements = ',sum
          ENDIF
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) 'Finished'
    IF (paral%io_parent)&
         WRITE(6,*)
    ! Copying the rotated Wannier functions to the original C0
    CALL dcopy(2*ncpw%ngw*nstate,BCrot,1,c0,1)

    ! Stop the clock!!
    CALL tihalt('  main_rot',isub)

    ! Freeing the memory used by the routine
    DEALLOCATE(BCrot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE rotatemywannier_para

  ! ==================================================================
  SUBROUTINE save_wann_restart_p(nstate,c0,etot)
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! == Saves the Wannier functions in a restart file named:         ==
    ! == RESTART.ref                                                  ==
    ! ==                                                              ==
    ! == NOTE: This works in parallel! (taken from RESTART_P)         ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: etot

    CHARACTER(len=128)                       :: filename
    INTEGER                                  :: controlbyte, fileunit, ia, &
                                                ie, ierror, irecord, isub
    LOGICAL                                  :: Error, ferror

    CALL tiset('SAVE_WAN_P',isub)

    ! setup filename
    fileunit = 2606
    filename = 'RESTART.ref'
    CALL xstring(filename,ia,ie)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(fileunit,filename,fo_def+fo_ufo,ferror)
       IF (ferror) GOTO 401
       IF (paral%io_parent)&
            REWIND(fileunit,err=401)
       Error=.FALSE.
       GOTO 402
401    Error=.TRUE.
402    CONTINUE

       IF (Error) THEN
          CALL stopgm('SAVE_WAN_P',&
               'Error OPENING RESTART FILE ('&
               //filename(ia:ie)//') for writing.',& 
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,*) 'Saving ',nstate,' states...'
       IF (paral%io_parent)&
            WRITE(6,*) 'Total energy for this system:'
       IF (paral%io_parent)&
            WRITE(6,*) 'Etot = ',etot,'A. U.'
       controlbyte = 730909
       IF (paral%io_parent)&
            WRITE (fileunit) controlbyte
       IF (paral%io_parent)&
            WRITE (fileunit) parm%a1,parm%a2,parm%a3
       IF (paral%io_parent)&
            WRITE (fileunit) nstate,ncpw%ngw
    ENDIF               ! parent
    ! ==--------------------------------------------------------------==
    CALL wr30wfn(fileunit,ierror,nstate,c0,tau0,'NIL',irecord)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (fileunit) etot
       IF (paral%io_parent)&
            WRITE (6,'(A)')&
            ' WANNIER wavefunction written to RESTART file '&
            //filename(ia:ie)//'.'
       IF (paral%io_parent)&
            CALL fileclose(fileunit)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    CALL tihalt('SAVE_WAN_P',isub)
    RETURN
  END SUBROUTINE save_wann_restart_p
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE read_wann_restart_p(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! == Reads the Wannier functions from the restart file named      ==
    ! == RESTART.ref                                                  ==
    ! ==                                                              ==
    ! == NOTE: This works in parallel! (taken from RESTART_P)         ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)

    CHARACTER(*), PARAMETER :: procedureN = 'read_wann_restart_p'

    CHARACTER(len=128)                       :: filename
    COMPLEX(real_8), ALLOCATABLE             :: Ct(:,:)
    INTEGER :: controlbyte, fileunit, ia, ie, ierr, info, irecord, istate, &
      isub, kstate, mol, ngw_, nstate_, states_to_read
    INTEGER(int_8)                           :: fpos
    LOGICAL                                  :: Error, ferror
    REAL(real_8)                             :: a1_(3), a2_(3), a3_(3), etot

    CALL tiset('READ_WAN_P',isub)

    ! setup filename
    fileunit = 2606
    filename = 'RESTART.ref'
    CALL xstring(filename,ia,ie)
    ! Calculating the number of states to read
    states_to_read=nstate/dmbi%nmol

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'WLOAD needs',states_to_read,&
            'states for each molecule'

       IF (paral%io_parent)&
            CALL fileopen(fileunit,filename,fo_old+fo_ufo,ferror)
       IF (ferror) GOTO 501

       IF (paral%io_parent)&
            REWIND(fileunit,err=501)
       Error=.FALSE.
       GOTO 502
501    Error=.TRUE.
502    CONTINUE

       IF (Error) THEN
          CALL stopgm('READ_WAN_P',&
               'Error OPENING RESTART FILE ('&
               //filename(ia:ie)//') for reading.',& 
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! GENERAL PART of restart file:
       IF (paral%io_parent)&
            READ(fileunit) controlbyte
       IF (paral%io_parent)&
            READ(fileunit) a1_,a2_,a3_
       IF (paral%io_parent)&
            READ(fileunit) nstate_,ngw_
       Error=.FALSE.
       ! TEST our controlbyte... is it a wannier restart
       IF (controlbyte .NE. 730909)&
            CALL stopgm('READ_WAN_P',&
            'Error READING RESTART FILE ('&
            //filename(ia:ie)//'): wrong parameters.',& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,*) 'The file contains ',nstate_,' states'
       IF (paral%io_parent)&
            WRITE(6,*) 'Each state has ',ngw_,' G vectors'
       ! TEST number of plane waves: has to be the same
       ! if (NGW.NE.ngw_) then
       ! error = .true.
       ! call stopgm('READ_WAN_P',
       !           'The number of G vectors found in the file '
       !              //'do not match NGW in the simulation! '//
       !              'At present this is not allowed...')
       ! endif
       ! TEST number of states: enough for each molecule
       IF (states_to_read.NE.nstate_) THEN
          CALL stopgm('READ_WAN_P',&
               'Need more states for each molecules...exiting',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF               ! parent stuff.
    CALL mp_bcast(nstate_,parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==         
    ! Everything is OK...allocating the Ct array (2*dims for CPLX*16)
    ALLOCATE(Ct(ncpw%ngw,nstate_),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! reading the wannier functions
    IF (paral%io_parent)&
         WRITE(6,*) 'Proceeding... loading the data for ',&
         dmbi%nmol,' molecule(s)'
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) READ(fileunit) irecord,fpos! dummy read.
    CALL rd30wfn(fileunit,Ct,nstate_,tau0,info,&
         .FALSE.,nstate_,1,'NIL',fpos)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            READ(fileunit) etot
       IF (paral%io_parent)&
            WRITE(6,*) 'Total energy of the saved state'
       IF (paral%io_parent)&
            WRITE(6,*) 'Etot = ',etot,'A. U.'
       CALL xstring(filename,ia,ie)
       IF (paral%io_parent)&
            WRITE (6,'(A)')&
            ' WANNIER wavefunction read from RESTART file '&
            //filename(ia:ie)//'.'
       IF (paral%io_parent)&
            CALL fileclose(fileunit)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*) 'Building',nstate,' Wannier states'
    ! Now filling up the C0 array with the required ammount of data
    ! initialising the global state counter
    kstate=1
    ! Loop on number of molecules
    DO mol=1,dmbi%nmol
       ! Storing the energy of the monomers on parent
       IF (paral%parent) emono(mol)=etot
       ! Storing starting state for each monomer
       STATEmin(mol)=kstate
       ! Loop on number of read states
       DO istate=1,nstate_
          CALL dcopy(2*ncpw%ngw,Ct(1,istate),1,c0(1,kstate),1)
          ! updating the global state counter
          kstate=kstate+1
       ENDDO
       STATEmax(mol)=kstate-1
    ENDDO

    IF (paral%parent) THEN
       DO mol=1,dmbi%nmol
          IF (paral%io_parent)&
               WRITE(6,*) 'molecule: ',mol,' STATEmin = ',&
               STATEmin(mol),' STATEmax = ',STATEmax(mol)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) 'Done! '
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('READ_WAN_P',isub)
    ! Freeing memory
    DEALLOCATE(Ct,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE read_wann_restart_p
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE read_multi_wann_restart_p(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! == Reads a Wannier function for each molecule from a series     ==
    ! == of restart file named:                                       ==
    ! == RESTART.ref1                                                 ==
    ! == RESTART.ref2                                                 ==
    ! == RESTART.ref3                                                 ==
    ! ==   ...                                                        ==
    ! ==                                                              ==
    ! == NOTE: This works in parallel! (taken from RESTART_P)         ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)

    CHARACTER(*), PARAMETER :: procedureN = 'read_multi_wann_restart_p'

    CHARACTER(len=128)                       :: filename
    CHARACTER(len=8)                         :: number
    COMPLEX(real_8), ALLOCATABLE             :: Ct(:,:)
    INTEGER                                  :: controlbyte, fileunit, ia, &
                                                ie, ierr, info, irecord, &
                                                istate, isub, kstate, mol, &
                                                ngw_, nstate_
    INTEGER(int_8)                           :: fpos
    LOGICAL                                  :: Error
    REAL(real_8)                             :: a1_(3), a2_(3), a3_(3), etot

    CALL tiset('READMULW_P',isub)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*) '|   PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL  |'
       IF (paral%io_parent)&
            WRITE(6,*) '| Multi wannier reference file reading subroutine |'
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*) '(This  routine will read ',dmbi%nmol,' reference files)'
    ENDIF

    ! initialising the global state counter
    kstate=1
    ! main loop on the molecules...
    DO mol=1,dmbi%nmol

       ! setup filename
       fileunit = 2606

       ! First creating the file name for the molecule mol
       IF (paral%io_parent)&
            WRITE(number,'(I8)')mol
       ! removing the spaces before and after the number (CPMD routine)
       CALL xstring(number,ia,ie)

       filename = 'RESTART.ref'//number(ia:ie)
       CALL xstring(filename,ia,ie)

       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,*) '======== Reading data for molecule ',mol,&
               ' ========'
          IF (paral%io_parent)&
               WRITE(6,*) 'Opening file '//filename(ia:ie)
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          Error=.FALSE.
          IF (paral%io_parent)&
               CALL fileopen(fileunit,filename,fo_old+fo_ufo,Error)
          IF (Error) CALL stopgm('READMULW_P',&
               'RESTART FILE NOT FOUND: '//filename(ia:ie)//'.',& 
               __LINE__,__FILE__)
          IF (paral%io_parent)&
               REWIND(fileunit,err=501)
          Error=.FALSE.
          GOTO 502
501       Error=.TRUE.
502       CONTINUE
          IF (Error) THEN
             CALL stopgm('READMULW_P',&
                  'Error OPENING RESTART FILE ('&
                  //filename(ia:ie)//') for reading.',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ==--------------------------------------------------------------==
          ! GENERAL PART of restart file:
          IF (paral%io_parent)&
               READ(fileunit) controlbyte
          IF (paral%io_parent)&
               READ(fileunit) a1_,a2_,a3_
          IF (paral%io_parent)&
               READ(fileunit) nstate_,ngw_
          Error=.FALSE.
          ! TEST our controlbyte... is it a wannier restart
          IF (controlbyte .NE. 730909) Error = .TRUE.
          IF (paral%io_parent)&
               WRITE(6,*) 'The file contains ',nstate_,' states'
          IF (paral%io_parent)&
               WRITE(6,*) 'Each state has ',ngw_,' G vectors'
          ! TEST number of plane waves: has to be the same
          ! if (NGW.NE.ngw_) then
          ! error = .true.
          ! call stopgm('READMULW_P',
          !           'The number of G vectors found in the file '
          !              //'do not match NGW in the simulation! '//
          !              'At present this is not allowed...')
          ! endif
          IF (Error) THEN
             CALL stopgm('READMULW_P',&
                  'Error READING RESTART FILE ('&
                  //filename(ia:ie)//'): wrong parameters.',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF             ! parent stuff.
       CALL mp_bcast(nstate_,parai%source,parai%allgrp)
       ! ==--------------------------------------------------------------==         
       ! Everything is OK...allocating the Ct array (2*dims for CPLX*16)
       ALLOCATE(Ct(ncpw%ngw,nstate_),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Clear memory just in case...
       CALL zeroing(Ct)!,SIZE(Ct))
       ! ==--------------------------------------------------------------==
       ! reading the wannier functions
       IF (paral%io_parent) READ(fileunit) irecord,fpos! dummy read.
       CALL rd30wfn(fileunit,Ct,nstate_,tau0,info,&
            .FALSE.,nstate_,1,'NIL',fpos)
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               READ(fileunit) etot
          ! Storing monomer energies on the parent
          emono(mol)=etot
          IF (paral%io_parent)&
               WRITE(6,*) 'Total energy of the saved state'
          IF (paral%io_parent)&
               WRITE(6,*) 'Etot = ',etot,'A. U.'
          CALL xstring(filename,ia,ie)
          IF (paral%io_parent)&
               WRITE (6,'(A)')&
               ' WANNIER wavefunction read from RESTART file '&
               //filename(ia:ie)//'.'
          IF (paral%io_parent)&
               CALL fileclose(fileunit)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Now filling up the C0 array with the data from this file
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Building',nstate_,' Wannier states for ',&
               'molecule',mol
          IF (paral%io_parent)&
               WRITE(6,*) 'Assigning states',kstate,'to',&
               kstate+nstate_-1
       ENDIF
       ! Assigning first and last state number for this molecule      
       STATEmin(mol)=kstate
       STATEmax(mol)=STATEmin(mol)+nstate_-1

       ! copy each read state to C0 state by state
       DO istate=1,nstate_
          CALL dcopy(2*ncpw%ngw,Ct(1,istate),1,c0(1,kstate),1)
          ! updating the global state counter
          kstate=kstate+1
       ENDDO
       ! Freeing memory
       DEALLOCATE(Ct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDDO

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '======== Summary of assigned states  ========'
       DO mol=1,dmbi%nmol
          IF (paral%io_parent)&
               WRITE(6,*) 'molecule: ',mol,' STATEmin = ',&
               STATEmin(mol),' STATEmax = ',STATEmax(mol)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*) '|   PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL  |'
       IF (paral%io_parent)&
            WRITE(6,*) '| Multi wannier subroutine FINISHED (no errors)   |'
       IF (paral%io_parent)&
            WRITE(6,*) '+-------------------------------------------------+'
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('READMULW_P',isub)
    RETURN
  END SUBROUTINE read_multi_wann_restart_p
  ! ==================================================================


  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! This subroutine shifts the wannier function by a vector     C
  ! defined by Vx,Vy,Vz.                                        C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE shift_my_wannier(st_state,fi_state,nstate,c0,Vx,Vy,Vz)

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

    IF (paral%io_parent)&
         WRITE(6,*) 'DOING STATE ',st_state,' TO ',fi_state

    ! Main loop on the states
    DO istate=st_state,fi_state
       ! Loop on the G vectors
       DO Gindx=1,ncpw%ngw
          ! Calculating the product G.r
          GbyR=gk(1,Gindx)*Vx+gk(2,gindx)*Vy+gk(3,gindx)*Vz

          ! Adjusting the units
          GbyR=-gbyr*parm%tpiba

          ! Calculating E^i(G.r) using cos(x)+i*sin(x)
          EXPiGR=CMPLX(COS(GbyR),SIN(gbyr),kind=real_8)

          ! multiplying the Wannier coeffs
          c0(Gindx,istate)=c0(gindx,istate)*EXPiGR
       ENDDO
    ENDDO
  END SUBROUTINE shift_my_wannier

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! This subroutine calculates the energy of a system   C
  ! based on the value of the coefficients (C0) of      C
  ! the G vectors. (Thanks to Aldo)                     C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  SUBROUTINE calc_wannier_energy(nstate,c0)
    ! so from now on NGW is the number of plane waves, GK(3,NGW) is the array
    ! of G vectors, HG(NGW) is the array containing the square of the norm of
    ! each G vector. (according to Daniel...)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate,1)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_wannier_energy'

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
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    IF (paral%io_parent)&
         WRITE(6,*) 'Rho size (',il_rhoe_1d,il_rhoe_2d,&
         ') psi size (',il_psi_1d, il_psi_2d,')'

    ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
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

    IF (paral%parent) THEN
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
    ENDIF

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
  END SUBROUTINE calc_wannier_energy


END MODULE rotate_my_wannier_para_p_utils
