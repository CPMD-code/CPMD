MODULE geofile_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fbohr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_nofp,&
                                             fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: tshl
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             cpat,&
                                             cpsp,&
                                             mm_go_mm,&
                                             mm_revert,&
                                             nat_grm
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: ipcurr
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: cntl,&
                                             maxsys
  use bicanonicalCpmd, only: biCanonicalEnsembleDo,&
    bicanonicalCpmdConfig, getNameGeometryTape,&
    getNameGeometryXyzTape, getNameGeometryScaleTape
  USE velupi_utils,                    ONLY: c_to_s

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: geofile

CONTAINS

  ! ==================================================================
  SUBROUTINE geofile(tau0,velp,my_tag)
    ! ==--------------------------------------------------------------==
    ! == READ AND WRITE THE GEOMETRY FILE                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:)
    CHARACTER(len=*)                         :: my_tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'geofile'
    INTEGER, PARAMETER                       :: iunit = 12 

    CHARACTER(len=100)                       :: filen, filen2, filen3
    CHARACTER(len=12)                        :: cflbod
    CHARACTER(len=5)                         :: tag
    INTEGER                                  :: i, i1, i2, ia, ierr, ip, is, k
    INTEGER, ALLOCATABLE                     :: gr_iat(:)
    LOGICAL                                  :: ferror, status
    REAL(real_8), ALLOCATABLE                :: gr_tau(:,:), gr_vel(:,:), &
                                                taus(:,:,:)

! ==--------------------------------------------------------------==

    tag=my_tag
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(gr_tau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_vel(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tpath.OR.tmw) THEN
       IF (cntl%tpath)ip=ipcurr
       IF (tmw)ip=mwi%walker_id
       IF (ip.EQ.0) THEN
          filen='GEOMETRY'
          filen2='GEOMETRY.xyz'
          filen3='GEOMETRY.scale'
       ELSE
          cflbod='GEOMETRY_'
          CALL mw_filename(cflbod,filen,ip)
          CALL xstring(filen,i1,i2)
          filen2=filen(i1:i2)//'.xyz'
          filen3=filen(i1:i2)//'.scale'
       ENDIF
     else if (biCanonicalEnsembleDo) then
       filen = getNameGeometryTape(bicanonicalCpmdConfig)
       filen2 = getNameGeometryXyzTape(bicanonicalCpmdConfig)
       filen3 = getNameGeometryScaleTape(bicanonicalCpmdConfig)
    ELSE
       filen='GEOMETRY'
       filen2='GEOMETRY.xyz'
       filen3='GEOMETRY.scale'
    ENDIF
    IF (tshl%txfmqc.AND.tmw) THEN
       cflbod='GEOMETRY_'
       CALL mw_filename(cflbod,filen,ip)
       CALL xstring(filen,i1,i2)
       filen2=filen(i1:i2)//'.xyz'
    ENDIF
    IF (INDEX(tag,'READ').NE.0) THEN
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_old+fo_nofp,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' GEOFILE! FILE GEOMETRY DOES NOT EXIST!!'
          CALL stopgm('GEOFILE','DO NOT READ OLD COORDINATES',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            REWIND(iunit)
       IF (cntl%tqmmm)THEN
          DO i=1,ions1%nat
             IF (paral%io_parent)&
                  READ(iunit,err=998,fmt=*)&
                  (gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
             DO k=1,3
                tau0(k,cpat(i),cpsp(i))=gr_tau(k,i)+clsaabox%mm_c_trans(k)
                velp(k,cpat(i),cpsp(i))=gr_vel(k,i)
             ENDDO
          ENDDO
       ELSE
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                IF (paral%io_parent)&
                     READ(iunit,err=998,fmt=*)&
                     (tau0(k,ia,is),k=1,3),(velp(k,ia,is),k=1,3)
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,A,A,/)') ' ATOMIC COORDINATES AND VELOCITIES READ ',&
            'FROM FILE GEOMETRY'
       IF (cntl%tprcp) THEN
          DO i=1,3
             IF (paral%io_parent)&
                  READ(iunit,err=999,fmt=*) metr_com%ht(i,1),metr_com%ht(i,2),metr_com%ht(i,3),&
                  metr_com%htvel(i,1),metr_com%htvel(i,2),metr_com%htvel(i,3)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(A,A,/)') ' CELL PARAMETERS AND VELOCITIES READ ',&
               'FROM FILE GEOMETRY'
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(iunit)
    ELSEIF (INDEX(tag,'WRITE').NE.0) THEN
       i=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             i=i+1
             gr_iat(nat_grm(i))=ions0%iatyp(is)
             DO k=1,3
                gr_tau(k,nat_grm(i))=tau0(k,ia,is)+clsaabox%mm_c_trans(k)
                gr_vel(k,nat_grm(i))=velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_def+fo_nofp,ferror)
       IF (paral%io_parent)&
            REWIND(iunit)
       DO i=1,ions1%nat
          IF (paral%io_parent)&
               WRITE(iunit,'(3F20.12,8X,3F20.12)')&
               (gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
       ENDDO
       IF (cntl%tprcp) THEN
          DO i=1,3
             IF (paral%io_parent)&
                  WRITE(iunit,*) metr_com%ht(i,1),metr_com%ht(i,2),metr_com%ht(i,3),&
                  metr_com%htvel(i,1),metr_com%htvel(i,2),metr_com%htvel(i,3)
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(iunit)

       ! WRITE ALSO A GEOMETRY.xyz file:
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen2,fo_def+fo_nofp,ferror)
       IF (paral%io_parent)&
            REWIND(iunit)
       IF (paral%io_parent)&
            WRITE(iunit,'(I8)') ions1%nat
       IF (paral%io_parent)&
            WRITE(iunit,'(A)') 'GEOMETRY FILE / created by CPMD'
       DO i=1,ions1%nat
          IF (paral%io_parent)&
               WRITE(iunit,'(A3,3F20.12,8X,3F20.12)') elem%el(gr_iat(i)),&
               (gr_tau(k,i)/fbohr,k=1,3),(gr_vel(k,i)/fbohr,k=1,3)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(iunit)

       IF (cntl%tscale) THEN
          ! WRITE ALSO A GEOMETRY.scale file:
          IF (paral%io_parent)&
               CALL fileopen(iunit,filen3,fo_def+fo_nofp,ferror)
          IF (paral%io_parent)&
               REWIND(iunit)
          IF (paral%io_parent)&
               WRITE(iunit,'(A)') "CELL MATRIX (cntl%bohr)"
          DO i=1,3
             IF (paral%io_parent)&
                  WRITE(iunit,'(3F20.12)') metr_com%ht(i,1),metr_com%ht(i,2),metr_com%ht(i,3)
          ENDDO
          IF (paral%io_parent)&
               WRITE(iunit,'(/,A)') "CELL MATRIX (ANGSTROM)"
          DO i=1,3
             IF (paral%io_parent)&
                  WRITE(iunit,'(3F20.12)')&
                  metr_com%ht(i,1)/fbohr,metr_com%ht(i,2)/fbohr,metr_com%ht(i,3)/fbohr
          ENDDO
          ALLOCATE(taus(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL c_to_s(tau0,taus)
          IF (paral%io_parent)&
               WRITE(iunit,'(/,A)') "SCALED ATOMIC COORDINATES"
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                IF (paral%io_parent)&
                     WRITE(iunit,'(3F20.12,8X,A3,I7)')&
                     (taus(k,ia,is),k=1,3),elem%el(ions0%iatyp(is)),ia
             ENDDO
          ENDDO
          DEALLOCATE(taus,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          IF (paral%io_parent)&
               CALL fileclose(iunit)
       ENDIF

    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' GEOFILE| UNKNOWN TAG: ',tag
    ENDIF
    CALL mm_dim(mm_revert,status)
    DEALLOCATE(gr_tau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_vel,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_iat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
    ! Error reading TAU0
998 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' GEOFILE ! ERROR WHILE READING COORDINATES'
    CALL stopgm('GEOFILE','READ ERROR',& 
         __LINE__,__FILE__)
999 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' GEOFILE ! ERROR WHILE READING VELOCITIES'
    CALL stopgm('GEOFILE','READ ERROR',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE geofile
  ! ==================================================================

END MODULE geofile_utils
