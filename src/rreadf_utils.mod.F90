MODULE rreadf_utils
  USE cnst,                            ONLY: fbohr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: grandparent,&
                                             pimd3,&
                                             repfname,&
                                             supergroup,&
                                             supersource,&
                                             trep
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE velupi_utils,                    ONLY: s_to_c

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rreadf

CONTAINS

  ! ==================================================================
  SUBROUTINE rreadf
    ! ==--------------------------------------------------------------==
    ! ==   Read replica coordinates from file                         ==
    ! ==--------------------------------------------------------------==
    ! DM1
    ! DM2
    CHARACTER(*), PARAMETER                  :: procedureN = 'rreadf'
    INTEGER, PARAMETER                       :: iunit = 34 

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ia, iat, ib, ie, ierr, ip, &
                                                is, k
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: dis, dismax, dismin, &
                                                rmean(3), rnp, var
    REAL(real_8), ALLOCATABLE                :: taus(:,:,:)

! DM1
! DM2
! ==--------------------------------------------------------------==

    IF (grandparent) THEN
       CALL xstring(repfname,ib,ie)
       rnp=REAL(pimd3%np_total,kind=real_8)
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(iunit,repfname(ib:ie),fo_old,ferror)
       IF (ferror) CALL stopgm('RREADF','CANNOT OPEN REPLICA FILE',& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            REWIND(iunit)
       ! DM1
       ! NOTE EMPTY LINE IN FRONT OF EVERY REPLICA BLOCK
       ! DM2
       ! We use these empty lines to specify some keywords:
       ! ANGSTROM cntl%bohr SCALE
       DO ip=1,pimd3%np_total
          IF (paral%io_parent)&
               READ(iunit,'(A)') line
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                IF (paral%io_parent)&
                     READ(iunit,*) (trep(i,ia,is,ip),i=1,3)
             ENDDO
          ENDDO
          IF (INDEX(line,'SCALE').NE.0) THEN
             ! Scale options: atomic coordinates in direct lattice basis
             ALLOCATE(taus(3,maxsys%nax,ions1%nsp),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ip),1,taus,1)
             CALL s_to_c(taus,trep(:,:,:,ip))
             DEALLOCATE(taus,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ELSEIF (INDEX(line,'ANGST').NE.0) THEN
             ! Conversion into cntl%bohr
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   DO i=1,3
                      trep(i,ia,is,ip)=trep(i,ia,is,ip)*fbohr
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(iunit)
       IF (paral%io_parent)&
            WRITE(6,'(A,A,A,A)') ' REPLICA GEOMETRIES READ FROM FILE:',&
            '"',repfname(ib:ie),'"'
       ! DM1
       ! Only for path integral calculations
       IF (cntl%tpimd) THEN
          ! ..check the radius of gyration and compare to free particle
          ! ..(Gillan in Computer Modelling of Fluids Polymers and Solids 1990)
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,63("*"))')
          IF (paral%io_parent)&
               WRITE(6,'(" *",61X,"*")')
          IF (paral%io_parent)&
               WRITE(6,'(" *",23X,A,23X,"*")') ' R OF GYRATION '
          IF (paral%io_parent)&
               WRITE(6,'(" *",61X,"*")')
          IF (paral%io_parent)&
               WRITE(6,'(1X,63("*"),/)')
          IF (paral%io_parent)&
               WRITE(6,'(A,A,A,A)') '   ATOM   R OF GYRATION',&
               '  MAXDIS',&
               '   MINDIS',&
               '   (ANGSTROM)'
          ! 
          iat=0
          DO is=1,ions1%nsp
             ! DM         VARFP=(315777.0_real_8/(4.0_real_8*PMAR(IS)*TEMPB))
             ! DM         VARFP=SQRT(VARFP)
             DO ia=1,ions0%na(is)
                DO k=1,3
                   rmean(k)=0.0_real_8
                   DO ip=1,pimd3%np_total
                      rmean(k)=rmean(k)+trep(k,ia,is,ip)
                   ENDDO
                   rmean(k)=rmean(k)/rnp
                ENDDO
                var=0.0_real_8
                DO ip=1,pimd3%np_total
                   DO k=1,3
                      var=var+(trep(k,ia,is,ip)-rmean(k))**2
                   ENDDO
                ENDDO
                var=SQRT(var/rnp)
                dismax=-10._real_8
                dismin=100000.0_real_8
                DO ip=1,pimd3%np_total
                   IF (ip.EQ.pimd3%np_total) THEN
                      dis=(trep(1,ia,is,ip)-trep(1,ia,is,1))**2+(trep(2,ia,&
                           is,ip)-trep(2,ia,is,1))**2+(trep(3,ia,is,ip)-trep(3,ia,is,1))&
                           **2
                   ELSE
                      dis=(trep(1,ia,is,ip)-trep(1,ia,is,ip+1))**2+(trep(2,&
                           ia,is,ip)-trep(2,ia,is,ip+1))**2+(trep(3,ia,is,ip)-trep(3,ia,&
                           is,ip+1))**2
                   ENDIF
                   dis=SQRT(dis)
                   IF (dis.GT.dismax) dismax=dis
                   IF (dis.LT.dismin) dismin=dis
                ENDDO
                ! 
                iat=iat+1
                IF (paral%io_parent)&
                     WRITE(6,'(I4,2X,A2,3X,F8.5,4X,F8.5,1X,F8.5)')
                ! DM           WRITE(6,'(I4,2X,A2,3X,F8.5,7X,F8.5,4X,F8.5,1X,F8.5)')IAT,EL(IATYP(IS)),VAR/FBOHR,DISMAX/FBOHR,DISMIN/FBOHR
                ! DM  *        VAR/FBOHR,VARFP0/FBOHR,DISMAX/FBOHR,DISMIN/FBOHR
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,*) '  '
       ENDIF                 ! IF(cntl%tpimd)
       ! DM2
    ENDIF                     ! IF(GRANDPARENT)
    CALL mp_bcast(trep,SIZE(trep),supersource,supergroup)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rreadf
  ! ==================================================================

END MODULE rreadf_utils
