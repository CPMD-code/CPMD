! ==================================================================
PROGRAM plotband
  ! ==--------------------------------------------------------------==
  USE inscan_utils,                    ONLY: inscan
  USE machine,                         ONLY: m_cputime,&
                                             m_datum,&
                                             m_getarg,&
                                             m_iargc,&
                                             m_walltime

  IMPLICIT NONE
  INTEGER, PARAMETER :: maxpanels=20 
  REAL(real_8), PARAMETER :: rydberg=13.606_real_8 
  CHARACTER (len=26) :: date
  CHARACTER (len=80) :: line,lineold
  INTEGER :: ipanel(maxpanels)
  REAL(real_8) :: rki(3,maxpanels),rkf(3,maxpanels),rk(3)

  REAL(real_8), ALLOCATABLE :: eband(:)
  INTEGER :: iunit,icarg,ierr,nkpts,npanel,nkpbd,i,&
       ibrav,nkpts0,nbands,iunitout,nkpoint,ip,iik,ikind,&
       iband,ii
  REAL(real_8) :: time1,wclk1,efermi,eshift,alat,wk
  ! ==--------------------------------------------------------------==
  time1 = m_cputime()
  wclk1 =m_walltime()
  CALL m_datum(date)
  IF (paral%io_parent)&
       WRITE(6,'(A,A)') ' PROGRAM PLOTBAND STARTED AT: ',date
  ! ==--------------------------------------------------------------==
  iunit = 9
  icarg=m_iargc()
  IF (icarg.LT.1) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' NO INPUT FILE NAME SPECIFIED '
     CALL stopgm('PLOTBAND',' ',& 
          __LINE__,__FILE__)
  ENDIF
  CALL m_getarg(1,line)
  IF (paral%io_parent)&
       OPEN(unit=iunit,file=line,status='OLD',err=200)
  ierr=inscan(iunit,'&SYSTEM')
  IF (ierr.NE.0) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' PLOTBAND| INPUT SECTION &SYSTEM NOT FOUND'
     CALL stopgm('CONTROL',' ',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  ! Read the input file
  line=''
  DO WHILE (INDEX(line,'&END').EQ.0)
     lineold=line
     IF (paral%io_parent)&
          READ(iunit,err=99,END=99,fmt='(A80)') line
     IF ( (INDEX(line,'KPOINT').NE.0).AND.&
          (INDEX(line,'BAND').NE.0) ) THEN
        nkpts=0
        npanel=0
        nkpbd=1
        DO WHILE (nkpbd.NE.0)
           npanel=npanel+1
           IF (paral%io_parent)&
                READ(iunit,err=99,END=99,fmt=*) nkpbd,&
                (rki(i,npanel),i=1,3), (rkf(i,npanel),i=1,3)
           IF (nkpbd.NE.0) THEN
              IF (nkpbd.LE.1) CALL stopgm('  SYSIN',&
                   'WRONG NUMBER OF K POINTS PER BAND',& 
                   __LINE__,__FILE__)
              ipanel(npanel)=nkpbd
              nkpts=nkpts+nkpbd
           ELSE
              npanel=npanel-1
           ENDIF
        ENDDO
        GOTO 100
     ENDIF
  ENDDO
  IF (paral%io_parent)&
       WRITE(6,*) 'NO KPOINT BAND INFORMATION IN INPUT FILE'
  CALL stopgm('PLOTBAND',' ',& 
       __LINE__,__FILE__)    
  ! ==--------------------------------------------------------------==
100 CONTINUE
  IF (paral%io_parent)&
       WRITE(6,'(A,$)') 'FERMI ENERGY: '
  IF (paral%io_parent)&
       READ(*,*) efermi
  IF (paral%io_parent)&
       WRITE(6,'(A,$)') 'SHIFT IN ENERGY: '
  IF (paral%io_parent)&
       READ(*,*) eshift
  ! We have information about panel, we read ENERGYBANDS.
  IF (paral%io_parent)&
       CLOSE(iunit)
  IF (paral%io_parent)&
       OPEN(unit=iunit,file='ENERGYBANDS',status='OLD',err=999)
  ! Read ALAT
  IF (paral%io_parent)&
       READ(iunit,err=999,END=999,fmt='(1X,A8,I12,F20.15)')&
       line(1:8), ibrav,alat
  IF (paral%io_parent)&
       READ(iunit,err=999,END=999,fmt='(1X,2I12)') nkpts0,nbands
  IF (nkpts0.NE.nkpts) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' PLOTBAND| WRONG NUMBER IN KPOINTS IN ENERGYBANDS'
     IF (paral%io_parent)&
          WRITE(6,*) '           NKPTS=',nkpts0,' INSTEAD OF',nkpts
     CALL stopgm('PLOTBAND',' ',& 
          __LINE__,__FILE__)
  ENDIF
  ALLOCATE(eband(nbands),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  iunitout=10
  IF (paral%io_parent)&
       OPEN(unit=iunitout,file='bnds.ext',status='UNKNOWN')
  IF (paral%io_parent)&
       WRITE(iunitout,'(I5,2F10.5)') nbands, efermi/rydberg, alat
  nkpoint=0
  DO ip=1,npanel
     IF (paral%io_parent)&
          WRITE(iunitout,'(I5)') ipanel(ip)
     DO iik=1,ipanel(ip)
        nkpoint=nkpoint+1
        IF (paral%io_parent)&
             READ(iunit,err=999,END=999,fmt='(I5,4(F12.6))')&
             ikind,(rk(i),i=1,3),wk
        IF (nkpoint.NE.ikind) GOTO 999
        DO iband=1,nbands
           IF (paral%io_parent)&
                READ(iunit,*) ii,eband(iband)
        ENDDO
        IF (paral%io_parent)&
             WRITE(iunitout,'(3F10.5/,(10F8.4))')&
             (rk(i),i=1,3),((eband(i)-eshift)/rydberg,i=1,nbands)
     ENDDO
  ENDDO
  IF (paral%io_parent)&
       WRITE(iunitout,'(I5,/,I5)') 0,0
  DEALLOCATE(eband,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  GOTO 99999
  ! ==--------------------------------------------------------------==
99 CONTINUE
  IF (paral%io_parent)&
       WRITE(6,*) ' PLOTBAND| ERROR IN READING INPUT FILE'
  IF (paral%io_parent)&
       WRITE(6,*) ' THE LAST TWO LINES READ WERE '
  IF (paral%io_parent)&
       WRITE(6,*) lineold
  IF (paral%io_parent)&
       WRITE(6,*) line
  CALL stopgm('PLOTBAND',' ',& 
       __LINE__,__FILE__)
200 CONTINUE
  IF (paral%io_parent)&
       WRITE(6,*) ' PLOTBAND| INPUT FILE NOT FOUND'
  CALL stopgm('PLOTBAND',' ',& 
       __LINE__,__FILE__)
999 CONTINUE
  IF (paral%io_parent)&
       WRITE(6,*) ' PLOTBAND| ERROR IN READING FILE ENERGYBANDS'
  CALL stopgm('PLOTBAND',' ',& 
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  CALL m_datum(date)
  IF (paral%io_parent)&
       WRITE(6,'(A,A)') ' PROGRAM PLOTBAND ENDED AT: ',date
99999 CONTINUE
END PROGRAM plotband
! ==================================================================
