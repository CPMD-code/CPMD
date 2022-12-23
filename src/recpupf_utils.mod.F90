MODULE recpupf_utils
  USE adat,                            ONLY: atwt,&
                                             defrag,&
                                             elem
  USE atom,                            ONLY: atom_common,&
                                             ecpfiles,&
                                             gnl,&
                                             rv,&
                                             rw,&
                                             vr
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corecg,&
                                             corei,&
                                             corel,&
                                             corer,&
                                             rcgrid
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE recpnew_utils,                   ONLY: ckgrid,&
                                             get_pplib,&
                                             tgrid
  USE rmas,                            ONLY: rmass
  USE sgpp,                            ONLY: mpro,&
                                             sgpp1,&
                                             sgpp2
  USE system,                          ONLY: lmaxx,&
                                             maxsys
  USE vdbt,                            ONLY: itmax,&
                                             vdbti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: recpupf

CONTAINS

  ! ==================================================================
  SUBROUTINE recpupf(isp,ecpnam)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Pseudopotential Input (UPF Format)                    ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! == ISP     Species index                                        ==
    ! == ECPNAM  Filename of pseudopotential file                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    CHARACTER(len=*)                         :: ecpnam

    CHARACTER(*), PARAMETER                  :: procedureN = 'recpupf'
    INTEGER, PARAMETER                       :: iunit = 21 

    CHARACTER(len=1)                         :: endc
    CHARACTER(len=120)                       :: ecplib
    CHARACTER(len=200)                       :: fnames
    CHARACTER(len=5000)                      :: longline
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ia, ibc, ie, ierr, iout, &
                                                ir, it, iv, j, jv, k, l, &
                                                lenecp, lp, m, meshv, nproje, &
                                                nwfn
    LOGICAL                                  :: erread, exists, upf2
    REAL(real_8)                             :: dij
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: c, dij_tmp, fc, temp

! ==--------------------------------------------------------------==

    sgpp1%tsgp(isp)=.FALSE.
    pslo_com%tnum(isp)=.FALSE.
    pslo_com%tlog(isp)=.FALSE.
    CALL get_pplib(ecplib,lenecp)
    CALL xstring(ecpnam,ia,ie)
    fnames=ecplib(1:lenecp)//ecpnam(ia:ie)
    INQUIRE(file=fnames,exist=exists)
    IF (.NOT.exists) THEN
       fnames=ecpnam(ia:ie)
       INQUIRE(file=fnames,exist=exists)
       IF (.NOT.exists) THEN
          WRITE(6,*) ' RECPUPF| ECPFILE NOT FOUND ',fnameS
          CALL stopgm('RECPUPF',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ecpfiles(isp)=fnames
    OPEN(unit=iunit,file=fnames,status='OLD')
    REWIND(iunit)
    ! ==--------------------------------------------------------------==
    ! ..Determine which UPF version we are using
    upf2=.TRUE.
    ierr=inscan(iunit,'<UPF ')
    IF (ierr.NE.0) THEN
       upf2=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ..The info section: no difference between V1 and V2
    ierr=inscan(iunit,'<PP_INFO>')
    WRITE(vdbti(1,isp),'(66("+"))')
    vdbti(2,isp) = ecpnam(ia:ie)
    WRITE(vdbti(3,isp),'(66("+"))')
    it=3
    DO i=1,60
       READ(iunit,END=20,err=20,fmt='(A)') line
       IF (INDEX(line,'</PP_INFO>').NE.0) GOTO 10
       it=it+1
       vdbti(it,isp)=line(1:66)
    ENDDO
10  CONTINUE
    it=it+1
    WRITE(vdbti(it,isp),'(66("+"))')
    itmax(isp)=it
    ! ==--------------------------------------------------------------==
    ! ..General info about PP and atom
    IF (upf2) THEN
       CALL get_section_info(iunit,'<PP_HEADER',longline,&
            LEN(longline))
       ! Element
       CALL get_key(longline,"element",line)
       CALL xstring(line,ia,ie)
       CALL find_element(line(ia:ie),ions0%iatyp(isp))
       ! Ultrasoft|Normconserving
       CALL get_key(longline,"is_ultrasoft",line)
       IF (INDEX(line,"T").NE.0) pslo_com%tvan(isp)=.TRUE.
       IF (pslo_com%tvan(isp)) THEN
          CALL stopgm('RECPUPF','ULTRASOFT UPF FILES NOT SUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF
       ! NLCC
       CALL get_key(longline,"core_correction",line)
       CALL xstring(line,ia,ie)
       IF (line(ia:ia).EQ."T") corel%tnlcc(isp)=.TRUE.
       ! Valence charge
       CALL get_key(longline,"z_valence",line)
       CALL readsr(line,1,iout,ions0%zv(isp),erread)
       ! LMAX
       CALL get_key(longline,"l_max",line)
       CALL readsi(line,1,iout,dpot_mod%lmax(isp),erread)
       dpot_mod%lmax(isp)=dpot_mod%lmax(isp)+1
       dpot_mod%lloc(isp)=dpot_mod%lmax(isp)
       dpot_mod%lskip(isp)=dpot_mod%lmax(isp)+5
       sgpp1%tsgp(isp)=.TRUE.
       pslo_com%tnum(isp)=.TRUE.
       dpot_mod%tkb(isp)=.FALSE.
       ions0%igau(isp)=0
       ! Mesh points
       CALL get_key(longline,"mesh_size",line)
       CALL readsi(line,1,iout,meshv,erread)
       IF (meshv.GT.maxsys%mmaxx) THEN
          CALL xstring(ecpnam,ia,ie)
          WRITE(6,*) ' RECPUPF! MESH FOR ',ECPNAM(IA:IE),' IS ',&
               meshV
          WRITE(6,*) ' RECPUPF! MAX NUMBER OF SPLINE POINTS:',maxsys%mmaxx
          WRITE(6,*) ' RECPUPF! INCREASE SPLINE POINTS NUMBER'
          CALL stopgm('RECPUPF','MESH TOO BIG',& 
               __LINE__,__FILE__)
       ENDIF
       atom_common%meshvp(isp)=meshv
       ! Number of wavefunctions, projectors
       CALL get_key(longline,"number_of_wfc",line)
       CALL readsi(line,1,iout,nwfn,erread)
       CALL get_key(longline,"number_of_proj",line)
       CALL readsi(line,1,iout,nproje,erread)
    ELSE
       ierr=inscan(iunit,'<PP_HEADER>')
       IF (ierr.NE.0) THEN
          WRITE(6,*) ' RECPUPF| PP_HEADER SECTION NOT FOUND'
          CALL stopgm('RECPUPF',' ',& 
               __LINE__,__FILE__)
       ENDIF
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Version number (not used)
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Element
       CALL xstring(line,ia,ie)
       CALL find_element(line(ia:ie),ions0%iatyp(isp))
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Ultrasoft|Normconserving
       IF (INDEX(line,"US").NE.0) pslo_com%tvan(isp)=.TRUE.
       IF (pslo_com%tvan(isp)) THEN
          CALL stopgm('RECPUPF','ULTRASOFT UPF FILES NOT SUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! NLCC
       CALL xstring(line,ia,ie)
       IF (line(ia:ia).EQ."T") corel%tnlcc(isp)=.TRUE.
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! DFT FUNCTIONAL (not used)
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Valence charge
       CALL readsr(line,1,iout,ions0%zv(isp),erread)
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Total energy (not used)
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Suggested cutoffs (not used)
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! LMAX
       CALL readsi(line,1,iout,dpot_mod%lmax(isp),erread)
       dpot_mod%lmax(isp)=dpot_mod%lmax(isp)+1
       dpot_mod%lloc(isp)=dpot_mod%lmax(isp)
       dpot_mod%lskip(isp)=dpot_mod%lmax(isp)+5
       sgpp1%tsgp(isp)=.TRUE.
       pslo_com%tnum(isp)=.TRUE.
       dpot_mod%tkb(isp)=.FALSE.
       ions0%igau(isp)=0
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Mesh points
       CALL readsi(line,1,iout,meshv,erread)
       IF (meshv.GT.maxsys%mmaxx) THEN
          CALL xstring(ecpnam,ia,ie)
          WRITE(6,*) ' RECPUPF! MESH FOR ',ECPNAM(IA:IE),' IS ',&
               meshV
          WRITE(6,*) ' RECPUPF! MAX NUMBER OF SPLINE POINTS:',maxsys%mmaxx
          WRITE(6,*) ' RECPUPF! INCREASE SPLINE POINTS NUMBER'
          CALL stopgm('RECPUPF','MESH TOO BIG',& 
               __LINE__,__FILE__)
       ENDIF
       atom_common%meshvp(isp)=meshv
       READ(iunit,END=20,err=20,fmt='(A)') line
       ! Number of wavefunctions, projectors
       CALL readsi(line,1,iout,nwfn,erread)
       CALL readsi(line,iout+1,iout,nproje,erread)
       ! Further info (not used)
    ENDIF
    ! 
    rmass%pma0(isp)=atwt(ions0%iatyp(isp))
    raggio(isp)=defrag(ions0%iatyp(isp))
    ! ==--------------------------------------------------------------==
    ! ..Radial Mesh
    IF (upf2) THEN
       endc=" "
    ELSE
       endc=">"
    ENDIF
    ierr=inscan(iunit,'<PP_MESH'//endc)
    IF (ierr.NE.0) THEN
       WRITE(6,*) ' RECPUPF| PP_MESH SECTION NOT FOUND'
       CALL stopgm('RECPUPF',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ierr=inscan(iunit,'<PP_R'//endc)
    IF (ierr.NE.0) THEN
       WRITE(6,*) ' RECPUPF| PP_R SECTION NOT FOUND'
       CALL stopgm('RECPUPF',' ',& 
            __LINE__,__FILE__)
    ENDIF
    READ(iunit,END=20,err=20,fmt=*) (rv(ir,isp),ir=1,meshv)
    ierr=inscan(iunit,'<PP_RAB'//endc)
    IF (ierr.NE.0) THEN
       WRITE(6,*) ' RECPUPF| PP_R SECTION NOT FOUND'
       CALL stopgm('RECPUPF',' ',& 
            __LINE__,__FILE__)
    ENDIF
    READ(iunit,END=20,err=20,fmt=*) (rw(ir,isp),ir=1,meshv)
    ! ==--------------------------------------------------------------==
    ! ..Local Potential section
    IF (upf2) THEN
       endc=" "
    ELSE
       endc=">"
    ENDIF
    ierr=inscan(iunit,'<PP_LOCAL'//endc)
    IF (ierr.NE.0) THEN
       WRITE(6,*) ' RECPUPF| PP_LOCAL SECTION NOT FOUND'
       CALL stopgm('RECPUPF',' ',& 
            __LINE__,__FILE__)
    ENDIF
    READ(iunit,END=20,err=20,fmt=*) (vr(ir,isp,1),ir=1,meshv)
    CALL dscal(meshv,0.5_real_8,vr(1,isp,1),1)
    ! ==--------------------------------------------------------------==
    ! ..Nonlocal Potential section
    IF (nproje.GT.0) THEN
       ierr=inscan(iunit,'<PP_NONLOCAL')
       IF (ierr.NE.0) THEN
          WRITE(6,*) ' RECPUPF| PP_NONLOCAL SECTION NOT FOUND'
          CALL stopgm('RECPUPF',' ',& 
               __LINE__,__FILE__)
       ENDIF
       DO i=1,lmaxx
          sgpp2%npro(i,isp)=0
       ENDDO
       ibc=0
100    CONTINUE
       READ(iunit,END=20,err=20,fmt='(A)') line
       IF (INDEX(line,'</PP_NONLOCAL>').NE.0) THEN
          CALL stopgm('RECPUPF','ERROR READING NL PROJECTORS',& 
               __LINE__,__FILE__)
       ENDIF
       IF (INDEX(line,'<PP_BETA').NE.0) THEN
          ibc=ibc+1
          IF (upf2) THEN
             WRITE(endc,'(I0)')ibc
             CALL get_section_info(iunit,'<PP_BETA.'//endc,longline,&
                  LEN(longline))
             CALL get_key(longline,"index",line)
             CALL readsi(line,1,iout,i,erread)
             CALL get_key(longline,"angular_momentum",line)
             CALL readsi(line,1,iout,l,erread)
             CALL get_key(longline,"size",line)
             CALL readsi(line,1,iout,k,erread)
          ELSE
             READ(iunit,END=20,err=20,fmt=*) i,l
             READ(iunit,END=20,err=20,fmt=*) k
          ENDIF
          sgpp2%pplist(i,1,isp)=l
          sgpp2%npro(l+1,isp)=sgpp2%npro(l+1,isp)+1
          IF (sgpp2%npro(l+1,isp).GT.mpro) THEN
             CALL stopgm('RECPUPF','TOO MANY PROJECTORS',& 
                  __LINE__,__FILE__)
          ENDIF
          sgpp2%pplist(i,2,isp)=sgpp2%npro(l+1,isp)
          READ(iunit,END=20,err=20,fmt=*) (gnl(ir,isp,ibc),ir=1,k)
          DO ir=k+1,meshv
             gnl(ir,isp,ibc)=0._real_8
          ENDDO
       ELSE
          GOTO 100
       ENDIF
       IF (ibc.LT.nproje) GOTO 100
       ! ==--------------------------------------------------------------==
       ! ..set up 
       iv=0
       lp=0
       DO l=1,dpot_mod%lmax(isp)
          DO m=1,2*l-1
             lp=lp+1
             DO k=1,sgpp2%npro(l,isp)
                iv=iv+1
                nghtol(iv,isp)=l-1
                sgpp2%lpval(iv,isp)=lp
                sgpp2%lfval(iv,isp)=k
             ENDDO
          ENDDO
       ENDDO
       nlps_com%ngh(isp)=iv
       ! ==--------------------------------------------------------------==
       ! ..DIJ section
       ierr=inscan(iunit,'<PP_DIJ')
       IF (ierr.NE.0) THEN
          WRITE(6,*) ' RECPUPF| PP_DIJ SECTION NOT FOUND'
          CALL stopgm('RECPUPF',' ',& 
               __LINE__,__FILE__)
       ENDIF
       DO l=1,dpot_mod%lmax(isp)
          DO i=1,mpro
             DO j=1,mpro
                sgpp2%hlsg(i,j,l,isp)=0._real_8
             ENDDO
          ENDDO
       ENDDO
       IF (upf2) THEN
          ALLOCATE(dij_tmp(nproje*nproje),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          READ(iunit,END=20,err=20,fmt=*)&
               (dij_tmp(ir),ir=1,nproje*nproje)
          ir=0
          DO i=1,nproje
             DO j=1,nproje
                ir=ir+1
                IF (i/=j) THEN
                   IF (dij_tmp(ir)/=0._real_8) THEN
                      CALL stopgm('RECPUPF','NON-ZERO DIJ IS INVALID! ',& 
                           __LINE__,__FILE__)
                   ENDIF
                   CYCLE
                ENDIF
                IF (sgpp2%pplist(i,1,isp).NE.sgpp2%pplist(j,1,isp)) THEN
                   CALL stopgm('RECPUPF','INVALID DIJ',& 
                        __LINE__,__FILE__)
                ENDIF
                l=sgpp2%pplist(i,1,isp)+1
                iv=sgpp2%pplist(i,2,isp)
                jv=sgpp2%pplist(j,2,isp)
                sgpp2%hlsg(iv,jv,l,isp)=0.5_real_8*dij_tmp(ir)
                sgpp2%hlsg(jv,iv,l,isp)=0.5_real_8*dij_tmp(ir)
             ENDDO
          ENDDO
          DEALLOCATE(dij_tmp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ELSE
          READ(iunit,END=20,err=20,fmt=*) m
          DO k=1,m
             READ(iunit,END=20,err=20,fmt=*) i,j,dij
             IF (sgpp2%pplist(i,1,isp).NE.sgpp2%pplist(j,1,isp)) THEN
                CALL stopgm('RECPUPF','INVALID DIJ',& 
                     __LINE__,__FILE__)
             ENDIF
             l=sgpp2%pplist(i,1,isp)+1
             iv=sgpp2%pplist(i,2,isp)
             jv=sgpp2%pplist(j,2,isp)
             sgpp2%hlsg(iv,jv,l,isp)=0.5_real_8*dij
             sgpp2%hlsg(jv,iv,l,isp)=0.5_real_8*dij
          ENDDO
       ENDIF
    ENDIF
    IF (corel%tnlcc(isp)) THEN
       ! ==--------------------------------------------------------------==
       ! ..NLCC section
       IF (upf2) THEN
          endc=" "
       ELSE
          endc=">"
       ENDIF
       corei%nlcct(isp)=2
       ierr=inscan(iunit,'<PP_NLCC'//endc)
       IF (ierr.NE.0) THEN
          WRITE(6,*) ' RECPUPF| PP_NLCC SECTION NOT FOUND'
          CALL stopgm('RECPUPF',' ',& 
               __LINE__,__FILE__)
       ENDIF
       READ(iunit,END=20,err=20,fmt=*) (corecg(ir,isp),ir=1,meshv)
       corei%meshcc(isp)=meshv
       ! 
       ! vw need to set a logarithmic grid...
       CALL ckgrid(rv(1,isp),rv(meshv,isp),rcgrid(1,isp),meshv,&
            corer%clogcc(isp))
       ALLOCATE(c(meshv),fc(meshv),temp(meshv),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL tgrid(rv(1,isp),meshv,rcgrid(1,isp),meshv,corecg(1,isp),&
            meshv,c,fc,temp)
       DEALLOCATE(c,fc,temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    WRITE(6,*) ' ERROR WHILE READING PP DEFINITIONS'
    CALL stopgm('RECPUPF',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE recpupf
  ! ==================================================================
  SUBROUTINE get_section_info(iunit, header, line, lline)
    INTEGER                                  :: iunit
    CHARACTER(len=*)                         :: header
    INTEGER                                  :: lline
    CHARACTER(len=lline)                     :: line

    INTEGER                                  :: ierr, ind

! ==--------------------------------------------------------------==

    ierr=inscan(iunit,TRIM(header))
    IF (ierr.NE.0) THEN
       WRITE(6,*) ' RECPUPF| '//TRIM(header)//' SECTION NOT FOUND'
       CALL stopgm('RECPUPF (GET_SECTION_INFO)',' ',& 
            __LINE__,__FILE__)
    ENDIF
    BACKSPACE(iunit)

    ind=1
    READ(iunit,END=20,err=20,fmt='(A)')line(ind:)
    ind=LEN_TRIM(line)+1
    DO WHILE ((INDEX(line,'>').EQ.0).AND.(ind.LE.lline))
       READ(iunit,END=20,err=20,fmt='(A)')line(ind:)
       ind=LEN_TRIM(line)+1
    ENDDO
    IF (ind.GT.lline) THEN
       WRITE(6,*) ' SECTION TOO LONG! INCREASE INTERNAL BUFFER!'
       CALL stopgm('RECPUPF (GET_SECTION_INFO)',' ',& 
            __LINE__,__FILE__)
    ENDIF
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    WRITE(6,*) ' ERROR WHILE READING: GET_SECTION_INFO'
    CALL stopgm('RECPUPF',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_section_info
  ! ==================================================================
  SUBROUTINE get_key(line,field,VALUE)
    CHARACTER(len=*)                         :: line, field, VALUE

    INTEGER                                  :: ind, indf

! ==--------------------------------------------------------------==

    ind=INDEX(line,TRIM(field)//'="')
    IF (ind==0) THEN
       WRITE(6,*) ' NO FIELD IN SECTION! '
       CALL stopgm('RECPUPF (GET_KEY)',' ',& 
            __LINE__,__FILE__)
    ENDIF
    indf=INDEX(line(ind+LEN_TRIM(field)+2:),'"')
    VALUE=line(ind+LEN_TRIM(field)+2:ind+LEN_TRIM(field)+2+indf-2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_key
  ! ==================================================================
  SUBROUTINE find_element(symb,ia)
    CHARACTER(len=*)                         :: symb
    INTEGER                                  :: ia

    INTEGER                                  :: i, i1, i2

! ==--------------------------------------------------------------==

    ia=0
    DO i=1,99
       CALL xstring(elem%el(i),i1,i2)
       IF (elem%el(i)(i1:i2).EQ.symb) THEN
          ia=i
          GOTO 100
       ENDIF
    ENDDO
    WRITE(6,"(A,A)") 'FIND_ELEMENT| Could not find element ',symB
    CALL stopgm('FIND_ELEMENT','INPUT ERROR',& 
         __LINE__,__FILE__)
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE find_element
  ! ==================================================================

END MODULE recpupf_utils
