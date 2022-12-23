#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE detsp_utils
  USE array_utils,                     ONLY: array_alloc,&
                                             array_realloc
  USE atom,                            ONLY: gnl,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE clas,                            ONLY: &
       clab, clas3, clas4, clasc, clasc0, clasf, clasfold, clasv, cltyp, &
       delclasc, maxclt, tclas
  USE coor,                            ONLY: lvelini,&
                                             tau0,&
                                             velp
  USE cotr,                            ONLY: duat
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_new,&
                                             fo_old
  USE if_parallel,                     ONLY: ifparai
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: al,&
                                             bl,&
                                             ions0,&
                                             maxgau,&
                                             rcl
  USE kinds,                           ONLY: real_8
  USE mm_input,                        ONLY: addh,&
                                             gqmmm,&
                                             igqmmm,&
                                             lqmmm,&
                                             sppph
  USE mm_parallel,                     ONLY: gparai,&
                                             gparal,&
                                             irslspar,&
                                             lrslspar
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_group,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE nlcc,                            ONLY: corecg,&
                                             corei,&
                                             corer,&
                                             rcgrid
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             rgh,&
                                             wgh,&
                                             wsg
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE readsr_utils,                    ONLY: xstring
  USE sgpp,                            ONLY: mpro
  USE system,                          ONLY: cntl,&
                                             lmaxx,&
                                             maxsp,&
                                             maxsys,&
                                             nhx
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: detsp

CONTAINS

  ! ==================================================================
  SUBROUTINE detsp
    ! ==--------------------------------------------------------------==
    ! ==  DETERMINES THE NUMBER OF SPECIES                            ==
    ! ==  AND THE MAXIMUM NUMBER OF ATOMS PER SPECIES                 ==
    ! ==  DOES DYNAMIC ALLOCATION RELATED TO THESE QUANTITIES         ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'detsp'

    CHARACTER(len=127)                       :: line
    INTEGER                                  :: i, ia, ie, ierr, is, isub, &
                                                iunit, nasp

    CALL tiset(procedureN,isub)


    lqmmm%qmmm=.FALSE.
#if defined (__GROMOS)
    IF (cntl%tqmmm)CALL mm_read_qmmm_input
    CALL array_alloc(gqmmm%gr_nasp,0)
    CALL array_alloc(gqmmm%is_dummy,0)
    CALL array_alloc(gqmmm%gr_atom,[0,0])
    CALL zeroing(gqmmm%gr_atom)
    CALL zeroing(gqmmm%gr_nasp)
#else
    IF (cntl%tqmmm) CALL stopgm('DETSP','GROMOS QM/MM CODE IS MISSING',& 
         __LINE__,__FILE__)
#endif

    ! Number of DUMMY ATOMS to zero
    duat%ndat=0
    IF (paral%io_parent) THEN
       iunit = 5
       ierr=inscan(iunit,'&ATOMS')
       IF (ierr.NE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,2A,I3)') 'DETSP: COULD NOT FIND SECTION &ATOMS',&
               ' ON UNIT',iunit
          CALL stopgm('DETSP','INPUT FILE PROBLEM',& 
               __LINE__,__FILE__)
       ENDIF
       maxsys%nsx=0
       maxsys%nax=0
       clas3%nclatom=0
       clas3%ncltyp=0
10     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt='(A)') line
       IF (INDEX(line,'&END').NE.0) GOTO 30
       ! Get the number of DUMMY ATOMS defined ATOMS section (if any)
       IF (INDEX(line,'DUMMY').NE.0 .AND.&
            INDEX(line,'ATOMS').NE.0)CALL getdan(iunit)
       IF (line(1:1).NE.'*') GOTO 10

       ! ..   NEW SPECIES FOUND, now we need to disregard
       ! the file name when testing for additional flags
       CALL xstring(line,ia,ie)
       ia=ie+1
       ie=LEN(line)
       IF (ie.LT.ia) ie=ia

       IF (INDEX(line(ia:ie),'CLASSIC').NE.0) THEN
          clas3%ncltyp=clas3%ncltyp+1
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt=*) nasp
          clas3%nclatom=clas3%nclatom+nasp
       ELSE
#if defined(__GROMOS)
          maxsys%nsx=maxsys%nsx+1
          CALL array_realloc ( gqmmm%is_dummy, maxsys%nsx, keep_data=.TRUE. )
          CALL array_realloc ( gqmmm%gr_nasp, maxsys%nsx, keep_data=.TRUE. )
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt=*)
          IF (INDEX(line(ia:ie),'ADD_H').NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,END=20,err=20,fmt=*) nasp
             IF (lqmmm%qmmm)THEN
                addh%n_added_h = nasp
                gqmmm%is_dummy(maxsys%nsx)=.FALSE.
             ENDIF
             igqmmm%sp_h_added=maxsys%nsx
             gqmmm%gr_nasp(maxsys%nsx)=nasp
             ions0%na(maxsys%nsx)=nasp
             maxsys%nax=MAX(nasp,maxsys%nax)
             clas3%nclatom=clas3%nclatom+nasp
             CALL array_realloc ( gqmmm%gr_atom, [maxsys%nsx, MAX(nasp,SIZE(gqmmm%gr_atom,2))], keep_data=.TRUE. )
          ELSE
             IF (paral%io_parent) THEN
                READ(iunit,END=20,err=20,fmt=*) nasp
                CALL array_realloc ( gqmmm%gr_atom, [maxsys%nsx, MAX(nasp,SIZE(gqmmm%gr_atom,2))], keep_data=.TRUE. )
             ENDIF
             IF (lqmmm%qmmm)THEN
                gqmmm%is_dummy(maxsys%nsx)=.FALSE.
                IF (paral%io_parent)&
                     READ(iunit,END=20,err=20,fmt=*)&
                     (gqmmm%gr_atom(maxsys%nsx,i),i=1,NASP)
                IF (INDEX(line(ia:ie),'DUM').NE.0)&
                     gqmmm%is_dummy(maxsys%nsx)=.TRUE.
             ENDIF
             gqmmm%gr_nasp(maxsys%nsx)=nasp
             ions0%na(maxsys%nsx)=nasp
             maxsys%nax=MAX(nasp,maxsys%nax)
             clas3%nclatom=clas3%nclatom+nasp
          ENDIF
#else
          maxsys%nsx=maxsys%nsx+1
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt=*)
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt=*) nasp
          ions0%na(maxsys%nsx)=nasp
          maxsys%nax=MAX(nasp,maxsys%nax)
          clas3%nclatom=clas3%nclatom+nasp
#endif
       ENDIF
       GOTO 10
20     CONTINUE
       IF (paral%io_parent)&
            WRITE(6,*) ' DETSP: ERROR WHILE READING ON UNIT ',iunit
       CALL stopgm('DETSP','PARSER ERROR',& 
            __LINE__,__FILE__)
30     CONTINUE
       IF (maxsys%nsx.GT.maxsp) THEN
          CALL stopgm('DETSP','TOO MANY ATOMIC SPECIES SPECIFIED',& 
               __LINE__,__FILE__)
       ENDIF
       maxsys%ncorx = (maxsys%nsx*(maxsys%nsx+1))/2
    ENDIF

    ! ==--------------------------------------------------------------==
    CALL mp_bcast_byte(maxsys, size_in_bytes_of(maxsys),parai%io_source, parai%cp_grp)
    ! ==--------------------------------------------------------------==
#if defined(__GROMOS)
    ! ==---------------- qmmm broadcast of data ----------------------==
    IF (lqmmm%Qmmm)THEN
       IF(.NOT.paral%io_parent) THEN
          CALL array_realloc(gqmmm%gr_nasp,maxsys%nsx)
          CALL array_realloc(gqmmm%is_dummy,maxsys%nsx)
          CALL array_realloc(gqmmm%gr_atom,[maxsys%nsx,maxsys%nax])
       ENDIF
       CALL mp_bcast(gqmmm%gr_nasp,SIZE(gqmmm%gr_nasp),parai%io_source,parai%cp_grp)
       CALL mp_bcast(gqmmm%is_dummy,SIZE(gqmmm%is_dummy),parai%io_source,parai%cp_grp)
       CALL mp_bcast(gqmmm%gr_atom,SIZE(gqmmm%gr_atom),parai%io_source,parai%cp_grp)
       CALL mp_bcast(addh%n_added_h,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(igqmmm, size_in_bytes_of(igqmmm),parai%io_source,parai%cp_grp)
    ENDIF
#endif
    CALL qmmmgrps_init
    ! ==--------------------------------------------------------------==
    CALL mp_sync(parai%qmmmgrp)
    CALL mp_bcast_byte(clas3, size_in_bytes_of(clas3),parai%io_source,parai%qmmmgrp)
    ! NDAT (number of DUMMY ATOMS) 
    CALL mp_bcast(duat%ndat,parai%io_source,parai%qmmmgrp)
    ! ==--------------------------------------------------------------==
    ALLOCATE(tau0(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! To avoid floating invalid when dived by SCSALE
    CALL zeroing(tau0)!,3*maxsys%nax*maxsys%nsx)
    ! Velocities
    ALLOCATE(velp(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    ALLOCATE(lvelini(0:maxsys%nax+1, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! lvelini has to be 0:
    !$omp parallel do private(IS,IA)
    DO is=1,maxsys%nsx
       DO ia=0,maxsys%nax+1
          lvelini(ia, is)=.FALSE.
       ENDDO
    ENDDO

#if defined(__GROMOS)
    ! if doing only the mm calculation nothing here is needed.
    IF (.NOT.paral%qmnode) RETURN
#endif

    ! Analytic PP
    ALLOCATE(al(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bl(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rcl(maxgau,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Gauss-Hermit integration
    ALLOCATE(wsg(maxsys%nsx,nhx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rgh(nhx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(wgh(nhx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Mapping function to l-value
    ALLOCATE(nghtol(nhx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! FIXME deallocate missing
    ALLOCATE(nghcom(nhx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! FIXME deallocate missing
    ! Real space function and grids for Kleinman-Bylander + Potentials
    ALLOCATE(gnl(maxsys%mmaxx,maxsys%nsx,lmaxx*mpro),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rps(maxsys%mmaxx,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rw(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rv(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vr(maxsys%mmaxx,maxsys%nsx,lmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Non Linear Core Correction
    ALLOCATE(rcgrid(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(corecg(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rcgrid)!,maxsys%mmaxx*maxsys%nsx)
    CALL zeroing(corecg)!,maxsys%mmaxx*maxsys%nsx)
    !$omp parallel do private(IS)
    DO is=1,maxsp
       corer%anlcc(is)=0.0_real_8
       corer%bnlcc(is)=0.0_real_8
       corer%enlcc(is)=0.0_real_8
       corer%clogcc(is)=0.0_real_8
       corei%nlcct(is)=0
       corei%meshcc(is)=0
    ENDDO
    ! Classical atoms
    IF (clas3%ncltyp.GT.0) THEN
       tclas=.TRUE.
       clas3%ncltyp=clas3%ncltyp+maxsys%nsx
       ALLOCATE(clasc(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(clasv(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(clasf(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(clasc0(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(clasfold(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(delclasc(3,clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(clasc)!,3*clas3%nclatom)
       CALL zeroing(clasv)!,3*clas3%nclatom)
       CALL zeroing(clasf)!,3*clas3%nclatom)
       CALL zeroing(clasc0)!,3*clas3%nclatom)
       CALL zeroing(delclasc)!,3*clas3%nclatom)
       CALL zeroing(clasfold)!,3*clas3%nclatom)
       ALLOCATE(cltyp(3*clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (clas3%ncltyp.GT.maxclt) CALL stopgm('DETSP',&
            'PARAMETER MAXCLT TOO SMALL',& 
            __LINE__,__FILE__)
       DO i=1,maxclt
          clas3%is_qm(i)=0
          clab(i)='    '
          clas4%clmas(i)=0._real_8
          clas3%ncrang(1,i)=0
          clas3%ncrang(2,i)=0
       ENDDO
    ELSE
       tclas=.FALSE.
    ENDIF
    IF (paral%parent) CALL prmem('     DETSP')
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE detsp
  ! ==================================================================
  SUBROUTINE qmmmgrps_init
    ! ==--------------------------------------------------------------==
    ! ==  Setup the groups for parallel MM and QMMM                   ==
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: grouplist(parai%cp_nproc), ip

! ==--------------------------------------------------------------==

#ifndef __GROMOS
    ! not defined gromos
    ! for compatibility with QM/MM
    parai%qmmmgrp=parai%cp_grp
    paral%qmnode=.TRUE.
    RETURN
#endif
    ! set to defaults for serial runs
    parai%qmmmnproc=parai%cp_nproc
    parai%qmmmsource=parai%io_source !check
    parai%qmmmgrp=parai%cp_grp
    gparai%mmnproc=parai%cp_nproc

    ifparai%ifgrp=parai%qmmmgrp
    paral%qmmmparent=paral%io_parent !check
    IF (lqmmm%qmmm) THEN
       ! set up mm group
       gparal%mmnode=.TRUE.
       gparai%mmmepos=parai%cp_me
       gparal%mmparent=paral%qmmmparent !check
       gparai%mmsource=parai%qmmmsource !check
       gparai%mmnproc=parai%qmmmnproc
       gparai%mmgrp=parai%qmmmgrp
       ! set up real space nodes for mm 
       lrslspar%rsnode=gparal%mmnode
       irslspar%nrsnodes=gparai%mmnproc
       lrslspar%rsparent=gparal%mmparent
       irslspar%rspos=gparai%mmmepos
       irslspar%rsgrp=gparai%mmgrp
       ! set up lattice sum nodes for mm
       lrslspar%lsnode=gparal%mmnode
       irslspar%nlsnodes=gparai%mmnproc
       lrslspar%lsparent=gparal%mmparent
       irslspar%lspos=gparai%mmmepos
       irslspar%lsgrp=gparai%mmgrp
    ENDIF
    IF (parai%nproc.GT.1 .AND. sppph%mm_split) THEN
       ! split off last node for the lattice sum part   
       DO ip=1,parai%cp_nproc
          grouplist(ip)=0
       ENDDO
       irslspar%nrsnodes=parai%qmmmnproc-1
       irslspar%nlsnodes=1
       DO ip=1,parai%cp_nproc
          IF(ip.LE.irslspar%nrsnodes) THEN
             IF (parai%cp_me.EQ.(ip-1)) THEN
                lrslspar%rsnode=.TRUE.
                irslspar%rspos=parai%cp_me
                lrslspar%lsnode=.FALSE.
             END IF
          ELSE
             IF(parai%cp_me.EQ.(ip-1)) THEN
                lrslspar%lsnode=.TRUE.
                lrslspar%rsnode=.FALSE.
                irslspar%lspos=(ip-1)-irslspar%nrsnodes
                lrslspar%lsparent=.FALSE.
                IF(irslspar%lspos.EQ.0) lrslspar%lsparent=.TRUE.
             ENDIF
          END IF
          grouplist(ip)=ip-1
       END DO
       ! the rs and ls groups are not actually used at the moment
       ! creation is done for completeness just in case it becomes 
       ! nescessary.
       irslspar%rsgrp=HUGE(0)
       irslspar%lsgrp=HUGE(0)
       CALL mp_group(irslspar%nrsnodes,grouplist(1:irslspar%nrsnodes),irslspar%rsgrp,parai%qmmmgrp)
       !     grouplist(1)=nproc-1
       CALL mp_group(irslspar%nlsnodes,grouplist(irslspar%nrsnodes+1),irslspar%lsgrp,parai%qmmmgrp)
    ENDIF
    ! split nodes to planes for the electrostatic interactions on the 
    ! has to be done after loadpa has been executed      
    ifparai%ifgrp=-1

    CALL mp_sync(parai%qmmmgrp)
    RETURN
  END SUBROUTINE qmmmgrps_init
  ! ==================================================================
  SUBROUTINE getdan(iunit)
    ! Reading the number of DUMMY ATOMS from &ATOMS section.
    ! The size of the dummy atoms is necessary for allocating their
    ! indices with NAT_cpmd array. But type of DUMMY ATOMS is read later.
    INTEGER                                  :: iunit

! Reading the next line after DUMMY ATOMS (number of dummies)

    IF (paral%io_parent)&
         READ(iunit,END=21,err=21,fmt=*)duat%ndat
    RETURN
21  CONTINUE
    CALL stopgm('DETSP','ERROR READING DUMMY ATOMS',& 
         __LINE__,__FILE__)
  END SUBROUTINE getdan

  ! ==================================================================
  SUBROUTINE my_syncfile(fname,fromnode,tonode,comm)
    ! transfer the text file 'fname' to all nodes that done have
    ! it locally. e.g. required in the NEC-SX.
    ! ==================================================================
    ! ARGS
    CHARACTER(len=*)                         :: fname
    INTEGER                                  :: fromnode, tonode, comm

    INTEGER, PARAMETER                       :: iunit = 72 , maxline = 120 

    CHARACTER(len=maxline)                   :: line
    INTEGER                                  :: ia, ie, inmsg
    LOGICAL                                  :: ferr, fexists

! try to figure out whether we have to receive the file.

    IF (fromnode.EQ.tonode) RETURN

    CALL mp_sync(comm)
    fexists=.FALSE.
    inmsg=0
    CALL xstring(fname,ia,ie)
    IF (paral%io_parent)&
         INQUIRE(file=fname(ia:ie),exist=fexists)
    IF (fexists) inmsg=1
    CALL mp_bcast(inmsg,tonode,comm)
    IF (inmsg.EQ.1) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'FILE ',fname(ia:ie),' AVAILABLE'
       RETURN
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) 'SYNCING FILE ', fname(ia:ie),&
            ' FROM NODE',fromnode,' TO NODE', tonode
    ENDIF


    ! send file       
    IF (parai%me.EQ.fromnode) THEN
       ferr=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(iunit,fname(ia:ie),fo_old,ferr)
       IF (ferr) GOTO 200
100    CONTINUE
       IF (paral%io_parent)&
            READ(iunit,fmt='(A)',err=200,END=200) line
       inmsg=1
       CALL mp_send(inmsg,tonode,0,comm)
       CALL mp_send(line,tonode,1,comm)
       GOTO 100
200    CONTINUE
       inmsg=0
       CALL mp_send(inmsg,tonode,0,comm)
       IF (paral%io_parent)&
            CALL fileclose(iunit)
       RETURN
    ENDIF
    IF (parai%me.EQ.tonode) THEN
       IF (paral%io_parent)&
            CALL fileopen(iunit,fname(ia:ie),fo_new,ferr)
300    CONTINUE
       inmsg=0
       CALL mp_recv(inmsg,fromnode,0,comm)
       IF (inmsg.EQ.0) THEN
          IF (paral%io_parent)&
               CALL fileclose(iunit)
          RETURN
       ENDIF
       CALL mp_recv(line,fromnode,1,comm)
       CALL xstring(line,ia,ie)
       IF (paral%io_parent)&
            WRITE(iunit,'(A)') line(1:ie)
       GOTO 300
    ENDIF
    RETURN
  END SUBROUTINE my_syncfile
  ! ==================================================================

END MODULE detsp_utils
