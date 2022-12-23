MODULE ratom_utils
  USE atom,                            ONLY: atom_common,&
                                             gnl,&
                                             patom1,&
                                             patom2,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE atwf,                            ONLY: atchg
  USE clas,                            ONLY: clab,&
                                             clas3,&
                                             clas4,&
                                             clasc,&
                                             cltyp,&
                                             tclas
  USE cnst_dyn,                        ONLY: &
       ekincv, icv_spin, imeta, inter_hill, inter_hill_max, kharm, lmeta, &
       ncolvar, nsubsys, rcc0, rmeta, vharm
  USE coninp_utils,                    ONLY: coninp
  USE coor,                            ONLY: tau0,&
                                             velp
  USE cotr,                            ONLY: &
       boad, cnpar, cnsval, cnsval_dest, cotc0, cotr007, duat, grate, gsrate, &
       lskcor, lskptr, maxdum, ntcnst, ntrest, resfor, respar, respos, &
       resval, resval_dest
  USE dpot,                            ONLY: dpot_mod
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ghermit_utils,                   ONLY: ghermit
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: tshl
  USE meta_cell_utils,                 ONLY: meta_cell_inp
  USE meta_colvar_inp_utils,           ONLY: meta_colvar_inp
  USE meta_dyn_def_utils,              ONLY: meta_dyn_def
  USE mm_dimmod,                       ONLY: mmdim
  USE mm_input,                        ONLY: g96_vel,&
                                             lqmmm
  USE movi,                            ONLY: imtyp
  USE mp_interface,                    ONLY: mp_bcast
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nlcc,                            ONLY: corel
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             nlps_com,&
                                             rgh,&
                                             wgh
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE recpnew_utils,                   ONLY: allelec,&
                                             recpeam,&
                                             recpnew
  USE recpupf_utils,                   ONLY: recpupf
  USE rmas,                            ONLY: rmass
  USE symm,                            ONLY: iun,&
                                             symmt
  USE system,                          ONLY: maxsp,&
                                             maxsys
  USE tst2min_inp_utils,               ONLY: tst2min_inp
  USE velocitinp_utils,                ONLY: velocitinp
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ratom

CONTAINS

  ! ==================================================================
  SUBROUTINE ratom
    ! ==--------------------------------------------------------------==
    ! ==  Reads Section &ATOMS from Unit iunit                        ==
    ! ==--------------------------------------------------------------==
    ! ==  DESCRIPTION OF INPUT                                        ==
    ! ==                                                              ==
    ! ==  &ATOMS                                                      ==
    ! ==                                                              ==
    ! ==  *ECPNAME [OPTIONS]                                          ==
    ! ==  {LMAX LOC SKIP , LMAX=L [LOC=L,SKIP=L]}                     ==
    ! ==  NATOMS                                                      ==
    ! ==    X     Y     Z                                             ==
    ! ==    X     Y     Z                                             ==
    ! ==    ....                                                      ==
    ! ==                                                              ==
    ! ==  ATOMIC CHARGES                                              ==
    ! ==   c1 c2 ... cnsp                                             ==
    ! ==                                                              ==
    ! ==  MOVIE TYPE                                                  ==
    ! ==   nt1 nt2 ...                                                ==
    ! ==                                                              ==
    ! ==  DUMMY ATOMS                                                 ==
    ! ==   ndat                                                       ==
    ! ==   TYPEx   x    y    z                                        ==
    ! ==   TYPEy   n    n1 n2 n3 ...                                  ==
    ! ==  ....                                                        ==
    ! ==                                                              ==
    ! ==  GENERATE COORDINATES                                        ==
    ! ==   iun1  iun2  ...  iunn                                      ==
    ! ==                                                              ==
    ! ==  ISOTOPE                                                     ==
    ! ==   mass1                                                      ==
    ! ==   ....                                                       ==
    ! ==                                                              ==
    ! ==  CONFINEMENT POTENTIAL                                       ==
    ! ==   a1 rc1                                                     ==
    ! ==   ....                                                       ==
    ! ==                                                              ==
    ! ==  CHANGE BONDS                                                ==
    ! ==    NBOAD                                                     ==
    ! ==    I  J   +/- 1                                              ==
    ! ==    ....                                                      ==
    ! ==                                                              ==
    ! ==  CONSTRAINTS                                                 ==
    ! ==    ....                                                      ==
    ! ==  END CONSTRAINTS                                             ==
    ! ==                                                              ==
    ! ==  META DYNAMICS {COLLECTIVE VARIABLE; FIX ; MULTI}            ==
    ! ==                                                              ==
    ! ==  SADDLE POINT                                                ==
    ! ==                                                              ==
    ! ==  VELOCITIES                                                  ==
    ! ==    ....                                                      ==
    ! ==  END VELOCITIES                                              ==
    ! ==                                                              ==
    ! ==  &END                                                        ==
    ! ==   OPTIONS: GAUSS-HERMIT=n  |  KLEINMAN-BYLANDER              ==
    ! ==            NLCC                                              ==
    ! ==            RAGGIO=ragg                                       ==
    ! ==            BINARY | FORMATTED                                ==
    ! ==            UPF                                               ==
    ! ==            NONE                                              ==
    ! ==            CLASSIC                                           ==
    ! ==            _EAM_                                             ==
    ! ==            FRAC                                              ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'ratom'

    CHARACTER(len=127)                       :: line
    CHARACTER(len=40)                        :: ecpnam
    INTEGER :: i, i1, i11, i2, i3, i4, ia, iat, ie, ierr, igh, il, ilm, iloc, &
      iout, is, iskip, iunit, iv, j, k, lmaxg, nasp, ncla, ncls, nghl, nghln, &
      NSX_q
    LOGICAL                                  :: erread, tfrac, tnone, tupf
    REAL(real_8)                             :: nel_eam, nel_sum, raggnew

! ==--------------------------------------------------------------==
! Initialization

    IF (lqmmm%qmmm)THEN
       NSX_q=mmdim%nspq
       maxsys%nsx=mmdim%nspm
    ELSE
       NSX_q=maxsys%nsx
       mmdim%nspm=maxsys%nsx
    ENDIF
    ! 
    patom1%pconf=.FALSE.
    ! 
    ALLOCATE(lskcor(3,maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lskptr(3,maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(atchg(NSX_q),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(atchg)!,NSX_q)
    CALL zeroing(lskptr)!,3*maxsys%nax*maxsys%nsx)
    DO i=1,maxsys%nax*maxsys%nsx
       DO j=1,3
          lskcor(j,i)=1
       ENDDO
    ENDDO
    ! Velocities
    ! McB
    IF ( g96_vel%ntx_vel.NE.1 ) THEN
       CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    ! McB
    ! Real space function and grids for Kleinman-Bylander + Potentials
    CALL zeroing(gnl)!,maxsys%mmaxx*NSX_q*lmaxx)
    CALL zeroing(rps)!,maxsys%mmaxx*NSX_q*lmaxx)
    CALL zeroing(vr)!,maxsys%mmaxx*NSX_q*lmaxx)
    CALL zeroing(rw)!,maxsys%mmaxx*NSX_q)
    CALL zeroing(rv)!,maxsys%mmaxx*NSX_q)
    ! ==--------------------------------------------------------------==
    IF (.NOT.paral%io_parent) GOTO 9999
    ! ==--------------------------------------------------------------==
    iunit = 5
    ierr=inscan(iunit,'&ATOMS')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A,I3)')&
            ' RATOM| COULD NOT FIND SECTION &ATOMS ON  UNIT',iunit
       CALL stopgm('RATOM',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    DO i=1,maxsp
       imtyp(i)=0
       dpot_mod%team(i)=.FALSE.
       dpot_mod%tkb(i)=.FALSE.
    ENDDO
    cotc0%lfcom=.FALSE.
    symmt%tgenc=.FALSE.
    ions1%nsp=0
    crge%nel=0._real_8
    nel_eam=0.0_real_8
    nel_sum=0.0_real_8
    duat%ndat  = 0
    duat%ndat1 = 0
    duat%ndat2 = 0
    duat%ndat3 = 0
    duat%ndat4 = 0
    cotc0%mcnstr = 0
    cotr007%mrestr = 0
    ncls = 0
    ncla = 0
    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    IF (INDEX(line,'&END').NE.0) GOTO 30
    IF (line(1:1).NE.'*') GOTO 10
    ! NEW SPECIES found
    CALL xstring(line,ia,ie)
    ecpnam=line(2:ie)
    ia=ie+1
    ie=LEN(line)
    IF (ie.LT.ia) ie=ia

    ncls=ncls+1
    IF (INDEX(line(ia:ie),'CLASSIC').NE.0) THEN
       i1=ia-1+INDEX(line(ia:ie),'LABEL')
       IF (i1.GE.ia) THEN
          i11=i1+INDEX(line(i1:127),'=')
          clab(ncls)(1:4)=line(i11:i11+3)
       ELSE
          CALL stopgm('RATOM','NO LABEL FOR CLASSICAL ATOM',& 
               __LINE__,__FILE__)
       ENDIF
       i1=ia-1+INDEX(line(ia:ie),'MASS')
       IF (i1.GE.ia) THEN
          i11=i1+INDEX(line(i1:127),'=')
          CALL readsr(line(i11:127),1,iout,clas4%clmas(ncls),erread)
       ELSE
          CALL stopgm('RATOM','NO MASS FOR CLASSICAL ATOM',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt=*) nasp
       clas3%ncrang(1,ncls)=ncla+1
       ncla=ncla+nasp
       clas3%ncrang(2,ncls)=ncla
       DO i=clas3%ncrang(1,ncls),clas3%ncrang(2,ncls)
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt=*) (clasc(k,i),k=1,3)
          cltyp(i)=ncls
       ENDDO
       clas3%is_qm(ncls)=0
    ELSEIF (INDEX(line(ia:ie),'_EAM_').NE.0) THEN
       ions1%nsp=ions1%nsp+1
       dpot_mod%team(ions1%nsp)=.TRUE.
       tnone=.FALSE.
       pslo_com%tvan(ions1%nsp)=.FALSE.
       pslo_com%tbin(ions1%nsp)=.FALSE.
       corel%tnlcc(ions1%nsp)=.FALSE.
       dpot_mod%tkb(ions1%nsp)=.FALSE.
       nghl=0
       dpot_mod%lmax(ions1%nsp)=0
       dpot_mod%lloc(ions1%nsp)=0
       dpot_mod%lskip(ions1%nsp)=1
       nlps_com%ngh(ions1%nsp)=0
       CALL recpeam(ions1%nsp,ecpnam)
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt='(A)') line
       ! ATOMIC COORDINATES
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt=*) ions0%na(ions1%nsp)
       DO i=1,ions0%na(ions1%nsp)
          READ(iunit,END=20,err=20,fmt='(A)') line
          ia=1
          DO k=1,3
             CALL readsr(line,ia,ie,tau0(k,i,ions1%nsp),erread)
             IF (erread) GOTO 20
             ia=ie+1
          ENDDO
       ENDDO
    ELSE
       ions1%nsp=ions1%nsp+1
       raggnew=-1._real_8
       tnone=.FALSE.
       tfrac=.FALSE.
       IF (INDEX(line(ia:ie),'NONE').NE.0) tnone=.TRUE.
       IF (INDEX(line(ia:ie),'FRAC').NE.0) tfrac=.TRUE.
       tupf=.FALSE.
       IF (INDEX(line(ia:ie),'UPF').NE.0) tupf=.TRUE.
       pslo_com%tvan(ions1%nsp)=.FALSE.
       pslo_com%tbin(ions1%nsp)=.FALSE.
       IF (INDEX(line(ia:ie),'FORMAT').NE.0) pslo_com%tvan(ions1%nsp)=.TRUE.
       IF (INDEX(line(ia:ie),'BINARY').NE.0) THEN
          IF (pslo_com%tvan(ions1%nsp))&
               CALL stopgm('RATOM','BINARY or FORMATTED uspp file?',& 
               __LINE__,__FILE__)
          pslo_com%tvan(ions1%nsp)=.TRUE.
          pslo_com%tbin(ions1%nsp)=.TRUE.
       ENDIF
       corel%tnlcc(ions1%nsp)=.FALSE.
       IF (INDEX(line(ia:ie),'NLCC').NE.0) corel%tnlcc(ions1%nsp)=.TRUE.
       nghl=20
       i1=ia-1+INDEX(line(ia:ie),'GAUSS')
       i2=ia-1+INDEX(line(ia:ie),'KLEINMAN')
       IF (i1.GE.ia) THEN
          dpot_mod%tkb(ions1%nsp)=.FALSE.
          i11=i1+INDEX(line(i1:127),'=')
          IF (i11.GT.i1) CALL readsi(line(i11:127),1,iout,nghl,erread)! TODO check if io_parent is needed
       ENDIF
       IF (i2.GT.ia) dpot_mod%tkb(ions1%nsp)=.TRUE.
       IF (i1.EQ.(ia-1) .AND. i2.EQ.(ia-1)) THEN
          dpot_mod%tkb(ions1%nsp)=.FALSE.
       ENDIF
       i3=ia-1+INDEX(line(ia:ie),'RAGGIO')
       IF (i3.GE.ia) THEN
          i4=i3+INDEX(line(i3:127),'=')
          CALL readsr(line(i4:127),1,iout,raggnew,erread)
       ENDIF
       CALL lread(iunit,dpot_mod%lmax(ions1%nsp),iloc,iskip)
       dpot_mod%lloc(ions1%nsp)=dpot_mod%lmax(ions1%nsp)+1
       dpot_mod%lskip(ions1%nsp)=dpot_mod%lmax(ions1%nsp)+2
       nlps_com%ngh(ions1%nsp)=0
       IF (dpot_mod%tkb(ions1%nsp)) THEN
          IF (iloc.GE.0.AND.iloc.LE.dpot_mod%lmax(ions1%nsp)) dpot_mod%lloc(ions1%nsp)=iloc+1
          IF (iskip.GE.0.AND.iskip.LE.dpot_mod%lmax(ions1%nsp)) dpot_mod%lskip(ions1%nsp)=iskip+1
       ELSEIF (i1.NE.0) THEN
          CALL ghermit(nghl,rgh(1,ions1%nsp),wgh(1,ions1%nsp))
          nghln=0
          DO i=1,nghl
             IF (wgh(i,ions1%nsp).GT.1.e-10_real_8) THEN
                nghln=nghln+1
                wgh(nghln,ions1%nsp)=wgh(i,ions1%nsp)
                rgh(nghln,ions1%nsp)=rgh(i,ions1%nsp)
             ENDIF
          ENDDO
          nghl=nghln
          nlps_com%rmaxn(ions1%nsp)=REAL(nghl,kind=real_8)
          IF (iloc.GE.0.AND.iloc.LE.dpot_mod%lmax(ions1%nsp)) dpot_mod%lloc(ions1%nsp)=iloc+1
          IF (iskip.GE.0.AND.iskip.LE.dpot_mod%lmax(ions1%nsp)) dpot_mod%lskip(ions1%nsp)=iskip+1
          lmaxg=dpot_mod%lmax(ions1%nsp)+1
          IF (lmaxg.EQ.dpot_mod%lloc(ions1%nsp)) lmaxg=dpot_mod%lmax(ions1%nsp)
          iv=0
          DO il=1,lmaxg
             DO ilm=1,2*il-1
                DO igh=1,nghl
                   iv=iv+1
                   nghtol(iv,ions1%nsp)=il-1
                ENDDO
             ENDDO
          ENDDO
          nlps_com%ngh(ions1%nsp)=iv
       ENDIF
       dpot_mod%lmax(ions1%nsp)=dpot_mod%lmax(ions1%nsp)+1
       IF (tclas) THEN
          i1=ia-1+INDEX(line(ia:ie),'LABEL')
          IF (i1.GE.ia) THEN
             i11=i1+INDEX(line(i1:127),'=')
             clab(ncls)(1:4)=line(i11:i11+3)
          ELSE
             CALL stopgm('RATOM','NO LABEL FOR QUANTUM ATOM',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! PSEUDOPOTENTIAL
       IF (lqmmm%qmmm)maxsys%nsx=mmdim%nspq
       IF (tnone) THEN
          CALL allelec(ions1%nsp,ecpnam)
       ELSEIF (tupf) THEN
          CALL recpupf(ions1%nsp,ecpnam)
       ELSE
          CALL recpnew(ions1%nsp,ecpnam)
       ENDIF
       IF (lqmmm%qmmm)maxsys%nsx=mmdim%nspm
       IF (raggnew.GT.0._real_8) raggio(ions1%nsp)=raggnew
       ! ATOMIC COORDINATES
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt=*) iat
       IF (iat.NE.ions0%na(ions1%nsp)) CALL stopgm('RATOM',&
            'INCONSISTENT READ OF &ATOMS SECTION',& 
            __LINE__,__FILE__)
       IF (.NOT.dpot_mod%team(ions1%nsp)) THEN
          ! FRACTIONAL CORE CHARGE?
          IF (tfrac) THEN
             nel_sum=nel_sum+REAL(ions0%na(ions1%nsp),kind=real_8)*ions0%zv(ions1%nsp)
          ELSE
             crge%nel=crge%nel+REAL(ions0%na(ions1%nsp)*NINT(ions0%zv(ions1%nsp)),kind=real_8)
          ENDIF
       ELSE
          nel_eam=nel_eam+REAL(ions0%na(ions1%nsp),kind=real_8)*ions0%zv(ions1%nsp)
       ENDIF
       IF (.NOT. lqmmm%qmmm)THEN
          DO i=1,ions0%na(ions1%nsp)
             IF (paral%io_parent)&
                  READ(iunit,END=20,err=20,fmt='(A)') line
             ia=1
             DO k=1,3
                CALL readsr(line,ia,ie,tau0(k,i,ions1%nsp),erread)
                IF (erread) GOTO 20
                ia=ie+1
             ENDDO
          ENDDO
       ENDIF
       IF (tclas) THEN
          clas3%ncrang(1,ncls)=ncla+1
          ncla=ncla+ions0%na(ions1%nsp)
          clas3%ncrang(2,ncls)=ncla
          DO i=clas3%ncrang(1,ncls),clas3%ncrang(2,ncls)
             cltyp(i)=ncls
          ENDDO
          clas3%is_qm(ncls)=ions1%nsp
       ENDIF
    ENDIF
    GOTO 10
    ! ==--------------------------------------------------------------==
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' RATOM| ERROR WHILE READING ON UNIT ',iunit
    CALL stopgm('RATOM',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
30  CONTINUE
    ! ==--------------------------------------------------------------==
    ! HANDLE FRACTIONAL CORE CHARGES
    IF (ABS(REAL(NINT(nel_sum),kind=real_8)-nel_sum).GT.1.0e-10_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,/,A,F20.10,/,A,I9,/)')&
            ' RATOM| WARNING! SUM OF CORE CHARGES IS NOT INTEGER',&
            ' RATOM| SUM OF FRACTIONAL CORE CHARGES IS:', nel_sum,&
            ' RATOM| NUMBER OF ELECTRONS ADDED IS:     ', NINT(nel_sum)
    ENDIF
    crge%nel=crge%nel+REAL(NINT(nel_sum),kind=real_8)
    ! ==--------------------------------------------------------------==
    cotc0%nboad=0
    ierr=inscan(iunit,'&ATOMS')
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ==  Metadynamics Defaults
    ![ EXACT FACTORIZATION
    IF (.NOT.tshl%txfmqc) THEN
       CALL meta_dyn_def
    ENDIF
    !] EXACT FACTORIZATION
    ! ==--------------------------------------------------------------==
40  CONTINUE

    IF (paral%io_parent)&
         READ(iunit,END=50,err=50,fmt='(A)') line
    IF (INDEX(line,'&END').NE.0) GOTO 60
    IF (INDEX(line,'ISOTOPE').NE.0) THEN
       ! READ ISOTOPIC MASSES
       DO i=1,ions1%nsp
          IF (paral%io_parent)&
               READ(iunit,END=50,err=50,fmt=*) rmass%pma0(i)
       ENDDO
       GOTO 40
    ELSEIF(INDEX(line,'ATOMIC').NE.0 .AND.&
         INDEX(line,'CHARGE').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,END=50,err=50,fmt=*) (atchg(ia),ia=1,NSX_q)
       GOTO 40
    ELSEIF(INDEX(line,'CONFIN').NE.0 .AND.&
         INDEX(line,'POTENT').NE.0) THEN
       ! Confinement potentials
       patom1%pconf=.TRUE.
       DO i=1,ions1%nsp
          IF (paral%io_parent)&
               READ(iunit,END=50,err=50,fmt=*) patom2%palpha(i),patom2%prc(i)
       ENDDO
       GOTO 40
    ELSEIF(INDEX(line,'MOVIE').NE.0 .AND.&
         INDEX(line,'TYPE').NE.0) THEN
       ! CHANGE DEFAULT MOVIE TYPE
       IF (paral%io_parent)&
            READ(iunit,END=50,err=50,fmt=*) (imtyp(is),is=1,ions1%nsp)
       GOTO 40
    ELSEIF(INDEX(line,'GENER').NE.0 .AND.&
         INDEX(line,'COORD').NE.0) THEN
       ! GENERATE ATOMIC COORDINATES FROM SYMMETRY UNIQUE ONES
       IF (paral%io_parent)&
            READ(iunit,END=50,err=50,fmt=*) (iun(is),is=1,ions1%nsp)
       symmt%tgenc=.TRUE.
       GOTO 40
    ELSEIF(INDEX(line,'DUMMY').NE.0 .AND.&
         INDEX(line,'ATOM').NE.0) THEN
       ! READ DUMMY ATOM DEFINITION
       IF (paral%io_parent)&
            READ(iunit,END=50,err=50,fmt=*) duat%ndat
       IF (duat%ndat.GT.maxdum) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' TOO MANY DUMMY ATOMS, CHANGE MAXDUM',&
               ' IN INCLUDE FILE COTR.INC'
          CALL stopgm('RATOM',' ',& 
               __LINE__,__FILE__)
       ENDIF
       DO iat=1,duat%ndat
          IF (paral%io_parent)&
               READ(iunit,END=50,err=50,fmt='(A)') line
          CALL rdumat(line,iat)
       ENDDO
       GOTO 40
    ELSEIF(INDEX(line,'CHANGE').NE.0 .AND.&
         INDEX(line,'BONDS').NE.0) THEN
       ! Add or delete bonds for the empirical Hessian
       IF (paral%io_parent) READ(iunit,END=50,err=50,fmt=*) cotc0%nboad
       ALLOCATE(boad(3,cotc0%nboad),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,cotc0%nboad
          IF (paral%io_parent)&
               READ(iunit,END=50,err=50,fmt=*) boad(1,i),boad(2,i),&
               boad(3,i)
       ENDDO
       GOTO 40
    ELSEIF (INDEX(line,'CONSTRAINT').NE.0) THEN
       ! Input constraints on the geometry
       CALL coninp(iunit)
       GOTO 40
    ELSEIF(INDEX(line,'META').NE.0 .AND.&
         INDEX(line,'DYNAMICS').NE.0) THEN
       IF (INDEX(line,'COLLE').NE.0 .AND.&
            INDEX(line,'VAR').NE.0) THEN
          ! Input for the meta dynamics in the space of the Collective Variables
          CALL meta_colvar_inp(iunit)
       ELSEIF (INDEX(line,'MULTI').NE.0) THEN
          ! Input for the meta dynamics with multiple walkers
          IF (INDEX(line,'WALKERS').NE.0)THEN
             tmw=.TRUE.
             il=INDEX(line,'NW=')
             IF (il.NE.0) THEN
                ia = il+3
                CALL readsi(line,ia,ie,mwi%nwalk,erread)
                IF (mwi%nwalk > parai%cp_nproc) THEN
                   WRITE(6,*)'WARNING: NWALK must be at most CP_NPROC, ' &
                        // 'reset NWALK = CP_NPROC'
                   mwi%nwalk = parai%cp_nproc
                ENDIF
             ELSE
                CALL stopgm('RATOM','NUMBER OF WALKERS NOT SPECIFIED',& 
                     __LINE__,__FILE__)
             ENDIF
             CALL meta_colvar_inp(iunit)
          ELSE
             ! Input for the meta dynamics in the space of the Collective Variables
             ! more than one set of CV is activated
             lmeta%tmulti = .TRUE.
             il=INDEX(line,'NS=')
             IF (il.NE.0) THEN
                ia = il+3
                CALL readsi(line,ia,ie,nsubsys,   erread)
             ENDIF
             CALL meta_colvar_inp(iunit)
          ENDIF
       ELSEIF (INDEX(line,'CELL').NE.0) THEN
          ! Input for the metadynamics of the cell (only if cntl%tprcp)
          CALL meta_cell_inp(iunit,line)
       ENDIF
       GOTO 40
    ELSEIF(INDEX(line,'SADDLE').NE.0 .AND.&
         INDEX(line,'POINT') .NE. 0) THEN
       CALL tst2min_inp(iunit)
       GOTO 40
    ELSEIF (INDEX(line,'VELOC').NE.0) THEN
       ! Input velocities
       CALL velocitinp(iunit)
       GOTO 40
    ELSE
       GOTO 40
    ENDIF
50  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' RATOM| ERROR WHILE READING KEYWORDS'
    IF (paral%io_parent)&
         WRITE(6,*) line
    CALL stopgm('RATOM',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
60  CONTINUE
    ! ==--------------------------------------------------------------==
    IF (ABS(nel_eam).GT.1.e-10_real_8) THEN
       CALL stopgm('RATOM','EAM charges do not add up to zero',& 
            __LINE__,__FILE__)
    ENDIF
9999 CONTINUE
    ! ==--------------------------------------------------------------==
    !cnst_dyn.mod.F90
    CALL mp_bcast(ekincv,parai%io_source,parai%cp_grp)
    CALL mp_bcast(vharm,parai%io_source,parai%cp_grp)
    ! COTR
    CALL mp_bcast_byte(cotc0, size_in_bytes_of(cotc0),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(cotr007, size_in_bytes_of(cotr007),parai%io_source,parai%cp_grp)
    ! TODO refactor all such common block implicit broadcasts by explicit: 
    CALL mp_bcast(lskcor,SIZE(lskcor),parai%io_source,parai%cp_grp)
    CALL mp_bcast(lskptr,SIZE(lskptr),parai%io_source,parai%cp_grp)
    IF (cotc0%mcnstr.GT.0) THEN 
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(ntcnst(6,cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cnsval(cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cnpar(2,cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(grate(cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cnsval_dest(cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! avoid 'not allocated' runtime error
          IF (.NOT.ALLOCATED(ntcnst)) THEN
             ALLOCATE(ntcnst(1,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cnsval(1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cnpar(1,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(grate(1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cnsval_dest(1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF ! if (.not.io_parent)
       CALL mp_bcast(ntcnst,SIZE(ntcnst),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cnsval,SIZE(cnsval),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cnpar,SIZE(cnpar),parai%io_source,parai%cp_grp)
       CALL mp_bcast(grate,SIZE(grate),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cnsval_dest,SIZE(cnsval_dest),parai%io_source,parai%cp_grp)
    ENDIF ! if (mcnstr.gt.0)
    IF (cotr007%mrestr.GT.0) THEN
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(ntrest(6,cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(resval(cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(resfor(cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(respar(2,cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gsrate(cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(resval_dest(cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(respos(3,cotr007%mrestr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__) ! cmb-kk
       ENDIF
       CALL mp_bcast(ntrest,SIZE(ntrest),parai%io_source,parai%cp_grp)
       CALL mp_bcast(resval,SIZE(resval),parai%io_source,parai%cp_grp)
       CALL mp_bcast(resfor,SIZE(resfor),parai%io_source,parai%cp_grp)
       CALL mp_bcast(gsrate,SIZE(gsrate),parai%io_source,parai%cp_grp)
       CALL mp_bcast(resval_dest,SIZE(resval_dest),parai%io_source,parai%cp_grp)
       CALL mp_bcast(respar,SIZE(respar),parai%io_source,parai%cp_grp)
       CALL mp_bcast(respos,SIZE(respos),parai%io_source,parai%cp_grp)
    ENDIF ! (mrestr.gt.0)
    ! DUMMY ATOMS
    CALL mp_bcast_byte(duat, size_in_bytes_of(duat),parai%io_source,parai%cp_grp)
    ! ATOM
    CALL mp_bcast_byte(atom_common, size_in_bytes_of(atom_common),parai%io_source,parai%cp_grp)
    ! NLPS
    CALL mp_bcast_byte(nlps_com, size_in_bytes_of(nlps_com),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rgh,SIZE(rgh),parai%io_source,parai%cp_grp)
    CALL mp_bcast(wgh,SIZE(wgh),parai%io_source,parai%cp_grp)
    ! NGH/NGHTOL
    CALL mp_bcast(nghtol,SIZE(nghtol),parai%io_source,parai%cp_grp)
    ! NGH/NGHCOM 
    CALL mp_bcast(nghcom,SIZE(nghcom),parai%io_source,parai%cp_grp)
    ! ATCHG
    CALL mp_bcast(atchg,SIZE(atchg),parai%io_source,parai%cp_grp)
    ! CLASSICAL ATOMS
    IF (tclas) THEN
       CALL mp_bcast_byte(clas3, size_in_bytes_of(clas3),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cltyp,SIZE(cltyp),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(clas4, size_in_bytes_of(clas4),parai%io_source,parai%cp_grp)
       CALL mp_bcast(clasc,SIZE(clasc),parai%io_source,parai%cp_grp)
    ENDIF
    ! CONFINEMENT POTENTIAL
    CALL mp_bcast_byte(patom1, size_in_bytes_of(patom1),parai%io_source,parai%cp_grp)
    IF (patom1%pconf) THEN
       CALL mp_bcast_byte(patom2,size_in_bytes_of(patom2),parai%io_source,parai%cp_grp)
    ENDIF
    ! META DYNAMICS OF ORDER PARAMETERS
    CALL mp_bcast_byte(lmeta,size_in_bytes_of(lmeta),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(imeta,size_in_bytes_of(imeta),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(rmeta,size_in_bytes_of(rmeta),parai%io_source,parai%cp_grp)
    CALL mp_bcast(ncolvar,parai%io_source,parai%cp_grp)
    CALL mp_bcast(inter_hill,parai%io_source,parai%cp_grp)
    CALL mp_bcast(inter_hill_max,parai%io_source,parai%cp_grp)
    IF (lmeta%tlocalizespin) THEN
       IF (.NOT. paral%parent) THEN
          ALLOCATE(icv_spin(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(icv_spin)!,ncolvar)
          ALLOCATE(rcc0(3,ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(rcc0)!,3*ncolvar)
          ALLOCATE(kharm(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(kharm)!,ncolvar)
       ENDIF

       CALL mp_bcast(icv_spin,SIZE(icv_spin),parai%io_source,parai%cp_grp)
       CALL mp_bcast(rcc0,SIZE(rcc0),parai%io_source,parai%cp_grp)
       CALL mp_bcast(kharm,SIZE(kharm),parai%io_source,parai%cp_grp)
    ENDIF
    ! MULTIPLE WALKER METADYNAMICS
    CALL mp_bcast(tmw,parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(mwi, size_in_bytes_of(mwi),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm)THEN
       maxsys%nsx=mmdim%nspq
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ratom
  ! ==================================================================
  SUBROUTINE lread(iunit,lmax,loc,iskip)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit, lmax, loc, iskip

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i1, i2, iin, iout
    LOGICAL                                  :: erread

    IF (paral%io_parent)&
         READ(iunit,err=20,END=20,fmt='(A)') line
    IF (INDEX(line,'LMAX').NE.0) THEN
       i1=INDEX(line,'LMAX')
       i2=i1+INDEX(line(i1:80),'=')
       lmax=lval(line(i2:i2))
       IF (INDEX(line,'LOC').NE.0) THEN
          i1=INDEX(line,'LOC')
          i2=i1+INDEX(line(i1:80),'=')
          loc=lval(line(i2:i2))
       ELSE
          loc=-1
       ENDIF
       IF (INDEX(line,'SKIP').NE.0) THEN
          i1=INDEX(line,'SKIP')
          i2=i1+INDEX(line(i1:80),'=')
          iskip=lval(line(i2:i2))
       ELSE
          iskip=-1
       ENDIF
    ELSE
       CALL readsi(line,1,iout,lmax,erread)
       iin=iout
       CALL readsi(line,iin,iout,loc,erread)
       iin=iout
       CALL readsi(line,iin,iout,iskip,erread)
    ENDIF
    GOTO 30
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' RATOM: ERROR WHILE READING ON UNIT ',iunit
    CALL stopgm('RATOM',' ',& 
         __LINE__,__FILE__)
30  CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lread
  ! ==================================================================
  FUNCTION lval(l)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: l
    INTEGER                                  :: lval

! ==--------------------------------------------------------------==

    lval=-1
    IF (l.EQ.'S') lval=0
    IF (l.EQ.'P') lval=1
    IF (l.EQ.'D') lval=2
    IF (l.EQ.'F') lval=3
    IF (l.EQ.'G') lval=4
    IF (l.EQ.'H') lval=5
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION lval
  ! ==================================================================
  SUBROUTINE rdumat(line,iat)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: iat

    INTEGER                                  :: i1, iout, k
    LOGICAL                                  :: erread

    IF (INDEX(line,'TYPE1').NE.0) THEN
       i1=INDEX(line,'TYPE1')+5
       duat%ndat1=duat%ndat1+1
       duat%listda(iat,1)=1
       duat%listda(iat,2)=duat%ndat1
       DO k=1,3
          CALL readsr(line,i1,iout,duat%dummy1(k,duat%ndat1),erread)
          i1=iout
       ENDDO
    ELSEIF (INDEX(line,'TYPE2').NE.0) THEN
       i1=INDEX(line,'TYPE2')+5
       duat%ndat2=duat%ndat2+1
       duat%listda(iat,1)=2
       duat%listda(iat,2)=duat%ndat2
       CALL readsi(line,i1,iout,duat%listd2(1,duat%ndat2),erread)
       DO k=1,duat%listd2(1,duat%ndat2)
          i1=iout
          CALL readsi(line,i1,iout,duat%listd2(k+1,duat%ndat2),erread)
       ENDDO
    ELSEIF (INDEX(line,'TYPE3').NE.0) THEN
       i1=INDEX(line,'TYPE3')+5
       duat%ndat3=duat%ndat3+1
       duat%listda(iat,1)=3
       duat%listda(iat,2)=duat%ndat3
       CALL readsi(line,i1,iout,duat%listd3(1,duat%ndat3),erread)
       DO k=1,duat%listd3(1,duat%ndat3)
          i1=iout
          CALL readsi(line,i1,iout,duat%listd3(k+1,duat%ndat3),erread)
       ENDDO
    ELSEIF (INDEX(line,'TYPE4').NE.0) THEN
       i1=INDEX(line,'TYPE4')+5
       duat%ndat4=duat%ndat4+1
       duat%listda(iat,1)=4
       duat%listda(iat,2)=duat%ndat4
       CALL readsi(line,i1,iout,duat%listd4(1,duat%ndat4),erread)
       DO k=1,duat%listd4(1,duat%ndat4)
          i1=iout
          CALL readsi(line,i1,iout,duat%listd4(k+1,duat%ndat4),erread)
       ENDDO
       DO k=1,duat%listd4(1,duat%ndat4)
          i1=iout
          CALL readsr(line,i1,iout,duat%weigd4(k,duat%ndat4),erread)
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' ERROR IN DUMMY ATOM INPUT LINE ',iat
       CALL stopgm('RATOM',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rdumat
  ! ==================================================================

END MODULE ratom_utils
