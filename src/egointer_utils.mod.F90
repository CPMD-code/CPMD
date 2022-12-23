MODULE egointer_utils
  USE atomc_utils,                     ONLY: atomc
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp
  USE cnst,                            ONLY: au_deb,&
                                             fbohr
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cppt,                            ONLY: indz,&
                                             nzh
  USE difrho_utils,                    ONLY: difrho
  USE dipo_utils,                      ONLY: dipo
  USE dipomod,                         ONLY: moment
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf,&
                                             textfld
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE elstpo_utils,                    ONLY: elstpo
  USE ener,                            ONLY: ener_com
  USE epot_types,                      ONLY: &
       epot1, epot2, epot3, epot4, epot5, maxind, maxintm, maxnear, maxoutr
  USE error_handling,                  ONLY: stopgm
  USE espchg_utils,                    ONLY: atfield,&
                                             espsolv,&
                                             selectp
  USE exdipo_utils,                    ONLY: exdipo
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fileopenmod,                     ONLY: fo_info
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE hip_utils,                       ONLY: give_qphi
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lodipo_utils,                    ONLY: lodipo
  USE lodp,                            ONLY: dmomlo,&
                                             extd,&
                                             focc
  USE machine,                         ONLY: m_getarg,&
                                             m_iargc,&
                                             m_sleep,&
                                             m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE mulliken_utils,                  ONLY: give_scr_mulliken,&
                                             mulliken
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE propin_utils,                    ONLY: propin
  USE proppt_utils,                    ONLY: give_scr_espc,&
                                             prdip
  USE pslo,                            ONLY: pslo_com
  USE purge_utils,                     ONLY: purge
  USE readsr_utils,                    ONLY: xstring
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rrane_utils,                     ONLY: rrane
  USE setirec_utils,                   ONLY: write_irec
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: cprint,&
                                             restart1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, ncpw, nkpt, parap, parm, spar
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrener
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INTERFACE

CONTAINS

  ! ***************************************************************************
  ! ***************************************************************************
  ! .....Interface to the cntl%md simulation program EGO.
  ! .....by M. Eichinger and J. Hutter
  ! .....extended by P.K.Biswas to support Gromacs.
  ! .....
  ! .....Functions in this file:
  ! .....INTERFACE(...)    The main loop for the interface mode.
  ! .....INTERFACE_WAIT(.) Waits for the QMCONTINUE file created by EGO.
  ! .....INTERFACE_READ(.) Reads the new data from EGO.
  ! .....RDEGOF(..)        Reads the external field.
  ! .....EXT_FIELD         Calculates the external field at grid points.
  ! .....CLCPOT(..)        Calculates the electrostatic field at any point.
  ! .....CLCFLD(..)        Calculates the same for a s-wave expanded point charge.
  ! .....INTERFACE_WRITE() Writes the forces on the QM atoms for EGO.
  ! .....EGOSTOP           Creates a file for stopping EGO in case of problems.
  ! .....MYPROPPT(...)     Calculates the dipole moment
  ! .....GMX_EL_FORCE(...) Calculates QM-2-MM electrostatic forces on MM atoms.
  ! ***************************************************************************
  ! ***************************************************************************
  ! ==================================================================
  SUBROUTINE INTERFACE(c0,c2,sc0,pme,gde,vpp,eigv)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,crge%n), &
                                                c2(nkpt%ngwk,crge%n), &
                                                sc0(nkpt%ngwk,crge%n), &
                                                pme(*), gde(*)
    REAL(real_8)                             :: vpp(*), eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'INTERFACE'

    CHARACTER(len=80)                        :: int_filen, line
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER                                  :: i, ierr, il_psi_1d, &
                                                il_psi_2d, il_rhoe_1d, &
                                                il_rhoe_2d, irec(100)
    LOGICAL                                  :: eofrun, fnowf
    REAL(real_8)                             :: detot, dummy(1), ekin1, &
                                                ekin2, ekincp, ekinh1, &
                                                ekinh2, etoto, tcpu, temp1, &
                                                temp2, time1, time2
    REAL(real_8), ALLOCATABLE                :: rhoe(:,:)

    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d,&
         il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CLOSE(5)
    time1 =m_walltime()
    ! ..+1 for esp charges
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(velp(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    CALL zeroing(velp)!,SIZE(velp))
    ! ..added: biswas
    ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(extf)!,kr1*kr2s*kr3s)
    ! ..addition ends
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(epot2%boxqm)!,6)
    CALL zeroing(epot2%boxdum)!,6)
    CALL zeroing(epot2%koloc)!,3*maxind)
    CALL zeroing(epot2%locexp)!,10*maxind)
    CALL zeroing(epot2%datnear)!,5*maxnear)
    CALL zeroing(epot2%konear)!,3*maxnear)
    CALL zeroing(epot1%locid)!,maxind)
    CALL zeroing(epot1%nincl)!,maxind)
    CALL zeroing(epot1%idnear)!,maxnear)
    CALL zeroing(epot1%lincl)!,maxind*maxnear)
    DO i=1,100
       epot2%myrag(i)=1.2_real_8
    ENDDO
    epot2%myrag(1)=0.85_real_8
    epot2%myrag(8)=1.23_real_8

    epot1%nnear = 0
    epot1%nloc  = 0

    restart1%restart=.FALSE.
    restart1%rwf    =.TRUE.
    eofrun=.FALSE.
    iteropt%nfi   = 0
    ener_com%ecnstr= 0.0_real_8
    ener_com%erestr= 0.0_real_8
    etoto = 0.0_real_8

    IF ((cnti%iftype.EQ.1).AND.paral%io_parent)&
         WRITE(6,'(A)') ' INTERFACE| ENTERING EGO-INTERFACE'
    IF ((cnti%iftype.EQ.2).AND.paral%io_parent)&
         WRITE(6,'(A)') ' INTERFACE| ENTERING GROMACS-INTERFACE'

    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! ==--------------------------------------------------------------==
    ! ==      THE BASIC MOLECULAR DYNAMICS LOOP                       ==
    ! ==--------------------------------------------------------------==
1   CONTINUE
    iteropt%nfi=iteropt%nfi+1
    ropt_mod%convwf=.FALSE.
    ropt_mod%sdiis =.TRUE.
    ropt_mod%spcg  =.TRUE.
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! 
    ! .......wait for the interface file
    ! 
    CALL interface_wait(int_filen)
    IF (.NOT.eofrun) CALL testex(eofrun)
    IF (eofrun) GOTO 2
    ! 
    ! .......read and process interface file
    ! 
    CALL interface_read(tau0)
    CALL phfac(tau0)
    IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)
    CALL mp_sync(parai%allgrp)
    ! 
    ! .......Randomize Wavefunction
    ! 
    IF (cntl%trane) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.8)')&
            ' RANDOMIZE INITIAL WAVEFUNCTION,  AMPLITUDE=',cntr%ampre
       CALL rrane(c0,c2,crge%n)
    ENDIF
    ! 
    fnowf=.FALSE.
    DO infi=1,cnti%nomore
       time1=m_walltime()
       IF (paral%parent) THEN
          ropt_mod%engpri=MOD(infi-1,cprint%iprint_step).EQ.0
       ELSE
          ropt_mod%engpri=.FALSE.
       ENDIF
       ! 
       ! .......update the wavefunctions 
       ! 
       CALL updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,&
            eigv,rhoe,psi,crge%n,fnowf,.TRUE.)
       ! 
       ! .......printout evolution of the progress
       ! 
       IF (paral%parent) THEN
          detot=ener_com%etot+ener_com%eext-etoto
          IF (infi.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF (ropt_mod%engpri) CALL wrener
          IF ((ropt_mod%engpri).AND.paral%io_parent)&
               WRITE(6,'(A)')&
               ' NFI     GEMAX     CNORM          ETOT         DETOT        TCPU'
          IF (paral%io_parent)&
               WRITE(6,'(I4,F10.6,F10.6,F14.5,F14.6,F12.2) ')&
               infi,gemax,cnorm,ener_com%etot+ener_com%eext,detot,tcpu
          etoto=ener_com%etot+ener_com%eext
       ENDIF
       IF (ropt_mod%convwf.AND.fnowf) GOTO 100
       IF (ropt_mod%convwf) THEN
          ropt_mod%convwf=.FALSE.
          fnowf=.TRUE.
       ENDIF
    ENDDO
100 CONTINUE
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)'EXTERNAL ENERGY    =',ener_com%eext,' AU'
       IF (paral%io_parent)&
            WRITE(6,*)'REAL TOTAL ENERGY  =',ener_com%eext+ener_com%etot,' AU'
    ENDIF
    ! ..
    ! ..Added by biswas: dec/2004.
    ! ..   For calculating electrostatic forces on Gromacs (MM) atoms
    IF (cnti%iftype.EQ.2) THEN
       CALL gmx_el_forces(c0,tau0,rhoe,psi,crge%n,1)
       ! CALL ELF(C0,TAU0,RHOE,PSI,SCR,LSCR,N,1)
    ENDIF
    ! ..addition ends here.
    ! ..
    ! ..Calculate gradient
    ! 
    CALL forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.TRUE.,.TRUE.)
    ! 
    ! 
    ! ..Calculate the dipole moment
    ! 
    IF (cntl%proper) THEN
       CALL m_getarg(1,line)
       IF (paral%io_parent) OPEN(unit=5,file=line,status='OLD',err=2)
       CALL myproppt(c0,tau0,rhoe,psi)
       IF (paral%io_parent) CLOSE(5)
    ENDIF
    ! 
    ! ..impose constraints on the gradient
    ! 
    CALL purge(tau0,fion)
    ! 
    ! ..process data and write the interface file
    ! 
    CALL interface_write(int_filen,tau0,fion,c0,taup,rhoe,psi)
    ! 
    ! ..write the restart file
    ! 
    CALL write_irec(irec)
    CALL zhwwf(2,irec,c0,c2,crge%n,dummy,tau0,velp,taup,iteropt%nfi)
    ! 
    ! ..test for exit-file
    ! 
    IF (.NOT.eofrun) CALL testex(eofrun)
    IF (eofrun) GOTO 2
    GOTO 1
2   CONTINUE
    ! 
    CALL write_irec(irec)
    CALL zhwwf(2,irec,c0,c2,crge%n,dummy,tau0,velp,taup,iteropt%nfi)
    IF (paral%parent) CALL finalp(tau0,fion,velp,dummy)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ..added: biswas.
    DEALLOCATE(extf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE INTERFACE
  ! ==================================================================
  ! 
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE interface_read(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ierr, imdum, isp, iunit, &
                                                nait

    IF (paral%parent) THEN
       iunit=5
       imdum=m_iargc()
       IF (imdum.LT.1) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' INTERFACE| NO INPUT FILE NAME SPECIFIED '
          CALL stopgm('INTERFACE',' ',& 
               __LINE__,__FILE__)
       ENDIF
       CALL m_getarg(1,line)
       IF (paral%io_parent)&
            WRITE(6,*) ' INTERFACE| READING NEW COORDINATES FROM FILE ',line
       IF (paral%io_parent)&
            OPEN(unit=iunit,file=line,status='OLD',err=300)
       ! 
       ! .......read atomic coordinates
       ! 
       ierr=inscan(iunit,'&ATOM')
       isp=0
10     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,END=200,err=300,fmt='(A)') line
       IF (INDEX(line,'&END').NE.0) GOTO 30
       IF (line(1:1).NE.'*') GOTO 10
       isp=isp+1
       IF (paral%io_parent)&
            READ(iunit,*)
       IF (paral%io_parent)&
            READ(iunit,*) nait
       DO i=1,nait
          IF (paral%io_parent)&
               READ(iunit,*) tau0(1,i,isp),tau0(2,i,isp),tau0(3,i,isp)
       ENDDO
       GOTO 10
30     CONTINUE
       ! 
       ! .......read external field
       ! 
       ierr=inscan(iunit,'&EXTE')
       IF (ierr.EQ.0) THEN
          CALL rdegof(iunit)
          textfld=.TRUE.
       ELSE
          textfld=.FALSE.
       ENDIF
       IF (paral%io_parent)&
            CLOSE(iunit)
       ! 
       ! .......terminate reading normal
       ! 
200    CONTINUE
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERFACE| READING NEW COORDINATES ... DONE'
       GOTO 400
300    CONTINUE
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERFACE| CAN NOT READ FROM FILE ',line
       CALL stopgm('INTERFACE',' ',& 
            __LINE__,__FILE__)
400    CONTINUE
    ENDIF
    ! .....broadcast information to all processors
    CALL mp_sync(parai%allgrp)
    ! .....tau0
    CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
    ! .....textfld
    CALL mp_bcast(textfld,parai%source,parai%allgrp)
    CALL mp_sync(parai%allgrp)
    IF (textfld) THEN
       ! .....epot1
       CALL mp_bcast_byte(epot1, size_in_bytes_of(epot1),parai%source,parai%allgrp)
       ! .....epot2
       CALL mp_bcast_byte(epot2, size_in_bytes_of(epot2),parai%source,parai%allgrp)
       ! ..added: biswas
       IF (cnti%iftype.EQ.2) THEN
          ! .....epot3
          CALL mp_bcast_byte(epot3, size_in_bytes_of(epot3),parai%source,parai%allgrp)
          ! .....epot4
          CALL mp_bcast_byte(epot4, size_in_bytes_of(epot4),parai%source,parai%allgrp)
          ! .....epot5
          CALL mp_bcast_byte(epot5, size_in_bytes_of(epot5),parai%source,parai%allgrp)
       ENDIF
       ! ..addition ends
    ENDIF
    ! .....calculate external field at grid points
    IF (textfld) CALL ext_field
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interface_read
  ! ==================================================================
  ! 
  ! 
  ! 
  ! ##############################################
  SUBROUTINE rdegof(iunit)
    ! ##############################################
    INTEGER                                  :: iunit

    CHARACTER(len=100)                       :: string
    CHARACTER(len=7)                         :: substring
    INTEGER                                  :: i, ilines, j, jj, k

! 
! ***** Read bounding box **********************************

    IF (paral%io_parent)&
         READ(iunit,'(A)') string
    IF (paral%io_parent)&
         READ(iunit,*) epot2%boxdum(1),epot2%boxdum(2),&
         epot2%boxdum(3),epot2%boxdum(4),&
         epot2%boxdum(5),epot2%boxdum(6)
    PRINT'(A10,F11.5,F11.5,F11.5,F11.5,F11.5,F11.5)',string,&
         epot2%boxdum(1),epot2%boxdum(2),&
         epot2%boxdum(3),epot2%boxdum(4),&
         epot2%boxdum(5),epot2%boxdum(6)
    ! 
    ! ***** Read QM-bounding box ************************************
    IF (paral%io_parent)&
         READ(iunit,'(A)') string
    IF (paral%io_parent)&
         READ(iunit,*) epot2%boxqm(1),epot2%boxqm(2),&
         epot2%boxqm(3),epot2%boxqm(4),&
         epot2%boxqm(5),epot2%boxqm(6)
    PRINT'(A10,F11.5,F11.5,F11.5,F11.5,F11.5,F11.5)',string,&
         epot2%boxqm(1),epot2%boxqm(2),&
         epot2%boxqm(3),epot2%boxqm(4),&
         epot2%boxqm(5),epot2%boxqm(6)
    ! 
    ! ..added: biswas // to implement multiply MM layer; 
    IF (cnti%iftype.EQ.2) THEN
       ! ***** Read Step-No *****************************************
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%stpno
       PRINT'(A15,I5)',string, epot3%stpno
       ! 
       ! ***** Read LMAX for partial wave expansion of MM-charge*******
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%pwlmax
       PRINT'(A15,I5)',string, epot3%pwlmax
       ! 
       ! ***** Read Intm-Freq *****************************************
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%intmf
       PRINT'(A15,I5)',string, epot3%intmf
       ! 
       ! ***** Read Outm-Freq *****************************************
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%outmf
       PRINT'(A15,I5)',string, epot3%outmf
       ! 
    ELSEIF (cnti%iftype.EQ.1) THEN
       ! .....Read local expansions ************************************
       ! 
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot1%nloc
       PRINT'(A11,I8)',string,epot1%nloc
       IF (epot1%nloc.GE.maxind) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' INTERFACE| Error: MAXIND not big enough! '
          IF (paral%io_parent)&
               WRITE(6,*) '            Increase this value in epot.inc'
          CALL stopgm('INTERFACE','To many QM-atoms.',& 
               __LINE__,__FILE__)
       ENDIF
       ! .....
       ! .....Loop over all local expansions
       ! 
       DO 30 i=1,epot1%nloc
          ! ...........Read position of local expansion
          IF (paral%io_parent)&
               READ(iunit,*) epot1%locid(i),epot2%koloc(1,i),epot2%koloc(2,i),epot2%koloc(3,i)
          IF (paral%io_parent)&
               READ(iunit,*) epot2%locexp(1,i)
          IF (paral%io_parent)&
               READ(iunit,*) epot2%locexp(2,i),epot2%locexp(3,i),epot2%locexp(4,i)
          IF (paral%io_parent)&
               READ(iunit,*) epot2%locexp(5,i),epot2%locexp(6,i),epot2%locexp(7,i),&
               epot2%locexp(8,I),epot2%locexp(9,I),epot2%locexp(10,I)
          IF (paral%io_parent)&
               READ(iunit,*) epot1%nincl(i)
          ! .....
          IF (i.LE.3) THEN
             PRINT'(I8,TR2,F11.5,F11.5,F11.5,TR2,A,I3)',epot1%locid(i),&
                  epot2%koloc(1,I),epot2%koloc(2,I),epot2%koloc(3,I),'INCL: ',epot1%nincl(I)
          ENDIF
          IF (epot1%nincl(i).GE.maxnear) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' INTERFACE| Error: MAXNEAR not big enough! '
             IF (paral%io_parent)&
                  WRITE(6,*) '            Increase this value in epot.inc'
             CALL stopgm('INTERFACE','To many near atoms.',& 
                  __LINE__,__FILE__)
          ENDIF
          ! 
          ! ...........Read the near inclusion list *********
          ! 
          ilines = INT(epot1%nincl(i)/10) + 1
          DO 35 k=0,ilines-1
             IF (paral%io_parent)&
                  READ(iunit,'(A)') string
             DO 40 j=0,9
                jj = 10*k + j
                IF (jj.GE.epot1%nincl(i)) GOTO 30
                substring = string(1+j*7:j*7+7)
                IF (paral%io_parent)&
                     READ(substring,*) epot1%lincl(i,jj+1)
                epot1%lincl(i,jj+1) = epot1%lincl(i,jj+1) + 1
40           CONTINUE
35        CONTINUE
30     CONTINUE
    ENDIF
    ! 
    ! .....Read list of near atoms ***************************
    ! 
    IF (paral%io_parent)&
         READ(iunit,'(A)') string
    IF (paral%io_parent)&
         READ(iunit,*) epot1%nnear
    PRINT'(A11,I5)',string,epot1%nnear
    IF (epot1%nnear.GE.maxnear) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' INTERFACE| Error: MAXNEAR not big enough! '
       IF (paral%io_parent)&
            WRITE(6,*) '            Increase this value in epot.inc'
       CALL stopgm('INTERFACE','To many near atoms.',& 
            __LINE__,__FILE__)
    ENDIF
    DO 50 i=1,epot1%nnear
       IF (paral%io_parent)&
            READ(iunit,*) epot1%idnear(i),&
            epot2%datnear(1,I),      & ! Partial Charge
            epot2%datnear(4,I),      & ! CutOffWidth
            epot2%datnear(5,I),      & ! CutOffValue
            epot2%datnear(2,I),      & ! Core charge
            epot2%datnear(3,I),      & ! Sigma
            epot2%konear(1,I),epot2%konear(2,I),epot2%konear(3,I) ! the position of the atoms
       ! ..modified: biswas/2005.
       epot4%zmmnear(I) = epot2%datnear(2,I) ! !nuclear charges are stored: biswas.
       ! ..
       IF (cnti%iftype.EQ.1) THEN
          IF (epot2%datnear(3,i).LE.0.0_real_8) THEN
             epot2%datnear(2,i) =0.0_real_8
             epot2%datnear(3,i) =-epot2%datnear(3,i)
          ENDIF
       ENDIF
       ! ..modification ends here.
       IF (i.LE.3) THEN
          PRINT'(I8,TR2,F11.5,F11.5,F11.5,F11.5,F11.5,F11.5)',&
               epot1%idnear(I),&
               epot2%datnear(1,I),&
               epot2%datnear(2,I),epot2%datnear(3,I),&
               epot2%konear(1,I),epot2%konear(2,I),epot2%konear(3,I)            
       ENDIF
50  CONTINUE
    ! .....Read list of intermediate layer MM atoms *************************C
    ! 
    IF (cnti%iftype.EQ.2) THEN
       ! 
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%nintm
       PRINT'(A11,I5)',string,epot3%nintm
       IF (epot3%nintm.GE.maxintm) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' INTERFACE| Error: MAXINTM not big enough! '
          IF (paral%io_parent)&
               WRITE(6,*) '            Increase this value in epot.inc'
          CALL stopgm('INTERFACE','To many intermediate MM atoms.',& 
               __LINE__,__FILE__)
       ENDIF
       IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
          DO 51 i=1,epot3%nintm
             IF (paral%io_parent)&
                  READ(iunit,*) epot3%idintm(i),&
                  epot4%qintm(1,I), & ! ......Partial Charge
                  epot4%qintm(2,I), & ! ......Sigma
                  epot4%kointm(1,I),epot4%kointm(2,I),epot4%kointm(3,I) ! position of the atoms
             ! ..   
51           CONTINUE
       ENDIF
       ! 
       ! .....Read list of outermost layer MM atoms *************************C
       ! 
       IF (paral%io_parent)&
            READ(iunit,'(A)') string
       IF (paral%io_parent)&
            READ(iunit,*) epot3%noutr
       PRINT'(A11,I5)',string,epot3%noutr
       IF (epot3%noutr.GE.maxoutr) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' INTERFACE| Error: MAXOUTR not big enough! '
          IF (paral%io_parent)&
               WRITE(6,*) '            Increase this value in epot.inc'
          CALL stopgm('INTERFACE','To many outermost MM atoms.',& 
               __LINE__,__FILE__)
       ENDIF
       ! ..
       IF (MOD(epot3%stpno,epot3%outmf).EQ.0.AND.epot3%noutr.GT.0) THEN
          DO 52 i=1,epot3%noutr
             IF (paral%io_parent)&
                  READ(iunit,*) epot3%idoutr(i),&
                  epot4%qoutr(1,I), & ! ..............Partial Charge
                  epot4%qoutr(2,I), & ! ...........Sigma, there is no QOUTR(2,I); 
                  epot4%kooutr(1,I),epot4%kooutr(2,I),epot4%kooutr(3,I)
             ! addl. dimension is to facilitate calc.
             ! ....................the position of the atoms
             ! ..
52           CONTINUE
       ENDIF
       ! ..
    ENDIF
    ! 
    ! *****
  END  SUBROUTINE RDEGOF
  ! *****END of SUBROUTINE RDEGOF(FNAME)
  ! 
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE ext_field
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Variables
    REAL(real_8) :: k0(0:2),di,dj,dk
    INTEGER :: isub,i,j,k,ikr1,ir
    ! ==--------------------------------------------------------------==
    CALL tiset(' EXT_FIELD',isub)
    ! .. memory for EXTF already allocated in INTERFACE(): biswas
    IKR1 = MOD(parm%nr1+1,2)   ! !added: biswas

    dk=(epot2%boxdum(6)-epot2%boxdum(5))/REAL(spar%nr3s,kind=real_8)
    dj=(epot2%boxdum(4)-epot2%boxdum(3))/REAL(spar%nr2s,kind=real_8)
    di=(epot2%boxdum(2)-epot2%boxdum(1))/REAL(spar%nr1s,kind=real_8)

    IF (cnti%iftype.EQ.1) THEN
       DO k=1,fpar%kr3s
          k0(2)=epot2%boxdum(5)+REAL(k-1,kind=real_8)*dk
          DO j=1,fpar%kr2s
             k0(1)=epot2%boxdum(3)+REAL(j-1,kind=real_8)*dj
             DO i=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)+ikr1! IKR1 added: biswas
                k0(0)=epot2%boxdum(1)+REAL(i-1,kind=real_8)*di
                extf((i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s) &
                     =clcpot(k0)
             ENDDO
          ENDDO
       ENDDO
       ! ..
    ELSEIF (cnti%iftype.EQ.2) THEN
       ! ..
       ir = 0
       DO k=1,fpar%kr3s
          k0(2)=epot2%boxdum(5)+REAL(k-1,kind=real_8)*dk
          DO j=1,fpar%kr2s
             k0(1)=epot2%boxdum(3)+REAL(j-1,kind=real_8)*dj
             DO i=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)+ikr1! IKR1 added: biswas
                ir = ir +1
                k0(0)=epot2%boxdum(1)+REAL(i-1,kind=real_8)*di
                extf(i-parap%nrxpl(parai%mepos,1)+1 + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s) &
                     =clcfld(ir,k0)
             ENDDO
          ENDDO
       ENDDO
       ! ..   
    ENDIF
    ! ..
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,'(A)') ' INTERFACE| EXTERNAL FIELD CALCULATED '
    ENDIF
    CALL tihalt(' EXT_FIELD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ext_field
  ! ==================================================================
  ! 
  ! 
  ! 
  ! ##############################################
  FUNCTION clcpot(ko)
    ! ##############################################
    REAL(real_8)                             :: ko(0:2), clcpot

    INTEGER                                  :: i, it, li
    REAL(real_8)                             :: addepot, d, dmin, dx, dy, dz

! Variables
! 
! .... find the nearest local expansion ********
! 

    clcpot = 0._real_8
    dmin   = 999999._real_8
    li     = -1
    ! 
    IF (epot1%nloc .EQ. 0) GOTO 150
    ! 
    DO 100 i=1,epot1%nloc
       dx=ko(0)-epot2%koloc(1,i)
       dy=ko(1)-epot2%koloc(2,i)
       dz=ko(2)-epot2%koloc(3,i)
       d=dx*dx+dy*dy+dz*dz
       IF (d .LT. dmin) THEN
          dmin=d
          li=i
       ENDIF
100    CONTINUE
    ! 
    ! .... evaluate potential from local expansion *****
    ! 
    dx=ko(0)-epot2%koloc(1,li)
    dy=ko(1)-epot2%koloc(2,li)
    dz=ko(2)-epot2%koloc(3,li)
    ! 
    clcpot =   -( epot2%locexp(1,li) +&
         DX*epot2%locexp(2,LI) + &
         DY*epot2%locexp(3,LI) + &
         DZ*epot2%locexp(4,LI) +&
         0.5_real_8*DX*DX*epot2%locexp(5,LI) +&
         0.5_real_8*DY*DY*epot2%locexp(6,LI) +&
         0.5_real_8*DZ*DZ*epot2%locexp(7,LI) +&
         DX*DY*epot2%locexp(8,LI) +&
         DX*DZ*epot2%locexp(9,LI) +&
         DY*DZ*epot2%locexp(10,LI) )
    ! 
150 CONTINUE
    ! 
    IF (epot1%nincl(li) .EQ. 0) GOTO 300
    ! 
    ! .....add El.Pot from all near atoms      
    ! 
    DO 200 i=1,epot1%nincl(li)
       it = epot1%lincl(li,i)
       ! 
       dx=ko(0)-epot2%konear(1,it)
       dy=ko(1)-epot2%konear(2,it)
       dz=ko(2)-epot2%konear(3,it)

       ! ............Distance in cntl%bohr. 
       d = SQRT(dx*dx + dy*dy + dz*dz)
       IF (d.LT.0.1_real_8) d=0.1_real_8
       ! 
       IF (d .GT. 6.0_real_8) THEN
          ! If we a far away take simple Coulumb-Potential
          ! Phi = q / d
          addepot = epot2%datnear(1,it) / d
       ELSE IF (epot2%datnear(2,it) .GT. 0.5_real_8) THEN
          ! If core-charges are provided use the more
          ! advanced atom-potential.
          dy=epot2%myrag(NINT(epot2%datnear(2,it)))
          dx = (epot2%datnear(2,it) - epot2%datnear(1,it))*&
               (1.0_real_8-EXP(-d/epot2%datnear(3,it)))
          addepot = cp_erf(d/dy)*(epot2%datnear(2,it)-dx)/d
       ELSE
          ! Use the simple atom-potential given 
          ! just by partial charges.            
          ! .              Phi = q * erf(d/sigma) / d
          dx = cp_erf(d/epot2%datnear(3,it))
          addepot = epot2%datnear(1,it) * dx / d
       ENDIF

       ! ............we have an excluded atom here
       IF (epot2%datnear(4,it) .GT. 0.01_real_8) THEN
          IF (d .GT. (epot2%datnear(5,it)+epot2%datnear(4,it))) THEN
             ! do nothing in that case
             ! because we are far away from it 
          ELSE IF (d .LT. epot2%datnear(5,it)) THEN
             ! we are to close, dont add any contribution from it
             addepot = 0.0_real_8
          ELSE
             ! Use a switching function: No typo here!
             dx = (d - epot2%datnear(5,it)) / epot2%datnear(4,it)
             dx = 1.0_real_8 - dx*dx
             addepot = (1.0_real_8 - dx*dx) * addepot
          ENDIF
       ENDIF

       clcpot = clcpot - addepot
200    CONTINUE
    IF (clcpot .LT. -0.3_real_8) clcpot = -0.3_real_8 - LOG(0.7_real_8-clcpot)
    ! 
300 CONTINUE
    ! 
  END FUNCTION clcpot

  ! 
  ! ..
  ! ..added by biswas for GMX interface.
  ! ##############################################
  FUNCTION clcfld(ir,ko)
    ! ##############################################
    INTEGER                                  :: ir
    REAL(real_8)                             :: ko(0:2), clcfld

    CHARACTER(*), PARAMETER                  :: procedureN = 'clcfld'

    INTEGER                                  :: i
    REAL(real_8)                             :: addepot, alam, alfa, beta, d, &
                                                dx, dy, dz, potintm, potoutr

! Variables
! 

    clcfld = 0._real_8
    alam = 1.3_real_8
    ! 
    ! ..!! For atoms in the innermost layer
    IF (epot1%nnear.GT.0) THEN
       DO 200 i=1,epot1%nnear
          ! 
          dx=ko(0)-epot2%konear(1,i)
          dy=ko(1)-epot2%konear(2,i)
          dz=ko(2)-epot2%konear(3,i)

          ! ............Distance in cntl%bohr
          d = SQRT(dx*dx + dy*dy + dz*dz)
          ! ..             
          alfa = alam*epot2%datnear(3,i)
          beta = 2._real_8*alfa
          ! ............biswas-modified Coulomb
          addepot = epot2%datnear(1,i)*&
               (1._real_8/d - EXP(-beta*d)*(1._real_8/d + alfa))
          ! 
          clcfld = clcfld - addepot
          ! 
200       CONTINUE
    ENDIF
    ! 
    ! ..   !! For atoms in the intermediate layer
    IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
       potintm = 0._real_8
       DO 300 i=1,epot3%nintm
          ! 
          dx=ko(0)-epot4%kointm(1,i)
          dy=ko(1)-epot4%kointm(2,i)
          dz=ko(2)-epot4%kointm(3,i)
          ! ..   
          d = SQRT(dx*dx + dy*dy + dz*dz)
          ! ..   
          alfa = alam*epot4%qintm(3,i)
          beta = 2._real_8*alfa
          ! .............simple Coulomb
          ! POTINTM = POTINTM + QINTM(I)/d  !! 
          ! ............biswas-modified Coulomb
          potintm = potintm + epot4%qintm(1,i)*&
               (1._real_8/d - EXP(-beta*d)*(1._real_8/d + alfa))
          ! 
300       CONTINUE
       ! 
          !vw commented next line (not allocated) and added stop
          !vw extfl1(ir) = potintm
          CALL stopgm(procedureN,'array not allocated', &
               __LINE__,__FILE__)

    ENDIF
    ! 
    ! ..   !! For atoms in the outermost layer
    IF (MOD(epot3%stpno,epot3%outmf).EQ.0.AND.epot3%noutr.GT.0) THEN
       potoutr = 0._real_8
       DO 400 i=1,epot3%noutr
          ! 
          dx=ko(0)-epot4%kooutr(1,i)
          dy=ko(1)-epot4%kooutr(2,i)
          dz=ko(2)-epot4%kooutr(3,i)
          ! ..   
          d = SQRT(dx*dx + dy*dy + dz*dz)
          ! ..   
          alfa = alam*epot4%qoutr(3,i)
          beta = 2._real_8*alfa
          ! .............simple Coulomb
          ! POTOUTR = POTOUTR + QOUTR(1,I)/d  !! 
          ! ............biswas-modified Coulomb
          potoutr = potoutr + epot4%qoutr(1,i)*&
               (1._real_8/d - EXP(-beta*d)*(1._real_8/d + alfa))
          ! 
400       CONTINUE
       !vw commented next line (not allocated) and added stop
       !vw extfl2(ir) = potoutr
          CALL stopgm(procedureN,'array not allocated', &
               __LINE__,__FILE__)
    ENDIF
    ! ..
    !vw commented next line (not allocated) and added stop
    !vw clcfld = clcfld - extfl1(ir) - extfl2(ir)
    CALL stopgm(procedureN,'array not allocated', &
         __LINE__,__FILE__)

    ! ..
    RETURN
  END FUNCTION clcfld
  ! 
  ! ==================================================================
  SUBROUTINE interface_wait(int_filen)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: int_filen

    CHARACTER(len=80)                        :: filen
    INTEGER                                  :: isub
    LOGICAL                                  :: test

    CALL tiset('INTER_WAIT',isub)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,63("="))')
       filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'QMCONTINUE'
       IF (paral%io_parent)&
            WRITE(6,*) ' INTERFACE| WAIT FOR CONTINUE-FILE ', filen
1      CONTINUE
       IF (paral%io_parent)&
            INQUIRE(file=filen,exist=test)
       IF (test) THEN
          int_filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'qmoutput.out'
          ! .......Now delete the QMCONTINUE file
          IF (paral%io_parent)&
               OPEN(5,file=filen,status='UNKNOWN')
          IF (paral%io_parent)&
               CLOSE(5,status='DELETE')
          IF (paral%io_parent)&
               WRITE(6,*) ' INTERFACE| CONTINUE QM CALCULATION ',int_filen
          GOTO 3
       ENDIF
       IF (paral%io_parent)&
            INQUIRE(file='EXIT',exist=test)
       IF (test) GOTO 3
       CALL m_sleep(1)
       GOTO 1
3      CONTINUE
    ENDIF
    CALL mp_sync(parai%allgrp)
    CALL tihalt('INTER_WAIT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interface_wait
  ! ==================================================================
  ! 
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE interface_write(int_filen,tau0,fion,c0,achrg,rhoe,v)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*)                         :: int_filen
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: achrg(ions1%nat), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v(maxfftn)

    CHARACTER(*), PARAMETER                  :: procedureN = 'interface_write'

    CHARACTER(len=80)                        :: line, str
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), qphi(:), &
                                                vtemp(:)
    INTEGER                                  :: i, ia, iat, idnr, ie, ierr, &
                                                ig, il_qphi, inr, ippc, is, &
                                                ish, l, mlen, nstate
    INTEGER, ALLOCATABLE                     :: isel(:)
    REAL(real_8)                             :: d1, d2, d3, dd, dx, dy, dz, &
                                                eef, ef(0:2), k0(0:2), k1(0:2)
    REAL(real_8), ALLOCATABLE                :: efield(:), rbuff(:), reirop(:)

    CALL setfftn(0)
    CALL zeroing(achrg)!,ions1%nat)
    ! 
    ! ..calculate atomic charges
    ! 
    IF (textfld.AND.cnti%iftype.EQ.1) THEN
       IF (cnti%icmet.EQ.2 .OR. cnti%icmet.EQ.3) THEN

          ! TODO check stat
          ! TODO align for BG 
          ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(reirop(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          CALL give_qphi(il_qphi)
          ALLOCATE(qphi(il_qphi),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

       ENDIF
       IF (cnti%icmet.EQ.1) THEN
          ! ..Mulliken Charges
          atwp%nattot=0
          DO is=1,ions1%nsp
             atwf_mod%numaor(is)=0
             DO ish=1,atwf_mod%nshell(is)
                l=atwf_mod%lshell(ish,is)
                atwp%nattot=atwp%nattot+ions0%na(is)*(2*l+1)
                atwf_mod%numaor(is)=atwf_mod%numaor(is)+(2*l+1)
             ENDDO
          ENDDO
          CALL mulliken(1,c0,nstate,tau0,achrg)
       ELSEIF (cnti%icmet.EQ.4) THEN
          ! ..Lowdin Charges
          atwp%nattot=0
          DO is=1,ions1%nsp
             atwf_mod%numaor(is)=0
             DO ish=1,atwf_mod%nshell(is)
                l=atwf_mod%lshell(ish,is)
                atwp%nattot=atwp%nattot+ions0%na(is)*(2*l+1)
                atwf_mod%numaor(is)=atwf_mod%numaor(is)+(2*l+1)
             ENDDO
          ENDDO
          CALL mulliken(2,c0,nstate,tau0,achrg)
       ELSEIF (cnti%icmet.EQ.2) THEN
          ! ..Hirshfeld Charges
          IF (pslo_com%tivan) CALL rnlsm(c0(:,1:crge%n),crge%n,1,1,.FALSE.)
          CALL rhoofr(c0(:,1:crge%n),rhoe,v,crge%n)


          ALLOCATE(rbuff(fpar%nnr1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          CALL atomc(rhoe(:,1),rbuff,tau0,achrg)

          DEALLOCATE(rbuff,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)



       ELSEIF (cnti%icmet.EQ.3) THEN
          ! ..ESP Charges
          IF (pslo_com%tivan) CALL rnlsm(c0(:,1:crge%n),crge%n,1,1,.FALSE.)
          CALL rhoofr(c0(:,1:crge%n),rhoe,v,crge%n)
          CALL eicalc(eivps,eirop)
          DO i=1,fpar%nnr1
             v(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(v,.FALSE.,parai%allgrp)
          DO ig=1,ncpw%nhg
             vtemp(ig) = v(nzh(ig))
          ENDDO
          CALL elstpo(vtemp,eirop,eivps,v,qphi)
          CALL dcopy(2*ncpw%nhg,v,1,vtemp,1)
          CALL zeroing(v)!,maxfft)
          !CDIR NODEP
          DO ig=1,ncpw%nhg
             v(indz(ig)) = CONJG(vtemp(ig))
             v(nzh(ig))  = vtemp(ig)
          ENDDO
          IF (geq0.AND.isos1%tclust) THEN
             v(nzh(1)) = vtemp(1)
          ELSEIF (geq0) THEN
             v(nzh(1)) = CMPLX(0._real_8,0._real_8,kind=real_8)
          ENDIF
          CALL  invfftn(v,.FALSE.,parai%allgrp)
          DO i=1,fpar%nnr1
             rhoe(i,1)=REAL(v(i))
          ENDDO
          mlen=fpar%nnr1/2 + 1
          ippc=0
          ALLOCATE(isel(mlen),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL selectp(isel,tau0,ippc)
          mlen=(ions1%nat+1)*ippc
          ALLOCATE(efield(mlen),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i=1,ippc
             efield(i)=rhoe(isel(i),1)
          ENDDO
          CALL atfield(efield,v,vtemp,reirop,qphi,isel,ippc)
          CALL espsolv(efield,achrg,ippc)
          DEALLOCATE(isel,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(efield,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSE
          CALL stopgm('INTERFACE_WRITE','INVALID ICMET',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! 
    ! ..writeout all data to interface-file
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            OPEN(5,file=int_filen,status='UNKNOWN')
       IF (paral%io_parent)&
            REWIND(5)
       iat=0
       IF (paral%io_parent)&
            WRITE(5,'(G15.6,TR1,G15.6)') ener_com%etot,ener_com%eext
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(5,'(G15.6,TR1,G15.6,TR1,G15.6,TR1,G15.6)')&
                  fion(1,ia,is),fion(2,ia,is),fion(3,ia,is),achrg(iat)
          ENDDO
       ENDDO
       ! 
       ! ..Calculate mean E-field and writeout dipole-moment
       ! 
       dz=epot2%boxdum(6)-epot2%boxdum(5)
       dy=epot2%boxdum(4)-epot2%boxdum(3)
       dx=epot2%boxdum(2)-epot2%boxdum(1)
       ! ..E-Field X-Komponente
       k0(0) = epot2%boxdum(1) + dx/3.0_real_8
       k0(1) = epot2%boxdum(3) + dy/2.0_real_8
       k0(2) = epot2%boxdum(5) + dz/2.0_real_8
       k1(0) = epot2%boxdum(1) + 2.0_real_8*dx/3.0_real_8
       k1(1) = epot2%boxdum(3) + dy/2.0_real_8
       k1(2) = epot2%boxdum(5) + dz/2.0_real_8
       ef(0) = 3.0_real_8*(clcpot(k1) - clcpot(k0))/dx
       ! ..E-Field Y-Komponente
       k0(0) = epot2%boxdum(1) + dx/2.0_real_8
       k0(1) = epot2%boxdum(3) + dy/3.0_real_8
       k0(2) = epot2%boxdum(5) + dz/2.0_real_8
       k1(0) = epot2%boxdum(1) + dx/2.0_real_8
       k1(1) = epot2%boxdum(3) + 2.0_real_8*dy/3.0_real_8
       k1(2) = epot2%boxdum(5) + dz/2.0_real_8
       ef(1) = 3.0_real_8*(clcpot(k1) - clcpot(k0))/dy
       ! ..E-Field Z-Komponente
       k0(0) = epot2%boxdum(1) + dx/2.0_real_8
       k0(1) = epot2%boxdum(3) + dy/2.0_real_8
       k0(2) = epot2%boxdum(5) + dz/3.0_real_8
       k1(0) = epot2%boxdum(1) + dx/2.0_real_8
       k1(1) = epot2%boxdum(3) + dy/2.0_real_8
       k1(2) = epot2%boxdum(5) + 2.0_real_8*dz/3.0_real_8
       ef(2) = 3.0_real_8*(clcpot(k1) - clcpot(k0))/dz
       eef   = SQRT(ef(0)*ef(0) + ef(1)*ef(1) + ef(2)*ef(2))
       ! ..Dipole Moment
       IF (prop1%locd) THEN
          d1=-au_deb*dmomlo(1,1)
          d2=-au_deb*dmomlo(2,1)
          d3=-au_deb*dmomlo(3,1)
          dd=SQRT(d1*d1+d2*d2+d3*d3)
       ELSE IF (prop1%ldip) THEN
          d1=-au_deb*moment%dmom(1)
          d2=-au_deb*moment%dmom(2)
          d3=-au_deb*moment%dmom(3)
          dd=SQRT(d1*d1+d2*d2+d3*d3)
       ELSE
          d1=-99999.9_real_8
          d2=-99999.9_real_8
          d3=-99999.9_real_8
          dd=-99999.9_real_8
       ENDIF
       ! ..Write it to file
       IF (paral%io_parent)&
            WRITE(5,'(F12.5,TR1,F12.5,TR1,F12.5,TR1,F12.5,A)')&
            EF(0),EF(1),EF(2),EEF,' E-Field'
       IF (paral%io_parent)&
            WRITE(5,'(F12.5,TR1,F12.5,TR1,F12.5,TR1,F12.5,A)')&
            D1,D2,D3,DD,' Debye'
       ! ..
       ! ..added by biswas (dec/2004) to write the forces on MM atoms by QM in qmoutput.out
       IF (cnti%iftype.EQ.2) THEN
          ! ..
          DO inr=1,epot1%nnear
             idnr=epot1%idnear(inr)
             IF (paral%io_parent)&
                  WRITE(5,'(I6,TR1,E15.6,TR1,E15.6,TR1,E15.6)')&
                  IDNR,epot5%gmx_fnear(1,INR),epot5%gmx_fnear(2,INR),epot5%gmx_fnear(3,INR)
          ENDDO
          ! ..
          IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
             DO inr=1,epot3%nintm
                idnr=epot3%idintm(inr)
                IF (paral%io_parent)&
                     WRITE(5,'(I6,TR1,E15.6,TR1,E15.6,TR1,E15.6)')&
                     IDNR,epot5%gmx_fintm(1,INR),epot5%gmx_fintm(2,INR),epot5%gmx_fintm(3,INR)
             ENDDO
          ENDIF
          ! ..
          IF (MOD(epot3%stpno,epot3%outmf).EQ.0.AND.epot3%noutr.GT.0) THEN
             DO inr=1,epot3%noutr
                idnr=epot3%idoutr(inr)
                IF (paral%io_parent)&
                     WRITE(5,'(I6,TR1,E15.6,TR1,E15.6,TR1,E15.6)')&
                     IDNR,epot5%gmx_foutr(1,INR),epot5%gmx_foutr(2,INR),epot5%gmx_foutr(3,INR)
             ENDDO
          ENDIF
          ! ..       
       ENDIF
       IF (paral%io_parent)&
            CLOSE(5)
       ! 
       ! ..Create the EGO/GMXCONTINUE file
       ! 
       IF (cnti%iftype.EQ.1) THEN
          str=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'EGOCONTINUE'
          CALL xstring(str,ia,ie)
          line=str(ia:ie)
          IF (paral%io_parent)&
               OPEN(5,file=line,status='UNKNOWN')
          IF (paral%io_parent)&
               REWIND(5)
          IF (paral%io_parent)&
               WRITE(5,'(A)') 'EGOCONTINUE'
          IF (paral%io_parent)&
               CLOSE(5)
          ! ..added
       ELSEIF (cnti%iftype.EQ.2) THEN
          str=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'GMXCONTINUE'
          CALL xstring(str,ia,ie)
          line=str(ia:ie)
          IF (paral%io_parent)&
               OPEN(5,file=line,status='UNKNOWN')
          IF (paral%io_parent)&
               REWIND(5)
          IF (paral%io_parent)&
               WRITE(5,'(A)') 'GMXCONTINUE'
          IF (paral%io_parent)&
               CLOSE(5)
       ELSE
          CALL stopgm('INTERFACE','NOT PROGRAMMED',& 
               __LINE__,__FILE__)
       ENDIF
       CALL wrgeof(tau0,fion)
       IF (paral%io_parent)&
            WRITE(6,*) ' INTERFACE| FORCES WRITTEN TO FILE '
    ENDIF
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(reirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(qphi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interface_write
  ! ==================================================================
  ! 
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE egostop
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(len=80)                        :: line, str
    INTEGER                                  :: ia, ie

! ==--------------------------------------------------------------==

    IF (cnti%iftype.EQ.1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(//,A)') ' Creating EGOEXIT file'
       str=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'EGOEXIT'
       CALL xstring(str,ia,ie)
       line=str(ia:ie)
       IF (paral%io_parent)&
            OPEN(5,file=line,status='UNKNOWN')
       IF (paral%io_parent)&
            REWIND(5)
       IF (paral%io_parent)&
            WRITE(5,'(A)') 'EGOEXIT'
       IF (paral%io_parent)&
            CLOSE(5)
    ELSEIF (cnti%iftype.EQ.2) THEN
       IF (paral%io_parent)&
            WRITE(6,'(//,A)') ' Creating GMXEXIT file'
       str=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'GMXEXIT'
       CALL xstring(str,ia,ie)
       line=str(ia:ie)
       IF (paral%io_parent)&
            OPEN(5,file=line,status='UNKNOWN')
       IF (paral%io_parent)&
            REWIND(5)
       IF (paral%io_parent)&
            WRITE(5,'(A)') 'GMXEXIT'
       IF (paral%io_parent)&
            CLOSE(5)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE egostop
  ! ==================================================================
  SUBROUTINE myproppt(c0,tau0,rhoe,v)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: tau0(:,:,:), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v(maxfftn)

    CHARACTER(*), PARAMETER                  :: procedureN = 'myproppt'

    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:)
    INTEGER                                  :: i, ierr, ir, j, ji
    LOGICAL                                  :: tinfo
    REAL(real_8), ALLOCATABLE                :: fback(:)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,63("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",9X,A,19X,"*")')&
            ' INTERFACE| PROPERTY CALCULATIONS'
       IF (paral%io_parent)&
            WRITE(6,'(1X,63("*"),/)')
    ENDIF

    tinfo=.FALSE.
    CALL propin(tinfo)
    ! ==--------------------------------------------------------------==
    prop2%numorb = MAX(crge%n,prop2%numorb)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE DIPOLE MOMENTS                                    ==
    ! ==--------------------------------------------------------------==
    IF (prop1%ldip.OR.prop1%locd.OR.prop1%lext) THEN
       CALL phfac(tau0)
       CALL rnlsm(c0(:,1:prop2%numorb),prop2%numorb,1,1,.FALSE.)
       ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF (prop1%lext) THEN
          ALLOCATE(fback(prop2%numorb),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(prop2%numorb,crge%f,1,fback,1)
          DO i=1,extd%nexdip
             DO j=1,prop2%numorb
                ji=(i-1)*prop2%numorb+j
                crge%f(j,1)=fback(j)-focc(ji)
             ENDDO
             CALL difrho(c0,rhoe,v,prop2%numorb)
             CALL zeroing(eirop)!CALL  fwfftn(v,.FALSE.)
             CALL exdipo(i,tau0,v,eirop)
          ENDDO
          CALL dcopy(prop2%numorb,fback,1,crge%f,1)
       ENDIF
       IF (prop1%ldip.OR.prop1%locd) THEN
          CALL rhoofr(c0(:,1:crge%n),rhoe,v,crge%n)
          CALL eicalc(eivps,eirop)
          DO ir=1,fpar%nnr1
             v(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          ENDDO
          CALL  fwfftn(v,.FALSE.,parai%allgrp)
          IF (prop1%ldip) CALL dipo(tau0,eirop,v)
          IF (prop1%locd) CALL lodipo(eirop,v)
       ENDIF
       IF (paral%parent) THEN
          CALL prdip
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE myproppt
  ! ==================================================================
  SUBROUTINE give_scr_interface_write(linterw,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: linterw
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lesp, lmulliken, lrhoofr, &
                                                lrnlsm, nstate

    nstate=crge%n
    lrhoofr=0
    lrnlsm=0
    lesp=0
    lmulliken=0
    IF (cnti%icmet.EQ.2 .OR. cnti%icmet.EQ.3) THEN
       CALL give_scr_rhoofr(lrhoofr,tag)
       CALL give_scr_espc(lesp,tag)
    ENDIF
    IF (cnti%icmet.EQ.1 .OR. cnti%icmet.EQ.4) THEN
       CALL give_scr_mulliken(lmulliken,tag,nstate)
    ELSEIF (cnti%icmet.EQ.2 .OR. cnti%icmet.EQ.3) THEN
       IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    ENDIF
    linterw=MAX(lrhoofr,lesp,lrnlsm,lmulliken)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_interface_write
  ! ==================================================================
  SUBROUTINE give_scr_myproppt(lmyproppt,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmyproppt
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lrhoofr, lrhopri, lrnlsm, num

    lrhopri=0
    lrnlsm=0
    lrhoofr=0
    lmyproppt=0
    IF (prop1%ldip.OR.prop1%locd.OR.prop1%lext) THEN
       num=MAX(crge%n,prop2%numorb)
       CALL give_scr_rhopri(lrhopri,tag,num)
       CALL give_scr_rnlsm(lrnlsm,tag,num,.FALSE.)
       ! EIVPS (2*NHG) EIROP (2*NHG)
       lmyproppt=4*ncpw%nhg
       IF (prop1%ldip.OR.prop1%locd) THEN
          CALL give_scr_rhoofr(lrhoofr,tag)
       ENDIF
    ENDIF
    lmyproppt=MAX(lmyproppt,lrnlsm,lrhoofr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_myproppt
  ! ==================================================================
  SUBROUTINE give_scr_interface(linter,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: linter
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcopot, lforcedr, linterw, &
                                                lmyproppt, lupd, nstate

    nstate=crge%n
    lcopot=0
    lmyproppt=0
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    CALL give_scr_updwf(lupd,tag,nstate,.FALSE.)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.TRUE.,.TRUE.)
    IF (cntl%proper) CALL give_scr_myproppt(lmyproppt,tag)
    CALL give_scr_interface_write(linterw,tag)
    linter=MAX(lcopot,lupd,lforcedr,lmyproppt,linterw)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_interface
  ! ==================================================================

  ! ==================================================================
  ! ..added in dec/2004: biswas; calculates forces on  gmx atoms due to QM.
  ! ==================================================================
  SUBROUTINE gmx_el_forces(c0,tau0,rhoe,v,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! Variables
    REAL(real_8)                             :: tau0(:,:,:), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v(maxfft)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gmx_el_forces'

    INTEGER                                  :: ikind, inr, ndim
    REAL(real_8)                             :: gmx_fmm(3,maxoutr)

! ==--------------------------------------------------------------==
! ..To create rho(r):

    DO ikind=1,nkpoint
       CALL rnlsm(c0(:,:,ikind),nstate,1,ikind,.FALSE.)
    ENDDO
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,v,nstate)
    ELSE
       CALL rhoofr(c0(:,:,1),rhoe,v,nstate)
    ENDIF
    ! ..End of creation of rho(r)
    ! ==--------------------------------------------------------------==
    ! ..Calls to calculate forces for three different layers: near, intermediate, outer
    ! ==--------------------------------------------------------------==

    ! NOTE: The storage extent of the dummy argument exceeds that of the
    ! actual argument.  [KONEAR]
    ! KONEAR is defined in epot.mod as (:,:), but in GMX_EL_FORCE dummy
    ! argument is TAU0(:,:,:)
    ! can compile only with this error suppressed [ifort]

    CALL stopgm(procedureN,'fix me',&
         __LINE__,__FILE__)
    !vw this is buggus      call gmx_el_force(epot1%nnear,epot2%konear,epot2%datnear,tau0,rhoe,gmx_fmm)

    DO inr=1,epot1%nnear
       DO ndim=1,3
          !vw commented next line and added stop
          epot5%gmx_fnear(ndim,inr)=gmx_fmm(ndim,inr)
          CALL stopgm(procedureN,'array not initialized', &
               __LINE__,__FILE__)
       ENDDO
    ENDDO
    ! ..
    IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
       CALL stopgm(procedureN,'fix me',&
            __LINE__,__FILE__)
       !vw this is buggus        call gmx_el_force(epot3%nintm,epot4%kointm,epot4%qintm,tau0,rhoe,gmx_fmm)
       DO inr=1,epot3%nintm
          DO ndim=1,3
             epot5%gmx_fintm(ndim,inr)=gmx_fmm(ndim,inr)
          ENDDO
       ENDDO
    ENDIF
    ! ..
    IF (MOD(epot3%stpno,epot3%intmf).EQ.0.AND.epot3%nintm.GT.0) THEN
       CALL gmx_el_force(epot3%noutr,epot4%kooutr,epot4%qoutr,tau0,rhoe,gmx_fmm)
       DO inr=1,epot1%nnear
          DO ndim=1,3
             epot5%gmx_foutr(ndim,inr)=gmx_fmm(ndim,inr)
          ENDDO
       ENDDO
    ENDIF
    ! ..
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gmx_el_forces
  ! 
  ! ==================================================================
  ! .....Calculates QM-2-MM electrostatic forces on MM atoms.
  ! ==================================================================
  SUBROUTINE gmx_el_force(nmm,komm,datmm,tau0,rhoe,gmx_fmm)
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: nmm
    REAL(real_8) :: komm(3,maxoutr), datmm(5,maxoutr), tau0(:,:,:), &
      rhoe(fpar%nnr1,clsd%nlsd), gmx_fmm(3,maxoutr)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gmx_el_force'
    REAL(real_8), PARAMETER :: force_const = 138.935485_real_8 

    INTEGER                                  :: i, ia, idmm(maxoutr), idnr, &
                                                ikr1, inr, ip, ipp, ir, is, &
                                                j, k, msgid, my_len, ndim
    REAL(real_8) :: alam = 1.3_real_8, alfa, beta, bohr_to_nm, di, dj, dk, &
      gfac, gmx_f_ip(3,maxoutr), gmx_f_ips(3,maxoutr), ko(0:2), qj_rho, &
      rffac, rij, rijq, rijs, rrj, rrjq, rrjs, xi, xj, yi, yj, zfac, zi, zj

! ==--------------------------------------------------------------==

    bohr_to_nm=0.1_real_8/fbohr
    ! ..
    dk=(epot2%boxdum(6)-epot2%boxdum(5))/REAL(spar%nr3s,kind=real_8)
    dj=(epot2%boxdum(4)-epot2%boxdum(3))/REAL(spar%nr2s,kind=real_8)
    di=(epot2%boxdum(2)-epot2%boxdum(1))/REAL(spar%nr1s,kind=real_8)
    ! ..   
    gfac=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ikr1 = MOD(parm%nr1+1,2)
    ! ..
    DO inr=1,nmm
       DO ndim=1,3
          gmx_fmm(ndim,inr)=0._real_8
          gmx_f_ip(ndim,inr)=0._real_8
          gmx_f_ips(ndim,inr)=0._real_8
       ENDDO
    ENDDO
    ! ..   
    ! ----------------------------------------------------------------------C
    ! ..Calculate contributions from densities on the grid.
    ! ----------------------------------------------------------------------C
    DO inr=1,nmm

       !vw commented next line and added stop 
       idnr=idmm(inr)
       CALL stopgm(procedureN,'array not initialized', &
            __LINE__,__FILE__)


       xj=komm(1,inr)
       yj=komm(2,inr)
       zj=komm(3,inr)
       ! ..   
       alfa = alam*datmm(3,inr)
       beta = 2._real_8*alfa
       ! ..   
       ir = 0
       DO k=1,fpar%kr3s
          ko(2)=epot2%boxdum(5)+REAL(k-1,kind=real_8)*dk
          DO j=1,fpar%kr2s
             ko(1)=epot2%boxdum(3)+REAL(j-1,kind=real_8)*dj
             DO i=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)+ikr1
                ko(0)=epot2%boxdum(1)+REAL(i-1,kind=real_8)*di
                ! ..   
                ir=ir+1
                ! ..   
                IF ((ir.GT.fpar%nnr1).AND.paral%io_parent)&
                     WRITE(6,*)'Error in looping in GMX_EL'
                ! ..   
                rrjs = (xj-ko(0))**2+(yj-ko(1))**2+(zj-ko(2))**2
                rrj = SQRT(rrjs)
                rrjq = rrj*rrjs
                qj_rho = datmm(1,inr)*rhoe(ir,1)
                ! ..   
                rffac=1._real_8/rrjq-EXP(-beta*rrj)*(1._real_8/rrjq+beta/rrjs+&
                     alfa*beta/rrj)
                ! RVFAC=1./RRJ - EXP(-BETA*RRJ)*(1./RRJ + ALFA) 
                ! ----------------------------------------------------------------------C
                gmx_f_ip(1,inr)=gmx_f_ip(1,inr)-qj_rho*(xj-ko(0))*rffac
                gmx_f_ip(2,inr)=gmx_f_ip(2,inr)-qj_rho*(yj-ko(1))*rffac
                gmx_f_ip(3,inr)=gmx_f_ip(3,inr)-qj_rho*(zj-ko(2))*rffac
                ! ----------------------------------------------------------------------C       
             ENDDO
          ENDDO
       ENDDO
       ! == Normalization for rho === MM-charge is already multiplied.
       gmx_f_ip(1,inr)=gmx_f_ip(1,inr)*gfac
       gmx_f_ip(2,inr)=gmx_f_ip(2,inr)*gfac
       gmx_f_ip(3,inr)=gmx_f_ip(3,inr)*gfac
       ! ..   
    ENDDO
    ! ----------------------------------------------------------------------C
    ! ..
    DO ipp=1,parai%nproc
       ip=parap%pgroup(ipp)
       msgid = ip
       my_len=nmm*3
       IF (paral%parent) THEN
          IF (ip.EQ.parai%me) THEN
             DO inr=1,nmm
                DO ndim=1,3
                   gmx_f_ips(ndim,inr) = gmx_f_ips(ndim,inr) +&
                        gmx_f_ip(ndim,inr)
                ENDDO
             ENDDO
          ELSE
             msgid=1
             CALL mp_recv(gmx_f_ip,my_len,ip,msgid,parai%allgrp)
             DO inr=1,nmm
                DO ndim=1,3
                   gmx_f_ips(ndim,inr) = gmx_f_ips(ndim,inr) +&
                        gmx_f_ip(ndim,inr)
                ENDDO
             ENDDO
          ENDIF
          ! ..   
       ELSE
          IF (ip.EQ.parai%me) THEN
             msgid=1
             CALL mp_send(gmx_f_ip,my_len,parap%pgroup(1),msgid,parai%allgrp)
          ENDIF
       ENDIF
    ENDDO
    ! ..   
    IF (paral%parent) THEN
       DO inr=1,nmm
          DO ndim=1,3
             gmx_fmm(ndim,inr) = gmx_f_ips(ndim,inr)
          ENDDO
       ENDDO
    ENDIF
    ! ..   
    ! ..
    ! ----------------------------------------------------------------------C
    ! ..Calculate & add contributions from QM-core ions.
    ! ----------------------------------------------------------------------C
    ! ..
    IF (paral%parent) THEN
       ! ..
       DO inr=1,nmm
          idnr=idmm(inr)
          xj=komm(1,inr)
          yj=komm(2,inr)
          zj=komm(3,inr)
          ! ..   
          alfa = alam*datmm(3,inr)
          beta = 2._real_8*alfa
          ! ..   
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                xi=tau0(1,ia,is)
                yi=tau0(2,ia,is)
                zi=tau0(3,ia,is)
                rijs = (xj-xi)**2+(yj-yi)**2+(zj-zi)**2
                rij = SQRT(rijs)
                rijq = rij*rijs
                zfac=1._real_8/rijq-EXP(-beta*rij)*&
                     (1._real_8/rijq+beta/rijs+alfa*beta/rij)
                ! ZVAC=1./RIJ - EXP(-BETA*RIJ)*(1./RIJ + ALFA) 
                ! ..
                gmx_fmm(1,inr)=gmx_fmm(1,inr)+datmm(1,inr)*ions0%zv(is)&
                     *(xj-xi)*zfac
                gmx_fmm(2,inr)=gmx_fmm(2,inr)+datmm(1,inr)*ions0%zv(is)&
                     *(yj-yi)*zfac
                gmx_fmm(3,inr)=gmx_fmm(3,inr)+datmm(1,inr)*ions0%zv(is)&
                     *(zj-zi)*zfac
                ! ..   
             ENDDO
          ENDDO
          ! ..   
          gmx_fmm(1,inr)=gmx_fmm(1,inr)*force_const/bohr_to_nm**2
          gmx_fmm(2,inr)=gmx_fmm(2,inr)*force_const/bohr_to_nm**2
          gmx_fmm(3,inr)=gmx_fmm(3,inr)*force_const/bohr_to_nm**2
          ! ..   
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gmx_el_force


END MODULE egointer_utils
