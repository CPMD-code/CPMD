#include "cpmd_global.h"

MODULE wv30_utils
  USE benc,                            ONLY: ibench
  USE cdft_utils,                      ONLY: wcdft_restart
  USE cdftmod,                         ONLY: cdfthda,&
                                             cdftlog
  USE cell,                            ONLY: cell_com
  USE clas,                            ONLY: clas3,&
                                             clasc,&
                                             clasv
  USE cotr,                            ONLY: cnsval,&
                                             cotc0
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_fpath_max,&
                                             fo_info
  USE filnmod,                         ONLY: filbod
  USE glemod,                          ONLY: glep,&
                                             glepar
  USE io_utils,                        ONLY: file_get_position,&
                                             file_open
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE lscal,                           ONLY: adtolr,&
                                             hesscr,&
                                             lallhs,&
                                             lvlhes,&
                                             nshess,&
                                             nvar
  USE machine,                         ONLY: m_system
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mm_dimmod,                       ONLY: clsaabox
  USE mm_extrap,                       ONLY: cold,&
                                             nnow,&
                                             numcold
  USE mp_interface,                    ONLY: mp_sum
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nose,                            ONLY: &
       eta, etadot, etap, etap1, etap1dot, etapdot, etapm, etapmdot, etc, &
       etcdot, lcttab, nchain, nchc, nche, nchx, nedof, nosl
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: ipcurr
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE prng,                            ONLY: prng_com
  USE readsr_utils,                    ONLY: xstring
  USE rlbfgs_io,                       ONLY: wr30lbfgs,&
                                             wr30prfo
  USE rw_linres_utils,                 ONLY: w_linres
  USE shop_rest,                       ONLY: sh03
  USE shop_rest_2,                     ONLY: c0old
  USE spin,                            ONLY: lspin2
  USE store_types,                     ONLY: &
       iface1, irec_ac, irec_cell, irec_clas, irec_co, irec_cons, &
       irec_ctrans, irec_cut, irec_eigv, irec_gle, irec_he, irec_ico, &
       irec_kpt, irec_lbfgs, irec_lrwf, irec_lrwv, irec_nas, irec_noc, &
       irec_noe, irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_nose, &
       irec_occ, irec_phes, irec_pot, irec_prfo, irec_prng, irec_pws, &
       irec_pwv, irec_rdtl, irec_rho, irec_sdp, irec_shop, irec_shopbo, &
       irec_vel, irec_wf, irec_xtrp, restart1, store1
  USE string_utils,                    ONLY: int2str
  USE symm,                            ONLY: symmi
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, dual00, maxsys, nacc, nkpt, parm, restf, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wfnio_utils,                     ONLY: w_wfnio
  USE wrintf_utils,                    ONLY: wintf
  use bicanonicalCpmd, only: biCanonicalEnsembleDo,&
    bicanonicalCpmdConfig, getNameLatestTape

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: zhwwf
  !public :: wv30
  !public :: file_update_record
  !public :: wr30pot

CONTAINS

  ! ==================================================================
  SUBROUTINE zhwwf (nw,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)
    ! ==--------------------------------------------------------------==
    ! == Write wavefunction in Version 3.0 format on unit nw          ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nw, irec(:)
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                cm(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,*), tau0(:,:,:), &
                                                taum(:,:,:), taui(:,:,:)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'zhwwf'

    CHARACTER(len=2*fo_fpath_max)            :: cpstring
    CHARACTER(len=20)                        :: fformat, latfile
    CHARACTER(len=fo_fpath_max)              :: filnam, filnam2
    INTEGER                                  :: i, i1, i2, ia, ie, iffa, &
                                                iffe, ifil, iostat, isub, j1, &
                                                j2
    INTEGER, SAVE                            :: icalls = 0
    LOGICAL                                  :: fexist

! ==--------------------------------------------------------------==

    IF (ibench(1).EQ.1) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==

    ! Construct the filename
    ! ONLY the IO_PARENT does the I/O
    IF (paral%io_parent) THEN

       CALL xstring(filbod,ia,ie)
       IF (cntl%tpath) THEN
          filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filbod(ia:ie)
       ELSE
          IF (tmw)THEN
             CALL mw_filename('LATEST_',latfile,mwi%walker_id)
          elseif (biCanonicalEnsembleDo) then
            latfile = getNameLatestTape(bicanonicalCpmdConfig)
            !JF call bicanonicalSetFileNameLatest(latfile)
          ELSE
             latfile='LATEST'
          ENDIF
          INQUIRE(file=latfile,exist=fexist)
          IF (restart1%rlate.AND.fexist) THEN
             OPEN(unit=nw,file=latfile,status='UNKNOWN',iostat=iostat)
             IF (iostat /= 0 ) CALL stopgm(procedureN,&
                  'Problem opening FILE = '//TRIM(latfile)//&
                  ' IOSTAT = '//TRIM(int2str(iostat)),& 
                  __LINE__,__FILE__)
             REWIND(unit=nw)
             READ(nw,'(A)') filnam
             READ(nw,*) icalls
             CLOSE(nw)
          ENDIF
          IF (.NOT.cdftlog%thda)THEN
             ifil=MOD(icalls,cnti%nrestf)+1
             icalls = icalls + 1
             WRITE(fformat,'(I10)') ifil
          ELSE
             IF (cdfthda%hdafirst)THEN
                fformat='STATE_1'
             ELSE
                fformat='STATE_2'
             ENDIF
          ENDIF
          CALL xstring(fformat,iffa,iffe)
          WRITE(filnam,'(A,A,A)')&
               fo_info%fpath(fo_info%iapath:fo_info%iepath),filbod(ia:ie),fformat(iffa:iffe)
       ENDIF

       CALL file_open(filnam,'UNFORMATTED','UNKNOWN',&
            'READWRITE',cntl%is_out_stream,.FALSE.,0,nw)

    ENDIF ! IO_PARENT

    CALL wv30(nw,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)

    IF (paral%io_parent) THEN
       CLOSE(unit=nw,iostat=iostat)
       IF (iostat /= 0 ) CALL stopgm(procedureN,&
            'Problem closing FILE = '//TRIM(filnam)//&
            ' IOSTAT = '//TRIM(int2str(iostat)),& 
            __LINE__,__FILE__)
    ENDIF

    IF (paral%io_parent) THEN
       CALL xstring(filnam,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(38,65-(i2-i1)),',A)'
       WRITE(6,fformat)&
            ' RESTART INFORMATION WRITTEN ON FILE ',filnam(i1:i2)
       IF (restf%nrcopy.GT.1) THEN
          DO i=2,restf%nrcopy
             ifil=MOD(icalls,cnti%nrestf)+1
             icalls = icalls + 1
             WRITE(fformat,'(I10)') ifil
             CALL xstring(fformat,iffa,iffe)
             WRITE(filnam2,'(A,A,A)')&
                  fo_info%fpath(fo_info%iapath:fo_info%iepath),filbod(ia:ie),fformat(iffa:iffe)
             CALL xstring(filnam2,j1,j2)
             cpstring='cp '//filnam(i1:i2)//'  '//filnam2(j1:j2)
             CALL m_system(cpstring)
             WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(38,65-(j2-j1)),',A)'
             WRITE(6,fformat)&
                  ' RESTART INFORMATION WRITTEN ON FILE ',filnam2(j1:j2)
          ENDDO
       ENDIF
       IF (.NOT.cntl%tpath) THEN
          OPEN(unit=nw,file=latfile,status='UNKNOWN',iostat=iostat)
          IF (iostat /= 0 ) CALL stopgm(procedureN,&
               'Problem opening FILE = '//TRIM(filnam)//&
               ' IOSTAT = '//TRIM(int2str(iostat)),& 
               __LINE__,__FILE__)
          REWIND(unit=nw)
          WRITE(nw,'(A)') filnam
          WRITE(nw,*) icalls
          CLOSE(nw)
       ENDIF
    ENDIF ! IO_PARENT
    ! INTERFACE FILE
    IF (iface1%intwrite) THEN
       CALL wintf(nw,c0,tau0)
    ENDIF
    IF (paral%parent)THEN
       IF (cntl%cdft)CALL wcdft_restart()
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zhwwf
  ! ==================================================================
  SUBROUTINE wv30(nw,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nw, irec(:)
    COMPLEX(real_8)                          :: c0(*), cm(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,*), tau0(:,:,:), &
                                                taum(:,:,:), taui(:,:,:)
    INTEGER                                  :: nfi

    CHARACTER(len=80)                        :: str
    INTEGER                                  :: i, ia, ierror, ik, ip, irec0, &
                                                irecord, is, isection, j, k, &
                                                l, nkpts0
    INTEGER(int_8)                           :: fpos_beg, fpos_end, fpos_null
    LOGICAL                                  :: lsavc0
    REAL(real_8)                             :: prngin(parai%nproc,2,3), &
                                                prngout(parai%nproc,2,3)

! ==--------------------------------------------------------------==
! ==                     GENERAL FILE FORMAT                      ==
! ==--------------------------------------------------------------==
! ==  Section  1: Header                                          ==
! ==  Section  2: Symmetry and Cell info                          ==
! ==  Section  3: Number of atomic species and atoms per species  ==
! ==  Section  4: Atomic coordinates                              ==
! ==  Section  5: Atomic velocities                               ==
! ==  Section  6: Initial atomic coordinates                      ==
! ==  Section  7: Cutoff, # of electrons, grid                    ==
! ==  Section  8: States, dual, plane waves                       ==
! ==  Section  9: PW Coefficients                                 ==
! ==  Section 10: PW Velocities                                   ==
! ==  Section 11: Accumulators                                    ==
! ==  Section 12: Nose thermostats general info                   ==
! ==  Section 13: Nose thermostats electrons                      ==
! ==  Section 14: Nose thermostats ions                           ==
! ==  Section 15: Nose thermostats ions (ULTRA)                   ==
! ==  Section 16: Nose thermostats ions (MASSIVE)                 ==
! ==  Section 17: Nose thermostats cell                           ==
! ==  Section 18: Potential                                       ==
! ==  Section 19: PW single states                                ==
! ==  Section 20: H Matrices (cell parameters)                    ==
! ==  Section 21: K-Points                                        ==
! ==  Section 22: Electron density                                ==
! ==  Section 23: Occupation numbers                              ==
! ==  Section 24: Eigenvalues                                     ==
! ==  Section 25: Classical Particles (coordinates and velocities)==
! ==  Section 26: LinRes PW Coefficients                          ==
! ==  Section 27: LinRes PW Velocities                            ==
! ==  Section 28: Partial Hessian (for microiterative TS search)  ==
! ==  Section 29: P-cntl%rfo status and work arrays                    ==
! ==  Section 30: L-cntl%bfgs history and status                       ==
! ==  Section 31: Adaptive tolerance status                       ==
! ==  Section 32: Constraints values                              ==
! ==  Section 33: Cell translation vector in QM/MM runs           ==
! ==  Section 34: Local temperature parameters                    ==
! ==  Section 35: Surface Hopping I                               ==
! ==  Section 36: Surface Hopping II                              ==
! ==  Section 37: Pseudo random number generator                  ==
! ==  Section 38: (generalized) Langevin equation                 ==
! ==  Section 39 - 98 : empty                                     ==
! ==  Section 99: Wavefunction history for cntl%textrap (must be last) ==
! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

    ierror=0
    ! ==--------------------------------------------------------------==
    IF (cntl%tpath.AND.cntl%tpimd) THEN
       ip=ipcurr
    ELSEIF (tmw)THEN
       ip=mwi%walker_id
    ELSE
       ip=1
    ENDIF
    irec0=0
    fpos_null=0
    IF (paral%io_parent) THEN
       CALL file_get_position(nw,.FALSE.,fpos_beg)
       ! ==------------------------------------------------------------==
       ! Section 1
       isection=1
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       str=' STUTTGART VERSION 3.0'
       IF (cntl%is_out_stream) str=TRIM(str)//' STREAM'
       WRITE (nw) str
       CALL file_get_position(nw,.FALSE.,fpos_end)
       fpos_beg=fpos_end
       ! ==------------------------------------------------------------==
       ! Section 2
       isection=2
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_cell).NE.0) THEN
          irecord=2
          WRITE (nw) irecord,fpos_null
          WRITE (nw) parm%ibrav,symmi%indpg
          WRITE (nw) (cell_com%celldm(i),i=1,6)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 3
       isection=3
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_nas).NE.0) THEN
          irecord=2
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) ions1%nsp
          WRITE (unit=nw,iostat=ierror)&
               (ions0%na(i),i=1,ions1%nsp)! ,(EL(IATYP(I)),I=1,NSP)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 4
       isection=4
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_co).NE.0) THEN
          irecord=ions1%nat
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO is = 1 , ions1%nsp
             DO ia = 1 , ions0%na(is)
                WRITE (unit=nw,iostat=ierror) (tau0(j,ia,is),j=1,3)
             ENDDO
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 5
       isection=5
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_vel).NE.0) THEN
          irecord=ions1%nat
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO is = 1 , ions1%nsp
             DO ia = 1 , ions0%na(is)
                WRITE (unit=nw,iostat=ierror) (taum(j,ia,is),j=1,3)
             ENDDO
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 6
       isection=6
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_ico).NE.0) THEN
          irecord=ions1%nat
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO is = 1 , ions1%nsp
             DO ia = 1 , ions0%na(is)
                WRITE (unit=nw,iostat=ierror) (taui(j,ia,is),j=1,3)
             ENDDO
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 7
       isection=7
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_cut).NE.0) THEN
          irecord=1
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               cntr%ecut,dual00%cdual,dual00%dual,NINT(crge%nel),spar%nr1s,spar%nr2s,spar%nr3s
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 8
       isection=8
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_sdp).NE.0) THEN
          irecord=1
          ! If gamma point nkpts0=0
          IF (tkpts%tkpnt) THEN
             nkpts0=nkpt%nkpts
          ELSE
             nkpts0=0
          ENDIF
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) crge%n,nkpts0,spar%ngwks,spar%ngwls,spar%nhgs,spar%nhgls
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    ! End of IF(PARENT). Test of IERROR
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 9
    isection=9
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION',isectioN
    ENDIF
    IF (irec(irec_wf).NE.0) THEN
       CALL w_wfnio(nw,ierror,nstate,c0,tau0,'C0',irecord)
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 10
    isection=10
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_pwv).NE.0) THEN
       CALL wr30wfn(nw,ierror,nstate,cm,tau0,'NIL',irecord)
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 11
    isection=11
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_ac).NE.0) THEN
       irecord=2
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) nfi,nacc
          WRITE (unit=nw,iostat=ierror) ( acc ( i ) , i = 1 , nacc )
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 12
    isection=12
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_nose).NE.0) THEN
       irecord=2
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) nche,nedof,nchain,nchc
          WRITE (unit=nw,iostat=ierror) nosl%tultra,nosl%tmnose
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 13
    isection=13
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_noe).NE.0) THEN
       IF (cntl%bsymm)THEN
          irecord=2
          IF (paral%io_parent) THEN
             WRITE (unit=nw,iostat=ierror) irecord,fpos_null
             WRITE (unit=nw,iostat=ierror)&
                  (eta(i,1),etadot(i,1),i=1,nchx)
             WRITE (unit=nw,iostat=ierror)&
                  (eta(i,2),etadot(i,2),i=1,nchx)
          ENDIF
       ELSE
          irecord=1
          IF (paral%io_parent) THEN
             WRITE (unit=nw,iostat=ierror) irecord,fpos_null
             WRITE (unit=nw,iostat=ierror)&
                  (eta(i,ip),etadot(i,ip),i=1,nchx)
          ENDIF
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 14
    isection=14
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_nop1).NE.0) THEN
       irecord=1
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               (etap1(i,ip),etap1dot(i,ip),i=1,nchx)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 15
    isection=15
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_nop2).NE.0) THEN
       irecord=1
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               ((etap(j,i,ip),etapdot(j,i,ip),i=1,nchx),j=1,maxsys%nsx)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 16
    isection=16
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_nop3).NE.0) THEN
       irecord=1
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               ((etapm(k,l,ip),etapmdot(k,l,ip),k=1,3*maxsys%nax*maxsys%nsx),l=1,nchx)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 17
    isection=17
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_noc).NE.0) THEN
       irecord=1
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               (etc(i,ip),etcdot(i,ip),i=1,nchx)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 18
    isection=18
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_pot).NE.0) THEN
       IF (cntl%tlsd) THEN
          irecord=2*spar%nr1s
       ELSEIF (lspin2%tlse) THEN
          irecord=4*spar%nr1s
       ELSE
          irecord=spar%nr1s
       ENDIF
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       CALL wr30pot(nw,ierror,potr(1,1))
       IF (cntl%tlsd) THEN
          CALL wr30pot(nw,ierror,potr(1,2))
       ENDIF
       IF (lspin2%tlse) THEN
          CALL wr30pot(nw,ierror,potr(1,2))
          CALL wr30pot(nw,ierror,potr(1,3))
          CALL wr30pot(nw,ierror,potr(1,4))
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 19
    isection=19
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_pws).NE.0) THEN
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       ! deb    CALL WR30WFN(NW,IERROR,1,CSTATE,TAU0,'NIL')
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       ! Section 20
       isection=20
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_he).NE.0) THEN
          irecord=6
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) metr_com%ht
          WRITE (unit=nw,iostat=ierror) metr_com%htm1
          WRITE (unit=nw,iostat=ierror) metr_com%htvel
          WRITE (unit=nw,iostat=ierror) metr_com%htfor
          WRITE (unit=nw,iostat=ierror) metr_com%ht0
          WRITE (unit=nw,iostat=ierror) metr_com%htm10
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 21
       isection=21
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_kpt).NE.0) THEN
          irecord=nkpt%nkpts
          IF (.NOT.tkpts%tkpnt) THEN
             irecord=0
          ENDIF
          IF (tkpts%tmonkp) THEN
             ! We add the indexes of the Monkhorst-Pack mesh.
             WRITE (unit=nw,iostat=ierror) irecord,fpos_null,kpts_com%nk1,kpts_com%nk2,kpts_com%nk3
          ELSE
             WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          ENDIF
          DO ik=1,irecord
             WRITE(unit=nw,iostat=ierror) (rk(j,ik),j=1,3),wk(ik)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 22
    isection=22
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_rho).NE.0) THEN
       IF (cntl%tlsd) THEN
          irecord=2*spar%nr1s
       ELSEIF (lspin2%tlse) THEN
          irecord=4*spar%nr1s
       ELSE
          irecord=spar%nr1s
       ENDIF
       IF (paral%io_parent) THEN
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       CALL wr30pot(nw,ierror,rhoo(1,1))
       IF (cntl%tlsd) THEN
          CALL wr30pot(nw,ierror,rhoo(1,2))
       ENDIF
       IF (lspin2%tlse) THEN
          CALL wr30pot(nw,ierror,rhoo(1,2))
          CALL wr30pot(nw,ierror,rhoo(1,3))
          CALL wr30pot(nw,ierror,rhoo(1,4))
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Section 23
    IF (paral%io_parent) THEN
       isection=23
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_occ).NE.0) THEN
          irecord=nkpt%nkpts
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO ik=1,nkpt%nkpts
             WRITE (unit=nw,iostat=ierror) (crge%f(i,ik),i=1,nstate)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 24
       isection=24
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_eigv).NE.0) THEN
          irecord=nkpt%nkpts+1
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) ener_com%amu
          DO ik=1,nkpt%nkpts
             WRITE (unit=nw,iostat=ierror) (eigv(i,ik),i=1,nstate)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 25
       isection=25
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_clas).NE.0) THEN
          irecord=clas3%nclatom
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO i=1,clas3%nclatom
             WRITE (unit=nw,iostat=ierror)&
                  (clasc(j,i),j=1,3),(clasv(j,i),j=1,3)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 26
    isection=26
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_lrwf).NE.0) THEN
       CALL w_linres (nw,ierror,"WF",irecord)
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 27
    isection=27
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_lrwv).NE.0) THEN
       CALL w_linres (nw,ierror,"WV",irecord)
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    IF (paral%io_parent) THEN
       ! ==------------------------------------------------------------==
       ! Section 28
       isection=28
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_phes).NE.0.AND..NOT.lallhs) THEN
          irecord=nvar+1
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          DO i=1,nvar
             WRITE (unit=nw,iostat=ierror) (hesscr(j,i),j=1,nvar)
          ENDDO
          WRITE (unit=nw,iostat=ierror) nshess, lvlhes
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 29
       isection=29
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_prfo).NE.0) THEN
          irecord=8
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          CALL wr30prfo (nw, ierror,irec(irec_lbfgs))
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 30
       isection=30
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_lbfgs).NE.0) THEN
          lsavc0 = .FALSE. ! Do not save previous wavefunction
          IF (lsavc0) THEN
             irecord=7
          ELSE
             irecord=7
          ENDIF
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          CALL wr30lbfgs (nw, ierror,lsavc0)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 31
       isection=31
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_rdtl).NE.0) THEN
          irecord=1
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               adtolr%delast, adtolr%demin, adtolr%tolmin, cntr%tolog, cntr%tolad, cntr%tolene
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 32
       isection=32
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_cons).NE.0 .AND. cotc0%mcnstr.GT.0) THEN
          irecord=2
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) cotc0%mcnstr
          WRITE (unit=nw,iostat=ierror) (cnsval(i),i=1,cotc0%mcnstr)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 33
       isection=33
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_ctrans).NE.0) THEN
          irecord=1
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) (clsaabox%mm_c_trans(i),i=1,3)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 34
       ! FIXME LOCT : WRITE GROUP SELECTIONS AND PARAMETERS TO FILE
       isection=34
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_nop4).NE.0) THEN
          irecord=2
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               ((lcttab(j,i),i=1,maxsys%nax),j=1,maxsys%nsx)
          WRITE (unit=nw,iostat=ierror)&
               ((etap(j,i,ip),etapdot(j,i,ip),i=1,nchx),j=1,maxsys%nax*maxsys%nsx)
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 35
       isection=35
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_shop).NE.0) THEN
          irecord=5
          WRITE(unit=nw,iostat=ierror) irecord,fpos_null
          WRITE(unit=nw,iostat=ierror) sh03%isurf,(sh03%pop(i),i=1,6),sh03%det,sh03%detold
          DO i=1,2
             WRITE(unit=nw,iostat=ierror) (sh03%coupl(i,k),sh03%couplold(i,k),k=1,2)
             WRITE(unit=nw,iostat=ierror) sh03%ec(i),sh03%eold(i)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    ! ==------------------------------------------------------------==
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 36
    isection=36
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF
    IF (irec(irec_shopbo).NE.0) THEN
       CALL w_wfnio(nw,ierror,nstate,c0old,tau0,'C0',irecord)
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 37
    isection=37
    prngin=0.0_real_8
    prngout=0.0_real_8
    IF (cntl%tpath.OR.tmw) THEN
       prngin(parai%mepos+1,:,:)=prng_com%paraseed
    ELSE
       prngin(parai%me+1,:,:)=prng_com%paraseed
    ENDIF
    !msglen = parai%nproc * 6 * 8
    CALL mp_sum(prngin,prngout,parai%nproc*6,parai%allgrp)
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_prng).NE.0) THEN
          irecord=1+parai%nproc
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror)&
               (prng_com%repseed(1,i),i=1,3),(prng_com%repseed(2,i),i=1,3)
          DO j=1,parai%nproc
             WRITE (unit=nw,iostat=ierror)&
                  (prngout(j,1,i),i=1,3),(prngout(j,2,i),i=1,3)
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 38
    isection=38
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION ',isectioN
       ENDIF
       IF (irec(irec_gle).NE.0) THEN
          irecord=1+ions1%nat*(glepar%gle_ns+1)
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          WRITE (unit=nw,iostat=ierror) glepar%egle
          DO is = 1 , ions1%nsp
             DO ia = 1, ions0%na(is)
                DO i= 1, glepar%gle_ns+1
                   WRITE (unit=nw,iostat=ierror) (glep(j,ia,is,i),j=1,3)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
       IF (cntl%is_out_stream) THEN
          CALL file_get_position(nw,.FALSE.,fpos_end)
          CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
       ENDIF
    ENDIF
    IF (ierror.NE.0) CALL stopgm('WV30','ERROR IN WRITING RESTART'//&
         ' FILE section='//TRIM(int2str(isection)),& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ! Section 39- 98
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'WRITE SECTION 38-98'
       ENDIF
       DO i=39,98
          isection=i
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
          IF (cntl%is_out_stream) THEN
             CALL file_get_position(nw,.FALSE.,fpos_end)
             CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 99
    isection=99
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'WRITE SECTION ',isectioN
    ENDIF

    IF (irec(irec_xtrp).NE.0) THEN
       ! NUMCOLD is the last 'official' record.
       ! the wavefunction history is then added.
       irecord=1
       IF (paral%io_parent) THEN
          WRITE(unit=nw,iostat=ierror) irecord,fpos_null
          WRITE(unit=nw,iostat=ierror) numcold
       ENDIF
       DO i=numcold,1,-1
          is=MOD(nnow+cnti%mextra-i,cnti%mextra)+1
          IF (store1%tdebio.AND.paral%io_parent) THEN
             WRITE(6,*) 'STORING WFN ',is,' OF ',numcolD
          ENDIF
          CALL w_wfnio(nw,ierror,nstate,cold(1,1,1,is),tau0,&
               'C0',irecord)
       ENDDO
    ELSE
       IF (paral%io_parent) THEN
          irecord=irec0
          WRITE (unit=nw,iostat=ierror) irecord,fpos_null
       ENDIF
    ENDIF
    IF (paral%io_parent.AND.cntl%is_out_stream) THEN
       CALL file_get_position(nw,.FALSE.,fpos_end)
       CALL file_update_record(nw,irecord,fpos_beg,fpos_end)
    ENDIF
    ! ==------------------------------------------------------------==
    RETURN
  END SUBROUTINE wv30
  SUBROUTINE file_update_record(u,irecord,fpos_beg,fpos_end)
    INTEGER                                  :: u, irecord
    INTEGER(int_8)                           :: fpos_beg, fpos_end

    CHARACTER(*), PARAMETER :: procedureN = 'file_update_record'

#if ! defined(_HASNT_BF_STREAM_IO)
    WRITE (u,pos=fpos_beg) irecord,fpos_end
    WRITE (u,pos=fpos_end)
#else
    CALL stopgm(procedureN,'we shouldnt be here',& 
         __LINE__,__FILE__)
#endif
    fpos_beg=fpos_end
    RETURN
  END SUBROUTINE file_update_record
  ! ==================================================================

END MODULE wv30_utils


SUBROUTINE wr30pot(nw,ierror,potr)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send,mp_recv
  USE system , ONLY: fpar,spar, parap
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror
  REAL(real_8)                               :: potr(fpar%kr1,*)

  CHARACTER(*), PARAMETER                    :: procedureN = 'wr30pot'

  INTEGER                                    :: i, ierr, ip, ipp, ipx, ir, &
                                                irr, kk, msgid
  REAL(real_8), ALLOCATABLE                  :: pscr(:)

! ==--------------------------------------------------------------==

  kk=fpar%kr2s*fpar%kr3s
  ALLOCATE(pscr(kk),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL mp_sync(parai%allgrp)
  ! For each plane
  DO ir=1,spar%nr1s
     ! Find the processor which has the given plane.
     DO ipp=1,parai%nproc
        IF (ir.GE.parap%nrxpl(ipp-1,1).AND.ir.LE.parap%nrxpl(ipp-1,2)) THEN
           ipx=ipp
           GOTO 10
        ENDIF
     ENDDO
     CALL stopgm('WR30POT','ERROR IN NRXPL',& 
          __LINE__,__FILE__)
10   CONTINUE
     ip=parap%pgroup(ipx)
     ! IP has the plane.
     IF (parai%me.EQ.ip) THEN
        irr=ir-parap%nrxpl(parai%mepos,1)+1
        CALL dcopy(kk,potr(irr,1),fpar%kr1,pscr(1),1)
        IF (.NOT.paral%parent) THEN
           ! IP sends array.
           msgid=1
           !msglen=kk * 8
           CALL mp_send(pscr,kk,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
     ELSEIF (paral%parent) THEN
        ! PARENT receives array.
        msgid=1
        !msglen=kk * 8
        CALL mp_recv(pscr,kk,ip,msgid,parai%allgrp)
     ENDIF
     IF (paral%parent) THEN
        IF (paral%io_parent)&
             WRITE(nw,iostat=ierror) (pscr(i),i=1,kk)
     ENDIF
     CALL mp_sync(parai%allgrp)
  ENDDO
  DEALLOCATE(pscr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wr30pot
! ==================================================================
