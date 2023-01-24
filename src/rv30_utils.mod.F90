MODULE rv30_utils
  USE cdftmod,                         ONLY: cdfthda,&
                                             cdftlog,&
                                             projlog
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
  USE filnmod,                         ONLY: filn
  USE geq0mod,                         ONLY: geq0
  USE glemod,                          ONLY: glep,&
                                             glepar
  USE io_utils,                        ONLY: file_open,&
                                             file_seek
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE lscal,                           ONLY: adtolr,&
                                             hesscr,&
                                             lallhs,&
                                             lundop,&
                                             lvlhes,&
                                             nshess,&
                                             nvar
  USE machine,                         ONLY: m_flush
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mimic_wrapper,                   ONLY: mimic_save_dim,&
                                             mimic_switch_dim,&
                                             mimic_revert_dim
  USE mm_dimmod,                       ONLY: clsaabox
  USE mm_extrap,                       ONLY: cold,&
                                             nnow,&
                                             numcold
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nose,                            ONLY: &
       eta, etadot, etap, etap1, etap1dot, etapdot, etapm, etapmdot, etc, &
       etcdot, lcttab, nchx
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: ipcurr
  USE poin,                            ONLY: cstate,&
                                             potr,&
                                             rhoo
  USE prng,                            ONLY: prng_com
  USE prng_utils,                      ONLY: prngparaskip
  USE readsr_utils,                    ONLY: xstring
  USE rlbfgs_io,                       ONLY: rd30lbfgs,&
                                             rd30prfo
  USE rw_linres_utils,                 ONLY: r_linres
  USE shop_rest,                       ONLY: sh03
  USE shop_rest_2,                     ONLY: c0old
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE store_types,                     ONLY: &
       irec_ac, irec_cell, irec_clas, irec_co, irec_cons, irec_ctrans, &
       irec_cut, irec_eigv, irec_gle, irec_he, irec_ico, irec_kpt, &
       irec_lbfgs, irec_lrwf, irec_lrwv, irec_nas, irec_noc, irec_noe, &
       irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_nose, irec_occ, &
       irec_phes, irec_pot, irec_prfo, irec_prng, irec_pws, irec_pwv, &
       irec_rdtl, irec_rho, irec_sdp, irec_shop, irec_shopbo, irec_vel, &
       irec_wf, irec_xtrp, restart1, store1
  USE string_utils,                    ONLY: int2str
  USE symm,                            ONLY: symmi
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, dual00, fpar, kpbeg, maxsp, maxsys, nacc, ncpw, &
       nkpbl, nkpt, parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: fskip
  USE wfnio_utils,                     ONLY: r_wfnio
  use bicanonicalCpmd, only: biCanonicalEnsembleDo,&
    bicanonicalCpmdConfig, getNameLatestTape
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: zhrwf
  PUBLIC :: rv30
  !public :: rd30pot
  !public :: rd30kpnt
  PUBLIC :: rdkpoints

CONTAINS

  ! ==================================================================
  SUBROUTINE zhrwf(nr,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)
    ! ==--------------------------------------------------------------==
    ! == Read wavefunction in Version 3.0 format on unit NR           ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, irec(:)
    COMPLEX(real_8)                          :: c0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(nkpt%ngwk,nstate)
    REAL(real_8)                             :: eigv(nstate,*), tau0(:,:,:), &
                                                taum(:,:,:), taui(:,:,:)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'zhrwf'

    CHARACTER(len=20)                        :: fformat, latfile
    CHARACTER(len=fo_fpath_max)              :: filnam
    INTEGER                                  :: ia, ie, iostat, isub
    LOGICAL                                  :: fexist

    CALL tiset(proceduren,isub)
    IF (cntl%mimic) THEN
       CALL mimic_save_dim()
       CALL mimic_switch_dim(go_qm=.FALSE.)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Construct the filename
    IF (paral%io_parent) THEN
       CALL xstring(filn,ia,ie)
       IF (cntl%tpath) THEN
          filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
       ELSE
          IF (restart1%rlate) THEN
             IF (tmw)THEN
                CALL mw_filename('LATEST_',latfile,mwi%walker_id)
             elseif (biCanonicalEnsembleDo) then
               latfile = getNameLatestTape(bicanonicalCpmdConfig)
               !JF call bicanonicalSetFileNameLatest(latfile)
             ELSE
                latfile='LATEST'
             ENDIF
             INQUIRE(file=latfile,exist=fexist)
             IF (fexist) THEN
                OPEN(unit=nr,file=latfile,status='UNKNOWN')
                READ(nr,'(A)') filnam
                CLOSE(nr)
                filnam=ADJUSTL(filnam)
             ELSE
                filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
             ENDIF
          ELSE
             filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
          ENDIF
       ENDIF
       IF (cdftlog%rcdft.AND.cdftlog%thda.AND..NOT.projlog%projnow)THEN
          IF (cdfthda%hdafirst)THEN
             filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//"RESTART.STATE_1"
          ELSE
             filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//"RESTART.STATE_2"
          ENDIF
       ENDIF
       INQUIRE(file=filnam,exist=fexist)
       IF (.NOT.fexist) THEN
          CALL xstring(filnam,ia,ie)
          WRITE(fformat,'(A,I2,A)')&
               '(/,A,T',MAX(35,65-(ie-ia)),',A)'
          WRITE(6,fformat) ' ZHRWF| RESTART FILE NOT FOUND:',filnam
          CALL stopgm('ZHRWF',' FILE NOT FOUND ',& 
               __LINE__,__FILE__)
       ENDIF

       CALL file_open(filnam,'UNFORMATTED','OLD','READ',&
            cntl%is_in_stream,.FALSE.,0,nr)

    ENDIF

    CALL rv30(nr,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)

    IF (paral%io_parent) THEN
       CLOSE(nr,iostat=iostat)
       IF (iostat /= 0 ) CALL stopgm(proceduren,&
            'Problem closing FILE = '//TRIM(filnam)//&
            ' IOSTAT = '//TRIM(int2str(iostat)),& 
            __LINE__,__FILE__)
    ENDIF

    IF (paral%io_parent) THEN
       CALL xstring(filnam,ia,ie)
       IF (cntl%tpath) THEN
          WRITE(fformat,'(A,I2,A)') '(A,T',MAX(35,65-(ie-ia)),',A)'
       ELSE
          WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(35,65-(ie-ia)),',A)'
       ENDIF
       WRITE(6,fformat)&
            ' RESTART INFORMATION READ ON FILE ',filnam(ia:ie)
       CLOSE(nr)
    ENDIF
    IF (cntl%mimic) THEN
       CALL mimic_revert_dim()
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(proceduren,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zhrwf
  ! ==================================================================
  SUBROUTINE rv30(nr,irec,c0,cm,nstate,eigv,tau0,taum,taui,nfi)
    ! ==--------------------------------------------------------------==
    ! == IREC: input -> flag array to indicate the read of record     ==
    ! ==       output-> 0 if not read or no information               ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, irec(:), nstate
    COMPLEX(real_8) :: cm(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate*nkpt%nkpnt)
    REAL(real_8)                             :: eigv(nstate,nkpt%nkpts)
    REAL(real_8), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: tau0, taum, taui
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'rv30'
    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8 

    CHARACTER(len=80)                        :: str
    INTEGER :: bcast(4:17), i, ia, ibrav0, ierr, ig, ik, indpg0, info, ip, &
      irec0, irecord, is, j, k, l, mcnstr0, n0, na0(maxsp), nacc0, nacc1, &
      nel0, ngwks0, ngwls0, nhgls0, nhgs0, nkmin, nkpts0, nmin, nr1s0, nr2s0, &
      nr3s0, nsp0
    INTEGER(int_8)                           :: fpos
    LOGICAL                                  :: dual0, tkpnt0, tlsd0, tlse0
    REAL(real_8)                             :: buff(3), cdual0, cenew(6), &
                                                ecut0, &
                                                prngin(parai%nproc,2,3), &
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
! ==  Section 24: Fermi energy and eigenvalues                    ==
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
! ==  Section 99: Wavefunction history for extrapolation          ==
! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

    bcast(:)=0
    IF (cntl%tpath.AND.cntl%tpimd) THEN
       ip=ipcurr
    ELSEIF (tmw)THEN
       ip=mwi%walker_id
    ELSE
       ip=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) REWIND (nr)
    IF (paral%io_parent) THEN
       ! Section 1
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 1'
       ENDIF
       READ(nr) str
       IF (INDEX(str,'STUTTGART VERSION 3.0').EQ.0) THEN
          CALL m_flush(6)
          WRITE(6,*) ' WRONG FILE FORMAT'
          WRITE(6,'(A)') str
          CALL stopgm(proceduren,'HEADER',& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%is_in_stream) THEN
          IF (INDEX(str,'STREAM').EQ.0) THEN
             WRITE(6,*) ' WRONG FILE FORMAT ',str
             CALL stopgm(proceduren,'The restart file'&
                  //' is not in STREAM format',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 2
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 2'
       ENDIF
       IF (irec(irec_cell).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.2) THEN
             READ(nr) ibrav0,indpg0
             READ(nr) (cenew(i),i=1,6)
             IF (ibrav0.NE.parm%ibrav)&
                  WRITE(6,'(A,T56,2I5)')&
                  ' RV30! SUPERCELL SYMMETRY   HAS CHANGED ',&
                  ibrav0,parm%ibrav
             IF (indpg0.NE.symmi%indpg)&
                  WRITE(6,'(A,T56,2I5)')&
                  ' RV30! POINT GROUP SYMMETRY HAS CHANGED ',&
                  indpg0,symmi%indpg
             DO i=1,6
                IF (ABS(cenew(i)-cell_com%celldm(i)).GT.small)&
                     WRITE(6,'(A,T65,I1)')&
                     ' RV30! SUPERCELL DIMENSIONS HAVE CHANGED: NUMBER=',I
             ENDDO
          ELSE
             WRITE(6,'(A)') ' RV30| WARNING! CANNOT READ SECTION 2 '
             irec(irec_cell)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 3
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 3'
       ENDIF
       IF (irec(irec_nas).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.2) THEN
             READ(nr) nsp0
             READ(nr) (na0(i),i=1,nsp0)
             IF (nsp0.NE.ions1%nsp)&
                  WRITE(6,'(A,T56,2I5)')&
                  ' RV30! NUMBER OF SPECIES HAS CHANGED ',&
                  nsp0,ions1%nsp
             DO is=1,MIN(nsp0,ions1%nsp)
                IF (ions0%na(is).NE.na0(is))&
                     WRITE(6,'(A,I3,A,T56,2I5)')&
                     ' RV30! NUMBER OF ATOMS FOR SPECIES ',IS,&
                     ' HAS CHANGED ',na0(is),ions0%na(is)
             ENDDO
          ELSE
             WRITE(6,'(A)') ' RV30| WARNING! CANNOT READ SECTION 3 '
             irec(irec_nas)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 4
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 4'
       ENDIF
       IF (irec(irec_co).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ( (irec(irec_nas).EQ.0).OR.(irecord.NE.ions1%nat) ) THEN
             WRITE(6,'(A)') ' RV30! INCONSISTENT NUMBER OF ATOMS.'
             WRITE(6,'(A)')&
                  ' RV30| WARNING! CANNOT READ COORDINATES FROM SECTION 4'
             irec(irec_co)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             bcast(4) = 1
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   READ(nr) (tau0(j,ia,is),j=1,3)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 5
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 5'
       ENDIF
       IF (irec(irec_vel).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ( (irec(irec_nas).EQ.0).OR.(irecord.NE.ions1%nat) ) THEN
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO IONIC VELOCITIES'
             irec(irec_vel)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             bcast(5) = 1
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   READ(nr) (taum(j,ia,is),j=1,3)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 6
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 6'
       ENDIF
       IF (irec(irec_ico).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ( (irec(irec_nas).EQ.0).OR.(irecord.NE.ions1%nat) ) THEN
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO INITIAL GEOMETRY'
             irec(irec_ico)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             bcast(6) = 1
             DO is = 1 , ions1%nsp
                DO ia = 1 , ions0%na ( is )
                   READ(nr) (taui(j,ia,is),j=1,3)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 7
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 7'
       ENDIF
       IF (irec(irec_cut).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             READ(nr) ecut0,cdual0,dual0,nel0,nr1s0,nr2s0,nr3s0
             IF (ABS(ecut0-cntr%ecut).GT.small) WRITE(6,'(A,T46,2F10.2)')&
                  ' RV30! CUTOFF HAS CHANGED ',ECUT0,cntr%ecut
             IF (ABS(cdual0-dual00%cdual).GT.small) WRITE(6,'(A,T46,2F10.2)')&
                  ' RV30! DUAL HAS CHANGED ',CDUAL0,dual00%cdual
             IF (nel0.NE.NINT(crge%nel)) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! NUMBER OF ELECTRONS HAS CHANGED ',&
                  nel0,NINT(crge%nel)
             IF (nr1s0.NE.spar%nr1s) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! X REAL SPACE MESH  HAS CHANGED ',nr1s0,spar%nr1s
             IF (nr2s0.NE.spar%nr2s) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! Y REAL SPACE MESH  HAS CHANGED ',nr2s0,spar%nr2s
             IF (nr3s0.NE.spar%nr3s) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! Z REAL SPACE MESH  HAS CHANGED ',nr3s0,spar%nr3s
          ELSE
             WRITE(6,'(A)') ' RV30| WARNING! CANNOT READ SECTION 7 '
             irec(irec_cut)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 8
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 8'
       ENDIF
       IF (irec(irec_sdp).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             READ(nr) n0,nkpts0,ngwks0,ngwls0,nhgs0,nhgls0
             IF (n0.NE.crge%n) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! NUMBER OF STATES   HAS CHANGED ',n0,crge%N
             IF (tkpts%tkpnt) THEN
                IF (nkpts0.EQ.0) THEN
                   ! WRITE(6,'(A)') ' RV30! TYPE   OF K-POINTS HAS CHANGED '
                ELSE IF (nkpts0.NE.nkpt%nkpts) THEN
                   WRITE(6,'(A,T56,2I5)')&
                        ' RV30! NUMBER OF K-POINTS HAS CHANGED ',&
                        nkpts0,nkpt%nkpts
                ENDIF
             ELSE IF (nkpts0.NE.0) THEN
                ! WRITE(6,'(A)') ' RV30! TYPE   OF K-POINTS HAS CHANGED '
             ENDIF
             IF (nkpts0.EQ.0) THEN
                tkpnt0=.FALSE.
                nkpts0=1
             ELSE
                tkpnt0=.TRUE.
             ENDIF
          ELSE
             WRITE(6,'(A)') ' RV30| WARNING! CANNOT READ SECTION 8 '
             irec(irec_sdp)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
             tkpnt0=.FALSE.
             nkpts0=1
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
          irec(irec_sdp)=0
          tkpnt0=.FALSE.
          nkpts0=1
       ENDIF
    ENDIF
    CALL mp_bcast(n0,parai%io_source,parai%cp_grp)
    CALL mp_bcast(nkpts0,parai%io_source,parai%cp_grp)
    CALL mp_bcast(tkpnt0,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! Section 9
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 9'
    ENDIF
    IF (irec(irec_wf).NE.0) THEN
       IF (irec(irec_sdp).NE.0.OR.projlog%projnow.OR.cntl%tsyscomb) THEN
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          IF (cntl%tsyscomb)THEN
             n0=nstate
             nkpts0=1
          ENDIF
          IF (paral%io_parent) THEN
             IF (cntl%is_in_stream)THEN
                READ(nr) irecord,fpos
             ELSE
                READ(nr) irecord
             ENDIF
          ENDIF
          CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
          IF (irecord.EQ.0) THEN
             info=1
          ELSEIF (irecord.LT.0) THEN
             CALL r_wfnio(nr,c0,nstate,info,tkpnt0,n0,nkpts0,&
                  'C0',fpos)
          ELSE
             CALL rd30wfn(nr,c0,nstate,tau0,info,&
                  tkpnt0,n0,nkpts0,'C0',fpos)
          ENDIF
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       ELSE
          info=1
       ENDIF
       IF (info.NE.0) THEN
          IF (paral%io_parent) WRITE(6,'(A)')&
               ' RV30| WARNING! NO INITIAL WAVEFUNCTIONS'
          irec(irec_wf)=0
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF

    ! CALL FILE_SEEK(NR,.FALSE.,FPOS)

    ! ==--------------------------------------------------------------==
    ! Section 10
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 10'
    ENDIF
    IF (irec(irec_pwv).NE.0) THEN
       IF (irec(irec_sdp).NE.0) THEN
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          IF (paral%io_parent) THEN
             IF (cntl%is_in_stream)THEN
                READ(nr) irecord,fpos
             ELSE
                READ(nr) irecord
             ENDIF
          ENDIF
          CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
          IF (irecord.EQ.0) THEN
             info=1
          ELSEIF (irecord.LT.0) THEN
             CALL r_wfnio(nr,cm,nstate,info,tkpnt0,n0,nkpts0,&
                  'NIL',fpos)
          ELSE
             CALL rd30wfn(nr,cm,nstate,tau0,info,&
                  tkpnt0,n0,nkpts0,'NIL',fpos)
          ENDIF
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       ELSE
          info=1
       ENDIF
       IF (info.NE.0) THEN
          IF (paral%io_parent) WRITE(6,'(A)')&
               ' RV30| WARNING! NO WAVEFUNCTION VELOCITIES '//&
               '(SET TO ZERO)'
          irec(irec_pwv)=0
          CALL zeroing(cm)!,nkpt%ngwk*nstate)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 11
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 11'
       ENDIF
       IF (irec(irec_ac).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.2) THEN
             WRITE(6,'(A,2I5)')&
                  ' RV30| WARNING! NO ACCUMULATOR INFORMATION'
             irec(irec_ac)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             READ(nr) nfi,nacc0
             IF (nacc0.NE.nacc) WRITE(6,'(A,T56,2I5)')&
                  ' RV30! NUMBER OF ACCUMULATOR ENTRIES HAS CHANGED ',&
                  nacc0,nacc
             nacc1=MIN(nacc0,nacc)
             READ(nr) ( acc ( i ) , i = 1 , nacc1 )
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 12
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 12'
       ENDIF
       IF (irec(irec_nose).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.2) THEN
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 12'
             irec(irec_nose)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             ! ..Nose parameters are fixed by input only!
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                READ(nr)
                READ(nr)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 13
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 13'
       ENDIF
       IF (irec(irec_noe).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ((irecord.EQ.2).AND.cntl%bsymm) THEN
             bcast(13) = 1
             READ(nr) (eta(i,1),etadot(i,1),i=1,nchx)
             READ(nr) (eta(i,2),etadot(i,2),i=1,nchx)
          ELSEIF ((irecord.EQ.2).AND.(.NOT.cntl%bsymm)) THEN
             bcast(13) = 1
             READ(nr) (eta(i,ip),etadot(i,ip),i=1,nchx)
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSEIF (irecord.EQ.1) THEN
             bcast(13) = 1
             READ(nr) (eta(i,ip),etadot(i,ip),i=1,nchx)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 13'
             irec(irec_noe)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 14
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 14'
       ENDIF
       IF (irec(irec_nop1).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             bcast(14) = 1
             READ(nr) (etap1(i,ip),etap1dot(i,ip),i=1,nchx)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 14'
             irec(irec_nop1)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 15
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 15'
       ENDIF
       IF (irec(irec_nop2).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             bcast(15) = 1
             READ(nr) ((etap(j,i,ip),etapdot(j,i,ip),&
                  i=1,nchx),j=1,maxsys%nsx)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 15'
             irec(irec_nop2)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! Section 16
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 16'
       ENDIF
       IF (irec(irec_nop3).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             bcast(16) = 1
             READ(nr) ((etapm(k,l,ip),etapmdot(k,l,ip),&
                  k=1,3*maxsys%nax*maxsys%nsx),&
                  l=1,nchx)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 16'
             irec(irec_nop3)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
       ! ==------------------------------------------------------------==
       ! Section 17
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 17'
       ENDIF
       IF (irec(irec_noc).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             bcast(17) = 1
             READ(nr) (etc(i,ip),etcdot(i,ip),i=1,nchx)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 17'
             irec(irec_noc)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF

    ! All this garbage should be bcast
    CALL mp_bcast(bcast,SIZE(bcast),parai%io_source,parai%cp_grp)

    IF (bcast(4)==1) THEN
       CALL mp_bcast(tau0,SIZE(tau0),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(5)==1) THEN
       CALL mp_bcast(taum,SIZE(taum),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(6)==1) THEN
       !vw bcasting taui causes segfault with gcc and intel compilers
       !CALL mp_bcast(taui,SIZE(taui),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(13)==1) THEN
       CALL mp_bcast(eta,SIZE(eta),parai%io_source,parai%cp_grp)
       CALL mp_bcast(etadot,SIZE(etadot),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(14)==1) THEN
       CALL mp_bcast(etap1,SIZE(etap1),parai%io_source,parai%cp_grp)
       CALL mp_bcast(etap1dot,SIZE(etap1dot),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(15)==1) THEN
       CALL mp_bcast(etap,SIZE(etap),parai%io_source,parai%cp_grp)
       CALL mp_bcast(etapdot,SIZE(etapdot),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(16)==1) THEN
       CALL mp_bcast(etapm,SIZE(etapm),parai%io_source,parai%cp_grp)
       CALL mp_bcast(etapmdot,SIZE(etapmdot),parai%io_source,parai%cp_grp)
    ENDIF
    IF (bcast(17)==1) THEN
       CALL mp_bcast(etc,SIZE(etc),parai%io_source,parai%cp_grp)
       CALL mp_bcast(etcdot,SIZE(etcdot),parai%io_source,parai%cp_grp)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! We communicate some information for the other processors
    CALL mp_bcast(nr1s0,parai%io_source,parai%cp_grp)
    CALL mp_bcast(nr2s0,parai%io_source,parai%cp_grp)
    CALL mp_bcast(nr3s0,parai%io_source,parai%cp_grp)
    CALL mp_bcast(nfi,parai%io_source,parai%cp_grp)
    ! Section 18
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 18'
    ENDIF
    IF (irec(irec_pot).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ (nr) irecord,fpos
          ELSE
             READ (nr) irecord
          ENDIF
       ENDIF
       CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
       IF ( (irecord.NE.0).AND.(irec(irec_cut).NE.0) ) THEN
          tlsd0=irecord.EQ.2*nr1s0
          IF (tlsd0)irecord=irecord/2
          tlse0=irecord.EQ.4*nr1s0
          IF (tlse0)irecord=irecord/4
          CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr(1,1))
          IF (cntl%tlsd) THEN
             IF (tlsd0) THEN
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr(1,2))
             ELSE
                CALL dcopy(fpar%nnr1,potr(1,1),1,potr(1,2),1)
             ENDIF
          ELSEIF (lspin2%tlse) THEN
             IF (tlse0) THEN
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr(1,2))
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr(1,3))
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr(1,4))
             ELSE
                CALL dcopy(fpar%nnr1,potr(1,1),1,potr(1,2),1)
                CALL dcopy(fpar%nnr1,potr(1,1),1,potr(1,3),1)
                CALL dcopy(fpar%nnr1,potr(1,1),1,potr(1,4),1)
             ENDIF
          ELSE
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                IF (tlsd0) CALL fskip(nr,irecord)
                IF (tlse0) CALL fskip(nr,irecord*3)
             ENDIF
          ENDIF
       ELSEIF (paral%io_parent) THEN
          WRITE(6,'(A)')&
               ' RV30| WARNING! NO POTENTIAL INFORMATION'
          irec(irec_pot)=0
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irecord)
          ENDIF
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 19
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 19'
    ENDIF
    IF (irec(irec_pws).NE.0) THEN
       IF (irec(irec_sdp).NE.0) THEN
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          IF (paral%io_parent) THEN
             IF (cntl%is_in_stream)THEN
                READ(nr) irecord,fpos
             ELSE
                READ(nr) irecord
             ENDIF
          ENDIF
          CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
          IF (irecord.EQ.0) THEN
             info=1
          ELSEIF (irecord.LT.0) THEN
             CALL r_wfnio(nr,cstate,1,info,tkpnt0,n0,nkpts0,&
                  'NIL',fpos)
          ELSE
             CALL rd30wfn(nr,cstate,1,tau0,info,&
                  tkpnt0,n0,nkpts0,'NIL',fpos)
          ENDIF
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       ELSE
          info=1
       ENDIF
       IF (info.NE.0.AND.paral%io_parent) THEN
          WRITE(6,'(A)')&
               ' RV30| WARNING! CANNOT READ SECTION 19'
          irec(irec_pws)=0
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 20
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 20'
    ENDIF
    IF (paral%io_parent) THEN
       IF (irec(irec_he).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.6) THEN
             WRITE(6,'(A)')&
                  ' RV30| WARNING! CANNOT READ CELL INFORMATION'
             irec(irec_he)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             READ(nr) metr_com%ht
             READ(nr) metr_com%htm1
             READ(nr) metr_com%htvel
             READ(nr) metr_com%htfor
             READ(nr) metr_com%ht0
             READ(nr) metr_com%htm10
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 21
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 21'
    ENDIF
    IF (irec(irec_kpt).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ (nr) irecord,fpos
          ELSE
             READ (nr) irecord
          ENDIF
       ENDIF
       CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
       IF (irecord.NE.0) THEN
          IF (tkpts%tkpnt.AND.(irecord.EQ.nkpts0).AND.&
               (irecord.NE.0).AND.(irec(irec_sdp).NE.0)) THEN
             CALL rd30kpnt(nr,nstate,nkpts0,tkpnt0,c0,cm,tau0,irec)
          ELSE
             IF ((irecord.NE.nkpts0).OR.(irec(irec_sdp).EQ.0)) THEN
                IF (paral%io_parent) THEN
                   WRITE(6,'(A)')&
                        ' RV30| WARNING! NO K-POINTS INFORMATION'
                   irec(irec_he)=0
                ENDIF
             ELSE
                IF (tkpts%tkpnt.AND.nkpt%nkpts.EQ.1.AND.(rk(1,1).EQ.0._real_8)&
                     .AND.(rk(2,1).EQ.0._real_8).AND.(rk(3,1).EQ.0._real_8)) THEN
                   DO i=1,nstate
                      DO ig=1,ncpw%ngw
                         c0(ig+ncpw%ngw,i)=CONJG(c0(ig,i))
                      ENDDO
                      IF (geq0) c0(ncpw%ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
                   ENDDO
                ELSE
                   IF (cntl%tlsd) THEN
                      CALL putconj(c0(1,1),spin_mod%nsup)
                      CALL putconj(c0(1,spin_mod%nsup+1),spin_mod%nsdown)
                   ELSE
                      CALL putconj(c0,nstate)
                   ENDIF
                ENDIF
                IF (paral%io_parent)&
                     WRITE(6,'(A)')&
                     ' RV30! TYPE   OF K-POINTS HAS CHANGED'
             ENDIF
             IF (paral%io_parent) THEN
                IF (cntl%is_in_stream) THEN
                   CALL file_seek(nr,.FALSE.,fpos)
                ELSE
                   CALL fskip(nr,irecord)
                ENDIF
             ENDIF
          ENDIF
       ELSE
          IF (tkpts%tkpnt) THEN
             DO i=1,nstate
                DO ig=1,ncpw%ngw
                   c0(ig+ncpw%ngw,i)=CONJG(c0(ig,i))
                ENDDO
                IF (geq0) c0(ncpw%ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,'(A)')&
                  ' RV30! TYPE   OF K-POINTS HAS CHANGED'
          ENDIF
          IF (paral%io_parent) THEN
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ENDIF
       ! EHR[
    ELSEIF (cntl%cmplx_wf.AND.cntl%start_real) THEN
       DO i=1,nstate
          DO ig=1,ncpw%ngw
             c0(ig+ncpw%ngw,i)=CONJG(c0(ig,i))
          ENDDO
          IF (geq0) c0(ncpw%ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            ' RV30! TYPE   OF K-POINTS HAS CHANGED 3'
       ! EHR]
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
          IF (( (tkpts%tkpnt.AND.(.NOT.tkpnt0)).OR.(.NOT.tkpts%tkpnt.AND.tkpnt0) )&
               )WRITE(6,'(A)') ' RV30! TYPE   OF K-POINTS HAS CHANGED'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 22
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 22'
    ENDIF
    IF (irec(irec_rho).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ (nr) irecord,fpos
          ELSE
             READ (nr) irecord
          ENDIF
       ENDIF
       CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
       IF ( (irecord.NE.0).AND.(irec(irec_cut).NE.0) ) THEN
          tlsd0=irecord.EQ.2*nr1s0
          IF (tlsd0)irecord=irecord/2
          tlse0=irecord.EQ.4*nr1s0
          IF (tlse0)irecord=irecord/4
          ! For LSD RHOO(1,1)=alpha+beta density, RHOO(1,2)=beta density
          CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,rhoo(1,1))
          IF (cntl%tlsd) THEN
             IF (tlsd0) THEN
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,rhoo(1,2))
             ELSE
                CALL zeroing(rhoo(:,2))!,nnr1)
                CALL daxpy(fpar%nnr1,0.5_real_8,rhoo(1,1),1,rhoo(1,2),1)
             ENDIF
             ! For LSE RHOO(1,1)=total density, RHOO(1,2)=beta density, m state
             ! For LSE RHOO(1,3)=total density, RHOO(1,4)=beta density, t state
          ELSEIF (lspin2%tlse) THEN
             IF (tlse0) THEN
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,rhoo(1,2))
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,rhoo(1,3))
                CALL rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,rhoo(1,4))
             ELSE
                CALL zeroing(rhoo(:,2))!,nnr1)
                CALL daxpy(fpar%nnr1,0.5_real_8,rhoo(1,1),1,rhoo(1,2),1)
                CALL zeroing(rhoo(:,3))!,nnr1)
                CALL dcopy(fpar%nnr1,rhoo(1,1),1,rhoo(1,3),1)
                CALL zeroing(rhoo(:,4))!,nnr1)
                CALL daxpy(fpar%nnr1,0.5_real_8,rhoo(1,3),1,rhoo(1,4),1)
             ENDIF
          ELSE
             IF (paral%io_parent) THEN
                IF (cntl%is_in_stream) THEN
                   CALL file_seek(nr,.FALSE.,fpos)
                ELSE
                   IF (tlsd0) CALL fskip(nr,irecord)
                   IF (tlse0) CALL fskip(nr,irecord*3)
                ENDIF
             ENDIF
          ENDIF
       ELSEIF (paral%io_parent) THEN
          WRITE(6,'(A)')&
               ' RV30| WARNING! NO DENSITY INFORMATION'
          irec(irec_rho)=0
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irecord)
          ENDIF
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 23
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 23'
    ENDIF
    IF (irec(irec_occ).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ( (irecord.NE.0).AND.(irecord.EQ.nkpts0).AND.&
               (irec(irec_sdp).NE.0) ) THEN
             nkmin=MIN(nkpt%nkpts,nkpts0)
             nmin=MIN(crge%n,n0)
             CALL zeroing(crge%f)!,n*nkpt%nkpts)
             DO ik=1,nkmin
                READ(nr) (crge%f(i,ik),i=1,nmin)
             ENDDO
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord-nkmin)
             ENDIF
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO OCCUPATION NUMBERS INFORMATION'
             irec(irec_occ)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ENDIF
       CALL mp_bcast(crge%f,SIZE(crge%f),parai%io_source,parai%cp_grp)
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 24
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 24'
    ENDIF
    IF (irec(irec_eigv).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF ( (irecord.NE.0).AND.(irecord.EQ.nkpts0+1).AND.&
               (irec(irec_sdp).NE.0) ) THEN
             READ(nr) ener_com%amu
             irecord=irecord-1
             nkmin=MIN(nkpt%nkpts,nkpts0)
             nmin=MIN(crge%n,n0)
             CALL zeroing(eigv)!,n*nkpt%nkpts)
             DO ik=1,nkmin
                READ(nr) (eigv(i,ik),i=1,nmin)
             ENDDO
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord-nkmin)
             ENDIF
          ELSE
             WRITE(6,'(A,A)') ' RV30| WARNING! ',&
                  'NO FERMI ENERGY AND EIGENVALUES INFORMATION'
             irec(irec_eigv)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ENDIF
       CALL mp_bcast(ener_com%amu,parai%io_source,parai%cp_grp)
       CALL mp_bcast(eigv,SIZE(eigv),parai%io_source,parai%cp_grp)
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 25
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 25'
    ENDIF
    IF (irec(irec_clas).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.clas3%nclatom) THEN
             WRITE(6,'(A,A)') ' RV30| WARNING! ',&
                  'CANNOT READ CORRECTLY CLASSICAL PARTICLE INFORMATION'
             irec(irec_clas)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ELSE
             DO i=1,clas3%nclatom
                READ(nr) (clasc(j,i),j=1,3),(clasv(j,i),j=1,3)
             ENDDO
          ENDIF
       ENDIF
       CALL mp_bcast(clasc,SIZE(clasc),parai%io_source,parai%cp_grp)
       CALL mp_bcast(clasv,SIZE(clasv),parai%io_source,parai%cp_grp)
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 26
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 26'
    ENDIF
    IF (irec(irec_lrwf).NE.0) THEN
       CALL r_linres (nr,"WF")
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 27
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 27'
    ENDIF
    IF (irec(irec_lrwv).NE.0) THEN
       CALL r_linres (nr,"WV")
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==------------------------------------------------------------==
    ! Section 28
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 28'
       ENDIF
       IF (irec(irec_phes).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord-1.NE.nvar .AND. nvar.NE.0) THEN
             WRITE(6,'(A,A,2I3)') ' RV30| WARNING! ',&
                  'CANNOT READ CORRECTLY PARTIAL HESSIAN INFORMATION'
             WRITE(6,'(A,2I4)') '   NVAR / IRECORD: ',nvar,irecorD
             irec(irec_phes)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
             irecord = 0
          ELSE IF (irecord-1.GT.0 .AND. nvar.EQ.0) THEN
             nvar = irecord-1
             WRITE(6,'(A,A,I3)') ' RV30| WARNING! ',&
                  'SETTING THE CORE SIZE FROM THE HESSIAN SIZE: ',nvaR
          ELSE IF (irecord-1.LE.0) THEN
             WRITE(6,'(A,A)') ' RV30| WARNING! ',&
                  'NO VALID PARTIAL HESSIAN RECORD IN RESTART FILE'
             irec(irec_phes)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
             irecord = 0
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
          irecord = 0
       ENDIF
       IF (irecord.NE.0) THEN
          IF(ALLOCATED(hesscr)) THEN
             DEALLOCATE(hesscr,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          ALLOCATE(hesscr(nvar,nvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          nshess = 0
          lallhs = .FALSE.
          IF (paral%io_parent) THEN
             k=0
             DO i=1,nvar
                READ(nr) (hesscr(j+k,1),j=1,nvar)! if NVAR is set here
                k=k+nvar
             ENDDO
             READ(nr) nshess, lvlhes
             IF (.NOT.lvlhes) THEN
                IF (nshess.GT.0 .AND. MOD(nshess,2).EQ.0) THEN
                   nshess = nshess - 1
                   lundop = .TRUE.
                   WRITE(6,'(A,A,I3)') ' RV30| WARNING! ',&
                        'INCOMPLETE HESSIAN, REDOING DIMENSION ',&
                        nshess/2+1
                ELSE
                   WRITE(6,'(A,A,I3)') ' RV30| WARNING! ',&
                        'INCOMPLETE HESSIAN, DOING DIMENSION ',&
                        nshess/2+1
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 29
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 29'
       ENDIF
       IF (irec(irec_prfo).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.8) THEN
             CALL rd30prfo (nr, irecord, irec(irec_lbfgs))
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO P-RFO INFORMATION'
             irec(irec_lbfgs) = 0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 30
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 30'
       ENDIF
       IF (irec(irec_lbfgs).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.GE.7) THEN
             CALL rd30lbfgs (nr, irecord)
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO L-BFGS INFORMATION'
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 31
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 31'
       ENDIF
       IF (irec(irec_rdtl).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.1) THEN
             READ(nr) adtolr%delast, adtolr%demin, adtolr%tolmin, cntr%tolog, cntr%tolad, cntr%tolene
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO ADAPTIVE TOLERANCE INFORMATION'
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(adtolr%delast,parai%io_source,parai%cp_grp)
    CALL mp_bcast(adtolr%tolmin,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntr%tolog,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntr%tolad,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntr%tolene,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! Section 32
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 32'
       ENDIF
       IF (irec(irec_cons).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.2) THEN
             READ(nr) mcnstr0
             IF (mcnstr0.EQ.cotc0%mcnstr) THEN
                READ(nr) (cnsval(i),i=1,cotc0%mcnstr)
             ELSE
                READ(nr)
                WRITE(6,'(A)')&
                     ' RV30| WARNING! INCONSISTENT CONSTRAINT INFORMATION'
             ENDIF
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO CONSTRAINT INFORMATION'
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 33
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 33'
       ENDIF
       IF (irec(irec_ctrans).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.1) THEN
             WRITE(6,'(2A)') ' RV30| NO CELL TRANSLATION DATA. ',&
                  'INITIALIZING TO (0.0,0.0,0.0).'
             irec(irec_ctrans)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
             clsaabox%mm_c_trans(1)=0.0_real_8
             clsaabox%mm_c_trans(2)=0.0_real_8
             clsaabox%mm_c_trans(3)=0.0_real_8
          ELSE
             READ(nr) (clsaabox%mm_c_trans(i),i=1,3)
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(clsaabox%mm_c_trans,SIZE(clsaabox%mm_c_trans),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! Section 34
    ! FIXME LOCT: READ GROUPS AND PARAMETERS, CHECK PARALELL
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 34'
       ENDIF
       IF (irec(irec_nop4).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.2) THEN
             READ(nr) ((lcttab(j,i),i=1,maxsys%nax),j=1,maxsys%nsx)
             READ(nr) ((etap(j,i,ip),etapdot(j,i,ip),i=1,nchx),&
                  j=1,maxsys%nax*maxsys%nsx)
          ELSE
             ! vw here we need to stop as the arrays are not initiaized at all
             CALL stopgm(proceduren,'NO NOSE INFORMATION 34',& 
                  __LINE__,__FILE__)
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO NOSE INFORMATION 34'
             irec(irec_nop2)=0
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irecord)
             ENDIF
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 35
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 35'
       ENDIF
       IF (irec(irec_shop).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.EQ.5) THEN
             READ(nr) sh03%isurf,(sh03%pop(i),i=1,6),sh03%det,sh03%detold
             DO i=1,2
                READ(nr) (sh03%coupl(i,k),sh03%couplold(i,k),k=1,2)
                READ(nr) sh03%ec(i),sh03%eold(i)
             ENDDO
          ELSE
             WRITE(6,'(A)')&
                  ' RV30| WARNING! NO SURFACE HOPPING INFORMATION'
             irec(irec_shop)=0
             CALL stopgm(proceduren,'ERROR READING SHOP',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(sh03%isurf,parai%io_source,parai%cp_grp)
    CALL mp_bcast(sh03%pop,SIZE(sh03%pop),parai%io_source,parai%cp_grp)
    CALL mp_bcast(sh03%coupl,SIZE(sh03%coupl),parai%io_source,parai%cp_grp)
    CALL mp_bcast(sh03%couplold,SIZE(sh03%couplold),parai%io_source,parai%cp_grp)
    CALL mp_bcast(sh03%ec,SIZE(sh03%ec),parai%io_source,parai%cp_grp)
    CALL mp_bcast(sh03%eold,SIZE(sh03%eold),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! Section 36
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 36'
    ENDIF
    IF (irec(irec_shopbo).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
       ENDIF
       CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
       IF (irecord.EQ.0) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(A)') 'RV30| WARNING! NO SHOP WAVEFUNCTIONS'
             CALL stopgm(proceduren,'ERROR READING PREVIOUS'&
                  //' WAVEFUNCTION',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF (irecord.LT.0) THEN
          CALL r_wfnio(nr,c0old,nstate,info,tkpnt0,n0,nkpts0,&
               'C0',fpos)
       ELSE
          CALL rd30wfn(nr,c0old,nstate,tau0,info,&
               TKPNT0,N0,NKPTS0,'C0',FPOS)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 37
    prngin=0.0_real_8
    prngout=0.0_real_8
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 37'
       ENDIF
       IF (irec(irec_prng).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.LT.2) THEN
             WRITE(6,'(A)') 'RV30| WARNING! NO PRNG STATE'
             CALL stopgm(proceduren,'ERROR READING PREVIOUS PRNG',& 
                  __LINE__,__FILE__)
          ELSE IF (irecord.NE.(1+parai%nproc)) THEN
             WRITE(6,'(A)')&
                  'RV30| WARNING! NUMPROC CHANGED, PARALLEL PRNG WILL BE RESETTED'
             READ(nr) (prng_com%repseed(1,i),i=1,3),(prng_com%repseed(2,i),i=1,3)
             READ(nr) (prng_com%paraseed(1,i),i=1,3),(prng_com%paraseed(2,i),i=1,3)
             DO j=1,irecord-2
                IF (j<=parai%nproc) THEN
                   READ(nr) (prngin(j,1,i),i=1,3), (prngin(j,2,i),i=1,3)
                ELSE
                   READ(nr) (buff(i),i=1,3), (buff(i),i=1,3)
                ENDIF
             ENDDO
          ELSE
             READ(nr) (prng_com%repseed(1,i),i=1,3),(prng_com%repseed(2,i),i=1,3)
             DO j=1,parai%nproc
                READ(nr) (prngin(j,1,i),i=1,3),(prngin(j,2,i),i=1,3)
             ENDDO
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF

    CALL mp_bcast(prng_com%repseed,6,parai%io_source,parai%cp_grp)
    CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)
    IF (irecord.NE.(1+parai%nproc)) THEN
       CALL mp_bcast(prng_com%paraseed,6,parai%io_source,parai%cp_grp)
       CALL prngparaskip(prng_com%paraseed,prng_com%paraseed)! skips so as to reconstruct the paraseed on individual nodes
    ELSE
       CALL mp_sum(prngin,prngout,parai%nproc*6,parai%cp_grp)
       IF (cntl%tpath.OR.tmw) THEN
          prng_com%paraseed=prngout(parai%mepos+1,:,:)
       ELSE
          prng_com%paraseed=prngout(parai%me+1,:,:)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 38
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 38'
       ENDIF
       IF (irec(irec_gle).NE.0) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.(1+ions1%nat*(glepar%gle_ns+1))) THEN
             WRITE(6,'(A)') 'RV30| WARNING! NO GLE DATA OR NS CHANGED'
             CALL stopgm(proceduren,'ERROR READING PREVIOUS GLE DATA',& 
                  __LINE__,__FILE__)
          ELSE
             READ(nr) glepar%egle
             DO is = 1 , ions1%nsp
                DO ia = 1, ions0%na(is)
                   DO i= 1, glepar%gle_ns+1
                      READ(nr) (glep(j,ia,is,i),j=1,3)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ELSE
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          irec0 = ABS(irec0)
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(irec(irec_gle),parai%io_source,parai%cp_grp)

    IF (irec(irec_gle).NE.0) THEN
       CALL mp_bcast(glepar%egle,parai%io_source,parai%cp_grp)

       CALL mp_bcast(glep,SIZE(glep),parai%io_source,parai%cp_grp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 39 - 98
    IF (paral%io_parent) THEN
       IF (store1%tdebio) THEN
          WRITE(6,*) 'READ SECTION 38-98'
       ENDIF
       DO i=39,98
          IF (cntl%is_in_stream)THEN
             READ(nr) irec0,fpos
          ELSE
             READ(nr) irec0
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Section 99
    IF (store1%tdebio.AND.paral%io_parent) THEN
       WRITE(6,*) 'READ SECTION 99'
    ENDIF
    numcold=0
    IF (irec(irec_xtrp).NE.0) THEN
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.NE.1) THEN
             WRITE(6,'(2A)') ' RV30| INCONSISTENT WF HISTORY DATA. ',&
                  ' SKIPPING AHEAD...'
             irec(irec_xtrp)=0
             numcold=0
          ELSE
             READ(nr) numcold
          ENDIF
       ENDIF
       nnow=numcold
       CALL mp_bcast(numcold,parai%io_source,parai%cp_grp)

       IF (numcold.GT.0) THEN
          DO i=1,numcold
             IF (paral%io_parent) THEN
                IF (cntl%is_in_stream)THEN
                   READ(nr) irecord,fpos
                ELSE
                   READ(nr) irecord
                ENDIF
             ENDIF
             CALL mp_bcast(irecord,parai%io_source,parai%cp_grp)

             IF (i.LE.cnti%mextra) THEN
                IF (store1%tdebio.AND.paral%io_parent) THEN
                   WRITE(6,'(A,I2,A,I2)') ' RV30| READING OLD WFN   ',&
                        i,' OF ',numcolD
                ENDIF
                IF (irecord.LT.0) THEN
                   CALL r_wfnio(nr,cold(1,1,1,i),nstate,info,tkpnt0,&
                        n0,nkpts0,'C0',fpos)
                ELSE
                   info=1
                ENDIF
             ELSE
                IF (paral%io_parent) THEN
                   WRITE(6,'(A,I2,A,I2)') ' RV30| SKIPPING OLD WFN  ',&
                        i,' OF ',numcolD
                   irec0 = ABS(irecord)
                   IF (cntl%is_in_stream) THEN
                      CALL file_seek(nr,.FALSE.,fpos)
                   ELSE
                      CALL fskip(nr,irec0)
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
          numcold=MIN(numcold,cnti%mextra)
       ELSE
          info=1
       ENDIF

       IF (info.NE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' RV30| WARNING! NO OR CORRUPT WAVEFUNCTION HISTORY'
          irec(irec_xtrp)=0
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          IF (irecord.GT.1) THEN
             WRITE(6,'(2A)') ' RV30| INCONSISTENT WF HISTORY DATA. ',&
                  ' SKIPPING AHEAD...'
             irec(irec_xtrp)=0
             numcold=0
          ELSEIF (irecord.EQ.0) THEN
             numcold=0
             irec(irec_xtrp)=0
          ELSE
             READ(nr) numcold
          ENDIF

          DO i=1,numcold
             IF (cntl%is_in_stream)THEN
                READ(nr) irec0,fpos
             ELSE
                READ(nr) irec0
             ENDIF
             irec0 = ABS(irec0)
             IF (store1%tdebio) THEN
                WRITE(6,*) 'SKIPPING WFN ',i,' OF ',numcold
             ENDIF
             IF (cntl%is_in_stream) THEN
                CALL file_seek(nr,.FALSE.,fpos)
             ELSE
                CALL fskip(nr,irec0)
             ENDIF
          ENDDO
          numcold=0
       ENDIF
    ENDIF
    nnow=numcold
    CALL mp_bcast(numcold,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! Broadcast information of IREC
    CALL mp_bcast(irec,SIZE(irec),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rv30
  ! ==================================================================
  SUBROUTINE rd30kpnt(nr,nstate,nkpts0,tkpnt0,c0,cm,tau0,irec)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, nstate, nkpts0
    LOGICAL                                  :: tkpnt0
    COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt), &
      cm(nkpt%ngwk,nstate,nkpt%nkpnt)
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: irec(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rd30kpnt'
    REAL(real_8), PARAMETER                  :: eps = 1.e-3_real_8 

    INTEGER :: icompr, ierr, ik, ikind, ikpoint, ikpt, index, iperm, ir, &
      irecord, istate, j, len, len1, n0, ngwa, ngwks0, nn, nread, nstate0
    CHARACTER(len=80)                        :: str
    INTEGER, ALLOCATABLE                     :: icmp(:), itabkpt(:), &
                                                itable(:), mapw(:)
    REAL(real_8)                             :: scale, wk0
    REAL(real_8), ALLOCATABLE                :: ca(:), cr(:), rk0(:,:)

! ==--------------------------------------------------------------==
! If options restart kpoints, the number of kpoints has to be equal.

    IF (restart1%rkpt) THEN
       IF (nkpt%nkpts.NE.nkpts0) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)&
                  ' RD30KPNT| WITH RESTART KPOINTS, THE NUMBER OF KPOINTS'
             WRITE(6,*)&
                  '           HAS TO BE THE SAME (NKPTS-OLD=',nkpts0,')'
          ENDIF
          CALL stopgm(' RD30KPNT|',&
               'THE NUMBER OF KPOINTS ARE DIFFERENT ! ',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent) CALL fskip(nr,nkpt%nkpts)
       RETURN
    ENDIF
    ALLOCATE(rk0(3,nkpts0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(itabkpt(nkpts0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(itable(nkpts0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! No options restart.
    IF (paral%io_parent) THEN
       DO ik=1,nkpts0
          READ(nr) (rk0(j,ik),j=1,3),wk0
       ENDDO
    ENDIF
    CALL mp_bcast(rk0,SIZE(rk0),parai%io_source,parai%cp_grp)
    ! Compare each kpoints in order to use the maximum information.
    ikpoint=0
    iperm=0
    CALL zeroing(itable)!,nkpts0)
    DO ikpt=1,nkpt%nblkp
       DO ikind=1,nkpbl(ikpt)
          ik=kpbeg(ikpt)+ikind
          index=1
          DO WHILE(index.LE.nkpts0)
             IF ((ABS(rk0(1,index)-rk(1,ik)).LT.eps).AND.&
                  (ABS(rk0(2,index)-rk(2,ik)).LT.eps).AND.&
                  (ABS(rk0(3,index)-rk(3,ik)).LT.eps)) THEN
                ikpoint=ikpoint+1
                IF ((index.NE.ik)) THEN
                   IF ( (irec(irec_wf).NE.0).OR.(irec(irec_wf).NE.0) ) THEN
                      iperm=iperm+1
                      itable(index)=ikind
                      itabkpt(index)=ikpt
                   ENDIF
                ENDIF
                index=nkpts0+1
             ELSE
                index=index+1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! Permutation if necessary
    IF ( (iperm.NE.0).AND.(irec(irec_wf).NE.0) ) THEN
       ! We need to reread the file
       IF (paral%io_parent) THEN
          REWIND(nr)
          READ(nr) str
          DO ir=2,8
             READ(nr) irecord
             CALL fskip(nr,irecord)
          ENDDO
          READ(nr) irecord
          READ(nr) n0,ngwks0,icompr,scale
       ENDIF
       CALL mp_bcast(n0,parai%io_source,parai%cp_grp)
       cnti%rcompb=icompr
       IF (cnti%rcompb.LT.0) THEN
          CALL aotopw(nr,c0,nstate,tau0,icompr,scale)
       ELSE
          IF (paral%io_parent) THEN
             len=2*nkpt%ngwk
             ALLOCATE(cr(4*len),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             len1=2*ncpw%ngw + 1
             ALLOCATE(mapw(len1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ngwa=MAX(ngwks0,spar%ngwks)
             ALLOCATE(ca(2*ngwa+100),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)! TODO why +100?
             ALLOCATE(icmp(2*ngwa+100),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(ca)!,2*ngwa+100)
             CALL zeroing(icmp)!,2*ngwa+100)
          ENDIF
          nstate0=n0/nkpts0
          nread=0
          ikpoint=0
          DO index=1,nkpts0
             ikind=itable(index)
             IF (ikind.NE.0) THEN
                ikpoint=ikpoint+1
                IF (tkpts%tkblock) THEN
                   ikpt=itabkpt(index)
                   CALL rkpt_swap(c0,nstate,ikpt,'C0')
                ENDIF
                DO istate=1,nstate0
                   nread=nread+1
                   IF (nread.LE.n0) THEN
                      CALL rdwfn(nr,c0(1,istate,ikind),ca,cr,icmp,mapw,&
                           ngwks0,icompr,scale,tkpnt0,.FALSE.)
                   ELSE
                      GOTO 30
                   ENDIF
                ENDDO
                IF (tkpts%tkblock) CALL wkpt_swap(c0,nstate,ikpt,'C0')
             ELSE
                nn=MIN(n0-nread,nstate0)
                IF (nn.LE.0) GOTO 30
                IF (paral%parent) CALL fskip(nr,2*nn)
                nread=nread+nstate
             ENDIF
          ENDDO
          DEALLOCATE(ca,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(icmp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(mapw,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
30  CONTINUE
    IF (ikpoint.NE.nkpt%nkpts) THEN
       IF (nkpt%nkpts-ikpoint.EQ.1) THEN
          WRITE(6,'(A,I3,A)')&
               ' RV30! ',nkpt%nkpts-ikpoint,' KPOINT IS DIFFERENT'
       ELSE
          WRITE(6,'(A,I3,A)')&
               ' RV30! ',nkpt%nkpts-ikpoint,' KPOINTS ARE DIFFERENT'
       ENDIF
    ENDIF
    DEALLOCATE(rk0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(itabkpt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(itable,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rd30kpnt
  ! ==================================================================
  SUBROUTINE rdkpoints
    ! ==--------------------------------------------------------------==
    ! == READ NKPTS, RK AND WK IN VERSION 3.0                         ==
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'rdkpoints'

    CHARACTER(len=20)                        :: latfile
    CHARACTER(len=80)                        :: str
    CHARACTER(len=fo_fpath_max)              :: filnam
    INTEGER                                  :: ia, ie, ierr, ik, ir, &
                                                irecord, j, n0, ngwks0, &
                                                ngwls0, nhgls0, nhgs0, &
                                                nkpts0, nr
    INTEGER(int_8)                           :: fpos
    LOGICAL                                  :: fexist

! ==--------------------------------------------------------------==

    nr=1
    CALL xstring(filn,ia,ie)
    IF (cntl%tpath) THEN
       filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
    ELSE
       IF (restart1%rlate) THEN
          IF (tmw)THEN
             CALL mw_filename('LATEST_',latfile,mwi%walker_id)
          ELSE
             latfile='LATEST'
          ENDIF
          INQUIRE(file=latfile,exist=fexist)
          IF (fexist) THEN
             OPEN(unit=nr,file=latfile,status='UNKNOWN')
             READ(nr,'(A)') filnam
             CLOSE(nr)
          ELSE
             filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
          ENDIF
       ELSE
          filnam=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filn(ia:ie)
       ENDIF
    ENDIF
    filnam=ADJUSTL(filnam)
    INQUIRE(file=filnam,exist=fexist)
    IF (.NOT.fexist) THEN
       WRITE(6,*) ' RDKPOINTS| RESTART FILE NOT FOUND:',filnam
       CALL stopgm('RDKPOINTS',' FILE NOT FOUND ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL file_open(filnam,'UNFORMATTED','OLD','READ',cntl%is_in_stream,&
         .FALSE.,0,nr)
    READ(nr) str
    IF (INDEX(str,'STUTTGART VERSION 3.0').EQ.0) THEN
       CALL stopgm('RDKPOINTS',&
            'NOT RESTART WITH KPOINTS WITH THIS FORMAT',& 
            __LINE__,__FILE__)
    ENDIF
    DO ir=2,7
       IF (cntl%is_in_stream)THEN
          READ(nr) irecord,fpos
       ELSE
          READ(nr) irecord
       ENDIF
       IF (cntl%is_in_stream) THEN
          CALL file_seek(nr,.FALSE.,fpos)
       ELSE
          CALL fskip(nr,irecord)
       ENDIF
    ENDDO
    IF (cntl%is_in_stream)THEN
       READ(nr) irecord,fpos
    ELSE
       READ(nr) irecord
    ENDIF
    IF (irecord.EQ.1) THEN
       READ(nr) n0,nkpts0,ngwks0,ngwls0,nhgs0,nhgls0
       IF (nkpts0.EQ.0) THEN
          nkpt%nkpts=1
          tkpts%tkpnt=.FALSE.
       ELSE
          nkpt%nkpts=nkpts0
          tkpts%tkpnt=.TRUE.
       ENDIF
       DO ir=9,20
          IF (cntl%is_in_stream)THEN
             READ(nr) irecord,fpos
          ELSE
             READ(nr) irecord
          ENDIF
          ! New format to read PW coefficients.
          IF (irecord.LT.0) THEN
             irecord = -irecord
          ENDIF
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irecord)
          ENDIF
       ENDDO
       ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (cntl%is_in_stream)THEN
          READ(nr) irecord,fpos
       ELSE
          READ(nr) irecord
       ENDIF
       IF (irecord.EQ.0.AND.nkpt%nkpts.EQ.1.OR.nkpts0.EQ.0) THEN
          tkpts%tkpnt=.FALSE.
          DO j=1,3
             rk(j,1)=0._real_8
          ENDDO
          wk(1)=1._real_8
       ELSEIF (irecord.EQ.nkpt%nkpts) THEN
          DO ik=1,nkpt%nkpts
             READ(nr) (rk(j,ik),j=1,3),wk(ik)
          ENDDO
       ELSE
          WRITE(6,'(A)')&
               ' RDKPOINTS| WARNING! CANNOT READ KPOINTS SECTION'
          CALL stopgm('RDKPOINTS','NOT RESTART WITH KPOINTS',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       WRITE(6,'(A)')&
            ' RDKPOINTS| WARNING! CANNOT READ NUMBER OF KPOINTS'
       CALL stopgm('RDKPOINTS','NOT RESTART WITH KPOINTS',& 
            __LINE__,__FILE__)
    ENDIF
    CLOSE(nr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rdkpoints
  ! ==================================================================


END MODULE rv30_utils


! ==================================================================
SUBROUTINE rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr)
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  ! == Read Potential or density array in RESTART file              ==
  ! == Works if the dimensions of reading array are not             ==
  ! == the same as POTR (but no guarantee)                          ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send, mp_recv
  USE system , ONLY:fpar,spar, parap
  USE parac, ONLY : paral,parai
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nr, irecord, nr1s0, nr2s0, &
                                                nr3s0
  REAL(real_8) :: potr(fpar%kr1,fpar%kr2s*fpar%kr3s)

  CHARACTER(*), PARAMETER                    :: procedureN = 'rd30pot'

  INTEGER                                    :: i, ierr, ip, ipp, ipx, ir, &
                                                irr, kk, kkmax, kks0, kr1s0, &
                                                kr2s0, kr3s0, msgid
  REAL(real_8), ALLOCATABLE                  :: pscr(:)

! Variables
! (kr2s*kr3s)
! ==--------------------------------------------------------------==

  kk=fpar%kr2s*fpar%kr3s
  ! Added to avoid a bug in T3E (Thierry Deutsch 28/04/97)
  CALL zeroing(potr)!,kr1*kk)
  ! Dimensions for POT in reading file.
  kr1s0=nr1s0+MOD(nr1s0+1,2)
  kr2s0=nr2s0+MOD(nr2s0+1,2)
  kr3s0=nr3s0+MOD(nr3s0+1,2)
  kks0=kr2s0*kr3s0
  kkmax=MAX(kks0,kk)
  ALLOCATE(pscr(kkmax),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  DO ir=1,irecord
     IF (paral%io_parent) READ(nr) (pscr(i),i=1,kks0)
     IF (kks0.LT.kk) CALL zeroing(pscr(kks0+1:kk))!,kk-kks0)
     IF (irecord.GT.spar%nr1s) THEN
        GOTO 20
     ELSE
        DO ipp=1,parai%nproc
           IF (ir.GE.parap%nrxpl(ipp-1,1).AND.ir.LE.parap%nrxpl(ipp-1,2)) THEN
              ipx=ipp
              GOTO 10
           ENDIF
        ENDDO
        CALL stopgm('RD30POT','ERROR IN NRXPL',& 
             __LINE__,__FILE__)
     ENDIF
10   CONTINUE
     ip=parap%pgroup(ipx)
     CALL mp_sync(parai%allgrp)
     IF (ip.NE.parai%source) THEN
        IF (paral%parent) THEN
           msgid=1
           !msglen=kk * 8
           CALL mp_send(pscr,kk,ip,msgid,parai%allgrp)
        ELSEIF (parai%me.EQ.ip) THEN
           msgid=1
           !msglen=kk * 8
           CALL mp_recv(pscr,kk,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
        CALL mp_sync(parai%allgrp)
     ENDIF
     IF (ip.EQ.parai%me) THEN
        irr=ir-parap%nrxpl(parai%mepos,1)+1
        CALL dcopy(kk,pscr(1),1,potr(irr,1),fpar%kr1)
     ENDIF
20   CONTINUE
  ENDDO
  DEALLOCATE(pscr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE rd30pot
