MODULE wrintf_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_ufo
  USE func,                            ONLY: func1
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE parac,                           ONLY: paral
  USE ragg,                            ONLY: raggio
  USE readsr_utils,                    ONLY: xstring
  USE store_types,                     ONLY: intfn
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             dual00,&
                                             nkpt,&
                                             parm,&
                                             spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wintf
  !public :: wintf36

CONTAINS

  ! ==================================================================
  SUBROUTINE wintf (nw,c0,taux)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nw
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*)
    REAL(real_8)                             :: taux(:,:,:)

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, ie
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       CALL xstring(intfn,ia,ie)
       IF (cntl%tpath) CALL stopgm("WINTF","Not implemented",& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            CALL fileopen(nw,intfn(ia:ie),fo_def+fo_ufo,ferror)
       IF (paral%io_parent)&
            REWIND(nw)
    ENDIF
    CALL wintf36(nw,c0,taux)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(nw)
       IF (paral%io_parent)&
            WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(38,65-(ie-ia)),',A)'
       IF (paral%io_parent)&
            WRITE(6,fformat)&
            ' RESTART INFORMATION WRITTEN ON FILE ',intfn(ia:ie)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wintf
  ! ==================================================================
  SUBROUTINE wintf36 (nw,c0,taux)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nw
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*)
    REAL(real_8)                             :: taux(:,:,:)

    INTEGER                                  :: ik, ipx, ipy, ipz, is, lines

    IF (paral%parent) THEN
       ! FILE TAG
       IF (paral%io_parent)&
            WRITE(nw) "CPMD3.6"
    ENDIF
    IF (paral%parent) THEN
       ! UNIT CELL PARAMETERS
       lines=3
       IF (paral%io_parent)&
            WRITE(nw) "CELL PARAMETER ATOMIC"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) parm%a1(1),parm%a2(1),parm%a3(1)
       IF (paral%io_parent)&
            WRITE(nw) parm%a1(2),parm%a2(2),parm%a3(2)
       IF (paral%io_parent)&
            WRITE(nw) parm%a1(3),parm%a2(3),parm%a3(3)
    ENDIF
    IF (paral%parent) THEN
       ! PERIODICITY
       lines=1
       ipx=1
       ipy=1
       ipz=1
       IF (isos1%tclust) THEN
          IF (isos1%toned) THEN
             ipy=0
             ipz=0
          ELSEIF (isos1%ttwod) THEN
             ipz=0
          ELSE
             ipx=0
             ipy=0
             ipz=0
          ENDIF
       ENDIF
       IF (paral%io_parent)&
            WRITE(nw) "PERIODIC"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) ipx,ipy,ipz
    ENDIF
    IF (paral%parent) THEN
       ! WAVEFUNCTION CUTOFF
       lines=1
       IF (paral%io_parent)&
            WRITE(nw) "CUTOFF WAVEFUNCTION RYDBERG"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) cntr%ecut
    ENDIF
    IF (paral%parent) THEN
       ! DENSITY CUTOFF
       lines=1
       IF (paral%io_parent)&
            WRITE(nw) "CUTOFF DENSITY RYDBERG"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) dual00%cdual*cntr%ecut
    ENDIF
    IF (paral%parent) THEN
       ! SPIN POLARISATION
       lines=1
       IF (paral%io_parent)&
            WRITE(nw) "SPIN POLARISATION"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) cntl%tlsd
    ENDIF
    IF (paral%parent) THEN
       ! FUNCTIONAL
       lines=1
       IF (paral%io_parent)&
            WRITE(nw) "FUNCTIONAL"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) func1%mfxcx,func1%mfxcc,func1%mgcx,func1%mgcc,func1%mhfx
    ENDIF
    IF (paral%parent) THEN
       ! KPOINTS
       IF (paral%io_parent)&
            WRITE(nw) "KPOINTS"
       IF (tkpts%tkpnt) THEN
          lines=2+nkpt%nkpnt
          IF (paral%io_parent)&
               WRITE(nw) lines
          IF (paral%io_parent)&
               WRITE(nw) "COMPLEX"
          IF (paral%io_parent)&
               WRITE(nw) nkpt%nkpnt
          DO ik=1,nkpt%nkpnt
             IF (paral%io_parent)&
                  WRITE(nw) wk(ik),rk(1,ik),rk(2,ik),rk(3,ik)
          ENDDO
       ELSE
          lines=1
          IF (paral%io_parent)&
               WRITE(nw) lines
          IF (paral%io_parent)&
               WRITE(nw) "GAMMA"
       ENDIF
    ENDIF
    IF (paral%parent) THEN
       ! ATOM TYPES
       lines=1+ions1%nsp
       IF (paral%io_parent)&
            WRITE(nw) "ATOM TYPES"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) ions1%nsp
       DO is=1,ions1%nsp
          IF (paral%io_parent)&
               WRITE(nw) ions0%iatyp(is),NINT(ions0%zv(is)),raggio(is)
       ENDDO
    ENDIF
    IF (paral%parent) THEN
       ! G-VECTORS
       lines=1
       IF (paral%io_parent)&
            WRITE(nw) "PLANE WAVES"
       IF (paral%io_parent)&
            WRITE(nw) lines
       IF (paral%io_parent)&
            WRITE(nw) spar%ngws
       ! WAVEFUNCTION CUTOFF
       DO is=1,spar%ngws
          IF (paral%io_parent)&
               WRITE(nw)
       ENDDO
       ! DENSITY CUTOFF
       IF (paral%io_parent)&
            WRITE(nw) spar%nhgs
       DO is=1,spar%nhgs
          IF (paral%io_parent)&
               WRITE(nw)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wintf36
  ! ==================================================================

END MODULE wrintf_utils
