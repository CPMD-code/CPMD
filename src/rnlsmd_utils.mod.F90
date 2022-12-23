MODULE rnlsmd_utils
  USE cppt,                            ONLY: twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mm_dimmod,                       ONLY: mmdim
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: eigr,&
                                             fnl,&
                                             fnl2,&
                                             tfnl2
  USE sumfnl_utils,                    ONLY: sumfnld
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsmd

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsmd(c0,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY FNL (ALSO CALLED BEC IN SOME VANDERBILT ROUTINES) ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), INTENT(in)              :: c0(:,:)
    INTEGER, INTENT(in)                      :: nstate, ikind

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsmd'
    INTEGER, PARAMETER                       :: nspar = 50

    COMPLEX(real_8), ALLOCATABLE             :: scra(:,:)
    INTEGER                                  :: i, ierr, is, isa0, isub, iv, &
                                                l, nstat
    INTEGER, SAVE                            :: ilfnlx = 0
    REAL(real_8), ALLOCATABLE, SAVE          :: fnlx(:)

    IF (tkpts%tkpnt) CALL stopgm('RNLSMD','TSDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tiset('    RNLSMD',isub)
    IF (ilfnlx.EQ.0) THEN
       ilfnlx=MAX(1,2*mmdim%naxq*nspar) !vw 2* complex -> real
       ALLOCATE(fnlx(ilfnlx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(scra(nkpt%ngwk,mmdim%naxq),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(scra)
    CALL zeroing(fnl(:,:,:,:,ikind))!,imagp*ldf1*maxsys%nhxs*nstate)
    IF (tfnl2) CALL zeroing(fnl2(:,:,:,:,ikind))!,imagp*ions1%nat*maxsys%nhxs*ldf2)
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          l=nghtol(iv,is)
          CALL scrph1(eigr(:,:,1),scra,twnl(:,iv,is,1),l,ions0%na(is),isa0)
          DO i=1,nstate,nspar
             nstat=MIN(nspar,nstate-i+1)
             IF (ncpw%ngw>0) THEN
                CALL dgemm('T','N',ions0%na(is),nstat,2*ncpw%ngw,2.0_real_8,&
                     scra,2*ncpw%ngw,c0(1,i),2*ncpw%ngw,0.0_real_8,fnlx,ions0%na(is))
             ENDIF
             CALL sumfnld(fnlx,isa0,ions0%na(is),iv,i,nstat,nstate,1,1,parai%allgrp)
          ENDDO
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO

    DEALLOCATE(scra,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('    RNLSMD',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rnlsmd
  ! ==================================================================
  SUBROUTINE scrph1(eigr,scra,twnl,l,natx,isa0)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), INTENT(in)              :: eigr(:,:)
    COMPLEX(real_8), INTENT(out)             :: scra(:,:)
    REAL(real_8)                             :: twnl(:)
    INTEGER, INTENT(in)                      :: l, natx, isa0

    COMPLEX(real_8)                          :: ci
    INTEGER                                  :: ia, ig, isa
    REAL(real_8)                             :: arg, cii, cir, ei, er

! ==--------------------------------------------------------------==

    ci=(0.0_real_8,-1.0_real_8)**l
    cir=REAL(ci)
    cii=AIMAG(ci)
    !$omp parallel do private(IA,IG,ISA,ARG,ER,EI)
    DO ia=1,natx
       isa=isa0+ia
       ! ..MAKE USE OF THE SPECIAL STRUCTURE OF CI
       IF (ABS(cir).GT.0.5_real_8) THEN
          ! ..CI IS REAL
          DO ig=1,ncpw%ngw
             arg=twnl(ig)*cir
             er=REAL(eigr(ig,isa))
             ei=AIMAG(eigr(ig,isa))
             scra(ig,ia) = CMPLX(arg*er,arg*ei,kind=real_8)
          ENDDO
       ELSE
          ! ..CI IS IMAGINARY
          DO ig=1,ncpw%ngw
             arg=twnl(ig)*cii
             er=REAL(eigr(ig,isa))
             ei=AIMAG(eigr(ig,isa))
             scra(ig,ia) = CMPLX(-arg*ei,arg*er,kind=real_8)
          ENDDO
       ENDIF
       IF (geq0) scra(1,ia)=0.5_real_8*scra(1,ia)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE scrph1
  ! ==================================================================

END MODULE rnlsmd_utils
