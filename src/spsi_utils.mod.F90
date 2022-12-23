MODULE spsi_utils
  USE cppt,                            ONLY: twnl
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: eigr,&
                                             fnl
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spsi
  PUBLIC :: give_scr_spsi

CONTAINS

  ! ==================================================================
  SUBROUTINE spsi(nstate,sc0)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate)
    INTEGER,   allocatable                   :: nivmap(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'spsi'

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: work(:,:)
    REAL(real_8)    ,ALLOCATABLE             :: aa(:,:)
    INTEGER                                  :: ia, ierr, ig, is, isa, isa0, &
                                                isub, iv, j, jv, idx, niv
    REAL(real_8)                             :: t1, t2, t3

    CALL tiset('      SPSI',isub)
    IF (cntl%tfdist) CALL stopgm('SPSI','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2) CALL stopgm('SPSI','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ALLOCATE(aa(ions1%nat*maxsys%nhxs,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(ncpw%ngw,ions1%nat*maxsys%nhxs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(nivmap(2,ions1%nat*maxsys%nhxs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    idx=0
    DO isa=1,ions1%nat
       is=iatpt(2,isa)
       IF (pslo_com%tvan(is)) THEN
          DO iv=1,nlps_com%ngh(is)
             idx=idx+1
             nivmap(1,idx)=isa
             nivmap(2,idx)=iv
          ENDDO
       ENDIF
    ENDDO
    niv=idx
!$OMP parallel private(ISA,IS,IV,IA,CI,J,JV)
!$OMP do private(IS,IV,JV)
    DO idx=1,niv
        isa=nivmap(1,idx)
        is=iatpt(2,isa)
        iv=nivmap(2,idx)
        DO j=1,nstate
            AA(idx,j)=0.0_real_8
            DO jv=1,nlps_com%ngh(is)
                aa(idx,j)=aa(idx,j)+qq(iv,jv,is)*fnl(1,isa,jv,j,1)
            ENDDO
        ENDDO
    ENDDO
!$OMP end do    
    DO idx=1,niv
        isa=nivmap(1,idx)
        is=iatpt(2,isa)
        iv=nivmap(2,idx)
        ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
!       MAKE USE OF THE SPECIAL STRUCTURE OF CI
        IF (ABS(REAL(ci)).GT.0.5_real_8) THEN
!           CI IS REAL
!$OMP do SCHEDULE(STATIC) private(t1,t2,t3)
            DO ig=1,ncpw%ngw
               t1=REAL(eigr(ig,isa,1))
               t2=AIMAG(eigr(ig,isa,1))
               t3=twnl(ig,iv,is,1)*REAL(ci)
               work(ig,idx)=CMPLX(t1*t3,t2*t3,kind=real_8)
            ENDDO
!$OMP end do
        ELSE
!          CI IS IMAGINARY
!$OMP do SCHEDULE(STATIC) private(t1,t2,t3)
           DO ig=1,ncpw%ngw
              t1=REAL(eigr(ig,isa,1))
              t2=AIMAG(eigr(ig,isa,1))
              t3=twnl(ig,iv,is,1)*AIMAG(ci)
              work(ig,idx)=CMPLX(-t2*t3,t1*t3,kind=real_8)
           ENDDO
!$OMP end do
        ENDIF
    ENDDO
!$OMP end parallel
    IF(ncpw%ngw /= 0) THEN
        CALL DGEMM('n','n',2*ncpw%ngw,nstate, &
            niv, 1.0_real_8,work(1,1),2*ncpw%ngw,&
            aa(1,1),ions1%nat*maxsys%nhxs,&
            1.0_real_8,sc0(1,1),2*ncpw%ngw)
    ENDIF
    DEALLOCATE(aa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(nivmap,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('      SPSI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE spsi
  ! ==================================================================
  SUBROUTINE give_scr_spsi(lspsi,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lspsi
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lspsi=ncpw%ngw*2
    tag ='NGW*2'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_spsi
  ! ==================================================================

END MODULE spsi_utils
