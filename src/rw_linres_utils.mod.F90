MODULE rw_linres_utils
  USE error_handling,                  ONLY: stopgm
  USE io_utils,                        ONLY: file_seek
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE linres,                          ONLY: &
       clrv, clrwf, lractive, nlinr, nlinw, nolr, nua, nub, td01, td03, urot
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt,&
                                             spar
  USE utils,                           ONLY: fskip,&
                                             izamax,&
                                             zgive
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: w_linres
  PUBLIC :: r_linres

CONTAINS

  ! ==================================================================
  SUBROUTINE w_linres(nw,ierror,tag,irecord)
    ! ==--------------------------------------------------------------==




    INTEGER                                  :: nw, ierror
    CHARACTER(len=2)                         :: tag
    INTEGER                                  :: irecord

    CHARACTER(*), PARAMETER                  :: procedureN = 'w_linres'

    INTEGER                                  :: i, icompr, ierr, irec0, isca, &
                                                j, k, len, len1, mtda, murec
    INTEGER(int_8)                           :: fpos_null
    INTEGER, ALLOCATABLE                     :: icmp(:), mapw(:)
    LOGICAL                                  :: ostda
    REAL(real_8)                             :: scale
    REAL(real_8), ALLOCATABLE                :: ca(:), cr(:)

    fpos_null=0
    IF (td03%tda) THEN
       mtda=1
    ELSE
       mtda=2
    ENDIF
    IF ( lractive ) THEN
       IF (paral%io_parent) THEN
          len=2*ncpw%ngw
          len1=ncpw%ngw + 1
          ALLOCATE(cr(4*len),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(mapw(2*len1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ca(2*spar%ngws),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(icmp(4*spar%ngws),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! To avoid Fortran runtime 'not allocated'
          ALLOCATE(cr(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(mapw(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ca(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(icmp(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (tag.EQ."WF") THEN
          !>>>
          !vw isca=idamax(2*nlinw*nolr*ngw*mtda,clrwf,1)
          !vw scale=abs(dgive(clrwf,isca))
          isca=izamax(nlinw*nolr*ncpw%ngw*mtda,clrwf,1)
          scale=ABS(zgive(clrwf,isca))
          !<<<
       ELSEIF (tag.EQ."WV") THEN
          !>>>
          !vw isca=idamax(2*nlinw*nolr*ngw*mtda,clrv,1)
          !vw scale=abs(dgive(clrv,isca))
          isca=izamax(nlinw*nolr*ncpw%ngw*mtda,clrv,1)
          scale=ABS(zgive(clrv,isca))
          !<<<
       ENDIF
       CALL mp_max(scale,parai%allgrp)
       icompr=cnti%wcompb
       IF (icompr.LT.1) THEN
          icompr=1
       ENDIF
       IF (td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon) THEN
          ostda=.TRUE.
          murec=nlinw
       ELSE
          ostda=.FALSE.
          murec=0
       ENDIF
       irecord = 2 + 2*nlinw*nolr*mtda + murec
       IF (paral%io_parent) THEN
          WRITE(unit=nw,iostat=ierror) irecord,fpos_null
          WRITE(unit=nw,iostat=ierror) spar%ngws,nolr,nlinw,mtda,icompr,&
               scale
          WRITE(unit=nw,iostat=ierror) ostda,nua,nub
       ENDIF
       DO j=1,mtda
          DO k=1,nlinw
             DO i=1,nolr
                IF (tag.EQ."WF") THEN
                   CALL wrwfn(nw,ierror,clrwf(1,i,k,j),ca,cr,icmp,mapw,&
                        icompr,scale,.FALSE.)
                ELSEIF (tag.EQ."WV") THEN
                   CALL wrwfn(nw,ierror,clrv(1,i,k,j),ca,cr,icmp,mapw,&
                        icompr,scale,.FALSE.)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       IF (ostda.AND.paral%io_parent) THEN
          DO k=1,nlinw
             WRITE(unit=nw,iostat=ierror) ((urot(i,j,k),i=1,nua),j=1,nub)
          ENDDO
       ENDIF
       IF (paral%io_parent) THEN
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
    ELSE
       irec0=0
       IF (paral%io_parent) WRITE(nw,iostat=ierror) irec0,fpos_null
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE w_linres
  ! ==================================================================
  SUBROUTINE r_linres(nr,tag)
    ! ==--------------------------------------------------------------==



    INTEGER                                  :: nr
    CHARACTER(len=2)                         :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'r_linres'

    INTEGER                                  :: i, icompr, ierr, irec0, &
                                                irecord, j, k, lca, len, &
                                                len1, licmp, mtda, mtda0, &
                                                ngwx, nolrx, nread, nuax, &
                                                nubx, nx
    INTEGER(int_8)                           :: fpos
    INTEGER, ALLOCATABLE                     :: icmp(:), mapw(:)
    LOGICAL                                  :: ostda, ostdax
    REAL(real_8)                             :: scale
    REAL(real_8), ALLOCATABLE                :: ca(:), cr(:)

    IF (td03%tda) THEN
       mtda=1
    ELSE
       mtda=2
    ENDIF
    IF ( lractive ) THEN
       IF (td03%tda.AND.td01%msubta.NE.0.AND..NOT.td03%tdacanon) THEN
          ostda=.TRUE.
       ELSE
          ostda=.FALSE.
       ENDIF
       IF (paral%io_parent) THEN
          READ(nr) irecord,fpos
          READ(nr) ngwx,nolrx,nlinr,mtda0,icompr,scale
          IF (mtda.NE.mtda0)&
               CALL stopgm("R_LINRES","cntl%tddft APPROXIMATION DIFFERS",& 
               __LINE__,__FILE__)
          IF (ngwx.NE.spar%ngws)&
               CALL stopgm("R_LINRES","NUMBER OF PW DIFFERS",& 
               __LINE__,__FILE__)
          IF (nolrx.NE.nolr)&
               CALL stopgm("R_LINRES","NUMBER OF STATES DIFFERS",& 
               __LINE__,__FILE__)
          READ(nr) ostdax,nuax,nubx
       ENDIF
       CALL mp_bcast(ngwx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nolrx,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nlinr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(icompr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(scale,parai%io_source,parai%cp_grp)
       CALL mp_bcast(ostdax,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nuax,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nubx,parai%io_source,parai%cp_grp)
       IF (paral%io_parent) THEN
          len=2*nkpt%ngwk
          ALLOCATE(cr(4*len),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          len1=2*ncpw%ngw + 1
          ALLOCATE(mapw(len1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          lca=2*spar%ngws+100
          ALLOCATE(ca(lca),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(ca)!,lca)
          ! Warning: ICMP is allocated as INTEGER*8
          licmp=2*spar%ngws+100
          ALLOCATE(icmp(2*licmp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! ICMP is integer*8 (length LICMP).
          CALL zeroing(icmp)!,2*licmp)
       ELSE
          ! Avoid 'not allocated' runtime errors
          ALLOCATE(ca(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(icmp(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(cr(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(mapw(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
       nread=MIN(nlinw,nlinr)
       DO j=1,mtda
          DO k=1,nread
             DO i=1,nolr
                IF (tag.EQ."WF") THEN
                   CALL rdwfn(nr,clrwf(1,i,k,j),ca,cr,icmp,mapw,&
                        spar%ngws,icompr,scale,.TRUE.,.FALSE.)
                ELSEIF (tag.EQ."WV") THEN
                   CALL rdwfn(nr,clrv(1,i,k,j),ca,cr,icmp,mapw,spar%ngws,icompr,&
                        scale,.TRUE.,.FALSE.)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       nx=nlinr-nread
       IF (paral%io_parent) THEN
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,2*mtda*nx*nolr)
          ENDIF
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
       IF (ostda.AND..NOT.ostdax) THEN
          CALL zeroing(urot)!,nua*nub*nlinw)
       ELSEIF (.NOT.ostda.AND.ostdax) THEN
          IF (paral%io_parent) THEN
             DO i=1,nlinr
                READ(nr)
             ENDDO
          ENDIF
       ELSEIF (ostda.AND.ostdax) THEN
          IF (paral%io_parent) THEN
             IF (nuax.NE.nua)&
                  CALL stopgm("R_LINRES","NUMBER OF OCC ORBITALS DIFFER",& 
                  __LINE__,__FILE__)
             IF (nubx.NE.nub)&
                  CALL stopgm("R_LINRES","NUMBER OF ACT ORBITALS DIFFER",& 
                  __LINE__,__FILE__)
             CALL zeroing(urot)!,nua*nub*nlinw)
             nread=MIN(nlinw,nlinr)
             DO k=1,nread
                READ(nr) ((urot(i,j,k),i=1,nua),j=1,nub)
             ENDDO
          ENDIF
          CALL mp_bcast(urot,nua*nub*nlinw,parai%io_source,parai%cp_grp)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          READ(nr) irec0,fpos
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos)
          ELSE
             CALL fskip(nr,irec0)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE r_linres
  ! ==================================================================

END MODULE rw_linres_utils
