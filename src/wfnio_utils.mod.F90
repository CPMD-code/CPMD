MODULE wfnio_utils
  USE bsym,                            ONLY: bsclcs,&
                                             bsfac,&
                                             hsdown,&
                                             hsup,&
                                             restbs
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gs_ortho,&
                                             gs_ortho_c
  USE io_utils,                        ONLY: file_seek
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE randtowf_utils,                  ONLY: randtowf
  USE readmod,                         ONLY: israndom,&
                                             isrlength
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: fskip,&
                                             izamax,&
                                             zgive
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: w_wfnio
  PUBLIC :: r_wfnio
  PUBLIC :: queryrands
  PUBLIC :: hasrands
  !public :: initrandp
  !public :: putrandp

CONTAINS

  ! ==================================================================
  SUBROUTINE w_wfnio(nw,ierror,nstate,c,tau0,tag_in,irecord)
    ! ==--------------------------------------------------------------==
    ! == WRITE THE WAVEFUNCTIONS IN RESTART FILE                      ==
    ! ==--------------------------------------------------------------==
    !!use wr30wfn_utils, only: wrwfns
    INTEGER                                  :: nw, ierror, nstate
    COMPLEX(real_8)                          :: c(nkpt%ngwk,*)
    REAL(real_8)                             :: tau0(:,:,:)
    CHARACTER(len=*)                         :: tag_in
    INTEGER                                  :: irecord

    CHARACTER(*), PARAMETER                  :: procedureN = 'w_wfnio'

    INTEGER :: first_state, icompr, ierr, ikpt, isca, isub, iwfntype, kbeg, &
      kend, kinc, len, len1, nkpoint, nstate_to_write
    CHARACTER(len=10)                        :: tag
    INTEGER(int_8)                           :: pos
    INTEGER, ALLOCATABLE                     :: icmp(:), mapw(:)
    REAL(real_8)                             :: scale, scale0
    REAL(real_8), ALLOCATABLE                :: ca(:), cr(:)

    CALL tiset(procedureN,isub)
    pos = -1
    irecord = 0
    tag = tag_in
    ! ..the type of wfn output
    iwfntype = 1
    ! CB: BS calculation needs to store two WFs
    IF (cntl%bsymm) THEN
       iwfntype=9
    ENDIF
    ! 
    IF (cnti%wcompb.LT.0) THEN

       icompr=cnti%wcompb
       CALL pwtoao(nw,ierror,nstate,nkpt%nkpnt,c,tau0,icompr,irecord)

    ELSE
       IF (paral%parent) THEN
          len=2*nkpt%ngwk
          len1=ncpw%ngw + 1
          ALLOCATE(cr(4*len),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(mapw(2*len1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)! TODO check length: previously mapw was allocated as *8
          ALLOCATE(ca(2*(spar%ngwks+MOD(spar%ngwks,8))),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)! mod(ngwks,8) is used in COMPRESS
          ALLOCATE(icmp(2*(2*spar%ngwks+100)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          icompr=cnti%wcompb
          IF (icompr.LT.1) THEN
             icompr=1
          ENDIF
       ELSE
          ! To avoid Fortran runtime error 'not allocated'
          ALLOCATE(cr(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ca(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(mapw(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(icmp(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL zeroing(cr)
       CALL zeroing(ca)
       CALL zeroing(mapw)
       CALL zeroing(icmp)
       IF (tag(1:3).EQ.'NIL') THEN
          irecord=-(2*nstate+1)
          !>>>
          !vw isca=idamax(2*nstate*nkpt%ngwk,c,1)
          !vw scale=abs(dgive(c(1,1),isca))
          isca=izamax(nstate*nkpt%ngwk,c,1)
          scale=ABS(zgive(c(1,1),isca))
          !<<<
          CALL mp_max(scale,parai%allgrp)
          IF (iwfntype.EQ.1) THEN
             IF (paral%io_parent) THEN
                WRITE(unit=nw,iostat=ierror) irecord,pos
                WRITE(unit=nw,iostat=ierror) iwfntype
                WRITE(unit=nw,iostat=ierror) nstate,spar%ngwks,icompr,scale
                WRITE(unit=nw,iostat=ierror) cntl%tlsd,spin_mod%nsup,spin_mod%nsdown
             ENDIF
             CALL wrwfns(nw,nstate,ierror,c,ca,cr,icmp,mapw,&
                  icompr,scale)
          ELSE
             CALL stopgm(procedureN,'2 NOT IMPLEMENTED',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          irecord=-(2*nstate*nkpt%nkpts*bsfac+3)
          scale=0._real_8
          CALL inq_swap(kbeg,kend,kinc)
          DO ikpt=kbeg,kend,kinc
             nkpoint=nkpbl(ikpt)
             IF (tkpts%tkblock) THEN
                CALL rkpt_swap(c,nstate,ikpt,tag)
             ENDIF
             !vw>>>
             !vw isca=idamax(2*nkpt%ngwk*nstate*nkpoint,c,1)
             !vw scale0=abs(dgive(c(1,1),isca))
             isca=izamax(nkpt%ngwk*nstate*nkpoint,c,1)
             scale0=ABS(zgive(c(1,1),isca))
             !vw<<<
             IF (scale0.GT.scale) THEN
                scale=scale0
             ENDIF
          ENDDO
          CALL mp_max(scale,parai%allgrp)
          IF (iwfntype.EQ.1) THEN
             IF (paral%io_parent) THEN
                WRITE(unit=nw,iostat=ierror) irecord,pos
                WRITE(unit=nw,iostat=ierror) iwfntype
                WRITE(unit=nw,iostat=ierror)&
                     nstate*nkpt%nkpts,spar%ngwks,icompr,scale
                WRITE(unit=nw,iostat=ierror) cntl%tlsd,spin_mod%nsup,spin_mod%nsdown
             ENDIF
             first_state=1
             DO ikpt=1,nkpt%nblkp
                IF (tkpts%tkblock) THEN
                   CALL rkpt_swap(c,nstate,ikpt,tag)
                ENDIF
                nstate_to_write=nstate*nkpbl(ikpt)
                CALL wrwfns(nw,nstate_to_write,ierror,&
                     c(1,first_state),ca,cr,icmp,mapw,icompr,scale)
                first_state=first_state+nstate_to_write
             ENDDO
             ! CB: Writing the two BS wavefunctions (BS, HS state)            
          ELSE IF (iwfntype.EQ.9) THEN
             IF (paral%io_parent) THEN
                WRITE(unit=nw,iostat=ierror) irecord,pos
                WRITE(unit=nw,iostat=ierror) iwfntype
                WRITE(unit=nw,iostat=ierror)&
                     nstate*nkpt%nkpts,spar%ngwks,icompr,scale
                WRITE(unit=nw,iostat=ierror)&
                     cntl%tlsd,INT(nstate/2),INT(nstate/2),hsup,hsdown
             ENDIF
             first_state=1
             DO ikpt=1,nkpt%nblkp
                IF (tkpts%tkblock) THEN
                   CALL rkpt_swap(c,nstate,ikpt,tag)
                ENDIF
                nstate_to_write=nstate*nkpbl(ikpt)*bsfac
                CALL wrwfns(nw,nstate*nkpbl(ikpt)*bsfac,ierror,&
                     c(1,first_state),ca,cr,icmp,mapw,icompr,scale)
                first_state=first_state+nstate_to_write
             ENDDO
          ELSE
             CALL stopgm(procedureN,'3 NOT IMPLEMENTED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF (paral%parent) THEN
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

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE w_wfnio
  ! ==================================================================
  SUBROUTINE r_wfnio(nr,c,nstate,info,&
       tkpnt0,nstate0,nkpts0,tag_in,fpos_end)
    ! ==--------------------------------------------------------------==
    ! == Read wavefunctions                                           ==
    ! == INPUT:                                                       ==
    ! ==   NR      Logical number of the file                         ==
    ! ==   NSTATE  Number of states                                   ==
    ! ==   TKPNT0  .TRUE. (has k points)                              ==
    ! ==   NSTATE0 Number of states  for the RESTART file             ==
    ! ==   NKPTS0  Number of kpoints for the RESTART file             ==
    ! ==   TAG_IN  'NIL' wavefunctions                                ==
    ! ==           'C0'  electronic wavefunctions                     ==
    ! == OUTPUT:                                                      ==
    ! ==   C      Wavefunctions                                       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, nstate
    COMPLEX(real_8)                          :: c(nkpt%ngwk,nstate,nkpt%nkpnt)
    INTEGER                                  :: info
    LOGICAL                                  :: tkpnt0
    INTEGER                                  :: nstate0, nkpts0
    CHARACTER(len=*)                         :: tag_in
    INTEGER(int_8)                           :: fpos_end

    CHARACTER(*), PARAMETER                  :: procedureN = 'r_wfnio'

    INTEGER :: i, icompr, ierr, ikind, ikpt, ikpt0, isub, iwfntype, lca, len, &
      len1, licmp, n0, na, nah, nb, nbh, ngwa, ngwks0, noortho, northo, &
      nread, nx
    COMPLEX(real_8), ALLOCATABLE             :: smat(:)
    CHARACTER(len=10)                        :: tag
    INTEGER, ALLOCATABLE                     :: icmp(:), mapw(:)
    LOGICAL                                  :: l_tlsd
    REAL(real_8)                             :: scale
    REAL(real_8), ALLOCATABLE                :: ca(:), cr(:)

!real(real_8), target :: c(2*nkpt%ngwk,nstate,nkpt%nkpnt)
! Variables
!complex(real_8), pointer :: c_complex(:,:,:)
! CB:Renamed from TLOG to avoid confusion w/ USPP
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
    ! ==--------------------------------------------------------------==
    tag = tag_in
    IF (tag(1:3).NE.'NIL') THEN
       CALL initrandp(nstate)
    ENDIF
    ! ==--------------------------------------------------------------==
    nread=0
    info=0
    IF (paral%io_parent) READ(nr) iwfntype
    CALL mp_bcast(iwfntype,parai%io_source,parai%cp_grp)
    IF (iwfntype.EQ.1.OR.iwfntype.EQ.9) THEN
       IF (paral%io_parent) THEN
          READ(nr) n0,ngwks0,icompr,scale
          IF (iwfntype.EQ.1) THEN
             READ(nr) l_tlsd,na,nb
          ELSE
             READ(nr) l_tlsd,na,nb,nah,nbh
             restbs=.FALSE.
             IF (bsclcs.EQ.2) THEN
                na=nah
                nb=nbh
             ENDIF
          ENDIF
       ENDIF
       CALL mp_bcast(n0,parai%io_source,parai%cp_grp)
       CALL mp_bcast(ngwks0,parai%io_source,parai%cp_grp)
       CALL mp_bcast(icompr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(l_tlsd,parai%io_source,parai%cp_grp)
       CALL mp_bcast(na,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nb,parai%io_source,parai%cp_grp)
       CALL mp_bcast(restbs,parai%io_source,parai%cp_grp)
       CALL mp_bcast(scale,parai%io_source,parai%cp_grp)
       cnti%rcompb=icompr
       IF (cnti%rcompb.EQ.-1) THEN
          CALL stopgm(procedureN,'1 NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
          GOTO 400
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          len=2*nkpt%ngwk
          ALLOCATE(cr(4*len),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          len1=ncpw%ngw + 1
          ALLOCATE(mapw(2*len1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ngwa=MAX(ngwks0,spar%ngwks)
          lca=MAX(2*ngwa+100,nstate0*MAX(nstate0,nstate))
          ALLOCATE(ca(lca),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(ca)!,lca)
          ! Warning: ICMP is allocated as INTEGER*8
          licmp=MAX(2*ngwa+100,(nstate-nstate0)/2+1)
          ALLOCATE(icmp(2*licmp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)! TODO maybe define ICMP as I*8?
          ! ICMP is integer*8 (length LICMP).
          CALL zeroing(icmp)!,2*licmp)
       ELSE
          ! avoid fortran runtime error 'not allocated'
          ALLOCATE(cr(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(mapw(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(ca(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(icmp(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF

       ALLOCATE(smat(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ! ==--------------------------------------------------------------==
       IF (tag(1:3).EQ.'NIL') THEN
          nread=0
          CALL rdwfns(nr,MIN(nstate,n0),c,ca,cr,icmp,mapw,&
               ngwks0,icompr,scale,tkpnt0,nread)
          nx=nstate-nread
          IF(nx.GT.0) CALL repprngu_vec_cmplx(nkpt%ngwk*nx,c(1,nread+1,1))
          GOTO 200
       ENDIF
       ! ==--------------------------------------------------------------==
       ! == Case for C0 (electronic wavefunctions)                       ==
       ! ==--------------------------------------------------------------==
       IF (tkpts%tkpnt.OR.tkpnt0) THEN
          ! Check consistency of N0.
          IF (n0.NE.nstate0*nkpts0) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) procedureN//'! N0=',N0,' NSTATE0=',NSTATE0,&
                  'NKPST0=',nkpts0
             CALL stopgm(procedureN,&
                  'THE NUMBER OF WAVEFUNCTIONS IS NOT CORRECT',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (iwfntype.EQ.9) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) procedureN//'! Cannot read BS RESTART w/ k-points.'
             CALL stopgm(procedureN,&
                  'BS RESTART cannot be read in k-point run.',& 
                  __LINE__,__FILE__)
          ENDIF
          nread=0
          ikpt0=0
          DO ikpt=1,nkpt%nblkp
             DO ikind=1,nkpbl(ikpt)
                ikpt0=ikpt0+1
                IF (ikpt0.LE.nkpts0) THEN
                   ! We load a new set of nstates for a new kpoint.
                   noortho=0
                   DO i=1,nstate
                      IF (i.LE.nstate0) THEN
                         ! We load a new state.
                         nread=nread+1
                         CALL rdwfn(nr,c(1,i,ikind),ca,cr,icmp,mapw,ngwks0,&
                              icompr,scale,tkpnt0,.FALSE.)
                         ! We adapt the format of C0.
                         IF (tkpts%tkpnt) THEN
                            IF (.NOT.tkpnt0) THEN
                               c(ncpw%ngw+1:2*ncpw%ngw,i,ikind)=CONJG(c(1:ncpw%ngw,i,ikind))
                               !call dcopy(2*ngw,c(1      ,i,ikind),1,&
                               !     c(1+2*ngw,i,ikind),1)
                               !call dscal(ngw,-1._real_8,c(2+2*ngw,i,ikind),2)
                               IF (geq0) THEN
                                  !c(2*ngw+1,i,ikind)=0._real_8
                                  !c(2*ngw+2,i,ikind)=0._real_8
                                  c(ncpw%ngw+1,i,ikind)=(0.0_real_8,0.0_real_8)
                               ENDIF
                            ENDIF
                            !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                            CALL c_clean(c(:,i,ikind),1,ikind)
                         ENDIF
                      ELSE
                         ! Use a randomized state.
                         !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                         CALL randtowf(c(1:,i:,ikind),1,ikind,ikpt0)
                         CALL putrandp(i)
                      ENDIF
                      ! C has to be real.
                      IF (tkpnt0.AND..NOT.tkpts%tkpnt) THEN
                         IF (cntl%tlsd) THEN
                            CALL putconj(c,spin_mod%nsup)
                            CALL putconj(c(1,spin_mod%nsup+1,1),spin_mod%nsdown)
                         ELSE
                            CALL putconj(c,nstate)
                         ENDIF
                      ENDIF
                   ENDDO
                   IF (nstate.GT.nstate0) THEN
                      noortho=nstate-nstate0
                      ! We orthogonalize.
                      IF (tkpts%tkpnt) THEN
                         IF (cntl%tlsd) THEN
                            IF (spin_mod%nsup.GT.nstate0) THEN
                               CALL gs_ortho_c(c(1,1,ikind),nstate0,&
                                    c(1,nstate0+1,ikind),spin_mod%nsup,smat)
                               ! We copy UP states  to DOWN states.
                               CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                                    c(1,spin_mod%nsup+1,ikind),1)
                            ELSEIF (spin_mod%nsup.EQ.nstate0) THEN
                               ! We copy UP states  to DOWN states.
                               CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                                    c(1,spin_mod%nsup+1,ikind),1)
                            ELSE
                               ! We orthogonalize only the down states.
                               northo=nstate0-spin_mod%nsup
                               CALL gs_ortho_c(c(1,spin_mod%nsup+1,ikind),northo,&
                                    c(1,spin_mod%nsup+northo+1,ikind),noortho,smat)
                            ENDIF
                         ELSE
                            CALL gs_ortho_c(c(1,1,ikind),nstate0,&
                                 c(1,nstate0+1,ikind),noortho,smat)
                         ENDIF
                      ELSE
                         IF (cntl%tlsd) THEN
                            IF (spin_mod%nsup.GT.nstate0) THEN
                               CALL gs_ortho(c(1,1,ikind),nstate0,&
                                    c(1,nstate0+1,ikind),spin_mod%nsup)
                               ! We copy UP states  to DOWN states.
                               CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                                    c(1,spin_mod%nsup+1,ikind),1)
                            ELSEIF (spin_mod%nsup.EQ.nstate0) THEN
                               ! We copy UP states  to DOWN states.
                               CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                                    c(1,spin_mod%nsup+1,ikind),1)
                            ELSE
                               ! We orthogonalize only the down states.
                               northo=nstate0-spin_mod%nsup
                               CALL gs_ortho(c(1,spin_mod%nsup+1,ikind),northo,&
                                    c(1,spin_mod%nsup+northo+1,ikind),noortho)
                            ENDIF
                         ELSE
                            CALL gs_ortho(c(1,1,ikind),nstate0,&
                                 c(1,nstate0+1,ikind),noortho)
                         ENDIF
                      ENDIF ! IF(TKPNT)
                   ELSEIF (paral%io_parent) THEN
                      ! NSTATE0.GE.NSTATE
                      nx=nstate0-nstate
                      DO i=1,nx
                         nread=nread+1
                         READ(nr)
                         READ(nr)
                      ENDDO
                   ENDIF     ! IF(NSTATE.GT.NSTATE0)
                ELSE
                   ! We have more k points -> we use the set of nstate
                   ! for ikind=1.
                   IF (ikind.NE.1) THEN
                      CALL dcopy(2*nkpt%ngwk*nstate,c(1,1,1),1,c(1,1,ikind),1)
                   ELSE
                      ! All the block has to be copied. Do nothing.
                      GOTO 100
                   ENDIF
                ENDIF         ! IF(IKPT0.LE.KPTS0)
             ENDDO             ! DO IKIND
100          CONTINUE
             IF (tkpts%tkblock) CALL wkpt_swap(c,nstate,ikpt,tag)
          ENDDO                 ! DO IKPT
       ELSE
          ! ..CASE with NO k-points (in read and write)
          nread=0
          IF (cntl%tlsd) THEN
             ! McB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             IF ( cntl%tshop ) THEN
                DO i=1,nstate
                   IF (i.LE.nstate0) THEN
                      ! We load a new state.
                      nread=nread+1
                      CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,&
                           icompr,scale,tkpnt0,.FALSE.)
                   ELSE
                      ! Use a randomized state.
                      !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                      CALL randtowf(c(1:,i:,1),1,1,1)
                      CALL putrandp(i)
                      CALL gs_ortho(c(1,1,1),i-1,c(1,i,1),1)
                   ENDIF
                ENDDO
                ! McB      surface hopping currently needs to read as if cntl%tlsd=L_TLSD=.F.
                ! mb-bug               GO TO 500
                ! McB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! mb-bug            ENDIF
             ELSE
                ! McB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF (l_tlsd) THEN
                   ! --> LSD:LSD
                   ! CB: We need to skip the BS wavefunction if it is not used.
                   IF (iwfntype.EQ.9.AND.bsclcs.EQ.2.AND.paral%io_parent)&
                        CALL fskip(nr,2*n0)
                   DO i=1,spin_mod%nsup
                      IF (i.LE.na) THEN
                         ! We load a new state.
                         nread=nread+1
                         CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,&
                              icompr,scale,tkpnt0,.FALSE.)
                      ELSE
                         ! Use a randomized state.
                         !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                         CALL randtowf(c(1:,i:,1),1,1,0)
                         CALL putrandp(i)
                         CALL gs_ortho(c(1,1,1),i-1,c(1,i,1),1)
                      ENDIF
                   ENDDO
                   IF (nread.LT.na.AND.paral%io_parent) THEN
                      nx=na-nread
                      DO i=1,nx
                         nread=nread+1
                         READ(nr)
                         READ(nr)
                      ENDDO
                   ENDIF
                   DO i=1,spin_mod%nsdown
                      IF (i.LE.nb) THEN
                         ! We load a new state.
                         nread=nread+1
                         CALL rdwfn(nr,c(1,spin_mod%nsup+i,1),ca,cr,icmp,mapw,ngwks0,&
                              icompr,scale,tkpnt0,.FALSE.)
                      ELSE
                         ! Use a randomized state.
                         !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                         CALL randtowf(c(1:,spin_mod%nsup+i:,1),1,1,0)
                         CALL putrandp(spin_mod%nsup+i)
                         CALL gs_ortho(c(1,spin_mod%nsup+1,1),i-1,c(1,spin_mod%nsup+i,1),1)
                      ENDIF
                   ENDDO
                   IF (nread.LT.(na+nb).AND.paral%io_parent) THEN
                      nx=na+nb-nread
                      DO i=1,nx
                         nread=nread+1
                         READ(nr)
                         READ(nr)
                      ENDDO
                   ENDIF
                   ! CB: We need to skip the HS wavefunction if it was not used.
                   IF (iwfntype.EQ.9.AND.bsclcs.EQ.1.AND.paral%io_parent)&
                        CALL fskip(nr,2*n0)
                ELSE
                   ! --> LSD:LDA
                   DO i=1,MAX(spin_mod%nsup,spin_mod%nsdown)
                      IF (i.LE.nstate0) THEN
                         ! We load a new state.
                         IF (i.LE.spin_mod%nsup) THEN
                            nread=nread+1
                            CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,&
                                 icompr,scale,tkpnt0,.FALSE.)
                            IF (i.LE.spin_mod%nsdown)&
                                 CALL dcopy(2*nkpt%ngwk,c(1,i,1),1,c(1,spin_mod%nsup+i,1),1)
                         ELSE
                            nread=nread+1
                            CALL rdwfn(nr,c(1,spin_mod%nsup+i,1),ca,cr,icmp,mapw,ngwks0,&
                                 icompr,scale,tkpnt0,.FALSE.)
                         ENDIF
                      ELSE
                         ! Use a randomized state.
                         IF (i.LE.spin_mod%nsup) THEN
                            !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                            CALL randtowf(c(1:,i:,1),1,1,1)
                            CALL putrandp(i)
                            CALL gs_ortho(c(1,1,1),i-1,c(1,i,1),1)
                         ELSE
                            !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                            CALL randtowf(c(1:,spin_mod%nsup+i:,1),1,1,1)
                            CALL putrandp(spin_mod%nsup+i)
                            CALL gs_ortho(c(1,spin_mod%nsup+1,1),i-1,c(1,spin_mod%nsup+i,1),1)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
                ! mb!!!!!!!!!!!!!!!!!!!!!!
             ENDIF
             ! mb!!!!!!!!!!!!!!!!!!!!!!
          ELSE
             IF (l_tlsd) THEN
                ! --> LDA:LSD (read only alpha states)
                DO i=1,nstate
                   IF (i.LE.na) THEN
                      ! We load a new state.
                      nread=nread+1
                      CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,&
                           icompr,scale,tkpnt0,.FALSE.)
                   ELSE
                      ! Use a randomized state.
                      !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                      CALL randtowf(c(1:,i:,1),1,1,1)
                      CALL putrandp(i)
                      CALL gs_ortho(c(1,1,1),i-1,c(1,i,1),1)
                   ENDIF
                ENDDO
             ELSE
                ! --> LDA:LDA
                ! McB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! McB          surface hopping currently needs to read as if cntl%tlsd=.F.
                ! mb-bug  500         continue
                ! McB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CALL rdwfns(nr,MIN(nstate0,nstate),c,ca,cr,icmp,mapw,&
                     ngwks0,icompr,scale,tkpnt0,nread)
                DO i=nstate0+1,nstate
                   !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                   CALL randtowf(c(1:,i:,1),1,1,1)
                   CALL putrandp(i)
                   CALL gs_ortho(c(1,1,1),i-1,c(1,i,1),1)
                ENDDO
             ENDIF
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
200    CONTINUE
       IF (n0.GT.nread.AND.paral%io_parent) THEN
          nx=n0-nread
          IF (cntl%is_in_stream) THEN
             CALL file_seek(nr,.FALSE.,fpos_end)
          ELSE
             CALL fskip(nr,2*nx)
          ENDIF
       ENDIF
       IF (paral%parent) THEN
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
       ELSE
          DEALLOCATE(ca,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(icmp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(cr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(mapw,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF

       DEALLOCATE(smat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       ! ==--------------------------------------------------------------==
400    CONTINUE
       ! distribute the wavefunction along the groups
       IF (parai%cp_nogrp.GT.1) THEN
          CALL mp_bcast(c,SIZE(c),0,parai%cp_inter_grp)
       ENDIF
    ELSE
       CALL stopgm(procedureN,'2 NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE r_wfnio
  ! ==================================================================
  SUBROUTINE initrandp(nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'initrandp'

    INTEGER                                  :: ierr, len

! ==--------------------------------------------------------------==

    IF (nstate.GT.isrlength) THEN
       IF (isrlength.NE.0) THEN
          DEALLOCATE(israndom,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       len = nstate * 4
       ALLOCATE(israndom(len),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       isrlength = nstate
    ENDIF
    CALL zeroing(israndom)!,isrlength)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initrandp
  ! ==================================================================
  SUBROUTINE putrandp(ns)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns

! ==--------------------------------------------------------------==

    IF (ns.GT.isrlength) THEN
       CALL stopgm('PUTRANDP','Out of range',& 
            __LINE__,__FILE__)
    ELSE
       israndom(ns)=1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putrandp
  ! ==================================================================
  LOGICAL FUNCTION hasrands()
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    INTEGER :: i
    ! ==--------------------------------------------------------------==
    hasrands = .FALSE.
    DO i=1,isrlength
       hasrands = hasrands .OR. israndom(i).NE.0
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION hasrands
  ! ==================================================================
  LOGICAL FUNCTION queryrands(ns)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    INTEGER :: ns
    ! ==--------------------------------------------------------------==
    queryrands = israndom(ns).NE.0
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION queryrands
  ! ==================================================================

END MODULE wfnio_utils
