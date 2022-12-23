MODULE mdpt_utils
  USE atwf,                            ONLY: atwp
  USE bsym,                            ONLY: bsfac
  USE cl_init_utils,                   ONLY: cl_init
  USE clas,                            ONLY: tclas
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE fusion_utils,                    ONLY: fusion,&
                                             separate
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: &
       c1k, ck, clrwf, nlinr, nlinw, nolr, nua, nub, td01, td03, tshl, urot
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE md_driver,                       ONLY: mddiag
  USE mdclas_utils,                    ONLY: mdclas
  USE mdfile_utils,                    ONLY: mdfile
  USE mdmain_utils,                    ONLY: mdmain
  USE mdshop_bo_utils,                 ONLY: mdshop_bo
  USE mdshop_cp_utils,                 ONLY: mdshop_cp
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_qm,&
                                             mm_revert
  USE mm_mddiag_utils,                 ONLY: mm_mddiag
  USE mm_mdmain_utils,                 ONLY: mm_mdmain
  USE mm_mdshop_bo_utils,              ONLY: mm_mdshop_bo
  USE mm_mdshop_cp_utils,              ONLY: mm_mdshop_cp
  USE mw,                              ONLY: tmw
  USE nabdy_md,                        ONLY: nabdy_dyn
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE shop,                            ONLY: sh02
  USE shop_rest_2,                     ONLY: c0old,&
                                             c0old2
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mdpt

CONTAINS


  ! ==================================================================
  SUBROUTINE mdpt
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'mdpt'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:,:), cm(:), &
                                                sc0(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER :: i, ierr, isub, iunit, mm, ms, nact, nc0, nc0_1d, nc0_2d, &
      nc0_3d, nc1, nc2, nc2_1d, nc2_2d, nc2_3d, ncm, nmin, nsc0, nstate, &
      nstx, nus, nxx
    LOGICAL                                  :: statusdummy
    REAL(real_8), ALLOCATABLE                :: gamx(:), gamy(:), vpp(:)

! TODO refactor this 
! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==  C0         WF POSITIONS                                     ==
! ==  CM         WF VELOCITIES                                    ==
! ==  C2         WF FORCES                                        ==
! ==  SC0        S*C0 (only Vanderbilt)                           ==
! ==  GAMX       CONSTRAINT MATRIX POSITIONS                      ==
! ==  GAMY       CONSTRAINT MATRIX VELOCITIES                     ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    IF ( paral%qmnode ) THEN

       IF (cntl%tddft.AND.(.NOT.tshl%txfmqc)) THEN
          CALL lr_in
          CALL tddft_input
          cntl%tsymrho=.FALSE.
       ENDIF

       nstate=crge%n
       ! McB
       IF ( cntl%tshop ) nstate=sh02%nst_s0+sh02%nst_s1
       ! McB
       nc1=1
       nua=1
       nub=1
       nlinw=0
       nlinr=0
       nus=1
       IF (cntl%tdiag) THEN
          IF (cntl%tddft) CALL stopgm('MDPT',&
               'CNTL%TDDFT.AND.TDIAG NOT POSSIBLE',& 
               __LINE__,__FILE__)
          IF (cntl%tdavi) THEN
             IF (cnti%ndavv.EQ.0) THEN
                cnti%ndavv=2*nstate
             ELSE
                cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
             ENDIF
             nstate=cnti%ndavv
          ENDIF
          nc0_1d = 2 * nkpt%ngwk 
          nc0_2d = nstate
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          IF (cntl%tfrsblk)THEN
             ncm=2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
          ELSEIF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
                  ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
          ELSE
             ncm=2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
          ENDIF
          nc2_1d = nkpt%ngwk
          nc2_2d = MAX(atwp%numaormax,nstate)
          nc2_3d = 1
          nc2 = nc2_1d * nc2_2d * nc2_3d
          !nc2=2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
          IF (cntl%tdavi)ncm=ncm*2
          ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! EHR[   NABDY[
       ELSEIF (cntl%tmdbo.OR.cntl%tmdfile.OR.cntl%tmdeh.OR.cntl%tnabdy) THEN
          ! EHR]   NABDY]
          IF (cntl%tsde) THEN
             ncm=8
          ELSE IF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
                  ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4+100
          ELSE IF (cntl%pcg) THEN
             ncm=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          ELSEIF (cntl%tmdeh) THEN
             ncm=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          ENDIF
          ! CB: Need memory for two wf in BS case          
          nc0_1d = nkpt%ngwk !<<<<<< 2 * nkpt%ngwk 
          nc0_2d = nstate * bsfac
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          nc2_1d = nkpt%ngwk
          nc2_2d = nstate
          nc2_3d = bsfac
          nc2 = nc2_1d * nc2_2d * nc2_3d
          !nc2=(2*nkpt%ngwk*nstate+8)*bsfac
          ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (cntl%tddft) THEN
             DO i=1,nstate
                IF (crge%f(i,1).LT.1.e-6_real_8)&
                     CALL stopgm('TDDFT','ALL STATES HAVE TO BE OCCUPIED',& 
                     __LINE__,__FILE__)
             ENDDO
             ms=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri)
             IF (td03%tda.AND.td01%ldiag.EQ.1) THEN
                mm=1
             ELSE
                mm=2
             ENDIF
             IF (td03%tda.AND.td01%msubta.GT.0) THEN
                nact=td01%msubta
             ELSE
                nact=nstate
             ENDIF
             nolr=nact
             nc1=2*ncpw%ngw*nact*ms*mm
             nua=nstate
             nub=nact
             nlinw=ms
             nlinr=0
             nus=2*nua*nub
             ! ..TSH[
             IF (td03%tdlz.OR.tshl%tdtully) THEN
                nc0_1d = nkpt%ngwk  !vw do we realy need this factor of !!! 2 * 
                nc0_2d = 2*(nstate+ms)
                nc0_3d = nkpt%nkpnt
                nc0=nc0_1d*nc0_2d*nc0_3d
                nc2_1d = nkpt%ngwk
                nc2_2d = 2*(nstate+ms)
                nc2_3d = nkpt%nkpnt
                nc2 = nc2_1d * nc2_2d * nc2_3d
                !nc2=2*nkpt%ngwk*(2*(nstate+ms))*nkpt%nkpnt+8
             ENDIF
             ! ..TSH]
          ENDIF
       ELSE
          IF (cntl%tddft) CALL stopgm('MDPT',&
               'TDDFT.AND.CP NOT POSSIBLE',& 
               __LINE__,__FILE__)
          ! NN  For the BS-CPMD, HS & BS wavefunctions are propagated at the 
          ! same time, thus the array sizes of C0,CM and C2 are multiplied by 2 
          nc0_1d = nkpt%ngwk !vw do we realy need this factor of !!! 2 *
          nc0_2d = nstate * bsfac
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          ncm=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
          nc2_1d = nkpt%ngwk
          nc2_2d = nstate
          nc2_3d = bsfac
          nc2 = nc2_1d * nc2_2d * nc2_3d
          !nc2=(2*nkpt%ngwk*nstate+8)*bsfac
       ENDIF
       nsc0=1
       IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%quenchb.AND.cntl%pcg.AND.cntl%pcgmin)) THEN
          nsc0=2*nkpt%ngwk*nstate*nkpt%nkpnt*bsfac
       ENDIF


       IF (cntl%tdipd.OR.vdwl%vdwd) THEN
          IF (cntl%tddft) CALL stopgm('MDPT',&
               'TDDFT AND DIPOLE DYNAMICS NOT POSSIBLE',& 
               __LINE__,__FILE__)
          lenbk=nxxfun(nstate)
          nxx=MAX(2*lenbk*parai%nproc,nc2)
          !vw here we need to fix one of the nc2_*d
          nc2_3d=FLOOR(REAL(nxx,real_8)/REAL(nc2_1d*nc2_2d,real_8))
          nc2 = nc2_1d * nc2_2d * nc2_3d
          nsc0=MAX(nxx,nsc0)
       ENDIF
       IF (cntl%tshop) THEN
          nc0=2*nc0
          nc0_2d = 2* nc0_2d
          ncm=2*ncm
          nc2_2d = 2*nc2_2d
          nc2=2*nc2
          ALLOCATE(c0old(nc0,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(c0old2(nc0,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          nmin=1
          ALLOCATE(c0old(nmin,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(c0old2(nmin,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(c0(nc0_1d,nc0_2d,nc0_3d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c0)
       ! ..TSH[
       IF (tshl%tdtully) THEN
          ! TODO determine correct dimensions for CK and C1K
          ALLOCATE(ck(ncpw%ngw,nc0/MAX(ncpw%ngw,1)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(c1k(ncpw%ngw,nstate,nc1/(MAX(ncpw%ngw,1)*nstate)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ..TSH]
       ALLOCATE(cm(ncm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cm)

       IF (nc1 /= 1) THEN
          ALLOCATE(c1(ncpw%ngw,nact,ms,mm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(c1(1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL zeroing(c1)
       ALLOCATE(c2(nc2_1d,nc2_2d,nc2_3d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c2)
       ALLOCATE(sc0(nsc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nus,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(urot)!,nus)
       clrwf => c1
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          ALLOCATE(gamx(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gamy(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(gamx(nstate*nstate*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gamy(nstate*nstate*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! FNL and DFNL
       CALL mm_dim(mm_go_qm,statusdummy)
       ! 
       IF (td03%sh_diag) THEN
          CALL fnlalloc(nstate+ms+2,.TRUE.,cntl%tpres)
       ELSE
          CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
       ENDIF
       CALL mm_dim(mm_revert,statusdummy)
    ELSE
       nmin = 10             ! TODO what is this???
       ! 
       ALLOCATE(c0(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c0)
       ALLOCATE(cm(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cm)

       ALLOCATE(c1(nmin,1,1,1),STAT=ierr) ! TODO understand what is nmin and which dimensions are correct
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c1)

       ALLOCATE(c2(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c2)
       ALLOCATE(sc0(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamx(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamy(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! 
       ALLOCATE(c0old(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0old2(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! 
    ENDIF

    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('      MDPT')
    ! ==--------------------------------------------------------------==

    IF ( cntl%tshop.AND.(cntl%tfusi.OR.cntl%tsep) ) THEN
       ! CALL STOPGM('MDPT','CNTL%TSHOP.AND.(CNTL%TFUSI.OR.TSEP) NOT POSSIBLE')
       IF (cntl%tfusi) THEN
          CALL fusion(c0,cm)
       ELSEIF (cntl%tsep) THEN
          CALL separate(c0,cm)
       ENDIF
       go to 500
    ENDIF

    IF ( cntl%tddft.AND.cntl%tshop ) THEN
       CALL stopgm('MDPT','CNTL%TDDFT.AND.TSHOP NOT POSSIBLE',& 
            __LINE__,__FILE__)
    ENDIF

    IF ( cntl%bsymm.AND.cntl%tshop ) THEN
       CALL stopgm('MDPT','BROKEN SYMMETRY WITH SURFACE HOPPING'&
            // 'NOT POSSIBLE',& 
            __LINE__,__FILE__)
    ENDIF
    !
    IF ((cntl%tshop.OR.cntl%tqmmm.OR.cntl%tddft.OR.cntl%tmdfile.OR.tclas)&
         .AND.tmw.AND.(.NOT.tshl%txfmqc)) &
         CALL stopgm('MDPT',' MULTIPLE WALKER IS NOT '&
         //'IMPLEMENTED FOR THE REQUIRED MD TYPE',& 
         __LINE__,__FILE__)
    !
    IF (tclas) THEN
       IF (tkpts%tkpnt)CALL stopgm('MDPT','K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (cntl%tddft) CALL stopgm('MDPT','TDDFT.AND.TCLAS NOT POSSIBLE',& 
            __LINE__,__FILE__)
       iunit=5
       CALL cl_init(iunit)
       CALL mdclas(c0,cm,c2,sc0,gamx,gamy)
       ! EHR[
    ELSEIF (cntl%tmdbo.OR.cntl%tmdeh) THEN
       ! EHR]
       IF ( cntl%tqmmm ) THEN
          IF ( cntl%tshop ) THEN
             ! CALL STOPGM('MDPT','SURFACE HOPPING WITH BO-MD AND QM/MM '
             ! .                    // 'NOT IMPLEMENTED')
             CALL mm_mdshop_bo(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
          ELSE
             CALL mm_mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
          ENDIF
       ELSE
          IF ( cntl%tshop ) THEN
             ! CALL STOPGM('MDPT','CNTL%TMDBO.AND.TSHOP NOT POSSIBLE')
             CALL mdshop_bo(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
          ELSE
             CALL mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
          ENDIF
       ENDIF
    ELSEIF (cntl%tmdfile) THEN
       CALL mdfile(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
       !NABDY[
    ELSEIF (cntl%tnabdy) THEN
       IF (cntl%tqmmm) THEN
          CALL stopgm('MDPT','NABDY QM/MM NOT IMPLEMENTED',&
               __LINE__,__FILE__)
          !call mm_nabdy_md(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
       ELSE
          CALL nabdy_dyn(c0,cm,c1,c2(:,:,1),sc0,vpp,gamx,gamy)
       ENDIF
       !NABDY]
    ELSE
       IF (tkpts%tkpnt)CALL stopgm('MDPT','K-POINTS NOT IMPLEMENTED',&  
            __LINE__,__FILE__)
       IF (cntl%tddft) CALL stopgm('MDPT','TDDFT.AND.CP NOT POSSIBLE',& 
            __LINE__,__FILE__)
       IF ( cntl%tqmmm ) THEN
          IF ( cntl%tshop ) THEN
             ! CALL STOPGM('MDPT','CNTL%TQMMM.AND.TSHOP NOT POSSIBLE')
             CALL mm_mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
          ELSE
             CALL mm_mdmain(c0,cm,c2,sc0,gamx,gamy)
          ENDIF
       ELSE
          IF ( cntl%tshop ) THEN
             ! CALL STOPGM('MDPT','SURFACE HOPPING NOT POSSIBLE')
             CALL mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
          ELSE
             CALL mdmain(c0,cm,c2,sc0,gamx,gamy)
          ENDIF
       ENDIF
    ENDIF
500 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tshop) THEN
       DEALLOCATE(c0old,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c0old2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF ( paral%qmnode ) THEN
       CALL fnldealloc(.TRUE.,cntl%tpres)
    ELSE
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ..TSH[
    IF (tshl%tdtully) THEN
       DEALLOCATE(ck,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c1k,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ..TSH]
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

  END SUBROUTINE mdpt
  ! ==================================================================

END MODULE mdpt_utils
