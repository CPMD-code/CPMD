MODULE td_nhdav_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr02,&
                                             td01,&
                                             td02,&
                                             td03
  USE lr_force_utils,                  ONLY: give_scr_lr_force,&
                                             lr_force
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_ortho
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE soft,                            ONLY: soft_com
  USE sort_utils,                      ONLY: sort2
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE td_dav_utils,                    ONLY: dotps
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_nhdav
  PUBLIC :: give_scr_td_nhdav

CONTAINS

  ! ==================================================================
  SUBROUTINE td_nhdav(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,&
       teig,nstate,nroot,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==  Non-Hermitian Davidson diagonalisation                      ==
    ! ==  Stratman et al. JCP 1098218 (1998)                         ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(ncpw%ngw,*), &
                                                sc0(*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:), drhoe(:,:), teig(:)
    INTEGER                                  :: nstate, nroot
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,nroot,2)
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_nhdav'
    INTEGER, PARAMETER                       :: mroot = 20, nd = 5*mroot 

    CHARACTER(len=5)                         :: hamb, hapb
    COMPLEX(real_8), ALLOCATABLE             :: cb(:,:,:), cr(:,:,:), &
                                                csm(:,:,:), csp(:,:,:)
    INTEGER                                  :: i, ierr, ig, info, is, istep, &
                                                isub, j, lwork, namat, nconv, &
                                                ndim, ndmax, ndold, nsigma, &
                                                ntogo
    INTEGER, ALLOCATABLE                     :: ind(:)
    REAL(real_8)                             :: e2(5), fee, so, somax, t1, &
                                                t2, tcpu, tx, vp
    REAL(real_8), ALLOCATABLE                :: ammat(:,:), apmat(:,:), &
                                                aux(:), bmmat(:,:), deig(:), &
                                                vpp(:), work(:)
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset(procedureN,isub)
    ndmax = MAX(2*nroot+2,td01%ndavspace)
    ndmax = 2*((ndmax+1)/2)
    IF (nroot.GT.mroot) CALL stopgm("TD_NHDAV","MAX NUMBER OF ROOTS EXCEEDED",& 
         __LINE__,__FILE__)
    IF (ndmax.GT.nd) CALL stopgm("TD_NHDAV","MAX SUBSPACE SIZE EXCEEDED",& 
         __LINE__,__FILE__)
    IF (INDEX(orbital,"CANON").EQ.0) CALL stopgm("TD_DAV","ONLY CANONICAL ORBITALS ALLOWED",& 
         __LINE__,__FILE__)
    ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(csp(ncpw%ngw,nstate,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(csm(ncpw%ngw,nstate,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cb(ncpw%ngw,nstate,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr(ncpw%ngw,nstate,4),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(apmat(nd,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ammat(nd,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bmmat(nd,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(deig(ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ind(ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(apmat)!,nd*ndmax)
    CALL zeroing(ammat)!,nd*ndmax)
    CALL zeroing(bmmat)!,nd*ndmax)
    ! ==--------------------------------------------------------------==
    ALLOCATE(aux(ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    lwork=4*ndmax
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (td03%tda) THEN
       hamb="TDA"
       hapb="TDA"
    ELSE
       hamb="TDAMB"
       hapb="TDAPB"
    ENDIF
    ntogo=td01%ndavmax
    ! ==--------------------------------------------------------------==
    cntl%prec=.TRUE.
    CALL ksdiag(vpp)
    DO ig=1,ncpw%ngw
       vp=vpp(ig)
       vpp(ig)=MAX(vp,lr02%lr_hthrs)
    ENDDO
    ! ==--------------------------------------------------------------==
    IF ( paral%parent .AND. kprint.GE.0 ) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T6,A,T64,"==")')&
            ' NON-HERMITIAN DAVIDSON DIAGONALISATION OF TDDFT MATRIX '
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    ! ..initilize vectors
1000 CONTINUE
    t1=m_walltime()
    CALL dcopy(2*ncpw%ngw*nstate*nroot,c1,1,cb,1)
    ndim=nroot
    IF (.NOT.td03%tda) THEN
       DO i = 1, nroot
          ndim=ndim+1
          IF (paral%parent.AND.ndim.GT.ndmax) THEN
             CALL stopgm("TD_NHDAV","SUBSPACE TOO SMALL",& 
                  __LINE__,__FILE__)
          ENDIF
          vp=dasum(2*ncpw%ngw*nstate,c1(1,1,i,2),1)
          CALL mp_sum(vp,parai%allgrp)
          IF (vp.GT.1.e-10_real_8) THEN
             ! Left eigenvector is available
             CALL dcopy(2*ncpw%ngw*nstate,c1(1,1,i,2),1,cb(1,1,ndim),1)
          ELSE
             ! Generate left eigenvector
             j=ndim-nroot
             CALL lr_ortho(nstate,c0,cb(:,:,j))
             CALL lr_force(c0,cb(:,:,j),cb(:,:,ndim),sc0,e2,&
                  ddxc,psi,drhoe,eigv,&
                  nstate,hapb,orbital)
          ENDIF
       ENDDO
    ENDIF
    DO i = 1, ndim
       CALL lr_ortho(nstate,c0,cb(1,1,i))
       DO j=1,i-1
          so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,j))
          CALL mp_sum(so,parai%allgrp)
          CALL daxpy(2*ncpw%ngw*nstate,-so,cb(:,:,j),1,cb(:,:,i),1)
       ENDDO
       so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,i))
       CALL mp_sum(so,parai%allgrp)
       CALL dscal(2*ncpw%ngw*nstate,1._real_8/SQRT(so),cb(:,:,i),1)
    ENDDO
    t2=m_walltime()
    IF ( paral%parent  .AND. kprint.GE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)') ' Number of states initialized ',ndim
       tcpu=(t2-t1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,F10.2,/)') ' Time for initialization ',tcpu
       IF (paral%io_parent)&
            WRITE(6,'(A,T43,A)')&
            '  ITER      STATES      SUBSPACE','   RESIDUAL        TCPU'
    ENDIF
    nsigma=0
    nconv=0
    DO istep=1,ntogo
       t1=m_walltime()
       ! ..initialize sigma vectors
       DO i = nsigma+1,ndim
          CALL lr_force(c0,cb(:,:,i),csp(:,:,i),sc0,e2,&
               ddxc,psi,drhoe,eigv,&
               nstate,hapb,orbital)
          CALL dscal(2*ncpw%ngw*nstate,-1._real_8,csp(:,:,i),1)
          CALL lr_ortho(nstate,c0,csp(:,:,i))
          CALL lr_force(c0,cb(:,:,i),csm(:,:,i),sc0,e2,&
               ddxc,psi,drhoe,eigv,&
               nstate,hamb,orbital)
          CALL dscal(2*ncpw%ngw*nstate,-1._real_8,csm(:,:,i),1)
          CALL lr_ortho(nstate,c0,csm(:,:,i))
       ENDDO
       ! ..Davidson matrices
       namat=ndim
       CALL zeroing(apmat)!,ndmax*nd)
       CALL zeroing(ammat)!,ndmax*nd)
       DO i=1,ndim
          DO j=i,ndim
             apmat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),csp(:,:,j))
             apmat(j,i)=apmat(i,j)
             ammat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),csm(:,:,j))
             ammat(j,i)=ammat(i,j)
          ENDDO
       ENDDO
       CALL mp_sum(apmat,nd*ndim,parai%allgrp)
       CALL mp_sum(ammat,nd*ndim,parai%allgrp)
       ! ..Calculate BMMAT = AMMAT * APMAT
       IF (paral%parent) THEN
          CALL dgemm('N','N',ndim,ndim,ndim,1._real_8,ammat,nd,apmat,nd,&
               0._real_8,bmmat,nd)
          ! ..diagonalise
          CALL dgeev('V','V',ndim,bmmat,nd,deig,aux,ammat,nd,apmat,&
               nd,work,lwork,info)
          IF (info.NE.0) CALL stopgm('TD_NHDAV','INFO DGEEV',& 
               __LINE__,__FILE__)
          tx=0._real_8
          DO i=1,ndim
             tx=tx+ABS(aux(i))
          ENDDO
          IF (tx.GT.1.e-10_real_8) CALL stopgm('TD_NHDAV',&
               'EIGNEVALUES NOT REAL',& 
               __LINE__,__FILE__)
          CALL sort2(deig,ndim,ind)
          CALL dcopy(nd*ndim,ammat,1,bmmat,1)
          DO i=1,ndim
             CALL dcopy(ndim,bmmat(1,ind(i)),1,ammat(1,i),1)
          ENDDO
          CALL dcopy(nd*ndim,apmat,1,bmmat,1)
          DO i=1,ndim
             CALL dcopy(ndim,bmmat(1,ind(i)),1,apmat(1,i),1)
          ENDDO
          DO i=1,ndim
             deig(i)=SQRT(deig(i))
          ENDDO
       ENDIF
       CALL mp_bcast(apmat,SIZE(apmat),parai%source,parai%allgrp)
       CALL mp_bcast(ammat,SIZE(ammat),parai%source,parai%allgrp)
       CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
       IF (ndim+1.GT.ndmax) THEN
          ! ..reduce trial vectors
          CALL zeroing(c1)!,ngw*nstate*nroot)
          DO i=1,nroot
             DO j=1,ndim
                CALL daxpy(2*ncpw%ngw*nstate,apmat(j,i),cb(:,:,j),1,&
                     c1(:,:,i,1),1)
                CALL daxpy(2*ncpw%ngw*nstate,ammat(j,i),cb(:,:,j),1,&
                     c1(:,:,i,2),1)
             ENDDO
          ENDDO
          ntogo=ntogo-ndim/2
          GOTO 1000
       ENDIF
       ! ..Residuals
       nsigma=ndim
       somax=0._real_8
       ndold=ndim
       DO i = nconv+1, nroot
          CALL zeroing(cr)!,2*ngw*nstate)
          DO j=1,ndim
             CALL daxpy(2*ncpw%ngw*nstate,ammat(j,i),csm(:,:,j),1,cr(:,:,1),1)
             fee=-deig(i)*apmat(j,i)
             CALL daxpy(2*ncpw%ngw*nstate,fee,cb(:,:,j),1,cr(:,:,1),1)
             CALL daxpy(2*ncpw%ngw*nstate,apmat(j,i),csp(:,:,j),1,cr(:,:,2),1)
             fee=-deig(i)*ammat(j,i)
             CALL daxpy(2*ncpw%ngw*nstate,fee,cb(:,:,j),1,cr(:,:,2),1)
          ENDDO
          so=dotps(ncpw%ngw,nstate,cr(:,:,1),cr(:,:,1)) +&
               dotps(ncpw%ngw,nstate,cr(:,:,2),cr(:,:,2))
          CALL mp_sum(so,parai%allgrp)
          ! ..Check for convergence
          IF (so.LT.td02%epstdav) nconv=nconv+1
          IF (so.GT.somax) somax=so
          IF (td03%tda) THEN
             IF (so.GT.td02%epstdav.AND.ndim.LT.ndmax) THEN
                ndim=ndim+1
                DO is=1,nstate
                   DO ig=1,ncpw%ngw
                      tx=1._real_8/(vpp(ig)-eigv(is)-deig(i))
                      cb(ig,is,ndim)=tx*cr(ig,is,1)
                   ENDDO
                ENDDO
             ENDIF
          ELSE
             IF (so.GT.td02%epstdav.AND.ndim.LT.ndmax-1) THEN
                ndim=ndim+2
                DO is=1,nstate
                   DO ig=1,ncpw%ngw
                      tx=1._real_8/(vpp(ig)-eigv(is)-deig(i))
                      cb(ig,is,ndim-1)=tx*cr(ig,is,1)
                      cb(ig,is,ndim)=tx*cr(ig,is,2)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDDO
       t2=m_walltime()
       tcpu=(t2-t1)*0.001_real_8
       IF ((paral%parent .AND. kprint.GE.0) .AND.paral%io_parent)&
            WRITE(6,'(I6,7X,I5,9X,I5,T34,G20.8,T58,F8.2)')&
            istep,nconv,namat,somax,tcpu
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (nconv.EQ.nroot .OR. istep.EQ.td01%ndavmax .OR. soft_com%exsoft) GOTO 999
       ! ..orthogonalize
       DO i=ndold+1,ndim
          CALL lr_ortho(nstate,c0,cb(:,:,i))
          DO j=1,i-1
             so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,j))
             CALL mp_sum(so,parai%allgrp)
             CALL daxpy(2*ncpw%ngw*nstate,-so,cb(:,:,j),1,cb(:,:,i),1)
          ENDDO
          so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,i))
          CALL mp_sum(so,parai%allgrp)
          CALL dscal(2*ncpw%ngw*nstate,1._real_8/SQRT(so),cb(:,:,i),1)
       ENDDO
    ENDDO
999 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ..Davidson matrices
    CALL zeroing(apmat)!,ndmax*nd)
    CALL zeroing(ammat)!,ndmax*nd)
    ndim=MIN(ndim,nsigma)
    DO i=1,ndim
       DO j=i,ndim
          apmat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),csp(:,:,j))
          apmat(j,i)=apmat(i,j)
          ammat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),csm(:,:,j))
          ammat(j,i)=ammat(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(apmat,nd*ndim,parai%allgrp)
    CALL mp_sum(ammat,nd*ndim,parai%allgrp)
    ! ..Calculate BMMAT = AMMAT * APMAT
    IF (paral%parent) THEN
       CALL dgemm('N','N',ndim,ndim,ndim,1._real_8,ammat,nd,apmat,nd,&
            0._real_8,bmmat,nd)
       ! ..diagonalise
       CALL dgeev('V','V',ndim,bmmat,nd,deig,aux,ammat,nd,apmat,&
            nd,work,lwork,info)
       IF (info.NE.0) CALL stopgm('TD_NHDAV','INFO DGEEV',& 
            __LINE__,__FILE__)
       tx=0._real_8
       DO i=1,ndim
          tx=tx+ABS(aux(i))
       ENDDO
       IF (tx.GT.1.e-10_real_8) CALL stopgm('TD_NHDAV',&
            'EIGNEVALUES NOT REAL',& 
            __LINE__,__FILE__)
       CALL sort2(deig,ndim,ind)
       CALL dcopy(nd*ndim,ammat,1,bmmat,1)
       DO i=1,ndim
          CALL dcopy(ndim,bmmat(1,ind(i)),1,ammat(1,i),1)
       ENDDO
       CALL dcopy(nd*ndim,apmat,1,bmmat,1)
       DO i=1,ndim
          CALL dcopy(ndim,bmmat(1,ind(i)),1,apmat(1,i),1)
       ENDDO
       DO i=1,ndim
          deig(i)=SQRT(deig(i))
       ENDDO
    ENDIF
    CALL mp_bcast(apmat,SIZE(apmat),parai%source,parai%allgrp)
    CALL mp_bcast(ammat,SIZE(ammat),parai%source,parai%allgrp)
    CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
    CALL zeroing(c1)!,ngw*nstate*nroot)
    DO i=1,nroot
       DO j=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,apmat(j,i),cb(:,:,j),1,&
               c1(:,:,i,1),1)
          CALL daxpy(2*ncpw%ngw*nstate,ammat(j,i),cb(:,:,j),1,&
               c1(:,:,i,2),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL dcopy (nroot,deig,1,teig,1)
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.nconv.LT.nroot) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,T64,"!!")')&
            ' DAVIDSON| NOT ALL ROOTS ARE CONVERGED'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(csp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(csm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(apmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ammat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bmmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(deig,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_nhdav
  ! ==================================================================
  SUBROUTINE give_scr_td_nhdav(ltd_nhdav,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltd_nhdav
    CHARACTER(len=30)                        :: tag

    INTEGER, PARAMETER                       :: mroot = 20, nd = 5*mroot 

    INTEGER                                  :: llr_force, llr_ortho, local

! ==--------------------------------------------------------------==

    local = nd*nd + nd
    CALL give_scr_lr_force(llr_force,"TDA",tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    ltd_nhdav=MAX(local,llr_force,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_td_nhdav
  ! ==================================================================

END MODULE td_nhdav_utils
