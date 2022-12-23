MODULE td_dav_utils
  USE cnst,                            ONLY: ry
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr02,&
                                             td01,&
                                             td02,&
                                             td03,&
                                             tshl
  USE lr_force_utils,                  ONLY: give_scr_lr_force,&
                                             lr_force,&
                                             lr_force_sub
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_ortho,&
                                             lr_orthos
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pw_hfx_resp_types,               ONLY: hfx_resp_env
  USE randtowf_utils,                  ONLY: randtowf
  USE soft,                            ONLY: soft_com
  USE sort_utils,                      ONLY: sort2
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_dav
  PUBLIC :: give_scr_td_dav
  PUBLIC :: td_sub
  PUBLIC :: setvpp
  PUBLIC :: dotps

CONTAINS

  ! ==================================================================
  SUBROUTINE td_dav(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,teig,&
       nstate,nroot,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==            Davidson diagonalization of cntl%tddft matrix          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:), drhoe(:,:), teig(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*)
    INTEGER                                  :: nroot
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_dav'
    INTEGER, PARAMETER                       :: mroot = 20, nd = 5*mroot

    COMPLEX(real_8), ALLOCATABLE             :: cb(:,:,:), cs(:,:,:)
    INTEGER :: i, ib, ierr, ig, indexr(nd), info, is, istep, isub, j, k, &
      lwork, namat, nconv, ndim, ndmax, ndold, nsigma, nzero
    LOGICAL                                  :: trdiis, tzero
    REAL(real_8)                             :: deij, e2(5), fee, sig, so, &
                                                soi(mroot), somax, t1, t2, &
                                                tcpu, vppmin
    REAL(real_8), ALLOCATABLE                :: amat(:,:), deig(:), &
                                                soback(:), vpp(:), work(:)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (.NOT.td03%tda) THEN
       CALL stopgm("TD_DAV","ONLY TAMM-DANCOFF APPROXIMATION ALLOWED",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (nroot.GT.mroot) THEN
       CALL stopgm("TD_DAV","MAX NUMBER OF ROOTS EXCEEDED",& 
            __LINE__,__FILE__)
    ENDIF
    ndmax = MAX(nroot+1,td01%ndavspace)
    IF (ndmax.GT.nd) THEN
       CALL stopgm("TD_DAV","MAX SUB SPACE SIZE EXCEEDED",& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(cs(ncpw%ngw,nstate,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cb(ncpw%ngw,nstate,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(amat(nd,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(deig(ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(amat)!,nd*ndmax)
    ! ==--------------------------------------------------------------==
    ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(soback(td01%ndavmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF ( paral%parent ) THEN
       IF ( kprint.GE.0 ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T13,A,T64,"==")')&
               '   DAVIDSON DIAGONALISATION OF TDDFT MATRIX '
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
       ENDIF
       IF (paral%io_parent.AND.kprint.EQ.0) THEN
          WRITE(6,'(A,T43,A)')&
               '  ITER      STATES      SUBSPACE','   RESIDUAL        TCPU'
       ENDIF
    ENDIF
    ! 
    trdiis=.FALSE.
    cntl%prec=.TRUE.
    CALL setvpp(vpp,ncpw%ngw,lr02%lr_hthrs,vppmin)
    ! ..initilize vectors
    CALL dcopy(2*ncpw%ngw*nstate*nroot,c1,1,cb,1)
    DO i = 1, nroot
       CALL lr_ortho(nstate,c0,cb(:,:,i))
       DO j=1,i-1
          so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,j))
          CALL mp_sum(so,parai%allgrp)
          CALL daxpy(2*ncpw%ngw*nstate,-so,cb(:,:,j),1,cb(:,:,i),1)
       ENDDO
       so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,i))
       CALL mp_sum(so,parai%allgrp)
       so=1._real_8/SQRT(so)
       CALL dscal(2*ncpw%ngw*nstate,so,cb(:,:,i),1)
    ENDDO
    nzero=0
    ndim=nroot
    nsigma=0
    nconv=0
    IF (paral%parent) THEN
       lwork=3*ndmax! TODO optimal LWORK
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF

    IF (hfx_resp_env%switch_lin2) THEN
       IF (paral%io_parent) WRITE(*,*) 'switch on lin lin'
       hfx_resp_env%use_lin_lin=.TRUE.
       hfx_resp_env%switch_lin2=.FALSE.
    ENDIF

    DO istep=1,td01%ndavmax
       t1=m_walltime()
1100   CONTINUE
       ! ..initialize sigma vectors
       DO i = nsigma+1,ndim
          CALL lr_force(c0,cb(:,:,i),cs(:,:,i),sc0,e2,ddxc,psi,drhoe,&
               eigv,nstate,"TDA",orbital)
          CALL dscal(2*ncpw%ngw*nstate,-1._real_8,cs(:,:,i),1)
          CALL lr_ortho(nstate,c0,cs(:,:,i))
       ENDDO
       ! ..Davidson matrix
       namat=ndim
       CALL zeroing(amat)!,ndmax*nd)
       DO i=1,ndim
          DO j=i,ndim
             amat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),cs(:,:,j))
             amat(j,i)=amat(i,j)
          ENDDO
       ENDDO
       CALL mp_sum(amat,nd*ndim,parai%allgrp)
       IF (paral%parent) THEN
          CALL dsyev('V','L',ndim,amat,nd,deig,work,lwork,info)
          IF (info.NE.0) CALL stopgm('TD_DAV','INFO! = 0 AFTER DSYEV',& 
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(amat,SIZE(amat),parai%source,parai%allgrp)
       CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
       ! ..Check for zero eigenvalues
       tzero=.FALSE.
       ! ..TSH[
       IF (.NOT.tshl%tdtully) THEN
          DO i=1,nroot
             IF (ABS(deig(i))*2._real_8*ry .LT. 0.01_real_8) THEN
                nzero=nzero+1
                tzero=.TRUE.
             ENDIF
          ENDDO
       ENDIF
       ! ..TSH]
       IF (tzero) THEN
          ! ..Remove zero eigenvectors
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' WARNING| Remove zero eigenvalues, reinitialize states'
          CALL davred(c0,c1,cb,cs,amat,ndim,nroot,nstate,nd)
          DO i=1,nroot
             IF (ABS(deig(i))*2._real_8*ry .LT. 0.01_real_8) THEN
                CALL randtowf(c1(:,:,i),nstate,0,0)
                CALL lr_ortho(nstate,c0,c1(:,:,i))
                deig(i)=i*1.e10_real_8
             ENDIF
          ENDDO
          ! sort states, new states last
          CALL sort2(deig,nroot,indexr)
          DO i = 1, nroot
             j=indexr(i)
             CALL dcopy(2*ncpw%ngw*nstate,c1(:,:,i),1,cb(:,:,j),1)
          ENDDO
          ! orthogonalize
          DO i = 1, nroot
             DO j=1,i-1
                so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,j))
                CALL mp_sum(so,parai%allgrp)
                CALL daxpy(2*ncpw%ngw*nstate,-so,cb(:,:,j),1,cb(:,:,i),1)
             ENDDO
             so=dotps(ncpw%ngw,nstate,cb(:,:,i),cb(:,:,i))
             CALL mp_sum(so,parai%allgrp)
             so=1._real_8/SQRT(so)
             CALL dscal(2*ncpw%ngw*nstate,so,cb(:,:,i),1)
          ENDDO
          ndim=nroot
          nsigma=0
          nconv=0
          GOTO 1100
       ELSEIF (ndim+1.GT.ndmax) THEN
          ! ..Reduce vector subspace
          CALL davred(c0,c1,cb,cs,amat,ndim,nroot,nstate,nd)
       ENDIF
       ! ..Residuals
       nsigma=ndim
       somax=0._real_8
       ndold=ndim
       DO i = nconv+1, nroot
          CALL zeroing(c2(:,1:nstate))!,ngw*nstate)
          DO j=1,ndim
             CALL daxpy(2*ncpw%ngw*nstate,amat(j,i),cs(:,:,j),1,c2,1)
             fee=-deig(i)*amat(j,i)
             CALL daxpy(2*ncpw%ngw*nstate,fee,cb(:,:,j),1,c2,1)
          ENDDO
          so=dotps(ncpw%ngw,nstate,c2,c2)
          CALL mp_sum(so,parai%allgrp)
          soi(i)=so
          ! ..Check for convergence
          IF (so.LT.td02%epstdav.AND.i.EQ.nconv+1) nconv=nconv+1
          IF (so.GT.somax) somax=so
          IF (so.GT.td02%epstdav.AND.ndim.LT.ndmax) THEN
             ndim=ndim+1
             DO is=1,nstate
                IF (INDEX(orbital,"CANON").NE.0) THEN
                   ib=is
                   sig=-1._real_8
                ELSE
                   ib=(is-1)*nstate+is
                   sig=1._real_8
                ENDIF
                deij=sig*eigv(ib)-deig(i)
                DO ig=1,ncpw%ngw
                   cb(ig,is,ndim)=c2(ig,is)/(vpp(ig)+deij)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       t2=m_walltime()
       tcpu=(t2-t1)*0.001_real_8
       IF (paral%parent) THEN
          IF (kprint.EQ.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(I6,7X,I5,9X,I5,T34,G20.8,T58,F8.2)')&
                  istep,nconv,namat,somax,tcpu
          ELSEIF (kprint.EQ.1) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,64("-"))')
             IF (paral%io_parent)&
                  WRITE(6,'(A,T20,I6,T34,A,T60,I6)') ' STEP  ',istep,&
                  'NUMBER OF CONVERGED STATES   ',nconv
             IF (paral%io_parent)&
                  WRITE(6,'(A,T20,I6,T34,A,T54,G12.4)')&
                  ' DAVIDSON MATRIX ',namat,'MAXIMUM RESIDIUM  ',somax
             IF (paral%io_parent)&
                  WRITE(6,'(T34,A,T58,F8.1)') 'WALL TIME  ',tcpu
             DO i = 1, nconv
                IF (paral%io_parent)&
                     WRITE(6,'(A,T20,I6,T34,F14.6," eV",T54,G12.4)')&
                     ' STATE ',i,deig(i)*ry*2._real_8,soi(i)
             ENDDO
             DO i = nconv+1, nroot
                IF (paral%io_parent)&
                     WRITE(6,'(A,T20,I6," NC",T34,F14.6," eV",T54,G12.4)')&
                     ' STATE ',i,deig(i)*ry*2._real_8,soi(i)
             ENDDO
          ENDIF
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (nconv.EQ.nroot .OR. istep.EQ.td01%ndavmax .OR. soft_com%exsoft) GOTO 999
       soback(istep)=somax
       IF (istep.GT.4) THEN
          IF (somax.LT.td02%rdiistin) THEN
             IF (somax.GT.0.9_real_8*soback(istep-4)) GOTO 999
          ENDIF
       ENDIF
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
    IF (nconv.LT.nroot) trdiis=.TRUE.
    ! ==--------------------------------------------------------------==
    ! ..Davidson matrix
    CALL zeroing(amat)!,ndmax*nd)
    ndim=MIN(ndim,nsigma)
    DO i=1,ndim
       DO j=i,ndim
          amat(i,j)=dotps(ncpw%ngw,nstate,cb(:,:,i),cs(:,:,j))
          amat(j,i)=amat(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(amat,nd*ndim,parai%allgrp)
    IF (paral%parent) THEN
       lwork=3*ndim-1
       CALL dsyev('V','L',ndim,amat,nd,deig,work,lwork,info)
       IF (info.NE.0) CALL stopgm('TD_DAV','INFO ! = 0 AFTER DSYEV',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) THEN
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(amat,SIZE(amat),parai%source,parai%allgrp)
    CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
    CALL zeroing(c1(:,:,1:nroot))!,ngw*nstate*nroot)
    DO i=1,nroot
       DO j=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,amat(j,i),cb(:,:,j),1,c1(:,:,i),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (soft_com%exsoft) trdiis=.FALSE.
    IF (trdiis) THEN
       j=nconv+1
       DO i=j,nroot
          DO k=1,i-1
             so=dotps(ncpw%ngw,nstate,c1(:,:,i),c1(:,:,k))
             CALL mp_sum(so,parai%allgrp)
             CALL daxpy(2*ncpw%ngw*nstate,-so,c1(:,:,k),1,c1(:,:,i),1)
          ENDDO
          so=dotps(ncpw%ngw,nstate,c1(:,:,i),c1(:,:,i))
          CALL mp_sum(so,parai%allgrp)
          so=1._real_8/SQRT(so)
          CALL dscal(2*ncpw%ngw*nstate,so,c1(:,:,i),1)
          CALL td_rdiis(c0,c1,c2,sc0,cb,cs,i,ndmax,nconv,vpp,ddxc,psi,&
               eigv,drhoe,deig(i),nstate,orbital,kprint)
       ENDDO
    ENDIF
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
    DEALLOCATE(soback,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(amat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(deig,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_dav
  ! ==================================================================
  SUBROUTINE give_scr_td_dav(ltd_dav,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltd_dav
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: llr_force, llr_ortho

    CALL give_scr_lr_force(llr_force,"TDA",tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    ltd_dav=MAX(llr_force,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_td_dav
  ! ==================================================================
  FUNCTION dotps(ngw,nstate,ca,cb)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngw, nstate
    COMPLEX(real_8)                          :: ca(:,:), cb(:,:)
    REAL(real_8)                             :: dotps

    INTEGER                                  :: i

    dotps=0._real_8
    DO i=1,nstate
       dotps=dotps+dotp(ngw,ca(:,i),cb(:,i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dotps
  ! ==================================================================
  SUBROUTINE td_sub(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,teig,&
       msub,nstate,nroot,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==               updates the wavefunctions                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), &
                                                c2(ncpw%ngw,*), sc0(*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:), drhoe(:,:), teig(:)
    INTEGER                                  :: msub
    COMPLEX(real_8)                          :: c1(ncpw%ngw,msub,*)
    INTEGER                                  :: nstate, nroot
    CHARACTER(len=30)                        :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_sub'
    INTEGER, PARAMETER                       :: mroot = 20, nd = 5*mroot

    COMPLEX(real_8), ALLOCATABLE             :: cb(:,:,:), cs(:,:,:)
    INTEGER                                  :: i, ib, id, ierr, ig, info, &
                                                is, istep, isub, j, lwork, &
                                                namat, nconv, ndim, ndmax, &
                                                ndold, nsigma
    REAL(real_8)                             :: e2(5), fee, sig, so, somax, &
                                                t1, t2, tcpu, vp
    REAL(real_8), ALLOCATABLE                :: amat(:,:), deig(:), vpp(:), &
                                                work(:)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (.NOT.td03%tda) THEN
       CALL stopgm("TD_SUB","ONLY TAMM-DANCOFF APPROXIMATION ALLOWED",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (nroot.GT.mroot) THEN
       CALL stopgm("TD_SUB","MAX NUMBER OF ROOTS EXCEEDED",& 
            __LINE__,__FILE__)
    ENDIF
    ndmax = MAX(nroot+1,td01%ndavspace)
    IF (ndmax.GT.nd) THEN
       CALL stopgm("TD_SUB","MAX SUB SPACE SIZE EXCEEDED",& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cs(ncpw%ngw,msub,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cb(ncpw%ngw,msub,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(amat(nd,ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(deig(ndmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(amat)!,nd*ndmax)
    ! ==--------------------------------------------------------------==
    cntl%prec=.TRUE.
    CALL ksdiag(vpp)
    !$omp parallel do private(IG,VP)
#ifdef __SR11000
    !poption parallel, tlocal(IG,VP)
#endif
    DO ig=1,ncpw%ngw
       vp=vpp(ig)
       vpp(ig)=MAX(vp,lr02%lr_hthrs)
    ENDDO
    ! ==--------------------------------------------------------------==
    IF ( paral%io_parent .AND. kprint.GE.0 ) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T13,A,T64,"==")')'   DAVIDSON DIAGONALISATION OF TDDFT MATRIX '
       WRITE(6,'(1X,"==",T6,A,T64,"==")')' TAMM-DANCOFF SUBSPACE APPROXIMATION (CANONICAL ORBITALS) '
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(A,T43,A)') '  ITER      STATES      SUBSPACE',&
            'RESIDUAL        TCPU'
    ENDIF
    ! ..initilize vectors
    CALL dcopy(2*ncpw%ngw*msub*nroot,c1,1,cb,1)
    DO i = 1, nroot
       CALL lr_orthos(nstate,msub,c0,cb(:,:,i))
       so=dotps(ncpw%ngw,msub,cb(:,:,i),cb(:,:,i))
       CALL mp_sum(so,parai%allgrp)
       CALL dscal(2*ncpw%ngw*msub,1._real_8/SQRT(so),cb(:,:,i),1)
    ENDDO
    ndim=nroot
    nsigma=0
    nconv=0
    IF (paral%parent) THEN
       lwork=3*ndmax
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    DO istep=1,td01%ndavmax
       t1=m_walltime()
       ! ..initialize sigma vectors
       IF (INDEX(orbital,"CANON").NE.0) THEN
          id=nstate-msub+1
          ib=nstate-msub+1
       ELSE
          id=nstate-msub+1
          ib=1
       ENDIF
       DO i = nsigma+1,ndim
          CALL lr_force_sub(c0(1,id),cb(:,:,i),cs(:,:,i),sc0,e2,ddxc,&
               psi,drhoe,eigv(ib:),nstate,msub,"TDA",orbital)
          CALL dscal(2*ncpw%ngw*msub,-1._real_8,cs(:,:,i),1)
          CALL lr_orthos(nstate,msub,c0,cs(:,:,i))
       ENDDO
       ! ..Davidson matrix
       namat=ndim
       CALL zeroing(amat)!,ndmax*nd)
       DO i=1,ndim
          DO j=i,ndim
             amat(i,j)=dotps(ncpw%ngw,msub,cb(:,:,i),cs(:,:,j))
             amat(j,i)=amat(i,j)
          ENDDO
       ENDDO
       CALL mp_sum(amat,nd*ndim,parai%allgrp)
       IF (paral%parent) THEN
          CALL dsyev('V','L',ndim,amat,nd,deig,work,lwork,info)
          IF (info.NE.0) CALL stopgm('TD_DAV','INFO! = 0 AFTER DSYEV',& 
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(amat,SIZE(amat),parai%source,parai%allgrp)
       CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
       IF (ndim+1.GT.ndmax) THEN
          ! ..reduce trial vectors
          CALL zeroing(c1(:,:,1:nroot))!,ngw*msub*nroot)
          DO i=1,nroot
             DO j=1,ndim
                CALL daxpy(2*ncpw%ngw*msub,amat(j,i),cb(:,:,j),1,c1(:,:,i),1)
             ENDDO
          ENDDO
          CALL zeroing(cb(:,:,1:nroot))!,ngw*msub*nroot)
          DO i=1,nroot
             DO j=1,ndim
                CALL daxpy(2*ncpw%ngw*msub,amat(j,i),cs(:,:,j),1,cb(:,:,i),1)
             ENDDO
          ENDDO
          CALL dcopy(2*msub*ncpw%ngw*nroot,cb,1,cs,1)
          CALL dcopy(2*msub*ncpw%ngw*nroot,c1,1,cb,1)
          ndim=nroot
          CALL zeroing(amat)!,ndmax*nd)
          DO i=1,nroot
             amat(i,i)=1._real_8
          ENDDO
       ENDIF
       ! ..Residuals
       nsigma=ndim
       somax=0._real_8
       ndold=ndim
       DO i = nconv+1, nroot
          CALL zeroing(c2(:,1:msub))!,ngw*msub)
          DO j=1,ndim
             CALL daxpy(2*ncpw%ngw*msub,amat(j,i),cs(:,:,j),1,c2,1)
             fee=-deig(i)*amat(j,i)
             CALL daxpy(2*ncpw%ngw*msub,fee,cb(:,:,j),1,c2,1)
          ENDDO
          so=dotps(ncpw%ngw,msub,c2(:,1:msub),c2(:,1:msub))
          CALL mp_sum(so,parai%allgrp)
          ! ..Check for convergence
          IF (so.LT.td02%epstdav) nconv=nconv+1
          IF (so.GT.somax) somax=so
          IF (so.GT.td02%epstdav.AND.ndim.LT.ndmax) THEN
             ndim=ndim+1
             DO is=1,msub
                IF (INDEX(orbital,"CANON").NE.0) THEN
                   id=nstate-msub+is
                   sig=-1._real_8
                ELSE
                   id=(is-1)*msub+is
                   sig=1._real_8
                ENDIF
                DO ig=1,ncpw%ngw
                   cb(ig,is,ndim)=c2(ig,is)/(vpp(ig)+sig*eigv(id)-deig(i))
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       t2=m_walltime()
       tcpu=(t2-t1)*0.0001_real_8
       IF (paral%io_parent.AND.kprint.GE.0) THEN
          WRITE(6,'(I6,7X,I5,9X,I5,T34,1PG20.8,T58,F8.2)')&
               istep,nconv,namat,somax,tcpU
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (nconv.EQ.nroot .OR. istep.EQ.td01%ndavmax .OR. soft_com%exsoft) GOTO 999
       ! ..orthogonalize
       DO i=ndold+1,ndim
          CALL lr_orthos(nstate,msub,c0,cb(:,:,i))
          DO j=1,i-1
             so=dotps(ncpw%ngw,msub,cb(:,:,i),cb(:,:,j))
             CALL mp_sum(so,parai%allgrp)
             CALL daxpy(2*ncpw%ngw*msub,-so,cb(:,:,j),1,cb(:,:,i),1)
          ENDDO
          so=dotps(ncpw%ngw,msub,cb(:,:,i),cb(:,:,i))
          CALL mp_sum(so,parai%allgrp)
          CALL dscal(2*ncpw%ngw*msub,1._real_8/SQRT(so),cb(:,:,i),1)
       ENDDO
    ENDDO
999 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ..Davidson matrix
    CALL zeroing(amat)!,ndmax*nd)
    ndim=MIN(ndim,nsigma)
    DO i=1,ndim
       DO j=i,ndim
          amat(i,j)=dotps(ncpw%ngw,msub,cb(:,:,i),cs(:,:,j))
          amat(j,i)=amat(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(amat,nd*ndim,parai%allgrp)
    IF (paral%parent) THEN
       CALL dsyev('V','L',ndim,amat,nd,deig,work,lwork,info)
       IF (info.NE.0) CALL stopgm('TD_DAV','INFO ! = 0 AFTER DSYEV',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL mp_bcast(amat,SIZE(amat),parai%source,parai%allgrp)
    CALL mp_bcast(deig,SIZE(deig),parai%source,parai%allgrp)
    CALL zeroing(c1(:,:,1:nroot))!,ngw*msub*nroot)
    DO i=1,nroot
       DO j=1,ndim
          CALL daxpy(2*ncpw%ngw*msub,amat(j,i),cb(:,:,j),1,c1(:,:,i),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL dcopy (nroot,deig,1,teig,1)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND.nconv.LT.nroot) THEN
       WRITE(6,'(1X,64("!"))')
       WRITE(6,'(" !!",A,T64,"!!")')&
            ' DAVIDSON| NOT ALL ROOTS ARE CONVERGED'
       WRITE(6,'(1X,64("!"))')
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(amat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(deig,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_sub
  ! ==================================================================
  SUBROUTINE setvpp(vpp,ngw,thres,vppmin)
    REAL(real_8)                             :: vpp(:)
    INTEGER                                  :: ngw
    REAL(real_8)                             :: thres, vppmin

    INTEGER                                  :: ig
    REAL(real_8)                             :: vp

! ==--------------------------------------------------------------==

    CALL ksdiag(vpp)
    vppmin=1.e9_real_8
    DO ig=1,ngw
       vp=vpp(ig)
       vpp(ig)=MAX(vp,thres)
       IF (vpp(ig).LT.vppmin) vppmin=vpp(ig)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setvpp
  ! ==================================================================
  SUBROUTINE davred(c0,c1,cb,cs,amat,ndim,nroot,nstate,nd)
    ! ==--------------------------------------------------------------==
    ! ==            Reduce Davidson vectors to nroots                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), cb(:,:,:), cs(:,:,:)
    REAL(real_8)                             :: amat(:,:)
    INTEGER                                  :: ndim, nroot, nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*)
    INTEGER                                  :: nd

    INTEGER                                  :: i, j

! 
! ==--------------------------------------------------------------==
! ..reduce trial vectors

    CALL zeroing(c1(:,:,1:nroot))!,ngw*nstate*nroot)
    DO i=1,nroot
       DO j=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,amat(j,i),cb(:,:,j),1,c1(:,:,i),1)
       ENDDO
    ENDDO
    CALL zeroing(cb(:,:,1:nroot))!,ngw*nstate*nroot)
    DO i=1,nroot
       DO j=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,amat(j,i),cs(:,:,j),1,cb(:,:,i),1)
       ENDDO
    ENDDO
    CALL dcopy(2*nstate*ncpw%ngw*nroot,cb,1,cs,1)
    CALL dcopy(2*nstate*ncpw%ngw*nroot,c1,1,cb,1)
    ndim=nroot
    CALL zeroing(amat)!,nroot*nd)
    DO i=1,nroot
       amat(i,i)=1._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE davred
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE td_rdiis(c0,c1,c2,sc0,cb,cs,nroot,nbdim,nconv,&
       vpp,ddxc,psi,eigv,drhoe,&
       teig,nstate,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==               updates the wavefunctions                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(ncpw%ngw,*), &
                                                sc0(*), cb(:,:,:), cs(:,:,:)
    INTEGER                                  :: nroot, nbdim, nconv
    REAL(real_8)                             :: vpp(:), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:), drhoe(:,:), teig
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*)
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_rdiis'
    INTEGER, PARAMETER                       :: matmax = 5 

    INTEGER                                  :: ib, ierr, ig, is, istep, &
                                                isub, ndiis, ndim, nnew, &
                                                nnow, nold, nrestdiis
    REAL(real_8) :: bc_array(matmax,matmax), deij, dismat(matmax,matmax), &
      e2(5), res1, res2, sig, so, t1, t2, tcpu, vc(matmax), vppmin

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (.NOT.td03%tda) CALL stopgm("TD_DAV",&
         "ONLY TAMM-DANCOFF APPROXIMATION ALLOWED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF ( paral%parent .AND. kprint.GE.0 ) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T62,I4)')&
            ' >>> RESIDUAL DIIS OPTIMIZATION OF ROOT ',nroot
       IF (paral%io_parent)&
            WRITE(6,'(A,T43,A)')&
            '  ITER      ENERGY      SUBSPACE','   RESIDUAL        TCPU'
    ENDIF
    ! 
    res1=1.e10_real_8
    res2=1.e10_real_8
    ndiis=MIN(nbdim,matmax-1)
    CALL zeroing(dismat)!,matmax*matmax)
    t1=m_walltime()
    cntl%prec=.TRUE.
    CALL setvpp(vpp,ncpw%ngw,lr02%lr_hthrs,vppmin)
    ! ..initilize vectors
    CALL dcopy(2*ncpw%ngw*nstate,c1(:,:,nroot),1,cb,1)
    CALL lr_ortho(nstate,c0,cb)
    so=dotps(ncpw%ngw,nstate,cb(:,:,1),cb(:,:,1))
    CALL mp_sum(so,parai%allgrp)
    so=1._real_8/SQRT(so)
    CALL dscal(2*ncpw%ngw*nstate,so,cb,1)
    nrestdiis=0
    nnew=1
    DO istep=1,td01%ntdiismax
       nnow=nnew
       ndim=MIN(ndiis,istep)
       ! ..calculate error vector (residual)
       CALL lr_force(c0,cb(:,:,nnow),cs(:,:,nnow),sc0,e2,&
            ddxc,psi,drhoe,eigv,&
            nstate,"TDA",orbital)
       CALL lr_ortho(nstate,c0,cs(:,:,nnow))
       teig=-dotps(ncpw%ngw,nstate,cb(:,:,nnow),cs(:,:,nnow))
       CALL mp_sum(teig,parai%allgrp)
       CALL daxpy(2*ncpw%ngw*nstate,teig,cb(:,:,nnow),1,cs(:,:,nnow),1)
       so=dotps(ncpw%ngw,nstate,cs(:,:,nnow),cs(:,:,nnow))
       CALL mp_sum(so,parai%allgrp)
       t2=m_walltime()
       tcpu=(t2-t1)*0.001_real_8
       IF ((paral%parent.AND.kprint.GE.0) .AND.paral%io_parent)&
            WRITE(6,'(I6,3X,F9.5,9X,I5,T34,G20.8,T58,F8.2)')&
            istep,teig*2._real_8*ry,ndim,so,tcpu
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (so.LT.td02%epstdav) nconv=nconv+1
       IF (so.LT.td02%epstdav .OR. soft_com%exsoft) GOTO 1000
       IF (so.GT.td02%rdiistin) GOTO 1000
       IF (nrestdiis.GT.td01%nrestdmax) GOTO 1000
       IF (res2.LT.so) THEN
          nold=MOD(ndiis+nnow-2,ndiis)+1
          CALL dcopy(2*ncpw%ngw*nstate,cb(:,:,nold),1,c1(:,:,nroot),1)
          DO is=1,nroot-1
             so=dotps(ncpw%ngw,nstate,c1(:,:,is),c1(:,:,nroot))
             CALL mp_sum(so,parai%allgrp)
             CALL daxpy(2*ncpw%ngw*nstate,-so,c1(:,:,is),1,c1(:,:,nroot),1)
          ENDDO
          so=dotps(ncpw%ngw,nstate,c1(:,:,nroot),c1(:,:,nroot))
          CALL mp_sum(so,parai%allgrp)
          so=1._real_8/SQRT(so)
          CALL dscal(2*ncpw%ngw*nstate,so,c1(:,:,nroot),1)
          res1=1.e10_real_8
          res2=1.e10_real_8
          nnew=1
          CALL dcopy(2*ncpw%ngw*nstate,c1(:,:,nroot),1,cb,1)
          nrestdiis=nrestdiis+1
          GOTO 999
       ENDIF
       res2=res1
       res1=so
       t1=m_walltime()
       ! Apply preconditioner
       DO is=1,nstate
          IF (INDEX(orbital,"CANON").NE.0) THEN
             ib=is
             sig=-1._real_8
          ELSE
             ib=(is-1)*nstate+is
             sig=1._real_8
          ENDIF
          deij=sig*eigv(ib)+teig
          DO ig=1,ncpw%ngw
             cs(ig,is,nnow)=cs(ig,is,nnow)/(vpp(ig)+deij)
          ENDDO
       ENDDO
       ! cntl%diis matrix
       DO is=1,ndim
          so=dotps(ncpw%ngw,nstate,cs(:,:,is),cs(:,:,nnow))
          CALL mp_sum(so,parai%allgrp)
          dismat(is,nnow)=so
          dismat(nnow,is)=so
       ENDDO
       CALL dcopy(matmax*matmax,dismat,1,bc_array,1)
       DO is=1,ndim+1
          vc(is)=0._real_8
          bc_array(is,ndim+1)=-1._real_8
          bc_array(ndim+1,is)=-1._real_8
       ENDDO
       vc(ndim+1)=-1._real_8
       ! Solve System of Linear Equations
       IF (ndim.EQ.1) THEN
          ierr=0
          vc(1)=1._real_8
       ELSE
          ierr=0
          CALL dissol(bc_array,matmax,ndim+1,vc,ierr)
       ENDIF
       IF (ierr.NE.0) THEN
          DO is=1,ndim
             vc(is)=0._real_8
          ENDDO
          vc(nnow)=1._real_8
          ndiis=MAX(ndiis-1,1)
       ENDIF
       ! Setup new vector
       nnew=MOD(nnow,ndiis)+1
       CALL zeroing(c2(:,1:nstate))!,ngw*nstate)
       DO is=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,vc(is),cb(:,:,is),1,c2,1)
       ENDDO
       CALL dcopy(2*ncpw%ngw*nstate,c2,1,cb(:,:,nnew),1)
       CALL zeroing(c2(:,1:nstate))!,ngw*nstate)
       DO is=1,ndim
          CALL daxpy(2*ncpw%ngw*nstate,vc(is),cs(:,:,is),1,c2,1)
       ENDDO
       CALL lr_ortho(nstate,c0,c2)
       CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c2,1,cb(:,:,nnew),1)
       so=dotps(ncpw%ngw,nstate,cb(:,:,nnew),cb(:,:,nnew))
       CALL mp_sum(so,parai%allgrp)
       so=1._real_8/SQRT(so)
       CALL dscal(2*ncpw%ngw*nstate,so,cb(:,:,nnew),1)
999    CONTINUE
    ENDDO
    nnow=nnew
1000 CONTINUE
    CALL dcopy(2*ncpw%ngw*nstate,cb(:,:,nnow),1,c1(:,:,nroot),1)
    DO is=1,nroot-1
       so=dotps(ncpw%ngw,nstate,c1(:,:,is),c1(:,:,nroot))
       CALL mp_sum(so,parai%allgrp)
       CALL daxpy(2*ncpw%ngw*nstate,-so,c1(:,:,is),1,c1(:,:,nroot),1)
    ENDDO
    so=dotps(ncpw%ngw,nstate,c1(:,:,nroot),c1(:,:,nroot))
    CALL mp_sum(so,parai%allgrp)
    so=1._real_8/SQRT(so)
    CALL dscal(2*ncpw%ngw*nstate,so,c1(:,:,nroot),1)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_rdiis
  ! ==================================================================
  SUBROUTINE give_scr_td_rdiis(ltd_rdiis,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltd_rdiis
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: llr_force, llr_ortho

    CALL give_scr_lr_force(llr_force,"TDA",tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    ltd_rdiis=MAX(llr_force,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_td_rdiis
  ! ==================================================================
  SUBROUTINE dissol(b,ldb,ndim,v,ierr)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: b(:,:)
    INTEGER                                  :: ldb, ndim
    REAL(real_8)                             :: v(:)
    INTEGER                                  :: ierr

    INTEGER, PARAMETER                       :: maxdis = 20 

    INTEGER                                  :: lw, nr
    REAL(real_8)                             :: scr1(maxdis+1,maxdis+1), &
                                                scr2(maxdis+1,maxdis+1), &
                                                toleig

    IF (ndim.GT.maxdis+1) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'DISSOL! NDIM=',ndim,' MAXDIS+1=',maxdis+1
       CALL stopgm('DISSOL','NDIM GREATER THAN MAXDIS+1',& 
            __LINE__,__FILE__)
    ENDIF
    toleig=1.e-8_real_8
    lw=(maxdis+1)*(maxdis+1)
    CALL dgelss(ndim,ndim,1,b,ldb,v,ldb,scr1,toleig,nr,scr2,lw,ierr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dissol
  ! ==================================================================


END MODULE td_dav_utils
