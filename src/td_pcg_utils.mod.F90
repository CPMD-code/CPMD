MODULE td_pcg_utils
  USE cnst,                            ONLY: ry
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr02,&
                                             td01,&
                                             td02,&
                                             td03,&
                                             urot
  USE lr_force_utils,                  ONLY: give_scr_lr_force,&
                                             give_scr_lru_force,&
                                             lr_force_sub,&
                                             lru_force
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_orthos
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE orbrot_utils,                    ONLY: get_u,&
                                             xgrad
  USE parac,                           ONLY: parai,&
                                             paral
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE td_dav_utils,                    ONLY: dotps
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_pcg
  PUBLIC :: give_scr_td_pcg

CONTAINS

  ! ==================================================================
  SUBROUTINE td_pcg(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,teig,&
       msub,nstate,nroot,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==               updates the wavefunctions                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(*), drhoe(:,:), teig(:)
    INTEGER                                  :: msub
    COMPLEX(real_8)                          :: c1(ncpw%ngw,msub,*)
    INTEGER                                  :: nstate, nroot
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_pcg'

    COMPLEX(real_8), ALLOCATABLE             :: cb(:,:), ch(:,:), cs(:,:)
    INTEGER                                  :: i, id, ierr, ig, iroot, irot, &
                                                iroteff, is, isub, iter, j
    REAL(real_8)                             :: dgg, e2(3), gam, gg, lam, &
                                                ores, rdg, rg, rotb, so, t1, &
                                                t2, tcpu, teigold, vp
    REAL(real_8), ALLOCATABLE :: cgrad(:,:), fmat(:,:), peig(:), pmat(:,:), &
      ufmat(:,:), ugrad(:,:), uh(:,:), umat(:,:,:), vpp(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cs(ncpw%ngw,msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cb(ncpw%ngw,msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ch(ncpw%ngw,msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (.NOT.td03%tdacanon) THEN
       ALLOCATE(fmat(msub,msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ufmat(msub,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pmat(msub,msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(peig(msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(umat(nstate,msub,nroot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ugrad(nstate,msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(uh(nstate,msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cgrad(nstate,msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
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
            WRITE(6,'(1X,"==",T15,A,T64,"==")')&
            '   CONJUGATED GRADIENTS OPTIMIZATION '
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T16,A,T64,"==")')&
            ' TAMM-DANCOFF SUBSPACE APPROXIMATION '
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    DO iroot = 1, nroot
       IF ( paral%parent .AND. kprint.GE.0 ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T25,A,I2,T64,"==")')&
               '   OPTIMIZING ROOT ', iroot
       ENDIF
       ! ..initilize vectors
       CALL dcopy(2*ncpw%ngw*msub,c1(:,:,iroot),1,cb,1)
       CALL lr_orthos(nstate,msub,c0,cb)
       IF (td03%tdacanon) THEN
          DO i = 1, iroot-1
             so=dotps(ncpw%ngw,msub,c1(:,:,i),cb)
             CALL mp_sum(so,parai%allgrp)
             CALL daxpy(2*ncpw%ngw*msub,-so,c1(:,:,i),1,cb,1)
          ENDDO
       ELSE
          CALL get_u(urot(:,:,iroot),umat(:,:,iroot),pmat,peig,nstate,&
               msub)
          CALL state_ortho(c1,cb,c2,umat,iroot,nstate,msub)
       ENDIF
       so=dotps(ncpw%ngw,msub,cb,cb)
       CALL mp_sum(so,parai%allgrp)
       CALL dscal(2*ncpw%ngw*msub,1._real_8/SQRT(so),cb,1)
       teigold=1.e30_real_8
       iroteff=1
       DO irot=1,td01%rotit
          IF ((paral%parent.AND. kprint.GE.0 ).AND.paral%io_parent)&
               WRITE(6,'(T25,A)')&
               'Iter            Residual             TCPU'
          IF (.NOT.td03%tdacanon) THEN
             CALL get_u(urot(:,:,iroot),umat(:,:,iroot),pmat,peig,nstate,&
                  msub)
             CALL dgemm('N','N',2*ncpw%ngw,msub,nstate,1._real_8,c0,2*ncpw%ngw,umat(:,:,&
                  iroot),nstate,0._real_8,sc0,2*ncpw%ngw)
             IF (INDEX(orbital,"CANON").NE.0) THEN
                DO i=1,msub
                   DO j=1,nstate
                      ufmat(i,j)=-umat(j,i,iroot)*eigv(j)
                   ENDDO
                ENDDO
             ELSE
                CALL dgemm('T','N',msub,nstate,nstate,1._real_8,umat(:,:,iroot)&
                     ,nstate,eigv,nstate,0._real_8,ufmat,msub)
             ENDIF
             CALL dgemm('N','N',msub,msub,nstate,1._real_8,ufmat,msub,umat(1,&
                  1,iroot),nstate,0._real_8,fmat,msub)
          ENDIF
          t1=m_walltime()
          DO iter=1,td01%tdpcgiter
             id=nstate-msub+1
             IF (td03%tdacanon) THEN
                CALL lr_force_sub(c0(:,id:),cb,cs,sc0,e2,ddxc,psi,drhoe,&
                     eigv(id),nstate,msub,"TDA",orbital)
             ELSE
                CALL lr_force_sub(sc0,cb,cs,c2,e2,ddxc,psi,drhoe,fmat,&
                     nstate,msub,"TDA","GENERAL")
             ENDIF
             CALL dscal(2*ncpw%ngw*msub,-1._real_8,cs,1)
             CALL lr_orthos(nstate,msub,c0,cs)
             teig(iroot)=dotps(ncpw%ngw,msub,cb,cs)
             CALL mp_sum(teig(iroot),parai%allgrp)
             CALL daxpy(2*ncpw%ngw*msub,-teig(iroot),cb,1,cs,1)
             CALL dcopy(2*ncpw%ngw*msub,cs,1,c2,1)
             IF (td03%tdacanon) THEN
                DO i = 1, iroot-1
                   so=dotps(ncpw%ngw,msub,c1(:,:,i),c2(:,:))
                   CALL mp_sum(so,parai%allgrp)
                   CALL daxpy(2*ncpw%ngw*msub,-so,c1(:,:,i),1,c2(:,:),1)
                ENDDO
             ELSE
                CALL state_ortho(c1,c2,cs,umat,iroot,nstate,msub)
             ENDIF
             ores=dotps(ncpw%ngw,msub,c2,c2)
             CALL mp_sum(ores,parai%allgrp)
             CALL dcopy(2*ncpw%ngw*msub,cb,1,c1(:,:,iroot),1)
             IF (paral%parent.AND.td03%tdacanon) THEN
                t2=m_walltime()
                tcpu=(t2-t1)*1.e-3_real_8
                IF (paral%io_parent)&
                     WRITE(6,'(T25,I4,5X,G15.6,5X,F12.2)') iter,ores,tcpu
             ENDIF
             IF (ores.LT.td02%epstdav) GOTO 1000
             IF (soft_com%exsoft) GOTO 1001
             t1=m_walltime()
             ! Precondtitoning of gradients
             DO is=1,msub
                id=nstate-msub+is
                DO ig=1,ncpw%ngw
                   c2(ig,is)=c2(ig,is)/(vpp(ig)-eigv(id))
                ENDDO
             ENDDO
             ! Update
             dgg=dotps(ncpw%ngw,msub,c2,c2)
             CALL mp_sum(dgg,parai%allgrp)
             IF (iter.EQ.1) THEN
                CALL dcopy(2*ncpw%ngw*msub,c2,1,ch,1)
             ELSE
                gam=dgg/gg
                CALL daxpy(2*ncpw%ngw*msub,gam,ch,1,c2,1)
                CALL dcopy(2*ncpw%ngw*msub,c2,1,ch,1)
             ENDIF
             gg=dgg
             IF (td03%tdpcgmin) THEN
                lam=td02%tdpcgstep
                IF (paral%io_parent)&
                     WRITE(6,*) "WARNING:not implemented"
             ELSE
                lam=td02%tdpcgstep
             ENDIF
             CALL daxpy(2*ncpw%ngw*msub,-lam,ch,1,cb,1)
             CALL lr_orthos(nstate,msub,c0,cb(:,:))
             IF (td03%tdacanon) THEN
                DO i = 1, iroot-1
                   so=dotps(ncpw%ngw,msub,c1(:,:,i),cb(:,:))
                   CALL mp_sum(so,parai%allgrp)
                   CALL daxpy(2*ncpw%ngw*msub,-so,c1(:,:,i),1,cb,1)
                ENDDO
             ELSE
                CALL state_ortho(c1,cb,cs,umat,iroot,nstate,msub)
             ENDIF
             so=dotps(ncpw%ngw,msub,cb,cb)
             CALL mp_sum(so,parai%allgrp)
             CALL dscal(2*ncpw%ngw*msub,1._real_8/SQRT(so),cb,1)
          ENDDO
1000      CONTINUE
          IF (paral%parent.AND..NOT.td03%tdacanon.AND. kprint.GE.0 ) THEN
             t2=m_walltime()
             tcpu=(t2-t1)*1.e-3_real_8
             IF (paral%io_parent)&
                  WRITE(6,'(T25,I4,5X,G15.6,5X,F12.2)') iter,ores,tcpu
          ENDIF
          IF (td03%tdacanon) GOTO 1002
          IF (teig(iroot).LT.teigold) THEN
             teigold=teig(iroot)
             CALL zeroing(ugrad)!,nstate*msub)
             CALL zeroing(cgrad)!,nstate*msub)
             CALL lru_force(c0,c1(:,:,iroot),c2,cgrad,psi(:,1),drhoe,ufmat,&
                  nstate,msub)
             CALL dscal(nstate*msub,-1._real_8,cgrad,1)
             CALL xgrad(urot(:,:,iroot),ugrad,cgrad,pmat,peig,nstate,&
                  msub)
             so=ddot(nstate*msub,ugrad,1,ugrad,1)
             IF (paral%parent.AND. kprint.GE.0 ) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,I4,T48,A)') ' STEP :',irot,&
                     ' SUBSPACE ROTATION'
                IF (paral%io_parent)&
                     WRITE(6,'(A,G15.6,T31,A,F12.4)') ' GRADIENT',so,&
                     ' EXCITATION ENERGY (eV)',2._real_8*teig(iroot)*ry
             ENDIF
             IF (so.LT.td02%epsrot) GOTO 1002
             IF (irot.NE.1.AND.iter.EQ.1) GOTO 1002
             IF (soft_com%exsoft) GOTO 1001
             IF (paral%parent) THEN
                rg=so
                IF (iroteff.EQ.1) THEN
                   CALL dcopy(nstate*msub,ugrad,1,uh,1)
                ELSE
                   gam=rg/rdg
                   CALL daxpy(nstate*msub,gam,uh,1,ugrad,1)
                   CALL dcopy(nstate*msub,ugrad,1,uh,1)
                ENDIF
                rdg=rg
                CALL daxpy(nstate*msub,-td02%rotstep,uh,1,urot(:,:,iroot),1)
             ENDIF
             iroteff=iroteff+1
          ELSE
             IF (soft_com%exsoft) GOTO 1001
             IF (paral%parent.AND. kprint.GE.0 ) THEN
                so=0._real_8
                IF (paral%io_parent)&
                     WRITE(6,'(A,I4,T48,A)') ' STEP :',irot,&
                     ' SUBSPACE ROTATION'
                IF (paral%io_parent)&
                     WRITE(6,'(A,G15.6,T31,A,F12.4)') ' GRADIENT',so,&
                     ' EXCITATION ENERGY (eV)',2._real_8*teig(iroot)*ry
                IF (paral%io_parent)&
                     WRITE(6,'(A,I4,T48,A)') ' STEP BACK, REINITIALIZE PCG '
                IF (iroteff.EQ.1) THEN
                   rotb=0.5_real_8*rotb
                ELSE
                   rotb=0.5_real_8*td02%rotstep
                ENDIF
                CALL daxpy(nstate*msub,rotb,uh,1,urot(:,:,iroot),1)
             ENDIF
             iroteff=1
          ENDIF
          CALL mp_bcast(urot(:,:,iroot),nstate*msub,parai%source,parai%allgrp)
       ENDDO
1002   CONTINUE
    ENDDO
1001 CONTINUE
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ch,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (.NOT.td03%tdacanon) THEN
       DEALLOCATE(fmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(umat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ufmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ugrad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(uh,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(pmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(peig,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cgrad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_pcg
  ! ==================================================================
  SUBROUTINE give_scr_td_pcg(ltd_pcg,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltd_pcg
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: llr_force, llr_ortho, llr_sub

    CALL give_scr_lr_force(llr_force,"TDA",tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    IF (td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon) THEN
       CALL give_scr_lru_force(llr_sub,td01%msubta,tag)
    ELSE
       llr_sub=0
    ENDIF
    ltd_pcg=MAX(llr_force,llr_ortho,llr_sub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_td_pcg
  ! ==================================================================
  SUBROUTINE state_ortho(c1,cb,cs,umat,iroot,nstate,msub)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cb(:,:), cs(:,:)
    INTEGER                                  :: iroot, nstate, msub
    REAL(real_8)                             :: umat(nstate,msub,*)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,msub,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'state_ortho'

    INTEGER                                  :: ierr, ir
    REAL(real_8)                             :: so
    REAL(real_8), ALLOCATABLE                :: aux(:)

    ALLOCATE(aux(msub*msub),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO ir=1,iroot-1
       CALL dgemm('T','N',msub,msub,nstate,1._real_8,umat(:,:,ir),nstate,&
            umat(:,:,iroot),nstate,0._real_8,aux,msub)
       CALL dgemm('N','N',2*ncpw%ngw,msub,msub,1._real_8,c1(:,:,ir),2*ncpw%ngw,aux,&
            msub,0._real_8,cs,2*ncpw%ngw)
       so=dotps(ncpw%ngw,msub,cs,cb)
       CALL mp_sum(so,parai%allgrp)
       CALL daxpy(2*ncpw%ngw*msub,-so,cs,1,cb,1)
    ENDDO
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE state_ortho
  ! ==================================================================

END MODULE td_pcg_utils
