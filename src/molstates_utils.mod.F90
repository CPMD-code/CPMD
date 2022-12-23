MODULE molstates_utils
  USE cnst,                            ONLY: ry
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: wcent
  USE localize_utils,                  ONLY: localize
  USE mm_input,                        ONLY: lqmmm
  USE mols,                            ONLY: msnum,&
                                             mstat,&
                                             numol
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc3
  USE poin,                            ONLY: potr
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE utils,                           ONLY: nxxfun,&
                                             rmatmov
  USE wann,                            ONLY: wannl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: molecular_states

CONTAINS

  ! ==================================================================
  SUBROUTINE molecular_states(c0,eigv,numorb,tau0)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: eigv(*)
    INTEGER                                  :: numorb
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,numorb)
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'molecular_states'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), psi(:), sc0(:,:)
    INTEGER                                  :: i, ia, iat, ierr, il_psi, &
                                                info, is, j, k, lscr, nxx
    INTEGER, ALLOCATABLE                     :: order(:,:)
    REAL(real_8)                             :: d1, d2, d3, dd, dm
    REAL(real_8), ALLOCATABLE                :: hmat(:,:), scr(:), umat(:)

    IF (cntl%tlsd) CALL stopgm('MOLSTATE','TLSD NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lqmmm%qmmm) CALL stopgm('MOLSTATE','QMMM NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    wannl%twann=.TRUE.
    CALL give_scr_ddipo(lscr,tag)
    lscr=MAX(lscr,2*ncpw%nhg)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    lenbk=nxxfun(numorb)
    nxx=MAX(2*lenbk*parai%nproc,2*nkpt%ngwk*numorb)
    ALLOCATE(sc0(nkpt%ngwk,nxx/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cs(nkpt%ngwk,nxx/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL localize(tau0,c0,cs,sc0,numorb)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(order(2,numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       DO i=1,numorb
          dm=1.e30_real_8
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                d1=wcent(1,i)-tau0(1,ia,is)
                d2=wcent(2,i)-tau0(2,ia,is)
                d3=wcent(3,i)-tau0(3,ia,is)
                CALL pbc3(d1,d2,d3,d1,d2,d3,1,parm%apbc,parm%ibrav)
                dd=SQRT(d1*d1+d2*d2+d3*d3)
                IF (dd.LT.dm) THEN
                   k=iat
                   dm=dd
                ENDIF
             ENDDO
          ENDDO
          order(1,i)=msnum(k)
          IF (paral%io_parent)&
               WRITE(6,'(A,I7,A,I7,A,I7)') '  STATE:',i,&
               '           ATOM:',k,'           MOLECULE:',msnum(k)
       ENDDO
       k=0
       DO i=1,numol
          is=0
          DO j=1,numorb
             IF (i.EQ.order(1,j)) is=is+1
          ENDDO
          IF (is.NE.mstat(i)) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' WARNING: NUMBER OF STATES IS INCONSISTENT'
             IF (paral%io_parent)&
                  WRITE(6,*) ' EXPECTED:',mstat(i),'   FOUND:',is
          ENDIF
          DO j=1,numorb
             IF (order(1,j).EQ.i) THEN
                k=k+1
                order(2,j)=k
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    CALL mp_bcast(order,SIZE(order),parai%source,parai%allgrp)
    DO i=1,numorb
       CALL dcopy(2*ncpw%ngw,c0(1,i),1,cs(1,order(2,i)),1)
    ENDDO
    ! 
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_hpsi(lscr,tag,numorb)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rhoe_psi_size(i,il_psi)
    ALLOCATE(psi(il_psi),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL hpsi(cs,c0,sc0,potr,psi,numorb,1,1)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hmat(numorb,numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(umat(numorb*numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL ovlap(numorb,hmat,c0,cs)
    CALL mp_sum(hmat,numorb*numorb,parai%allgrp)
    CALL dscal(numorb*numorb,-1._real_8,hmat,1)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(T30,A)') 'MOLECULAR STATES'
    ENDIF
    is=1
    DO i=1,numol
       ia=mstat(i)
       IF (paral%parent) THEN
          CALL rmatmov(ia,ia,hmat(is,is),numorb,umat,ia)
          CALL dsyev('V','U',ia,umat,ia,eigv(is),scr,lscr,info)
          IF (info.NE.0) CALL stopgm('MOLSTAT','INFO! = 0 AFTER DSYEV',& 
               __LINE__,__FILE__)
          IF (paral%io_parent)&
               WRITE(6,'(A,I6,A)') ' MOLECULE:',i,'    EIGENVALUES [eV]'
          IF (paral%io_parent)&
               WRITE(6,'(T16,4(F11.4,2X))') (eigv(is+j-1)*2._real_8*ry,j=1,ia)
       ENDIF
       CALL mp_bcast(umat,ia*ia,parai%source,parai%allgrp)
       CALL rotate(1._real_8,cs(:,is:),0._real_8,c0(:,is:),umat,ia,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
            spin_mod%nsdown)
       is=is+mstat(i)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    CALL mp_bcast(eigv,numorb,parai%source,parai%allgrp)
    ! 
    DEALLOCATE(cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(umat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(order,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE molecular_states
  ! ==================================================================

END MODULE molstates_utils
