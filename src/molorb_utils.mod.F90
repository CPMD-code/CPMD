MODULE molorb_utils
  USE adat,                            ONLY: covrad
  USE cnst,                            ONLY: fbohr
  USE empf,                            ONLY: c,&
                                             ibind,&
                                             naty,&
                                             zan
  USE empfor_utils,                    ONLY: cpmdbonds,&
                                             smolec
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc3
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parm
  USE wann,                            ONLY: wannr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: molorb

CONTAINS

  ! ==================================================================
  SUBROUTINE molorb(c0,c2,tau0,nstate,centers)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: centers(4,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'molorb'

    INTEGER                                  :: i, i1, i2, iat, ierr, info, &
                                                j, k, laux, lcenters, lhmat, &
                                                lnwann, lwmol, minat, n1, n2, &
                                                ncc, nma, nmpa1, nst, nummol
    INTEGER, ALLOCATABLE                     :: imol(:), nbonds(:), &
                                                nwann(:,:), wmol(:)
    LOGICAL                                  :: debug
    REAL(real_8)                             :: d1, d2, d3, dist, fac, &
                                                mindist, rcov(0:104), xmfac
    REAL(real_8), ALLOCATABLE                :: aux(:), hmat(:,:), &
                                                hmat2(:,:), molchg(:), &
                                                rotmat(:,:), w(:)

    debug=.FALSE.

    lcenters=4*nstate
    lwmol=nstate/2+1
    ALLOCATE(wmol(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    lnwann=(ions1%nat*clsd%nlsd)/2+1
    ALLOCATE(nwann(ions1%nat, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)


    IF (paral%parent) THEN
       nma=ions1%nat
       nmpa1=(ions1%nat*(ions1%nat+1))/2

       ! TODO check this 
       ALLOCATE(zan(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ibind(nmpa1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(nbonds(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(imol(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(molchg(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL zeroing(imol)!,ions1%nat)
       CALL zeroing(ibind)!,nmpa1)
       CALL zeroing(nbonds)!,ions1%nat)
       naty=ions1%nat
       k=0
       DO i=1,ions1%nsp
          DO j=1,ions0%na(i)
             k=k+1
             zan(k)=REAL(ions0%iatyp(i),kind=real_8)
             c(1,k)=tau0(1,j,i)-wannr%w_ref(1)
             c(2,k)=tau0(2,j,i)-wannr%w_ref(2)
             c(3,k)=tau0(3,j,i)-wannr%w_ref(3)
          ENDDO
       ENDDO
       DO i=1,99
          rcov(i)=covrad(i)*fbohr
       ENDDO
       ncc=3*ions1%nat
       xmfac=1.35_real_8
       CALL cpmdbonds(rcov,xmfac,nbonds)

       ! Search for molecules
       CALL smolec(imol,nummol)

       DO i=1,nstate
          d1=centers(1,i)-c(1,1)
          d2=centers(2,i)-c(2,1)
          d3=centers(3,i)-c(3,1)
          CALL pbc3(d1,d2,d3,d1,d2,d3,1,parm%apbc,parm%ibrav)
          mindist = d1*d1 + d2*d2 + d3*d3
          minat=1
          DO j=2,naty
             d1=centers(1,i)-c(1,j)
             d2=centers(2,i)-c(2,j)
             d3=centers(3,i)-c(3,j)
             CALL pbc3(d1,d2,d3,d1,d2,d3,1,parm%apbc,parm%ibrav)
             dist = d1*d1 + d2*d2 + d3*d3
             IF (dist.LT.mindist) THEN
                minat=j
                mindist=dist
             ENDIF
          ENDDO
          wmol(i)=imol(minat)
       ENDDO

       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '          ORB.      MOLECULE '
       DO i=1,nstate
          IF (paral%io_parent)&
               WRITE(6,*) i,wmol(i)
       ENDDO


       CALL zeroing(nwann)!,ions1%nat*clsd%nlsd)
       CALL zeroing(molchg)!,ions1%nat)
       IF (cntl%tlsd) THEN
          fac=1.0_real_8
          DO i=1,spin_mod%nsup
             j=wmol(i)
             nwann(j,1)=nwann(j,1)+1
             molchg(j)=molchg(j)-fac
          ENDDO
          DO i=spin_mod%nsup+1,nstate
             j=wmol(i)
             nwann(j,2)=nwann(j,2)+1
             molchg(j)=molchg(j)-fac
          ENDDO
       ELSE
          fac=2.0_real_8
          DO i=1,nstate
             j=wmol(i)
             nwann(j,1)=nwann(j,1)+1
             molchg(j)=molchg(j)-fac
          ENDDO
       ENDIF

       iat=0
       DO i=1,ions1%nsp
          DO j=1,ions0%na(i)
             iat=iat+1
             k=imol(iat)
             molchg(k)=molchg(k)+ions0%zv(i)
          ENDDO
       ENDDO

       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '       MOLECULE      CHARGE '
       DO i=1,nummol
          IF (paral%io_parent)&
               WRITE(6,*) i,molchg(i)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)

       DEALLOCATE(zan,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ibind,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       DEALLOCATE(nbonds,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(imol,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(molchg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

    ENDIF

    CALL mp_bcast(centers,lcenters,parai%source,parai%allgrp)
    CALL mp_bcast(wmol,lwmol,parai%source,parai%allgrp)
    CALL mp_bcast(nwann,lwmol,parai%source,parai%allgrp)
    CALL mp_bcast(nummol,parai%source,parai%allgrp)



    lhmat=nstate*nstate
    laux=12*nstate
    ALLOCATE(hmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hmat2(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rotmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(12*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL ovlap(nstate,hmat,c0,c2)
    CALL mp_sum(hmat,nstate*nstate,parai%allgrp)

    CALL zeroing(rotmat)!,nstate*nstate)
    IF (cntl%tlsd) THEN
       i1=1
       i2=spin_mod%nsup+1
       DO i=1,nummol
          DO j=1,spin_mod%nsup
             IF (wmol(j).EQ.i) THEN
                rotmat(j,i1)=1.0_real_8
                i1=i1+1
             ENDIF
          ENDDO
          DO j=spin_mod%nsup+1,nstate
             IF (wmol(j).EQ.i) THEN
                rotmat(j,i2)=1.0_real_8
                i2=i2+1
             ENDIF
          ENDDO
       ENDDO
    ELSE
       i1=1
       DO i=1,nummol
          DO j=1,nstate
             IF (wmol(j).EQ.i) THEN
                rotmat(j,i1)=1.0_real_8
                i1=i1+1
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    CALL rotate(1.0_real_8,c0,0.0_real_8,c2,rotmat,nstate,&
         2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    ! debug      CALL DCOPY(2*NGW*NSTATE,C2,1,C0,1)

    IF (cntl%tlsd) THEN
       n2=spin_mod%nsup+1
       CALL dgemm('N','N',spin_mod%nsup,spin_mod%nsup,spin_mod%nsup,1.0_real_8,hmat,nstate,rotmat,&
            nstate,0.0_real_8,hmat2,nstate)
       CALL dgemm('T','N',spin_mod%nsup,spin_mod%nsup,spin_mod%nsup,1.0_real_8,rotmat,nstate,hmat2,&
            nstate,0.0_real_8,hmat,nstate)

       CALL dgemm('N','N',spin_mod%nsdown,spin_mod%nsdown,spin_mod%nsdown,1.0_real_8,hmat(n2,n2),&
            nstate,rotmat(n2,n2),nstate,0.0_real_8,hmat2(n2,n2),nstate)
       CALL dgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,spin_mod%nsdown,1.0_real_8,rotmat(n2,n2),&
            nstate,hmat2(n2,n2),nstate,0.0_real_8,hmat(n2,n2),nstate)
    ELSE
       CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,hmat,nstate,&
            rotmat,nstate,0.0_real_8,hmat2,nstate)
       CALL dgemm('T','N',nstate,nstate,nstate,1.0_real_8,rotmat,nstate,&
            hmat2,nstate,0.0_real_8,hmat,nstate)
    ENDIF


    IF (cntl%tlsd) THEN
       n1=1
       n2=spin_mod%nsup+1
       DO i=1,nummol
          nst=nwann(i,1)
          CALL dsyev('V','U',nst,hmat(n1,n1),nstate,w,aux,laux,info)
          CALL dgemm('N','N',2*ncpw%ngw,nst,nst,1.0_real_8,c2(1,n1),2*ncpw%ngw,&
               hmat(n1,n1),nstate,0.0_real_8,c0(1,n1),2*ncpw%ngw)
          n1=n1+nst

          nst=nwann(i,2)
          CALL dsyev('V','U',nst,hmat(n2,n2),nstate,w,aux,laux,info)
          CALL dgemm('N','N',2*ncpw%ngw,nst,nst,1.0_real_8,c2(1,n2),2*ncpw%ngw,&
               hmat(n2,n2),nstate,0.0_real_8,c0(1,n2),2*ncpw%ngw)
          n2=n2+nst
       ENDDO
    ELSE
       n1=1
       DO i=1,nummol
          nst=nwann(i,1)
          CALL dsyev('V','U',nst,hmat(n1,n1),nstate,w,aux,laux,info)
          CALL dgemm('N','N',2*ncpw%ngw,nst,nst,1.0_real_8,c2(1,n1),2*ncpw%ngw,&
               hmat(n1,n1),nstate,0.0_real_8,c0(1,n1),2*ncpw%ngw)
          n1=n1+nst
       ENDDO
    ENDIF


    DEALLOCATE(hmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hmat2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rotmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE molorb
  ! ==--------------------------------------------------------------==

END MODULE molorb_utils
