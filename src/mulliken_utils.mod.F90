MODULE mulliken_utils
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             loadc_foc_array_size,&
                                             xsmat,&
                                             xxmat
  USE augchg_utils,                    ONLY: augchg
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop2
  USE pslo,                            ONLY: pslo_com
  USE setbasis_utils,                  ONLY: loadc
  USE sfac,                            ONLY: fnl
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mulliken
  PUBLIC :: give_scr_mulliken

CONTAINS

  ! ==================================================================
  SUBROUTINE mulliken(ityp,c0,nstate,tau0,achrg)
    ! ==--------------------------------------------------------------==
    ! == PROJECT WAVEFUNCTIONS ON PSEUDO ATOMIC ORBITALS              ==
    ! == AND CALCULATE MULLIKEN POPULATION ANALYSIS                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ityp
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(*), achrg(ions1%nat)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mulliken'

    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cscr(:,:)
    INTEGER                                  :: i, ia, iao, iaorb, iat, ierr, &
                                                is, isub, ixx, izero, j, k, &
                                                nao, natst, nomax
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: tlsd2
    REAL(real_8)                             :: fac, &
                                                foc(loadc_foc_array_size), &
                                                ql, qm, rlen, sfc
    REAL(real_8), ALLOCATABLE                :: a(:), aux(:), aux2(:), w(:), &
                                                z(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: comp(:), qa(:,:), xtmat(:,:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('  MULLIKEN',isub)
    ! ..Scratch space
    ! ==--------------------------------------------------------------==
    nomax = MAX(nstate,atwp%nattot)
    IF (ifirst.EQ.0) THEN
       ALLOCATE(comp(prop2%numorb),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xtmat(atwp%nattot,nomax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xsmat(atwp%nattot,nomax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xxmat(atwp%nattot,nomax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(catom(ncpw%ngw,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cscr(ncpw%ngw,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(qa(ions1%nat,3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(qa)!,3*ions1%nat)
    ! ..Load AO in PW basis
    iaorb=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL loadc(cscr,foc,ncpw%ngw,ncpw%ngw,atwp%nattot,SIZE(foc),is,iat,natst)
          DO ixx=1,natst
             sfc=dotp(ncpw%ngw,cscr(:,ixx),cscr(:,ixx))
             CALL mp_sum(sfc,parai%allgrp)
             sfc=1._real_8/SQRT(sfc)
             CALL dscal(2*ncpw%ngw,sfc,cscr(1,ixx),1)
          ENDDO
          CALL dcopy(2*ncpw%ngw*natst,cscr,1,catom(1,iaorb),1)
          iaorb=iaorb+natst
       ENDDO
    ENDDO
    ! ..Overlap matrix
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    CALL ovlap(atwp%nattot,xsmat,catom,catom)
    cntl%tlsd=tlsd2
    CALL summat(xsmat,atwp%nattot)
    CALL dcopy(atwp%nattot*atwp%nattot,xsmat,1,xtmat,1)
    ! ..O**(-1/2)
    ! TODO align for BG
    ALLOCATE(a(atwp%nattot*atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(z(atwp%nattot, atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    k=0
    DO j=1,atwp%nattot
       DO i=j,atwp%nattot
          k=k+1
          a(k)=xsmat(i,j)
       ENDDO
    ENDDO
    CALL dspevy(1,a,w,z,atwp%nattot,atwp%nattot,aux,3*atwp%nattot)
    izero=0
    DO i=1,atwp%nattot
       IF (ABS(w(i)).LT.1.e-5_real_8) THEN
          w(i)=1.0e30_real_8
          izero=izero+1
          DO j=1,atwp%nattot
             z(j,i)=0.0_real_8
          ENDDO
       ENDIF
    ENDDO
    IF ((izero.NE.0).AND.paral%io_parent)&
         WRITE(6,'(A,I3)')&
         ' NUMBER OF LINEARLY DEPENDENT ATOMIC BASIS FUNCTIONS :',izero
    DO i=1,atwp%nattot
       w(i)=1._real_8/SQRT(w(i))
    ENDDO
    CALL zeroing(xsmat)!,atwp%nattot*atwp%nattot)
    DO k=1,atwp%nattot
       DO j=1,atwp%nattot
          fac=w(k)*z(j,k)
          CALL daxpy(atwp%nattot,fac,z(1,k),1,xsmat(1,j),1)
       ENDDO
    ENDDO
    ! ..Orthonormalize atomic orbitals
    CALL dgemm('N','N',2*ncpw%ngw,atwp%nattot,atwp%nattot,1.0_real_8,catom,2*ncpw%ngw,xsmat,&
         atwp%nattot,0.0_real_8,cscr,2*ncpw%ngw)
    CALL dcopy(2*atwp%nattot*ncpw%ngw,cscr,1,catom,1)
    ! ..Overlap AO with wavefunctions
    CALL ovlap2(ncpw%ngw,atwp%nattot,prop2%numorb,xxmat,catom,c0,.TRUE.)
    CALL mp_sum(xxmat,atwp%nattot*prop2%numorb,parai%allgrp)
    IF (pslo_com%tivan) THEN
       ! ..Calculate augmentation charges
       CALL augchg(fnl,crge%f,qa(1,3),prop2%numorb)
    ENDIF
    DO i=1,prop2%numorb
       rlen=dotp(ncpw%ngw,c0(:,i),c0(:,i))
       CALL mp_sum(rlen,parai%allgrp)
       comp(i)=ddot(atwp%nattot,xxmat(1,i),1,xxmat(1,i),1)/rlen
    ENDDO
    IF (paral%parent) THEN
       ! ..Lowdin Population Analysis
       iat=0
       iao=0
       DO is=1,ions1%nsp
          nao=atwf_mod%numaor(is)
          DO ia=1,ions0%na(is)
             iat=iat+1
             ! ..QA(.,3) is the true augmentation charge!
             qa(iat,2)=ions0%zv(is)-qa(iat,3)
             DO i=iao+1,iao+nao
                DO j=1,nstate
                   qa(iat,2)=qa(iat,2)-crge%f(j,1)*xxmat(i,j)**2
                ENDDO
             ENDDO
             iao=iao+nao
          ENDDO
       ENDDO
    ENDIF
    ! ..Rotate back to AO basis
    ALLOCATE(aux2(atwp%nattot*prop2%numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL dgemm('N','N',atwp%nattot,prop2%numorb,atwp%nattot,1.0_real_8,xsmat,atwp%nattot,xxmat,&
         atwp%nattot,0.0_real_8,aux2,atwp%nattot)
    CALL dcopy(atwp%nattot*prop2%numorb,aux2,1,xxmat,1)
    IF (paral%parent) THEN
       ! ..Mulliken Population Analysis
       DO i=1,atwp%nattot
          DO j=1,atwp%nattot
             xsmat(i,j)=0.0_real_8
             DO k=1,nstate
                xsmat(i,j)=xsmat(i,j)+crge%f(k,1)*xxmat(i,k)*xxmat(j,k)
             ENDDO
          ENDDO
       ENDDO
       CALL zeroing(aux2)!,atwp%nattot)
       DO j=1,atwp%nattot
          DO i=1,atwp%nattot
             aux2(i)=aux2(i)+xsmat(i,j)*xtmat(i,j)
          ENDDO
       ENDDO
       iat=0
       iao=0
       DO is=1,ions1%nsp
          nao=atwf_mod%numaor(is)
          DO ia=1,ions0%na(is)
             iat=iat+1
             ! ..QA(.,3) is the true augmentation charge!
             qa(iat,1)=ions0%zv(is)-qa(iat,3)
             DO i=iao+1,iao+nao
                qa(iat,1)=qa(iat,1)-aux2(i)
             ENDDO
             iao=iao+nao
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(aux2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ..unassigned charge
    qm=-crge%charge
    ql=-crge%charge
    DO iat=1,ions1%nat
       qm=qm+qa(iat,1)
       ql=ql+qa(iat,2)
    ENDDO
    qm=qm/REAL(ions1%nat,kind=real_8)
    ql=ql/REAL(ions1%nat,kind=real_8)
    ! ..add unassigned charge
    DO iat=1,ions1%nat
       qa(iat,1)=qa(iat,1)-qm
       qa(iat,2)=qa(iat,2)-ql
    ENDDO
    ! ..ITYP=1 --> MULLIKEN
    ! ..ITYP=2 --> LOWDIN
    CALL dcopy(ions1%nat,qa(1,ityp),1,achrg,1)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(z,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  MULLIKEN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mulliken
  ! ==================================================================
  SUBROUTINE give_scr_mulliken(lmulliken,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmulliken
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lmax, lsum, n2

    CALL give_scr_summat(lsum,tag,nstate)
    n2=atwp%nattot*atwp%nattot
    lmax=2*n2+4*atwp%nattot
    lmulliken=MAX(lsum,lmax,atwp%nattot*prop2%numorb)
    tag = 'MAX(LSUM,LMAX,NATTOT*NUMORB)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mulliken
  ! ==================================================================

END MODULE mulliken_utils
