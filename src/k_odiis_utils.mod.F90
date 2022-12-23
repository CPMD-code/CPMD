MODULE k_odiis_utils
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !public :: k_odiis
  !public :: k_grimax

  !contains

END MODULE k_odiis_utils


! ==================================================================
SUBROUTINE k_odiis(c0,c2,vpp,ik,nstate,pme,gde,svar2,&
     reinit,calcrho)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum,mp_bcast
  USE system , ONLY:cnti,cntl,maxdis,nkpt
  USE parac, ONLY : paral,parai
  USE ener , ONLY:ener_com
  USE elct , ONLY:crge
  USE geq0mod , ONLY:geq0
  USE odiis_utils, ONLY : updis, solve
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  COMPLEX(real_8)                            :: c0(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
  REAL(real_8)                               :: vpp(*)
  INTEGER                                    :: ik, nstate
  REAL(real_8) :: pme(nkpt%ngwk*nstate+8,cnti%mdiis,*), &
      gde((nkpt%ngwk*nstate+8)/4,cnti%mdiis,*), svar2
  LOGICAL                                    :: reinit, calcrho

  CHARACTER(*), PARAMETER                    :: procedureN = 'k_odiis'

  EXTERNAL                                   :: ddot
  INTEGER                                    :: i, ierr, ig, isub, j, k, nowm
  INTEGER, SAVE                              :: ifirst = 0, istate = 0, &
                                                ndiis, nowv, nsize
  REAL(real_8)                               :: bc(maxdis+1,maxdis+1), ddot, &
                                                de1, ff, g1, g2, raim, ratio, &
                                                ration, vc(maxdis+1)
  REAL(real_8), ALLOCATABLE, SAVE            :: diism(:,:,:), GAMMA(:), &
                                                gimax(:,:), gnorm(:,:), &
                                                grmax(:,:)
  REAL(real_8), SAVE                         :: eold

  CALL tiset(procedureN,isub)

  IF (ifirst .EQ. 0 ) THEN
     ALLOCATE(diism(maxdis,maxdis,nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(grmax(maxdis,nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(gimax(maxdis,nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(gnorm(maxdis,nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(GAMMA(nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ifirst = 1
     nowv=0
  ENDIF
  IF (reinit.OR.istate.EQ.0 .AND. ik .EQ. 1 ) THEN
     ! Reinitialization of cntl%diis procedure
     istate=0
     eold=9999._real_8
     CALL zeroing(gamma)!,nkpt%nkpnt)
     reinit=.FALSE.
     ndiis=0
  ENDIF
  ! Empirical rules to assure convergence
  IF (istate.NE.0) THEN
     ! Check for an increase in energy

     de1=ener_com%etot-eold
     IF (calcrho) THEN
        IF (geq0) THEN
           IF (de1.GT.0.0_real_8) THEN
              GAMMA(ik)=GAMMA(ik)+1._real_8/vpp(1)
           ELSE
              GAMMA(ik)=0.75_real_8*GAMMA(ik)
           ENDIF
        ENDIF
     ENDIF
     CALL mp_bcast(GAMMA(ik),parai%io_source,parai%cp_grp)
     reinit=.FALSE.
     IF (ndiis.GT.cnti%nreset.AND.cnti%nreset.GE.3) THEN
        ! ..restart if there was no progress over the last steps 
        nowm=nowv+1
        IF (nowv.EQ.cnti%mdiis) nowm=1
        g1=SQRT(grmax(nowv,ik)**2 + gimax(nowv,ik)**2)
        g2=SQRT(grmax(nowm,ik)**2 + gimax(nowm,ik)**2)
        ratio=g1/g2
        ration=gnorm(nowv,ik)/gnorm(nowm,ik)
        raim=1.0_real_8-REAL(cnti%mdiis,kind=real_8)/100._real_8
        IF (ratio.GT.raim.AND.ration.GT.raim .AND. g1.GT.1.e-5_real_8)&
             reinit=.TRUE.
        IF ((reinit .AND. paral%parent).AND.paral%io_parent)&
             WRITE(6,'(A,I8,A)')&
             ' ODIIS| Insufficient progress for K = ',ik,'; reset! '
        IF (reinit) THEN
           istate=0
           GAMMA(ik)=0.0_real_8
           reinit=.FALSE.
        ENDIF
        IF (calcrho) THEN
           IF (de1 .GE. 1.0_real_8 .AND. ik .EQ. 1) THEN
              ! Reinitialization of cntl%diis procedure
              istate=0
              eold=9999._real_8
              CALL zeroing(gamma)!,nkpt%nkpnt)
              reinit=.FALSE.
              ndiis=0
           ELSEIF (ABS(de1) .LE. 1.0e-8_real_8 .AND. ik .EQ. 1) THEN
              istate=0
              eold=9999._real_8
              CALL zeroing(gamma)!,nkpt%nkpnt)
              reinit=.TRUE.
              ndiis=0
           ENDIF
           IF ((reinit .AND. paral%parent).AND.paral%io_parent)&
                WRITE(6,'(A)')&
                ' ODIIS| Insufficient progress; reset! '
        ENDIF
     ENDIF
  ENDIF

  IF (istate.EQ.0 .AND. ik .EQ. 1)&
       CALL updis(ndiis,nowv,nsize,cnti%mdiis,0)
  IF (ik .EQ. 1)      istate=2
  ! Perform an electronic cntl%diis step
  IF (ik .EQ. 1) CALL updis(ndiis,nowv,nsize,cnti%mdiis,1)

  DO ig=1,nkpt%ngwk
     vpp(ig)=vpp(ig)/(1.0_real_8+vpp(ig)*GAMMA(ik))
  ENDDO
  ! Update cntl%diis buffers
  CALL dscopy(nkpt%ngwk*nstate,c0,pme(1,nowv,ik))
  CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,c2(1,1),1)
  CALL k_grimax(c2,nstate,grmax(nowv,ik),gimax(nowv,ik),&
       gnorm(nowv,ik))
  CALL dicopy(nkpt%ngwk*nstate,c2,gde(1,nowv,ik),&
       grmax(nowv,ik),gimax(nowv,ik))
  ! Scale Gradient with Preconditioner(**2)
  DO k=1,nstate
     IF (cntl%prec .AND. crge%f(k,ik).GT.0.001_real_8) THEN
        ff=1.0_real_8/crge%f(k,ik)
     ELSE
        ff=1.0_real_8
     ENDIF
     !$omp parallel do private(IG) shared(FF)
     DO ig=1,nkpt%ngwk
        c2(ig,k)=c2(ig,k)*2.0_real_8*ff*ff*vpp(ig)*vpp(ig)
     ENDDO
  ENDDO
  ! Update cntl%diis matrix
  DO i=1,nsize-1
     diism(i,nowv,ik)=0.0_real_8
  ENDDO
  DO i=1,nsize-1
     CALL idcopy(nkpt%ngwk*nstate,c0,gde(1,i,ik),grmax(i,ik),gimax(i,ik))
     DO k=1,nstate
        diism(i,nowv,ik)=diism(i,nowv,ik)+&
             ddot(2*nkpt%ngwk,c0(1,k),1,c2(1,k),1)
     ENDDO
  ENDDO
  CALL mp_sum(diism(:,nowv,ik),nsize-1,parai%allgrp)
  DO i=1,nsize-1
     diism(nowv,i,ik)=diism(i,nowv,ik)
  ENDDO
  ! Set up cntl%diis Matrix
  CALL zeroing(bc)!,(maxdis+1)*(maxdis+1))
  !$omp parallel do private(I,J)
  DO i=1,nsize-1
     DO j=1,nsize-1
        bc(i,j)=diism(i,j,ik)
     ENDDO
  ENDDO
  DO i=1,nsize-1
     vc(i)=0._real_8
     bc(i,nsize)=-1._real_8
     bc(nsize,i)=-1._real_8
  ENDDO
  vc(nsize)=-1._real_8
  bc(nsize,nsize)=0._real_8
  ! Solve System of Linear Equations
  CALL solve(bc,maxdis+1,nsize,vc)
#if defined(__VECTOR)
  ! Compute Interpolated Coefficient Vectors
  CALL zeroing(c0(:,1:nstate))!,nkpt%ngwk*nstate)
  DO i=1,nsize-1
     CALL sdcopy(nkpt%ngwk*nstate,pme(1,i,ik),c2)
     CALL daxpy(2*nkpt%ngwk*nstate,vc(i),c2(1,1),1,c0(1,1),1)
  ENDDO
#else
  ! do the above in one go and avoid having to use C2.
  CALL intpcoef(nkpt%ngwk*nstate,nsize,nkpt%ngwk*nstate+8,pme(1,1,ik),vc,c0)
#endif
  ! Estimate New Parameter Vectors 
  DO i=1,nsize-1
     CALL idcopy(nkpt%ngwk*nstate,c2,gde(1,i,ik),grmax(i,ik),gimax(i,ik))
     ff=1.0_real_8
     DO k=1,nstate
        IF (cntl%prec.AND.crge%f(k,ik).GT.0.001_real_8) THEN
           ff=1.0_real_8/crge%f(k,ik)
        ELSE
           ff=1.0_real_8
        ENDIF
        !$omp parallel do private(IG) shared(FF)
        DO ig=1,nkpt%ngwk
           c0(ig,k)=c0(ig,k)-vc(i)*ff*vpp(ig)*c2(ig,k)
        ENDDO
     ENDDO
  ENDDO
  IF ( ik .EQ. nkpt%nkpnt ) eold=ener_com%etot
  ! DO IG=1,NGWK
  ! VPP(IG)=1._real_8/(1.0_real_8/VPP(IG)-GAMMA)
  ! ENDDO
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE k_odiis
! ==================================================================


! ==================================================================
SUBROUTINE k_grimax(c2,n,grmax,gimax,gnorm)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum,mp_max
  USE system , ONLY:nkpt,spar
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: n
  REAL(real_8)                               :: c2(2,nkpt%ngwk*n), grmax, &
                                                gimax, gnorm

  EXTERNAL                                   :: ddot
  INTEGER                                    :: i, idamax, ii, ir
  REAL(real_8)                               :: ddot

! ==--------------------------------------------------------------==

  ir=idamax(n*nkpt%ngwk,c2(1,1),2)
  ii=idamax(n*nkpt%ngwk,c2(2,1),2)
  grmax=c2(1,ir)
  gimax=c2(2,ii)
  gnorm=0.0_real_8
  DO i=1,n
     ii=(i-1)*nkpt%ngwk+1
     gnorm=gnorm+ddot(2*nkpt%ngwk,c2(1,ii),1,c2(1,ii),1)
  ENDDO
  CALL mp_sum(gnorm,parai%allgrp)
  grmax=ABS(grmax);gimax=ABS(gimax)
  CALL mp_max(grmax,parai%allgrp)
  CALL mp_max(gimax,parai%allgrp)
  gnorm=SQRT(gnorm/(n*spar%ngwks))
  RETURN
END SUBROUTINE k_grimax
! ==================================================================
