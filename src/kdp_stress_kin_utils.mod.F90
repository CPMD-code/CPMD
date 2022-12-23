MODULE kdp_stress_kin_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             dekin
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kdp_stress_kin

CONTAINS

  SUBROUTINE kdp_stress_kin(nstate,c,gagk,wkdp,xkdp,pkdp,bkdp,nkdp)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c(nkpt%ngwk,nstate)
    REAL(real_8)                             :: gagk(ncpw%nhg,6), &
                                                pkdp(3,nstate,nstate)
    INTEGER                                  :: nkdp
    COMPLEX(real_8)                          :: bkdp(nstate,nstate,nkdp)
    REAL(real_8)                             :: xkdp(3,nkdp), wkdp(nkdp)

    CHARACTER(*), PARAMETER                  :: procedureN = 'kdp_stress_kin'

    COMPLEX(real_8), ALLOCATABLE, SAVE       :: brsum(:,:), brsum3(:,:,:)
    INTEGER                                  :: i, ic, ierr, ig, ikdp, j, kk
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: aux, cgcg, xkkk
    REAL(real_8), ALLOCATABLE, SAVE          :: xkk(:,:)

! 
! ==--------------------------------------------------------------==
! == ALLOCATE ARRAYS FOR KDP ROUTINES                             ==
! ==--------------------------------------------------------------==

    IF (ifirst.EQ.0) THEN
       ALLOCATE(brsum(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(brsum3(nstate,nstate,6),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xkk(6,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    ! 
    ! ==  Initialization
    ! 
    DO kk=1,6
       DO ikdp=1,nkdp
          xkk(kk,ikdp)=xkdp(alpha(kk),ikdp)*xkdp(beta(kk),ikdp)
       ENDDO
       dekin(kk) = 0._real_8
    ENDDO
    ! 
    CALL zeroing(brsum)!,SIZE(brsum))
    CALL zeroing(brsum3)!,SIZE(brsum3))
    DO i=1,nstate
       DO j=1,i
          DO ikdp=1,nkdp
             brsum(i,j)=brsum(i,j)+0.5_real_8*wkdp(ikdp)*(bkdp(i,j,ikdp)+bkdp(&
                  j,i,ikdp))
          ENDDO
          brsum(j,i)=brsum(i,j)
          DO ic=1,3
             DO ikdp=1,nkdp
                brsum3(i,j,ic)=brsum3(i,j,ic)+wkdp(ikdp)*xkdp(ic,ikdp)*(&
                     bkdp(i,j,ikdp)-bkdp(j,i,ikdp))
             ENDDO
             brsum3(j,i,ic)=-brsum3(i,j,ic)
          ENDDO
       ENDDO
    ENDDO
    ! 
    ! Calculation of the kinetic energy contribution to the stress (k.p)
    ! 
    ! ==   1) g.g + k.k term

    DO kk=1,6
       ! 
       ! Off-diagonal contributions (different states), only g.g
       DO i=1,nstate
          DO j=1,i-1
             cgcg=0._real_8
             DO ig=1,ncpw%ngw
                cgcg=cgcg+REAL(CONJG(c(ig,i))*c(ig,j))*gagk(ig,kk)
             ENDDO
             cgcg=4._real_8*cgcg
             dekin(kk)=dekin(kk)+REAL(brsum(i,j)*cgcg)
          ENDDO
       ENDDO
       ! 
       ! Diagonal contributions (same state)
       DO i=1,nstate
          ! g.g
          cgcg=0._real_8
          DO ig=1,ncpw%ngw
             cgcg = cgcg + REAL(CONJG(c(ig,i))*c(ig,i))*gagk(ig,kk)
          ENDDO
          cgcg=2._real_8*cgcg
          dekin(kk)=dekin(kk)+REAL(brsum(i,i)*cgcg)
          ! k.k
          xkkk=0._real_8
          DO ikdp=1,nkdp
             xkkk=xkkk+wkdp(ikdp)*REAL(bkdp(i,i,ikdp))*xkk(kk,ikdp)
          ENDDO
          xkkk=xkkk/REAL(parai%nproc,kind=real_8)
          dekin(kk)=dekin(kk)+xkkk
          ! 
       ENDDO
       ! 
    ENDDO
    ! 
    ! ==   2) g.k + k.g term
    ! 
    ! Only off-diagonal contributions (different states).
    DO kk=1,6
       aux=0._real_8
       DO i=1,nstate
          DO j=1,i-1
             aux=aux-AIMAG(brsum3(i,j,alpha(kk)))*pkdp(beta(kk),i,j)-&
                  AIMAG(brsum3(i,j,beta(kk)))*pkdp(alpha(kk),i,j)
          ENDDO
       ENDDO
       aux=aux/REAL(parai%nproc,kind=real_8)
       dekin(kk)=dekin(kk)+aux
    ENDDO
    ! 
    DO kk=1,6
       dekin(kk) = - dekin(kk)
    ENDDO
    ! 
    RETURN
  END SUBROUTINE kdp_stress_kin

END MODULE kdp_stress_kin_utils
