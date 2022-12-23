MODULE conduct_utils
  USE cnst,                            ONLY: pi,&
                                             ry
  USE condu,                           ONLY: condavg,&
                                             condavg2,&
                                             condpa,&
                                             conduct,&
                                             conduct2,&
                                             normcon
  USE cppt,                            ONLY: gk
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conductivity

CONTAINS

  ! ==================================================================
  SUBROUTINE conductivity(c0,we,f,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == Compute optical conductivity using Greenwood-Kubo formula    ==
    ! == Sigma(w) = 2 pi e^2/(3 m^2 Vw) sum_i,j (f_i-fj) |<j|p|i>|^2  ==
    ! == * delta(Ej-Ei-h_bar w)                                       ==
    ! ==                                                              ==
    ! == Original version - P. L. Silvestrelli (18/10/1995)           ==
    ! == Updated version - M. Boero (22/1/2002)                       ==
    ! ==                   M. Boero & P. L. Silvestrelli (10/ 4/2002) ==
    ! ==                   M. Boero & F. Gervasio        ( 8/ 5/2002) ==
    ! ==                   M. Boero                      (02/06/2005) ==
    ! ==     (printout <i|x|j> ~ (f_i-f_j)*|M(i,j)|^2/(eps_j - eps_i) ==
    ! ==                                                              ==
    ! == N. B.:  If TSEGHE is .true. additional quantities are printed==
    ! == on the FORT.89 file, in the following form:                  ==
    ! ==                                                              ==
    ! == i              f_i    eps_i                                  ==
    ! == i    j=i+1     (f_i - f_j)   (eps_j - eps_i)   M_ij          ==
    ! == i    j=i+2     (f_i - f_j)   (eps_j - eps_i)   M_ij          ==
    ! == ...................................................          ==
    ! == ...................................................          ==
    ! == ...................................................          ==
    ! == i    j=N       (f_i - f_j)   (eps_j - eps_i)   M_ij          ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                f(nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'conductivity'
    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8) 
    REAL(real_8), PARAMETER                  :: ftoll = 1.e-10_real_8 

    COMPLEX(real_8)                          :: caux, csumx, csumy, csumz
    INTEGER                                  :: i, ibin, ierr, ig, ik, j
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror
    REAL(real_8) :: acon, aiaux, aux, auxx, auxy, auxz, csri(2), deng, df, &
      factor, fracon, ryby2, sum, sumx, sumx2, sumy, sumy2, sumz, sumz2

! ==--------------------------------------------------------------==
! == FACTOR = (2._real_8*PI/3._real_8/OMEGA)                                ==
! == FACTOR is divided by 2 because this factor is already        ==
! == contained in f_i: 4.6E+04 is (?) the conversion factor       ==
! == from a.u. to (ohm cm)**(-1)                                  ==
! ==--------------------------------------------------------------==

    ryby2=2._real_8*ry
    factor=pi/3._real_8/parm%omega*ryby2**2*4.6e04_real_8
    acon=parm%tpiba2
    IF (ifirst.EQ.0) THEN
       ifirst=1

       ! [irina] TODO free is missing for CONDUCT, CONDUCT2, ...

       ALLOCATE(conduct(condpa%nconduct),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(conduct2(condpa%nconduct),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(condavg(condpa%nconduct),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(condavg2(condpa%nconduct),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(normcon(condpa%nconduct),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(condavg)!,condpa%nconduct)
       CALL zeroing(condavg2)!,condpa%nconduct)
       CALL zeroing(normcon)!,condpa%nconduct)
    ENDIF
    CALL zeroing(conduct)!,condpa%nconduct)
    ! mb-> write the transition matrix elements into the file MATRIX.DAT
    ! mb-> Note: r=(x,y,z); only NON-ZERO transitions are written
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(59,'MATRIX.DAT',fo_def,ferror)
       IF (paral%io_parent)&
            WRITE(59,'(A,A,A)')&
            "  |k>  |i>  |j>  f(i) f(j)",&
            "  E(i) [eV]  E(j) [eV]     |<i|x|j>|**2",&
            "     |<i|y|j>|**2     |<i|z|j>|**2     |<i|r|j>|**2"
       IF (paral%io_parent)&
            WRITE(59,*)
    ENDIF
    ! mb-> uncomment the following section to check and print out
    ! mb-> eigenvalues and occupation numbers in debugging
    ! IF(PARENT) THEN
    ! DO IK=1,NKPOINT
    ! WRITE(6,'(1x,''K-POINT ='',i6)')IK
    ! WRITE(6,'(1x,'' STATE  OCCUPATION  EIGENVALUE (eV)'')')
    ! DO I=1,NSTATE
    ! WRITE(6,'(1x,i5,2x,f7.3,6x,f10.3)')
    ! .              I,F(I,IK),RYBY2*WE(I,IK)
    ! ENDDO
    ! WRITE(6,*)
    ! ENDDO
    ! ENDIF
    ! mb-> end of eigenvalues and occupation section
    DO ik=1,nkpoint
       DO i=1,nstate-1
          DO j=i+1,nstate
             df=f(i,ik)-f(j,ik)
             IF (df.LT.ftoll) GOTO 999
             deng=ryby2*(we(j,ik)-we(i,ik))
             ibin=INT(deng/condpa%condstep)+1
             IF (ibin.LE.condpa%nconduct) THEN
                IF (tkpts%tkpnt) THEN
                   csumx=zzero
                   csumy=zzero
                   csumz=zzero
                   ! mb-> first G > 0 with the k-points
#ifdef __SR8000
                   !poption parallel
#endif
                   !$omp parallel do private(IG,CAUX) reduction(+:CSUMX,CSUMY,CSUMZ)
                   DO ig=1,ncpw%ngw
                      caux=CONJG(c0(ig,i,ik))*c0(ig,j,ik)
                      csumx=csumx+(rk(1,ik)+gk(1,ig))*caux
                      csumy=csumy+(rk(2,ik)+gk(2,ig))*caux
                      csumz=csumz+(rk(3,ik)+gk(3,ig))*caux
                   ENDDO
                   ! mb-> then G < 0 with the k-points
#ifdef __SR8000
                   !poption parallel
#endif
                   !$omp parallel do private(IG,CAUX) reduction(+:CSUMX,CSUMY,CSUMZ)
                   DO ig=1,ncpw%ngw
                      caux=CONJG(c0(ig+ncpw%ngw,i,ik))*c0(ig+ncpw%ngw,j,ik)
                      csumx=csumx+(rk(1,ik)-gk(1,ig))*caux
                      csumy=csumy+(rk(2,ik)-gk(2,ig))*caux
                      csumz=csumz+(rk(3,ik)-gk(3,ig))*caux
                   ENDDO
                   csri(1)=REAL(csumx)
                   csri(2)=AIMAG(csumx)
                   CALL mp_sum(csri,2,parai%allgrp)
                   csumx=CMPLX(csri(1),csri(2),kind=real_8)
                   csri(1)=REAL(csumy)
                   csri(2)=AIMAG(csumy)
                   CALL mp_sum(csri,2,parai%allgrp)
                   csumy=CMPLX(csri(1),csri(2),kind=real_8)
                   csri(1)=REAL(csumz)
                   csri(2)=AIMAG(csumz)
                   CALL mp_sum(csri,2,parai%allgrp)
                   csumz=CMPLX(csri(1),csri(2),kind=real_8)
                   sumx2=REAL(csumx*CONJG(csumx))
                   sumy2=REAL(csumy*CONJG(csumy))
                   sumz2=REAL(csumz*CONJG(csumz))
                   sum=sumx2+sumy2+sumz2
                ELSE
                   sumx=0._real_8
                   sumy=0._real_8
                   sumz=0._real_8
                   ! mb->  Sum_{G} <c(j)|G|c(i)> = 2 * Sum_{G>0} G Im{<c(j)|c(i)>} at Gamma
#ifdef __SR8000
                   !poption parallel
#endif
                   !$omp parallel do private(IG,AIAUX) reduction(+:SUMX,SUMY,SUMZ)
                   DO ig=1,ncpw%ngw
                      aiaux=AIMAG(CONJG(c0(ig,i,ik))*c0(ig,j,ik))
                      sumx=sumx+gk(1,ig)*aiaux
                      sumy=sumy+gk(2,ig)*aiaux
                      sumz=sumz+gk(3,ig)*aiaux
                   ENDDO
                   CALL mp_sum(sumx,parai%allgrp)
                   CALL mp_sum(sumy,parai%allgrp)
                   CALL mp_sum(sumz,parai%allgrp)
                   sumx2=4.0_real_8*sumx*sumx
                   sumy2=4.0_real_8*sumy*sumy
                   sumz2=4.0_real_8*sumz*sumz
                   sum=sumx2+sumy2+sumz2
                ENDIF
                ! AK 2005/03/26: not used
                ! FREQU=(real(IBIN-1,kind=real_8)+0.5_real_8)*CONDSTEP
                fracon=factor*df*acon/deng/condpa%condstep
                auxx=fracon*sumx2
                auxy=fracon*sumy2
                auxz=fracon*sumz2
                aux=fracon*sum
                IF (tkpts%tkpnt) THEN
                   auxx=auxx*wk(ik)
                   auxy=auxy*wk(ik)
                   auxz=auxz*wk(ik)
                   aux=aux*wk(ik)
                ENDIF
                conduct(ibin)=conduct(ibin)+aux
                normcon(ibin)=normcon(ibin)+1
                IF (paral%io_parent)&
                     WRITE(59,'(3i5,1x,2f5.2,2f11.3,4e17.8)')&
                     ik,i,j,f(i,ik),f(j,ik)&
                     ,ryby2*we(i,ik),ryby2*we(j,ik)&
                     ,auxx,auxy,auxz,aux
             ENDIF
999          CONTINUE
          ENDDO
       ENDDO
    ENDDO
    ! mb-> and now fill the remaining bins
    DO ibin=1,condpa%nconduct
       conduct2(ibin)=conduct(ibin)*conduct(ibin)
       condavg(ibin) =condavg(ibin)+conduct(ibin)
       condavg2(ibin)=condavg2(ibin)+conduct2(ibin)
    ENDDO
    ! ==--------------------------------------------------------------==
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(59)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE conductivity
  ! ==================================================================

END MODULE conduct_utils
