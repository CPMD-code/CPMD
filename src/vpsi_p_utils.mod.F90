MODULE vpsi_p_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: hg,&
                                             indzs,&
                                             nzhs
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE parac,                           ONLY: parai
  USE prcp,                            ONLY: prcp_com
  USE reshaper,                        ONLY: reshape_inplace
  USE response_pmod,                   ONLY: c0real,&
                                             response1
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: v0psi1_p
  PUBLIC :: v1psi0_p
  PUBLIC :: allocate_c0real_p
  PUBLIC :: release_c0real_p

CONTAINS

  ! ==================================================================
  SUBROUTINE v0psi1_p(c0,c1,c12,vpot,psi,drhoe,nstate)
    ! ==================================================================
    ! == applies the ground state potential (V0) to the perturbation  ==
    ! == orbitals. TUNED VERSION, under construction!                 ==
    ! == ************************************************************ ==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == C0,C1:  IN ORBITALS                                          ==
    ! == DRHOE:  OUT FIRST ORDER DENSITY in REAL SPACE                ==
    ! == C12:    V[ n0 ] applied to C1                                ==
    ! ==--  d.s. 3/2001  ---------------------------------------------==
    ! ==--  d.s.11/2002  (lsd) ---------------------------------------==
    REAL(real_8)                             :: vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8), TARGET                  :: psi(:)
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c12(ncpw%ngw,nstate), &
                                                c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'v0psi1_p'

    COMPLEX(real_8)                          :: fm, fp, psii, psin
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: aux(:)
    INTEGER                                  :: ierr, ig, ilsd1, ilsd2, ir, &
                                                is, is1, is2, iss, isub
    REAL(real_8)                             :: density_factor, ff1, ff2, g2, &
                                                rho1, rho2
    REAL(real_8), POINTER                    :: auxre(:,:), psire(:,:)

    CALL tiset('    V0PSI1',isub)
    IF (group%nogrp.NE.1) THEN
       CALL stopgm('V0psi1','ORBITAL-GROUPS not implemented.',& 
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(aux(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(aux, (/2, maxfft/), auxre)
    CALL reshape_inplace(psi, (/2, maxfft/), psire)

    IF (response1%trho1) CALL zeroing(drhoe)!,nnr1*clsd%nlsd)
    IF (tkpts%tkpnt) THEN
       CALL stopgm('V0psi1','K-POINTS not implemented.',& 
            __LINE__,__FILE__)
    ENDIF

    DO is=1,nstate,2
       is1 = is
       is2 = is+1
       iss = is / 2 + 1
       ilsd1=1
       ilsd2=1
       IF (cntl%tlsd) THEN
          IF (is1 .GT. spin_mod%nsup) ilsd1 = 2
          IF (is2 .GT. spin_mod%nsup) ilsd2 = 2
       ENDIF
       ff1 = crge%f(is1,1)
       ff2 = 0._real_8
       IF (is2 .LE. nstate) ff2=crge%f(is2,1)
       ! ==------------------------------------------------------------==
       ! == FIRST, the GROUND STATE ORBITALS are FFTed to real space   ==
       ! ==------------------------------------------------------------==
       IF (response1%trho1) THEN
          IF (response1%tkeeprealspacewfn) THEN
             CALL dcopy(2*fpar%nnr1,c0real(1,iss),1,aux,1)
          ELSE
             CALL zeroing(aux)!,maxfft)
             IF (is2.GT.nstate) THEN
                !CDIR NODEP
                !ocl NOVREC
#ifdef __SR8000
                !poption parallel
#endif
                !$omp parallel do private(IG) shared(IS1)
                DO ig=1,ncpw%ngw
                   aux(nzhs(ig))      = c0(ig,is1)
                   aux(indzs(ig))     = CONJG(c0(ig,is1))
                ENDDO
                IF (geq0) aux(nzhs(1))=c0(1,is1)
             ELSE
                !CDIR NODEP
                !ocl NOALIAS 
                !ocl NOVREC  
#ifdef __SR8000
                !poption parallel
#endif
                !$omp parallel do private(IG) shared(IS1,IS2)
                DO ig=1,ncpw%ngw
                   auxre(1,nzhs(ig))  = REAL(c0(ig,is1))-AIMAG(c0(ig,is2))
                   auxre(2,nzhs(ig))  = AIMAG(c0(ig,is1))+REAL(c0(ig,is2))
                   auxre(1,indzs(ig)) = REAL(c0(ig,is1))+AIMAG(c0(ig,is2))
                   auxre(2,indzs(ig)) =-AIMAG(c0(ig,is1))+REAL(c0(ig,is2))
                ENDDO
                IF (geq0) aux(nzhs(1))= c0(1,is1) + uimag*c0(1,is2)
             ENDIF         ! odd number of states and IS=last state...
             CALL  invfftn(aux,.TRUE.,parai%allgrp)! DO IT.
          ENDIF             ! keep real space wfns (C0)
       ENDIF                 ! if (trho1)
       ! AUX now contains the real-space ground state wavefunctions.        
       ! ==------------------------------------------------------------==
       ! == SECOND, the PERTURBATION ORBITALS are FFTed to real space  ==
       ! ==------------------------------------------------------------==
       CALL zeroing(psi)!,maxfft)
       IF (is2.GT.nstate) THEN
          !CDIR NODEP
          !ocl NOVREC
#ifdef __SR8000
          !poption parallel
#endif
          !$omp parallel do private(IG) shared(IS1)
          DO ig=1,ncpw%ngw
             psi(nzhs(ig))=c1(ig,is1)
             psi(indzs(ig))=CONJG(c1(ig,is1))
          ENDDO
          IF (geq0) psi(nzhs(1))=c1(1,is1)
       ELSE
          !CDIR NODEP
          !ocl NOALIAS 
          !ocl NOVREC  
#ifdef __SR8000
          !poption parallel
#endif
          !$omp parallel do private(IG) shared(IS1,IS2)
          DO ig=1,ncpw%ngw
             psire(1,nzhs(ig))  = REAL(c1(ig,is1)) - AIMAG(c1(ig,is2))
             psire(2,nzhs(ig))  = AIMAG(c1(ig,is1)) + REAL(c1(ig,is2))
             psire(1,indzs(ig)) = REAL(c1(ig,is1)) + AIMAG(c1(ig,is2))
             psire(2,indzs(ig)) =-AIMAG(c1(ig,is1)) + REAL(c1(ig,is2))
          ENDDO
          IF (geq0) psi(nzhs(1))=c1(1,is1)+uimag*c1(1,is2)
       ENDIF                 ! odd number of states and IS=last state...
       CALL  invfftn(psi,.TRUE.,parai%allgrp)    ! DO IT.
       ! NB: If IS2 .GT. NSTATE, then this FFT will result in an array
       ! PSI which is purely real. Therefore all the calculations below
       ! will be done with zero numbers and thus not contribute.


       ! ==------------------------------------------------------------==
       ! == Compute the PERTURBATION DENSITY n^1 (r), and              ==
       ! == Apply the potential (V) to the unperturbed orbitals        ==
       ! ==------------------------------------------------------------==
       IF (response1%trho1) THEN
          !$omp parallel do private(ir,rho1,rho2)
          DO ir=1,fpar%nnr1
             rho1 = auxre(1,ir)*psire(1,ir)
             rho2 = auxre(2,ir)*psire(2,ir)
             drhoe(ir,ilsd1)   = drhoe(ir,ilsd1)   + ff1*rho1
             drhoe(ir,ilsd2)   = drhoe(ir,ilsd2)   + ff2*rho2
             psire(1,ir)   = vpot(ir,ilsd1)  *psire(1,ir)
             psire(2,ir)   = vpot(ir,ilsd2)  *psire(2,ir)
          ENDDO
       ELSE
          !$omp parallel do private(ir)
          DO ir=1,fpar%nnr1
             psire(1,ir)   = vpot(ir,ilsd1)  *psire(1,ir)
             psire(2,ir)   = vpot(ir,ilsd2)  *psire(2,ir)
          ENDDO
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions and add the        ==
       ! == kinetic energy term.                                       ==
       ! ==------------------------------------------------------------==
       IF (prcp_com%akin.GT.1.e-10_real_8) THEN
          CALL stopgm('V0psi1','CONSTANT CUTOFF not implemented.',& 
               __LINE__,__FILE__)
       ENDIF
       !$omp parallel do private(IG,G2,psin,psii,FP,FM)
       DO ig=1,ncpw%ngw
          g2 = - 0.5_real_8*parm%tpiba2*hg(ig)
          psin=psi(nzhs(ig))! access only once
          psii=psi(indzs(ig))! these mem locations
          fp=-0.5_real_8*(psin+psii)
          fm=-0.5_real_8*(psin-psii)
          ! mb           FP = - 0.5_real_8*(PSI(NZHS(IG))+PSI(INDZS(IG)))
          ! mb           FM = - 0.5_real_8*(PSI(NZHS(IG))-PSI(INDZS(IG)))
          c12(ig,is1)= ff1 * (g2*c1(ig,is1)&
               + CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
          IF (is2.LE.nstate)&
               c12(ig,is2)= ff2 * (g2*c1(ig,is2)&
               +  CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       ENDDO
    ENDDO                     ! Loop over IS

    IF (response1%trho1) THEN
       IF (cntl%tlsd) CALL daxpy(fpar%nnr1,1._real_8,drhoe(1,2),1,drhoe(1,1),1)
       density_factor = 2._real_8 / parm%omega
       CALL dscal(clsd%nlsd*fpar%nnr1,density_factor,drhoe,1)
    ENDIF

    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    V0PSI1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE v0psi1_p
  ! ==================================================================
  SUBROUTINE v1psi0_p(c0,c12,vpot,psi,nstate)
    ! ==================================================================
    ! == VPOT:   IN INPUT POTENTIAL V_Hxc [n1]                        ==
    ! == C12:    IN partial gradient  dE/dC1                         ==
    ! == C12:    OUT C12 + V_Hxc [n1] * C0                            ==
    ! ==--  d.s. 3/2001  ---------------------------------------------==
    ! ==--  d.s.11/2002  (lsd)  --------------------------------------==
    REAL(real_8)                             :: vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8), TARGET                  :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c12(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    COMPLEX(real_8)                          :: fm, fp, psii, psin
    INTEGER                                  :: ig, ilsd1, ilsd2, ir, is, &
                                                is1, is2, iss, isub
    REAL(real_8)                             :: ff1, ff2
    REAL(real_8), POINTER                    :: psire(:,:)

    CALL tiset('    V0PSI1',isub)

    IF (group%nogrp.NE.1) THEN
       CALL stopgm('V0psi1','ORBITAL-GROUPS not implemented.',& 
            __LINE__,__FILE__)
    ENDIF
    IF (tkpts%tkpnt) THEN
       CALL stopgm('V0psi1','K-POINTS not implemented.',& 
            __LINE__,__FILE__)
    ENDIF

    CALL reshape_inplace(psi, (/2, maxfft/), psire)

    DO is=1,nstate,2
       is1 = is
       is2 = is+1
       iss = is / 2 + 1
       ilsd1=1
       ilsd2=1
       IF (cntl%tlsd) THEN
          IF (is1 .GT. spin_mod%nsup) ilsd1 = 2
          IF (is2 .GT. spin_mod%nsup) ilsd2 = 2
       ENDIF
       ff1 = crge%f(is1,1)
       ff2 = 0._real_8
       IF (is2 .LE. nstate) ff2=crge%f(is2,1)
       ! ==------------------------------------------------------------==
       ! == FIRST, the GROUND STATE ORBITALS are FFTed to real space   ==
       ! ==------------------------------------------------------------==
       IF (response1%tkeeprealspacewfn) THEN
          CALL dcopy(2*fpar%nnr1,c0real(1,iss),1,psi,1)
       ELSE
          CALL zeroing(psi)!,maxfft)
          IF (is2.GT.nstate) THEN
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IG) shared(IS1)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig))=c0(ig,is1)
                psi(indzs(ig))=CONJG(c0(ig,is1))
             ENDDO
             IF (geq0) psi(nzhs(1))=c0(1,is1)
          ELSE
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IG) shared(IS1,IS2)
             DO ig=1,ncpw%ngw
                psire(1,nzhs(ig))  = REAL(c0(ig,is1)) - AIMAG(c0(ig,is2))
                psire(2,nzhs(ig))  = AIMAG(c0(ig,is1)) + REAL(c0(ig,is2))
                psire(1,indzs(ig)) = REAL(c0(ig,is1)) + AIMAG(c0(ig,is2))
                psire(2,indzs(ig)) =-AIMAG(c0(ig,is1)) + REAL(c0(ig,is2))
             ENDDO
             IF (geq0) psi(nzhs(1))=c0(1,is1)+uimag*c0(1,is2)
          ENDIF             ! odd number of states and IS=last state...
          CALL  invfftn(psi,.TRUE.,parai%allgrp)! DO IT.
       ENDIF                 ! REAL SPACE c0
       ! ==------------------------------------------------------------==
       ! == Apply the potential (VPOT) to the orbitals                 ==
       ! ==------------------------------------------------------------==
       !$omp parallel do private(ir) shared(ilsd1,ilsd2)
       DO ir=1,fpar%nnr1
          psire(1,ir)   = vpot(ir,ilsd1)  *psire(1,ir)
          psire(2,ir)   = vpot(ir,ilsd2)  *psire(2,ir)
       ENDDO

       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions                    ==
       ! ==------------------------------------------------------------==
       IF (prcp_com%akin.GT.1.e-10_real_8) THEN
          CALL stopgm('V0psi1','CONSTANT CUTOFF not implemented.',& 
               __LINE__,__FILE__)
       ELSE
          !CDIR NODEP
          DO ig=1,ncpw%ngw
             psii=psi(nzhs(ig))! access only once
             psin=psi(indzs(ig))! these mem locations
             fp=0.5_real_8*(psii+psin)
             fm=0.5_real_8*(psii-psin)
             ! mb            FP = 0.5_real_8*(PSI(NZHS(IG))+PSI(INDZS(IG)))
             ! mb            FM = 0.5_real_8*(PSI(NZHS(IG))-PSI(INDZS(IG)))
             c12(ig,is1)=c12(ig,is1)-ff1*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             IF (is2.LE.nstate)&
                  c12(ig,is2)=c12(ig,is2)-ff2*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ENDIF
    ENDDO                     ! Loop over IS

    CALL tihalt('    V0PSI1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE v1psi0_p
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE allocate_c0real_p(c0,psi,nstate)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'allocate_c0real_p'

    INTEGER                                  :: ierr, ig, il_c0real, is, is1, &
                                                is2, iss

! ==--------------------------------------------------------------==

    il_c0real = INT((nstate)/2+1)*2 * fpar%nnr1
    ALLOCATE(c0real(fpar%nnr1,il_c0real/fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0real)!, nnr1*nstate/2)

    DO is=1,nstate,2
       CALL zeroing(psi)!,maxfft)
       is1 = is
       is2 = is+1
       IF (is2.GT.nstate) THEN
          !CDIR NODEP
#ifdef __SR8000
          !poption parallel
#endif
          !$omp parallel do private(IG) shared(IS1)
          DO ig=1,ncpw%ngw
             psi(nzhs(ig))   = c0(ig,is1)
             !psi(1,nzhs(ig))   = real(c0(ig,is1))
             !psi(2,nzhs(ig))   = aimag(c0(ig,is1))

             psi(indzs(ig))  = CONJG(c0(ig,is1))

             !psi(1,indzs(ig))  = real(c0(ig,is1))
             !psi(2,indzs(ig))  =-aimag(c0(ig,is1))
          ENDDO
          IF (geq0) THEN
             psi(nzhs(1)) = CMPLX( REAL( c0(1,is1), KIND=real_8 ), 0.0_real_8, KIND=real_8 )
             !psi(1,nzhs(1))=real(c0(1,is1))
             !psi(2,nzhs(1))=0._real_8
          ENDIF
       ELSE
          !CDIR NODEP
#ifdef __SR8000
          !poption parallel
#endif
          !$omp parallel do private(IG) shared(IS1,IS2)
          DO ig=1,ncpw%ngw

             psi(nzhs(ig)) = CMPLX(REAL(c0(ig,is1),KIND=real_8)-AIMAG(c0(ig,is2)), &
                  AIMAG(c0(ig,is1))+REAL(c0(ig,is2),KIND=real_8), KIND=real_8 )

             !psi(1,nzhs(ig))  = real(c0(ig,is1)) - aimag(c0(ig,is2))
             !psi(2,nzhs(ig))  = aimag(c0(ig,is1)) + real(c0(ig,is2))

             psi(indzs(ig)) = CMPLX(REAL(c0(ig,is1),KIND=real_8)+AIMAG(c0(ig,is2)), &
                  -AIMAG(c0(ig,is1))+REAL(c0(ig,is2),KIND=real_8), KIND=real_8 )

             !psi(1,indzs(ig)) = real(c0(ig,is1)) + aimag(c0(ig,is2))
             !psi(2,indzs(ig)) =-aimag(c0(ig,is1)) + real(c0(ig,is2))
          ENDDO
          IF (geq0) THEN
             psi(nzhs(1)) = REAL(c0(1,is1),KIND=real_8)
             !psi(1,nzhs(1)) = real(c0(1,is1))
             !psi(2,nzhs(1)) = real(c0(1,is2))
          ENDIF
       ENDIF                 ! odd number of states and IS=last state...
       CALL  invfftn(psi,.TRUE.,parai%allgrp)    ! DO IT.
       iss = is / 2 + 1
       CALL dcopy(2*fpar%nnr1,psi,1,c0real(1,iss),1)
    ENDDO                     ! IS
    RETURN
  END SUBROUTINE allocate_c0real_p
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE release_c0real_p
    CHARACTER(*), PARAMETER :: procedureN = 'release_c0real_p'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==

    DEALLOCATE(c0real,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE release_c0real_p
  ! ==================================================================

END MODULE vpsi_p_utils
