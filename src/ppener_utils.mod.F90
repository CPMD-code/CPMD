#include "cpmd_global.h"

MODULE ppener_utils
  USE cppt,                            ONLY: nzh,&
                                             scg
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE nvtx_utils
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ppener

CONTAINS

  ! ==================================================================
  SUBROUTINE ppener(eh,ei,ee,eps,v,vtemp,eivps,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == ELECTROSTATIC AND PSEUDOPOTENTIAL ENERGIES AND THE           ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    ! == Input:                                                       ==
    ! ==   V     Electronic density in G-space (input)                ==
    ! ==   EIVPS phase factor times local pseudopotential  (VPS)      ==
    ! ==   EIROP phase factor times Gaussian charge distributions     ==
    ! ==         which replaced ionic point charges (RHOPS)           ==
    ! == Output:                                                      ==
    ! ==   EH    Hartree term                                         ==
    ! ==         (e-e, ion-ion and e-ion [part] interactions)         ==
    ! ==         Used by Car-Parrinello scheme                        ==
    ! ==   EI    Gaussian charge dist. (RHOPS) (Ion-Ion interaction)  ==
    ! ==   EE    e-e interaction                                      ==
    ! ==         Both used by diagonalisation schemes (no G=0 term)   ==
    ! ==   EPS   Integration of electronic density*local pp           ==
    ! ==         EPSEU=2._real_8*real(EPS)*OMEGA)                         ==
    ! ==   VPLOC G=0 pp part (given from pp data) (Energy=NEL*VPLOC)  ==
    ! ==   VTEMP Stored the potential in G-space                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eh, ei, ee, eps, v(:), &
                                                vtemp(:), eivps(:), eirop(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ppener'

    COMPLEX(real_8)                          :: rhet, rhets, rhog, rhogs, rp, &
                                                vcg, vp
    INTEGER                                  :: ig, ig1, isub

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    IF (geq0) THEN
       vp=eivps(1)
       vploc=REAL(vp)
       eps=0.5_real_8*vp*CONJG(v(nzh(1)))
       ig1=2
       rp=eirop(1)
       rhet=v(nzh(1))
       rhog=rhet+rp
       eh=0.5_real_8*scg(1)*REAL(rhog)*REAL(rhog)
       ei=0.5_real_8*scg(1)*rp*rp
       ee=0.5_real_8*scg(1)*rhet*rhet
       vtemp(1)=scg(1)*rhog
    ELSE
       ig1=1
       eps=(0.0_real_8,0.0_real_8)
       vploc=0.0_real_8
       eh=(0.0_real_8,0.0_real_8)
       ei=(0.0_real_8,0.0_real_8)
       ee=(0.0_real_8,0.0_real_8)
    ENDIF
    !$omp  parallel do private(IG,VP,RP,RHET,RHOG,RHETS,RHOGS,VCG) &
    !$omp  reduction(+:EH,EI,EE,EPS)
#ifdef __SR8000
    !poption parallel
    !poption tlocal(IG,VP,RP,RHET,RHOG,RHETS,RHOGS,VCG)
    !poption psum(EH,EI,EE,EPS)
#endif 
    DO ig=ig1,ncpw%nhg
       ! ==------------------------------------------------------------==
       ! == In VTEMP it is stored temporarily the potential in G-space.==
       ! ==------------------------------------------------------------==
       vp=eivps(ig)
       rp=eirop(ig)
       rhet=v(nzh(ig))
       rhog=rhet+rp
       rhets=CONJG(rhet)
       rhogs=CONJG(rhog)
       vcg=scg(ig)*rhog
       vtemp(ig)=vcg+vp
       ! Hartree term (e-e, ion-ion and e-ion [part] interactions)
       eh=eh+vcg*rhogs
       ! Gaussian charge dist. (RHOPS) Ion-Ion interaction
       ei=ei+scg(ig)*rp*CONJG(rp)
       ! e-e interaction
       ee=ee+scg(ig)*rhet*CONJG(rhet)
       eps=eps+rhets*vp
    ENDDO
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ppener
  ! ==================================================================

END MODULE ppener_utils
