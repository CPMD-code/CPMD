MODULE strs
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == Stress quantities                                            ==
  ! ==--------------------------------------------------------------==
  ! == DEHC  Calculated but not used in HTRSTR                      ==
  ! == DEHT  Hartree energy contribution to the stress tensor       ==
  ! == DEPST Calculated but not used in DEPST                       ==
  ! == DEKIN Kinetic energy contribution (to be symmetrized)        ==
  ! == DEXC  Exchange-Correlation contribution                      ==
  ! == DESR  Residual part of the Ion-Ion interaction contribution  ==
  ! == DEPS  Calculated in VLOCST                                   ==
  ! == DENL  Non-local pseudopotential contribution                 ==
  ! == DLAM  Orthogonality constrained in Vanderbilt pp (NLSL)      ==
  ! == DEGC  Gradient Correction                                    ==
  ! == DECC  Non-linear core correction                             ==
  ! == DELTA ???                                                    ==
  ! == PAIU  Stress tensor in cartesian coordinates                 ==
  ! == PAIL  Stress tensor in crystal coordinates                   ==
  ! ==================================================================
  COMPLEX(real_8), SAVE :: dehc(6)=(0.0_real_8,0.0_real_8),&
       depst(6)=(0.0_real_8,0.0_real_8)
  REAL(real_8), SAVE :: dekin(6)=0.0_real_8,&
       deht(6)=0.0_real_8,&
       dexc(6)=0.0_real_8,&
       desr(6)=0.0_real_8,&
       deps(6)=0.0_real_8,&
       denl(6)=0.0_real_8,&
       dlam(6)=0.0_real_8,&
       degc(6)=0.0_real_8,&
       decc(6)=0.0_real_8,&
       delta(3,3)=0.0_real_8,&
       paiu(3,3)=0.0_real_8,&
       pail(3,3)=0.0_real_8
  ! ==--------------------------------------------------------------==
  ! == Indexes fro Voigt notation                                   ==
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE :: alpha(6)= 0,beta(6)= 0
  TYPE :: presi_t
     INTEGER :: alpha(6) = 0
     INTEGER :: beta(6) = 0
  END TYPE presi_t
  !  TYPE(presi_t), SAVE :: presi
  ! ==================================================================

END MODULE strs
