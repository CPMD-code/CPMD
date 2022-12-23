MODULE sphe
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE TO HAVE SPHERICAL CUTOFF FOR K-POINTS I.E.      ==
  ! == | G + K |^2 < Ecutoff (GCUTW)                                ==
  ! ==--------------------------------------------------------------==
  ! == TSPHERE = .TRUE. Use spherical cutoff (default)              ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: tsphere
  ! ==--------------------------------------------------------------==
  ! == GCUTKA Maximum K point amplitude (Kmax) in Brillouin zone.   ==
  ! == For each k and each G, |k+G|^2 < GCUTW                       ==
  ! == so we define GCUTWMIN and GCUTWMAX as:                       ==
  ! ==    GCUTWMIN < |G|^2 < GCUTWMAX, we have to apply a mask.     ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: gcutka,gcutwmin,gcutwmax
  ! ==--------------------------------------------------------------==
  ! == MASKL Number of G vectors                                    ==
  ! ==       with GCUTW - GCUTKAMIN <= |G|^2 < GCUTW + GCUTKAMAX    ==
  ! == NGWKPS(1:NKPTS) Number of G-vectors for each K-points        ==
  ! == NGW_GAMMA       Number of G-vectors for Gamma point          ==
  ! == MASKGW(2,MASKL,NKPNT) Mask                                   ==
  ! ==                (1,:,:) for  G                                ==
  ! ==                (2,:,:) for -G                                ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: ngwkps(:)
  INTEGER :: maskl,ngw_gamma
  REAL(real_8), ALLOCATABLE :: maskgw(:,:,:)

  ! ==================================================================

END MODULE sphe
