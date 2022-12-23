MODULE symm4
  USE kinds,                           ONLY: int_8

  IMPLICIT NONE

  ! ==================================================================
  ! == Use in parallel version                                      ==
  ! == NUMSHEL   Number of shells of equivalent G-vectors           ==
  ! ==           equivalent means symmetry links all G-vectors      ==
  ! == INDSHEL(NHG) For each G-vector index of the equivalent       ==
  ! ==          shell, i.e. 2 G-vectors related by symmetry has the ==
  ! ==          same index                                          ==
  ! ==          If TSYMMORPHIC the sign of INDSHEL indicates that   ==
  ! ==          -IG is taken instead of IG                          ==
  ! == ISYSHEL(2,NHG)                                               ==
  ! ==          If TSYMMORPHIC                                      ==
  ! ==          give for IG  if ISYSHEL>0                           ==
  ! ==          or   for -IG if ISYSHEL<0 the multiplicative factor ==
  ! == NINSHEL(NHG) !!!WARNING INTEGER*8                            ==
  ! ==          If TSYMMORPHIC unused                               ==
  ! ==          If .NOT.TSYMMORPHIC                                 ==
  ! ==          For equivalent shell,                               ==
  ! ==          give for each G-vector, the related rotations.      ==
  ! ==          For the rotations IR, the IRth bit sets to 1        ==
  ! == ISISHEL(NHG) !!!WARNING INTEGER*8                            ==
  ! ==          If TSYMMORPHIC unused                               ==
  ! ==          If .NOT.TSYMMORPHIC                                 ==
  ! ==          For equivalent shell, give for each G-vector        ==
  ! ==          if IG or -IG has to take. The IRth bits sets to 1   ==
  ! ==          the component IG is taken for the IRth rotation     ==
  ! == WARNING: you need to use INTEGER*8 in order ot have 64 bits  ==
  ! ==--------------------------------------------------------------==

  ! ==================================================================

#if defined(__NOINT8)
  INTEGER, ALLOCATABLE :: ninshel(:,:)
  INTEGER, ALLOCATABLE :: isishel(:,:)
#else
  INTEGER(int_8), ALLOCATABLE :: ninshel(:)
  INTEGER(int_8), ALLOCATABLE :: isishel(:)
#endif

END MODULE symm4
