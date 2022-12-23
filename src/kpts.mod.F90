MODULE kpts
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INFORMATION ABOUT K POINTS MESH (MONKHORST-PACK)             ==
  ! ==   NK1,NK2,NK3 : DIMENSIONS OF MONKHORST-PACK MESH            ==
  ! == LOGICAL FLAGS                                                ==
  ! ==   TMONKP      : WE USE A MONKHORST-PACK MESH                 ==
  ! ==   TKPNT       : K POINTS USED (COMPLEX WAVEFUNCTIONS)        ==
  ! ==   TKSCALE     : EACH THREE COMPONENTS (INDEX I) ARE DIVIDED  ==
  ! ==                 NOT BY 2PI/ALAT BUT 2PI/ALAT(I)              ==
  ! ==   TSYMKP      : SYMMETRIZE SPECIAL K-POINTS                  ==
  ! ==   TKFULL      : USE THE FULL MONKHORST-PACK MESH             ==
  ! ==                 WITH ONLY INVERSION SYMMETRY                 ==
  ! ==   TKBLOCK     : CALCULATE K-POINTS PER BLOCK (LOW MEMORY)    ==
  ! ==   TKALL       : SWAP ALL ARRAYS DEPENDING ON NKPNT (BLOCK)   ==
  ! ==   TKNOSWAP    : USE WITH OPTION TKBLOCK                      ==
  ! ==                 USE TO CALCULATE K POINT PER K POINT         ==
  ! ==                 C0 IS NOT SWAPPED                            ==
  ! ==   TKBCALC     : USE WITH OPTION TKBLOCK                      ==
  ! ==                 CALCULATE EACH TIME ALL ARRAYS DEPENDING ON  ==
  ! ==                 K POINTS                                     ==
  ! ==   NKPTALL = NKPTS IF NOT TKALL OTHERWISE NKPNT               ==
  ! ==   TONLYDIAG   : ONLY ONE DIAGONALIZATION PER EACH K POINT    ==
  ! ==                 WITHOUT SELF-CONSISTENCE CYCLE               ==
  ! ==================================================================


  ! ==--------------------------------------------------------------==
  ! == WVK0(3) MAC DONALD SHIT FOR GENERATION OF MONKHORST-PACK MESH==
  ! ==         SEE MACDONALD PRB (1978) 18, 5897-5899               ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: wvk0(3)
  ! ==================================================================

  TYPE :: kpts_com_t
     INTEGER :: nk1
     INTEGER :: nk2
     INTEGER :: nk3
     INTEGER :: nkptall
  END TYPE kpts_com_t
  TYPE(kpts_com_t) :: kpts_com
  TYPE :: tkpts_t
     LOGICAL :: tkpnt
     LOGICAL :: tkscale
     LOGICAL :: tmonkp
     LOGICAL :: tsymkp
     LOGICAL :: tkfull
     LOGICAL :: tkblock
     LOGICAL :: tkall
     LOGICAL :: tknoswap
     LOGICAL :: tkbcalc
     LOGICAL :: tonlydiag
  END TYPE tkpts_t
  TYPE(tkpts_t) :: tkpts

END MODULE kpts
