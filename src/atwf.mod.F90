MODULE atwf
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR GENERATION OF INITIAL ATOMIC WAVEFUNCTION   ==
  ! == (ATOMWF,ATRHO,SETBASIS,...)                                  ==
  ! ==--------------------------------------------------------------==
  ! == NSHELL(MAXSP)         Number of shell n per species          ==
  ! == LSHELL(M1SHLX,MAXSP)  Number of orbital per species          ==
  ! == NUMAOR(MAXSP)         Number of  orbitals per species        ==
  ! == NBCUT                 0:WFN cutoff  1:DENSITY cutoff         ==
  ! ==================================================================
  INTEGER, PARAMETER :: m1shlx=20,m1spli=256,loadc_foc_array_size=256
  ! ==--------------------------------------------------------------==
  ! == M1SHL                                                        ==
  ! == NATTOT                Total number of  orbitals              ==
  ! == NUMAORMAX             Maximum of NUMAOR(MAXSP)               ==
  ! ==--------------------------------------------------------------==
  TYPE :: atwp_t
     INTEGER :: m1shl=HUGE(0)
     INTEGER :: nattot=HUGE(0)
     INTEGER :: numaormax=HUGE(0)
  END TYPE atwp_t
  TYPE(atwp_t), SAVE :: atwp
  ! ==--------------------------------------------------------------==
  ! == OC(M1SHLX,MAXSP)                                             ==
  ! == STOEXP(M1SHLX,MAXSP)                                         ==
  ! == NQSTO(M1SHLX,MAXSP)                                          ==
  ! == LSHELL(M1SHLX,MAXSP)                                         ==
  ! == NUMAOR(MAXSP)                                                ==
  ! == NDSHELL(MAXSP)                                               ==
  ! == NBCUT                                                        ==
  ! ==--------------------------------------------------------------==
  TYPE :: atwf_mod_t
     REAL(real_8) :: oc(m1shlx,maxsp)=HUGE(0.0_real_8)
     REAL(real_8) :: stoexp(m1shlx,maxsp)=HUGE(0.0_real_8)
     INTEGER :: nqsto(m1shlx,maxsp)=HUGE(0)
     INTEGER :: nshell(maxsp)=HUGE(0)
     INTEGER :: lshell(m1shlx,maxsp)=HUGE(0)
     INTEGER :: numaor(maxsp)=HUGE(0)
     INTEGER :: ndshell(maxsp)=HUGE(0)
     INTEGER :: nbcut=HUGE(0)
  END TYPE atwf_mod_t
  TYPE(atwf_mod_t), SAVE :: atwf_mod
  ! ==--------------------------------------------------------------==
  ! == CAT                                                          ==
  ! == XXMAT(NATTOT,NATTOT)                                         ==
  ! == XSMAT(NATTOT,NATTOT)                                         ==
  ! == DBAS                                                         ==
  ! == CATOM                                                        ==
  ! == ZXMAT(NATTOT,NATTOT)                                         ==
  ! == ZSMAT(NATTOT,NATTOT)                                         ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: cat(:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: xxmat(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: xsmat(:,:)

  COMPLEX(real_8), ALLOCATABLE, TARGET, SAVE :: catom(:,:)
  COMPLEX(real_8), ALLOCATABLE, SAVE :: zxmat(:,:)
  COMPLEX(real_8), ALLOCATABLE, SAVE :: zsmat(:,:)
  ! ==--------------------------------------------------------------==
  ! == MESHAT                                                       ==
  ! == CLOGAT                                                       ==
  ! == ATRG                                                         ==
  ! == ATWFR                                                        ==
  ! == ATCHG(NSX) Atomic charge                                     ==
  ! ==--------------------------------------------------------------==

  REAL(real_8), ALLOCATABLE, SAVE :: atrg(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: atwfr(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: atchg(:)
  REAL(real_8), ALLOCATABLE, SAVE :: atrg_epr(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: atwfr_epr(:,:,:)

  TYPE :: atwr_t
     INTEGER :: meshat(maxsp)=HUGE(0)
     INTEGER :: mesh_epr(maxsp)=HUGE(0)
     REAL(real_8) :: clogat(maxsp)=HUGE(0.0_real_8)
  END TYPE atwr_t
  TYPE(atwr_t), SAVE :: atwr
  ! ==--------------------------------------------------------------==
  ! == DMOVRMIX Parameter for MOVERHO (mixing when density moves)   ==
  ! == TMOVR    .TRUE. uses MOVERHO to move density when moving     ==
  ! ==          atoms (projects density to atomic orbital basis set)==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: dmovrmix=HUGE(0.0_real_8)
  LOGICAL, SAVE :: tmovr=.FALSE.

  ! ==================================================================
  ! CDFT
  ! ------------------------------------------------------------------
  TYPE :: cdftoc_t
     REAL(real_8) :: ocun(m1shlx,maxsp)=HUGE(0.0_real_8)
     REAL(real_8) :: oc2(m1shlx,maxsp)=HUGE(0.0_real_8)
  END TYPE cdftoc_t
  TYPE(cdftoc_t), SAVE :: cdftoc
  TYPE :: cdftocl_t
     LOGICAL :: ocset=.FALSE.
  END TYPE cdftocl_t
  TYPE(cdftocl_t), SAVE :: cdftocl
  ! == CATORD(M1SHLX,NAX)  Order of atomic orbital shells in CAT ==
  INTEGER, ALLOCATABLE, SAVE :: catord(:,:)
  INTEGER, SAVE :: natsave=HUGE(0)
  ! ==================================================================

END MODULE atwf
