MODULE ropt
  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR PRINTP (TRAJECTORY, STRESS, MOVIE FILES)    ==
  ! ==================================================================
  ! == TESR    .TRUE. Ewald sum for electrostatic energy            ==
  ! == SDIIS   .TRUE. Initialize or Reinitialize DIIS scheme        ==
  ! == SPCG    .TURE. Call of HESELE                                ==
  ! == CONVWF  CONVergence of WaveFunctions                         ==
  ! == CONVGE  CONVergence of GEometry                              ==
  ! == PRTEIG  PRinT EIGenvalues                                    ==
  ! == ENGPRI  ENerGy PRInt                                         ==
  ! == MOVIE   Write movie film                                     ==
  ! == TXYZ    Write xyz trajectory film                            ==
  ! == TDCD    Write dcd trajectory film                            ==
  ! == CALSTE  .TRUE. calculations of stress tensor                 ==
  ! ==         for the iteration                                    ==
  ! == CALSTC  .TRUE. calculations of stress tensor                 ==
  ! ==         for the given iteration                              ==
  ! == RPRINT  Save DENSITY file for the given iteration            ==
  ! == MODENS  Use atomic orbital for updated electronic density    ==
  ! ==         for the next iteration                               ==
  ! ==--------------------------------------------------------------==

  TYPE :: ropt_mod_t
     LOGICAL :: tesr
     LOGICAL :: sdiis
     LOGICAL :: spcg
     LOGICAL :: convwf
     LOGICAL :: convge
     LOGICAL :: prteig
     LOGICAL :: engpri
     LOGICAL :: txyz
     LOGICAL :: tdcd
     LOGICAL :: movie
     LOGICAL :: calste
     LOGICAL :: calstc
     LOGICAL :: rprint
     LOGICAL :: modens
  END TYPE ropt_mod_t
  TYPE(ropt_mod_t) :: ropt_mod
  ! ==================================================================
  ! == NFI  Total number of iteration (for wavefunction)            ==
  ! == INFI Total number of self-iteration (GEO or MD)              ==
  ! == IINFI Total number of self-iteration QUENCHBO
  ! == INFW Number of iterations for wavefunction (during 1 SC iter)==
  ! == IESR Number of lattice point (IESRxIESRxIESR) for Ewald sum  ==
  ! == NDISRS Number of DIIS resets since last diagonalization      ==
  ! == NDISTP Number of DIIS steps since last DIIS reset            ==
  ! ==--------------------------------------------------------------==
  TYPE :: iteropt_t
     INTEGER :: nfi
     INTEGER :: infi
     INTEGER :: iinfi
     INTEGER :: infw
     INTEGER :: iesr
     INTEGER :: ndisrs
     INTEGER :: ndistp
  END TYPE iteropt_t
  TYPE(iteropt_t), TARGET :: iteropt
  INTEGER, POINTER :: infi, infw ! inf[iw] cannot be DO loop iterator
  ! ==================================================================
  ! == BSNFI  Total number of iteration (for BS wavefunction)       ==
  ! ==--------------------------------------------------------------==
  INTEGER :: bsnfi
  TYPE :: bsiopt_t
     INTEGER :: bsnfi
  END TYPE bsiopt_t
  TYPE(bsiopt_t) :: bsiopt
  ! ==================================================================
CONTAINS
  SUBROUTINE init_pinf_pointers()

    infi => iteropt%infi
    infw => iteropt%infw
  END SUBROUTINE init_pinf_pointers
END MODULE ropt
