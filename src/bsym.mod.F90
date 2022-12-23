MODULE bsym
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == BROKEN SYMMETRY                                              ==
  ! ==================================================================
  ! == HSPIN: MULTIPLICITY OF THE HIGHSPIN STATE                    ==
  ! == HUPEL HDOEL: number of up and down electrons, HS state       ==
  ! == BSCLCS: BS CaLCulation State (1:BS WF, 2:HS WF)              ==
  ! == BSFAC: Factor to account for to wfncts in memory allocations ==
  ! == RESTBS: True if only BS state in RESTART file                ==
  ! == SMIN,SMAX: Lowest and highest spins in the BS system         ==
  ! == SA,SB: Spins on coupled centers A and B                      ==
  ! == RTSASB,CNSTWGT: Derived values                               ==
  ! == EKINC_HS,EKINC_BS: Electronic kinetic energy of HS and BS WF ==
  ! ==--------------------------------------------------------------==
  INTEGER :: hspin
  INTEGER :: hupel
  INTEGER :: hdoel
  INTEGER :: hsup
  INTEGER :: hsdown
  INTEGER :: bsclcs
  INTEGER :: bsfac
  LOGICAL :: restbs
  LOGICAL :: resthswf
  REAL(real_8) :: smin
  REAL(real_8) :: smax
  REAL(real_8) :: sa
  REAL(real_8) :: sb
  REAL(real_8) :: rtsasb
  REAL(real_8) :: cnstwgt
  REAL(real_8) :: ekinc_hs
  REAL(real_8) :: ekinc_bs
  ! ==================================================================
  ! == BROKEN SYMMETRY CONSTANTS                                    ==
  ! ==================================================================
  ! == AUTOCM: Conversion factor a.u. to cm^-1 (AUTOCM = 219474.6)  ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: autocm
  ! ==================================================================

END MODULE bsym
