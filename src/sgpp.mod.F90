MODULE sgpp
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: lmaxx,&
                                             maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == STEFAN GOEDECKER PSEUDOPOTENIAL                              ==
  ! ==================================================================
  ! == See Goedecker et al., PRB 54 (1996) pp 1703--1710            ==
  ! ==================================================================
  ! == MPRO  Maximum of projectors for each angular momentum        ==
  ! == MCFUN Maximum number of C coefficients for the local part    ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: mpro=4,mcfun=4 
  ! ==================================================================
  ! == TSGPALL      True if ALL PP are SG PP                        ==
  ! == TSGP(MAXSP)  For each specie, true if SG PP                  ==
  ! == NCLSG(MAXSP) Number of C coefficient for local part          ==
  ! ==--------------------------------------------------------------==


  ! ==================================================================
  ! == RCSG(MAXSP) R_{loc} local radius                             ==
  ! == CLSG(MCFUN,MAXSP) C coefficients for the local potential part==
  ! == RCNL(LMAXX,MAXSP) Radius for each ang. mom. quantum number   ==
  ! == HLSG(MPRO,MPRO,LMAXX,MAXSP): For each specie give the matrix ==
  ! ==     h^l_(i,j) i=1:NPRO(LMAXX,MAXSP)                          ==
  ! == NPRO(LMAXX,MAXSP) Number of projectors                       ==
  ! == LPVAL(NHX,MAXSP) For each projector gives                    ==
  ! ==      the l*m number                                          ==
  ! == LFVAL(NHX,MAXSP) For each projector gives the                ==
  ! ==      projector number (from 1 to NPRO(LMAXX,MAXSP))          ==
  ! == PPLIST(LMAX*MPRO,2,MAXSP) For each projector gives the       ==
  ! ==      value l and number of projector within l                ==
  ! ==--------------------------------------------------------------==


  ! ==================================================================
  TYPE :: sgpp1_t
     LOGICAL :: tsgpall
     LOGICAL :: tsgp(maxsp)
     INTEGER :: nclsg(maxsp)
  END TYPE sgpp1_t
  TYPE(sgpp1_t) :: sgpp1
  TYPE :: sgpp2_t
     REAL(real_8) :: rcsg(maxsp)
     REAL(real_8) :: clsg(mcfun,maxsp)
     REAL(real_8) :: rcnl(lmaxx,maxsp)
     REAL(real_8) :: hlsg(mpro,mpro,lmaxx,maxsp)
     INTEGER :: npro(lmaxx,maxsp)
     INTEGER :: lpval(lmaxx*lmaxx*mpro,maxsp)
     INTEGER :: lfval(lmaxx*lmaxx*mpro,maxsp)
     INTEGER :: pplist(lmaxx*mpro,2,maxsp)
  END TYPE sgpp2_t
  TYPE(sgpp2_t) :: sgpp2

END MODULE sgpp
