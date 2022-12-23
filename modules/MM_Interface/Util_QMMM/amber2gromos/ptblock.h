C $Id: ptblock.h,v 1.1 2006-12-27 11:30:34 itavern Exp $ -*-fortran-*-
C define block types for perturbation files



      INTEGER NPTTIT, NPTATM,NPTATP,NPTBHG,NPTBNG,NPTBAH
      INTEGER NPTBAG,NPTQDH,NPTQDG,NPTPDH,NPTPDG
      INTEGER NPTPIN,NPTEND


      PARAMETER (NPTTIT =  1)
      PARAMETER (NPTATM =  2)
      PARAMETER (NPTATP =  3)

      PARAMETER (NPTBHG =  4)
      PARAMETER (NPTBNG =  5)

      PARAMETER (NPTBAH =  6)
      PARAMETER (NPTBAG =  7)

      PARAMETER (NPTQDH =  8)
      PARAMETER (NPTQDG =  9)

      PARAMETER (NPTPDH = 10)
      PARAMETER (NPTPDG = 11)

      PARAMETER (NPTPIN = 12)
      PARAMETER (NPTEND = 13)

C add other block types here.
C remember to adjust NPTMAX in accordance

      INTEGER NPTERR,NPTMIN, NPTMAX
      PARAMETER (NPTERR = 0,NPTMIN = NPTTIT, NPTMAX = NPTEND)


C maximum length of identifier strings
      INTEGER MXPTST
      PARAMETER (MXPTST = 30)
      LOGICAL LPTGOT
      COMMON /PTBLCK/ LPTGOT(NPTMAX)

      CHARACTER *(MXPTST) PTBNAM
      COMMON /PTCHAR/PTBNAM(NPTMAX)






