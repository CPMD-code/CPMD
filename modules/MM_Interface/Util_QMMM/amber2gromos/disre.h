C$Id: disre.h,v 1.1 2006-12-27 11:20:27 itavern Exp $ -*-fortran-*-
C constants defining behaviour of SUBR. DISRE


COMMVAR IHTREL,IHLCH1,IHRCH1,IHLCH2,IVLCH2,IH1CH3,IH2CH3,MINIHT,MAXIHT
C     Here Hydrogen type codes 
      INTEGER IHTREL,IHLCH1,IHRCH1,IHLCH2,IVLCH2,IH1CH3,IH2CH3
      INTEGER MINIHT,MAXIHT

C real hydrogen atom
      PARAMETER (IHTREL = 0)
C virtual H, aliphatic CH1
      PARAMETER (IHLCH1 = 1)
C virtual H, aromatic CH1
      PARAMETER (IHRCH1 = 2)
C pseudo H, aliphatic CH2
      PARAMETER (IHLCH2 = 3)
C virtual H, aliphatic  CH2
      PARAMETER (IVLCH2 = 4)
C pseudo H, one CH3 group
      PARAMETER (IH1CH3 = 5)
C pseudo H, two CH3 groups
      PARAMETER (IH2CH3 = 6)
      PARAMETER (MINIHT = IHTREL, MAXIHT = IH2CH3)
COMMEND



COMMVAR NDRALL,NDREF,NDREL,NDRE,NDRNDL,NDRNDE,NDRMMI,NDRMMA
C     values for NDRCTR in subroutine L<DISRE>

      INTEGER NDRALL,NDREF,NDREL,NDRE,NDRNDL,NDRNDE,NDRMMI,NDRMMA
      PARAMETER (NDRALL =-1)
      PARAMETER (NDREF  = 0)
      PARAMETER (NDREL  = 1)
      PARAMETER (NDRE   = 2)
      PARAMETER (NDRNDL = 3)
      PARAMETER (NDRNDE = 4)

C NDRMIN and NDRMAX are already defined...
      PARAMETER (NDRMMI = NDRALL,NDRMMA = NDRNDE)

C     NDRCTR is   NB means     NM means  
C     -1          num of d.r.  num of molecules    EB,F,XB0,EB0
C     0           num of d.r.  num of molecules    EB,F
C     1           num of d.r.  num of molecules    EB,XB0,EB0
C     2           num of d.r.  num of molecules    EB
C     3           seq number   seq number molecule calc length
C     4           seq number   seq number molecule calc energy
COMMEND



