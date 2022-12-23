C$Id: nbpml.h,v 1.1 2006-12-27 11:27:56 itavern Exp $  -*-fortran-*-
C common blocks for SUBR. NBPML

C-------------------------------------------
C INCLUDE forcesz.h BEFORE THIS FILE!
C-------------------------------------------


C     
C     COMMON BLOCK FOR NBPML
C     
      LOGICAL LPERTL
      LOGICAL LDOTRA,LMONO,LOCTO,LVAC,LDOVIR,L4D
      LOGICAL LPIDOP

      INTEGER NSPT

      REAL COSB,COSB2,BOXOH,BOXOQ,RCUTP2,RCUTL2
      REAL RFF, RFE, RFC, C1, RCRF2, RMULOC
      REAL PININV

      COMMON /PMLINT/
     $     LDOTRA,LMONO,LOCTO,LVAC,LDOVIR,L4D,LPIDOP,
     $     LPERTL

      COMMON /PMINT/
     $     NSPT(MAXINB)

      COMMON /PMLREL/
     $     RFF, RFE, RFC, C1, RCRF2, 
     $     COSB,COSB2,BOXOH,BOXOQ,RCUTP2,RCUTL2,PININV,
     $     RMULOC

