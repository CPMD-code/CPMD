C  $Id: ffar.h,v 1.1 2006-12-27 11:22:27 itavern Exp $ -*-fortran-*-
C  force field arrays
C

C only define things here which are not defined in toposz.h
C----------------------------------------------------------
C remember to include ffsz.h and toposz.h before this one
C----------------------------------------------------------

COMMVAR FFTIT,NFFTLN
C     force field title
      CHARACTER FFTIT*(MAXTIT)
      COMMON /FFCHAR/FFTIT(MAXLNS)
      INTEGER NFFTLN
      COMMON /FFTLN/NFFTLN
COMMEND


COMMVAR NUMIMC,IMC,ATMAS,ATMASN
C     mass codes, names and values
      INTEGER NUMIMC,IMC
      REAL ATMAS
      CHARACTER ATMASN*(MAXTLE)

      COMMON /IMCINT/NUMIMC,IMC(MAXIMC)
      COMMON /IMCREL/ATMAS(MAXIMC)
      COMMON /IMCCHR/ATMASN(MAXIMC)
COMMEND



COMMVAR C612,C1212,CS612,CS1212,LPAIR
C     for atom interaction parameters
      REAL C612,C1212
      REAL CS612,CS1212

      INTEGER LPAIR
      COMMON /C12REL/C612(MAXATT),C1212(MAXC12,MAXATT),
     $     CS612(MAXATT),CS1212(MAXATT)
      COMMON /C12INT/LPAIR(MAXATT,MAXATT)
COMMEND

COMMVAR IMIX,JMIX,NRMIX,C6MIX,C12MIX,CS6MIX,CS12MI
C     for mix atom pairs
      INTEGER IMIX,JMIX,NRMIX
      REAL C6MIX,C12MIX,CS6MIX,CS12MI

      COMMON /MIXINT/IMIX(MAXMIX),JMIX(MAXMIX),NRMIX
      COMMON /MIXREL/ C6MIX(MAXMIX),C12MIX(MAXMIX),
     $     CS6MIX(MAXMIX),CS12MI(MAXMIX)
COMMEND

