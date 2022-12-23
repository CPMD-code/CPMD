C $Id: fileio.h,v 1.1 2006-12-27 11:23:49 itavern Exp $ -*-fortran-*-

COMMVAR MAXDLN
C     maximum length of a line in a formatted file
      INTEGER MAXDLN
      PARAMETER (MAXDLN = 150)
COMMEND

COMMVAR NDXRI,NDXLE
C     used by the chopping routines
      INTEGER NDXRI,NDXLE
      COMMON /FIOINT/NDXRI,NDXLE
COMMEND

COMMVAR FIOLIN
C     the line read in by L<RAWGET>
      CHARACTER*(MAXDLN) FIOLIN
      COMMON /FIOCHR/FIOLIN
COMMEND

COMMVAR FMRMAX,FMIMAX,FMTLEN
C number of real and integer format statements
      INTEGER FMRMAX,FMIMAX
      PARAMETER (FMRMAX = 15, FMIMAX = 10)
C length of format statements
      INTEGER FMTLEN
      PARAMETER (FMTLEN = 8)
COMMEND

COMMVAR FMRTAB,FMITAB
C     these strings are initialized in a BLOCKDATA
C     statement in fileio.f
C
      CHARACTER*(FMTLEN) FMRTAB
      CHARACTER*(FMTLEN) FMITAB
      COMMON /FIOFMT/FMRTAB(FMRMAX),FMITAB(FMIMAX)
COMMEND


