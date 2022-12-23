C $Id: feear.h,v 1.1 2006-12-27 11:22:09 itavern Exp $ -*-fortran-*-
C     Arrays needed in program PROFEE and associated subroutines

C-----------------------------------------------------------------------
C     remember to include toposz.h, md.h and forcesz.h before this one
C-----------------------------------------------------------------------


C title of fee files (read in inifee)

      CHARACTER FEETIT*(MAXTIT)
      INTEGER NTITF
      COMMON /FTIT/FEETIT(MAXLNS)
      COMMON /FTITL/NTITF

C atoms to be created
C maxcrt is the maxmum number of atoms to be created
      INTEGER MAXCRT
      PARAMETER (MAXCRT = 10)

      INTEGER NCRT,IA,JA,KA,LA
      REAL BOND,ALPHA,PHI
      COMMON /COMCRT/BOND(MAXCRT),ALPHA(MAXCRT),
     $               PHI(MAXCRT),
     $               NCRT,IA(MAXCRT),JA(MAXCRT),
     $               KA(MAXCRT),LA(MAXCRT)
            


C Common blocks used to transfer variables energies

C   nonbonded atom group energy components
      REAL EPLJ0,EPLJAV0,EPLJFL0,
     $     EPLJ1,EPLJAV1,EPLJFL1,
     $     EPLJ2,EPLJAV2,EPLJFL2,
     $     EPLJ3,EPLJAV3,EPLJFL3

      COMMON /FEPLJ/EPLJ0(MXNRE2),EPLJAV0(MXNRE2),EPLJFL0(MXNRE2),
     $              EPLJ1(MXNRE2),EPLJAV1(MXNRE2),EPLJFL1(MXNRE2),
     $              EPLJ2(MXNRE2),EPLJAV2(MXNRE2),EPLJFL2(MXNRE2),
     $              EPLJ3(MXNRE2),EPLJAV3(MXNRE2),EPLJFL3(MXNRE2)

      DOUBLE PRECISION EPLJSU0,EPLJSQ0,
     $                 EPLJSU1,EPLJSQ1,
     $                 EPLJSU2,EPLJSQ2,
     $                 EPLJSU3,EPLJSQ3

      COMMON /FEPLJD/EPLJSU0(MXNRE2),EPLJSQ0(MXNRE2),
     $               EPLJSU1(MXNRE2),EPLJSQ1(MXNRE2),
     $               EPLJSU2(MXNRE2),EPLJSQ2(MXNRE2),
     $               EPLJSU3(MXNRE2),EPLJSQ3(MXNRE2)



      REAL EPEL0,EPELAV0,EPELFL0,
     $     EPEL1,EPELAV1,EPELFL1,
     $     EPEL2,EPELAV2,EPELFL2,
     $     EPEL3,EPELAV3,EPELFL3

      COMMON /FEPEL/EPEL0(MXNRE2),EPELAV0(MXNRE2),EPELFL0(MXNRE2),
     $              EPEL1(MXNRE2),EPELAV1(MXNRE2),EPELFL1(MXNRE2),
     $              EPEL2(MXNRE2),EPELAV2(MXNRE2),EPELFL2(MXNRE2),
     $              EPEL3(MXNRE2),EPELAV3(MXNRE2),EPELFL3(MXNRE2)

      DOUBLE PRECISION EPELSU0,EPELSQ0,
     $                 EPELSU1,EPELSQ1,
     $                 EPELSU2,EPELSQ2,
     $                 EPELSU3,EPELSQ3

      COMMON /FEPELD/EPELSU0(MXNRE2),EPELSQ0(MXNRE2),
     $               EPELSU1(MXNRE2),EPELSQ1(MXNRE2),
     $               EPELSU2(MXNRE2),EPELSQ2(MXNRE2),
     $               EPELSU3(MXNRE2),EPELSQ3(MXNRE2)


      REAL EPRF0,EPRFAV0,EPRFFL0,
     $     EPRF1,EPRFAV1,EPRFFL1,
     $     EPRF2,EPRFAV2,EPRFFL2,
     $     EPRF3,EPRFAV3,EPRFFL3

      COMMON /FEPRF/EPRF0(MXNRE2),EPRFAV0(MXNRE2),EPRFFL0(MXNRE2),
     $              EPRF1(MXNRE2),EPRFAV1(MXNRE2),EPRFFL1(MXNRE2),
     $              EPRF2(MXNRE2),EPRFAV2(MXNRE2),EPRFFL2(MXNRE2),
     $              EPRF3(MXNRE2),EPRFAV3(MXNRE2),EPRFFL3(MXNRE2)


      DOUBLE PRECISION EPRFSU0,EPRFSQ0,
     $                 EPRFSU1,EPRFSQ1,
     $                 EPRFSU2,EPRFSQ2,
     $                 EPRFSU3,EPRFSQ3

      COMMON /FEPRFD/EPRFSU0(MXNRE2),EPRFSQ0(MXNRE2),
     $               EPRFSU1(MXNRE2),EPRFSQ1(MXNRE2),
     $               EPRFSU2(MXNRE2),EPRFSQ2(MXNRE2),
     $               EPRFSU3(MXNRE2),EPRFSQ3(MXNRE2)


      REAL EPRC0,EPRCAV0,EPRCFL0,
     $     EPRC1,EPRCAV1,EPRCFL1,
     $     EPRC2,EPRCAV2,EPRCFL2,
     $     EPRC3,EPRCAV3,EPRCFL3

      COMMON /FEPRC/EPRC0(MXNRE2),EPRCAV0(MXNRE2),EPRCFL0(MXNRE2),
     $              EPRC1(MXNRE2),EPRCAV1(MXNRE2),EPRCFL1(MXNRE2),
     $              EPRC2(MXNRE2),EPRCAV2(MXNRE2),EPRCFL2(MXNRE2),
     $              EPRC3(MXNRE2),EPRCAV3(MXNRE2),EPRCFL3(MXNRE2)


      DOUBLE PRECISION EPRCSU0,EPRCSQ0,
     $                 EPRCSU1,EPRCSQ1,
     $                 EPRCSU2,EPRCSQ2,
     $                 EPRCSU3,EPRCSQ3

      COMMON /FEPRCD/EPRCSU0(MXNRE2),EPRCSQ0(MXNRE2),
     $               EPRCSU1(MXNRE2),EPRCSQ1(MXNRE2),
     $               EPRCSU2(MXNRE2),EPRCSQ2(MXNRE2),
     $               EPRCSU3(MXNRE2),EPRCSQ3(MXNRE2)


C     energy, dE/dLamda and dE/dmu (original + triplicate arrays)

      REAL ENER0,EPAVE0,EPFLC0,
     $     ENER1,EPAVE1,EPFLC1,
     $     ENER2,EPAVE2,EPFLC2,
     $     ENER3,EPAVE3,EPFLC3

      COMMON /FENER/ENER0(MXETBL),EPAVE0(MXETBL),EPFLC0(MXETBL),
     $              ENER1(MXETBL),EPAVE1(MXETBL),EPFLC1(MXETBL),
     $              ENER2(MXETBL),EPAVE2(MXETBL),EPFLC2(MXETBL),
     $              ENER3(MXETBL),EPAVE3(MXETBL),EPFLC3(MXETBL)


      DOUBLE PRECISION EPSUM0,EPSQ0,
     $                 EPSUM1,EPSQ1,
     $                 EPSUM2,EPSQ2,
     $                 EPSUM3,EPSQ3

      COMMON /FENERD/EPSUM0(MXETBL),EPSQ0(MXETBL),
     $               EPSUM1(MXETBL),EPSQ1(MXETBL),
     $               EPSUM2(MXETBL),EPSQ2(MXETBL),
     $               EPSUM3(MXETBL),EPSQ3(MXETBL)



      REAL DEDLAM0,DEDLAV0,DEDLFL0,
     $     DEDLAM1,DEDLAV1,DEDLFL1,
     $     DEDLAM2,DEDLAV2,DEDLFL2,
     $     DEDLAM3,DEDLAV3,DEDLFL3,
     $     DEDMU0,DEDMAV0,DEDMFL0,
     $     DEDMU1,DEDMAV1,DEDMFL1,
     $     DEDMU2,DEDMAV2,DEDMFL2,
     $     DEDMU3,DEDMAV3,DEDMFL3


      COMMON /FDED/DEDLAM0(MXETBL),DEDLAV0(MXETBL),DEDLFL0(MXETBL),
     $             DEDLAM1(MXETBL),DEDLAV1(MXETBL),DEDLFL1(MXETBL),
     $             DEDLAM2(MXETBL),DEDLAV2(MXETBL),DEDLFL2(MXETBL),
     $             DEDLAM3(MXETBL),DEDLAV3(MXETBL),DEDLFL3(MXETBL),
     $             DEDMU0(MXETBL),DEDMAV0(MXETBL),DEDMFL0(MXETBL),
     $             DEDMU1(MXETBL),DEDMAV1(MXETBL),DEDMFL1(MXETBL),
     $             DEDMU2(MXETBL),DEDMAV2(MXETBL),DEDMFL2(MXETBL),
     $             DEDMU3(MXETBL),DEDMAV3(MXETBL),DEDMFL3(MXETBL)

      DOUBLE PRECISION DEDLSU0,DEDLSQ0,
     $                 DEDLSU1,DEDLSQ1,
     $                 DEDLSU2,DEDLSQ2,
     $                 DEDLSU3,DEDLSQ3,
     $                 DEDMSU0,DEDMSQ0,
     $                 DEDMSU1,DEDMSQ1,
     $                 DEDMSU2,DEDMSQ2,
     $                 DEDMSU3,DEDMSQ3


      COMMON /FDEDD/DEDLSU0(MXETBL),DEDLSQ0(MXETBL),
     $              DEDLSU1(MXETBL),DEDLSQ1(MXETBL),
     $              DEDLSU2(MXETBL),DEDLSQ2(MXETBL),
     $              DEDLSU3(MXETBL),DEDLSQ3(MXETBL),
     $              DEDMSU0(MXETBL),DEDMSQ0(MXETBL),
     $              DEDMSU1(MXETBL),DEDMSQ1(MXETBL),
     $              DEDMSU2(MXETBL),DEDMSQ2(MXETBL),
     $              DEDMSU3(MXETBL),DEDMSQ3(MXETBL)


C     arrays for additional force field terms and averages thereof

      REAL ENERES0, EREAVE0, EREFLC0,
     $     ENERES1, EREAVE1, EREFLC1,
     $     ENERES2, EREAVE2, EREFLC2,
     $     ENERES3, EREAVE3, EREFLC3

      COMMON /FENERE/ENERES0(MXCTBL),EREAVE0(MXCTBL),EREFLC0(MXCTBL),
     $               ENERES1(MXCTBL),EREAVE1(MXCTBL),EREFLC1(MXCTBL),
     $               ENERES2(MXCTBL),EREAVE2(MXCTBL),EREFLC2(MXCTBL),
     $               ENERES3(MXCTBL),EREAVE3(MXCTBL),EREFLC3(MXCTBL)

      DOUBLE PRECISION ERESUM0,ERESQ0,
     $                 ERESUM1,ERESQ1,
     $                 ERESUM2,ERESQ2,
     $                 ERESUM3,ERESQ3

      COMMON /FENERED/ERESUM0(MXCTBL),ERESQ0(MXCTBL),
     $                ERESUM1(MXCTBL),ERESQ1(MXCTBL),
     $                ERESUM2(MXCTBL),ERESQ2(MXCTBL),
     $                ERESUM3(MXCTBL),ERESQ3(MXCTBL)


C VOLPRT: some sytem parameters saved to file.
      REAL VOLPRT0,VPAVE0,VPFLC0,
     $     VOLPRT1,VPAVE1,VPFLC1,
     $     VOLPRT2,VPAVE2,VPFLC2,
     $     VOLPRT3,VPAVE3,VPFLC3

      COMMON /FVOLPR/VOLPRT0(MXVTBL),VPAVE0(MXVTBL),VPFLC0(MXVTBL),
     $               VOLPRT1(MXVTBL),VPAVE1(MXVTBL),VPFLC1(MXVTBL),
     $               VOLPRT2(MXVTBL),VPAVE2(MXVTBL),VPFLC2(MXVTBL),
     $               VOLPRT3(MXVTBL),VPAVE3(MXVTBL),VPFLC3(MXVTBL)

      DOUBLE PRECISION VPSUM0,VPSQ0,
     $                 VPSUM1,VPSQ1,
     $                 VPSUM2,VPSQ2,
     $                 VPSUM3,VPSQ3

      COMMON /FVOLPRD/VPSUM0(MXVTBL),VPSQ0(MXVTBL),
     $                VPSUM1(MXVTBL),VPSQ1(MXVTBL),
     $                VPSUM2(MXVTBL),VPSQ2(MXVTBL),
     $                VPSUM3(MXVTBL),VPSQ3(MXVTBL)



