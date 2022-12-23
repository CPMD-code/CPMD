C $Id: topommt2.h,v 1.1 2006-12-27 11:34:22 itavern Exp $
C

C     COMMON BLOCK DEFINITIONS FOR THE 2. GROMOS TOPOLOGY IN PROMMT.F
C
C----------------------------------------------------------
C     NOTE: INCLUDE THE FILE TOPOSZ.H BEFORE INCLUDING
C     THIS ONE.
C----------------------------------------------------------


      CHARACTER*(MAXNLE) PANM2
      CHARACTER*(MAXRLE) AANM2

      INTEGER IAC2,NRP2,NRAT2,NRAT22,
     $        NRAA22,
     $        NBTY2,
     $        IB2,JB2,ICB2,NBON2,
     $        IBH2,JBH2,ICBH2,NBONH2,
     $        NTTY2,
     $        ITH2,JTH2,KTH2,ICTH2,NTHEH2,
     $        IT2,JT2,KT2,ICT2,NTHE2,
     $        NQTY2,
     $        IQ2,JQ2,KQ2,LQ2,ICQ2,NQHI2,
     $        IQH2,JQH2,KQH2,LQH2,ICQH2,NQHIH2,
     $        NPTY2,
     $        IP2,JP2,KP2,LP2,ICP2,NPHI2,
     $        IPH2,JPH2,KPH2,LPH2,ICPH2,NPHIH2,
     $        MRES2,INC2,NCAG2,INE2,KNE2,INE24,KNE24,
     $        JSNE2,NAEX2,JSNE24,NAEX42

      REAL    CG2,WINV2,WMAS2

      COMMON /CTOP2/
     $        PANM2(MAXNRP),
     $        AANM2(MAXAA2)

      COMMON /ITOP2/
     $       IAC2(MAXNRP),NRP2,
     $       NRAT2,NRAT22,
     $       NRAA22,
     $       NBTY2,
     $       IB2(MAXBON),JB2(MAXBON),ICB2(MAXBON),NBON2,
     $       IBH2(MAXBNH),JBH2(MAXBNH),ICBH2(MAXBNH),NBONH2,
     $       ITH2(MXQHEH),JTH2(MXQHEH),KTH2(MXQHEH),
     $       ICTH2(MXQHEH),NTHEH2,
     $       IT2(MAXTHE),JT2(MAXTHE),KT2(MAXTHE),
     $       ICT2(MAXTHE),NTHE2,NTTY2,
     $       NQTY2
      
      COMMON /IITOP2/
     $       IQ2(MAXQHI),JQ2(MAXQHI),KQ2(MAXQHI),LQ2(MAXQHI),
     $       ICQ2(MAXQHI),NQHI2,
     $       IQH2(MAXHIH),JQH2(MAXHIH),
     $       KQH2(MAXHIH),LQH2(MAXHIH),ICQH2(MAXHIH),NQHIH2,
     $       NPTY2,
     $       IP2(MAXPHI),JP2(MAXPHI),KP2(MAXPHI),LP2(MAXPHI),
     $       ICP2(MAXPHI),NPHI2,
     $       IPH2(MXPHIH),JPH2(MXPHIH),KPH2(MXPHIH),LPH2(MXPHIH),
     $       ICPH2(MXPHIH),NPHIH2,
     $       MRES2(MAXNRP),
     $       INC2(MAXCAG),NCAG2,
     $       INE2(MAXNRP),KNE2(MAXNRP),INE24(MAXNRP),KNE24(MAXNRP),
     $       JSNE2(MAXAEX),NAEX2,
     $       JSNE24(MXEX14),NAEX42


      COMMON /RTOP2/
     $        CG2(MAXNRP),WINV2(MAXNRP),WMAS2(MAXNRP)
