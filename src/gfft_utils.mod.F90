#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE gfft_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE parac,                           ONLY: paral
  USE timer,                           ONLY: tihalt,&
                                             tiset

IMPLICIT NONE
PRIVATE

PUBLIC :: ctrig

CONTAINS

! ==================================================================
        SUBROUTINE ctrig(n,trig,after,before,now,isign,ic)
! ==--------------------------------------------------------------==
! Copyright by Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
! modified by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
! modified by Stefan Goedecker, Stuttgart, Germany, October 6, 1995
! Commercial use is prohibited 
! without the explicit permission of the author.
! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: trig(2,1024)
    INTEGER                                  :: after(7), before(7), now(7), &
                                                isign, ic

    INTEGER                                  :: i, itt, j
    INTEGER, DIMENSION(7, 100) :: idata = RESHAPE((/3,   3, 1, 1, 1, 1, 1,    &
      4,   4, 1, 1, 1, 1, 1,5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1&
      ,8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,12,   4, 3, 1, 1, 1,&
      1,      15,   5, 3, 1, 1, 1, 1,16,   4, 4, 1, 1, 1, 1,      18,   6, 3, &
      1, 1, 1, 1,20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,25,   5, &
      5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,30,   6, 5, 1, 1, 1, 1,      &
      32,   8, 4, 1, 1, 1, 1,36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1,&
      1,45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,54,   6, 3, 3, 1, &
      1, 1,      60,   5, 4, 3, 1, 1, 1,64,   4, 4, 4, 1, 1, 1,      72,   8, &
      3, 3, 1, 1, 1,75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,81,   &
      3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,96,   8, 4, 3, 1, 1, 1,   &
      100,   5, 5, 4, 1, 1, 1,108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, &
      1, 1,125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,135,   5, 3, 3&
      , 3, 1, 1,     144,   4, 4, 3, 3, 1, 1,150,   6, 5, 5, 1, 1, 1,     160,&
      8, 5, 4, 1, 1, 1,162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,&
      192,   4, 4, 4, 3, 1, 1,     200,   8, 5, 5, 1, 1, 1,216,   8, 3, 3, 3, &
      1, 1,     225,   5, 5, 3, 3, 1, 1,240,   5, 4, 4, 3, 1, 1,     243,   3,&
      3, 3, 3, 3, 1,256,   4, 4, 4, 4, 1, 1,     270,   6, 5, 3, 3, 1, 1,288, &
      8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,320,   5, 4, 4, 4, 1, 1,  &
      324,   4, 3, 3, 3, 3, 1,360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, &
      1, 1,384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,405,   5, 3, 3&
      , 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,450,   6, 5, 5, 3, 1, 1,     480,&
      8, 5, 4, 3, 1, 1,486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,&
      512,   8, 4, 4, 4, 1, 1,     540,   5, 4, 3, 3, 3, 1,576,   4, 4, 4, 3, &
      3, 1,     600,   8, 5, 5, 3, 1, 1,625,   5, 5, 5, 5, 1, 1,     640,   8,&
      5, 4, 4, 1, 1,648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,720, &
      5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,750,   6, 5, 5, 5, 1, 1,  &
      768,   4, 4, 4, 4, 3, 1,800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, &
      3, 1,864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,960,   5, 4, 4&
      , 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,1000,   8, 5, 5, 5, 1, 1,    1024&
      ,   4, 4, 4, 4, 4, 1,1080,   8, 5, 3, 3, 3, 1,1152,   8, 4, 4, 3, 3, 1, &
      1200,   5, 5, 4, 4, 3, 1,1280,   5, 4, 4, 4, 4, 1,1296,   4, 4, 3, 3, 3,&
      3,    1350,   6, 5, 5, 3, 3, 1,1440,   6, 5, 4, 4, 3, 1,    1458,   6, 3&
      , 3, 3, 3, 3,1500,   5, 5, 5, 4, 3, 1,    1536,   8, 4, 4, 4, 3, 1,1600,&
      8, 8, 5, 5, 1, 1,    1620,   5, 4, 3, 3, 3, 3,1728,   8, 8, 3, 3, 3, 1, &
      1800,   8, 5, 5, 3, 3, 1,1920,   8, 5, 4, 4, 3, 1,    1944,   8, 3, 3, 3&
      , 3, 3,2000,   5, 5, 5, 4, 4, 1,    2048,   8, 4, 4, 4, 4, 1 /),(/7,100&
      /))
    REAL(real_8)                             :: angle, twopi

      DO 111,i=1,100
        IF (n.eq.idata(1,i)) THEN
          ic=0
          DO 11,j=1,6
            itt=idata(1+j,i)
            IF (itt.gt.1) THEN
              ic=ic+1
              now(j)=idata(1+j,i)
            ELSE
              GOTO 1000
            ENDIF
 11       CONTINUE
          GOTO 1000
        ENDIF
 111  CONTINUE
      IF (paral%io_parent)&
      WRITE(6,*)'VALUE OF',n,&
           'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
 37     FORMAT(15(i5))
      IF (paral%io_parent)&
      WRITE(6,37) (idata(1,i),i=1,100)
        CALL stopgm('CTRIG',' ',& 
 __LINE__,__FILE__)
 1000   CONTINUE

        after(1)=1
        before(ic)=1
        DO 22,i=2,ic
          after(i)=after(i-1)*now(i-1)
          before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
 22     CONTINUE

! 12      format(6(i3))
! write(6,12) (after(i),i=1,ic)
! write(6,12) (now(i),i=1,ic)
! write(6,12) (before(i),i=1,ic)

        twopi=8._real_8*ATAN(1._real_8)
        angle=isign*twopi/n
        trig(1,1)=1._real_8
        trig(2,1)=0._real_8
!$omp parallel do private(i) shared(angle)
        DO i=1,n-1
          trig(1,i+1)=COS(i*angle)
          trig(2,i+1)=SIN(i*angle)
        ENDDO
! heck
! if (mod(n,2).eq.0) then
! trig(2,n/2+1)=-100000._real_8
! endif

        RETURN
        END SUBROUTINE ctrig

END MODULE gfft_utils


         SUBROUTINE fftpre(mm,nfft,m,nn,n,zin,zout,trig,now,after,before,isign)
! Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
! modified by Stefan Goedecker, Stuttgart, Germany, October 15, 1995
! Commercial use is prohibited without the explicit permission of the author.
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: mm, nfft, m, nn, n
  REAL(real_8)                               :: zin(2,m,mm), zout(2,nn,n), &
                                                trig(2,1024)
  INTEGER                                    :: now, after, before, isign

  INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
  REAL(real_8) :: am, ap, bb, bm, bp, ci2, ci3, ci4, ci5, cm, cos2, cos4, cp, &
      cr2, cr3, cr4, cr5, dm, dp, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, &
      r8, rt2i, s, s1, s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, &
      ui2, ui3, ur1, ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

        atn=after*now
        atb=after*before

! sqrt(0.5_real_8)
        rt2i=0.7071067811865475_real_8
        IF (now.eq.4) THEN
        IF (isign.eq.1) THEN
                ia=1
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                DO 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4001,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r2=zin(1,nin2,j)
                s2=zin(2,nin2,j)
                r3=zin(1,nin3,j)
                s3=zin(2,nin3,j)
                r4=zin(1,nin4,j)
                s4=zin(2,nin4,j)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r - s
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r + s
                zout(2,j,nout4) = r - s
4001                CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
                DO 4000,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4010,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=-zin(2,nin3,j)
                        s3=zin(1,nin3,j)
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
4010                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4020,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,nin3,j)
                        s=zin(2,nin3,j)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
4020                        CONTINUE
                ENDIF
4000                CONTINUE
        ELSE
                ia=1
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4101,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r2=zin(1,nin2,j)
                s2=zin(2,nin2,j)
                r3=zin(1,nin3,j)
                s3=zin(2,nin3,j)
                r4=zin(1,nin4,j)
                s4=zin(2,nin4,j)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
4101                CONTINUE
#ifdef __SR8000 
!poption indep(zout)
!voption pvfunc(3)
#endif 
                DO 4100,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                        DO 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4110,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,nin3,j)
                        s3=-zin(1,nin3,j)
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
4110                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                        DO 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4120,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,nin3,j)
                        s=zin(2,nin3,j)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
4120                        CONTINUE
                ENDIF
4100                CONTINUE
        ENDIF
        ELSE IF (now.eq.8) THEN
        IF (isign.eq.-1) THEN
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                        DO 8120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8120,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r2=zin(1,nin2,j)
                        s2=zin(2,nin2,j)
                        r3=zin(1,nin3,j)
                        s3=zin(2,nin3,j)
                        r4=zin(1,nin4,j)
                        s4=zin(2,nin4,j)
                        r5=zin(1,nin5,j)
                        s5=zin(2,nin5,j)
                        r6=zin(1,nin6,j)
                        s6=zin(2,nin6,j)
                        r7=zin(1,nin7,j)
                        s7=zin(2,nin7,j)
                        r8=zin(1,nin8,j)
                        s8=zin(2,nin8,j)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = (-cp + dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp = ( cm - dp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dp
8120                        CONTINUE
        ELSE
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                        DO 8121,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8121,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r2=zin(1,nin2,j)
                        s2=zin(2,nin2,j)
                        r3=zin(1,nin3,j)
                        s3=zin(2,nin3,j)
                        r4=zin(1,nin4,j)
                        s4=zin(2,nin4,j)
                        r5=zin(1,nin5,j)
                        s5=zin(2,nin5,j)
                        r6=zin(1,nin6,j)
                        s6=zin(2,nin6,j)
                        r7=zin(1,nin7,j)
                        s7=zin(2,nin7,j)
                        r8=zin(1,nin8,j)
                        s8=zin(2,nin8,j)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= (-cm + dp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dp
8121                        CONTINUE
        ENDIF
        ELSE IF (now.eq.3) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
        DO 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3001,j=1,nfft
        r1=zin(1,nin1,j)
        s1=zin(2,nin1,j)
        r2=zin(1,nin2,j)
        s2=zin(2,nin2,j)
        r3=zin(1,nin3,j)
        s3=zin(2,nin3,j)
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
3001        CONTINUE
#ifdef __SR8000 
!poption indep(zout)
!voption pvfunc(3)
#endif 
        DO 3000,ia=2,after
        ias=ia-1
        IF (4*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
#ifdef ___SR8000 
!poption indep(zout)
#endif 
                DO 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3010,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r2=-zin(2,nin2,j)
                s2=zin(1,nin2,j)
                r3=-zin(1,nin3,j)
                s3=-zin(2,nin3,j)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3010                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                DO 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3020,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r2=zin(2,nin2,j)
                s2=-zin(1,nin2,j)
                r3=-zin(1,nin3,j)
                s3=-zin(2,nin3,j)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3020                CONTINUE
        ENDIF
        ELSE IF (8*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                DO 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3030,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r=zin(1,nin2,j)
                s=zin(2,nin2,j)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=-zin(2,nin3,j)
                s3=zin(1,nin3,j)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3030                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                DO 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3040,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r=zin(1,nin2,j)
                s=zin(2,nin2,j)
                r2=(r + s)*rt2i
                s2=(-r + s)*rt2i
                r3=zin(2,nin3,j)
                s3=-zin(1,nin3,j)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3040                CONTINUE
        ENDIF
        ELSE
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
        DO 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3090,j=1,nfft
        r1=zin(1,nin1,j)
        s1=zin(2,nin1,j)
        r=zin(1,nin2,j)
        s=zin(2,nin2,j)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,nin3,j)
        s=zin(2,nin3,j)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
3090        CONTINUE
        ENDIF
3000        CONTINUE
        ELSE IF (now.eq.5) THEN
! cos(2._real_8*pi/5._real_8)
        cos2=0.3090169943749474_real_8
! cos(4._real_8*pi/5._real_8)
        cos4=-0.8090169943749474_real_8
! sin(2._real_8*pi/5._real_8)
        sin2=isign*0.9510565162951536_real_8
! sin(4._real_8*pi/5._real_8)
        sin4=isign*0.5877852522924731_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
        DO 5001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        DO 5001,j=1,nfft
        r1=zin(1,nin1,j)
        s1=zin(2,nin1,j)
        r2=zin(1,nin2,j)
        s2=zin(2,nin2,j)
        r3=zin(1,nin3,j)
        s3=zin(2,nin3,j)
        r4=zin(1,nin4,j)
        s4=zin(2,nin4,j)
        r5=zin(1,nin5,j)
        s5=zin(2,nin5,j)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,j,nout1) = r1 + r25 + r34
        r = cos2*r25 + cos4*r34 + r1
        s = sin2*s25 + sin4*s34
        zout(1,j,nout2) = r - s
        zout(1,j,nout5) = r + s
        r = cos4*r25 + cos2*r34 + r1
        s = sin4*s25 - sin2*s34
        zout(1,j,nout3) = r - s
        zout(1,j,nout4) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,j,nout1) = s1 + s25 + s34
        r = cos2*s25 + cos4*s34 + s1
        s = sin2*r25 + sin4*r34
        zout(2,j,nout2) = r + s
        zout(2,j,nout5) = r - s
        r = cos4*s25 + cos2*s34 + s1
        s = sin4*r25 - sin2*r34
        zout(2,j,nout3) = r + s
        zout(2,j,nout4) = r - s
5001        CONTINUE
#ifdef __SR8000 
!poption indep(zout)
!voption pvfunc(3)
#endif 
        DO 5000,ia=2,after
        ias=ia-1
        IF (8*ias.eq.5*after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 5010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5010,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=-zin(2,nin3,j)
                        s3=zin(1,nin3,j)
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=-zin(1,nin5,j)
                        s5=-zin(2,nin5,j)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
5010                        CONTINUE
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                        DO 5020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5020,j=1,nfft
                        r1=zin(1,nin1,j)
                        s1=zin(2,nin1,j)
                        r=zin(1,nin2,j)
                        s=zin(2,nin2,j)
                        r2=(r + s)*rt2i
                        s2=(-r + s)*rt2i
                        r3=zin(2,nin3,j)
                        s3=-zin(1,nin3,j)
                        r=zin(1,nin4,j)
                        s=zin(2,nin4,j)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r5=-zin(1,nin5,j)
                        s5=-zin(2,nin5,j)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
5020                        CONTINUE
                ENDIF
        ELSE
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
                DO 5100,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                DO 5100,j=1,nfft
                r1=zin(1,nin1,j)
                s1=zin(2,nin1,j)
                r=zin(1,nin2,j)
                s=zin(2,nin2,j)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,nin3,j)
                s=zin(2,nin3,j)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,nin4,j)
                s=zin(2,nin4,j)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,nin5,j)
                s=zin(2,nin5,j)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,j,nout1) = r1 + r25 + r34
                r = cos2*r25 + cos4*r34 + r1
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = cos4*r25 + cos2*r34 + r1
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 + s34
                r = cos2*s25 + cos4*s34 + s1
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = cos4*s25 + cos2*s34 + s1
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
5100                CONTINUE
        ENDIF
5000        CONTINUE
       ELSE IF (now.eq.6) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8

        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout)
#endif 
        DO 6120,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        DO 6120,j=1,nfft
        r2=zin(1,nin3,j)
        s2=zin(2,nin3,j)
        r3=zin(1,nin5,j)
        s3=zin(2,nin5,j)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,nin1,j)
        s1=zin(2,nin1,j)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,nin6,j)
        s2=zin(2,nin6,j)
        r3=zin(1,nin2,j)
        s3=zin(2,nin2,j)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,nin4,j)
        s1=zin(2,nin4,j)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,j,nout1)=ur1+vr1
        zout(2,j,nout1)=ui1+vi1
        zout(1,j,nout5)=ur2+vr2
        zout(2,j,nout5)=ui2+vi2
        zout(1,j,nout3)=ur3+vr3
        zout(2,j,nout3)=ui3+vi3
        zout(1,j,nout4)=ur1-vr1
        zout(2,j,nout4)=ui1-vi1
        zout(1,j,nout2)=ur2-vr2
        zout(2,j,nout2)=ui2-vi2
        zout(1,j,nout6)=ur3-vr3
        zout(2,j,nout6)=ui3-vi3
6120        CONTINUE

        ELSE
          CALL stopgm('fftpre','error',& 
 __LINE__,__FILE__)
        ENDIF

        RETURN
        END SUBROUTINE fftpre

! ==================================================================
      SUBROUTINE fftstp(mm,nfft,m,nn,n,zin,zout&
                       ,trig,now,after,before,isign)
! ==--------------------------------------------------------------==
! Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
! modified by Stefan Goedecker, Stuttgart, Germany, October 15, 1995
! Commercial use is prohibited without the explicit permission of 
! the author.
! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: mm, nfft, m, nn, n
  REAL(real_8)                               :: zin(2,mm,m), zout(2,nn,n), &
                                                trig(2,1024)
  INTEGER                                    :: now, after, before, isign

  INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
  REAL(real_8) :: am, ap, bb, bm, bp, ci2, ci3, ci4, ci5, cm, cos2, cos4, cp, &
      cr2, cr3, cr4, cr5, dm, dp, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, &
      r8, rt2i, s, s1, s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, &
      ui2, ui3, ur1, ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

      atn=after*now
      atb=after*before

! sqrt(0.5_real_8)
      rt2i=0.7071067811865475_real_8
      IF (now.eq.4) THEN
      IF (isign.eq.1) THEN
         ia=1
         nin1=ia-after
         nout1=ia-atn
#ifdef __SR8000
!poption indep(zout)
!voption pvfunc(3)
#endif 
                DO 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4001,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r - s
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r + s
                zout(2,j,nout4) = r - s
4001                CONTINUE
#ifdef __SR8000
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
                DO 4000,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=-zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
4010                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
4020                        CONTINUE
                ENDIF
4000                CONTINUE
        ELSE
                ia=1
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000
!poption indep(zout) 
#endif 
                DO 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4101,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
4101                CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
                DO 4100,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4110,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=-zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
4110                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
4120                        CONTINUE
                ENDIF
4100                CONTINUE
        ENDIF
        ELSE IF (now.eq.8) THEN
        IF (isign.eq.-1) THEN
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 8120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = (-cp + dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp = ( cm - dp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dp
8120                        CONTINUE
        ELSE
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 8121,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8121,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= (-cm + dp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dp
8121                        CONTINUE
        ENDIF
        ELSE IF (now.eq.3) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
3001        CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
        DO 3000,ia=2,after
        ias=ia-1
        IF (4*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
                DO 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3010,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=-zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=-zin(1,j,nin3)
                s3=-zin(2,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3010                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
                DO 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3020,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=-zin(1,j,nin2)
                r3=-zin(1,j,nin3)
                s3=-zin(2,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3020                CONTINUE
        ENDIF
        ELSE IF (8*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3030,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=-zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3030                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3040,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(-r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=-zin(1,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
3040                CONTINUE
        ENDIF
        ELSE
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3090,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r=zin(1,j,nin2)
        s=zin(2,j,nin2)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,j,nin3)
        s=zin(2,j,nin3)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
3090        CONTINUE
        ENDIF
3000        CONTINUE
        ELSE IF (now.eq.5) THEN
! cos(2._real_8*pi/5._real_8)
        cos2=0.3090169943749474_real_8
! cos(4._real_8*pi/5._real_8)
        cos4=-0.8090169943749474_real_8
! sin(2._real_8*pi/5._real_8)
        sin2=isign*0.9510565162951536_real_8
! sin(4._real_8*pi/5._real_8)
        sin4=isign*0.5877852522924731_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 5001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        DO 5001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r4=zin(1,j,nin4)
        s4=zin(2,j,nin4)
        r5=zin(1,j,nin5)
        s5=zin(2,j,nin5)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,j,nout1) = r1 + r25 + r34
        r = cos2*r25 + cos4*r34 + r1
        s = sin2*s25 + sin4*s34
        zout(1,j,nout2) = r - s
        zout(1,j,nout5) = r + s
        r = cos4*r25 + cos2*r34 + r1
        s = sin4*s25 - sin2*s34
        zout(1,j,nout3) = r - s
        zout(1,j,nout4) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,j,nout1) = s1 + s25 + s34
        r = cos2*s25 + cos4*s34 + s1
        s = sin2*r25 + sin4*r34
        zout(2,j,nout2) = r + s
        zout(2,j,nout5) = r - s
        r = cos4*s25 + cos2*s34 + s1
        s = sin4*r25 - sin2*r34
        zout(2,j,nout3) = r + s
        zout(2,j,nout4) = r - s
5001        CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
        DO 5000,ia=2,after
        ias=ia-1
        IF (8*ias.eq.5*after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 5010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=-zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=-zin(1,j,nin5)
                        s5=-zin(2,j,nin5)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
5010                        CONTINUE
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
                        DO 5020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(-r + s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=-zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r5=-zin(1,j,nin5)
                        s5=-zin(2,j,nin5)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
5020                        CONTINUE
                ENDIF
        ELSE
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 5100,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                DO 5100,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,j,nin5)
                s=zin(2,j,nin5)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,j,nout1) = r1 + r25 + r34
                r = cos2*r25 + cos4*r34 + r1
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = cos4*r25 + cos2*r34 + r1
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 + s34
                r = cos2*s25 + cos4*s34 + s1
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = cos4*s25 + cos2*s34 + s1
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
5100                CONTINUE
        ENDIF
5000        CONTINUE
       ELSE IF (now.eq.6) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8

        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000
!poption indep(zout) 
#endif 
        DO 6120,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        DO 6120,j=1,nfft
        r2=zin(1,j,nin3)
        s2=zin(2,j,nin3)
        r3=zin(1,j,nin5)
        s3=zin(2,j,nin5)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,j,nin6)
        s2=zin(2,j,nin6)
        r3=zin(1,j,nin2)
        s3=zin(2,j,nin2)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin4)
        s1=zin(2,j,nin4)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,j,nout1)=ur1+vr1
        zout(2,j,nout1)=ui1+vi1
        zout(1,j,nout5)=ur2+vr2
        zout(2,j,nout5)=ui2+vi2
        zout(1,j,nout3)=ur3+vr3
        zout(2,j,nout3)=ui3+vi3
        zout(1,j,nout4)=ur1-vr1
        zout(2,j,nout4)=ui1-vi1
        zout(1,j,nout2)=ur2-vr2
        zout(2,j,nout2)=ui2-vi2
        zout(1,j,nout6)=ur3-vr3
        zout(2,j,nout6)=ui3-vi3
6120        CONTINUE

        ELSE
          CALL stopgm('FFTSTP','error',& 
 __LINE__,__FILE__)
        ENDIF

        RETURN
     END SUBROUTINE fftstp

     SUBROUTINE fftrot(mm,nfft,m,nn,n,zin,zout,&
                       trig,now,after,before,isign)
! Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
! modified by Stefan Goedecker, Stuttgart, Germany, October 15, 1995
! Commercial use is prohibited without the explicit permission of the author.
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: mm, nfft, m, nn, n
  REAL(real_8)                               :: zin(2,mm,m), zout(2,n,nn), &
                                                trig(2,1024)
  INTEGER                                    :: now, after, before, isign

  INTEGER :: atb, atn, ia, ias, ib, itrig, itt, j, nin1, nin2, nin3, nin4, &
      nin5, nin6, nin7, nin8, nout1, nout2, nout3, nout4, nout5, nout6, &
      nout7, nout8
  REAL(real_8) :: am, ap, bb, bm, bp, ci2, ci3, ci4, ci5, cm, cos2, cos4, cp, &
      cr2, cr3, cr4, cr5, dm, dp, r, r1, r2, r25, r3, r34, r4, r5, r6, r7, &
      r8, rt2i, s, s1, s2, s25, s3, s34, s4, s5, s6, s7, s8, sin2, sin4, ui1, &
      ui2, ui3, ur1, ur2, ur3, vi1, vi2, vi3, vr1, vr2, vr3

        atn=after*now
        atb=after*before

! sqrt(0.5_real_8)
        rt2i=0.7071067811865475_real_8
        IF (now.eq.4) THEN
        IF (isign.eq.1) THEN
                ia=1
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4001,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r - s
                zout(1,nout4,j) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r + s
                zout(2,nout4,j) = r - s
4001                CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
                DO 4000,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=-zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r + s
                        zout(2,nout4,j) = r - s
4010                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r + s
                        zout(2,nout4,j) = r - s
4020                        CONTINUE
                ENDIF
4000                CONTINUE
        ELSE
                ia=1
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO 4101,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r + s
                zout(1,nout4,j) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r - s
                zout(2,nout4,j) = r + s
4101                CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
                DO 4100,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4110,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=-zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4110                        CONTINUE
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO 4120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4120                        CONTINUE
                ENDIF
4100                CONTINUE
        ENDIF

        ELSE IF (now.eq.8) THEN
        IF (isign.eq.-1) THEN
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 8120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am + dm
                        zout(2,nout3,j) = cm - bm
                        zout(1,nout7,j) = am - dm
                        zout(2,nout7,j) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = (-cp + dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp = ( cm - dp)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp
8120                        CONTINUE
        ELSE
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 8121,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO 8121,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am - dm
                        zout(2,nout3,j) = cm + bm
                        zout(1,nout7,j) = am + dm
                        zout(2,nout7,j) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= (-cm + dp)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp
8121                        CONTINUE
        ENDIF
        ELSE IF (now.eq.3) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2
        zout(2,nout3,j) = s1 - r2
3001        CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
        DO 3000,ia=2,after
        ias=ia-1
        IF (4*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3010,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=-zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=-zin(1,j,nin3)
                s3=-zin(2,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2
                zout(2,nout3,j) = s1 - r2
3010                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3020,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=-zin(1,j,nin2)
                r3=-zin(1,j,nin3)
                s3=-zin(2,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2
                zout(2,nout3,j) = s1 - r2
3020                CONTINUE
        ENDIF
        ELSE IF (8*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3030,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=-zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2
                zout(2,nout3,j) = s1 - r2
3030                CONTINUE
        ELSE
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO 3040,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(-r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=-zin(1,j,nin3)
                r=r2 + r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - 0.5_real_8*r
                s1=s1 - 0.5_real_8*s
                r2=bb*(r2-r3)
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2
                zout(2,nout3,j) = s1 - r2
3040                CONTINUE
        ENDIF
        ELSE
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO 3090,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r=zin(1,j,nin2)
        s=zin(2,j,nin2)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,j,nin3)
        s=zin(2,j,nin3)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2
        zout(2,nout3,j) = s1 - r2
3090        CONTINUE
        ENDIF
3000        CONTINUE
        ELSE IF (now.eq.5) THEN
! cos(2._real_8*pi/5._real_8)
        cos2=0.3090169943749474_real_8
! cos(4._real_8*pi/5._real_8)
        cos4=-0.8090169943749474_real_8
! sin(2._real_8*pi/5._real_8)
        sin2=isign*0.9510565162951536_real_8
! sin(4._real_8*pi/5._real_8)
        sin4=isign*0.5877852522924731_real_8
        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 5001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        DO 5001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r4=zin(1,j,nin4)
        s4=zin(2,j,nin4)
        r5=zin(1,j,nin5)
        s5=zin(2,j,nin5)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,nout1,j) = r1 + r25 + r34
        r = cos2*r25 + cos4*r34 + r1
        s = sin2*s25 + sin4*s34
        zout(1,nout2,j) = r - s
        zout(1,nout5,j) = r + s
        r = cos4*r25 + cos2*r34 + r1
        s = sin4*s25 - sin2*s34
        zout(1,nout3,j) = r - s
        zout(1,nout4,j) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,nout1,j) = s1 + s25 + s34
        r = cos2*s25 + cos4*s34 + s1
        s = sin2*r25 + sin4*r34
        zout(2,nout2,j) = r + s
        zout(2,nout5,j) = r - s
        r = cos4*s25 + cos2*s34 + s1
        s = sin4*r25 - sin2*r34
        zout(2,nout3,j) = r + s
        zout(2,nout4,j) = r - s
5001        CONTINUE
#ifdef __SR8000 
!poption indep(zout) 
!voption pvfunc(3) 
#endif 
        DO 5000,ia=2,after
        ias=ia-1
        IF (8*ias.eq.5*after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO 5010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=-zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=-(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=-zin(1,j,nin5)
                        s5=-zin(2,j,nin5)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,nout1,j) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5010                        CONTINUE
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                        DO 5020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO 5020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(-r + s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=-zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=-(r + s)*rt2i
                        r5=-zin(1,j,nin5)
                        s5=-zin(2,j,nin5)
                        r25 = r2 + r5
                        r34 = r3 + r4
                        s25 = s2 - s5
                        s34 = s3 - s4
                        zout(1,nout1,j) = r1 + r25 + r34
                        r = cos2*r25 + cos4*r34 + r1
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = cos4*r25 + cos2*r34 + r1
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 - r5
                        r34 = r3 - r4
                        s25 = s2 + s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 + s34
                        r = cos2*s25 + cos4*s34 + s1
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = cos4*s25 + cos2*s34 + s1
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5020                        CONTINUE
                ENDIF
        ELSE
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
                DO 5100,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                DO 5100,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,j,nin5)
                s=zin(2,j,nin5)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,nout1,j) = r1 + r25 + r34
                r = cos2*r25 + cos4*r34 + r1
                s = sin2*s25 + sin4*s34
                zout(1,nout2,j) = r - s
                zout(1,nout5,j) = r + s
                r = cos4*r25 + cos2*r34 + r1
                s = sin4*s25 - sin2*s34
                zout(1,nout3,j) = r - s
                zout(1,nout4,j) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,nout1,j) = s1 + s25 + s34
                r = cos2*s25 + cos4*s34 + s1
                s = sin2*r25 + sin4*r34
                zout(2,nout2,j) = r + s
                zout(2,nout5,j) = r - s
                r = cos4*s25 + cos2*s34 + s1
                s = sin4*r25 - sin2*r34
                zout(2,nout3,j) = r + s
                zout(2,nout4,j) = r - s
5100                CONTINUE
        ENDIF
5000        CONTINUE
       ELSE IF (now.eq.6) THEN
! 0.5_real_8*sqrt(3._real_8)
        bb=isign*0.8660254037844387_real_8

        ia=1
        nin1=ia-after
        nout1=ia-atn
#ifdef __SR8000 
!poption indep(zout) 
#endif 
        DO 6120,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        DO 6120,j=1,nfft
        r2=zin(1,j,nin3)
        s2=zin(2,j,nin3)
        r3=zin(1,j,nin5)
        s3=zin(2,j,nin5)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,j,nin6)
        s2=zin(2,j,nin6)
        r3=zin(1,j,nin2)
        s3=zin(2,j,nin2)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin4)
        s1=zin(2,j,nin4)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - 0.5_real_8*r
        s1=s1 - 0.5_real_8*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,nout1,j)=ur1+vr1
        zout(2,nout1,j)=ui1+vi1
        zout(1,nout5,j)=ur2+vr2
        zout(2,nout5,j)=ui2+vi2
        zout(1,nout3,j)=ur3+vr3
        zout(2,nout3,j)=ui3+vi3
        zout(1,nout4,j)=ur1-vr1
        zout(2,nout4,j)=ui1-vi1
        zout(1,nout2,j)=ur2-vr2
        zout(2,nout2,j)=ui2-vi2
        zout(1,nout6,j)=ur3-vr3
        zout(2,nout6,j)=ui3-vi3
6120        CONTINUE

       ELSE
        CALL stopgm('FFTROT','error',& 
 __LINE__,__FILE__)
       ENDIF

        RETURN
     END SUBROUTINE fftrot
