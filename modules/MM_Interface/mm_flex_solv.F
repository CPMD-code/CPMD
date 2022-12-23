C     ==================================================================
      SUBROUTINE mm_flex_solv(NRPT,NSOLV)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
! GROMOSDIP
cmb - Revised on 21 June 2006
      use system
      use mm_input
      use store_types
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen
      use coordar
      implicit none
      !include 'Gromos/coordsz.h'
      !include 'Gromos/coordar.h'
      include 'Gromos/toposz.h'
      include 'Gromos/topoar.h'
      include 'Gromos/md.h'
      integer NRPT,NSOLV
!
      integer i,ig0,ires,ic,ige,ine0,is,ns,NSM_fin,iw,ig0s
      logical dostop,ferror
      logical, allocatable ::  to_solute(:)
      real*8 xs
!
      IF(restart1%restart) call stopgm('mm_flex_solv',
     c  'NO RESTART ALLOWED IF THE SOLVENT IS TURNED INTO SOLUTE',
     & __LINE__,"__file__")
      kne(1)=0
      DO i=2,NRP
        kne(i)=kne(i-1)+INE(i-1)
      ENDDO
      ig0=NRP
      ires=NRAA2
!
      ALLOCATE(to_solute(NSM))
      IF(solqmmm%all_water)THEN
!$OMP parallel do private(ns)
#ifdef __SR11000
*poption parallel, tlocal(ns)
#endif
        DO ns=1,NSM 
          to_solute(ns)=.true.
        ENDDO
      ELSE
        ig0s=NRP
        DO ns=1,NSM
          to_solute(ns)=.false.
!$OMP parallel do private(iw)
          DO iw=1,solqmmi%nwq
            if((ig0+1).eq.solqmmi%ig_wq(iw)) to_solute(ns)=.true.
          ENDDO
          IF(to_solute(ns)) THEN
            DO i=1,NRAM*3
              xs=x(3*ig0+i)
              x(3*ig0+i)=x(3*ig0s+i)
              x(3*ig0s+i)=xs
            ENDDO
            ig0s=ig0s+NRAM
          ENDIF
          ig0=ig0+NRAM
        ENDDO
      ENDIF
!
      NSM_fin=NSM
      ig0=NRP
      DO ns=1,NSM
        IF(to_solute(ns)) THEN
          NSM_fin=NSM_fin-1
          ires=ires+1
          NCAG=NCAG+1
          AANM(ires)='H2O'
! add bonds
          DO ic=1,NCONS
            NBONH=NBONH+1
            ICBH(NBONH)=solqmmi%solv_bond(ic)
            IBH(NBONH)=ig0+ICONS(ic)
            JBH(NBONH)=ig0+JCONS(ic)
          ENDDO
! add exclusions
          ine0=NRAM-1
          DO i=ig0+1,ig0+NRAM
            kne(i)=kne(i-1)+INE(i-1)
            INE14(i)=0
            INE(i)=ine0
            DO ige=1,INE(i)
              JSNE(kne(i)+ige)=i+ige
            ENDDO
            ine0=max(0,ine0-1)
          ENDDO
! charges
!$OMP parallel do private(i,is)
          DO i=ig0+1,ig0+NRAM
            MRES(i)=ires
            is=i-ig0
            CG(i)=CGS(is)
            IAC(i)=IACS(is)
            PANM(i)=ANMS(is)
            WMAS(i)=WMASS(is)
          ENDDO
          INC(NCAG)=ig0+NRAM
          ig0=ig0+NRAM
        ENDIF
      ENDDO
!
      NRP=ig0
      NRAA2=ires
      NRPT=NRP
      NSM=NSM_fin
      IF(solqmmm%all_water .and. (NSM.ne.0))
     &     call stopgm('MM_FLEX_SOLV','more water than expected',
     & __LINE__,"__file__")
      NSOLV=NRAM*NSM
!
      DEALLOCATE(to_solute)
      dostop=.false.
!
      IF(NRP.ge.MAXNRP)THEN
        write(6,*)'MAXNRP should be greater than',NRP
        dostop=.true.
      ENDIF
      IF(kne(ig0).ge.MAXAEX)THEN
        write(6,*)'MAXAEX should be greater than',kne(ig0)
        dostop=.true.
      ENDIF
      IF(ires.ge.MAXAA2)THEN
        write(6,*)'MAXAA2 should be greater than',ires
        dostop=.true.
      ENDIF
      IF(NBONH.ge.MAXBNH)THEN
        write(6,*)'MAXBNH should be greater than',NBONH
        dostop=.true.
      ENDIF
      IF(NCAG.ge.MAXCAG)THEN
        write(6,*)'MAXCAG should be greater than',NCAG
        dostop=.true.
      ENDIF
!
      IF(dostop)call stopgm('MM_FLEX_SOLV','DIMENSIONING PROBLEM',
     & __LINE__,"__file__")
      CALL FILEOPEN(1111,'MM_TOPOLOGY',FO_DEF,FERROR)
      CALL wrtopo(1111)
      CALL FILECLOSE(1111)
      RETURN
      END
