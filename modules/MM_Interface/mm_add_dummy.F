C     ==================================================================
       SUBROUTINE mm_add_dummy(NRPTp,NSOLVp)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
! GROMOSDIP
       use mm_input
       use coordar
       implicit none

       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h'
       !include 'Gromos/coordsz.h'
       !include 'Gromos/coordar.h'
       include 'Gromos/md.h'
       include 'Gromos/forcesz.h'
       include 'Gromos/forcear.h'
       INTEGER NRPT,NSOLV,NATTOT,NUM4DA,NRP_old
       INTEGER IUTRJX,IUTRJV,IUTRJE,IUTRJG,IFNCFG
       common/arg_i/NRPT,NSOLV,NATTOT,NUM4DA,IUTRJX,IUTRJV,IUTRJE,
     c             IUTRJG,IFNCFG 
!
       integer     NRPTp,NSOLVp
!
       integer mr,ic,iac_dummy,i,ipp,j,ii
       character CHARO

       addh%igdummy=NRP+1
C (GM) determine non vdW interacting atom type
C note that all added hydrogens show up in the
C vdW+electrostatic interaction lists (explicit exclusions are not constructed)
       do i=1,NRATT
         ic=i*(i+1)/2
         if((DABS(C12(ic)).lt.1.d-24).and.(DABS(C6(ic)).lt.1.d-24))
     $          goto 111
       enddo
       call stopgm('MM_ADD_DUMMY','NO SUITABLE VDW TYPE FOR DUMMY',
     & __LINE__,"__file__")
111    continue
       iac_dummy=i
       NRAA2=NRAA2+addh%n_added_h
!$OMP parallel do private(i)
       DO i=NRAA2-addh%n_added_h+1,NRAA2
         AANM(i)="DUM"
       ENDDO
       mr=MRES(NRP)
       ii=0
       do i=NRP+1,NRP+addh%n_added_h
         ii=ii+1
         gqmmm%gr_atom(igqmmm%sp_H_added,ii)=i
         mr=mr+1
         MRES(i)=mr
         PANM(i)="DUM"
         IAC(i)=iac_dummy
         CG(i)=0.d0
         INE(i)=0
         INE14(i)=0
         WMAS(i)= 1.008d0
         WINV(i)= 1.D0 / WMAS(i) 
       enddo
       ic=NRP
       do i=NCAG+1,NCAG+addh%n_added_h
         ic=ic+1
         INC(i)=ic
       enddo
       NCAG=NCAG+addh%n_added_h
       NRP=NRP+addh%n_added_h
       NATTOT=NATTOT+addh%n_added_h

       nsp_t%NSP(nsp_t%NSPM+1)=NRP
       nsp_t%NSPM=nsp_t%NSPM+1

       CHARO = CHAR(1)
!$OMP parallel do private(I)
       DO I=1,coordsz_t%MAXNAT
          IAGRP(I) = CHARO
       ENDDO
       NUSNRE = 0
       NUNRE2 = 1

!$OMP parallel do private(I,J)
       DO I=1,MAXNRE
         DO J=I,MAXNRE
           NRELKP(I,J) = 1
           NRELKP(J,I) = 1
         ENDDO
       ENDDO

       NRPT = NRP*NPM
       NSOLV= NRAM*NSM
       NRP_old=NRP-addh%n_added_h
       do i=1,3*NSOLV
         ipp=3*(NRP_old+NSOLV)-i+1
         x(ipp+addh%n_added_h*3)=x(ipp)
         v(ipp+addh%n_added_h*3)=v(ipp)
       enddo
!$OMP parallel do private(i)
       DO i=3*(NRP_old)+1,3*NRP
         x(i)=DBLE(i)*1.d-6
         v(i)=0.0D0
       ENDDO
       NRPTp=NRPT
       NSOLVp=NSOLV

       RETURN
       END
