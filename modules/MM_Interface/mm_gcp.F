C      -*- fortran -*-
C  note the following stuff is to revert in case this is really needed.
C if you compile with -DMM_GCP_NOSHIFT the coordinates are not shifted
C when converting between gromos and cpmd.
C
C HOWEVER, when looking at the CRD_INI.grm and CRD_FIN.grm files, they
C have often (small) parts of the structure wrapped back into the box,
C but mirrored.
C this happens most frequently, when importing structures from AMBER.
C note, that the initial gromos file created with amber2gromos is ok,
C but the ones from CPMD are not, and that excludes them from being useful.
C when shifting all atoms in the box by half a box diameter, this effect
C seems to go away without really changing the forces and energies
C (at least not for me), so this seems much more useful to me.
C 
C please let me know, if you find a case, where this does not work.
C axel kohlmeyer <akohlmey@cmm.upenn.edu> 2005/10/26.
C
#if defined(MM_GCP_NOSHIFT)
C     ==================================================================
       SUBROUTINE mm_gcp(nax,nsx,epot_mm_gr,epot_mm_cp,
     c                   g_to_cp,convtype,xgr,xcp)
C     ==--------------------------------------------------------------==
! GROMOSDIP
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
       use mm_dimmod
       use cnst
       use coordsz
       implicit none

       !include 'Gromos/coordsz.h'
       integer   nax,nsx
       REAL*8    xcp(3,nax,*),xgr(*)
       real*8    epot_mm_gr,epot_mm_cp
       character g_to_cp,convtype
c
c locals
c
       integer i,NRTOT
       real*8 fconv
c
       NRTOT=solsolv%NRPT+solsolv%NSOLV

       fconv=0.D0
       if (g_to_cp.eq.'C')then  !convert to cp array
         if(convtype.eq.'P')fconv=10.d0*FBOHR
         if(convtype.eq.'F')fconv=0.1d0/FBOHR/AU_KJM
         DO i=1,NRTOT
c          write(6,*)i,cpat(i),cpsp(i)
            xcp(1,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+1)
            xcp(2,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+2)
            xcp(3,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+3)
         ENDDO
         epot_mm_cp=epot_mm_gr/AU_KJM
       else if (g_to_cp.eq.'G') then  !convert to gromos array
         if(convtype.eq.'P') fconv=0.1d0/FBOHR
         if(convtype.eq.'F') fconv=AU_KJM*10.d0*FBOHR
         DO i=1,NRTOT
            xgr(3*(i-1)+1) =fconv*xcp(1,cpat(i),cpsp(i))
            xgr(3*(i-1)+2) =fconv*xcp(2,cpat(i),cpsp(i))
            xgr(3*(i-1)+3) =fconv*xcp(3,cpat(i),cpsp(i))
         ENDDO
c        epot_mm_gr=epot_mm_cp*ht_to_kj
       else
         call stopgm('MM_GCP','UNKNOWN CONVERSION FLAG',
     & __LINE__,"__file__")
       endif
       
       RETURN
       END
#else
C     ==================================================================
       SUBROUTINE mm_gcp(nax,nsx,epot_mm_gr,epot_mm_cp,
     c                   g_to_cp,convtype,xgr,xcp)
C     ==--------------------------------------------------------------==
! GROMOSDIP
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
       use mm_dimmod
       use cnst
       use coordsz
       implicit none

       !include 'Gromos/coordsz.h'
       integer   nax,nsx
       REAL*8    xcp(3,nax,nsx),xgr(3*coordsz_t%MAXNAT)
       real*8    epot_mm_gr,epot_mm_cp
       character g_to_cp,convtype
c
c locals
c
       integer i,NRTOT
       real*8 fconv,boxby2(3)
c
       NRTOT=solsolv%NRPT+solsolv%NSOLV

       fconv=0.D0
       if (g_to_cp.eq.'C')then  !convert to cp array
         if(convtype.eq.'P') then 
           fconv=10.d0*FBOHR
           boxby2(1)=clsaabox%BOX_au(1)*0.5d0
           boxby2(2)=clsaabox%BOX_au(2)*0.5d0
           boxby2(3)=clsaabox%BOX_au(3)*0.5d0
           DO i=1,NRTOT
c             write(6,*)i,cpat(i),cpsp(i)
             xcp(1,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+1)-boxby2(1)
             xcp(2,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+2)-boxby2(2)
             xcp(3,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+3)-boxby2(3)
           ENDDO
         elseif(convtype.eq.'F') then
           fconv=0.1d0/FBOHR/AU_KJM
           DO i=1,NRTOT
c             write(6,*)i,cpat(i),cpsp(i)
             xcp(1,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+1)
             xcp(2,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+2)
             xcp(3,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+3)
           ENDDO
         else
           call stopgm('MM_GCP','UNKNOWN CONVERSION TYPE',
     & __LINE__,"__file__")
         endif
         epot_mm_cp=epot_mm_gr/AU_KJM
       else if (g_to_cp.eq.'G') then  !convert to gromos array
         if(convtype.eq.'P') then
           fconv=0.1d0/FBOHR
           boxby2(1)=clsaabox%BOX_au(1)*0.5d0*fconv
           boxby2(2)=clsaabox%BOX_au(2)*0.5d0*fconv
           boxby2(3)=clsaabox%BOX_au(3)*0.5d0*fconv
           DO i=1,NRTOT
             xgr(3*(i-1)+1) =fconv*xcp(1,cpat(i),cpsp(i))+boxby2(1)
             xgr(3*(i-1)+2) =fconv*xcp(2,cpat(i),cpsp(i))+boxby2(2)
             xgr(3*(i-1)+3) =fconv*xcp(3,cpat(i),cpsp(i))+boxby2(3)
           ENDDO
         elseif(convtype.eq.'F') then
           fconv=AU_KJM*10.d0*FBOHR
           DO i=1,NRTOT
             xgr(3*(i-1)+1) =fconv*xcp(1,cpat(i),cpsp(i))
             xgr(3*(i-1)+2) =fconv*xcp(2,cpat(i),cpsp(i))
             xgr(3*(i-1)+3) =fconv*xcp(3,cpat(i),cpsp(i))
           ENDDO
         else
           call stopgm('MM_GCP','UNKNOWN CONVERSION TYPE',
     & __LINE__,"__file__")
         endif
c        epot_mm_gr=epot_mm_cp*ht_to_kj
       else
         call stopgm('MM_GCP','UNKNOWN CONVERSION DIRECTION',
     & __LINE__,"__file__")
       endif
       
       RETURN
       END
#endif
C     ==================================================================
       SUBROUTINE mm_gcp_vel(nax,nsx,g_to_cp,convtype,xgr,xcp)
C     ==--------------------------------------------------------------==
! GROMOSDIP
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
       use cnst
       use mm_dimmod
       use coordsz
       implicit none

       !include 'Gromos/coordsz.h'
       integer   nax,nsx
       REAL*8    xcp(3,nax,nsx),xgr(3*coordsz_t%MAXNAT)
       character g_to_cp,convtype
c
c locals
c
       integer i,NRTOT
       real*8 fconv

       if ( convtype.ne.'V' ) then
          call stopgm('MM_GCP_VEL','CONVERSION FLAG MUST BE: V',
     & __LINE__,"__file__")
       endif

       NRTOT=solsolv%NRPT+solsolv%NSOLV

       fconv=0.D0
       if (g_to_cp.eq.'C')then  !convert to cp array
         write (6,*) ' '
         write (6,*) 'LOADING GROMOS VELOCITIES.'
         write (6,*) ' '
         if(convtype.eq.'V')fconv=1.0d-2*FBOHR*AU_FS
c        fconv=1.0D0
         DO i=1,NRTOT
c          write(6,*)i,cpat(i),cpsp(i)
            xcp(1,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+1)
            xcp(2,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+2)
            xcp(3,cpat(i),cpsp(i))=fconv*xgr(3*(i-1)+3)
         ENDDO
       else if (g_to_cp.eq.'G') then  !convert to gromos array
         if(convtype.eq.'V')fconv=1.0d2/FBOHR/AU_FS
         DO i=1,NRTOT
            xgr(3*(i-1)+1) =fconv*xcp(1,cpat(i),cpsp(i))
            xgr(3*(i-1)+2) =fconv*xcp(2,cpat(i),cpsp(i))
            xgr(3*(i-1)+3) =fconv*xcp(3,cpat(i),cpsp(i))
         ENDDO
       else
         call stopgm('MM_GCP_VEL','UNKNOWN CONVERSION FLAG',
     & __LINE__,"__file__")
       endif
       
       RETURN
       END
