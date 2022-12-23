! GROMOSDIP
C     ==================================================================
       SUBROUTINE mm_readgromos(nax,nsx,TAU)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 23 Feb. 2005
       use coordar
       implicit none

       integer nax,nsx
       REAL*8     TAU(3,nax,nsx)
       !include 'Gromos/coordsz.h'
       !include 'Gromos/coordar.h'
c
c locals
c
       REAL*8 dum1,dum2
       dum1=0.D0
       dum2=0.D0
       CALL mm_gcp(nax,nsx,dum1,dum2,'C','P',x,TAU)

       RETURN
       END
C     ==================================================================
       SUBROUTINE mm_readgromos_vel(nax,nsx,VELP)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==

       use coordar
       implicit none

       integer i

       integer nax,nsx
       REAL*8     VELP(*)
       !include 'Gromos/coordsz.h'
       !include 'Gromos/coordar.h'

       CALL mm_gcp_vel(nax,nsx,'C','V',v,VELP)
       RETURN
       END

