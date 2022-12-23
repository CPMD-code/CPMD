! GROMOSDIP
C     ==================================================================
      SUBROUTINE mm_write_gromos_coord(filename,TAU,VELP,nax,nsx)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 23 Feb. 2005
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen
      use machine, only: m_datum
      use coordar
      implicit none
      integer nax,nsx
      real*8  TAU(*),VELP(*)
      character*11 filename
      CHARACTER DATX*26

      include 'Gromos/md.h'
      !include 'Gromos/coordsz.h'
      !include 'Gromos/coordar.h'
      include 'Gromos/toposz.h'
      include 'Gromos/topoar.h'
      include 'Gromos/box.h'

      integer NAT
      real*8 dum1,dum2
      LOGICAL FERROR

      CALL mm_gcp(nax,nsx,dum1,dum2,'G','P',x,tau)
      CALL mm_gcp_vel(nax,nsx,'G','V',v,VELP)

      CALL FILEOPEN(512,filename,FO_DEF,FERROR)
      WRITE(512,'(a5)')'TITLE'
      CALL m_datum(datx)
      WRITE(512,*)datx
      WRITE(512,'(a3)')'END'
      NAT=NRP+NSM*NRAM
      CALL  WRXPOS(512,.true.,.false.,NAT,1,NSM,3,X) 
      CALL  WRXVEL(512,.true.,.false.,NAT,1,NSM,3,V) 
      CALL  WRBOX(512,.true.,BOX) 
      CALL FILECLOSE(512)

      RETURN
      END
     
