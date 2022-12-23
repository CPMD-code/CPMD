! GROMOSDIP
C     ==================================================================
       SUBROUTINE mm_add_hydrogen
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
cmb - Revised on 23 Feb. 2005
C     ==--------------------------------------------------------------==
       use mm_input
       use coordar
       implicit none

       !include 'Gromos/coordsz.h'
       !include 'Gromos/coordar.h'
       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h'
!
       character*3 hybr
       real*8      rr(3,5)
       integer     NBC,ibc(4),i,igc
       integer     NBTOT,NH_DUM,ic
       integer Nadd

!
! calculate the number of bonds of carbon igc
!
       Nadd = 0
       do ic=1,addh%n_carbon
         igc=addh%ig_carbon(ic)
         NBC=0
         do i=1,NBON
           if(ib(i).eq.igc .or. jb(i).eq.igc) then 
             NBC=NBC+1
             if(ib(i).eq.igc)ibc(NBC)=jb(i)
             if(jb(i).eq.igc)ibc(NBC)=ib(i)
           endif
         enddo
         do i=1,NBONH
           if(ibh(i).eq.igc .or. jbh(i).eq.igc) then 
             NBC=NBC+1
             if(ibh(i).eq.igc)ibc(NBC)=jbh(i)
             if(jbh(i).eq.igc)ibc(NBC)=ibh(i)
           endif
         enddo
!
! calculate the hybridization
!
         if(IAC(igc).eq.11) NH_DUM=0
         if(IAC(igc).eq.12) NH_DUM=1
         if(IAC(igc).eq.13) NH_DUM=2
         if(IAC(igc).eq.14) NH_DUM=3
         if(IAC(igc).eq.15) NH_DUM=4
         if(IAC(igc).eq.16) NH_DUM=1

         NBTOT=NBC + NH_DUM
         if(NBTOT .eq. 3) then
           hybr = 'sp2'
         else if(NBTOT .eq. 4) then
           hybr = 'sp3'
         else
           write(6,*)'wrong hybridization',NBTOT,NBC,NH_DUM
           call stopgm('MM_ADD_HYDROGEN','HYBRIDIZATION ERROR',
     & __LINE__,"__file__")
         endif
!
!  save the coordinetes
!
         rr(1,1)=x(1+(igc-1)*3)
         rr(2,1)=x(2+(igc-1)*3)
         rr(3,1)=x(3+(igc-1)*3)
!$OMP parallel do private(i)
         DO i=1,NBC
           rr(1,1+i)=x(1+(ibc(i)-1)*3)
           rr(2,1+i)=x(2+(ibc(i)-1)*3)
           rr(3,1+i)=x(3+(ibc(i)-1)*3)
         ENDDO

         call add_H(hybr,NBC,rr)
         WMAS(igc)= 12.011D0
         WINV(igc)= 1.D0/WMAS(igc)

         DO i=1,NH_DUM
           IF(iac(addh%igdummy).ne.19.and.iac(addh%igdummy).ne.18
     .                   .and.iac(addh%igdummy).ne.17) THEN
           call stopgm('mm_add_hydrogen','igdummy is not a dummy',
     & __LINE__,"__file__")
           ENDIF
           Nadd=Nadd+1
           x(1+(addh%igdummy-1)*3)=rr(1,1+NBC+i)
           x(2+(addh%igdummy-1)*3)=rr(2,1+NBC+i)
           x(3+(addh%igdummy-1)*3)=rr(3,1+NBC+i)
           WMAS(addh%igdummy)= 1.008D0
           WINV(addh%igdummy)= 1.D0/WMAS(igc)
           addh%igdummy = addh%igdummy+1
         ENDDO
       enddo

       RETURN
       END
C     ==================================================================
      SUBROUTINE add_H(hybrid,NBC,rr)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use constr_utils, only : funcr,funcd,funcdd,funct,funco,
     &     normalize,funcp,diffd,diffdd,diffr,difft,diffo,diffp,
     &     getscal,getnorm,normvec,vecmul,addvec,copyvec,zerovec,
     &     vecprod,vecrotz,vecroty,vecrotx,vecrotvec
C     ==--------------------------------------------------------------==
      IMPLICIT NONE

      character*3   hybrid
      integer       NBC
      real*8        rr(3,5)
      real*8        vec(3,10)  ! 10 dummy vectors
      real*8        hblength  ! our hydrogen bond length for all cases  
      real*8        pi,d1,d2

      pi=DACOS(-1.0D0)
C     hblength=1.1D0 ! angstrom
!     hblength=2.1D0 ! atomic units
      hblength=0.105D0  ! nm

      IF (hybrid.eq.'sp2') THEN
         IF (NBC.eq.1) THEN
            ! hard coded random vector
            vec(1,1)=0.91084371024D0
            vec(2,1)=0.23843383392D0
            vec(3,1)=0.17363649163D0
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,2),vec(1,2))
            call vecprod(vec(1,1),vec(1,2),vec(1,3))
            call normvec(vec(1,3),vec(1,3))
            call vecrotvec(vec(1,3),2*pi/3,vec(1,2),vec(1,4))
            call normvec(vec(1,4),vec(1,4))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,4),rr(1,3))
            call vecrotvec(vec(1,3),-2*pi/3,vec(1,2),vec(1,5))
            call normvec(vec(1,5),vec(1,5))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,5),rr(1,4))
         ELSE IF (NBC.eq.2) THEN  
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,2),vec(1,1))
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,3),vec(1,2))
            call vecrotvec(vec(1,1),pi,vec(1,2),vec(1,3))  
            call normvec(vec(1,3),vec(1,3))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,3),rr(1,4))
         ELSE
            write(6,*) "not able to put hydrogens back in " ,hybrid,NBC
            call stopgm('MM_ADD_HYDROGEN','NO SPACE TO ADD HYDROGEN',
     & __LINE__,"__file__")
         ENDIF
      ELSE IF (hybrid.eq.'sp3') THEN
         IF (NBC.eq.1) THEN
            ! hard coded random vector
            vec(1,1)=0.91084371024D0
            vec(2,1)=0.23843383392D0
            vec(3,1)=0.17363649163D0
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,2),vec(1,2))
            call vecprod(vec(1,1),vec(1,2),vec(1,3))
            call normvec(vec(1,3),vec(1,3))

            call vecrotvec(vec(1,3),1.9D0,vec(1,2),vec(1,4))
            call normvec(vec(1,4),vec(1,4))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,4),rr(1,3))

            call vecrotvec(vec(1,2),pi*2/3,vec(1,4),vec(1,5))
            call normvec(vec(1,5),vec(1,5))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,5),rr(1,4))

            call vecrotvec(vec(1,2),-pi*2/3,vec(1,4),vec(1,6))
            call normvec(vec(1,6),vec(1,6))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,6),rr(1,5))
            
         ELSE IF (NBC.eq.2) THEN
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,2),vec(1,1))
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,3),vec(1,2))

            call vecrotvec(vec(1,1),pi*2/3,vec(1,2),vec(1,3))
            call normvec(vec(1,3),vec(1,3))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,3),rr(1,4))

            call vecrotvec(vec(1,1),-pi*2/3,vec(1,2),vec(1,4))
            call normvec(vec(1,4),vec(1,4))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,4),rr(1,5))
         ELSE IF (NBC.eq.3) THEN
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,2),vec(1,1))
            call addvec(-1.0D0,rr(1,1),1.0D0,rr(1,3),vec(1,2))

            call vecrotvec(vec(1,1),pi*2/3,vec(1,2),vec(1,3))
            call normvec(vec(1,3),vec(1,3))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,3),vec(1,5))
            call addvec(1.0D0,vec(1,5),-1.0D0,rr(1,4),vec(1,7)) ! distance to 4
            call getnorm(vec(1,7),d1)

            call vecrotvec(vec(1,1),-pi*2/3,vec(1,2),vec(1,4))
            call normvec(vec(1,4),vec(1,4))
            call addvec(1.0D0,rr(1,1),hblength,vec(1,4),vec(1,6))
            call addvec(1.0D0,vec(1,6),-1.0D0,rr(1,4),vec(1,8)) ! distance to 4
            call getnorm(vec(1,8),d2)

            if (d1.lt.d2) then
               call copyvec(vec(1,6),rr(1,5))
            else
               call copyvec(vec(1,5),rr(1,5))
            endif
         ELSE
            write(6,*) "not able to put hydrogens back in ",hybrid,NBC
            call stopgm('MM_ADD_HYDROGEN','ERROR IN HYBRIDIZATION',
     & __LINE__,"__file__")
         ENDIF
      ELSE
        call stopgm('MM_ADD_HYDROGEN','WRONG HYBRID CARBON',
     & __LINE__,"__file__")
      ENDIF
 
      END SUBROUTINE
C     ==--------------------------------------------------------------==
