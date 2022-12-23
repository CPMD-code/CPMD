C     ==================================================================
       SUBROUTINE mm_get_NSX(NSPq,NRPT,NSOLV)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
! GROMOSDIP
! determines new maxsys%nax and maxsys%nsx
       use system
       use ions
       use mm_input
       use mm_parallel
       use coordsz
       implicit none

       integer NRPT,NSOLV
       integer NSPq            
!
       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h'
       !include 'Gromos/coordsz.h'
!
       logical, allocatable ::  mm_q(:)
       integer niac(MAXATT),sn
       integer i,is,ia,ik,NRTOT

       NRTOT=NRPT+NSOLV
       allocate(mm_q(NRTOT))
       NSPq = maxsys%nsx
       IF (gparal%mmparent) THEN
         write(6,*)
         write(6,*)'  QUANTUM SYSTEM:'
         write(6,*)'  NAX:',maxsys%nax
         write(6,*)'  NSX:',maxsys%nsx
         write(6,*)
         IF (MAXSP-maxsys%nsx-MAXATT.le.0) THEN
           write(6,*)'MAXSP=',MAXSP
           write(6,*)'nsx=',maxsys%nsx,'  MAXATT=',MAXATT
           write(6,*)'MAXSP-nsx-MAXATT=',MAXSP-maxsys%nsx-MAXATT
           write(6,*)'MAXSP is too small, I am sorry...'
           call stopgm('mm_init','  ',
     & __LINE__,"__file__")
        ENDIF
       ENDIF
c   get a good estimate for maxsys%nax
       maxsys%nax=max(maxsys%nax,NRTOT/(MAXSP-maxsys%nsx-MAXATT))            
c
c   find quantum atoms
c
       do i=1,NRTOT
         mm_q(i)=.false.
         do is=1,NSPq
           do ia=1,gqmmm%gr_nasp(is)
             if(gqmmm%gr_atom(is,ia).eq.i)then
               if(i.gt.NRPT) call stopgm('mm_get_NSX',
     c          'solvent moleules can not be quantum',
     & __LINE__,"__file__")
               mm_q(i)=.true.
             endif
           enddo
         enddo
       enddo
c
c   given the iac code of the mm atoms as defined in gromos topology, 
c   defines the cp species and atom number
c
!$OMP parallel do private(i)
       DO i=1,MAXATT
         niac(i)=0
       ENDDO
c   solute
       do i=1,NRPT
         if(iac(i).gt.MAXATT)  then
           write(6,*)i,iac(i),MAXATT
           call stopgm('mm_get_NSX','WRONG IAC',
     & __LINE__,"__file__")
         endif
         if(.not. mm_q(i)) niac(iac(i))=niac(iac(i))+1
       enddo
       do i=NRPT+1,NRTOT
         sn=mod(i-1-NRPT,NRAM)+1
         if(iacs(sn).gt.MAXATT)  then
           write(6,*)sn,iacs(sn),MAXATT
           call stopgm('mm_get_NSX','WRONG IAC',
     & __LINE__,"__file__")
         endif
         if(.not. mm_q(i)) niac(iacs(sn))=niac(iacs(sn))+1
       enddo

 100   continue
       do ik=1,MAXATT
c        write(6,*)'niac',niac(ik)
         if(niac(ik).gt.0) then
            maxsys%nsx=maxsys%nsx+(niac(ik)-1)/maxsys%nax+1
         endif
       enddo
c   we may have an unfortunate distribution of the
c   of the MM atom types. crank up maxsys%nax and try again.
       IF(maxsys%nsx.GT.MAXSP) THEN
         maxsys%nax=maxsys%nax+1
         maxsys%nsx=NSPq
         GOTO 100
       ENDIF

       IF (gparal%mmparent) THEN
          write(6,*)'  FULL SYSTEM:'
          write(6,*)'  NAX:',maxsys%nax
          write(6,*)'  NSX:',maxsys%nsx
          write(6,*)
       END IF

       deallocate(mm_q)
       RETURN
       END
