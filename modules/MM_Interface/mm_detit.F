C     ==================================================================
       SUBROUTINE mm_detit(cpat,cpsp,NAT_cpmd,NAT_grm,gratom,
     .            NAq,NSPq,NAm,NSPm,NATq,NAXq,NATm,NRPT,NSOLV)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 23 Feb. 2013
! GROMOS_DEPENDENCES
       use system
       use parac
       use ions
       use mm_input
       use fileopenmod
       use fileopen_utils, only : fileclose,fileopen

       use cotr
       use coordsz
       implicit none

       integer NAT_cpmd(0:*)
       integer NAT_grm(0:*)
       integer cpat(*)
       integer cpsp(*)
       integer gratom(*)
       integer NRPT,NSOLV
       integer NSPq            !number of quantum species
       integer NSPm            !number of total species
       integer NAq(*)      !number of quantum atoms per species
       integer NAm(*)      !number of total atoms per species
       integer NATq           !number of quantum atoms
       integer NAXq           !max number of quantum atoms per species
       integer NATm           !number of atoms
       
       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h'
       !include 'Gromos/coordsz.h'
c locals
c
       integer niac(MAXATT)
       integer SP(MAXATT)
       integer ik,spt,sn,ig
       integer i,im,is,j,ia
       integer NRTOT
       logical ferror
       logical , allocatable :: mm_q(:)
c
c  quantum atoms 
c
       NAXq=1
       NRTOT=NRPT+NSOLV
       IF(NRP.ne.NRPT) THEN
         call stopgm('mm_detit',
     .        'MULTIPLE SOLUTE MOLECULES NOT IMPLEMENTED',
     & __LINE__,"__file__")
       ENDIF
       allocate( mm_q(NRTOT) )
       do i=1,NRTOT
         mm_q(i)=.false.
         do is=1,NSPq
           NAXq=MAX(NAXq,gqmmm%gr_nasp(is))
           NAq(is)=gqmmm%gr_nasp(is)
           do ia=1,gqmmm%gr_nasp(is)
             if(gqmmm%gr_atom(is,ia).eq.i)then
               mm_q(i)=.true. 
               cpat(i)=ia
               cpsp(i)=is
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

       do i=1,NRPT
         if(.not. mm_q(i)) niac(iac(i))=niac(iac(i))+1
       enddo

       do i=NRPT+1,NRPT+NSOLV
         sn=mod(i-1-NRPT,NRAM)+1
         if(.not. mm_q(i)) niac(iacs(sn))=niac(iacs(sn))+1
       enddo

       NSPm=NSPq
       do ik=1,MAXATT
         if(niac(ik).gt.0) then
            SP(ik)=NSPm+1
            NSPm=NSPm+(niac(ik)-1)/maxsys%nax+1
         endif
       enddo
       if(NSPm.ne.maxsys%nsx)
     &   call stopgm('mm_detit','error defining nsx',
     & __LINE__,"__file__")
       

c   gromos atom i corresponds to  cp atom cpat(i) of species cpsp(i) 
!$OMP parallel do private(i)
       DO i=1,NSPq
         NAm(i)=NAq(i)
       ENDDO
!$OMP parallel do private(i)
       DO i=NSPq+1,NSPm
         NAm(i)=0
       ENDDO

       do i=1,NRPT
         if(.not. mm_q(i))then
           spt=sp(iac(i))
           NAm(spt)=NAm(spt)+1
           if(NAm(spt).gt.maxsys%nax)then
             NAm(spt)=maxsys%nax
             SP(iac(i))=SP(iac(i))+1
             spt=sp(iac(i))
             NAm(spt)=1
           endif
           cpsp(i)=spt
           cpat(i)=NAm(spt)
         endif
       enddo

       do i=NRPT+1,NRPT+NSOLV
         sn=mod(i-1-NRPT,NRAM)+1
         if(.not. mm_q(i))then
           spt=sp(iacs(sn))
           NAm(spt)=NAm(spt)+1
           if(NAm(spt).gt.maxsys%nax)then
             NAm(spt)=maxsys%nax
             SP(iacs(sn))=SP(iacs(sn))+1
             spt=sp(iacs(sn))
             NAm(spt)=1
           endif
           cpsp(i)=spt
           cpat(i)=NAm(spt)
         endif
       enddo

       NATq=0
!$OMP parallel do private(is) reduction(+:NATq)
       DO is=1,NSPq
         NATq=NATq+NAq(is)
       ENDDO 
       NATm=0
!$OMP parallel do private(is) reduction(+:NATm)
       DO is=NSPq+1,NSPm
         NATm=NATm+NAm(is)
       ENDDO
       NATm=NATm+NATq

       if (NATm.ne.NRTOT) then
         write(6,*)NATq,NATm,NRTOT
         call stopgm('mm_detit','error defining NATm',
     & __LINE__,"__file__")
       endif

C NAT_cpmd translates Gromos ordering into CPMD ordering
       do ig=1,NRTOT
         NAT_cpmd(ig)=0
         do is=1,cpsp(ig)-1
           NAT_cpmd(ig)= NAT_cpmd(ig) + NAm(is)
         enddo
         NAT_cpmd(ig)= NAT_cpmd(ig) + cpat(ig)
       enddo
C      Fills the DUMMY ATOM indices in NAT_cpmd
       if(duat%NDAT.gt.0)then
         WRITE(6,'(/,A)')' DUMMY ATOMS INDICIES: #DUM,  INDEX'
         do i=NRTOT+1,NRTOT+duat%NDAT
           NAT_cpmd(i)=i
           WRITE(6,'(20X,I5,I8)') IA, I
         end do
       end if

!$OMP parallel do private(i) shared(maxsys)
       DO i=1,NRTOT
         gratom((cpsp(i)-1)*maxsys%nax+cpat(i))=i
       ENDDO
C     write table of the internal ordering.
       IF(paral%io_parent) THEN  ! only the parent must write on 17 !
       CALL FILEOPEN(17,'QMMM_ORDER',FO_DEF,FERROR)
       WRITE(17,'(i8,1x,i8)')NRTOT,NATq
       ENDIF
       i=0
       NAT_grm(0)=0
       NAT_cpmd(0)=0
       DO is=1,NSPm
         DO ia=1,NAm(is)
           i=i+1
           IF(paral%io_parent) THEN
           WRITE(17,'(4(i8,1x))')gratom((is-1)*maxsys%nax+ia),i,is,ia
           ENDIF
C NAT_grm translates CPMD ordering into Gromos ordering
           NAT_grm(i)=gratom((is-1)*maxsys%nax+ia)
         ENDDO
       ENDDO
C      Fills the DUMMY ATOM indices in NAT_grm
       IF(duat%NDAT.gt.0) THEN
         do i=NRTOT+1,NRTOT+duat%NDAT
           NAT_grm(i)=i
         ENDDO
       ENDIF
       IF(paral%io_parent) CALL FILECLOSE(17)

       deallocate(mm_q)
C     ==--------------------------------------------------------------==
       RETURN
       END
