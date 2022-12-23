c     =================================================================
       SUBROUTINE mm_CAP_H(NRPT,NSOLV)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==-------------------------------------------------------------==
C       Saturate dangling bonds with H. The maximum number of dummy
C       atoms is set to 200 in the present release.
! GROMOSDIP
C     ==-------------------------------------------------------------==
cmb - Revised on 20 Apr. 2005
       use mm_input
       use system
       use mm_parallel
       use coordar
       IMPLICIT NONE

       integer NRPT, NSOLV
       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h' 
       !include 'Gromos/coordsz.h'
       !include 'Gromos/coordar.h' 

cmb    paramdummy
       INTEGER npax
       PARAMETER(npax=200) ! we might change this. FIXME: add this to ARRAYSIZES.AK
cmb
       integer n_capping
       integer dang_qm(npax)
       integer dang_mm(npax)
       integer i, m1, m2, ia, is
       logical, allocatable ::  mm_q(:)
       real*8  h_bond_d, qmx, qmy, qmz, mmx, mmy, mmz
       real*8  vecx, vecy, vecz, norm

       h_bond_d = 0.109d0

       allocate(mm_q(NRPT+NSOLV))
!$OMP parallel do private(i)
       DO i = 1,NRPT+NSOLV
         mm_q(i) = .false.
       ENDDO

!$OMP parallel do private(is,ia)
       DO is = 1, maxsys%nsx !MAXSPQ
         DO ia=1,gqmmm%gr_nasp(is) 
           mm_q(gqmmm%gr_atom(is,ia))=.true.
         ENDDO
       ENDDO

!$OMP parallel do private(i)
       DO i=1,npax
         dang_qm(i)=0
         dang_mm(i)=0
       ENDDO
!
! We identify the atoms that belongs to the QM systems and are 
! bonded to MM atoms
!

       n_capping = 0
       DO i = 1,NBON   
         m1=IB(i)
         m2=JB(i)
         if(mm_q(M1) .and. (.not. mm_q(M2))) then
           n_capping = n_capping + 1
           dang_qm(n_capping) = m1
           dang_mm(n_capping) = m2
           if(gparal%mmparent) write(6,'(A31,2I6)')
     c     "CAPPING  DANGLING BOND BETWEEN: ",m1,m2
!          write(6,'(A18,I5,A2,F7.4,A13)') "MM charge of atom ", m2,
!    c     " (",CG(m2),") set to 0.00"
!          CG(m1) = CG(m1)+CG(m2)
!          CG(m2) = 0.0d0
         else if(mm_q(M2) .and. (.not. mm_q(M1))) then
           n_capping = n_capping + 1
           dang_qm(n_capping) = m2
           dang_mm(n_capping) = m1
           if(gparal%mmparent)write(6,'(A31,2I6)')
     c    "CAPPING  DANGLING BOND BETWEEN: ",m2,m1
!          write(6,'(A18,I5,A2,F7.4,A13)') "MM charge of atom ", m1,
!    c     " (",CG(m1),") set to 0.00"
!          CG(m2) = CG(m1)+CG(m2)
!          CG(m1) = 0.0d0
         endif
       ENDDO 
       IF(gparal%mmparent.AND.n_capping.ne.(NRP-addh%igdummy+1)) THEN
          write(6,*)NRP-addh%igdummy+1,n_capping,addh%igdummy
          call stopgm('MM_CAP_H','Wrong number of dummy atoms',
     & __LINE__,"__file__")
       ENDIF
       IF(gparal%mmparent.AND.n_capping.gt.npax) THEN
         call stopgm('MM_CAP_H','Too many dummy atoms',
     & __LINE__,"__file__")
       ENDIF

!
! Finally we determine the starting coordinates for these new atoms
! from the versor that connects each quantum atom to the
! classical counterpart. The assumed distance is 0.109 nm
!

      DO i = 1, n_capping
!
!  extract coordinates of the qm  and mm atom in nanometers
!  and assign the coordinates to the capping dummy atoms
!  it assumes that the dummy atoms are located at the end of
!  the solute atoms
!
        qmx = x((dang_qm(i)-1)*3 + 1) 
        qmy = x((dang_qm(i)-1)*3 + 2) 
        qmz = x((dang_qm(i)-1)*3 + 3) 
        mmx = x((dang_mm(i)-1)*3 + 1) 
        mmy = x((dang_mm(i)-1)*3 + 2) 
        mmz = x((dang_mm(i)-1)*3 + 3) 
        vecx = qmx - mmx
        vecy = qmy - mmy
        vecz = qmz - mmz
        norm = h_bond_d/
     .         DSQRT(vecx*vecx+vecy*vecy+vecz*vecz)
        vecx = vecx*norm
        vecy = vecy*norm
        vecz = vecz*norm
        
        x(3*(addh%igdummy-1)+1)=qmx-vecx
        x(3*(addh%igdummy-1)+2)=qmy-vecy
        x(3*(addh%igdummy-1)+3)=qmz-vecz

        addh%igdummy=addh%igdummy+1

      ENDDO

      deallocate(mm_q)
      RETURN
      END
