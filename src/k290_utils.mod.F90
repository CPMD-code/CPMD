MODULE k290_utils
  USE error_handling,                  ONLY: stopgm
  USE k290_2_utils,                    ONLY: group1,&
                                             sppt2
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k290


CONTAINS

  ! ==================================================================
  SUBROUTINE k290(iout,nat,nkpoint,nsp,iq1,iq2,iq3,istriz,&
       a1,a2,a3,alat,strain,xkapa,rx,tvec,&
       ty,isc,f0,ntvec,wvk0,wvkl,lwght,lrot,&
       nhash,includ,list,rlist,delta)
    ! ==================================================================
    ! WRITTEN ON SEPTEMBER 12TH, 1979.
    ! IBM-RETOUCHED ON OCTOBER 27TH, 1980.
    ! Tsukuba-retouched on March 19th, 2008.
    ! GENERATION OF SPECIAL POINTS MODIFIED ON 26-MAY-82 BY OHN.
    ! RETOUCHED ON JANUARY 8TH, 1997
    ! INTEGRATION IN CPMD-FEMD PROGRAM BY THIERRY DEUTSCH
    ! ==--------------------------------------------------------------==
    ! PLAYING WITH SPECIAL POINTS AND CREATION OF 'CRYSTALLOGRAPHIC'
    ! FILE FOR BAND STRUCTURE CALCULATIONS.
    ! GENERATION OF SPECIAL POINTS IN K-SPACE FOR AN ARBITRARY LATTICE,
    ! FOLLOWING THE METHOD MONKHORST,PACK, PHYS. REV. B13 (1976) 5188
    ! MODIFIED BY MACDONALD, PHYS. REV. B18 (1978) 5897
    ! MODIFIED ALSO BY OLE HOLM NIELSEN ("SYMMETRIZATION")
    ! ==--------------------------------------------------------------==
    ! TESTING THEIR EFFICIENCY AND PREPARATION OF THE
    ! "STRUCTURAL" FILE FOR RUNNING THE
    ! SELF-CONSISTENT BAND STRUCTURE PROGRAMS.
    ! IN THE CASES WHERE THE POINT GROUP OF THE CRYSTAL DOES NOT
    ! CONTAIN INVERSION, THE LATTER IS ARTIFICIALLY ADDED, IN ORDER
    ! TO MAKE USE OF THE HERMITICITY OF THE HAMILTONIAN
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   IOUT    LOGIC FILE NUMBER                                  ==
    ! ==   NAT     NUMBER OF ATOMS                                    ==
    ! ==   NKPOINT MAXIMAL NUMBER OF K POINTS                         ==
    ! ==   NSP     NUMBER OF SPECIES                                  ==
    ! ==   IQ1,IQ2,IQ3 THE MONKHORST-PACK MESH PARAMETERS             ==
    ! ==   ISTRIZ  SWITCH FOR SYMMETRIZATION                          ==
    ! ==   A1(3),A2(3),A3(3) LATTICE VECTORS                          ==
    ! ==   ALAT    LATTICE CONSTANT                                   ==
    ! ==   STRAIN(3,3)  STRAIN APPLIED TO LATTICE IN ORDER            ==
    ! ==           TO HAVE K POINTS WITH SYMMETRY OF STRAINED LATTICE ==
    ! ==   XKAPA(3,NAT)   ATOMS COORDINATES                           ==
    ! ==   TY(NAT)      TYPES OF ATOMS                                ==
    ! ==   WVK0(3) SHIFT FOR K POINTS MESh (MACDONALD ARTICLE)        ==
    ! ==   NHASH   SIZE OF THE HASH TABLES (LIST)                     ==
    ! ==   DELTA   REQUIRED ACCURACY (1.e-6_real_8 IS A GOOD VALUE)          ==
    ! ==           K-VECTOR < DELTA IS CONSIDERED ZERO                ==
    ! == OUTPUT:                                                      ==
    ! ==   RX(3,NAT) SCRATCH ARRAY USED BY GROUP1 ROUTINE             ==
    ! ==   TVEC(1:3,1:NTVEC) TRANSLATION VECTORS (SEE NTVEC)          ==
    ! ==   ISC(NAT)  SCRATCH ARRAY USED BY GROUP1 ROUTINE             ==
    ! ==   F0(49,NAT) ATOM TRANSFORMATION TABLE                       ==
    ! ==              IF NTVEC/=1 THE 49TH GIVES INEQUIVALENT ATOMS   ==
    ! ==   NTVEC NUMBER OF TRANSLATION VECTORS (IF NOT PRIMITIVE CELL)==
    ! ==   WVKL(3,NKPOINT) SPECIAL KPOINTS GENERATED                  ==
    ! ==   LWGHT(NKPOINT)  WEIGHT FOR EACH K POINT                    ==
    ! ==   LROT(48,NKPOINT) SYMMETRY OPERATION FOR EACh K POINTS      ==
    ! ==   INCLUD(NKPOINT)  SCRATCH ARRAY USED BY SPPT2               ==
    ! ==   LIST(NKPOINT+NHASH) HASH TABLE USED BY SPPT2               ==
    ! ==   RLIST(3,NKPOINT) SCRATCH ARRAY USED BY SPPT2               ==
    ! ==--------------------------------------------------------------==
    ! SUBROUTINES NEEDED: 
    ! SPPT2, GROUP1, PGL1, ATFTM1, ROT1, STRUCT,
    ! BZRDUC, INBZ, MESH, BZDEFI
    ! (GROUP1, PGL1, ATFTM1, ROT1 FROM THE
    ! "COMPUTER PHYSICS COMMUNICATIONS" PACKAGE "ACMI" - (1971,1974)
    ! WORLTON-WARREN).
    ! ==================================================================
    INTEGER                                  :: iout, nat, nkpoint, nsp, iq1, &
                                                iq2, iq3, istriz
    REAL(real_8)                             :: a1(3), a2(3), a3(3), alat, &
                                                strain(6), xkapa(3,nat), &
                                                rx(3,nat), tvec(3,nat)
    INTEGER                                  :: ty(nat), isc(nat), &
                                                f0(49,nat), ntvec
    REAL(real_8)                             :: wvk0(3), wvkl(3,nkpoint)
    INTEGER                                  :: lwght(nkpoint), &
                                                lrot(48,nkpoint), nhash, &
                                                includ(nkpoint), &
                                                list(nkpoint+nhash)
    REAL(real_8)                             :: rlist(3,nkpoint), delta

    CHARACTER(len=10), DIMENSION(48), PARAMETER :: rname_cubic = (/&
      ' 1        ',' 2[ 10 0] ',' 2[ 01 0] ',' 2[ 00 1] ',' 3[-1-1-1]',&
      ' 3[ 11-1] ',' 3[-11 1] ',' 3[ 1-11] ',' 3[ 11 1] ',' 3[-11-1] ',&
      ' 3[-1-11] ',' 3[ 1-1-1]',' 2[-11 0] ',' 4[ 00 1] ',' 4[ 00-1] ',&
      ' 2[ 11 0] ',' 2[ 0-11] ',' 2[ 01 1] ',' 4[ 10 0] ',' 4[-10 0] ',&
      ' 2[-10 1] ',' 4[ 0-10] ',' 2[ 10 1] ',' 4[ 01 0] ','-1        ',&
      '-2[ 10 0] ','-2[ 01 0] ','-2[ 00 1] ','-3[-1-1-1]','-3[ 11-1] ',&
      '-3[-11 1] ','-3[ 1-11] ','-3[ 11 1] ','-3[-11-1] ','-3[-1-11] ',&
      '-3[ 1-1-1]','-2[-11 0] ','-4[ 00 1] ','-4[ 00-1] ','-2[ 11 0] ',&
      '-2[ 0-11] ','-2[ 01 1] ','-4[ 10 0] ','-4[-10 0] ','-2[-10 1] ',&
      '-4[ 0-10] ','-2[ 10 1] ','-4[ 01 0] '/)
    CHARACTER(len=11), DIMENSION(24), PARAMETER :: rname_hexai = (/&
      ' 1         ',' 6[ 00  1] ',' 3[ 00  1] ',' 2[ 00  1] ',' 3[ 00 -1] ',&
      ' 6[ 00 -1] ',' 2[ 01  0] ',' 2[-11  0] ',' 2[ 10  0] ',' 2[ 21  0] ',&
      ' 2[ 11  0] ',' 2[ 12  0] ','-1         ','-6[ 00  1] ','-3[ 00  1] ',&
      '-2[ 00  1] ','-3[ 00 -1] ','-6[ 00 -1] ','-2[ 01  0] ','-2[-11  0] ',&
      '-2[ 10  0] ','-2[ 21  0] ','-2[ 11  0] ','-2[ 12  0] '/)
    CHARACTER(len=12), DIMENSION(7), PARAMETER :: icst = (/'TRICLINIC   ',&
      'MONOCLINIC  ','ORTHORHOMBIC','TETRAGONAL  ','CUBIC       ',&
      'TRIGONAL    ','HEXAGONAL   '/)

    INTEGER :: i, ib(48), ib0(48), ihc, ihc0, ihg, ihg0, indpg, indpg0, &
      invadd, istrin, iswght, isy, isy0, itype, j, k, l, li, li0, lmax, n, &
      nc, nc0, ntot, ntvec0
    INTEGER, DIMENSION(49, 1)                :: f00 = 0
    REAL(real_8) :: a01(3), a02(3), a03(3), b01(3), b02(3), b03(3), b1(3), &
      b2(3), b3(3), dtotstr, origin(3), origin0(3), proj1, proj2, proj3, &
      r(3,3,48), r0(3,3,48), totstr, tvec0(3,1), volum, vv0(3)
    REAL(real_8), DIMENSION(3, 1)            :: x0 = 0._real_8 
    REAL(real_8), DIMENSION(3, 48)           :: v = 0._real_8, v0 = 0._real_8

! ==--------------------------------------------------------------==
! READ IN LATTICE STRUCTURE
! ==--------------------------------------------------------------==

    DO i = 1,3
       a01(i) = a1(i)/alat
       a02(i) = a2(i)/alat
       a03(i) = a3(i)/alat
    ENDDO
    IF (paral%io_parent)&
         WRITE(iout,'(" A1",3F10.5,/," A2",3F10.5,/," A3",3F10.5)')&
         a01,a02,a03
    IF (paral%io_parent)&
         WRITE(iout,'(/," NUMBER OF ATOMS (STRUCT):",I6)') nat
    IF (paral%io_parent)&
         WRITE(iout,'(5X,"K   TYPE",14X,"X(K)")')
    itype = 0
    DO i = 1,nat
       ! Assign an atomic type (for internal purposes)
       IF (i .NE. 1) THEN
          DO j = 1, (i - 1)
             IF (ty(j) .EQ. ty(i)) THEN
                ! Type located
                GOTO 178
             ENDIF
          ENDDO
          ! New type
       ENDIF
       itype = itype + 1
       IF (itype .GT. nsp) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I4,")")')&
               ' NUMBER OF ATOMIC TYPES EXCEEDS DIMENSION (NSP=)',&
               nsp
          IF (paral%io_parent)&
               WRITE(6,'(" THE ARRAY TY IS:",/,9(1X,10I7,/))')&
               (ty(j),j=1,nat)
          CALL stopgm('K290','FATAL ERROR',& 
               __LINE__,__FILE__)
       ENDIF
178    CONTINUE
       IF (paral%io_parent)&
            WRITE(iout,'(1X,I5,I6,3F10.5)')&
            i,ty(i),(xkapa(j,i),j=1,3)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! IS THE STRAIN SIGNIFICANT ?
    ! ==--------------------------------------------------------------==
    dtotstr = delta*delta
    totstr = 0._real_8
    istrin = 0
    DO i = 1,6
       totstr  = totstr + ABS(strain(i))
    ENDDO
    IF (totstr .GT. dtotstr) istrin = 1
    ! ==--------------------------------------------------------------==
    ! Volume of the cell.
    volum  = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) +&
         a3(1)*a1(2)*a2(3) - a1(3)*a2(2)*a3(1) -&
         A2(3)*A3(2)*A1(1) - A3(3)*A1(2)*A2(1)
    volum  = ABS(volum)
    b1(1) = (a2(2)*a3(3) - a2(3)*a3(2)) / volum
    b1(2) = (a2(3)*a3(1) - a2(1)*a3(3)) / volum
    b1(3) = (a2(1)*a3(2) - a2(2)*a3(1)) / volum
    b2(1) = (a3(2)*a1(3) - a3(3)*a1(2)) / volum
    b2(2) = (a3(3)*a1(1) - a3(1)*a1(3)) / volum
    b2(3) = (a3(1)*a1(2) - a3(2)*a1(1)) / volum
    b3(1) = (a1(2)*a2(3) - a1(3)*a2(2)) / volum
    b3(2) = (a1(3)*a2(1) - a1(1)*a2(3)) / volum
    b3(3) = (a1(1)*a2(2) - a1(2)*a2(1)) / volum
    ! ==--------------------------------------------------------------==
    DO i=1,3
       b01(i)=b1(i)*alat
       b02(i)=b2(i)*alat
       b03(i)=b3(i)*alat
    ENDDO
    IF (paral%io_parent)&
         WRITE(iout,'(/," B1",3F10.5,/," B2",3F10.5,/," B3",3F10.5,/)')&
         b01,b02,b03
    ! ==--------------------------------------------------------------==
    ! == GROUP-THEORY ANALYSIS OF LATTICE                             ==
    ! ==--------------------------------------------------------------==
    CALL group1(iout,a1,a2,a3,nat,ty,xkapa,b1,b2,b3,&
         ihg,ihc,isy,li,nc,indpg,ib,ntvec,&
         v,f0,r,tvec,origin,rx,isc,delta)
    ! ==--------------------------------------------------------------==
    DO n = nc+1,48
       ib(n) = 0
    ENDDO
    ! ==--------------------------------------------------------------==
    invadd=0
    IF (li .EQ. 0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(A,/,A,/,A)')&
            ' ALTHOUGH THE POINT GROUP OF THE CRYSTAL DOES NOT',&
            ' CONTAIN INVERSION, THE SPECIAL POINT GENERATION ALGORITHM',&
            ' WILL CONSIDER IT AS A SYMMETRY OPERATION'
       invadd=1
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == CRYSTALLOGRAPHIC DATA                                        ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(iout,'(/," CRYSTALLOGRAPHIC DATA:")')
    IF (paral%io_parent)&
         WRITE(iout,'(4X,"A1",3F10.5/,4X,"A2",3F10.5/,4X,"A3",3F10.5,/)')&
         a1,a2,a3
    IF (paral%io_parent)&
         WRITE(iout,'(4X,"B1",3F10.5/,4X,"B2",3F10.5/,4X,"B3",3F10.5)')&
         b1,b2,b3
    ! ==--------------------------------------------------------------==
    ! == GROUP-THEORETICAL INFORMATION                                ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(iout,'(/," GROUP-THEORETICAL INFORMATION:")')
    ! IHG .... Point group of the primitive lattice, holohedral
    IF (paral%io_parent)&
         WRITE(iout,&
         '("    POINT GROUP OF THE PRIMITIVE LATTICE: ",A," SYSTEM")')&
         icst(ihg)
    ! IHC .... Code distinguishing between hexagonal and cubic groups
    ! ISY .... Code indicating whether the space group is symmorphic
    IF (isy.EQ.0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"NONSYMMORPHIC GROUP")')
    ELSEIF (isy.EQ.1) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"SYMMORPHIC GROUP")')
    ELSEIF (isy.EQ.-1) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"SYMMORPHIC GROUP WITH NON-STANDARD ORIGIN")')
    ELSEIF (isy.EQ.-2) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"NONSYMMORPHIC GROUP???")')
    ENDIF
    ! LI ..... Inversions symmetry
    IF (li.EQ.0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"NO INVERSION SYMMETRY")')
    ELSEIF (li.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(4X,"INVERSION SYMMETRY")')
    ENDIF
    ! NC ..... Total number of elements in the point group
    IF (paral%io_parent)&
         WRITE(iout,&
         '(4X,"TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP:",I3)') nc
    IF (paral%io_parent)&
         WRITE(iout,'(4X,"TO SUM UP: (",I1,5I3,")")')&
         ihg,ihc,isy,li,nc,indpg
    ! IB ..... List of the rotations constituting the point group
    IF (paral%io_parent)&
         WRITE(iout,'(/,4X,"LIST OF THE ROTATIONS:")')
    IF (paral%io_parent)&
         WRITE(iout,'(7X,12I4)') (ib(i),i=1,nc)
    ! V ...... Nonprimitive translations (for nonsymmorphic groups)
    IF (isy.LE.0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(/,4X,"NONPRIMITIVE TRANSLATIONS:")')
       IF (paral%io_parent)&
            WRITE(iout,'(A,A)')&
            ' ROT    V  IN THE BASIS A1, A2, A3      ',&
            'V  IN CARTESIAN COORDINATES'
       ! Cartesian components of nonprimitive translation.
       DO i=1,nc
          DO j=1,3
             vv0(j) = v(1,i)*a1(j) + v(2,i)*a2(j) + v(3,i)*a3(j)
          ENDDO
          IF (paral%io_parent)&
               WRITE(iout,'(1X,I3,3F10.5,3X,3F10.5)')&
               ib(i),(v(j,i),j=1,3),vv0
       ENDDO
    ENDIF
    ! F0 ..... The function defined in Maradudin, Ipatova by
    ! eq. (3.2.12): atom transformation table.
    IF (paral%io_parent)&
         WRITE(iout,&
         '(/,4X,"ATOM TRANSFORMATION TABLE (MARADUDIN,VOSKO):")')
    IF (paral%io_parent)&
         WRITE(iout,'(5(4X,"R  AT->AT"))')
    IF (paral%io_parent)&
         WRITE(iout,'(I5," [Identity]")') 1
    DO k=2,nc
       DO j=1,nat
          IF (paral%io_parent)&
               WRITE(iout,'(I5,2I4)',advance="no") ib(k),j,f0(k,j)
          IF ((MOD(j,5).EQ.0).AND.paral%io_parent)&
               WRITE(iout, *)
       ENDDO
       IF ((MOD(j-1,5).NE.0).AND.paral%io_parent)&
            WRITE(iout, *)
    ENDDO
    ! R ...... List of the 3 x 3 rotation matrices
    IF (paral%io_parent)&
         WRITE(iout,'(/,4X,"LIST OF THE 3 X 3 ROTATION MATRICES:")')
    IF (ihc.EQ.0) THEN
       DO k=1,nc
          IF (paral%io_parent)&
               WRITE(iout,&
               '(4X,I3," (",I2,": ",A11,")",2(3F14.6,/,25X),3F14.6)')&
               k,ib(k),rname_hexai(ib(k)),((r(i,j,ib(k)),j=1,3),i=1,3)
       ENDDO
    ELSE
       DO k=1,nc
          IF (paral%io_parent)&
               WRITE(iout,&
               '(4X,I3," (",I2,": ",A10,") ",2(3F14.6,/,25X),3F14.6)')&
               k,ib(k),rname_cubic(ib(k)),((r(i,j,ib(k)),j=1,3),i=1,3)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! GENERATE THE BRAVAIS LATTICE
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(iout,'(/,1X,20("-"),A,20("-"),/,A)')&
         ' THE (UNSTRAINED) BRAVAIS LATTICE ',&
         ' (USED FOR GENERATING THE LARGEST POSSIBLE MESH IN THE B.Z.)'
    ! ==--------------------------------------------------------------==
    CALL group1(iout,a01,a02,a03,1,ty,x0,b01,b02,b03,&
         ihg0,ihc0,isy0,li0,nc0,indpg0,ib0,ntvec0,&
         v0,f00,r0,tvec0,origin0,rx,isc,delta)
    ! ==--------------------------------------------------------------==
    ! It is assumed that the same 'type' of symmetry operations
    ! (cubic/hexagonal) will apply to the crystal as well as the Bravais
    ! lattice.
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(iout,'(/,1X,19("*"),A,25("*"))')&
         ' GENERATION OF SPECIAL POINTS '
    ! Parameter Q of Monkhorst and Pack, generalized for 3 axes B1,2,3
    IF (paral%io_parent)&
         WRITE(iout,'(A,/,1X,3I5)')&
         ' MONKHORST-PACK PARAMETERS (GENERALIZED) IQ1,IQ2,IQ3:',&
         iq1,iq2,iq3
    ! WVK0 is the shift of the whole mesh (see Macdonald)
    IF (paral%io_parent)&
         WRITE(iout,'(A,/,1X,3F10.5)')&
         ' CONSTANT VECTOR SHIFT (MACDONALD) OF THIS MESH:',wvk0
    IF (iabs(iq1) + iabs(iq2) + iabs(iq3) .EQ. 0) GOTO 710
    IF (ABS(istriz) .NE. 1) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(" INVALID SWITCH FOR SYMMETRIZATION",I10)') istriz
       IF (paral%io_parent)&
            WRITE(6,'(" INVALID SWITCH FOR SYMMETRIZATION",I10)') istriz
       CALL stopgm('K290','ISTRIZ WRONG ARGUMENT',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(iout,'(" SYMMETRIZATION SWITCH: ",I3)',advance="no") istriz
    IF (istriz.EQ.1) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(" (SYMMETRIZATION OF MONKHORST-PACK MESH)")')
    ELSE
       IF (paral%io_parent)&
            WRITE(iout,'(" (NO SYMMETRIZATION OF MONKHORST-PACK MESH)")')
    ENDIF
    ! Set to 0.
    DO i=1, nkpoint
       lwght(i) = 0
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == Generation of the points (they are not multiplied            ==
    ! == by 2*Pi because B1,2,3  were not,either)                     ==
    ! ==--------------------------------------------------------------==
    IF (nc.GT.nc0) THEN
       ! Due to non-use of primitive cell, the crystal has more
       ! rotations than Bravais lattice.
       ! We use only the rotations for Bravais lattices
       IF (ntvec.EQ.1) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'NUMBER OF ROTATIONS FOR BRAVAIS LATTICE',nc0
          IF (paral%io_parent)&
               WRITE(6,*) 'NUMBER OF ROTATIONS FOR CRYSTAL LATTICE',nc
          IF (paral%io_parent)&
               WRITE(6,*) 'NO DUPLICATION FOUND'
          CALL stopgm('ERROR',&
               'SOMETHING IS WRONG IN GROUP DETERMINATION',& 
               __LINE__,__FILE__)
       ENDIF
       nc=nc0
       DO i=1,nc0
          ib(i)=ib0(i)
       ENDDO
       IF (paral%io_parent)&
            WRITE(iout,'(/,1X,20("! "),"WARNING",20("!"))')
       IF (paral%io_parent)&
            WRITE(iout,'(A)')&
            ' THE CRYSTAL HAS MORE SYMMETRY THAN THE BRAVAIS LATTICE'
       IF (paral%io_parent)&
            WRITE(iout,'(A)')&
            ' BECAUSE THIS IS NOT A PRIMITIVE CELL'
       IF (paral%io_parent)&
            WRITE(iout,'(A)')&
            ' USE ONLY SYMMETRY FROM BRAVAIS LATTICE'
       IF (paral%io_parent)&
            WRITE(iout,'(1X,20("! "),"WARNING",20("!"),/)')
    ENDIF
    CALL sppt2 (iout,iq1,iq2,iq3,wvk0,nkpoint,&
         a01,a02,a03,b01,b02,b03,&
         invadd,nc,ib,r,ntot,wvkl,lwght,lrot,nc0,ib0,istriz,&
         nhash,includ,list,rlist,delta)
    ! ==--------------------------------------------------------------==
    ! == Check on error signals                                       ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(iout,'(/,1X,I5," SPECIAL POINTS GENERATED")') ntot
    IF (ntot .EQ. 0) THEN
       GOTO 710
    ELSE IF (ntot .LT. 0) THEN
       IF (paral%io_parent)&
            WRITE(iout,'(A,I5,/,A,/,A)') ' DIMENSION NKPOINT =',nkpoint,&
            ' INSUFFICIENT FOR ACCOMMODATING ALL THE SPECIAL POINTS',&
            ' WHAT FOLLOWS IS AN INCOMPLETE LIST'
       ntot = iabs(ntot)
    ENDIF
    ! Before using the list WVKL as wave vectors, they have to be
    ! multiplied by 2*Pi
    ! The list of weights LWGHT is not normalized
    iswght = 0
    DO i=1,ntot
       iswght = iswght + lwght(i)
    ENDDO
    IF (paral%io_parent)&
         WRITE(iout,'(8X,A,T33,A,4X,A)')&
         'WAVEVECTOR K','WEIGHT','UNFOLDING ROTATIONS'
    ! Set near-zeroes equal to zero:
    DO l=1,ntot
       DO i = 1,3
          IF (ABS(wvkl(i,l)) .LT. delta) wvkl(i,l) = 0._real_8
       ENDDO
       IF (istrin .NE. 0) THEN
          ! Express special points in (unstrained) basis.
          proj1 = 0._real_8
          proj2 = 0._real_8
          proj3 = 0._real_8
          DO i = 1,3
             proj1 = proj1 + wvkl(i,l)*a01(i)
             proj2 = proj2 + wvkl(i,l)*a02(i)
             proj3 = proj3 + wvkl(i,l)*a03(i)
          ENDDO
          DO i = 1,3
             wvkl(i,l) = proj1*b1(i) + proj2*b2(i) + proj3*b3(i)
          ENDDO
       ENDIF
       lmax=lwght(l)
       IF (paral%io_parent)&
            WRITE(iout,fmt='(1X,I5,3F8.4,I8,T42,12I3)')&
            l,(wvkl(i,l),i=1,3),lwght(l),(lrot(i,l),i=1,MIN(lmax,12))
       DO j=13,lmax,12
          IF (paral%io_parent)&
               WRITE(iout,fmt='(T42,12I3)')&
               (lrot(i,l),i=j,MIN(lmax,j-1+12))
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(iout,'(24X,"TOTAL:",I8)') iswght
    ! ==--------------------------------------------------------------==
710 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k290
  ! ==================================================================

END MODULE k290_utils
