C     ==================================================================
         PROGRAM Amber2Gromos
C     ==--------------------------------------------------------------==
C     ==  Converts Amber topology and coordinates files into Gromos   ==
C     ==  WARNING: AMBER-6, -7 and -12 have different formats !       ==
C     ==  THE PRESENT VERSION WORKS ALSO WITH AMBER v7, 8, 9 and 12   ==
C     ==  FORMATS AND THE CORRESPONDING MAKEFILE AUTOMATICALLY USES   ==
C     ==  THE SELECTED KEYWORD                                        ==
C     ==  Revision cmb Tokyo/Strasbourg june 2014                     ==
C     ==--------------------------------------------------------------==
         IMPLICIT NONE
         CHARACTER*60 amber_crd,gromos_crd
         CHARACTER*60 amber_topofile,gromos_topofile
         CHARACTER*60 input_gromos
         REAL*8  box(3),Rcut_alpha,Rcut
         INTEGER NRP,NT_UPDATE_NLIST,ngp_x,ngp_y,ngp_z,NTB,NSM,NRAM
         LOGICAL fexist,flexwater,readpdb
         INTEGER iargc,nargs,argoffs

         WRITE(6,*) '==----------------------------------------------=='
#if defined(_AMBER6)
         WRITE(6,*) ' Converting an AMBER-6 style topology to GROMOS'
#elif defined(_AMBER7)
         WRITE(6,*) ' Converting an AMBER-7 style topology to GROMOS'
#elif defined(_AMBER12)
         WRITE(6,*) ' Converting an AMBER-12 style topology to GROMOS'
#else
         PRINT*, 'Unsupported topology format. Check compilation flags'
         STOP ' COMPILATION PROBLEM: No AMBER version specified !'
#endif
         PRINT*
         PRINT*, '   Program version 1.00 (5/06/2014)'
         PRINT*
         WRITE(6,*) '==----------------------------------------------=='

C AK: FIXME: check this with cpmd2cube and cpmd for portability.         
         argoffs=0
         flexwater=.false.
         readpdb=.false.
         nargs=iargc()
 999     CONTINUE
         IF(nargs.lt.2) THEN
#if defined(_AMBER6)
           PRINT 1, 6
#elif defined(_AMBER7)
           PRINT 1, 7
#elif defined(_AMBER12)
           PRINT 2, 12
#endif
           STOP
         ENDIF
 1       FORMAT (' usage: amber',i1,'togromos.x [-flexible|-rigid]',
     .        ' [-pdb] <topology> <coordinates>',/)
 2       FORMAT (' usage: amber',i2,'togromos.x [-flexible|-rigid]',
     .        ' [-pdb] <topology> <coordinates>',/)
         
         CALL getarg(1+argoffs,amber_topofile)
C     parse potential flags
         IF(amber_topofile(1:3).eq.'-fl') THEN
           PRINT*,'FLEXIBLE water molecules requested'
           nargs=nargs-1
           argoffs=argoffs+1
           flexwater=.true.
           GOTO 999
         ENDIF
         IF(amber_topofile(1:3).eq.'-ri') THEN
           PRINT*,'RIGID water molecules requested'
           nargs=nargs-1
           argoffs=argoffs+1
           flexwater=.false.
           GOTO 999
         ENDIF
         IF(amber_topofile(1:4).eq.'-pdb') THEN
           PRINT*,'reading coordinates from PDB file requested'
           nargs=nargs-1
           argoffs=argoffs+1
           readpdb=.true.
           GOTO 999
         ENDIF

         CALL getarg(2+argoffs,amber_crd)
         gromos_topofile='gromos.top'
         input_gromos='gromos.inp'
         gromos_crd='gromos.crd'

C     test for existence of files.
         fexist=.false.
         INQUIRE(file=amber_topofile,exist=fexist)
         IF(.not.fexist) THEN
           PRINT*, 'AMBER topology file ', amber_topofile, ' not found.'
           STOP
         ENDIF
         INQUIRE(file=amber_crd,exist=fexist)
         IF(.not.fexist) THEN
           PRINT*, 'AMBER coordinates file ', amber_crd, ' not found.'
           STOP
         ENDIF
         NSM=0
         NRAM=1
         CALL amber_to_gromos(amber_topofile,gromos_topofile,
     c              BOX,NRP,NTB,NSM,NRAM,flexwater)
         CALL amber_to_gromos_coord(amber_crd,gromos_crd,box,
     c        NTB,readpdb)
         

! HARD CODED INPUT VALUES
         NT_UPDATE_NLIST=10
         Rcut=1.0
         Rcut_alpha=0.8d0
         ngp_x=64
         ngp_y=64
         ngp_z=64
         call write_gromos_inp(input_gromos,
     c   ngp_x,ngp_y,ngp_z,NT_UPDATE_NLIST,NRP,
     c   Rcut_alpha,Rcut,box,NTB,NSM,NRAM) 
         STOP
         END
C     ==================================================================
         SUBROUTINE amber_to_gromos(amber_topofile,gromos_topofile,
     &     BOX,NTOTAT,NTB,NSM_r,NRAM_r,flexwater)
C     ==--------------------------------------------------------------==
! gromos
         IMPLICIT NONE
         CHARACTER*30 amber_topofile,gromos_topofile
         LOGICAL flexwater
         INTEGER NTB
         REAL*8  BOX(3)
         INCLUDE 'toposz.h'
         INCLUDE 'topoar.h'
! amber
         CHARACTER*4 ITITL
         INTEGER NTOTAT,NTYPES,NBONH_a,NBONA,NTHETH,NTHETS
         INTEGER NPHIH_a,NPHIA,JHPARM,NTTYP,NTYPES_d,NTTYP_d
         INTEGER JPARM,NEXT, NTOTRS,MBONA,MTHETS,MPHIA,MUMBND,MUMANG   
         INTEGER MPTRA,NATYP,NHB,IFPERT,NBPER,NGPER,NSM_r,NRAM_r
         INTEGER NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP      
         INTEGER NATCAP                  ! Water cap parameters
         REAL*8  CUTCAP,XCAP,YCAP,ZCAP   ! Water cap parameters
         CHARACTER*4, ALLOCATABLE :: IGRAPH(:)
         REAL*8, ALLOCATABLE :: CHRG(:)
         REAL*8, ALLOCATABLE :: AMASS(:)
#if defined(_AMBER12)
         INTEGER, ALLOCATABLE :: ZATOM(:)
#endif
         INTEGER, ALLOCATABLE :: IAC_a(:)
         INTEGER, ALLOCATABLE :: NUMEX(:)
         INTEGER, ALLOCATABLE :: ICO(:)
         CHARACTER*4,ALLOCATABLE :: LABRES(:)
         INTEGER, ALLOCATABLE :: IPRES(:)
         REAL*8, ALLOCATABLE :: RK(:)
         REAL*8, ALLOCATABLE :: REQ(:)
         REAL*8, ALLOCATABLE :: TK(:)
         REAL*8, ALLOCATABLE :: TEQ(:)
         REAL*8, ALLOCATABLE :: PK(:)
         REAL*8, ALLOCATABLE :: PN(:)
         REAL*8, ALLOCATABLE :: PHASE(:)
#if defined(_AMBER12)
         REAL*8, ALLOCATABLE :: SCEE(:)
         REAL*8, ALLOCATABLE :: SCNB(:)
#endif
         REAL*8, ALLOCATABLE :: SOLTY(:)
         REAL*8, ALLOCATABLE :: CN1(:)
         REAL*8, ALLOCATABLE :: CN2(:)
         INTEGER, ALLOCATABLE :: IBH_a(:),JBH_a(:),ICBH_a(:)
         INTEGER, ALLOCATABLE :: IB_a(:),JB_a(:),ICB_a(:)
         INTEGER, ALLOCATABLE :: ITH_a(:),JTH_a(:)
         INTEGER, ALLOCATABLE :: KTH_a(:),ICTH_a(:)
         INTEGER, ALLOCATABLE :: IT_a(:),JT_a(:),KT_a(:),ICT_a(:)
         INTEGER, ALLOCATABLE :: IPH_a(:),JPH_a(:),KPH_a(:)
         INTEGER, ALLOCATABLE :: LPH_a(:),ICPH_a(:)
         INTEGER, ALLOCATABLE :: IP_a(:),JP_a(:)
         INTEGER, ALLOCATABLE :: KP_a(:),LP_a(:),ICP_a(:)
         INTEGER, ALLOCATABLE :: NATEX(:)
         REAL*8, ALLOCATABLE :: AG(:)
         REAL*8, ALLOCATABLE :: BG(:)
         REAL*8, ALLOCATABLE :: HBCUT(:)
         CHARACTER*4, ALLOCATABLE :: ISYMBL(:)
         CHARACTER*4, ALLOCATABLE :: ITREE(:)
         REAL*8  BOX_a(3)
         INTEGER, ALLOCATABLE :: JOIN(:)
         INTEGER, ALLOCATABLE :: IROTAT(:)
         INTEGER, ALLOCATABLE :: NSP(:)
         INTEGER IPTRES,NSPM,NSPSOL
!  locals
         REAL*8 conv_cg,conv_mas,conv_cb,conv_t0
         REAL*8 conv_c12,conv_c6,fac_14,fac_14vdw
         REAL*8 beta_a,to_rad,pi,twoekbt
         REAL*8 tol_cg
         REAL*8,parameter :: kcal_to_kj=4.184d0      !kcal/mol to kj/mol
         REAL*8,parameter :: ekbt=2.5d0              !E*kb*T in kJ/mol, adapted to the value used in the gromos code
         REAL*8,parameter :: fconst_rescale=2.0d0    !Amber uses k(X-X0)**2 instead of 0.5*k(X-X0)**2, 
                                                     !as done in GROMOS, for harmonic potentials, so we multiply
                                                     !the corresponding force constants by fconst_rescale=2.0 
         REAL*8,parameter :: aa_to_nm=0.1d0          !Angstrom to nanometer (see cell sizes)


         INTEGER i,j,i1,i2,NF,ibbb
         INTEGER, ALLOCATABLE :: EXCL14(:,:)
         INTEGER, ALLOCATABLE :: EXCL(:,:)
         INTEGER sav(40),jsmall,jf,smaller

         REAL*8 cht, ch1, ch2, tol_ch_1, tol_ch_2, d_ch1,d_ch12
         INTEGER ic, max_bond
         INTEGER, ALLOCATABLE :: n_bond(:)
         INTEGER, ALLOCATABLE :: bond(:,:)
         INTEGER, ALLOCATABLE :: stai_con(:)
         REAL*8, ALLOCATABLE :: ch_grp(:)

         INTEGER NPM,NSM,NTX,IG,NTXO
         COMMON /MSTART/ NPM,NSM,NTX,IG,NTXO

         tol_cg=0.1
!
! 
! WARNING we assume that the scalign factor for the non bonding vdw interactions is
! 2.0 this should work  fine for the AMBER force field.
!
         OPEN(12,file=amber_topofile,status='old')
         OPEN(10,file=gromos_topofile,status='unknown')
         NF=12
! Warning: Amber 7 and later vesions print out many more comment
!          lines that we just skip with READ(NF,FMT=*,ERR=666,END=666) ! cmb
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9108,ERR=666,END=666) ITITL
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) NTOTAT,
     +            NTYPES,   !(# vdw types)
! bond
     +            NBONH_a,
     +            NBONA,
!angles
     +            NTHETH,
     +            NTHETS,
! dihedrals
     +            NPHIH_a,
     +            NPHIA,
     +            JHPARM,    ! read but not used
     +            JPARM,     ! read but not used
     +            NEXT,      ! number of excluded atoms
     +            NTOTRS,    ! number of residues
     +            MBONA,     ! number of bons + number of constrained distances
     +            MTHETS,    ! number of angles + number of constrained angles
     +            MPHIA,     ! number of dihedrals + number of constrained dihedrals
     +            MUMBND,    ! number of bond types 
     +            MUMANG,    ! number of angles type
     +            MPTRA,     ! number of dihedrals parameters types
     +            NATYP,     ! number of Force Field atom types (??)
     +            NHB,       ! number of Hydrogen bond types (??)
     +            IFPERT,    ! flag for perturbation
     +            NBPER,     ! flag for perturbation
     +            NGPER,     ! flag for perturbation
     +            NDPER,     ! flag for perturbation
     +            MBPER,     ! flag for perturbation
     +            MGPER,     ! flag for perturbation
     +            MDPER,     ! flag for perturbation
     +            IFBOX,     ! periodic bondary conditions (1 box, 2 octa)
     +            NMXRS,     ! number of atoms in the largest residue
     +            IFCAP      ! logical flag for water CAP
#else
         READ(NF,9108,ERR=666,END=666) ITITL
         READ(NF,9118,ERR=666,END=666) NTOTAT,
     +            NTYPES,   !(# vdw types)
! bond
     +            NBONH_a,
     +            NBONA,
!angles
     +            NTHETH,
     +            NTHETS,
! dihedrals
     +            NPHIH_a,
     +            NPHIA,
     +            JHPARM,    ! read but not used
     +            JPARM,     ! read but not used
     +            NEXT,      ! number of excluded atoms
     +            NTOTRS,    ! number of residues
     +            MBONA,     ! number of bons + number of constrained distances
     +            MTHETS,    ! number of angles + number of constrained angles
     +            MPHIA,     ! number of dihedrals + number of constrained dihedrals
     +            MUMBND,    ! number of bond types 
     +            MUMANG,    ! number of angles type
     +            MPTRA,     ! number of dihedrals parameters types
     +            NATYP,     ! number of Force Field atom types (??)
     +            NHB,       ! number of Hydrogen bond types (??)
     +            IFPERT,    ! flag for perturbation
     +            NBPER,     ! flag for perturbation
     +            NGPER,     ! flag for perturbation
     +            NDPER,     ! flag for perturbation
     +            MBPER,     ! flag for perturbation
     +            MGPER,     ! flag for perturbation
     +            MDPER,     ! flag for perturbation
     +            IFBOX,     ! periodic bondary conditions (1 box, 2 octa)
     +            NMXRS,     ! number of atoms in the largest residue
     +            IFCAP      ! logical flag for water CAP
#endif
         write(6,*) 'NTOTAT=',NTOTAT ! cmb
         write(6,*) 'IFCAP =',IFCAP  ! cmb
         NTTYP  = (NTYPES*(NTYPES+1))/2
         ALLOCATE (IGRAPH(NTOTAT))
         ALLOCATE (CHRG(NTOTAT))
#if defined(_AMBER12)
         ALLOCATE (ZATOM(NTOTAT))
#endif
         ALLOCATE (AMASS(NTOTAT))
         ALLOCATE (IAC_a(NTOTAT))
         ALLOCATE (NUMEX(NTOTAT))
         ALLOCATE (ICO(NTYPES*NTYPES))
         ALLOCATE (LABRES(NTOTRS))
         ALLOCATE (IPRES(NTOTRS))
         ALLOCATE (RK(MUMBND))
         ALLOCATE (REQ(MUMBND))
         ALLOCATE (TK(MUMANG))
         ALLOCATE (TEQ(MUMANG))
         ALLOCATE (PK(MPTRA))
         ALLOCATE (PN(MPTRA))
         ALLOCATE (PHASE(MPTRA))
#if defined(_AMBER12)
         ALLOCATE (SCEE(MPTRA))
         ALLOCATE (SCNB(MPTRA))
#endif
         ALLOCATE (SOLTY(NATYP))
         ALLOCATE (CN1(NTTYP))
         ALLOCATE (CN2(NTTYP))
         ALLOCATE (IBH_a(NBONH_a),JBH_a(NBONH_a),ICBH_a(NBONH_a))
         ALLOCATE (IB_a(MBONA),JB_a(MBONA),ICB_a(MBONA))
         ALLOCATE (ITH_a(NTHETH),JTH_a(NTHETH))
         ALLOCATE (KTH_a(NTHETH),ICTH_a(NTHETH))
         ALLOCATE (IT_a(MTHETS),JT_a(MTHETS),KT_a(MTHETS),ICT_a(MTHETS))
         ALLOCATE (IPH_a(NPHIH_a),JPH_a(NPHIH_a),KPH_a(NPHIH_a))
         ALLOCATE (LPH_a(NPHIH_a),ICPH_a(NPHIH_a))
         ALLOCATE (IP_a(MPHIA),JP_a(MPHIA),KP_a(MPHIA))
         ALLOCATE (LP_a(MPHIA),ICP_a(MPHIA))
         ALLOCATE (NATEX(NEXT))
         ALLOCATE (AG(NHB))
         ALLOCATE (BG(NHB))
         ALLOCATE (HBCUT(NHB))
         ALLOCATE (ISYMBL(NTOTAT))

#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9108,ERR=666,END=666) (IGRAPH(I),I = 1,NTOTAT)   ! atom names
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (CHRG(I),  I = 1,NTOTAT)   ! charges on atoms
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
#if defined(_AMBER12)
         READ(NF,9117,ERR=666,END=666) (ZATOM(I), I = 1,NTOTAT)   ! atomic numbers on Amber12
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (AMASS(I), I = 1,NTOTAT)   ! atom masses
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) (IAC_a(I),   I = 1,NTOTAT)   ! Integer atom Code
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) (NUMEX(I), I = 1,NTOTAT)   ! Number of atoms excluded from NB interactions 
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) (ICO(I), I = 1,NTYPES*NTYPES) ! index to the nonbonded parameters
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#else
         READ(NF,9118,ERR=666,END=666) (IAC_a(I),   I = 1,NTOTAT)   ! Integer atom Code
         READ(NF,9118,ERR=666,END=666) (NUMEX(I), I = 1,NTOTAT)   ! Number of atoms excluded from NB interactions 
         READ(NF,9118,ERR=666,END=666) (ICO(I), I = 1,NTYPES*NTYPES) ! index to the nonbonded parameters
#endif
         READ(NF,9108,ERR=666,END=666) (LABRES(I),I = 1,NTOTRS)   ! name of residues
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) (IPRES(I), I = 1,NTOTRS)   ! number of first atom in residue I
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#else
         READ(NF,9118,ERR=666,END=666) (IPRES(I), I = 1,NTOTRS)   ! number of first atom in residue I
#endif
         READ(NF,9128,ERR=666,END=666) (RK(I),    I = 1,MUMBND)   ! bond constant
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (REQ(I),   I = 1,MUMBND)   ! bond eq distance
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (TK(I),    I = 1,MUMANG)   ! angle constant
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (TEQ(I),   I = 1,MUMANG)   ! angle Theta equilibrium
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (PK(I),    I = 1,MPTRA)    ! Dihedral force constant
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (PN(I),    I = 1,MPTRA)    ! Dihedral periodicity
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (PHASE(I), I = 1,MPTRA)    ! Dihedral phase
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
#if defined(_AMBER12)
         READ(NF,9128,ERR=666,END=666) (SCEE(I),  I = 1,MPTRA)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9128,ERR=666,END=666) (SCNB(I),  I = 1,MPTRA)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         WRITE(6,'(''AMBER-12 ONLY: SCEE(1),SCNB(1) = '',2E16.8)') 
     .         SCEE(1),SCNB(1)
#endif
         READ(NF,9128,ERR=666,END=666) (SOLTY(I), I = 1,NATYP)    ! Atom type hydrophylicity (NOT USED)
         NTTYP   = (NTYPES*(NTYPES+1))/2          ! Number of non bonded types
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (CN1(I),   I = 1,NTTYP)    ! L-J coefficient 6     
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (CN2(I),   I = 1,NTTYP)    ! L-J coefficient 12      
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666)
     &        (IBH_a(I),JBH_a(I),ICBH_a(I),I = 1,NBONH_a) ! parameters for H-bonds
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666)
     &        (IB_a(I),JB_a(I),ICB_a(I),   I = 1,MBONA) ! parameters for non-H-bonds 
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) 
     c   (ITH_a(I),JTH_a(I),KTH_a(I),ICTH_a(I),I = 1,NTHETH) ! parameters for angles of H
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) 
     c   (IT_a(I),JT_a(I),KT_a(I),ICT_a(I),    I = 1,MTHETS) ! parameters for angles of non-H
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) 
     c   (IPH_a(I),JPH_a(I),KPH_a(I),LPH_a(I),ICPH_a(I),I = 1,NPHIH_a) ! parameters for diedrals of H
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) 
     c   (IP_a(I),JP_a(I),KP_a(I),LP_a(I),ICP_a(I),     I = 1,MPHIA) ! parameters for dihedrals of non-H
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,9117,ERR=666,END=666) (NATEX(I), I = 1,NEXT)  ! Vector of excluded atoms
#else
         READ(NF,9118,ERR=666,END=666)
     &        (IBH_a(I),JBH_a(I),ICBH_a(I),I = 1,NBONH_a) ! parameters for H-bonds
         READ(NF,9118,ERR=666,END=666)
     &        (IB_a(I),JB_a(I),ICB_a(I),   I = 1,MBONA) ! parameters for non-H-bonds
         READ(NF,9118,ERR=666,END=666) 
     c   (ITH_a(I),JTH_a(I),KTH_a(I),ICTH_a(I),I = 1,NTHETH) ! parameters for angles of H
         READ(NF,9118,ERR=666,END=666) 
     c   (IT_a(I),JT_a(I),KT_a(I),ICT_a(I),    I = 1,MTHETS) ! parameters for angles of non-H
         READ(NF,9118,ERR=666,END=666) 
     c   (IPH_a(I),JPH_a(I),KPH_a(I),LPH_a(I),ICPH_a(I),I = 1,NPHIH_a) ! parameters for diedrals of H
         READ(NF,9118,ERR=666,END=666) 
     c   (IP_a(I),JP_a(I),KP_a(I),LP_a(I),ICP_a(I),     I = 1,MPHIA) ! parameters for dihedrals of non-H
         READ(NF,9118,ERR=666,END=666) (NATEX(I), I = 1,NEXT)  ! Vector of excluded atoms
#endif
!
! End of the useful stuff, the rest is not strictly necessary
!         
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (AG(I),    I = 1,NHB)   ! parametrers for the Hydrogen bond interaction
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (BG(I),    I = 1,NHB)   ! parametrers for the Hydrogen bond interaction
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9128,ERR=666,END=666) (HBCUT(I), I = 1,NHB)   ! parametrers for the Hydrogen bond interaction
#if defined(_AMBER7) || defined(_AMBER12)
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
         READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
         READ(NF,9108,ERR=666,END=666) (ISYMBL(I), I=1,NTOTAT)    ! connectivity informations
         IFBOX=0
         if(IFBOX.gt.0)then
           if(IFBOX.ne.1)stop 'only cubic cells'
         
           ALLOCATE(ITREE(NTOTAT))
           ALLOCATE(JOIN(NTOTAT))
           ALLOCATE(IROTAT(NTOTAT))
           ALLOCATE(NSP(NSPM))
#if defined(_AMBER7) || defined(_AMBER12)
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#endif
           READ(NF,9108,ERR=666,END=666) (ITREE(I),  I=1,NTOTAT)
#if defined(_AMBER7) || defined(_AMBER12)
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,9117,ERR=666,END=666) (JOIN(I),   I=1,NTOTAT)
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,9117,ERR=666,END=666) (IROTAT(I), I=1,NTOTAT)    
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,9117,ERR=666,END=666) IPTRES,NSPM,NSPSOL      ! Periodic bondary conditions parameters
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,9117,ERR=666,END=666) (NSP(I),I=1,NSPM)       ! total number of atoms in each molecule
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
           READ(NF,FMT=*,ERR=666,END=666)  ! cmb
#else
           READ(NF,9118,ERR=666,END=666) (JOIN(I),   I=1,NTOTAT)
           READ(NF,9118,ERR=666,END=666) (IROTAT(I), I=1,NTOTAT)    
           READ(NF,9118,ERR=666,END=666) IPTRES,NSPM,NSPSOL      ! Periodic bondary conditions parameters
           READ(NF,9118,ERR=666,END=666) (NSP(I),I=1,NSPM)       ! total number of atoms in each molecule
#endif
           READ(NF,9128,ERR=666,END=666) beta_a,(BOX_a(I),I=1,3) ! cell angle and cell size in Angstrom             
           if(DABS(beta_a-90.).gt.1.d-3) stop 'only cubic box'
           do i=1,3
             BOX(i)=0.1d0*BOX_a(i) ! conversion picometers/Angstroms
           enddo
         else
           NTB=0
           do i=1,3
             BOX(i)=0.d0
           enddo
         endif
         IF(IFCAP.gt.0) THEN
#if defined(_AMBER7) || defined(_AMBER12)
          READ(NF,FMT=*,ERR=666,END=666)  ! cmb
          READ(NF,FMT=*,ERR=666,END=666)  ! cmb
          READ(NF,9117,ERR=666,END=666) NATCAP                  ! Water cap parameters
#else
          READ(NF,9118,ERR=666,END=666) NATCAP                  ! Water cap parameters
#endif
          READ(NF,9128,ERR=666,END=666) CUTCAP,XCAP,YCAP,ZCAP   ! Water cap parameters
         ENDIF
!
! Here there can be other stuff if ones wishes to use Polarizability or to
! perform Free energy perturbation calculations
! These parameters are not currently transferred to the Gromos Topology
!
! *********************************************************************************************
!
!  HERE BEGINS THE CONVERSION OF THE ARRAYS IN GROMOS FORMAT
!
! *********************************************************************************************

         TOPTIT=ITITL
CAK  we add one atom type for the dummy atom potential.
         NRATT=NTYPES+1
         NRATT2=NRATT*(NRATT+1)/2
         NBONH=NBONH_a
         NBON=NBONA
         NTHEH=NTHETH
         NTHE=NTHETS
         NRAA2=NTOTRS
         NSM=0
         NRAM=1
         if(flexwater) then
           NRP=NTOTAT
         else
CAK search for water solvent (residue name WAT)
CAK it has to be at the end of the topology, so we search
CAK backwards until we hit the first non-water residue.
           do i=NTOTRS,1,-1
             if(LABRES(I).eq.'WAT') then
               NSM=NSM+1
               NRAM=3
             else
CAK          no more water, end loop.
               if(NSM.gt.0) then
                 NRAA2=I
CAK          we assume, that the water is kept rigid
CAK          with three distance constraints.
                 NCONS=3
                 ICONS(1)=1
                 JCONS(1)=2
                 ICONS(2)=1
                 JCONS(2)=3
                 ICONS(3)=2
                 JCONS(3)=3
                 goto 333
               endif
             endif
           end do
 333       continue
         endif
         NSM_r=NSM
         NRAM_r=NRAM
         NRP=NTOTAT-(NSM*NRAM) ! subtract solvent from solute
         write(6,*) 'NRP   =',NRP
         write(6,*) 'NSM   =',NSM
         NPTY=MPTRA
         NQTY=0   ! impropers are treated as standard dihedrals in Amber (we hope...)
         conv_cg=1.d0/18.2223d0
         conv_mas=1.d0

!$omp parallel do private(i)
         do i=1,NRP+NSM*NRAM
           PANM(i)=IGRAPH(I)
           CG(i)=conv_cg*CHRG(I)
           WMAS(i)=conv_mas*AMASS(I)
           WINV(i)=1.d0/WMAS(i)
           IAC(i)=IAC_a(I)
         enddo
!$omp parallel do private(i)
         do i=1,NRAA2
           AANM(i)=LABRES(I)
         enddo
! atom names
!$omp parallel do private(i)
        do i=1,MAXATT
          FFTYPE(I)='DUM'
        enddo
!$omp parallel do private(i)
        do i=1,NRP
          FFTYPE(IAC(I))=ISYMBL(I)
        enddo
! add atom names for solvent
        if(NSM.gt.0) THEN
          do i=NRP+1,NRP+NRAM
            FFTYPE(IAC(I))=ISYMBL(I)
          enddo
        endif
! residue number
        J=1
        do i=1,NRP
          if(J.ne.NTOTRS)then
            if(i.ge.IPRES(J+1)) J=J+1
          endif
          MRES(I)=J
        enddo

!    bonds

         conv_cb=fconst_rescale*(kcal_to_kj*(1/aa_to_nm)**2) !kcal/mol/A^2  to kj/mol/nm^2
         NBTY=MUMBND
!$omp parallel do private(i)
         do i=1,NBTY
           B0(i)=aa_to_nm*REQ(I)
           CB(i)=0.5d0*conv_cb*RK(i)/B0(i)**2
         enddo
         NBONH=0
         do i=1,NBONH_a
           IBH(i)=IBH_a(i)/3+1
           JBH(i)=JBH_a(i)/3+1
           ICBH(I)=ICBH_a(I)
C   pack NBONH list
           if((IBH(i).le.NRP).and.(JBH(i).le.NRP)) then
             NBONH=NBONH+1
             IBH(NBONH)=IBH(i)
             JBH(NBONH)=JBH(i)
             ICBH(NBONH)=ICBH(i)
           endif
C   pick water contraints for parameters of first water molecule
           do j=1,3
             if (IBH(I).eq.(NRP+ICONS(j))
     &          .and.JBH(I).eq.(NRP+JCONS(j))) then
               CONS(J)=B0(ICBH_a(I))
             endif
             if (IBH(I).eq.(NRP+JCONS(j))
     &          .and.JBH(I).eq.(NRP+ICONS(j))) then
               CONS(J)=B0(ICBH_a(I))
             endif
           enddo
         enddo
         NBON=0
         do i=1,MBONA
           IB(i)=IB_a(i)/3+1
           JB(i)=JB_a(i)/3+1
           ICB(I)=ICB_a(I)
C   pack NBON list
           if((IB(i).le.NRP).and.(JB(i).le.NRP)) then
             NBON=NBON+1
             IB(NBON)=IB(i)
             JB(NBON)=JB(i)
             ICB(NBON)=ICB(i)
           endif
         enddo
!  angles
 
         conv_t0=1.d0
         to_rad=DACOS(-1.d0)/180.d0
         twoekbt=2.0d0*ekbt
         NTTY=MUMANG
         do i=1,NTTY
           T0(i)=TEQ(I)/to_rad
           TK(i)=fconst_rescale*TK(i)*kcal_to_kj
           if(TK(i).gt.0.1d-4)then
             CT(i)=twoekbt/
     c          ((DCOS(TEQ(I)+DSQRT(ekbt/TK(i)))-DCOS(TEQ(I)))**2+
     c           (DCOS(TEQ(I)-DSQRT(ekbt/TK(i)))-DCOS(TEQ(I)))**2)
           else
             CT(i) = 0.0d0
           end if
         enddo
         NTHEH=NTHETH
!$omp parallel do private(i)
         do i=1,NTHEH
           ITH(i)=ITH_a(i)/3+1
           JTH(i)=JTH_a(i)/3+1
           KTH(i)=KTH_a(i)/3+1
           ICTH(I)=ICTH_a(I)
         enddo
         NTHE=MTHETS
!$omp parallel do private(i)
         do i=1,NTHE
           IT(i)=IT_a(i)/3+1
           JT(i)=JT_a(i)/3+1
           KT(i)=KT_a(i)/3+1
           ICT(I)=ICT_a(I)
         enddo
CAK we also need to remove solvent water angles from the parameters (if present).
CAK since we assume above, that solvent waters are added in a continous block
CAK from the bottom, we can search for them here backwards, too.
CAK i.e. any references atom indices higher than NRP have to be water.
CAK in fact, _all_ three indices have to be larger than NRP or none.
CAK otherwise we have an inconsistent situation and stop with an error.
         NTHETH=0
         do i=NTHEH,1,-1
           if ((ITH(i).gt.NRP).or.(JTH(i).gt.NRP)
     &           .or.(KTH(i).gt.NRP)) then
             if ((ITH(i).gt.NRP).and.(JTH(i).gt.NRP)
     &           .and.(KTH(i).gt.NRP)) then
               NTHETH=i
             else
               write(6,*)'ERROR:angle ',i,' includes solvent and solute'
               stop 'INCONSISTENT Amber topology file. Please check.'
             endif
           endif
         enddo
CAK NTHETH is now the first solven water angle potential (or zero).
         if (NTHETH.gt.0) then
           write(6,*) 'Solvent angle parameters detected'
           write(6,*) 'Adapting NTHEH from ', NTHEH, ' to ', NTHETH-1
           NTHEH=NTHETH-1
         endif

! dihedrals
         pi=DACOS(-1.d0)
         do i =1, NPTY
           CP(I) = kcal_to_kj*PK(I)
           if(DABS(PHASE(I)-pi).lt.1.d-4)then
             PD(I)=-1.d0
           else if(DABS(PHASE(I)).lt.1.d-8)then
             PD(I)= 1.d0
           else
             write(6,*)'not allowed phase for dihedral',i,PHASE(I),pi
             stop  'only 180 and 0 allowed'
           endif
           NP(I) = nint(DABS(PN(I)))
         end do  
         NQHIH=0
         NPHIH=0
         do I=1,NPHIH_a
           NPHIH=NPHIH+1
           IPH(NPHIH)=IPH_a(I)/3+1
           JPH(NPHIH)=JPH_a(I)/3+1
           KPH(NPHIH)=ABS(KPH_a(I))/3+1
           LPH(NPHIH)=ABS(LPH_a(I))/3+1
           ICPH(NPHIH)=ICPH_a(I)
         enddo
         NQHI=0
         NPHI=0
         do I=1,NPHIA
           NPHI=NPHI+1
           IP(NPHI)=IP_a(I)/3+1
           JP(NPHI)=JP_a(I)/3+1
           KP(NPHI)=ABS(KP_a(I))/3+1
           LP(NPHI)=ABS(LP_a(I))/3+1
           ICP(NPHI)=ICP_a(I)
           if(LP_a(I).lt.0)then
             if(NP(ICP(NPHI)).ne.2 .and. PD(ICP(NPHI)).ne.-1.d0)then
               write(6,*)I,NP(ICP(NPHI)), PD(ICP(NPHI))
               stop 'improper dihedral with wrong phase or multiplicity'
             endif
           endif
         enddo                                     
! lennard-jones
         conv_c12=4.184d0*10.d0**(-12)
         conv_c6=4.184d0*10.d0**(-6)
         fac_14=1.0d0/1.2d0
         fac_14vdw=0.5d0
!$omp parallel do private(i)
         do i=1,NTTYP
           C12(I)=conv_c12*CN1(I)
           C6(I)=conv_c6*CN2(I)
           CS12(I)=fac_14vdw*C12(I)
           CS6(I)=fac_14vdw*C6(I)
         enddo
! add LJ parameters for dummy atoms
         NTYPES_d = NTYPES+1
         NTTYP_d  = (NTYPES_d*(NTYPES_d+1))/2
!$omp parallel do private(i)
         do i=NTTYP+1,NTTYP_d
           C12(I)=0.0d0
           C6(I)=0.0d0
           CS12(I)=0.0d0
           CS6(I)=0.0d0
         end do

! exclusions
!$omp parallel do private(i)
         do i=1,NRP
           INE(i)=NUMEX(I)
         enddo
!$omp parallel do private(i)
         do i=1,NEXT
           JSNE(i)=NATEX(I)
         enddo
         allocate(EXCL14(NRP,20))
         do i=1,NRP
           INE14(i)=0
           do j=1,20
             EXCL14(i,j)=99999999
           enddo
         enddo

         do I=1,NPHIA
           if(KP_a(I).ge.0.and.LP_a(I).ge.0) THEN
             i1=IP_a(I)/3+1
             i2=LP_a(I)/3+1
             if(i2 .gt. i1) then
               INE14(i1)=INE14(i1)+1
               EXCL14(i1,INE14(i1))=i2
             else
               INE14(i2)=INE14(i2)+1
               EXCL14(i2,INE14(i2))=i1
             endif
           END IF
         enddo

         do I=1,NPHIH_a
           if(KPH_a(I).ge.0.and.LPH_a(I).ge.0) THEN
             i1=IPH_a(I)/3+1
             i2=LPH_a(I)/3+1
             if(i2 .gt. i1) then
               INE14(i1)=INE14(i1)+1
               EXCL14(i1,INE14(i1))=i2
             else
               INE14(i2)=INE14(i2)+1
               EXCL14(i2,INE14(i2))=i1
             endif
           END IF
         enddo

!  sort EXCL14
         do i=1,NRP
           do j=1,INE14(i)
             sav(j)=EXCL14(i,j)
           enddo
           do jf=1,INE14(i)
             jsmall=0
             smaller=99999999
             do j=1,INE14(i)
               if(sav(j).lt.smaller)then
                 smaller=sav(j)
                 jsmall=j
               endif
             enddo
             EXCL14(i,jf)=sav(jsmall)
             sav(jsmall)=99999999
           enddo
         enddo

        allocate(EXCL(NRP,40))
        KNE(1)=0
        i1=0
        do i=1,NRP
          KNE(i+1)=KNE(i)+INE(i)
          do j=1,INE(i)
            EXCL(i,j)=JSNE(KNE(i)+j)
            sav(j)=EXCL(i,j)
          enddo
          do jf=1,INE(i)
            do j=1,INE14(i)
              if(EXCL14(I,j).eq.sav(jf))sav(jf)=-1
            enddo
            if(EXCL(I,jf).eq.0)then
              if(jf.eq.1)then
                sav(jf)=-1
              else
                stop 'please contact kollmann for further discussions.'
C               hmmm, please suggest a suitable medium!
              endif
            endif
          enddo
          jf=0
          do j=1,INE(i)
            if(sav(j).ne.-1)then
              jf=jf+1
              EXCL(i,jf)=sav(j)
            endif
          enddo
          INE(i)=jf
          do j=1,INE(i)
            i1=i1+1
            JSNE(i1)=EXCL(i,j)
          enddo
        enddo

         i1=0
         do i=1,NRP
           do j=1,INE14(i)
             i1=i1+1
             JSNE14(i1)=EXCL14(i,j)
           enddo
         enddo
!
         NAEX=NEXT
         NAEX14=i1
!
         KNE(1)=0
         do i=1,NRP-1
           KNE(i+1)=KNE(i)+INE(i)
         enddo
!
         KNE14(1)=0
         do i=1,NRP-1
           KNE14(i+1)=KNE14(i)+INE14(i)
         enddo

! charge groups

         max_bond = 8

         allocate(n_bond(NRP))
         allocate(bond(NRP,MAX_BOND))
         allocate(stai_con(NRP))
         allocate(ch_grp(NRP))

         tol_ch_1 = 0.1d0
         tol_ch_2 = 0.7d0
 
         do i=1,NRP
           n_bond(i)=0
           do j=1,MAX_BOND
             bond(i,j)=0
           enddo
         enddo
         do i=1,NBONH
           i1=IBH(i)
           i2=JBH(i)
C make sure we keep only bonds from the solute.
           if((i1.le.NRP).and.(i2.le.NRP)) then
             n_bond(i1)=n_bond(i1)+1
             n_bond(i2)=n_bond(i2)+1
             bond(i1,n_bond(i1))=i2
             bond(i2,n_bond(i2))=i1
           endif
         enddo
         do i=1,NBON
           i1=IB(i)
           i2=JB(i)
C make sure we keep only bonds from the solute.
           if((i1.le.NRP).and.(i2.le.NRP)) then
             n_bond(i1)=n_bond(i1)+1
             n_bond(i2)=n_bond(i2)+1
             bond(i1,n_bond(i1))=i2
             bond(i2,n_bond(i2))=i1
           endif
         enddo

         do i=1,NRP
           stai_con(i)=i
           if(n_bond(i).eq.1) stai_con(i)=bond(i,1)
         enddo

         do i=2,NRP
           if(stai_con(i).lt.stai_con(i-1)) then
             write(6,*) 'WARNING: 
     c       Non contiguous terminal atom bound: ', i, stai_con(i)
           end if         
         end do

!$omp parallel do private(i)
         do i=1,NRP
           ch_grp(i) = 0.0d0 
         enddo                  

!$omp parallel do private(i)
         do i=1,NRP
           ch_grp(stai_con(i)) = ch_grp(stai_con(i)) + CG(i)
         enddo
 
         do i=1,NRP
           if(stai_con(i).ne.i) cycle
           ch1=ch_grp(i)
           d_ch1=DABS(ch1-ANINT(ch1))
           if(d_ch1.lt.tol_ch_1)cycle
           do ic=1,n_bond(i)
             ibbb=bond(i,ic)
             ch2=ch_grp(ibbb)
             cht=ch1+ch2
             d_ch12=DABS(cht-ANINT(cht))
             if(d_ch12/d_ch1.lt.tol_ch_2)then
               if(i.gt.ibbb)then
                 ch_grp(ibbb)=cht
                 ch_grp(i)=0.d0
                 do j=1,NRP
                   if(stai_con(j).eq.i)stai_con(j)=ibbb
                 enddo
               endif
               if(i.lt.ibbb)then
                 ch_grp(ibbb)=0.d0
                 ch_grp(i)=cht
                 do j=1,NRP
                   if(stai_con(j).eq.ibbb)stai_con(j)=i
                 enddo
               endif
             endif
           enddo
         enddo

        NCAG=0
        do i=2,NRP              ! handle flexible water
          if(aanm(mres(i)).eq.'WAT') then
            stai_con(i)=i
            if(mres(i).eq.mres(i-1))stai_con(i)=stai_con(i-1)
            ch_grp(i)=0.d0
          endif
        enddo
        do i=2,NRP
          if(stai_con(i).ne.stai_con(i-1)) then
            NCAG=NCAG+1
            inc(NCAG)=i-1
          endif
        enddo
        if(inc(NCAG).ne.NRP) then
          NCAG=NCAG+1
          inc(NCAG)=NRP
        endif

!        do i=1,4000
!          write(6,333)i,panm(i),aanm(mres(i)),ch_grp(i),stai_con(i)
!        enddo
!333     format(i6,a6,1x,a6,f10.5,1x,i6)

! hard coded variables
        NTPLNS=1
! solvent
        if(NSM.gt.0) then
          do i=1,NRAM
            IACS(i)=IAC_a(NRP+I)
            ANMS(i)=IGRAPH(NRP+I)
            CGS(i)=conv_cg*CHRG(NRP+I)
            WMASS(i)=conv_mas*AMASS(NRP+I)
            WINVS(i)=1.d0/WMAS(NRP+I)
         enddo
        else
          IACS(1)=1
          NRAM=1
          CGS(1)=0.d0
          WINVS(1)=1.d0
          WMASS(1)=1.d0
          ANMS(1)='???'
          NCONS=0
        endif
!
        FPEPSI=0.1389354E+03
        HBAR=0.6350780E-01
!
        NPIA=0
        NPID=0
!$omp parallel do private(i)
        do i=1,NRP
          IPIC(I)=0
        enddo
        call wrtopo(10)
        close(10)
        close(12)
 9108 FORMAT(20A4)
#if defined(_AMBER7) || defined(_AMBER12)
 9117 FORMAT(10I8)  ! new Amber7 format
#else
 9118 FORMAT(12I6)  ! old Amber format
#endif
 9128 FORMAT(5E16.8) 
      RETURN

 666  CONTINUE
      PRINT*,'PROBLEM READING AMBER DATA. ABORTING...'
      STOP
      END
C     ==================================================================
      SUBROUTINE amber_to_gromos_coord(amber_crd,gromos_crd,box_r,NTB_r,
     $     readpdb)
C     ==--------------------------------------------------------------==
      implicit none
C     arguments
      character*(*) amber_crd,gromos_crd
      real*8 box_r(3)
      integer  NTB_r
      logical readpdb
C     local variables
      character*80 aline, dummy
      real*8 maxx,maxy,maxz,minx,miny,minz,ang_r(3)
      integer NAT_l,i
      include 'md.h'
      include 'coordsz.h'
      include 'coordar.h'
      include 'toposz.h'
      include 'topoar.h' 

      LOGICAL LREDTO
      COMMON /TOPOK/LREDTO

      LREDTO=.true.
      NTB_r=0
C     initialize box/angle settings to zero.
      do i=1,3
        box_r(i)=0.d0
        ang_r(i)=0.d0
        v(i)=0.d0
        v(i+3)=0.d0
      enddo

      open(10,file=amber_crd,status='old')

C  PBB format for records that we recognize. version 2.3 (1998).
C
C     The CRYST1 record presents the unit cell parameters, space group,
C     and Z value. If the structure was not determined by
C     crystallographic means, CRYST1 simply defines a unit cube.
C
C Record Format
C
C COLUMNS      DATA TYPE            FIELD        DEFINITION
C ----------------------------------------------------------
C  1 -  6      Record name          "CRYST1"
C  7 - 15      Real(9.3)            a            a (Angstroms).
C 16 - 24      Real(9.3)            b            b (Angstroms).
C 25 - 33      Real(9.3)            c            c (Angstroms).
C 34 - 40      Real(7.2)            alpha        alpha (degrees).
C 41 - 47      Real(7.2)            beta         beta (degrees).
C 48 - 54      Real(7.2)            gamma        gamma (degrees).
C 56 - 66      LString              sGroup       Space group.
C 67 - 70      Integer              z            Z value.

C     The ATOM records present the atomic coordinates for standard
C     residues (see
C     http://deposit.pdb.org/public-component-erf.cif). They also
C     present the occupancy and temperature factor for each
C     atom. Heterogen coordinates use the HETATM record type. The
C     element symbol is always present on each ATOM record; segment
C     identifier and charge are optional.  Record Format
C
C COLUMNS      DATA TYPE        FIELD      DEFINITION
C ------------------------------------------------------
C  1 -  6      Record name      "ATOM    "
C  7 - 11      Integer          serial     Atom serial number.
C 13 - 16      Atom             name       Atom name.
C 17           Character        altLoc     Alternate location indicator.
C 18 - 20      Residue name     resName    Residue name.
C 22           Character        chainID    Chain identifier.
C 23 - 26      Integer          resSeq     Residue sequence number.
C 27           AChar            iCode      Code for insertion of residues.
C 31 - 38      Real(8.3)        x          Orthogonal coordinates for X in
C                                          Angstroms
C 39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in
C                                          Angstroms
C 47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in
C                                          Angstroms
C 55 - 60      Real(6.2)        occupancy  Occupancy.
C 61 - 66      Real(6.2)        tempFactor Temperature factor.
C 77 - 78      LString(2)       element    Element symbol, right-justified.
C 79 - 80      LString(2)       charge     Charge on the atom.

C     The HETATM records present the atomic coordinate records for atoms
C     within "non-standard" groups. These records are used for water
C     molecules and atoms presented in HET groups (see
C     http://deposit.pdb.org/public-component-erf.cif).
C Record Format
C
C COLUMNS     DATA TYPE        FIELD         DEFINITION
C --------------------------------------------------------------
C  1 - 6      Record name      "HETATM"
C  7 - 11     Integer          serial        Atom serial number.
C 13 - 16     Atom             name          Atom name.
C 17          Character        altLoc        Alternate location indicator.
C 18 - 20     Residue name     resName       Residue name.
C 22          Character        chainID       Chain identifier.
C 23 - 26     Integer          resSeq        Residue sequence number.
C 27          AChar            iCode         Code for insertion of residues.
C 31 - 38     Real(8.3)        x             Orthogonal coordinates for X.
C 39 - 46     Real(8.3)        y             Orthogonal coordinates for Y.
C 47 - 54     Real(8.3)        z             Orthogonal coordinates for Z.
C 55 - 60     Real(6.2)        occupancy     Occupancy.
C 61 - 66     Real(6.2)        tempFactor    Temperature factor.
C 77 - 78     LString(2)       element       Element symbol; right-justified.
C 79 - 80     LString(2)       charge        Charge on the atom.
C

      IF (readpdb) THEN
C     parse coordinate and box segment(s) from .pdb file
        NAT_l=0
 99     CONTINUE
        READ(10,'(A80)',END=110,ERR=100) aline
        IF (len(aline).LT.3) GOTO 99
        IF (aline(1:3).EQ.'END') GOTO 110
        IF (len(aline).LT.4) goto 99
        IF (aline(1:4).EQ.'ATOM') THEN
 101      FORMAT (30X,3F8.3)
          READ(aline,101) x(3*NAT_l+1),x(3*NAT_l+2),x(3*NAT_l+3)
          NAT_l=NAT_l+1
        ENDIF
        IF (len(aline).lt.6) goto 99
        IF (aline(1:6).eq.'HETATM') THEN
          READ(aline,101) x(3*NAT_l+1),x(3*NAT_l+2),x(3*NAT_l+3)
          NAT_l=NAT_l+1
        ENDIF
        IF (aline(1:6).eq.'CRYST1') THEN
          READ(aline,'(6X,3F9.3,3F7.2)')
     $        box_r(1),box_r(2),box_r(3),ang_r(1),ang_r(2),ang_r(3)
          NTB_r=1
        ENDIF
        GOTO 99
C     error exit from parsing
 100    CONTINUE
        stop 'error reading .pdb file'
C     regular exit from parsing
 110    CONTINUE
      else
C     parse coordinate and box segment(s) from amber .crd file
        READ(10,*)
        READ(10,*) NAT_l
        if(NAT_l.ne.(NRP+(NSM*NRAM)))
     &       stop 'inconsistent amber topology and coord file'
 113    format(6f12.7)
C     read coordinates and convert to nanometers
        read(10,113,err=116,end=116) (x(i),i=1,3*NAT_l)

C     read velocities, if present. if not v(1-6) will contain the box/angle settings.
        read(10,113,ERR=114,END=114) (v(i),i=1,3*NAT_l)
C     read box and angle in case of a restart.
        read(10,113,err=111,end=111)(box_r(i),i=1,3),(ang_r(i),i=1,3)
        goto 115
 116    continue
        stop 'error reading amber coordinate file'
C     no velocities, copy values back from vel array.
 114    continue
        do i=1,3
          box_r(i)=v(i)
          ang_r(i)=v(3+i)
        enddo

 115    continue
        NTB_r=1
cmb
 111    continue
      endif

      PRINT*,'amber coordinates read successfully'
      DO i=1,3*NAT_l
        x(i)=0.1d0*x(i) ! conversion picometers/Angstrom
      ENDDO
      WRITE(6,*) 'box/angle settings from coordinates file:'
      WRITE(6,*) (box_r(i),i=1,3),(ang_r(i),i=1,3)

cmb
      do i=1,3
        box_r(i)=0.1d0*box_r(i)
      enddo

C     AK 2005/07/24: the input may not have box dimensions.
C     we try to salvage the situation by guessing an orthorhombic
C     box from minmax values of the coordinates.
C     coordinates are already converted to nm.
      if(DABS(box_r(1))+DABS(box_r(2))+DABS(box_r(3)).lt.1.d0) then
        write(6,*)
        write(6,*) 'no reasonable box data in coordinate file.'
        write(6,*) 'trying to guess it. check your files.'
        write(6,*)
        maxx=-9999.d0
        maxy=-9999.d0
        maxz=-9999.d0
        minx= 9999.d0
        miny= 9999.d0
        minz= 9999.d0
        do i=1,3*NAT_l,3
          if(x(i)  .gt.maxx) maxx=x(i)
          if(x(i)  .lt.minx) minx=x(i)
          if(x(i+1).gt.maxy) maxy=x(i)
          if(x(i+1).lt.miny) miny=x(i)
          if(x(i+2).gt.maxz) maxz=x(i)
          if(x(i+2).lt.minz) minz=x(i)
        enddo
        box_r(1)=maxx-minx
        box_r(2)=maxy-miny
        box_r(3)=maxz-minz
        ang_r(1)=90.0d0
        ang_r(2)=90.0d0
        ang_r(3)=90.0d0
        write(6,*) 'box/angle settings from guess:'
        write(6,*) (box_r(i),i=1,3),(ang_r(i),i=1,3)
      endif

      do i=1,3
        if(DABS(ang_r(i)-90.0d0).gt.1.d-6)then
          write(6,*) 'WARNING: angle(',i,')=',ang_r(i) 
          write(6,*) 'only orthorombic cells allowed'
          stop
        endif
      enddo

      open(12,file=gromos_crd,status='unknown')
      write(12,'(a5)')'TITLE'
      write(12,*)'coordinates generated from amber coord file'
      write(12,*)amber_crd
      write(12,'(a3)')'END'
C                  formatted, no4d, ...,NPM,...,NDIM,...
      CALL WRXPOS(12,.true.,.false.,NAT_l,1,NSM,3,X)
      if(NTB_r.eq.1)then
        write(12,'(a3)')'BOX'
        write(12,'(3(f15.8,1x))')(box_r(i),i=1,3)
        write(12,'(a3)')'END'
      endif
      close(10)            
      close(12)            
      RETURN
      END
C     ==================================================================
      SUBROUTINE write_gromos_inp(input_gromos,
     c   ngp_x,ngp_y,ngp_z,NT_UPDATE_NLIST,NAT,
     c   Rcut_alpha,Rcut,box,NTB,NSM,NRAM)
C     ==--------------------------------------------------------------==
      implicit none
      character*30 input_gromos
      integer ngp_x,ngp_y,ngp_z,NT_UPDATE_NLIST,NAT,NTB,NSM,NRAM
      real*8 Rcut_alpha,Rcut,box(3)
 
      open(10,file=input_gromos,status='unknown')

 1010 format ('TITLE',/,A,/,'END')
      write(10,1010)' input generated by amber2gromos'

 1020 format ('SYSTEM',/,
     &     '# NPM >= 0 number of (identical) solute molecules',/,
     &     '# NSM >= 0 number of (identical) solvent molecules',/,
     &     '#      NPM     NSM',/,2I8,/,'END')
      write(10,1020) 1, NSM

 1030 format ('START',/,
     & '# NTX 1..3 reading/generation of initial 3D coordinates (X),',/,
     & '# 3D velocities (V), and stochastic integrals (SX)',/,
     & '# action                                  1    2    3  ',/,
     & '# read X from IOXVI (21)                  yes  yes  yes',/,
     & '# set V to zero if TEMPI = 0.0            yes  no   no ',/,
     & '# read V from IOXVI (21) if  TEMPI = 0.0  no   yes  yes',/,
     & '# read SX from IOXVI (21)                 no   no   yes',/,
     & '# ',/,
     & '# INIT 1..4 startup configuration',/,
     & '# action          1      2    3    4',/,
     & '# shake X        yes     no   no   no',/,
     & '# shake V        yes     yes  no   no',/,
     & '# C.O.M removal  yes     yes  yes  no',/,
     & '# if NTCM =1',/,
     & '#     NTX     INIT      IG     TEMPI      HEAT  NTXO   BOLTZ',/,
     &     3I9,F9.1,F8.5,I8,X,G11.5,/,'END')
      write(10,1030) 1,1, 210185,300.0,0.0,1,8.31441d-3

 1040 format ('STEP',/,
     &     '# NSTLIM > 0 total number of steps',/,
     &     '# T >= 0.0 time at beginning of simulation',/,
     &     '# DT >= 0.0 time step of simulations',/,
     &     '#    NSTLIM         T        DT',/,I10,F10.2,F10.5,/,'END')
      write(10,1040)  10, 0.0, 0.002

 1050 format ('BOUNDARY',/,
     & '# NTB -2..2 define boundary conditions',/,
     & '# <0 truncated octahedron periodic boundary conditions',/,
     & '# >0 rectangular/monoclinic periodic boundary conditions',/,
     & '# 0 vacuum',/,
     & '# abs(NTB) = 2 the virial (pressure) is calculated',/,
     & '# BOX(1..3) > 0.0 periodic box dimensions',/,
     & '# BETA > 0.0 angle of box (inclination between x- ',
     & 'and z-axes) in degrees',/,
     & '# BETA must be 90.0 degrees if abs(NTB) = 2 or if NTB<0',/,
     & '# NRDBOX 0,1 controls reading of BOX dimensions',/,
     & '# 0 use the box dimensions from the BOUNDARY block',/,
     & '# 1 read box dimensions from startup file (BOX) IOXVI',/,
     & '#    NTB   BOX(1)       BOX(2)       BOX(3)',
     & '         BETA   NRDBOX',/,
     &  4X,I4,3(F12.7,X),F10.2,I4,/,'END')
      write(10,1050) NTB,box(1),box(2),box(3), 90.0, 0

 1060 format ('SUBMOLECULES',/,
     & '# NSPM > 0 number of (sub)molecules the solute consists of',/,
     & '# NSP(1.. NSPM) > 0 atom sequence number in topology',/,
     & '# of last atom in submolecule I',/,
     & '#      NSPM  NSP(1.. NSPM)',/,2I8,/,'END')
      write(10,1060) 1,NAT-NSM*NRAM

 1070 format ('TCOUPLE',/,
     & '# NTT(1) -3..3 controls temperature coupling of solute',/,
     & '# internal and rotational degrees of freedom ( d.o.f.)',
     & ' in three dimensions',/,
     & '# NTT(2) -3..3 controls temperature coupling of solute ',/,
     & '# centre of mass translational d.o.f. in three dimensions',/,
     & '# NTT(3) -3..3 controls temperature coupling of ',
     & 'solvent d.o.f.',/,'# in three dimensions',/,
     & '# 0 no temperature coupling for set of d.o.f.',/,
     & '# 1 couple one set of d.o.f. to one bath',/,
     & '# 2,-2 couple two sets of d.o.f. to one bath (+ sign)',/,
     & '# -3,3 couple three sets of d.o.f. to one bath (+ sign)',/,
     & '# TEMP0 >= 0.0 bath reference temperature',/,
     & '# TAUT >= 0.0 coupling time ',/,
     & '#      NTT     TEMP0      TAUT',/,
     &  3(I10,F10.1,F10.3,/),'END') 
      write(10,1070) 0, 300.0, 0.1,  0, 300.0, 0.1, 0, 300.0, 0.1 ! cmb

 1080 format ('CENTREOFMASS',/,
     & '# NDFMIN >= 0 number of degrees of freedom to subtract ',
     & 'from total',/,'# when calculating the system temperature',/,
     & '# NTCM 0,1 controls initial centre of mass motion removal',/,
     & '# 0 no initial centre of mass motion removal',/,
     & '# 1 inital centre of mass motion is removed',/,
     & '# NSCM >= 0 controls centre of mass motion removal during',/,
     & '# simulation',/,
     & '# 0 no centre of mass motion removal',/,
     & '# > 0 centre of mass motion removal every NSCM steps',/,
     & '#    NDFMIN      NTCM      NSCM',/,
     &  3I10,/,'END')
      write(10,1080) 0,0,100000

 1090 format ('PRINT',/,
     &     '# NTPR: print out energies, etc. every NTPR steps',/,
     &     '# NTPL: print out C.O.M motion and total energy fitting',
     &     ' every NTPL steps',/,
     &     '# NTPP: =1 perform dihedral angle transition monitoring',/,
     &     '#      NTPR      NTPL      NTPP',/,3I10,/,'END')
      write(10,1090) 1000,100,0

 1100 format ('WRITE',/,
     &     '# NTPW = 0 : binary',/,
     &     '# NTPW = 1 : formatted',/,
     &     '# NTWSE = configuration selection parameter',/,
     &     '# =0: write normal trajectory',/,
     &     '# >0: chose min energy for writing configurations',/,
     &     '#     NTWX     NTWSE      NTWV      NTWE      NTWG',
     &     '      NTPW',/,6I10,/,'END')
      write(10,1100) 20,0,5000,1,0,1

 1110 format ('SHAKE',/,
     &     '# NTC 1,2,3 controls use of SHAKE for the solute',/,
     &     '# 1 no shake is performed (solute)',/,
     &     '# 2 SHAKE the bonds involving hydrogens (solute)',/,
     &     '# 3 SHAKE all bonds (solute)',/,
     &     '# TOL > 0.0 SHAKE tolerance (relative geometric ',/,
     &     '# precision), for solute and solvent',/,
     &     '#      NTC       TOL',/,I10,F14.7,/,'END')
      write(10,1110) 1,0.0000010
      
 1120 format ('FORCE',/,
     &     '# NTF(I)=0 do not use term I',/,
     &     '# NTF(I)=1 use term I',/,
     &     '# NTF(1) bonds involving hydrogens',/,
     &     '# NTF(2) bonds not involving hydrogens',/,
     &     '# NTF(3) bond angles involving hydrogens',/,
     &     '# NTF(4) bond angles not involving hydrogens',/,
     &     '# NTF(5) improper dihedrals involving hydrogens',/,
     &     '# NTF(6) improper dihedrals not involving hydrogens',/,
     &     '# NTF(7) dihedrals involving hydrogens',/,
     &     '# NTF(8) dihedrals not involving hydrogens',/,
     &     '# NTF(9) nonbonded charge interactions',/,
     &     '# NTF(10) nonbonded interactions',/,
     &     '# NEGR >= 0 number of energy groups',/,
     &     '# NRE > 0 last atoms in each energy group',/,
     &     '# NRE values must be in ascending order',/,
     &     '#  bonds    angles    imp.     dihe     charge nonbonded',/,
     &     '#  H        H         H        H',/,
     &     2X,2I3,3X,2I3,3X,2I3,3X,2I3,3X,2I3,/,
     &     '#       NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)',/,
     &     4X,2I8,/,'END')
      write(10,1120) 1,1,1,1,1,1,1,1,1,1,1,NAT

 1130 format ('PLIST',/,
     &     '# NTNB 0,1 controls pairlist construction before the',
     &     ' first step',/,
     &     '# 0 no pairlist is constructed before the first step',/,
     &     '# 1 a pairlist is constructed before the first step',/,
     &     '# NSNB > 0 controls frequency (number of steps) a ',
     &     'pairlist is',/,'#     constructed',/,
     &     '# RCUTP > 0.0 cut-off used in pairlist construction',/,
     &     '# RCUTL > 0.0 cut-off used in long range interaction',/,
     &     '#      NTNB      NSNB     RCUTP     RCUTL',/,
     &     3X,2I10,2F10.2,/,'END')
      write(10,1130) 1,NT_UPDATE_NLIST,Rcut,Rcut

 1140 format ('LONGRANGE',/,
     &     '# EPSRF =0.0,or >= 1.0 eps used in reaction-field calc.',/,
     &     '# APPAK >= 0.0 kappa used in reaction-field calc.',/,
     &     '# RCRF > 0.0 cut-off radius used in reaction-field calc.',/,
     &     '#     EPSRF     APPAK      RCRF',/,
     &     2F12.1,E12.1,/,'END')
      write(10,1140)50.0 ,0.0,0.7d10

 1150 format ('LATSUM',/,
     &     '# 32 is the ewald grid, 64 is a more conservative choice',/,
     &     '# 0.7 is the size of the head function for ewald summation'
     &     ,/,4I6,I8,2F12.5,I12,/,'END')
     &     
      write(10,1150) 2,ngp_x,ngp_y,ngp_z,0,Rcut_alpha,1.33,100000

 1160 format ('POSREST',/,
     &     '# values for NTR',/,
     &     '# 0: no position re(con)straining',/,
     &     '# 1: use CHO',/,
     &     '# 2: use CHO/ ATOMIC B-FACTORS',/,
     &     '# 3: position constraining',/,
     &     '#      NTR       CHO     NRDRX',/,
     &     3I10,/,'END')
      write(10,1160) 0,1000,0

      close(10)
      RETURN
      END
