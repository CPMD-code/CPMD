C $Id: box.h,v 1.1 2006-12-27 11:15:19 itavern Exp $ -*-fortran-*-

C---------------------------------------------------
C remember to include coordsz.h before this file
C---------------------------------------------------

COMMVAR BOXW
C     the box size in the fourth dimension is
C     set to this value when running in 4D.
C     It can be any large value.
      REAL BOXW
      PARAMETER (BOXW = 1000.0)
COMMEND

COMMVAR NDRMIN,NDRMAX,NDHMIN,NDHMAX,NDIM
      INTEGER NDRMIN,NDRMAX
      PARAMETER (NDRMIN = 1, NDRMAX = 3)

C     NDHMIN and NDHMAX are used to loop over the
C     "harmonic degrees of freedom", and the "real degrees of freedom".
C     In the current implementation the former is the fourth dim,
C     and the real degrees of freedom are dims 1..3.
C     
C     Simulation in three dims:
C     NDIM = 3: NDRMIN = 1, NDRMAX = 3,NDHMIN = 4, NHDMAX = 3
C     Simulation in four dims:
C     NDIM = 4: NDRMIN = 1, NDRMAX = 3,NDHMIN = 4, NDHMAX = 4

C     L<NDO4D> determines the number of dimensions a given term is
C     calculated in, i.e. if NDO4D(N4DNBD) = 4, then the nonbonded
C     calculations are performed in 4D.
COMMEND


      INTEGER NTB,NDIM,NDHMIN,NDHMAX
      COMMON /MBOUND/ NTB,NDIM,NDHMIN,NDHMAX


COMMVBAR BOX,BOXH,BOXINV,BETA
C     BOX   : the dimensions of the periodic box
C     BOXH  : half the box dimensions
C     BOXINV: the inverse box lengths
C     BETA  : the angle of the monoclinic box in degrees
C     See also <LNTB>,L<MAXDIM>,<NDIM>

      REAL BOX,BOXH,BOXINV,BETA
      COMMON /MBNDFP/BOX(MAXDIM),BOXH(MAXDIM),
     $     BOXINV(MAXDIM),BETA
COMMEND


COMMVAR NTB,NTBVAC,NTBVIR,NTBMIN,NTBMAX
C     NTB < 0 : PERIODICITY IS APPLIED
C               BOX IS TRUNCATED OCTAHEDRON (BETA=90)
C         = 0 : NO PERIODICITY IS APPLIED
C         > 0 : PERIODICITY IS APPLIED
C               BOX IS RECTANGULAR OR MONOCLINIC, DEPENDING ON BETA
C     IF ABS(NTB)=2, THE VIRIAL IS CALCULATED (BETA=90)

      INTEGER NTBVAC
      PARAMETER (NTBVAC = 0)
C   remember to compare ABSOLUTE value of NTB
      INTEGER NTBVIR
      PARAMETER (NTBVIR = 2)

      INTEGER NTBMIN,NTBMAX
      PARAMETER (NTBMIN  = -NTBVIR, NTBMAX = NTBVIR)
COMMEND
