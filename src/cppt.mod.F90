MODULE cppt
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! ==           DYNAMIC ALLOCATION OF PERMANENT ARRAYS             ==
  ! ==================================================================
  ! == TSHELS = FALSE if all TSHEL(IS) false                        ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: tshels,tshel(maxsp)
  ! ==================================================================
  ! == INYH(3,NHG) coordinates +NHI (I=1,2,3) for G-vectors         ==
  ! == IGL(NHG)    the index of the shell for each G-vector         ==
  ! == NZH(NHG)    index in PSI or RHO array (IG1>=0)               ==
  ! == INDZ(NHG)   index in PSI or RHO array (IG1<=0)               ==
  ! == ISPTR(NHGL+1) last index IG for the shell                    ==
  ! == INDZS(NGW)  index for G-compon. of wavefunction in PSI (I1<0)==
  ! == NZH(NGW)    index for G-compon. of wavefunction in PSI (I1>0)==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, TARGET :: nzh(:)
  INTEGER, ALLOCATABLE, TARGET :: nzhs(:)
  INTEGER, ALLOCATABLE, TARGET :: indz(:)
  INTEGER, ALLOCATABLE, TARGET :: indzs(:)
  INTEGER, ALLOCATABLE, TARGET :: inyh(:,:)

  INTEGER, ALLOCATABLE :: igl(:)
  INTEGER, ALLOCATABLE :: isptr(:)

  ! ==================================================================
  ! == HG(1:NHG) Square norm of G-vectors                           ==
  ! == GK(1:3,1:NHG) Components of G-vectors                        ==
  ! == GL(1:NHGL) Square norm of G-vectors for each shell           ==
  ! == VPS(1:NSP,1:NHG) Local pseudopotential per species in G space==
  ! == RHOPS: smeared ionic charge density in G space               ==
  ! ==        ionic point charges replaced by Gaussian charge       ==
  ! ==        distributions (see ener.inc) calculated in PUTPS      ==
  ! == TWNL(1:NGW,1:NGH(IS),1:NSP) Non-Local projectors array       ==
  ! ==        for each G-components (Kleinman-Bylander form)        ==
  ! == USED BY VANDERBILT PSEUDOPOTENTIALS                          ==
  ! == QRAD                                                         ==
  ! == YLMB(NHGK,LPMAX,NKPNT) Initialized in PUTWNL                 ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: hg(:)
  REAL(real_8), ALLOCATABLE :: gk(:,:)
  REAL(real_8), ALLOCATABLE :: gl(:)
  REAL(real_8), ALLOCATABLE :: vps(:,:)
  REAL(real_8), ALLOCATABLE :: rhops(:,:)

  REAL(real_8), ALLOCATABLE :: twnl(:,:,:,:) ! (:)??(:,:,:)??

  REAL(real_8), ALLOCATABLE :: qrad(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: twnls(:,:,:)
  REAL(real_8), ALLOCATABLE :: ylmb(:,:,:)


  ! ==================================================================
  ! == Dimension of HGPOT (for isolated system -- HIP)              ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nr1h,nr2h,nr3h,nr3pl
  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: hgpot(:,:,:)
  REAL(real_8), ALLOCATABLE :: hipz(:)

  COMPLEX(real_8), ALLOCATABLE, TARGET :: scg(:)
  COMPLEX(real_8), POINTER :: scgx(:)



  ! ==================================================================

END MODULE cppt
