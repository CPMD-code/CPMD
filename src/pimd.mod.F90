MODULE pimd
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == PATH INTEGRAL PARAMETERS                                     ==
  ! ==--------------------------------------------------------------==
  ! == !! NOTE: IF YOU ADD ITEMS TO THE COMMON BLOCKS INCREASE THE  ==
  ! ==          RESPECTIVE PARAMETERS PI_LOG,PI_REL,PI_INT          ==
  ! ==================================================================
  INTEGER, PARAMETER :: pi_log=10
  INTEGER, PARAMETER :: pi_rel=5 
  INTEGER, PARAMETER :: pi_int=4 
  ! ==--------------------------------------------------------------==
  ! == PILOG       Only for broadcasting                            ==
  ! == RREP        See TINIT                                        ==
  ! == TREAD_CENT  Classical coordinates are read in as centroid    ==
  ! ==             positions and normal modes are sampled.          ==
  ! == TSTAGE      Use the staging representation                   ==
  ! ==             for the path integrals                           ==
  ! == TPINM       Use the normal mode representation               ==
  ! ==             for the path integrals                           ==
  ! == TINIT       Initialization run:                              ==
  ! ==             Supply full start configuration either by        ==
  ! ==             (i)  generation -> with RREP=.TRUE. or by        ==
  ! ==             (ii) reading it -> with REPREAD=.TRUE            ==
  ! == REPREAD     Read replica coordinates from external file      ==
  ! ==             to be named 'repfname' in the loop format        ==
  ! ==             "DO IP, DO IS, DO IA ... (K=1,3)"                ==
  ! ==             (see 'rread.F') .with an empty line              ==
  ! ==             Before EVERY set of replicas (also the first!)   ==
  ! == TCENTRO     Do centroid path integral MD and                 ==
  ! ==             do not thermostat centroid mode,                 ==
  ! ==             i.e., first Trotter slice                        ==
  ! == TRINGP      Do ring-polymer path integral MD                 ==
  ! == TESTPI      Classical path integral test                     ==
  ! ==--------------------------------------------------------------==

  ! ==--------------------------------------------------------------==
  ! == PIREL       Only for broadcasting                            ==
  ! == TEMPB       Use TEMPB as the reference temperature (in K)    ==
  ! ==             to generate the De Broglie displaced NP-1 other  ==
  ! ==             particles create free particle paths             ==
  ! ==             using a Gaussian Levy walk at temperature TEMPB  ==
  ! ==             according to the masses of the species           ==
  ! ==             starting from the first bead from the classical  ==
  ! ==             configuration that is obtained from input        ==
  ! == WMASS       Ficticious ionic masses are WMASS*M              ==
  ! ==             No specification of WMASS -> WMASS=1.0           ==
  ! ==             (choice of Berne et al.)                         ==
  ! == FACSTAGE    The ficticious normal mode mass is               ==
  ! ==             (FACSTAGE*WMASS) times heavier than the corres-  ==
  ! ==             -ponding real physical mass of the same species  ==
  ! ==             The fictitious non-centroid masses are           ==
  ! ==             FACSTAGE times lighter than the centroid mass    ==
  ! ==             (which has the physical mass as fictitious mass) ==
  ! == GLE_LAMBDA  Friction coefficient of GLE thermostat is        ==
  ! ==             GLE_LAMBDA times larger than the frequency       ==
  ! ==             of the harmonic potential                        ==
  ! ==--------------------------------------------------------------==

  ! ==--------------------------------------------------------------==
  ! == PIINT    Only for broadcasting                               ==
  ! == NP_TOTAL Number of replica                                   ==
  ! == LEVPRT   Printing level: minimal output for LEVPRT < 5       ==
  ! == LOUTFN   Output files for each processor, group              ==
  ! ==          or only grandparent                                 ==
  ! ==--------------------------------------------------------------==

  ! ==================================================================
  ! == REPFNAME Name of the external file to read replica coord.    ==
  ! ==--------------------------------------------------------------==
  CHARACTER(len=80) :: repfname
  ! ==================================================================
  ! == VARIABLES FOR PARALLEL WORK                                  ==
  ! ==--------------------------------------------------------------==
  ! == SUPERGROUP                                                   ==
  ! == SUPERSOURCE                                                  ==
  ! == PARENTGROUP                                                  ==
  ! == PC_GRP                                                       ==
  ! == PC_GROUPS   Number of processor goups                        ==
  ! == PCG_POS                                                      ==
  ! == NPROC_TOT                                                    ==
  ! == IPCURR                                                       ==
  ! == GRANDPARENT                                                  ==
  ! ==--------------------------------------------------------------==
  INTEGER :: supergroup,supersource,parentgroup
  INTEGER, ALLOCATABLE, DIMENSION(:) :: pc_grp
  INTEGER :: pc_groups,pcg_pos,nproc_tot,ipcurr
  LOGICAL :: grandparent,PARENT_PARENTGROUP
  ! ==================================================================
  ! == MAXNP  Maximum number of replica                             ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: maxnp=64 
  ! ==--------------------------------------------------------------==
  ! == NUM_NP                                                       ==
  ! == IPP1                                                         ==
  ! == IPM1                                                         ==
  ! == NP_LOW                                                       ==
  ! == NP_HIGH                                                      ==
  ! == NP_LOCAL                                                     ==
  ! ==--------------------------------------------------------------==
  INTEGER :: num_np(maxnp),ipp1(maxnp),ipm1(maxnp)
  INTEGER :: np_low,np_high,np_local
  ! ==--------------------------------------------------------------==
  ! == ENERGIES PER REPLICA                                         ==
  ! ==--------------------------------------------------------------==
  ! == EHAR                                                         ==
  ! == ETOTV(MAXNP)                                                 ==
  ! == EKINV(MAXNP)                                                 ==
  ! == EHTV(MAXNP)                                                  ==
  ! == EPSEUV(MAXNP)                                                ==
  ! == ENLV(MAXNP)                                                  ==
  ! == EXCV(MAXNP)                                                  ==
  ! == ESELFV(MAXNP)                                                ==
  ! == ESRV(MAXNP)                                                  ==
  ! == EGCV(MAXNP)                                                  ==
  ! == EEIGV(MAXNP)                                                 ==
  ! == ECNSTRV(MAXNP)                                               ==
  ! == EHARV(MAXNP)    Harmonic energy                              ==
  ! == CSUMGV(MAXNP)                                                ==
  ! == CSUMRV(MAXNP)                                                ==
  ! == VDBCHGV(MAXNP)                                               ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: ehar
  REAL(real_8) :: etotv(maxnp),ekinv(maxnp),ehtv(maxnp),epseuv(maxnp)&
       ,enlv(maxnp),excv(maxnp),eselfv(maxnp),esrv(maxnp),egcv(&
       maxnp),eeigv(maxnp),ecnstrv(maxnp),eharv(maxnp),csumgv(maxnp)&
       ,csumrv(maxnp),vdbchgv(maxnp,maxnp)
  ! ==================================================================
  ! == TREP        Atomic coordinates for each replica              ==
  ! == TREPNM      Normal modes for each replica                    ==
  ! == FIONKS      Atomic force coordinates for each replica        ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: trep(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: trepnm(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: fionks(:,:,:,:)


  ! ==================================================================
  ! == PMAR(MAXSP)                                                  ==
  ! == PMAR0(MAXSP)                                                 ==
  ! == FLNM(MAXNP)                                                  ==
  ! == PMARS(MAXSP,MAXNP)                                           ==
  ! == PMA0S(MAXSP,MAXNP)                                           ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: pmar(maxsp),pmar0(maxsp),flnm(maxnp)
  REAL(real_8) :: pmars(maxsp,maxnp),pma0s(maxsp,maxnp)
  ! ==--------------------------------------------------------------==
  ! == TNM          Forward normal mode transformation matrix       ==
  ! == TNMI                                                         ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: tnm(:,:)
  REAL(real_8), ALLOCATABLE :: tnmi(:,:)


  ! ==================================================================
  ! == FPGYR(MAXSP)                                                 ==
  ! == AVGYR(MAXSP)                                                 ==
  ! == AVGYRA(MAXSP,2)                                              ==
  ! == FPSUS(MAXSP)                                                 ==
  ! == AVSUS(MAXSP)                                                 ==
  ! == AVSUSA(MAXSP,2)                                              ==
  ! == FPCOR(MAXSP)                                                 ==
  ! == AVCOR(MAXSP)                                                 ==
  ! == AVCORA(MAXSP,2)                                              ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: fpgyr(maxsp),avgyr(maxsp),avgyra(maxsp,2),fpsus(&
       maxsp),avsus(maxsp),avsusa(maxsp,2),fpcor(maxsp),avcor(maxsp)&
       ,avcora(maxsp,2)
  ! ==--------------------------------------------------------------==
  ! == RGYR(NAX,NSX)                                                ==
  ! == RGYRA(NAX,NSX,2)                                             ==
  ! == RSUS(NAX,NSX)                                                ==
  ! == RSUSA(NAX,NSX,2)                                             ==
  ! == RCOR(NAX,NSX)                                                ==
  ! == RCORA(NAX,NSX,2)                                             ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: rgyr(:,:)
  REAL(real_8), ALLOCATABLE :: rgyra(:,:,:)
  REAL(real_8), ALLOCATABLE :: rsus(:,:)
  REAL(real_8), ALLOCATABLE :: rsusa(:,:,:)
  REAL(real_8), ALLOCATABLE :: rcor(:,:)
  REAL(real_8), ALLOCATABLE :: rcora(:,:,:)


  ! ==================================================================
  ! == IRCORR(NAX,NSX,NP_TOTAL)                                     ==
  ! == IFCORR(NAX,NSX,NP_TOTAL)                                     ==
  ! == RCORR (NAX,NSX,NP_TOTAL)                                     ==
  ! == FCORR (NAX,NSX,NP_TOTAL)                                     ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: ircorr(:,:,:)
  INTEGER, ALLOCATABLE :: ifcorr(:,:,:)

  REAL(real_8), ALLOCATABLE :: rcorr(:,:,:)
  REAL(real_8), ALLOCATABLE :: fcorr(:,:,:)


  REAL(real_8) :: pi_omega(maxnp)
  REAL(real_8), ALLOCATABLE :: pi_egle(:)
  REAL(real_8), ALLOCATABLE :: pi_glec(:)
  REAL(real_8), ALLOCATABLE :: pi_gles(:,:)
  REAL(real_8), ALLOCATABLE :: pi_glet(:,:)
  REAL(real_8), ALLOCATABLE :: pi_glep(:,:,:,:,:)

  ! ==================================================================

  TYPE :: pimd1_t
     LOGICAL :: pilog
     LOGICAL :: rrep
     LOGICAL :: tread_cent
     LOGICAL :: tstage
     LOGICAL :: tpinm
     LOGICAL :: tinit
     LOGICAL :: repread
     LOGICAL :: tcentro
     LOGICAL :: tringp
     LOGICAL :: testpi
  END TYPE pimd1_t
  TYPE(pimd1_t) :: pimd1
  TYPE :: pimd2_t
     REAL(real_8) :: pirel
     REAL(real_8) :: tempb
     REAL(real_8) :: wmass
     REAL(real_8) :: facstage
     REAL(real_8) :: gle_lambda
  END TYPE pimd2_t
  TYPE(pimd2_t) :: pimd2
  TYPE :: pimd3_t
     INTEGER :: piint
     INTEGER :: np_total
     INTEGER :: levprt
     INTEGER :: loutfn
  END TYPE pimd3_t
  TYPE(pimd3_t) :: pimd3

END MODULE pimd
