MODULE iffi_elstat_utils
  USE efld,                            ONLY: extf,&
                                             textfld
  USE iffi_types                       
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nkpt,  &
       parm, spar, parap
  USE elct,                            ONLY: crge
  USE ragg,                            ONLY: raggio
  USE zeroing_utils,                   ONLY: zeroing
  USE iffi_comm_utils,                 ONLY: SHIFT_MULTIPOLES_WRAP,INHERIT_LOCAL_EXP_WRAP
  USE rhoofr_utils,                    ONLY: rhoofr
  USE coor,                            ONLY: tau0
  USE iffi_wfct
  USE timer, ONLY: tiset, tihalt
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: QMCOREFORCE,CALC_IMPORT,CALC_EXPORT
  
CONTAINS
!....Interface to the MD simulation program      IPHIGENIE
!....Original version by M. Eichinger and J. Hutter
!....rewritten by B. Breitenfeld and M. Schwoerer
!..... INTERFACE_INIT              : initially assign memory
!..... CALC_IMPORT               : calculate electrostatic potential from near atoms, voxels and taylor expansions
!..... OLE_IMPORT                  : build up nearlist and olelist, calculate ole coefficients
!..... CALC_EXPORT                     : calculate electrostatics at near atoms, multipole moments of voxels and qm atoms
!..... OLE_EXPORT                  : calculate electrostatics at the mm atoms in nearlist and olelist
!..... MMCOREFORCE                 : calculate electrostatics at near atoms due to QM cores
!..... QMCOREFORCE                 : calculate force on QM atoms caused by external potential
!***************************************************************************
! 1/sqrt(pi) 
#define M_SQRTPI_1 0.56418958354775628695  


!==================================================================
SUBROUTINE CALC_IMPORT
!==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER                    :: procedureN = 'CALC_IMPORT'
      
      REAL*8  KO_X,KO_Y,KO_Z,DI,DJ,DK
      REAL*8 DELTA 

!     control variables to determine what quantities will be called in this function
      LOGICAL CALC_MM,CALC_PFF                              

      REAL*8  DX,DY,DZ,DX2,DY2,DZ2                              
      REAL*8  DXXX,DYYY,DZZZ,DXXY,DYYX,DZZX,DZZY,DXXZ,DYYZ                              
      REAL*8  XXXX,YYYY,ZZZZ,XXYZ,YYXZ,ZZXY                              
!     Variables
      INTEGER IM,N,VO                    
      REAL*8  NEARPOT_MM,NEARPOT_PFF                            
      REAL*8  D,DERF,ERRORF                              
      REAL*8  pdotr_j,D_1,D_2,arg                               
      REAL*8  SIGMA_1,SIGMA_2,LAMBDA,MU      

#ifdef __DERF                              
      EXTERNAL DERF                              
#endif                              
      REAL*8  CLCPOT_MM, CLCPOT_TE, CLCPOT_PFF
      INTEGER ISUB,I,J,K
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/
      
      INTEGER ikr1
      INTEGER idx !to access extf
!     ==--------------------------------------------------------------==

      IKR1 = MOD(parm%nr1+1,2)   ! !added: biswas

      CALL TISET(procedureN,ISUB)
      
      IF(paral%parent) WRITE (6,'(A)') 'INTERFACE| CALCULATING EXTERNAL POTENTIAL...'

!    zeroing(ECHRG)
      
!.....build cpmd near and ole lists, calculate expansion coefficients at voxels
      CALL OLE_IMPORT

!........decide what to calculate
      IF (runinfo%pff) THEN
          CALC_MM  = .FALSE.
          CALC_PFF = .TRUE.
        IF (.NOT.epot1_var%tcalconlypotfrompff) THEN
          CALC_MM  = .TRUE.
        ENDIF
      ELSE
          CALC_MM  = .TRUE.
          CALC_PFF = .FALSE.
      ENDIF

!     reset EXTFM
      CALL zeroing(EXTF)!,KR1*KR2S*KR3S)
!     reset EXTFMM if this is the first run after an integration step
      IF (CALC_MM) THEN
        CALL zeroing(EXTFMM)!,KR1*KR2S*KR3S)
      ENDIF
      
      
      IF (runinfo%meanfieldmode.AND.IFIRST.EQ.1) THEN
        CALL zeroing(EXTFMEAN)!,KR1*KR2S*KR3S)
        CALL zeroing(EXTFLAST)!,KR1*KR2S*KR3S)
      ENDIF

!.....lattice constant
      DK=(runinfo%boxdum(6)-runinfo%boxdum(5))/DBLE(spar%nr3s)
      DJ=(runinfo%boxdum(4)-runinfo%boxdum(3))/DBLE(spar%nr2s)
      DI=(runinfo%boxdum(2)-runinfo%boxdum(1))/DBLE(spar%nr2s)

!$omp parallel
!$omp do private(KO_X,KO_Y,KO_Z,CLCPOT_MM,CLCPOT_TE,CLCPOT_PFF,DELTA, &
!$omp  DX,DY,DZ,VO,DX2,DY2,DZ2,                                       &
!$omp  DXXX,DYYY,DZZZ,DXXY,DYYX,DZZX,DZZY,DXXZ,DYYZ,                  &
!$omp XXXX,YYYY,ZZZZ,XXYZ,YYXZ,ZZXY,                                  &
!$omp  IM,NEARPOT_MM,NEARPOT_PFF,                                     &
!$omp  D,D_1,ERRORF,pdotr_j,D_2,SIGMA_1,SIGMA_2,ARG,LAMBDA,MU,        &
!$omp  I,J,K,idx)                                                  
!CCCCCCCCCCCCCCCCCCreduction (+:GRIDNEARWW)                           
!.....LOOP OVER GRIDPOINTS
      DO K=1,fpar%kr3s
       KO_Z=runinfo%boxdum(5)+REAL(k-1,kind=real_8)*dk
       DO J=1,fpar%kr2s
        KO_Y=runinfo%boxdum(3)+REAL(j-1,kind=real_8)*dj
        DO I=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)  + ikr1! IKR1 added: biswas
         KO_X=runinfo%boxdum(1)+REAL(i-1,kind=real_8)*di

         idx = (i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s
         
         CLCPOT_MM  = 0.D0
         CLCPOT_TE  = 0.D0
         CLCPOT_PFF = 0.D0
         
         IF (epot1_var%updatemeanfield.AND.epot1_var%tdoextrap) THEN !epot1_var%tdoextrap==FIRSTRUNAFTERINT
!           update mean potential..(from last update step)
            DELTA = EXTFLAST(idx)&
                   - EXTFMEAN(idx)
            EXTFMEAN(idx) = &
                   EXTFMEAN(idx)+ DELTA/epot1_var%nmean
!           ..and write it back later
           !now calculate a new ext pot and save it into EXTFLAST
         ENDIF
        

!.......calculate the potential
!       first evaluate potential from local expansion 

!.......find the voxel VO to which this gridpoint is assigned to
        VO=VOXPART(idx)

!       vector from voxel VO to gridpoint
        DX=KO_X-KOVOX(1,VO)
        DY=KO_Y-KOVOX(2,VO)
        DZ=KO_Z-KOVOX(3,VO)

        DX2=DX*DX
        DY2=DY*DY
        DZ2=DZ*DZ 

! CPMD wants negative potential!
! all expansions coefficients from iffi save the potential carry a minus sign!
        CLCPOT_TE = -VOXEXP(KPHI,VO) + &
                DX*VOXEXP(KX,VO) + &
                DY*VOXEXP(KY,VO) + &
                DZ*VOXEXP(KZ,VO) + &
         0.5D0*DX2*VOXEXP(KXX,VO) + &
         0.5D0*DY2*VOXEXP(KYY,VO) + &
         0.5D0*DZ2*(-VOXEXP(KXX,VO)-VOXEXP(KYY,VO)) + &
             DX*DY*VOXEXP(KXY,VO) + &
             DX*DZ*VOXEXP(KXZ,VO) + &
             DY*DZ*VOXEXP(KYZ,VO) 

!       taylor expansion to third order
!        IF (runinfo%sammtype.GE.4) THEN
          DXXX=DX2*DX
          DYYY=DY2*DY
          DZZZ=DZ2*DZ
          DXXY=DX2*DY
          DYYX=DY2*DX
          DZZX=DZ2*DX
          DZZY=DZ2*DY
          DXXZ=DX2*DZ
          DYYZ=DY2*DZ

          CLCPOT_TE = CLCPOT_TE +       ( &
                   VOXEXP(KXXX,VO) *DXXX+ &
                   VOXEXP(KYYY,VO) *DYYY+ &
                 (-VOXEXP(KZXX,VO) &
                  -VOXEXP(KZYY,VO))*DZZZ+ &
           3.0D0*  VOXEXP(KXYY,VO) *DYYX+ &
           3.0D0*(-VOXEXP(KXXX,VO) &
                  -VOXEXP(KXYY,VO))*DZZX+ &
           3.0D0*  VOXEXP(KYXX,VO) *DXXY+ &
           3.0D0*(-VOXEXP(KYYY,VO) &
                  -VOXEXP(KYXX,VO))*DZZY+ &
           3.0D0*  VOXEXP(KZXX,VO) *DXXZ+ &
           3.0D0*  VOXEXP(KZYY,VO) *DYYZ+ &
           6.0D0*  VOXEXP(KXYZ,VO) *DX*DY*DZ &
                                 ) / 6.0D0
!        ENDIF

!       taylor expansion to forth order
!        IF (runinfo%sammtype.GE.5) THEN
          XXXX = -VOXEXP(KXXYY,VO)-VOXEXP(KXXZZ,VO)
          YYYY = -VOXEXP(KXXYY,VO)-VOXEXP(KYYZZ,VO)
          ZZZZ = -VOXEXP(KXXZZ,VO)-VOXEXP(KYYZZ,VO)
          XXYZ = -VOXEXP(KYYYZ,VO)-VOXEXP(KZZZY,VO)
          YYXZ = -VOXEXP(KXXXZ,VO)-VOXEXP(KZZZX,VO)
          ZZXY = -VOXEXP(KXXXY,VO)-VOXEXP(KYYYX,VO)

          CLCPOT_TE = CLCPOT_TE +    (         & 
                           XXXX    *DX2*DX2 +  &
                           YYYY    *DY2*DY2 +  &
                           ZZZZ    *DZ2*DZ2 +  &
            6.0D0* VOXEXP(KXXYY,VO)*DXXY*DY +  &
            6.0D0* VOXEXP(KXXZZ,VO)*DXXZ*DZ +  &
            6.0D0* VOXEXP(KYYZZ,VO)*DYYZ*DZ +  &
            4.0D0* VOXEXP(KXXXY,VO)*DXXX*DY +  &
            4.0D0* VOXEXP(KXXXZ,VO)*DXXX*DZ +  &
            4.0D0* VOXEXP(KYYYX,VO)*DYYY*DX +  &
            4.0D0* VOXEXP(KYYYZ,VO)*DYYY*DZ +  &
            4.0D0* VOXEXP(KZZZX,VO)*DZZZ*DX +  &
            4.0D0* VOXEXP(KZZZY,VO)*DZZZ*DY +  &
           12.0D0*         XXYZ    *DXXY*DZ +  &
           12.0D0*         YYXZ    *DYYX*DZ +  &
           12.0D0*         ZZXY    *DZZX*DY    &
                          ) / 24.0D0

!.....CLCPOT_TE now carries the contribution from the Taylorexpansion

!.....Now add potential from all NEAR atoms (class C_-3)

!.....loop over all near mmatoms of the current voxel
      DO 200 N=1,NEARLISTLEN(VO)
        IM=NEARLIST(VO,N)
        NEARPOT_MM  = 0.D0
        NEARPOT_PFF = 0.D0
        
!       if runinfo%pff: skip atoms without dipole in all but first import for each step
        IF (CALC_MM.OR.(CALC_PFF.AND.          &
           epot1_fix%neardata(SIGDIPOLE,IM).GT.0.000001D0))  THEN

!            vector from near mmatom to gridpoint
             DX=KO_X-epot1_fix%neardata(N_KX,IM)
             DY=KO_Y-epot1_fix%neardata(N_KY,IM)
             DZ=KO_Z-epot1_fix%neardata(N_KZ,IM)
!
             D = DSQRT(DX*DX + DY*DY + DZ*DZ)

!            important, as mm atoms can be very close to or even on grid points!
             IF (D .LT. 0.1) D = 0.1D0
             
             D_1 = 1.0D0 / D

!+-+-+-+-+-+- POTENTIAL FROM MM CHARGES+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
            IF (CALC_MM) THEN
!.............check for nonzero charge
              IF (epot1_fix%neardata(PCHARGE,IM).GT. 0.0000001D0.OR. &
                 epot1_fix%neardata(PCHARGE,IM).LT.-0.0000001D0) THEN
     
!              decide whether we are in the gauss- or point-region:
               IF (D.GT.(runinfo%olensig*epot1_fix%neardata(SIGCHARGE,IM))) THEN
!.................use point charges
!                 Phi = q / d
                  NEARPOT_MM = epot1_fix%neardata(PCHARGE,IM) * D_1
               ELSE
!.................use gaussian charges
!.                Phi = q * erf(d/sigma) / d
                  ERRORF = DERF(D/epot1_fix%neardata(SIGCHARGE,IM))
                  NEARPOT_MM = epot1_fix%neardata(PCHARGE,IM)*ERRORF*D_1
!              end of if gaussian or point
               ENDIF
!            end of check for nonzero charge
             ENDIF
!           end of if CALC_MM
            ENDIF
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 

!+-+-+-+-+-+- POTENTIAL FROM PMM DIPOLES +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!...........If runinfo%pff is on and this mmatom is polarizable (sigma = epot1_fix%neardata(SIGDIPOLE,IM) > 0.0, then calc dipole potenial)
            IF (runinfo%pff.AND.epot1_fix%neardata(SIGDIPOLE,IM).GT.0.D0) THEN
!             decide whether we are in the gauss- or point-region:
              IF (D.GT.(runinfo%olensig*epot1_fix%neardata(SIGDIPOLE,IM))) THEN
!.............use point dipoles
                 pdotr_j = epot1_var%neardata_var(N_PX,IM)*DX + &
                          epot1_var%neardata_var(N_PY,IM)*DY +  &
                          epot1_var%neardata_var(N_PZ,IM)*DZ
                 NEARPOT_PFF = pdotr_j / (D*D*D)
              ELSE
!.................use gaussian dipoles
                 !formulas taken from appendix B of sebastian bauer's diss
                 D_2 = D_1 * D_1
                 SIGMA_1 = 1.0D0/epot1_fix%neardata(SIGDIPOLE,IM) !sqrt2 implied
                 SIGMA_2 = SIGMA_1 * SIGMA_1
                 ARG=D*SIGMA_1  !sqrt2 implied
                 LAMBDA=(2.0D0 * M_SQRTPI_1 )      &
                         * SIGMA_1 * DEXP(-ARG*ARG) !sqrt2 implied
                 MU=D_1 * DERF(ARG)
                 NEARPOT_PFF = (MU-LAMBDA) * D_2 *       &
                              (epot1_var%neardata_var(N_PX,IM)*DX +   &
                               epot1_var%neardata_var(N_PY,IM)*DY +   &
                               epot1_var%neardata_var(N_PZ,IM)*DZ)
 
!             end of if gaussian or point
              ENDIF
!            end of if runinfo%pff
             ENDIF
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 

!............add NEAR mm and pff contributions to CLCPOT (from other NEAR atoms)                        
             CLCPOT_MM  = CLCPOT_MM  - NEARPOT_MM                   
             CLCPOT_PFF = CLCPOT_PFF - NEARPOT_PFF

!       end of if skip this atom if runinfo%pff and not CALC_MM and not polarizable
        ENDIF

!.....end loop over all near mmatoms
200   ENDDO

!........store the potential in the extf structures
         IF (runinfo%pff) THEN
           IF (CALC_MM) THEN
             EXTFMM(idx) = CLCPOT_MM
           ENDIF

           extf(idx)=EXTFMM(idx)+ &
            CLCPOT_TE + CLCPOT_PFF

         ELSE
            extf(idx)= CLCPOT_TE + CLCPOT_MM
         ENDIF

!........uncomment for printing association gridpoints-voxels to cube file through RHOOUT
!        EXTF(I-parap%nrxpl(parai%mepos,1)+1,J,K)=VO

         IF (epot1_var%updatemeanfield) THEN 
          !save most recent potential
          EXTFLAST(idx)=extf(idx)
          
          IF(epot1_var%nmean.EQ.1)  & !at the very first step do update mean every time
            EXTFMEAN(idx)= EXTFLAST(idx)
          
          !..and use mean pot for wavefunction
          extf(idx)=EXTFMEAN(idx)

         ENDIF

!       end of loop over all grid points
        ENDDO
       ENDDO
      ENDDO
!$omp end do 
!$omp end parallel 
      
      IF(paral%parent) WRITE (6,'(A)') 'INTERFACE| ...EXTERNAL POTENTIAL CALCULATED'
      
      IF(IFIRST.EQ.1) IFIRST=0

      CALL TIHALT(procedureN,ISUB)
!==--------------------------------------------------------------==
 RETURN
 END SUBROUTINE
!==================================================================
!




!==================================================================
SUBROUTINE OLE_IMPORT
!==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER                    :: procedureN = 'OLE_IMPORT'
      
      INTEGER ISUB
      INTEGER VO, IQ, L,J, IM, N, th
      REAL*8 D,DMIN, DX,DY,DZ

      REAL*8 f,f1,f2,f3,DIP_X,DIP_Y,DIP_Z
      REAL*8 dR2eX2,dR2eY2,dR2eZ2,dR2eR2
      REAL*8 distMDXY,distMDXZ,distMDYZ,distMDXYZ
      REAL*8 rpowPR1,rpowPR2,rpowPR3,rpowPR5
      REAL*8 FivePR2,One_FivePR2_X2,One_FivePR2_Y2,tmp35PR4,tmp35PR2
      REAL*8 tmp15PR2,tmpYYY,tmpXXX,tmpX2,tmpY2, tmp45_rpowPR2
      REAL*8 tmpZ2,pYdrY,pXdrX

      LOGICAL OLEALLOWED
!$    INTEGER   OMP_GET_THREAD_NUM
!$    EXTERNAL  OMP_GET_THREAD_NUM

      if (paral%parent.AND.runinfo%tverbose) print *,"ENTER OLE_IMPORT"
      
      CALL TISET(procedureN,ISUB)

!.....reset expansions and lists
      CALL zeroing(NEARLISTLEN)!,MYNVOX)
      CALL zeroing(OLELISTLEN)!,MYNVOX)
      CALL zeroing(VOXEXP)!,dimlocexp*MYNVOX)
      IF (.NOT.runinfo%pff.OR..NOT.epot1_var%tcalconlypotfrompff) THEN
          CALL zeroing(VOXEXP_MM)!,dimlocexp*MYNVOX)
      ENDIF

!.....assign voxels to QM atoms (VO2IQ)
      IF (epot1_var%tgridpart) THEN
       CALL zeroing(VO2IQ)!,MYNVOX)
!     important for non-openmp runs
      th=1
!$omp parallel private(VO,DMIN,IQ,L,DX,DY,DZ,D,th)
!$     th=OMP_GET_THREAD_NUM()+1
!$omp do 
       DO VO=1,MYNVOX
           DMIN   = 999999D0
           IQ     = -1
           DO  L=1,epot1_list%nloc
             DX=KOVOX(1,VO)-epot1_fix%koloc(1,L)
             DY=KOVOX(2,VO)-epot1_fix%koloc(2,L)
             DZ=KOVOX(3,VO)-epot1_fix%koloc(3,L)
             D=DX*DX+DY*DY+DZ*DZ
             IF (D .LT. DMIN) THEN
              DMIN=D
              IQ=L
             END IF
           ENDDO
!           LI2VOX(IQ)=VO
           VO2IQ(VO)=IQ
!          save maximum distance qmatom-voxel
       ENDDO !end loop over voxels
!$omp end do 
!$omp end parallel
      ENDIF
!.....end of assign voxels to QM atoms

!.....VOXEXP_MM (from charges) is calculated only once each integration step, VOXEXP (from dipoles) each cpmd call

!.....now handle coefficients and near/ole for all my voxels
!     no need to thread-parallelize VOXEXP-structure as parallelization is over VO
!$omp parallel private(VO,IQ,DX,DY,DZ,N,IM,D,                             &
!$omp dR2eX2,dR2eY2,dR2eZ2,dR2eR2,                                        &
!$omp OLEALLOWED,distMDXY,distMDXZ,distMDYZ,distMDXYZ,                    &
!$omp rpowPR1,rpowPR2,f,FivePR2,One_FivePR2_X2,One_FivePR2_Y2,tmp35PR2,   &
!$omp tmp35PR4,f1,tmpX2,tmpY2,tmpZ2,                                      &
!$omp DIP_X,DIP_Y,DIP_Z,pXdrX,pYdrY,rpowPR3,rpowPR5,f2,f3,                &
!$omp tmp45_rpowPR2,tmp15PR2,tmpYYY,tmpXXX)
!$omp do 
      DO VO=1,MYNVOX
!.......inherit coefficients from QM atom to voxel
        IQ=VO2IQ(VO)
!       vector from qmatom to voxel
        DX=KOVOX(1,VO)-epot1_fix%koloc(1,IQ)
        DY=KOVOX(2,VO)-epot1_fix%koloc(2,IQ)
        DZ=KOVOX(3,VO)-epot1_fix%koloc(3,IQ)

!.......inherit expansions from level C_-1 (qm atoms) to level C_-2 (voxels)
        CALL INHERIT_LOCAL_EXP_WRAP(VO, &
                        DX,DY,DZ, IQ) !reads epot1_var%locexp writes to VOXEXP
 
!.......for all iffi near atoms (failed C_-1-IAC) of qmatom IQ belonging to voxel VO:
!       decide between near (C_-3) and ole (C_^-2) (and calc coefficients)
        DO N=1,epot1_list%nincl(IQ)
           IM = epot1_list%lincl(IQ,N)
!          vector from nearatom to voxel
           DX=KOVOX(1,VO)-epot1_fix%neardata(N_KX,IM)
           DY=KOVOX(2,VO)-epot1_fix%neardata(N_KY,IM)
           DZ=KOVOX(3,VO)-epot1_fix%neardata(N_KZ,IM)

!          DR2 (vectorElementMultiplication)
           dR2eX2 = DX*DX           
           dR2eY2 = DY*DY           
           dR2eZ2 = DZ*DZ   
           dR2eR2 = dR2eX2+dR2eY2+dR2eZ2 
           D = DSQRT(dR2eR2)

           OLEALLOWED=.FALSE.
! 1st crit: check if the mm object is farer away than n sigma:  d >? nsig * max(sig_q,sig_p) !IAC Eq. (6) in Schwoerer/JCP2015
           IF (D.GT.runinfo%olensig*epot1_fix%neardata(SIGCHARGE,IM)) THEN !1/sqrt(2) implied in runinfo%olensig
             IF(D.GT.runinfo%olensig*epot1_fix%neardata(SIGDIPOLE,IM)) THEN
! 2nd crit: oleradius calculated from voxel size and epsilon
               IF (D.GT.runinfo%oleradius) THEN    !IAC Eq. (5) in Schwoerer/JCP2015
                 OLEALLOWED=.TRUE.
               ENDIF
             ENDIF
           ENDIF

           IF (OLEALLOWED) THEN
!.....helping quantities
!     distm (calcDistanceMatrix)
      IF ((.NOT.runinfo%pff).OR.(.NOT.epot1_var%tcalconlypotfrompff).OR. &
          epot1_fix%neardata(SIGDIPOLE,IM).GT.0.0D0) THEN
      distMDXY = DX * DY 
      distMDXZ = DX * DZ 
      distMDYZ = DY * DZ 
      distMDXYZ= distMDXY* DZ 
!     distance matrix   (calcPowersOfDistance)             (scalefac = 1)
      rpowPR1      =   1.0D0     / D     
      rpowPR2      = rpowPR1 * rpowPR1
      ENDIF

! CHARGE      (setMonopoleContributionsToExpansion)
      IF (.NOT.runinfo%pff.OR..NOT.epot1_var%tcalconlypotfrompff) THEN
       f = epot1_fix%neardata(PCHARGE,IM) * rpowPR1 
       VOXEXP_MM(KPHI,VO) = VOXEXP_MM(KPHI,VO)+f

       f = f * rpowPR2
       VOXEXP_MM(KX,VO) = VOXEXP_MM(KX,VO) + f * dX
       VOXEXP_MM(KY,VO) = VOXEXP_MM(KY,VO) + f * dY
       VOXEXP_MM(KZ,VO) = VOXEXP_MM(KZ,VO) + f * dZ
    
       f1 = -3.0D0 * f * rpowPR2
       VOXEXP_MM(KXX,VO) = VOXEXP_MM(KXX,VO) + (f + dR2eX2 * f1)
       VOXEXP_MM(KYY,VO) = VOXEXP_MM(KYY,VO) + (f + dR2eY2 * f1)
       VOXEXP_MM(KXY,VO) = VOXEXP_MM(KXY,VO) +  distMDXY * f1
       VOXEXP_MM(KXZ,VO) = VOXEXP_MM(KXZ,VO) +  distMDXZ * f1
       VOXEXP_MM(KYZ,VO) = VOXEXP_MM(KYZ,VO) +  distMDYZ * f1

       FivePR2 = 5.0D0 * rpowPR2
       One_FivePR2_X2 =  1.0D0 - FivePR2 * dR2eX2
       One_FivePR2_Y2 =  1.0D0 - FivePR2 * dR2eY2

       VOXEXP_MM(KXXX,VO) = VOXEXP_MM(KXXX,VO) +                   &
                         f1 * ( 3.0D0 - FivePR2 * dR2eX2 ) * dX
       VOXEXP_MM(KYYY,VO) = VOXEXP_MM(KYYY,VO) +                   &
                         f1 * ( 3.0D0 - FivePR2 * dR2eY2 ) * dY
       VOXEXP_MM(KXYY,VO) = VOXEXP_MM(KXYY,VO) +                   &
                         f1 * ( One_FivePR2_Y2 ) * dX
       VOXEXP_MM(KYXX,VO) = VOXEXP_MM(KYXX,VO) +                   &
                         f1 * ( One_FivePR2_X2 ) * dY            
       VOXEXP_MM(KZXX,VO) = VOXEXP_MM(KZXX,VO) +                   &
                         f1 * ( One_FivePR2_X2 ) * dZ
       VOXEXP_MM(KZYY,VO) = VOXEXP_MM(KZYY,VO) +                   &
                         f1 * ( One_FivePR2_Y2 ) * dZ
       VOXEXP_MM(KXYZ,VO) = VOXEXP_MM(KXYZ,VO)                     &
                         -f1 *   FivePR2 * distMDXYZ

       tmp35PR2 = 7.0D0 * FivePR2
       tmp35PR4 = tmp35PR2 * rpowPR2

       VOXEXP_MM(KXXYY,VO) = VOXEXP_MM(KXXYY,VO)                    &
                     + f1 * ( dR2eX2 * dR2eY2 * tmp35PR4            &
                     - ( dR2eX2 + dR2eY2 ) * FivePR2 + 1.0D0 )
       VOXEXP_MM(KXXZZ,VO) = VOXEXP_MM(KXXZZ,VO)                    &
                     + f1 * ( dR2eX2 * dR2eZ2 * tmp35PR4            &
                     - ( dR2eX2 + dR2eZ2 ) * FivePR2 + 1.0D0 )
       VOXEXP_MM(KYYZZ,VO) = VOXEXP_MM(KYYZZ,VO)                    &
                     + f1 * ( dR2eY2 * dR2eZ2 * tmp35PR4            &
                     - ( dR2eY2 + dR2eZ2 ) * FivePR2 + 1.0D0 )

       f1 = f1 * rpowPR2
       tmpX2 = f1 * ( dR2eX2 * tmp35PR2 - 15.0D0)
       tmpY2 = f1 * ( dR2eY2 * tmp35PR2 - 15.0D0)
       tmpZ2 = f1 * ( dR2eZ2 * tmp35PR2 - 15.0D0)

       VOXEXP_MM(KXXXY,VO) = VOXEXP_MM(KXXXY,VO) + tmpX2 * distMDXY
       VOXEXP_MM(KXXXZ,VO) = VOXEXP_MM(KXXXZ,VO) + tmpX2 * distMDXZ
       VOXEXP_MM(KYYYX,VO) = VOXEXP_MM(KYYYX,VO) + tmpY2 * distMDXY
       VOXEXP_MM(KYYYZ,VO) = VOXEXP_MM(KYYYZ,VO) + tmpY2 * distMDYZ
       VOXEXP_MM(KZZZX,VO) = VOXEXP_MM(KZZZX,VO) + tmpZ2 * distMDXZ
       VOXEXP_MM(KZZZY,VO) = VOXEXP_MM(KZZZY,VO) + tmpZ2 * distMDYZ
      ENDIF

! dipole    (addDipoleContributionsToExpansion)
      IF (epot1_fix%neardata(SIGDIPOLE,IM).GT.0.0D0) THEN
       DIP_X =epot1_var%neardata_var(N_PX,IM) 
       DIP_Y =epot1_var%neardata_var(N_PY,IM) 
       DIP_Z =epot1_var%neardata_var(N_PZ,IM) 

       pXdrX = DIP_X*dX
       pYdrY = DIP_Y*dY
       rpowPR3      = rpowPR2 * rpowPR1
       rpowPR5      = rpowPR3 * rpowPR2

       f1 = pXdrX+pYdrY+DIP_Z*DZ

       VOXEXP(KPHI,VO) = VOXEXP(KPHI,VO) + f1 * rpowPR3

       f2 = 3.0D0 * f1
       f3 = f2 * rpowPR2
       VOXEXP(KX,VO) = VOXEXP(KX,VO) + (f3 * dX -DIP_X) *rpowPR3
       VOXEXP(KY,VO) = VOXEXP(KY,VO) + (f3 * dY -DIP_Y) *rpowPR3
       VOXEXP(KZ,VO) = VOXEXP(KZ,VO) + (f3 * dZ -DIP_Z) *rpowPR3

       f3 = f3 * 5.0D0

       VOXEXP(KXX,VO) = VOXEXP(KXX,VO) +                        &
                 (f2 + 6.0D0 * pXdrX - dR2eX2*f3)*rpowPR5
       VOXEXP(KYY,VO) = VOXEXP(KYY,VO) +                        &
                 (f2 + 6.0D0 * pYdrY - dR2eY2*f3)*rpowPR5
       VOXEXP(KXY,VO) = VOXEXP(KXY,VO) +                        &
                 (3.0D0 * (DIP_X*dY + DIP_Y*dX)                 &
                     - distMDXY * f3) * rpowPR5
       VOXEXP(KXZ,VO) = VOXEXP(KXZ,VO) +                        &
                 (3.0D0 * (DIP_X*dZ + DIP_Z*dX)                 &
                    - distMDXZ * f3) * rpowPR5
       VOXEXP(KYZ,VO) = VOXEXP(KYZ,VO) +                        &
                 (3.0D0 * (DIP_Y*dZ + DIP_Z*dY)                 &
                    - distMDYZ * f3) * rpowPR5

       tmp45_rpowPR2 = 45.0D0 * rpowPR2
       f3 = f3 * 7.0D0

       VOXEXP(KXYZ,VO) = VOXEXP(KXYZ,VO) +                       &
                ( f3 * distMDXYZ - 15.0D0 * ( DIP_X * distMDYZ   &
                    + DIP_Y * distMDXZ + DIP_Z * distMDXY)  )    &
                    * rpowPR5 * rpowPR2 

            f3 = f3* rpowPR2

       VOXEXP(KXXX,VO) = VOXEXP(KXXX,VO) +                       &
                 ( ( f3 * dR2eX2                                 &
                      - tmp45_rpowPR2 * ( pXdrX + f1 ) ) * dX    &
                     + 9.0D0 * DIP_X ) * rpowPR5 
       VOXEXP(KYYY,VO) = VOXEXP(KYYY,VO) +                       &
                 ( ( f3 * dR2eY2                                 &
                      - tmp45_rpowPR2 * ( pYdrY + f1 ) ) * dY    &
                     + 9.0D0 * DIP_Y ) * rpowPR5 

       f1 = f1 *  15.0D0

       tmp15PR2 = 15.0 * rpowPR2
       tmpYYY =   f3 * dR2eY2 - ( 30.0D0 * pYdrY + f1 ) * rpowPR2
       tmpXXX =   f3 * dR2eX2 - ( 30.0D0 * pXdrX + f1 ) * rpowPR2
       tmpX2  = dR2eX2 * tmp15PR2 - 3.0D0 
       tmpY2  = dR2eY2 * tmp15PR2 - 3.0D0 

       VOXEXP(KXYY,VO) = VOXEXP(KXYY,VO)                        &
                + (tmpYYY * dX - tmpY2 * DIP_X ) * rpowPR5
       VOXEXP(KYXX,VO) = VOXEXP(KYXX,VO)                        &
                + (tmpXXX * dY - tmpX2 * DIP_Y ) * rpowPR5
       VOXEXP(KZXX,VO) = VOXEXP(KZXX,VO)                        &
                + (tmpXXX * dZ - tmpX2 * DIP_Z ) * rpowPR5
       VOXEXP(KZYY,VO) = VOXEXP(KZYY,VO)                        &
                + (tmpYYY * dZ - tmpY2 * DIP_Z ) * rpowPR5

      ENDIF !end if dipole
      
!           save atom in ole list
            OLELISTLEN(VO) = OLELISTLEN(VO) + 1
            OLELIST(VO,OLELISTLEN(VO)) = IM
            IF (OLELISTLEN(VO).GT.maxolepervo)         & 
               CALL STOPGM('OLE_IMPORT','Increase maxolepervo!',&
                   __LINE__,__FILE__)

      ELSE !not OLE
!           save atom in near list
            NEARLISTLEN(VO) = NEARLISTLEN(VO) + 1
            NEARLIST(VO,NEARLISTLEN(VO)) = IM
            IF (NEARLISTLEN(VO).GT.maxnearpervo)       &
               CALL STOPGM('OLE_IMPORT','Increase maxnearpervo!',&
                   __LINE__,__FILE__)

           ENDIF
         ENDDO !loop over mm atoms
      ENDDO  !loop over voxels
!$omp end do 
!$omp end parallel 

!$omp parallel private(VO,J)
!$omp do 
      DO VO=1,MYNVOX
        DO J=1,dimlocexp
          VOXEXP(J,VO) = VOXEXP(J,VO) + VOXEXP_MM(J,VO)
        ENDDO
      ENDDO
!$omp end do 
!$omp end parallel 

      CALL TIHALT(procedureN,ISUB)

 END SUBROUTINE
!==================================================================



!==================================================================
SUBROUTINE CALC_EXPORT
!==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER                    :: procedureN = 'CALC_EXPORT'
!     Variables      
      REAL*8  KO_X,KO_Y,KO_Z,DI,DJ,DK
      REAL*8  D,D_1,D_2,D_3
      REAL*8  DX,DY,DZ,DX2,DY2,DZ2,DXY,DXZ,DYZ
      
      REAL*8  sigmaT,neg_inv_2_sigmaT_sqr,inv_sigmaT, &
             prefac_sigmaT
      REAL*8  factor1,factor2
      REAL*8  FktGauss,FktErf_D
      REAL*8  Q,QFAC,Q3,Q5
      REAL*8  RRQ,RR
      INTEGER ISUB,I,J,K,VO,IQ,A,IM,WARN

      REAL*8 MU, LAMBDA, FAC1, FAC2, SIGMA_1, SIGMA_2, ARG
      LOGICAL CALC_MM
      REAL*8 ONEOVERNUMGRIDPOINTS,NUMGRIDPOINTS
      REAL*8 WEIGHT

      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/
      
      REAL*8  DERF
#ifdef __DERF
      EXTERNAL DERF
#endif

      INTEGER ikr1,ir
      INTEGER idx !to access extf
            
!     openmp variables
      INTEGER th,mom
!$    INTEGER   OMP_GET_THREAD_NUM
!$    EXTERNAL  OMP_GET_THREAD_NUM

!!
      IF (paral%io_parent) WRITE (6,'(3A)') 'INTERFACE| CALCULATING ELECTROSTATICS AT MM-ATOMS...'

!     (re)calculate RHOE
      CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
   
     
      CALL TISET(procedureN,ISUB)

!     zero potential expansion ad MM atoms
      CALL zeroing(MMPOTEXP)!,dimpotexp*maxnear*mxno_threads)  
      CALL zeroing(MMPOTEXP_VAR)!,DIMPOTEXP_VAR*maxnear*mxno_threads)  
      IF (paral%parent) THEN
        CALL zeroing(MMPOTEXPGLOB)!,dimpotexp*maxnear)
        CALL zeroing(MMPOTEXPGLOB_VAR)!,DIMPOTEXP_VAR*maxnear)
      ENDIF
      
!     zero QM multipole moments
      CALL zeroing(VOMULTEXP)!,dimmultexp*MYNVOX*mxno_threads)
      CALL zeroing(ATMULTEXP)!,dimmultexp*maxind*mxno_threads)
      IF (paral%parent) CALL zeroing(ATMULTEXPGLOB)!,dimmultexp*maxind) !there VOMULTEXPGLOB-structure

!     zero QM rgyr moments
      IF (epot1_var%tupdatergyr) THEN
        CALL zeroing(VOXRGYREXP)!,dimrgyrexp*MYNVOX*mxno_threads)
        CALL zeroing(ATRGYREXP)!, dimrgyrexp*maxind*mxno_threads)
        IF (paral%parent) CALL zeroing(ATRGYREXPGLOB)!,dimrgyrexp*maxind)
      ENDIF


!.....lattice constant
      DK=(runinfo%boxdum(6)-runinfo%boxdum(5))/DBLE(spar%nr3s)
      DJ=(runinfo%boxdum(4)-runinfo%boxdum(3))/DBLE(spar%nr2s)
      DI=(runinfo%boxdum(2)-runinfo%boxdum(1))/DBLE(spar%nr2s)

      WARN=0
      QFAC=-parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)     
      
      NUMGRIDPOINTS=DBLE(spar%nr3s*spar%nr2s*spar%nr2s)
      ONEOVERNUMGRIDPOINTS=1.0D0/NUMGRIDPOINTS 

      ir = 0
      IKR1 = MOD(parm%nr1+1,2)   ! !added: biswas
      
       
!     important for non-openmp runs
      th=1
      
      if (.NOT.runinfo%pff.OR.(.NOT.epot1_var%tcalconlypfffield)) then
        CALC_MM = .TRUE.
      else
        CALC_MM = .FALSE.
      endif

!$omp parallel private(I,J,K,                                &
!$omp  KO_Z,KO_Y,KO_X,VO,IQ,Q,WEIGHT,                        &
!$omp  Q3,Q5,DX,DY,DZ,D,DX2,DY2,DZ2,RR,RRQ,DXY,DXZ,DYZ,      &
!$omp  A,IM,D_1,D_2,D_3,factor1,sigmaT,neg_inv_2_sigmaT_sqr, &
!$omp  prefac_sigmaT,FktGauss,FktErf_D,inv_sigmaT,           &
!$omp  factor2,SIGMA_1,SIGMA_2,ARG,LAMBDA,MU,FAC1,FAC2,th,idx)
!$      th=OMP_GET_THREAD_NUM()+1
!$omp do 
!.....LOOP OVER GRIDPOINTS
      DO K=1,fpar%kr3s
       KO_Z=runinfo%boxdum(5)+REAL(k-1,kind=real_8)*dk
       DO J=1,fpar%kr2s
         KO_Y=runinfo%boxdum(3)+REAL(j-1,kind=real_8)*dj
         DO I=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)  + ikr1! IKR1 added: biswas
           KO_X=runinfo%boxdum(1)+REAL(i-1,kind=real_8)*di
     
           idx = (i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s
! +++++++++ THE SAMM PART of CALC_EXPORT +++++++++++++++++++++++++++++
!..........calculate multipole expansions of voxels

!..........look up the  voxel associated to the current grid point
           VO = VOXPART(idx)

!..........get grid charge
           Q=RHOE(idx,clsd%nlsd)*QFAC

           Q3=3.D0*Q
           Q5=5.D0*Q
           
!..........handle gyration tensors
           IF (epot1_var%tupdatergyr.OR.runinfo%meanfieldmode) THEN
             WEIGHT=1.0D0
             WEIGHT=-Q/crge%nel !number of electrons   !Eq. (8) in Schwoerer/JCP2015 
             IF (WEIGHT.GT.ONEOVERNUMGRIDPOINTS ) THEN    !scaling factor c would be here  
               WEIGHT=1.0D0
             ELSE
               WEIGHT=WEIGHT*NUMGRIDPOINTS                !scaling factor c would be here           
             ENDIF
             
             IF (epot1_var%tupdatergyr) THEN
!..............calculate contribution to gyration moments (normalization is done later)
               VOXRGYREXP(NGAMMA,VO,th)= VOXRGYREXP(NGAMMA,VO,th)+1.0D0
                
            VOXRGYREXP(RGZERO,VO,th)= VOXRGYREXP(RGZERO,VO,th) + WEIGHT
            VOXRGYREXP(RGX,VO,th) = VOXRGYREXP(RGX,VO,th) + WEIGHT*KO_X
            VOXRGYREXP(RGY,VO,th) = VOXRGYREXP(RGY,VO,th) + WEIGHT*KO_Y
            VOXRGYREXP(RGZ,VO,th) = VOXRGYREXP(RGZ,VO,th) + WEIGHT*KO_Z
               VOXRGYREXP(RGXX,VO,th) = VOXRGYREXP(RGXX,VO,th)        &
                                               + WEIGHT*KO_X*KO_X
               VOXRGYREXP(RGYY,VO,th) = VOXRGYREXP(RGYY,VO,th)        &
                                               + WEIGHT*KO_Y*KO_Y
               VOXRGYREXP(RGZZ,VO,th) = VOXRGYREXP(RGZZ,VO,th)        &
                                               + WEIGHT*KO_Z*KO_Z
               VOXRGYREXP(RGXY,VO,th) = VOXRGYREXP(RGXY,VO,th)        &
                                               + WEIGHT*KO_X*KO_Y
               VOXRGYREXP(RGXZ,VO,th) = VOXRGYREXP(RGXZ,VO,th)        &
                                               + WEIGHT*KO_X*KO_Z
               VOXRGYREXP(RGYZ,VO,th) = VOXRGYREXP(RGYZ,VO,th)        &
                                               + WEIGHT*KO_Y*KO_Z
               ENDIF
            ENDIF
            
!..........vector from voxel (=refpoint) to gridpoint
           DX=KO_X-KOVOX(1,VO)
           DY=KO_Y-KOVOX(2,VO)
           DZ=KO_Z-KOVOX(3,VO)
           
           DX2=DX*DX
           DY2=DY*DY
           DZ2=DZ*DZ
           DXY=DX*DY
           DXZ=DX*DZ
           DYZ=DY*DZ   
           
           RR=DX2+DY2+DZ2
           RRQ=Q*RR

!..........now calculate and add up all multipole moments
!          charge
           VOMULTEXP(LAD,VO,th)=VOMULTEXP(LAD,VO,th)+Q

!          dipolemoment 
           VOMULTEXP(PX,VO,th)=VOMULTEXP(PX,VO,th)+Q*DX
           VOMULTEXP(PY,VO,th)=VOMULTEXP(PY,VO,th)+Q*DY
           VOMULTEXP(PZ,VO,th)=VOMULTEXP(PZ,VO,th)+Q*DZ

!          quadrupolemoment 
           VOMULTEXP(QXX,VO,th)=VOMULTEXP(QXX,VO,th)+(Q3*DX2-RRQ)
           VOMULTEXP(QYY,VO,th)=VOMULTEXP(QYY,VO,th)+(Q3*DY2-RRQ)
           VOMULTEXP(QXY,VO,th)=VOMULTEXP(QXY,VO,th)+(Q3*DXY)
           VOMULTEXP(QXZ,VO,th)=VOMULTEXP(QXZ,VO,th)+(Q3*DXZ)
           VOMULTEXP(QYZ,VO,th)=VOMULTEXP(QYZ,VO,th)+(Q3*DYZ)

!          octupolemoment
!           IF (runinfo%sammtype.GE.4) THEN
           VOMULTEXP(OXYZ,VO,th)=VOMULTEXP(OXYZ,VO,th)     &
            +Q5*DXY*DZ
            VOMULTEXP(OXXY,VO,th)=VOMULTEXP(OXXY,VO,th)    &
                                 +(Q5*DX2-RRQ)*DY
           VOMULTEXP(OXXZ,VO,th)=VOMULTEXP(OXXZ,VO,th)     &
                                 +(Q5*DX2-RRQ)*DZ
           VOMULTEXP(OYYX,VO,th)=VOMULTEXP(OYYX,VO,th)     &
                                 +(Q5*DY2-RRQ)*DX
           VOMULTEXP(OYYZ,VO,th)=VOMULTEXP(OYYZ,VO,th)     &
                                 +(Q5*DY2-RRQ)*DZ
           VOMULTEXP(OZZX,VO,th)=VOMULTEXP(OZZX,VO,th)     &
                                 +(Q5*DZ2-RRQ)*DX
           VOMULTEXP(OZZY,VO,th)=VOMULTEXP(OZZY,VO,th)     &
                                 +(Q5*DZ2-RRQ)*DY
!           ENDIF
!          hexadecapolemoment
!           IF (runinfo%sammtype.GE.5) THEN
            VOMULTEXP(HXXYY,VO,th) = VOMULTEXP(HXXYY,VO,th)+       &
                                 Q5*(7.0D0* DX2*DY2 -              &
                                       RR*(DX2+DY2)) + RR*RRQ
            VOMULTEXP(HXXZZ,VO,th) = VOMULTEXP(HXXZZ,VO,th)+       &
                                 Q5*(7.0D0*DX2*DZ2 -               &
                                       RR*(DX2+DZ2)) + RR*RRQ
            VOMULTEXP(HYYZZ,VO,th) = VOMULTEXP(HYYZZ,VO,th)+       &
                                 Q5*(7.0D0*DY2*DZ2-                &
                                       RR*(DY2+DZ2)) + RR*RRQ
            VOMULTEXP(HXXXY,VO,th) = VOMULTEXP(HXXXY,VO,th)+       &
                                 Q5*( DXY*(7.0D0*DX2-3.0D0*RR) )
            VOMULTEXP(HXXXZ,VO,th) = VOMULTEXP(HXXXZ,VO,th)+       &
                                 Q5*DXZ*(7.0D0*DX2-3.0D0*RR) 
            VOMULTEXP(HYYYX,VO,th) = VOMULTEXP(HYYYX,VO,th)+       &
                                 Q5*DXY*(7.0D0*DY2-3.0D0*RR) 
            VOMULTEXP(HYYYZ,VO,th) = VOMULTEXP(HYYYZ,VO,th)+       &
                                 Q5* DYZ*(7.0D0*DY2-3.0D0*RR)         
            VOMULTEXP(HZZZX,VO,th) = VOMULTEXP(HZZZX,VO,th)+       &
                                 Q5*DXZ*(7.0D0*DZ2-3.0D0*RR) 
            VOMULTEXP(HZZZY,VO,th) = VOMULTEXP(HZZZY,VO,th)+       &
                                 Q5* DYZ*(7.0D0*DZ2-3.0D0*RR) 
!           ENDIF
!


! +++++++++ THE NEAR PART of CALC_EXPORT +++++++++++++++++++++++++++++

!..........for near mm atoms calculate potential and derivatives generated by current gridpoint

!..........loop over all NEAR mmatoms of voxel VO
           DO A=1,NEARLISTLEN(VO)
!           atomid IM of nearatom A 
            IM=NEARLIST(VO,A)
             
!           vector from gridpoint to near mm-atom
            DX=epot1_fix%neardata(N_KX,IM)-KO_X
            DY=epot1_fix%neardata(N_KY,IM)-KO_Y
            DZ=epot1_fix%neardata(N_KZ,IM)-KO_Z
            DX2=DX*DX
            DY2=DY*DY
            DZ2=DZ*DZ
!           Distance
            D = DSQRT(DX2+DY2+DZ2)
!           important, as mm atoms can be very close to or even on grid points!
            IF (D .LT. 0.1) D = 0.1D0

            D_1 = 1.0D0/D
            D_2 = D_1 * D_1
            D_3 = D_2 * D_1

! +++++++++ MM REGION of CALC_EXPORT++++++++++++++++++
!...........note: all gaussian widths carry a factor of sqrt(2)
!           potential and field at charges      
            IF (CALC_MM) THEN
!            decide whether we are in the gauss- or point-region
             IF (D.GT.(runinfo%olensig*epot1_fix%neardata(SIGCHARGE,IM))) THEN
!              use point charges
!              potential  Phi=q/d
               MMPOTEXP(MPOT,IM,th)=MMPOTEXP(MPOT,IM,th)+Q*D_1
!              field at charges E(r)=q/r^3*r
               factor1=Q*D_3
               MMPOTEXP(MEX,IM,th)=MMPOTEXP(MEX,IM,th)+factor1*DX
               MMPOTEXP(MEY,IM,th)=MMPOTEXP(MEY,IM,th)+factor1*DY
               MMPOTEXP(MEZ,IM,th)=MMPOTEXP(MEZ,IM,th)+factor1*DZ
   
              ELSE
!               use gaussian charges                 
                 sigmaT=epot1_fix%neardata(SIGCHARGE,IM)  
                 neg_inv_2_sigmaT_sqr = -1.D0 /(sigmaT*sigmaT)
                 inv_sigmaT = 1.D0 / sigmaT
                 prefac_sigmaT = 2.D0*M_SQRTPI_1                &
                                   / sigmaT
                 FktGauss  = DEXP(neg_inv_2_sigmaT_sqr*D*D)
                 FktErf_D    = DERF(inv_sigmaT*D)*D_1

!                potential             Phi = q * erf(d/sigmaT) / d
                 MMPOTEXP(MPOT,IM,th)=MMPOTEXP(MPOT,IM,th)      &
                                       +Q*FktErf_D
!                field     E(r)= q*r/r^3*erf(r/(sqrt(r)*sigmaT) - sqrt(2/pi)*q*r/r^2*1/sigmaT * exp(-r^2/(sigmaT^2)
                 factor1 = Q*D_2*                               &
                       (FktErf_D-prefac_sigmaT*FktGauss)
                 MMPOTEXP(MEX,IM,th)=MMPOTEXP(MEX,IM,th)+factor1*DX
                 MMPOTEXP(MEY,IM,th)=MMPOTEXP(MEY,IM,th)+factor1*DY
                 MMPOTEXP(MEZ,IM,th)=MMPOTEXP(MEZ,IM,th)+factor1*DZ
!             end of gaussian or point MM
              ENDIF
!           end of if CALC_MM
            ENDIF
! ++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++ PMM REGION of CALC_EXPORT++++++++++++++++++
!           if runinfo%pff and we have a MM atom carrying an induced dipole...
            IF (runinfo%pff.AND.epot1_fix%neardata(SIGDIPOLE,IM).GT.0.D0)   THEN
!             decide whether we are in the gauss- or point-region
              IF (D.GT.(runinfo%olensig*epot1_fix%neardata(SIGDIPOLE,IM))) THEN
!               use point dipoles
                factor1=Q*D_3
                MMPOTEXP_VAR(MEXPFF,IM,th)=MMPOTEXP_VAR(MEXPFF,IM,th)  &
                                          +factor1*DX
                MMPOTEXP_VAR(MEYPFF,IM,th)=MMPOTEXP_VAR(MEYPFF,IM,th)  &
                                          +factor1*DY
                MMPOTEXP_VAR(MEZPFF,IM,th)=MMPOTEXP_VAR(MEZPFF,IM,th)  &
                                          +factor1*DZ

!               calculate the field gradient if needed       
                IF (.NOT.epot1_var%tcalconlypfffield) THEN  
                  factor2 = 3.0 * factor1 * D_2
!                 fieldgradient        gradE(r)=q/r^3*(delta_ij-3*r_i*r_j/r^2)
                  MMPOTEXP(MEXX,IM,th)= MMPOTEXP(MEXX,IM,th)          &
                                            +(factor1-DX2*factor2)   
                  MMPOTEXP(MEYY,IM,th)= MMPOTEXP(MEYY,IM,th)          &
                                            +(factor1-DY2*factor2)   
                  MMPOTEXP(MEZZ,IM,th)= MMPOTEXP(MEZZ,IM,th)          &
                                            +(factor1-DZ2*factor2)   
                  MMPOTEXP(MEXY,IM,th)= MMPOTEXP(MEXY,IM,th)          &
                                         -DX*DY*factor2             
                  MMPOTEXP(MEXZ,IM,th)= MMPOTEXP(MEXZ,IM,th)          &
                                         -DX*DZ*factor2             
                  MMPOTEXP(MEYZ,IM,th)= MMPOTEXP(MEYZ,IM,th)          &
                                         -DY*DZ*factor2
                ENDIF

            ELSE
!             use gaussian dipoles
              sigmaT=epot1_fix%neardata(SIGDIPOLE,IM) 

              SIGMA_1 = 1.0D0/sigmaT !sqrt2 implied
              SIGMA_2 = SIGMA_1 * SIGMA_1 !sqrt2 implied
              ARG=D*SIGMA_1  !sqrt2 implied
              LAMBDA=(2.0D0 * M_SQRTPI_1 )                       &
                      * SIGMA_1 * DEXP(-ARG*ARG) !sqrt2 implied
              MU=D_1 * DERF(ARG)
              
              FAC1=(MU-LAMBDA)*D_2
              FAC2=Q*(3.0D0 * FAC1 - 2.0D0 *LAMBDA*SIGMA_2) * D_2  !factor 2 due to sqrt2 implied
              FAC1=Q*FAC1
              
!             the field at the PMM atoms
              MMPOTEXP_VAR(MEXPFF,IM,th)= MMPOTEXP_VAR(MEXPFF,IM,th) &
                                        +FAC1*DX
              MMPOTEXP_VAR(MEYPFF,IM,th)= MMPOTEXP_VAR(MEYPFF,IM,th) &
                                        +FAC1*DY
              MMPOTEXP_VAR(MEZPFF,IM,th)= MMPOTEXP_VAR(MEZPFF,IM,th) &
                                        +FAC1*DZ

!           calculate the field gradient at the PMM atoms if needed
              IF (.NOT.epot1_var%tcalconlypfffield) THEN
                MMPOTEXP(MEXX,IM,th)= MMPOTEXP(MEXX,IM,th)          &
                                    -FAC2*DX2+FAC1
                MMPOTEXP(MEYY,IM,th)= MMPOTEXP(MEYY,IM,th)          &
                                    -FAC2*DY2+FAC1
                MMPOTEXP(MEZZ,IM,th)= MMPOTEXP(MEZZ,IM,th)          &
                                    -FAC2*DZ2+FAC1
                MMPOTEXP(MEXY,IM,th)= MMPOTEXP(MEXY,IM,th)          &
                                    -FAC2*DX*DY 
                MMPOTEXP(MEXZ,IM,th)= MMPOTEXP(MEXZ,IM,th)          &
                                    -FAC2*DX*DZ                 
                MMPOTEXP(MEYZ,IM,th)= MMPOTEXP(MEYZ,IM,th)          &
                                    -FAC2*DY*DZ 
              ENDIF 

!............end of gaussian or point runinfo%pff
             ENDIF
!...........end of if runinfo%pff.AND.SIGDIP > 0
            ENDIF
            
!..........end of loop over all near mmatoms of voxel VO
           ENDDO
          CONTINUE

!.....end of loop over all gridpoints:
          ENDDO
         ENDDO
       ENDDO
!$omp end do 
!$omp end parallel 



!     compression of threaded variables to the unthreaded structures - part 1
!     voxel multipole moments: they have to be compressed before OLE_EXPORT!
      IF (parai%ncpus.GT.1) THEN !needed for intel
       DO th=2,parai%ncpus !parai%ncpus is the number of omp threads per mpi process (NPROC is number of mpi processes)
!       for all qm atoms sum up all moments
        DO VO=1,MYNVOX
          DO mom=LAD,dimmultexp
           VOMULTEXP(mom,VO,1)=VOMULTEXP(mom,VO,1)+VOMULTEXP(mom,VO,th)
          ENDDO
          DO mom=RGZERO,dimrgyrexp
            VOXRGYREXP(mom,VO,1)=VOXRGYREXP(mom,VO,1) &
             +VOXRGYREXP(mom,VO,th)
          ENDDO
        ENDDO
       ENDDO
      ENDIF

!     calculate ole contributions from VOMULTEXP to MMPOTEXP and shift VOMULTEXP to ATMULTEXP
      CALL OLE_EXPORT
      

      IF (parai%ncpus.GT.1) THEN !needed for intel
!     compression of threaded variables to the unthreaded structures  - part 2
       DO th=2,parai%ncpus !parai%ncpus is the number of omp threads per mpi process (NPROC is number of mpi processes)
!       for all qm atoms sum up all moments
        DO IQ=1,epot1_list%nloc
          DO mom=LAD,dimmultexp
           ATMULTEXP(mom,IQ,1)=ATMULTEXP(mom,IQ,1)+ATMULTEXP(mom,IQ,th)
          ENDDO
          DO mom=RGZERO,dimrgyrexp
           ATRGYREXP(mom,IQ,1)=ATRGYREXP(mom,IQ,1)+ATRGYREXP(mom,IQ,th)
          ENDDO
        ENDDO !IQ

!       for all near atoms sum up the potential expansions
        DO A=1,epot1_list%nnear
          DO mom=MPOT,dimpotexp
              MMPOTEXP(mom,A,1)  =MMPOTEXP(mom,A,1)  +MMPOTEXP(mom,A,th)
          ENDDO
          DO mom=MEXPFF,DIMPOTEXP_VAR
              MMPOTEXP_VAR(mom,A,1) = MMPOTEXP_VAR(mom,A,1) &
                                    +MMPOTEXP_VAR(mom,A,th)
          ENDDO
        ENDDO !A
       ENDDO !th
      END IF


!.....now add contributions from core charges to ATMULTEXP and MMPOTEXP......
      IF(paral%parent) CALL MMCOREFORCE(TAU0)
!..................................................

      IF (paral%io_parent)  WRITE (6,'(3A)') 'INTERFACE| ...ELECTROSTATICS AT MM-ATOMS CALCULATED'      
      CALL TIHALT(procedureN,ISUB)

!==--------------------------------------------------------------==
 RETURN
 END SUBROUTINE
!==================================================================



!==================================================================
SUBROUTINE OLE_EXPORT
!==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER                    :: procedureN = 'OLE_EXPORT'
      
      INTEGER ISUB,th
      INTEGER VO,O,OL,IQ,J
      REAL*8 DX,DY,DZ
      REAL*8 TAR(dimmultexp)

      REAL*8 OLEEXP(1+3+6),OLEEXPO(1+3),OLEEXPHMPOT

      REAL*8 f,f1,f2,f3,r5,D
      REAL*8 PdotR, PhiQ, QxLine, QyLine, QzLine,PR7_5_2PhiQ,PhiO
      REAL*8 myXXXX,myYYYY,myZZZZ
      REAL*8 ThreeR2X,ThreeR2Y,ThreeR2Z
      REAL*8 dR2eX2,dR2eY2,dR2eZ2,dR2eR2
      REAL*8 distMDXY,distMDXZ,distMDYZ,distMDXYZ
      REAL*8 rpowPR1,rpowPR2,rpowPR4,rpowPR3,rpowPR5,rpowPR7
      REAL*8 rpowPR2_7,rpowPR7_5,rpowPR7_10,rpowPR7_5_2

      LOGICAL ATOMHASDIPOLE

!$    INTEGER   OMP_GET_THREAD_NUM
!$    EXTERNAL  OMP_GET_THREAD_NUM


      CALL TISET(procedureN,ISUB)
      CALL zeroing(ATMULTEXP)!,dimmultexp*maxind*mxno_threads)

      if (paral%parent.AND.runinfo%tverbose) print *,"ENTER OLE_EXPORT"

!     important for non-openmp runs
      th=1

!$omp parallel private(VO,IQ,DX,DY,DZ,J,TAR,O,OL,OLEEXP,                 &
!$omp  OLEEXPO,dR2eX2,dR2eY2,dR2eZ2,dR2eR2,                              &
!$omp  D,distMDXY,distMDXZ,distMDYZ,distMDXYZ,ATOMHASDIPOLE,             &
!$omp  rpowPR1,rpowPR2,rpowPR4,rpowPR3,rpowPR5,rpowPR7,rpowPR2_7,        &
!$omp  rpowPR7_5,rpowPR7_10,rpowPR7_5_2,                                 &
!$omp  f,f1,PdotR,f2,f3,r5,PhiQ,PR7_5_2PhiQ,QxLine,QyLine,QzLine,PhiO,   &
!$omp  myXXXX,myYYYY,myZZZZ,ThreeR2X,ThreeR2Y,ThreeR2Z,OLEEXPHMPOT,th)   
!$     th=OMP_GET_THREAD_NUM()+1
!$omp do 
!.....handle multipoles and ole electrostatics for all my voxels
!.....note: the parai%sources are VOMULTEXP(*,VO,1), threads have been compressed before
      DO VO=1,MYNVOX
        IQ=VO2IQ(VO)

!.......shift multipole moments from voxels to qm atoms
!       vector from voxel to qm atom (=refpoint)
        DX=epot1_fix%koloc(1,IQ)-KOVOX(1,VO)
        DY=epot1_fix%koloc(2,IQ)-KOVOX(2,VO)
        DZ=epot1_fix%koloc(3,IQ)-KOVOX(3,VO)
        
!.......shift expansions from level C_-2 (voxels) to C_-1 (qm atoms)
        call SHIFT_MULTIPOLES_WRAP(TAR,DX,DY,DZ,VO) !reads VOMULTEXP,write TAR
        DO J=1,dimmultexp
          ATMULTEXP(J,IQ,th) = ATMULTEXP(J,IQ,th) + TAR(J) 
        ENDDO

!       add RGYR tensors
        DO J=1,dimrgyrexp
          ATRGYREXP(J,IQ,th) = ATRGYREXP(J,IQ,th) + VOXRGYREXP(J,VO,1) 
        ENDDO

!.......calculate potential, field, and field gradient at ole (C^-2) MM atoms from voxels
        DO O=1,OLELISTLEN(VO)
          OL=OLELIST(VO,O)
          ATOMHASDIPOLE=(epot1_fix%neardata(SIGDIPOLE,OL).GT.0.0D0)

!         vector from voxel to ole atom
          DX=epot1_fix%neardata(N_KX,OL)-KOVOX(1,VO)
          DY=epot1_fix%neardata(N_KY,OL)-KOVOX(2,VO)
          DZ=epot1_fix%neardata(N_KZ,OL)-KOVOX(3,VO)

!     no reset of OLEEXP,OLEEXP,OLEEXPHMPOT as they are overwritten 

!.....helping quantities
!     DR2 (vectorElementMultiplication)
      dR2eX2 = DX*DX           
      dR2eY2 = DY*DY           
      dR2eZ2 = DZ*DZ          
      dR2eR2 = dR2eX2+dR2eY2+dR2eZ2 
      D   = DSQRT(dR2eR2)
!     distm (calcDistanceMatrix)
      distMDXY = DX * DY; 
      distMDXZ = DX * DZ; 
      distMDYZ = DY * DZ; 
      distMDXYZ= distMDXY* DZ; 
!     distance matrix   (calcPowersOfDistance)             (scalefac = 1)
      rpowPR1      =   1.0D0     / D;     
      rpowPR2      = rpowPR1 * rpowPR1
      rpowPR4      = rpowPR2 * rpowPR2
      rpowPR3      = rpowPR2 * rpowPR1
      rpowPR5      = rpowPR3 * rpowPR2
      rpowPR7      = rpowPR5 * rpowPR2
      rpowPR2_7    =   7.0D0  * rpowPR2
      rpowPR7_5    =   5.0D0  * rpowPR7
      rpowPR7_10   =  10.0D0  * rpowPR7
      rpowPR7_5_2  =   2.5D0  * rpowPR7

! CHARGE      (setMonopoleContributionsToExpansion)
      f = VOMULTEXP(LAD,VO,1) * rpowPR1 
      OLEEXP(MPOT) = f

      f = f * rpowPR2
      OLEEXP(MEX) = f * dX
      OLEEXP(MEY) = f * dY
      OLEEXP(MEZ) = f * dZ
    
      IF (ATOMHASDIPOLE.AND..NOT.epot1_var%tcalconlypfffield) then
        f1 = -3.0D0 * f * rpowPR2
        OLEEXP(MEXX) = (f + dR2eX2 * f1)
        OLEEXP(MEYY) = (f + dR2eY2 * f1)
        OLEEXP(MEXY) =  distMDXY * f1
        OLEEXP(MEXZ) =  distMDXZ * f1
        OLEEXP(MEYZ) =  distMDYZ * f1
      ENDIF

! dipole    (addDipoleContributionsToExpansion)
      PdotR = VOMULTEXP(PX,VO,1)*DX+VOMULTEXP(PY,VO,1)*DY+VOMULTEXP(PZ,VO,1)*DZ

      OLEEXP(MPOT) = OLEEXP(MPOT) + PdotR * rpowPR3

      f2 = 3.0D0 * PdotR
      f3 = f2 * rpowPR2
      OLEEXP(MEX)= OLEEXP(MEX) + (f3 * dX -VOMULTEXP(PX,VO,1)) *rpowPR3
      OLEEXP(MEY)= OLEEXP(MEY) + (f3 * dY -VOMULTEXP(PY,VO,1)) *rpowPR3
      OLEEXP(MEZ)= OLEEXP(MEZ) + (f3 * dZ -VOMULTEXP(PZ,VO,1)) *rpowPR3

      IF (ATOMHASDIPOLE.AND..NOT.epot1_var%tcalconlypfffield) then   
        r5 = rpowPR5
        f3 = f3 * 5.0D0                                                                                                                                                                               

        OLEEXP(MEXX) = OLEEXP(MEXX)+(f2 + 6.0D0 * VOMULTEXP(PX,VO,1)*dX - dR2eX2*f3)*r5
        OLEEXP(MEYY) = OLEEXP(MEYY)+(f2 + 6.0D0 * VOMULTEXP(PY,VO,1)*dY - dR2eY2*f3)*r5
        OLEEXP(MEXY) = OLEEXP(MEXY)+(3.0D0 * (VOMULTEXP(PX,VO,1)*dY + VOMULTEXP(PY,VO,1)*dX) - distMDXY * f3) * r5
        OLEEXP(MEXZ) = OLEEXP(MEXZ)+(3.0D0 * (VOMULTEXP(PX,VO,1)*dZ + VOMULTEXP(PZ,VO,1)*dX) - distMDXZ * f3) * r5
        OLEEXP(MEYZ) = OLEEXP(MEYZ)+(3.0D0 * (VOMULTEXP(PY,VO,1)*dZ + VOMULTEXP(PZ,VO,1)*dY) - distMDYZ * f3) * r5
      ENDIF

! quadrupole (addQuadrupoleContributionsToExpansion)
      PhiQ =       VOMULTEXP(QXX,VO,1) * dR2eX2 +                      &
                   VOMULTEXP(QYY,VO,1) * dR2eY2 +                      &
                 (-VOMULTEXP(QXX,VO,1)-VOMULTEXP(QYY,VO,1)) * dR2eZ2 + & 
                   2.0D0 * VOMULTEXP(QXY,VO,1) * distMDXY +            &
                   2.0D0 * VOMULTEXP(QXZ,VO,1) * distMDXZ +            &
                   2.0D0 * VOMULTEXP(QYZ,VO,1) * distMDYZ

      OLEEXP(MPOT) = OLEEXP(MPOT) + 0.5 * PhiQ * rpowPR5

      PR7_5_2PhiQ = rpowPR7_5_2 * PhiQ
      QxLine = VOMULTEXP(QXX,VO,1) * dX + VOMULTEXP(QXY,VO,1)    * dY  &
                    + VOMULTEXP(QXZ,VO,1) * dZ
      QyLine = VOMULTEXP(QXY,VO,1) * dX + VOMULTEXP(QYY,VO,1)    * dY  &
                    + VOMULTEXP(QYZ,VO,1) * dZ
      QzLine = VOMULTEXP(QXZ,VO,1) * dX + VOMULTEXP(QYZ,VO,1)    * dY  &
                    + (-VOMULTEXP(QXX,VO,1)-VOMULTEXP(QYY,VO,1)) * dZ

      OLEEXP(MEX) = OLEEXP(MEX)+PR7_5_2PhiQ * dX - QxLine * rpowPR5
      OLEEXP(MEY) = OLEEXP(MEY)+PR7_5_2PhiQ * dY - QyLine * rpowPR5
      OLEEXP(MEZ) = OLEEXP(MEZ)+PR7_5_2PhiQ * dZ - QzLine * rpowPR5
    
      IF (ATOMHASDIPOLE.AND..NOT.epot1_var%tcalconlypfffield) then

        OLEEXP(MEXX) = OLEEXP(MEXX) + PR7_5_2PhiQ *  &
                     (1.0D0 - rpowPR2_7 * dR2eX2)    &
                      + rpowPR7_10 * QxLine * dX     &
                      - VOMULTEXP(QXX,VO,1) * rpowPR5
        OLEEXP(MEYY) = OLEEXP(MEYY) + PR7_5_2PhiQ *  &
                     (1.0D0 - rpowPR2_7 * dR2eY2)    &
                      + rpowPR7_10 * QyLine * dY     &
                      - VOMULTEXP(QYY,VO,1) * rpowPR5

        f = PR7_5_2PhiQ * rpowPR2_7

        OLEEXP(MEXY) = OLEEXP(MEXY) - f * distMDXY + rpowPR7_5 * &
                       (QxLine * dY +  QyLine * dX)              &
                       - VOMULTEXP(QXY,VO,1) * rpowPR5
        OLEEXP(MEXZ) = OLEEXP(MEXZ) - f * distMDXZ + rpowPR7_5 * &
                       (QxLine * dZ +  QzLine * dX)              &
                       - VOMULTEXP(QXZ,VO,1) * rpowPR5
        OLEEXP(MEYZ) = OLEEXP(MEYZ) - f * distMDYZ + rpowPR7_5 * &
                       (QyLine * dZ +  QzLine * dY)              &
                       - VOMULTEXP(QYZ,VO,1) * rpowPR5
      ENDIF

! octupole     (addOctupoleContributionsToExpansion)
      PhiO = ( 6.0D0 * VOMULTEXP(OXYZ,VO,1) * distmDXYZ +                &
      3.0D0*((VOMULTEXP(OXXY,VO,1)*dY+ VOMULTEXP(OXXZ,VO,1) * dZ )      &
                                                             *dR2eX2    &
       + ( VOMULTEXP(OYYX,VO,1) * dX + VOMULTEXP(OYYZ,VO,1) * dZ )      &
                                                             *dR2eY2    &
       + ( VOMULTEXP(OZZX,VO,1) * dX + VOMULTEXP(OZZY,VO,1) * dY )      &
                                                             *dR2eZ2)   &
       + (-VOMULTEXP(OYYX,VO,1) -VOMULTEXP(OZZX,VO,1)) * dR2eX2 * dX    &
       + (-VOMULTEXP(OXXY,VO,1) -VOMULTEXP(OZZY,VO,1)) * dR2eY2 * dY    &
       + (-VOMULTEXP(OYYZ,VO,1) -VOMULTEXP(OXXZ,VO,1)) * dR2eZ2 * dZ )  &
       *rpowPR2

       OLEEXPO(MPOT) = (0.5D0*PhiO*rpowPR5)

       OLEEXPO(MEX) =  ( 3.5D0 * dX * PhiO                             &
      -1.5D0 * (VOMULTEXP(OYYX,VO,1) * ( dR2eY2 - dR2eX2 ) +           &
                VOMULTEXP(OZZX,VO,1) * ( dR2eZ2 - dR2eX2 ) +           &
       2.0D0 * (VOMULTEXP(OXXY,VO,1) * distMDXY +                      &
                VOMULTEXP(OXXZ,VO,1) * distMDXZ +                      &
                VOMULTEXP(OXYZ,VO,1) * distMDYZ   ) ) ) * rpowPR7
       OLEEXPO(MEY) =  ( 3.5D0 * dY * PhiO                             &
      -1.5D0 * (VOMULTEXP(OXXY,VO,1) * ( dR2eX2 - dR2eY2 ) +           &
                VOMULTEXP(OZZY,VO,1) * ( dR2eZ2 - dR2eY2 ) +           &
       2.0D0 * (VOMULTEXP(OYYX,VO,1) * distMDXY +                      &
                VOMULTEXP(OXYZ,VO,1) * distMDXZ +                      &
                VOMULTEXP(OYYZ,VO,1) * distMDYZ ) ) ) * rpowPR7
       OLEEXPO(MEZ) =  ( 3.5D0 * dZ * PhiO                             &
      -1.5D0 * (VOMULTEXP(OXXZ,VO,1) * ( dR2eX2 - dR2eZ2 ) +           &
                VOMULTEXP(OYYZ,VO,1) * ( dR2eY2 - dR2eZ2 ) +           &
       2.0D0 * (VOMULTEXP(OXYZ,VO,1) * distMDXY +                      &
                VOMULTEXP(OZZX,VO,1) * distMDXZ +                      &
                VOMULTEXP(OZZY,VO,1) * distMDYZ ) ) ) * rpowPR7

      IF (ATOMHASDIPOLE.AND..NOT.epot1_var%tcalconlypfffield) then
          OLEEXP(MEXX) = OLEEXP(MEXX)                                   &
            + 7.0D0*OLEEXPO(MPOT)*rpowPR2*(1.0D0+5.0D0*rpowPR2*dR2eX2)  & 
            - 7.0D0*rpowPR2*(            2.0D0*OLEEXPO(MEX) * dX)       &
            - 1.0D0/6.0D0 *rpowPR5*rpowPR2 *                            &
            18.0D0*(dX*(- VOMULTEXP(OYYX,VO,1) - VOMULTEXP(OZZX,VO,1))  &
            + dY* VOMULTEXP(OXXY,VO,1)                                  &
            + dZ* VOMULTEXP(OXXZ,VO,1))                                 
          OLEEXP(MEYY) = OLEEXP(MEYY)                                   &
            + 7.0D0*OLEEXPO(MPOT)*rpowPR2*(1.0D0+5.0D0*rpowPR2* dR2eY2) &  
            - 7.0D0*rpowPR2*(            2.0D0*OLEEXPO(MEY) * dY)       &
          - 1.0D0/6.0D0*rpowPR5*rpowPR2*18.0D0*(dX*VOMULTEXP(OYYX,VO,1) &        
            + dY*(- VOMULTEXP(OXXY,VO,1) - VOMULTEXP(OZZY,VO,1))        &
            + dZ* VOMULTEXP(OYYZ,VO,1))                                 
          OLEEXP(MEXY) = OLEEXP(MEXY)                                   &
            + 7.0D0*OLEEXPO(MPOT)*rpowPR2*(    5.0D0*rpowPR2* distMDXY) &   
            - 7.0D0*rpowPR2*(OLEEXPO(MEX) * dY + OLEEXPO(MEY) * dX)     &
          - 1.0D0/6.0D0*rpowPR5*rpowPR2*18.0D0*(dZ*VOMULTEXP(OXYZ,VO,1) &        
            + dY* VOMULTEXP(OYYX,VO,1)                                  &
            + dX* VOMULTEXP(OXXY,VO,1))                                 
          OLEEXP(MEXZ) = OLEEXP(MEXZ)                                   &
            +7.0D0*OLEEXPO(MPOT)*rpowPR2*(     5.0D0*rpowPR2* distMDXZ) &   
            - 7.0D0*rpowPR2*(OLEEXPO(MEX) * dZ + OLEEXPO(MEZ) * dX)     &
          - 1.0D0/6.0D0*rpowPR5*rpowPR2*18.0D0*(dY*VOMULTEXP(OXYZ,VO,1) &        
            + dZ* VOMULTEXP(OZZX,VO,1)                                  &
            + dX* VOMULTEXP(OXXZ,VO,1))
          OLEEXP(MEYZ) = OLEEXP(MEYZ)                                   &
            +7.0D0*OLEEXPO(MPOT)*rpowPR2*(     5.0D0*rpowPR2* distMDYZ) &   
            - 7.0D0*rpowPR2*(OLEEXPO(MEY) * dZ + OLEEXPO(MEZ) * dY)     &
          - 1.0D0/6.0D0*rpowPR5*rpowPR2*18.0D0*(dX*VOMULTEXP(OXYZ,VO,1) &        
            + dZ* VOMULTEXP(OZZY,VO,1)                                  &
            + dY* VOMULTEXP(OYYZ,VO,1))
      ENDIF

! hexadecapole (addHexadecapoleContributionsToExpansion)
      myXXXX = dR2eX2 * dR2eX2
      myYYYY = dR2eY2 * dR2eY2
      myZZZZ = dR2eZ2 * dR2eZ2
      ThreeR2X = 3.0 * dR2eX2
      ThreeR2Y = 3.0 * dR2eY2
      ThreeR2Z = 3.0 * dR2eZ2
        
      OLEEXPHMPOT = 0.125D0 * (                                          &
         (6.0D0* dR2eX2 * dR2eY2 -myXXXX -myYYYY)*VOMULTEXP(HXXYY,VO,1)  &
       + (6.0D0* dR2eX2 * dR2eZ2 -myXXXX -myZZZZ)*VOMULTEXP(HXXZZ,VO,1)  &
       + (6.0D0* dR2eY2 * dR2eZ2 -myYYYY -myZZZZ)*VOMULTEXP(HYYZZ,VO,1)  &
       + 4.0D0 * (distMDXY *((dR2eX2 - ThreeR2Z)*VOMULTEXP(HXXXY,VO,1)   &
       + ( dR2eY2 - ThreeR2Z ) * VOMULTEXP(HYYYX,VO,1) )                 &
       + distMDXZ * ( ( dR2eX2 - ThreeR2Y ) * VOMULTEXP(HXXXZ,VO,1)      &
       + ( dR2eZ2 - ThreeR2Y ) * VOMULTEXP(HZZZX,VO,1) )                 &
       + distMDYZ * ( ( dR2eY2 - ThreeR2X ) * VOMULTEXP(HYYYZ,VO,1)      &
       + ( dR2eZ2 - ThreeR2X ) * VOMULTEXP(HZZZY,VO,1) ))   )            &
        * rpowPR5 * rpowPR4

        OLEEXP(MEX) = OLEEXP(MEX)+ dX * 9.0D0*OLEEXPHMPOT*rpowPR2       &
               - 0.5D0 * rpowPR5*rpowPR4 * (                            &
               + (3.0D0* dR2eY2 -   dR2eX2)* dX * VOMULTEXP(HXXYY,VO,1) &
               + (3.0D0* dR2eZ2 -   dR2eX2)* dX * VOMULTEXP(HXXZZ,VO,1) &
               +  3.0D0*(dR2eX2 -   dR2eZ2)* dY * VOMULTEXP(HXXXY,VO,1) &
               +    (dR2eY2 - 3.0D0*dR2eZ2)* dY * VOMULTEXP(HYYYX,VO,1) &
               +  3.0D0*(dR2eX2 -   dR2eY2)* dZ * VOMULTEXP(HXXXZ,VO,1)  &
               +   (dR2eZ2 - 3.0D0*dR2eY2)* dZ * VOMULTEXP(HZZZX,VO,1)   &
       -6.0D0*distmDXYZ *(VOMULTEXP(HYYYZ,VO,1)+VOMULTEXP(HZZZY,VO,1)))  
        OLEEXP(MEY) = OLEEXP(MEY) + dY * 9.0D0*OLEEXPHMPOT*rpowPR2       &
               - 0.5D0 * rpowPR5*rpowPR4 * (                             &
               + (3.0D0* dR2eX2 -  dR2eY2)* dY * VOMULTEXP(HXXYY,VO,1)   &
               + (3.0D0* dR2eZ2 -  dR2eY2)* dY * VOMULTEXP(HYYZZ,VO,1)   &
               +   (dR2eX2 - 3.0D0*dR2eZ2)* dX * VOMULTEXP(HXXXY,VO,1)   &
               +  3.0D0*(dR2eY2 -  dR2eZ2)* dX * VOMULTEXP(HYYYX,VO,1)   &
               +  3.0D0*(dR2eY2 -  dR2eX2)* dZ * VOMULTEXP(HYYYZ,VO,1)   &
               +   (dR2eZ2 - 3.0D0*dR2eX2)* dZ * VOMULTEXP(HZZZY,VO,1)   &
       - 6.0D0*distmDXYZ*(VOMULTEXP(HXXXZ,VO,1)+VOMULTEXP(HZZZX,VO,1)))
        OLEEXP(MEZ) = OLEEXP(MEZ)+ dZ * 9.0D0*OLEEXPHMPOT*rpowPR2       &
                - 0.5D0 * rpowPR5*rpowPR4 * (                           &
                + (3.0D0*dR2eX2 -   dR2eZ2)* dZ * VOMULTEXP(HXXZZ,VO,1) &
                + (3.0D0*dR2eY2 -   dR2eZ2)* dZ * VOMULTEXP(HYYZZ,VO,1) &
                +    (dR2eX2 -3.0D0*dR2eY2)* dX * VOMULTEXP(HXXXZ,VO,1) &
                + 3.0D0*(dR2eZ2 -   dR2eY2)* dX * VOMULTEXP(HZZZX,VO,1) &
                +    (dR2eY2 -3.0D0*dR2eX2)* dY * VOMULTEXP(HYYYZ,VO,1) &
                + 3.0D0*(dR2eZ2 -   dR2eX2)* dY * VOMULTEXP(HZZZY,VO,1) &
       -6.0D0*distmDXYZ *(VOMULTEXP(HXXXY,VO,1)+VOMULTEXP(HYYYX,VO,1)))


       IF (.NOT.epot1_var%tcalconlypfffield) THEN
         MMPOTEXP(MPOT,OL,th)=MMPOTEXP(MPOT,OL,th)+ OLEEXP(MPOT)+OLEEXPO(MPOT)+OLEEXPHMPOT

         MMPOTEXP(MEX,OL,th)=MMPOTEXP(MEX,OL,th)  + OLEEXP(MEX)+OLEEXPO(MEX)
         MMPOTEXP(MEY,OL,th)=MMPOTEXP(MEY,OL,th)  + OLEEXP(MEY)+OLEEXPO(MEY)
         MMPOTEXP(MEZ,OL,th)=MMPOTEXP(MEZ,OL,th)  + OLEEXP(MEZ)+OLEEXPO(MEZ)
       ENDIF

       IF (epot1_fix%neardata(SIGDIPOLE,OL).GT.0.0D0) THEN
         MMPOTEXP_VAR(MEXPFF,OL,th)=MMPOTEXP_VAR(MEXPFF,OL,th) +OLEEXP(MEX)+OLEEXPO(MEX)
         MMPOTEXP_VAR(MEYPFF,OL,th)=MMPOTEXP_VAR(MEYPFF,OL,th) +OLEEXP(MEY)+OLEEXPO(MEY)
         MMPOTEXP_VAR(MEZPFF,OL,th)=MMPOTEXP_VAR(MEZPFF,OL,th) +OLEEXP(MEZ)+OLEEXPO(MEZ)

         IF (.NOT.epot1_var%tcalconlypfffield) THEN
           MMPOTEXP(MEXX,OL,th)=MMPOTEXP(MEXX,OL,th)+OLEEXP(MEXX)
           MMPOTEXP(MEYY,OL,th)=MMPOTEXP(MEYY,OL,th)+OLEEXP(MEYY)
           MMPOTEXP(MEZZ,OL,th)=MMPOTEXP(MEZZ,OL,th)-OLEEXP(MEXX)-OLEEXP(MEYY)
           MMPOTEXP(MEXY,OL,th)=MMPOTEXP(MEXY,OL,th)+OLEEXP(MEXY)
           MMPOTEXP(MEXZ,OL,th)=MMPOTEXP(MEXZ,OL,th)+OLEEXP(MEXZ)
           MMPOTEXP(MEYZ,OL,th)=MMPOTEXP(MEYZ,OL,th)+OLEEXP(MEYZ)
         ENDIF
       ENDIF
!       end of loop OL over ole atoms of voxel VO
        ENDDO
!     end of loop VO over voxels
      ENDDO
!$omp end do
!$omp end parallel

      CALL TIHALT(procedureN,ISUB)
END SUBROUTINE
!==================================================================





!==================================================================
SUBROUTINE MMCOREFORCE(TAU0)
!==--------------------------------------------------------------==
      IMPLICIT NONE
!     Arguments
      REAL*8  TAU0(:,:,:)
!     Variables
      REAL*8  sigmaT,neg_inv_2_sigmaT_sqr,inv_sigmaT,prefac_sigmaT
      REAL*8  factor1,factor2,PSSIG,Q
      REAL*8  FktGauss,FktErf
      REAL*8  D,D_1,D_2,D_3
      REAL*8  D_X,D_Y,D_Z
      REAL*8  DX,DY,DZ,DX2,DY2,DZ2
      INTEGER IA,IS,IQ,A,IM,TMP
      REAL*8  DERF

      REAL*8  SIGMA_1,SIGMA_2,ARG,LAMBDA,MU,FAC1,FAC2

      IF(.NOT.paral%parent) RETURN
      
      TMP = 0
!.....loop over NSP atomspecies IS
      DO IS=1,ions1%nsp
!.......store charge and gaussian width of current species
        Q=ions0%zv(is)
        PSSIG=RAGGIO(IS)
!       RAGGIO gives the radius of the core charge,
!       per default, RAGGIO(IS) is 1.2 for all atom species IS,
!       the data is stored in the structure DEFRAG2 in atoms.F .
!       A factor of sqrt(2) is already implicated in the value of
!       RAGGIO, cf. the use of RAGGIO in the calculation of
!       energy and force on the ions due to the external potential in
!       the function VEXTNU in eextern.F. There, the sqrt(2) is omitted
!       in the gaussian function, too.

!.......for each species IS loop over NA qm-atoms 
        DO IA=1,ions0%na(IS)
          TMP = TMP + 1
          IF (epot1_list%nloc.EQ.0) THEN
             CALL STOPGM('MMCOREF','epot1_list%nloc is zero! This cannot be!',&
                   __LINE__,__FILE__)
          ENDIF
!.........loop over all qm atoms
          DO A=1,epot1_list%nloc
            D_X=epot1_fix%koloc(1,A)-TAU0(1,IA,IS)
            D_Y=epot1_fix%koloc(2,A)-TAU0(2,IA,IS)
            D_Z=epot1_fix%koloc(3,A)-TAU0(3,IA,IS)

!           calculate distance between qm-atom and atom in species-list.
!           if this is small, identify them as the same
            IF(D_X*D_X+D_Y*D_Y+D_Z*D_Z.LT.0.1D0) THEN
               IQ=A
               GOTO 400
            ENDIF
          ENDDO
          CALL STOPGM('MMCOREF','NEAREST EXPANSION NOT FOUND',&
                   __LINE__,__FILE__)
 400      CONTINUE

!       print *,"QMATOMS",IA,IS,TMP,IQ
        if (TMP.NE.IQ) print *,"STRANGE DAYS"
          
!.......write out core charges WRITECHARGE
        IF (0.EQ.1) THEN
        WRITE(*,'(F21.14,TR1,F21.14,TR1,F21.14,TR1,F21.14,TR1,F21.14,A)') &
             epot1_fix%koloc(1,IQ),epot1_fix%koloc(2,IQ),epot1_fix%koloc(3,IQ),Q,PSSIG,' QC'
        ENDIF

!++++++++ THE SAMM PART   OF  MM  C O R E  FORCE
!.........add core charge of current qm-atom to its monopole moment (no voxels for core-mm interactions)
          ATMULTEXP(LAD,IQ,1)=ATMULTEXP(LAD,IQ,1)+Q

!++++++++ THE NEAR PART   OF  MM  C O R E  FORCE

!.........loop over all mm-atoms A which are near qm-atom IQ  
          DO A=1,epot1_list%nincl(IQ)
! ...........atomid IM of nearatom A from qmatom IQ
             IM=epot1_list%lincl(IQ,A)           
!............vector from qm-core to near mm-atom
             DX=epot1_fix%neardata(N_KX,IM)-epot1_fix%koloc(1,IQ)
             DY=epot1_fix%neardata(N_KY,IM)-epot1_fix%koloc(2,IQ)
             DZ=epot1_fix%neardata(N_KZ,IM)-epot1_fix%koloc(3,IQ)

             DX2 = DX*DX
             DY2 = DY*DY
             DZ2 = DZ*DZ

             D = DSQRT(DX2+DY2+DZ2)
             D_1 = 1.0D0/D
             D_2 = D_1 * D_1
             D_3 = D_2 * D_1

! ++++++++++ MM REGION of MMCOREFORCE
!           potential and field at charges (not at beginning of integration step)
            IF (.NOT.runinfo%pff.OR.(.NOT.epot1_var%tcalconlypfffield)) THEN
!             gaussian or point MM
              sigmaT=DSQRT(epot1_fix%neardata(SIGCHARGE,IM)**2+PSSIG**2) !two gauss objects: use mixed gaussian width
!             factor sqrt(2) is already in epot1_fix%neardata(SIGCHARGE,IM) and also in PSSIG

              IF (D.GT.(8.60D0/DSQRT(2.0D0)*sigmaT)) THEN
!               use point charges
                MMPOTEXP(MPOT,IM,1)=MMPOTEXP(MPOT,IM,1)+Q*D_1

                factor1=Q*D_3
                MMPOTEXP(MEX,IM,1)=MMPOTEXP(MEX,IM,1)+factor1*DX
                MMPOTEXP(MEY,IM,1)=MMPOTEXP(MEY,IM,1)+factor1*DY
                MMPOTEXP(MEZ,IM,1)=MMPOTEXP(MEZ,IM,1)+factor1*DZ
              ELSE
!               use gaussian charges
                neg_inv_2_sigmaT_sqr = -1.0D0/(sigmaT*sigmaT)
                inv_sigmaT = 1.0D0 / sigmaT
                prefac_sigmaT = 2.0D0*M_SQRTPI_1/ sigmaT
                FktGauss  = DEXP(neg_inv_2_sigmaT_sqr*D*D)
                FktErf    = DERF(inv_sigmaT*D)

                MMPOTEXP(MPOT,IM,1)=MMPOTEXP(MPOT,IM,1)+Q*FktErf / D

                factor1 = Q*(D_3*FktErf-prefac_sigmaT*D_2*FktGauss)
                MMPOTEXP(MEX,IM,1)=MMPOTEXP(MEX,IM,1)+factor1*DX
                MMPOTEXP(MEY,IM,1)=MMPOTEXP(MEY,IM,1)+factor1*DY
                MMPOTEXP(MEZ,IM,1)=MMPOTEXP(MEZ,IM,1)+factor1*DZ
!             end of gaussian or point MM
              ENDIF
!           end of if not runinfo%pff or runinfo%pff iteration has converged:
            ENDIF


! ++++++++++ PMM REGION of MMCOREFORCE
!           if runinfo%pff and we have a MM atom carrying an induced dipole...
            IF (runinfo%pff.AND.epot1_fix%neardata(SIGDIPOLE,IM).GT.0.D0)   THEN
              sigmaT=DSQRT(epot1_fix%neardata(SIGDIPOLE,IM)**2+PSSIG**2)  !two gauss objects: use mixed gaussian width
!             gaussian or point PMM
              IF (D.GT.(8.60D0/DSQRT(2.0D0)*sigmaT)) THEN
!               use point dipoles
!               for runinfo%pff calculate field and fieldgradient at polarizable atoms
                factor1=Q*D_3
                MMPOTEXP_VAR(MEXPFF,IM,1)=MMPOTEXP_VAR(MEXPFF,IM,1) +factor1*DX
                MMPOTEXP_VAR(MEYPFF,IM,1)=MMPOTEXP_VAR(MEYPFF,IM,1) +factor1*DY
                MMPOTEXP_VAR(MEZPFF,IM,1)=MMPOTEXP_VAR(MEZPFF,IM,1) +factor1*DZ

!               calculate the field gradient at the PMM atoms if needed     
                IF (.NOT.epot1_var%tcalconlypfffield) THEN
                 factor2 = 3.0 * factor1 * D_2
                 MMPOTEXP(MEXX,IM,1)= MMPOTEXP(MEXX,IM,1) +(factor1-DX*DX*factor2)
                 MMPOTEXP(MEYY,IM,1)= MMPOTEXP(MEYY,IM,1) +(factor1-DY*DY*factor2)
                 MMPOTEXP(MEZZ,IM,1)= MMPOTEXP(MEZZ,IM,1) +(factor1-DZ*DZ*factor2)
                 MMPOTEXP(MEXY,IM,1)= MMPOTEXP(MEXY,IM,1) -DX*DY*factor2
                 MMPOTEXP(MEXZ,IM,1)= MMPOTEXP(MEXZ,IM,1) -DX*DZ*factor2
                 MMPOTEXP(MEYZ,IM,1)= MMPOTEXP(MEYZ,IM,1) -DY*DZ*factor2
                ENDIF
              ELSE
!               use gaussian dipoles  
                SIGMA_1 = 1.0D0/sigmaT !sqrt2 implied
                SIGMA_2 = SIGMA_1 * SIGMA_1 !sqrt2 implied
                ARG=D*SIGMA_1  !sqrt2 implied
                LAMBDA=(2.0D0 * M_SQRTPI_1 ) * SIGMA_1 * DEXP(-ARG*ARG) !sqrt2 implied
                MU=D_1 * DERF(ARG)
                 
                FAC1=(MU-LAMBDA)*D_2
                FAC2=Q*(3.0D0 * FAC1 - 2.0D0 *LAMBDA*SIGMA_2) * D_2  !factor 2 due to sqrt2 implied

                FAC1=Q*FAC1
                 
!               save field at dipoles in VAR structures                 
                MMPOTEXP_VAR(MEXPFF,IM,1)= MMPOTEXP_VAR(MEXPFF,IM,1) +FAC1*DX
                MMPOTEXP_VAR(MEYPFF,IM,1)= MMPOTEXP_VAR(MEYPFF,IM,1) +FAC1*DY
                MMPOTEXP_VAR(MEZPFF,IM,1)= MMPOTEXP_VAR(MEZPFF,IM,1) +FAC1*DZ

!               calculate the field gradient at the PMM atoms if needed
                IF (.NOT.epot1_var%tcalconlypfffield) THEN
                  MMPOTEXP(MEXX,IM,1)= MMPOTEXP(MEXX,IM,1) -FAC2*DX*DX+FAC1
                  MMPOTEXP(MEYY,IM,1)= MMPOTEXP(MEYY,IM,1) -FAC2*DY*DY+FAC1
                  MMPOTEXP(MEZZ,IM,1)= MMPOTEXP(MEZZ,IM,1) -FAC2*DZ*DZ+FAC1
                  MMPOTEXP(MEXY,IM,1)= MMPOTEXP(MEXY,IM,1) -FAC2*DX*DY 
                  MMPOTEXP(MEXZ,IM,1)= MMPOTEXP(MEXZ,IM,1) -FAC2*DX*DZ 
                  MMPOTEXP(MEYZ,IM,1)= MMPOTEXP(MEYZ,IM,1) -FAC2*DY*DZ 
               ENDIF 
!............end of gaussian or point runinfo%pff
             ENDIF

!...........end of runinfo%pff
            ENDIF

!......end of loop over all mm-atoms A which are near qm-atom IQ
       ENDDO

!.......end of loop over NA qm-atoms of species IS
        ENDDO
!.....end of loop over NSP atomspecies IS
      ENDDO
!==--------------------------------------------------------------==
 RETURN
 END SUBROUTINE
!==================================================================



!==================================================================
SUBROUTINE QMCOREFORCE(TAU0,EEXT,FION)
!==--------------------------------------------------------------==
      IMPLICIT NONE
!     Arguments
      REAL*8  TAU0(:,:,:),EEXT,FION(:,:,:)
!     Variables
      REAL*8  TMP_POT,TMP_FIELD(3)
      REAL*8  D,D_1,D_2,D_3
      REAL*8  D_X,D_Y,D_Z
      REAL*8  DX,DY,DZ,DX2,DY2,DZ2
      INTEGER IA,IS,IQ,A,IM
      REAL*8  DERF

      REAL*8  sigmaT,neg_inv_2_sigmaT_sqr,inv_sigmaT,prefac_sigmaT
      REAL*8  factor1,PSSIG,Q
      REAL*8  FktGauss,FktErf
      REAL*8  NEARPOT_PFF
      REAL*8  SIGMA_1,SIGMA_2,ARG,LAMBDA,MU,FAC1,FAC2
      
      IF(.NOT.paral%parent) RETURN
      
!.....this is symetric to MMCOREFORCE

!.....loop over NSP atomspecies IS
      DO IS=1,ions1%nsp
!.......store charge and gaussian width of current species
        Q=ions0%zv(is)
        PSSIG=RAGGIO(IS)
!       RAGGIO gives the radius of the core charge, cf. comment im MMCOREFORCE

!.......for each species IS loop over NA qm-atoms 
        DO IA=1,ions0%na(IS)
          IF (epot1_list%nloc.EQ.0) GOTO 400
!.........loop over all qm atoms
          DO A=1,epot1_list%nloc
            D_X=epot1_fix%koloc(1,A)-TAU0(1,IA,IS)
            D_Y=epot1_fix%koloc(2,A)-TAU0(2,IA,IS)
            D_Z=epot1_fix%koloc(3,A)-TAU0(3,IA,IS)

!           calculate distance between qm-atom and atom in species-list.
!           if this is small, identify them as the same
            IF(D_X*D_X+D_Y*D_Y+D_Z*D_Z.LT.0.1D0) THEN
               IQ=A
               GOTO 400
            ENDIF
          ENDDO
          PRINT *,"ERROR! if you want to use QMHOMFIELD, ", &
                   "set #define NEWCOREFORCE 0 in eextern.F"
          CALL STOPGM('QMCOREF','NEAREST QM ATOM NOT FOUND',&
                   __LINE__,__FILE__)
 400      CONTINUE

!++++++++ THE SAMM PART OF  QM  C O R E  FORCE
!       evalute taylor expansion at reference point
!       epot1_var%locexp(K?,*) already has negative sign from iffi!
        TMP_POT=epot1_var%locexp(KPHI,IQ)
        TMP_FIELD(1)=epot1_var%locexp(KX,IQ)
        TMP_FIELD(2)=epot1_var%locexp(KY,IQ)
        TMP_FIELD(3)=epot1_var%locexp(KZ,IQ)

!++++++++ THE NEAR PART   OF  QM  C O R E  FORCE
!.........loop over all mm-atoms A which are near qm-atom IQ  
          DO A=1,epot1_list%nincl(IQ)
! ...........atomid IM of nearatom A from qmatom IQ
             IM=epot1_list%lincl(IQ,A)           
!............vector from near mm-atom to qm-core 
             DX=epot1_fix%koloc(1,IQ)-epot1_fix%neardata(N_KX,IM)
             DY=epot1_fix%koloc(2,IQ)-epot1_fix%neardata(N_KY,IM)
             DZ=epot1_fix%koloc(3,IQ)-epot1_fix%neardata(N_KZ,IM)

             DX2 = DX*DX
             DY2 = DY*DY
             DZ2 = DZ*DZ

             D = DSQRT(DX2+DY2+DZ2)
             D_1 = 1.0D0/D
             D_2 = D_1 * D_1
             D_3 = D_2 * D_1

! ++++++++++ MM REGION of QMCOREFORCE
           sigmaT=DSQRT(epot1_fix%neardata(SIGCHARGE,IM)**2+PSSIG**2) !two gauss objects: use mixed gaussian width
!          decide if we are far enough away to use point objects
           IF (D.GT.(8.60D0/DSQRT(2.0D0)*sigmaT)) THEN
             TMP_POT=TMP_POT + epot1_fix%neardata(PCHARGE,IM) * D_1
             factor1=epot1_fix%neardata(PCHARGE,IM) * D_3
             TMP_FIELD(1)=TMP_FIELD(1)+factor1*DX
             TMP_FIELD(2)=TMP_FIELD(2)+factor1*DY
             TMP_FIELD(3)=TMP_FIELD(3)+factor1*DZ
           ELSE
             neg_inv_2_sigmaT_sqr = -1.0D0/(sigmaT*sigmaT)
             inv_sigmaT = 1.0D0 / sigmaT
             prefac_sigmaT = 2.0D0*M_SQRTPI_1 / sigmaT
             FktGauss  = DEXP(neg_inv_2_sigmaT_sqr*D*D)
             FktErf    = DERF(inv_sigmaT*D)

             TMP_POT=TMP_POT + epot1_fix%neardata(PCHARGE,IM)*FktErf * D_1

             factor1 = epot1_fix%neardata(PCHARGE,IM)*           &
                      (D_3*FktErf-prefac_sigmaT*D_2*FktGauss)
             TMP_FIELD(1)=TMP_FIELD(1)+factor1*DX
             TMP_FIELD(2)=TMP_FIELD(2)+factor1*DY
             TMP_FIELD(3)=TMP_FIELD(3)+factor1*DZ
           ENDIF

! ++++++++++ PMM REGION of QMCOREFORCE
!           if runinfo%pff and we have a MM atom carrying an induced dipole...
            IF (runinfo%pff.AND.epot1_fix%neardata(SIGDIPOLE,IM).GT.0.D0)   THEN 
              sigmaT=DSQRT(epot1_fix%neardata(SIGDIPOLE,IM)**2+PSSIG**2)  !two gauss objects: use mixed gaussian width
!             decide if we are far enough away to use point objects
              IF (D.GT.(8.60D0/DSQRT(2.0D0)*sigmaT)) THEN
                 NEARPOT_PFF = D_3 *                                  &
                               (epot1_var%neardata_var(N_PX,IM)*DX +  &
                                epot1_var%neardata_var(N_PY,IM)*DY +  &
                                epot1_var%neardata_var(N_PZ,IM)*DZ)
 
                 TMP_POT=TMP_POT + NEARPOT_PFF

!               the field by a gaussian dipole                  
                 FAC1=D_3
                 FAC2=(3.0D0 * FAC1) * D_2 

                 TMP_FIELD(1) = TMP_FIELD(1) + (                          &
                    (FAC2*DX*DX - FAC1) * epot1_var%neardata_var(N_PX,IM) &
                   +(FAC2*DX*DY       ) * epot1_var%neardata_var(N_PY,IM) &
                   +(FAC2*DX*DZ       ) * epot1_var%neardata_var(N_PZ,IM))

                 TMP_FIELD(2) = TMP_FIELD(2) + (                          &
                    (FAC2*DY*DX       ) * epot1_var%neardata_var(N_PX,IM) &
                   +(FAC2*DY*DY - FAC1) * epot1_var%neardata_var(N_PY,IM) &
                   +(FAC2*DY*DZ       ) * epot1_var%neardata_var(N_PZ,IM))

                 TMP_FIELD(3) = TMP_FIELD(3) + (                          &
                    (FAC2*DZ*DX       ) * epot1_var%neardata_var(N_PX,IM) &
                   +(FAC2*DZ*DY       ) * epot1_var%neardata_var(N_PY,IM) &
                   +(FAC2*DZ*DZ - FAC1) * epot1_var%neardata_var(N_PZ,IM))
              ELSE
                 !formulas taken from appendix B of sebastian bauer's diss
                 SIGMA_1 = 1.0D0/sigmaT !sqrt2 implied
                 SIGMA_2 = SIGMA_1 * SIGMA_1
                 ARG=D*SIGMA_1  !sqrt2 implied
                 LAMBDA=(2.0D0 * M_SQRTPI_1 ) * SIGMA_1 * DEXP(-ARG*ARG) !sqrt2 implied
                 MU=D_1 * DERF(ARG)
                 NEARPOT_PFF = (MU-LAMBDA) * D_2 *                   &
                               (epot1_var%neardata_var(N_PX,IM)*DX + &
                                epot1_var%neardata_var(N_PY,IM)*DY + &
                                epot1_var%neardata_var(N_PZ,IM)*DZ)
 
                 TMP_POT=TMP_POT + NEARPOT_PFF

!               the field by a gaussian dipole                  
                 FAC1=(MU-LAMBDA) * D_2
                 FAC2=(3.0D0 * FAC1 - 2.0D0 *LAMBDA*SIGMA_2) * D_2  !factor 2 due to sqrt2 implied

                 TMP_FIELD(1) = TMP_FIELD(1) + (                           &
                    (FAC2*DX*DX - FAC1) * epot1_var%neardata_var(N_PX,IM)  &
                   +(FAC2*DX*DY       ) * epot1_var%neardata_var(N_PY,IM)  &
                   +(FAC2*DX*DZ       ) * epot1_var%neardata_var(N_PZ,IM))

                 TMP_FIELD(2) = TMP_FIELD(2) + (                           &
                    (FAC2*DY*DX       ) * epot1_var%neardata_var(N_PX,IM)  &
                   +(FAC2*DY*DY - FAC1) * epot1_var%neardata_var(N_PY,IM)  &
                   +(FAC2*DY*DZ       ) * epot1_var%neardata_var(N_PZ,IM))

                 TMP_FIELD(3) = TMP_FIELD(3) + (                           &
                    (FAC2*DZ*DX       ) * epot1_var%neardata_var(N_PX,IM)  &
                   +(FAC2*DZ*DY       ) * epot1_var%neardata_var(N_PY,IM)  &
                   +(FAC2*DZ*DZ - FAC1) * epot1_var%neardata_var(N_PZ,IM))
              ENDIF

!...........end of runinfo%pff
            ENDIF

!......end of loop over all mm-atoms A which are near qm-atom IQ: 
       ENDDO

       
!.....save energies and forces on qm atoms 
       EEXT = EEXT + TMP_POT * Q 

       FION(1,IA,IS)=FION(1,IA,IS)+Q*TMP_FIELD(1) 
       FION(2,IA,IS)=FION(2,IA,IS)+Q*TMP_FIELD(2)
       FION(3,IA,IS)=FION(3,IA,IS)+Q*TMP_FIELD(3)

!......save potential and field for transfer to iffi (only for output)
       QMPOTANDFLD(1,IQ) = TMP_POT
       QMPOTANDFLD(2,IQ) = TMP_FIELD(1)
       QMPOTANDFLD(3,IQ) = TMP_FIELD(2)
       QMPOTANDFLD(4,IQ) = TMP_FIELD(3)

!.......end of loop IA over NA qm-atoms of species IS:
        ENDDO
!.....end of loop IS over NSP atomspecies
      ENDDO

!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE
!     ==================================================================


END MODULE iffi_elstat_utils