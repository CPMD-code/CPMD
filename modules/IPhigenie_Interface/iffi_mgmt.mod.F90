MODULE iffi_mgmt_utils
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE efld,                            ONLY: extf,&
                                             textfld
  USE elct,                            ONLY: crge
  USE elstpo_utils,                    ONLY: elstpo
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE mm_dimmod,                       ONLY: clsaabox
  USE mm_extrap,                       ONLY: cold
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE pslo,                            ONLY: pslo_com
  USE purge_utils,                     ONLY: purge
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rrane_utils,                     ONLY: rrane
  USE setirec_utils,                   ONLY: write_irec
  USE store_types,                     ONLY: cprint,&
                                             restart1,&
                                             rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, maxsys, ncpw, ncpw, nkpt,  &
       parm, spar, parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrener
  USE wv30_utils,                      ONLY: zhwwf

  USE machine,                         ONLY:  m_walltime,m_flush
  
  USE iffi_types
  USE iffi_wfct
  USE mp_interface,                    ONLY: mp_bcast,&
                                                mp_recv,&
                                                mp_send,&
                                                mp_sync
  USE phfac_utils, ONLY : phfac

  USE iffi_elstat_utils,               ONLY: QMCOREFORCE
#include "sizeof.h"

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC  :: ITERATE_WFCT,iffi_give_scr_interface



CONTAINS
!....Interface to the MD simulation program IPHIGENIE
!.....
!..... ITERATE_WFCT             : KS iteration
!..... CALC_EXPORT              : calls export routines 
!..... WRITE_EXTF               : write cube files
!..... SAVEWF                   : save wavefunction (for normal mode analysis)
!..... TRANSLATEWF              : translate wfct on the grid
!..... MYPROPPT_IFFI            : properties calculation
!..... GIVE_SCR_MYPROPPT        : reserve scr memory
!..... GIVE_SCR_INTERFACE_WRITE :
!..... GIVE_SCR_INTERFACE       :
!***************************************************************************

!==================================================================
SUBROUTINE ITERATE_WFCT
!==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER(*), PARAMETER                  :: procedureN = 'ITERATE_WFCT'
      
      LOGICAL    FNOWF
      INTEGER    IREC(100)
      REAL*8     TIME1,TIME2,TIMEF
      REAL*8     ETOTO, &
                 DETOT,TCPU,dummy(1)
      REAL*8     eext_elec
      INTEGER ISUB

      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/


      EXTERNAL TIMEF

      CALL TISET(procedureN,ISUB)
      
        iteropt%nfi=iteropt%nfi+1
        ropt_mod%convwf=.FALSE.
        ropt_mod%sdiis =.TRUE.
        ropt_mod%spcg  =.TRUE.
        ropt_mod%modens=.FALSE.
        ropt_mod%engpri=.FALSE.
        ropt_mod%calste=.FALSE.

        etoto = 0.0_real_8

        CALL phfac(tau0)
!       TINLCC: True if Non Linear Core Correction is used for a species (nlcc.inc)
!       COPOT: COMPUTES THE POTENTIAL AND ENERGIES FOR EXCHANGE AND        
!       CORRELATION FOR THE CORE CHARGES (copot.F)                           
        IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)
        CALL mp_sync(parai%allgrp)
!
!.......Randomize Wavefunction
    IF (cntl%trane) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.8)')&
            ' RANDOMIZE INITIAL WAVEFUNCTION,  AMPLITUDE=',cntr%ampre
       CALL rrane(c0,c2,crge%n)
    ENDIF


!........translate the wavefunction if necessary
         IF (epot1_fix%boxtrans(1).NE.0.0D0.OR.epot1_fix%boxtrans(2).NE.0.0D0.OR.epot1_fix%boxtrans(3).NE.0.0D0) THEN
           IF (paral%io_parent)  WRITE(6,'(A)') 'INTERFACE| TRANSLATING WAVEFCT'

           CALL mp_sync(parai%allgrp)

           CALL translatewf(tau0,c0(:,:,1),c2,crge%n)  !FIXME c2 correct??

!          reset boxtrans vector
           epot1_fix%boxtrans(1)=0.0D0
           epot1_fix%boxtrans(2)=0.0D0
           epot1_fix%boxtrans(3)=0.0D0
        ENDIF

!.......if this is the first run after an integration step: extrapolate the wavefunction
        IF(cntl%textrap.AND.epot1_var%tdoextrap.AND..NOT.runinfo%normalmodeanalysis) THEN
            CALL mp_sync(parai%allgrp)
            IF (paral%io_parent) WRITE(6,*) 'DO WAVEFUNCTION EXTRAPOLATION..'
!           note that NFI is the total number of CPMD runs, not only
!           the runs after integration. But NFI is not used in EXTRAPWF() anyway.
            CALL extrapwf(infi,c0(:,:,1),scr,cold,crge%n,cnti%mextra)  !n was nstate
        ENDIF


! periodic output of density/wavefunction etc.
        IF(rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
          cntl%tepot=.TRUE. 
          IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
!.....set the box translation vector in order to align cpmd output with iffi trajectory
            clsaabox%mm_c_trans(1) = epot1_fix%boxoffs(1)
            clsaabox%mm_c_trans(2) = epot1_fix%boxoffs(2)
            clsaabox%mm_c_trans(3) = epot1_fix%boxoffs(3)
            CALL STOPGM('ITERATE_WFCT','RHOOUT NOT SUPPORTED IN INTERFACE MODE YET',0,"")!FIXME
            !CALL WRITE_EXTF 
 
!...........just to be sure
            clsaabox%mm_c_trans(1)=0.0d0
            clsaabox%mm_c_trans(2)=0.0d0
            clsaabox%mm_c_trans(3)=0.0d0
          ENDIF
        ENDIF

!.....for iffi normalmodeanalysis FIXME
!      IF (runinfo%normalmodeanalysis) THEN
!        IF (IFIRST.EQ.1) THEN
!         save the gas phase wavefunction that cpmd calculates at the beginning..
!          CALL SAVEWF(TAU0,C0,N,1)
!          IFIRST=0
!        ELSE
!         ..and restore it before each noma step to have identical starting conditions
!         (noma in cpmd does it this way, too)
!          CALL SAVEWF(TAU0,C0,N,0) 
!        ENDIF
!      ENDIF

!.................................................................
!.......ITERATE THE WAVE FUNCTION.................................
    fnowf=.FALSE.
    DO infi=1,cnti%nomore
       time1=m_walltime()
       IF (paral%parent) THEN
          ropt_mod%engpri=MOD(infi-1,cprint%iprint_step).EQ.0
       ELSE
          ropt_mod%engpri=.FALSE.
       ENDIF
       ! 
       ! .......update the wavefunctions 
       ! 
       CALL updwf(c0(:,:,1),c2,sc0(:,:,1),tau0,fion,pme,gde,vpp,&
            eigv,rhoe,psi,crge%n,fnowf,.TRUE.)
       ! 
       ! .......printout evolution of the progress
       ! 
       IF (paral%parent) THEN
          detot=ener_com%etot+ener_com%eext-etoto
          IF (infi.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF (ropt_mod%engpri) CALL wrener
          IF ((ropt_mod%engpri).AND.paral%io_parent)&
               WRITE(6,'(A)')&
               ' NFI     GEMAX     CNORM          ETOT         DETOT        TCPU'
          IF (paral%io_parent)&
               !WRITE(6,'(I4,F10.6,F10.6,F14.5,F14.6,F12.2) ')&
               WRITE(6,'(I4,1PE11.3,1PE12.3,0PF15.6,1PE13.3,0PF10.2) ')&
               infi,gemax,cnorm,ener_com%etot+ener_com%eext,detot,tcpu
          etoto=ener_com%etot+ener_com%eext
       ENDIF
       IF (ropt_mod%convwf.AND.fnowf) GOTO 100
       IF (ropt_mod%convwf) THEN
          ropt_mod%convwf=.FALSE.
          fnowf=.TRUE.
       ENDIF
!      if this is the MAXSTEPS'th iteration: throw an error!
       IF(INFI.EQ.cnti%nomore) GOTO 2
    ENDDO
100 CONTINUE
!.......end of ITERATE THE WAVE FUNCTION.................................


!..Calculate gradient
      IF (.NOT.epot1_var%tcalconlypfffield) THEN
        CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.TRUE.,.TRUE.)         
      ENDIF
!
    eext_elec = ener_com%eext
    
    IF (paral%parent) CALL QMCOREFORCE(tau0,ener_com%eext,fion)

    IF (paral%parent) THEN
       IF (paral%io_parent) THEN
            WRITE(6,*)'ELECTRONIC EXTERNAL ENERGY  =',eext_elec,' AU'
            WRITE(6,*)'IONIC EXTERNAL ENERGY       =',ener_com%eext-eext_elec,' AU'
            WRITE(6,*)'TOTAL EXTERNAL ENERGY       =',ener_com%eext,' AU'
            WRITE(6,*)'REAL TOTAL ENERGY           =',ener_com%eext+ener_com%etot,' AU'
       ENDIF
    ENDIF

                        



!.......Calculate the dipole moment if requested (PROPER)
!       and if pff: only if pffmode=onlywriteout or writeall
        prop1%espc=.FALSE.
        if (epot1_var%calcespcharges.EQ.1) prop1%espc=.TRUE.

        !FIXME properties
        !IF(PROPER.AND.(                                               &
        !               (.NOT.runinfo%pff).OR.                         &
        !               (.NOT.epot1_var%tcalconlypfffield).OR.         &
        !                (epot1_var%tcalconlypotfrompff.AND.           &
        !                 .NOT.epot1_var%tnopotandwfctupd.AND.         &
        !                 .NOT.epot1_var%tcalconlypfffield) ) ) THEN
        !    CALL MYPROPPT_IFFI(C0,TAU0,RHOE,PSI,SCR,LSCR)
        !   CALL myproppt(c0,tau0,rhoe,psi)
        !NDIF

! periodic output of density/wavefunction etc.
        IF(rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
          cntl%tepot=.TRUE. 
          IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
!.....set the box translation vector in order to align cpmd output with iffi trajectory
            clsaabox%mm_c_trans(1) = epot1_fix%boxoffs(1)
            clsaabox%mm_c_trans(2) = epot1_fix%boxoffs(2)
            clsaabox%mm_c_trans(3) = epot1_fix%boxoffs(3)

!          CALL RHOPRI(C0,TAU0,RHOE,PSI,SCR,LSCR,NSTATE,NKPNT)
           CALL rhopri(c0(:,:,1),tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
           
!...........just to be sure
            clsaabox%mm_c_trans(1)=0.0d0
            clsaabox%mm_c_trans(2)=0.0d0
            clsaabox%mm_c_trans(3)=0.0d0
          ENDIF
        ENDIF

!..impose constraints on the gradient
        CALL purge(tau0,fion)
                  
!..write the restart file
        if (epot1_var%twriterestartfile.AND..NOT.epot1_var%tcalconlypfffield) THEN
          CALL write_irec(irec)
          CALL zhwwf(2,irec,c0(:,:,1),c2,crge%n,dummy,tau0,velp,taup,iteropt%nfi)
        ENDIF

       GOTO 1
    2  CONTINUE


       PRINT*,'TOO MANY ITERATIONS NEEDED IN THIS STEP! ', &
               'SOMETHING IS WRONG OR INCREASE MAXSTEPS'
       CALL STOPGM('ITERATE_WFCT','TOO MANY ITERATIONS NEEDED',0,"")

    1 CONTINUE
      
      CALL TIHALT(procedureN,ISUB)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE
!     ==================================================================



! ==================================================================
SUBROUTINE iffi_give_scr_interface(linter,tag)
! ==--------------------------------------------------------------==
    
    IMPLICIT NONE
    
    INTEGER                                  :: linter
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcopot, lforcedr, linterw, &
                                                lmyproppt, lupd, nstate
    INTEGER                                  :: lrhopri                                              
    INTEGER                                  :: lesp, lmulliken, lrhoofr, &
                                                lrnlsm
                                            
   
    !from give_scr_interface(linterw,tag)
    nstate=crge%n
    lcopot=0
    lmyproppt=0

    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    CALL give_scr_updwf(lupd,tag,nstate,.FALSE.)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.TRUE.,.TRUE.)
    !FIXME properties IF (cntl%proper) CALL give_scr_myproppt(lmyproppt,tag)
    
    lrhopri=0
    !CALL give_scr_rhopri(lrhopri,tag,nstate) FIXME 
        
    !from give_scr_interface_write(linterw,tag)
    lrhoofr=0
    lrnlsm=0
    lesp=0
    lmulliken=0
    
    
    CALL give_scr_rhoofr(lrhoofr,tag)
    !FIXME ESP   CALL give_scr_espc(lesp,tag)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    
    
    linterw=MAX(lrhoofr,lesp,lrnlsm,lmulliken,lrhopri)    
    linter=MAX(lcopot,lupd,lforcedr,lmyproppt,linterw)
! ==--------------------------------------------------------------==
    RETURN
END SUBROUTINE iffi_give_scr_interface
! ==================================================================




#if 0
!.....write external potential to file MMPOT.*
!     ==================================================================
      SUBROUTINE WRITE_EXTF
!     ==--------------------------------------------------------------==

      IMPLICIT NONE
!       INCLUDE 'system.h'
!       INCLUDE 'ropt.inc'
!       INCLUDE 'store.inc'
!       INCLUDE 'efld.inc'
!       INCLUDE 'coor.inc'
!       INCLUDE 'cppt.inc'
!       INCLUDE 'iffi.inc'
!       INCLUDE 'spin.inc'

      INTEGER I,J,K,I1,I2,N1,N2,IG
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/
      CHARACTER  FILEN*30,CIPNUM*30

      REAL*8  EXTF(fpar%NNR1)
      DIMENSION EXTFCOMPLEX(fpar%NNR1)
      DIMENSION MYVTEMP(fpar%NNR1)


      IF(IFIRST.EQ.1) THEN
        IFIRST=0
        CALL MEMORY(IP_EXTFCOMPLEX,fpar%NNR1*2,'EXTFCOMPLEX')
        CALL MEMORY(IP_MYVTEMP,fpar%NNR1*2,'MYVTEMP')
      ENDIF

!     sign change to get phi(r)
!$OMP parallel do private(I)
      DO I=1,fpar%NNR1
        EXTFCOMPLEX(I)=DCMPLX(-EXTF(I),0.D0)
      ENDDO
!$OMP end parallel do

      CALL FWFFTN(EXTFCOMPLEX,.FALSE.)

!$OMP parallel do private(IG)
      DO IG=1,ncpw%NHG
        MYVTEMP(IG) = EXTFCOMPLEX(NZH(IG)) 
      ENDDO

      FILEN='MMPOT'
      CALL XSTRING(FILEN,I1,I2)
      
      IF(rout1%nrhoout.GT.0) THEN
        WRITE(CIPNUM,'(I14)') iteropt%nfi
        CALL XSTRING(CIPNUM,N1,N2)
        FILEN=FILEN(I1:I2) // '.' // CIPNUM(N1:N2)
        CALL XSTRING(FILEN,I1,I2)
      ENDIF

      CALL DENSTO(MYVTEMP,TAU0,FILEN(I1:I2))

!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
#endif


#if 0
!     IF SVE==1 wavefunction is saved to COLD(*,*,1,1) for a more accurate normal mode analysis with IPHIGENIE/CPMD
!     (mis)uses structures used for wavefunction extrapolation, which therefore has to be switched on (is not used then)
!     ==================================================================
      SUBROUTINE SAVEWF(TAU,C0,NSTATE,SVE)
!     ==--------------------------------------------------------------==
      IMPLICIT NONE 
!       INCLUDE 'system.h' 
!       INCLUDE 'parac.inc'
!       INCLUDE 'cppt.inc'
!       INCLUDE 'gvec.inc'
!       INCLUDE 'isos.inc'
!       INCLUDE 'mm_dim.inc'
!       INCLUDE 'iffi.inc'

      integer NSTATE
      COMPLEX*16 C0(ncpw%NGW,NSTATE,*) !changed for BSYMM 
      real*8 TAU(3,NAX,NSX)

!       INCLUDE 'elct.inc'      ! for N
!       INCLUDE 'mm_extrap.inc' ! for COLD()

      integer IG,IS
      INTEGER SVE
      
      IF(SVE.EQ.1) THEN
        if (PARENT) WRITE(6,'(A)') 'SAVING WAVEFUNCTION..'
        DO IS=1,NSTATE
          DO IG=1,ncpw%NGW
            COLD(IG,IS,1,1)=C0(IG,IS,1)
          ENDDO
        ENDDO
      ELSE
        if (PARENT) WRITE(6,'(A)') 'RESTORING WAVEFUNCTION..'
        DO IS=1,NSTATE
          DO IG=1,ncpw%NGW
            C0(IG,IS,1)=COLD(IG,IS,1,1)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
!     ==================================================================      
#endif



#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE4
#define __COLLAPSE3
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#define __COLLAPSE3 collapse(3)
#define __COLLAPSE4 collapse(4)
#endif

!     adapted from  mm_translate_qmmm in MM_Interface/mm_transl_c0.F
!     ==================================================================
      SUBROUTINE translatewf(TAU,C0,CM,NSTATE)
!     ==================================================================
      !use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
!     ==--------------------------------------------------------------==
!cmb - Revised on 25 October 2013
      use system 
      use parac
      use cppt
      use gvec
      use isos
      use mm_dimmod
      use geq0mod
      use elct      ! for N
      use mm_extrap ! for COLD()
      IMPLICIT NONE 

      integer NSTATE
      COMPLEX*16 C0(nkpt%NGWK,NSTATE,*) !changed for BSYMM 
      COMPLEX*16 CM(nkpt%NGWK,NSTATE,*) !changed for BSYMM

      real*8 TAU(3,maxsys%nax,maxsys%nsx)

      integer IG,IS,IH,IK
      integer NH1,NH2,NH3
      real*8 Gt1,Gt2,Gt3,rI,rJ,rK,gdotr,c_trans(3)
      complex*16, allocatable ::  eigr_l(:)

!      Translation vector in real space
       c_trans(1) = -epot1_fix%boxtrans(1)
       c_trans(2) = -epot1_fix%boxtrans(2)
       c_trans(3) = -epot1_fix%boxtrans(3)
       
      !clsaabox%mm_c_trans(1)= clsaabox%mm_c_trans(1)+c_trans(1)
      !clsaabox%mm_c_trans(2)= clsaabox%mm_c_trans(2)+c_trans(2)
      !clsaabox%mm_c_trans(3)= clsaabox%mm_c_trans(3)+c_trans(3)

      if (.FALSE.) print *,TAU(1,1,1) !to prevent compiler warning
      
      allocate (eigr_l(ncpw%NGW))
      NH1=spar%nr1s/2+1
      NH2=spar%nr2s/2+1
      NH3=spar%nr3s/2+1 

      if (paral%io_parent) WRITE(6,'(A,F21.14,TR1,F21.14,TR1,F21.14)')  &
       'TRANSLATING WFCT:',c_trans(1),c_trans(2),c_trans(3)
     
!$OMP parallel do private(IG,rI,rJ,rK,Gt1,Gt2,Gt3,gdotr)
      DO IG=1,ncpw%NGW
        rI=DBLE(INYH(1,IG)-NH1)
        rJ=DBLE(INYH(2,IG)-NH2)
        rK=DBLE(INYH(3,IG)-NH3)
        Gt1=rI*gvec_com%b1(1)+rJ*gvec_com%b2(1)+rK*gvec_com%b3(1)
        Gt2=rI*gvec_com%b1(2)+rJ*gvec_com%b2(2)+rK*gvec_com%b3(2)
        Gt3=rI*gvec_com%b1(3)+rJ*gvec_com%b2(3)+rK*gvec_com%b3(3)
        gdotr=parm%tpiba*(Gt1*c_trans(1) &
                   +Gt2*c_trans(2)      &
                   +Gt3*c_trans(3))
        eigr_l(IG)=DCMPLX(DCOS(gdotr),-DSIN(gdotr)) 
      ENDDO
!
!$OMP parallel do private(IS,IG) __COLLAPSE2
        DO IS=1,NSTATE
          DO IG=1,ncpw%NGW
            C0(IG,IS,1)=C0(IG,IS,1)*eigr_l(IG)
            CM(IG,IS,1)=CM(IG,IS,1)*eigr_l(IG)
          ENDDO
        ENDDO

      IF(CNTL%TMDBO.AND.CNTL%TEXTRAP) THEN
!     if using wavefunction extrapolation, apply the same
!     shift also to the wavefunction history.

!$OMP parallel do private(IH,IK,IS,IG) __COLLAPSE4
        DO IH=1,cnti%MEXTRA
          DO IK=1,nkpt%nkpnt
            DO IS=1,NSTATE
              DO IG=1,ncpw%NGW
                COLD(IG,IS,IK,IH)=COLD(IG,IS,IK,IH)*eigr_l(IG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      deallocate(eigr_l)

      RETURN
      END SUBROUTINE
      
END MODULE iffi_mgmt_utils
