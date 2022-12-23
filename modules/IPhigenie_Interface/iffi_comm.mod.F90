MODULE iffi_comm_utils
  USE cnst,                            ONLY: au_deb,&
                                             fbohr
  USE coor,                            ONLY: fion,&
                                             tau0
  USE dipomod,                         ONLY: moment
   USE efld,                            ONLY: textfld
  USE ener,                            ONLY: ener_com
  USE iffi_types                       
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE lodp,                            ONLY: dmomlo
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop1,&
                                             prop2
  USE ropt,                            ONLY: infi,&
                                             iteropt
  USE error_handling,                  ONLY: stopgm
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys,  nkpt,  &
       parm, spar, parap
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync, &
                                             mp_sum

#include "sizeof.h"

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: IFFI_READ, IFFI_WRITE,SHIFT_MULTIPOLES_WRAP,INHERIT_LOCAL_EXP_WRAP
  
CONTAINS
!....Interface to the MD simulation program IPHIGENIE
!    here the data exchange with iffi and mpi communication is performed
!..... IFFI_READ                : get interface data from iphigenie
!..... READ_RUNINFO             :  "
!..... READ_QMATOMDATA          :  "
!..... BCAST_IFFIDATA           : broadcast interface data to nodes
!..... IFFI_WRITE               : put interface data to iphigenie
!..... GATHER_IFFIDATA          : gather interface data at master

!..... INIT_VOXELS              : init voxels
!..... INHERIT_LOCAL_EXP_WRAP   : wrapper function
!..... SHIFT_MULTIPOLES_WRAP    : wrapper function
!***************************************************************************
#define DBG_IMPEXP 3
#define defPRINTONLYATSTARTUP .true.

!     ==================================================================
      SUBROUTINE IFFI_READ
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/

      IF(paral%parent.AND.runinfo%tverbose) WRITE(6,'(A)') "IFFI_READ"

      cnti%icmet = 7 !charge method for symmetric interface

      IF(paral%parent) THEN
          CALL READ_RUNINFO
          CALL READ_QMATOMDATA
      ENDIF
    
      !MPI communication: broadcast interface data from master
      CALL BCAST_IFFIDATA(IFIRST)

      IF(IFIRST.EQ.1) THEN
        IFIRST=0
        !every node inits the voxels once
        CALL INIT_VOXELS
      ENDIF


!.........print some information about what will be done in this run:
      IF (runinfo%pff.AND.paral%parent.AND.runinfo%tverbose) THEN
        WRITE(6,*)                                                                &
        'INTERFACE| runinfo%pff',runinfo%pff,cntr%tolog,epot1_var%tcalconlypotfrompff, &
                         epot1_var%tnopotandwfctupd,epot1_var%tcalconlypfffield

        IF (epot1_var%tcalconlypotfrompff) THEN
          WRITE(6,'(2A)') 'INTERFACE| INFO: UPDATE ONLY POTENTIAL FROM DIPOLES'
        ELSE
          WRITE(6,'(2A)') 'INTERFACE| INFO:  UPDATE POTENTIAL FROM CHARGES AND DIPOLES'
        ENDIF

        IF (epot1_var%tcalconlypfffield) THEN
          WRITE(6,'(2A)') 'INTERFACE| INFO: WRITE OUT ONLY E_PFF'
        ELSE
          WRITE(6,'(2A)') 'INTERFACE| INFO: WRITE OUT FULL ELECTROSTATICS'
        ENDIF

        IF (epot1_var%tnopotandwfctupd) THEN
          WRITE(6,'(2A)') 'INTERFACE| INFO: NO WAVEFCT UPDATE'
        ELSE
          WRITE(6,'(2A)') 'INTERFACE| INFO: DO WAVEFCT UPDATE'
        ENDIF
 
      ENDIF

      END SUBROUTINE
!     ==================================================================



!     ==================================================================
      SUBROUTINE BCAST_IFFIDATA(IFIRST)
!     ==--------------------------------------------------------------==

      INTEGER :: IFIRST

      !CALL mp_sync(parai%allgrp)

!.....now broadcast information to all processors
!     tau0
      CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)

!     cntr%tolog contains the convergence criterion for "convergence orbitals".
!     uncomment the 2 lines after that, if multiple criterions are needed
!     for any other criterion than "convergence orbitals" (cf. above)
      CALL mp_bcast(cntr%tolog,parai%source,parai%allgrp)
!     CALL MY_BCAST(CNTR1,MREAL*8,SOURCE,ALLGRP)
!     CALL MY_BCAST(CNTI1,MINTE*8/IRAT,SOURCE,ALLGRP)

      CALL mp_bcast(textfld,parai%source,parai%allgrp)
      
!.....now broadcast the import data
!     broadcast the variable contributions first, as they steer the boradcasting
!     of the fixed quantities
        IF(IFIRST.EQ.1) THEN
          IF (paral%parent.AND.runinfo%tverbose) WRITE(6,*) "INTERFACE | BROADCAST RUNINFO DATA"
!.........broadcasting runinfo
          CALL mp_bcast_byte(runinfo, size_in_bytes_of(runinfo),parai%source,parai%allgrp)
!         NOTE: IFIRST is set to 1 below 
        ENDIF
        
      IF (paral%parent.AND.runinfo%tverbose) WRITE(6,*) "INTERFACE | BROADCAST VAR DATA"
!.....broadcasting epot1_var
      CALL mp_bcast_byte(epot1_var, size_in_bytes_of(epot1_var),parai%source,parai%allgrp)
        
!.....broadcast fixed data only if not runinfo%pff or right after an integration step
      IF(.NOT.epot1_var%tcalconlypotfrompff) THEN     !no need to do this during pff iteration 

!.........broadcast lists only if updated
          IF (epot1_var%tupdateqmintlists) THEN
            IF (paral%parent.AND.runinfo%tverbose) WRITE(6,*) "INTERFACE | BROADCAST QMMM LISTS"
!...........broadcasting epot1_list
            CALL mp_bcast_byte(epot1_list, size_in_bytes_of(epot1_list),parai%source,parai%allgrp)
          ENDIF
!.........broadcasting epot1_fix
          IF (paral%parent.AND.runinfo%tverbose) WRITE(6,*) "INTERFACE | BROADCAST FIXED DATA"
          CALL mp_bcast_byte(epot1_fix, size_in_bytes_of(epot1_fix),parai%source,parai%allgrp)
      ENDIF

!     ==--------------------------------------------------------------==
      END SUBROUTINE
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE READ_RUNINFO
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 CONVCRIT
      REAL*8 BOXDUM_TMP(3)
      INTEGER SKIPUPDATEGRIDPART,FIRSTRUNAFTERINT,IWRITERESTARTFILE
      INTEGER ICALCONLYPOTFROMPFF,INOPOTANDWFCTUPD,ICALCONLYPFFFIELD
      INTEGER IUPDATEQMWWLISTS,INOMA,IUPDATERGYR,IMEANFIELDMODE
      INTEGER IVOXNUM(3),IPFF,ISAMMTYPE,INRQMATOMSIFFI,IUPDATEMEANFIELD
      REAL*8  VERBOSE
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/


      IF(.NOT.paral%parent) RETURN

!.....set default values for every step

      epot1_var%tcalconlypotfrompff = .FALSE.
      epot1_var%tnopotandwfctupd = .FALSE.
      epot1_var%tcalconlypfffield = .FALSE.
      epot1_var%tgridpart = .TRUE.
      epot1_var%tdoextrap=.FALSE.
      epot1_var%tupdateqmintlists=.TRUE.
      epot1_var%calcespcharges = 0
      epot1_var%tupdatergyr = .FALSE.
      epot1_var%updatemeanfield = .FALSE.

      IF (IFIRST.EQ.1) THEN
!.......set default values once
        TEXTFLD=.TRUE. !is always true and does not change!
        runinfo%normalmodeanalysis=.FALSE.
        runinfo%meanfieldmode=.FALSE.
        epot1_var%nmean=0

        call get_runinfo_first_int(IPFF,ISAMMTYPE,IVOXNUM,INRQMATOMSIFFI,INOMA,IMEANFIELDMODE)

       IF (INRQMATOMSIFFI.NE.runinfo%nrqmatoms) THEN
         WRITE(6,*) 'WRONG NR OF QMATOMS IN CPMD!'
         WRITE(6,*) 'IFFI:',INRQMATOMSIFFI,'CPMD:',runinfo%nrqmatoms
         WRITE(6,*)'Maybe you have forgotten to include all atoms types'
         WRITE(6,*) 'in you template file.'

         CALL STOPGM('READ QMATOMDATA','WRONG NR OF QMATOMS',0,"")
       endif

        IF (INOMA.GT.0) THEN
          runinfo%normalmodeanalysis=.TRUE.
!.........wfct extrapolation needs to be switched on to have the memory available
          IF (.NOT.cntl%textrap) THEN
            WRITE(6,*) 'runinfo%normalmodeanalysis-MODE NEEDS'
            WRITE(6,*) 'WAVEFUNCTION EXTRAPOLATION TO BE SWITCHED ON.'
            CALL STOPGM('READ QMATOMDATA','runinfo%normalmodeanalysis-MODE NEEDS WAVEFUNCTION EXTRAPOLATION',0,"")
          ENDIF
        ENDIF

        if (IPFF.GT.0) runinfo%pff = .TRUE.

        runinfo%sammtype  = ISAMMTYPE
        runinfo%voxnum(1) = IVOXNUM(1)
        runinfo%voxnum(2) = IVOXNUM(2)
        runinfo%voxnum(3) = IVOXNUM(3)
        
        IF (IMEANFIELDMODE.EQ.1) runinfo%meanfieldmode=.TRUE.

!......ole parameter
        call get_runinfo_first_dbl3(runinfo%oleeps,runinfo%olensig,VERBOSE)
        runinfo%olensig=runinfo%olensig/DSQRT(2.0D0) !to divide out factor of sqrt(2) that all sigmas carry
        
        IF (VERBOSE.GT.0.0D0) THEN
          runinfo%tverbose=.TRUE.
          WRITE(6,*) 'VERBOSE OUTPUT ON!'
        ENDIF
         
!.......box dimensions
!.......FROMIFFI:
        call get_boxdum(BOXDUM_TMP)
        runinfo%boxdum(1) = 0.0D0
        runinfo%boxdum(2) = BOXDUM_TMP(1)
        runinfo%boxdum(3) = 0.0D0
        runinfo%boxdum(4) = BOXDUM_TMP(2)
        runinfo%boxdum(5) = 0.0D0
        runinfo%boxdum(6) = BOXDUM_TMP(3)

        IFIRST = 0
        
         WRITE (6,*) "INTERFACE | RUNINFO",                                             &
           runinfo%normalmodeanalysis,runinfo%pff,runinfo%sammtype,                     &
           runinfo%voxnum(1),runinfo%voxnum(2),runinfo%voxnum(3),runinfo%meanfieldmode, &
           runinfo%oleeps,runinfo%olensig,                                              & 
           runinfo%boxdum(2),runinfo%boxdum(4),runinfo%boxdum(6)

      ENDIF !end if FIRST

!.....FROMIFFI:
!.....important: get_runinfo_int needs to be called as the first get_* function, for after that is it known what to calculate
      call get_runinfo_int(FIRSTRUNAFTERINT,                  &
      &                 ICALCONLYPOTFROMPFF,INOPOTANDWFCTUPD, &
                        ICALCONLYPFFFIELD,                    &
                        SKIPUPDATEGRIDPART,                   &
                        epot1_var%calcespcharges,             &
                        IWRITERESTARTFILE,                    &
                        IUPDATEQMWWLISTS,                     &
                        IUPDATERGYR,                          &
                        IUPDATEMEANFIELD)

      call get_runinfo_dbl(CONVCRIT)

!.....convcrit
!.....write convcrit in the "CONVERGENCE ORBITALS" variable
      cntr%tolog = CONVCRIT
!     note: this implies that "CONVERGENCE ORBITALS" is used.
!           for different criteria types search for "Convergence criteria" in control.F .
!           if adding new types, don't forget to broadcast


!.....FROMIFFI:
      call get_boxtrans(epot1_fix%boxtrans)
!.....boxtrans
      IF (runinfo%tverbose)                                                    &
            WRITE(6,*) 'INTERFACE| SETTING BOX TRANSLATION VECTOR ',           &
            epot1_fix%boxtrans(1),epot1_fix%boxtrans(2),epot1_fix%boxtrans(3)

!.....FROMIFFI:
      call get_boxoffs(epot1_fix%boxoffs)

!     decide whether to extrapolate or not if runinfo%pff is done according to PFFSTATUS 
      IF (FIRSTRUNAFTERINT.EQ.1) THEN
          epot1_var%tdoextrap=.TRUE.
          IF (runinfo%tverbose) WRITE(6,*)  'INTERFACE| FIRST RUN AFTER INT'
      ENDIF

      IF (IWRITERESTARTFILE.EQ.1) THEN
          epot1_var%twriterestartfile=.TRUE.
          IF (runinfo%tverbose) WRITE(6,*)  'INTERFACE| writing RESTART file'
      ENDIF

!.....determine what to calculate in this run
      IF (runinfo%pff) THEN
       IF (ICALCONLYPOTFROMPFF.EQ.1) epot1_var%tcalconlypotfrompff = .TRUE.
       IF (INOPOTANDWFCTUPD.EQ.1)    epot1_var%tnopotandwfctupd    = .TRUE.
       IF (ICALCONLYPFFFIELD.EQ.1)   epot1_var%tcalconlypfffield   = .TRUE.
         
       IF (runinfo%tverbose) THEN
         IF(epot1_var%tcalconlypotfrompff) THEN
           WRITE(6,'(A)',advance='no') 'INTERFACE| LIGHT IMPORT'
         ELSE
           WRITE(6,'(A)',advance='no') 'INTERFACE| FULL IMPORT'
         ENDIF
         IF(epot1_var%tcalconlypfffield) THEN
          WRITE(6,'(A)',advance='yes') ' / LIGHT EXPORT'
         ELSE
          WRITE(6,'(A)',advance='yes') ' / FULL EXPORT'
         ENDIF
       ENDIF
       IF(epot1_var%tnopotandwfctupd) THEN
         WRITE(6,*) 'INTERFACE| NO WAVEFUNCTION UPDATE'
         call STOPGM("READ_RUNINFO","DISCONTINUED",0,"")
!       ELSE
!        WRITE(6,*) 'INTERFACE| WAVEFUNCTION UPDATE'
       ENDIF
      ENDIF !end of if runinfo%pff

      
      IF (SKIPUPDATEGRIDPART.EQ.1.AND.iteropt%nfi.GT.0) THEN
        epot1_var%tgridpart = .FALSE.
        IF (runinfo%tverbose) WRITE(6,*) 'INTERFACE| NO GRIDPART UPDATE '
      ELSE
        IF (runinfo%tverbose) WRITE(6,*) 'INTERFACE| DO GRIDPART UPDATE '
      ENDIF
      IF (IUPDATERGYR.EQ.1.AND..NOT.epot1_var%tcalconlypotfrompff) THEN
        IF (runinfo%tverbose) WRITE(6,*) 'INTERFACE| UPDATE RGYR'
        epot1_var%tupdatergyr=.TRUE.
      ENDIF


      IF (IUPDATEQMWWLISTS.EQ.1.OR.iteropt%nfi.EQ.0) THEN
        IF (runinfo%tverbose) WRITE(6,*) 'INTERFACE| DO QMINTLISTS UPDATE '
      ELSE
        epot1_var%tupdateqmintlists = .FALSE.
        IF (runinfo%tverbose) WRITE(6,*) 'INTERFACE| NO QMINTLISTS UPDATE'
      ENDIF
      IF (epot1_var%calcespcharges.EQ.1.AND.runinfo%tverbose) WRITE(6,*) 'INTERFACE| CALC ESP CHARGES'
      
   
      IF (runinfo%meanfieldmode) THEN
        IF(IUPDATEMEANFIELD.GT.0) THEN
          IF (FIRSTRUNAFTERINT.EQ.1) THEN
              epot1_var%nmean=epot1_var%nmean+1 !increment only at beginning of integration step
          ENDIF
          IF (epot1_var%nmean.GT.0) THEN
            epot1_var%updatemeanfield=.TRUE.
            IF (runinfo%tverbose)  WRITE (6,*) 'INTERFACE| UPDATING EXTF AND RHOE',epot1_var%nmean
          ENDIF
          
          ener_com%eext=0.0D0
        ENDIF
        epot1_var%tupdatergyr=epot1_var%updatemeanfield
        epot1_var%tgridpart = epot1_var%updatemeanfield
        
        IF (epot1_var%nmean.EQ.0) THEN 
          epot1_var%tgridpart=.TRUE.
          epot1_var%tupdatergyr=.TRUE.
        ENDIF
      ENDIF
      

       WRITE (6,*) "INTERFACE | INFO ",                                                                &
       cntr%tolog,epot1_fix%boxtrans(1),epot1_fix%boxtrans(2),epot1_fix%boxtrans(3),                   &
       epot1_fix%boxoffs(1),epot1_fix%boxoffs(2),epot1_fix%boxoffs(3),                                 &
       FIRSTRUNAFTERINT,epot1_var%tdoextrap,epot1_var%twriterestartfile,                                 &
       epot1_var%tcalconlypotfrompff, epot1_var%tnopotandwfctupd, epot1_var%tcalconlypfffield,         &
       epot1_var%tgridpart,epot1_var%tupdatergyr,epot1_var%tupdateqmintlists,epot1_var%calcespcharges, &
       runinfo%meanfieldmode, IUPDATEMEANFIELD, epot1_var%updatemeanfield

      END SUBROUTINE
!     ==================================================================


!     ==================================================================
      SUBROUTINE READ_QMATOMDATA
!     ==--------------------------------------------------------------==
      IMPLICIT NONE

      REAL*8 QMATOMPOS(3)
      REAL*8 LOCEXP_TMP(dimlocexp)
      REAL*8 NEARDATA_TMP(sizeneardata+sizeneardata_var)
      REAL*8 KCMOL_H
      INTEGER I,J,K,LOCID_TMP,NINCL_TMP
      INTEGER LINCL_TMP(maxnear)
      INTEGER IA,IS
#if DBG_IMPEXP
!     counter for printing import/export only at startup
      INTEGER IFIRSTDBG
      SAVE    IFIRSTDBG
      DATA    IFIRSTDBG /DBG_IMPEXP/
#endif
      INTEGER IFIRST
      SAVE    IFIRST
      DATA    IFIRST /1/

      IF(.NOT.paral%parent) RETURN

!......change this if pseudoatoms are used
       epot1_list%nloc=runinfo%nrqmatoms

       I=0
       DO IS=1,ions1%nsp
         DO IA=1,ions0%na(IS)
            I=I+1
            call get_qmatom_position(I,QMATOMPOS)

!...........security check
!           check if the qm atom coordinates written to the cpmd.inp file
!           are identical to those in the iffi-structure to validate the mapping
            if (IFIRST.EQ.1.AND.                                          &
                (ABS(TAU0(1,IA,IS)-QMATOMPOS(1)) > 0.0000001D0.OR.        &
                 ABS(TAU0(2,IA,IS)-QMATOMPOS(2)) > 0.0000001D0.OR.        &
                 ABS(TAU0(3,IA,IS)-QMATOMPOS(3)) > 0.0000001D0)) THEN     
              WRITE (*,'(A,F12.6,F12.6,F12.6)') "TAU0  ", TAU0(1,IA,IS),  &
                       TAU0(2,IA,IS),TAU0(3,IA,IS)                        
              WRITE (*,'(A,F12.6,F12.6,F12.6)') "QMPOS ", QMATOMPOS(1),   &
                       QMATOMPOS(2),QMATOMPOS(3)                          
              CALL m_flush(6)
              CALL STOPGM('READ QMATOMDATA','Wrong qmmap!',0,"")
            endif
            IFIRST=0
!...........end security check

!.......now the mapping between the (NA,NSP)-counting and the I=1,epot1_list%nloc counting is being fixed
!       set coordinates that cpmd uses
        TAU0(1,IA,IS)=QMATOMPOS(1)
        TAU0(2,IA,IS)=QMATOMPOS(2)
        TAU0(3,IA,IS)=QMATOMPOS(3)
!       set coordinates that interface uses for computation of electrostatics
        epot1_fix%koloc(1,I)=QMATOMPOS(1)
        epot1_fix%koloc(2,I)=QMATOMPOS(2)
        epot1_fix%koloc(3,I)=QMATOMPOS(3)

!.......EXPANSION COEFFICIENTS...................................
!......FROMIFFI:
        call get_qmatom_expansion(I,LOCID_TMP,NINCL_TMP,LOCEXP_TMP)

        epot1_list%locid(I)=LOCID_TMP
!       length of near inclusion list   
        epot1_list%nincl(I)=NINCL_TMP

!       SAMM coefficients
        DO J=1,dimlocexp
            epot1_var%locexp(J,I)=LOCEXP_TMP(J)
        END DO

!.......security check
        IF (epot1_list%nincl(I).GE.maxnear) THEN
           WRITE(6,*) 'INTERFACE| Error: maxnear not big enough!'
           WRITE(6,*) '           Increase this value in iffi.inc'
           CALL STOPGM('INTERFACE','Too many near atoms.',0,"")
        ENDIF

        IF (.NOT.epot1_var%tcalconlypotfrompff) THEN     !no need to do this during pff iteration
!.........IDS OF NEAR MM ATOMS OF THAT QM ATOM...............................
!.........FROMIFFI:
          call get_qmatom_nearlist(I,LINCL_TMP)

!.........epot1_list%lincl(I,J) gives the position of the J'th NEAR atom of
!                  QM atom I in the NEAR list  epot1_fix%neardata(?,J)
!                  so we loop for each QM atom LI 
!                  over N=1,epot1_list%nincl(LI) and get  IT = epot1_list%lincl(LI,N) and e.g. epot1_fix%neardata(PCHARGE,IT)
          DO J=1,epot1_list%nincl(I)
            epot1_list%lincl(I,J)=LINCL_TMP(J)
          ENDDO
        ENDIF



#if DBG_IMPEXP
      IF (runinfo%tverbose.AND.IFIRSTDBG.GT.0.OR..NOT.defPRINTONLYATSTARTUP) THEN
!    myprint LOCAL EXPANSIONS
        WRITE(*,'(A,I6,F12.8,TR1,F12.8,TR1,F12.8,TR1,I6)') &
             "QMAT",epot1_list%locid(I),epot1_fix%koloc(1,I),epot1_fix%koloc(2,I),epot1_fix%koloc(3,I),epot1_list%nincl(I)
        WRITE(*,'(25(F21.13,TR1))')  (epot1_var%locexp(J,I),J=1,25)

        DO J=1,epot1_list%nincl(I)
          if (MOD(J,10).EQ.0) THEN
             WRITE(6,'(I6)',advance='yes') epot1_list%lincl(I,J)
          else
             WRITE(6,'(I6)',advance='no') epot1_list%lincl(I,J)
          endif  
        ENDDO

        IF(paral%parent) THEN 
        print *," "        
        ENDIF
      ENDIF
#endif       

         ENDDO !IA
       ENDDO   !IS
!      end of loop over DFT atoms


!.....LIST OF NEAR MM ATOMS...................................
!.....FROMIFFI:
       IF (.NOT.epot1_var%tcalconlypotfrompff) THEN     !no need to do this during pff iteration
         call get_nearatoms_len(epot1_list%nnear)

         IF (epot1_list%nnear.GE.maxnear) THEN
            WRITE(6,*) 'INTERFACE| Error: maxnear not big enough!'
            WRITE(6,*) '           Increase this value in iffi.inc' 
            CALL STOPGM('INTERFACE','Too many near atoms.',0,"")
         ENDIF

       ENDIF

       DO J=1,epot1_list%nnear
!........FROMIFFI:         
         IF (.NOT.epot1_var%tcalconlypotfrompff) THEN     !no need to do this during pff iteration
           call get_mergednearlistitem(J,epot1_list%idnear(J),NEARDATA_TMP)
           DO K=1,sizeneardata
             epot1_fix%neardata(K,J) = NEARDATA_TMP(K)
           ENDDO
         ENDIF

         call get_mergednearlistitem_dipole(J,NEARDATA_TMP)         

         epot1_var%neardata_var(N_PX,J) = NEARDATA_TMP(1)
         epot1_var%neardata_var(N_PY,J) = NEARDATA_TMP(2)
         epot1_var%neardata_var(N_PZ,J) = NEARDATA_TMP(3)


         KCMOL_H = 0.00159360144639910963D0*332.0636D0*FBOHR
!                   = KCAL_TO_HART * f_elstat / BOHR_TO_ANGS =  1.0003461691537221959

         IF (.NOT.epot1_var%tcalconlypotfrompff) THEN     !no need to do this during pff iteration
!........Multiply conversion factor from kcal/mol to Hartree
           epot1_fix%neardata(PCHARGE,J)    = epot1_fix%neardata(PCHARGE,J)*KCMOL_H
           epot1_fix%neardata(CORECHARGE,J) = epot1_fix%neardata(CORECHARGE,J)*KCMOL_H 

!          to avoid multiplication with sqrt(2) in all formulas later
           epot1_fix%neardata(SIGCHARGE,J)=epot1_fix%neardata(SIGCHARGE,J)*DSQRT(2.0D0)
           epot1_fix%neardata(SIGDIPOLE,J)=epot1_fix%neardata(SIGDIPOLE,J)*DSQRT(2.0D0)
         ENDIF


!        p = q*d            
         epot1_var%neardata_var(N_PX,J) = epot1_var%neardata_var(N_PX,J)*KCMOL_H
         epot1_var%neardata_var(N_PY,J) = epot1_var%neardata_var(N_PY,J)*KCMOL_H
         epot1_var%neardata_var(N_PZ,J) = epot1_var%neardata_var(N_PZ,J)*KCMOL_H


#if DBG_IMPEXP
      IF (runinfo%tverbose.AND.IFIRSTDBG.GT.0.OR..NOT.defPRINTONLYATSTARTUP) THEN     
        IF(paral%parent) THEN   
          WRITE(*,'(I8,TR2,F11.5,F11.5,F11.5,F11.5,F11.5,                                           &
                  &F11.5,F11.5,F11.5,F21.13,F21.13,F21.13,F11.5)' )                                 &
               epot1_list%idnear(J),epot1_fix%neardata(PCHARGE,J),                                  &
               epot1_fix%neardata(CO_WIDTH,J),epot1_fix%neardata(CO_VALUE,J),                       &
               epot1_fix%neardata(CORECHARGE,J),epot1_fix%neardata(SIGCHARGE,J),                    &
               epot1_fix%neardata(N_KX,J),epot1_fix%neardata(N_KY,J),epot1_fix%neardata(N_KZ,J),    &        
               epot1_var%neardata_var(N_PX,J),epot1_var%neardata_var(N_PY,J),                                                 &
               epot1_var%neardata_var(N_PZ,J),epot1_fix%neardata(SIGDIPOLE,J)      
        ENDIF
      ENDIF
#endif

      ENDDO
!     end loop over all NEAR mm atoms

#if DBG_IMPEXP
!     set counter for debug information
      IF (IFIRSTDBG.GT.0) THEN
           IFIRSTDBG = IFIRSTDBG - 1
      ENDIF
#endif

      END SUBROUTINE
!     ==================================================================




!     ==================================================================
      SUBROUTINE IFFI_WRITE
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 QMATOMFOR(3),QMATOMPOS(3),LATTICEVEC(3)
      REAL*8 MULTEXP_TMP(dimmultexp),QMPOTANDFLD_TMP(4)
      REAL*8 MMPOTEXP_TMP(dimpotexp+DIMPOTEXP_VAR)
      REAL*8 DIPOLE_TMP(3),ECHRG
      REAL*8 RGYREXP_TMP(dimrgyrexp)

      INTEGER I,J,k,L,LOCID_TMP,IA,IS
#if DBG_IMPEXP
!     counter for printing import/export only at startup
      INTEGER IFIRSTDBG
      SAVE    IFIRSTDBG
      DATA    IFIRSTDBG /DBG_IMPEXP/
#endif

      !MPI communication: gather interface data at master
      CALL GATHER_IFFIDATA

      IF(.NOT.paral%parent) RETURN

      if (paral%parent.AND.runinfo%tverbose) WRITE(6,'(A)') "IFFI_WRITE"
      
!.......dipole moment calculated by properties
      IF (.NOT.epot1_var%tcalconlypfffield) THEN  !no need to do until end of pff iteration 
         IF(prop1%locd) THEN
           DIPOLE_TMP(1)=au_deb*dmomlo(1,1)
           DIPOLE_TMP(2)=au_deb*dmomlo(2,1)
           DIPOLE_TMP(3)=au_deb*dmomlo(3,1)
         ELSE IF(prop1%ldip) THEN
           DIPOLE_TMP(1)=au_deb*moment%dmom(1)
           DIPOLE_TMP(2)=au_deb*moment%dmom(2)
           DIPOLE_TMP(3)=au_deb*moment%dmom(3)
         ELSE
!          these values trigger that iffi uses the dipole calculated from multipoles expansions at qm atoms
           DIPOLE_TMP(1)=-99999.9
           DIPOLE_TMP(2)=-99999.9
           DIPOLE_TMP(3)=-99999.9
         ENDIF

!........TOIFFI:
         call set_qm_energies(ener_com%etot,ener_com%eext)

!........TOIFFI:
         call set_qmfrag_dipole(DIPOLE_TMP)
      ENDIF

!.....lattice vector
      LATTICEVEC(1)=(runinfo%boxdum(2)-runinfo%boxdum(1))/DBLE(spar%nr1s)
      LATTICEVEC(2)=(runinfo%boxdum(4)-runinfo%boxdum(3))/DBLE(spar%nr2s)
      LATTICEVEC(3)=(runinfo%boxdum(6)-runinfo%boxdum(5))/DBLE(spar%nr3s)
!.....TOIFFI:
      call set_latticevec(LATTICEVEC)

#if DBG_IMPEXP
        IF (runinfo%tverbose.AND.IFIRSTDBG.GT.0) THEN
            WRITE(*,'(A,F12.8,TR1,F12.8,TR1,F12.8)') "LATTICEVEC",LATTICEVEC(1),LATTICEVEC(2),LATTICEVEC(3)   
        ENDIF
#endif


       I=0
       DO IS=1,ions1%nsp
         DO IA=1,ions0%na(IS)
            I=I+1
!.......POSITION and FORCES
!           note: positions are only passed to enable to check if the data is stored
!           in the right structure
             IF (.NOT.epot1_var%tcalconlypfffield) THEN  !no need to do until end of pff iteration 
              QMATOMPOS(1)=TAU0(1,IA,IS)
              QMATOMPOS(2)=TAU0(2,IA,IS)
              QMATOMPOS(3)=TAU0(3,IA,IS)

!.............FORCES................................................
              QMATOMFOR(1)=FION(1,IA,IS)
              QMATOMFOR(2)=FION(2,IA,IS)
              QMATOMFOR(3)=FION(3,IA,IS)
!.............TOIFFI:
              call set_qmatom_forces(I,QMATOMFOR,QMATOMPOS)
            ENDIF

!...........MULTIPOLE MOMENTS...................................
            LOCID_TMP=epot1_list%locid(I)
            DO J=1,dimmultexp
              MULTEXP_TMP(J)=ATMULTEXPGLOB(J,I)
            ENDDO
!...........TOIFFI:
            call set_qmatom_multipoles(I,LOCID_TMP,MULTEXP_TMP)

!...........POTENTIAL AND FIELD...................................
            DO J=1,4
              QMPOTANDFLD_TMP(J) = QMPOTANDFLD(J,I)
            ENDDO

!...........TOIFFI:
            ECHRG=0.0D0 !FIXME ESP
            call set_qmatom_potandfield(I,LOCID_TMP,QMPOTANDFLD_TMP,ECHRG)

!...........RGYR EXPANSION...................................
            DO J=1,dimrgyrexp
              RGYREXP_TMP(J)=ATRGYREXPGLOB(J,I)
            ENDDO
!...........TOIFFI:
            call set_qmatom_rgyrexp(I,epot1_list%locid(I),RGYREXP_TMP,runinfo%oleradius*runinfo%oleeps/2.0D0)

#if DBG_IMPEXP        
      IF (runinfo%tverbose.AND.IFIRSTDBG.GT.0.OR..NOT.defPRINTONLYATSTARTUP) THEN
!       myprint multipoles
        IF(paral%parent) THEN  
        WRITE(*,'(I6,F12.8,TR1,F12.8,TR1,F12.8,TR1,F12.8,TR1,   &
                   &F12.8,TR1,F12.8,TR1,A)') LOCID_TMP,          &
                TAU0(1,IA,IS), TAU0(2,IA,IS), TAU0(3,IA,IS),    &
                FION(1,IA,IS), FION(2,IA,IS), FION(3,IA,IS),    &
        "(VOMULTEXP)"
        WRITE(*,'(25(F12.8,TR1))')  (ATMULTEXPGLOB(J,I),J=1,25)
        ENDIF
      ENDIF
#endif
         ENDDO !IA
       ENDDO !IS
!      end of loop over DFT atoms

!.......ELECTROSTATICS AT NEAR MM ATOMS...................................
!.......loop over all MM atoms in NEARLIST
        DO J=1,epot1_list%nnear
           MMPOTEXP_TMP(dimpotexp+MEXPFF)=MMPOTEXPGLOB_VAR(MEXPFF,J)
           MMPOTEXP_TMP(dimpotexp+MEYPFF)=MMPOTEXPGLOB_VAR(MEYPFF,J)
           MMPOTEXP_TMP(dimpotexp+MEZPFF)=MMPOTEXPGLOB_VAR(MEZPFF,J)

          IF (.NOT.epot1_var%tcalconlypfffield) THEN  !no need to do until end of pff iteration 
            DO K=1,dimpotexp
              MMPOTEXP_TMP(K)=MMPOTEXPGLOB(K,J)
            ENDDO
!           set whole elstat
            call set_mmatom_elstat(J,epot1_list%idnear(J),MMPOTEXP_TMP,0) ! 3rd argument 0 means: whole elstat
          ELSE
!           set only elstat at dipole
            call set_mmatom_elstat(J,epot1_list%idnear(J),MMPOTEXP_TMP,1) ! 3rd argument 1 means: only elstat at dipole
          ENDIF

        ENDDO

#if DBG_IMPEXP
      IF (runinfo%tverbose.AND.IFIRSTDBG.GT.0.OR..NOT.defPRINTONLYATSTARTUP) THEN
        DO J=1,epot1_list%nnear
!      myprint NEAR expansion
        IF(paral%parent) THEN    
          WRITE(*,'(I6,F12.8,TR1,F12.8,TR1,F12.8,TR1,A)')                                        &
            J,epot1_fix%neardata(N_KX,J),epot1_fix%neardata(N_KY,J),epot1_fix%neardata(N_KZ,J),    &
            "(MMPOTEXP)"
          WRITE(*,'(25(F12.8,TR1))')  (MMPOTEXPGLOB(L,J),L=1,10),(MMPOTEXPGLOB_VAR(L,J),L=1,3) 
        ENDIF
        ENDDO
      ENDIF

      
!     set counter for debug information
      IF (IFIRSTDBG.GT.0) THEN
         IFIRSTDBG = IFIRSTDBG - 1
      ENDIF
#endif
!     ==--------------------------------------------------------------==
      END SUBROUTINE
!     ==================================================================


!     ==================================================================
      SUBROUTINE GATHER_IFFIDATA
!     ==--------------------------------------------------------------==
       
      IMPLICIT NONE
      INTEGER :: msglen
      
      !no need to do this until end of pff iteration
      IF (.NOT.epot1_var%tcalconlypfffield) THEN     
        msglen=dimpotexp * epot1_list%nnear
        CALL mp_sum(MMPOTEXP(:,:,1),MMPOTEXPGLOB,msglen,parai%source,parai%allgrp)
        !(data_in,data_out,n, root, comm)
      ENDIF
      
      msglen=DIMPOTEXP_VAR * epot1_list%nnear
      CALL mp_sum(MMPOTEXP_VAR(:,:,1),MMPOTEXPGLOB_VAR,msglen,parai%source,parai%allgrp)


      msglen=dimmultexp * epot1_list%nloc
      CALL mp_sum(ATMULTEXP(:,:,1),ATMULTEXPGLOB,msglen,parai%source,parai%allgrp)

      IF (epot1_var%tupdatergyr) THEN
        msglen=dimrgyrexp * maxind
        CALL mp_sum(ATRGYREXP(:,:,1),ATRGYREXPGLOB,msglen,parai%source,parai%allgrp)
      ENDIF
!     ==--------------------------------------------------------------==
      END SUBROUTINE
!     ==================================================================


!     ==================================================================
      SUBROUTINE INIT_VOXELS
!     ==--------------------------------------------------------------==
      USE kinds, ONLY: real_8
      
      USE system
      
      IMPLICIT NONE
      CHARACTER(len=*), PARAMETER              :: procedureN = "INIT_VOXELS"
      
      INTEGER V,VG,VO, I,J,K
      REAL*8 DK,DJ,DI, KO_X,KO_Y,KO_Z, DX,DY,DZ
      REAL*8 D,DMIN
      REAL*8 LBOX(3),LTARGET,DELTAL,DELTALMIN,LVOXTMP,LVOX(3)
      INTEGER NVOX(3)
      LOGICAL TFOUND
      REAL*8 VOXRGYR

      INTEGER NVOXTOT

      INTEGER GPOINTS(3)
      INTEGER ierr

      INTEGER ikr1
      INTEGER idx !to access extf
      REAL(real_8), ALLOCATABLE  :: KOVOXG(:,:)
      INTEGER, ALLOCATABLE  :: MYVOX(:),MYVG2VO(:)
      
      MYNVOX  = 0 !number of voxels on this node
      
      IKR1 = MOD(parm%nr1+1,2)   ! !added: biswas

!     box extension
      LBOX(1)=runinfo%boxdum(2)-runinfo%boxdum(1)
      LBOX(2)=runinfo%boxdum(4)-runinfo%boxdum(3)
      LBOX(3)=runinfo%boxdum(6)-runinfo%boxdum(5)

!     grid constant
      DI=LBOX(1)/DBLE(spar%nr1s)
      DJ=LBOX(2)/DBLE(spar%nr2s)
      DK=LBOX(3)/DBLE(spar%nr3s)

      
      !IF (paral%parent) THEN
      ! print *,"GRIDINFO NRxS :",spar%nr1s,spar%nr2s,spar%nr3s !whole real space (rs) mesh (cf. system.h)
      ! print *,"GRIDINFO NRx  :",parm%nr1,parm%nr2,parm%nr3    !node-local real space mesh
      ! print *,"GRIDINFO KRxS :",kr1s,kr2s,kr3s !whole fft mesh
      ! print *,"GRIDINFO KRx  :",kr1,kr2,kr3    !node-local fft mesh
      !ENDIF
      !print *,"GRIDINFO NRXPL:",MEPOS,NRXPL(MEPOS,1),NRXPL(MEPOS,2)
      !note: fft mesh is always odd: its equal to rs mesh if the latter is odd, else its rs mesh + 1 (cf. LEADIM in loadpa.F)

!.......set voxel size
        GPOINTS(1) = spar%nr1s
        GPOINTS(2) = spar%nr2s
        GPOINTS(3) = spar%nr3s
        LTARGET=0.35 * 2.0D0 * FBOHR !target: half side length of 0.35 A

        DO I=1,3  !three dimensions
          DELTALMIN=9999999.0D0

          IF(paral%parent) WRITE(6,'(A,I1,A)',advance='no') "VOXINFO: POSSIBLE NO. OF VOXELS IN DIM. ",I, ": "
          DO J=1,GPOINTS(I)
            K = MOD(GPOINTS(I), J) !find integer divisors of grid 
            IF (K.EQ.0 .AND. J.GT.1 .AND. J.LT.GPOINTS(I)) THEN ! at least two grid points per voxel and at least two voxels per direction
              IF(paral%parent) WRITE(6,'(I5)',advance='no') J

              LVOXTMP=LBOX(I)/DBLE(J)
              DELTAL=ABS(LVOXTMP-LTARGET)
              IF (DELTAL.LT.DELTALMIN) THEN
                NVOX(I)=J
                LVOX(I)=LVOXTMP
                DELTALMIN=DELTAL
                IF(paral%parent) WRITE(6,'(A,F5.2,A)',advance='no') " (",LVOX(I)/FBOHR,") "
              ENDIF

            ENDIF
          ENDDO
          IF(paral%parent) WRITE(6,'(A)',advance='yes') " ."
       ENDDO

      IF (runinfo%voxnum(1).GT.0.AND.runinfo%voxnum(2).GT.0.AND.runinfo%voxnum(3).GT.0) THEN
        NVOX(1) = runinfo%voxnum(1)
        NVOX(2) = runinfo%voxnum(2)
        NVOX(3) = runinfo%voxnum(3)
!       voxel side length
        LVOX(3)=LBOX(3)/DBLE(NVOX(3))
        LVOX(2)=LBOX(2)/DBLE(NVOX(2))
        LVOX(1)=LBOX(1)/DBLE(NVOX(1))
      ENDIF
      
       VOXRGYR=0.5D0 * DSQRT( (LVOX(1)**2+LVOX(2)**2+LVOX(3)**2)/3.0D0 ) !Eq. S36
       
!.....print some information
      IF (paral%parent) THEN
        print *,"VOXINFO: NUM GRIDPOINTS  :",spar%nr1s,spar%nr2s,spar%nr3s
        print *,"VOXINFO: GRIDPOINTS/VOXEL:",DBLE(spar%nr1s)/DBLE(NVOX(1)), &
                     DBLE(spar%nr2s)/DBLE(NVOX(2)),DBLE(spar%nr3s)/DBLE(NVOX(3))

        WRITE(6,*)  "VOXINFO: ND, NVOX, LVOX (ANGS), RGYR",parai%mepos, &
         NVOX(1),NVOX(2),NVOX(3),LVOX(1)/FBOHR,LVOX(2)/FBOHR,LVOX(3)/FBOHR,VOXRGYR/FBOHR
      ENDIF
!.....end of print some information


      IF (NVOX(1).EQ.0.OR.NVOX(2).EQ.0.OR.NVOX(3).EQ.0) THEN
        CALL STOPGM('INIT_VOXELS','runinfo%voxnum MUST BE GREATER ZERO',0,"")
      ENDIF

      
      
!     allocate temporary memory for voxel placement (will be freed later)
      NVOXTOT = NVOX(1)*NVOX(2)*NVOX(3)
      ALLOCATE(KOVOXG(3,NVOXTOT),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      ALLOCATE(MYVOX(NVOXTOT),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      ALLOCATE(MYVG2VO(NVOXTOT),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)   
           
      NVOXTOT = 0     
           
!     position voxels. in this subroutine, each node knows all voxel coordinates KOVOXG (G for global)
      DO I=1,NVOX(1)
        DO J=1,NVOX(2)
          DO K=1,NVOX(3)
            NVOXTOT = NVOXTOT + 1
            KOVOXG(1,NVOXTOT) = runinfo%boxdum(1)-DI*0.50D0 + LVOX(1) * (DBLE(I)-0.5)
            KOVOXG(2,NVOXTOT) = runinfo%boxdum(3)-DJ*0.50D0 + LVOX(2) * (DBLE(J)-0.5)
            KOVOXG(3,NVOXTOT) = runinfo%boxdum(5)-DK*0.50D0 + LVOX(3) * (DBLE(K)-0.5)
          ENDDO
        ENDDO
      ENDDO
      
      
      ALLOCATE(VOXPART(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr) !kr1 not kr1s!
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)


!     loop over my part of the grid and assign gridpoints to nearest voxel
      IF (paral%parent)  WRITE(6,'(A)') "ASSIGNING VOXELS TO GRIDPOINTS..."
!$omp parallel private(K,J,I,KO_X,KO_Y,KO_Z, &
!$omp   DMIN,V,VG,DX,DY,DZ,D,idx)
!$omp do 
!.....LOOP OVER GRIDPOINTS
      DO K=1,fpar%kr3s
       KO_Z=runinfo%boxdum(5)+REAL(k-1,kind=real_8)*dk
       DO J=1,fpar%kr2s
         KO_Y=runinfo%boxdum(3)+REAL(j-1,kind=real_8)*dj
         DO I=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2) + ikr1! IKR1 added: biswas
           KO_X=runinfo%boxdum(1)+REAL(i-1,kind=real_8)*di
     
           idx = (i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s
           
           DMIN   = 999999D0
           VG     = -1
           DO V=1,NVOXTOT
             DX=KO_X-KOVOXG(1,V)
             DY=KO_Y-KOVOXG(2,V)
             DZ=KO_Z-KOVOXG(3,V)
             D=DX*DX+DY*DY+DZ*DZ
             IF (D .LT. DMIN) THEN
               DMIN=D
               VG=V
             END IF
           ENDDO
           VOXPART(idx)=VG
          ENDDO
        ENDDO
      ENDDO
!$omp end do 
!$omp end parallel


!.....now build voxel list of this node.(this is seperated from the above loop due to avoid cumbersome threadsave programming)
!     note that one voxel can live on multiple nodes
!.....loop over my gridpoints
!.....LOOP OVER GRIDPOINTS
      DO K=1,fpar%kr3s
       DO J=1,fpar%kr2s
         DO I=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2) + ikr1! IKR1 added: biswas
         
           idx = (i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s
                              
!          get voxel this gridpoint belongs to
           VG=VOXPART(idx)
           TFOUND = .FALSE.
!          loop over my voxellist...
           DO VO=1,MYNVOX
!            ..and continue with next gridpoint if voxel is already in my voxellist
             IF (MYVOX(VO).EQ.VG) THEN
                TFOUND=.TRUE.
                EXIT
             ENDIF
           ENDDO
!          ..if voxel is not in my voxellist, add it
           IF (.NOT.TFOUND) THEN
             MYNVOX=MYNVOX + 1
!            save mapping for resorting later in this function
             MYVG2VO(VG)    = MYNVOX
             MYVOX(MYNVOX)  = VG 
           ENDIF
          ENDDO
        ENDDO
      ENDDO
!.....now each node knows its number of voxels MYNVOX

!.....now we know MYNVOX so we can assign memory
      ALLOCATE(KOVOX(3,MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
           
      ALLOCATE(VO2IQ(MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      
!.....copy coordinates of voxel from global voxel list to my voxel list
!$omp parallel private(VO,VG)
!$omp do 
      DO VO=1,MYNVOX
        VG=MYVOX(VO)
        KOVOX(1,VO) = KOVOXG(1,VG)
        KOVOX(2,VO) = KOVOXG(2,VG)
        KOVOX(3,VO) = KOVOXG(3,VG)
      ENDDO
!$omp end do 
!$omp end parallel

     
!     loop once more over grid and remap VOXPART to the node-internal voxel numbering
!$omp parallel private(K,J,I,VO,VG,idx)
!$omp do 
      DO K=1,fpar%kr3s
       DO J=1,fpar%kr2s
         DO I=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)  + ikr1 ! IKR1 added: biswas
           idx = (i-parap%nrxpl(parai%mepos,1)+1) + (j-1)*fpar%kr1 + (k-1)*fpar%kr1*fpar%kr2s
           
           VG=VOXPART(idx)
           VO=MYVG2VO(VG)
           VOXPART(idx)=VO
          ENDDO
        ENDDO
      ENDDO
!$omp end do 
!$omp end parallel

     


!      calc oleradius from 2*rgyr / eps
       runinfo%oleradius = 2.0D0 * VOXRGYR / runinfo%oleeps   
      !note: each node does this calculation, i.e. runinfo%oleradius does not have to be broadcasted

       IF (paral%parent) THEN
         print *,"VOXINFO: runinfo%oleeps ", runinfo%oleeps,                 &
                        " runinfo%olensig ", runinfo%olensig*DSQRT(2.0D0),   &
               "runinfo%oleradius (BOHR,A)", runinfo%oleradius,runinfo%oleradius/FBOHR
       ENDIF
       
      
!     deallocate temporary memory for voxel placement (will be freed later)
      DEALLOCATE(KOVOXG,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
      DEALLOCATE(MYVOX,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
      DEALLOCATE(MYVG2VO,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)


      END SUBROUTINE
!     ==================================================================


!     ==================================================================
      SUBROUTINE INHERIT_LOCAL_EXP_WRAP(VO,DX,DY,DZ, LI)
!     ==--------------------------------------------------------------==
      INTEGER VO, LI
      REAL*8  DX,DY,DZ

      CALL inherit_local_exp(VOXEXP,VO,dimlocexp,DX,DY,DZ, epot1_var%locexp, LI)

      END SUBROUTINE
!     ==================================================================


      
!     ==================================================================
      SUBROUTINE SHIFT_MULTIPOLES_WRAP(TAR,DX,DY,DZ,VO) 
!     ==--------------------------------------------------------------==
      INTEGER VO
      REAL*8 DX,DY,DZ
      REAL*8 TAR(dimmultexp)

      
      call shift_multipoles(TAR,DX,DY,DZ,VOMULTEXP,VO,dimmultexp)
     
      END SUBROUTINE
!     ==================================================================

END MODULE iffi_comm_utils
