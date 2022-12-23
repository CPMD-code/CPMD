MODULE rinit_utils
  USE broy,                            ONLY: broy1
  USE cell,                            ONLY: cell_com
  USE clas,                            ONLY: clas8,&
                                             tclas
  USE cnst,                            ONLY: au_kb
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE gvec,                            ONLY: gvec_com
  USE isos,                            ONLY: isos1,&
                                             isos2,&
                                             isos3
  USE kdpc,                            ONLY: tkdp
  USE kdpoints_utils,                  ONLY: kdpoints
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE metr,                            ONLY: metr_com
  USE nlps,                            ONLY: imagp
  USE parac,                           ONLY: paral
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: response1,&
                                             rrk,&
                                             wwk
  USE rggen_utils,                     ONLY: rggen
  USE rkpnt_utils,                     ONLY: rkpnt
  USE sphe,                            ONLY: ngw_gamma,&
                                             ngwkps,&
                                             tsphere
  USE store_types,                     ONLY: restart1,&
                                             rout1,&
                                             store1
  USE symm,                            ONLY: stag,&
                                             symmi,&
                                             symmr,&
                                             symmt
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             dual00,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rinit

CONTAINS


  ! ==================================================================
  SUBROUTINE rinit
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'rinit'
    CHARACTER(len=10), DIMENSION(32), PARAMETER :: pgrp = (/'    1 (c1)',&
      '   <1>(ci)','    2 (c2)','   m (c1h)','  2/m(c2h)','    3 (c3)',&
      '  <3>(c3i)','   32 (d3)','  3m (c3v)',' <3>m(d3d)','    4 (c4)',&
      '   <4>(s4)','  4/m(c4h)','   422(d4)','  4mm(c4v)','<4>2m(d2d)',&
      '4/mmm(d4h)','    6 (c6)','  <6>(c3h)','  6/m(c6h)','   622(d6)',&
      '  6mm(c6v)','<6>m2(d3h)','6/mmm(d6h)','   222(d2)','  mm2(c2v)',&
      '  mmm(d2h)','    23 (t)','   m3 (th)','   432 (o)',' <4>3m(td)',&
      '   m3m(oh)'/)
    CHARACTER(len=12), DIMENSION(7), PARAMETER :: ihgp = (/'   TRICLINIC',&
      '  MONOCLINIC','ORTHORHOMBIC','  TETRAGONAL','       CUBIC',&
      '    TRIGONAL','   HEXAGONAL'/)
    REAL(real_8), PARAMETER                  :: deltakin = 1.e-10_real_8 

    CHARACTER(len=14)                        :: cihgp
    INTEGER                                  :: i, i1, i2, ierr, ik, isub
    REAL(real_8)                             :: rcerf1, rcerf2, rcerf3, &
                                                scr1(3,3)

    CALL tiset(procedureN,isub)
    IF (tkpts%tkpnt) THEN
       nkpt%ngwk=2*ncpw%ngw
       nkpt%nhgk=2*ncpw%nhg
       ! ..FNL is complex
       imagp=2
       ! G-VECTORS
       CALL rggen
       ! INITIALIZE K-POINTS
       CALL rkpnt
    ELSEIF (tkdp) THEN
       nkpt%ngwk=ncpw%ngw
       nkpt%nhgk=ncpw%nhg
       imagp=1
       rout1%teband=.FALSE.
       restart1%rocc=restart1%rocc.AND.cntl%tdiag
       ! G-VECTORS
       CALL rggen
       ! Initialize k.p k-points
       CALL kdpoints(kpts_com%nk1*kpts_com%nk2*kpts_com%nk3)
    ELSEIF (response1%tkpert) THEN
       nkpt%ngwk=ncpw%ngw
       nkpt%nhgk=ncpw%nhg
       imagp=1
       rout1%teband=.FALSE.
       restart1%rocc=restart1%rocc.AND.cntl%tdiag
       ! G-VECTORS
       CALL rggen
       ALLOCATE(rrk(3,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wwk(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(3*nkpt%nkpts,rk,1,rrk,1)
       CALL dcopy(nkpt%nkpts,wk,1,wwk,1)
       CALL zeroing(rk)!,3*nkpt%nkpts)
       DO ik = 1,nkpt%nkpts
          wk(ik) = 1.0_real_8
       ENDDO
    ELSE
       nkpt%ngwk=ncpw%ngw
       nkpt%nhgk=ncpw%nhg
       ! For non parallel set to one in SYSIN.
       nkpt%nkpts=1
       nkpt%nkpnt=1
       nkpt%nblkp=1
       imagp=1
       IF(.NOT.ALLOCATED(wk)) THEN
          ALLOCATE(wk(1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       wk(1)=1._real_8
       rout1%teband=.FALSE.
       restart1%rocc=restart1%rocc.AND.cntl%tdiag
       ! G-VECTORS
       CALL rggen
    ENDIF
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1) THEN
       IF (isos2%rcerf.LT.1.e-4_real_8) THEN
          rcerf1=cell_com%celldm(1)*0.5_real_8
          rcerf2=cell_com%celldm(2)*rcerf1
          rcerf3=cell_com%celldm(3)*rcerf1
          isos2%rcerf=(rcerf1+rcerf2+rcerf3)/3._real_8
          isos2%rcerf=isos2%rcerf/7._real_8
          isos2%rcerf=MIN(isos2%rcerf,2._real_8)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE(6,'(/,A,A)')' ************************** SUPERCELL',&
            ' ***************************'
       IF (isos1%tclust) THEN
          IF (isos1%toned) THEN
             WRITE(6,'(A)') ' THIS IS A POLYMER CALCULATION'
          ELSEIF (isos1%ttwod) THEN
             WRITE(6,'(A)') ' THIS IS A SURFACE CALCULATION'
          ELSE
             WRITE(6,'(A)') ' THIS IS AN ISOLATED SYSTEM CALCULATION'
          ENDIF
          IF (isos3%ps_type.EQ.1) THEN
             WRITE(6,'(A,T54,A)') ' POISSON EQUATION SOLVER  : ',&
                  ' HOCKNEY'
             WRITE(6,'(A,T54,F12.3)') ' COULOMB SMOOTHING RADIUS : ',&
                  isos2%rcerf
          ELSEIF (isos3%ps_type.EQ.2) THEN
             WRITE(6,'(A,T46,A)') ' POISSON EQUATION SOLVER  : ',&
                  ' TUCKERMAN & MARTYNA'
             WRITE(6,'(A,T54,F12.3)')&
                  ' SHORT RANGE POTENTIAL LENGTH * BOX LENGTH ',isos2%alphal
          ELSEIF (isos3%ps_type.EQ.3) THEN
             WRITE(6,'(A,T54,A)') ' POISSON EQUATION SOLVER  : ',&
                  ' MORTENSEN'
          ENDIF
       ENDIF
       IF ((parm%ibrav.EQ.1))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','       SIMPLE CUBIC'
       IF ((parm%ibrav.EQ.2))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','FACE CENTERED CUBIC'
       IF ((parm%ibrav.EQ.3))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','BODY CENTERED CUBIC'
       IF ((parm%ibrav.EQ.4))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','          HEXAGONAL'
       IF ((parm%ibrav.EQ.5))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','       RHOMBOHEDRAL'
       IF ((parm%ibrav.EQ.6))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','         TETRAGONAL'
       IF ((parm%ibrav.EQ.7))&
            WRITE(6,'(A,T43,A23)')&
            ' SYMMETRY:','BODY CENTRED TETRAGONAL'
       IF ((parm%ibrav.EQ.8))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','       ORTHORHOMBIC'
       IF ((parm%ibrav.EQ.12))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','         MONOCLINIC'
       IF ((parm%ibrav.EQ.14))&
            WRITE(6,'(A,T47,A19)') ' SYMMETRY:','          TRICLINIC'
       IF (symmi%indpg.NE.0) THEN
          IF (symmt%tmsym) THEN
             WRITE(6,'(A,A,I3)') ' MOLECULAR POINT GROUP : ',stag,symmi%naxis
             IF ((symmi%indpg.LT.0))&
                  WRITE(6,'(A,A,5X,A,I3)')&
                  ' POINT GROUP : SUBGROUP OF ',stag,'GROUP ORDER=',symmi%nrot
             symmi%indpg=ABS(symmi%indpg)
             DO i=1,symmi%nrot
                CALL dgemm('N','N',3,3,3,1._real_8,symmr%xtable(1,1,i),3,metr_com%ht(1,1),3,&
                     0._real_8,scr1(1,1),3)
                CALL dgemm('N','N',3,3,3,1._real_8,metr_com%htm1(1,1),3,scr1(1,1),3,&
                     0._real_8,symmr%xtable(1,1,i),3)
             ENDDO
          ELSE
             IF (symmi%indpg.GT.0) THEN
                IF (symmi%ihg.EQ.0) THEN
                   WRITE(cihgp,'(14(" "))')
                ELSE
                   CALL xstring(ihgp(symmi%ihg),i1,i2)
                   WRITE(cihgp,'("[",A,"]")') ihgp(symmi%ihg)(i1:i2)
                ENDIF
                WRITE(6,'(" POINT GROUP : ",T42,A14,A10)')cihgp,pgrp(&
                     symmi%indpg)
             ELSE
                WRITE(6,'(A,A10,T51,A12,I3)')' POINT GROUP : SUBGROUP OF'&
                     ,pgrp(-symmi%indpg),'GROUP ORDER=',symmi%nrot
                symmi%indpg=ABS(symmi%indpg)
             ENDIF
          ENDIF
       ENDIF
       IF (cntl%tprcp .OR. cntl%tpres) THEN
          WRITE(6,'(A,T54,F12.5)')' REFERENCE LATTICE CONSTANT(a.u.):',&
               parm%alat
          WRITE(6,'(A,T28,F8.3,5F6.3)') ' REFERENCE CELL DIMENSION:',&
               (prcp_com%cellrf(i),i=1,6)
          WRITE(6,'(A,T54,F12.4)')' REFERENCE VOLUME '//&
               '(OMEGA0 IN BOHR^3):', prcp_com%omega0
          WRITE(6,'(A,T28,F8.3,5F6.3)') ' INITIAL CELL DIMENSION:',&
               (cell_com%celldm(i),i=1,6)
          WRITE(6,'(A,T54,F12.4)')' INITIAL VOLUME(OMEGA IN BOHR^3):',&
               parm%omega
       ELSE
          WRITE(6,'(A,29X,F12.5)')' LATTICE CONSTANT(a.u.):',parm%alat
          WRITE(6,'(A,T18,6F8.4)')' CELL DIMENSION: ',(cell_com%celldm(i),i=1,6)
          WRITE(6,'(A,19X,F20.5)')' VOLUME(OMEGA IN BOHR^3): ',parm%omega
       ENDIF
       IF (cntl%tprcp) THEN
          IF ((prcpl%tzflex))&
               WRITE(6,'(A)') ' Z-DIRECTION FLEXIBLE CELL '
          IF ((prcpl%tisot))&
               WRITE(6,'(A)') ' ISOTROPICALLY FLEXIBLE CELL '

          IF (prcpl%tisot.OR.prcpl%tzflex) THEN
             WRITE(6,'(A,T46,F20.4)') ' REFERENCE PRESSURE(KBAR):',prcp_com%druck*&
                  au_kB
          ELSE
             WRITE(6,'(A,T42,3F8.2)') ' REFERENCE STRESS TENSOR(KBAR):',&
                  (prcp_com%stens(1,i)*au_kb,i=1,3)
             WRITE(6,'(T42,3F8.2)')(prcp_com%stens(2,i)*au_kb,i=1,3)
             WRITE(6,'(T42,3F8.2)')(prcp_com%stens(3,i)*au_kb,i=1,3)
          ENDIF
       ENDIF
       WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A1(BOHR): ',parm%a1
       WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A2(BOHR): ',parm%a2
       WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A3(BOHR): ',parm%a3
       WRITE(6,'(A,3F11.4)') ' RECIP. LAT. VEC. B1(2Pi/BOHR): ',gvec_com%b1(1)/&
            parm%alat,gvec_com%b1(2)/parm%alat,gvec_com%b1(3)/parm%alat
       WRITE(6,'(A,3F11.4)') ' RECIP. LAT. VEC. B2(2Pi/BOHR): ',gvec_com%b2(1)/&
            parm%alat,gvec_com%b2(2)/parm%alat,gvec_com%b2(3)/parm%alat
       WRITE(6,'(A,3F11.4)') ' RECIP. LAT. VEC. B3(2Pi/BOHR): ',gvec_com%b3(1)/&
            parm%alat,gvec_com%b3(2)/parm%alat,gvec_com%b3(3)/parm%alat
       IF (tclas.AND.isos1%tclust) THEN
          WRITE(6,'(A,T18,6F8.4)')' CLASSICAL CELL: ',(clas8%cellcl(i),i=1,6)
          WRITE(6,'(A,1X,T46,F20.5)')' VOLUME(OMEGA IN BOHR^3): ',&
               clas8%clomega
          WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A1(BOHR): ',clas8%cla1
          WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A2(BOHR): ',clas8%cla2
          WRITE(6,'(A,6X,3F11.4)') ' LATTICE VECTOR A3(BOHR): ',clas8%cla3
       ENDIF
       IF (prcp_com%akin.GT.deltakin) THEN
          WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)')&
               ' KINETIC ENERGY SMOOTHING(Ry) : A=',prcp_com%akin,'  S=',prcp_com%skin,&
               '  GC=',prcp_com%eckin
       ENDIF
       WRITE(6,'(A,T27,3(8X,I5))') ' REAL SPACE MESH:',spar%nr1s,spar%nr2s,spar%nr3s
       WRITE(6,'(A,T54,F12.5)')' WAVEFUNCTION CUTOFF(RYDBERG):',cntr%ecut
       WRITE(6,'(A,T36,A6,F5.2,A1,T54,F12.5)')&
            ' DENSITY CUTOFF(RYDBERG):','(DUAL=',dual00%cdual,')',dual00%cdual*cntr%ecut
       IF (tkpts%tkpnt.AND.tsphere) THEN
          WRITE(6,'(A,T54,I12)')&
               ' NUMBER OF PLANE WAVES PER WAVEFUNCTION ',spar%ngwks-1
          WRITE(6,'(A,T54,I12)')&
               ' NUMBER OF PLANE WAVES AT GAMMA POINT ',ngw_gammA
       ELSE
          IF (tkpts%tkpnt) THEN
             WRITE(6,'(A,T54,I12)')&
                  ' NUMBER OF PLANE WAVES FOR WAVEFUNCTION CUTOFF: ',&
                  spar%ngwks-1
          ELSE
             WRITE(6,'(A,T54,I12)')&
                  ' NUMBER OF PLANE WAVES FOR WAVEFUNCTION CUTOFF: ',&
                  spar%ngwks
          ENDIF
       ENDIF
       WRITE(6,'(A,T54,I12)')&
            ' NUMBER OF PLANE WAVES FOR DENSITY CUTOFF: ',spar%nhgs
       IF ((broy1%tgbroy) )&
            WRITE(6,'(A,T54,I12)')&
            ' NUMBER OF PLANE WAVES FOR BROYDEN MIXING: ',broy1%ngbroy
       IF ((fint1%ttrot))&
            WRITE(6,'(A,T54,1PE12.5)') ' TROTTER FACTOR :',fint1%betap
       IF (tkpts%tkpnt) THEN
          WRITE(6,'(1X,A,T58,I8)')&
               'KPOINTS (IN CARTESIAN COORDINATES AS INPUT):',nkpt%nkpts
          ! In input and output k points are in cartesian coordinates.
          WRITE(6,'(A,A)') '  NKP       KX           KY',&
               '           KZ         WEIGHT       NGW'

          DO ik=1,nkpt%nkpts
             WRITE(6,'(1X,I4,3F13.5,F11.5,I10)') ik,(rk(i,ik),i=1,3),wk(&
                  ik),ngwkps(ik)
          ENDDO
          IF ((tkpts%tkscale))&
               WRITE(6,'(1X,A)')&
               '[K-POINTS IN INPUT FILE WERE IN RECIPROCAL COORDINATES]'
       ENDIF
       WRITE(6,'(1X,64("*"))')
       WRITE(6,*)
       ! ==------------------------------------------------------------==
       ! Check if the number of states is not bigger that the number of
       ! components for each set of k points.
       ! If it is, the orthogonalisation is impossible.
       IF (tkpts%tkpnt) THEN
          DO ik=1,nkpt%nkpts
             IF (ngwkps(ik).LT.crge%n) THEN
                WRITE(6,'(A,I4,3F13.5,/,A,I6,/,A,I6,/,A)')&
                     ' RINIT! For the K points',ik,(rk(i,ik),i=1,3),&
                     ' RINIT! the number of components NGWK',ngwkps(ik),&
                     ' RINIT! is not greater than the number of states',crge%n,&
                     ' RINIT! The orthogonalisation is impossible !'
                CALL stopgm('RINIT','TOO MANY STATES',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDDO
       ELSE
          IF (2*spar%ngws-1.LT.crge%n) THEN
             WRITE(6,'(A,I6,/,A,I6,/,A)')&
                  ' RINIT! The number of components (2*NGWS-1)',&
                  2*spar%ngws-1,&
                  ' RINIT! is not greater than the number of states',crge%n,&
                  ' RINIT! The orthogonalisation is impossible !'
             CALL stopgm('RINIT','TOO MANY STATES',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ==------------------------------------------------------------==
       ! Check for the option TKNOSWAP.
       IF (tkpts%tknoswap) THEN
          IF (.NOT.cntl%wfopt) THEN
             WRITE(6,'(A,A)')&
                  ' RINIT| WITH YOUR OPTIONS, THE WAVEFUNCTIONS ARE ',&
                  'IRRELEVANT'
             CALL stopgm('RINIT',&
                  'USED ONLY WITH WAVEFUNCTIONS OPTMIZATION',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cnti%nomore.NE.1) THEN
             WRITE(6,'(A,A)')&
                  ' RINIT| WITH YOUR OPTIONS, THE WAVEFUNCTIONS ARE ',&
                  'IRRELEVANT'
             CALL stopgm('RINIT',&
                  'THE MAXIMUM NUMBER OF STEP HAS TO BE EQUAL TO 1',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (store1%swf.OR.store1%spot.OR.store1%srho) THEN
             WRITE(6,'(A,A)')&
                  ' RINIT| WITH YOUR OPTIONS, THE WAVEFUNCTIONS ARE ',&
                  'IRRELEVANT'
             CALL stopgm('RINIT',&
                  'USE: STORE WAVEFUNCTIONS DENSITY POTENTIAL OFF',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cntl%tpres.OR.cntl%tprcp) THEN
             WRITE(6,'(A,A)')&
                  ' RINIT| WITH YOUR OPTIONS, THE WAVEFUNCTIONS ARE ',&
                  'IRRELEVANT'
             CALL stopgm('RINIT','NO CALCULATION OF STRESS TENSOR',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
    ENDIF                     ! IF(IO_PARENT)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinit
  ! ==================================================================
END MODULE rinit_utils
