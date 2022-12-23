MODULE tst2min_inp_utils
  USE cnst_dyn,                        ONLY: &
       atcvar, cv_min2, cv_min_tol, cv_tol_ext, cvpar, file_str_ab, iangcv, &
       iatdlmn, icv_rmsd_ab, imincheck, itol_type, lmeta, max_minchk, &
       max_nnvar, max_pippo, max_search, ncolvar, ncvmin, nrmsd_ab, &
       specindex, trmsd_ab, tycvar
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert,&
                                             mmdim,&
                                             nat_cpmd
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tst2min_inp

CONTAINS

  ! ==================================================================
  SUBROUTINE tst2min_inp(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Collective Variables Input for Minima Finding         ==
    ! ==--------------------------------------------------------------==
    ! 
    ! NCOLVAR   :  number of Collective Variables
    ! TYCVAR    :  type of Collective Variables
    ! ATCVAR    :  atoms which define the Collective Variables
    ! CVPAR     :  other parameters that are required for the C.V.
    ! 
    ! ==--------------------------------------------------------------==
    ! how C.V. definition is read
    ! 
    ! DEFINE VARIABLE
    ! ==   nfix                                                       ==
    ! ==  DIST      n1  n2                                            ==
    ! ==  STRETCH   n1  n2                                            ==
    ! ==  BEND      n1  n2  n3                                        ==
    ! ==  TORSION   n1  n2  n3  n4                                    ==
    ! ==  OUTP      n1  n2  n3  n4                                    ==
    ! ==  COORD     n1   K  RC                                        ==
    ! ==  DIFFER    n1  n2  n3                                        ==
    ! ==  COORSP
    ! ==  COOR_RF  {2SHELL}
    ! ==  BNSWT
    ! ==  TOT_COOR  {2SHELL}
    ! ==  PLNANG
    ! ==  DIFCOOR
    ! ==  ...  
    ! END DEFINITION
    ! ==--------------------------------------------------------------==
    ! NCVMIN    :  number of known possible minima
    ! CV_MIN2   :  set of CV values for the all the known possible    ==
    ! minuma that could be found
    ! dimension  NCOLVAR*NCVMIN
    ! CV_MIN_TOL:  tolerance for accepting the ions configuration as  ==
    ! a minimum configuration (%)
    ! 
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'tst2min_inp'
    CHARACTER(len=10), DIMENSION(17), PARAMETER :: styp = (/'  STRETCH',&
      '     BEND','  TORSION',' DISTANCE','     OUTP',' COORDIN.','DIFFEREN.',&
      '   COORSP','  COOR_RF','BN.SWITCH',' TOT C_RF','  RMSD_AB','   DISPL1',&
      '  COORDIS','   PLNANG','  HBONDCH','  DIFCOOR'/)

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    CHARACTER(len=120)                       :: line2
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, i1, ia, iaa, icv, idummy, &
                                                ie, iee, ierr, ii, imin, &
                                                iout, ityp, k, my_nat, nnvar, &
                                                numspec
    INTEGER, SAVE                            :: firstdlmn = 0
    LOGICAL                                  :: erread, status
    REAL(real_8)                             :: tol_typ(13)

    lmeta%lsadpnt  = .TRUE.

    ! Default tolerance
    tol_typ(1) =  0.01_real_8  ! strech
    tol_typ(2) =  1.00_real_8  ! bend
    tol_typ(3) =  1.00_real_8  ! torsion
    tol_typ(4) =  0.01_real_8  ! distance
    tol_typ(5) =  1.00_real_8  ! outp
    tol_typ(6) =  0.20_real_8  ! coord
    tol_typ(7) =  0.01_real_8  ! differ
    tol_typ(8) =  0.03_real_8  ! coord species dependent
    tol_typ(9) =  0.03_real_8  ! coord sp. dep. with rational function
    tol_typ(10)=  0.03_real_8  ! bond switch  with rat. func.
    tol_typ(11)=  0.03_real_8  ! sp. dep. total coord  with rat. func.
    tol_typ(12)=  1.00_real_8
    tol_typ(13)=  0.03_real_8

    ! Default
    ! CV_MIN_TOL = 0.1_real_8
    imincheck  = 10
    max_minchk = 1
    max_search = 1000
    ncolvar    = 0
    itol_type = 1

    ! ==--------------------------------------------------------------==
    ! get total number of atoms and set indexing.
    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
    ELSE
       my_nat=mmdim%natm
    ENDIF

    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=20,END=20,fmt='(A)') line

    IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'SADDLE').NE.0)&
         GOTO 30
    IF (INDEX(line,'DEF').NE.0.AND.&
         INDEX(line,'VARIABLE').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) nnvar

       ! Allocate Memory
       ALLOCATE(tycvar(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(tycvar)!,nnvar)
       ALLOCATE(atcvar(15,nnvar/15),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(atcvar)!,15*nnvar)
       ALLOCATE(cvpar(10,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvpar)!,4*nnvar)
       ALLOCATE(cv_min_tol(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(specindex(nnvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(specindex)!,nnvar)


       GOTO 12
       ! Read Type and Atoms for the Collective Variables to be set 
11     READ(iunit,err=20,END=20,fmt='(A)') line
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'DEF').NE.0)&
            GOTO 10
12     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'DEF').NE.0) THEN
          ! Problem the number of coll. variables NCOLVAR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               '! # OF GIVEN C.V. =',nnvaR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               '! # OF FOUND C.V.=',ncolvaR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')&
               '! END C.V. IS REACHED.'
          CALL stopgm('META_COLCAR_INP',&
               'ERROR WHILE READING FIX STRUCTURES',& 
               __LINE__,__FILE__)
          ! ==------------------- DISTANCE --------------------------------==
       ELSEIF (INDEX(line,'DIST').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 4
          i1=INDEX(line,'DIST')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(4)
          ! ==--------------- DIFFERENCE AMONG DISTANCES ------------------==
       ELSEIF (INDEX(line,'DIFFER').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 7
          i1=INDEX(line,'DIFFER')+6
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(7)
          ! ==------------------ STRETCH ----------------------------------==
       ELSEIF (INDEX(line,'STRETCH').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 1
          i1=INDEX(line,'STRETCH')+7
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(1)
          ! ==------------------  BEND ------------------------------------==
       ELSEIF (INDEX(line,'BEND').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 2
          i1=INDEX(line,'BEND')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(2)
          ! ==------------------ TORSION ANGLE  ----------------------------==
       ELSEIF (INDEX(line,'TORSION').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 3
          i1=INDEX(line,'TORSION')+7
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(3)
          ! ==------------------ OUT OF LPANE ANGLE ------------------------==
       ELSEIF (INDEX(line,'OUTP').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 5
          i1=INDEX(line,'OUTP')+4
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          cv_min_tol(ncolvar) = tol_typ(5)
          ! ==------------------ COORDINATION NUMBER ------------------------==
       ELSEIF (INDEX(line,'COORD').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 6
          i1=INDEX(line,'COORD')+5
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
          cv_min_tol(ncolvar) = tol_typ(6)
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER  -----------==
          ! ==                 f = 1/(1+exp(k*(R-R_0)))                      ==
       ELSEIF (INDEX(line,'COORSP').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 8
          i1=INDEX(line,'COORSP')+6
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER -------------==
          ! ==            f =sum_i (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'COOR_RF').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 9
          i1=INDEX(line,'COOR_RF')+7
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') + 5
          ELSE
             specindex(ncolvar) = 1
          ENDIF
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF

             IF (paral%io_parent)&
                  READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
                  ii=1,atcvar(2,ncolvar))

          ENDIF
          ii=INDEX(line,'2SHELL')
          IF (ii.NE.0) THEN
             !>>> vw does it read an int or dble????
             !call readsi(line,i1,iout,cvpar(2,ncolvar),erread)
             CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
             !<<<
             ia=atcvar(3,ncolvar)
             ie=atcvar(4,ncolvar)
             IF (MOD(ia,2).NE.0 .OR. MOD(ie,2).NE.0 ) THEN
                CALL stopgm('TST2MIN_INP',&
                     'for COOR_RF with 2nd shell use even exponents',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          ! ==--------------- BOND SWITCH WITH RATIONAL F -------------------==
          ! ==            f = (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'BNSWT').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 10
          i1=INDEX(line,'BNSWT')+5
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          ! ==--------- TOTAL SPECIES DEPENDENT COORDINATION NUMBER ---------==
          ! ==            f =SUM_j[sum_i (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))]    ==
       ELSEIF (INDEX(line,'TOT_COOR').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 11
          i1=INDEX(line,'TOT_COOR') +8
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') + 5
          ELSE
             specindex(ncolvar) = 1
          ENDIF
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF

             IF (paral%io_parent)&
                  READ(iunit,*)  (iatdlmn(ii+my_nat*(ncolvar-1)),&
                  ii=1,atcvar(2,ncolvar))
          ENDIF
          ii=INDEX(line,'2SHELL')
          IF (ii.NE.0) THEN
             !>>> vw does it read an int or dble????
             !call readsi(line,i1,iout,cvpar(2,ncolvar),erread)
             CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
             !<<<
             ia=atcvar(3,ncolvar)
             ie=atcvar(4,ncolvar)
             IF (MOD(ia,2).NE.0 .OR. MOD(ie,2).NE.0 ) THEN
                CALL stopgm('M_COLVAR_INP',&
                     'for TOT_COOR with 2nd shell use even exponents',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          ! ==------- RMSD wrt 2 configurations A and B (read from file)-----==
          ! ==                f = (RMSD_A-RMSD_B)/(RMSD_A+RMSD_B)            ==
       ELSEIF (INDEX(line,'RMSD_AB').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 12
          IF (.NOT. trmsd_ab) THEN
             IF (max_nnvar.LT.nnvar) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) "NNVAR:",nnvar,&
                     "  MAX_NNVAR:",max_nnvar
                CALL stopgm("FILE_STR_AB",&
                     "PARAMETER MAX_NNVAR too small",& 
                     __LINE__,__FILE__)
             ENDIF

             ALLOCATE(icv_rmsd_ab(nnvar),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(icv_rmsd_ab)!,nnvar)
             trmsd_ab=.TRUE.
             nrmsd_ab = 0
          ENDIF
          icv_rmsd_ab(ncolvar)=ncolvar
          nrmsd_ab  = nrmsd_ab + 1

          i1=INDEX(line,'RMSD_AB') +7
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          IF (atcvar(1,ncolvar) .GT. 13) CALL stopgm('M_COLVAR_INP',&
               'RMSD_AB! too many species (>13)',& 
               __LINE__,__FILE__)
          numspec = atcvar(1,ncolvar)
          DO i = 1,numspec
             CALL readsi(line,i1,iout,atcvar(i+1,ncolvar),erread)
             i1=iout
          ENDDO
          ii=INDEX(line,'FILEAB')
          IF (ii.NE.0) THEN
             ia = ii+6
             line2  = line(ia:ia+20)
             CALL xstring(line2,ia,ie)
             file_str_ab(ncolvar) = line2(ia:ie)
          ELSE
             file_str_ab(ncolvar) = 'STRUCTURE_AB'
          ENDIF
          atcvar(numspec+2,ncolvar) = nrmsd_ab
          ! ==---- Angle between 2 planes, each defined by 3 given pnts  -----==
       ELSEIF (INDEX(line,'PLNANG').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 15
          ! IANGCV(NCOLVAR) = 2
          i1=INDEX(line,'PLNANG')+6
          iangcv(ncolvar) = 1
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(1,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(2,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(3,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(4,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(5,ncolvar)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          atcvar(6,ncolvar)=NAT_cpmd(idummy)
          ! ==------------ DIFFERENCE OF COORDINATION NUMBERS -------------==
          ! ==  f =sum_i (1-(Rbi/R_0)^n)/(1-(Rbi/R_0)^(n+m))-
          ! ==     sum_i (1-(Rai/R_0)^n)/(1-(Rai/R_0)^(n+m))
       ELSEIF (INDEX(line,'DIFCOOR').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 17
          i1=INDEX(line,'DIFCOOR')+7
          IF (INDEX(line,'INDAT') .NE. 0) THEN
             specindex(ncolvar) = 0
             i1 = INDEX(line,'INDAT') + 5
          ELSE
             specindex(ncolvar) = 1
          ENDIF
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          IF (specindex(ncolvar) .EQ. 0) THEN
             IF (firstdlmn .EQ. 0) THEN
                ALLOCATE(iatdlmn(my_nat*nnvar),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(iatdlmn)!,my_nat*nnvar)
                firstdlmn = 1
             ENDIF
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line2
             i1=1
             DO ii = 1,atcvar(3,ncolvar)
                CALL readsi(line2,i1,iout,&
                     iatdlmn(ii+my_nat*(ncolvar-1)),erread)
                i1=iout
             ENDDO
          ENDIF
          ! ==------------ COORDINATION OF SECOND NEIGHBORS -------------==
          ! ==    f =sum_jik [(1+(Rij/R_0)^n)/(1+(Rij/R_0)^(n+m))*       ==
          ! ==                1+(Rjk/R_0)^n)/(1+(Rjk/R_0)^(n+m))]/NA/NB  == 

       ELSEIF (INDEX(line,'COOR_CHAIN').NE.0) THEN
          ncolvar = ncolvar + 1
          tycvar(ncolvar) = 18
          i1=INDEX(line,'COOR_CHAIN')+10
          CALL readsi(line,i1,iout,atcvar(1,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(2,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(3,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(4,ncolvar),erread)
          i1=iout
          CALL readsi(line,i1,iout,atcvar(5,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(1,ncolvar),erread)
          i1=iout
          CALL readsr(line,i1,iout,cvpar(2,ncolvar),erread)
       ENDIF

       IF (ncolvar.GT.nnvar) THEN
          ! Problem the number of coll.variables NCOLVAR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'TST2MIN_INP! # OF GIVEN COL. VAR= ',nnvaR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)')&
               'TST2MIN_INP! # OF FOUND COL. VAR=',ncolvaR
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')&
               'THE NUMBER OF GIVEN COL. VAR. IS TOO SMALL.'
          CALL stopgm('TST2MIN_INP',&
               'ERROR WHILE READING FIX STRUCTURES',& 
               __LINE__,__FILE__)
       ELSEIF (ncolvar.GT. max_pippo) THEN
          CALL stopgm('TST2MIN_INP','NCOLVAR.GT. MAX_PIPPO',& 
               __LINE__,__FILE__)
       ELSEIF (ncolvar.LT.nnvar) THEN
          GOTO 12
       ELSE
          GOTO 11
       ENDIF
    ELSEIF (INDEX(line,'KNOWN_MINIMA').NE.0) THEN
       itol_type = 1
       IF (INDEX(line,'EXTENDED').NE.0) itol_type =2

       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) ncvmin

       ! Allocate Memory
       ALLOCATE(cv_min2(max_pippo,ncvmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_min2)!,ncvmin*max_pippo)

       IF (ncolvar.EQ.0) CALL stopgm('TST2MIN_INP',&
            'SET of KNOWN MIN.: # of CV=0, CV to be defined before ',& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,*) 'ITOL_TYPE ',itol_type
       IF (itol_type .EQ. 1) THEN
          DO imin = 1,ncvmin
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line
             i1 = 1
             ! read minima
             DO icv = 1,ncolvar
                CALL readsr(line,i1,iout,cv_min2(icv,imin),erread)
                i1 = iout
             ENDDO

          ENDDO
       ELSEIF (itol_type .EQ. 2) THEN
          ALLOCATE(cv_tol_ext(max_pippo,ncvmin),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(cv_tol_ext)!,ncvmin*max_pippo)

          DO imin = 1,ncvmin
             IF (paral%io_parent)&
                  READ(iunit,err=20,END=20,fmt='(A)') line
             i1 = 1
             ! read minima
             DO icv = 1,ncolvar
                CALL readsr(line,i1,iout,cv_min2(icv,imin),erread)
                i1 = iout
             ENDDO
             ! read tolerance
             DO icv = 1,ncolvar
                CALL readsr(line,i1,iout,cv_tol_ext(icv,imin),erread)
                i1 = iout
             ENDDO
             ! write(6,*) imin, (CV_MIN2(ICV,IMIN),ICV = 1,NCOLVAR),
             ! &                       (CV_TOL_EXT(ICV,IMIN),ICV = 1,NCOLVAR)
          ENDDO

       ENDIF
       GOTO 10
    ELSEIF(INDEX(line,'SADDLE').NE.0 .AND.&
         INDEX(line,'TOLERANCE') .NE.0) THEN

       IF (ncolvar.EQ.0) CALL stopgm('TST2MIN_INP',&
            'SET of TOLERANCE.: # of CV=0, CV to be defined before ',& 
            __LINE__,__FILE__)

       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*)&
            (cv_min_tol(icv),icv = 1,ncolvar)
       ! write(6,*) 'tol ',(CV_MIN_TOL(ICV),ICV = 1,NCOLVAR)
       GOTO 10
    ELSEIF (INDEX(line,'STEPCHECK').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) line
       CALL readsi(line,1,ie,imincheck,erread)
       IF (erread) GOTO 20
       GOTO 10
    ELSEIF (INDEX(line,'MAXCHECK').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) line
       CALL readsi(line,1,ie,max_minchk,erread)
       IF (erread) GOTO 20
       GOTO 10
    ELSEIF (INDEX(line,'MAXSEARCH').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) line
       CALL readsi(line,1,ie,max_search,erread)
       IF (erread) GOTO 20
       GOTO 10
    ELSE
       GOTO 10
    ENDIF  ! file biginning

20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING METADYNAMICS VARIABLES '
    CALL stopgm('TST2MIN_INP',' ',& 
         __LINE__,__FILE__)   
30  CONTINUE
    ! ==--------------------------------------------------------------==
    ! TEST consistencies

    lmeta%lmdreinit   = .FALSE.
    IF (ncvmin .EQ. 0) CALL stopgm('TST2MIN_INP',' NO MINIMA DEFINED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Print initialization
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(1x,A,A)')&
         '***           FROM SADDLE TO MINIMA',&
         '               ***'
    IF (paral%io_parent)&
         WRITE(6,'(5x,A,I5)')   '- Number of known minima = ',ncvmin
    IF (paral%io_parent)&
         WRITE(6,'(5x,A,I5,A)') '- Check performed each ',&
         imincheck,' MD STEP'
    IF (paral%io_parent)&
         WRITE(6,'(5x,A,I8)') '- Max. # of checks per search ', max_minchk
    IF (paral%io_parent)&
         WRITE(6,'(5x,A,I8)') '- Number of searching runs ',max_search
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(5X,A,10X,A)') 'TYPE','PARAMETERS'
    DO i=1,ncolvar
       ityp=tycvar(i)
       IF (ityp .EQ. 6) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,1X,I5,1X,2f8.4)') styp(ityp),&
               atcvar(1,i),cvpar(1,i),cvpar(2,i)
       ELSEIF (ityp.EQ.8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,1X,2I5,10X,2f8.4)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),cvpar(1,i),cvpar(2,i)
       ELSEIF (ityp .GE.9 .AND. ityp.LE.11) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,1X,2I4,6X,2I4,f8.5)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               cvpar(1,i)
       ELSEIF (ityp.EQ.15) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,1X,6I5)') styp(ityp),&
               (atcvar(k,i),k=1,6)
       ELSEIF (ityp .EQ. 17) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,1X,3I4,2x,2I3,f7.3)') styp(ityp),&
               atcvar(1,i),atcvar(2,i),atcvar(3,i),atcvar(4,i),&
               atcvar(5,i),cvpar(1,i)
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,3X,4I6)') styp(ityp),atcvar(1,i),&
               atcvar(2,i),atcvar(3,i),atcvar(4,i)
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)
    lineform = '(1X,I5,'//chnum(iaa:iee)//'f16.8)'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '       KNOWN MINIMA     '
    DO imin = 1,ncvmin
       IF (paral%io_parent)&
            WRITE(6,lineform) imin, (cv_min2(icv,imin),icv = 1,ncolvar)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (itol_type .EQ. 2) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') '       TOLERANCE PER MINIMUM       '
       DO imin = 1,ncvmin
          IF (paral%io_parent)&
               WRITE(6,lineform) imin, (cv_tol_ext(icv,imin),icv = 1,ncolvar)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,status)
    RETURN
  END SUBROUTINE tst2min_inp
  ! ==================================================================

END MODULE tst2min_inp_utils
