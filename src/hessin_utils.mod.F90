MODULE hessin_utils
  USE coor,                            ONLY: tau0
  USE cotr,                            ONLY: cotc0,&
                                             hess,&
                                             lskcor,&
                                             sdpl
  USE empfor_utils,                    ONLY: empfor
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE lscal,                           ONLY: hesscr,&
                                             icore,&
                                             lmap_p,&
                                             lnomap,&
                                             map_p,&
                                             ncore,&
                                             nvar
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE symtrz_utils,                    ONLY: symmat
  USE system,                          ONLY: cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: unitmx

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hessin
  PUBLIC :: gethes
  PUBLIC :: phessci
  PUBLIC :: clcmap

CONTAINS

  ! ==================================================================
  SUBROUTINE hessin(tau0,readhe)
    ! ==--------------------------------------------------------------==
    ! ==  Initialization of the nuclear hessian                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    LOGICAL                                  :: readhe

    CHARACTER(*), PARAMETER                  :: procedureN = 'hessin'

    CHARACTER(LEN=100)                       :: fname
    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: i, i1, i2, ierr, isub, iunit, &
                                                j, lhscr, nat00
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror
    REAL(real_8), ALLOCATABLE                :: hscr(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    fname='HESSIAN'
    iunit=21
    IF (ifirst.EQ.0) THEN
       ! For compatibility with secder (HESS(NODIM,NODIM)).
       ALLOCATE(hess(3*ions1%nat,3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    lhscr=2*3*ions1%nat*3*ions1%nat+10
    ALLOCATE(hscr(lhscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == CHECK IF FILE EXISTS                                         ==
    ! ==--------------------------------------------------------------==
    IF (readhe) THEN
       IF (paral%io_parent) THEN
          CALL fileopen(iunit,fname,fo_old,ferror)
          IF (ferror) THEN
             CALL xstring(fname,i1,i2)
             WRITE(6,*) 'HESSIN| HESSIAN FILE NOT FOUND:', fname(i1:i2)
             WRITE(6,*) 'HESSIN| INITIALIZE HESSIAN MATRIX'
             readhe=.FALSE.
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  READ OLD NUCLEAR HESSIAN FROM FILE                          ==
    ! ==--------------------------------------------------------------==
    IF (readhe) THEN
       IF (paral%io_parent) THEN
          REWIND(iunit)
          READ(iunit,err=999,END=999,fmt=*) nat00
          IF (cotc0%nodim.NE.nat00) THEN
             WRITE(6,'(A,2I5)')&
                  ' HESSIN| NUMBER OF PARAMETERS NOT CONSISTENT ',&
                  nat00,cotc0%nodim
             GOTO 999
          ENDIF
          DO i=1,cotc0%nodim
             READ(iunit,err=999,END=999,fmt=*) (hess(i,j),j=1,cotc0%nodim)
          ENDDO
          CALL fileclose(iunit)
       ENDIF
       GOTO 2000
    ELSE
       GOTO 1000
    ENDIF
999 CONTINUE
    CALL xstring(fname,i1,i2)
    IF (paral%io_parent)&
         WRITE(fformat,'(A,I2,A)')&
         '(A,T',MAX(33,65-(i2-i1)),',A)'
    IF (paral%io_parent)&
         WRITE(6,fformat) ' HESSIN| BAD HESSIAN FILE:',fname
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' HESSIN| INITIALIZE HESSIAN MATRIX'
    readhe=.FALSE.
    IF (paral%io_parent)&
         CALL fileclose(iunit)
1000 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  INITIAL HESSIAN IS UNIT MATRIX                              ==
    ! ==--------------------------------------------------------------==
    IF (cnti%npara.EQ.-1) THEN
       CALL unitmx(hess,cotc0%nodim)
    ELSEIF (cnti%npara.EQ.0 .OR. cnti%npara.EQ.1) THEN
       ! ==--------------------------------------------------------------==
       ! ==  EMPIRICAL INITIAL HESSIAN                                   ==
       ! ==--------------------------------------------------------------==
       CALL empfor(tau0,hscr)
       CALL gethes(hscr)
    ELSEIF (cnti%npara.GE.2) THEN
       ! ==--------------------------------------------------------------==
       ! ==  INITIAL HESSIAN FROM TRANSITION STATE SEARCH PARTIAL HESS.  ==
       ! ==--------------------------------------------------------------==
       CALL phessci(hess,cotc0%nodim)
    ENDIF
2000 CONTINUE
    ! WE CANNOT SYMMETRIZE THE MATRIX IF THE DEGREES NUMBER IS 
    ! SMALLER THAN 3 NAT.
    IF (cotc0%nodim.GE.3*ions1%nat) THEN
       CALL symmat(hess,0)
    ENDIF
    ! We use inverse Hessian matrix.
    IF (sdpl%tinvbfgs) THEN
       DO i=1,cotc0%nodim
          IF (hess(i,i).NE.0._real_8) THEN
             hess(i,i)=1._real_8/hess(i,i)
          ENDIF
       ENDDO
    ENDIF
    DEALLOCATE(hscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! we need to bcast the hessian if cp_groups > 1
    IF (parai%cp_nogrp.GT.1) THEN
       CALL mp_bcast(hess,cotc0%nodim**2,parai%cp_inter_io_source,parai%cp_inter_grp)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hessin
  ! ==================================================================
  SUBROUTINE gethes(rhess)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION OF THE NUCLEAR HESSIAN                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhess(3*ions1%nat,3*ions1%nat)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gethes'

    INTEGER                                  :: i1, i2, isub, j1, j2, k1, k2, &
                                                l1, l2, m1, m2

    CALL tiset(procedureN,isub)

    CALL unitmx(hess,cotc0%nodim)
    k1=0
    DO i1=1,ions1%nat
       DO j1=1,3
          l1=lskcor(j1,i1)
          IF (l1.NE.0) THEN
             k1=k1+1
             m1=(i1-1)*3+j1
             k2=0
             DO i2=1,ions1%nat
                DO j2=1,3
                   l2=lskcor(j2,i2)
                   IF (l2.NE.0) THEN
                      k2=k2+1
                      m2=(i2-1)*3+j2
                      hess(k1,k2) = rhess(m1,m2)
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

    RETURN
  END SUBROUTINE gethes
  ! ==================================================================

  SUBROUTINE phessci (hess, nodim)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nodim
    REAL(real_8)                             :: hess(nodim,nodim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'phessci'

    INTEGER                                  :: i, ierr, im, isub, j, jm
    REAL(real_8), ALLOCATABLE                :: hscr(:)

    CALL tiset(procedureN,isub)

    IF (.NOT.paral%parent) RETURN
    IF (cnti%npara.GE.3)  THEN
       ! CALL MEMORY (IP_HSCR, 2*3*NAT*3*NAT, 'HSCR')
       ALLOCATE(hscr(2*3*ions1%nat*3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Create map core <-> active if not yet there
    IF (.NOT.lnomap .AND. lmap_p.EQ.0) THEN
       ! CALL MEMORY (IP_MAP_P, NODIM, 'MAP_P')
       ALLOCATE(map_p(nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       lmap_p = nodim
       CALL clcmap (map_p, nodim, nvar, icore, ncore)

       ! CALL FREEM (IP_ICORE)
       DEALLOCATE(icore,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Init the full Hessian with unit matrix
    IF (cnti%npara.EQ.2) THEN
       CALL unitmx(hess,nodim)
       ! Init the full Hessian with empirical Hessian
    ELSE IF (cnti%npara.GE.3) THEN
       cnti%npara = cnti%npara - 3
       CALL empfor(tau0,hscr)
       CALL gethes(hscr)
       cnti%npara = cnti%npara + 3
    ENDIF
    ! Partial Hessian is usually stored as upper triangular matrix
    DO i = 1,nvar-1
       CALL dcopy (nvar-i, hesscr(i,i+1), nvar, hesscr(i+1,i), 1)
    ENDDO
    ! Check in partial Hessian
    IF (nvar.EQ.nodim) THEN
       CALL dcopy (nvar*nvar, hesscr, 1, hess, 1)
    ELSE IF (lnomap) THEN
       DO i = 1,nvar
          DO j = 1,nodim
             hess(i,j) = 0.0_real_8
             hess(j,i) = 0.0_real_8
          ENDDO
       ENDDO
       DO i = 1,nvar
          CALL dcopy (nvar, hesscr(1,i), 1, hess(1,i), 1)
       ENDDO
    ELSE
       DO i = 1,nvar
          im = map_p(i)
          DO j = 1,nodim
             hess(im,j) = 0.0_real_8
             hess(j,im) = 0.0_real_8
          ENDDO
       ENDDO
       DO i = 1,nvar
          im = map_p(i)
          DO j = 1,nvar
             jm = map_p(j)
             hess(jm,im) = hesscr(j,i)
          ENDDO
       ENDDO
    ENDIF
    IF (cnti%npara.GE.3) THEN
       ! CALL FREEM (IP_HSCR)
       DEALLOCATE(hscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

    RETURN
  END SUBROUTINE phessci


  SUBROUTINE clcmap (map_p, ndim, nvar, icore, ncore)
    ! ==--------------------------------------------------------------==
    ! ==  Calculate the mapping MAP_P between the degrees of freedom  ==
    ! ==  to be optimized (NODIM) and the ones of the reaction core   ==
    ! ==  (NVAR) and the environment (NODIM-NVAR) according to the    ==
    ! ==  atom list ICORE of the reaction core.                       ==
    ! ==                                                              ==
    ! ==  MAP_P(i_in_coord) = i_in_xpar (from RGMOPT)                 ==
    ! ==  i_in_coord      = 1.. NVAR:       index in COORD (core)    ==
    ! ==  i_in_coord-NVAR = 1.. NODIM-NVAR: index in VARS  (env.)    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim, map_p(ndim), nvar, &
                                                icore(nvar/3), ncore

    INTEGER                                  :: i, ia, iat, ic, ie, is, k, l
    LOGICAL                                  :: lenv

    iat = 0
    ic  = 0
    ie  = nvar
    l   = 0
    DO is = 1,ions1%nsp
       DO ia = 1,ions0%na(is)
          iat = iat + 1
          ! 
          ! ==  CHECK IF ATOM IAT IS CONTAINED IN ICORE(1..NCORE)           ==
          ! 
          i = 1
          lenv = .TRUE.
10        IF (lenv.AND.i.LE.ncore) THEN
             lenv = (iat.NE.icore(i))
             i = i + 1
             go to 10
          ENDIF
          ! 
          ! ==  REGISTER DEGREE OF FREEDOM EITHER FOR CORE OR ENVIRONMENT   ==
          ! 
          DO k = 1,3
             IF (lskcor(k,iat).NE.0) THEN
                l = l + 1! L is index in XPAR
                IF (lenv) THEN
                   ie = ie + 1
                   map_p(ie) = l
                ELSE
                   ic = ic + 1
                   map_p(ic) = l
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! 
    ! ==  FLAG INCONSISTENT INPUT                                     ==
    ! 
    IF (ic.NE.nvar .OR. ie.NE.cotc0%nodim) THEN
       IF (paral%io_parent)&
            WRITE (6,*) '  INCONSISTENT INPUT: ', ncore, ' ATOMS IN CORE'
       IF (paral%io_parent)&
            WRITE (6,*) '  NUMBER OF DEGREES OF FREEDOM SPECIFIED: ', nvar
       IF (paral%io_parent)&
            WRITE (6,*) '  NUMBER OF DEGREES OF FREEDOM FOUND:     ', ic
       CALL stopgm ('CLCMAP', 'INCONSISTENT INPUT',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE clcmap

END MODULE hessin_utils
