MODULE tst2min_utils
  USE cnst_dyn,                        ONLY: &
       cv_ist, cv_min2, cv_min_tol, cv_tol_ext, imincheck, itol_type, lmeta, &
       max_minchk, max_search, ncolvar, ncvmin
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE kinds,                           ONLY: real_8
  USE meta_colvar_inp_utils,           ONLY: colvarofr
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE soft,                            ONLY: soft_com

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tst2min

CONTAINS

  ! ==================================================================
  SUBROUTINE tst2min(taup,tscr,infi)
    ! ==--------------------------------------------------------------==
    ! ==  Checks whether the a minimum is reached
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taup(:,:,:), tscr(:,:,:)
    INTEGER                                  :: infi

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform1, lineform2
    CHARACTER(len=20)                        :: file1 = 'sadsrch.out', &
                                                file2 = 'cv_trj.out'
    INTEGER                                  :: iaa, icv, iee, imin, min_found
    INTEGER, SAVE                            :: i_checks = 0, i_search = 1, &
                                                if1 = fo_verb, if2 = fo_verb
    LOGICAL                                  :: ferror, lminreached
    REAL(real_8)                             :: diff_dv

! ==--------------------------------------------------------------==

    lmeta%lmdreinit   = .FALSE.
    lminreached = .FALSE.


    ! DO IMIN = 1,NCVMIN
    ! write(6,'(I5,3f16.8)') imin, (CV_MIN2(ICV,IMIN),ICV = 1,NCOLVAR)
    ! ENDDO
    ! STOP  'tst2'
    ! DO IMIN = 1,NCVMIN
    ! write(6,*) CV_MIN(1,IMIN),CV_MIN(2,IMIN)
    ! ENDDO

    ! stop 'cv_min'
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)
    lineform1 = '(1X,4I8,'//chnum(iaa:iee)//'f16.8)'
    lineform2 = '(1X,2I6,'//chnum(iaa:iee)//'f16.8)'

    IF (.NOT.paral%io_parent) GOTO 9999

    IF (MOD(infi,imincheck) .EQ. 0 ) THEN
       i_checks = i_checks+1


       ! ==--------------------------------------------------------------==
       ! Calculate the Collective Variables and their Derivatives
       CALL colvarofr(taup,tscr)

       IF (paral%io_parent)&
            CALL fileopen(52,file2,fo_app+if2,ferror)
       ! Write output
       IF (paral%io_parent)&
            WRITE(52,lineform2) i_search,i_checks,&
            (cv_ist(icv),icv=1,ncolvar)
       IF (paral%io_parent)&
            CALL fileclose(52)
       if2=0

       DO imin = 1,ncvmin

          DO icv = 1,ncolvar
             IF (paral%io_parent)&
                  WRITE(6,*)  'CV_MIN ',imin,cv_min2(icv,imin)

             diff_dv = ABS(cv_ist(icv)-cv_min2(icv,imin))

             IF (itol_type .EQ. 1) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(/A,2I5,2f16.8/)')&
                     'chk ', imin,icv,diff_dv,cv_min_tol(icv)
                IF (diff_dv .GT. cv_min_tol(icv)) GOTO 100
             ELSEIF (itol_type .EQ. 2) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(/A,2I5,2f16.8/)')&
                     'chk ', imin,icv,diff_dv,cv_tol_ext(icv,imin)

                IF (diff_dv .GT. cv_tol_ext(icv,imin)) GOTO 100
             ENDIF
          ENDDO
          lminreached  = .TRUE.
          min_found  = imin

          GOTO 200

100       CONTINUE
       ENDDO
200    CONTINUE

       IF (lminreached .OR. i_checks .EQ. max_minchk) THEN
          lmeta%lmdreinit = .TRUE.
          IF (paral%io_parent)&
               WRITE(6,'(/10x,A,I6,A/)') ' Search # ',i_search,' ENDED'
          IF (paral%io_parent)&
               WRITE(6,'(/10x,A,I6)') ' TOTAL # OF MD STEP =  ', infi


          ! Open output
          IF (paral%io_parent)&
               CALL fileopen(51,file1,fo_app+if1,ferror)

          IF (lminreached) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6/)')&
                  ' Minimum found corrisponding to # ',min_found
             IF (paral%io_parent)&
                  WRITE(51,lineform1) i_search, i_checks, infi,min_found,&
                  (cv_ist(icv),icv=1,ncolvar)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A/)') 'Minimum not found, max # checks reached'
             IF (paral%io_parent)&
                  WRITE(51,lineform1) i_search, i_checks, infi, 0 ,&
                  (cv_ist(icv),icv=1,ncolvar)
          ENDIF

          IF (paral%io_parent)&
               CALL fileclose(51)
          if2=0

          i_checks = 0

          IF (i_search .EQ. max_search) THEN
             soft_com%exsoft = .TRUE.
             lmeta%lmdreinit = .FALSE.
             IF (paral%io_parent)&
                  WRITE(6,'(/10x,A)') 'Maximun # of searchs reached'
             IF (paral%io_parent)&
                  WRITE(6,'(/10x,A)') 'Iteration stops'
          ENDIF
          i_search  = i_search + 1
       ENDIF  ! LMINREACHED

    ENDIF      ! IMINCHECK

9999 CONTINUE   ! PARENT

    ! IF(I_CHECKS .EQ. 5) STOP '5 checks'

    CALL mp_sync(parai%cp_grp)

    CALL mp_bcast(lmeta%lmdreinit,parai%io_source,parai%cp_grp)
    CALL mp_bcast(soft_com%exsoft,parai%io_source,parai%cp_grp)


    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE tst2min
  ! ==================================================================

END MODULE tst2min_utils
