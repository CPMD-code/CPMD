MODULE wfn_print_utils
  USE cppt,                            ONLY: inyh
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: gloc_list,&
                                             gloci
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wfn_print
  !!public :: wfnz_print

CONTAINS

  ! ==================================================================
  SUBROUTINE wfn_print(infi,c0,nstate,psi)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: infi
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8) :: psi(-fpar%kr1s/2:fpar%kr1s/2,-fpar%kr2s/2:fpar%kr2s/2,&
      -fpar%kr3s/2:fpar%kr3s/2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wfn_print'

    CHARACTER(len=15)                        :: filbod, filen, filtr
    CHARACTER(len=8)                         :: numb
    INTEGER                                  :: i, i1, i2, ia, ie, ierr, ig, &
                                                istate, j, k, nh1, nh2, nh3, &
                                                no
    INTEGER, ALLOCATABLE, SAVE               :: orbpri(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror

! ==--------------------------------------------------------------==

    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1

    IF (ifirst.EQ.0) THEN
       no=nstate/2+1
       ALLOCATE(orbpri(no),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(orbpri)!,nstate)
       IF (gloci%gloc_all.NE.0) THEN
          DO i=1,nstate
             orbpri(i)=1
          ENDDO
       ELSEIF (gloci%gloc_first.GT.0) THEN
          DO i=gloci%gloc_first,gloci%gloc_last
             orbpri(i)=1
          ENDDO
       ELSEIF (gloci%gloc_orb.GT.0) THEN
          DO i=1,gloci%gloc_orb
             orbpri(gloc_list(i))=1
          ENDDO
       ENDIF
       ifirst=1
    ENDIF
    IF (paral%parent) THEN
       filtr='WFN_'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') infi
       CALL xstring(numb,ia,ie)
       filbod=filtr(1:4)//numb(ia:ie)
    ENDIF
    DO istate=1,nstate
       IF (orbpri(istate).NE.0) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(numb,'(I8)') istate
             CALL xstring(numb,ia,ie)
             CALL xstring(filbod,i1,i2)
             filen=filbod(i1:i2)//'.'//numb(ia:ie)
          ENDIF
          CALL zeroing(psi)!,kr1s*kr2s*kr3s)
          DO ig=1,ncpw%ngw
             i=inyh(1,ig)-nh1
             j=inyh(2,ig)-nh2
             k=inyh(3,ig)-nh3
             psi(i,j,k)=REAL(c0(ig,istate))*REAL(c0(ig,istate))+&
                  AIMAG(c0(ig,istate))*AIMAG(c0(ig,istate))
             psi(-i,-j,-k)=REAL(c0(ig,istate))*REAL(c0(ig,istate))+&
                  AIMAG(c0(ig,istate))*AIMAG(c0(ig,istate))
          ENDDO
          CALL mp_sum(psi,fpar%kr1s*fpar%kr2s*fpar%kr3s,parai%allgrp)

          IF (paral%io_parent) THEN
             CALL fileopen(12,filen,fo_def,ferror)
             WRITE(12,'(A,/)') filen
             WRITE(12,'(I4,3f10.4)') 1,0.0_real_8,0.0_real_8,0.0_real_8
             WRITE(12,'(I4,3f10.4)') spar%nr1s+1, 1.0_real_8,0.0_real_8,0.0_real_8
             WRITE(12,'(I4,3f10.4)') spar%nr2s+1, 0.0_real_8,1.0_real_8,0.0_real_8
             WRITE(12,'(I4,3f10.4)') spar%nr3s+1, 0.0_real_8,0.0_real_8,1.0_real_8

             WRITE(12,'(I4,4f10.4)') 8, 0.0_real_8,REAL(spar%nr1s,kind=real_8)/2._real_8,REAL(spar%nr2s,kind=real_8)/&
                  2._real_8,REAL(spar%nr3s,kind=real_8)/2._real_8
             DO i = -spar%nr1s/2,spar%nr1s/2
                DO j =-spar%nr2s/2,spar%nr2s/2
                   WRITE (12, '(6F15.7 )')  (psi(i,j,k), k = -spar%nr3s/2,spar%nr3s/&
                        2 )
                ENDDO
             ENDDO
             IF (paral%io_parent)&
                  CALL fileclose(12)
          ENDIF
       ENDIF
    ENDDO



    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wfn_print
  ! ==================================================================

  ! ==================================================================

END MODULE wfn_print_utils
SUBROUTINE wfnz_print(infi,ikind,c0,nstate,psi,nstate_loc)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:fpar,ncpw,spar
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:inyh
  USE g_loc , ONLY:glocal,gloci,gloc_list,ind_st,lostate
  USE fileopen_utils, ONLY : fileopen,fileclose
  USE fileopenmod , ONLY:fo_def
  USE readsr_utils, ONLY : readsr, xstring
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: infi, ikind, nstate
  COMPLEX(real_8)                            :: c0(2*ncpw%ngw,nstate)
  REAL(real_8) :: psi(-fpar%kr1s/2:fpar%kr1s/2,-fpar%kr2s/2:fpar%kr2s/2,-fpar%&
      kr3s/2:fpar%kr3s/2)
  INTEGER                                    :: nstate_loc

  CHARACTER(*), PARAMETER                    :: procedureN = 'wfnz_print'

  CHARACTER(len=15)                          :: filbod, filen, filtr
  CHARACTER(len=8)                           :: numb
  INTEGER                                    :: i, i1, i2, ia, icont, ie, &
                                                ierr, ig, istate, j, k, nh1, &
                                                nh2, nh3, no
  INTEGER, ALLOCATABLE, SAVE                 :: orbpri(:)
  INTEGER, SAVE                              :: ifirst = 0
  LOGICAL                                    :: ferror

! ==--------------------------------------------------------------==

  nh1=spar%nr1s/2+1
  nh2=spar%nr2s/2+1
  nh3=spar%nr3s/2+1

  IF (ifirst.EQ.0) THEN
     no=nstate/2+1
     ALLOCATE(orbpri(no),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(orbpri)!,nstate)
     IF (gloci%gloc_all.NE.0) THEN
        DO i=1,nstate
           orbpri(i)=1
        ENDDO
     ELSEIF (gloci%gloc_first.GT.0) THEN
        DO i=gloci%gloc_first,gloci%gloc_last
           orbpri(i)=1
        ENDDO
     ELSEIF (gloci%gloc_orb.GT.0) THEN
        DO i=1,gloci%gloc_orb
           orbpri(gloc_list(i))=1
        ENDDO
     ENDIF
     ifirst=1
  ENDIF
  IF (paral%parent) THEN
     filtr='WFN_'
     IF (paral%io_parent)&
          WRITE(numb,'(I8)') infi
     CALL xstring(numb,ia,ie)
     filbod=filtr(1:4)//numb(ia:ie)
     IF (glocal%tg_complex) THEN
        CALL xstring(filbod,i1,i2)
        filtr = filbod(i1:i2)
        IF (paral%io_parent)&
             WRITE(numb,'(I8)') ikind
        CALL xstring(numb,ia,ie)
        CALL xstring(filtr,i1,i2)
        filbod=filtr(i1:i2)//'_k'//numb(ia:ie)
     ENDIF
  ENDIF
  DO istate=1,nstate
     IF (orbpri(istate).NE.0) THEN
        IF (paral%parent) THEN
           IF (paral%io_parent)&
                WRITE(numb,'(I8)') istate
           CALL xstring(numb,ia,ie)
           CALL xstring(filbod,i1,i2)
           filen=filbod(i1:i2)//'.'//numb(ia:ie)
        ENDIF
        CALL zeroing(psi)!,kr1s*kr2s*kr3s)
        DO ig=1,ncpw%ngw
           i=inyh(1,ig)-nh1
           j=inyh(2,ig)-nh2
           k=inyh(3,ig)-nh3
           psi(i,j,k)=REAL(c0(ig,istate))*REAL(c0(ig,istate))+&
                AIMAG(c0(ig,istate))*AIMAG(c0(ig,istate))
           psi(-i,-j,-k)=REAL(c0(ig+ncpw%ngw,istate))*REAL(c0(ig+ncpw%ngw,&
                istate))+AIMAG(c0(ig+ncpw%ngw,istate))*AIMAG(c0(ig+ncpw%ngw,&
                istate))
        ENDDO
        CALL mp_sum(psi,fpar%kr1s*fpar%kr2s*fpar%kr3s,parai%allgrp)
        IF (paral%io_parent) THEN
           CALL fileopen(12,filen,fo_def,ferror)
           WRITE(12,'(A,/)') filen
           WRITE(12,'(I4,3f10.4)') 1,0.0_real_8,0.0_real_8,0.0_real_8
           WRITE(12,'(I4,3f10.4)') spar%nr1s+1, 1.0_real_8,0.0_real_8,0.0_real_8
           WRITE(12,'(I4,3f10.4)') spar%nr2s+1, 0.0_real_8,1.0_real_8,0.0_real_8
           WRITE(12,'(I4,3f10.4)') spar%nr3s+1, 0.0_real_8,0.0_real_8,1.0_real_8

           WRITE(12,'(I4,4f10.4)') 8, 0.0_real_8,REAL(spar%nr1s,kind=real_8)/2._real_8,REAL(spar%nr2s,kind=real_8)/&
                2._real_8,REAL(spar%nr3s,kind=real_8)/2._real_8

           DO i = -spar%nr1s/2,spar%nr1s/2
              DO j =-spar%nr2s/2,spar%nr2s/2
                 WRITE (12, '(6F15.7 )')  (psi(i,j,k), k = -spar%nr3s/2,spar%nr3s/&
                      2 )
              ENDDO
           ENDDO
           IF (paral%io_parent)&
                CALL fileclose(12)
        ENDIF
     ENDIF
  ENDDO


  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(numb,'(A3)') 'TOT'
     CALL xstring(numb,ia,ie)
     CALL xstring(filbod,i1,i2)
     filen=filbod(i1:i2)//'.'//numb(ia:ie)
  ENDIF
  CALL zeroing(psi)!,kr1s*kr2s*kr3s)
  DO icont = 1,nstate_loc
     IF (lostate%state_all) THEN
        istate = icont
     ELSE
        istate = ind_st(icont)
     ENDIF
     DO ig=1,ncpw%ngw
        i=inyh(1,ig)-nh1
        j=inyh(2,ig)-nh2
        k=inyh(3,ig)-nh3
        psi(i,j,k)=psi(i,j,k)+REAL(c0(ig,istate))*REAL(c0(ig,istate)&
             )+AIMAG(c0(ig,istate))*AIMAG(c0(ig,istate))
        psi(-i,-j,-k)=psi(-i,-j,-k)+REAL(c0(ig+ncpw%ngw,istate))*REAL(c0(&
             ig+ncpw%ngw,istate))+AIMAG(c0(ig+ncpw%ngw,istate))*AIMAG(c0(ig+ncpw%ngw,&
             istate))
     ENDDO
  ENDDO
  CALL mp_sum(psi,fpar%kr1s*fpar%kr2s*fpar%kr3s,parai%allgrp)

  IF (paral%io_parent) THEN
     CALL fileopen(12,filen,fo_def,ferror)
     WRITE(12,'(A,/)') fileN
     WRITE(12,'(I4,3f10.4)') 1,0.0_real_8,0.0_real_8,0.0_real_8
     WRITE(12,'(I4,3f10.4)') spar%nr1s+1, 1.0_real_8,0.0_real_8,0.0_real_8
     WRITE(12,'(I4,3f10.4)') spar%nr2s+1, 0.0_real_8,1.0_real_8,0.0_real_8
     WRITE(12,'(I4,3f10.4)') spar%nr3s+1, 0.0_real_8,0.0_real_8,1.0_real_8

     WRITE(12,'(I4,4f10.4)') 8, 0.0_real_8,REAL(spar%nr1s,kind=real_8)/2._real_8,REAL(spar%nr2s,kind=real_8)/&
          2._real_8,REAL(spar%nr3s,kind=real_8)/2._real_8
     DO i = -spar%nr1s/2,spar%nr1s/2
        DO j =-spar%nr2s/2,spar%nr2s/2
           WRITE (12, '(6F15.7 )')  (psi(i,j,k), k = -spar%nr3s/2,spar%nr3s/2 )
        ENDDO
     ENDDO
     IF (paral%io_parent)&
          CALL fileclose(12)
  ENDIF



  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wfnz_print
! ==================================================================
