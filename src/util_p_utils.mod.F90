#define warning_level 5.e-3_real_8

MODULE util_p_utils
  USE cppt,                            ONLY: gk
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm,&
                                             spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dot_ga

CONTAINS

  SUBROUTINE clean_rspace(array)
    REAL(real_8) :: array(fpar%kr1,fpar%kr2s,fpar%kr3s)

    INTEGER                                  :: i1, i2, i3

    i1 = fpar%kr1
    !$omp parallel do private(i3,i2)
    DO i3=1,spar%nr3s
       DO i2=1,spar%nr2s
          array(i1,i2,i3) = 0._real_8
       ENDDO
    ENDDO
    i2 = fpar%kr2s
    !$omp parallel do private(i3,i1)
    DO i3=1,spar%nr3s
       DO i1=1,parm%nr1
          array(i1,i2,i3) = 0._real_8
       ENDDO
    ENDDO
    i3 = fpar%kr3s
    !$omp parallel do private(i2,i1)
    DO i2=1,spar%nr2s
       DO i1=1,parm%nr1
          array(i1,i2,i3) = 0._real_8
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE clean_rspace
  ! ==================================================================
  FUNCTION dot_ga(ig,which_b) RESULT(res)
    ! computes and returns the dotproduct
    ! G . B_i
    ! where G is the reciprocal space vector defined by the index ig,
    ! and B_i is the reciprocal lattice vector of the system, 
    ! indexed by which_b. which_b=1, 2 or 3 corresponding to the system
    ! vectors B1, B2, B3.
    ! =----------------------------------------------------------------=
    INTEGER                                  :: ig, which_b
    REAL(real_8)                             :: res

    IF (ig .GT. ncpw%nhg)&
         CALL stopgm('DOT_GB','INTERNAL ERROR: wrong G vector index',& 
         __LINE__,__FILE__)

    IF (which_b .EQ. 1) THEN
       res = gk(1,ig)*parm%a1(1) + gk(2,ig)*parm%a1(2) + gk(3,ig)*parm%a1(3)
    ELSEIF (which_b .EQ. 2) THEN
       res = gk(1,ig)*parm%a2(1) + gk(2,ig)*parm%a2(2) + gk(3,ig)*parm%a2(3)
    ELSEIF (which_b .EQ. 3) THEN
       res = gk(1,ig)*parm%a3(1) + gk(2,ig)*parm%a3(2) + gk(3,ig)*parm%a3(3)
    ELSE
       CALL stopgm('DOT_GB','INTERNAL ERROR: wrong direction index',& 
            __LINE__,__FILE__)
    ENDIF
    res = parm%tpiba * res
    RETURN
  END FUNCTION dot_ga
  ! ==================================================================
  ! ==================================================================


END MODULE util_p_utils


! 
! WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
! TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
! THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
! 
! LINE   FORMAT      CONTENTS
! ===============================================================
! 1     A           TITLE
! 2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
! 3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
! 4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
! #ATOMS LINES OF ATOM COORDINATES:
! ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
! REST: 6E13.5      CUBE DATA
! 
! ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
! 
! ==================================================================
SUBROUTINE cubefile(filename,array,center,psi,thalfmesh)
  ! ==--------------------------------------------------------------==
  USE mp_interface, ONLY: mp_sync, mp_bcast, mp_send, mp_recv
  USE readsr_utils, ONLY : xstring
  USE system , ONLY:fpar,parm,spar,parap
  USE parac, ONLY : paral,parai
  USE ions , ONLY:ions0,ions1
  USE coor , ONLY:tau0
  USE metr , ONLY:metr_com
  USE response_pmod , ONLY:nris,real_8
  USE fileopen_utils, ONLY : fileopen
  USE fileopen_utils, ONLY : fileclose
  USE fileopenmod , ONLY:fo_def
  IMPLICIT NONE
  CHARACTER(len=*)                           :: filename
  REAL(real_8) :: array(fpar%kr1,fpar%kr2s,fpar%kr3s), center(3), &
      psi(fpar%kr2s,fpar%kr3s)
  LOGICAL                                    :: thalfmesh

  INTEGER                                    :: i1, i2, i3, ia, iat, ie, ii1, &
                                                ii2, ilowerleft(3), ir, isp, &
                                                msgid, source_pe, step
  INTEGER, SAVE                              :: fileunit = 99
  LOGICAL                                    :: ferror
  REAL(real_8)                               :: lowerleft(3), s(3)

  CALL mp_bcast(center,SIZE(center),parai%source,parai%allgrp)
  step = 1
  IF (thalfmesh) step=2
  nris(1)=spar%nr1s
  nris(2)=spar%nr2s
  nris(3)=spar%nr3s

  DO ir=1,3
     lowerleft(ir) = center(ir)&
          - 0.5_real_8*(parm%a1(ir)+parm%a2(ir)+parm%a3(ir))
  ENDDO
  CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,lowerleft,1,0._real_8,s,1)
  DO i1=1,3
     ilowerleft(i1) = NINT(s(i1)*nris(i1))+1
     ! generates 1..NR1S + n*NR1S
     lowerleft(i1) = 0._real_8
     DO i2=1,3
        lowerleft(i1)  = lowerleft(i1) +&
             metr_com%ht(i1,i2)*(ilowerleft(i2)-1)/REAL(nris(i2),kind=real_8)
     ENDDO
     ilowerleft(i1) = MOD((ilowerleft(i1)-1)+&
          1000*nris(i1),nris(i1))+1
  ENDDO

  IF (paral%parent) THEN
     ia=1
     CALL xstring(filename,ia,ie)
     IF (paral%io_parent)&
          CALL fileopen(fileunit,filename(ia:ie),fo_def,ferror)
     IF (ferror) GOTO 100

     ! make the file identifiable by its name in the title string.
     ! the second line tells other programs, that this is a
     ! 'density-style' cubefile (in case they actually check it).
     IF (paral%io_parent)&
          WRITE (fileunit,'(A)') ' CPMD CUBE FILE: '//filename(ia:ie)
     IF (paral%io_parent)&
          WRITE (fileunit,'(A)') ' Total SCF Density'
     IF (paral%io_parent)&
          WRITE (fileunit,200) ions1%nat, (lowerleft(ir),ir=1,3)
     ! note that we print the first mesh point in each direction
     ! twice to have smooth transitions for periodic visualizations.
     IF (paral%io_parent)&
          WRITE (fileunit,200) (spar%nr1s/step)+1,&
          (REAL(step,kind=real_8)*parm%a1(ir)/REAL(spar%nr1s,kind=real_8),ir=1,3)
     IF (paral%io_parent)&
          WRITE (fileunit,200) (spar%nr2s/step)+1,&
          (REAL(step,kind=real_8)*parm%a2(ir)/REAL(spar%nr2s,kind=real_8),ir=1,3)
     IF (paral%io_parent)&
          WRITE (fileunit,200) (spar%nr3s/step)+1,&
          (REAL(step,kind=real_8)*parm%a3(ir)/REAL(spar%nr3s,kind=real_8),ir=1,3)

     ! write atom coordinates
     DO isp=1,ions1%nsp
        DO iat=1,ions0%na(isp)
           IF (paral%io_parent)&
                WRITE (fileunit,201) ions0%iatyp(isp),REAL(ions0%iatyp(isp),kind=real_8),&
                (tau0(ir,iat,isp),ir=1,3)
        ENDDO
     ENDDO
  ENDIF

  ! NOTE: the '+step' is needed to print the first meshpoint 
  ! again at the end of the loop.
  DO i1=ilowerleft(1),spar%nr1s+ilowerleft(1)-1+step,step
     ii1=MOD(i1-1,spar%nr1s)+1
     ! Now: Copy array(ii1,*,*) into psi. NB: This is one YZ-plane,
     ! so it is NOT consecutively stored in memory...
     CALL mp_sync(parai%allgrp)! to avoid long transfer queues
     ! find the PE who has the plane i1:
     source_pe = 0
     DO WHILE (.NOT.&
          (ii1.GE.parap%nrxpl(source_pe,1) .AND.&
          ii1.LE.parap%nrxpl(source_pe,2)))
        source_pe = source_pe + 1
     ENDDO
     source_pe = parap%pgroup(source_pe+1)

     IF (parai%mepos .EQ. source_pe) THEN
        CALL dcopy(fpar%kr2s*fpar%kr3s,&
             array(ii1-parap%nrxpl(parai%mepos,1)+1,1,1),fpar%kr1,&
             psi,1)
        IF (.NOT. paral%parent) THEN
           !msglen = 8*kr2s*kr3s! one X-plane.
           msgid = 2! MPI message tag.
           CALL mp_send(psi,fpar%kr2s*fpar%kr3s,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
     ELSEIF (paral%parent) THEN
        !msglen = 8*kr2s*kr3s! one X-plane.
        msgid = 2     ! MPI message tag.
        CALL mp_recv(psi,fpar%kr2s*fpar%kr3s,source_pe,msgid,parai%allgrp)
     ENDIF

     ! Now, the parent should have the plane in PSI: write to cubefile.
     ! again, add step to get one more meshpoint in y- and z-direction 
     IF (paral%parent) THEN
        DO i2=ilowerleft(2),spar%nr2s+ilowerleft(2)-1+step,step
           ii2=MOD(i2-1,spar%nr2s)+1
           IF (paral%io_parent)&
                WRITE (fileunit,'(6E13.5)') (psi(ii2,MOD(i3-1,spar%nr3s)+1),&
                i3=ilowerleft(3),spar%nr3s+ilowerleft(3)-1+step,step)
        ENDDO
     ENDIF
  ENDDO                     ! i1

  IF ((paral%parent).AND.paral%io_parent)&
       CALL fileclose(fileunit)
  RETURN

  IF (paral%io_parent)&
       WRITE(6,*)'ERROR FINALLY CLOSING FILE.'
  RETURN

200 FORMAT (i5,3f12.6)
201 FORMAT (i5,4f12.6)

100 CONTINUE
  IF (paral%io_parent)&
       WRITE(6,*)'ERROR OPENING FILE:',filename
  IF (paral%io_parent)&
       CALL fileclose(fileunit)
  IF (paral%io_parent)&
       WRITE(6,*)'FILE SUCCESSFULLY RE-CLOSED:',filename
  RETURN
END SUBROUTINE cubefile
! ==================================================================

FUNCTION eigr_dot(isa,func) RESULT(my_res)
  ! ==--------------------------------------------------------------==
  ! INPUT:   One function f(G), packed storage, dim=NHG
  ! \        Species/Atom index ISA.
  ! OUTPUT:  function value returns
  ! \        Sum_G  exp(i G R_ion[isa])  f(G)
  ! \        i.e. the Fourier transform of f(G) at the position of
  ! \        the nucleus ISA.
  ! \        ATTENTION: The G=0 component of func is set to zero!!
  ! NOT PARALLEL -- corresponds to a DOTP-call
  ! ==--------------------------------------------------------------==
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,ncpw
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:inyh
  USE sfac , ONLY:ei1,ei2,ei3,eigrb,real_8
  USE geq0mod , ONLY:geq0
  IMPLICIT NONE
  INTEGER                                    :: isa
  REAL(real_8)                               :: func(2,*), my_res

  COMPLEX(real_8)                            :: expigr
  INTEGER                                    :: ig, ig1, isub
  REAL(real_8)                               :: dot
  REAL(real_8), EXTERNAL                     :: dotp

  CALL tiset('EIGR_DOT       ',isub)
  ig1=1
  IF (geq0) THEN
     ig1=2
     func(1,1) = 0._real_8
     func(2,1) = 0._real_8
  ENDIF
  IF (cntl%bigmem) THEN
     my_res = dotp(ncpw%nhg,eigrb(1,isa),func)
  ELSE
     my_res = 0._real_8
     !$omp parallel do private(ig,expigr,dot) reduction(+:my_res)
     DO ig=ig1,ncpw%nhg
        expigr = ei1(isa,inyh(1,ig))&
             * ei2(isa,inyh(2,ig))&
             * ei3(isa,inyh(3,ig))
        dot    = REAL(expigr)*func(1,ig)&
             + AIMAG(expigr)*func(2,ig)
        my_res = my_res + dot
     ENDDO
     my_res = 2._real_8 * my_res
  ENDIF
  CALL tihalt('EIGR_DOT       ',isub)
END FUNCTION eigr_dot
! ==================================================================

FUNCTION eigr_dot2(isa,func) RESULT(my_res)
  ! ==--------------------------------------------------------------==
  ! DIFFERENCE TO PREVIOUS ROUTINE (eigr_dot):
  ! THIS ONE DOES NOT SET THE IG=0 COMPONENT TO ZERO!!
  ! INPUT:   One function f(G), packed storage, dim=NHG
  ! \        Species/Atom index ISA.
  ! OUTPUT:  function value returns
  ! \        Sum_G  exp(i G R_ion[isa])  f(G)
  ! \        i.e. the Fourier transform of f(G) at the position of
  ! \        the nucleus ISA.
  ! \        ATTENTION: The G=0 component of func is set to zero!!
  ! NOT PARALLEL -- corresponds to a DOTP-call
  ! ==--------------------------------------------------------------==
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,ncpw
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:inyh
  USE sfac , ONLY:ei1,ei2,ei3,eigrb,real_8
  USE geq0mod , ONLY:geq0
  IMPLICIT NONE
  INTEGER                                    :: isa
  REAL(real_8)                               :: func(2,*), my_res

  COMPLEX(real_8)                            :: expigr
  INTEGER                                    :: ig, ig1, isub
  REAL(real_8)                               :: dot
  REAL(real_8), EXTERNAL                     :: dotp

  CALL tiset('EIGR_DOT2      ',isub)
  ig1=1
  my_res=0.0_real_8
  IF (geq0) THEN
     ig1=2
     expigr = ei1(isa,inyh(1,1))&
          * ei2(isa,inyh(2,1))&
          * ei3(isa,inyh(3,1))
     dot    = REAL(expigr)*func(1,1)&
          + AIMAG(expigr)*func(2,1)
     my_res = 0.5_real_8 * dot
  ENDIF
  IF (cntl%bigmem) THEN
     my_res = dotp(ncpw%nhg,eigrb(1,isa),func)
  ELSE
     !$omp parallel do private(ig,expigr,dot) reduction(+:my_res)
     DO ig=ig1,ncpw%nhg
        expigr = ei1(isa,inyh(1,ig))&
             * ei2(isa,inyh(2,ig))&
             * ei3(isa,inyh(3,ig))
        dot    = REAL(expigr)*func(1,ig)&
             + AIMAG(expigr)*func(2,ig)
        my_res = my_res + dot
     ENDDO
     my_res = 2._real_8 * my_res
  ENDIF
  CALL tihalt('EIGR_DOT2       ',isub)
END FUNCTION eigr_dot2

! ==================================================================
