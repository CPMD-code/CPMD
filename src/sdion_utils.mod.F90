MODULE sdion_utils
  USE cotr,                            ONLY: cotc0,&
                                             dmax,&
                                             dtm,&
                                             pxpar,&
                                             sdpl
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fixcom_utils,                    ONLY: fixcom
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE sort_utils,                      ONLY: sort2
  USE tpar,                            ONLY: dt2bym
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sdion
  PUBLIC :: minlin
  !public :: sort
  !public :: xlag2
  PUBLIC :: give_scr_sdion
  PUBLIC :: scalar

CONTAINS

  ! ==================================================================
  SUBROUTINE sdion(xpar,dxpar)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE IONIC POSITIONS WITH STEEPEST DESCENT OR             ==
    ! ==         WITH CONJUGATE GRADIENT                              ==
    ! ==  ADD MINIMIZATION ALONG THE LINE                             ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   XPAR(NODIM)                                                ==
    ! ==  DXPAR(NODIM)                                                ==
    ! == OUTPUT:                                                      ==
    ! ==   XPAR(NODIM)                                                ==
    ! ==  (NODIM Number of Freedom degres)                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim), &
                                                dxpar(cotc0%nodim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'sdion'
    INTEGER, PARAMETER                       :: nstepmax = 20

    INTEGER                                  :: i, ierr, info
    INTEGER, SAVE                            :: ifirst = 0, nstep
    REAL(real_8)                             :: beta, dr, f2, fac1, fpoldb, &
                                                hlen, p2, pf1
    REAL(real_8), ALLOCATABLE                :: aux(:)
    REAL(real_8), EXTERNAL                   :: ddot
    REAL(real_8), SAVE                       :: dd, ddold, de(nstepmax), &
                                                dstep, dx(nstepmax), &
                                                ee(nstepmax), f2old, fpold, &
                                                p2old

    ALLOCATE(aux(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (ifirst.EQ.0) THEN
       ! Determines the first step
       dstep=0._real_8
       DO i=1,ions1%nsp
          dstep=dstep+dt2bym(i)
       ENDDO
       dstep=dstep/REAL(ions1%nsp,kind=real_8)
       dd=dstep
       ALLOCATE(pxpar(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(pxpar)!,cotc0%nodim)
       IF (sdpl%tcgp) THEN
          fac1=0.91_real_8
          f2old = 1._real_8
          p2old = 0._real_8
       ELSE
          fac1=0._real_8
       ENDIF
       nstep=0
       ifirst=1
    ENDIF
    IF (sdpl%tsdpline.OR.sdpl%tcgp) THEN
       ! ==------------------------------------------------------------==
       ! == MINIMIZATION ALONG THE LINE BY STEEPEST DESCENT            ==
       ! ==------------------------------------------------------------==
       nstep=nstep+1
       ee(nstep)=ener_com%etot
       de(nstep)=scalar(cotc0%nodim,pxpar,dxpar)
       IF (nstep.EQ.1) THEN
          dx(nstep)=0._real_8
          IF (dstep.EQ.0._real_8) THEN
             dstep=dd
          ELSE
             dd=dstep
          ENDIF
       ELSE
          dx(nstep)=dstep
       ENDIF
       CALL minlin(nstep,dx,ee,de,dstep,info)
       IF (info.NE.0) THEN
          ! First step
          ddold=0._real_8
          dx(nstep)=0._real_8
          IF (sdpl%tcgp) THEN
             ! Conjugate gradient
             f2=scalar(cotc0%nodim,dxpar,dxpar)
             IF (info.EQ.2) THEN
                ! We reset the conjugate gradient direction
                DO i=1,cotc0%nodim
                   pxpar(i)=0._real_8
                ENDDO
                p2old=0._real_8
             ENDIF
             fpold=0._real_8
             DO i=1,cotc0%nodim
                fpold=dxpar(i)*pxpar(i)
             ENDDO
             IF (f2old.EQ.0._real_8) THEN
                beta=0._real_8
             ELSE
                beta = f2 / f2old * fac1
             ENDIF
             fpoldb = fpold * beta
             p2 = f2 + p2old*beta*beta + 2._real_8*fpoldb
             ! If cos( (P,F) ) < 0.1 , we reinitialize P = F .
             IF (f2.NE.0._real_8.AND.p2.NE.0._real_8) THEN
                IF ( (f2 + fpoldb)/SQRT(f2*p2).LT.0.1_real_8) THEN
                   beta=0._real_8
                   p2=f2
                ENDIF
             ENDIF
             f2old = f2
             p2old = p2
             ! First step
             DO i=1,cotc0%nodim
                pxpar(i)=dxpar(i)+beta*pxpar(i)
             ENDDO
             de(1)=scalar(cotc0%nodim,pxpar,dxpar)
          ELSE
             ! Steepest descent
             beta=0._real_8
             ! First step
             DO i=1,cotc0%nodim
                pxpar(i)=dxpar(i)
             ENDDO
             de(1)=scalar(cotc0%nodim,dxpar,dxpar)
          ENDIF
       ENDIF
       dr=dstep-ddold
       ! Compute the ionic steps
       DO i=1,cotc0%nodim
          aux(i)=-dr*pxpar(i)
       ENDDO
       ddold=dstep
    ELSE
       ! ==------------------------------------------------------------==
       ! == OLD VERSION OF SDION                                       ==
       ! ==------------------------------------------------------------==
       ! Compute the ionic steps
       DO i=1,cotc0%nodim
          aux(i)=-dtm(i)*dxpar(i)
       ENDDO
       ! Check if displacements are not too big
       hlen=SQRT(ddot(cotc0%nodim,aux(1),1,aux(1),1))
       IF (hlen.GT.dmax) THEN
          pf1=dmax/hlen
          CALL dscal(cotc0%nodim,pf1,aux(1),1)
       ENDIF
    ENDIF
    ! Fix the center of mass (if LFCOM=.TRUE.)
    CALL fixcom(aux)
    CALL daxpy(cotc0%nodim,1.0_real_8,aux(1),1,xpar(1),1)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sdion
  ! ==================================================================
  SUBROUTINE minlin(nstep,dx,ee,de,dstep,info)
    ! ==--------------------------------------------------------------==
    ! == MINIMISATION ALONG LINE                                      ==
    ! == The NSTEPth point is the last calculated one                 ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   NSTEP   Number of steps                                    ==
    ! ==   DX(NSTEP) Abscissae of points                              ==
    ! ==   EE(NSTEP) Energy of points                                 ==
    ! ==   DE(NSTEP  First derivative of points                       ==
    ! == OUTPUT:                                                      ==
    ! ==   DSTEP   Estimated abscissae for root                       ==
    ! ==   INFO    0 O.K. we continue                                 ==
    ! ==           1 We need to have initialization                   ==
    ! ==             (new direction)                                  ==
    ! ==   NSTEP   Number of points stored for the next call          ==
    ! ==   DX,EE, and EE are sorted                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstep
    REAL(real_8)                             :: dx(nstep), ee(nstep), &
                                                de(nstep), dstep
    INTEGER                                  :: info

    INTEGER, PARAMETER                       :: nstepmax = 20 
    REAL(real_8), PARAMETER                  :: golden = 1.618034_real_8 

    INTEGER                                  :: i, INDEX(nstepmax)
    REAL(real_8)                             :: dmax, dmin, dnew, &
                                                xnew(nstepmax)

    IF (nstep.LE.1) THEN
       info=1
       RETURN
    ENDIF
    ! Determine min et max
    dmin=dx(1)
    dmax=dx(1)
    DO i=1,nstep
       IF (dmin.GT.dx(i)) dmin=dx(i)
       IF (dmax.LT.dx(i)) dmax=dx(i)
    ENDDO
    ! We sort the points
    CALL sort2(ee,nstep,index)
    CALL sort(nstep,dx,xnew,index)
    CALL sort(nstep,de,xnew,index)
    CALL xlag2(nstep,ee,dx,de,dmin,dmax,dnew,info)
    IF (info.EQ.0) THEN
       ! Good
       dstep=dnew
       IF (nstep.EQ.2) THEN
          ! We check if the approximation is O.K.
          info=0
       ELSE
          ! Check is the last point has the lowest EE
          IF (INDEX(1).EQ.nstep) THEN
             ! We have finished line minimisation
             info=1
             nstep=1
             ! DSTEP a little bit bigger
             dstep=dstep*1.1_real_8
          ELSE
             ! Use the lowest
             dstep=dx(1)
             nstep=0
             info=0
          ENDIF
       ENDIF
       RETURN
    ELSEIF (info.EQ.-1) THEN
       ! DNEW smaller than DMIN
       ! Error
       info=2
       nstep=1
    ELSEIF (info.EQ.1) THEN
       ! DNEW bigger than DMAX
       ! Check if DNEW not too big
       IF (dnew.GT.dmax*golden*golden*golden*golden) THEN
          dstep=dmax*golden*golden*golden*golden
       ELSE
          dstep=dnew
       ENDIF
       info=0
       nstep=1
    ELSEIF (info.EQ.2) THEN
       ! No interpolation
       info=1
       nstep=1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE minlin
  ! ==================================================================
  SUBROUTINE sort(n,x,xnew,index)
    ! ==--------------------------------------------------------------==
    ! == SORT X ARRAY WITH INDEX information                          ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==  N Dimension of X(1:N), X(1:N) and INDEX(1:N)                ==
    ! == OUTPUT:   X(1:N)                                             ==
    ! ==  XNEW(1:N) is a scratch array                                ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), xnew(n)
    INTEGER                                  :: INDEX(n)

    INTEGER                                  :: i

    DO i = 1 , n
       xnew(i)=x(INDEX(i))
    ENDDO
    CALL dcopy(n,xnew,1,x,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sort
  ! ==================================================================
  SUBROUTINE xlag2(nstep,e,d,de,dmin,dmax,dstep,info)
    ! ==--------------------------------------------------------------==
    ! == GIVE ROOT OF THE SECOND POLYNOMIAL GIVEN BY TWO POINTS       ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   E(1:NSTEP)    Energy for each point                        ==
    ! ==   D(1:NSTEP)    Coordinate                                   ==
    ! ==   DE(1:NSTEP)   First derivative of energy along the line    ==
    ! ==   DMIN          Minimal abscissae (bounds)                   ==
    ! ==   DMAX          maximal abscissae                            ==
    ! == OUTPUT:                                                      ==
    ! ==   DSTEP the abscissae of the new point                       ==
    ! ==   INFO  0 DSTEP is between DMIN and DMAX                     ==
    ! ==        -1 DSTEP < DMIN                                       ==
    ! ==         1 DSTEP > DMAX                                       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstep
    REAL(real_8)                             :: e(nstep), d(nstep), &
                                                de(nstep), dmin, dmax, dstep
    INTEGER                                  :: info

    INTEGER                                  :: i, iother

    IF (nstep.LE.1) RETURN
    ! ==--------------------------------------------------------------==
    ! The first one has the lowest values
    iother=0
    DO i=2,nstep
       IF (de(1)*de(i).LT.0._real_8) THEN
          iother=i
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    IF (iother.NE.0) THEN
       ! Good , we have bounds
       info=0
       IF (de(iother)-de(1).EQ.0) THEN
          ! No interpolation (normally impossible)
          info=2
       ELSE
          dstep=d(1)-de(1)*(d(iother)-d(1))/(de(iother)-de(1))
       ENDIF
    ELSE
       ! Use the two lowest values
       IF (de(2)-de(1).EQ.0) THEN
          ! No interpolation (normally impossible)
          info=2
       ELSE
          dstep=d(1)-de(1)*(d(2)-d(1))/(de(2)-de(1))
       ENDIF
       ! Check if the root is inside the bounds
       IF (dstep.LT.dmin) THEN
          info=-1
       ELSEIF (dstep.GT.dmax) THEN
          info=1
       ELSE
          info=0
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xlag2
  ! ==================================================================
  FUNCTION scalar(n,a,b)
    ! ==--------------------------------------------------------------==
    ! == SCALAR PRODUCT                                               ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n), b(n), scalar

    INTEGER                                  :: i

    scalar=0._real_8
    DO i=1,n
       scalar=scalar+a(i)*b(i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION scalar
  ! ==================================================================
  SUBROUTINE give_scr_sdion(lsdion,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsdion
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       lsdion=cotc0%nodim
       tag ='NODIM'
    ELSE
       lsdion=0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_sdion
  ! ==================================================================

END MODULE sdion_utils
