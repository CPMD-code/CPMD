MODULE td_os_berry_utils
  USE cnst,                            ONLY: ry
  USE ddip,                            ONLY: lenbk,&
                                             ngwmax
  USE ddipo_utils,                     ONLY: setdip
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum,&
                                             mp_sync
  USE opeigr_p_utils,                  ONLY: opeigr_p
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_os_berry

CONTAINS

  ! ======================================================================
  SUBROUTINE td_os_berry(c0,c1,eigenvalue,nstate)
    ! ======================================================================
    REAL(real_8)                             :: eigenvalue
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'td_os_berry'

    COMPLEX(real_8)                          :: cone, czero, dd
    COMPLEX(real_8), ALLOCATABLE             :: cw(:), ddmat_1(:,:), &
                                                ddmat_p(:,:), f_raman(:,:,:), &
                                                sc0(:), work(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cwork(:,:,:)
    INTEGER                                  :: i, ierr, ig, il_cw, il_cw1, &
                                                il_framan, info, is1, isub, &
                                                k, n1, nmcol, nop1, nop2, &
                                                nop3, nxx
    INTEGER, ALLOCATABLE                     :: ipiv(:)
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: a, b, c, d, f3, f4, osx, osy, &
                                                osz, ox, spinfac

    CALL tiset('    td_os_raman_p',isub)

    cone=CMPLX(1._real_8,0._real_8,kind=real_8)
    czero=CMPLX(0._real_8,0._real_8,kind=real_8)
    ! ----------------------------------------------------------------------
    ! Memory allocations
    ! ----------------------------------------------------------------------
    il_framan=6*ncpw%ngw*nstate+12

    ! TODO check if +8 was needed for alignment only;
    ! have to remove now because it leads to segfault due to different allocation shape
    il_cw=2*ncpw%ngw*nstate        ! +8

    il_cw1=4*ncpw%ngw*nstate

    ALLOCATE(f_raman(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(f_raman)!,3*ngw*nstate)

    ALLOCATE(ddmat_p(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddmat_1(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipiv(50*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(nstate*50),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    lenbk=nxxfun(nstate)
    nxx=MAX(2*lenbk*parai%nproc,il_cw)
    il_cw=nxx
    ALLOCATE(sc0(nxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cw(il_cw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cw)!,SIZE(cw))
    IF (ifirst.EQ.0) THEN
       n1=0
       DO i=0,parai%nproc-1
          n1=MAX(n1,parap%sparm(3,i))
       ENDDO
       ngwmax=n1
       ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nmcol=(parai%nproc*ngwmax)
       ALLOCATE(mapcol(nmcol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cwork(ncpw%ngw,nstate,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setdip(mapful,mapcol)
       ifirst=1
    ENDIF
    ! ----------------------------------------------------------------------
    ! End of memory allocations
    ! ----------------------------------------------------------------------

    ! .....Spin degeneracy factor
    IF (cntl%tlsd) THEN
       spinfac=1._real_8
    ELSE
       spinfac=2._real_8
    ENDIF

    ! ----------------------------------------------------------------------
    ! Set up phase operator matrices (k=1,3) and f_raman.
    ! ----------------------------------------------------------------------
    osx = 0._real_8
    osy = 0._real_8
    osz = 0._real_8
    DO k=1,3

       CALL zeroing(cw)!,SIZE(cw))
       CALL zeroing(ddmat_p)!,nstate*nstate)
       CALL zeroing(ddmat_1)!,nstate*nstate)
       CALL zeroing(cwork)!,2*ngw*nstate)
       CALL zeroing(sc0)!,SIZE(sc0))
       CALL zeroing(f_raman)!,SIZE(f_raman))

       ! .......Compute cwork = opeigr*c(g) and Q(i,j) = <phi_i | opeigr | phi_j>.
       IF (parm%ibrav.EQ.1 .OR. parm%ibrav.EQ.8) THEN
          nop1=k
          nop2=0
          nop3=0
       ELSE IF (parm%ibrav.EQ.2) THEN
          IF (k.EQ.1) THEN
             nop1=1
             nop2=2
             nop3=0
          ELSE IF (k.EQ.2) THEN
             nop1=1
             nop2=3
             nop3=0
          ELSE
             nop1=2
             nop2=3
             nop3=0
          ENDIF
       ELSE IF (parm%ibrav.EQ.3) THEN
          IF (k.EQ.1) THEN
             nop1=1
             nop2=-2
             nop3=-3
          ELSE IF (k.EQ.2) THEN
             nop1=1
             nop2=2
             nop3=-3
          ELSE
             nop1=1
             nop2=2
             nop3=3
          ENDIF
       ELSE
          IF (paral%io_parent)&
               WRITE (6,*)&
               'TD_OS_BERRY: ONLY SUPERCELLS OF TYPE 1,2,3,8 allowed.'
       ENDIF
       CALL mp_sync(parai%allgrp)
       CALL opeigr_p(c0,cw,sc0,nstate,mapful,mapcol,ddmat_p,&
            nop1,nop2,nop3,dd,cwork)
       CALL mp_sync(parai%allgrp)

       ! ........Make a copy of ddmat_p (ddmat_p need be inverted).
       CALL zcopy(nstate*nstate,ddmat_p,1,ddmat_1,1)

       ! ----------------------------------------------------------------------
       ! Compute inverse of ddmat_p (Q^-1) for mu=k in Putrino et al.
       ! Note ddmat_p is overwrittwen.
       ! ----------------------------------------------------------------------
       CALL zgetrf(nstate,nstate,ddmat_p,nstate,ipiv,info)
       IF (info.EQ.0) THEN
          CALL zgetri(nstate,ddmat_p,nstate,ipiv,work,nstate,info)
       ELSE
          CALL stopgm('td_os_berry','error in matrix inversion',& 
               __LINE__,__FILE__)
       ENDIF
       ! ----------------------------------------------------------------------
       ! End.
       ! ----------------------------------------------------------------------


       ! ----------------------------------------------------------------------
       ! Build f_raman.
       ! ----------------------------------------------------------------------
       ! ........cwork(1,*,*) = [opeigr*c(g)] * Q^-1;
       ! ........cwork(2,*,*) = [opeigr*c(-g)] * Q^-1.
       DO i=1,2
          CALL zeroing(cw)!,SIZE(cw))
          CALL zgemm('N','N',ncpw%ngw,nstate,nstate,cone,cwork(1,1,i),ncpw%ngw,&
               ddmat_p(1,1),nstate,czero,cw,ncpw%ngw)
          CALL zcopy(ncpw%ngw*nstate,cw,1,cwork(1,1,i),1)
       ENDDO
       CALL mp_sync(parai%allgrp)

       IF (geq0) THEN
          is1=2
          DO i=1,nstate
             b=AIMAG(cwork(1,i,1))
             f_raman(1,i,k)=CMPLX(b,0._real_8,kind=real_8)
          ENDDO
       ELSE
          is1=1
       ENDIF

       ! This is to have a sum running on ngw below.
       DO i=1,nstate
          DO ig=is1,ncpw%ngw
             a=REAL(cwork(ig,i,1))
             b=AIMAG(cwork(ig,i,1))
             c=REAL(cwork(ig,i,2))
             d=AIMAG(cwork(ig,i,2))
             f3=a-c
             f4=b+d
             f_raman(ig,i,k)=0.5_real_8*CMPLX(f4,-f3,kind=real_8)
          ENDDO
       ENDDO
       ! ----------------------------------------------------------------------
       ! End.
       ! ----------------------------------------------------------------------


       ! ----------------------------------------------------------------------
       ! Compute oscillator strengths.
       ! ----------------------------------------------------------------------
       ! .......Compute os = sum_k sum_i <chi_i| opeigr(k) | phi_i> / Q + c.c.
       IF (geq0) THEN
          is1=2
          DO i=1,nstate
             ox = REAL( CONJG(c1(1,i))*f_raman(1,i,k) ) +&
                  REAL( CONJG(f_raman(1,i,k))*c1(1,i) )

             IF (k.EQ.1) osx = osx + ox
             IF (k.EQ.2) osy = osy + ox
             IF (k.EQ.3) osz = osz + ox
          ENDDO
       ELSE
          is1=1
       ENDIF

       DO i=1,nstate
          DO ig=is1,ncpw%ngw
             ox = 2._real_8*REAL( CONJG(c1(ig,i))*f_raman(ig,i,k)  ) +&
                  2._real_8*REAL( CONJG(f_raman(ig,i,k))*c1(ig,i)  )

             IF (k.EQ.1) osx = osx + ox
             IF (k.EQ.2) osy = osy + ox
             IF (k.EQ.3) osz = osz + ox
          ENDDO
       ENDDO
       ! ----------------------------------------------------------------------
       ! End.
       ! ----------------------------------------------------------------------
    ENDDO                                        ! On cartesian components
    CALL mp_sync(parai%allgrp)
    CALL mp_sum(osx,parai%allgrp)
    CALL mp_sum(osy,parai%allgrp)
    CALL mp_sum(osz,parai%allgrp)

    IF (paral%io_parent)&
         WRITE(6,'(A,T25,A,F10.5," eV",T54,A,F10.5)')&
         ' TD_OS_BERRY|',' dE=',&
         eigenvalue*RY*2._real_8,'f=',&
         2._real_8 *&
         (2._real_8/3._real_8)* eigenvalue * ( (spinfac*osx)**2&
         +   (spinfac*osy)**2&
         +   (spinfac*osz)**2 )
    IF (paral%io_parent)&
         WRITE(6,'(A,T26,"x:",F10.5,"  y:",F10.5,"  z:",F10.5)')&
         ' TD_OS_BERRY|',2._real_8*osx, 2._real_8*osy, 2._real_8*osz

    IF (paral%parent) THEN
       WRITE(6,'(A,T26,"x:",F10.5,"  y:",F10.5,"  z: ", F10.5)')&
            ' TD_OS_BERRY|',2._real_8*osx, 2._real_8*osy, 2._real_8*osz
    ENDIF

    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(f_raman,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ipiv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddmat_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddmat_1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    td_os_raman_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE td_os_berry


END MODULE td_os_berry_utils
