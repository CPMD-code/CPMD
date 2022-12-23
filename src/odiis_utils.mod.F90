#include "cpmd_global.h"

MODULE odiis_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: maxdis
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !!public :: odiis
  PUBLIC :: updis
  PUBLIC :: solve
  !public :: grimax
  !public :: dicopy
  !public :: idcopy
  !public :: sdcopy
  !public :: dscopy
  !public :: daxpy
  !public :: dscal
  !public :: intpcoef

CONTAINS

  SUBROUTINE updis(ndiis,nowv,nsize,mdiis,iact)
    ! ==--------------------------------------------------------------==
    ! == Update variables for cntl%diis if IACT=1                          ==
    ! == Set to 0                  if IACT/=1                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiis, nowv, nsize, mdiis, &
                                                iact

! ==--------------------------------------------------------------==

    IF (iact.EQ.1) THEN
       nowv=MOD(ndiis,mdiis)+1
       ndiis=ndiis+1
       nsize=ndiis+1
       IF (nsize.GT.mdiis) nsize=mdiis+1
    ELSE
       ndiis=0
       nowv=0
       nsize=0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updis
  ! ==================================================================
  SUBROUTINE solve(b,ldb,ndim,v)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldb
    REAL(real_8)                             :: b(ldb,*)
    INTEGER                                  :: ndim
    REAL(real_8)                             :: v(*)

    INTEGER                                  :: info, lw, nr
    REAL(real_8)                             :: scr1(maxdis+1,maxdis+1), &
                                                scr2(maxdis+1,maxdis+1), &
                                                toleig

    IF (ndim.GT.maxdis+1) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'SOLVE! NDIM=',ndim,' MAXDIS+1=',maxdis+1
       CALL stopgm('SOLVE','MDIM GREATER THAN MAXDIS+1',& 
            __LINE__,__FILE__)
    ENDIF
    toleig=0.2e-15_real_8
    lw=(maxdis+1)*(maxdis+1)
    CALL dgelss(ndim,ndim,1,b,ldb,v,ldb,scr1,toleig,nr,scr2,lw,info)
    IF (info.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'SOLVE! INFO=',infO
       CALL stopgm('SOLVE','COULD NOT SOLVE DIIS EQUATION',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE solve

#if defined(__ES)
  ! mb   Please, do NOT remove the part here below: It is vital for ES !
  ! ==================================================================
  SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
    ! ==--------------------------------------------------------------==
    ! constant times a vector plus a vector.
    ! uses unrolled loops for increments equal to one.
    ! jack dongarra, linpack, 3/11/78.
    ! modified 12/3/93, array(1) declarations changed to array(*)
    ! mb   Updated 7/2/2004 for ES
    ! ==--------------------------------------------------------------==
    DOUBLE PRECISION dx(*),dy(*),da
    INTEGER :: i,incx,incy,ix,iy,m,mp1,n
    ! ==--------------------------------------------------------------==
    IF (n.LE.0) RETURN
    IF (da .EQ. 0.0_real_8) RETURN
    IF (incx.EQ.1.AND.incy.EQ.1) GOTO 20
    ! 
    ! code for unequal increments or equal increments
    ! not equal to 1
    ! 
    ix = 1
    iy = 1
    IF (incx.LT.0) ix = (-n+1)*incx + 1
    IF (incy.LT.0) iy = (-n+1)*incy + 1
    DO i = 1,n
       dy(iy) = dy(iy) + da*dx(ix)
       ix = ix + incx
       iy = iy + incy
    ENDDO
    RETURN
    ! 
    ! code for both increments equal to 1
20  CONTINUE
    !$omp parallel do private(i)
    DO i = 1,n
       dy(i) = dy(i) + da*dx(i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE daxpy
  ! ==================================================================
  SUBROUTINE dscal(n,da,dx,incx)
    ! ==--------------------------------------------------------------==
    ! scales a vector by a constant.
    ! uses unrolled loops for increment equal to one.
    ! jack dongarra, linpack, 3/11/78.
    ! modified 3/93 to return if incx .le. 0.
    ! modified 12/3/93, array(1) declarations changed to array(*)
    ! mb   Updated 7/2/2004 for ES
    ! ==--------------------------------------------------------------==
    DOUBLE PRECISION da,dx(*)
    INTEGER :: i,incx,m,mp1,n,nincx
    ! ==--------------------------------------------------------------==
    IF (incx.EQ.1) THEN
       ! 
       ! code for increment equal to 1
       !$omp parallel do private(i)
       DO i = 1,n
          dx(i) = da*dx(i)
       ENDDO
    ELSE
       ! 
       ! code for increment not equal to 1
       ! 
       nincx = n*incx
       DO i = 1,nincx,incx
          dx(i) = da*dx(i)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dscal
  ! ==================================================================
#endif



END MODULE odiis_utils




! ==================================================================
SUBROUTINE intpcoef(num,nsize,ldpme,pme,vc,c0)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: num, nsize, ldpme
  COMPLEX(real_4)                            :: pme(ldpme,*)
  REAL(real_8)                               :: vc(*)
  COMPLEX(real_8)                            :: c0(*)

  COMPLEX(real_8)                            :: ctmp
  INTEGER                                    :: i, j

  !$omp parallel do private(I,J,CTMP) schedule(static)
  DO i=1,num
     ctmp=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
     DO j=1,nsize-1
        ctmp=ctmp+vc(j)*CMPLX(pme(i,j),kind=real_8)
     ENDDO
     c0(i)=ctmp
  ENDDO
END SUBROUTINE intpcoef
! ==================================================================
SUBROUTINE idcopy(n,c,x,grmax,gimax)
  ! ==--------------------------------------------------------------==
  ! == UNPACK N complex(8) :: NUMBERS FROM TWO 1 BYTE INTEGERS         ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
#if defined (__NEC)
  ! Arguments
  INTEGER :: n,x(n/4+3)
  REAL(real_8) :: c(2,n),grmax,gimax
  ! Variables
  REAL(real_8) :: fr,fi
  INTEGER :: i,j,nrest
  ! ==--------------------------------------------------------------==
  fr=grmax/127._real_8
  fi=gimax/127._real_8
  !$omp parallel do private(I,J)
  DO i=1,n/4
     j=(i-1)*4
     c(1,j+1)=ISHFT(x(i),-56)
     c(2,j+1)=ISHFT(ISHFT(x(i), 8),-56)
     c(1,j+2)=ISHFT(ISHFT(x(i),16),-56)
     c(2,j+2)=ISHFT(ISHFT(x(i),24),-56)
     c(1,j+3)=ISHFT(ISHFT(x(i),32),-56)
     c(2,j+3)=ISHFT(ISHFT(x(i),40),-56)
     c(1,j+4)=ISHFT(ISHFT(x(i),48),-56)
     c(2,j+4)=ISHFT(ISHFT(x(i),56),-56)
  ENDDO
  nrest=MOD(n,4)
  IF (nrest.EQ.1) THEN
     c(1,n)=ISHFT(x(n/4+1),-56)
     c(2,n)=ISHFT(ISHFT(x(n/4+1), 8),-56)
  ELSEIF (nrest.EQ.2) THEN
     c(1,n-1)=ISHFT(x(n/4+1),-56)
     c(2,n-1)=ISHFT(ISHFT(x(n/4+1), 8),-56)
     c(1,n)=ISHFT(ISHFT(x(n/4+1),16),-56)
     c(2,n)=ISHFT(ISHFT(x(n/4+1),24),-56)
  ELSEIF (nrest.EQ.3) THEN
     c(1,n-2)=ISHFT(x(n/4+1),-56)
     c(2,n-2)=ISHFT(ISHFT(x(n/4+1), 8),-56)
     c(1,n-1)=ISHFT(ISHFT(x(n/4+1),16),-56)
     c(2,n-1)=ISHFT(ISHFT(x(n/4+1),24),-56)
     c(1,n)=ISHFT(ISHFT(x(n/4+1),32),-56)
     c(2,n)=ISHFT(ISHFT(x(n/4+1),40),-56)
  ENDIF
  !$omp parallel do private(I)
  DO i=1,n
     c(1,i)=(c(1,i)-128._real_8)*fr
     c(2,i)=(c(2,i)-128._real_8)*fi
  ENDDO
#else
  ! INTEGER = INTEGER*4
  INTEGER(int_4) :: n
  REAL(real_8) :: c(2,n), grmax,gimax, fr,fi
  INTEGER(int_4) :: x(*), i,j,nrest
  ! ==--------------------------------------------------------------==
  fr=grmax/127._real_8
  fi=gimax/127._real_8
#if defined(__SR11000) || defined(__PRIMEHPC)
  !$omp parallel do private(I,J)
#endif 
  DO i=1,n/2
     j=(i-1)*2
     c(1,j+1)=ISHFT(x(i),-24)
     c(2,j+1)=ISHFT(ISHFT(x(i), 8),-24)
     c(1,j+2)=ISHFT(ISHFT(x(i),16),-24)
     c(2,j+2)=ISHFT(ISHFT(x(i),24),-24)
  ENDDO
  nrest=MOD(n,2)
  IF (nrest.EQ.1) THEN
     c(1,n)=ISHFT(x(n/2+1),-24)
     c(2,n)=ISHFT(ISHFT(x(n/2+1), 8),-24)
  ENDIF
  !$omp parallel do private(I)
  DO i=1,n
     c(1,i)=(c(1,i)-128._real_8)*fr
     c(2,i)=(c(2,i)-128._real_8)*fi
  ENDDO
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE idcopy
! ==================================================================
SUBROUTINE sdcopy(n,x,y)
  ! ==--------------------------------------------------------------==
  ! == UNPACK N complex(8) :: NUMBERS FROM HALF PRECISION              ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
#if defined (__NEC)
  ! Arguments
  INTEGER :: n,x(*)
  REAL :: y(*)
  ! Variables
  INTEGER, PARAMETER :: nbig=10**9 
  INTEGER :: i,i1,i2
  ! ==--------------------------------------------------------------==
  !$omp parallel do private(I,I1,I2)
  DO i=1,n
     i1=2*i-1
     i2=i1+1
     y(i1)=ISHFT(x(i),-32)
     y(i2)=ISHFT(ISHFT(x(i),32),-32)
  ENDDO
  !$omp parallel do private(I)
  DO i=1,2*n
     y(i)=(y(i)-nbig)*1.e-8_real_8
  ENDDO
#else
  INTEGER :: n, i
  COMPLEX(real_4) :: x(n)
  COMPLEX(real_8) :: y(n)
  ! ==--------------------------------------------------------------==
  !$omp parallel do private(i)
  DO i=1,n
     y(i)=x(i)
  ENDDO
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sdcopy
! ==================================================================
SUBROUTINE dicopy(n,c,x,grmax,gimax)
  ! ==--------------------------------------------------------------==
  ! == PACK N complex(8) :: NUMBERS TO TWO 1 BYTE INTEGERS             ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
#if defined (__NEC) 
  ! Arguments
  INTEGER :: n,x(n/4+3)
  REAL(real_8) :: c(2,n),grmax,gimax
  ! Variables
  REAL(real_8) :: fr,fi
  INTEGER :: i,j,i11,i21,i12,i22,i13,i23,i14,i24,nrest
  ! ==--------------------------------------------------------------==
  fr=127._real_8/grmax
  fi=127._real_8/gimax
  DO i=1,n/4
     j=(i-1)*4
     i11=NINT(c(1,j+1)*fr)+128
     i21=NINT(c(2,j+1)*fi)+128
     i12=NINT(c(1,j+2)*fr)+128
     i22=NINT(c(2,j+2)*fi)+128
     i13=NINT(c(1,j+3)*fr)+128
     i23=NINT(c(2,j+3)*fi)+128
     i14=NINT(c(1,j+4)*fr)+128
     i24=NINT(c(2,j+4)*fi)+128
     x(i)=ISHFT(i11,56)+ISHFT(i21,48)+ISHFT(i12,40)+&
          ISHFT(i22,32)+ISHFT(i13,24)+ISHFT(i23,16)+&
          ISHFT(i14,8)+i24
  ENDDO
  nrest=MOD(n,4)
  IF (nrest.EQ.1) THEN
     i11=NINT(c(1,n)*fr)+128
     i21=NINT(c(2,n)*fi)+128
     x(n/4+1)=ISHFT(i11,56)+ISHFT(i21,48)
  ELSEIF (nrest.EQ.2) THEN
     i11=NINT(c(1,n-1)*fr)+128
     i21=NINT(c(2,n-1)*fi)+128
     i12=NINT(c(1,n)*fr)+128
     i22=NINT(c(2,n)*fi)+128
     x(n/4+1)=ISHFT(i11,56)+ISHFT(i21,48)+ISHFT(i12,40)+&
          ISHFT(i22,32)
  ELSEIF (nrest.EQ.3) THEN
     i11=NINT(c(1,n-2)*fr)+128
     i21=NINT(c(2,n-2)*fi)+128
     i12=NINT(c(1,n-1)*fr)+128
     i22=NINT(c(2,n-1)*fi)+128
     i13=NINT(c(1,n)*fr)+128
     i23=NINT(c(2,n)*fi)+128
     x(n/4+1)=ISHFT(i11,56)+ISHFT(i21,48)+ISHFT(i12,40)+&
          ISHFT(i22,32)+ISHFT(i13,24)+ISHFT(i23,16)
  ENDIF
#else
  INTEGER(int_4) :: n
  REAL(real_8) :: c(2,n),grmax,gimax, fr,fi
  INTEGER(int_4) :: x(*), i,j,nrest,i11,i21,i12,i22
  ! ==--------------------------------------------------------------==
  fr=127._real_8/grmax
  fi=127._real_8/gimax
#ifdef __SR11000
  !$omp parallel do private(I,J,I11,I21,I12,I22)
#endif 
  DO i=1,n/2
     j=(i-1)*2
     i11=NINT(c(1,j+1)*fr)+128
     i21=NINT(c(2,j+1)*fi)+128
     i12=NINT(c(1,j+2)*fr)+128
     i22=NINT(c(2,j+2)*fi)+128
     x(i)=ISHFT(i11,24)+ISHFT(i21,16)+&
          ISHFT(i12,8) +i22
  ENDDO
  nrest=MOD(n,2)
  IF (nrest.EQ.1) THEN
     i11=NINT(c(1,n)*fr)+128
     i21=NINT(c(2,n)*fi)+128
     x(n/2+1)=ISHFT(i11,24)+ISHFT(i21,16)
  ENDIF
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dicopy
! ==================================================================
SUBROUTINE grimax(c2,ngw_l,n,grmax,gimax,gnorm)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum,mp_max
  USE system , ONLY:spar
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: ngw_l, n
  REAL(real_8)                               :: c2(2,ngw_l*n), grmax, gimax, &
                                                gnorm

  EXTERNAL                                   :: dotp
  INTEGER                                    :: i, idamax, ii, ir
  REAL(real_8)                               :: dotp

! ==--------------------------------------------------------------==

  IF (ngw_l.GT.0.AND.n.GT.0) THEN
     ir=idamax(n*ngw_l,c2(1,1),2)
     ii=idamax(n*ngw_l,c2(2,1),2)
     grmax=c2(1,ir)
     gimax=c2(2,ii)
     gnorm=0.0_real_8
     DO i=1,n
        ii=(i-1)*ngw_l+1
        gnorm=gnorm+dotp(ngw_l,c2(1,ii),c2(1,ii))
     ENDDO
  ELSE
     grmax=0.0_real_8
     gimax=0.0_real_8
     gnorm=0.0_real_8
  ENDIF
  CALL mp_sum(gnorm,parai%cp_grp)
  GRMAX=ABS(GRMAX);GIMAX=ABS(GIMAX)
  CALL mp_max(grmax,parai%cp_grp)
  CALL mp_max(gimax,parai%cp_grp)
  ! vw protect agains division by zero
  IF ( gimax == 0.0_real_8 ) gimax = EPSILON( 0.0_real_8 )
  gnorm=SQRT(gnorm/REAL(n*spar%ngws,kind=real_8))
  RETURN
END SUBROUTINE grimax
! ==================================================================
SUBROUTINE dscopy(n,x,y)
  ! ==--------------------------------------------------------------==
  ! == PACK N complex(8) :: NUMBERS IN HALF PRECISION                  ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
#if defined (__NEC)
  ! Arguments
  INTEGER :: n,y(n)
  REAL :: x(2*n)
  ! Variables
  INTEGER, PARAMETER :: nbig=10**9 
  INTEGER :: i,i1,i2,ix1,ix2
  ! ==--------------------------------------------------------------==
  DO i=1,n
     i1=2*i-1
     i2=i1+1
     ix1=NINT(x(i1)*1.e8_real_8)+nbig
     ix2=NINT(x(i2)*1.e8_real_8)+nbig
     y(i)=ISHFT(ix1,32)+ix2
  ENDDO
#else
  INTEGER :: n, i
  COMPLEX(real_8) :: x(n)
  COMPLEX(real_4) :: y(n)
  ! ==--------------------------------------------------------------==
  !$omp parallel do private(i)
  DO i=1,n
     y(i)=x(i)
  ENDDO
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dscopy
! ==================================================================


! ==================================================================
SUBROUTINE odiis(c0,c2,vpp,nstate,pme,gde,svar2,reinit)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum, mp_bcast
  USE system , ONLY:cnti,cntl,maxdis,ncpw
  USE parac, ONLY : paral,parai
  USE ener , ONLY:ener_com
  USE elct , ONLY: crge
  USE cp_grp_utils, ONLY : cp_grp_redist
  USE cp_grp_utils, ONLY : cp_grp_get_sizes
  USE odiis_utils, ONLY : updis, solve
  USE zeroing_utils,                   ONLY: zeroing
  USE nvtx_utils
  IMPLICIT NONE
  REAL(real_8)                               :: vpp(*)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
  REAL(real_8)                               :: pme(ncpw%ngw*nstate+8,*), &
                                                gde((ncpw%ngw*nstate+8)/4,*), &
                                                svar2
  LOGICAL                                    :: reinit

  CHARACTER(*), PARAMETER                    :: procedureN = 'odiis'

  COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: c0_local, c2_local
  INTEGER                                    :: i, ibeg_c0, iend_c0, ierr, &
                                                ig, isub, isub3, j, k, &
                                                nempty, ngw_local, nowm, nsize
  INTEGER, SAVE                              :: istate = 0, ndiis, nocc, nowv
  LOGICAL                                    :: einc1, geq0_local
  REAL(real_8)                               :: bc(maxdis+1,maxdis+1), de1, &
                                                e1thr, ff, g1, g2, raim, &
                                                ratio, ration, vc(maxdis+1)
  REAL(real_8), EXTERNAL                     :: dotp
  REAL(real_8), SAVE                         :: diism(maxdis,maxdis), eold, &
                                                gamma, gimax(maxdis), &
                                                gnorm(maxdis), grmax(maxdis)

  CALL tiset(proceduren,isub)
  __NVTX_TIMER_START ( procedureN )

  IF (reinit.OR.istate.EQ.0) THEN
     ! Reinitialization of cntl%diis procedure
     istate=0
     eold=9999._real_8
     gamma=0._real_8
     reinit=.FALSE.
     nocc=0
     ndiis=0
     DO i=1,nstate
        IF (crge%f(i,1).GT.1.e-5_real_8) THEN
           nocc=nocc+1
        ENDIF
     ENDDO
  ENDIF
  ! >>>>>>> cp_grp trick
  CALL cp_grp_get_sizes(ngw_l=ngw_local,geq0_l=geq0_local,&
       first_g=ibeg_c0,last_g=iend_c0)
  ALLOCATE(c0_local(ngw_local,nstate),c2_local(ngw_local,nstate),&
       stat=ierr)
  IF (ierr.NE.0) CALL stopgm(proceduren,'Allocation problem',& 
       __LINE__,__FILE__)
  CALL cp_grp_copy_wfn_to_local(c0,ncpw%ngw,c0_local,ngw_local,&
       ibeg_c0,ngw_local,nstate)
  CALL cp_grp_copy_wfn_to_local(c2,ncpw%ngw,c2_local,ngw_local,&
       ibeg_c0,ngw_local,nstate)
  ! <<<<<<<
  ! Empirical rules to assure convergence
  IF (istate.NE.0) THEN
     ! Check for an increase in energy
     e1thr=1.0_real_8
     de1=ener_com%etot-eold
     einc1=.FALSE.
     IF (de1.GT.e1thr) einc1=.TRUE.
     IF (.NOT. cntl%tfrho_upw) THEN
        ! Trust region parameter
        IF (geq0_local) THEN
           IF (de1.GT.0.0_real_8) THEN
              gamma=gamma+1._real_8/vpp(1)
           ELSE
              gamma=0.75_real_8*gamma
           ENDIF
        ENDIF
     ENDIF
     CALL mp_bcast(gamma,parai%io_source,parai%cp_grp)
     reinit=.FALSE.
     IF (ndiis.GT.cnti%nreset.AND.cnti%nreset.GE.3) THEN
        ! ..restart if there was no progress over the last steps 
        nowm=nowv+1
        IF (nowv.EQ.cnti%mdiis) nowm=1
        g1=SQRT(grmax(nowv)**2 + gimax(nowv)**2)
        g2=SQRT(grmax(nowm)**2 + gimax(nowm)**2)
        ratio=g1/g2
        ration=gnorm(nowv)/gnorm(nowm)
        raim=1.0_real_8-REAL(cnti%mdiis,kind=real_8)/100._real_8
        IF (ratio.GT.raim.AND.ration.GT.raim) reinit=.TRUE.
        IF (reinit.AND.paral%io_parent)&
             WRITE(6,'(A)') ' ODIIS| Insufficient progress; reset! '
        IF (reinit) THEN
           istate=0
           gamma=0.0_real_8
           ! REINIT=.FALSE. ! we need this info in UPDWF, reset there
        ENDIF
     ENDIF
  ENDIF
  IF (istate.EQ.0) CALL updis(ndiis,nowv,nsize,cnti%mdiis,0)
  istate=2
  ! Perform an electronic cntl%diis step
  CALL updis(ndiis,nowv,nsize,cnti%mdiis,1)
  IF (gamma.NE.0.0_real_8) THEN
     !$omp parallel do private(ig)
     DO ig=ibeg_c0,iend_c0
        vpp(ig)=vpp(ig)/(1.0_real_8+vpp(ig)*gamma)
     ENDDO
     !$omp end parallel do
  ENDIF
  ! Update cntl%diis buffers
  CALL dscopy(ngw_local*nocc,c0_local,pme(1,nowv))
  CALL dscal(2*ngw_local*nocc,-1.0_real_8,c2_local,1)
  CALL grimax(c2_local,ngw_local,nocc,&
       grmax(nowv),gimax(nowv),gnorm(nowv))
  CALL dicopy(ngw_local*nocc,c2_local,gde(1,nowv),&
       grmax(nowv),gimax(nowv))

  !$omp parallel do private(K,FF,IG)
#ifdef __SR8000
  !poption parallel, tlocal(K,FF,IG)
#endif
  DO k=1,nocc
     IF (cntl%prec.AND.crge%f(k,1).GT.0.1_real_8) THEN
        ff=1.0_real_8/crge%f(k,1)
     ELSE
        ff=1.0_real_8
     ENDIF
     DO ig=1,ngw_local
        c2_local(ig,k)=c2_local(ig,k)*&
             2.0_real_8*ff*ff*vpp(ig+ibeg_c0-1)*vpp(ig+ibeg_c0-1)
     ENDDO
  ENDDO
  !$omp end parallel do

  ! Update cntl%diis matrix
  DO i=1,nsize-1
     diism(i,nowv)=0.0_real_8
  ENDDO
  IF (ngw_local.GT.0) THEN
     DO i=1,nsize-1
        CALL idcopy(ngw_local*nocc,c0_local,gde(1,i),&
             grmax(i),gimax(i))
        DO k=1,nocc
           diism(i,nowv)=diism(i,nowv)+&
                dotp(ngw_local,c0_local(1,k),c2_local(1,k))
        ENDDO
     ENDDO
  ENDIF
  CALL mp_sum(diism(:,nowv),nsize-1,parai%cp_grp)
  DO i=1,nsize-1
     diism(nowv,i)=diism(i,nowv)
  ENDDO
  ! Set up cntl%diis Matrix
  CALL zeroing(bc)!,(maxdis+1)*(maxdis+1))

  DO i=1,nsize-1
     DO j=1,nsize-1
        bc(i,j)=diism(i,j)
     ENDDO
  ENDDO
  DO i=1,nsize-1
     vc(i)=0._real_8
     bc(i,nsize)=-1._real_8
     bc(nsize,i)=-1._real_8
  ENDDO
  vc(nsize)=-1._real_8
  bc(nsize,nsize)=0._real_8
  ! Solve System of Linear Equations
  CALL solve(bc,maxdis+1,nsize,vc)
  ! Compute Interpolated Coefficient Vectors
#if defined(__VECTOR)
  CALL zeroing(c0_local)!,ngw_local*nocc)
  DO i=1,nsize-1
     CALL sdcopy(ngw_local*nocc,pme(1,i),c2_local)
     CALL daxpy(2*ngw_local*nocc,vc(i),c2_local(1,1),1,&
          c0_local(1,1),1)
  ENDDO
#else
  ! do the above in one go and avoid having to use C2.
  CALL intpcoef(ngw_local*nocc,nsize,&
       ncpw%ngw*nstate+8,pme,vc,c0_local)
#endif
  ! Estimate New Parameter Vectors 
  DO i=1,nsize-1
     CALL idcopy(ngw_local*nocc,c2_local,gde(1,i),&
          grmax(i),gimax(i))
     ff=1.0_real_8
     !$omp parallel do private(K,FF,IG)
#ifdef __SR8000
     !poption parallel, tlocal(K,FF,IG)
#endif 
     DO k=1,nocc
        IF (cntl%prec.AND.crge%f(k,1).GT.0.1_real_8) THEN
           ff=1.0_real_8/crge%f(k,1)
        ELSE
           ff=1.0_real_8
        ENDIF
        DO ig=1,ngw_local
           c0_local(ig,k)=c0_local(ig,k)-&
                vc(i)*ff*vpp(ig+ibeg_c0-1)*c2_local(ig,k)
        ENDDO
     ENDDO
     !$omp end parallel do
  ENDDO

  IF (nocc.LT.nstate) THEN
     ! Use steepest descent for empty states
     nempty=nstate-nocc
     CALL daxpy(2*ngw_local*nempty,svar2,c2_local(1,nocc+1),1,&
          c0_local(1,nocc+1),1)
  ENDIF
  eold=ener_com%etot
  IF (gamma.NE.0.0_real_8) THEN
     !$omp parallel do private(ig)
     DO ig=ibeg_c0,iend_c0
        vpp(ig)=1._real_8/(1.0_real_8/vpp(ig)-gamma)
     ENDDO
     !$omp end parallel do
  ENDIF

  ! >>>>>>> cp_grp trick
  CALL tiset(proceduren//'_grps_b',isub3)
  ! we need to zero the C0 that we can do the reduce
  CALL cp_grp_copy_local_to_wfn(c0_local,ngw_local,c0,ncpw%ngw,&
       ibeg_c0,ngw_local,nstate)
  ! we reduce so that we get back the cp_grp distribution of C0
  CALL cp_grp_zero_g(c0,ncpw%ngw,nstate,ibeg_c0,iend_c0)
  CALL cp_grp_redist(c0,ncpw%ngw,nstate)
  DEALLOCATE(c0_local,c2_local,stat=ierr)
  IF (ierr.NE.0) CALL stopgm(proceduren,'Deallocation problem',& 
       __LINE__,__FILE__)
  CALL tihalt(proceduren//'_grps_b',isub3)
  ! <<<<<<<

  __NVTX_TIMER_STOP
  CALL tihalt(proceduren,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE odiis
! ==================================================================
