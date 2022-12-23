MODULE g_loc_opeigr_utils

  IMPLICIT NONE

  PRIVATE

  !!public :: g_loc_opeigr

  !contains

END MODULE g_loc_opeigr_utils


SUBROUTINE g_loc_opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
     nop1,nop2)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE opeigr_utils, ONLY : opapp
  USE system , ONLY: ncpw,spar,parap
  USE parac, ONLY : paral,parai
  USE ddip , ONLY:lenbk,ngg1,ngg2,ngg3,ngwmax
  USE g_loc , ONLY:lostate
  USE reshaper , ONLY:reshape_inplace
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  COMPLEX(real_8)                            :: c0(2*ncpw%ngw,*)
  COMPLEX(real_8), TARGET                    :: c2(:,:), sc0(:,:)
  INTEGER                                    :: nstate, mapful(2,*), &
                                                mapcol(ngwmax,*)
  COMPLEX(real_8)                            :: ddmat(nstate,nstate)
  INTEGER                                    :: nop1, nop2

  CHARACTER(*), PARAMETER                    :: procedureN = 'g_loc_opeigr'

  COMPLEX(real_8), ALLOCATABLE               :: aux(:)
  COMPLEX(real_8), POINTER                   :: c2f(:,:), c2s(:,:), sc0s(:,:)
  INTEGER                                    :: i, ierr, ig, igii, ii, ip, &
                                                n1, nggp, nn

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

  CALL reshape_inplace(c2, (/2*spar%ngws, nstate/), c2f)
  CALL zeroing(ddmat)!,nstate*nstate)
  IF (lostate%state_all) THEN
     CALL reshape_inplace(c2, (/2*lenbk, parai%nproc/), c2s)
     CALL reshape_inplace(sc0, (/2*lenbk, parai%nproc/), sc0s)

     DO ip=1,parai%nproc
        nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
        n1=parap%nst12(ip-1,1)
        CALL dcopy(2*2*ncpw%ngw*nn,c0(1,n1),1,c2s(1,ip),1)
     ENDDO

     CALL my_trans(c2s,sc0s,16*2*lenbk,1)
     nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
     DO ip=1,parai%nproc
        nggp = parap%sparm(3,ip-1)
        igii = 0
        DO ii=1,nn
           DO ig=1,nggp
              igii=igii+1
              c2f(mapcol(ig,ip)     ,ii) = sc0s(igii     ,ip)
              c2f(mapcol(ig,ip)+spar%ngws,ii) = sc0s(igii+nggp,ip)
           ENDDO
           igii = igii+nggp
        ENDDO
     ENDDO
  ELSE
     CALL dcopy(2*2*ncpw%ngw*nstate,c0,1,c2f,1)
  ENDIF

  ALLOCATE(aux(ngg1*ngg2*ngg3),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  IF (lostate%state_all) THEN
     DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
        ii=i-parap%nst12(parai%mepos,1)+1
        CALL zeroing(aux)!,SIZE(aux))

        DO ig=1,spar%ngws
           aux(mapful(1,ig))=c2f(ig     ,ii)
           aux(mapful(2,ig))=c2f(ig+spar%ngws,ii)
        ENDDO
        CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
        CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
        DO ig=1,spar%ngws
           sc0(ig     ,ii) = aux(mapful(1,ig))
           sc0(ig+spar%ngws,ii) = aux(mapful(2,ig))
        ENDDO
        sc0(1+spar%ngws,ii)  = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
     ENDDO
     nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
     DO ip=1,parai%nproc
        nggp = parap%sparm(3,ip-1)
        igii = 0
        DO ii=1,nn
           DO ig=1,nggp
              igii=igii+1
              c2s(igii     ,ip) = sc0(mapcol(ig,ip)     ,ii)
              c2s(igii+nggp,ip) = sc0(mapcol(ig,ip)+spar%ngws,ii)
           ENDDO
           igii=igii+nggp
        ENDDO
     ENDDO
     CALL my_trans(c2s,sc0s,16*2*lenbk,1)
     DO ip=1,parai%nproc
        nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
        n1=parap%nst12(ip-1,1)
        CALL dcopy(2*2*ncpw%ngw*nn,sc0s(1,ip),1,c2(1,n1),1)
     ENDDO
  ELSE
     DO ii = 1,nstate
        CALL zeroing(aux)!,SIZE(aux))

        DO ig=1,spar%ngws
           aux(mapful(1,ig))=c2f(ig     ,ii)
           aux(mapful(2,ig))=c2f(ig+spar%ngws,ii)
        ENDDO
        CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
        CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
        DO ig=1,spar%ngws
           sc0(ig     ,ii) = aux(mapful(1,ig))
           sc0(ig+spar%ngws,ii) = aux(mapful(2,ig))
        ENDDO
        sc0(1+spar%ngws,ii)  = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
     ENDDO

     CALL dcopy(2*2*ncpw%ngw*nstate,sc0,1,c2,1)

  ENDIF
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,2*ncpw%ngw,&
       c2,2*ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)

  CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)

  RETURN
END SUBROUTINE g_loc_opeigr

