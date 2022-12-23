MODULE epr_efg_utils
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             inyh,&
                                             rhops
  USE eicalc_utils,                    ONLY: eicalc
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_mark,&
                                             fo_verb
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop5,&
                                             rho_save_efg,&
                                             rho_save_epr
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_efg
  !!public :: calc_gradient
  PUBLIC :: apply_fpibg2
  PUBLIC :: apply_p_d
  PUBLIC :: save_rho
  PUBLIC :: epr_efg_outp
  PUBLIC :: trace
  PUBLIC :: epr_efg_outp_ev
  PUBLIC :: calcion

CONTAINS

  ! ============================================================
  SUBROUTINE epr_efg (rhoe,psi,nfi)
    ! ============================================================

    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: psi(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: nfi

    IF (prop5%teprefg) THEN
       CALL calc_gradient(rho_save_epr,rhoe,psi,1,nfi)
    ENDIF
    IF (prop5%tefg) THEN
       CALL calc_gradient(rho_save_efg,rhoe,psi,2,nfi)
    ENDIF
    RETURN
  END SUBROUTINE epr_efg


  ! ==========================================================      
  SUBROUTINE apply_fpibg2 (src)
    ! ==========================================================

    COMPLEX(real_8)                          :: src(ncpw%nhg)

    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: x, xx

    xx = fpi / parm%tpiba2
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       x = xx / hg(ig)
       src(ig) = src(ig) * x
    ENDDO
    RETURN
  END SUBROUTINE apply_fpibg2

  ! ==========================================================      
  SUBROUTINE apply_p_d (src,dest,coordinate)
    ! ==========================================================

    COMPLEX(real_8)                          :: src(ncpw%nhg), dest(ncpw%nhg)
    INTEGER                                  :: coordinate

    INTEGER                                  :: ig
    REAL(real_8)                             :: coeff, im, re

    DO ig=1,ncpw%nhg
       coeff    =   parm%tpiba  * gk(coordinate,ig)
       re       =   -coeff*AIMAG(src(ig))
       im       =   coeff*REAL(src(ig))
       dest(ig) =   CMPLX(re,im,kind=real_8)
       ! dest(ig) =   CMPLX(0.0_real_8,coeff)  * src(ig)
    ENDDO

    RETURN
  END SUBROUTINE apply_p_d


  ! ==============================================================
  SUBROUTINE save_rho (rhoe)
    ! ==============================================================

    REAL(real_8)                             :: rhoe(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'save_rho'

    INTEGER                                  :: ierr
    INTEGER, SAVE                            :: ifirst = 0

    IF (prop5%teprefg.OR.prop5%tefg) THEN
       IF (ifirst .EQ. 0) THEN
          ifirst = 1
          IF (prop5%teprefg) THEN
             ALLOCATE(rho_save_epr(fpar%nnr1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(rho_save_epr)!,nnr1)
          ENDIF
          IF (prop5%tefg) THEN
             ALLOCATE(rho_save_efg(fpar%nnr1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(rho_save_efg)!,nnr1)
          ENDIF
       ENDIF
       IF (prop5%teprefg) THEN
          CALL dcopy(fpar%nnr1,rhoe(1,1),1,rho_save_epr,1)
          CALL daxpy(fpar%nnr1,-2._real_8,rhoe(1,2),1,rho_save_epr,1)
       ENDIF
       IF (prop5%tefg) THEN
          CALL dcopy(fpar%nnr1,rhoe(1,1),1,rho_save_efg,1)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE save_rho

  ! ============================================================
  SUBROUTINE epr_efg_outp (m,f,p,i,nfi)
    ! ============================================================ 

    REAL(real_8)                             :: m(3,3,maxsys%nax,*), &
                                                f(3,maxsys%nax,*), &
                                                p(maxsys%nax,*)
    INTEGER                                  :: i, nfi

    CHARACTER(len=3)                         :: filename
    INTEGER                                  :: iat, isa, isp
    INTEGER, SAVE                            :: ifirst = -1
    LOGICAL                                  :: err

    IF (paral%io_parent) THEN
       IF (i.EQ.1) THEN
          filename='EPR'
       ELSEIF (i.EQ.2) THEN
          filename='EFG'
       ENDIF
       IF (ifirst.LT.0) ifirst=fo_verb+fo_mark
       CALL fileopen(95,filename,fo_app+ifirst,err)
       ifirst=0
       IF (i.EQ.1) THEN
10        FORMAT(a3,1x,'matrix elements '//&
               '(1,1) (2,2) (3,3) (1,2) (1,3) (2,3) rho(0)',/,1x)
          WRITE(95,10) filename
       ELSEIF (i.EQ.2) THEN
15        FORMAT(a3,1x,'matrix elements '//&
               '(1,1) (2,2) (3,3) (1,2) (1, 3) (2,3)  Ex Ey Ez',/,1x)
          WRITE(95,15) filename
       ENDIF
       WRITE(95,*) nfi
       isa = 0
       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             isa = isa + 1
             IF (i.EQ.1) THEN
                WRITE(95,20) isa,m(1,1,iat,isp),m(2,2,iat,isp),&
                     m(3,3,iat,isp),m(2,1,iat,isp),m(3,1,iat,isp),&
                     m(3,2,iat,isp),p(iat,isp)
20              FORMAT(i3,3x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,&
                     4x,f13.10)
             ELSEIF (i.EQ.2) THEN
                WRITE(95,25) isa,m(1,1,iat,isp),m(2,2,iat,isp),m(3,3,iat,&
                     isp),m(2,1,iat,isp),m(3,1,iat,isp),m(3,2,iat,isp),&
                     f(1,iat,isp),f(2,iat,isp),f(3,iat,isp),p(iat,isp)
25              FORMAT(i3,3x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,&
                     4x,f8.5,2x,f8.5,2x,f8.5,3x,f13.10)
             ENDIF
          ENDDO
       ENDDO
       CALL fileclose(95)
    ENDIF
    RETURN
  END SUBROUTINE epr_efg_outp

  ! ===============================================================
  SUBROUTINE trace (matrix)
    ! ===============================================================

    REAL(real_8)                             :: matrix(3,3,maxsys%nax,*)

    INTEGER                                  :: i, iat, isp
    REAL(real_8)                             :: sum

    IF (paral%parent) THEN
       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             sum=0._real_8
             DO i=1,3
                sum=sum+matrix(i,i,iat,isp)
             ENDDO
             DO i=1,3
                matrix(i,i,iat,isp)=matrix(i,i,iat,isp)-sum/3
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE trace

  ! =============================================================
  SUBROUTINE epr_efg_outp_ev (eigenvalues,i,nfi)
    ! =============================================================

    REAL(real_8)                             :: eigenvalues(3,maxsys%nax,*)
    INTEGER                                  :: i, nfi

    CHARACTER(len=6)                         :: filename
    INTEGER                                  :: iat, isa, isp
    INTEGER, SAVE                            :: ifirst = -1
    LOGICAL                                  :: err

    IF (paral%parent) THEN
       IF (i.EQ.1) THEN
          filename='EV_EPR'
       ELSEIF (i.EQ.2) THEN
          filename='EV_EFG'
       ENDIF

       IF (ifirst.LT.0) ifirst=fo_verb+fo_mark
       IF (paral%io_parent) CALL fileopen(95,filename,fo_app+ifirst,err)
       ifirst=0
       IF (i.EQ.1) THEN
          IF (paral%io_parent) WRITE(95,10) filename
10        FORMAT(a3,1x,'Anisotropic Spin Dipole Couplings '//&
               ' Eigenvalues (principal axis system) ',/,1x)
       ELSEIF (i.EQ.2) THEN
          IF (paral%io_parent) WRITE(95,15) filename
15        FORMAT(a3,1x,'Electric field gradient '//&
               'Eigenvalues (principal axis system) ',/,1x)
       ENDIF
    ENDIF
    IF (paral%io_parent) WRITE(95,*) nfi
    isa = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          isa = isa + 1
          IF (paral%io_parent) WRITE(95,20) isa,eigenvalues(1,iat,isp),&
               eigenvalues(2,iat,isp),eigenvalues(3,iat,isp)
20        FORMAT(i3,3x,f8.4,2x,f8.4,2x,f8.4)
       ENDDO
    ENDDO
    IF (paral%io_parent) CALL fileclose(95)
    RETURN
  END SUBROUTINE epr_efg_outp_ev

  ! =============================================================
  SUBROUTINE calcion (eirop,isa,is)
    ! =============================================================
    COMPLEX(real_8)                          :: eirop(ncpw%nhg)
    INTEGER                                  :: isa, is

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: ig
    REAL(real_8)                             :: ei, er

    CALL zeroing(eirop)!,nhg)
    DO ig=1,ncpw%nhg
       ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
            ei3(isa,inyh(3,ig))
       er=REAL(ei123)
       ei=AIMAG(ei123)
       eirop(ig)=eirop(ig)+CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)
    ENDDO

    RETURN
  END SUBROUTINE calcion


  ! =========================================================================


END MODULE epr_efg_utils


! =======================================================
SUBROUTINE calc_gradient (rhor,rhoe,psi,i,nfi)
  ! =======================================================  
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:maxsys,ncpw,fpar
  USE parac, ONLY : paral,parai
  USE ions , ONLY:ions0,ions1
  USE spin , ONLY:clsd
  USE epr_efg_utils, ONLY : calcion,apply_fpibg2,apply_p_d,trace,&
       & epr_efg_outp,epr_efg_outp_ev
  USE eicalc_utils, ONLY : eicalc
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  REAL(real_8)                               :: rhor(fpar%nnr1)
  COMPLEX(real_8)                            :: rhoe(ncpw%nhg), &
                                                psi(fpar%nnr1,clsd%nlsd)
  INTEGER                                    :: i, nfi

  CHARACTER(*), PARAMETER                    :: procedureN = 'calc_gradient'

  COMPLEX(real_8), ALLOCATABLE               :: scr(:,:)
  INTEGER                                    :: dir, dir1, iat, ierr, ig, &
                                                info, isa, isp
  INTEGER, SAVE                              :: ifirst = 0
  REAL(real_8)                               :: eigr_dot, work(3*3)
  REAL(real_8), ALLOCATABLE, SAVE            :: efield(:,:,:), &
                                                eigenvalues(:,:,:), &
                                                gradient_matrix(:,:,:,:), &
                                                pot(:,:)

  ALLOCATE(scr(ncpw%nhg, 3),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)

  IF (ifirst.EQ.0) THEN
     ifirst=1
     ! TODO deallocate - ?
     ALLOCATE(gradient_matrix(3,3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(efield(3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(pot(maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(eigenvalues(3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(gradient_matrix)!,3*3*maxsys%nax*maxsys%nsx)
     CALL zeroing(efield)!,3*maxsys%nax*maxsys%nsx)
     CALL zeroing(pot)!,maxsys%nax*maxsys%nsx)
     CALL zeroing(eigenvalues)!,3*maxsys%nax*maxsys%nsx)
  ENDIF
  CALL ffttog(rhor,rhoe,psi,ncpw%nhg,.TRUE.)
  IF (i.EQ.2) THEN
     CALL eicalc(scr(:,1),scr(:,2))
     DO ig=1,ncpw%nhg
        rhoe(ig)=rhoe(ig)+scr(ig,2)
     ENDDO
  ENDIF
  isa = 0
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        isa = isa + 1
        DO ig=1,ncpw%nhg
           scr(ig,1)=rhoe(ig)
        ENDDO
        IF (i.EQ.2) THEN
           CALL calcion(scr(1,2),isa,isp)
           DO ig=1,ncpw%nhg
              scr(ig,1)=scr(ig,1)-scr(ig,2)
           ENDDO
        ENDIF
        CALL apply_fpibg2 (scr(1,1))
        IF (i.EQ.1) THEN
           pot(iat,isp)=eigr_dot(isa,rhoe)
        ELSE
           pot(iat,isp)=(-1._real_8)*eigr_dot(isa,scr(1,1))
        ENDIF
        DO dir=1,3
           CALL apply_p_d (scr(1,1),scr(1,2),dir)
           IF (i.EQ.2) THEN
              efield(dir,iat,isp)=eigr_dot(isa,scr(1,2))
           ENDIF
           DO dir1=dir,3
              CALL apply_p_d (scr(1,2),scr(1,3),dir1)
              IF ((i.EQ.2).AND.(dir1.EQ.dir)) THEN
                 efield(dir,iat,isp)=eigr_dot(isa,scr(1,2))
              ENDIF
              gradient_matrix(dir1,dir,iat,isp) &
                   =eigr_dot(isa,scr(1,3))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(scr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  CALL mp_sum(gradient_matrix,3*3*maxsys%nax*maxsys%nsx,parai%allgrp)
  CALL mp_sum(efield,3*maxsys%nax*maxsys%nsx,parai%allgrp)
  CALL mp_sum(pot,maxsys%nax*maxsys%nsx,parai%allgrp)
  IF (i.EQ.2) CALL trace(gradient_matrix)
  CALL epr_efg_outp(gradient_matrix,efield,pot,i,nfi)
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        CALL dsyev('N','L',3,gradient_matrix(1,1,iat,isp),3,&
             eigenvalues(1,iat,isp),work,9,info)
        IF (paral%io_parent.AND.(info.NE.0)) THEN
           WRITE (6, *)'ERROR IN DSYEV, code ',info
        ENDIF
     ENDDO
  ENDDO
  CALL epr_efg_outp_ev(eigenvalues,i,nfi)
  RETURN
END SUBROUTINE calc_gradient
