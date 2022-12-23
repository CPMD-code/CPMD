MODULE pi_stress_utils
  USE cnst,                            ONLY: au_kb,&
                                             factem
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_vmark
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: paral,&
                                             parai
  USE pbc_utils,                       ONLY: pbc
  USE pimd,                            ONLY: pimd3,&
                                             ipp1,&
                                             pmar,&
                                             pmars,&
                                             pma0s,&
                                             grandparent
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE store_types,                     ONLY: rout1
  USE symtrz_utils,                    ONLY: symstress
  USE system,                          ONLY: cnti,&
                                             cntr,&
                                             iatpt,&
                                             parm,&
                                             cntl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_stress_vir
  PUBLIC :: pi_stress_nmm
  PUBLIC :: pi_stress_pri
  PUBLIC :: wr_pi_stress

CONTAINS

  ! ==================================================================
  FUNCTION pi_stress_vir(paiuall,stagep,vstage,pitaup,fionks,tinfo) RESULT(paiu)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: paiuall(:,:,:),&
                                                stagep(:,:,:),&
                                                vstage(:,:,:),&
                                                pitaup(:,:,:,:),&
                                                fionks(:,:,:,:),&
                                                paiu(3,3)
    LOGICAL                                  :: tinfo
    INTEGER                                  :: npx, ip, ia, is, i, j
    REAL(real_8)                             :: fact, rnp, sum, xlm, ylm, zlm, xlm_, ylm_, zlm_, &
      paiks(3,3), paifr(3,3), paikin(3,3)
    REAL(real_8), PARAMETER                  :: epsilon = 1.e-15_real_8

    IF (.NOT.paral%io_parent) GOTO 999
    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    ! Electronic term
    CALL zeroing(paiks)
    DO ip=1,npx
       paiks(:,:)=paiks(:,:)+paiuall(:,:,ip)
    ENDDO
    ! Force term other than centroids 
    CALL zeroing(paifr)
    !$omp parallel do private(i,ia,is,ip,xlm,ylm,zlm) reduction(+:paifr)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       DO ip=1,npx
          xlm_=pitaup(1,ia,is,ip)-stagep(1,ia,is)
          ylm_=pitaup(2,ia,is,ip)-stagep(2,ia,is)
          zlm_=pitaup(3,ia,is,ip)-stagep(3,ia,is)
          CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
          paifr(1,1)=paifr(1,1)-fionks(1,ia,is,ip)*xlm
          paifr(2,1)=paifr(2,1)-fionks(2,ia,is,ip)*xlm
          paifr(3,1)=paifr(3,1)-fionks(3,ia,is,ip)*xlm
          paifr(1,2)=paifr(1,2)-fionks(1,ia,is,ip)*ylm
          paifr(2,2)=paifr(2,2)-fionks(2,ia,is,ip)*ylm
          paifr(3,2)=paifr(3,2)-fionks(3,ia,is,ip)*ylm
          paifr(1,3)=paifr(1,3)-fionks(1,ia,is,ip)*zlm
          paifr(2,3)=paifr(2,3)-fionks(2,ia,is,ip)*zlm
          paifr(3,3)=paifr(3,3)-fionks(3,ia,is,ip)*zlm
       ENDDO
    ENDDO
    !$omp end parallel do
    ! Symmetric matrix and set to 0 if < EPSILON
    DO i=1,3
       DO j=i,3
          IF (ABS(paifr(j,i)).LT.epsilon) THEN
             paifr(j,i)=0._real_8
             paifr(i,j)=0._real_8
          ELSE
             sum=0.5_real_8*(paifr(j,i)+paifr(i,j))
             paifr(j,i)=sum
             paifr(i,j)=sum
          ENDIF
       ENDDO
    ENDDO
    ! Kinetic term of centroids only
    CALL zeroing(paikin)
    !$omp parallel do private(i,ia,is,fact) reduction(+:paikin)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       fact=pma0s(is,1)
       paikin(1,1)=paikin(1,1)+fact*vstage(1,ia,is)*vstage(1,ia,is)
       paikin(2,1)=paikin(2,1)+fact*vstage(2,ia,is)*vstage(1,ia,is)
       paikin(3,1)=paikin(3,1)+fact*vstage(3,ia,is)*vstage(1,ia,is)
       paikin(1,2)=paikin(1,2)+fact*vstage(1,ia,is)*vstage(2,ia,is)
       paikin(2,2)=paikin(2,2)+fact*vstage(2,ia,is)*vstage(2,ia,is)
       paikin(3,2)=paikin(3,2)+fact*vstage(3,ia,is)*vstage(2,ia,is)
       paikin(1,3)=paikin(1,3)+fact*vstage(1,ia,is)*vstage(3,ia,is)
       paikin(2,3)=paikin(2,3)+fact*vstage(2,ia,is)*vstage(3,ia,is)
       paikin(3,3)=paikin(3,3)+fact*vstage(3,ia,is)*vstage(3,ia,is)
    ENDDO
    !$omp end parallel do
    paiu=paiks+paifr   ! stress tensor minus velocity-dependent part
    ! Symmetrization
    CALL symstress(paiu)
    ! htfp
    CALL dcopy(9,paiu(1,1),1,metr_com%htfp(1,1),1)
    IF (prcpl%tzflex) THEN
       metr_com%htfp(3,3)=metr_com%htfp(3,3)-prcp_com%druck*parm%omega
    ELSE IF (prcpl%tisot.OR.cntl%tshock) THEN
       DO i=1,3
          metr_com%htfp(i,i)=metr_com%htfp(i,i)-prcp_com%druck*parm%omega
       ENDDO
    ELSE
       CALL daxpy(9,-parm%omega,prcp_com%stens(1,1),1,metr_com%htfp(1,1),1)
    ENDIF
    ! htfor
    CALL dgemm('N','N',3,3,3,1.0_real_8,metr_com%htfp(1,1),3,metr_com%htm1(1,1),3,&
         0.0_real_8,metr_com%htfor(1,1),3)
    paiu=paiu+paikin   ! true stress tensor
    ! Print out the components if requested
    IF (grandparent.AND.tinfo.AND.pimd3%levprt>=3) THEN
       WRITE(6,'(3(4X,8("-"),A,8("-")))') &
       '   PAIKS  ','   PAIFR  ','  PAIKIN  '
       DO i=1,3
         WRITE(6,'(9F10.3)') paiks(i,:)*au_kb/parm%omega,&
                             paifr(i,:)*au_kb/parm%omega,&
                             paikin(i,:)*au_kb/parm%omega
       ENDDO
    ENDIF
999 CONTINUE
    CALL mp_bcast(paiu,9,parai%io_source,parai%cp_grp)
    CALL mp_bcast(metr_com%htfp,9,parai%io_source,parai%cp_grp)
    CALL mp_bcast(metr_com%htfor,9,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION pi_stress_vir
  ! ==================================================================
  FUNCTION pi_stress_nmm(paiuall,stagep,vstage,tinfo) RESULT(paiu)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: paiuall(:,:,:),&
                                                stagep(:,:,:,:),&
                                                vstage(:,:,:,:),&
                                                paiu(3,3)
    LOGICAL                                  :: tinfo
    INTEGER                                  :: ip, ia, is, i, j, npx
    REAL(real_8)                             :: fact, rnp, sum, omegap, omp2, &
      paiks(3,3), paihar(3,3), paikin(3,3)
    REAL(real_8), PARAMETER                  :: epsilon = 1.e-15_real_8

    IF (.NOT.paral%io_parent) GOTO 999
    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    ! Electronic term
    CALL zeroing(paiks)
    DO ip=1,npx
       paiks(:,:)=paiks(:,:)+paiuall(:,:,ip)
    ENDDO
    ! Harmonic term
    omegap = SQRT(rnp)*cntr%tempw/factem
    omp2 = omegap*omegap
    CALL zeroing(paihar)
    !$omp parallel do private(i,ia,is,ip,fact) reduction(+:paihar)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       DO ip=2,npx
          fact=-pmars(is,ip)*omp2
          paihar(1,1)=paihar(1,1)+fact*stagep(1,ia,is,ip)*stagep(1,ia,is,ip)
          paihar(2,1)=paihar(2,1)+fact*stagep(2,ia,is,ip)*stagep(1,ia,is,ip)
          paihar(3,1)=paihar(3,1)+fact*stagep(3,ia,is,ip)*stagep(1,ia,is,ip)
          paihar(1,2)=paihar(1,2)+fact*stagep(1,ia,is,ip)*stagep(2,ia,is,ip)
          paihar(2,2)=paihar(2,2)+fact*stagep(2,ia,is,ip)*stagep(2,ia,is,ip)
          paihar(3,2)=paihar(3,2)+fact*stagep(3,ia,is,ip)*stagep(2,ia,is,ip)
          paihar(1,3)=paihar(1,3)+fact*stagep(1,ia,is,ip)*stagep(3,ia,is,ip)
          paihar(2,3)=paihar(2,3)+fact*stagep(2,ia,is,ip)*stagep(3,ia,is,ip)
          paihar(3,3)=paihar(3,3)+fact*stagep(3,ia,is,ip)*stagep(3,ia,is,ip)
       ENDDO
    ENDDO
    !$omp end parallel do
    ! Kinetic term 
    CALL zeroing(paikin)
    !$omp parallel do private(i,ia,is,ip,fact) reduction(+:paikin)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       DO ip=1,npx
          fact=pma0s(is,ip)
          paikin(1,1)=paikin(1,1)+fact*vstage(1,ia,is,ip)*vstage(1,ia,is,ip)
          paikin(2,1)=paikin(2,1)+fact*vstage(2,ia,is,ip)*vstage(1,ia,is,ip)
          paikin(3,1)=paikin(3,1)+fact*vstage(3,ia,is,ip)*vstage(1,ia,is,ip)
          paikin(1,2)=paikin(1,2)+fact*vstage(1,ia,is,ip)*vstage(2,ia,is,ip)
          paikin(2,2)=paikin(2,2)+fact*vstage(2,ia,is,ip)*vstage(2,ia,is,ip)
          paikin(3,2)=paikin(3,2)+fact*vstage(3,ia,is,ip)*vstage(2,ia,is,ip)
          paikin(1,3)=paikin(1,3)+fact*vstage(1,ia,is,ip)*vstage(3,ia,is,ip)
          paikin(2,3)=paikin(2,3)+fact*vstage(2,ia,is,ip)*vstage(3,ia,is,ip)
          paikin(3,3)=paikin(3,3)+fact*vstage(3,ia,is,ip)*vstage(3,ia,is,ip)
       ENDDO
    ENDDO
    !$omp end parallel do 
    paiu=paiks+paihar+paikin
    ! Symmetrization
    CALL symstress(paiu)
    ! Print out the components if requested
    IF (grandparent.AND.tinfo.AND.pimd3%levprt>=3) THEN
       WRITE(6,'(3(4X,8("-"),A,8("-")))') &
       '   PAIKS  ','  PAIHAR  ','  PAIKIN  '
       DO i=1,3
         WRITE(6,'(9F10.3)') paiks(i,:)*au_kb/parm%omega,&
                             paihar(i,:)*au_kb/parm%omega,&
                             paikin(i,:)*au_kb/parm%omega
       ENDDO
    ENDIF
999 CONTINUE
    CALL mp_bcast(paiu,9,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION pi_stress_nmm
  ! ==================================================================
  FUNCTION pi_stress_pri(paiuall,pitaup,pivelp,tinfo) RESULT(paiu)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: paiuall(:,:,:),&
                                                pitaup(:,:,:,:),&
                                                pivelp(:,:,:,:),&
                                                paiu(3,3)
    LOGICAL                                  :: tinfo
    INTEGER                                  :: ip, ia, is, i, j, npx
    REAL(real_8)                             :: fact, rnp, xlm, ylm, zlm, sum, &
      omegap, omp2, paiks(3,3), paihar(3,3), paikin(3,3), xlm_, ylm_, zlm_
    REAL(real_8), PARAMETER                  :: epsilon = 1.e-15_real_8

    IF (.NOT.paral%io_parent) GOTO 999
    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    ! Electronic term
    CALL zeroing(paiks)
    DO ip=1,npx
       paiks(:,:)=paiks(:,:)+paiuall(:,:,ip)
    ENDDO
    ! Harmonic term
    omegap = SQRT(rnp)*cntr%tempw/factem
    omp2 = omegap*omegap
    CALL zeroing(paihar)
    !$omp parallel do private(i,ia,is,fact,ip,xlm,ylm,zlm) reduction(+:paihar)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       fact=-pmar(is)*omp2
       DO ip=1,npx
          xlm_=pitaup(1,ia,is,ip)-pitaup(1,ia,is,ipp1(ip))
          ylm_=pitaup(2,ia,is,ip)-pitaup(2,ia,is,ipp1(ip))
          zlm_=pitaup(3,ia,is,ip)-pitaup(3,ia,is,ipp1(ip))
          CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
          paihar(1,1)=paihar(1,1)+fact*xlm*xlm
          paihar(2,1)=paihar(2,1)+fact*ylm*xlm
          paihar(3,1)=paihar(3,1)+fact*zlm*xlm
          paihar(1,2)=paihar(1,2)+fact*xlm*ylm
          paihar(2,2)=paihar(2,2)+fact*ylm*ylm
          paihar(3,2)=paihar(3,2)+fact*zlm*ylm
          paihar(1,3)=paihar(1,3)+fact*xlm*zlm
          paihar(2,3)=paihar(2,3)+fact*ylm*zlm
          paihar(3,3)=paihar(3,3)+fact*zlm*zlm
       ENDDO
    ENDDO
    !$omp end parallel do
    ! Kinetic term 
    CALL zeroing(paikin)
    !$omp parallel do private(i,ia,is,fact,ip) reduction(+:paikin)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       fact=rmass%pma(is)
       DO ip=1,npx
          paikin(1,1)=paikin(1,1)+fact*pivelp(1,ia,is,ip)*pivelp(1,ia,is,ip)
          paikin(2,1)=paikin(2,1)+fact*pivelp(2,ia,is,ip)*pivelp(1,ia,is,ip)
          paikin(3,1)=paikin(3,1)+fact*pivelp(3,ia,is,ip)*pivelp(1,ia,is,ip)
          paikin(1,2)=paikin(1,2)+fact*pivelp(1,ia,is,ip)*pivelp(2,ia,is,ip)
          paikin(2,2)=paikin(2,2)+fact*pivelp(2,ia,is,ip)*pivelp(2,ia,is,ip)
          paikin(3,2)=paikin(3,2)+fact*pivelp(3,ia,is,ip)*pivelp(2,ia,is,ip)
          paikin(1,3)=paikin(1,3)+fact*pivelp(1,ia,is,ip)*pivelp(3,ia,is,ip)
          paikin(2,3)=paikin(2,3)+fact*pivelp(2,ia,is,ip)*pivelp(3,ia,is,ip)
          paikin(3,3)=paikin(3,3)+fact*pivelp(3,ia,is,ip)*pivelp(3,ia,is,ip)
       ENDDO
    ENDDO
    !$omp end parallel do
    paiu=paiks+paihar+paikin
    ! Symmetrization
    CALL symstress(paiu)
    ! Print out the components if requested
    IF (grandparent.AND.tinfo.AND.pimd3%levprt>=3) THEN
       WRITE(6,'(3(4X,8("-"),A,8("-")))') &
       '   PAIKS  ','  PAIHAR  ','  PAIKIN  '
       DO i=1,3
         WRITE(6,'(9F10.3)') paiks(i,:)*au_kb/parm%omega,&
                             paihar(i,:)*au_kb/parm%omega,&
                             paikin(i,:)*au_kb/parm%omega
       ENDDO
    ENDIF
999 CONTINUE
    CALL mp_bcast(paiu,9,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION pi_stress_pri
  ! ==================================================================
  SUBROUTINE wr_pi_stress(paiu)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: paiu(3,3)
    CHARACTER(len=10), PARAMETER             :: f4 = 'STRESS    ', &
                                                f5 = 'CELL      '
    INTEGER                                  :: i
    INTEGER, SAVE                            :: if2 = fo_vmark
    LOGICAL                                  :: ferror

    CALL fileopen(33,f4,fo_app+if2,ferror)
    WRITE(33,*) '  TOTAL STRESS TENSOR (kB): Step:', iteropt%nfi
    DO i=1,3
       WRITE(33,'(5X,3(F20.8))') (paiu(i,:)/parm%omega)*au_kb
    ENDDO
    CALL fileclose(33)

    IF (cntl%tprcp.AND.(ropt_mod%rprint.OR.&
       (ropt_mod%txyz.AND.rout1%xtout).OR.(ropt_mod%tdcd.AND.rout1%dcout))) THEN
       CALL fileopen(32,f5,fo_app+if2,ferror)
       WRITE(32,*) '  CELL PARAMETERS at Step:', iteropt%nfi
       DO i=1,3
          WRITE(32,'(3(1X,F14.6),8X,3(1X,F12.6))')&
          metr_com%ht(i,:),metr_com%htvel(i,:)
       ENDDO
       CALL fileclose(32)
    ENDIF
    if2=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr_pi_stress
  ! ==================================================================

END MODULE pi_stress_utils
