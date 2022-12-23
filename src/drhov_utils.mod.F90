MODULE drhov_utils
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

!!$private
!!$
!!$public :: drhov
!!$
!!$contains

END MODULE drhov_utils


! ==================================================================
! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- DRHOV_C WILL BE WRITTEN   ==
! == NO CALL TO THIS ROUTINE                                      ==
! ==================================================================
SUBROUTINE drhov(nstate,qg,ctmp)
  ! ==--------------------------------------------------------------==
  ! ==    THIS ROUTINE CALCULATES THE DERIVATIVE VANDERBILT DENSITY ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE prmem_utils, ONLY: prmem
  USE reshaper, ONLY: reshape_inplace
  USE system , ONLY:cntl,ncpw,parap
  USE parac, ONLY : paral,parai
  USE ions , ONLY:ions0,ions1
  USE pslo , ONLY:pslo_com
  USE nlps , ONLY:nlps_com
  USE elct , ONLY:crge
  USE cppt , ONLY:inyh
  USE sfac , ONLY:ei1,ei2,ei3,eigrb,fnl
  USE spin , ONLY:spin_mod
  USE str2 , ONLY:becs,dqg,drhovg,gagk,qrada
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: qg(ncpw%nhg), ctmp(ncpw%nhg), &
                                                ctmps(ncpw%nhg,2)

  CHARACTER(*), PARAMETER                    :: procedureN = 'drhov'

  INTEGER                                    :: i, ia, ierr, ig, ii, is, isa, &
                                                isa0, isub, iv, jv, kk
  INTEGER, SAVE                              :: ifirst = 0
  REAL(real_8)                               :: dsum, sum, &
                                                dsum_up, dsum_down, &
                                                sum_up, sum_down
  REAL(real_8), ALLOCATABLE, SAVE            :: spline(:)
  COMPLEX(real_8), POINTER                   :: dqg_2d(:,:)

  CALL tiset('     DRHOV',isub)
  ! ==--------------------------------------------------------------==
  IF (cntl%tfdist) CALL stopgm('DRHOV','TFDIST NOT IMPLEMENTED',& 
       __LINE__,__FILE__)
  ! ==-------------------------------------------------------------==
  ! ==  ALLOCATE MEMORY                                            ==
  ! ==-------------------------------------------------------------==
  IF (ifirst.EQ.0) THEN
     ifirst=1
     ALLOCATE(spline(ncpw%nhg),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     IF (paral%parent) CALL prmem('     DRHOV')
  ENDIF
  ! ==--------------------------------------------------------------==
  CALL reshape_inplace(dqg, (/ncpw%nhg,6/), dqg_2d)
  CALL zeroing(drhovg)!,nhg*6)
  isa0=0
  DO is=1,ions1%nsp
     IF (pslo_com%tvan(is)) THEN
        DO iv=1,nlps_com%ngh(is)
           DO jv=iv,nlps_com%ngh(is)
              CALL dqvan2(iv,jv,is,qrada,ctmp,gagk,qg,dqg_2d,spline)
              IF (cntl%bigmem) THEN
                 IF (cntl%tlsd) THEN
                    CALL zeroing(ctmps)
                    DO ia=1,ions0%na(is)
                       isa=isa0+ia
                       sum=0.0_real_8
                       DO i=1,spin_mod%nsup
                          sum = sum + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum=2._real_8*sum
                       CALL daxpy(2*ncpw%nhg,sum,eigrb(1,isa),1,ctmps(1,1),1)
                       sum=0.0_real_8
                       DO i=spin_mod%nsup+1,spin_mod%nsup+spin_mod%nsdown
                          sum = sum + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum=2._real_8*sum
                       CALL daxpy(2*ncpw%nhg,sum,eigrb(1,isa),1,ctmps(1,2),1)
                    ENDDO
                    !$omp parallel do private(kk,ig)
                    DO kk=1,6
                       DO ig=1,ncpw%nhg
                          drhovg(ig,kk)=drhovg(ig,kk)+dqg_2d(ig,kk)*ctmps(ig,1)
                          drhovg(ig,6+kk)=drhovg(ig,6+kk)+dqg_2d(ig,kk)*ctmps(ig,2)
                       ENDDO
                    ENDDO
                 ELSE
                    CALL zeroing(ctmp)!,nhg)
                    DO ia=1,ions0%na(is)
                       isa=isa0+ia
                       sum=0.0_real_8
                       DO i=1,nstate
                          sum = sum + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum=2._real_8*sum
                       CALL daxpy(2*ncpw%nhg,sum,eigrb(1,isa),1,ctmp(1),1)
                    ENDDO
                    !$omp parallel do private(KK,IG)
                    DO kk=1,6
                       DO ig=1,ncpw%nhg
                          drhovg(ig,kk)=drhovg(ig,kk)+dqg_2d(ig,kk)*ctmp(ig)
                       ENDDO
                    ENDDO
                 ENDIF
              ELSE
                 DO ia=1,ions0%na(is)
                    isa=isa0+ia
                    !$omp parallel do private(IG)
#ifdef __SR8000
                    !poption parallel
#endif
                    DO ig=1,ncpw%nhg
                       ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                            ei3(isa,inyh(3,ig))
                    ENDDO
                    IF (cntl%tlsd) THEN
                       sum_up=0.0_real_8
                       DO i=1,spin_mod%nsup
                          sum_up = sum_up + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum_up=2._real_8*sum_up
                       sum_down=0.0_real_8
                       DO i=spin_mod%nsup+1,spin_mod%nsup+spin_mod%nsdown
                          sum_down = sum_down + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum_down=2._real_8*sum_down
                       !$omp parallel do private(kk,ig)
                       DO kk=1,6
                          DO ig=1,ncpw%nhg
                             drhovg(ig,kk)=drhovg(ig,kk)+dqg_2d(ig,kk)*sum_up*ctmp(ig)
                             drhovg(ig,6+kk)=drhovg(ig,6+kk)+dqg_2d(ig,kk)*sum_down*ctmp(ig)
                          ENDDO
                       ENDDO
                    ELSE
                       sum=0.0_real_8
                       DO i=1,nstate
                          sum = sum + crge%f(i,1)*&
                               fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                       ENDDO
                       IF (iv.NE.jv) sum=2._real_8*sum
                       !$omp parallel do private(KK,IG)
                       DO kk=1,6
                          DO ig=1,ncpw%nhg
                             drhovg(ig,kk)=drhovg(ig,kk)+dqg_2d(ig,kk)*sum*ctmp(ig)
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
              IF (cntl%bigmem) THEN
                 DO kk=1,6
                    IF (cntl%tlsd) THEN
                       CALL zeroing(ctmps)!,nhg)
                       DO ia=1,ions0%na(is)
                          isa=isa0+ia
                          dsum_up=0.0_real_8; dsum_down=0.0_real_8
                          DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                             ii=i-parap%nst12(parai%mepos,1)+1
                             IF (i>=1.AND.i<=spin_mod%nsup) THEN
                                dsum_up=dsum_up +&
                                   crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                   crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                             ENDIF
                             IF (i>=spin_mod%nsup+1.AND.i<=spin_mod%nsup+spin_mod%nsdown) THEN
                                dsum_down=dsum_down +&
                                   crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                   crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                             ENDIF
                          ENDDO
                       ENDDO
                       CALL mp_sum(dsum_up,parai%allgrp)
                       CALL mp_sum(dsum_down,parai%allgrp)
                       IF (iv.NE.jv) dsum_up=2._real_8*dsum_up
                       IF (iv.NE.jv) dsum_down=2._real_8*dsum_down
                       CALL daxpy(2*ncpw%nhg,dsum_up,eigrb(1,isa),1,ctmps(1,1),1)
                       CALL daxpy(2*ncpw%nhg,dsum_down,eigrb(1,isa),1,ctmps(1,2),1)
                       DO ig=1,ncpw%nhg
                          drhovg(ig,kk)=drhovg(ig,kk)+qg(ig)*ctmps(ig,1)
                          drhovg(ig,6+kk)=drhovg(ig,6+kk)+qg(ig)*ctmps(ig,2)
                       ENDDO
                    ELSE
                       CALL zeroing(ctmp)!,nhg)
                       DO ia=1,ions0%na(is)
                          isa=isa0+ia
                          dsum=0.0_real_8
                          DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                             ii=i-parap%nst12(parai%mepos,1)+1
                             dsum=dsum +&
                                  crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                  crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                          ENDDO
                          CALL mp_sum(dsum,parai%allgrp)
                          IF (iv.NE.jv) dsum=2._real_8*dsum
                          CALL daxpy(2*ncpw%nhg,dsum,eigrb(1,isa),1,ctmp(1),1)
                       ENDDO
                       DO ig=1,ncpw%nhg
                          drhovg(ig,kk)=drhovg(ig,kk)+qg(ig)*ctmp(ig)
                       ENDDO
                    ENDIF
                 ENDDO
              ELSE
                 DO ia=1,ions0%na(is)
                    isa=isa0+ia
                    !$omp parallel do private(IG)
#ifdef __SR8000
                    !poption parallel
#endif
                    DO ig=1,ncpw%nhg
                       ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                            ei3(isa,inyh(3,ig))
                    ENDDO
                    DO kk=1,6
                       IF (cntl%tlsd) THEN
                          dsum_up=0.0_real_8; dsum_down=0.0_real_8
                          DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                             ii=i-parap%nst12(parai%mepos,1)+1
                             IF (i>=1.AND.i<=spin_mod%nsup) THEN
                                dsum_up=dsum_up +&
                                   crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                   crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                             ENDIF
                             IF (i>=spin_mod%nsup+1.AND.i<=spin_mod%nsup+spin_mod%nsdown) THEN
                                dsum_down=dsum_down +&
                                   crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                   crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                             ENDIF
                          ENDDO
                          CALL mp_sum(dsum_up,parai%allgrp)
                          CALL mp_sum(dsum_down,parai%allgrp)
                          IF (iv.NE.jv) dsum_up=2._real_8*dsum_up
                          IF (iv.NE.jv) dsum_down=2._real_8*dsum_down
                          DO ig=1,ncpw%nhg
                             drhovg(ig,kk)=drhovg(ig,kk)+qg(ig)*dsum_up*ctmp(ig)
                             drhovg(ig,6+kk)=drhovg(ig,6+kk)+qg(ig)*dsum_down*ctmp(ig)
                          ENDDO
                       ELSE
                          dsum=0.0_real_8
                          DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                             ii=i-parap%nst12(parai%mepos,1)+1
                             dsum=dsum +&
                                  crge%f(i,1)*fnl(1,isa,iv,i,1)*becs(1,isa,jv,kk,ii,1)+&
                                  crge%f(i,1)*becs(1,isa,iv,kk,ii,1)*fnl(1,isa,jv,i,1)
                          ENDDO
                          CALL mp_sum(dsum,parai%allgrp)
                          IF (iv.NE.jv) dsum=2._real_8*dsum
                          DO ig=1,ncpw%nhg
                             drhovg(ig,kk)=drhovg(ig,kk)+qg(ig)*dsum*ctmp(ig)
                          ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDIF
     isa0=isa0+ions0%na(is)
  ENDDO
  CALL tihalt('     DRHOV',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE drhov
! ==================================================================

