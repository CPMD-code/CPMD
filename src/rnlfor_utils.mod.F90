MODULE rnlfor_utils
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,&
                                             fnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ipept,&
                                             maxsys,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlfor
  !public :: rcasfor

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlfor(fion,f,wk,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: f(nstate,nkpoint), wk(nkpoint)

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                ispin, isub, iv, jv, ki, kj, &
                                                l, l2, li, lj
    REAL(real_8)                             :: tdbl, temp, tt, weight, &
                                                wk1_1, wk1_2, wk1_3, wk2_1, &
                                                wk2_2, wk2_3

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('    RNLFOR',isub)
    DO ik=1,nkpoint
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ! Vanderbild pp
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   tdbl=1._real_8
                   IF (iv.NE.jv) tdbl=2.0_real_8
                   DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                      ii=i-parap%nst12(parai%mepos,1)+1
                      ispin=1
                      weight=wk(ik)*f(i,ik)
                      IF (ABS(weight).GT.1.e-12_real_8) THEN
                         IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
                         IF (cntl%tfdist) THEN
                            IF (imagp.EQ.2) THEN
                               !$omp parallel do private(IA,ISA,TEMP,wk1_1,wk1_2,wk1_3) &
                               !$omp             private(wk2_1,wk2_2,wk2_3)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  temp=dvan(iv,jv,is)+deeq(isa,jv,iv,ispin)
                                  wk1_1=dfnl(1,isa,jv,1,ii,ik)*fnl2(1,isa,iv,&
                                       ii,ik)+dfnl(2,isa,jv,1,ii,ik)*fnl2(2,isa,&
                                       iv,  ii,ik)
                                  wk1_2=dfnl(1,isa,jv,2,ii,ik)*fnl2(1,isa,iv,&
                                       ii,ik)+dfnl(2,isa,jv,2,ii,ik)*fnl2(2,isa,&
                                       iv,  ii,ik)
                                  wk1_3=dfnl(1,isa,jv,3,ii,ik)*fnl2(1,isa,iv,ii,&
                                       ik)+dfnl(2,isa,jv,3,ii,ik)*fnl2(2,isa,iv,&
                                       ii,ik)
                                  wk2_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,jv,&
                                       ii,ik)+dfnl(2,isa,iv,1,ii,ik)*fnl2(2,isa,&
                                       jv,  ii,ik)
                                  wk2_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,jv,&
                                       ii,ik)+dfnl(2,isa,iv,2,ii,ik)*fnl2(2,isa,&
                                       jv,  ii,ik)
                                  wk2_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,jv,&
                                       ii,ik)+dfnl(2,isa,iv,3,ii,ik)*fnl2(2,isa,&
                                       jv,  ii,ik)
                                  fion(1,ia,is)=fion(1,ia,is)-weight*tdbl*&
                                       temp*(wk1_1+wk2_1)
                                  fion(2,ia,is)=fion(2,ia,is)-weight*tdbl*&
                                       temp*(wk1_2+wk2_2)
                                  fion(3,ia,is)=fion(3,ia,is)-weight*tdbl*&
                                       temp*(wk1_3+wk2_3)
                               ENDDO
                            ELSE
                               !$omp parallel do private(IA,ISA,TEMP,wk1_1,wk1_2,wk1_3) &
                               !$omp             private(wk2_1,wk2_2,wk2_3)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  temp=dvan(iv,jv,is)+deeq(isa,jv,iv,ispin)
                                  wk1_1=dfnl(1,isa,jv,1,ii,ik) *fnl2(1,isa,iv,&
                                       ii,ik)
                                  wk1_2=dfnl(1,isa,jv,2,ii,ik) *fnl2(1,isa,iv,&
                                       ii,ik)
                                  wk1_3=dfnl(1,isa,jv,3,ii,ik) *fnl2(1,isa,iv,&
                                       ii,ik)
                                  wk2_1=dfnl(1,isa,iv,1,ii,ik) *fnl2(1,isa,jv,&
                                       ii,ik)
                                  wk2_2=dfnl(1,isa,iv,2,ii,ik) *fnl2(1,isa,jv,&
                                       ii,ik)
                                  wk2_3=dfnl(1,isa,iv,3,ii,ik) *fnl2(1,isa,jv,&
                                       ii,ik)
                                  fion(1,ia,is)=fion(1,ia,is)-weight*tdbl*&
                                       temp*(wk1_1+wk2_1)
                                  fion(2,ia,is)=fion(2,ia,is)-weight*tdbl*&
                                       temp*(wk1_2+wk2_2)
                                  fion(3,ia,is)=fion(3,ia,is)-weight*tdbl*&
                                       temp*(wk1_3+wk2_3)
                               ENDDO
                            ENDIF
                         ELSE
                            IF (imagp.EQ.2) THEN
                               !$omp parallel do private(IA,ISA,TEMP,wk1_1,wk1_2,wk1_3) &
                               !$omp             private(wk2_1,wk2_2,wk2_3)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  temp=dvan(iv,jv,is)+deeq(isa,jv,iv,ispin)
                                  wk1_1=dfnl(1,isa,jv,1,ii,ik)*fnl(1,isa,iv, i,&
                                       ik)+dfnl(2,isa,jv,1,ii,ik)*fnl(2,isa,iv,&
                                       i,ik)
                                  wk1_2=dfnl(1,isa,jv,2,ii,ik)*fnl(1,isa,iv, i,&
                                       ik)+dfnl(2,isa,jv,2,ii,ik)*fnl(2,isa,iv,&
                                       i,ik)
                                  wk1_3=dfnl(1,isa,jv,3,ii,ik)*fnl(1,isa,iv, i,&
                                       ik)+dfnl(2,isa,jv,3,ii,ik)*fnl(2,isa,iv,&
                                       i,ik)
                                  wk2_1=dfnl(1,isa,iv,1,ii,ik)*fnl(1,isa,jv, i,&
                                       ik)+dfnl(2,isa,iv,1,ii,ik)*fnl(2,isa,jv,&
                                       i,ik)
                                  wk2_2=dfnl(1,isa,iv,2,ii,ik)*fnl(1,isa,jv, i,&
                                       ik)+dfnl(2,isa,iv,2,ii,ik)*fnl(2,isa,jv,&
                                       i,ik)
                                  wk2_3=dfnl(1,isa,iv,3,ii,ik)*fnl(1,isa,jv, i,&
                                       ik)+dfnl(2,isa,iv,3,ii,ik)*fnl(2,isa,jv,&
                                       i,ik)
                                  fion(1,ia,is)=fion(1,ia,is)-weight*tdbl*&
                                       temp*(wk1_1+wk2_1)
                                  fion(2,ia,is)=fion(2,ia,is)-weight*tdbl*&
                                       temp*(wk1_2+wk2_2)
                                  fion(3,ia,is)=fion(3,ia,is)-weight*tdbl*&
                                       temp*(wk1_3+wk2_3)
                               ENDDO
                            ELSE
                               !$omp parallel do private(IA,ISA,TEMP,wk1_1,wk1_2,wk1_3) &
                               !$omp             private(wk2_1,wk2_2,wk2_3)
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  temp=dvan(iv,jv,is)+deeq(isa,jv,iv,ispin)
                                  wk1_1=dfnl(1,isa,jv,1,ii,ik)*fnl(1,isa,iv,&
                                       i,ik)
                                  wk1_2=dfnl(1,isa,jv,2,ii,ik)*fnl(1,isa,iv,&
                                       i,ik)
                                  wk1_3=dfnl(1,isa,jv,3,ii,ik)*fnl(1,isa,iv,&
                                       i,ik)
                                  wk2_1=dfnl(1,isa,iv,1,ii,ik)*fnl(1,isa,jv,&
                                       i,ik)
                                  wk2_2=dfnl(1,isa,iv,2,ii,ik)*fnl(1,isa,jv,&
                                       i,ik)
                                  wk2_3=dfnl(1,isa,iv,3,ii,ik)*fnl(1,isa,jv,&
                                       i,ik)
                                  fion(1,ia,is)=fion(1,ia,is)-weight*tdbl*&
                                       temp*(wk1_1+wk2_1)
                                  fion(2,ia,is)=fion(2,ia,is)-weight*tdbl*&
                                       temp*(wk1_2+wk2_2)
                                  fion(3,ia,is)=fion(3,ia,is)-weight*tdbl*&
                                       temp*(wk1_3+wk2_3)
                               ENDDO
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker pp
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         weight=wk(ik)*f(i,ik)
                         IF (ABS(weight).GT.1.e-12_real_8) THEN
                            tt=2.0_real_8*weight*sgpp2%hlsg(ki,kj,l,is)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            IF (imagp.EQ.2) THEN
                               IF (cntl%tfdist) THEN
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,1,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,2,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,3,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                  ENDDO
                               ELSE
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,1,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,2,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,3,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                  ENDDO
                               ENDIF
                            ELSE
                               IF (cntl%tfdist) THEN
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                  ENDDO
                               ELSE
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                  ENDDO
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! Other pp (numeric)
             DO iv=1,nlps_com%ngh(is)
                temp=2._real_8*wsg(is,iv)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   weight=wk(ik)*f(i,ik)
                   IF (ABS(weight).GT.1.e-12_real_8) THEN
                      ii=i-parap%nst12(parai%mepos,1)+1
                      IF (cntl%tfdist) THEN
                         IF (imagp.EQ.2) THEN
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,1,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,2,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,3,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ELSE
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ENDIF
                      ELSE
                         IF (imagp.EQ.2) THEN
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,1,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,2,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,3,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ELSE
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0 = isa0 + ions0%na(is)
       ENDDO
    ENDDO
    CALL tihalt('    RNLFOR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlfor
  ! ==================================================================
  SUBROUTINE rcasfor(fion)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rcasfor'

    INTEGER                                  :: ia, iatl, ierr, ii, is, isa, &
                                                isa0, isub, iv, jv, k, ki, &
                                                kj, l, l2, li, lj
    REAL(real_8)                             :: tt
    REAL(real_8), ALLOCATABLE                :: dfab(:,:,:)

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('   RCASFOR',isub)
    ALLOCATE(dfab(ions1%nat,maxsys%nhxs,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO k=1,3
       CALL zeroing(dfab)!,2*ions1%nat)
       IF (clsd%ialpha.GE.parap%nst12(parai%mepos,1).AND.clsd%ialpha.LE.parap%nst12(parai%mepos,2) ) THEN
          ii=clsd%ialpha-parap%nst12(parai%mepos,1)+1
          CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,1),1)
       ENDIF
       IF (clsd%ibeta.GE.parap%nst12(parai%mepos,1).AND.clsd%ibeta.LE.parap%nst12(parai%mepos,2) ) THEN
          ii=clsd%ibeta-parap%nst12(parai%mepos,1)+1
          CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,2),1)
       ENDIF
       CALL mp_sum(dfab,2*ions1%nat*maxsys%nhxs,parai%allgrp)
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ! Vanderbild pp
             CALL stopgm("RCASFOR","VDB NOT IMPLEMENTED",& 
                  __LINE__,__FILE__)
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker pp
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      tt=sgpp2%hlsg(ki,kj,l,is)
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         IF (cntl%tfdist) THEN
                            iatl=isa-ipept(1,parai%mepos)+1
                            fion(k,ia,is)=fion(k,ia,is)-0.5_real_8*tt*(dfab(isa,iv,&
                                 1)*fnl(1,iatl,jv,clsd%ibeta,1)+dfab(isa,iv,2)*&
                                 fnl(1,iatl,jv,clsd%ialpha,1)+ dfab(isa,jv,1)*&
                                 fnl(1,iatl,iv,clsd%ibeta,1)+dfab(isa,jv,2)*&
                                 fnl(1,iatl, iv, clsd%ialpha,1))
                         ELSE
                            fion(k,ia,is)=fion(k,ia,is)-tt*(dfab(isa,iv,1)*&
                                 fnl(1,isa,jv,clsd%ibeta,1)+dfab(isa,iv,2)*fnl(1,&
                                 isa, jv, clsd%ialpha, 1)+dfab(isa, jv,1)*fnl(1,&
                                 isa,iv, clsd%ibeta,1)+ dfab(isa,jv,2)*fnl(1,isa,&
                                 iv,clsd%ialpha,1))
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! Other pp (numeric)
             DO iv=1,nlps_com%ngh(is)
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   IF (isa.GE.ipept(1,parai%mepos).AND.isa.LE.ipept(2,parai%mepos))&
                        THEN
                      IF (cntl%tfdist) THEN
                         iatl=isa-ipept(1,parai%mepos)+1
                         fion(k,ia,is)=fion(k,ia,is)-wsg(is,iv)*(dfab(isa,iv,&
                              1)*fnl(1,iatl,iv,clsd%ibeta,1)+fnl(1,iatl,iv,clsd%ialpha,&
                              1)* dfab(isa,iv,2))
                      ELSE
                         fion(k,ia,is)=fion(k,ia,is)-wsg(is,iv)*(dfab(isa,iv,&
                              1)*fnl(1,isa,iv,clsd%ibeta,1)+fnl(1,isa,iv,clsd%ialpha,&
                              1)* dfab(isa,iv,2))
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0 = isa0 + ions0%na(is)
       ENDDO
    ENDDO
    DEALLOCATE(dfab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('   RCASFOR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rcasfor
  ! ==================================================================

END MODULE rnlfor_utils
