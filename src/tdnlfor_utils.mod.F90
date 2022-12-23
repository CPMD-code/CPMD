MODULE tdnlfor_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tdnlfor

CONTAINS

  ! ==================================================================
  SUBROUTINE tdnlfor(fion,f,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)

    INTEGER                                  :: i, ia, ii, is, isa, isa0, &
                                                isub, iv, jv, k, ki, kj, l, &
                                                l2, li, lj
    REAL(real_8)                             :: temp, tt

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('   TDNLFOR',isub)
    DO k=1,3
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             CALL stopgm("TDNLFOR","Not supported",& 
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
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         tt=f(i)*sgpp2%hlsg(ki,kj,l,is)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         IF (cntl%tfdist) THEN
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               fion(k,ia,is)=fion(k,ia,is)-tt*&
                                    dfnl(1,isa,iv,k,ii,1)*fnl2(1,isa,jv,ii,1)
                            ENDDO
                         ELSE
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               fion(k,ia,is)=fion(k,ia,is)-tt*&
                                    dfnl(1,isa,iv,k,ii,1)*fnl2(1,isa,jv,i,1)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! Other pp (numeric)
             DO iv=1,nlps_com%ngh(is)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   temp=wsg(is,iv)*f(i)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   IF (cntl%tfdist) THEN
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         fion(k,ia,is)=fion(k,ia,is)-temp*&
                              dfnl(1,isa,iv,k,ii,1)*fnl2(1,isa,iv,ii,1)
                      ENDDO
                   ELSE
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         fion(k,ia,is)=fion(k,ia,is)-temp*&
                              dfnl(1,isa,iv,k,ii,1)*fnl2(1,isa,iv,i,1)
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0 = isa0 + ions0%na(is)
       ENDDO
    ENDDO
    CALL tihalt('   TDNLFOR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tdnlfor
  ! ==================================================================

END MODULE tdnlfor_utils
