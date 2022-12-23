MODULE sd_nl_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             parap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sd_nl

CONTAINS

  ! ==================================================================
  SUBROUTINE sd_nl(fion,f,fnl,fnl2,dfnl,fnl1,fnl3,dfnl1,nstate)
    ! ==---------------------------------------------------------------==
    ! ==                        computes                               ==
    ! ==  the non-local potential contribution to the dynamical matrix ==
    ! ==---------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), fnl2(ions1%nat,&
                                                maxsys%nhxs,*), &
                                                fnl3(ions1%nat,maxsys%nhxs,*)
    INTEGER                                  :: nstate
    REAL(real_8) :: dfnl1(ions1%nat,maxsys%nhxs,3,nstate), &
      fnl1(ions1%nat,maxsys%nhxs,nstate), &
      dfnl(ions1%nat,maxsys%nhxs,3,nstate), &
      fnl(ions1%nat,maxsys%nhxs,nstate), f(nstate)

    INTEGER                                  :: i, ia, iat, ii, is, iv, jv, &
                                                k, ki, kj, l, l2, li, lj
    REAL(real_8)                             :: temp, tt, tx

! variables
! ==--------------------------------------------------------------==

    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             tt=0._real_8
             IF (pslo_com%tvan(is)) THEN
                CALL stopgm('SD_NL','VDB not implemented',& 
                     __LINE__,__FILE__)
             ELSEIF (sgpp1%tsgp(is)) THEN
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
                            tx=2.0_real_8*f(i)*sgpp2%hlsg(ki,kj,l,is)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            IF (cntl%tfdist) THEN
                               tt=tt+tx*(&
                                    fnl3(iat,iv,ii)*dfnl(iat,jv,k,ii)+&
                                    dfnl1(iat,iv,k,ii)*fnl2(iat,jv,ii))
                            ELSE
                               tt=tt+tx*(&
                                    fnl1(iat,iv,i)*dfnl(iat,jv,k,ii)+&
                                    dfnl1(iat,iv,k,ii)*fnl(iat,jv,i))
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   temp=2._real_8*wsg(is,iv)
                   DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                      ii=i-parap%nst12(parai%mepos,1)+1
                      IF (cntl%tfdist) THEN
                         tt=tt+temp*f(i)*(&
                              fnl3(iat,iv,ii)*dfnl(iat,iv,k,ii)+&
                              dfnl1(iat,iv,k,ii)*fnl2(iat,iv,ii))
                      ELSE
                         tt=tt+temp*f(i)*(&
                              fnl1(iat,iv,i)*dfnl(iat,iv,k,ii)+&
                              dfnl1(iat,iv,k,ii)*fnl(iat,iv,i))
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             fion(k,ia,is)=fion(k,ia,is)+tt
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sd_nl
  ! ==================================================================

END MODULE sd_nl_utils
