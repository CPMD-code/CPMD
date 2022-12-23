MODULE eind_nl_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             parap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eind_nl

CONTAINS

  ! ==================================================================
  SUBROUTINE eind_nl(eind,is,isa,ka,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE DYNAMICAL MATRIX ==
    ! == THERE ARE ONLY DIAGONAL CONTRIBUTIONS                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eind
    INTEGER                                  :: is, isa, ka, nstate

    INTEGER                                  :: i, ii, iv, j, jv, ki, kj, kk, &
                                                l, l2, li, lj
    REAL(real_8)                             :: temp, tt, weight

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    IF (ka.EQ.1) THEN
       kk = 1
    ELSEIF (ka.EQ.2) THEN
       kk = 4
    ELSE
       kk = 6
    ENDIF
    IF (pslo_com%tvan(is)) THEN
       ! ..Vanderbild pp
       CALL stopgm('SD_NL2','VDB PP NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ELSEIF (sgpp1%tsgp(is)) THEN
       ! ..Stefan Goedecker pp
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
                   tt=2.0_real_8*crge%f(i,1)*sgpp2%hlsg(ki,kj,l,is)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   IF (cntl%tfdist) THEN
                      j=ii
                   ELSE
                      j=i
                   ENDIF
                   eind=eind+tt*&
                        (ddfnl(isa,iv,kk,ii)*fnl2(1,isa,jv, j,1)+&
                        dfnl(1,isa,iv,ka,ii,1)*dfnl(1,isa,jv,ka,ii,1))
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ELSE
       ! ..Other pp (numeric)
       DO iv=1,nlps_com%ngh(is)
          temp=2._real_8*wsg(is,iv)
          DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
             weight=crge%f(i,1)*temp
             ii=i-parap%nst12(parai%mepos,1)+1
             IF (cntl%tfdist) THEN
                j=ii
             ELSE
                j=i
             ENDIF
             eind=eind+weight*&
                  (ddfnl(isa,iv,kk,ii)*fnl2(1,isa,iv, j,1)+&
                  dfnl(1,isa,iv,ka,ii,1)*dfnl(1,isa,iv,ka,ii,1))
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eind_nl
  ! ==================================================================

END MODULE eind_nl_utils
