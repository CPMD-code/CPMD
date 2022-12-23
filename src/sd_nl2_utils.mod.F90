MODULE sd_nl2_utils
  USE elct,                            ONLY: crge
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
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             parap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sd_nl2

CONTAINS

  ! ==================================================================
  SUBROUTINE sd_nl2(sder,iato,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE DYNAMICAL MATRIX ==
    ! == THERE ARE ONLY DIAGONAL CONTRIBUTIONS                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder(3*ions1%nat,3*ions1%nat)
    INTEGER                                  :: iato(maxsys%nax,maxsys%nsx), &
                                                nstate

    INTEGER                                  :: i, ia, ii, is, isa, isa0, iv, &
                                                j, jj, jv, ki, kj, l, l2, li, &
                                                lj
    REAL(real_8)                             :: temp, tt, weight

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    isa0=0
    DO is=1,ions1%nsp
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
                      DO ia=1,ions0%na(is)
                         jj=iato(ia,is)*3
                         isa=isa0+ia
                         sder(jj+1,jj+1)=sder(jj+1,jj+1)+tt*&
                              (ddfnl(isa,iv,1,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,jv,1,ii,1))
                         sder(jj+1,jj+2)=sder(jj+1,jj+2)+tt*&
                              (ddfnl(isa,iv,2,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,jv,2,ii,1))
                         sder(jj+1,jj+3)=sder(jj+1,jj+3)+tt*&
                              (ddfnl(isa,iv,3,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,jv,3,ii,1))
                         sder(jj+2,jj+2)=sder(jj+2,jj+2)+tt*&
                              (ddfnl(isa,iv,4,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,jv,2,ii,1))
                         sder(jj+2,jj+3)=sder(jj+2,jj+3)+tt*&
                              (ddfnl(isa,iv,5,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,jv,3,ii,1))
                         sder(jj+3,jj+3)=sder(jj+3,jj+3)+tt*&
                              (ddfnl(isa,iv,6,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,jv,3,ii,1))
                         sder(jj+2,jj+1)=sder(jj+2,jj+1)+tt*&
                              (ddfnl(isa,iv,2,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,jv,1,ii,1))
                         sder(jj+3,jj+1)=sder(jj+3,jj+1)+tt*&
                              (ddfnl(isa,iv,3,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,jv,1,ii,1))
                         sder(jj+3,jj+2)=sder(jj+3,jj+2)+tt*&
                              (ddfnl(isa,iv,5,ii)*fnl2(1,isa,jv,j,1)+&
                              dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,jv,2,ii,1))
                      ENDDO
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
                DO ia=1,ions0%na(is)
                   jj=iato(ia,is)*3
                   isa=isa0+ia
                   sder(jj+1,jj+1)=sder(jj+1,jj+1)+weight*&
                        (ddfnl(isa,iv,1,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,iv,1,ii,1))
                   sder(jj+1,jj+2)=sder(jj+1,jj+2)+weight*&
                        (ddfnl(isa,iv,2,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,iv,2,ii,1))
                   sder(jj+1,jj+3)=sder(jj+1,jj+3)+weight*&
                        (ddfnl(isa,iv,3,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,1,ii,1)*dfnl(1,isa,iv,3,ii,1))
                   sder(jj+2,jj+2)=sder(jj+2,jj+2)+weight*&
                        (ddfnl(isa,iv,4,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,iv,2,ii,1))
                   sder(jj+2,jj+3)=sder(jj+2,jj+3)+weight*&
                        (ddfnl(isa,iv,5,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,iv,3,ii,1))
                   sder(jj+3,jj+3)=sder(jj+3,jj+3)+weight*&
                        (ddfnl(isa,iv,6,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,iv,3,ii,1))
                   sder(jj+2,jj+1)=sder(jj+2,jj+1)+weight*&
                        (ddfnl(isa,iv,2,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,2,ii,1)*dfnl(1,isa,iv,1,ii,1))
                   sder(jj+3,jj+1)=sder(jj+3,jj+1)+weight*&
                        (ddfnl(isa,iv,3,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,iv,1,ii,1))
                   sder(jj+3,jj+2)=sder(jj+3,jj+2)+weight*&
                        (ddfnl(isa,iv,5,ii)*fnl2(1,isa,iv,j,1)+&
                        dfnl(1,isa,iv,3,ii,1)*dfnl(1,isa,iv,2,ii,1))
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sd_nl2
  ! ==================================================================

END MODULE sd_nl2_utils
