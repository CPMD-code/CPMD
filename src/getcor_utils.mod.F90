MODULE getcor_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pbc_utils,                       ONLY: pbc
  USE pimd,                            ONLY: fcorr,&
                                             fionks,&
                                             ifcorr,&
                                             ircorr,&
                                             pimd3,&
                                             rcorr
  USE system,                          ONLY: parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getcor

CONTAINS

  ! ==================================================================
  SUBROUTINE getcor(taup)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taup(:,:,:,:)

    INTEGER                                  :: ia, iipp, ip, ipp, is, itau
    REAL(real_8)                             :: var, xlm, ylm, zlm , xlm_, ylm_, zlm_

! Variables
! ==--------------------------------------------------------------==
! ==  CALCULATE IMAGINARY TIME CORRELATION OF POSITIONS           == 
! ==--------------------------------------------------------------==

    CALL zeroing(rcorr)!,maxsys%nax*maxsys%nsx*pimd3%np_total)
    CALL zeroing(ircorr)!,maxsys%nax*maxsys%nsx*pimd3%np_total)
    DO ip=1,pimd3%np_total
       itau=0
       DO ipp=ip,ip+pimd3%np_total-1
          iipp=ipp-pimd3%np_total*(ipp/(pimd3%np_total+1))
          itau=itau+1
          IF (itau.LE.pimd3%np_total/2+1) THEN
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   xlm_=taup(1,ia,is,ip)-taup(1,ia,is,iipp)
                   ylm_=taup(2,ia,is,ip)-taup(2,ia,is,iipp)
                   zlm_=taup(3,ia,is,ip)-taup(3,ia,is,iipp)
                   CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
                   var=xlm**2+ylm**2+zlm**2
                   rcorr(ia,is,itau)=rcorr(ia,is,itau)+var
                   ircorr(ia,is,itau)=ircorr(ia,is,itau)+1
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE IMAGINARY TIME CORRELATION OF KS-FORCES           ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(fcorr)!,maxsys%nax*maxsys%nsx*pimd3%np_total)
    CALL zeroing(ifcorr)!,maxsys%nax*maxsys%nsx*pimd3%np_total)
    DO ip=1,pimd3%np_total
       itau=0
       DO ipp=ip,ip+pimd3%np_total-1
          iipp=ipp-pimd3%np_total*(ipp/(pimd3%np_total+1))
          itau=itau+1
          IF (itau.LE.pimd3%np_total/2+1) THEN
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   xlm=fionks(1,ia,is,ip)-fionks(1,ia,is,iipp)
                   ylm=fionks(2,ia,is,ip)-fionks(2,ia,is,iipp)
                   zlm=fionks(3,ia,is,ip)-fionks(3,ia,is,iipp)
                   var=xlm**2+ylm**2+zlm**2
                   fcorr(ia,is,itau)=fcorr(ia,is,itau)+var
                   ifcorr(ia,is,itau)=ifcorr(ia,is,itau)+1
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getcor
  ! ==================================================================

END MODULE getcor_utils
