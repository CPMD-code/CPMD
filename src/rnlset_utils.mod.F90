MODULE rnlset_utils
  USE atom,                            ONLY: atom_common,&
                                             gnl,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE dcacp_utils,                     ONLY: dcacp
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             rgh,&
                                             wgh,&
                                             wsg
  USE sgpp,                            ONLY: sgpp1
  USE system,                          ONLY: lmaxx,&
                                             nhx
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: simpsn

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlset

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlset(is,gmax,fint,wsgtmp)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE FORM FACTORS GNL AND WSG OF THE NON-LOCAL PSEUDOPOTENTIAL==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: gmax(*), fint(*), wsgtmp(*)

    INTEGER                                  :: i, ierr, ir, isub, iv, k, l, &
                                                lm, m, nghl, lmax
    REAL(real_8)                             :: rcm2, rcmgh, rgh2, &
                                                vect(nhx,lmaxx)

    IF (sgpp1%tsgp(is).OR.dpot_mod%lloc(is).LE.0) RETURN
    CALL tiset('    RNLSET',isub)
    CALL dcopy(atom_common%meshvp(is),vr(1,is,dpot_mod%lloc(is)),1,gmax(1),1)
    IF (dpot_mod%tkb(is)) THEN
       ! Kleinman-Bylander scheme
       IF (atom_common%meshw(is).NE.atom_common%meshvp(is))&
            CALL stopgm('RNLSET','INCOMPATIBLE MESHES',& 
            __LINE__,__FILE__)
       IF (atom_common%clogw(is).NE.atom_common%clogvp(is))&
            CALL stopgm('RNLSET','INCOMPATIBLE MESHES',& 
            __LINE__,__FILE__)
       IF (dcacp%active) THEN
          lmax = 4
       ELSE
          lmax = dpot_mod%lmax(is)
       ENDIF
       DO l=1,lmax
          IF (l == dpot_mod%lskip(is) .or. l == dpot_mod%lloc(is) .or. &
              dcacp%skip(l,is) ) THEN
             wsgtmp(l)=0._real_8
          ELSE
             CALL dcopy(atom_common%meshvp(is),vr(1,is,l),1,gnl(1,is,l),1)
             DO ir=1,atom_common%meshw(is)
                gnl(ir,is,l)=gnl(ir,is,l)-gmax(ir)
                IF (ABS(gnl(ir,is,l)).LT.1.e-10_real_8) gnl(ir,is,l) = 0._real_8
             ENDDO
             DO ir=1,atom_common%meshw(is)
                fint(ir)=rps(ir,is,l)**2*gnl(ir,is,l)*rw(ir,is)
             ENDDO
             CALL simpsn(atom_common%meshw(is),fint,wsg(is,l))
             wsgtmp(l)=1._real_8/(atom_common%clogw(is)*wsg(is,l))
          ENDIF
       ENDDO
       ! ..aa
       k=0
       DO l=1,lmax
          IF ( l /= dpot_mod%lloc(is) .AND. l /= dpot_mod%lskip(is) .AND. &
               .not. dcacp%skip(l,is) ) THEN
             DO lm=1,2*l-1
                k=k+1
                wsg(is,k)=wsgtmp(l)
             ENDDO
          ENDIF
       ENDDO
       ! ..aa
    ELSE
       ! Gauss-Hermit integration
       ! scaling parameter for integration points
       rcmgh=0.72_real_8
       rcm2=rcmgh*rcmgh
       nghl=NINT(nlps_com%rmaxn(is))
       DO iv=1,nghl
          rgh(iv,is)=rgh(iv,is)*rcmgh
       ENDDO
       DO l=1,dpot_mod%lmax(is)
          CALL dcopy(atom_common%meshvp(is),vr(1,is,l),1,gnl(1,is,l),1)
          DO ir=1,atom_common%meshvp(is)
             gnl(ir,is,l)=gnl(ir,is,l)-gmax(ir)
             IF (ABS(gnl(ir,is,l)).LT.1.e-10_real_8) gnl(ir,is,l) = 0._real_8
          ENDDO
          CALL curv1(atom_common%meshvp(is),rv(1,is),gnl(1,is,l),0.0_real_8,0.0_real_8,3,&
               wsgtmp,fint,0.0_real_8,ierr)
          DO i=1,nghl
             vect(i,l)=curv2(rgh(i,is),atom_common%meshvp(is),rv(1,is),gnl(1,is,l),&
                  wsgtmp,0.0_real_8)
          ENDDO
       ENDDO
       DO iv=1,nlps_com%ngh(is)
          m=MOD(iv-1,nghl)+1
          l=nghtol(iv,is)+1
          rgh2=rgh(m,is)**2
          wsg(is,iv)=rcmgh*wgh(m,is)*rgh2*vect(m,l)*EXP(rgh2/rcm2)
       ENDDO
    ENDIF
    CALL tihalt('    RNLSET',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlset
  ! ==================================================================

END MODULE rnlset_utils
