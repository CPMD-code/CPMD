MODULE hnlmat_utils
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hnlmat

CONTAINS

  ! ==================================================================
  SUBROUTINE hnlmat(hmat,f,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate), hmat(nstate,nstate)

    INTEGER                                  :: chunk, i, ia, ind, ind1, is, &
                                                isa, isa0, ispin, isub, iv, &
                                                j, jmax, jv, ki, kj, l, l2, &
                                                li, lj, mxmypos, mxrank, mypos
    REAL(real_8)                             :: dd, fdd, ffi

! Variables
! ==--------------------------------------------------------------==
! == Compute the non-local contribution to the Hamilton matrix    ==
! ==--------------------------------------------------------------==

    CALL tiset('    HNLMAT',isub)
    IF (cntl%tfdist) CALL stopgm('HNLMAT','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('HNLMAT','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO i=1,nstate
       ffi=f(i)
       IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
       DO j=1,nstate
          hmat(i,j)=hmat(i,j)/ffi
       ENDDO
    ENDDO

    ! AK   FIXME: 20080131 the HNLMAT() optimization does not work for LSD.
    ! AK   FIXME: we drop back to old algorithm.

    IF (cntl%tlsd) THEN
       ! 1.. NSTATE for serial !
       DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
          ispin=1
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
          DO j=i,jmax
             isa0=0
             DO is=1,ions1%nsp
                IF (pslo_com%tvan(is)) THEN
                   ! VANDERBILT PP
                   DO iv=1,nlps_com%ngh(is)
                      DO jv=1,nlps_com%ngh(is)
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            dd=deeq(isa,iv,jv,ispin)+dvan(iv,jv,is)
                            fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,jv,j,1)
                            hmat(i,j)=hmat(i,j)-fdd
                            IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                         ENDDO
                      ENDDO
                   ENDDO
                ELSEIF (sgpp1%tsgp(is)) THEN
                   ! Stefan Goedecker PP
                   DO iv=1,nlps_com%ngh(is)
                      l=nghtol(iv,is)+1
                      ki=sgpp2%lfval(iv,is)
                      li=sgpp2%lpval(iv,is)
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l2.EQ.l.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               fdd=sgpp2%hlsg(ki,kj,l,is)*fnl(1,isa,iv,i,1)*&
                                    fnl(1,isa,jv,j,1)
                               hmat(i,j)=hmat(i,j)-fdd
                               IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ELSE
                   ! BACHELET HAMANN SCHLUTER
                   DO iv=1,nlps_com%ngh(is)
                      dd=wsg(is,iv)
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,iv,j,1)
                         hmat(i,j)=hmat(i,j)-fdd
                         IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                      ENDDO
                   ENDDO
                   ! ==-------------------------------------------------------==
                ENDIF
                isa0=isa0+ions0%na(is)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       chunk = 0
       mxrank = 0
       IF (parap%nst12(parai%mepos,1).LE.parap%nst12(parai%mepos,2)) THEN
          ind1 = 0
          DO i=0, parai%nproc-1
             IF (parap%nst12(i,1).LE.parap%nst12(i,2)) THEN
                ind1 = ind1 + 1
                IF (i.GT.mxrank) mxrank=i
                IF (parai%mepos.EQ.i) mxmypos = ind1
             ENDIF
          ENDDO

          chunk = (nstate*(nstate+1)/2) / ind1
          IF (parai%mepos.EQ.mxrank) chunk=(nstate*(nstate+1)/2)-(ind1-1)*&
               chunk
          ind = 1
          mypos = (mxmypos-1)*chunk+1
          IF (parai%mepos.EQ.mxrank) THEN
             mypos=(ind1-1)*((nstate*(nstate+1)/2)/ind1)+1
          ENDIF
          DO i=1, nstate
             DO j=i, nstate
                IF (ind.EQ.mypos) GOTO 999
                ind = ind +1
             ENDDO
          ENDDO
999       CONTINUE
       ENDIF

       ! 1.. NSTATE for serial !
       isa0 = 0
       DO ind = 1, chunk

          ispin=1
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                ! VANDERBILT PP
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         dd=deeq(isa,iv,jv,ispin)+dvan(iv,jv,is)
                         fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,jv,j,1)
                         hmat(i,j)=hmat(i,j)-fdd
                         IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF (sgpp1%tsgp(is)) THEN
                ! Stefan Goedecker PP
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l2.EQ.l.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            fdd=sgpp2%hlsg(ki,kj,l,is)*fnl(1,isa,iv,i,1)*&
                                 fnl(1,isa,jv,j,1)
                            hmat(i,j)=hmat(i,j)-fdd
                            IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                ! BACHELET HAMANN SCHLUTER
                DO iv=1,nlps_com%ngh(is)
                   dd=wsg(is,iv)
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,iv,j,1)
                      hmat(i,j)=hmat(i,j)-fdd
                      IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                   ENDDO
                ENDDO
                ! ==-------------------------------------------------------==
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
          j = j+1
          isa0 = 0
          IF (j.GT.nstate) THEN
             j = i+1
             i = i+1
             isa0 = 0
          ENDIF
       ENDDO
       ! AK: end of .NOT.cntl%tlsd branch. FIXME.
    ENDIF

    DO i=1,nstate
       ffi=f(i)
       IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
       DO j=1,nstate
          hmat(i,j)=hmat(i,j)*ffi
       ENDDO
    ENDDO
    CALL tihalt('    HNLMAT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hnlmat
  ! ==================================================================

END MODULE hnlmat_utils
