MODULE ldos_utils
  USE cnst,                            ONLY: ry
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: ngrm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: invfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_mark,&
                                             fo_verb
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE kpts,                            ONLY: tkpts
  USE ldosmod,                         ONLY: cldos
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: fpar,&
                                             group,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ldos
  PUBLIC :: give_scr_ldos

CONTAINS

  ! ==================================================================
  SUBROUTINE ldos(c0,psi,we,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  LAYER PROJECTED DENSITY OF STATES                           ==
    ! ==--------------------------------------------------------------==
    ! ==  WARNING: IF YOU USE SPECIAL K-POINTS FOR A SPECIAL STRUCTURE==
    ! ==           YOU NEED TO SYMMETRIZE CHARGE DENSITY              ==
    ! ==           FOR THAT -> SPECIFY THE POINT GROUP                ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:,:), &
                                                psi(maxfft*group%nogrp)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: we(nstate,nkpt%nkpnt)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ldos'

    CHARACTER(len=100)                       :: fileldos
    INTEGER :: i, i1, i2, i3, ia, ib, ibb, ic, id, ierr, ifft, ig, ii, ik, &
      ikind, ikk, ikpt, is1, isub, iuldos, kbeg, kend, kinc, l, lead, leadx, &
      nkpoint, nnrx, nsta
    INTEGER, SAVE                            :: ifirst = 1
    LOGICAL                                  :: ferror, tfcal
    REAL(real_8)                             :: coef3, r1, r2
    REAL(real_8), ALLOCATABLE                :: dldos(:,:,:)

!(nkpt%ngwk,nstate,nkpt%nkpnt)

    CALL tiset('      LDOS',isub)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("<")/20X,A16/1X,64(">"))')&
            'CALCULATING LDOS'
    ENDIF
    IF (ifirst.EQ.1) THEN
       ifirst=0
       ALLOCATE(dldos(cldos%nlayer,nstate,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.GT.1) THEN
       CALL stopgm(' LDOS',' DOESNT WORK FOR NOGRP GT 1',& 
            __LINE__,__FILE__)
    ENDIF
    IF (group%nogrp.EQ.1) THEN
       lead  = fpar%kr1s*parai%ngrays
       leadx = fpar%nnr1
       ifft=1
       nnrx=fpar%nnr1
    ELSE
       lead  = fpar%kr1s*ngrm
       leadx = fpar%krx*fpar%kr2s*fpar%kr3s
       ifft=2
       nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
    ENDIF
    ! Initialize
    CALL zeroing(dldos)!,cldos%nlayer*nstate*nkpt%nkpnt)
    ! Accumulate the charge 
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,'HGKP HGKM C0')
       DO ikind = 1, nkpoint
          ikk=kpbeg(ikpt)+ikind
          ! Loop over the electronic states
          DO id=1,nstate,group%nogrp
             nsta=MIN((nstate-id+1),group%nogrp)
             tfcal=.TRUE.
             DO i=id,id+(nsta-1)
                tfcal=tfcal.AND.(crge%f(i,ikk).NE.0._real_8)
             ENDDO
             IF (tfcal) THEN
                CALL zeroing(psi)!,maxfft*group%nogrp)
                DO ib=1,nsta
                   is1=id+(ib-1)
                   ibb=(ib-1)*lead
                   !CDIR NODEP
                   DO ig=1,ncpw%ngw
                      psi(nzhs(ig)+ibb)=c0(ig,is1,ikind)
                      IF (tkpts%tkpnt) THEN
                         psi(indzs(ig)+ibb)=c0(ig+ncpw%ngw,is1,ikind)
                      ELSE
                         psi(indzs(ig)+ibb)=CONJG(c0(ig,is1,ikind))
                      ENDIF
                   ENDDO
                   IF (geq0) psi(nzhs(1)+ibb)=c0(1,is1,ikind)
                ENDDO
                ! Fourier transform the wave functions to real space.
                IF (ifft.EQ.1) THEN
                   CALL invfftn(psi,.TRUE.,parai%allgrp)
                ELSE
                   CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
                        __LINE__,__FILE__)
                ENDIF
                ! Compute the charge density from the wave functions
                ! in real space 
                IF (group%nogrp.GT.1) THEN
                   DO ib=1,nsta
                      IF (group%nolist(ib).EQ.parai%me) THEN
                         is1=id+(ib-1)
                      ENDIF
                   ENDDO
                ELSE
                   is1=id
                ENDIF
                coef3=wk(ikk)*crge%f(is1,ikk)/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
                IF (cldos%ldosn.EQ.3) THEN
                   l=0
                   DO  i3=1,fpar%kr3s
                      DO i2=1,fpar%kr2s
                         DO i1=1,fpar%kr1
                            l=l+1
                            ii=INT(i3*(cldos%nlayer+1)/fpar%kr3s)
                            IF (ii.LE.cldos%nlayer) THEN
                               r1=REAL(psi(l))
                               r2=AIMAG(psi(l))
                               dldos(ii,is1,ikk)= dldos(ii,is1,ikk)+&
                                    coef3*(r1*r1+r2*r2)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSEIF (cldos%ldosn.EQ.2) THEN
                   l=0
                   DO  i3=1,fpar%kr3s
                      DO i2=1,fpar%kr2s
                         DO i1=1,fpar%kr1
                            l=l+1
                            ii=INT(i2*(cldos%nlayer+1)/fpar%kr2s)
                            IF (ii.LE.cldos%nlayer) THEN
                               r1=REAL(psi(l))
                               r2=AIMAG(psi(l))
                               dldos(ii,is1,ikk)= dldos(ii,is1,ikk)+&
                                    coef3*(r1*r1+r2*r2)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ELSEIF (cldos%ldosn.EQ.1) THEN
                   l=0
                   DO  i3=1,fpar%kr3s
                      DO i2=1,fpar%kr2s
                         DO i1=1,fpar%kr1
                            ia=parap%nrxpl(parai%me,1)
                            ib=parap%nrxpl(parai%me,2)
                            ic=ib-ia
                            l=l+1
                            ii=INT((ia+i1-1)*(cldos%nlayer+1)/fpar%kr1s)
                            IF (ii.LE.cldos%nlayer.AND.ic.GE.0) THEN
                               r1=REAL(psi(l))
                               r2=AIMAG(psi(l))
                               dldos(ii,is1,ikk)= dldos(ii,is1,ikk)+&
                                    coef3*(r1*r1+r2*r2)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF         ! Endif for TFCAL
          ENDDO             ! End loop over electronic states
       ENDDO                 ! End loop over k points (IKIND)
    ENDDO                     ! End loop over IKPT
    CALL mp_sum(dldos,cldos%nlayer*nstate*nkpt%nkpnt,parai%allgrp)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       fileldos='LDOS'
       iuldos=35
       IF (paral%io_parent)&
            CALL fileopen(iuldos,fileldos,fo_app+fo_verb+fo_mark,ferror)
       IF (paral%io_parent)&
            WRITE(iuldos,*) cldos%nlayer,nstate,nkpt%nkpnt
       IF (paral%io_parent)&
            WRITE(iuldos,'(A13,1X,F10.5)') 'FERMI ENERGY:',ener_com%amu*2._real_8*ry
       DO ik=1,nkpt%nkpnt
          DO i=1,nstate
             IF (paral%io_parent)&
                  WRITE(iuldos,'(I4,I4,1X,2F10.4,1X,F10.5)')&
                  ik,i,we(i,ik)*2._real_8*ry,crge%f(i,ik),wk(ik)
             IF (paral%io_parent)&
                  WRITE(iuldos,'(8E10.4)') (dldos(ii,i,ik),ii=1,cldos%nlayer)
          ENDDO
       ENDDO
       DO i=1,cldos%nlayer
          IF (paral%io_parent)&
               WRITE(57,*)  i,(dldos(i,ii,1),ii=1,nstate)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(iuldos)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SYMMETRIZE DENSITY IF POINT GROUP SPECIFIED 
    ! (NEED FOR SPECIAL K-POINTS.
    ! deb  IF(cntl%tlsd) THEN
    ! deb    CALL SYMRHO(RHOE(1,1),PSI)
    ! deb    CALL SYMRHO(RHOE(1,2),PSI)
    ! deb  ELSE
    ! deb    CALL SYMRHO(RHOE(1,1),PSI)
    ! deb  ENDIF
    CALL tihalt('      LDOS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ldos
  ! ==================================================================
  SUBROUTINE give_scr_ldos(xxldos,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: xxldos
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    xxldos=1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ldos
  ! ==================================================================

END MODULE ldos_utils
