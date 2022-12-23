MODULE velupa_utils
  USE harm,                            ONLY: dtan2c,&
                                             dtan2w,&
                                             freq,&
                                             xmu
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_elec,&
                                             dtb2me

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: velupa

CONTAINS

  ! ==================================================================
  SUBROUTINE velupa(c0,cm,c2,nstate,nnlst)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                cm(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
    INTEGER                                  :: nstate, nnlst

    INTEGER                                  :: i, ig, isub
    REAL(real_8)                             :: dtx, hfo, hgi

    CALL tiset('    VELUPA',isub)
    ! WAVEFUNCTION VELOCITY UPDATE
    IF (cntl%tharm) THEN
       !$omp parallel do private(I,IG,HFO,HGI)
       DO i=1,nstate
          DO ig=1,ncpw%ngw
             IF (cntl%tmass) THEN
                hfo=-xmu(ig)*freq(ig)**2
             ELSE
                hfo=-cntr%emass*freq(ig)**2
             ENDIF
             IF (nnlst.EQ.1) THEN
                hgi=dtan2w(ig)
                cm(ig,i)=cm(ig,i)+hgi*(c2(ig,i)-hfo*c0(ig,i))
             ELSE
                hgi=dtan2c(ig)
                cm(ig,i)=cm(ig,i)+hgi*c2(ig,i)
             ENDIF
          ENDDO
       ENDDO
    ELSE
       IF (cntl%tmass) THEN
          !$omp parallel do private(I,IG,HGI)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                hgi=REAL(nnlst,kind=real_8)*0.5_real_8*dt_elec/xmu(ig)
                cm(ig,i)=cm(ig,i)+hgi*c2(ig,i)
             ENDDO
          ENDDO
       ELSE
          dtx=REAL(nnlst,kind=real_8)*dtb2me
          CALL daxpy(2*nstate*nkpt%ngwk,dtx,c2(1,1),1,cm(1,1),1)
       ENDIF
    ENDIF
    CALL tihalt('    VELUPA',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupa
  ! ==================================================================

END MODULE velupa_utils
