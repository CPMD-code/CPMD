MODULE opt_lr_utils
  USE cppt,                            ONLY: indz,&
                                             nzh,&
                                             scg
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: kr1m
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_ortho
  USE lr_upd_utils,                    ONLY: give_scr_lr_upd,&
                                             lr_upd
  USE machine,                         ONLY: m_walltime
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: ropt_mod
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE utils,                           ONLY: zclean
  USE vpsi_utils,                      ONLY: vpsimt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: opt_lr
  PUBLIC :: give_scr_opt_lr
  !public :: lr_print
  !public :: lr_fpri

CONTAINS

  ! ==================================================================
  SUBROUTINE opt_lr(c0,c1,c2,sc0,eind,eigv,drhoe,h1nl,eivps1,eirop1,&
       ddxc,vpp,pme,gde,psi,nstate,TYPE,orbital)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eind, eigv(:), drhoe(:,:)
    COMPLEX(real_8)                          :: eivps1(*), eirop1(*)
    REAL(real_8)                             :: ddxc(:,:), vpp(:)
    COMPLEX(real_8)                          :: pme(*), gde(*), psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: h1nl(nkpt%ngwk,nstate), sc0(nkpt%ngwk,nstate), &
      c2(nkpt%ngwk,nstate), c1(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate)
    CHARACTER(len=*)                         :: TYPE, orbital

    CHARACTER(*), PARAMETER                  :: procedureN = 'opt_lr'

    COMPLEX(real_8), ALLOCATABLE             :: v1(:)
    INTEGER                                  :: ierr, ig, ir, is, isub, knfi
    LOGICAL                                  :: debug
    REAL(real_8)                             :: detot, e2(5), etot0, etot2, &
                                                fi, tcpu, time1, time2
    REAL(real_8), ALLOCATABLE                :: vpot(:,:)

! VARIABLES
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    CALL setfftn(0)
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    ! ==--------------------------------------------------------------==
    ! 
    debug=.FALSE.
    etot2 = 0._real_8
    ! ==--------------------------------------------------------------==
    ! == CONTINUATION: STEP USING PREVIOUS FORCES INSTEAD OF INIT     ==
    ! ==--------------------------------------------------------------==
    IF (INDEX(TYPE,"NOINIT").NE.0) THEN
       TYPE="PHONON"
       ! 
       ropt_mod%convwf=.FALSE.
       CALL lr_upd(c0,c1,c2,sc0,e2,ddxc,psi,eigv,drhoe,h1nl,vpp,pme,&
            gde,nstate,TYPE,orbital,0,.FALSE.)
       GOTO 50
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    IF (INDEX(TYPE,"PHONON").NE.0) THEN
       ALLOCATE(v1(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(vpot(fpar%nnr1,MAX(fpar%krx,kr1m*group%nogrp)*fpar%kr2s*fpar%kr3s*clsd%nlsd/fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)! TODO CHECK
       ! ALLOCATE(VPOT(NNR1,NLSD))
       DO ig=1,ncpw%nhg
          v1(ig)   = eivps1(ig) + scg(ig) * eirop1(ig)
       ENDDO
       CALL zeroing(psi(:,1))!,maxfft)
       DO ig=1,ncpw%nhg
          psi(indz(ig),1) = CONJG(v1(ig))
          psi(nzh(ig),1)  = v1(ig)
       ENDDO
       IF (geq0) psi(nzh(1),1) = v1(1)
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
       DO ir=1,fpar%nnr1
          vpot(ir,1)=REAL(psi(ir,1))
       ENDDO
       IF (cntl%tlsd.OR.lspin2%tlse) CALL dcopy(fpar%nnr1,vpot(1,1),1,vpot(1,2),1)
       IF (lspin2%tlse) CALL dcopy(fpar%nnr1,vpot(1,1),1,vpot(1,3),1)
       CALL vpsimt(c0,h1nl,crge%f(:,1),vpot,psi(:,1),nstate,clsd%nlsd,.FALSE.)
       ! ..H1NL now contains the total constant part of the gradient
       ! ..H1NL = dVPS + dVNL + dVH
       DEALLOCATE(v1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(vpot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ..initialize C1
    DO is=1,nstate
       IF (crge%f(is,1).EQ.0.0_real_8) THEN
          fi = 1.0_real_8
       ELSE
          fi = 1.0_real_8/crge%f(is,1)
       ENDIF
       DO ig=1,ncpw%ngw
          c1(ig,is)=dt2bye*vpp(ig)*h1nl(ig,is)*fi
       ENDDO
    ENDDO
    IF (geq0) CALL zclean(c1,nstate,ncpw%ngw)
    CALL lr_ortho(nstate,c0,c1)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    etot0=0._real_8
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F8.2,A8)')&
            ' TIME FOR INITIALIZATION OF LINEAR RESPONSE:',&
            tcpu,' SECONDS'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==      THE BASIC LOOP FOR WAVEFUNCTION OPTIMIZATION            ==
    ! ==--------------------------------------------------------------==
50  ropt_mod%convwf=.FALSE.
    DO knfi=1,lr01%nmaxlr
       time1=m_walltime()
       CALL zeroing(e2)!,5)
       e2(5)=eind
       CALL lr_upd(c0,c1,c2,sc0,e2,ddxc,psi,eigv,drhoe,h1nl,vpp,pme,&
            gde,nstate,TYPE,orbital,knfi,.TRUE.)
       etot2=e2(1)+e2(2)+e2(3)+e2(4)+e2(5)
       IF (knfi.EQ.1) THEN
          detot=0._real_8
       ELSE
          detot=etot2-etot0
       ENDIF
       IF (paral%parent) THEN
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
       ENDIF
       etot0 = etot2
       IF (paral%parent) CALL lr_print(knfi,etot2,detot,tcpu)
       IF (ropt_mod%convwf) GOTO 100
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==     end of main loop                                         ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,T64,"!!")')&
            ' OPT_SDWFN1| the maximum number of steps is reached'
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,F10.6,A,T64,"!!")')&
            '         but no convergence (d E / d phi_1 =',gemax,')'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
100 CONTINUE
    IF (paral%parent) CALL lr_fpri(e2,etot2)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE opt_lr
  ! ==================================================================
  SUBROUTINE give_scr_opt_lr(lopt_lr,TYPE,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lopt_lr
    CHARACTER(len=*)                         :: TYPE
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: llr_ortho, llr_upd

! ==--------------------------------------------------------------==

    lopt_lr=2*ncpw%nhg + MAX(fpar%krx,kr1m*group%nogrp)*fpar%kr2s*fpar%kr3s * clsd%nlsd +100
    CALL give_scr_lr_upd(llr_upd,TYPE,tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    lopt_lr=MAX(lopt_lr,llr_upd,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_opt_lr
  ! ==================================================================
  SUBROUTINE lr_print(knfi,etot2,detot,tcpu)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: knfi
    REAL(real_8)                             :: etot2, detot, tcpu

! ==--------------------------------------------------------------==

    IF ((knfi.EQ.1).AND.paral%io_parent)&
         WRITE(6,'(A,A)')&
         '   NFI      GEMAX     CNORM',&
         '         ETOT(2)        DETOT     TCPU'
    IF (paral%io_parent)&
         WRITE(6,'(I6,F11.6,F10.6,F16.8,G13.4,F9.3)')&
         knfi,gemax,cnorm,etot2,detot,tcpu
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_print
  ! ==================================================================
  SUBROUTINE lr_fpri(e2,etot2)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: e2(5), etot2

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(/)')
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ======== 2nd Order Perturbation Energy =',etot2,' ========='
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ========              Equad(2) Hartree =',e2(1),' ========='
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ========              Equad(2)      XC =',e2(2),' ========='
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ========              Equad(2)   (H-e) =',e2(3),' ========='
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ========              Elin(2)          =',e2(4),' ========='
    IF (paral%io_parent)&
         WRITE(6,'(A,T40,F16.8,A)')&
         ' ========              Eind(2)          =',e2(5),' ========='
    IF (paral%io_parent)&
         WRITE(6,'(/)')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_fpri
  ! ==================================================================

END MODULE opt_lr_utils
