MODULE inr_dr_utils
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cotr,                            ONLY: cotc0,&
                                             hess
  USE do_perturbation_p_utils,         ONLY: give_scr_perturbation
  USE eicalc_utils,                    ONLY: eicalc
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gnodim
  USE hess_eta_p_utils,                ONLY: hess_eta_p
  USE hessup_utils,                    ONLY: hessup
  USE implhv,                          ONLY: rs_v,&
                                             sd0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE lanc_phon_p_utils,               ONLY: no_mode,&
                                             setsd0
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nlps,                            ONLY: imagp,&
                                             ndfnl
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: lag_mult
  USE phonons_p_utils,                 ONLY: give_scr_phonon,&
                                             settras
  USE puttau_utils,                    ONLY: gettau,&
                                             puttau
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00,&
                                             rho0,&
                                             vofrho0
  USE rgdiis_utils,                    ONLY: rgdiis
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_p_utils,                   ONLY: rnlsm3
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: infi
  USE sdion_utils,                     ONLY: sdion
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             fnl
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tiset
  USE utils,                           ONLY: dspevy
  USE vofrho_utils,                    ONLY: vofrho
  USE wrgeo_utils,                     ONLY: wrgeof_inr
  USE xinr,                            ONLY: &
       correct, direct, gnx_inr, hess_1, inr_integer, inr_logical, rmixsd, &
       tol_inr, tolx_inr, zprec
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: inr_dr
  PUBLIC :: inr
  !public :: inrsolve
  !public :: hessininr
  PUBLIC :: give_scr_inr
  !public :: cont_inr
  !public :: delta_rcm

CONTAINS

  ! =====================================================================
  SUBROUTINE inr_dr(c0,c1,nstate,xpar,dxpar,psi,rhoe,&
       tras,tmp,eirop,eivps,told,gold)
    ! ==-----------------------------------------------------------------==
    ! ==  Calls INR routine or SDION/RGDIIS, if mixed scheme is adopted  ==
    ! ==  SOLVE THE EXACT QUADRATIC EQUATION AND DETERMINE AN IONIC STEP ==
    ! ==-----------------------------------------------------------------==
    ! == INPUT:                                                          ==
    ! ==   XPAR(NODIM)                                                   ==
    ! ==   DXPAR(NODIM)                                                  ==
    ! ==   TOL_INR                                                       ==
    ! == OUTPUT:                                                         ==
    ! ==   XPAR(NODIM)                                                   ==
    ! ==   (NODIM Number of degrees of freedom                           ==
    ! ==-----------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: xpar(*), dxpar(*)
    COMPLEX(real_8)                          :: psi(maxfftn)
    REAL(real_8)                             :: rhoe(:,:), &
                                                tras(3*ions1%nat,*), tmp(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    REAL(real_8)                             :: told(*), gold(*)

    INTEGER                                  :: i
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL, SAVE                            :: done

! ==--------------------------------------------------------------==
! ==  SOLVE HESS*S = -FION FOR HESS POSITIVE DIFINITE             ==
! ==--------------------------------------------------------------==
! CALL TISET('    INR_DR',ISUB)
! ==--------------------------------------------------------------==

    IF (ifirst.EQ.0) THEN
       ifirst=1
       IF ((infi.EQ.1).AND.(inr_logical%tmixsd.OR.inr_logical%tmixgdiis)) THEN
          done=.FALSE.
       ELSE
          done=.TRUE.
       ENDIF
    ENDIF

    IF (.NOT.done) THEN
       IF (paral%parent) THEN
          IF (inr_logical%tmixsd) CALL sdion(xpar,dxpar)
          IF (inr_logical%tmixgdiis) CALL rgdiis(xpar,dxpar,told,gold)
          CALL gnodim(dxpar,cotc0%nodim,gnmax,gnorm)

          IF (gnmax.LT.rmixsd) THEN
             done=.TRUE.
             ! ... changing tmixgdiis for hessup in rgmopt
             inr_logical%tmixgdiis=.FALSE.
             IF (paral%io_parent)&
                  WRITE(6,'(/,4x,a,/)') 'THE MIXED PROCEDURE HAS REACHED'&
                  //' THE PRESCRIBED THRESHOLD'
             IF (paral%io_parent)&
                  WRITE(6,'(4x,a,/)') '******* STARTING NOW INR *******'
             IF (paral%io_parent)&
                  WRITE(6,*) ' '
          ENDIF
       ENDIF
       CALL mp_bcast(done,parai%source,parai%allgrp)
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) GOTO 541
    ENDIF

    ! ... at this point "done" can be true because we are not using a mixed 
    ! ... optimization OR we reached the condition to start INR. The double 
    ! ... "if" statement avoids a second, useless, c0 optimization 

    IF (done) THEN
       IF (paral%parent) THEN
          CALL gnodim(dxpar,cotc0%nodim,gnmax,gnorm)
          DO i=1,inr_integer%nreg
             IF (gnmax.GT.gnx_inr(i)) THEN
                tol_inr=tolx_inr(i)
                GOTO 312
             ENDIF
          ENDDO
312       CONTINUE
       ENDIF
       CALL mp_bcast(tol_inr,parai%source,parai%allgrp)

       CALL inr(c0,c1,nstate,xpar,dxpar,psi,rhoe,&
            tras,tmp,eirop,eivps)
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) GOTO 541
    ENDIF

541 CONTINUE

    ! CALL TIHALT('    INR_DR',ISUB)

    RETURN
  END SUBROUTINE inr_dr


  ! *********************************************************************

  ! =====================================================================
  SUBROUTINE inr(c0,c1,nstate,xpar,dxpar,psi,rhoe,&
       tras,tmp,eirop,eivps)
    ! ==-----------------------------------------------------------------==
    ! ==  SOLVE THE EXACT QUADRATIC EQUATION AND DETERMINE AN IONIC STEP ==
    ! ==-----------------------------------------------------------------==
    ! == INPUT:                                                          ==
    ! ==   XPAR(NODIM)                                                    ==
    ! ==   DXPAR(NODIM)                                                  ==
    ! ==   TOL_INR                                                       ==
    ! == OUTPUT:                                                         ==
    ! ==   XPAR(NODIM)                                                   ==
    ! ==   (NODIM Number of Freedom degres                               ==
    ! ==-----------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: xpar(*), dxpar(*)
    COMPLEX(real_8)                          :: psi(maxfftn,1)
    REAL(real_8)                             :: rhoe(:,:), &
                                                tras(3*ions1%nat,*), tmp(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'inr'

    COMPLEX(real_8), ALLOCATABLE             :: pippo(:)
    INTEGER                                  :: dummy, i, ierr, iter, ldfnl, &
                                                lfnl
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: occurred
    REAL(real_8)                             :: ak, ak1, ak2, bk, bk1, bk2, &
                                                deq_acc, dequad, dnrm2, dx, &
                                                dxmax, er1, gnorm, norm
    REAL(real_8), ALLOCATABLE                :: drhoe(:), dypar(:), &
                                                fscr(:,:,:), tscr(:,:,:), &
                                                ypar(:), z11(:,:)
    REAL(real_8), EXTERNAL                   :: ddot

! ... HESSUP

    IF (paral%parent)  THEN
       ALLOCATE(ypar(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent)  THEN
       ALLOCATE(dypar(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(pippo(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(z11(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(drhoe(2*maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! ... init
    occurred=.FALSE.
    ! sconvgd=.false.
    dxmax=0.2_real_8

    ! ==--------------------------------------------------------------==
    ! ==  SOLVE HESS*S = -FION FOR HESS POSITIVE DIFINITE             ==
    ! ==--------------------------------------------------------------==
    ! CALL TISET('      INR',ISUB)

    IF (inr_logical%inr_prec.AND.paral%parent) THEN
       ALLOCATE(hess_1(cotc0%nodim,cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL hessininr
    ENDIF

    ! the potential v[n0]:
    CALL rhoofr(c0,rho0,psi(:,1),nstate)
    IF (cntl%ttau) CALL stopgm("inr_dr","no tau functionals",& 
         __LINE__,__FILE__)
    CALL dcopy(fpar%nnr1, rho0, 1, vofrho0, 1)
    CALL vofrho(tau0,fion,vofrho0,psi,.FALSE.,.FALSE.)

    ! calculate the lagrange parameters z11
    CALL lag_mult(c0,c1,psi,rhoe,z11,nstate)
    CALL eicalc(eivps,eirop)

    CALL rnlsm3(c0,nstate,ddfnl)

    ! ==--------------------------------------------------------------==
    ! FNL,DFNL,DDFNL...
    ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    lfnl = imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt

    CALL rnlsm(c0,nstate,1,1,.TRUE.)
    CALL dcopy(lfnl,   fnl,1,  fnl00,1)
    CALL dcopy(ldfnl, dfnl,1, dfnl00,1)

    CALL zeroing(sd0)!,3*3*ions1%nat*ions1%nat)
    CALL setsd0(eirop,ddfnl,nstate,psi(:,1))
    CALL zeroing(rs_v)!,2*3*ions1%nat)

    IF (paral%parent) CALL settras(tras)

    iter=0
    deq_acc=0._real_8
    IF (paral%parent) THEN
       IF (inr_logical%inr_cont.AND.(ifirst.EQ.0)) THEN
          CALL cont_inr(cotc0%nodim,iter,bk1,tol_inr,deq_acc,&
               correct,xpar,'READ')
          IF (paral%io_parent)&
               WRITE(6,'(/,2x,a,3x,i6,/)') 'RESTARTING  A PREVIOUS JOB --'&
               //' WITH ITERATION',iter+1
          ifirst=1
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,2x,a,e10.4)') 'STARTING  INR:  '//&
            'RESIDUAL THRESHOLD FOR CONVERGENCE =',tol_inr
    ENDIF

    CALL mp_bcast(ifirst,parai%source,parai%allgrp)
    CALL mp_bcast(iter,parai%source,parai%allgrp)
    CALL mp_bcast(tol_inr,parai%source,parai%allgrp)
    CALL mp_bcast(correct,cotc0%nodim,parai%source,parai%allgrp)


199 FORMAT(18x,'RESIDUAL AT ',i4,3x,'ITERATION ',6x,'=',e10.4)
100 IF (iter.LT.inr_integer%itmax_inr) THEN

       iter=iter+1

       IF (iter.EQ.1) THEN

          ! .. r(0)=gradient
          IF (paral%parent) THEN

             CALL dcopy(cotc0%nodim,dxpar,1,correct,1)
             CALL dscal(cotc0%nodim,-1._real_8,correct,1)
             CALL no_mode(tras,correct,ions1%nat)

             ! ..   compute z(0)=A(cntl%prec)^(-1)*r(0)
             IF (inr_logical%inr_prec) THEN
                CALL inrsolve(correct,zprec)
             ELSE
                CALL dcopy(cotc0%nodim,correct,1,zprec,1)
             ENDIF
             CALL no_mode(tras,zprec,ions1%nat)

             CALL dcopy(cotc0%nodim,zprec,1,direct,1)
             bk1=ddot(cotc0%nodim,correct,1,zprec,1)

             CALL gnodim(dxpar,cotc0%nodim,er1,gnorm)

             IF (paral%io_parent)&
                  WRITE(6,199)iter,er1
          ENDIF

          CALL mp_bcast(direct,cotc0%nodim,parai%source,parai%allgrp)
       ELSE
          IF (paral%parent) THEN
             ! .. compute z(0)=A(cntl%prec)^(-1)*r(0)
             IF (inr_logical%inr_prec) THEN
                CALL inrsolve(correct,zprec)
             ELSE
                CALL dcopy(cotc0%nodim,correct,1,zprec,1)
             ENDIF
             CALL no_mode(tras,zprec,ions1%nat)

             ! ..   calculation of the coefficients necessary for the linear combination
             bk2=bk1

             bk1=ddot(cotc0%nodim,correct,1,zprec,1)

             bk=bk1/bk2
             ! ..   Z_k=Z_k+bk*p_(k-1)
             CALL daxpy(cotc0%nodim,bk,direct,1,zprec,1)
             ! ..   P_K=Z=K
             CALL dcopy(cotc0%nodim,zprec,1,direct,1)
          ENDIF

          CALL mp_bcast(direct,cotc0%nodim,parai%source,parai%allgrp)

       ENDIF

       ! ...  compute Hp(iter) and store it in z 
       CALL dcopy(cotc0%nodim,direct,1,rs_v(1,1),1)

       ! ... normalizing to 1 the perturbation vector for numerical reasons
       norm=dnrm2(3*ions1%nat,rs_v(1,1),1)
       CALL dscal(3*ions1%nat, 1._real_8/norm,rs_v(1,1),1)
       CALL dcopy(cotc0%nodim,rs_v(1,1),1,tmp,1)

       ! ==--------------------------------------------------------------==
       ! ... H.p product
       ! ==--------------------------------------------------------------==
       CALL hess_eta_p(c0,c1,psi(:,1),rhoe,drhoe,&
            eirop,eivps,z11,nstate,1,dummy,ddfnl,tmp)
       ! ==--------------------------------------------------------------==

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)

       IF (soft_com%exsoft) THEN
          CALL puttau(tau0,xpar)
          GOTO 214
       ENDIF

       ! ... restoring the original norm in the H.p vector
       CALL dscal(3*ions1%nat,norm,rs_v(1,2),1)
       CALL mp_sum(rs_v(:,2),cotc0%nodim,parai%allgrp)

       IF (paral%parent) THEN
          CALL no_mode(tras,rs_v(1,2),ions1%nat)

          IF (inr_logical%inr_step) THEN
             ak1=ddot(cotc0%nodim,rs_v(1,2),1,correct,1)
             ak2=(dnrm2(cotc0%nodim,rs_v(1,2),1))**2
          ELSE
             ak2=ddot(cotc0%nodim,direct,1,rs_v(1,2),1)
             ak1=bk1
          ENDIF

          ak=ak1/ak2
          IF (ak.LT.0._real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' '
             IF (paral%io_parent)&
                  WRITE(6,559) iter,ak
             IF (paral%io_parent)&
                  WRITE(6,*) ' !!!!!!!!!  WARNING WARNING WARNING  !!!!!!!'
             IF (paral%io_parent)&
                  WRITE(6,*) ' ALPHA(K) IS NEGATIVE: BREAKING INR CYCLE'
             IF (paral%io_parent)&
                  WRITE(6,*) ' '
             ! ... ak is set to tolng, minimal displacement
             ak=cntr%tolng
             occurred=.TRUE.
          ENDIF
559       FORMAT(1x,' ALPHA(',i3,') =',f10.5)

          CALL gettau(tscr,xpar)
          CALL gettau(fscr,direct)
          CALL wrgeof_inr(tscr,fscr)
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("*"),/)')

          ! ccc  .... Hessian....
          CALL dcopy (cotc0%nodim,xpar,1,ypar,1)
          CALL dcopy (cotc0%nodim,correct,1,dypar,1)

          ! ... x_k=x_(k-1)+alfa_k.p_k
          CALL daxpy(cotc0%nodim,ak,direct,1,xpar,1)

          ! ... change in quadratic energy
          dequad = ddot(cotc0%nodim,direct,1,correct,1)
          dequad = -dequad*ak+ak*ak*ak2*0.5_real_8
          deq_acc=deq_acc+dequad

          ! ... r_k=r_(k-1)-alfa_k.A.p_k
          CALL daxpy(cotc0%nodim,-ak,rs_v(1,2),1,correct,1)

          ! ... HESSIAN UPDATE
          CALL hessup(xpar,ypar,dypar,correct)

          dxmax=0.0_real_8
          ! Check convergence on atomic positions (for the next it.)
          ! (DXMAX .LT. 0.01*TOLNG)
          DO i=1,cotc0%nodim
             dx=ABS(xpar(i)-ypar(i))
             dxmax=MAX(dx,dxmax)
          ENDDO
          IF (dxmax.LT.1.e-2_real_8*cntr%tolng) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1x,a,/)') 'NEGLIGIBLE DISPLACEMENT;'&
                  //' BREAKING INR CYCLE'
             occurred=.TRUE.
          ENDIF

          CALL gnodim(correct,cotc0%nodim,er1,gnorm)

          IF (paral%io_parent)&
               WRITE(6,'(1X,64("*"))')
          IF (paral%io_parent)&
               WRITE(6,561) iter,iter,ak
          IF (paral%io_parent)&
               WRITE(6,198) er1,dequad,dxmax
          IF (inr_logical%inr_verbose) CALL delta_rcm(direct,cotc0%nodim)

          IF (paral%io_parent)&
               WRITE(6,'(1X,64("*"),/)')
198       FORMAT(1x,'RESIDUAL=',2x,e10.4,3x,'DEQUAD.=',1x,e10.4,4x,&
               'DXMAX=',1x,e10.4)
561       FORMAT(1x,'ITERATION',3x,i3,4x,' ALPHA(',i3,') =',f9.5)

       ENDIF

       CALL mp_bcast(er1,parai%source,parai%allgrp)
       CALL mp_bcast(occurred,parai%source,parai%allgrp)

       IF (paral%parent) THEN
          CALL gettau(tscr,xpar)
          CALL geofile(tscr,fscr,'WRITE')
          CALL cont_inr(3*ions1%nat,iter,bk1,tol_inr,deq_acc&
               ,correct,xpar,'WRITE')
       ENDIF

       IF (occurred) GOTO 214

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) THEN
          CALL puttau(tau0,xpar)
          GOTO 214
       ENDIF
       IF (er1.GT.tol_inr) GOTO 100
    ENDIF

214 CONTINUE

    IF (paral%parent) THEN
       IF (er1.LE.tol_inr) THEN
          IF (paral%io_parent)&
               WRITE(6,215) deq_acc
       ELSE
          IF (paral%io_parent)&
               WRITE(6,216) deq_acc
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"),/)')
       IF (inr_logical%inr_prec) DEALLOCATE(hess_1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
215 FORMAT(1x,'INR CONVERGED:        CYCLE CHANGE IN ENERGY',3x,e10.4,&
         2x,'A.U.')
216 FORMAT(1x,'INR DID NOT CONVERGE: CYCLE CHANGE IN ENERGY',3x,e10.4,&
         2x,'A.U.')

    IF (paral%parent) DEALLOCATE(ypar,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) DEALLOCATE(dypar,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pippo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! CALL TIHALT('      INR',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE inr
  ! ==================================================================


  ! **************************** SUBROUTINE **********************************



  ! ==================================================================
  SUBROUTINE inrsolve(x,y)
    ! ==================================================================
    ! == Y=A^(-1)*X                                                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x(cotc0%nodim), y(cotc0%nodim)

    CALL zeroing(y)!,cotc0%nodim)
    CALL dgemv ('N',cotc0%nodim,cotc0%nodim,1._real_8,hess_1,cotc0%nodim,x,1,0._real_8,y,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE inrsolve
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hessininr
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'hessininr'

    INTEGER                                  :: i, ierr, j, k
    REAL(real_8) :: vect(cotc0%nodim,cotc0%nodim), &
      vect_1(cotc0%nodim,cotc0%nodim), w(cotc0%nodim), work(3*cotc0%nodim), &
      wwd, wwds, wwu, wwus
    REAL(real_8), ALLOCATABLE                :: hess_lpack(:)

! variables

    IF (paral%io_parent)&
         WRITE(6,'(2x,a)') 'INR: CONSTRUCTION OF THE PRECONDITIONER'
    ! allocations
    ALLOCATE(hess_lpack(cotc0%nodim*(cotc0%nodim+1)/2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! ... parameters
    ! ...1.5e-1_real_8 a.u. ~ 65 cm^-1

    wwd=1.e-2_real_8
    wwds=0.9_real_8
    wwu=1.05_real_8
    wwus=1.05_real_8

    k = 0
    DO  i=1,cotc0%nodim
       DO  j=i,cotc0%nodim
          k = k+1
          hess_lpack(k)=hess(j,i)
       ENDDO
    ENDDO
    CALL dspevy(1,hess_lpack,w,vect,cotc0%nodim,cotc0%nodim,work,3*cotc0%nodim)

    CALL zeroing(vect_1)!,cotc0%nodim*cotc0%nodim)

    DO i=1,cotc0%nodim
       IF (w(i).LE.wwd) THEN

          CALL dcopy(cotc0%nodim,vect(1,i),1,vect_1(1,i),1)
          IF ((inr_logical%inr_verbose).AND.paral%io_parent)&
               WRITE(6,890)i,w(i),&
               SIGN(5140.487_real_8*SQRT(ABS(w(i))),w(i)),wwds
          CALL dscal(cotc0%nodim,1/wwds,vect_1(1,i),1)

       ELSE IF (w(i).GT.wwu) THEN

          CALL dcopy(cotc0%nodim,vect(1,i),1,vect_1(1,i),1)
          CALL dscal(cotc0%nodim,1/wwus,vect_1(1,i),1)
          IF ((inr_logical%inr_verbose).AND.paral%io_parent)&
               WRITE(6,890)i,w(i),&
               SIGN(5140.487_real_8*SQRT(ABS(w(i))),w(i)),wwus

       ELSE

          IF ((inr_logical%inr_verbose).AND.paral%io_parent)&
               WRITE(6,889)i,w(i),&
               SIGN(5140.487_real_8*SQRT(ABS(w(i))),w(i))
          CALL dcopy(cotc0%nodim,vect(1,i),1,vect_1(1,i),1)
          CALL dscal(cotc0%nodim,1/w(i),vect_1(1,i),1)

       ENDIF
    ENDDO

889 FORMAT(2x,'EIGENV. ', i3,' (', e11.4,')','(', f12.3,')',' SCALED')
890 FORMAT(2x,'EIGENV. ', i3,' (', e11.4,')','(', f12.3,')',&
         ' SCALED BY', f10.6)
    CALL dgemm('N','T',cotc0%nodim,cotc0%nodim,cotc0%nodim,1._real_8,vect,&
         cotc0%nodim,vect_1,cotc0%nodim,0._real_8,hess_1,cotc0%nodim)


    DEALLOCATE(hess_lpack,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hessininr
  ! ==================================================================



  ! ==================================================================
  SUBROUTINE give_scr_inr(lrinr,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrinr
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lphonon

! NODIM MCNSTR
! ==--------------------------------------------------------------==

    lrinr=cotc0%nodim*4
    CALL give_scr_phonon(lphonon,tag,nstate)
    lrinr=MAX(lphonon,lrinr)
    CALL give_scr_perturbation(lphonon,tag,nstate)
    lrinr=MAX(lphonon,lrinr)
    tag ='NODIM'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_inr
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cont_inr(ndim,iter,bk1,tol_inr,deq_acc,&
       correct,xpar,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim, iter
    REAL(real_8)                             :: bk1, tol_inr, deq_acc, &
                                                correct(*), xpar(*)
    CHARACTER(len=*)                         :: tag

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: iunit, k
    LOGICAL                                  :: ferror

    iunit=21

    IF (INDEX(tag,'READ').NE.0) THEN
       filen='INR_CONTINUE'
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_old,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' CONT_INR! NO FILE INR_CONTINUE WAS FOUND'
          CALL stopgm('CONT_INR','CANNOT CONTINUE PREVIOUS JOB',& 
               __LINE__,__FILE__)
       ENDIF

       IF (paral%io_parent)&
            REWIND(iunit)
       IF (paral%io_parent)&
            READ(iunit,'(i12,3f15.10)') iter,bk1,tol_inr,deq_acc
       DO k=1,ndim
          IF (paral%io_parent)&
               READ(iunit,*) correct(k),xpar(k)
       ENDDO
    ELSEIF (INDEX(tag,'WRITE').NE.0) THEN
       filen='INR_CONTINUE.1'
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_def,ferror)

       IF (paral%io_parent)&
            REWIND(iunit)
       IF (paral%io_parent)&
            WRITE(iunit,'(i12,3f15.10)') iter,bk1,tol_inr,deq_acc
       DO k=1,ndim
          IF (paral%io_parent)&
               WRITE(iunit,*) correct(k),xpar(k)
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' CONT_INR| UNKNOWN TAG: ',tag
    ENDIF

    IF (paral%io_parent)&
         CALL fileclose(iunit)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cont_inr
  ! ==================================================================


  ! ==================================================================

  ! ==--------------------------------------------------------------==
  SUBROUTINE delta_rcm(mode,nrighe)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: mode(*)
    INTEGER                                  :: nrighe

    INTEGER                                  :: i, j, k, l
    REAL(real_8)                             :: norm, SUM(3)

! ... arguments
! ... internal variables

    norm=0._real_8
    DO j=1,ions1%nsp
       norm=norm+ions0%na(j)*rmass%pma0(j)**2
    ENDDO
    l=0
    DO k=1,3
       SUM(k)=0._real_8
    ENDDO
    DO j=1,ions1%nsp
       DO i=1,ions0%na(j)
          DO k=1,3
             l=l+1
             SUM(k)=SUM(k)+mode(l)*rmass%pma0(j)/norm
             ! sum=sum+mode(l)*SQRT(masse(i))
          ENDDO
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,198) (SUM(k),k=1,3)
198 FORMAT(1x,'COM',3x,'DX=  ',e10.4,8x,'DY= ',e10.4,7x,'DZ= ',e10.4)

    RETURN
  END SUBROUTINE delta_rcm

  ! ==--------------------------------------------------------------==



END MODULE inr_dr_utils
