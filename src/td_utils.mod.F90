MODULE td_utils
  USE calc_pij_utils,                  ONLY: calc_pij
  USE cnst,                            ONLY: fpi,&
                                             pi,&
                                             uimag
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cppt,                            ONLY: gk,&
                                             gl,&
                                             hg,&
                                             indzs,&
                                             isptr,&
                                             nzh,&
                                             nzhs
  USE densrd_utils,                    ONLY: densrd
  USE dg,                              ONLY: ipooldg,&
                                             tdgcomm
  USE dipo_utils,                      ONLY: dipo
  USE dipomod,                         ONLY: moment
  USE efld,                            ONLY: extf
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: llr1
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fileopen_utils,                  ONLY: fileopen
  USE fileopenmod,                     ONLY: FO_NEW
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE hpsi_utils,                      ONLY: hpsi
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpnt,                            ONLY: rk
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE readsr_utils,                    ONLY: xstring
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhopri_utils,                    ONLY: rhopri_eh
  USE rnlrh_utils,                     ONLY: rnlrhg
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: infi
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod,&
                                             tdsp1
  USE system,                          ONLY: &
       cntl, cntr, fpar, group, mapgp, maxsys, ncpw, nkpt, parap, parm, spar
  USE td_input,                        ONLY: &
       gampl, gaugefield, genl, gndir, gpot, gpotv, itermax, maskreal, &
       pointcharge, ppcorr, td_prop
  USE utils,                           ONLY: zclean_k
  USE vofrho_utils,                    ONLY: vofrho
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pos_oper
  PUBLIC :: getnorm_k
  PUBLIC :: kick
  PUBLIC :: derivs
  PUBLIC :: gaugepot_laser
  PUBLIC :: gaugepot_laser_circ_pol
  PUBLIC :: dext_pot
  PUBLIC :: dext_pipot
  PUBLIC :: dext_ch
  PUBLIC :: dyn_analysis
  PUBLIC :: currentoper
  PUBLIC :: load_ex_states
  PUBLIC :: eh_initialize_data
  PUBLIC :: dipolox
  PUBLIC :: applymask
  PUBLIC :: tmprd_prop
  PUBLIC :: tmpwr_prop
  PUBLIC :: initialize_ehrenfest_dyn

CONTAINS

  ! ==================================================================
  SUBROUTINE scalar_product_r(phi,hphi,nstate,neigv)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: hphi(2*ncpw%ngw,nstate), &
                                                phi(2*ncpw%ngw,nstate)
    REAL(real_8)                             :: neigv(nstate)

    INTEGER                                  :: lst
    REAL(real_8)                             :: ddot

    DO lst=1,nstate
       neigv(lst)=ddot(2*ncpw%ngw*2,phi(1,lst),1,hphi(1,lst),1)
       CALL mp_sum(neigv(lst),parai%allgrp)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE scalar_product_r
  ! ==================================================================
  SUBROUTINE pos_oper(pos,dir)
    ! COMPUTES the position operator referred to a Cartesian coordinate
    ! system placed at the center of the box
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: pos(fpar%nnr1)
    INTEGER                                  :: dir

    INTEGER                                  :: ia, ir, ir1, ir2, ir3, is
    REAL(real_8)                             :: dx, dy, dz, help, r2, rmax, &
                                                rmax2, x, x0, y, y0, z, z0

    CALL zeroing(pos)!,nnr1)
    ! -- Origin --
    ! call phfac(tau0)
    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    help=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          x0=x0+tau0(1,ia,is)*ions0%zv(is)
          y0=y0+tau0(2,ia,is)*ions0%zv(is)
          z0=z0+tau0(3,ia,is)*ions0%zv(is)
          help=help+ions0%zv(is)
       ENDDO
    ENDDO
    x0=x0/help
    y0=y0/help
    z0=z0/help
    ! 
    ! -- Mesh parameters --
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)

    ! -- Position oprtator --
    CALL zeroing(pos)!,nnr1)
    rmax=0.5_real_8*parm%a1(1)
    rmax2=rmax**2
    DO ir3=1,parm%nr3
       DO ir2=1,parm%nr2
          DO ir1=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
             x=(ir1-1)*dx-x0+dx/2._real_8
             y=(ir2-1)*dy-y0+dy/2._real_8
             z=(ir3-1)*dz-z0+dz/2._real_8
             ir=(ir3-1)*parm%nr1*parm%nr2   +&
                  (ir2-1)*parm%nr1       +&
                  ir1               -&
                  parap%nrxpl(parai%mepos,1) + 1
             r2=x**2+y**2+z**2
             ! if (r2.lt.rmax2) then
             IF (dir.EQ.1) THEN
                pos(ir)=x
             ELSEIF (dir.EQ.2) THEN
                pos(ir)=y
             ELSEIF (dir.EQ.3) THEN
                pos(ir)=z
             ELSEIF (dir.EQ.4) THEN
                pos(ir)=x+y
             ELSEIF (dir.EQ.5) THEN
                pos(ir)=x+y+z
             ENDIF
             ! endif
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pos_oper
  ! ================================================================== 
  SUBROUTINE getnorm_k(phi,nstate,norms)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: phi(nkpt%ngwk,nstate)
    REAL(real_8)                             :: norms(nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'getnorm_k'

    COMPLEX(real_8)                          :: zdotc
    INTEGER                                  :: lst

    DO lst=1,nstate
       norms(lst)=REAL(zdotc(nkpt%ngwk,phi(1,lst),1,phi(1,lst),1))
    ENDDO
    CALL mp_sum(norms,nstate,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getnorm_k
  ! ==================================================================
  SUBROUTINE maxnorm_k(phil,dist)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: phil(nkpt%ngwk,crge%n,td_prop%nzeros,2)
    REAL(real_8)                             :: dist

    CHARACTER(*), PARAMETER                  :: procedureN = 'maxnorm_k'

    COMPLEX(real_8), ALLOCATABLE             :: phi1(:,:), phi2(:,:)
    INTEGER                                  :: ierr, ig, lst, msglen, &
                                                nstate, t_n
    REAL(real_8)                             :: help

    nstate=crge%n
    dist=0._real_8
    msglen= 8
    ALLOCATE(phi1(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(phi2(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    DO t_n=1,td_prop%nzeros
       ! ---------------
       DO lst=1,nstate
          DO ig=1,nkpt%ngwk
             phi1(ig,lst)=phil(ig,lst,t_n,1)
             phi2(ig,lst)=phil(ig,lst,t_n,2)
             help=REAL((phi1(ig,lst)-phi2(ig,lst)) *&
                  CONJG(phi1(ig,lst)-phi2(ig,lst)) )
             IF (help.GT.dist) dist=help
          ENDDO
       ENDDO
       CALL mp_max(dist,parai%allgrp)
    ENDDO
    ! -----
    CALL mp_bcast(dist,parai%source,parai%allgrp)
    DEALLOCATE(phi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(phi2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE maxnorm_k
  ! ==================================================================
  SUBROUTINE kick(c2,nstate,norms,pos,phi,psi)
    ! -- Momentum kick to the wavefunctions                           --
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    REAL(real_8)                             :: norms(nstate), pos(fpar%nnr1)
    COMPLEX(real_8)                          :: phi(nkpt%ngwk,nstate), psi(:)

    COMPLEX(real_8)                          :: expp, fm, fp
    INTEGER                                  :: fi, i, ia, ib, ibb, ifft, ig, &
                                                ig1, ik1, ik2, ik3, ir, lead, &
                                                leadx, lst, njump, nsta
    REAL(real_8)                             :: help

    lead  = fpar%kr1s*parai%ngrays
    leadx = fpar%nnr1
    ifft  = 1
    ! 
    IF (geq0) THEN
       ! The first component is for G=0.
       ig1=2
    ELSE
       ig1=1
    ENDIF
    ! 
    DO lst=1,nstate
       IF (geq0) THEN
          phi(1,lst)=CMPLX(REAL(c2(1,lst)),0.0_real_8,kind=real_8)
          phi(ncpw%ngw+1,lst)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDIF
       DO ig=ig1,nkpt%ngwk
          phi(ig,lst)=c2(ig,lst)
       ENDDO
    ENDDO
    ! 
    njump=group%nogrp
    ! 
    DO ia=1,nstate,njump
       ! 
       CALL zeroing(psi)!,maxfft*group%nogrp)
       nsta=MIN(nstate-ia+1,group%nogrp)
       DO ib=1,nsta
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          DO ig=1,ncpw%ngw
             ! the +g component of the state i
             psi(nzhs(ig)+ibb)=phi(ig,i)
             ! the -g component of the state i
             psi(indzs(ig)+ibb)=phi(ig+ncpw%ngw,i)
          ENDDO
          IF (geq0) psi(nzhs(1)+ibb)=phi(1,i)
       ENDDO
       ! 
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       ! ig=1
       ! do ir=1,nnr1
       ! if (psi(ir).ne.0._real_8) then
       ! psi(ir)=(1.0_real_8+uimag*ampl*pos(ig))*psi(ir)
       ! help=pos(ig)*ampl
       ! expp=cmplx(cos(help),-sin(help),kind=real_8)
       ! psi(ir)=psi(ir)*expp
       ! ig=ig+1
       ! endif
       ! enddo

       ig=0
       ir=0
       DO ik3=1,fpar%kr3
          DO ik2=1,fpar%kr2
             DO ik1=1,fpar%kr1
                ig=ig+1
                IF ((ik3.LE.parm%nr3).AND.(ik2.LE.parm%nr2).AND.(ik1.LE.parm%nr1)) THEN
                   ir=ir+1
                   help=pos(ir)*td_prop%ampl
                   expp=CMPLX(COS(help),-SIN(help),kind=real_8)
                   psi(ig)=psi(ig)*expp
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       ! 
       IF (ifft.EQ.1) THEN
          CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       DO ib=1,nsta
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          fi=crge%f(i,1)
          IF (fi.EQ.0._real_8) fi=2._real_8
          DO ig=1,ncpw%ngw
             fp=psi(nzhs(ig)+ibb)
             fm=psi(indzs(ig)+ibb)
             phi(ig,    i)=fp
             phi(ig+ncpw%ngw,i)=fm
          ENDDO
          IF (geq0)phi(1+ncpw%ngw,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
       ! 
    ENDDO
    ! 
    CALL getnorm_k(phi,nstate,norms)
    ! 
    DO lst=1,nstate
       DO ig=1,ncpw%ngw
          phi(ig,lst)=phi(ig,lst)*(1._real_8/SQRT(norms(lst)))
          phi(ig+ncpw%ngw,lst)=phi(ig+ncpw%ngw,lst)*(1._real_8/SQRT(norms(lst)))
       ENDDO
    ENDDO
    ! 
    CALL getnorm_k(phi,nstate,norms)
    ! 
    DO lst=1,nstate
       IF (paral%io_parent)&
            WRITE(6,'(1x,A,I5,A,F12.6)')&
            'Norm of the perturbed orbital',&
            lst,'   :',norms(lst)
    ENDDO
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kick
  ! ==================================================================
  SUBROUTINE gaugepot_ind(gpot0,gpot1,t_np,time,c0,c2)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: gpot0, gpot1
    INTEGER                                  :: t_np
    REAL(real_8)                             :: time(td_prop%nzeros)
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,crge%n), &
                                                c2(nkpt%ngwk,crge%n)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gaugepot_ind'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: cpr, cpx, cpy, cpz
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:)
    INTEGER                                  :: i, ierr, lrnlsm, nstate
    LOGICAL                                  :: debug, update_d, update_p
    REAL(real_8)                             :: b, d, d1, d2, dt, dydt(2), p, &
                                                y(2), ynew(2)
    REAL(real_8), ALLOCATABLE                :: rhoe(:,:)
    REAL(real_8), SAVE                       :: abst = 0.0_real_8

    nstate=crge%n
    ! call memory(ip_auxc,2*ngwk*maxsys%nax,'auxc')
    ALLOCATE(auxc(MAX((fpar%nnr1+100),nkpt%ngwk*maxsys%nax)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoe(fpar%nnr1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    lrnlsm=lrnlsm+100
    ! 
    debug   =.FALSE.
    update_p=.FALSE.
    update_d=.FALSE.
    ! 
    IF (gaugefield%tkgauge)THEN
       CALL stopgm(' gaugepot_cayl',' tkgauge not implemented ',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    dt=time(1)
    ! 
    abst=abst+dt
    IF ((paral%parent.AND.debug).AND.paral%io_parent)&
         WRITE(6,*) dt,abst
    ! ----------------------------------------------------
    ! CALCULATION OF d (constant term, not dependent on A)
    ! d=0._real_8
    ! do i=1,nstate
    ! call calc_pij(c0(1,i),c0(1,i),cpx,cpy,cpz,1)
    ! if (pertdir.eq.1) then
    ! cpr=cpx
    ! elseif(pertdir.eq.2) then
    ! cpr=cpy
    ! elseif(pertdir.eq.3) then
    ! cpr=cpz
    ! endif
    ! call mp_sum(cpr,2,allgrp)
    ! d=d+f(i)*real(cpr)*tpiba
    ! enddo
    d=0._real_8
    DO i=1,nstate
       CALL calc_pij(c0(1,i),c0(1,i),cpx,cpy,cpz,1)
       IF (td_prop%pertdir.EQ.1) THEN
          cpr=cpx
       ELSEIF (td_prop%pertdir.EQ.2) THEN
          cpr=cpy
       ELSEIF (td_prop%pertdir.EQ.3) THEN
          cpr=cpz
       ENDIF
       CALL mp_sum(cpr,parai%allgrp)
       d=d+crge%f(i,1)*REAL(cpr)*parm%tpiba
    ENDDO
    d1=d
    d=0._real_8
    DO i=1,nstate
       CALL calc_pij(c2(1,i),c2(1,i),cpx,cpy,cpz,1)
       IF (td_prop%pertdir.EQ.1) THEN
          cpr=cpx
       ELSEIF (td_prop%pertdir.EQ.2) THEN
          cpr=cpy
       ELSEIF (td_prop%pertdir.EQ.3) THEN
          cpr=cpz
       ENDIF
       CALL mp_sum(cpr,parai%allgrp)
       d=d+crge%f(i,1)*REAL(cpr)*parm%tpiba
    ENDDO
    d2=d
    d=(d1+d2)/2._real_8
    ! ----------------------------------------------------
    ! CALCULATION OF p (PP dependent part)
    IF (ppcorr) THEN
       IF (update_p) THEN
          ! call rnlsmg(c0,nstate,auxc,lrnlsm,1,0) ! fnl rotated
          ! call rnlsmg(c0,nstate,auxc,lrnlsm,1,1) ! fnlgp
          genl=0._real_8
          CALL rnlrhg(genl,nstate,1)
          CALL mp_sum(genl,parai%allgrp)
       ELSE
          genl=0._real_8
       ENDIF
    ELSE
       genl=0._real_8
    ENDIF
    p=genl
    ! -------------------------------------------
    b=(fpi/parm%omega)*crge%nel
    d=(fpi/parm%omega)*d
    p=(fpi/parm%omega**2)*p
    ! p=(fpi/omega)*p
    ! 29.03.05 ----------------------------------
    d=-1._real_8*d
    b=-1._real_8*b
    p=-1._real_8*p
    IF (paral%io_parent)&
         WRITE(118,'(5f18.12)') abst,gpot,gpot1,ener_com%enl,genl
    ! 29.03.05 ----------------------------------
    y(1)=gpot0
    y(2)=gpot1
    dydt(1)=y(2)
    dydt(2)=d+b*y(1)+p
    CALL rk4_ind(d,b,p,y,dydt,2,abst,dt,ynew,c0,c2,crge%n,update_d)
    gpot0=ynew(1)
    gpot1=ynew(2)
    ! gpot=gpot0
    ! 
    IF (debug) THEN
       IF (paral%io_parent)&
            WRITE(6,'(i4,7e14.5)') parai%me,b,d,genl,p,gpot
    ENDIF
    ! 
    CALL mp_sync(parai%allgrp)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gaugepot_ind
  ! ==================================================================
  SUBROUTINE rk4_ind(d,b,p,y,dydt,ny,t,dt,yout,c0,c2,nstate,&
       update_d)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: d, b, p
    INTEGER                                  :: ny
    REAL(real_8)                             :: dydt(ny), y(ny), t, dt, &
                                                yout(ny)
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    LOGICAL                                  :: update_d

    INTEGER, PARAMETER                       :: nmax = 2 

    INTEGER                                  :: i
    REAL(real_8)                             :: dym(nmax), dyt(nmax), h6, hh, &
                                                th, yt(nmax)

    hh=dt*0.5_real_8
    h6=dt/6._real_8
    th=t+hh
    DO i=1,ny
       yt(i)=y(i)+hh*dydt(i)
    ENDDO
    CALL derivs(1,c0,c2,d,b,p,nstate,update_d,yt,dyt)
    DO i=1,ny
       yt(i)=y(i)+hh*dyt(i)
    ENDDO
    CALL derivs(1,c0,c2,d,b,p,nstate,update_d,yt,dym)
    DO i=1,ny
       yt(i)=y(i)+dt*dym(i)
       dym(i)=dyt(i)+dym(i)
    ENDDO
    CALL derivs(2,c0,c2,d,b,p,nstate,update_d,yt,dyt)
    DO i=1,ny
       yout(i)=y(i)+h6*(dydt(i)+dyt(i)+2._real_8*dym(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rk4_ind
  ! ==================================================================
  SUBROUTINE derivs(step,c0,c2,d,b,p,nstate,update_d,yt,dyt)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: step
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
    REAL(real_8)                             :: d, b, p
    INTEGER                                  :: nstate
    LOGICAL                                  :: update_d
    REAL(real_8)                             :: yt(2), dyt(2)

    COMPLEX(real_8)                          :: cpr, cpx, cpy, cpz
    INTEGER                                  :: i
    REAL(real_8)                             :: d1, d2

    IF (update_d) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'rk4 update d'
       IF (step.EQ.1) THEN
          d=0._real_8
          DO i=1,nstate
             CALL calc_pij(c0(1,i),c0(1,i),cpx,cpy,cpz,1)
             IF (td_prop%pertdir.EQ.1) THEN
                cpr=cpx
             ELSEIF (td_prop%pertdir.EQ.2) THEN
                cpr=cpy
             ELSEIF (td_prop%pertdir.EQ.3) THEN
                cpr=cpz
             ENDIF
             CALL mp_sum(cpr,parai%allgrp)
             d=d+crge%f(i,1)*REAL(cpr)*parm%tpiba
          ENDDO
          d=(fpi/parm%omega)*d
          d1=-1._real_8*d
          d=0._real_8
          DO i=1,nstate
             CALL calc_pij(c2(1,i),c2(1,i),cpx,cpy,cpz,1)
             IF (td_prop%pertdir.EQ.1) THEN
                cpr=cpx
             ELSEIF (td_prop%pertdir.EQ.2) THEN
                cpr=cpy
             ELSEIF (td_prop%pertdir.EQ.3) THEN
                cpr=cpz
             ENDIF
             CALL mp_sum(cpr,parai%allgrp)
             d=d+crge%f(i,1)*REAL(cpr)*parm%tpiba
          ENDDO
          d=(fpi/parm%omega)*d
          d2=-1._real_8*d
          d=(d1+d2)/2._real_8
       ELSEIF (step.EQ.2) THEN
          d=0._real_8
          DO i=1,nstate
             CALL calc_pij(c2(1,i),c2(1,i),cpx,cpy,cpz,1)
             CALL mp_sum(cpz,parai%allgrp)
             d=d+crge%f(i,1)*REAL(cpz)*parm%tpiba
          ENDDO
          d=(fpi/parm%omega)*d
          d=-1._real_8*d
       ENDIF
    ENDIF
    ! update p
    ! ----------------------------------------------------------------==
    dyt(1)=yt(2)
    dyt(2)=d+b*yt(1)+p
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE derivs
  ! ==================================================================
  SUBROUTINE gaugepot_laser(gpot0,gpot1,t)
    ! Updates the vector potential A that drives the electronic dynamics
    ! Coupling given by p \dot A
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: gpot0, gpot1, t

    REAL(real_8)                             :: ampli, osc, swtf

    IF (.NOT. gaugefield%pi_pulse) THEN
       ! pulse shape: E(t)= E0 sin(pi/T t)^2 sin(freq t)
       ! use ampl    for E0
       ! gdfreq  for 1/T
       ! gpot is the gauge pot A=-dE(t)/dt
       ! A=-1/4 * E0 (    Cos[(2 pi/T + freq)t]/(2 pi/T + freq)  
       ! - 2 Cos[freq t]/freq 
       ! - Cos[(2 pi/T - freq)t]/(2 pi/T - freq) 
       ! )
       ampli=td_prop%ampl
       IF (t .GT. gaugefield%gswtt) THEN
          swtf=1._real_8
       ELSE
          swtf=t/gaugefield%gswtt
       ENDIF
       gpot0=-0.25*td_prop%ampl*&
            (&
            COS((2.0*pi*gaugefield%gdfreq+cntr%gfreq)*t)/(2.0*pi*gaugefield%gdfreq+cntr%gfreq) -&
            COS((2.0*pi*gaugefield%gdfreq-cntr%gfreq)*t)/(2.0*pi*gaugefield%gdfreq-cntr%gfreq) -&
            2.0*COS(cntr%gfreq*t)/(cntr%gfreq)&
            )
       ! gpot=gpot0
    ELSE
       osc=gaugefield%muij
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'PI-PULSE, \mu_AB',gaugefield%muij
          CALL m_flush(6)
       ENDIF
       ! F(t)=ampl*sin^2(pi/TP*t)
       ! E(t)=F(t)*cos(w*t)
       ! Tp=n 2 pi/omega
       ! F(t)=ampl*sin^2(w*t/n)
       ! A(t)=\int E(t)=....
       ! tp=nperiod*(2._real_8*Pi/gfreq)
       ! gpot0= ampl*(gfreq/nperiod)*
       ! &       ((-1._real_8+nperiod**2-(nperiod**2)*cos(t*gfreq/nperiod))
       ! &       *sin(t*gfreq)+nperiod*cos(t*gfreq)*sin(t*gfreq/nperiod))/
       ! &        (2*gfreq*(nperiod**2 -1.0))
       gpot0=td_prop%ampl*(&
            -SIN((4*td_prop%ampl*osc-cntr%gfreq)*t)/(4*(4*td_prop%ampl*osc-cntr%gfreq)) +&
            SIN(cntr%gfreq*t)/(2*cntr%gfreq)                           -&
            SIN((4*td_prop%ampl*osc+cntr%gfreq)*t)/(4*(4*td_prop%ampl*osc+cntr%gfreq))&
            )
       ! gpot=gpot0
    ENDIF
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gaugepot_laser
  ! ==================================================================
  SUBROUTINE gaugepot_laser_circ_pol(t)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t

    INTEGER                                  :: id
    LOGICAL                                  :: isin
    REAL(real_8)                             :: ci, nt, twosigsq, w

! if gndir=1: linearily polarized field (along a direction given
! by the only nonzero amplitude). if gndir=2: elliptically pol. 
! isin is a logical variable used to assign a sinusoidal
! behavior to the first nonzero component.

    w=2._real_8*pi*cntr%gfreq
    twosigsq=2._real_8*gaugefield%gswtt*gaugefield%gswtt
    nt=t-gaugefield%gdfreq
    ci=1.0_real_8 !137._real_8
    isin=.FALSE.
    ! write(6,*) 'eli',gdfreq,gfreq,gswtt,gampl(1),gampl(2)
    DO id=1,3
       gpotv(id)=0._real_8
       IF (gampl(id).NE.0._real_8) THEN
          IF (.NOT.isin) THEN
             gpotv(id)=-ci*gampl(id)*EXP(-nt*nt/twosigsq)*SIN(w*nt)/w
             isin=.TRUE.
             ! write(6,*) 'eli',ci,id,gampl(id)
          ELSE
             gpotv(id)=ci*gampl(id)*EXP(-nt*nt/twosigsq)*COS(w*nt)/w
             ! write(6,*) 'eli',ci,id,gampl(id)
          ENDIF
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gaugepot_laser_circ_pol
  ! ==================================================================
  SUBROUTINE dext_pot(t,elfield,icall)
    ! Sets a time-dependent potential to couple with the electronic 
    ! dynamics. 
    ! The shape (frequency) of the field is hardcoded.
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t, elfield
    INTEGER                                  :: icall

    INTEGER                                  :: i, imax, ir, ir1, ir2, ir3
    REAL(real_8)                             :: beta, dfreq, dx, dy, dz, &
                                                freq, fx, fy, fz, nampl, &
                                                swtf, swtt, x, x0, xm, y, y0, &
                                                ym, z, z0, zm

! == Origin ==================================================
! TODO refactor EXTF: defined as 1D, but here is used as 3D... 

    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    ! ------------------------------------------------------------
    dfreq=0.0002
    imax=0

    ! beta=(imax)**2/5.0
    beta=1.0

    swtt=50.0
    ! ------------------------------------------------------------
    ! 
    ! == Mesh parameters =========================================
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)
    ! ============================================================
    ! 
    CALL zeroing(extf)!,kr1*kr2s*kr3s)
    ! 
    ! time dependent switching function
    ! -------------------------------------------------------------
    IF (t .GT. swtt) THEN
       swtf=1._real_8
    ELSE
       swtf=t/swtt
    ENDIF
    ! -------------------------------------------------------------
    ! 
    IF (paral%io_parent)&
         WRITE(6,*) '  swtf     dfreq    tdfreq   ampl'
    IF (paral%io_parent)&
         WRITE(6,'(4F9.4)') swtf,dfreq,td_prop%tdfreq,td_prop%ampl


    xm=moment%dmom(1)
    ym=moment%dmom(2)
    zm=moment%dmom(3)
    ! 
    ! 
    DO ir3=1,fpar%kr3s
       DO ir2=1,fpar%kr2s
          DO ir1=1,fpar%kr1

             x=(parap%nrxpl(parai%mepos,1)+ir1-2)*dx-x0! +dx/2._real_8
             y=(ir2-1)*dy-y0    ! +dy/2._real_8
             z=(ir3-1)*dz-z0    ! +dz/2._real_8
             ir=(ir3-1)*parm%nr1*parm%nr2   +&
                  (ir2-1)*parm%nr1       +&
                  ir1               -&
                  parap%nrxpl(parai%mepos,1) + 1

             ! kill the boundaries
             ! -------------------------------------------------------------
             fx=-1.0+1._real_8/(EXP(-2._real_8*parm%a1(1)*(x-0.1*parm%a1(1)))+1._real_8)+&
                  1._real_8/(EXP(2._real_8*parm%a1(1)*(x-0.9*parm%a1(1)))+1._real_8)
             fy=-1.0+1._real_8/(EXP(-2._real_8*parm%a2(2)*(y-0.1*parm%a2(2)))+1._real_8)+&
                  1._real_8/(EXP(2._real_8*parm%a2(2)*(y-0.9*parm%a2(2)))+1._real_8)
             fz=-1.0+1._real_8/(EXP(-2._real_8*parm%a3(3)*(z-0.1*parm%a3(3)))+1._real_8)+&
                  1._real_8/(EXP(2._real_8*parm%a3(3)*(z-0.9*parm%a3(3)))+1._real_8)

             DO i=-imax,imax
                ! set the frequency and amplitude
                ! -------------------------------------------------------------
                freq=td_prop%tdfreq+i*dfreq
                ! nampl=ampl*exp(-(i**2)/beta)
                ! nampl=nampl/(SQRT(pi*beta))
                ! nampl=nampl*10.0
                nampl=td_prop%ampl
                elfield=nampl*swtf*COS(freq*t)
                ! -------------------------------------------------------------
                IF (td_prop%pertdir.EQ.1) THEN
                   extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        -nampl*xm*swtf*COS(freq*t)*fx
                ELSEIF (td_prop%pertdir.EQ.2) THEN
                   extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        -nampl*ym*swtf*COS(freq*t)*fy
                ELSEIF (td_prop%pertdir.EQ.3) THEN
                   extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        -nampl*zm*swtf*COS(freq*t)*fz
                ELSEIF (td_prop%pertdir.EQ.4) THEN
                   extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        -nampl*(xm*fx+ym*fy)*swtf*COS(freq*t)
                ELSEIF (td_prop%pertdir.EQ.5) THEN
                   extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                        -nampl*(xm*fx+ym*fy+zm*fz)*swtf*COS(freq*t)
                ENDIF
             ENDDO
             ! 
          ENDDO
       ENDDO
    ENDDO
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dext_pot
  ! ==================================================================
  SUBROUTINE dext_pipot(t,elfield)
    ! Sets a time-dependent potential to couple with the electronic 
    ! dynamics. The shape of the field is hardcoded.
    ! Pi-pulse version
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t, elfield

    INTEGER                                  :: ir, ir1, ir2, ir3
    REAL(real_8)                             :: dx, dy, dz, nampl, x, x0, y, &
                                                y0, z, z0

    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    ! ------------------------------------------------------------
    ! == Mesh parameters =========================================
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)
    ! ============================================================
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'computing the pi_pulse'
    ENDIF
    ! 
    CALL zeroing(extf)!,kr1*kr2s*kr3s)
    ! 
    DO ir3=1,fpar%kr3s
       DO ir2=1,fpar%kr2s
          DO ir1=1,fpar%kr1

             x=(parap%nrxpl(parai%mepos,1)+ir1-2)*dx-x0
             y=(ir2-1)*dy-y0+dy/2._real_8
             z=(ir3-1)*dz-z0+dz/2._real_8

             ir=(ir3-1)*parm%nr1*parm%nr2   +&
                  (ir2-1)*parm%nr1       +&
                  ir1               -&
                  parap%nrxpl(parai%mepos,1) + 1

             ! ---------------------------------------------------
             nampl=td_prop%ampl*((SIN(2.0*td_prop%ampl*gaugefield%muij*t))**2)*COS(td_prop%tdfreq*t)
             elfield=nampl
             ! ---------------------------------------------------
             IF (td_prop%pertdir.EQ.1) THEN
                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     =extf((ir3-1)*fpar%kr1*fpar%kr2+(ir2-1)*fpar%kr1+ir1)+nampl*x
             ELSEIF (td_prop%pertdir.EQ.2) THEN
                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     =extf((ir3-1)*fpar%kr1*fpar%kr2+(ir2-1)*fpar%kr1+ir1)+nampl*y
             ELSEIF (td_prop%pertdir.EQ.3) THEN
                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     =extf((ir3-1)*fpar%kr1*fpar%kr2+(ir2-1)*fpar%kr1+ir1)+nampl*z
             ELSEIF (td_prop%pertdir.EQ.4) THEN
                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     =extf((ir3-1)*fpar%kr1*fpar%kr2+(ir2-1)*fpar%kr1+ir1)&
                     +nampl*(x+y)
             ELSEIF (td_prop%pertdir.EQ.5) THEN
                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     =extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)&
                     +nampl*(x+y+z)
             ENDIF
             ! 
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dext_pipot
  ! ==================================================================
  SUBROUTINE dext_ch(t)
    ! Adds external potential generated by a point cahrge in position 
    ! charg1(*) active in the time interval defined by the variables
    ! pcht0 and pchdt
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t

    INTEGER                                  :: ir, ir1, ir2, ir3
    REAL(real_8)                             :: charg1(3), dist1, dx, dy, dz, &
                                                fr, fx, fy, fz, gaussb, r, &
                                                temp, x, x0, y, y0, z, z0

    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    ! -- Mesh parameters --
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)
    charg1(1)=parm%a1(1)*pointcharge%pchfx+dx/2._real_8
    charg1(2)=parm%a2(2)*pointcharge%pchfy+dy/2._real_8
    charg1(3)=parm%a3(3)*pointcharge%pchfz+dz/2._real_8
    ! 
    temp=0.1
    gaussb=0.1
    CALL zeroing(extf)!,kr1*kr2s*kr3s)
    ! ------------------------
    IF ((t.GT.(pointcharge%pcht0-pointcharge%pchdt/2.0)).AND.(t.LT.(pointcharge%pcht0+pointcharge%pchdt/2.0))) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'ADD EXTERNAL PULSE'
       ! ------------------------
       DO ir3=1,fpar%kr3s
          DO ir2=1,fpar%kr2s
             DO ir1=1,fpar%kr1
                x=(parap%nrxpl(parai%mepos,1)+ir1-2)*dx-x0! +dx/2._real_8
                y=(ir2-1)*dy-y0     ! +dy/2._real_8
                z=(ir3-1)*dz-z0     ! +dz/2._real_8
                ir=(ir3-1)*parm%nr1*parm%nr2   +&
                     (ir2-1)*parm%nr1       +&
                     ir1               -&
                     parap%nrxpl(parai%mepos,1) + 1
                ! 
                fx=-1.0+1._real_8/(EXP(-temp*parm%a1(1)*(x-0.15*parm%a1(1)))+1._real_8)+&
                     1._real_8/(EXP(temp*parm%a1(1)*(x-0.85*parm%a1(1)))+1._real_8)
                fy=-1.0+1._real_8/(EXP(-temp*parm%a2(2)*(y-0.15*parm%a2(2)))+1._real_8)+&
                     1._real_8/(EXP(temp*parm%a2(2)*(y-0.85*parm%a2(2)))+1._real_8)
                fz=-1.0+1._real_8/(EXP(-temp*parm%a3(3)*(z-0.15*parm%a3(3)))+1._real_8)+&
                     1._real_8/(EXP(temp*parm%a3(3)*(z-0.85*parm%a3(3)))+1._real_8)
                r=SQRT(x**2+y**2+z**2)
                fr=-1.0+1._real_8/(EXP(-temp*parm%a1(1)*(r-0.1*parm%a1(1)))+1._real_8)+&
                     1._real_8/(EXP(temp*parm%a1(1)*(r-0.9*parm%a1(1)))+1._real_8)

                dist1=SQRT((x-charg1(1))**2+&
                     (y-charg1(2))**2+&
                     (z-charg1(3))**2 )
                ! 
                ! extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)= extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)
                ! &                    +1.0_real_8                       *
                ! &                    (2._real_8*gaussb/Pi)**(3._real_8/4._real_8) *
                ! &                     EXP(-gaussb*dist1)
                ! extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)= extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)
                ! &                    -1.0_real_8                        *
                ! &                    (2._real_8*gaussb/Pi)**(3._real_8/4._real_8) *
                ! &                     EXP(-gaussb*dist2)

                extf((ir3-1)*fpar%kr1*fpar%kr2 + (ir2-1)*fpar%kr1 + ir1)=-pointcharge%pchint *& ! (cos(pfreq*t))
                     1._real_8/dist1
                ! &                    (exp(-dist1/1.0))/dist1
                ! extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)=extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)*fx*fy*fz
                ! extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)=extf((ir3-1)*kr1*kr2 + (ir2-1)*kr1 + ir1)*fr
                ! 
             ENDDO
          ENDDO
       ENDDO
       ! ------------------------
    ENDIF
    ! ------------------------
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dext_ch
  ! ==================================================================
  SUBROUTINE dyn_analysis(first,c0,rhoe,rhoer,rhog,psi,&
       nstate,respo,time)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: first
    COMPLEX(real_8)                          :: c0(*)
    REAL(real_8)                             :: rhoe(:,:), rhoer(*)
    COMPLEX(real_8)                          :: rhog(*), psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: respo(*)
    REAL(real_8)                             :: time

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    INTEGER                                  :: i, i1, i2, ig, ig1, ind, &
                                                indk, ir, ish, maxish, n1, n2
    INTEGER, SAVE                            :: icall = 0
    LOGICAL                                  :: fexist
    REAL(real_8)                             :: ga, gb, gc, gmax, gmin, &
                                                isusc(5000), lindh(5000), &
                                                meangl(5000), rcount, &
                                                rsusc(5000)

    IF (icall.EQ.0) THEN
       IF (paral%parent) THEN
          filen="suscept.dat"
          IF (paral%io_parent)&
               INQUIRE(file=filen,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=177,file=filen,status='unknown')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=177,file=filen,status='unknown',&
                  position='append')
          ENDIF
       ENDIF
    ENDIF
    icall=icall+1
    ! 
    CALL rhoofr_c(c0,rhoe,psi,nstate)
    ! 
    ! if (icall.ne.0) then
    ! -- compute and write density difference against reference density (rhoer)
    ! -- result in rhoe
    DO i=1,fpar%nnr1
       rhoe(i,1)=rhoe(i,1)-rhoer(i)
    ENDDO
    CALL rhopri_eh(c0,tau0,rhoe,psi,nstate,1)
    ! endif
    ! -- transform rhoe in G space
    CALL zeroing(psi(1:fpar%nnr1))!,nnr1)
    CALL zeroing(rhog(1:ncpw%nhg))!,nhg)
    !$omp parallel do private(ir)
    DO ir=1,fpar%nnr1
       psi(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    !$omp parallel do private(ig)
    DO ig=1,ncpw%nhg
       rhog(ig) = psi(nzh(ig))
    ENDDO
    IF (geq0) THEN
       ig1=2
    ELSE
       ig1=1
    ENDIF
    ga=0.1_real_8
    gb=0.1_real_8
    gc=0.1_real_8
    CALL zeroing(respo(1:ncpw%nhg))!,nhg)
    gmax=0.0
    gmin=100000.0
    DO ig=ig1,ncpw%nhg
       IF (hg(ig).GT.gmax) gmax=hg(ig)
       IF (hg(ig).LT.gmin) gmin=hg(ig)
       ! respo(ig)=rhog(ig)/(cos(pfreq*time))
    ENDDO
    gmin=-gmin
    CALL mp_max(gmin,parai%allgrp)
    gmin=-gmin
    gmax=-gmax
    CALL mp_max(gmax,parai%allgrp)
    gmax=-gmax
    IF ((paral%parent.AND.(icall.EQ.1)).AND.paral%io_parent)&
         WRITE(6,*) 'gmin, gmax :', gmin,gmin*parm%tpiba,gmax,gmax*parm%tpiba
    ! write real and imaginary part of the suscettibility
    IF (geq0)THEN
       ig1=2
    ELSE
       ig1=1
    ENDIF
    ! run over the K vector shells
    CALL zeroing(rsusc)!,5000)
    CALL zeroing(isusc)!,5000)
    CALL zeroing(lindh)!,5000)
    rcount=0._real_8
    ind=0
    maxish=100
    DO ish=ig1,ncpw%nhgl
       IF (ish.LE.maxish) THEN
          ind=ind+1
          indk=NINT(gl(ish))
          meangl(indk)=gl(ish)
          rcount=0._real_8
          DO ig=isptr(ish),isptr(ish+1)-1
             rsusc(indk)=rsusc(indk)+REAL(psi(nzh(ig)))
             isusc(indk)=isusc(indk)+AIMAG(psi(nzh(ig)))
             rcount=rcount+1._real_8
          ENDDO
          rsusc(indk)=rsusc(indk)/rcount
          isusc(indk)=isusc(indk)/rcount
       ENDIF
    ENDDO
    CALL mp_sum(rsusc,5000,parai%allgrp)
    CALL mp_sum(isusc,5000,parai%allgrp)
    rsusc=rsusc/parai%nproc
    isusc=isusc/parai%nproc
    CALL m_flush(6)
    IF (paral%parent) THEN
       DO i=1,maxish
          IF (paral%io_parent)&
               WRITE(177,'(5F29.16)')&
               SQRT(REAL(i)),SQRT(REAL(i))*parm%tpiba,&
               rsusc(i),isusc(i)
       ENDDO
       CALL m_flush(177)
    ENDIF
    ! 
    GOTO 111
    CALL mp_max(gmax,parai%allgrp)
    !cmb-bugfix    CALL my_min_d(gmin,1,parai%allgrp) !! this sbr does not exist !!
    IF (paral%io_parent)&
         WRITE(6,*) 'gmax,gmin',gmax,gmin
    cflbod='response.'
    IF (paral%io_parent)&
         WRITE(cipnum,'(i4)') icall
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filen=cflbod(n1:n2)//cipnum(i1:i2)
    CALL densto(respo,tau0,filen)
111 CONTINUE
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dyn_analysis
  ! ==================================================================
  SUBROUTINE currentoper(c0,rhoe,nstate,tau0,psi,nkpoint)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*)
    REAL(real_8)                             :: rhoe(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nkpoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'currentoper'

    CHARACTER(len=14)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: filename
    COMPLEX(real_8)                          :: uimag
    COMPLEX(real_8), ALLOCATABLE             :: c0r(:), dc0(:,:), dc0r(:,:)
    INTEGER                                  :: i, i1, i2, iat, ierr, ig, &
                                                ikind, ir, isp, k, n1, n2, &
                                                nnat
    REAL(real_8)                             :: center(3), maxcurr
    REAL(real_8), ALLOCATABLE                :: curr_r(:,:)

!(nnr1,*)

    ALLOCATE(c0r(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dc0r(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dc0(nkpt%ngwk,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(curr_r(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    center(1) = 0._real_8
    center(2) = 0._real_8
    center(3) = 0._real_8
    nnat = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          center(1) = center(1) + tau0(1,iat,isp)
          center(2) = center(2) + tau0(2,iat,isp)
          center(3) = center(3) + tau0(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    center(1) =  center(1) / REAL(nnat,kind=real_8)
    center(2) =  center(2) / REAL(nnat,kind=real_8)
    center(3) =  center(3) / REAL(nnat,kind=real_8)
    ! 
    ikind=1
    CALL zeroing(curr_r)!,nnr1*3)
    uimag=CMPLX(0._real_8,1._real_8,kind=real_8)
    ! 
    ! see for instance SUBROUTINE CALC_PI
    ! 
    DO i=1,nstate
       DO k=1,3
          DO ig=1,ncpw%ngw
             dc0(ig,k)=uimag*(rk(k,1)+gk(k,ig))*c0(ig,i)*parm%tpiba
             dc0(ig+ncpw%ngw,k)=uimag*(rk(k,1)-gk(k,ig))*c0(ig+ncpw%ngw,i)*parm%tpiba
          ENDDO
          IF (geq0) dc0(1+ncpw%ngw,k)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
       ! 
       CALL zeroing(dc0r)!,nnr1*3)
       DO k=1,3
          CALL zeroing(psi)!,maxfft)
          DO ig=1,ncpw%ngw
             psi(nzhs(ig))=dc0(ig,k)
             psi(indzs(ig))=dc0(ig+ncpw%ngw,k)
          ENDDO
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
          DO ir=1,fpar%nnr1
             dc0r(ir,k)=psi(ir)/SQRT(parm%omega)
          ENDDO
       ENDDO
       CALL zeroing(c0r)!,nnr1)
       CALL zeroing(psi)!,maxfft)
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))=c0(ig,i)
          psi(indzs(ig))=c0(ig+ncpw%ngw,i)
       ENDDO
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       !$omp parallel do private(ir)
       DO ir=1,fpar%nnr1
          c0r(ir)=psi(ir)/SQRT(parm%omega)
       ENDDO
       ! 
       DO k=1,3
          DO ir=1,fpar%nnr1
             curr_r(ir,k)=curr_r(ir,k)+0.5_real_8 * crge%f(i,1) *&
                  AIMAG(&
                  CONJG(c0r(ir))   *       dc0r(ir,k) -&
                  c0r(ir)    *CONJG(dc0r(ir,k))&
                  )
             ! curr_r(ir,k)=curr_r(ir,k)+ f(i,1) *
             ! &             aimag(conjg(c0r(ir))* dc0r(ir,k))
          ENDDO
       ENDDO
    ENDDO
    ! 
    maxcurr=0._real_8
    DO ir=1,fpar%nnr1
       DO k=1,3
          IF (curr_r(ir,k).GT.maxcurr) maxcurr=curr_r(ir,k)
       ENDDO
    ENDDO
    CALL mp_max(maxcurr,parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,*) 'maxcurr :',maxcurr
    IF ((parai%me.EQ.10).AND.paral%io_parent)&
         WRITE(6,*) 'maxcurr :',maxcurr

    ! do ir=1,nnr1
    ! do k=1,3
    ! curr_r(ir,k)=curr_r(ir,k)/maxcurr
    ! enddo
    ! enddo
    ! 
    ! do k=1,3
    ! call currpri(tau0,curr_r(1,1),psi,scr,lscr,nstate,nkpoint,k)
    ! enddo
    ! 
    cflbod="currx.cube_"
    IF (paral%io_parent) WRITE(cipnum,'(I5)') infi
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filename=cflbod(n1:n2)//cipnum(i1:i2)
    CALL cubefile(filename,curr_r(1,1),center,psi,.FALSE.)
    cflbod="curry.cube_"
    IF (paral%io_parent) WRITE(cipnum,'(I5)') infi
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filename=cflbod(n1:n2)//cipnum(i1:i2)
    CALL cubefile(filename,curr_r(1,2),center,psi,.FALSE.)
    cflbod="currz.cube_"
    IF (paral%io_parent) WRITE(cipnum,'(I5)') infi
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filename=cflbod(n1:n2)//cipnum(i1:i2)
    CALL cubefile(filename,curr_r(1,3),center,psi,.FALSE.)
    ! 
    ! compute the velocity field
    DO ir=1,fpar%nnr1
       DO k=1,3
          curr_r(ir,k)=curr_r(ir,k)/rhoe(ir,1)
       ENDDO
    ENDDO
    ! 
    filename="velx.cube"
    CALL cubefile(filename,curr_r(1,1),center,psi,.FALSE.)
    filename="vely.cube"
    CALL cubefile(filename,curr_r(1,2),center,psi,.FALSE.)
    filename="velz.cube"
    CALL cubefile(filename,curr_r(1,3),center,psi,.FALSE.)
    ! 
    ! scale current
    maxcurr=0._real_8
    DO ir=1,fpar%nnr1
       curr_r(ir,1)=SQRT(curr_r(ir,1)**2+&
            curr_r(ir,2)**2+&
            curr_r(ir,3)**2)
       IF (curr_r(ir,1).GT.maxcurr) maxcurr=curr_r(ir,1)
    ENDDO
    ! 
    CALL mp_max(maxcurr,parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,'(A,F10.5)') 'max abs current :', maxcurr
    !$omp parallel do private(ir)
    DO ir=1,fpar%nnr1
       curr_r(ir,1)=curr_r(ir,1)/maxcurr
    ENDDO
    ! 
    ! call currpri(tau0,curr_r(1,1),psi,scr,lscr,nstate,nkpoint,0)
    cflbod="curr_scaled.cube"
    IF (paral%io_parent) WRITE(cipnum,'(I5)') infi
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filename=cflbod(n1:n2)//cipnum(i1:i2)
    CALL cubefile(filename,curr_r(1,1),center,psi,.FALSE.)
    ! 
    DEALLOCATE(c0r,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dc0r,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(curr_r,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE currentoper
  ! ==================================================================
  SUBROUTINE load_ex_states(c0)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,crge%n)

    CHARACTER(*), PARAMETER                  :: procedureN = 'load_ex_states'

    CHARACTER(len=1)                         :: xll
    CHARACTER(len=10)                        :: numer
    CHARACTER(len=100)                       :: file_name
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, ido, ierr, ist, itemp, j, &
                                                nstate, state_number

! -- Arguments
! -- Variables
! (nhg)
! 

    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nstate=crge%n
    DO itemp=1,nstate
       CALL zeroing(vtemp)!,nhg)
       ! state_number=n-(itemp-1)   
       state_number=itemp
       ist=state_number
       IF (ist.LT.0) THEN
          ido=-1
          ist=ABS(ist)
       ELSE
          ido=1
       ENDIF
       IF (ABS(ido).EQ.10) THEN
          xll='a'
       ELSEIF (ABS(ido).EQ.100) THEN
          xll='b'
       ELSE
          xll=' '
       ENDIF
       IF (paral%io_parent)&
            WRITE(numer,'(i6,a1)') ist,xll
       DO i=1,LEN(numer)
          IF (numer(i:i).NE.' ') THEN
             file_name='WAVEFUNCTION.'//numer(i:LEN(numer))
             GOTO 10
          ENDIF
       ENDDO
10     CONTINUE
       CALL densrd(vtemp,file_name)
       DO j=1,ncpw%ngw
          c0(j,state_number)=vtemp(j)
          c0(j+ncpw%ngw,state_number)=CONJG(vtemp(j))
       ENDDO
       IF (geq0)&
            c0(ncpw%ngw+1,state_number)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    ENDDO
    ! 
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE load_ex_states
  ! ==================================================================
  SUBROUTINE  eh_initialize_data(rhoe,rhoer,rhog,v,psi,eirop,elfield,&
       nstate,norms,do_dipole,dmom0,work_time,last_time,&
       gpot_i0,gpot_i1,gpot_l0,gpot_l1,ispin,icall,&
       f_dipole,f_time,f_gauge_i,f_gauge_l,f_el)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:), rhoer(*)
    COMPLEX(real_8)                          :: rhog(*), v(:,:), psi(:), &
                                                eirop(ncpw%nhg)
    REAL(real_8)                             :: elfield
    INTEGER                                  :: nstate
    REAL(real_8)                             :: norms(*)
    LOGICAL                                  :: do_dipole
    REAL(real_8)                             :: dmom0(3), work_time, &
                                                last_time, gpot_i0, gpot_i1, &
                                                gpot_l0, gpot_l1
    INTEGER                                  :: ispin, icall
    CHARACTER(len=30)                        :: f_dipole, f_time, f_gauge_i, &
                                                f_gauge_l, f_el

    CHARACTER(*), PARAMETER :: procedureN = 'eh_initialize_data'

    CHARACTER(len=30)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: eivps(:)
    INTEGER                                  :: ierr, ig, ir, k
    LOGICAL                                  :: debug, fexist
    REAL(real_8)                             :: rhoemax

!(nnr1,*)
! -- Arguments
! -- Variables
! 

    debug=.FALSE.

    ! === Initial dipol
    ! -----------------
    ! 
    IF (icall.EQ.0) THEN
       IF (td_prop%tpointch) THEN
          IF (pointcharge%pchwr.EQ.0) THEN
             CALL dcopy(fpar%nnr1*ispin,rhoe,1,rhoer,1)
             CALL zeroing(psi)!,maxfft*clsd%nlsd)
             DO ir=1,fpar%nnr1
                psi(ir)=CMPLX(rhoer(ir),0._real_8,kind=real_8)
             ENDDO
             CALL  fwfftn(psi,.FALSE.,parai%allgrp)
             DO ig=1,ncpw%nhg
                rhog(ig)=psi(nzh(ig))
             ENDDO
             filen='DENSITY_REF'
             CALL densto(rhog,tau0,filen)
          ELSE
             filen='DENSITY_REF'
             CALL densrd(v,filen)
             CALL ffttor(v,rhoer,psi,ncpw%nhg,.FALSE.)
          ENDIF
       ENDIF
       CALL zeroing(v)!,maxfft*clsd%nlsd)
       CALL zeroing(eirop)!,nhg)
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL eicalc(eivps,eirop)
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       rhoemax=0._real_8
       !$omp parallel do private(ir)
       DO ir=1,fpar%nnr1
          v(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
          IF (rhoe(ir,1).GT.rhoemax) rhoemax=rhoe(ir,1)
       ENDDO
       CALL mp_max(rhoemax,parai%allgrp)
       IF ((paral%parent.AND.debug).AND.paral%io_parent)&
            WRITE(6,*)'MAX DENSITY',rhoemax
       IF (do_dipole) THEN
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL dipo(tau0,eirop,v(:,1))
          CALL dcopy(3,moment%dmom,1,dmom0,1)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1x,a,3f12.6,/)') 'Initial dipole components:',&
                  dmom0(1),dmom0(2),dmom0(3)
             IF (paral%io_parent)&
                  INQUIRE(file=f_dipole,exist=fexist)
             IF (.NOT.fexist) THEN
                IF (paral%io_parent)&
                     OPEN(unit=50,file=f_dipole,status='UNKNOWN')
                work_time=0._real_8
                IF (paral%io_parent)&
                     CLOSE(50)
             ELSE
                IF (paral%io_parent)&
                     OPEN(unit=50,file=f_dipole,status='UNKNOWN')
                IF (paral%io_parent)&
                     REWIND(50)
                IF (paral%io_parent)&
                     READ(50,*) work_time,dmom0(1),dmom0(2),dmom0(3)
                IF (paral%io_parent)&
                     CLOSE(50)
                IF (paral%io_parent)&
                     OPEN(unit=50,file=f_dipole,status='UNKNOWN',&
                     position='APPEND')
                IF (paral%io_parent)&
                     BACKSPACE(50)
                IF (paral%io_parent)&
                     READ(50,*) work_time,moment%dmom(1),moment%dmom(2),moment%dmom(3)
                IF (paral%io_parent)&
                     CLOSE(50)
             ENDIF
          ENDIF
          CALL mp_bcast(work_time,parai%source,parai%allgrp)
          CALL mp_bcast(dmom0,SIZE(dmom0),parai%source,parai%allgrp)
          CALL mp_bcast(moment%dmom,SIZE(moment%dmom),parai%source,parai%allgrp)
          last_time=work_time
       ENDIF
    ENDIF
    ! 
    IF (cntl%tgaugep.OR.cntl%tgaugef.OR.td_prop%td_extpot) THEN
       ! ----------------------------------------
       IF (icall.EQ.0) THEN
          IF (paral%parent) THEN

             IF (cntl%tgaugep) THEN
                IF (paral%io_parent)&
                     INQUIRE(file=f_gauge_i,exist=fexist)
                IF (.NOT.fexist) THEN
                   IF (paral%io_parent)&
                        OPEN(unit=72,file=f_gauge_i,status='UNKNOWN')
                   IF (.NOT.cntl%tgaugef) THEN
                      gpot_i0=td_prop%ampl
                      gpot=td_prop%ampl
                   ENDIF
                   work_time=0._real_8
                   IF (paral%io_parent)&
                        CLOSE(72)
                ELSE
                   IF (paral%io_parent)&
                        OPEN(unit=72,file=f_gauge_i,status='UNKNOWN',&
                        position='APPEND')
                   IF (paral%io_parent)&
                        BACKSPACE(72)
                   IF (paral%io_parent)&
                        READ(72,'(3F20.12)') work_time,gpot_i0,gpot_i1
                   IF (paral%io_parent)&
                        WRITE(6,*) ' Starting induced field ',work_time,gpot_i0,gpot_i1
                   IF (paral%io_parent)&
                        CLOSE(72)
                ENDIF
             ENDIF

             IF (cntl%tgaugef)  THEN
                IF (paral%io_parent)&
                     INQUIRE(file=f_gauge_l,exist=fexist)
                IF (.NOT.fexist.AND.(gndir.EQ.0)) THEN
                   gpot_l0=0._real_8
                   gpot=0._real_8
                   work_time=0._real_8
                ELSEIF (.NOT.fexist.AND.(gndir.GT.0)) THEN
                   DO k=1,3
                      gpotv(k)=0._real_8
                   ENDDO
                   work_time=0._real_8
                ELSE
                   IF (paral%io_parent)&
                        OPEN(unit=73,file=f_gauge_l,status='UNKNOWN',&
                        position='APPEND')
                   IF (paral%io_parent)&
                        BACKSPACE(73)
                   IF (gndir.EQ.0) THEN
                      IF (paral%io_parent)&
                           READ(73,'(3F20.12)') work_time,gpot_l0,gpot_l1
                      IF (paral%io_parent)&
                           PRINT*, ' Starting laser field ',work_time,gpot_l0,gpot_l1
                   ELSEIF (gndir.GT.0) THEN
                      IF (paral%io_parent)&
                           READ(73,'(4F20.12)') work_time,(gpotv(k),k=1,3)
                      IF (paral%io_parent)&
                           PRINT*, ' Starting laser field ',work_time,(gpotv(k),k=1,3)
                   ENDIF
                   IF (paral%io_parent)&
                        CLOSE(73)
                ENDIF
             ENDIF

             IF (td_prop%td_extpot) THEN
                IF (paral%io_parent)&
                     INQUIRE(file=f_el,exist=fexist)
                IF (.NOT.fexist) THEN
                   IF (paral%io_parent)&
                        OPEN(unit=74,file=f_el,status='UNKNOWN')
                   gpot_l0=0._real_8
                   gpot=0._real_8
                   work_time=0._real_8
                   IF (paral%io_parent)&
                        CLOSE(74)
                ELSE
                   IF (paral%io_parent)&
                        OPEN(unit=74,file=f_el,status='UNKNOWN',&
                        position='APPEND')
                   IF (paral%io_parent)&
                        BACKSPACE(74)
                   IF (paral%io_parent)&
                        READ(74,'(2F20.12)') td_prop%ttime,elfield
                   IF (paral%io_parent)&
                        WRITE(6,*) ' Starting field ',td_prop%ttime,elfield
                   work_time=td_prop%ttime
                   IF (paral%io_parent)&
                        CLOSE(74)
                ENDIF
             ENDIF
             ! 
          ENDIF
          IF (cntl%tgaugep.OR.cntl%tgaugef) THEN
             IF (cntl%tgaugef) THEN
                IF (gndir.EQ.0) THEN
                   CALL mp_bcast(gpot_l0,parai%source,parai%allgrp)
                   CALL mp_bcast(gpot_l1,parai%source,parai%allgrp)
                ELSE
                   CALL mp_bcast(gpotv,SIZE(gpotv),parai%source,parai%allgrp)
                ENDIF
             ELSE
                CALL mp_bcast(gpot_i0,parai%source,parai%allgrp)
                CALL mp_bcast(gpot_i1,parai%source,parai%allgrp)
             ENDIF
             CALL mp_bcast(work_time,parai%source,parai%allgrp)
          ELSEIF (td_prop%td_extpot) THEN
             CALL mp_bcast(elfield,parai%source,parai%allgrp)
          ENDIF
          last_time=work_time
          gpot=gpot_l0+gpot_i0
       ENDIF
       ! 
    ELSEIF (td_prop%tpointch.AND.(icall.EQ.0)) THEN
       ! ---------------------------------------
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               INQUIRE(file=f_time,exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=53,file=f_time,status='UNKNOWN')
             work_time=0._real_8
             IF (paral%io_parent)&
                  CLOSE(53)
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=53,file=f_time,status='UNKNOWN',&
                  position='APPEND')
             IF (paral%io_parent)&
                  BACKSPACE(53)
             IF (paral%io_parent)&
                  READ(53,'(3F20.12)') work_time
             IF (paral%io_parent)&
                  CLOSE(53)
          ENDIF
       ENDIF
       CALL mp_bcast(work_time,parai%source,parai%allgrp)
       last_time=work_time
    ENDIF
    ! -----
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eh_initialize_data
  ! ==================================================================
  SUBROUTINE dipolox(posx,rhoe,dipox)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: posx(*), rhoe(:), dipox

    INTEGER                                  :: ig, ir

    ig=1
    dipox=0._real_8
    DO ir=1,fpar%nnr1
       IF (rhoe(ir).NE.0.0_real_8) THEN
          dipox=dipox+rhoe(ir)*posx(ig)
          ig=ig+1
       ENDIF
    ENDDO
    dipox=dipox*(parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s))
    CALL mp_sum(dipox,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dipolox
  ! ==================================================================
  ! SUBROUTINES FOR cheby INTERPOLATION
  ! ==================================================================
  SUBROUTINE chbcf_zeros(n,t,chbcfzx,chbcfzt)
    ! inputs: N, total number of points in time domain [0,...,T]
    ! T, time intervall in a.u.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: t, chbcfzx(n), chbcfzt(n)

    INTEGER                                  :: i
    REAL(real_8)                             :: pi

! -- Arguments:
! -- Variables:
! 

    pi=3.141592653589793_real_8
    DO i=1,n
       chbcfzx(i)=-COS(pi*((i-1)+0.5)/n)
       chbcfzt(i)=0.5*t*(1._real_8+chbcfzx(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chbcf_zeros
  ! ==================================================================
  SUBROUTINE chbcf_delta(chbcfzt,nzeros,tintrvll,deltat)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nzeros
    REAL(real_8)                             :: chbcfzt(nzeros), tintrvll, &
                                                deltat(nzeros)

    INTEGER                                  :: i

! -- Arguments:
! -- Variables:
! 

    deltat(1)=chbcfzt(1)+(tintrvll-chbcfzt(nzeros))
    DO i=2,nzeros
       deltat(i)=chbcfzt(i)-chbcfzt(i-1)
    ENDDO
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chbcf_delta
  ! ==================================================================
  SUBROUTINE chbcf_expan_coeff_i(imat,n,chbcfzx)
    ! inputs:  N      , total number of points in time domain [0,...,T]
    ! chbcfzx, positions of the zeroes in the domain [-1,...,1]
    ! outputs: IMAT      , coefficients matrix
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: imat(n,n), chbcfzx(n)

    INTEGER                                  :: i, j, k
    REAL(real_8)                             :: prec

! -- Arguments:
! -- Variables:
! 

    prec=0.001_real_8
    ! 
    DO i=1,n
       DO j=1,n
          ! imat(i,j)=(1._real_8/(2*n))*chbcff(1,chbcfzx(i))*sfkt(1,chbcfzx(j))
          ! imat(i,j)=(1._real_8/(2*n))*chbcff(1,chbcfzx(i))
          ! &                        *sfkt_num(1,chbcfzx(j),prec)
          imat(i,j)=(1._real_8/(2*n))*chbcff(1,chbcfzx(i))&
               *sfkt_rec(0,chbcfzx(j))
          DO k=2,n
             ! imat(i,j)=imat(i,j)+(1._real_8/n)*
             ! &          chbcff(k,chbcfzx(i))*sfkt(k,chbcfzx(j))
             ! imat(i,j)=imat(i,j)+(1._real_8/n)*
             ! &          chbcff(k,chbcfzx(i))*sfkt_num(k,chbcfzx(j),prec)
             imat(i,j)=imat(i,j)+(1._real_8/n)*&
                  chbcff(k,chbcfzx(i))*sfkt_rec(k-1,chbcfzx(j))
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chbcf_expan_coeff_i
  ! ==================================================================
  FUNCTION chbcff(i_order,x)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_order
    REAL(real_8)                             :: x, chbcff

    INTEGER                                  :: i_value

! -- Arguments:
! -- Variables:
! 

    i_value=i_order-1
    chbcff=COS(i_value*dacos(x))
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION chbcff
  ! ==================================================================
  FUNCTION sfkt(i_order,x)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_order
    REAL(real_8)                             :: x, sfkt

    INTEGER                                  :: i_value
    REAL(real_8)                             :: theta

! -- Arguments:
! -- Variables:
! 

    i_value=i_order-1
    theta=dacos(x)
    IF ((i_value).EQ.1) THEN
       sfkt=(SIN(theta)**2)/2._real_8
    ELSE
       sfkt=(x*COS(i_value*theta) - 1&
            +  (i_value)*SIN(i_value*theta)*SIN(theta))&
            /   ((i_value)**2-1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION sfkt
  ! ==================================================================
  RECURSIVE FUNCTION sfkt_rec(i_value,x) RESULT(res)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_value
    REAL(real_8)                             :: x, res

    REAL(real_8)                             :: theta

! -- Arguments:
! -- Variables:
! 

    theta=dacos(x)
    IF ((i_value).EQ.0) THEN
       res=x+1.0_real_8
    ELSEIF ((i_value).EQ.1) THEN
       res=(x**2-1._real_8)/2._real_8
    ELSE
       res=(2._real_8*(x**2-1) * COS((i_value-1)*theta)  +&
            (i_value-3)*sfkt_rec(i_value-2,x))/(i_value+1)
    ENDIF
    ! ==--------------------------------------------------------------==
  END FUNCTION sfkt_rec
  ! ==================================================================
  FUNCTION sfkt_num(i_order,x,dz)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_order
    REAL(real_8)                             :: x, dz, sfkt_num

    INTEGER                                  :: i, i_value, n_points
    REAL(real_8)                             :: interv, z

! -- Arguments:
! -- Variables:
! 

    i_value=i_order-1
    interv=x+1._real_8
    n_points=NINT(interv/dz)+1

    sfkt_num=0._real_8
    DO i=1,n_points
       z=-1.0_real_8+(i-1)*dz
       sfkt_num=sfkt_num + chbcff(i_order,z) * dz
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION sfkt_num
  ! ==================================================================
  SUBROUTINE initial_guess_eigv(&
       phi2,phil,num,nstate,neigv,chbcfzt,ciclo)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: num, nstate
    COMPLEX(real_8) :: phil(2*ncpw%ngw,nstate,td_prop%nzeros,*), &
      phi2(2*ncpw%ngw,nstate)
    REAL(real_8)                             :: neigv(nstate), &
                                                chbcfzt(td_prop%nzeros)
    INTEGER                                  :: ciclo

    CHARACTER(*), PARAMETER :: procedureN = 'initial_guess_eigv'

    INTEGER                                  :: ierr, ig, lst, t_n
    REAL(real_8)                             :: shift_t
    REAL(real_8), ALLOCATABLE                :: coeff(:)

! -- Arguments:
! -- Variables:
! (nzeros)
! 

    ALLOCATE(coeff(td_prop%nzeros),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    IF (ciclo.EQ.1) THEN
       shift_t=0.0_real_8
    ELSE
       shift_t=chbcfzt(1)
    ENDIF
    ! 
    DO lst=1,nstate
       DO t_n=1,td_prop%nzeros
          coeff(t_n)=(shift_t+chbcfzt(t_n))*neigv(lst)
          DO ig=1,2*ncpw%ngw
             phil(ig,lst,t_n,num)=phi2(ig,lst)&
                  * CMPLX(COS(coeff(t_n)),-SIN(coeff(t_n)),kind=real_8)
          ENDDO
       ENDDO
    ENDDO
    ! 
    DEALLOCATE(coeff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initial_guess_eigv
  ! ==================================================================
  SUBROUTINE initial_guess_rec(&
       phi2,num,nstate,neigv,chbcfzt,t_n,ciclo)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: num, nstate
    COMPLEX(real_8)                          :: phi2(2*ncpw%ngw,nstate)
    REAL(real_8)                             :: neigv(nstate), &
                                                chbcfzt(td_prop%nzeros)
    INTEGER                                  :: t_n, ciclo

    CHARACTER(*), PARAMETER :: procedureN = 'initial_guess_rec'

    INTEGER                                  :: ierr, ig, lst
    REAL(real_8), ALLOCATABLE                :: coeff(:)

! -- Arguments:
! -- Variables:
! (nzeros)
! 

    ALLOCATE(coeff(td_prop%nzeros),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    DO lst=1,nstate
       IF (t_n.EQ.1.AND.ciclo.EQ.1) THEN
          coeff(1)=chbcfzt(1)*neigv(lst)
       ELSEIF (t_n.EQ.1.AND.ciclo.GT.1) THEN
          coeff(1)=chbcfzt(1)*neigv(lst)*2
       ELSE
          coeff(t_n)=(chbcfzt(t_n)-chbcfzt(t_n-1))*neigv(lst)
       ENDIF
       DO ig=1,2*ncpw%ngw
          phi2(ig,lst)=phi2(ig,lst)&
               * CMPLX(COS(coeff(t_n)),-SIN(coeff(t_n)),kind=real_8)
       ENDDO
    ENDDO
    ! 
    DEALLOCATE(coeff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initial_guess_rec
  ! ==================================================================
  SUBROUTINE initial_guess_taylor(&
       phi,phil,hphi,num,psi,nstate,&
       rhoe,v,sc0,&
       chbcfzt,ciclo)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: num
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: hphi(nkpt%ngwk,nstate), &
      phil(nkpt%ngwk,nstate,td_prop%nzeros,*), phi(nkpt%ngwk,nstate)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(maxfftn,1), sc0(nkpt%ngwk,*)
    REAL(real_8)                             :: chbcfzt(td_prop%nzeros)
    INTEGER                                  :: ciclo

    CHARACTER(*), PARAMETER :: procedureN = 'initial_guess_taylor'

    INTEGER                                  :: ierr, ig, ispin, lst, t_n
    REAL(real_8), ALLOCATABLE                :: coeff(:)

!(nnr1,clsd%nlsd)
!rhoe(*)
! -- Arguments:
! -- Variables:
! (nzeros)
! 

    ALLOCATE(coeff(td_prop%nzeros),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ispin=clsd%nlsd
    ! 
    DO t_n=1,td_prop%nzeros

       IF (t_n.EQ.1) THEN
          coeff(1)=chbcfzt(1)
          IF (ciclo.GT.1) coeff(1)=2*coeff(1)
       ELSE
          coeff(t_n)=(chbcfzt(t_n)-chbcfzt(t_n-1))
       ENDIF

       DO lst=1,nstate
          IF (t_n.EQ.1) THEN
             DO ig=1,nkpt%ngwk
                phil(ig,lst,1,num)=phi(ig,lst)
             ENDDO
          ELSE
             DO ig=1,nkpt%ngwk
                phil(ig,lst,t_n,num)=phil(ig,lst,t_n -1,num)
                phi(ig,lst)=phil(ig,lst,t_n,num)
             ENDDO
          ENDIF
       ENDDO

       CALL rhoofr_c(phi,rhoe,psi,nstate)
       CALL vofrho(tau0,fion,rhoe,v,.FALSE.,.FALSE.)
       CALL hpsi(phi,hphi,sc0,rhoe,psi,nstate,1,ispin)

       DO lst=1,nstate
          DO ig=1,nkpt%ngwk
             phil(ig,lst,t_n,num)=phil(ig,lst,t_n,num)-&
                  uimag*hphi(ig,lst)*coeff(t_n)
          ENDDO
       ENDDO

       CALL hpsi(hphi,phi,sc0,rhoe,psi,nstate,1,ispin)

       DO lst=1,nstate
          DO ig=1,nkpt%ngwk
             phil(ig,lst,t_n,num)=phil(ig,lst,t_n,num)-0.5_real_8*&
                  phi(ig,lst)*coeff(t_n)**2
          ENDDO
       ENDDO
       ! 
    ENDDO
    DEALLOCATE(coeff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initial_guess_taylor
  ! ==================================================================
  SUBROUTINE initial_guess_cayley(&
       phi,phil,hphi,num,psi,nstate,&
       rhoe,v,sc0,norms,&
       chbcfzt,ciclo)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: num
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: hphi(nkpt%ngwk,nstate), &
      phil(nkpt%ngwk,nstate,td_prop%nzeros,*), phi(nkpt%ngwk,nstate)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(*), sc0(nkpt%ngwk,*)
    REAL(real_8)                             :: norms(*), &
                                                chbcfzt(td_prop%nzeros)
    INTEGER                                  :: ciclo

    CHARACTER(*), PARAMETER :: procedureN = 'initial_guess_cayley'

    COMPLEX(real_8)                          :: phiold, zdotc
    COMPLEX(real_8), ALLOCATABLE             :: hphi0(:,:), phi0(:,:)
    INTEGER                                  :: ierr, ig, ispin, iter, lst, &
                                                t_n
    REAL(real_8)                             :: delta, deltamax, deltaold, &
                                                tinc2
    REAL(real_8), ALLOCATABLE                :: coeff(:)

! -- Arguments:
! -- Variables:
! (ngwk,*)
! (ngwk,*)
! 

    ALLOCATE(coeff(td_prop%nzeros),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(phi0(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hphi0(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    itermax=1
    ispin=clsd%nlsd

    DO t_n=1,td_prop%nzeros

       IF (t_n.EQ.1) THEN
          coeff(1)=chbcfzt(1)
          IF (ciclo.GT.1) coeff(1)=2*coeff(1)
       ELSE
          coeff(t_n)=(chbcfzt(t_n)-chbcfzt(t_n-1))
       ENDIF
       tinc2=coeff(t_n)/2.0_real_8
       DO lst=1,nstate
          IF (t_n.EQ.1) THEN
             DO ig=1,nkpt%ngwk
                phil(ig,lst,1,num)=phi(ig,lst)
                phi0(ig,lst)=phi(ig,lst)
             ENDDO
          ELSE
             DO ig=1,nkpt%ngwk
                phil(ig,lst,t_n,num)=phil(ig,lst,t_n -1,num)
                phi0(ig,lst)=phil(ig,lst,t_n,num)
             ENDDO
          ENDIF
       ENDDO

       ! CALL rhoofr_c(phi,rhoe,psi,nstate)
       ! CALL VOFRHO(tau0,fion,rhoe,v,.false.,.false.)

       CALL hpsi(phi0,hphi0,sc0,rhoe,psi,nstate,1,ispin)
       CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,hphi0(1,1),1)
       ! 
       DO lst=1,nstate
          CALL dcopy(2*nkpt%ngwk,phi0(1,lst),1,phi(1,lst),1)
       ENDDO

       iter=0
       deltamax=1._real_8
       ! ---------------------------------------------------------------
700    CONTINUE
       deltaold=deltamax
       CALL hpsi(phi,hphi,sc0,rhoe,psi,nstate,1,ispin)
       CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,hphi(1,1),1)

       deltamax=0._real_8
       DO lst=1,nstate
          DO ig=1,nkpt%ngwk
             phiold=phi(ig,lst)

             phi(ig,lst)=phi0(ig,lst)&
                  - uimag * tinc2 * hphi0(ig,lst)&
                  - uimag * tinc2 * hphi(ig,lst)

             delta=ABS(phiold-phi(ig,lst))
             IF (delta.GT.deltamax) deltamax=delta
          ENDDO
       ENDDO
       CALL mp_max(deltamax,parai%allgrp)
       ! write(6,*) t_n,deltamax
       iter=iter+1
       IF ((deltamax.GT.0.001).AND.(deltamax.LE.deltaold)) GOTO 700
       ! ---------------------------------------------------------------

       IF (geq0) CALL zclean_k(phi,nstate,ncpw%ngw)
       DO lst=1,nstate
          norms(lst)=REAL(zdotc(nkpt%ngwk,phi(1,lst),1,phi(1,lst),1))
          CALL mp_sum(norms(lst),parai%allgrp)
       ENDDO
       ! 
       DO lst=1,nstate
          CALL dscal(2*nkpt%ngwk,1._real_8/SQRT(norms(lst)),phi(1,lst),1)
          CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,phil(1,lst,t_n,num),1)
       ENDDO

       IF (iter.GT.itermax) itermax=iter
    ENDDO
    ! 
    IF (paral%io_parent)&
         WRITE(6,*) '  -- cayley max iteration steps',itermax
    DEALLOCATE(coeff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(phi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hphi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initial_guess_cayley

  ! ==================================================================
  SUBROUTINE applymask(c0,psi,nstate,icall,ikind)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    INTEGER                                  :: icall, ikind

    CHARACTER(*), PARAMETER                  :: procedureN = 'applymask'

    COMPLEX(real_8), ALLOCATABLE             :: pdi(:)
    INTEGER                                  :: i, ia, ib, ibb, ierr, ifft, &
                                                ig, ii, ix, iy, iz, lead, &
                                                leadx, lst, njump, nnrx, &
                                                nostat, nsta
    COMPLEX(real_8)                          :: phi(nkpt%ngwk,nstate)
    REAL(real_8), ALLOCATABLE                :: maskfun(:)

    ALLOCATE(maskfun(fpar%nnr1),stat=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pdi(maxfft),stat=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! Checks
    IF (tdgcomm%tdg) THEN
       IF (lspin2%tlse) CALL stopgm(procedureN, &
            'tdg together with lse is not yet implemented ', &
            __LINE__,__FILE__)
       IF(group%nogrp.GT.1)CALL stopgm(procedureN, &
            'tdg together with task groups not yet implemented ', & 
            __LINE__,__FILE__)
       CALL setfftn(ipooldg)
    ELSE
       CALL setfftn(0)
    ENDIF
    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',&
         __LINE__,__FILE__)
    !
    lead  = fpar%kr1s*parai%ngrays
    leadx = fpar%nnr1
    nnrx  = llr1 
    ifft  = 1

    njump=group%nogrp
    nostat = nstate

    DO ia=1,nostat,njump
       ! 
       CALL zeroing(psi)!,maxfft*group%nogrp)
       CALL zeroing(pdi)!,maxfft*group%nogrp)

       nsta=MIN(nstate-ia+1,group%nogrp)
       DO ib=1,nsta
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          DO ig=1,ncpw%ngw
             ! the +g component of the state i
             psi(nzhs(ig)+ibb)=c0(ig,i)
             ! the -g component of the state i
             psi(indzs(ig)+ibb)=c0(ig+ncpw%ngw,i)
          ENDDO
          IF (geq0) psi(nzhs(1)+ibb)=c0(1,i)
       ENDDO
       ! 
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",&
               __LINE__,__FILE__)
       ENDIF
       CALL get_mask_fun(maskfun,1,2)
       !    skipping only for cluster
       ig=0
       ii=0
       DO iz=1,parm%nr3
          DO iy=1,fpar%kr2
             DO ix=1,fpar%kr1
                ig=ig+1
                IF((ix.LE.parm%nr1).AND.(iy.LE.parm%nr2)) THEN
                   ii=ii+1
                   psi(ig)=psi(ig)*maskfun(ii)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !    backward fourier transform psi --> "phi"    
       CALL fwfftn(psi,.TRUE.,parai%allgrp)
       nsta=MIN(nostat-ia+1,group%nogrp)
       DO ib=1,nsta
          !if(tkpnt) then
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          DO ig=1,ncpw%ngw
             phi(ig,i)= psi(nzhs(ig)+ibb)
             phi(ig+ncpw%ngw,i)= psi(indzs(ig)+ibb)
          ENDDO
          IF(geq0) phi(1+ncpw%ngw,i)=CMPLX(0._real_8,0._real_8)
       ENDDO
    ENDDO ! ia  
    ! call getnorm_k(phi,nstate,nos)
    CALL c_clean(phi,nstate,1)
    DO lst=1,nstate
       CALL dcopy(2*nkpt%ngwk,phi(1,lst),1,c0(1,lst),1)
    ENDDO
    !
    DEALLOCATE(maskfun,stat=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pdi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE applymask
  !
  ! ==================================================================
  SUBROUTINE get_mask_fun(pos,switch1,switch2)
    !     COMPUTES the absorbing boundary mask
    !     pos:     output the mask
    !     switch1: switch for momentum distribution (not active)
    !     switch2: switch for centering the mask
    !              1: centering around the COM of the molecule
    !              2: centering at the center of the box
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: pos(fpar%nnr1)
    INTEGER                                  :: switch1, switch2

    INTEGER                                  :: ia, ir, ir1, ir2, ir3, is
    REAL(real_8)                             :: abx, aby, abz, dx, dy, dz, &
                                                help, r2, rmax, rmax2, x, x0, &
                                                xa, xc, xxa, xxc, y, y0, yxa, &
                                                yxc, z, z0, zxa, zxc

! Arguments:
! Variables
!

    xa=maskreal%maskpar1 ! 0.42_real_8
    xc=maskreal%maskpar2 ! 0.48_real_8
    !
    CALL zeroing(pos)!,nnr1)
    ! -- Origin --
    ! call phfac(tau0)
    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    help=0._real_8

    IF (switch2.EQ.1) THEN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             x0=x0+tau0(1,ia,is)*ions0%zv(is)
             y0=y0+tau0(2,ia,is)*ions0%zv(is)
             z0=z0+tau0(3,ia,is)*ions0%zv(is)
             help=help+ions0%zv(is)
          ENDDO
       ENDDO
       x0=x0/help
       y0=y0/help
       z0=z0/help
    ELSE
       x0=parm%a1(1)/2._real_8
       y0=parm%a2(2)/2._real_8
       z0=parm%a3(3)/2._real_8
    ENDIF
    !
    ! -- Mesh parameters --
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)
    ! -- Mask parameters
    xxa=xa*parm%a1(1)!/2._real_8
    xxc=xc*parm%a1(1)!/2._real_8
    yxa=xa*parm%a2(2)!/2._real_8
    yxc=xc*parm%a2(2)!/2._real_8
    zxa=xa*parm%a3(3)!/2._real_8
    zxc=xc*parm%a3(3)!/2._real_8
    !
    ! if (parent) then
    !  write(*,*) 'mask boundaries r1r2(xyz)'
    !  write(*,*) xxa,xxc,yxa,yxc,zxa,zxc
    ! endif
    ! -- Position operator --
    CALL zeroing(pos)!,nnr1)
    rmax=0.5_real_8*parm%a1(1)
    rmax2=rmax**2
    DO ir3=1,parm%nr3
       DO ir2=1,parm%nr2
          DO ir1=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
             x=(ir1-1)*dx-x0+dx/2._real_8
             y=(ir2-1)*dy-y0+dy/2._real_8
             z=(ir3-1)*dz-z0+dz/2._real_8
             ir=(ir3-1)*parm%nr1*parm%nr2   + &
                  (ir2-1)*parm%nr1       + &
                  ir1              - &
                  parap%nrxpl(parai%mepos,1) + 1
             r2=X**2+Y**2+Z**2
             pos(ir)=1._real_8
             ! eq 7 = mom distr
             IF (switch1.EQ.2) THEN
                xxa=xxa*maskreal%maskpar3 ! default 0.8_real_8
                xxc=xxc*maskreal%maskpar3
                yxa=yxa*maskreal%maskpar3
                yxc=yxc*maskreal%maskpar3
                zxa=zxa*maskreal%maskpar3
                zxc=zxc*maskreal%maskpar3
             ENDIF
             ! end eq 7
             abx=ABS(x)
             IF (abx.GE.xxa) THEN
                IF (xa.EQ.xc) THEN
                   pos(ir)=0._real_8
                   GOTO 11
                ENDIF
                IF (abx.LT.xxc) THEN
                   pos(ir)=pos(ir)* &
                        COS((abx-xxa)*pi/(xxc-xxa)/2._real_8)**2
                ELSE
                   pos(ir)=0._real_8
                   GOTO 11
                ENDIF
             ENDIF
             aby=ABS(y)
             IF (aby.GE.yxa) THEN
                IF (xa.EQ.xc) THEN
                   pos(ir)=0._real_8
                   GOTO 11
                ENDIF
                IF (aby.LT.yxc) THEN
                   pos(ir)=pos(ir)* &
                        COS((aby-yxa)*pi/(yxc-yxa)/2._real_8)**2
                ELSE
                   pos(ir)=0._real_8
                   GOTO 11
                ENDIF
             ENDIF
             abz=ABS(z)
             IF (abz.GE.zxa) THEN
                IF (xa.EQ.xc) THEN
                   pos(ir)=0._real_8
                   GOTO 11
                ENDIF
                IF (abz.LT.zxc) THEN
                   pos(ir)=pos(ir)* &
                        COS((abz-zxa)*pi/(zxc-zxa)/2._real_8)**2
                ELSE
                   pos(ir)=0._real_8
                ENDIF
             ENDIF
11           CONTINUE
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE  get_mask_fun
  ! ==================================================================
  ! SUBROUTINES FOR OUTPUT
  ! ==================================================================
  SUBROUTINE wrdens(filename,array,psi)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=256)                       :: filename
    REAL(real_8) :: array(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      psi(fpar%kr2s,fpar%kr3s)

    CHARACTER(len=256)                       :: filename_b
    INTEGER :: fileunit, fileunit_b, i1, i2, i3, ia, iat, ie, ig, ii1, ii2, &
      ii3, ii3p, ilowerleft(3), ir, isp, msgid, msglen, source_pe
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: center(3), int_dens, &
                                                lowerleft(3)

! ==--------------------------------------------------------------==
! Arguments:
! array: the R-space array to plot
! filename: obvious
! tau0: the ionic position vectors
! center: the position (in the CPMD system) where I should place
! \       the center of the plot
! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       ! center(1)=celldm(1)/2._real_8
       ! center(2)=celldm(1)*celldm(2)/2._real_8
       ! center(3)=celldm(1)*celldm(3)/2._real_8
       DO ir=1,3
          center(ir)=(parm%a1(ir)+parm%a2(ir)+parm%a3(ir))/2._real_8
       ENDDO
    ENDIF
    CALL mp_bcast(center,SIZE(center),parai%source,parai%allgrp)

    DO ir=1,3
       lowerleft(ir) = center(ir)&
            - 0.5_real_8*(parm%a1(ir)+parm%a2(ir)+parm%a3(ir))
    ENDDO
    DO ir=1,3
       ilowerleft(ir) = rgrid_m(lowerleft,ir)
    ENDDO
    CALL putongrid_m(lowerleft)
    ! 
    IF (paral%parent) THEN
       fileunit=2604
       CALL xstring(filename,ia,ie)
       ! ori         call fileopen(fileunit,filename(ia:ie),'UNKNOWN')
       IF (paral%io_parent)&
            CALL fileopen(fileunit,filename(ia:ie),FO_NEW,ferror)
       fileunit_b=2605
       filename_b=filename(ia:ie)//'_1d'
       ig=ie+3
       IF (paral%io_parent)&
            CALL fileopen(fileunit_b,filename_b(ia:ig),FO_NEW,ferror)
       IF (paral%io_parent)&
            WRITE (2604,*) 'CPMD CUBE FILE.'
       IF (paral%io_parent)&
            WRITE (2604,*) 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
       IF (paral%io_parent)&
            WRITE (2604,200) ions1%nat, (lowerleft(ir),ir=1,3)
       IF (paral%io_parent)&
            WRITE (2604,200) spar%nr1s,(parm%a1(ir)/REAL(spar%nr1s,kind=real_8),ir=1,3)
       IF (paral%io_parent)&
            WRITE (2604,200) spar%nr2s,(parm%a2(ir)/REAL(spar%nr2s,kind=real_8),ir=1,3)
       IF (paral%io_parent)&
            WRITE (2604,200) spar%nr3s,(parm%a3(ir)/REAL(spar%nr3s,kind=real_8),ir=1,3)

       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             IF (paral%io_parent)&
                  WRITE (2604,201) ions0%iatyp(isp),0.0_real_8,&
                  (tau0(ir,iat,isp),ir=1,3)
          ENDDO
       ENDDO
    ENDIF
    ! 
    DO i1=ilowerleft(1),spar%nr1s+ilowerleft(1)-1
       CALL mp_sync(parai%allgrp)! to avoid long transfer queues
       ii1=i1
       IF (ii1.GT.spar%nr1s) ii1=ii1-spar%nr1s
       ! find the PE who has the plane i1:
       source_pe = 0
       DO WHILE (.NOT.&
            (ii1.GE.parap%nrxpl(source_pe,1) .AND.&
            ii1.LE.parap%nrxpl(source_pe,2)))
          source_pe = source_pe + 1
       ENDDO
       source_pe = parap%pgroup(source_pe+1)

       IF (parai%mepos .EQ. source_pe) THEN
          CALL dcopy(fpar%kr2s*fpar%kr3s,&
               array(ii1-parap%nrxpl(parai%mepos,1)+1,1,1),fpar%kr1,&
               psi,1)
          IF (.NOT. paral%parent) THEN
             msglen = 8*fpar%kr2s*fpar%kr3s! one X-plane.
             msgid = 2! MPI message tag.
             CALL mp_send(psi,fpar%kr2s*fpar%kr3s,parap%pgroup(1),msgid,parai%allgrp)
          ENDIF
       ELSEIF (paral%parent) THEN
          !msglen = 8*kr2s*kr3s! one X-plane.
          msgid = 2     ! MPI message tag.
          CALL mp_recv(psi,fpar%kr2s*fpar%kr3s,source_pe,msgid,parai%allgrp)
       ENDIF
       ! Now, the parent should have the plane in PSI.
       IF (paral%parent) THEN
          int_dens=0._real_8
          DO i2=ilowerleft(2),spar%nr2s+ilowerleft(2)-1
             ii2=i2
             IF (ii2.GT.spar%nr2s) ii2=ii2-spar%nr2s
             DO i3=ilowerleft(3),spar%nr3s+ilowerleft(3)-1-1,2
                ii3=i3
                ii3p=i3+1
                IF (ii3 .GT.spar%nr3s) ii3 =ii3 -spar%nr3s
                IF (ii3p.GT.spar%nr3s) ii3p=ii3p-spar%nr3s
                int_dens=int_dens+psi(ii2,ii3)+psi(ii2,ii3p)
                IF (paral%io_parent)&
                     WRITE (2604,202) psi(ii2,ii3),psi(ii2,ii3p)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE (2605,203) int_dens*(parm%omega)/(spar%nr1s*spar%nr2s*spar%nr3s)
       ENDIF
    ENDDO                     ! i1

    IF ((paral%parent).AND.paral%io_parent)&
         CLOSE(2604)
    IF ((paral%parent).AND.paral%io_parent)&
         CLOSE(2605)

200 FORMAT (5x,i4,6x,3f13.5)
201 FORMAT (5x,i4,f6.2,3f13.8)
202 FORMAT (5x,2f12.3)
203 FORMAT (5x,1f15.7)
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrdens
  ! ==================================================================
  SUBROUTINE putongrid_m(position)
    ! CHANGES the position such that it is placed on the nearest
    ! real space gridpoint. No unit-cell-wrapping is done.
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: position(3)

    INTEGER                                  :: i1, i2, i3

! ==--------------------------------------------------------------==

    i1 = NINT((position(1)*gvec_com%b1(1)&
         +position(2)*gvec_com%b1(2)+position(3)*gvec_com%b1(3))*spar%nr1s/parm%alat)
    i2 = NINT((position(1)*gvec_com%b2(1)&
         +position(2)*gvec_com%b2(2)+position(3)*gvec_com%b2(3))*spar%nr2s/parm%alat)
    i3 = NINT((position(1)*gvec_com%b3(1)&
         +position(2)*gvec_com%b3(2)+position(3)*gvec_com%b3(3))*spar%nr3s/parm%alat)

    position(1) = i1*parm%a1(1)/REAL(spar%nr1s,kind=real_8)&
         + i2*parm%a2(1)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(1)/REAL(spar%nr3s,kind=real_8)
    position(2) = i1*parm%a1(2)/REAL(spar%nr1s,kind=real_8)&
         + i2*parm%a2(2)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(2)/REAL(spar%nr3s,kind=real_8)
    position(3) = i1*parm%a1(3)/REAL(spar%nr1s,kind=real_8)&
         + i2*parm%a2(3)/REAL(spar%nr2s,kind=real_8) + i3*parm%a3(3)/REAL(spar%nr3s,kind=real_8)
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putongrid_m
  ! ==================================================================
  INTEGER FUNCTION rgrid_m(position,direction)
    ! returns the integer grid index ix/iy/iz of an arbitrary position
    ! in space (also outside the box). Convention: 1 <= rgrid <= NR#S
    ! where    # = (direction).
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    REAL(real_8) :: position(3), projection
    INTEGER :: direction, nris, iprojection

    REAL(real_8) :: projector(3)
    ! ==--------------------------------------------------------------==
    IF (direction .EQ. 1) THEN
       projector(:) = gvec_com%b1(:)
       nris=spar%nr1s
    ELSEIF (direction .EQ. 2) THEN
       projector(:) = gvec_com%b2(:)
       nris=spar%nr2s
    ELSEIF (direction .EQ. 3) THEN
       projector(:) = gvec_com%b3(:)
       nris=spar%nr3s
    ELSE
       CALL stopgm('RGRID_M','INVALID DIRECTION',& 
            __LINE__,__FILE__)
    ENDIF

    projection = (position(1)*projector(1) +&
         position(2)*projector(2) +&
         position(3)*projector(3)) /parm%alat

    iprojection = NINT( (projection - NINT(projection - 0.5_real_8))&
         * REAL(nris,kind=real_8)) + 1
    IF (iprojection .GT. nris) iprojection=1

    rgrid_m = iprojection
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION rgrid_m
  ! ==================================================================
  SUBROUTINE tmprd_prop(phi,nstate,itime,tag)
    ! -- Reads complex wavefunction to restart dynamics               --
    ! -- Input: PHI: complex wavefunction in Fourier space            --
    ! --     NSTATE: nb of states                                     --
    ! --     ITIME: index of time in propagation                      --
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: phi(2*ncpw%ngw,nstate)
    INTEGER                                  :: itime
    CHARACTER(len=*)                         :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'tmprd_prop'

    INTEGER                                  :: i, ierr, ig, incr, ip, ipp, &
                                                ipw, len, lst, msgid
    COMPLEX(real_8), ALLOCATABLE             :: phi2(:), phisx(:)
    INTEGER, ALLOCATABLE                     :: mapw(:)

! -- Arguments --
! -- Variables --

    IF (paral%parent) THEN
       len=spar%ngws
    ELSE
       len=ncpw%ngw
    ENDIF
    ALLOCATE(phisx(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       ALLOCATE(phi2(spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! call memory90(phi2,(/(ngw+100)/),'phi2')
       ALLOCATE(mapw(2*(ncpw%ngw+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            OPEN(unit=120,file='wavefunctions',status='old',&
            form='unformatted')
       IF (paral%io_parent)&
            READ(120) itime
    ENDIF
    ! msglen=8
    ! call my_bcast(itime,msglen,source,allgrp)

    DO lst=1,nstate
       DO i=0,1
          incr=i*ncpw%ngw
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  READ(120) (phisx(ig),ig=1,spar%ngws)
          ENDIF

          DO ipp=1,parai%nproc
             ip=parap%pgroup(ipp)
             msgid=ip
             IF (paral%parent) THEN
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      mapw(ig)=mapgp(ig)
                   ENDDO
                ELSE
                   msgid=2
                   !msglen=parap%sparm(3,ipp-1)*8/irat
                   CALL mp_recv(mapw,parap%sparm(3,ipp-1),ip,msgid,parai%allgrp)
                ENDIF
                DO ipw=1,parap%sparm(3,ipp-1)
                   phi2(ipw)=phisx(mapw(ipw))
                ENDDO
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      phi(ig+incr,lst)=phi2(ig)
                   ENDDO
                ELSE
                   msgid=1
                   !msglen=2*parap%sparm(3,ipp-1)*8
                   CALL mp_send(phi2,parap%sparm(3,ipp-1),ip,msgid,parai%allgrp)
                ENDIF
             ELSE
                IF (ip.EQ.parai%me) THEN
                   msgid=2
                   !msglen=ngw*8/irat
                   CALL mp_send(mapgp,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   msgid=1
                   !msglen=2*ngw*8
                   CALL mp_recv(phisx,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   DO ig=1,ncpw%ngw
                      phi(ig+incr,lst)=phisx(ig)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    IF (paral%parent .AND. paral%io_parent) CLOSE(120)
    IF (paral%parent) THEN
       DEALLOCATE(mapw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(phi2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(phisx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tmprd_prop
  ! ==================================================================
  SUBROUTINE tmpwr_prop(phi,nstate,itime,tag)
    ! -- Writes temporary wave function on file TMPPHI                --
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: phi(2*ncpw%ngw,nstate)
    INTEGER                                  :: itime
    CHARACTER(len=*)                         :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'tmpwr_prop'

    COMPLEX(real_8), ALLOCATABLE             :: phi2(:), phisx(:)
    INTEGER                                  :: i, ierr, ig, incr, ip, ipp, &
                                                ipw, len_g, lst, msgid
    INTEGER, ALLOCATABLE                     :: mapw(:)

! -- Arguments --
! -- Variables --

    IF (paral%parent) THEN
       len_g=spar%ngws
    ELSE
       len_g=ncpw%ngw
    ENDIF
    ALLOCATE(phisx(len_g),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       ALLOCATE(phi2(spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(mapw(2*(ncpw%ngw+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            OPEN(unit=120,file=tag,status='unknown',form='unformatted')
       IF (paral%io_parent)&
            WRITE(120) itime
    ENDIF

    DO lst=1,nstate
       DO i=0,1
          incr=i*ncpw%ngw
          DO ipp=1,parai%nproc
             ip=parap%pgroup(ipp)
             msgid=ip
             IF (paral%parent) THEN
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      phi2(ig)=phi(ig+incr,lst)
                      mapw(ig)=mapgp(ig)
                   ENDDO
                ELSE
                   msgid=1
                   !msglen=2*parap%sparm(3,ip)*8
                   CALL mp_recv(phi2,parap%sparm(3,ip),ip,msgid,parai%allgrp)
                   msgid=2
                   !msglen=parap%sparm(3,ip)*8/irat
                   CALL mp_recv(mapw,parap%sparm(3,ip),ip,msgid,parai%allgrp)
                ENDIF
                DO ipw=1,parap%sparm(3,ipp-1)
                   phisx(mapw(ipw))=phi2(ipw)
                ENDDO
             ELSE
                IF (ip.EQ.parai%me) THEN
                   DO ig=1,ncpw%ngw
                      phisx(ig)=phi(ig+incr,lst)
                   ENDDO
                   msgid=1
                   !msglen=2*ngw*8
                   CALL mp_send(phisx,ncpw%ngw,parai%source,msgid,parai%allgrp)
                   msgid=2
                   !msglen=ngw*8/irat
                   CALL mp_send(mapgp,ncpw%ngw,parai%source,msgid,parai%allgrp)
                ENDIF
             ENDIF
          ENDDO
          CALL mp_sync(parai%allgrp)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE (120) (phisx(ig),ig=1,spar%ngws)
          ENDIF
       ENDDO
    ENDDO
    IF ((paral%parent).AND.paral%io_parent)&
         CLOSE(120)
    IF (paral%parent) THEN
       DEALLOCATE(mapw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(phi2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(phisx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tmpwr_prop
  ! ==================================================================
  SUBROUTINE drawsection(rhoe,tag)
    ! -- Draws cross-section of RHOE at Z=0, output in a file         --
    ! -- Input: RHOE density, wavefunction, etc... in real space      --
    ! --        TAG file name                                         --
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:)
    CHARACTER(len=*)                         :: tag

    INTEGER                                  :: ia, ir, ir1, ir2, ir3, is
    REAL(real_8)                             :: charge, dx, dy, dz, x, x0, y, &
                                                y0, z, z0

! nnr1,nrx,nrxs,krx,nrxpl,mepos
! tau0
! na,nsp,zv
! nnr1,nrx,nrxs,krx,nrxpl,mepos
! tau0
! na,nsp,zv
! -- Arguments --
! -- Variables --
! -- Origin --

    CALL phfac(tau0)
    x0=0._real_8
    y0=0._real_8
    z0=0._real_8
    charge=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          x0=x0+tau0(1,ia,is)*ions0%zv(is)
          y0=y0+tau0(2,ia,is)*ions0%zv(is)
          z0=z0+tau0(3,ia,is)*ions0%zv(is)
          charge=charge+ions0%zv(is)
       ENDDO
    ENDDO
    x0=x0/charge
    y0=y0/charge
    z0=z0/charge

    ! -- mesh parameters --
    dx=parm%a1(1)/REAL(spar%nr1s,KIND=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,KIND=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,KIND=real_8)

    ! -- Output --
    IF (paral%io_parent)&
         OPEN(125,file=tag)
    ir3=INT(parm%nr3*z0/parm%a1(1))+1
    DO ir2=1,parm%nr2
       DO ir1=parap%nrxpl(parai%mepos,1),parap%nrxpl(parai%mepos,2)
          x=(ir1-1)*dx-x0
          y=(ir2-1)*dy-y0
          z=(ir3-1)*dz-z0
          ir=(ir3-1)*fpar%kr1*fpar%kr2+(ir2-1)*fpar%kr1+ir1-parap%nrxpl(parai%mepos,1)+1
          IF (paral%io_parent)&
               WRITE(125,'(3f24.18)') x,y,rhoe(ir)
       ENDDO
       IF (paral%io_parent)&
            WRITE (125,*)
    ENDDO
    IF (paral%io_parent)&
         CLOSE(125)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE drawsection
  ! ==================================================================
  SUBROUTINE initialize_ehrenfest_dyn(c0,ch,rhoe,psi,norms)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), INTENT(inout)           :: c0(:,:,:), ch(:,:,:)
    REAL(real_8), INTENT(inout)              :: rhoe(:,:)
    COMPLEX(real_8), INTENT(inout)           :: psi(:,:)
    REAL(real_8), INTENT(inout)              :: norms(:)

    CHARACTER(*), PARAMETER :: procedureN = 'initialize_ehrenfest_dyn'

    COMPLEX(real_8), ALLOCATABLE             :: phi(:)
    INTEGER                                  :: i, ierr, itemp, j

!

    ALLOCATE(phi(nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (cntl%start_real) THEN
       IF (td_prop%read_wf.EQ.1) THEN
          ! LOAD WAVEFUNCTIONS FORM FILE DENSITY.-NM
          CALL load_ex_states(c0)
          CALL getnorm_k(c0,crge%n,norms)
          DO i=1,crge%n
             CALL dscal(2*nkpt%ngwk,1._real_8/SQRT(norms(i)),c0(1,i,1),1)
          ENDDO
       ELSEIF (td_prop%read_wf.EQ.2) THEN
          ! TAKE PROP. WFS FROM THE OPT. WFS IN RESTART
          ! 1,....,tdsp1%nupel, tdsp1%nupel+1,...,crge%n
          ! a,....,a          , b, ...., b
          ! #nupel     <=         #ndownel          
          !IF (paral%io_parent) THEN
          !   write(6,*)
          !   'nupel,ndownel',tdsp1%nupel,tdsp1%ndoel,spin_mod%nsup,spin_mod%nsdown,crge%n
          !ENDIF
          DO i=1,crge%n
             CALL dcopy(2*nkpt%ngwk,ch(1,i,1),1,c0(1,i,1),1)
          ENDDO
       ELSEIF (td_prop%ionize_from_state) THEN
          DO i=1,crge%n
             CALL dcopy(2*nkpt%ngwk,ch(1,i,1),1,c0(1,i,1),1)
          ENDDO
          ! need to reorder orbitals
          DO i=1,nkpt%ngwk
             phi(i)=c0(i,spin_mod%nsup+1,1)
          ENDDO
          DO j=spin_mod%nsup+2,crge%n
             DO i=1,nkpt%ngwk 
                ch(i,j-1,1)=ch(i,j,1)
             ENDDO
          ENDDO
          DO i=1,nkpt%ngwk
             ch(i,crge%n,1)=phi(i)
          ENDDO
          ! do ionization by removing orbital in td_prop%ionized_state
          j=1
          DO i=1,tdsp1%ndoel
             IF (i.EQ.td_prop%ionized_state) THEN
                IF (paral%io_parent) WRITE(6,*) 'skip orbital ',tdsp1%nupel+i
                CYCLE
             ENDIF
             IF (paral%io_parent) WRITE(6,*) 'copy orbital ', tdsp1%nupel+i, 'into ',j
             CALL dcopy(2*nkpt%ngwk,ch(1,tdsp1%nupel+i,1),1,c0(1,j,1),1)
             j=j+1
          ENDDO
       ELSEIF (td_prop%ionize_from_rdm) THEN
          IF (.NOT.cntl%tlsd) THEN
             CALL stopgm(procedureN,'Ionization from gamma matrix only with LSD', &
                  __LINE__,__FILE__)
          ENDIF
          CALL restart_from_ionic_state(ch,c0,psi,norms)
          IF (paral%io_parent) WRITE(6,*) '   Initiated from ionic state'
       ENDIF
    ELSE
       ! READ FROM THE EXTRA INPUT FILE 'wavefunctions'
       CALL tmprd_prop(c0,crge%n,itemp,'wavefunctions')
    ENDIF
    CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
    CALL getnorm_k(c0,crge%n,norms)
    IF (paral%io_parent) THEN
       WRITE(6,*) " NORM OF THE STARTING KS STATES"
       WRITE(6,'(2X,63("-"))')
       WRITE(6,'(7(1X,F8.3))') (norms(i),i=1,crge%n)
       WRITE(6,'(2X,63("-"))')
    ENDIF

    DEALLOCATE(phi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initialize_ehrenfest_dyn
  ! ==================================================================
  SUBROUTINE restart_from_ionic_state(ch,c0,psi,norms)

    COMPLEX(real_8), INTENT(inout)           :: ch(:,:,:), c0(:,:,:), psi(:,:)
    REAL(real_8)                             :: norms(:)

    CHARACTER(*), PARAMETER :: procedureN = 'restart_from_ionic_state'

    COMPLEX(real_8)                          :: new_c1, new_c2
    COMPLEX(real_8), ALLOCATABLE             :: chr(:,:), chr2(:,:), &
                                                chr3(:,:), gamma_mat(:,:), &
                                                psi1(:), xx0(:,:), yy0(:,:)
    INTEGER                                  :: i, ierr, ig, ir, j
    INTEGER, ALLOCATABLE                     :: quadrant(:)
    REAL(real_8)                             :: dist1, dist2, new_p, old_p, &
                                                rr, total_norm, xx, xx1, xx2, &
                                                yy, yy1, yy2

    ALLOCATE(gamma_mat(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi1(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(chr(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xx0(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(yy0(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(chr2(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(chr3(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(quadrant(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! START WITH WAVEFUNCTIONS IN CK AND PUT RESULT IN C0
    ! GO TO REAL TO DO OPERATIONS
    IF (paral%io_parent) THEN
       OPEN(unit=111,file='gamma.dat',status='old')
       DO i=1,crge%n
          READ(111,*) (gamma_mat(i,j),j=1,crge%n)
          !write(*,*)  (gamma_mat(i,j),j=1,crge%n)
       ENDDO
       CLOSE(111)
    ENDIF
    CALL mp_bcast(gamma_mat,SIZE(gamma_mat),parai%source,parai%allgrp)

    ! convert ch orbitals to real space
    DO i=1,crge%n
       CALL zeroing(psi1)!,maxfft)
       DO ig=1,ncpw%ngw
          psi1(nzhs(ig))=ch(ig,i,1)
          psi1(indzs(ig))=ch(ig+ncpw%ngw,i,1)
       ENDDO
       CALL  invfftn(psi1,.TRUE.,parai%allgrp)
       DO ir=1,fpar%nnr1
          chr(ir,i)=psi1(ir)!/SQRT(parm%omega)
       ENDDO
    ENDDO
    ! MAKE A SPECIFIC COMBINATION STARTING FROM THE ORBITALS IN C0 
    DO i=1,crge%n
       DO ir=1,fpar%nnr1
          xx0(ir,i)=REAL(chr(ir,i))
          yy0(ir,i)=AIMAG(chr(ir,i))
       ENDDO
    ENDDO

    CALL zeroing(chr2) 

    DO i=1,crge%n
       DO ir=1,fpar%nnr1
          ! generate the new complex starting orbitals -> chr2(ir,i)
          chr2(ir,i)=chr2(ir,i)+ (1.0_real_8-gamma_mat(i,i)) *chr(ir,i)**2.d0
       ENDDO

       DO j=1,crge%n
          IF (j.EQ.i) CYCLE   !novedad
          DO ir=1,fpar%nnr1
             chr2(ir,i)=chr2(ir,i) - gamma_mat(i,j) * chr(ir,i) * chr(ir,j)
          ENDDO
       ENDDO
    ENDDO

    DO i=1,crge%n
       DO ir=1,fpar%nnr1
          xx=REAL(chr2(ir,i))
          yy=AIMAG(chr2(ir,i))
          rr=SQRT(xx**2+yy**2)
          IF (xx.GT.0.0d0) THEN
             old_p=ATAN(yy/xx)
          ELSEIF(xx.LT.0.0d0 .AND. yy .GE. 0.0d0) THEN
             old_p=ATAN(yy/xx) + pi
          ELSEIF (xx.LT. 0.0d0 .AND. yy .LT. 0.0d0) THEN
             old_p=ATAN(yy/xx) - pi
          ELSEIF(xx.EQ.0.0d0 .AND. yy .GT. 0.0d0) THEN
             old_p=pi/2.0d0
          ELSEIF(xx.EQ.0.0d0 .AND. yy.LT.0.0d0) THEN
             old_p=-pi/2.0d0
          ENDIF
          new_p=old_p/2.0d0
          new_c1=SQRT(rr) * dcmplx(dcos(new_p),dsin(new_p))
          new_p=(old_p+2.0d0*pi)/2.0d0
          new_c2=SQRT(rr) * dcmplx(dcos(new_p),dsin(new_p))
          !
          xx1=REAL(new_c1)
          yy1=AIMAG(new_c1)
          dist1=SQRT( (xx0(ir,i)-xx1)**2 + (yy0(ir,i)-yy1)**2)
          xx2=REAL(new_c2)
          yy2=AIMAG(new_c2)
          dist2=SQRT( (xx0(ir,i)-xx2)**2 + (yy0(ir,i)-yy2)**2)
          IF (dist1.LT.dist2) THEN
             chr3(ir,i)=new_c1
          ELSE
             chr3(ir,i)=new_c2
          ENDIF
       ENDDO
    ENDDO

    ! back to G-space
    DO i=1,crge%n 
       CALL zeroing(psi1)!,maxfft*clsd%nlsd)
       DO ir=1,fpar%nnr1
          psi1(ir)=chr3(ir,i)
       ENDDO
       CALL fwfftn(psi1,.TRUE.,parai%allgrp)
       DO ig=1,ncpw%ngw
          c0(ig,i,1)=psi1(nzhs(ig))
          c0(ig+ncpw%ngw,i,1)=psi1(indzs(ig))
       ENDDO
       IF (geq0) c0(1+ncpw%ngw,i,1)=CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDDO

    CALL getnorm_k(c0,crge%n,norms)
    IF (paral%io_parent) WRITE(6,*) ' RDM norms'
    total_norm=0._real_8
    DO i=1,crge%n
       IF (paral%io_parent) THEN
          WRITE(6,'(I4,1X,F12.8)') i, norms(i)
       ENDIF
       total_norm=total_norm+norms(i)
    ENDDO
    IF (paral%io_parent) WRITE(6,'(A,F12.8)') 'RDM total norm',total_norm

    DEALLOCATE(chr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(chr2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(chr3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamma_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(quadrant,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xx0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(yy0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE restart_from_ionic_state

  ! ==================================================================
END MODULE td_utils
