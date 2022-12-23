MODULE raman_p_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fbohr,&
                                             uimag
  USE d_mat_p_utils,                   ONLY: d_mat_nonloc
  USE ddip,                            ONLY: lenbk,&
                                             ngwmax,&
                                             pdipole,&
                                             pdipolt
  USE ddipo_utils,                     ONLY: setdip
  USE dipomod,                         ONLY: moment
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE opeigr_p_utils,                  ONLY: opeigr_p
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: determ,&
                                             invmat,&
                                             nxxfun
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: raman_p

CONTAINS

  ! ==================================================================
  ! .PAGLIAI/CINECA ADDED TAU0
  SUBROUTINE raman_p(tau0,c0,c1,psi,rhoe,drhoe,&
       eirop,eivps,&
       z11,nstate)
    ! .PAGLIAI/CINECA END ADDED TAU0
    ! ==--------------------------------------------------------------==
    ! .PAGLIAI/CINECA ADDED dipo.inc (DIPOLE MOMENT, see ddipo.F)
    ! .PAGLIAI/CINECA END ADDED dipo.inc (DIPOLE MOMENT, see ddipo.F)

    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'raman_p'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/'x','y','z'/)
    CHARACTER(len=14), PARAMETER             :: fborn = 'APT', &
                                                fdip = 'POLARIZATION', &
                                                fpol = 'POLARIZABILITY'

    COMPLEX(real_8)                          :: dd, det, ephase, one, zero
    COMPLEX(real_8), ALLOCATABLE             :: cw(:), ddmat_p(:,:), &
                                                eirop1(:), f_raman(:,:,:), &
                                                sc0(:), v1(:), v1_loc(:), &
                                                work(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cwork(:,:,:)
    INTEGER :: i, ia, iat, icol, ierr, ig, ig1, il_eirop1, il_framan, &
      il_v1_loc, info, isa, isa0, isp, isub, iv, ix, j, k, k1, k2, nmcol, &
      nop1, nop2, nop3, nxx
    COMPLEX(real_8)                          :: ddmat_temp(nstate,nstate)
    INTEGER, ALLOCATABLE                     :: ipiv(:)
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:,:), mapful(:,:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: aux(3,3), cd_tot, convchem, &
                                                cost, d(3), deb, fac, omegon, &
                                                p(3), phase, tpi, tpzz, &
                                                x0(3), zvtot
    REAL(real_8), ALLOCATABLE                :: born(:,:), born1(:,:), &
                                                charges(:), gmetr(:,:), &
                                                pol(:,:), pol2(:,:)

    CALL tiset('    raman_p',isub)
    one=CMPLX(1._real_8,0._real_8,kind=real_8)
    zero=CMPLX(0._real_8,0._real_8,kind=real_8)
    il_framan=6*ncpw%ngw*nstate ! +12 <-- ??? [IF]
    il_v1_loc=2*ncpw%nhg
    il_eirop1=2*ncpw%nhg
    ALLOCATE(f_raman(ncpw%ngw,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(f_raman)!,SIZE(f_raman))
    ALLOCATE(ddmat_p(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ipiv(20*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(nstate*20),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(v1_loc(il_v1_loc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(v1_loc)!,SIZE(v1_loc))
    ALLOCATE(eirop1(il_eirop1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eirop1)!,SIZE(eirop1))
    ALLOCATE(pol(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(pol)!,9)
    ALLOCATE(pol2(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(pol2)!,9)
    ALLOCATE(gmetr(3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gmetr)!,9)

    lenbk=nxxfun(nstate)
    nxx=MAX(lenbk*parai%nproc,ncpw%ngw*nstate)
    ALLOCATE(sc0(nxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cw(nxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cw)!,nxx)

    IF (ifirst.EQ.0) THEN
       ngwmax=0
       DO i=0,parai%nproc-1
          ngwmax=MAX(ngwmax,parap%sparm(3,i))
       ENDDO
       ALLOCATE(mapful(2,spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nmcol=parai%nproc*ngwmax
       ALLOCATE(mapcol(ngwmax,nmcol/ngwmax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cwork(ncpw%ngw,nstate,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setdip(mapful,mapcol)

       ifirst=1
    ENDIF

    ! born charges
    ALLOCATE(v1(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(born(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(born)!,9*ions1%nat)
    ALLOCATE(born1(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(born1)!,9*ions1%nat)
    ALLOCATE(charges(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(charges)!,ions1%nat)

    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dfnl00(imagp,ions1%nat,maxsys%nhxs,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL rnlsm(c0,nstate,1,1,.TRUE.)
    ! save  fnl and  dfnl in fnl00 and dfnl00      
    isa0=0
    DO isp=1,ions1%nsp
       DO iv=1,nlps_com%ngh(isp)
          DO iat=1,ions0%na(isp)
             isa=isa0+iat
             DO i=1,nstate
                fnl00(1,isa,iv,i,1)=fnl(1,isa,iv,i,1)
             ENDDO
             DO i=1,parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
                DO k=1,3
                   dfnl00(1,isa,iv,k,i,1)=dfnl(1,isa,iv,k,i,1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       isa0=isa0+ions0%na(isp)
    ENDDO

    ! .PAGLIAI/CINECA DIPOLE MOMENT CALCULATION (see ddipo.F)

    IF (cntl%tlsd) THEN
       fac=1._real_8
    ELSE
       fac=2._real_8
    ENDIF
    omegon=1.0_real_8/parm%omega

    ! .....ionic contribution
    ! .....center of charge

    p(1)=0._real_8
    p(2)=0._real_8
    p(3)=0._real_8
    x0(1)=0._real_8
    x0(2)=0._real_8
    x0(3)=0._real_8
    zvtot=0._real_8
    DO isp=1,ions1%nsp
       DO ia=1,ions0%na(isp)
          x0(1)=x0(1)+ions0%zv(isp)*tau0(1,ia,isp)
          x0(2)=x0(2)+ions0%zv(isp)*tau0(2,ia,isp)
          x0(3)=x0(3)+ions0%zv(isp)*tau0(3,ia,isp)
       ENDDO
       zvtot=zvtot+ions0%na(isp)*ions0%zv(isp)
    ENDDO
    x0(1)=x0(1)/zvtot
    x0(2)=x0(2)/zvtot
    x0(3)=x0(3)/zvtot
    DO isp=1,ions1%nsp
       DO ia=1,ions0%na(isp)
          p(1)=p(1)+ions0%zv(isp)*(tau0(1,ia,isp)-x0(1))
          p(2)=p(2)+ions0%zv(isp)*(tau0(2,ia,isp)-x0(2))
          p(3)=p(3)+ions0%zv(isp)*(tau0(3,ia,isp)-x0(3))
       ENDDO
    ENDDO
    p(1)=p(1)*omegon
    p(2)=p(2)*omegon
    p(3)=p(3)*omegon

    IF (ABS(p(1)+p(2)+p(3)) .GT. 1.e-10_real_8)&
         CALL stopgm("DDIPO|","REFERENCE POINT",& 
         __LINE__,__FILE__)

    ! .PAGLIAI/CINECA END DIPOLE MOMENT CALCULATION (see ddipo.F)


    ! ==--------------------------------------------------------------==
    ! ==         the basic loops for perturbations                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/," ",22("*"),a,22("*"),/)')&
         '   perturbations    '
    ! no symmetrization of gradients in this part!

    CALL zeroing(f_raman)!,SIZE(f_raman))

    ! .PAGLIAI/CINECA PERTURBATIONS
    ! THIS SECTION CONTAINS AN UPDATED VERSION OF PREVIOUS 
    ! PERTURBATIONS CALCULATION. NOW IT IS INDEPENDENT
    ! FROM THE CELL SHAPE

    ! .....metric tensor

    gmetr(1,1) = gvec_com%b1(1)
    gmetr(2,1) = gvec_com%b1(2)
    gmetr(3,1) = gvec_com%b1(3)
    gmetr(1,2) = gvec_com%b2(1)
    gmetr(2,2) = gvec_com%b2(2)
    gmetr(3,2) = gvec_com%b2(3)
    gmetr(1,3) = gvec_com%b3(1)
    gmetr(2,3) = gvec_com%b3(2)
    gmetr(3,3) = gvec_com%b3(3)

    ! .....computes position operator (as in ddipo.F, shape independent)

    DO k=1,3
       CALL zeroing(cw)!,ngw*nstate)
       CALL zeroing(ddmat_p)!,nstate*nstate)
       CALL zeroing(cwork)!,2*ngw*nstate)
       CALL zeroing(sc0)!,SIZE(sc0))
       IF (k.EQ.1) THEN
          nop1=1
          nop2=0
          nop3=0
       ELSE IF (k.EQ.2) THEN
          nop1=0
          nop2=2
          nop3=0
       ELSE IF (k.EQ.3) THEN
          nop1=0
          nop2=0
          nop3=3
       ENDIF

       CALL opeigr_p(c0,cw,sc0,nstate,mapful,mapcol,ddmat_p,&
            nop1,nop2,nop3,dd,cwork)

       ! ........dipole moment calculation (as in ddipo.F)
       ! 
       ! ........det(ddmat_p)
       CALL dcopy(2*nstate*nstate,ddmat_p,1,ddmat_temp,1)
       CALL determ(ddmat_temp,nstate,nstate,dd)

       ! .....phase factors due to center of charge motion

       tpzz=parm%tpiba*nstate
       phase=(gmetr(1,k)*x0(1)+gmetr(2,k)*x0(2)+gmetr(3,k)*x0(3))*tpzz
       ephase=CMPLX(COS(phase),SIN(phase),kind=real_8)
       det=dd*ephase
       d(k)=ATAN(AIMAG(det)/REAL(det))

       ! ........polarizability tensor calculation

       CALL zgetrf(nstate,nstate,ddmat_p,nstate,ipiv,info)
       IF (info.NE.0)&
            CALL stopgm('RAMAN_P','error in matrix inversion (1)',& 
            __LINE__,__FILE__)
       CALL zgetri(nstate,ddmat_p,nstate,ipiv,work,nstate,info)
       IF (info.NE.0)&
            CALL stopgm('RAMAN_P','error in matrix inversion (2)',& 
            __LINE__,__FILE__)

       DO i=1,2
          CALL zeroing(cw)!,ngw*nstate)
          CALL zgemm('N','N',ncpw%ngw,nstate,nstate,one,cwork(1,1,i),&
               ncpw%ngw,ddmat_p(1,1),nstate,zero,cw,ncpw%ngw)
          CALL zcopy(ncpw%ngw*nstate,cw,1,cwork(1,1,i),1)
       ENDDO

       ig1=1
       IF (geq0) THEN
          ig1=2
          DO i=1,nstate
             f_raman(1,i,k) = crge%f(i,1) * AIMAG(cwork(1,i,1))
          ENDDO
       ENDIF
       DO i=1,nstate
          DO ig=ig1,ncpw%ngw
             f_raman(ig,i,k) = crge%f(i,1) * 0.5_real_8*uimag*(CONJG(cwork(ig,i,&
                  2))-cwork(ig,i,1))
          ENDDO
       ENDDO
    ENDDO

    ! .PAGLIAI/CINECA END PERTURBATIONS

    DO k=1,3                  ! the direction of the E-field
       ! ==--------------------------------------------------------------==
       CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,v1_loc,f_raman(1,1,k),&
            z11,nstate,eirop1)
       ! ==--------------------------------------------------------------==
       ! calculate the polarizability
       DO k1=1,3
          pol(k,k1)=0._real_8
          ig1=1
          DO i=1,nstate
             pol(k,k1) = pol(k,k1)&
                  +  2._real_8*dotp(ncpw%ngw,c1(:,i),f_raman(:,i,k1))
          ENDDO
       ENDDO


       ! born charges
       CALL ffttog(drhoe,drhoe,psi,ncpw%nhg,.TRUE.)
       CALL rnlsm(c1,nstate,1,1,.TRUE.)
       iat=0
       DO isp=1,ions1%nsp
          DO ia=1,ions0%na(isp)
             iat=iat+1
             DO k2=1,3
                ix=3*(iat-1)+k2
                icol=k
                CALL d_mat_locps(born(ix,icol),drhoe,iat,k2,isp)
                CALL d_mat_nonloc(born(ix,icol),fnl,dfnl,fnl00,dfnl00,&
                     crge%f,wk,isp,k2,k2,iat,nstate,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO                    ! k

    ! -----> DIPOLE MOMENT CALCULATION (see ddipo.F)
    ! .....Sum up dipoles
    tpi=2._real_8*ACOS(-1._real_8)
    cost=fac*omegon/tpi
    DO i=1,3
       pdipole(i)=cost*(d(1)*parm%a1(i)+d(2)*parm%a2(i)+d(3)*parm%a3(i))
       pdipolt(i)=pdipole(i) + p(i)
       moment%dmom(i)=parm%omega*pdipolt(i)! cmb
    ENDDO
    ! -----> END DIPOLE MOMENT CALCULATION (see ddipo.F)

    CALL summat(pol,3)
    CALL mp_sum(born,9*ions1%nat,parai%allgrp)
    CALL invmat(3,gmetr,aux,info)
    CALL dgemm('T','N',3,3,3,1._real_8,gmetr,3,pol,3,0._real_8,aux,3)
    CALL dgemm('N','N',3,3,3,1._real_8/parm%tpiba2,aux,3,gmetr,3,0._real_8,pol2,3)

    ! born charges
    iat=0
    DO isp=1,ions1%nsp
       DO ia=1,ions0%na(isp)
          iat=iat+1
          DO i=1,3
             icol=3*(iat-1)+i
             DO j=1,3
                DO k=1,3
                   born1(icol,j)=born1(icol,j)&
                        - born(icol,k)*gmetr(k,j) / parm%tpiba
                ENDDO
             ENDDO
             charges(iat)=charges(iat)+born1(icol,i)/3._real_8
          ENDDO
       ENDDO
    ENDDO

    ! .PAGLIAI/CINECA BORN CHARGES (standard output printing)
    ! MINOR CHANGES IN THE LAYOUT
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'BORN CHARGES'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) '        ATOM               CHARGE'
       iat=0
       DO isp=1,ions1%nsp
          DO ia=1,ions0%na(isp)
             iat=iat+1
             charges(iat)=charges(iat)+ions0%zv(isp)
             IF (paral%io_parent)&
                  WRITE(6,'(I8,6X,A2,4X,F15.6)')&
                  iat,elem%el(ions0%iatyp(isp)),charges(iat)
             DO i=1,3
                icol=3*(iat-1)+i
                born1(icol,i)=born1(icol,i)+ions0%zv(isp)
             ENDDO
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*)
       ! .PAGLIAI/CINECA END BORN CHARGES (standard output printing)

       ! .PAGLIAI/CINECA DIPOLE MOMENT (standard output printing)
       ! NEW OUTPUT SECTION
       ! -----> DIPOLE MOMENT IN ATOMIC AND CONVENTIONAL UNITS
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' DIPOLE MOMENT (atomic units)'
       deb=2.541765_real_8
       cd_tot=SQRT(moment%dmom(1)**2+moment%dmom(2)**2+moment%dmom(3)**2)
       IF (paral%io_parent)&
            WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
       IF (paral%io_parent)&
            WRITE(6,'(4F12.5)') moment%dmom(1),moment%dmom(2),moment%dmom(3),cd_tot
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' DIPOLE MOMENT (Debye)'
       IF (paral%io_parent)&
            WRITE(6,'(3(10X,A),6X,A)') ' X',' Y',' Z',' TOTAL'
       IF (paral%io_parent)&
            WRITE(6,'(4F12.5)') deb*moment%dmom(1),deb*moment%dmom(2),deb*moment%dmom(3),deb*&
            cd_tot
       IF (paral%io_parent)&
            WRITE(6,*)
       ! .PAGLIAI/CINECA END DIPOLE MOMENT (standard output printing)

       ! .PAGLIAI/CINECA POLARIZABILITY (standard output printing)
       ! NOW PRINTED ALSO IN Angstrom^3
       ! -----> POLARIZABILITY TENSOR IN ATOMIC AND CONVENTIONAL UNITS        
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE (6,*) 'POLARIZABILITY TENSOR: (atomic units)'
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(1,i),i=1,3)
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(2,i),i=1,3)
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(3,i),i=1,3)
       WRITE (6,'(A,E11.5)') ' ChkSum(POLARI) = ',SUM(ABS(pol2))
       convchem=(1.0_real_8/fbohr)**3
       IF (paral%io_parent)&
            WRITE (6,*) 'POLARIZABILITY TENSOR: (Angstrom^3)'
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(1,i)*convchem,i=1,3)
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(2,i)*convchem,i=1,3)
       IF (paral%io_parent)&
            WRITE (6,'(3F15.5)') (pol2(3,i)*convchem,i=1,3)
       ! .PAGLIAI/CINECA END POLARIZABILITY (standard output printing)

       ! .PAGLIAI/CINECA NEW I/O UNITS USED FOR:
       ! .PAGLIAI/CINECA PRINT ON UNIT 91 = 'POLARIZABILITY' 
       ! .PAGLIAI/CINECA PRINT ON UNIT 92 = 'APT' 
       ! .PAGLIAI/CINECA PRINT ON UNIT 93 = 'POLARIZATION' 
       ! .PAGLIAI/CINECA WARNING: IN VERSION 3.11 I/O UNIT 44 USED BY 
       ! .PAGLIAI/CINECA          wannier_center.F AND raman_p.F 

       IF (paral%io_parent)&
            CALL fileopen(91,fpol,fo_app+fo_verb,ferror)
       IF (paral%io_parent)&
            CALL fileopen(92,fborn,fo_app+fo_verb,ferror)
       IF (paral%io_parent)&
            CALL fileopen(93,fdip,fo_app+fo_verb,ferror)
       IF (paral%io_parent)&
            WRITE(91,'(9(2X,1PE16.6E2))')&
            pol2(1,1),pol2(2,2),pol2(3,3),&
            pol2(1,2),pol2(1,3),pol2(2,1),&
            pol2(2,3),pol2(3,1),pol2(3,2)
       IF (paral%io_parent)&
            CALL fileclose(91)
       DO iat=1,ions1%nat
          DO i=1,3
             icol=3*(iat-1)+i
             IF (paral%io_parent)&
                  WRITE(92,*) (born1(icol,j),j=1,3)
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(92)
       IF (paral%io_parent)&
            WRITE(93,'(3F15.5)') moment%dmom(1),moment%dmom(2),moment%dmom(3)
       IF (paral%io_parent)&
            CALL fileclose(93)
    ENDIF                     ! parent

    DEALLOCATE(dfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(charges,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(born1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(born,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gmetr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pol2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ipiv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddmat_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(f_raman,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    raman_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE raman_p


END MODULE raman_p_utils
