MODULE soc
  USE cnst,                            ONLY: fpi,&
                                             pi,&
                                             uimag
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indzs,&
                                             nzhs
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftutil_utils,                   ONLY: phase
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: c_sh,&
                                             lrsym,&
                                             td01,&
                                             tshi,&
                                             tshl
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE lr_tddft_utils,                  ONLY: lr_tddft
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE response_pmod,                   ONLY: nmr_options
  USE soc_types,                       ONLY: &
       c01, cs1, ct1, do_soc_because_sh, ene_sin, ene_tri, etot_isc_md, &
       facthcm, md_on_sin, md_on_tri, norm_s, norm_t, socskip, socvar, tausoc
  USE specpt_utils,                    ONLY: specpt
  USE system,                          ONLY: cnti,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: soc_calc
  PUBLIC :: soc_get_rks
  PUBLIC :: do_tsh_isc
  PUBLIC :: do_tsh_isc_lz
  PUBLIC :: isc_read_start
  PUBLIC :: isc_write_start

CONTAINS

  !==================================================================
  SUBROUTINE soc_calc 
    !==================================================================


    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_calc'

    CHARACTER(len=3)                         :: s_lab(3)
    COMPLEX(real_8)                          :: soc_el(3)
    INTEGER                                  :: i, ierr, is, isub, js, nroot, &
                                                nstate
    LOGICAL                                  :: debug, do_all_socs
    REAL(real_8)                             :: norm_l, soc_mat(100,100), &
                                                soc_tot, soc_vec_gs(100)

    CALL tiset(procedureN,isub)
    s_lab(1)='Sx:'
    s_lab(2)='Sy:'
    s_lab(3)='Sz:'
    debug=.FALSE.
    do_all_socs=.TRUE.

    CALL lr_in
    CALL tddft_input

    nroot=td01%ns_sin
    nstate=crge%n

    IF((cnti%isocs .NE. 0) .OR. (do_all_socs))THEN
       ALLOCATE(cs1(ncpw%ngw,nstate,nroot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ct1(ncpw%ngw,nstate,nroot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c01(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tausoc(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(norm_s(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(norm_t(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       socvar%do_sing=.TRUE.
       CALL specpt
       socvar%do_sing=.FALSE.
       CALL specpt

    ELSEIF((cnti%isocs .EQ. 0) .AND. (.NOT. do_all_socs))THEN
       ALLOCATE(ct1(ncpw%ngw,nstate,nroot+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c01(ncpw%ngw,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tausoc(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(norm_s(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(norm_t(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       td01%ns_sin=0
       socvar%do_sing=.FALSE.
       CALL specpt

    ENDIF

    IF(debug)THEN
       IF(paral%io_parent) WRITE(6,*) 'norms of the lr orbitals'
       IF(paral%io_parent) WRITE(6,*) ' '
       DO i=1,nstate
          norm_l=dotp(ncpw%ngw,c01(:,i),c01(:,i))
          CALL mp_sum(norm_l,parai%allgrp)
          norm_s(i)=norm_l
          IF (paral%io_parent) WRITE(6,'(a,i4,f12.5)') 'norm c0 ',i,norm_l
       ENDDO

       DO i=1,nstate
          norm_l=dotp(ncpw%ngw,cs1(:,i,cnti%isocs),cs1(:,i,cnti%isocs))
          CALL mp_sum(norm_l,parai%allgrp)
          norm_s(i)=norm_l
          IF (paral%io_parent) WRITE(6,'(a,i4,f12.5)') 'norm lr singl',i,norm_l
       ENDDO
       IF(paral%io_parent) WRITE(6,*) ' '
       DO i=1,nstate
          norm_l=dotp(ncpw%ngw,ct1(:,i,cnti%jsoct),ct1(:,i,cnti%jsoct))
          CALL mp_sum(norm_l,parai%allgrp)
          norm_t(i)=norm_l
          IF(paral%io_parent) WRITE(6,'(A,I4,F12.5)') 'norm LR TRIPL',i,norm_l
       ENDDO
       IF(paral%io_parent) WRITE(6,*) ' '
    ENDIF

    CALL zeroing(soc_el)!,3)


    IF(cnti%isocs .EQ. 0)THEN
       CALL soc_get_rks_gs(c01,cnti%jsoct,nstate,soc_el)
    ELSE
       CALL soc_get_rks(c01,cnti%isocs,cnti%jsoct,nstate,soc_el)
    ENDIF

    !IF (paral%io_parent) THEN
    !   DO i=1,3
    !      WRITE(6,'(2F15.7)') real(soc_el(i),kind=real_8),AIMAG(soc_el(i))
    !   ENDDO
    !ENDIF

    IF (paral%io_parent) THEN
       WRITE(6,'(/,A,2I5)') ' SOC components in cm-1 (s/t)',&
            cnti%isocs,cnti%jsoct
       DO i=1,3
          WRITE(6,'(5X,A3,1X,2F15.7)') s_lab(i),REAL(soc_el(i),kind=real_8)*FACTHCM, &
               AIMAG(soc_el(i))*FACTHCM
       ENDDO
    ENDIF

    !vw this is buggus call  zeroing(soc_tot,3)
    soc_tot=0.0_real_8

    DO i=1,3
       soc_tot=soc_tot+REAL(soc_el(i)*CONJG(soc_el(i)),kind=real_8)
    ENDDO
    soc_tot=SQRT(soc_tot)

    IF(paral%io_parent) &
         WRITE(6,'(A,4X,f15.7)') ' TOTAL SOC in cm-1: ',soc_tot*facthcm


    ! -------------------------------------------------------------

    IF (do_all_socs) THEN
       IF (paral%io_parent) THEN 
          WRITE(6,'(/,A31)') &
               ' START SOC LOOP OVER ALL STATES '
          WRITE(6,'(1X,30("-"))')
       ENDIF
       DO is=0,nroot
          DO js=1,nroot
             CALL zeroing(soc_el) !,3)
             IF(is .EQ. 0)THEN
                CALL soc_get_rks_gs(c01,js,nstate,soc_el)
             ELSE
                CALL soc_get_rks(c01,is,js,nstate,soc_el)
             ENDIF

             IF (paral%io_parent) THEN
                WRITE(6,'(/,A,2I5)')' SOC components in cm-1 (s/t): ',&
                     is,js
                DO i=1,3
                   WRITE(6,'(5X,A3,1X,2f15.7)') s_lab(i),REAL(soc_el(i),kind=real_8)*facthcm, &
                        AIMAG(soc_el(i))*facthcm
                ENDDO
             ENDIF
             soc_tot=0.0_real_8
             DO i=1,3
                soc_tot=soc_tot+REAL(soc_el(i)*CONJG(soc_el(i)),kind=real_8)
             ENDDO
             soc_tot=SQRT(soc_tot)

             IF(paral%io_parent)THEN
                WRITE(6,'(A,4X,f15.7)') ' TOTAL SOC in cm-1: ', soc_tot*facthcm
             ENDIF
             CALL m_flush(6)
             IF(is .NE. 0)THEN
                soc_mat(is,js)=soc_tot*facthcm
             ELSE
                soc_vec_gs(js)=soc_tot*facthcm
             ENDIF
          ENDDO
       ENDDO
       IF (paral%io_parent) THEN
          WRITE(6,'(/,A)') ' MATRIX OF SOC COUPLINGS in cm-1'
          WRITE(6,'(A,/)') ' rows: singlets, columns: triplets '
          WRITE(6,'(7X,100(I15,2X))') (js,js=1,nroot)
          DO is=1,nroot
             WRITE(6,'(I5,2X,100(F15.7,2X))') is, &
                  (soc_mat(is,js),js=1,nroot)
          ENDDO
          WRITE(6,'(I5,2X,100(F15.7,2X))') 0, &
               (soc_vec_gs(js),js=1,nroot)
       ENDIF
    ENDIF
    !----------------------------------------------------------------

    IF((cnti%isocs .NE. 0) .OR. (do_all_socs))THEN
       DEALLOCATE(cs1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ct1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c01,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tausoc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(norm_s,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(norm_t,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF((cnti%isocs .EQ. 0) .AND. (.NOT. do_all_socs))THEN
       DEALLOCATE(ct1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(c01,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tausoc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(norm_t,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==-------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE soc_calc
  !==================================================================
  SUBROUTINE soc_get_rks(c0,is,jt,nstate,so)
    !==================================================================
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: is, jt, nstate
    COMPLEX(real_8)                          :: so(3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_get_rks'

    COMPLEX(real_8)                          :: soc_el(3)
    INTEGER                                  :: is1, isub, jt1
    COMPLEX(real_8)                          :: B(3,nstate,nstate), &
                                                A(3,nstate,nstate)
    LOGICAL                                  :: debug
    REAL(real_8)                             :: scal(nstate,nstate), scalb, &
                                                sqtwo

!   COMPLEX(real_8)                          :: scal(nstate,nstate)

    CALL tiset(procedureN,isub)
    sqtwo=1.41421356237
    debug=.FALSE.

    DO is1=1,nstate
       CALL soc_matelement(cs1(:,:,is),ct1(:,:,jt),nstate,A(:,is1,is1),is1,is1)
       DO jt1=1,nstate
          CALL soc_matelement(c0,c0,nstate,B(:,is1,jt1),is1,jt1)
          scalb=dotp(ncpw%ngw,cs1(:,is1,is),ct1(:,jt1,jt))
          CALL mp_sum(scalb,parai%allgrp)
          CALL mp_bcast(scalb,parai%source,parai%allgrp)
          scal(is1,jt1)=scalb
       ENDDO
    ENDDO

    CALL zeroing(soc_el)!,3)

    DO is1=1,nstate
       so(1)=so(1)+A(1,is1,is1)+uimag*A(2,is1,is1)
       so(2)=so(2)+A(1,is1,is1)-uimag*A(2,is1,is1)
       so(3)=so(3)+2._real_8*(A(3,is1,is1)-B(3,is1,is1)*scal(is1,is1))
       DO jt1=1,nstate
          so(1)=so(1)-scal(is1,jt1)*(B(1,jt1,is1)+uimag*B(2,is1,is1))
          so(2)=so(2)-scal(is1,jt1)*(B(1,jt1,is1)-uimag*B(2,is1,is1))
          IF(is1 .NE. jt1)THEN
             so(3)=so(3)-2._real_8*B(3,jt1,is1)*scal(is1,jt1)
          ENDIF
       ENDDO
    ENDDO

    so(1)=so(1)*sqtwo  
    so(2)=so(2)*sqtwo


    CALL tihalt(procedureN,isub)
  END SUBROUTINE soc_get_rks

  SUBROUTINE soc_get_rks_gs(c0,jt,nstate,so)
    !====================================================================
    !               GROUND-STATE-EXCITED-STATE COUPLING
    !====================================================================
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: jt, nstate
    COMPLEX(real_8)                          :: so(3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_get_rks_gs'

    COMPLEX(real_8)                          :: matel(3)
    INTEGER                                  :: is1, is1r, isub, jt1r
    LOGICAL                                  :: debug
    REAL(real_8)                             :: sqtwo

    CALL tiset(procedureN,isub)
    sqtwo=1.41421356237
    debug=.FALSE.

    DO is1=1,nstate
       is1r=is1
       jt1r=is1r
       CALL soc_matelement(c0,ct1(:,:,jt),nstate,matel,is1r,jt1r)
       so(3)=so(3)+2._real_8*matel(3)

       so(1)=so(1)+sqtwo*matel(1) &
            -sqtwo*matel(2)*uimag

    ENDDO

    so(2)=CONJG(so(1))

    CALL tihalt(procedureN,isub)
  END SUBROUTINE soc_get_rks_gs
  !==================================================================
  SUBROUTINE soc_matelement(cl,cr,nstate,matel,is1r,jt1r)
    !==================================================================
    !
    !     | x1/r3 |       | p1 |     | s23= x2/r3 p3  - s32=x3/r3 p2 |
    !     | x2/r3 |   x   | p2 |  =  | s31= x3/r3 p1  - s13=x1/r3 p3 |
    !     | x3/r3 |       | p3 |     | s12= x1/r3 p2  - s21=x2/r3 p1 |
    !
    COMPLEX(real_8)                          :: cl(:,:), cr(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: matel(:)
    INTEGER                                  :: is1r, jt1r

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_matelement'

    CHARACTER(LEN=20)                        :: filename
    COMPLEX(real_8), ALLOCATABLE :: cr12(:), cr13(:), cr21(:), cr23(:), &
      cr31(:), cr32(:), crx(:), crxrg(:), cry(:), cryrg(:), crz(:), crzrg(:), &
      psi(:)
    INTEGER                                  :: i, ia, iat, ierr, ig, ig1, &
                                                ir, is, isa, isp, isub, nnat
    LOGICAL                                  :: debug, dophase, &
                                                multiplybyi_f, multiplybyi_l, &
                                                print_cube, skip_contrib
    REAL(real_8)                             :: alpha, center(3), eigr_dot2, &
                                                fact, norm_l, norm_r, omtp, &
                                                s12r, s13r, s21r, s23r, s31r, &
                                                s32r
    REAL(real_8), ALLOCATABLE                :: clr(:), crxr(:), cryr(:), &
                                                crzr(:)

    CALL tiset(procedureN,isub)
    skip_contrib=.FALSE.
    nmr_options%tfast=.TRUE.
    debug=.FALSE.
    print_cube=.FALSE.
    center(1) = 0_real_8
    center(2) = 0_real_8
    center(3) = 0_real_8
    nnat = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          center(1) = center(1) + tausoc(1,iat,isp)
          center(2) = center(2) + tausoc(2,iat,isp)
          center(3) = center(3) + tausoc(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    center(1) =  center(1) / REAL(nnat,real_8)
    center(2) =  center(2) / REAL(nnat,real_8)
    center(3) =  center(3) / REAL(nnat,real_8)
    !
    multiplybyi_f=.TRUE.
    multiplybyi_l=.FALSE.
    dophase=.TRUE.
    !
    ! skip calculation if the norm of the LR orbitals is less than 10D-6
    norm_r=dotp(ncpw%ngw,cl(:,is1r),cl(:,is1r))
    norm_l=dotp(ncpw%ngw,cr(:,jt1r),cr(:,jt1r))
    CALL mp_sum(norm_r,parai%allgrp)
    CALL mp_sum(norm_l,parai%allgrp)
    !IF (paral%parent.AND.((norm_l*norm_r).LT.1.0D-8)) skip_contrib=.TRUE.
    !CALL mp_bcast(skip_contrib,parai%source,parai%allgrp)
    skip_contrib = (norm_l*norm_r).LT.1.0e-8_real_8
    IF (skip_contrib) THEN
       DO i=1,3
          matel(i)=CMPLX(0.0,0.0)
       ENDDO
       CALL tihalt(procedureN,isub)
       RETURN
    ENDIF
    !
    !   if (parent) write(*,*) 'multiplybyi,dophase',multiplybyi,dophase
    !   if (parent) write(*,*) 'omega'              ,omega,SQRT(omega)
    !   
    !   PRINT *,'basile', NSTATE,MATEL(3),IS1R,JT1R,SIG,TAU
    !   
    !   --------------------------------------------------
    !   norm=dotp(ngw,cl(1,nstate),cl(1,nstate)) 
    !   call mp_sum(norm,parai%allgrp)
    !   if (paral%io_parent) write(*,*) 'norm in soc 3: NSTATE',norm
    !   norm=dotp(ngw,cl(1,1),cl(1,1)) 
    !   call mp_sum(norm,parai%allgrp)
    !   if (paral%io_parent) write(*,*) 'norm in soc 3: 1',norm
    !   ---------------------------------------------------
    !   
    ALLOCATE(crx(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cry(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)  
    ALLOCATE(crz(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !vw XL with 4 mpi -> 1 mpi hangs in the following allocate
    ALLOCATE(crxr(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cryr(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)  
    ALLOCATE(crzr(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(clr(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(crxrg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cryrg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(crzrg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr12(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr21(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr13(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr31(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr23(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cr32(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !   
    !   fine structure cst.!
    alpha = (1._real_8/137.0359895_real_8)
    CALL zeroing(matel)!,3)
    CALL zeroing(cr12)!,nhg)
    CALL zeroing(cr21)!,nhg)
    CALL zeroing(cr31)!,nhg)
    CALL zeroing(cr13)!,nhg)
    CALL zeroing(cr23)!,nhg)
    CALL zeroing(cr32)!,nhg)
    CALL zeroing(crx)!,ngw)
    CALL zeroing(cry)!,ngw)
    CALL zeroing(crz)!,ngw)
    !   
    DO ig=1,ncpw%ngw
       crx(ig)=gk(1,ig)*cr(ig,jt1r)*parm%tpiba
       cry(ig)=gk(2,ig)*cr(ig,jt1r)*parm%tpiba
       crz(ig)=gk(3,ig)*cr(ig,jt1r)*parm%tpiba
    ENDDO
    !   this call returns  crx=p_1.cr=-i d/dx cr
    !   
    !   make crx real (multiplication by -i)
    IF (multiplybyi_f) THEN
       CALL zscal(ncpw%ngw,-uimag,crx(1),1)
       CALL zscal(ncpw%ngw,-uimag,cry(1),1)
       CALL zscal(ncpw%ngw,-uimag,crz(1),1)
    ENDIF
    !   
    CALL zeroing(crxr)!,nnr1)
    CALL zeroing(cryr)!,nnr1)
    CALL zeroing(crzr)!,nnr1)
    CALL zeroing(clr) !,nnr1)
    !   
    CALL zeroing(psi)!,maxfft)
    !   
    CALL fftTOr(crx,       crxr,psi,ncpw%ngw,dophase)
    CALL fftTOr(cry,       cryr,psi,ncpw%ngw,dophase)
    CALL fftTOr(crz,       crzr,psi,ncpw%ngw,dophase)
    CALL fftTOr(cl(1,is1r),clr, psi,ncpw%ngw,dophase)

    !   --------------------------------------------------
    !   norm=0._real_8
    !   do ir=1,nnr1
    !    norm=norm+crxr(ir)*crxr(ir)/parm%omega
    !   enddo
    !   norm=norm*parm%omega/dble(spar%nr1s*spar%nr2s*spar%nr3s)
    !   call mp_sum(norm,parai%allgrp)
    !   if (paral%io_parent) write(*,*) 'nxorm',norm
    !   --------------------------------------------------
    !   
    DO ir=1,fpar%nnr1
       crxr(ir)=crxr(ir)*clr(ir) 
       cryr(ir)=cryr(ir)*clr(ir)
       crzr(ir)=crzr(ir)*clr(ir)
    ENDDO
    !   
    !   --------------------------------------------------
    !   note: real space orbitals need to be normalized
    !         (/SQRT(omega)) when used as such.
    !     For internal manipulation (FFT/INVFFT) 
    !     within CPMD they don t need normalization.
    !   --------------------------------------------------
    !   norm=0._real_8
    !   if (parent) write(*,*) 'is1r,jt1r',is1r,jt1r
    !   do ir=1,nnr1
    !    norm=norm+crxr(ir)*crxr(ir)/parm%omega
    !   enddo
    !   norm=norm*parm%omega/dble(spar%nr1s*spar%nr2s*spar%nr3s)
    !   call mp_sum(norm,parai%allgrp)
    !   if (paral%io_parent) write(*,*) 'nxorm',norm
    !   ---------------------------------------------------
    !   
    CALL zeroing(crxrg)!,nhg)
    CALL zeroing(cryrg)!,nhg)
    CALL zeroing(crzrg)!,nhg)
    !   
    CALL fftTOg(crxr,crxrg,psi,ncpw%nhg,dophase)
    CALL fftTOg(cryr,cryrg,psi,ncpw%nhg,dophase)
    CALL fftTOg(crzr,crzrg,psi,ncpw%nhg,dophase)
    !   
    ig1=1
    IF (geq0) THEN
       ig1=2
       cr12(1)=CMPLX(0.0_real_8,0.0_real_8)
       cr21(1)=CMPLX(0.0_real_8,0.0_real_8)
       cr23(1)=CMPLX(0.0_real_8,0.0_real_8)
       cr32(1)=CMPLX(0.0_real_8,0.0_real_8)
       cr13(1)=CMPLX(0.0_real_8,0.0_real_8)
       cr31(1)=CMPLX(0.0_real_8,0.0_real_8)
    ENDIF
    !   Computer Physics Communications, 147, (2002) 707–710
    DO ig=ig1,ncpw%nhg
       cr12(ig)=-fpi*uimag*(gk(1,ig)/hg(ig))*cryrg(ig)/parm%tpiba ! s12
       cr21(ig)=-fpi*uimag*(gk(2,ig)/hg(ig))*crxrg(ig)/parm%tpiba ! s21
       cr23(ig)=-fpi*uimag*(gk(2,ig)/hg(ig))*crzrg(ig)/parm%tpiba ! s23
       cr32(ig)=-fpi*uimag*(gk(3,ig)/hg(ig))*cryrg(ig)/parm%tpiba ! s32
       cr13(ig)=-fpi*uimag*(gk(1,ig)/hg(ig))*crzrg(ig)/parm%tpiba ! s13
       cr31(ig)=-fpi*uimag*(gk(3,ig)/hg(ig))*crxrg(ig)/parm%tpiba ! s31
    ENDDO
    !   
    IF (multiplybyi_l) THEN
       CALL zscal(ncpw%nhg,uimag,cr12(1),1)
       CALL zscal(ncpw%nhg,uimag,cr21(1),1)     
       CALL zscal(ncpw%nhg,uimag,cr23(1),1)
       CALL zscal(ncpw%nhg,uimag,cr32(1),1)
       CALL zscal(ncpw%nhg,uimag,cr13(1),1)
       CALL zscal(ncpw%nhg,uimag,cr31(1),1)
    ENDIF

    IF (print_cube) THEN
       CALL fftTOr(cl(1,nstate),crxr,psi,ncpw%ngw,.TRUE.)
       filename="homo.cube"
       CALL cubefile(filename,crxr,center,psi,.FALSE.)
       CALL fftTOr(cl(1,nstate-1),crxr,psi,ncpw%ngw,.TRUE.)
       filename="homo-1.cube"
       CALL cubefile(filename,crxr,center,psi,.FALSE.)

       !    call fftTOr(cr12(1),crxr,psi,nhg,dophase)
       !    if ((is1r.eq.nstate).and.(jt1r.eq.(nstate-1))) then
       !     filename="last12.cube"
       !     call cubefile(filename,crxr,center,psi,.false.)
       !    endif
       !        --------------------------------------------------
       !    call fftTOr(cr21(1),crxr,psi,nhg,dophase)
       !    if ((is1r.eq.nstate).and.(jt1r.eq.(nstate-1))) then
       !     filename="last21.cube"
       !     call cubefile(filename,crxr,center,psi,.false.)
       !    endif
    ENDIF
    !   call phfac(tausoc)
    !   
    omtp=1.0_real_8/(parm%omega)


    isa=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          isa=isa+1
          s12r=eigr_dot2(isa,cr12)
          s21r=eigr_dot2(isa,cr21)
          s13r=eigr_dot2(isa,cr13)
          s31r=eigr_dot2(isa,cr31)
          s23r=eigr_dot2(isa,cr23)
          s32r=eigr_dot2(isa,cr32)
          CALL mp_sum(s12r,parai%allgrp)
          CALL mp_sum(s21r,parai%allgrp)
          CALL mp_sum(s13r,parai%allgrp)
          CALL mp_sum(s31r,parai%allgrp)
          CALL mp_sum(s23r,parai%allgrp)
          CALL mp_sum(s32r,parai%allgrp)
          s12r=s12r*omtp
          s21r=s21r*omtp
          s13r=s13r*omtp
          s31r=s31r*omtp
          s23r=s23r*omtp
          s32r=s32r*omtp

          !          fact=(alpha**2)/2._real_8 * ions0%zv(is)
          fact=((alpha**2)/2._real_8)*DBLE(ions0%iatyp(is))
          !
          matel(1)=matel(1)+(s23r-s32r)*fact*0.5_real_8
          matel(2)=matel(2)+(s31r-s13r)*fact*0.5_real_8
          matel(3)=matel(3)+(s12r-s21r)*fact*0.5_real_8


       ENDDO
    ENDDO
    !   --------------------------------------------------------------
    DEALLOCATE(crx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cry,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(crz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(crxr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cryr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(crzr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(crxrg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cryrg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(crzrg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(clr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr12,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr21,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr13,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr31,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr23,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cr32,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
  END SUBROUTINE soc_matelement

  !==================================================================
  SUBROUTINE do_tsh_isc(c0,c1,c2,sc0,rhoe,psi,eigv,etot_isc,& 
       soc_el,soc_array,infi,nstate,nroot)
    !==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), &
                                                c1(nkpt%ngwk,crge%n,*), &
                                                c2(:,:), sc0(:)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:,:), etot_isc
    COMPLEX(real_8)                          :: soc_el(:)
    REAL(real_8)                             :: soc_array(:)
    INTEGER                                  :: infi, nstate, nroot

    CHARACTER(*), PARAMETER                  :: procedureN = 'do_tsh_isc'

    INTEGER                                  :: i, isoc, isoc1, isoc2, isub, &
                                                j, jsoc, jsoc1, jsoc2
    LOGICAL                                  :: fexist
    REAL(real_8)                             :: minene

    CALL tiset(procedureN,isub)
    IF (tshl%isc) THEN
       IF (md_on_sin) THEN
          DO i=1,nstate
             DO j=1,nroot
                CALL dcopy(2*ncpw%ngw,c1(1,i,j),1,cs1(1,i,j),1)
             ENDDO
          ENDDO
          td01%ns_tri=td01%ns_sin
          td01%ns_sin=0
          lrsym=3
       ELSE
          DO i=1,nstate
             DO j=1,nroot
                CALL dcopy(2*ncpw%ngw,c1(1,i,j),1,ct1(1,i,j),1)
             ENDDO
          ENDDO
          td01%ns_sin=td01%ns_tri
          td01%ns_tri=0
          lrsym=1
       ENDIF
       !
       tshl%isc_2ndlr=.TRUE.
       CALL lr_tddft(c0(:,:),c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,.FALSE.,td01%ioutput)
       tshl%isc_2ndlr=.FALSE.
       !
       IF (md_on_sin) THEN
          DO i=1,nstate
             DO j=1,nroot
                CALL dcopy(2*ncpw%ngw,c1(1,i,j),1,ct1(1,i,j),1)
                CALL dcopy(2*ncpw%ngw,cs1(1,i,j),1,c1(1,i,j),1)
             ENDDO
          ENDDO
       ELSE
          DO i=1,nstate
             DO j=1,nroot
                CALL dcopy(2*ncpw%ngw,c1(1,i,j),1,cs1(1,i,j),1)
                CALL dcopy(2*ncpw%ngw,ct1(1,i,j),1,c1(1,i,j),1)
             ENDDO
          ENDDO
       ENDIF

       !Do SOC calculation
       !goto 111
       IF (MOD(infi,socskip).EQ.0.OR.do_soc_because_sh) THEN
          IF (paral%io_parent) WRITE(6,'(/,A)') &
               ' Start SOC calculation ... '
          CALL zeroing(soc_el)!,3)
          IF (md_on_sin) THEN
             !the reference state is the running state
             cnti%isocs=td01%fstate
             DO jsoc=1,nroot
                soc_array(jsoc)=0._real_8
             ENDDO
             ! find the two closest triplet states to the running state
             minene=10000._real_8
             DO jsoc=1,nroot
                IF (dabs(ene_sin(1,td01%fstate+1)-ene_tri(1,jsoc+1)).LT.minene) THEN
                   minene=dabs(ene_sin(1,td01%fstate+1)-ene_tri(1,jsoc+1))
                   jsoc1=jsoc
                ENDIF
             ENDDO
             minene=10000._real_8
             DO jsoc=1,nroot
                IF (jsoc.EQ.jsoc1) CYCLE
                IF (dabs(ene_sin(1,td01%fstate+1)-ene_tri(1,jsoc+1)).LT.minene) THEN
                   minene=dabs(ene_sin(1,td01%fstate+1)-ene_tri(1,jsoc+1))
                   jsoc2=jsoc
                ENDIF
             ENDDO
             DO jsoc=1,nroot
                IF ((jsoc.EQ.jsoc1).OR.(jsoc.EQ.jsoc2)) THEN
                   IF (paral%io_parent) WRITE(6,'(A,I5,A,I5)') &
                        ' singlet:', td01%fstate, ' with triplet:', jsoc
                   CALL soc_get_rks(c0,cnti%isocs,jsoc,nstate,soc_el)
                   soc_array(jsoc)=0._real_8
                   DO i=1,3
                      soc_array(jsoc)=soc_array(jsoc)+REAL((soc_el(i)*CONJG(soc_el(i))),kind=real_8)
                   ENDDO
                   soc_array(jsoc)=SQRT(soc_array(jsoc))
                ENDIF
             ENDDO
          ELSEIF (md_on_tri) THEN
             !the reference state is the running state
             DO isoc=1,nroot
                soc_array(isoc)=0._real_8
             ENDDO
             ! find the two closest triplet states to the running state
             minene=10000._real_8
             DO isoc=1,nroot
                IF (dabs(ene_tri(1,td01%fstate+1)-ene_sin(1,isoc+1)).LT.minene) THEN
                   minene=dabs(ene_tri(1,td01%fstate+1)-ene_sin(1,isoc+1))
                   isoc1=isoc
                ENDIF
             ENDDO
             minene=10000._real_8
             DO isoc=1,nroot
                IF (isoc.EQ.isoc1) CYCLE
                IF (dabs(ene_tri(1,td01%fstate+1)-ene_sin(1,isoc+1)).LT.minene) THEN
                   minene=dabs(ene_tri(1,td01%fstate+1)-ene_sin(1,isoc+1))
                   isoc2=isoc
                ENDIF
             ENDDO
             cnti%jsoct=td01%fstate
             DO isoc=1,nroot
                IF ((isoc.EQ.isoc1).OR.(isoc.EQ.isoc2)) THEN
                   IF (paral%io_parent) WRITE(6,'(A,I5,A,I5)') &
                        ' triplt:', td01%fstate, ' with singlet:', isoc
                   CALL soc_get_rks(c0,isoc,cnti%jsoct,nstate,soc_el)
                   soc_array(isoc)=0._real_8
                   DO i=1,3
                      soc_array(isoc)=soc_array(isoc)+REAL((soc_el(i)*CONJG(soc_el(i))),kind=real_8)
                   ENDDO
                   soc_array(isoc)=SQRT(soc_array(isoc))
                ENDIF
             ENDDO
          ENDIF
          IF (paral%io_parent) WRITE(6,'(A,/)') &
               ' End SOC calculation'
          !
          !write soc
          IF (md_on_sin) THEN
             IF (paral%io_parent) THEN
                WRITE(6,'(/,1x,"ISC:",60("-"))')
                WRITE(6,'(A)') ' SOCs with the triplets(cm-1 ):'
                WRITE(6,'(100F12.6,/)') (soc_array(jsoc)*facthcm,jsoc=1,nroot)
                WRITE(6,'(1x,"ISC:",60("-"),/)')
             ENDIF
          ELSEIF (md_on_tri) THEN
             IF (paral%io_parent) THEN
                WRITE(6,'(/,1x,"ISC:",60("-"))')
                WRITE(6,'(A)') ' SOCs with the singlets (cm-1):'
                WRITE(6,'(100F12.6,/)') (soc_array(isoc)*facthcm,isoc=1,nroot)
                WRITE(6,'(1x,"ISC:",60("-"),/)')
             ENDIF
          ENDIF
          fexist=.FALSE.
          IF (paral%io_parent)&
               INQUIRE(file="ISC_SOC.dat",exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent)&
                  OPEN(unit=50,file="ISC_SOC.dat",status='NEW')
          ELSE
             IF (paral%io_parent)&
                  OPEN(unit=50,file="ISC_SOC.dat",status='UNKNOWN',position='APPEND')
          ENDIF
          IF (paral%io_parent) THEN
             IF (md_on_sin) THEN
                WRITE(50,'(3I6,100F12.6,/)')&
                     infi,0,td01%fstate,(soc_array(isoc)*facthcm,isoc=1,nroot)
                !   step,mult,state,(SOCs(i),i=1,all other states with different multiplicity)
             ELSEIF (md_on_tri) THEN
                WRITE(50,'(3I6,100F12.6,/)')&
                     infi,3,td01%fstate,(soc_array(isoc)*facthcm,isoc=1,nroot)
                !   step,mult,state,(SOCs(i),i=1,all other states with different multiplicity)
             ENDIF
          ENDIF
          IF (paral%io_parent) CLOSE(50)
       ENDIF
       !111  continue
       !
       IF (md_on_sin) THEN
          td01%ns_sin=td01%ns_tri
          td01%ns_tri=0
          lrsym=1
       ELSE
          td01%ns_tri=td01%ns_sin
          td01%ns_sin=0
          lrsym=3
       ENDIF
    ENDIF


    CALL tihalt(procedureN,isub)
  END SUBROUTINE do_tsh_isc

  !==================================================================
  SUBROUTINE do_tsh_isc_lz(soc_array,dt,nroot,nskip)
    !==--------------------------------------------------------------==
    !-- Arguments:
    REAL(real_8)                             :: soc_array(:), dt
    INTEGER                                  :: nroot, nskip

    CHARACTER(*), PARAMETER                  :: procedureN = 'do_tsh_isc_lz'

    INTEGER                                  :: i, isub, j
    INTEGER, SAVE                            :: ifirst = 0, iswitch = 0
    LOGICAL                                  :: fexist, isc_constr, switch
    REAL(real_8)                             :: de1, de2, dedt, &
                                                delta_ene(2,nroot), massey, &
                                                prob, prod_delta, test

    CALL tiset(procedureN,isub)
    ifirst=ifirst+1
    ! if (paral%io_parent) write(6,'(A,I5,/)') 'ifirst',ifirst 
    !---------------------------------------------------------------

    switch=.FALSE.
    ! running state: td01%fstate

    IF ((ene_tri(2,td01%fstate+1).EQ.0._real_8) .OR. &
         (ene_sin(2,td01%fstate+1).EQ.0._real_8) ) THEN
       !---------------------------------------------------------------
       ! not enough history to do derivatives
       IF (paral%io_parent) WRITE(6,'(/,A,I5)') &
            ' fstate in do_tsh_isc_lz', td01%fstate
       IF (paral%io_parent) WRITE(6,'(10F8.3)') (ene_sin(1,i),i=1,nroot)
       IF (paral%io_parent) WRITE(6,'(10F8.3)') (ene_sin(2,i),i=1,nroot)
       IF (paral%io_parent) WRITE(6,'(10F8.3)') (ene_tri(1,i),i=1,nroot)
       IF (paral%io_parent) WRITE(6,'(10F8.3)') (ene_tri(2,i),i=1,nroot)
       !----------------------------------------------------------------
       CALL tihalt(procedureN,isub)
       RETURN  
    ENDIF




    IF (paral%io_parent) THEN
       WRITE(101,*) '###################'
       WRITE(102,*) '###################'
       WRITE(6,'(A)') &
            ' ISC Landau Zener test:    state    LZ prob   random no'
       DO i=1,nroot
          IF (md_on_sin) THEN
             ! test to see if it is a crossing point
             delta_ene(1,i)=ene_sin(2,td01%fstate+1)-ene_tri(2,i+1)
             delta_ene(2,i)=ene_sin(1,td01%fstate+1)-ene_tri(1,i+1)
             prod_delta=delta_ene(1,i)*delta_ene(2,i)  
             IF(prod_delta.GE.0.0_real_8)THEN
                isc_constr=.FALSE.
             ELSEIF(prod_delta.LT.0.0_real_8)THEN
                isc_constr=.TRUE.
             ENDIF

             de2=dabs(ene_sin(2,td01%fstate+1)-ene_tri(2,i+1))
             de1=dabs(ene_sin(1,td01%fstate+1)-ene_tri(1,i+1))
             dedt=dabs(de2-de1)/(dt*nskip)
             massey=(dabs(soc_array(i))**2)/(dedt)
             prob=1.0_real_8-EXP(-2._real_8*pi*massey)
             test=repprngu()
             !---------------------------------------------------------------
             IF (paral%io_parent) WRITE(6,'(27X,I5,4X,F7.2,5X,F7.2)') i,prob,test
             !---------------------------------------------------------------
             IF (iswitch.LT.(socskip+2))      isc_constr=.FALSE.
             IF (test.LT.prob.AND.isc_constr) THEN
                switch=.TRUE.
                td01%fstate=i
                iswitch=0
                GOTO 111
             ENDIF
          ELSE
             ! test to see if it is a crossing point
             delta_ene(1,i)=ene_tri(2,td01%fstate+1)-ene_sin(2,i+1)
             delta_ene(2,i)=ene_tri(1,td01%fstate+1)-ene_sin(1,i+1)
             prod_delta=delta_ene(1,i)*delta_ene(2,i)  
             IF(prod_delta.GE.0.0_real_8)THEN
                isc_constr=.FALSE.
             ELSEIF(prod_delta.LT.0.0_real_8)THEN
                isc_constr=.TRUE.
             ENDIF

             de2=dabs(ene_tri(2,td01%fstate+1)-ene_sin(2,i+1))
             de1=dabs(ene_tri(1,td01%fstate+1)-ene_sin(1,i+1))
             dedt=dabs(de2-de1)/(dt*nskip)
             massey=(dabs(soc_array(i))**2)/(dedt)
             prob=1.0_real_8-EXP(-2._real_8*pi*massey)
             test=repprngu()
             !---------------------------------------------------------------
             IF (paral%io_parent) WRITE(6,'(27X,I5,4X,F7.2,5X,F7.2)') i,prob,test
             !---------------------------------------------------------------
             IF (iswitch.LT.(socskip+2))      isc_constr=.FALSE.
             IF (test.LT.prob.AND.isc_constr) THEN
                switch=.TRUE.
                td01%fstate=i
                iswitch=0
                GOTO 111
             ENDIF
          ENDIF
          WRITE(101,'(I5,I2,L2)') tshi%shstep,i,isc_constr
          WRITE(102,'(I5,I2,3F15.9)') tshi%shstep,i,delta_ene(1,i),delta_ene(2,i),prod_delta
       ENDDO
       WRITE(6,*) ' '
    ENDIF
    iswitch=iswitch+1
111 CONTINUE
    !bcast decision
    CALL mp_bcast(switch,parai%source,parai%allgrp)
    IF (switch) THEN
       CALL mp_bcast(td01%fstate,parai%source,parai%allgrp)
       CALL mp_bcast(iswitch,parai%source,parai%allgrp)
    ENDIF
    IF (paral%io_parent.AND.switch) THEN
       IF (md_on_sin) THEN
          WRITE(6,'(/,A)') ' ISC from singlet to triplet'
          WRITE(6,'(A,I6,/)') ' NEW state: triplet ',td01%fstate
       ELSEIF (md_on_tri) THEN
          WRITE(6,'(/,A)') ' ISC from triplet to singlet'
          WRITE(6,'(A,I6,/)') ' NEW state: singlet ',td01%fstate
       ENDIF
    ENDIF
    ! 
    IF (md_on_sin.AND.switch) THEN
       md_on_sin=.FALSE.
       md_on_tri=.TRUE.
       td01%ns_tri=td01%ns_sin
       td01%ns_sin=0
       lrsym=3
       CALL zeroing(c_sh)!,2*(nroot+1))
       c_sh(1,td01%fstate+1)=(1._real_8,0._real_8)
       c_sh(2,td01%fstate+1)=(1._real_8,0._real_8)
    ELSEIF(md_on_tri.AND.switch) THEN
       md_on_tri=.FALSE.
       md_on_sin=.TRUE.
       td01%ns_sin=td01%ns_tri
       td01%ns_tri=0
       lrsym=1
       CALL zeroing(c_sh)!,2*(nroot+1))
       c_sh(1,td01%fstate+1)=(1._real_8,0._real_8)
       c_sh(2,td01%fstate+1)=(1._real_8,0._real_8)
    ENDIF

    fexist=.FALSE.
    IF (paral%io_parent)&
         INQUIRE(file="ISC_ENERGIES.dat",exist=fexist)
    IF (.NOT.fexist) THEN
       IF (paral%io_parent)&
            OPEN(unit=50,file="ISC_ENERGIES.dat",status='NEW')
    ELSE
       IF (paral%io_parent)&
            OPEN(unit=50,file="ISC_ENERGIES.dat",status='UNKNOWN',position='APPEND')
    ENDIF

    IF (md_on_sin) THEN
       IF (paral%io_parent)&
            WRITE(50,'(i10,40f15.9)') tshi%shstep,&
            (etot_isc_md+ene_sin(1,j),j=1,nroot+1),&
            (etot_isc_md+ene_tri(1,j),j=1,nroot+1),&
            etot_isc_md+ene_sin(1,td01%fstate+1)
    ELSEIF (md_on_tri) THEN
       IF (paral%io_parent)&
            WRITE(50,'(i10,40f15.9)') tshi%shstep,&
            (etot_isc_md+ene_sin(1,j),j=1,nroot+1),&
            (etot_isc_md+ene_tri(1,j),j=1,nroot+1),&
            etot_isc_md+ene_tri(1,td01%fstate+1)
    ENDIF
    IF (paral%io_parent) CLOSE(50)

    fexist=.FALSE.
    IF (paral%io_parent)&
         INQUIRE(file="ISC_ENERGIES_ex.dat",exist=fexist)
    IF (.NOT.fexist) THEN
       IF (paral%io_parent)&
            OPEN(unit=50,file="ISC_ENERGIES_ex.dat",status='NEW')
    ELSE
       IF (paral%io_parent)&
            OPEN(unit=50,file="ISC_ENERGIES_ex.dat",status='UNKNOWN',position='APPEND')
    ENDIF
    IF (paral%io_parent.AND.md_on_sin)&
         WRITE(50,'(i10,40f15.9)') tshi%shstep,&
         (ene_sin(1,j),j=2,nroot+1),&
         (ene_tri(1,j),j=2,nroot+1),ene_sin(1,td01%fstate+1)
    IF (paral%io_parent.AND.md_on_tri)&
         WRITE(50,'(i10,40f15.9)') tshi%shstep,&
         (ene_sin(1,j),j=2,nroot+1),&
         (ene_tri(1,j),j=2,nroot+1),ene_tri(1,td01%fstate+1)
    IF (paral%io_parent) CLOSE(50)

    CALL tihalt(procedureN,isub)
  END SUBROUTINE do_tsh_isc_lz

  !==================================================================
  SUBROUTINE isc_read_start
    !==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'isc_read_start'

    INTEGER                                  :: isub, multipl
    LOGICAL                                  :: fexist

    CALL tiset(procedureN,isub)
    fexist=.FALSE.

    IF (paral%io_parent) THEN
       INQUIRE(file="ISC_RESTART_STATE.dat",exist=fexist)
       IF (.NOT.fexist) THEN
          !multiplicity and running state read from input
          OPEN(unit=50,file="ISC_RESTART_STATE.dat",status='new')
          IF (td01%ns_sin.NE.0) THEN
             WRITE(50,'(I2,2X,I5)')  0,td01%fstate
          ELSE
             WRITE(50,'(I2,2X,I5)')  3,td01%fstate
          ENDIF
       ELSE
          OPEN(unit=50,file="ISC_RESTART_STATE.dat",status='old')
          READ(50,*) multipl,td01%fstate
          IF (multipl.EQ.0) THEN
             md_on_sin=.TRUE.
             md_on_tri=.FALSE. 
             td01%ns_sin=MAX(td01%ns_sin,td01%ns_tri)
             td01%ns_tri=0
          ELSEIF(multipl.EQ.3) THEN
             md_on_sin=.FALSE.
             md_on_tri=.TRUE.
             td01%ns_tri=MAX(td01%ns_sin,td01%ns_tri)
             td01%ns_sin=0
          ELSE
             CALL stopgm(proceduren,'error reading ISC_RESTART_STATE.dat',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       CLOSE(50)
    ENDIF
    !
    CALL mp_bcast(fexist,parai%source,parai%allgrp)
    !
    IF (fexist) THEN
       CALL mp_bcast(td01%ns_sin,parai%source,parai%allgrp) 
       CALL mp_bcast(td01%ns_tri,parai%source,parai%allgrp)
       CALL mp_bcast(md_on_sin,parai%source,parai%allgrp)
       CALL mp_bcast(md_on_tri,parai%source,parai%allgrp)
       CALL mp_bcast(td01%fstate,parai%source,parai%allgrp)
    ENDIF
    CALL tihalt(procedureN,isub)
  END SUBROUTINE isc_read_start

  !==================================================================
  SUBROUTINE isc_write_start
    !==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'isc_write_start'

    IF (paral%io_parent) THEN
       OPEN(unit=50,file="ISC_RESTART_STATE.dat",status='replace')
       IF (td01%ns_sin.NE.0) THEN
          WRITE(50,'(I2,2X,I5)')  0,td01%fstate
       ELSE
          WRITE(50,'(I2,2X,I5)')  3,td01%fstate
       ENDIF
       CLOSE(50)
    ENDIF
    RETURN
  END SUBROUTINE isc_write_start
  !==================================================================


  SUBROUTINE soc_real_space_mult(cin,pos,k,cout)
    !==--------------------------------------------------------------==
    !-- Arguments:
    COMPLEX(real_8)                          :: cin(ncpw%ngw)
    REAL(real_8)                             :: pos(*)
    INTEGER                                  :: k
    COMPLEX(real_8)                          :: cout(ncpw%ngw)

    CHARACTER(*), PARAMETER :: procedureN = 'soc_real_space_mult'

    COMPLEX(real_8), ALLOCATABLE             :: psi(:)
    INTEGER                                  :: ierr, ig, ir, isub
    LOGICAL                                  :: do_phase
    REAL(real_8), ALLOCATABLE                :: dest(:), dest_real(:), &
                                                real_psi(:)

    CALL tiset(procedureN,isub)
    ALLOCATE(dest(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dest_real(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(real_psi(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    !=----------------------------------------------------------------=
    do_phase=.FALSE.

    CALL zeroing(psi)!,maxfft)
    !
    DO ig=1,ncpw%ngw
       psi(nzhs(ig))  = cin(ig)*CMPLX(0.0_real_8,1.0_real_8)
       psi(indzs(ig)) = CONJG(cin(ig)*CMPLX(0.0_real_8,1.0_real_8))
    ENDDO
    CALL invfftn(psi,.TRUE.,parai%allgrp)
    IF (do_phase) CALL phase(psi)

    DO ir=1,fpar%nnr1
       psi(ir)=psi(ir)/SQRT(parm%omega)
       real_psi(ir)=REAL(psi(ir),kind=real_8)
       ! real_psi(ir)=dimag(psi(ir))
    ENDDO

    CALL soc_mult_psi_pos(pos,k,real_psi,dest_real)

    DO ir=1,fpar%nnr1
       dest_real(ir)=dest_real(ir)*SQRT(parm%omega)
    ENDDO

    CALL zeroing(psi)!,maxfft)
    !parallel do private(ir)
    DO ir=1,fpar%nnr1
       psi(ir)  = CMPLX(dest_real(ir),0.0_real_8)
    ENDDO

    IF (do_phase) CALL phase(psi)
    CALL fwfftn(psi,.FALSE.,parai%allgrp)
    !parallel do private(ig)
    DO ig=1,ncpw%ngw
       cout(ig)  = psi(nzhs(ig))
    ENDDO
    !
    DEALLOCATE(dest,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dest_real,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(real_psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
  END SUBROUTINE soc_real_space_mult
  !================================================================== 
  SUBROUTINE soc_mult_psi_pos(pos,k,psi,dest)
    !==================================================================
    REAL(real_8) :: pos(fpar%kr1,fpar%kr2,fpar%kr3,3)
    INTEGER                                  :: k
    REAL(real_8) :: psi(fpar%kr1,fpar%kr2,fpar%kr3), &
      dest(fpar%kr1,fpar%kr2,fpar%kr3)

    CHARACTER(*), PARAMETER :: procedureN = 'soc_mult_psi_pos'

    INTEGER                                  :: ix, iy, iz

    DO iz=1,parm%nr3
       DO iy=1,parm%nr2
          !$OMP parallel do private(ix)
          DO ix=1,parm%nr1
             dest(ix,iy,iz)=psi(ix,iy,iz)*pos(ix,iy,iz,k)
          ENDDO
       ENDDO
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_mult_psi_pos
  !==================================================================
  SUBROUTINE soc_mod_pos(pos,vect)
    !==--------------------------------------------------------------==
    !Arguments --
    REAL(real_8) :: pos(fpar%kr1,fpar%kr2,fpar%kr3,3), vect(3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_mod_pos'

    INTEGER                                  :: i, ix, iy, iz
    REAL(real_8)                             :: min_delta, r1, r3, xtol

    xtol=1.0e-16_real_8 ! machine precision
    !
    min_delta=1000.0
    min_delta=MIN(parm%a1(1)/spar%nr1s,min_delta)
    min_delta=MIN(parm%a2(2)/spar%nr2s,min_delta)
    min_delta=MIN(parm%a3(3)/spar%nr3s,min_delta)
    !
    DO iz=1,parm%nr3
       DO iy=1,parm%nr2
          DO ix=1,parm%nr1
             r1=SQRT((vect(1)-pos(ix,iy,iz,1))**2+ &
                  (vect(2)-pos(ix,iy,iz,2))**2+ &
                  (vect(3)-pos(ix,iy,iz,3))**2)
             r3=r1**3
             !      if (r1.le.(min_delta/20._real_8)) then
             IF (r3.LE.xtol) THEN
                DO i=1,3
                   pos(ix,iy,iz,i)=0._real_8
                ENDDO
                WRITE(*,*)'one point excluded from the integral'
             ELSE
                DO i=1,3
                   pos(ix,iy,iz,i)=(pos(ix,iy,iz,i)-vect(i))/r3
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_mod_pos
  !==================================================================
  SUBROUTINE soc_pos(taus,pos)
    !==--------------------------------------------------------------==
    REAL(real_8) :: taus(:,:,:), pos(fpar%kr1,fpar%kr2,fpar%kr3,3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_pos'

    INTEGER                                  :: i, ia, is, ix, iy, iz, nx0
    REAL(real_8)                             :: dx(3,3), rr(3), x(3), zz

    DO i=1,3
       rr(i)=0.0_real_8
    ENDDO
    zz=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO i=1,3
             rr(i)=rr(i)+taus(i,ia,is)*ions0%zv(is)
          ENDDO
          zz=zz+ions0%zv(is)
       ENDDO
    ENDDO
    DO i=1,3
       rr(i)=rr(i)/zz
       socvar%rr_soc(i)=rr(i)
    ENDDO
    !if(paral%io_parent) write (6,'(/,a,3f9.5,/)')
    !     ' center of integration (core charge): ', (rr(i), i=1,3)
    !Calculate the dipole moment in real space
    DO i=1,3
       dx(i,1)=parm%a1(i)/spar%nr1s
       dx(i,2)=parm%a2(i)/spar%nr2s
       dx(i,3)=parm%a3(i)/spar%nr3s
    ENDDO
    nx0=parap%nrxpl(parai%mepos,1)
    DO iz=1,parm%nr3
       DO iy=1,parm%nr2
          DO i=1,3
             x(i)=(iz-1)*dx(i,3) + (iy-1)*dx(i,2) + (nx0-1)*dx(i,1)
          ENDDO
          DO ix=1,parm%nr1
             DO i=1,3
                pos(ix,iy,iz,i)=x(i)-rr(i)
                x(i)=x(i)+dx(i,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_pos
  !==================================================================
  SUBROUTINE soc_dotp(c0,c1,res)
    !==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw), c1(ncpw%ngw)
    REAL(real_8)                             :: res

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_dotp'

    COMPLEX(real_8)                          :: caux
    INTEGER                                  :: ig

    res = 0._real_8
    DO ig=1,ncpw%ngw
       caux=CONJG(c0(ig))*c1(ig)
       res=res+REAL(caux,kind=real_8)
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_dotp
  !==================================================================
  !     SUBROUTINES FOR TESTS
  !==================================================================
  SUBROUTINE soc_dodipole(rhoe,pos,dipole)
    !==--------------------------------------------------------------==
    !Arguments --
    REAL(real_8) :: rhoe(fpar%kr1,fpar%kr2,fpar%kr3), &
      pos(fpar%kr1,fpar%kr2,fpar%kr3,3), dipole(3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_dodipole'

    INTEGER                                  :: i, ix, iy, iz
    REAL(real_8)                             :: rho, scl

    scl=parm%omega/DBLE(spar%nr1s*spar%nr2s*spar%nr3s)
    !
    CALL zeroing(dipole)!,3)
    DO iz=1,parm%nr3
       DO iy=1,parm%nr2
          DO ix=1,parm%nr1
             rho=-rhoe(ix,iy,iz)
             DO i=1,3
                dipole(i)=dipole(i) + pos(ix,iy,iz,i)*rho
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL dscal(3,scl,dipole,1)
    CALL mp_sum(dipole,3,parai%allgrp)
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_dodipole
  !================================================================== 
  SUBROUTINE soc_dipole_real(cin,taus,nstate)
    !==--------------------------------------------------------------==
    !-- Arguments:
    COMPLEX(real_8)                          :: cin(ncpw%ngw,*)
    REAL(real_8)                             :: taus(:,:,:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'soc_dipole_real'

    COMPLEX(real_8), ALLOCATABLE             :: psi(:)
    INTEGER                                  :: i, ierr, ig, ir
    REAL(real_8)                             :: dipole(3), norm_l
    REAL(real_8), ALLOCATABLE                :: dest(:), pos(:,:), rhoe(:)

    ALLOCATE(pos(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dest(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rhoe(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !=----------------------------------------------------------------=
    CALL zeroing(rhoe)!,nnr1)
    CALL soc_pos(taus,pos)

    !if (paral%io_parent) write(6,*) 'nstate ', nstate
    !
    DO i=1,nstate
       CALL zeroing(psi)!,maxfft)
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))  = cin(ig,i)
          psi(indzs(ig)) = CONJG(cin(ig,i))
       ENDDO
       CALL invfftn(psi,.TRUE.,parai%allgrp)
       !
       DO ir=1,fpar%nnr1
          psi(ir)=psi(ir)/SQRT(parm%omega)
       ENDDO
       !
       DO ir=1,fpar%nnr1
          rhoe(ir)  = rhoe(ir)+ REAL(psi(ir)*CONJG(psi(ir)),kind=real_8)
       ENDDO
    ENDDO

    norm_l=0._real_8
    DO ir=1,fpar%nnr1
       norm_l=norm_l+rhoe(ir)* parm%omega/DBLE(spar%nr1s*spar%nr2s*spar%nr3s)
    ENDDO
    CALL mp_sum(norm_l,parai%allgrp)

    CALL soc_dodipole(rhoe,pos,dipole)

    IF (paral%io_parent) WRITE(6,'(A10,4(1X,E12.6))') 'dipole ',&
         norm_l,dipole(1),dipole(2),dipole(3)

    DEALLOCATE(pos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dest,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__) 
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE soc_dipole_real
  !==================================================================

END MODULE soc

