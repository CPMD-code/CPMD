MODULE g_loc_dr_utils
  USE cnst,                            ONLY: fbohr
  USE cppt,                            ONLY: gk
  USE ddipo_utils,                     ONLY: set_operator
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: &
       filcen, filmat, filpen, filspr, g2g_mem, gifirst, glocal, gloci, &
       glocr, ind_st, indstate, ist_list, lostate, shiftkx, shiftky, shiftkz, &
       wan_mem
  USE g_loc_exp_sum_utils,             ONLY: g_loc_exp_sum
  USE g_loc_optim_utils,               ONLY: g_loc_exp_ide,&
                                             g_loc_xyzmat
  USE g_loc_realspace_utils,           ONLY: p_real_space
  USE g_loc_spread_ide_utils,          ONLY: g_loc_spread_ide
  USE g_loc_spread_sum_utils,          ONLY: g_loc_spread_sum
  USE g_loc_util_utils,                ONLY: r_matrix
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE setirec_utils,                   ONLY: write_irec
  USE system,                          ONLY: fpar,&
                                             mapgp,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE utils,                           ONLY: zclean_k
  USE wann,                            ONLY: wannc
  USE wannier_center_utils,            ONLY: wannier_center
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
!!use rotate_utils, only : rotate_c
!!use wfn_print_utils, only : wfnz_print

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: g_loc_dr_real
  PUBLIC :: g_loc_dr_comp

CONTAINS

  ! ==================================================================
  SUBROUTINE g_loc_dr_real(c0,c2,nstate,eigv,tau0,taup,nfi)
    ! ==--------------------------------------------------------------==
    ! == PREPARE ENVIRONMENT FOR LOCALIZATION IN G SPACE              ==
    ! ==--------------------------------------------------------------==
    ! input
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(2*ncpw%ngw,nstate)
    REAL(real_8)                             :: eigv(*), tau0(:,:,:), &
                                                taup(:,:,:)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_loc_dr_real'

    CHARACTER(len=15)                        :: filorb, filtr
    CHARACTER(len=8)                         :: numb
    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: c1(:,:), psi(:), rotmat(:,:), &
                                                xyzmat(:,:,:)
    INTEGER                                  :: cont1, cont2, i, ia, ie, &
                                                ierr, ig, ikind, il_psi, is, &
                                                js
    LOGICAL                                  :: ferror, test
    REAL(real_8)                             :: im, mod, mod2, over_sup, re
    REAL(real_8), ALLOCATABLE                :: cenr(:,:), center(:,:)

! TEST
! ==--------------------------------------------------------------==

    test = .FALSE.
    IF (tkpts%tkpnt)  CALL stopgm('G_LOC_DR',&
         'KPOINT AND REAL WF INCONSISTENT',& 
         __LINE__,__FILE__)

    ALLOCATE(c1(2*ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1)!,2*ngw*nstate)

    DO is = 1,nstate
       DO ig = 1,ncpw%ngw
          re  = REAL(c0(ig,is))
          im  = AIMAG(c0(ig,is))
          c1(ig    ,is)  = CMPLX(re, im,kind=real_8)
          c1(ig+ncpw%ngw,is)  = CMPLX(re,-im,kind=real_8)
       ENDDO
    ENDDO

    IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)


    nkpt%ngwk  = 2*ncpw%ngw
    spar%ngwks = 2*spar%ngws

    IF (test) THEN
       DO is = 1,nstate
          filtr='ORB_in_'
          IF (paral%io_parent)&
               WRITE(numb,'(I8)') is
          CALL xstring(numb,ia,ie)
          filorb=filtr(1:7)//numb(ia:ie)//'.dat'
          IF (paral%io_parent)&
               CALL fileopen(17,filorb,fo_def,ferror)
          DO ig = nkpt%ngwk,ncpw%ngw,-1
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
             IF (paral%io_parent)&
                  WRITE(17,'(I8,2(2x,1PE12.5))') -ig+ncpw%ngw,mod,LOG10(mod)
          ENDDO
          DO ig = 1,ncpw%ngw
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
             IF (paral%io_parent)&
                  WRITE(17,'(I8,2(2x,1PE12.5))') ig,mod,LOG10(mod)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(17)
       ENDDO
    ENDIF

    IF (glocal%tglocprint) THEN
       ! print selected localized functions
       il_psi = fpar%kr1s*fpar%kr2s*fpar%kr3s
       ALLOCATE(psi(il_psi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL wfnz_print(1,1,c1,nstate,psi,nstate)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (glocal%tglocrealp) THEN
       ! print selected localized functions
       il_psi = fpar%kr1*fpar%kr2s*fpar%kr3s
       ALLOCATE(psi(il_psi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL p_real_space(1,1,c1,nstate,psi,nstate)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF


    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64(1H*))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)')&
            'LOCALIZATION OF WAVEFUNCTION IN RECIPROCAL SPACE '
       IF (gloci%gloc_type .EQ. 1) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       FUNCTIONAL:  Sum_(n,i)[ <n|Gi^2|n> - <n|Gi|n>^2 ]'
       ELSEIF (gloci%gloc_type .EQ. 2) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       FUNCTIONAL:  Sum_(n,i)[1 - |<n| EXP(iGir_o) |n>|^2]'
       ENDIF
       IF ((gloci%gloc_opt.EQ.1).AND.paral%io_parent)&
            WRITE(6,'(A)')&
            '       OPTIMISATION: U UPDATE BY STEEPEST DESCENT'
       IF ((gloci%gloc_opt.EQ.2) .AND.paral%io_parent)&
            WRITE(6,'(A,/,A)')&
            '       OPTIMISATION: U NEAR TO 1, WFN UPDATE EACH STEP',&
            '       SMALL ROTATIONS ASSUMED'
       IF (gloci%gloc_const .EQ. 1) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,/,A,1PE10.3,A,I10)')&
               '       UNITARITY  BY LAGRANGIAN MULTIPLIERS',&
               '       TOLERANCE = ',glocr%gepslag,' MAX n ITERATION = ',&
               gloci%gloc_maxit
       ELSEIF (gloci%gloc_const .EQ. 2) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       UNITARITY  BY ORTHOGONALIZATION [S^(-1/2)]'
       ELSEIF (gloci%gloc_const .EQ. 3) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       UNITARITY  APPROXIMATION U = 1+iA (A HERMITIAN)'
       ENDIF
       IF (glocal%tg_linesearch) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               '       TIME STEP TUNING BY SEARCH ALONG THE LINE'
       ENDIF
       IF (glocal%tg_read_matrix) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               '       FIRST UNITARY TRANSFORMATION READ FROM FILE'
       ENDIF
       IF (glocal%tgwannier) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/A,A)')&
               '       THE SPREAD FUNCTIONAL IN REAL',&
               ' SPACE IS ALSO CALCULATED'
          IF (paral%io_parent)&
               WRITE(6,'(A,f6.3/A,f6.3)')&
               '            optimization: W_GFUN = ',glocr%g2g_weight,&
               '                        : W_RFUN = ',glocr%wan_weight
          IF (gloci%gloc_opt.NE.1 .AND. gloci%gloc_type .NE. 2)&
               CALL stopgm('G_LOC_DR',&
               'REAL SPACE LOC. ONLY WITH GLOC_OPT=1, GLOC_TYPE=2',& 
               __LINE__,__FILE__)
       ENDIF
       IF (glocal%tg_antisymm) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,A)')&
               'PENALTY FUNCTION FOR SYMMETRIC WFN CALCULATED',&
               '  FOR FIRST 100 STEPS'
          IF ((glocal%tg_penalty).AND.paral%io_parent)&
               WRITE(6,'(A,/)')&
               '     PENALTY IMPOSED IN THE GRADIENT (ADDED TERM)'
       ENDIF
       ! IF(TCGRAD) THEN
       ! IF (GLOC_OPT.EQ.2) CALL STOPGM('G_LOC_DR_REAL',
       ! &          'CONJUGATE GRADIENT and U ALMOST 1 NOT POSSIBLE')
       ! WRITE(6,'(A)') '       CONJUGATE GRADIENT USED'
       ! ENDIF

       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         CONVERGENCE CRITERIA:',glocr%gloc_eps
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            '         MAXIMUM # OF STEPS:',gloci%gloc_maxs
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         RANDOMIZATION AMPLITUDE:',glocr%gloc_ran
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         STEP SIZE:',glocr%gloc_step
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            '         INITIAL STEP',gloci%gloc_init+1
       IF (glocal%tglocprint) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1x,A)')&
               'STARTING  AND OPTIMIZED WFN PRINTED IN SEPARETED FILES'
          IF (gloci%gloc_all .EQ. 1) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>>  all wfn printed  <<<")')
          ELSEIF (gloci%gloc_first .GT. 0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>> wfn from ",I5," to ",I6," printed  <<<<")')
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(">>>>  printed wfn list number = ",I6)') gloci%gloc_orb
          ENDIF
       ENDIF
    ENDIF

    gifirst = 1
    ikind = 0
    wan_mem = glocr%wan_weight
    g2g_mem = glocr%g2g_weight

    filtr='SPR_'
    IF (paral%io_parent)&
         WRITE(numb,'(I8)') ikind
    CALL xstring(numb,ia,ie)
    filspr=filtr(1:4)//numb(ia:ie)//'.dat'
    filtr='CEN_'
    IF (paral%io_parent)&
         WRITE(numb,'(I8)') ikind
    CALL xstring(numb,ia,ie)
    filcen=filtr(1:4)//numb(ia:ie)//'.xyz'
    filtr='PENALTY_'
    IF (paral%io_parent)&
         WRITE(numb,'(I8)') ikind
    CALL xstring(numb,ia,ie)
    filpen=filtr(1:8)//numb(ia:ie)//'.dat'
    filtr='MATRIX_'
    IF (paral%io_parent)&
         WRITE(numb,'(I8)') ikind
    CALL xstring(numb,ia,ie)
    filmat=filtr(1:7)//numb(ia:ie)

    over_sup = 0.0_real_8
    DO is = 1,nstate
       DO ig = 1,nkpt%ngwk
          mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
               AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
          DO js = is+1,nstate
             mod2= REAL(c1(ig,js))*REAL(c1(ig,js))+&
                  AIMAG(c1(ig,js))*AIMAG(c1(ig,js))
             over_sup = over_sup +  mod*mod2
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(over_sup,parai%allgrp)
    IF (paral%io_parent)&
         WRITE(6,'(A,2(1PE12.4))') '   >>>>>>> INITIAL SUPPORT OVERLAP ',&
         over_sup,over_sup*2.0_real_8/(nstate*(nstate-1))


    IF (test) THEN
       ! TEST -----------------------
       ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL r_matrix(rotmat,nstate)
       zone = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
       zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
       CALL rotate_c(zone,c1,zzero,c2,rotmat,nstate)
       CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
       IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)
       CALL set_operator(.FALSE.)
       ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cenr(3,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL g_loc_xyzmat(c1,c2,xyzmat,nstate,tau0)

       ALLOCATE(center(4,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
       cont1 = 1
       DO i = 1,nstate
          cenr(1,i) = center(1,i)
          cenr(2,i) = center(2,i)
          cenr(3,i) = center(3,i)
          IF (cenr(3,i) .LT. -0.5_real_8*parm%a3(3)) THEN
             cenr(3,i) = cenr(3,i) + parm%a3(3)
          ENDIF
          cont1 = cont1+4
          ! write(6,'(I5,3(2x,f10.5))') I,CENR(1,I),CENR(2,I),CENR(3,I)
       ENDDO
       DEALLOCATE(center,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       cont1 = 0
       cont2 = 0
       DO i = 1,nstate
          IF (ABS(cenr(3,i)) .GT. 3.0_real_8) THEN
             CALL dcopy(2*2*ncpw%ngw,c1(1,i),1,c2(1,nstate-cont2),1)
             ! write(6,*) I , ' in ', NSTATE-CONT2,CENR(3,I)
             cont2 = cont2 + 1
          ELSE
             cont1 = cont1 + 1
             CALL dcopy(2*2*ncpw%ngw,c1(1,i),1,c2(1,cont1),1)
             ! write(6,*) I , ' in ', CONT1,CENR(3,I)
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) cont1,cont2, cont1+cont2

       CALL zeroing(c1)!,2*ngw*nstate)
       CALL dcopy(2*2*ncpw%ngw*nstate,c2(1,1),1,c1(1,1),1)
       CALL zeroing(c2)!,2*ngw*nstate)


       DEALLOCATE(xyzmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cenr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rotmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       ! STOP 'CENTERS'
    ENDIF
    ! END TEST --------------------




    IF (gloci%gloc_type .EQ. 1)  THEN

       IF (gloci%gloc_opt .EQ. 1) THEN
          CALL g_loc_spread_sum&
               (c1,c2,nstate,nfi)
       ELSEIF (gloci%gloc_opt .EQ. 2) THEN
          CALL g_loc_spread_ide&
               (c1,c2,nstate,nfi)
       ELSE
          CALL stopgm('G_OC_DR_REAL',&
               'OPTIMIZATION ALGORITHM NOT FOUND',& 
               __LINE__,__FILE__)
       ENDIF

    ELSEIF (gloci%gloc_type .EQ.2) THEN

       IF (gloci%gloc_opt .EQ. 1) THEN
          IF (test) THEN
             CALL g_loc_exp_sum&
                  (c1(1,1),c2,nstate/2,nfi,tau0)
          ELSE
             CALL g_loc_exp_sum&
                  (c1,c2,nstate,nfi,tau0)
          ENDIF
       ELSEIF (gloci%gloc_opt .EQ. 2) THEN
          CALL g_loc_exp_ide&
               (c1,c2,nstate,nfi)
       ELSE
          CALL stopgm('G_OC_DR_REAL',&
               'OPTIMIZATION ALGORITHM NOT FOUND',& 
               __LINE__,__FILE__)
       ENDIF

    ENDIF

    IF (glocal%tglocprint) THEN
       ! print selected localized functions
       il_psi = fpar%kr1s*fpar%kr2s*fpar%kr3s
       ALLOCATE(psi(il_psi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL wfnz_print(2,1,c1,nstate,psi,nstate)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (glocal%tglocrealp) THEN
       ! print selected localized functions
       il_psi = fpar%kr1*fpar%kr2s*fpar%kr3s
       ALLOCATE(psi(il_psi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       CALL p_real_space(2,1,c1,nstate,psi,nstate)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    IF (test) THEN
       DO is = 1,nstate
          filtr='ORB_fi_'
          IF (paral%io_parent)&
               WRITE(numb,'(I8)') is
          CALL xstring(numb,ia,ie)
          filorb=filtr(1:7)//numb(ia:ie)//'.dat'
          IF (paral%io_parent)&
               CALL fileopen(17,filorb,fo_def,ferror)
          DO ig = nkpt%ngwk,ncpw%ngw,-1
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
             IF (paral%io_parent)&
                  WRITE(17,'(I8,2(2x,1PE12.5))') -ig+ncpw%ngw,mod,LOG10(mod)
          ENDDO
          DO ig = 1,ncpw%ngw
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
             IF (paral%io_parent)&
                  WRITE(17,'(I8,2(2x,1PE12.5))') ig,mod,LOG10(mod)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(17)
       ENDDO

       over_sup = 0.0_real_8
       DO is = 1,nstate
          DO ig = 1,nkpt%ngwk
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             DO js = is+1,nstate
                mod2= REAL(c1(ig,js))*REAL(c1(ig,js))+&
                     AIMAG(c1(ig,js))*AIMAG(c1(ig,js))
                over_sup = over_sup +  mod*mod2
             ENDDO
          ENDDO
       ENDDO
       CALL mp_sum(over_sup,parai%allgrp)
       IF (paral%io_parent)&
            WRITE(6,'(A,2(1PE12.4))') '   >>>>>>> FINAL SUPPORT OVERLAP ',&
            over_sup,over_sup*2.0_real_8/(nstate*(nstate-1))
    ENDIF
    ! stop
    ! CALL WRITE_IREC(IREC)
    ! IREC(10) = 0
    ! CALL ZHWWF(2,IREC,C1,C2,NSTATE,EIGV,TAU0,TAUP,TAUP,NFI)

    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE g_loc_dr_real
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE g_loc_dr_comp(c0,c2,nstate,eigv,tau0,taup,nfi)
    ! ==--------------------------------------------------------------==
    ! == PREPARE ENVIRONMENT FOR LOCALIZATION IN G SPACE              ==
    ! ==--------------------------------------------------------------==

    ! input
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c2(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate,nkpt%nkpts)
    REAL(real_8)                             :: eigv(*), tau0(:,:,:), &
                                                taup(:,:,:)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_loc_dr_comp'

    CHARACTER(len=15)                        :: filorb, filtr
    CHARACTER(len=8)                         :: numb
    COMPLEX(real_8), ALLOCATABLE             :: c1(:,:), psi(:)
    INTEGER                                  :: ia, icont, ie, ierr, ig, &
                                                ikind, il_psi, ip, ipw, &
                                                irec(100), is, istate, js, &
                                                len1, msgid, ngwip, nstate_loc
    INTEGER, ALLOCATABLE                     :: mapw(:)
    LOGICAL                                  :: ferror, tgrid
    REAL(real_8)                             :: fac, mod, mod2, over_sup

! ==--------------------------------------------------------------==

    len1=2*ncpw%ngw
    ALLOCATE(mapw(len1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    IF (lostate%state_all) THEN
       nstate_loc = nstate
    ELSEIF (lostate%state_range) THEN
       IF (indstate%ist_last .GT. nstate) CALL stopgm('G_LOC_DR',&
            'WRONG STATES  RANGE, IST_LAST > NSTATE',& 
            __LINE__,__FILE__)
       nstate_loc = indstate%ist_last - indstate%ist_first + 1
       ALLOCATE(ind_st(nstate_loc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       icont = 0
       DO istate =  indstate%ist_first,indstate%ist_last
          icont = icont + 1
          ind_st(icont) = istate
       ENDDO
    ELSE
       nstate_loc = indstate%nst_list
       IF (indstate%nst_list .GT. nstate) CALL stopgm('G_LOC_DR',&
            'TOO MANY STATES, NST_LIST > NSTATE',& 
            __LINE__,__FILE__)
       ALLOCATE(ind_st(nstate_loc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO istate =  1,indstate%nst_list
          IF (ist_list(istate) .GT. nstate) CALL stopgm('G_LOC_DR',&
               'WRONG STATES INDEX IN THE LIST',& 
               __LINE__,__FILE__)
          ind_st(istate) = ist_list(istate)
       ENDDO
    ENDIF

    ALLOCATE(c1(2*ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1)!,SIZE(c1))


    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64(1H*))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)')&
            'LOCALIZATION OF WAVEFUNCTION IN RECIPROCAL SPACE '
       IF (gloci%gloc_type .EQ. 1) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       FUNCTIONAL:  Sum_(n,i)[ <n|Gi^2|n> - <n|Gi|n>^2 ]'
       ELSEIF (gloci%gloc_type .EQ. 2) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       FUNCTIONAL:  Sum_(n,i)[1 - |<n| EXP(iGir_o) |n>|^2]'
       ENDIF
       IF ((gloci%gloc_opt.EQ.1).AND.paral%io_parent)&
            WRITE(6,'(A)')&
            '       OPTIMISATION: U UPDATE BY STEEPEST DESCENT'
       IF ((gloci%gloc_opt.EQ.2) .AND.paral%io_parent)&
            WRITE(6,'(A,/,A)')&
            '       OPTIMISATION: U NEAR TO 1, WFN UPDATE EACH STEP',&
            '       SMALL ROTATIONS ASSUMED'
       IF (gloci%gloc_const .EQ. 1) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,/,A,1PE10.3,A,I10)')&
               '       UNITARITY  BY LAGRANGIAN MULTIPLIERS',&
               '       TOLERANCE = ',glocr%gepslag,' MAX n ITERATION = ',&
               gloci%gloc_maxit
       ELSEIF (gloci%gloc_const .EQ. 2) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       UNITARITY  BY ORTHOGONALIZATION [S^(-1/2)]'
       ELSEIF (gloci%gloc_const .EQ. 3) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '       UNITARITY  APPROXIMATION U = 1+iA (A HERMITIAN)'
       ENDIF
       IF (glocal%tg_linesearch) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               '       TIME STEP TUNING BY SEARCH ALONG THE LINE'
       ENDIF
       IF (lostate%state_all) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,I8)')&
               '       ALL THE STATES ARE CONSIDERED, NSTATE = ',nstate
       ELSEIF (lostate%state_range) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,I8,A,I8,A,I8)')&
               '      ONLY STATES BETWEEN ',indstate%ist_first,' AND ',indstate%ist_last,&
               ' ARE CONSIDERED, NSTATE = ',indstate%ist_last - indstate%ist_first + 1
       ELSEIF (lostate%state_list) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/A)')&
               '      ONLY THE FOLLOWING STATES ARE CONSIDERED: '
          IF (paral%io_parent)&
               WRITE(6,*) (ist_list(istate),istate = 1,indstate%nst_list)
       ENDIF

       IF (glocal%tg_read_matrix) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               '       FIRST UNITARY TRANSFORMATION READ FROM FILE'
       ENDIF
       IF (glocal%tgwannier) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/A,A)')&
               '       THE SPREAD FUNCTIONAL IN REAL',&
               ' SPACE IS ALSO CALCULATED'
          IF (paral%io_parent)&
               WRITE(6,'(A,f6.3/A,f6.3)')&
               '            optimization: W_GFUN = ',glocr%g2g_weight,&
               '                        : W_RFUN = ',glocr%wan_weight
          IF (gloci%gloc_opt.NE.1 .AND. gloci%gloc_type .NE. 2)&
               CALL stopgm('G_LOC_DR',&
               'REAL SPACE LOC. ONLY WITH GLOC_OPT=1, GLOC_TYPE=2',& 
               __LINE__,__FILE__)
       ENDIF
       IF (glocal%tg_antisymm) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,A)')&
               'PENALTY FUNCTION FOR SYMMETRIC WFN CALCULATED',&
               '  FOR FIRST 100 STEPS'
          IF ((glocal%tg_penalty).AND.paral%io_parent)&
               WRITE(6,'(A,/)')&
               '     PENALTY IMPOSED IN THE GRADIENT (ADDED TERM)'
       ENDIF
       ! IF(TCGRAD) THEN
       ! IF (GLOC_OPT.EQ.2) CALL STOPGM('G_LOC_DR_REAL',
       ! &          'CONJUGATE GRADIENT and U ALMOST 1 NOT POSSIBLE')
       ! WRITE(6,'(A)') '       CONJUGATE GRADIENT USED'
       ! ENDIF

       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         CONVERGENCE CRITERIA:',glocr%gloc_eps
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            '         MAXIMUM # OF STEPS:',gloci%gloc_maxs
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         RANDOMIZATION AMPLITUDE:',glocr%gloc_ran
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.4)')&
            '         STEP SIZE:',glocr%gloc_step
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            '         INITIAL STEP',gloci%gloc_init+1
       IF (glocal%tglocprint) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1x,A)')&
               'STARTING  AND OPTIMIZED WFN PRINTED IN SEPARETED FILES'
          IF (gloci%gloc_all .EQ. 1) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>>  all wfn printed  <<<")')
          ELSEIF (gloci%gloc_first .GT. 0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(">>>> wfn from ",I5," to ",I6," printed  <<<<")')
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(">>>>  printed wfn list number = ",I6)') gloci%gloc_orb
          ENDIF
       ENDIF
    ENDIF

    wan_mem = glocr%wan_weight
    g2g_mem = glocr%g2g_weight


    ! NSTATE = NSTATE - 6

    DO ikind = 1,nkpt%nkpts

       IF (tkpts%tkpnt)THEN
          shiftkx = rk(1,ikind)
          shiftky = rk(1,ikind)
          shiftkz = rk(1,ikind)
       ELSE
          shiftkx = 0.0_real_8
          shiftky = 0.0_real_8
          shiftkz = 0.0_real_8
       ENDIF

       IF (lostate%state_all) THEN
          CALL dcopy(4*ncpw%ngw*nstate_loc,c0(1,1,ikind),1,c1,1)
       ELSEIF (lostate%state_range) THEN
          CALL dcopy(4*ncpw%ngw*nstate_loc,c0(1,indstate%ist_first,ikind),1,c1,1)
       ELSE
          DO istate = 1,nstate_loc
             CALL dcopy(4*ncpw%ngw,c0(1,ind_st(istate),ikind),1,&
                  c1(1,istate),1)
          ENDDO
       ENDIF

       gifirst = 1

       filtr='SPR_K'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') ikind
       CALL xstring(numb,ia,ie)
       filspr=filtr(1:5)//numb(ia:ie)//'.dat'
       filtr='CEN_K'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') ikind
       CALL xstring(numb,ia,ie)
       filcen=filtr(1:5)//numb(ia:ie)//'.xyz'
       filtr='PENALTY_K'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') ikind
       CALL xstring(numb,ia,ie)
       filpen=filtr(1:9)//numb(ia:ie)//'.dat'
       filtr='MATRIX_K'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') ikind
       CALL xstring(numb,ia,ie)
       filmat=filtr(1:8)//numb(ia:ie)

       IF (glocal%tglocprint) THEN
          ! print selected localized functions
          il_psi = fpar%kr1s*fpar%kr2s*fpar%kr3s
          ALLOCATE(psi(il_psi),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL wfnz_print(1,ikind,c0(1,1,ikind),nstate,psi,nstate_loc)
          DEALLOCATE(psi,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF

       IF (glocal%tglocrealp) THEN
          ! print selected localized functions
          il_psi = fpar%kr1*fpar%kr2s*fpar%kr3s
          ALLOCATE(psi(il_psi),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL p_real_space(1,ikind,c0(1,1,ikind),nstate,psi,&
               nstate_loc)
          DEALLOCATE(psi,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF

       glocal%torbdist = .FALSE.
       tgrid    = .TRUE.

       IF (ikind .EQ. 1 .AND. glocal%torbdist) THEN
          DO is = 1,nstate_loc
             filtr='ORB_in_'
             IF (paral%io_parent)&
                  WRITE(numb,'(I8)') is
             CALL xstring(numb,ia,ie)
             filorb=filtr(1:7)//numb(ia:ie)//'.dat'
             IF (paral%io_parent)&
                  CALL fileopen(17,filorb,fo_def,ferror)
             il_psi = fpar%kr1*fpar%kr2s*fpar%kr3s
             ALLOCATE(psi(il_psi),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DO ip=0,parai%nproc-1
                IF (paral%parent) THEN
                   IF (parap%pgroup(ip+1).EQ.parai%me) THEN
                      DO ig=1,ncpw%ngw
                         c2(ig,is)=c1(ig,is)
                         mapw(ig)=mapgp(ig)
                      ENDDO
                      DO ig=ncpw%ngw+1,2*ncpw%ngw
                         c2(ig,is)=c1(ig,is)
                      ENDDO
                   ELSE
                      msgid=1
                      !msglen = 2 * 2 * parap%sparm(3,ip) * 8
                      CALL mp_recv(c2(:,is),2 * parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
                      msgid=2
                      !msglen = parap%sparm(3,ip) * 8/irat
                      CALL mp_recv(mapw,parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
                   ENDIF
                   DO ipw=1,parap%sparm(3,ip)
                      psi(mapw(ipw))=c2(ipw,is)
                   ENDDO
                   ngwip=parap%sparm(3,ip)
                   DO ipw=1,ngwip
                      psi(mapw(ipw)+spar%ngws)=c2(ngwip+ipw,is)
                   ENDDO
                ELSE
                   IF (parap%pgroup(ip+1).EQ.parai%me) THEN
                      msgid=1
                      !msglen = 2 * 2 * ngw * 8
                      CALL mp_send(c1(:,is),2 * ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
                      msgid=2
                      !msglen = ngw * 8/irat
                      CALL mp_send(mapgp,ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
                   ENDIF
                ENDIF
             ENDDO

             IF (paral%parent) THEN
                DO ig = spar%ngwks,spar%ngws,-1
                   mod = REAL(psi(ig))*REAL(psi(ig))+&
                        AIMAG(psi(ig))*AIMAG(psi(ig))
                   IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
                   IF (paral%io_parent)&
                        WRITE(17,'(I8,2(2x,1PE12.5))') -ig+spar%ngws,mod,LOG10(mod)
                ENDDO
             ENDIF
             DEALLOCATE(psi,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             IF (paral%io_parent) CALL fileclose(17)
          ENDDO

       ENDIF

       over_sup = 0.0_real_8
       DO is = 1,nstate_loc
          DO ig = 1,nkpt%ngwk
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             DO js = is+1,nstate_loc
                mod2= REAL(c1(ig,js))*REAL(c1(ig,js))+&
                     AIMAG(c1(ig,js))*AIMAG(c1(ig,js))
                over_sup = over_sup +  mod*mod2
             ENDDO
          ENDDO
       ENDDO
       CALL mp_sum(over_sup,parai%allgrp)
       IF (paral%io_parent)&
            WRITE(6,'(A,2(1PE12.4))') '   >>>>>>> INITIAL SUPPORT OVERLAP ',&
            over_sup,over_sup*2.0_real_8/(nstate_loc*(nstate_loc-1))




       IF (gloci%gloc_type .EQ. 1)  THEN

          IF (gloci%gloc_opt .EQ. 1) THEN
             CALL g_loc_spread_sum&
                  (c1,c2,nstate_loc,nfi)
          ELSEIF (gloci%gloc_opt .EQ. 2) THEN
             CALL g_loc_spread_ide&
                  (c1,c2,nstate_loc,nfi)
          ELSE
             CALL stopgm('G_OC_DR_REAL',&
                  'OPTIMIZATION ALGORITHM NOT FOUND',& 
                  __LINE__,__FILE__)
          ENDIF

       ELSEIF (gloci%gloc_type .EQ.2) THEN

          IF (gloci%gloc_opt .EQ. 1) THEN
             CALL g_loc_exp_sum&
                  (c1,c2,nstate_loc,nfi,tau0)
          ELSEIF (gloci%gloc_opt .EQ. 2) THEN
             CALL g_loc_exp_ide&
                  (c1,c2,nstate_loc,nfi)
          ELSE
             CALL stopgm('G_OC_DR_REAL',&
                  'OPTIMIZATION ALGORITHM NOT FOUND',& 
                  __LINE__,__FILE__)
          ENDIF

       ENDIF

       IF (lostate%state_all) THEN
          CALL dcopy(4*ncpw%ngw*nstate_loc,c1,1,c0(1,1,ikind),1)
       ELSEIF (lostate%state_range) THEN
          CALL dcopy(4*ncpw%ngw*nstate_loc,c1,1,c0(1,indstate%ist_first,ikind),1)
       ELSE
          DO istate = 1,nstate_loc
             CALL dcopy(4*ncpw%ngw,c1(1,istate),1,&
                  c0(1,ind_st(istate),ikind),1)
          ENDDO
       ENDIF


       IF (glocal%tglocprint) THEN
          ! print selected localized functions
          il_psi = fpar%kr1s*fpar%kr2s*fpar%kr3s
          ALLOCATE(psi(il_psi),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL wfnz_print(2,ikind,c0(1,1,ikind),nstate,psi,nstate_loc)
          DEALLOCATE(psi,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF

       IF (glocal%tglocrealp) THEN
          ! print selected localized functions
          il_psi = fpar%kr1*fpar%kr2s*fpar%kr3s
          ALLOCATE(psi(il_psi),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL p_real_space(2,ikind,c0(1,1,ikind),nstate,psi,&
               nstate_loc)
          DEALLOCATE(psi,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF


       IF (ikind .EQ. 1) THEN
          DO is = 1,nstate_loc

             glocal%torbdist = .FALSE.
             IF ( glocal%torbdist) THEN
                filtr='ORB_fi_'
                IF (paral%io_parent)&
                     WRITE(numb,'(I8)') is
                CALL xstring(numb,ia,ie)
                filorb=filtr(1:7)//numb(ia:ie)//'.dat'
                IF (paral%io_parent)&
                     CALL fileopen(17,filorb,fo_def,ferror)
                il_psi = 2*spar%ngws 
                ALLOCATE(psi(il_psi),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                DO ip=0,parai%nproc-1
                   IF (paral%parent) THEN
                      IF (parap%pgroup(ip+1).EQ.parai%me) THEN
                         DO ig=1,ncpw%ngw
                            c2(ig,is)=c1(ig,is)
                            mapw(ig)=mapgp(ig)
                         ENDDO
                         DO ig=ncpw%ngw+1,nkpt%ngwk
                            c2(ig,is)=c1(ig,is)
                         ENDDO
                      ELSE
                         msgid=1
                         !msglen = 2 * 2 * parap%sparm(3,ip) * 8
                         CALL mp_recv(c2(:,is),2 * parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
                         msgid=2
                         !msglen = parap%sparm(3,ip) * 8/irat
                         CALL mp_recv(mapw,parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
                      ENDIF
                      DO ipw=1,parap%sparm(3,ip)
                         psi(mapw(ipw))=c2(ipw,is)
                      ENDDO
                      ngwip=parap%sparm(3,ip)
                      DO ipw=1,ngwip
                         psi(mapw(ipw)+spar%ngws)=c2(ngwip+ipw,is)
                      ENDDO
                   ELSE
                      IF (parap%pgroup(ip+1).EQ.parai%me) THEN
                         msgid=1
                         !msglen = 2 * 2 * ngw * 8
                         CALL mp_send(c1(:,is),2 * ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
                         msgid=2
                         !msglen = ngw * 8/irat
                         CALL mp_send(mapgp,ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
                      ENDIF
                   ENDIF
                ENDDO
                IF (paral%parent) THEN
                   DO ig = spar%ngwks,spar%ngws,-1
                      mod = REAL(psi(ig))*REAL(psi(ig))+&
                           AIMAG(psi(ig))*AIMAG(psi(ig))
                      IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
                      IF (paral%io_parent)&
                           WRITE(17,'(I8,2(2x,1PE12.5))') -ig+spar%ngws,mod,LOG10(mod)
                   ENDDO
                   DO ig = 1,spar%ngws
                      mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                           AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
                      IF (mod .LT. 1.e-12_real_8) mod = 1.e-12_real_8
                      IF (paral%io_parent)&
                           WRITE(17,'(I8,2(2x,1PE12.5))') ig,mod,LOG10(mod)
                   ENDDO
                ENDIF
                DEALLOCATE(psi,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
                IF (paral%io_parent)CALL fileclose(17)
             ENDIF         ! TORBDIST

             IF (tgrid) THEN
                filtr='GRID_fi_'
                IF (paral%io_parent)&
                     WRITE(numb,'(I8)') is
                CALL xstring(numb,ia,ie)
                filorb=filtr(1:8)//numb(ia:ie)//'.xyz'
                IF (paral%io_parent)&
                     CALL fileopen(17,filorb,fo_def,ferror)
                fac = parm%tpiba * fbohr
                DO ip = 0,parai%nproc-1
                   IF (parap%pgroup(ip+1).EQ.parai%me) THEN
                      DO ig = nkpt%ngwk,ncpw%ngw,-1
                         mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                              AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
                         IF (mod .GT. 0.5_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'N',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.1_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'O',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.01_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'F',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.005_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'C',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.001_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'P',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.0005_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'Fe',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                         ELSEIF (mod .GT. 0.0001_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)')&
                                 'Ni',-gk(1,ig-ncpw%ngw)*fac,&
                                 -gk(2,ig-ncpw%ngw)*fac,-gk(3,ig-ncpw%ngw)*fac,ig,mod
                            ! ELSEIF(MOD .GT. 0.00001_real_8) THEN
                            ! write(17,'(A,3f12.5,I10,1PE12.5)')
                            ! &                   'Co',-GK(1,IG-NGW)*FAC,
                            ! &                 -GK(2,IG-NGW)*FAC,-GK(3,IG-NGW)*FAC,IG,MOD
                            ! ELSEIF(MOD .GT. 0.000001_real_8) THEN
                            ! write(17,'(A,3f12.5,I10,1PE12.5)')
                            ! &                   'He',-GK(1,IG-NGW)*FAC,
                            ! &                 -GK(2,IG-NGW)*FAC,-GK(3,IG-NGW)*FAC,IG,MOD
                         ENDIF
                      ENDDO
                      DO ig = 1,ncpw%ngw
                         mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                              AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
                         IF (mod .GT. 0.5_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'N',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.1_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'O',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.01_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'F',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.005_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'C',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.001_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'P',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.0005_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'Fe',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                         ELSEIF (mod .GT. 0.0001_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(17,'(A,3f12.5,I10,1PE12.5)') 'Ni',gk(1,ig)*fac,&
                                 gk(2,ig)*fac,gk(3,ig)*fac,ig,mod
                            ! ELSEIF(MOD .GT. 0.00001_real_8) THEN
                            ! write(17,'(A,3f12.5,I10,1PE12.5)') 'Co',GK(1,IG)*FAC,
                            ! &                              GK(2,IG)*FAC,GK(3,IG)*FAC,IG,MOD
                            ! ELSEIF(MOD .GT. 0.000001_real_8) THEN
                            ! write(17,'(A,3f12.5,I10,1PE12.5)') 'He',GK(1,IG)*FAC,
                            ! &                              GK(2,IG)*FAC,GK(3,IG)*FAC,IG,MOD
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL mp_sync(parai%allgrp)
                ENDDO
                IF ((paral%parent).AND.paral%io_parent)&
                     CALL fileclose(17)
             ENDIF       ! TGRID
          ENDDO
       ENDIF

       over_sup = 0.0_real_8
       DO is = 1,nstate_loc
          DO ig = 1,nkpt%ngwk
             mod = REAL(c1(ig,is))*REAL(c1(ig,is))+&
                  AIMAG(c1(ig,is))*AIMAG(c1(ig,is))
             DO js = is+1,nstate_loc
                mod2= REAL(c1(ig,js))*REAL(c1(ig,js))+&
                     AIMAG(c1(ig,js))*AIMAG(c1(ig,js))
                over_sup = over_sup +  mod*mod2
             ENDDO
          ENDDO
       ENDDO
       CALL mp_sum(over_sup,parai%allgrp)
       IF (paral%io_parent)&
            WRITE(6,'(A,2(1PE12.4))') '   >>>>>>> FINAL SUPPORT OVERLAP ',&
            over_sup,over_sup*2.0_real_8/(nstate_loc*(nstate_loc-1))


    ENDDO

    CALL write_irec(irec)
    irec(10) = 0
    CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,taup,taup,nfi)

    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE g_loc_dr_comp

  ! ==================================================================
  ! ==================================================================




END MODULE g_loc_dr_utils
