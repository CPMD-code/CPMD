MODULE nmr_current_p_utils
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: gk,&
                                             hg
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nmr_position_p_utils,            ONLY: set_origin
  USE nmr_util_p_utils,                ONLY: make_123
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE response_pmod,                   ONLY: chi,&
                                             chi_si_to_shift_ppm,&
                                             nmr_options,&
                                             response1,&
                                             shift_factor
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_current
  PUBLIC :: give_scr_nmr_current

CONTAINS

  ! ==================================================================
  SUBROUTINE calc_current(c0,c1,nstate,psi)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,9), &
                                                c0(ncpw%ngw,nstate), psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_current'
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: labels = (/'x','y','z'/)

    CHARACTER(len=1)                         :: label
    COMPLEX(real_8), ALLOCATABLE             :: BIND(:,:)
    INTEGER                                  :: i_offs_for_psiL, iat, iB, &
                                                iBB, iBBB, ierr, il_psi, &
                                                il_rho, ip, ir, is, isa, isp, &
                                                isub, nattot
    REAL(real_8)                             :: center(3), chi_scr(3,3), &
                                                chifactor, j2
    REAL(real_8), ALLOCATABLE                :: gpsi0(:,:), gpsi1L(:,:), &
                                                j(:,:), psi0(:), psi1L(:), &
                                                psi1p(:,:), this_shift(:,:,:)
    REAL(real_8), EXTERNAL                   :: eigr_dot2

! ==--------------------------------------------------------------==

    CALL tiset('   CURRENT',isub)
    CALL rhoe_psi_size(il_rho,il_psi)
    IF (il_rho .LT. fpar%kr1*fpar%kr2*fpar%kr3)&
         CALL stopgm('NMR_CURRENT','il_rho .lt. kr1*kr2*kr3',& 
         __LINE__,__FILE__)
    IF (fpar%nnr1 .NE. fpar%kr1*fpar%kr2*fpar%kr3)&
         CALL stopgm('NMR_CURRENT','nnr1 .ne. kr1*kr2*kr3',& 
         __LINE__,__FILE__)

    ! BUGFIX:
    i_offs_for_psiL = 3
    IF (nmr_options%tcurrent .AND. nmr_options%tnmr_full) i_offs_for_psiL = 6
    ! BUGFIX ends.



    ALLOCATE(BIND(ncpw%nhg, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(j(il_rho, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(psi0(il_rho),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi1p(il_rho, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi1L(il_rho),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(gpsi0(il_rho, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gpsi1L(il_rho, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(this_shift(3,3,maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   this_shift)!,3*3*maxsys%nax*maxsys%nsx)
    ! ==--------------------------------------------------------------==
    CALL zeroing(     psi)!,     il_psi)
    CALL zeroing(     psi0)!,     il_rho)
    CALL zeroing(    gpsi0)!,   3*il_rho)
    CALL zeroing(    psi1L)!,     il_rho)
    CALL zeroing(   gpsi1L)!,   3*il_rho)
    CALL zeroing(    psi1p)!,   3*il_rho)
    CALL zeroing(        j)!,   3*il_rho)
    IF (paral%parent) CALL prmem('NMR (current)   ')
    ! ==--------------------------------------------------------------==
    ! Find the center of the molecule, used for the .xyz file
    nattot=0
    DO isp=1,ions1%nsp
       nattot=nattot+ions0%na(isp)
    ENDDO
    DO ir=1,3
       center(ir) = 0._real_8
       DO isp=1,ions1%nsp
          DO iat=1,ions0%na(isp)
             center(ir) = center(ir) + tau0(ir,iat,isp)
          ENDDO
       ENDDO
       center(ir) = center(ir)/REAL(nattot,kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    DO iB=1,3
       label = labels(iB)
       CALL make_123(iB,iBB,iBBB)
       CALL zeroing(j)!,3*il_rho)
       DO is=1,nstate
          CALL set_origin(is)
          CALL fft2tor(c0(1,is), psi0, c1(1,is,i_offs_for_psiL+iB),&
               psi1L,psi,ncpw%ngw,.TRUE.)
          DO ip=1,3
             CALL apply_op_p(c0(1,is),      gpsi0 (1,ip),ip, ncpw%ngw)
             CALL apply_op_p(c1(1,is,iB+i_offs_for_psiL),&
                  gpsi1L(1,ip),ip, ncpw%ngw)
             CALL fft2tor(gpsi0(1,ip),gpsi0(1,ip),&
                  gpsi1L(1,ip),gpsi1l(1,ip), psi, ncpw%ngw,.TRUE.)
          ENDDO



          ! Magnetic field B in direction a
          ! current components j_b
          ! => j_ba = current component b for a field in direction a
          ! TERM 1:    t1  =  (p_b psi0) (r) * psi1^(L_a) (r)
          ! TERM 2:    t2  =  psi0(r) *  (p_b psi1^(L_a)) (r)

          DO ir=1,fpar%nnr1
             DO ip=1,3
                j(ir,ip) = j(ir,ip)&
                                ! !   &                 + gpsi0(ir,ip) * psi1L(ir) ! t1
                                ! !   &                 - psi0(ir) * gpsi1L(ir,ip) ! t2
                     + crge%f(is,1)*( gpsi0(ir,ip) * psi1L(ir)   & ! t1
                     -psi0(ir) * gpsi1L(ir,ip) ) ! t2
             ENDDO
          ENDDO
          ! these two terms result in 
          ! identical shift values, by the way.



          ! TERM 3:   t3 = (p_b psi0) (ir)  *  ( (r-d)_mu  phi1^(p_nu) )(ir) 
          ! \                               * eps_a_nu_mu
          ! TERM 4:   t4 = phi0(ir) * ( p_b (r-d)_mu  psi1^(p_nu) ) (ir)
          ! \                       * eps_a_nu_mu
          ! \            = phi0(ir) * ( (r-d)_mu  p_b psi1^(p_nu) ) (ir)
          ! \                       * eps_a_nu_my
          ! BECAUSE the (r-d)_mu term is NOT an operator, it comes from
          ! the CSGT-"r-prime". So it can only be an operator when next
          ! to the   |r> <r|  of the current density operator, i.e. when
          ! already having applied the p_b to psi1^(p_nu).


          CALL fft2tor(c1(1,is,iBB), psi1p(1,1),&
               c1(1,is,iBBB), psi1p(1,2),psi,ncpw%ngw,.TRUE.)
          CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,iBBB)
          CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,iBB)
          ! t3. In fact, we calculate (for iB=3, i.e. B || Oz)
          ! psi1p(*,1) =  (r-d)_y psi1^px  -  (r-d)_x psi1^py


          DO ir=1,fpar%nnr1
             DO ip=1,3
                j(ir,ip) = j(ir,ip)&
                                ! !   &                 + gpsi0(ir,ip) * psi1p(ir,1) ! t3
                     + crge%f(is,1)*gpsi0(ir,ip) * psi1p(ir,1) ! t3
             ENDDO
          ENDDO


          ! t4.
          ! \            = phi0(ir) * ( (r-d)_mu  p_b psi1^(p_nu) ) (ir)
          ! \                       * eps_a_nu_my


          DO ip=1,3
             CALL apply_op_p(c1(1,is,iBB),psi1p(1,1),ip,ncpw%ngw)
             CALL apply_op_p(c1(1,is,iBBB),psi1p(1,2),ip,ncpw%ngw)
             CALL fft2tor(psi1p(1,1),psi1p(1,1),&
                  psi1p(1,2),psi1p(1,2),psi,ncpw%ngw,.TRUE.)
             CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,iBBB)
             CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,iBB)
             DO ir=1,fpar%nnr1
                j(ir,ip) = j(ir,ip)&
                                ! !   &                 - psi0(ir) * psi1p(ir,1) ! t4
                     - crge%f(is,1)*psi0(ir) * psi1p(ir,1) ! t4
             ENDDO
          ENDDO

       ENDDO               ! is (states)



       CALL cubefile('jx_B'//label//'.cube ',&
            j(1,1),center,psi,.FALSE.)
       CALL cubefile('jy_B'//label//'.cube ',&
            j(1,2),center,psi,.FALSE.)
       CALL cubefile('jz_B'//label//'.cube ',&
            j(1,3),center,psi,.FALSE.)

       DO ir=1,fpar%nnr1
          j2 = j(ir,1)*j(ir,1) + j(ir,2)*j(ir,2)&
               + j(ir,3)*j(ir,3)
          psi1L(ir) = SQRT(j2)! abuse of psi1L
       ENDDO
       CALL cubefile('jabs_B'//label//'.cube     ',&
            psi1L,center,psi,.FALSE.)


       ! Current density in G space
       CALL fft2tog(j(1,1),j(1,1),&
            j(1,2),j(1,2),psi,ncpw%nhg,.TRUE.)
       CALL ffttog(j(1,3),j(1,3),psi,ncpw%nhg,.TRUE.)


#if 0
       ! DIVERGENCE of current density
       IF (iB .EQ.3) THEN
          CALL apply_op_p(j(1,1),gpsi0(1,1),1,nhg)
          CALL apply_op_p(j(1,2),gpsi0(1,2),2,nhg)
          CALL apply_op_p(j(1,3),gpsi0(1,3),3,nhg)

          CALL fft2tor(gpsi0(1,1),gpsi0(1,1),&
               gpsi0(1,2),gpsi0(1,2),psi,nhg,.TRUE.)
          CALL ffttor(gpsi0(1,3),gpsi0(1,3),psi,nhg,.TRUE.)
987       FORMAT ("div j(ir=",i7,")  =  ",f10.6," + ",&
               f10.6," + ",f10.6," = ",f10.6)
          DO ir=243550,243600
             IF (paral%io_parent)&
                  WRITE (6,987)&
                  ir,gpsi0(ir,1),gpsi0(ir,2),gpsi0(ir,3),&
                  gpsi0(ir,1)+gpsi0(ir,2)+gpsi0(ir,3)
          ENDDO
       ENDIF
#endif




       ! INDUCED FIELD:
       ! B_ind = Oz:
       CALL apply_op_p(j(1,1),gpsi0(1,1),2,ncpw%nhg)
       CALL apply_op_p(j(1,2),BIND(1,3),1,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,BIND(1,3),1)
       ! B_ind = Oy:
       CALL apply_op_p(j(1,3),gpsi0(1,1),1,ncpw%nhg)
       CALL apply_op_p(j(1,1),BIND(1,2),3,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,BIND(1,2),1)
       ! B_ind = Ox:
       CALL apply_op_p(j(1,2),gpsi0(1,1),3,ncpw%nhg)
       CALL apply_op_p(j(1,3),BIND(1,1),2,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,BIND(1,1),1)

       CALL nmr_interaction3(BIND(1,1),BIND(1,2),BIND(1,3))


       ! Add the susceptibility contribution to B(G=0)
       CALL dcopy(9,chi,1,chi_scr,1)
       CALL mp_sum(chi_scr,3*3,parai%allgrp)
       IF (geq0) THEN
          chifactor = - chi_SI_to_shift_ppm/shift_factor
          DO ir=1,3
             BIND(1,ir) = chifactor * chi_scr(ir,iB)
          ENDDO
       ENDIF

       CALL fft2tor(BIND(1,1),j(1,1),BIND(1,2),j(1,2),&
            psi,ncpw%nhg,.TRUE.)   ! abuse of j
       CALL ffttor(BIND(1,3),j(1,3),psi,ncpw%nhg,.TRUE.)
       CALL dscal(fpar%nnr1,shift_factor,j(1,1),1)
       CALL dscal(fpar%nnr1,shift_factor,j(1,2),1)
       CALL dscal(fpar%nnr1,shift_factor,j(1,3),1)
       CALL cubefile('Bindx_B'//label//'.cube ',&
            j(1,1),center,psi,.FALSE.)
       CALL cubefile('Bindy_B'//label//'.cube ',&
            j(1,2),center,psi,.FALSE.)
       CALL cubefile('Bindz_B'//label//'.cube ',&
            j(1,3),center,psi,.FALSE.)


       ! SHIFTS:
       isa = 0
       DO isp = 1,maxsys%nsx
          DO iat = 1, ions0%na(isp)
             isa = isa+1
             this_shift(1,iB,isa) = eigr_dot2(isa,BIND(1,1))
             this_shift(2,iB,isa) = eigr_dot2(isa,BIND(1,2))
             this_shift(3,iB,isa) = eigr_dot2(isa,BIND(1,3))
          ENDDO
       ENDDO
    ENDDO                     ! iB

    CALL mp_sum(this_shift,3*3*isa,parai%allgrp)

    IF (paral%parent) THEN
       isa = 0
       DO isp = 1,maxsys%nsx
          DO iat = 1, ions0%na(isp)
             isa = isa+1
             DO iB=1,3
                IF (paral%io_parent)&
                     WRITE (6,'("SHIELDINGS B_ext=",I1,'//&
                     '", atom ",I3," --> ",3F15.5)')&
                     iB,isa,&
                     (shift_factor * this_shift(ir,iB,isa),ir=1,3)
             ENDDO
          ENDDO
       ENDDO
    ENDIF


    ! ==--------------------------------------------------------------==
    DEALLOCATE(BIND,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(j,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(psi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi1p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi1L,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(gpsi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gpsi1L,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('   CURRENT',isub)
    RETURN
  END SUBROUTINE calc_current
  ! ==================================================================





  ! ==================================================================
  SUBROUTINE remove_divergence(j1,j2,iB)
    ! expects INPUT: j1,j2 in G space (packed format)
    ! where j1 = j_(iB+1), j2 = j_(iB+2).
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: j1(ncpw%nhg), j2(ncpw%nhg)
    INTEGER                                  :: iB

    COMPLEX(real_8)                          :: coeff
    INTEGER                                  :: iBB, iBBB, ig
    REAL(real_8)                             :: dotp, normafter(2), &
                                                normbefore(2)

    CALL stopgm('REMOVE_DIV','NO MORE SUPPORTED',& 
         __LINE__,__FILE__)

    IF (.NOT. nmr_options%tfast) THEN
       normbefore(1) = dotp(ncpw%nhg,j1,j1)
       normbefore(2) = dotp(ncpw%nhg,j2,j2)
       CALL mp_sum(normbefore,2,parai%allgrp)
       normbefore(1) = SQRT(normbefore(1))
       normbefore(2) = SQRT(normbefore(2))
    ENDIF
    CALL make_123(iB,iBB,iBBB)


    DO ig=1,ncpw%nhg
       coeff = j1(ig) * gk(iBB,ig)&
            + j2(ig) * gk(iBBB,ig)
       IF (hg(ig) .GT. 1.e-5_real_8) coeff = coeff / hg(ig)
       j1(ig) = j1(ig) - coeff * gk(iBB,ig)
       j2(ig) = j2(ig) - coeff * gk(iBBB,ig)
    ENDDO

    IF (.NOT. nmr_options%tfast) THEN
       normafter(1) = dotp(ncpw%nhg,j1,j1)
       normafter(2) = dotp(ncpw%nhg,j2,j2)
       CALL mp_sum(normafter,2,parai%allgrp)
       normafter(1) = SQRT(normafter(1))
       normafter(2) = SQRT(normafter(2))
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'REMOVE DIV: Norm 1: ',&
               normbefore(1),' -> ',normafter(1)
          IF (paral%io_parent)&
               WRITE(6,*)'REMOVE DIV: Norm 2: ',&
               normbefore(2),' -> ',normafter(2)
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE remove_divergence
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE give_scr_nmr_current(lscr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lscr
    CHARACTER(len=30)                        :: tag

    lscr = 1
    IF (.NOT. response1%tnmr) RETURN
    IF (.NOT. nmr_options%tcurrent) RETURN
    lscr = 9*fpar%nnr1
    tag  = 'nmr-response          '
    RETURN
  END SUBROUTINE give_scr_nmr_current
  ! ==================================================================



END MODULE nmr_current_p_utils
