MODULE epr_current_p_utils
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nmr_position_p_utils,            ONLY: set_origin
  USE nmr_util_p_utils,                ONLY: make_123
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE response_pmod,                   ONLY: response1
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_current
  PUBLIC :: give_scr_epr_current

CONTAINS

  ! ==================================================================
  SUBROUTINE epr_current(c0,c1,c1b,c1t,nstate,BIND,j_a,j_b,&
       psi,epr_chi)
    ! ==================================================================
    ! ==--------------------------------------------------------------==
    ! Subroutine adapted from nmr_current_p
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c1t(ncpw%ngw,nstate,1), c1b(ncpw%ngw,nstate,3), &
      c1(ncpw%ngw,nstate,3), c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: BIND(fpar%nnr1,9), &
                                                j_a(fpar%nnr1,9), &
                                                j_b(fpar%nnr1,9)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: epr_chi(3,3)

    CHARACTER(*), PARAMETER                  :: procedureN = 'epr_current'

    COMPLEX(real_8), ALLOCATABLE             :: Bindt(:,:)
    INTEGER                                  :: iat, ib, ibb, ibbb, ierr, &
                                                il_psi, il_rho, ip, ir, is, &
                                                isp, isub, nattot
    REAL(real_8)                             :: center(3)
    REAL(real_8), ALLOCATABLE                :: gpsi0(:,:), gpsi1L(:,:), &
                                                j(:,:), psi0(:), psi1L(:), &
                                                psi1p(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('   CURRENT',isub)
    CALL rhoe_psi_size(il_rho,il_psi)
    IF (il_rho .LT. fpar%kr1*fpar%kr2*fpar%kr3) THEN
       CALL stopgm('EPR_CURRENT','IL_RHO .LT. KR1*KR2*KR3',& 
            __LINE__,__FILE__)
    ENDIF
    IF (fpar%nnr1 .NE. fpar%kr1*fpar%kr2*fpar%kr3) THEN
       CALL stopgm('EPR_CURRENT','nnr1 .ne. kr1*kr2*kr3',& 
            __LINE__,__FILE__)
    ENDIF

    ! TODO check stat
    ALLOCATE(psi0(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gpsi0(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi1L(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psi1p(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gpsi1L(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(j(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(Bindt(ncpw%nhg,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! TODO check and remove this  
    IF (il_rho > fpar%nnr1) THEN
       ! need to understand what should be size of these arrays and why do
       ! we azzero using il_rho??? 
       IF (paral%io_parent)&
            WRITE(6,*) 'epr_current_p: REFACTORING ERROR','il_rho > nnr1'
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(      psi)!,     il_psi)
    CALL zeroing(     psi0)!,     il_rho)
    CALL zeroing(    gpsi0)!,   3*il_rho)
    CALL zeroing(    psi1L)!,     il_rho)
    CALL zeroing(   gpsi1L)!,   3*il_rho)
    CALL zeroing(    psi1p)!,   3*il_rho)
    CALL zeroing(        j)!,   3*il_rho)
    IF (paral%parent) CALL prmem('EPR (current)   ')
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
    DO ib=1,3
       CALL make_123(ib,ibb,ibbb)
       CALL zeroing(j)!,3*il_rho)

       DO ip=1,3
          epr_chi(ib,ip)=0._real_8
       ENDDO

       ! alfa current
       DO is=1,spin_mod%nsup
          CALL set_origin(is)
          CALL fft2tor(c0(1,is),psi0,c1b(1,is,ib),psi1L, psi,ncpw%ngw,.TRUE.)
          DO ip=1,3
             CALL apply_op_p(c0(1,is),gpsi0(1,ip),ip,ncpw%ngw)
             CALL apply_op_p(c1b(1,is,ib),gpsi1L(1,ip),ip,ncpw%ngw)
             CALL fft2tor(gpsi0(1,ip),gpsi0(1,ip),gpsi1L(1,ip),&
                  gpsi1L(1,ip),psi,ncpw%ngw,.TRUE.)
          ENDDO

          DO ir=1,fpar%nnr1
             DO ip=1,3
                j_a(ir,(ib-1)*3+ip) = j_a(ir,(ib-1)*3+ip)+ gpsi0(ir,&
                     ip)*psi1L(ir)- psi0(ir)*gpsi1L(ir,ip)
             ENDDO
          ENDDO

          CALL fft2tor(c1(1,is,ibb),psi1p(1,1),c1(1,is,ibbb),&
               psi1p(1,2),psi,ncpw%ngw,.TRUE.)
          CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,ibbb)
          CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,ibb)

          DO ir=1,fpar%nnr1
             DO ip=1,3
                j_a(ir,(ib-1)*3+ip) = j_a(ir,(ib-1)*3+ip) +&
                     gpsi0(ir,ip) * psi1p(ir,1)
             ENDDO
          ENDDO

          DO ip=1,3
             CALL apply_op_p(c1(1,is,ibb),psi1p(1,1),ip,ncpw%ngw)
             CALL apply_op_p(c1(1,is,ibbb),psi1p(1,2),ip,ncpw%ngw)
             CALL fft2tor(psi1p(1,1),psi1p(1,1),psi1p(1,2),psi1p(1,2),&
                  psi,ncpw%ngw,.TRUE.)
             CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,ibbb)
             CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,ibb)
             DO ir=1,fpar%nnr1
                j_a(ir,(ib-1)*3+ip) = j_a(ir,(ib-1)*3+ip) -&
                     psi0(ir)*psi1p(ir,1)
             ENDDO
          ENDDO
          ! Recover current density for last state
          DO ir=1,fpar%nnr1
             DO ip=1,3
                j(ir,ip)=j_a(ir,(ib-1)*3+ip)-j(ir,ip)
             ENDDO
          ENDDO
          CALL apply_op_rx(j(1,3),1._real_8,gpsi0(1,1),0._real_8,2)
          CALL apply_op_rx(j(1,2),-1._real_8,gpsi1L(1,1),0._real_8,3)

          CALL apply_op_rx(j(1,1),1._real_8,gpsi0(1,2),0._real_8,3)
          CALL apply_op_rx(j(1,3),-1._real_8,gpsi1L(1,2),0._real_8,1)

          CALL apply_op_rx(j(1,2),1._real_8,gpsi0(1,3),0._real_8,1)
          CALL apply_op_rx(j(1,1),-1._real_8,gpsi1L(1,3),0._real_8,2)

          DO ir=1,fpar%nnr1
             DO ip=1,3
                epr_chi(ib,ip)=epr_chi(ib,ip)+2._real_8*(gpsi0(ir,ip)+&
                     gpsi1L(ir,ip))
             ENDDO
          ENDDO
       ENDDO                 ! is (alfa states)

       DO is=spin_mod%nsup+1,nstate
          CALL set_origin(is)
          CALL fft2tor(c0(1,is), psi0, c1b(1,is,ib),&
               psi1L, psi,ncpw%ngw,.TRUE.)
          DO ip=1,3
             CALL apply_op_p(c0(1,is),gpsi0(1,ip),ip,ncpw%ngw)
             CALL apply_op_p(c1b(1,is,ib),gpsi1L(1,ip),ip,ncpw%ngw)
             CALL fft2tor(gpsi0(1,ip),gpsi0(1,ip),gpsi1L(1,ip),&
                  gpsi1L(1,ip),psi,ncpw%ngw,.TRUE.)
          ENDDO

          DO ir=1,fpar%nnr1
             DO ip=1,3
                j_b(ir,(ib-1)*3+ip) = j_b(ir,(ib-1)*3+ip)+&
                     gpsi0(ir,ip)*psi1L(ir)- psi0(ir)*gpsi1L(ir,ip)
             ENDDO
          ENDDO

          CALL fft2tor(c1(1,is,ibb),psi1p(1,1),c1(1,is,ibbb),&
               psi1p(1,2),psi,ncpw%ngw,.TRUE.)
          CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,ibbb)
          CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,ibb)

          DO ir=1,fpar%nnr1
             DO ip=1,3
                j_b(ir,(ib-1)*3+ip) = j_b(ir,(ib-1)*3+ip) +&
                     gpsi0(ir,ip)*psi1p(ir,1)
             ENDDO
          ENDDO

          DO ip=1,3
             CALL apply_op_p(c1(1,is,ibb),psi1p(1,1),ip,ncpw%ngw)
             CALL apply_op_p(c1(1,is,ibbb),psi1p(1,2),ip,ncpw%ngw)
             CALL fft2tor(psi1p(1,1),psi1p(1,1),psi1p(1,2),psi1p(1,2),&
                  psi,ncpw%ngw,.TRUE.)
             CALL apply_op_rx(psi1p(1,1),1._real_8,psi1p(1,1),0._real_8,ibbb)
             CALL apply_op_rx(psi1p(1,2),-1._real_8,psi1p(1,1),1._real_8,ibb)
             DO ir=1,fpar%nnr1
                j_b(ir,(ib-1)*3+ip) = j_b(ir,(ib-1)*3+ip)-&
                     psi0(ir)*psi1p(ir,1)
             ENDDO
          ENDDO
       ENDDO                 ! is (beta states)

       ! j = alfa-current + beta-current
       ! - (beta-current - alfa-current) = 2*alfa-current

       DO ir=1,fpar%nnr1
          DO ip=1,3
             j(ir,ip) = 2.0_real_8*j_a(ir,(ib-1)*3+ip)
          ENDDO
       ENDDO

       ! Current density in G space
       CALL fft2tog(j(1,1),j(1,1),j(1,2),j(1,2),psi,ncpw%nhg,.TRUE.)
       CALL ffttog(j(1,3),j(1,3),psi,ncpw%nhg,.TRUE.)

       ! INDUCED FIELD:
       ! B_ind = Oz:
       CALL apply_op_p(j(1,1),gpsi0(1,1),2,ncpw%nhg)
       CALL apply_op_p(j(1,2),Bindt(1,3),1,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,Bindt(1,3),1)
       ! B_ind = Oy:
       CALL apply_op_p(j(1,3),gpsi0(1,1),1,ncpw%nhg)
       CALL apply_op_p(j(1,1),Bindt(1,2),3,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,Bindt(1,2),1)
       ! B_ind = Ox:
       CALL apply_op_p(j(1,2),gpsi0(1,1),3,ncpw%nhg)
       CALL apply_op_p(j(1,3),Bindt(1,1),2,ncpw%nhg)
       CALL daxpy(2*ncpw%nhg,-1._real_8,gpsi0(1,1),1,Bindt(1,1),1)

       CALL nmr_interaction3(Bindt(1,1),bindt(1,2),bindt(1,3))
       ! Remove G=0 value
       IF (geq0) THEN
          Bindt(1,1)=0
          Bindt(1,2)=0
          Bindt(1,3)=0
       ENDIF

       CALL fft2tor(Bindt(1,1),BIND(1,(ib-1)*3+1),bindt(1,2),&
            BIND(1,(ib-1)*3+2),psi,ncpw%nhg,.TRUE.)
       CALL ffttor(Bindt(1,3),BIND(1,(ib-1)*3+3),psi,ncpw%nhg,.TRUE.)
       ! * FPI/(fine structure constant)^2
       CALL dscal(fpar%nnr1,0.0006691760526_real_8,BIND(1,(ib-1)*3+1),1)
       CALL dscal(fpar%nnr1,0.0006691760526_real_8,BIND(1,(ib-1)*3+2),1)
       CALL dscal(fpar%nnr1,0.0006691760526_real_8,BIND(1,(ib-1)*3+3),1)

       ! CALL cubefile('Bindx_B'//label//'.cube ',
       ! &                 Bind(1,(iB-1)*3+1),center,psi,.FALSE.)
       ! CALL cubefile('Bindy_B'//label//'.cube ',
       ! &                 Bind(1,(iB-1)*3+2),center,psi,.FALSE.)
       ! CALL cubefile('Bindz_B'//label//'.cube ',
       ! &                 Bind(1,(iB-1)*3+3),center,psi,.FALSE.)

    ENDDO                     ! IB

    ! ==--------------------------------------------------------------==
    DEALLOCATE(psi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gpsi0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi1L,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psi1p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gpsi1L,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(j,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(Bindt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('   CURRENT',isub)
    RETURN
  END SUBROUTINE epr_current
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE give_scr_epr_current(lscr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lscr
    CHARACTER(len=30)                        :: tag

    lscr = 1
    IF (.NOT. response1%tepr) RETURN
    lscr = 12*fpar%nnr1
    tag  = 'epr-response          '
    RETURN
  END SUBROUTINE give_scr_epr_current
  ! ==================================================================



END MODULE epr_current_p_utils
