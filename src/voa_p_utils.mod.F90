MODULE voa_p_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: scmass
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             velp
  USE ddip,                            ONLY: lenbk
  USE dipo_utils,                      ONLY: dipo
  USE dipomod,                         ONLY: moment
  USE dotp_utils,                      ONLY: dotp
  USE eicalc_utils,                    ONLY: eicalc1
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE fnonloc_p_utils,                 ONLY: fnonloc_p
  USE fnonloc_utils,                   ONLY: fnonloc
  USE h0psi1_p_utils,                  ONLY: h0psi1_p
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: localize
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             ndfnl
  USE nmr_position_p_utils,            ONLY: calc_lower_left_new,&
                                             print_llc,&
                                             set_origin
  USE nmr_util_p_utils,                ONLY: epsi,&
                                             make_123
  USE parac,                           ONLY: parai,&
                                             paral
  USE prop,                            ONLY: prop2
  USE response_pmod,                   ONLY: &
       dfnl00, ener1, fnl00, lower_left, lower_left_value, nmr_options, nris, &
       voa_data, voa_options
  USE rhoofr_p_utils,                  ONLY: rhoofr_p
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: infi
  USE rotate_utils,                    ONLY: rotate
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE spin,                            ONLY: clsd,&
                                             spin_mod
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
  USE utils,                           ONLY: nxxfun
  USE wann,                            ONLY: wannr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: voa_p

  INTERFACE print_tensor
     MODULE PROCEDURE print_tensor_rank2
     MODULE PROCEDURE print_tensor_rank3
  END INTERFACE print_tensor

CONTAINS

  ! ==================================================================
  SUBROUTINE voa_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,nstate)
    ! ==--------------------------------------------------------------==
    ! == NVPT calculations of VOA/VCD main routine                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:), &
                                                drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'voa_p'

    CHARACTER(30)                            :: filen
    INTEGER                                  :: iatom, idir, isub

650 FORMAT (1X,63("*"),1X)
    ! ==--------------------------------------------------------------==
    CALL tiset('       voa',isub)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND.voa_options%tverbose) THEN
       WRITE(6,*) ''
       WRITE(6,*) ''
       WRITE(6,650)
       WRITE(6,*) '*                                                             *'
       WRITE(6,*) '*            Nuclear Velocity Perturbation Theory (NVPT)      *'
       WRITE(6,*) '*                     calculations of                         *'
       WRITE(6,*) '*             Vibrational Optical Activity (VOA)              *'
       WRITE(6,*) '*                                                             *'
       WRITE(6,650)
       WRITE(6,*) 'L. A. Nafie'
       WRITE(6,*) 'J. Chem. Phys. 79, 4950 (1983)'
       WRITE(6,*) 'http://dx.doi.org/10.1063/1.445588'
       WRITE(6,*) ''
       WRITE(6,*) 'Scherrer, A; Sebastiani, D; Vuilleumier, R'
       WRITE(6,*) 'J. Chem. Theory Comput. 9, 5305 (2013)'
       WRITE(6,*) 'http://dx.doi.org/10.1021/ct400700c'
       WRITE(6,*) ''
       WRITE(6,*) 'Scherrer, A; Agostini, F; Sebastiani, D; Gross, E K U; Vuilleumier, R'
       WRITE(6,*) 'J. Chem. Phys. 143, 074106 (2015)'
       WRITE(6,*) 'http://dx.doi.org/10.1063/1.4928578'
       WRITE(6,650)
       WRITE(6,*) 'Selected Options:'
       IF(voa_options%tat)       WRITE(6,*) ' - Atomic Axial/Polar Tensor (AAT/APT) in velocity form'
       IF(voa_options%tmd)       WRITE(6,*) ' - Dipole Moments along Molecular Dynamics'
       IF(voa_options%tpf)       WRITE(6,*) ' - Atomic Polar Tensor (APT) in position form'
       IF(voa_options%tcurrent)  WRITE(6,*) ' - Electronic Probability Density Current'
       IF(voa_options%tverbose)  WRITE(6,*) ' - Verbose Output'
       IF(voa_options%tdebug)    WRITE(6,*) ' - Debug Modus'
       IF(voa_options%tatomlist) WRITE(6,*) ' - Calculaton only for Subset in Atomlist'
       IF(voa_options%tamat)     WRITE(6,*) ' - A-Matrix'
       WRITE(6,650)
    ENDIF
    IF (paral%io_parent.AND.voa_options%tdebug) THEN
       WRITE(6,*) 'DEBUG OUTPUT:'
       WRITE(6,*) 'voa_options%tat', voa_options%tat
       WRITE(6,*) 'voa_options%tmd', voa_options%tmd
       WRITE(6,*) 'voa_options%tpf', voa_options%tpf
       WRITE(6,*) 'voa_options%tcurrent', voa_options%tcurrent
       WRITE(6,*) 'voa_options%tdensity', voa_options%tdensity
       WRITE(6,*) 'voa_options%tverbose', voa_options%tverbose
       WRITE(6,*) 'voa_data%fmt ', voa_data%fmt
       WRITE(6,*) 'voa_data%w_eps', voa_data%w_eps
       WRITE(6,*) 'voa_options%tatomlist', voa_options%tatomlist
       WRITE(6,*) 'voa_options%atomlistline', TRIM(voa_options%atomlistline)
       WRITE(6,*) 'voa_options%tamat', voa_options%tamat
       WRITE(6,650)
    ENDIF
    nmr_options%tfast=.NOT.voa_options%tdebug
    ! ==--------------------------------------------------------------==
    IF(.NOT.((parm%ibrav.EQ.0).OR.(parm%ibrav.EQ.1).OR.(parm%ibrav.EQ.8)))&
         CALL stopgm(procedureN,'Only IBRAV 0/1/8 were tested',__LINE__,__FILE__)
    IF(.NOT.(voa_options%tat.OR.voa_options%tmd)) &
         CALL stopgm(procedureN,'Either AT or MD has to be set',__LINE__,__FILE__)
    IF(voa_options%tat.AND.voa_options%tmd) &
         CALL stopgm(procedureN,'AT and MD are not compatible',__LINE__,__FILE__)
    IF(voa_options%tcurrent.AND.(.NOT.voa_options%tmd)) &
         CALL stopgm(procedureN,'Current only available for MD',__LINE__,__FILE__)
    IF(voa_options%tdensity.AND.(.NOT.voa_options%tcurrent)) &
         CALL stopgm(procedureN,'Density only available with current',__LINE__,__FILE__)
    IF(voa_options%tpf.AND.(.NOT.voa_options%tat)) &
         CALL stopgm(procedureN,'apt_pf only available for AT',__LINE__,__FILE__)
    IF(voa_options%tamat.AND.(.NOT.voa_options%tat)) &
         CALL stopgm(procedureN,'A-Matrix only available for AT',__LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL initialization(nstate)
    CALL localization(c0,nstate)
    CALL commutator(voa_data%rc0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    IF(voa_options%tat) THEN
       CALL rhoofr(c0,rhoe,psi(:,1),nstate)
       CALL dipole_gs(psi,rhoe,eirop)
       DO iatom=1,ions1%nat
          IF(voa_data%atomlist(iatom)) THEN
             DO idir=1,3
                CALL pertham(c0,c1,psi(:,1),rhoe,drhoe,eirop,eivps,z11,nstate,iatom,idir)
                IF(voa_options%tpf) CALL dipole_pf(c0,psi,drhoe,eirop,nstate,iatom,idir)
                CALL sternheimer(c0,c1,psi(:,1),rhoe,drhoe,eirop,eivps,z11,nstate)
                CALL rotate(1._real_8,c1,0._real_8,voa_data%rc1,voa_data%rotmat,nstate,&
                     2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
                CALL moments(voa_data%rc1,rhoe,psi,nstate,idir,iatom)
                IF(voa_options%tamat) &
                     CALL dcopy(2*ncpw%ngw*nstate,c1,1,voa_data%c1src(1,1,idir,iatom),1)
             ENDDO
          ENDIF
       ENDDO
       IF(paral%io_parent) CALL atomic_tensors()
       IF(voa_options%tamat) CALL amatrix(c0,psi,drhoe,z11,nstate)
    ELSE IF(voa_options%tmd) THEN
       IF(voa_options%tdensity) THEN
          CALL rhoofr(c0,rhoe,psi(:,1),nstate)
          WRITE(filen,'(A,I6.6,A)') 'DENSITY-',infi,'.cube'
          CALL cubefile(filen,rhoe,voa_data%box_center,psi,.FALSE.)
       ENDIF
       CALL pertham(c0,c1,psi(:,1),rhoe,drhoe,eirop,eivps,z11,nstate,1,1)
       CALL sternheimer(c0,c1,psi(:,1),rhoe,drhoe,eirop,eivps,z11,nstate)
       CALL rotate(1._real_8,c1,0._real_8,voa_data%rc1,voa_data%rotmat,nstate,&
            2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       CALL moments(voa_data%rc1,rhoe,psi,nstate,1,1)
       CALL print_moments(nstate)
       IF(voa_options%tcurrent) CALL current(voa_data%rc0,voa_data%rc1,psi,nstate)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL finalization()
    ! ==--------------------------------------------------------------==
    CALL tihalt('       voa',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE voa_p
  ! ==================================================================
  SUBROUTINE initialization(nstate)
    ! ==--------------------------------------------------------------==
    ! == allocation of memory and initialization                      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'initialization'

    CHARACTER(24)                            :: form
    CHARACTER(4)                             :: atom_index
    INTEGER                                  :: i1, ia, iat, ierr, is

650 FORMAT (1X,63("*"),1X)
651 FORMAT (1X,I4.4,1X,L1,1X,A2,3X,F10.5,3X,F3.1) ! atom bool symb mass charge


    WRITE(*,'(A,I1,A,A,A)') '(',3,'(',TRIM(voa_data%fmt),'))'

    WRITE(form,'(A,I1,A,A,A)')'(',3,'(',TRIM(voa_data%fmt),'))'

    ALLOCATE(voa_data%h1psi0(ncpw%ngw,nstate),&
         voa_data%vnlc0(ncpw%ngw,nstate),&
         voa_data%pic0(ncpw%ngw,nstate,3),&
         voa_data%ric0(ncpw%ngw,nstate,3),&
         voa_data%vnlric0(ncpw%ngw,nstate,3),&
         voa_data%rxpic0(ncpw%ngw,nstate,3),&
         voa_data%atomlist(ions1%nat),&
         voa_data%p_vf(3,nstate),&
         voa_data%m_vf(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(voa_data%h1psi0)
    CALL zeroing(voa_data%vnlc0)
    CALL zeroing(voa_data%pic0)
    CALL zeroing(voa_data%ric0)
    CALL zeroing(voa_data%vnlric0)
    CALL zeroing(voa_data%rxpic0)
    CALL zeroing(voa_data%p_vf)
    CALL zeroing(voa_data%m_vf)

    IF(voa_options%tatomlist) THEN
       DO iat=1,ions1%nat
          WRITE(atom_index,'(I4.4)') iat
          voa_data%atomlist(iat) = INDEX(voa_options%atomlistline,atom_index).NE.0
       ENDDO
    ELSE
       DO iat=1,ions1%nat
          voa_data%atomlist(iat) = .TRUE.
       ENDDO
    ENDIF

    IF(voa_options%tat) THEN
       ALLOCATE(voa_data%apt_pf(3,3,ions1%nat),&
            voa_data%apt_vf(3,3,ions1%nat),&
            voa_data%aat_cog(3,3,ions1%nat),&
            voa_data%aat_dog(3,3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
       CALL zeroing(voa_data%apt_pf)
       CALL zeroing(voa_data%apt_vf)
       CALL zeroing(voa_data%aat_cog)
       CALL zeroing(voa_data%aat_dog)
    ENDIF

    IF(voa_options%tamat) THEN
       ALLOCATE(voa_data%c1src(ncpw%ngw,nstate,3,ions1%nat),&
            voa_data%amat_nl(3,ions1%nat,3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
       CALL zeroing(voa_data%c1src)
       CALL zeroing(voa_data%amat_nl)
    ENDIF

    IF(paral%io_parent.AND.voa_options%tverbose) THEN
       WRITE(6,'(A)') ' ATOMIC MASSES AND VALENCE CHARGES'
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat = iat + 1
             WRITE(6,651) iat,voa_data%atomlist(iat),elem%el(ions0%iatyp(is)),&
                  rmass%pma(is)/scmass,ions0%zv(is)
          ENDDO
       ENDDO
       WRITE(6,650)
       WRITE(6,'(A)') ' NUCLEAR POSITIONS'
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat = iat + 1
             WRITE(6,form) (tau0(i1,ia,is),i1=1,3)
          ENDDO
       ENDDO
       WRITE(6,650)
    ENDIF
    RETURN
  END SUBROUTINE initialization
  ! ==================================================================
  SUBROUTINE localization(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! == localization to max. localized Wannier orbitals              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'localization'

    COMPLEX(real_8), ALLOCATABLE             :: cs(:), sc0(:)
    INTEGER                                  :: i1, ierr, nc
    REAL(real_8)                             :: w_eps_save

650 FORMAT (1X,63("*"),1X)

    ALLOCATE(voa_data%rc0(ncpw%ngw,nstate),&
         voa_data%rc1(ncpw%ngw,nstate),&
         voa_data%rotmat(nstate,nstate),&
         lower_left(3,nstate),&
         lower_left_value(3,nstate),&
         voa_data%r_wc(3,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL dcopy(2*ncpw%ngw*nstate,c0,1,voa_data%rc0,1)
    CALL zeroing(voa_data%rc1)
    CALL zeroing(voa_data%rotmat)
    CALL zeroing(lower_left)
    CALL zeroing(lower_left_value)
    CALL zeroing(voa_data%r_wc)

    nris(1) = spar%nr1s
    nris(2) = spar%nr2s
    nris(3) = spar%nr3s
    prop2%numorb=MAX(prop2%numorb,nstate)
    nc=2*nkpt%ngwk*prop2%numorb
    lenbk=nxxfun(prop2%numorb)
    nc=MAX(2*lenbk*parai%nproc,nc)
    ALLOCATE(sc0(nc),cs(nc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    w_eps_save  = wannr%w_eps
    wannr%w_eps = voa_data%w_eps
    CALL localize(tau0,voa_data%rc0,cs,sc0,nstate)
    wannr%w_eps = w_eps_save
    DEALLOCATE(sc0,cs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)

    CALL calc_lower_left_new(voa_data%rc0,nstate)
    CALL print_llc(nstate)
    IF(paral%io_parent) THEN
       CALL print_tensor(voa_data%r_wc,'R_WC',3,nstate)
       WRITE(6,650)
    ENDIF
    DO i1=1,3
       voa_data%box_center(i1)=0.5_real_8*(parm%a1(i1)+parm%a2(i1)+parm%a3(i1))
    ENDDO
    RETURN
  END SUBROUTINE localization
  ! ==================================================================
  SUBROUTINE commutator(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! == eval. of: p_i|c0>, r_i|c0>, Vnl|c0>, Vnlr_i|c0>, rxp|c0>     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'commutator'

    INTEGER                                  :: i1, i2, i3, ierr, is
    REAL(real_8), ALLOCATABLE                :: tmpscr(:,:)

    ALLOCATE(tmpscr(fpar%nnr1,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(tmpscr)

    CALL rnlsm(c0,nstate,1,1,.FALSE.)
    CALL fnonloc(voa_data%vnlc0,crge%f,nstate,1,clsd%nlsd,.FALSE.)
    DO i1=1,3
       DO is=1,nstate
          CALL apply_op_p(c0(1,is),voa_data%pic0(1,is,i1),i1,ncpw%ngw)
       ENDDO
    ENDDO
    DO i1=1,3
       CALL make_123(i1,i2,i3)
       DO is=1,nstate
          CALL set_origin(is)
          CALL fft2tor(voa_data%pic0(1,is,i2),tmpscr(1,1),&
               voa_data%pic0(1,is,i3),tmpscr(1,2),&
               psi,ncpw%ngw,.FALSE.)
          CALL apply_op_rx(tmpscr(1,1),-1._real_8,rhoe,0._real_8,i3)
          CALL apply_op_rx(tmpscr(1,2),+1._real_8,rhoe,1._real_8,i2)
          CALL ffttog(rhoe,voa_data%rxpic0(1,is,i1),psi,ncpw%ngw,.FALSE.)

          CALL ffttor(c0(1,is),rhoe,psi,ncpw%ngw,.FALSE.)
          CALL apply_op_rx(rhoe,1._real_8,tmpscr(1,1),0._real_8,i1)
          CALL ffttog(tmpscr(1,1),voa_data%ric0(1,is,i1),psi,ncpw%ngw,.FALSE.)
       ENDDO
       CALL rnlsm(voa_data%ric0(:,:,i1),nstate,1,1,.FALSE.)
       CALL fnonloc(voa_data%vnlric0(:,:,i1),crge%f,nstate,1,clsd%nlsd,.FALSE.)
    ENDDO

    DEALLOCATE(tmpscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__) 

    RETURN
  END SUBROUTINE commutator
  ! ==================================================================
  SUBROUTINE pertham(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,nstate,iatom,idir)
    ! ==--------------------------------------------------------------==
    ! == calculation of the pert. Ham. <d(phi_0)/dR_nu,dot(R_nu)>     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), &
                                                drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: iatom, idir

    CHARACTER(*), PARAMETER                  :: procedureN = 'pertham'

    COMPLEX(real_8), ALLOCATABLE             :: eirop1(:), peirop1(:), &
                                                pv1_loc(:,:), &
                                                pv1_nonloc(:,:), v1_loc(:,:), &
                                                v1_nonloc(:,:)
    INTEGER                                  :: i1, ia, iat, ierr, il_eirop1, &
                                                il_v1_loc, il_v1_nonloc, is, &
                                                is0, ldfnl, lfnl
    REAL(real_8)                             :: velnorm, velnormtmp

    IF (paral%io_parent.AND.voa_options%tverbose)&
         WRITE(6,'(A)')' CALCULATION OF PERTURBATION HAMILTONIAN'
    CALL zeroing(voa_data%h1psi0)
    voa_data%timagpert = .FALSE.

    ! Potentials and related:
    il_v1_loc=2*ncpw%nhg*clsd%nlsd
    il_v1_nonloc=2*nkpt%ngwk*nstate*nkpt%nkpnt
    il_eirop1=2*ncpw%nhg
    ALLOCATE(v1_loc(ncpw%nhg,clsd%nlsd),eirop1(il_eirop1),&
         v1_nonloc(nkpt%ngwk,il_v1_nonloc/nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(v1_loc)
    CALL zeroing(v1_nonloc)
    CALL zeroing(eirop1)

    ! FNL,DFNL
    ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    lfnl = imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt
    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),&
         dfnl00(imagp,ions1%nat,maxsys%nhxs,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL rnlsm(c0,nstate,1,1,.TRUE.)
    CALL dcopy(lfnl,fnl,1,fnl00,1)
    CALL dcopy(ldfnl,dfnl,1,dfnl00,1)

    IF(voa_options%tat) THEN
       IF(paral%io_parent) WRITE(6,'(A,I4.4,A,I1)') &
            ' DISPLACEMENT OF ATOM ',iatom,' IN DIRECTION ',idir
       iat=0
       is0=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF(iatom.EQ.iat) is0=is
          ENDDO
       ENDDO
       CALL eicalc1(idir,is0,iatom,v1_loc,eirop1)
       CALL fnonloc_p(c0,psi,v1_nonloc,crge%f,nstate,idir,is0,iatom,1)
       CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,v1_loc,&
            v1_nonloc,z11,nstate,eirop1)
       CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c1,1,voa_data%h1psi0,1)
    ELSE IF(voa_options%tmd) THEN
       ALLOCATE(pv1_loc(ncpw%nhg,clsd%nlsd),peirop1(il_eirop1),&
            pv1_nonloc(nkpt%ngwk,il_v1_nonloc/nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
       CALL zeroing(pv1_loc)
       CALL zeroing(pv1_nonloc)
       CALL zeroing(peirop1)

       iat=0
       velnorm=0._real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF(voa_data%atomlist(iat)) THEN
                velnormtmp=0._real_8
                DO i1=1,3
                   velnormtmp=velnormtmp+velp(i1,ia,is)**2
                ENDDO
                velnorm=velnorm+dsqrt(velnormtmp)
                DO i1=1,3 ! projection on velocity
                   CALL eicalc1(i1,is,iat,v1_loc,eirop1)
                   CALL fnonloc_p(c0,psi,v1_nonloc,crge%f,nstate,i1,is,iat,1)
                   CALL daxpy(il_eirop1,velp(i1,ia,is),eirop1,1,peirop1,1)
                   CALL daxpy(il_v1_loc,velp(i1,ia,is),v1_loc,1,pv1_loc,1)
                   CALL daxpy(il_v1_nonloc,velp(i1,ia,is),v1_nonloc,1,pv1_nonloc,1)
                ENDDO
             ENDIF
          ENDDO
       ENDDO

       IF(velnorm.LE.1.E-9_real_8) &
            CALL stopgm(procedureN,'Too small velocity!',__LINE__,__FILE__)
       velnormtmp=1._real_8/velnorm
       IF(paral%io_parent) &
            WRITE(6,'(A,1PE16.6E2)') ' SCALING WITH 1/|v| =',velnormtmp
       CALL dscal(il_eirop1,velnormtmp,peirop1,1)
       CALL dscal(il_v1_loc,velnormtmp,pv1_loc,1)
       CALL dscal(il_v1_nonloc,velnormtmp,pv1_nonloc,1)
       CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,&
            pv1_loc,pv1_nonloc,z11,nstate,peirop1)
       CALL daxpy(2*ncpw%ngw*nstate,velnorm,c1,1,voa_data%h1psi0,1)

       DEALLOCATE(pv1_loc,pv1_nonloc,peirop1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)
    ENDIF

    DEALLOCATE(v1_loc,v1_nonloc,eirop1,fnl00,dfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)

    RETURN
  END SUBROUTINE pertham
  ! ==================================================================
  SUBROUTINE sternheimer(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,nstate)
    ! ==--------------------------------------------------------------==
    ! == calculation of phi_1=(H-eps)^(-1)<d(phi_0)/dR_nu,dot(R_nu)>  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), &
                                                drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'sternheimer'

    COMPLEX(real_8)                          :: eirop1(1), v1_loc(1,1)
    INTEGER                                  :: is

    IF (paral%io_parent.AND.voa_options%tverbose) &
         WRITE(6,'(A)')' CALCULATION OF THE ADIABATAIC CORRECTION'

    DO is=1,nstate
       CALL dscal(2*ncpw%ngw,crge%f(is,1),voa_data%h1psi0(1,is),1)
    ENDDO
    ener1%eloc1   = 0.0_real_8
    ener1%eht1    = 0.0_real_8
    voa_data%timagpert = .TRUE.
    CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,v1_loc,&
         voa_data%h1psi0,z11,nstate,eirop1)
    voa_data%timagpert = .FALSE.

    RETURN
  END SUBROUTINE sternheimer
  ! ==================================================================
  SUBROUTINE moments(c1,rhoe,psi,nstate,idir,iatom)
    ! ==--------------------------------------------------------------==
    ! == calculation of dipole moment derivatives (velocity forms)    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate)
    INTEGER                                  :: idir, iatom

    CHARACTER(*), PARAMETER                  :: procedureN = 'moments'

    COMPLEX(real_8), ALLOCATABLE             :: ric1(:,:,:), vnlc1(:,:)
    INTEGER                                  :: i1, i2, i3, ia, ia0, iat, &
                                                ierr, is, is0
    REAL(real_8)                             :: val_loc, val_nlc, vec(3)
    REAL(real_8), ALLOCATABLE                :: tmpscr(:,:)

650 FORMAT (A,1X,I4.4,1X,I1,1X,3(E18.10)) ! label atom dir moment(3)
651 FORMAT (A,1X,I1,1X,I4.4,1X,3(E18.10)) ! label dir state moment(3)

    ALLOCATE(vnlc1(ncpw%ngw,nstate),ric1(ncpw%ngw,nstate,3),&
         tmpscr(fpar%nnr1,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(vnlc1)
    CALL zeroing(ric1)
    CALL zeroing(tmpscr)

    CALL rnlsm(c1,nstate,1,1,.FALSE.)
    CALL fnonloc(vnlc1,crge%f,nstate,1,clsd%nlsd,.FALSE.)

    CALL zeroing(voa_data%p_vf)
    CALL zeroing(voa_data%m_vf)

    DO i1=1,3
       DO is=1,nstate
          CALL set_origin(is)
          CALL ffttor(c1(1,is),rhoe,psi,ncpw%ngw,.FALSE.)
          CALL apply_op_rx(rhoe,1._real_8,tmpscr(1,1),0._real_8,i1)
          CALL ffttog(tmpscr(1,1),ric1(1,is,i1),psi,ncpw%ngw,.FALSE.)
          !val_loc=(<C0|p_i1)(|C1>)
          val_loc=crge%f(is,1)*dotp(ncpw%ngw,c1(:,is),voa_data%pic0(:,is,i1))
          !val_nlc=(<C1|Vnl)(r_i1|C0>)-(<C0|Vnl)(r_i1|C1>)
          val_nlc=dotp(ncpw%ngw,vnlc1(:,is),voa_data%ric0(:,is,i1))-&
               dotp(ncpw%ngw,voa_data%vnlc0(:,is),ric1(:,is,i1))
          IF(voa_options%tamat) voa_data%amat_nl(i1,iatom,idir)=&
               voa_data%amat_nl(i1,iatom,idir)+2._real_8*val_nlc
          voa_data%p_vf(i1,is)=2._real_8*(val_loc+val_nlc)
       ENDDO
    ENDDO
    CALL mp_sum(voa_data%p_vf,nstate*3,parai%allgrp)
    IF(voa_options%tamat)&
         CALL mp_sum(voa_data%amat_nl(:,iatom,idir),3,parai%allgrp)

    IF(voa_options%tdebug.AND.paral%io_parent) THEN
       WRITE(6,*) ''
       WRITE(6,*) ' LABEL DIRECTION STATE VALUE'
       DO i1=1,3
          DO is=1,nstate
             WRITE(6,651) 'P_VF ', i1, is, voa_data%p_vf(i1,is)
          ENDDO
       ENDDO
    ENDIF

    DO i1=1,3
       CALL make_123(i1,i2,i3)
       DO is=1,nstate
          !val_loc=<C0|r_i2 p_i3-r_i3 p_i2)(|C1>)
          val_loc=crge%f(is,1)*dotp(ncpw%ngw,c1(:,is),voa_data%rxpic0(:,is,i1))
          !val_nlc=(<C0|r_i3 Vnl)(r_i2|C1>)-(<C0|r_i2 Vnl)(r_i3|C1>)
          val_nlc=dotp(ncpw%ngw,voa_data%vnlric0(:,is,i3),ric1(:,is,i2))-&
               dotp(ncpw%ngw,voa_data%vnlric0(:,is,i2),ric1(:,is,i3))
          voa_data%m_vf(i1,is)=val_loc+val_nlc
       ENDDO
    ENDDO
    CALL mp_sum(voa_data%m_vf,nstate*3,parai%allgrp)

    IF(voa_options%tdebug.AND.paral%io_parent) THEN
       WRITE(6,*) ''
       WRITE(6,'(A)') 'LABEL DIRECTION STATE VALUE'
       DO i1=1,3
          DO is=1,nstate
             WRITE(6,651) 'M_VF ', i1, is, voa_data%m_vf(i1,is)
          ENDDO
       ENDDO
    ENDIF

    IF(voa_options%tat) THEN
       iat=0
       ia0=0
       is0=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF(iatom.EQ.iat) THEN
                ia0=ia
                is0=is
             ENDIF
          ENDDO
       ENDDO
       DO is=1,nstate
          DO i1=1,3
             voa_data%apt_vf(i1,idir,iatom)=voa_data%apt_vf(i1,idir,iatom)+voa_data%p_vf(i1,is)
             voa_data%aat_dog(i1,idir,iatom)=voa_data%aat_dog(i1,idir,iatom)+voa_data%m_vf(i1,is)
          ENDDO
       ENDDO
       CALL dcopy(3,voa_data%aat_dog(1,idir,iatom),1,voa_data%aat_cog(1,idir,iatom),1)
       DO is=1,nstate
          DO i1=1,3
             vec(i1) = voa_data%r_wc(i1,is)-tau0(i1,ia0,is0)
          ENDDO
          vec(1)=vec(1)-NINT(vec(1)/parm%a1(1))*parm%a1(1)
          vec(2)=vec(2)-NINT(vec(2)/parm%a2(2))*parm%a2(2)
          vec(3)=vec(3)-NINT(vec(3)/parm%a3(3))*parm%a3(3)
          DO i1=1,3
             DO i2=1,3
                DO i3=1,3
                   voa_data%aat_dog(i3,idir,iatom)=voa_data%aat_dog(i3,idir,iatom)+0.5_real_8*&
                        epsi(i1,i2,i3)*vec(i1)*voa_data%p_vf(i2,is)
                   voa_data%aat_cog(i3,idir,iatom)=voa_data%aat_cog(i3,idir,iatom)+0.5_real_8*&
                        epsi(i1,i2,i3)*(voa_data%r_wc(i1,is)-voa_data%box_center(i1))*voa_data%p_vf(i2,is)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(voa_options%tverbose.AND.paral%io_parent.AND.voa_options%tat) THEN
       WRITE(6,*) ''
       WRITE(6,'(A)') 'LABEL ATOM DIRECTION VALUES (1-3)'
       IF(voa_options%tpf) &
            WRITE(6,650) 'E_PF ', iatom, idir, (voa_data%apt_pf(i1,idir,iatom),i1=1,3)
       WRITE(6,650) 'E_VF ', iatom, idir, (voa_data%apt_vf(i1,idir,iatom),i1=1,3)
       WRITE(6,650) 'I_COG', iatom, idir, (voa_data%aat_cog(i1,idir,iatom),i1=1,3)
       WRITE(6,650) 'I_DOG', iatom, idir, (voa_data%aat_dog(i1,idir,iatom),i1=1,3)
       WRITE(6,*) ''
    ENDIF

    DEALLOCATE(vnlc1,ric1,tmpscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)

    RETURN
  END SUBROUTINE moments
  ! ==================================================================
  SUBROUTINE dipole_gs(psi,rhoe,eirop)
    ! ==--------------------------------------------------------------==
    ! == electric dipole moment (ground state)                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dipole_gs'

    INTEGER                                  :: i1, ir

    DO ir=1,fpar%nnr1
       psi(ir,1)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL fwfftn(psi(:,1),.FALSE.,parai%allgrp)
    CALL dipo(tau0,eirop,psi)
    DO i1=1,3
       voa_data%p_gs(i1)=moment%dmom(i1)
    ENDDO

    RETURN
  END SUBROUTINE dipole_gs
  ! ==================================================================
  SUBROUTINE dipole_pf(c0,psi,drhoe,eirop,nstate,iatom,idir)
    ! ==--------------------------------------------------------------==
    ! == APT for nuclear displacement (position form)                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)
    INTEGER                                  :: iatom, idir

    CHARACTER(*), PARAMETER                  :: procedureN = 'dipole_pf'

    INTEGER                                  :: i1, ir

    CALL rhoofr_p(c0,voa_data%h1psi0,drhoe,psi,nstate)
    DO ir=1,fpar%nnr1
       psi(ir,1)=CMPLX(drhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL fwfftn(psi(:,1),.FALSE.,parai%allgrp)
    CALL dipo(tau0,eirop,psi)
    DO i1=1,3
       voa_data%apt_pf(i1,idir,iatom)=moment%dmom(i1)
    ENDDO

    RETURN
  END SUBROUTINE dipole_pf
  ! ==================================================================
  SUBROUTINE atomic_tensors
    ! ==--------------------------------------------------------------==
    ! == nuclear contributios / printing / sumrules                   ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'atomic_tensors'

    INTEGER                                  :: i1, i2, i3, i4, ia, iat, is
    REAL(real_8)                             :: sum_rules(3,3,6)

650 FORMAT (1X,63("*"),1X)

    CALL zeroing(sum_rules)
    iat = 0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat = iat + 1
          IF(voa_data%atomlist(iat)) THEN
             DO i1=1,3
                ! N_12^nu = Z_nu delta_12
                voa_data%apt_pf(i1,i1,iat)=voa_data%apt_pf(i1,i1,iat)+ions0%zv(is)
                voa_data%apt_vf(i1,i1,iat)=voa_data%apt_vf(i1,i1,iat)+ions0%zv(is)
                DO i2=1,3
                   DO i3=1,3
                      ! I_12^nu = 0.5*Z_nu sum_3 eps_123 R^nu_3
                      voa_data%aat_cog(i2,i1,iat)=voa_data%aat_cog(i2,i1,iat)+0.5_real_8*&
                           ions0%zv(is)*epsi(i1,i2,i3)*(tau0(i3,ia,is)-voa_data%box_center(i3))
                   ENDDO
                ENDDO
                DO i2=1,3
                   ! SIGMA^0_12 = sum_nu P^nu_12
                   sum_rules(i2,i1,1)=sum_rules(i2,i1,1)+voa_data%apt_pf(i2,i1,iat)
                   sum_rules(i2,i1,2)=sum_rules(i2,i1,2)+voa_data%apt_vf(i2,i1,iat)
                   DO i3=1,3
                      DO i4=1,3
                         ! SIGMA^2_12 = sum_nu34 eps_234 P^nu_41 R^nu_3
                         sum_rules(i2,i1,4)=sum_rules(i2,i1,4)+epsi(i2,i3,i4)*voa_data%apt_pf(i1,i4,iat)*&
                              (tau0(i3,ia,is)-voa_data%box_center(i3))
                         sum_rules(i2,i1,5)=sum_rules(i2,i1,5)+epsi(i2,i3,i4)*voa_data%apt_vf(i1,i4,iat)*&
                              (tau0(i3,ia,is)-voa_data%box_center(i3))
                      ENDDO
                   ENDDO
                   ! SIGMA^3 = 2c sum_nu M^nu_12
                   sum_rules(i2,i1,6)=sum_rules(i2,i1,6)+2._real_8*voa_data%aat_cog(i2,i1,iat)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    DO i1=1,3
       DO i2=1,3
          DO i3=1,3
             ! SIGMA^1_12 = sum_c eps_123 dip_c
             sum_rules(i2,i1,3) = sum_rules(i2,i1,3)+epsi(i1,i2,i3)*voa_data%p_gs(i3)
          ENDDO
       ENDDO
    ENDDO

    WRITE(6,650)
    IF(voa_options%tpf) &
         CALL print_tensor(voa_data%apt_pf,'APT_PF',3,3,ions1%nat)
    CALL print_tensor(voa_data%apt_vf,'APT_VF',3,3,ions1%nat)
    CALL print_tensor(voa_data%aat_cog,'AAT_COG',3,3,ions1%nat)
    CALL print_tensor(voa_data%aat_dog,'AAT_DOG',3,3,ions1%nat)
    WRITE(6,650)
    IF(voa_options%tpf) &
         CALL print_tensor(sum_rules(:,:,1),'SIGMA_0_PF',3,3)
    CALL print_tensor(sum_rules(:,:,2),'SIGMA_0_VF',3,3)
    CALL print_tensor(sum_rules(:,:,3),'SIGMA_1_GS',3,3)
    IF(voa_options%tpf) &
         CALL print_tensor(sum_rules(:,:,4),'SIGMA_2_PF',3,3)
    CALL print_tensor(sum_rules(:,:,5),'SIGMA_2_VF',3,3)
    CALL print_tensor(sum_rules(:,:,6),'SIGMA_3_COG',3,3)
    WRITE(6,650)
    WRITE(6,'(A,F15.7)') 'VOA CHECKSUM=' , SUM(ABS(sum_rules)) 
    RETURN
  END SUBROUTINE atomic_tensors
  ! ==================================================================
  SUBROUTINE print_tensor_rank2(tensor,label,n1,n2)
    CHARACTER(*)                             :: label
    INTEGER                                  :: n1, n2
    REAL(real_8)                             :: tensor(n1,n2)

    CHARACTER(*), PARAMETER :: procedureN = 'print_tensor_rank2'

    INTEGER                                  :: i1, i2
    LOGICAL                                  :: ferror

    IF (.NOT.paral%io_parent) RETURN
    !WRITE(form,'(A,I1,A,A,A)')'(',n1,'(',TRIM(voa_data%fmt),'))'
    CALL fileopen(44,label,fo_app+fo_verb,ferror)
    WRITE(6,'(A)') label
    DO i2=1,n2
       WRITE(6,'(100(E18.10))') (tensor(i1,i2),i1=1,n1)
       WRITE(44,'(100(E18.10))') (tensor(i1,i2),i1=1,n1)
    ENDDO
    CALL fileclose(44)

    RETURN
  END SUBROUTINE print_tensor_rank2
  ! ==================================================================
  SUBROUTINE print_tensor_rank3(tensor,label,n1,n2,n3)
    CHARACTER(*)                             :: label
    INTEGER                                  :: n1, n2, n3
    REAL(real_8)                             :: tensor(n1,n2,n3)

    CHARACTER(*), PARAMETER :: procedureN = 'print_tensor_rank3'

    INTEGER                                  :: i1, i2, i3
    LOGICAL                                  :: ferror

    IF (.NOT.paral%io_parent) RETURN
    !WRITE(form,'(A,I1,A,A,A)')'(',n1,'(',TRIM(voa_data%fmt),'))'
    CALL fileopen(44,label,fo_app+fo_verb,ferror)
    WRITE(6,'(A)') label
    DO i3=1,n3
       DO i2=1,n2
          WRITE(6,'(100(E18.10))') (tensor(i1,i2,i3),i1=1,n1)
          WRITE(44,'(100(E18.10))') (tensor(i1,i2,i3),i1=1,n1)
       ENDDO
    ENDDO
    CALL fileclose(44)

    RETURN
  END SUBROUTINE print_tensor_rank3
  ! ==================================================================
  SUBROUTINE print_moments(nstate)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'print_moments'

    CHARACTER(12)                            :: filen
    INTEGER                                  :: i1, is
    LOGICAL                                  :: ferror

657 FORMAT (I7,9(2X,1PE16.6E2)) ! TSTEP R(3) P(3) M(3)

    IF (.NOT.paral%io_parent) RETURN
    WRITE(filen,'(A)') 'MOMENTS'
    WRITE(6,'(A)') filen
    CALL fileopen(44,filen,fo_app+fo_verb,ferror)
    DO is=1,nstate
       WRITE(6,657) infi,(voa_data%r_wc(i1,is),I1=1,3),&
            (voa_data%p_vf(i1,is),i1=1,3),(voa_data%m_vf(i1,is),i1=1,3)
       WRITE(44,657) infi,(voa_data%r_wc(i1,is),I1=1,3),&
            (voa_data%p_vf(i1,is),i1=1,3),(voa_data%m_vf(i1,is),i1=1,3)
    ENDDO
    CALL fileclose(44)

    RETURN
  END SUBROUTINE print_moments
  ! ==================================================================
  SUBROUTINE finalization
    ! ==--------------------------------------------------------------==
    ! == deallocation of memory                                       ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'finalization'

    INTEGER                                  :: ierr

    DEALLOCATE(voa_data%rotmat,voa_data%rc0,voa_data%rc1,voa_data%h1psi0,&
         voa_data%vnlc0,voa_data%pic0,voa_data%ric0,voa_data%vnlric0,&
         voa_data%rxpic0,lower_left,lower_left_value,voa_data%atomlist,&
         voa_data%r_wc,voa_data%p_vf,voa_data%m_vf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)
    IF(voa_options%tat) THEN
       DEALLOCATE(voa_data%apt_pf,voa_data%apt_vf,voa_data%aat_cog,&
            voa_data%aat_dog,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)
    ENDIF
    IF(voa_options%tamat) THEN
       DEALLOCATE(voa_data%c1src,voa_data%amat_nl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE finalization
  ! ==================================================================
  SUBROUTINE current(c0,c1,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! == calculation of electronic current                            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'current'

    CHARACTER(30)                            :: filen
    COMPLEX(real_8), ALLOCATABLE             :: pic1(:,:)
    INTEGER                                  :: i1, ierr, ir, is
    REAL(real_8)                             :: pref
    REAL(real_8), ALLOCATABLE                :: crnt(:), tmpscr(:,:)

    IF (cntl%tlsd) &
         CALL stopgm(procedureN,'voa: LSD not yet tested for current',__LINE__,__FILE__)

    ALLOCATE(tmpscr(fpar%nnr1,4),pic1(ncpw%ngw,nstate),crnt(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(tmpscr)
    CALL zeroing(pic1)

    pref=2._real_8/parm%omega
    DO i1=1,3
       CALL zeroing(crnt)
       DO is=1,nstate
          CALL apply_op_p(c1(1,is),pic1(1,is),i1,ncpw%ngw)
          CALL fft2tor(c0(1,is),tmpscr(1,1),c1(1,is),tmpscr(1,3),psi,ncpw%ngw,.FALSE.)
          CALL fft2tor(voa_data%pic0(1,is,i1),tmpscr(1,4),pic1(1,is),tmpscr(1,2),psi,ncpw%ngw,.FALSE.)
          DO ir=1,fpar%nnr1
             crnt(ir)=crnt(ir)+pref*(tmpscr(ir,1)*tmpscr(ir,2)-tmpscr(ir,3)*tmpscr(ir,4))
          ENDDO
       ENDDO
       WRITE(filen,'(A,I6.6,A,I1,A)') 'CURRENT-',infi,'-',i1,'.cube'
       CALL cubefile(filen,crnt,voa_data%box_center,psi,.FALSE.)
    ENDDO

    DEALLOCATE(tmpscr,pic1,crnt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__) 

    RETURN
  END SUBROUTINE current
  ! ==================================================================
  SUBROUTINE amatrix(c0,psi,drhoe,z11,nstate)
    ! ==--------------------------------------------------------------==
    ! == calculation of the vector potential matrix                   ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'amatrix'

    CHARACTER(30)                            :: filen
    COMPLEX(real_8), ALLOCATABLE             :: c12(:,:)
    INTEGER                                  :: i1, i2, iat1, iat2, ierr, is
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: corr
    REAL(real_8), ALLOCATABLE                :: amat(:,:,:,:)

675 FORMAT (1X,I3,1X,I1,1X,I3,1X,I1,1PE16.6E2,1PE16.6E2,1PE16.6E2)

    IF (cntl%tlsd) &
         CALL stopgm(procedureN,'voa: LSD not yet tested for amatrix',__LINE__,__FILE__)

    ALLOCATE(c12(nkpt%ngwk,nstate),amat(3,ions1%nat,3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',__LINE__,__FILE__)
    CALL zeroing(amat)

    DO iat1=1,ions1%nat
       IF(voa_data%atomlist(iat1)) THEN
          DO i1=1,3
             CALL zeroing(c12)
             ! c12 = - H_0 F | c1 (iat1,i1) >
             CALL h0psi1_p(voa_data%c1src(1,1,i1,iat1),c0,drhoe,crge%f,tau0,fion,psi(:,1),nstate,c12)
             ! c12 = - (H_0 - E_0) F | c1 (iat1,i1) >
             DO is=1,nstate
                CALL daxpy(2*ncpw%ngw,z11(is,is),voa_data%c1src(1,is,i1,iat1),1,c12(1,is),1)
             ENDDO
             DO iat2=1,ions1%nat
                IF(voa_data%atomlist(iat2)) THEN
                   DO i2=1,3
                      DO is=1,nstate
                         ! - 2 < c1 (iat2,i2) | (H_0 - E_0) F | c1 (iat1,i1) >
                         amat(i1,iat1,i2,iat2)=amat(i1,iat1,i2,iat2)-&
                              2._real_8*dotp(ncpw%ngw,voa_data%c1src(:,is,i2,iat2),c12(:,is))
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL mp_sum(amat,3*3*ions1%nat*ions1%nat,parai%allgrp)

    IF(paral%io_parent) THEN
       WRITE(filen,'(A)') 'A_MATRIX'
       CALL fileopen(44,filen,fo_app+fo_verb,ferror)
       WRITE(6,'(A)') 'A_MATRIX'
       DO iat1=1,ions1%nat
          DO i1=1,3
             DO iat2=1,ions1%nat
                DO i2=1,3
                   IF (iat1.EQ.iat2) THEN
                      corr=0.5_real_8*(voa_data%amat_nl(i2,iat1,i1)+voa_data%amat_nl(i1,iat1,i2))
                   ELSE
                      corr=0._real_8
                   ENDIF
                   WRITE(6,675) iat1,i1,iat2,i2,&
                        amat(i1,iat1,i2,iat2),corr,amat(i1,iat1,i2,iat2)-corr
                   WRITE(44,675) iat1,i1,iat2,i2,&
                        amat(i1,iat1,i2,iat2),corr,amat(i1,iat1,i2,iat2)-corr
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       WRITE(6,'(1X,63("*"),1X)')
       CALL fileclose(44)
    ENDIF

    DEALLOCATE(c12,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',__LINE__,__FILE__) 

    RETURN
  END SUBROUTINE amatrix
  ! ==================================================================
END MODULE voa_p_utils
