MODULE v1ofrho_p_utils
  USE cppt,                            ONLY: nr1h,&
                                             nr2h,&
                                             nr3pl,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE reshaper,                        ONLY: reshape_inplace
  USE response_pmod,                   ONLY: ener1,&
                                             response1
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1xc_p_utils,                    ONLY: give_scr_v1xc_p,&
                                             v1xc_p
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: v1ofrho_p
  PUBLIC :: give_scr_v1ofrho

CONTAINS

  ! ==================================================================
  SUBROUTINE v1ofrho_p(v1_loc,eirop1,drhoe,v,psi,first)
    ! ==--------------------------------------------------------------==
    ! calculate the implicit part of the first order potential, i.e.
    ! the part induced by the response of the density wrt the bare
    ! perturbation.
    ! input: first order density in
    ! drhoe(r) = <phi1|r><r|phi0> + cc
    ! \       local potential in v_loc
    ! output: first order potential in v(r):
    ! \       v(r0)  =  drhoe(r) / |r0-r|  +  d v_xc / d rho  * drhoe
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: v1_loc(ncpw%nhg), &
                                                eirop1(ncpw%nhg)
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd)
    REAL(real_8), TARGET                     :: v(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:)
    LOGICAL                                  :: first

    CHARACTER(*), PARAMETER                  :: procedureN = 'v1ofrho_p'

    COMPLEX(real_8)                          :: eirop1gg
    COMPLEX(real_8), ALLOCATABLE             :: aux(:)
    COMPLEX(real_8), POINTER                 :: vg(:)
    INTEGER                                  :: ierr, ig, ig1, isub

! ==--------------------------------------------------------------==
! ==         VARIABLE USAGE:
! == v:      electronic density in g-space (input)                ==
! == eivps : phase factor times local pseudopotential  (vps)      ==
! == eirop : phase factor times gaussian charge distributions     ==
! ==         which replaced ionic point charges (rhops)           ==
! ==--------------------------------------------------------------==

    CALL tiset('    v1ofrho',isub)

    CALL reshape_inplace(v, (/ncpw%nhg/), vg)

    CALL zeroing(vg)!,nhg)
    ig1 = 1
    IF (geq0 .AND. (parm%ibrav .NE. 0)) THEN
       ig1 = 2
    ENDIF
    ! ==--------------------------------------------------------------==
    ener1%eht1  = 0._real_8
    IF (.NOT. first) THEN     ! rho1 - dependent stuff
       ! Hartree potential.
       ALLOCATE(aux(maxfftn),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL ffttog(drhoe,psi,aux,ncpw%nhg,.TRUE.)
       DEALLOCATE(aux,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       ! The Coulomb potential is done from the total density, up+down
       IF (geq0 .AND. (parm%ibrav .NE. 0)) THEN
          psi(1) = CMPLX(0._real_8,0._real_8,kind=real_8)! This should be satisfied anyway.
          vg(1)  = CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDIF

       ! Coulomb potential:
       DO ig=ig1,ncpw%nhg
          vg(ig) = psi(ig) * scg(ig)
          ener1%eht1   = ener1%eht1 + REAL(CONJG(vg(ig))*psi(ig))
       ENDDO
    ELSE                      ! first=.true.
       CALL zeroing(psi)!,maxfft)
    ENDIF

    IF (response1%phonon .OR. response1%tlanphon .OR. response1%tvoa) THEN
       DO ig=ig1,ncpw%nhg
          eirop1gg = eirop1(ig) * scg(ig)
          vg(ig)   = vg(ig) + eirop1gg
          ener1%eht1     = ener1%eht1 + 2._real_8*REAL(CONJG(eirop1gg) * psi(ig))
       ENDDO
    ENDIF


    ! Local external potential.
    ! TO BE ELIMINATED. Use   V1_NONLOC * PSI0   instead.
    ener1%eloc1 = 0.0_real_8
    IF (.NOT. (response1%tinteraction.OR.response1%tdummyatom)) THEN
       DO ig=ig1,ncpw%nhg
          vg(ig) = vg(ig) + v1_loc(ig)
          ener1%eloc1  = ener1%eloc1 + 2._real_8 * REAL(CONJG(psi(ig))*v1_loc(ig))
       ENDDO
    ENDIF

    ener1%eht1   =  ener1%eht1  * parm%omega
    ener1%eloc1  =  ener1%eloc1 * parm%omega

    ! ==--  v(G) -> v(R)  --------------------------------------------==
    CALL ffttor(vg,v,psi,ncpw%nhg,.TRUE.)
    ! NB: The arrays V and VG are actually pointing to identical 
    ! memory addresses.

    ! ==  XC potential  ----------------------------------------------==
    CALL v1xc_p(drhoe,v,psi)
    ! Adds the xc - contribution to V (creates spin split if 
    ! required, V(up), V(down)) and calculates also the XC-energy EXC1

    ! ==--------------------------------------------------------------==
    ! == v contains the total potential in r-space                    ==

    CALL tihalt('    v1ofrho',isub)
    RETURN
  END SUBROUTINE v1ofrho_p
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE give_scr_v1ofrho(lvofrho,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lvofrho
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: l1, l2, l4

! variables

    IF (isos1%tclust) THEN
       l1       = 2*ncpw%nhg*2+2*(nr1h+1)*nr2h*nr3pl
    ELSE
       l1       = 2*ncpw%nhg*2
    ENDIF
    CALL give_scr_v1xc_p(l2,tag)
    l4   = MAX(ncpw%nhg*2*clsd%nlsx,2*fpar%nnr1,2*ncpw%nhg*clsd%nlsd)


    lvofrho  =  MAX (l1,l2,l4) + 8
    tag      =  'v1 of rho'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_v1ofrho
  ! ==================================================================

END MODULE v1ofrho_p_utils
