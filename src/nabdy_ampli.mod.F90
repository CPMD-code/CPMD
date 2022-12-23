MODULE nabdy_ampli
  USE coor,                            ONLY: tau0,&
                                             velp
  USE cppt,                            ONLY: nzh,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nabdy_initialize,                ONLY: mom_to_vel
  USE nabdy_types,                     ONLY: nabdyvar,&
                                             nabdyvec
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE rmas,                            ONLY: rmass
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntr,&
                                             fpar,&
                                             maxsys,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nabdy_initialize_action
  PUBLIC :: nabdy_action
  PUBLIC :: nabdy_repartition_ampli
  PUBLIC :: nabdy_keep_classical

CONTAINS

  !==================================================================
  SUBROUTINE nabdy_ampli_calc
    !==--------------------------------------------------------------==
    ! loop over atoms
    ! (no atom index requires from here below)
    !  call get_grid
    !  call prepare_spline
    !  loop over atomic space (in R3)
    !    call action_fct (to get S nabla^2 S)
    !  endloop
    ! endloop
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_ampli_calc
  !==================================================================
  SUBROUTINE nabdy_initialize_action
    INTEGER                                  :: i, ia, is, k, method

    method=2
    !use spline to reconstruct the axtion from nabla s
    IF (method.EQ.1) THEN
       IF (paral%io_parent) &
            WRITE(6,*) 'nabdy_initialize_action; method not implemented'
       ! loop over atoms
       ! (no atom index requires fro here below)
       !  call get_grid
       !  call prepare_spline
       !  loop over atomic space (in R3)
       !    call action_fct (to get S nabla^2 S)
       !  endloop
       ! endloop
       ! use definition of sction and integration over the lagrangian
    ELSEIF (method.EQ.2) THEN
       !   initialize S
       !   initialize dds
       DO i=1,nabdyvar%ntrajbd
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                nabdyvec%naact(ia,is,i)=0._real_8
                DO k=1,3
                   nabdyvec%naact(ia,is,i)=nabdyvec%naact(ia,is,i)+ &
                        nabdyvec%nacoor(k,ia,is,i)*nabdyvec%namom(k,ia,is,i)
                ENDDO
                nabdyvec%ddnaact(ia,is,i)=0._real_8
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_initialize_action
  !==================================================================
  SUBROUTINE nabdy_action(rhoe,velpmem,itraj,step)
    !==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:), velpmem(:,:,:,:)
    INTEGER                                  :: itraj, step

    INTEGER                                  :: coinc(1000,4), i, ia, ia1, &
                                                ia2, is, is1, is2, j, k, &
                                                method, ncoinc, scheme
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: check_collisions
    REAL(real_8) :: disp(3), dist, distcutoff, help, mindispl, minvel, &
      v(3,2), vect(3), versor(3), vpar(3,2), vper(3,2)

    IF (ifirst.EQ.0) THEN
       CALL zeroing(nabdyvec%past_namom)!,3*maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd)
       CALL zeroing(nabdyvec%past_naact)!,maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd)
       CALL zeroing(nabdyvec%past_dds)!,maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                nabdyvec%nacoorp(k,ia,is,itraj)=nabdyvec%nacoor(k,ia,is,itraj)
                nabdyvec%past_namom(k,ia,is,itraj)=nabdyvec%namom(k,ia,is,itraj)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    !Initialization
    method=2
    scheme=1
    !if scheme .eq. 1
    mindispl=0.0001
    !if scheme .eq. 2
    minvel=0.0001

    check_collisions=.FALSE.
    distcutoff=0.05

    !use spline to reconstruct the axtion from nabla s
    IF (method.EQ.1) THEN
       ! loop over atoms
       ! (no atom index requires fro here below)
       !  call get_grid
       !  call prepare_spline
       !  loop over atomic space (in R3)
       !    call action_fct (to get S nabla^2 S)
       !  endloop
       ! endloop
       ! use definition of sction and integration over the lagrangian
    ELSEIF (method.EQ.2) THEN
       !ivano
       IF (paral%io_parent) THEN 
          IF (scheme.EQ.1) THEN
             !  ----------------------------------------------------------------
             !   update dds and dts
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   DO k=1,3
                      vect(k)=nabdyvec%namom(k,ia,is,itraj)-nabdyvec%past_namom(k,ia,is,itraj)
                      disp(k)=nabdyvec%nacoor(k,ia,is,itraj)-nabdyvec%nacoorp(k,ia,is,itraj)
                   ENDDO
                   !IF (((disp(1).LT.mindispl).AND. &
                   !     (disp(2).LT.mindispl).AND. &
                   !     (disp(3).LT.mindispl)      &
                   !     ).OR.(step.LT.100)) THEN
                   IF ((ifirst.EQ.0) .OR. (step.LE.2)) THEN

                      nabdyvec%ddnaact(ia,is,itraj)=nabdyvec%past_dds(ia,is,itraj)

                      WRITE(1111,'(I8,7E16.5)')                   &
                           itraj,(vect(k),k=1,3),(disp(k),k=1,3),  &
                           nabdyvec%ddnaact(ia,is,itraj)
                   ELSE
                      nabdyvec%ddnaact(ia,is,itraj)=vect(1)/disp(1)+       &
                                                    vect(2)/disp(2)+       &
                                                    vect(3)/disp(3)

                      WRITE(1111,'(I8,7E16.5)')                   &
                           itraj,(vect(k),k=1,3),(disp(k),k=1,3),  &
                           nabdyvec%ddnaact(ia,is,itraj)
                   ENDIF
                ENDDO
             ENDDO

          ELSEIF (scheme.EQ.2) THEN
             !  ----------------------------------------------------------------
             !   update dds and dts
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   DO k=1,3
                      vect(k)=nabdyvec%namom(k,ia,is,itraj)-nabdyvec%past_namom(k,ia,is,itraj)
                      disp(k)=nabdyvec%nacoor(k,ia,is,itraj)-nabdyvec%nacoorp(k,ia,is,itraj)
                   ENDDO

                   !IF (((ABS(nabdyvec%namom(1,ia,is,itraj)/rmass%pma(is)).LT.minvel).OR.  &
                   !     (ABS(nabdyvec%namom(2,ia,is,itraj)/rmass%pma(is)).LT.minvel).OR.  &
                   !     (ABS(nabdyvec%namom(3,ia,is,itraj)/rmass%pma(is)).LT.minvel)      &
                   !     ).OR.(step.LT.100)) THEN
                    IF ((ifirst.EQ.0) .OR. (step.LE.2)) THEN

                      nabdyvec%ddnaact(ia,is,itraj)=nabdyvec%past_dds(ia,is,itraj)

                      WRITE(1111,'(I8,7E16.5)')                                &
                           itraj,(nabdyvec%namom(k,ia,is,itraj),k=1,3),(vect(k),k=1,3),  &
                           nabdyvec%ddnaact(ia,is,itraj)
                   ELSE

                      nabdyvec%ddnaact(ia,is,itraj)=0._real_8
                      DO k=1,3
                         !         d_xd_x S= d_x p = (dp/dt) * (dt/dx) = (dp/dt) * v^(-1)
                         !         v^-1=1/3(v1^-1,v_2^-1,v_3^-1)
                         nabdyvec%ddnaact(ia,is,itraj)=nabdyvec%ddnaact(ia,is,itraj)+     &
                              (vect(k)/cntr%delt_ions)*                     &
                              (rmass%pma(is)/(3.0_real_8*nabdyvec%namom(k,ia,is,itraj)))
                      ENDDO

                      !        if (dabs(nabdyvec%ddnaact(ia,is,itraj)-nabdyvec%past_dds(ia,is,itraj)) &
                      !            .gt.100) then
                      !           nabdyvec%ddnaact(ia,is,itraj)=nabdyvec%past_dds(ia,is,itraj)
                      !        endif

                      WRITE(1111,'(I8,7E16.5)')    &
                           itraj,(nabdyvec%namom(k,ia,is,itraj),k=1,3),(vect(k),k=1,3),  &
                           nabdyvec%ddnaact(ia,is,itraj)

                   ENDIF
                ENDDO
             ENDDO
          ENDIF ! if scheme
       ENDIF ! if paral%io_parent
       ! ------------------------------------------------------------------
       ! check collisions

       IF (check_collisions) THEN
          IF (paral%io_parent) THEN
             ncoinc=0
             DO is1=1,ions1%nsp
                is2=is1
                DO ia1=1,ions0%na(is1)
                   DO ia2=ia1+1,ions0%na(is2)
                      dist=0._real_8
                      DO k=1,3
                         dist=dist+(nabdyvec%nacoor(k,ia1,is1,itraj)-       &
                              nabdyvec%nacoor(k,ia2,is2,itraj))**2
                      ENDDO
                      dist=SQRT(dist)
                      IF (dist.LT.distcutoff) THEN
                         ncoinc=ncoinc+1
                         coinc(ncoinc,1)=is1
                         coinc(ncoinc,2)=ia1
                         coinc(ncoinc,3)=is2
                         coinc(ncoinc,4)=ia2
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             DO i=1,ncoinc
                DO j=1,2
                   DO k=1,3
                      v(k,j)=(nabdyvec%namom(k,coinc(i,2*j),coinc(i,2*j-1),itraj)) &
                           /rmass%pma(coinc(i,2*j-1))
                   ENDDO
                ENDDO
                DO k=1,3
                   versor(k)=nabdyvec%nacoor(k,coinc(i,4),coinc(i,3),itraj)-     &
                        nabdyvec%nacoor(k,coinc(i,2),coinc(i,1),itraj)
                ENDDO
                help=0._real_8
                DO k=1,3
                   help=help+versor(k)**2
                ENDDO
                DO k=1,3
                   versor(k)=versor(k)/SQRT(help)
                ENDDO
                DO j=1,2
                   CALL get_components(v,j,versor,vpar,vper)
                ENDDO
                DO k=1,3
                   nabdyvec%namom(k,coinc(i,2),coinc(i,1),itraj)=       &
                        ( -vpar(k,2)+vper(k,1))*rmass%pma(coinc(i,1))
                   nabdyvec%namom(k,coinc(i,4),coinc(i,3),itraj)=       &
                        (vpar(k,1)+vper(k,2))*rmass%pma(coinc(i,3))
                ENDDO
             ENDDO
          ENDIF ! paral%io_parent
          !  bcast nabdyvec%namom and velp 
          IF (ncoinc.GE.1) THEN
             CALL mp_bcast(nabdyvec%namom,SIZE(nabdyvec%namom),parai%source,parai%allgrp)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   DO k=1,3
                      velpmem(k,ia,is,itraj)=nabdyvec%namom(k,ia,is,itraj)/rmass%pma(is)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF ! if_check_collisions

       ! ----------------------------------------------------------------
       ! update nabdyvec%past_namom and past_action
       IF (paral%io_parent) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   nabdyvec%past_namom(k,ia,is,itraj)=nabdyvec%namom(k,ia,is,itraj)
                   nabdyvec%nacoorp(k,ia,is,itraj)=nabdyvec%nacoor(k,ia,is,itraj)
                ENDDO
                nabdyvec%past_naact(ia,is,itraj)=nabdyvec%naact(ia,is,itraj)
                nabdyvec%past_dds(ia,is,itraj)=nabdyvec%ddnaact(ia,is,itraj)
             ENDDO
          ENDDO
       ENDIF ! if paral%io_parent
       ! ----------------------------------------------------------------
       !bcast
       CALL mp_bcast(nabdyvec%nacoorp,SIZE(nabdyvec%nacoorp),parai%source,parai%allgrp)
       CALL mp_bcast(nabdyvec%past_namom,SIZE(nabdyvec%past_namom),parai%source,parai%allgrp)
       CALL mp_bcast(nabdyvec%ddnaact,SIZE(nabdyvec%ddnaact),parai%source,parai%allgrp)
       CALL mp_bcast(nabdyvec%past_naact,SIZE(nabdyvec%past_naact),parai%source,parai%allgrp)
       CALL mp_bcast(nabdyvec%past_dds,SIZE(nabdyvec%past_dds),parai%source,parai%allgrp)
       !----------------------------------------------------------------
    ENDIF ! method
    !------------------------------------------------------------------
    ifirst=1
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_action
  !==================================================================
  SUBROUTINE get_components(v,j,delta,vpar,vper)
    !==--------------------------------------------------------------==
    REAL(real_8)                             :: v(3,2)
    INTEGER                                  :: j
    REAL(real_8)                             :: delta(3), vpar(3,2), vper(3,2)

    INTEGER                                  :: k
    REAL(real_8)                             :: scal

!scal = v * delta
!vpar = |v| cos(theta) * delta  ! delta is of unit length: |d|=1
!     = |v| (scal/|v| |d|) * delta
!     = scal * delta

    scal=0._real_8
    DO k=1,3
       scal=scal+delta(k)*v(k,j)
    ENDDO
    DO k=1,3
       vpar(k,j)=delta(k)*scal
    ENDDO
    DO k=1,3
       IF (j.EQ.1) THEN
          vper(k,j)=v(k,j)-vpar(k,j)
       ELSE
          vper(k,j)=v(k,j)+vpar(k,j)
       ENDIF
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_components
  !==================================================================
  SUBROUTINE nabdy_ppener(eii,eei,eps,rhoe,eivps,eirop)
    !==--------------------------------------------------------------==
    !==                        computes                              ==
    !== electrostatic and pseudopotential energies and the           ==
    !== potential energy contributions to the forces on the ions     ==
    !==--------------------------------------------------------------==
    !== input:                                                       ==
    !==   rhoe  electronic density in real space
    !==   eivps phase factor times local pseudopotential  (vps)      ==
    !==   eirop phase factor times gaussian charge distributions     ==
    !==         which replaced ionic point charges (rhops)           ==
    !== output:                                                      ==
    !==   eii   gaussian charge dist. (rhops) (ion-ion interaction)  ==
    !==   eei   e-e interaction                                      ==
    !==         both used by diagonalisation schemes (no g=0 term)   ==
    !==   eps   integration of electronic density*local pp           ==
    !==         epseu=2._real_8*dreal(eps)*omega)                         ==
    !==   vploc g=0 pp part (given from pp data) (energy=nel*vploc)  ==
    !==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eii, eei, eps
    REAL(real_8)                             :: rhoe(:)
    COMPLEX(real_8)                          :: eivps(ncpw%nhg), &
                                                eirop(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'nabdy_ppener'

    COMPLEX(real_8)                          :: rhet, rhog, rhoi, rp, vp
    COMPLEX(real_8), ALLOCATABLE             :: vect(:)
    INTEGER                                  :: ierr, ig, ig1, ir

    ALLOCATE(vect(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vect)!,maxfft)
    !==--------------------------------------------------------------==
    !transform the density to g space
    DO ir=1,fpar%nnr1
       vect(ir)=CMPLX(rhoe(ir),0._real_8)
    ENDDO
    CALL fwfftn(vect,.FALSE.,parai%allgrp)

    IF(geq0) THEN
       vp=eivps(1)
       vploc=REAL(vp,kind=real_8)
       eps=0.5_real_8*vp*CONJG(vect(nzh(1)))
       ig1=2
       rhoi=eirop(1)
       rhet=vect(nzh(1))
       rhog=rhet+rhoi
       eei=0.5_real_8*scg(1)*REAL(rhog,kind=real_8)*REAL(rhoi,kind=real_8)
       !vw rp never set before, HUGEified
       rp = HUGE(0.0_real_8)
       eii=0.5_real_8*scg(1)*rp*rp
    ELSE
       ig1=1
       eps=(0.0_real_8,0.0_real_8)
       vploc=0.0_real_8
       eei=(0.0_real_8,0.0_real_8)
       eii=(0.0_real_8,0.0_real_8)
    ENDIF
    DO ig=ig1,ncpw%nhg
       !  rho e
       rhet=vect(nzh(ig))
       !  rho i
       rhoi=eirop(ig)
       !  rho e+i
       rhog=rhet+rhoi

       eei=eei+(scg(ig)*rhog)*CONJG(rhoi)
       eii=eii+scg(ig)*rhoi*CONJG(rhoi)
       eps=eps+CONJG(rhet)*eivps(ig)
    ENDDO
    !
    DEALLOCATE(vect,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_ppener
  !==================================================================
  SUBROUTINE nabdy_redistribute_ampli
    !==--------------------------------------------------------------==
    !Variables
    !call nabdy_check_low_ampli(n_points_to_redistr)
    !do i=1,n_points_to_redistr
    ! call nabdy_find_ampli_com
    ! call nabdy_put_grid
    ! call nabdy_find_fe_pos
    ! call nabdy_distr_ampli
    ! call nabdy_normalize
    !enddo
    !==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE nabdy_redistribute_ampli
  !==================================================================
  SUBROUTINE nabdy_repartition_ampli
    !==--------------------------------------------------------------==
    !Variables
    CHARACTER(*), PARAMETER :: procedureN = 'nabdy_repartition_ampli'

    INTEGER                                  :: i, ia, ierr, is, itraj, k, &
                                                n_points
    INTEGER, ALLOCATABLE                     :: index_rm(:), index_split(:)
    LOGICAL                                  :: contact
    REAL(real_8)                             :: coor(2,3), d, d2, factor, &
                                                min_ampli, xn(3)
    REAL(real_8), ALLOCATABLE                :: ndelta(:)

    IF (paral%io_parent) THEN
       ALLOCATE(index_split(nabdyvar%ntrajbd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(index_rm(nabdyvar%ntrajbd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ndelta(ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       min_ampli=1.0e-6_real_8
       contact=.FALSE.

       CALL nabdy_check_ampli(n_points,index_split,index_rm,ndelta)

       IF (n_points.GT.0) THEN
          WRITE(6,'(A,2X,I5)')'Number amplitude redistributions:',n_points
       ENDIF

       DO i=1,n_points

          DO is=1,ions1%nsp
             IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
             DO ia=1,ions0%na(is)
                !    -------------------------------
                factor=1.0
111             CONTINUE
                DO k=1,3
                   xn(k)=repprngu()
                ENDDO
                d=SQRT(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
                DO k=1,3
                   xn(k)=xn(k)/d
                ENDDO
                d2=ndelta(is)/(2.0*factor)
                DO k=1,3
                   coor(1,k)=nabdyvec%nacoor(k,ia,is,index_split(i))+d2*xn(k)
                   coor(2,k)=nabdyvec%nacoor(k,ia,is,index_split(i))-d2*xn(k)
                ENDDO
                CALL check_contacts(coor,is,ia,contact,d2)
                IF (contact) THEN
                   factor=factor+0.01
                   GOTO 111
                ENDIF
                !    --------------------------------
                DO k=1,3
                   nabdyvec%nacoor(k,ia,is,index_rm(i))=coor(1,k)
                   nabdyvec%nacoor(k,ia,is,index_split(i))=coor(2,k)
                   nabdyvec%nacoorp(k,ia,is,index_rm(i))= &
                        nabdyvec%nacoorp(k,ia,is,index_split(i))+d2*xn(k)
                   nabdyvec%nacoorp(k,ia,is,index_split(i))= &
                        nabdyvec%nacoorp(k,ia,is,index_split(i))-d2*xn(k)
                   nabdyvec%namom(k,ia,is,index_rm(i))=nabdyvec%namom(k,ia,is,index_split(i))/SQRT(2._real_8)
                   nabdyvec%namom(k,ia,is,index_split(i))= &
                        nabdyvec%namom(k,ia,is,index_split(i))/SQRT(2._real_8)
                ENDDO
                nabdyvec%naampl(ia,is,index_rm(i))=nabdyvec%naampl(ia,is,index_split(i))/SQRT(2._real_8)
                nabdyvec%naampl(ia,is,index_split(i))=nabdyvec%naampl(ia,is,index_split(i))/SQRT(2._real_8)
                nabdyvec%naampl_temp(ia,is,index_rm(i))= &
                     nabdyvec%naampl_temp(ia,is,index_split(i))/SQRT(2._real_8)
                nabdyvec%naampl_temp(ia,is,index_split(i))= &
                     nabdyvec%naampl_temp(ia,is,index_split(i))/SQRT(2._real_8)
                nabdyvec%naact(ia,is,index_rm(i))=nabdyvec%naact(ia,is,index_split(i))/2.0
                nabdyvec%naact(ia,is,index_split(i))=nabdyvec%naact(ia,is,index_split(i))/2.0

                nabdyvec%ddnaact(ia,is,index_rm(i))=nabdyvec%ddnaact(ia,is,index_split(i))
                nabdyvec%naomega(ia,is,index_rm(i))=nabdyvec%naomega(ia,is,index_split(i))
                nabdyvec%nanormalf(ia,is,index_rm(i))=nabdyvec%nanormalf(ia,is,index_split(i))
             ENDDO
          ENDDO

       ENDDO

       DEALLOCATE(index_split,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(index_rm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ndelta,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

    ENDIF

    CALL mp_bcast(n_points,parai%source,parai%allgrp)
    IF (n_points.EQ.0) RETURN

    !bcast
    CALL mp_bcast(nabdyvec%nacoor,SIZE(nabdyvec%nacoor),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%nacoorp,SIZE(nabdyvec%nacoorp),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%namom,SIZE(nabdyvec%namom),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%naampl,SIZE(nabdyvec%naampl),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%naampl_temp,SIZE(nabdyvec%naampl_temp),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%naact,SIZE(nabdyvec%naact),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%ddnaact,SIZE(nabdyvec%ddnaact),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%naomega,SIZE(nabdyvec%naomega),parai%source,parai%allgrp)
    CALL mp_bcast(nabdyvec%nanormalf,SIZE(nabdyvec%nanormalf),parai%source,parai%allgrp)

    !and transfer nabdyvec%nacoor and nabdyvec%namom to TAU0 and VELP
    DO itraj=1,nabdyvar%ntrajbd
       CALL dcopy(3*maxsys%nax*maxsys%nsx,nabdyvec%nacoor(1,1,1,itraj),1,tau0(1,1,1),1)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,nabdyvec%namom(1,1,1,itraj),1,velp(1,1,1),1)
       CALL mom_to_vel(velp)
       ! call dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,velpmem(1,1,1,itraj),1)
    ENDDO
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_repartition_ampli
  !==================================================================
  SUBROUTINE  nabdy_keep_classical(velpmem)
    !==--------------------------------------------------------------==
    REAL(real_8)                             :: velpmem(:,:,:,:)

    INTEGER                                  :: i, ia, is, k

!Variables
!real(real_8)     ::     nacoor_mem(3,maxsys%nax,maxsys%nsx,*),namom_mem(3,maxsys%nax,maxsys%nsx,*)
!pointer    (ip_nacoor_mem,nacoor_mem),(ip_namom_mem,namom_mem)
!integer          ::    icall
!data       icall/0/
!save       icall
!save       ip_nacoor_mem,ip_namom_mem
!if (icall.eq.0) then
! call memory(ip_nacoor_mem,3*maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd,'nacoor_mem')
! call memory(ip_namom_mem,3*maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd,'namom_mem')
! do is=1,ions1%nsp
!   do ia=1,ions0%na(is)
!     do i=1,nabdyvar%ntrajbd
!      do k=1,3
!       nacoor_mem(k,ia,is,i)=nacoor(k,ia,is,i)
!       namom_mem(k,ia,is,i)=namom(k,ia,is,i)
!      enddo
!     enddo
!   enddo
! enddo
!else

    DO is=1,ions1%nsp
       IF (ions0%iatyp(is).LE.nabdyvar%nabdy_zmax) CYCLE
       DO ia=1,ions0%na(is)
          DO i=2,nabdyvar%ntrajbd
             DO k=1,3
                !       nacoor(k,ia,is,i)=nacoor_mem(k,ia,is,i)
                !       namom(k,ia,is,i)=namom_mem(k,ia,is,i)
                nabdyvec%namom(k,ia,is,i)=nabdyvec%namom(k,ia,is,1)
                velpmem(k,ia,is,i)=velpmem(k,ia,is,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !endif
    !icall=1
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE  nabdy_keep_classical
  !==================================================================
  SUBROUTINE  check_contacts(coor,is,ia,contact,d2)
    !==--------------------------------------------------------------==
    REAL(real_8)                             :: coor(2,*)
    INTEGER                                  :: is, ia
    LOGICAL                                  :: contact
    REAL(real_8)                             :: d2

    INTEGER                                  :: i, j
    REAL(real_8)                             :: dd, mindd

!Variables

    mindd=100.0
    contact=.FALSE.

    DO i=1,nabdyvar%ntrajbd
       DO j=i+1,nabdyvar%ntrajbd
          dd=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 +  &
               (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 +  &
               (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2)
          IF (dd.LT.mindd)  mindd=dd
       ENDDO
    ENDDO

    IF (mindd.LT.d2) contact=.TRUE.

    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE  check_contacts
  !==================================================================
  SUBROUTINE nabdy_check_ampli(n_points,index_split,index_rm,ndelta)
    !==--------------------------------------------------------------==
    INTEGER                                  :: n_points, index_split(*), &
                                                index_rm(*)
    REAL(real_8)                             :: ndelta(*)

    CHARACTER(*), PARAMETER :: procedureN = 'nabdy_check_ampli'

    INTEGER                                  :: i, ia, ierr, is, j
    REAL(real_8)                             :: avampli, dd, min_ampli, &
                                                min_temp, mindd
    REAL(real_8), ALLOCATABLE                :: totampli(:)

!Variables

    ALLOCATE(totampli(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    !determine average ampli
    DO i=1,nabdyvar%ntrajbd
       totampli(i)=1._real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             totampli(i)=totampli(i)*nabdyvec%naampl(ia,is,i)
          ENDDO
       ENDDO
    ENDDO

    !find the molecular amplitudes larger than 2.0_real_8*avampli
    !they are in number of n_points

    n_points=0
    avampli=0.0_real_8
    DO i=1,nabdyvar%ntrajbd
       avampli=avampli+totampli(i)
    ENDDO
    avampli=avampli/nabdyvar%ntrajbd
    DO i=1,nabdyvar%ntrajbd
       IF ((totampli(i).GT.(2.0_real_8*avampli))              &
            .OR.(totampli(i).GT.0.7)) THEN
          n_points=n_points+1
          index_split(n_points)=i
       ENDIF
       !if (paral%io_parent) write(6,*) 'xxxxx',i,totampli(i),3.0_real_8*avampli
    ENDDO

    !find the 'n_points' minima 

    IF (n_points.GT.0) THEN
       DO i=1,n_points
          IF (i.EQ.1) THEN
             min_ampli=0._real_8
          ELSE
             min_ampli=totampli(index_rm(i-1))
          ENDIF

          min_temp=100.0
          DO j=1,nabdyvar%ntrajbd 
             IF ((totampli(j).LT.min_temp).AND.(totampli(j).GT.min_ampli)) THEN
                index_rm(i)=j
                min_temp=totampli(j)
             ENDIF
          ENDDO

       ENDDO
    ENDIF

    !find minimum distance between fe atoms (per species)
    DO is=1,ions1%nsp

       mindd=1000.0
       DO ia=1,ions0%na(is)
          DO i=1,nabdyvar%ntrajbd
             DO j=i+1,nabdyvar%ntrajbd
                dd=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 +  &
                     (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 +  &
                     (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2)
                IF (dd.LT.mindd) THEN
                   mindd=dd
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       ndelta(is)=mindd
    ENDDO

    DEALLOCATE(totampli,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_check_ampli
  !==================================================================

END MODULE nabdy_ampli
