MODULE nabdy_initialize
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fbohr,&
                                             pi
  USE coor,                            ONLY: tau0,&
                                             velp
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             cpat,&
                                             cpsp,&
                                             mm_go_mm,&
                                             mm_revert,&
                                             nat_grm
  USE mp_interface,                    ONLY: mp_bcast
  USE nabdy_types,                     ONLY: nabdyfric,&
                                             nabdyvar,&
                                             nabdyvec
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nabdy_mem, nabdy_rel, nabdy_load, mom_to_vel, vel_to_mom


CONTAINS

  !==================================================================
  SUBROUTINE nabdy_mem
    !==--------------------------------------------------------------==
    !call memory for all arrays in nabdy.inc

    CHARACTER(*), PARAMETER                  :: procedureN = 'nabdy_mem'

    INTEGER                                  :: ierr

    ALLOCATE(nabdyvec%nacoor(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%nacoorp(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%namom(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%nafor(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naampl(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naqp_pat(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naampl_temp(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naomega(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naact(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%ddnaact(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%gridshep(3,1000),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%ddelta(3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%nanormalf(maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%nabdy_mean(maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%nabdy_dx(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%nabdy_dx2(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%nabdy_dp(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%nabdy_dp2(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%natotqp(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%natotcl(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%natotk(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%naampfe(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%dispx_ref(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyfric%dispp_ref(maxsys%nsx,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%past_namom(3,maxsys%nax,maxsys%nax,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%past_naact(maxsys%nax,maxsys%nax,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nabdyvec%past_dds(maxsys%nax,maxsys%nax,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    !call memory arrays for spline done in NABDY_SPLINE_MEM

    CALL nabdy_azzero_at_init
    !==--------------------------------------------------------------==
    RETURN 
  END SUBROUTINE nabdy_mem
  !==================================================================
  SUBROUTINE nabdy_azzero_at_init
    !==--------------------------------------------------------------==
    CALL zeroing(nabdyfric%nabdy_dx)!,maxsys%nsx*maxsys%nax)
    CALL zeroing(nabdyfric%nabdy_dx2)!,maxsys%nsx*maxsys%nax)  
    CALL zeroing(nabdyfric%nabdy_dp)!,maxsys%nsx*maxsys%nax)
    CALL zeroing(nabdyfric%nabdy_dp2)!,maxsys%nsx*maxsys%nax)
    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_azzero_at_init
  !==================================================================
  SUBROUTINE nabdy_rel
    !==--------------------------------------------------------------==
    !call memory for all arrays in nabdy.inc

    CHARACTER(*), PARAMETER                  :: procedureN = 'nabdy_rel'

    INTEGER                                  :: ierr

    DEALLOCATE(nabdyvec%nacoor,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%nacoorp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%namom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%nafor,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naampl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naqp_pat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naampl_temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naomega,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naact,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%ddnaact,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%gridshep,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%ddelta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%nanormalf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%nabdy_mean,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%nabdy_dx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%nabdy_dx2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%nabdy_dp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%nabdy_dp2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%natotqp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%natotcl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%natotk,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%naampfe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%dispx_ref,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyfric%dispp_ref,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%past_namom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%past_naact,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nabdyvec%past_dds,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    !==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_rel
  !==================================================================
  SUBROUTINE nabdy_load(temp_inst)
    !==--------------------------------------------------------------==
    !To be called in nabby_md.F because it requires the temperature
    !If temp_inst.eq.0 then it takes tempw from system.h

    REAL(real_8), INTENT(in)                 :: temp_inst

    CHARACTER(*), PARAMETER                  :: procedureN = 'nabdy_load'

    CHARACTER(LEN=14)                        :: filen1, filen4
    CHARACTER(LEN=28)                        :: filen2, filen3
    INTEGER                                  :: i, ia, ierr, is, j, jj, k
    INTEGER, ALLOCATABLE                     :: gr_iat(:)
    LOGICAL                                  :: fexist, first, status
    REAL(real_8)                             :: a, cutdist_h, cutdist_l, &
                                                dist, mass, mindist, s, temp
    REAL(real_8), ALLOCATABLE                :: ampli(:), coord(:,:), &
                                                lomega(:), veloc(:,:)

!Variables
!
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(coord(3,ions1%nat*nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(veloc(3,ions1%nat*nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ampli(ions1%nat*nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lomega(ions1%nat*nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !
    !from geofile.F
    i=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          i=i+1
          gr_iat(nat_grm(i))=ions0%iatyp(is)
       ENDDO
    ENDDO
    !
    temp=temp_inst
    IF (temp.EQ.0.0) temp=cntr%tempw
    IF (temp.EQ.0.0) temp=300
    cutdist_h=0.001
    cutdist_l=0.01

    filen1="nabdy_geometry"
    filen2="nabdy_geometry.xyz"
    filen3="nabdy_gaussian.dat"
    filen4='nabdy_time.dat'

    IF (paral%io_parent) THEN
       INQUIRE(file=filen1,exist=fexist)
       IF(.NOT.fexist) THEN
          OPEN(unit=777,file=filen1,status='new')
       ELSE
          OPEN(unit=777,file=filen1,status='old')
       ENDIF
    ENDIF

    IF (paral%io_parent.AND..NOT.fexist) THEN
       ! First run:
       ! copy all coordinates in nabdyvec%nacoor(3,nax,nsx,1)
       ! spread the remaning ntrajbd-1 coordinates
       ! the same for the velocities


       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)

             DO k=1,3
                nabdyvec%nacoor(k,ia,is,1)=tau0(k,ia,is)
                nabdyvec%namom(k,ia,is,1)=velp(k,ia,is)
                coord(k,1)=tau0(k,ia,is)
                veloc(k,1)=velp(k,ia,is)
             ENDDO
             mass=rmass%pma(is)

             CALL distribute_atoms(is,mass,temp,coord,veloc,a,s,.TRUE.)
             nabdyvec%naampl(ia,is,1)=a
             nabdyvec%naomega(ia,is,1)=s

             mindist=HUGE(0.0_real_8)
             DO j=2,nabdyvar%ntrajbd
                !     call distribute_atoms(is,mass,temp,coord,veloc,a,s,first)
100             CONTINUE
                WRITE(*,*) 'j, mindist',j,mindist
                CALL distribute_atoms(is,mass,temp,coord,veloc,a,s,.FALSE.)
                mindist=1000._real_8
                DO jj=1,j-1
                   dist=0._real_8
                   DO k=1,3
                      dist=dist+(coord(k,2)-nabdyvec%nacoor(k,ia,is,jj))**2
                   ENDDO
                   dist=SQRT(dist)
                   IF (dist.LT.mindist) mindist=dist
                ENDDO
                IF ( ((ions0%iatyp(is).LE.nabdyvar%nabdy_zmax) .AND. &
                      (mindist.LT.cutdist_l))                        .OR. & 
                     ((ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) .AND. &
                      (mindist.LT.cutdist_h))) THEN
                   ! refuse and retry
                   GOTO 100
                ELSE
                   ! accept and go ahead
                ENDIF
                DO k=1,3
                   nabdyvec%nacoor(k,ia,is,j)=coord(k,2)
                   nabdyvec%namom(k,ia,is,j)=veloc(k,2)*mass
                ENDDO
                nabdyvec%naampl(ia,is,j)=a
                nabdyvec%naomega(ia,is,j)=s
             ENDDO
          ENDDO
       ENDDO
       ! write files nabdy_geometry and nabdy_geometry.xyz
       OPEN(unit=888,file=filen2,status='new')
       OPEN(unit=999,file=filen3,status='new')
       WRITE(888,'(I6)') ions1%nat*nabdyvar%ntrajbd
       WRITE(888,'(A8,I5)') 'GEOMETRY'
       DO j=1,nabdyvar%ntrajbd
          !  write(888,'(A8,I5)') 'GEOMETRY',j
          !  from geofile.F: nabdyvec%nacoor > coord
          i=0
          DO is=1,ions1%nsp
             mass=rmass%pma(is)
             DO ia=1,ions0%na(is)
                i=i+1
                DO k=1,3
                   coord(k,nat_grm(i))=nabdyvec%nacoor(k,ia,is,j)
                   veloc(k,nat_grm(i))=nabdyvec%namom(k,ia,is,j)/mass
                ENDDO
                ampli(nat_grm(i))=nabdyvec%naampl(ia,is,j)
                lomega(nat_grm(i))=nabdyvec%naomega(ia,is,j)
             ENDDO
          ENDDO
          DO i=1,ions1%nat
             WRITE(888,'(a3,3f20.12,8x,3f20.12)') elem%el(gr_iat(i)),&
                  (coord(k,i)/fbohr,k=1,3),(veloc(k,i)/fbohr,k=1,3)
             WRITE(777,'(3f20.12,8x,3f20.12)') &
                  (coord(k,i),k=1,3),(veloc(k,i),k=1,3)
             WRITE(999,'(2f20.12)') ampli(i),lomega(i)
          ENDDO
       ENDDO
       OPEN(unit=666,file=filen4,status='new')
       nabdyvar%nabdy_time=0._real_8
       WRITE(666,'(F12.4)') nabdyvar%nabdy_time
       CLOSE(666)
       CLOSE(888)
       CLOSE(999)
    ELSEIF (paral%io_parent.AND.fexist) THEN
       ! is a continuation run
       ! read all coordinates and velocities from nabdy_geometry
       OPEN(unit=999,file=filen3,status='old')
       OPEN(unit=666,file=filen4,status='old')
       DO j=1,nabdyvar%ntrajbd
          IF(cntl%tqmmm)THEN
             DO i=1,ions1%nat
                is=cpsp(i)
                mass=rmass%pma(is)
                READ(777,fmt=*)(coord(k,i),k=1,3),(veloc(k,i),k=1,3)
                READ(999,fmt=*) ampli(i),lomega(i)
                DO k=1,3
                   nabdyvec%nacoor(k,cpat(i),cpsp(i),j)=coord(k,i)+clsaabox%mm_c_trans(k)
                   nabdyvec%namom(k,cpat(i),cpsp(i),j)=veloc(k,i)*mass
                END DO
                nabdyvec%naampl(cpat(i),cpsp(i),j)=ampli(i)
                nabdyvec%naomega(cpat(i),cpsp(i),j)=lomega(i)
             ENDDO
          ELSE
             DO is=1,ions1%nsp
                mass=rmass%pma(is)
                DO ia=1,ions0%na(is)
                   READ(777,fmt=*) (nabdyvec%nacoor(k,ia,is,j),k=1,3),(nabdyvec%namom(k,ia,is,j),k=1,3)
                   READ(999,fmt=*) nabdyvec%naampl(ia,is,j),nabdyvec%naomega(ia,is,j)
                   DO k=1,3
                      nabdyvec%namom(k,ia,is,j)=nabdyvec%namom(k,ia,is,j)*mass
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       CLOSE(999)
       READ(666,fmt=*) nabdyvar%nabdy_time
       CLOSE(666)
    ENDIF
    !
    IF (paral%io_parent) THEN 
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO i=1,nabdyvar%ntrajbd
                nabdyvec%nanormalf(ia,is,i)=(2.0*pi)**(1.5) * (nabdyvec%naomega(ia,is,i))**3
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !
    IF (paral%io_parent) CLOSE(777)
    !
    !bcast to all nodes
    CALL mp_bcast(nabdyvec%nacoor,SIZE(nabdyvec%nacoor),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvec%namom,SIZE(nabdyvec%namom),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvec%naampl,SIZE(nabdyvec%naampl),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvec%naomega,SIZE(nabdyvec%naomega),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvec%nanormalf,SIZE(nabdyvec%nanormalf),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvar%nabdy_time,parai%io_source,parai%cp_grp)
    !==--------------------------------------------------------------==
    !format of the file nabdy_geometry.xyz
    !NATOMS
    !GEOMETRY          1
    !ATOM NAME      X   Y   Z   V1   V2   V3
    !....
    !GEOMETRY          2
    !ATOM NAME     X   Y   Z   V1   V2   V3
    !....
    !GEOMETRY          3
    !...
    !
    !==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,status)
    DEALLOCATE(coord,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)    
    DEALLOCATE(veloc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ampli,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lomega,STAT=ierr) 
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_iat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE nabdy_load
  !==================================================================
  SUBROUTINE distribute_atoms(is,mass,temp,coord,veloc,ampli,sigma,first)
    !==--------------------------------------------------------------==
    !use http://hyperphysics.phy-astr.gsu.edu/hbase/quantum/debrog2.html
    !with mass of the element and kinetic energy corresponding to the 
    !temperature T
    !KB     = 8.6173324 E-5 eV/K [Wikipedia]
    !3/2 KB = .000129259986 eV/K
    !mp c^2 = 9.382734 E8 eV
    !hc     = 12398.4 eV A
    !mp     = proton mass
    !lambda[A]= hc/pc
    !"Particles can be described by wave packets. The particle velocity
    ! equals the group velocity of the wavepacket and the particle can 
    ! be localized not better than within the de Broglie wavelength"
    ! Atoms, Molecules and Photons: An Introduction to Atomic-, Molecular-, 
    ! and  Quantum Physics (second edition) by Wolfgang Demtr\"oder, 
    ! Springer Verlag, pp 104-105.
    !
    INTEGER                                  :: is
    REAL(real_8)                             :: mass, temp, coord(3,*), &
                                                veloc(3,*), ampli, sigma
    LOGICAL                                  :: first

    REAL(real_8), PARAMETER                  :: hc = 12398.4 , &
                                                kb32 = 0.000129259986_real_8, &
                                                mpc2 = 9.382734e8_real_8

    INTEGER                                  :: k
    LOGICAL                                  :: homog_distr
    REAL(real_8)                             :: dist2, ke, lambda, lomega, &
                                                mc2, pc, segno, sigma_p, &
                                                sigma_x, x0, x1, zi1, zi2

!vw real(real_8), dimension(:,:), intent(inout) :: coord,veloc
!Variables
! [eV A]

    homog_distr=.TRUE.

    ke=kb32*temp
    mc2=mpc2*mass
    pc=SQRT(2.0_real_8*ke*mc2)   
    lambda=hc/pc * 1.889726
    WRITE(6,*) 'lambda',lambda,ke,mc2,pc,temp
    !1.889726 = conversion factor A -> au
    !generate random numbers distributed according to a gaussian
    !distribution: M.P. Allem and D.J. Tildesley, "Computer Simulation 
    !of Liquids", Appendix G.3, p347. 
    !For a proton at room temperature lamda=.14682328 nm in agreement
    !with the value in Cohen-Tannoudji, diu, laloe, p42, lambda = 1.4 A
    WRITE(6,*) 'mass',mass
    IF (mass.LT.2000) THEN
       ! sigma_x=lambda/1.0
       sigma_x=lambda*2.0
    ELSE
       sigma_x=lambda
    ENDIF
    ! the sqrt(variance) of a SINGLE gaussian is setted to 1/n (n=2,...,10) of
    ! the value of sigma_x 
    !NORMALIZED GAUSSIAN in 3D:
    ! 1/(sigma^3 (2 pi)^3/2) exp^(-1/2 x^2/(sigma^2))
    IF (mass.LT.2000) THEN
       !  sigma=sigma_x/1.5 !  used for H2
       sigma=sigma_x/2.0
    ELSE
       sigma=sigma_x/2.0
    ENDIF
    lomega=(1._real_8/(2.0*sigma**2))
    IF (paral%io_parent) WRITE(*,*) 'sigma, omega',sigma,lomega

    IF (first) THEN
       ! FE at the central position (maximal amplitude for it)
       DO k=1,3
          coord(k,2)=coord(k,1)
          veloc(k,2)=veloc(k,1)
       ENDDO
       IF (.NOT.homog_distr) THEN
          ampli=1.0/((sigma_x)**3*(2.0*PI)**1.5_real_8)
       ELSE
          ampli=sqrt(1.0_real_8/nabdyvar%ntrajbd)
       ENDIF
    ELSE
       IF (.NOT.homog_distr) THEN
          DO k=1,3
             segno=1._real_8
             zi1=repprngu()
             zi2=repprngu()
             IF (zi2.GT.0.5) segno=-1._real_8
             x0=coord(k,1)
             coord(k,2)=x0+segno*2*sigma_x*zi1
          ENDDO
          dist2=0.0
          DO k=1,3
             dist2=dist2+(coord(k,2)-coord(k,1))**2
          ENDDO
          ampli=1.0/((sigma_x)**3*(2.0*PI)**1.5_real_8)*& 
               EXP((-0.5_real_8 * dist2)/(sigma_x**2))
       ELSE
          DO k=1,3
             zi1=repprngu()
             zi2=repprngu()
             x0=coord(k,1)
             x1=(-2._real_8*LOG(zi1))**0.5 * COS(2._real_8*pi*zi2)
             !   x2=(-2._real_8*LOG(zi1))**0.5 * dsin(2._real_8*pi*zi2)
             coord(k,2)=x0+sigma_x*x1
          ENDDO
          ampli=sqrt(1.0_real_8/nabdyvar%ntrajbd)
       ENDIF
       ! distribute velocities
       IF (ions0%iatyp(is).LE.nabdyvar%nabdy_zmax) THEN
          DO k=1,3
             sigma_p=0.1_real_8 * veloc(k,1)  !1/(2*sigma_x)
             zi1=repprngu()
             zi2=repprngu()
             x0=veloc(k,1)
             x1=(-2._real_8*LOG(zi1))**0.5 * COS(2._real_8*pi*zi2)
             !   x2=(-2._real_8*LOG(zi1))**0.5 * dsin(2._real_8*pi*zi2)
             veloc(k,2)=x0+sigma_p*x1
          ENDDO
       ELSE
          DO k=1,3
             veloc(k,2)=veloc(k,1)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN 
  END SUBROUTINE distribute_atoms
  ! ==================================================================
  SUBROUTINE mom_to_vel(vect)
    REAL(real_8)                             :: vect(:,:,:)

    INTEGER                                  :: ia, is, k

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             vect(k,ia,is)=vect(k,ia,is)/rmass%pma(is)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mom_to_vel
  ! ==================================================================
  SUBROUTINE vel_to_mom(vect,itraj)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vect(:,:,:,:)
    INTEGER                                  :: ITRAJ

    INTEGER                                  :: ia, is, k

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             vect(k,ia,is,itraj)=vect(k,ia,is,itraj)*rmass%pma(is)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vel_to_mom
  ! ==================================================================

END MODULE nabdy_initialize

