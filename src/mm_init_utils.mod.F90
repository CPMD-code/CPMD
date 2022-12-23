MODULE mm_init_utils
  USE adat,                            ONLY: elem
  USE cell,                            ONLY: cell_com
  USE coor,                            ONLY: lvelini,&
                                             tau0,&
                                             velp
  USE cotr,                            ONLY: duat
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE ions,                            ONLY: ions0
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: &
       clsaabox, cpat, cpsp, gratom, inc_l, mm_charge, mm_stat, mmdim, nam, &
       naq, nat_cpmd, nat_grm, ncag_l, solsolv, solvvv
  USE mm_input,                        ONLY: addh,&
                                             capp,&
                                             clc,&
                                             eqm_r,&
                                             excl_comm,&
                                             g96_vel,&
                                             lqmmm,&
                                             solqmmm
  USE mm_ion_dens,                     ONLY: mm_raggio
  USE mm_parallel,                     ONLY: gparai,&
                                             gparal
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE rmas,                            ONLY: rmass
  USE store_types,                     ONLY: cprint
  USE system,                          ONLY: maxsys
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_init

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_init
    ! ==--------------------------------------------------------------==

#if defined (__GROMOS)
    ! McB
#endif      
    IMPLICIT NONE

    ! QM/MM. see. below for compatibility init for normal QM runs.
#if defined (__GROMOS)
    ! McB
    INTEGER :: npm,nsm,ntx,ig,ntxo
    COMMON /mstart/ npm,nsm,ntx,ig,ntxo
    ! locals
    INTEGER :: NSX_s
    INTEGER :: i,ia,is,iis,nrtot
    INTEGER :: diffe1,diffe2
    REAL(real_8) :: xm(3)
    REAL(real_8) :: diffevol
    INTEGER :: msglen,ierr
    LOGICAL :: ferror
    CHARACTER(*),PARAMETER::procedureN='mm_init'
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==

    mm_stat=.TRUE.

    IF (.NOT. lqmmm%qmmm) THEN
       CALL mm_compat_init
       RETURN
    ENDIF

    NSX_s=maxsys%nsx
    CALL allocate_gromos

    IF (gparal%mmnode) THEN
       CALL mm_setup(solsolv%nrpt,solsolv%nsolv)
       IF (solqmmm%tflexsol)THEN
          CALL mm_flex_solv(solsolv%nrpt,solsolv%nsolv)
       ENDIF
       IF (capp%cap_h .OR. addh%n_carbon.GT.0)THEN
          CALL mm_add_dummy(solsolv%nrpt,solsolv%nsolv)
       ENDIF
       IF (addh%n_carbon.GT.0)THEN
          CALL mm_add_hydrogen
       ENDIF
       IF (capp%cap_h)THEN
          CALL mm_cap_h( solsolv%nrpt, solsolv%nsolv)
       ENDIF
       nrtot=solsolv%nrpt+solsolv%nsolv
       ! define new maxsys%nax and maxsys%nsx
       CALL mm_get_nsx(mmdim%nspq,solsolv%nrpt,solsolv%nsolv)
       IF (.NOT.clc%classical) THEN
          CALL mm_quantum_topo(mmdim%nspq,solsolv%nrpt,solsolv%nsolv)
       ENDIF
    ENDIF
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL mp_bcast(maxsys%nax,gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast_byte(maxsys, size_in_bytes_of(maxsys),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast_byte(solsolv, size_in_bytes_of(solsolv),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(solsolv%nsolv,gparai%mmsource,parai%qmmmgrp)
    nrtot=solsolv%nrpt+solsolv%nsolv
    ! FIXME: AK 2005/05/24. memory allocates real(8) :: words, if we would
    ! use (NRTOT+1)/IRAT to significantly reduce the memory requirements for 
    ! all integer arrays based on NRTOT on systems with IRAT=2 (i.e. most current
    ! platforms). this will be especially useful for systems with a large MM part.
    ALLOCATE(cpat(nrtot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cpsp(nrtot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(NAm(maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(NAq(maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! FIMXE: AK 2005/05/24 gratom can become especially big.
    ! we should try to get rid of it. it is rarely used, too.
    ! -> mm_nlist.F/mm_short_range_classic.F
    ALLOCATE(gratom(maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nat_cpmd(0:nrtot+duat%ndat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO implement 0: though memory90
    ALLOCATE(nat_grm(0:nrtot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO implement 0: through memory90
    IF (gparal%mmnode) THEN
       CALL mm_detit(cpat,cpsp,nat_cpmd,nat_grm,gratom,&
            NAq,mmdim%nspq,NAm,mmdim%nspm,mmdim%natq,mmdim%naxq,mmdim%natm,solsolv%nrpt,solsolv%nsolv)
    ENDIF
    ! check wether the exclusions are valid
    ferror=.FALSE.
    DO i=1,excl_comm%nce
       ia=excl_comm%atom_qm_excl(i)
       IF ((ia.GT.mmdim%natm).OR.(ia.LT.1))THEN
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(6,200) i,'QM ATOM',ia,'IS NOT A VALID ATOM'
             ferror=.TRUE.
          ENDIF
       ENDIF
       IF (nat_cpmd(ia).GT.mmdim%natq) THEN
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(6,200) i,'   ATOM',ia,'IS NOT A QM ATOM'
             ferror=.TRUE.
          ENDIF
       ENDIF
       ia=excl_comm%atom_mm_excl(i)
       IF ((ia.GT.mmdim%natm).OR.(ia.LT.1))THEN
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(6,200) i,'MM ATOM',ia,'IS NOT A VALID ATOM'
             ferror=.TRUE.
          ENDIF
       ENDIF
       IF (nat_cpmd(ia).LE.mmdim%natq) THEN
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(6,200) i,'   ATOM',ia,'IS NOT AN MM ATOM'
             ferror=.TRUE.
          ENDIF
       ENDIF
    ENDDO
200 FORMAT (1x,'MM_INIT| EXCLUSION ',i5,':',2x,a,i8,1x,a)
    IF (ferror) CALL stopgm('MM_INIT','ERROR IN EXCLUSION LIST',& 
         __LINE__,__FILE__)
    IF (cprint%maxwriteatom.GT.mmdim%natm) cprint%maxwriteatom=mmdim%natm

    CALL mp_bcast(NAq,SIZE(naq),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(NAm,SIZE(nam),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast_byte(mmdim,size_in_bytes_of(mmdim),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(cprint%maxwriteatom,gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(cpat,SIZE(cpat),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(cpsp,SIZE(cpsp),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(nat_cpmd,SIZE(nat_cpmd),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(nat_grm,SIZE(nat_grm),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(gratom,SIZE(gratom),gparai%mmsource,parai%qmmmgrp)
    ! 
    ! mb      ALLOCATE(mm_charge(maxsys%nax,maxsys%nsx))
    ALLOCATE(mm_charge(maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(INC_l(maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (gparal%mmnode) THEN
       CALL  mm_pardef(rmass%pma,rmass%pma0,rmass%pmat0,rmass%pmatot,&
            mm_charge,clsaabox%box_au,INC_l,NCAG_l,&
            solvvv%nram_gr,solvvv%ncons_gr,&
            ions0%iatyp,mm_RAGGIO,&
            solsolv%nrpt,solsolv%nsolv,maxsys%nax,cpat,cpsp,mmdim%nspm,NAm)
    ENDIF
    CALL mp_bcast_byte(rmass,size_in_bytes_of(rmass),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(mm_charge,SIZE(mm_charge),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(clsaabox%box_au,SIZE(clsaabox%box_au),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(INC_l,SIZE(INC_l),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(NCAG_l,gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast_byte(solvvv,size_in_bytes_of(solvvv),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(mm_RAGGIO,SIZE(mm_RAGGIO),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(ions0%iatyp,SIZE(ions0%iatyp),gparai%mmsource,parai%qmmmgrp)
    ! 
    ! set initial positions; no initial velocities allowed
    ! 
    DEALLOCATE(tau0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tau0(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(velp(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (gparal%mmparent) THEN
       CALL mm_readgromos(maxsys%nax,maxsys%nsx,tau0)
       ! the quantum system can be split in two pieces due to pbc. 
       ! First take the minimal image with respect the position of the first atom

       IF (isos1%tcent)THEN
          ! Beware, changed (F Gervasio): now it loops on all atoms
          ! To look for the minimal image with the smaller rmax-rmin (r=x,y,z)
          CALL mm_best(tau0,maxsys%nax,maxsys%nsx,xm,velp,diffe1,diffe2,diffevol)
          CALL mm_min_im(tau0,velp,maxsys%nax,maxsys%nsx,xm)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,tau0(1,1,1),1)
          ! The cdm of the quantum system is translated to the center of the quantum cell
          CALL mm_center(tau0,xm,.FALSE.)
          IF (gparal%mmparent.AND.paral%io_parent)&
               WRITE(6,*) "Cell Volume",cell_com%volcel
          ! The minimal image with respect to the center of the quantum cell is taken
          xm(1)=0.5_real_8*cell_com%celldm(1)
          xm(2)=0.5_real_8*cell_com%celldm(2)*cell_com%celldm(1)
          xm(3)=0.5_real_8*cell_com%celldm(3)*cell_com%celldm(1)
          CALL mm_min_im(tau0,velp,maxsys%nax,maxsys%nsx,xm)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,tau0(1,1,1),1)
       ENDIF
    ENDIF

    CALL mp_bcast(tau0,SIZE(tau0),gparai%mmsource,parai%qmmmgrp)

    ! McB
    CALL zeroing(VELP)!,3*maxsys%nax*maxsys%nsx)  ! no initial velocity!!!

    IF ( gparal%mmparent ) THEN
       g96_vel%ntx_vel=0
       IF ( ntx.EQ.2 ) THEN
          g96_vel%ntx_vel=1
          ! ... load GROMOS velocities ...
          CALL mm_readgromos_vel(maxsys%nax,maxsys%nsx,velp)
       ENDIF
    ENDIF

    CALL mp_bcast_byte(g96_vel, size_in_bytes_of(g96_vel),gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(velp,SIZE(velp),gparai%mmsource,parai%qmmmgrp)
    ! McB
    DEALLOCATE(lvelini,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! TODO LVELINI has to be 0:, not sure how to implement this with memory90
    ALLOCATE(lvelini(0:maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO is=1,maxsys%nsx
       iis=(is-1)*(maxsys%nax+1)
       DO ia=1,maxsys%nax+1
          ! LVELINI(IIS+IA,1)=.FALSE.
       ENDDO
    ENDDO
    lvelini = .FALSE. ! use F90 array operations feature
    maxsys%nsx=NSX_s
    ! McB
    ! write table of the internal ordering.
    IF (paral%io_parent)&
         CALL fileopen(17,'QMMM_ORDER',fo_def,ferror)
    IF (paral%io_parent)&
         WRITE(17,'(''#'')')
    IF (paral%io_parent)&
         WRITE(17,'(''# NRTOT ='',i8,3x,''NATq ='',i8)') nrtot,mmdim%natq
    IF (paral%io_parent)&
         WRITE(17,'(''#'')')
    IF (paral%io_parent)&
         WRITE(17,'(''#__GROMOS_____CPMD_______is_______ia___SYMBOL_'')')
    i=0
    DO is=1,mmdim%nspm
       DO ia=1,NAm(is)
          i=i+1
          IF (paral%io_parent)&
               WRITE(17,'(4(1x,i8),7x,a2)') nat_grm(i),i,is,ia,elem%el(ions0%iatyp(is))
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(17,'(''#---------------------------------------------'')')
    IF (paral%io_parent)&
         CALL fileclose(17)
    ! McB
    RETURN
#else /* ! __GROMOS */
    ! directly call stub version for non-QM/MM compile.
    CALL mm_compat_init
    RETURN
#endif
  END SUBROUTINE mm_init
#if defined (__GROMOS)
  ! ==================================================================
  SUBROUTINE mm_best(tau0,nax,nsx,xm,velp,diffe1,diffe2,diffevol)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nax
    REAL(real_8)                             :: tau0(3,nax,*)
    INTEGER                                  :: nsx
    REAL(real_8)                             :: xm(3), velp(3,nax,*)
    INTEGER                                  :: diffe1, diffe2
    REAL(real_8)                             :: diffevol

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: diffesto

    diffesto=1.0e09_real_8
    diffe1=-HUGE(0)
    diffe2=-HUGE(0)
    DO is=1,mmdim%nspq
       DO ia=1,NAq(is)
          xm(1)=tau0(1,ia,is)
          xm(2)=tau0(2,ia,is)
          xm(3)=tau0(3,ia,is)
          CALL mm_min_im(tau0,velp,nax,nsx,xm)
          CALL dcopy(3*nax*nsx,velp(1,1,1),1,tau0(1,1,1),1)
          ! The cdm of the quantum system is translated to the center of the quantum cell
          CALL mm_center(tau0,xm,.FALSE.)
          ! Here it checks if the new minimal image has
          ! the smallest rmax-rmin
          diffevol=cell_com%volcel
          IF (diffesto.GT.diffevol.AND.diffevol.NE.0) THEN
             diffe1=ia
             diffe2=is
             diffesto=diffevol
             IF (gparal%mmparent.AND.paral%io_parent)&
                  WRITE(6,*) 'best',ia,is,diffevol
          ENDIF
       ENDDO
    ENDDO
    xm(1)=tau0(1,diffe1,diffe2)
    xm(2)=tau0(2,diffe1,diffe2)
    xm(3)=tau0(3,diffe1,diffe2)
    RETURN
  END SUBROUTINE mm_best
  ! ==================================================================
#endif

  ! this subroutine sets up some arrays and variable introduced with the
  ! QM/MM code to normal-QM compatible stub versions. the purpose of this
  ! is to limit the amount of special case code and thus improve manageability.
  ! AK 2005/05/25.
  SUBROUTINE mm_compat_init

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_compat_init'

    INTEGER                                  :: i, ia, ierr, is, NATm

    mm_stat=.TRUE.       ! we are (and stay) in QM dimensions.
    lqmmm%qmmm_verbose=.FALSE. ! no verbose QM/MM output.

    ! Initialize QM/MM energies
    eqm_r%eqm=0.0_real_8
    eqm_r%eqmmm=0._real_8
    eqm_r%eqmmm0=0._real_8
    eqm_r%eqmmm_cl=0._real_8
    eqm_r%emm=0._real_8
    eqm_r%eext0=0._real_8
    eqm_r%eexcl=0._real_8

    ! initialize NAT_cpmd to have an 1:1 atom index mapping 
    ! for normal inputs and outputs.
    mmdim%natm=0
    natm = 0
    !$omp parallel do private(IS) reduction(+:NATm)
    DO is=1,maxsys%nsx
       natm=natm+ions0%na(is)
    ENDDO
    mmdim%natm = natm
    mmdim%natq=mmdim%natm
    mmdim%naxq=maxsys%nax
    IF (cprint%maxwriteatom.GT.mmdim%natm) cprint%maxwriteatom=mmdim%natm

    CALL mp_bcast(mmdim%natm,parai%io_source,parai%cp_grp)
    CALL mp_bcast(mmdim%natq,parai%io_source,parai%cp_grp)
    CALL mp_bcast(mmdim%naxq,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cprint%maxwriteatom,parai%io_source,parai%cp_grp)

    i=(mmdim%natm+duat%ndat+2)
    ALLOCATE(nat_cpmd(0:i-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)   ! TODO nat_cpmd has to be 0:, not sure how to implement with memory90
    ! more: in subsequent do i=0, NATm loop we address to nat_cpmd(i), not (i-1)!
    ALLOCATE(nat_grm(0:i-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO nat_grm has to be 0:, not sure how to implement with memory90
    ALLOCATE(cpsp(i),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cpat(i),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gratom(maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(gratom)!,maxsys%nax*maxsys%nsx)
    !$omp parallel do private(I)
#ifdef __SR11000
    !poption parallel, tlocal(I)
#endif
    DO i=0,mmdim%natm
       nat_cpmd(i)=i
       nat_grm(i)=i
    ENDDO
    ! Dummy atoms are included in NAT_cpmd
    IF (duat%ndat.GT.0)THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')' DUMMY ATOMS INDICES: #DUM,  INDEX'
       ENDIF
       ia=0
       DO i=mmdim%natm+1,mmdim%natm+duat%ndat
          ia=ia+1
          nat_cpmd(i)=i
          nat_grm(i)=i
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(20X,I5,I8)') ia, i
          ENDIF
       ENDDO
    ENDIF
    i=1
    DO is=1,maxsys%nsx
       DO ia=1,ions0%na(is)
          cpat(i)=ia
          cpsp(i)=is
          gratom((is-1)*maxsys%nax+ia)=i
          i=i+1
       ENDDO
    ENDDO
    CALL mp_bcast(nat_cpmd,SIZE(nat_cpmd),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nat_grm,SIZE(nat_grm),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cpsp,SIZE(cpsp),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cpat,SIZE(cpat),parai%io_source,parai%cp_grp)
    CALL mp_bcast(gratom,SIZE(gratom),parai%io_source,parai%cp_grp)
    RETURN
  END SUBROUTINE mm_compat_init

END MODULE mm_init_utils
