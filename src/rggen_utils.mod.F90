#include "cpmd_global.h"

MODULE rggen_utils
  USE broy,                            ONLY: broy1
  USE cell,                            ONLY: cell_com,&
                                             lcell
  USE cppt,                            ONLY: gk,&
                                             gl,&
                                             hg,&
                                             igl,&
                                             indz,&
                                             inyh,&
                                             isptr,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE glopar_utils,                    ONLY: glopar
  USE gvec,                            ONLY: epsg,&
                                             gvec_com
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: latgen
  USE loadpa_utils,                    ONLY: leadim,&
                                             loadpa
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast, mp_sum,&
                                             mp_get_node_env,&
                                             mp_get_processor_name,&
                                             mp_max_processor_name
  USE mp_interface,                    ONLY: mp_max, mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE setsc_utils,                     ONLY: ihmat
  USE sphe,                            ONLY: gcutwmax
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: numcpus

#ifdef _HASNT_OMP_SET_NESTED
  !$ USE omp_lib, ONLY: omp_get_dynamic
#else
  !$ USE omp_lib, ONLY: omp_get_nested, omp_get_dynamic
#endif
  !$ USE omp_lib, ONLY: omp_get_max_active_levels

#ifndef _HASNT_OMP_45
  !$ USE omp_lib, ONLY: omp_get_proc_bind, omp_get_num_places, omp_get_place_proc_ids
  !$ USE omp_lib, ONLY: omp_get_place_num_procs
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rggen
  PUBLIC :: gvector
  PUBLIC :: recips

CONTAINS

  ! ==================================================================
  SUBROUTINE rggen
    ! ==--------------------------------------------------------------==
    ! == GENERATES THE RECIPROCAL LATTICE VECTORS WITH LENGTH SQUARED ==
    ! == LESS THAN GCUT, AND RETURNS THEM IN ORDER OF INCREASING      ==
    ! == LENGTH.                                                      ==
    ! == G=I*B1+J*B2+K*B3, WHERE B1,B2,B3 ARE THE VECTORS DEFINING    ==
    ! == THE RECIPROCAL LATTICE, AND I,J,K EACH GO FROM -(NR-1) TO    ==
    ! == +(NR-1).                                                     ==
    ! == N1,N2, AND N3 ARE THE FAST-FOURIER TRANSFORM INDICES OF G.   ==
    ! == IF I.GE.0, N1=I+1; IF I.LT.0, N1=I+1+NR, WITH SIMILAR        ==
    ! == DEFINITIONS FOR N2 AND N3. THIS CONNECTS G S WITH NEGATIVE   ==
    ! == COEFFICIENTS C TO THOSE WITH POSITIVE COEFFICIENTS BY        ==
    ! == INTRODUCING A FACTOR EXP(M*2PI*I)=1 IN THE FOURIER TRANSFORM.==
    ! == GCUT SHOULD BE .GE. GCUTW. NG IS THE TOTAL NUMBER OF VECTORS ==
    ! == WITH LENGTH SQUARED LESS THAN GCUT. NGW IS THE TOTAL NUMBER  ==
    ! == OF VECTORS WITH LENGTH SQUARED LESS THAN GCUTW.              ==
    ! == THE G S ARE IN UNITS OF 2PI/A.                               ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'rggen'
    REAL(real_8), PARAMETER                  :: delta = 1.e-6_real_8 

    INTEGER                                  :: i, ierr, ig, ig1, igs, indy1, &
                                                indy2, indy3, info, ish, &
                                                isub, j, k, nh1, nh2, nh3, &
                                                nshl
    LOGICAL                                  :: debug, ferror
    REAL(real_8)                             :: diff, g2, gcutbroy, glen, &
                                                hgold, t1, t2, t3

! ==--------------------------------------------------------------==

    debug=.FALSE.
    CALL tiset(procedureN,isub)
    ! Leading dimensions
    CALL leadim(parm%nr1,parm%nr2,parm%nr3,fpar%kr1,fpar%kr2,fpar%kr3)
    fpar%nnr1=fpar%kr1*fpar%kr2*fpar%kr3
    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    ! Global System Parameters
    CALL glopar
    ! Generate Taskgroups
    ! vw      CALL GROUPS
    ! Distribution of Work and Information to all Nodes
    CALL loadpa
    CALL aliasing
#if defined(__VECTOR)
    !$omp parallel do private(IG,I,J,K,T1,T2,T3)
#else
    !$omp parallel do private(IG,I,J,K,T1,T2,T3) schedule(static)
#endif
    DO ig=1,ncpw%nhg
       i=inyh(1,ig)-nh1
       j=inyh(2,ig)-nh2
       k=inyh(3,ig)-nh3
       t1=REAL(i,kind=real_8)*gvec_com%b1(1)+REAL(j,kind=real_8)*gvec_com%b2(1)+REAL(k,kind=real_8)*gvec_com%b3(1)
       t2=REAL(i,kind=real_8)*gvec_com%b1(2)+REAL(j,kind=real_8)*gvec_com%b2(2)+REAL(k,kind=real_8)*gvec_com%b3(2)
       t3=REAL(i,kind=real_8)*gvec_com%b1(3)+REAL(j,kind=real_8)*gvec_com%b2(3)+REAL(k,kind=real_8)*gvec_com%b3(3)
       hg(ig)=t1*t1+t2*t2+t3*t3
    ENDDO
    ! Broyden mixing cutoff and number of plane waves 
    broy1%ngbroy=0
    IF (broy1%ecutbroy.LE.0._real_8) THEN
       broy1%ecutbroy=gvec_com%gcut*parm%tpiba2
       gcutbroy=gvec_com%gcut
    ELSE
       gcutbroy=broy1%ecutbroy/parm%tpiba2
    ENDIF
    IF (broy1%tgbroy) THEN
       DO ig=1,ncpw%nhg
          g2=hg(ig)
          IF (g2.LT.gcutbroy) THEN
             broy1%ngbroy=broy1%ngbroy+1
          ENDIF
       ENDDO
    ENDIF
    IF (debug) THEN
       ! Test of reordering
       info=0
       DO ig=2,ncpw%nhg
          IF (hg(ig)-hg(ig-1).LT.-epsg) THEN
             IF (paral%io_parent) WRITE(6,*) 'PROC=',parai%me,' IG=',ig,&
                  ' HG(IG)=',hg(ig),' HG(IG-1)=',hg(ig-1)
             info=1
          ENDIF
       ENDDO
       IF (info.EQ.1) CALL stopgm(procedureN,'BAD REORDERING',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(nzh(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(indz(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ig=1,ncpw%nhg
       indy1=inyh(1,ig)
       indy2=inyh(2,ig)
       indy3=inyh(3,ig)
       nzh(ig) = indy1 + (indy2-1)*fpar%kr1s + (indy3-1)*fpar%kr1s*fpar%kr2s
       indy1=-indy1+nh1*2
       indy2=-indy2+nh2*2
       indy3=-indy3+nh3*2
       indz(ig) = indy1 + (indy2-1)*fpar%kr1s + (indy3-1)*fpar%kr1s*fpar%kr2s
    ENDDO
    ! Number of Shells
    ncpw%nhgl=1
    hgold=hg(1)
    DO ig=2,ncpw%nhg
       IF (ABS(hg(ig)-hgold).GT.delta) THEN
          ncpw%nhgl=ncpw%nhgl+1
          hgold=hg(ig)
          IF (hg(ig).LT.gcutwmax) ncpw%ngwl=ncpw%nhgl
       ENDIF
    ENDDO
    IF (ncpw%ngwl.EQ.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)')&
            'INCREASE THE FOURIER MESH OR THE ENERGY CUTOFF'
       CALL stopgm(procedureN,'MESH TOO SMALL',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(gl(ncpw%nhgl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nshl=1
    gl(1)=hg(1)
    DO ig=2,ncpw%nhg
       IF (ABS(hg(ig)-gl(nshl)).GT.delta) THEN
          nshl=nshl+1
          IF (nshl.GT.ncpw%nhgl) CALL stopgm(procedureN,&
               'NSHL BIGGER THAN NHGL',& 
               __LINE__,__FILE__)
          gl(nshl)=hg(ig)
       ENDIF
    ENDDO
    ALLOCATE(gk(3,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(igl(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(isptr(ncpw%nhgl+4),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Pointers PW --> Shell
    nshl=1
    DO ig=1,ncpw%nhg
100    CONTINUE
       IF (ABS(hg(ig)-gl(nshl)).GT.delta) THEN
          nshl=nshl+1
          GOTO 100
       ELSE
          igl(ig)=nshl
       ENDIF
    ENDDO
    IF (debug) THEN
       ! Test of shell reordering
       info=0
       DO ish=2,nshl
          IF (gl(ish)-gl(ish-1).LT.-epsg) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) 'PROC=',parai%me,' ISH=',ish,&
                  ' GL(ISH)=',gl(ish),' GL(ISH-1)=',gl(ish-1)
             info=1
          ENDIF
       ENDDO
       IF (info.EQ.1) CALL stopgm(procedureN,'BAD SHELL ORDERING',& 
            __LINE__,__FILE__)
    ENDIF
    ! Pointers to Shells
    DO ish=1,ncpw%nhgl
       glen=gl(ish)
       ig1=1
       IF (ish.GT.1) ig1=isptr(ish-1)
       DO ig=ig1,ncpw%nhg
          diff=hg(ig)-glen
          IF (ABS(diff).LE.delta .OR. diff.GT.delta) THEN
             isptr(ish)=ig
             GOTO 200
          ENDIF
       ENDDO
       isptr(ish)=ncpw%nhg+1
200    CONTINUE
    ENDDO
    isptr(ncpw%nhgl+1)=ncpw%nhg+1
    ! AK: for the sake of consistency. CELLDM() (the CELL keyword)
    ! AK: should _always_ be the actual cell, regardless of whether
    ! AK: we specify a reference cell or not. ...and particularly
    ! AK: it should not depend on whether we want to calculate
    ! AK: the stress tensor or not.
    ! AK: 2008/05/30.   
    ! IF(cntl%tprcp .OR. cntl%tpres) THEN
    IF (.NOT. lcell%tcellvectors) THEN
       CALL latgen(parm%ibrav,cell_com%celldm,parm%a1,parm%a2,parm%a3,parm%omega)
    ENDIF
    DO i=1,3
       metr_com%ht(1,i) = parm%a1(i)
       metr_com%ht(2,i) = parm%a2(i)
       metr_com%ht(3,i) = parm%a3(i)
    ENDDO
    CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
    ! AK:      ENDIF
    CALL gvector
    ! ==--------------------------------------------------------------==
    CALL numcpus(parai%ncpus)
    CALL print_omp_info( )

    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CALL prmem(procedureN)
    ! ==--------------------------------------------------------------==
    ! mb...write the file of the G-shell (to be used for S(q))
    IF (paral%io_parent .AND. cntl%gshellout) THEN
       CALL fileopen(71,"GSHELL",fo_def,ferror)
       WRITE(71,'(1x,2f15.6,i7)') 0.0_real_8,0.0_real_8,0! dummy origin
       WRITE(71,'(1x,2f15.6,i7)')&
            gl(1),SQRT(gl(1))*parm%tpiba,1
       DO igs=2,ncpw%nhgl
          IF (gl(igs).NE.gl(igs-1)) THEN
             WRITE(71,'(1x,2f15.6,i7)')&
                  GL(IGS),SQRT(GL(IGS))*parm%tpiba,IGS
          ENDIF
       ENDDO
       CALL fileclose(71)
    ENDIF
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE rggen

  SUBROUTINE print_omp_info()

    integer :: i, place_num, num_places, place_num_procs, partition_num_places
    integer, allocatable, dimension(:) :: place_proc_ids,partition_place_nums

    CHARACTER(mp_max_processor_name)         :: node_io_name, proc_name
    integer :: node_numtasks, node_taskid, iproc, from
    integer :: proc_bind
    LOGICAL                                  :: node_io_parent
    INTEGER, DIMENSION(:), ALLOCATABLE :: io_node_taskid_list

    IF (paral%io_parent) THEN
       WRITE(6,'(/," ",10("OPENMP"),"OPEN")')
       WRITE(6,'(A,T60,I6)') " OMP: NUMBER OF CPUS PER TASK",parai%ncpus
#ifndef _HASNT_OMP_SET_NESTED
       !$ WRITE(6,'(A,T60,L6)') " OMP: omp_get_nested",omp_get_nested( )
#endif
       !$ WRITE(6,'(A,T60,L6)') " OMP: omp_get_dynamic",omp_get_dynamic( )
       !$ WRITE(6,'(A,T54,I12)') " OMP: omp_get_max_active_levels",omp_get_max_active_levels( )
    ENDIF


    CALL mp_get_node_env ( parai%cp_grp, node_numtasks, node_taskid )
    !>vw get which node is the io node
    CALL mp_get_processor_name ( proc_name )
    IF( paral%io_parent ) node_io_name = proc_name
    CALL mp_bcast( node_io_name, parai%io_source, parai%cp_grp )
    node_io_parent = proc_name == node_io_name
    ALLOCATE(io_node_taskid_list(0:node_numtasks-1))
    io_node_taskid_list(:) = 0
    IF(node_io_parent) io_node_taskid_list(node_taskid) = parai%cp_me
    CALL mp_sum( io_node_taskid_list, SIZE(io_node_taskid_list), parai%cp_grp )

    DO iproc = 0,node_numtasks-1
       from = io_node_taskid_list( iproc )
       proc_bind = -1
#ifndef _HASNT_OMP_45
       !$ proc_bind = omp_get_proc_bind()
#endif
       num_places = -1
#ifndef _HASNT_OMP_45
       !$ num_places = omp_get_num_places()
#endif
       CALL mp_bcast( proc_bind, from, parai%cp_grp )
       CALL mp_bcast( num_places, from, parai%cp_grp )

       IF( paral%io_parent ) THEN
          IF( proc_bind > -1 .AND. num_places > -1 ) THEN
             !$ WRITE(6,'(1X,2(A,I0))') 'OMP: MPI local ', iproc,' global ', from
             !$ WRITE(6,'(1X,A,1X,I0)') 'OMP:   proc_bind',proc_bind
             !$ WRITE(6,'(1X,A,1X,I0)') 'OMP:   num_places', num_places
          ENDIF
       ENDIF

       DO place_num = 0, num_places-1
          place_num_procs = 0
#ifndef _HASNT_OMP_45
          !$ place_num_procs = omp_get_place_num_procs(place_num)
#endif
          allocate(place_proc_ids(place_num_procs))
#ifndef _HASNT_OMP_45
          !$ call omp_get_place_proc_ids(place_num,place_proc_ids)
#endif
          CALL mp_bcast( place_proc_ids, SIZE(place_proc_ids), from, parai%cp_grp )
          IF( paral%io_parent ) THEN
             !$ WRITE(*,'(1X,2(A,1X,I0))')  'OMP:     place_num',place_num,' place_num_procs', place_num_procs
             !$ WRITE(*,'(1X,A,99(1X,I0))') 'OMP:     place_proc_ids',place_proc_ids
          ENDIF
          deallocate(place_proc_ids)
       enddo
    enddo

    DEALLOCATE(io_node_taskid_list)

    IF (paral%io_parent) THEN
       WRITE(6,'(" ",10("OPENMP"),"OPEN",/)')
    ENDIF

  END SUBROUTINE print_omp_info
  ! ==================================================================
  SUBROUTINE gvector
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ig, iri1, iri2, iri3, ish, &
                                                nh1, nh2, nh3

! ==--------------------------------------------------------------==

    CALL recips(parm%alat,parm%a1,parm%a2,parm%a3,gvec_com%b1,gvec_com%b2,gvec_com%b3)
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    DO ig=1,ncpw%nhg
       iri1=inyh(1,ig)-nh1
       iri2=inyh(2,ig)-nh2
       iri3=inyh(3,ig)-nh3
       gk(1,ig)=iri1*gvec_com%b1(1)+iri2*gvec_com%b2(1)+iri3*gvec_com%b3(1)
       gk(2,ig)=iri1*gvec_com%b1(2)+iri2*gvec_com%b2(2)+iri3*gvec_com%b3(2)
       gk(3,ig)=iri1*gvec_com%b1(3)+iri2*gvec_com%b2(3)+iri3*gvec_com%b3(3)
       hg(ig)=gk(1,ig)**2+gk(2,ig)**2+gk(3,ig)**2
    ENDDO
    DO ish=1,ncpw%nhgl
       ig=isptr(ish)
       gl(ish)=hg(ig)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gvector
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE aliasing
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'aliasing'

    INTEGER                                  :: if1, if2, if3, iff, ig, iri1, &
                                                iri2, iri3, isub, nh1, nh2, &
                                                nh3
    REAL(real_8)                             :: xif1, xif2, xif3

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    if1=0
    if2=0
    if3=0
    DO ig=1,ncpw%nhg
       iri1=inyh(1,ig)
       iri2=inyh(2,ig)
       iri3=inyh(3,ig)
       IF (iri1.LT.1) if1=MAX(if1,ABS(iri1+1))
       IF (iri2.LT.1) if2=MAX(if2,ABS(iri2+1))
       IF (iri3.LT.1) if3=MAX(if3,ABS(iri3+1))
       IF (iri1.GT.spar%nr1s) if1=MAX(if1,iri1-spar%nr1s)
       IF (iri2.GT.spar%nr2s) if2=MAX(if2,iri2-spar%nr2s)
       IF (iri3.GT.spar%nr3s) if3=MAX(if3,iri3-spar%nr3s)
    ENDDO
    xif1=if1
    xif2=if2
    xif3=if3
    CALL mp_max(xif1,parai%allgrp)
    CALL mp_max(xif2,parai%allgrp)
    CALL mp_max(xif3,parai%allgrp)
    if1=xif1
    if2=xif2
    if3=xif3
    IF (paral%io_parent) THEN
       IF (if1.GT.0) WRITE(6,'(A,I3,A,I3,A)')&
            ' SET MESH PARAMETER NR1 TO ',spar%nr1s+2*if1,&
            ' (CURRENT ',spar%nr1s,')'
       IF (if2.GT.0) WRITE(6,'(A,I3,A,I3,A)')&
            ' SET MESH PARAMETER NR2 TO ',spar%nr2s+2*if2,&
            ' (CURRENT ',spar%nr2s,')'
       IF (if3.GT.0) WRITE(6,'(A,I3,A,I3,A)')&
            ' SET MESH PARAMETER NR3 TO ',spar%nr3s+2*if3,&
            ' (CURRENT ',spar%nr3s,')'
       iff=if1+if2+if3
       IF (iff.GT.0) THEN
          WRITE(6,'(A,3I4)') ' CURRENT MESH: ',spar%nr1s,spar%nr2s,spar%nr3s
          WRITE(6,'(A)')&
               ' USE OPTION MESH IN SYSTEM SECTION'
          CALL stopgm(procedureN,'MESH TOO SMALL',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE aliasing
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE recips(a,a1,a2,a3,b1,b2,b3)
    ! ==--------------------------------------------------------------==
    ! == GENERATES THE RECIPROCAL LATTICE VECTORS B1,B2,B3 GIVEN THE  ==
    ! == REAL SPACE VECTORS A1,A2,A3.  B IN UNITS OF 2PI/A.           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a, a1(3), a2(3), a3(3), &
                                                b1(3), b2(3), b3(3)

    INTEGER                                  :: i, iperm, ir, j, k, l
    REAL(real_8)                             :: den, s

    den=0._real_8
    i=1
    j=2
    k=3
    s=1._real_8
1   CONTINUE
    DO iperm=1,3
       den=den+s*a1(i)*a2(j)*a3(k)
       l=i
       i=j
       j=k
       k=l
    ENDDO
    i=2
    j=1
    k=3
    s=-s
    IF (s.LT.0._real_8) go to 1
    i=1
    j=2
    k=3
    den=a/ABS(den)
    DO ir=1,3
       b1(ir)=den*(a2(j)*a3(k)-a2(k)*a3(j))
       b2(ir)=den*(a3(j)*a1(k)-a3(k)*a1(j))
       b3(ir)=den*(a1(j)*a2(k)-a1(k)*a2(j))
       l=i
       i=j
       j=k
       k=l
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE recips
  ! ==================================================================

END MODULE rggen_utils
