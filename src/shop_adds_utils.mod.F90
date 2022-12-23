MODULE shop_adds_utils
  USE coor,                            ONLY: fion,&
                                             velp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_cputime,&
                                             m_walltime
  USE mm_dimmod,                       ONLY: mmdim,&
                                             naq
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE shop,                            ONLY: fs0,&
                                             fs1,&
                                             sh02,&
                                             tlsd0
  USE shop_ekinqm,                     ONLY: ekinqmsh
  USE shop_rest,                       ONLY: prob,&
                                             prob1,&
                                             sh03,&
                                             zz
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             maxsys,&
                                             ncpw,&
                                             norbpe,&
                                             parap
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_shmd
  PUBLIC :: state_select
  PUBLIC :: s0_s1_overlap
  PUBLIC :: decide

  ! Period parameter
  INTEGER, PARAMETER, PRIVATE :: n = 624, n1 = n + 1
  TYPE rng_block_type
     INTEGER                                  :: mt(0:n-1)
     INTEGER                                  :: mti=n1
  END TYPE rng_block_type
  TYPE (rng_block_type), SAVE, PRIVATE              :: rng_block

CONTAINS

  ! ==================================================================
  ! == ADDITIONAL ROUTINES FOR SURFACE HOPPING                      ==
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE write_shmd(iwahl,iunit,infi,tempp)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iwahl, iunit, infi
    REAL(real_8)                             :: tempp

    IF (iwahl.EQ.0) THEN
       IF (paral%io_parent)&
            WRITE(iunit,'(5A)') '#___STEP______T[K]',&
            '_________E(1)[a.u.]_________E(2)[a.u.]',&
            '___ISH_________DE____________PROB',&
            '____________RAND_________________NORM',&
            '_'

    ELSEIF (iwahl.EQ.1) THEN
       IF (paral%io_parent)&
            WRITE(iunit,'(1X,I7,F10.3,2(1X,F18.10),1X,I5,1X,F10.5, 2(1X,D15.7),1X,F20.15)')&
            INFI,TEMPP,sh03%ec(1),sh03%ec(2),sh03%ioldsurf,sh03%ec(2)-sh03%ec(1),&
            PROB(sh03%ioldsurf),ZZ,prob1%d1sq+prob1%d2sq

    ELSEIF (iwahl.EQ.-1) THEN
       IF (paral%io_parent)&
            WRITE(iunit,'(5A)') '#-----------------',&
            '--------------------------------------',&
            '---------------------------------',&
            '-------------------------------------',&
            '-'

    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,I5)') ' INVALID OPTION:',iwahl
       CALL stopgm('WRITE_SHMD',' ',& 
            __LINE__,__FILE__)

    ENDIF

    RETURN
  END SUBROUTINE write_shmd

  ! ==================================================================
  SUBROUTINE state_select(tag)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: tag

    INTEGER                                  :: i
    REAL(real_8)                             :: xsaim, xsnow, xstates

    IF (tag.EQ."S0" ) THEN
       IF ( paral%parent ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,'' =='',5X,A16,39X,''==''/)') 'STATE SELECT: S0'
       ENDIF
       lspin2%tlse=.FALSE.
       IF ( tlsd0 ) THEN
          cntl%tlsd=.TRUE.
          clsd%nlsd=2
          clsd%nlsx=3
          crge%n=sh02%nst_s0
       ELSE
          clsd%nlsd=1
          clsd%nlsx=1
          crge%n=sh02%nst_s0
       ENDIF
       ! -PARALLEL
       xstates=REAL(crge%n,kind=real_8)
       xsnow=0.0_real_8
       DO i=parai%nproc,1,-1
          xsaim = xsnow + xstates/parai%nproc
          parap%nst12(i-1,1)=NINT(xsnow)+1
          parap%nst12(i-1,2)=NINT(xsaim)
          IF ( NINT(xsaim).GT.crge%n ) parap%nst12(i-1,2)=crge%n
          IF ( i.EQ.1 ) parap%nst12(i-1,2)=crge%n
          xsnow = xsaim
       ENDDO
       norbpe=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
       ! WRITE(6,*)'S0',me,NORBPE,MEPOS,NST12(MEPOS,2),NST12(MEPOS,1)
       ! -ENDPARALLEL
       CALL dcopy(sh02%nst_s0,fs0,1,crge%f,1)
    ELSEIF (tag.EQ."S1") THEN
       ! ...    SELECT S1 
       IF ( paral%parent ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,'' =='',5X,A16,39X,''==''/)') 'STATE SELECT: S1'
       ENDIF
       cntl%tlsd=.FALSE.
       lspin2%tlse=.TRUE.
       lspin2%troks=.TRUE.
       clsd%nlsd=4
       clsd%nlsx=3
       crge%n=sh02%nst_s1
       ! -PARALLEL---------------------------------
       xstates=REAL(crge%n,kind=real_8)
       xsnow=0.0_real_8
       DO i=parai%nproc,1,-1
          xsaim = xsnow + xstates/parai%nproc
          parap%nst12(i-1,1)=NINT(xsnow)+1
          parap%nst12(i-1,2)=NINT(xsaim)
          IF (NINT(xsaim).GT.crge%n) parap%nst12(i-1,2)=crge%n
          IF (i.EQ.1) parap%nst12(i-1,2)=crge%n
          xsnow = xsaim
       ENDDO
       norbpe=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
       ! WRITE(6,*)'S1',me,NORBPE,MEPOS,NST12(MEPOS,2),NST12(MEPOS,1)
       ! -ENDPARALLEL   
       ! McB   copy NST_S0-NST_S1 ZEROS, TOO !
       CALL dcopy(sh02%nst_s0,fs1,1,crge%f,1)
       ! CALL DCOPY(NST_S1,FS1,1,F,1)
       ! McB
    ELSE
       CALL stopgm('STATE_SELECT','UNKNOWN TAG',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE state_select
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE s0_s1_overlap(c0,cm,det,c12,c21)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), cm(ncpw%ngw,*)
    REAL(real_8)                             :: det, c12, c21

    CHARACTER(*), PARAMETER                  :: procedureN = 's0_s1_overlap'

    INTEGER                                  :: idstate, ierr, ig, iistate, &
                                                info, istate, iustate, &
                                                jdstate, jjstate, jstate, &
                                                justate
    REAL(real_8)                             :: det2, detdm, detum
    REAL(real_8), ALLOCATABLE :: d2dotm(:,:), d2dotm21(:,:), ddotm(:,:), &
      ddotm21(:,:), dm(:,:), dm2(:,:), u2dotm(:,:), u2dotm21(:,:), &
      udotm(:,:), udotm21(:,:), um(:,:), um2(:,:), work(:,:)

    ALLOCATE(um(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dm(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(um2(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dm2(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(udotm(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddotm(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(u2dotm(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(d2dotm(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(udotm21(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddotm21(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(u2dotm21(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(d2dotm21(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(sh02%nelb2,sh02%nelb2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (tlsd0)THEN

       ! if ( parent ) WRITE(6,*)'IN OVERLAP',NELB2

       CALL zeroing(um)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(dm)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(um2)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(dm2)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(udotm)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(ddotm)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(u2dotm)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(d2dotm)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(udotm21)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(ddotm21)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(u2dotm21)!,sh02%nelb2*sh02%nelb2)
       CALL zeroing(d2dotm21)!,sh02%nelb2*sh02%nelb2)

       DO istate=1,sh02%nelb2
          iistate=sh02%nelb2+istate
          iustate=sh02%nst_s0+istate
          idstate=iustate
          IF (istate.EQ.sh02%nelb2)idstate=iustate+1

          DO jstate=1,sh02%nelb2
             jjstate=sh02%nelb2+jstate
             justate=jstate+sh02%nst_s0
             jdstate=justate
             IF (jstate.EQ.sh02%nelb2)jdstate=justate+1
             ! 
             DO ig=1,ncpw%ngw
                IF (ig.EQ.1.AND.geq0) THEN
                   ! ...<S0|m1>
                   um(istate,jstate)=um(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JUSTATE)

                   dm(istate,jstate)=dm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        C0(IG,JDSTATE)
                   ! ...<S0|d/dt|m1>
                   udotm(istate,jstate)=udotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JUSTATE)

                   ddotm(istate,jstate)=ddotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JDSTATE)
                   ! ...<m1|d/dt|S0>
                   udotm21(istate,jstate)=udotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JSTATE)

                   ddotm21(istate,jstate)=ddotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JJSTATE)
                   ! ...<S0|m2>
                   um2(istate,jstate)=um2(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JDSTATE)

                   dm2(istate,jstate)=dm2(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        C0(IG,JUSTATE)
                   ! ...<S0|d/dt|m2>
                   u2dotm(istate,jstate)=u2dotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JDSTATE)

                   d2dotm(istate,jstate)=d2dotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JUSTATE)
                   ! ...<m2|d/dt|S0>
                   u2dotm21(istate,jstate)=u2dotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JJSTATE)

                   d2dotm21(istate,jstate)=d2dotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JJSTATE)
                ELSE
                   ! ...<S0|m1>
                   um(istate,jstate)=um(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JUSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(C0(IG,JUSTATE))

                   dm(istate,jstate)=dm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        C0(IG,JDSTATE)+&
                        C0(IG,IISTATE)*&
                        CONJG(C0(IG,JDSTATE))
                   ! ...<S0|d/dt|m1>
                   udotm(istate,jstate)=udotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JUSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(CM(IG,JUSTATE))

                   ddotm(istate,jstate)=ddotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JDSTATE)+&
                        C0(IG,IISTATE)*&
                        CONJG(CM(IG,JDSTATE))
                   ! ...<m1|d/dt|S0>
                   udotm21(istate,jstate)=udotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JSTATE)+&
                        C0(IG,IUSTATE)*&
                        CONJG(CM(IG,JSTATE))

                   ddotm21(istate,jstate)=ddotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JJSTATE)+&
                        C0(IG,IDSTATE)*&
                        CONJG(CM(IG,JJSTATE))
                   ! ...<S0|m2>
                   um2(istate,jstate)=um2(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JDSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(C0(IG,JDSTATE))


                   dm2(istate,jstate)=dm2(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        C0(IG,JUSTATE)+&
                        C0(IG,IISTATE)*&
                        CONJG(C0(IG,JUSTATE))
                   ! ...<S0|d/dt|m2>
                   u2dotm(istate,jstate)=u2dotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JDSTATE)+&
                        C0(IG,IISTATE)*&
                        CONJG(CM(IG,JDSTATE))

                   d2dotm(istate,jstate)=d2dotm(istate,jstate) +&
                        CONJG(C0(IG,IISTATE))*&
                        CM(IG,JUSTATE)+&
                        C0(IG,IISTATE)*&
                        CONJG(CM(IG,JUSTATE))
                   ! ...<m2|d/dt|S0>
                   u2dotm21(istate,jstate)=u2dotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JJSTATE)+&
                        C0(IG,IDSTATE)*&
                        CONJG(CM(IG,JJSTATE))

                   d2dotm21(istate,jstate)=d2dotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JJSTATE)+&
                        C0(IG,IUSTATE)*&
                        CONJG(CM(IG,JJSTATE))

                ENDIF
             ENDDO
             ! if ( parent ) then 
             ! WRITE(6,*)'UM',istate,jstate,udotm(istate,jstate)
             !                  ,u2dotm(istate,jstate)
             ! WRITE(6,*)'DM',istate,jstate,ddotm(istate,jstate)
             !                  ,d2dotm(istate,jstate)
             ! endif
          ENDDO
       ENDDO

       CALL mp_sum(um,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(dm,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(um2,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(dm2,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(udotm,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(ddotm,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(u2dotm,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(d2dotm,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(udotm21,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(ddotm21,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(u2dotm21,sh02%nelb2*sh02%nelb2,parai%allgrp)
       CALL mp_sum(d2dotm21,sh02%nelb2*sh02%nelb2,parai%allgrp)

       ! NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
       ! ...<S0|m1>
       CALL det_and_inv(sh02%nelb2,um,work,info,detum)
       CALL det_and_inv(sh02%nelb2,dm,work,info,detdm)
       det=detum*detdm/SQRT(2.0_real_8)
       ! WRITE(6,*)'DET 1',det
       ! .................................................................
       c12=0.0_real_8
       c21=0.0_real_8

       DO istate=1,sh02%nelb2
          DO jstate=1,sh02%nelb2
             c12=c12+det*( udotm(istate,jstate)*um(jstate,istate)+&
                  DDOTM(ISTATE,JSTATE)*DM(JSTATE,ISTATE) )
             c21=c21+det*( udotm21(istate,jstate)*um(istate,jstate)+&
                  DDOTM21(ISTATE,JSTATE)*DM(ISTATE,JSTATE) )
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)'D_21, D_12 [1]',c21,c12
       ! ...<S0|m2>
       CALL det_and_inv(sh02%nelb2,um2,work,info,detum)
       CALL det_and_inv(sh02%nelb2,dm2,work,info,detdm)
       det2=detum*detdm/SQRT(2.0_real_8)
       ! WRITE(6,*)'DET 2',det2
       ! .................................................................

       DO istate=1,sh02%nelb2
          DO jstate=1,sh02%nelb2
             c12=c12+det2*( u2dotm(istate,jstate)*um2(jstate,istate)+&
                  D2DOTM(ISTATE,JSTATE)*DM2(JSTATE,ISTATE) )
             c21=c21+det2*(&
                  U2DOTM21(ISTATE,JSTATE)*UM2(ISTATE,JSTATE)+&
                  D2DOTM21(ISTATE,JSTATE)*DM2(ISTATE,JSTATE) )
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)'D_21, D_12 [2]',c21,c12
       det=det+det2
       ! NNNNNN UM AND DM ARE NOW THE INVERSE MATRICES !!!!!!!!!!! 
       ! STOP
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*)'ACHTUNG, KEIN LSD! !'

       CALL zeroing(um)!,sh02%nst_s0*sh02%nst_s0)
       CALL zeroing(dm)!,sh02%nst_s0*sh02%nst_s0)
       CALL zeroing(udotm)!,sh02%nst_s0*sh02%nst_s0)
       CALL zeroing(ddotm)!,sh02%nst_s0*sh02%nst_s0)
       CALL zeroing(udotm21)!,sh02%nst_s0*sh02%nst_s0)
       CALL zeroing(ddotm21)!,sh02%nst_s0*sh02%nst_s0)

       DO istate=1,sh02%nst_s0
          iustate=sh02%nst_s0+istate
          idstate=iustate
          IF ( istate.EQ.sh02%nelb2 ) idstate=iustate+1

          DO jstate=1,sh02%nst_s0
             justate=jstate+sh02%nst_s0
             jdstate=justate
             IF ( jstate.EQ.sh02%nelb2 ) jdstate=justate+1

             DO ig=1,ncpw%ngw
                IF ( (ig.EQ.1).AND.geq0 ) THEN

                   um(istate,jstate)=um(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JUSTATE)

                   dm(istate,jstate)=dm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JDSTATE)

                   udotm(istate,jstate)=udotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JUSTATE)

                   ddotm(istate,jstate)=ddotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JDSTATE)

                   udotm21(istate,jstate)=udotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JSTATE)

                   ddotm21(istate,jstate)=ddotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JSTATE)
                ELSE

                   um(istate,jstate)=um(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JUSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(C0(IG,JUSTATE))

                   dm(istate,jstate)=dm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        C0(IG,JDSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(C0(IG,JDSTATE))

                   udotm(istate,jstate)=udotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JUSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(CM(IG,JUSTATE))

                   ddotm(istate,jstate)=ddotm(istate,jstate) +&
                        CONJG(C0(IG,ISTATE))*&
                        CM(IG,JDSTATE)+&
                        C0(IG,ISTATE)*&
                        CONJG(CM(IG,JDSTATE))

                   udotm21(istate,jstate)=udotm21(istate,jstate) +&
                        CONJG(C0(IG,IUSTATE))*&
                        CM(IG,JSTATE)+&
                        C0(IG,IUSTATE)*&
                        CONJG(CM(IG,JSTATE))

                   ddotm21(istate,jstate)=ddotm21(istate,jstate) +&
                        CONJG(C0(IG,IDSTATE))*&
                        CM(IG,JSTATE)  +&
                        C0(IG,IDSTATE)*&
                        CONJG(CM(IG,JSTATE))

                ENDIF
             ENDDO

          ENDDO
       ENDDO

       CALL mp_sum(um,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       CALL mp_sum(dm,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       CALL mp_sum(udotm,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       CALL mp_sum(ddotm,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       CALL mp_sum(udotm21,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       CALL mp_sum(ddotm21,sh02%nst_s0*sh02%nst_s0,parai%allgrp)
       ! NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
       CALL det_and_inv(sh02%nst_s0,um,work,info,detum)
       CALL det_and_inv(sh02%nst_s0,dm,work,info,detdm)
       ! NNNNNN UM AND DM ARE NOW THE INVERSE MATRICES !!!!!!!!!!! 
       det=SQRT(2.0_real_8)*detum*detdm
       ! .................................................................
       c12=0.0_real_8
       c21=0.0_real_8
       DO istate=1,sh02%nelb2
          DO jstate=1,sh02%nelb2
             c12=c12+det*( udotm(istate,jstate)*um(jstate,istate)+&
                  DDOTM(ISTATE,JSTATE)*DM(JSTATE,ISTATE) )
             c21=c21+det*( udotm21(istate,jstate)*um(istate,jstate)+&
                  DDOTM21(ISTATE,JSTATE)*DM(ISTATE,JSTATE) )
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)'D_21, D_12',c21,c12
    ENDIF

    ! ...FREE MEMORY
    DEALLOCATE(um,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(um2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dm2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(udotm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddotm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(u2dotm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(d2dotm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(udotm21,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddotm21,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(u2dotm21,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(d2dotm21,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s0_s1_overlap
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE det_and_inv(n,a,b,info,det)
    ! ==--------------------------------------------------------------==
    ! == Inverse Matrix A(N,N)                                        ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:  N Dimension                                          ==
    ! ==         A(N,N) Matrix                                        ==
    ! == OUTPUT: A(N,N) inverse matrix of A                           ==
    ! ==         B is a work array                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n), b(n,n)
    INTEGER                                  :: info
    REAL(real_8)                             :: det

    INTEGER                                  :: i, lb, ipivot(n)
    REAL(real_8)                             :: sign, work(n)

    IF (n.EQ.1) THEN
       a(1,1)=1._real_8/a(1,1)
    ELSE
       lb=(n-1)*n
       ! Compute an LU factorization
       ! WRITE(6,*)'A',A(1,1),A(1,2)
       ! WRITE(6,*)'A',A(2,1),A(2,2)
       CALL dgetrf(n,n,a,n,ipivot,info)
       ! WRITE(6,*)'INFO',INFO
       det=1.0_real_8
       DO i=1,n
          sign=1.0_real_8
          IF (ipivot(i).NE.i)sign=-1.0_real_8
          ! WRITE(6,*)'A(I,I)',A(I,I),IPIVOT(I)
          det=det*sign*a(i,i)
          ! WRITE(6,*)'det',DET,SIGN
       ENDDO
       IF (info.EQ.0) THEN
          ! Compute the inverse of a matrix using the LU factorization
          CALL dgetri(n,a,n,ipivot,work,n,info)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE det_and_inv
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE decide(e,fion0,fion1,jump)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: e(2), fion0(:,:,:), &
                                                fion1(:,:,:)
    LOGICAL                                  :: jump

    CHARACTER(len=2)                         :: tag(2)
    INTEGER                                  :: i, ia, is, j, jsurf, k
    REAL(real_8)                             :: arg, delt, ekinp, ekinsh, &
                                                ppp(2), skal

    COMMON /ekin/ekinp

    delt=cntr%delt_ions
    ! DELT=DELT_ELEC
    ! McB
    tag(1)='S0'
    tag(2)='S1'
    IF (sh03%isurf.EQ.1) jsurf=2
    IF (sh03%isurf.EQ.2) jsurf=1
    sh03%ioldsurf=sh03%isurf
    IF (.NOT.jump) THEN
       CALL state_select(tag(sh03%isurf))
       ener_com%etot=e(sh03%isurf)
       IF (sh03%isurf.EQ.2) THEN
          DO i=1,3
             DO j=1,maxsys%nax
                DO k=1,maxsys%nsx
                   fion(i,j,k)= fion1(i,j,k)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO i=1,3
             DO j=1,maxsys%nax
                DO k=1,maxsys%nsx
                   fion(i,j,k)= fion0(i,j,k)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ELSE
       ! NNNN IF SURFACE JUMP ALLOWED:
       ppp(2)=-delt*prob1%d22dot/prob1%d2sq
       ppp(1)=-delt*prob1%d11dot/prob1%d1sq
       prob(1)=ppp(1)
       prob(2)=ppp(2)
       ! ---- WRITE COEFFICIENTS TO FILE --------------------------
       ! if (parent) then
       ! call fileopen(789,'prob.dat',FO_APP,FERROR)
       ! write(789,12)infi,d1sq,d2sq
       ! 12     format(i6,4(1x,f12.8))     
       ! call fileclose(789)
       ! endif
       ! --------------------------------------------------------      
       CALL shop_grnd(zz)
       ! McB
       ! ZZ=ZZ*1.0e-5_real_8
       ! McB
       IF (paral%io_parent)&
            WRITE(6,*)'tot',prob1%d1sq+prob1%d2sq
       ! IF(PARENT) write(6,*)'d1,d2',d1sq,d2sq
       ! IF(PARENT) write(6,*)'d11dot,d22dot',d11dot,d22dot
       ! IF(PARENT) WRITE(6,*)'ZZ',ZZ,'PROB',P(ISURF),ISURF
       IF (paral%io_parent)&
            WRITE(6,*)'ZZ',zz,'PROB',ppp(sh03%isurf),sh03%isurf
       ! NNNN IF RANDOM NUMBER SMALLER THAN TRANSITION PROBABILITY:
       ! IF (ZZ.LT.P(ISURF)) THEN
       IF (zz.LT.ppp(sh03%isurf)) THEN
          IF (cntl%tqmmm) THEN
             ekinsh=ekinqmsh
          ELSE
             ekinsh=ekinp
          ENDIF
          arg=(e(sh03%isurf)-e(jsurf))/ekinsh+1.0_real_8
          IF (arg.GT.0.0_real_8) THEN! classically allowed transition
             CALL state_select(tag(jsurf))
             ener_com%etot=e(jsurf)
             skal=SQRT(arg)
             IF (paral%io_parent)&
                  WRITE(6,*)'SKAL',skal,ekinsh
             DO i=1,3
                ! scale QM centres only!
                ! DO J=1,maxsys%nax
                ! DO K=1,maxsys%nsx
                DO is=1,mmdim%nspq
                   DO ia=1,NAq(is)
                      velp(i,ia,is)= skal*velp(i,ia,is)! rescale velos
                   ENDDO
                ENDDO
             ENDDO
             IF (jsurf.EQ.2)THEN
                DO i=1,3
                   DO j=1,maxsys%nax
                      DO k=1,maxsys%nsx
                         fion(i,j,k)= fion1(i,j,k)
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO i=1,3
                   DO j=1,maxsys%nax
                      DO k=1,maxsys%nsx
                         fion(i,j,k)= fion0(i,j,k)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             sh03%isurf=jsurf
          ELSE! classically forbidden transition (arg<0)
             ! NNNN IF RANDOM NUMBER NOT SMALLER THAN TRANSITION PROBABILITY:
             CALL state_select(tag(sh03%isurf))
             ener_com%etot=e(sh03%isurf)
             IF (sh03%isurf.EQ.2)THEN! isurf=2
                DO i=1,3
                   DO j=1,maxsys%nax
                      DO k=1,maxsys%nsx
                         fion(i,j,k)= fion1(i,j,k)
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO i=1,3
                   DO j=1,maxsys%nax
                      DO k=1,maxsys%nsx
                         fion(i,j,k)= fion0(i,j,k)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF        ! isurf=2
          ENDIF             ! arg<0
       ENDIF                 ! zz>p
    ENDIF                 ! jump=true
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE decide
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE shop_grnd(rnv)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rnv

    INTEGER, PARAMETER :: lmask = 2147483647 , m = 397, mata = -1727483681 , &
      tmaskb = -1658038656, tmaskc = -272236544, umask = -2147483647

    CHARACTER(len=15)                        :: rng_state_file
    INTEGER                                  :: iunit, kk, y
    INTEGER, DIMENSION(0:1), SAVE            :: mag01 = (/0, mata/)
    LOGICAL                                  :: ferror

! -----------------------------------------------------------------------
! generate N words at one time

    IF (rng_block%mti.GE.n) THEN

       ! initialize or restore random number generator on first call.
       IF (rng_block%mti.EQ.n+1) THEN
          CALL shop_sgrnd
       ENDIF
       ! 
       DO kk=0,n-m-1
          y=IOR(IAND(rng_block%mt(kk),umask),IAND(rng_block%mt(kk+1),lmask))
          rng_block%mt(kk)=IEOR(IEOR(rng_block%mt(kk+m),ISHFT(y,-1)),mag01(IAND(y,1)))
       ENDDO

       DO kk=n-m,n-2
          y=IOR(IAND(rng_block%mt(kk),umask),IAND(rng_block%mt(kk+1),lmask))
          rng_block%mt(kk)=IEOR(IEOR(rng_block%mt(kk+(m-n)),ISHFT(y,-1)),mag01(IAND(y,1)))
       ENDDO

       y=IOR(IAND(rng_block%mt(n-1),umask),IAND(rng_block%mt(0),lmask))
       rng_block%mt(n-1)=IEOR(IEOR(rng_block%mt(m-1),ISHFT(y,-1)),mag01(IAND(y,1)))
       rng_block%mti = 0
    ENDIF
    ! 
    y=rng_block%mt(rng_block%mti)
    rng_block%mti=rng_block%mti+1
    ! 
    y=IEOR(y,ISHFT(y,-11))
    y=IEOR(y,IAND(ISHFT(y,7),tmaskb))
    y=IEOR(y,IAND(ISHFT(y,15),tmaskc))
    y=IEOR(y,ISHFT(y,-18))
    ! 
    IF (y.LT.0) THEN
       rnv=(REAL(y,kind=real_8)+2.0_real_8**32)/(2.0_real_8**32-1.0_real_8)
    ELSE
       rnv=REAL(y,kind=real_8)/(2.0_real_8**32-1.0_real_8)
    ENDIF

    IF (paral%parent) THEN
       ! record the state of the random number generator in the file
       ! SHOP_RNG_STATE. this should not pose a performance problem,
       ! since this is only done once per cntl%md step.
       ferror=.FALSE.
       iunit=69
       rng_state_file='SHOP_RNG_STATE'
       IF (paral%io_parent)&
            CALL fileopen(iunit,rng_state_file,fo_def,ferror)
       IF (.NOT.ferror) THEN
          IF (paral%io_parent)&
               WRITE(iunit,'(I14)')   rng_block%mti
          IF (paral%io_parent)&
               WRITE(iunit,'(5I14)') (rng_block%mt(kk),kk=0,n-1)
          IF (paral%io_parent)&
               CALL fileclose(iunit)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE shop_grnd
  ! ==================================================================

  ! ==================================================================
  ! seed or restore the random number generator.
  ! for surface hopping the seed is supposed to be new on every
  ! new run, so we try to read the recorded state of the last run first.
  SUBROUTINE shop_sgrnd
    ! ==--------------------------------------------------------------==

    CHARACTER(len=15)                        :: rng_state_file
    INTEGER                                  :: i, iunit, mtitmp, mttmp, seed
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,200)
       ferror=.FALSE.
       iunit=69
       rng_state_file='SHOP_RNG_STATE'
       IF (paral%io_parent)&
            CALL fileopen(iunit,rng_state_file,fo_old,ferror)
       IF (.NOT.ferror) THEN
          ! read state vector from restart in a fault tolerant way.
          IF (paral%io_parent)&
               WRITE(6,'(1X,2A)',advance="no")'READING OLD STATE FROM: ',&
               rng_state_file
          IF (paral%io_parent)&
               READ(iunit,'(I14)',err=100,END=100)   rng_block%mti
          IF (paral%io_parent)&
               READ(iunit,'(5I14)',err=100,END=100) (rng_block%mt(i),i=0,n-1)
          ferror=.FALSE.
          IF (paral%io_parent)&
               WRITE(6,'('' OK.'',/)')
          GOTO 110
100       CONTINUE
          IF (paral%io_parent)&
               WRITE(6,'('' ERROR.'',/)')
          ferror=.TRUE.
110       CONTINUE
          IF (paral%io_parent)&
               CALL fileclose(iunit)
       ENDIF

       ! no previous random number state was available or readable. 
       ! setting initial seeds for the rng_block%mt() array using 
       ! a simple random number generator.
       IF (ferror) THEN
          ! SEED=ANINT(abs(TIMEC()+TIMEF()+repprngu())*1.0e8_real_8)+42
          ! McB
          seed=ANINT(ABS(m_cputime()+m_walltime()+repprngu()*1.0e8_real_8))+42
          ! SEED=ANINT(abs(TIMEC()+TIMEF()+repprngu())*1.0e2_real_8)+42
          ! McB

          ! McB
          ! if (parent) then
          ! write (6,*) 'sgrnd> ',' TIMEC() = ',TIMEC()
          ! write (6,*) 'sgrnd> ',' TIMEF() = ',TIMEF()
          ! write (6,*) 'sgrnd> ',' repprngu() = ',repprngu()
          ! endif
          ! McB
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I12,/)') 'SEEDING RNG WITH SEED:',seed
          ENDIF
          rng_block%mt(0)= IAND(seed,-1)
          DO mtitmp =1, n-1
             ! McB
             mttmp=69069 * rng_block%mt(mtitmp-1)
             ! mt(mti) = IAND(69069 * mt(mti-1),-1)
             rng_block%mt(mtitmp) = IAND(mttmp,-1)

             ! if (parent) write (6,*) 'sgrnd> ',mti,mttmp,mt(mti)
             ! McB
          ENDDO
          rng_block%mti=mtitmp
          ! McB
          ! stop 
          ! McB
       ELSE
       ENDIF
    ENDIF

    ! distribute the RNG state vector across all nodes
    ! so we get the same sequence of numbers on all of them.
    CALL mp_bcast(rng_block%mti,parai%source,parai%allgrp)
    CALL mp_bcast(rng_block%mt,SIZE(rng_block%mt),parai%source,parai%allgrp)

200 FORMAT(/' SEEDING RANDOM NUMBER GENERATOR FOR SURFACE HOPPING.',&
         /,' REF: M. Matsumoto and T. Nishimura,',&
         /,' ACM Transactions on Modeling and Computer Simulation,',&
         /,' Vol. 8, No. 1, January 1998, pp 3--30.',/)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE shop_sgrnd
  ! ==================================================================

END MODULE shop_adds_utils
