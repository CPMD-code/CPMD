MODULE cl_init_utils
  USE clas,                            ONLY: &
       clab, clas3, clas4, clas7, defy, imovc, myclat1, myclat2, ndefo, &
       ntrac, pote1, pr12, rcr12
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readff_utils,                    ONLY: readff
  USE readsr_utils,                    ONLY: xstring
  USE rmas,                            ONLY: rmass
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cl_init
  !public :: clean_label

CONTAINS

  ! ==================================================================
  SUBROUTINE cl_init(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &CLASSIC &END                ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &CLASSIC                                                 ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    FORCE FIELD                                               ==
    ! ==      ......                                                  ==
    ! ==    END FORCE FIELD                                           ==
    ! ==                                                              ==
    ! ==    FREEZE QUANTUM                                            ==
    ! ==    PRINT COORDINATES                                         ==
    ! ==    PRINT FF                                                  ==
    ! ==    FULL TRAJECTORY                                           ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'cl_init'

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ierr, is, j, n1, n2
    REAL(real_8)                             :: xat, xsaim, xsnow

! ==--------------------------------------------------------------==
! ==  DEFAULTS                                                    ==
! ==--------------------------------------------------------------==

    clas7%tfreeze=.FALSE.
    clas7%twcoor=.FALSE.
    clas7%tpff=.FALSE.
    clas7%tftra=.FALSE.
    clas7%tdefo=.FALSE.
    clas7%tmovc=.FALSE.
    ! NTRAC frequency of writing FF trajectory
    ! NDEFO frequency of applying deformation
    ntrac=1
    ndefo=1
    imovc=1
    ! ==--------------------------------------------------------------==
    ! ==  FORCE FIELDS (DEFAULTS)                                     ==
    ! ==--------------------------------------------------------------==
    pote1%t_r12=.FALSE.
    pote1%initr12=0
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       ! ===--------------------------------------------------------------==
       CALL clean_label
       ierr=inscan(iunit,'&CLAS')
       IF (ierr.NE.0) GOTO 100
       ! ==--------------------------------------------------------------==
10     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,err=99,END=99,fmt='(A80)') line
       IF (INDEX(line,'&END').NE.0) GOTO 100
       IF (INDEX(line,'FORCE').NE.0 .AND. INDEX(line,'FIELD').NE.0)THEN
          ! ..Classical Force Field
          CALL readff(iunit)
          GOTO 10
       ELSE IF (INDEX(line,'FREEZE').NE.0&
            .AND. INDEX(line,'QUANT').NE.0) THEN
          ! ..freeze quantum particles
          clas7%tfreeze=.TRUE.
          GOTO 10
       ELSE IF (INDEX(line,'PRINT').NE.0&
            .AND. INDEX(line,'COORD').NE.0) THEN
          ! ..print all coordinates
          clas7%twcoor=.TRUE.
          GOTO 10
       ELSE IF (INDEX(line,'PRINT').NE.0&
            .AND. INDEX(line,'FF').NE.0) THEN
          ! ..print force fields
          clas7%tpff=.TRUE.
          GOTO 10
       ELSE IF (INDEX(line,'FULL').NE.0&
            .AND. INDEX(line,'TRAJ').NE.0) THEN
          ! ..save full trajectory

          clas7%tftra=.TRUE.
          IF (INDEX(line,'SAMPLE').NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*) ntrac
          ENDIF
          GOTO 10
          ! ==--------------------------------------------------------------==
       ELSE IF (INDEX(line,'MOVIE').NE.0) THEN
          clas7%tmovc=.TRUE.
          IF (paral%io_parent)&
               READ(iunit,err=99,END=99,fmt=*) imovc
          GOTO 10
          ! ==--------------------------------------------------------------==
       ELSE
          IF (INDEX(line,'DEFORM').NE.0) THEN
             clas7%tdefo=.TRUE.
             IF (INDEX(line,'SAMPLE').NE.0) THEN
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) ndefo
             ENDIF
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*) defy%expx, defy%expy, defy%expz
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*) defy%xcen, defy%ycen, defy%zcen
          ENDIF
          ! ..Dummy line
          GOTO 10
       ENDIF
       ! ==--------------------------------------------------------------==
99     CONTINUE
       IF (paral%io_parent)&
            WRITE(6,*) ' CL_INIT : ERROR IN READING INPUT FILE'
       CALL stopgm('CL_INIT',' ',& 
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
100    CONTINUE
       ! ==--------------------------------------------------------------==
       ! ==  Test of options                                             ==
       ! ==--------------------------------------------------------------==
       ! ==--------------------------------------------------------------==
       ! ==  WRITE INFO TO OUTPUT                                        ==
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (clas7%tfreeze.AND.paral%io_parent)&
            WRITE(6,'(A)') ' QUANTUM PARTICLES WILL NOT MOVE'
       IF (clas7%tftra.AND.paral%io_parent)&
            WRITE(6,'(A,I6,A)')&
            ' FULL TRAJECTORY WILL BE SAVED EVERY',&
            ntrac,' STEPS'
       IF (clas7%tmovc.AND.paral%io_parent)&
            WRITE(6,'(A,I6,A)')&
            ' FULL TRAJECTORY WILL BE SAVED EVERY',&
            imovc,' STEPS IN MOVIE FORMAT'
       IF (clas7%tdefo.AND.paral%io_parent)&
            WRITE(6,'(A,I6,A)')&
            ' DEFORMATION WILL BE APPLIED EVERY',&
            ndefo,' STEPS'
       IF (clas7%tpff) THEN
          IF (pote1%t_r12) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' Force field of type: (a/R)^12'
             DO i=1,clas3%ncltyp
                DO j=1,i
                   IF (paral%io_parent)&
                        WRITE(6,'(6X,A,A,5X,A,A,5X,A,F10.4,A,F8.2)')&
                        'ATOM=',clab(i),'ATOM=',clab(j),&
                        'PARAM=',pr12(i,j),'   RC=',rcr12(i,j)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*)
       ! ==--------------------------------------------------------------==
    ENDIF
    ! ..distribute data
    CALL mp_bcast(clas7%tfreeze,parai%source,parai%allgrp)
    CALL mp_bcast(clas7%twcoor,parai%source,parai%allgrp)
    CALL mp_bcast(clas7%tpff,parai%source,parai%allgrp)
    CALL mp_bcast(clas7%tftra,parai%source,parai%allgrp)
    CALL mp_bcast(clas7%tdefo,parai%source,parai%allgrp)
    CALL mp_bcast(ndefo,parai%source,parai%allgrp)
    CALL mp_bcast_byte(defy, size_in_bytes_of(defy),parai%source,parai%allgrp)
    CALL mp_bcast(pote1%t_r12,parai%source,parai%allgrp)
    IF (pote1%t_r12) THEN
       IF (pote1%initr12.EQ.0) THEN
          ALLOCATE(pr12(clas3%ncltyp,clas3%ncltyp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rcr12(clas3%ncltyp,clas3%ncltyp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          pote1%initr12=1
       ENDIF
       CALL mp_bcast(pr12,SIZE(pr12),parai%source,parai%allgrp)
       CALL mp_bcast(rcr12,SIZE(rcr12),parai%source,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ..distribute work for parallel jobs
    xat=REAL(clas3%nclatom,kind=real_8)
    xsnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xsaim = xsnow + xat/parai%nproc
       n1=NINT(xsnow)+1
       n2=NINT(xsaim)
       IF (NINT(xsaim).GT.clas3%nclatom) n2=clas3%nclatom
       IF (i.EQ.1) n2=clas3%nclatom
       xsnow = xsaim
       IF (i-1.EQ.parai%mepos) THEN
          myclat1=n1
          myclat2=n2
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ..update mass information
    DO i=1,clas3%ncltyp
       is=clas3%is_qm(i)
       IF (is.GT.0) clas4%clmas(i)=rmass%pma0(is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cl_init
  ! ==================================================================
  SUBROUTINE clean_label

    INTEGER                                  :: i, i1, i2, len
    CHARACTER(len=4)                         :: temp

! ==--------------------------------------------------------------==

    DO i=1,clas3%ncltyp
       temp='    '
       CALL xstring(clab(i),i1,i2)
       len=i2-i1+1
       IF (len.GT.0) THEN
          temp(1:len)=clab(i)(i1:i2)
       ENDIF
       clab(i)=temp
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE clean_label
  ! ==================================================================

END MODULE cl_init_utils
