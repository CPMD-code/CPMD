MODULE wannier_center_utils
  USE adat,                            ONLY: elem
  USE cdftmod,                         ONLY: cdfthda,&
                                             cdftlog
  USE cnst,                            ONLY: fbohr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE hfxmod,                          ONLY: hfxc3,&
                                             wcentx
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mm_dimmod,                       ONLY: clsaabox
  USE mp_interface,                    ONLY: mp_bcast
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc3
  USE pimd,                            ONLY: ipcurr,&
                                             np_low
  USE response_pmod,                   ONLY: nmr_para,&
                                             response1
  USE ropt,                            ONLY: iteropt
  USE store_types,                     ONLY: cprint,&
                                             iprint_wann,&
                                             rout1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             parm
  USE vdw_utils,                       ONLY: wwannier
  USE vdwcmod,                         ONLY: nfragx,&
                                             nwfcx,&
                                             rwann,&
                                             swann,&
                                             tauref,&
                                             twannupx,&
                                             vdwl,&
                                             vdwwfl
  USE wann,                            ONLY: wannc,&
                                             wanni,&
                                             wannr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wannier_center
  !public :: wc_trajectory

CONTAINS

  ! ==================================================================
  SUBROUTINE wannier_center(xyzmat,ldx,nstate,center,tau0)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: xyzmat(:,:,:)
    INTEGER                                  :: ldx, nstate
    REAL(real_8)                             :: center(:,:), tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wannier_center'
    REAL(real_8), PARAMETER                  :: zero = 0.0_real_8

    INTEGER                                  :: i, ierr, ipx, istep, k, l
    INTEGER, SAVE                            :: lcounter = 0
    REAL(real_8) :: funcv, phasex, phasey, phasez, pi, pi24, &
      quadrupole(6,nstate), s(3), wsum, x, xmatim, xmatre, y, ymatim, ymatre, &
      z, zmatim, zmatre

! ==--------------------------------------------------------------==

    lcounter = lcounter +1
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       pi=dacos(-1._real_8)
       pi24=4._real_8*pi*pi
       DO i=1,nstate
          IF (wanni%w_type.EQ.1) THEN
             funcv=0._real_8
             wsum=0._real_8
             DO k=1,wannc%nwanopt
                funcv=funcv+&
                     wannc%wwei(k)*(REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,k)))
                wsum=wsum+wannc%wwei(k)
             ENDDO
             center(4,i)=SQRT((wsum-funcv)/pi24)
          ELSE
             funcv=0._real_8
             DO k=1,wannc%nwanopt
                funcv=funcv-&
                     wannc%wwei(k)*LOG(REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,k)))
             ENDDO
             center(4,i)=SQRT(funcv/pi24)
          ENDIF
          ! RV----quadrupole------------
          DO l=1,6
             funcv=0._real_8
             DO k=1,wannc%nwanopt
                funcv=funcv-&
                     wannc%wquadi(l,k)*LOG(REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,k)))
             ENDDO
             quadrupole(l,i)=funcv/pi24
          ENDDO
          ! RV----quadrupole------------
          xmatre=REAL(xyzmat(i,i,1))
          xmatim=AIMAG(xyzmat(i,i,1))
          phasex=ATAN(xmatim/xmatre)
          IF (xmatre.LT.0._real_8) phasex=phasex+pi
          s(1)=-phasex/2._real_8/pi
          ymatre=REAL(xyzmat(i,i,2))
          ymatim=AIMAG(xyzmat(i,i,2))
          phasey=ATAN(ymatim/ymatre)
          IF (ymatre.LT.0._real_8) phasey=phasey+pi
          s(2)=-phasey/2._real_8/pi
          zmatre=REAL(xyzmat(i,i,3))
          zmatim=AIMAG(xyzmat(i,i,3))
          phasez=ATAN(zmatim/zmatre)
          IF (zmatre.LT.0._real_8) phasez=phasez+pi
          s(3)=-phasez/2._real_8/pi
          ! PBC the centers in the unit box
          s(1)=s(1)-dnint(s(1))
          s(2)=s(2)-dnint(s(2))
          s(3)=s(3)-dnint(s(3))
          ! hit with HMAT
          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,s,1,0.0_real_8,center(1,i),1)
          ! SUBTRACT OFF REFERENCE
          center(1,i)=center(1,i)-wannr%w_ref(1)
          center(2,i)=center(2,i)-wannr%w_ref(2)
          center(3,i)=center(3,i)-wannr%w_ref(3)
       ENDDO
       istep=cnti%npdip*lcounter
       CALL wc_trajectory(istep,center,quadrupole,tau0,nstate)
       istep=MOD(cnti%npdip*lcounter,cprint%iprint_step)
       IF (cprint%iprint(iprint_wann).EQ.1.AND.istep.LT.cnti%npdip) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(/,1X,64("*"))')
             WRITE(6,'(" *",12X,A,18X,A,1X,"*")')&
                  ' WANNIER CENTERS  ','<X^2> - <X>^2'
             WRITE(6,'(1X,64("*"))')
             DO i=1,nstate
                WRITE(6,'(3F12.4,17X,F12.4)')&
                     center(1,i),center(2,i),center(3,i),center(4,i)
             ENDDO
             WRITE(6,'(1X,64("*"),/)')
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (hfxc3%twscr) THEN
       IF (lcounter.EQ.1) THEN
          ALLOCATE(wcentx(4,nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(4*nstate,zero,0,wcentx,1)
       ENDIF
       IF (paral%io_parent) CALL dcopy(4*nstate,center,1,wcentx,1)
       CALL mp_bcast(wcentx,SIZE(wcentx),parai%io_source,parai%cp_grp)
    ENDIF
    IF (vdwl%vdwd) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr-np_low+1
       ELSE
          ipx=1
       ENDIF
       IF (cntl%tpath) vdwwfl%twannup=twannupx(ipx)
       IF (vdwwfl%twannup) THEN
          IF (paral%io_parent) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,tauref(1,1,1,ipx),1)
             DO i=1,nstate
                CALL pbc3(center(1,i),center(2,i),center(3,i),&
                     x,y,z,1,parm%apbc,parm%ibrav)
                rwann(1,i,ipx)=x+wannr%w_ref(1)-clsaabox%mm_c_trans(1)
                rwann(2,i,ipx)=y+wannr%w_ref(2)-clsaabox%mm_c_trans(2)
                rwann(3,i,ipx)=z+wannr%w_ref(3)-clsaabox%mm_c_trans(3)
                swann(i,ipx)=center(4,i)
             ENDDO
             CALL wwannier(nstate,tau0,rwann(:,:,ipx),swann(:,ipx))
          ENDIF
          CALL mp_bcast(tauref(:,:,:,ipx),3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
          CALL mp_bcast(rwann(:,:,ipx),3*nwfcx*nfragx,parai%io_source,parai%cp_grp)
          CALL mp_bcast(swann(:,ipx),nwfcx*nfragx,parai%io_source,parai%cp_grp)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wannier_center
  ! ==================================================================
  SUBROUTINE wc_trajectory(istep,center,quadrupole,tau0,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: istep
    REAL(real_8)                             :: center(:,:), quadrupole(:,:), &
                                                tau0(:,:,:)
    INTEGER                                  :: nstate

    CHARACTER(len=30)                        :: filen
    INTEGER                                  :: i, ia, ipx, is, j
    INTEGER, SAVE                            :: if1 = fo_verb, if2 = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: x, xx, y, yy, z, zz

! ==--------------------------------------------------------------==
! dsebasti

    IF (response1%tnmr) THEN
       IF (.NOT. nmr_para%nmr_superparent) RETURN
    ENDIF
    IF (response1%tepr) RETURN

    ! write only if commensurate with SAMPLE value.
    IF (MOD(iteropt%nfi-1,cnti%npdip).NE.0) RETURN

    IF (.NOT.paral%io_parent) RETURN

    ! suffix for replica
    IF (tmw) THEN
       ipx=mwi%walker_id
    ELSE IF (cntl%tpath) THEN
       ipx=ipcurr
    ENDIF
    ! dsebasti
    ! AK 2005/06/23: follow the output format selection scheme
    ! for the trajectory files. i.e. write a
    ! WANNIER_CENTER file, only if a TRAJECTORY is written...
    IF (rout1%rout.OR.rout1%xtout.OR.rout1%dcout) THEN
       filen='WANNIER_CENTER'
       IF (tmw.OR.cntl%tpath) THEN
          CALL mw_filename('WANNIER_CENTER_',filen,ipx)
       ENDIF
       CALL fileopen(4,filen,fo_app+if2,ferror)
       DO i=1,nstate
          WRITE(4,'(I7,3F18.8,F12.6)') iteropt%nfi,(center(j,i),j=1,4)
       ENDDO
       CALL fileclose(4)
       ! stores the quadrupole
       filen='WC_QUAD'
       IF (tmw.OR.cntl%tpath) THEN
          CALL mw_filename('WC_QUAD_',filen,ipx)
       ENDIF
       CALL fileopen(84,filen,fo_app+if2,ferror)
       DO i=1,nstate
          WRITE(84,'(I7,1x,6F14.6)') iteropt%nfi,(quadrupole(j,i),j=1,6)
       ENDDO
       CALL fileclose(84)
       if2=0
    ENDIF
    ! AK 2005/06/23: FIXME: for the time being always write IONS+CENTERS.xyz 
    ! in the long run this should be coupled to XTOUT.
    ! 
    ! write file "IONS+CENTERS" like in version 3.0 
    ! change of format to XYZ
    IF (cdftlog%thda)THEN
       IF (cdfthda%hdafirst)THEN
          filen='IONS+CENTERS-S1.xyz'
          IF (tmw.OR.cntl%tpath) THEN
             CALL mw_filename('IONS+CENTERS-S1_',filen,ipx)
             filen=TRIM(filen)//'.xyz'
          ENDIF
       ELSE
          filen='IONS+CENTERS-S2.xyz'
          IF (tmw.OR.cntl%tpath) THEN
             CALL mw_filename('IONS+CENTERS-S2_',filen,ipx)
             filen=TRIM(filen)//'.xyz'
          ENDIF
       ENDIF
    ELSE
       filen='IONS+CENTERS.xyz'
       IF (tmw.OR.cntl%tpath) THEN
          CALL mw_filename('IONS+CENTERS_',filen,ipx)
          filen=TRIM(filen)//'.xyz'
       ENDIF
    ENDIF
    CALL fileopen(44,filen,fo_app+if1,ferror)
    WRITE(44,*) ions1%nat+nstate
    WRITE(44,*) iteropt%nfi

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          ! mb Set the ions according to WFCs reference and apply pbc
          ! jh always apply pbc, even for isolated systems
          ! ak ...and then shift back to the original position.
          xx=tau0(1,ia,is)-wannr%w_ref(1)
          yy=tau0(2,ia,is)-wannr%w_ref(2)
          zz=tau0(3,ia,is)-wannr%w_ref(3)
          CALL pbc3(xx,yy,zz,x,y,z,1,parm%apbc,parm%ibrav)
          x=x+wannr%w_ref(1)-clsaabox%mm_c_trans(1)
          y=y+wannr%w_ref(2)-clsaabox%mm_c_trans(2)
          z=z+wannr%w_ref(3)-clsaabox%mm_c_trans(3)
          WRITE(44,'(A4,3F18.8)') elem%el(ions0%iatyp(is)),&
               x/fbohr, y/fbohr, z/fbohr
       ENDDO
    ENDDO
    DO i=1,nstate
       CALL pbc3(center(1,i),center(2,i),center(3,i),&
            x,y,z,1,parm%apbc,parm%ibrav)
       x=x+wannr%w_ref(1)-clsaabox%mm_c_trans(1)
       y=y+wannr%w_ref(2)-clsaabox%mm_c_trans(2)
       z=z+wannr%w_ref(3)-clsaabox%mm_c_trans(3)
       WRITE(44,'(A4,3F18.8)') "X", x/fbohr, y/fbohr, z/fbohr
    ENDDO
    CALL fileclose(44)

    ! store the spread for the case we have only the .xyz files.
    filen='WC_SPREAD'
    IF (tmw.OR.cntl%tpath) THEN
       CALL mw_filename('WC_SPREAD_',filen,ipx)
    ENDIF
    CALL fileopen(4,filen,fo_app+if1,ferror)
    WRITE (4,'(2I8)') iteropt%nfi, nstatE
    WRITE(4,'(6F12.6)') (center(4,i),i=1,nstate)
    CALL fileclose(4)
    if1=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wc_trajectory
  ! ==================================================================

END MODULE wannier_center_utils
