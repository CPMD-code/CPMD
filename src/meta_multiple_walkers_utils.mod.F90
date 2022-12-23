MODULE meta_multiple_walkers_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_dyn, cv_ist, cv_path, cv_vel, ekincv, &
       ekincv_walk, fmtdres, hllh_val, hllw_val, iangcv, imeta, inter_hill, &
       inter_hill_max, kharm, lmeta, ncolvar, rcc, rmeta, vharm_walk
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def,&
                                             fo_old,&
                                             fo_ufo,&
                                             fo_verb
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_group,&
                                             mp_split,&
                                             mp_sync
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             nproc_tot,&
                                             parentgroup,&
                                             pc_groups,&
                                             pc_grp,&
                                             pcg_pos,&
                                             supergroup,&
                                             supersource
  USE prng_utils,                      ONLY: repprngu
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: iteropt
  USE set_cp_grp_utils,                ONLY: reset_cp_grp
  USE system,                          ONLY: cntr,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wann,                            ONLY: wan05,&
                                             wannl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mw_init
  PUBLIC :: mw_filename
  PUBLIC :: rmtdres_mw
  PUBLIC :: calc_pos_dyn_mw
  PUBLIC :: wmtdres_mw
  PUBLIC :: rinvelcv_mw
  PUBLIC :: cv_exlagr_out_mw
  PUBLIC :: mw_assign_filenames

CONTAINS

  ! ================================================================== 
  SUBROUTINE calc_pos_dyn_mw(tollm,i_cvst_mw,i_temp_mw,i_meta,hc_last,hill_add_mw)

    REAL(real_8)                             :: tollm
    INTEGER                                  :: i_cvst_mw(mwi%nwalk), &
                                                i_temp_mw(mwi%nwalk), i_meta
    REAL(real_8)                             :: hc_last(ncolvar*mwi%nwalk)
    LOGICAL                                  :: hill_add_mw(mwi%nwalk)

    CHARACTER(len=50)                        :: decmw(3)
    INTEGER                                  :: icv, ipw0, ipw1, iu, iwalk0, &
                                                iwalk1, walknear
    REAL(real_8)                             :: av_disp, &
                                                av_disp_rel(mwi%nwalk), &
                                                cvdist, diff

    decmw(1)=' MIN.MTD.STEP '
    decmw(2)=' DISP.GT.TOL  '
    decmw(3)=' MAX.MTD.STEP '

    DO iwalk0=1,mwi%nwalk
       ! position of first cv of a walker
       ipw0=(iwalk0-1)*ncolvar
       hill_add_mw(iwalk0)=.FALSE.
       ! check the inter meta step is >= minimum mtd step
       IF (i_cvst_mw(iwalk0).GE.inter_hill)THEN
          i_temp_mw(iwalk0)=i_temp_mw(iwalk0)+1
          ! Displacement of a walker is the minimum displacement from the
          ! last hills put by all walkers             
          av_disp_rel(iwalk0)=99999.0! a big value
          walknear=0
          DO iwalk1=1,mwi%nwalk
             ! position of first cv of a walker
             ipw1=(iwalk1-1)*ncolvar
             av_disp=0.0_real_8
             DO icv=1,ncolvar
                ! Displacement of a cv
                diff=(cv_dyn(icv+ipw0)-hc_last(icv+ipw1))
                IF (iangcv(icv) .EQ. 1 .AND. diff .GT. pi) THEN
                   diff = diff -2.0_real_8*pi
                ELSEIF (iangcv(icv) .EQ. 1 .AND. diff .LT. -pi) THEN
                   diff = diff +2.0_real_8*pi
                ENDIF
                diff=diff/cscl_fac(1,icv)
                ! Displacement square
                av_disp=av_disp+diff*diff
             ENDDO
             ! average Displacement
             av_disp=SQRT(av_disp)
             IF (av_disp.LT.av_disp_rel(iwalk0))THEN
                av_disp_rel(iwalk0)=av_disp
                walknear=iwalk1
             ENDIF
          ENDDO
          IF (walknear.EQ.0)CALL stopgm('CALC_POS_DYN_MW',&
               'WALKNEAR.EQ.0',& 
               __LINE__,__FILE__)
          ! deciding if hill need to be added due to any conditions
          iu=0
          ! (a) add hill for the first walker, if the first time 
          ! minimum mtd step is reached
          IF ((i_cvst_mw(iwalk0).EQ.inter_hill).AND.&
               (i_meta.EQ.1).AND.&
               ((iwalk0).EQ.1) )THEN
             hill_add_mw(iwalk0)=.TRUE.
             iu=1
             ! (b) if max mtd step reached
          ELSE IF (i_temp_mw(iwalk0).EQ.inter_hill_max)THEN
             hill_add_mw(iwalk0)=.TRUE.
             iu=3
             ! (c) if displacement need to be checked and displacement >= tolerance 
             ! plus the constraint that distance between CVs at time t is also >=tolerance
          ELSE IF((MOD(i_temp_mw(iwalk0),imeta%icheck).EQ.0).AND.&
               (av_disp_rel(iwalk0).GE.tollm))THEN
             hill_add_mw(iwalk0)=.TRUE.
             iu=2
             IF (iwalk0.GT.1)THEN
                DO iwalk1=1,iwalk0-1
                   ipw1=(iwalk1-1)*ncolvar
                   cvdist=0.0_real_8
                   DO icv=1,ncolvar
                      cvdist=cvdist+&
                           (cv_dyn(icv+ipw0)-cv_dyn(icv+ipw1))**2
                   ENDDO
                   cvdist=SQRT(cvdist)
                   IF (cvdist.LE.tollm) THEN
                      hill_add_mw(iwalk0)=.FALSE.
                      iu=0
                      IF (paral%io_parent)&
                           WRITE(6,'(A,I8,A,I8,A,F16.6,A)')&
                           '**Walker',iwalk0,' is close to',iwalk1,&
                           ' Distance:',cvdist,' **'
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF (iu.GT.0)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(T2,A,I8,A,I8,A,F16.8,A,F16.8,A)')&
                  ' WALKER:',iwalk0,&
                  ' NEIGHBOR:',walknear,&
                  ' ||DISP.tot|| = ', av_disp_rel(iwalk0),&
                  '  Tolerance ',tollm,&
                  decmw(iu)
          ENDIF
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE calc_pos_dyn_mw

  ! ==================================================================
  SUBROUTINE rinvelcv_mw(vel,mass,cscl_fac,nvar,temp,ekincv)
    INTEGER                                  :: nvar
    REAL(real_8)                             :: cscl_fac(3,nvar), mass(nvar), &
                                                vel(nvar), temp, ekincv

    INTEGER                                  :: icv
    REAL(real_8)                             :: alfa, const, rnr, sigma, &
                                                tempp, tscal, vscale

    CALL zeroing(vel)!,nvar)
    IF (temp.LT.1.e-5_real_8) GOTO 100
    DO icv=1,nvar
       sigma=SQRT(temp/(mass(icv)*factem))
       rnr = repprngu()
       alfa=2.0_real_8*pi*rnr
       vel(icv)&
            =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa)*sigma
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE KINETIC ENERGY OF THE CV                          ==
    ! ==--------------------------------------------------------------==
    ekincv= 0._real_8
    tempp = 0._real_8
    DO icv=1,nvar
       const=0.5_real_8*mass(icv)*vel(icv)*vel(icv)
       ekincv=ekincv+const
       tempp = tempp + const/(cscl_fac(1,icv)*cscl_fac(1,icv))
    ENDDO
    tempp=tempp*factem*2.0_real_8/REAL(nvar,kind=real_8)
    IF (tempp.GT.1.e-5_real_8) THEN
       tscal=temp/tempp
       vscale=SQRT(tscal)
       ! !$OMP parallel do private(ICV)
       DO icv=1,nvar
          vel(icv)=vel(icv)*vscale
       ENDDO
    ENDIF

100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinvelcv_mw

  ! ==================================================================
  SUBROUTINE  cv_exlagr_out_mw(i_meta,hc_last,f_harm,f_hill,f_wall,&
       cvd_ist,Displacement,ekinc,ekinp,iwalk1)
    INTEGER                                  :: i_meta
    REAL(real_8) :: hc_last(ncolvar), f_harm(ncolvar), f_hill(ncolvar), &
      f_wall(ncolvar), cvd_ist(ncolvar), Displacement, ekinc, ekinp
    INTEGER                                  :: iwalk1

    CHARACTER(len=20), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file3 = 'istvar_mtd', file4 = 'forfac_mtd', &
      file5 = 'disvar_mtd', file6 = 'velvar_mtd', file7 = 'enevar_mtd', &
      file8 = 'spinpo_mtd'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    INTEGER                                  :: iaa, icv, iee
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: t_ion

    IF (paral%io_parent)&
         CALL fileopen(51,file1,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(52,file2,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(53,file3,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(54,file4,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(55,file5,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(56,file6,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(57,file7,fo_app+ifirst,ferror)
    IF (lmeta%tlocalizespin) THEN
       IF (paral%io_parent)&
            CALL fileopen(58,file8,fo_app+ifirst,ferror)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Verbose open only on the first run
    ifirst=0
    ! ==--------------------------------------------------------------==

    ! Write output
    ! Check if rmeta%hllh /= 0.0 - In case not, we do not write the 
    ! hill information, because it was not added..
    IF (rmeta%hllh /=0.0_real_8) THEN
       IF (paral%io_parent)&
            WRITE(chnum,'(I5)') ncolvar
       CALL xstring(chnum,iaa,iee)
       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E13.3,2I8)'
       IF (paral%io_parent)&
            WRITE(51,lineform) iteropt%nfi,(hc_last(icv),icv=1,ncolvar),&
            (cscl_fac(1,icv),icv=1,ncolvar),iwalk1,i_meta

       lineform =&
            '(1X,I7,3F14.6,2I8)'
       IF (paral%io_parent)&
            WRITE(52,lineform) iteropt%nfi,&
            Displacement,rmeta%hllw,rmeta%hllh,iwalk1,i_meta

       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
       IF (paral%io_parent)&
            WRITE(53,lineform) iteropt%nfi,(cv_ist(icv),icv=1,ncolvar),&
            (cv_ist(icv)-hc_last(icv),icv=1,ncolvar)

       ! Print force factors
       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'&
            //chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E16.4)'
       IF (paral%io_parent)&
            WRITE(54,lineform) iteropt%nfi,(f_harm(icv),icv=1,ncolvar),&
            (f_hill(icv),icv=1,ncolvar),(f_wall(icv),icv=1,ncolvar)

       ! Print last displacements
       lineform =&
            '(I6,'//chnum(iaa:iee)//'f11.6,'//chnum(iaa:iee)//&
            'f10.6,'//chnum(iaa:iee)//'f10.6)'
       IF (paral%io_parent)&
            WRITE(55,lineform) iteropt%nfi,&
            ((cv_dyn(icv)-hc_last(icv)),icv=1,ncolvar),&
            (cvd_ist(icv),icv=1,ncolvar),&
            (kharm(icv),icv=1,ncolvar)

       ! Print last velocities, temperature and kinetic energy
       ! EHAM_HILL=EKINP+ETOT+ENOSE+ENOSP+ECNSTR+EKINC+VHARM+EKINC
       IF (ALLOCATED(cv_vel)) THEN
          lineform =&
               '(I6,'//chnum(iaa:iee)//'f11.6,f15.6,6f16.8)'
          IF (paral%io_parent)&
               WRITE(56,lineform) iteropt%nfi,(cv_vel(icv),icv=1,ncolvar)
       END IF


       t_ion = ekinp*factem*2._real_8/glib
       IF (ALLOCATED(ekincv_walk)) THEN
          lineform = '(I6,8E16.6)'
          IF (paral%io_parent)&
               WRITE(57,lineform) iteropt%nfi,&
               t_ion,ekinc,ekincv_walk(iwalk1),&
               vharm_walk(iwalk1),rmeta%gausspot,&
               ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
       ELSE
          lineform = '(I6,4E16.6)'
          IF (paral%io_parent)&
               WRITE(57,lineform) iteropt%nfi,&
               t_ion,ekinc, rmeta%gausspot,&
               ener_com%etot
       END IF


       IF (lmeta%tlocalizespin) THEN
          ! Print Position of the center of the spin density
          IF (paral%io_parent)&
               WRITE(58,'(I6,3f12.6)') iteropt%nfi, rcc(1),rcc(2),rcc(3)
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') " >> Hills in Metadynamics have an height of zero!"
          WRITE(6,'(A)') " >> No Hills will be spawned. Metadynamics files may be empty!"
       END IF
    END IF

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(51)
    IF (paral%io_parent)&
         CALL fileclose(52)
    IF (paral%io_parent)&
         CALL fileclose(53)
    IF (paral%io_parent)&
         CALL fileclose(54)
    IF (paral%io_parent)&
         CALL fileclose(55)
    IF (paral%io_parent)&
         CALL fileclose(56)
    IF (paral%io_parent)&
         CALL fileclose(57)
    IF (lmeta%tlocalizespin.AND.paral%io_parent)&
         CALL fileclose(58)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cv_exlagr_out_mw

  ! ==================================================================
  SUBROUTINE wmtdres_mw(nw,ntot_iter,i_meta,iwalk1)

    INTEGER                                  :: nw, ntot_iter, i_meta, iwalk1

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, ipw, it
    LOGICAL                                  :: ferror

    ipw=(iwalk1-1)*ncolvar
    IF (paral%parent) THEN
       CALL xstring(fmtdres,ia,ie)
       IF (paral%io_parent)&
            CALL fileopen(nw,fmtdres(ia:ie),fo_def+fo_ufo,ferror)
       IF (paral%io_parent)&
            REWIND(nw)
       IF (paral%io_parent)&
            WRITE(nw) i_meta-1
       IF (paral%io_parent)&
            WRITE(nw) ncolvar

       DO icv = 1,ncolvar
          ! full path of the cv and the scaling factors are written
          IF (paral%io_parent)&
               WRITE(nw) (cv_path(it,icv),it=1,i_meta-1)
          IF (paral%io_parent)&
               WRITE(nw) (cscl_val(it,icv),it=1,i_meta-1)
       ENDDO

       IF (lmeta%lextlagrange) THEN
          ! each walker writes the difference of s(t) with the s(hill_last) written
          IF (paral%io_parent)&
               WRITE(nw)(cv_dyn(icv+ipw)-cv_path(i_meta-1,icv),&
               icv=1,ncolvar)
          ! each walker writes its own velocity
          IF (paral%io_parent)&
               WRITE(nw)(cv_vel(icv+ipw),icv=1,ncolvar)
       ENDIF
       IF (paral%io_parent)&
            WRITE(nw) (hllw_val(it,1),it=1,i_meta-1)
       IF (paral%io_parent)&
            WRITE(nw) (hllh_val(it,1),it=1,i_meta-1)

       IF (paral%io_parent)&
            CALL fileclose(nw)
       IF (paral%io_parent)&
            WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
       IF (paral%io_parent)&
            WRITE(6,fformat)&
            ' MTD RESTART INFO WRITTEN ON FILE ',fmtdres(ia:ie)

    ENDIF  ! parent
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wmtdres_mw
  ! ==--------------------------------------------------------------==

  SUBROUTINE rmtdres_mw(nr,cv_disp,tfull,ntot_iter,iwalk1)
    INTEGER                                  :: nr
    REAL(real_8)                             :: cv_disp(ncolvar)
    LOGICAL                                  :: tfull
    INTEGER                                  :: ntot_iter, iwalk1

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, it, nc
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       CALL xstring(fmtdres,ia,ie)
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(nr,fmtdres(ia:ie),fo_old+fo_ufo,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)')&
               '(/,A,T',MAX(36,65-(ie-ia)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat)&
               'RMTDRES| MTD RESTART FILE NOT FOUND:',fmtdres(ia:ie)
          CALL stopgm('RMTDRES',' FILE NOT FOUND ',& 
               __LINE__,__FILE__)
       ENDIF

       IF (paral%io_parent)&
            REWIND(nr)
       IF (paral%io_parent)&
            READ(nr) imeta%i_meta_res
       IF (paral%io_parent)&
            READ(nr) nc

       IF (nc .NE. ncolvar) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I3,A,I3)') 'RMTDRES! DIFFERENT # of CV ',&
               nc,' vs',ncolvar
          CALL stopgm('RMTDRES','USE SAME # of CV',& 
               __LINE__,__FILE__)
       ENDIF

       IF (tfull) THEN
          DO icv = 1,ncolvar
             ! MW         each walker is reading the FULL hills; indivisual HC_LAST can be reconstructed from the diffusion
             IF (paral%io_parent)&
                  READ(nr) (cv_path(it,icv),it=1,imeta%i_meta_res)
             ! MW         each walker is reading the FULL scaling factors
             IF (paral%io_parent)&
                  READ(nr) (cscl_val(it,icv),it=1,imeta%i_meta_res)
          ENDDO

          IF (lmeta%lextlagrange) THEN
             ! MW         each walker is reading his own part of the displacement
             IF (paral%io_parent)&
                  READ(nr)(cv_disp(icv),icv=1,ncolvar)
             ! MW         each walker is reading his own velocity part only
             IF (paral%io_parent)&
                  READ(nr)(cv_vel(icv+(iwalk1-1)*ncolvar),icv=1,ncolvar)
          ENDIF

          IF (paral%io_parent)&
               READ(nr) (hllw_val(it,1),it=1,imeta%i_meta_res)
          IF (paral%io_parent)&
               READ(nr) (hllh_val(it,1),it=1,imeta%i_meta_res)
          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat)&
               ' MTD RESTART INFO READ FROM FILE ',fmtdres(ia:ie)
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(nr)

    ENDIF  ! parent
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmtdres_mw

  ! ==--------------------------------------------------------------==
  SUBROUTINE mw_filename(filenamein,filenameout,ipl)
    ! ==--------------------------------------------------------------==
    ! == Some multiple walker metadynamics output files have names    ==
    ! == in the format NAME_#id(.EXTENTION)                           ==
    ! ==--------------------------------------------------------------==
    ! 
    CHARACTER(len=*)                         :: filenamein, filenameout
    INTEGER                                  :: ipl

    CHARACTER(len=50)                        :: cipnum
    INTEGER                                  :: i1, i2, n1, n2

    filenameout=' '
    IF (paral%io_parent)&
         WRITE(cipnum,'(I4)') ipl
    CALL xstring(filenamein,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filenameout=filenamein(n1:n2)//cipnum(i1:i2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mw_filename

  ! ==--------------------------------------------------------------==
  SUBROUTINE mw_assign_filenames(ipcurr,filen,filenbs)
    ! ==--------------------------------------------------------------==
    ! == Create filenames: RESTART_#ip.  RESTART_#ip 
    ! ==                  ENERGIES_#ip  BS_ENERG_#ip 
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: ipcurr
    CHARACTER(len=*)                         :: filen, filenbs

    INTEGER                                  :: n1, n2

    IF (.NOT.tmw)RETURN
    CALL mw_filename('ENERGIES_',filen,ipcurr)
    CALL mw_filename('BS_ENERG_',filenbs,ipcurr)
    CALL mw_filename('RESTART_',filn,ipcurr)
    CALL xstring(filn,n1,n2)
    filbod=filn(n1:n2)//'.'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mw_assign_filenames

  ! ==--------------------------------------------------------------==
  SUBROUTINE mw_init
    ! ==--------------------------------------------------------------==
    ! ==    INITIALIZE MULTIPLE WALKER METADYNAMICS                   ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'mw_init'

    CHARACTER(len=12)                        :: fileon, fileout
    INTEGER                                  :: color, i, i1, i2, ip, ipp, &
                                                isub, metmp, nhigh, nlow, &
                                                pg_me
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: rproc

    IF (.NOT.tmw)RETURN
    CALL tiset(procedureN,isub)

    CALL mp_bcast(cntr%ecut,parai%io_source,parai%cp_grp) ! TODO is it needed?

    ! ==--------------------------------------------------------------==
    ! Reassigning PARENTS and SOURCE plus define GRANDPARENT and SUPERSOURCE
    grandparent=paral%io_parent
    supergroup=parai%cp_grp
    supersource=parai%io_source
    nproc_tot=parai%cp_nproc
    pg_me=parai%cp_me
    ! ==--------------------------------------------------------------==
    ! Generate processor groups
    CALL zeroing(parap%pgroup)!,maxcpu)
    CALL zeroing(parap%nlink)!,maxcpu)
    ! Always 1 walker per processor groups
    pc_groups=mwi%nwalk
    ! 
    rproc=REAL(parai%cp_nproc,kind=real_8)/mwi%nwalk
    DO i=1,mwi%nwalk
       nlow=NINT((i-1)*rproc)
       pc_grp(i)=nlow
       nhigh=NINT(i*rproc)-1
       IF (parai%cp_me.GE.nlow .AND. parai%cp_me.LE.nhigh) THEN
          parai%cp_nproc=nhigh-nlow+1! number of processors in a group is redifined on all corresponding processors
          parai%mepos=parai%cp_me-nlow! mepos is redifined on all processors of a group
          pcg_pos=i
          mwi%walker_id=i
          DO ip=1,parai%cp_nproc
             parap%pgroup(ip)=ip-1
             ipp=parap%pgroup(ip)
             parap%nlink(ipp)=ip-1
          ENDDO
       ENDIF
    ENDDO

    color=parai%cp_me/parai%cp_nproc
    CALL mp_split(supergroup,color,parai%cp_me,parai%cp_grp,metmp)
    IF (parai%mepos.NE.metmp) THEN
       CALL stopgm('MW_INIT','! MY_SPLIT RETURNED ERROR!',& 
            __LINE__,__FILE__)
    ENDIF
    parai%cp_me=metmp
    CALL mp_sync(supergroup)
    CALL mp_group(mwi%nwalk,pc_grp,parentgroup,supergroup)
    !
    !     Reset the cp group communicator
    !
    CALL reset_cp_grp()

    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND. pg_me+1>parai%cp_nproc) THEN
       WRITE(fileon,'(I3)') pcg_poS
       CALL xstring(fileon,i1,i2)
       fileout='OUTPUT_'//fileon(i1:i2)
       CALL fileclose(6)
       CALL fileOpen(6,fileout,fo_def,ferror)
#if defined(__Linux) || defined (__ALPHALINUX)
       ! avoid printmemsize messes
       CALL silentstdout
#endif
    ENDIF

    IF (grandparent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,'(1X,A4,A18,A38,A4)') '**  ','WALKERID', 'PROCESSORS',&
            '**'
       WRITE(6,'(1X,64("-"))')
       DO i=1,mwi%nwalk
          nlow=pc_grp(i)
          IF (i.EQ.mwi%nwalk) THEN
             nhigh=nproc_tot-1
          ELSE
             nhigh=pc_grp(i+1)-1
          ENDIF
          WRITE(6,'(1X,A4,I18,I17,A,I17,A4)') '**  ',i,&
               nlow,' -> ',nhigh, '**'
       ENDDO
       WRITE(6,'(1X,64("*"),/)')
    ENDIF

    IF (wannl%twann) THEN
       IF (wan05%loc_npgrp>parai%cp_nproc) wan05%loc_npgrp=parai%cp_nproc
    ENDIF

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE mw_init
  ! ==================================================================

END MODULE meta_multiple_walkers_utils
