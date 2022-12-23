MODULE setsys_utils
  USE adat,                            ONLY: elem
  USE atom,                            ONLY: gnl,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE bsym,                            ONLY: &
       autocm, bsclcs, bsfac, cnstwgt, hdoel, hsdown, hspin, hsup, hupel, &
       rtsasb, sa, sb, smax, smin
  USE cell,                            ONLY: cell_com
  USE chksym_utils,                    ONLY: chksym,&
                                             gentau,&
                                             givesym,&
                                             symtau
  USE clas,                            ONLY: tclas
  USE cnst,                            ONLY: fbohr,&
                                             ry,&
                                             scmass
  USE cnst_dyn,                        ONLY: atcvar,&
                                             cvpar,&
                                             ncolvar,&
                                             rccnga,&
                                             tycvar,&
                                             vbound
  USE coor,                            ONLY: lvelini,&
                                             tau0,&
                                             velp
  USE cotr,                            ONLY: &
       cnpar, cnsval, cnsval_dest, cotc0, cotr007, duat, grate, gsrate, &
       lskcor, ntcnst, ntrest, respar, resval, resval_dest
  USE cplngs_utils,                    ONLY: cpl_para
  USE dpot,                            ONLY: dpot_mod
  USE elct,                            ONLY: crge
  USE elct2,                           ONLY: tfixo
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE fint,                            ONLY: fint1,&
                                             fint4,&
                                             fint5
  USE hubbardu,                        ONLY: hubbu
  USE ions,                            ONLY: coord_fdiff,&
                                             ions0,&
                                             ions1,&
                                             r_fdiff
  USE isos,                            ONLY: isos1
  USE kdpc,                            ONLY: tkdp
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: kcnorm,&
                                             kgemax,&
                                             rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE metr,                            ONLY: metr_com
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_go_mm,&
                                             mm_go_qm,&
                                             mmdim,&
                                             cpsp
  USE mm_extrap,                       ONLY: nnow, numcold
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast
  USE multtb_utils,                    ONLY: multtb
  USE nlcc,                            ONLY: corecg,&
                                             corei,&
                                             corel,&
                                             corer,&
                                             rcgrid
  USE nlps,                            ONLY: nlm,&
                                             nlps_com
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE response_pmod,                   ONLY: response1
  USE rkpnt_utils,                     ONLY: bzmesh,&
                                             rkpnt
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: ropt_mod
  USE rv30_utils,                      ONLY: rdkpoints
  USE sfac,                            ONLY: natx
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE shock,                           ONLY: shock1
  USE shop,                            ONLY: fs0,&
                                             fs1,&
                                             sh02
  USE sphe,                            ONLY: tsphere
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2,&
                                             lspin3,&
                                             lspin4,&
                                             spin_mod,&
                                             tdsp1
  USE store_types,                     ONLY: restart1
  USE symm,                            ONLY: &
       irt, irtvec, isymu, numshel, sxscale, symmi, symmr, symmt, syscale, &
       szscale, tcartesian, tvec
  USE system,                          ONLY: &
       cnti, cntl, cntr, kpbeg, lx, maxsys, nbrx, nhx, nkpbl, nkpt, parm
  USE vdbp,                            ONLY: &
       betar, dion, ncpr1, qfunc, qqq, qrl, r, rab, rsatom, rscore, ru, rucore
  USE vdbt,                            ONLY: itmax,&
                                             vdbti
  USE velupi_utils,                    ONLY: s_to_c
  USE wann,                            ONLY: wannl,&
                                             wannr
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setsys

CONTAINS

  ! ==================================================================
  SUBROUTINE setsys
    ! ==--------------------------------------------------------------==
    ! ==  FINAL SETUP OF THE SYSTEM                                   ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'setsys'
    CHARACTER(len=1), DIMENSION(5), &
      PARAMETER                              :: cang = (/'S','P','D','F','G'/)
    CHARACTER(len=10), DIMENSION(2), PARAMETER :: &
      ctyp = (/'     LOCAL','  NONLOCAL'/)
    CHARACTER(len=12), DIMENSION(7), PARAMETER :: ischeme = (/'GAUSS-HERMIT',&
      '    KLEINMAN','  VANDERBILT','   GOEDECKER','         EAM',&
      '            ','     UPF  NC'/)
    CHARACTER(len=3), DIMENSION(2), &
      PARAMETER                              :: inlcc = (/'YES',' NO'/)
    CHARACTER(len=7), DIMENSION(10), PARAMETER :: mul = (/'SINGLET','DOUBLET',&
      'TRIPLET','QUARTET','QUINTET',' SEXTET',' SEPTET','  OCTET','  NONET',&
      ' DECTET'/)
    REAL(real_8), PARAMETER                  :: betaelmax = 1.e33_real_8 

    INTEGER :: i, ia, iat, ierr, ifac, ik, il, ipgrpa, is, isin, it, iu, ja, &
      js, k, m, mbl, n1, nelhalf, nkpoint, nnnx, nsic, NSX_q, nt, nx
    LOGICAL                                  :: fixedkblock, lfnotset, &
                                                lveltrue, status, statusdummy
    REAL(real_8)                             :: cx, cy, cz, dcx, dcy, dcz, &
                                                ddr, dmin, dr(3), dr_(3), ff, &
                                                htsave(3,3), rkt(3), sumwk, &
                                                xnel, ynel
    REAL(real_8), ALLOCATABLE                :: ft(:), taus(:,:,:)
    REAL(real_8), EXTERNAL                   :: dasum

    IF (lqmmm%qmmm)THEN
       NSX_q=mmdim%nspq
    ELSE
       NSX_q=maxsys%nsx
       mmdim%nspm =maxsys%nsx
    ENDIF
    IF (lqmmm%qmmm)CALL mm_dim(mm_go_mm,status)
    IF (.NOT.paral%io_parent) GOTO 9999
    IF (cntl%tscale) THEN
       ! ==--------------------------------------------------------------==
       ! == SCALE OPTIONS: ATOMIC COORDINATES IN DIRECT LATTICE BASIS    ==
       ! ==--------------------------------------------------------------==
       CALL dscal(maxsys%nax*maxsys%nsx,1._real_8/sxscale,tau0(1,1,1),3)
       CALL dscal(maxsys%nax*maxsys%nsx,1._real_8/syscale,tau0(2,1,1),3)
       CALL dscal(maxsys%nax*maxsys%nsx,1._real_8/szscale,tau0(3,1,1),3)
       IF (tcartesian) THEN
          CALL dscal(maxsys%nax*maxsys%nsx,cell_com%celldm(1),tau0(1,1,1),3)
          CALL dscal(maxsys%nax*maxsys%nsx,cell_com%celldm(1)*cell_com%celldm(2),tau0(2,1,1),3)
          CALL dscal(maxsys%nax*maxsys%nsx,cell_com%celldm(1)*cell_com%celldm(3),tau0(3,1,1),3)
       ELSE
          ALLOCATE(taus(3,maxsys%nax,ions1%nsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,taus,1)
          CALL s_to_c(taus,tau0)
          DEALLOCATE(taus,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ELSEIF (.NOT.cntl%bohr) THEN
       ! ==--------------------------------------------------------------==
       ! == ANGSTROM OPTION: ATOMIC COORDINATES IN cntl%bohr UNIT             ==
       ! ==--------------------------------------------------------------==
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             tau0(1,ia,is)=tau0(1,ia,is)*fbohr
             tau0(2,ia,is)=tau0(2,ia,is)*fbohr
             tau0(3,ia,is)=tau0(3,ia,is)*fbohr
             velp(1,ia,is)=velp(1,ia,is)*fbohr
             velp(2,ia,is)=velp(2,ia,is)*fbohr
             velp(3,ia,is)=velp(3,ia,is)*fbohr
          ENDDO
       ENDDO

       ! convert wannier reference coordinate to atomic units if set in input.
       IF (wannr%w_ref(1).GT.-999999.9_real_8) THEN
          DO i=1,3
             wannr%w_ref(i) = wannr%w_ref(i)*fbohr
          ENDDO
       ENDIF
       ! convert finite differences reference coordinate and radius to atomic units
       DO i=1,3
          coord_fdiff(i) = coord_fdiff(i)*fbohr
       ENDDO
       r_fdiff = r_fdiff*fbohr

       ! convert type 1 dummy atom coordinates to atomic units
       DO i=1,duat%ndat1
          duat%dummy1(1,i)=duat%dummy1(1,i)*fbohr
          duat%dummy1(2,i)=duat%dummy1(2,i)*fbohr
          duat%dummy1(3,i)=duat%dummy1(3,i)*fbohr
       ENDDO

       ! convert distance based constraints to atomic units
       IF (cotc0%mcnstr.GT.0) THEN
          DO i=1,cotc0%mcnstr
             IF (ntcnst(1,i).EQ.1 .OR.&
                  ntcnst(1,i).EQ.4 .OR.&
                  ntcnst(1,i).EQ.7 .OR.&
                  ntcnst(1,i).EQ.12) THEN
                IF (cnsval(i).NE.-999._real_8) cnsval(i)=cnsval(i)*fbohr
                grate(i)=grate(i)*fbohr
                IF (cnsval_dest(i).NE.-999._real_8)&
                     cnsval_dest(i)=cnsval_dest(i)*fbohr
             ENDIF
             IF (ntcnst(1,i).EQ.6 .OR. ntcnst(1,i).EQ.8 .OR.&
                  ntcnst(1,i).EQ.9 .OR. ntcnst(1,i).EQ.10) THEN
                IF (ntcnst(1,i).EQ.6.OR.ntcnst(1,i).EQ.8)THEN
                   cnpar(1,i) = cnpar(1,i)/fbohr
                ELSE
                   cnpar(1,i) = cnpar(1,i)*fbohr
                ENDIF
                cnpar(2,i) = cnpar(2,i)*fbohr
             ENDIF
          ENDDO
       ENDIF
       ! convert distance based restraints to atomic units
       IF (cotr007%mrestr.GT.0) THEN
          DO i=1,cotr007%mrestr
             IF (ntrest(1,i).EQ.1 .OR.&
                  ntrest(1,i).EQ.4 .OR.&
                  ntrest(1,i).EQ.7 .OR.&
                  ntrest(1,i).EQ.12) THEN
                IF (resval(i).NE.-999._real_8) resval(i)=resval(i)*fbohr
                gsrate(i)=gsrate(i)*fbohr
                IF (resval_dest(i).NE.-999._real_8)&
                     resval_dest(i)=resval_dest(i)*fbohr
             ENDIF
             IF (ntrest(1,i).EQ.6 .OR. ntrest(1,i).EQ.8 .OR.&
                  ntrest(1,i).EQ.9 .OR. ntrest(1,i).EQ.10) THEN
                IF (ntrest(1,i).EQ.6 .OR. ntrest(1,i).EQ.8)THEN
                   respar(1,i)=respar(1,i)/fbohr
                ELSE
                   respar(1,i)=respar(1,i)*fbohr
                ENDIF
                respar(2,i)=respar(2,i)*fbohr
             ENDIF
          ENDDO
       ENDIF

       ! convert metadynamics parameters.
       IF (ncolvar.GT.0) THEN
          ipgrpa=0
          DO i=1,ncolvar

             ! convert upper/lower bounds for distance based CVs
             IF ((tycvar(i).EQ.1) .OR. (tycvar(i).EQ.4)   & ! STRETCH, DIST
                  .OR. (tycvar(i).EQ.7) .OR. (tycvar(i).EQ.13)  & ! DIFFER, DISPL1
                  .OR. (tycvar(i).EQ.23)) THEN                  ! CELLSIDE
                vbound(1,i) = vbound(1,i)*fbohr
                vbound(3,i) = vbound(3,i)*fbohr
             ENDIF
             IF (tycvar(i).EQ.12) THEN   ! RMSD_AB
                vbound(1,i) = vbound(1,i)*fbohr*fbohr
                vbound(3,i) = vbound(3,i)*fbohr*fbohr
             ENDIF
             IF (tycvar(i).EQ.22) THEN   ! VOLVAR
                vbound(1,i) = vbound(1,i)*fbohr*fbohr*fbohr
                vbound(3,i) = vbound(3,i)*fbohr*fbohr*fbohr
             ENDIF

             ! convert length parameters in the definition of the CVs
             IF ((tycvar(i).EQ.6).OR.(tycvar(i).EQ.8)) THEN! COORD, COORSP
                cvpar(1,i) = cvpar(1,i)/fbohr
                cvpar(2,i) = cvpar(2,i)*fbohr
             ENDIF
             IF ((tycvar(i).EQ.10).OR.(tycvar(i).EQ.17)) THEN! BNSWT, DIFCOOR
                cvpar(1,i) = cvpar(1,i)*fbohr
             ENDIF
             IF ((tycvar(i).EQ.9).OR.(tycvar(i).EQ.11)) THEN! COOR_RF, TOT_COOR
                cvpar(1,i) = cvpar(1,i)*fbohr
                cvpar(2,i) = cvpar(2,i)*fbohr          ! 2SHELL
             ENDIF
             IF ((tycvar(i).EQ.16).OR.(tycvar(i).EQ.18)) THEN! HBONDCH, COOR_CHAIN
                cvpar(1,i) = cvpar(1,i)*fbohr
                cvpar(2,i) = cvpar(2,i)*fbohr
             ENDIF
             IF ((tycvar(i).EQ.29)) THEN! COORGROUP
                cvpar(1,i) = cvpar(1,i)/fbohr! Kappa is inverse of distance
                DO k=1,atcvar(1,i)
                   ipgrpa=ipgrpa+1
                   rccnga(ipgrpa)=rccnga(ipgrpa)*fbohr
                ENDDO
             ENDIF

             ! FIXME: conversion for HYDRONIUM and DIS_HYD still missing.
             ! is a bit tricky, since there are default parameters.
          ENDDO
       ENDIF

    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ATOMIC MASSES                                                ==
    ! ==--------------------------------------------------------------==
    IF (.NOT.lqmmm%qmmm) THEN
       rmass%pmatot=0._real_8
       rmass%pmat0=0._real_8
       ions1%nat=0
       DO is=1,ions1%nsp
          rmass%pma(is)=rmass%pma0(is)*scmass
          rmass%pmat0=rmass%pmat0+rmass%pma0(is)*ions0%na(is)
          rmass%pmatot=rmass%pmatot+rmass%pma(is)*ions0%na(is)
          ions1%nat=ions1%nat+ions0%na(is)
       ENDDO
       ! initialize cell mass to total mass, if not set in input
       ! variable cell is not available for qm/mm so we are safe.
       IF (cntr%cmass.LE.0.0_real_8) THEN
!!!          cntr%cmass=rmass%pmat0  ! wrong mass in chemical units
          cntr%cmass=rmass%pmatot  ! correct mass in atomic units
          IF (cntl%tsdc.OR.cntl%tprcp)&
             ! WRITE(6,'(A,T54,F12.2)')&
               WRITE(6,'(A,T54,E12.6)')&
               ' AUTOMATIC FICTITIOUS MD CELL MASS '&
               // 'SET TO [A.U.]:',cntr%cmass
       ENDIF
       IF (cntl%tshock) shock1%vshock=shock1%vshock*SQRT(rmass%pmatot)
       IF (ions1%nat.EQ.0) THEN
          CALL stopgm('SETSYS','THE NUMBER OF ATOMS IS NULL',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! safety check.
    IF ((cntl%tsdc.OR.cntl%tprcp).AND.(cntr%cmass.LT.1.0_real_8))&
         CALL stopgm('SETSYS','INVALID FICTITIOUS MD CELL MASS <1.0',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == MOVE COM TO CENTER OF BOX                                    ==
    ! ==--------------------------------------------------------------==
    clsaabox%mm_c_trans(1)=0.0_real_8
    clsaabox%mm_c_trans(2)=0.0_real_8
    clsaabox%mm_c_trans(3)=0.0_real_8
    ! EHR[
    IF (isos1%tclust.AND.isos1%tcent.AND..NOT.tclas&
         .AND..NOT.cntl%tmdeh.AND..NOT.cntl%tpspec) THEN
       ! EHR]
       IF (lqmmm%qmmm)THEN
#if defined (__GROMOS)
          CALL mm_center(tau0,clsaabox%mm_c_trans,.FALSE.)
#endif
       ELSE
          cx=0.0_real_8
          cy=0.0_real_8
          cz=0.0_real_8
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                cx=cx+tau0(1,ia,is)*rmass%pma0(is)/rmass%pmat0
                cy=cy+tau0(2,ia,is)*rmass%pma0(is)/rmass%pmat0
                cz=cz+tau0(3,ia,is)*rmass%pma0(is)/rmass%pmat0
             ENDDO
          ENDDO
          dcx=0.5_real_8*cell_com%celldm(1)-cx
          dcy=0.5_real_8*cell_com%celldm(2)*cell_com%celldm(1)-cy
          dcz=0.5_real_8*cell_com%celldm(3)*cell_com%celldm(1)-cz
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                tau0(1,ia,is)=tau0(1,ia,is)+dcx
                tau0(2,ia,is)=tau0(2,ia,is)+dcy
                tau0(3,ia,is)=tau0(3,ia,is)+dcz
             ENDDO
          ENDDO
          WRITE(6,'(/,A,A,/)')&
               ' >>>>>>>> CENTER OF MASS HAS BEEN MOVED',&
               ' TO CENTER OF BOX <<<<<<<<'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == SET WANNIER REFERENCE IF NOT SET IN INPUT                    ==
    ! ==--------------------------------------------------------------==
    IF (wannr%w_ref(1).LT.-999999.9_real_8) THEN
       ! default to middle of box if centering is active
       IF (isos1%tcent) THEN
          DO i=1,3
             wannr%w_ref(i)=0.5_real_8*(parm%a1(i)+parm%a2(i)+parm%a3(i))
          ENDDO
       ELSE
          DO i=1,3
             wannr%w_ref(i)=0.0_real_8
          ENDDO
       ENDIF

       IF (wannl%twann)&
            WRITE(6,'(A,3F11.5)')&
            ' AUTOMATIC WANNIER REF POINT [A.U.]:',(wannr%w_ref(i),i=1,3)

    ENDIF
    ! == SYMMETRIES (POINT GROUP)                                     ==
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm.AND.symmt%tgenc)CALL stopgm('SETSYS','symmetry not implemented'&
         ,& 
         __LINE__,__FILE__)
    IF (symmt%tpgauto) THEN
       WRITE(6,'(/," AUTOMATIC DETERMINATION OF THE POINT GROUP:")')
       CALL givesym(tau0)
    ENDIF
    CALL gentau(tau0)
    CALL symtau(tau0)
    CALL chksym(tau0)
    ! Construct table multiplication.
    CALL multtb(symmi%indpg,symmi%nrot,symmi%ntvec,ions1%nat,symmr%xtable,symmi%inve,symmi%multab,.TRUE.)
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm)THEN
       IF ( mmdim%natm.GT.500) GOTO 111
    ENDIF
    WRITE(6,'(/,1X,29("*"),A,28("*"))') ' ATOMS '
    WRITE(6,'(A,A)') '   NR   TYPE        X(BOHR)        Y(BOHR)',&
         '        Z(BOHR)     MBL'
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          mbl=lskcor(1,iat)+lskcor(2,iat)+lskcor(3,iat)
          WRITE(6,'(I5,5X,A2,3F15.6,I8)')&
               iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3),mbL
       ENDDO
    ENDDO
    WRITE(6,'(1X,64("*"))')
111 CONTINUE
    ! ==--------------------------------------------------------------==
    ! == CHECK IF TWO ATOMS HAVE THE SAME COORDINATES.                ==
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm)CALL mm_dim(mm_go_qm,statusdummy)
    IF (cntl%tnogeocheck) THEN
       dmin = 1.e-10_real_8
    ELSE
       dmin = 0.5_real_8
    ENDIF

    CALL dcopy(3*3,metr_com%ht,1,htsave,1)
    DO i=1,3
       metr_com%ht(1,i) = parm%a1(i)
       metr_com%ht(2,i) = parm%a2(i)
       metr_com%ht(3,i) = parm%a3(i)
    ENDDO

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO js=is,ions1%nsp
             DO ja=1,ions0%na(js)
                IF ( (is.NE.js).OR.(ia.NE.ja) ) THEN
                   dr_(1) = tau0(1,ia,is)-tau0(1,ja,js)
                   dr_(2) = tau0(2,ia,is)-tau0(2,ja,js)
                   dr_(3) = tau0(3,ia,is)-tau0(3,ja,js)
                   CALL pbc(dr_(1),dr_(2),dr_(3),dr(1),dr(2),dr(3),1,parm%apbc,parm%ibrav)
                   ddr = SQRT( dr(1)**2 + dr(2)**2 + dr(3)**2 )
                   IF (ddr .LE. dmin) THEN
                      WRITE(6,*) 'ATOM TYPE=',is,' NUM=',ia,&
                           (tau0(i,ia,is),i=1,3)
                      WRITE(6,*) 'ATOM TYPE=',js,' NUM=',ja,&
                           (tau0(i,ja,js),i=1,3)
                      CALL stopgm(' SETSYS','ATOMS ARE VERY CLOSE',& 
                           __LINE__,__FILE__)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL dcopy(3*3,htsave,1,metr_com%ht,1)
    ! ==--------------------------------------------------------------==
    ! == FILENAMES                                                    ==
    ! ==--------------------------------------------------------------==
    filn='RESTART'
    filbod='RESTART.'
    ! ==--------------------------------------------------------------==
    ! == SPECIFIED VELOCITIES                                         ==
    ! ==--------------------------------------------------------------==
    lveltrue=.FALSE.
    DO is=1,ions1%nsp
       lveltrue=lveltrue.OR.lvelini(0,is)
    ENDDO
    IF (lveltrue) THEN
       WRITE(6,'(/,1X,17("*"),A,18("*"))')&
            ' SPECIFIED ATOMIC VELOCITIES '
       WRITE(6,'(A,A)')&
            '   NR   TYPE           VX(a.u.)         VY(a.u.)',&
            '         VZ(a.u.)'
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (lvelini(ia,is))&
                  WRITE(6,'(I5,5X,A2,2X,3F17.6)')&
                  iat,elem%el(ions0%iatyp(is)),(velp(k,ia,is),k=1,3)
          ENDDO
       ENDDO
       WRITE(6,'(1X,64("*"))')
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ROKS AND LOW SPIN EXCITATION                                 ==
    ! ==--------------------------------------------------------------==
    IF (lspin2%tros.AND.lspin2%tlse) THEN
       WRITE(6,*) '                                                 '
       WRITE(6,*) '      *******************************************'
       WRITE(6,*) '      WARNING: USING ROKS IN THE &CPMD SECTION   '
       WRITE(6,*) '      AND LSE SIMULTANEOUSLY MIGHT CAUSE PROBLEMS'
       WRITE(6,*) '      ONLY DO THIS IF YOU KNOW WHAT YOU ARE DOING'
       WRITE(6,*) '      *******************************************'
    ELSE IF (lspin2%tros) THEN
       lspin2%tlse=.TRUE.
       lspin2%troks=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == SET F OCCUPATION NUMBER                                      ==
    ! ==--------------------------------------------------------------==
    ynel=crge%nel-crge%charge
    IF (cntl%tlsd) THEN
       ifac=2
    ELSE
       ifac=1
    ENDIF
    IF (crge%n.EQ.0) THEN
       ! NKSSTA is the number of unoccupied states
       ! (parameter of KOHN-SHAM ENERGIES)
       IF (cntl%tfint.AND.cnti%nkssta.EQ.0) THEN
          ! For free energy functional use 4 states + NEL/2 for
          ! the default behaviour.
          cnti%nkssta=4
       ENDIF
       nnnx=INT(ynel + 1)+ifac*cnti%nkssta
       ALLOCATE(crge%f(nnnx,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       DO i=1,nnnx
          crge%f(i,1) = -1._real_8
       ENDDO
    ENDIF
    IF (crge%f(1,1).NE.-1._real_8) THEN
       lfnotset=.FALSE.
       xnel=dasum(crge%n,crge%f(1,1),1)
       IF (ABS(ynel-xnel).GT.1.e-6_real_8) THEN
          WRITE(6,*) ' INCOMPATIBLE NUMBER OF ELECTRONS',ynel,xneL
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       lfnotset=.TRUE.
       IF (cntl%tlsd) THEN
          n1=idint(ynel)
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,n1
             crge%f(i,1)=1._real_8
          ENDDO
          IF (ABS(REAL(n1,kind=real_8)-ynel).GT.1.e-4_real_8) THEN
             WRITE(6,*)&
                  'FOR THE MOMENT USE INTEGRAL NUMBER OF ELECTRONS WITH LSD'&
                  ,n1
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ELSE
             ynel=REAL(n1,kind=real_8)
          ENDIF
          IF (crge%n.EQ.0) crge%n=n1
          IF (crge%n.LT.n1) THEN
             WRITE(6,*) ' SPECIFIED NUMBER OF STATES TOO SMALL ',crge%n,n1
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ELSE
             !$omp parallel do private(I)
#ifdef __SR11000
             !poption parallel, tlocal(I)
#endif
             DO i=n1+1,crge%n
                crge%f(i,1)=0.0_real_8
             ENDDO
          ENDIF
       ELSEIF (lspin2%tlse) THEN
          n1=idint(ynel)
          IF (ABS(REAL(n1,kind=real_8)-ynel).GT.1.e-4_real_8) THEN
             WRITE(6,*)&
                  'FOR THE MOMENT USE INTEGRAL NUMBER OF ELECTRONS WITH LSE'&
                  ,n1
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ELSE
             ynel=REAL(n1,kind=real_8)
          ENDIF
          IF (MOD(n1,2).NE.0) THEN
             WRITE(6,*) ' LSE NEEDS EVEN NUMBER OF ELECTRONS ',n1
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ENDIF
          crge%n=n1/2+1
          ! McB     ... surface hopping stuff ...
          IF ( cntl%tshop ) crge%n=crge%nel/2+1
          ! McB
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,crge%n-2
             crge%f(i,1)=2.0_real_8
          ENDDO
          crge%f(crge%n-1,1)=1.0_real_8
          crge%f(crge%n,1)=1.0_real_8
          ! Transition State Density
          IF (lspin2%tlsets) THEN
             crge%f(crge%n-1,1)=1.5_real_8
             crge%f(crge%n,1)=0.5_real_8
          ENDIF
       ELSE
          n1=idint(ynel/2._real_8)
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,n1
             crge%f(i,1)=2._real_8
          ENDDO
          IF (2._real_8*REAL(n1,kind=real_8).LT.ynel) THEN
             crge%f(n1+1,1)=ynel-2._real_8*REAL(n1,kind=real_8)
             n1=n1+1
          ENDIF
          IF (crge%n.NE.0) THEN
             IF (crge%n.LT.n1) THEN
                WRITE(6,*) ' SPECIFIED NUMBER OF STATES TOO SMALL ',crge%n,n1
                CALL stopgm('SETSYS',' ',& 
                     __LINE__,__FILE__)
             ELSE
                !$omp parallel do private(I)
#ifdef __SR11000
                !poption parallel, tlocal(I)
#endif
                DO i=n1+1,crge%n
                   crge%f(i,1)=0.0_real_8
                ENDDO
             ENDIF
          ELSE
             crge%n=n1
          ENDIF
       ENDIF
    ENDIF
    crge%nel=crge%nel-crge%charge
    ! Add NKSSTA into N
    IF (lfnotset) THEN
       !$omp parallel do private(I)
#ifdef __SR11000
       !poption parallel, tlocal(I)
#endif
       DO i=crge%n+1,crge%n+ifac*cnti%nkssta
          crge%f(i,1)=0._real_8
       ENDDO
    ENDIF
    crge%n=crge%n+ifac*cnti%nkssta
    ! .. set nkry_block if not specified in input  
    IF (cnti%nkry_block.LE.0) THEN
       IF (crge%n.GT.100) THEN
          ! NKRY_BLOCK=N/(INT(N/100)+1) 
          cnti%nkry_block=(crge%n+1)/(INT(crge%n/100)+1)
       ELSE
          cnti%nkry_block=crge%n
       ENDIF
       fixedkblock=.TRUE.
    ELSE
       fixedkblock=.FALSE.
    ENDIF
    bsfac=1
    bsclcs=1
    IF (cntl%bsymm) THEN
       ! CB: Some array sizes have to be doubled for BS
       bsfac=2
       ! CB: The multiplicity of the high spin state has to be given.
       ! CB: The spin of the low spin state is assumed to be 0.
       hspin=spin_mod%nspin
       spin_mod%nspin=0
       ! CB: Now the high spin system has to be set up.
       ! CB: Setting up of the F array has to be done at each change
       ! CB: between BS and HS state!
       nx=idint(crge%nel)-hspin+1
       IF (MOD(nx,2).NE.0) THEN
          WRITE(6,*) ' NUMBER OF ELECTRONS INCOMPATIBLE WITH ',&
               ' SPIN MULTIPLICITY'
          WRITE(6,*) ' NEL ',crge%nel,'    MULT ',hspiN
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
       hupel=REAL(nx/2,kind=real_8)
       hdoel=REAL(idint(crge%nel)-hupel,kind=real_8)
       hsup=crge%n/2-(hdoel-hupel)/2
       hsdown=crge%n-hsup
       ! CB: Setting constants for BS equations
       smin=0.0_real_8
       smax=(hspin-1.0_real_8)/2.0_real_8
       sa=0.5_real_8*smax
       sb=sa
       rtsasb=0.5_real_8/(sa*sb)
       cnstwgt=0.5_real_8*(smax-smin)*rtsasb
       ! CB: Setting general constants used in BS equations
       autocm=219474.6_real_8
    ENDIF
    ! .. 
    IF (cntl%tlsd.AND.lfnotset) THEN
       ! Spin up and down
       n1=ynel
       IF (spin_mod%nspin.EQ.0) THEN
          IF (cntl%tdiag) tdsp1%tfxsp=.FALSE.
          ! .. N is number of states
          tdsp1%nupel=n1/2
          tdsp1%ndoel=n1-tdsp1%nupel
          IF (spin_mod%nsup.EQ.0) spin_mod%nsup=crge%n/2
          spin_mod%nsup=MAX(spin_mod%nsup,INT(tdsp1%nupel))
          spin_mod%nsdown=crge%n-spin_mod%nsup
          CALL zeroing(crge%f)!,n)
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,tdsp1%nupel
             crge%f(i,1)=1._real_8
          ENDDO
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,tdsp1%ndoel
             crge%f(i+spin_mod%nsup,1)=1._real_8
          ENDDO
          spin_mod%nspin=ABS(tdsp1%nupel-tdsp1%ndoel)+1
       ELSE
          nx=idint(crge%nel)-spin_mod%nspin+1
          ! in principle we could allow fractional spin
          IF (MOD(nx,2).NE.0) THEN
             WRITE(6,*) ' NUMBER OF ELECTRONS INCOMPATIBLE WITH ',&
                  'SPIN MULTIPLICITY'
             WRITE(6,*) ' NEL ',crge%nel,'    MULT ',spin_mod%nspin
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (.NOT.cntl%tdiag) tdsp1%tfxsp=.TRUE.
          IF (cntl%tdiag.AND.cntl%ksener) tdsp1%tfxsp=.TRUE.
          ! NUPEL (NDOEL) is number of alpha (beta) electrons
          tdsp1%nupel=REAL(nx/2,kind=real_8)
          tdsp1%ndoel=REAL(idint(crge%nel)-tdsp1%nupel,kind=real_8)
          ! NSUP (NSDOWN) is number of alpha (beta) states
          IF (spin_mod%nsup.EQ.0) spin_mod%nsup=crge%n/2-(tdsp1%ndoel-tdsp1%nupel)/2
          spin_mod%nsup=MAX(spin_mod%nsup,INT(tdsp1%nupel))
          spin_mod%nsdown=crge%n-spin_mod%nsup
          CALL zeroing(crge%f)!,n)
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,INT(tdsp1%nupel)
             crge%f(i,1)=1._real_8
          ENDDO
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=spin_mod%nsup+1,spin_mod%nsup+INT(tdsp1%ndoel)
             crge%f(i,1) = 1._real_8
          ENDDO
       ENDIF
    ELSE IF (cntl%tlsd.AND..NOT.lfnotset) THEN
       IF (cnti%nkssta.NE.0) THEN
          WRITE(6,*) ' WARNING: OCCUPATION NUMBERS SET, ',&
               'NKSSTA IGNORED, NKSSTA=', cnti%nkssta
       ENDIF
       IF (spin_mod%nsup.EQ.0) spin_mod%nsup=crge%n/2
       spin_mod%nsdown=crge%n-spin_mod%nsup
       xnel=dasum(spin_mod%nsup,crge%f(1,1),1)
       tdsp1%nupel=INT(xnel)
       xnel=dasum(spin_mod%nsdown,crge%f(spin_mod%nsup+1,1),1)
       tdsp1%ndoel=INT(xnel)
       IF (tdsp1%nupel+tdsp1%ndoel.NE.INT(ynel)) THEN
          WRITE(6,*)&
               ' WARNING! FRACTIONAL ELECTRONS, GO ON AT YOUR OWN RISK'
          WRITE(6,*) ' ELECTRONS (UP/DOWN/TOTAL): ', tdsp1%nupel, tdsp1%ndoel,&
               INT(YNEL)
       ENDIF
       DO i=1,crge%n
          IF (crge%f(i,1).EQ.-1.0_real_8) crge%f(i,1) = 0.0_real_8
       ENDDO
    ENDIF
    IF (lspin2%tlsets) THEN
       IF (cntl%tshop) THEN
          WRITE(6,*) ' TRANSITION STATE ROKS NOT YET IMPLEMENTED WITH ',&
               'SURFACE HOPPING'
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
       clsd%ialpha=0
       clsd%ibeta=0
       DO i=1,crge%n
          ff=ABS(crge%f(i,1)-1.5_real_8)
          IF (ff.LT.1.e-5_real_8) clsd%ialpha=i
          IF (ff-0.5_real_8.GT.1.e-5_real_8) clsd%ibeta=i
       ENDDO
       IF (ABS(clsd%ibeta-clsd%ialpha-1.0_real_8).GT.1.e-5_real_8) THEN
          WRITE(6,*) 'HIGHER EXCITED STATES NOT PROGRAMMED'
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE IF (lspin2%tlse) THEN
       IF (cntl%tshop) THEN
          ! ..       initiate open shell ground state 
          isin=0
          nelhalf=crge%nel/2! FIXME: AK 2005/07/19  is NEL/2.eq.real(NEL,kind=real_8)/2.0?
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,nelhalf-1
             crge%f(i,1)=2.0_real_8
          ENDDO
          crge%f(nelhalf,1)=1.0_real_8
          crge%f(nelhalf+1,1)=1.0_real_8
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=nelhalf+2,crge%n
             crge%f(i,1)=0.0_real_8
          ENDDO
       ELSE
          isin=0
          DO i=1,crge%n
             IF (ABS(crge%f(i,1)-1._real_8).LT.1.e-5_real_8) isin=isin+1
          ENDDO
          IF (isin.NE.2) THEN
             WRITE(6,*) ' LSE NEEDS TWO SINGLY OCCUPIED STATES',isiN
             CALL stopgm('SETSYS',' ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ..find the singly occupied orbitals in LSE
       clsd%ialpha=0
       clsd%ibeta=0
       DO i=1,crge%n
          ff=ABS(crge%f(i,1)-1._real_8)
          IF (ff.LT.1.e-5_real_8) THEN
             IF (clsd%ialpha.EQ.0) THEN
                clsd%ialpha=i
             ELSE
                clsd%ibeta=i
             ENDIF
          ENDIF
       ENDDO
       IF (ABS(clsd%ibeta-clsd%ialpha-1.0_real_8).GT.1.e-5_real_8) THEN
          WRITE(6,*) 'HIGHER EXCITED STATES NOT PROGRAMMED'
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ..some setup for surface hopping (get all the occupations right)
    IF (cntl%tshop) THEN
       sh02%eaddsh=sh02%eaddsh/(ry*2._real_8)
       IF (cntl%tlsd) THEN
          sh02%nst_s0=crge%n
          ALLOCATE(fs0(sh02%nst_s0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          sh02%nst_s1=crge%nel/2+1
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,sh02%nst_s0
             fs0(i)=1._real_8
          ENDDO
          ! McB     ... surface hopping stuff ...
          ! McB       store at least NST_S0-NST_S1 ZEROS !
          ALLOCATE(fs1(sh02%nst_s0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(sh02%nst_s0,crge%f,1,fs1,1)
          ! McB
       ELSE
          sh02%nst_s0=crge%n-1
          ALLOCATE(fs0(sh02%nst_s0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          sh02%nst_s1=crge%n
          ALLOCATE(fs1(sh02%nst_s1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          !$omp parallel do private(I)
#ifdef __SR11000
          !poption parallel, tlocal(I)
#endif
          DO i=1,sh02%nst_s0
             fs0(i)=2._real_8
          ENDDO
          CALL dcopy(sh02%nst_s1,crge%f,1,fs1,1)
       ENDIF
       sh02%nelb2=sh02%nst_s1-1
    ENDIF
    ! mb_ssic Check the number of spin-up/spin-down electrons
    IF (cntl%tsic.AND.cntl%tlsd) THEN
       nsic=ABS(spin_mod%nsup-spin_mod%nsdown)
       WRITE(6,*)
       WRITE(6,*) ' SIC: NUMBER OF UNPAIRED ELECTRONS =',nsiC
    ENDIF
    ! mb_ssic
    ! ==--------------------------------------------------------------==
    ! == KPOINTS                                                      ==
    ! ==--------------------------------------------------------------==
    IF (restart1%rkpt) THEN
       ! Restart k points
       CALL rdkpoints
       ALLOCATE(kgemax(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kgemax)!,nkpt%nkpts)
       ALLOCATE(kcnorm(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kcnorm)!,nkpt%nkpts)
    ELSE IF (tkpts%tkpnt .OR. response1%tkpert) THEN
       IF (tkpts%tmonkp) THEN
          ! Generates k points
          nkpoint=kpts_com%nk1*kpts_com%nk2*kpts_com%nk3
          IF (tkpts%tsymkp) THEN
             ! 2^3 image per k points
             nkpoint=8*nkpoint
          ENDIF
          CALL bzmesh(nkpoint)
       ELSE
          ! In input and output k points are in cartesian coordinates.
          DO ik=1,nkpt%nkpts
             CALL dcopy(3,rk(1,ik),1,rkt(1),1)
             ! We have always inversion symmetry.
             ! By convention RK(1,IK)>=0
             IF (rkt(1).LT.0._real_8) CALL dscal(3,-1._real_8,rkt(1),1)
             CALL dcopy(3,rkt(1),1,rk(1,ik),1)
          ENDDO
          sumwk=0._real_8
          DO ik=1,nkpt%nkpts
             sumwk=sumwk+wk(ik)
          ENDDO
          IF (sumwk.EQ.0)CALL stopgm('RKPNT',&
               'THE WEIGHT OF KPOINTS ARE EQUAL TO ZERO',& 
               __LINE__,__FILE__)
          CALL dscal(nkpt%nkpts,1._real_8/sumwk,wk(1),1)
       ENDIF
       ALLOCATE(kgemax(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kgemax)!,nkpt%nkpts)
       ALLOCATE(kcnorm(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kcnorm)!,nkpt%nkpts)
    ELSE IF (.NOT.tkdp) THEN
       ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Definition of F occupation numbers for all kpoints
    IF (tkpts%tkpnt.AND.nkpt%nkpts.GT.1) THEN
       ALLOCATE(ft(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(crge%n,crge%f,1,ft(1),1)
       DEALLOCATE(crge%f,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(crge%f(crge%n,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO ik=1,nkpt%nkpts
          CALL dcopy(crge%n,ft(1),1,crge%f(1,ik),1)
       ENDDO
       DEALLOCATE(ft,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tshop) THEN
       WRITE(6,'(/,A,T60,I6,/,A,T60,I6,/,A,T54,F12.5)')&
            ' NUMBER OF STATES FOR S0 STATE:',sh02%nst_s0,&
            ' NUMBER OF STATES FOR S1 STATE:',sh02%nst_s1,&
            ' NUMBER OF ELECTRONS:',crge%nel
    ELSE
       WRITE(6,'(/,A,T60,I6,/,A,T54,F12.5)')&
            ' NUMBER OF STATES:',crge%n,&
            ' NUMBER OF ELECTRONS:',crge%nel
    ENDIF
    WRITE(6,'(A,T54,F12.5)') ' CHARGE:  ',crge%charge
    WRITE(6,'(A,T54,F12.5)')&
         ' ELECTRON TEMPERATURE(KELVIN):  ',fint1%betael
    IF (fint1%betael.GT.0._real_8) THEN
       fint1%betael=1._real_8/(fint1%betael*1.380662e-23_real_8/1.60219e-19_real_8/13.6058_real_8/2._real_8)
    ELSE
       fint1%betael=betaelmax
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tshop) THEN
       ! McB    ... surface hopping has to treat ground  a n d  excited state
       IF (cntl%tlsd) THEN
          WRITE(6,'(/2A)') ' SURFACE HOPPING: ',&
               ' GROUND STATE (UNRESTRICTED)'
          IF (spin_mod%nspin.LE.10) THEN
             IF(spin_mod%nspin==0) THEN
                WRITE(6,'(A,T47,A19)')  ' SPIN MULTIPLICITY:','NOT DEFINE IN INPUT'
             ELSE
                WRITE(6,'(A,T59,A7)')  ' SPIN MULTIPLICITY:',mul(spin_mod%nspin)
             ENDIF
          ELSE
             WRITE(6,'(A,T59,I7)')  ' SPIN MULTIPLICITY:',spin_mod%nspin
          ENDIF
          WRITE(6,'(A,T61,I5)') ' NUMBER OF ALPHA STATES:',spin_mod%nsup
          WRITE(6,'(A,T61,I5)') ' NUMBER OF BETA STATES:',spin_mod%nsdown
          WRITE(6,'(A)') ' ALPHA OCCUPATION'
          WRITE(6,'(13F5.1)') (crge%f(i,1),i=1,spin_mod%nsup)
          WRITE(6,'(A)') ' BETA OCCUPATION'
          WRITE(6,'(13F5.1)') (crge%f(i,1),i=spin_mod%nsup+1,crge%n)
       ELSE
          WRITE(6,'(/2A)') ' SURFACE HOPPING: ',&
               ' GROUND STATE (RESTRICTED)'
          WRITE(6,'(A)') ' OCCUPATION'
          WRITE(6,'(13F5.1)') (fs0(i),i=1,sh02%nst_s0)
       ENDIF
       ! McB    ... set default method for LSE ...
       lspin2%troks=.TRUE.
       WRITE(6,'(/2A)') ' SURFACE HOPPING: ',&
            ' LOW SPIN EXCITATION (RESTRICTED OPEN SHELL)'
       ! McB    ... repeat print block from below ...
       IF (lspin2%troks) THEN
          WRITE(6,'(A)') ' ROKS: I. FRANK ET AL, JCP 108 p4060 (1998)'
          WRITE(6,'(A,T50,2F8.3)') ' LSE PARAMETERS: ',lspin1%lsea,lspin1%lseb
          IF ( (lspin1%lsea.EQ.2).AND.(lspin1%lseb.EQ.-1) ) THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(1)
          ELSEIF ( (lspin1%lsea.EQ.1).AND.(lspin1%lseb.EQ.0) ) THEN
             WRITE(6,'(A,T42,A25)') ' SPIN MULTIPLICITY:',&
                  ' MIXTURE SINGLET/TRIPLET'
          ELSEIF ( (lspin1%lsea.EQ.0).AND.(lspin1%lseb.EQ.1) ) THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(3)
          ELSE
             WRITE(6,'(A)') ' CALCULATE THE MULTIPLICITY YOURSELF'
          ENDIF
       ENDIF
       WRITE(6,'(A)') ' OCCUPATION'
       WRITE(6,'(13F5.1)') (crge%f(i,1),i=1,crge%n)
    ELSEIF (cntl%tlsd) THEN
       IF (spin_mod%nspin.LE.10) THEN
          IF( spin_mod%nspin == 0 ) THEN
             WRITE(6,'(A,T47,A19)')  ' SPIN MULTIPLICITY:','NOT DEFINE IN INPUT'
          ELSE
             WRITE(6,'(A,T59,A7)')  ' SPIN MULTIPLICITY:',mul(spin_mod%nspin)
          ENDIF
       ELSE
          WRITE(6,'(A,T59,I7)')  ' SPIN MULTIPLICITY:',spin_mod%nspin
       ENDIF
       IF (cntl%bsymm) THEN
          IF (hspin.LE.10) THEN
             WRITE(6,'(A,T59,A7)')' HIGH SPIN MULTIPLICITY:',mul(hspin)
          ELSE
             WRITE(6,'(A,T59,I7)')' HIGH SPIN MULTIPLICITY:',hspiN
          ENDIF
       ENDIF
       WRITE(6,'(A,T61,I5)') ' NUMBER OF ALPHA STATES:',spin_mod%nsup
       WRITE(6,'(A,T61,I5)') ' NUMBER OF BETA STATES:',spin_mod%nsdown
       IF (cntl%bsymm) THEN
          WRITE(6,'(A,T61,I5)')' NUMBER OF HS ALPHA STATES:',hupeL
          WRITE(6,'(A,T61,I5)')' NUMBER OF HS BETA STATES:',hdoeL
       ENDIF
       WRITE(6,'(A)') ' ALPHA OCCUPATION'
       WRITE(6,'(13F5.1)') (crge%f(i,1),i=1,spin_mod%nsup)
       WRITE(6,'(A)') ' BETA OCCUPATION'
       WRITE(6,'(13F5.1)') (crge%f(i,1),i=spin_mod%nsup+1,crge%n)
    ELSEIF (lspin2%tlse) THEN
       IF (lspin2%tros) THEN
          WRITE(6,'(A)') ' RESTRICTED OPEN SHELL KOHN-SHAM (ROKS)'
       ELSE
          WRITE(6,'(A)') ' LOW SPIN EXCITATION (RESTRICTED OPEN SHELL)'
       ENDIF
       IF (lspin2%troks) THEN
          WRITE(6,'(A)') ' ROKS: I. FRANK ET AL, JCP 108 p4060 (1998)'
          WRITE(6,'(A,T50,2F8.3)') ' LSE PARAMETERS: ',lspin1%lsea,lspin1%lseb
          IF (lspin1%lsea.EQ.2.AND.lspin1%lseb.EQ.-1)THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(1)
          ELSEIF (lspin1%lsea.EQ.1.AND.lspin1%lseb.EQ.0)THEN
             WRITE(6,'(A,T42,A25)') ' SPIN MULTIPLICITY:',&
                  ' MIXTURE SINGLET/TRIPLET'
          ELSEIF (lspin1%lsea.EQ.0.AND.lspin1%lseb.EQ.1)THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(3)
          ELSE
             WRITE(6,'(A)') ' CALCULATE THE MULTIPLICITY YOURSELF'
          ENDIF
       ELSEIF (lspin2%troot) THEN
          WRITE(6,'(A)') ' ROKS: ROOTHAAN METHOD'
          WRITE(6,'(A)')' M. Filatov and S. Shaik, JCP 110 p116 (1999)'
          WRITE(6,'(A,T50,2F8.3)') ' LSE PARAMETERS: ',lspin1%lsea,lspin1%lseb
          IF (lspin1%lsea.EQ.2.AND.lspin1%lseb.EQ.-1)THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(1)
          ELSEIF (lspin1%lsea.EQ.1.AND.lspin1%lseb.EQ.0)THEN
             WRITE(6,'(A,T42,A25)') ' SPIN MULTIPLICITY:',&
                  ' MIXTURE SINGLET/TRIPLET'
          ELSEIF (lspin1%lsea.EQ.0.AND.lspin1%lseb.EQ.1)THEN
             WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(3)
          ELSE
             WRITE(6,'(A)') ' CALCULATE THE MULTIPLICITY YOURSELF'
          ENDIF
       ELSEIF (lspin2%tross) THEN
          WRITE(6,'(A)')&
               ' ROSS: J. GRAEFENSTEIN ET AL, CPL 288 p593 (1998) '
          WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(1)
       ELSEIF (lspin2%tcas22) THEN
          WRITE(6,'(A)')&
               ' CAS22: unpublished method '
          WRITE(6,'(A,T58,A8)') ' SPIN MULTIPLICITY:',mul(1)
       ENDIF
       IF (lspin2%tpenal) THEN
          WRITE(6,'(A)') ' PENALTY FUNCTION ADDED TO ENERGY'
          WRITE(6,'(A,T58,F8.3)') '         PROP. FACTOR',lspin4%apenal
       ENDIF
       IF (lspin2%tmgoe) THEN
          WRITE(6,'(A)')&
               ' MODIFIED GOEDECKER: GRIMM ET AL JCP 119 p11574 (2003) '
          WRITE(6,'(A,T50,2F8.3)') ' PARAMETERS ', lspin3%mgmab, lspin3%mgmba
       ENDIF
       IF (lspin2%tros) THEN
          WRITE(6,'(A)')&
               ' MODIFIED GOEDECKER: FRIEDRICHS ET AL'
          WRITE(6,'(A,/,A)') ' CHEM PHYS (2007) (in press)',&
               ' PARAMETERS'
          WRITE(6,'(A, F5.2,A,F5.2)')' A_AC:',lspin3%mgab(1), ' B_AC:',lspin3%mgab(2)
          WRITE(6,'(A, F5.2,A,F5.2)')' A_BC:',lspin3%mgab(3), ' B_BC:',lspin3%mgab(4)
          WRITE(6,'(A, F5.2,A,F5.2)')' A_CA:',lspin3%mgab(5), ' B_CA:',lspin3%mgab(6)
          WRITE(6,'(A, F5.2,A,F5.2)')' A_CB:',lspin3%mgab(7), ' B_CB:',lspin3%mgab(8)
          WRITE(6,'(A, F5.2,A,F5.2)')' A_BA:',lspin3%mgab(9), ' B_BA:',lspin3%mgab(10)
          WRITE(6,'(A, F5.2,A,F5.2)')' A_AB:',lspin3%mgab(11),' B_AB:',lspin3%mgab(12)
       ENDIF
       IF (lspin2%tlsets) THEN
          WRITE(6,'(A)') ' SLATER TRANSITION-STATE DENSITY WITH ROKS'
          WRITE(6,'(A)') ' BILLETER AND EGLI PRB SUBMITTED (2006)'
       ENDIF
       WRITE(6,'(A)') ' OCCUPATION'
       WRITE(6,'(13F5.1)') (crge%f(i,1),i=1,crge%n)
    ELSE
       WRITE(6,'(A)') ' OCCUPATION'
       WRITE(6,'(13F5.1)') (crge%f(i,1),i=1,crge%n)
    ENDIF
    IF(cntl%THUBB)THEN
       WRITE(*,'(/,1X,64("*"))')
       WRITE(*,'(A)')' DFT+U CORRECTIONS ARE INCLUDED'
       !IF(DFTUP(IUP_NORM_ORTHO))THEN
       IF(hubbu%pnorm)THEN
          IF(hubbu%portho) THEN
            WRITE(*,'(A,A)')' PROJECTORS ARE ORTHO-NORMALIZED',' ATOMIC ORBITALS'
          ELSE
            WRITE(*,'(A,A)')' PROJECTORS ARE NORMALIZED ',' ATOMIC ORBITALS'
          ENDIF
       ELSE
          IF(hubbu%portho) THEN
             WRITE(*,'(A,A)')' PROJECTORS ARE ORTHOGONALIZED ',' ATOMIC ORBITALS'
          ELSE
             WRITE(*,'(A,A)')' PROJECTORS ARE TRUE',' ATOMIC ORBITALS'
          ENDIF
       ENDIF
       IF(cntl%BSYMM)THEN
         WRITE(*,'(A4,2A8,2A9,A10,2A8)')'NR.','ATM.NR.','ELEMENT','U_BS(eV)','U_HS(eV)','ALPHA(eV)','SHELL','L-VALUE'
       ELSE
         WRITE(*,'(A4,2A8,A9,A10,2A8)') 'NR.','ATM.NR.','ELEMENT','U(eV)','ALPHA(eV)','SHELL','L-VALUE'
       END IF
       DO I=1,hubbu%nuatm
         DO K=1,hubbu%nl(I)
           IF(cntl%bsymm)THEN
             WRITE(*,'(I4,I8,A8,2F9.3,F10.3,2I8)')&
              I,hubbu%uatm(I),elem%el(ions0%iatyp(cpsp(hubbu%uatm(I)))),hubbu%U(I),hubbu%hs(I),&
              hubbu%a(I),hubbu%s(I,K),hubbu%l(I,K)
           ELSE
             WRITE(*,'(I4,I8,A8,F9.3,F10.3,2I8)')&
              I,hubbu%uatm(I),elem%el(ions0%iatyp(cpsp(hubbu%uatm(I)))),hubbu%u(I),hubbu%a(I),&
              hubbu%s(I,K),hubbu%l(I,K)
           END IF
         END DO
!  alph and U parameters for +U is in a.u.
!  EVTOAU = 1/(2*Ry) 
         hubbu%a(I)=hubbu%a(I)/(2*Ry)
         hubbu%u(I)=hubbu%u(I)/(2*Ry)
         hubbu%bs(I)=hubbu%u(I)
         hubbu%hs(I)=hubbu%hs(I)/(2*Ry)
       ENDDO  
       WRITE(*,'(1X,64("*"),/)')
       hubbu%firstcall=0
      ENDIF
    IF (tfixo.AND.tkpts%tkpnt.AND.nkpt%nkpts.NE.1) THEN
       WRITE(6,'(A,/,A)')&
            ' FIXED OCCUPATION NUMBERS IN COMBINATION WITH',&
            ' MORE THAN ONE K POINT NOT IMPLEMENTED'
       CALL stopgm('SETSYS','FIXED OCCUPATION & KPOINTS',& 
            __LINE__,__FILE__)
    ENDIF
    IF (tfixo.AND.cntl%tfint.AND.fint1%betael.NE.betaelmax) THEN
       WRITE(6,'(A,/,A)')&
            ' FIXED OCCUPATION NUMBERS IN COMBINATION WITH',&
            ' ELECTRONIC TEMPERATURE /= 0: INCOMPATIBLE'
       CALL stopgm('SETSYS',&
            'FIXED OCCUPATION & ELECTRON TEMPERATURE /= 0',& 
            __LINE__,__FILE__)
    ENDIF
    IF (tfixo) THEN
       WRITE(6,'(A)') ' OCCUPATION NUMBERS FIXED'
    ENDIF
    IF (cntl%tlanc.AND.fixedkblock) THEN
       ! NKRY_BLOCK is fixed in function of N.
       WRITE(6,'(/,A,/,A,T63,I3)')&
            ' LANCZOS DIAGONALIZATION (KRYLOV SUBSPACE)',&
            '    MAX. KRYLOV BLOCK SIZE   ',cnti%nkry_block
    ENDIF
    IF (tkpts%tonlydiag) THEN
       ! ONLY ONE DIAGONALIZATION, NO SELF CONSISTENCY
       WRITE(6,'(/,A,/,A,/)')&
            ' >>>>>>     ONLY ONE DIAGONALIZATION',&
            '    THE KS ENERGY ARE COMPUTED PER EACH K-POINT   '
       IF (tsphere) WRITE(6,'(/,A,/,A,A)')&
            '! !!!!!!!!!!!!  ACHTUNG ACHTUNG  !!!!!!!!!! ',&
            'TSPHERE is TRUE, BE SURE OF BUILDING CORRECT',&
            'WAVEFUNCTIONS, ESPECIALLY IF FROM RESTART'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == PSEUDOPOTENTIALS                                             ==
    ! ==--------------------------------------------------------------==
    ! IS THERE A NLCC
    IF (lqmmm%qmmm)CALL mm_dim(mm_go_qm,statusdummy)
    corel%tinlc=.FALSE.
    DO is=1,ions1%nsp
       corel%tinlc=corel%tinlc.OR.corel%tnlcc(is)
    ENDDO
    ! IS THERE A VANDERBILT PSEUDOPOTENTIAL?
    pslo_com%tivan=.FALSE.
    DO is=1,ions1%nsp
       pslo_com%tivan=pslo_com%tivan.OR.pslo_com%tvan(is)
    ENDDO
    IF (lspin2%tlse) THEN
       IF (pslo_com%tivan) THEN
          WRITE(6,*) ' VANDERBILT NOT IMPLEMENTED WITH LSE'
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
       IF (corel%tinlc) THEN
          WRITE(6,*) ' NLCC NOT IMPLEMENTED WITH LSE'
          CALL stopgm('SETSYS',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! MAXIMAL VALUE OF NGH (GAUSS-HERMIT POINT OR VDB CHANNEL)
    maxsys%nhxs=1
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nt=ncpr1%nvales(is)*ncpr1%nbeta(is)
          maxsys%nhxs=MAX(nt,maxsys%nhxs)
       ELSE
          maxsys%nhxs=MAX(nlps_com%ngh(is),maxsys%nhxs)
       ENDIF
    ENDDO
    IF (maxsys%nhxs.GT.nhx) THEN
       WRITE(6,*) ' PARAMETER NHX TOO SMALL :',maxsys%nhxs,nhX
       CALL stopgm('SETSYS',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    WRITE(6,*)
    DO is=1,ions1%nsp
       DO it=1,itmax(is)
          WRITE(6,'(72A)') vdbti(it,is)
       ENDDO
       WRITE(6,*)
    ENDDO
    WRITE(6,'(1X,64("*"))')
    WRITE(6,9900)
9900 FORMAT(1x,"*  "," ATOM ","      MASS ","  RAGGIO",&
         " NLCC ","             PSEUDOPOTENTIAL ","*")
    sgpp1%tsgpall=.TRUE.
    DO is=1,ions1%nsp
       it=1
       IF (dpot_mod%tkb(is)) it=2
       IF (pslo_com%tvan(is)) it=3
       IF (sgpp1%tsgp(is)) THEN
          IF (pslo_com%tnum(is)) THEN
             it=7
          ELSE
             it=4
          ENDIF
       ELSE
          sgpp1%tsgpall=.FALSE.
       ENDIF
       IF (dpot_mod%team(is)) it=5
       iu=2
       IF (corel%tnlcc(is)) iu=1
       IF (pslo_com%tvan(is)) THEN
          WRITE(6,9901)&
               elem%el(ions0%iatyp(is)),rmass%pma0(is),raggio(is),inlcc(iu),ischeme(it)
9901      FORMAT(1x,"* ",3x,a2,2x,f10.4,1x,f8.4,1x,a3,13x,a12,6x,"*")
       ELSEIF (it.EQ.7) THEN
          WRITE(6,9903)&
               elem%el(ions0%iatyp(is)),rmass%pma0(is),raggio(is),inlcc(iu),ischeme(it),&
               dpot_mod%lmax(is)
9903      FORMAT(1x,"* ",4x,a2,1x,f10.4,1x,f8.4,&
               1x,a3,1x,a12,3x,i3,' PROJECTORS',1x,"*")
       ELSE
          k=2
          IF (dpot_mod%lloc(is).EQ.1.OR.dpot_mod%lskip(is).EQ.1) k=1
          IF (k.EQ.1.AND.dpot_mod%lmax(is).LT.2) it=6
          WRITE(6,9902)&
               elem%el(ions0%iatyp(is)),rmass%pma0(is),raggio(is),inlcc(iu),ischeme(it),&
               cang(1),ctyp(k)
9902      FORMAT(1x,"* ",4x,a2,1x,f10.4,1x,f8.4,&
               1x,a3,1x,a12,6x,a1,a10,1x,"*")
          DO il=2,dpot_mod%lmax(is)
             k=2
             IF (dpot_mod%lloc(is).EQ.il.OR.dpot_mod%lskip(is).EQ.il) k=1
             WRITE(6,'(1X,"* ",49X,A1,A10,1X,"*")') cang(il),ctyp(k)
          ENDDO
          IF (it.EQ.1)&
               WRITE(6,'(1X,"* ",31X,A,2X,I3,1X,"*")')&
               'GH INTEGRATION POINTS:  ',INT(nlps_com%rmaxn(is))
       ENDIF
       IF ((raggio(is).GT.parm%alat/5._real_8).AND.(.NOT.ropt_mod%tesr)) THEN
          WRITE(6,'(" * ",A,F8.3,A,T65,"*")')&
               'WARNING| RAGGIO HAS TO BE LESS THAN ',parm%alat/5._real_8,&
               ' (ALAT/5)'
          WRITE(6,'(" * ",A,T65,"*")')&
               'WARNING| EWALD SUM WILL BE WRONG! '
          WRITE(6,'(" * ",A,T65,"*")')&
               'WARNING| DECREASE RAGGIO OR USE TESR'
       ENDIF
    ENDDO
    WRITE(6,'(1X,64("*"))')
    WRITE(6,*)
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    ! Initialisation for all processors
    ! Number of equivalent shell (0 for serial version)
    IF (symmi%indpg.EQ.0) THEN
       symmi%nrot=1
       symmi%inversion=0
       symmi%ntvec=0
       symmt%torigin=.FALSE.
       DO i=1,3
          symmr%origin(i)=0._real_8
       ENDDO
    ENDIF
    numshel=0
    numcold=0
    nnow=0
    ! ==--------------------------------------------------------------==
    ! cntl%tshock
    CALL mp_bcast_byte(shock1, size_in_bytes_of(shock1),parai%io_source,parai%cp_grp)
    CALL mp_bcast(maxsys%nhxs,parai%io_source,parai%cp_grp)
    CALL mp_bcast(maxsys%lpmax,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cntr%cmass,parai%io_source,parai%cp_grp)
    ! RAGGIO
    CALL mp_bcast(raggio,SIZE(raggio),parai%io_source,parai%cp_grp)
    ! RMASS
    CALL mp_bcast_byte(rmass, size_in_bytes_of(rmass),parai%io_source,parai%cp_grp)
    ! HUBBARDU
    CALL mp_bcast_byte(hubbu,size_in_bytes_of(hubbu),parai%io_source,parai%allgrp)
    ! KPOINTS
    CALL mp_bcast_byte(tkpts, size_in_bytes_of(tkpts),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nkpt%nkpts,parai%io_source,parai%cp_grp)
    IF (.NOT.paral%io_parent)  THEN
       ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(rk,SIZE(rk),parai%io_source,parai%cp_grp)
    IF (.NOT.paral%io_parent)  THEN
       ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(wk,SIZE(wk),parai%io_source,parai%cp_grp)

    IF (.NOT.paral%io_parent) THEN
       ALLOCATE(kgemax(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kgemax)!,nkpt%nkpts)
       ALLOCATE(kcnorm(nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(kcnorm)!,nkpt%nkpts)
    ENDIF

    ! ELCT (NEL,CHARGE,F occupation number)
    CALL mp_bcast(crge%n,parai%io_source,parai%cp_grp)
    CALL mp_bcast(crge%charge,parai%io_source,parai%cp_grp)
    CALL mp_bcast(crge%nel,parai%io_source,parai%cp_grp)
    IF (response1%tkpert) THEN
       IF (.NOT.paral%io_parent)  THEN
          ALLOCATE(crge%f(crge%n,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(crge%f,SIZE(crge%f),parai%io_source,parai%cp_grp)
    ELSE
       IF (.NOT.paral%io_parent)  THEN
          IF (ALLOCATED(crge%f)) THEN
             DEALLOCATE(crge%f,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
          ALLOCATE(crge%f(crge%n,nkpt%nkpts),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(crge%f,crge%n * nkpt%nkpts,parai%io_source,parai%cp_grp)
    ENDIF
    ! SPIN
    CALL mp_bcast_byte(spin_mod, size_in_bytes_of(spin_mod),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(tdsp1, size_in_bytes_of(tdsp1),parai%io_source,parai%cp_grp)
    ! IONS
    CALL mp_bcast_byte(ions0, size_in_bytes_of(ions0),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(ions1, size_in_bytes_of(ions1),parai%io_source,parai%cp_grp)
    ! TAU0
    CALL mp_bcast(tau0,SIZE(tau0),parai%io_source,parai%cp_grp)
    ! NLCC
    CALL mp_bcast_byte(corel, size_in_bytes_of(corel),parai%io_source,parai%cp_grp)
    IF (corel%tinlc) THEN
       CALL mp_bcast_byte(corer, size_in_bytes_of(corer),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(corei, size_in_bytes_of(corei),parai%io_source,parai%cp_grp)
       CALL mp_bcast(rcgrid,SIZE(rcgrid),parai%io_source,parai%cp_grp)
       CALL mp_bcast(corecg,SIZE(corecg),parai%io_source,parai%cp_grp)
    ENDIF
    ! NORM CONSERVING PP IN REAL SPACE
    CALL mp_bcast(gnl,SIZE(gnl),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rps,SIZE(rps),parai%io_source,parai%cp_grp)
    CALL mp_bcast(vr,SIZE(vr),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rw,SIZE(rw),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rv,SIZE(rv),parai%io_source,parai%cp_grp)
    ! PSLO
    CALL mp_bcast_byte(pslo_com, size_in_bytes_of(pslo_com),parai%io_source,parai%cp_grp)
    ! NKRY_BLOCK   
    CALL mp_bcast_byte(cnti, size_in_bytes_of(cnti), parai%io_source, parai%cp_grp)
    IF (pslo_com%tivan) THEN
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(rscore(maxsys%mmaxx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(dion(nbrx,nbrx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(betar(maxsys%mmaxx,nbrx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qqq(nbrx,nbrx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qfunc(maxsys%mmaxx,nbrx,nbrx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qrl(maxsys%mmaxx,nbrx,nbrx,lx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(r(maxsys%mmaxx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rucore(maxsys%mmaxx,nbrx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ru(maxsys%mmaxx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rab(maxsys%mmaxx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rsatom(maxsys%mmaxx,NSX_q),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(rscore,SIZE(rscore),parai%io_source,parai%cp_grp)
       CALL mp_bcast(dion,SIZE(dion),parai%io_source,parai%cp_grp)
       CALL mp_bcast(qqq,SIZE(qqq),parai%io_source,parai%cp_grp)
       CALL mp_bcast(betar,SIZE(betar),parai%io_source,parai%cp_grp)
       CALL mp_bcast(qfunc,SIZE(qfunc),parai%io_source,parai%cp_grp)
       CALL mp_bcast(qrl,SIZE(qrl),parai%io_source,parai%cp_grp)
       CALL mp_bcast(r,SIZE(r),parai%io_source,parai%cp_grp)
       CALL mp_bcast(rucore,SIZE(rucore),parai%io_source,parai%cp_grp)
       CALL mp_bcast(rab,SIZE(rab),parai%io_source,parai%cp_grp)
       CALL mp_bcast(rsatom,SIZE(rsatom),parai%io_source,parai%cp_grp)
       ! VDBP
       CALL mp_bcast_byte(ncpr1, size_in_bytes_of(ncpr1),parai%io_source,parai%cp_grp)
    ENDIF
    ! SGPP
    CALL mp_bcast_byte(sgpp1, size_in_bytes_of(sgpp1),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(sgpp2, size_in_bytes_of(sgpp2),parai%io_source,parai%cp_grp)
    ! DPOT
    CALL mp_bcast_byte(dpot_mod, size_in_bytes_of(dpot_mod),parai%io_source,parai%cp_grp)
    ! SYMMETRY
    CALL mp_bcast_byte(symmt, size_in_bytes_of(symmt),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(symmi, size_in_bytes_of(symmi),parai%io_source,parai%cp_grp)
    IF (symmi%indpg.NE.0) THEN
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(isymu(ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(irt(120,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (symmi%ntvec.NE.0) THEN
             ALLOCATE(tvec(3,symmi%ntvec),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(irtvec(ions1%nat,symmi%ntvec),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! XTABLE,FTABLE,FTAU,ORIGIN,DELTASYM
       CALL mp_bcast_byte(symmr, size_in_bytes_of(symmr),parai%io_source,parai%cp_grp)
       ! INVE,MULTAB
       ! ISYMU
       CALL mp_bcast(isymu,SIZE(isymu),parai%io_source,parai%cp_grp)
       ! IRT
       CALL mp_bcast(irt,SIZE(irt),parai%io_source,parai%cp_grp)
       IF (symmi%ntvec.NE.0) THEN
          ! TVEC
          CALL mp_bcast(tvec,SIZE(tvec),parai%io_source,parai%cp_grp)
          CALL mp_bcast(irtvec,SIZE(irtvec),parai%io_source,parai%cp_grp)
       ENDIF
    ELSE
       ! Number of rotations.
       symmi%nrot=1
       ! Number of translation vectors
       symmi%ntvec=0
    ENDIF
    ! FDIFF
    CALL mp_bcast(coord_fdiff,SIZE(coord_fdiff),parai%io_source,parai%cp_grp)
    ! FINT
    CALL mp_bcast_byte(fint1, size_in_bytes_of(fint1),parai%io_source,parai%cp_grp)
    ! DENSTROT
    CALL mp_bcast_byte(fint4, size_in_bytes_of(fint4),parai%io_source,parai%cp_grp)
    ! DENSBETAP
    CALL mp_bcast_byte(fint5, size_in_bytes_of(fint5),parai%io_source,parai%cp_grp)
    ! ROKS
    CALL mp_bcast_byte(lspin2, size_in_bytes_of(lspin2),parai%io_source,parai%cp_grp)
    ! TLSE orbitals
    IF (lspin2%tlse) THEN
       CALL mp_bcast(clsd%ialpha,parai%io_source,parai%cp_grp)
       CALL mp_bcast(clsd%ibeta,parai%io_source,parai%cp_grp)
    ENDIF
    ! SHOP
    IF (cntl%tshop) THEN
       CALL mp_bcast_byte(sh02, size_in_bytes_of(sh02),parai%io_source,parai%cp_grp)
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(fs0(sh02%nst_s0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(fs1(sh02%nst_s1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(fs0,SIZE(fs0),parai%io_source,parai%cp_grp)
       CALL mp_bcast(fs1,SIZE(fs1),parai%io_source,parai%cp_grp)
    ENDIF
    ! CB: cntl%bsymm
    CALL mp_bcast(bsfac,parai%io_source,parai%cp_grp)
    CALL mp_bcast(bsclcs,parai%io_source,parai%cp_grp)
    ! WANNIER
    CALL mp_bcast(wannr%w_ref,3,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! NLM is used by CALC_ALM, EHPSI, SUMFNL, RNLSM1, RNLSM2.
    nlm=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nlm=nlm+ions0%na(is)*ncpr1%nvales(is)*ncpr1%nbeta(is)
       ELSE
          nlm=nlm+ions0%na(is)*nlps_com%ngh(is)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    ! define an odd leading dimension for EI1 etc
    natx = 2*(ions1%nat/2) + 1
    IF (natx.LT.ions1%nat) CALL stopgm("SETSYS","NATX.LT.NAT",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! kpt
    IF (.NOT.tkpts%tkpnt) THEN
       ALLOCATE(nkpbl(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nkpbl(1)=1
       ALLOCATE(kpbeg(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       kpbeg(1)=0
       kpts_com%nkptall=1
    ENDIF
    IF (tkpts%tkpnt.AND.pslo_com%tivan) CALL stopgm('SETSYS',&
         'K-POINTS WITH VANDERBILT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (tkpts%tkpnt.AND.cntl%tlowd) CALL stopgm('SETSYS',&
         'K-POINTS WITH LOWDIN ORTHO NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! kpt
    ! ==--------------------------------------------------------------==
    CALL mp_bcast_byte(clsaabox, size_in_bytes_of(clsaabox), parai%io_source, parai%cp_grp)
    ! ==--------------------------------------------------------------==
    IF (lqmmm%qmmm)CALL mm_dim(mm_go_qm,status) ! because it is needed
    CALL cpl_para
    RETURN
  END SUBROUTINE setsys
  ! ==================================================================

END MODULE setsys_utils
