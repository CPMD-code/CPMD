MODULE setsc_utils
  USE bc,                              ONLY: bc_com
  USE cell,                            ONLY: cell_com,&
                                             lcell
  USE clas,                            ONLY: clas8
  USE cnst,                            ONLY: fbohr,&
                                             pi
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: gvec_com
  USE isos,                            ONLY: isos1,&
                                             isos2
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: genlat,&
                                             latgen,&
                                             omegagen
  USE metr,                            ONLY: metr_com
  USE molsym_utils,                    ONLY: molsym
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE sort_utils,                      ONLY: sort2
  USE symm,                            ONLY: stable,&
                                             stag,&
                                             symmi,&
                                             symmr,&
                                             symmt
  USE symmetry_utils,                  ONLY: alfasym,&
                                             bccsym,&
                                             bctsym,&
                                             fccsym,&
                                             pgsym,&
                                             trgsym
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             dual00,&
                                             parm
  USE utils,                           ONLY: invmat
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setsc
  PUBLIC :: ihmat

CONTAINS

  ! ==================================================================
  SUBROUTINE setsc
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    ! Variables
    REAL(real_8), PARAMETER                  :: deltakin = 1.e-10_real_8 

    INTEGER                                  :: i, i1, i2, indx(4), kk
    REAL(real_8)                             :: bbbb, cccc, eedif, fact, &
                                                sinn, xl(4)

    IF (.NOT.paral%io_parent) GOTO 9999
    ! ==--------------------------------------------------------------==
    ! ==  SETUP DIMENSIONS OF THE SUPERCELL FOR APPLYING PBC          ==
    ! ==--------------------------------------------------------------==
    IF (parm%ibrav.EQ.-2 .OR. lcell%tcellvectors) THEN
       ! LATTICE VECTORS GIVEN
       IF (.NOT.cntl%bohr) THEN
          DO i=1,3
             parm%a1(i) = parm%a1(i)*fbohr
             parm%a2(i) = parm%a2(i)*fbohr
             parm%a3(i) = parm%a3(i)*fbohr
          ENDDO
       ENDIF
       CALL genlat(parm%a1,parm%a2,parm%a3,cell_com%celldm,parm%omega)
    ELSE
       IF (.NOT.cntl%bohr) cell_com%celldm(1)=cell_com%celldm(1)*fbohr
       IF (.NOT.cntl%bohr) prcp_com%cellrf(1)=prcp_com%cellrf(1)*fbohr
       IF (.NOT.cntl%bohr) isos2%rcerf=isos2%rcerf*fbohr
    ENDIF

    ! ==--------------------------------------------------------------==
    IF (lcell%tcellrefvec) THEN
       ! LATTICE VECTORS GIVEN
       IF (.NOT.cntl%bohr) THEN
          DO i=1,3
             cell_com%ar1(i) = cell_com%ar1(i)*fbohr
             cell_com%ar2(i) = cell_com%ar2(i)*fbohr
             cell_com%ar3(i) = cell_com%ar3(i)*fbohr
          ENDDO
       ENDIF
       CALL genlat(cell_com%ar1,cell_com%ar2,cell_com%ar3,prcp_com%cellrf,parm%omega)
    ELSEIF (prcp_com%cellrf(1).LE.0.0_real_8) THEN
       CALL dcopy(6,cell_com%celldm(1),1,prcp_com%cellrf(1),1)
       IF (lcell%tcellvectors) THEN
          lcell%tcellrefvec = .TRUE.
          DO i=1,3
             cell_com%ar1(i) = parm%a1(i)
             cell_com%ar2(i) = parm%a2(i)
             cell_com%ar3(i) = parm%a3(i)
          ENDDO
       ENDIF
    ENDIF
    IF (cntl%tprcp.AND.(.NOT.prcpl%tisot))THEN
       ! Triclinic cell
       parm%alat=prcp_com%cellrf(1)
    ELSEIF (parm%ibrav.EQ.-2 .OR. lcell%tcellvectors) THEN
       ! Triclinic cell
       parm%alat=prcp_com%cellrf(1)
    ELSE
       sinn=SQRT(1._real_8-cell_com%celldm(4)**2)
       parm%alat=prcp_com%cellrf(1)
       IF (parm%ibrav.EQ.1) THEN
          ! Simple cubic
          parm%apbc(1)=parm%alat/2._real_8
       ELSEIF (parm%ibrav.EQ.2) THEN
          ! Face centred cubic
          parm%apbc(1)=parm%alat*SQRT(2._real_8)/4._real_8
       ELSEIF (parm%ibrav.EQ.4) THEN
          ! Hexagonal
          parm%apbc(1)=parm%alat/2._real_8
          parm%apbc(2)=parm%alat*SQRT(3._real_8)/4._real_8
          parm%apbc(3)=parm%alat/2._real_8*prcp_com%cellrf(3)
          parm%apbc(4)=SQRT(3._real_8)/3._real_8
       ELSEIF (parm%ibrav.EQ.8) THEN
          ! Orthorhombic
          parm%apbc(1)=parm%alat/2._real_8
          parm%apbc(2)=parm%alat/2._real_8*prcp_com%cellrf(2)
          parm%apbc(3)=parm%alat/2._real_8*prcp_com%cellrf(3)
       ELSEIF (parm%ibrav.EQ.12) THEN
          ! JOK      NEW GEOMETRY PROGRAMMED (JOK)
          bbbb=parm%alat*prcp_com%cellrf(2)*prcp_com%cellrf(4)
          cccc=parm%alat*prcp_com%cellrf(2)*sinn
          ! LOOK FOR THE CLOSEST IMAGES TO BUILD THE WIGNER-SEITZ CELL.
          DO kk=1,4
             xl(kk)=(parm%alat-REAL(kk-1,kind=real_8)*bbbb)**2+(REAL(kk-1,kind=real_8)*cccc)**2
          ENDDO
          ! ORDER THE IMAGES WITH INCREASING DISTANCES
          CALL sort2(xl,4,indx)
          ! HERE WE DEFINE THE DIMENSIONS OF THE W-S CELL. APBC(1),APBC(2)
          ! AND APBC(4) ARE IN THE (X,Y) PLANE, WHERE THE CELL IS A NON-
          ! REGULAR HEXAGON; WHILE APBC(3) CORRESPONDS TO THE Z-DIRECTION.
          ! THE HEXAGON IS DEFINED BY THE 3 CLOSEST IMAGES IN THE (X,Y)
          ! PLANE. TWO OF THEM WERE CALCULATED BEFORE, AND THE 3RD (THE
          ! SMALLEST ONE) IS THAT IN THE DIRECTION OF THE LATTICE VECTOR
          ! NUMBER 2 (CELLDM(2)).
          parm%apbc(1)=SQRT(xl(1))/2._real_8
          parm%apbc(2)=parm%alat/2._real_8*prcp_com%cellrf(2)
          parm%apbc(3)=parm%alat/2._real_8*prcp_com%cellrf(3)
          parm%apbc(4)=SQRT(xl(2))/2._real_8
       ELSEIF (parm%ibrav.EQ.14) THEN
          ! Special case of triclinic: monoclinic with gamma<>90 degrees
          IF (ABS(prcp_com%cellrf(4)).LT.1.e-10_real_8.AND.ABS(prcp_com%cellrf(5)).LT.1.e-10_real_8) &
               THEN
             parm%apbc(1)=parm%alat/2._real_8
             parm%apbc(2)=parm%alat/2._real_8*prcp_com%cellrf(2)*SQRT(1._real_8-prcp_com%cellrf(6)**2)
             parm%apbc(3)=parm%alat/2._real_8*prcp_com%cellrf(3)
             parm%apbc(4)=-prcp_com%cellrf(6)/SQRT(1._real_8-prcp_com%cellrf(6)**2)
          ELSE
             parm%apbc(2)=0.0_real_8
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  CONSTANTS RELATED TO THE DIMENSIONS OF THE RECIPROCAL CELL  ==
    ! ==--------------------------------------------------------------==
    parm%tpiba=2._real_8*pi/parm%alat
    parm%tpiba2=parm%tpiba*parm%tpiba
    gvec_com%gcutw=cntr%ecut/parm%tpiba2
    IF (prcp_com%akin.GT.deltakin) THEN
       prcp_com%gckin=prcp_com%eckin/parm%tpiba2
       prcp_com%gakin=prcp_com%akin/parm%tpiba2
       prcp_com%gskin=prcp_com%skin/parm%tpiba2
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DUAL APPROXIMATION  ***  BE CAREFUL ***                     ==
    ! ==--------------------------------------------------------------==
    IF (dual00%dual) THEN
       IF (dual00%cdual.LT.1._real_8) dual00%cdual=1._real_8
    ELSE
       dual00%cdual=4.0_real_8
    ENDIF
    gvec_com%gcut=gvec_com%gcutw*dual00%cdual
    ! ==--------------------------------------------------------------==
    ! ==  METRIC TENSOR OF SUPERCELL                                  ==
    ! ==--------------------------------------------------------------==
    IF (lcell%tcellvectors) THEN
       DO i=1,3
          metr_com%ht(1,i) = parm%a1(i)
          metr_com%ht(2,i) = parm%a2(i)
          metr_com%ht(3,i) = parm%a3(i)
       ENDDO
    ELSEIF (parm%ibrav.GE.0) THEN
       IF (cntl%tprcp) THEN
          CALL latgen(parm%ibrav,cell_com%celldm,parm%a1,parm%a2,parm%a3,parm%omega)
       ELSE
          CALL latgen(parm%ibrav,cell_com%celldm,parm%a1,parm%a2,parm%a3,parm%omega)
          CALL cry(parm%ibrav,parm%a1,parm%a2,parm%a3,indx)
       ENDIF
       DO i=1,3
          metr_com%ht(1,i) = parm%a1(i)
          metr_com%ht(2,i) = parm%a2(i)
          metr_com%ht(3,i) = parm%a3(i)
       ENDDO
    ENDIF
    CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
    CALL dcopy(9,  metr_com%ht(1,1),1,  metr_com%ht0(1,1),1)
    CALL dcopy(9,metr_com%htm1(1,1),1,metr_com%htm10(1,1),1)

    DO i=1,3
       fact=1.0_real_8/SQRT(metr_com%ht(i,1)**2+metr_com%ht(i,2)**2+metr_com%ht(i,3)**2)
       prcp_com%hunit(i,1)=fact*metr_com%ht(i,1)
       prcp_com%hunit(i,2)=fact*metr_com%ht(i,2)
       prcp_com%hunit(i,3)=fact*metr_com%ht(i,3)
    ENDDO

    prcp_com%omega0=parm%omega
    ! ==--------------------------------------------------------------==
    ! ==  GENERATE LATTICE VECTORS                                    ==
    ! ==--------------------------------------------------------------==
    IF (lcell%tcellrefvec) THEN
       DO i = 1,3
          parm%a1(i)=cell_com%ar1(i)
          parm%a2(i)=cell_com%ar2(i)
          parm%a3(i)=cell_com%ar3(i)
       ENDDO
       CALL omegagen(parm%a1,parm%a2,parm%a3,prcp_com%omega0)
    ELSEIF (parm%ibrav.GE.0) THEN
       IF (cntl%tprcp) THEN
          CALL latgen(parm%ibrav,prcp_com%cellrf,parm%a1,parm%a2,parm%a3,prcp_com%omega0)
       ELSE
          CALL latgen(parm%ibrav,prcp_com%cellrf,parm%a1,parm%a2,parm%a3,prcp_com%omega0)
          CALL cry(parm%ibrav,parm%a1,parm%a2,parm%a3,indx)
       ENDIF
    ENDIF
    IF (parm%ibrav.EQ.-2) THEN
       ! LATTICE VECTORS ARE GIVEN
       ! WE USE TRICLINIC SYSTEM
       parm%ibrav=14
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==  SET POINT GROUP TRANSFORMATION MATRICES                     ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(stable)!,9*48)
    CALL zeroing(symmr%xtable)!,9*120)
    IF (symmi%indpg.NE.0.AND..NOT.symmt%tmsym) THEN
       IF (parm%ibrav.EQ.1) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.2) THEN
          CALL fccsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.3) THEN
          CALL bccsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.4) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.5) THEN
          CALL trgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.6) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.7) THEN
          CALL bctsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.8) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.10) THEN
          CALL alfasym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.12) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ELSEIF (parm%ibrav.EQ.14) THEN
          CALL pgsym(symmi%indpg,stable,symmi%nrot)
       ENDIF
       DO i=1,symmi%nrot
          DO i1=1,3
             DO i2=1,3
                ! (T. Deutsch 29/1/98 change to transpose)
                ! XTABLE(I1,I2,I)=real(STABLE(I1,I2,I),kind=real_8)
                symmr%xtable(i1,i2,i)=REAL(stable(i2,i1,i),kind=real_8)
             ENDDO
          ENDDO
       ENDDO
    ELSEIF (symmi%indpg.NE.0.AND.symmt%tmsym) THEN
       CALL molsym(stag,symmi%naxis,symmr%xtable,symmi%nrot)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  CLASSICAL BOX                                               ==
    ! ==--------------------------------------------------------------==
    IF (clas8%cellcl(1).EQ.0.0_real_8) THEN
       CALL dcopy(6,cell_com%celldm(1),1,clas8%cellcl(1),1)
    ENDIF
    eedif=0._real_8
    DO i=1,6
       eedif=eedif+(cell_com%celldm(i)-clas8%cellcl(i))**2
    ENDDO
    IF (eedif.GT.1.e-12_real_8) THEN
       IF (.NOT.isos1%tclust) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' ISOLATED SYSTEM OPTION REQUIRED FOR '
          IF (paral%io_parent)&
               WRITE(6,*) ' NONEQUAL CLASSICAL AND QUANTUM CELLS'
          CALL stopgm('SETSC','SPECIFY IBRAV=0',& 
               __LINE__,__FILE__)
       ELSEIF (clas8%cellcl(1).LT.cell_com%celldm(1)) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' CLASSICAL CELL < QUANTUM CELL '
          CALL stopgm('SETSC','CL CELL < QM CELL',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    CALL latgen(parm%ibrav,clas8%cellcl,clas8%cla1,clas8%cla2,clas8%cla3,clas8%clomega)
    DO i=1,3
       clas8%clht(1,i) = clas8%cla1(i)
       clas8%clht(2,i) = clas8%cla2(i)
       clas8%clht(3,i) = clas8%cla3(i)
    ENDDO
    CALL ihmat(clas8%clht,clas8%clhtm1,clas8%clomega)
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    ! ..PARM
    CALL mp_bcast_byte(parm, size_in_bytes_of(parm),parai%io_source,parai%cp_grp)
    ! ..PRCP
    CALL mp_bcast_byte(prcp_com, size_in_bytes_of(prcp_com),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(prcpl, size_in_bytes_of(prcpl),parai%io_source,parai%cp_grp)
    ! ..METR
    CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
    ! ..BC
    CALL mp_bcast_byte(bc_com, size_in_bytes_of(bc_com),parai%io_source,parai%cp_grp)
    ! ..GVEC
    CALL mp_bcast_byte(gvec_com, size_in_bytes_of(gvec_com),parai%io_source,parai%cp_grp)
    ! ..CELLDM
    CALL mp_bcast_byte(cell_com, size_in_bytes_of(cell_com),parai%io_source,parai%cp_grp)
    ! ..TCELLVECTORS
    CALL mp_bcast_byte(lcell, size_in_bytes_of(lcell),parai%io_source,parai%cp_grp)
    ! ..CELLCL
    CALL mp_bcast_byte(clas8, size_in_bytes_of(clas8),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setsc
  ! ==================================================================
  SUBROUTINE ihmat(ht,htm1,omega)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ht(3,3), htm1(3,3), omega

    INTEGER                                  :: i, info, iperm, j, k, l
    REAL(real_8)                             :: aux(1000), s

    CALL dcopy(9,ht(1,1),1,htm1(1,1),1)
    CALL invmat(3,htm1,aux,info)
    ! Compute the volume of the Supercell (OMEGA) 
    omega=0._real_8
    s=1._real_8
    i=1
    j=2
    k=3
101 CONTINUE
    DO iperm=1,3
       omega=omega+s*ht(1,i)*ht(2,j)*ht(3,k)
       l=i
       i=j
       j=k
       k=l
    ENDDO
    i=2
    j=1
    k=3
    s=-s
    IF (s.LT.0._real_8) go to 101
    omega=ABS(omega)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ihmat
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cry(ibrav,a1,a2,a3,indx)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE CALCULATES THE VECTORS V1CRYS, V2CRYS AND V3CRYS, WHICH ARE ==
    ! == THEN USED BY ROUTINE PBC.                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ibrav
    REAL(real_8)                             :: a1(3), a2(3), a3(3)
    INTEGER                                  :: indx(4)

    INTEGER                                  :: i
    REAL(real_8)                             :: aux

    IF (ibrav.EQ.12) THEN
       bc_com%v1crys(1)=a1(1)-REAL(indx(1)-1,kind=real_8)*a2(1)
       bc_com%v1crys(2)=-REAL(indx(1)-1,kind=real_8)*a2(2)
       aux=SQRT(bc_com%v1crys(1)**2+bc_com%v1crys(2)**2)
       bc_com%v1crys(3)=bc_com%v1crys(1)/aux
       bc_com%v1crys(4)=bc_com%v1crys(2)/aux
       bc_com%v2crys(1)=a2(1)
       bc_com%v2crys(2)=a2(2)
       aux=SQRT(bc_com%v2crys(1)**2+bc_com%v2crys(2)**2)
       bc_com%v2crys(3)=bc_com%v2crys(1)/aux
       bc_com%v2crys(4)=bc_com%v2crys(2)/aux
       bc_com%v3crys(1)=a1(1)-REAL(indx(2)-1,kind=real_8)*a2(1)
       bc_com%v3crys(2)=-REAL(indx(2)-1,kind=real_8)*a2(2)
       aux=SQRT(bc_com%v3crys(1)**2+bc_com%v3crys(2)**2)
       bc_com%v3crys(3)=bc_com%v3crys(1)/aux
       bc_com%v3crys(4)=bc_com%v3crys(2)/aux
    ELSE
       DO i=1,4
          bc_com%v1crys(i)=0._real_8
          bc_com%v2crys(i)=0._real_8
          bc_com%v3crys(i)=0._real_8
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cry
  ! ==================================================================

END MODULE setsc_utils
