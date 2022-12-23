MODULE clas_force_utils
  USE clas,                            ONLY: &
       clas3, clas8, clasc, clasc0, clasf, claspres, cltyp, myclat1, myclat2, &
       pote1, pr12, rcr12, rlrc
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai
  USE pbc_utils,                       ONLY: pbc_clas
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: icopy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: clas_force
  PUBLIC :: check
  !public :: save0
  !public :: cell_listr12
  !public :: check_cell

CONTAINS

  ! ==================================================================
  ! mdebug      SUBROUTINE CLAS_FORCE(ENERGY,UPDATE)
  SUBROUTINE clas_force(energy,update,eint)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: energy
    LOGICAL                                  :: update
    REAL(real_8)                             :: eint

    CHARACTER(*), PARAMETER                  :: procedureN = 'clas_force'
    INTEGER, PARAMETER                       :: maxnab = 100 

    INTEGER                                  :: ia1, ia2, ierr, is1, is2, &
                                                isub, jbeg, jend, jnab, nlist
    INTEGER, ALLOCATABLE                     :: list2(:)
    INTEGER, ALLOCATABLE, SAVE               :: cellr12(:,:), list(:), &
                                                point(:)
    INTEGER, DIMENSION(3), SAVE              :: cdimr12
    INTEGER, SAVE                            :: ifirst = 0, listlen
    REAL(real_8) :: fx, fy, fz, rcut, rcutsq, rij, rijsq, rlistsq, rxia1, &
      ryia1, rzia1, sigmasq, sr12, sr2, src12, vij, wij, x1, x2, y1, y2, z1, &
      z2

    CALL tiset('CLAS_FORCE',isub)
    IF (ifirst.EQ.0) THEN
       ALLOCATE(point(clas3%nclatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       listlen=maxnab*(myclat2-myclat1+1)
       ALLOCATE(list(listlen),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(point)!,clas3%nclatom)
       CALL zeroing(list)!,listlen)
       IF (pote1%t_r12) THEN
          ALLOCATE(cellr12(3,(3*clas3%nclatom)/3),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(cellr12)!,3*clas3%nclatom)
       ENDIF
       ifirst=1
    ENDIF

    ! PARAMETER RLRC is the difference R_list - R_cutoff
    ! it has the same value for all pairs of particle types
    rlrc=2.0_real_8

    energy = 0._real_8
    eint=0._real_8
    CALL zeroing(claspres)!,9)
    CALL zeroing(clasf)!,3*clas3%nclatom)
    ! ==--------------------------------------------------------------==

    IF (pote1%t_r12) THEN

       IF (update) THEN
          ! Verlet list has to be updated

          CALL save0
          nlist = 0
          CALL cell_listr12(cdimr12,cellr12)

          DO ia1=myclat1,myclat2
             is1=cltyp(ia1)

             point(ia1)=nlist+1

             rxia1 = clasc(1,ia1)
             ryia1 = clasc(2,ia1)
             rzia1 = clasc(3,ia1)

             DO ia2 = ia1+1,clas3%nclatom
                is2=cltyp(ia2)

                IF (clas3%is_qm(is1)*clas3%is_qm(is2).EQ.0) THEN
                   IF (check_cell(cellr12(1,ia1),cellr12(1,ia2),cdimr12))&
                        THEN
                      sigmasq = pr12(is1,is2)*pr12(is1,is2)
                      rcutsq = rcr12(is1,is2)**2
                      rcut=rcr12(is1,is2)
                      rlistsq = (rcr12(is1,is2)+rlrc)**2

                      x1 = rxia1 - clasc(1,ia2)
                      y1 = ryia1 - clasc(2,ia2)
                      z1 = rzia1 - clasc(3,ia2)

                      CALL pbc_clas(x1,y1,z1,x2,y2,z2)

                      rijsq = x2**2 + y2**2 + z2**2
                      rij=SQRT(rijsq)

                      IF (rijsq.LT.rlistsq) THEN
                         nlist=nlist+1
                         IF (nlist.EQ.listlen) THEN
                            ALLOCATE(list2(listlen),STAT=ierr)
                            IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                                 __LINE__,__FILE__)
                            CALL icopy(nlist,list,1,list2,1)
                            DEALLOCATE(list,STAT=ierr)
                            IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                                 __LINE__,__FILE__)
                            listlen=3*listlen/2
                            ALLOCATE(list(listlen),STAT=ierr)
                            IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                                 __LINE__,__FILE__)
                            CALL icopy(nlist,list2,1,list,1)
                            DEALLOCATE(list,STAT=ierr)
                            IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                                 __LINE__,__FILE__)
                         ENDIF
                         list(nlist)=ia2
                         IF (rijsq.LT.rcutsq) THEN

                            sr2   = sigmasq/rijsq
                            sr12  = sr2**6
                            src12= (sigmasq/rcutsq)**6
                            ! m                           VIJ   = SR12
                            vij   = sr12 +src12*(12._real_8*rij-13._real_8*rcut)/rcut
                            wij   = 12.0_real_8 * sr12 / rijsq-12._real_8*src12/rij/&
                                 rcut
                            fx = wij * x2
                            fy = wij * y2
                            fz = wij * z2

                            clasf(1,ia1) = clasf(1,ia1) + fx
                            clasf(2,ia1) = clasf(2,ia1) + fy
                            clasf(3,ia1) = clasf(3,ia1) + fz

                            clasf(1,ia2) = clasf(1,ia2) - fx
                            clasf(2,ia2) = clasf(2,ia2) - fy
                            clasf(3,ia2) = clasf(3,ia2) - fz

                            energy = energy + vij
                            IF ((clas3%is_qm(is1)+clas3%is_qm(is2)).NE.0)eint=eint+vij

                            claspres(1,1) = claspres(1,1) + fx * x2
                            claspres(2,2) = claspres(2,2) + fy * y2
                            claspres(3,3) = claspres(3,3) + fz * z2
                            claspres(1,2) = claspres(1,2) + fx * y2
                            claspres(1,3) = claspres(1,3) + fx * z2
                            claspres(2,3) = claspres(2,3) + fy * z2
                            claspres(2,1) = claspres(1,2)
                            claspres(3,1) = claspres(1,3)
                            claspres(3,2) = claspres(2,3)

                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
          IF (myclat2.LT.clas3%nclatom) point(myclat2+1)=nlist+1

       ELSE
          ! we use the list to find the neighbours

          DO ia1=myclat1,myclat2
             IF (ia1.LT.clas3%nclatom) THEN
                is1=cltyp(ia1)

                jbeg=point(ia1)
                jend=point(ia1+1)-1

                IF (jbeg.LE.jend) THEN

                   rxia1 = clasc(1,ia1)
                   ryia1 = clasc(2,ia1)
                   rzia1 = clasc(3,ia1)

                   DO jnab=jbeg,jend
                      ia2=list(jnab)
                      is2=cltyp(ia2)

                      sigmasq = pr12(is1,is2)*pr12(is1,is2)
                      rcutsq = rcr12(is1,is2)**2
                      rcut=rcr12(is1,is2)
                      x1 = rxia1 - clasc(1,ia2)
                      y1 = ryia1 - clasc(2,ia2)
                      z1 = rzia1 - clasc(3,ia2)

                      CALL pbc_clas(x1,y1,z1,x2,y2,z2)

                      rijsq = x2**2 + y2**2 + z2**2
                      rij=SQRT(rijsq)

                      IF (rijsq.LT.rcutsq) THEN

                         sr2   = sigmasq/rijsq
                         sr12  = sr2**6
                         src12 = (sigmasq/rcutsq)**6
                         vij   = sr12+src12*(12._real_8*rij-13._real_8*rcut)/rcut
                         wij   = 12.0_real_8 * sr12 / rijsq-12._real_8*src12/rij/&
                              rcut
                         ! VIJ   = SR12
                         ! WIJ   = 12.0 * SR12 / RIJSQ
                         ! 
                         fx = wij * x2
                         fy = wij * y2
                         fz = wij * z2

                         clasf(1,ia1) = clasf(1,ia1) + fx
                         clasf(2,ia1) = clasf(2,ia1) + fy
                         clasf(3,ia1) = clasf(3,ia1) + fz

                         clasf(1,ia2) = clasf(1,ia2) - fx
                         clasf(2,ia2) = clasf(2,ia2) - fy
                         clasf(3,ia2) = clasf(3,ia2) - fz

                         energy = energy + vij
                         IF ( (clas3%is_qm(is1)*clas3%is_qm(is2).EQ.0)  .AND.    ( (clas3%is_qm&
                              (is1)+clas3%is_qm(is2)).NE.0 ) )eint=eint+vij

                         claspres(1,1) = claspres(1,1) + fx * x2
                         claspres(2,2) = claspres(2,2) + fy * y2
                         claspres(3,3) = claspres(3,3) + fz * z2
                         claspres(1,2) = claspres(1,2) + fx * y2
                         claspres(1,3) = claspres(1,3) + fx * z2
                         claspres(2,3) = claspres(2,3) + fy * z2
                         claspres(2,1) = claspres(1,2)
                         claspres(3,1) = claspres(1,3)
                         claspres(3,2) = claspres(2,3)

                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
          ENDDO

       ENDIF

    ELSE
       CALL stopgm('CLAS_FORCE','UNDEFINED FORCE FIELD',& 
            __LINE__,__FILE__)
    ENDIF
    CALL mp_sync(parai%allgrp)
    CALL mp_sum(energy,parai%allgrp)
    CALL mp_sum(eint,parai%allgrp)
    CALL mp_sum(claspres,9,parai%allgrp)
    CALL mp_sum(clasf,3*clas3%nclatom,parai%allgrp)
    CALL tihalt('CLAS_FORCE',isub)
    RETURN
  END SUBROUTINE clas_force


  SUBROUTINE check ( update )

    ! *******************************************************************
    ! ** DECIDES WHETHER THE LIST NEEDS TO BE RECONSTRUCTED.           **
    ! **                                                               **
    ! ** PRINCIPAL VARIABLES:                                          **
    ! **                                                               **
    ! ** CLASC(J,IA)                    ATOM POSITIONS                 **
    ! ** CLASC0(J,IA)                   COORDINATES AT LAST UPDATE     **
    ! ** REAL     RLIST                 RADIUS OF VERLET LIST          **
    ! ** REAL     RCUT                  CUTOFF DISTANCE FOR FORCES     **
    ! ** REAL     DISPMX                LARGEST DISPLACEMENT           **
    ! ** LOGICAL  UPDATE                IF TRUE THE LIST IS UPDATED    **
    ! **                                                               **
    ! ** USAGE:                                                        **
    ! **                                                               **
    ! ** CHECK IS CALLED TO SET UPDATE BEFORE EVERY CALL TO FORCE.     **
    ! *******************************************************************

    LOGICAL                                  :: update

    INTEGER                                  :: ia
    REAL(real_8)                             :: dispmx

! *******************************************************************
! ** CALCULATE MAXIMUM DISPLACEMENT SINCE LAST UPDATE **

    dispmx = 0.0_real_8

    DO ia=myclat1,myclat2
       dispmx = MAX ( ABS ( clasc0(1,ia)-clasc(1,ia) ), dispmx )
       dispmx = MAX ( ABS ( clasc0(2,ia)-clasc(2,ia) ), dispmx )
       dispmx = MAX ( ABS ( clasc0(3,ia)-clasc(3,ia) ), dispmx )
    ENDDO
    CALL mp_max(dispmx,parai%allgrp)
    ! ** A CONSERVATIVE TEST OF THE LIST SKIN CROSSING **
    dispmx = 2.0_real_8 * SQRT  ( 3.0_real_8 * dispmx * dispmx )

    update = dispmx.GT.rlrc

    RETURN
  END SUBROUTINE check


  SUBROUTINE save0


    ! *******************************************************************
    ! ** SAVES CURRENT CONFIGURATION FOR FUTURE CHECKING.              **
    ! **                                                               **
    ! ** PRINCIPAL VARIABLES:                                          **
    ! **                                                               **
    ! ** CLASC(J,IA)                    ATOM POSITIONS                 **
    ! ** CLASC0(J,IA)                   STORAGE LOCATIONS FOR SAVE     **
    ! **                                                               **
    ! ** SAVE IS CALLED WHENEVER THE NEW VERLET LIST IS CONSTRUCTED.   **
    ! *******************************************************************

    INTEGER                                  :: ia

! *******************************************************************

    DO ia=1,clas3%nclatom
       clasc0(1,ia)=clasc(1,ia)
       clasc0(2,ia)=clasc(2,ia)
       clasc0(3,ia)=clasc(3,ia)
    ENDDO

    RETURN
  END SUBROUTINE save0

  ! *******************************************************************
  SUBROUTINE cell_listr12(cdimr12,cellr12)

    INTEGER                                  :: cdimr12(3), cellr12(3,*)

    INTEGER                                  :: ia, imx, imy, imz, is1, is2
    REAL(real_8)                             :: cellm, dx, dy, dz, xdim, xp, &
                                                ydim, yp, zdim, zp

    cellm =0._real_8
    DO is1=1,clas3%ncltyp
       DO is2=is1,clas3%ncltyp
          cellm = MAX(cellm,(rcr12(is1,is2)+rlrc)**2)
       ENDDO
    ENDDO
    cellm=SQRT(cellm)

    xdim=clas8%cellcl(1)
    ydim=clas8%cellcl(2)*clas8%cellcl(1)
    zdim=clas8%cellcl(3)*clas8%cellcl(1)
    imx=MAX(INT(xdim/cellm),1)
    imy=MAX(INT(ydim/cellm),1)
    imz=MAX(INT(zdim/cellm),1)
    dx=xdim/REAL(imx,kind=real_8)
    dy=xdim/REAL(imy,kind=real_8)
    dz=xdim/REAL(imz,kind=real_8)
    cdimr12(1)=imx
    cdimr12(2)=imy
    cdimr12(3)=imz

    dx=1._real_8/dx
    dy=1._real_8/dy
    dz=1._real_8/dz
    DO ia=1,clas3%nclatom
       CALL pbc_clas(clasc(1,ia),clasc(2,ia),clasc(3,ia),xp,yp,zp)
       cellr12(1,ia)=INT((xp+xdim*0.5_real_8)*dx)+1
       cellr12(2,ia)=INT((yp+ydim*0.5_real_8)*dy)+1
       cellr12(3,ia)=INT((zp+zdim*0.5_real_8)*dz)+1
       IF (cellr12(1,ia).LT.1.OR.cellr12(1,ia).GT.cdimr12(1))&
            CALL stopgm('CELL_LISTR12','CELL NUMBER 1',& 
            __LINE__,__FILE__)
       IF (cellr12(2,ia).LT.1.OR.cellr12(2,ia).GT.cdimr12(2))&
            CALL stopgm('CELL_LISTR12','CELL NUMBER 2',& 
            __LINE__,__FILE__)
       IF (cellr12(3,ia).LT.1.OR.cellr12(3,ia).GT.cdimr12(3))&
            CALL stopgm('CELL_LISTR12','CELL NUMBER 3',& 
            __LINE__,__FILE__)
    ENDDO

    RETURN
  END SUBROUTINE cell_listr12

  FUNCTION check_cell(c1,c2,Cm)
    INTEGER                                  :: c1(3), c2(3), Cm(3)
    LOGICAL                                  :: check_cell

    INTEGER                                  :: dc, i, ima(3)

    DO i=1,3
       ima(i)=0
       dc=ABS(c1(i)-c2(i))
       IF (dc.LE.1) ima(i)=1
       IF (dc.GE.Cm(i)-1) ima(i)=1
    ENDDO
    check_cell=ima(1)*ima(2)*ima(3).EQ.1

    RETURN
  END FUNCTION check_cell

END MODULE clas_force_utils
