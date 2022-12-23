MODULE interp3d_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interp3d

CONTAINS

  ! ==================================================================
  SUBROUTINE interp3d(grid0,kx0,ky0,kz0,grid1,kx1,ky1,kz1)
    ! ==--------------------------------------------------------------==
    ! == SIMPLE INTERPOLATION 3D. (not efficient - T. Deutsch)        ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: kx0, ky0, kz0
    REAL(real_8)                             :: grid0(kx0,ky0,kz0)
    INTEGER                                  :: kx1, ky1, kz1
    REAL(real_8)                             :: grid1(kx1,ky1,kz1)

    REAL(real_8), PARAMETER                  :: densmax = 1.e1_real_8, &
                                                toll = -1.e-15_real_8 

    INTEGER :: incx, incy, incz, ix1, iy1, iz1, maxx1, maxy1, maxz1, minx1, &
      miny1, minz1, nnr0, nnr1, nx0, nx1, ny0, ny1, nz0, nz1
    REAL(real_8)                             :: px, py, pz, rx01, ry01, rz01, &
                                                sumgrid0, sumgrid1, VALUE, &
                                                x1, y1, z1
    REAL(real_8), EXTERNAL                   :: dasum

    CALL zeroing(grid1)!,kx1*ky1*kz1)
    nx0=kx0-MOD(kx0+1,2)
    ny0=ky0-MOD(ky0+1,2)
    nz0=kz0-MOD(kz0+1,2)
    nx1=kx1-MOD(kx1+1,2)
    ny1=ky1-MOD(ky1+1,2)
    nz1=kz1-MOD(kz1+1,2)
    nnr0=nx0*ny0*nz0
    nnr1=nx1*ny1*nz1
    rx01=REAL(nx0,kind=real_8)/REAL(nx1,kind=real_8)
    IF ((nx1-nx0).EQ.0) THEN
       incx=0
    ELSE
       incx=1
    ENDIF
    ry01=REAL(ny0,kind=real_8)/REAL(ny1,kind=real_8)
    IF ((ny1-ny0).EQ.0) THEN
       incy=0
    ELSE
       incy=1
    ENDIF
    rz01=REAL(nz0,kind=real_8)/REAL(nz1,kind=real_8)
    IF ((nz1-nz0).EQ.0) THEN
       incz=0
    ELSE
       incz=1
    ENDIF
    sumgrid1=0._real_8
    DO iz1=1,nz1
       z1=rz01*iz1
       minz1=INT(z1)
       maxz1=minz1+incz
       pz=z1-minz1
       IF (minz1.LT.1)   minz1=nz0
       IF (maxz1.GT.nz0) maxz1=1
       DO iy1=1,ny1
          y1=ry01*iy1
          miny1=INT(y1)
          maxy1=miny1+incy
          py=y1-miny1
          IF (miny1.LT.1)   miny1=ny0
          IF (maxy1.GT.ny0) maxy1=1
          DO ix1=1,nx1
             x1=rx01*ix1
             minx1=INT(x1)
             maxx1=minx1+incx
             px=x1-minx1
             IF (minx1.LT.1)   minx1=nx0
             IF (maxx1.GT.nx0) maxx1=1
             VALUE=(1._real_8-px)*(1._real_8-py)*(1._real_8-pz)*grid0(minx1,miny1,minz1)&
                  +(1._real_8-px)*      py *(1._real_8-pz)*grid0(minx1,maxy1,minz1)&
                  +(1._real_8-px)*(1._real_8-py)*      pz *grid0(minx1,miny1,maxz1)&
                  +(1._real_8-px)*      py *      pz *grid0(minx1,maxy1,maxz1)&
                  +      px *(1._real_8-py)*(1._real_8-pz)*grid0(maxx1,miny1,minz1)&
                  +      px *      py *(1._real_8-pz)*grid0(maxx1,maxy1,minz1)&
                  +      px *(1._real_8-py)*      pz *grid0(maxx1,miny1,maxz1)&
                  +      px *      py *      pz *grid0(maxx1,maxy1,maxz1)
             IF (VALUE.LT.toll) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) ' INTERP3D! ',ix1,iy1,iz1, VALUE
                CALL stopgm('INTERP3D','NEGATIVE DENSITY',& 
                     __LINE__,__FILE__)
             ENDIF
             grid1(ix1,iy1,iz1)=VALUE
             sumgrid1=sumgrid1+VALUE
          ENDDO
       ENDDO
    ENDDO
    sumgrid0=dasum(nnr0,grid0,1)
    sumgrid0=sumgrid0/REAL(nnr0,kind=real_8)
    sumgrid1=sumgrid1/REAL(nnr1,kind=real_8)
    IF (sumgrid1.GT.densmax*sumgrid0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' INTERP3D! ','SUM OF INTERPOLATED MESH=',sumgrid1
       CALL stopgm('INTERP3D','WRONG INTERPOLATION',& 
            __LINE__,__FILE__)
    ENDIF
    ! Renormalisation
    CALL dscal(nnr1,sumgrid0/sumgrid1,grid1(1,1,1),1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interp3d
  ! ==================================================================

END MODULE interp3d_utils
