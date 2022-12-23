C     ==================================================================
      SUBROUTINE mm_center(TAU,c_trans,be_verbose)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
C     == Centering the QM system                                      ==
C     ==--------------------------------------------------------------==
      use system
      use parac
      use cell
      use rmas
      use mm_dimmod
      use isos
      use mm_input
      use bsym
      IMPLICIT NONE

      real*8  tau(3,maxsys%nax,*)
      real*8  c_trans(3)
      logical be_verbose
c
      real*8 xmin,ymin,zmin
      real*8 xmax,ymax,zmax
      real*8 cx,cy,cz
      real*8 PMATq
      integer IS,IA
      integer count
      save count
      data count /0/

      if(.not.isos1%tcent)return
      xmax=-1.D9
      ymax=-1.D9
      zmax=-1.D9
      xmin= 1.D9
      ymin= 1.D9
      zmin= 1.D9
      DO IS=1,mmdim%NSPq
        DO IA=1,NAq(IS)
          xmax=max(xmax,TAU(1,IA,IS))
          ymax=max(ymax,TAU(2,IA,IS))
          zmax=max(zmax,TAU(3,IA,IS))
          xmin=min(xmin,TAU(1,IA,IS))
          ymin=min(ymin,TAU(2,IA,IS))
          zmin=min(zmin,TAU(3,IA,IS))
        ENDDO
      ENDDO
      IF(be_verbose) count=count+1
cMcB
c     VOLCEL= (xmax-xmin)*(ymax-ymin)*(zmax-zmin) 
      CELL_COM%VOLCEL=1.0D0
      if ( (xmax-xmin).gt.0.0D0 ) CELL_COM%VOLCEL=
     &     CELL_COM%VOLCEL*(xmax-xmin)
      if ( (ymax-ymin).gt.0.0D0 ) CELL_COM%VOLCEL=
     &     CELL_COM%VOLCEL*(ymax-ymin)
      if ( (zmax-zmin).gt.0.0D0 ) CELL_COM%VOLCEL=
     &     CELL_COM%VOLCEL*(zmax-zmin)
cMcB
      if(paral%io_parent .and. (count.eq.1))then
        if(lqmmm%qmmm_reflex) then
         write(6,'(/1X,''USE REFLECTING WALLS TO CONFINE QM SYSTEM'')')
         write(6,'(/1X,''RE-CENTERING AT EACH TIME STEP TURNED OFF''/)')
        else
         write(6,'(/1X,''RE-CENTERING QM SYSTEM AT EVERY TIME STEP''/)')
        endif
      endif
      if(paral%io_parent .and. (MOD(count,100).eq.1))then
        write(6,*)'  BOX TOLERANCE [a.u.]  ', rqmmm%box_toll
        write(6,'(/1X,''BOX SIZE [a.u.]'',21X,
     .                                      ''QM SYSTEM SIZE [a.u.]'')')
        write(6,61)CELL_COM%CELLDM(1),xmax-xmin
        write(6,62)CELL_COM%CELLDM(2)*CELL_COM%CELLDM(1),ymax-ymin
        write(6,63)CELL_COM%CELLDM(3)*CELL_COM%CELLDM(1),zmax-zmin
      endif
61    format('  X DIRECTION:  CELLDIM = ',f10.4,'; XMAX-XMIN= ',f10.4)
62    format('  Y DIRECTION:  CELLDIM = ',f10.4,'; YMAX-YMIN= ',f10.4)
63    format('  Z DIRECTION:  CELLDIM = ',f10.4,'; ZMAX-ZMIN= ',f10.4)

      IF(paral%io_parent.and.be_verbose)THEN
        if(xmax-xmin .ge. cell_com%celldm(1)-rqmmm%box_toll) then
            write(6,*)'WARNING!!!  Box too small in x direction'
            write(6,61)CELL_COM%CELLDM(1),xmax-xmin
        endif
        if(ymax-ymin .ge. cell_com%celldm(2)*CELL_COM%CELLDM(1)
     &       -rqmmm%box_toll) then
            write(6,*)'WARNING!!!  Box too small in y direction'
            write(6,62)CELL_COM%CELLDM(2)*CELL_COM%CELLDM(1),ymax-ymin
        endif
        if(zmax-zmin .ge. cell_com%celldm(3)*CELL_COM%CELLDM(1)
     &       -rqmmm%box_toll) then
            write(6,*)'WARNING!!!  Box too small in z direction'
            write(6,63)CELL_COM%CELLDM(3)*CELL_COM%CELLDM(1),zmax-zmin
        endif
      ENDIF
 
      CX=(CELL_COM%CELLDM(1)-xmax-xmin)*0.5D0
      CY=(CELL_COM%CELLDM(2) * CELL_COM%CELLDM(1)-ymax-ymin)*0.5D0
      CZ=(CELL_COM%CELLDM(3) * CELL_COM%CELLDM(1)-zmax-zmin)*0.5D0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     PMATq =0.D0
!     do is=1,NSPq
!       PMATq = PMATq + NAq(is)
!     enddo
!     CX=0.D0
!     CY=0.D0
!     CZ=0.D0
!     DO IS=1,NSPq
!       DO IA=1,NAq(IS)
!         CX=CX+TAU(1,IA,IS)/PMATq
!         CY=CY+TAU(2,IA,IS)/PMATq
!         CZ=CZ+TAU(3,IA,IS)/PMATq
!       ENDDO
!     ENDDO
!     CX=0.5D0*CELLDM(1)-CX
!     CY=0.5D0*CELLDM(2)*CELLDM(1)-CY
!     CZ=0.5D0*CELLDM(3)*CELLDM(1)-CZ
!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      c_trans(1)=cx
      c_trans(2)=cy
      c_trans(3)=cz
      IF(CNTL%BSYMM.AND.(BSCLCS.EQ.2))THEN 
       IF(paral%io_parent)WRITE(6,'(/,A,/)')
     &   'MM_CENTER| WARNING! COORDINATES NOT SHIFTED FOR HS STATE'
       RETURN
      END IF
!$OMP parallel do private(IS,IA) shared(cx,cy,cz)
      DO IS=1,mmdim%NSPm
        DO IA=1,NAm(IS)
          TAU(1,IA,IS)=TAU(1,IA,IS)+cx
          TAU(2,IA,IS)=TAU(2,IA,IS)+cy
          TAU(3,IA,IS)=TAU(3,IA,IS)+cz
        ENDDO
      ENDDO
      RETURN
      END
