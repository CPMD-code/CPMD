MODULE cnstpr_utils
  USE coninp_utils,                    ONLY: raddeg
  USE cotr,                            ONLY: &
       cnpar, cnsval, cnsval_dest, cotc0, cotr007, fv, grate, gsrate, ntcnst, &
       ntrest, resfor, respar, respos, resval, resval_dest
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert,&
                                             nat_grm
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cnstpr

CONTAINS

  ! ==================================================================
  SUBROUTINE cnstpr
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(len=10), DIMENSION(-1:1) :: &
      ltyp = (/'     MINUS','     FIXED','      PLUS'/)
    CHARACTER(len=10), DIMENSION(13) :: styp = (/'   STRETCH','      BEND',&
      '   TORSION','  DISTANCE','      OUTP','COORDINAT.','DIFFERENCE',&
      '    COORSP','   COOR_RF','BN. SWITCH','  TOT_COOR','   DISAXIS',&
      '    RESPOS'/)
    INTEGER                                  :: i, ia1, ia2, ia3, ia4, ityp
    LOGICAL                                  :: status
    REAL(real_8)                             :: cval, fval

! ==--------------------------------------------------------------==
! ==  PRINT SOME INFO                                             ==
! ==--------------------------------------------------------------==

    IF (cotc0%mcnstr.EQ.0.AND.cotr007%mrestr.EQ.0) RETURN
    CALL mm_dim(mm_go_mm,status)
    IF (cotc0%lshove) THEN
       IF (cotc0%mcnstr.GT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,20X,A)') ' <<<<< CONSTRAINTS >>>>>'
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' FIXED STRUCTURE ELEMENTS'
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') '      TYPE     ATOM ATOM ATOM ATOM',&
               '       VALUE       DIRECTION'
       ENDIF
       DO i=1,cotc0%mcnstr
          cval=cnsval(i)
          fval=fv(i)
          ityp=ntcnst(1,i)
          ia1=NAT_grm(ntcnst(2,i))
          ia2=NAT_grm(ntcnst(3,i))
          ia3=NAT_grm(ntcnst(4,i))
          ia4=NAT_grm(ntcnst(5,i))
          IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) THEN
             CALL raddeg(cval,1)
             CALL raddeg(fval,1)
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(A,4X,4I5,F12.5,6X,A10)') styp(ityp),ia1,ia2,ia3,ia4,&
               fval,ltyp(ntcnst(6,i))
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
    ELSE
       IF (cotc0%mcnstr.GT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,20X,A)') ' <<<<< CONSTRAINTS >>>>>'
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' FIXED STRUCTURE ELEMENTS'
          IF (paral%io_parent)&
               WRITE(6,'(A,7x,A)') '      TYPE     ATOM ATOM ATOM ATOM',&
               '         VALUE      DIFFERENCE'
       ENDIF
       DO i=1,cotc0%mcnstr
          cval=cnsval(i)
          fval=fv(i)
          ityp=ntcnst(1,i)
          ia1=NAT_grm(ntcnst(2,i))
          ia2=NAT_grm(ntcnst(3,i))
          ia3=NAT_grm(ntcnst(4,i))
          ia4=NAT_grm(ntcnst(5,i))
          IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) THEN
             CALL raddeg(cval,1)
             CALL raddeg(fval,1)
          ENDIF
          IF (ityp.EQ.6) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,I5,5X,2F10.5,2X,F10.5,1PE16.6)') styp(ityp),&
                  ia1,cnpar(1,i),cnpar(2,i),fval,fval-cval
          ELSEIF (ityp.EQ.8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,2I5,3F10.5,1PE16.6)') styp(ityp),ia1,ia2,&
                  cnpar(1,i),cnpar(2,i),fval,fval-cval
          ELSEIF (ityp.EQ.9 .OR. ityp.EQ. 11) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,4I5,2X,2F10.5,1PE16.6)') styp(ityp),&
                  ia1,ia2,ia3,ia4,cnpar(1,i),fval,fval-cval
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,4I5,10X,0PF12.5,1PE16.6)') styp(ityp),&
                  ia1,ia2,ia3,ia4,fval,fval-cval
          ENDIF

          IF (ABS(grate(i)).GT.1.e-10_real_8) THEN
             cval=cnsval_dest(i)
             IF (cval.NE.-999._real_8) THEN
                IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) THEN
                   CALL raddeg(cval,1)
                ENDIF
                grate(i)=SIGN(grate(i),cval-fval)
             ENDIF 
             IF (paral%io_parent)&
                  WRITE(6,'(A,10X,T22,A,10X,0PF18.5,A)') styp(ityp),&
                  " Growth Rate: ",grate(i),"/at.unit"
             IF (cval.NE.-999._real_8) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,10X,T22,A,10X,0PF18.5)') styp(ityp),&
                     " Destination: ",cval
             ENDIF
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (cotr007%mrestr.GT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,20X,A)') ' <<<<< RESTRAINTS >>>>> '
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' RESTRAINED STRUCTURE ELEMENTS'
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') '      TYPE     ATOM ATOM ATOM ATOM',&
               '     VALUE   FORCE CONSTANT'
       ENDIF
       DO i=1,cotr007%mrestr
          cval=resval(i)
          ityp=ntrest(1,i)
          ia1=NAT_grm(ntrest(2,i))
          ia2=NAT_grm(ntrest(3,i))
          ia3=NAT_grm(ntrest(4,i))
          ia4=NAT_grm(ntrest(5,i))

          IF (ityp.EQ.13) THEN! cmb-kk
             IF (paral%io_parent)&
                  WRITE(6,'(A,4X,I6,3F6.2,F10.5,2X,G15.4)') styp(ityp),&
                  ia1,respos(1,i),respos(2,i),respos(3,i),&
                  cval,resfor(i)
          ELSEIF (ityp.EQ.6) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,I5,5X,2F10.5,2X,F10.5,2X,G15.4)') styp(ityp),&
                  ia1,respar(1,i),respar(2,i),cval,resfor(i)
          ELSEIF (ityp.EQ.8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,2I5,3F10.5,G15.4)') styp(ityp),ia1,ia2,&
                  cnpar(1,i),cnpar(2,i),cval,resfor(i)
          ELSEIF (ityp.EQ.9.OR.ityp.EQ.11) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,3X,4I5,2X,2F10.5,G15.4)') styp(ityp),&
                  ia1,ia2,ia3,ia4,&
                  respar(1,i),cval,resfor(i)
          ELSE
             IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) CALL raddeg(cval,1)
             IF (paral%io_parent)&
                  WRITE(6,'(A,4X,4I6,F10.5,2X,G15.4)') styp(ityp),ia1,&
                  ia2,ia3,ia4,cval,resfor(i)
          ENDIF

          IF (ABS(gsrate(i)).GT.1.e-10_real_8) THEN
             fval=cval
             cval=resval_dest(i)
             IF (cval.NE.-999._real_8) THEN
                IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) CALL raddeg(cval,1)
                gsrate(i)=SIGN(gsrate(i),cval-fval)
             ENDIF
             IF (paral%io_parent)&
                  WRITE(6,'(A,10X,T22,A,10X,0PF18.5,A)') styp(ityp),&
                  " Growth Rate: ",gsrate(i),"/at.unit"
             IF (cval.NE.-999._real_8) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,10X,T22,A,10X,0PF18.5)') styp(ityp),&
                     " Destination: ",cval
             ENDIF
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,status)
    RETURN
  END SUBROUTINE cnstpr
  ! ==================================================================

END MODULE cnstpr_utils
