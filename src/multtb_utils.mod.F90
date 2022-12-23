MODULE multtb_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE symm,                            ONLY: irt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: multtb

CONTAINS

  ! ==================================================================
  ! IRT generates not-allocated usage error if INDPG=0
  SUBROUTINE multtb(indpg,nrot,ntvec,nat,xtable,&
       inve,multab,prt)
    ! ==--------------------------------------------------------------==
    ! ==  CONSTRUCT GROUP MULTIPLICATION TABLE                        ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   INDPG Group number (if 0 return otherwise not use)         ==
    ! ==   NROT  Number of rotations                                  ==
    ! ==   NTVEC Number of translation vectors associated with        ==
    ! ==         Identity rotation                                    ==
    ! ==         NTVEC /= 1 if the cell is not the primitive one      ==
    ! ==   NAT   Number of atoms                                      ==
    ! ==   XTABLE(3,3,NROT) Rotations                                 ==
    ! ==   INVE(120) Number of inverse operation for each rotation    ==
    ! ==   IRT(120,NAT)   Atom transformation table for rotations     ==
    ! ==   PRT .TRUE. PRINT SOME INFORMATION                          ==
    ! == OUTPUT:                                                      ==
    ! ==   MULTAB(120,120) multiplication table                       ==
    ! ==   INVE(NROT) inverse operation for each rotation             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, nrot, ntvec, nat
    REAL(real_8)                             :: xtable(3,3,nrot)
    INTEGER                                  :: inve(120), multab(120,120)
    LOGICAL                                  :: prt

    CHARACTER(len=15)                        :: fformat
    INTEGER :: i, i1, i2, i3, ia, iclass, id3max, idgen(12), idiff, ir, ira, &
      irb, irc, irtm(200), ix, iy, j, tclass(24,120), tcn(24)
    LOGICAL                                  :: lmtest
    REAL(real_8)                             :: xab(3,3)

! ,IRT(120,NAT)
! Variables
! ==--------------------------------------------------------------==

    IF (indpg.EQ.0) RETURN
    lmtest=.TRUE.
    ! ==--------------------------------------------------------------==
    ! Group multiplication table
    DO ira=1,nrot
       DO irb=1,nrot
          DO j=1,3
             DO i=1,3
                xab(i,j)=xtable(i,1,irb)*xtable(1,j,ira)+&
                     xtable(i,2,irb)*xtable(2,j,ira)+&
                     xtable(i,3,irb)*xtable(3,j,ira)
             ENDDO
          ENDDO
          DO irc=1,nrot
             idiff=0
             DO j=1,3
                DO i=1,3
                   idiff=idiff+NINT(ABS(xtable(i,j,irc)-xab(i,j)))
                ENDDO
             ENDDO
             IF (idiff.EQ.0) THEN
                multab(ira,irb)=irc
                GOTO 100
             ENDIF
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,*)&
               'MULTTB| CANNOT FIND PROPER MULTIPLICATION OF ',&
               ira,' AND ',irb
          CALL stopgm('MULTTB',' ',& 
               __LINE__,__FILE__)
100       CONTINUE
          ! ==---------------------------------------------------------==
          ! Test if the cell is a primary cell (or no replication).
          irc=multab(ira,irb)
          DO ia=1,nat
             irtm(ia)=irt(irb,irt(ira,ia))
          ENDDO
          idiff=0
          DO i=1,nat
             idiff=idiff+ABS(irt(irc,i)-irtm(i))
          ENDDO
          IF (idiff.NE.0) THEN
             IF (lmtest) THEN
                IF ((prt).AND.paral%io_parent)&
                     WRITE(6,*)&
                     'SUPERCELL ATOMS FORM NO BASIS FOR ',&
                     'THIS POINT GROUP'
             ENDIF
             lmtest=.FALSE.
          ENDIF
       ENDDO
    ENDDO
    IF (ntvec.NE.1) THEN
       IF (lmtest) THEN
          IF ((prt).AND.paral%io_parent)&
               WRITE(6,*)&
               'MULTTB| SUPERCELL ATOMS FORM NO BASIS FOR ',&
               'THIS POINT GROUP'
       ENDIF
       lmtest=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Look for the inverse of a symmetry operation
    DO ir=1,nrot
       DO ix=1,nrot
          IF (multab(ir,ix).EQ.1) THEN
             inve(ir)=ix
             GOTO 200
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)&
            'MULTTB| CANNOT FIND INVERSE OPERATION OF ',ir
       CALL stopgm('MULTTB',' ',& 
            __LINE__,__FILE__)
200    CONTINUE
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Generate the class structure of the group
    iclass=0
    DO ir=1,nrot
       DO i=1,iclass
          DO j=1,tcn(i)
             IF (ir.EQ.tclass(i,j)) GOTO 300
          ENDDO
       ENDDO
       iclass=iclass+1
       tcn(iclass)=0
       DO ix=1,nrot
          iy=multab(ix,ir)
          iy=multab(iy,inve(ix))
          DO i=1,tcn(iclass)
             IF (tclass(iclass,i).EQ.iy) GOTO 400
          ENDDO
          tcn(iclass)=tcn(iclass)+1
          tclass(iclass,tcn(iclass))=iy
400       CONTINUE
       ENDDO
300    CONTINUE
    ENDDO
    IF ((prt).AND.paral%io_parent)&
         WRITE(6,'(A,T63,I3)')&
         ' NUMBER OF IRREDUCIBLE REPRESENTATIONS:',iclass
    ! ==--------------------------------------------------------------==
    ! Calculation of degre of irreducible representation
    ! Sum d_{IR}^2=NROT with d=1,2 or 3
    ! Maximum number of IR of degree 3
    id3max=nrot/3**2
    DO i3=id3max,0,-1
       ! ID2 and ID1 are automatically determined
       i=nrot-i3*3**2-(iclass-i3)
       IF (MOD(i,3).EQ.0) THEN
          ! We got a solution
          i2=i/3
          i1=iclass-i3-i2
          GOTO 500
       ENDIF
    ENDDO
    ! Error
    CALL stopgm('MULTTB','DEGRES OF IR NO DETERMINED',& 
         __LINE__,__FILE__)
500 CONTINUE
    DO ir=1,iclass
       IF (ir.LE.i1) THEN
          idgen(ir)=1
       ELSEIF (ir.LE.(i1+i2)) THEN
          idgen(ir)=2
       ELSE
          idgen(ir)=3
       ENDIF
    ENDDO
    IF (prt) THEN
       IF (paral%io_parent)&
            WRITE(fformat,'(A,I2,A,I2,A)')&
            '(A,T',66-iclass*3,',',iclass,'I3)'
       IF (paral%io_parent)&
            WRITE(6,fformat)&
            ' DIMENSION OF IR:',(idgen(i),i=1,iclass)
    ENDIF
    iy=0
    DO i=1,iclass
       iy=iy+idgen(i)*idgen(i)
    ENDDO
    IF (iy.NE.nrot) CALL stopgm('MULTTB',&
         'INCORRECT CALUCLATION OF DIMENSIONS OF IR',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE multtb
  ! ==================================================================

END MODULE multtb_utils
