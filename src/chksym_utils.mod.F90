MODULE chksym_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE k290_2_utils,                    ONLY: group1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE molsym_utils,                    ONLY: molsym
  USE multtb_utils,                    ONLY: multtb
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring
  USE rmas,                            ONLY: rmass
  USE symm,                            ONLY: &
       irt, irtvec, isymu, iun, numshel, stag, symmi, symmr, symmt, tvec
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: maxsys,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dgive
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: updatsym
  PUBLIC :: symtau
  PUBLIC :: gentau 
  PUBLIC :: chksym
  PUBLIC :: givesym

CONTAINS

  ! ==================================================================
  SUBROUTINE givesym(tau0)
    ! ==--------------------------------------------------------------==
    ! == Find automatically symmetry from atomic coordinates (TAU0)   ==
    ! == and lattice vectors                                          ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'givesym'

    INTEGER                                  :: i, i1, i2, ia, ib(48), ierr, &
                                                ifirst = 0, ihc, iout, is, &
                                                isub, isy, li, nat0
    INTEGER, ALLOCATABLE                     :: f0(:,:), isc(:), ty(:)
    REAL(real_8)                             :: a0(3,3), b0(3,3), r(3,3,48), &
                                                temp(3,3)
    REAL(real_8), ALLOCATABLE                :: rx(:,:), tvec0(:,:), &
                                                xkapa(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('   GIVESYM',isub)
    ! ==--------------------------------------------------------------==
    ! == Allocations of arrays                                        ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(xkapa(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rx(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tvec0(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ty(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(isc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(f0(49,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nat0=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          nat0=nat0+1
          ty(nat0) = is
          DO i=1,3
             xkapa(i,nat0) = tau0(i,ia,is)
          ENDDO
       ENDDO
    ENDDO
    IF (paral%parent) THEN
       ! Message from GROUP1 to standard output
       iout=6
    ELSE
       ! No message for other processors
       iout=0
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL group1(iout,parm%a1,parm%a2,parm%a3,ions1%nat,ty,xkapa,b0(1,1),b0(1,2),b0(1,3),&
         symmi%ihg,ihc,isy,li,symmi%nrot,symmi%indpg,ib,symmi%ntvec,symmr%ftau,f0,r,tvec0,symmr%origin,rx,&
         isc,symmr%deltasym)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(xkapa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ty,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(isc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Save Atomic table transformation
    IF (ifirst.EQ.0) THEN
       IF(ALLOCATED(irt)) THEN
          DEALLOCATE(irt,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(irt(120,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF(ALLOCATED(tvec)) THEN
          DEALLOCATE(tvec,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(tvec(3,symmi%ntvec),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    CALL zeroing(irt)!,120*ions1%nat)
    DO i=1,symmi%nrot
       DO ia=1,ions1%nat
          irt(i,ia)=f0(i,ia)
       ENDDO
    ENDDO
    IF (symmi%ntvec.NE.0) THEN
       ! Save 49th row of F0 containing inequivalent atoms for transla.
       DO ia=1,ions1%nat
          irt(49,ia)=f0(49,ia)
       ENDDO
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'TRANSLATION VECTORS:'
          DO i=1,symmi%ntvec
             IF (paral%io_parent)&
                  WRITE(6,'(" TVEC(",I3,"): [",F6.3,", ",F6.3,", ",F6.3,"]")')&
                  i,tvec0(1,i),tvec0(2,i),tvec0(3,i)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF
       DEALLOCATE(tvec,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(tvec(3,(3*symmi%ntvec+1)/3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(3*symmi%ntvec,tvec0,1,tvec,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(f0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tvec0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF ( (symmi%indpg.NE.0).AND.(.NOT.symmt%tmsym) ) THEN
       IF (isy.EQ.0.OR.isy.EQ.-2) THEN
          symmt%tsymmorphic=.FALSE.
       ELSEIF (isy.EQ.1) THEN
          symmt%tsymmorphic=.TRUE.
          symmt%torigin=.TRUE.
       ELSEIF (isy.EQ.-1) THEN
          symmt%tsymmorphic=.TRUE.
          symmt%torigin=.FALSE.
       ENDIF
       DO i = 1,3
          a0(i,1) = parm%a1(i)/parm%alat
          a0(i,2) = parm%a2(i)/parm%alat
          a0(i,3) = parm%a3(i)/parm%alat
       ENDDO
       DO i2=1,3
          DO i1=1,3
             b0(i1,i2) = parm%alat*b0(i1,i2)
          ENDDO
       ENDDO
       IF (li.GT.0) THEN
          symmi%inversion=li
       ELSE
          symmi%inversion=0
       ENDIF
       ! Transpose(B).A0 = 1
       ! Put in crystal coordinates
       DO i=1,symmi%nrot
          CALL dgemm('N','N',3,3,3,1._real_8,r(1,1,ib(i)),3,a0(1,1),3,0._real_8,&
               temp(1,1),3)
          CALL dgemm('T','N',3,3,3,1._real_8,b0(1,1),3,temp(1,1),3,0._real_8,&
               symmr%xtable(1,1,i),3)
       ENDDO
    ELSEIF (symmi%indpg.NE.0.AND.symmt%tmsym) THEN
       CALL molsym(stag,symmi%naxis,symmr%xtable,symmi%nrot)
    ENDIF
    CALL tihalt('   GIVESYM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE givesym
  ! ==================================================================
  SUBROUTINE chksym(tau0)
    ! ==--------------------------------------------------------------==
    ! == Check symmetry operations and                                ==
    ! == find fractional translation vectors, translation vectors     ==
    ! == associated with Identity symmetry                            ==
    ! == Built some arrays used for symmetrisation of some quantities ==
    ! == (forces, density)                                            ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'chksym'

    CHARACTER(len=12), DIMENSION(7) :: ihgp = (/'   TRICLINIC','  MONOCLINIC',&
      'ORTHORHOMBIC','  TETRAGONAL','       CUBIC','    TRIGONAL',&
      '   HEXAGONAL'/)
    INTEGER                                  :: i, ia, iat, ib, ic, ierr, &
                                                indpg0, invers, ir, is, iv, &
                                                j, k, lxau, nsym
    INTEGER, ALLOCATABLE                     :: irt0(:), isc(:), ityp(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: nodupli, oksym, sym(120), &
                                                tsymmor2
    REAL(real_8)                             :: aa(3,3), bb(3,3), cm(3), &
                                                deltasym2, ft(3), sum, &
                                                temp(3,3), vs
    REAL(real_8), ALLOCATABLE                :: rau(:,:), xau(:,:)

! ==--------------------------------------------------------------==

    IF (symmi%indpg.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    deltasym2=1.e-4_real_8*symmr%deltasym
    lxau=3*ions1%nat*2              ! Used as scratch array for symvec
    ALLOCATE(xau(3,lxau/3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ityp(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(isc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(irt0(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (ifirst.EQ.0) THEN
       IF (.NOT.symmt%tpgauto.AND..NOT.symmt%tsymm) THEN
          ALLOCATE(irt(120,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF(ALLOCATED(isymu)) THEN
          DEALLOCATE(isymu,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(isymu(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (symmt%tmsym) THEN
       ! MOVE CENTER OF MASS TO ORIGIN
       CALL zeroing(cm)!,3)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             cm(1)=cm(1)+tau0(1,ia,is)*rmass%pma(is)
             cm(2)=cm(2)+tau0(2,ia,is)*rmass%pma(is)
             cm(3)=cm(3)+tau0(3,ia,is)*rmass%pma(is)
          ENDDO
       ENDDO
       cm(1)=cm(1)/rmass%pmatot
       cm(2)=cm(2)/rmass%pmatot
       cm(3)=cm(3)/rmass%pmatot
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             tau0(1,ia,is)=tau0(1,ia,is)-cm(1)
             tau0(2,ia,is)=tau0(2,ia,is)-cm(2)
             tau0(3,ia,is)=tau0(3,ia,is)-cm(3)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ATOMIC COORDINATES IN DIRECT LATTICE VECTOR BASIS
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          ityp(iat)=is
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,tau0(1,ia,is),1,0.0_real_8,&
               xau(1,iat),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Check translation vectors
    IF (symmt%tpgauto) THEN
       ! We know the number of translation vectors:
       ! Calculation of the atom transf. table
    ELSE
       symmi%ntvec=1
       rau(1,1)=0._real_8
       rau(2,1)=0._real_8
       rau(3,1)=0._real_8
       ! We do not know the number of translation vectors:
       DO ia=2,ions1%nat
          IF (ityp(1).EQ.ityp(ia)) THEN
             ft(1)=xau(1,ia)-xau(1,1)-NINT(xau(1,ia)-xau(1,1))
             ft(2)=xau(2,ia)-xau(2,1)-NINT(xau(2,ia)-xau(2,1))
             ft(3)=xau(3,ia)-xau(3,1)-NINT(xau(3,ia)-xau(3,1))
             CALL checksym(ions1%nat,ityp,xau,xau,ft,.TRUE.,&
                  oksym,irt0,isc,symmr%deltasym)
             IF (oksym) THEN
                symmi%ntvec=symmi%ntvec+1
                rau(1,symmi%ntvec)=ft(1)
                rau(2,symmi%ntvec)=ft(2)
                rau(3,symmi%ntvec)=ft(3)
             ENDIF
          ENDIF
       ENDDO
       IF (symmi%ntvec.NE.0) THEN
          IF (ifirst.EQ.0) THEN
             ALLOCATE(tvec(3,symmi%ntvec),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          DO ir=1,symmi%ntvec
             tvec(1,ir)=rau(1,ir)
             tvec(2,ir)=rau(2,ir)
             tvec(3,ir)=rau(3,ir)
          ENDDO
       ENDIF
    ENDIF
    ! Calculation of IRTVEC(NTVEC,NAT):
    ! Atomic tranformation table for translation vectors
    IF (ifirst.EQ.0) THEN
       IF(ALLOCATED(irtvec)) THEN
          DEALLOCATE(irtvec,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(irtvec(ions1%nat,symmi%ntvec),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(irtvec,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(irtvec(ions1%nat,symmi%ntvec),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DO ir=1,symmi%ntvec
       CALL checksym(ions1%nat,ityp,xau,xau,tvec(1,ir),&
            .TRUE.,oksym,irtvec(1,ir),isc,symmr%deltasym)
       IF (.NOT.oksym.AND.paral%parent) THEN
          CALL stopgm('CHKSYM','WRONG TRANSLATION VECTORS',& 
               __LINE__,__FILE__)
       ENDIF
    ENDDO
    nodupli=symmi%ntvec.EQ.1
    ! ==--------------------------------------------------------------==
    ! LOOP OVER SYMMETRY OPERATIONS
    IF (symmt%tpgauto) tsymmor2=symmt%tsymmorphic
    symmt%tsymmorphic=.TRUE.
    DO ir=1,symmi%nrot
       symmt%tsymmor(ir)=.TRUE.
    ENDDO
    vs=0._real_8
    DO ir=1,symmi%nrot
       DO ia=1,ions1%nat
          DO k=1,3
             rau(k,ia)=symmr%xtable(k,1,ir)*xau(1,ia)+&
                  symmr%xtable(k,2,ir)*xau(2,ia)+&
                  symmr%xtable(k,3,ir)*xau(3,ia)
          ENDDO
       ENDDO
       IF (symmt%tpgauto) THEN
          ! TRY TO FIND A FRACTIONAL TRANSLATION VECTOR
          ib=irt(ir,1)
          IF (ityp(1).EQ.ityp(ib)) THEN
             ft(1)=xau(1,ib)-rau(1,1)-NINT(xau(1,ib)-rau(1,1))
             ft(2)=xau(2,ib)-rau(2,1)-NINT(xau(2,ib)-rau(2,1))
             ft(3)=xau(3,ib)-rau(3,1)-NINT(xau(3,ib)-rau(3,1))
             CALL checksym(ions1%nat,ityp,xau,rau,ft,nodupli,&
                  sym(ir),irt0,isc,symmr%deltasym)
             IF (sym(ir)) THEN
                IF (ABS(ft(1))+ABS(ft(2))+ABS(ft(3)).GE.symmr%deltasym) THEN
                   vs=vs+ABS(ft(1))+ABS(ft(2))+ABS(ft(3))
                   symmt%tsymmorphic=.FALSE.
                   symmt%tsymmor(ir)=.FALSE.
                ENDIF
                GOTO 100
             ELSE
                IF (paral%parent) CALL stopgm('CHKSYM',&
                     'WRONG TRANSLATION VECTOR FROM AUTO',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSE
             IF (paral%parent) CALL stopgm('CHKSYM',&
                  'WRONG ATOMIC TRANSFORMATION TABLE FROM AUTO',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          CALL zeroing(ft)!,3)
          ! CHECK FOR THIS SYMMETRY WITH FT=0
          ! (IMPORTANT IF NOT UNIQUE ATOM)
          CALL checksym(ions1%nat,ityp,xau,rau,ft,nodupli,&
               sym(ir),irt0,isc,symmr%deltasym)
          IF (.NOT.sym(ir).AND..NOT.symmt%tmsym) THEN
             DO ib=1,ions1%nat
                IF (ityp(1).EQ.ityp(ib)) THEN
                   ft(1)=xau(1,1)-rau(1,ib)-NINT(xau(1,1)-rau(1,ib))
                   ft(2)=xau(2,1)-rau(2,ib)-NINT(xau(2,1)-rau(2,ib))
                   ft(3)=xau(3,1)-rau(3,ib)-NINT(xau(3,1)-rau(3,ib))
                   CALL checksym(ions1%nat,ityp,xau,rau,ft,nodupli,sym(ir),&
                        irt0,isc,symmr%deltasym)
                   IF (sym(ir)) THEN
                      vs=vs+ABS(ft(1))+ABS(ft(2))+ABS(ft(3))
                      symmt%tsymmorphic=.FALSE.
                      symmt%tsymmor(ir)=.FALSE.
                      GOTO 100
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDIF
100    CONTINUE
       IF (sym(ir)) THEN
          IF (symmt%tpgauto) THEN
             ! FT and FTAU should be the same for the same atomic table
             ! multiplication.
             DO ia=1,ions1%nat
                IF (irt(ir,ia).NE.irt0(ia)) THEN
                   IF (paral%parent) CALL stopgm('CHKSYM',&
                        'WRONG ATOMIC TRANSFORMATION TABLE',& 
                        __LINE__,__FILE__)
                ENDIF
             ENDDO
             IF (  (ABS((symmr%ftau(1,ir)-ft(1)) - NINT(symmr%ftau(1,ir)-ft(1)))&
                  .GT.symmr%deltasym).OR.&
                  (ABS((symmr%ftau(2,ir)-ft(2)) - NINT(symmr%ftau(2,ir)-ft(2)))&
                  .GT.symmr%deltasym).OR.&
                  (ABS((symmr%ftau(3,ir)-ft(3)) - NINT(symmr%ftau(3,ir)-ft(3)))&
                  .GT.symmr%deltasym) ) THEN
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(" IR=",I3,T10,"FTAU=",3F12.6)')&
                        ir,(symmr%ftau(k,ir),k=1,3)
                   IF (paral%io_parent)&
                        WRITE(6,'(T10,"  FT=",3F12.6)') (ft(k),k=1,3)
                   ! FT and FTAU should be the same
                   IF (paral%parent)&
                        CALL stopgm('CHKSYM','WRONG TRANSLATION VECTOR',& 
                        __LINE__,__FILE__)
                ENDIF
             ENDIF
          ELSE
             DO ia=1,ions1%nat
                irt(ir,ia)=irt0(ia)
             ENDDO
             symmr%ftau(1,ir) = ft(1)
             symmr%ftau(2,ir) = ft(2)
             symmr%ftau(3,ir) = ft(3)
          ENDIF
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    nsym=0
    DO ir=1,symmi%nrot
       IF (sym(ir))nsym=nsym+1
       IF (.NOT.sym(ir)) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,"! !!",I3," OPERATION NOT CORRECT!!!")') iR
             IF (symmt%tpgauto) CALL stopgm('CHKSYM',&
                  'WRONG AUTOMATIC DETERMINATION OF SYMMETRY',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
    ENDDO
    IF (nsym.NE.symmi%nrot) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,64("! "))')
          IF (paral%io_parent)&
               WRITE(6,'(" ! !",A,T64,"!!")')&
               ' CHKSYM| NUMBER OF SYMMETRY OPERATIONS HAS BEEN CHANGED'
          IF (paral%io_parent)&
               WRITE(6,'(" ! !",T13,A,I2,A,I2,T64,"!!")')&
               'FROM ',symmi%nrot,' TO ',nsym
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("! "))')
       ENDIF
       is=0
       DO ir=1,symmi%nrot
          IF (sym(ir)) THEN
             is=is+1
             symmt%tsymmor(is)=symmt%tsymmor(ir)
             DO j=1,3
                DO i=1,3
                   symmr%xtable(i,j,is)=symmr%xtable(i,j,ir)
                ENDDO
             ENDDO
             DO i=1,3
                symmr%ftau(i,is)=symmr%ftau(i,ir)
             ENDDO
             DO i=1,ions1%nat
                irt(is,i)=irt(ir,i)
             ENDDO
          ENDIF
       ENDDO
       symmi%nrot=nsym
       symmi%indpg=-symmi%indpg
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Check if there is inversion symmetry
    invers=0
    DO ir=1,symmi%nrot
       IF (  (ABS(symmr%xtable(1,1,ir)+1._real_8).LT.deltasym2).AND.&
            (ABS(symmr%xtable(2,2,ir)+1._real_8).LT.deltasym2).AND.&
            (ABS(symmr%xtable(3,3,ir)+1._real_8).LT.deltasym2) ) THEN
          sum=   ABS(symmr%xtable(1,2,ir))&
               + ABS(symmr%xtable(1,3,ir))&
               + ABS(symmr%xtable(2,1,ir))&
               + ABS(symmr%xtable(2,3,ir))&
               + ABS(symmr%xtable(3,1,ir))&
               + ABS(symmr%xtable(3,2,ir))
          IF (sum.LT.deltasym2) THEN
             invers=ir
          ENDIF
       ENDIF
    ENDDO
    IF (symmt%tpgauto) THEN
       IF (invers.NE.symmi%inversion) THEN
          IF (paral%parent) CALL stopgm('CHKSYM',&
               'INCORRECT DETERMINATION OF INVERSION',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    symmi%inversion=invers
    ! ==--------------------------------------------------------------==
    ! Look for center of symmetry (ORIGIN)
    IF (symmt%tpgauto) THEN
       IF (symmt%tsymmorphic.AND.symmt%torigin) THEN
          IF (  (symmt%tsymmorphic.AND..NOT.tsymmor2).OR.&
               (.NOT.symmt%tsymmorphic.AND.tsymmor2) ) THEN
             IF (paral%parent) CALL stopgm('CHKSYM',&
                  'INCORRECT DETERMINATION OF SPACE GROUP TYPE',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
    ELSE
       IF (symmt%tsymmorphic) THEN
          symmt%torigin=.TRUE.
       ELSE
          symmt%torigin=.FALSE.
       ENDIF
       DO i=1,3
          symmr%origin(i)=0._real_8
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! CONSTRUCTION OF FTABLE (ROTATION MATRICES FOR FOURIER SPACE)
    ! `B.A0 = 1 (`B=Transpose of B) [HT=A=`B¯¹, HTM1= B=`A¯¹]
    ! Put in reciprocal coordinates : 
    ! FTABLE = `B¯¹.`A.XTABLE.A¯¹.B = A.`A.XTABLE.`B.B
    CALL dgemm('N','T',3,3,3,1._real_8/(parm%alat*parm%alat),&
         metr_com%ht(1,1),3,metr_com%ht(1,1),3,0._real_8,aa(1,1),3)
    CALL dgemm('T','N',3,3,3,(parm%alat*parm%alat),&
         metr_com%htm1(1,1),3,metr_com%htm1(1,1),3,0._real_8,bb(1,1),3)
    DO ir=1,symmi%nrot
       CALL dgemm('N','N',3,3,3,1._real_8,&
            symmr%xtable(1,1,ir),3,bb(1,1),3,0._real_8,temp(1,1),3)
       CALL dgemm('N','N',3,3,3,1._real_8,&
            aa(1,1),3,temp(1,1),3,0._real_8,symmr%ftable(1,1,ir),3)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! LOOK FOR SYMMETRY UNIQUE ATOMS
    DO ia=1,ions1%nat
       isymu(ia)=1
    ENDDO
    ! Unique atoms for rotations
    DO ia=1,ions1%nat
       IF (isymu(ia).NE.0) THEN
          DO ir=1,symmi%nrot
             ib=irt(ir,ia)
             IF (ib.GT.ia) THEN
                IF (isymu(ib).NE.0) THEN
                   isymu(ia)=isymu(ia)+1
                   isymu(ib)=0
                ENDIF
             ENDIF
             ! Unique atoms for translation vectors
             IF (symmi%ntvec.NE.1) THEN
                DO iv=1,symmi%ntvec
                   ic=irtvec(ib,iv)
                   IF (ic.GT.ia) THEN
                      IF (isymu(ic).NE.0) THEN
                         isymu(ia)=isymu(ia)+1
                         isymu(ic)=0
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    symmi%iunique=0
    DO ia=1,ions1%nat
       IF (isymu(ia).NE.0) symmi%iunique=symmi%iunique+1
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (symmt%tmsym) THEN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             tau0(1,ia,is)=tau0(1,ia,is)+cm(1)
             tau0(2,ia,is)=tau0(2,ia,is)+cm(2)
             tau0(3,ia,is)=tau0(3,ia,is)+cm(3)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (.NOT.symmt%tpgauto) THEN
       indpg0=ABS(symmi%indpg)
       IF (indpg0.GE.1.AND.indpg0.LE.2) THEN
          symmi%ihg=1
       ELSEIF (indpg0.GE.3.AND.indpg0.LE.5) THEN
          symmi%ihg=2
       ELSEIF (indpg0.GE.6.AND.indpg0.LE.10) THEN
          symmi%ihg=6
       ELSEIF (indpg0.GE.11.AND.indpg0.LE.17) THEN
          symmi%ihg=4
       ELSEIF (indpg0.GE.18.AND.indpg0.LE.24) THEN
          symmi%ihg=7
       ELSEIF (indpg0.GE.25.AND.indpg0.LE.27) THEN
          symmi%ihg=3
       ELSEIF (indpg0.GE.28.AND.indpg0.LE.32) THEN
          symmi%ihg=5
       ELSEIF (indpg0.EQ.100) THEN
          symmi%ihg=100
       ELSE
          CALL stopgm('CHKSYM','BAD GROUP NUMBER',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%parent) THEN
          IF (parm%ibrav.EQ.4 .AND. symmi%nrot.EQ.24) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,A)')&
                  ' THE POINT GROUP OF THE CRYSTAL IS THE FULL ',&
                  'HEXAGONAL GROUP'
          ELSEIF (parm%ibrav.EQ.1 .AND. symmi%nrot.EQ.48) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,A)')&
                  ' THE POINT GROUP OF THE CRYSTAL IS THE FULL ',&
                  'CUBIC GROUP'
          ELSEIF (symmi%ihg.EQ.100) THEN
             CALL xstring(ihgp(symmi%ihg),ia,ib)
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,I2,A)')&
                  ' MOLECULAR SYSTEM  WITH ',symmi%nrot,' OPERATIONS'
          ELSE
             CALL xstring(ihgp(symmi%ihg),ia,ib)
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,A,A,I2,A)')&
                  ' THE CRYSTAL SYSTEM IS ',&
                  ihgp(symmi%ihg)(ia:ib),' WITH ',symmi%nrot,' OPERATIONS'
          ENDIF
          IF (symmi%ihg.NE.100) THEN
             IF (symmt%tsymmorphic) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A)')&
                     ' THE SPACE GROUP OF THE CRYSTAL IS SYMMORPHIC'
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A,/,A,A,/,3X,A,F15.6,A)')&
                     ' THE SPACE GROUP IS NON-SYMMORPHIC,',&
                     ' OR ELSE A NON STANDARD ORIGIN OF COORDINATES WAS',&
                     ' USED.',&
                     ' (SUM OF TRANSLATION VECTORS=',vs,')'
             ENDIF
             IF (symmi%ntvec.EQ.1) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,T60,I6)')&
                     ' NUMBER OF PRIMITIVE CELL:',symmi%ntvec
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A,T60,I6)')&
                     ' NUMBER OF PRIMITIVE CELLS:',symmi%ntvec
             ENDIF
             IF (symmi%inversion.NE.0.AND.paral%io_parent)&
                  WRITE(6,'(1X,A)')&
                  'THE POINT GROUP OF THE CRYSTAL CONTAINS THE INVERSION'
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF
    ENDIF
    IF (paral%parent) THEN
       IF (symmi%iunique.EQ.1) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,T61,I5)')&
               ' SYMMETRY UNIQUE (INEQUIVALENT) ATOM:',symmi%iunique
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,T61,I5)')&
               ' SYMMETRY UNIQUE (INEQUIVALENT) ATOMS:',symmi%iunique
          IF (symmi%iunique.NE.ions1%nat) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,6X)',advance="no") ' INDEXES:'
             j=0
             DO i=1,ions1%nat
                IF (isymu(i).NE.0) THEN
                   j=j+1
                   IF (j.LE.10) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(I4)',advance="no") i
                   ELSE
                      IF (paral%io_parent)&
                           WRITE(6,'(/,15X,I4)',advance="no") i
                      j=0
                   ENDIF
                ENDIF
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,*)
          ENDIF
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,1PE10.2)')&
            ' REQUIRED PRECISION FOR SYMMETRY:',symmr%deltasym
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       ifirst=1
    ENDIF
    DEALLOCATE(xau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ityp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(isc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(irt0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE chksym
  ! ==================================================================
  SUBROUTINE symtau(tau0)
    ! ==--------------------------------------------------------------==
    ! == Symmetrisation of atomic coordinates from symmetry operations==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'symtau'

    INTEGER                                  :: ia, iat, ib, idamax, ierr, &
                                                imax, ir, is, k, lxau, nsym
    INTEGER, ALLOCATABLE                     :: irt0(:), isc(:), ityp(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: sym(120)
    REAL(real_8)                             :: cm(3), disp, ft(3), vs
    REAL(real_8), ALLOCATABLE                :: rau(:,:), xau(:,:)

! ==--------------------------------------------------------------==

    IF (.NOT.symmt%tsymm) RETURN
    IF (symmi%indpg.EQ.0) RETURN
    IF (paral%io_parent)&
         WRITE(6,'(A)')&
         ' SYMTAU| SYMMETRIZATION OF ATOMIC COORDINATES '
    ! ==--------------------------------------------------------------==
    lxau=6*maxsys%nax*maxsys%nsx
    ALLOCATE(xau(3,2*maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rau(3,2*maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ityp(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(isc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(irt0(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       IF (.NOT.symmt%tpgauto) THEN
          ALLOCATE(irt(120,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ifirst=1
    ! ==--------------------------------------------------------------==
    ! MOVE CENTER OF MASS TO ORIGIN
    CALL zeroing(cm)!,3)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          cm(1)=cm(1)+tau0(1,ia,is)*rmass%pma(is)
          cm(2)=cm(2)+tau0(2,ia,is)*rmass%pma(is)
          cm(3)=cm(3)+tau0(3,ia,is)*rmass%pma(is)
       ENDDO
    ENDDO
    cm(1)=cm(1)/rmass%pmatot
    cm(2)=cm(2)/rmass%pmatot
    cm(3)=cm(3)/rmass%pmatot
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          tau0(1,ia,is)=tau0(1,ia,is)-cm(1)
          tau0(2,ia,is)=tau0(2,ia,is)-cm(2)
          tau0(3,ia,is)=tau0(3,ia,is)-cm(3)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ATOMIC COORDINATES IN DIRECT LATTICE VECTOR BASIS
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          ityp(iat)=is
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,tau0(1,ia,is),1,0.0_real_8,&
               xau(1,iat),1)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! LOOP TO FIND MATCHING TOLERANCE
200 CONTINUE
    ! ==--------------------------------------------------------------==
    ! LOOP OVER SYMMETRY OPERATIONS
    CALL zeroing(symmr%ftau)!,3*symmi%nrot)
    vs=0._real_8
    symmt%tsymmorphic=.TRUE.
    DO ir=1,symmi%nrot
       symmt%tsymmor(ir)=.TRUE.
    ENDDO
    DO ir=1,symmi%nrot
       DO ia=1,ions1%nat
          DO k=1,3
             rau(k,ia)=symmr%xtable(k,1,ir)*xau(1,ia)+&
                  symmr%xtable(k,2,ir)*xau(2,ia)+&
                  symmr%xtable(k,3,ir)*xau(3,ia)
          ENDDO
       ENDDO
       CALL zeroing(ft)!,3)
       ! CHECK FOR THIS SYMMETRY WITH FT=0 (IMPORTANT IF NOT UNIQUE ATOM)
       CALL checksym(ions1%nat,ityp,xau,rau,ft,.TRUE.,sym(ir),irt0,&
            isc,symmr%deltasym)
       IF (.NOT.sym(ir).AND..NOT.symmt%tmsym) THEN
          ! TRY TO FIND A FRACTIONAL TRANSLATION VECTOR
          DO ib=1,ions1%nat
             IF (ityp(1).EQ.ityp(ib)) THEN
                ft(1)=xau(1,1)-rau(1,ib)-NINT(xau(1,1)-rau(1,ib))
                ft(2)=xau(2,1)-rau(2,ib)-NINT(xau(2,1)-rau(2,ib))
                ft(3)=xau(3,1)-rau(3,ib)-NINT(xau(3,1)-rau(3,ib))
                CALL checksym(ions1%nat,ityp,xau,rau,ft,.TRUE.,&
                     sym(ir),irt0,isc,symmr%deltasym)
                IF (sym(ir)) THEN
                   vs=vs+ABS(ft(1))+ABS(ft(2))+ABS(ft(3))
                   symmt%tsymmorphic=.FALSE.
                   symmt%tsymmor(ir)=.FALSE.
                   GOTO 100
                ENDIF
             ENDIF
          ENDDO
       ENDIF
100    CONTINUE
       IF (sym(ir)) THEN
          DO iat=1,ions1%nat
             irt(ir,iat)=irt0(iat)
          ENDDO
          symmr%ftau(1,ir) = ft(1)
          symmr%ftau(2,ir) = ft(2)
          symmr%ftau(3,ir) = ft(3)
       ENDIF
    ENDDO
    nsym=0
    DO ir=1,symmi%nrot
       IF (sym(ir))nsym=nsym+1
    ENDDO
    IF (nsym.EQ.symmi%nrot) GOTO 300
    IF (parm%alat.GT.2._real_8*symmr%deltasym) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' SYMTAU! FAILED TO SYMMETRIZE ATOMIC COORDINATES '
          CALL stopgm('SYMTAU',' ERROR ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    GOTO 200
300 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,'(A,G10.3)')&
         ' SYMTAU! TOLERANCE NEEDED FOR FULL SYMMETRY      :',&
         symmr%deltasym
    CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,rau(1,1),1)
    CALL symvec(tau0)
    disp=0.0_real_8
    CALL daxpy(3*maxsys%nax*maxsys%nsx,-1.0_real_8,tau0(1,1,1),1,rau(1,1),1)
    imax=idamax(3*maxsys%nax*maxsys%nsx,rau(1,1),1)
    disp=ABS(dgive(rau(1,1),imax))
    IF (paral%io_parent)&
         WRITE(6,'(A,G10.3)')&
         ' SYMTAU! LARGEST DISPLACEMENT FOR SYMMETRIZATION :',disp
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          tau0(1,ia,is)=tau0(1,ia,is)+cm(1)
          tau0(2,ia,is)=tau0(2,ia,is)+cm(2)
          tau0(3,ia,is)=tau0(3,ia,is)+cm(3)
       ENDDO
    ENDDO
    DEALLOCATE(xau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ityp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(irt0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE symtau
  ! ==================================================================
  SUBROUTINE gentau(tau0)
    ! ==--------------------------------------------------------------==
    ! == Generate atomic coordinates from unique atoms                ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gentau'

    INTEGER                                  :: ia, iat, iatn, ib, ierr, ir, &
                                                is, k
    REAL(real_8)                             :: d1, d2, d3, dd
    REAL(real_8), ALLOCATABLE                :: rau(:,:), xau(:,:)

    IF (.NOT.symmt%tgenc) RETURN
    IF (symmi%indpg.EQ.0) RETURN
    IF (paral%io_parent) THEN
       WRITE(6,'(/,A,A)')' GENERATE ATOMIC COORDINATES FROM ',&
            'SYMMETRY UNIQUE ONES'
    ENDIF
    ALLOCATE(xau(3,maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rau(3,maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! LOOP OVER SPECIES AND SYMMETRY UNIQUE ATOMS
    iat=0
    DO is=1,ions1%nsp
       ! ATOMIC COORDINATES IN DIRECT LATTICE VECTOR BASIS
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,tau0(1,ia,is),1,0.0_real_8,&
               xau(1,ia),1)
          DO k=1,3
             xau(k,ia)=xau(k,ia)-NINT(xau(k,ia))
          ENDDO
       ENDDO
       iat=iun(is)
       DO ir=1,symmi%nrot
          DO ia=1,iat
             DO k=1,3
                rau(k,ia)=symmr%xtable(k,1,ir)*xau(1,ia)+&
                     symmr%xtable(k,2,ir)*xau(2,ia)+&
                     symmr%xtable(k,3,ir)*xau(3,ia)
                rau(k,ia)=rau(k,ia)-NINT(rau(k,iat))
             ENDDO
          ENDDO
          iatn=iat
          DO ia=1,iat
             DO ib=1,iatn
                d1=ABS(rau(1,ia)-xau(1,ib))
                d2=ABS(rau(2,ia)-xau(2,ib))
                d3=ABS(rau(3,ia)-xau(3,ib))
                dd=SQRT(d1*d1+d2*d2+d3*d3)
                IF (dd.LT.symmr%deltasym) GOTO 100
             ENDDO
             ! NEW ATOM
             iatn=iatn+1
             CALL dcopy(3,rau(1,ia),1,xau(1,iatn),1)
100          CONTINUE
          ENDDO
          iat=iatn
       ENDDO
       ! CHECK THAT ALL ATOMS HAVE BEEN GENERATED
       IF (ions0%na(is).NE.iat) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,T61,I5)')&
                  ' GENTAU| INCONSISTENCY FOR SPECIES:',is
             IF (paral%io_parent)&
                  WRITE(6,'(A,T51,I5,A5,I5)')&
                  ' GENTAU| NUMBER OF ATOMS:',iat,' NOT ',ions0%na(is)
             CALL stopgm('GENTAU',' CHECK INPUT ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ATOMIC COORDINATES IN CARTESIAN COORDINATES
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,xau(1,ia),1,0.0_real_8,tau0(1,ia,is),&
               1)
       ENDDO
    ENDDO
    DEALLOCATE(xau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gentau
  ! ==================================================================
  SUBROUTINE updatsym(tau0)
    ! ==--------------------------------------------------------------==
    ! == UPDATE SYMMETRY OPERATIONS                                   ==
    ! == (GEOMETRY OPTIMISATION OR MOLECULAR DYNAMICS)                ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   TAU0(1:3,maxsys%nax,maxsys%nsx) Atomic coordinates                       ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: i

    IF (symmt%tpgauto) CALL givesym(tau0)
    ! Initialisation for all processors
    IF (symmi%indpg.EQ.0) THEN
       symmi%nrot=1
       symmi%inversion=0
       symmt%torigin=.FALSE.
       DO i=1,3
          symmr%origin(i)=0._real_8
       ENDDO
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL chksym(tau0)
    ! Construct table multiplication.
    ! IRT generates 'non-allocated usage' error when INDPG=0, should not be passed
    CALL multtb(symmi%indpg,symmi%nrot,symmi%ntvec,ions1%nat,symmr%xtable,symmi%inve,symmi%multab,.FALSE.)
    ! Number of equivalent shell (0 for serial version)
    ! Use for symmetrisation of density (force call to BUILDSHELL)
    numshel=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updatsym
  ! ==================================================================
  SUBROUTINE checksym(nat,ityp,xau,rau,ft,nodupli,sym,irt,isc,delta)
    ! ==--------------------------------------------------------------==
    ! == CHECK FOR EACH ATOM IA=1:NAT THAT THE IMAGE ATOM RAU(3,1:NAT)==
    ! == GIVEN BY THE IR ROTATION OPERATION WITH THE TRANSLATION      ==
    ! == VECTOR FT(3) CORRESPONDS TO AN ATOM                          ==
    ! == BUILD ALSO THE ATOMIC TRANSFORMATION TABLE                   ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   NAT         Number of atoms                                ==
    ! ==   ITYP(NAT)   Type of atoms                                  ==
    ! ==   XAU(3,NAT)  Atomic coordinates                             ==
    ! ==   RAU(3,NAT)  Images from pure rotation of atomic coordinates==
    ! ==   FT(3)       Translation vector                             ==
    ! ==   NODUPLI     .TRUE. the cell is duplicated, we can speed up ==
    ! ==   DELTA       Precision required                             ==
    ! == OUTPUT:                                                      ==
    ! ==   SYM         .TRUE. if everything is o.k.                   ==
    ! ==   IRT(NAT)    Atomic transformation table                    ==
    ! ==   ISC(NAT)    Scratch array                                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nat, ityp(nat)
    REAL(real_8)                             :: xau(3,nat), rau(3,nat), ft(3)
    LOGICAL                                  :: nodupli, sym
    INTEGER                                  :: irt(nat), isc(nat)
    REAL(real_8)                             :: delta

    INTEGER                                  :: ia, ib

    DO ia=1,nat
       isc(ia)=0
    ENDDO
    ext_loop: DO ia = 1, nat
       sym = .FALSE.
       DO ib = 1, nat
          IF (ityp(ia).EQ.ityp(ib).AND.isc(ib).EQ.0) THEN
             sym = eqvect(rau(1,ia),xau(1,ib),ft,delta)
             IF (sym) THEN
                ! THE ROTATED ATOM DOES COINCIDE WITH ONE OF THE LIKE ATOMS
                ! KEEP TRACK OF WHICH ATOM THE ROTATED ATOM COINCIDES WITH
                irt(ia)=ib
                IF (nodupli) isc(ib)=1
                CYCLE ext_loop
             ENDIF
          ENDIF
       ENDDO
       ! THE ROTATED ATOM DOES NOT COINCIDE WITH ANY OF THE LIKE ATOMS
       ! S(IR) + FT IS NOT A SYMMETRY OPERATION
       RETURN
    END DO ext_loop
    ! S(IR) + FT IS A SYMMETRY OPERATION
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE checksym
  ! ==================================================================
  FUNCTION eqvect(x,y,f,delta)
    ! ==--------------------------------------------------------------==
    ! == TEST IF Y = X + F (VECTOR)                                   ==
    ! == DELTA GIVES THE REQUIRED PRECISION                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x(3), y(3), f(3), delta
    LOGICAL                                  :: eqvect

    eqvect = ABS(x(1)-y(1)+f(1)-NINT(x(1)-y(1)+f(1))).LT. delta  .AND.&
         ABS(x(2)-y(2)+f(2)-NINT(x(2)-y(2)+f(2))).LT. delta  .AND.&
         ABS(x(3)-y(3)+f(3)-NINT(x(3)-y(3)+f(3))).LT. delta
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION eqvect
  ! ==================================================================

END MODULE chksym_utils
