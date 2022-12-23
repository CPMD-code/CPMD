MODULE znum_mat_utils
  USE cppt,                            ONLY: gk
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: znum_mat
  PUBLIC :: give_scr_znum_mat

CONTAINS

  ! ==================================================================
  SUBROUTINE znum_mat(c0,c2,nstate,ddmat,index)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(2*ncpw%ngw,*)
    COMPLEX(real_8), TARGET                  :: c2(2*ncpw%ngw,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ddmat(nstate,nstate)
    INTEGER                                  :: index

    CHARACTER(*), PARAMETER                  :: procedureN = 'znum_mat'

    INTEGER                                  :: i, ierr, ig, ii, ip1, ip2, &
                                                ip3, nh(3)
    REAL(real_8)                             :: im, numsec(3), prod, re, x0(3)
    REAL(real_8), ALLOCATABLE                :: eii(:), eir(:), sc0(:,:)

    CALL zeroing(ddmat)!,nstate*nstate)

    ip1 = 1
    ip2 = ip1 + 2*2*ncpw%ngw
    ip3 = ip2 + 2*ncpw%ngw

    ALLOCATE(sc0(2, 2*ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eir(2*ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eii(2*ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    numsec(1) = REAL(spar%nr1s,kind=real_8)
    numsec(2) = REAL(spar%nr2s,kind=real_8)
    numsec(3) = REAL(spar%nr3s,kind=real_8)

    nh(1)=spar%nr1s/2+1
    nh(2)=spar%nr2s/2+1
    nh(3)=spar%nr3s/2+1


    ! .... MINIMUM R-VECTOR in the grid in along the axes INDEX 

    ! X0(1) = ipos(1)*A1(1)/NUMSEC(1)+
    ! &        ipos(2)*A2(1)/NUMSEC(2)+ipos(3)*A3(1)/NUMSEC(3)
    ! X0(2) = ipos(1)*A1(2)/NUMSEC(1)+
    ! &        ipos(2)*A2(2)/NUMSEC(2)+ipos(3)*A3(2)/NUMSEC(3)
    ! X0(3) = ipos(1)*A1(3)/NUMSEC(1)+
    ! &        ipos(2)*A2(3)/NUMSEC(2)+ipos(3)*A3(3)/NUMSEC(3)

    x0(1) = metr_com%ht(index,1)/numsec(index)
    x0(2) = metr_com%ht(index,2)/numsec(index)
    x0(3) = metr_com%ht(index,3)/numsec(index)
    IF (paral%io_parent)&
         WRITE(6,'("X0",3f12.6)') x0(1),x0(2),x0(3)

    DO ig=1,ncpw%ngw

       prod = gk(1,ig)*x0(1)+gk(2,ig)*x0(2)+gk(3,ig)*x0(3)
       prod = parm%tpiba*prod
       eir(ig)     =  COS(prod)
       eii(ig)     =  SIN(prod)

       ! EIR(IG)     =  cos(GK(INDEX,IG)*X0(INDEX))
       ! EII(IG)     =  sin(GK(INDEX,IG)*X0(INDEX))

       eir(ig+ncpw%ngw) =  eir(ig)
       eii(ig+ncpw%ngw) = -eii(ig)

       ! write(6,'(2(I8,2f12.6))') IG,EIR(IG),EII(IG),
       ! &                          IG+NGW,EIR(IG+NGW),EII(IG+NGW)
       ! write(6,'(8x,3f12.6)') GK(1,IG),GK(3,IG),GK(2,IG)


    ENDDO
    ! stop

    DO i= 1,nstate ! NST12(MEPOS,1),NST12(MEPOS,2)

       ii=i! -NST12(MEPOS,1)+1

       CALL dcopy(4*ncpw%ngw,c0(1,ii),1,sc0,1)

       DO ig=1,2*ncpw%ngw
          re = eir(ig)*sc0(1,ig)-eii(ig)*sc0(2,ig)
          im = eir(ig)*sc0(2,ig)+eii(ig)*sc0(1,ig)
          c2(ig,ii)=CMPLX(re,im,kind=real_8)
       ENDDO
       IF (geq0)&
            c2(1+ncpw%ngw,ii)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ! DO IG=1,10
       ! write(6,'(a,2f12.6)') 'EIR',EIR(IG),EII(IG)
       ! write(6,'(i8,4f12.6)') IG,C0(IG,II),C2F(IG,II)
       ! write(6,'(i8,4f12.6)') IG+NGW,C0(IG+NGW,II),C2F(IG+NGW,II)
       ! enddo
    ENDDO

    IF (cntl%tlsd) THEN
       CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,2*ncpw%ngw,&
            c2,2*ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
       CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),&
            c0(1,spin_mod%nsup+1),2*ncpw%ngw,c2(1,spin_mod%nsup+1),2*ncpw%ngw,&
            CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
    ELSE


       CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),&
            c0,2*ncpw%ngw,c2,2*ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)


       ! write(6,'(2I4,2f12.6,2x,2f12.6)') 18,16,DDMAT(18,16),DDMAT(16,18)

       ! stop
    ENDIF
    CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)


    ! DO I = 10,20
    ! IF(PARENT) WRITE(6,'(I4,11(2x,2f7.4))') I,( DDMAT(I,II),II=10,20)
    ! ENDDO
    ! 
    ! stop
    ! ==--------------------------------------------------------------==
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eir,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eii,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE znum_mat
  ! ==================================================================
  SUBROUTINE give_scr_znum_mat(lznum,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lznum
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lznum=2*nstate*nstate+4*ncpw%ngw+4*ncpw%ngw+1
    tag='2*NSTATE*NSTATE+2*NGW+4*NGW'
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_znum_mat
  ! ==================================================================

END MODULE znum_mat_utils
