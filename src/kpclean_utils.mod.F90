MODULE kpclean_utils
  USE kinds,                           ONLY: real_8
  USE sphe,                            ONLY: maskgw,&
                                             maskl,&
                                             tsphere
  USE system,                          ONLY: ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: c_clean
  PUBLIC :: r_clean

CONTAINS

  ! ==================================================================
  SUBROUTINE c_clean(c,nstate,ik)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c(ncpw%ngw,2,nstate)
    INTEGER                                  :: ik

    INTEGER                                  :: ig, igg, is, nn

    IF (tsphere) THEN
       nn=ncpw%ngw-maskl
       DO is=1,nstate
          DO ig=nn+1,ncpw%ngw
             igg=ig-nn
             IF (maskgw(1,igg,ik).LT.0.5_real_8) c(ig,1,is)=CMPLX(0._real_8,0._real_8,kind=real_8)
          ENDDO
          DO ig=nn+1,ncpw%ngw
             igg=ig-nn
             IF (maskgw(2,igg,ik).LT.0.5_real_8) c(ig,2,is)=CMPLX(0._real_8,0._real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE c_clean
  ! ==================================================================
  SUBROUTINE r_clean(t,nstate,ik)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: t(ncpw%ngw,2,nstate)
    INTEGER                                  :: ik

    INTEGER                                  :: ig, igg, is, nn

    IF (tsphere) THEN
       nn=ncpw%ngw-maskl
       DO is=1,nstate
          DO ig=nn+1,ncpw%ngw
             igg=ig-nn
             IF (maskgw(1,igg,ik).LT.0.5_real_8) t(ig,1,is)=0._real_8
          ENDDO
          DO ig=nn+1,ncpw%ngw
             igg=ig-nn
             IF (maskgw(2,igg,ik).LT.0.5_real_8) t(ig,2,is)=0._real_8
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE r_clean
  ! ==================================================================


END MODULE kpclean_utils
