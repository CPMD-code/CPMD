MODULE calc_pij_utils
  USE cppt,                            ONLY: gk
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk
  USE kpts,                            ONLY: tkpts
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_pij
  PUBLIC :: calc_pi

CONTAINS

  ! ==================================================================
  SUBROUTINE calc_pij(c0,c1,cpx,cpy,cpz,ikind)
    ! ==--------------------------------------------------------------==
    ! == GIVEN C0 AND C1 IN THE PLANE WAVE REPRESENTATION, THIS       ==
    ! == ROUTINE RETURNS  <C0|P|C1>                                   ==
    ! == IN UNITS OF 2 PI/A                                           ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk), c1(nkpt%ngwk), &
                                                cpx, cpy, cpz
    INTEGER                                  :: ikind

    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: caux, caux2
    INTEGER                                  :: ig
    REAL(real_8)                             :: px, py, pz

! ==--------------------------------------------------------------==

    IF (tkpts%tkpnt) THEN
       cpx = zzero
       cpy = zzero
       cpz = zzero
       DO ig=1,ncpw%ngw
          caux=c1(ig)*CONJG(c0(ig))
          caux2=c1(ig+ncpw%ngw)*CONJG(c0(ig+ncpw%ngw))
          cpx=cpx+rk(1,ikind)*(caux+caux2)&
               +GK(1,IG)*(CAUX-CAUX2)
          cpy=cpy+rk(2,ikind)*(caux+caux2)&
               +GK(2,IG)*(CAUX-CAUX2)
          cpz=cpz+rk(3,ikind)*(caux+caux2)&
               +GK(3,IG)*(CAUX-CAUX2)
       ENDDO
    ELSE
       px = 0._real_8
       py = 0._real_8
       pz = 0._real_8
       DO ig=1,ncpw%ngw
          caux=CONJG(c0(ig))*c1(ig)
          px=px+gk(1,ig)*AIMAG(caux)
          py=py+gk(2,ig)*AIMAG(caux)
          pz=pz+gk(3,ig)*AIMAG(caux)
       ENDDO
       cpx=CMPLX(px*2._real_8,0._real_8,kind=real_8)
       cpy=CMPLX(py*2._real_8,0._real_8,kind=real_8)
       cpz=CMPLX(pz*2._real_8,0._real_8,kind=real_8)
    ENDIF
    RETURN
  END SUBROUTINE calc_pij
  ! ==================================================================
  SUBROUTINE calc_pi(c0,c1,j,ikind)
    ! ==--------------------------------------------------------------==
    ! == GIVEN C0 IN THE PLANE WAVE REPRESENTATION, THIS              ==
    ! == ROUTINE RETURNS  C1=P_J.C0=-id/dx C0                         ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk), c1(nkpt%ngwk)
    INTEGER                                  :: j, ikind

    INTEGER                                  :: ig

    IF (tkpts%tkpnt) THEN
       DO ig=1,ncpw%ngw
          c1(ig)    =(rk(j,ikind)+gk(j,ig))*c0(ig)*parm%tpiba
          c1(ig+ncpw%ngw)=(rk(j,ikind)-gk(j,ig))*c0(ig+ncpw%ngw)*parm%tpiba
       ENDDO
       IF (geq0)c1(1+ncpw%ngw)=CMPLX(0._real_8,0._real_8,kind=real_8)
    ELSE
       DO ig=2,ncpw%ngw
          c1(ig)=gk(j,ig)*c0(ig)*parm%tpiba
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_pi
  ! ==================================================================

END MODULE calc_pij_utils
