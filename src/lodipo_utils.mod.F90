MODULE lodipo_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: gk,&
                                             nzh
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE lodp,                            ONLY: &
       chrld, dmomlo, numbld, rcc, xmaxld, xminld, ymaxld, yminld, zmaxld, &
       zminld
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lodipo

CONTAINS

  ! IF   * SUBROUTINE FOR CALCULATION OF LOCAL DIPOLE MOMENTS *
  ! IF   * (I. FRANK, OCT. 94)
  ! ==================================================================
  SUBROUTINE lodipo(eirop,rhog)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCAL DIPOLE MOMENTS;                             ==
    ! == DOES THE INTEGRATION USING CARTESIAN COORDINATES;            ==
    ! == FOR EXPLANATION OF SOME VARIABLES SEE proppt.f               ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                rhog(fpar%nnr1)

    COMPLEX(real_8), PARAMETER               :: cone = (1._real_8,0._real_8) 

    COMPLEX(real_8)                          :: cg, cg1, cg2, cg3, dr(3), &
                                                emax(3), emin(3), gmax(3), &
                                                gmin(3), rdr(3), rsum(3), sum
    INTEGER                                  :: i, ig, ipm, isub, j, k
    REAL(real_8)                             :: g(3), g1, g2, g3, rmax(3), &
                                                rmin(3)

! ==--------------------------------------------------------------==

    CALL tiset('    LODIPO',isub)
    ! Loop over the NUMBLD areas given in the input file
    DO i=1,numbld
       rmin(1)=xminld(i)
       rmax(1)=xmaxld(i)
       rmin(2)=yminld(i)
       rmax(2)=ymaxld(i)
       rmin(3)=zminld(i)
       rmax(3)=zmaxld(i)
       ! RCC: coordinates of the center of charge
       DO k=1,3
          rcc(k,i)=(rmin(k)+rmax(k))/2
          rsum(k)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
       sum=CMPLX(0._real_8,0._real_8,kind=real_8)
       ! Calculate the center of positive charge
       ! sum over all components of the wave function
       DO ig=1,ncpw%nhg
          ! sum twice (for positive and negative g values)
          DO ipm=-1,1,2
             ! do not sum twice if |g|=0
             IF (ABS(gk(1,ig)).GT.0.001_real_8&
                  .OR.ABS(gk(2,ig)).GT.0.001_real_8&
                  .OR.ABS(gk(3,ig)).GT.0.001_real_8&
                  .OR.ipm.GT.0) THEN
                cg=CMPLX(REAL(eirop(ig)),&
                     AIMAG(eirop(ig))*ipm,kind=real_8)
                DO k=1,3
                   g(k)=gk(k,ig)*parm%tpiba*ipm
                   gmin(k)=uimag*g(k)*rmin(k)
                   gmax(k)=uimag*g(k)*rmax(k)
                   emin(k)=EXP(gmin(k))
                   emax(k)=EXP(gmax(k))
                   IF (ABS(gk(k,ig)).GT.0.001_real_8) THEN
                      rdr(k)=-((gmax(k)-cone)*emax(k)&
                           -(gmin(k)-cone)*emin(k))&
                           /g(k)/g(k)
                      dr(k)=-uimag*(emax(k)-emin(k))/g(k)
                   ELSE
                      rdr(k)=(rmax(k)*rmax(k)&
                           -rmin(k)*rmin(k))/2
                      dr(k)=rmax(k)-rmin(k)
                   ENDIF
                ENDDO
                rsum(1)=rsum(1)-cg*rdr(1)*dr(2)*dr(3)
                rsum(2)=rsum(2)-cg*dr(1)*rdr(2)*dr(3)
                rsum(3)=rsum(3)-cg*dr(1)*dr(2)*rdr(3)
                sum=sum-cg*dr(1)*dr(2)*dr(3)
             ENDIF
          ENDDO
       ENDDO
       CALL mp_sum(sum,parai%allgrp)
       CALL mp_sum(rsum,3,parai%allgrp)
       IF (ABS(REAL(sum)).GT.0.001_real_8) THEN
          DO k=1,3
             rcc(k,i)=REAL(rsum(k))/REAL(sum)
          ENDDO
       ENDIF
       chrld(i)=0._real_8
       DO ig=1,ncpw%nhg
          g1=gk(1,ig)*parm%tpiba
          g2=gk(2,ig)*parm%tpiba
          g3=gk(3,ig)*parm%tpiba
          IF (ABS(g1).LT.0.001_real_8) THEN
             cg1=rmax(1)-rmin(1)
          ELSE
             cg1=-uimag/g1*(EXP(uimag*rmax(1)*g1)-EXP(uimag*rmin(1)*g1))
          ENDIF
          IF (ABS(g2).LT.0.001_real_8) THEN
             cg2=rmax(2)-rmin(2)
          ELSE
             cg2=-uimag/g2*(EXP(uimag*rmax(2)*g2)-EXP(uimag*rmin(2)*g2))
          ENDIF
          IF (ABS(g3).LT.0.001_real_8) THEN
             cg3=rmax(3)-rmin(3)
          ELSE
             cg3=-uimag/g3*(EXP(uimag*rmax(3)*g3)-EXP(uimag*rmin(3)*g3))
          ENDIF
          chrld(i)=chrld(i)+2._real_8*REAL(cg1*cg2*cg3*rhog(nzh(ig)))
       ENDDO
       IF (geq0) THEN
          chrld(i)=chrld(i)-REAL(rhog(nzh(1)))*&
               (rmax(1)-rmin(1))*(rmax(2)-rmin(2))*(rmax(3)-rmin(3))
       ENDIF
       CALL mp_sum(chrld(i),parai%allgrp)
       ! start calculating the dipole moment
       ! DMOMLO(J,I): component J of the dipole moment for area I
       DO j=1,3
          dmomlo(j,i)=0._real_8
       ENDDO
       ! calculate integration range relative to center of charge
       ! sum over all components of the wave function
       DO ig=1,ncpw%nhg
          ! sum twice (for positive and negative g values)
          DO ipm=-1,1,2
             ! do not sum twice if |g|=0
             IF (ABS(gk(1,ig)).GT.0.001_real_8&
                  .OR.ABS(gk(2,ig)).GT.0.001_real_8&
                  .OR.ABS(gk(3,ig)).GT.0.001_real_8&
                  .OR.ipm.GT.0) THEN
                cg=CMPLX(REAL(rhog(nzh(ig))+eirop(ig)),&
                     AIMAG(rhog(nzh(ig))+eirop(ig))*ipm,kind=real_8)
                DO k=1,3
                   g(k)=gk(k,ig)*parm%tpiba*ipm
                   gmin(k)=uimag*g(k)*rmin(k)
                   gmax(k)=uimag*g(k)*rmax(k)
                   emin(k)=EXP(gmin(k))
                   emax(k)=EXP(gmax(k))
                   IF (ABS(gk(k,ig)).GT.0.001_real_8) THEN
                      rdr(k)=-((gmax(k)-cone)*emax(k)&
                           -(gmin(k)-cone)*emin(k))&
                           /g(k)/g(k)
                      dr(k)=-uimag*(emax(k)-emin(k))/g(k)
                   ELSE
                      rdr(k)=(rmax(k)*rmax(k)&
                           -rmin(k)*rmin(k))/2
                      dr(k)=rmax(k)-rmin(k)
                   ENDIF
                ENDDO
                dmomlo(1,i)=dmomlo(1,i)&
                     -REAL(cg*(rdr(1)-dr(1)*rcc(1,i))*dr(2)*dr(3))
                dmomlo(2,i)=dmomlo(2,i)&
                     -REAL(cg*dr(1)*(rdr(2)-dr(2)*rcc(2,i))*dr(3))
                dmomlo(3,i)=dmomlo(3,i)&
                     -REAL(cg*dr(1)*dr(2)*(rdr(3)-dr(3)*rcc(3,i)))
             ENDIF
          ENDDO
       ENDDO
       CALL mp_sum(dmomlo(:,i),3,parai%allgrp)
    ENDDO
    CALL tihalt('    LODIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lodipo
  ! ==================================================================

END MODULE lodipo_utils
