MODULE kin_energy_utils
  USE cppt,                            ONLY: hg
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE prcp,                            ONLY: prcp_com
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kin_energy

CONTAINS

  ! ==================================================================
  SUBROUTINE kin_energy(c0,nstate,rsum)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE KINETIC ENERGY EKIN. IT IS DONE IN RECIPROCAL SPACE     ==
    ! ==  WHERE THE ASSOCIATED OPERATORS ARE DIAGONAL.                ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rsum

    REAL(real_8), PARAMETER                  :: deltakin = 1.e-10_real_8 

    INTEGER                                  :: i, ig, is1, isub
    REAL(real_8)                             :: arg, g2, ima, imb, ra, rb, &
                                                sk1, sk2, xkin, xskin

#if defined(__SR8000)
    REAL(real_8), ALLOCATABLE :: g2w(:)

    INTEGER, SAVE :: ifirst=0
#endif
    ! ==--------------------------------------------------------------==
    CALL tiset('KIN_ENERGY',isub)
    ! ==--------------------------------------------------------------==
#if defined(__SR8000)
    IF (ifirst.EQ.0) THEN
       ALLOCATE(g2w(ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    ! Accumulate the charge and kinetic energy
    rsum=0._real_8
    xkin=0._real_8
    !$omp parallel do private(I,SK1,XSKIN,IS1,ARG,G2,IG) &
    !$omp  reduction(+:RSUM,XKIN)
    DO i=1,nstate
       IF (crge%f(i,1).NE.0._real_8) THEN! TODO check F(I,1) === F(I)
          rsum=rsum+crge%f(i,1)*dotp(ncpw%ngw,c0(:,i),c0(:,i))
          sk1=0.0_real_8
          IF (prcp_com%akin.GT.deltakin) THEN
             xskin=1._real_8/prcp_com%gskin
             is1=1
             IF (geq0) THEN
                is1=2
                arg=-prcp_com%gckin*xskin
                g2=0.5_real_8*prcp_com%gakin*(1._real_8+cp_erf(arg))
                sk1=sk1+REAL(g2*CONJG(c0(1,i))*c0(1,i))
             ENDIF
#if defined(__SR8000)
             !poption parallel
             !poption tlocal(ARG)
             DO ig=is1,ngw
                arg=(hg(ig)-prcp_com%gckin)*xskin
                g2w(ig)=hg(ig)+prcp_com%gakin*(1._real_8+cp_erf(arg))
             ENDDO
#ifdef __SR8000
             !poption parallel, tlocal(IG), psum(SK1)
#endif 
             DO ig=is1,ngw
                sk1=sk1+REAL(g2w(ig)*CONJG(c0(ig,i))*c0(ig,i))
             ENDDO
#else
             DO ig=is1,ncpw%ngw
                arg=(hg(ig)-prcp_com%gckin)*xskin
                g2=hg(ig)+prcp_com%gakin*(1._real_8+cp_erf(arg))
                sk1=sk1+REAL(g2*CONJG(c0(ig,i))*c0(ig,i))
             ENDDO
#endif
          ELSE
#ifdef __SR8000
             !poption parallel, tlocal(IG), psum(SK1)
#endif 
             DO ig=1,ncpw%ngw
                sk1=sk1+REAL(hg(ig)*CONJG(c0(ig,i))*c0(ig,i))
             ENDDO
          ENDIF
          xkin=xkin+crge%f(i,1)*sk1
       ENDIF
    ENDDO
    ener_com%ekin=xkin*parm%tpiba2
    ! Other kinetic energies for CAS22 method
    IF (lspin2%tlse .AND. (lspin2%tcas22.OR.lspin2%tpenal)) THEN
       IF (prcp_com%akin.GT.deltakin) CALL stopgm('RHOOFR',&
            'CAS22 and AKIN not implemented',& 
            __LINE__,__FILE__)
       sk1=0._real_8
       DO ig=1,ncpw%ngw
          ra=REAL(c0(ig,clsd%ialpha))
          rb=REAL(c0(ig,clsd%ibeta))
          ima=AIMAG(c0(ig,clsd%ialpha))
          imb=AIMAG(c0(ig,clsd%ibeta))
          sk1=sk1+hg(ig)*(ra*rb+ima*imb)
       ENDDO
       ener_c%ekin_ab = sk1 * parm%tpiba2
    ENDIF
    IF (lspin2%tlse .AND. lspin2%tcas22) THEN
       sk1=0._real_8
       sk2=0._real_8
       DO ig=1,ncpw%ngw
          sk1=sk1+REAL(hg(ig)*CONJG(c0(ig,clsd%ialpha))*c0(ig,clsd%ialpha))
          sk2=sk2+REAL(hg(ig)*CONJG(c0(ig,clsd%ibeta))*c0(ig,clsd%ibeta))
       ENDDO
       ener_c%ekin_a = ener_com%ekin - parm%tpiba2 * ( sk2 - sk1 )
       ener_c%ekin_2 = ener_com%ekin + parm%tpiba2 * ( sk2 - sk1 )
    ENDIF
    ! 
    CALL tihalt('KIN_ENERGY',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kin_energy
  ! ==================================================================

END MODULE kin_energy_utils
