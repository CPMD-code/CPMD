#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#define __COLLAPSE4
#else
#define __COLLAPSE2 collapse(2)
#define __COLLAPSE4 collapse(4)
#endif
!     ==================================================================
      SUBROUTINE mm_translate_qmmm(TAU,C0,CM,NSTATE)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
!     ==--------------------------------------------------------------==
!cmb - Revised on 25 October 2013
      use system 
      use parac
      use cppt
      use gvec
      use isos
      use mm_dimmod
      use geq0mod
      use elct      ! for N
      use mm_extrap ! for COLD()
      IMPLICIT NONE 

      integer NSTATE
CEHR[
      COMPLEX*16 C0(nkpt%NGWK,NSTATE,*) !changed for BSYMM 
      COMPLEX*16 CM(nkpt%NGWK,NSTATE,*) !changed for BSYMM
CEHR]
      real*8 TAU(3,maxsys%nax,maxsys%nsx)

      integer IG,IS,IH,IK
      integer NH1,NH2,NH3
      real*8 Gt1,Gt2,Gt3,rI,rJ,rK,gdotr,c_trans(3)
      complex*16, allocatable ::  eigr_l(:)

      if(.not.isos1%tcent)return
      CALL  mm_center(TAU,c_trans,.TRUE.)
      clsaabox%mm_c_trans(1)= clsaabox%mm_c_trans(1)+c_trans(1)
      clsaabox%mm_c_trans(2)= clsaabox%mm_c_trans(2)+c_trans(2)
      clsaabox%mm_c_trans(3)= clsaabox%mm_c_trans(3)+c_trans(3)

      if(.not.paral%qmnode) return

      allocate (eigr_l(ncpw%ngw))
      NH1=spar%nr1s/2+1
      NH2=spar%nr2s/2+1
      NH3=spar%nr3s/2+1 

!$OMP parallel do private(IG,rI,rJ,rK,Gt1,Gt2,Gt3,gdotr)
#ifdef __SR11000
*poption parallel, tlocal(IG,rI,rJ,rK,Gt1,Gt2,Gt3,gdotr)
#endif
      DO IG=1,ncpw%ngw
        rI=DBLE(INYH(1,IG)-NH1)
        rJ=DBLE(INYH(2,IG)-NH2)
        rK=DBLE(INYH(3,IG)-NH3)
        Gt1=rI*gvec_com%b1(1)+rJ*gvec_com%b2(1)+rK*gvec_com%b3(1)
        Gt2=rI*gvec_com%b1(2)+rJ*gvec_com%b2(2)+rK*gvec_com%b3(2)
        Gt3=rI*gvec_com%b1(3)+rJ*gvec_com%b2(3)+rK*gvec_com%b3(3)
        gdotr=parm%tpiba*(Gt1*c_trans(1)
     .              +Gt2*c_trans(2)
     .              +Gt3*c_trans(3))
        eigr_l(IG)=DCMPLX(DCOS(gdotr),-DSIN(gdotr)) 
      ENDDO
!
!     If a Broken symmetry run..the HS wf and its velocities are
!     scaled simultaneously (NN)
      IF(CNTL%BSYMM)THEN
#ifdef __SR11000
*poption parallel, tlocal(IS,IG)
#endif
!$OMP parallel do private(IS,IG) __COLLAPSE2
        DO IS=1,NSTATE
          DO IG=1,ncpw%ngw
            C0(IG,IS,1)=C0(IG,IS,1)*eigr_l(IG) !BS
            C0(IG,IS,2)=C0(IG,IS,2)*eigr_l(IG) !HS
            CM(IG,IS,1)=CM(IG,IS,1)*eigr_l(IG) !BS
            CM(IG,IS,2)=CM(IG,IS,2)*eigr_l(IG) !HS
          ENDDO
        ENDDO
      ELSE
#ifdef __SR11000
*poption parallel, tlocal(IS,IG)
#endif
!$OMP parallel do private(IS,IG) __COLLAPSE2
        DO IS=1,NSTATE
          DO IG=1,ncpw%ngw
            C0(IG,IS,1)=C0(IG,IS,1)*eigr_l(IG)
            CM(IG,IS,1)=CM(IG,IS,1)*eigr_l(IG)
          ENDDO
        ENDDO
CEHR[
        IF (CNTL%TMDEH) THEN
!$OMP parallel do private(IS,IG)
         DO IS=1,NSTATE
          DO IG=1,ncpw%ngw
            C0(IG+ncpw%ngw,IS,1)=C0(IG+ncpw%ngw,IS,1)*DCONJG(eigr_l(IG))
            CM(IG+ncpw%ngw,IS,1)=CM(IG+ncpw%ngw,IS,1)*DCONJG(eigr_l(IG))
          ENDDO
          IF(GEQ0) C0(ncpw%ngw+1,IS,1)=DCMPLX(0.D0,0.D0)
         ENDDO
        ENDIF
CEHR]
      END IF

      IF(CNTL%TMDBO.AND.CNTL%TEXTRAP) THEN
C     if using wavefunction extrapolation, apply the same
C     shift also to the wavefunction history.
#ifdef __SR11000
*poption parallel, tlocal(IH,IK,IS,IG)
#endif
!$OMP parallel do private(IH,IK,IS,IG) __COLLAPSE4
        DO IH=1,cnti%MEXTRA
          DO IK=1,nkpt%nkpnt
            DO IS=1,NSTATE
              DO IG=1,ncpw%ngw
                COLD(IG,IS,IK,IH)=COLD(IG,IS,IK,IH)*eigr_l(IG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      deallocate(eigr_l)

      RETURN
      END
