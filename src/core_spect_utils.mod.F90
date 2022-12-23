MODULE core_spect_utils
  USE adat,                            ONLY: elem
  USE atwf,                            ONLY: atrg,&
                                             atwr
  USE cnst,                            ONLY: fpi,&
                                             ry
  USE cores,                           ONLY: core_atwfr,&
                                             core_c0,&
                                             core_cat,&
                                             coresi,&
                                             coresr
  USE cppt,                            ONLY: gk,&
                                             hg
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk
  USE kpts,                            ONLY: tkpts
  USE lsfbtr_utils,                    ONLY: lsfbtr
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE qspl,                            ONLY: ggng,&
                                             nsplpo
  USE recpnew_utils,                   ONLY: tgrid
  USE setbasis_utils,                  ONLY: storps
  USE sfac,                            ONLY: eigr
  USE sphe,                            ONLY: gcutka
  USE system,                          ONLY: iatpt,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE ylmr_utils,                      ONLY: ylmr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: core_spectra

CONTAINS

  ! ==================================================================
  SUBROUTINE core_spectra(tau0,c0,we,f,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == Compute the X-ray core spectra within the framework of       ==
    ! == M. Cavalleri et al. J. Chem. Phys. 121, 10065 (2004)         ==
    ! == Present release: Tokai-mura/Tsukuba, 14 June 2005            ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                f(nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'core_spectra'

    COMPLEX(real_8)                          :: caux, ci, csumx, csumy, csumz
    COMPLEX(real_8), ALLOCATABLE             :: work(:,:)
    INTEGER                                  :: i, ia, iat, ierr, ig, ik, il, &
                                                ir, is, iv, k, l, ll, ly, &
                                                mmax, n2, n22, nwork
    INTEGER, DIMENSION(5)                    :: lpp = (/0,1,4,9,16/)
    LOGICAL                                  :: ferror, saved
    REAL(real_8) :: aiaux, cc, csri(2), disc, dln2, gmin, rmin, ryby2, &
      sum2(9), sumx, sumx2(9), sumy, sumy2(9), sumz, sumz2(9), vol, xmax
    CHARACTER(len=1), DIMENSION(5)           :: cang = (/'S','P','D','F','G'/)
    REAL(real_8), ALLOCATABLE                :: bsint(:), fint(:), gg(:), &
                                                temp(:,:)

    ryby2=2._real_8*ry
    dln2=1._real_8/LOG(2._real_8)
    iat=coresi%core_atom
    ia=iatpt(1,iat)
    is=iatpt(2,iat)
    l=coresi%core_lshell
    ll=2*l+1
    nwork=2*maxsys%mmaxx
    ! ..Allocation of local arrays
    ALLOCATE(work(maxsys%mmaxx,5),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fint(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bsint(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gg(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(temp(maxsys%mmaxx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(core_atwfr(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(core_cat(nsplpo,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(core_c0(nkpt%ngwk,ll),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..Construct core AO 
    CALL storps(atrg(1,is),core_atwfr,atwr%meshat(is),&
         coresi%core_nqsto,l,coresr%core_stoexp)
    ! ..Transform core AO in R-space to G-space
    CALL zeroing(core_cat)!,nsplpo*2)
    mmax=atwr%meshat(is)
    xmax=REAL(mmax,kind=real_8)
    n2=NINT(dln2*LOG(xmax)+0.499999_real_8)
    n22=2**n2
    rmin=LOG(atrg(1,is))
    gmin=LOG(SQRT(gvec_com%gcutw+gcutka)*parm%tpiba)-(mmax-1)*atwr%clogat(is)
    !$omp parallel do private(IL)
#ifdef __SR8000
    !poption parallel
#endif
    DO il=1,mmax
       gg(il)=(EXP(gmin+(il-1)*atwr%clogat(is))/parm%tpiba)**2
    ENDDO
    CALL zeroing(fint(1:n22))!,n22)
    !$omp parallel do private(IR)
#ifdef __SR8000
    !poption parallel
#endif
    DO ir=1,mmax
       fint(ir)=fpi*core_atwfr(ir)/atrg(ir,is)
    ENDDO
    ! ..Fourier transformation
    saved=.FALSE.
    CALL lsfbtr(fint,bsint,l,rmin,gmin,atwr%clogat(is),n2,saved,work(1,1),&
         work(1,2),work(1,3),work(1,5),nwork,disc)
    CALL tgrid(gg,mmax,ggng,nsplpo,bsint,maxsys%mmaxx,temp(1,1),temp(1,2),&
         temp(1,3))
    CALL dcopy(nsplpo,bsint(1),1,core_cat(1,1),1)
    IF (l.GT.0.AND.ggng(1).LT.1.e-12_real_8) core_cat(1,1)=0.0_real_8
    CALL curv1(nsplpo,ggng,core_cat(1,1),0.0_real_8,0.0_real_8,3,core_cat(1,2),&
         temp,0.0_real_8,ierr)
    ! ..Load core AO in plane wave basis
    vol=1._real_8/SQRT(parm%omega)
    ci=(0.0_real_8,1.0_real_8)**l
    CALL zeroing(core_c0(:,1:ll))!,nkpt%ngwk*ll)
    DO iv=1,ll
       ly=lpp(l+1)+iv
       DO ig=1,ncpw%ngw
          cc=curv2(hg(ig),nsplpo,ggng(1),core_cat(1,1),core_cat(1,2),&
               0.0_real_8)*vol
          core_c0(ig,iv)=ci*ylmr(ly,ig,gk(1,1))*cc*eigr(ig,iat,1)
       ENDDO
    ENDDO
    IF (paral%parent) THEN
       IF (tkpts%tkpnt) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,A)')'CORE_SPECTRA:',' K-POINTS CALCULATION'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(/,A,A)')'CORE_SPECTRA:',' GAMMA-POINT CALCULATION'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') ' EXCITED ATOM:'
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            '    NR   TYPE        X(BOHR)        Y(BOHR)        Z(BOHR)'
       IF (paral%io_parent)&
            WRITE(6,'(1X,I5,5X,A2,3F15.6)')&
            IAT,elem%EL(ions0%iatyp(IS)),(TAU0(K,IA,IS),K=1,3)
       IF (paral%io_parent)&
            WRITE(6,'(A)') '    EXCITED STATE:'
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            '                      n   l          Ec[Ha]         Ec[eV]'
       IF (paral%io_parent)&
            WRITE(6,'(22X,I1,3X,A,1X,2F15.6,/)')&
            coresi%core_nqsto,CANG(L+1),coresr%core_level,RYBY2*coresr%core_level
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            '  k    i   Ei[eV]  Ei-Ec[eV]   fi  ',&
            ' |<i|px|c>|**2  |<i|py|c>|**2  |<i|pz|c>|**2'
       ! mb-> write on file 60
       IF (paral%io_parent)&
            CALL fileopen(60,'XRAYSPEC.DAT',fo_def,ferror)
       IF (tkpts%tkpnt) THEN
          IF (paral%io_parent)&
               WRITE(60,'(A,A)')'CORE_SPECTRA:',' K-POINTS CALCULATION'
       ELSE
          IF (paral%io_parent)&
               WRITE(60,'(A,A)')'CORE_SPECTRA:',' GAMMA-POINT CALCULATION'
       ENDIF
       IF (paral%io_parent)&
            WRITE(60,'(/,A)') ' EXCITED ATOM:'
       IF (paral%io_parent)&
            WRITE(60,'(A,A)')&
            '    NR   TYPE        X(BOHR)        Y(BOHR)        Z(BOHR)'
       IF (paral%io_parent)&
            WRITE(60,'(1X,I5,5X,A2,3F15.6)')&
            IAT,elem%EL(ions0%iatyp(IS)),(TAU0(K,IA,IS),K=1,3)
       IF (paral%io_parent)&
            WRITE(60,'(A)') '    EXCITED STATE:'
       IF (paral%io_parent)&
            WRITE(60,'(A)')&
            '                      n   l          Ec[Ha]         Ec[eV]'
       IF (paral%io_parent)&
            WRITE(60,'(22X,I1,3X,A,1X,2F15.6,/)')&
            coresi%core_nqsto,CANG(L+1),coresr%core_level,RYBY2*coresr%core_level
       IF (paral%io_parent)&
            WRITE(60,'(A,A)')&
            '  k    i   Ei[eV]  Ei-Ec[eV]   fi  ',&
            ' |<i|px|c>|**2  |<i|py|c>|**2  |<i|pz|c>|**2  |<i|p|c>|**2'
    ENDIF
    ! ..Compute the transition matrix 
    IF (tkpts%tkpnt) THEN            ! K-points
       DO ik=1,nkpoint
          DO i=1,nstate
             DO iv=1,ll
                csumx=(0.0_real_8,0.0_real_8)
                csumy=(0.0_real_8,0.0_real_8)
                csumz=(0.0_real_8,0.0_real_8)
                ! mb-> first G > 0 with the k-points
#ifdef __SR8000
                !poption parallel
#endif
                !$omp parallel do private(IG,CAUX) reduction(+:CSUMX,CSUMY,CSUMZ)
                DO ig=1,ncpw%ngw
                   caux=CONJG(c0(ig,i,ik))*core_c0(ig,iv)
                   csumx=csumx+(rk(1,ik)+gk(1,ig))*caux
                   csumy=csumy+(rk(2,ik)+gk(2,ig))*caux
                   csumz=csumz+(rk(3,ik)+gk(3,ig))*caux
                ENDDO
                ! mb-> then G < 0 with the k-points
#ifdef __SR8000
                !poption parallel
#endif
                !$omp parallel do private(IG,CAUX) reduction(+:CSUMX,CSUMY,CSUMZ)
                DO ig=1,ncpw%ngw
                   caux=CONJG(c0(ig+ncpw%ngw,i,ik))*core_c0(ig,iv)
                   csumx=csumx+(rk(1,ik)-gk(1,ig))*caux
                   csumy=csumy+(rk(2,ik)-gk(2,ig))*caux
                   csumz=csumz+(rk(3,ik)-gk(3,ig))*caux
                ENDDO
                csri(1)=REAL(csumx)
                csri(2)=AIMAG(csumx)
                CALL mp_sum(csri,2,parai%allgrp)
                csumx=CMPLX(csri(1),csri(2),kind=real_8)
                csri(1)=REAL(csumy)
                csri(2)=AIMAG(csumy)
                CALL mp_sum(csri,2,parai%allgrp)
                csumy=CMPLX(csri(1),csri(2),kind=real_8)
                csri(1)=REAL(csumz)
                csri(2)=AIMAG(csumz)
                CALL mp_sum(csri,2,parai%allgrp)
                csumz=CMPLX(csri(1),csri(2),kind=real_8)
                sumx2(iv)=REAL(CONJG(csumx)*csumx)
                sumy2(iv)=REAL(CONJG(csumy)*csumy)
                sumz2(iv)=REAL(CONJG(csumz)*csumz)
                sum2(iv)=sumx2(iv)+sumy2(iv)+sumz2(iv)
             ENDDO
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(I3,I5,F9.3,F11.3,F6.2,27E15.7)')&
                     IK,I,RYBY2*WE(I,IK),RYBY2*(WE(I,IK)-coresr%core_level),F(I,IK),&
                     (SUMX2(IV),SUMY2(IV),SUMZ2(IV),IV=1,LL)
                IF (paral%io_parent)&
                     WRITE(60,'(I3,I5,F9.3,F11.3,F6.2,36E15.7)')&
                     IK,I,RYBY2*WE(I,IK),RYBY2*(WE(I,IK)-coresr%core_level),F(I,IK),&
                     (SUMX2(IV),SUMY2(IV),SUMZ2(IV),SUM2(IV),IV=1,LL)
             ENDIF
          ENDDO
       ENDDO
    ELSE                      ! Gamma-point only
       ik=1
       DO i=1,nstate
          DO iv=1,ll
             sumx=0.0_real_8
             sumy=0.0_real_8
             sumz=0.0_real_8
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IG,AIAUX) reduction(+:SUMX,SUMY,SUMZ)
             DO ig=1,ncpw%ngw
                aiaux=AIMAG(CONJG(c0(ig,i,ik))*core_c0(ig,iv))
                sumx=sumx+gk(1,ig)*aiaux
                sumy=sumy+gk(2,ig)*aiaux
                sumz=sumz+gk(3,ig)*aiaux
             ENDDO
             CALL mp_sum(sumx,parai%allgrp)
             CALL mp_sum(sumy,parai%allgrp)
             CALL mp_sum(sumz,parai%allgrp)
             sumx2(iv)=4.0_real_8*sumx*sumx
             sumy2(iv)=4.0_real_8*sumy*sumy
             sumz2(iv)=4.0_real_8*sumz*sumz
             sum2(iv)=sumx2(iv)+sumy2(iv)+sumz2(iv)
          ENDDO
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(I3,I5,F9.3,F11.3,F6.2,27E15.7)')&
                  IK,I,RYBY2*WE(I,IK),RYBY2*(WE(I,IK)-coresr%core_level),F(I,IK),&
                  (SUMX2(IV),SUMY2(IV),SUMZ2(IV),IV=1,LL) 
             IF (paral%io_parent)&
                  WRITE(60,'(I3,I5,F9.3,F11.3,F6.2,36E15.7)')&
                  IK,I,RYBY2*WE(I,IK),RYBY2*(WE(I,IK)-coresr%core_level),F(I,IK),&
                  (SUMX2(IV),SUMY2(IV),SUMZ2(IV),SUM2(IV),IV=1,LL) 
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Deallocation of local arrays
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fint,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bsint,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(core_atwfr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(core_cat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(core_c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            CALL fileclose(60)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE core_spectra
  ! ==================================================================

END MODULE core_spect_utils
