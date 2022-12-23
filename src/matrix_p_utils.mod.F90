MODULE matrix_p_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE response_pmod,                   ONLY: &
       ddfnl_ddk, ddtwnl_ddk, dfnl_dk, dtwnl_dk, fnl00, h0_11, h1_00, h1_10, &
       h2_00, vofrho0
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE sfac,                            ONLY: eigr,&
                                             fnl
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean_k
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: matrix_p
  PUBLIC :: give_scr_matrix_p

CONTAINS

  ! ==================================================================
  SUBROUTINE matrix_p(c0,cu1,nstate,psi)
    ! ==================================================================
    ! ==                        COMPUTES                              ==
    ! ==        The  matrixes NSTATExNSTATE not dependent on k vector ==
    ! ==    The final Hamiltonian matrix is constitued by combinations==
    ! ==    of these matrixes. The computed matrixes depend on the    ==
    ! ==    x,y and z directions and on the overlapping wavefunctions ==
    ! ==    At last they are:
    ! ==    <c0|H1x|c0>,<c0|H1y|c0>,<c0|H1z|c0>
    ! ==    <c0|H2xx|c0>,<c0|H2xy|c0>,<c0|H2xz|c0>
    ! ==    <c0|H2yx|c0>,<c0|H2yy|c0>,<c0|H2yz|c0>
    ! ==    <c0|H2zx|c0>,<c0|H2zy|c0>,<c0|H2zz|c0>
    ! ==    <c1x|H1x|c0>,<c1x|H1y|c0>,<c1x|H1z|c0>
    ! ==    <c1y|H1x|c0>,<c1y|H1y|c0>,<c1y|H1z|c0>
    ! ==    <c1z|H1x|c0>,<c1z|H1y|c0>,<c1z|H1z|c0>
    ! ==    <c1x|H0|c1x>,<c1x|H0|c1y>,<c1x|H0|c1z>
    ! ==    <c1y|H0|c1x>,<c1y|H0|c1y>,<c1y|H0|c1z>
    ! ==    <c1z|H0|c1x>,<c1z|H0|c1y>,<c1z|H0|c1z>
    ! ==--------------------------------------------------------------==

    ! F(istate) = 1
    ! Matrix defined with reverse sign
    ! Matrix complex and hermitian
    ! Computed blocks: <c0|H|c0> <c1|H|c0> <c1|H|c1>

    ! forces always put in vector c2
    ! then called overlap routine to build  the matrix


    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cu1(ncpw%ngw,nstate,3), &
                                                c0(ncpw%ngw,nstate), &
                                                psi(maxfftn)

    CHARACTER(*), PARAMETER                  :: procedureN = 'matrix_p'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ctm, g, i_to_l
    COMPLEX(real_8), ALLOCATABLE             :: c00(:,:), c2(:,:)
    INTEGER                                  :: ia, iat, ierr, ig, il_auxc, &
                                                il_ddia, isp_, istate, isub, &
                                                iv, jv, k2_, k_, ki, kj, kk_, &
                                                l, l2, li, lj
    REAL(real_8)                             :: dd0, dd1(3), dd2(3), &
                                                foc(nstate), tcpu, time1, &
                                                time2
    REAL(real_8), ALLOCATABLE                :: auxc(:), ddia(:)

! ==--------------------------------------------------------------==

    CALL tiset('    MATRIX_P',isub)
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    ! ALLOCATE MEMORY

    ALLOCATE(c2(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c2)!,2*ngw*nstate)
    ALLOCATE(c00(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c00)!,2*ngw*nstate)
    ALLOCATE(h1_00(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(h1_00)!,3*nstate*nstate)
    ALLOCATE(h2_00(nstate,nstate,3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(h2_00)!,3*3*nstate*nstate)
    ALLOCATE(h1_10(nstate,nstate,3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(h1_10)!,3*3*nstate*nstate)
    ALLOCATE(h0_11(nstate,nstate,3,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(h0_11)!,3*3*nstate*nstate)


    DO istate = 1,nstate
       DO ig = 1,ncpw%ngw
          c00(ig,istate) = c0(ig,istate)
          c00(ig+ncpw%ngw,istate) = CONJG(c0(ig,istate))
       ENDDO
    ENDDO

    IF (geq0)   CALL zclean_k(c00,nstate,ncpw%ngw)

    ! M_0H0 = <c0|H|c0>
    ! devo calcolare:
    ! H1_x_00 = <c0|-i dx|c0>
    ! H1_y_00 = <c0|-i dy|c0>
    ! H1_z_00 = <c0|-i dz|c0>

    DO k_ = 1,3
       CALL zeroing(c2)!,2*ngw*nstate)
       DO istate = 1,nstate
          DO ig = 1,ncpw%ngw
             g=CMPLX(parm%tpiba*gk(k_,ig),0.0_real_8,kind=real_8)
             c2(ig,istate) =  g*c00(ig,istate)
             c2(ig+ncpw%ngw,istate)= -g*c00(ig+ncpw%ngw,istate)
          ENDDO
       ENDDO

       ! H1_x_00 =+ <c0|d_x P>w<P |c0>+<c0|P>w<d_x P|c0>
       ! H1_y_00 =+ <c0|d_y P>w<P |c0>+<c0|P>w<d_y P|c0>
       ! H1_z_00 =+ <c0|d_z P>w<P |c0>+<c0|P>w<d_z P|c0>

       iat = 0
       DO istate = 1,nstate
          DO isp_ = 1,ions1%nsp
             ! ==--------------------------------------------------------------==
             ! ==  VANDERBILT                                                  ==
             ! ==--------------------------------------------------------------==
             IF (pslo_com%tvan(isp_)) THEN
                CALL stopgm('K_PERT_PREP',&
                     'Vanderbilt PP not implemented',& 
                     __LINE__,__FILE__)

             ELSEIF (sgpp1%tsgp(isp_)) THEN
                ! ==--------------------------------------------------------------==
                ! ==  STEFAN GOEDECKER                                            ==
                ! ==--------------------------------------------------------------==
                DO ia=1,ions0%na(isp_)
                   iat=iat+1
                   DO iv=1,nlps_com%ngh(isp_)
                      l=nghtol(iv,isp_)+1
                      ki=sgpp2%lfval(iv,isp_)
                      li=sgpp2%lpval(iv,isp_)
                      dd0=0.0_real_8
                      DO k2_ = 1,3
                         dd2(k2_) = 0.0_real_8
                         dd1(k2_) = 0.0_real_8
                      ENDDO
                      DO jv=1,nlps_com%ngh(isp_)
                         l2=nghtol(jv,isp_)+1
                         lj=sgpp2%lpval(jv,isp_)
                         IF (l.EQ.l2.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,isp_)
                            dd0=dd0 +&
                                 fnl00(1,iat,jv,istate,1)*sgpp2%hlsg(ki,kj,l,isp_)
                            dd1(k_)=dd1(k_)+dfnl_dk(2,iat,jv,istate, k_)*&
                                 sgpp2%hlsg(ki,kj,l,isp_)
                         ENDIF
                      ENDDO
                      ctm=- CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,isp_)
                      DO ig=1,ncpw%ngw
                         c2(ig,istate)=c2(ig,istate)+ctm * eigr(ig, iat,1)*&
                              (twnl(ig,iv,isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*dd1(k_)+&
                              dtwnl_dk(ig,iv,isp_,k_)*dd0)

                         c2(ig+ncpw%ngw,istate)=c2(ig+ncpw%ngw,istate)+ctm *&
                              CONJG(eigr(ig,iat,1))*((-1)**(nghtol(iv,&
                              isp_))* twnl(ig,iv, isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*&
                              dd1(k_)+dtwnl_dk(ig+ncpw%ngw,iv,isp_, k_)*dd0)
                      ENDDO
                   ENDDO
                ENDDO
             ELSE

                DO ia = 1,ions0%na(isp_)
                   iat = iat + 1
                   DO iv = 1,nlps_com%ngh(isp_)

                      ! ==--------------------------------------------------------------==
                      ! ==  BACHELET HAMANN SCHLUTER                                    ==
                      ! ==--------------------------------------------------------------==
                      i_to_l = CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,isp_)
                      ctm    = - i_to_l *wsg(isp_,iv)
                      DO ig=1,ncpw%ngw
                         c2(ig,istate)=c2(ig,istate)+ctm * eigr(ig, iat,1)*&
                              (twnl(ig,iv,isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*&
                              dfnl_dk(2,iat,iv, istate,k_)+dtwnl_dk(ig,iv,&
                              isp_,k_)*fnl00(1,iat,iv,istate,1))

                         c2(ig+ncpw%ngw,istate)=c2(ig+ncpw%ngw,istate)+ctm *&
                              CONJG(eigr(ig,iat,1))*((-1)**(nghtol(iv,&
                              isp_))* twnl(ig,iv, isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*&
                              dfnl_dk(2,iat,iv,istate,k_)+ dtwnl_dk(ig+ncpw%ngw,&
                              iv,isp_,k_)*fnl00(1,iat,iv,istate,1))
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          iat = 0
       ENDDO
       IF (geq0)   CALL zclean_k(c2,nstate,ncpw%ngw)

       CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,zone,c00(1,1),2*ncpw%ngw,&
            c2(1,1),2*ncpw%ngw,zzero,h1_00(1,1,k_),nstate)

       CALL mp_sum(h1_00(:,:,k_),nstate*nstate,parai%allgrp)
    ENDDO

    ! 
    ! H2_ab_00 = <c0|d_a P>w<d_b P |c0>+
    ! <c0|d_ab P>w<P |c0>+<c0| P>w<d_ab P|c0>
    ! 
    DO k_ = 1,3
       DO kk_ = 1,3
          CALL zeroing(c2)!,2*ngw*nstate)
          iat = 0
          ! ==--------------------------------------------------------------==
          DO istate = 1,nstate
             iat = 0
             DO isp_ = 1,ions1%nsp
                IF (sgpp1%tsgp(isp_)) THEN
                   ! ==--------------------------------------------------------------==
                   ! ==  STEFAN GOEDECKER                                            ==
                   ! ==--------------------------------------------------------------==
                   DO ia=1,ions0%na(isp_)
                      iat=iat+1
                      DO iv=1,nlps_com%ngh(isp_)
                         l=nghtol(iv,isp_)+1
                         ki=sgpp2%lfval(iv,isp_)
                         li=sgpp2%lpval(iv,isp_)
                         dd0=0.0_real_8
                         DO k2_ = 1,3
                            dd2(k2_) = 0.0_real_8
                            dd1(k2_) = 0.0_real_8
                         ENDDO
                         DO jv=1,nlps_com%ngh(isp_)
                            l2=nghtol(jv,isp_)+1
                            lj=sgpp2%lpval(jv,isp_)
                            IF (l.EQ.l2.AND.li.EQ.lj) THEN
                               kj=sgpp2%lfval(jv,isp_)
                               dd0=dd0 +&
                                    fnl00(1,iat,jv,istate,1)*sgpp2%hlsg(ki,kj,l,isp_)
                               dd1(kk_)=dd1(kk_)+dfnl_dk(2,iat,jv,&
                                    istate,kk_)*sgpp2%hlsg(ki,kj,l,isp_)
                               dd2(kk_)=dd2(kk_)+ddfnl_ddk(1,iat,jv,&
                                    istate,k_,kk_)*sgpp2%hlsg(ki,kj,l,isp_)
                            ENDIF
                         ENDDO
                         ctm=- CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,isp_)
                         DO ig=1,ncpw%ngw
                            c2(ig,istate)= &
                                 c2(ig,istate) + &
                                 ctm*eigr(ig, iat,1)*(dtwnl_dk(ig,iv,isp_,k_)&
                                 * CMPLX(0._real_8,1._real_8,kind=real_8)*dd1(kk_) &
                                 + 0.5_real_8*ddtwnl_ddk(ig,iv,isp_,k_,kk_) &
                                 * dd0+twnl(ig,iv,isp_,1)*0.5_real_8*dd2(kk_))

                            c2(ig+ncpw%ngw,istate)= &
                                 c2(ig+ncpw%ngw,istate) + &
                                 ctm * CONJG(eigr(ig,iat,1))*&
                                 (dtwnl_dk(ig+ncpw%ngw,iv, isp_,k_)*&
                                 CMPLX(0._real_8,1._real_8,kind=real_8)*dd1(kk_) + &
                                 0.5_real_8* ddtwnl_ddk(ig+ncpw%ngw,iv,isp_,k_,kk_)* &
                                 dd0+0.5_real_8* (-1)**(nghtol(iv,isp_)) * &
                                 twnl(ig,iv,isp_,1)* dd2(kk_))
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE

                   DO ia = 1,ions0%na(isp_)
                      iat = iat + 1

                      DO iv = 1,nlps_com%ngh(isp_)
                         i_to_l = CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,&
                              isp_)
                         ctm    = - i_to_l *wsg(isp_,iv)

                         DO ig=1,ncpw%ngw
                            c2(ig,istate)=c2(ig,istate)+&
                                 ctm * eigr(ig, iat,1) *&
                                 (dtwnl_dk(ig,iv,isp_,k_)*CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dfnl_dk(2, iat,iv,istate,kk_) +&
                                 0.5_real_8*ddtwnl_ddk(ig,iv,isp_,k_,kk_)*&
                                 fnl00(1,iat,iv,istate,1) +&
                                 0.5_real_8*twnl(ig,iv,isp_,1)*&
                                 ddfnl_ddk(1,iat,iv,istate,k_,kk_))

                            c2(ig+ncpw%ngw,istate)= &
                                 c2(ig+ncpw%ngw,istate)+&
                                 ctm*CONJG(eigr(ig,iat,1))*&
                                 (dtwnl_dk(ig+ncpw%ngw,iv,isp_,k_)*&
                                 CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dfnl_dk(2,iat,iv,istate,kk_) + &
                                 0.5_real_8*ddtwnl_ddk(ig+ncpw%ngw,iv,isp_,k_,kk_)*&
                                 fnl00(1,iat,iv,istate,1)+&
                                 0.5_real_8*(-1)**(nghtol(iv,isp_))*&
                                 twnl(ig,iv,isp_,1)*&
                                 ddfnl_ddk(1,iat,iv,istate,k_,kk_))
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          IF (geq0)   CALL zclean_k(c2,nstate,ncpw%ngw)

          CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,zone,c00(1,1),2*ncpw%ngw,&
               c2(1,1),2*ncpw%ngw,zzero,h2_00(1,1,k_,kk_),nstate)

          CALL sumhmat(h2_00(1,1,k_,kk_),nstate)
       ENDDO
    ENDDO

    ! <c1|H|c0>

    DO k_ = 1,3
       CALL zeroing(c00)!,2*ngw*nstate)
       DO istate = 1,nstate
          DO ig = 1,ncpw%ngw
             c00(ig,istate) = cu1(ig,istate,k_)
             c00(ig+ncpw%ngw,istate) = CONJG(cu1(ig,istate,k_))
          ENDDO
       ENDDO

       IF (geq0)   CALL zclean_k(c00,nstate,ncpw%ngw)

       ! H1_x_10 = <c1|-i dx|c0>
       ! H1_y_10 = <c1|-i dy|c0>
       ! H1_z_10 = <c1|-i dz|c0>

       DO kk_ = 1,3
          CALL zeroing(c2)!,2*ngw*nstate)

          DO istate = 1,nstate

             DO ig = 1,ncpw%ngw
                g=CMPLX(parm%tpiba*gk(kk_,ig),0.0_real_8,kind=real_8)

                c2(ig,istate) = g*c0(ig,istate)
                c2(ig+ncpw%ngw,istate)=-g*CONJG(c0(ig,istate))
             ENDDO
          ENDDO

          ! H1_x_10 =+ <c0|d_x P>w<P |c0>+<c0|P>w<d_x P|c0>
          ! H1_y_10 =+ <c0|d_y P>w<P |c0>+<c0|P>w<d_y P|c0>
          ! H1_z_10 =+ <c0|d_z P>w<P |c0>+<c0|P>w<d_z P|c0>

          iat = 0
          DO istate = 1,nstate
             DO isp_ = 1,ions1%nsp
                IF (sgpp1%tsgp(isp_)) THEN
                   ! ==--------------------------------------------------------------==
                   ! ==  STEFAN GOEDECKER                                            ==
                   ! ==--------------------------------------------------------------==
                   DO ia=1,ions0%na(isp_)
                      iat=iat+1
                      DO iv=1,nlps_com%ngh(isp_)
                         l=nghtol(iv,isp_)+1
                         ki=sgpp2%lfval(iv,isp_)
                         li=sgpp2%lpval(iv,isp_)
                         dd0=0.0_real_8
                         DO k2_ = 1,3
                            dd2(k2_) = 0.0_real_8
                            dd1(k2_) = 0.0_real_8
                         ENDDO
                         DO jv=1,nlps_com%ngh(isp_)
                            l2=nghtol(jv,isp_)+1
                            lj=sgpp2%lpval(jv,isp_)
                            IF (l.EQ.l2.AND.li.EQ.lj) THEN
                               kj=sgpp2%lfval(jv,isp_)
                               dd0=dd0 +&
                                    fnl00(1,iat,jv,istate,1)*sgpp2%hlsg(ki,kj,l,isp_)
                               dd1(kk_)=dd1(kk_)+dfnl_dk(2,iat,jv,&
                                    istate,kk_)*sgpp2%hlsg(ki,kj,l,isp_)
                            ENDIF
                         ENDDO
                         ctm=- CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,isp_)
                         DO ig=1,ncpw%ngw
                            c2(ig,istate)=c2(ig,istate)+ctm * eigr(ig, iat,1)*&
                                 (twnl(ig,iv,isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dd1(kk_)+ dtwnl_dk(ig,iv,isp_,kk_)*dd0)

                            c2(ig+ncpw%ngw,istate)=c2(ig+ncpw%ngw,istate)+ctm *&
                                 CONJG(eigr(ig,iat,1))*((-1)**(nghtol(iv,&
                                 isp_))*twnl(ig,iv, isp_,1)*CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dd1(kk_)+dtwnl_dk(ig+ncpw%ngw,iv,isp_, kk_)*dd0)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO ia = 1,ions0%na(isp_)
                      iat = iat + 1
                      DO iv = 1,nlps_com%ngh(isp_)
                         i_to_l = CMPLX(0.0_real_8,-1.0_real_8,kind=real_8)**nghtol(iv,&
                              isp_)
                         ctm    = - i_to_l *wsg(isp_,iv)

                         DO ig=1,ncpw%ngw
                            c2(ig,istate)=c2(ig,istate)+ctm * eigr(ig, iat,1)*&
                                 (twnl(ig,iv,isp_,1) *CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dfnl_dk(2,iat,iv, istate,kk_)+dtwnl_dk(ig,iv,&
                                 isp_,kk_)*fnl00(1,iat,iv,istate,1))

                            c2(ig+ncpw%ngw,istate)=c2(ig+ncpw%ngw,istate)+ctm *&
                                 CONJG(eigr(ig,iat,1))*((-1)**(nghtol(iv,&
                                 isp_))*twnl(ig,iv, isp_,1)*CMPLX(0._real_8,1._real_8,kind=real_8)*&
                                 dfnl_dk(2,iat,iv,istate,kk_)+ dtwnl_dk(ig+&
                                 ncpw%ngw,iv,isp_,kk_)*fnl00(1,iat,iv,istate,1))
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             iat = 0
          ENDDO
          IF (geq0)   CALL zclean_k(c2,nstate,ncpw%ngw)

          CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,zone,c00(1,1),2*ncpw%ngw,&
               c2(1,1),2*ncpw%ngw,zzero,h1_10(1,1,k_,kk_),nstate)

          CALL mp_sum(h1_10(:,:,k_,kk_),nstate*nstate,parai%allgrp)
       ENDDO
    ENDDO

    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! <c1|H|c1>
    ! in this case  only H0 is applied
    ! the same structure of h0psi1 is used, but 
    ! H0_11 defined as real, since CU1 are real

    DO istate = 1,nstate
       foc(istate)=1.0_real_8
    ENDDO

    ALLOCATE(c2(ncpw%ngw,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c2)!,ngw*nstate)

    CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
    ALLOCATE(auxc(il_auxc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(auxc)!,il_auxc)
    ALLOCATE(ddia(il_ddia),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ddia)!,il_ddia)

    DO k_ = 1,3
       CALL zeroing(c2)!,ngw*nstate)
       CALL zeroing(fnl)!,ions1%nat*maxsys%nhxs*nstate*1)

       CALL rnlsm(cu1(:,:,k_),nstate,1,1,.FALSE.)

       ! c2(g)  =  <g|  v_Hxc(0)  |cu1>; vpsi also initializes c12
       CALL vpsi(cu1(:,:,k_),c2,foc,VofRho0,psi,nstate,1,clsd%nlsd,.TRUE.)

       ! c2(g) +=  <g|  | fnl[cu1] >
       CALL fnonloc(c2,foc,nstate,1,clsd%nlsd,.TRUE.)

       DO kk_ = 1,3
          ! calculate the hamiltonian matrix <c1|H0|c1> = <c1|c2>:
          CALL ovlap(nstate,h0_11(:,:,kk_,k_),cu1(:,:,kk_),c2)

          CALL mp_sum(h0_11(:,:,kk_,k_),nstate*nstate,parai%allgrp)
       ENDDO
    ENDDO

    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(a,t50,f8.2,a8)')&
            ' cpu time for matrix elements calculation:',&
            tcpu,' seconds'
    ENDIF

    CALL tihalt('    MATRIX_P',isub)
    RETURN
  END SUBROUTINE matrix_p
  ! ==============================================================
  SUBROUTINE give_scr_matrix_p(lmatrix_p,tag,nstate)
    INTEGER                                  :: lmatrix_p
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_auxc, l_ddia, l_rnlsm, lmat

    CALL give_scr_rnlsm(l_rnlsm,tag,nstate,.FALSE.)
    CALL give_scr_fnonloc(l_auxc,l_ddia,nstate)
    lmat = 2*nstate*(nstate+1)

    lmatrix_p = MAX(l_rnlsm,l_auxc+l_ddia,lmat)
    RETURN
  END SUBROUTINE give_scr_matrix_p
  ! ==================================================================

END MODULE matrix_p_utils
