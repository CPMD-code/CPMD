MODULE nlforce_utils
  USE cppt,                            ONLY: twnl
  USE cvan,                            ONLY: deeq,&
                                             dvan,&
                                             qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: eigr,&
                                             fnl
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE tbxc,                            ONLY: toldcode
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlforce
  PUBLIC :: give_scr_nlforce

CONTAINS

  ! ==================================================================
  SUBROUTINE nlforce(c2,f,gam,auxc,ddia,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(ncpw%ngw,*)
    REAL(real_8)                             :: f(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ddia(maxsys%nax,nstate)
    COMPLEX(real_8)                          :: auxc(ncpw%ngw,nstate)
    REAL(real_8)                             :: gam(nstate,*)

    COMPLEX(real_8)                          :: ct, ctm
    INTEGER                                  :: i, ia, ig, is, isa, isa0, &
                                                ispin, isub, iv, j, jv, ki, &
                                                kj, l, l2, li, lj, nmax, nmin
    REAL(real_8)                             :: dd, fac, ffi

! Variables
! ==--------------------------------------------------------------==
! == Compute the non-local contributions to the electron gradient ==
! ==--------------------------------------------------------------==

    CALL tiset('   NLFORCE',isub)
    IF (imagp.EQ.2) CALL stopgm('NLFORCE','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (cntl%tfdist) CALL stopgm('NLFORCE','FNL DIST. NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! FIXME: AK 2008/02/14: this was broken by changes from 2006/05/12.
    ! not sure whether those were the real cause.
    ! undoing those changes does not make it work
    ! again. so we stop here until somebody fixes it
    ! UPDATE AK 2008/05/24: it looks as if it can be worked around using 
    ! OLDCODE. it seems to be faster, too.
    IF (pslo_com%tivan.AND.cntl%tlsd) THEN
       IF (.NOT.toldcode)&
            CALL stopgm('NLFORCE','VANDERBILT WITH LSD REQUIRES USE OF'&
            // ' "OLDCODE" FLAG IN &DFT SECTION.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          ! ==--------------------------------------------------------==
          ! ==  VANDERBILT PP                                         ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(ddia)!,maxsys%nax*nstate)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ispin=1
                IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
                DO jv=1,nlps_com%ngh(is)
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      ddia(ia,i)=ddia(ia,i)-ffi*fnl(1,isa,jv,i,1)*&
                           (deeq(isa,jv,iv,ispin)+dvan(jv,iv,is))
                   ENDDO
                   IF (ABS(qq(jv,iv,is)).GT.1.e-5_real_8) THEN
                      IF (cnti%iproj.EQ.1.AND.paral%parent) THEN
                         fac=-gam(i,i)*qq(jv,iv,is)
                         CALL daxpy(ions0%na(is),fac,fnl(1,isa0+1,jv,i,1),1,&
                              ddia(1,i),1)
                      ELSE
                         IF (cntl%tlsd) THEN
                            IF (i.LE.spin_mod%nsup) THEN
                               nmin=1
                               nmax=spin_mod%nsup
                            ELSE
                               nmin=spin_mod%nsup+1
                               nmax=nstate
                            ENDIF
                         ELSE
                            nmin=1
                            nmax=nstate
                         ENDIF
                         DO j=nmin,nmax
                            fac=-gam(i,j)*qq(jv,iv,is)
                            CALL daxpy(ions0%na(is),fac,fnl(1,isa0+1,jv,j,1),1,&
                                 ddia(1,i),1)
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
             CALL mp_sum(ddia,nstate*maxsys%nax,parai%allgrp)
             IF (ncpw%ngw.GT.0)&
                  CALL dgemm('N','N',2*ncpw%ngw,nstate,ions0%na(is),1.0_real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,ddia,maxsys%nax,0.0_real_8,auxc,2*ncpw%ngw)
             ctm=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             !$omp parallel do private(IG,I)
             DO i=1,nstate
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! ==--------------------------------------------------------==
          ! ==  Stefan Goedecker PP                                   ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(auxc)!,ngw*nstate)
             CALL zeroing(ddia)!,maxsys%nax*nstate)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   dd=0.0_real_8
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         ddia(ia,i)=ddia(ia,i)+&
                              fnl(1,isa,jv,i,1)*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             CALL mp_sum(ddia,nstate*maxsys%nax,parai%allgrp)
             IF (ncpw%ngw.GT.0)&
                  CALL dgemm("N","N",2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,ddia(1,1),maxsys%nax,0._real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             !$omp parallel do private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          ! ==--------------------------------------------------------==
          ! ==  BACHELET HAMANN SCHLUTER                              ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(ddia)!,maxsys%nax*nstate)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                CALL dcopy(ions0%na(is),fnl(1,isa0+1,iv,i,1),1,ddia(1,i),1)
             ENDDO
             CALL mp_sum(ddia,nstate*maxsys%nax,parai%allgrp)
             IF (ncpw%ngw.GT.0)&
                  CALL dgemm('N','N',2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,ddia,maxsys%nax,0.0_real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)*wsg(is,iv)
             !$omp parallel do private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    CALL tihalt('   NLFORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlforce
  ! ==================================================================
  SUBROUTINE give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: il_gam, il_auxc, il_ddia, &
                                                nstate

! ==--------------------------------------------------------------==

    il_gam =imagp*nstate*nstate
    il_auxc=2*nkpt%ngwk*nstate+2
    il_ddia=2*maxsys%nax*nstate+2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_nlforce
  ! ==================================================================


END MODULE nlforce_utils
