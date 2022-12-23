MODULE fukui_p_utils
  USE cnst,                            ONLY: ry
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: indz,&
                                             nzh
  USE d_mat_p_utils,                   ONLY: d_mat_nonloc
  USE densrd_utils,                    ONLY: densrd
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00,&
                                             nf,&
                                             numf,&
                                             tweight,&
                                             wghtf
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rnlfor_utils,                    ONLY: rnlfor
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE sfac,                            ONLY: dfnl,&
                                             fnl
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho_p_utils,                 ONLY: v1ofrho_p
  USE vpsi_p_utils,                    ONLY: v1psi0_p
  USE wrgeo_utils,                     ONLY: wrgeof
!!use densto_utils, only : densto
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fukui_p

CONTAINS

  ! ==================================================================
  SUBROUTINE fukui_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,&
       nstate)
    ! ==--------------------------------------------------------------==
    ! nuclear fukui functions 


    COMPLEX(real_8)                          :: psi(maxfftn)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fukui_p'

    CHARACTER(len=15)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: filen
    COMPLEX(real_8)                          :: v1_nonloc(1)
    COMPLEX(real_8), ALLOCATABLE             :: cf(:,:), eirop1(:), psif(:), &
                                                v1_loc(:), vf(:), vtemp(:)
    INTEGER :: i, i1, i2, ia, iat, icol, ierr, ig, il_cf, il_psif, il_rhof, &
      il_v1_loc, ir, is, isa, isa0, ist, isub, iv, k, n1, n2, nbeg, nend, nft
    REAL(real_8)                             :: hardness, r1, scal, sumr
    REAL(real_8), ALLOCATABLE                :: ff(:), ffmat(:,:,:), &
                                                fion_p(:), hard(:,:), &
                                                rhof(:), xcf(:)

! variables
! ==--------------------------------------------------------------==

    CALL tiset('   fukui_p',isub)

    il_cf=2*ncpw%ngw*(numf+1)
    ALLOCATE(cf(ncpw%ngw,numf+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL rhoe_psi_size(il_rhof,il_psif)
    ALLOCATE(rhof(il_rhof),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xcf(il_rhof),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psif(il_psif),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vf(il_psif),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eirop1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vtemp)!,nhg)
    CALL zeroing(eirop1)!,nhg)
    il_v1_loc=2*ncpw%nhg
    ALLOCATE(v1_loc(il_v1_loc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hard(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion_p(3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ffmat(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ff(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ff)!,nstate)

    ! ==--------------------------------------------------------------==
    ! reads the wavefunction whose occupation number changes

    IF (cntl%tlsd) THEN   ! LSD
       DO i=1,numf
          ist=-nf(i)
          IF (ist.GT.spin_mod%nsup) THEN
             ist=ist-spin_mod%nsup
             cflbod='WAVEFUNCTION.B'
          ELSE
             cflbod='WAVEFUNCTION.A'
          ENDIF
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ist
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
          CALL densrd(vtemp,filen)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%ngw
             cf(ig,i)=vtemp(ig)
          ENDDO
       ENDDO
    ELSE
       DO i=1,numf! LDA
          ist=-nf(i)
          cflbod='WAVEFUNCTION.'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ist
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
          CALL densrd(vtemp,filen)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%ngw
             cf(ig,i)=vtemp(ig)
          ENDDO
       ENDDO
    ENDIF

    ! eventually compute linear combination

    IF (tweight) THEN
       nbeg=1
       nend=numf
       nft=numf
    ELSE
       scal=0.0_real_8
       DO i=1,numf
          scal=scal+wghtf(i)*wghtf(i)
       ENDDO
       scal=1._real_8/SQRT(scal)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'All coeffs are scaled by',scal
       ENDIF
       CALL dgemm('N','N',2*ncpw%ngw,1,numf,scal,cf(1,1),2*ncpw%ngw,wghtf,&
            numf+1,0._real_8,cf(1,numf+1),2*ncpw%ngw)
       nbeg=numf+1
       nend=numf+1
       nft=1
       wghtf(numf+1)=1.0_real_8
    ENDIF

    ! adds non-local contribution to the nuclear fukui functions
    CALL zeroing(c1)!,ngw*nstate)
    CALL dcopy(2*ncpw%ngw*nft,cf(1,nbeg),1,c1,1)
    CALL rnlsm(c1,nstate,1,1,.TRUE.)
    CALL zeroing(ffmat)!,3*maxsys%nax*maxsys%nsx)
    DO i=1,nft
       ff(i)=wghtf(nbeg+i-1)
    ENDDO
    !$omp parallel do private(i)
    DO i=nft+1,nstate
       ff(i)=0.0_real_8
    ENDDO
    CALL rnlfor(ffmat,ff,wk,nstate,1)
    CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,ffmat,1)

    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dfnl00(imagp,ions1%nat,maxsys%nhxs,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL rnlsm(c0,nstate,1,1,.TRUE.)

    ! save  fnl and  dfnl in fnl00 and dfnl00      
    isa0=0
    DO is=1,ions1%nsp           ! parallel impossible
       DO iv=1,nlps_com%ngh(is)
          DO ia=1,ions0%na(is)! parallel possible
             isa=isa0+ia
             DO i=1,nstate
                fnl00(1,isa,iv,i,1)=fnl(1,isa,iv,i,1)
             ENDDO
             DO i=1,parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
                dfnl00(1,isa,iv,1,i,1)=dfnl(1,isa,iv,1,i,1)
                dfnl00(1,isa,iv,2,i,1)=dfnl(1,isa,iv,2,i,1)
                dfnl00(1,isa,iv,3,i,1)=dfnl(1,isa,iv,3,i,1)
             ENDDO
          ENDDO
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO

    ! ==--------------------------------------------------------------==
    ! calculates the perturbating density

    CALL zeroing(rhof)!,nnr1)
    DO i=nbeg,nend
       CALL zeroing(psif)!,maxfftn)
       DO ig=1,ncpw%ngw
          psif(nzh(ig))=cf(ig,i)
          psif(indz(ig))=CONJG(cf(ig,i))
       ENDDO
       IF (geq0) psif(nzh(1))=cf(1,i)
       CALL  invfftn(psif,.FALSE.,parai%allgrp)
       !$omp parallel do private(IR,R1)
       DO ir=1,fpar%nnr1
          r1=REAL(psif(ir))
          rhof(ir)=rhof(ir)+r1*r1/parm%omega*wghtf(i)
       ENDDO
    ENDDO

    sumr=0.0_real_8
    !$omp parallel do private(IR) reduction(+:SUMR)
    DO ir=1,fpar%nnr1
       sumr=sumr+rhof(ir)
    ENDDO
    CALL mp_sum(sumr,parai%allgrp)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Norm of the perturbating density:',&
            sumr*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! calculates the corresponding perturbating potential
    CALL zeroing(vtemp)!,nhg)
    CALL v1ofrho_p(vtemp,eirop1,rhof,rhoe,psif,.FALSE.)
    CALL zeroing(vf)!,maxfftn)
    !$omp parallel do private(ir)
    DO ir=1,fpar%nnr1
       vf(ir)=CMPLX(rhoe(ir,1),0.0_real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(vf,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       v1_loc(ig)=vf(nzh(ig))
    ENDDO

    CALL rhoofr(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==         the basic loops for perturbations                    ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",22("*"),a,22("*"),/)')'   perturbations    '
    ENDIF
    ! ...  optimization for c1
    CALL zeroing(c1)!,ngw*nstate)
    ! ==--------------------------------------------------------------==
    CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,&
         eivps,v1_loc,v1_nonloc,z11,nstate,eirop1)
    ! ==--------------------------------------------------------------==
    ! 
    ! ==--------------------------------------------------------------==
    ! calculate the linear perturbation in the density 
    ! \   drhoe=<c0|r><r|c1>+cc (from rwfopt_p) and rhof=rhof+drhoe
    CALL dcopy(fpar%nnr1,drhoe,1,rhof,1)
    DO i=nbeg,nend
       CALL zeroing(psif)!,maxfftn)
       DO ig=1,ncpw%ngw
          psif(nzh(ig))=cf(ig,i)
          psif(indz(ig))=CONJG(cf(ig,i))
       ENDDO
       IF (geq0) psif(nzh(1))=cf(1,i)
       CALL  invfftn(psif,.FALSE.,parai%allgrp)
       !$omp parallel do private(IR,R1)
       DO ir=1,fpar%nnr1
          r1=REAL(psif(ir))
          rhof(ir)=rhof(ir)+r1*r1/parm%omega*wghtf(i)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! computes global hardness
    ! ==--------------------------------------------------------------==
    CALL zeroing(psif)!,maxfftn)
    CALL  invfftn(vf,.FALSE.,parai%allgrp)
    hardness=0.0_real_8
    DO ir=1,fpar%nnr1
       r1=REAL(vf(ir))
       hardness=hardness+rhof(ir)*r1
       psif(ir)=CMPLX(rhof(ir),0.0_real_8,kind=real_8)
    ENDDO
    CALL mp_sum(hardness,parai%allgrp)
    hardness=hardness*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) 'HARDNESS (eV/e) = deps_n/dn = ',2*ry*hardness
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    CALL  fwfftn(psif,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtemp(ig)=psif(nzh(ig))
    ENDDO
    filen='FUKUI'
    CALL densto(vtemp,tau0,filen)

    ! ==--------------------------------------------------------------==
    ! computes nuclear Fukui functions

    CALL zeroing(fion_p)!,3*ions1%nat)
    CALL rnlsm(c1,nstate,1,1,.TRUE.)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             icol=3*(iat-1)+k
             CALL d_mat_locps(fion_p(icol),vtemp,&
                  iat,k,is)
             CALL d_mat_nonloc(fion_p(icol),fnl,dfnl,fnl00,dfnl00,&
                  crge%f,wk,is,k,k,iat,nstate,1)
          ENDDO
       ENDDO
    ENDDO

    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          ffmat(1,ia,is)=ffmat(1,ia,is)+fion_p(3*iat+1)
          ffmat(2,ia,is)=ffmat(2,ia,is)+fion_p(3*iat+2)
          ffmat(3,ia,is)=ffmat(3,ia,is)+fion_p(3*iat+3)
          iat=iat+1
       ENDDO
    ENDDO

    CALL mp_sum(ffmat,3*maxsys%nax*maxsys%nsx,parai%allgrp)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Nuclear Fukui functions:'
       CALL wrgeof(tau0,ffmat)
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! computes "per state" hardness
    CALL zeroing(vtemp)!,nhg)
    CALL v1ofrho_p(vtemp,eirop1,rhof,rhoe,psif,.FALSE.)

    ! potential * orbitals in c1
    CALL zeroing(c1)!,nstate*ngw)
    CALL v1psi0_p(c0,c1,rhoe,psi,nstate)

    ! overlap
    CALL ovlap(nstate,hard,c0,c1)
    CALL mp_sum(hard,nstate*nstate,parai%allgrp)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'Per state hardness (eV/e) '
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,*) ' STATE     deps_i/dn'
       DO i=1,nstate
          IF (paral%io_parent)&
               WRITE(6,'(I4,2x,F15.6)') nstate+1-i,-1.0_real_8*ry*hard(i,i)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    DEALLOCATE(cf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhof,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xcf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psif,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hard,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion_p,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ffmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('   fukui_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fukui_p

END MODULE fukui_p_utils
