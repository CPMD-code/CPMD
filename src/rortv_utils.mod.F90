MODULE rortv_utils
  USE crotwf_utils,                    ONLY: crotwf,&
                                             give_scr_crotwf
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: dtan2w,&
                                             xmu
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE linalg_utils,                    ONLY: symm_da
  USE mp_interface,                    ONLY: mp_sum
  USE nort,                            ONLY: nort_com
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             parap,&
                                             paraw
  USE tpar,                            ONLY: dtb2me
!!use ovlap_utils, only : ovlap2
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rortv
  !public :: rvharm
  PUBLIC :: give_scr_rortv

CONTAINS

  ! ==================================================================
  SUBROUTINE rortv(c0,cm,c2,sc0,gamy,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(ncpw%ngw,*), &
                                                sc0(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gamy(nstate,nstate)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rortv'

    INTEGER                                  :: i, ierr, ig, ip, j, nstx, nx
    REAL(real_8)                             :: ai, bi, pf4, s2, yi
    REAL(real_8), ALLOCATABLE                :: a1mat(:,:), a2mat(:,:)

    IF (cntl%nonort) THEN
       ! Norm constraints
       DO i=1,nstate
          IF (cntl%tharm.OR.cntl%tmass) THEN
             s2=0.0_real_8
             IF (cntl%tharm) THEN
                !$omp parallel do private(IG,PF4,AI,BI) reduction(+:S2)
                DO ig=1,ncpw%ngw
                   pf4=2.0_real_8*dtan2w(ig)/dtb2me
                   ai=REAL(c0(ig,i))
                   bi=AIMAG(c0(ig,i))
                   s2=s2+pf4*(ai*ai+bi*bi)
                ENDDO
             ELSE
                !$omp parallel do private(IG,PF4,AI,BI) reduction(+:S2)
                DO ig=1,ncpw%ngw
                   pf4=2.0_real_8*cntr%emass/xmu(ig)
                   ai=REAL(c0(ig,i))
                   bi=AIMAG(c0(ig,i))
                   s2=s2+pf4*(ai*ai+bi*bi)
                ENDDO
             ENDIF
             CALL mp_sum(s2,parai%allgrp)
             IF (geq0) THEN
                IF (cntl%tharm) THEN
                   pf4=dtan2w(1)/dtb2me
                ELSE
                   pf4=cntr%emass/xmu(1)
                ENDIF
             ENDIF
             yi=-dotp(ncpw%ngw,cm(:,i),c0(:,i))/s2
             CALL mp_sum(yi,parai%allgrp)
             IF (cntl%tharm) THEN
                !$omp parallel do private(IG,PF4)
                DO ig=1,ncpw%ngw
                   pf4=dtan2w(ig)/dtb2me
                   cm(ig,i)=cm(ig,i)+pf4*yi*c0(ig,i)
                ENDDO
             ELSE
                !$omp parallel do private(IG,PF4)
                DO ig=1,ncpw%ngw
                   pf4=cntr%emass/xmu(ig)
                   cm(ig,i)=cm(ig,i)+pf4*yi*c0(ig,i)
                ENDDO
             ENDIF
          ELSE
             yi=-dotp(ncpw%ngw,cm(:,i),c0(:,i))
             CALL mp_sum(yi,parai%allgrp)
             CALL daxpy(2*ncpw%ngw,yi,c0(1,i),1,cm(1,i),1)
          ENDIF
       ENDDO
       ! UNITARY ROTATION?
       IF (nort_com%scond.GT.nort_com%slimit) THEN
          CALL crotwf(c0,cm,c2,sc0,nstate,gamy)
       ENDIF
    ELSE
       ! Overlap constraints
       IF (cntl%tharm.OR.cntl%tmass) THEN
          CALL rvharm(c0,cm,gamy,nstate)
       ELSE
          IF (cntl%tdmal) THEN
             CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
             ALLOCATE(a1mat(nstate, nstx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)! TODO ALIGN FOR BG
             ALLOCATE(a2mat(nstate, nstx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             CALL zeroing(gamy(:,1:nstx))!,nstate*nstx)
             DO ip=0,parai%nproc-1
                nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
                CALL zeroing(a1mat)!,nstate*nstx)
                IF (nx.GT.0) THEN
                   CALL ovlap2(ncpw%ngw,nstate,nx,a1mat,cm,c0(1,paraw%nwa12(ip,1)),&
                        .TRUE.)
                   CALL mp_sum(a1mat,a2mat,&
                        nstate*nstx,parap%pgroup(ip+1),parai%allgrp)
                   IF (parai%mepos.EQ.parap%pgroup(ip+1)) THEN
                      CALL dcopy(nstate*nstx,a2mat,1,gamy,1)
                   ENDIF
                ENDIF
             ENDDO
             CALL dscal(nstate*nstx,-1.0_real_8,gamy,1)
             CALL symm_da(gamy,a1mat,a2mat,nstate,paraw%nwa12(0,1),paraw%nwa12(0,2),&
                  nstx,parai%mepos,parai%nproc,parai%allgrp)
             CALL rotate_da(1._real_8,c0(1,1),1._real_8,cm(1,1),gamy,2*ncpw%ngw,2*ncpw%ngw,&
                  nstate,paraw%nwa12(0,1),paraw%nwa12(0,2),nstx,parai%mepos,&
                  parap%pgroup,parai%nproc,parai%allgrp,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
             DEALLOCATE(a1mat,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(a2mat,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ELSE
             CALL ovlap(nstate,gamy,cm,c0)
             CALL mp_sum(gamy,nstate*nstate,parai%allgrp)
             !$omp parallel do private(I,J) schedule(static,1)
             DO i=1,nstate
                DO j=i+1,nstate
                   gamy(i,j)=-0.5_real_8*(gamy(i,j)+gamy(j,i))
                   gamy(j,i)=gamy(i,j)
                ENDDO
                gamy(i,i)=-gamy(i,i)
             ENDDO
             CALL rotate(1.0_real_8,c0,1.0_real_8,cm,gamy,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
                  spin_mod%nsdown)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rortv
  ! ==================================================================
#if defined(__VECTOR) && (! (__SR8000))
  ! ==================================================================
  SUBROUTINE rvharm(c0,cm,gamy,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gamy(nstate,*)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rvharm'

    INTEGER                                  :: i, ierr, ig, iter, j, jmax, &
                                                modi, modj
    REAL(real_8) :: ai, aj, aj0, aj1, aj2, aj3, aj4, aj5, aj6, aj7, bi, bj, &
      bj0, bj1, bj2, bj3, bj4, bj5, bj6, bj7, dgam, dgmax, fdtb2me, gij, &
      gij00, gij01, gij02, gij03, gij04, gij05, gij06, gij07, gij10, gij11, &
      gij12, gij13, gij14, gij15, gij16, gij17, gij20, gij21, gij22, gij23, &
      gij24, gij25, gij26, gij27, gij30, gij31, gij32, gij33, gij34, gij35, &
      gij36, gij37, gij40, gij41, gij42, gij43, gij44, gij45, gij46, gij47, &
      gij50, gij51, gij52, gij53, gij54, gij55, gij56, gij57, gij60, gij61, &
      gij62, gij63, gij64, gij65, gij66, gij67, gij70, gij71, gij72, gij73, &
      gij74, gij75, gij76, gij77, pf4, t0, t1
    REAL(real_8) :: t2, t3, t4, t5, t6, t7, temp, tol, wdtb2me
    REAL(real_8), ALLOCATABLE                :: s1(:,:), s2(:,:), s3(:,:), &
                                                s4(:,:)

    ALLOCATE(s1(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO align for BG
    ALLOCATE(s2(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(s3(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(s4(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    tol = cntr%epsog
    fdtb2me=1.0_real_8/dtb2me
    wdtb2me=2.0_real_8*fdtb2me
    CALL ovlap(nstate,s1,cm,c0)
    CALL mp_sum(s1,nstate*nstate,parai%allgrp)
    CALL zeroing(s2)!,nstate*nstate)
    IF (cntl%tharm) THEN
       pf4=dtan2w(1)*fdtb2me
       !$omp parallel do private(I,J,IG,JMAX,AI,BI,AJ,BJ,TEMP,MODJ, &
       !$omp  AJ0,AJ1,AJ2,AJ3,AJ4,AJ5,AJ6,AJ7, &
       !$omp  BJ0,BJ1,BJ2,BJ3,BJ4,BJ5,BJ6,BJ7, &
       !$omp  T0,T1,T2,T3,T4,T5,T6,T7) &
       !$omp  schedule(static,1)
       DO i=1,nstate
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          modj=MOD(jmax-i+1,8)
          DO j=i,i+modj-1
             temp=0.0_real_8
             DO ig=1,ncpw%ngw
                ai=REAL(c0(ig,i))
                bi=AIMAG(c0(ig,i))
                aj=REAL(c0(ig,j))
                bj=AIMAG(c0(ig,j))
                temp=temp+dtan2w(ig)*(ai*aj+bi*bj)
             ENDDO
             s2(i,j)=wdtb2me*temp
             IF (geq0) s2(i,j)=s2(i,j)-pf4*REAL(c0(1,i))*REAL(c0(1,j))
             s2(j,i)=s2(i,j)
          ENDDO
          DO j=i+modj,jmax,8
             t0=0.0_real_8
             t1=0.0_real_8
             t2=0.0_real_8
             t3=0.0_real_8
             t4=0.0_real_8
             t5=0.0_real_8
             t6=0.0_real_8
             t7=0.0_real_8
             DO ig=1,ncpw%ngw
                ai =REAL(c0(ig,i  ))
                bi =AIMAG(c0(ig,i  ))
                aj0=REAL(c0(ig,j  ))
                aj1=REAL(c0(ig,j+1))
                aj2=REAL(c0(ig,j+2))
                aj3=REAL(c0(ig,j+3))
                aj4=REAL(c0(ig,j+4))
                aj5=REAL(c0(ig,j+5))
                aj6=REAL(c0(ig,j+6))
                aj7=REAL(c0(ig,j+7))
                bj0=AIMAG(c0(ig,j  ))
                bj1=AIMAG(c0(ig,j+1))
                bj2=AIMAG(c0(ig,j+2))
                bj3=AIMAG(c0(ig,j+3))
                bj4=AIMAG(c0(ig,j+4))
                bj5=AIMAG(c0(ig,j+5))
                bj6=AIMAG(c0(ig,j+6))
                bj7=AIMAG(c0(ig,j+7))
                t0=t0+dtan2w(ig)*(ai*aj0+bi*bj0)
                t1=t1+dtan2w(ig)*(ai*aj1+bi*bj1)
                t2=t2+dtan2w(ig)*(ai*aj2+bi*bj2)
                t3=t3+dtan2w(ig)*(ai*aj3+bi*bj3)
                t4=t4+dtan2w(ig)*(ai*aj4+bi*bj4)
                t5=t5+dtan2w(ig)*(ai*aj5+bi*bj5)
                t6=t6+dtan2w(ig)*(ai*aj6+bi*bj6)
                t7=t7+dtan2w(ig)*(ai*aj7+bi*bj7)
             ENDDO
             s2(i,j  )=wdtb2me*t0
             s2(i,j+1)=wdtb2me*t1
             s2(i,j+2)=wdtb2me*t2
             s2(i,j+3)=wdtb2me*t3
             s2(i,j+4)=wdtb2me*t4
             s2(i,j+5)=wdtb2me*t5
             s2(i,j+6)=wdtb2me*t6
             s2(i,j+7)=wdtb2me*t7
             IF (geq0) THEN
                temp=pf4*REAL(c0(1,i))
                s2(i,j  )=s2(i,j  )-temp*REAL(c0(1,j  ))
                s2(i,j+1)=s2(i,j+1)-temp*REAL(c0(1,j+1))
                s2(i,j+2)=s2(i,j+2)-temp*REAL(c0(1,j+2))
                s2(i,j+3)=s2(i,j+3)-temp*REAL(c0(1,j+3))
                s2(i,j+4)=s2(i,j+4)-temp*REAL(c0(1,j+4))
                s2(i,j+5)=s2(i,j+5)-temp*REAL(c0(1,j+5))
                s2(i,j+6)=s2(i,j+6)-temp*REAL(c0(1,j+6))
                s2(i,j+7)=s2(i,j+7)-temp*REAL(c0(1,j+7))
             ENDIF
             s2(j  ,i)=s2(i,j  )
             s2(j+1,i)=s2(i,j+1)
             s2(j+2,i)=s2(i,j+2)
             s2(j+3,i)=s2(i,j+3)
             s2(j+4,i)=s2(i,j+4)
             s2(j+5,i)=s2(i,j+5)
             s2(j+6,i)=s2(i,j+6)
             s2(j+7,i)=s2(i,j+7)
          ENDDO
          DO j=jmax+1,nstate
             s2(i,j)=0.0_real_8
             s2(j,i)=0.0_real_8
          ENDDO
       ENDDO
    ELSE
       pf4=cntr%emass/xmu(1)
       !$omp parallel do private(I,J,IG,JMAX,AI,BI,AJ,BJ,TEMP,MODJ, &
       !$omp  AJ0,AJ1,AJ2,AJ3,AJ4,AJ5,AJ6,AJ7, &
       !$omp  BJ0,BJ1,BJ2,BJ3,BJ4,BJ5,BJ6,BJ7, &
       !$omp  T0,T1,T2,T3,T4,T5,T6,T7) &
       !$omp  schedule(static,1)
       DO i=1,nstate
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          modj=MOD(jmax-i+1,8)
          DO j=i,i+modj-1
             temp=0.0_real_8
             DO ig=1,ncpw%ngw
                ai=REAL(c0(ig,i))
                bi=AIMAG(c0(ig,i))
                aj=REAL(c0(ig,j))
                bj=AIMAG(c0(ig,j))
                temp=temp+(ai*aj+bi*bj)/xmu(ig)
             ENDDO
             s2(i,j)=2.0_real_8*cntr%emass*temp
             IF (geq0) s2(i,j)=s2(i,j)-pf4*REAL(c0(1,i))*REAL(c0(1,j))
             s2(j,i)=s2(i,j)
          ENDDO
          DO j=i+modj,jmax,8
             t0=0.0_real_8
             t1=0.0_real_8
             t3=0.0_real_8
             t4=0.0_real_8
             t5=0.0_real_8
             t6=0.0_real_8
             t7=0.0_real_8
             DO ig=1,ncpw%ngw
                ai =REAL(c0(ig,i  ))
                bi =AIMAG(c0(ig,i  ))
                aj0=REAL(c0(ig,j  ))
                aj1=REAL(c0(ig,j+1))
                aj2=REAL(c0(ig,j+2))
                aj3=REAL(c0(ig,j+3))
                aj4=REAL(c0(ig,j+4))
                aj5=REAL(c0(ig,j+5))
                aj6=REAL(c0(ig,j+6))
                aj7=REAL(c0(ig,j+7))
                bj0=AIMAG(c0(ig,j  ))
                bj1=AIMAG(c0(ig,j+1))
                bj2=AIMAG(c0(ig,j+2))
                bj3=AIMAG(c0(ig,j+3))
                bj4=AIMAG(c0(ig,j+4))
                bj5=AIMAG(c0(ig,j+5))
                bj6=AIMAG(c0(ig,j+6))
                bj7=AIMAG(c0(ig,j+7))
                t0=t0+(ai*aj0+bi*bj0)/xmu(ig)
                t1=t1+(ai*aj1+bi*bj1)/xmu(ig)
                t2=t2+(ai*aj2+bi*bj2)/xmu(ig)
                t3=t3+(ai*aj3+bi*bj3)/xmu(ig)
                t4=t4+(ai*aj4+bi*bj4)/xmu(ig)
                t5=t5+(ai*aj5+bi*bj5)/xmu(ig)
                t6=t6+(ai*aj6+bi*bj6)/xmu(ig)
                t7=t7+(ai*aj7+bi*bj7)/xmu(ig)
             ENDDO
             s2(i,j  )=2.0_real_8*cntr%emass*t0
             s2(i,j+1)=2.0_real_8*cntr%emass*t1
             s2(i,j+2)=2.0_real_8*cntr%emass*t2
             s2(i,j+3)=2.0_real_8*cntr%emass*t3
             s2(i,j+4)=2.0_real_8*cntr%emass*t4
             s2(i,j+5)=2.0_real_8*cntr%emass*t5
             s2(i,j+6)=2.0_real_8*cntr%emass*t6
             s2(i,j+7)=2.0_real_8*cntr%emass*t7
             IF (geq0) THEN
                s2(i,j  )=s2(i,j  )-pf4*REAL(c0(1,i))*REAL(c0(1,j  ))
                s2(i,j+1)=s2(i,j+1)-pf4*REAL(c0(1,i))*REAL(c0(1,j+1))
                s2(i,j+2)=s2(i,j+2)-pf4*REAL(c0(1,i))*REAL(c0(1,j+2))
                s2(i,j+3)=s2(i,j+3)-pf4*REAL(c0(1,i))*REAL(c0(1,j+3))
                s2(i,j+4)=s2(i,j+4)-pf4*REAL(c0(1,i))*REAL(c0(1,j+4))
                s2(i,j+5)=s2(i,j+5)-pf4*REAL(c0(1,i))*REAL(c0(1,j+5))
                s2(i,j+6)=s2(i,j+6)-pf4*REAL(c0(1,i))*REAL(c0(1,j+6))
                s2(i,j+7)=s2(i,j+7)-pf4*REAL(c0(1,i))*REAL(c0(1,j+7))
             ENDIF
             s2(j  ,i)=s2(i,j  )
             s2(j+1,i)=s2(i,j+1)
             s2(j+2,i)=s2(i,j+2)
             s2(j+3,i)=s2(i,j+3)
             s2(j+4,i)=s2(i,j+4)
             s2(j+5,i)=s2(i,j+5)
             s2(j+6,i)=s2(i,j+6)
             s2(j+7,i)=s2(i,j+7)
          ENDDO
          DO j=jmax+1,nstate
             s2(i,j)=0.0_real_8
             s2(j,i)=0.0_real_8
          ENDDO
       ENDDO
    ENDIF
    CALL mp_sum(s2,nstate*nstate,parai%allgrp)
    DO i=1,nstate
       s2(i,i)=s2(i,i)-1.0_real_8
    ENDDO
    ! Initial Guess for GAMY
    !$omp parallel do private(I,J)
    DO j=1,nstate
       DO i=1,nstate
          gamy(i,j)=-0.5_real_8*(s1(i,j)+s1(j,i))
          s4(i,j)=gamy(i,j)
       ENDDO
    ENDDO
    ! Iteration
    iter=0
10  CONTINUE
    iter=iter+1
    CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,s2,nstate,gamy,&
         nstate,0.0_real_8,s3,nstate)
    dgmax=0.0_real_8
    !$omp parallel do private(I,J,DGAM) reduction(max:DGMAX)
    DO j=1,nstate
       DO i=1,nstate
          gamy(i,j)=-0.5_real_8*(s1(i,j)+s1(j,i)+s3(i,j)+s3(j,i))
          dgam = (gamy(i,j)-s4(i,j))*(gamy(i,j)-s4(i,j))
          dgmax=MAX(dgmax,dgam)
          s4(i,j)=gamy(i,j)
       ENDDO
    ENDDO
    dgmax=SQRT(dgmax)
    IF (dgmax.GT.tol) GOTO 10
    ! Apply the constraint
    IF (cntl%tharm) THEN
       modi=MOD(nstate,8)
       modj=MOD(nstate,8)
       !$omp parallel do private(I,J,IG,GIJ, &
       !$omp  GIJ00,GIJ01,GIJ02,GIJ03,GIJ04,GIJ05,GIJ06,GIJ07)
       DO i=1,modi
          DO j=1,modj
             gij=gamy(i,j)*fdtb2me
             DO ig=1,ncpw%ngw
                cm(ig,i)=cm(ig,i)+dtan2w(ig)*gij*c0(ig,j)
             ENDDO
          ENDDO
          DO j=modj+1,nstate,8
             gij00=gamy(i,j  )*fdtb2me
             gij01=gamy(i,j+1)*fdtb2me
             gij02=gamy(i,j+2)*fdtb2me
             gij03=gamy(i,j+3)*fdtb2me
             gij04=gamy(i,j+4)*fdtb2me
             gij05=gamy(i,j+5)*fdtb2me
             gij06=gamy(i,j+6)*fdtb2me
             gij07=gamy(i,j+7)*fdtb2me
             DO ig=1,ncpw%ngw
                cm(ig,i  )=cm(ig,i  ) +dtan2w(ig)*(gij00*c0(ig,j  )&
                     +gij01*c0(ig,j+1)+gij02*c0(ig,j+2)+gij03*c0(ig,j+3)&
                     +gij04*c0(ig,j+4)+gij05*c0(ig,j+5)+gij06*c0(ig,j+6)&
                     +gij07*c0(ig,j+7))
             ENDDO
          ENDDO
       ENDDO
       !$omp parallel do private(I,J,IG, &
       !$omp  GIJ00,GIJ01,GIJ02,GIJ03,GIJ04,GIJ05,GIJ06,GIJ07, &
       !$omp  GIJ10,GIJ11,GIJ12,GIJ13,GIJ14,GIJ15,GIJ16,GIJ17, &
       !$omp  GIJ20,GIJ21,GIJ22,GIJ23,GIJ24,GIJ25,GIJ26,GIJ27, &
       !$omp  GIJ30,GIJ31,GIJ32,GIJ33,GIJ34,GIJ35,GIJ36,GIJ37, &
       !$omp  GIJ40,GIJ41,GIJ42,GIJ43,GIJ44,GIJ45,GIJ46,GIJ47, &
       !$omp  GIJ50,GIJ51,GIJ52,GIJ53,GIJ54,GIJ55,GIJ56,GIJ57, &
       !$omp  GIJ60,GIJ61,GIJ62,GIJ63,GIJ64,GIJ65,GIJ66,GIJ67, &
       !$omp  GIJ70,GIJ71,GIJ72,GIJ73,GIJ74,GIJ75,GIJ76,GIJ77)
       DO i=modi+1,nstate,8
          DO j=1,modj
             gij00=gamy(i  ,j)*fdtb2me
             gij10=gamy(i+1,j)*fdtb2me
             gij20=gamy(i+2,j)*fdtb2me
             gij30=gamy(i+3,j)*fdtb2me
             gij40=gamy(i+4,j)*fdtb2me
             gij50=gamy(i+5,j)*fdtb2me
             gij60=gamy(i+6,j)*fdtb2me
             gij70=gamy(i+7,j)*fdtb2me
             DO ig=1,ncpw%ngw
                cm(ig,i  )=cm(ig,i  )+dtan2w(ig)*gij00*c0(ig,j)
                cm(ig,i+1)=cm(ig,i+1)+dtan2w(ig)*gij10*c0(ig,j)
                cm(ig,i+2)=cm(ig,i+2)+dtan2w(ig)*gij20*c0(ig,j)
                cm(ig,i+3)=cm(ig,i+3)+dtan2w(ig)*gij30*c0(ig,j)
                cm(ig,i+4)=cm(ig,i+4)+dtan2w(ig)*gij40*c0(ig,j)
                cm(ig,i+5)=cm(ig,i+5)+dtan2w(ig)*gij50*c0(ig,j)
                cm(ig,i+6)=cm(ig,i+6)+dtan2w(ig)*gij60*c0(ig,j)
                cm(ig,i+7)=cm(ig,i+7)+dtan2w(ig)*gij70*c0(ig,j)
             ENDDO
          ENDDO
          DO j=modj+1,nstate,8
             gij00=gamy(i  ,j  )*fdtb2me
             gij10=gamy(i+1,j  )*fdtb2me
             gij20=gamy(i+2,j  )*fdtb2me
             gij30=gamy(i+3,j  )*fdtb2me
             gij40=gamy(i+4,j  )*fdtb2me
             gij50=gamy(i+5,j  )*fdtb2me
             gij60=gamy(i+6,j  )*fdtb2me
             gij70=gamy(i+7,j  )*fdtb2me
             gij01=gamy(i  ,j+1)*fdtb2me
             gij11=gamy(i+1,j+1)*fdtb2me
             gij21=gamy(i+2,j+1)*fdtb2me
             gij31=gamy(i+3,j+1)*fdtb2me
             gij41=gamy(i+4,j+1)*fdtb2me
             gij51=gamy(i+5,j+1)*fdtb2me
             gij61=gamy(i+6,j+1)*fdtb2me
             gij71=gamy(i+7,j+1)*fdtb2me
             gij02=gamy(i  ,j+2)*fdtb2me
             gij12=gamy(i+1,j+2)*fdtb2me
             gij22=gamy(i+2,j+2)*fdtb2me
             gij32=gamy(i+3,j+2)*fdtb2me
             gij42=gamy(i+4,j+2)*fdtb2me
             gij52=gamy(i+5,j+2)*fdtb2me
             gij62=gamy(i+6,j+2)*fdtb2me
             gij72=gamy(i+7,j+2)*fdtb2me
             gij03=gamy(i  ,j+3)*fdtb2me
             gij13=gamy(i+1,j+3)*fdtb2me
             gij23=gamy(i+2,j+3)*fdtb2me
             gij33=gamy(i+3,j+3)*fdtb2me
             gij43=gamy(i+4,j+3)*fdtb2me
             gij53=gamy(i+5,j+3)*fdtb2me
             gij63=gamy(i+6,j+3)*fdtb2me
             gij73=gamy(i+7,j+3)*fdtb2me
             gij04=gamy(i  ,j+4)*fdtb2me
             gij14=gamy(i+1,j+4)*fdtb2me
             gij24=gamy(i+2,j+4)*fdtb2me
             gij34=gamy(i+3,j+4)*fdtb2me
             gij44=gamy(i+4,j+4)*fdtb2me
             gij54=gamy(i+5,j+4)*fdtb2me
             gij64=gamy(i+6,j+4)*fdtb2me
             gij74=gamy(i+7,j+4)*fdtb2me
             gij05=gamy(i  ,j+5)*fdtb2me
             gij15=gamy(i+1,j+5)*fdtb2me
             gij25=gamy(i+2,j+5)*fdtb2me
             gij35=gamy(i+3,j+5)*fdtb2me
             gij45=gamy(i+4,j+5)*fdtb2me
             gij55=gamy(i+5,j+5)*fdtb2me
             gij65=gamy(i+6,j+5)*fdtb2me
             gij75=gamy(i+7,j+5)*fdtb2me
             gij06=gamy(i  ,j+6)*fdtb2me
             gij16=gamy(i+1,j+6)*fdtb2me
             gij26=gamy(i+2,j+6)*fdtb2me
             gij36=gamy(i+3,j+6)*fdtb2me
             gij46=gamy(i+4,j+6)*fdtb2me
             gij56=gamy(i+5,j+6)*fdtb2me
             gij66=gamy(i+6,j+6)*fdtb2me
             gij76=gamy(i+7,j+6)*fdtb2me
             gij07=gamy(i  ,j+7)*fdtb2me
             gij17=gamy(i+1,j+7)*fdtb2me
             gij27=gamy(i+2,j+7)*fdtb2me
             gij37=gamy(i+3,j+7)*fdtb2me
             gij47=gamy(i+4,j+7)*fdtb2me
             gij57=gamy(i+5,j+7)*fdtb2me
             gij67=gamy(i+6,j+7)*fdtb2me
             gij77=gamy(i+7,j+7)*fdtb2me
             DO ig=1,ncpw%ngw
                cm(ig,i  )=cm(ig,i  ) +dtan2w(ig)*(gij00*c0(ig,j  )&
                     +gij01*c0(ig,j+1)+gij02*c0(ig,j+2)+gij03*c0(ig,j+3)&
                     +gij04*c0(ig,j+4)+gij05*c0(ig,j+5)+gij06*c0(ig,j+6)&
                     +gij07*c0(ig,j+7))
                cm(ig,i+1)=cm(ig,i+1) +dtan2w(ig)*(gij10*c0(ig,j  )&
                     +gij11*c0(ig,j+1)+gij12*c0(ig,j+2)+gij13*c0(ig,j+3)&
                     +gij14*c0(ig,j+4)+gij15*c0(ig,j+5)+gij16*c0(ig,j+6)&
                     +gij17*c0(ig,j+7))
                cm(ig,i+2)=cm(ig,i+2) +dtan2w(ig)*(gij20*c0(ig,j  )&
                     +gij21*c0(ig,j+1)+gij22*c0(ig,j+2)+gij23*c0(ig,j+3)&
                     +gij24*c0(ig,j+4)+gij25*c0(ig,j+5)+gij26*c0(ig,j+6)&
                     +gij27*c0(ig,j+7))
                cm(ig,i+3)=cm(ig,i+3) +dtan2w(ig)*(gij30*c0(ig,j  )&
                     +gij31*c0(ig,j+1)+gij32*c0(ig,j+2)+gij33*c0(ig,j+3)&
                     +gij34*c0(ig,j+4)+gij35*c0(ig,j+5)+gij36*c0(ig,j+6)&
                     +gij37*c0(ig,j+7))
                cm(ig,i+4)=cm(ig,i+4) +dtan2w(ig)*(gij40*c0(ig,j  )&
                     +gij41*c0(ig,j+1)+gij42*c0(ig,j+2)+gij43*c0(ig,j+3)&
                     +gij44*c0(ig,j+4)+gij45*c0(ig,j+5)+gij46*c0(ig,j+6)&
                     +gij47*c0(ig,j+7))
                cm(ig,i+5)=cm(ig,i+5) +dtan2w(ig)*(gij50*c0(ig,j  )&
                     +gij51*c0(ig,j+1)+gij52*c0(ig,j+2)+gij53*c0(ig,j+3)&
                     +gij54*c0(ig,j+4)+gij55*c0(ig,j+5)+gij56*c0(ig,j+6)&
                     +gij57*c0(ig,j+7))
                cm(ig,i+6)=cm(ig,i+6) +dtan2w(ig)*(gij60*c0(ig,j  )&
                     +gij61*c0(ig,j+1)+gij62*c0(ig,j+2)+gij63*c0(ig,j+3)&
                     +gij64*c0(ig,j+4)+gij65*c0(ig,j+5)+gij66*c0(ig,j+6)&
                     +gij67*c0(ig,j+7))
                cm(ig,i+7)=cm(ig,i+7) +dtan2w(ig)*(gij70*c0(ig,j  )&
                     +gij71*c0(ig,j+1)+gij72*c0(ig,j+2)+gij73*c0(ig,j+3)&
                     +gij74*c0(ig,j+4)+gij75*c0(ig,j+5)+gij76*c0(ig,j+6)&
                     +gij77*c0(ig,j+7))
             ENDDO
          ENDDO
       ENDDO
    ELSE
       modi=MOD(nstate,8)
       modj=MOD(nstate,8)
       !$omp parallel do private(I,J,IG,GIJ, &
       !$omp  GIJ00,GIJ01,GIJ02,GIJ03,GIJ04,GIJ05,GIJ06,GIJ07)
       DO i=1,modi
          DO j=1,modj
             gij=cntr%emass*gamy(i,j)
             DO ig=1,ncpw%ngw
                cm(ig,i)=cm(ig,i)+gij*c0(ig,j)/xmu(ig)
             ENDDO
          ENDDO
          DO j=modj+1,nstate,8
             gij00=cntr%emass*gamy(i,j  )
             gij01=cntr%emass*gamy(i,j+1)
             gij02=cntr%emass*gamy(i,j+2)
             gij03=cntr%emass*gamy(i,j+3)
             gij04=cntr%emass*gamy(i,j+4)
             gij05=cntr%emass*gamy(i,j+5)
             gij06=cntr%emass*gamy(i,j+6)
             gij07=cntr%emass*gamy(i,j+7)
             DO ig=1,ncpw%ngw
                cm(ig,i)=cm(ig,i)+(gij00*c0(ig,j  )&
                     +gij01*c0(ig,j+1)+gij02*c0(ig,j+2)+gij03*c0(ig,j+3)&
                     +gij04*c0(ig,j+4)+gij05*c0(ig,j+5)+gij06*c0(ig,j+6)&
                     +gij07*c0(ig,j+7))/xmu(ig)
             ENDDO
          ENDDO
       ENDDO
       !$omp parallel do private(I,J,IG, &
       !$omp  GIJ00,GIJ01,GIJ02,GIJ03,GIJ04,GIJ05,GIJ06,GIJ07, &
       !$omp  GIJ10,GIJ11,GIJ12,GIJ13,GIJ14,GIJ15,GIJ16,GIJ17, &
       !$omp  GIJ20,GIJ21,GIJ22,GIJ23,GIJ24,GIJ25,GIJ26,GIJ27, &
       !$omp  GIJ30,GIJ31,GIJ32,GIJ33,GIJ34,GIJ35,GIJ36,GIJ37, &
       !$omp  GIJ40,GIJ41,GIJ42,GIJ43,GIJ44,GIJ45,GIJ46,GIJ47, &
       !$omp  GIJ50,GIJ51,GIJ52,GIJ53,GIJ54,GIJ55,GIJ56,GIJ57, &
       !$omp  GIJ60,GIJ61,GIJ62,GIJ63,GIJ64,GIJ65,GIJ66,GIJ67, &
       !$omp  GIJ70,GIJ71,GIJ72,GIJ73,GIJ74,GIJ75,GIJ76,GIJ77)
       DO i=modi+1,nstate,8
          DO j=1,modj
             gij00=cntr%emass*gamy(i  ,j)
             gij10=cntr%emass*gamy(i+1,j)
             gij20=cntr%emass*gamy(i+2,j)
             gij30=cntr%emass*gamy(i+3,j)
             gij40=cntr%emass*gamy(i+4,j)
             gij50=cntr%emass*gamy(i+5,j)
             gij60=cntr%emass*gamy(i+6,j)
             gij70=cntr%emass*gamy(i+7,j)
             DO ig=1,ncpw%ngw
                cm(ig,i  )=cm(ig,i  )+gij00*c0(ig,j)/xmu(ig)
                cm(ig,i+1)=cm(ig,i+1)+gij10*c0(ig,j)/xmu(ig)
                cm(ig,i+2)=cm(ig,i+2)+gij20*c0(ig,j)/xmu(ig)
                cm(ig,i+3)=cm(ig,i+3)+gij30*c0(ig,j)/xmu(ig)
                cm(ig,i+4)=cm(ig,i+4)+gij40*c0(ig,j)/xmu(ig)
                cm(ig,i+5)=cm(ig,i+5)+gij50*c0(ig,j)/xmu(ig)
                cm(ig,i+6)=cm(ig,i+6)+gij60*c0(ig,j)/xmu(ig)
                cm(ig,i+7)=cm(ig,i+7)+gij70*c0(ig,j)/xmu(ig)
             ENDDO
          ENDDO
          DO j=modj+1,nstate,8
             gij00=cntr%emass*gamy(i  ,j  )
             gij10=cntr%emass*gamy(i+1,j  )
             gij20=cntr%emass*gamy(i+2,j  )
             gij30=cntr%emass*gamy(i+3,j  )
             gij40=cntr%emass*gamy(i+4,j  )
             gij50=cntr%emass*gamy(i+5,j  )
             gij60=cntr%emass*gamy(i+6,j  )
             gij70=cntr%emass*gamy(i+7,j  )
             gij01=cntr%emass*gamy(i  ,j+1)
             gij11=cntr%emass*gamy(i+1,j+1)
             gij21=cntr%emass*gamy(i+2,j+1)
             gij31=cntr%emass*gamy(i+3,j+1)
             gij41=cntr%emass*gamy(i+4,j+1)
             gij51=cntr%emass*gamy(i+5,j+1)
             gij61=cntr%emass*gamy(i+6,j+1)
             gij71=cntr%emass*gamy(i+7,j+1)
             gij02=cntr%emass*gamy(i  ,j+2)
             gij12=cntr%emass*gamy(i+1,j+2)
             gij22=cntr%emass*gamy(i+2,j+2)
             gij32=cntr%emass*gamy(i+3,j+2)
             gij42=cntr%emass*gamy(i+4,j+2)
             gij52=cntr%emass*gamy(i+5,j+2)
             gij62=cntr%emass*gamy(i+6,j+2)
             gij72=cntr%emass*gamy(i+7,j+2)
             gij03=cntr%emass*gamy(i  ,j+3)
             gij13=cntr%emass*gamy(i+1,j+3)
             gij23=cntr%emass*gamy(i+2,j+3)
             gij33=cntr%emass*gamy(i+3,j+3)
             gij43=cntr%emass*gamy(i+4,j+3)
             gij53=cntr%emass*gamy(i+5,j+3)
             gij63=cntr%emass*gamy(i+6,j+3)
             gij73=cntr%emass*gamy(i+7,j+3)
             gij04=cntr%emass*gamy(i  ,j+4)
             gij14=cntr%emass*gamy(i+1,j+4)
             gij24=cntr%emass*gamy(i+2,j+4)
             gij34=cntr%emass*gamy(i+3,j+4)
             gij44=cntr%emass*gamy(i+4,j+4)
             gij54=cntr%emass*gamy(i+5,j+4)
             gij64=cntr%emass*gamy(i+6,j+4)
             gij74=cntr%emass*gamy(i+7,j+4)
             gij05=cntr%emass*gamy(i  ,j+5)
             gij15=cntr%emass*gamy(i+1,j+5)
             gij25=cntr%emass*gamy(i+2,j+5)
             gij35=cntr%emass*gamy(i+3,j+5)
             gij45=cntr%emass*gamy(i+4,j+5)
             gij55=cntr%emass*gamy(i+5,j+5)
             gij65=cntr%emass*gamy(i+6,j+5)
             gij75=cntr%emass*gamy(i+7,j+5)
             gij06=cntr%emass*gamy(i  ,j+6)
             gij16=cntr%emass*gamy(i+1,j+6)
             gij26=cntr%emass*gamy(i+2,j+6)
             gij36=cntr%emass*gamy(i+3,j+6)
             gij46=cntr%emass*gamy(i+4,j+6)
             gij56=cntr%emass*gamy(i+5,j+6)
             gij66=cntr%emass*gamy(i+6,j+6)
             gij76=cntr%emass*gamy(i+7,j+6)
             gij07=cntr%emass*gamy(i  ,j+7)
             gij17=cntr%emass*gamy(i+1,j+7)
             gij27=cntr%emass*gamy(i+2,j+7)
             gij37=cntr%emass*gamy(i+3,j+7)
             gij47=cntr%emass*gamy(i+4,j+7)
             gij57=cntr%emass*gamy(i+5,j+7)
             gij67=cntr%emass*gamy(i+6,j+7)
             gij77=cntr%emass*gamy(i+7,j+7)
             DO ig=1,ncpw%ngw
                cm(ig,i  )=cm(ig,i  )+(gij00*c0(ig,j  )&
                     +gij01*c0(ig,j+1)+gij02*c0(ig,j+2)+gij03*c0(ig,j+3)&
                     +gij04*c0(ig,j+4)+gij05*c0(ig,j+5)+gij06*c0(ig,j+6)&
                     +gij07*c0(ig,j+7))/xmu(ig)
                cm(ig,i+1)=cm(ig,i+1)+(gij10*c0(ig,j  )&
                     +gij11*c0(ig,j+1)+gij12*c0(ig,j+2)+gij13*c0(ig,j+3)&
                     +gij14*c0(ig,j+4)+gij15*c0(ig,j+5)+gij16*c0(ig,j+6)&
                     +gij17*c0(ig,j+7))/xmu(ig)
                cm(ig,i+2)=cm(ig,i+2)+(gij20*c0(ig,j  )&
                     +gij21*c0(ig,j+1)+gij22*c0(ig,j+2)+gij23*c0(ig,j+3)&
                     +gij24*c0(ig,j+4)+gij25*c0(ig,j+5)+gij26*c0(ig,j+6)&
                     +gij27*c0(ig,j+7))/xmu(ig)
                cm(ig,i+3)=cm(ig,i+3)+(gij30*c0(ig,j  )&
                     +gij31*c0(ig,j+1)+gij32*c0(ig,j+2)+gij33*c0(ig,j+3)&
                     +gij34*c0(ig,j+4)+gij35*c0(ig,j+5)+gij36*c0(ig,j+6)&
                     +gij37*c0(ig,j+7))/xmu(ig)
                cm(ig,i+4)=cm(ig,i+4)+(gij40*c0(ig,j  )&
                     +gij41*c0(ig,j+1)+gij42*c0(ig,j+2)+gij43*c0(ig,j+3)&
                     +gij44*c0(ig,j+4)+gij45*c0(ig,j+5)+gij46*c0(ig,j+6)&
                     +gij47*c0(ig,j+7))/xmu(ig)
                cm(ig,i+5)=cm(ig,i+5)+(gij50*c0(ig,j  )&
                     +gij51*c0(ig,j+1)+gij52*c0(ig,j+2)+gij53*c0(ig,j+3)&
                     +gij54*c0(ig,j+4)+gij55*c0(ig,j+5)+gij56*c0(ig,j+6)&
                     +gij57*c0(ig,j+7))/xmu(ig)
                cm(ig,i+6)=cm(ig,i+6)+(gij60*c0(ig,j  )&
                     +gij61*c0(ig,j+1)+gij62*c0(ig,j+2)+gij63*c0(ig,j+3)&
                     +gij64*c0(ig,j+4)+gij65*c0(ig,j+5)+gij66*c0(ig,j+6)&
                     +gij67*c0(ig,j+7))/xmu(ig)
                cm(ig,i+7)=cm(ig,i+7)+(gij70*c0(ig,j  )&
                     +gij71*c0(ig,j+1)+gij72*c0(ig,j+2)+gij73*c0(ig,j+3)&
                     +gij74*c0(ig,j+4)+gij75*c0(ig,j+5)+gij76*c0(ig,j+6)&
                     +gij77*c0(ig,j+7))/xmu(ig)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(s1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s4,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rvharm
  ! ==================================================================
#else 
  ! ==================================================================
  SUBROUTINE rvharm(c0,cm,gamy,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gamy(nstate,*)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rvharm'

    INTEGER                                  :: i, ierr, ig, iter, j, jmax
    REAL(real_8)                             :: ai, aj, bi, bj, dgam, dgmax, &
                                                gij, pf4, tol
    REAL(real_8), ALLOCATABLE                :: s1(:,:), s2(:,:), s3(:,:), &
                                                s4(:,:)

    ALLOCATE(s1(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO align for BG
    ALLOCATE(s2(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(s3(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(s4(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    tol = cntr%epsog
    CALL ovlap(nstate,s1,cm,c0)
    CALL mp_sum(s1,nstate*nstate,parai%allgrp)
    CALL zeroing(s2)!,nstate*nstate)
    !$omp parallel do private(I,J,IG,JMAX,PF4,AI,BI,AJ,BJ) &
    !$omp  schedule(static,1)
#ifdef __SR8000
    !poption parallel
    !poption tlocal(I,J,IG,JMAX,PF4,AI,BI,AJ,BJ), cyclic
#endif
    DO i=1,nstate
       jmax=nstate
       IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
       IF (cntl%tharm) THEN
          DO j=i,jmax
             DO ig=1,ncpw%ngw
                pf4=2.0_real_8*dtan2w(ig)/dtb2me
                ai=REAL(c0(ig,i))
                bi=AIMAG(c0(ig,i))
                aj=REAL(c0(ig,j))
                bj=AIMAG(c0(ig,j))
                s2(i,j)=s2(i,j)+pf4*(ai*aj+bi*bj)
             ENDDO
             IF (geq0) THEN
                pf4=dtan2w(1)/dtb2me
                s2(i,j)=s2(i,j)-pf4*REAL(c0(1,i))*REAL(c0(1,j))
             ENDIF
             s2(j,i)=s2(i,j)
          ENDDO
       ELSE
          DO j=i,jmax
             DO ig=1,ncpw%ngw
                pf4=2.0_real_8*cntr%emass/xmu(ig)
                ai=REAL(c0(ig,i))
                bi=AIMAG(c0(ig,i))
                aj=REAL(c0(ig,j))
                bj=AIMAG(c0(ig,j))
                s2(i,j)=s2(i,j)+pf4*(ai*aj+bi*bj)
             ENDDO
             IF (geq0) THEN
                pf4=cntr%emass/xmu(1)
                s2(i,j)=s2(i,j)-pf4*REAL(c0(1,i))*REAL(c0(1,j))
             ENDIF
             s2(j,i)=s2(i,j)
          ENDDO
       ENDIF
       DO j=jmax+1,nstate
          s2(i,j)=0.0_real_8
          s2(j,i)=0.0_real_8
       ENDDO
    ENDDO
    CALL mp_sum(s2,nstate*nstate,parai%allgrp)
    DO i=1,nstate
       s2(i,i)=s2(i,i)-1.0_real_8
    ENDDO
    ! Initial Guess for GAMY
    !$omp parallel do private(I,J)
    DO i=1,nstate
       DO j=1,nstate
          gamy(i,j)=-0.5_real_8*(s1(i,j)+s1(j,i))
          s4(i,j)=gamy(i,j)
       ENDDO
    ENDDO
    ! Iteration
    iter=0
10  CONTINUE
    iter=iter+1
    CALL dgemm('N','N',nstate,nstate,nstate,1.0_real_8,s2,nstate,gamy,&
         nstate,0.0_real_8,s3,nstate)
    dgmax=0.0_real_8
    !$omp parallel do private(I,J,DGAM) reduction(max:DGMAX)
#ifdef __SR8000
    !poption parallel, tlocal(I,J,DGAM), max(DGMAX)
#endif
    DO i=1,nstate
       DO j=1,nstate
          gamy(i,j)=-0.5_real_8*(s1(i,j)+s1(j,i)+s3(i,j)+s3(j,i))
          dgam = (gamy(i,j)-s4(i,j))*(gamy(i,j)-s4(i,j))
          dgmax=MAX(dgmax,dgam)
          s4(i,j)=gamy(i,j)
       ENDDO
    ENDDO
    dgmax=SQRT(dgmax)
    IF (dgmax.GT.tol) GOTO 10
    ! Apply the constraint
    IF (cntl%tharm) THEN
       !$omp parallel do private(I,J,IG,GIJ,PF4)
#ifdef __SR8000
       !poption parallel, tlocal(I,J,IG,GIJ,PF4)
#endif
       DO i=1,nstate
          DO j=1,nstate
             gij=gamy(i,j)
             DO ig=1,ncpw%ngw
                pf4=dtan2w(ig)/dtb2me
                cm(ig,i)=cm(ig,i)+gij*pf4*c0(ig,j)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(I,J,IG,GIJ,PF4)
#ifdef __SR8000
       !poption parallel, tlocal(I,J,IG,GIJ,PF4)
#endif
       DO i=1,nstate
          DO j=1,nstate
             gij=gamy(i,j)
             DO ig=1,ncpw%ngw
                pf4=cntr%emass/xmu(ig)
                cm(ig,i)=cm(ig,i)+gij*pf4*c0(ig,j)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(s1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s4,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rvharm
  ! ==================================================================
#endif 
  ! ==================================================================
  SUBROUTINE give_scr_rortv(lrortv,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrortv
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: nstx

! ==--------------------------------------------------------------==

    IF (cntl%nonort) THEN
       IF (nort_com%scond.GT.nort_com%slimit) THEN
          CALL give_scr_crotwf(lrortv,tag,nstate)
       ELSE
          lrortv=0
       ENDIF
    ELSE
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          lrortv=3*nstate*nstx
       ELSE
          lrortv=3*nstate*nstate+MAX(nstate*nstate,3*nstate)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rortv
  ! ==================================================================

END MODULE rortv_utils
