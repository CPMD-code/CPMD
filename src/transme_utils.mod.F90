MODULE transme_utils
  USE cdftmod,                         ONLY: cdftlog,&
                                             cdftpi,&
                                             projlog,&
                                             wd,&
                                             wdiff
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cppt,                            ONLY: inyh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filn
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE setirec_utils,                   ONLY: read_irec
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE transmemod,                      ONLY: biga,&
                                             fullp,&
                                             fullw,&
                                             gindex,&
                                             minyh,&
                                             mxhg,&
                                             transme_tol,&
                                             w
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: transme

CONTAINS

  ! ==================================================================
  ! Source file of the cntl%cdft transition matrix element calculation code
  ! H. Oberhofer (ho246@cam.ac.uk) 2009
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE overlap(pa,pb,c,sab)
    ! ==--------------------------------------------------------------==
    ! == Computes the overlap matrix S of PA and PB                   ==
    ! == and its determinant to get SAB and the transposed of its     ==
    ! == Cofactor Matrix                                              ==
    ! == IN: PA,PB coefficients of the Orbitals of States A and B     ==
    ! == OUT: SAB Overlap Matrix Element, C cofactor matrix           ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: pa(ncpw%ngw,*), pb(ncpw%ngw,*)
    REAL(real_8)                             :: c(crge%n,crge%n), sab

    INTEGER                                  :: gl, il, ilaenv, ipiv(crge%n), &
                                                jl, oblock, status
    LOGICAL                                  :: neg
    REAL(real_8)                             :: res, work(crge%n,crge%n)

    sab=1.0_real_8
    CALL zeroing(c)!,n*n)

    IF (cntl%tlsd)THEN
       ! UP
       DO il=1,spin_mod%nsup
          DO jl=1,spin_mod%nsup
             res=0.0_real_8
             !$omp parallel do private(GL) reduction(+:RES)
             DO gl=2,ncpw%ngw
                res=res+2.0_real_8*(REAL(pa(gl,il))*&
                     REAL(pb(gl,jl))+AIMAG(pa(gl,il))*AIMAG(pb(gl,jl)))
             ENDDO
             IF (geq0)THEN
                c(il,jl)=res+REAL(pa(1,il))*REAL(pb(1,jl))
             ELSE
                c(il,jl)=res+2.0_real_8*(REAL(pa(1,il))*REAL(pb(1,jl))+&
                     AIMAG(pa(1,il))*AIMAG(pb(1,jl)))
             ENDIF
          ENDDO
       ENDDO
       ! Down
       DO il=spin_mod%nsup+1,crge%n
          DO jl=spin_mod%nsup+1,crge%n
             res=0.0_real_8
             !$omp parallel do private(GL) reduction(+:RES)
             DO gl=2,ncpw%ngw
                res=res+2.0_real_8*(REAL(pa(gl,il))*&
                     REAL(pb(gl,jl))+AIMAG(pa(gl,il))*AIMAG(pb(gl,jl)))
             ENDDO
             IF (geq0)THEN
                c(il,jl)=res+REAL(pa(1,il))*REAL(pb(1,jl))
             ELSE
                c(il,jl)=res+2.0_real_8*(REAL(pa(1,il))*REAL(pb(1,jl))+&
                     AIMAG(pa(1,il))*AIMAG(pb(1,jl)))
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO il=1,crge%n
          DO jl=1,crge%n
             res=0.0_real_8
             !$omp parallel do private(GL) reduction(+:RES)
             DO gl=2,ncpw%ngw
                res=res+2.0_real_8*(REAL(pa(gl,il))*&
                     REAL(pb(gl,jl))+AIMAG(pa(gl,il))*AIMAG(pb(gl,jl)))
             ENDDO
             IF (geq0)THEN
                c(il,jl)=res+REAL(pa(1,il))*REAL(pb(1,jl))
             ELSE
                c(il,jl)=res+2.0_real_8*(REAL(pa(1,il))*REAL(pb(1,jl))+&
                     AIMAG(pa(1,il))*AIMAG(pb(1,jl)))
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    CALL mp_sum(c,crge%n*crge%n,parai%allgrp)
    IF (paral%io_parent.AND.cdftlog%tphio)THEN
       DO il=1,crge%n
          DO jl=1,crge%n
             WRITE(71,'( F13.8,"  ")',advance="no")c(il,jl)
          ENDDO
          WRITE(71,*) " "
       ENDDO
       WRITE(71,*) " "
    ENDIF

    ! Decompose Matrix and calculate Determinant
    CALL zeroing(ipiv)!,n)
    CALL dgetrf(crge%n,crge%n,c,crge%n,ipiv,status)

    IF (status.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*)"Output of Decomposition ",status
       CALL stopgm('OVERLAP','DECOMPOSITION FAILED',& 
            __LINE__,__FILE__)
    ENDIF
    neg=.FALSE.
    DO il=1,crge%n
       sab=sab*c(il,il)
       IF (ipiv(il).NE.il) neg=(.NOT.neg)
    ENDDO
    IF (neg)sab=sab*(-1.0_real_8)
    IF (.NOT.cntl%tlsd)sab=sab*sab

    ! Calculate the inverse
    oblock = ilaenv(1,"DGETRI"," ",crge%n,crge%n,-1,-1)
    CALL zeroing(work)!,n*n)
    CALL dgetri(crge%n,c,crge%n,ipiv,work,crge%n*oblock,status)
    IF (status.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*)"Output of Inversion ",status
       CALL stopgm('OVERLAP','INVERSION FAILED',& 
            __LINE__,__FILE__)
    ENDIF

    CALL zeroing(work)!,n*n)
    DO il=1,crge%n
       DO jl=1,crge%n
          work(il,jl)=sab*c(jl,il)
       ENDDO
    ENDDO
    CALL dcopy(crge%n*crge%n,work,1,c,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE overlap
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE mkfullp(pa)
    ! ==--------------------------------------------------------------==
    ! == Gathers the full wavefunction                                ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: pa(*)

    INTEGER                                  :: il, msgl

    CALL zeroing(w)!,mxhg)
    DO il=1,ncpw%ngw
       w(il)=pa(il)
    ENDDO

    msgl=2*mxhg*8
    CALL my_concat(w,fullp,msgl,parai%allgrp)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mkfullp
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE transme(pa,pb,ipar,ipbr,fa,fb,va,vb,nai,nbi)
    ! ==--------------------------------------------------------------==
    ! == Computes the transition matrix element                       ==
    ! == IN: PA,PB coefficients of the Orbitals of States A and B     ==
    ! == IN: FA,FB uncorrected energies for States A and B            ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: pa(ncpw%ngw,*), &
                                                pb(ncpw%ngw,*), &
                                                ipar(fpar%nnr1,*), &
                                                ipbr(fpar%nnr1,*)
    REAL(real_8)                             :: fa, fb, va(2), vb(2), nai(2), &
                                                nbi(2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'transme'

    INTEGER                                  :: i, ierr, il, isub, jl, kl, &
                                                lowi, lowi1, lowi2
    REAL(real_8) :: csm, diffa, diffb, e(2), h(2,2), hab, hab2, habs, hba, &
      hba2, hbas, hp(2,2), hproj(2,2), hrest, l(2), na, na2, nas, nb, nb2, &
      nbs, s(2), saa, sbb, t(2,2), up(2,2), v(2,2), vi(2,2), wab(2,2), x(2,2)
    REAL(real_8), ALLOCATABLE                :: caa(:,:), cab(:,:), cba(:,:), &
                                                cbb(:,:), par(:,:), pbr(:,:)

! ,TEST
! for projection

    IF (paral%parent)CALL tiset('   TRANSME',isub)
    hrest=0.0_real_8
    lowi=1
    lowi1=1
    lowi2=spin_mod%nsup+1
    IF (paral%io_parent)WRITE(6,*)"Calculating Overlap Matrix"
    ALLOCATE(cab(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cba(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(caa(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cbb(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(par(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pbr(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent.AND.cdftlog%tphio)THEN
       OPEN(71,file="PHI_MAT",status='UNKNOWN')
       REWIND(unit=71)
    ENDIF
    CALL overlap(pb,pb,cbb,sbb)
    CALL overlap(pa,pa,caa,saa)
    CALL overlap(pa,pb,cab,s(2))
    CALL overlap(pb,pa,cba,s(1))
    IF (paral%io_parent.AND.cdftlog%tphio) CLOSE(71)

    ! ==--------------------------------------------------------------==

    hab=hrest
    hba=hab
    CALL zeroing(par)!,nnr1*n)
    CALL zeroing(pbr)!,nnr1*n)
    IF (paral%io_parent) WRITE(6,*)"Constructing real space wavefunctions"
    kl=part_1d_nbr_el_in_blk(crge%n,parai%cp_inter_me,parai%cp_nogrp)
    kl=part_1d_get_el_in_blk(kl,crge%n,parai%cp_inter_me,parai%cp_nogrp)
    DO i=1,part_1d_nbr_el_in_blk(crge%n,parai%cp_inter_me,parai%cp_nogrp),2
       il=part_1d_get_el_in_blk(i,crge%n,parai%cp_inter_me,parai%cp_nogrp)
       !$omp parallel do private(JL)
       DO jl=1,fpar%nnr1
          par(jl,il)=REAL(ipar(jl,il/2+1))
          IF (il.LT.kl) par(jl,il+1)=AIMAG(ipar(jl,il/2+1))
       ENDDO
    ENDDO
    CALL cp_grp_redist(par,fpar%nnr1,crge%n) ! distribute over all groups
    DO i=1,part_1d_nbr_el_in_blk(crge%n,parai%cp_inter_me,parai%cp_nogrp),2
       il=part_1d_get_el_in_blk(i,crge%n,parai%cp_inter_me,parai%cp_nogrp)
       !$omp parallel do private(JL)
       DO jl=1,fpar%nnr1
          pbr(jl,il)=REAL(ipbr(jl,il/2+1))
          IF (il.LT.kl) pbr(jl,il+1)=AIMAG(ipbr(jl,il/2+1))
       ENDDO
    ENDDO
    CALL cp_grp_redist(pbr,fpar%nnr1,crge%n) ! distribute over all groups

    na=0.0_real_8
    nb=na
    hab=0.0_real_8
    hba=hab
    na2=0.0_real_8
    nb2=na2
    hab2=0.0_real_8
    hba2=hab2
    csm=1.0_real_8
    IF (cdftlog%tspinc)csm=-1.0_real_8

    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) WRITE(6,*)"Starting Summation"
    IF (cntl%tlsd)THEN
       ! Up
       DO il=lowi1,spin_mod%nsup
          DO jl=1,spin_mod%nsup
             hab=hab+matel_real(par(1,il),pbr(1,jl),wd)*cab(il,jl)
             nb=nb+matel_real(pbr(1,il),pbr(1,jl),wd)*cbb(il,jl)
             hba=hba+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na=na+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
       DO il=lowi2,crge%n
          DO jl=spin_mod%nsup+1,crge%n
             hab2=hab2+matel_real(par(1,il),pbr(1,jl),wd)*cab(il,jl)
             nb2=nb2+matel_real(pbr(1,il),pbr(1,jl),wd)*cbb(il,jl)
             hba2=hba2+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na2=na2+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
       IF (.NOT.cdftlog%tcall)THEN
          hab=hab+csm*hab2
          hba=hba+csm*hba2
          na=na+csm*na2
          nb=nb+csm*nb2
       ELSE
          hab=hab+hab2
          hba=hba+hba2
          na=na+na2
          nb=nb+nb2
          habs=hab-2.0_real_8*hab2
          hbas=hba-2.0_real_8*hba2
          nas=na-2.0_real_8*na2
          nbs=nb-2.0_real_8*nb2
       ENDIF
    ELSE
       DO il=lowi,crge%n
          DO jl=1,crge%n
             hab=hab+matel_real(par(1,il),pbr(1,jl),wd)*cab(il,jl)
             nb=nb+matel_real(pbr(1,il),pbr(1,jl),wd)*cbb(il,jl)
             hba=hab+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na=na+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
       hab=hab*2.0_real_8
       nb=nb*2.0_real_8
       hba=hba*2.0_real_8
       na=na*2.0_real_8
    ENDIF
    IF (paral%io_parent)THEN
       diffa=na+zvdiff()
       diffb=nb+zvdiff()
       IF (ABS(diffa-nai(1)).GT.transme_tol&
            .OR.ABS(diffb-nbi(1)).GT.transme_tol)THEN
          WRITE(6,'(1X,64("="),/)')
          WRITE(6,*)"WARNING! "
          WRITE(6,*)"INTEGRATED N VALUES DO NOT COINCIDE WITH N_C"
          WRITE(6,'(/,2X,"Integrated Values")')
          WRITE(6,'(3X,F13.8,F13.8,/)')diffa,diffB
          WRITE(6,'(/,2X,"N_C")')
          WRITE(6,'(3X,F13.8,F13.8,/)')nai(1),nbi(1)
          WRITE(6,'(1X,64("="),/)')
       ENDIF
       IF (cdftlog%tcall)THEN
          IF (ABS(nas-nai(2)).GT.transme_tol&
               .OR.ABS(nbs-nbi(2)).GT.transme_tol)THEN
             WRITE(6,'(1X,64("="),/)')
             WRITE(6,*)"WARNING! "
             WRITE(6,*)"INTEGRATED N VALUES DO NOT COINCIDE WITH N_S"
             WRITE(6,'(/,2X,"Integrated Values")')
             WRITE(6,'(3X,F13.8,F13.8,/)')nas,nbS
             WRITE(6,'(/,2X,"N_S")')
             WRITE(6,'(3X,F13.8,F13.8,/)')nai(2),nbi(2)
             WRITE(6,'(1X,64("="),/)')
          ENDIF
       ENDIF
       IF (cdftlog%tcall)THEN
          wab(1,2)=hba+hbas*vb(2)/vb(1)
          wab(2,1)=hab+habs*va(2)/va(1)
          wab(1,1)=nb+nbs*vb(2)/vb(1)
          wab(2,2)=na+nas*va(2)/va(1)
       ELSE
          wab(1,2)=hba
          wab(2,1)=hab
          wab(1,1)=nb
          wab(2,2)=na
       ENDIF

       h(1,1)=fb
       h(2,2)=fa
       IF (cdftlog%tcall)THEN
          h(2,1)=(fb+vb(1)*nb+vb(2)*nbs)*s(2)-vb(1)*wab(1,2)
          h(1,2)=(fa+va(1)*na+va(2)*nas)*s(1)-va(1)*wab(2,1)
       ELSE
          h(2,1)=(fb+vb(1)*nb)*s(2)-vb(1)*wab(1,2)
          h(1,2)=(fa+va(1)*na)*s(1)-va(1)*wab(2,1)
       ENDIF

       CALL geneigv(s,wab,v,l)
       CALL vinv(v,vi)
       CALL hprime(h,hp,s)
       CALL htrans2(t,v,vi,hp,h)
       CALL odiag(hp,x,e)
    ENDIF

    IF (cdftlog%thdaproj)THEN
       CALL dcopy(4,x,1,up,1)
       hproj(1,1)=e(1)
       hproj(2,2)=e(2)
       projlog%projnow=.TRUE.
       CALL hda_project(pa,pb,up,hproj)
    ENDIF

    IF (paral%io_parent)THEN
       WRITE(6,'(1X,64("-"),/)')
       WRITE(6,'(2X,"Matrix Element Calculation",/)')
       WRITE(6,'(2X,"Overlap Matrix (S)  ChkSum = ",F14.8)')&
            sbb+s(1)+s(2)+saA
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')sbb,s(1),&
            s(2),saA
       WRITE(6,'(2X,"Weight Matrix (W)  ChkSum = ",F14.8)')&
            SUM(wab)
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')wab(1,1),&
            wab(1,2),wab(2,1),wab(2,2)
       WRITE(6,'(2X,A,F14.8)')"Non-orthogonal diabatic H Matrix (H'')"&
            //'  ChkSum = ',SUM(t)
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')t(1,1),t(1,2),&
            t(2,1),t(2,2)
       WRITE(6,'(2X,A,F14.8)')"Full non-orthogonal Hamiltonian (H')"&
            //'  ChkSum = ',SUM(hp)
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')hp(1,1),hp(1,2),&
            hp(2,1),hp(2,2)
       WRITE(6,'(2X,"Generalised Eigenvalues")')
       WRITE(6,'(3X,F14.8,F14.8,/)')l(1),l(2)
       WRITE(6,'(2X,"Generalised Eigenvectors (V)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')v(1,1),v(1,2),&
            v(2,1),v(2,2)
       WRITE(6,'(2X,"Orthogonal diabatic Hamiltonian (H)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')h(1,1),h(1,2),&
            h(2,1),h(2,2)
       WRITE(6,'(2X,"Adiabatic Hamiltonian")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')e(1),0.0_real_8,&
            0.0_real_8,e(2)
       WRITE(6,'(2X,"Mixing Coefficients")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')x(1,1),x(1,2),&
            x(2,1),x(2,2)
       IF (cdftlog%thdaproj)THEN
          WRITE(6,*)
          WRITE(6,'(2X,"Diabatic Hamiltonian in orthogonalised"'//&
               ',/" dressed basis")')
          WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')hproj(1,1),&
               hproj(1,2),hproj(2,1),hproj(2,2)
          WRITE(6,'(2X,"Projection Matrix")')
          WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')&
               up(1,1),up(1,2),up(2,1),up(2,2)
       ENDIF
       WRITE(6,'(1X,64("-"),/)')
       CALL tihalt('   TRANSME',isub)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE transme
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE transmeext(pa,pb,ipar,ipbr,fa,fb)
    ! ==--------------------------------------------------------------==
    ! == Computes the transition matrix element for a given external  ==
    ! ==  potential, potentials for both states are stored in WA/WD   ==
    ! == IN: PA,PB coefficients of the Orbitals of States A and B     ==
    ! == IN: FA,FB uncorrected energies for States A and B            ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: pa(ncpw%ngw,*), &
                                                pb(ncpw%ngw,*), &
                                                ipar(fpar%nnr1,*), &
                                                ipbr(fpar%nnr1,*)
    REAL(real_8)                             :: fa, fb

    CHARACTER(*), PARAMETER                  :: procedureN = 'transmeext'

    INTEGER                                  :: ierr, il, isub, jl, kl, lowi, &
                                                lowi1, lowi2
    REAL(real_8)                             :: h(2,2), hab, hba, hp(2,2), &
                                                hrest, l(2), na, nb, s(2), &
                                                saa, sbb, t(2,2), v(2,2), &
                                                vi(2,2), wab(2,2)
    REAL(real_8), ALLOCATABLE                :: caa(:,:), cab(:,:), cba(:,:), &
                                                cbb(:,:), par(:,:), pbr(:,:)

! ,TEST

    IF (paral%parent)CALL tiset('   TRANSME',isub)
    hrest=0.0_real_8
    lowi=1
    lowi1=1
    lowi2=spin_mod%nsup+1

    IF (paral%io_parent)WRITE(6,*)"Calculating Overlap Matrix"
    ALLOCATE(cab(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cba(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(caa(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cbb(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(par(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pbr(fpar%nnr1,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent.AND.paral%io_parent.AND.cdftlog%tphio)THEN
       OPEN(71,file="PHI_MAT",status='UNKNOWN')
       REWIND(unit=71)
    ENDIF
    CALL overlap(pb,pb,cbb,sbb)
    CALL overlap(pa,pa,caa,saa)
    CALL overlap(pa,pb,cab,s(2))
    CALL overlap(pb,pa,cba,s(1))
    IF ((paral%parent.AND.cdftlog%tphio).AND.paral%io_parent)&
         CLOSE(71)

    ! ==--------------------------------------------------------------==

    hab=hrest
    hba=hab
    CALL zeroing(par)!,nnr1*n)
    CALL zeroing(pbr)!,nnr1*n)
    IF (paral%io_parent)&
         WRITE(6,*)"Constructing real space wavefunctions"
    kl=1
    DO il=1,crge%n,2
       !$omp parallel do private(JL)
       DO jl=1,fpar%nnr1
          par(jl,il)=REAL(ipar(jl,kl))
          IF (il.LT.crge%n) par(jl,il+1)=AIMAG(ipar(jl,kl))
       ENDDO
       kl=kl+1
    ENDDO
    kl=1
    DO il=1,crge%n,2
       !$omp parallel do private(JL)
       DO jl=1,fpar%nnr1
          pbr(jl,il)=REAL(ipbr(jl,kl))
          IF (il.LT.crge%n) pbr(jl,il+1)=AIMAG(ipbr(jl,kl))
       ENDDO
       kl=kl+1
    ENDDO
    na=0.0_real_8
    nb=na
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*)"Starting Summation"
    IF (cntl%tlsd)THEN
       ! Up
       DO il=lowi1,spin_mod%nsup
          DO jl=1,spin_mod%nsup
             hab=hab+matel_real(par(1,il),pbr(1,jl),wdiff)*cab(il,jl)
             nb=nb+matel_real(pbr(1,il),pbr(1,jl),wdiff)*cbb(il,jl)
             hba=hba+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na=na+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
       ! Down
       DO il=lowi2,crge%n
          DO jl=spin_mod%nsup+1,crge%n
             hab=hab+matel_real(par(1,il),pbr(1,jl),wdiff)*cab(il,jl)
             nb=nb+matel_real(pbr(1,il),pbr(1,jl),wdiff)*cbb(il,jl)
             hba=hba+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na=na+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
    ELSE
       DO il=lowi,crge%n
          DO jl=1,crge%n
             hab=hab+matel_real(par(1,il),pbr(1,jl),wdiff)*cab(il,jl)
             nb=nb+matel_real(pbr(1,il),pbr(1,jl),wdiff)*cbb(il,jl)
             hba=hab+matel_real(pbr(1,il),par(1,jl),wdiff)*cba(il,jl)
             na=na+matel_real(par(1,il),par(1,jl),wdiff)*caa(il,jl)
          ENDDO
       ENDDO
       hab=hab*2.0_real_8
       nb=nb*2.0_real_8
       hba=hba*2.0_real_8
       na=na*2.0_real_8
    ENDIF
    IF (paral%parent)THEN
       wab(1,2)=hba
       wab(2,1)=hab
       wab(1,1)=nb
       wab(2,2)=na

       h(1,1)=fb-nb
       h(2,2)=fa-na
       h(2,1)=fb*s(2)-wab(1,2)
       h(1,2)=fa*s(1)-wab(2,1)
       h(2,1)=(h(1,2)+h(2,1))/2.0_real_8
       h(1,2)=h(2,1)

       CALL geneigv(s,wab,v,l)
       CALL vinv(v,vi)
       CALL hprime(h,hp,s)
       CALL htrans2(t,v,vi,hp,h)
    ENDIF

    IF (paral%io_parent)THEN
       WRITE(6,'(1X,64("-"),/)')
       WRITE(6,'(2X,"Matrix Element Calculation",/)')
       WRITE(6,'(2X,"Overlap Matrix (S)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')sbb,s(1),&
            s(2),saA
       WRITE(6,'(2X,"Weight Matrix (W)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')wab(1,1),&
            wab(1,2),wab(2,1),wab(2,2)
       WRITE(6,'(2X,A)')"Non-orthogonal diabatic H Matrix (H'')"
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')t(1,1),t(1,2),&
            t(2,1),t(2,2)
       WRITE(6,'(2X,A)')"Full non-orthogonal Hamiltonian (H')"
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')hp(1,1),hp(1,2),&
            hp(2,1),hp(2,2)
       WRITE(6,'(2X,"Generalised Eigenvalues")')
       WRITE(6,'(3X,F14.8,F14.8,/)')l(1),l(2)
       WRITE(6,'(2X,"Generalised Eigenvectors (V)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')v(1,1),v(1,2),&
            v(2,1),v(2,2)
       WRITE(6,'(2X,"Orthogonal diabatic H Matrix (H)")')
       WRITE(6,'(3X,F14.8,F14.8,/,3X,F14.8,F14.8,/)')h(1,1),h(1,2),&
            h(2,1),h(2,2)
       WRITE(6,'(1X,64("-"),/)')

       CALL tihalt('   TRANSME',isub)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE transmeext
  ! ==================================================================


  ! ==================================================================
  FUNCTION matel(pa)
    ! ==--------------------------------------------------------------==
    ! ==Calculate the W Matrix Element in G space (currently not used)==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: pa(*)
    REAL(real_8)                             :: matel

    COMPLEX(real_8)                          :: buf, pa0con, pb0, res
    INTEGER                                  :: g, gminus(3), gp, gplus(3), &
                                                gpx(3), gs, gstart, gtm, gtp, &
                                                gx(3), ip, nh(3)
    LOGICAL                                  :: con

! GX : G vector; GPX : GP vector; GMINUS :  GPX-GX ; GPLUS : GPX+GX
! GP+G and GP-G
! ==--------------------------------------------------------------==

    nh(1)=spar%nr1s/2+1
    nh(2)=spar%nr2s/2+1
    nh(3)=spar%nr3s/2+1

    pa0con=CONJG(pa(1))
    pb0=fullp(1+parai%igeq0*mxhg)

    IF (geq0)THEN
       res=pa0con*pb0*fullw(1)
       gstart=2
    ELSE
       res=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
       gstart=1
    ENDIF


    DO g=gstart,ncpw%ngw
       res=res+2.0_real_8*pb0*REAL(pa(g)*fullw(g+parai%mepos*mxhg))

       gx(1)=inyh(1,g)-nh(1)
       gx(2)=inyh(2,g)-nh(2)
       gx(3)=inyh(3,g)-nh(3)
       DO ip=0,parai%nproc-1
          IF (ip.EQ.parai%igeq0)THEN
             gs=2
          ELSE
             gs=1
          ENDIF
          DO gp=gs+ip*mxhg,parap%sparm(3,ip)+ip*mxhg
             gpx(1)=minyh(1,gp)-nh(1)
             gpx(2)=minyh(2,gp)-nh(2)
             gpx(3)=minyh(3,gp)-nh(3)
             gminus(1)=gpx(1)-gx(1)
             IF (gminus(1).LT.0)THEN
                gminus(1)=nh(1)-gminus(1)
                con=.TRUE.
             ELSE
                gminus(1)=gminus(1)+nh(1)
                con=.FALSE.
             ENDIF
             gminus(2)=gpx(2)-gx(2)+nh(2)
             gminus(3)=gpx(3)-gx(3)+nh(3)
             gplus(1)=gpx(1)+gx(1)+nh(1)
             gplus(2)=gpx(2)+gx(2)+nh(2)
             gplus(3)=gpx(3)+gx(3)+nh(3)

             gtp=gindex(gplus(1),gplus(2),gplus(3))
             gtm=gindex(gminus(1),gminus(2),gminus(3))

             buf=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             IF (gtp.GT.0)THEN
                buf=buf+REAL(pa(g)*fullp(gp)*fullw(gtp))
             ENDIF
             IF (gtm.GT.0)THEN
                IF (con)THEN
                   buf=buf+REAL(pa(g)*CONJG(fullp(gp))*fullw(gtm))
                ELSE
                   buf=buf+REAL(pa(g)*CONJG(fullp(gp)*fullw(gtm)))
                ENDIF
             ENDIF
             res=res+2.0_real_8*buf
          ENDDO
       ENDDO
    ENDDO
    IF (geq0)THEN
       DO ip=0,parai%nproc-1
          IF (ip.EQ.parai%igeq0)THEN
             gs=2
          ELSE
             gs=1
          ENDIF
          DO gp=gs+ip*mxhg,parap%sparm(3,ip)+ip*mxhg
             res=res+2.0_real_8*pa0con*&
                  REAL(fullp(gp)*fullw(gp))
          ENDDO
       ENDDO
    ENDIF

    matel=REAL(res)
    CALL mp_sum(matel,parai%allgrp)
    RETURN
  END FUNCTION matel
  ! ==================================================================

  ! ==================================================================
  FUNCTION matel_real(pa,pb,w)
    ! ==--------------------------------------------------------------==
    ! == Calculate the W Matrix Element in Real Space (more efficient)==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: pa(*), pb(*), w(*), matel_real

    INTEGER                                  :: il
    REAL(real_8)                             :: res

    res=0.0_real_8
    !$omp parallel do private(IL) reduction(+:RES)
    DO il=1,fpar%nnr1
       res=res+pa(il)*pb(il)*w(il)
    ENDDO

    res=res/(spar%nr1s*spar%nr2s*spar%nr3s)

    CALL mp_sum(res,parai%allgrp)
    matel_real=res
    RETURN
  END FUNCTION matel_real
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE mkgindex()
    ! ==--------------------------------------------------------------==
    ! == Construct the G vector index, belongs to G space MATEL       ==
    ! ==--------------------------------------------------------------==

    CHARACTER(*), PARAMETER                  :: procedureN = 'mkgindex'

    INTEGER                                  :: check, gi, i, i2, ierr, ip, &
                                                ir, j, j2, k, k2, maxl, msgl
    INTEGER, ALLOCATABLE                     :: nibuf(:,:)
    REAL(real_8)                             :: g2, t

    ALLOCATE(gindex(fpar%kr1s,fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    mxhg=parap%sparm(1,0)
    DO ip=1,parai%nproc-1
       IF (parap%sparm(1,ip).GT.mxhg)mxhg=parap%sparm(1,ip)
    ENDDO
    mxhg=mxhg+1
    biga=mxhg*parai%nproc
    ALLOCATE(minyh(3,spar%nhgs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nibuf(3,mxhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Get G Vector Mapping from nodes
    CALL zeroing(nibuf)!,SIZE(nibuf))
    DO gi=1,ncpw%nhg
       nibuf(1,gi)=inyh(1,gi)
       nibuf(2,gi)=inyh(2,gi)
       nibuf(3,gi)=inyh(3,gi)
    ENDDO
    msgl=3*mxhg*8/2
    CALL my_concat(nibuf,minyh,msgl,parai%allgrp)
    ! ==--------------------------------------------------------------==
    maxl=fpar%kr1s*fpar%kr2s*fpar%kr3s
    CALL zeroing(gindex)!,SIZE(gindex))
    check=0
    DO gi=1,biga
       i=minyh(1,gi)
       j=minyh(2,gi)
       k=minyh(3,gi)
       i2=i-spar%nr1s/2-1
       j2=j-spar%nr2s/2-1
       k2=k-spar%nr3s/2-1
       g2=0.0_real_8
       DO ir=1,3
          t=REAL(i2,kind=real_8)*gvec_com%b1(ir)+REAL(j2,kind=real_8)*gvec_com%b2(ir)+REAL(k2,kind=real_8)*gvec_com%b3(ir)
          g2=g2+t*t
       ENDDO
       IF (i.GT.0.AND.j.GT.0.AND.k.GT.0.AND.g2.LT.gvec_com%gcut)THEN
          gindex(i,j,k)=gi
          check=check+1
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE mkgindex
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE geneigv(s,w,v,l)
    ! ==--------------------------------------------------------------==
    ! == Solve the generalised eigenvalue problem                     ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: s(2), w(2,2), v(2,2), l(2)

    REAL(real_8)                             :: a, b, c, d, norm

! (SBA,SAB), ((NB,WBA)(WAB,NA)), ((Eigenvector1)(Eigenvector2))

    a=1-s(1)*s(2)
    b=s(1)*w(2,1)+s(2)*w(1,2)-w(1,1)-w(2,2)
    c=w(1,1)*w(2,2)-w(1,2)*w(2,1)
    d=b*b-4.0_real_8*a*c

    l(1)=(-b-SQRT(d))/(2.0_real_8*a)
    l(2)=(-b+SQRT(d))/(2.0_real_8*a)

    ! V(1,1)=L(1)*S(1)-W(2,1)
    v(1,1)=l(1)*s(1)-w(1,2)
    v(2,1)=w(1,1)-l(1)
    norm=SQRT(v(1,1)*(v(1,1)+s(2)*v(2,1))+&
         v(2,1)*(s(1)*v(1,1)+v(2,1)))
    v(1,1)=v(1,1)/norm! normalise
    v(2,1)=v(2,1)/norm! normalise
    v(1,2)=l(2)*s(1)-w(1,2)
    ! V(1,2)=L(2)*S(1)-W(2,1)
    v(2,2)=w(1,1)-l(2)
    norm=SQRT(v(1,2)*(v(1,2)+s(2)*v(2,2))+&
         v(2,2)*(s(1)*v(1,2)+v(2,2)))
    v(1,2)=v(1,2)/norm! normalise
    v(2,2)=v(2,2)/norm! normalise

    RETURN
  END SUBROUTINE geneigv
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE htrans(t,v,h)
    ! ==--------------------------------------------------------------==
    ! == Transform the Hamiltonian                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t(2,2), v(2,2), h(2,2)

    t(1,1)=h(1,1)
    t(2,1)=h(2,1)! normal
    t(1,2)=h(1,2)! normal
    t(2,2)=h(2,2)

    h(1,1)=v(1,1)*v(1,1)*t(1,1)+v(1,1)*v(2,1)*t(1,2)+&
         v(1,1)*v(2,1)*t(2,1)+v(2,1)*v(2,1)*t(2,2)
    h(2,1)=v(1,1)*v(1,2)*t(1,1)+v(2,1)*v(1,2)*t(1,2)+&
         v(1,1)*v(2,2)*t(2,1)+v(2,1)*v(2,2)*t(2,2)
    h(1,2)=v(1,1)*v(1,2)*t(1,1)+v(1,1)*v(2,2)*t(1,2)+&
         v(2,1)*v(1,2)*t(2,1)+v(2,1)*v(2,2)*t(2,2)
    h(2,2)=v(1,2)*v(1,2)*t(1,1)+v(1,2)*v(2,2)*t(1,2)+&
         v(1,2)*v(2,2)*t(2,1)+v(2,2)*v(2,2)*t(2,2)

    h(1,2)=ABS(h(1,2))
    h(2,1)=ABS(h(2,1))

    RETURN
  END SUBROUTINE htrans
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE htrans2(t,v,vi,hp,h)
    ! ==--------------------------------------------------------------==
    ! == Transform the Hamiltonian                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t(2,2), v(2,2), vi(2,2), &
                                                hp(2,2), h(2,2)

    t(1,1)=h(1,1)
    t(2,1)=h(2,1)
    t(1,2)=h(1,2)
    t(2,2)=h(2,2)

    h(1,1)=v(1,1)*vi(1,1)*hp(1,1)+v(2,1)*vi(1,1)*hp(1,2)+&
         v(1,1)*vi(1,2)*hp(2,1)+v(2,1)*vi(1,2)*hp(2,2)
    h(2,1)=v(1,1)*vi(2,1)*hp(1,1)+v(2,1)*vi(2,1)*hp(1,2)+&
         v(1,1)*vi(2,2)*hp(2,1)+v(2,1)*vi(2,2)*hp(2,2)
    h(1,2)=v(1,2)*vi(1,1)*hp(1,1)+v(2,2)*vi(1,1)*hp(1,2)+&
         v(1,2)*vi(1,2)*hp(2,1)+v(2,2)*vi(1,2)*hp(2,2)
    h(2,2)=v(1,2)*vi(2,1)*hp(1,1)+v(2,2)*vi(2,1)*hp(1,2)+&
         v(1,2)*vi(2,2)*hp(2,1)+v(2,2)*vi(2,2)*hp(2,2)

    h(1,2)=ABS(h(1,2))
    h(2,1)=ABS(h(2,1))

    RETURN
  END SUBROUTINE htrans2
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE htrans3(h,v)
    ! ==--------------------------------------------------------------==
    ! == Transform the Hamiltonian  (V.H.V^T)                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: h(2,2), v(2,2)

    REAL(real_8)                             :: t(2,2)

    t(1,1)=h(1,1)
    t(2,1)=h(2,1)
    t(1,2)=h(1,2)
    t(2,2)=h(2,2)

    h(1,1)=v(1,1)*v(1,1)*t(1,1)+v(1,1)*v(1,2)*t(1,2)+&
         v(1,1)*v(1,2)*t(2,1)+v(1,2)*v(1,2)*t(2,2)
    h(2,1)=v(1,1)*v(2,1)*t(1,1)+v(1,2)*v(2,1)*t(1,2)+&
         v(1,1)*v(2,2)*t(2,1)+v(1,2)*v(2,2)*t(2,2)
    h(1,2)=v(1,1)*v(2,1)*t(1,1)+v(1,1)*v(2,2)*t(1,2)+&
         v(1,2)*v(2,1)*t(1,2)+v(1,2)*v(2,2)*t(2,2)
    h(2,2)=v(2,1)*v(2,1)*t(1,1)+v(2,1)*v(2,2)*t(1,2)+&
         v(2,1)*v(2,2)*t(2,1)+v(2,2)*v(2,2)*t(2,2)

    RETURN
  END SUBROUTINE htrans3
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE odiag(h,v,l)
    ! ==--------------------------------------------------------------==
    ! == Solve the eigenvalue problem                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: h(2,2), v(2,2), l(2)

    REAL(real_8)                             :: t

    t=SQRT(h(1,1)**2+h(2,2)**2-2.0_real_8*h(1,1)*h(2,2)+&
         4.0_real_8*h(1,2)*h(2,1))

    l(1)=(h(1,1)+h(2,2)-t)/2.0_real_8
    l(2)=(h(1,1)+h(2,2)+t)/2.0_real_8

    v(1,1)=-(h(2,2)-h(1,1)+t)/(2.0_real_8*h(2,1))
    v(1,2)=-(h(2,2)-h(1,1)-t)/(2.0_real_8*h(2,1))
    v(2,1)=1.0_real_8
    v(2,2)=1.0_real_8
    t=SQRT(v(1,1)*v(1,1)+v(2,1)*v(2,1))
    v(1,1)=v(1,1)/t
    v(2,1)=v(2,1)/t
    t=SQRT(v(1,2)*v(1,2)+v(2,2)*v(2,2))
    v(1,2)=v(1,2)/t
    v(2,2)=v(2,2)/t

    RETURN
  END SUBROUTINE odiag
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hprime(h,hp,s)
    ! ==--------------------------------------------------------------==
    ! == Construct the full Hamiltonian                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: h(2,2), hp(2,2), s(2)

    REAL(real_8)                             :: t

    t=1.0_real_8-s(1)**2
    h(2,1)=(h(2,1)+h(1,2))/2.0_real_8! symmetrised
    h(1,2)=h(2,1)! symmetrised

    hp(1,1)=(h(1,1)-s(1)*h(2,1))/t
    hp(1,2)=(h(1,2)-s(1)*h(2,2))/t
    hp(2,1)=(h(2,1)-s(2)*h(1,1))/t
    hp(2,2)=(h(2,2)-s(2)*h(1,2))/t

    RETURN
  END SUBROUTINE hprime
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE vinv(v,vi)
    ! ==--------------------------------------------------------------==
    ! == Calculate the inverse transformation matrix                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: v(2,2), vi(2,2)

    REAL(real_8)                             :: t

    t=v(1,1)*v(2,2)-v(1,2)*v(2,1)
    vi(1,1)=v(2,2)/t
    vi(1,2)=-v(1,2)/t
    vi(2,1)=-v(2,1)/t
    vi(2,2)=v(1,1)/t

    RETURN
  END SUBROUTINE vinv
  ! ==================================================================

  ! ==================================================================
  FUNCTION zvdiff()
    ! ==--------------------------------------------------------------==
    ! == Calculate the zv difference of Donor and Acceptor group      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: zvdiff

    INTEGER                                  :: it
    REAL(real_8)                             :: res

    res=0.0_real_8
    DO it=1,cdftpi%naccr
       res=res-ions0%zv(cdftpi%spa(it))
    ENDDO
    DO it=1,cdftpi%ndon
       res=res+ions0%zv(cdftpi%spd(it))
    ENDDO
    zvdiff=res
    RETURN
  END FUNCTION zvdiff

  ! ==================================================================
  SUBROUTINE hda_project(cc1,cc2,u,h)
    ! ==--------------------------------------------------------------==
    ! == Project constrained states on reference states which are     ==
    ! == read from RESTART.REF1 and RESTART.REF2                      ==
    ! == Parameters: CC1,CC2 ... constrained states                   ==
    ! ==             U       ... 2D transform to adiabtic Ham  (in)   ==
    ! ==             U       ... 2D transform to diabtic Ham   (out)  ==
    ! ==             H       ... adiabatic Hamiltonian  (in)          ==
    ! ==             H       ... diabatic Hamiltonian  (out)          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cc1(ncpw%ngw,*), &
                                                cc2(ncpw%ngw,*)
    REAL(real_8)                             :: u(2,2), h(2,2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hda_project'

    CHARACTER(len=20)                        :: filn_back
    COMPLEX(real_8), ALLOCATABLE             :: c01(:,:,:), c02(:,:,:)
    INTEGER                                  :: ierr, irec(100)
    LOGICAL                                  :: rl_save
    REAL(real_8)                             :: ga1, ga2, gb1, gb2, ka1, ka2, &
                                                kb1, kb2, n1, n2, s(2,2), sa, &
                                                sab, sb
    REAL(real_8), ALLOCATABLE                :: cab(:,:)

    ALLOCATE(c01(ncpw%ngw,crge%n,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c02(ncpw%ngw,crge%n,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cab(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    rl_save=restart1%rlate
    filn_back=filn
    restart1%rlate=.FALSE.

    CALL read_irec(irec)
    IF (paral%io_parent)&
         WRITE(6,*)"Reading reference states"
    filn = "RESTART.REF1"
    !vw commeted line (cm not allocated)+ stop
    !vw CALL zhrwf(1,irec,c01,cm,n,eigv,tau0,velp,taui,nfi)
    CALL stopgm(procedureN,'not allocated array',&
         __LINE__,__FILE__)
    filn = "RESTART.REF2"
    !vw commeted line (cm not allocated)+ stop
    !vw CALL zhrwf(1,irec,c02,cm,n,eigv,tau0,velp,taui,nfi)
    CALL stopgm(procedureN,'not allocated array',&
         __LINE__,__FILE__)

    restart1%rlate=rl_save
    filn=filn_back
    IF (paral%io_parent)&
         WRITE(6,*)"Projecting constrained wavefunctions to",&
         " reference states"

    ! Get normalisation for reference states
    CALL overlap(c01,c01,cab,n1)
    CALL overlap(c02,c02,cab,n2)
    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,"Normalisation of reference states")')
       IF (paral%io_parent)&
            WRITE(6,'(3X,F14.8,F14.8,/)')n1,n2
    ENDIF


    CALL overlap(cc2,c01,cab,sb)
    CALL overlap(cc1,c01,cab,sa)
    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,"Overlap with reference state 1")')
       IF (paral%io_parent)&
            WRITE(6,'(3X,F14.8,F14.8,/)')sa,sb
    ENDIF
    ka1=u(1,1)*sb/n1+u(1,2)*sa/n1
    ka2=u(2,1)*sb/n1+u(2,2)*sa/n1
    CALL overlap(cc2,c02,cab,sb)
    CALL overlap(cc1,c02,cab,sa)
    IF (paral%parent)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,"Overlap with reference state 2")')
       IF (paral%io_parent)&
            WRITE(6,'(3X,F14.8,F14.8,/)')sa,sb
    ENDIF
    kb1=u(1,1)*sb/n2+u(1,2)*sa/n2
    kb2=u(2,1)*sb/n2+u(2,2)*sa/n2
    CALL overlap(cc1,cc2,cab,sab)

    n1=SQRT(u(1,1)*u(1,1)+2.0_real_8*u(1,1)*u(1,2)*sab+u(1,2)*u(1,2))
    n2=SQRT(u(2,1)*u(2,1)+2.0_real_8*u(2,1)*u(2,2)*sab+u(2,2)*u(2,2))

    ga1=ka1*u(1,1)/n1+ka2*u(2,1)/n2
    ga2=ka1*u(1,2)/n1+ka2*u(2,2)/n2
    gb1=kb1*u(1,1)/n1+kb2*u(2,1)/n2
    gb2=kb1*u(1,2)/n1+kb2*u(2,2)/n2


    ! Construct S matrix
    s(1,1)=gb1*gb1+gb2*gb2+2.0_real_8*gb1*gb2*sab
    s(2,2)=ga1*ga1+ga2*ga2+2.0_real_8*ga1*ga2*sab
    s(1,2)=ga1*gb1+ga2*gb2+(ga1*gb2+ga2*gb1)*sab
    s(2,1)=s(1,2)

    CALL mat_isqrt(s)

    u(1,1)=s(1,1)*ka1+s(1,2)*kb1
    u(1,2)=s(1,1)*ka2+s(1,2)*kb2
    u(2,1)=s(2,1)*ka1+s(2,2)*kb1
    u(2,2)=s(2,1)*ka2+s(2,2)*kb2

    ! Transform the adiabatic Hamiltonian
    CALL htrans3(h,u)

    RETURN
  END SUBROUTINE hda_project
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE mat_isqrt(u)
    ! ==--------------------------------------------------------------==
    ! == Calculates the square root of the inverse of a symmetric     ==
    ! ==  2x2 Matrix                                                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: u(2,2)

    REAL(real_8)                             :: l(2), v(2,2)

    CALL odiag(u,v,l)

    l(1)=1.0_real_8/SQRT(l(1))
    l(2)=1.0_real_8/SQRT(l(2))

    u(1,1)=v(1,1)*v(1,1)*l(1)+v(1,2)*v(1,2)*l(2)
    u(2,1)=v(1,1)*v(2,1)*l(1)+v(1,2)*v(2,2)*l(2)
    u(1,2)=u(2,1)
    u(2,2)=v(2,1)*v(2,1)*l(1)+v(2,2)*v(2,2)*l(2)

    RETURN
  END SUBROUTINE mat_isqrt
  ! ==================================================================

END MODULE transme_utils
