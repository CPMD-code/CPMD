MODULE eam_pot_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: ry
  USE dpot,                            ONLY: dpot_mod
  USE eam,                             ONLY: eam2,&
                                             team_cross
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE system,                          ONLY: iatpe,&
                                             ipept,&
                                             maxsp,&
                                             maxsys,&
                                             natpe,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eamin
  PUBLIC :: eam_pot

CONTAINS

  ! ==================================================================
  SUBROUTINE eamin
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &EAM &END ON UNIT IUNIT      ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &EAM                                                     ==
    ! ==        Species 1 eamparam                                    ==
    ! ==        Species 2 eamparam                                    ==
    ! ==        ...                                                   ==
    ! ==        Species n eamparam                                    ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  rc re fe phi_e alpha beta gamma rho_e PHI_e E_c             ==
    ! ==  rc re : cntl%bohr                                                ==
    ! ==  fe    : dimensionless                                       ==
    ! ==  phi_e : eV                                                  ==
    ! ==  alpha beta gamma : dimensionless                            ==
    ! ==  rho_e : dimensionless                                       ==
    ! ==  PHI_e : eV                                                  ==
    ! ==  E_c   : eV                                                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i, ierr, is, is2, iunit

! ==--------------------------------------------------------------==
! ==  DEFAULTS                                                    ==
! ==--------------------------------------------------------------==

    IF (.NOT.paral%io_parent) GOTO 9999
    iunit = 5
    CALL zeroing(eam2%eamre )!,maxsp)
    CALL zeroing(eam2%eamfe )!,maxsp)
    CALL zeroing(eam2%eampe )!,maxsp)
    CALL zeroing(eam2%eama  )!,maxsp)
    CALL zeroing(eam2%eamb  )!,maxsp)
    CALL zeroing(eam2%eamc  )!,maxsp)
    CALL zeroing(eam2%eamrho)!,maxsp)
    CALL zeroing(eam2%eambpe)!,maxsp)
    CALL zeroing(eam2%eamec )!,maxsp)
    CALL zeroing(eam2%eamcut)!,maxsp)
    CALL zeroing(eam2%eamsmooth)!,maxsp)
    DO is2=1,ions1%nsp
       DO is=1,ions1%nsp
          team_cross(is,is2) = .FALSE.
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ierr=inscan(iunit,'&EAM')
    IF (ierr.NE.0) GOTO 200
    ! ==--------------------------------------------------------------==
    DO is=1,ions1%nsp
       IF (paral%io_parent)&
            READ(iunit,END=300,err=300,fmt=*)&
            eam2%eamre(is),eam2%eamfe(is),eam2%eampe(is),&
            eam2%eama(is),eam2%eamb(is),eam2%eamc(is),eam2%eamrho(is),&
            eam2%eambpe(is),eam2%eamec(is),eam2%eamcut(is),eam2%eamsmooth(is)
    ENDDO
    DO i=1,ions1%nsp*ions1%nsp
       IF (paral%io_parent)&
            READ(iunit,END=100,err=100,fmt=*) is,is2
       IF (dpot_mod%team(is).AND.dpot_mod%team(is2)) team_cross(is,is2)=.TRUE.
    ENDDO
    ! ==--------------------------------------------------------------==
100 CONTINUE
    DO is2=1,ions1%nsp
       DO is=1,ions1%nsp
          IF ( team_cross(is,is2) ) THEN
             team_cross(is2,is)=.TRUE.
             dpot_mod%team(is)=.TRUE.
             dpot_mod%team(is2)=.TRUE.
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  WRITE INFO TO OUTPUT                                        ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' EMBEDDED-ATOM-MODEL PARAMETER '
    DO is=1,ions1%nsp
       IF (paral%io_parent)&
            WRITE(6,'(A,A3,T20,A,F7.2,A,T48,A,F7.2,A)')&
            " ATOM= ",elem%el(ions0%iatyp(is)),"r_c = ",eam2%eamcut(is)," a.u.",&
            "r_e = ",eam2%eamre(is)," a.u. "
       IF (paral%io_parent)&
            WRITE(6,'(T19,A,F9.5,A)')&
            " r_s = ",eam2%eamsmooth(is)," a.u."
       IF (paral%io_parent)&
            WRITE(6,'(T20,A,F7.2,T48,A,F7.2,A)')&
            "f_e = ",eam2%eamfe(is),"phi_e = ",eam2%eampe(is)," eV"
       IF (paral%io_parent)&
            WRITE(6,'(T20,A,F7.2,T40,A,F7.2,T55,A,F7.2)')&
            "a = ",eam2%eama(is),"b = ",eam2%eamb(is),"c = ",eam2%eamc(is)
       IF (paral%io_parent)&
            WRITE(6,'(T20,A,F7.2,T48,A,F7.2,A)')&
            "rho_e = ",eam2%eamrho(is),"PHI_e = ",eam2%eambpe(is)," eV"
       IF (paral%io_parent)&
            WRITE(6,'(T33,A,F7.2,A)') "Cohesive energy: E_c = ",&
            eam2%eamec(is)," eV"
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
    DO is=1,ions1%nsp
       eam2%eampe(is)=eam2%eampe(is)/(2._real_8*ry)
       eam2%eambpe(is)=eam2%eambpe(is)/(2._real_8*ry)
       eam2%eamec(is)=eam2%eamec(is)/(2._real_8*ry)
    ENDDO
    ! ==--------------------------------------------------------------==
9999 CONTINUE

    CALL mp_bcast_byte(eam2, size_in_bytes_of(eam2),parai%io_source,parai%cp_grp)
    CALL mp_bcast(team_cross,maxsp * ions1%nsp,parai%io_source,parai%cp_grp)
    RETURN

200 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) "EAMIN: No &EAM section"
    CALL stopgm("EAMIN","No &EAM section",& 
         __LINE__,__FILE__)

300 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) "EAMIN: Error reading &EAM section"
    CALL stopgm("EAMIN","Error reading &EAM section",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE eamin
  ! ==================================================================
  SUBROUTINE eam_pot(esr,tau0,iesr,fion,tfor)
    ! nmaxnm:  ..maximum number of nearest neighbours on the list = NMAXNN
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: esr, tau0(:,:,:)
    INTEGER                                  :: iesr
    REAL(real_8)                             :: fion(:,:,:)
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'eam_pot'
    INTEGER, PARAMETER                       :: nmaxnn = 100 
    REAL(real_8), PARAMETER                  :: dens_eps = 1.0e-12_real_8 , &
                                                skin = 2.0_real_8, &
                                                smooth_eps = 1.0e-12_real_8

    INTEGER                                  :: atidx(maxsys%nax,ions1%nsp), &
                                                ia1, ia2, iat, idamax, ierr, &
                                                ii, ioff, is1, is2, isub, &
                                                jat, len
    INTEGER, ALLOCATABLE, SAVE               :: nliste(:,:,:)
    INTEGER, SAVE                            :: iflag = 0
    LOGICAL                                  :: lupdate, tprint
    REAL(real_8) :: b, b1, b2, c1, c2, cutoff, dcut, dcutoff, dd, desr, dfi, &
      dfj, dx, dx1, dy, dy1, dz, dz1, fa, fb, fe, fr12, paa, pbb, r1, r2, re, &
      rho, rr, rrho, rsmooth, x, xp, y, yp
    REAL(real_8), ALLOCATABLE, SAVE          :: dens(:), efio(:,:,:), &
                                                otau(:,:,:), transl(:,:,:)
    REAL(real_8), SAVE                       :: redist

! ==--------------------------------------------------------------==

    IF (iflag.EQ.0) THEN
       ALLOCATE(dens(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(otau(3,maxsys%nax,ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(efio(3,maxsys%nax,ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       natpe=ipept(2,parai%mepos)-ipept(1,parai%mepos)+1
       len=(nmaxnn*2*natpe+1)/2
       ALLOCATE(nliste(nmaxnn,2,len/(2*nmaxnn)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       len=(3*nmaxnn*natpe)
       ALLOCATE(transl(3,nmaxnn,len/(3*nmaxnn)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       lupdate=.TRUE.
       iat=idamax(ions1%nsp,eam2%eamcut,1)
       jat=idamax(ions1%nsp,eam2%eamsmooth,1)
       redist=eam2%eamcut(iat)+skin+eam2%eamsmooth(jat)*LOG(1._real_8/smooth_eps-1._real_8)
       redist=redist**2
       iflag=1
       tprint=.TRUE.
    ELSE
       CALL check_tau(tau0,otau,skin,lupdate)
       tprint=.FALSE.
    ENDIF
    IF (lupdate) CALL do_list(tau0,nliste,redist,iesr,nmaxnn,transl)
    ! ==--------------------------------------------------------------==
    CALL tiset('   EAM_POT',isub)

    iat=0
    DO is1=1,ions1%nsp
       DO ia1=1,ions0%na(is1)
          iat=iat+1
          atidx(ia1,is1) = iat
       ENDDO
    ENDDO

    ! ..calculate the density
    ioff=ipept(1,parai%mepos)-1
    CALL zeroing(dens)!,ions1%nat)
    iat=0
    DO is1=1,ions1%nsp
       IF (.NOT.dpot_mod%team(is1)) go to 1100
       DO ia1=1,ions0%na(is1)
          iat=iat+1
          ii=iat-ioff
          IF (iatpe(iat).NE.parai%mepos) GOTO 1000
          DO jat=1,nmaxnn
             IF (nliste(jat,1,ii).EQ.0) GOTO 1000
             ia2=nliste(jat,1,ii)
             is2=nliste(jat,2,ii)

             IF (.NOT.dpot_mod%team(is2).OR..NOT.team_cross(is1,is2)) go to 1200

             dx=tau0(1,ia1,is1)-tau0(1,ia2,is2)
             dy=tau0(2,ia1,is1)-tau0(2,ia2,is2)
             dz=tau0(3,ia1,is1)-tau0(3,ia2,is2)

             dcut=(eam2%eamcut(is1)+eam2%eamcut(is2))/2._real_8
             rsmooth=(eam2%eamsmooth(is1)+eam2%eamsmooth(is2))/2._real_8

             dx1=dx+transl(1,jat,ii)
             dy1=dy+transl(2,jat,ii)
             dz1=dz+transl(3,jat,ii)
             dd=SQRT(dx1*dx1+dy1*dy1+dz1*dz1)

             cutoff=1._real_8/(EXP((dd-dcut)/rsmooth)+1._real_8)
             IF (cutoff.GT.smooth_eps) THEN
                re=eam2%eamre(is2)
                fe=eam2%eamfe(is2)
                b =eam2%eamb (is2)
                rr=dd/re
                dens(atidx(ia1,is1))=dens(atidx(ia1,is1))&
                     + CUTOFF * FE*EXP(-B*(RR-1._real_8))
                re=eam2%eamre(is1)
                fe=eam2%eamfe(is1)
                b =eam2%eamb (is1)
                rr=dd/re
                dens(atidx(ia2,is2))=dens(atidx(ia2,is2))&
                     + CUTOFF * FE*EXP(-B*(RR-1._real_8))
             ENDIF

1200         CONTINUE
          ENDDO

1000      CONTINUE
       ENDDO

1100   CONTINUE
    ENDDO
    CALL mp_sum(dens,ions1%nat,parai%allgrp)

    IF ( tprint .AND. paral%parent ) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F20.16,/)') "DBG: rho(nat) = ", dens(ions1%nat)
    ENDIF

    ! ..calculate energy and forces
    IF (tfor) CALL zeroing(efio)!,3*maxsys%nax*ions1%nsp)
    iat=0
    DO is1=1,ions1%nsp
       IF (.NOT.dpot_mod%team(is1)) go to 2100

       ! *apsi* TMPTMPTMP
       ! EAMEC(:NSP) = 0.0
       ! EAMBPE(:NSP) = 0.0
       ! EAMPE(:NSP) = 0.0

       DO ia1=1,ions0%na(is1)
          iat=iat+1
          ii=iat-ioff
          IF (iatpe(iat).NE.parai%mepos) GOTO 2000
          rho=dens(atidx(ia1,is1))
          IF (rho.GT.dens_eps) THEN
             rrho=rho/eam2%eamrho(is1)
             x=rrho**(eam2%eama(is1)/eam2%eamb(is1))
             y=rrho**(eam2%eamc(is1)/eam2%eamb(is1))
             esr=esr-eam2%eamec(is1)*(1._real_8-LOG(x))*x-eam2%eambpe(is1)*y
             xp=eam2%eama(is1)/eam2%eamb(is1) * x/rho
             yp=eam2%eamc(is1)/eam2%eamb(is1) * y/rho
             dfi=eam2%eamec(is1)*LOG(x)*xp-eam2%eambpe(is1)*yp
          ELSE
             dfi=0._real_8
          ENDIF

          DO jat=1,nmaxnn
             IF (nliste(jat,1,ii).EQ.0) GOTO 2000
             ia2=nliste(jat,1,ii)
             is2=nliste(jat,2,ii)
             IF (.NOT.dpot_mod%team(is2).OR..NOT.team_cross(is1,is2)) go to 2200

             rho=dens(atidx(ia2,is2))
             IF (rho.GT.dens_eps) THEN
                rrho=rho/eam2%eamrho(is2)
                x=rrho**(eam2%eama(is2)/eam2%eamb(is2))
                y=rrho**(eam2%eamc(is2)/eam2%eamb(is2))

                xp=eam2%eama(is2)/eam2%eamb(is2) * x/rho
                yp=eam2%eamc(is2)/eam2%eamb(is2) * y/rho
                dfj=eam2%eamec(is2)*LOG(x)*xp-eam2%eambpe(is2)*yp
             ELSE
                dfj=0._real_8
             ENDIF

             dx=tau0(1,ia1,is1)-tau0(1,ia2,is2)
             dy=tau0(2,ia1,is1)-tau0(2,ia2,is2)
             dz=tau0(3,ia1,is1)-tau0(3,ia2,is2)

             dx1=dx+transl(1,jat,ii)
             dy1=dy+transl(2,jat,ii)
             dz1=dz+transl(3,jat,ii)
             dd=SQRT(dx1*dx1+dy1*dy1+dz1*dz1)

             fr12=0._real_8
             desr=0._real_8
             dcut=(eam2%eamcut(is1)+eam2%eamcut(is2))/2._real_8
             rsmooth=(eam2%eamsmooth(is1)+eam2%eamsmooth(is2))/2._real_8
             cutoff=1._real_8/(EXP((dd-dcut)/rsmooth)+1._real_8)
             dcutoff=-EXP((dd-dcut)/rsmooth)/rsmooth*cutoff**2

             IF (cutoff.GT.smooth_eps) THEN
                IF (is1.EQ.is2) THEN
                   rr=dd/eam2%eamre(is1)
                   fa=eam2%eampe(is1)*EXP(-eam2%eamc(is1)*(rr-1._real_8))
                   desr=desr+fa
                   IF (tfor) THEN
                      fr12=fr12-eam2%eamc(is1)/eam2%eamre(is1)*fa
                   ENDIF
                ELSE
                   r1=eam2%eamre(is1)
                   b1=eam2%eamb(is1)
                   c1=eam2%eamc(is1)
                   rr=dd/r1
                   fa=eam2%eamfe(is1)*EXP(-b1*(rr-1._real_8))
                   paa=eam2%eampe(is1)*EXP(-c1*(rr-1._real_8))

                   r2=eam2%eamre(is2)
                   b2=eam2%eamb(is2)
                   c2=eam2%eamc(is2)
                   rr=dd/eam2%eamre(is2)
                   fb=eam2%eamfe(is2)*EXP(-b2*(rr-1._real_8))
                   pbb=eam2%eampe(is2)*EXP(-c2*(rr-1._real_8))

                   desr=desr+0.5_real_8*(fb/fa*paa+fa/fb*pbb)
                   IF (tfor) THEN
                      fr12=fr12+0.5_real_8*&
                           ( (-b2/r2+b1/r1-c1/r1)*fb/fa*paa&
                           + (-b1/r1+b2/r2-c2/r2)*fa/fb*pbb )
                   ENDIF
                ENDIF

                IF (tfor) THEN
                   r1=dd/eam2%eamre(is1)
                   r2=dd/eam2%eamre(is2)
                   fr12=fr12&
                        +dfi*eam2%eamfe(is2)*EXP(-eam2%eamb(is2)*(r2-1._real_8))*&
                        (-eam2%eamb(is2)/eam2%eamre(is2))&
                        +dfj*eam2%eamfe(is1)*EXP(-eam2%eamb(is1)*(r1-1._real_8))*&
                        (-eam2%eamb(is1)/eam2%eamre(is1))
                   fr12=fr12*cutoff&
                        +(DFI*eam2%eamfe(IS2)&
                        *EXP(-eam2%eamb(IS2)*(R2-1._real_8))+DESR*0.5_real_8)*DCUTOFF&
                        +(DFJ*eam2%eamfe(IS1)&
                        *EXP(-eam2%eamb(IS1)*(R1-1._real_8))+DESR*0.5_real_8)*DCUTOFF
                   efio(1,ia1,is1)=efio(1,ia1,is1)-dx1/dd*fr12
                   efio(2,ia1,is1)=efio(2,ia1,is1)-dy1/dd*fr12
                   efio(3,ia1,is1)=efio(3,ia1,is1)-dz1/dd*fr12
                   efio(1,ia2,is2)=efio(1,ia2,is2)+dx1/dd*fr12
                   efio(2,ia2,is2)=efio(2,ia2,is2)+dy1/dd*fr12
                   efio(3,ia2,is2)=efio(3,ia2,is2)+dz1/dd*fr12
                ENDIF

                esr = esr + desr * cutoff
             ENDIF

2200         CONTINUE
          ENDDO

2000      CONTINUE
       ENDDO

2100   CONTINUE
    ENDDO

    IF ( tfor ) THEN
       DO is1=1,ions1%nsp
          DO ia1=1,ions0%na(is1)
             fion(1,ia1,is1)=fion(1,ia1,is1)+efio(1,ia1,is1)
             fion(2,ia1,is1)=fion(2,ia1,is1)+efio(2,ia1,is1)
             fion(3,ia1,is1)=fion(3,ia1,is1)+efio(3,ia1,is1)
          ENDDO
       ENDDO
    ENDIF

    CALL tihalt('   EAM_POT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eam_pot
  ! ==================================================================
  SUBROUTINE check_tau(tau0,otau,skin,lupdate)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), otau(:,:,:), skin
    LOGICAL                                  :: lupdate

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: dcut, dd, dx, dy, dz

! VARIABLES
! ==--------------------------------------------------------------==

    dcut=skin**2
    lupdate=.FALSE.
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          dx=tau0(1,ia,is)-otau(1,ia,is)
          dy=tau0(2,ia,is)-otau(2,ia,is)
          dz=tau0(3,ia,is)-otau(3,ia,is)
          CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=dx*dx+dy*dy+dz*dz
          lupdate=lupdate .OR. (dd.GT.dcut)
       ENDDO
    ENDDO
    IF (lupdate) CALL dcopy(3*maxsys%nax*ions1%nsp,tau0,1,otau,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE check_tau
  ! ==================================================================
  SUBROUTINE do_list(tau0,nliste,redist,iesr,nmaxnn,transl)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), redist
    INTEGER                                  :: iesr, nmaxnn, &
                                                nliste(nmaxnn,2,natpe)
    REAL(real_8)                             :: transl(3,nmaxnn,natpe)

    INTEGER                                  :: ia1, ia2, ias, iat1, ii, &
                                                ioff, is1, is2, ix, iy, iz, &
                                                jat, minx
    REAL(real_8)                             :: dd, dx, dx1, dx2, dy, dy1, &
                                                dy2, dz, dz1, dz2

    natpe=ipept(2,parai%mepos)-ipept(1,parai%mepos)+1
    CALL zeroing(nliste)!,nmaxnn*2*natpe)
    CALL zeroing(transl)!,3*nmaxnn*natpe)
    ioff=ipept(1,parai%mepos)-1
    iat1=0
    DO is1=1,ions1%nsp
       DO ia1=1,ions0%na(is1)
          iat1=iat1+1
          ii=iat1-ioff
          IF (iatpe(iat1).NE.parai%mepos) GOTO 1000
          jat=0
          DO is2=is1,ions1%nsp
             ias=1
             IF (is1.EQ.is2) ias=ia1
             DO ia2=ias,ions0%na(is2)
                dx=tau0(1,ia1,is1)-tau0(1,ia2,is2)
                dy=tau0(2,ia1,is1)-tau0(2,ia2,is2)
                dz=tau0(3,ia1,is1)-tau0(3,ia2,is2)
                CALL pbc(dx,dy,dz,dx1,dy1,dz1,1,parm%apbc,parm%ibrav)

                minx=-iesr
                IF (ia1.EQ.ia2.AND.is1.EQ.is2) minx=0

                DO iz=-iesr,iesr
                   DO iy=minx,iesr
                      DO ix=minx,iesr
                         dx2=dx1+ix*metr_com%ht(1,1)+iy*metr_com%ht(2,1)+iz*metr_com%ht(3,1)
                         dy2=dy1+ix*metr_com%ht(1,2)+iy*metr_com%ht(2,2)+iz*metr_com%ht(3,2)
                         dz2=dz1+ix*metr_com%ht(1,3)+iy*metr_com%ht(2,3)+iz*metr_com%ht(3,3)
                         dd=dx2*dx2+dy2*dy2+dz2*dz2

                         IF (dd.LT.redist.AND.dd.GT.1.0e-10_real_8) THEN
                            jat=jat+1
                            IF (jat.GT.nmaxnn) THEN
                               CALL stopgm('DO_LIST','TOO MANY NEIGHBOURS',& 
                                    __LINE__,__FILE__)
                            ENDIF
                            nliste(jat,1,ii)=ia2
                            nliste(jat,2,ii)=is2
                            transl(1,jat,ii)=dx2-dx
                            transl(2,jat,ii)=dy2-dy
                            transl(3,jat,ii)=dz2-dz
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO
          ENDDO

1000      CONTINUE
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE do_list
  ! ==================================================================

END MODULE eam_pot_utils
