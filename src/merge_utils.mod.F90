MODULE merge_utils
  USE coor,                            ONLY: tau0,&
                                             velp
  USE cppt,                            ONLY: gk
  USE elct,                            ONLY: crge
  USE fft_maxfft,                      ONLY: maxfftn
  USE filnmod,                         ONLY: filn
  USE gsortho_utils,                   ONLY: gsortho
  USE kinds,                           ONLY: real_8
  USE mergemod,                        ONLY: merge01,&
                                             merge02
  USE mp_interface,                    ONLY: mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: merge

CONTAINS

  ! ==================================================================
  SUBROUTINE MERGE(c0,cm,rhoe,psi)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), &
                                                cm(ncpw%ngw,*)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(maxfftn)

    CHARACTER(len=20)                        :: filn_back
    COMPLEX(real_8)                          :: phase
    INTEGER                                  :: i, ig, irec(100), nfi, nstate
    REAL(real_8)                             :: eigv(1), &
                                                taui(3,maxsys%nax,maxsys%nsx)

! variables
! 
! 
! ==================================================================
! nstate=6
! mnst1=10
! mnst2=5

    nstate = merge02%mnst1 + merge02%mnst2

    ! 
    CALL mp_sync(parai%allgrp)
    ! 
    filn_back = filn
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ', ' rco =',restart1%rco
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ', ' rwf =',restart1%rwf
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ',' nstate =',nstate
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ',' n =',crge%n
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ',' mshift()',merge01%mshift(1), merge01%mshift(2),&
            merge01%mshift(3)
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ',' mortho =',merge02%mortho
    ENDIF

    CALL read_irec(irec)
    ! ..   read lower states
    filn = merge01%mfiln1
    CALL zhrwf(1,irec,c0,cm,merge02%mnst1,eigv,tau0,velp,taui,nfi)
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ','filn: ',filn
       IF (paral%io_parent)&
            WRITE (6,*) 'merge1: ',1,' -',merge02%mnst1,' (',merge02%mnst1,')'
    ENDIF
    ! 
    CALL mp_sync(parai%allgrp)
    ! 
    ! ..   read uppers states

    filn = merge01%mfiln2
    CALL zhrwf(1,irec,c0(:,merge02%mnst1+1:,:),cm(1,merge02%mnst1+1),&
         merge02%mnst2,eigv,tau0,velp,taui,nfi)
    ! 
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'merge: ','filn: ',filn
       IF (paral%io_parent)&
            WRITE (6,*) 'merge2: ',merge02%mnst1+1,' -',merge02%mnst1+merge02%mnst2,&
            ' (',merge02%mnst2,')'
    ENDIF
    ! 
    DO i=merge02%mnst1+1,merge02%mnst1+merge02%mnst2
       DO ig=1,ncpw%ngw
          phase=CMPLX(0.0_real_8,-parm%tpiba*(gk(1,ig)*merge01%mshift(1)+&
               gk(2,ig)*merge01%mshift(2)+gk(3,ig)*merge01%mshift(3)),kind=real_8)
          c0(ig,i,1)=c0(ig,i,1)*EXP(phase)
       ENDDO
    ENDDO
    ! 
    IF (merge02%mortho.EQ.1) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'orthogonalizing wavefunctions...'
       ENDIF
       CALL gsortho(c0,ncpw%ngw,1,nstate)
    ENDIF
    CALL rhoofr(c0(:,:,1),rhoe,psi,nstate)
    RETURN
  END SUBROUTINE merge
  ! ==================================================================

END MODULE merge_utils
