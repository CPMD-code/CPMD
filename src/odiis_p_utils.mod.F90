MODULE odiis_p_utils
  USE dotp_utils,                      ONLY: dotp
  USE ener,                            ONLY: ener_com
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE odiis_utils,                     ONLY: solve
  USE parac,                           ONLY: parai
  USE pcgrad_p_utils,                  ONLY: precon_p
  USE system,                          ONLY: cnti,&
                                             maxdis,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: odiis_p
  !public :: updis_p

CONTAINS

  ! ==================================================================
  SUBROUTINE odiis_p(c0,c2,vpp,z11,nstate,pme,gde,svar2,reinit)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vpp(ncpw%ngw)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: pme(ncpw%ngw*nstate+8,*), &
                                                gde((ncpw%ngw*nstate+8)/4,*), &
                                                svar2
    LOGICAL                                  :: reinit

    INTEGER                                  :: i, isub, j, k, nsize
    INTEGER, SAVE                            :: istate = 0, ndiis, nowv
    REAL(real_8)                             :: bc_array(maxdis+1,maxdis+1), &
                                                de1, vc(maxdis+1)
    REAL(real_8), SAVE                       :: diism(maxdis,maxdis), eold, &
                                                gimax(maxdis), gnorm(maxdis), &
                                                grmax(maxdis)

    CALL tiset('   ODIIS_P',isub)
    de1 = 0._real_8
    IF (istate.NE.0) THEN
       ! Check for an increase in energy
       de1=ener_com%etot-eold
    ENDIF

    IF (reinit.OR.istate.EQ.0 .OR. de1 .GE. 1._real_8) THEN
       ! Reinitialization of cntl%diis procedure
       istate=0
       eold=9999._real_8
       reinit=.FALSE.
       ndiis=0
    ENDIF

    IF (istate.EQ.0) CALL updis_p(ndiis,nowv,nsize,cnti%mdiis,0)
    istate=2
    ! Perform an electronic cntl%diis step
    CALL updis_p(ndiis,nowv,nsize,cnti%mdiis,1)

    ! Update cntl%diis buffers
    CALL dscopy(ncpw%ngw*nstate,c0,pme(1,nowv))
    CALL dscal(2*ncpw%ngw*nstate,-1.0_real_8,c2(1,1),1)
    CALL grimax(c2,ncpw%ngw,nstate,grmax(nowv),gimax(nowv),gnorm(nowv))
    CALL dicopy(ncpw%ngw*nstate,c2,gde(1,nowv),&
         grmax(nowv),gimax(nowv))



    CALL precon_p(c2,c2,vpp,nstate,z11,1)
    CALL precon_p(c2,c2,vpp,nstate,z11,1)




    ! Update cntl%diis matrix
    DO i=1,nsize-1
       diism(i,nowv)=0.0_real_8
    ENDDO
    DO i=1,nsize-1
       CALL idcopy(ncpw%ngw*nstate,c0,gde(1,i),grmax(i),gimax(i))
       DO k=1,nstate
          diism(i,nowv)=diism(i,nowv)+dotp(ncpw%ngw,c0(:,k),c2(:,k))
       ENDDO
    ENDDO
    CALL mp_sum(diism(:,nowv),nsize-1,parai%allgrp)

    DO i=1,nsize-1
       diism(nowv,i)=diism(i,nowv)
    ENDDO
    ! Set up cntl%diis Matrix
    CALL zeroing(bc_array)!,(maxdis+1)*(maxdis+1))
    DO i=1,nsize-1
       DO j=1,nsize-1
          bc_array(i,j)=diism(i,j)
       ENDDO
    ENDDO
    DO i=1,nsize-1
       vc(i)=0._real_8
       bc_array(i,nsize)=-1._real_8
       bc_array(nsize,i)=-1._real_8
    ENDDO
    vc(nsize)=-1._real_8
    bc_array(nsize,nsize)=0._real_8

    ! Solve System of Linear Equations
    CALL solve(bc_array,maxdis+1,nsize,vc)
    ! Compute Interpolated Coefficient Vectors
    CALL zeroing(c0)!,SIZE(c0))
    DO i=1,nsize-1
       CALL sdcopy(ncpw%ngw*nstate,pme(1,i),c2)
       CALL daxpy(2*ncpw%ngw*nstate,vc(i),c2(1,1),1,c0(1,1),1)
    ENDDO
    ! Estimate New Parameter Vectors 


    DO i=1,nsize-1
       CALL idcopy(ncpw%ngw*nstate,c2,gde(1,i),grmax(i),gimax(i))
       CALL precon_p(c2,c2,vpp,nstate,z11,1)
       CALL daxpy(2*nstate*ncpw%ngw,-vc(i),c2,1,c0,1)
    ENDDO
    eold=ener_com%etot

    CALL tihalt('   ODIIS_P',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE odiis_p
  ! ==================================================================
  SUBROUTINE updis_p(ndiis,nowv,nsize,mdiis,iact)
    ! ==--------------------------------------------------------------==
    ! == Update variables for cntl%diis if IACT=1                          ==
    ! == Set to 0                  if IACT/=1                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiis, nowv, nsize, mdiis, &
                                                iact

! ==--------------------------------------------------------------==

    IF (iact.EQ.1) THEN
       nowv=MOD(ndiis,mdiis)+1
       ndiis=ndiis+1
       nsize=ndiis+1
       IF (nsize.GT.mdiis) nsize=mdiis+1
    ELSE
       ndiis=0
       nowv=0
       nsize=0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updis_p
  ! ==================================================================

END MODULE odiis_p_utils
