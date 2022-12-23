MODULE nl_res_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,&
                                             eigr,&
                                             fnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nl_res
  PUBLIC :: give_scr_nl_res

CONTAINS

  ! ==================================================================
  SUBROUTINE nl_res(is,isa,k,c2,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is, isa, k, nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'nl_res'

    COMPLEX(real_8)                          :: ctm, g
    INTEGER                                  :: i, ierr, ig, ii, isub, iv, &
                                                jv, ki, kj, l, l2, li, lj
    LOGICAL                                  :: debug
    REAL(real_8)                             :: t, tdft, tft
    REAL(real_8), ALLOCATABLE                :: dft(:,:), ft(:,:)

! Variables
! CHARACTER*30 TAG
! dimension EIGR(NGW,NAT)
! dimension TWNL(NGW,maxsys%nhxs,maxsys%nsx)
! dimension FNL2(NAT,maxsys%nhxs,LDF2)
! real(8), allocatable :: FNL(:,:,:)
! real(8), allocatable :: DFNL(:,:,:,:)
! INTEGER LNL_RES
! ==--------------------------------------------------------------==

    IF (ncpw%ngw.EQ.0) RETURN
    CALL tiset('    NL_RES',isub)
    ! 
    debug=.FALSE.
    ! 
    ALLOCATE(ft(maxsys%nhxs,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dft(maxsys%nhxs,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..collect all the FNL and DFNL elements used
    CALL zeroing(ft)!,maxsys%nhxs*nstate)
    CALL zeroing(dft)!,maxsys%nhxs*nstate)
    IF (cntl%tfdist) THEN
       DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
          ii=i-parap%nst12(parai%mepos,1)+1
          DO iv=1,nlps_com%ngh(is)
             ft(iv,i)=fnl2(1,isa,iv,ii,1)
          ENDDO
       ENDDO
       CALL mp_sum(ft,maxsys%nhxs*nstate,parai%allgrp)
    ELSE
       DO i=1,nstate
          DO iv=1,nlps_com%ngh(is)
             ft(iv,i)=fnl(1,isa,iv,i,1)
          ENDDO
       ENDDO
    ENDIF
    DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
       ii=i-parap%nst12(parai%mepos,1)+1
       DO iv=1,nlps_com%ngh(is)
          dft(iv,i)=dfnl(1,isa,iv,k,ii,1)
       ENDDO
    ENDDO
    CALL mp_sum(dft,maxsys%nhxs*nstate,parai%allgrp)
    ! 
    CALL zeroing(c2)!,SIZE(c2))
    ! 
    DO i=1,nstate
       IF (pslo_com%tvan(is)) THEN
          CALL stopgm('FNONLOC_P','VDB not implemented',& 
               __LINE__,__FILE__)
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! ==--------------------------------------------------------------==
          ! ==  Stefan Goedeckers PP                                        ==
          ! ==--------------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             l=nghtol(iv,is)+1
             ki=sgpp2%lfval(iv,is)
             li=sgpp2%lpval(iv,is)
             tft=0._real_8
             tdft=0._real_8
             DO jv=1,nlps_com%ngh(is)
                l2=nghtol(jv,is)+1
                lj=sgpp2%lpval(jv,is)
                IF (l.EQ.l2.AND.li.EQ.lj) THEN
                   kj=sgpp2%lfval(jv,is)
                   tft=tft+ft(jv,i)*sgpp2%hlsg(ki,kj,l,is)
                   tdft=tdft+dft(jv,i)*sgpp2%hlsg(ki,kj,l,is)
                ENDIF
             ENDDO
             ctm=-crge%f(i,1)*(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             DO ig=1,ncpw%ngw
                g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba
                c2(ig,i)=c2(ig,i)+&
                     ctm*eigr(ig,isa,1)*twnl(ig,iv,is,1)*(tdft+g*tft)
             ENDDO
          ENDDO
       ELSE
          ! ==--------------------------------------------------------------==
          ! ==  BACHELET HAMANN SCHLUTER                                    ==
          ! ==--------------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             t=-crge%f(i,1)*wsg(is,iv)
             ctm=t*(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             DO ig=1,ncpw%ngw
                g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba
                c2(ig,i)=c2(ig,i)+ ctm*eigr(ig,isa,1)*twnl(ig,iv,is,1)*&
                     (dft(iv,i)+g*ft(iv,i))
             ENDDO
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
    ENDDO
    DEALLOCATE(ft,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dft,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    NL_RES',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nl_res
  ! ==================================================================
  SUBROUTINE give_scr_nl_res(lnl_res,nstate,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnl_res, nstate
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lnl_res=maxsys%nhxs*nstate*2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_nl_res
  ! ==================================================================


END MODULE nl_res_utils
