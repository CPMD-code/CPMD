MODULE fnonloc_p_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE response_pmod,                   ONLY: dfnl00,&
                                             fnl00
  USE sfac,                            ONLY: eigr
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fnonloc_p

CONTAINS

  ! ==================================================================
  SUBROUTINE fnonloc_p(c0,psi,h1nl,f,nstate,k_,isp_,isa_,&
       ikind)
    ! ==--------------------------------------------------------------==
    ! == Computes H1NL (from FNL00), the perturbation needed for      ==
    ! == the phonon perturbations.                                    ==
    ! == the derivative of the nonlocal pseudopotential with respect  ==
    ! == to the ionic position R(is,isa), applied to c0.              ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: h1nl(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: k_, isp_, isa_, ikind

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnonloc_p'

    COMPLEX(real_8)                          :: ctm, expigr, g, i_to_l
    INTEGER                                  :: ierr, ig, istate, istate_rel, &
                                                isub, iv, jv
    LOGICAL                                  :: match
    REAL(real_8)                             :: Sum_dfnl, Sum_fnl, tw
    REAL(real_8), ALLOCATABLE                :: dloc(:,:)

    CALL tiset(' FNONLOC_P',isub)
    IF (ikind .NE. 1)CALL stopgm('F_NL_P','iKind not supported.',& 
         __LINE__,__FILE__)
    CALL zeroing(h1nl)!,nkpt%ngwk*nstate)
    ! ==--------------------------------------------------------------==
    ! Prepare the DFNL on all processors:
    ALLOCATE(dloc(maxsys%nhxs,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   dloc)!,nstate*maxsys%nhxs)
    DO istate=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
       istate_rel = istate - parap%nst12(parai%mepos,1)+1
       DO iv=1,nlps_com%ngh(isp_)
          dloc(iv,istate)=dfnl00(1,isa_,iv,k_,istate_rel,ikind)
       ENDDO
    ENDDO
    CALL mp_sum(dloc,nstate*maxsys%nhxs,parai%allgrp)

    ! ==--------------------------------------------------------------==
    ! ==  VANDERBILT                                                  ==
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tvan(isp_)) THEN
       CALL stopgm('F_NL_P','Vanderbilt PP not implemented',& 
            __LINE__,__FILE__)
    ELSEIF (sgpp1%tsgp(isp_)) THEN


       ! ==--------------------------------------------------------------==
       ! ==  STEFAN GOEDECKER                                            ==
       ! ==--------------------------------------------------------------==
       ! (ds 3/2000)
       IF (tkpts%tkpnt) CALL stopgm('F_NL_P','K-points are complicated.',& 
            __LINE__,__FILE__)
       IF (cntl%tfdist) CALL stopgm('F_NL_P',&
            'distributed DFNLs not supported yet.',& 
            __LINE__,__FILE__)

       DO istate=1,nstate
          DO iv=1,nlps_com%ngh(isp_)
             Sum_fnl  = 0._real_8
             Sum_dfnl = 0._real_8
             DO jv=1,nlps_com%ngh(isp_)
                match =    (nghtol(iv,isp_) .EQ. nghtol(jv,isp_))&
                     .AND. ( sgpp2%lpval(iv,isp_) .EQ.  sgpp2%lpval(jv,isp_))

                IF (match) THEN
                   Sum_fnl = sum_fnl&
                        + fnl00(1,isa_,jv,istate,ikind)&
                        * sgpp2%hlsg( sgpp2%lfval(iv,isp_), sgpp2%lfval(jv,isp_),&
                        nghtol(iv,isp_)+1, isp_)
                   Sum_dfnl = sum_dfnl&
                        + dloc(jv,istate)&
                        * sgpp2%hlsg( sgpp2%lfval(iv,isp_), sgpp2%lfval(jv,isp_),&
                        nghtol(iv,isp_)+1, isp_)
                ENDIF
             ENDDO         ! jv (projectors)
             ! if (cntl%tfdist) call mp_sum(dd,1,allgrp)!  for later

             ! Second, prepare i^^L, where L is that of the projector "iv":
             i_to_l = - CMPLX(0._real_8,-1._real_8,kind=real_8)**nghtol(iv,isp_)

             ! Finally, put everything together:
             ! H1nl += [ FNL*(-iG) + DFNL ]  * exp(iGR)  * i^L * TWnl * HLSM:
             !$omp parallel do private(ig,expigr,G,tw) shared(i_to_l)
             DO ig=1,ncpw%ngw
                expigr = eigr(ig,isa_,1)
                g      = CMPLX(0._real_8,-parm%tpiba*gk(k_,ig),kind=real_8)
                tw     = f(istate) * TWnl(ig,iv,isp_,ikind)
                h1nl(ig,istate) = h1nl(ig,istate)&
                     + (Sum_fnl * g + Sum_dfnl)&
                     * expigr * i_to_l * tw
             ENDDO         ! ig
          ENDDO             ! iv (projectors)
       ENDDO                 ! istate
    ELSE


       ! ==--------------------------------------------------------------==
       ! ==  BACHELET HAMANN SCHLUTER                                    ==
       ! ==--------------------------------------------------------------==
       DO istate=1,nstate
          DO iv=1,nlps_com%ngh(isp_)
             i_to_l = (0.0_real_8,-1.0_real_8)**nghtol(iv,isp_)
             ctm    = - i_to_l * f(istate)*wsg(isp_,iv)
             DO ig=1,nkpt%ngwk
                g=CMPLX(0._real_8,-parm%tpiba*gk(k_,ig),kind=real_8)

                h1nl(ig,istate)=h1nl(ig,istate)+&
                     TWnl(ig,iv,isp_,ikind) * ctm * eigr(ig,isa_,1)&
                     * (dloc(iv,istate)+ g*fnl00(1,isa_,iv,istate,ikind))
                IF (tkpts%tkpnt.AND.geq0) h1nl(1+ncpw%ngw,istate)=CMPLX(0._real_8,0._real_8,kind=real_8)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(dloc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(' FNONLOC_P',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnonloc_p
  ! ==================================================================

END MODULE fnonloc_p_utils
