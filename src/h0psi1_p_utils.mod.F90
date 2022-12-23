MODULE h0psi1_p_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE kinds,                           ONLY: real_8
  USE perturbation_p_utils,            ONLY: restrain_ngw_zero
  USE response_pmod,                   ONLY: dmbi,&
                                             response1,&
                                             vofrho0
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_p_utils,                    ONLY: v0psi1_p

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: h0psi1_p
  PUBLIC :: give_scr_h0psi1

CONTAINS

  ! ==================================================================
  SUBROUTINE h0psi1_p(c1,c0,drhoe,f,tau0,fion,psi,nstate,&
       c12)
    ! ==================================================================
    ! Arguments:
    COMPLEX(real_8)                          :: c1(nkpt%ngwk,*), &
                                                c0(nkpt%ngwk,*)
    REAL(real_8)                             :: drhoe(*), tau0(:,:,:), &
                                                fion(:,:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: c12(nkpt%ngwk,*)

    INTEGER                                  :: isub
    LOGICAL                                  :: tfor

! Local variables:

    CALL tiset('h0psi1           ',isub)

    ! calculate   fnl(g) = sum(l)      exp igR  *  i^l  proj[l]  c1(g)
    ! and in the case of PHONON perturbation also
    ! dfnl(g) = sum(l) ig * exp igR  *  i^l  proj[l]  c1(g):
    tfor = (response1%phonon .OR. response1%tlanphon .OR. response1%teigens .OR.\
    cntl%tinr .OR. response1%tvoa)
    IF (dmbi%cutoff_restr) CALL restrain_ngw_zero(c1,dmbi%ngw_zero,ncpw%ngw,nstate)
    CALL rnlsm(c1(:,1:nstate),nstate,1,1,tfor)

    ! c12(g)  =  <g|  v_Hxc(0)  |c1>; vpsi also initializes c12
    CALL v0psi1_p(c0,c1,c12,VofRho0,psi,drhoe,nstate)
    ! NB: V0PSI1_P also calculates the first order density response DRHOE

    ! c12(g) +=  <g|  | fnl[c1] >
    CALL fnonloc(c12,f,nstate,1,clsd%nlsd,.TRUE.)

    CALL tihalt('h0psi1          ',isub)
    RETURN
  END SUBROUTINE h0psi1_p
  ! ==================================================================
  SUBROUTINE give_scr_h0psi1(l_h0psi1,tag,nstate)
    INTEGER                                  :: l_h0psi1
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_auxc, l_ddia, l_rnlsm, &
                                                lv0psi1
    LOGICAL                                  :: tfor

    tfor = (response1%phonon .OR. response1%tlanphon .OR. response1%teigens .OR.\
    cntl%tinr .OR. response1%tvoa)

    CALL give_scr_rnlsm(l_rnlsm,tag,nstate,tfor)
    CALL give_scr_fnonloc(l_auxc,l_ddia,nstate)
    lv0psi1 = 2*maxfft

    l_h0psi1 = MAX(l_rnlsm,l_auxc+l_ddia,lv0psi1)
    RETURN
  END SUBROUTINE give_scr_h0psi1
  ! ==================================================================

END MODULE h0psi1_p_utils
