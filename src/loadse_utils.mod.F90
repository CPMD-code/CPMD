MODULE loadse_utils
  USE cppt,                            ONLY: hg,&
                                             inyh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: epsg,&
                                             epsgx,&
                                             gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE parac,                           ONLY: parai
  USE sphe,                            ONLY: gcutwmax
  USE system,                          ONLY: &
       iatpe, iatpt, ipept, natpe, ncpw, norbpe, parap, parm, spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: loadse

CONTAINS

  ! ==================================================================
  SUBROUTINE loadse
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'loadse'

    INTEGER                                  :: i, ia, iat, ierr, ig, ir, is, &
                                                j, jmax, jmin, k, kmax, kmin, &
                                                ngwy, nh1, nh2, nh3, nhgy
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: g2, gvsign, t

! ==--------------------------------------------------------------==

    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    parap%nrxpl(0,1)=1
    parap%nrxpl(0,2)=parm%nr1
    ! ALLOCATE MEMORY FOR THE PERMANENT ARRAYS
    ALLOCATE(hg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(inyh(3,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Set up the G-Vectors
    nhgy=ncpw%nhg
    ngwy=ncpw%ngw
    ! To Get unique ordering of G-Vectors, break the symmetry
    ! EPSG=EPSGX * ACOS(-1._real_8) * real(NHGS,kind=real_8)
    ! EPSG=GCUT*EPSGX
    epsg=SQRT(REAL(spar%nhgs-1,kind=real_8))*gvec_com%gcut*epsgx
    ! Construction of g-vectors
    ig=0
    ncpw%ngw=0
    ncpw%nhg=0
    DO i=0,parm%nr1-1
       jmin=-parm%nr2+1
       jmax=parm%nr2-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr3+1
          kmax=parm%nr3-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          !CDIR NOVECTOR
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                g2=g2+t*t
             ENDDO
             IF (g2.LT.gvec_com%gcut) THEN
                ig=ig+1
                ncpw%nhg=ncpw%nhg+1
                IF (g2.LT.gcutwmax) THEN
                   ncpw%ngw=ncpw%ngw+1
                   gvsign=-1._real_8
                ELSE
                   gvsign=+1._real_8
                ENDIF
                ! HG(IG)=G2+SQRT(FLOAT(IG-1))*EPSG*GVSIGN
                ! G2*EPSGX*GVSIGN is the epsilon added (>= epsilon(G2))
                ! SQRT(FLOAT(IG-1)) is to break the symmetry
                hg(ig)=g2*(1._real_8+SQRT(REAL(ig-1,kind=real_8))*epsgx*gvsign)
                inyh(1,ig)=nh1+i
                inyh(2,ig)=nh2+j
                inyh(3,ig)=nh3+k
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    IF (nhgy.NE.ncpw%nhg) THEN
       CALL stopgm('LOADSE','NHG.NE.NHGY',& 
            __LINE__,__FILE__)
    ENDIF
    IF (ngwy.NE.ncpw%ngw) THEN
       CALL stopgm('LOADSE','NGW.NE.NGWY',& 
            __LINE__,__FILE__)
    ENDIF
    ! REORDER THE G S IN ORDER OF INCREASING MAGNITUDE
    CALL gsort(hg,inyh,ncpw%nhg)
    parap%nst12(parai%mepos,1)=1
    parap%nst12(parai%mepos,2)=crge%n
    ! 
    norbpe=crge%n
    natpe=ions1%nat
    ipept(1,parai%mepos)=1
    ipept(2,parai%mepos)=ions1%nat
    ALLOCATE(iatpe(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO j=1,ions1%nat
       iatpe(j)=parai%mepos
    ENDDO

    CALL mm_dim(mm_go_mm,oldstatus)
    ALLOCATE(iatpt(2,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          iatpt(1,iat)=ia
          iatpt(2,iat)=is
       ENDDO
    ENDDO
    CALL mm_dim(mm_revert,oldstatus)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE loadse
  ! ==================================================================

END MODULE loadse_utils
