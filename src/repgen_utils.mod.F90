MODULE repgen_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: factem,&
                                             fbohr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE movi,                            ONLY: imtyp
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: &
       flnm, grandparent, maxnp, pimd1, pimd2, pimd3, pmar, supergroup, &
       supersource, trep, trepnm
  USE pinmtrans_utils,                 ONLY: pinmtrans
  USE prng_utils,                      ONLY: repprngg
  USE store_types,                     ONLY: rout1
  USE system,                          ONLY: maxsys,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: repgen

CONTAINS

  ! ==================================================================
  SUBROUTINE repgen(tau0)
    ! ==--------------------------------------------------------------==
    ! ==   creates the NP-1 replica for the coordinates based on      ==
    ! ==   the first replica (read in from standard input):           ==
    ! ==   uses the Levy-Brown procedure from                         ==
    ! ==   Fosdick, Jordan: PR 143 (1966) p.58                        ==
    ! ==--------------------------------------------------------------==
    ! DM1
    ! DM2
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: ia, iat, ip, is, it0, iz0, k
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: dis, dismax, dismin, gauss, &
                                                om2p, omegap, rmean(3), rnp, &
                                                sigma, tau(maxnp+1), var, &
                                                varfp

! Variables
! ==--------------------------------------------------------------==
! Trotter time displacements
! deb 
! write(6,*) ' me np     ',me,np_total
! deb    

    DO ip=1,pimd3%np_total+1
       tau(ip)=(REAL(ip,kind=real_8)-1._real_8)/REAL(pimd3%np_total,kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  if(tread_cent) then read in classical coordinates as        ==
    ! ==  centroids and sample randomly the normal modes              ==
    ! ==--------------------------------------------------------------==
    IF (grandparent) THEN
       rnp=REAL(pimd3%np_total,kind=real_8)
       IF (pimd1%tread_cent) THEN
          omegap = SQRT(rnp)*pimd2%tempb/factem
          om2p = omegap*omegap
          ! deb     
          ! write(6,*) ' me tempb om2p  ',me,tempb,om2p
          ! deb
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   trepnm(k,ia,is,1)=tau0(k,ia,is)
                   DO ip=2,pimd3%np_total
                      gauss = repprngg()
                      sigma = SQRT(pimd2%tempb/(factem*om2p*flnm(ip)*pmar(is)))
                      ! deb     
                      ! write(6,*) ' is, ia, k, ip ',is,ia,k,ip 
                      ! write(6,*) ' me gauss sigma flnm(ip) pmar(is)  ',
                      ! *       me, gauss, sigma, flnm(ip), pmar(is) 
                      ! deb
                      trepnm(k,ia,is,ip) = gauss*sigma
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          ! TRANSFORM FROM NORMAL MODES TO PRIMITIVE COORDINATES
          CALL pinmtrans(trep,trepnm,1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,trep,1,tau0,1)
       ELSE
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   trep(k,ia,is,1)=tau0(k,ia,is)
                   DO ip=1,pimd3%np_total-1
                      gauss = repprngg()
                      trep(k,ia,is,ip+1)=&
                           ( trep(k,ia,is,ip)*(tau(pimd3%np_total+1)-tau(ip+1))&
                           + trep(k,ia,is,1)*(tau(ip+1)-tau(ip)) ) /&
                           (tau(pimd3%np_total+1)-tau(ip))&
                           + gauss*SQRT(&
                           (factem/(pmar(is)*pimd2%tempb))*&
                           (tau(ip+1)-tau(ip))*(tau(pimd3%np_total+1)-tau(ip+1)) /&
                           (tau(pimd3%np_total+1)-tau(ip))&
                           )
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! DM1
       ! final check if all positions are inside cluster
       IF (isos1%tclust) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   DO ip=1,pimd3%np_total
                      IF ( (trep(k,ia,is,ip).GE.parm%alat) .OR.&
                           (trep(k,ia,is,ip).LE.0._real_8) ) THEN
                         IF (paral%io_parent)&
                              WRITE(6,*)&
                              ' ! ! STOP IN REPGEN: TAU(K,IA,IS,IP) OUTSIDE CLUSTER'
                         IF (paral%io_parent)&
                              WRITE(6,*) parm%alat,sngl(trep(k,ia,is,ip)),k,ia,is,ip
                         CALL stopgm('REPGEN','REPLICA OUTSIDE CLUSTER',& 
                              __LINE__,__FILE__)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! DM2
       ! 
       ! check the radius of gyration and compare to free particle
       ! (Gillan in Computer Modelling of Fluids Polymers and Solids 1990) 
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,63("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",61X,"*")')
       IF (paral%io_parent)&
            WRITE(6,'(" *",23X,A,23X,"*")') ' R OF GYRATION '
       IF (paral%io_parent)&
            WRITE(6,'(" *",61X,"*")')
       IF (paral%io_parent)&
            WRITE(6,'(1X,63("*"),/)')
       IF (pimd1%tread_cent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' CLASSICAL CONFIGURATION READ IN AS CENTROIDS'
          IF (paral%io_parent)&
               WRITE(6,*) ' NORMAL MODES SAMPLED DIRECTLY'
          IF (.NOT.pimd1%tpinm) THEN
             IF (pimd1%tstage) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) ' STAGING COORDINATES COMPUTED FROM NORMAL MODES'
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,*)&
                     ' PRIMITIVE COORDINATES COMPUTED FROM NORMAL MODES'
             ENDIF
          ENDIF
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) '  REPLICA HAVE GAUSSIAN DISTRIBUTION AROUND IP=1'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*) '  '
       IF (paral%io_parent)&
            WRITE(6,'(A,A,A,A,A)') '   ATOM   R OF GYRATION',&
            '  FREE PARTICLE',&
            '  MAXDIS',&
            '  MINDIS',&
            ' (ANGSTROM)'
       ! DM     *                 ' (ANGSTROEM)'
       ! 
       iat=0
       DO is=1,ions1%nsp
          varfp=(factem/(4.0_real_8*pmar(is)*pimd2%tempb))
          varfp=SQRT(varfp)
          DO ia=1,ions0%na(is)
             DO k=1,3
                rmean(k)=0.0_real_8
                DO ip=1,pimd3%np_total
                   rmean(k)=rmean(k)+trep(k,ia,is,ip)
                ENDDO
                rmean(k)=rmean(k)/rnp
             ENDDO
             var=0.0_real_8
             DO ip=1,pimd3%np_total
                DO k=1,3
                   var=var+(trep(k,ia,is,ip)-rmean(k))**2
                ENDDO
             ENDDO
             var=SQRT(var/rnp)
             dismax=-10._real_8
             dismin=100000.0_real_8
             DO ip=1,pimd3%np_total
                IF (ip.EQ.pimd3%np_total) THEN
                   dis=(trep(1,ia,is,ip)-trep(1,ia,is,1))**2&
                        +(trep(2,ia,is,ip)-trep(2,ia,is,1))**2&
                        +(trep(3,ia,is,ip)-trep(3,ia,is,1))**2
                ELSE
                   dis=(trep(1,ia,is,ip)-trep(1,ia,is,ip+1))**2&
                        +(trep(2,ia,is,ip)-trep(2,ia,is,ip+1))**2&
                        +(trep(3,ia,is,ip)-trep(3,ia,is,ip+1))**2
                ENDIF
                dis=SQRT(dis)
                IF (dis.GT.dismax) dismax=dis
                IF (dis.LT.dismin) dismin=dis
             ENDDO
             ! 
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(I4,2X,A2,3X,F8.5,7X,F8.5,4X,F8.5,1X,F8.5)')&
                  iat,elem%el(ions0%iatyp(is)),&
                  var/fbohr,varfp/fbohr,dismax/fbohr,dismin/fbohr
          ENDDO
       ENDDO
       ! deb
       ! write(6,*) ' me fbohr ',me,fbohr 
       ! deb
       IF (paral%io_parent)&
            WRITE(6,*) '  '
       ! 
       ! CHECK POSITIONS OF CENTROIDS
       ! 
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                rmean(k)=0.0_real_8
                DO ip=1,pimd3%np_total
                   rmean(k)=rmean(k)+trep(k,ia,is,ip)
                ENDDO
                rmean(k)=rmean(k)/rnp
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    CALL mp_bcast(trep,3*maxsys%nax*maxsys%nsx*pimd3%np_total,supersource,supergroup)
    IF (pimd1%tread_cent) THEN
       CALL mp_bcast(trepnm,3*maxsys%nax*maxsys%nsx*pimd3%np_total,supersource,supergroup)
    ENDIF
    ! DM1
    ! ..WRITE MOVIE FILE
    IF (grandparent.AND.rout1%mout) THEN
       IF (paral%io_parent)&
            CALL fileopen(11,'MOVIE',fo_def,ferror)
       DO ip=1,pimd3%np_total
          DO is=1,ions1%nsp
             iz0=imtyp(is)
             DO ia=1,ions0%na(is)
                IF (ions0%iatyp(is).EQ.4) THEN
                   ! ..Change Be -> Na ; problem with movie
                   it0=11
                   IF (paral%io_parent)&
                        WRITE(11,'(3(2X,F12.4),2I4)')&
                        (trep(k,ia,is,ip)/fbohr,k=1,3),it0,iz0
                ELSE
                   IF (paral%io_parent)&
                        WRITE(11,'(3(2X,F12.4),2I4)')&
                        (trep(k,ia,is,ip)/fbohr,k=1,3),ions0%iatyp(is),iz0
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(11)
    ENDIF
    ! DM2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE repgen
  ! ==================================================================

END MODULE repgen_utils
