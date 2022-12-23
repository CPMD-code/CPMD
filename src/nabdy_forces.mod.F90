MODULE nabdy_forces
  USE coor,                            ONLY: fion,&
                                             velp
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nabdy_types,                     ONLY: nabdyvar,&
                                             nabdyvec
  USE parac,                           ONLY: paral
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nabdy_forces_calc

CONTAINS

  ! ==================================================================
  SUBROUTINE nabdy_forces_calc(itraj)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: itraj

    CHARACTER(*), PARAMETER :: procedureN = 'nabdy_forces_calc'

    INTEGER                                  :: ia, is, j, k
    LOGICAL                                  :: add_friction, debug
    REAL(real_8)                             :: df(3), dg(3), dist2, f, &
                                                friction, ftot_sfe, g, &
                                                gtot_sfe, length1, length2, &
                                                perc, vect(3)

    debug=.TRUE.
    add_friction=.FALSE.
    friction=1.0e-9_real_8

    CALL modfor(fion,length1,itraj,1)
    IF (add_friction) THEN
       CALL stopgm(procedureN,'please fix me',& 
            __LINE__,__FILE__)
       !vw there is a bug here! call mp_bcast(velp,parai%source,parai%allgrp)
    ENDIF

    DO is=1,ions1%nsp
       IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
       DO ia=1,ions0%na(is)
          !
          IF (paral%io_parent.AND.debug) WRITE(6,'(a,2i5,3f22.7)') 'check-1', &
               ia,is,(fion(k,ia,is),k=1,3)

          f=0.0_real_8
          g=0.0_real_8
          DO k=1,3
             df(k)=0.0_real_8
             dg(k)=0.0_real_8
          ENDDO
          gtot_sfe=0.0_real_8
          ftot_sfe=0.0_real_8 
          DO j=1,nabdyvar%ntrajbd
             !      if (j.eq.itraj) cycle
             CALL dist_traj(is,ia,itraj,j,dist2,vect)

             !      if (paral%io_parent.and.debug) write(6,*) '1',itraj,j,dist2
             !      if (paral%io_parent.and.debug) write(6,*) '2',nabdyvec%naomega(ia,is,j),nabdyvec%naampl(ia,is,j)

             f=f+((1._real_8/(nabdyvec%naomega(ia,is,j))**2) * &
                  (nabdyvec%naampl(ia,is,j)/(nabdyvec%nanormalf(ia,is,j))) &
                  ) * &
                  EXP(-dist2/(2.0_real_8*(nabdyvec%naomega(ia,is,j)**2))) * &
                  (3._real_8-(1._real_8/(nabdyvec%naomega(ia,is,j))**2)*dist2)  
             g=g+(nabdyvec%naampl(ia,is,j))/(nabdyvec%nanormalf(ia,is,j))* &
                  EXP(-dist2/(2.0_real_8*(nabdyvec%naomega(ia,is,j)**2)))
             !       nabdyvec%naqp_pat(ia,is,j)=-0.5_real_8*(1.0/rmass%pma(is))*f/g
             gtot_sfe=gtot_sfe+g
             ftot_sfe=ftot_sfe+f
             DO k=1,3
                df(k)=df(k)+((1._real_8/(nabdyvec%naomega(ia,is,j))**4) * &
                     (nabdyvec%naampl(ia,is,j)/(nabdyvec%nanormalf(ia,is,j))) &
                     ) * &
                     EXP(-dist2/(2.0_real_8*(nabdyvec%naomega(ia,is,j)**2))) * &
                     vect(k) * &
                     ((1._real_8/(nabdyvec%naomega(ia,is,j))**2)*dist2-5._real_8) 
                dg(k)=dg(k)-((1._real_8/(nabdyvec%naomega(ia,is,j))**2) * &
                     (nabdyvec%naampl(ia,is,j)/(nabdyvec%nanormalf(ia,is,j))) &
                     ) * &
                     EXP(-dist2/(2.0_real_8*(nabdyvec%naomega(ia,is,j)**2))) * &
                     vect(k)
             ENDDO
          ENDDO

          nabdyvec%naqp_pat(ia,is)=-0.5_real_8*(1._real_8/rmass%pma(is))*(ftot_sfe/gtot_sfe)

          IF (.NOT.nabdyvar%naforce_screened) THEN
             DO k=1,3
                nabdyvec%nafor(k,ia,is,itraj)= &
                     -0.5_real_8*(1._real_8/rmass%pma(is))*(df(k)-(f/g)*dg(k))/g 
                !    -0.5_real_8*(1._real_8/rmass%pma(is))*(df(k)*g-f*dg(k))/(g**2)
                fion(k,ia,is)=fion(k,ia,is)+nabdyvec%nafor(k,ia,is,itraj)
             ENDDO
          ELSE
             DO k=1,3
                nabdyvec%nafor(k,ia,is,itraj)= &
                     -0.5_real_8*(1._real_8/rmass%pma(is))*(df(k)-(f/g)*dg(k)) &
                     /SQRT(g**2+nabdyvar%nasoftening)
                fion(k,ia,is)=fion(k,ia,is)+nabdyvec%nafor(k,ia,is,itraj)
             ENDDO
          ENDIF
          !
          !     if (paral%io_parent.and.debug) write(6,'(a,5f22.7)') '3',g,f,df(1),dg(1),
          !                                           NAFOR(1,IA,IS,ITRAJ)
          IF (paral%io_parent.AND.debug) WRITE(6,'(a,2i5,3f22.7)') 'check-2', &
               ia,is,(nabdyvec%nafor(k,ia,is,itraj),k=1,3)

          IF (add_friction) THEN
             DO k=1,3
                fion(k,ia,is)=fion(k,ia,is)-friction*rmass%pma(is)*velp(k,ia,is)
             ENDDO
          ENDIF

       ENDDO
    ENDDO
    
    CALL modfor(nabdyvec%nafor,length2,itraj,2)
    IF (length1.NE.0.0_real_8) THEN
     perc=length2/length1*100.0
    ELSE
     perc=1.0_real_8
    ENDIF
    IF (paral%io_parent.AND.debug) WRITE(6,'(a,2x,3e12.5)') ' force ratio fq/fc %: ', &
         length1,length2, perc
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_forces_calc
  ! ==================================================================
  SUBROUTINE dist_traj(is,ia,i,j,dist2,vect)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is, ia, i, j
    REAL(real_8)                             :: dist2, vect(3)

    INTEGER                                  :: k

! variables

    dist2=0.0_real_8
    DO k=1,3
       vect(k)=(nabdyvec%nacoor(k,ia,is,i)-nabdyvec%nacoor(k,ia,is,j))
       dist2=dist2+vect(k)**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_traj
  ! ==================================================================
  SUBROUTINE dist_grid(is,ia,j,rcoord,dist2,vect)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is, ia, j
    REAL(real_8)                             :: rcoord(3), dist2, vect(3)

    INTEGER                                  :: k

! variables

    dist2=0.0_real_8
    DO k=1,3
       vect(k)=(rcoord(k)-nabdyvec%nacoor(k,ia,is,j))
       dist2=dist2+vect(k)**2 
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_grid
  ! ==================================================================
  SUBROUTINE modfor(vect,length,itraj,flag)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: vect(3,maxsys%nax,maxsys%nsx,*), length
    INTEGER                                  :: itraj, flag

    INTEGER                                  :: i, ia, is, k
    REAL(real_8)                             :: temp

    length=0.0_real_8
    IF (flag.EQ.1) THEN
       DO is=1,ions1%nsp
          IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
          DO ia=1,ions0%na(is)
             DO k=1,3
                length=length+(vect(k,ia,is,1))**2
             ENDDO
          ENDDO
       ENDDO
    ELSEIF (flag.EQ.2) THEN
       !DO i=1,nabdyvar%ntrajbd
          temp=0._real_8
          DO is=1,ions1%nsp
             IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
             DO ia=1,ions0%na(is)
                DO k=1,3
                   temp=temp+(vect(k,ia,is,itraj))**2
                ENDDO
             ENDDO
          ENDDO
          length=length+temp
       !ENDDO
       length=length/nabdyvar%ntrajbd
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE modfor
  ! ==================================================================

END MODULE nabdy_forces
