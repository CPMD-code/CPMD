MODULE printave_utils
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: acc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pacca
  PUBLIC :: paccb
  PUBLIC :: paccc
  PUBLIC :: paccd
  PUBLIC :: pacce

CONTAINS

  ! ==================================================================
  SUBROUTINE pacca(ekinc,tempp,etot,econs,eham,enose,enosp,ecnstr,&
       erestr,disa,tcpu,nfi,it)
    ! ==--------------------------------------------------------------==
    ! == IT:                                                          ==
    ! ==   1  Update accumulators                                     ==
    ! ==      The arguments are used                                  ==
    ! == !=1  Print  averaged values of accumulators                  ==
    ! ==      The arguments are unused                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinc, tempp, etot, econs, &
                                                eham, enose, enosp, ecnstr, &
                                                erestr, disa, tcpu
    INTEGER                                  :: nfi, it

    REAL(real_8)                             :: ac1, ac10, ac2, ac3, ac4, &
                                                ac5, ac6, ac7, ac8, ac9, anor

    IF (it.EQ.1) THEN
       acc(1)=acc(1)+ekinc
       acc(2)=acc(2)+tempp
       acc(3)=acc(3)+etot
       acc(4)=acc(4)+econs
       acc(5)=acc(5)+eham
       acc(6)=acc(6)+enose
       acc(7)=acc(7)+enosp
       acc(8)=acc(8)+ecnstr
       acc(9)=acc(9)+disa
       acc(10)=acc(10)+tcpu
       acc(11)=acc(11)+ekinc*ekinc
       acc(12)=acc(12)+tempp*tempp
       acc(13)=acc(13)+etot*etot
       acc(14)=acc(14)+econs*econs
       acc(15)=acc(15)+eham*eham
       acc(16)=acc(16)+enose*enose
       acc(17)=acc(17)+enosp*enosp
       acc(18)=acc(18)+ecnstr*ecnstr
       acc(19)=acc(19)+disa*disa
       acc(20)=acc(20)+tcpu*tcpu
       acc(21)=acc(21)+erestr
       acc(22)=acc(22)+erestr*erestr
    ELSE
       anor=1._real_8/REAL(nfi,kind=real_8)
       ac1=ABS(acc(11)*anor-(acc(1)*anor)**2)
       ac2=ABS(acc(12)*anor-(acc(2)*anor)**2)
       ac3=ABS(acc(13)*anor-(acc(3)*anor)**2)
       ac4=ABS(acc(14)*anor-(acc(4)*anor)**2)
       ac5=ABS(acc(15)*anor-(acc(5)*anor)**2)
       ac6=ABS(acc(16)*anor-(acc(6)*anor)**2)
       ac7=ABS(acc(17)*anor-(acc(7)*anor)**2)
       ac8=ABS(acc(18)*anor-(acc(8)*anor)**2)
       ac9=ABS(acc(19)*anor-(acc(9)*anor)**2)
       ac10=ABS(acc(22)*anor-(acc(21)*anor)**2)
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",21X,A,21X,"*")') ' AVERAGED QUANTITIES'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '      MEAN VALUE       +/-  RMS DEVIATION'
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '             <x>     [<x^2>-<x>^2]**(1/2)'
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' ELECTRON KINETIC ENERGY  ',&
            acc(1)*anor,SQRT(ac1)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' IONIC TEMPERATURE        ',&
            acc(2)*anor,SQRT(ac2)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' DENSITY FUNCTIONAL ENERGY',&
            acc(3)*anor,SQRT(ac3)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CLASSICAL ENERGY         ',&
            acc(4)*anor,SQRT(ac4)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CONSERVED ENERGY         ',&
            acc(5)*anor,SQRT(ac5)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY ELECTRONS    ',&
            acc(6)*anor,SQRT(ac6)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY IONS         ',&
            acc(7)*anor,SQRT(ac7)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CONSTRAINTS ENERGY       ',&
            acc(8)*anor,SQRT(ac8)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' RESTRAINTS ENERGY        ',&
            acc(21)*anor,SQRT(ac10)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ION DISPLACEMENT         ',&
            acc(9)*anor,SQRT(ac9)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' CPU TIME                 ',&
            acc(10)*anor
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pacca
  ! ==================================================================
  SUBROUTINE paccb(ekinc,tempp,etot,econs,eham,disa,tcpu,&
       temph,volume,enose,enosp,nfi,it)
    ! ==--------------------------------------------------------------==
    ! == IT:                                                          ==
    ! ==   1  Update accumulators                                     ==
    ! ==      The arguments are used                                  ==
    ! == !=1  Print  averaged values of accumulators                  ==
    ! ==      The arguments are unused                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinc, tempp, etot, econs, &
                                                eham, disa, tcpu, temph, &
                                                volume, enose, enosp
    INTEGER                                  :: nfi, it

    REAL(real_8)                             :: ac1, ac10, ac11, ac2, ac3, &
                                                ac4, ac5, ac6, ac8, ac9, anor

    IF (it.EQ.1) THEN
       acc(1)=acc(1)+ekinc
       acc(2)=acc(2)+tempp
       acc(3)=acc(3)+etot
       acc(4)=acc(4)+econs
       acc(5)=acc(5)+eham
       acc(6)=acc(6)+disa
       acc(7)=acc(7)+tcpu
       acc(8)=acc(8)+temph
       acc(9)=acc(9)+volume
       acc(10)=acc(10)+enose
       acc(11)=acc(11)+enosp
       acc(12)=acc(12)+ekinc*ekinc
       acc(13)=acc(13)+tempp*tempp
       acc(14)=acc(14)+etot*etot
       acc(15)=acc(15)+econs*econs
       acc(16)=acc(16)+eham*eham
       acc(17)=acc(17)+disa*disa
       acc(18)=acc(18)+tcpu
       acc(19)=acc(19)+temph*temph
       acc(20)=acc(20)+volume*volume
       acc(21)=acc(21)+enose *enose
       acc(22)=acc(22)+enosp *enosp
    ELSE
       anor=1._real_8/REAL(nfi,kind=real_8)
       ac1=ABS(acc(12)*anor-(acc(1)*anor)**2)
       ac2=ABS(acc(13)*anor-(acc(2)*anor)**2)
       ac3=ABS(acc(14)*anor-(acc(3)*anor)**2)
       ac4=ABS(acc(15)*anor-(acc(4)*anor)**2)
       ac5=ABS(acc(16)*anor-(acc(5)*anor)**2)
       ac6=ABS(acc(17)*anor-(acc(6)*anor)**2)
       ac8=ABS(acc(19)*anor-(acc(8)*anor)**2)
       ac9=ABS(acc(20)*anor-(acc(9)*anor)**2)
       ac10=ABS(acc(21)*anor-(acc(10)*anor)**2)
       ac11=ABS(acc(22)*anor-(acc(11)*anor)**2)
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",21X,A,21X,"*")') ' AVERAGED QUANTITIES'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)') '   MEAN VALUE <x>   DEVIATION <x^2>-<x>^2'
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ELECTRON KINETIC ENERGY  ',&
            acc(1)*anor,SQRT(ac1)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,G16.6)') ' IONIC TEMPERATURE        ',&
            acc(2)*anor,SQRT(ac2)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,G16.6)') ' CELL TEMPERATURE         ',&
            acc(8)*anor,SQRT(ac8)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' DENSITY FUNCTIONAL ENERGY',&
            acc(3)*anor,SQRT(ac3)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CLASSICAL ENERGY         ',&
            acc(4)*anor,SQRT(ac4)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CONSERVED ENERGY         ',&
            acc(5)*anor,SQRT(ac5)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY ELECTRONS    ',&
            acc(10)*anor,SQRT(ac10)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY IONS         ',&
            acc(11)*anor,SQRT(ac11)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ION DISPLACEMENT         ',&
            acc(6)*anor,SQRT(ac6)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' CELL VOLUME              ',&
            acc(9)*anor,SQRT(ac9)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' CPU TIME                 ',&
            acc(7)*anor
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE paccb
  ! ==================================================================
  SUBROUTINE paccc(tempp,etot,econs,enose,enosp,ecnstr,erestr,&
       ebogo,disa,tcpu,nfi,it)
    ! ==--------------------------------------------------------------==
    ! == IT:                                                          ==
    ! ==   1  Update accumulators                                     ==
    ! ==      The arguments are used                                  ==
    ! == !=1  Print  averaged values of accumulators                  ==
    ! ==      The arguments are unused                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tempp, etot, econs, enose, &
                                                enosp, ecnstr, erestr, ebogo, &
                                                disa, tcpu
    INTEGER                                  :: nfi, it

    REAL(real_8)                             :: ac1, ac2, ac3, ac4, ac5, ac6, &
                                                ac7, ac8, ac9, anor

    IF (it.EQ.1) THEN
       acc(1)=acc(1)+tempp
       acc(2)=acc(2)+etot
       acc(3)=acc(3)+econs
       acc(4)=acc(4)+enose
       acc(5)=acc(5)+enosp
       acc(6)=acc(6)+ecnstr
       acc(7)=acc(7)+ebogo
       acc(8)=acc(8)+disa
       acc(9)=acc(9)+tcpu
       acc(10)=acc(10)+erestr
       acc(11)=acc(11)+tempp*tempp
       acc(12)=acc(12)+etot*etot
       acc(13)=acc(13)+econs*econs
       acc(14)=acc(14)+enose*enose
       acc(15)=acc(15)+enosp*enosp
       acc(16)=acc(16)+ecnstr*ecnstr
       acc(17)=acc(17)+ebogo*ebogo
       acc(18)=acc(18)+disa*disa
       acc(19)=acc(19)+tcpu
       acc(20)=acc(20)+erestr*erestr
    ELSE
       anor=1._real_8/REAL(nfi,kind=real_8)
       ac1=ABS(acc(11)*anor-(acc(1)*anor)**2)
       ac2=ABS(acc(12)*anor-(acc(2)*anor)**2)
       ac3=ABS(acc(13)*anor-(acc(3)*anor)**2)
       ac4=ABS(acc(14)*anor-(acc(4)*anor)**2)
       ac5=ABS(acc(15)*anor-(acc(5)*anor)**2)
       ac6=ABS(acc(16)*anor-(acc(6)*anor)**2)
       ac7=ABS(acc(17)*anor-(acc(7)*anor)**2)
       ac8=ABS(acc(18)*anor-(acc(8)*anor)**2)
       ac9=ABS(acc(20)*anor-(acc(10)*anor)**2)
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",21X,A,21X,"*")') ' AVERAGED QUANTITIES'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)') '   MEAN VALUE <x>   DEVIATION <x^2>-<x>^2'
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,G16.6)') ' IONIC TEMPERATURE        ',&
            acc(1)*anor,SQRT(ac1)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' DENSITY FUNCTIONAL ENERGY',&
            acc(2)*anor,SQRT(ac2)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CLASSICAL ENERGY         ',&
            acc(3)*anor,SQRT(ac3)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY ELECTRONS    ',&
            acc(4)*anor,SQRT(ac4)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' NOSE ENERGY IONS         ',&
            acc(5)*anor,SQRT(ac5)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' ENERGY OF CONSTRAINTS    ',&
            acc(6)*anor,SQRT(ac6)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' ENERGY OF RESTRAINTS     ',&
            acc(10)*anor,SQRT(ac9)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' BOGOLIUBOV CORRECTION    ',&
            acc(7)*anor,SQRT(ac7)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ION DISPLACEMENT         ',&
            acc(8)*anor,SQRT(ac8)
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' CPU TIME                 ',&
            acc(9)*anor
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE paccc
  ! ==================================================================
  SUBROUTINE paccd(ekinc,etotv,econs,eham,tempp,eharv,enose,enosp,&
       disa,tcpu,qpkinp,qvkinp,ekina,etota,econsa,&
       ehama,tempa,equant,accus,nfi,it)
    ! ==--------------------------------------------------------------==
    ! == IT:                                                          ==
    ! ==   1  Update accumulators                                     ==
    ! ==      The arguments are used                                  ==
    ! == !=1  Print  averaged values of accumulators                  ==
    ! ==      The arguments are unused                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: ekinc, etotv, econs, eham, tempp, eharv, enose, enosp, &
      disa, tcpu, qpkinp, qvkinp, ekina, etota, econsa, ehama, tempa, equant, &
      accus(*)
    INTEGER                                  :: nfi, it

    REAL(real_8)                             :: ac1, ac2, ac3, ac4, ac5, ac6, &
                                                ac7, ac8, ac9, anor

    IF (it.EQ.1) THEN
       accus(1)=accus(1)+ekinc
       accus(2)=accus(2)+etotv
       accus(3)=accus(3)+econs
       accus(4)=accus(4)+eham
       accus(5)=accus(5)+tempp
       accus(6)=accus(6)+eharv
       accus(7)=accus(7)+enose
       accus(8)=accus(8)+enosp
       accus(9)=accus(9)+disa
       accus(10)=accus(10)+tcpu
       accus(11)=accus(11)+qpkinp
       accus(12)=accus(12)+qvkinp
       accus(13)=accus(13)+ekina
       accus(14)=accus(14)+etota
       accus(15)=accus(15)+econsa
       accus(16)=accus(16)+ehama
       accus(17)=accus(17)+tempa
       accus(18)=accus(18)+equant
       accus(21)=accus(21)+qpkinp*qpkinp
       accus(22)=accus(22)+qvkinp*qvkinp
       accus(23)=accus(23)+ekina*ekina
       accus(24)=accus(24)+etota*etota
       accus(25)=accus(25)+econsa*econsa
       accus(26)=accus(26)+ehama*ehama
       accus(27)=accus(27)+tempa*tempa
       accus(28)=accus(28)+equant*equant
       accus(29)=accus(29)+disa*disa
    ELSE
       anor=1._real_8/REAL(nfi,kind=real_8)
       ac1=accus(21)*anor-(accus(11)*anor)**2
       ac2=accus(22)*anor-(accus(12)*anor)**2
       ac3=accus(23)*anor-(accus(13)*anor)**2
       ac4=accus(24)*anor-(accus(14)*anor)**2
       ac5=accus(25)*anor-(accus(15)*anor)**2
       ac6=accus(26)*anor-(accus(16)*anor)**2
       ac7=accus(27)*anor-(accus(17)*anor)**2
       ac8=accus(28)*anor-(accus(18)*anor)**2
       ac9=accus(29)*anor-(accus( 9)*anor)**2
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",21X,A,21X,"*")') ' AVERAGED QUANTITIES'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '      MEAN VALUE       +/-  RMS DEVIATION'
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '             <x>     [<x^2>-<x>^2]**(1/2)'
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ELECTRON KINETIC ENERGY  ',&
            accus(13)*anor,SQRT(ABS(ac3))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,F16.2)') ' IONIC TEMPERATURE        ',&
            accus(17)*anor,SQRT(ABS(ac7))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' DENSITY FUNCTIONAL ENERGY',&
            accus(14)*anor,SQRT(ABS(ac4))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' QUANTUM KINETIC ENERGY(P)',&
            accus(11)*anor,SQRT(ABS(ac1))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' QUANTUM KINETIC ENERGY(V)',&
            accus(12)*anor,SQRT(ABS(ac2))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' TOTAL QUANTUM ENERGY(P)  ',&
            accus(18)*anor,SQRT(ABS(ac8))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CLASSICAL ENERGY         ',&
            accus(15)*anor,SQRT(ABS(ac5))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CONSERVED ENERGY         ',&
            accus(16)*anor,SQRT(ABS(ac6))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ION DISPLACEMENT         ',&
            accus(9)*anor,SQRT(ABS(ac9))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' CPU TIME                 ',&
            accus(10)*anor
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE paccd
  ! ==================================================================
  SUBROUTINE pacce(ekinc,etotv,econs,eham,tempp,eharv,enose,enosp,&
       disa,tcpu,qpkinp,qvkinp,ekina,etota,econsa,&
       ehama,tempa,equant,temph,volume,accus,nfi,it)
    ! ==--------------------------------------------------------------==
    ! == IT:                                                          ==
    ! ==   1  Update accumulators                                     ==
    ! ==      The arguments are used                                  ==
    ! == !=1  Print  averaged values of accumulators                  ==
    ! ==      The arguments are unused                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: ekinc, etotv, econs, eham, tempp, eharv, enose, enosp, &
      disa, tcpu, qpkinp, qvkinp, ekina, etota, econsa, ehama, tempa, equant, &
      temph, volume, accus(*)
    INTEGER                                  :: nfi, it

    REAL(real_8)                             :: ac1, ac2, ac3, ac4, ac5, ac6, &
                                                ac7, ac8, ac9, ac10, ac11, anor

    IF (it.EQ.1) THEN
       accus(1)=accus(1)+ekinc
       accus(2)=accus(2)+etotv
       accus(3)=accus(3)+econs
       accus(4)=accus(4)+eham
       accus(5)=accus(5)+tempp
       accus(6)=accus(6)+eharv
       accus(7)=accus(7)+enose
       accus(8)=accus(8)+enosp
       accus(9)=accus(9)+disa
       accus(10)=accus(10)+tcpu
       accus(11)=accus(11)+qpkinp
       accus(12)=accus(12)+qvkinp
       accus(13)=accus(13)+ekina
       accus(14)=accus(14)+etota
       accus(15)=accus(15)+econsa
       accus(16)=accus(16)+ehama
       accus(17)=accus(17)+tempa
       accus(18)=accus(18)+equant
       accus(19)=accus(19)+temph
       accus(20)=accus(20)+volume
       accus(21)=accus(21)+qpkinp*qpkinp
       accus(22)=accus(22)+qvkinp*qvkinp
       accus(23)=accus(23)+ekina*ekina
       accus(24)=accus(24)+etota*etota
       accus(25)=accus(25)+econsa*econsa
       accus(26)=accus(26)+ehama*ehama
       accus(27)=accus(27)+tempa*tempa
       accus(28)=accus(28)+equant*equant
       accus(29)=accus(29)+disa*disa
       accus(30)=accus(30)+temph*temph
       accus(31)=accus(31)+volume*volume
    ELSE
       anor=1._real_8/REAL(nfi,kind=real_8)
       ac1=accus(21)*anor-(accus(11)*anor)**2
       ac2=accus(22)*anor-(accus(12)*anor)**2
       ac3=accus(23)*anor-(accus(13)*anor)**2
       ac4=accus(24)*anor-(accus(14)*anor)**2
       ac5=accus(25)*anor-(accus(15)*anor)**2
       ac6=accus(26)*anor-(accus(16)*anor)**2
       ac7=accus(27)*anor-(accus(17)*anor)**2
       ac8=accus(28)*anor-(accus(18)*anor)**2
       ac9=accus(29)*anor-(accus( 9)*anor)**2
       ac10=accus(30)*anor-(accus(19)*anor)**2
       ac11=accus(31)*anor-(accus(20)*anor)**2
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(" *",21X,A,21X,"*")') ' AVERAGED QUANTITIES'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '      MEAN VALUE       +/-  RMS DEVIATION'
       IF (paral%io_parent)&
            WRITE(6,'(24X,A,A)')&
            '             <x>     [<x^2>-<x>^2]**(1/2)'
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ELECTRON KINETIC ENERGY  ',&
            accus(13)*anor,SQRT(ABS(ac3))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,F16.2)') ' IONIC TEMPERATURE        ',&
            accus(17)*anor,SQRT(ABS(ac7))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.2,9X,G16.6)') ' CELL TEMPERATURE         ',&
            accus(19)*anor,SQRT(ABS(ac10))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' DENSITY FUNCTIONAL ENERGY',&
            accus(14)*anor,SQRT(ABS(ac4))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' QUANTUM KINETIC ENERGY(P)',&
            accus(11)*anor,SQRT(ABS(ac1))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' QUANTUM KINETIC ENERGY(V)',&
            accus(12)*anor,SQRT(ABS(ac2))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' TOTAL QUANTUM ENERGY(P)  ',&
            accus(18)*anor,SQRT(ABS(ac8))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CLASSICAL ENERGY         ',&
            accus(15)*anor,SQRT(ABS(ac5))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.6,9X,G16.6)') ' CONSERVED ENERGY         ',&
            accus(16)*anor,SQRT(ABS(ac6))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' ION DISPLACEMENT         ',&
            accus(9)*anor,SQRT(ABS(ac9))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,G12.6,9X,G16.6)') ' CELL VOLUME              ',&
            accus(20)*anor,SQRT(ABS(ac11))
       IF (paral%io_parent)&
            WRITE(6,'(A,2X,F12.4,9X,G16.6)') ' CPU TIME                 ',&
            accus(10)*anor
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pacce
  ! ==================================================================

END MODULE printave_utils
