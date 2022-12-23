! GROMOSDIP
C     ==================================================================
       subroutine mm_pardef(PMA,PMA0,PMAT0,PMATOT,
     c          mm_charge,box_au,INC_l,NCAG_l,
     c            NRAM_gr,NCONS_gr,
     c           IATYP,mm_RAGGIO,
     c           NRPT,NSOLV,nax,cpat,cpsp,NSPm,NAm)
C     ==--------------------------------------------------------------==
cmb - Revised on 01 Spetember 2006
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
       use adat
       use cnst
       use coordsz
       implicit none

       integer nax
       REAL*8   PMAT0,PMATOT,PMA(*),PMA0(*),mm_RAGGIO(*)
       REAL*8   mm_charge(nax,*),box_au(3)
       integer  INC_l(*),NCAG_l,NRAM_gr,NCONS_gr,IATYP(*)
       integer icg,cpat(*),cpsp(*),NRPT,NSOLV,NSPm,NAm(*)
c mass and charges of the atoms 
       include 'Gromos/toposz.h'
       include 'Gromos/topoar.h'
       !include 'Gromos/coordsz.h'
       include 'Gromos/box.h'
c
c locals
       integer i,is,sn
C Must be SCMASS !!!
C       REAL*8     mass_prot
C       parameter(mass_prot=1822.89D0)
c gromos charged are multiplied by 1./sqrt(\epsilon_0  4 \pi)=sqrt(138.9354)
c  to get back charges in units of e_0, we have to multiply by g_to_au.
       REAL*8     g_to_au
       parameter(g_to_au=8.483861D-2)
       real*8 nm_to_au
C       parameter(nm_to_au = 10.D0/0.529177248D0)
c
c locals
       integer NSM
       integer zatom
       external zatom
       character*(MAXTLE) name

       nm_to_au=10.D0*FBOHR
c
c gromos masses are assigned to the atoms
c solute
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif 
       DO i=1,NRPT
         PMA0(cpsp(i))=WMAS(i)
       ENDDO
c solvent
!$OMP parallel do private(i,sn)
#ifdef __SR11000
*poption parallel, tlocal(i,sn)
#endif 
       DO i=NRPT+1,NRPT+NSOLV
         sn=mod(i-1-NRPT,NRAM)+1
         PMA0(cpsp(i))=WMASS(sn)
       ENDDO
c
       PMAT0 =0.0d0
!$OMP parallel do private(is) reduction(+:PMAT0)
#ifdef __SR11000
*poption parallel, tlocal(is), psum(PMAT0)
#endif 
       do is=1,NSPm
         PMAT0 = PMAT0 + PMA(is)*NAm(is)
         PMA(is)=PMA0(is)*SCMASS
       enddo
       PMATOT=PMAT0*SCMASS
!
!$OMP parallel do private(i)
#ifdef __SR11000
*poption parallel, tlocal(i)
#endif 
       DO i=1,NRPT        !solute atoms
         mm_charge(cpat(i),cpsp(i))=g_to_au*CG(mod(i-1,NRP)+1)
       ENDDO
!$OMP parallel do private(i,sn)
#ifdef __SR11000
*poption parallel, tlocal(i,sn)
#endif 
       DO i=NRPT+1,NRPT+NSOLV  !sovent charges
         sn=mod(i-1-NRPT,NRAM)+1
         mm_charge(cpat(i),cpsp(i))=g_to_au*CGS(sn)
       ENDDO
!
      do i=1,3
        BOX_au(i)=BOX(i)*nm_to_au
      enddo
! fake charge groups for solvent molecule
       NSM=NSOLV/NRAM
!$OMP parallel do private(icg)
       DO icg=1,NCAG
         INC_l(icg)=INC(icg)
       ENDDO
!$OMP parallel do private(i,icg)
       DO i=1,NSM
         icg=NCAG+i
         INC_l(icg)=INC_l(NCAG)+i*NRAM
       ENDDO
       NCAG_l=NCAG+NSM

! constraints
       NRAM_gr=NRAM
       NCONS_gr=NCONS
! IATYP
!$OMP parallel do private(is)
       DO is=1,NSPm
         iatyp(is)=999
       ENDDO
       DO i=1,NRPT
         is=cpsp(i)
         if(iatyp(is).eq.999)then
           name=fftype(iac(i))
           iatyp(is)=zatom(name)
         else
           name=fftype(iac(i))
           IF(zatom(name).ne.iatyp(is))THEN
              WRITE(6,'(4(A,I6))') ' WARNING: for is=',is,' ia=',i,
     c          ',  zatom=',zatom(name),'  but iatyp=',iatyp(is)
              call stopgm('mm_pardef','wrong iatyp assignment',
     & __LINE__,"__file__")
           ENDIF
         ENDIF
       ENDDO
       do i=NRPT+1,NRPT+NSOLV
         sn=mod(i-1-NRPT,NRAM)+1
         is=cpsp(i)
         if(iatyp(is).eq.999)then
           name=fftype(iacs(sn))
           iatyp(is)=zatom(name)
         else
           name=fftype(iacs(sn))
           IF(zatom(name).ne.iatyp(is)) THEN
             call stopgm('mm_pardef','wrong iatyp assignment',
     & __LINE__,"__file__")
           ENDIF
         endif
       enddo
! radius of the core density as used for electrostatic interaction.
!$OMP parallel do private(is)
       DO is=1,NSPm
         if(iatyp(is).le.99) then
         mm_RAGGIO(IS)= mm_core_raggio(IATYP(IS)) 
         endif
       ENDDO

       RETURN
       END
C     ==================================================================
       integer function zatom(name)
C     ==--------------------------------------------------------------==
       implicit none

       character*(*) name
       character*1 a1
       character*2 a2
       character*3 a3
       character*4 a4
       character*5 a5 

        zatom=0
        a1=name
        a2=name
        a3=name
        a4=name
        a5=name

        if(a1.eq.'H')zatom=1
        if(a1.eq.'C')zatom=6
        if(a1.eq.'N')zatom=7
        if(a1.eq.'O')zatom=8
        if(a1.eq.'F')zatom=9
        if(a1.eq.'P')zatom=15
        if(a1.eq.'S')zatom=16
        if(a1.eq.'K')zatom=19
        if(a1.eq.'V')zatom=23
        if(a1.eq.'I')zatom=53

        if(a2.eq.'CA')zatom=6 
        if(a2.eq.'CB')zatom=6      
        if(a2.eq.'CC')zatom=6                 
        if(a2.eq.'CD')zatom=6                 
        if(a2.eq.'CK')zatom=6                 
        if(a2.eq.'CM')zatom=6                 
        if(a2.eq.'CN')zatom=6                 
        if(a2.eq.'CQ')zatom=6                 
        if(a2.eq.'CR')zatom=6                 
        if(a2.eq.'CT')zatom=6                 
        if(a2.eq.'CV')zatom=6                 
        if(a2.eq.'CW')zatom=6         
        if(a2.eq.'C*')zatom=6                 
        if(a2.eq.'CY')zatom=6                 
        if(a2.eq.'CZ')zatom=6                 
        if(a2.eq.'H1')zatom=1                 
        if(a2.eq.'H2')zatom=1                 
        if(a2.eq.'h2')zatom=1      
        if(a2.eq.'H3')zatom=1                 
        if(a2.eq.'HA')zatom=1                 
        if(a2.eq.'H4')zatom=1                 
        if(a2.eq.'H5')zatom=1                 
        if(a2.eq.'h5')zatom=1                 
        if(a2.eq.'HO')zatom=1      
        if(a2.eq.'ho')zatom=1      
        if(a2.eq.'HS')zatom=1      
        if(a2.eq.'HW')zatom=1      
        if(a2.eq.'HP')zatom=1      
        if(a2.eq.'HZ')zatom=1      
        if(a2.eq.'HX')zatom=1      
        if(a2.eq.'hx')zatom=1      
        if(a2.eq.'HC')zatom=1            
        if(a2.eq.'Cl')zatom=17     
        if(a2.eq.'CL')zatom=17   
        if(a2.eq.'IM')zatom=17     
        if(a2.eq.'IB')zatom=11     
        if(a2.eq.'MG')zatom=12     
        if(a2.eq.'NA')zatom=7      
        if(a2.eq.'NB')zatom=7  
        if(a2.eq.'NC')zatom=7  
        if(a2.eq.'N2')zatom=7  
        if(a2.eq.'N3')zatom=7  
        if(a2.eq.'NT')zatom=7  
        if(a2.eq.'N*')zatom=7  
        if(a2.eq.'NY')zatom=7  
        if(a2.eq.'NL')zatom=7  
        if(a2.eq.'NR')zatom=7  
        if(a2.eq.'NZ')zatom=7  
        if(a2.eq.'NE')zatom=7  
        if(a2.eq.'O2')zatom=8  
        if(a2.eq.'OW')zatom=8  
        if(a2.eq.'OH')zatom=8  
        if(a2.eq.'oh')zatom=8  
        if(a2.eq.'OS')zatom=8  
        if(a2.eq.'os')zatom=8  
        if(a2.eq.'OM')zatom=8  
        if(a2.eq.'OA')zatom=8      
        if(a2.eq.'P4')zatom=15
        if(a2.eq.'p4')zatom=15
        if(a2.eq.'P5')zatom=15
        if(a2.eq.'p5')zatom=15
        if(a2.eq.'SH')zatom=16     
        if(a2.eq.'CU')zatom=29        
        if(a2.eq.'Cu')zatom=29        
        if(a2.eq.'FE')zatom=26        
        if(a2.eq.'Fe')zatom=26        
        if(a2.eq.'IP')zatom=11     
        if(a2.eq.'Na')zatom=11       
        if(a2.eq.'SI')zatom=14   
        if(a2.eq.'Si')zatom=14   
        if(a2.eq.'CL')zatom=17
        if(a2.eq.'Cl')zatom=17
        if(a2.eq.'AR')zatom=18   
        if(a2.eq.'Ar')zatom=18   
        if(a2.eq.'K+')zatom=19   
        if(a2.eq.'Ca')zatom=20
        if(a2.eq.'Zn')zatom=30
        if(a2.eq.'ZN')zatom=30
        if(a2.eq.'BR')zatom=35   
        if(a2.eq.'Br')zatom=35     
        if(a2.eq.'Rb')zatom=37      
        if(a2.eq.'Cs')zatom=55    

        if(a3.eq.'CH1')zatom=6       
        if(a3.eq.'CH2')zatom=6      
        if(a3.eq.'CH3')zatom=6      
        if(a3.eq.'CH4')zatom=6      
        if(a3.eq.'CR1')zatom=6      
        if(a3.eq.'DUM')zatom=1      
        if(a3.eq.'NA+')zatom=11     
        if(a3.eq.'Na+')zatom=11     
        if(a3.eq.'Cl-')zatom=17     
        if(a3.eq.'CL-')zatom=17     

        if(a4.eq.'LCH2')zatom=6
        if(a4.eq.'LCH3')zatom=6
        if(a4.eq.'CU1+')zatom=29     
        if(a4.eq.'CU2+')zatom=29     
        if(a4.eq.'CD2+')zatom=48     
        if(a4.eq.'ZN2+')zatom=30     
        if(a4.eq.'MG2+')zatom=12     
        if(a4.eq.'CA2+')zatom=20     
        if(a4.eq.'CMET')zatom=6      
        if(a4.eq.'OMET')zatom=8      
        if(a4.eq.'CCHL')zatom=6      
        if(a4.eq.'HCHL')zatom=1      
        if(a4.eq.'CCL4')zatom=6      

        if(a5.eq.'CLCHL')zatom=17     
        if(a5.eq.'SDMSO')zatom=16        
        if(a5.eq.'CDMSO')zatom=6   
        if(a5.eq.'ODMSO')zatom=8   
        if(a5.eq.'CLCL4')zatom=17   

        if(zatom.eq.0)  then
          
          WRITE(6,*)
     .    'WARNING: atom ',a5,' : unknown Force Field type !'
          if(a2.eq.'He')zatom=2
          if(a2.eq.'Li')zatom=3     
          if(a2.eq.'Be')zatom=4     
          if(a2.eq.'B ')zatom=5     
          if(a2.eq.'Ne')zatom=10    
          if(a2.eq.'Al')zatom=13    
          if(a2.eq.'Sc')zatom=21    
          if(a2.eq.'Ti')zatom=22    
          if(a2.eq.'V ')zatom=23    
          if(a2.eq.'Cr')zatom=24    
          if(a2.eq.'Mn')zatom=25    
          if(a2.eq.'Co')zatom=27    
          if(a2.eq.'Ga')zatom=31    
          if(a2.eq.'Ge')zatom=32    
          if(a2.eq.'As')zatom=33    
          if(a2.eq.'Se')zatom=34    
          if(a2.eq.'Kr')zatom=36    
          if(a2.eq.'Sr')zatom=38    
          if(a2.eq.'Y ')zatom=39    
          if(a2.eq.'Zr')zatom=40    
          if(a2.eq.'Nb')zatom=41    
          if(a2.eq.'Mo')zatom=42    
          if(a2.eq.'Tc')zatom=43    
          if(a2.eq.'Ru')zatom=44    
          if(a2.eq.'Rh')zatom=45    
          if(a2.eq.'Pd')zatom=46    
          if(a2.eq.'Ag')zatom=47    
          if(a2.eq.'Cd')zatom=48    
          if(a2.eq.'In')zatom=49    
          if(a2.eq.'Sn')zatom=50    
          if(a2.eq.'Sb')zatom=51      
          if(a2.eq.'Te')zatom=52      
          if(a2.eq.'Xe')zatom=54
          if(a2.eq.'Ba')zatom=56      
          if(a2.eq.'La')zatom=57      
          if(a2.eq.'Hf')zatom=72      
          if(a2.eq.'Ta')zatom=73      
          if(a2.eq.'W ')zatom=74      
          if(a2.eq.'Re')zatom=75      
          if(a2.eq.'Os')zatom=76      
          if(a2.eq.'Ir')zatom=77      
          if(a2.eq.'Pt')zatom=78
          if(a2.eq.'Au')zatom=79      
          if(a2.eq.'Hg')zatom=80      
          if(a2.eq.'Tl')zatom=81      
          if(a2.eq.'Pb')zatom=82      
          if(a2.eq.'Bi')zatom=83      
          if(a2.eq.'Po')zatom=84      
          if(a2.eq.'At')zatom=85      
          if(a2.eq.'Rn')zatom=86      
          if(a2.eq.'U ')zatom=92

        end if

        if(zatom.eq.0)write(6,*)'zatom:detit unknown atom ',a5

        RETURN
        END
