#Defaults: add overrides below
LD=$(FC)
###########################################
#   SELECT HERE BELOW THE COMPILER YOU
#   NEED ACCORDING TO YOUR OWN MACHINE
###########################################
#- RECENT INTEL FORTRAN COMPILER --
FC= ifort -c -O2 -unroll -fixed
CPP= fpp -C -P 
LD= ifort -static-intel -O2
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
CPPFLAGS = -D_AMBER12
#-------INTEL FORTRAN COMPILER-------
#FC= ifort -c -pc64 -O2 -unroll
#CPP = /lib/cpp -P -C -traditional
#LD= ifort -static-intel -O2
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------PGF90--------------
#FC= pgf90 -fast -pc=64 -c
#CPP = /lib/cpp -P -C -traditional
#LD= pgf90 -pc=64
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#-------------GFORTRAN---------------
#FC= gfortran -Wall -O2 -c
#CPP = /lib/cpp -P -C -traditional
#LD= gfortran
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------IBM XLF------------
#FC= xlf -O2 -q64 -qarch=auto -qdpc -qstrict -c
#CPP = /usr/ccs/lib/cpp -P -C
#LD= xlf -O2 -q64 -qstrict
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#-------SGI-Altix-IFC--------------
#FC= ifort -c -O2 -r8 -align dcommons -fpp
#CPP = /usr/bin/cpp -E -P -traditional
#LD= ifort -i-static -O2
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------MACOSX-XLF---------
#FC= xlf -O2 -qstrict -qarch=auto -qdpc -c
#CPP = /opt/ibmcmp/xlf/8.1/exe/cpp -P -C
#LD= xlf -O2 -qstrict -qarch=auto -qdpc
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------Fujitsu-VPP500-----
#FC= frt -c -Aw -Sw -O5 -Kfast -X9 -Fixed -No -Os,-X,-l
#CPP = /lib/cpp -P -C
#LD= frt -Aw -Sw -O5 -Kfast -X9 -Fixed -No -Os,-X,-l
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------Fujitsu-FX10-----
#FC= frtpx -c -Aw -O3 -Kfast -X9 -Free -Cpp
#CPP = /lib/cpp -P -C
#LD= frtpx -Aw -O3 -Kfast -X9 -Free -Cpp
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------NEC-SX6-ES---------
##X-compiler
#FC= esf90 -c -ew -float0 -f0 -Nftrace -V -R5 -P static \
#   -Cvopt -Cssafe -pi -Wf"-pvctl noassume loopcnt=99999999 \
#   fullmsg vr256 vwork=static vworksz=4M novlchk noverrchk"
#CPP = /usr/ccs/lbin/cpp -P -T
#LD= esf90 -ew -float0 -f0 -Nftrace -V -R5 -P static \
#   -Cvopt -Cssafe -pi -Wf"-pvctl noassume loopcnt=99999999 \
#   fullmsg vr256 vwork=static vworksz=4M novlchk noverrchk"
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#---------------NEC-SX6------------
##Standard f90 compiler
#FC= f90 -c -ew -float0 -f0 -Nftrace -V -R5 -P static \
#   -Cvopt -Cssafe -pi -Wf"-pvctl noassume loopcnt=99999999 \
#   fullmsg vr256 vwork=static vworksz=4M novlchk noverrchk"
#CPP = /usr/ccs/lbin/cpp -P -T
#LD= f90 -ew -float0 -f0 -Nftrace -V -R5 -P static \
#   -Cvopt -Cssafe -pi -Wf"-pvctl noassume loopcnt=99999999 \
#   fullmsg vr256 vwork=static vworksz=4M novlchk noverrchk"
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6
##For AMBER-7,8,9: 
#CPPFLAGS = -D_AMBER7
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12
#----------Hitachi SR8000----------
#FC= xf90 -c -Oss -procnum=8 -fixed -conti199 -pvfunc=2 -noscope \
#    -noloopsplit -norapidcall -msg=w -intptr -nolimit -64
#CPP = xcpp -PC
#LD= xf90 -Oss -procnum=8 -fixed -conti199 -pvfunc=2 -noscope \
#    -noloopsplit -norapidcall -msg=w -intptr -nolimit -64
##For AMBER-6 and previous versions: 
#CPPFLAGS = -D_AMBER6 -D__SR8000
##For AMBER-7,8,9: 
#CPPFLAGS = -D_AMBER7 -D__SR8000
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12 -D__SR8000
#----------Hitachi SR11000-K1------
#FC= f90 -c -Oss -64 -model=K1 -fullmsg -intptr -nolimit -loglist
#CPP = /opt/optc/bin/cpp -PC
#LD= f90 -Oss -64 -model=K1 -fullmsg -intptr -nolimit -loglist
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6 -D__SR8000 -D__SR11000
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7 -D__SR8000 -D__SR11000
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12 -D__SR8000 -D__SR11000
#----------Hitachi SR11000-J1------
#FC= f90 -c -Oss -64 -model=J1 -fullmsg -intptr -nolimit -loglist
#CPP = /opt/optc/bin/cpp -PC
#LD= f90 -Oss -64 -model=J1 -fullmsg -intptr -nolimit -loglist
##For AMBER-6 and previous versions:
#CPPFLAGS = -D_AMBER6 -D__SR8000 -D__SR11000
##For AMBER-7,8,9:
#CPPFLAGS = -D_AMBER7 -D__SR8000 -D__SR11000
##For AMBER-12 and next version:
#CPPFLAGS = -D_AMBER12 -D__SR8000 -D__SR11000
###################################

#OBJFILE=blockio.o  conv.o  fileio.o  rdtopo.o  wrtopo.o 
OBJFILE=blockio.o  fileio.o  rdtopo.o  wrtopo.o 

default: ambertogromos.x

# pattern rules:
.SUFFIXES : .F
#
.F.f :
	$(CPP) $(CPPFLAGS) $< $@
#
.f.o:
	$(FC) $(FLAG) -c $<

conv.f: conv.F
	$(CPP) $(CPPFLAGS) $< $@

ambertogromos.x : $(OBJFILE) conv.o
	$(LD) -o  $@  conv.o $(OBJFILE) 

#
gtar: clean
	@( d=`date +%Y%m%d` ; b=`basename $$PWD` ; cd .. ;		   \
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude CVS \
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log	   \
        --exclude \*.out --exclude \*.prj --exclude \*.chk 		   \
	--exclude \*.orig --exclude \*.rej -zvvf $$b-$$d.tar.gz $$b &&	   \
        echo successfully created $$b-$$d-tar.gz ; cd $$b )

	
clean :
	rm -f $(OBJFILE) *.o *.f *.L *~
#
blockio.o:blockio.f dataid.h fileio.h \
                toposz.h topoar.h \
                coordsz.h disre.h \
                formats.h

fileio.o:fileio.f units.h dataid.h \
                fileio.h

rdtopo.o:rdtopo.f units.h toposz.h \
                topoar.h fileio.h \
                topblock.h formats.h 

wrtopo.o:wrtopo.f toposz.h topoar.h \
                topblock.h 

conv.o:conv.f toposz.h topoar.h \
       md.h coordsz.h coordar.h

