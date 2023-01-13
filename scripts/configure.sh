#!/bin/bash

CPMD_ROOT=${PWD}

SRC_DIR=${CPMD_ROOT}/src
MOD_DIR=${CPMD_ROOT}/modules

#--------------------------------------------------------------------------#
#Create Makefile for cpmd.x program.                                       #
#--------------------------------------------------------------------------#
# Options for CPPFLAG                                                      #
#                                                                          #
#  -D__HAS_FFT_DEFAULT Use default FFT routine in the code (from Goedecker)#
#  -D__HAS_FFT_ESSL    Use FFT ESSL (need ESSL library)                    #
#  -D__HAS_FFT_HP      Use FFT HP (Hewlet-Packard)                         #
#  -D__HAS_FFT_FFTW3   Use FFTW FFT routines                               #
#  -D__PARALLEL        For PARALLEL computers, you need to have a          #
#                      MPI parallel library installed                      #
#  -D__VECTOR          For vectorial computers                             #
#  -D__NOINT8          If the compiler does not use INTEGER*8 or           #
#                      if the function BTEST and IBTEST do not like it     #
#  -D__NO_MEMSIZE      to omit the detection of the current program size   #
#--------------------------------------------------------------------------#
#Create Makefile for cpmd.x program.

#Display all configurations
Message () {
cat << END >&2

./configure.sh  [options] Name_of_Configuration
where the name of the configuration is:
END
\ls configure >&2
cat << END >&2

To display additional information about a configuration 
type: ./configure.sh -info SUN
Use configure.sh in the directory where the SOURCE FILES are.
Ex: ./configure.sh SUN 
See the description of further options with './configure.sh -help'

To create a new configuration, add a new file in the "configure" directory.

END
}


#Help
Help () {
Message
cat << END >&2

Description of options:
   -debug      (-d) Give special options for debugging
   -help       (-h) Print this message
   -info       (-i) Print additional information about a machine configuration
   -verbose    (-v) Display the name of files for the dependencies
   -coverage   (-c) With GCC compiler only, allows for specific configuration files,
                    to generate the code coverage/profiling during an execution
   -qmmm            Generates a makefile for QMMM 
   -iphigenie       Support for external interface to iphigenie
   -omp             Enables the use of OMP instructions (if the config file allows that)
                    OMP3 is in general triggered by the config keyword OMP3_DISABLED,
                    which can either be true or false (or a script that sets true/false
                    according a certain compiler version)
   -disable_omp3    Overrides any specification (compiler/configuration file) and disables
                    OMP3 instructions
   -minpack=<PATH>  Compiles enabling MINPACK and links using the provided library
   -VAR=value       Use the value of VAR in the makefile
                    Ex: -DEST=destdir specifies where the compilation will be
   You can use:
     FFLAGS   Fortran flags
     LFLAGS   Link editor flags
     CFLAGS   C compiler flags
     CPP      Preprocessor program
     CPPFLAGS CPP flags
     FC       Fortran compiler
     LD       Link editor
     AR       Archive library maintainer
     RANLIB   Converts archives to random libraries
Note: cpmd.x will be compiled in a different directory
      than the source files.
      All compilation files will be available in the obj directory,
      while the executable will be place in the bin directory.
END
}

#By default, we try to use nawk
if [ -x /usr/bin/nawk ]; then
   AWK='/usr/bin/nawk'
else
   AWK='awk'
fi

#Check if ${AWK} does exist.
haveawk=0
echo 1 | ${AWK} '{ }' > /dev/null 2>&1 && haveawk=1 || haveawk=0
if [ ${haveawk} = 0 ]; then
  echo "The command ${AWK} does not exist!" >&2
  exit 3
fi

#Check for grep
GREP='grep'

#No argument: Error
if [ $# -eq 0 ]; then
  Message
 exit 2
fi

echo " " >&2
echo "** CPMD 4.1 configuration script" >&2
echo "** " >&2
#Is it help option or debug option
info=0
opt=0
i=1

while [ $i -le $# ];
do
  #Select the i-th argument
  eval option=\$$i
  case $option in
    -omp|-o)
      omp=1
      omp3=1
      echo "** Enabling OMP instructions (if the config file allows that)" >&2
      ;;
    -disable_omp3)
      omp3=0
      ;;
    -qmmm|-q)
      qmmm=1
      echo "** Enabling QM/MM (if Gromos and MM_Interface modules will be available)" >&2
      ;;
    -debug|-d)
      debug=1
      echo "** Debug option is used." >&2
      ;;
    -help|-h)
      Help
      exit 0
      ;;
    -coverage|-c)
      coverage=1
      echo "** Enabling coverage/profiling (if the configuration file supports it)" >&2
      ;;
    -verbose|-v)
      verbose=1
      ;;
    -info|-i)
      info=1
      ;;
    -iphigenie)
     echo "building interface to iphigenie" >&2
     iffi=1
     ;;
    -minpack=*)
      MINPACKCPP="-D__MINPACK"
      MINPACKLIB=`echo $option | cut -c2- | cut -d= -f 2`
      echo "** Linking MINPACK. Library available in ${MINPACKLIB}" >&2
      ;;
    -DEST=*)
      opt=1
      DEST=`echo $option | cut -c2- | cut -d= -f 2`
      ;;
    -*=*)
      opt=1
      ;;
    -*)
      echo "** Unknow option: $option" >&2
      Message
      exit 1
      ;;
    *)
      Configuration=$option
      echo "** Default configuration for $Configuration." >&2
      ;;
  esac
  i=`expr $i + 1`
done
echo "** " >&2
echo "** Starting now configuration.." >&2
echo " " >&2

#No configuration given
if [ -z "${Configuration}" ]; then
  echo "configure.sh: No configuration name" >&2
  Message
  exit 2
fi

#Check if the Configuration does exist.
if [ -f configure/${Configuration} ]; then
   . configure/${Configuration}
else
   echo "configure.sh: Unknown configuration '${Configuration}'" >&2
   Message
   exit 2
fi

#Print info message and exit
if [ $info -eq 1 ]
then
   echo
   echo "Additional information for configuration: ${Configuration}"
   echo
   ${GREP} '^#INFO#' configure/${Configuration} | cut -c7-
   exit 0
fi

#Print a warning for OMP3 support, when not officially tested in the configure file
CPPFLAGS_OMP3=""
if [ "${omp3}" == "1" ]
then
  echo "Compiler or Configure script reported OMP3_DISABLED with value: " ${OMP3_DISABLED} >&2
  if [ "${OMP3_DISABLED}" == "true" ]
  then
    CPPFLAGS_OMP3="-D__HASNT_OMP_COLLAPSE"
    echo " OMP3 support is disabled for this configuration file." >&2
    echo " WARNING: The clause <collapse(n)> may not work in old compilers and" >&2
    echo "          is known to have troubles on the Fujitsu compiler (K-computer" >&2
    echo "          and Fujitsu FX10 machine) and Intel Fortran Compiler version" >&2
    echo "          11 and former ones. Since ifort v.12.0.3 bugs have been fixed." >&2
    echo "          Please, refer to the discussion you see in" >&2
    echo "          http://openmp.org/forum/viewtopic.php?f=3&t=1023&start=10" >&2
    echo "          and check carefully your compiler and OS." >&2
    echo "" >&2
    echo "If you believe that your compiler fully supports OMP3, enable it by" >&2
    echo "changing the keyword OMP3_DISABLED in your configure file." >&2
    echo "" >&2
  else
    CPPFLAGS_OMP3=""
  fi
fi

#--------------------------------------------------------------------#
# Default Configurations                                             #
#--------------------------------------------------------------------#

CPPFLAGS_GROMOS='-DEWALD -DEWATCUT -DHAT_SHAPE -DUNPACKED_GRID'

#QM/MM compilation setup
if [ $qmmm ]; then
  if [ -f ${MOD_DIR}/QMMM_SOURCES ]; then
    QMMM_FLAGS='-D__GROMOS'
  else
    printf "\nThe file QMMM_SOURCES does not exist.\n" >&2
    echo "The file QMMM_SOURCES does not exist" >&2
    exit 1
  fi
else
  QMMM_FLAGS=' '
  CPPFLAGS_GROMOS=' '
  FFLAGS_GROMOS=' '
fi
if [ $iffi ]; then
  if [ -f ${MOD_DIR}/IPhigenie_Interface/IFFIINTER_SOURCES ]; then
    QMMM_FLAGS='-D__IFFIINTER'
  else
    printf "\nThe file IFFIINTER_SOURCES does not exist.\n" >&2
    echo "The file IFFIINTER_SOURCES does not exist" >&2
    exit 1
  fi
fi
      
if [ -n $DEST ]; then
    # Makefile requires DEST to be full path for correct work
    if [[ "$DEST" != /* ]]; then
        DEST=${CPMD_ROOT}/$DEST
    fi
else
    DEST=${CPMD_ROOT}
fi

OBJ_DIR=${DEST}/obj
BIN_DIR=${DEST}/bin
LIB_DIR=${DEST}/lib

mkdir -p ${BIN_DIR} ${OBJ_DIR} ${LIB_DIR}
#--------------------------------------------------------------------#
#End of Configurations                                               #
#--------------------------------------------------------------------#

# File descriptor usage:
# 0 standard input
# 1 file creation (standard output or makefile)
# 2 errors and warnings
# 3 Makefile if $makefile or &1

# Configure Makefile
echo "Creation of Makefile: ${DEST}/Makefile" >&2
exec 3>${DEST}/Makefile

cat << END >&3
#########################################################################
# Makefile for cpmd.x (plane wave electronic calculation)
# Configuration: ${Configuration}
# Creation of Makefile: `date '+%b %e %Y'`
# on `uname -a`
# Author: `who am i | cut -d' ' -f 1`
#
# To compile: go into the DEST directory (if using DEST) or stay here 
# and type make.
#########################################################################
#
SHELL = /bin/sh

#########################################################################
# project-independent stuff

CPMDROOT = ${CPMD_ROOT}
DEST     = ${DEST}
SRCDIR   = \${CPMDROOT}/src
MODDIR   = \${CPMDROOT}/modules

BINDIR   = \${DEST}/bin
OBJDIR   = \${DEST}/obj
LIBDIR   = \${DEST}/lib

MAKEFILE = \${DEST}/Makefile

TARGET        = \$(BINDIR)/cpmd.x
CPMD_LIB      = \$(LIBDIR)/libcpmd.a
END

if [ $qmmm ]; then
cat << END >&3
GROMOS_LIB    = \$(LIBDIR)/libgromos.a
INTERFACE_LIB = \$(LIBDIR)/libinterface.a
END
fi

if [ $iffi ]; then
cat << END >&3
IFFIINTER_LIB =  \$(LIBDIR)/libiffiinter.a
END
fi

cat << END >&3

.SUFFIXES: .F90 .f90 .c .o

#########################################################################
# platform-dependent stuff
# this section is built by the configure.sh script: no manual editing
# should be required.
#########################################################################
END

cat << END >&3
FFLAGS = ${FFLAGS} -I\${SRCDIR} -I\${OBJDIR}
LFLAGS = ${LFLAGS} ${MINPACKLIB}
CFLAGS = ${CFLAGS} -I\${SRCDIR}
NVCCFLAGS = ${NVCCFLAGS} -I\${SRCDIR}
CPP = ${CPP}
CPPFLAGS = ${CPPFLAGS} ${QMMM_FLAGS} ${CPPFLAGS_OMP3} ${MINPACKCPP}  -I\${SRCDIR} -D'SVN_REV="r\$(shell svnversion -n ${CPMD_ROOT})"'
NOOPT_FLAG = ${NOOPT_FLAG}
END

if [ $qmmm ]; then
cat << END >&3
CPPFLAGS_GROMOS = ${CPPFLAGS_GROMOS} ${CPPFLAGS_OMP3} 
FFLAGS_GROMOS = ${FFLAGS_GROMOS} -I\${MODDIR}
FFLAGS_GROMOS_MODULES = ${FFLAGS_GROMOS_MODULES} -I\${MODDIR}
FFLAGS += -I${MOD_DIR}
END
fi

if [ $iffi ]; then
echo "FFLAGS += -I\$(MODDIR)/IPhigenie_Interface" >&3
fi

cat << END >&3
CC = ${CC-cc}
FC = ${FC-f90 -c -O}
LD = ${LD-f90 -O}
NVCC = ${NVCC}
AR = ${AR-ar -r}
RANLIB = ${RANLIB-ranlib}
#########################################################################
END

#Check the options and personal variables 
cat << END >&3
# Personal Configuration
END
printf "# My_Conf: ">&3
i=1
while [ $i -le $# ];
do
  eval option=\$$i
  case $option in
    -*=*)
             printf "%s " $option >&3
          ;;
    *)
       ;;
  esac
  i=`expr $i + 1`
done
printf "\n">&3
echo "# All arguments: " $@ >&3
cat << END >&3
CONFARGS = $@
#########################################################################
END

#There is a Personal Configuration.
if [ $opt -eq 1 ]; then
  printf "Personal configurations..." >&2
  
  i=1
  while [ $i -le $# ];
  do
    eval option=\$$i
    case $option in
      -*=*)
	    var=`echo $option | cut -c2- | cut -d= -f 1`
	    val=`echo $option | cut -c2- | cut -d= -f 2`
	    echo "# $var = $val" >&3
	    eval $var='$val'
	    ;;
      *)
	 ;;
    esac
    i=`expr $i + 1`
  done
  
  if [ $qmmm ]; then
    echo "FC = $FC -I. -I\$(OBJDIR)/MM_Interface -I\$(OBJDIR)/Gromos" >&3
  elif [ $iffi ]; then
     echo "FC = $FC -I. -I\$(OBJDIR)/IPhigenie_Interface" >&3
  else
    echo "FC = $FC -I." >&3
  fi
fi #End of Personal Configuration.

cat << END >&3
#########################################################################
# End of Personal Configuration
#########################################################################
END

echo "done." >&2
if [ -n "${BIN_DIR}" ]; then
  echo "cpmd.x executable will be in: ${BIN_DIR}/cpmd.x ." >&2
fi

# record how we called configure for make reconfig
echo "CFGDEST = ${OBJ_DIR}" >&3
echo "CFGMACH = ${Configuration}" >&3 
if [ $qmmm ] ; then
  echo "CFGQMMM = -qmmm" >&3

  printf "QM/MM: Testing for directory modules.\n" >&2
  if [ ! -d ${MOD_DIR} ] ; then
    echo "Directory ${MOD_DIR} does not exist! Abort QM/MM installation."
    echo "Unable to configure the QM/MM module."
    echo "Contact developers at: cpmd@cpmd.org for external modules support."
    exit -9
  fi

  printf "QM/MM: Testing for directory MM_Interface.\n" >&2
  if [ ! -d ${SRC_DIR}/MM_Interface ] ; then
    if [ -d ${MOD_DIR}/MM_Interface ] ; then
       printf "QM/MM: MM_Interface exists. Continuing..\n" >& 2
    else
      echo "Directory ${MOD_DIR}/MM_Interface does not exist! Abort QM/MM installation."
      echo "Contact developers at: cpmd@cpmd.org for QM/MM support."
      exit -9 
    fi
  fi

  if [ ! -d ${OBJ_DIR}/MM_Interface ] ; then
    echo "Directory ${OBJ_DIR}/MM_Interface does not exist! Creating.."  >& 2
    mkdir ${OBJ_DIR}/MM_Interface
  fi

  printf "QM/MM: Testing for directory Gromos.\n" >& 2
  if [ ! -d ${SRC_DIR}/Gromos  ] ; then 
    if [ -d ${MOD_DIR}/Gromos ] ; then
       printf "QM/MM: GROMOS exists. Continuing..\n" >& 2
    else
      echo "Directory ${MOD_DIR}/Gromos does not exist! Abort QM/MM installation."
      echo "Contact developers at: cpmd@cpmd.org for QM/MM support."
      exit -9 
    fi
  fi

  if [ ! -d ${OBJ_DIR}/Gromos ] ; then
    echo "Directory ${OBJ_DIR}/Gromos does not exist! Creating.." >& 2
    mkdir ${OBJ_DIR}/Gromos
  fi

else 
  echo "CFGQMMM = " >&3
  rm -rf ${OBJ_DIR}/Gromos
  rm -rf ${OBJ_DIR}/MM_Interface

  rm -rf ${SRC_DIR}/Gromos
  rm -rf ${SRC_DIR}/MM_Interface
fi

if [ $iffi ]; then
    test -d ${OBJ_DIR}/IPhigenie_Interface || mkdir ${OBJ_DIR}/IPhigenie_Interface
else
   rm -rf  ${OBJ_DIR}/IPhigenie_Interface 
fi

echo "The source directory is: ${SRC_DIR}" >&2
echo "The object directory is: ${OBJ_DIR}" >&2  

#We change gromos.h if necessary
if [ -f ${OBJ_DIR}/gromos.h ]; then
  gromosfile="${OBJ_DIR}/gromos.h"
elif [ -f ${OBJ_DIR}/gromos.h ]; then
  gromosfile="${SRC_DIR}/gromos.h"
else
  comm=1
fi

if [ -n "$gromosfile" ]; then
  comm=`awk 'BEGIN{comm=0} 
     {if (index(toupper($0),"INCLUDE") != 0) {if (toupper($1)=="c" || $1=="!") comm=1}} 
     END{ print comm}' ${gromosfile}`   
fi

if [ $qmmm ]; then
  if [ $comm -eq 1 -o ! -f ${OBJ_DIR}/gromos.h ]; then
    #if gromos not writable give the good right.
    if [ -f ${OBJ_DIR}/gromos.h -a ! -w ${OBJ_DIR}/gromos.h ]; then
      echo "+ chmod u+w ${OBJ_DIR}/gromos.h" >&2
      chmod u+w ${OBJ_DIR}/gromos.h
    fi
    cat << END > ${OBJ_DIR}/gromos.h
      include 'Gromos/toposz.h'
      include 'Gromos/topoar.h'
      include 'Gromos/box.h'     
END
    echo "The file ${OBJ_DIR}/gromos.h is changed." >&2
  else
    echo "The file ${OBJ_DIR}/gromos.h is consistent with qmmm." >&2
  fi
else
  if [ $comm -eq 0 -o ! -f ${OBJ_DIR}/gromos.h ]; then
    #if gromos not writable give the good right.
    if [ -f ${SRC_DIR}/gromos.h -a ! -w ${OBJ_DIR}/gromos.h ]; then
      echo "+ chmod u+w ${OBJ_DIR}/gromos.h" >&2
      chmod u+w ${OBJ_DIR}/gromos.h
    fi
    cat << END > ${OBJ_DIR}/gromos.h
!     include 'Gromos/toposz.h'
!     include 'Gromos/topoar.h'
!     include 'Gromos/box.h'     
END
    echo "The file ${OBJ_DIR}/gromos.h is changed." >&2
  else
    echo "The file ${OBJ_DIR}/gromos.h is consistent with normal qm." >&2
  fi
fi

#Add OBJECTS to the makefile
if [ -f ${SRC_DIR}/SOURCES ]; then
  printf "Add OBJECTS (object files)..." >&2
  cat << END  >&3
################################################################################
# files 
# load SRC_ variables from SOURCES file: correct module compilation order
include \$(SRCDIR)/SOURCES

OBJ_MODULES = \$(SRC_MODULES:%.F90=%.o)
OBJ_F90     = \$(SRC_F90:%.F90=%.o)
OBJ_CC      = \$(SRC_CC:%.c=%.o)
ifneq (\$(NVCC),)
OBJ_CU      = \$(SRC_CU:%.cu=%.o)
else
OBJ_CU      = 
endif
OBJ_CPMD     = \$(OBJ_MODULES) \$(OBJ_F90) \$(OBJ_CC) \$(OBJ_CU)

END
  echo "done."  >&2
else
  echo "The file SOURCES does not exist" >&2
  exit 1
fi

#QMMM OBJECTS
if [ $qmmm ]; then
  printf "Add QMMM_OBJECTS (object files)..." >&2
  cat << END >&3
# load QMMM SRC_ variables from SOURCES file
include \$(MODDIR)/QMMM_SOURCES

OBJECTS_MODULES_GROMOS = \$(GROMOS_MODULES:%.F90=%.o) 
OBJECTS_GROMOS    = \$(GROMOS_SRC:%.F=%.o)
OBJECTS_INTERFACE = \$(INTERFACE_SRC:%.F=%.o)

END
  echo "done."  >&2
fi

if [ $iffi ]; then
    echo "include \$(MODDIR)/IPhigenie_Interface/IFFIINTER_SOURCES" >&3
    echo "OBJECTS_IFFIINTER = \$(IFFIINTER_SRC:%.F90=%.o)"     >&3
fi

printf "Add explicit rules..." >&2
cat << END >&3
################################################################################
# Explicit Rules
################################################################################


default, all:
	+\$(MAKE) -C \$(OBJDIR) -f \$(MAKEFILE) \$(TARGET)

lib:
	+\$(MAKE) -C \$(OBJDIR) -f \$(MAKEFILE) \$(CPMD_LIB)
END
if [ $qmmm ]; then
cat << END >&3
	+\$(MAKE) -C \$(OBJDIR) -f \$(MAKEFILE) \$(GROMOS_LIB)
	+\$(MAKE) -C \$(OBJDIR) -f \$(MAKEFILE) \$(INTERFACE_LIB)
END
fi
if [ $iffi ]; then
cat << END >&3
	+\$(MAKE) -C \$(OBJDIR) -f \$(MAKEFILE) \$(IFFIINTER_LIB)
END
fi

cat << END >&3

distclean:  realclean
	rm -f  \$(TARGET)
	rm -rf \$(DEST)/Makefile

realclean:  clean
	rm -rf \$(OBJDIR) \$(BINDIR) \$(LIBDIR)

clean:  purge
	rm -f \$(CPMD_LIB) \$(GROMOS_LIB) \$(INTERFACE_LIB) \$(OBJDIR)/*.o \$(OBJDIR)/*.mod \$(OBJDIR)/*.f90 \$(OBJDIR)/*.c \$(OBJDIR)/*.cu
	rm -f \$(OBJDIR)/Gromos/*.o \$(OBJDIR)/Gromos/*.f*
	rm -f \$(OBJDIR)/MM_Interface/*.o \$(OBJDIR)/MM_Interface/*.f
	rm -f \$(OBJDIR)/IPhigenie_Interface/*.{o,f}

purge: 
	rm -f \$(SRCDIR)/*~

reconfigure:
	\$(shell cd \$(CPMDROOT); ./configure.sh \$(CONFARGS) )

%.o:    \$(SRCDIR)/%.c
	( cd \$(SRCDIR); cp \$(<F) \$(OBJDIR)/\$(<F); cd \$(OBJDIR) )
	\$(CC) -c \$(CFLAGS) \$(CPPFLAGS) -o \$@ \$(<F)

ifneq (\$(NVCC),)
%.o:    \$(SRCDIR)/%.cu
	( cd \$(SRCDIR); cp \$(<F) \$(OBJDIR)/\$(<F); cd \$(OBJDIR) )
	\$(NVCC) -c \$(NVCCFLAGS) \$(CPPFLAGS) -o \$@ \$(<F)
endif

ifneq (\$(CPP),)
%.f90:  \$(SRCDIR)/%.F90
	( cd \$(SRCDIR); \$(CPP) \$(CPPFLAGS) \$(<F) \$(OBJDIR)/\$(@F); cd \$(OBJDIR) )

%.mod.o: %.mod.f90
	\$(FC) -c \$(FFLAGS) -o \$@ \$(@F:.o=.f90)

%.o:    \$(OBJ_MODULES) %.f90 
	\$(FC) -c \$(FFLAGS) -o \$@ \$(@F:.o=.f90)
else
%.mod.o: \$(SRCDIR)/%.mod.F90
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -o \$@ \$(SRCDIR)/\$(@F:.o=.F90)

%.o:    \$(OBJ_MODULES) %.f90 
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -o \$@ \$(SRCDIR)/\$(@F:.o=.F90)
endif
END
if [ $qmmm ]; then
cat << END >&3
\$(TARGET): \$(CPMD_LIB) \$(GROMOS_LIB) \$(INTERFACE_LIB) timetag.o cpmd.o
	\$(LD) \$(FFLAGS) -o \$(TARGET) timetag.o cpmd.o \$(CPMD_LIB) \$(GROMOS_LIB) \$(INTERFACE_LIB) \$(LFLAGS)
	@ls -l \$(TARGET)
	@echo "Compilation done."
END
elif [ $iffi ]; then
cat << END >&3
ifneq (\$(CPP),)
\$(IFFIINTER_SRC:.F90=.f90):  \$(MODDIR)/IPhigenie_Interface/\$(@F:.f90=.F90)
	( \$(CPP) \$(CPPFLAGS)   \$(MODDIR)/IPhigenie_Interface/\$(@F:.f90=.F90) \$(OBJDIR)/IPhigenie_Interface/\$(@F);  )

\$(IFFIINTER_SRC:.F90=.mod.o): %.mod.f90
	 \$(FC) -c \$(FFLAGS) -o \$@ \$(OBJDIR)/IPhigenie_Interface/\$(@F:.o=.f90)

\$(IFFIINTER_SRC:.F90=.o):    \$(OBJ_MODULES) \$(@F:.o=.f90)
	\$(FC) -c \$(FFLAGS) -o \$@  \$(OBJDIR)/IPhigenie_Interface/\$(@F:.o=.f90)
else
\$(IFFIINTER_SRC:.F90=.mod.o):  \$(MODDIR)/IPhigenie_Interface/\$(@F:.mod.o=.F90)
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -o \$@ \$(MODDIR)/IPhigenie_Interface/\$(@F:.o=.F90)
\$(IFFIINTER_SRC:.F90=.o):    \$(OBJ_MODULES) \$(@F:.o=.f90)
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -o \$@ \$(MODDIR)/IPhigenie_Interface/\$(@F:.o=.F90)
endif

# create cpmd.x and independent libiffiinter.a
\$(TARGET): \$(CPMD_LIB) \$(IFFIINTER_LIB) timetag.o cpmd.o
	\$(LD) \$(FFLAGS) -o \$(TARGET) timetag.o cpmd.o \$(CPMD_LIB) \$(LFLAGS)
	@ls -l \$(TARGET)
	@echo "Compilation done."
END
else
cat << END >&3
\$(TARGET): \$(CPMD_LIB) timetag.o cpmd.o
	\$(LD) \$(FFLAGS) -o \$(TARGET) timetag.o cpmd.o \$(CPMD_LIB) \$(LFLAGS)
	@ls -l \$(TARGET)
	@echo "Compilation done."
END
fi
cat << END >&3

\$(CPMD_LIB): \$(OBJ_CPMD)
	\$(AR) \$(CPMD_LIB) \$(OBJ_CPMD)
	\$(RANLIB) \$(CPMD_LIB)

################################################################################
#       SPECIAL RULES ONLY
################################################################################
# IRAT module is compiled separately to have -D_IRAT_
irat.mod.o: \$(SRCDIR)/irat.mod.F90
ifneq (\$(CPP),) 
	\$(CPP) \$(CPPFLAGS) -D_IRAT_=${IRAT} \$< \$(@F:.o=.f90)
	\$(FC) -c \$(FFLAGS) -o \$@ \$(@F:.o=.f90)
else
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -D_IRAT_=${IRAT} -o \$@ \$(SRCDIR)/\$(@F:.o=.F90)
endif

# timetag is compiled separately to have the compilation date always up2date
.PHONY: timetag.o
timetag.o:
ifneq (\$(CPP),) 
	\$(CPP) \$(CPPFLAGS) \$(SRCDIR)/timetag.F90 timetag.f90
	\$(FC) -c \$(FFLAGS) -o timetag.o timetag.f90
else
	\$(FC) -c \$(FFLAGS) \$(CPPFLAGS) -o timetag.o \$(SRCDIR)/timetag.F90 
endif

.SECONDARY: \$(SRC_ALL:%.F90=%.f90) cpmd.f90 # do not delete intermediate .f90 files
.SECONDARY: \$(GROMOS_SRC:%.F=%.f) \$(INTERFACE_SRC:%.F=%.f)
.SECONDARY: \$(GROMOS_MODULES:%.mod.F90=%.mod.f90) 

END

if [ $qmmm ]; then
cat << END >&3
################################################################################
#       QM/MM rules
################################################################################
# TODO write special rules for forcematch_ only
\$(OBJ_MODULES): \$(OBJECTS_MODULES_GROMOS)

ifneq (\$(CPP),)
\$(OBJECTS_GROMOS:.o=.f): 
	rm -f \$@
	\$(CPP) \$(CPPFLAGS_GROMOS) \$(MODDIR)/\$(@:.f=.F) \$@

\$(OBJECTS_MODULES_GROMOS:.o=.f90):
	rm -f \$@
	\$(CPP) \$(CPPFLAGS_GROMOS) \$(MODDIR)/\$(@:.f90=.F90) \$@

\$(OBJECTS_MODULES_GROMOS): \$(OBJECTS_MODULES_GROMOS:.o=.f90)

	\$(FC) -c \$(FFLAGS_GROMOS_MODULES) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ -o \$@ \$(@:.o=.f90) 

\$(OBJECTS_GROMOS): \$(OBJECTS_MODULES_GROMOS)
	\$(FC) -c \$(FFLAGS_GROMOS) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ -o \$@ \$(@:.o=.f)

else
\$(OBJECTS_MODULES_GROMOS):
	\$(FC) -c \$(FFLAGS_GROMOS_MODULES) \$(CPPFLAGS_GROMOS) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ -o \$@ \$(MODDIR)/\$(@:.o=.F90) 

\$(OBJECTS_GROMOS): \$(OBJECTS_MODULES_GROMOS) 
	\$(FC) -c \$(FFLAGS_GROMOS) \$(CPPFLAGS_GROMOS) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ -o \$@ \$(MODDIR)/\$(@:.o=.F) 
endif

\$(GROMOS_LIB): \$(OBJECTS_MODULES_GROMOS) \$(OBJECTS_GROMOS)
	\$(AR) \$(GROMOS_LIB) \$(OBJECTS_MODULES_GROMOS) \$(OBJECTS_GROMOS)
	\$(RANLIB) \$(GROMOS_LIB)

ifneq (\$(CPP),)
\$(OBJECTS_INTERFACE:.o=.f): 
	rm -f \$@
	\$(CPP) \$(CPPFLAGS) \$(MODDIR)/\$(@:.f=.F) \$@

\$(OBJECTS_INTERFACE): 
	\$(FC) -c \$(FFLAGS_GROMOS) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ \$< -o \$@
else
\$(OBJECTS_INTERFACE): 
	\$(FC) -c \$(FFLAGS_GROMOS) \$(CPPFLAGS_GROMOS) -I\$(MODDIR)/Gromos/ -I\$(MODDIR)/MM_Interface/ -o \$@ \$(MODDIR)/\$(@:.o=.F)
endif

\$(INTERFACE_LIB): \$(OBJECTS_INTERFACE)
	\$(AR) \$(INTERFACE_LIB) \$(OBJECTS_INTERFACE)
	\$(RANLIB) \$(INTERFACE_LIB)
END
fi
echo "done."  >&2

if [ $iffi ]; then
cat << END >&3
\$(IFFIINTER_LIB): \$(OBJ_CPMD) \$(OBJECTS_IFFIINTER) timetag.o
	\$(AR) \$(IFFIINTER_LIB) \$(OBJ_CPMD) \$(OBJECTS_IFFIINTER) timetag.o
	\$(RANLIB) \$(IFFIINTER_LIB)
END
fi

cat << END >&3
################################################################################
# Module dependencies
# 
END

printf "Create CPMD dependencies..." >&2
cd src >&2
FortranFiles=`ls *.F90  *.c `
${GREP} -i '^ *module ' *.F90 |\
  ${GREP} -i -v procedure |\
  ${GREP} -i -v iso_c_binding | \
                           ${AWK} -F":| " 'BEGIN{suffix=".o"}
                                   {modN=$3; fileN=substr($1,1,length($1)-4);
                                    line = sprintf("%s.mod:",modN);
                                    printf "%-16s%s%s\n	@true\n", line,fileN,suffix}' >&3


cat << END >&3
################################################################################
# Object dependencies: CPMD
# 
END

for name in ${FortranFiles}
do
  if [ $verbose ]; then
    printf "[%s]" $name
  fi
  ${AWK} -v qmmm=${qmmm}  '
       NR==1 { # Add here any exclusion to external modules
               SkipInclude["xc_f03_lib_m"] = 0;
               SkipInclude["mpif.h"] = 0;
               SkipInclude["mpi"] = 0;
               SkipInclude["rhjsx.inc"] = 0;
               SkipInclude["uhjsx.inc"] = 0;
               if (qmmm != 1) {
                  SkipInclude["coordsz"] = 0;
               }
               MaxLength=60;
               ll = length(FILENAME);
               #find suffix
               i=ll
               while (substr(FILENAME,i,1) != ".") {i-=1};
               prefix = substr(FILENAME,1,i-1);
               suffix = substr(FILENAME,i,ll-i+1);
               src  = "$(SRCDIR)";
               dest = "$(OBJDIR)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F90") {
                 line = sprintf("%s.f90:", prefix);
                 printf "%-16s%s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f90";
               }
               line = sprintf("%s.o:",prefix);
               if (suffix == ".c") {
                 line = sprintf("%-16s%s/%s%s",line,src,prefix,suffix);
               } else {
                 line = sprintf("%-16s%s%s",line,prefix,suffix);
               }
       }
       found = 0;
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               chck = Name;
               found = 1;
       }
       /^[ ]*USE / || /^[ ]*use / {
               split($2,a,",");
               chck = tolower(a[1]);
               Name = sprintf("%s.mod",a[1]);
               found = 1;
       }
       { if ( found == 1 ) {
               no = 0;
               for ( Include in SkipInclude ) {
                 if ( Include == chck ) {
                   no = 1;
                   break;
                 }
               }
               if((Name == "gromos.h") && (qmmm == 1)) 
                 Name=sprintf("$(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h")
               if (no != 1) {
                 SkipInclude[Name] = 0;
                 lline = length(line);
                 if (lline+2 >= MaxLength) {
                   printf "%s \\\n", line;
                   line = sprintf("%-15s"," ");
                 }
                 line=sprintf("%s %s",line,Name)
               }
         }
       }
       END   { printf "%s\n", line }' ${name} >&3

       # append compilation rules for objects with disabled optimization.
       if [ -n "$NOOPT_OBJS" ] 
       then
           for file in $NOOPT_OBJS
           do 
               if [ ${file%.o} == ${name%\.*} ]
               then
       	         echo "	\$(FC) -c \$(NOOPT_FLAG) -o \$@ \$(@F:.o=.f90)" >&3
                 echo "# Optimizations disabled for ${name}." >&3
               fi
           done
       fi
       echo "" >&3
done                               
cd ../ >&2
echo "done."  >&2


if [ $qmmm ]; then
  cd modules >&2
  FortranFiles=`ls Gromos/*.F Gromos/*.F90`
  printf "Create GROMOS dependencies..." >&2 
  cd Gromos >&2
  ${GREP} -i '^ *module ' *.F90 |\
  ${GREP} -i -v procedure |\
  ${GREP} -i -v iso_c_binding | \
                           ${AWK} -F":| " 'BEGIN{suffix=".o"; src  = "Gromos/"}
                                   {modN=$3; fileN=substr($1,1,length($1)-4);
                                    line = sprintf("%s.mod:",modN);
                                    printf "%-16s%s%s%s\n	@true\n", line,src,fileN,suffix}' >&3
  cd ../>&2
cat << END >&3
################################################################################
# Object dependencies: QM/MM modules
# 

END

  for name in ${FortranFiles}
  do
    if [ $verbose ]; then
      printf "[%s]" $name
    fi
    awk 'NR==1 { # Add here any exclusion to external modules
               SkipInclude["xc_f03_lib_m"] = 0;
               SkipInclude["mpif.h"] = 0;
               SkipInclude["mpi"] = 0;
               MaxLength=60;
               ll = length(FILENAME);
               #find suffix
               i=ll
               while (substr(FILENAME,i,1) != ".") {i-=1};
               prefix = substr(FILENAME,1,i-1);
               suffix = substr(FILENAME,i,ll-i+1);
               src  = "$(MODDIR)";
               dest = "$(OBJDIR)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F90") {
                 line = sprintf("%s.f90:", prefix);
                 printf "%-16s%s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f90";
               }
               if (suffix == ".F") {
                 line = sprintf("%s.f:", prefix);
                 printf "%-16s%s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f";
               }
               src    = "$(MODDIR)/Gromos";
               line = sprintf("%s.o:",prefix);
               line = sprintf("%-16s%s%s",line,prefix,suffix);
             }
       found = 0;
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               chck = Name;
               found = 1;
               dir=src;
       }
       /^[ ]*USE / || /^[ ]*use / {
               split($2,a,",");
               Name = sprintf("%s.mod",a[1]);
               chck = tolower(a[1]);
               found = 1;
               dir="";
       }
       { if ( found == 1 ) {
               no = 0;
               if(Name=="") no=1;
               if(Name==" ") no=1;
               for ( Include in SkipInclude ) {
                 if ( Include == chck ) { 
                   no = 1;
                   break;
                 }
               }
               if (no != 1) {
                 SkipInclude[Name] = 0;
                 lline = length(line);
                 if (lline+2 >= MaxLength) {
                   printf "%s \\\n", line;
                   line = sprintf("%-15s"," ");
                 }
                 if (Name == "irat.inc") {
                   line=sprintf("%s %s/%s",line,dest,Name)
                 } else {
                   src_pp=dir;
                   Name_pp[2]=Name;
                   if (match(Name,"/") != 0) {src_pp="$(MODDIR)"};
                   if (dir != "") {line=sprintf("%s %s/%s",line,src_pp,Name_pp[2])}
                   else {line=sprintf("%s %s",line,Name_pp[2])}
                 }
               }
         }
       }
       END   { printf "%s\n\n", line }' ${name} >&3
  done
  echo "done."  >&2


  FortranFiles=`ls MM_Interface/*.F`
  printf "Create MM_INTERFACE dependencies..." >&2 

  for name in ${FortranFiles}
  do
    if [ $verbose ]; then
      printf "[%s]" $name
    fi
    awk 'NR==1 { # Add here any exclusion to external modules
               SkipInclude["xc_f03_lib_m"] = 0;
               SkipInclude["mpif.h"] = 0;
               SkipInclude["mpi"] = 0;
               MaxLength=60;
               ll = length(FILENAME);
               prefix = substr(FILENAME,1,ll-2);
               suffix = substr(FILENAME,ll-1,2);
               src    = "$(MODDIR)/MM_Interface";
               dest = "$(OBJDIR)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F") {
                 src    = "$(MODDIR)";
                 line = sprintf("%s.f:", prefix);
                 printf "%-16s%s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f";
                 src    = "$(MODDIR)/MM_Interface";
               }
               line = sprintf("%s.o:", prefix);
               line = sprintf("%-16s%s%s",line,prefix,suffix);
             }
       found = 0;
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               chck = Name;
               found = 1;
               dir=src;
       }
       /^[ ]*USE / || /^[ ]*use / {
               split($2,a,",");
               Name = sprintf("%s.mod",a[1]);
               chck = tolower(a[1]);
               found = 1;
               dir="";
       }
       { if ( found == 1 ) {
               no = 0;
               for ( Include in SkipInclude ) {
                 if (Include == chck ) { 
                   no = 1;
                   break;
                 }
               }
               if (no != 1) {
                 SkipInclude[Name] = 0;
                 lline = length(line);
                 if (lline+2 >= MaxLength) {
                   printf "%s \\\n", line;
                   line = sprintf("%-15s"," ");
                 }
                 if (Name == "irat.inc") {
                   line=sprintf("%s %s/%s",line,dest,Name)
                 } else {
                   src_pp=dir;
                   Name_pp[2]=Name;
                   if (match(Name,"/") != 0) {src_pp="$(MODDIR)"};
                   if (dir != "") {line=sprintf("%s %s/%s",line,src_pp,Name_pp[2])}
                   else {line=sprintf("%s %s",line,Name_pp[2])}
                 }
               }
          }
       }
       END   { printf "%s\n\n", line }' ${name} >&3
  done
  echo "done."  >&2
fi
#end if [ $qmmm ]

if [ $iffi ] ; then
    printf "Create iphigenie interface dependencies..." >&2
    cd $MOD_DIR/IPhigenie_Interface >&2
    FortranFiles=`ls *.F90 `
  ${GREP} -i '^ *module ' *.F90 |\
  ${GREP} -i -v procedure |\
  ${GREP} -i -v iso_c_binding | \
                           ${AWK} -F":| " 'BEGIN{suffix=".o"}
                                   {modN=$3; fileN=substr($1,1,length($1)-4);
                                    line = sprintf("%s.mod:",modN);
                                    printf "%-16s%s%s\n	@true\n", line,fileN,suffix}' >&3
cat << END >&3
################################################################################
# Object dependencies: Iphigenie Interface
# 
END

for name in ${FortranFiles}
do
  if [ $verbose ]; then
    printf "[%s]" $name
  fi
  ${AWK} -v qmmm=${qmmm}  '
       NR==1 { # Add here any exclusion to external modules
               SkipInclude["xc_f03_lib_m"] = 0;
               SkipInclude["mpif.h"] = 0;
               SkipInclude["mpi"] = 0;
               SkipInclude["rhjsx.inc"] = 0;
               SkipInclude["uhjsx.inc"] = 0;
               if (qmmm != 1) {
                  SkipInclude["coordsz"] = 0;
               }
               MaxLength=60;
               ll = length(FILENAME);
               #find suffix
               i=ll
               while (substr(FILENAME,i,1) != ".") {i-=1};
               prefix = substr(FILENAME,1,i-1);
               suffix = substr(FILENAME,i,ll-i+1);
               src  = "$(MODDIR)/IPhigenie_Interface";
               dest = "$(OBJDIR)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F90") {
                 line = sprintf("%s.f90:", prefix);
                 printf "%-16s%s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f90";
               }
               line = sprintf("%s.o:",prefix);
               if (suffix == ".c") {
                 line = sprintf("%-16s%s/%s%s",line,src,prefix,suffix);
               } else {
                 line = sprintf("%-16s%s%s",line,prefix,suffix);
               }
       }
       found = 0;
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               chck = Name;
               found = 1;
       }
       /^[ ]*USE / || /^[ ]*use / {
               split($2,a,",");
               chck = tolower(a[1]);
               Name = sprintf("%s.mod",a[1]);
               found = 1;
       }
       { if ( found == 1 ) {
               no = 0;
               for ( Include in SkipInclude ) {
                 if ( Include == chck ) {
                   no = 1;
                   break;
                 }
               }
               if((Name == "gromos.h") && (qmmm == 1)) 
                 Name=sprintf("$(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h")
               if (no != 1) {
                 SkipInclude[Name] = 0;
                 lline = length(line);
                 if (lline+2 >= MaxLength) {
                   printf "%s \\\n", line;
                   line = sprintf("%-15s"," ");
                 }
                 line=sprintf("%s %s",line,Name)
               }
         }
       }
       END   { printf "%s\n", line }' ${name} >&3

       # append compilation rules for objects with disabled optimization.
       if [ -n "$NOOPT_OBJS" ] 
       then
           for file in $NOOPT_OBJS
           do 
               if [ ${file%.o} == ${name%\.*} ]
               then
                 echo " \$(FC) -c \$(NOOPT_FLAG) -o \$@ \$(@F:.o=.f90)" >&3
                 echo "# Optimizations disabled for ${name}." >&3
               fi
           done
       fi
       echo "" >&3
done                               
cd - >&2
echo "done."  >&2
fi

cat << END >&3

# create a regular timestamped cpmd package from the cvs tree
gtar-cpmd  :
	@( d=\`date +%Y%m%d\` ; b=\`basename \$\$PWD\` ; cd \${CPMDROOT}/.. ;	\\
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude .svn	\\
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log 		\\
	--exclude \*.out --exclude \*.prj --exclude \*.chk			\\
	--exclude modules/MM_Interface --exclude modules/Gromos 		\\
        --exclude modules/IPhigenie_Interface					\\
	--exclude \*.orig --exclude \*.rej -zvvf \$\$b-\$\$d.tar.gz \$\$b &&	\\
	echo successfully created \$\$b-\$\$d-tar.gz ; cd \$\$b )

# create a QM/MM timestamped cpmd package from the cvs tree
gtar-qmmm :
	@( d=\`date +%Y%m%d\` ; b=\`basename \$\$PWD\` ; cd \${CPMDROOT}/.. ;	\\
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude .svn	\\
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log 		\\
	--exclude \*.out --exclude \*.prj --exclude \*.chk			\\
	--exclude \*.orig --exclude \*.rej -zvvf \$\$b-\$\$d.tar.gz \$\$b &&	\\
	echo successfully created \$\$b-\$\$d-tar.gz ; cd \$\$b )

END

echo "Configuration successful!"  >&2
echo "Go to $DEST directory and type make!"  >&2
echo "O.K." >&2
exit 0
