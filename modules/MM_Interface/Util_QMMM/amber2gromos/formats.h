C $Id: formats.h,v 1.1 2006-12-27 11:25:59 itavern Exp $ -*-fortran-*-

C some strings used as format statements for writing out errors
C these strings are initialized in a BLOCKDATA statement in fileio.f

      CHARACTER*(40) FMNII
      CHARACTER*(40) FMNIR
      CHARACTER*(32) FMGEI
      CHARACTER*(32) FMGER
      CHARACTER*(32) FMGTI
      CHARACTER*(32) FMGTR
      CHARACTER*(40) FMBTI
      CHARACTER*(46) FMBTII
      CHARACTER*(44) FMBTR
      CHARACTER*(28) FMBOO
      CHARACTER*(64) FMSHKP
      CHARACTER*(64) FMSHKS
      COMMON /FMTSTR/FMNII,FMNIR,FMGEI,FMGER,FMGTI,FMGTR,
     $     FMBTI,FMBTR,FMBOO,FMBTII,FMSHKP,FMSHKS
