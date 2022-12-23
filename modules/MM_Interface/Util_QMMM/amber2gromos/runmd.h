C $Id: runmd.h,v 1.1 2006-12-27 11:32:24 itavern Exp $  -*-fortran -*-

C define some constants used in runmd.f

COMMVAR NFTTO,NFTPIR,NFTPCM,NFTSLV,NFT4D,NFTS4D,NFTP4D,NFTMAX
C*****indices into temperature, ekin and scaling related arrays
C     used in subroutine L<RUNMD>.
C
C     NFTTO :total kinetic energy in 4D
C     NFTPIR: solute internal and rotational (in 3D)
C     NFTPCM: solute c.o.m. translational
C     NFTSLV: solvent total 3D
C     NFT4D : 4th D (solvent and solute)
C     we also need (internally, never reported)
C     NFTS4D: 4th D solvent
C     NFTP4D: 4th D  solute
C     NFTMAX: the size of the temp array
C
      INTEGER NFTTO,NFTPIR,NFTPCM,NFTSLV,NFT4D,NFTS4D,NFTP4D,NFTMAX

      PARAMETER (NFTPIR = 1,NFTPCM = 2,NFTSLV = 3)
      PARAMETER (NFT4D = 4, NFTTO = 5,NFTS4D = 6, NFTP4D = 7)
      PARAMETER (NFTMAX = NFTP4D)
COMMEND

COMMVAR MXBATH
C     size of arrays used in unravelling the L<NTT> switches
C     in subroutine L<RUNMD>.
C
      INTEGER MXBATH
      PARAMETER (MXBATH = 4)
COMMEND


COMMVAR MAXFMT,FMNDIM
C     this is a format string which is used to
C     print out a host of things to stdout from
C     within runmd.
C     FMNDIM is initialized in a blockdata statement
C     in runmd.f
C
      INTEGER MAXFMT
      PARAMETER (MAXFMT = 20)

      CHARACTER FMNDIM*(MAXFMT)
      COMMON /RUNNY/FMNDIM
COMMEND
