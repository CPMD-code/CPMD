! Purpose: utilities for interface with external programs
! Author: Pablo Baudin
! Date: May 2018
module interface_utils

   use adat,                            only: elem
   use fileopenmod,                     only: fo_info
   use ions,                            only: ions0,ions1
   use machine,                         only: m_sleep, m_system
   use mp_interface,                    only: mp_bcast, mp_sync
   use parac,                           only: paral,parai
   use system,                          only: iatpt,maxsys
   use zeroing_utils,                   only: zeroing
   use kinds,                           only: real_8
   use timer,                           only: tiset, tihalt
   use cnst,                            only: fbohr
   use error_handling,                  only: stopgm

   implicit none

   private 

   public :: get_external_forces

contains   

   ! Purpose: Call external program to get ionic forces from current position
   ! Author: Pablo Baudin
   ! Date: May 2018 
   subroutine get_external_forces(prog, positions, forces, path, args)
      implicit none
      !> Name of the program or script to call
      character(len=*), intent(in) :: prog
      !> position of the ions to be used for the calculation of the forces
      real(real_8), intent(in) ::  positions(3,maxsys%nax,*)
      !> forces produced by the external program
      real(real_8), intent(out) ::  forces(:,:,:)
      !> path of the script and arguments to pass to the script
      character(len=*), intent(in), optional :: path, args

      character(len=100) :: mypath, myargs, title
      character(len=3) :: symbol
      logical :: exists
      integer :: funit, iat, ia, is, l, mat, msglen

      mypath = fo_info%fpath(fo_info%iapath:fo_info%iepath)
      if (present(path)) mypath = path
      myargs = ''
      if (present(args)) myargs = ' '//args
      funit = 211


      ! write the ionic positions to file to be read by program
      if (paral%io_parent) then
         open(funit,file='geometry.xyz',status='unknown')
         write(funit,*) ions1%nat 
         write(funit,*) 'new GEOMETRY from CPMD'
         do iat=1,ions1%nat
            ia=iatpt(1,iat)
            is=iatpt(2,iat)
            !look for atomic label
            write(funit,'(A3,3F20.12)') elem%el(ions0%iatyp(is)), (positions(l,ia,is)/fbohr,l=1,3)
         enddo
         close(funit)
      endif
      call mp_sync(parai%allgrp)


      ! check if external program exists and call it
      if(paral%io_parent) then
         inquire(file=trim(mypath)//prog,exist=exists)
         if (exists) then
            ! call program
            call m_system(trim(mypath)//prog//myargs)
         else
            call stopgm('call_external','External program '//trim(mypath)//prog//' not present!',&
            __LINE__,__FILE__)
         end if
      end if


      ! wait until program is done (check for existence of file: forces.xyz)
      call interface_f_wait('forces.xyz', 'EXT.error', mypath)

      ! read forces
      call zeroing(forces)
      if (paral%io_parent) then
         open(funit,file=trim(mypath)//'forces.xyz',status='unknown')
         read(funit,*) mat 
         read(funit,'(a)') title
         write(6,*) title
         do iat=1,mat
            ia=iatpt(1,iat)
            is=iatpt(2,iat)
            read(funit,*) symbol,(forces(l,ia,is),l=1,3)
         enddo
         close(funit)
      endif
      call mp_sync(parai%allgrp)
      msglen = 3*maxsys%nax*maxsys%nsx
      call mp_bcast(forces,msglen,parai%io_source,parai%allgrp)

   end subroutine get_external_forces

   ! Purpose: wait until a given file has been produced
   ! Author: Pablo Baudin
   ! Date: May 2018
   subroutine interface_f_wait(filename, error_file, path)
      implicit none
      !> name of the file to be produced
      character(len=*), intent(in) :: filename, error_file
      !> path of file
      character(len=*), intent(in), optional :: path

      logical ::  test, error
      integer ::  isub
      character(*), parameter :: procedureN = 'interface_f_wait'
      character(len=100) :: mypath

      mypath = fo_info%fpath(fo_info%iapath:fo_info%iepath)
      if (present(path)) mypath = path

      call tiset(procedureN,isub)
      if (paral%io_parent) then
         do while (.not. test)

            ! wait until required file is produced   
            inquire(file=trim(mypath)//filename,exist=test)

            ! wait 1 sec.
            call m_sleep(1)

            ! check for error file
            inquire(file=trim(mypath)//error_file,exist=error)
            if (error) call stopgm(proceduren,'External program error see '//error_file//' file',&
               __LINE__,__FILE__)

         end do
      end if
      call tihalt(procedureN,isub)

   end subroutine interface_f_wait

end module interface_utils
