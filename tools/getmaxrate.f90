      program main
!:
!: Find the maximum rate for a set of ION_COEFF-e1-di.dat files
!: and write the maximum rate to maxrate.dat
!:
!: Usage:
!:   getmaxrate dirstart dirend nstart nend
!:
!: Command line input:
!:    dirstart, dirend - starting and ending directions for search for the maximum rate
!:    nstart, nend - starting and ending time steps for search for the maximum rate
!:
!: Output
!:    maxrate.dat - file with the maximum rate
!:
      implicit none
      character(4) dirstr
      character(80) line
      character(80) ionfile,outfile
      integer(8) funit1,funit2,i,j,idata,idir,dirstart,dirend
      integer(8) istep,nstart,nend,ntimes,nstates
      real(8) rate,maxrate
      real(8), allocatable ::  scratch(:)

      funit1 = 10
      funit2 = 11
!:
!: get command line arguments
!:
      call get_command_argument(1,line) 
        read(line,*) dirstart
      call get_command_argument(2,line) 
        read(line,*) dirend
      call get_command_argument(3,line) 
        read(line,*) nstart
      call get_command_argument(4,line) 
        read(line,*) nend
      write(*,*) dirstart,dirend,nstart,nend
!:
!: output file
!:
      outfile = "maxrate.dat"
      open(funit1,file=trim(outfile))
!:
!: loop over directions
!:
      do idir = dirstart,dirend
!:
!: read ion_coefficient file
!:
        write( dirstr, '(i0)' )  idir
        ionfile = 'ION_COEFF-e1-d'//trim(dirstr)//'.dat'
          write(*,*) ionfile
        open(funit2,file=trim(ionfile))
        read(funit2,"(A80)") line
        read(funit2,"(7x,I10)") ntimes
        write(*,*) " ntimes ",ntimes
        read(funit2,"(9x,I11)") nstates
        write(*,*) " nstates ",nstates
        if(idir .eq. 1) allocate( scratch(2*nstates) )
!:
!: loop over time steps to find the maximum rate
!:
        maxrate = 0.d0
        do istep = 1, nend
          read(funit2,"(i5,10x,f16.10)") idata,rate
          write(*,*) istep,idata,rate
          if(istep.ge.nstart .and. abs(rate).gt.maxrate) maxrate = abs(rate)
!:          write(*,*) " maxrate ",maxrate
          do i = 1, 2*nstates
            read(funit2,*) (scratch(j),j=1,2*nstates)
          enddo
        enddo
      enddo
      write(funit1,"(f16.10)") maxrate
      close(funit1)
      close(funit2)
 
      end program main
