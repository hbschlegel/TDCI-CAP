

!:
!: Construct a formated checkpoint file from a template formatted chekpoint file 
!: and the MO density data at a selected step in an MO_densty-ei-dj.dat file
!: The ouput formatted checkpoint file contains the total density
!:
!: Usage:
!:   density2fchk infile MOdensityfile nstep diff iscale
!:
!: Command line input:
!:    infile -  full name of the template formatted checkpoint file (with ".fchk")
!:    MOdensityfile - full name of the ion coefficient file (e.g. MO_density-e1-d1.dat)
!:    nstep - step number to be used for the constructiion of the NTOs and hole density
!:    diff = 0 - calculate the density
!:         = 1 - calculate the density minus the density of the field-free neutral
!:         = 2 - calculate the density minus the density in the template fchk file
!:         = 3 - calculate the density minus the density at time step 1
!:    iscale = 0 - do not scale NTOs or hole density
!:           = 1 - scale by the ratio of the rate to the maximum rate
!:
!: Output
!:    infile-1.fchk - formatted checkpoint file with total density
!:    total SCF density replaced bt the density or density difference
!:


module density2fchk

contains

 
subroutine read_dbin(array, n, filename, ios)
  real(8), intent(out) :: array(n)
  integer(8), intent(in) :: n
  character(len=*) :: filename
  integer(8), intent(out) :: ios
  integer(8) :: unit

  ! Open the file
  open(newunit=unit, file=filename, form='unformatted', &
       access='stream', status='old', action='read', iostat=ios)

  if (ios /= 0) then
    print *, 'Error opening file for reading.'
    return
  end if

  ! Read the array
  read(unit, iostat=ios) array

  if (ios /= 0) then
    print *, 'Error reading from file.'
    return
  end if

  ! Close the file
  close(unit)
end subroutine read_dbin


function remove_extension(filename) result(basefilename)
    implicit none
    character(80), intent(in) :: filename
    character(255) :: basefilename
    integer :: pos

    pos = scan(filename, ".", back=.true.)

    if (pos == 0) then
        basefilename = filename
    else
        basefilename = filename(:pos-1)
    end if

end function remove_extension




end module density2fchk

      program main
      use density2fchk

      implicit none
      character(80) line
      character(80) infile,densityfile,outfile
      integer(8) funit1,funit2,funit3
      integer(8) nea,neb,nbasis,nbsuse,ntt,nocc,norb,nstep,nend,diff,iscale
      integer(8) i,j,ij,k,kk,n
      real(8) time,x,rate,maxrate,factor
      real(8), allocatable ::  cmo(:),scratch(:),density(:)


      character(255) :: filename
      integer(8) :: ios

      write(*,*) "Start density2fchk"

      funit1 = 10
      funit2 = 11
      funit3 = 12
!:
!: get command line arguments
!: template fchk file, ion_coeffiient file, step
!:
      call get_command_argument(1,infile) 
      call get_command_argument(2,densityfile) 
      call get_command_argument(3,line) 
        read(line,*) nstep
      call get_command_argument(4,line) 
        read(line,*) diff
      call get_command_argument(5,line) 
        read(line,*) iscale
!:
!: open fchk file and read info
!:
      write(*,"('  Input file:                ',A)") trim(infile)
      open(funit1,file=trim(infile))
      do 
        read(funit1,"(A80)",ERR=100) line
        if(trim(line(1:25))==trim("Number of alpha electrons")) then
          read(line,"(51x,i10)") nea
        endif
        if(trim(line(1:24))==trim("Number of beta electrons"))  then
          read(line,"(51x,i10)") neb
        endif
        if(trim(line(1:25))==trim("Number of basis functions"))  then
          read(line,"(51x,i10)") nbasis
          ntt = nbasis*(nbasis+1)/2
        endif
        if(trim(line(1:31))==trim("Number of independent functions")) then
          read(line,"(51x,i10)") nbsuse
!:          write(*,*) nea,nbasis,ntt,nbsuse
!:
!: allocate arrays
!:
          allocate( cmo(nbasis*nbsuse) )
          allocate( scratch(nbasis*nbasis) )
          allocate( density(nbasis*nbasis) )
        endif
        if(trim(line(1:21))==trim("Alpha MO coefficients")) then
          read(funit1,*) (cmo(j),j=1,nbasis*nbsuse)
          exit
        endif
      enddo
  100 continue
      close(funit1)
!:
!: open density file 
!:
      write(*,"('  Density file:              ',A)") trim(densityfile)
      open(funit2,file=trim(densityfile))
      read(funit2,"(7x,I10)") nend
!:      write(*,*) " nend ",nend
      read(funit2,"(11x,I10)") nocc
      read(funit2,"(10x,I10)") nocc
!:      write(*,*) " nocc ",nocc
      read(funit2,"(15x,I10)") norb
!:      write(*,*) " norb ",norb
!:
!: skip records in MO density file to get to nstep
!:
      density = 0.d0
      if(nstep.gt.1) then
        do n = 1, nstep-1
          if(n.eq.1) write(*,*) "1st loop "
          read(funit2,"(f13.9)") time
          do
            read(funit2,"(i5,i5,1x,f13.10)") i,j,x
            if(i.eq.0 .and. j.eq.0) go to 150
            if( diff.eq.3 .and. n.eq.1 ) then
!:              if(abs(x).gt.1.d-4) write(*,"(2I5,1x,f16.10)") i,j,x
              if(i.eq.j) then
                density((i-1)*norb+i) = -x
              else
                density((i-1)*norb+j) = -0.5d0 * x
                density((j-1)*norb+i) = -0.5d0 * x
              endif
            endif
          enddo
  150     continue
        enddo
      endif
!:
!: read elements of the density matrix
!:
        read(funit2,"(f13.9,f16.10)") time,rate
        do
          read(funit2,"(i5,i5,1x,f13.10)") i,j,x
          if(i.eq.0 .and. j.eq.0) go to 160
!:          if(abs(x).gt.1.d-4) write(*,"(2I5,1x,f16.10)") i,j,x
          if(i.eq.j) then
            density((i-1)*norb+i) = density((i-1)*norb+i) + x
!:            write(*,*) i,x
          else
            density((i-1)*norb+j) = density((i-1)*norb+j) + 0.5d0 * x
            density((j-1)*norb+i) = density((j-1)*norb+i) + 0.5d0 * x
          endif
!:          if(abs(x).gt.1.d-5) write(*,"(2I5,1x,f16.10)") i,j,density((i-1)*norb+j)
        enddo
  160   continue
        if( nstep.eq.1 .and. diff.eq.3 ) density = 0.d0
        close(funit2)

!:
!: read density matrix from bin file
!:

!: Keep the previous code for rate, but overwrite density
density = 0.d0

write( filename, '(A,A,I0,A)') trim(remove_extension(densityfile)), trim("."), &
                               nstep, "00.bin"

write(*,"(A,A)") "Reading file: ", filename

call read_dbin( density, norb*norb, trim(filename) , ios)

write(*,*) "ios: ", ios

write(*,"(A, E15.8, E15.8)") "density bin: ", density(1), density(2)
write(*,"(A,E15.8)") "sum of density: ", sum(density)


!: Let's make sure the dimension is actually norb by checking the diagonal
write(*,*) "density(1) = ", density(1)
write(*,*) "density((2-1)*norb+2) = ", density((2-1)*norb+2)
write(*,*) "density((3-1)*norb+3) = ", density((3-1)*norb+3)



!write(*,"(5F12.8)")  (density(i), i=1,25)

write(*,*) "first 7x7 block of matrix: "
do i=1,7
  do j=1,7
    write(*, "(E10.3,1X)", advance='no') density((i-1)*norb+j)
  end do
  write(*,*) " "
end do




!:
!: subtract density of the neutral if requested by diff = 1
!
        if(diff.eq.1) then
          do i = 1,nocc
            density((i-1)*norb+i) = density((i-1)*norb+i) - 2.d0
          enddo
        endif
!:
!: if scaling is requested, read maxrate from file maxrate.dat and calculate scale factor
!:
      if(iscale.eq.0) then
        factor = 1.d0
      else
        open(funit2,file=trim("maxrate.dat"))
        read(funit2,"(f16.10)") maxrate
        close(funit2)
        !factor = abs(rate) / maxrate
        factor = 1.0/maxrate
      endif
      write(*,"(A, E15.8)") "Scaling factor: ", factor

      !do i = 1,norb*norb
      !  if (density(i) < -1E-8) then
      !    write(*,"(A, I5, E15.8)") "MO density < 0 :  ",i, density(i)
      !  end if 
      !  density(i) =  abs(factor*density(i))
      !end do

!:
!: transform density from MO basis to AO basis
!:
        scratch = 0.d0
        do i = 1,nbasis
          do j = 1,norb
            do k = 1,norb
              kk = (k+nea-nocc-1)*nbasis
              scratch((i-1)*norb+j) = scratch((i-1)*norb+j) &
                + cmo(kk+i) * density((k-1)*norb+j)
            enddo
          enddo
        enddo
        density = 0.d0
        do i = 1,nbasis
          do j = 1,nbasis
            do k = 1,norb
              kk = (k+nea-nocc-1)*nbasis
              density((i-1)*nbasis+j) = density((i-1)*nbasis+j) &
                + scratch((i-1)*norb+k) * cmo(kk+j)
            enddo
          enddo
        enddo

  !do i=1,nbasis*nbasis
  !  if (density(i) < -1E-8) then
  !    write(*,"(A, I8, E15.8)") "AO density < 0 :  ",i, density(i)
  !  end if
  !end do
!:
!: copy template fchk file to output fchk file and insert density
!:
        open(funit1,file=trim(infile))
!:
!: build the output file name
!:
        i = INDEX(infile,".fch",.true.) - 1
        outfile = infile(1:i)//"-1.fchk"
        write(*,"('  Output file for step',I4,':  ',A)") &
          nstep,trim(outfile)
        open(funit3,file=trim(outfile))
!:
!: read template formatted checkpoint file
!:
        do
          read(funit1,"(A80)",ERR=200) line
          if(trim(line(1:16))==trim("Gaussian Version")) then
            write(funit3,"(A80)") line
            read(funit1,"(A80)",ERR=200) line
            write(funit3,"(A80)") line
            exit
          endif
!:
!: write density in place of the total SCF density
!:
          if(trim(line(1:17))==trim("Total SCF Density")) then
            write(funit3,"(A80)") line
            read(funit1,*) (scratch(i),i=1,ntt)
!:
!: subtract density from the template fchk file if requested by diff = 2
!:
            if(diff.eq.2) then
              ij = 0
              do i = 1,nbasis
                do j = 1,i
                  ij = ij + 1
                  density((i-1)*nbasis+j) = density((i-1)*nbasis+j) - scratch(ij)
                end do
              end do
            end if
            write(funit3,"(5(1x,E15.8))") &
              ((density((i-1)*nbasis+j),j=1,i),i=1,nbasis)
          else
            write(funit3,"(A80)") line
          endif
        enddo
  200   continue
        close(funit1)
        close(funit3)
      
      write(*,*) "End density2fchk"
 
      end program main
