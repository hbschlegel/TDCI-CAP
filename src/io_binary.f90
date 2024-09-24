module io_binary

  use variables_global !: Really only need iout=42
  use iso_fortran_env, only: system

  implicit none
contains


!: Remove directory if it exists, and the re-make the directory
subroutine cleanup_directory(dir_name)
  character(len=*), intent(in) :: dir_name
  integer :: stat
  character(len=255) :: cmd
  

  ! Check if the subdirectory exists
  inquire(file=trim(dir_name)//"/.", exist=stat)
  if (stat == 1) then
    ! Subdirectory exists, delete it recursively
    write(iout,*) "Overwriting subdirectory '"//trim(dir_name)//"'"
    cmd = "rm -rf "//trim(dir_name)
    call execute_command_line(cmd)
  end if

  ! Create the subdirectory
  cmd = "mkdir -p "//trim(dir_name)
  call execute_command_line(cmd)

  write(iout,*) "end of cleanup_directory()"

end subroutine cleanup_directory


!: Unit test for write_dbin and write_dbin
subroutine io_bin_test
  implicit none
  integer(8), parameter :: n = 45
  real(8) :: write_array(n), read_array(n)
  character(len=100) :: filename
  integer(8) :: i, ios

  ! Initialize an array with some test values
  do i = 1, n
     write_array(i) = i * 2.0d0
  end do

  ! Set the filename
  filename = 'test_array.bin'

  ! Write the array to a binary file
  call write_dbin(write_array, n, filename)

  ! Read the array from the binary file
  call read_dbin(read_array, n, filename, ios)

  ! Check for errors in reading the file
  if (ios /= 0) then
    print *, 'Error reading file. IOSTAT = ', ios
    stop
  end if

  ! Compare the arrays
  do i = 1, n
    if (write_array(i) /= read_array(i)) then
      print *, 'Test failed: Arrays do not match.'
      stop
    end if
  end do

  print *, 'Test passed: Arrays match.'

  ! Test the edge case: file does not exist
  call read_dbin(read_array, n, 'nonexistent.bin', ios)
  if (ios == 0) then
    print *, 'Edge case test failed: No error for non-existent file.'
  else
    print *, 'Edge case test passed: Error caught for non-existent file.'
  end if

end subroutine io_bin_test



subroutine write_dbin_safe(array, n, filename, funit)
  real(8), intent(in) :: array(n)
  integer(8), intent(in) :: n
  character(len=*) :: filename
  integer(8), intent(in), optional :: funit

  logical :: verify_correct, verify_exists
  integer(8) :: verify_attempts, ios
  integer(8), parameter :: max_verify_attempts = 10
  real(8) :: wait_time = 3.0
  real(8) :: verify_array(n)


  !: Maybe putting an extra critical region here will prevent our thread unsafety
  !$OMP CRITICAL (write_dbin_safe_lock)
  ! Call the original subroutine to attempt writing the file
  if (present(funit)) then
    call write_dbin(array, n, filename, funit)
  else
    call write_dbin(array, n, filename)
  end if

  ! Now, let's check if the file has been written
  verify_correct = .false.
  verify_attempts = 0
  do while (.not. verify_correct .and. verify_attempts < max_verify_attempts)
    verify_exists = .false.
    inquire(file=filename, exist=verify_exists)
    if (verify_exists) then
      ! Check that the first element is correct
      verify_array = 0.d0
      call read_dbin(verify_array, n, filename, ios)
      if (ios .ne. 0) then
        write(iout,*) "Error reading for verify: ", filename
      else 
        if (verify_array(1) .eq. array(1)) then
          verify_correct = .true.
        else
          write(iout,*) "FIRST ELEMENT OF ", filename, "IS NOT CORRECT: ", array(1), " ", verify_array(1)
        end if
      end if
    end if
    if (.not. verify_correct) then
      ! if something messed up, try again!
      call sleep(wait_time)
      if (present(funit)) then
        call write_dbin(array, n, filename, funit)
      else
        call write_dbin(array, n, filename)
      end if
      verify_attempts = verify_attempts + 1
    else
      exit ! File exists, let's break out of the loop
    end if
  end do

  if (.not. verify_correct) then
    ! If after all attempts the file still isn't there, let's inform the user
    write(iout, *) "ERROR: ", filename, " could not be verified after ", max_verify_attempts, " attempts :("
  else
    if (verify_attempts > 1) then
      write(iout, *) "SUCCESS: ", filename, " has been successfully written and verified! (", verify_attempts, " attempts)"
    end if
  end if
  !$OMP END CRITICAL (write_dbin_safe_lock)
end subroutine write_dbin_safe



subroutine write_dbin(array, n, filename, funit)
  real(8), intent(in) :: array(n)
  integer(8), intent(in) :: n
  character(len=*) :: filename
  integer(8), intent(in), optional :: funit
  integer(8) :: localUnit, iost
  logical    :: fileisopen
  
  integer(8) :: attempts
  integer(8), parameter :: max_attempts = 12
  real(8), parameter :: sleep_time = 2.0 ! in seconds

  attempts = 1
  fileisopen = .false.
  if (present(funit)) then
    inquire(unit=funit, opened=fileisopen)
    !: If the file unit we were given is already open, close it!
    do while (fileisopen)
      flush(funit)
      call sleep(sleep_time)
      close(funit)
      call sleep(sleep_time)
      inquire(unit=funit, opened=fileisopen)
    end do
  end if
  
  do while (.not. fileisopen .and. attempts < max_attempts)
    if (present(funit)) then
      !write(iout, *) "Provided unit: ", funit, " to write file ", filename ; flush(iout)
      localUnit = funit
      open(unit=localUnit, file=filename, form='unformatted', &
           access='stream', status='replace', action='write', iostat=iost)
    else
      !write(iout, *) "Getting new unit to write file ", filename ; flush(iout)
      open(newunit=localUnit, file=filename, form='unformatted', &
           access='stream', status='replace', action='write', iostat=iost)
    end if
    if (iost /= 0) then
      write(iout, *) "Error opening file: ", iost
      attempts = attempts+1
      call sleep(sleep_time)
      cycle ! skips to next loop, like python continue
    end if

    inquire(unit=localUnit, opened=fileisopen)

    if (.not. fileisopen) then
      write(iout, *) "Attempt ", attempts, ": ", filename, " DIDNT OPEN!!!"
      attempts = attempts + 1
      call sleep(sleep_time) ! wait a lil bit before trying again
    end if
  end do

  if (fileisopen) then
    if (attempts > 1) then
      write(iout, *) "Attempt ", attempts, ": ", filename, " properly opened!"
    end if
    ! Write the array
    write(localUnit, iostat=iost) array
    if (iost /= 0) then
      write(iout, *) "Error writing file: ", iost
    end if
    flush(localUnit)
    ! Close the file
    close(localUnit)
  else
    write(iout, *) "ERROR: Could not open ", filename, " after ", max_attempts, " attempts!"
  end if

end subroutine write_dbin


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


end module io_binary
