module io_binary
  implicit none
contains


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

subroutine write_dbin(array, n, filename)
  real(8), intent(in) :: array(n)
  integer(8), intent(in) :: n
  character(len=*) :: filename
  integer(8) :: unit, iost

  ! Open the file
  open(newunit=unit, file=filename, form='unformatted', &
       access='stream', status='replace', action='write', iostat=iost)

  if (iost /= 0) then
    print *, 'Error opening file for writing.'
    return
  end if

  ! Write the array
  write(unit) array

  ! Close the file
  close(unit)
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
