module sorting_module
  implicit none
contains

  subroutine quicksort_descending(A, B, key)
    real(8), intent(in)  :: A(:)
    real(8), intent(out) :: B(size(A))
    integer(8), intent(out) :: key(size(A))
    integer(8) :: i, lo, hi

    B = A
    do i = 1,size(A)
      key(i) = i
    end do

    lo = 1
    hi = size(A)
    call quicksort_helper(B, key, lo, hi)
  end subroutine quicksort_descending
  
  recursive subroutine quicksort_helper(arr, key, lo, hi)
    real(8), intent(inout) :: arr(:)
    integer(8), intent(inout) :: key(:)
    integer(8), intent(in) :: lo, hi
    integer(8) :: p
    if (lo < hi) then
      p = partition(arr, key, lo, hi)
      call quicksort_helper(arr, key, lo, p - 1)
      call quicksort_helper(arr, key, p + 1, hi)
    end if
  end subroutine quicksort_helper

  integer function partition(arr, key, lo, hi) result(p)
    real(8), intent(inout) :: arr(:)
    integer(8), intent(inout) :: key(:)
    integer(8), intent(in) :: lo, hi
    real(8) :: temp_val
    integer(8) :: temp_key, i, j

    i = lo
    do j = lo, hi - 1
      if (arr(j) > arr(hi)) then
        ! Swap array values
        temp_val = arr(i)
        arr(i) = arr(j)
        arr(j) = temp_val
        ! Swap indices in sortkey
        temp_key = key(i)
        key(i) = key(j)
        key(j) = temp_key

        i = i + 1
      end if
    end do
    temp_val = arr(i)
    arr(i) = arr(hi)
    arr(hi) = temp_val
    temp_key = key(i)
    key(i) = key(hi)
    key(hi) = temp_key
    p = i
  end function partition

  subroutine sort_eigenpairs_desc(eigenvalues, eigenvectors_1D, n)
    implicit none
    integer(8), intent(in) :: n
    real(8), dimension(n), intent(inout) :: eigenvalues
    real(8), dimension(n*n), intent(inout) :: eigenvectors_1D
    integer(8) :: i, j, k
    real(8) :: temp_val

    ! Bubble sort
    do i = 1, n-1
      do j = 1, n-i
        if (eigenvalues(j) < eigenvalues(j+1)) then
          ! Swap eigenvalues
          temp_val = eigenvalues(j)
          eigenvalues(j) = eigenvalues(j+1)
          eigenvalues(j+1) = temp_val

          ! Swap corresponding eigenvectors (columns)
          do k = 1, n
            temp_val = eigenvectors_1D((j-1)*n + k)
            eigenvectors_1D((j-1)*n + k) = eigenvectors_1D((j)*n + k)
            eigenvectors_1D((j)*n + k) = temp_val
          end do
        end if
      end do
    end do

  end subroutine sort_eigenpairs_desc





















end module sorting_module

