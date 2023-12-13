module sorting_module
  implicit none
contains

  subroutine quicksort_descending(A, B)
    real(8), intent(in)  :: A(:)
    real(8), intent(out) :: B(size(A))
    B = A
    call quicksort_helper(B, 1, size(A))
  end subroutine quicksort_descending
  
  recursive subroutine quicksort_helper(arr, lo, hi)
    real(8), intent(inout) :: arr(:)
    integer, intent(in) :: lo, hi
    integer :: p
    if (lo < hi) then
      p = partition(arr, lo, hi)
      call quicksort_helper(arr, lo, p - 1)
      call quicksort_helper(arr, p + 1, hi)
    end if
  end subroutine quicksort_helper

  integer function partition(arr, lo, hi) result(p)
    real(8), intent(inout) :: arr(:)
    integer, intent(in) :: lo, hi
    real(8) :: pivot, temp
    integer :: i, j

    pivot = arr(hi)
    i = lo
    do j = lo, hi - 1
      if (arr(j) > pivot) then
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        i = i + 1
      end if
    end do
    temp = arr(i)
    arr(i) = arr(hi)
    arr(hi) = temp
    p = i
  end function partition

  subroutine sort_eigenpairs_desc(eigenvalues, eigenvectors_1D, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: eigenvalues
    real(8), dimension(n*n), intent(inout) :: eigenvectors_1D
    integer :: i, j, k
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

