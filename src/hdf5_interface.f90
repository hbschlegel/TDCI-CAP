module hdf5_writer

use, intrinsic :: iso_fortran_env, only: int64, real64
use hdf5

implicit none

private

type, public :: HDF5Writer
  integer(HID_T) :: file_id = -1
  character(len=:), allocatable :: filename
  contains
    procedure, public :: open_file        => h5w_open_file
    procedure, public :: close_file       => h5w_close_file
    procedure, public :: create_group     => h5w_create_group
    procedure, public :: close_group      => h5w_close_group
    procedure, public :: write_attribute  => h5w_write_attribute
    procedure, public :: create_dataset_real    => h5w_create_dataset_real
    procedure, public :: create_dataset_int     => h5w_create_dataset_int
    procedure, public :: create_dataset_complex => h5w_create_dataset_complex
    procedure, public :: append_real           => h5w_append_real
    procedure, public :: append_int            => h5w_append_int
    procedure, public :: append_complex        => h5w_append_complex
end type HDF5Writer

contains

!------------------------------------------------------------------------------
! 1) Open an HDF5 file (create or overwrite) using modern Fortran bindings.
!------------------------------------------------------------------------------
subroutine h5w_open_file(this, fname)
   class(HDF5Writer), intent(inout) :: this
   character(len=*),  intent(in)    :: fname

   integer :: status

   this%filename = fname
   ! The modern routine for creating a file is h5fcreate_f
   call h5fcreate_f(fname, H5F_ACC_TRUNC_F, this%file_id, status)
   if (status /= 0) then
      write(*,*) "Error creating HDF5 file: ", trim(fname)
      stop
   end if
end subroutine h5w_open_file


!------------------------------------------------------------------------------
! 2) Close the HDF5 file
!------------------------------------------------------------------------------
subroutine h5w_close_file(this)
   class(HDF5Writer), intent(inout) :: this
   integer :: status

   if (this%file_id /= -1) then
      call h5fclose_f(this%file_id, status)
      this%file_id = -1
      write(*,*) "Closed HDF5 file:", trim(this%filename)
   end if
end subroutine h5w_close_file


!------------------------------------------------------------------------------
! 3) Create (or open) a group, e.g. "/direction_i/field_j".
!    group_id is returned for further actions on that group.
!------------------------------------------------------------------------------
subroutine h5w_create_group(this, group_path, group_id)
   class(HDF5Writer), intent(inout) :: this
   character(len=*),   intent(in)  :: group_path
   integer(HID_T),     intent(out) :: group_id

   integer :: status

   call h5gcreate_f(this%file_id, group_path, group_id, status)
   if (status /= 0) then
      ! Possibly the group already exists; try opening it
      call h5gopen_f(this%file_id, group_path, group_id, status)
      if (status /= 0) then
         write(*,*) "Error creating/opening group:", trim(group_path)
         stop
      end if
   end if
end subroutine h5w_create_group


subroutine h5w_close_group(this, group_id)
   class(HDF5Writer), intent(inout) :: this
   integer(HID_T),     intent(in) :: group_id

   integer :: status

   call h5gclose_f(group_id, status)
   if (status /= 0) then
     write(*,*) "Error closing group:", group_id
   end if
end subroutine h5w_close_group


!------------------------------------------------------------------------------
! 4) Write an attribute (e.g. a small real array) to a group
!------------------------------------------------------------------------------
subroutine h5w_write_attribute(this, group_id, attr_name, data)
   class(HDF5Writer), intent(in)    :: this
   integer(HID_T),     intent(in)   :: group_id
   character(len=*),   intent(in)   :: attr_name
   real(8),       intent(in)   :: data(:)

   integer :: status
   integer(HID_T) :: attr_id, aspace_id
   integer(HSIZE_T), dimension(1) :: dims

   !real(H5T_IEEE_F64LE), allocatable :: tempdata(:)

   dims(1) = size(data)
   !allocate(tempdata(size(data)))

   call h5screate_simple_f(1, dims, aspace_id, status)
   if (status /= 0) stop "Error creating attribute dataspace."

   call h5acreate_f(group_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, status)
   if (status /= 0) stop "Error creating attribute."

   call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, status)
   if (status /= 0) stop "Error writing attribute."

   call h5aclose_f(attr_id, status)
   call h5sclose_f(aspace_id, status)
end subroutine h5w_write_attribute


!------------------------------------------------------------------------------
! 5) Create an UNLIMITED dataset for real(8), shape [time, dim2, dim3, ...].
!------------------------------------------------------------------------------
subroutine h5w_create_dataset_real(this, group_id, dset_name, dims, dset_id)
   class(HDF5Writer),                intent(inout) :: this
   integer(HID_T),                   intent(in)    :: group_id
   character(len=*),                 intent(in)    :: dset_name
   integer(HSIZE_T), dimension(:),   intent(in)    :: dims
   integer(HID_T),                   intent(out)   :: dset_id

   integer :: status, rank, i
   integer(HID_T) :: space_id, dcpl_id
   integer(HSIZE_T), allocatable :: current_dims(:), maxdims(:), chunk_dims(:)

   rank = size(dims) + 1   ! +1 for the unlimited (time) dimension
   allocate(current_dims(rank), maxdims(rank), chunk_dims(rank))

   ! First dimension is "time", start at 0, unlimited
   current_dims(1) = 0
   maxdims(1)      = H5S_UNLIMITED_F

   ! The rest are fixed
   do i = 2, rank
      current_dims(i) = dims(i-1)
      maxdims(i)      = dims(i-1)
   end do

   ! Create the dataspace
   call h5screate_simple_f(rank, current_dims, space_id, status, maxdims=maxdims)
   if (status /= 0) stop "Error creating dataspace for real dataset."

   ! Create property list for chunking
   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, status)
   if (status /= 0) stop "Error creating DCPL."

   ! Example chunk dims: chunk along time in blocks of 8
   chunk_dims(1) = 8
   do i = 2, rank
      chunk_dims(i) = current_dims(i)
   end do

   call h5pset_chunk_f(dcpl_id, rank, chunk_dims, status)
   if (status /= 0) stop "Error setting chunk dims."

   ! Create the dataset
   call h5dcreate_f(                          &
        loc_id   = group_id,                  &
        name     = dset_name,                 &
        type_id  = H5T_NATIVE_DOUBLE,            &
        space_id = space_id,                  &
        dset_id  = dset_id,                   &
        hdferr   = status,                    &
        lcpl_id  = H5P_DEFAULT_F,             &
        dcpl_id  = dcpl_id,                   &
        dapl_id  = H5P_DEFAULT_F)
   if (status /= 0) stop "Error creating real dataset."

   call h5pclose_f(dcpl_id, status)
   call h5sclose_f(space_id, status)
end subroutine h5w_create_dataset_real


!------------------------------------------------------------------------------
! 6) Create an unlimited dataset for integer(8).
!------------------------------------------------------------------------------
subroutine h5w_create_dataset_int(this, group_id, dset_name, dims, dset_id)
   class(HDF5Writer),                intent(inout) :: this
   integer(HID_T),                   intent(in)    :: group_id
   character(len=*),                 intent(in)    :: dset_name
   integer(HSIZE_T), dimension(:),   intent(in)    :: dims
   integer(HID_T),                   intent(out)   :: dset_id

   integer :: status, rank, i
   integer(HID_T) :: space_id, dcpl_id
   integer(HSIZE_T), allocatable :: current_dims(:), maxdims(:), chunk_dims(:)

   rank = size(dims) + 1
   allocate(current_dims(rank), maxdims(rank), chunk_dims(rank))

   current_dims(1) = 0
   maxdims(1)      = H5S_UNLIMITED_F

   do i = 2, rank
      current_dims(i) = dims(i-1)
      maxdims(i)      = dims(i-1)
   end do

   call h5screate_simple_f(rank, current_dims, space_id, status, maxdims=maxdims)
   if (status /= 0) stop "Error creating dataspace for int dataset."

   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, status)
   if (status /= 0) stop "Error creating DCPL."

   ! Example chunk dims
   chunk_dims(1) = 8
   do i = 2, rank
      chunk_dims(i) = current_dims(i)
   end do
   call h5pset_chunk_f(dcpl_id, rank, chunk_dims, status)

   ! For integer(8):
   call h5dcreate_f(                       &
        loc_id   = group_id,               &
        name     = dset_name,              &
        type_id  = H5T_STD_I64LE,          &
        space_id = space_id,               &
        dset_id  = dset_id,                &
        hdferr   = status,                 &
        lcpl_id  = H5P_DEFAULT_F,          &
        dcpl_id  = dcpl_id,                &
        dapl_id  = H5P_DEFAULT_F)
   if (status /= 0) stop "Error creating int dataset."

   call h5pclose_f(dcpl_id, status)
   call h5sclose_f(space_id, status)
end subroutine h5w_create_dataset_int


!------------------------------------------------------------------------------
! 7) Create an unlimited dataset for COMPLEX(8).
!    We store an extra dimension of size 2 for real+imag.
!------------------------------------------------------------------------------
subroutine h5w_create_dataset_complex(this, group_id, dset_name, dims, dset_id)
   class(HDF5Writer),                intent(inout) :: this
   integer(HID_T),                   intent(in)    :: group_id
   character(len=*),                 intent(in)    :: dset_name
   integer(HSIZE_T), dimension(:),   intent(in)    :: dims
   integer(HID_T),                   intent(out)   :: dset_id

   integer :: status, rank, i
   integer(HID_T) :: space_id, dcpl_id
   integer(HSIZE_T), allocatable :: current_dims(:), maxdims(:), chunk_dims(:)

   ! rank = user dims + 1 (time) + 1 (real/imag)
   rank = size(dims) + 2
   allocate(current_dims(rank), maxdims(rank), chunk_dims(rank))

   ! time dimension
   current_dims(1) = 0
   maxdims(1)      = H5S_UNLIMITED_F

   ! user dims
   do i = 2, rank-1
      current_dims(i) = dims(i-1)
      maxdims(i)      = dims(i-1)
   end do

   ! final dimension = 2
   current_dims(rank) = 2
   maxdims(rank)      = 2

   call h5screate_simple_f(rank, current_dims, space_id, status, maxdims=maxdims)
   if (status /= 0) stop "Error creating dataspace for complex dataset."

   call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, status)

   chunk_dims(1) = 8
   do i = 2, rank
      chunk_dims(i) = current_dims(i)
   end do
   call h5pset_chunk_f(dcpl_id, rank, chunk_dims, status)

   call h5dcreate_f(                       &
        loc_id   = group_id,               &
        name     = dset_name,              &
        type_id  = H5T_NATIVE_DOUBLE,         & ! store double in last dimension
        space_id = space_id,               &
        dset_id  = dset_id,                &
        hdferr   = status,                 &
        lcpl_id  = H5P_DEFAULT_F,          &
        dcpl_id  = dcpl_id,                &
        dapl_id  = H5P_DEFAULT_F)
   if (status /= 0) stop "Error creating complex dataset."

   call h5pclose_f(dcpl_id, status)
   call h5sclose_f(space_id, status)
end subroutine h5w_create_dataset_complex


!------------------------------------------------------------------------------
! 8) Append a 1D real(8) array to an unlimited dataset of shape [time, Nx].
!    Extends dimension 1 by 1, then writes data(:).
!------------------------------------------------------------------------------
subroutine h5w_append_real(this, dset_id, data)
   class(HDF5Writer), intent(inout) :: this
   integer(HID_T),    intent(in)    :: dset_id
   real(real64), dimension(:), intent(in) :: data

   integer :: status
   integer(HID_T) :: filespace_id, memspace_id
   integer(HSIZE_T), allocatable :: dset_dims(:), dset_dims_max(:), newsize(:), start(:), count(:)
   integer :: rank
   integer(HSIZE_T) :: tempdims(2)

   ! 1) Get current dataspace & rank
   call h5dget_space_f(dset_id, filespace_id, status)
   if (status /= 0) stop "Error getting file dataspace."

   call h5sget_simple_extent_ndims_f(filespace_id, rank, status)
   allocate(dset_dims(rank), dset_dims_max(rank), newsize(rank), start(rank), count(rank))

   call h5sget_simple_extent_dims_f(filespace_id, dset_dims, dset_dims_max, status)
   ! 2) Extend time dimension by 1
   newsize = dset_dims
   newsize(1) = newsize(1) + 1

   call h5dextend_f(dset_id, newsize, status)
   if (status /= 0) stop "Error extending dataset."

   ! Reacquire the (new) dataspace after extension
   call h5sclose_f(filespace_id, status)
   call h5dget_space_f(dset_id, filespace_id, status)
   if (status /= 0) stop "Error getting extended dataspace."

   ! 3) Select hyperslab for the new row
   start(1) = newsize(1) - 1
   count(1) = 1
   if (rank /= 2) stop "Append real: expected rank=2 dataset."

   start(2) = 0
   count(2) = size(data)

   call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, start, count, status)
   if (status /= 0) stop "Error selecting hyperslab."

   ! 4) Create memspace for writing [1 x Nx]
   tempdims = [1, size(data)]
   call h5screate_simple_f(2, tempdims, memspace_id, status)
   if (status /= 0) stop "Error creating memspace."

   ! 5) Write
   call h5dwrite_f(                             &
        dset_id       = dset_id,                &
        mem_type_id   = H5T_NATIVE_DOUBLE,       & ! Or H5T_NATIVE_DOUBLE
        buf           = data,                    &
        dims          = tempdims,         &
        hdferr        = status,                  &
        file_space_id = filespace_id,            &
        mem_space_id  = memspace_id)
   if (status /= 0) stop "Error writing real data."

   ! Cleanup
   call h5sclose_f(memspace_id, status)
   call h5sclose_f(filespace_id, status)
end subroutine h5w_append_real


!------------------------------------------------------------------------------
! 9) Append a 1D integer(int64) array to an unlimited dataset of shape [time, Nx].
!------------------------------------------------------------------------------
subroutine h5w_append_int(this, dset_id, data)
   class(HDF5Writer), intent(inout) :: this
   integer(HID_T),    intent(in)    :: dset_id
   integer(int64),    intent(in)    :: data(:)

   integer :: status
   integer(HID_T) :: filespace_id, memspace_id
   integer(HSIZE_T), allocatable :: dset_dims(:), dset_dims_max(:), newsize(:), start(:), count(:)
   integer :: rank
   integer(HSIZE_T) :: tempdims(2)

   ! 1) Get current dataspace & rank
   call h5dget_space_f(dset_id, filespace_id, status)
   if (status /= 0) stop "Error getting file dataspace."

   call h5sget_simple_extent_ndims_f(filespace_id, rank, status)
   allocate(dset_dims(rank), dset_dims_max(rank), newsize(rank), start(rank), count(rank))
   call h5sget_simple_extent_dims_f(filespace_id, dset_dims, dset_dims_max, status)

   ! 2) Extend along time dimension
   newsize = dset_dims
   newsize(1) = newsize(1) + 1
   call h5dextend_f(dset_id, newsize, status)
   if (status /= 0) stop "Error extending int dataset."

   call h5sclose_f(filespace_id, status)
   call h5dget_space_f(dset_id, filespace_id, status)

   ! 3) Select the new row
   start(1) = newsize(1) - 1
   count(1) = 1
   if (rank /= 2) stop "Append int: expected rank=2 dataset."

   start(2) = 0
   count(2) = size(data)

   call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, start, count, status)
   if (status /= 0) stop "Error selecting hyperslab."

   ! 4) Memspace
   tempdims = [1, size(data)]
   call h5screate_simple_f(2, tempdims, memspace_id, status)
   if (status /= 0) stop "Error creating memspace."

   ! 5) Write
   call h5dwrite_f(                         &
        dset_id       = dset_id,            &
        mem_type_id   = H5T_STD_I64LE,   &
        buf           = data,               &
        dims          = tempdims,    &
        hdferr        = status,             &
        file_space_id = filespace_id,       &
        mem_space_id  = memspace_id)
   if (status /= 0) stop "Error writing int data."

   call h5sclose_f(memspace_id, status)
   call h5sclose_f(filespace_id, status)
end subroutine h5w_append_int


!------------------------------------------------------------------------------
! 10) Append a 2D real array for COMPLEX(8) data to shape [time, Nx, 2].
!------------------------------------------------------------------------------
subroutine h5w_append_complex(this, dset_id, data)
   class(HDF5Writer),           intent(inout) :: this
   integer(HID_T),              intent(in)    :: dset_id
   real(real64), dimension(:,:),intent(in)    :: data
   ! data(:,1) = real part, data(:,2) = imag

   integer :: status
   integer(HID_T) :: filespace_id, memspace_id
   integer(HSIZE_T), allocatable :: dset_dims(:), dset_dims_max(:), newsize(:), start(:), count(:)
   integer :: rank
   integer(HSIZE_T) :: tempdims(3)

   ! 1) Get the current dataspace & rank
   call h5dget_space_f(dset_id, filespace_id, status)
   if (status /= 0) stop "Error getting dataspace (complex)."

   call h5sget_simple_extent_ndims_f(filespace_id, rank, status)
   allocate(dset_dims(rank), dset_dims_max(rank), newsize(rank), start(rank), count(rank))
   call h5sget_simple_extent_dims_f(filespace_id, dset_dims, dset_dims_max, status)

   ! 2) Extend time dimension by 1
   newsize = dset_dims
   newsize(1) = newsize(1) + 1
   call h5dextend_f(dset_id, newsize, status)
   if (status /= 0) stop "Error extending complex dataset."

   call h5sclose_f(filespace_id, status)
   call h5dget_space_f(dset_id, filespace_id, status)

   ! 3) Select hyperslab for the new "time" slice
   start(1) = newsize(1) - 1
   count(1) = 1
   if (rank /= 3) stop "Append complex: expected rank=3 dataset (time x Nx x 2)."

   start(2) = 0
   count(2) = size(data,1)

   start(3) = 0
   count(3) = 2

   call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, start, count, status)
   if (status /= 0) stop "Error selecting hyperslab for complex."

   ! 4) Create a 3D memspace [1, Nx, 2]
   tempdims = [1, size(data,1), 2]
   call h5screate_simple_f(3, tempdims, memspace_id, status)
   if (status /= 0) stop "Error creating memspace for complex data."

   ! 5) Write
   call h5dwrite_f(                                 &
        dset_id       = dset_id,                    &
        mem_type_id   = H5T_NATIVE_DOUBLE,           &
        buf           = data,                       &
        dims          = tempdims,        &
        hdferr        = status,                     &
        file_space_id = filespace_id,               &
        mem_space_id  = memspace_id)
   if (status /= 0) stop "Error writing complex data."

   call h5sclose_f(memspace_id, status)
   call h5sclose_f(filespace_id, status)
end subroutine h5w_append_complex

end module hdf5_writer
