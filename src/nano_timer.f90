!: integer(8) nano_time() function returns a high resolution wall-clock time
module nano_timer
   use iso_c_binding
   implicit none
   private
   public :: nano_time, elapsed_ns

   type, bind(C) :: timespec
      integer(c_long) :: tv_sec
      integer(c_long) :: tv_nsec
   end type timespec

   interface
      function clock_gettime(clk_id, tp) bind(C, name="clock_gettime")
         import :: c_int, timespec
         integer(c_int),  value :: clk_id
         type(timespec)        :: tp
         integer(c_int)        :: clock_gettime
      end function clock_gettime
   end interface
   integer(c_int), parameter :: CLOCK_MONOTONIC = 1

contains
   function nano_time() result(ns)
      integer(8) :: ns
      integer(c_int) :: cstat
      type(timespec)     :: ts
      !call clock_gettime(CLOCK_MONOTONIC, ts)
      cstat = clock_gettime(CLOCK_MONOTONIC, ts)
      ns = int(ts%tv_sec,  8) * 1000000000_8 + &
           int(ts%tv_nsec, 8)
   end function nano_time
   function elapsed_ns(t_start, t_end) result(dt)
      integer(c_int64_t), intent(in) :: t_start
      integer(c_int64_t), intent(in) :: t_end
      integer(c_int64_t)             :: dt
      dt = t_end - t_start
   end function elapsed_ns

end module nano_timer
