module propagate
  
use variables_global
use analysis
use util
use sort
use io_binary  ! io_bin_test, write_dbin, read_dbin
use hdf5
use hdf5_interface
use nano_timer
  
implicit none


type PropagationPrivate

  !: field info
  real(8) :: dirx1, diry1, dirz1, emax1, efield1
  real(8) :: dirx2, diry2, dirz2, emax2, efield2
  real(8) :: efieldx, efieldy, efieldz
  real(8) :: temp, temp1, temp2

  !: results and psi stuff
  real(8)    :: norm, normV, rate, mux, muy, muz
  !complex(8) :: psi_j, psi_k
  !complex(8) :: psi(nstuse), psi1(nstates)
  real(8),allocatable :: rate_a(:),rate_b(:),rate_aa(:),rate_ab(:),rate_ba(:),rate_bb(:)
  complex(8), allocatable :: density_complex(:)
  real(8), allocatable :: transition_rates(:)

  !: file stuff
  integer(8)   :: funit(7) !: Output data file units
  character(4) :: dirstr, emaxstr
  character(100) :: cifile, datafile
 !: lapack stuff and misc
  integer(8) :: info2, info3
  integer(8) :: wallstart1, wallstart2, wallfinish1, wallfinish2
  real(8)    :: start1, start2, start3, finish1, finish2, finish3
  real(8)    :: h5cumtime, plaincumtime
  !integer(8) :: info1, lscratch, liwork
  !integer(8), allocatable :: iwork(:)
  !real(8), allocatable    :: scratch(:)

  integer(8) :: &
    nva95_MO, nva99_MO, nva95_NO, nva99_NO, nva95_debug, &
    nva99_debug, nva95_MO_sort, nva99_MO_sort, nva95_NO_sort, nva99_NO_sort, &
    nva95_direct, nva99_direct

  real(8) :: rate_density, rate_direct, rate_debug, rate_noabs

  integer(HID_T) :: file_id
  integer(HID_T) :: dir_grpid, field_grpid
  !: real(8) datasets
  integer(HID_T) :: time_dsid, efield1_dsid, efield2_dsid 
  integer(HID_T) :: efieldx_dsid, efieldy_dsid, efieldz_dsid
  integer(HID_T) :: norm2_dsid, mux_dsid, muy_dsid, muz_dsid
  integer(HID_T) :: rate_dsid, pop1_dsid, rate_aa_dsid, rate_ab_dsid, &
    rate_bb_dsid, rate_a_dsid, rate_b_dsid, ion_dsid, ion_coeff_dsid
  integer(HID_T) :: dens_dsid, psi_dsid, psi_det0_dsid
  !: integer(8) datasets
  integer(HID_T) :: nva99MO_dsid, nva99MOsort_dsid, nva99NO_dsid, nva99NOsort_dsid

  contains
    procedure :: initialize => initialize_PropPriv
    procedure :: deconstruct => deconstruct_PropPriv
    !procedure :: write_outdata => write_outdata_PropPriv


end type PropagationPrivate

    

type PropagationShared

  !: Natural Orbital generation
  real(8), allocatable :: opdm_avg(:)  !: Averaged one-particle reduced density matrix (1-RDM)
  real(8), allocatable :: opdm_avg_abs(:)  !: average(abs(opdm))
  real(8), allocatable :: natorb_occ(:) !: Natural orbital occupations (eigenvalues)
  real(8), allocatable :: natorb_occ_abs(:) !: Natural orbital occupations (eigenvalues)
  real(8), allocatable :: U_NO(:) !: MO->NO transformation matrix
  real(8), allocatable :: U_NO_abs(:) !: MO->NO_abs transformation matrix
  integer(8) :: opdm_avg_N, datetime_start(8)
  real(8), allocatable :: U_NO_input(:) !: U_NO read from file

  real(8) :: ratemax_density, ratemax_direct, ratemax_debug
  

  contains
    procedure :: initialize => initialize_PropShared
    procedure :: deconstruct => deconstruct_PropShared

end type PropagationShared


contains !: MODULE SUBROUTINES BELOW

!==================================================================!
! PropagationPrivate member subroutines start
!==================================================================!
subroutine initialize_PropPriv(this)
  implicit none

  class(PropagationPrivate), intent(inout) :: this

  !write(iout, *) "Allocating Priv!"
  allocate( this%rate_aa(noa*noa), this%rate_ab(noa*nob), this%rate_ba(nob*noa), this%rate_bb(nob*nob) )
  allocate( this%rate_a(noa+nva), this%rate_b(nob+nvb+2) )
  allocate( this%density_complex(nrorb*nrorb) )
  allocate( this%transition_rates(nrorb*nrorb) )
  this%rate_aa = 0.d0 ; this%rate_ab = 0.d0 ; this%rate_ba = 0.d0 ; this%rate_bb = 0.d0
  this%rate_a = 0.d0 ; this%rate_b = 0.d0

  this%dirx1 = 0.d0 ; this%diry1 = 0.d0 ; this%dirz1 = 0.d0
  this%dirx2 = 0.d0 ; this%diry2 = 0.d0 ; this%dirz2 = 0.d0
  this%emax1 = 0.d0; this%efield1 = 0.d0; this%emax2 = 0.d0; this%efield2 = 0.d0
  this%efieldx = 0.d0 ; this%efieldy = 0.d0; this%efieldz = 0.d0
  this%temp = 0.d0; this%temp1 = 0.d0; this%temp2 = 0.d0
  this%norm = 0.d0; this%normV = 0.d0; this%rate = 0.d0
  this%mux = 0.d0; this%muy = 0.d0 ; this%muz = 0.d0

  this%nva95_MO = 0 ; this%nva99_MO = 0 ; this%nva95_NO = 0 ; this%nva99_NO = 0
  this%nva95_direct = 0 ; this%nva99_direct = 0
  this%nva95_debug = 0 ; this%nva99_debug = 0
  this%nva95_MO_sort = 0
  this%nva99_MO_sort = 0
  this%nva95_NO_sort = 0
  this%nva99_NO_sort = 0    
  this%rate_density = 0.d0 ; this%rate_direct = 0.d0 ; this%rate_debug = 0.d0
  this%rate_noabs = 0.d0

  this%h5cumtime = 0.d0 ; this%plaincumtime = 0.d0

  !if( QeigenDC ) then
  !  allocate( this%iwork(3+5*nstuse) )
  !  allocate( this%scratch(1+8*nstuse+2*nstuse*nstuse) )
  !else
  !  allocate( this%iwork(2) )
  !  allocate( this%scratch(nstuse*nstuse) )
  !end if

  !allocate( this%psi(nstuse) )
  !allocate( this%psi1(nstates) )



end subroutine initialize_PropPriv

subroutine deconstruct_PropPriv(this)
  implicit none
  class(PropagationPrivate), intent(inout) :: this

  if(h5inc_enable) then
    !$OMP CRITICAL (IO_LOCK)
    call h5_close_group(this%dir_grpid)
    call h5_close_file(this%file_id)
    !$OMP END CRITICAL (IO_LOCK)
  end if

  deallocate( this%density_complex, this%transition_rates )
  deallocate( this%rate_aa(noa*noa), this%rate_ab(noa*nob) )
  deallocate( this%rate_ba(nob*noa), this%rate_bb(nob*nob) )
  deallocate( this%rate_a(noa+nva), this%rate_b(nob+nva+2) )

end subroutine deconstruct_PropPriv

subroutine write_h5metadata(datetime_start)
  implicit none
  integer(8), intent(in) :: datetime_start(8)

  integer(HID_T) :: file_id, group_id
  logical :: file_exists
  character(len=2048) :: opath, group_path
  character(len=128) :: str_id
  integer(8) :: datetime_end(8)

  write(iout, *) "Writing metadata.h5"
  write(opath, '( "metadata.h5" )')
  write(group_path, '( "/metadata" )')
  inquire(file=opath, exist=file_exists)
  if (file_exists) then
    call execute_command_line("rm -f " // opath) ! I don't think this works.
  end if
  call h5_open_file(trim(opath), file_id)
  call h5_create_group(file_id, group_path, group_id)
  call h5_write_attribute_int(group_id, "nemax", [nemax])
  call h5_write_attribute_int(group_id, "ndir", [ndir])
  call h5_write_attribute_int(group_id, "nbasis", [nbasis])
  call h5_write_attribute_int(group_id, "nrorb", [nrorb])
  call h5_write_attribute_int(group_id, "noa", [noa])
  call h5_write_attribute_int(group_id, "nob", [nob])
  call h5_write_attribute_int(group_id, "nva", [nva])
  call h5_write_attribute_int(group_id, "nvb", [nvb])
  call h5_write_attribute_int(group_id, "nstates", [nstates])
  call h5_write_attribute_int(group_id, "nstuse", [nstuse])
  call h5_write_attribute_int(group_id, "ip_states", [ip_states])
  call DATE_AND_TIME(VALUES=datetime_end)
  call h5_write_attribute_int(group_id, "datetime_start", datetime_start)
  call h5_write_attribute_int(group_id, "datetime_end", datetime_end)

  call h5_write_static_real(group_id, "orben", Mol%orben)
  call h5_write_static_real(group_id, "vabsmoa", Mol%vabsmoa)
  call h5_write_static_real(group_id, "cmo_a", Mol%cmo_a)
  call h5_write_static_real(group_id, "dipxmoa", Mol%dipxmoa)
  call h5_write_static_real(group_id, "dipymoa", Mol%dipymoa)
  call h5_write_static_real(group_id, "dipzmoa", Mol%dipzmoa)

  call h5_write_static_int(group_id, "hole_index", hole_index)
  call h5_write_static_int(group_id, "part_index", part_index)
  call h5_write_static_int(group_id, "hole_ip_index", hole_ip_index)
  call h5_write_static_int(group_id, "part_ip_index", part_ip_index)


  call h5_close_group(group_id)
  call h5_close_file(file_id)

end subroutine write_h5metadata



!: Create new h5 file for the direction thread
subroutine init_h5dir(Priv, idir)
  implicit none
  class(PropagationPrivate), intent(inout) :: Priv
  integer(8), intent(in) ::  idir

  integer(HID_T) :: file_id, grp_id, dset_id_real
  character(len=128) :: gpath, fpath
  logical :: file_exists

  write(fpath, '( "thread", i0, ".h5" )') idir
  inquire(file=fpath, exist=file_exists)
  if (file_exists) then
    call execute_command_line("rm -f " // fpath)
  end if
  file_id = -1
  call h5_open_file(fpath, file_id)
  Priv%file_id = file_id

  write(gpath, '( "/direction_", i0 )') idir
  call h5_create_group(file_id, gpath, grp_id)
  call h5_write_attribute_real(grp_id, &
         "field_direction", [Priv%dirx1,Priv%diry1,Priv%dirz1])
  Priv%dir_grpid = grp_id

end subroutine init_h5dir

!: Create the /directionX/fieldY direction
subroutine init_h5emax(Priv, iemax, idir)
  class(PropagationPrivate), intent(inout) :: Priv
  integer(8), intent(in) :: iemax, idir
  integer(HSIZE_T), dimension(1) :: one_dim, nrorb_dim, nrorb2_dim, norb_dim
  integer(HSIZE_T), dimension(1) :: noa_dim, nob_dim, noa2_dim, noanob_dim, nob2_dim, &
    ip_states_dim, rate_b_dim, ion_coeff_dim, nstuse_dim, nstates_dim
  
  character(len=128) :: gpath
  integer(HID_T) :: grp_id

  grp_id = -1
  write(gpath, '( "/direction_", i0, "/field_", i0 )') idir, iemax
  !write(iout, *) "before create emax group ", Priv%file_id, gpath ; flush(iout)

  call h5_create_group(Priv%dir_grpid, gpath, grp_id)
  Priv%field_grpid = grp_id
  !write(iout, *) "emax group created: ", grp_id, Priv%field_grpid ; flush(iout)

  !: We have to declare these as integer(HSIZE_T) arrays
  one_dim = [1]
  nrorb_dim = [nrorb]
  norb_dim = [norb]
  nrorb2_dim = [nrorb*nrorb]
  noa2_dim = [noa*noa]
  noanob_dim = [noa*nob]
  nob2_dim = [nob*nob]
  noa_dim = [noa]
  nob_dim = [nob]
  ip_states_dim = [ip_states]
  rate_b_dim = [nrorb+2]
  ion_coeff_dim = [max(2*ip_states*(nva+nvb),ip_states*ip_states)]
  nstuse_dim = [nstuse]
  nstates_dim = [nstates]

  call h5_create_dataset_real(grp_id, "time", one_dim, Priv%time_dsid)
  call h5_create_dataset_real(grp_id, "efield1", one_dim, Priv%efield1_dsid)
  call h5_create_dataset_real(grp_id, "efield2", one_dim, Priv%efield2_dsid)
  call h5_create_dataset_real(grp_id, "efieldx", one_dim, Priv%efieldx_dsid)
  call h5_create_dataset_real(grp_id, "efieldy", one_dim, Priv%efieldy_dsid)
  call h5_create_dataset_real(grp_id, "efieldz", one_dim, Priv%efieldz_dsid)
  call h5_create_dataset_real(grp_id, "norm2", one_dim, Priv%norm2_dsid)
  call h5_create_dataset_real(grp_id, "mux", one_dim, Priv%mux_dsid)
  call h5_create_dataset_real(grp_id, "muy", one_dim, Priv%muy_dsid)
  call h5_create_dataset_real(grp_id, "muz", one_dim, Priv%muz_dsid)
  call h5_create_dataset_real(grp_id, "rate", one_dim, Priv%rate_dsid)

  call h5_create_dataset_int(grp_id, "nva99MO", one_dim, Priv%nva99MO_dsid)
  call h5_create_dataset_int(grp_id, "nva99MOsort", one_dim, Priv%nva99MOsort_dsid)
  call h5_create_dataset_int(grp_id, "nva99NO", one_dim, Priv%nva99NO_dsid)
  call h5_create_dataset_int(grp_id, "nva99NOsort", one_dim, Priv%nva99NOsort_dsid)
  if(h5inc_density) then
    call h5_create_dataset_real(grp_id, "density", nrorb2_dim, Priv%dens_dsid)
  end if
  if(h5inc_psi) then
    call h5_create_dataset_complex(grp_id, "psi", nstuse_dim, Priv%psi_dsid)
  end if
  if(h5inc_psi_det0) then
    call h5_create_dataset_complex(grp_id, "psi_det0", nstates_dim, Priv%psi_det0_dsid)
  end if

  !: funit(3)
  call h5_create_dataset_real(grp_id, "pop1", norb_dim, Priv%pop1_dsid)
  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then

    call h5_create_dataset_real(grp_id, "rate_aa", noa2_dim, Priv%rate_aa_dsid)
    call h5_create_dataset_real(grp_id, "rate_ab", noanob_dim, Priv%rate_ab_dsid)
    call h5_create_dataset_real(grp_id, "rate_bb", nob2_dim, Priv%rate_bb_dsid)

  else !: flag_cis, flag_soc
    call h5_create_dataset_real(grp_id, "rate_a", nrorb_dim, Priv%rate_a_dsid)
    call h5_create_dataset_real(grp_id, "rate_b", rate_b_dim, Priv%rate_b_dsid)
  end if

  !: funit(4)
  call h5_create_dataset_real(grp_id, "ion", norb_dim, Priv%ion_dsid)
  call h5_create_dataset_complex(grp_id, "ion_coeff", ion_coeff_dim, Priv%ion_coeff_dsid)


end subroutine init_h5emax


subroutine write_h5_step(Priv, psi, psi1, psi_det0, Zion_coeff, &
                                  ion_coeff, scratch, itime, idata, pop1, ion)
  implicit none

  class(PropagationPrivate), intent(inout) :: Priv
  
  complex(8), intent(inout) :: psi(nstuse), psi1(nstates), psi_det0(:)
  complex(8), intent(inout) :: Zion_coeff(:), ion_coeff(:)
  real(8), intent(inout) :: scratch(:)

  integer(8), intent(inout) :: itime, idata
  real(8), intent(inout) :: pop1(:), ion(:)
  integer(8) :: ndim2

  ndim2 = (noa+nva)*(noa+nva)

  call h5_append_real(Priv%time_dsid , [dble(itime)*dt*au2fs])
  call h5_append_real(Priv%efield1_dsid , [Priv%efield1])
  call h5_append_real(Priv%efield2_dsid , [Priv%efield2])
  call h5_append_real(Priv%efieldx_dsid , [Priv%efieldx])
  call h5_append_real(Priv%efieldy_dsid , [Priv%efieldy])
  call h5_append_real(Priv%efieldz_dsid , [Priv%efieldz])
  call h5_append_real(Priv%norm2_dsid , [Priv%norm**2])
  call h5_append_real(Priv%mux_dsid , [Priv%mux])
  call h5_append_real(Priv%muy_dsid , [Priv%muy])
  call h5_append_real(Priv%muz_dsid , [Priv%muz])
  call h5_append_real(Priv%rate_dsid , [Priv%rate/au2fs])

  call h5_append_int(Priv%nva99MO_dsid , [Priv%nva99_MO])
  call h5_append_int(Priv%nva99MOsort_dsid , [Priv%nva99_MO_sort])
  call h5_append_int(Priv%nva99NO_dsid , [Priv%nva99_NO])
  call h5_append_int(Priv%nva99NOsort_dsid , [Priv%nva99_NO_sort])
  if(h5inc_density) then
    call h5_append_real(Priv%dens_dsid, scratch(:ndim2))
  end if
  if(h5inc_psi) then
    call h5_append_complex(Priv%psi_dsid, psi(:nstuse))
  end if
  if(h5inc_psi_det0) then
    call h5_append_complex(Priv%psi_det0_dsid, psi_det0)
  end if

  call h5_append_real(Priv%pop1_dsid, pop1)
  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    call h5_append_real(Priv%rate_aa_dsid, Priv%rate_aa)
    call h5_append_real(Priv%rate_ab_dsid, Priv%rate_ab)
    call h5_append_real(Priv%rate_bb_dsid, Priv%rate_bb)
  else !: flag_cis, flag_soc
    call h5_append_real(Priv%rate_a_dsid, Priv%rate_a)
    call h5_append_real(Priv%rate_b_dsid, Priv%rate_b)
  end if
  call h5_append_real(Priv%ion_dsid, ion)
  !write(*,*) "test", max(2*ip_states*(nva+nvb),ip_states*ip_states), size(ion_coeff)
  call h5_append_complex(Priv%ion_coeff_dsid, ion_coeff)



end subroutine write_h5_step


!: Run once for each thread at start of propagation
!: Sets up the funit values and writes headers to .dat files.
subroutine PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, psi_det0, Zion_coeff, MODE )
  implicit none

  class(PropagationPrivate), intent(inout) :: Priv
  integer(8), intent(in) :: iemax, idir
  type(tdcidat), allocatable, intent(inout) :: tdciresults(:)
  complex(8), intent(inout) :: psi0(:), psi_det0(:), Zion_coeff(:)
  !: MODE 0=trotter_linear,  1=trotter_circular
  !:      2=Ztrotter_linear, 3=Ztrotter_circular
  integer, intent(in) :: MODE

  !: Local variables
  integer(8) :: i,j

  !: Aliases
  !: Apparently the PGI compiler doesn't like this...
  !integer(8), pointer :: funit(:) => Priv%funit
  !character(4), pointer :: dirstr => Priv%dirstr
  !character(4), pointer :: emaxstr => Priv%emaxstr
  !character(100), pointer :: cifile => Priv%cifile
  !character(100), pointer :: datafile => Priv%datafile

  !: for writing out files later
  write( Priv%emaxstr, '(i0)' ) iemax
  write( Priv%dirstr, '(i0)' )  idir

  Priv%funit(1) = 0    +1000*iemax + idir
  Priv%funit(2) = 10000+1000*iemax + idir
  Priv%funit(3) = 20000+1000*iemax + idir
  Priv%funit(4) = 30000+1000*iemax + idir
  Priv%funit(5) = 40000+1000*iemax + idir
  Priv%funit(6) = 50000+1000*iemax + idir
  Priv%funit(7) = 60000+1000*iemax + idir


  !: cifile binary
  if( Qci_save ) then
    !Priv%cifile ='CI-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.bin'
    !open( unit=Priv%funit(1), file=trim(Priv%cifile), form='unformatted' )
    !write(Priv%funit(1)) ndata, nstuse, nstates
    Priv%cifile ='CI-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
    open( unit=Priv%funit(1), file=trim(Priv%cifile) )
    write(Priv%funit(1), "(A,i6,A,i6,A,i6)") "ndata=",ndata,"  nstuse=",nstuse, &
     "  nstates=",nstates
    do j=1,nstuse
      write(Priv%funit(1),"(A,i6,A)") "CI state",j,"   Eigenvector:"
      write(Priv%funit(1),"(5d18.10)") (cis_vec((j-1)*nstates+i),i=1,nstates)
      end do
    write(Priv%funit(1), "(A)") "time,field1,field2,dirx1,diry1,dirz1,dirx2,diry2,dirz2,norm2"
    select case (MODE)
      case(0, 2) !: trotter_linear, Ztrotter_linear
        write(Priv%funit(1), "(10f10.6)") 0.d0, 0.d0, 0.d0, Priv%dirx1, Priv%diry1, &
          Priv%dirz1, 0.d0, 0.d0, 0.d0, 1.d0
      case(1, 3) !: trotter_circular, Ztrotter_circular
        write(Priv%funit(1),"(10f10.6)") 0.d0, 0.d0, Priv%dirx1, Priv%diry1, &
          Priv%dirz1, Priv%dirx2, Priv%diry2, Priv%dirz2, 1.d0
    end select
  !:  write(Priv%funit(1)) real(psi0)
  !:  write(Priv%funit(1)) aimag(psi0) 
    write(Priv%funit(1), "(A)") &
      "State,Real(psi),Imag(psi),hole_index(1),hole_index(2),particle_index(1),Real(psi_det),Imag(Psi_det)"
    do j=1,nstates
      if(j.le.nstuse) then
        write(Priv%funit(1),"(i6,2d18.10,3i4,2d18.10)") &
          j,real(psi0(j)),aimag(psi0(j)), &
          hole_index(j,1),hole_index(j,2),part_index(j,1), &
          real(psi_det0(j)),aimag(psi_det0(j))
      else
        write(Priv%funit(1),"(42x,3i4,2d18.10)") &
          hole_index(j,1),hole_index(j,2),part_index(j,1), &
          real(psi_det0(j)),aimag(psi_det0(j))
      end if
    end do
    flush(Priv%funit(1))      
  end if    
  
  !: RESULTS datafile
  Priv%datafile = 'RESULTS-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(2),file=trim(Priv%datafile) )
  write( Priv%funit(2), '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'

  select case (MODE)
    case(0, 2) !: trotter_linear, Ztrotter_linear
      write( Priv%funit(2), "( '#',  20(f16.10,1x) )" ) Priv%emax1, 0.d0, &
             tdciresults(idir+(iemax-1)*ndir)%theta0, tdciresults(idir+(iemax-1)*ndir)%phi0,&
             0.d0,0.d0,  0.d0,0.d0, &
             Priv%dirx1, Priv%diry1, Priv%dirz1,  0.d0,0.d0,0.d0,  0.d0,0.d0,0.d0
    case(1, 3) !: trotter_circular, Ztrotter_circular
      write( Priv%funit(2), "( '#',  20(f16.10,1x) )" ) Priv%emax1, Priv%emax2, &
           tdciresults(1+(iemax-1)*ndir)%theta0, tdciresults(1+(iemax-1)*ndir)%phi0, &
           tdciresults(1+(iemax-1)*ndir)%theta1, tdciresults(1+(iemax-1)*ndir)%phi1, &
           tdciresults(1+(iemax-1)*ndir)%theta2, tdciresults(1+(iemax-1)*ndir)%phi2, &
           tdciresults(1+(iemax-1)*ndir)%x0, tdciresults(1+(iemax-1)*ndir)%y0, tdciresults(1+(iemax-1)*ndir)%z0, & 
           Priv%dirx1, Priv%diry1, Priv%dirz1,  Priv%dirx2, Priv%diry2, Priv%dirz2
  end select

  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
    write( Priv%funit(2),"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)', & 
      'NO99 MO99 ','field1','field2','fieldx','fieldy','fieldz', &
      'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)','6 x |psi(i)|**2'
  else
    write( Priv%funit(2),"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)', &
      'NO99 MO99 ','field1','field2','fieldx','fieldy','fieldz', &
      'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
  end if

  !: POP datafile
  Priv%datafile = 'POP-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(3),file=trim(Priv%datafile) )
  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    write( Priv%funit(3),"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
      'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
  else
    write( Priv%funit(3),"(a5,12(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
      'rate(fs-1)','rate_a(occ)','rate_b(occ)','rate_a(virt)','rate_b(virt)','c0* c0 Vabs00','rate'
  end if

  !: ION datafile
  Priv%datafile = 'ION-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(4),file=trim(Priv%datafile) )
  write( Priv%funit(4),"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'

  !: ION_COEFF datafile
  if( Qread_ion_coeff .or. Qwrite_ion_coeff ) then
    Priv%datafile = 'ION_COEFF-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.bin'
    open( unit=Priv%funit(5),file=trim(Priv%datafile),form='unformatted' )
  else
    Priv%datafile = 'ION_COEFF-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
    open( unit=Priv%funit(5),file=trim(Priv%datafile) )
    write( Priv%funit(5),"(a5,2(a14,1x),3(a24,1x))" ) '#','time(fs)','rate/norm2', &
      's(i)(i=1,ip_states)','proj Zion_coeff(j,i)','Zion_coeff(j,i)'
    write(Priv%funit(5),"('ntimes ',i0)") int(nstep/outstep)
    write(Priv%funit(5),"('ip_states ',i0)") ip_states
  end if
  if( Qread_ion_coeff ) then 
    read(Priv%funit(5)) Zion_coeff
  end if

  !: MO density datafile
  if( Qmo_dens ) then
    Priv%datafile = 'MO_density-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
    open( unit=Priv%funit(6),file=trim(Priv%datafile) )
    !open( newunit=Priv%funit(6),file=trim(Priv%datafile) )
    write(Priv%funit(6),"('ntimes ',i0)") int(nstep/outstep)
    write(Priv%funit(6),"('alpha_homo ',i0)") noa
    write(Priv%funit(6),"('beta_homo ',i0)")  nob
    write(Priv%funit(6),"('alpha_orbitals ',i0)") noa + nva
  end if 


end subroutine PropWriteDataHeaders


subroutine PropWriteData(Priv, psi, psi1, psi_det0, Zion_coeff, ion_coeff, scratch, &
                         itime, idata, pop1, ion)
  implicit none

  class(PropagationPrivate), intent(inout) :: Priv
  
  complex(8), intent(inout) :: psi(nstuse), psi1(nstates), psi_det0(:)
  complex(8), intent(inout) :: Zion_coeff(:), ion_coeff(:)
  real(8), intent(inout) :: scratch(:)

  integer(8), intent(inout) :: itime, idata
  real(8), intent(inout) :: pop1(:), ion(:)

  integer(8) :: i, j, kk
  !real(8) :: tmprate
  logical :: fileopen

  real(8), allocatable :: density_AO(:), density_MO(:)
  character(len=1024) :: density_filename
  

  if( Qci_save) then
  !:  write(Priv%funit(1)) dble(itime)*dt*au2fs, Priv%efield1, 0.d0, &
  !:   Priv%dirx1, Priv%diry1, Priv%dirz1, 0.d0, 0.d0, 0.d0, Priv%norm**2
    write(Priv%funit(1), "(A)") "time,field1,field2,dirx1,diry1,dirz1,dirx2,diry2,dirz2,norm2"
    write(Priv%funit(1), "(10f10.6)") dble(itime)*dt*au2fs, Priv%efield1, Priv%efield2, &
     Priv%dirx1, Priv%diry1, Priv%dirz1, Priv%dirx2, Priv%diry2, Priv%dirz2, Priv%norm**2
 !:   write(Priv%funit(1)) real( psi )
 !:   write(Priv%funit(1)) aimag( psi )
    write(Priv%funit(1), "(A)") &
      "State,Real(psi),Imag(psi),hole_index(1),hole_index(2),particle_index,Real(psi_det),Imag(Psi_det)"
    do j=1,nstates
      if(j.le.nstuse) then
        write(Priv%funit(1),"(i6,2d18.10,3i4,2d18.10)") &
          j,real(psi(j)),aimag(psi(j)), &
          hole_index(j,1),hole_index(j,2),part_index(j,1), &
          real(psi_det0(j)),aimag(psi_det0(j))
      else
        write(Priv%funit(1),"(42x,3i4,2d18.10)") &
          hole_index(j,1),hole_index(j,2),part_index(j,1), &
          real(psi_det0(j)),aimag(psi_det0(j))
      end if
    end do
    flush(Priv%funit(1))
  end if

  !tmprate = 0.d0
  !do i=0, (noa+nva)*(noa+nva)
  !  tmprate = tmprate + scratch(i)*Mol%vabsmoa(i)
  !end do

  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    if(Priv%norm.ne.0) then
      write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs,Priv%nva99_NO,Priv%nva99_MO,Priv%efield1,Priv%efield2,&
             Priv%efieldx,Priv%efieldy,Priv%efieldz, &
             Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz, &
             !Priv%norm**2, tmprate, Priv%rate, Priv%rate/au2fs, Priv%muz, &
             (dble(dconjg(psi(i))*psi(i))/Priv%norm**2,i=1,20)
    else !: norm is zero? when does this happen?
      write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs,Priv%nva99_NO,Priv%nva99_MO,Priv%efield1,Priv%efield2,&
             Priv%efieldx,Priv%efieldy,Priv%efieldz, &
             Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz,&
             0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
    end if
  else
    write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs,Priv%nva99_NO,Priv%nva99_MO,Priv%efield1,Priv%efield2,&
           Priv%efieldx,Priv%efieldy,Priv%efieldz, &
           Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz 
           !Priv%norm**2, tmprate, Priv%rate, Priv%rate/au2fs, Priv%muz 
  end if
  flush(Priv%funit(2))

  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    write( Priv%funit(3),"( i5,f10.4,500(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs, Priv%norm**2, &
           (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),Priv%rate/au2fs,&
           (Priv%rate_aa(i),i=1,noa*noa),&
           (Priv%rate_ab(i),i=1,noa*nob),(Priv%rate_ba(i),i=1,nob*noa),(Priv%rate_bb(i),i=1,nob*nob)
  else !: flag_cis, flag_soc
    write( Priv%funit(3),"( i5,f10.4,5000(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs, Priv%norm**2, &
           (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob), &
           (Priv%rate_a(i),i=1,noa),(Priv%rate_b(i),i=1,nob), &
           (Priv%rate_a(i+noa),i=1,nva),(Priv%rate_b(i+nob),i=1,nvb)
  end if
  flush(Priv%funit(3))

  write( Priv%funit(4),"( i5,f10.4,500(1x,f15.10))") &
         idata, dble(itime)*dt*au2fs, Priv%rate/Priv%norm**2, Priv%normV, &
         (ion(i),i=1,noa),(ion(noa+nva+i),i=1,nob),(ion_coeff(i),i=1,ip_states)
  flush(Priv%funit(4))
              
  if( Qmo_dens) then
    fileopen = .false.
    inquire( unit=Priv%funit(6), opened=fileopen )
    if (.not. fileopen) then
      write(iout, *) "funit(6) isn't opened! ", fileopen
    end if
    call write_density_difference( Priv%funit(6), dble(itime)*dt*au2fs, &
                                   Priv%rate, noa, nva, scratch, Mol%vabsmoa )
  end if

  if ( write_binaries ) then
    write( density_filename, '(A,A,A,A,A,I0,A)') "matrices/MO_density-e", trim(Priv%emaxstr), &
           "-d", trim(Priv%dirstr), ".", itime, ".bin"
    call write_density_bin( density_filename, dble(itime)*dt*au2fs, &
                            Priv%rate, noa, nva, scratch, Mol%vabsmoa, Priv%funit(7) )

  end if ! write_binaries


  if (.not. Qwrite_ion_coeff) then
    select case (trim(jobtype))
    case( flag_cis, flag_soc )
      call get_ion_coeff(hole_index,part_index,psi_det0,psi1,&
             Priv%norm,Mol%vabsmoa,Mol%vabsmob, &
             Priv%rate,Zion_coeff,ion_coeff,scratch)
      if( trim(jobtype).eq.flag_cis ) then
        call get_proj_ion(iout,noa+nob,ip_vec,&
               Zion_coeff(noa+nob+1),Zproj_ion)
      else !: flag_soc
        call get_Zproj_ion(iout,noa+nob,Zip_vec,&
               Zion_coeff(noa+nob+1),Zproj_ion)
      end if
      !write(iout,*) " ip_states",ip_states
      !do i=1,noa+nob
      !  write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
      !  write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
      !  write(iout,"('ipvec',i3,16f13.7)") i,(ip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
      !  write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
      !end do
      write( Priv%funit(5),"(i5,f10.4,50(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs, Priv%rate,(scratch(i),i=1,noa+nob)
      do i = 1,noa+nob
        write( Priv%funit(5),"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
      end do
      do i = 1,noa+nob
        write( Priv%funit(5),"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
      end do
    case( flag_ip, flag_socip )
      call get_ion_coeff_ip(hole_index,part_index,&
             state_ip_index,psi_det0,psi1,Priv%norm,&
             Mol%vabsmoa,Mol%vabsmob,&
             Priv%rate,Zion_coeff,ion_coeff,scratch)
      if( trim(jobtype).eq.flag_ip ) then
        call get_proj_ion(iout,ip_states,ip_vec,&
               Zion_coeff(ip_states+1),Zproj_ion)
      else !: flag_socip
        call get_Zproj_ion(iout,ip_states,Zip_vec,&
          Zion_coeff(noa+nob+1),Zproj_ion)
      end if
      !write(iout,*) " ip_states",ip_states
      !write(iout,"('s    ',14f13.7/5x,14f13.7/5x,14f13.7/5x,14f13.7)") &
      !            (scratch(j),j=1,ip_states)
      !write(iout,"('rates',14f13.7/5x,14f13.7/5x,14f13.7/5x,14f13.7)") &
      !            (abs(Zion_coeff(j)),j=1,ip_states)
      !flush(iout)
      !do i=1,ip_states
      !write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
      !write(iout,"('coeff',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") &
      !            i,(Zion_coeff(j+i*(ip_states)),j=1,ip_states)
      !write(iout,"('ipvec',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") &
      !            i,(ip_vec(j+(i-1)*(ip_states)),j=1,ip_states)
      !write(iout,"('proj ',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") &
      !            i,(Zproj_ion(j+(i-1)*(ip_states)),j=1,ip_states)
      !end do
      !flush(iout)
      write( Priv%funit(5),"(i5,f10.4,50(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs, Priv%rate,(scratch(i),i=1,ip_states)
      do i = 1,ip_states
        write( Priv%funit(5),"(60(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
      end do
      do i = 1,ip_states
        write( Priv%funit(5),"(60(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
      end do
      flush(Priv%funit(5))
    end select
  end if !: (.not. Qwrite_ion_coeff)

end subroutine PropWriteData

subroutine PrivSetField(Priv, itime, iemax)
  implicit none
  class(PropagationPrivate), intent(inout) :: Priv
  integer(8), intent(in) :: itime, iemax
  integer(8) :: i 
  
  
  !: modified midpoint 
  Priv%efield1 = 0.5d0 * Priv%emax1*(fvect1(itime)+fvect1(itime+1)) !: modified midpoint
  Priv%efield2 = 0.5d0 * Priv%emax2*(fvect2(itime)+fvect2(itime+1)) !: modified midpoint
  if( read_shift(iemax).ne.0 ) then
    if(itime.eq.1) write(iout,"(' Pulse has been shifted by ',i6,' steps')") read_shift(iemax)
    i = itime-read_shift(iemax)
    if( i.gt.0 .and. i.lt.nstep ) then
      Priv%efield1 = 0.5d0 * Priv%emax1 * ( fvect1(i) + fvect1(i+1) )
      Priv%efield2 = 0.5d0 * Priv%emax2 * ( fvect2(i) + fvect2(i+1) )
    else
      Priv%efield1 = 0.d0
      Priv%efield2 = 0.d0
    end if
  end if
  Priv%efieldx = Priv%dirx1 * Priv%efield1
  Priv%efieldy = Priv%diry1 * Priv%efield1
  Priv%efieldz = Priv%dirz1 * Priv%efield1
  if(.not.linear) then
    Priv%efieldx = Priv%efieldx + Priv%dirx2 * Priv%efield2 
    Priv%efieldy = Priv%efieldy + Priv%diry2 * Priv%efield2
    Priv%efieldz = Priv%efieldz + Priv%dirz2 * Priv%efield2
    Priv%temp1 = dt * Priv%efield1 * 0.5d0 
    Priv%temp2 = dt * Priv%efield2
  end if
    
end subroutine PrivSetField

!==================================================================!
!  PropagationPrivate subroutines end
!==================================================================!

!==================================================================!
!  PropagationShared subroutines start
!==================================================================!
subroutine initialize_PropShared(Prop)
  implicit none

  class(PropagationShared), intent(inout) :: Prop

  integer(8) :: info

  real(8) :: min_temp, max_temp, dip_temp
  integer(8) :: i, j, i_min, i_max, j_min, j_max
  
  !: Natural Orbital initialization
  allocate( Prop%opdm_avg(nrorb*nrorb) )
  allocate( Prop%opdm_avg_abs(nrorb*nrorb) )
  allocate( Prop%natorb_occ(nrorb))
  allocate( Prop%natorb_occ_abs(nrorb))
  allocate( Prop%U_NO(nrorb*nrorb) )
  allocate( Prop%U_NO_abs(nrorb*nrorb) )
  allocate( Prop%U_NO_input(nrorb*nrorb) )
  Prop%opdm_avg = 0.d0
  Prop%opdm_avg_abs = 0.d0
  Prop%natorb_occ = 0.d0
  Prop%natorb_occ_abs = 0.d0
  Prop%U_NO = 0.d0
  Prop%U_NO_abs = 0.d0
  Prop%opdm_avg_N = ndir*(nstep/outstep)
  !write(iout, '(A,i5)') "opdm_avg_N: ", Prop%opdm_avg_N

  !if ( .true. ) then
  if ( flag_ReadU_NO ) then
    call read_dbin( Prop%U_NO_input, (nrorb*nrorb), "U_NO_input.bin", info)
  end if

  if(datfile_enable) then
    call write_dbin_safe( Mol%vabsmoa, nrorb*nrorb, "matrices/Vabs_MO.bin" )
    call write_dbin_safe( Mol%cmo_a, nbasis*nrorb, "matrices/CMO.bin" )
    call write_dbin_safe( Mol%dipxmoa, nbasis*nrorb, "matrices/dipxmoa.bin" )
    call write_dbin_safe( Mol%dipymoa, nbasis*nrorb, "matrices/dipymoa.bin" )
    call write_dbin_safe( Mol%dipzmoa, nbasis*nrorb, "matrices/dipzmoa.bin" )
  end if

  if(verbosity.gt.1) then
    write(iout, "('Maximum(cis_eig)=',f10.4)") maxval(cis_eig, 1)
    write(iout, "('Minimum(cis_eig)=',f10.4)") minval(cis_eig, 1)
  end if
  min_temp = 0.0 ; max_temp = 0.0
  i_min = 0 ; i_max = 0 ; j_min = 0 ; j_max = 0
  do i = 1,nrorb
    do j = 1,nrorb
      if (Mol%vabsmoa((i-1)*nrorb+j) > max_temp) then
        max_temp = Mol%vabsmoa((i-1)*nrorb+j)
        i_max = i
        j_max = j
      end if
      if (Mol%vabsmoa((i-1)*nrorb+j) < min_temp) then
        min_temp = Mol%vabsmoa((i-1)*nrorb+j)
        i_min = i
        j_min = j
      end if
    end do
  end do
  if(verbosity.gt.1) then
    write(iout, "('Maximum(vabsmoa)=',i5,i5,f10.4)") i_max, j_max, max_temp
    write(iout, "('Minimum(vabsmoa)=',i5,i5,f10.4)") i_min, j_min, min_temp
  end if
  ! Dipole
  min_temp = 0.0 ; max_temp = 0.0
  i_min = 0 ; i_max = 0 ; j_min = 0 ; j_max = 0
  do i = 1,nrorb
    do j = 1,nrorb
      flush(iout)
      dip_temp = sqrt( (Mol%dipxmoa((i-1)*nrorb+j))**2 + (Mol%dipymoa((i-1)*nrorb+j))**2 + (Mol%dipzmoa((i-1)*nrorb+j))**2 )
      if (dip_temp > max_temp) then
        max_temp = dip_temp
        i_max = i
        j_max = j
      end if
      if ((i.eq.j) .and. (dip_temp > min_temp)) then
        min_temp = dip_temp
        i_min = i
      end if
    end do
  end do
  if(verbosity.gt.1) then
    write(iout, "('Maximum(MO dipole norm                 )=',i5,i5,f10.4)") i_max, j_max, max_temp
    write(iout, "('Maximum(MO dipole norm (diagonals only))=',i5,i5,f10.4)") i_min, i_min, min_temp
  end if
  call DATE_AND_TIME(VALUES=Prop%datetime_start)



end subroutine initialize_PropShared


subroutine deconstruct_PropShared(this)
  class(PropagationShared), intent(inout) :: this

  deallocate( this%opdm_avg , this%opdm_avg_abs )
  deallocate( this%natorb_occ, this%natorb_occ_abs )
  deallocate( this%U_NO, this%U_NO_abs )
  deallocate( this%U_NO_input )
end subroutine deconstruct_PropShared




!==================================================================!
!==================================================================!
subroutine get_expVabs
  
  !: get ready to propagate.  Get exp(-V_abs)
  
  implicit none
  
  real(8) :: hdt
  integer(8) :: i,j,ij
  

  call write_header( 'get_expVabs', 'propagate','enter' )
  write(iout,'(A)') ' before propagation, getting exp( -V_abs * dt/2 )'
  flush(iout)    

  !: If Heuristic.gt.0, add Heuristic*Sqrt(CISig(I))
  !: to diagonal of abp for unbound states
  !if( heuristic.gt.0.d0 ) then
  !  write(iout,"(' heuristic lifetime =', f8.4,' sqrt(cis_eig-ionization)')") heuristic
  !  do i = 1, nstuse
  !    if( cis_eig(i).gt.ionization ) abp(nstuse*(i-1)+i) = abp(nstuse*(i-1)+i)+heuristic*sqrt(cis_eig(i)-ionization)
  !  end do
  !else
  !  write(iout,'(A)') " heuristic = 0.0 "
  !end if
  
  hdt = -0.5d0*dt
  call cpu_time(start)
  if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
    call Zexp_mat(iout,nstuse,Zabp(1:nstuse*nstuse),Zexp_abp(1:nstuse*nstuse),hdt,QeigenDC)
    call testherm(nstuse,Zexp_abp,'ZexpVabs')
  else
    call exp_mat(iout,nstuse,abp(1:nstuse*nstuse),exp_abp(1:nstuse*nstuse),hdt,QeigenDC)
  end if
  call cpu_time(finish)
  write(iout,"(' exponential of V_abs time:',f12.4,' s')") finish-start       
      
  call write_header( 'get_expVabs','propagate','leave' )
  flush(iout)    
  
  
end subroutine get_expVabs
!==================================================================!
!==================================================================!
subroutine trotter_linear
  
  use omp_lib    
  implicit none
  
  !: thread shared variables
  class(PropagationShared), allocatable :: Prop
  integer(8) :: nstuse2 
  complex(8) :: exphel(nstuse), psi0(nstuse)
  real(8)    :: tdvals1(nstuse)  !: eigenvectors of hp1
  real(8)    :: hp1(nstuse*nstuse)
  
  !: thread private variables
  class(PropagationPrivate), allocatable :: Priv    
  !: Psi arrays must be allocated on the stack for performance
  complex(8) :: psi_j, psi_k, psi(nstuse), psi1(nstates)

  real(8) :: norm0
  real(8),allocatable :: pop0(:),pop1(:),ion(:)
  complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)

  integer(8) :: i, j, ii, jj, k, kk
  integer(8) :: itime, idir, iemax
  integer(8) :: ithread, idata
  complex(8) :: cdum

  real(8) :: start1, finish1, start2, finish2, start3, finish3
  real(8) :: times(100)

  !: lapack stuff
  integer(8) :: info1, lscratch, liwork
  integer(8), allocatable :: iwork(:)
  real(8), allocatable    :: scratch(:)

  !: variables not used in OMP
  real(8) :: temp

  !=------------------------=!

  call write_header( 'trotter_linear','propagate','enter' )    
  call writeme_propagate( 'trot_lin', 'equation' ) 
  call cpu_time(start)
  
  nstuse2 = nstuse * nstuse

  !: allocation and exphel = exp(-iH*dt/2)
  call trotter_init(Prop, exphel, psi0, norm0, pop0, pop1, ion, psi_det0, ion_coeff, Zion_coeff, iwork, scratch)

  !: write psi0
  call writeme_propagate( 'trot_lin', 'psi0' )    

  !: Start loop over directions.

  !$OMP PARALLEL DEFAULT(NONE),&
  !$OMP PRIVATE(Priv, i, idata, idir, iemax, ii, itime, j, jj, k, kk, ithread, &
  !$OMP norm0, lscratch, liwork, scratch, iwork, info1, start1, finish1, &
  !$OMP psi_j, psi_k, psi, psi1, start2, finish2, start3, finish3, times, &
  !$OMP pop1, ion, ion_coeff, psi_det0,  &
  !$OMP hp1, tdvals1, Zion_coeff ),  &
  !$OMP SHARED( Mol, Prop, jobtype, nbasis, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
  !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
  !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
  !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
  !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
  !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
  !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
  !$OMP flag_ReadU_NO, write_orb_transitionrates, h5inc_enable, datfile_enable, verbosity )
  
  !$OMP DO  
  dir_loop : do idir=1, ndir
     
    ithread = omp_get_thread_num()
    call cpu_time( start1 )

    psi = 0 ; psi1 = 0
    tdvals1 = 0.d0 ;

    allocate(Priv)
    call Priv%initialize()

    !: get directions stored in TDCItdciresults
    Priv%dirx1 = tdciresults(idir)%x0  
    Priv%diry1 = tdciresults(idir)%y0  
    Priv%dirz1 = tdciresults(idir)%z0  

    if(h5inc_enable) then
      !$OMP CRITICAL (IO_LOCK)
      !: Create the /direction_idir group in the h5 file
      call init_h5dir(Priv, idir)
      !$OMP END CRITICAL (IO_LOCK)
    end if

    !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
    hp1(:) = Priv%dirx1*tdx(1:nstuse2) &
           + Priv%diry1*tdy(1:nstuse2) + Priv%dirz1*tdz(1:nstuse2)

    !: diagonalize mu dot E
    info1 = 10  
    if( QeigenDC ) then
      lscratch = 1+6*nstuse+2*nstuse*nstuse
      liwork = 3+5*nstuse
      call dsyevd('v','u',nstuse, hp1, nstuse, tdvals1, &
        scratch, lscratch, iwork, liwork, info1)
    else
      lscratch = nstuse*nstuse
      call dsyev('v','u',nstuse, hp1, nstuse, tdvals1, &
        scratch, lscratch, info1)
    end if

    !: hp = W * exp(-Vabs dt/2)
    call dgemm('t','n',nstuse,nstuse,nstuse,1.d0,hp1,nstuse, & 
      exp_abp,nstuse,0.d0, scratch,nstuse)
    hp1 = scratch(1:nstuse*nstuse)

    call cpu_time(finish1)
     
    !: loop over intensities
    emax_loop : do iemax=1, nemax

      !: get emax
      Priv%emax1 = tdciresults(1+(iemax-1)*ndir)%fstrength0

      !$OMP CRITICAL (IO_LOCK)

      if(datfile_enable) then
        call cpu_time(start2)
        call PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, psi_det0, Zion_coeff, 0)
        call cpu_time(finish2)
        Priv%plaincumtime = Priv%plaincumtime + (finish2-start2)
      end if

      if(h5inc_enable) then
        call cpu_time(start2)
        call init_h5emax(Priv, iemax, idir)
        call cpu_time(finish2)
        Priv%h5cumtime = Priv%h5cumtime + (finish2-start2)
      end if

      if(verbosity.gt.1) then
        if(iemax.eq.1) write(iout,"(12x,'Init, TD diag, and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
        write(iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0)" ) idir, iemax, ithread
      end if
      if(verbosity.gt.2) then
        write( iout, "('  Maximum(tdvals1) = ',f12.4)") maxval(tdvals1)
        write( iout, "('  Minimum(tdvals1) = ',f12.4)") minval(tdvals1)
      end if
      flush( iout )

      !$OMP END CRITICAL (IO_LOCK)
                  
      !: initialize psi
      call cpu_time(start3)
      psi = psi0
      if( read_state1(iemax).ne.0 ) then
        psi = dcmplx(0.d0,0.d0)
        write(iout,"(' *** Reset initial wavefunction ***')")
        if( read_state1(iemax).ne.0 ) then
          i = read_state1(iemax)
          psi(i) = read_coeff1(iemax)
          write(iout,"(' Initial coefficient for state',i4,' is ',2f12.8)") i,psi(i)
        end if
        if( read_state2(iemax).ne.0 ) then
          i = read_state2(iemax)
          psi(i) = read_coeff2(iemax)
          write(iout,"(' Initial coefficient for state',i4,' is ',2f12.8)") i,psi(i)
        end if
        call get_norm( norm0, nstuse, psi )
        psi = psi / norm0 
        psi0 = psi
      end if
      if( Qread_ion_coeff ) psi = dcmplx(0.d0,0.d0)
      call cpu_time(finish3)
      times(1) = start3 - finish3

      !: begin Looping over time
      call cpu_time( Priv%start2 )
      timestep_loop : do itime=1, nstep-1

        call cpu_time(start3)
        call PrivSetField(Priv, itime, iemax) 

        !: How can we factor out this ion sampling stuff?
        if( itime.lt.ion_sample_start(iemax) ) psi = dcmplx(0.d0,0.d0)
        if( Qread_ion_coeff ) then
          if( itime.ge.ion_sample_start(iemax) .and. &
              itime.le.ion_sample_start(iemax)+ion_sample_width(iemax) ) then
            ii = (itime-1)*(noa+nob)*(noa+nob+1)
            !: tranform ion_coeff to CI basis and add to psi
            !write(iout,"('R0',i5,i3,16F10.6)") itime,ion_sample_state(iemax)
            call get_ion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                   dt,cis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
          end if
        else
          if( itime.eq.ion_sample_start(iemax) ) psi = psi0
        end if
        call cpu_time(finish3)
        times(2) = finish3 - start3

        !: Below follows Equation (6) from https://doi.org/10.3389/fchem.2022.866137 
        call cpu_time(start3)
        !: psi contains C(t)
        !: psi = exp(-iHel dt/2 ) * C(t)
        do i=1, nstuse
          psi(i) = exphel(i) * psi(i)
        end do
        call cpu_time(finish3)
        times(3) = finish3 - start3

         
        call cpu_time(start3)
        !: psi contains exp(-iHel dt/2 ) * C(t)
        !: psi1 = W * exp(-Vabs dt/2) * psi
        psi1 = dcmplx(0.d0,0.d0)
        do j = 1, nstuse
          jj = ( j-1) * nstuse  
          psi_j = psi(j)
          do k=1, nstuse
            psi1(k) = psi1(k) + hp1( jj+k ) * psi_j
          end do
        end do
        call cpu_time(finish3)
        times(4) = finish3 - start3

         
        call cpu_time(start3)
        !: psi1 contains W * exp(-Vabs dt/2) * exp(-iHel dt/2 ) * C(t)
        !: psi1 = exp(-E(t+dt/2)*mu*dt) * psi1
        do j = 1, nstuse
          Priv%temp = dt * Priv%efield1 * tdvals1(j)
          psi1(j) = dcmplx( dcos(Priv%temp),dsin(Priv%temp) ) * psi1(j)
        end do
        call cpu_time(finish3)
        times(5) = finish3 - start3

         
        call cpu_time(start3)
        !: psi1 contains exp(-E(t+dt/2)*mu*dt) * W 
        !:               * exp(-Vabs dt/2) * exp(-iHel dt/2 ) * C(t)
        !: psi = exp(-iHel dt/2) * exp(-Vabs dt/2) * W.T * psi1
        do j = 1, nstuse
          jj = nstuse * ( j-1 )
          psi_j = dcmplx( 0.d0, 0.d0 )
          do k=1, nstuse
             psi_j  = psi_j + hp1( jj+k ) * psi1(k)
          end do
            psi(j) = exphel(j) * psi_j
        end do
        !: psi now contains the result of Equation (6), C(t+dt) = 
        !:   exp(-iHel dt/2) * exp(-Vabs dt/2)
        !:    * W.T * exp(-E(t+dt/2)*mu*dt) * W
        !:    * exp(-Vabs dt/2) * exp(-iHel dt/2 ) * C(t)
        call cpu_time(finish3)
        times(6) = finish3 - start3
         
        call cpu_time(start3)
        if( Qwrite_ion_coeff ) then
          call get_norm( Priv%norm, nstuse, psi )
          call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )
          ii = (itime-1)*(noa+nob)*(noa+nob+1)
          call get_ion_coeff(hole_index,part_index,psi_det0,psi1,Priv%norm,&
                 Mol%vabsmoa,Mol%vabsmob,Priv%rate,&
                 Zion_coeff( ii+1 : ii+(noa+nob)*(noa+nob+1) ),ion_coeff,scratch)
          !write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,scratch(1:noa+nob)
          do i = 1,noa+nob
            do j =1,noa+nob
              Zion_coeff(ii+i*(noa+nob)+j) =  &
                  dsqrt(Priv%rate) * scratch(i) * Zion_coeff(ii+i*(noa+nob)+j)
            end do
          end do
          !if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
          !                                       itime,norm,(scratch(i),i=1,noa+nob)               
        end if
        call cpu_time(finish3)
        times(7) = finish3 - start3

        analysis : if ( mod(itime,outstep).eq.0 ) then                
            
          idata = int( itime/outstep)  
            
          call get_norm( Priv%norm, nstuse, psi )
          call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )

          !$OMP CRITICAL (IO_LOCK)

          if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
            call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                 hole_index,part_index,state_ip_index,ip_states, &
                 pop1,ion,ion_coeff,Priv%rate_aa,Priv%rate_ab,Priv%rate_ba,Priv%rate_bb, &
                 psi_det0,psi1,Priv%normV,Mol%vabsmoa,Mol%vabsmob,scratch,au2fs,Priv%rate_direct)
          else
            call pop_rate(Mol,jj,kk,hole_index,part_index,state_ip_index, &
                 pop1,ion,ion_coeff,Priv%rate_a,Priv%rate_b,psi_det0,psi1,Priv%normV, &
                 Mol%vabsmoa,Mol%vabsmob,scratch,Priv%rate_direct)
          end if 
          !: 1-RDM (opdm) should be stored in scratch now.

          !: Update nva
          call update_maxnva( Priv%nva95_MO, Priv%nva99_MO, Priv%nva95_NO, Priv%nva99_NO, &
            Priv%nva95_MO_sort, Priv%nva99_MO_sort, Priv%nva95_NO_sort, Priv%nva99_NO_sort, &
            Priv%nva95_debug, Priv%nva99_debug, Priv%rate_debug, Priv%rate_noabs, &
            scratch, Prop%U_NO_input, Mol%vabsmoa, .false. ) !: last arg verbosity


          call add_opdm_average( Prop%opdm_avg, scratch, Prop%opdm_avg_N )
          call add_opdm_average( Prop%opdm_avg_abs, scratch, &
                                 Prop%opdm_avg_N, .false., .true. )


          !: Should we rip this out? It never worked.
          ! Calculate complex density matrix and transition rates
          if ( write_orb_transitionrates ) then
            write(iout,*) "before make_complex_density  ", ithread, " ", itime 
            flush(iout)
            call make_complex_density(hole_index(:,1),part_index(:,1),psi_det0,noa,nva,nstates,Priv%density_complex)
            write(iout,*) "before make_transition_rate  ", ithread, " ", itime 
            flush(iout)
            call make_transition_rate((noa+nva), Mol%ham_mo, Mol%dipxmoa, Mol%dipymoa, Mol%dipzmoa, Priv%efieldx, &
                   Priv%efieldy, Priv%efieldz, Mol%vabsmoa, Priv%density_complex, Mol%orben, Priv%transition_rates)

            write( Priv%datafile, '(A,A,A,A,A,I0,A)') "matrices/MO_TransitionRate-e",&
              trim(Priv%emaxstr),"-d", trim(Priv%dirstr), ".", itime, ".bin"

            call write_dbin_safe(Priv%transition_rates, (noa+nva)*(noa+nva), Priv%datafile )
          end if
          !: 1-RDM stored in scratch now.

          call get_norm( Priv%norm,nstuse, psi )
          call get_expectation( nstuse, Priv%norm, psi, abp, Priv%rate) !: rate expectation value
          call get_expectation( nstuse, Priv%norm, psi, tdx, Priv%mux ) !: mux  expectation value
          call get_expectation( nstuse, Priv%norm, psi, tdy, Priv%muy ) !: muy  expectation value
          call get_expectation( nstuse, Priv%norm, psi, tdz, Priv%muz ) !: muz  expectation value
          !write(iout,*) " expectation Priv%rate", itime,Priv%norm,Priv%rate
          Priv%rate = -2.d0 * Priv%rate * Priv%norm**2

          if(datfile_enable) then
            call cpu_time(start2)
            call PropWriteData(Priv, psi, psi1, psi_det0, Zion_coeff,ion_coeff,&
                               scratch, itime, idata, pop1, ion)
            call cpu_time(finish2)
            Priv%plaincumtime = Priv%plaincumtime + (finish2-start2)
          end if

          if(h5inc_enable) then
            call cpu_time(start2)
            call write_h5_step(Priv, psi, psi1, psi_det0, Zion_coeff, ion_coeff,&
                               scratch, itime, idata, pop1, ion)
            call cpu_time(finish2)
            Priv%h5cumtime = Priv%h5cumtime + (finish2-start2)
          end if

          !$OMP END CRITICAL (IO_LOCK)

        end if analysis
         
      end do timestep_loop
      
      call cpu_time(Priv%finish2)

      if( Qci_save ) close(Priv%funit(1))
      close(Priv%funit(2))
      close(Priv%funit(3))
      close(Priv%funit(4))

      if( Qwrite_ion_coeff ) write(Priv%funit(5)) Zion_coeff
      if( Qread_ion_coeff .or. Qwrite_ion_coeff ) close(Priv%funit(5))

      if( Qmo_dens ) close(Priv%funit(6))

      !$OMP CRITICAL (IO_LOCK)
      if(h5inc_enable) then
        !write(iout,*) "Closing h5 emax group ", idir, iemax, Priv%dir_grpid, Priv%field_grpid
        call h5_close_group(Priv%field_grpid)
        !write(iout,*) "Done closing h5 emax group ", idir, iemax, Priv%dir_grpid, Priv%field_grpid
      end if
      !: record data at last timestep
      tdciresults( idir + (iemax-1)*ndir)%norm0 = Priv%norm**2
      tdciresults( idir + (iemax-1)*ndir)%dipx  = Priv%mux
      tdciresults( idir + (iemax-1)*ndir)%dipy  = Priv%muy
      tdciresults( idir + (iemax-1)*ndir)%dipz  = Priv%muz
      
      if(verbosity.gt.1) then
        write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
        write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           Priv%dirx1, Priv%diry1, Priv%dirz1, Priv%emax1
        if(h5inc_enable) then
          write(iout,"(12x,'h5 calls time:            ',f12.4,' s')") Priv%h5cumtime
        end if
        if(datfile_enable) then 
          write(iout,"(12x,'plain output calls time:  ',f12.4,' s')") Priv%plaincumtime
        end if
        write(iout,  "(12x,'propagation time:         ',f12.4,' s')") Priv%finish2 - Priv%start2
        write(iout,  "(12x,'final norm = ',f10.5)")  Priv%norm**2
        write(iout,  "(12x,'final rate = ',f18.10)") Priv%rate
      end if !: verbosity.gt.1

      !: debug times
      if(verbosity.gt.2) then
        do i=1,7
          write(iout, "('time ',i5,': ',f14.6,' s')") i, times(i) 
        end do
      end if

      !$OMP END CRITICAL (IO_LOCK)
     end do emax_loop
    call Priv%deconstruct()
    deallocate(Priv)

  end do dir_loop

  !$OMP END DO
  !$OMP END PARALLEL

  call write_dbin_safe(Prop%opdm_avg, (noa+nva)*(noa+nva), "opdm_avg.bin")
  call write_dbin_safe(Prop%opdm_avg_abs, (noa+nva)*(noa+nva), "opdm_avg_abs.bin")

  call generate_natural_orbitals( Prop%opdm_avg, Prop%U_NO, Prop%natorb_occ, .false. )
  call generate_natural_orbitals( Prop%opdm_avg_abs, Prop%U_NO_abs, Prop%natorb_occ_abs, .false. )
  call write_dbin_safe(Prop%U_NO, (noa+nva)*(noa+nva), "U_NO_out.bin")
  call write_dbin_safe(Prop%U_NO_abs, (noa+nva)*(noa+nva), "U_NO_abs_out.bin")



  if(h5inc_enable) then
    call write_h5metadata(Prop%datetime_start)
  end if
  call Prop%deconstruct()
  call write_header( 'trotter_linear','propagate','leave' )
  call cpu_time(finish)
  write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start  
  flush(iout)

end subroutine trotter_linear
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
! SUBROUTINE TROTTER_CIRCULAR
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!  
subroutine trotter_circular


  ! <C> propagation using circularly polarized lights.  Takes in fvect1 and fvect2
  ! <C> C(t+dt) = exp(-iHel dt/2)exp(-Vabs dt/2) * W1exp(-iE1(t) dt/2)W1* W2*exp(-iE2(t+dt/2)mu dt/2)*W2
  ! <C>           *W1exp(-iE1(t) dt/2)W1 * exp(-Vabs dt/2)exp(-iHel dt/2)*C(t)


  use omp_lib
  implicit none

  !: thread shared variables
  class(PropagationShared), allocatable :: Prop
  integer(8) :: nstuse2 
  complex(8) :: exphel(nstuse), psi0(nstuse)
  real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse) 
  real(8) :: hp2(nstuse*nstuse), tdvals2(nstuse)        
  

  !: thread private variables
  class(PropagationPrivate), allocatable :: Priv    
  !: Psi arrays must be allocated on the stack for performance
  complex(8) :: psi_j, psi_k, psi(nstuse), psi1(nstates)

  real(8) :: norm0
  real(8),allocatable :: pop0(:),pop1(:),ion(:)
  complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)

  integer(8) :: i, j, ii, jj, k, kk
  integer(8) :: itime, idir, iemax
  integer(8) :: ithread, idata
  complex(8) :: cdum

  real(8) :: start1, finish1, start2, finish2, start3, finish3
  real(8) :: times(100)

  !: lapack stuff
  integer(8) :: info1, info2, lscratch, liwork
  integer(8), allocatable :: iwork(:)
  real(8), allocatable    :: scratch(:)

  call write_header( 'trotter_circular','propagate','enter' )    
  call writeme_propagate( 'trot_cir', 'equation' ) 
  call cpu_time(start)

  nstuse2 = nstuse*nstuse

  !: allocation and exphel = exp(-iH*dt/2)
  call trotter_init(Prop, exphel, psi0, norm0, pop0, pop1, ion, psi_det0, ion_coeff, Zion_coeff, iwork, scratch)

  !: write psi0
  call writeme_propagate( 'trot_lin', 'psi0' )


  !: Start loop over directions

  !$OMP PARALLEL DEFAULT(NONE),&
  !$OMP PRIVATE(Priv, i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
  !$OMP ithread, cdum, &
  !$OMP norm0, lscratch, liwork, scratch, iwork, info1, info2, start1, finish1, &
  !$OMP start2, finish2, start3, finish3, &
  !$OMP psi_j, psi_k, psi, psi1, times, &
  !$OMP pop1, ion, ion_coeff, psi_det0,  &
  !$OMP hp1, hp2, tdvals1, tdvals2, Zion_coeff ),  &
  !$OMP SHARED( Mol, Prop, jobtype, nbasis, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
  !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
  !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
  !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
  !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
  !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
  !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
  !$OMP flag_ReadU_NO, write_orb_transitionrates, h5inc_enable, datfile_enable, verbosity )

  !$OMP DO
  dir_loop : do idir=1, ndir 

    ithread = omp_get_thread_num()
    call cpu_time( start1 ) 

    psi = 0 ; psi1 = 0
    tdvals1 = 0.d0 ; tdvals2 = 0.d0

    !write(iout, '(A, I5, L1)') "Priv allocated at start of thread ", ithread , allocated(Priv)
    allocate(Priv)
    call Priv%initialize()

    !: get directions stored in TDCItdciresults
    Priv%dirx1 = tdciresults(idir)%x1 ; Priv%dirx2 = tdciresults(idir)%x2
    Priv%diry1 = tdciresults(idir)%y1 ; Priv%diry2 = tdciresults(idir)%y2
    Priv%dirz1 = tdciresults(idir)%z1 ; Priv%dirz2 = tdciresults(idir)%z2

    if(h5inc_enable) then
      !$OMP CRITICAL (IO_LOCK)
      !: Create the /direction_idir group in the h5 file
      call init_h5dir(Priv, idir)
      !$OMP END CRITICAL (IO_LOCK)
    end if

    !: Form the transition dipole in (dirx1,diry1,dirz1) and (dirx2,diry2,dirz2) directions,
    hp1(:) = Priv%dirx1*tdx(1:nstuse2) + Priv%diry1*tdy(1:nstuse2) + Priv%dirz1*tdz(1:nstuse2)
    hp2(:) = Priv%dirx2*tdx(1:nstuse2) + Priv%diry2*tdy(1:nstuse2) + Priv%dirz2*tdz(1:nstuse2)

    !: diagonalize hp1 and hp2
    info1=10
    info2=10
    if( QeigenDC ) then
      lscratch = 1+6*nstuse+2*nstuse*nstuse
      liwork = 3+5*nstuse
      call dsyevd('v','u',nstuse, hp1, nstuse, tdvals1, scratch, lscratch, iwork, liwork, info1)
      call dsyevd('v','u',nstuse, hp2, nstuse, tdvals2, scratch, lscratch, iwork, liwork, info2)
    else
      lscratch = nstuse*nstuse
      call dsyev('v','u',nstuse, hp1, nstuse, tdvals1, scratch, lscratch, info1)
      call dsyev('v','u',nstuse, hp2, nstuse, tdvals2, scratch, lscratch, info2)
    end if
     
    !: hp2 = hp2 * hp1 ; hp2 = W2 * W1
    call dgemm('t','n',nstuse,nstuse,nstuse,1.d0,hp2,nstuse,hp1,nstuse,0.d0,scratch,nstuse)
    hp2 = scratch(1:nstuse*nstuse)
    !: hp1 = W1 * exp(-Vabs dt/2)
    call dgemm('t','n',nstuse,nstuse,nstuse,1.d0,hp1,nstuse,exp_abp,nstuse,0.d0,scratch,nstuse)
    hp1 = scratch(1:nstuse*nstuse)

    call cpu_time( finish1 )

    !: loop over intensities
    emax_loop : do iemax=1, nemax

      !: get emax ; emax1==emax2 in this version
      Priv%emax1 = tdciresults((iemax-1)*ndir+1)%fstrength1
      Priv%emax2 = tdciresults((iemax-1)*ndir+1)%fstrength2

      !$OMP CRITICAL (IO_LOCK)

      if(datfile_enable) then
        call cpu_time(start2)
        call PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, psi_det0, Zion_coeff, 0)
        call cpu_time(finish2)
        Priv%plaincumtime = Priv%plaincumtime + (finish2-start2)
      end if

      if(h5inc_enable) then
        call cpu_time(start2)
        call init_h5emax(Priv, iemax, idir)
        call cpu_time(finish2)
        Priv%h5cumtime = Priv%h5cumtime + (finish2-start2)
      end if

      if(verbosity.gt.1) then
        if(iemax.eq.1) write(iout,"(12x,'Init, TD diag, and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
        write(iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0)" ) idir, iemax, ithread
      end if
      if(verbosity.gt.2) then
        write( iout, "('  Maximum(tdvals1) = ',f12.4)") maxval(tdvals1)
        write( iout, "('  Minimum(tdvals1) = ',f12.4)") minval(tdvals1)
      end if
      flush( iout )
      !$OMP END CRITICAL (IO_LOCK)
        
      !: initialize psi
      call cpu_time(start1)
      psi = psi0
      if( read_state1(iemax).ne.0 ) then
        psi = dcmplx(0.d0,0.d0)
        write(iout,"(' *** Reset initial wavefunction ***')")
        if( read_state1(iemax).ne.0 ) then
          i = read_state1(iemax)
          psi(i) = read_coeff1(iemax)
          write(iout,"(' Initial coefficient for state',i4,' is ',2f12.8)") i,psi(i)
        end if
        if( read_state2(iemax).ne.0 ) then
          i = read_state2(iemax)
          psi(i) = read_coeff2(iemax)
          write(iout,"(' Initial coefficient for state',i4,' is ',2f12.8)") i,psi(i)
        end if
        call get_norm( norm0, nstuse, psi )
        psi = psi / norm0
        psi0 = psi
        flush(iout)
      end if
      if( Qread_ion_coeff ) psi = dcmplx(0.d0,0.d0)

      call cpu_time(finish1)
      times(1) = start1 - finish1 !: initialize psi time

      !: begin Looping over time
      call cpu_time( start2 )
      timestep_loop : do itime=1, nstep-1
           
        call cpu_time(start3)
        call PrivSetField(Priv, itime, iemax)

        !: How can we factor out this ion sampling stuff?
        if( itime.lt.ion_sample_start(iemax) ) psi = dcmplx(0.d0,0.d0)
        if( Qread_ion_coeff ) then
          if( itime.ge.ion_sample_start(iemax) .and. &
              itime.le.ion_sample_start(iemax)+ion_sample_width(iemax) ) then
            ii = (itime-1)*(noa+nob)*(noa+nob+1)
            !: tranform ion_coeff to CI basis and add to psi
            !:write(iout,"('R0',i5,i3,16F10.6)") itime,ion_sample_state(iemax)
            call get_ion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                  dt,cis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
          end if
          else
            if( itime.eq.ion_sample_start(iemax) ) psi = psi0
          end if
          call cpu_time(finish3)
          times(2) = finish3 - start3
           
          !: Propagation code below. TODO: add reference.
          !: psi contains C(t)
          !: psi = exp(-iHel dt/2 ) * C(t)
          do j = 1, nstuse
            psi(j) = exphel(j) * psi(j)
          end do
           
          !: psi contains exp(-iHel dt/2 ) * C(t)
          !: hp1 is W1 := 
          !: psi1 = W1 * exp(-Vabs dt/2) * psi
          psi1 = dcmplx( 0.d0, 0.d0 )
          do j = 1, nstuse
            jj = (j-1) * nstuse
            psi_j = psi(j)
            do k = 1 , nstuse
              psi1(k) = psi1(k) + hp1(jj+k) * psi_j
            end do
          end do

          !: psi1 contains W1 * exp(-Vabs dt/2) * exp(-iHel dt/2 ) * C(t)
          !: psi1 = exp( iE1(t+dt/2)*mu * dt/2) * psi1
          do j = 1, nstuse
            Priv%temp = Priv%temp1 * tdvals1(j)
            psi1(j) = dcmplx( dcos(Priv%temp),dsin(Priv%temp) ) * psi1(j)
          end do
           
          !: W2*W1 * psi
          psi = dcmplx( 0.d0, 0.d0 )
          do j = 1, nstuse
            jj = (j-1)*nstuse
            psi_j  = psi1(j)
            do k=1, nstuse
              psi(k) = psi(k) + hp2(jj+k) * psi_j
            end do
          end do
           
          !: exp( iE2(t+dt/2)*mu* dt ) * psi
          do j = 1, nstuse
            Priv%temp = Priv%temp2 * tdvals2(j)
            psi(j) = dcmplx(dcos(Priv%temp),dsin(Priv%temp)) * psi(j)
          end do


          !: W2*W1 * psi
          do j = 1, nstuse
            jj = nstuse*(j-1)
            psi_j = dcmplx( 0.d0, 0.d0 )
            do k=1, nstuse
              psi_j = psi_j + hp2(jj+k) * psi(k)
            end do
            psi1(j) = psi_j
          end do
           
          !: exp( iE1(t+dt/2)*mu*dt/2) * psi
          do j = 1, nstuse
            Priv%temp = Priv%temp1 * tdvals1(j)
            psi1(j) = dcmplx(dcos(Priv%temp),dsin(Priv%temp)) * psi1(j)
          end do
           
          !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W1 * psi
          do j = 1, nstuse
            jj = nstuse*(j-1)
            psi_j = dcmplx( 0.d0, 0.d0 )
            do k=1, nstuse
              psi_j = psi_j + hp1(jj+k) * psi1(k)
            end do
            psi(j) = exphel(j) * psi_j
          end do
           
          if( Qwrite_ion_coeff ) then
            call get_norm( Priv%norm, nstuse, psi )
            call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )
            ii = (itime-1)*(noa+nob)*(noa+nob+1)
            call get_ion_coeff(hole_index,part_index, &
              psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob, &
              Priv%rate, Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,scratch)
            !:write(iout,"('W0',i5,i7,16f12.8)") itime,ii,Priv%rate,scratch(1:noa+nob)
            do i = 1,noa+nob
              do j =1,noa+nob
                Zion_coeff(ii+i*(noa+nob)+j) =  &
                  dsqrt(Priv%rate) * scratch(i) * Zion_coeff(ii+i*(noa+nob)+j)
              end do
            end do
            !:if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
            !:  itime,norm,(scratch(i),i=1,noa+nob)
          end if
           
          analysis : if ( mod(itime,outstep).eq.0 ) then

            idata = int(itime/outstep)

            call get_norm( Priv%norm, nstuse, psi )
            call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )

            !$OMP CRITICAL (IO_LOCK)

            if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
              call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                 hole_index,part_index,state_ip_index,ip_states, &
                 pop1,ion,ion_coeff,Priv%rate_aa,Priv%rate_ab,Priv%rate_ba,Priv%rate_bb, &
                 psi_det0,psi1,Priv%normV,Mol%vabsmoa,Mol%vabsmob,scratch,au2fs,Priv%rate_direct)
            else
              call pop_rate(Mol,jj,kk,hole_index,part_index,state_ip_index, &
                 pop1,ion,ion_coeff,Priv%rate_a,Priv%rate_b,psi_det0,psi1,Priv%normV, &
                 Mol%vabsmoa,Mol%vabsmob,scratch,Priv%rate_direct)
            end if
            !: MO-based density matrix stored in scratch now.

            call get_norm( Priv%norm,nstuse, psi )
            call get_expectation( nstuse, Priv%norm, psi, abp, Priv%rate) !: rate expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdx, Priv%mux ) !: mux  expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdy, Priv%muy ) !: muy  expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdz, Priv%muz ) !: muz  expectation value
            !write(iout,*) " expectation rate", itime,Priv%rate
            Priv%rate = -2.d0 * Priv%rate * Priv%norm**2

            if(datfile_enable) then
              call cpu_time(start2)
              call PropWriteData(Priv, psi, psi1, psi_det0, Zion_coeff,ion_coeff,&
                                 scratch, itime, idata, pop1, ion)
              call cpu_time(finish2)
              Priv%plaincumtime = Priv%plaincumtime + (finish2-start2)
            end if

            if(h5inc_enable) then
              call cpu_time(start2)
              call write_h5_step(Priv, psi, psi1, psi_det0, Zion_coeff, ion_coeff,&
                                 scratch, itime, idata, pop1, ion)
              call cpu_time(finish2)
              Priv%h5cumtime = Priv%h5cumtime + (finish2-start2)
            end if

            !$OMP END CRITICAL (IO_LOCK)

          end if analysis

        end do timestep_loop
        call cpu_time(finish2)

        if( Qci_save ) close(Priv%funit(1))
        close(Priv%funit(2))
        close(Priv%funit(3))
        close(Priv%funit(4))
        if( Qwrite_ion_coeff ) write(Priv%funit(5)) Zion_coeff
        if( Qread_ion_coeff .or. Qwrite_ion_coeff ) close(Priv%funit(5))

        if( Qmo_dens ) close(Priv%funit(6))

        !$OMP CRITICAL (IO_LOCK)
        !: record data at last timestep
        tdciresults( idir + (iemax-1)*ndir)%norm0 = Priv%norm**2
        tdciresults( idir + (iemax-1)*ndir)%dipx  = Priv%mux
        tdciresults( idir + (iemax-1)*ndir)%dipy  = Priv%muy
        tdciresults( idir + (iemax-1)*ndir)%dipz  = Priv%muz

        if(verbosity.gt.1) then
          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread,idir,iemax
          write(iout,"(12x,'dir1 = (',f8.5,',',f8.5,',',f8.5,')  emax1 = ',f8.5,' au')") Priv%dirx1,Priv%diry1,Priv%dirz1,Priv%emax1
          write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") Priv%dirx2,Priv%diry2,Priv%dirz2,Priv%emax2
          write(iout,"(12x,'(iemax=1) TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
          write(iout,"(12x,'(iemax=1) LAPACK dysev TD diagonalization INFO=',i0)") info1
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2-start2
          write(iout,"(12x,'final norm = ',f10.5)")          Priv%norm**2
          write(iout,"(12x,'final rate = ',f10.5)")          Priv%rate
        end if !: verbosity.gt.1

        !: debug times
        if(verbosity.gt.2) then
          do i=1,7
            write(iout, "('time ',i5,': ',f14.6,' s')") i, times(i) 
          end do
        end if

        !$OMP END CRITICAL (IO_LOCK)
     end do emax_loop

    call Priv%deconstruct()
    deallocate(Priv)

  end do dir_loop

  !$OMP END DO
  !$OMP END PARALLEL

  if(h5inc_enable) then
    call write_h5metadata(Prop%datetime_start)
  end if
  call Prop%deconstruct()
  call write_header( 'trotter_circular','propagate','leave' )
  call cpu_time(finish)
  write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start  
  flush(iout)

end subroutine trotter_circular
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
! SUBROUTINE WRITEME_PROPAGATE
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
subroutine writeme_propagate(myroutine,option,pop,norb)

  implicit none

  character(*), intent(in) :: myroutine
  character(*), intent(in) :: option

  integer(8), optional, intent(in) :: norb
  real(8),    optional, intent(in) :: pop(norb)

  integer(8) :: i, iorb


  select case( trim(myroutine) )
  case( 'trot_lin' )
     trotter_linear : select case( trim(option) )
     case('equation')
        write(iout,'(A)') " C(t+dt) = exp[-iHel dt/2] * exp[-Vabs dt/2] * "
        write(iout,'(A)') "           W^T * exp[iE(t+dt/2)*mu dt] * W * "
        write(iout,'(A)') "           exp[-Vabs dt/2] * exp[-iHel dt/2] * C(t)"                 
        write(iout,"(' start propagation for ',i10,' timesteps')")        nstep
        write(iout,"('                   for ',i10,' field directions')") ndir
        write(iout,"('                   for ',i10,' field strengths')")  nemax    
     case('psi0')
        write(iout,'(A)') ' '
        write(iout,"(A)") " INITIALIZED STATE: "
        write(iout,"(A)",advance="no") '     |psi(0)> = '
        write(iout,100) ( '(',init_coeffs(i),')', init_states(i), i=1, init_states(0) )
     case('pop0')
        write(iout,'(A)') ' '
        write(iout,"(A)") " INITIALIZED MO POPULATION"
        if( unrestricted ) then
           write(iout,101) ( pop(iorb), iorb=1, noa )
           write(iout,300) ( pop(iorb), iorb=(noa+1), nrorb )
           write(iout,200) ( pop(iorb), iorb=(nrorb+1),(nrorb+nob) )
           write(iout,400) ( pop(iorb), iorb=(nrorb+nob+1),norb )
        else
           write(iout,500) ( pop(iorb), iorb=1, noa )
           write(iout,600) ( pop(iorb), iorb=(noa+1),nrorb )
        end if
     end select trotter_linear
  case( 'trot_circ' )
     trotter_circ : select case( option ) 
     case('equation')
        write(iout,'(A)') " C(t+dt) = exp[-iHel dt/2] * exp[-Vabs dt/2] * W1^T * exp[iE1(t+dt/2)*mu dt/2] * W1 *"
        write(iout,'(A)') "           W2^T * exp[iE2(t+dt/2)*mu dt] * W2 * "
        write(iout,'(A)') "           W1^T * exp[iE1(t+dt/2)*mu dt/2]  * W1 * exp[-Vabs dt/2] * exp[-iHel dt/2] * C(t)"
        
        write(iout,"(' start propagation for ',i10,' timesteps')")  nstep
        write(iout,"('                   for ',i10,' directions')") ndir
        write(iout,"('                   for ',i10,' strengths')")  nemax
     end select trotter_circ
  end select
  

  flush(iout)

100 format( 10(a1,f7.5,','f7.5,a1,'|',i0,'>  ') )
101 format( '  occ_a:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
300 format( '  vir_a:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
200 format( '  occ_b:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
400 format( '  vir_b:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
500 format( '    occ:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
600 format( '    vir:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )


end subroutine writeme_propagate
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

!: Accepts a one-particle reduced density matrix, and adds it to the averaged
!:  matrix. This is run on every step we print to outfiles, 
!:  so opdm_avg_N must be nstep/outstep to ensure normalization.
subroutine add_opdm_average( opdm_avg, opdm, opdm_avg_N, verbosity, useabs )
  
  implicit none

  integer(8), intent(in) :: opdm_avg_N
  real(8), intent(inout) :: opdm_avg(norb*norb)
  real(8), intent(in)    :: opdm(norb*norb)
  logical, intent(in), optional :: verbosity
  logical, intent(in), optional :: useabs

  integer(8) :: i, j, ndim
  logical :: verbose, useabs_

  !: default to nonverbose unless optional argument is true
  verbose = .false.
  if (present(verbosity)) verbose = verbosity
  useabs_ = .false.
  if (present(useabs)) useabs_ = useabs


  !: LATER ADD CASES FOR UHF
  ndim = noa+nvb

  if (useabs_) then
    do i = 1,ndim*ndim
      opdm_avg(i) = opdm_avg(i) + abs(opdm(i)/real(opdm_avg_N))
    end do
  else
    do i = 1,ndim*ndim
      opdm_avg(i) = opdm_avg(i) + opdm(i)/real(opdm_avg_N)
    end do
  end if    


end subroutine add_opdm_average

!: Diagonalize opdm_avg for natural orbitals.
!: Output: Transformation matrix from MO -> NO.
subroutine generate_natural_orbitals( opdm, U_NO, evals, verbosity )

  implicit none
  real(8), intent(in) :: opdm(norb*norb)
  real(8), intent(out) :: U_NO(norb*norb)
  real(8), intent(out) :: evals(norb)
  logical, intent(in), optional :: verbosity

  integer :: info_, lwork, liwork
  integer, allocatable :: iwork(:)
  real(8), allocatable ::  work(:)

  integer(8) :: ndim, i, j, idx1, idx2
  real(8) :: temp

  logical :: verbose

  !: default to nonverbose unless optional argument is true
  verbose = .false.
  if (present(verbosity)) verbose = verbosity

  !: LATER ADD CASES FOR UHF
  ndim = noa+nvb

  !: Write opdm
  if (verbose) then
    write(iout, '(A)') 'OPDM at start of generate_natural_orbitals:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, '(F8.4)', advance='no') opdm((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Copy opdm to U_NO
  call dcopy(ndim*ndim, opdm, 1, U_NO, 1)

  !: Symmetrize
  do i = 1,ndim
    do j = 1,i-1
      idx1 = (i-1)*ndim+j
      idx2 = (j-1)*ndim+i
      temp = 0.5d0 * (U_NO(idx1) + U_NO(idx2))
      U_NO(idx1) = temp
      U_NO(idx2) = temp
    end do
  end do
  
  !: workspace query to determine size for work arrays
  lwork = -1
  liwork = -1
  allocate(work(1))
  allocate(iwork(1))
  call dsyevd('v','u', ndim, U_NO, ndim, evals, work, lwork, iwork, liwork, info_)
  lwork = nint(work(1))
  liwork = iwork(1)
  deallocate(work, iwork)

  allocate(work(lwork))
  allocate(iwork(liwork))

  !: Diagonalize!
  call dsyevd('v','u', ndim, U_NO, ndim, evals, work, lwork, iwork, liwork, info_)

  write(iout,'(A)') 'After generate_natural_orbitals'
  if (info_) then
    write(iout, *) 'Natural Orb Diag error!!: ', info_
  end if

  !: Write U_NO
  if (verbose) then
    write(iout, '(A)') 'Natural Orb Occ unsorted: '
    do i=1,ndim
      write(iout, '(F12.10)') evals(i)
    end do
    write(iout, '(A)') 'U_NO unsorted:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, '(F8.4)', advance='no') U_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if


  deallocate(work, iwork)

  !: abs() all eigenvalues:
  !:  In standard case all eigenvalues are positive
  !:  In abs case, some small elements maybe negative, and we want to sort by magnitude.
  do i=1,ndim
    evals(i) = abs(evals(i))
  end do
  !: dsyevd arranges eigenpairs so that evals are ascending, we want descending.
  call sort_eigenpairs_desc(evals, U_NO, ndim)


  !: Write U_NO
  if (verbose) then
    write(iout, '(A)') 'Natural Orb Occ sorted: '
    do i=1,ndim
      write(iout, '(F12.10)') evals(i)
    end do
    write(iout, '(A)') 'U_NO sorted:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, '(F8.4)', advance='no') U_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if


end subroutine generate_natural_orbitals



subroutine update_maxnva( nva95_MO, nva99_MO, nva95_NO, nva99_NO, &
                          nva95_MO_sort, nva99_MO_sort, nva95_NO_sort, nva99_NO_sort, &
                          nva95_debug, nva99_debug, rate_debug, rate_noabs, &
                          opdm, U_NO, Vabs, verbosity)
  implicit none
  integer(8), intent(inout) :: nva95_MO, nva99_MO
  integer(8), intent(inout) :: nva95_NO, nva99_NO
  integer(8), intent(inout) :: nva95_MO_sort, nva99_MO_sort
  integer(8), intent(inout) :: nva95_NO_sort, nva99_NO_sort
  integer(8), intent(inout) :: nva95_debug, nva99_debug 
  real(8), intent(inout) :: rate_debug, rate_noabs  !, rate_density
  real(8), intent(in) :: opdm(:) !: One-particle density matrix in MO basis
  real(8), intent(in) :: U_NO(:) !: Transformation matrix from MO -> NO.
  !: Interesting, vabsmo is a 2D array (:,:), so I can't do Vabs(:) here.
  !:  but I can do Vabs(norb*norb)...
  real(8), intent(in) :: Vabs(nrorb*nrorb) !: <n|Vabs|m> in MO basis
  logical, intent(in), optional :: verbosity

  logical :: verbose
  integer(8) :: i, j, k, ndim, idx
  integer(8) :: nva95, nva99
  real(8) :: rate, tmp
  real(8) :: tmprate(nrorb*nrorb) !: way too big but doesnt matter
  real(8), allocatable :: opdm_NO(:), Vabs_NO(:)
  real(8), allocatable :: temp(:), temp2(:)
  real(8), allocatable :: sort_vals(:), vals_sorted(:)
  integer(8), allocatable :: sort_key(:)


  !: default to nonverbose unless optional argument is true
  verbose = .false.
  if (present(verbosity)) verbose = verbosity

  ndim = noa + nva
  allocate( temp(ndim*ndim) )
  allocate( temp2(ndim*ndim) )
  allocate( opdm_NO(ndim*ndim) )
  allocate( Vabs_NO(ndim*ndim) )
  allocate( sort_vals(ndim) )
  allocate( vals_sorted(ndim) )
  allocate( sort_key(ndim) )
  temp = 0.d0
  temp2 = 0.d0
  opdm_NO = 0.d0
  Vabs_NO = 0.d0
  sort_vals = 0.d0
  sort_key = 0
  
  !: Transform Vabs to NO basis
  !: U_NO^T * Vabs * U_NO = Vabs_NO
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, Vabs, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, Vabs_NO, ndim)

  !: Transform opdm to NO basis
  !: U_NO^T * opdm * U_NO = opdm_NO
  temp = 0.d0
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, opdm, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, opdm_NO, ndim)

  !write(iout, "(A)"), "Ensuring matrices are symmetric..."
  !if ( .not. matrix_is_symmetric( opdm, ndim ) ) then
  !  write(iout, "(A)"), "WARNING: opdm is not symmetric! (delta=1E-4)"
  !end if
  !if ( .not. matrix_is_symmetric( Vabs, ndim ) ) then
  !  write(iout, "(A)"), "WARNING: Vabs is not symmetric! (delta=1E-4)"
  !end if
  !if ( .not. matrix_is_symmetric( opdm_NO, ndim ) ) then
  !  write(iout, "(A)"), "WARNING: opdm_NO is not symmetric! (delta=1E-4)"
  !end if
  !if ( .not. matrix_is_symmetric( Vabs_NO, ndim ) ) then
  !  write(iout, "(A)"), "WARNING: Vabs_NO is not symmetric! (delta=1E-4)"
  !end if

  !: For NOs (UNSORTED)
  rate = 0.d0
  nva95_NO = 0
  nva99_NO = 0
  do i=1,ndim
    do j=1,ndim
      !: 2*Tr(opdm_NO * Vabs_NO) componentwise
      rate = rate + abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j))
    end do
  end do
  do k=1,ndim
    tmp = 0.d0
    do i=1,k
      do j=1,k
        tmp = tmp + abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j))
      end do
    end do
    if (tmp/rate .lt. 0.95d0) nva95_NO = k
    if (tmp/rate .lt. 0.99d0) nva99_NO = k
  end do


  !: For MOs ( DEBUG, UNSORTED, ABS(OPDM*VABS) )
  rate = 0.d0
  tmp = 0.d0
  nva95_debug = 0
  nva99_debug = 0
  tmprate = 0.d0
  tmprate(1) = abs(opdm(1)*Vabs(1))
  do i = 2,noa+nva
    tmprate(i) = tmprate(i-1) + abs(opdm(i+(i-1)*ndim)*Vabs((i-1)*ndim+i))
    do j = 1,i-1
      tmprate(i) = tmprate(i) + abs(opdm(i+(j-1)*ndim)*Vabs((i-1)*ndim+j))
      tmprate(i) = tmprate(i) + abs(opdm(j+(i-1)*ndim)*Vabs((j-1)*ndim+i))
    end do
  end do
  do i = 1,nva
    if (tmprate(i).lt.0.95d0*tmprate(noa+nva)) nva95_debug = i
    if (tmprate(i).lt.0.99d0*tmprate(noa+nva)) nva99_debug = i
  end do
  rate_debug = tmprate(noa+nva)

  !: For MOs (no abs())
  rate = 0.d0
  temp = 0.d0
  rate_noabs = 0.d0
  do i=1,ndim*ndim
    rate_noabs = rate_noabs + opdm(i)*Vabs(i)
  end do

  !: For MOs (UNSORTED)
  rate = 0.0d0
  tmp = 0.0d0
  nva95_MO = 0
  nva99_MO = 0
  do i=1,ndim
    do j=1,ndim
      !: 2*Tr(opdm * Vabs) componentwise
      rate = rate + abs(opdm((i-1)*ndim+j) * Vabs((i-1)*ndim+j))
    end do
  end do
  do k=1,ndim
    tmp = 0.d0
    do i=1,k
      do j=1,k
        tmp = tmp + abs(opdm((i-1)*ndim+j) * Vabs((i-1)*ndim+j))
      end do
    end do
    if (tmp/rate .lt. 0.95d0) nva95_MO = k
    if (tmp/rate .lt. 0.99d0) nva99_MO = k
  end do

  
  !: For NOs (SORTED)
  rate = 0.0d0
  tmp = 0.0d0
  nva95_NO_sort = 0
  nva99_NO_sort = 0
  sort_vals = 0.d0
  vals_sorted = 0.d0
  sort_key = 0

  !: Contribution from Orb i approximated by:
  !:  Sum_j{ opdm(i,j)*Vabs(i,j) + opdm(j,i)*Vabs(j,i) }
  do i=1,ndim
    do j=1,ndim
      sort_vals(i) = sort_vals(i) + abs(opdm_NO((i-1)*ndim+j)*Vabs_NO((i-1)*ndim+j))
      !sort_vals(i) = sort_vals(i) + abs(opdm_NO((j-1)*ndim+i)*Vabs_NO((j-1)*ndim+i))
    end do
  end do
  call quicksort_descending(sort_vals, vals_sorted, sort_key)

  if (verbose) then 
    write(iout, *) "NO Contribution Analysis"
    write(iout, *) "i, i_sort, rate, rate_sort"
    write(iout, *) "=========================="
    do i=1,ndim
      write(iout, "(i4, i4, f12.8, f12.8)") i, sort_key(i), sort_vals(i), vals_sorted(i)
    end do
  end if

  do i=1,ndim
    do j=1,ndim
      !: 2*Tr(opdm * Vabs) componentwise
      rate = rate + abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j))
    end do
  end do
  do k=1,ndim
    tmp = 0.d0
    do i=1,k
      do j=1,k
        !: index (i-1)*ndim+j, but i,j go through sort_key
        idx = (sort_key(i)-1)*ndim+sort_key(j)
        tmp = tmp + abs(opdm_NO(idx) * Vabs_NO(idx))
      end do
    end do
    if (tmp/rate .lt. 0.95d0) nva95_NO_sort = k
    if (tmp/rate .lt. 0.99d0) nva99_NO_sort = k
  end do


  !: For MOs (SORTED)
  rate = 0.0d0
  tmp = 0.0d0
  nva95_MO_sort = 0
  nva99_MO_sort = 0
  sort_vals = 0.d0
  vals_sorted = 0.d0
  sort_key = 0

  !: Contribution from Orb i approximated by:
  !:  Sum_j{ opdm(i,j)*Vabs(i,j) + opdm(j,i)*Vabs(j,i) }
  do i=1,ndim
    do j=1,ndim
      sort_vals(i) = sort_vals(i) + abs(opdm((i-1)*ndim+j)*Vabs((i-1)*ndim+j))
      !sort_vals(i) = sort_vals(i) + abs(opdm((j-1)*ndim+i)*Vabs((j-1)*ndim+i))
    end do
  end do
  call quicksort_descending(sort_vals, vals_sorted, sort_key)

  if (verbose) then
    write(iout, *) "MO Contribution Analysis"
    write(iout, *) "i, i_sort, rate, rate_sort"
    write(iout, *) "=========================="
    do i=1,ndim
      write(iout, "(i4, i4, f12.8, f12.8)") i, sort_key(i), sort_vals(i), vals_sorted(i)
    end do
  end if

  do i=1,ndim
    do j=1,ndim
      !: 2*Tr(opdm * Vabs) componentwise
      rate = rate + abs(opdm((i-1)*ndim+j) * Vabs((i-1)*ndim+j))
    end do
  end do
  do k=1,ndim
    tmp = 0.d0
    do i=1,k
      do j=1,k
        !: index (i-1)*ndim+j, but i,j go through sort_key
        idx = (sort_key(i)-1)*ndim+sort_key(j)
        tmp = tmp + abs(opdm(idx) * Vabs(idx))
      end do
    end do
    if (tmp/rate .lt. 0.95d0) nva95_MO_sort = k
    if (tmp/rate .lt. 0.99d0) nva99_MO_sort = k
  end do

  if (verbose) then
    write(iout,"('nva99_MO=',i5,' , nva99_MO_sort=',i5,' , nva99_NO=',i5,' , nva99_NO_sort=',i5)") nva99_MO, nva99_MO_sort, nva99_NO, nva99_NO_sort
  end if

  deallocate( temp )
  deallocate( temp2 )
  deallocate( opdm_NO )
  deallocate( Vabs_NO )

end subroutine update_maxnva


function matrix_is_symmetric(A, n) result(output)
  implicit none
  integer(8), intent(in) :: n
  real(8), intent(in) :: A(n*n)
  logical, intent(out) :: output

  integer(8) :: i,j
  real(8) :: delta

  output = .true.
  do i=1,n
    do j=1,n
      if (i .ne. j) then
        delta = abs( A((i-1)+j) - A((j-1)+i) )
        if ( delta .ge. 0.0001 ) then
          output = .false.
        end if
      end if
    end do
  end do 
end function matrix_is_symmetric


subroutine write_density_bin(filename, time, rate, noa, nva, density, vabs, funit)
  implicit none
  character(len=*), intent(in) :: filename
  real(8), intent(in)    :: time, rate
  integer(8), intent(in)    :: noa, nva
  real(8), intent(in)    :: density(:)
  real(8), allocatable, intent(in)    :: vabs(:)
  integer(8), intent(in), optional :: funit

  integer(8) :: i,j, idx, ndim, ndim2
  real(8) :: temp_density( (noa+nva)*(noa+nva) )
  real(8) :: temp, temptrace
  logical :: isSymm


  ndim = (noa+nva)
  ndim2 = ndim*ndim
  
  !: Prepare weighted density difference
  temp_density = density(:ndim2)

  if (present(funit)) then
    call write_dbin_safe(temp_density, nrorb*nrorb, trim(filename), funit)
    !call write_dbin(temp_density, nrorb*nrorb, trim(filename), funit)
  else
    call write_dbin_safe(temp_density, nrorb*nrorb, trim(filename))
    !call write_dbin(temp_density, nrorb*nrorb, trim(filename))
  end if

end subroutine write_density_bin

subroutine write_density_difference_bin(filename, time, rate, noa, nva, density, vabs, funit)
  implicit none
  character(len=*), intent(in) :: filename
  real(8), intent(in)    :: time, rate
  integer(8), intent(in)    :: noa, nva
  real(8), intent(in)    :: density(:)
  real(8), allocatable, intent(in)    :: vabs(:)
  integer(8), intent(in), optional :: funit

  integer(8) :: i,j, idx, ndim, ndim2
  real(8) :: temp_density( (noa+nva)*(noa+nva) )
  real(8) :: temp, temptrace
  logical :: isSymm

  ndim = (noa+nva)
  ndim2 = ndim*ndim

  temp_density = 0.d0
  do i=1, ndim2
    temp_density(i) = density(i)
  end do

  do i=1, noa !: Density difference from HF
    temp_density(i+(i-1)*ndim) = temp_density(i+(i-1)*ndim) - 2.d0
  end do

  if (present(funit)) then
    call write_dbin_safe(temp_density, ndim2, trim(filename), funit)
    !call write_dbin(temp_density, ndim2, trim(filename), funit)
  else
    call write_dbin_safe(temp_density, ndim2, trim(filename))
    !call write_dbin(temp_density, ndim2, trim(filename))
  end if

end subroutine write_density_difference_bin


subroutine write_vdens_diff_bin(filename, time, rate, noa, nva, density, vabs, funit)
  implicit none
  character(len=*), intent(in) :: filename
  real(8), intent(in)    :: time, rate
  integer(8), intent(in)    :: noa, nva
  real(8), intent(in)    :: density(:)
  real(8), allocatable, intent(in)    :: vabs(:,:)
  integer(8), intent(in), optional :: funit

  integer(8) :: i,j, idx, ndim, ndim2
  real(8) :: temp_density( (noa+nva)*(noa+nva) )
  real(8) :: temp, temptrace
  logical :: isSymm

  ndim = (noa+nva)
  ndim2 = ndim*ndim

  temp_density = 0.d0
  do i=1, ndim2
    temp_density(i) = density(i)
  end do

  do i=1, noa !: Density difference from HF
    temp_density(i+(i-1)*ndim) = temp_density(i+(i-1)*ndim) - 2.d0
  end do

  do i=1,ndim
    do j=1,ndim
      temp_density((i-1)*ndim+j) = temp_density((i-1)*ndim+j) * vabs(i,j)
    end do
  end do

  if (present(funit)) then
    call write_dbin_safe(temp_density, ndim2, trim(filename), funit)
    !call write_dbin(temp_density, ndim2, trim(filename), funit)
  else
    call write_dbin_safe(temp_density, ndim2, trim(filename))
    !call write_dbin(temp_density, ndim2, trim(filename))
  end if

end subroutine write_vdens_diff_bin




subroutine write_density_difference(funit, time, rate, noa, nva, density, vabs)
  implicit none
  integer(8), intent(in)              :: funit
  real(8), intent(in)                 :: time, rate
  integer(8), intent(in)              :: noa, nva
  real(8), intent(in)                 :: density(:)
  real(8), allocatable, intent(in)    :: vabs(:,:)

  integer(8) :: i,j, idx, ndim, ndim2
  real(8) :: temp_density( (noa+nva)*(noa+nva) )
  real(8) :: temp, temptrace
  logical :: isSymm

  ndim = (noa+nva)
  ndim2 = ndim*ndim

  !: Prepare weighted density difference
  temp_density = density(:ndim2)
  do i=1, noa !: Density difference from HF
    temp_density(i+(i-1)*ndim) = temp_density(i+(i-1)*ndim) - 2.d0
  end do

  !: Weight density by vabs
  !do i=1,ndim2
  !  temp_density(i) = temp_density(i) * vabs(i)
  !end do
  do i=1,ndim
    do j=1,ndim
      temp_density(j+(i-1)*ndim) = temp_density(j+(i-1)*ndim) * vabs(i,j)
    end do
  end do

  !: Re-normalize density matrix
  !temptrace =  sum([(temp_density(i+(i-1)*ndim), i=1,ndim)])
  !temp_density = temp_density / temptrace

  !: Scale so max value is 2.0
  !temptrace = 0.5*maxval(abs(temp_density))
  !temptrace = 0.0020426 !: scaling in ch3br at end of 045 field
  temptrace = 1.00
  temp_density = temp_density/temptrace

  !write(iout, "('density scaling factor: ', f16.10)") temptrace

  !write(funit,"(f13.9,f16.10,f16.10, L1)") time, rate, temptrace, isSymm
  write(funit,"(f13.9,f16.10)") time, rate
  do i=1, ndim !: Diagonal
    idx = i+(i-1)*ndim
    if ( abs(density(idx)).gt.1.d-9 ) then
      write(funit,"(i5,i5,1x,f13.10)") &
            i, i, temp_density(idx)
    end if
    do j =1, i-1 !: Off-Diagonals
      idx = i+(j-1)*ndim
      if ( abs(density(idx)).gt.1.d-9 ) then 
        write(funit,"(i5,i5,1x,f13.10)") &
              i, j, temp_density(idx)
      end if
    end do
  end do
  write(funit,"(i5,i5,1x,f13.10)") 0,0,0.d0

end subroutine write_density_difference


subroutine write_vabsmo_table(nrorb, Vabs_MO, orben)
  implicit none
  integer(8), intent(in) :: nrorb
  real(8), intent(in) :: Vabs_MO(nrorb*nrorb), orben(2*nrorb)

  integer(8), parameter :: iout = 42
  integer(8) :: i,j
  real(8) :: vabs_sum

  write(iout,*) "Vabs MO summary:"
  write(iout,*) "  i,     Energy (eV), Vabs(i,i), sum(Vabs(i)) "
  write(iout,*) "====================================="
  do i=1,nrorb
    vabs_sum = 0.0
    do j=1,nrorb
      vabs_sum = vabs_sum + Vabs_MO((i-1)*nrorb+j)
    end do
    write(iout,'(I4,A2,F14.4,A2,F8.4,A2,F8.4)') i, ", ", orben(i)*au2ev, ", ", Vabs_MO((i-1)*nrorb+i), ", ", vabs_sum
    !write(iout,*) i, ",", orben(i), ",", Vabs_MO((i-1)*nrorb+i), ",", vabs_sum
  end do
  write(iout,*) ""


end subroutine write_vabsmo_table



!: Duplicate code for trotter_linear and trotter_circular
!: before the propagation loop
subroutine trotter_init(Prop, exphel, psi0, norm0, pop0, pop1, ion, psi_det0, ion_coeff, Zion_coeff, iwork, scratch, cwork)
  implicit none
  class(PropagationShared), allocatable, intent(inout) :: Prop
  complex(8), intent(inout) :: exphel(nstuse), psi0(nstuse)

  real(8), intent(inout) :: norm0
  real(8),allocatable, intent(inout) :: pop0(:),pop1(:),ion(:)
  complex(8), allocatable, intent(inout) :: psi_det0(:),ion_coeff(:),Zion_coeff(:)
  integer(8), allocatable, intent(inout)  :: iwork(:)
  real(8), allocatable, intent(inout)     :: scratch(:)
  complex(8), allocatable, optional, intent(inout)     :: cwork(:)

  real(8) :: temp
  integer(8) :: i


  !: initialize psi0
  psi0 = dcmplx(0.d0,0.d0)
  do i=1, init_states(0)
    psi0( init_states(i) ) = init_coeffs(i)
  end do

  !: normalize
  call get_norm( norm0, nstuse, psi0 )
  psi0 = psi0 / norm0
  !: Initialize shared variables
  allocate(Prop)
  call Prop%initialize()

  !: get initial population
  allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
  i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
  allocate(ion_coeff(i))
  if( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
    allocate(Zion_coeff(nstep*ip_states*(ip_states+1)))
  else
    allocate(Zion_coeff(ip_states*(ip_states+1)))
    allocate(Zproj_ion(ip_states*(ip_states+1)))
  end if

  norm0 = 1.d0


  if( (trim(jobtype).eq.flag_soc) ) then
    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
  else !: linear trotter 
    call get_psid( nstuse, nstates, cis_vec, norm0, psi0, psi_det0 )
  end if

  select case ( trim(jobtype) )
  case( flag_cis )
    if ( unrestricted ) then
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
    else
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
    end if
  case( flag_soc )
    if ( unrestricted ) then
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
    else
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
    end if
  case( flag_tda )
    if ( unrestricted ) then
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
    else
      call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
    end if
  case( flag_ip)
    call pop_ip( hole_index, noa,nob,norb,nstates,nva,nvb, part_index(:,1), pop0, psi_det0 )
  case( flag_socip)
    call pop_ip( hole_index, noa,nob,norb,nstates,nva,nvb, part_index(:,1), pop0, psi_det0 )
  end select
  call writeme_propagate( 'trot_lin','pop0', pop0, norb )

  deallocate( pop0 )


  if( QeigenDC ) then
    allocate( iwork(3+5*nstuse) )
    allocate( scratch(1+8*nstuse+2*nstuse*nstuse) )
    if(trim(jobtype).eq.flag_soc) then !: Ztrotter
      if(linear) then
        allocate( cwork(1+8*nstuse+2*nstuse*nstuse) ) !: Ztrotter_linear
      else
        allocate( cwork(nstuse*nstuse+2*nstuse) ) !: Ztrotter_circular
      end if
    end if
  else !: .not. QeigenDC
    if((.not.linear).and.(trim(jobtype).eq.flag_soc)) then !: Ztrotter_circular
      allocate( iwork(2) )
      i = max(3*nstuse-2,(ip_states)*(ip_states))
      allocate( scratch(i) ) !: weird that rwork (scratch) and cwork are swapped for circular. TODO: double-check array bounds in Ztrotter_circular
      allocate( cwork(nstuse*nstuse) )
    else !: trotter_circular, trotter_linear, Ztrotter_linear
      allocate( iwork(2) )
      if(trim(jobtype).eq.flag_soc) then !: Ztrotter_linear
        i = max(3*nstuse,(ip_states)*(ip_states))
        allocate( cwork(i) )
      end if
      allocate( scratch(nstuse*nstuse) )
    end if
  end if

  !: exphel = exp(-iH*dt/2)
  do i=1, nstuse
    temp = -0.5d0 * cis_eig(i) * dt
    exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
  end do

  call write_vabsmo_table(nrorb, Mol%vabsmoa, Mol%orben )
  ! cis_vec stores hamiltonian
  if ( write_orb_transitionrates ) then
    allocate(Mol%ham_mo(nrorb*nrorb))
    !call ham_cis_to_mo(hole_index,part_index,psi_det0, noa,nva,nstates,nrorb,cis_vec,Mol%ham_mo,nob,nvb)
    call ham_cis_to_mo(hole_index(:,1),part_index(:,1),psi_det0, noa,nva,nstates,nrorb,cis_vec,Mol%ham_mo,nob,nvb)
  end if

end subroutine trotter_init

end module propagate

