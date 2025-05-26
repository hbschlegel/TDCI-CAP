module variables_global


  !: contains main variables 


  use variables_units
  use variables_setup
  use variables_control

  
  implicit none



  !: default file parameters
  integer(8),   parameter :: iout = 42            !: unit number for outfile 'tdci.log'
  character(5), parameter :: inputfile  = 'input' !: see make_input
  character(5), parameter :: fieldfile  = 'enput' !: see make_enput
 
  !: pretty write
  character(48), parameter :: divide      = ' _______________________________________________'
  character(48), parameter :: smalldivide = ' ______'
  
  !: jobtype flags
  character(3), parameter :: flag_cis  = 'cis'
  character(3), parameter :: flag_soc  = 'soc'   !: for cis with soc
  character(5), parameter :: flag_socip= 'socip' !: for cisd-ip with soc
  character(3), parameter :: flag_tda  = 'tda'
  character(2), parameter :: flag_ip   = 'ip'  
  character(4), parameter :: flag_cisd = 'cisd'

  
  !: MolInfo class instance, contains everything from variables_setup
  class(MolInfo), allocatable :: Mol
  


  !: data type for recording tdciresults
  type tdcidat 
     real(8) :: fstrength0, fstrength1, fstrength2
     real(8) :: theta0, phi0
     real(8) :: x0, y0, z0
     real(8) :: theta1, phi1, theta2, phi2
     real(8) :: x1,y1,z1, x2,y2,z2
     real(8) :: norm0
     real(8) :: dipx, dipy, dipz
  end type tdcidat

  
  integer(8)     :: tot_use_mem
  
  
  !: /FIELD/ variables
  integer(8) ::   ncyc     !: total number of field cycles
  real(8) ::      ellipt   !: ellipticity for elliptically polarized light
  real(8) ::      euler    !: the other Euler angle for circularly polarized light
  character(4) :: envelope !: flag to set the shape of pulse
  character(5) :: dirform  !: directions are specified in 'cart' or 'polar'  
  
  !: /FIELD_UNITS/ variables
  character(10) ::  &
       omega_units     !: 'au' or 'nm' 
       !angle_units    !: 'rad' or 'deg'
  
  !: /FIELD_strengths/ variables
  integer(8) :: nemax         !: total number of different field strength
  real(8)    :: read_emax(40) !: change if want to have more than 40 different field strengths 
  integer(8) :: read_state1(40), read_state2(40) !: initial states to superimpose
  complex(8) :: read_coeff1(40), read_coeff2(40) !: initial coefficients for superimposed states
  integer(8) :: read_shift(40) !: shift pulse by n timesteps
  integer(8) :: ion_sample_start(40),ion_sample_width(40),ion_sample_state(40)
             !: sequential double ionization sampling of ion_coeff read in ip or socip simulations

  !: /FIELD_directions/ variables
  integer(8) :: ndir !: number of polarization direction per field strength
  real(8) ::            &       
       read_theta(100), & !: stores theta.  Change if have more than 100 directions for each field strength
       read_phi(100),   & !: stores phi.    Change if have more than 100 directions for each field strength 
       read_x(100), read_y(100), read_z(100) !: or could read the directions in x,y,z directions.  not automated
        

 !: /SYSTEM/ variables
  integer(8) ::    &
       nstep  ,    & !: total number of timesteps 
       outstep,    & !: write to outfiles every outstep
       nactive,    & !: number of active orbitals in cisd-ip calculations
       nvirtual      !: number of virtual orbitals in calculations if less thaan nva
  real(8)    ::    &
       omega ,     & !: field-frequency in au units
       phase         !: phase of the field in radians
  real(8)    ::    &
       eigmax    , & !: ignore all cis states with energy above eigmax threshold
       ionization, & !: used to setup eigmax
       dt        , & !: size of timestep in au
       heuristic , & !: for harmonic generation; J. Chem. Phys. 131, 114304 (2009); https://doi.org/10.1063/1.3218847
       socfac    , & !: scaling factor for spin-orbit coupling (default socfac = 1.d0)
       socfacz   , & !: scaling factor for z component of spin-orbit coupling (default socfacz = 1.d0)
       ffieldx   , & !: finite field in the x direction added to Hamiltonian before diagonalizing (default 0.d0)
       ffieldy   , & !: finite field in the y direction added to Hamiltonian before diagonalizing (default 0.d0)
       ffieldz       !: finite field in the z direction added to Hamiltonian before diagonalizing (default 0.d0)
  
  
  !: /SYSTEM_units/ variables
  character(10) ::    &
       eigmax_units , & !: 'au' or 'eV'
       dt_units         !: 'au' or 'fs' or 'as'       

  
  !: /DYNAMICS/ variables
  integer(8) :: init_states(0:10) !: read in initial states.  
                                  !: init_states(0) = # states in superposition at t=0
                                  !: if starting with superposition of more than 10 states, change array size
  complex(8) :: init_coeffs(1:10) !: values CI coefficients at t=0
  logical    :: restart
  

  !: /InOutFILES/ variables
  character(1000) ::    & 
       tdcidatfile,     & !: defaults to 'TDCI.dat' if not specified.  Will die if 'TDCI.dat' not found 
       outputfile,      & !: defaults to 'OUTPUT' if not specified
       restartbinfile     !: restart bin file
  
  !: /hdf5/ variables
  logical :: h5inc_density


  !: data from tdci.dat
  integer(8) :: &
       nbasis, & !: total number of basis functions
       nrorb,  & !: total number of occupied+virtual  MOs for alpha = for beta
       noa,    & !: total number of occupied alpha MOs
       nva,    & !: total number of virtual alpha MOs
       nob,    & !: total number of occupied beta MOs
       nvb       !: total number of virtual beta MOs
  

  !: computed info for electronic system
  integer(8) :: &
       norb,     & !: if unrestricted, norb=2*nrorb ; else norb=nrorb
       noanva,   & !: noa*nva
       nobnvb,   & !: nob*nvb
       noanvb,   & !: noa*nvb
       nobnva,   & !: nob*nva
       nstates,  & !: total number of CIS states
       nstuse,   & !: number of CIS states that will be used; set from ionization+eigmax threshold
       ip_states,& !: total number of ionized states
       read_states !: number of states to read for SDI
  
  
  
  !: computed info for e-field
  real(8)    ::        &
       period       ,  & !: total number of timesteps per period of E-field
       field_duration    !: total number of timesteps during the field duration
  
  
  real(8), allocatable :: &
       fvect1(:) ,        & !: contains shaped e-field ; for lin and circ pulse
       fvect2(:)            !: contains shaped e-field ; only used for circ pulse


  integer(8), allocatable :: hole_index(:,:), part_index(:,:)        !: holes and particle indices for states
  integer(8), allocatable :: hole_ip_index(:,:)                      !: holes and particle indices for ionized states
  integer(8), allocatable :: part_ip_index(:,:)                      !: holes and particle indices for ionized states
  integer(8), allocatable :: state_index(:,:), cisd_indices(:,:,:,:) !: cisd indices
  integer(8), allocatable :: state_ip_index(:,:)                     !: ionized state indices

  
  real(8), allocatable :: &
       cis_vec(:),        & !: CIS eigenstates 
       cisd_vec(:,:),     & !: cisd 
       cis_eig(:),        & !: CIS eigenvalue energies
       ip_vec(:),         & !: ionized state eigenvectors
       ip_eig(:),         & !: ionized state energies
       proj_ion(:)          !: project ion coefficients from wavefunction
  
  complex(8), allocatable :: &
       Zcis_vec(:),&  !: COMPLEX CIS eigenstates for soc
       Zip_vec(:), &  !: ionized state eigenvectors
       ion_coeff(:),& !: coefficients of ionized states in absorbed wavefunction
       Zproj_ion(:)   !: project ion coefficients from wavefunction

  real(8), allocatable ::       &
       abp(:), exp_abp(:),      & !: Vabs matrix elements in CIS basis and [R][exp(abp')][R]^T
       tdx(:), tdy(:), tdz(:),  & !: transition dipole moments in CIS basis  
       polar(:,:)                 !: dipole polarizability

  complex(8), allocatable ::    &
       Zabp(:), Zexp_abp(:),    & !: COMPLEX Vabs matrix elements in CIS basis and [R][exp(abp')][R]^T
       Ztdx(:), Ztdy(:), Ztdz(:)  !: COMPLEX transition dipole moments in CIS basis

       
  type(tdcidat), allocatable :: tdciresults(:)

  
  !: arrays to write files
  integer(8)    :: ndata !: number of times to record dat



  !: misc_variables 
  integer(8) :: torbitals = 0 !: total # of orbitals including core orbitals (read in from checkpoint file)
  
  
  !: timestamp variables
  real(8) :: start, finish
  

contains
  !:-----------------!
  ! subroutine dnt 
  !:-----------------!
  subroutine dnt(iout)


    !: Retrieve Date and Time
    !: value(1)=year...value(2)=month...value(3)=day...value(4)=timezone
    !: value(5)=hour...value(6)=mins....value(7)=sec...value(8)=millisec


    implicit none

    integer(8), intent(in) :: iout

    integer(8) :: values(8)

    character(3) :: month, day 
    character(5) :: zone
    character(8) :: date
    character(10):: time
    character(2) :: ampm
    

    values=0
    call DATE_AND_TIME(date,time,zone,values)
    select case ( values(2) )
    case (1)  ; month = 'Jan'
    case (2)  ; month = 'Feb'
    case (3)  ; month = 'Mar'
    case (4)  ; month = 'Apr'
    case (5)  ; month = 'May'
    case (6)  ; month = 'Jun'
    case (7)  ; month = 'Jul'
    case (8)  ; month = 'Aug'
    case (9)  ; month = 'Sep'
    case (10) ; month = 'Oct'
    case (11) ; month = 'Nov'
    case (12) ; month = 'Dec'
    end select

    if ( values(5).ge.12 ) then
       ampm = "PM"
       values(5) = mod( values(5),12 )
    else if ( values(5).eq.0 ) then
       ampm = "AM"
       values(5) = 12
    else
       ampm = "AM"
    end if

    !: WRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITE
    write(iout,300) month, values(3), values(1), values(5), values(6), values(7), values(8), ampm
    
300 format(1x,a3,1x,i2,1x,i4,3x,i2,':',i2,':',i2,':',i3,1x,a2)

  end subroutine dnt
  !:--------------------!
  ! subroutine stopme
  !:--------------------!
  subroutine stopme

    implicit none


    write(iout,'(A)') divide
    write(iout,'(A)') ' '
    write(iout,'(A)')  ' STOP IN THE NAME OF LOVE'
    call dnt(iout) 
    flush(iout)
  

    stop

    
  end subroutine stopme
  !:---------------------------!
  ! subroutine write_header
  !:---------------------------!
  subroutine write_header( curr_sub, curr_mod, enter_or_leave )
    
    implicit none

    character(*), intent(in) :: curr_sub, curr_mod, enter_or_leave
    
    select case( trim(enter_or_leave) )
    case( 'enter' )
       write( iout,'(A)' ) ' in subroutine '//trim(curr_sub)//' in MODULE '//trim(curr_mod)
       write( iout,'(A)' ) ' '
    case( 'leave' )
       write( iout,'(A)' ) ' '
       write( iout,'(A)' ) ' leaving subroutine '//trim(curr_sub)
       write( iout,'(A)' ) divide
    end select
    flush(iout)

  end subroutine write_header
  !:----------------------------!
  ! subroutine track_mem
  !:----------------------------!
  subroutine track_mem( nelements )
  
    integer(8), intent(in) :: nelements
    

    tot_use_mem = tot_use_mem + nelements 
    write(iout,"(' total GB memory in use',f10.5)") dble(tot_use_mem * 8) / giga2byte
    
    
  end subroutine track_mem
  !:----------------------------!
  !:----------------------------!
end module variables_global
