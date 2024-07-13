module initialize
  

  use variables_global
  implicit none


  !: contains subroutine read_input
  !: contains subroutine set_default
  !: contains subroutine read_tdcihead
  !: contains subroutine set_variables
  !: contains subroutine allocate_main( option )
  !: contains subroutine deallocate_main( option )
  !: contains subroutine read_hamdata
  !: contains subroutine write_input( option )
  !: contains subroutine writeme_modreadin( myroutine, option )
  !: contains subroutine erorrs_modreadin( myroutine, option )
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_INPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_input  
    

    !: reads in file 'input'    
    !: sets default parameters if file 'input' is not found
    implicit none

    
    integer(8), parameter :: myfile=100
    
    logical :: iamhere
    integer(8) :: i, mystat
    

    !: /DYNAMICS/         init_coeffs(0:10), init_states(0:10).  to change, mod_variables_global.f90
    !: /FIELD_strengths/  read_emax(20).  to change, mod_variables_global.f90
    !: /FIELD_directions/ read_*(100).  to change, mod_variables_global.f90
    namelist /DYNAMICS/         init_coeffs, init_states, restart, Qsave
    namelist /FIELD/            dirform, ellipt, envelope, ncyc, omega, phase, euler
    namelist /FIELD_units/      omega_units !:, angle_units
    namelist /FIELD_strengths/  nemax, read_emax, read_state1, read_state2, read_coeff1, read_coeff2, &
                                read_shift, ion_sample_start, ion_sample_width, ion_sample_state
    namelist /FIELD_directions/ ndir, read_theta, read_phi, read_x, read_y, read_z 
    namelist /SYSTEM/           dt, eigmax, heuristic, ionization, jobtype, nstep, outstep, &
                                nactive, nvirtual, IP_alpha_beta, socfac, socfacz, ffieldx, &
                                ffieldy, ffieldz, QsocA2B, QeigenDC, Qserial
    namelist /SYSTEM_units/     dt_units, eigmax_units
    namelist /InOutFILES/       tdcidatfile, outputfile, restartbinfile, &
                                Qread_TDCIdata, Qwrite_ion_coeff, Qread_ion_coeff, &
                                Qmo_dens, Qci_save, write_binaries
    namelist /DAVIDSON/         flag_davidson
    namelist /ReadU_NO/         flag_ReadU_NO
    

    !: set default optional parameters
    call set_default

    !: in input file does not exist, set default parameters
    inquire( file=trim(inputfile), exist=iamhere )        
    if ( .not.iamhere ) call write_input(myfile)
    

    !: read input.  see namelists.  not all WARNINGs are fatal
    open( unit=10, file=trim(inputfile) )

    read( 10, nml=FIELD, iostat=mystat,err=40 ) 
40  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD, iostat= ',i0)") mystat 
    rewind(10)

    read( 10, nml=FIELD_units, iostat=mystat,err=41 )
41  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_units, iostat= ',i0)") mystat 
    rewind(10)

    read( 10, nml=FIELD_strengths, iostat=mystat, err=47 )
47  if ( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_strengths, iostat= ',i0)") mystat
    rewind(10)

    read( 10, nml=FIELD_directions, iostat=mystat, err=48 )
48  if ( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_directions, iostat= ',i0)") mystat
    rewind(10)    
    
    read( 10, nml=SYSTEM, iostat=mystat, err=42 ) 
42  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=SYSTEM, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10, nml=SYSTEM_units, iostat=mystat, err=43 )
43  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=SYSTEM_units, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10,nml=DYNAMICS,iostat=mystat,err=44 )
44  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=DYNAMICS, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10, nml=InOutFILES, iostat=mystat, err=45 )
45  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=InOutFILES, iostat= ',i0)") mystat
    rewind(10)

    read( 10, nml=ReadU_NO, iostat=mystat, err=49 )
49  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=ReadU_NO, iostat= ',i0)") mystat
    rewind(10)

    read( 10, nml=DAVIDSON, iostat=mystat, err=46 )
46  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=DAVIDSON, iostat= ',i0)") mystat 
    continue
    
    close(10)


    !: officially opening outputfile.  
    open( iout,file=trim(outputfile) )


    call writeme_modreadin( 'read_input', 'greeting' ) 
    call writeme_modreadin( 'read_input', 'date' )
    
    if ( .not.iamhere ) write(iout,'(A)') " WARNING: Could not find file 'input'. Using default input settings"
    write(iout,'(A)') " finished reading Field and System variables from file 'input' "
    call write_input(iout)
    
    
    call write_header( 'read_input', 'initialize', 'leave' )
    

  end subroutine read_input
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SET_DEAFULTINPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine set_default
    

    !: sets default values for variables read in from file 'input'
    implicit none
    

    !: default for /DYNAMICS/
    !: default init staate = ground state
    !: # of superposition states = init_stats(0) = 1
    init_coeffs    = dcmplx( 0.d0,0.d0 )
    init_coeffs(1) = dcmplx( 1.d0,0.d0 )
    init_states(0) = 1  
    init_states(1) = 1
    restart = .False.

    !: default for /FIELD/
    dirform  = 'polar'
    ellipt   = 1.d0
    envelope = 'trap'
    ncyc     = 7
    omega    = 0.057d0
    phase    = 90.d0
    euler    = 0.d0
    
    !: default for /FIELD_units/
    omega_units = 'au'
    !angle_units = 'deg'
    
    !: default for /FIELD_strengths/
    nemax  = 3 
    read_emax = -1000.d0  !: will assign later with MO energies if not read-in
    read_state1 = 0
    read_state2 = 0
    read_coeff1 = dcmplx(0.d0, 0.d0)
    read_coeff2 = dcmplx(0.d0, 0.d0)
    read_shift  = 0
    ion_sample_start = 0
    ion_sample_width = 16000
    ion_sample_state = 0
    
    !: default for /FIELD_directions/   
    ndir     = 62
    read_theta = -1000.d0
    read_phi   = -1000.d0
    
    !: default for /InOutFiles/
    tdcidatfile     = 'TDCI.dat'
    outputfile      = 'OUTPUT'
    restartbinfile  = 'none'

    !: default for /SYSTEM/
    dt         = 0.05d0
    eigmax     = 10.d0
    ionization = -1.0 
    jobtype    = flag_cis
    nstep      = 16000
    outstep    = 50
    nactive    = -1
    nvirtual   = -1
    socfac     = 1.D0
    socfacz    = 1.D0
    ffieldx    = 0.d0
    ffieldy    = 0.d0
    ffieldz    = 0.d0

    !: default for /SYSTEM_units/
    dt_units     = 'au'
    eigmax_units = 'au'    
    
    !: default for davidson diagonalization
    flag_davidson = .False.

    flag_ReadU_NO = .False.
    write_binaries = .false.

    
  end subroutine set_default
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_TDCIHEAD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_tdcihead
    

    !: read in 'TDCI.dat'.  Exits if file does not exist
    implicit none
    
    logical        :: iamhere
    character(100) :: cskip
    integer(4)     :: i
  
  
    call write_header( 'read_tdcihead','initialize','enter' )
    
    
    !: stop if TDCI.dat file does not exist
    inquire( file=trim(tdcidatfile), exist=iamhere )
    if ( .not.iamhere ) call errors_modreadin( 'read_tdci1', 'no_tdcidat' )    
    
    !: read TDCI.dat
    open( unit=10, file=trim(tdcidatfile) )    
    read(10,*)     nbasis, nrorb, noa, nva, nob, nvb, Mol%unrstedflag, Mol%vabsflag
    read(10,'(A)') Mol%job_title
    read(10,*)     Mol%ICharg, MoL%Multip, Mol%natoms    
   
    if(nactive.gt.noa) nactive = noa
    if(nvirtual.gt.nva) nvirtual = nva 

    !: allocate coordinate arrays
    allocate( Mol%xcoord(Mol%natoms), &
         Mol%ycoord(Mol%natoms), &
         Mol%zcoord(Mol%natoms), &
         Mol%myatom(Mol%natoms) )    
    read(10,*) ( (Mol%myatom(i), &
         Mol%xcoord(i), &
         Mol%ycoord(i), &
         Mol%zcoord(i)) , i=1, Mol%natoms )
    read(10,*) cskip, cskip, Mol%dipx00, Mol%dipy00, Mol%dipz00
    
    if( restart ) then
       close(10) 
       write(iout,'(A)') " finished reading heading in '"//trim(tdcidatfile)//"'"
    else
       write(iout,'(A)') " finished reading heading in '"//trim(tdcidatfile)//"'"
       write(iout,'(A)') ' '//trim(tdcidatfile)//" still open for reading "
    end if


    call write_header( 'read_tdcihead','initialize','leave' )
    
    
  end subroutine read_tdcihead
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SET_VARIABLES
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine set_variables


    !: set global variables.  convert units to atomic units
    implicit none
    
    integer(8) :: itheta, iphi, idir, nphi, ntheta, natrim, nbtrim, nva1, nvb1
    real(8) :: norm, dumtheta, dumphi, dtheta, dphi
    logical :: iamhere, gen_direction
    
    
    call write_header( 'set_variables','initialize','enter' )
    
    
    !: check units, convert to AU
    if ( trim(omega_units).ne.'au' ) then
       if ( trim(omega_units).ne.'nm' ) call errors_modreadin( 'read_input', 'omega_error' )
       omega = convert_freq(omega)
    end if
    
    if ( trim(dt_units).ne.'au' ) then
       dt = convert_time( dt, trim(dt_units), 'au' )
       if ( dt.eq.-1.d20 ) call errors_modreadin( 'read_input', 'dt_error' )
    end if
    
    !if ( trim(angle_units).ne.'rad' ) then
    !   if ( trim(angle_units).ne.'deg' ) call errors_modreadin( 'read_input', 'euler_error' )
    !   euler = convert_angles( euler )
    !   phase = convert_angles( phase )
    !end if
    
    if ( trim(eigmax_units).ne.'au' ) then
       eigmax = convert_energy( eigmax, trim(eigmax_units), 'au' )
       if ( eigmax.eq.0.d0 ) call errors_modreadin( 'read_input', 'eigmax_error' )
    end if
    
    
    !: set up flag for linear or circular pulse
    linear = .True.
    if ( envelope.eq.'cirl' .or. envelope.eq.'cirr' ) linear = .False.
    
    !: if unreasonable envelope, set to default
    select case ( trim(envelope) ) 
    case( 'none' ) ; go to 150
    case( 'cos2' ) ; go to 150
    case( 'trap' ) ; go to 150
    case( 'gaus' ) ; go to 150
    case( 'stat' ) ; go to 150
    case( 'band' ) ; go to 150
    case( 'ramp' ) ; go to 150
    case( 'sin2' ) ; go to 150
    case( 'cirl' ) ; go to 150
    case( 'cirr' ) ; go to 150
    case default 
       envelope = 'trap'
       linear   = .true.
       call errors_modreadin( 'set_variables', 'no_envelope' )
    end select
    
150 continue

    
    !: if unreasonable eigmax, set eigmax to default
    if(eigmax.eq.0.d0) then
       eigmax = 10.d0
       call errors_modreadin( 'read_input', 'eigmax_error' )
    end if
    

    !: if unreasonable nstep, set nstep to 1.5 * number of cycles + 2
    period = int( 2.d0 * pi / omega / dt )    !: total number of timestep per period
    if ( nstep.le.1 ) then       
       nstep  = int( 1.5d0 * dble(ncyc) * dble(period) + 2.d0 ) !: total number of propagation steps
       call errors_modreadin( 'set_variables', 'nstep_error' )
    end if
    nstep = nstep + 1 
    field_duration = dble(ncyc) * period + 1
    

    if ( ionization.lt.0.d0 ) write(iout,"(' ionization < 0, setting energy cutoff automatically')")
    if ( ionization.gt.0.d0 ) write(iout,"(' ionization > 0, energy cutoff set to ', f10.4, 'au')") ionization
    
    write(iout,'(A)') ''
    if ( trim(dirform).eq.'polar' ) write(iout,"(' field/light polarization/prop directions in polar coordinates')")
    if ( trim(dirform).eq.'cart'  ) write(iout,"(' field/light polarization/prop directions in Cartesian coordinates')") 

    
    !: set default directions if not specified     
    gen_direction = .False.
    do idir=1, ndir 
       if ( read_theta(idir).lt.-100.0 .or. read_phi(idir).lt.-100.0 ) gen_direction = .True.
    end do

    if ( gen_direction ) then
       dirform = 'polar'
       read_theta = 0.d0
       read_phi   = 0.d0
       dtheta = 30.d0     ;    ntheta = int(180.d0/dtheta)
       dphi   = 30.d0     ;    nphi   = int(360.d0/dphi) 
       write(iout,*) ntheta, dphi
       write(iout,'(A)') ' WARNING:  setting ndir = 62 for automatic field generation' ; write(iout,'(A)') ''
       write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") 1, read_theta(1), 1, read_phi(1)
       idir = 1
       do itheta=2, ntheta
          dumtheta = dble(itheta-1)*dtheta
          do iphi=1, nphi
             idir = idir + 1
             read_theta(idir) = dumtheta
             read_phi(idir)   = dble(iphi-1)*dphi
             write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") idir, read_theta(idir), idir, read_phi(idir)
          end do
       end do
       idir = idir + 1 
       read_theta(idir) = 180.d0
       read_phi(idir)   = 0.d0
       write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") idir, read_theta(idir), idir, read_phi(idir)
    end if
    
    
    if( read_emax(1).lt.-100.0 ) nemax = 3 

    !: outsteps
    ndata = int(nstep/outstep) 
    write(iout,'(A)') '' ; write(iout,"(' write norm and wavefunction out every ',i0,' steps, total of ' i0,' steps')") outstep, ndata
    
    
    !: DATA FROM TDCI.DAT    

    if( trim(jobtype) .eq. flag_ip ) Qalloc_indices = .True.
    if( Mol%vabsflag .eq. 0 ) Qalloc_vabs = .False.
    if( trim(jobtype) .eq. flag_soc .or.  trim(jobtype) .eq. flag_socip ) then
       Qalloc_Zcomplex = .True.
       Qalloc_socmo    = .True.
       Qread_socx_ao   = .True.
       Qread_socy_ao   = .True.
       Qread_socz_ao   = .True.
       Qwrite_ham0     = .False.
    end if
    

    !: restart, set job route
    if( restart ) then
       if ( trim(restartbinfile) .eq. 'none' ) call errors_modreadin( 'set_variables','norestart' )
       Qallocate_main = .True.
       Qread_hamdata  = .False.
       Qgetfield     = .True.
       Qgetham0      = .False.
       Qget_nstuse    = .False.
       Qalloc_vabsmo  = .True.
       Qalloc_dipmo   = .True.       
       Qget_1eham     = .False.
       Qget_expVabs   = .False.
       Qread_binaries = .True.
       Qsave          = .False.
       Qdealloc       = .False.
       Qpropagate     = .True.
       Qserial        = .False.
    end if

    !: davidson flag
    if ( flag_davidson ) then
       Qallocate_main = .False.
       Qread_hamdata  = .True.
       Qgetfield     = .False.
       Qwrite_specifics1 = .False.
       Qwrite_ham0    = .False.
       Qwrite_specifics2 = .False.
       Qgetham0      = .False.
       Qget_nstuse    = .False.
       Qget_1eham     = .False.
       Qget_expVabs   = .False.
       Qread_binaries = .False.
       Qsave          = .False.
       Qdealloc       = .False.
       Qpropagate     = .False.
       Qserial        = .False.
       !: alloc variables
       Qalloc_vabsmo  = .False.
       Qalloc_dipmo   = .False.
       Qalloc_socmo   = .False.
       Qread_vabs_ao  = .False.
       Qread_dipx_ao  = .True.
       Qread_dipy_ao  = .True.
       Qread_dipz_ao  = .True.
       Qread_socx_ao  = .False.
       Qread_socy_ao  = .False.
       Qread_socz_ao  = .False.
       Qread_orben    = .True.
       Qread_cmo      = .True.
    end if
    call writeme_modreadin( 'set_variables','routine' )

    
    !: parameters for unrestricted or restricted
    nva1 = nva
    if(nvirtual.gt.0) nva1 = nvirtual
    nvb1 = nvb
    if(nvirtual.gt.0) nvb1 = nvirtual + noa - nob
    select case( trim(jobtype) ) 
    case( flag_cis )
       noanva = noa*nva
       nobnvb = nob*nvb
       if(Mol%unrstedflag.ne.0) then
          unrestricted  = .True.
          norb          = noa + nva + nob + nvb
          nstates       = noa*nva1 + nob*nvb1 + 1
          ip_states = noa + nob
       else
          unrestricted = .False.
          norb         = nrorb
          nstates      = noa*nva1 + 1
          ip_states = noa
       end if
    case( flag_soc )
       unrestricted = .True.
       noanva = noa*nva
       nobnvb = nob*nvb
       noanvb = noa*nvb
       nobnva = nob*nva
       nstates = noa*nva1 + nob*nvb1 + noa*nvb1 + nob*nva1 + 1 
       If(.not.QsocA2B) nstates = nstates - noa*nvb1 -nob*nva1
       norb = noa + nva + nob + nvb
       ip_states = noa + nob
    case( flag_tda )
       noanva = noa*nva
       nobnvb = nob*nvb
       if(Mol%unrstedflag.ne.0) then
          unrestricted  = .True.
          norb          = noa + nva + nob + nvb
          nstates       = noa*nva1 + nob*nvb1 + 1
          ip_states = noa + nob
       else
          unrestricted = .False.
          nstates      = noa*nva1 + 1
          norb         = nrorb
          ip_states = noa
       end if
    case ( flag_ip ) 
       unrestricted = .True.
       if ( .not. IP_alpha_beta ) then 
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = nob*(nob-1)/2*nvb1 + nob*noa*nva1 + nob
          ip_states    = nob*(nob-1)/2 + nob*noa
          read_states  = nob
          if(nactive.gt.0) then 
             nstates=nactive*(noa*nva1+nob*nvb1)+nactive- &
                     ((nactive+1)*nactive/2)*nvb1
             read_states = nactive
          end If
       else
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = nob*(nob-1)/2*nvb1 + nob*noa*nva1 + nob +  &
                         noa*(noa-1)/2*nva1 + noa*nob*nvb1 + noa           
          ip_states    = noa*(noa-1)/2 + nob*(nob-1)/2 + noa*nob           
          read_states  = noa + nob
          if(nactive.gt.0) then
             nstates=2*nactive*(noa*nva1+nob*nvb1)+2*nactive- &
                     ((nactive+1)*nactive/2)*(nva1+nvb1)
             read_states = 2*nactive
          end If
       end if
    case ( flag_socip ) 
       unrestricted = .True.
       if ( .not. IP_alpha_beta ) then 
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = (nob*(nob-1)/2)*(nva1+nvb1) + nob*noa*nva1 + nob
          if(.not.QsocA2B) nstates = nstates - (nob*(nob-1)/2)*nva1
          ip_states    = noa*(noa-1)/2 + noa*nob           
          read_states  = nob
          if(nactive.gt.0) then
             nbtrim = nob - nactive
             nstates=nstates-nbtrim-(nbtrim*(nbtrim-1)/2)*(nva1+nvb1)-natrim*nbtrim*(nva1+nvb1)
!:             nstates=nstates-nbtrim-(nbtrim*(nbtrim-1)/2)*(nva1+nvb1)-nbtrim*noa*nva1
             if(.not.QsocA2B) nstates = nstates + (nbtrim*(nbtrim-1)/2)*nva1
             ip_states=ip_states-nbtrim*(nbtrim-1)/2-nbtrim*noa
             read_states = nactive
          end if  
       else
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = noa*(noa-1)/2*(nva1+nvb1) + noa + nob +  &
                         nob*(nob-1)/2*(nva1+nvb1) + noa*nob*(nva1+nvb1)           
          if(.not.QsocA2B) nstates = nstates - (noa*(noa-1)/2)*nvb1 - (nob*(nob-1)/2)*nva1
          ip_states    = noa*(noa-1)/2 + nob*(nob-1)/2 + noa*nob           
          read_states  = noa + nob
          if(nactive.gt.0) then
             natrim = noa - nactive
             nbtrim = nob - nactive
             nstates=nstates-natrim*(natrim-1)/2*(nva1+nvb1)-natrim-nbtrim- &
               nbtrim*(nbtrim-1)/2*(nva1+nvb1)-natrim*nbtrim*(nva1+nvb1)
             if(.not.QsocA2B) nstates = nstates + (natrim*(natrim-1)/2)*nvb1 &
                                                + (nbtrim*(nbtrim-1)/2)*nva1
             ip_states=ip_states-natrim*(natrim-1)/2+nbtrim*(nbtrim-1)/2-natrim*nbtrim
             read_states = 2*nactive
          end if
       end if

    case( flag_cisd ) 
       if(Mol%unrstedflag.ne.0) then
          unrestricted = .True.
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates = noanva*nobnvb + (noa)*(noa-1)*(nva)*(nva-1)/4 + (nob)*(nob-1)*(nvb)*(nvb-1)/4 + &
               noanva + nobnvb + 1 !: for singles and ground state
       else
          unrestricted = .False.
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = (noa)*(noa-1)*(nva)*(nva-1)/4 + noanva + 1
       end if
    end select
    
    
    !: initialize to CIS parameters
    nstuse  = nstates


    !: keep track of mem
    tot_use_mem = 0

    
    call write_header( 'set_variables','initialize','leave' )
    
    
  end subroutine set_variables
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ALLOCATE_MAIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine allocate_main( option )


    
    !: allocate main arrays
    implicit none
    
    character(*), intent(in) :: option
    integer(8), parameter :: npart = 2 ! nparticle
    

100 format(' allocated rank 1 array',a15,' of length = ',i11)


    call write_header( 'allocate_main','initialize','enter' )
    
    
    select case( trim(option) )
    case( 'main' )
       
       !: arrays for field 
       if( Qalloc_field ) then
          allocate( fvect1(nstep) )        ; fvect1 = 0.d0
          allocate( Mol%field_env(nstep) )           ; Mol%field_env = 0.d0
          write(iout,100) 'fvect1', nstep  ; call track_mem( nstep )
          write(iout,100) 'Mol%field_env',    nstep  ; call track_mem( nstep )
          !: circular stuff
          if ( .not.linear ) then
             allocate( fvect2(nstep) )        ; fvect2 = 0.d0
             write(iout,100) 'fvect2', nstep  ; call track_mem(nstep)
          end if
          allocate( tdciresults(nemax*ndir) ) 
          write(iout,"(' allocated tdciresults ')")
       end if

       !: for eigen-vectors and eigen-values
       select case( trim(jobtype) )
       case( flag_cis ) 
          allocate( cis_vec(nstates*nstates) )     ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )             ; cis_eig = 0.d0
          allocate( ip_vec(ip_states*ip_states) )  ; ip_vec = 0.d0
       case( flag_soc ) 
          allocate( Zcis_vec(nstates*nstates) )    ; Zcis_vec = dcmplx(0.d0,0.d0)
          allocate( cis_eig(nstates) )             ; cis_eig = 0.d0 
          allocate( Zip_vec(ip_states*ip_states) ) ; Zip_vec = dcmplx(0.d0,0.d0)
       case( flag_ip ) 
          allocate( cis_vec(nstates*nstates) )     ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )             ; cis_eig = 0.d0
          allocate( ip_vec(ip_states*ip_states) )  ; ip_vec = 0.d0
       case( flag_socip ) 
          allocate( Zcis_vec(nstates*nstates) )    ; Zcis_vec = dcmplx(0.d0,0.d0)
          allocate( cis_eig(nstates) )             ; cis_eig = 0.d0 
          allocate( Zip_vec(ip_states*ip_states) ) ; Zip_vec = dcmplx(0.d0,0.d0)
       case ( flag_cisd )
          allocate( cis_vec(nstates*nstates) )     ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )             ; cis_eig = 0.d0
       end select
       allocate( ip_eig(ip_states) )               ; ip_eig = 0.d0 
       allocate( ion_coeff(ip_states) )            ; ip_eig = 0.d0 

       if ( Qalloc_Zcomplex ) then
          write(iout,100) 'cis_vec', 2*nstates*nstates      ;  call track_mem( 2*nstates*nstates)
          write(iout,100) 'cis_eig', 2*nstates              ;  call track_mem( 2*nstates )
          write(iout,100) 'Zip_vec', 2*ip_states*ip_states  ;  call track_mem( 2*ip_states*ip_states)
          write(iout,100) 'ip_eig',  2*ip_states            ;  call track_mem( 2*ip_states )
       else
          write(iout,100) 'cis_vec', nstates*nstates        ;  call track_mem( nstates*nstates)
          write(iout,100) 'cis_eig', nstates                ;  call track_mem( nstates )
          write(iout,100) 'Zip_vec', 2*ip_states*ip_states  ;  call track_mem( 2*ip_states*ip_states)
          write(iout,100) 'ip_eig',  2*ip_states            ;  call track_mem( 2*ip_states )
       end if

       !: allocate arrays to keep track of indices .  optional
       if ( Qalloc_indices ) then
          select case( trim(jobtype) ) 
          case( flag_cis ) 
             allocate( hole_index(nstates,2) )
             allocate( part_index(nstates,1) )
             write(iout,100) 'hole_index', 2*nstates
             write(iout,100) 'part_index', nstates
             call track_mem( 2*nstates )
             hole_index = 0; part_index = 0
          case( flag_soc ) 
             allocate( hole_index(nstates,2) )
             allocate( part_index(nstates,1) ) 
             write(iout,100) 'hole_index', 2*nstates 
             write(iout,100) 'part_index', nstates 
             call track_mem( 2*nstates )
             hole_index = 0; part_index = 0
          case( flag_tda ) 
             allocate( hole_index(nstates,2) ) 
             allocate( part_index(nstates,1) )
             write(iout,100) 'hole_index', 2*nstates 
             write(iout,100) 'part_index', nstates             
             call track_mem( 2*nstates )
             hole_index = 0; part_index = 0
          case( flag_ip )
             allocate( hole_index(nstates,2) )
             allocate( part_index(nstates,1) )
             write(iout,100) 'hole_index', 2*nstates
             write(iout,100) 'part_index', nstates 
             call track_mem( 3*nstates )
             hole_index = 0; part_index = 0
          case( flag_socip )
             allocate( hole_index(nstates,2) ) 
             allocate( part_index(nstates,1) )
             write(iout,100) 'hole_index', 2*nstates 
             write(iout,100) 'part_index', nstates 
             call track_mem( 3*nstates )
             hole_index = 0; part_index = 0
          end select
          allocate( hole_ip_index(ip_states,2) ) 
          allocate( part_ip_index(ip_states,1) )
          allocate( state_ip_index(noa+nob,noa+nob) )
          write(iout,100) 'hole_ip_index', 2*ip_states
          write(iout,100) 'part_ip_index', ip_states
          write(iout,100) 'state_ip_index', (noa+nob)*(noa+nob)
          call track_mem( 3*ip_states+(noa+nob)*(noa+nob) )
          hole_ip_index = 0; part_ip_index = 0; state_ip_index = 0
       end if
              
       !: for Vabs elements 
       if( Qalloc_vabs ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Zabp(nstates*nstates) )         ; Zabp = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Zabp', 2*nstates*nstates ; call track_mem( 2*nstates*nstates )
          else 
             allocate( abp(nstates*nstates) )         ;  abp = 0.d0
             write(iout,100) 'abp',  nstates*nstates  ;  call track_mem( nstates*nstates )
          end if
       end if
       
       !: for [R][e^(-Vabs)][R]^T elements
       if( Qalloc_exp_vabs ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Zexp_abp(nstates*nstates) )          ; Zexp_abp = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Zexp_abp', 2*nstates*nstates  ; call track_mem( 2*nstates*nstates )
          else
             allocate( exp_abp(nstates*nstates) )        ; exp_abp = 0.d0
             write(iout,100) 'exp_abp', nstates*nstates  ; call track_mem( nstates*nstates )
          end if
       end if
       
       !: for x, y, z dipoles in CIS basis
       if( Qalloc_tdxyz ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Ztdx(nstates*nstates) )  ;  tdx = dcmplx( 0.d0,0.d0 )
             allocate( Ztdy(nstates*nstates) )  ;  tdy = dcmplx( 0.d0,0.d0 )
             allocate( Ztdz(nstates*nstates) )  ;  tdz = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Ztdx', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
             write(iout,100) 'Ztdy', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
             write(iout,100) 'Ztdz', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
          else             
             allocate( tdx(nstates*nstates) )        ;  tdx = 0.d0
             allocate( tdy(nstates*nstates) )        ;  tdy = 0.d0
             allocate( tdz(nstates*nstates) )        ;  tdz = 0.d0
             write(iout,100) 'tdx', nstates*nstates  ;  call track_mem( nstates*nstates )
             write(iout,100) 'tdy', nstates*nstates  ;  call track_mem( nstates*nstates )
             write(iout,100) 'tdz', nstates*nstates  ;  call track_mem( nstates*nstates )
          end if
          allocate( polar(nstates,6) )               ;  polar = 0.d0
          write(iout,100) 'polar',nstates*6          ;  call track_mem( nstates*6 )
       end if
       
       
    case( 'mo' )
       
       !: Vabs MO elements
       if( Qalloc_vabsmo ) then
          allocate( Mol%vabsmoa(nrorb,nrorb) )            ;  Mol%vabsmoa = 0.d0
          write(iout,100) 'vabs_moa', nrorb*nrorb     ;  call track_mem( nrorb*nrorb )
          if( unrestricted ) then
             allocate( Mol%vabsmob(nrorb,nrorb) )         ;  Mol%vabsmob = 0.d0
             write(iout,100) 'vabs_mob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          end if
       end if

       !: dipx, dipy, dipz MO elements
       if( Qalloc_dipmo ) then
          allocate( Mol%dipxmoa(nrorb,nrorb) )        ;  Mol%dipxmoa = 0.d0
          allocate( Mol%dipymoa(nrorb,nrorb) )        ;  Mol%dipymoa = 0.d0
          allocate( Mol%dipzmoa(nrorb,nrorb) )        ;  Mol%dipzmoa = 0.d0
          write(iout,100) 'Mol%dipxmoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          write(iout,100) 'Mol%dipymoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          write(iout,100) 'Mol%dipzmoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          if( unrestricted ) then
             allocate( Mol%dipxmob(nrorb,nrorb) )        ;  Mol%dipxmob = 0.d0
             allocate( Mol%dipymob(nrorb,nrorb) )        ;  Mol%dipymob = 0.d0
             allocate( Mol%dipzmob(nrorb,nrorb) )        ;  Mol%dipzmob = 0.d0
             write(iout,100) 'Mol%dipxmob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
             write(iout,100) 'Mol%dipymob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
             write(iout,100) 'Mol%dipzmob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          end if
       end if

       if ( Qalloc_socmo ) then
          allocate( Mol%socmoAA(nrorb,nrorb) )  ;  Mol%socmoAA = 0.d0
          allocate( Mol%socmoBB(nrorb,nrorb) )  ;  Mol%socmoBB = 0.d0
          allocate( Mol%socmoAB(nrorb,nrorb) )  ;  Mol%socmoAB = 0.d0
          write(iout,100) 'Mol%socmoAA', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
          write(iout,100) 'Mol%socmoBB', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
          write(iout,100) 'Mol%socmoAB', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
       end if

    end select
    
    
    call write_header( 'allocate_main','initialize','leave' )


  end subroutine allocate_main
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DEALLOCATE_MAIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine deallocate_main( option )
    
    
    use readintegrals
    
    implicit none

    character(*), intent(in) :: option
    

    call write_header( 'deallocate_main','initialize','enter' )

    select case( trim(option) )
    case( '2e_int' )
       if( allocated(Mol%dabcdBB) ) then
          deallocate( Mol%dabcdBB )
          write(iout,"( deallocated Mol%dabcdBB )" )  ;  call track_mem( -nvb3*nvb3 )
       end if
       if( allocated(Mol%dabcdAB) ) then
          deallocate( Mol%dabcdAB )
          write(iout,"(' deallocated Mol%dabcdAB ')")  ;  call track_mem( -nva2*nvb2 )
       end if
       if( allocated(Mol%dabcdAA) ) then
          deallocate( Mol%dabcdAA )  
          write(iout,"(' deallocated Mol%dabcdAA ')")  ;  call track_mem( -nva3*nva3 )
       end if
       if( allocated(Mol%diabcBB) ) then
          deallocate( Mol%diabcBB )
          write(iout,"(' deallocated Mol%diabcBB ')")  ;  call track_mem( -nobnvb*nvb3 )
       end if
       if( allocated(Mol%diabcBA) ) then
          deallocate( Mol%diabcBA )
          write(iout,"(' deallocated Mol%diabcBA ')")  ;  call track_mem( -nobnvb*nva2 )
       end if
       if( allocated(Mol%diabcAB) ) then
          deallocate( Mol%diabcAB )
          write(iout,"(' deallocated Mol%diabcAB ')")  ;  call track_mem( -noanva*nvb2 )
       end if
       if( allocated(Mol%diabcAA) ) then
          deallocate( Mol%diabcAA )
          write(iout,"(' deallocated Mol%diabcAA ')")  ;  call track_mem( -noanva*nva3 )
       end if
       if( allocated(Mol%dijkaBB) ) then
          deallocate( Mol%dijkaBB )
          write(iout,"(' deallocated Mol%dijkaBB ')")  ;  call track_mem( -nob3*nobnvb )
       end if
       if( allocated(Mol%dijkaBA) ) then
          deallocate( Mol%dijkaBA )
          write(iout,"(' deallocated Mol%dijkaBA ')")  ;  call track_mem( -nob2*noanva )
       end if
       if( allocated(Mol%dijkaAB) ) then
          deallocate( Mol%dijkaAB )
          write(iout,"(' deallocated Mol%dijkaAB ')")  ;  call track_mem( -noa2*nobnvb )
       end if
       if( allocated(Mol%dijkaAA) ) then
          deallocate( Mol%dijkaAA )
          write(iout,"(' deallocated Mol%dijkaAA ')")  ;  call track_mem( -noa3*noanva )
       end if
       if( allocated(Mol%diajbBB) ) then
          deallocate( Mol%diajbBB ) 
          write(iout,"(' deallocated Mol%diajbBB ')")  ;  call track_mem( -nobnvb*nobnvb )
       end if
       if ( allocated(Mol%dijklBB) ) then
          deallocate( Mol%dijklBB )
          write(iout,'(" deallocated Mol%dijklBB ")')  ;  call track_mem( -nob3*nob3 )
       end if
       if( allocated(Mol%diajbBA) ) then
          deallocate( Mol%diajbBA )
          write(iout,'(" deallocated Mol%diajbBA ")')  ;  call track_mem( -nva2*nob2 )
       end if
       if( allocated(Mol%diajbAB) ) then
          deallocate( Mol%diajbAB )
          write(iout,'(" deallocated Mol%diajbAB ")')  ;  call track_mem( -nvb*noa2 )
       end if
       if( allocated(Mol%dijklAB) ) then
          deallocate( Mol%dijklAB )
          write(iout,'(" deallocated Mol%dijklAB ")')  ;  call track_mem( -nob2*noa2 )
       end if
       if( allocated(Mol%diajbAA) ) then
          deallocate( Mol%diajbAA )
          write(iout,'(" deallocated Mol%diajbAA ")')  ;  call track_mem( -noanva*noanva )
       end if
       if( allocated(Mol%dijklAA) ) then
          deallocate( Mol%dijklAA )
          write(iout,'(" deallocated Mol%dijklAA ")')  ;  call track_mem( -noa3*noa3 )
       end if
       if( allocated(Mol%dijabBB) ) then
          deallocate( Mol%dijabBB )
          write(iout,'(" deallocated Mol%dijabBB ")')  ;  call track_mem( -nvb3*nob3 )
       end if
       if( allocated(Mol%dijabAB) ) then
          deallocate( Mol%dijabAB ) 
          write(iout,'(" deallocated Mol%dijabAB ")')  ;  call track_mem( -nobnvb*noanva )
       end if
       if( allocated(Mol%dijabAA) ) then
          deallocate( Mol%dijabAA )
          write(iout,'(" deallocated Mol%dijabAB ")')  ;  call track_mem( -nva3*noa3 )
       end if
       
    case( '1e_int' )
    
       ntt = nbasis*(nbasis+1)/2

       if( allocated(Mol%soczao)  ) then
          deallocate(Mol%socxao, Mol%socyao, Mol%soczao)
          write(iout,'(A)') ' deallocated Mol%soczao, Mol%socyao, Mol%socxao '  
          call track_mem( -3*ntt )
       end if
       if( allocated(Mol%cmo_b) ) then
          deallocate(Mol%cmo_b)
          write(iout,'(A)') ' deallocated Mol%cmo_b '  
          call track_mem( -nrorb*nbasis )
       end if
       !: Need this to convert density back to AO basis
       !if( allocated(Mol%cmo_a) ) then
       !   deallocate(Mol%cmo_a) 
       !   write(iout,'(A)') ' deallocated Mol%cmo_a '
       !   call track_mem( -nrorb*nbasis )
       !end if
       if( allocated(Mol%orben) ) then
          deallocate(Mol%orben)
          write(iout,'(A)') ' deallocated Mol%orben '
          call track_mem( -norb )
       end if
       if( allocated(Mol%dipzao)) then
          deallocate(Mol%dipzao, Mol%dipyao, Mol%dipxao)
          write(iout,'(A)') ' deallocated Mol%dipzao, Mol%dipyao, Mol%dipxao '
          call track_mem( -3*ntt )
       end if
       if( allocated(Mol%vabsao)) then
          deallocate(Mol%vabsao)   
          write(iout,'(A)') ' deallocated Mol%vabsao '
          call track_mem( -ntt )
       end if
       if( allocated(Mol%dipzmob) ) then
          deallocate(Mol%dipzmob, Mol%dipymob, Mol%dipxmob) 
          write(iout,'(A)') ' deallocated Mol%dipzmob, Mol%dipymob, Mol%dipxmob '
          call track_mem( -3*nrorb*nrorb )
       end if
       if( allocated(Mol%dipzmoa) ) then
          write(iout,'(A)') ' DEALLOCATE Mol%dipzmoa, Mol%dipymoa, Mol%dipxmoa '
          flush(iout)
          deallocate(Mol%dipzmoa,Mol%dipymoa,Mol%dipxmoa) 
          write(iout,'(A)') ' deallocated Mol%dipzmoa, Mol%dipymoa, Mol%dipxmoa '
          flush(iout)
          call track_mem( -3*nrorb*nrorb )
       end if
       if( allocated(Mol%vabsmob) ) then
!:          deallocate(Mol%vabsmob)
!:          write(iout,'(A)') ' deallocated Mol%vabsmob '
!:          call track_mem( -nrorb*nrorb )
       end if
       if( allocated(Mol%vabsmoa) ) then
!:          deallocate(Mol%vabsmoa)
!:          write(iout,'(A)') ' deallocated Mol%vabsmoa '
!:          call track_mem( -nrorb*nrorb )
       end if
       if( allocated(Mol%socmoAB) ) then
          write(iout,*) " deallocating Mol%socmoAB,Mol%socmoBB,Mol%socmoAAMol%socmoAB,Mol%socmoBB,Mol%socmoAA"
          deallocate(Mol%socmoAB,Mol%socmoBB,Mol%socmoAA) 
          write(iout,'(A)') ' deallocated Mol%socmoAB, Mol%socmoBB, Mol%socmoAA'
          call track_mem( -3*nrorb*nrorb )
          flush(iout)
       end if
       if( allocated(Mol%field_env) ) then  
          deallocate(Mol%field_env)
          write(iout,'(A)') ' deallocated Mol%field_env'
          flush(iout)
       end if

    end select

    call write_header( 'deallocate_main','initialize','leave' )
          flush(iout)


  end subroutine deallocate_main
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_HAMDATA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_hamdata
    

    !: read in rest of tdci.dat to build the Hamiltonian


    use readintegrals
    implicit none


    integer(8) :: ii, jj
    
    
    call write_header( 'read_hamdata','initialize','enter' )
    

    call cpu_time(start)    
        

    call read_vabs_ao( Qread_vabs_ao ) !: read CAP elements in AO basis
    call read_dipx_ao( Qread_dipx_ao ) !: read Dipole-X in AO basis        
    call read_dipy_ao( Qread_dipy_ao ) !: read Dipole-Y in AO basis        
    call read_dipz_ao( Qread_dipz_ao ) !: read Dipole-Z in AO basis     
    call read_socx_ao( Qread_socx_ao ) !: read del x terms in AO basis    
    call read_socy_ao( Qread_socy_ao ) !: read del y terms in AO basis    
    call read_socz_ao( Qread_socz_ao ) !: read del z terms in AO basis    
    call read_orben( Qread_orben )     !: read in MO orbital energies     
    call read_cmo( Qread_cmo )         !: read in mo-lcao coefficients 


    Qread_buckets(2) = .True.
    Qread_buckets(5) = .True.
    if( unrestricted ) Qread_buckets(10) = .True. 
    if( trim(jobtype).eq.flag_soc )    Qread_buckets(1:8)  = .True.
    if( trim(jobtype).eq.flag_ip )     Qread_buckets(1:14) = .True.
    if( trim(jobtype).eq.flag_socip )  Qread_buckets(1:14) = .True.
    if( trim(jobtype).eq.flag_cisd )   Qread_buckets(1:21) = .True.


    if( trim(jobtype).eq.flag_tda ) then
       Qread_buckets = .False.
       call read_custom_ham
    end if
    
    
    !: read in integrals
    if( Qread_buckets(2) )  call read_bucket_2  !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB
    if( Qread_buckets(5) )  call read_bucket_5  !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB       
    if( Qread_buckets(10) ) call read_bucket_10 !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB
    if( Qread_buckets(1) )  call read_bucket_1  !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB
    if( Qread_buckets(3) )  call read_bucket_3  !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB          
    if( Qread_buckets(7) )  call read_bucket_7  !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB       
    if( Qread_buckets(8) )  call read_bucket_8  !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB       
    if( Qread_buckets(4) )  call read_bucket_4  !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L,IJ.LE.KL  IJKL       
    if( Qread_buckets(6) )  call read_bucket_6  !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL       
    if( Qread_buckets(9) )  call read_bucket_9  !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL       
    if( Qread_buckets(11) ) call read_bucket_11 !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA       
    if( Qread_buckets(12) ) call read_bucket_12 !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA       
    if( Qread_buckets(13) ) call read_bucket_13 !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA       
    if( Qread_buckets(14) ) call read_bucket_14 !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA       
    if( Qread_buckets(15) ) call read_bucket_15 !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC       
    if( Qread_buckets(16) ) call read_bucket_16 !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC
    if( Qread_buckets(17) ) call read_bucket_17 !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC       
    if( Qread_buckets(18) ) call read_bucket_18 !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC       
    if( Qread_buckets(19) ) call read_bucket_19 !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD       
    if( Qread_buckets(20) ) call read_bucket_20 !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD       
    if( Qread_buckets(21) ) call read_bucket_21 !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD
    
    
    call cpu_time(finish)    
    
    
    close(10)
    write(iout,"(' TDCI.dat read time:',f12.4,' seconds')") finish-start
    write(iout,'(A)') ' closed file '//trim(tdcidatfile)

    
    call write_header( 'read_hamdata','initialize','leave' )
    
    
  end subroutine read_hamdata
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITE_INPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine write_input(myoption)


    implicit none
    integer(8), intent(in) :: myoption

    integer(8) :: i
    

    if( myoption.ne.iout ) open(unit=myoption, file='input')

    write(myoption,'(A)') ' &DYNAMICS'
    write(myoption,"(' init_states(0) =',i2)") init_states(0)
    do i=1, init_states(0) 
       write(myoption,"(' init_states(', i1 , ') =',i2)") i, init_states(i)
    end do
    do i=1, init_states(0)
       write(myoption,"( ' init_coeffs(', i1, ') = ','(',f7.5,',',f7.5,')' )") i, real(init_coeffs(i)), aimag(init_coeffs(i)) 
    end do
    if ( restart )      write(myoption,'(A)') ' restart        = .True.  '
    if ( .not.restart ) write(myoption,'(A)') ' restart        = .False. '    
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') ' &FIELD_units'
    write(myoption,'(A)') " omega_units = "//"'"//trim(omega_units)//"'"
    !write(myoption,'(A)') " angle_units = "//"'"//trim(angle_units)//"'"
    write(myoption,'(A)') ' /'
    
    write(myoption,'(A)') ' &FIELD'
    write(myoption,"(' dirform  =', a8  )") adjustr("'"//trim(dirform)//"'")
    write(myoption,"(' ellipt   =', f7.3)") ellipt
    write(myoption,"(' envelope =', a7  )") adjustr("'"//trim(envelope)//"'")
    write(myoption,"(' ncyc     =', i7  )") ncyc
    write(myoption,"(' omega    =', f7.3)") omega
    write(myoption,"(' phase    =', f7.3)") phase
    write(myoption,"(' euler    =', f7.3)") euler
    write(myoption,'(A)') ' /'
    
    write(myoption,'(A)') '&FIELD_strengths '
    write(myoption,"(' nemax = ', i7)") nemax
    do i=1, nemax
       write(myoption,"(' read_emax(',i1,') = ',f10.4 )") i, read_emax(i)
    end do
    do i=1, nemax
      if(read_state1(i).ne.0)  write(myoption,"(' read_state1(',i1,') = ',i5 )") i,read_state1(i)
      if(read_state1(i).ne.0)  write(myoption,"(' read_coeff1(',i1,') = ',2f12.8 )") i,read_coeff1(i)
      if(read_state2(i).ne.0)  write(myoption,"(' read_state2(',i1,') = ',i5 )") i,read_state2(i)
      if(read_state2(i).ne.0)  write(myoption,"(' read_coeff2(',i1,') = ',2f12.8 )") i,read_coeff2(i)
    end do
    do i=1, nemax
      if(ion_sample_start(i).ne.0)  write(myoption,"(' ion_sample_start(',i1,') = ',i5 )") i,ion_sample_start(i)
      if(ion_sample_start(i).ne.0)  write(myoption,"(' ion_sample_width(',i1,') = ',i5 )") i,ion_sample_width(i)
      if(ion_sample_state(i).ne.0)  write(myoption,"(' ion_sample_state(',i1,') = ',i5 )") i,ion_sample_state(i)
    end do
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') '&FIELD_directions '
    write(myoption,"(' ndir = ', i7 )") ndir
    if( trim(dirform).eq.'cart' ) then
       do i=1, ndir
          write(myoption,"(' read_x(',i2,',') = ',f8.4 )") i, read_x(i)
          write(myoption,"(' read_y(',i2,',') = ',f8.4 )") i, read_y(i)
          write(myoption,"(' read_z(',i2,',') = ',f8.4 )") i, read_z(i)
       end do
    else 
       do i=1, ndir
          write(myoption,"(' read_theta(',i2,') = ',f10.4)") i, read_theta(i)
       end do
       do i=1, ndir
          write(myoption,"(' read_phi(',i2,')   = ',f10.4)") i, read_phi(i)
       end do
    end if
    
    write(myoption,'(A)') ' &SYSTEM_units'
    write(myoption,'(A)') " dt_units     = 'au' "
    write(myoption,'(A)') " eigmax_units = 'au' "
    write(myoption,'(A)') ' /'
    
    write(myoption,'(A)') ' &SYSTEM'
    write(myoption,"(' dt          =', f7.3)") dt
    write(myoption,"(' eigmax      =', f7.3)") eigmax
    write(myoption,"(' ionization  =', f7.3)") ionization
    write(myoption,"(' jobtype     =', a7 )") adjustr("'"//trim(jobtype)//"'")
    write(myoption,"(' nstep       =', i7  )") nstep
    write(myoption,"(' outstep     =', i7  )") outstep
    write(myoption,"(' nactive     =', i7  )") nactive
    write(myoption,"(' nvirtual    =', i7  )") nvirtual
    write(myoption,"(' socfac      =', f12.8)") socfac
    write(myoption,"(' socfacz     =', f12.8)") socfacZ
    write(myoption,"(' ffieldx     =', f9.5)") ffieldx
    write(myoption,"(' ffieldy     =', f9.5)") ffieldy
    write(myoption,"(' ffieldz     =', f9.5)") ffieldz
    if ( IP_alpha_beta )      write(myoption,'(A)') ' IP_alpha_beta  = .True.  '
    if ( .not.IP_alpha_beta ) write(myoption,'(A)') ' IP_alpha_beta  = .False.  '
    if ( QsocA2B )      write(myoption,'(A)') ' QsocA2B  = .True.  '
    if ( .not.QsocA2B ) write(myoption,'(A)') ' QsocA2B  = .False.  '
    if ( QeigenDC )      write(myoption,'(A)') ' QeigenDC  = .True.  '
    if ( .not.QeigenDC ) write(myoption,'(A)') ' QeigenDC  = .False.  '
    write(myoption,'(A)') ' /'        

    write(myoption,'(A)') ' &InOutFILES'
    if ( Qread_tdcidata )      write(myoption,'(A)') ' Qread_TDCIdata  = .True.  '
    if ( .not.Qread_tdcidata ) write(myoption,'(A)') ' Qread_TDCIdata  = .False.  '
    write(myoption,'(A)') ' tdcidatfile     = '//"'"//trim(tdcidatfile)//"'"
    write(myoption,'(A)') ' outputfile      = '//"'"//trim(outputfile)//"'"
    write(myoption,'(A)') ' restartbinfile  = '//"'"//trim(restartbinfile)//"'"
    if ( Qwrite_ion_coeff )   write(myoption,'(A)') ' Qwrite_ion_coeff  = .True.  '
    if ( Qread_ion_coeff )    write(myoption,'(A)') ' Qread_ion_coeff  = .True.  '
    if ( Qmo_dens )    write(myoption,'(A)') ' Qmo_dens  = .True.  '
    if ( Qci_save )    write(myoption,'(A)') ' Qci_save  = .True.  '
    if ( write_binaries )  write(myoption,'(A)') ' write_binaries = .True.   '
    write(myoption,'(A)') ' /'    

    write(myoption,'(A)') '&Davidson'
    if( flag_davidson ) write(myoption,'(A)') " flag_davidson = .True. "
    if ( .not. flag_davidson ) write(myoption,'(A)') " flag_davidson = .False. "    
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') ' &ReadU_NO'
    if( flag_ReadU_NO )      write(myoption,'(A)') " flag_ReadU_NO = .True. "
    if( .not.flag_ReadU_NO ) write(myoption,'(A)') " flag_ReadU_NO = .False. "
    write(myoption,'(A)') ' /'

    if( myoption.ne.iout ) close(myoption)
    
    
  end subroutine write_input
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_MODREADIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_modreadin(myroutine,option)
    
    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option

    
    select case(trim(myroutine))
    !: subroutine read_input
    case('read_input')
       select case( trim(option) ) 
       case( 'greeting' ) 
          write(iout,'(A)') ' THIS IS A WONDERFUL DAY, I HAVE NEVER SEEN THIS ONE BEFORE'
          write(iout,'(A)') ' - Maya Angelous'
          write(iout,'(A)') divide
       case( 'date' )          
          write(iout,'(A)') ' I WAS COMPILED ON Tue Jul  9 02:31:42 PM EDT 2024 '
          write(iout,'(A)') ' I AM A REVISED CODE FOR CW PULSE GENERATION '
          write(iout,'(A)') ' RAMPING PARAMETER SET TO RAMP_STEP=16000, NOT NSTEP'
          call dnt(iout)          
          call write_header( 'read_input','initialize','enter' )
       end select

       !: subroutine set_variables
    case('set_variables')
       select case( trim(option) )
       case('routine') 
          write(iout,'(A)') ' job route'
          write(iout,"(' read_input            ',l1)") Qread_input
          write(iout,"(' read_tdcidata         ',l1)") Qread_tdcidata
          write(iout,"(' set_variables         ',l1)") Qset_variables
          write(iout,"(' get_1eham             ',l1)") Qget_1eham
          write(iout,"(' get_expVabs           ',l1)") Qget_expVabs
          write(iout,"(' read_restart_binaries ',l1)") Qread_binaries
          write(iout,"(' save_ham              ',l1)") Qsave
          write(iout,"(' deallocate_main       ',l1)") Qdealloc
          write(iout,"(' propagate             ',l1)") Qpropagate
          write(iout,"(' serial                ',l1)") Qserial
       end select
       
    end select
    
    flush(iout)
    
100 format(' allocated rank 1 array',a15,' of length = ',i11)

    
  end subroutine writeme_modreadin
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ERRORS_MODREADIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_modreadin(myroutine,option)

    implicit none
    
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option


    select case( trim(myroutine) )

    !: subroutine read_tdcihead
    case( 'read_tdci1' )
       select case( trim(option) )
       case('no_tdcidat')
          write(iout,'(A)') " ERROR: Cannot find setup file 'TDCI.dat' "
          go to 100
       end select

    !: subroutine set_variables
    case( 'read_input' )
       select case( trim(option) ) 
       case('omega_error') 
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         omega_units must be in 'au' or 'nm' "
          go to 100 
       case('dt_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         dt_units must be in 'au' or 'ps' or 'fs' or 'as' "
          go to 100
       case('euler_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         euler_units must be in 'rad' or 'deg' "
          go to 100
       case('phase_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         euler_units must be in 'rad' or 'deg' "
          go to 100          
       case('eigmax_error')
          write(iout,'(A)') " WARNING:  OVERWRITING eigmax to 10.0 au "
          eigmax  = 10.d0
          go to 200
       end select

    !: subroutine set_variables
    case( 'set_variables' )
       select case( trim(option) )
       case('no_envelope')
          write(iout,'(A)') ' WARNING: envelope type not found'
          write(iout,'(A)') ' WARNING: overriding user-defined pulse shape to static pulse'
          go to 200 
       case('nstep_error') 
          write(iout,'(A)') ' WARNING: nstep is less than 1'
          write(iout,"(' WARNING: overriding user-defined total number of timesteps to: ',i0)") nstep
          go to 200
       case( 'norestart' )
          write(iout,'(A)') ' ERROR:  could not find restart binary file'
          go to 100
       end select
       
    end select
    
100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " DON'T BE DISMAYED BY GOODBYES - Richard Bach "
    stop
    
200 continue
    
  end subroutine errors_modreadin
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module initialize
