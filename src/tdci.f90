!--------------------------------------------------------------------------------!
! FORTRAN code to propagate a molecule in an electric field                      !
! using Absorbing Boundary Potential to calculate ionization rates               !
!                                                                                !
! Initial version written by J.Sonk, Sept 2012                                   !
! Extensions/modifications by P.Krause, H.B.Schlegel                             !
! Restricted and unrestricted CIS, parallel implementation by H.B.Schlegel       !
!--------------------------------------------------------------------------------!
program main

  
  use variables_global !: inherits units, variables_setup, variables_control
  use initialize
  use readintegrals
  use write_info
  use getfield
  use getham          !: inherits getham0
  use propagate
  use Zpropagate 
  use davidson_ip  !: only to check Hamiltonian

  
  implicit none
  integer(8) :: i,j,a,b,ia, istate, iout2 !: temporary cisd variables
  real(8)  :: masterstart, masterfinish
  

  !: output opened in read_input with iout=42   
  call cpu_time(masterstart)
  
  !: Allocate MolInfo class instance
  allocate(Mol)
  
  !: read and set input parameters (module initialize) 
  if( Qread_input       ) call read_input  

  !: read matrix elements file to get all the integrals and data 
  if(.not.Qread_tdcidata) call read_matrixelements

  !: alteratively, read TDCI.dat file 
  if( Qread_tdcidata    ) call read_tdcihead
  if( Qset_variables    ) call set_variables
  

  !: allocate main arrays before allocating temporary arrays (module initialize)
  if ( Qallocate_main ) then
     call allocate_main( 'main' )
     call allocate_main( 'mo' )
  end if
  
  
  !: read 1-electron and 2-electron matrix elements from TDCI.dat file (module initialize)
  if( Qread_hamdata ) then
     if( Qread_tdcidata     ) call read_hamdata
     if( Qwrite_mo_energies ) call write_mo_energies
  end if
  
  !: generate field 
  !: read polarization directions, dirform<0 in polar ; dirform>0, in Cartesian
  !: for circularly polarized pulses, two vectors perpendicular to propdirection   
  if( Qgetfield ) then
     call check_emax
     call get_lindirection
     if (.not.linear) call get_circdirection
     call shape_field
     if ( Qwrite_fshape ) call write_field_shape
  end if
  

  !: write out useful information I
  if( Qwrite_specifics1 ) call write_specifics1
 
  !: form CIS Hamiltonian and diagonalize 
  !: matrix elements stored in cis_vec since cis_vec will be fed into dysev diagonalization
  !: dysev will spit out eigenvectors into cis_vec.
  if ( Qgetham0 ) then
     select case ( trim(jobtype) ) 
     case( flag_cis ) 
        call get_cisN
        cis_vec(1)=-500.d0
        call diagonalize
        cis_eig(1)=0.d0     
        call get_sip_cis
     case( flag_soc ) 
        call get_soc_cisN
        call Zdiagonalize
        call get_sip_soc
     case( flag_tda )
        call get_cis_index
        cis_vec(1)=-500.d0
        call diagonalize
        cis_eig(1)=0.d0   
     case( flag_ip  ) 
        call get_ip_cisd
        call diagonalize
        call get_dip_cis
     case( flag_socip ) 
        call get_socip_cisd
        call Zdiagonalize
        call get_dip_soc
     case( flag_cisd ) 
        call get_cisd_index
        call diagonalize         
     end select

     ! We need the 2e integrals for IP version of eq18
     !call deallocate_main( '2e_int' )     
     if( Qwrite_ham0 ) call write_ham0 
  end if

  !: davidson
  if ( flag_davidson ) then
     select case( trim(jobtype) ) 
     case( flag_ip )
        allocate( hole_index(nstates,2) ) ; hole_index = 0
        allocate( part_index(nstates,1) ) ; part_index = 0
        call get_ip_index
        call davidson_ip
     case( flag_cisd )
        write(iout,'(A)') " Getting there "
        stop
     end select
  end if

  !: set NstUse = number of states to use <= nstates
  if( Qget_nstuse ) call get_nstuse  
  
  !: Transform ONE electron integrals from AO to MO basis (VabsAO,Dip*AO-->VabsMO,Dip*MO)
  !: Transform from MO basis to CIS states (abp,TDX,TDY,TDZ)    
  !: Calls get_ao2mo --> get_form1det --> form1cis, OMP parallelized
  !: Or read in the binaries if restart = .True.
  if( Qget_1eham ) then
     call get_1eham
     if( Qwrite_mo_elements ) call write_mo_elements
  end if

 
  if( Qread_binaries ) call read_restart_bin
  
  !: write out useful information II
  if( Qwrite_specifics2 ) call write_specifics2
 
  !: compute exp(-Vabs *dt/2) if restart = .False.
  if( Qget_expVabs ) call get_expVabs

  
  !: save for later restart.  
  if( Qsave ) call save_restart_bin

  !: clean matrices folder, save Vabs AO 
  call cleanup_directory("matrices")
  call write_dbin( Mol%vabsao, nbasis*(nbasis+1)/2, "matrices/Vabs_AO.bin")
  

  write(iout,*) ' Qdealloc',Qdealloc
  !: deallocate un-used arrays
  if( Qdealloc ) call deallocate_main( '1e_int' )
  
  !: PROPAGATE 
  if ( Qpropagate ) then
     if( linear ) then
        select case ( trim(jobtype) )
        case( flag_cis ) ;    call trotter_linear
        case( flag_ip  ) ;    call trotter_linear
        case( flag_soc ) 
          if(Qserial) then
            call Ztrott_serial_lin
          else
            call Ztrotter_linear
          end if
        case( flag_socip ) 
          if(Qserial) then
            call Ztrott_serial_lin
          else
            call Ztrotter_linear
          end if
        end select
     else
        select case ( trim(jobtype) )
        case( flag_cis ) ;    call trotter_circular
        case( flag_ip  ) ;    call trotter_circular
        case( flag_soc ) 
          if(Qserial) then
            call Ztrott_serial_cir
          else
            call Ztrotter_circular
          end if
        case( flag_socip ) 
          if(Qserial) then
            call Ztrott_serial_cir
          else
            call Ztrotter_circular
          end if
        end select
     end if
     flush(iout)
     call write_summary
  end if
  

  call cpu_time(masterfinish)  

  
  write(iout,'(A)') divide
  write(iout,'(A)') ' '
  write(iout,"(' ---> total elapsed time:',f12.4,' s')") masterfinish-masterstart
  write(iout,'(A)') ' '
  call dnt(iout)
  write(iout,'(A)') ' '
  write(iout,'(A)') " MY MISSION IN LIFE IS NOT MERELY TO SURVIVE, BUT TO THRIVE;"
  write(iout,'(A)') " AND TO DO SO WITH SOME PASSION,"
  write(iout,'(A)') " SOME COMPASSION, SOME HUMOR, AND SOME STYLE"
  write(iout,'(A)') ' -Maya Angelou'
  flush(iout)

  close(iout)
  

end program main
