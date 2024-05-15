module variables_control 


  !: module variables_control
  !: controls the flow of the program 
  

  implicit none


  !: method 
  character(10) :: jobtype       !: 'cis', 'read', 'tda', 'ip'
  logical       :: unrestricted  !: restricted or unrestricted flag 
  logical       :: flag_davidson !: turn davidson on
  logical       :: flag_ReadU_NO      !: Read in test orbitals before propagation
  
  
  !: if envelope > 6 then circular pulse
  logical    :: linear = .true.

  
  !: subroutine route control
  logical :: &
       QRead_TDCIdata = .false., &
       Qread_input    = .true. , &
       Qread_tdcihead = .true. , &
       Qset_variables = .true. , &
       Qallocate_main = .true. , &
       Qread_hamdata  = .true. , &
       Qget_nstuse    = .true. , &
       Qgetfield     = .true.,  &
       Qgetham0      = .true. , &
       Qget_1eham     = .true. , &
       Qread_binaries = .false., &
       Qget_expVabs   = .true. , &
       Qsave          = .true. , &
       Qdealloc       = .true. , &
       Qpropagate     = .true. , &
       Qserial        = .false., &
       QeigenDC       = .true.

  
  !: for IP and SOCIP
!#  logical :: IP_alpha_beta = .true.
!# Attempted fix for less states
!# Include onlu beta ionizations
  logical :: IP_alpha_beta = .false.
!# Include alpha to beta and beta to apla excitations for SOC
  logical :: QsocA2B = .true.

  
  !: write files/specifics control
  logical :: &
       Qwrite_fshape      = .true., & 
       Qwrite_mo_energies = .true., &
       Qwrite_specifics1  = .true., &
       Qwrite_ham0        = .false., &
       Qwrite_mo_elements = .true., &
       Qwrite_specifics2  = .true. 
       
  !: write binaries control
  logical :: &
       Qwrite_vec = .true., &      !: write CI vectors
       Qwrite_tdx = .true., &      !: write x dipole matrix in CI basis
       Qwrite_tdy = .true., &      !: write y dipole matrix in CI basis
       Qwrite_tdz = .true., &      !: write z dipole matrix in CI basis
       Qwrite_abp = .true., &      !: write Vabs matrix in CI basis
       Qwrite_expabp = .true., &   !: write [R][exp(-Vabs_diag*dt/2)][R]^T 
       Qwrite_ion_coeff=.false., & !: write coefficients for singly ionized states from cis and soc simulations
       Qread_ion_coeff=.false., &  !: read coefficients from singly ionized states from cis and soc simulations
       Qmo_dens = .false., &       !: write density matrix in the MO basis
       Qci_save = .false.          !: write CI vectors to a binary file

  logical :: write_binaries !: Dump many binary arrays in matrices/ folder 
                            !: (This should be turned off for multithreaded calculations)
 
  
  !: allocate control
  logical :: &
       Qalloc_indices  = .True., &
       Qalloc_vabs     = .True., &
       Qalloc_exp_vabs = .True., &
       Qalloc_tdxyz    = .True., &
       Qalloc_field    = .True., &
       Qalloc_vabsmo   = .True., &
       Qalloc_dipmo    = .True., &
       Qalloc_socmo    = .False., &
       Qalloc_Zcomplex = .False.

  
  !: allocate building arrays
  logical :: &
       Qread_vabs_ao = .True.,  &
       Qread_dipx_ao = .True.,  &
       Qread_dipy_ao = .True.,  &
       Qread_dipz_ao = .True.,  &
       Qread_socx_ao = .False., &
       Qread_socy_ao = .False., &
       Qread_socz_ao = .False., &
       Qread_orben   = .True.,  &
       Qread_cmo     = .True.
  

  logical :: Qread_buckets(21) = .False.
       
  
end module variables_control
