module propagate
  
use variables_global
use analysis
use util
use sort
use io_binary  ! io_bin_test, write_dbin, read_dbin
  
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

  !: file stuff
  integer(8)   :: funit(6) !: Output data file units
  character(4) :: dirstr, emaxstr
  character(100) :: cifile, datafile
  !: lapack stuff and misc
  integer(8) :: info2, info3
  real(8)    :: start1, start2, finish1, finish2
  !integer(8) :: info1, lscratch, liwork
  !integer(8), allocatable :: iwork(:)
  !real(8), allocatable    :: scratch(:)

  integer(8) :: &
    nva95max, nva99max, nva95maxMO, nva99maxMO, nva95max_debug, &
    nva99max_debug, nva95MO_sort, nva99MO_sort, nva95NO_sort, nva99NO_sort, &
    nva95max_direct, nva99max_direct, MO99, NO99

  real(8) :: rate_density, rate_direct, rate_debug

  contains
    procedure :: initialize => initialize_PropPriv
    procedure :: deconstruct => deconstruct_PropPriv
    !procedure :: write_outdata => write_outdata_PropPriv


end type PropagationPrivate

    

type PropagationShared

  integer(8) :: &
    nva95maxmax, nva99maxmax, nva95maxmaxMO, nva99maxmaxMO, &
    nva95maxmax_direct,nva99maxmax_direct,nva95maxMO_sort,nva99maxMO_sort, &
    nva95maxNO_sort, nva99maxNO_sort

  !: Natural Orbital generation
  real(8), allocatable :: opdm_avg(:)  !: Averaged one-particle reduced density matrix (1-RDM)
  real(8), allocatable :: opdm_avg_abs(:)  !: average(abs(opdm))
  real(8), allocatable :: natorb_occ(:) !: Natural orbital occupations (eigenvalues)
  real(8), allocatable :: natorb_occ_abs(:) !: Natural orbital occupations (eigenvalues)
  real(8), allocatable :: U_NO(:) !: MO->NO transformation matrix
  real(8), allocatable :: U_NO_abs(:) !: MO->NO_abs transformation matrix
  integer(8) :: opdm_avg_N, ndim
  real(8), allocatable :: U_NO_input(:) !: U_NO read from file

  real(8) :: ratemax_density, ratemax_direct, ratemax_debug
  

  contains
    procedure :: initialize => initialize_PropShared

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
  allocate( this%rate_a(noa), this%rate_b(nob) )
  this%rate_aa = 0.d0 ; this%rate_ab = 0.d0 ; this%rate_ba = 0.d0 ; this%rate_bb = 0.d0
  this%rate_a = 0.d0 ; this%rate_b = 0.d0

  this%dirx1 = 0.d0 ; this%diry1 = 0.d0 ; this%dirz1 = 0.d0
  this%dirx2 = 0.d0 ; this%diry2 = 0.d0 ; this%dirz2 = 0.d0
  this%emax1 = 0.d0; this%efield1 = 0.d0; this%emax2 = 0.d0; this%efield2 = 0.d0
  this%efieldx = 0.d0 ; this%efieldy = 0.d0; this%efieldz = 0.d0
  this%temp = 0.d0; this%temp1 = 0.d0; this%temp2 = 0.d0
  this%norm = 0.d0; this%normV = 0.d0; this%rate = 0.d0
  this%mux = 0.d0; this%muy = 0.d0 ; this%muz = 0.d0

  this%nva95max = 0 ; this%nva99max = 0 ; this%nva95maxMO = 0 ; this%nva99maxMO = 0
  this%nva95max_direct = 0 ; this%nva99max_direct = 0
  this%nva95max_debug = 0 ; this%nva99max_debug = 0
  this%nva95MO_sort = 0
  this%nva99MO_sort = 0
  this%nva95NO_sort = 0
  this%nva99NO_sort = 0    
  this%rate_density = 0.d0 ; this%rate_direct = 0.d0 ; this%rate_debug = 0.d0

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

  !write(iout, *) "Deallocating Priv!"
  deallocate( this%rate_aa(noa*noa), this%rate_ab(noa*nob), this%rate_ba(nob*noa), this%rate_bb(nob*nob) )
  deallocate( this%rate_a(noa), this%rate_b(nob) )

end subroutine deconstruct_PropPriv

!: Run once for each thread at start of propagation
!: Sets up the funit values and writes headers to .dat files.
subroutine PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, Zion_coeff, MODE )
  implicit none

  class(PropagationPrivate), intent(inout) :: Priv
  integer(8), intent(in) :: iemax, idir
  type(tdcidat), allocatable, intent(inout) :: tdciresults(:)
  complex(8), intent(inout) :: psi0(:), Zion_coeff(:)
  !: MODE 0=trotter_linear, 1=trotter_circular
  integer, intent(in) :: MODE


  !: Local variables

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

  !: cifile binary
  if( Qci_save ) then
    Priv%funit(1)  = iemax*100 + idir
    Priv%cifile ='CI-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.bin'
    open( unit=Priv%funit(1), file=trim(Priv%cifile), form='unformatted' )
    write(Priv%funit(1)) ndata, nstuse, nstates
    select case (MODE)
      case(0) !: trotter_linear
        write(Priv%funit(1)) 0.d0, 0.d0, 0.d0, Priv%dirx1, Priv%diry1, &
          Priv%dirz1, 0.d0, 0.d0, 0.d0, 1.d0
      case(1) !: trotter_circular
        write(Priv%funit(1)) 0.d0, 0.d0, Priv%dirx1, Priv%diry1, Priv%dirz1, Priv%dirx2, Priv%diry2, Priv%dirz2, 1.d0
    end select
    write(Priv%funit(1)) real(psi0)
    write(Priv%funit(1)) aimag(psi0) 
    flush(Priv%funit(1))      
  end if    
  
  !: RESULTS datafile
  Priv%funit(2) = 1000+100*iemax + idir
  Priv%datafile = 'RESULTS-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(2),file=trim(Priv%datafile) )
  write( Priv%funit(2), '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'

  select case (MODE)
    case(0) !: trotter_linear
      write( Priv%funit(2), "( '#',  20(f16.10,1x) )" ) Priv%emax1, 0.d0, &
             tdciresults(idir+(iemax-1)*ndir)%theta0, tdciresults(idir+(iemax-1)*ndir)%phi0,&
             0.d0,0.d0,  0.d0,0.d0, &
             Priv%dirx1, Priv%diry1, Priv%dirz1,  0.d0,0.d0,0.d0,  0.d0,0.d0,0.d0
    case(1) !: trotter_circular
      write( Priv%funit(2), "( '#',  20(f16.10,1x) )" ) Priv%emax1, Priv%emax2, &
           tdciresults(1+(iemax-1)*ndir)%theta0, tdciresults(1+(iemax-1)*ndir)%phi0, &
           tdciresults(1+(iemax-1)*ndir)%theta1, tdciresults(1+(iemax-1)*ndir)%phi1, &
           tdciresults(1+(iemax-1)*ndir)%theta2, tdciresults(1+(iemax-1)*ndir)%phi2, &
           tdciresults(1+(iemax-1)*ndir)%x0, tdciresults(1+(iemax-1)*ndir)%y0, tdciresults(1+(iemax-1)*ndir)%z0, & 
           Priv%dirx1, Priv%diry1, Priv%dirz1,  Priv%dirx2, Priv%diry2, Priv%dirz2
  end select

  if( trim(jobtype).eq.flag_ip) then
    write( Priv%funit(2),"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)', & 
      'NO99 MO99 ','field1','field2','fieldx','fieldy','fieldz', &
      'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)','6 x |psi(i)|**2'
  else
    write( Priv%funit(2),"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)', &
      'NO99 MO99 ','field1','field2','fieldx','fieldy','fieldz', &
      'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
  end if

  !: POP datafile
  Priv%funit(3) = 2000+100*iemax + idir
  Priv%datafile = 'POP-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(3),file=trim(Priv%datafile) )
  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    write( Priv%funit(3),"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
      'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
  else
    write( Priv%funit(3),"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
      'rate(fs-1)','rate_a(occ)','rate_b(occ)'
  end if

  !: ION datafile
  Priv%funit(4) = 3000+100*iemax + idir
  Priv%datafile = 'ION-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
  open( unit=Priv%funit(4),file=trim(Priv%datafile) )
  write( Priv%funit(4),"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'

  !: ION_COEFF datafile
  If( Qread_ion_coeff .or. Qwrite_ion_coeff ) then
    Priv%funit(5) = 4000+100*iemax + idir
    Priv%datafile = 'ION_COEFF-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.bin'
    open( unit=Priv%funit(5),file=trim(Priv%datafile),form='unformatted' )
  else
    Priv%funit(5) = 4000+100*iemax + idir
    Priv%datafile = 'ION_COEFF-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
    open( unit=Priv%funit(5),file=trim(Priv%datafile) )
    write( Priv%funit(5),"(a5,2(a14,1x),3(a24,1x))" ) '#','time(fs)','rate/norm2', &
      's(i)(i=1,ip_states)','proj Zion_coeff(j,i)','Zion_coeff(j,i)'
    write(Priv%funit(5),"('ntimes ',i0)") int(nstep/outstep)
    write(Priv%funit(5),"('ip_states ',i0)") ip_states
  end If
  If( Qread_ion_coeff ) then 
    read(Priv%funit(5)) Zion_coeff
  end if

  !: MO density datafile
  If( Qmo_dens ) then
    Priv%funit(6) = 5000+100*iemax + idir
    Priv%datafile = 'MO_density-e'//trim(Priv%emaxstr)//'-d'//trim(Priv%dirstr)//'.dat'
    open( unit=Priv%funit(6),file=trim(Priv%datafile) )
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
  

  if( Qci_save) then
    write(Priv%funit(1)) dble(itime)*dt*au2fs, Priv%efield1, 0.d0, &
     Priv%dirx1, Priv%diry1, Priv%dirz1, 0.d0, 0.d0, 0.d0, Priv%norm**2
    write(Priv%funit(1)) real( psi )
    write(Priv%funit(1)) aimag( psi )
    flush(Priv%funit(1))
  end if


  if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
    if(Priv%norm.ne.0) then
      write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs,Priv%nva99max,Priv%nva99max_direct,Priv%efield1,0.d0,&
             Priv%efieldx,Priv%efieldy,Priv%efieldz, &
             Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz, &
             dble(dconjg(psi(1))*psi(1))/Priv%norm**2, &
             dble(dconjg(psi(2))*psi(2))/Priv%norm**2, &
             dble(dconjg(psi(3))*psi(3))/Priv%norm**2, &
             dble(dconjg(psi(4))*psi(4))/Priv%norm**2, &
             dble(dconjg(psi(5))*psi(5))/Priv%norm**2, &
             dble(dconjg(psi(6))*psi(6))/Priv%norm**2, &
             dble(dconjg(psi(7))*psi(7))/Priv%norm**2, &
             dble(dconjg(psi(8))*psi(8))/Priv%norm**2
    else
      write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
             idata, dble(itime)*dt*au2fs,Priv%nva99max,Priv%nva99max_direct,Priv%efield1,0.d0,&
             Priv%efieldx,Priv%efieldy,Priv%efieldz, &
             Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz,&
             0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
    end if
  else
    write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs,Priv%nva99max,Priv%nva99max_direct,Priv%efield1,0.d0,&
           Priv%efieldx,Priv%efieldy,Priv%efieldz, &
           Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz
  end if
  flush(Priv%funit(2))



  if( trim(jobtype).eq.flag_ip) then
    write( Priv%funit(3),"( i5,f10.4,500(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs, Priv%norm**2, &
           (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),Priv%rate/au2fs,&
           (Priv%rate_aa(i),i=1,noa*noa),&
           (Priv%rate_ab(i),i=1,noa*nob),(Priv%rate_ba(i),i=1,nob*noa),(Priv%rate_bb(i),i=1,nob*nob)
  else
    write( Priv%funit(3),"( i5,f10.4,500(1x,f15.10))") &
           idata, dble(itime)*dt*au2fs, Priv%norm**2, &
           (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob), &
           (Priv%rate_a(i),i=1,noa),(Priv%rate_b(i),i=1,nob)
  end if
  flush(Priv%funit(3))

  write( Priv%funit(4),"( i5,f10.4,500(1x,f15.10))") &
         idata, dble(itime)*dt*au2fs, Priv%rate/Priv%norm**2, Priv%normV, &
         (ion(i),i=1,noa),(ion(noa+nva+i),i=1,nob)
  flush(Priv%funit(4))
              
  if( Qmo_dens) then
    call write_density_difference( Priv%funit(6), dble(itime)*dt*au2fs, &
                                   Priv%rate, noa, nva, scratch, Mol%vabsmoa )
    !write(Priv%funit(6),"(f13.9,f16.10)") dble(itime)*dt*au2fs,Priv%rate
    !do i=1, noa+nva
    !  if ( abs(scratch(i+(i-1)*(noa+nva))).gt.1.d-6 ) then
    !    !write(Priv%funit(6),"(i5,i5,1x,f13.10,',')",advance='no')
    !    write(Priv%funit(6),"(i5,i5,1x,f13.10)") &
    !          i, i, scratch(i+(i-1)*(noa+nva))
    !  end if
    !  do j=(i+1), noa+nva
    !    if ( abs(scratch(i+(j-1)*(noa+nva))).gt.1.d-6 ) then 
    !      !write(Priv%funit(6),"(i5,i5,1x,f13.10,',')",advance='no') 
    !      write(Priv%funit(6),"(i5,i5,1x,f13.10)") &
    !            i, j, 2.d0*scratch(i+(j-1)*(noa+nva))
    !    end if
    !  end do
    !end do
    !write(Priv%funit(6),"(i5,i5,1x,f13.10)") 0,0,0.d0
  end if



  if(trim(jobtype).eq.flag_cis .and. (.not. Qwrite_ion_coeff) ) then
    call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                       psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                       Priv%rate,Zion_coeff,ion_coeff,scratch)
    call get_proj_ion(iout,noa+nob,ip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
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
  end if



  if(trim(jobtype).eq.flag_ip .and. (.not. Qwrite_ion_coeff) ) then
    call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                ip_states,state_ip_index,psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob, &
                Priv%rate,Zion_coeff,ion_coeff,scratch)
    call get_proj_ion(iout,ip_states,ip_vec,Zion_coeff(ip_states+1),Zproj_ion)
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
  end if


end subroutine PropWriteData


!==================================================================!
!  PropagationPrivate subroutines end
!==================================================================!

!==================================================================!
!  PropagationShared subroutines start
!==================================================================!
subroutine initialize_PropShared(Prop, ndim)
  implicit none

  class(PropagationShared), intent(inout) :: Prop

  integer(8), intent(in) :: ndim

  integer(8) :: info
  
  !: Natural Orbital initialization
  allocate( Prop%opdm_avg(ndim*ndim) )
  allocate( Prop%opdm_avg_abs(ndim*ndim) )
  allocate( Prop%natorb_occ(ndim))
  allocate( Prop%natorb_occ_abs(ndim))
  allocate( Prop%U_NO(ndim*ndim) )
  allocate( Prop%U_NO_abs(ndim*ndim) )
  allocate( Prop%U_NO_input(ndim*ndim) )
  Prop%opdm_avg = 0.d0
  Prop%opdm_avg_abs = 0.d0
  Prop%natorb_occ = 0.d0
  Prop%natorb_occ_abs = 0.d0
  Prop%U_NO = 0.d0
  Prop%U_NO_abs = 0.d0
  Prop%opdm_avg_N = ndir*(nstep/outstep)
  write(iout, '(A,i5)') "opdm_avg_N: ", Prop%opdm_avg_N

  !if ( .true. ) then
  if ( flag_ReadU_NO ) then
    call read_dbin( Prop%U_NO_input, (ndim*ndim), "U_NO_input.bin", info)
  end if

  Prop%nva95maxmax = 0 ; Prop%nva99maxmax = 0 ; Prop%nva95maxmaxMO = 0 ; Prop%nva99maxmaxMO = 0
  Prop%nva95maxmax_direct = 0 ; Prop%nva99maxmax_direct = 0 ; Prop%nva95maxMO_sort = 0
  Prop%nva99maxMO_sort = 0 ; Prop%nva95maxNO_sort = 0 ; Prop%nva99maxNO_sort = 0

end subroutine initialize_PropShared



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
  
  !: shared variables
  class(PropagationShared), allocatable :: Prop
  integer(8) :: nstuse2 
  complex(8) :: exphel(nstuse), psi0(nstuse)
  real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse) 
  

  !: private variables
  class(PropagationPrivate), allocatable :: Priv    
  !: Psi arrays must be not allocatable for performance
  complex(8) :: psi_j, psi_k, psi(nstuse), psi1(nstates)

  real(8) :: norm0
  real(8),allocatable :: pop0(:),pop1(:),ion(:)
  complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)

  integer(8) :: i, j, ii, jj, k, kk
  integer(8) :: itime, idir, iemax, ndim
  integer(8) :: ithread, idata
  complex(8) :: cdum

  real(8) :: start1, finish1, start3, finish3
  real(8) :: times(100)

  !: lapack stuff
  integer(8) :: info1, lscratch, liwork
  integer(8), allocatable :: iwork(:)
  real(8), allocatable    :: scratch(:)



  call write_header( 'trotter_linear','propagate','enter' )    
  call writeme_propagate( 'trot_lin', 'equation' ) 
  call cpu_time(start)
  

  !: nstuse2
  nstuse2 = nstuse * nstuse
  
  !: initialize psi0
  psi0 = dcmplx(0.d0,0.d0)
  do i=1, init_states(0)
    psi0( init_states(i) ) = init_coeffs(i)
  end do
  
  !: normalize
  call get_norm( norm0, nstuse, psi0 )
  psi0 = psi0 / norm0 

  !: write psi0
  call writeme_propagate( 'trot_lin', 'psi0' )    

  ndim = noa+nva
  !: Initialize shared variables
  allocate(Prop)
  call Prop%initialize(ndim)


  !: get initial population
  allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )    
  i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
  allocate(ion_coeff(i))
  If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
    allocate(Zion_coeff(nstep*ip_states*(ip_states+1)))
  else
    allocate(Zion_coeff(ip_states*(ip_states+1)))
    allocate(Zproj_ion(ip_states*ip_states))
  end If

  norm0 = 1.d0
  
  call get_psid( nstuse, nstates, cis_vec, norm0, psi0, psi_det0 )
  select case ( trim(jobtype) ) 
  case( flag_cis ) 
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
  end select
  call writeme_propagate( 'trot_lin','pop0', pop0, norb )

  deallocate( pop0 )
  

  if( QeigenDC ) then
    allocate( iwork(3+5*nstuse) )
    allocate( scratch(1+8*nstuse+2*nstuse*nstuse) )
  else
    allocate( iwork(2) )
    allocate( scratch(nstuse*nstuse) )
  end if

  !: exphel = exp(-iH*dt/2)
  do i=1, nstuse
     Priv%temp = -0.5d0 * cis_eig(i) * dt
     exphel(i) = dcmplx( dcos(Priv%temp) , dsin(Priv%temp) )
  end do

  !: Start loop over directions.  Counters need to be passed in as non-derived datatype

  !$OMP PARALLEL DEFAULT(NONE),&
  !$OMP PRIVATE(Priv, i, idata, idir, iemax, ii, itime, j, jj, k, kk, ithread, &
  !$OMP norm0, lscratch, liwork, scratch, iwork, info1, start1, finish1, &
  !$OMP psi_j, psi_k, psi, psi1, start3, finish3, times, &
  !$OMP pop1, ion, ion_coeff, psi_det0,  &
  !$OMP hp1, tdvals1, Zion_coeff ),  &
  !$OMP SHARED( Mol, Prop, jobtype, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
  !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
  !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
  !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
  !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
  !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
  !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
  !$OMP flag_ReadU_NO )
  
  !$OMP DO  

  dir_loop : do idir=1, ndir
     
    ithread = omp_get_thread_num()
    call cpu_time( start1 )

    psi = 0 ; psi1 = 0
    tdvals1 = 0.d0 ;

    !write(iout, '(A, I5, L1)') "Priv allocated at start of thread ", ithread , allocated(Priv)
    allocate(Priv)
    call Priv%initialize()

    !: get directions stored in TDCItdciresults
    Priv%dirx1 = tdciresults(idir)%x0  
    Priv%diry1 = tdciresults(idir)%y0  
    Priv%dirz1 = tdciresults(idir)%z0  

    !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
    hp1(:) = Priv%dirx1*tdx(1:nstuse2) &
           + Priv%diry1*tdy(1:nstuse2) + Priv%dirz1*tdz(1:nstuse2)

    !: diagonalize mu dot E
    info1 = 10  
    If( QeigenDC ) then
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

      call cpu_time(start1)
      !: get emax
      Priv%emax1 = tdciresults(1+(iemax-1)*ndir)%fstrength0

      call PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, Zion_coeff, 0)

      !$OMP CRITICAL
      !: all thread execute the code, but only one at a time
      if(iemax.eq.1) write(iout,"(12x,'TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
      write( iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0)" ) idir, iemax, ithread
      flush( iout )
      !$OMP END CRITICAL
                  
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
          !call get_norm( norm0, nstuse, psi )
          !write(iout,"(i5,16F10.6)") itime,norm0
          call cpu_time(start3)
           
          !: modified midpoint 
          Priv%efield1 = 0.5d0 * Priv%emax1 * ( fvect1(itime) + fvect1(itime+1) )
          if( read_shift(iemax).ne.0 ) then
            if(itime.eq.1) write(iout,"(' Pulse has been shifted by ',i6,' steps')") read_shift(iemax)
            i = itime-read_shift(iemax)
            if( i.gt.0 .and. i.lt.nstep ) then
              Priv%efield1 = 0.5d0 * Priv%emax1 * ( fvect1(i) + fvect1(i+1) )
            else
              Priv%efield1 = 0.d0
            end if
          end if
          Priv%efield2 = 0.d0
          Priv%efieldx = Priv%dirx1 * Priv%efield1
          Priv%efieldy = Priv%diry1 * Priv%efield1
          Priv%efieldz = Priv%dirz1 * Priv%efield1
           
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


          call cpu_time(start3)
          !: exp(-iHel dt/2 ) * psi
          do i=1, nstuse
            psi(i) = exphel(i) * psi(i)
          end do
          call cpu_time(finish3)
          times(3) = finish3 - start3

           
          call cpu_time(start3)
          !: W * exp(-Vabs dt/2) * psi
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
          !: exp(-E(t+dt/2)*mu*dt) * psi
          do j = 1, nstuse
            Priv%temp = dt * Priv%efield1 * tdvals1(j)
            psi1(j) = dcmplx( dcos(Priv%temp),dsin(Priv%temp) ) * psi1(j)
          end do
          call cpu_time(finish3)
          times(5) = finish3 - start3

           
          call cpu_time(start3)
          !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W * psi
          do j = 1, nstuse
            jj = nstuse * ( j-1 )
            psi_j = dcmplx( 0.d0, 0.d0 )
            do k=1, nstuse
               psi_j  = psi_j + hp1( jj+k ) * psi1(k)
            end do
              psi(j) = exphel(j) * psi_j
          end do
          call cpu_time(finish3)
          times(6) = finish3 - start3
           
          call cpu_time(start3)
          if( Qwrite_ion_coeff ) then
            call get_norm( Priv%norm, nstuse, psi )
            call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )
            ii = (itime-1)*(noa+nob)*(noa+nob+1)
            call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                     psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                     Priv%rate, Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,scratch)
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

            if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
              call pop_rate_ip(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                   hole_index,part_index,state_ip_index,ip_states, &
                   pop1,ion,ion_coeff,Priv%rate_aa,Priv%rate_ab,Priv%rate_ba,Priv%rate_bb, &
                   psi_det0,psi1,Priv%normV,Mol%vabsmoa,Mol%vabsmob,scratch,au2fs,Priv%rate_direct)
            else
              call pop_rate(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                   hole_index,part_index,state_ip_index,ip_states, &
                   pop1,ion,ion_coeff,Priv%rate_a,Priv%rate_b,psi_det0,psi1,Priv%normV, &
                   Mol%vabsmoa,Mol%vabsmob,unrestricted,scratch,au2fs,Priv%rate_direct)
            end if 
            !: 1-RDM (opdm) should be stored in scratch now.

            !: Check max nva for input NOs at this step.
            Priv%nva95max_direct = jj
            Priv%nva99max_direct = kk
            !: hmm why is this here
            !Priv%ratemax_direct = Priv%rate_direct
            call update_maxnva( Priv%nva95maxMO, Priv%nva99maxMO, Priv%nva95max, Priv%nva99max, &
                    Priv%nva95MO_sort, Priv%nva99MO_sort, Priv%nva95NO_sort, Priv%nva99NO_sort, &
                    Priv%nva95max_debug, Priv%nva99max_debug, Priv%rate_debug, &
                    Prop%opdm_avg, Prop%U_NO_input, Mol%vabsmoa )

            !$OMP CRITICAL
            !: Add 1-RDM (opdm) to 1-RDM average for Natural Orbital generation.
              
            call add_opdm_average( Prop%opdm_avg, scratch, Prop%opdm_avg_N )
            call add_opdm_average( Prop%opdm_avg_abs, scratch, &
                                   Prop%opdm_avg_N, .false., .true. )

            !call generate_natural_orbitals( scratch, U_NO, natorb_occ )
            !call NO_rate_sanity( scratch, U_NO, natorb_occ, Mol%vabsmoa, Mol%cmo_a)
            !flush(iout)    
            !call dgemm_sanity 

            !$OMP END CRITICAL

            call get_norm( Priv%norm,nstuse, psi )
            call get_expectation( nstuse, Priv%norm, psi, abp, Priv%rate) !: rate expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdx, Priv%mux ) !: mux  expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdy, Priv%muy ) !: muy  expectation value
            call get_expectation( nstuse, Priv%norm, psi, tdz, Priv%muz ) !: muz  expectation value
            !write(iout,*) " expectation Priv%rate", itime,Priv%norm,Priv%rate
            Priv%rate = -2.d0 * Priv%rate * Priv%norm**2


            call PropWriteData(Priv, psi, psi1, psi_det0, Zion_coeff,ion_coeff,&
                               scratch, itime, idata, pop1, ion)

          end if analysis
           
        end do timestep_loop
        call cpu_time(Priv%finish2)

        if( Qci_save ) close(Priv%funit(1))
        close(Priv%funit(2))
        close(Priv%funit(3))
        close(Priv%funit(4))

        if( Qwrite_ion_coeff ) write(Priv%funit(5)) Zion_coeff
        close(Priv%funit(5))

        if( Qmo_dens ) close(Priv%funit(6))

        !scratch(1:3*(noa+nob)) = 0.d0
        !psi(1:(noa+nob)**2) = dcmplx(0.d0,0.d0)
        !Priv%normV = 0.d0
        !Zion_coeff = Zion_coeff*dt
        !do itime = 1,nstep
        !  ii = (itime-1)*(noa+nob)**2
        !  write(iout,"('coeff2',16f10.6)") 1.d+2*Zion_coeff(ii+1+noa+nob:ii+2*(noa+nob))
        !  psi(1:(noa+nob)**2) = psi(1:(noa+nob)**2) + Zion_coeff(ii+1:ii+(noa+nob)**2)
        !  psi(1+noa+nob:2*(noa+nob)) = psi(1+noa+nob:2*(noa+nob)) &
        !      + Zion_coeff(ii+1+noa+nob:ii+2*(noa+nob))
        !  write(iout,"('psi2  ',16f10.6)") 1.d+2*psi(1+noa+nob:2*(noa+nob))
        !  scratch(1:2*(noa+nob)) = 0.d0
        !  do j = 1,noa+nob
        !    jj = ii + (j-1)*(noa+nob) 
        !    call get_norm( scratch(j),noa+nob,Zion_coeff(jj+1:jj+(noa+nob)) )
        !    call get_norm( scratch(j+noa+nob),noa+nob,psi((j-1)*(noa+nob)+1:j*(noa+nob)) )
        !    scratch(j+2*(noa+nob)) = scratch(j+2*(noa+nob))+ scratch(j)
        !  end do
        !  Priv%norm = 0.d0
        !  do j = 1,noa+nob
        !    Priv%norm = Priv%norm + scratch(j)**2
        !  end do
        !  Priv%norm = sqrt(Priv%norm)
        !  Priv%normV = Priv%normV + Priv%norm
        !  write(iout,"(i5,9f12.8)")itime,Priv%norm,(scratch(j),j=1,noa+nob)
        !  write(iout,"(' ion ',9f12.8)") Priv%normV,(scratch(j+2*(noa+nob)),j=1,noa+nob)
        !  Priv%norm = 0.d0
        !  do j = 1,noa+nob
        !    Priv%norm = Priv%norm + scratch(j+noa+nob)**2
        !  end do
        !  Priv%norm = sqrt(Priv%norm)
        !  write(iout,"(' psi ',9f12.8)") Priv%norm,(scratch(j+(noa+nob)),j=1,noa+nob)
        !end do

        !$OMP CRITICAL                    
        !: record data at last timestep
        tdciresults( idir + (iemax-1)*ndir)%norm0 = Priv%norm**2
        tdciresults( idir + (iemax-1)*ndir)%dipx  = Priv%mux
        tdciresults( idir + (iemax-1)*ndir)%dipy  = Priv%muy
        tdciresults( idir + (iemax-1)*ndir)%dipz  = Priv%muz
        
        write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
        write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           Priv%dirx1, Priv%diry1, Priv%dirz1, Priv%emax1
        write(iout,"(12x,'propagation time:',f12.4,' s')") Priv%finish2 - Priv%start2  
        write(iout,"(12x,'final norm = ',f10.5)")          Priv%norm**2
        !: debug times
        do i=1,7
          write(iout, "('time ',i5,': ',f14.6,' s')") i, times(i) 
        end do


        write(iout,"(12x,'For input NOs, nva95max=',i5,', nva99max=',i5)") Priv%nva95max, Priv%nva99max
        write(iout,"('For MOs (direct): rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
              Priv%rate_direct, Priv%nva95max_direct, Priv%nva99max_direct
        write(iout,"('For MOs (density): rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
              Priv%rate_debug, Priv%nva95max_debug, Priv%nva99max_debug
        !write(iout,"('DEBUG:   rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
        !      rate_debug, nva95max_debug, nva99max_debug

        !: Put current nva maxs into overall nva max
        if (Prop%nva95maxmaxMO .lt. Priv%nva95maxMO) Prop%nva95maxmaxMO = Priv%nva95maxMO
        if (Prop%nva99maxmaxMO .lt. Priv%nva99maxMO) then
          Prop%nva99maxmaxMO = Priv%nva99maxMO
          !ratemax_density = rate_density
          Prop%ratemax_density = Priv%rate_debug
        end if
        if (Prop%nva95maxmax_direct .lt. Priv%nva95max_direct) Prop%nva95maxmax_direct = Priv%nva95max_direct
        if (Prop%nva99maxmax_direct .lt. Priv%nva99max_direct) then
          Prop%nva99maxmax_direct = Priv%nva99max_direct
          Prop%ratemax_direct = Priv%rate_direct
        end if
        !: Put current nva maxs into overall nva max
        if (Prop%nva95maxmax .lt. Priv%nva95max) Prop%nva95maxmax = Priv%nva95max
        if (Prop%nva99maxmax .lt. Priv%nva99max) Prop%nva99maxmax = Priv%nva99max

        if (Prop%nva95maxMO_sort .lt. Priv%nva95MO_sort) Prop%nva95maxMO_sort = Priv%nva95MO_sort
        if (Prop%nva99maxMO_sort .lt. Priv%nva99MO_sort) Prop%nva99maxMO_sort = Priv%nva99MO_sort
        if (Prop%nva95maxNO_sort .lt. Priv%nva95NO_sort) Prop%nva95maxNO_sort = Priv%nva95NO_sort
        if (Prop%nva99maxNO_sort .lt. Priv%nva99NO_sort) Prop%nva99maxNO_sort = Priv%nva99NO_sort

        flush(iout)

        !$OMP END CRITICAL
     end do emax_loop

    call Priv%deconstruct()
    deallocate(Priv)
    !write(iout, '(A, I5, L1)') "Priv Deallocated at end of thread ", ithread , allocated(Priv)

  end do dir_loop

  !$OMP END DO
  ithread = omp_get_thread_num()
  !$OMP END PARALLEL

  !: We use Priv variables a couple times as placeholders so lets reallocate it
  allocate(Priv)
  call Priv%initialize()

  write(iout,*) "POST-PROPAGATION ANALYSIS:"
  flush(iout)

  if (flag_ReadU_NO) then
    write(iout,"('For input NOs, all directions (UNSORTED): nva95max=',i5,', nva99max=',i5)") Prop%nva95maxmax, Prop%nva99maxmax
    write(iout,"('For input NOs, all directions (  SORTED): nva95max=',i5,', nva99max=',i5)") Prop%nva95maxNO_sort, Prop%nva99maxNO_sort

    write(iout,"('For MOs, all directions: rate(density,   SORTED)=',F10.7,', nva95max=',i5, &
                ', nva99max=',i5)") Prop%ratemax_density, Prop%nva95maxMO_sort, Prop%nva99maxMO_sort
    write(iout,"('For MOs, all directions: rate(density, UNSORTED)=',F10.7,', nva95max=',i5, &
                ', nva99max=',i5)") Prop%ratemax_density, Prop%nva95maxmaxMO, Prop%nva99maxmaxMO
    write(iout,"('For MOs, all directions: rate(direct, UNSORTED)=',F10.7,', nva95max=',i5, & 
                ', nva99max=',i5)") Prop%ratemax_direct, Prop%nva95maxmax_direct, Prop%nva99maxmax_direct
    Priv%nva95max = 0
    Priv%nva99max = 0
    Priv%nva95maxMO = 0
    Priv%nva99maxMO = 0
    Priv%nva95max_direct = 0
    Priv%nva99max_direct = 0
    Priv%nva95max_debug = 0
    Priv%nva99max_debug = 0
    Priv%nva95MO_sort = 0
    Priv%nva99MO_sort = 0
    Priv%nva95NO_sort = 0
    Priv%nva99NO_sort = 0
    call update_maxnva( Priv%nva95maxMO, Priv%nva99maxMO, Priv%nva95max, Priv%nva99max, &
                        Priv%nva95MO_sort, Priv%nva99MO_sort, Priv%nva95MO_sort,Priv%nva99NO_sort, &
                        Priv%nva95max_debug, Priv%nva99max_debug, Priv%rate_debug, &
                        !Priv%nva95max_density, Priv%nva99max_density, Priv%rate_density, &
                        Prop%opdm_avg, Prop%U_NO_input, Mol%vabsmoa, .True. )
  end if

  call generate_natural_orbitals( Prop%opdm_avg, Prop%U_NO, Prop%natorb_occ, .false. )
  call generate_natural_orbitals( Prop%opdm_avg_abs, Prop%U_NO_abs, Prop%natorb_occ_abs, .false. )
  call NO_rate_sanity2( Prop%opdm_avg, Prop%U_NO, Prop%natorb_occ, Prop%opdm_avg_abs, Prop%U_NO_abs, Prop%natorb_occ_abs, Mol%vabsmoa, Mol%cmo_a)

  !call io_bin_test
  call write_dbin(Prop%U_NO, (noa+nva)*(noa+nva), "U_NO_out.bin")

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

  !: shared variables
  class(PropagationShared), allocatable :: Prop
  integer(8) :: nstuse2 
  complex(8) :: exphel(nstuse), psi0(nstuse)
  real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse) 
  real(8) :: hp2(nstuse*nstuse), tdvals2(nstuse)        
  

  !: private variables
  class(PropagationPrivate), allocatable :: Priv    
  !: Psi arrays must be not allocatable for performance
  complex(8) :: psi_j, psi_k, psi(nstuse), psi1(nstates)

  real(8) :: norm0
  real(8),allocatable :: pop0(:),pop1(:),ion(:)
  complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)

  integer(8) :: i, j, ii, jj, k, kk
  integer(8) :: itime, idir, iemax, ndim
  integer(8) :: ithread, idata
  complex(8) :: cdum

  real(8) :: start1, finish1, start2, finish2, start3, finish3
  real(8) :: times(100)

  !: lapack stuff
  integer(8) :: info1, info2, lscratch, liwork
  integer(8), allocatable :: iwork(:)
  real(8), allocatable    :: scratch(:)

  call write_header( 'trotter_circular','propagate','enter' )    
  write(iout, *) "Before writeme_propagate equation" ; flush(iout)
  call writeme_propagate( 'trot_cir', 'equation' ) 
  write(iout, *) "After writeme_propagate equation" ; flush(iout)
  call cpu_time(start)

  !: nstuse
  nstuse2 = nstuse*nstuse

  !: initialize psi0
  psi0 = dcmplx(0.d0,0.d0)
  do i=1, init_states(0)
    psi0( init_states(i) ) = init_coeffs(i)
  end do

  !: normalize
  call get_norm( norm0, nstuse, psi0 )
  psi0 = psi0 / norm0 
  
  !: write psi0
  write(iout, *) "Before writeme_propagate psi0" ; flush(iout)
  call writeme_propagate( 'trot_lin', 'psi0' )    
  write(iout, *) "After writeme_propagate psi0" ; flush(iout)
  

  ndim = noa+nva
  !: Initialize shared variables
  allocate(Prop)
  call Prop%initialize(ndim)

  !: get initial population
  allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
  i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
  allocate(ion_coeff(i))
  If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
    allocate(Zion_coeff(nstep*ip_states*(ip_states+1)))
  else
    allocate(Zion_coeff(ip_states*(ip_states+1)))
    allocate(Zproj_ion(ip_states*ip_states))
  end If
  
  norm0 = 1.d0
  
  call get_psid( nstuse, nstates, cis_vec, norm0, psi0, psi_det0 )
  select case ( trim(jobtype) ) 
  case( flag_cis ) 
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
  end select
  write(iout, *) "Before writeme_propagate pop0" ; flush(iout)
  call writeme_propagate( 'trot_lin','pop0', pop0, norb )
  write(iout, *) "After writeme_propagate pop0" ; flush(iout)

  deallocate( pop0 )
 
  allocate(Priv)
  call Priv%initialize()

  If( QeigenDC ) then
    allocate( iwork(3+5*nstuse) )
    allocate( scratch(1+8*nstuse+2*nstuse*nstuse) )
  else
    allocate( iwork(2) )
    allocate( scratch(nstuse*nstuse) )
  end if

  !: exphel = exp(-iH*dt/2)
  do i=1, nstuse
    Priv%temp = -0.5d0 * cis_eig(i) * dt
    exphel(i) = dcmplx( dcos(Priv%temp) , dsin(Priv%temp) )
  end do


  !: Start loop over directions

  !$OMP PARALLEL DEFAULT(NONE),&
  !$OMP PRIVATE(Priv, i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
  !$OMP ithread, cdum, &
  !$OMP norm0, lscratch, liwork, scratch, iwork, info1, info2, start1, finish1, &
  !$OMP start2, finish2, start3, finish3, &
  !$OMP psi_j, psi_k, psi, psi1, times, &
  !$OMP pop1, ion, ion_coeff, psi_det0,  &
  !$OMP hp1, hp2, tdvals1, tdvals2, Zion_coeff ),  &
  !$OMP SHARED( Mol, Prop, jobtype, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
  !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
  !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
  !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
  !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
  !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
  !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
  !$OMP flag_ReadU_NO )

  
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

     !: Form the transition dipole in (dirx1,diry1,dirz1) and (dirx2,diry2,dirz2) directions,
     hp1(:) = Priv%dirx1*tdx(1:nstuse2) + Priv%diry1*tdy(1:nstuse2) + Priv%dirz1*tdz(1:nstuse2)
     hp2(:) = Priv%dirx2*tdx(1:nstuse2) + Priv%diry2*tdy(1:nstuse2) + Priv%dirz2*tdz(1:nstuse2)

     !: diagonalize hp1 and hp2
     info1=10  
     info2=10  
     If( QeigenDC ) then
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

        call PropWriteDataHeaders(Priv, iemax, idir, tdciresults, psi0, Zion_coeff, 1)

        !$OMP CRITICAL
        write(iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0 )") idir, iemax, ithread
        flush(iout)
        !$OMP END CRITICAL
        
        !: initialize psi
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
        If( Qread_ion_coeff ) psi = dcmplx(0.d0,0.d0)

        !: begin Looping over time
        call cpu_time( start2 )
        timestep_loop : do itime=1, nstep-1
           
           call cpu_time(start3)

           Priv%efield1 = 0.5d0 * Priv%emax1*(fvect1(itime)+fvect1(itime+1)) !: modified midpoint
           Priv%efield2 = 0.5d0 * Priv%emax2*(fvect2(itime)+fvect2(itime+1)) !: modified midpoint
           if( read_shift(iemax).ne.0 ) then
             If(itime.eq.1) write(iout,"(' Pulse has been shifted by ',i6,' steps')") read_shift(iemax)
             i = itime-read_shift(iemax)
             if( i.gt.0 .and. i.lt.nstep ) then
               Priv%efield1 = 0.5d0 * Priv%emax1 * ( fvect1(i) + fvect1(i+1) )
               Priv%efield2 = 0.5d0 * Priv%emax2 * ( fvect2(i) + fvect2(i+1) )
             else
               Priv%efield1 = 0.d0
               Priv%efield2 = 0.d0
             end if
           end if

           Priv%efieldx = Priv%dirx1 * Priv%efield1 + Priv%dirx2 * Priv%efield2
           Priv%efieldy = Priv%diry1 * Priv%efield1 + Priv%diry2 * Priv%efield2
           Priv%efieldz = Priv%dirz1 * Priv%efield1 + Priv%dirz2 * Priv%efield2
           Priv%temp1 = dt * Priv%efield1 * 0.5d0 
           Priv%temp2 = dt * Priv%efield2
           
           If( itime.lt.ion_sample_start(iemax) ) psi = dcmplx(0.d0,0.d0)
           If( Qread_ion_coeff ) then
             If( itime.ge.ion_sample_start(iemax) .and. &
               itime.le.ion_sample_start(iemax)+ion_sample_width(iemax) ) then
               ii = (itime-1)*(noa+nob)*(noa+nob+1)
               !: tranform ion_coeff to CI basis and add to psi
!:                  write(iout,"('R0',i5,i3,16F10.6)") itime,ion_sample_state(iemax)
               call get_ion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                    dt,cis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
             end If
           else
             If( itime.eq.ion_sample_start(iemax) ) psi = psi0
           end If
           
           !: exp(-iHel dt/2 ) * psi
           do j = 1, nstuse
              psi(j) = exphel(j) * psi(j)
           end do
           
           !: W1 * exp(-Vabs dt/2) * psi
           psi1 = dcmplx( 0.d0, 0.d0 )
           do j = 1, nstuse
              jj = (j-1) * nstuse
              psi_j = psi(j)
              do k = 1 , nstuse
                 psi1(k) = psi1(k) + hp1(jj+k) * psi_j
              end do
           end do

           !: exp( iE1(t+dt/2)*mu * dt/2) * psi
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
           
           If( Qwrite_ion_coeff ) then
              call get_norm( Priv%norm, nstuse, psi )
              call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )
              ii = (itime-1)*(noa+nob)*(noa+nob+1)
              call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                Priv%rate, Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,scratch)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,Priv%rate,scratch(1:noa+nob)
                do i = 1,noa+nob
                   do j =1,noa+nob
                     Zion_coeff(ii+i*(noa+nob)+j) =  &
                       dsqrt(Priv%rate) * scratch(i) * Zion_coeff(ii+i*(noa+nob)+j)
                   end do
                end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
           end If
           
           analysis : if ( mod(itime,outstep).eq.0 ) then

              idata = int(itime/outstep)

              call get_norm( Priv%norm, nstuse, psi )
              call get_psid( nstuse, nstates, cis_vec, Priv%norm, psi, psi_det0 )
              if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                call pop_rate_ip(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                  hole_index,part_index,state_ip_index,ip_states, &
                  pop1,ion,ion_coeff,Priv%rate_aa,Priv%rate_ab,Priv%rate_ba,Priv%rate_bb, &
                  psi_det0,psi1,Priv%normV,Mol%vabsmoa,Mol%vabsmob,scratch,au2fs,Priv%rate_direct)
              else
                call pop_rate(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                  hole_index,part_index,state_ip_index,ip_states, &
                  pop1,ion,ion_coeff,Priv%rate_a,Priv%rate_b,psi_det0,psi1,Priv%normV, &
                  Mol%vabsmoa,Mol%vabsmob,unrestricted,scratch,au2fs,Priv%rate_direct)
              end if
              !: 1-RDM should be stored in scratch now.
              !: Why does scratch have size (nstuse,nstuse) while density in
              !:  pop_rate has size (norb,norb)?

              !: Check max nva for input NOs at this step.
              Priv%nva95max_direct = jj
              Priv%nva99max_direct = kk
              Prop%ratemax_direct = Priv%rate_direct
              call update_maxnva( Priv%nva95maxMO, Priv%nva99maxMO, Priv%nva95max, Priv%nva99max, &
                      Prop%nva95maxMO_sort, Prop%nva99maxMO_sort, Prop%nva95maxMO_sort,Prop%nva99maxMO_sort, &
                      Priv%nva95max_debug, Priv%nva99max_debug,Priv%rate_debug, &
                      !Priv%nva95max_density, Priv%nva99max_density, Priv%rate_density, &
                      Prop%opdm_avg, Prop%U_NO_input, Mol%vabsmoa )

              !$OMP CRITICAL
              !: set critical for sum so we dont have multiple threads
              !:   modifying opdm_avg at a time... is this going to destroy
              !:   performance?
              !: Add 1-RDM (opdm) to 1-RDM average for Natural Orbital generation.
              
              call add_opdm_average( Prop%opdm_avg, scratch, Prop%opdm_avg_N )
              call add_opdm_average( Prop%opdm_avg_abs, scratch, Prop%opdm_avg_N, .false., .true. )

              !call generate_natural_orbitals( scratch, U_NO, natorb_occ )
              !call NO_rate_sanity( scratch, U_NO, natorb_occ, Mol%vabsmoa, Mol%cmo_a)
              !flush(iout)    
              !call dgemm_sanity 

              !$OMP END CRITICAL

              ion_coeff = dcmplx( 0.d0,0.d0 )
              do i=1,ip_states
                cdum = dcmplx( 0.d0,0.d0 )
                do j=1,nstuse
!:                    cdum = cdum + proj_ion(j+(i-1)*nstuse)*psi(j)
                end do
!:                write(iout,"(' ion_coeff_det',I4,4f12.7)") i,cdum/Priv%norm,abs(cdum/Priv%norm)
                do k=1,ip_states
                  ion_coeff(k) = ion_coeff(k) + ip_vec(i+(k-1)*ip_states)*cdum
                end do
              end do

              call get_norm( Priv%norm,nstuse, psi )                
              call get_expectation( nstuse, Priv%norm, psi, abp, Priv%rate) !: rate expectation value
              call get_expectation( nstuse, Priv%norm, psi, tdx, Priv%mux ) !: mux  expectation value
              call get_expectation( nstuse, Priv%norm, psi, tdy, Priv%muy ) !: muy  expectation value
              call get_expectation( nstuse, Priv%norm, psi, tdz, Priv%muz ) !: muz  expectation value
!:                write(iout,*) " expectation rate", itime,Priv%rate
              Priv%rate = -2.d0 * Priv%rate * Priv%norm**2

              if( Qci_save ) then
                write(Priv%funit(1)) dble(itime)*dt*au2fs, Priv%efield1, Priv%efield2, Priv%dirx1, Priv%diry1, Priv%dirz1, &
                  Priv%dirx2, Priv%diry2, Priv%dirz2, Priv%norm**2
                write(Priv%funit(1)) real( psi )
                write(Priv%funit(1)) aimag( psi )
                flush(Priv%funit(1))
              end if
              
              if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                  if(Priv%norm.ne.0) then
                    write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs,Priv%nva99max,kk,Priv%efield1,0.d0,Priv%efieldx,Priv%efieldy,Priv%efieldz, &
                    Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz, &
                    dble(dconjg(psi(1))*psi(1))/Priv%norm**2,dble(dconjg(psi(2))*psi(2))/Priv%norm**2, &
                    dble(dconjg(psi(3))*psi(3))/Priv%norm**2,dble(dconjg(psi(4))*psi(4))/Priv%norm**2, &
                    dble(dconjg(psi(5))*psi(5))/Priv%norm**2,dble(dconjg(psi(6))*psi(6))/Priv%norm**2, &
                    dble(dconjg(psi(7))*psi(7))/Priv%norm**2,dble(dconjg(psi(8))*psi(8))/Priv%norm**2
                 else
                    write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs,Priv%nva99max,kk,Priv%efield1,0.d0,Priv%efieldx,Priv%efieldy,Priv%efieldz, &
                    Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
                 end if
              else
                  write( Priv%funit(2),"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs,Priv%nva99max,kk,Priv%efield1,Priv%efield2,Priv%efieldx,Priv%efieldy,Priv%efieldz, &
                    Priv%norm**2, Priv%rate/au2fs, Priv%mux, Priv%muy, Priv%muz
              end if
              flush(Priv%funit(2))

              if( trim(jobtype).eq.flag_ip) then
                  write( Priv%funit(3),"( i5,f10.4,500(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, Priv%norm**2, &
                    (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),Priv%rate/au2fs,(Priv%rate_aa(i),i=1,noa*noa),&
                    (Priv%rate_ab(i),i=1,noa*nob),(Priv%rate_ba(i),i=1,nob*noa),(Priv%rate_bb(i),i=1,nob*nob)
              else
                  write( Priv%funit(3),"( i5,f10.4,500(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, Priv%norm**2, &
                    (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob), &
                    (Priv%rate_a(i),i=1,noa),(Priv%rate_b(i),i=1,nob)
              end if
              flush(Priv%funit(3))

              write( Priv%funit(4),"( i5,f10.4,500(1x,f15.10))") &
                idata, dble(itime)*dt*au2fs, Priv%rate/Priv%norm**2, Priv%normV, &
                (ion(i),i=1,noa),(ion(noa+nva+i),i=1,nob),(ion_coeff(i),i=1,ip_states)
              flush(Priv%funit(4))

              if( Qmo_dens) then
              !: write real part of total density
                write(Priv%funit(6),"(f13.9,f16.10)") dble(itime)*dt*au2fs,Priv%rate
                do i=1, noa+nva
                  if ( abs(scratch(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(Priv%funit(6),"(i5,i5,1x,f13.10,',')",advance='no') &
                    write(Priv%funit(6),"(i5,i5,1x,f13.10)") &
                      i, i, scratch(i+(i-1)*(noa+nva))
                  do j=(i+1), noa+nva
                    if ( abs(scratch(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(Priv%funit(6),"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(Priv%funit(6),"(i5,i5,1x,f13.10)") &
                        i, j, 2.d0*scratch(i+(j-1)*(noa+nva))
                    end do
                  end do
                  write(Priv%funit(6),"(i5,i5,1x,f13.10)") 0,0,0.d0
              end if

              If(trim(jobtype).eq.flag_cis .and. (.not. Qwrite_ion_coeff) ) then
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                  Priv%rate,Zion_coeff,ion_coeff,scratch)
                call get_proj_ion(iout,noa+nob,ip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(ip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                write( Priv%funit(5),"(i5,f10.4,50(1x,f15.10))") &
                  idata, dble(itime)*dt*au2fs, Priv%rate,(scratch(i),i=1,noa+nob)
                do i = 1,noa+nob
                  write( Priv%funit(5),"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                end do
                do i = 1,noa+nob
                  write( Priv%funit(5),"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                end do
              end If
              If(trim(jobtype).eq.flag_ip .and. (.not. Qwrite_ion_coeff) ) then
                call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  ip_states,state_ip_index,psi_det0,psi1,Priv%norm,Mol%vabsmoa,Mol%vabsmob, &
                  Priv%rate,Zion_coeff,ion_coeff,scratch)
                call get_proj_ion(iout,ip_states,ip_vec,Zion_coeff(ip_states+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  write(iout,"('s    ',14f13.7/5x,14f13.7/5x,14f13.7/5x,14f13.7)") (scratch(j),j=1,ip_states)
!:                  write(iout,"('rates',14f13.7/5x,14f13.7/5x,14f13.7/5x,14f13.7)") (abs(Zion_coeff(j)),j=1,ip_states)
!:                  flush(iout)
!:                  do i=1,ip_states
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
!:                    write(iout,"('coeff',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") i,(Zion_coeff(j+i*(ip_states)),j=1,ip_states)
!:                    write(iout,"('ipvec',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") i,(ip_vec(j+(i-1)*(ip_states)),j=1,ip_states)
!:                    write(iout,"('proj ',i3,16f13.7/8x,16f13.7/8x,16f13.7/8x,16f13.7)") i,(Zproj_ion(j+(i-1)*(ip_states)),j=1,ip_states)
!:                  end do
!:                  flush(iout)
                write( Priv%funit(5),"(i5,f10.4,50(1x,f15.10))") &
                  idata, dble(itime)*dt*au2fs, Priv%rate,(scratch(i),i=1,ip_states)
                do i = 1,ip_states
                  write( Priv%funit(5),"(60(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                end do
                do i = 1,ip_states
                  write( Priv%funit(5),"(60(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                end do
                flush(Priv%funit(5))
              end If

           end if analysis

        end do timestep_loop
        call cpu_time(finish2)

        if( Qci_save ) close(Priv%funit(1))
        close(Priv%funit(2))
        close(Priv%funit(3))
        close(Priv%funit(4))
        If( Qwrite_ion_coeff ) write(Priv%funit(5)) Zion_coeff
        If( Qread_ion_coeff .or. Qwrite_ion_coeff ) close(Priv%funit(5))

        If( Qmo_dens ) close(Priv%funit(6))

        !$OMP CRITICAL
        !: record data at last timestep
        tdciresults( idir + (iemax-1)*ndir)%norm0 = Priv%norm**2
        tdciresults( idir + (iemax-1)*ndir)%dipx  = Priv%mux
        tdciresults( idir + (iemax-1)*ndir)%dipy  = Priv%muy
        tdciresults( idir + (iemax-1)*ndir)%dipz  = Priv%muz

        write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread,idir,iemax
        write(iout,"(12x,'dir1 = (',f8.5,',',f8.5,',',f8.5,')  emax1 = ',f8.5,' au')") Priv%dirx1,Priv%diry1,Priv%dirz1,Priv%emax1
        write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") Priv%dirx2,Priv%diry2,Priv%dirz2,Priv%emax2
        write(iout,"(12x,'(iemax=1) TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
        write(iout,"(12x,'(iemax=1) LAPACK dysev TD diagonalization INFO=',i0)") info1
        write(iout,"(12x,'propagation time:',f12.4,' s')") finish2-start2
        write(iout,"(12x,'final norm = ',f10.5)")          Priv%norm**2


        if (flag_ReadU_NO) then
          write(iout,"('For MOs (direct): rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
                  Priv%rate_direct, Priv%nva95max_direct, Priv%nva99max_direct
          write(iout,"('For MOs (density): rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
                  Priv%rate_debug, Priv%nva95max_debug, Priv%nva99max_debug
          !write(iout,"('DEBUG:   rate=',E24.16,' nva95max=',i5,', nva99max=',i5)") &
          !        Priv%rate_debug, nva95max_debug, nva99max_debug
          !: Put current nva maxs into overall nva max
          if (Prop%nva95maxmaxMO .lt. Priv%nva95maxMO) Prop%nva95maxmaxMO = Priv%nva95maxMO
          if (Prop%nva99maxmaxMO .lt. Priv%nva99maxMO) then
            Prop%nva99maxmaxMO = Priv%nva99maxMO
            !ratemax_density = Priv%rate_density
            Prop%ratemax_density = Priv%rate_debug
          end if
          if (Prop%nva95maxmax_direct .lt. Priv%nva95max_direct) Prop%nva95maxmax_direct = Priv%nva95max_direct
          if (Prop%nva99maxmax_direct .lt. Priv%nva99max_direct) then
            Prop%nva99maxmax_direct = Priv%nva99max_direct
            Prop%ratemax_direct = Priv%rate_direct
          end if

          write(iout,"(12x,'For input NOs, nva95max=',i5,', nva99max=',i5)") Priv%nva95max, Priv%nva99max
          !: Put current nva maxs into overall nva max
          if (Prop%nva95maxmax .lt. Priv%nva95max) Prop%nva95maxmax = Priv%nva95max
          if (Prop%nva99maxmax .lt. Priv%nva99max) Prop%nva99maxmax = Priv%nva99max

          if (Prop%nva95maxMO_sort .lt. Priv%nva95MO_sort) Prop%nva95maxMO_sort =Priv%nva95MO_sort
          if (Prop%nva99maxMO_sort .lt. Priv%nva99MO_sort) Prop%nva99maxMO_sort =Priv%nva99MO_sort
          if (Prop%nva95maxNO_sort .lt. Priv%nva95NO_sort) Prop%nva95maxNO_sort =Priv%nva95NO_sort
          if (Prop%nva99maxNO_sort .lt. Priv%nva99NO_sort) Prop%nva99maxNO_sort =Priv%nva99NO_sort

        end if
        flush(iout)

        !$OMP END CRITICAL

     end do emax_loop

    call Priv%deconstruct()
    deallocate(Priv)
    !write(iout, '(A, I5, L1)') "Priv Deallocated at end of thread ", ithread , allocated(Priv)

  end do dir_loop

  !$OMP END DO
  !$OMP END PARALLEL

  write(iout,*) "POST-PROPAGATION ANALYSIS:"
  flush(iout)

  if (flag_ReadU_NO) then
    write(iout,"('For input NOs, all directions (UNSORTED): nva95max=',i5,', nva99max=',i5)") Prop%nva95maxmax, Prop%nva99maxmax
    write(iout,"('For input NOs, all directions (  SORTED): nva95max=',i5,', nva99max=',i5)") Prop%nva95maxmax, Prop%nva99maxmax

    write(iout,"('For MOs, all directions: rate(density,   SORTED)=',F10.7,', nva95max=',i5, &
                ', nva99max=',i5)") Prop%ratemax_density, Prop%nva95maxMO_sort, Prop%nva99maxmaxMO
    write(iout,"('For MOs, all directions: rate(density, UNSORTED)=',F10.7,', nva95max=',i5, &
                ', nva99max=',i5)") Prop%ratemax_density, Prop%nva95maxmaxMO, Prop%nva99maxmaxMO
    write(iout,"('For MOs, all directions: rate(direct, UNSORTED)=',F10.7,', nva95max=',i5, & 
                ', nva99max=',i5)") Prop%ratemax_direct, Prop%nva95maxmax_direct, Prop%nva99maxmax_direct
    Priv%nva95max = 0
    Priv%nva99max = 0
    Priv%nva95maxMO = 0
    Priv%nva99maxMO = 0
    Priv%nva95max_direct = 0
    Priv%nva99max_direct = 0
    Priv%nva95max_debug = 0
    Priv%nva99max_debug = 0
    Prop%nva95maxMO_sort = 0
    Prop%nva99maxMO_sort = 0
    Prop%nva95maxNO_sort = 0
    Prop%nva99maxNO_sort = 0
    call update_maxnva( Priv%nva95maxMO, Priv%nva99maxMO, Priv%nva95max, Priv%nva99max, &
                        Prop%nva95maxMO_sort, Prop%nva99maxMO_sort, Prop%nva95maxMO_sort, Prop%nva99maxNO_sort, &
                        Priv%nva95max_debug, Priv%nva99max_debug, Priv%rate_debug, &
                        Prop%opdm_avg, Prop%U_NO_input, Mol%vabsmoa, .True. )
  end if

  call generate_natural_orbitals( Prop%opdm_avg, Prop%U_NO, Prop%natorb_occ, .false. )
  call generate_natural_orbitals( Prop%opdm_avg_abs, Prop%U_NO_abs, Prop%natorb_occ_abs, .false. )
  call NO_rate_sanity2( Prop%opdm_avg, Prop%U_NO, Prop%natorb_occ, Prop%opdm_avg_abs, Prop%U_NO_abs, &
                        Prop%natorb_occ_abs, Mol%vabsmoa, Mol%cmo_a)

  !call io_bin_test
  call write_dbin(Prop%U_NO, (noa+nva)*(noa+nva), "U_NO_out.bin")

  call write_header( 'trotter_circular','propagate','leave' )
  call cpu_time(finish)
  write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start  


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


subroutine NO_rate_sanity( opdm, U_NO, NO_occ, Vabs, C_AO )
  implicit none

  real(8), intent(in) :: opdm(norb*norb) !: One-particle density matrix in MO basis
  real(8), intent(in) :: U_NO(norb*norb) !: Transformation matrix from MO -> NO.
  real(8), intent(in) :: NO_occ(norb)    !: NO occupations (eigenvalues)
  real(8), intent(in) :: Vabs(norb*norb) !: Expectation value <n|Vabs|m> in MO basis
  real(8), intent(in) :: C_AO(norb*norb) !: AO -> MO Coefficients

  real(8) :: rate_NO, rate_MO, rate_NO_abs, rate_MO_abs, rate_NO2 
  integer :: nva95_NO, nva99_NO, nva95_MO, nva99_MO
  real(8) :: tmp_NO, tmp_MO
 
  integer :: i, j, k, ndim
  real(8) :: temp(norb*norb), temp2(norb*norb) !: Temporary matrix for multiplication
  real(8) :: tmprate_NO(norb), tmprate_MO(norb)
  real(8) :: opdm_NO(norb*norb) !: opdm in NO basis
  real(8) :: Vabs_NO(norb*norb) !: Expectation value <n|Vabs|m> in NO basis
  logical :: verbose

  temp = 0.d0
  opdm_NO = 0.d0
  Vabs_NO = 0.d0
  tmprate_NO = 0.d0
  tmprate_MO = 0.d0
  verbose = .false.
  !: LATER ADD CASES FOR UHF
  ndim = noa+nvb

  !: Write U_NO
  if (verbose) then
    write(iout, '(A)') 'U_NO:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') U_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Write Vabs (MO)
  if (verbose) then
    write(iout, '(A)') 'Vabs (MO):'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') Vabs((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Transform Vabs to NO basis
  !: U_NO^T * Vabs * U_NO = Vabs_NO
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, Vabs, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, Vabs_NO, ndim)

  !: Transform opdm to NO basis
  !: U_NO^T * opdm * U_NO = opdm_NO
  temp = 0.d0
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, opdm, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, opdm_NO, ndim)


  !: rate = 2*Tr(opdm_NO * Vabs_NO)
  !: temp = opdm_NO * Vabs_NO
  temp = 0.d0
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, opdm_NO, ndim, Vabs_NO, ndim, 0.d0, temp, ndim)
  !: temp2 = opdm_MO * Vabs_MO
  temp2 = 0.d0
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, opdm, ndim, Vabs, ndim, 0.d0, temp2, ndim)


  !: Write opdm_NO
  if (verbose) then
    write(iout, '(A)') 'opdm_NO:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') opdm_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Calculate rate (Matrix multiply method)
  rate_NO = 0.d0
  rate_MO = 0.d0
  rate_NO_abs = 0.d0
  rate_MO_abs = 0.d0
  do i = 1, ndim
    rate_NO = rate_NO + 2.d0*(temp((i-1)*ndim+i))  !: temp[i][i]
    rate_MO = rate_MO + 2.d0*(temp2((i-1)*ndim+i)) !: temp2[i][i]
    rate_NO_abs = rate_NO_abs + 2.d0*abs(temp((i-1)*ndim+i))  !: temp[i][i]
    rate_MO_abs = rate_MO_abs + 2.d0*abs(temp2((i-1)*ndim+i)) !: temp2[i][i]
  end do

  !: Calculate rate as dot(occ,Vabs_NO) since Vabs_NO is diagonal
  rate_NO2 = 0.d0
  do i=1,ndim
    rate_NO2 = rate_NO2 + 2.d0*NO_occ(i)*Vabs_NO((i-1)*ndim+i)
  end do

  !write(iout, '(A)') 'Natural Orbital Rate Analysis (Matrix Multiply):'
  !write(iout, '(A)') '  Orb#  Occ            Vabs_NO        Rate'
  !write(iout, '(A)') '=========================================='
  !do i = ndim,1, -1
  !  write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ(i), &
  !   ",  ", Vabs_NO((i-1)*ndim+i), ",  ", abs(2.d0*temp((i-1)*ndim+i))
  !end do

  write(iout, '(A, E24.16)') '2*Tr(    opdm_NO * Vabs_NO)):   ', rate_NO
  write(iout, '(A, E24.16)') '2*Tr(    opdm_MO * Vabs_MO)):   ', rate_MO
  write(iout, '(A, E24.16)') '2*Tr(abs(opdm_NO * Vabs_NO)):   ', rate_NO_abs
  write(iout, '(A, E24.16)') '2*Tr(abs(opdm_MO * Vabs_MO)):   ', rate_MO_abs
  write(iout, '(A, E24.16)') 'dot(NO_occ, diag(Vabs_NO)   : ', rate_NO2

  !: Calculate rate for Matrix Multiply rate = 2*Tr(abs(opdm*Vabs))
  tmp_NO = 0.d0
  tmp_MO = 0.d0
  nva95_NO = 0
  nva99_NO = 0
  nva95_MO = 0
  nva99_MO = 0
  ! U_NO should be sorted by descending evals
  !do i=ndim,1, -1
  do i=1,ndim
    !write(iout, '(I5, A, F12.10)' ) i, ' NO: ', 2*temp((i-1)*ndim+i)
    tmp_NO = tmp_NO + abs(2.d0*temp((i-1)*ndim+i))
    if (tmp_NO/rate_NO .lt. 0.95d0) nva95_NO = nva95_NO + 1
    if (tmp_NO/rate_NO .lt. 0.99d0) nva99_NO = nva99_NO + 1
  end do
  do i=1,ndim
    !write(iout, '(I5, A, F12.10)' ) i, ' MO: ', 2*temp2((i-1)*ndim+i)
    tmp_MO = tmp_MO + abs(2.d0*temp2((i-1)*ndim+i))
    if (tmp_MO/rate_MO .lt. 0.95d0) nva95_MO = i
    if (tmp_MO/rate_MO .lt. 0.99d0) nva99_MO = i
  end do

  write(iout, '(A)') "Rate with matrix multiply rate = 2*Tr(opdm*Vabs):"
  !write(iout, '(A, I5, I5, I5, I5, I5)') "nva, nva95_NO, nva99_NO, nva95_MO, nva99_MO: ", nvb, &
  !  nva95_NO, nva99_NO, nva95_MO, nva99_MO

  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",rate_NO,", nva=",nva, &
    ", nva95_NO= ",nva95_NO,", nva99_NO= ",nva99_NO

  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",rate_MO,", nva=",nva, &
    ", nva95_MO= ",nva95_MO,", nva99_MO= ",nva99_MO


  write(iout, '(A)') "Rate with rate = sum_{i,j} abs(opdm_{i,j} * Vabs_{i,j}):"

  !: Calculate rate component-wise way rate = sum_{i,j} abs(opdm_{i,j} * Vabs_{i,j})
  tmp_NO = 0.d0
  tmp_MO = 0.d0
  nva95_NO = 0
  nva99_NO = 0
  nva95_MO = 0
  nva99_MO = 0
  tmprate_MO(1) = 2*abs( opdm(1) * Vabs(1) )
  do i=2,ndim
    tmprate_MO(i) = tmprate_MO(i-1) + 2*abs(opdm((i-1)*ndim+i) * Vabs((i-1)*ndim+i) )
    do j=1,i-1
      tmprate_MO(i) = tmprate_MO(i) + 2*abs(opdm((i-1)*ndim+j) * Vabs((i-1)*ndim+j) )
      tmprate_MO(i) = tmprate_MO(i) + 2*abs(opdm((j-1)*ndim+i) * Vabs((j-1)*ndim+i) )
    end do
  end do
  do i=1,ndim
    if(tmprate_MO(i).lt.0.95d0*tmprate_MO(ndim)) nva95_MO = i
    if(tmprate_MO(i).lt.0.99d0*tmprate_MO(ndim)) nva99_MO = i
  end do
  
  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmprate_MO(ndim),", nva=",nva, &
    ", nva95_MO= ",nva95_MO,", nva99_MO= ",nva99_MO


  tmprate_NO(1) = 2*abs( opdm_NO(1) * Vabs_NO(1) )
  do i=2,ndim
    tmprate_NO(i) = tmprate_NO(i-1) + 2*abs(opdm_NO((i-1)*ndim+i) * Vabs_NO((i-1)*ndim+i) )
    do j=1,i-1
      tmprate_NO(i) = tmprate_NO(i) + 2*abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j) )
      tmprate_NO(i) = tmprate_NO(i) + 2*abs(opdm_NO((j-1)*ndim+i) * Vabs_NO((j-1)*ndim+i) )
    end do
  end do
  do i=1,ndim
    if(tmprate_NO(i).lt.0.95d0*tmprate_NO(ndim)) nva95_NO = i
    if(tmprate_NO(i).lt.0.99d0*tmprate_NO(ndim)) nva99_NO = i
  end do
  
  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmprate_NO(ndim),", nva=",nva, &
    ", nva95_NO= ",nva95_NO,", nva99_NO= ",nva99_NO


  !write(iout, '(A)') 'Natural Orbital Rate Analysis (Component-wise):'
  !write(iout, '(A)') '  Orb#  Occ          Vabs_NO(i,i)   Cumulative Rate'
  !write(iout, '(A)') '====================================================='
  !!do i = ndim,1, -1
  !do i = 1,ndim
  !  write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ(i), &
  !   ",  ", Vabs_NO((i-1)*ndim+i), ",  ", tmprate_NO(i)
  !end do

  !write(iout, '(A)') 'Molecular Orbital Rate Analysis (Component-wise):'
  !write(iout, '(A)') '  Orb#            Vabs_MO(i,i)   Cumulative Rate'
  !write(iout, '(A)') '================================================='
  !!do i = ndim,1, -1
  !do i = 1,ndim
  !  write(iout,'(I5, A, F12.10, A, F12.10)') i,  &
  !   ",  ", Vabs((i-1)*ndim+i), ",  ", tmprate_MO(i)
  !end do
  flush(iout)    

end subroutine NO_rate_sanity

!: Sorry to code replicate!
!: This one includes analysis for the abs, keeping the old one for each-step analysis.
subroutine NO_rate_sanity2( opdm, U_NO, NO_occ, opdm_abs, U_absNO, NO_occ_abs, Vabs, C_AO )
  implicit none

  real(8), intent(in) :: opdm(norb*norb) !: One-particle density matrix in MO basis
  real(8), intent(in) :: U_NO(norb*norb) !: Transformation matrix from MO -> NO.
  real(8), intent(in) :: NO_occ(norb)    !: NO occupations (eigenvalues)
  real(8), intent(in) :: opdm_abs(norb*norb) !: One-particle density matrix in MO basis
  real(8), intent(in) :: U_absNO(norb*norb) !: Transformation matrix from MO -> NO.
  real(8), intent(in) :: NO_occ_abs(norb)    !: NO occupations (eigenvalues)
  real(8), intent(in) :: Vabs(norb*norb) !: <n|Vabs|m> in MO basis
  real(8), intent(in) :: C_AO(norb*norb) !: AO -> MO Coefficients

  real(8) :: rate_NO, rate_MO, rate_NO_abs, rate_MO_abs, rate_NO2 
  real(8) :: rate_absNO, rate_absNO2
  integer :: nva95_NO, nva99_NO, nva95_MO, nva99_MO
  integer :: nva95_absNO, nva99_absNO
  real(8) :: tmp_NO, tmp_MO, tmp_absNO
 
  integer :: i, j, k, ndim
  real(8) :: temp(norb*norb), temp2(norb*norb), temp3(norb*norb) !: Temporary matrix for multiplication
  real(8) :: tmprate_NO(norb), tmprate_MO(norb), tmprate_absNO(norb)
  real(8) :: opdm_NO(norb*norb), opdm_absNO(norb*norb) !: opdm in NO basis
  real(8) :: Vabs_NO(norb*norb), Vabs_absNO(norb*norb) !: Expectation value <n|Vabs|m> in NO basis
  logical :: verbose

  temp = 0.d0
  opdm_NO = 0.d0
  Vabs_NO = 0.d0
  tmprate_NO = 0.d0
  tmprate_MO = 0.d0
  verbose = .false.
  !: LATER ADD CASES FOR UHF
  ndim = noa+nvb

  !: Write U_NO
  if (verbose) then
    write(iout, '(A)') 'U_NO:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') U_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Write Vabs (MO)
  if (verbose) then
    write(iout, '(A)') 'Vabs (MO):'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') Vabs((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Transform Vabs to NO basis
  !: U_NO^T * Vabs * U_NO = Vabs_NO
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, Vabs, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, Vabs_NO, ndim)

  !: Transform opdm to NO basis
  !: U_NO^T * opdm * U_NO = opdm_NO
  temp = 0.d0
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, opdm, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, opdm_NO, ndim)


  !: Transform Vabs to absNO basis
  !: U_absNO^T * Vabs * U_absNO = Vabs_absNO
  temp = 0.d0
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_absNO, ndim, Vabs, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_absNO, ndim, 0.d0, Vabs_absNO, ndim)

  !: Transform opdm_abs to absNO basis
  !: U_absNO^T * opdm_abs * U_absNO = opdm_absNO
  temp = 0.d0
  call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_absNO, ndim, opdm_abs, ndim, 0.d0, temp, ndim)
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_absNO, ndim, 0.d0, opdm_absNO, ndim)


  !: rate = 2*Tr(opdm_NO * Vabs_NO)
  !: temp = opdm_NO * Vabs_NO
  temp = 0.d0
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, opdm_NO, ndim, Vabs_NO, ndim, 0.d0, temp, ndim)
  !: temp2 = opdm_MO * Vabs_MO
  temp2 = 0.d0
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, opdm, ndim, Vabs, ndim, 0.d0, temp2, ndim)
  !: temp3 = opdm_absNO * Vabs_absNO
  temp3 = 0.d0
  call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, opdm_absNO, ndim, Vabs_absNO, ndim, 0.d0, temp3, ndim)


  !: Write opdm_NO
  if (verbose) then
    write(iout, '(A)') 'opdm_NO:'
    do i = 1, ndim
      do j = 1, ndim
        write(iout, *, advance='no') opdm_NO((i-1)*ndim+j)
      end do
      write(iout, *) ! New line after each row
    end do
  end if

  !: Calculate rate (Matrix multiply method)
  rate_NO = 0.d0
  rate_MO = 0.d0
  rate_NO_abs = 0.d0
  rate_MO_abs = 0.d0
  rate_absNO = 0.d0

  do i = 1, ndim
    rate_NO = rate_NO + 2.d0*(temp((i-1)*ndim+i))  !: temp[i][i]
    rate_MO = rate_MO + 2.d0*(temp2((i-1)*ndim+i)) !: temp2[i][i]
    rate_NO_abs = rate_NO_abs + 2.d0*abs(temp((i-1)*ndim+i))  !: temp[i][i]
    rate_MO_abs = rate_MO_abs + 2.d0*abs(temp2((i-1)*ndim+i)) !: temp2[i][i]
    rate_absNO = rate_absNO + 2.d0*abs(temp3((i-1)*ndim+i))
  end do

  !: Calculate rate as dot(occ,Vabs_NO) since Vabs_NO is diagonal
  rate_NO2 = 0.d0
  do i=1,ndim
    rate_NO2 = rate_NO2 + 2.d0*NO_occ(i)*Vabs_NO((i-1)*ndim+i)
  end do
  rate_absNO2 = 0.d0
  do i=1,ndim
    rate_absNO2 = rate_absNO2 + 2.d0*NO_occ_abs(i)*Vabs_absNO((i-1)*ndim+i)
  end do

  !write(iout, '(A)') 'Natural Orbital Rate Analysis (Matrix Multiply):'
  !write(iout, '(A)') '  Orb#  Occ            Vabs_NO        Rate'
  !write(iout, '(A)') '=========================================='
  !do i = ndim,1, -1
  !  write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ(i), &
  !   ",  ", Vabs_NO((i-1)*ndim+i), ",  ", abs(2.d0*temp((i-1)*ndim+i))
  !end do

  write(iout, '(A, E24.16)') '2*Tr(    opdm_NO * Vabs_NO) :   ', rate_NO
  write(iout, '(A, E24.16)') '2*Tr(    opdm_MO * Vabs_MO) :   ', rate_MO
  write(iout, '(A, E24.16)') '2*Tr(abs(opdm_NO * Vabs_NO)):   ', rate_NO_abs
  write(iout, '(A, E24.16)') '2*Tr(abs(opdm_MO * Vabs_MO)):   ', rate_MO_abs
  write(iout, '(A, E24.16)') '2*Tr(opdm_absNO * Vabs_absNO):   ', rate_absNO
  write(iout, '(A, E24.16)') 'dot(NO_occ, diag(Vabs_NO)   : ', rate_NO2
  write(iout, '(A, E24.16)') 'dot(NO_occ_abs, diag(Vabs_absNO)   : ', rate_absNO2

  !: Calculate rate for Matrix Multiply rate = 2*Tr(abs(opdm*Vabs))
  tmp_NO = 0.d0
  tmp_MO = 0.d0
  tmp_absNO = 0.d0
  nva95_NO = 0
  nva99_NO = 0
  nva95_MO = 0
  nva99_MO = 0
  nva95_absNO = 0
  nva99_absNO = 0
  ! U_NO should be sorted by descending evals
  !do i=ndim,1, -1
  do i=1,ndim
    !write(iout, '(I5, A, F12.10)' ) i, ' NO: ', 2*temp((i-1)*ndim+i)
    tmp_NO = tmp_NO + abs(2.d0*temp((i-1)*ndim+i))
    if (tmp_NO/rate_NO .lt. 0.95d0) nva95_NO = nva95_NO + 1
    if (tmp_NO/rate_NO .lt. 0.99d0) nva99_NO = nva99_NO + 1
  end do
  do i=1,ndim
    !write(iout, '(I5, A, F12.10)' ) i, ' MO: ', 2*temp2((i-1)*ndim+i)
    tmp_MO = tmp_MO + abs(2.d0*temp2((i-1)*ndim+i))
    if (tmp_MO/rate_MO .lt. 0.95d0) nva95_MO = i
    if (tmp_MO/rate_MO .lt. 0.99d0) nva99_MO = i
  end do
  do i=1,ndim
    !write(iout, '(I5, A, F12.10)' ) i, ' NO: ', 2*temp((i-1)*ndim+i)
    tmp_absNO = tmp_absNO + abs(2.d0*temp3((i-1)*ndim+i))
    if (tmp_absNO/rate_absNO .lt. 0.95d0) nva95_absNO = nva95_absNO + 1
    if (tmp_absNO/rate_absNO .lt. 0.99d0) nva99_absNO = nva99_absNO + 1
  end do

  write(iout, '(A)') "Rate with matrix multiply rate = 2*Tr(opdm*Vabs):"
  !write(iout, '(A, I5, I5, I5, I5, I5)') "nva, nva95_NO, nva99_NO, nva95_MO, nva99_MO: ", nvb, &
  !  nva95_NO, nva99_NO, nva95_MO, nva99_MO

  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",rate_NO,", nva=",nva, &
    ", nva95_NO= ",nva95_NO,", nva99_NO= ",nva99_NO

  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",rate_MO,", nva=",nva, &
    ", nva95_MO= ",nva95_MO,", nva99_MO= ",nva99_MO

  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",rate_absNO,", nva=",nva, &
    ", nva95_absNO= ",nva95_absNO,", nva99_absNO= ",nva99_absNO

  write(iout, '(A)') "Rate with rate = sum_{i,j} abs(opdm_{i,j} * Vabs_{i,j}):"

  !: Calculate rate component-wise way rate = sum_{i,j} abs(opdm_{i,j} * Vabs_{i,j})
  tmp_NO = 0.d0
  tmp_MO = 0.d0
  tmp_absNO = 0.d0
  nva95_NO = 0
  nva99_NO = 0
  nva95_MO = 0
  nva99_MO = 0
  nva95_absNO = 0
  nva99_absNO = 0
  tmprate_MO(1) = 2*abs( opdm(1) * Vabs(1) )
  do i=2,ndim
    tmprate_MO(i) = tmprate_MO(i-1) + 2*abs(opdm((i-1)*ndim+i) * Vabs((i-1)*ndim+i) )
    do j=1,i-1
      tmprate_MO(i) = tmprate_MO(i) + 2*abs(opdm((i-1)*ndim+j) * Vabs((i-1)*ndim+j) )
      tmprate_MO(i) = tmprate_MO(i) + 2*abs(opdm((j-1)*ndim+i) * Vabs((j-1)*ndim+i) )
    end do
  end do
  do i=1,ndim
    if(tmprate_MO(i).lt.0.95d0*tmprate_MO(ndim)) nva95_MO = i
    if(tmprate_MO(i).lt.0.99d0*tmprate_MO(ndim)) nva99_MO = i
  end do
  
  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmprate_MO(ndim),", nva=",nva, &
    ", nva95_MO= ",nva95_MO,", nva99_MO= ",nva99_MO


  tmprate_NO(1) = 2*abs( opdm_NO(1) * Vabs_NO(1) )
  do i=2,ndim
    tmprate_NO(i) = tmprate_NO(i-1) + 2*abs(opdm_NO((i-1)*ndim+i) * Vabs_NO((i-1)*ndim+i) )
    do j=1,i-1
      tmprate_NO(i) = tmprate_NO(i) + 2*abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j) )
      tmprate_NO(i) = tmprate_NO(i) + 2*abs(opdm_NO((j-1)*ndim+i) * Vabs_NO((j-1)*ndim+i) )
    end do
  end do
  do i=1,ndim
    if(tmprate_NO(i).lt.0.95d0*tmprate_NO(ndim)) nva95_NO = i
    if(tmprate_NO(i).lt.0.99d0*tmprate_NO(ndim)) nva99_NO = i
  end do
  
  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmprate_NO(ndim),", nva=",nva, &
    ", nva95_NO= ",nva95_NO,", nva99_NO= ",nva99_NO

  tmprate_absNO(1) = 2*abs( opdm_absNO(1) * Vabs_absNO(1) )
  do i=2,ndim
    tmprate_absNO(i) = tmprate_absNO(i-1) + 2*abs(opdm_absNO((i-1)*ndim+i) * Vabs_absNO((i-1)*ndim+i) )
    do j=1,i-1
      tmprate_absNO(i) = tmprate_absNO(i) + 2*abs(opdm_absNO((i-1)*ndim+j) * Vabs_absNO((i-1)*ndim+j) )
      tmprate_absNO(i) = tmprate_absNO(i) + 2*abs(opdm_absNO((j-1)*ndim+i) * Vabs_absNO((j-1)*ndim+i) )
    end do
  end do
  do i=1,ndim
    if(tmprate_absNO(i).lt.0.95d0*tmprate_absNO(ndim)) nva95_absNO = i
    if(tmprate_absNO(i).lt.0.99d0*tmprate_absNO(ndim)) nva99_absNO = i
  end do
  
  write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmprate_absNO(ndim),", nva=",nva, &
    ", nva95_absNO= ",nva95_absNO,", nva99_absNO= ",nva99_absNO

  !write(iout, '(A)') 'Natural Orbital Rate Analysis (Component-wise):'
  !write(iout, '(A)') '  Orb#  Occ          Vabs_NO(i,i)   Cumulative Rate'
  !write(iout, '(A)') '====================================================='
  !!do i = ndim,1, -1
  !do i = 1,ndim
  !  write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ(i), &
  !   ",  ", Vabs_NO((i-1)*ndim+i), ",  ", tmprate_NO(i)
  !end do

  !write(iout, '(A)') 'Molecular Orbital Rate Analysis (Component-wise):'
  !write(iout, '(A)') '  Orb#            Vabs_MO(i,i)   Cumulative Rate'
  !write(iout, '(A)') '================================================='
  !!do i = ndim,1, -1
  !do i = 1,ndim
  !  write(iout,'(I5, A, F12.10, A, F12.10)') i,  &
  !   ",  ", Vabs((i-1)*ndim+i), ",  ", tmprate_MO(i)
  !end do


  write(iout, '(A)') 'ABS Natural Orbital Rate Analysis (Component-wise):'
  write(iout, '(A)') '  Orb#  Occ       Vabs_absNO(i,i)   Cumulative Rate'
  write(iout, '(A)') '====================================================='
  !do i = ndim,1, -1
  do i = 1,ndim
    write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ_abs(i), &
     ",  ", Vabs_absNO((i-1)*ndim+i), ",  ", tmprate_absNO(i)
  end do

  flush(iout)    

end subroutine NO_rate_sanity2

subroutine update_maxnva( nva95maxMO, nva99maxMO, nva95maxNO, nva99maxNO, &
                          nva95maxMO_sort, nva99maxMO_sort, nva95maxNO_sort, nva99maxNO_sort, &
                          nva95max_debug, nva99max_debug, rate_debug, &
                          !nva95max_density, nva99max_density, rate_density, &
                          opdm, U_NO, Vabs, verbosity)
  implicit none
  integer(8), intent(inout) :: nva95maxMO, nva99maxMO
  integer(8), intent(inout) :: nva95maxNO, nva99maxNO
  integer(8), intent(inout) :: nva95maxMO_sort, nva99maxMO_sort
  integer(8), intent(inout) :: nva95maxNO_sort, nva99maxNO_sort
  integer(8), intent(inout) :: nva95max_debug, nva99max_debug 
  !integer(8), intent(inout) :: nva95max_density, nva99max_density 
  real(8), intent(inout) :: rate_debug  !, rate_density
  real(8), intent(in) :: opdm(:) !: One-particle density matrix in MO basis
  real(8), intent(in) :: U_NO(:) !: Transformation matrix from MO -> NO.
  !: Interesting, vabsmo is a 2D array (:,:), so I can't do Vabs(:) here.
  !:  but I can do Vabs(norb*norb)...
  real(8), intent(in) :: Vabs(norb*norb) !: <n|Vabs|m> in MO basis
  logical, intent(in), optional :: verbosity

  logical :: verbose
  integer(8) :: i, j, k, ndim, idx
  integer(8) :: nva95_NO, nva99_NO, nva95_MO, nva99_MO
  integer(8) :: nva95_NO_sort, nva99_NO_sort, nva95_MO_sort, nva99_MO_sort
  integer(8) :: nva95, nva99
  real(8) :: rate, tmp
  real(8) :: tmprate(norb*norb) !: way too big but doesnt matter
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

  if ( nva95_NO .gt. nva95maxNO) nva95maxNO = nva95_NO
  if ( nva99_NO .gt. nva99maxNO) nva99maxNO = nva99_NO


  !: For MOs ( DEBUG, UNSORTED )
  rate = 0.d0
  tmp = 0.d0
  nva95 = 0
  nva99 = 0
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
    if (tmprate(i).lt.0.95d0*tmprate(noa+nva)) nva95 = i
    if (tmprate(i).lt.0.99d0*tmprate(noa+nva)) nva99 = i
  end do
  if ( nva95 .gt. nva95maxMO) nva95max_debug = nva95
  if ( nva99 .gt. nva99maxMO) then
    nva99max_debug = nva99
    rate_debug = tmprate(noa+nva)
  end if

  !: For MOs (UNSORTED)
  rate = 0.0d0
  tmp = 0.0d0
  nva95_NO = 0
  nva99_NO = 0
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
    if (tmp/rate .lt. 0.95d0) nva95_NO = k
    if (tmp/rate .lt. 0.99d0) nva99_NO = k
  end do

  if ( nva95_NO .gt. nva95maxMO) nva95maxMO = nva95_NO
  if ( nva99_NO .gt. nva99maxMO) nva99maxMO = nva99_NO

  
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

  !if (verbose) then 
  !  write(iout, *) "NO Contribution Analysis"
  !  write(iout, *) "i, i_sort, rate, rate_sort"
  !  write(iout, *) "=========================="
  !  do i=1,ndim
  !    write(iout, "(i4, i4, f12.8, f12.8)") i, sort_key(i), sort_vals(i), vals_sorted(i)
  !  end do
  !end if

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
  if ( nva95_NO_sort .gt. nva95maxNO_sort) nva95maxNO_sort = nva95_NO_sort
  if ( nva99_NO_sort .gt. nva99maxNO_sort) nva99maxNO_sort = nva99_NO_sort



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

  !if (verbose) then
    !write(iout, *) "MO Contribution Analysis"
    !write(iout, *) "i, i_sort, rate, rate_sort"
    !write(iout, *) "=========================="
    !do i=1,ndim
    !  write(iout, "(i4, i4, f12.8, f12.8)") i, sort_key(i), sort_vals(i), vals_sorted(i)
    !end do
  !end if

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
  if ( nva95_MO_sort .gt. nva95maxMO_sort) nva95maxMO_sort = nva95_MO_sort
  if ( nva99_MO_sort .gt. nva99maxMO_sort) nva99maxMO_sort = nva99_MO_sort


  deallocate( temp )
  deallocate( temp2 )
  deallocate( opdm_NO )
  deallocate( Vabs_NO )

  if (verbose) then
    write(iout,"('For input NOs, averaged opdm (UNSORTED) : nva95max=',i5,', nva99max=',i5)") nva95maxNO, nva99maxNO
    write(iout,"('For       MOs, averaged opdm (UNSORTED) : nva95max=',i5,', nva99max=',i5)") nva95maxMO, nva99maxMO
    write(iout,"('For input NOs, averaged opdm (  SORTED) : nva95max=',i5,', nva99max=',i5)") nva95_NO_sort, nva99_NO_sort
    write(iout,"('For       MOs, averaged opdm (  SORTED) : nva95max=',i5,', nva99max=',i5)") nva95_MO_sort, nva99_MO_sort
  end if

end subroutine


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


subroutine write_density_difference(funit, time, rate, noa, nva, density, vabs)
  integer(8), intent(in) :: funit
  real(8), intent(in)    :: time, rate
  integer(8), intent(in)    :: noa, nva
  real(8), intent(in)    :: density(:)
  real(8), allocatable, intent(in)    :: vabs(:)

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
  !temp_density = temp_density * vabs(:ndim2)
  temp_density = temp_density * vabs
  !: Re-normalize density matrix
  !temptrace =  sum([(temp_density(i+(i-1)*ndim), i=1,ndim)])
  !temp_density = temp_density / temptrace

  !: Scale so max value is 2.0
  !temptrace = 0.5*maxval(abs(temp_density))
  temptrace = 0.0020426 !: scaling in ch3br at end of 045 field
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










end module propagate

