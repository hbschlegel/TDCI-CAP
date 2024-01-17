module get_ham
    
  use global_variables
  use get_ham0
  use get_ham0_cisd
  use sorting_module
  
  implicit none


  !: contains subroutine get_nstuse
  !: contains subroutine get_1eham
  !: contains subroutine get_ao2mo( mythread, myoption )
  !: contains subroutine get_form1det( mythread, myoption )
  !: contains subroutine read_restart_bin
  !: contains writeme_ham( myroutine, option)
  !: contains errors_modgetham( myroutine, option )
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_NSTUSE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_nstuse

    
    !: get the number of states to use for dynamics.  NStUse<=NStates
    
    implicit none    

    integer(8) :: i
    real(8)    :: cutoff

    
    call write_header( 'get_nstuse','get_ham','enter' )    
    
    
    !: set ionization to Koopman's value, determine number of
    if ( ionization.lt.0.d0 ) then       
       ionization = min( orben(noa),orben(nob) )
       ionization = abs(ionization)
       write(iout,'(A)') ' ionization parameter < 0.  Setting value to HOMO MO energy'
    else
       write(iout,'(A)') " ionization parameter from 'input' file "
    end if
    
    write(iout,"(' eigmax     = ',f6.2,' au', f8.2,' eV')") eigmax, eigmax*au2eV
    write(iout,"(' ionization = ',f6.2,' au', f8.2,' eV')") ionization, ionization*au2eV       
    write(iout,'(A)') ' cutoff energy = ionization + eigmax '
    write(iout,'(A)') ' all eigenstates above cutoff energy will be discarded'
    

    nstuse = nstates
    cutoff = ionization + eigmax
    
    
    !: get nstuse
    if( trim(jobtype).eq.flag_soc ) then
       do i=1, nstates
          if ( cis_eig(i).lt.cutoff ) nstuse = i
       end do
    else     
       do i=1, nstates
          if ( cis_eig(i).lt.cutoff ) nstuse = i
       end do
    end if
    
    
    call write_header( 'get_nstuse','get_ham','leave' )
    
    
  end subroutine get_nstuse
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_1EHAM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_1eham
    

    use omp_lib
    use util
    implicit none

    
    integer(8) :: i, j, k
    integer(8) :: nrorbA, nrorbB
    integer(8) :: thread1, thread2, thread3, thread4   
    real(8)    :: start1, finish1
    real(8),    allocatable :: work1(:),work2(:),work3(:),work4(:)
    complex(8), allocatable :: Zwork1(:),Zwork2(:),Zwork3(:),Zwork4(:)

    
    call write_header( 'get_1eham','get_ham','enter' )
    write(iout,'(A)') ' starting AO-->MO-->det-->CIS transformation, OMP parallelization'    
    call cpu_time(start1)
    
    nrorbA = noa + nva
    nrorbB = nob + nvb
    
    if ( Qalloc_Zcomplex ) then
       allocate( Zwork1(nstates*nstuse) )
       allocate( Zwork2(nstates*nstuse) )
       allocate( Zwork3(nstates*nstuse) )
       allocate( Zwork4(nstates*nstuse) )
    else
       allocate( work1(nstates*nstuse) )
       allocate( work2(nstates*nstuse) )
       allocate( work3(nstates*nstuse) )
       allocate( work4(nstates*nstuse) )
    end if

    !: Start OMP parallelization start non-iterative parallelization
    !$OMP PARALLEL
    !$OMP SECTIONS
    
    !$OMP SECTION  !: Vabs
    thread1 = OMP_get_thread_num()       
    write(iout,"(' Vabs Thread # ',i0)") thread1       
    call get_ao2mo( thread1, 'vabs' )   !: <AO|Vabs|AO> --> <MO|Vabs|MO>     
    call get_form1det( thread1, 'vabs') !: <MO|Vabs|MO> --> <psi_ia|Vabs|psi_jb> 
    if ( Qalloc_Zcomplex ) then
       call testherm(nstates,Zabp,'VabsMO')
!:       call Zform1cis(nstates,nstuse,Zabp,Zcis_vec)
       call Zform1cis0(nstates,nstuse,Zabp,Zcis_vec,Zwork1)
       call testherm(nstuse,Zabp,'VabsCI')
    else
!:       call form1cis(nstates,nstuse,abp,cis_vec)
       call form1cis0(nstates,nstuse,abp,cis_vec,work1)
    end if
    
    
    !$OMP SECTION  !: DipX
    thread2 = OMP_get_thread_num()
    write(iout,"(' DipX Thread # ',i0)") thread2
    call get_ao2mo( thread2, 'dipx' ) 
    call get_form1det( thread2, 'dipx') 
    if ( Qalloc_Zcomplex ) then
       call testherm(nstates,Ztdx,'TDxMO')
!:       call Zform1cis(nstates,nstuse,Ztdx,Zcis_vec)
       call Zform1cis0(nstates,nstuse,Ztdx,Zcis_vec,Zwork2)
       call testherm(nstuse,Ztdx,'TDxCI')
    else 
!:       call form1cis(nstates,nstuse,tdx,cis_vec)        
       call form1cis0(nstates,nstuse,tdx,cis_vec,work2)
    end if

    !$OMP SECTION  !: DipY
    thread3 = OMP_get_thread_num()
    write(iout,"(' DipY Thread # ',i0)") thread3       
    call get_ao2mo( thread3, 'dipy' ) 
    call get_form1det( thread3, 'dipy')
    if ( Qalloc_Zcomplex ) then
       call testherm(nstates,Ztdy,'TDyMO')
!:       call Zform1cis(nstates,nstuse,Ztdy,Zcis_vec)
       call Zform1cis0(nstates,nstuse,Ztdy,Zcis_vec,Zwork3)
       call testherm(nstuse,Ztdy,'TDyCI')
    else 
!:       call form1cis(nstates,nstuse,tdy,cis_vec)  
       call form1cis0(nstates,nstuse,tdy,cis_vec,work3)
    end if
    
    !$OMP SECTION  !: DipZ
    thread4 = OMP_get_thread_num()
    write(iout,"(' DipZ Thread # ',i0)") thread4
    call get_ao2mo( thread4, 'dipz' )
    call get_form1det( thread4, 'dipz' )
    if ( Qalloc_Zcomplex ) then
       call testherm(nstates,Ztdz,'TDzMO')
!:       call Zform1cis(nstates,nstuse,Ztdz,Zcis_vec)
       call Zform1cis0(nstates,nstuse,Ztdz,Zcis_vec,Zwork4)
       call testherm(nstuse,Ztdz,'TDzCI')
    else 
!:       call form1cis(nstates,nstuse,tdz,cis_vec)
       call form1cis0(nstates,nstuse,tdz,cis_vec,work4)
    end if
       
    !$OMP END SECTIONS
    !$OMP END PARALLEL
   
    !: AD
    write(iout, *) "Eq 18 (Field Perturbed Orbitals) OmniDirectional:"
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.030d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.035d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.040d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.045d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.050d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.055d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.060d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.065d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.070d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.075d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.080d0, 1.d0, 1.d0, 1.d0)
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.090d0, 1.d0, 1.d0, 1.d0)
    write(iout, *) "Eq 18 (Field Perturbed Orbitals) Directional:"
    write(iout, *) "Emax = 0.050 au"
    write(iout, *) "Theta=0"
    call generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, 0.050d0, 1.d0, 0.d0, 0.d0)
    write(iout, *) "Theta=45"
    call generate_field_perturbed_orbitals( &
      norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, &
      0.050d0, 1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0), 0.d0)
    write(iout, *) "Theta=90"
    call generate_field_perturbed_orbitals( &
      norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, &
      0.050d0, 0.d0, 1.d0, 0.d0)
    write(iout, *) "Theta=135"
    call generate_field_perturbed_orbitals( &
      norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, &
      0.050d0, -1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0), 0.d0)
    write(iout, *) "Theta=180"
    call generate_field_perturbed_orbitals( &
      norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, &
      0.050d0, -1.d0, 1.d0/sqrt(2.d0), 0.d0)

    write(iout, *) "Emax = 0.035"
    write(iout, *) "Field directon (0, 1/sqrt(2), 1/sqrt(2))"
    call generate_field_perturbed_orbitals( &
      norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabsmoa, &
      0.0350d0, 0.d0, 1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0))
    write(iout, *) "Diagonalize EFock::"
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.002d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.005d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.010d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.030d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.035d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.040d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.070d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.080d0 )
    call diag_EFock( dipxmoa, dipymoa, dipzmoa, vabsmoa,  orben, 0.090d0 )

    call write_header( 'get_1eham','get_ham','leave' )
    call cpu_time(finish1)
    write(iout,"(' One electorn integral transformation time:',f12.4,' seconds')") finish1-start1 

    
    
  end subroutine get_1eham

  subroutine generate_field_perturbed_orbitals( norb, noa, nva, orben, dipxmoa, dipymoa, dipzmoa, vabs_a, Emax, xscale, yscale, zscale)

    implicit none

    !: dummy variables
    integer(8),        intent(in) :: norb, noa, nva
    real(8),              intent(in) :: orben(norb)
    real(8),              intent(in) :: dipxmoa(noa+nva, noa+nva)
    real(8),              intent(in) :: dipymoa(noa+nva, noa+nva)
    real(8),              intent(in) :: dipzmoa(noa+nva, noa+nva)
    real(8),              intent(in) :: vabs_a(noa+nva, noa+nva)
    real(8),              intent(in) :: Emax
    real(8),              intent(in) :: xscale, yscale, zscale
    !: Should I be multiplying C_out with cmo_a ? 
    !real,              intent(in) :: cmo_a(nbasis, noa+nva) !: AO -> MO 

    !: function result
    !: local variables
    integer(8) :: i,a,p, nva95, nva99
    !real(8)    :: Emax = 0.05d0 !: 0.005339 = ~1E+16 W/m^2
    real(8)    :: tmpval = 0.d0
    real(8)    :: tmpsum = 0.d0
    real(8)    :: dipole_norm = 0.d0
    real(8)    :: C_out(noa+nva, noa+nva)
    real(8)    :: mo_rate(noa+nva)
    real(8)    :: mo_rate_sorted(noa+nva)
    integer(8) :: sort_key(noa+nva)

    !: Make sure array is zero'd
    C_out = 0.d0
    mo_rate = 0.d0
    mo_rate_sorted = 0.d0
    nva95 = 0
    nva99 = 0

    write(iout, *) "(noa, nva): ", noa, nva

    do i=1,(noa+nva) !: c_ij = delta_ij, c_ab = delta_ab. 
      C_out(i,i) = 1.d0
    end do

    !: off-diagonals
    do i=1,noa
      do a=noa+1,noa+nva
        dipole_norm = sqrt( (xscale*dipxmoa(a,i))**2 &
                          + (yscale*dipymoa(a,i))**2 + (zscale*dipzmoa(a,i))**2  )
        tmpval = Emax*dipole_norm/( orben(a)-orben(i) )
        C_out(a,i) = tmpval
        C_out(i,a) = tmpval
      end do
    end do

    do a=1,nva
      !: Calculate ionization rate for orbital a
      mo_rate(a) = 0.d0
      do i=1,noa
        do p=1,noa+nva
          !: debug prints
          !write(iout, '(A, I4, A, I4, A, I4,)') '(a,i,p) = ', a, ", ", i, ", ", p
          !write(iout, '(F10.7, A, F10.7, A, F10.7, A, F10.7)') mo_rate(a), ", ", C_out(p,i), &
          !  ", ",C_out(a,i), ", ",vabs_a(p,a)
          mo_rate(a) = mo_rate(a) + abs(C_out(p,i)*C_out(a,i)*vabs_a(p,a))
        end do
      end do
    end do

    !write(iout, '(A)') NEW_LINE('A')
    !write(iout, '(A)') 'Perturbative Rate Analysis:'
    !write(iout, '(A)') 'Orb#   Rate   <i|Vabs|i>'
    !write(iout, '(A)') '========================'
    !do i = 1,nva
    !  write(iout,'(I5, A, F10.7, A, F10.7)') i, ",  ", mo_rate(i), ",  ", vabs_a(noa+i,noa+i)
    !end do

    tmpval = 0.d0
    tmpsum = 0.d0
    !: Sum rate for percent denominator
    !do i = 1,noa+nva
    write(iout, '(A)') 'Perturbative Rate Analysis:'
    write(iout, '(A)') 'i   tmpsum'
    write(iout, '(A)') '=========='
    do i = 1,nva
      tmpsum = tmpsum + mo_rate(i)
      write(iout, '(I5, A, F15.10)') i, ", ", tmpsum
    end do
    !: Calculate cumulative rate by unsorted MOs
    !do i = 1,noa+nva
    do i = 1,nva
      tmpval = tmpval + mo_rate(i)/tmpsum
      if(tmpval.lt.0.95d0) nva95 = i
      if(tmpval.lt.0.99d0) nva99 = i
    end do
    write(iout, '(A, F10.7)') 'Perturbative Rate Summary:  Emax = ', Emax
    write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmpsum,", nva=",nva, &
      ", (UNSORTED) nva95= ",nva95,", nva99= ",nva99

    call quicksort_descending(mo_rate, mo_rate_sorted, sort_key)

    !: Calculate cumulative rate by sorted MOs
    tmpval = 0.d0
    do i = 1,noa+nva
      tmpval = tmpval + mo_rate_sorted(i)/tmpsum
      if(tmpval.lt.0.95d0) nva95 = i
      if(tmpval.lt.0.99d0) nva99 = i
    end do
    write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmpsum,", nva=",nva, &
      ",   (SORTED) nva95= ",nva95,", nva99= ",nva99


  end subroutine generate_field_perturbed_orbitals

  !: Fock is diagonal in HF MO basis, so just use orbital energies.
  subroutine diag_EFock( dipxmoa, dipymoa, dipzmoa,  Vabs, orb_eng, Emax )
    implicit none

    real(8), intent(in) :: dipxmoa(noa+nva, noa+nva)
    real(8), intent(in) :: dipymoa(noa+nva, noa+nva)
    real(8), intent(in) :: dipzmoa(noa+nva, noa+nva)
    real(8), intent(in) :: Vabs((noa+nva),(noa+nva))
    real(8), intent(in) :: orb_eng(ndim)
    real(8), intent(in) :: Emax

    integer :: ndim, i,j,k,a,p, nva95, nva99
    logical :: Vabs_isSymmetric
    !real(8) :: Emax = 0.001d0 !: 0.005339 = ~1E+16 W/m^2

    real(8) :: FieldFock((noa+nva)*(noa+nva))
    real(8) :: evals((noa+nva))
    real(8)    :: tmpval = 0.d0
    real(8)    :: tmpsum = 0.d0
    real(8)    :: dipole_norm = 0.d0
    real(8)    :: mo_rate(noa+nva)
    real(8)    :: mo_rate_sorted(noa+nva)

    !: dsyevd stuff
    integer :: info_, lwork, liwork
    integer, allocatable :: iwork(:)
    real(8), allocatable ::  work(:)
    integer(8) :: sort_key((noa+nva))

    ndim = noa+nva

    !: Sanity checks
    !call isSymmetric( Vabs, ndim, Vabs_isSymmetric )
    !if (Vabs_isSymmetric) then
    !  write(iout, '(A)') 'Vabs_moa is symmetric'
    !else
    !  write(iout, '(A)') 'Vabs_moa is NOT symmetric'
    !end if
    !call matrixIndexSanity( Vabs, Vabs, ndim )

    !: Construct FieldFock to be diagonalized
    FieldFock = 0.d0
    evals = 0.d0
    !:  Copy orb energies
    do i=1,ndim
      FieldFock( (i-1)*ndim+i ) = orb_eng(i)
    end do
    !: Add field
    do i=1,ndim
      do j=1,ndim
        tmpval = Emax*sqrt( dipxmoa(j,i)**2 + dipymoa(j,i)**2 + dipzmoa(j,i)**2 )
        FieldFock( (i-1)*ndim+j ) = FieldFock( (i-1)*ndim+j ) - tmpval
      end do
    end do

    !call isSymmetric( FieldFock, ndim, Vabs_isSymmetric )
    !if (Vabs_isSymmetric) then
    !  write(iout, '(A)') 'FieldFock is symmetric'
    !else
    !  write(iout, '(A)') 'FieldFock is NOT symmetric'
    !end if

    !: Diagonalize
    !: workspace query to determine size for work arrays
    lwork = -1
    liwork = -1
    allocate(work(1))
    allocate(iwork(1))
    call dsyevd('v','u', ndim, FieldFock, ndim, evals, work, lwork, iwork, liwork, info_)
    lwork = nint(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)

    allocate(work(lwork))
    allocate(iwork(liwork))
    call dsyevd('v','u', ndim, FieldFock, ndim, evals, work, lwork, iwork, liwork, info_)

    write(iout,'(A)') 'After diagonalzie FieldFock'
    if (info_) then
      write(iout, *) 'dsyevd error!!: ', info_
    end if
    deallocate(work, iwork)

    !: Calculate rate
    do a=1,noa+nva
      mo_rate(a) = 0.d0
      do i=1,noa
        do p=1,noa+nva
          !mo_rate(a) = mo_rate(a) + abs( FieldFock  (p,i)*FieldFock(a,i) * Vabs(p,a) )
          mo_rate(a) = mo_rate(a) + abs( FieldFock( (i-1)*ndim+p )*FieldFock( (i-1)*ndim+a ) * Vabs(p,a) )
        end do
      end do
    end do

    !write(iout, '(A)') NEW_LINE('A')
    !write(iout, '(A)') 'diag(FieldFock) Rate Analysis:'
    !write(iout, '(A)') 'Orb#   Rate   <i|Vabs|i>'
    !write(iout, '(A)') '========================'
    !do i = 1,noa+nva
    !  write(iout,'(I5, A, F10.7, A, F10.7)') i, ",  ", mo_rate(i), ",  ", Vabs(i,i)
    !end do

    tmpval = 0.d0
    tmpsum = 0.d0
    !: Sum rate for percent denominator
    do i = 1,noa+nva
      tmpsum = tmpsum + mo_rate(i)
    end do
    !: Calculate cumulative rate by unsorted MOs
    do i = 1,noa+nva
      tmpval = tmpval + mo_rate(i)/tmpsum
      if(tmpval.lt.0.95d0) nva95 = i
      if(tmpval.lt.0.99d0) nva99 = i
    end do
    write(iout, '(A, F10.7)') 'diag(FieldFock) Rate Summary:  Emax = ', Emax
    write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmpsum,", nva=",nva, &
      ", (UNSORTED) nva95= ",nva95,", nva99= ",nva99

    call quicksort_descending(mo_rate, mo_rate_sorted, sort_key)

    !: Calculate cumulative rate by sorted MOs
    tmpval = 0.d0
    do i = 1,noa+nva
      tmpval = tmpval + mo_rate_sorted(i)/tmpsum
      if(tmpval.lt.0.95d0) nva95 = i
      if(tmpval.lt.0.99d0) nva99 = i
    end do
    write(iout, '(A, F10.7, A, I5, A, I5, A, I5)') " rate= ",tmpsum,", nva=",nva, &
      ",   (SORTED) nva95= ",nva95,", nva99= ",nva99
    
  end subroutine

  subroutine isSymmetric( A, n, output )
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n)
    logical, intent(out) :: output

    integer :: i,j

    output = .true.
    do i=1,n
      do j=1,i-1
        if ( abs(A(i,j)-A(j,i)) .gt. 0.0001 ) then
          output = .false.
        end if
      end do
    end do

  end subroutine isSymmetric

  !: Feed a 2d matrix into both A_1D and A_2D
  subroutine matrixIndexSanity( A_1D, A_2D, n )
    implicit none

    integer, intent(in) :: n
    real(8), intent(in) :: A_1D(n*n), A_2D(n,n)

    integer :: i,j
    logical :: insane
    insane = .false.
    do i=1,n
      do j=1,n
        if ( abs( A_2D(i,j)-A_1D((i-1)*n+j) ) .gt. 0.0001 ) then
          insane = .true.
        end if
      end do
    end do

  if (insane) then
    write(iout, '(A)') 'MATRIX INDEX USAGE NOT SANE!!'
    write(iout, '(A)') 'A_1D:'
    do i = 1, n
      do j = 1, n
        write(iout, *, advance='no') A_1D((i-1)*n+j)
      end do
      write(iout, *) ! New line after each row
    end do
    write(iout, '(A)') 'A_2D:'
    do i = 1, n
      do j = 1, n
        write(iout, *, advance='no') A_2D(i,j)
      end do
      write(iout, *) ! New line after each row
    end do
  else
    write(iout, '(A)') 'MATRIX INDEX USAGE SANITY CHECKED.'
  end if

  end subroutine matrixIndexSanity


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_AO2MO
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_ao2mo( mythread, myoption )

    !: intermediary routine to call ao2mo 
    !: no real need, just makes the code look cleaner
    

    use util
    implicit none

    integer(8),   intent(in) :: mythread
    character(*), intent(in) :: myoption
    real(8)    :: start, finish


    call cpu_time(start)
 
    select case( trim(myoption) )
    case( 'vabs' )  ; call ao2mo(nbasis,nrorb,vabsao,vabsmoa,cmo_a)
    case( 'dipx' )  ; call ao2mo(nbasis,nrorb,dipxao,dipxmoa,cmo_a)
    case( 'dipy' )  ; call ao2mo(nbasis,nrorb,dipyao,dipymoa,cmo_a)
    case( 'dipz' )  ; call ao2mo(nbasis,nrorb,dipzao,dipzmoa,cmo_a)
    end select
    
    if (unrestricted ) then    
       select case ( trim(myoption) )
       case( 'vabs' )  ; call ao2mo(nbasis,nrorb,vabsao,vabsmob,cmo_b) 
       case( 'dipx' )  ; call ao2mo(nbasis,nrorb,dipxao,dipxmob,cmo_b)       
       case( 'dipy' )  ; call ao2mo(nbasis,nrorb,dipyao,dipymob,cmo_b)
       case( 'dipz' )  ; call ao2mo(nbasis,nrorb,dipzao,dipzmob,cmo_b)
       end select
    end if

    !: < 0 | V | 0 > element
    if( unrestricted ) then
       call det00(noa,nva,vabsmoa,vabs00,nob,nvb,vabsmob) 
    else
       call det00(noa,nva,vabsmoa,vabs00)  
       vabs00 = 2.d0 * vabs00
    end if

    call cpu_time(finish)
    write(iout,"(' Thread # ',i0, a6,' AO-->MO  time: ',f12.4,' s')") mythread, trim(myoption), finish-start
    flush(iout)
    

  end subroutine get_ao2mo
  !:---------------------------!
  !: SUBROUTINE GET_SOC_AO2MO
  !:---------------------------!
!  subroutine get_soc_ao2mo
!
!    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS
!    !: Form complex spin-orbit matrices VPM,VZA,VZB in MO basis
!    !: real VX, VY and VZ integrals in AO basis
!    
!    !: 
!    !: AO basis, lower triangle, real
!    !: VX=(mu/SOCx/nu)/i, VY=(mu/SOCy/nu)/i, VZ=(mu/SOCz/nu)/i
!    !:
!    !: MO basis, full matrix, complex (Z, raising and lowering)
!    !:    VZA = (p/SOCz/q) for alpha,alpha
!    !:    VZB = (p/SOCz/q) for beta,beta
!    !:    VPM = (p/SOCx/q)+i(p/SOCy/q) for alpha,beta
!    !:    VM(p,q) = VPM(p,q), VP(p,q) = VPM(q,p)*
!    !:
!    !: conversion factor for spin-orbit integrals
!    !: (see Salem, Angnew Chem Internat - Vol 11 (1972),No 2,pp92-111)
!    !: constant=((h/2pi)**2*e**2)/4*m***2*c**2 = 2.9217cm**-1
!    !: = 0.00001331224 au
!    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS
!
!    use util
!
!    implicit none
!    real(8), parameter :: constant = 0.00001331224d0
!    
!    integer(8) :: i,j,ij, ntt
!    complex(8), allocatable :: scratch1(:,:)
!
!    ntt = nbasis*(nbasis+1)/2
!
!    !: convert to au units
!    socxao(:) = constant * socxao(:)
!    socyao(:) = constant * socyao(:)
!    soczao(:) = constant * soczao(:)
!    
!
!    allocate( scratch1(nbasis,nbasis) )
!    scratch1 = 0.d0    
!    
!    !: socmoAA
!    do i=1, nbasis
!       do j=1, i
!          ij = i*(i-1)/2 + j
!          scratch1(i,j) = dcmplx( 0.d0,-soczao(ij) )
!          scratch1(j,i) = dconjg( scratch1(i,j) )
!       end do
!    end do
!    call ao2mo_complex(nbasis, nrorb, scratch1, socmoAA, cmo_a, cmo_a )
!
!    open( unit=100,file='VZA.OUT' )
!    do i=1, nrorb
!       do j=1, nrorb
!          write(100,"(f15.10,f15.10)") socmoAA(j,i)
!       end do
!    end do
!    close(100)
!    
!    
!    !: socmoBB
!    scratch1 = dcmplx( 0.d0, 0.d0 )
!    do i=1, nbasis
!       do j=1, i
!          ij = i*(i-1)/2 + j
!          scratch1(i,j) = dcmplx( 0.d0,soczao(ij) )
!          scratch1(j,i) = dconjg( scratch1(i,j) )
!       end do
!    end do
!    call ao2mo_complex(nbasis, nrorb, scratch1, socmoBB, cmo_b, cmo_b )
!
!    open( unit=100,file='VZB.OUT' )
!    do i=1, nrorb
!       do j=1, nrorb
!          write(100,"(f15.10,f15.10)") socmoBB(j,i)
!       end do
!    end do
!    close(100)    
!
!
!    !: socmoAB
!    do i=1, nbasis
!       do j=1, i
!          ij = i*(i-1)/2 + j 
!          scratch1(i,j) = dcmplx( -socyao(ij), -socxao(ij) )
!          scratch1(j,i) = dcmplx( socyao(ij), socxao(ij) )
!       end do
!    end do
!    call ao2mo_complex(nbasis, nrorb, scratch1, socmoAB, cmo_a, cmo_b )
!
!
!    open( unit=100,file='VPM.OUT' )
!    do i=1, nrorb
!       do j=1, nrorb
!          write(100,"(f15.10,f15.10)") socmoAB(j,i)
!       end do
!    end do
!    close(100)    
!    
!
!    deallocate( scratch1 )
!
!
!
!  end subroutine get_soc_ao2mo
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_FORM1H
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_form1det( mythread, myoption)
    
    !: intermediary routine to call form1h
    !: no real need, just makes the code look cleaner
    
    use util
    implicit none

    integer(8),   intent(in) :: mythread
    character(*), intent(in) :: myoption
    real(8)    :: start, finish
    

    call cpu_time(start)

    select case( trim(jobtype) )
    case( flag_cis )     
       if ( unrestricted ) then
          select case ( trim(myoption) )
          case( 'vabs' ) ; call form1det(noa,nva,nstates,abp,vabsmoa,vabs00,hole_index(:,1),part_index(:,1),nob,nvb,vabsmob)
          case( 'dipx' ) ; call form1det(noa,nva,nstates,tdx,dipxmoa,dipx00,hole_index(:,1),part_index(:,1),nob,nvb,dipxmob)
          case( 'dipy' ) ; call form1det(noa,nva,nstates,tdy,dipymoa,dipy00,hole_index(:,1),part_index(:,1),nob,nvb,dipymob)
          case( 'dipz' ) ; call form1det(noa,nva,nstates,tdz,dipzmoa,dipz00,hole_index(:,1),part_index(:,1),nob,nvb,dipzmob)
          end select
       else
          select case ( trim(myoption) )
          case( 'vabs' ) ; call form1det(noa,nva,nstates,abp,vabsmoa,vabs00,hole_index(:,1),part_index(:,1))
          case( 'dipx' ) ; call form1det(noa,nva,nstates,tdx,dipxmoa,dipx00,hole_index(:,1),part_index(:,1))
          case( 'dipy' ) ; call form1det(noa,nva,nstates,tdy,dipymoa,dipy00,hole_index(:,1),part_index(:,1))
          case( 'dipz' ) ; call form1det(noa,nva,nstates,tdz,dipzmoa,dipz00,hole_index(:,1),part_index(:,1))
          end select
       end if
    case( flag_soc ) 
       select case( trim(myoption) )
       case( 'vabs' ) ; call Zform1det(noa,nva,nstates,Zabp,vabsmoa,vabs00,hole_index(:,1),part_index(:,1),nob,nvb,vabsmob)
       case( 'dipx' ) ; call Zform1det(noa,nva,nstates,Ztdx,dipxmoa,dipx00,hole_index(:,1),part_index(:,1),nob,nvb,dipxmob)
       case( 'dipy' ) ; call Zform1det(noa,nva,nstates,Ztdy,dipymoa,dipy00,hole_index(:,1),part_index(:,1),nob,nvb,dipymob)
       case( 'dipz' ) ; call Zform1det(noa,nva,nstates,Ztdz,dipzmoa,dipz00,hole_index(:,1),part_index(:,1),nob,nvb,dipzmob)
       end select
    case( flag_tda )
       select case ( trim(myoption) )
       case( 'vabs' ) ; call form1det(noa,nva,nstates,abp,vabsmoa,vabs00,hole_index(:,1),part_index(:,1),nob,nvb,vabsmob)
       case( 'dipx' ) ; call form1det(noa,nva,nstates,tdx,dipxmoa,dipx00,hole_index(:,1),part_index(:,1),nob,nvb,dipxmob)
       case( 'dipy' ) ; call form1det(noa,nva,nstates,tdy,dipymoa,dipy00,hole_index(:,1),part_index(:,1),nob,nvb,dipymob)
       case( 'dipz' ) ; call form1det(noa,nva,nstates,tdz,dipzmoa,dipz00,hole_index(:,1),part_index(:,1),nob,nvb,dipzmob)
       end select
    case( flag_ip ) 
       select case ( trim(myoption) )
       case( 'vabs' ) ; call form1det_ip(noa,nva,nstates,abp,vabsmoa,vabs00,hole_index,part_index(:,1),nob,nvb,vabsmob)
       case( 'dipx' ) ; call form1det_ip(noa,nva,nstates,tdx,dipxmoa,dipx00,hole_index,part_index(:,1),nob,nvb,dipxmob)
       case( 'dipy' ) ; call form1det_ip(noa,nva,nstates,tdy,dipymoa,dipy00,hole_index,part_index(:,1),nob,nvb,dipymob)
       case( 'dipz' ) ; call form1det_ip(noa,nva,nstates,tdz,dipzmoa,dipz00,hole_index,part_index(:,1),nob,nvb,dipzmob)
       end select
    case( flag_socip ) 
       select case( trim(myoption) )
       case( 'vabs' ) ; call Zform1det_ip(noa,nva,nstates,Zabp,vabsmoa,vabs00,hole_index,part_index(:,1),nob,nvb,vabsmob)
       case( 'dipx' ) ; call Zform1det_ip(noa,nva,nstates,Ztdx,dipxmoa,dipx00,hole_index,part_index(:,1),nob,nvb,dipxmob)
       case( 'dipy' ) ; call Zform1det_ip(noa,nva,nstates,Ztdy,dipymoa,dipy00,hole_index,part_index(:,1),nob,nvb,dipymob)
       case( 'dipz' ) ; call Zform1det_ip(noa,nva,nstates,Ztdz,dipzmoa,dipz00,hole_index,part_index(:,1),nob,nvb,dipzmob)
       end select
    end select
    
    
    call cpu_time(finish)
    write(iout,"(' Thread # ',i0, a6,' MO-->det time: ',f12.4,' s')") mythread, trim(myoption), finish-start
    flush(iout)


  end subroutine get_form1det
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_RESTART_BIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_restart_bin
    
    implicit none

    logical    :: iamhere
    integer(8) :: read_nstates, i, j


    call write_header( 'read_restart_bin','get_ham','enter' )
    
    inquire( file=trim(restartbinfile), exist=iamhere )
    if ( .not.iamhere ) call errors_modgetham( 'restart', 'nofile' )

    open( unit=50, file=trim(restartbinfile), form='unformatted' )

    read(50) read_nstates, nstuse
    write( iout,"(' nstates = ',i0,' nstuse = ',i0)" ) read_nstates, nstuse
    if( read_nstates .ne. nstates ) call errors_modgetham( 'restart', 'mismatch' )
    if( nstates .lt. nstuse )       call errors_modgetham( 'restart', 'nstuse' )
    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
      read(50) Zcis_vec
      read(50) cis_eig
      read(50) Ztdx
      read(50) Ztdy
      read(50) Ztdz
      read(50) Zabp
      read(50) Zexp_abp
      read(50) Zip_vec
    else
      read(50) cis_vec
      read(50) cis_eig
      read(50) tdx
      read(50) tdy
      read(50) tdz
      read(50) abp
      read(50) exp_abp
      read(50) ip_vec
    end if
    close(50)
  
    write(iout,*) " ip_states",ip_states 
    write(iout,*) (Zip_vec(i),i=1,10)
    if(Qalloc_vabsmo.and.Qalloc_dipmo) then
      open( unit=50, file=trim(outputfile)//'_RESTART_MO.bin', form='unformatted' )

      read(50) i,i,i,i 
      read(50) ( vabsmoa(:,i), i=1, nrorb )
!:      write(42,*) ( vabsmoa(:,i), i=1, nrorb )
      read(50) ( dipxmoa(:,i), i=1, nrorb )
      read(50) ( dipymoa(:,i), i=1, nrorb )
      read(50) ( dipzmoa(:,i), i=1, nrorb )
      if ( unrestricted ) then
        read(50) ( vabsmob(:,i), i=1, nrorb )
        read(50) ( dipxmob(:,i), i=1, nrorb )
        read(50) ( dipymob(:,i), i=1, nrorb )
        read(50) ( dipzmob(:,i), i=1, nrorb )
      end if
      if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
        read(50) ( socmoAA(:,i), i=1, nrorb )
        read(50) ( socmoBB(:,i), i=1, nrorb )
        read(50) ( socmoAB(:,i), i=1, nrorb )
      end if
    end if
    close(50)

    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       read(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
!:    if( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
      do i=1, ip_states
        read(100,"(i7,i7,i7)") hole_ip_index(i,1), hole_ip_index(i,2), part_ip_index(i,1)
      end do
      do i=1, noa+nob
        read(100,"(400i5)") (state_ip_index(i,j), j=1,noa+nob)
      end do
!:    end if
    close(100)

    write( iout,"(' nstates = ',i0,' nstuse = ',i0)" ) read_nstates, nstuse
    write( iout,"(' finished reading cis_vec, cis_eig, tdx, tdy, tdz, abp, exp_abp')")
    
    !select case( trim(jobtype) ) 
    !case( flag_cis ) ; call get_cis_index
    !case( flag_tda ) ; call get_cis_index
    !case( flag_ip )  ; call get_ip_index       
    !end select

    call write_header( 'read_restart_bin','get_ham','leave')


  end subroutine read_restart_bin
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_GETHAM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_ham(myroutine,option)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option


    select case( trim(myroutine) )
    end select
    
    flush(iout)
    

  end subroutine writeme_ham
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE errors_modget_ham
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_modgetham(myroutine,option)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option


    select case( trim(myroutine) )
    case( 'restart' )
       select case( trim(option) )
       case('nofile')
          write(iout,'(A)' ) ' ERROR:  could not find restart binary file'
          go to 100
       case('mismatch') 
          write(iout,'(A)' ) ' ERROR:  unexpected value for nstates'
          go to 100
       case('nstuse') 
          write(iout,'(A)' ) ' ERROR:  nstuse > nstates'
          go to 100
       end select
    end select

100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " DON'T BE DISMAYED BY GOODBYES - Richard Bach "
    stop

200 continue
          

  end subroutine errors_modgetham
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module get_ham
