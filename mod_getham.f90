module get_ham
    
  use global_variables
  use get_ham0
  use get_ham0_cisd
  
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
   
    call write_header( 'get_1eham','get_ham','leave' )
    call cpu_time(finish1)
    write(iout,"(' One electorn itegral transformation time:',f12.4,' seconds')") finish1-start1 

    
    
  end subroutine get_1eham
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
