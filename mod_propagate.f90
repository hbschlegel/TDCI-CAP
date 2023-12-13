module propagate
  
  use global_variables
  use analyze_psi
  use util
  use sorting_module
  use io_binary  ! io_bin_test, write_dbin, read_dbin
  
  implicit none


contains
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
    !   write(iout,"(' heuristic lifetime =', f8.4,' sqrt(cis_eig-ionization)')") heuristic
    !   do i = 1, nstuse
    !      if( cis_eig(i).gt.ionization ) abp(nstuse*(i-1)+i) = abp(nstuse*(i-1)+i)+heuristic*sqrt(cis_eig(i)-ionization)
    !   end do
    !else
    !   write(iout,'(A)') " heuristic = 0.0 "
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
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse) 
    real(8), allocatable :: hp2(:), tdvals2(:)        
    
    !: temporary arrays to be deallocated
    real(8) :: norm0
    real(8),allocatable :: pop0(:),pop1(:),ion(:)
    real(8),allocatable :: rate_a(:),rate_b(:),rate_aa(:),rate_ab(:),rate_ba(:),rate_bb(:)
    complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)


    !: private variables
    integer(8) :: i, j, ii, jj, k, kk
    integer(8) :: itime, idir, iemax
    integer(8) :: ithread, idata
    complex(8) :: cdum
    
    !: field info
    real(8) :: dirx1, diry1, dirz1, emax1, efield1
    real(8) :: dirx2, diry2, dirz2, emax2, efield2
    real(8) :: efieldx, efieldy, efieldz
    real(8) :: temp, temp1, temp2

    !: results and psi stuff
    real(8)    :: norm, normV, rate, mux, muy, muz
    complex(8) :: psi_j, psi_k
    complex(8) :: psi(nstuse), psi1(nstates)

    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, info3, lscratch, liwork
    real(8)    :: start1, start2, finish1, finish2
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: scratch(:)

    !: Natural Orbital generation
    real(8), allocatable :: opdm_avg(:)  !: Averaged one-particle reduced density matrix (1-RDM)
    real(8), allocatable :: opdm_avg_abs(:)  !: average(abs(opdm))
    real(8), allocatable :: natorb_occ(:) !: Natural orbital occupations (eigenvalues)
    real(8), allocatable :: natorb_occ_abs(:) !: Natural orbital occupations (eigenvalues)
    real(8), allocatable :: U_NO(:) !: MO->NO transformation matrix
    real(8), allocatable :: U_NO_abs(:) !: MO->NO_abs transformation matrix
    integer(8) :: opdm_avg_N, ndim
    real(8), allocatable :: U_NO_input(:) !: U_NO read from file
    integer(8) :: nva95max, nva99max, nva95maxmax, nva99maxmax



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

    !: Natural Orbital initialization
    ndim = noa+nva
    allocate( opdm_avg(ndim*ndim) )
    allocate( opdm_avg_abs(ndim*ndim) )
    allocate(natorb_occ(ndim))
    allocate(natorb_occ_abs(ndim))
    allocate( U_NO(ndim*ndim) )
    allocate( U_NO_abs(ndim*ndim) )
    allocate( U_NO_input(ndim*ndim) )
    opdm_avg = 0.d0
    opdm_avg_abs = 0.d0
    natorb_occ = 0.d0
    natorb_occ_abs = 0.d0
    U_NO = 0.d0
    U_NO_abs = 0.d0
    opdm_avg_N = ndir*(nstep/outstep)
    write(iout, '(A,i5)') "opdm_avg_N: ", opdm_avg_N
    nva95max = 0
    nva99max = 0
    nva95maxmax = 0
    nva99maxmax = 0
    !if ( .true. ) then
    if ( flag_ReadU_NO ) then
      call read_dbin( U_NO_input, (ndim*ndim), "U_NO_input.bin", info3)
    end if

    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )    
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*ip_states*(ip_states+1)))
    else
      allocate(Zion_coeff(ip_states*(ip_states+1)))
      allocate(Zproj_ion(ip_states*ip_states))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( scratch(1+8*nstuse+2*nstuse*nstuse) )
    else
      allocate( iwork(2) )
      allocate( scratch(nstuse*nstuse) )
    end if

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
    
    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do


    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    !$OMP PARALLEL DEFAULT(NONE),&
    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
    !$OMP cifile,datafile,dirstr,emaxstr,finish1,finish2,funit1,funit2,funit3,funit4,funit5,funit6, &
    !$OMP info1,info2,ithread,lscratch,liwork,start1,start2, &
    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, &
    !$OMP psi_j, psi_k, temp, temp1, temp2, cdum,           &
    !$OMP norm, norm0, normV, mux, muy, muz, rate, efieldx, efieldy, efieldz, &
    !$OMP pop1, ion, ion_coeff, rate_a, rate_b, rate_aa, rate_ab, rate_ba, rate_bb, psi_det0,  &
    !$OMP hp1, hp2, psi, psi1, scratch, iwork, tdvals1, tdvals2, Zion_coeff, &
    !$OMP nva95max, nva99max ),  &
    !$OMP SHARED( jobtype, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
    !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
    !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
    !$OMP vabsmoa, vabsmob, noa, nob, nva, nvb, norb, hole_index, part_index, &
    !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
    !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
    !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
    !$OMP opdm_avg, opdm_avg_abs, opdm_avg_N, U_NO, U_NO_abs, natorb_occ, natorb_occ_abs, cmo_a, &
    !$OMP U_NO_input, flag_ReadU_NO, nva95maxmax, nva99maxmax )
    
    !$OMP DO  

    dir_loop : do idir=1, ndir
       
       ithread = omp_get_thread_num()
       call cpu_time( start1 )

       !: reset nva max
       nva95max = 0
       nva99max = 0

       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x0  
       diry1 = tdciresults(idir)%y0  
       dirz1 = tdciresults(idir)%z0  

       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*tdx(1:nstuse2) + diry1*tdy(1:nstuse2) + dirz1*tdz(1:nstuse2)
       
       !: diagonalize mu dot E
       info1 = 10  
       If( QeigenDC ) then
         lscratch = 1+6*nstuse+2*nstuse*nstuse
         liwork = 3+5*nstuse
         call dsyevd('v','u',nstuse, hp1, nstuse, tdvals1, scratch, lscratch, iwork, liwork, info1)
       else
         lscratch = nstuse*nstuse
         call dsyev('v','u',nstuse, hp1, nstuse, tdvals1, scratch, lscratch, info1)
       end if
       
       !: hp = W * exp(-Vabs dt/2)
       call dgemm('t','n',nstuse,nstuse,nstuse,1.d0,hp1,nstuse,exp_abp,nstuse,0.d0,scratch,nstuse)
       hp1 = scratch(1:nstuse*nstuse)
       
       call cpu_time(finish1)
       
       !: loop over intensities
       emax_loop : do iemax=1, nemax
          
          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' )  idir

          !: get emax
          emax1 = tdciresults(1+(iemax-1)*ndir)%fstrength0
          

          !: cifile binary
          If( Qci_save ) then
          funit1  = iemax*100 + idir
          cifile ='CI-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
          open( unit=funit1, file=trim(cifile), form='unformatted' )
            write(funit1) ndata, nstuse, nstates
            write(funit1) 0.d0, 0.d0, 0.d0, dirx1, diry1, dirz1, 0.d0, 0.d0, 0.d0, 1.d0
            write(funit1) real(psi0)
            write(funit1) aimag(psi0) 
            flush(funit1)      
          end if    
          
          !: RESULTS datafile
          funit2 = 1000+100*iemax + idir
          datafile = 'RESULTS-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit2,file=trim(datafile) )
          write( funit2, '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'
          write( funit2, "( '#',  20(f16.10,1x) )" ) emax1, 0.d0, &
               tdciresults(idir+(iemax-1)*ndir)%theta0, tdciresults(idir+(iemax-1)*ndir)%phi0,  0.d0,0.d0,  0.d0,0.d0, &
               dirx1, diry1, dirz1,  0.d0,0.d0,0.d0,  0.d0,0.d0,0.d0

          if( trim(jobtype).eq.flag_ip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)','6 x |psi(i)|**2'
          else
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          end if

          !: POP datafile
          funit3 = 2000+100*iemax + idir
          datafile = 'POP-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit3,file=trim(datafile) )
          if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
          else
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_a(occ)','rate_b(occ)'
          end if

          !: ION datafile
          funit4 = 3000+100*iemax + idir
          datafile = 'ION-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit4,file=trim(datafile) )
          write( funit4,"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'
 
          !: ION_COEFF datafile
          If( Qread_ion_coeff .or. Qwrite_ion_coeff ) then
            funit5 = 4000+100*iemax + idir
            datafile = 'ION_COEFF-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
            open( unit=funit5,file=trim(datafile),form='unformatted' )
          else
            funit5 = 4000+100*iemax + idir
            datafile = 'ION_COEFF-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
            open( unit=funit5,file=trim(datafile) )
            write( funit5,"(a5,2(a14,1x),3(a24,1x))" ) '#','time(fs)','rate/norm2', &
              's(i)(i=1,ip_states)','proj Zion_coeff(j,i)','Zion_coeff(j,i)'
            write(funit5,"('ntimes ',i0)") int(nstep/outstep)
            write(funit5,"('ip_states ',i0)") ip_states
          end If
          If( Qread_ion_coeff ) then 
            read(funit5) Zion_coeff
          end if

          !: MO density datafile
          If( Qmo_dens ) then
            funit6 = 5000+100*iemax + idir
            datafile = 'MO_density-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
            open( unit=funit6,file=trim(datafile) )
            write(funit6,"('ntimes ',i0)") int(nstep/outstep)
            write(funit6,"('alpha_homo ',i0)") noa
            write(funit6,"('beta_homo ',i0)")  nob
            write(funit6,"('alpha_orbitals ',i0)") noa + nva
          end if 

          !$OMP CRITICAL
          !: all thread execute the code, but only one at a time
          if(iemax.eq.1) write(iout,"(12x,'TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
          write( iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0)" ) idir, iemax, ithread
          flush( iout )
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
          end if
          If( Qread_ion_coeff ) psi = dcmplx(0.d0,0.d0)

          !: begin Looping over time
          call cpu_time( start2 )
          timestep_loop : do itime=1, nstep-1
!:             call get_norm( norm0, nstuse, psi )
!:             write(iout,"(i5,16F10.6)") itime,norm0
             
             !: modified midpoint 
             efield1 = 0.5d0 * emax1 * ( fvect1(itime) + fvect1(itime+1) )
             if( read_shift(iemax).ne.0 ) then
               If(itime.eq.1) write(iout,"(' Pulse has been shifted by ',i6,' steps')") read_shift(iemax)
               i = itime-read_shift(iemax)
               if( i.gt.0 .and. i.lt.nstep ) then
                 efield1 = 0.5d0 * emax1 * ( fvect1(i) + fvect1(i+1) )
               else
                 efield1 = 0.d0
               end if
             end if
             efield2 = 0.d0
             efieldx = dirx1 * efield1
             efieldy = diry1 * efield1
             efieldz = dirz1 * efield1
             
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
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do

             
             !: W * exp(-Vabs dt/2) * psi
             psi1 = dcmplx(0.d0,0.d0)
             do j = 1, nstuse
                jj = ( j-1) * nstuse  
                psi_j = psi(j)
                do k=1, nstuse
                   psi1(k) = psi1(k) + hp1( jj+k ) * psi_j
                end do
             end do

             
             !: exp(-E(t+dt/2)*mu*dt) * psi
             do j = 1, nstuse
                temp = dt * efield1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do

             
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W * psi
             do j = 1, nstuse
                jj = nstuse * ( j-1 )
                psi_j = dcmplx( 0.d0, 0.d0 )
                do k=1, nstuse
                   psi_j  = psi_j + hp1( jj+k ) * psi1(k)
                end do
                psi(j) = exphel(j) * psi_j
             end do
             
             If( Qwrite_ion_coeff ) then
                call get_norm( norm, nstuse, psi )
                call get_psid( nstuse, nstates, cis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,vabsmoa,vabsmob,unrestricted, &
                  rate, Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,scratch)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,scratch(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * scratch(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
             
             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
                
                call get_norm( norm, nstuse, psi )
                call get_psid( nstuse, nstates, cis_vec, norm, psi, psi_det0 )

                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,vabsmoa,vabsmob,scratch,au2fs)
                else
                  call pop_rate(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    vabsmoa,vabsmob,unrestricted,scratch,au2fs)
                end if 
                !: 1-RDM (opdm) should be stored in scratch now.

                !: Check max nva for input NOs at this step.
                if (flag_ReadU_NO) then
                  call update_maxnva( nva95max, nva99max, scratch, U_NO_input, vabsmoa )
                end if

                !$OMP CRITICAL
                !: set critical for sum so we dont have multiple threads
                !:   modifying opdm_avg at a time... is this going to destroy
                !:   performance?
                !: Add 1-RDM (opdm) to 1-RDM average for Natural Orbital generation.
                
                call add_opdm_average( opdm_avg, scratch, opdm_avg_N )
                call add_opdm_average( opdm_avg_abs, scratch, opdm_avg_N, .false., .true. )

                !call generate_natural_orbitals( scratch, U_NO, natorb_occ )
                !call NO_rate_sanity( scratch, U_NO, natorb_occ, vabsmoa, cmo_a)
                !flush(iout)    
                !call dgemm_sanity 

                !$OMP END CRITICAL
 
                call get_norm( norm,nstuse, psi )
                call get_expectation( nstuse, norm, psi, abp, rate) !: rate expectation value
                call get_expectation( nstuse, norm, psi, tdx, mux ) !: mux  expectation value
                call get_expectation( nstuse, norm, psi, tdy, muy ) !: muy  expectation value
                call get_expectation( nstuse, norm, psi, tdz, muz ) !: muz  expectation value
                !write(iout,*) " expectation rate", itime,norm,rate
                rate = -2.d0 * rate * norm**2

                if( Qci_save) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, 0.d0, &
                   dirx1, diry1, dirz1, 0.d0, 0.d0, 0.d0, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
                
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    if(norm.ne.0) then
                      write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, &
                      dble(dconjg(psi(1))*psi(1))/norm**2,dble(dconjg(psi(2))*psi(2))/norm**2, &
                      dble(dconjg(psi(3))*psi(3))/norm**2,dble(dconjg(psi(4))*psi(4))/norm**2, &
                      dble(dconjg(psi(5))*psi(5))/norm**2,dble(dconjg(psi(6))*psi(6))/norm**2, &
                      dble(dconjg(psi(7))*psi(7))/norm**2,dble(dconjg(psi(8))*psi(8))/norm**2
                   else
                      write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
                   end if
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)

                if( trim(jobtype).eq.flag_ip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),rate/au2fs,(rate_aa(i),i=1,noa*noa),&
                      (rate_ab(i),i=1,noa*nob),(rate_ba(i),i=1,nob*noa),(rate_bb(i),i=1,nob*nob)
                else
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob), &
                      (rate_a(i),i=1,noa),(rate_b(i),i=1,nob)
                end if
                flush(funit3)

                write( funit4,"( i5,f10.4,500(1x,f15.10))") &
                  idata, dble(itime)*dt*au2fs, rate/norm**2, normV, &
                  (ion(i),i=1,noa),(ion(noa+nva+i),i=1,nob)
                flush(funit4)
                
                if( Qmo_dens) then
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(scratch(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, scratch(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(scratch(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*scratch(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_cis .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,vabsmoa,vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,scratch)
                  call get_proj_ion(iout,noa+nob,ip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(ip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(scratch(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If
                If(trim(jobtype).eq.flag_ip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,vabsmoa,vabsmob, &
                    rate,Zion_coeff,ion_coeff,scratch)
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
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(scratch(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(60(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(60(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
                  flush(funit5)
                end If

             end if analysis
             
          end do timestep_loop
          call cpu_time(finish2)

          If( Qci_save ) close(funit1)
          close(funit2)
          close(funit3)
          close(funit4)

          If( Qwrite_ion_coeff ) write(funit5) Zion_coeff
          close(funit5)
 
          If( Qmo_dens ) close(funit6)

!:          scratch(1:3*(noa+nob)) = 0.d0
!:          psi(1:(noa+nob)**2) = dcmplx(0.d0,0.d0)
!:          normV = 0.d0
!:          Zion_coeff = Zion_coeff*dt
!:          do itime = 1,nstep
!:             ii = (itime-1)*(noa+nob)**2
!:             write(iout,"('coeff2',16f10.6)") 1.d+2*Zion_coeff(ii+1+noa+nob:ii+2*(noa+nob))
!:             psi(1:(noa+nob)**2) = psi(1:(noa+nob)**2) + Zion_coeff(ii+1:ii+(noa+nob)**2)
!:             psi(1+noa+nob:2*(noa+nob)) = psi(1+noa+nob:2*(noa+nob)) &
!:                 + Zion_coeff(ii+1+noa+nob:ii+2*(noa+nob))
!:             write(iout,"('psi2  ',16f10.6)") 1.d+2*psi(1+noa+nob:2*(noa+nob))
!:             scratch(1:2*(noa+nob)) = 0.d0
!:             do j = 1,noa+nob
!:               jj = ii + (j-1)*(noa+nob) 
!:               call get_norm( scratch(j),noa+nob,Zion_coeff(jj+1:jj+(noa+nob)) )
!:               call get_norm( scratch(j+noa+nob),noa+nob,psi((j-1)*(noa+nob)+1:j*(noa+nob)) )
!:               scratch(j+2*(noa+nob)) = scratch(j+2*(noa+nob))+ scratch(j)
!:             end do
!:             norm = 0.d0
!:             do j = 1,noa+nob
!:               norm = norm + scratch(j)**2
!:             end do
!:             norm = sqrt(norm)
!:             normV = normV + norm
!:             write(iout,"(i5,9f12.8)")itime,norm,(scratch(j),j=1,noa+nob)
!:             write(iout,"(' ion ',9f12.8)") normV,(scratch(j+2*(noa+nob)),j=1,noa+nob)
!:             norm = 0.d0
!:             do j = 1,noa+nob
!:               norm = norm + scratch(j+noa+nob)**2
!:             end do
!:             norm = sqrt(norm)
!:             write(iout,"(' psi ',9f12.8)") norm,(scratch(j+(noa+nob)),j=1,noa+nob)
!:          end do
 
          !$OMP CRITICAL                    
          !: record data at last timestep
          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz
          
          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2

          if (flag_ReadU_NO) then
            write(iout,"(12x,'For input NOs, nva95max=',i5,', nva99max=',i5)") nva95max, nva99max
            !: Put current nva maxs into overall nva max
            if (nva95maxmax .lt. nva95max) nva95maxmax = nva95max
            if (nva99maxmax .lt. nva99max) nva99maxmax = nva99max
          end if
          flush(iout)
          

          !$OMP END CRITICAL
       end do emax_loop
    end do dir_loop

    !$OMP END DO
    ithread = omp_get_thread_num()
    !$OMP END PARALLEL

    write(iout,*) "AVERAGED NATURAL ORBITALS START"
    write(iout,"('For input NOs, all directions: nva95maxmax=',i5,', nva99maxmax=',i5)") nva95maxmax, nva99maxmax

    call generate_natural_orbitals( opdm_avg, U_NO, natorb_occ, .false. )
    call generate_natural_orbitals( opdm_avg_abs, U_NO_abs, natorb_occ_abs, .false. )
    call NO_rate_sanity2( opdm_avg, U_NO, natorb_occ, opdm_avg_abs, U_NO_abs, natorb_occ_abs, vabsmoa, cmo_a)

    !call io_bin_test
    call write_dbin(U_NO, (noa+nva)*(noa+nva), "U_NO_out.bin")

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
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse)
    real(8) :: hp2(nstuse*nstuse), tdvals2(nstuse)
    
    !: temporary arrays to be deallocated
    real(8)    :: norm0
    real(8),allocatable :: pop0(:),pop1(:),ion(:)
    real(8),allocatable :: rate_a(:),rate_b(:),rate_aa(:),rate_ab(:),rate_ba(:),rate_bb(:)
    complex(8), allocatable :: psi_det0(:), ion_coeff(:),Zion_coeff(:)
    

    !: private variables
    integer(8) :: i, j, ii, jj, k, kk 
    integer(8) :: itime, idir, iemax
    integer(8) :: ithread, idata
    complex(8) :: cdum

    !: field info
    real(8) :: dirx1, diry1, dirz1, emax1, efield1
    real(8) :: dirx2, diry2, dirz2, emax2, efield2
    real(8) :: efieldx, efieldy, efieldz
    real(8) :: temp, temp1, temp2

    !: results and psi stuff
    real(8)    :: norm, normV, rate, mux, muy, muz
    complex(8) :: psi_j, psi_k
    complex(8) :: psi(nstuse), psi1(nstates) 

    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch, liwork
    real(8)    :: start1, start2, finish1, finish2
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: scratch(:)

    !: Natural Orbital generation
    real(8), allocatable :: opdm_avg(:)  !: Averaged one-particle reduced density matrix (1-RDM)
    real(8), allocatable :: opdm_avg_abs(:)  !: average(abs(opdm))
    real(8), allocatable :: natorb_occ(:) !: Natural orbital occupations (eigenvalues)
    real(8), allocatable :: natorb_occ_abs(:) !: Natural orbital occupations (eigenvalues)
    real(8), allocatable :: U_NO(:) !: MO->NO transformation matrix
    real(8), allocatable :: U_NO_abs(:) !: MO->NO_abs transformation matrix
    integer(8) :: opdm_avg_N

    call write_header( 'trotter_circular','propagate','enter' )    
    call writeme_propagate( 'trot_cir', 'equation' ) 
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
    call writeme_propagate( 'trot_lin', 'psi0' )    
    
    !: Natural Orbital initialization
    allocate( opdm_avg(norb*norb) )
    allocate( opdm_avg_abs(norb*norb) )
    allocate(natorb_occ(norb))
    allocate(natorb_occ_abs(norb))
    allocate( U_NO(norb*norb) )
    allocate( U_NO_abs(norb*norb) )
    opdm_avg = 0.d0
    opdm_avg_abs = 0.d0
    natorb_occ = 0.d0
    natorb_occ_abs = 0.d0
    U_NO = 0.d0
    U_NO_abs = 0.d0
    opdm_avg_N = ndir*(nstep/outstep)
    write(iout, '(A,i5)') "opdm_avg_N: ", opdm_avg_N
    
    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*ip_states*(ip_states+1)))
    else
      allocate(Zion_coeff(ip_states*(ip_states+1)))
      allocate(Zproj_ion(ip_states*ip_states))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( scratch(1+8*nstuse+2*nstuse*nstuse) )
    else
      allocate( iwork(2) )
      allocate( scratch(nstuse*nstuse) )
    end if
    
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
   
    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do
    

    !: Start loop over directions

    !$OMP PARALLEL DEFAULT(NONE),&
    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
    !$OMP cifile,datafile,dirstr,emaxstr,finish1,finish2,funit1,funit2,funit3,funit4,funit5,funit6, &
    !$OMP info1,info2,ithread,lscratch,liwork, start1,start2, &
    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, &
    !$OMP psi_j, psi_k, temp, temp1, temp2, cdum,           &
    !$OMP norm, norm0, normV, mux, muy, muz, rate, efieldx, efieldy, efieldz, &
    !$OMP pop1, ion, ion_coeff, rate_a, rate_b, rate_aa, rate_ab, rate_ba, rate_bb, psi_det0,  &
    !$OMP hp1, hp2, psi, psi1, scratch, iwork, tdvals1, tdvals2, Zion_coeff ),        &
    !$OMP SHARED( jobtype, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
    !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
    !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz, &
    !$OMP vabsmoa, vabsmob, noa, nob, nva, nvb, norb, hole_index, part_index, &
    !$OMP state_ip_index, ip_states, read_states, ip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
    !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
    !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save, &
    !$OMP opdm_avg, opdm_avg_abs, opdm_avg_N, U_NO, U_NO_abs, natorb_occ, natorb_occ_abs, cmo_a)

    
    !$OMP DO
    dir_loop : do idir=1, ndir 

       ithread = omp_get_thread_num()
       call cpu_time( start1 ) 

       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x1 ; dirx2 = tdciresults(idir)%x2
       diry1 = tdciresults(idir)%y1 ; diry2 = tdciresults(idir)%y2
       dirz1 = tdciresults(idir)%z1 ; dirz2 = tdciresults(idir)%z2

       !: Form the transition dipole in (dirx1,diry1,dirz1) and (dirx2,diry2,dirz2) directions,
       hp1(:) = dirx1*tdx(1:nstuse2) + diry1*tdy(1:nstuse2) + dirz1*tdz(1:nstuse2)
       hp2(:) = dirx2*tdx(1:nstuse2) + diry2*tdy(1:nstuse2) + dirz2*tdz(1:nstuse2)

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

          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' ) idir

          !: get emax ; emax1==emax2 in this version
          emax1 = tdciresults((iemax-1)*ndir+1)%fstrength1
          emax2 = tdciresults((iemax-1)*ndir+1)%fstrength2

          !: cifile binary
          if( Qci_save ) then
            funit1  = iemax*100 + idir
            cifile ='CI-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
            open( unit=funit1, file=trim(cifile), form='unformatted' )
            write(funit1) ndata, nstuse, nstates
            write(funit1) 0.d0, 0.d0, dirx1, diry1, dirz1, dirx2, diry2, dirz2, 1.d0
            write(funit1) real(psi0)
            write(funit1) aimag(psi0) 
            flush(funit1)      
          end if    
          
          !: RESULTS datafile
          funit2 = 1000+100*iemax + idir
          datafile = 'RESULTS-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit2,file=trim(datafile) )
          write( funit2, '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'
          write( funit2, "( '#',  20(f16.10,1x) )" ) emax1, emax2, &
               tdciresults(1+(iemax-1)*ndir)%theta0, tdciresults(1+(iemax-1)*ndir)%phi0, &
               tdciresults(1+(iemax-1)*ndir)%theta1, tdciresults(1+(iemax-1)*ndir)%phi1, &
               tdciresults(1+(iemax-1)*ndir)%theta2, tdciresults(1+(iemax-1)*ndir)%phi2, &
               tdciresults(1+(iemax-1)*ndir)%x0, tdciresults(1+(iemax-1)*ndir)%y0, tdciresults(1+(iemax-1)*ndir)%z0, & 
               dirx1, diry1, dirz1,  dirx2, diry2, dirz2
 
          if( trim(jobtype).eq.flag_ip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nav95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)','6 x |psi(i)|**2'
          else
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          end if

          !: POP datafile
          funit3 = 2000+100*iemax + idir
          datafile = 'POP-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit3,file=trim(datafile) )
          if( trim(jobtype).eq.flag_ip ) then
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
          else
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_a(occ)','rate_b(occ)'
          end if

          !: ION datafile
          funit4 = 3000+100*iemax + idir
          datafile = 'ION-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit4,file=trim(datafile) )
          write( funit4,"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'

          !: ION_COEFF datafile
          If( Qread_ion_coeff .or. Qwrite_ion_coeff ) then
            funit5 = 4000+100*iemax + idir
            datafile = 'ION_COEFF-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
            open( unit=funit5,file=trim(datafile),form='unformatted' )
          else
            funit5 = 4000+100*iemax + idir
            datafile = 'ION_COEFF-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
            open( unit=funit5,file=trim(datafile) )
            write( funit5,"(a5,2(a14,1x),3(a24,1x))" ) '#','time(fs)','rate/norm2', &
              's(i)(i=1,ip_states)','proj Zion_coeff(j,i)','Zion_coeff(j,i)'
            write(funit5,"('ntimes ',i0)") int(nstep/outstep)
            write(funit5,"('ip_states ',i0)") ip_states
          end If
          If( Qread_ion_coeff ) then
            read(funit5) Zion_coeff
          end if

          !: MO density datafile
          If( Qmo_dens ) then
            funit6 = 5000+100*iemax + idir
            datafile = 'MO_density-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
            open( unit=funit6,file=trim(datafile) )
            write(funit6,"('ntimes ',i0)") int(nstep/outstep)
            write(funit6,"('alpha_homo ',i0)") noa
            write(funit6,"('beta_homo ',i0)")  nob
            write(funit6,"('alpha_orbitals ',i0)") noa + nva
          end if 

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
             
             efield1 = 0.5d0 * emax1*(fvect1(itime)+fvect1(itime+1)) !: modified midpoint
             efield2 = 0.5d0 * emax2*(fvect2(itime)+fvect2(itime+1)) !: modified midpoint
             if( read_shift(iemax).ne.0 ) then
               If(itime.eq.1) write(iout,"(' Pulse has been shifted by ',i6,' steps')") read_shift(iemax)
               i = itime-read_shift(iemax)
               if( i.gt.0 .and. i.lt.nstep ) then
                 efield1 = 0.5d0 * emax1 * ( fvect1(i) + fvect1(i+1) )
                 efield2 = 0.5d0 * emax2 * ( fvect2(i) + fvect2(i+1) )
               else
                 efield1 = 0.d0
                 efield2 = 0.d0
               end if
             end if

             efieldx = dirx1 * efield1 + dirx2 * efield2
             efieldy = diry1 * efield1 + diry2 * efield2
             efieldz = dirz1 * efield1 + dirz2 * efield2
             temp1 = dt * efield1 * 0.5d0 
             temp2 = dt * efield2
             
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
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
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
                temp = temp2 * tdvals2(j)
                psi(j) = dcmplx(dcos(temp),dsin(temp)) * psi(j)
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
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx(dcos(temp),dsin(temp)) * psi1(j)
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
                call get_norm( norm, nstuse, psi )
                call get_psid( nstuse, nstates, cis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,vabsmoa,vabsmob,unrestricted, &
                  rate, Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,scratch)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,scratch(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * scratch(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
             
             analysis : if ( mod(itime,outstep).eq.0 ) then

                idata = int(itime/outstep)

                call get_norm( norm, nstuse, psi )
                call get_psid( nstuse, nstates, cis_vec, norm, psi, psi_det0 )
                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,vabsmoa,vabsmob,scratch,au2fs)
                else
                  call pop_rate(iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    vabsmoa,vabsmob,unrestricted,scratch,au2fs)
                end if
                !: 1-RDM should be stored in scratch now.
                !: Why does scratch have size (nstuse,nstuse) while density in
                !:  pop_rate has size (norb,norb)?

                !$OMP CRITICAL
                !: set critical for sum so we dont have multiple threads
                !:   modifying opdm_avg at a time... is this going to destroy
                !:   performance?
                !: Add 1-RDM (opdm) to 1-RDM average for Natural Orbital generation.
                
                call add_opdm_average( opdm_avg, scratch, opdm_avg_N )
                call add_opdm_average( opdm_avg_abs, scratch, opdm_avg_N, .false., .true. )

                call generate_natural_orbitals( scratch, U_NO, natorb_occ )
                call NO_rate_sanity( scratch, U_NO, natorb_occ, vabsmoa, cmo_a)
                flush(iout)    
                !call dgemm_sanity 

                !$OMP END CRITICAL

                ion_coeff = dcmplx( 0.d0,0.d0 )
                do i=1,ip_states
                  cdum = dcmplx( 0.d0,0.d0 )
                  do j=1,nstuse
!:                    cdum = cdum + proj_ion(j+(i-1)*nstuse)*psi(j)
                  end do
!:                write(iout,"(' ion_coeff_det',I4,4f12.7)") i,cdum/norm,abs(cdum/norm)
                  do k=1,ip_states
                    ion_coeff(k) = ion_coeff(k) + ip_vec(i+(k-1)*ip_states)*cdum
                  end do
                end do

                call get_norm( norm,nstuse, psi )                
                call get_expectation( nstuse, norm, psi, abp, rate) !: rate expectation value
                call get_expectation( nstuse, norm, psi, tdx, mux ) !: mux  expectation value
                call get_expectation( nstuse, norm, psi, tdy, muy ) !: muy  expectation value
                call get_expectation( nstuse, norm, psi, tdz, muz ) !: muz  expectation value
!:                write(iout,*) " expectation rate", itime,rate
                rate = -2.d0 * rate * norm**2

                if( Qci_save ) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, efield2, dirx1, diry1, dirz1, &
                    dirx2, diry2, dirz2, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
                
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    if(norm.ne.0) then
                      write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, &
                      dble(dconjg(psi(1))*psi(1))/norm**2,dble(dconjg(psi(2))*psi(2))/norm**2, &
                      dble(dconjg(psi(3))*psi(3))/norm**2,dble(dconjg(psi(4))*psi(4))/norm**2, &
                      dble(dconjg(psi(5))*psi(5))/norm**2,dble(dconjg(psi(6))*psi(6))/norm**2, &
                      dble(dconjg(psi(7))*psi(7))/norm**2,dble(dconjg(psi(8))*psi(8))/norm**2
                   else
                      write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
                   end if
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,efield2,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)

                if( trim(jobtype).eq.flag_ip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),rate/au2fs,(rate_aa(i),i=1,noa*noa),&
                      (rate_ab(i),i=1,noa*nob),(rate_ba(i),i=1,nob*noa),(rate_bb(i),i=1,nob*nob)
                else
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob), &
                      (rate_a(i),i=1,noa),(rate_b(i),i=1,nob)
                end if
                flush(funit3)

                write( funit4,"( i5,f10.4,500(1x,f15.10))") &
                  idata, dble(itime)*dt*au2fs, rate/norm**2, normV, &
                  (ion(i),i=1,noa),(ion(noa+nva+i),i=1,nob),(ion_coeff(i),i=1,ip_states)
                flush(funit4)

                if( Qmo_dens) then
                !: write real part of total density
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(scratch(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, scratch(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(scratch(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*scratch(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_cis .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,vabsmoa,vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,scratch)
                  call get_proj_ion(iout,noa+nob,ip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,scratch(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(ip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(scratch(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If
                If(trim(jobtype).eq.flag_ip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,vabsmoa,vabsmob, &
                    rate,Zion_coeff,ion_coeff,scratch)
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
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(scratch(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(60(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(60(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
                  flush(funit5)
                end If

             end if analysis

          end do timestep_loop
          call cpu_time(finish2)

          if( Qci_save ) close(funit1)
          close(funit2)
          close(funit3)
          close(funit4)
          If( Qwrite_ion_coeff ) write(funit5) Zion_coeff
          If( Qread_ion_coeff .or. Qwrite_ion_coeff ) close(funit5)

          If( Qmo_dens ) close(funit6)

          !$OMP CRITICAL
          !: record data at last timestep
          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz

          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread,idir,iemax
          write(iout,"(12x,'dir1 = (',f8.5,',',f8.5,',',f8.5,')  emax1 = ',f8.5,' au')") dirx1,diry1,dirz1,emax1
          write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") dirx2,diry2,dirz2,emax2
          write(iout,"(12x,'(iemax=1) TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
          write(iout,"(12x,'(iemax=1) LAPACK dysev TD diagonalization INFO=',i0)") info1
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2-start2
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)
          !$OMP END CRITICAL


       end do emax_loop
    end do dir_loop

    !$OMP END DO
    !$OMP END PARALLEL

    write(iout,*) "AVERAGED NATURAL ORBITALS START"

    call generate_natural_orbitals( opdm_avg, U_NO, natorb_occ, .false. )
    call generate_natural_orbitals( opdm_avg_abs, U_NO_abs, natorb_occ_abs, .false. )
    call NO_rate_sanity2( opdm_avg, U_NO, natorb_occ, opdm_avg_abs, U_NO_abs, natorb_occ_abs, vabsmoa, cmo_a)

    call write_header( 'trotter_linear','propagate','leave' )
    call cpu_time(finish)
    write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start  
    flush(iout)


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

    integer :: ndim, i, j, idx1, idx2
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

    write(iout, '(A)') 'Natural Orbital Rate Analysis (Component-wise):'
    write(iout, '(A)') '  Orb#  Occ          Vabs_NO(i,i)   Cumulative Rate'
    write(iout, '(A)') '====================================================='
    !do i = ndim,1, -1
    do i = 1,ndim
      write(iout,'(I5, A, F12.10, A, F12.10, A, F12.10)') i, ",  ", NO_occ(i), &
       ",  ", Vabs_NO((i-1)*ndim+i), ",  ", tmprate_NO(i)
    end do

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

  subroutine update_maxnva( nva95max, nva99max, opdm, U_NO, Vabs)
    implicit none
    integer(8), intent(inout) :: nva95max, nva99max
    real(8), intent(in) :: opdm(:) !: One-particle density matrix in MO basis
    real(8), intent(in) :: U_NO(:) !: Transformation matrix from MO -> NO.
    !: Interesting, vabsmo is a 2D array (:,:), so I can't do Vabs(:) here.
    !:  but I can do Vabs(norb*norb)...
    real(8), intent(in) :: Vabs(norb*norb) !: <n|Vabs|m> in MO basis

    integer(8) :: i, j, k, ndim
    integer(8) :: nva95_NO, nva99_NO, nva95_MO, nva99_MO
    real(8) :: rate_NO, rate_MO
    real(8) :: tmp_NO, tmp_MO
    real(8), allocatable :: opdm_NO(:), Vabs_NO(:)
    real(8), allocatable :: temp(:), temp2(:)

    ndim = noa + nva
    allocate( temp(ndim*ndim) )
    allocate( temp2(ndim*ndim) )
    allocate( opdm_NO(ndim*ndim) )
    allocate( Vabs_NO(ndim*ndim) )
    temp = 0.d0
    temp2 = 0.d0
    opdm_NO = 0.d0
    Vabs_NO = 0.d0
    
    !: Transform Vabs to NO basis
    !: U_NO^T * Vabs * U_NO = Vabs_NO
    call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, Vabs, ndim, 0.d0, temp, ndim)
    call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, Vabs_NO, ndim)

    !: Transform opdm to NO basis
    !: U_NO^T * opdm * U_NO = opdm_NO
    temp = 0.d0
    call dgemm('T', 'N', ndim, ndim, ndim, 1.d0, U_NO, ndim, opdm, ndim, 0.d0, temp, ndim)
    call dgemm('N', 'N', ndim, ndim, ndim, 1.d0, temp, ndim, U_NO, ndim, 0.d0, opdm_NO, ndim)

    do i=1,ndim
      do j=1,ndim
        !: 2*Tr(opdm_NO * Vabs_NO) componentwise
        rate_NO = rate_NO + 2.d0*abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j))
      end do
    end do
    do k=1,ndim
      tmp_NO = 0.d0
      do i=1,k
        do j=1,k
          tmp_NO = tmp_NO + 2.d0*abs(opdm_NO((i-1)*ndim+j) * Vabs_NO((i-1)*ndim+j))
        end do
      end do
      if (tmp_NO/rate_NO .lt. 0.95d0) nva95_NO = k
      if (tmp_NO/rate_NO .lt. 0.99d0) nva99_NO = k
    end do

    if ( nva95_NO .gt. nva95max) nva95max = nva95_NO
    if ( nva99_NO .gt. nva99max) nva99max = nva99_NO

    deallocate( temp )
    deallocate( temp2 )
    deallocate( opdm_NO )
    deallocate( Vabs_NO )

  end subroutine



end module propagate

