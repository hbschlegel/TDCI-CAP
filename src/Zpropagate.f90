module Zpropagate
  
  use variables_global
  use analysis
  use util
  
  implicit none


contains
  !==================================================================!
  !==================================================================!
  subroutine Ztrotter_linear
    
    use omp_lib    
    implicit none
    
    !: shared variables
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    real(8)    :: tdvals1(nstuse)  !: eigenvectors of hp1
    complex(8) :: hp1(nstuse*nstuse)
    complex(8), allocatable :: hp2(:), tdvals2(:)  
    
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
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: psi2(:), scratch(:)

    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch, lrwork, liwork
    real(8)    :: start1, start2, finish1, finish2
    
    
    call write_header( 'Ztrotter_linear','propagate','enter' )    
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
    

    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
!:    write(iout,*) " allocate ion_coeff",i
!:    flush(iout)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*(ip_states)*(ip_states+1)))
    else
      allocate(Zion_coeff(ip_states*(ip_states+1)))
      allocate(Zproj_ion(ip_states*(ip_states+1)))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( rwork(1+8*nstuse+2*nstuse*nstuse) )
      allocate( scratch(nstuse*nstuse+5*nstuse) )
    else
      allocate( iwork(2) )
      i = max(3*nstuse,(ip_states)*(ip_states))
      allocate( rwork(i) )
      allocate( scratch(nstuse*nstuse) )
    end if
 
!: test
          psi = dcmplx(0.d0,0.d0)
          psi(1) = dcmplx(1.d0,0.d0)
          call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),Zexp_abp,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm Vabs 0 ", norm
!: test

    norm0 = 1.d0
    
    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
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

    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do

    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
    !$OMP cifile,datafile,dirstr,emaxstr,finish1,finish2,funit1,funit2,funit3,funit4,funit5,funit6, &
    !$OMP info1,info2,ithread,lscratch,lrwork,liwork,start1,start2, &
    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, &
    !$OMP psi_j, psi_k, temp, temp1, temp2, cdum,           &
    !$OMP norm, norm0, normV, mux, muy, muz, rate, efieldx, efieldy, efieldz, &
    !$OMP pop1, ion, ion_coeff, rate_a, rate_b, rate_aa, rate_ab, rate_ba, rate_bb, psi_det0,  &
    !$OMP hp1, hp2, psi, psi1, scratch, rwork, iwork, tdvals1, tdvals2, Zion_coeff ),  &
    !$OMP SHARED( Mol, jobtype, nbasis, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
    !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
    !$OMP Zabp, Zcis_vec, Zexp_abp, exphel, fvect1, fvect2, psi0, tdciresults, Ztdx, Ztdy, Ztdz, &
    !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
    !$OMP state_ip_index, ip_states, read_states, Zip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
    !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
    !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save)
    
    !$OMP DO  

    dir_loop : do idir=1, ndir
 
       ithread = omp_get_thread_num()
       call cpu_time( start1 )
       
       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x0  
       diry1 = tdciresults(idir)%y0  
       dirz1 = tdciresults(idir)%z0  

       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*Ztdx(1:nstuse2) + diry1*Ztdy(1:nstuse2) + dirz1*Ztdz(1:nstuse2)
!: test       
          psi = dcmplx(0.d0,0.d0)
          psi(1) = dcmplx(1.d0,0.d0)
         call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),Zexp_abp,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm Vabs 0 ", norm
                flush(iout)
!: test       
      !: diagonalize mu dot E
      call testherm(nstuse,hp1,'hp1')
      If( QeigenDC ) then
         liwork = -1  
         call zheevd('v','u',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,lrwork, &
                     iwork,liwork,info1)
!: test       
                lscratch = scratch(1)
                lrwork = rwork(1)
                liwork = iwork(1)
                write(iout,*) " zheevd 0", lscratch, lrwork, liwork
                flush(iout)
!: test       
         lscratch = nstuse*nstuse + 2*nstuse
         lrwork = 1+5*nstuse+2*nstuse*nstuse
         liwork = 3+5*nstuse
                write(iout,*) " zheevd 1", lscratch, lrwork, liwork
                flush(iout)
         info1 = 10  
         call zheevd('v','u',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,lrwork, &
                     iwork,liwork,info1)
       else
         lscratch = nstuse*nstuse
         lrwork = 3*nstuse-2
         info1 = 10  
         call zheev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,info1)
       end if
       call Zfix_phase(nstuse,hp1)

!: test       
          psi = dcmplx(0.d0,0.d0)
          psi(1) = dcmplx(1.d0,0.d0)
         call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " info ", info1
                write(iout,*) " norm hp1 1 ", norm
         call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm hp1 1c ", norm
                flush(iout)
!: test       
       !: hp = W * exp(-Vabs dt/2)
       call testherm(nstuse,Zexp_abp,'Zexp_abs')
       call zgemm('t','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp1,nstuse,Zexp_abp,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
!: test       
          psi = dcmplx(0.d0,0.d0)
          psi(1) = dcmplx(1.d0,0.d0)
         call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm hp1 2 ", norm
         call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),Zexp_abp,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm Vabs 1 ", norm
         call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)             
                call get_norm( norm, nstuse, psi )
                write(iout,*) " norm hp1 Zexo_abp ", norm
          psi = dcmplx(0.d0,0.d0)
          psi(1) = dcmplx(1.d0,0.d0)
         call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),scratch,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm scratch 2 ", norm
         hp1 = scratch(1:nstuse*nstuse)
         call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm hp1 3 ", norm
                flush(iout)
!: test       
       
       call cpu_time(finish1)
       
       !: loop over intensities
       emax_loop : do iemax=1, nemax
          
          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' )  idir

          !: get emax
          emax1 = tdciresults(1+(iemax-1)*ndir)%fstrength0
!:          write(iout,*) " open unit 1"
          !: cifile binary
          if( Qci_save ) then
            funit1  = iemax*100 + idir
            cifile ='CI-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
            open( unit=funit1, file=trim(cifile), form='unformatted' )
            write(funit1) ndata, nstuse, nstates
            write(funit1) 0.d0, 0.d0, 0.d0, dirx1, diry1, dirz1, 0.d0, 0.d0, 0.d0, 1.d0
            write(funit1) real(psi0)
            write(funit1) aimag(psi0) 
            flush(funit1)      
          end if    
          
!:          write(iout,*) " open unit 2"
          !: RESULTS datafile
          funit2 = 1000+100*iemax + idir
          datafile = 'RESULTS-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit2,file=trim(datafile) )
          write( funit2, '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'
         
          write( funit2, "( '#',  20(f16.10,1x) )" ) emax1, 0.d0, &
               tdciresults(idir+(iemax-1)*ndir)%theta0, tdciresults(idir+(iemax-1)*ndir)%phi0,  0.d0,0.d0,  0.d0,0.d0, &
               dirx1, diry1, dirz1,  0.d0,0.d0,0.d0,  0.d0,0.d0,0.d0

          if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)', '6 x |psi(i)|**2'
          else
              write( funit2,"(a5,7(a10,1x),2(1x,a15),7(1x,a15) )" ) '#','time(fs)','nva95/99','field1','field2','fieldx','fieldy','fieldz', &
                'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          end if
 
!:          write(iout,*) " open unit 3"
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
 
!:          write(iout,*) " open unit 4"
          !: ION datafile
          funit4 = 3000+100*iemax + idir
          datafile = 'ION-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit4,file=trim(datafile) )
          write( funit4,"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'

!:          write(iout,*) " open unit 5"
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
            flush(funit5)
          end If
          If( Qread_ion_coeff ) then
            read(funit5) Zion_coeff
          end if

!:          write(iout,*) " open unit 6"
          !: MO density datafile
          If( Qmo_dens ) then
            funit6 = 5000+100*iemax + idir
            datafile = 'MO_density-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
            open( unit=funit6,file=trim(datafile) )
            write(funit6,"('ntimes ',i0)") int(nstep/outstep)
            write(funit6,"('alpha_homo ',i0)") noa
            write(funit6,"('beta_homo ',i0)")  nob
            write(funit6,"('alpha_orbitals ',i0)") noa + nva
            flush(funit6)
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
!: test
          write(iout,*) " start time loop"
                call get_norm( norm, nstuse, psi )
                write(iout,*) " norm 0 ", norm
          call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),Zexp_abp,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm Vabs 2 ", norm
                flush(iout)
!: test
          !: begin Looping over time
          call cpu_time( start2 )
          timestep_loop : do itime=1, nstep-1
 
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
                 call get_Zion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                      dt,Zcis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
               end If
             else
               If( itime.eq.ion_sample_start(iemax) ) psi = psi0
             end If

!: test
                call get_norm( norm, nstuse, psi )
                write(iout,*) " itime ",itime
                write(iout,*) " norm 1 ", norm
                flush(iout)
!: test
             !: exp(-iHel dt/2 ) * psi
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do
             
!: test
                call get_norm( norm, nstuse, psi )
                write(iout,*) " norm 2 ", norm
!: test
             !: W * exp(-Vabs dt/2) * psi
             psi1 = dcmplx(0.d0,0.d0)
             do j = 1, nstuse
                jj = ( j-1) * nstuse  
                psi_j = psi(j)
                do k=1, nstuse
                   psi1(k) = psi1(k) + hp1( jj+k ) * psi_j
                end do
             end do
!:             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
             
!: test
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm 3 ", norm
                flush(iout)
!: test
             !: exp(-E(t+dt/2)*mu*dt) * psi
             do j = 1, nstuse
                temp = dt * efield1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do
             
!: test
                call get_norm( norm, nstuse, psi1 )
                write(iout,*) " norm 4 ", norm
                flush(iout)
!: test
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W * psi
             do j = 1, nstuse
                jj = nstuse * ( j-1 )
                psi_j = dcmplx( 0.d0, 0.d0 )
                do k=1, nstuse
                   psi_j  = psi_j + dconjg( hp1(jj+k) ) * psi1(k)
                end do
!:                psi(j) = exphel(j) * psi_j
                psi(j) = psi_j
             end do
!:             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)             
 
!: test
                call get_norm( norm, nstuse, psi )
                write(iout,*) " norm 5 ", norm
                flush(iout)
!: test
             do j = 1, nstuse
                psi(j) = exphel(j) * psi(j)
             end do
 
!: test
                call get_norm( norm, nstuse, psi )
                write(iout,*) " norm 6 ", norm
                flush(iout)
!: test
             If( Qwrite_ion_coeff ) then
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                  rate,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,rwork)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,rwork(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * rwork(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
             
             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
   
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,Mol%vabsmoa,Mol%vabsmob,rwork,au2fs)
                else
                  call pop_rate(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    Mol%vabsmoa,Mol%vabsmob,unrestricted,rwork,au2fs)
                end if

                call get_norm( norm,nstuse, psi )                
                call get_Zexpectation( nstuse, norm, psi, Zabp, psi1, rate) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdx, psi1, mux ) !: mux  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdy, psi1, muy ) !: muy  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdz, psi1, muz ) !: muz  expectation value
!:                write(iout,*) " expect rate", itime,norm,rate 
!:                flush(iout)
                rate = -2.d0 * rate * norm**2

                if( Qci_save ) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, efield2, &
                   dirx1, diry1, dirz1, dirx2, diry2, dirz2, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
 
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, (dble(dconjg(psi(i))*psi(i)),i=1,20)
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)
                      
                if( trim(jobtype).eq.flag_socip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),(rate_aa(i),i=1,noa*noa),rate/au2fs,&
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
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(rwork(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, rwork(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(rwork(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*rwork(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_soc .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,noa+nob,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  flush(iout)
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,rwork(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(Zip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If

                If(trim(jobtype).eq.flag_socip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,ip_states,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  flush(iout)
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
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
          
          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)

          !$OMP END CRITICAL
          
       end do emax_loop
    end do dir_loop

    !$OMP END DO
    !$OMP END PARALLEL
    
    call write_header( 'Ztrotter_linear','propagate','leave' )
    call cpu_time(finish)
    write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start

    
  end subroutine Ztrotter_linear

  !==================================================================!
  !==================================================================!

  subroutine Ztrotter_circular

    ! <C> propagation using circularly polarized lights.  Takes in fvect1 and fvect2
    ! <C> C(t+dt) = exp(-iHel dt/2) * [exp(-Vabs dt/2)*W1] * exp(-iE1(t) dt/2) * [W1*W2] * exp(-iE2(t+dt/2)mu dt/2) * {W2*W1]T
    ! <C>           * exp(-iE1(t) dt/2) * [exp(-Vabs dt/2)*W1]T * exp(-iHel dt/2)*C(t)

    use omp_lib    
    implicit none
    
    !: shared variables
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    complex(8) :: hp1(nstuse*nstuse), hp2(nstuse*nstuse)
    real(8)    :: tdvals1(nstuse), tdvals2(nstuse)
    
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
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: psi2(:), scratch(:)


    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch, lrwork, liwork
    real(8)    :: start1, start2, finish1, finish2
    
    
    call write_header( 'Ztrotter_circular','propagate','enter' )    
    call writeme_propagate( 'trot_cir', 'equation' ) 
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
    

    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*(ip_states)*(ip_states+1)))
    else
      allocate(Zion_coeff((ip_states)*(ip_states+1)))
      allocate(Zproj_ion((ip_states)*(ip_states)))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( rwork(1+8*nstuse+2*nstuse*nstuse) )
      allocate( scratch(nstuse*nstuse+2*nstuse) )
    else
      allocate( iwork(2) )
      i = max(3*nstuse-2,(ip_states)*(ip_states))
      allocate( rwork(i) )
      allocate( scratch(nstuse*nstuse) )
    end if

    norm0 = 1.d0
    
    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
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
!:    write(iout,*) "pop0",norb,nrorb,pop0
!:    write(iout,*) "Vabsmo=",( Mol%vabsmoa(:,i), i=1, nrorb )

    call writeme_propagate( 'trot_lin','pop0', pop0, norb )
    
    deallocate( pop0 )

    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do
    

    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
    !$OMP cifile,datafile,dirstr,emaxstr,finish1,finish2,funit1,funit2,funit3,funit4,funit5,funit6, &
    !$OMP info1,info2,ithread,lscratch,lrwork,liwork,start1,start2, &
    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, &
    !$OMP psi_j, psi_k, temp, temp1, temp2, cdum,           &
    !$OMP norm, norm0, normV, mux, muy, muz, rate, efieldx, efieldy, efieldz, &
    !$OMP pop1, ion, ion_coeff, rate_a, rate_b, rate_aa, rate_ab, rate_ba, rate_bb, psi_det0,  &
    !$OMP hp1, hp2, psi, psi1, scratch, rwork, iwork, tdvals1, tdvals2, Zion_coeff ),  &
    !$OMP SHARED( Mol, jobtype, nbasis, flag_cis, flag_tda, flag_ip, flag_soc, flag_socip, &
    !$OMP au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
    !$OMP Zabp, Zcis_vec, Zexp_abp, exphel, fvect1, fvect2, psi0, tdciresults, Ztdx, Ztdy, Ztdz, &
    !$OMP noa, nob, nva, nvb, norb, hole_index, part_index, &
    !$OMP state_ip_index, ip_states, read_states, Zip_vec, Zproj_ion, Qread_ion_coeff, Qwrite_ion_coeff, &
    !$OMP read_state1, read_state2, read_coeff1, read_coeff2, read_shift, &
    !$OMP ion_sample_start, ion_sample_width, ion_sample_state, unrestricted, QeigenDC, Qmo_dens, Qci_save)

    !$OMP DO  

    dir_loop : do idir=1, ndir
       
       ithread = omp_get_thread_num()
       call cpu_time( start1 )
       
       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x1 ; dirx2 = tdciresults(idir)%x2
       diry1 = tdciresults(idir)%y1 ; diry2 = tdciresults(idir)%y2
       dirz1 = tdciresults(idir)%z1 ; dirz2 = tdciresults(idir)%z2
       
       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*Ztdx(1:nstuse2) + diry1*Ztdy(1:nstuse2) + dirz1*Ztdz(1:nstuse2)
       hp2(:) = dirx2*Ztdx(1:nstuse2) + diry2*Ztdy(1:nstuse2) + dirz2*Ztdz(1:nstuse2)

       !: diagonalize mu dot E
       info1 = 10  
       info2 = 10  
       If( QeigenDC ) then
          lscratch = nstuse*nstuse + 2*nstuse
          lrwork = 1+5*nstuse+2*nstuse*nstuse
          liwork = 3+5*nstuse
          call zheevd('v','u',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,lrwork, &
                      iwork,liwork,info1)
          call zheevd('v','u',nstuse,hp2,nstuse,tdvals2,scratch,lscratch,rwork,lrwork, &
                      iwork,liwork,info2)
        else
          lscratch = nstuse*nstuse
          lrwork = 3*nstuse-2
          call zheev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,info1)
          call zheev('vectors','upper',nstuse,hp2,nstuse,tdvals2,scratch,lscratch,rwork,info2)
        end if
        call Zfix_phase(nstuse,hp1)
        call Zfix_phase(nstuse,hp2)

       !: hp2 = hp2 * hp1 ; hp2 = W2 * W1
       call zgemm('c','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp2,nstuse,hp1,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
       hp2 = scratch(1:nstuse*nstuse)
       !: hp1 = W1 * exp(-Vabs dt/2) 
       call zgemm('c','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp1,nstuse,Zexp_abp,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
       hp1 = scratch(1:nstuse*nstuse)
       
       call cpu_time(finish1)
       
       !: loop over intensities
       emax_loop : do iemax=1, nemax
          
          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' )  idir

          !: get emax
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
         
          if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)', '6 x |psi(i)|**2'
          else
              write( funit2,"(a5,7(a10,1x),2(1x,a15),7(1x,a15) )" ) '#','time(fs)','nva95/99','field1','field2','fieldx','fieldy','fieldz', &
                'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          end if
          flush(funit2)
 
          !: POP datafile
          funit3 = 2000+100*iemax + idir
          datafile = 'POP-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit3,file=trim(datafile) )
          if( trim(jobtype).eq.flag_socip) then
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
          else
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_a(occ)','rate_b(occ)'
          end if
          flush(funit3)
 
          !: ION datafile
          funit4 = 3000+100*iemax + idir
          datafile = 'ION-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit4,file=trim(datafile) )
          write( funit4,"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'
          flush(funit4)

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

             !: modified midpoint
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
                 call get_Zion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                      dt,Zcis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
               end If
               If( itime.eq.ion_sample_start(iemax) ) psi = psi0
             end If

             !: exp(-iHel dt/2 ) * psi
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do
             
             !: W1 * exp(-Vabs dt/2) * psi
!:             psi1 = dcmplx(0.d0,0.d0)
!:             do j = 1, nstuse
!:                jj = ( j-1) * nstuse  
!:                psi_j = psi(j)
!:                do k=1, nstuse
!:                   psi1(k) = psi1(k) + hp1( jj+k ) * psi_j
!:                end do
!:             end do
             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)  

             !: exp( iE1(t+dt/2)*mu * dt/2) * psi
             do j = 1, nstuse
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do
             
             !: W2*W1 * psi 
!:             psi = dcmplx( 0.d0, 0.d0 )
!:             do j = 1, nstuse
!:                jj = (j-1)*nstuse
!:                psi_j  = psi1(j)
!:                do k=1, nstuse
!:                   psi(k) = psi(k) + hp2(jj+k) * psi_j
!:                end do
!:             end do
             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp2,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)  
             
             !: exp( iE2(t+dt/2)*mu* dt ) * psi
             do j = 1, nstuse
                temp = temp2 * tdvals2(j)
                psi(j) = dcmplx(dcos(temp),dsin(temp)) * psi(j)
             end do
 
             !: W1*W2 * psi
!:             do j = 1, nstuse
!:                jj = nstuse*(j-1)
!:                psi_j = dcmplx( 0.d0, 0.d0 )
!:                do k=1, nstuse 
!:                   psi_j = psi_j + dconjg( hp2(jj+k) ) * psi(k)
!:                end do
!:                psi1(j) = psi_j
!:             end do
             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp2,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)

             !: exp( iE1(t+dt/2)*mu*dt/2) * psi
             do j = 1, nstuse
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx(dcos(temp),dsin(temp)) * psi1(j)
             end do
             
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W1 * psi
!:             do j = 1, nstuse
!:                jj = nstuse * ( j-1 )
!:                psi_j = dcmplx( 0.d0, 0.d0 )
!:                do k=1, nstuse
!:                   psi_j  = psi_j + dconjg( hp1(jj+k) ) * psi1(k)
!:                end do
!:                psi(j) = exphel(j) * psi_j
!:             end do
             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)
             
             do j = 1, nstuse
                psi(j) = exphel(j) * psi(j)
             end do

             If( Qwrite_ion_coeff ) then
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                  rate,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,rwork)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,rwork(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * rwork(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
            
             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
        
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,Mol%vabsmoa,Mol%vabsmob,rwork,au2fs)
                else
                  call pop_rate(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    Mol%vabsmoa,Mol%vabsmob,unrestricted,rwork,au2fs)
                end if
        
                call get_norm( norm,nstuse, psi )                
                call get_Zexpectation( nstuse, norm, psi, Zabp, psi1, rate) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdx, psi1, mux ) !: mux  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdy, psi1, muy ) !: muy  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdz, psi1, muz ) !: muz  expectation value
!:                write(iout,*) " expect rate", itime,norm,rate 
                rate = -2.d0 * rate * norm**2

                if( Qci_save ) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, efield2, &
                   dirx1, diry1, dirz1, dirx2, diry2, dirz2, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
 
                if( trim(jobtype).eq.flag_socip) then
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,efield2,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, (dble(dconjg(psi(i))*psi(i)),i=1,20)
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,efield2,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)
                      
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),(rate_aa(i),i=1,noa*noa),rate/au2fs,&
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
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(rwork(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, rwork(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(rwork(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*rwork(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_soc .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,noa+nob,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,rwork(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(Zip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If

                If(trim(jobtype).eq.flag_socip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,ip_states,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  flush(iout)
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
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
          
          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") dirx2,diry2,dirz2,emax2
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)

          !$OMP END CRITICAL
          
       end do emax_loop
    end do dir_loop

    !$OMP END DO
    !$OMP END PARALLEL
    
    call write_header( 'Ztrotter_circular','propagate','leave' )
    call cpu_time(finish)
    write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start

    
  end subroutine Ztrotter_circular

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
  !==================================================================!
  !==================================================================!
  subroutine Ztrott_serial_lin
    
    implicit none
    
    !: shared variables
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    real(8)    :: tdvals1(nstuse)  !: eigenvectors of hp1
    complex(8) :: hp1(nstuse*nstuse)
    complex(8), allocatable :: hp2(:), tdvals2(:)  
    
    !: temporary arrays to be deallocated
    real(8) :: norm0
    real(8),allocatable :: pop0(:),pop1(:),ion(:)
    real(8),allocatable :: rate_a(:),rate_b(:),rate_aa(:),rate_ab(:),rate_ba(:),rate_bb(:)
    complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)
    

    !: private variables
    integer(8) :: i, j, ii, jj, k, kk
    integer(8) :: itime, idir, iemax
    integer(8) :: idata
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
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: psi2(:), scratch(:)

    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch, lrwork, liwork
    real(8)    :: start1, start2, finish1, finish2
    
    
    call write_header( 'Ztrotter_linear','propagate','enter' )    
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
    

    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*(noa+nob)*(noa+nob+1)))
    else
      allocate(Zion_coeff((noa+nob)**2))
      allocate(Zproj_ion((noa+nob)*(noa+nob)))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( rwork(1+8*nstuse+2*nstuse*nstuse) )
      allocate( scratch(nstuse*nstuse+2*nstuse) )
    else
      allocate( iwork(2) )
      i = max(3*nstuse-2,(noa*nva)*(noa+nva))
      allocate( rwork(i) )
      allocate( scratch(nstuse*nstuse) )
    end if

    norm0 = 1.d0
    
    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
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

    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do

    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    dir_loop : do idir=1, ndir
       
       call cpu_time( start1 )
       

       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x0  
       diry1 = tdciresults(idir)%y0  
       dirz1 = tdciresults(idir)%z0  
       

       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*Ztdx(1:nstuse2) + diry1*Ztdy(1:nstuse2) + dirz1*Ztdz(1:nstuse2)
       
       !: diagonalize mu dot E
       info1 = 10  
      If( QeigenDC ) then
         lscratch = nstuse*nstuse + 2*nstuse
         lrwork = 1+5*nstuse+2*nstuse*nstuse
         liwork = 3+5*nstuse
         call zheevd('v','u',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,lrwork, &
                     iwork,liwork,info1)
       else
         lscratch = nstuse*nstuse
         lrwork = 3*nstuse-2
         call zheev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,info1)
       end if
       call Zfix_phase(nstuse,hp1)
       
       !: hp = W * exp(-Vabs dt/2)
!:       call Zmatmult( hp1, Zexp_abp(1:nstuse2), nstuse )
       call zgemm('c','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp1,nstuse,Zexp_abp,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
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
          if( Qci_save ) then
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

          if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95,99','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)', '6 x |psi(i)|**2'
          else
            write( funit2,"(a5,7(a10,1x),2(1x,a15),7(1x,a15) )" ) '#','time(fs)','nva95/99','field1','field2','fieldx','fieldy','fieldz', &
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

          if(iemax.eq.1) write(iout,"(12x,'TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
          write( iout,"(' start propagation for direction',i4,' intensity',i4)" ) idir, iemax
          flush( iout )
                    
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
                 call get_Zion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                      dt,Zcis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
               end If
               If( itime.eq.ion_sample_start(iemax) ) psi = psi0
             end If

             !: exp(-iHel dt/2 ) * psi
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do
             
             !: W * exp(-Vabs dt/2) * psi
             psi1 = dcmplx(0.d0,0.d0)
             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)             
             
             !: exp(-E(t+dt/2)*mu*dt) * psi
             do j = 1, nstuse
                temp = dt * efield1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do
             
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W * psi
             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)             
             do j = 1, nstuse
                psi(j) = exphel(j) * psi(j)
             end do
             
             If( Qwrite_ion_coeff ) then
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                  rate,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,rwork)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,rwork(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * rwork(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
            
             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
   
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,Mol%vabsmoa,Mol%vabsmob,rwork,au2fs)
                else
                  call pop_rate(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    Mol%vabsmoa,Mol%vabsmob,unrestricted,rwork,au2fs)
                end if

                call get_norm( norm,nstuse, psi )                
                call get_Zexpectation( nstuse, norm, psi, Zabp, psi1, rate) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdx, psi1, mux ) !: mux  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdy, psi1, muy ) !: muy  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdz, psi1, muz ) !: muz  expectation value
!:                write(iout,*) " expectation rate", itime,rate 
                rate = -2.d0 * rate * norm**2

                if( Qci_save ) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, efield2, &
                   dirx1, diry1, dirz1, dirx2, diry2, dirz2, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
 
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, (dble(dconjg(psi(i))*psi(i)),i=1,20)
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,0.d0,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)
                      
                if( trim(jobtype).eq.flag_socip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),(rate_aa(i),i=1,noa*noa),rate/au2fs,&
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
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(rwork(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, rwork(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(rwork(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*rwork(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_soc .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,noa+nob,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,rwork(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(Zip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If

                If(trim(jobtype).eq.flag_socip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,ip_states,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  flush(iout)
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
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

          !: record data at last timestep
          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz
          
          write(iout,"(' propagation done for direction',i4,' and intensity',i4)") idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)

       end do emax_loop
    end do dir_loop

    call write_header( 'Ztrotter_linear','propagate','leave' )
    call cpu_time(finish)
    write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start

    
  end subroutine Ztrott_serial_lin

  !==================================================================!
  !==================================================================!

  subroutine Ztrott_serial_cir

    ! <C> propagation using circularly polarized lights.  Takes in fvect1 and fvect2
    ! <C> C(t+dt) = exp(-iHel dt/2) * [exp(-Vabs dt/2)*W1] * exp(-iE1(t) dt/2) * [W1*W2] * exp(-iE2(t+dt/2)mu dt/2) * {W2*W1]T
    ! <C>           * exp(-iE1(t) dt/2) * [exp(-Vabs dt/2)*W1]T * exp(-iHel dt/2)*C(t)

    implicit none
    
    !: shared variables
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    complex(8) :: hp1(nstuse*nstuse), hp2(nstuse*nstuse)
    real(8)    :: tdvals1(nstuse), tdvals2(nstuse)
    
    !: temporary arrays to be deallocated
    real(8) :: norm0
    real(8),allocatable :: pop0(:),pop1(:),ion(:)
    real(8),allocatable :: rate_a(:),rate_b(:),rate_aa(:),rate_ab(:),rate_ba(:),rate_bb(:)
    complex(8), allocatable :: psi_det0(:),ion_coeff(:),Zion_coeff(:)
    

    !: private variables
    integer(8) :: i, j, ii, jj, k, kk
    integer(8) :: itime, idir, iemax
    integer(8) :: idata
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
    integer(8), allocatable :: iwork(:)
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: psi2(:), scratch(:)


    !: file stuff
    integer(8)   :: funit1, funit2, funit3, funit4, funit5, funit6
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch, lrwork, liwork
    real(8)    :: start1, start2, finish1, finish2
    
    
    call write_header( 'Ztrotter_circular','propagate','enter' )    
    call writeme_propagate( 'trot_cir', 'equation' ) 
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
    

    !: get initial population
    allocate( pop0(norb), pop1(norb), ion(norb), psi_det0(nstates) )
    allocate( rate_a(noa), rate_b(nob) )
    allocate( rate_aa(noa*noa), rate_ab(noa*nob), rate_ba(nob*noa), rate_bb(nob*nob) )
    i = max(2*ip_states*(nva+nvb),ip_states*ip_states)
    allocate(ion_coeff(i))
    If( Qwrite_ion_coeff .or. Qread_ion_coeff ) then
      allocate(Zion_coeff(nstep*(noa+nob)*(noa+nob+1)))
    else
      allocate(Zion_coeff((noa+nob)**2))
      allocate(Zproj_ion((noa+nob)*(noa+nob)))
    end If
    If( QeigenDC ) then
      allocate( iwork(3+5*nstuse) )
      allocate( rwork(1+8*nstuse+2*nstuse*nstuse) )
      allocate( scratch(nstuse*nstuse+2*nstuse) )
    else
      allocate( iwork(2) )
      i = max(3*nstuse-2,(noa*nva)*(noa+nva))
      allocate( rwork(i) )
      allocate( scratch(nstuse*nstuse) )
    end if

    norm0 = 1.d0
    
    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
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
!:    write(iout,*) "pop0",norb,nrorb,pop0
!:    write(iout,*) "Vabsmo=",( Mol%vabsmoa(:,i), i=1, nrorb )

    call writeme_propagate( 'trot_lin','pop0', pop0, norb )
    
    deallocate( pop0 )

    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do
    

    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    dir_loop : do idir=1, ndir
       
       call cpu_time( start1 )
       
       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x1 ; dirx2 = tdciresults(idir)%x2
       diry1 = tdciresults(idir)%y1 ; diry2 = tdciresults(idir)%y2
       dirz1 = tdciresults(idir)%z1 ; dirz2 = tdciresults(idir)%z2
       
       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*Ztdx(1:nstuse2) + diry1*Ztdy(1:nstuse2) + dirz1*Ztdz(1:nstuse2)
       hp2(:) = dirx2*Ztdx(1:nstuse2) + diry2*Ztdy(1:nstuse2) + dirz2*Ztdz(1:nstuse2)

       !: diagonalize mu dot E
       info1 = 10  
       info2 = 10  
       If( QeigenDC ) then
          lscratch = nstuse*nstuse + 2*nstuse
          lrwork = 1+5*nstuse+2*nstuse*nstuse
          liwork = 3+5*nstuse
          call zheevd('v','u',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,lrwork, &
                      iwork,liwork,info1)
          call zheevd('v','u',nstuse,hp2,nstuse,tdvals2,scratch,lscratch,rwork,lrwork, &
                      iwork,liwork,info2)
        else
          lscratch = nstuse*nstuse
          lrwork = 3*nstuse-2
          call zheev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,info1)
          call zheev('vectors','upper',nstuse,hp2,nstuse,tdvals2,scratch,lscratch,rwork,info2)
        end if
        call Zfix_phase(nstuse,hp1)
        call Zfix_phase(nstuse,hp2)

       !: hp2 = hp2 * hp1 ; hp2 = W2 * W1
       call zgemm('c','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp2,nstuse,hp1,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
       hp2 = scratch(1:nstuse*nstuse)
       !: hp1 = W1 * exp(-Vabs dt/2) 
       call zgemm('c','n',nstuse,nstuse,nstuse,dcmplx(1.d0,0.d0), &
         hp1,nstuse,Zexp_abp,nstuse,dcmplx(0.d0,0.d0),scratch,nstuse)
       hp1 = scratch(1:nstuse*nstuse)
       
       call cpu_time(finish1)
       
       !: loop over intensities
       emax_loop : do iemax=1, nemax
          
          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' )  idir

          !: get emax
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
         
          if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
            write( funit2,"(a5,7(a10,1x),2(1x,a15),10(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
              'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)', '6 x |psi(i)|**2'
          else
              write( funit2,"(a5,7(a10,1x),2(1x,a15),7(1x,a15) )" ) '#','time(fs)','nva95/99 ','field1','field2','fieldx','fieldy','fieldz', &
                'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          end if
          flush(funit2)
 
          !: POP datafile
          funit3 = 2000+100*iemax + idir
          datafile = 'POP-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit3,file=trim(datafile) )
          if( trim(jobtype).eq.flag_socip) then
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_aa(occ)','rate_ab(occ)','rate_ba(occ)','rate_bb(occ)'
          else
            write( funit3,"(a5,8(a14,1x))" ) '#','time(fs)','norm2','pop_a(occ)','pop_b(occ)', &
              'rate(fs-1)','rate_a(occ)','rate_b(occ)'
          end if
          flush(funit3)
 
          !: ION datafile
          funit4 = 3000+100*iemax + idir
          datafile = 'ION-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.dat'
          open( unit=funit4,file=trim(datafile) )
          write( funit4,"(a5,8(a14,1x))" ) '#','time(fs)','rate/norm2','normV','ion_a(occ)','ion_b(occ)','ion_coeff'
          flush(funit4)

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

          if(iemax.eq.1) write(iout,"(12x,'TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
          write( iout,"(' start propagation for direction',i4,' intensity',i4)" ) idir, iemax
          flush( iout )
                    
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

             !: modified midpoint
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
                 call get_Zion_psi1(iout,nstates,nstuse,noa+nob,ion_sample_state(iemax), &
                      dt,Zcis_vec,psi,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),psi_det0)
               end If
               If( itime.eq.ion_sample_start(iemax) ) psi = psi0
             end If

             !: exp(-iHel dt/2 ) * psi
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do
             
             !: W1 * exp(-Vabs dt/2) * psi
             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)  

             !: exp( iE1(t+dt/2)*mu * dt/2) * psi
             do j = 1, nstuse
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do
             
             !: W2*W1 * psi 
             call zgemv('n',nstuse,nstuse,dcmplx(1.d0,0.d0),hp2,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)  
             
             !: exp( iE2(t+dt/2)*mu* dt ) * psi
             do j = 1, nstuse
                temp = temp2 * tdvals2(j)
                psi(j) = dcmplx(dcos(temp),dsin(temp)) * psi(j)
             end do
 
             !: W1*W2 * psi
             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp2,nstuse,psi,1,dcmplx(0.d0,0.d0),psi1,1)

             !: exp( iE1(t+dt/2)*mu*dt/2) * psi
             do j = 1, nstuse
                temp = temp1 * tdvals1(j)
                psi1(j) = dcmplx(dcos(temp),dsin(temp)) * psi1(j)
             end do
             
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W1 * psi
             call zgemv('c',nstuse,nstuse,dcmplx(1.d0,0.d0),hp1,nstuse,psi1,1,dcmplx(0.d0,0.d0),psi,1)
             
             do j = 1, nstuse
                psi(j) = exphel(j) * psi(j)
             end do

             If( Qwrite_ion_coeff ) then
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                ii = (itime-1)*(noa+nob)*(noa+nob+1)
                call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                  psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                  rate,Zion_coeff(ii+1:ii+(noa+nob)*(noa+nob+1)),ion_coeff,rwork)
!:                  write(iout,"('W0',i5,i7,16f12.8)") itime,ii,rate,rwork(1:noa+nob)
                  do i = 1,noa+nob
                     do j =1,noa+nob
                       Zion_coeff(ii+i*(noa+nob)+j) =  &
                         dsqrt(rate) * rwork(i) * Zion_coeff(ii+i*(noa+nob)+j)
                     end do
                  end do
!:                  if ( mod(itime,outstep).eq.0 )  write(iout,"(i5,' SVD ',10(5F12.8/))") &
!:                                                  itime,norm,(scratch(i),i=1,noa+nob)               
             end If
            
             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
        
                call get_norm( norm, nstuse, psi )
                call get_Zpsid( nstuse, nstates, Zcis_vec, norm, psi, psi_det0 )
                if ( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip ) then
                  call pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb, &
                    psi_det0,psi1,normV,Mol%vabsmoa,Mol%vabsmob,rwork,au2fs)
                else
                  call pop_rate(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,jj,kk, &
                    hole_index,part_index,state_ip_index,ip_states, &
                    pop1,ion,ion_coeff,rate_a,rate_b,psi_det0,psi1,normV, &
                    Mol%vabsmoa,Mol%vabsmob,unrestricted,rwork,au2fs)
                end if
        
                call get_norm( norm,nstuse, psi )                
                call get_Zexpectation( nstuse, norm, psi, Zabp, psi1, rate) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdx, psi1, mux ) !: mux  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdy, psi1, muy ) !: muy  expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdz, psi1, muz ) !: muz  expectation value
!:                write(iout,*) " expectation rate", itime,rate 
                rate = -2.d0 * rate * norm**2

                if( Qci_save) then
                  write(funit1) dble(itime)*dt*au2fs, efield1, efield2, &
                   dirx1, diry1, dirz1, dirx2, diry2, dirz2, norm**2
                  write(funit1) real( psi )
                  write(funit1) aimag( psi )
                  flush(funit1)
                end if
 
                if( trim(jobtype).eq.flag_socip) then
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,efield2,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz, (dble(dconjg(psi(i))*psi(i)),i=1,20)
                else
                    write( funit2,"( i5,f10.4,2i5,1x,5(f10.7,1x),2(1x,f15.10),500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs,jj,kk,efield1,efield2,efieldx,efieldy,efieldz, &
                      norm**2, rate/au2fs, mux, muy, muz
                end if
                flush(funit2)
                      
                if( trim(jobtype).eq.flag_ip .or. trim(jobtype).eq.flag_socip) then
                    write( funit3,"( i5,f10.4,500(1x,f15.10))") &
                      idata, dble(itime)*dt*au2fs, norm**2, &
                      (pop1(i),i=1,noa),(pop1(noa+nva+i),i=1,nob),(rate_aa(i),i=1,noa*noa),rate/au2fs,&
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
                  write(funit6,"(f13.9,f16.10)") dble(itime)*dt*au2fs,rate
                  do i=1, noa+nva
                    if ( abs(rwork(i+(i-1)*(noa+nva))).gt.1.d-6 ) &
!:                      write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                      write(funit6,"(i5,i5,1x,f13.10)") &
                        i, i, rwork(i+(i-1)*(noa+nva))
                    do j=(i+1), noa+nva
                      if ( abs(rwork(i+(j-1)*(noa+nva))).gt.1.d-6 ) &
!:                        write(funit6,"(i5,i5,1x,f13.10,',')",advance='no') &
                        write(funit6,"(i5,i5,1x,f13.10)") &
                          i, j, 2.d0*rwork(i+(j-1)*(noa+nva))
                      end do
                    end do
                    write(funit6,"(i5,i5,1x,f13.10)") 0,0,0.d0
                end if

                If(trim(jobtype).eq.flag_soc .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob,unrestricted, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,noa+nob,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  do i=1,noa+nob
!:                    write(iout,"('SVD s',i8,16f13.7)") idata,rwork(i)
!:                    write(iout,"('coeff',i3,16f13.7)") i,(Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('ipvec',i3,16f13.7)") i,(Zip_vec(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                    write(iout,"('proj ',i3,16f13.7)") i,(Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
!:                  end do
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,noa+nob)
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*(noa+nob)),j=1,noa+nob)
                  end do
                  do i = 1,noa+nob
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*(noa+nob)),j=1,noa+nob)
                  end do
                end If

                If(trim(jobtype).eq.flag_socip .and. (.not. Qwrite_ion_coeff) ) then
                  call get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole_index,part_index, &
                    ip_states,state_ip_index,psi_det0,psi1,norm,Mol%vabsmoa,Mol%vabsmob, &
                    rate,Zion_coeff,ion_coeff,rwork)
                  call get_Zproj_ion(iout,ip_states,Zip_vec,Zion_coeff(noa+nob+1),Zproj_ion)
!:                  write(iout,*) " ip_states",ip_states
!:                  flush(iout)
                  write( funit5,"(i5,f10.4,50(1x,f15.10))") &
                    idata, dble(itime)*dt*au2fs, rate,(rwork(i),i=1,ip_states)
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zproj_ion(j+(i-1)*ip_states),j=1,ip_states)
                  end do
                  do i = 1,ip_states
                    write( funit5,"(50(1x,f15.10))") (Zion_coeff(j+i*ip_states),j=1,ip_states)
                  end do
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
 
          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz
          
          write(iout,"(' propagation done for direction',i4,' and intensity',i4)") idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") dirx2,diry2,dirz2,emax2
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)

       end do emax_loop
    end do dir_loop

    call write_header( 'Ztrotter_circular','propagate','leave' )
    call cpu_time(finish)
    write(iout,"(2x,'Total propagation time:',f12.4,' s')") finish - start

    
  end subroutine Ztrott_serial_cir

end module Zpropagate
