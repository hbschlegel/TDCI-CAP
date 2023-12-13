module mod_davidson_cisd

  
  !: Davidson-Liu diagonalization for CISD-IP
  !: REF:  M.L. Leininger, C.D. Sherrill, W.D. Allen, H.F. Schaeffer III
  !:       J. Comp. Chem. 22, 1574-1589 (2001)
  !:       Systematic Study of Selected Diagonalization Methods for Configuration Interaction Matrices
  !: REF:  E.R. Davidson and W.J. Thompson, Computers in Physics, 7, 519-522 (1993)
  !:       Monster Matrices:  Their Eigenvalues and Eigenvectors
  
  use omp_lib
  use global_variables
  use util

  implicit none
  
  
contains
  !--------------------------------------------------------------------------------------! 
  !--------------------------------------------------------------------------------------!
  subroutine davidson_cisd
    
    use get_ham0_cisd

    implicit none
    
    integer(8), parameter :: mroot = 5
    real(8),    parameter :: energy_tol    = 1.d-9
    real(8),    parameter :: residue_tol   = 1.d-7
    real(8),    parameter :: expand_cutoff = 0.01d0    
    
    integer(8) :: ibasis, jbasis, iroot, ia
    integer(8) :: lbasis0, lbasis, max_basis, start_basis
    integer(8) :: i, j, a, b, x
    
    real(8) :: alpha_k, rho_iroot, norm2, projection
    real(8) :: rtmp, rtmpx, rtmpy, rtmpz
    logical :: collapse
    
    !: big arrays
    real(8) :: diag(nstates), delta(nstates), residue(nstates), tmpmat(nstates)
    real(8), allocatable :: bigG(:,:), b_i(:,:), sigma(:,:), eig(:)
    
    !: convergence criteria
    integer(8) :: iiter
    real(8)    :: max_endiff, max_resdiff, energy_diff(mroot), res_diff(mroot), store_eig(mroot)

    !: dipole 
    real(8)    :: tdipx(mroot), tdipy(mroot), tdipz(mroot)
    
    !: LAPACK variables
    integer(8) :: lscratch, info, lwork
    real(8), allocatable :: work(:)
    
    !: timestamp
    real(8) :: lstart, lfinish
    
    
    
    call write_header( 'davidson_ucisdip', 'davidson', 'enter' )
    write(iout,"(' convergence on energy :  ',e10.4,' au ',e10.4,' eV ')") energy_tol, energy_tol*au2eV
    write(iout,"(' convergence on residue:  ',e10.4)") residue_tol
    write(iout,"(A)") '' 
    
    
    call cpu_time(lstart)           
    
    
    !: allocate guess space
    max_basis = 10*mroot
    allocate( b_i(nstates,max_basis) )    ; b_i   = 0.d0  
    allocate( sigma(nstates,max_basis) )  ; sigma = 0.d0 
    allocate( bigG(max_basis,max_basis) ) ; bigG  = 0.d0 
    allocate( eig(max_basis) )            ; eig   = 0.d0
    
    
    !: get guess vectors, all single excitation from occupied AND two lowest double excitations
    b_i = 0.d0
    !: ground state
    b_i(1) = 1.d0  
    !: singly excited states
    ia = state_index(-noa,0)     ; b_i(ia) = 1.d0 
    ia = state_index(-(noa-1),0) ; b_i(ia) = 1.d0
    ia = state_index( nob,0)     ; b_i(ia) = 1.d0
    ia = state_index( (nob-1),0) ; b_i(ia) = 1.d0
    !: doubly excited states
    ia = state_index(-(noa-1),-noa) ; b_i(ia) = 1.d0
    ia = state_index(-(noa-2),-noa) ; b_i(ia) = 1.d0
    ia = state_index( (nob-1), nob) ; b_i(ia) = 1.d0
    ia = state_index( (nob-2), nob) ; b_i(ia) = 1.d0

    
    !: get diag
    !$OMP PARALLEL DEFAULT(NONE), PRIVATE( ia, tmpmat ), SHARED( nstates, diag )
    !$OMP DO
    do ia=1, nstates
       call get_h_column(ia,tmpmat)
       diag(ia) = tmpmat(ia)
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    write(iout,"(' ')") 
    write(iout,"(' NSTATES = ',i0)") nstates
    write(iout,"(' NUMBER OF ROOTS to seek = ',i0)") mroot
    write(iout,"(' MAX vector dimension    = ',i0)") max_basis
    write(iout,"(' INITIAL guess dimension = ',i0)") lbasis
    write(iout,"(A)") ' '
    flush(iout)


    start_basis = 1    
    iiter    = 0   
    collapse = .false.
    max_resdiff = 10.d0 ; max_endiff = 10.d0 ; store_eig=0.d0 
    
    
    !: ITERATEITERATEITERATEITERATEITERATEITERATEITERATEITERATEITERATE
    converge : do while ( max_resdiff.ge.residue_tol .and. max_endiff.ge.energy_tol )        
       
       lbasis0 = lbasis 
       iiter   = iiter + 1        
       
       !: G_ji = < b_j | H | b_i > = < b_j | sigma_i >
       get_sigma : do ibasis=start_basis, lbasis
          
          !: use residue as temporary vector for OMP
          !$OMP PARALLEL DEFAULT(NONE), PRIVATE(ia, i,a,j,b, rtmp, tmpmat, residue), SHARED(b_i, ibasis, nstates, sigma, state_index, noa, nva, nob, nvb)
          residue = 0.d0

          !: singles alpha
          !$OMP DO
          do i=1, noa
             do a=1, nva                
                ia = state_index(-i,0) + a - 1 
                rtmp = b_i(ia,ibasis)
                if ( rtmp.ne.0.d0 ) then
                   call get_cisd_ia_AA(i,a,tmpmat)
                   residue(:) = residue(:) + rtmp*tmpmat(:)
                end if
             end do
          end do
          !$OMP END DO
          !$OMP CRITICAL
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL

          !: singles beta
          !$OMP DO
          do i=1, nob
             do a=1, nvb                
                ia = state_index(i,0) + a - 1 
                rtmp = b_i(ia,ibasis)
                if ( rtmp.ne.0.d0 ) then
                   call get_cisd_ia_BB(i,a,tmpmat)
                   residue(:) = residue(:) + rtmp*tmpmat(:)
                end if
             end do
          end do
          !$OMP END DO
          !$OMP CRITICAL
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL


          !: doubles alpha alpha
          !$OMP DO
          do i=1, noa
             do j=(i+1), noa
                do a=1, nva
                   do b=(a+1), nva
                      ia = state_index(-i,-j) - 1 + (a-1)*nva + a*(a-1)/2 + (b-a)
                      rtmp = b_i(ia,ibasis)
                      if ( rtmp.ne.0.d0 ) then
                         call get_cisd_iajb_AAAA(i,j,a,b,tmpmat)
                         residue(:) = residue(:) + rtmp*tmpmat(:)
                      end if
                   end do
                end do
             end do
          end do
          !$OMP END DO
          !$OMP CRITICAL
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL


          !: doubles beta beta
          !$OMP DO
          do i=1, nob
             do j=(i+1), nob
                do a=1, nvb
                   do b=(a+1), nvb
                      ia = state_index(i,j) - 1 + (a-1)*nvb + a*(a-1)/2 + (b-a)
                      rtmp = b_i(ia,ibasis)
                      if ( rtmp.ne.0.d0 ) then
                         call get_cisd_iajb_BBBB(i,j,a,b,tmpmat)
                         residue(:) = residue(:) + rtmp*tmpmat(:)
                      end if
                   end do
                end do
             end do
          end do
          !$OMP END DO
          !$OMP CRITICAL
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL


          !: doubles alpha beta
          !$OMP DO
          do i=1, noa
             do j=1, nob
                do a=1, noa
                   do b=1, nvb
                      ia = state_index(-i,j) - 1 + (a-1)*nvb + b
                      rtmp = b_i(ia,ibasis)
                      if ( rtmp.ne.0.d0 ) then
                         call get_cisd_iajb_ABAB(i,j,a,b,tmpmat)
                         residue(:) = residue(:) + rtmp*tmpmat(:)
                      end if
                   end do
                end do
             end do
          end do
          !$OMP END DO
          !$OMP CRITICAL
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL
          
          
          !$OMP END PARALLEL
          
       end do get_sigma


       
       bigG = 0.d0
       get_bigG : do ibasis=1, lbasis
          do jbasis=1, lbasis             
             bigG(jbasis,ibasis) = dot_product( b_i(:,jbasis), sigma(:,ibasis) )
          end do
       end do get_bigG
       

       
       !: write bigG
       if( iiter.eq.1 ) then
          do ibasis=1, lbasis
             write(iout,"(100(f10.5))") ( bigG(jbasis,ibasis), jbasis=1, lbasis ) 
          end do
          write(iout,'(A)') ' ----------------------------'
          flush(iout)
       end if
       

       !: diagonalize bigG
       lwork = 3*lbasis ; allocate( work(lwork) ) ; info = 10 ; eig=0.d0
       !: dysev(  'v'   ,  'u'  ,   N  ,  A , LDA  , W ,WORK,LWORK,INFO)
       call dsyev('vectors','upper',lbasis,bigG(1:lbasis,1:lbasis),lbasis,eig(1:lbasis),work,lwork,info)
       deallocate( work) 
       !: bigG contains eigenvectors (column vectors) ; eig contains the eigenvalues       
       

       !: write bigG
       if( iiter.eq.1 ) then
          do ibasis=1, lbasis
             write(iout,"(100(f10.5))") eig(ibasis), ( bigG(jbasis,ibasis), jbasis=1, lbasis )
          end do
       end if


       
       !: get residue vectors + correction vectors
       expand_space : do iroot=1, mroot
          
          rho_iroot = eig(iroot)
          residue   = 0.d0

          do ibasis=1, lbasis
             alpha_k    = bigG(ibasis,iroot)
             residue(:) = residue(:) + alpha_k * ( sigma(:,ibasis) - rho_iroot * b_i(:,ibasis) )
          end do
          
          !: standard Davidson-Liu Jacobi (preconditioner) algorithm
          delta = 0.d0
          do ia=1, nstates
             if ( diag(ia) - rho_iroot .eq. 0.d0 ) then
                delta(ia) = 0.d0
             else
                delta(ia) = - residue(ia) / ( diag(ia) - rho_iroot )
             end if
          end do
             
          !: normalize correction vectors
          norm2 = dot_product( delta, delta )
          delta = delta / sqrt(norm2)
          
          !: Gram-Schmidt orthogonalize
          tmpmat = delta
          do ibasis=1, lbasis
             projection = dot_product( delta(:), b_i(:,ibasis) )
             tmpmat(:) = tmpmat(:) - projection * b_i(:,ibasis)
          end do
          delta(:) = tmpmat(:)
          
             
          !: expand subspace
          norm2 = dot_product( delta, delta ) 
          if ( sqrt(norm2).ge.expand_cutoff ) then             
             lbasis = lbasis + 1 
             if ( lbasis.gt.max_basis ) then
                collapse = .true.
                go to 200                 
             end if
             b_i(:,lbasis) = delta(:) / sqrt(norm2)
          end if
          
          
          !: record residue and delta_energy for convergence
          res_diff(iroot)     = dot_product( residue, residue )
          energy_diff(iroot)  = abs( eig(iroot) - store_eig(iroot) )
          store_eig(iroot)    = eig(iroot)

       end do expand_space
       
       !-----------------------------!
       max_resdiff = maxval( res_diff )
       max_endiff  = maxval( energy_diff )
       !-----------------------------!              
       
       start_basis = lbasis0 + 1 

       
       !WRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITE
85     format(' root ',i4,' : ',f20.15,' eV || Delta_E ',f20.15,'   Residue',f20.15 )
       write(iout,"(' ITERATION ',i0,'    SUBSPACE ',i0)") iiter, lbasis0
       write(iout,85) (iroot, eig(iroot)*au2ev, energy_diff(iroot)*au2ev, res_diff(iroot), iroot=1, mroot )
       flush(iout)
       !WRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITEWRITE
       

200    continue
       
       if ( collapse ) then
          write(iout,"(' ITERATION ',i0,'    SUBSPACE ',i0)") iiter, lbasis0
          sigma = 0.d0 
          sigma = b_i 
          b_i   = 0.d0
          do i=1, mroot
             do j=1, lbasis0
                b_i(:,i) = b_i(:,i) + bigG(j,i)*sigma(:,j)
             end do
          end do
          
          start_basis = 1
          bigG   = 0.d0
          sigma  = 0.d0              
          lbasis = mroot
          write(iout,"(' COLLAPSING SUBSPACE to ',i0)") lbasis
          collapse = .false.
          flush(iout)            
       end if
       
       
    end do converge
    

    !: get transition dipoles. Use sigma as dummy array   
    call get_mo_elements
    call get_trans_elements(mroot, lbasis0, max_basis, sigma, b_i, bigG, tdipx, tdipy, tdipz)
    !call get_trans_elements2(mroot, lbasis0, max_basis, sigma, b_i, bigG, tdipx, tdipy, tdipz)
    

    write(iout,'(A)') divide
    write(iout,"(' CONVERGED STATES')")    
    
    do iroot=1, mroot
       
       write(iout,"(' State ',i4,' : ',f20.15,' eV ',3x,' trans_dip= ',f20.15 )") &
            iroot, eig(iroot)*au2ev, tdipx(iroot)**2 + tdipy(iroot)**2 + tdipz(iroot)**2
       
       write_states: do ia=1, nstates
          rtmp = sigma(ia,iroot)
          if ( abs(rtmp).gt.0.1d0 ) then
             x = hole_index(ia,1)
             i = hole_index(ia,2)
             a = part_index(ia,1)
             if ( i .eq. 0 ) write(iout,"('    x = ',i3,17x,f10.5)") x, rtmp
             if ( i .lt. 0 ) write(iout,"('    x = ',i3,' iA= ',i3,' aA = ',i3, f10.5)") x, -i, -a, rtmp
             if ( i .gt. 0 ) write(iout,"('    x = ',i3,' iB= ',i3,' aB = ',i3, f10.5)") x, i, a, rtmp
          end if
       end do write_states

    end do
    

    call cpu_time(lfinish)
    write(iout,'(" ")') 
    write(iout,"(' total diagonalization time:'f12.4,' seconds')") lfinish-lstart

    
    
  end subroutine davidson_cisd
  !--------------------------------------------------------------------------------------! 
  !--------------------------------------------------------------------------------------!
end module mod_davidson_cisd
