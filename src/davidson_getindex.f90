module davidson_getindex

  
  !: Davidson-Liu diagonalization for CISD-IP
  !: REF:  M.L. Leininger, C.D. Sherrill, W.D. Allen, H.F. Schaeffer III
  !:       J. Comp. Chem. 22, 1574-1589 (2001)
  !:       Systematic Study of Selected Diagonalization Methods for Configuration Interaction Matrices
  !: REF:  E.R. Davidson and W.J. Thompson, Computers in Physics, 7, 519-522 (1993)
  !:       Monster Matrices:  Their Eigenvalues and Eigenvectors
  
  use omp_lib
  use variables_global
  use util
  implicit none
  
  
contains
  !--------------------------------------------------------------------------------------! 
  !--------------------------------------------------------------------------------------!
  subroutine davidson_ip
    
    implicit none
    
    integer(8), parameter :: mroot = 5
    real(8),    parameter :: energy_tol    = 1.d-9
    real(8),    parameter :: residue_tol   = 1.d-7
    real(8),    parameter :: expand_cutoff = 0.01d0    
    
    integer(8) :: ibasis, jbasis, iroot, ia
    integer(8) :: lbasis0, lbasis, max_basis, start_basis
    integer(8) :: i, j, a, x

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
    real(8) :: start, finish
    
    
    
    call write_header( 'davidson_ucisdip', 'davidson', 'enter' )
    write(iout,"(' convergence on energy :  ',e10.4,' au ',e10.4,' eV ')") energy_tol, energy_tol*au2eV
    write(iout,"(' convergence on residue:  ',e10.4)") residue_tol
    write(iout,"(A)") '' 
    
    
    call cpu_time(start)           
    
    
    !: allocate guess space
    max_basis = 10*mroot
    allocate( b_i(nstates,max_basis) )    ; b_i   = 0.d0  
    allocate( sigma(nstates,max_basis) )  ; sigma = 0.d0 
    allocate( bigG(max_basis,max_basis) ) ; bigG  = 0.d0 
    allocate( eig(max_basis) )            ; eig   = 0.d0
    
    
    !: get guess vectors, all single excitation from occupied AND two lowest double excitations
    call ip_initial_guess( nstates, max_basis, lbasis, b_i, tmpmat )
    
    
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
    iiter  = 0   
    collapse = .false.
    max_resdiff = 10.d0 ; max_endiff = 10.d0 ; store_eig=0.d0 
    
    
    !: ITERATEITERATEITERATEITERATEITERATEITERATEITERATEITERATEITERATE
    converge : do while ( max_resdiff.ge.residue_tol .and. max_endiff.ge.energy_tol )        
       
       lbasis0 = lbasis 
       iiter   = iiter + 1        
       
       !: G_ji = < b_j | H | b_i > = < b_j | sigma_i >
       get_sigma : do ibasis=start_basis, lbasis
          
          !: use residue as temporary vector for OMP
          !$OMP PARALLEL DEFAULT(NONE), PRIVATE(ia, rtmp, tmpmat, residue), SHARED(b_i, ibasis, nstates, sigma)
          residue = 0.d0
          !$OMP DO
          do ia=1, nstates
             rtmp = b_i(ia,ibasis)
             if ( rtmp.ne.0.d0 ) then
                call get_h_column(ia,tmpmat)
                residue(:) = residue(:) + rtmp*tmpmat(:)
             end if
          end do
          !$OMP END DO
          !$OMP CRITICAL (form_sigma_lock)
          sigma(:,ibasis) = sigma(:,ibasis) + residue(:)
          !$OMP END CRITICAL (form_sigma_lock)
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
    

    call cpu_time(finish)
    write(iout,'(" ")') 
    write(iout,"(' total diagonalization time:'f12.4,' seconds')") finish-start

    
    
  end subroutine davidson_ip
  !--------------------------------------------------------------------------------------! 
  !--------------------------------------------------------------------------------------!
  subroutine ip_initial_guess(nstates, max_basis, lbasis, b_i, tmpmat)
    
    implicit none
    integer(8), intent(in) :: nstates, max_basis
    integer(8), intent(inout) :: lbasis
    real(8),    intent(inout) :: b_i( nstates,max_basis ), tmpmat(nstates)
    
    integer(8) :: x, i, a, ia
    integer(8) :: ibasis, xstart, mystate
    integer(8) :: store_ia(100)
    real(8)    :: curr_max

    !: initialize
    ibasis = 0
    xstart = nob-5
    if ( nob.le.10 ) xstart=2
    
    b_i = 0.d0
    store_ia = 0
    
    do x=xstart, nob

       !: single excitation
       ia = x
       ibasis = ibasis + 1
       store_ia(ibasis) = ia
       b_i( ia,ibasis ) = 1.d0
       write( iout,"(5x,' x= ',i0,' i= ',i0,' a= ',i0)") hole_index(ia,1), hole_index(ia,2), part_index(ia,1)

       !: get state with largest coupling to x
       call get_h_column( ia,tmpmat )
       do i=1, ibasis
          tmpmat( store_ia(i) ) = 0.d0
       end do
       do i=1, 2
          curr_max = 0.d0
          do ia=nob+1, nstates
             if( abs(tmpmat(ia)).gt.curr_max ) then
                curr_max = abs(tmpmat(ia))
                mystate = ia
             end if
          end do
          ia = mystate
          ibasis = ibasis + 1
          store_ia(ibasis) = ia
          b_i( ia,ibasis ) = 1.d0
          write( iout,"(5x,' x= ',i0,' i= ',i0,' a= ',i0)") hole_index(ia,1), hole_index(ia,2), part_index(ia,1)
          tmpmat(ia) = 0.d0
       end do

    end do

    tmpmat = 0.d0
    lbasis=ibasis
    

  end subroutine ip_initial_guess
  !:--------------------------!
  !:--------------------------!
  subroutine get_h_column(ia, hcol)

    use readintegrals
    implicit none

    integer(8), intent(in) :: ia
    real(8), intent(inout) :: hcol(nstates)
    
    integer(8) :: istate, xorb, iorb, aorb
    integer(8) :: xx, ii, aa, yy, jj, bb, x, y, i, a, j, b
    integer(8) :: jb, ia2, jb2, yx, xb, yj, ji, ix, jy, ab, ya, ib, ja, yb, xa, xi
    real(8) :: sign, storeme


    hcol = 0.d0

    xx = hole_index(ia,1) ; x = abs(xx)
    ii = hole_index(ia,2) ; i = abs(ii)
    aa = part_index(ia,1) ; a = abs(aa)

    !: orbital indices
    if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
    if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii
    if ( aa.lt.0 ) aorb = -aa + noa ; if ( aa.gt.0 ) aorb = nrorb + nob + aa

    do jb=1, nstates
       
       yy = hole_index(jb,1) ; y = abs(yy)
       jj = hole_index(jb,2) ; j = abs(jj)
       bb = part_index(jb,1) ; b = abs(bb)
       
       storeme = 0.d0
       
       SS : if( ii.eq.0 .and. jj.eq.0 ) then
          if( xx.eq.yy ) storeme = -Mol%orben(xorb)
          go to 78
       end if SS
       
       !: need to consider ii.ne.0 and jj.eq.0 if do loops go from jb=1
       SD : if( ii.eq.0 .and. jj.ne.0 ) then                 
          if ( jj.lt.0 ) storeme = - get_dijkaBA(y,j,x,b) !: <X|H|jb_Y> = -<Yj||Xb> 
          if ( jj.gt.0 ) storeme = - get_dijkaBB(y,j,x,b) !: <X|H|JB_Y> = -<YJ||XB> 
          go to 78
       else if ( jj.eq.0 .and. ii.ne.0 ) then          
          if ( ii.lt.0 ) storeme = - get_dijkaBA(x,i,y,a) !: <Y|H|ia_X> = -<Xi||Ya> 
          if ( ii.gt.0 ) storeme = - get_dijkaBB(x,i,y,a) !: <Y|H|IA_X> = -<XI||YA> 
          go to 78
       end if SD
       
       DD : if( ii.ne.0 .and. jj.ne.0 ) then
          
          kdelta_ab : if ( aa.eq.bb ) then
             if ( aa.lt.0 ) storeme = storeme + get_dijklAB(j,y,i,x) !: <jY||iX> 
             if ( aa.gt.0 ) storeme = storeme + get_dijklBB(j,y,i,x) !: <JY||IX> 
          end if kdelta_ab
          
          kdelta_ij : if ( ii.eq.jj ) then
             if ( aa.lt.0 ) storeme = storeme - get_diajbBA(y,a,x,b) !: -<Ya||Xb> 
             if ( aa.gt.0 ) storeme = storeme - get_diajbBB(y,a,x,b) !: -<YA||XB> 
          end if kdelta_ij
          
          kdelta_xy : if ( XX.eq.YY ) then
             if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme - get_diajbAA(j,a,i,b) !: -<ja||ib> 
             if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme - get_diajbBB(j,a,i,b) !: -<JA||IB> 
             if ( ii.lt.0 .and. JJ.gt.0 ) storeme = get_dijabAB(i,j,a,b) !: <iJ||aB> 
             if ( II.gt.0 .and. jj.lt.0 ) storeme = get_dijabAB(j,i,b,a) !: <jI||bA> 
          end if kdelta_xy
          
          kdelta_jx : if ( JJ.eq.XX ) then
             if ( ii.gt.0 ) storeme = storeme + get_diajbBB(y,a,i,b) !: <YA||IB> 
             if ( ii.lt.0 ) storeme = -get_dijabAB(i,y,a,b)          !: -<iY||aB> 
          end if kdelta_jx
          
          kdelta_iy : if ( II.eq.YY ) then
             if ( JJ.gt.0 ) storeme = storeme + get_diajbBB(j,a,x,b) !: <JA||XB> 
             if ( jj.lt.0 ) storeme = -get_dijabAB(j,x,b,a)          !: -<jX||bA>
          end if kdelta_iy
          
          diagonal : if ( jb.eq.ia ) then
             if ( aa.lt.0 ) storeme = storeme - Mol%orben(i) - Mol%orben(x) + Mol%orben(noa+a)
             if ( AA.gt.0 ) storeme = storeme - Mol%orben(nrorb+I) - Mol%orben(nrorb+X) + Mol%orben(nrorb+nob+A)
          end if diagonal
          
          special_off : if ( ii.eq.yy .and. xx.eq.jj .and. aa.eq.bb ) then
             storeme = storeme - Mol%orben(nrorb+I) - Mol%orben(nrorb+X) + Mol%orben(nrorb+nob+A)
          end if special_off
          
       end if DD
       
78     continue
       
       hcol(jb) = storeme
       
    end do !jb
    

  end subroutine get_h_column
  !:--------------------------!
  !:--------------------------!
  subroutine get_trans_elements(mroot, lbasis, max_basis, eig_states, b_i, bigG, dipx, dipy, dipz)
    
    implicit none
    
    integer(8), intent(in) :: mroot, max_basis, lbasis
    real(8), intent(in)    :: bigG(max_basis, max_basis), b_i(nstates,max_basis)
    real(8), intent(inout) :: eig_states(nstates,max_basis), dipx(mroot), dipy(mroot), dipz(mroot)
    
    integer(8) :: ia, iroot, ibasis
    
    real(8) :: coeff_ia
    real(8) :: densA(nrorb,nrorb), densB(nrorb,nrorb)
    real(8) :: dipx_ia(nstates), dipy_ia(nstates), dipz_ia(nstates)
    

    dipx = 0.d0
    dipy = 0.d0
    dipz = 0.d0

    
    !: sigma holds the eigenstates of the mroots
    eig_states = 0.d0
    do iroot=1, mroot
       do ibasis=1, lbasis
          eig_states(:,iroot) = eig_states(:,iroot) + bigG(ibasis,iroot)*b_i(:,ibasis)
       end do
    end do


    !$OMP PARALLEL DEFAULT(NONE), SHARED( eig_states, nstates, mroot, &
    !$OMP Mol, &
    !$OMP dipx_ia, dipy_ia, dipz_ia, dipx, dipy, dipz ),&
    !$OMP PRIVATE( ia, iroot, coeff_ia )
    !$OMP SECTIONS
    
    !$OMP SECTION
    do ia=1, nstates       
       coeff_ia = eig_states(ia,1)
       dipx_ia  = 0.d0
       call get_dip_col( ia, Mol%dipx00, Mol%dipxmoa, Mol%dipxmob, dipx_ia )       
       do iroot=1, mroot
          dipx(iroot) = dipx(iroot) + coeff_ia * dot_product( dipx_ia(:), eig_states(:,iroot) )
       end do
    end do
    
    !$OMP SECTION
    do ia=1, nstates       
       coeff_ia = eig_states(ia,1)  
       dipy_ia  = 0.d0
       call get_dip_col( ia, Mol%dipy00, Mol%dipymoa, Mol%dipymob, dipy_ia )       
       do iroot=1, mroot
          dipy(iroot) = dipy(iroot) + coeff_ia * dot_product( dipy_ia(:), eig_states(:,iroot) )
       end do
    end do
    
    !$OMP SECTION
    do ia=1, nstates       
       coeff_ia = eig_states(ia,1) 
       dipz_ia  = 0.d0
       call get_dip_col( ia, Mol%dipz00, Mol%dipzmoa, Mol%dipzmob, dipz_ia )       
       do iroot=1, mroot
          dipz(iroot) = dipz(iroot) + coeff_ia * dot_product( dipz_ia(:), eig_states(:,iroot) )
       end do
    end do

    !$OMP END SECTIONS
    !$OMP END PARALLEL
    
    
  end subroutine get_trans_elements
  !:--------------------------!
  !:--------------------------!
  subroutine get_trans_elements2(mroot, lbasis, max_basis, eig_states, b_i, bigG, dipx, dipy, dipz)
    
    implicit none
    
    integer(8), intent(in) :: mroot, max_basis, lbasis
    real(8), intent(in)    :: bigG(max_basis, max_basis), b_i(nstates,max_basis)
    real(8), intent(inout) :: eig_states(nstates,max_basis), dipx(mroot), dipy(mroot), dipz(mroot)
    
    integer(8) :: i, ia, iroot, ibasis
    
    real(8) :: rtmp
    real(8) :: densA(nrorb,nrorb), densB(nrorb,nrorb)
    real(8) :: diptmp(nrorb,nrorb)
    

    dipx = 0.d0
    dipy = 0.d0
    dipz = 0.d0


    !: sigma holds the eigenstates of the mroots
    eig_states = 0.d0
    do iroot=1, mroot
       do ibasis=1, lbasis
          eig_states(:,iroot) = eig_states(:,iroot) + bigG(ibasis,iroot)*b_i(:,ibasis)
       end do
    end do
    

    !$OMP PARALLEL DEFAULT(NONE), SHARED( eig_states, mroot, nrorb, &
    !$OMP Mol, dipx, dipy, dipz ),&
    !$OMP PRIVATE( iroot, i, rtmp, diptmp, densA, densB )
    !$OMP SECTIONS
    
    !$OMP SECTION
    do iroot=1, mroot
       call get_density( eig_states(:,1), eig_states(:,iroot), densA, densB )
       diptmp = matmul( densA, Mol%dipxmoa )
       rtmp = 0.d0
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       diptmp = matmul( densB, Mol%dipxmob )
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       dipx(iroot) = rtmp
    end do

    !$OMP SECTION
    do iroot=1, mroot
       call get_density( eig_states(:,1), eig_states(:,iroot), densA, densB )
       diptmp = matmul( densA, Mol%dipymoa )
       rtmp = 0.d0
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       diptmp = matmul( densB, Mol%dipymob )
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       dipy(iroot) = rtmp
    end do

    !$OMP SECTION
    do iroot=1, mroot
       call get_density( eig_states(:,1), eig_states(:,iroot), densA, densB )
       diptmp = matmul( densA, Mol%dipzmoa )
       rtmp = 0.d0
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       diptmp = matmul( densB, Mol%dipzmob )
       do i=1, nrorb
          rtmp = rtmp + diptmp(i,i)
       end do
       dipz(iroot) = rtmp
    end do
    
    !$OMP END SECTIONS
    !$OMP END PARALLEL
    
    
  end subroutine get_trans_elements2
  !:--------------------------!
  !:--------------------------!
  subroutine get_mo_elements

    implicit none


    allocate( Mol%dipxmoa(nrorb,nrorb), Mol%dipymoa(nrorb,nrorb), Mol%dipzmoa(nrorb,nrorb) )
    allocate( Mol%dipxmob(nrorb,nrorb), Mol%dipymob(nrorb,nrorb), Mol%dipzmob(nrorb,nrorb) )


    Mol%dipxmoa = 0.d0 ; Mol%dipymoa = 0.d0 ; Mol%dipzmoa = 0.d0
    Mol%dipxmob = 0.d0 ; Mol%dipymob = 0.d0 ; Mol%dipzmob = 0.d0

    !$OMP PARALLEL
    !$OMP SECTIONS

    !$OMP SECTION
    call ao2mo(nbasis,nrorb,Mol%dipxao,Mol%dipxmoa,Mol%cmo_a)
    call ao2mo(nbasis,nrorb,Mol%dipxao,Mol%dipxmob,Mol%cmo_b)
    
    !$OMP SECTION
    call ao2mo(nbasis,nrorb,Mol%dipyao,Mol%dipymoa,Mol%cmo_a)
    call ao2mo(nbasis,nrorb,Mol%dipyao,Mol%dipymob,Mol%cmo_b)

    !$OMP SECTION
    call ao2mo(nbasis,nrorb,Mol%dipzao,Mol%dipzmoa,Mol%cmo_a)
    call ao2mo(nbasis,nrorb,Mol%dipzao,Mol%dipzmob,Mol%cmo_b)
    
    !$OMP END SECTIONS
    !$OMP END PARALLEL


  end subroutine get_mo_elements
  !:--------------------------!
  !:--------------------------! 
  subroutine get_dip_col(ia, mo00, mo_mataa, mo_matbb, dip_col)

    implicit none

    integer(8), intent(in) :: ia
    real(8), intent(in)    :: mo00, mo_mataa(nrorb,nrorb), mo_matbb(nrorb,nrorb)
    real(8), intent(inout) :: dip_col(nstates)
    
    
    real(8) :: rdum
    integer(8) :: ii, jj, aa, bb, xx, yy, i, j, a, b, x, y
    integer(8) :: jb
    
    
    xx = hole_index(ia,1)
    ii = hole_index(ia,2)
    aa = part_index(ia,1)
    
    do jb=1, nstates
       
       yy = hole_index(jb,1)
       jj = hole_index(jb,2)
       bb = part_index(jb,1)
       
       rdum = 0.d0
       
       SS: if ( ii.eq.0 .and. jj.eq.0 ) then
          if ( xx*yy .gt. 0 ) then !: same sign
             if ( xx.lt.0 ) then
                x = -xx
                y = -yy
                rdum = - mo_matAA(y,x)
             else if ( xx.gt.0 ) then
                x = xx
                y = yy
                rdum = - mo_matBB(y,x)
             end if
          end if
          go to 78
       end if SS


       !: <y|mu|ia_y>
       SD1: if ( yy.eq.xx .and. jj.eq.0 ) then
          if ( ii.lt.0 ) then
             a = -aa + noa
             i = -ii
             rdum = mo_matAA(a,i)
          else if ( ii.gt.0 ) then
             a = aa + nob
             i = ii
             rdum = mo_matBB(a,i)
          end if
          go to 78
       end if SD1

       !: <jb_y|mu|x>
       SD2 : if ( yy.eq.xx .and. ii.eq.0 ) then
          if ( jj.lt.0 ) then
             b = -bb + noa
             j = -jj
             rdum = mo_matAA(b,j)
          else if ( jj.gt.0 ) then
             b = bb + nob
             j = jj
             rdum = mo_matBB(b,j)
          end if
          go to 78          
       end if SD2
       
       !: <y|mu|ya_x>
       SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
          same_spin0 : if ( xx*aa .gt. 0 ) then
             if ( xx.lt.0 ) then
                a = -aa + noa
                x = -xx
                rdum = -mo_matAA(a,x)
             else if ( xx.gt.0 ) then
                a = aa + nob
                x = xx
                rdum = -mo_matBB(a,x)
             end if
          end if same_spin0
          go to 78
       end if SD3

       !: <xb_y|mu|x>
       SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
          same_spin00 : if ( yy*bb .gt. 0 ) then
             if ( yy.lt.0 ) then
                b = -bb + noa
                y = -yy
                rdum = -mo_matAA(b,y)
             else if ( yy.gt.0 ) then
                b = bb + nob
                y = yy
                rdum = -mo_matBB(b,y)
             end if
          end if same_spin00
          go to 78
       end if SD4
       
       !: doubles
       
       kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
          if ( ii.lt.0 ) then
             i = -ii
             j = -jj
             rdum = -mo_matAA(i,j)
          else if ( ii.gt.0 ) then
             i = ii
             j = jj
             rdum = -mo_matBB(i,j)
          end if
       end if kdelta_xy_ab
       
       kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
          if ( aa.lt.0 ) then
             a = -aa + noa
             b = -bb + noa
             rdum = rdum + mo_matAA(a,b)
          else if ( aa.gt.0 ) then
             a = aa + nob
             b = bb + nob
             rdum = rdum + mo_matBB(a,b)
          end if
       end if kdelta_xy_ij
       
       kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
          same_spin1 : if ( xx*yy .gt. 0 ) then
             if ( xx.lt.0 ) then
                x = -xx
                y = -yy
                rdum = rdum - mo_matAA(x,y)
             else
                x = xx
                y = yy
                rdum = rdum - mo_matBB(x,y)
             end if
          end if same_spin1
       end if kdelta_ij_ab
       
       kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
          same_spin2 : if ( ii*yy .gt. 0 ) then
             if ( ii.lt.0 ) then
                i = - ii
                y = - yy
                rdum = rdum + mo_matAA(i,y)
             else
                i = ii
                y = yy
                rdum = rdum + mo_matBB(i,y)
             end if
          end if same_spin2
       end if kdelta_jx_ab
       
       kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
          same_spin3 : if ( xx*jj .gt. 0 ) then
             if ( jj.lt.0 ) then
                j = -jj
                x = -xx
                rdum = rdum + mo_matAA(j,x)
             else
                j = jj
                x = xx
                rdum = rdum + mo_matBB(j,x)
             end if
          end if same_spin3
       end if kdelta_yi_ab
       
       kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
          same_spin4 : if ( aa*bb .gt. 0 ) then
             if ( aa.lt.0 ) then
                a = -aa + noa
                b = -bb + noa
                rdum = rdum - mo_matAA(a,b)
             else if ( aa.gt.0 ) then
                a = aa + nob
                b = bb + nob
                rdum = rdum - mo_matBB(a,b)
             end if
          end if same_spin4
       end if kdelta_iy_xj
       
78     continue
       if( ia.eq.jb ) rdum = rdum + mo00
       
       dip_col(jb) = rdum
       
    end do
    
    
  end subroutine get_dip_col
  !:--------------------------!
  !:--------------------------! 
  subroutine get_density( psi1, psi2, rhoA, rhoB )

    implicit none

    real(8), intent(in) :: psi1(nstates), psi2(nstates)
    real(8), intent(inout) :: rhoA(nrorb,nrorb), rhoB(nrorb,nrorb)

    real(8) :: coeff_ia, coeff_jb
    integer(8) :: ia, jb, xx, ii, aa, yy, jj, bb
    integer(8) :: i, j, a, b, x, y

    
    rhoA = 0.d0
    rhoB = 0.d0

    do ia=1, nstates
       
       xx = hole_index(ia,1)
       ii = hole_index(ia,2)
       aa = part_index(ia,1)

       coeff_ia = psi1(ia)
       
       do jb=1, nstates
          
          yy = hole_index(jb,1)
          jj = hole_index(jb,2)
          bb = part_index(jb,1)

          coeff_jb = psi2(jb)
          
          SS : if ( jj.eq.0 .and. ii.eq.0 ) then
             if ( xx.gt.0 ) then
                x = xx
                y = yy
                rhoB(y,x) = rhoB(y,x) - coeff_ia * coeff_jb
             end if
             go to 78
          end if SS
          
          SD1 : if ( xx.eq.yy .and. ii.eq.0 ) then
             if ( jj.lt.0 ) then
                j = -jj
                b = -bb + noa
                rhoA(b,j) = rhoA(b,j) + coeff_ia * coeff_jb
             else if ( jj.gt.0 ) then
                j = jj
                b = bb + nob
                rhoB(b,j) = rhoB(b,j) + coeff_ia * coeff_jb
             end if
             go to 78
          end if SD1
          
          SD2 : if ( xx.eq.yy .and. jj.eq.0 ) then
             if ( ii.lt.0 ) then
                i = -ii
                a = -aa + noa
                rhoA(a,i) = rhoA(a,i) + coeff_ia * coeff_jb
             else if ( ii.gt.0 ) then
                i = ii
                a = aa + nob
                rhoB(a,i) = rhoB(a,i) + coeff_ia * coeff_jb
             end if
             go to 78
          end if SD2
          
          SD3 : if ( xx.eq.jj .and. ii.eq.0 ) then
             if ( yy*bb .gt. 0 ) then
                if ( yy.lt.0 ) then
                   y = -yy
                   b = -bb + noa
                   rhoA(b,y) = rhoA(b,y) - coeff_ia * coeff_jb
                else if ( yy.gt.0 ) then
                   y = yy
                   b = bb + nob
                   rhoB(b,y) = rhoB(b,y) - coeff_ia * coeff_jb
                end if
             end if
          end if SD3
          
          SD4:  if ( yy.eq.ii .and. jj.eq.0 ) then
             if ( xx*aa .gt. 0 ) then
                if ( xx.lt.0 ) then
                   x = -xx
                   a = -aa + noa
                   rhoA(a,x) = rhoA(a,x) - coeff_ia * coeff_jb
                else if ( xx.gt.0 ) then
                   x = xx
                   a = aa + nob
                   rhoB(a,x) = rhoB(a,x) - coeff_ia * coeff_jb
                end if
             end if
          end if SD4

          !: doubles
          xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             if ( ii.lt.0 ) then
                i = -ii
                j = -jj
                rhoA(j,i) = rhoA(j,i) - coeff_ia * coeff_jb
             else if  ( ii.gt.0 ) then
                i = ii
                j = jj
                rhoB(j,i) = rhoB(j,i) - coeff_ia* coeff_jb
             end if
          end if xy_ab

          xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             if ( aa.lt.0 ) then
                a = -aa + noa
                b = -bb + noa
                rhoA(b,a) = rhoA(b,a) + coeff_ia * coeff_jb
             else if ( aa.gt.0 ) then
                a = aa + nob
                b = bb + nob
                rhoB(b,a) = rhoB(b,a) + coeff_ia * coeff_jb
             end if
          end if xy_ij

          ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             if ( xx.lt.0 ) then
                x = -xx 
                y = -yy
                rhoA(y,x) = rhoA(y,x) - coeff_ia * coeff_jb
             else if ( xx.gt.0 ) then
                x = xx
                y = yy
                rhoB(y,x) = rhoB(y,x) - coeff_ia * coeff_jb
             end if
          end if ij_ab

          jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             if ( ii*yy.gt.0 ) then
                if ( ii.lt.0 ) then
                   i = -ii
                   y = -yy
                   rhoA(y,i) = rhoA(y,i) + coeff_ia * coeff_jb
                else if ( ii.gt.0 ) then
                   i = ii
                   y = yy
                   rhoB(y,i) = rhoB(y,i) + coeff_ia * coeff_jb
                end if
             end if
          end if jx_ab

          yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             if ( jj*xx.gt.0 ) then
                if ( jj.lt.0 ) then
                   j = -jj
                   x = -xx
                   rhoA(x,j) = rhoA(x,j) + coeff_ia * coeff_jb
                else if ( jj.gt.0 ) then
                   j = jj
                   x = xx
                   rhoB(x,j) = rhoB(x,j) + coeff_ia * coeff_jb
                end if
             end if
          end if yi_ab

          iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             if ( aa.lt.0 ) then
                a = -aa + noa
                b = -bb + noa
                rhoA(b,a) = rhoA(b,a) - coeff_ia * coeff_jb
             else if ( aa.gt.0 ) then
                a = aa + nob
                b = bb + nob
                rhoB(b,a) = rhoB(b,a) - coeff_ia * coeff_jb
             end if
          end if iy_xj
          
78        continue

          !: diagonal
          if ( ia.eq.jb ) then
             do i=1, noa
                rhoA(i,i) = rhoA(i,i) + coeff_ia * coeff_jb
             end do
             do i=1, nob
                rhoB(i,i) = rhoB(i,i) + coeff_ia * coeff_jb
             end do
          end if
          
       end do
    end do
  !:--------------------------!
  !:--------------------------! 
  end subroutine get_density
end module davidson_getindex
