module analysis
  
  !use variables_global ! <C> use at your risk, beware of private & shared OMP variables
  use variables_setup ! MolInfo class definition
  use util
  implicit none

  
contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_NORM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_norm( norm, nstuse, psi )
    
    ! <C> Get norm of a complex vector
    
    real(8), intent(inout) :: norm    
    integer(8), intent(in) :: nstuse
    complex(8), intent(in) :: psi(nstuse)
    
    integer(8) :: i
    real(8)    :: rdum, cdum
    
     
    norm = dot_product( psi, psi )
    norm = dsqrt(norm)

    
  end subroutine get_norm
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_EXPECTATION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_expectation(nstuse,norm,psi,obsv,expect_value)
    
    ! <C> get expectation value of obsv

    implicit none

    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: norm
    real(8)   , intent(in) :: obsv(nstuse*nstuse)
    complex(8), intent(in) :: psi(nstuse)
    real(8), intent(inout) :: expect_value
    
    integer(8) :: i, ii, j, ij
    complex(8) :: psi_i

    
    expect_value = 0.d0
    
    
    ! <C> off-diagonal
    do i = 2, nstuse
       ii = (i-1)*nstuse
       psi_i = dconjg( psi(i) )
       do j = 1, (i-1)
          ij = ii + j
          expect_value = expect_value + obsv(ij) * real( psi_i*psi(j) )
       end do
    end do

    expect_value = 2.d0 * expect_value
    
    ! <C> diagonal
    do i = 1, nstuse
       ii = (i-1)*nstuse + i
       expect_value = expect_value + obsv(ii) * dconjg(psi(i))*psi(i) 
    end do
    
    
    if(norm.ne.0.d0) expect_value = expect_value / norm**2
    
    
  end subroutine get_expectation
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_ZEXPECTATION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_Zexpectation(nstuse,norm,psi,Zobsv,psi1,expect_value)
    
    ! <C> get expectation value of obsv

    implicit none

    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: Zobsv(nstuse*nstuse)
    complex(8), intent(in) :: psi(nstuse), psi1(nstuse)
    real(8), intent(inout) :: expect_value
    
    integer(8) :: i, ii, j, ij
    complex(8) :: psi_i, cdum, zero, zdotc
    real(8)    :: expect_value1

    zero = dcmplx(0.d0,0.d0)
    cdum = dcmplx(1.d0,0.d0)
    if(norm.ne.0) cdum = dcmplx(1.d0/norm**2,0.d0)
    call zgemv('n',nstuse,nstuse,cdum,Zobsv,nstuse,psi,1,zero,psi1,1)
    expect_value = dble(zdotc(nstuse,psi,1,psi1,1))
 
!:    expect_value1 = 0.d0
!:    cdum = zero
!:    
!:    do i=1, nstuse
!:       ii = (i-1)*nstuse
!:       psi_i = dconjg( psi(i) )
!:       do j=1, nstuse
!:          ij = i + (j-1)*nstuse
!:          cdum = cdum + Zobsv(ij) * psi(j) * psi_i
!:       end do
!:    end do
!:
!:    !if ( aimag(cdum).ne.0.d0 ) then
!:    !   write(*,*) cdum
!:    !end if
!:    
!:    expect_value1 = dble(cdum) / norm**2
!:    write(42,*) expect_value1,expect_value
    
  end subroutine get_Zexpectation
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET PSI_DET
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_psid(nstuse,nstates,cis_vec,norm,psi,psi_det)

    implicit none

    integer(8), intent(in) :: nstuse, nstates
    real(8),    intent(in) :: cis_vec(nstates*nstates), norm
    complex(8), intent(in) :: psi(nstuse)
    complex(8), intent(inout) :: psi_det(nstates)

    integer(8) :: k, kk
    complex(8) :: psi_k

    psi_det = dcmplx(0.d0,0.d0)

    do k = 1, nstuse
       kk = nstates * (k-1)
       psi_k = psi(k) 
       psi_det(:) = psi_det(:) + cis_vec(kk+1:kk+nstates)*psi_k
    end do
    
    psi_det = psi_det / norm
    
  end subroutine get_psid
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET ZPSID
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_Zpsid(nstuse,nstates,Zcis_vec,norm,psi,psi_det)

    implicit none

    integer(8), intent(in) :: nstuse, nstates
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: Zcis_vec(nstates*nstates)
    complex(8), intent(in) :: psi(nstuse)
    complex(8), intent(inout) :: psi_det(nstates)

    integer(8) :: k, kk
    complex(8) :: psi_k

!:    write(42,*) "Zpsid",nstuse,norm

    psi_det = dcmplx(0.d0,0.d0)

    do k = 1, nstuse
       kk = nstates * (k-1)
       psi_k = psi(k) 
!:       psi_det(:) = psi_det(:) + dconjg(Zcis_vec(kk+1:kk+nstates))*psi_k
       psi_det(:) = psi_det(:) + Zcis_vec(kk+1:kk+nstates)*psi_k
    end do
    
    psi_det = psi_det / norm
    
  end subroutine get_Zpsid
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_CIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_cis( hole,noa,norb,nstates,nva,part,pop,psi_det,nob,nvb)
    
    ! <C> Calculate MO occupancies for a CIS vector

    implicit none
    
    integer(8), intent(in) :: noa, norb, nstates, nva
    integer(8), intent(in) :: hole(nstates), part(nstates)
    complex(8), intent(in) :: psi_det(nstates)
    real(8), intent(inout) :: pop(norb)
    integer(8), optional, intent(in) :: nob, nvb

    integer(8) :: i, a, ia, ii, aa
    real(8)    :: psi2
    
    
    pop = 0.d0
    
    if ( present(nob) ) then
       pop(1:noa) = 1.d0
       pop(noa+nva+1:noa+nva+nob) = 1.d0
    else
       pop(1:noa) = 2.d0
    end if
    
    do ia = 2, nstates
       ii   = hole(ia) 
       aa   = part(ia)
       psi2 = dble(dconjg(psi_det(ia))*psi_det(ia))
       if ( ii.lt.0 ) then
          i = -ii
          pop(i) = pop(i) - psi2
       else
          i = ii + noa + nva
          pop(i) = pop(i) - psi2
       end if
       if ( aa.lt.0 ) then
          a = -aa + noa          
          pop(a) = pop(a) + psi2
       else
          a = aa + noa + nva + nob
          pop(a) = pop(a) + psi2
       end if
    end do

    
  end subroutine pop_cis
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_IP
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ip( hole,noa,nob,norb,nstates,nva,nvb,part,pop,psi_det )
    
    implicit none
    
    integer(8), intent(in) :: noa, nob, norb, nstates, nva, nvb
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    real(8), intent(inout) :: pop(norb)
    complex(8), intent(in) :: psi_det(nstates)

    integer(8) :: ia, ii, aa, xx, i, a, x, j
    integer(8) :: nrorb
    real(8)    :: psi2
    

    pop   = 0.d0
    nrorb = noa + nva
    

    pop(1:noa) = 1.d0
    pop( (nrorb+1):(nrorb+nob) ) = 1.d0

    do ia=1, nstates
       
       xx   = hole(ia,1) 
       ii   = hole(ia,2)
       aa   = part(ia)

       x = xx + nrorb       
       i = ii + nrorb
       a = aa + nrorb + nob 
       if ( xx.lt.0 ) x = -xx
       if ( ii.lt.0 ) i = -ii
       if ( aa.lt.0 ) a = -aa + noa
       
       psi2 = dble(dconjg(psi_det(ia)) * psi_det(ia))
       
       S_or_D : if ( ii.eq.0 ) then
          
          pop(x) = pop(x) - psi2
          
       else

          pop(x) = pop(x) - psi2
          pop(i) = pop(i) - psi2
          pop(a) = pop(a) + psi2
          
       end if S_OR_D
    
    end do


  end subroutine pop_ip
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_UION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ion( noa,nva,norb,nstates,psi_det,vabsmoa,pop,nob,nvb,vabsmob)
    
    implicit none
    
    integer(8), intent(in) :: noa, nva, norb, nstates
    real(8),    intent(in) :: vabsmoa(noa+nva,noa+nva)
    complex(8), intent(in) :: psi_det(nstates)
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: vabsmob(nob+nvb,nob+nvb)
    real(8),    intent(inout) :: pop(norb)
    
    integer(8) :: i, a, b, ia, ib
    real(8)    :: psi2, rdum
    complex(8) :: psi0
    
    
    pop = 0.d0
    
    ! <C> diagonal <ia|V|ia> alpha
    do i=1, noa
       rdum = 0.d0
       do a=1, nva
          ia = (i-1)*nva + a + 1
          psi0 = dconjg(psi_det(ia))
          rdum = rdum + dconjg(psi0)*psi0 * vabsmoa(noa+a,noa+a)
          do b=(a+1), nva
             ib = (i-1)*nva + b + 1 
             rdum = rdum + 2.d0*real(psi0*psi_det(ib)) * vabsmoa(noa+a,noa+b) 
          end do
       end do
       pop(i) = rdum
    end do
    

    ! <C> diagonal <ia|V|ia> beta
    if ( present(nob) ) then
       do i=1, nob
          rdum = 0.d0
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             psi0 = dconjg(psi_det(ia))
             rdum = rdum + dconjg(psi0)*psi0 * vabsmob(nob+a,nob+a)
             do b=(a+1), nvb
                ib = (i-1)*nvb + b + noa*nva + 1
                rdum = rdum + 2.d0*real(psi0*psi_det(ib)) * vabsmob(nob+a,nob+b)
             end do
          end do
          pop(noa+nva+i) = rdum
       end do
    end if
    
    
  end subroutine pop_ion
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_ION2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ion2( noa,nva,norb,nstates,psi_det,vabsmoa,pop,nob,nvb,vabsmob )
    
    implicit none

    integer(8), intent(in) :: noa,nva,norb,nstates
    real(8)   , intent(in) :: vabsmoa(noa+nva,noa+nva)
    complex(8), intent(in) :: psi_det(nstates)
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: vabsmob(nob+nvb,nob+nvb)
    real(8), intent(inout) :: pop(norb)
    
    integer(8) :: i,j,k,a,ia,b,jb
    real(8) :: tmp(norb,norb), dens(norb,norb), const
    complex(8) :: psi2, psi0
    
    
    dens = 0.d0 ; pop = 0.d0
    tmp(1:noa+nva,1:noa+nva) = vabsmoa(:,:)

    
    if ( present(nob) ) then
       tmp( noa+nva+1:norb, noa+nva+1:norb ) = vabsmob(:,:)
       const = 1.d0 
       do i=1, noa
          dens(i,i) = 1.d0
       end do
       do i=1, nob
          dens(noa+nva+i,noa+nva+i) = 1.d0
       end do
    else
       const = dsqrt(2.d0)
       do i=1, noa
          dens(i,i) = 2.d0
       end do
    end if
    

    ! <C> get one-electron density matrix < 0 | ia >_(N-1) alpha
    psi0 = dconjg(psi_det(1))
    do i=1, noa
       do a=1, nva
          ia = (i-1)*nva + a + 1 
          dens(i,noa+a) = const * real( psi0 * psi_det(ia) )
          dens(noa+a,i) = const * real( psi0 * psi_det(ia) )
       end do
    end do
    
    ! <C> < ia | jb > alpha
    do ia=1, noa*nva
       i = (ia-1)/nva + 1 
       a = ia - (i-1)*nva + noa
       psi2 = dconjg(psi_det(ia+1))
       do jb=1, noa*nva
          j = (jb-1)/nva + 1 
          b = jb - (j-1)*nva + noa
          if( i.eq.j ) dens(b,a) = dens(b,a) + real(psi2*psi_det(jb+1))
          if( a.eq.b ) dens(j,i) = dens(j,i) - real(psi2*psi_det(jb+1))
       end do
    end do

    
    ! <C> < ia | jb > beta
    if ( present(nob) ) then
       ! <C> get one-electron density matrix < 0 | ia >_(N-1) beta
       psi0 = dconjg(psi_det(1))
       do i=1, nob
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             dens(noa+nva+i,noa+nva+nob+a) = real( psi0 * psi_det(ia) )
             dens(noa+nva+nob+a,noa+nva+i) = real( psi0 * psi_det(ia) )
          end do
       end do
       ! <C> < ia | jb > beta
       do ia=1, nob*nvb
          i = (ia-1)/nvb + 1 
          a = ia - (i-1)*nvb + noa + nva + nob
          i = i + noa + nva
          psi2 = dconjg(psi_det(ia+1+noa*nva))
          do jb=1, nob*nvb
             j = (jb-1)/nvb + 1 
             b = jb - (j-1)*nvb + noa + nva + nob
             j = j + noa + nva
             if( i.eq.j ) dens(b,a) = dens(b,a) + real(psi2*psi_det(jb+1+noa*nva))
             if( a.eq.b ) dens(j,i) = dens(j,i) - real(psi2*psi_det(jb+1+noa*nva))
          end do
       end do
    end if
    

    tmp = matmul( dens, tmp ) 
    
    do i=1, norb
       pop(i) = tmp(i,i) 
    end do
    

  end subroutine pop_ion2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_DYSON
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_dyson( noa,nva,norb,torbitals,nstuse,nstates,nbasis,cis_vec, &
                        cmo_a,psi,psi0,dirstr,emaxstr,nob,nvb,cmo_b )
    
    ! <C> use at your risk.  need to be double-checked

    implicit none
    
    integer(8),  intent(in) :: noa, nva, norb, torbitals, nstates, nstuse, nbasis
    real(8),     intent(in) :: cis_vec(nstates*nstates), cmo_a(nbasis,torbitals)
    complex(8),  intent(in) :: psi(nstuse), psi0(nstuse)
    character(4),intent(in) :: dirstr, emaxstr
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: cmo_b(nbasis,torbitals)
    
    integer(8) :: i,j,k,a,ia,b,jb, ibasis, jbasis, iorb, jorb, i0, den0
    real(8)    :: const, c, norm, norm0, dens_mo(torbitals,torbitals), dens_ao(nbasis*(nbasis+1)/2)
    complex(8) :: psi_det(nstates), psi_det0(nstates), psi2
    character(100) :: outfile
    
    
    norm0 = 1.d0
    norm = dot_product( psi,psi ) ; norm = dsqrt(norm)
    call get_psid(nstuse,nstates,cis_vec,norm0,psi0,psi_det0)
    call get_psid(nstuse,nstates,cis_vec,norm,psi,psi_det)
    
   
    dens_mo = 0.d0 

    den0 = 2.d0 ; const = dsqrt(2.d0)
    if ( present(nob) ) then
       den0 = 1.d0
       const = 1.d0
    end if

    i0 = torbitals - norb/2

    ! <C> ALPHA DENSITY 
    do i=1, i0
       dens_mo(i,i) = den0
    end do
    do i=1, noa
       dens_mo(i+i0,i+i0) = den0
    end do
    
    ! <C>  < 0 | ia >_(N-1) alpha
    psi2 = dconjg(psi_det0(1))
    do i=1, noa
       do a=1, nva
          ia = (i-1)*nva + a + 1 
          dens_mo(i+i0,noa+a+i0) = const * real( psi2 * psi_det(ia) )
          dens_mo(noa+a+i0,i+i0) = const * real( psi2 * psi_det(ia) )
       end do
    end do
    
    ! <C> < ia | jb > alpha
    do ia=1, noa*nva
       i = (ia-1)/nva + 1 
       a = ia - (i-1)*nva + noa + i0
       i = i + i0
       psi2 = dconjg(psi_det0(ia+1))
       do jb=1, noa*nva
          j = (jb-1)/nva + 1 
          b = jb - (j-1)*nva + noa + i0
          j = j + i0
          if( i.eq.j ) dens_mo(b,a) = dens_mo(b,a) + real(psi2*psi_det(jb+1))
          if( a.eq.b ) dens_mo(j,i) = dens_mo(j,i) - real(psi2*psi_det(jb+1))
       end do
    end do
    
    dens_ao = 0.d0 ; k = 0
    ! <C> get rho_ij to rho_uv
    do ibasis=1, nbasis
       do jbasis=1, ibasis
          k = k + 1 ; c = 0.d0
          do iorb=1, torbitals
             do jorb=1, torbitals
                c = c + cmo_a(ibasis,iorb) * cmo_a(jbasis,jorb) * dens_mo(jorb,iorb)
             end do
          end do
          if ( abs(c).lt.1.d-10 ) c = 0.d0
          dens_ao(k) = c
       end do
    end do


    
    ! <C> BETA DENSITY
    beta_den : if( present(nob) ) then

       dens_mo = 0.d0
       do i=1, i0
          dens_mo(i,i) = den0
       end do
       do i=1, nob
          dens_mo(i+i0,i+i0) = den0
       end do
       
       ! <C> get one-electron dens_moity matrix < 0 | ia >_(N-1) beta
       psi2 = dconjg(psi_det0(1))
       do i=1, nob
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             dens_mo(i+i0,nob+a+i0) = real( psi2 * psi_det(ia) )
             dens_mo(nob+a+i0,i+i0) = real( psi2 * psi_det(ia) )
          end do
       end do
       
       ! <C> < ia | jb > beta
       do ia=1, nob*nvb
          i = (ia-1)/nvb + 1 
          a = ia - (i-1)*nva + nob + i0
          i = i + i0
          psi2 = dconjg(psi_det(ia+1+noa*nva))
          do jb=1, nob*nvb
             j = (jb-1)/nvb + 1 
             b = jb - (j-1)*nvb + nob + i0
             j = j + i0
             if( i.eq.j ) dens_mo(b,a) = dens_mo(b,a) + real(psi2*psi_det(jb+1+noa*nva))
             if( a.eq.b ) dens_mo(j,i) = dens_mo(j,i) - real(psi2*psi_det(jb+1+noa*nva))
          end do
       end do
       
       k = 0
       ! <C> get rho_ij to rho_uv
       do ibasis=1, nbasis
          do jbasis=1, ibasis
             k = k + 1 ; c = 0.d0
             do iorb=1, torbitals
                do jorb=1, torbitals
                   c = c + cmo_b(ibasis,iorb) * cmo_b(jbasis,jorb) * dens_mo(jorb,iorb)
                end do
             end do
             if ( abs(c).lt.1.d-10 ) c = 0.d0
             dens_ao(k) = dens_ao(k) + c
          end do
       end do
       
    end if beta_den

    outfile='DYSON'//'-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.out'
    open( unit=434,file=trim(outfile) )
    write(434,"(5es16.8)") ( dens_ao(k), k=1, nbasis*(nbasis+1)/2 )
    close(434)


    
  end subroutine get_dyson
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_MODENSITY
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_modensity( noa,nva,norb,nstates,hole,part,psi_det,norm,dens,nob,nvb )

    implicit none 
    integer(8), intent(in) :: noa, nva, norb, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: psi_det(nstates)

    real(8), intent(inout) :: dens(norb,norb)
    integer(8), optional, intent(in) :: nob, nvb
    
    integer(8) :: ia, jb, i, a, j, b, ii, aa, jj, bb
    real(8)    :: const
    complex(8) :: c0, cdum

    
    const = dsqrt(2.d0)
    if( present(nob) ) const = 1.d0

    dens = 0.d0
    do i=1, noa
       dens(i,i) = norm * norm
    end do

    if ( present(nob) ) then
       do i = (noa+nva+1) , (noa+nva+nob )
          dens(i,i) =  norm * norm
       end do
    end if

    
    ! <C> 0,ia
    c0 = dconjg( psi_det(1) ) 
    do ia = 2, nstates
       
       ii = hole(ia)
       aa = part(ia)
       
       i = -ii
       a = -aa + noa
       if ( ii.gt.0 ) i = ii + noa + nva
       if ( aa.gt.0 ) a = aa + noa + nva + nob
       
       cdum = c0 * psi_det(ia)
       dens(a,i) = abs( cdum ) **2
       dens(i,a) = abs( cdum ) **2
       
    end do
    
    
    do ia = 2, nstates

       ii = hole(ia)
       aa = part(ia)
       i = -ii
       a = -aa + noa
       if ( ii.gt.0 ) i = ii + noa + nva
       if ( aa.gt.0 ) a = aa + noa + nva + nob

       c0 = dconjg(psi_det(ia))
       do jb = 2, nstates
          
          jj = hole(jb)
          bb = part(jb)
          j = -jj
          b = -bb + noa 
          if ( jj.gt.0 ) j = jj + noa + nva
          if ( bb.gt.0 ) b = bb + noa + nva  + nob
          
          if ( i.eq.j ) dens(a,b) = dens(a,b) + abs(c0 * psi_det(jb))**2
          if ( a.eq.b ) dens(j,i) = dens(j,i) - abs(c0 * psi_det(jb))**2
          
       end do
    end do


  end subroutine get_modensity
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FOURIER TRANSFORM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine fourier_transform( ntime, dt, corrt, outfile )
  
    implicit none

    real(8), parameter :: pi = 2.d0*acos(0.d0)      ! pi
    real(8), parameter :: au2eV = 27.21138602d0     ! hartrees to eV
    real(8), parameter :: au2fs = 2.418884326509d-2 ! au time to fs 
    
    integer(8), intent(in) :: ntime
    real(8),    intent(in) :: dt
    complex(8), intent(in) :: corrt(0:ntime)
    character(100) , intent(in) :: outfile
    

    ! <C> blackman window
    real(8), parameter :: a0 = 7938.d0/18608.d0
    real(8), parameter :: a1 = 9240.d0/18608.d0
    real(8), parameter :: a2 = 1430.d0/18608.d0

    integer(8) :: it, iw, i
    real(8)    :: pic, window(0:ntime), dw
    complex(8) :: corrw(0:ntime), c, expk0, expk
    

    corrw = dcmplx(0.d0,0.d0)
    dw = 2.d0 * pi / dt / dble(ntime)
    pic   = 2.d0 * pi / dble(ntime)

    ! <C> assign window function
    !do it=0, ntime
    !   window(it) = a0 - a1*cos(pic*dble(it)) + a2*cos(2.d0*pic*dble(it))
    !end do

    
    do iw=0, int(ntime/2)
       c = corrt(0)
       !expk0 = exp( eye * dble(iw) * pic )
       expk0 = dcmplx( cos(-dble(iw)*pic), sin(-dble(iw)*pic) )
       expk  = dcmplx( 1.d0,0.d0 )
       do it=1, ntime
          expk = expk * expk0
          c = c + corrt(it) * expk
       end do
       corrw(iw) = c
    end do

    open( unit=100,file=trim(outfile) )
    write(100,"( 2a20,'|',2a20 )") 'omega(eV)', 'corr(w)', 'time(fs)' ,'corr(t)'
    do i=0, ntime
       write(100,"(2f20.10,'|',2f20.10)")  dw*dble(i)*au2eV, abs( corrw(i) ), dt*dble(i)*au2fs, abs( corrt(i) )
    end do
    close(100)
    
    
  end subroutine fourier_transform
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITE_DATA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine write_data(nstep,outstep,norb,ndata,dt,noa,nva,nob,nvb,normdata,ratedata,&
       popdata,iondata,mudata,normfile,popfile,ionfile,dipolefile,tdx,tdy,tdz,        &
       iout,writenorm,writepop,writeion,writemu,unrestricted)
    

    implicit none

    real(8), parameter :: au2fs = 2.418884326509d-2

    integer(8), intent(in) :: nstep, outstep, norb, ndata, noa, nva, nob, nvb
    integer(8), intent(in) :: iout
    real(8),    intent(in) :: dt, tdx, tdy, tdz
    real(8),    intent(in) :: normdata(ndata), ratedata(ndata)
    real(8),    intent(in) :: popdata(norb,ndata), iondata(norb,ndata), mudata(3,ndata)
    character(60), intent(in) :: normfile, popfile, ionfile, dipolefile
    logical, intent(in) :: writenorm, writepop, writeion, writemu, unrestricted

    integer(8) :: idata, j
    
    
    if( writenorm ) then
       open( 18, file=trim(normfile) )
       write(18,"( 2(1x,a9),2(1x,a12) )") 'time(au)','time(fs)','norm', 'rate(fs)'
       write(18,"( 2(1x,f9.3),2(1x,f13.10) )") 0.d0, 0.d0, 1.d0, 0.d0
       do idata = 1, ndata
          write(18,"( 2(1x,f9.3),2(1x,f13.10) )")   &
               dble(outstep)*dble(idata)*dt,       &
               dble(outstep)*dble(idata)*dt*au2fs, &
               normdata(idata), ratedata(idata)/au2fs
       end do
       close(18)
       write(iout,100) 'norm written to                         ',trim(normfile)
    end if
    

    if ( writepop ) then
       open( 19, file=trim(popfile) )
       do idata = 1, ndata
          write(19,"(' TIME(fs)= ',f9.2)") (outstep)*dble(idata)*dt*au2fs
          if ( unrestricted ) then
             write(19,101) ( popdata(j,idata), j=1,noa )
             write(19,300) ( popdata(j,idata), j=noa+1,noa+nva )
             write(19,200) ( popdata(j,idata), j=noa+nva+1,noa+nva+nob )
             write(19,400) ( popdata(j,idata), j=noa+nva+nob+1,noa+nva+nob+nvb )
          else
             write(19,500) ( popdata(j,idata), j=1,noa )
             write(19,600) ( popdata(j,idata), j=noa+1,noa+nva )
          end if
       end do
       close(19)
       write(iout,100) 'MO populations written to               ',trim(popfile)
    end if

101 format( '   occ_a:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
200 format( '   occ_b:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
300 format( '   vir_a:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
400 format( '   vir_b:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
500 format( '     occ:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
600 format( '     vir:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )

    if ( writeion ) then
       open( 20,file=trim(ionfile) )
       do idata = 1, ndata
          write(20,"(' TIME(fs)=',f9.2,5x,'NORM=',f13.10,5x,'RATE(1/fs)=',f13.10,5x,'CHECK(1/fs)=',f13.10)") &
               (outstep)*dble(idata)*dt*au2fs, normdata(idata), ratedata(idata)/au2fs, sum(iondata(:,idata))/au2fs
          if ( unrestricted ) then
             write(20,101) ( iondata(j,idata)/au2fs, j=1,noa )
             write(20,300) ( iondata(j,idata)/au2fs, j=noa+1,noa+nva )
             write(20,200) ( iondata(j,idata)/au2fs, j=noa+nva+1,noa+nva+nob )
             write(20,400) ( iondata(j,idata)/au2fs, j=noa+nva+nob+1,noa+nva+nob+nvb )
          else
             write(20,500) ( iondata(j,idata)/au2fs, j=1,noa )
             write(20,600) ( iondata(j,idata)/au2fs, j=noa+1,noa+nva )
          end if
       end do
       close(20)
       write(iout,100) 'Instantaneous rate & ion pop written to ',trim(ionfile)
       write(iout,'(A)') '                                integrated MO populations for the ion'
    end if


    if ( writemu ) then
       open( 28,file=trim(dipolefile) )
       write(28,"(a9,3(1x,a12))") 'time(fs)','mu_x (au)','mu_y (au)','mu_z (au)'
       write(28,"(f9.2,3(1x,f13.10))") 0.d0, tdx, tdy, tdz
       do idata = 1, ndata
          write(28,"(f9.2,3(1x,f13.10))")           &
               dble(outstep)*dble(idata)*dt*au2fs, &
               mudata(1,idata), mudata(2,idata), mudata(3,idata)
       end do
       close(28)
       write(iout,100) 'Dipole moments written to               ',trim(dipolefile)
    end if
    

100 format(12x,a39,a18)    


    flush(iout)


  end subroutine write_data
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine pop_rate(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,nva95,nva99,hole,part,state_ip_index,ip_states, &
      pop,ion,ion_coeff,rate_a,rate_b,psi_det,psiV,normV,vabs_a,vabs_b,unrestricted,density,au2fs,nva_rate)

    implicit none

    class(MolInfo), intent(in)       :: Mol
    integer(8), intent(in)    :: nbasis,iout, noa, nob, norb, nstates, nva, nvb, ip_states
    integer(8), intent(in)    :: hole(nstates,2), part(nstates), state_ip_index(noa+nob,noa+nob)
    real(8), intent(in)       :: au2fs, vabs_a(noa+nva,noa+nva), vabs_b(nob+nvb,nob+nvb)
    real(8), intent(inout)    :: pop(norb), ion(norb)
    real(8), allocatable, intent(inout) :: rate_a(:), rate_b(:)
    real(8), intent(inout)    :: normV, density(norb*norb)
    complex(8),intent(inout)  :: psiv(nstates),ion_coeff(ip_states)
    complex(8), intent(in)    :: psi_det(nstates)
    logical, intent(in)       :: unrestricted
    integer(8),intent(inout):: nva95, nva99
    real(8), intent(inout), optional :: nva_rate

    integer(8) :: i, i1, j, j1, a, a1, b, b1, ia, jb, ii, jj, aa, bb, istate, ndim, k
    integer(8) :: ij, iselect
    real(8)    :: const, const2, psi2, vabs00_, rdum, rate
    real(8)    :: density_AO(nbasis*nbasis)
    real(8)    :: tmprate(80000), tmpsum, tmp2, tmprate_direct
    integer(8) :: nva95_direct, nva99_direct, nva95_density, nva99_density
    complex(8) :: psi0, psi_ia
    logical    :: total_dens
    
     total_dens = .true.
     ndim = norb
     if( total_dens ) ndim = noa + nva
     const  = dsqrt(2.d0) 
     const2 = 2.d0
     if ( unrestricted ) then
       const  = 1.d0
       const2 = 1.d0
     end if
    
     pop = 0.d0
     ion = 0.d0
     density(1:ndim*ndim) = 0.d0
     ion_coeff = dcmplx(0.d0,0.d0)
     rate = 0.d0
     rate_a = 0.d0
     rate_b = 0.d0
     iselect = 8

     !write(iout, *) "pop_rate ndim: ", ndim
       
!### code assumes psi_det is normalized ###

       call get_norm(normV,nstates,psi_det)
!:       write(iout,*) " norm-det =",normV
!:       psi_det = psi_det/normV

     !: ground state
!:     write(iout,*) " starting pop_rate"
     vabs00_ = 0.d0
     do i=1, noa
       pop(i) = const2
       ion(i) = const2
       density(i+(i-1)*ndim) = const2
       rate_a(i) = rate_a(i) + const2 * vabs_a(i,i)
       vabs00_ = vabs00_ + const2 * vabs_a(i,i)
     end do
     if ( unrestricted ) then
       do i=1, nob
         pop(noa+nva+i) = const2
         ion(noa+nva+i) = const2
         if( total_dens ) then
           density(i+(i-1)*ndim) = density(i+(i-1)*ndim) + const2
         else
           density(i+noa+nva+(i+noa+nva-1)*ndim) = const2
         end if
         rate_b(i) = rate_b(i) + const2 * vabs_b(i,i)
         vabs00_ = vabs00_ + const2 * vabs_b(i,i)
       end do
     end if
!:     write(iout,*)  " VAbs00",vabs00_
!:     write(iout,*) rate_a
!:     write(iout,*) rate_b
     
     psi0 = psi_det(1)
     psiV = dcmplx(0.d0,0.d0)
     rate = dconjg(psi_det(1))*psi_det(1)*vabs00_
     psiV(1) = psi_det(1) * vabs00_
              
     do ia=2, nstates
!:     write(iout,"(' pop ',i5,6f10.5)")ia,(density(i+(i-1)*ndim),i=1,noa+1)
          
       psi_ia = dconjg( psi_det(ia) )
          
       ii = hole(ia,1) 
       aa = part(ia) 

       i1 = ii + noa + nva
       a1 = aa + noa + nva + nob
       if ( ii.lt.0 ) i1 = -ii
       if ( aa.lt.0 ) a1 = -aa + noa

       psi2 = dble(dconjg(psi_ia)*psi_ia)

       pop(i1) = pop(i1) - psi2
       pop(a1) = pop(a1) + psi2
 
       i = abs(ii)
       a = aa+nob
       if(aa.lt.0) a = -aa+noa

!:     Indices for density
 
       i1 = abs(ii)
       a1 = abs(aa) + noa
       if( .not. total_dens ) then
         if ( ii.gt.0 ) i1 = ii + noa + nva
         if ( aa.gt.0 ) a1 = aa + noa + nva + nob
       end if
!:       write(iout,"(12I5)") ia,ii,aa,i1,a1

       psi2 = dble(psi_ia*psi_det(1))
  
       rdum = 0.d0
       if ( ii.lt.0 .and. aa.lt.0 ) then
         density(i1+(a1-1)*ndim) = density(i1+(a1-1)*ndim) + const * psi2
         density(a1+(i1-1)*ndim) = density(a1+(i1-1)*ndim) + const * psi2
         rdum =  const * vabs_a(a,i)
         rate_a(i) = rate_a(i) + 2.d0 * psi2 * rdum
         if(iselect.eq.0.or.iselect.eq.i) rate_a(a) = rate_a(a) + 2.d0 * psi2 * rdum
         rate = rate + 2.d0 * psi2 * rdum
         psiV(1) = psiV(1) + psi_det(ia) * rdum
         psiV(ia) = psiV(ia) + psi_det(1) * rdum
       end if
       if ( ii.gt.0 .and. aa.gt.0 ) then
         density(i1+(a1-1)*ndim) = density(i1+(a1-1)*ndim) + const * psi2
         density(a1+(i1-1)*ndim) = density(a1+(i1-1)*ndim) + const * psi2
         rdum =  const * vabs_b(a,i)
         rate = rate + 2.d0 * psi2 * rdum
         rate_b(i) = rate_b(i) + 2.d0 * psi2 * rdum
         if(iselect.eq.0.or.iselect.eq.i) rate_b(a) = rate_b(a) + 2.d0 * psi2 * rdum
         psiV(1) = psiV(1) + psi_det(ia) * rdum
         psiV(ia) = psiV(ia) + psi_det(1) * rdum
       end if

       do jb=2, nstates
             
         jj = hole(jb,1)  
         bb = part(jb) 

         j = abs(jj)
         b = bb+nob
         if(bb.lt.0) b = -bb+noa
 
!:     Indices for density
 
         j1 = abs(jj)
         b1 = abs(bb) + noa
         if( .not. total_dens ) then
           if ( jj.gt.0 ) j1 = jj + noa + nva
           if ( bb.gt.0 ) b1 = bb + noa + nva + nob
         end if
!:         write(iout,"(12I5)") ia,ii,aa,jb,jj,bb,i1,a1,j1,b1
         psi2 = dble(psi_ia*psi_det(jb))
         if(aa.eq.bb .and. ii*jj.gt.0) density(j1+(i1-1)*ndim) = density(j1+(i1-1)*ndim) - psi2
         if(ii.eq.jj .and. aa*bb.gt.0) density(a1+(b1-1)*ndim) = density(a1+(b1-1)*ndim) + psi2

         rdum = 0.d0
             
         if ( ii.eq.jj .and. ii.lt.0 ) then
           if ( aa.lt.0 .and. bb.lt.0 ) then
             rdum = vabs_a(a,b)
             if(iselect.eq.0.or.iselect.eq.i) rate_a(a) = rate_a(a) + psi2 * rdum
           end if
           if ( aa.gt.0 .and. bb.gt.0 ) then
             rdum = vabs_b(a,b)
             if(iselect.eq.0.or.iselect.eq.i) rate_b(a) = rate_b(a) + psi2 * rdum
           end if
           rate = rate + psi2 * rdum
           rate_a(i) = rate_a(i) + psi2 * rdum
           psiV(ia) = psiV(ia) + psi_det(jb) * rdum
         end if
             
         if ( ii.eq.jj .and. ii.gt.0 ) then
           if ( aa.lt.0 .and. bb.lt.0 ) then
             rdum = vabs_a(a,b)
             if(iselect.eq.0.or.iselect.eq.i) rate_a(a) = rate_a(a) + psi2 * rdum
           end if
           if ( aa.gt.0 .and. bb.gt.0 ) then
             rdum = vabs_b(a,b)
             if(iselect.eq.0.or.iselect.eq.i) rate_b(a) = rate_b(a) + psi2 * rdum
           end if
           rate = rate + psi2 * rdum
           rate_b(i) = rate_b(i) + psi2 * rdum
           psiV(ia) = psiV(ia) + psi_det(jb) * rdum
         end if
             
         if ( aa.eq.bb ) then
           if ( ii.lt.0 .and. jj.lt.0 ) then
             rdum = - vabs_a(j,i)
             rate = rate + psi2 * rdum
             rate_a(i) = rate_a(i) + psi2 * rdum
             if(iselect.eq.0.or.iselect.eq.i) rate_a(a) = rate_a(a) + psi2 * rdum
             psiV(ia) = psiV(ia) + psi_det(jb) * rdum
           end if
           if ( ii.gt.0 .and. jj.gt.0 ) then
             rdum = - vabs_b(j,i)
             rate = rate + psi2 * rdum
             rate_b(i) = rate_b(i) + psi2 * rdum
             if(iselect.eq.0.or.iselect.eq.i) rate_b(a) = rate_b(a) + psi2 * rdum
             psiV(ia) = psiV(ia) + psi_det(jb) * rdum
           end if
         end if

         if( ia.eq.jb ) then
           rdum = vabs00_
           if ( ii.lt.0 .and. jj.lt.0 ) then
             if(iselect.eq.0.or.iselect.eq.i) rate_a(a) = rate_a(a) + psi2 * rdum
           end if
           if ( ii.gt.0 .and. jj.gt.0 ) then
             if(iselect.eq.0.or.iselect.eq.i) rate_b(a) = rate_b(a) + psi2 * rdum
           end if
           rate = rate + psi2 * rdum
           psiV(ia) = psiV(ia) + psi_det(jb) * rdum
         end if
!:       write(iout,998) ia,jb,ii,aa,jj,bb,rate_a,rate_b
998    format(6I4/4F20.15/4F20.15)
       end do
     end do
      rate_b(noa+nva+1) = dconjg(psi_det(1))*psi_det(1)*vabs00_
      rate_b(noa+nva+2) = rate
      !write(iout,*) " rate0 = ", rate
      rate = 0.d0
      do ia = 1,nstates
        rate = rate + dconjg(psi_det(ia))*psiV(ia)
      end do
      !write(iout,*) " rate1 = ", rate


       rate = 0.d0
       do i = 1,ndim
         do j = 1,ndim
           rate = rate + density(i+(j-1)*ndim) * vabs_a(i,j)
         end do
       end do
      !write(iout,*) " rate2 = ", rate

      !: Check rate in AO basis
      !call mo2ao_full(nbasis, ndim, density, density_AO, Mol%cmo_a)
      !rate = 0.d0
      !do i = 1,nbasis
      !  do j = 1,nbasis
      !    ij  = i*(i-1)/2 + j ! vabsao has triangular indexing
      !    if (j.gt.i) ij = j*(j-1)/2 + i
      !    rate = rate + density_AO(i+(j-1)*nbasis) * Mol%vabsao(ij)
      !  end do
      !end do
      !write(iout,*) " rate3 = ", rate

       


     
!:     Population of absorbed wavefunction
 
!:       write(iout,*) "Populations"
!:       write(iout,"(10f8.4)") (pop(i),i=1,noa+nva+nob+nvb)
!: ********** HBS
       tmprate = 0.D0
       do ia = 2,nstates
         aa = IAbs(part(ia))
         if(aa.lt.1.or.aa.gt.nva) write(iout,*) " out of range ",aa
         tmprate(aa)=tmprate(aa)+abs(dconjg(psi_det(ia))*psiV(ia))
       end do
       tmpsum=0.d0
       do aa = 1,nva
         tmpsum=tmpsum+tmprate(aa)
       end do 
       tmprate_direct = tmpsum
!:       write(iout,8888) (tmprate(aa)/tmpsum,aa=1,nva)
       tmprate(1) = tmprate(1)/tmpsum
       do aa = 2,nva
         tmprate(aa) = tmprate(aa-1)+tmprate(aa)/tmpsum
         if(tmprate(aa).lt.0.95d0) nva95 = aa
         if(tmprate(aa).lt.0.99d0) nva99 = aa
       end do
       nva95_direct = nva95
       nva99_direct = nva99
       if (present(nva_rate)) nva_rate = tmpsum
       
       !write(iout,8889) rate,nva,nva95,nva99
       nva95 = 0
       nva99 = 0
       tmprate = 0.D0
       tmprate(1) = abs(density(1)*vabs_a(1,1))
       do i = 2,noa+nva
         tmprate(i) = tmprate(i-1) + abs(density(i+(i-1)*ndim)*vabs_a(i,i))
         do j = 1,i-1
           tmprate(i) = tmprate(i) + abs(density(i+(j-1)*ndim)*vabs_a(i,j))
           tmprate(i) = tmprate(i) + abs(density(j+(i-1)*ndim)*vabs_a(j,i))
         end do
         !write(iout,*) "tmprate ",i,tmprate(i)
       end do
       do aa = 1,nva
         if(tmprate(aa).lt.0.95d0*tmprate(noa+nva)) nva95 = aa
         if(tmprate(aa).lt.0.99d0*tmprate(noa+nva)) nva99 = aa
       end do
       nva95_density = nva95
       nva99_density = nva95

       !: Make sure RESULTS files have eq12/direct ci vector method.
       nva95 = nva95_direct
       nva99 = nva99_direct

       !write(iout,8889) tmprate(noa+nva),nva,nva95,nva99
 8888  format(10F10.7)
 8889  format(" (direct) rate= ",F10.7,"  nva= ",I5,"  nva95= ",I5,"  nva99= ",I5)

       !if ( abs(tmprate_direct - tmprate(nva)) .gt. 1e-6 ) then
       !  write(iout, *) "Discrepency detected between direct and density derived rates: "
       !  write(iout, "(A, F10.7, A, I5, A, I5)") " (direct ) rate= ",tmprate_direct,", nva95= ", &
       !              nva95_direct, ", nva99= ", nva99_direct
       !  write(iout, "(A, F10.7, A, I5, A, I5)") " (density) rate= ",tmprate(nva),", nva95= ", &
       !              nva95, ", nva99= ", nva99
       !end if
 

       flush(iout)
!: ********** HBS
       call get_norm(normV,nstates,psiV)
       psiV = psiV/normV
!:      write(iout,*) " normV = ", normV

       do ia = 2, nstates

        psi2 = dble( dconjg(psiV(ia)) * psiV(ia) )

       ii = hole(ia,1) 
       aa = part(ia) 

       i = ii + noa + nva
       a = aa + noa + nva + nob
       if ( ii.lt.0 ) i = -ii
       if ( aa.lt.0 ) a = -aa + noa

        ion(i) = ion(i) - psi2
        ion(a) = ion(a) + psi2
      end do
 
!:       write(iout,*) "Ion Populations"
!:       write(iout,"(10f8.4)") (ion(i),i=1,noa+nva+nob+nvb)
 
     rate_a = rate_a * 2.d0 / au2fs
     rate_b = rate_b * 2.d0 / au2fs
  
  end subroutine pop_rate
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine pop_rate_ip(Mol,nbasis,iout,noa,nob,norb,nstates,nva,nvb,nva95,nva99,hole,part,state_ip_index,ip_states, &
      pop,ion,ion_coeff,rate_aa,rate_ab,rate_ba,rate_bb,psi_det,psiV,normV,vabs_a,vabs_b,density,au2fs,nva_rate)

    implicit none

    class(MolInfo), intent(in)     :: Mol
    integer(8), intent(in)  :: nbasis,iout,noa, nob, norb, nstates, nva, nvb, ip_states
    integer(8), intent(in)  :: hole(nstates,2), part(nstates), state_ip_index(noa+nob,noa+nob)
    real(8), intent(in)     :: au2fs, vabs_a(noa+nva,noa+nva), vabs_b(nob+nvb,nob+nvb)
    real(8), intent(inout)  :: pop(norb),ion(norb),rate_aa(noa,noa),rate_ab(noa,nob), &
                               rate_ba(nob,noa),rate_bb(nob,nob),normV,density(norb*norb)
    complex(8),intent(inout):: psiv(nstates),ion_coeff(ip_states)
    complex(8), intent(in)  :: psi_det(nstates)
    integer(8),intent(inout):: nva95, nva99
    real(8), intent(inout), optional :: nva_rate

    integer(8) :: i,i1,j,j1,a,a1,b,b1,ia,jb,ii,jj,x,x1,y,y1,aa,bb,xx,yy,istate,ndim
    real(8)    :: psi2, vabs00_, rate, rdum
    real(8)    :: tmprate(1000), tmpsum, tmprate_direct
    integer(8) :: nva95_direct, nva99_direct, nva95_density, nva99_density
    complex(8) :: psi_ia
    logical    :: total_dens
    
    total_dens = .true.
    ndim = norb
    if( total_dens ) ndim = noa + nva
    pop = 0.d0
    ion = 0.d0
    density(1:ndim*ndim) = 0.d0
    ion_coeff = dcmplx(0.d0,0.d0)
    rate = 0.d0
    rate_aa = 0.d0
    rate_ab = 0.d0           
    rate_ba = 0.d0
    rate_bb = 0.d0           


!### code assumes psi_det is normalized ###
  
     !: unrestricted ground state
 
     vabs00_ = 0.d0
     do i=1, noa
       pop(i) = 1.d0
       ion(i) = 1.d0
       density(i+(i-1)*ndim) = 1.d0
       vabs00_ = vabs00_ + vabs_a(i,i)
     end do
     do i=1, nob
       pop(noa+nva+i) = 1.d0
       ion(noa+nva+i) = 1.d0
       if( total_dens ) then
         density(i+(i-1)*ndim) = density(i+(i-1)*ndim) + 1.d0
       else
         density(i+noa+nva+(i+noa+nva-1)*ndim) = 1.d0
       end if
       vabs00_ = vabs00_ + vabs_b(i,i)
     end do

     psiV = 0.d0
     do ia=1, nstates

       psi_ia = dconjg(psi_det(ia))
       psi2 = dble( dconjg(psi_ia) * psi_ia )

       xx = hole(ia,1) 
       ii = hole(ia,2) 
       aa = part(ia)  
 
       x = abs(xx)
       i = abs(ii)
       a = abs(aa) + noa
       if( aa.gt.0 ) a = aa + nob

       x1 = x
       i1 = i
       a1 = a
       if ( xx.gt.0 ) x1 = xx + noa + nva
       if ( ii.gt.0 ) i1 = ii + noa + nva
       if ( aa.gt.0 ) a1 = aa + noa + nva + nob

!:       if(psi2.gt.1.d-7) write(iout,"(6i5,4F12.6)") xx,ii,aa,x1,i1,a1,psi_det(ia),psi2
 
!:     Population of wavefunction

       if ( ii.eq.0 ) then
         If(x1.ne.0) pop(x1) = pop(x1) - psi2
       else
         If(x1.ne.0) pop(x1) = pop(x1) - psi2
         pop(i1) = pop(i1) - psi2
         pop(a1) = pop(a1) + psi2
       end if 

!:     Indices for density
 
       x1 = abs(xx)
       i1 = abs(ii)
       a1 = abs(aa) + noa
       if( .not. total_dens ) then
         if ( xx.gt.0 ) x1 = xx + noa + nva
         if ( ii.gt.0 ) i1 = ii + noa + nva
         if ( aa.gt.0 ) a1 = aa + noa + nva + nob
       end if

         do jb=1, nstates
             
           yy = hole(jb,1)
           jj = hole(jb,2) 
           bb = part(jb)   

           y = abs(yy)
           j = abs(jj)
           b = abs(bb) + noa
           if (bb.gt.0) b = bb + nob

!:     Indices for density
 
           y1 = abs(yy)
           j1 = abs(jj)
           b1 = abs(bb) + noa
           if( .not. total_dens ) then
             if ( yy.gt.0 ) y1 = yy + noa + nva
             if ( jj.gt.0 ) j1 = jj + noa + nva
             if ( bb.gt.0 ) b1 = bb + noa + nva + nob
           end if

           psi2 = dble(psi_ia*psi_det(jb))

           rdum = 0.d0


          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
              if(xx*yy.gt.0) density(y1+(x1-1)*ndim) = density(y1+(x1-1)*ndim) - psi2
              !: < x |A| y >
              if ( xx.lt.0 .and. yy.lt.0 ) rdum = - vabs_a(y,x)
              !: < X |A| Y >
              if ( xx.gt.0 .and. yy.gt.0 ) rdum = - vabs_b(Y,X)
             go to 78
          end if SS


          !: singles doubles
          
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 ) then             
             if(ii*aa.gt.0) density(a1+(i1-1)*ndim) = density(a1+(i1-1)*ndim) + psi2
             !: < i->a , x |A| x >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) rdum = vabs_a(a,i) 
             !: < I->A, x |A| x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) rdum = vabs_b(A,I) 
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 ) then
             if(jj*bb.gt.0) density(b1+(j1-1)*ndim) = density(b1+(j1-1)*ndim) + psi2
             !: < x |A| j->b, x >  =  < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) rdum = vabs_a(b,j)
             !: < x |A| J->B, x >  =  < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) rdum = vabs_b(b,j)
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
              if(xx*aa.gt.0) density(a1+(x1-1)*ndim) = density(a1+(x1-1)*ndim) - psi2
              !: < y->a, x |A| y >  =  - < a |A| x >
              if ( xx.lt.0 .and. aa.lt.0 ) rdum = - vabs_a(a,x)
              !: < Y->A, X |A| Y >  =  - < A |A| X >
              if ( xx.gt.0 .and. aa.gt.0 ) rdum = - vabs_b(A,X)
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
              if(yy*bb.gt.0) density(b1+(y1-1)*ndim) = density(b1+(y1-1)*ndim) - psi2
              !: < x |A| x->b, y >  =  - < y | A | b >
              if ( yy.lt.0 .and. bb.lt.0 ) rdum = - vabs_a(b,y)
              !: < X |A| X->B, Y >  =  - < Y |A| B >
              if ( yy.gt.0 .and. bb.gt.0 ) rdum = - vabs_b(B,Y)
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             if(ii*jj.gt.0) density(i1+(j1-1)*ndim) = density(i1+(j1-1)*ndim) - psi2
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) rdum = - vabs_a(i,j)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( ii.gt.0 .and. jj.gt.0 ) rdum = - vabs_b(I,J)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             if(aa*bb.gt.0) density(a1+(b1-1)*ndim) = density(a1+(b1-1)*ndim) + psi2
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum + vabs_a(a,b)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum + vabs_b(A,B)
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             if(xx*yy.gt.0) density(x1+(y1-1)*ndim) = density(x1+(y1-1)*ndim) - psi2
             !: < i->a, x |A| i->a, y >  =  - < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) rdum = rdum - vabs_a(x,y)
             !: < I->A, X |A| I->A, Y >  =  - < Y |A| X >
             if ( xx.gt.0 .and. yy.gt.0 ) rdum = rdum - vabs_b(X,Y)
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             if(ii*yy.gt.0) density(i1+(y1-1)*ndim) = density(i1+(y1-1)*ndim) + psi2
             !: < i->a, x |A| x->a, y >  =  < y |A| i >
             if ( ii.lt.0 .and. yy.lt.0 ) rdum = rdum + vabs_a(i,y)
             !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
             if ( ii.gt.0 .and. yy.gt.0 ) rdum = rdum + vabs_b(I,Y) 
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             if(jj*xx.gt.0) density(j1+(x1-1)*ndim) = density(j1+(x1-1)*ndim) + psi2
             !: < y->a, x |A| j->a, y >  =  < j |A| x >
             if ( jj.lt.0 .and. xx.lt.0 ) rdum = rdum + vabs_a(j,x)
             !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
             if ( jj.gt.0 .and. xx.gt.0 ) rdum = rdum + vabs_b(J,X)
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             if(aa*bb.gt.0) density(a1+(b1-1)*ndim) = density(a1+(b1-1)*ndim) - psi2
             !: < y->a, x |A| x->b, y >  =  - < a |A| b >
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum - vabs_a(a,b)
             !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum - vabs_b(A,B)
          end if kdelta_iy_xj
          
78        continue
 
          if( ia.eq.jb ) rdum = rdum + vabs00_
          psiV(ia) = psiV(ia) + rdum * psi_det(jb)

995       format(i5,4e20.8)
996       format(2i5,4e20.8)
997       format(4e20.8)
998       format(8I4,6f20.14)
!:          write(iout,998) ia,jb,ii,xx,aa,jj,yy,bb,rdum
!:          write(iout,*) psi2,psiV(ia)
!:          flush(iout)
          rate = rate + psi2 * rdum
          if( ii.ne.0 ) then 
            if( ii.lt.0 .and. xx.lt.0 ) rate_aa(-ii,-xx) = rate_aa(-ii,-xx) + psi2 * rdum
            if( ii.lt.0 .and. xx.gt.0 ) rate_ab(-ii, xx) = rate_ab(-ii, xx) + psi2 * rdum
            if( ii.gt.0 .and. xx.lt.0 ) rate_ba( ii,-xx) = rate_ba( ii,-xx) + psi2 * rdum
            if( ii.gt.0 .and. xx.gt.0 ) rate_bb( ii, xx) = rate_bb( ii, xx) + psi2 * rdum
          end if
          if( ii.eq.0 .and. jj.ne.0 ) then 
            if( jj.lt.0 .and. yy.lt.0 ) rate_aa(-jj,-yy) = rate_aa(-jj,-yy) + psi2 * rdum
            if( jj.lt.0 .and. yy.gt.0 ) rate_ab(-jj, yy) = rate_ab(-jj, yy) + psi2 * rdum
            if( jj.gt.0 .and. yy.lt.0 ) rate_ba( jj,-yy) = rate_ba( jj,-yy) + psi2 * rdum
            if( jj.gt.0 .and. yy.gt.0 ) rate_bb( jj, yy) = rate_bb( jj, yy) + psi2 * rdum
          end if
          if( ii.eq.0 .and. jj.eq.0 ) then 
            if( xx.lt.0 .and. yy.lt.0 .and. -xx.ge.-yy) rate_aa(-xx,-yy) = rate_aa(-xx,-yy) + psi2 * rdum
            if( xx.lt.0 .and. yy.lt.0 .and. -xx.lt.-yy) rate_aa(-yy,-xx) = rate_aa(-yy,-xx) + psi2 * rdum
            if( xx.lt.0 .and. yy.gt.0 .and. -xx.ge. yy) rate_ab(-xx, yy) = rate_ab(-xx, yy) + psi2 * rdum
            if( xx.lt.0 .and. yy.gt.0 .and. -xx.lt. yy) rate_ab( yy,-xx) = rate_ab( yy,-xx) + psi2 * rdum
            if( xx.gt.0 .and. yy.lt.0 .and.  xx.ge.-yy) rate_ba( xx,-yy) = rate_ba( xx,-yy) + psi2 * rdum
            if( xx.gt.0 .and. yy.lt.0 .and.  xx.lt.-yy) rate_ba(-yy, xx) = rate_ba(-yy, xx) + psi2 * rdum
            if( xx.gt.0 .and. yy.gt.0 .and.  xx.ge. yy) rate_bb( xx, yy) = rate_bb( xx, yy) + psi2 * rdum
            if( xx.gt.0 .and. yy.gt.0 .and.  xx.lt. yy) rate_bb( yy, xx) = rate_bb( yy, xx) + psi2 * rdum
          end if
          end do
 
79        continue
 
       end do

!:      write(iout,*) " rate0 = ", rate
!:      rate = 0.d0
!:      do ia = 1,nstates
!:        rate = rate + dconjg(psi_det(ia))*psiV(ia)
!:      end do
!:      write(iout,*) " rate1 = ", rate
!:
!:       rate = 0.d0
!:       do i = 1,ndim
!:         do j = 1,ndim
!:           rate = rate + density(i+(j-1)*ndim) * vabs_a(i,j)
!:         end do
!:       end do
!:      write(iout,*) " rate2 = ", rate
!:     
!:       write(iout,*) "Populations"
!:       write(iout,"(10f8.4)") (pop(i),i=1,noa+nva+nob+nvb)
!:       write(iout,*) "Ion Populations"
!:       write(iout,"(10f8.4)") (ion(i),i=1,noa+nva+nob+nvb)

!:     Population of absorbed wavefunction

!: ********** HBS
       tmprate = 0.D0
       do ia = 1,nstates
         aa = IAbs(part(ia))
         if(aa.ge.1.and.aa.le.nva) then
           tmprate(aa)=tmprate(aa)+abs(dconjg(psi_det(ia))*psiV(ia))
           endif
       end do
       tmpsum=0.d0
       do aa = 1,nva
         tmpsum=tmpsum+tmprate(aa)
       end do
       tmprate_direct = tmpsum
!:       write(iout,8888) (tmprate(aa)/tmpsum,aa=1,nva)
!:       flush(iout)
       tmprate(1) = tmprate(1)/tmpsum
       do aa = 2,nva
         tmprate(aa) = tmprate(aa-1)+tmprate(aa)/tmpsum
         if(tmprate(aa).lt.0.95d0) nva95 = aa
         if(tmprate(aa).lt.0.99d0) nva99 = aa
       end do
       nva95_direct = nva95
       nva99_direct = nva99
       if (present(nva_rate)) nva_rate = tmpsum

       !write(iout,8889) rate,nva,nva95,nva99
       nva95 = 0
       nva99 = 0
       tmprate = 0.D0
       tmprate(1) = abs(density(1)*vabs_a(1,1))
       do i = 2,noa+nva
         tmprate(i) = tmprate(i-1) + abs(density(i+(i-1)*ndim)*vabs_a(i,i))
         do j = 1,i-1
           tmprate(i) = tmprate(i) + abs(density(i+(j-1)*ndim)*vabs_a(i,j))
           tmprate(i) = tmprate(i) + abs(density(j+(i-1)*ndim)*vabs_a(j,i))
         end do
       end do
       do aa = 1,nva
         if(tmprate(aa).lt.0.95d0*tmprate(noa+nva)) nva95 = aa
         if(tmprate(aa).lt.0.99d0*tmprate(noa+nva)) nva99 = aa
       end do
       nva95_density = nva95
       nva99_density = nva95

       !: Make sure RESULTS files have eq12/direct ci vector method.
       nva95 = nva95_direct
       nva99 = nva99_direct

!:       flush(iout)
 8887  format(2I5,F10.7)
 8888  format(10F10.7)
 8889  format(" (density) rate= ",F10.7,"  nva= ",I5,"  nva95= ",I5,"  nva99= ",I5)

       !if ( (abs(tmprate_direct-tmprate(nva)).gt.1e-6) .or. (nva99_direct .ne. nva99 )) then
       !  write(iout, *) "Discrepency detected between direct and density derived rates: "
       !  write(iout, "(A, F10.7, A, I5, A, I5)") " (direct ) rate= ",tmprate_direct,", nva95= ", &
       !              nva95_direct, ", nva99= ", nva99_direct
       !  write(iout, "(A, F10.7, A, I5, A, I5)") " (density) rate= ",tmprate(nva),", nva95= ", &
       !              nva95_density, ", nva99= ", nva99_density
       !end if

!: ********** HBS
 
       call get_norm(normV,nstates,psiV)
       psiV = psiV/normV

       do ia = 1, nstates

        xx = hole(ia,1) 
        ii = hole(ia,2) 
        aa = part(ia)  

        x = abs(xx)
        i = abs(ii)
        a = abs(aa) + noa
        if ( xx.gt.0 ) x = xx + noa + nva
        if ( ii.gt.0 ) i = ii + noa + nva
        if ( aa.gt.0 ) a = aa + noa + nva + nob

        x1 = x
        i1 = i
        if ( xx.gt.0 ) x1 = xx + noa
        if ( ii.gt.0 ) i1 = ii + noa

        psi2 = dble( dconjg(psiV(ia)) * psiV(ia) )

          if ( ii.eq.0 ) then
            If(x.ne.0) ion(x) = ion(x) - psi2
          else
            If(x.ne.0) ion(x) = ion(x) - psi2
            ion(i) = ion(i) - psi2
            ion(a) = ion(a) + psi2
          end if 
       end do
 
    rate_aa = rate_aa * 2.d0 / au2fs
    rate_ab = rate_ab * 2.d0 / au2fs
    rate_ba = rate_ba * 2.d0 / au2fs
    rate_bb = rate_bb * 2.d0 / au2fs
 
  end subroutine pop_rate_ip
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_ion_coeff_ip(iout,noa,nob,nva,nvb,nstates,hole,part,ip_states,state_ip_index, &
      psi_det,psiV,norm,vabs_a,vabs_b,rate,Zion_coeff,work,s)

    use util
    implicit none

    integer(8), intent(in)  :: iout, noa, nob, nva, nvb, nstates, ip_states, &
                               hole(nstates,2), part(nstates), state_ip_index(noa+nob,noa+nob)
    real(8), intent(in)     :: vabs_a(noa+nva,noa+nva), vabs_b(nob+nvb,nob+nvb)
    real(8), intent(inout)  :: norm, rate, s(nstates*nstates)
    complex(8), intent(in)  :: psi_det(nstates)
    complex(8),intent(inout):: psiv(nstates),Zion_coeff(ip_states*(ip_states+1)),work(2*ip_states*(nva+nvb))

    integer(8) :: i,i1,ix,j,j1,ij,a,a1,b,b1,ia,jb,ii,jj,x,x1,y,y1,aa,bb,xx,yy,istate,ndim,info
    real(8)    :: psi2, vabs00_, normV, rdum
!:    real(8)    :: rate_aa(noa,noa),rate_ab(noa,nob),rate_ba(nob,noa),rate_bb(nob,nob), &
!:                  pop(noa+nva+nob+nvb), ion(noa+nva+nob+nvb)

    complex(8) :: psi_ia, u(2), vt(2)
    
!:    do i=1,noa+nob
!:      write(iout,"(16i4)") (state_ip_index(i,j),j=1,noa+nob)
!:    end do
    Zion_coeff = dcmplx(0.d0,0.d0)
    rate = 0.d0
!:    pop = 0.d0
!:    ion = 0.d0
!:    rate_aa = 0.d0
!:    rate_ab = 0.d0           
!:    rate_ba = 0.d0
!:    rate_bb = 0.d0           


!### code assumes psi_det is normalized ###
  
     !: unrestricted ground state
 
     vabs00_ = 0.d0
     do i=1, noa
!:       pop(i) = 1.d0
!:       ion(i) = 1.d0
       vabs00_ = vabs00_ + vabs_a(i,i)
     end do
     do i=1, nob
!:       pop(noa+nva+i) = 1.d0
!:       ion(noa+nva+i) = 1.d0
       vabs00_ = vabs00_ + vabs_b(i,i)
     end do

     psiV = 0.d0
     do ia=1, nstates

       psi_ia = dconjg(psi_det(ia))
       psi2 = dble( dconjg(psi_ia) * psi_ia )

       xx = hole(ia,1) 
       ii = hole(ia,2) 
       aa = part(ia)  
 
       x = abs(xx)
       i = abs(ii)
       a = abs(aa) + noa
       if( aa.gt.0 ) a = aa + nob

       x1 = x
       i1 = i
       a1 = a
       if ( xx.gt.0 ) x1 = x + noa + nva
       if ( ii.gt.0 ) i1 = i + noa + nva
       if ( aa.gt.0 ) a1 = a + noa + nva + nob

!:       if(psi2.gt.1.d-7) write(iout,"(6i5,4F12.6)") xx,ii,aa,x1,i1,a1,psi_det(ia),psi2
!: 
!:     Population of wavefunction
!:
!:       if ( ii.eq.0 ) then
!:         If(x1.ne.0) pop(x1) = pop(x1) - psi2
!:       else
!:         If(x1.ne.0) pop(x1) = pop(x1) - psi2
!:         pop(i1) = pop(i1) - psi2
!:         pop(a1) = pop(a1) + psi2
!:       end if 

       x1 = x
       i1 = i
       if ( xx.gt.0 ) x1 = x + noa
       if ( ii.gt.0 ) i1 = i + noa

         do jb=1, nstates
             
           yy = hole(jb,1)
           jj = hole(jb,2) 
           bb = part(jb)   

           y = abs(yy)
           j = abs(jj)
           b = abs(bb) + noa
           if (bb.gt.0) b = bb + nob

           y1 = y
           j1 = j
           if ( yy.gt.0 ) y1 = y + noa
           if ( jj.gt.0 ) j1 = j + noa

           psi2 = dble(psi_ia*psi_det(jb))

           rdum = 0.d0


          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
              !: < x |A| y >
              if ( xx.lt.0 .and. yy.lt.0 ) rdum = - vabs_a(y,x)
              !: < X |A| Y >
              if ( xx.gt.0 .and. yy.gt.0 ) rdum = - vabs_b(Y,X)
             go to 78
          end if SS


          !: singles doubles
          
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 ) then             
             !: < i->a , x |A| x >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) rdum = vabs_a(a,i) 
             !: < I->A, x |A| x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) rdum = vabs_b(A,I) 
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 ) then
             !: < x |A| j->b, x >  =  < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) rdum = vabs_a(b,j)
             !: < x |A| J->B, x >  =  < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) rdum = vabs_b(b,j)
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
              !: < y->a, x |A| y >  =  - < a |A| x >
              if ( xx.lt.0 .and. aa.lt.0 ) rdum = - vabs_a(a,x)
              !: < Y->A, X |A| Y >  =  - < A |A| X >
              if ( xx.gt.0 .and. aa.gt.0 ) rdum = - vabs_b(A,X)
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
              !: < x |A| x->b, y >  =  - < y | A | b >
              if ( yy.lt.0 .and. bb.lt.0 ) rdum = - vabs_a(b,y)
              !: < X |A| X->B, Y >  =  - < Y |A| B >
              if ( yy.gt.0 .and. bb.gt.0 ) rdum = - vabs_b(B,Y)
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) rdum = - vabs_a(i,j)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( ii.gt.0 .and. jj.gt.0 ) rdum = - vabs_b(I,J)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum + vabs_a(a,b)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum + vabs_b(A,B)
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             !: < i->a, x |A| i->a, y >  =  - < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) rdum = rdum - vabs_a(x,y)
             !: < I->A, X |A| I->A, Y >  =  - < Y |A| X >
             if ( xx.gt.0 .and. yy.gt.0 ) rdum = rdum - vabs_b(X,Y)
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             !: < i->a, x |A| x->a, y >  =  < y |A| i >
             if ( ii.lt.0 .and. yy.lt.0 ) rdum = rdum + vabs_a(i,y)
             !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
             if ( ii.gt.0 .and. yy.gt.0 ) rdum = rdum + vabs_b(I,Y) 
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             !: < y->a, x |A| j->a, y >  =  < j |A| x >
             if ( jj.lt.0 .and. xx.lt.0 ) rdum = rdum + vabs_a(j,x)
             !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
             if ( jj.gt.0 .and. xx.gt.0 ) rdum = rdum + vabs_b(J,X)
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             !: < y->a, x |A| x->b, y >  =  - < a |A| b >
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum - vabs_a(a,b)
             !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum - vabs_b(A,B)
          end if kdelta_iy_xj
          
78        continue
 
          if( ia.eq.jb ) rdum = rdum + vabs00_
          psiV(ia) = psiV(ia) + rdum * psi_det(jb)

!:          write(iout,"(8i4,f20.14)") ia,jb,ii,xx,aa,jj,yy,bb,rdum

          rate = rate + psi2 * rdum
          if(ii.ne.0) istate = state_ip_index(i1,x1) 
          if(ii.eq.0.and.jj.ne.0) istate = state_ip_index(j1,y1) 
          if(ii.eq.0.and.jj.eq.0) istate = state_ip_index(x1,y1) 

!:          write(iout,*) psi2,psiV(ia)
!:          write(iout,"(5i5)") i1,x1,j1,y1,istate
!:          flush(iout)

          Zion_coeff(istate) = Zion_coeff(istate) + dcmplx(psi2*rdum,0.d0)

!:          if( ii.ne.0 ) then 
!:            if( ii.lt.0 .and. xx.lt.0 ) rate_aa(-ii,-xx) = rate_aa(-ii,-xx) + psi2 * rdum
!:            if( ii.lt.0 .and. xx.gt.0 ) rate_ab(-ii, xx) = rate_ab(-ii, xx) + psi2 * rdum
!:            if( ii.gt.0 .and. xx.lt.0 ) rate_ba( ii,-xx) = rate_ba( ii,-xx) + psi2 * rdum
!:            if( ii.gt.0 .and. xx.gt.0 ) rate_bb( ii, xx) = rate_bb( ii, xx) + psi2 * rdum
!:          end if
!:          if( ii.eq.0 .and. jj.ne.0 ) then 
!:            if( jj.lt.0 .and. yy.lt.0 ) rate_aa(-jj,-yy) = rate_aa(-jj,-yy) + psi2 * rdum
!:            if( jj.lt.0 .and. yy.gt.0 ) rate_ab(-jj, yy) = rate_ab(-jj, yy) + psi2 * rdum
!:            if( jj.gt.0 .and. yy.lt.0 ) rate_ba( jj,-yy) = rate_ba( jj,-yy) + psi2 * rdum
!:            if( jj.gt.0 .and. yy.gt.0 ) rate_bb( jj, yy) = rate_bb( jj, yy) + psi2 * rdum
!:          end if
!:          if( ii.eq.0 .and. jj.eq.0 ) then 
!:            if( xx.lt.0 .and. yy.lt.0 ) rate_aa(-xx,-yy) = rate_aa(-xx,-yy) + psi2 * rdum
!:            if( xx.lt.0 .and. yy.gt.0 ) rate_ab(-xx, yy) = rate_ab(-xx, yy) + psi2 * rdum
!:            if( xx.gt.0 .and. yy.lt.0 ) rate_ba( xx,-yy) = rate_ba( xx,-yy) + psi2 * rdum
!:            if( xx.gt.0 .and. yy.gt.0 ) rate_bb( xx, yy) = rate_bb( xx, yy) + psi2 * rdum
!:          end if
          end do
       end do

       call get_norm(normV,nstates,psiV)
       psiV = psiV/normV

!:       do i = 1, noa+nob
!:         do j = 1, noa+nob
!: fix for nactive .ne. 0
!:           ij = i*(i-1)/2 + j
!:           if(j.gt.i) ij = j*(j-1)/2 + i
!:           if(i.le.noa .and. j.le.noa) Zion_coeff(ij) = Zion_coeff(ij) = dcmplx(rate_aa(i,j),0.d0)
!:           if(i.le.noa .and. j.gt.noa) Zion_coeff(ij) = Zion_coeff(ij) = dcmplx(rate_ab(i,j),0.d0)
!:           if(i.gt.noa .and. j.le.noa) Zion_coeff(ij) = Zion_coeff(ij) = dcmplx(rate_ba(i,j),0.d0)
!:           if(i.gt.noa .and. j.gt.noa) Zion_coeff(ij) = Zion_coeff(ij) = dcmplx(rate_bb(i,j),0.d0)
!:           end do
!:         end do
!:
!:     Population of absorbed wavefunction
!: 
!:       do ia = 1, nstates
!:
!:        xx = hole(ia,1) 
!:        ii = hole(ia,2) 
!:        aa = part(ia)  
!:
!:        x1 = abs(xx)
!:        i1 = abs(ii)
!:        a1 = abs(aa) + noa
!:        if ( xx.gt.0 ) x1 = xx + noa + nva
!:        if ( ii.gt.0 ) i1 = ii + noa + nva
!:        if ( aa.gt.0 ) a1 = aa + noa + nva + nob
!:
!:        psi2 = dble( dconjg(psiV(ia)) * psiV(ia) )
!:
!:          if ( ii.eq.0 ) then
!:            If(x1.ne.0) ion(x1) = ion(x1) - psi2
!:          else
!:            If(x1.ne.0) ion(x1) = ion(x1) - psi2
!:            ion(i1) = ion(i1) - psi2
!:            ion(a1) = ion(a1) + psi2
!:          end if 
!:       end do

       ndim = ip_states*(nva+nvb)
       work(1:ndim) = dcmplx(0.d0,0.d0)
       do ia = 1, nstates

        xx = hole(ia,1) 
        ii = hole(ia,2) 
        aa = part(ia)  

        x1 = abs(xx)
        i1 = abs(ii)
        a1 = abs(aa)
        if ( xx.gt.0 ) x1 = xx + noa
        if ( ii.gt.0 ) i1 = ii + noa
        if ( aa.gt.0 ) a1 = aa + nva
 
        if(ii.ne.0) then
          istate = state_ip_index(x1,i1)
          if(istate.gt.0.and.istate.le.ip_states) work(istate+(a1-1)*ip_states) = psiV(ia)
!:          write(iout,"('a2',5i5,2f13.7)") ia,xx,ii,aa,istate,psiV(ia)
        end if
!:        ix = i*(i-1)/2 + x
!:        if ( x.gt.i ) ix = x*(x-1)/2 + i
!:         work(ix+(a-1)*ip_states) = psiV(ia)
       end do

! SVD of ion_coeff
  
       s = 0.d0 
!:       do i =1,ip_states
!:         write(iout,"('a2 ',i3,16F12.8)")i,(work(i+(j-1)*ip_states),j=1,nva+nvb)
!:       end do
       info = 10
       call zgesvd('O','N',ip_states,nva+nvb,work,ip_states, &
         s(ip_states+1:2*ip_states),u,2,vt,2,work(1+ndim:2*ndim), &
         ndim,s(1+2*ip_states:7*ip_states),info)
!:       write(iout,"('a2a0',28F12.8)")(s(i+ip_states),i=1,ip_states)
!:       do i =1,ip_states
!:         write(iout,"('a2 ',i3,16F12.8)")i,(work(i+(j-1)*ip_states),j=1,nva+nvb)
!:       end do
!:       return
!:       write(iout,"(' after return')")
       call Zfix_phase(ip_states,work)
! add code to fix ordering of degenerate eigenvectors

       rdum = 0.d0
       do i = 1, ip_states
         ii = (i-1)*ip_states
         jj = i*ip_states
         Zion_coeff(jj+1:jj+ip_states) = work(ii+1:ii+ip_states)
         s(i) = s(i+ip_states)
         rdum = rdum + s(i+ip_states)**2
       end do
!:       write(iout,"('a2b',2F12.8)") rate,rdum
!:       flush(iout)
 
  end subroutine get_ion_coeff_ip
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_ion_coeff(iout,noa,nob,nva,nvb,nstates,hole,part,psi_det,psiV,norm, &
                 vabs_a,vabs_b,unrestricted,rate,Zion_coeff,work,s)
    use util
    implicit none

    integer(8), intent(in)    :: iout, noa, nob, nva, nvb, nstates, &
                                 hole(nstates,2), part(nstates)
    real(8), intent(in)       :: vabs_a(noa+nva,noa+nva), vabs_b(nob+nvb,nob+nvb)
    real(8), intent(inout)    :: norm, rate, s(nstates*nstates)
    complex(8), intent(in)    :: psi_det(nstates),work(:)
    complex(8), intent(inout) :: psiv(nstates),Zion_coeff((noa+nob)*(nva+nvb+1))
    logical, intent(in)       :: unrestricted

    integer(8) :: i, j, a, b, ia, jb, ii, jj, aa, bb, info, noab, nvab
    integer(8) :: index(50)
    real(8)    :: const, const2, psi2, vabs00_, rdum, rdum1, norm1
    complex(8) :: psi_ia, u(2), vt(2)
    
!    get_psid produces psi_det that is normalized 
!    norm is the norm of psi_det before normalization
!    psiV = Vabs*psi_det is returned normalized and in the determinantal basis
!    norm is returned as the norm of psiV before normalization
!    ion_coeff are the coefficients of the ionic wavefunction obtained from normalized psiV by SVD
!    s are the weights of the ionic wavefunctions obtained by SVD (sum s**2 = 1)
!    ion_coeff is multiplied by unnormalized rate to yield the raw rate 
!      at which the ionic wavefunctions are produced by the simulation

     noab = noa + nob
     nvab = nva + nvb

     norm1 = 0.d0
     do ia = 1,nstates
       norm1 = norm1 + dble(dconjg(psi_det(ia))*psi_det(ia))
     end do
     s(1:noab) = 0.d0
     const  = dsqrt(2.d0)
     const2 = 2.d0
     if ( unrestricted ) then
       const  = 1.d0
       const2 = 1.d0
     end if
     !: ground state
     vabs00_ = 0.d0
     do i=1, noa
       vabs00_ = vabs00_ + const2 * vabs_a(i,i)
       s(i) = s(i) +  vabs_a(i,i)
     end do
     do i=1, nob
       if ( unrestricted ) then
         vabs00_ = vabs00_ + const2 * vabs_b(i,i)
         s(i+noa) = s(i+noa) + vabs_b(i,i)
       else
         s(i+noa) = s(i+noa) + vabs_a(i,i)
       end if
     end do

     psiV = dcmplx(0.d0,0.d0)
     psi2 = dreal(dconjg(psi_det(1))*psi_det(1))
     rate = psi2 * vabs00_
     psiV(1) = psi_det(1) * vabs00_
     rdum=0.d0
     do i=1,noa+nob
       rdum=rdum+s(i)
       s(i) = psi2 * s(i)
     end do
!:     write(iout,*) " vabs00_",vabs00_,rdum

     do ia=2, nstates

       ii = hole(ia,1)
       aa = part(ia)

       i =  abs(ii)
       a =  aa+nob
       if ( aa.lt.0 ) a = -aa+noa
       psi2 = dreal(dconjg(psi_det(ia))*psi_det(1)+dconjg(psi_det(1))*psi_det(ia))

       rdum = 0.d0
       if ( ii.lt.0 .and. aa.lt.0 ) then
         rdum =  const * vabs_a(a,i)
         s(i) = s(i) + psi2 * rdum
         rate = rate + psi2 * rdum
         psiV(1)  = psiV(1)  + psi_det(ia) * rdum
         psiV(ia) = psiV(ia) + psi_det(1)  * rdum
       end if
       if ( ii.gt.0 .and. aa.gt.0 ) then
         rdum =  const * vabs_b(a,i)
         s(i+noa) = s(i+noa) + psi2 * rdum
         rate = rate + psi2 * rdum
         psiV(1)  = psiV(1)  + psi_det(ia) * rdum
         psiV(ia) = psiV(ia) + psi_det(1)  * rdum
       end if

       do jb=2, nstates

         jj = hole(jb,1)
         bb = part(jb)
         j  = abs(jj)
         b  = bb+nob
         if(bb.lt.0) b = -bb+noa
         psi2 = dreal(dconjg(psi_det(ia))*psi_det(jb))

         if ( ii.eq.jj ) then
           rdum = 0.d0
           if ( aa.lt.0 .and. bb.lt.0 ) rdum = vabs_a(a,b)
           if ( aa.gt.0 .and. bb.gt.0 ) rdum = vabs_b(a,b)
           if ( ii.lt.0 ) then
              s(i) = s(i) + psi2 * rdum
           else
              s(i+noa) = s(i+noa) + psi2 * rdum
           end if
           rate = rate + psi2 * rdum
           psiV(ia) = psiV(ia) + psi_det(jb) * rdum
         end if

         if ( aa.eq.bb ) then
           rdum = 0.d0
           if ( ii.lt.0 .and. jj.lt.0 ) rdum = -vabs_a(j,i)
           if ( ii.gt.0 .and. jj.gt.0 ) rdum = -vabs_b(j,i)
           if ( ii.lt.0 ) then
              s(i) = s(i) + psi2 * rdum
           else
              s(i+noa) = s(i+noa) + psi2 * rdum
           end if
           rate = rate + psi2 * rdum
           psiV(ia) = psiV(ia) + psi_det(jb) * rdum
         end if
         if( ia.eq.jb ) then
           if ( ii.lt.0 ) then
              s(i) = s(i) + psi2 * vabs00_
           else
              s(i+noa) = s(i+noa) + psi2 * vabs00_
           end if
           rate = rate + psi2 * vabs00_
           psiV(ia) = psiV(ia) + psi_det(jb) * vabs00_
         end if
       end do
     end do

!:     write(iout,"('a1a',14f12.8)") norm1,norm,rate
     rate = 0.d0
     do ia = 1, nstates
       rate = rate +dble(dconjg(psi_det(ia))*psiV(ia))
     end do
     rate = 2.0d0*rate*norm**2
     s(1:noab) = 2.D0*s(1:noab)*norm**2
     rdum = 0.d0
     do i = 1, noab
       Zion_coeff(i) = dcmplx(s(i),0.d0)
       rdum = rdum + s(i)
     end do
!:     write(iout,"('a1b',14f12.8)") norm1,norm,0.5d0*rate/norm**2,0.5d0*rdum/norm**2,s(1:noab)

       call get_norm(norm1,nstates,psiV)
       psiV = psiV/norm1

       work(1:noab*nvab) = dcmplx(0.d0,0.d0)
       do ia = 2, nstates

         ii = hole(ia,1)
         aa = part(ia)

         i = ii + noa
         a = aa + nva
         if ( ii.lt.0 ) i = -ii
         if ( aa.lt.0 ) a = -aa

         work(i+(a-1)*noab) = psiV(ia)
       end do

! SVD of ion_coeff

       info = 10
       call zgesvd('O','N',noab,nvab,work,noab, &
         s(noab+1:2*noab),u,2,vt,2,work(1+noab*nvab:2*noab*nvab),noab*nvab,s(1+2*noab:7*noab),info)
!:       write(iout,"('a2a',17F12.8)")(s(i+noab),i=1,noab)
!:     fix phase and alpha beta ordering
       call Zfix_phase(noab,work)
       do i = 2, noab
         if(abs(s(i-1+noab)-s(i+noab)).lt.1.d-7) &
           call fix_order(noa,nob,noa+nob,work(1+(i-2)*noab),work(1+(i-1)*noab))
         end do
       do i = 2, noab
         if(abs(s(i-1+noab)-s(i+noab)).lt.1.d-7) &
           call fix_order(noa,nob,noa+nob,work(1+(i-2)*noab),work(1+(i-1)*noab))
         end do
       do i = 2, noab
         if(abs(s(i-1+noab)-s(i+noab)).lt.1.d-7) &
           call fix_order(noa,nob,noa+nob,work(1+(i-2)*noab),work(1+(i-1)*noab))
         end do

       rdum = 0.d0
       do i = 1, noab
         ii = (i-1)*noab
         jj = i*noab
         Zion_coeff(jj+1:jj+noab) = work(ii+1:ii+noab)
         s(i) = s(i+noab)
!:         write(iout,"(i3,17f10.6)") i,s(i),Zion_coeff(jj+1:jj+noab)
         rdum = rdum + s(i+noab)**2
       end do

!:       write(iout,"(i3,16f10.6)") 0,Zion_coeff(1:noab)
!:       do i = 1, noab
!:         ii = i*noab
!:         write(iout,"(i3,16f10.6)") i,Zion_coeff(ii+1:ii+noab)
!:       end do
!:       write(iout,"('a2b',17F12.8)") norm,0.5d0*rate/norm**2,norm1,rdum,(s(i),i=1,noab)

  end subroutine get_ion_coeff
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_ion_psi(iout,noa,nob,nva,nvb,nstates,nst_ion,hole,part, &
                 ion_sample_state,Zion_coeff,cis_vec,psi)

    implicit none

    integer(8), intent(in)    :: iout, noa, nob, nva, nvb, nstates, nst_ion, &
                                 hole(nstates,2), part(nstates), ion_sample_state
    real(8), intent(in)       :: cis_vec(nstates*nstates)
    complex(8), intent(in)    :: Zion_coeff((noa+nob)*(nva+nvb))
    complex(8), intent(inout) :: psi(nstates)

    integer(8) :: i, j, k, x, a, ia, noab
    complex(8) :: cdum, zdum

     noab = noa + nob

!: assume the lowest nst_ion states include all the singly ionized ground states
     do k = 1, nst_ion
       cdum = dcmplx(0.d0,0.d0)
       do ia = 1, nstates
         a = part(ia)
         if( a.eq.0 ) then
           x = -hole(ia,1)
           i = -hole(ia,2)
           if( hole(ia,1).gt.0 ) x = hole(ia,1) + noa
             if( ion_sample_state.gt.0 ) then
               zdum = Zion_coeff(x+(ion_sample_state-1)*noab)
             else
               !: for ion_sample_state = 0, sum all of the ion states
               zdum = dcmplx(0.d0,0.d0)
               do j = 1, noab
                 zdum = zdum + Zion_coeff(x+(j-1)*noab)
               end do
             end if
             cdum = cdum + cis_vec(ia+(k-1)*nstates)*zdum
!:             write(iout,"(4i5,4f10.6)") ia,x,i,a,zdum,cis_vec(ia+(k-1)*nstates)
         end if
       end do
       write(iout,"(i5,9f12.7)") k,cdum,psi(k),psi(k)+cdum
!:       psi(k) = psi(k) + cdum
       psi(k) = psi(k) + cabs(cdum) * psi(k)/cabs(psi(k))
     end do
!:     write(iout,"(9f12.7)") (psi(k),k=1,noab)
  end subroutine get_ion_psi
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_Zion_psi(iout,noa,nob,nva,nvb,nstates,nst_ion,hole,part, &
                 ion_sample_state,Zion_coeff,Zcis_vec,psi)

    implicit none

    integer(8), intent(in)    :: iout, noa, nob, nva, nvb, nstates, nst_ion, &
                                 hole(nstates,2), part(nstates), ion_sample_state
    complex(8), intent(in)    :: Zion_coeff((noa+nob)*(nva+nvb)),Zcis_vec(nstates*nstates)
    complex(8), intent(inout) :: psi(nstates)

    integer(8) :: i, j, k, x, a, ia, noab
    complex(8) :: cdum, zdum
    
     noab = noa + nob

!: assume the lowest nst_ion states include all the singly ionized ground states
     do k = 1, nst_ion
       cdum = dcmplx(0.d0,0.d0)
       do ia = 1, nstates
         a = part(ia)
         if( a.eq.0 ) then
           x = -hole(ia,1)
           i = -hole(ia,2)
           if( hole(ia,1).gt.0 ) x = hole(ia,1) + noa
             if( ion_sample_state.gt.0 ) then
               zdum = Zion_coeff(x+(ion_sample_state-1)*noab)
             else
               !: for ion_sample_state = 0, sum all of the ion states
               zdum = dcmplx(0.d0,0.d0)
               do j = 1, noab
                 zdum = zdum + Zion_coeff(x+(j-1)*noab)
               end do 
             end if
             cdum = cdum + Zcis_vec(ia+(k-1)*nstates)*zdum
!:             write(iout,"(4i5,4f10.6)") ia,x,i,a,zdum,cis_vec(ia+(k-1)*nstates)
         end if
       end do
!:       write(iout,"(i5,9f12.7)") k,cdum,psi(k),psi(k)+cdum*dt
       psi(k) = psi(k) + cdum
     end do
  end subroutine get_Zion_psi
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FIXPHASE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine fixphase(n,c)
 
    !: multiply eigenvectors so that the largest coefficient is real and positive  
    
    implicit none
    
    integer(8), intent(in) :: n
    complex(8), intent(inout) :: c(n)
 
    integer(8) :: i, imax
    real(8) :: cmax, cmax1
    complex(8) :: phase

    cmax = dble(dconjg(c(1))*c(1))
    imax = 1
    do i = 1, n
      cmax1 = dble(dconjg(c(i))*c(i))
      if(cmax1 .gt. cmax+1.d-9) then
        imax = i
        cmax = cmax1
      end if
    end do
    if(imax.gt.0) then
      phase = dconjg(c(imax))/dsqrt(cmax)
    else
      phase = dcmplx(1.d0,0.d0)
    end if
    do i = 1, n
      c(i) = phase*c(i)
    end do

  end subroutine fixphase 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_ion_psi1(iout,nstates,nstuse,nstion,isample, &
               dt,cis_vec,psi,Zion_coeff,psi_det)

    implicit none

    integer(8), intent(in)    :: iout, nstates, nstuse, nstion, isample
    real(8), intent(in)       :: dt,cis_vec(nstates*nstates)
    complex(8), intent(in)    :: Zion_coeff(nstion*(nstion+1)), psi_det(nstion)
    complex(8), intent(inout) :: psi(nstuse)

    integer(8) :: i, j, k, kk
    real(8)    :: rdum, acdum, rate, norm0, norm1, norm2, norm3
    complex(8) :: cdum, zdum, factor

!: transform psi to psi_det for the singly ionized configurations
     
   psi_det(1:nstates) = dcmplx(0.d0,0.d0)
!:   norm0 = 0.d0
   do k = 1, nstuse
!:     norm0 = norm0 + dble(dconjg(psi(k))*psi(k))
     kk = (k-1)*nstates
!::     psi_det(1:nstates) = psi_det(1:nstates) + cis_vec(kk+1:kk+nstates)*psi(k)
     psi_det(1:nstion) = psi_det(1:nstion) + cis_vec(kk+1:kk+nstion)*psi(k)
   end do
!:     norm1 = 0.d0
!:     do k = 1, nstion
!:       norm1 = norm1 + dble(dconjg(psi_det(k))*psi_det(k))
!:     end do
     
!: add contribution from Zion_coeff; the factor is chosen so that the change
!: in norm**2 is equal to the rate*dt; the phase is chosen to match psi

!:   write(iout,"(i3,16f12.8)") isample,(Zion_coeff(i+isample*nstion),i=1,nstion)
   if(isample.ne.0) then

!;   uses Zion_coeff = sqrt(rate) * s(isample) * u(isample) from SVD of absorbed wfn
     rate = 0.d0
     cdum = dcmplx(0.d0,0.d0)
     do i = 1, nstion
       zdum = Zion_coeff(i+isample*nstion)
       cdum = cdum + dconjg(psi_det(i)) * zdum
       rate = rate + dble(dconjg(zdum) * zdum)
     end do
     rate = rate * dt
     acdum = abs(cdum)
     if(acdum/rate.gt.1.d-7) then
       factor= (dsqrt((acdum/rate)**2+1.d0)-acdum/rate)*dconjg(cdum)/acdum
     else
       factor = dcmplx(1.d0,0.d0)
     end if
!:     write(iout,"('rate',8f16.10)") rate,factor,acdum/rate,zdum,cdum
!:     write(iout,"(17f12.8)") factor,(Zion_coeff(i+isample*nstion),i=1,nstion)
     do i = 1, nstion
!::       psi_det(i) = psi_det(i) + factor*Zion_coeff(i+isample*nstion)
       psi_det(i) = factor*Zion_coeff(i+isample*nstion)
     end do

   else

!:   uses Zion_coeff(i) = rate(i) from partitioning of the rate into orbital contributions
     rate = 0.d0
     do i = 1, nstion
       rate = rate + cabs(Zion_coeff(i)) * dt
       rdum = cabs(psi_det(i))
       if(rdum.gt.1.d-7) then
!::         psi_det(i) = psi_det(i) + (dsqrt(rdum**2+cabs(Zion_coeff(i))*dt)-rdum)*psi_det(i)/rdum
         psi_det(i) = (dsqrt(rdum**2+cabs(Zion_coeff(i))*dt)-rdum)*psi_det(i)/rdum
       else
!::         psi_det(i) = psi_det(i) + dsqrt(cabs(Zion_coeff(i))*dt)
         psi_det(i) = dsqrt(cabs(Zion_coeff(i))*dt)
       end if
     end do
!:     write(iout,"('cabs(Zion)',9f12.8)") (cabs(Zion_coeff(i)),i=1,nstion)
!:     write(iout,"('cabs(psid)',9f12.8)") (cabs(psi_det(i)),i=1,nstion)

   end if
!:     norm2 = 0.d0
!:     do k = 1, nstion
!:       norm2 = norm2 + dble(dconjg(psi_det(k))*psi_det(k))
!:     end do

!: transform psi_det to psi for the singly ionized configurations
    
!:   norm3 = 0.d0 
!:: do not zero psi
!::   psi(1:nstuse) = dcmplx(0.d0,0.d0)
   do k = 1, nstuse
     kk = (k-1)*nstates
     cdum = dcmplx(0.d0,0.d0)
!::     do i = 1, nstates
     do i = 1, nstion
       cdum = cdum +cis_vec(i+kk)*psi_det(i)
     end do
     psi(k) = psi(k) + cdum
!:     norm3 = norm3 + dble(dconjg(psi(k))*psi(k))
   end do
!:       write(iout,"('norm',8f16.12)") norm0,norm1,norm2,norm3,norm2-norm1,norm3-norm0,rate
!:     write(iout,"(9f12.7)") (psi(i),i=1,nstuse)

  end subroutine get_ion_psi1
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    subroutine get_Zion_psi1(iout,nstates,nstuse,nstion,isample, &
               dt,Zcis_vec,psi,Zion_coeff,psi_det)

    implicit none

    integer(8), intent(in)    :: iout, nstates, nstuse, nstion, isample
    real(8), intent(in)       :: dt
    complex(8), intent(in)    :: Zcis_vec(nstates*nstates), Zion_coeff(nstion*(nstion+1)), psi_det(nstion)
    complex(8), intent(inout) :: psi(nstuse)

    integer(8) :: i, j, k, kk
    real(8)    :: rdum, acdum, rate, norm0, norm1, norm2, norm3
    complex(8) :: cdum, zdum, factor

!: transform psi to psi_det for the singly ionized configurations
     
   psi_det(1:nstates) = dcmplx(0.d0,0.d0)
!:   norm0 = 0.d0
   do k = 1, nstuse
!:     norm0 = norm0 + dble(dconjg(psi(k))*psi(k))
     kk = (k-1)*nstates
     psi_det(1:nstion) = psi_det(1:nstion) + Zcis_vec(kk+1:kk+nstion)*psi(k)
   end do
!:     norm1 = 0.d0
!:     do k = 1, nstates
!:       norm1 = norm1 + dble(dconjg(psi_det(k))*psi_det(k))
!:     end do
     
!: add contribution from Zion_coeff; the factor is chosen so that the change
!: in norm**2 is equal to the rate*dt; the phase is chosen to match psi

!:   write(iout,"(i3,16f12.8)") isample,(Zion_coeff(i+isample*nstion),i=1,nstion)
   if(isample.ne.0) then

!;   uses Zion_coeff = sqrt(rate) * s(isample) * u(isample) from SVD of absorbed wfn
     rate = 0.d0
     cdum = dcmplx(0.d0,0.d0)
     do i = 1, nstion
       zdum = Zion_coeff(i+isample*nstion)
       cdum = cdum + dconjg(psi_det(i)) * zdum
       rate = rate + dble(dconjg(zdum) * zdum)
     end do
     rate = rate * dt
     acdum = abs(cdum)
     if(acdum/rate.gt.1.d-7) then
       factor= (dsqrt((acdum/rate)**2+1.d0)-acdum/rate)*dconjg(cdum)/acdum
     else
       factor = dcmplx(1.d0,0.d0)
     end if
!:     write(iout,"('rate',8f16.10)") rate,factor,acdum/rate,zdum,cdum
!:     write(iout,"(17f12.8)") factor,(Zion_coeff(i+isample*nstion),i=1,nstion)
     do i = 1, nstion
       psi_det(i) = factor*Zion_coeff(i+isample*nstion)
     end do

   else

!:   uses Zion_coeff(i) = rate(i) from partitioning of the rate into orbital contributions
     rate = 0.d0
     do i = 1, nstion
       rate = rate + cabs(Zion_coeff(i)) * dt
       rdum = cabs(psi_det(i))
       if(rdum.gt.1.d-7) then
         psi_det(i) = (dsqrt(rdum**2+cabs(Zion_coeff(i))*dt)-rdum)*psi_det(i)/rdum
       else
         psi_det(i) = dsqrt(cabs(Zion_coeff(i))*dt)
       end if
     end do
!:     write(iout,"('cabs(Zion)',9f12.8)") (cabs(Zion_coeff(i)),i=1,nstion)
!:     write(iout,"('cabs(psid)',9f12.8)") (cabs(psi_det(i)),i=1,nstion)

   end if
!:     norm2 = 0.d0
!:     do k = 1, nstates
!:       norm2 = norm2 + dble(dconjg(psi_det(k))*psi_det(k))
!:     end do

!: transform psi_det to psi for the singly ionized configurations
    
!:   norm3 = 0.d0 
   do k = 1, nstuse
     kk = (k-1)*nstates
     cdum = dcmplx(0.d0,0.d0)
     do i = 1, nstion
       cdum = cdum + dconjg(Zcis_vec(i+kk))*psi_det(i)
     end do
     psi(k) = psi(k) + cdum
!:     norm3 = norm3 + dble(dconjg(psi(k))*psi(k))
   end do
!:       write(iout,"('norm',8f16.12)") norm0,norm1,norm2,norm3,norm2-norm1,norm3-norm0,rate
!:     write(iout,"(9f12.7)") (psi(i),i=1,nstuse)

  end subroutine get_Zion_psi1




  subroutine ham_cis_to_mo(hole, part, psi_det, noa, nva, nstates, nrorb, ham_cis, ham_mo, nob, nvb)

      ! Arguments
      implicit none
      integer(8), intent(in) :: noa, nva, nstates, nrorb
      integer(8), intent(in) :: hole(nstates), part(nstates)
      complex(8), intent(in) :: psi_det(nstates)
      real(8), intent(in) :: ham_cis(nstates * nstates)  ! Hamiltonian in CIS state basis (1D array)
      real(8), intent(inout) :: ham_mo(nrorb * nrorb)      ! Hamiltonian in MO basis (1D array)
      integer(8), optional, intent(in) :: nob, nvb

      ! Local variables
      integer(8) :: i, j, a, b, ia, ib, ii, aa, jj, bb
      integer(8) :: cis_index, mo_index
      real(8) :: psi2_i, psi2_j, ham_contrib

      ! Initialize MO Hamiltonian
      ham_mo = 0.d0

      ! Loop over all CIS states for contributions
      do ia = 1, nstates
          ii = hole(ia)
          aa = part(ia)

          ! Only calculate for alpha spin orbitals.
          if ( ii.lt.0 .and. aa.lt.0 ) then
            psi2_i = dble(dconjg(psi_det(ia)) * psi_det(ia))  ! Norm of CIS coefficient
            i = ii
            a = aa + noa

            ! Loop over other states to compute matrix element contributions
            do ib = 1, nstates
                jj = hole(ib)
                bb = part(ib)

                ! Only calculate for alpha spin orbitals.
                if ( jj.lt.0 .and. bb.lt.0 ) then
                  psi2_j = dble(dconjg(psi_det(ib)) * psi_det(ib))  ! Norm of CIS coefficient
                  j = jj
                  b = bb + noa

                  cis_index = (ia - 1) * nstates + ib

                  ! Contribution to MO Hamiltonian
                  ham_contrib = ham_cis(cis_index) * psi2_i * psi2_j

                  ham_mo((i - 1) * nrorb + j) = ham_mo((i - 1) * nrorb + j) + ham_contrib
                  ham_mo((a - 1) * nrorb + b) = ham_mo((a - 1) * nrorb + b) + ham_contrib
                  ham_mo((i - 1) * nrorb + b) = ham_mo((i - 1) * nrorb + b) - ham_contrib
                  ham_mo((a - 1) * nrorb + j) = ham_mo((a - 1) * nrorb + j) - ham_contrib
                end if
            end do
          end if
      end do

  end subroutine ham_cis_to_mo

  !: Generate the complex reduced MO density
  subroutine make_complex_density(hole, part, psi_det, noa, nva, nstates, density_complex)
      ! Arguments
      implicit none
      integer(8), intent(in) :: noa, nva, nstates
      integer(8), intent(in) :: hole(nstates), part(nstates)
      complex(8), intent(in) :: psi_det(nstates)            ! CI vector in state basis
      complex(8), allocatable, intent(inout) :: density_complex(:)  ! Complex MO density matrix

      ! Local variables
      integer(8) :: i, j, i1, j1, a1,b1, ii, jj, aa, bb, ia, jb, idx
      complex(8) :: psi_ia, psi_jb
      integer(8) :: nrorb
      integer(8), parameter :: iout = 42

      nrorb = noa+nva

      ! Initialization
      density_complex = dcmplx(0.d0, 0.d0)

      ! Populate the density matrix
      !: Loop from 2 because 1 represents HF determinant
      !:   and causes segfault bc hole(1,1) is 0.
      do ia = 2, nstates
          psi_ia = psi_det(ia)

          ii = hole(ia)
          aa = part(ia)

          i1 = abs(ii)
          !if (ii > 0) i1 = ii + noa + nva

          a1 = abs(aa) + noa
          !if (aa > 0) a1 = aa + noa + nva + noa

          ! Add diagonal contributions
          idx = (i1-1)*nrorb + i1
          density_complex(idx) = density_complex(idx) + dconjg(psi_ia) * psi_ia
          idx = (a1-1)*nrorb + a1
          density_complex(idx) = density_complex(idx) + dconjg(psi_ia) * psi_ia

          ! Add off-diagonal contributions
          do jb = 1, nstates
              psi_jb = psi_det(jb)

              jj = hole(jb)
              bb = part(jb)

              j1 = abs(jj)
              !if (jj > 0) j1 = jj + noa + nva

              b1 = abs(bb) + noa
              !if (bb > 0) b1 = bb + noa + nva + noa

              ! Off-diagonal element between states (ia, jb)
              idx = (i1-1)*nrorb + j1
              density_complex(idx) = density_complex(idx) + dconjg(psi_ia) * psi_jb
              idx = (a1-1)*nrorb + b1
              density_complex(idx) = density_complex(idx) + dconjg(psi_ia) * psi_jb
          end do
      end do
  end subroutine make_complex_density



  subroutine make_transition_rate(nrorb, ham_mo, dipxmoa, dipymoa, dipzmoa, efieldx, efieldy, efieldz, vabsmoa, density_MO, transition_rates)

    ! Arguments
    implicit none
    integer(8), intent(in) :: nrorb
    real(8), intent(in) :: ham_mo(nrorb*nrorb)
    real(8), intent(in) :: dipxmoa(nrorb*nrorb)
    real(8), intent(in) :: dipymoa(nrorb*nrorb)
    real(8), intent(in) :: dipzmoa(nrorb*nrorb)
    real(8), intent(in) :: efieldx, efieldy, efieldz
    real(8), intent(in) :: vabsmoa(nrorb*nrorb)
    complex(8), intent(in) :: density_MO(nrorb*nrorb)
    real(8), intent(out) :: transition_rates(nrorb*nrorb)

    ! Local variables
    integer(8) :: i, j, idx
    real(8) :: field_contribution, absorbing_contribution
    real(8) :: hamT_idx
    complex(8) :: density_element
    real(8) :: hbar, rate_contribution

    ! Constants
    hbar = 1.0d0  ! Planck's constant (arbitrary units)

    ! Initialize transition rates to zero
    transition_rates = 0.d0

    ! Loop over all alpha orbitals (restricted wavefunction: nrorb x nrorb block)
    do i = 1, nrorb
        do j = 1, nrorb
            idx = (i-1)*nrorb+j

            ! Add field contributions: E(t) * mu
            field_contribution = efieldx * dipxmoa(idx) + efieldy * dipymoa(idx) + efieldz * dipzmoa(idx)

            ! Total Hamiltonian element including time-dependent contributions
            hamT_idx = ham_mo(idx) + field_contribution + vabsmoa(idx)

            ! Compute transition rate contribution
            !: The 2.0 factor is from the Liouville-von Neumann expansion, not the restricted density matrix.
            transition_rates(idx) = 2.d0*aimag(hamT_idx * conjg(density_MO(idx)))

        end do
    end do

  end subroutine make_transition_rate












end module analysis
