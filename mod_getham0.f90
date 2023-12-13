  module get_ham0
  
  use global_variables
  use util
  
  implicit none

  integer(8) :: double_startAA, double_startBB, double_startAB, double_startBA


  !: contains subroutine get_cis
  !: contains subroutine get_cis_index
  !: contains subroutine get_ip_cisd
  !: contains subroutine get_sip_cis
  !: contains subroutine get_sip_soc
  !: contains subroutine get_dip_cis
  !: contains subroutine get_dip_soc
  !: contains subroutine get_ip_index
  !: contains subroutine diagonalize
  !: contains subroutine writeme_ham0( myroutine, option )
  !: contains subrotuine errors_getham0( myroutine, option, mystring )
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_RCIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cis
    
    !: Form restricted or unrestricted  CIS Hamiltonian                                       
    !: ia = (i-1)*nva + a ; ib = (i-1)*nva + b
    !: jb = (j-1)*nvb + b ; ja = (j-1)*nvb + a     
    
    use read_integrals
    implicit none
    
    integer(8) :: istate, i, j, k, l, a, b
    integer(8) :: ia2, jb2, ia, ib, ja, jb
    real(8)    :: hmat(nstates,nstates)
    

    call write_header( 'get_cis','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'cis','form' ) 
    call get_cis_index
    
    !: temporarily store ham.  will transfer ham(rank2) to cis_vec(rank1)
    hmat   = 0.d0    
    
    if ( .not.unrestricted ) then    
       
       restricted : do ia = 2, nstates
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          ia2 = ia - 1 !: dijabAB index
          do jb = 2, nstates
             j = -hole_index(jb,1)
             b = -part_index(jb,1)   
             jb2 = jb - 1        !: dijabAB index
             ib  = (i-1)*nva + b !: diajbAA index
             ja  = (j-1)*nva + a !: diajbAA index
             hmat(jb,ia) = dijabAB(jb2,ia2) - diajbAA(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(noa+a) - orben(i)
       end do restricted
       
    else 
       
       aa : do ia = 2 , 1+noanva
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          do jb = 2 , 1+noanva
             j = -hole_index(jb,1)
             b = -part_index(jb,1)             
             ib = (i-1)*nva + b !: diajbAA index
             ja = (j-1)*nva + a !: diajbAA index
             hmat(jb,ia) = -diajbAA(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(noa+a) - orben(i)
       end do aa

       !: beta, beta block 
       bb : do ia = 2+noanva, nstates
          i = hole_index(ia,1)
          a = part_index(ia,1)
          do jb = 2+noanva, nstates
             j  = hole_index(jb,1)
             b  = part_index(jb,1)       
             ib = nvb*(i-1) + b !: diajbBB index
             ja = nvb*(j-1) + a !: diabBB index
             hmat(jb,ia) = -diajbBB(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(a+nob+nrorb) - orben(i+nrorb)
       end do bb
       
       !: alpha,beta block ;; beta,alpha block
       ab : do ia = 2, 1+noanva
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          ia2 = (i-1)*nva + a
          do jb = 2+noanva, nstates
             j  = hole_index(jb,1)
             b  = part_index(jb,1)                          
             jb2 = (j-1)*nvb + b
             hmat(jb,ia) = dijabAB(jb2,ia2)
             hmat(ia,jb) = hmat(jb,ia)
          end do
       end do ab
       
    end if
    
    !: transfer to cis_vec
    istate=0
    do i=1, nstates
       do j=1, nstates
          istate = istate + 1 
          cis_vec( istate ) = hmat(j,i)
       end do
    end do
 
    if ( Qwrite_ham0 ) then
      open(unit=100,file='HAM.OUT')
      do i=1, nstates
        do j=1, nstates
          write(100,"(i5,i5,f15.10,f15.10)") i,j, hmat(j,i)
        end do
      end do
      close(100)    
    end if

    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start
    
    
    call write_header( 'get_cis','get_ham0','leave' )

    
  end subroutine get_cis
  !:-------------------------!
  !: subroutine get_soc_cis  !
  !:-------------------------!
  subroutine get_soc_cis


    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS
    !:     Form CIS Hamiltonian with spin-orbit coupling
    !:     (RHF or ROHF not not UHF)
    !:
    !:             BKT NO  INTEGRAL  SPIN-CASE    SYM           SEQ
    !:      --------------------------------------------------------
    !:      Dijab     2    (IJ//AB)    ABAB       ALL           IAJB
    !:      DiajbA    5    (IA//JB)    AAAA       ALL           IAJB
    !:      DiajbAB   7    (IA//JB)    ABAB   I.LE.J,A.LE.B     IJAB
    !:      DiajbBA   8    (IA//JB)    BABA   I.LE.J,A.LE.B     IJAB
    !:      DiajbB   10    (IA//JB)    BBBB       ALL           IAJB
    !:
    !:     CIS state order and two electron matrix elements
    !:
    !:       ia \ jb  |    0    alpha,alpha; beta,beta; alpha,beta; beta,alpha
    !:     ----------------------------------------------------------
    !:          0     |    0         0          0           0          0
    !:     alpha,alpha|    0      -DiajbA      Dijab        0          0
    !:     beta,beta  |    0       Dijab      -DiajbB       0          0
    !:     alpha,beta |    0         0          0        DiajbAB       0
    !:     beta,alpha |    0         0          0           0        DiajbBA
    !:
    !:
    !:     One electron spin-orbit matrix elements
    !:       Vx = (mu/SOCx/nu)/i, Vy, Vz in AO basis from Gaussian
    !:       VZA = (p/SOCz/q) in MO basis for alpha,alpha
    !:       VZB = (p/SOCz/q) in MO basis for beta,beta
    !:       VPM = (p/SOCx/q)+i(p/SOCy/q) for alpha,beta
    !:       VM(p,q) = VPM(p,q), VP(p,q) = VPM(q,p)*
    !:
    !: 
    !:       ia \ jb  |  0    alpha,alpha; beta,beta; alpha,beta; beta,alpha
    !:     -----------------------------------------------------------------
    !:          0     |  0       VZAjb      VZBjb       VPMjb      -VPMbj*
    !:     alpha,alpha|VZAia* VZAij,VZAab     0         VPMab       VPMij
    !:     beta,beta  |VZBia*      0     VZBij,VZBab   -VPMji*     -VPMba*
    !:     alpha,beta |VPMia*    VPMba*    -VPMij     VZAij,VZBab     0
    !:     beta,alpha |-VPMai    VPMji*    -VPMab          0     VZBij,VZAab
    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS
    

    use read_integrals
    implicit none

    integer(8) :: istate, ii, jj, aa, bb, i, j, a, b
    integer(8) :: ia, jb, aorb, iorb
    complex(8) :: cdum, hmat(nstates,nstates)
  
    
    call write_header( 'get_soc_cis','get_ham0','enter' )

    call cpu_time(start)
    call get_soc_ao2mo
    call get_cis_index
    

    hmat = dcmplx( 0.d0,0.d0 )
    
    
    ia = 1 
    do i=1, noa
       do a=1, nva
          ia = ia + 1 
          iorb = i
          aorb = noa + a          
!:          if(-i.ne.hole_index(ia,1).or.-a.ne.part_index(ia,1)) write(42,*) "ia=",ia,i,a
          
          !: < 0 |H| i->a >
          hmat(1,ia) = socmoAA(i,noa+a) 
          hmat(ia,1) = dconjg( hmat(1,ia) ) 
          
          !: < i->a |H| j->b >
          do jb=2, ia
             j = abs( hole_index(jb,1) )
             b = abs( part_index(jb,1) )
             !: - <ib||ja>
             cdum = dcmplx( -get_diajbAA(i,b,j,a), 0.d0 )
             if ( i.eq.j ) cdum = cdum + socmoAA(noa+a,noa+b)
             if ( a.eq.b ) cdum = cdum - socmoAA(j,i)
             hmat(ia,jb) = cdum
             hmat(jb,ia) = dconjg(cdum)
          end do
          hmat(ia,ia) = hmat(ia,ia) + dcmplx( orben(aorb)-orben(iorb), 0.d0 )
          
          !: < i->a |H| J->B >
          JB = 1 + noanva 
          do J=1, NOB
             do B=1, NVB
                JB = JB + 1 
!:          if(J.ne.hole_index(JB,1).or.B.ne.part_index(JB,1)) &
!:            write(42,"(A,5I5)") "(1) jb=",JB,J,B,hole_index(JB,1),part_index(JB,1)
                !: <iJ||aB>
                cdum = dcmplx( get_dijabAB(i,J,a,B), 0.d0 )
                hmat(ia,JB) = cdum
                hmat(JB,ia) = dconjg( cdum ) 
             end do
          end do
          
          !: < i->a |H| i->B >
          j = i
          jB = 1 + noanva + NOBNVB + (j-1)*NVB
          do B=1, NVB
             jB = jB + 1
!:          if(-J.ne.hole_index(JB,1).or.B.ne.part_index(JB,1))&
!:            write(42,"(A,5I5)") "(2) jb=",JB,J,B,hole_index(JB,1),part_index(JB,1)
             hmat(ia,jB) = socmoAB(a+noa,B+NOB)
             hmat(jB,ia) = dconjg( hmat(ia,jB) ) !: < i->B |H| i->a >
          end do
          
          !: < i->a |H| J->a >
          b = a 
          do J=1, NOB
             Jb = 1 + noanva + NOBNVB + noaNVB + (J-1)*nva + b
!:          if(J.ne.hole_index(JB,1).or.-B.ne.part_index(JB,1)) &
!:            write(42,"(A,5I5)") "(3) jb=",JB,J,B,hole_index(JB,1),part_index(JB,1)
             hmat(Jb,ia) = - socmoAB(i,J) !: < J->a |H| i->a >
             hmat(ia,Jb) = dconjg( hmat(Jb,ia) )
          end do
          
       end do
    end do

    
    IA = 1 + noanva
    do I=1, NOB
       do A=1, NVB
          IA = IA + 1 
!:          if(i.ne.hole_index(ia,1).or.a.ne.part_index(ia,1)) write(42,*) "ia=",ia,i,a
          AORB = nrorb + NOB + A 
          IORB = nrorb + I
          
          !: < 0 |H| I->A >
          hmat(1,IA) = socmoBB(I,NOB+A)
          hmat(IA,1) = dconjg( hmat(1,IA) )
          
          !: < I->A |H| J->B >
          do JB=(1+NOANVA+1), IA
             J = hole_index(JB,1)
             B = part_index(JB,1)
             !: - <ib||ja>
             cdum = dcmplx( -get_diajbBB(I,B,J,A), 0.d0 )
             if ( I.eq.J ) cdum = cdum + socmoBB(A+NOB,B+NOB)
             if ( A.eq.B ) cdum = cdum - socmoBB(J,I)
             hmat(IA,JB) = cdum
             hmat(JB,IA) = dconjg( cdum )
          end do
          hmat(IA,IA) = hmat(IA,IA) + dcmplx( orben(AORB)-orben(IORB), 0.d0 )
          
          !: < I->A |H| j->A >
          B = A
          do j=1, noa
             JB = 1 + noanva + NOBNVB + (j-1)*NVB + B
          if(-J.ne.hole_index(JB,1).or.B.ne.part_index(JB,1)) &
!:            write(42,"(A,5I5)") "(4) jb=",JB,J,B,hole_index(JB,1),part_index(JB,1)
!:             hmat(IA,jB) = - socmoAB(j,I)
             hmat(jB,IA) = dconjg( hmat(IA,jB) )
          end do
          
          !: < I->A |H| I->b >
          J = I
          Jb = 1 + noanva + NOBNVB + noaNVB + (J-1)*nva
          do b=1, nva
             Jb = Jb + 1
!:          if(J.ne.hole_index(JB,1).or.-B.ne.part_index(JB,1)) &
!:            write(42,"(A,5I5)") "(5) jb=",JB,J,B,hole_index(JB,1),part_index(JB,1)
             hmat(Jb,IA) = socmoAB(noa+b,NOB+A) !: < I->b |H| I->A >
             hmat(IA,Jb) = dconjg( hmat(Jb,IA) )
          end do
          
       end do
    end do

    
    iA = 1 + noanva + NOBNVB
    do i=1, noa
       do A=1, NVB
          iA = iA + 1 
          if(-i.ne.hole_index(ia,1).or.a.ne.part_index(ia,1)) write(42,*) "ia=",ia,i,a
          iorb = i
          AORB = nrorb + NOB + A 

          !: < i->A |H| 0 >
          hmat(1,iA) = socmoAB(i,NOB+A) !: < 0 |H| i->A > 
          hmat(iA,1) = dconjg( hmat(1,iA) )
          
          !: < i->A |H| j->B >
          do jB = ( 1+noanva+NOBNVB+1 ), iA
             j = abs( hole_index(jB,1) )
             B = part_index(jB,1)
             !: - <iB||jA>
             cdum = dcmplx( -get_diajbAB(i,B,j,A), 0.d0 )
             if ( i.eq.j ) cdum = cdum + socmoBB(NOB+A,NOB+B)
             if ( A.eq.B ) cdum = cdum - socmoAA(j,i)
             hmat(iA,jB) = cdum
             hmat(jB,iA) = dconjg( cdum )
          end do
          hmat(iA,iA) = hmat(iA,iA) + dcmplx( orben(AORB)-orben(iorb), 0.d0 )
       end do
    end do

    
    Ia = 1 + noanva + NOBNVB + noaNVB
    do I=1, NOB
       do a=1, nva
          Ia = Ia + 1 
          if(i.ne.hole_index(ia,1).or.-a.ne.part_index(ia,1)) write(42,*) "ia=",ia,i,a
          IORB = nrorb + I
          aorb = noa + a 
          
          !: < I->a |H| 0 >
          hmat(Ia,1) = socmoAB(noa+a,I)
          hmat(1,Ia) = dconjg( hmat(Ia,1) ) !: < 0 |H| I->a >

          !: < I->a |H| J->b >
          do Jb=(1+noanva+NOBNVB+noaNVB+1), ia 
             J = hole_index(Jb,1)
             b = abs(part_index(Jb,1))
             !: -<Ib||Ja>
             cdum = dcmplx( -get_diajbBA(I,b,J,a), 0.d0 )
             if ( I.eq.J ) cdum = cdum + socmoAA(noa+a,noa+b)
             if ( a.eq.b ) cdum = cdum - socmoBB(J,I)
             hmat(Ia,Jb) = cdum
             hmat(Jb,Ia) = dconjg( cdum )
          end do
          hmat(Ia,Ia) = hmat(Ia,Ia) + dcmplx( orben(aorb)-orben(IORB), 0.d0 )
       end do
    end do

    
    ia = 0
    Zcis_vec = dcmplx(0.d0,0.d0)
    do i=1, nstates
       do j=1, nstates
          ia = ia + 1 
          Zcis_vec(ia) = hmat(j,i)
       end do
    end do
    
    if ( Qwrite_ham0 ) then
      open(unit=100,file='HAM.OUT')
      do i=1, nstates
        do j=1, nstates
          write(100,"(i5,i5,f15.10,f15.10)") i,j, hmat(j,i)
        end do
      end do
      close(100)
    end if

    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start


    call write_header( 'get_soc_cis','get_ham0','leave' )


  end subroutine get_soc_cis
  !:---------------------------!
  !: SUBROUTINE GET_SOC_AO2MO
  !:---------------------------!
  subroutine get_soc_ao2mo

    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS
    !: Form complex spin-orbit matrices VPM,VZA,VZB in MO basis
    !: real VX, VY and VZ integrals in AO basis

    !:
    !: AO basis, lower triangle, real
    !: VX=(mu/SOCx/nu)/i, VY=(mu/SOCy/nu)/i, VZ=(mu/SOCz/nu)/i
    !:
    !: MO basis, full matrix, complex (Z, raising and lowering)
    !:    VZA = (p/SOCz/q) for alpha,alpha
    !:    VZB = (p/SOCz/q) for beta,beta
    !:    VPM = (p/SOCx/q)+i(p/SOCy/q) for alpha,beta
    !:    VM(p,q) = VPM(p,q), VP(p,q) = VPM(q,p)*
    !:
    !: conversion factor for spin-orbit integrals
    !: (see Salem, Angnew Chem Internat - Vol 11 (1972),No 2,pp92-111)
    !: constant=((h/2pi)**2*e**2)/4*m***2*c**2 = 2.9217cm**-1
    !: = 0.00001331224 au
    !: HBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBSHBS

    use util

    implicit none
    real(8), parameter :: constant = 0.00001331224d0
    
    integer(8) :: i, j, ij
    complex(8), allocatable :: scratch1(:,:)


    write(iout,'(A)') ' in subroutine get_soc_ao2mo in MODULE get_ham0 '

    !: convert to au units
    socxao(:,:) = socfac * constant * socxao(:,:)
    socyao(:,:) = socfac * constant * socyao(:,:)
    soczao(:,:) = socfacz *  constant * soczao(:,:)
!:    soczao(:,:) = socfacz * socfac * constant * soczao(:,:)

    
    allocate( scratch1(nbasis,nbasis) )
    scratch1 = 0.d0
    
    !: socmoAA
    scratch1(:,:) = - soczao(:,:)
    call ao2mo_complex(nbasis, nrorb, scratch1, socmoAA, cmo_a, cmo_a )

    open( unit=100,file='VZA.OUT' )
    do i=1, nrorb
       do j=1, nrorb
          write(100,"(f15.10,f15.10)") socmoAA(j,i)
       end do
    end do
    close(100)

    
    !: socmoBB
    scratch1(:,:) = soczao(:,:)
    call ao2mo_complex(nbasis, nrorb, scratch1, socmoBB, cmo_b, cmo_b )

    open( unit=100,file='VZB.OUT' )
    do i=1, nrorb
       do j=1, nrorb
          write(100,"(f15.10,f15.10)") socmoBB(j,i)
       end do
    end do
    close(100)
    
    
    !: socmoAB
    do i=1, nbasis
       do j=1, nbasis
          scratch1(j,i) =  - ( socxao(j,i) - eye * socyao(j,i) )
       end do
    end do
    call ao2mo_complex(nbasis, nrorb, scratch1, socmoAB, cmo_a, cmo_b )
    
    
    open( unit=100,file='VPM.OUT' )
    do i=1, nrorb
       do j=1, nrorb
          write(100,"(f15.10,f15.10)") socmoAB(j,i)
       end do
    end do
    close(100)


    deallocate( scratch1 )

    write(iout,'(A)') ' finished assigning socmoAA, socmoBB, socmoAB ' 
    write(iout,'(A)') ' leaving subroutine get_soc_ao2mo'
    

  end subroutine get_soc_ao2mo
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_RCIS_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cis_index


    !: assign hole and particle excitation indices for each state
    implicit none
    
    integer(8) :: i, j, a, b, istate
    
    
    call writeme_ham0( 'cis_index', 'imhere' )

    !: HF reference state
    istate = 1
    hole_index(istate,1) = 0
    part_index(istate,1) = 0
    
    do i=1, noa
       do a=1, nva
          istate = istate + 1 
          hole_index(istate,1) = -i
          part_index(istate,1) = -a ! noa + a 
       end do
    end do    

    if ( unrestricted ) then
       do I=1, NOB
          do A=1, NVB
             istate = istate + 1 
             hole_index(istate,1) = I ! noa+nva+i
             part_index(istate,1) = A ! noa+nva+nob+a
          end do
       end do
    end if

    
    !: with spin orbit coupling
    if (QsocA2B .and. trim(jobtype) .eq. flag_soc ) then
       !: alpha -> beta
       do i=1, noa
          do A=1, NVB
             istate = istate + 1 
             hole_index(istate,1) = -i
             part_index(istate,1) = A
          end do
       end do
       !: beta -> alpha
       do I=1, NOB
          do a=1, nva
             istate = istate + 1 
             hole_index(istate,1) = I
             part_index(istate,1) = -a
          end do
       end do
    end if

    if ( istate.ne.nstates) write(iout,'(A)')  ' ERROR ERROR ERROR IN INDEX ASSIGNMENT IN CIS' 
    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    close(100)    
    
    call writeme_ham0( 'cis_index','imdone' )
    

  end subroutine get_cis_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_IP_CISD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_ip_cisd

    !: get IP-CISD Hamiltonian

    use read_integrals
    implicit none

    integer(8) :: istate, xorb, iorb, aorb
    integer(8) :: xx, ii, aa, yy, jj, bb, x, y, i, a, j, b
    integer(8) :: ia, jb, ia2, jb2, yx, xb, yj, ji, ix, xi, jy, ab, ya, ib, ja, yb, xa
    real(8) :: sign, storeme
    
    
    call write_header( 'get_ip_cisd','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'ip', 'form' ) 
    call get_ip_index 
    
    ia : do ia=1, nstates
       
       xx = hole_index(ia,1) ; x = abs(xx)
       ii = hole_index(ia,2) ; i = abs(ii)
       aa = part_index(ia,1) ; a = abs(aa)
       
       !: orbital indices
       if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii
       if ( aa.lt.0 ) aorb = -aa + noa ; if ( aa.gt.0 ) aorb = nrorb + nob + aa
       
       jb : do jb=ia, nstates
          
          yy = hole_index(jb,1) ; y = abs(yy)
          jj = hole_index(jb,2) ; j = abs(jj)
          bb = part_index(jb,1) ; b = abs(bb)
          
          storeme = 0.d0
          
          SS : if( ii.eq.0 .and. jj.eq.0 ) then
             if( xx.eq.yy ) storeme = -orben(xorb)
             go to 78
          end if SS
          
          
          SD : if ( ii*jj .eq. 0 ) then

             if ( ii.eq.0 ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then                   
                   !: <x|H|jb_y> = -<yj||xb> = -<ij||ka> : DijkaAA(ka,i<j) = DijkaAA(xb,y<j)
                   if ( jj.lt.0 ) storeme = - get_dijkaAA(y,j,x,b)
                   !: <x|H|JB_y> = -<yJ||xB> = -<iJ||kA> : DijkaAB(i<=k,JA) = DijkaAB(y<=x,JB)
                   if ( JJ.gt.0 ) storeme = - get_dijkaAB(y,J,x,B)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <X|H|jb_Y> = -<Yj||Xb> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(jb2,Y<=X)
                   if ( jj.lt.0 ) storeme = - get_dijkaBA(Y,j,X,b)
                   !: <X|H|JB_Y> = -<YJ||XB> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(XB,Y<J)
                   if ( JJ.gt.0 ) storeme = - get_dijkaBB(Y,J,X,B)
                end if
             else if ( jj.eq.0 ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <y|H|ia_x> = -<xi||ya> = -<ij||ka> :  DijkaAA(ka,i<j) = DijkaBB(ya,x<i)                   
                   if ( ii.lt.0 ) storeme = - get_dijkaAA(x,i,y,a)
                   !: <y|H|IA_x> = -<xI||yA> = -<iJ||kA> : DijkaAB(JA,i<=k) = DijkaAB(IA,x<=y)
                   if ( II.gt.0 ) storeme = - get_dijkaAB(x,I,y,A)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <Y|H|ia_X> = -<Xi||Ya> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(ia, X<=Y)
                   if ( ii.lt.0 ) storeme = - get_dijkaBA(X,i,Y,a)
                   !: <Y|H|IA_X> = -<XI||YA> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(YA,X<I)
                   if ( II.gt.0 ) storeme = - get_dijkaBB(X,I,Y,A)
                end if
             end if

             go to 78
             
          end if SD
          
          
          DD : if( ii.ne.0 .and. jj.ne.0 ) then
             
             !: < ix->a |H| jy->a >
             kdelta_ab : if ( aa.eq.bb ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <jy||ix> = <ij||kl> : DijklAA(k<l,i<j)=DijklAA(i<x,j<y)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAA(j,y,i,x)
                   !: <Jy||Ix> = <yJ||xI> = <iJ||kL> : DijklAB(J<=L,i<=k) = DijklAB(J<=I,y<=x)                   
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklAB(y,J,x,I)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <jY||iX> = <iJ||kL> : DijklAB(J<=L,i<=k)=DijklAB(Y<=X,j<=i)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAB(j,Y,i,X)
                   !: <JY||IX> = <IJ||KL> : DijklBB(k<l,i<j)=DijklBB(I<X,J<Y)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklBB(J,Y,I,X)
                end if
             end if kdelta_ab
             
             !: < ix->a |H| iy->b >
             kdelta_ij : if ( ii.eq.jj ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: -<ya||xb> = <ia||jb> : DiajbAA(jb,ia) = DiajbAA(xb,ya)
                   if ( aa.lt.0 .and. bb.lt.0 ) storeme = storeme - get_diajbAA(y,a,x,b)
                   !: -<yA||xB> = <iA||jB> : DiajbAB(A<=B,i<=j) = DiajbAB(A<=B,y<=x)
                   if ( AA.gt.0 .and. BB.gt.0 ) storeme = storeme - get_diajbAB(y,A,x,B)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: -<Ya||Xb> = <Ia||Jb> : DiajbBA(a<=b,I<=J)=DiajbBA(a<=b,Y<=X)
                   if ( aa.lt.0 .and. bb.lt.0 ) storeme = storeme - get_diajbBA(Y,a,X,b)
                   !: -<YA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,YA)                   
                   if ( AA.gt.0 .and. BB.gt.0 ) storeme = storeme - get_diajbBB(Y,A,X,B)
                end if
             end if kdelta_ij
             
             
             !: < ix->a |H| jx->b >
             kdelta_xy : if ( XX.eq.YY ) then
                !: -<ja||ib> = -<ia||jb> : DiajbAA(jb,ia) = DiajbAA(ib,ja)
                if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme - get_diajbAA(j,a,i,b)
                !: -<JA||IB> = -<ia||jb> : DiajbBB(jb,ia) = DiajbBB(IB,JA)
                if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme - get_diajbBB(J,A,I,B)
                !: <iJ||aB>  : DijabAB(JB,ia)
                if ( ii.lt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijabAB(i,J,a,B)
                !: <Ij||Ab> = <jI||bA> : DijabAB(IA,jb)
                if ( II.gt.0 .and. jj.lt.0 ) storeme = storeme + get_dijabAB(j,I,b,A)
             end if kdelta_xy


             !: < ix->a |H| xy->b >
             kdelta_jx : if ( JJ.eq.XX ) then
                aa_or_bb : if ( jj.lt.0 ) then
                   !: < ix->a |H| xy->b > :  <ya||ib> = <ia||jb> : DiajbAA(jb,ia)=DiajbAA(ib,ya)
                   if ( ii.lt.0 .and. yy.lt.0 ) storeme = storeme + get_diajbAA(y,a,i,b)
                   !: < Ix->A |H| xy->b > : - <Iy||Ab> = -<yI||bA>
                   if ( II.gt.0 .and. yy.lt.0 ) storeme = storeme - get_dijabAB(y,I,b,A)
                else if ( JJ.gt.0 ) then
                   !: < IX->A |H| XY->B > :  <YA||IB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(IB,YA)
                   if ( II.gt.0 .and. YY.gt.0 ) storeme = storeme + get_diajbBB(Y,A,I,B)
                   !: < iX->a |H| XY->B > : -<iY||aB>
                   if ( ii.lt.0 .and. YY.gt.0 ) storeme = storeme - get_dijabAB(i,Y,a,B)                   
                end if aa_or_bb
             end if kdelta_jx
             

             !: < yx->a |H| jy->b >
             kdelta_iy : if ( II.eq.YY ) then
                aa_or_bb2:  if ( ii.lt.0 ) then
                   !: < yx->a |H| jy->b > :  <ja||xb> = <ia||jb> : DiajbAA(jb,ia)=DiajbAA(xb,ja)
                   if ( jj.lt.0 .and. xx.lt.0 ) storeme = storeme + get_diajbAA(j,a,x,b)
                   !: < yx->a |H| Jy->B > :  - <xJ||aB> 
                   if ( JJ.gt.0 .and. xx.lt.0 ) storeme = storeme - get_dijabAB(x,J,a,B)
                else if ( II.gt.0 ) then
                   !: < YX->A |H| JY->B > :  <JA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,JA)
                   if ( JJ.gt.0 .and. XX.gt.0 ) storeme = storeme + get_diajbBB(J,A,X,B)
                   !: < YX->A |H| jY->b > :  - <jX||bA>
                   if ( jj.lt.0 .and. XX.gt.0 ) storeme = storeme - get_dijabAB(j,X,b,A)                   
                end if aa_or_bb2
             end if kdelta_iy                         
             
          end if DD

          if ( jb.eq.ia ) storeme = storeme - orben(iorb) - orben(xorb) + orben(aorb)
          
78        continue

          istate = (ia-1)*nstates + jb
          cis_vec(istate) = storeme

          istate = (jb-1)*nstates + ia
          cis_vec(istate) = storeme
          
       end do jb
    end do ia
   
    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start    

    call write_header( 'get_ip_cisd','get_ham0','leave' )
    

  end subroutine get_ip_cisd
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_SIP_CIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_sip_cis

    !: assign hole and particle excitation indices for singly ionized states
    implicit none

    integer(8) :: lscratch, info, lwork
    integer(8) :: i, j, ii, ia, x, xx, xorb, istate, xstart, xstart1
    real(8), allocatable :: work(:)
    character(20) :: cinfo

    call writeme_ham0( 'sip_cis','imhere' )
   
    part_ip_index(:,:) = 0
    hole_ip_index(:,:) = 0
    state_ip_index = 0
    istate = 0
    xstart = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
       singlesA : do x=xstart, noa
          istate = istate + 1 
          hole_ip_index(istate,1) = -x
          state_ip_index(x,1) = istate
          state_ip_index(1,x) = istate
       end do singlesA
    end if

    singlesB : do x=xstart, nob
       istate = istate + 1 
       hole_ip_index(istate,1) = x
       state_ip_index(x+noa,1) = istate
       state_ip_index(1,x+noa) = istate
    end do singlesB
 
!:    write(iout,*) " sip_states",ip_states,istate
!:    
!:      do i=1, ip_states 
!:        write(iout,123) i,hole_ip_index(i,1),hole_ip_index(i,2),part_ip_index(i,1) 
!:      end do 
!:      do i=1, noa+nob
!:        write(iout,'(4i5)') i, state_ip_index(i,1)
!:      end do 

    ip_vec = 0.d0

    do ia=1, ip_states
       ii = ia + ip_states*(ia-1)
       xx = hole_ip_index(ia,1) 
       xorb = -xx
       if( xx.gt.0 ) xorb = nrorb + xx
       ip_vec(ii) = -orben(xorb)
!:       write(iout,"(4i5,f12.6)") ia,ii,xx,xorb,orben(xorb)
    end do
!:    write(iout,"(8f12.6)") ip_vec

    info = 10
    lwork = ip_states**2
    allocate( work(lwork) )
    call dsyev('vectors','upper',ip_states,ip_vec,ip_states,ip_eig,work,lwork,info)
    call fix_phase(ip_states,ip_vec)
    deallocate( work )

    write(iout,'(A)') ' ENERGIES for SINGLY IONIZED STATES'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
    end do

    write(iout,'(A)') ' EIGENVECTORS for SINGLY IONIZED STATES (coeff and |coeff|**2)'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
      do ia=1, ip_states
        ii = ia+ip_states*(i-1)
        if(abs(ip_vec(ii)).gt.1.d-1) write(iout,60) hole_ip_index(ia,1), ip_vec(ii)
      end do
    end do

    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    do i=1, ip_states
       write(100,"(i7,i7,i7)") hole_ip_index(i,1), hole_ip_index(i,2), part_ip_index(i,1)
    end do
    do i=1, noa+nob
       write(100,"(400i5)") (state_ip_index(i,j),j = 1,noa+nob)
    end do
    close(100)    
    
50  format( i5,1x,f15.10 )
60  format( 7x,i3,2x,4f17.10 )

    call writeme_ham0( 'sip_cis','imdone' )
    
  end subroutine get_sip_cis
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_SIP_SOC
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_sip_soc

    !: assign hole and particle excitation indices for singly ionized states
    implicit none

    integer(8) :: lscratch, info, lwork
    integer(8) :: i, j, ii, ia, x, xx, xorb, istate, xstart, xstart1
    real(8), allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    character(20) :: cinfo

    call writeme_ham0( 'sip_soc','imhere' )
   
    part_ip_index(:,:) = 0
    hole_ip_index(:,:) = 0
    state_ip_index = 0
    istate = 0
    xstart = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
       singlesA : do x=xstart, noa
          istate = istate + 1 
          hole_ip_index(istate,1) = -x
          state_ip_index(x,1) = istate
          state_ip_index(1,x) = istate
       end do singlesA
    end if

    singlesB : do x=xstart, nob
       istate = istate + 1 
       hole_ip_index(istate,1) = x
       state_ip_index(x+noa,1) = istate
       state_ip_index(1,x+noa) = istate
    end do singlesB
!: 
!:    write(iout,*) " sip_states",ip_states,istate
!:    
!:      do i=1, ip_states 
!:        write(iout,123) i,hole_ip_index(i,1),hole_ip_index(i,2),part_ip_index(i,1) 
!:      end do 
!:      do i=1, noa+nob
!:        write(iout,123) i, state_ip_index(i,1)
!:      end do 
!: 123  format(4i5)

    Zip_vec = dcmplx(0.d0,0.d0)
    call Zform1det_socip(noa,nva,ip_states,Zip_vec,socmoAA, &
      hole_ip_index,part_ip_index,nob,nvb,socmoBB,socmoAB)

    do ia=1, ip_states
       ii = ia + ip_states*(ia-1)
       xx = hole_ip_index(ia,1) 
       xorb = -xx
       if( xx.gt.0 ) xorb = nrorb + xx
       Zip_vec(ii) = dcmplx(-orben(xorb),0.d0)
    end do

    info = 10
    lwork = ip_states**2
    allocate( work(lwork) )
    allocate( rwork(3*ip_states-2) ) ;  rwork = 0.d0
    call zheev('vectors','upper',ip_states,Zip_vec,ip_states,ip_eig,work,lwork,rwork,info)
    call Zfix_phase(ip_states,Zip_vec)
    deallocate( rwork )
    deallocate( work )

    write(iout,'(A)') ' ENERGIES for SINGLY IONIZED STATES'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
    end do

    write(iout,'(A)') ' EIGENVECTORS for SINGLY IONIZED STATES (coeff and |coeff|**2)'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
      do ia=1, ip_states
        ii = ia+ip_states*(i-1)
        if(abs(Zip_vec(ii)).gt.1.d-1) write(iout,60) hole_ip_index(ia,1), &
           Zip_vec(ii),dble(dconjg(Zip_vec(ii))*Zip_vec(ii))
      end do
    end do

    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    do i=1, ip_states
       write(100,"(i7,i7,i7)") hole_ip_index(i,1), hole_ip_index(i,2), part_ip_index(i,1)
    end do
    do i=1, noa+nob
       write(100,"(400i5)") (state_ip_index(i,j),j = 1,noa+nob)
    end do
    close(100)    
    
50  format( i5,1x,f15.10 )
60  format( 7x,i3,2x,4f17.10 )

    call writeme_ham0( 'sip_soc','imdone' )
    
  end subroutine get_sip_soc
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_DIP_CIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_dip_cis

    !: assign hole and particle excitation indices for doubly ionized states

    use read_integrals
    use util

    implicit none

    integer(8) :: lscratch, info, lwork
    integer(8) :: i, j, x, y, istate, xstart, xstart1
    integer(8) :: ia, jb, xx, ii, yy, jj, xorb, iorb
    real(8) :: sign, storeme
    real(8), allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    character(20) :: cinfo

    call writeme_ham0( 'dip_cis','imhere' )
    
    part_ip_index(:,:) = 0
    hole_ip_index(:,:) = 0
    state_ip_index = 0
    istate = 0
    xstart  = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
 
!:  alpha,alpha 
        do x=xstart1, noa
          if(xstart.le.x-1) then
            do i=xstart, x-1
              istate = istate + 1 
              hole_ip_index(istate,1) = -x
              hole_ip_index(istate,2) = -i
              state_ip_index(i,x) = istate
              state_ip_index(x,i) = istate
            end do
          end if
       end do 
    end if
 
!:  alpha,beta 
        do x=xstart, noa
          do i=xstart, nob
            istate = istate + 1 
            hole_ip_index(istate,1) = -x
            hole_ip_index(istate,2) = i
            state_ip_index(i+noa,x) = istate
            state_ip_index(x,i+noa) = istate
         end do
       end do 
    
!:  beta,beta
      do x=xstart1, nob
        if(xstart.le.x-1) then
          do i=xstart, x-1
            istate = istate + 1 
            hole_ip_index(istate,1) = x
            hole_ip_index(istate,2) = i
            state_ip_index(i+noa,x+noa) = istate
            state_ip_index(x+noa,i+noa) = istate
          end do
        end if
      end do 
!:
!:      write(iout,*) " dip_states",ip_states,istate
!:      do i=1, ip_states
!:        write(iout,123) i,hole_ip_index(i,1),hole_ip_index(i,2),part_ip_index(i,1)
!:      end do 
!:      do i=1, noa+nob
!:        do x=1, noa+nob
!:          write(iout,123) i,x,state_ip_index(i,x)
!:        end do 
!:      end do 
!:      flush(iout)
 123  format(4i5)
    
    ip_vec = 0.d0

    ia : do ia=1, ip_states
       xx = hole_ip_index(ia,1) ; x = abs(xx)
       ii = hole_ip_index(ia,2) ; i = abs(ii)
       !: orbital indices
       if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii

       jb : do jb=ia, ip_states
          yy = hole_ip_index(jb,1) ; y = abs(yy)
          jj = hole_ip_index(jb,2) ; j = abs(jj)

          storeme = 0.d0

             !: < ix->a |H| jy->a >
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <jy||ix> = <ij||kl> : DijklAA(k<l,i<j)=DijklAA(i<x,j<y)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAA(j,y,i,x)
                   !: <Jy||Ix> = <yJ||xI> = <iJ||kL> : DijklAB(J<=L,i<=k) = DijklAB(J<=I,y<=x)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklAB(y,J,x,I)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <jY||iX> = <iJ||kL> :
                   !DijklAB(J<=L,i<=k)=DijklAB(Y<=X,j<=i)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAB(j,Y,i,X)
                   !: <JY||IX> = <IJ||KL> : DijklBB(k<l,i<j)=DijklBB(I<X,J<Y)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklBB(J,Y,I,X)
                end if
 
          if ( jb.eq.ia ) storeme = storeme - orben(iorb) - orben(xorb)
          ip_vec((ia-1)*ip_states+jb) = ip_vec((ia-1)*ip_states+jb) + storeme
          if( ia.ne.jb) ip_vec((jb-1)*ip_states+ia) = ip_vec((jb-1)*ip_states+ia) + storeme

       end do jb
    end do ia
 
    info = 10
    lwork = ip_states**2
    allocate( work(lwork) )
    call dsyev('vectors','upper',ip_states,ip_vec,ip_states,ip_eig,work,lwork,info)
    call fix_phase(ip_states,ip_vec)
    deallocate( work )

    write(iout,'(A)') ' ENERGIES for DOUBLY IONIZED STATES'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
    end do
    flush(iout)

    write(iout,'(A)') ' EIGENVECTORS for DOUBLY IONIZED STATES '
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
      do ia=1, ip_states
        ii = ia+ip_states*(i-1)
        if(abs(ip_vec(ii)).gt.1.d-1) write(iout,60) hole_ip_index(ia,1),hole_ip_index(ia,2),ip_vec(ii)
      end do
    end do
    flush(iout)

    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    do i=1, ip_states
       write(100,"(i7,i7,i7)") hole_ip_index(i,1), hole_ip_index(i,2), part_ip_index(i,1)
    end do
    do i=1, noa+nob
       write(100,"(400i5)") (state_ip_index(i,j),j=1,noa+nob) 
    end do
    close(100)    
    
50  format( i5,1x,f15.10 )
60  format( 7x,2i5,2x,4f17.10 )

    call writeme_ham0( 'dip_cis','imdone' )
    
  end subroutine get_dip_cis
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_DIP_SOC
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_dip_soc

    !: assign hole and particle excitation indices for doubly ionized states

    use read_integrals
    use util

    implicit none

    integer(8) :: lscratch, info, lwork
    integer(8) :: i, j, x, y, istate, xstart, xstart1
    integer(8) :: ia, jb, xx, ii, yy, jj, xorb, iorb
    real(8) :: sign, storeme
    real(8), allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    character(20) :: cinfo

    call writeme_ham0( 'dip_soc','imhere' )
    
    part_ip_index(:,:) = 0
    hole_ip_index(:,:) = 0
    state_ip_index = 0
    istate = 0
    xstart  = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
 
!:  alpha,alpha 
        do x=xstart1, noa
          if(xstart.le.x-1) then
            do i=xstart, x-1
              istate = istate + 1 
              hole_ip_index(istate,1) = -x
              hole_ip_index(istate,2) = -i
              state_ip_index(i,x) = istate
              state_ip_index(x,i) = istate
            end do
          end if
       end do 
    end if
 
!:  alpha,beta 
        do x=xstart, noa
          do i=xstart, nob
            istate = istate + 1 
            hole_ip_index(istate,1) = -x
            hole_ip_index(istate,2) = i
            state_ip_index(i+noa,x) = istate
            state_ip_index(x,i+noa) = istate
         end do
       end do 
    
!:  beta,beta
      do x=xstart1, nob
        if(xstart.le.x-1) then
          do i=xstart, x-1
            istate = istate + 1 
            hole_ip_index(istate,1) = x
            hole_ip_index(istate,2) = i
            state_ip_index(i+noa,x+noa) = istate
            state_ip_index(x+noa,i+noa) = istate
          end do
        end if
      end do 

!:      write(iout,*) " dip_states",ip_states,istate
!:      do i=1, ip_states
!:        write(iout,123) i,hole_ip_index(i,1),hole_ip_index(i,2),part_ip_index(i,1)
!:      end do 
!:      do i=1, noa+nob
!:        do x=1, noa+nob
!:          write(iout,123) i,x,state_ip_index(i,x)
!:        end do 
!:      end do 
!: 123  format(5i5)
    
    Zip_vec = dcmplx(0.d0,0.d0)
    call Zform1det_socip(noa,nva,ip_states,Zip_vec,socmoAA, &
       hole_ip_index,part_ip_index,nob,nvb,socmoBB,socmoAB)

    ia : do ia=1, ip_states
       xx = hole_ip_index(ia,1) ; x = abs(xx)
       ii = hole_ip_index(ia,2) ; i = abs(ii)
       !: orbital indices
       if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii

       jb : do jb=ia, ip_states
          yy = hole_ip_index(jb,1) ; y = abs(yy)
          jj = hole_ip_index(jb,2) ; j = abs(jj)

          storeme = 0.d0

             !: < ix->a |H| jy->a >
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <jy||ix> = <ij||kl> : DijklAA(k<l,i<j)=DijklAA(i<x,j<y)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAA(j,y,i,x)
                   !: <Jy||Ix> = <yJ||xI> = <iJ||kL> : DijklAB(J<=L,i<=k) = DijklAB(J<=I,y<=x)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklAB(y,J,x,I)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <jY||iX> = <iJ||kL> :
                   !DijklAB(J<=L,i<=k)=DijklAB(Y<=X,j<=i)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAB(j,Y,i,X)
                   !: <JY||IX> = <IJ||KL> : DijklBB(k<l,i<j)=DijklBB(I<X,J<Y)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklBB(J,Y,I,X)
                end if
 
          if ( jb.eq.ia ) storeme = storeme - orben(iorb) - orben(xorb)
          Zip_vec((ia-1)*ip_states+jb) = Zip_vec((ia-1)*ip_states+jb) + dcmplx(storeme,0.d0)
          if( ia.ne.jb) Zip_vec((jb-1)*ip_states+ia) = Zip_vec((jb-1)*ip_states+ia) + dcmplx(storeme,0.d0)

       end do jb
    end do ia
 
    call testherm(ip_states,Zip_vec,'Zip_vec DIP')
    info = 10
    lwork = ip_states**2
    allocate( work(lwork) )
    allocate( rwork(3*ip_states-2) ) ;  rwork = 0.d0
    call zheev('vectors','upper',ip_states,Zip_vec,ip_states,ip_eig,work,lwork,rwork,info)
    call Zfix_phase(ip_states,Zip_vec)
    deallocate( rwork )
    deallocate( work )

    write(iout,'(A)') ' ENERGIES for DOUBLY IONIZED STATES'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
    end do

    write(iout,'(A)') ' EIGENVECTORS for DOUBLY IONIZED STATES (coeff and |coeff|**2)'
    do i=1, ip_states
      write(iout,50) i, ip_eig(i)*au2eV
      do ia=1, ip_states
        ii = ia+ip_states*(i-1)
        if(abs(Zip_vec(ii)).gt.1.d-1) write(iout,60) hole_ip_index(ia,1),hole_ip_index(ia,2), &
          Zip_vec(ii),dble(dconjg(Zip_vec(ii))*Zip_vec(ii))
      end do
    end do

    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    do i=1, ip_states
       write(100,"(i7,i7,i7)") hole_ip_index(i,1), hole_ip_index(i,2), part_ip_index(i,1)
    end do
    do i=1, noa+nob
       write(100,"(400i5)") (state_ip_index(i,j),j=1,noa+nob) 
    end do
    close(100)    
    
50  format( i5,1x,f15.10 )
60  format( 7x,2i5,2x,4f17.10 )

    call writeme_ham0( 'dip_soc','imdone' )
    
  end subroutine get_dip_soc
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_IP_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_ip_index

    !: assign hole and particle excitation indices for each state
    implicit none

    integer(8) :: i, j, a, b, x, istate, xstart, xstart1


    call writeme_ham0( 'ip_index','imhere' )
    
    part_index(:,:) = 0
    hole_index(:,:) = 0
    istate = 0
    xstart = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
       singlesA : do x=xstart, noa
          istate = istate + 1 
          hole_index(istate,1) = -x
          hole_index(istate,2) = 0
          part_index(istate,1) = 0
       end do singlesA
    end if

    singlesB : do x=xstart, nob
       istate = istate + 1 
       hole_index(istate,1) = x
       hole_index(istate,2) = 0
       part_index(istate,1) = 0
    end do singlesB
    
    if ( IP_alpha_beta ) then
 
!:  alpha,alpha -> alpha
        do x=xstart1, noa
          do i=1, x-1
             do a=1, nva
                istate = istate + 1 
                hole_index(istate,1) = -x
                hole_index(istate,2) = -i
                part_index(istate,1) = -a 
             end do
          end do
       end do 
 
!:  alpha,beta -> beta
        do x=xstart, noa
          do i=1, nob
             do a=1, nvb
                istate = istate + 1 
                hole_index(istate,1) = -x
                hole_index(istate,2) = i
                part_index(istate,1) = a 
             end do
         end do
       end do 
    end if

!:  beta,alpha -> alpha
      do x=xstart, nob
        do i=1, noa
          do a=1, nva
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = -i
             part_index(istate,1) = -a
          end do
        end do 
      end do 
    
!:  beta,beta -> beta
      do x=xstart1, nob
        do i=1, x-1
          do a=1, nvb
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = i
             part_index(istate,1) = a
          end do
        end do
      end do 

!:      do i=1, istate 
!:        write(42,123) i,hole_index(i,1),hole_index(i,2),part_index(i,1) 
!:      end do 
!: 123  format(4i5)
    
    if ( istate.ne.nstates) write(iout,'(A)')  ' ERROR ERROR ERROR IN INDEX ASSIGNMENT IN IP' 
    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    close(100)    

    call writeme_ham0( 'ip_index','imdone' )

    
  end subroutine get_ip_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_SOC_CIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_soc_cisN

    !: get SOC-CIS Hamiltonian

    use read_integrals
    use util
 
    implicit none

    integer(8) :: istate, xorb, iorb, aorb
    integer(8) :: xx, ii, aa, yy, jj, bb, x, y, i, a, j, b
    integer(8) :: ia, jb, ia2, jb2, yx, xb, yj, ji, ix, xi, jy, ab, ya, ib, ja, yb, xa
    real(8) :: sign, storeme
    
    
    call write_header( 'get_soc_cis','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'soc', 'form' ) 

    call get_cis_index 
    
    call get_soc_ao2mo
    call Zform1det_soc(noa,nva,nstates,Zcis_vec,socmoAA,hole_index,part_index,nob,nvb,socmoBB,socmoAB)
     
    ia : do ia=2, nstates
       
       ii = hole_index(ia,1) ; i = abs(ii)
       aa = part_index(ia,1) ; a = abs(aa)
       
       !: orbital indices
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii
       if ( aa.lt.0 ) aorb = -aa + noa ; if ( aa.gt.0 ) aorb = nrorb + nob + aa
       
       jb : do jb=ia, nstates
          
          jj = hole_index(jb,1) ; j = abs(jj)
          bb = part_index(jb,1) ; b = abs(bb)

          storeme = 0.d0

          !: < i->a |H| j->b >
          !: -<ja||ib> = -<ia||jb> : DiajbAA(jb,ia) = DiajbAA(ib,ja)
          if ( ii.lt.0 .and. jj.lt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_diajbAA(j,a,i,b)
          !: -<jA||iB> = -<iA||jB> : DiajbAB(jB,iA) = DiajbAB(iB,jA)
          if ( ii.lt.0 .and. jj.lt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_diajbAB(j,A,i,B)
          !: -<JA||IB> = -<ia||jb> : DiajbBB(jb,ia) = DiajbBB(IB,JA)
          if ( II.gt.0 .and. JJ.gt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_diajbBB(J,A,I,B)
          !: -<Ja||Ib> = -<Ia||Jb> : DiajbBA(jb,ia) = DiajbBA(IB,JA)
          if ( II.gt.0 .and. JJ.gt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_diajbBA(J,a,I,b)
          !: <iJ||aB>  : DijabAB(JB,ia)
          if ( ii.lt.0 .and. JJ.gt.0 .and. aa.lt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_dijabAB(i,J,a,B)
          !: <Ij||Ab> = <jI||bA> : DijabAB(IA,jb)
          if ( II.gt.0 .and. jj.lt.0 .and. AA.gt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_dijabAB(j,I,b,A)

          if ( jb.eq.ia ) storeme = storeme - orben(iorb) + orben(aorb)
          
          istate = (ia-1)*nstates + jb
!:          write(42,"(6i5,4f16.10)") ii,aa,jj,bb,ia,jb,storeme,Zcis_vec(istate)
          Zcis_vec(istate) = Zcis_vec(istate) + dcmplx(storeme,0.d0)

          istate = (jb-1)*nstates + ia
          if(ia.ne.jb) Zcis_vec(istate) = Zcis_vec(istate) + dcmplx(storeme,0.d0)
         
       end do jb
    end do ia
!:     Dump Hamiltonian for analysis with Mathematica
!:     call math_ham
   
    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start    

    call write_header( 'get_soc_cis','get_ham0','leave' )
    

  end subroutine get_soc_cisN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_SOCIP_CISD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_socip_cisd

    !: get SOCIP-CISD Hamiltonian

    use read_integrals
    use util
 
    implicit none

    integer(8) :: istate, xorb, iorb, aorb
    integer(8) :: xx, ii, aa, yy, jj, bb, x, y, i, a, j, b
    integer(8) :: ia, jb, ia2, jb2, yx, xb, yj, ji, ix, xi, jy, ab, ya, ib, ja, yb, xa
    real(8) :: sign, storeme
    
    
    call write_header( 'get_socip_cisd','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'socip', 'form' ) 

    call get_socip_index 
    
    call get_soc_ao2mo
    call Zform1det_socip(noa,nva,nstates,Zcis_vec,socmoAA,hole_index,part_index,nob,nvb,socmoBB,socmoAB)
     
    ia : do ia=1, nstates
       
       xx = hole_index(ia,1) ; x = abs(xx)
       ii = hole_index(ia,2) ; i = abs(ii)
       aa = part_index(ia,1) ; a = abs(aa)
       
       !: orbital indices
       if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii
       if ( aa.lt.0 ) aorb = -aa + noa ; if ( aa.gt.0 ) aorb = nrorb + nob + aa
       
       jb : do jb=ia, nstates
          
          yy = hole_index(jb,1) ; y = abs(yy)
          jj = hole_index(jb,2) ; j = abs(jj)
          bb = part_index(jb,1) ; b = abs(bb)

          storeme = 0.d0
          
          SS : if( ii.eq.0 .and. jj.eq.0 ) then
             if( xx.eq.yy ) storeme = -orben(xorb)
             go to 78
          end if SS
          
          
          SD : if ( ii*jj .eq. 0 ) then

             if ( ii.eq.0 ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then                   
                   !: <x|H|jb_y> = -<yj||xb> = -<ij||ka> : DijkaAA(ka,i<j) = DijkaAA(xb,y<j)
                   if ( jj.lt.0 .and. bb.lt.0 ) storeme = - get_dijkaAA(y,j,x,b)
                   !: <x|H|JB_y> = -<yJ||xB> = -<iJ||kA> : DijkaAB(i<=k,JA) = DijkaAB(y<=x,JB)
                   if ( JJ.gt.0 .and. BB.gt.0 ) storeme = - get_dijkaAB(y,J,x,B)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <X|H|jb_Y> = -<Yj||Xb> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(jb2,Y<=X)
                   if ( jj.lt.0 .and. bb.lt.0 ) storeme = - get_dijkaBA(Y,j,X,b)
                   !: <X|H|JB_Y> = -<YJ||XB> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(XB,Y<J)
                   if ( JJ.gt.0 .and. BB.gt.0 ) storeme = - get_dijkaBB(Y,J,X,B)
                end if
             else if ( jj.eq.0 ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <y|H|ia_x> = -<xi||ya> = -<ij||ka> :  DijkaAA(ka,i<j) = DijkaBB(ya,x<i) 
                   if ( ii.lt.0 .and. aa.lt.0 ) storeme = - get_dijkaAA(x,i,y,a)
                   !: <y|H|IA_x> = -<xI||yA> = -<iJ||kA> : DijkaAB(JA,i<=k) = DijkaAB(IA,x<=y)
                   if ( II.gt.0 .and. AA.gt.0 ) storeme = - get_dijkaAB(x,I,y,A)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <Y|H|ia_X> = -<Xi||Ya> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(ia, X<=Y)
                   if ( ii.lt.0 .and. aa.lt.0 ) storeme = - get_dijkaBA(X,i,Y,a)
                   !: <Y|H|IA_X> = -<XI||YA> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(YA,X<I)
                   if ( II.gt.0 .and. AA.gt.0 ) storeme = - get_dijkaBB(X,I,Y,A)
                end if
             end if

             go to 78
             
          end if SD
          
          
          DD : if( ii.ne.0 .and. jj.ne.0 ) then
             
             !: < ix->a |H| jy->a >
             kdelta_ab : if ( aa.eq.bb ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: <jy||ix> = <ij||kl> : DijklAA(k<l,i<j)=DijklAA(i<x,j<y)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAA(j,y,i,x)
                   !: <Jy||Ix> = <yJ||xI> = <iJ||kL> : DijklAB(J<=L,i<=k) = DijklAB(J<=I,y<=x) 
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklAB(y,J,x,I)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: <jY||iX> = <iJ||kL> : DijklAB(J<=L,i<=k)=DijklAB(Y<=X,j<=i)
                   if ( ii.lt.0 .and. jj.lt.0 ) storeme = storeme + get_dijklAB(j,Y,i,X)
                   !: <JY||IX> = <IJ||KL> : DijklBB(k<l,i<j)=DijklBB(I<X,J<Y)
                   if ( II.gt.0 .and. JJ.gt.0 ) storeme = storeme + get_dijklBB(J,Y,I,X)
                end if
             end if kdelta_ab
 
             
             !: < ix->a |H| iy->b >
             kdelta_ij : if ( ii.eq.jj ) then
                if ( xx.lt.0 .and. yy.lt.0 ) then
                   !: -<ya||xb> = <ia||jb> : DiajbAA(jb,ia) = DiajbAA(xb,ya)
                   if ( aa.lt.0 .and. bb.lt.0 ) storeme = storeme - get_diajbAA(y,a,x,b)
                   !: -<yA||xB> = <iA||jB> : DiajbAB(A<=B,i<=j) = DiajbAB(A<=B,y<=x)
                   if ( AA.gt.0 .and. BB.gt.0 ) storeme = storeme - get_diajbAB(y,A,x,B)
                else if ( XX.gt.0 .and. YY.gt.0 ) then
                   !: -<Ya||Xb> = <Ia||Jb> : DiajbBA(a<=b,I<=J)=DiajbBA(a<=b,Y<=X)
                   if ( aa.lt.0 .and. bb.lt.0 ) storeme = storeme - get_diajbBA(Y,a,X,b)
                   !: -<YA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,YA)                   
                   if ( AA.gt.0 .and. BB.gt.0 ) storeme = storeme - get_diajbBB(Y,A,X,B)
                end if
             end if kdelta_ij
             
             
             !: < ix->a |H| jx->b >
             kdelta_xy : if ( XX.eq.YY ) then
                !: -<ja||ib> = -<ia||jb> : DiajbAA(jb,ia) = DiajbAA(ib,ja)
                if ( ii.lt.0 .and. jj.lt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_diajbAA(j,a,i,b)
                !: -<jA||iB> = -<iA||jB> : DiajbAB(jB,iA) = DiajbAB(iB,jA)
                if ( ii.lt.0 .and. jj.lt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_diajbAB(j,A,i,B)
                !: -<JA||IB> = -<ia||jb> : DiajbBB(jb,ia) = DiajbBB(IB,JA)
                if ( II.gt.0 .and. JJ.gt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_diajbBB(J,A,I,B)
                !: -<Ja||Ib> = -<Ia||Jb> : DiajbBA(jb,ia) = DiajbBA(IB,JA)
                if ( II.gt.0 .and. JJ.gt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_diajbBA(J,a,I,b)
                !: <iJ||aB>  : DijabAB(JB,ia)
                if ( ii.lt.0 .and. JJ.gt.0 .and. aa.lt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_dijabAB(i,J,a,B)
                !: <Ij||Ab> = <jI||bA> : DijabAB(IA,jb)
                if ( II.gt.0 .and. jj.lt.0 .and. AA.gt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_dijabAB(j,I,b,A)
             end if kdelta_xy


             !: < ix->a |H| xy->b >
             kdelta_jx : if ( JJ.eq.XX ) then
                !: < ix->a |H| xy->b > :  <ya||ib> = <ia||jb> : DiajbAA(jb,ia)=DiajbAA(ib,ya)
                if ( ii.lt.0 .and. yy.lt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_diajbAA(y,a,i,b)
                !: < ix->A |H| xy->B > :  <yA||iB> = <iA||jB> : DiajbAB(jb,ia)=DiajbAB(ib,ya)
                if ( ii.lt.0 .and. yy.lt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_diajbAB(y,A,i,B)
                !: < Ix->A |H| xy->b > : - <Iy||Ab> = -<yI||bA>
                if ( II.gt.0 .and. yy.lt.0 .and. AA.gt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_dijabAB(y,I,b,A)
                !: < IX->A |H| XY->B > :  <YA||IB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(IB,YA)
                if ( II.gt.0 .and. YY.gt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_diajbBB(Y,A,I,B)
                !: < IX->a |H| XY->b > :  <Ya||Ib> = <Ia||Jb> : DiajbBA(Jb,Ia)=DiajbBA(Ib,Ya)
                if ( II.gt.0 .and. YY.gt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_diajbBA(Y,a,I,b)
                !: < iX->a |H| XY->B > : -<iY||aB>
                if ( ii.lt.0 .and. YY.gt.0 .and. aa.lt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_dijabAB(i,Y,a,B)                   
             end if kdelta_jx
             

             !: < yx->a |H| jy->b >
             kdelta_iy : if ( II.eq.YY ) then
                !: < yx->a |H| jy->b > :  <ja||xb> = <ia||jb> : DiajbAA(jb,ia)=DiajbAA(xb,ja)
                if ( jj.lt.0 .and. xx.lt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_diajbAA(j,a,x,b)
                !: < yx->A |H| jy->B > :  <jA||xB> = <iA||jB> : DiajbAB(jB,iA)=DiajbAB(xB,jA)
                if ( jj.lt.0 .and. xx.lt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_diajbAB(j,A,x,B)
                !: < yx->a |H| Jy->B > :  - <xJ||aB> 
                if ( JJ.gt.0 .and. xx.lt.0 .and. aa.lt.0 .and. BB.gt.0 ) &
                                   storeme = storeme - get_dijabAB(x,J,a,B)
                !: < YX->A |H| JY->B > :  <JA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,JA)
                if ( JJ.gt.0 .and. XX.gt.0 .and. AA.gt.0 .and. BB.gt.0 ) &
                                   storeme = storeme + get_diajbBB(J,A,X,B)
                !: < YX->a |H| JY->b > :  <Ja||Xb> = <Ia||Jb> : DiajbBA(Jb,Ia)=DiajbBA(Xb,Ja)
                if ( JJ.gt.0 .and. XX.gt.0 .and. aa.lt.0 .and. bb.lt.0 ) &
                                   storeme = storeme + get_diajbBA(J,a,X,b)
                !: < YX->A |H| jY->b > :  - <jX||bA>
                if ( jj.lt.0 .and. XX.gt.0 .and. AA.gt.0 .and. bb.lt.0 ) &
                                   storeme = storeme - get_dijabAB(j,X,b,A)                   
             end if kdelta_iy                         
             
          end if DD

          if ( jb.eq.ia ) storeme = storeme - orben(iorb) - orben(xorb) + orben(aorb)
          
78        continue
          
          istate = (ia-1)*nstates + jb
!:          write(42,123) ii,aa,xx,jj,bb,yy,ia,jb,storeme,Zcis_vec(istate)
 123 format(8i5,4f16.10)
          Zcis_vec(istate) = Zcis_vec(istate) + dcmplx(storeme,0.d0)

          istate = (jb-1)*nstates + ia
          if(ia.ne.jb) Zcis_vec(istate) = Zcis_vec(istate) + dcmplx(storeme,0.d0)
         
       end do jb
    end do ia
!:     Dump Hamiltonian for analysis with Mathematica
!:     call math_ham
   
    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start    

    call write_header( 'get_socip_cisd','get_ham0','leave' )
    

  end subroutine get_socip_cisd
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_SOCIP_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_socip_index

    !: assign hole and particle excitation indices for each state
    implicit none

    integer(8) :: i, j, a, b, x, istate, xstart, xstart1, natrim, nbtrim


    call writeme_ham0( 'socip_index','imhere' )
    
    part_index(:,:) = 0
    hole_index(:,:) = 0
    if(nactive.gt.0) then
      natrim = noa - nactive
      nbtrim = noa - nactive
    else
      natrim = 0
      nbtrim = 0
    end if
    istate = 0
    xstart = 1
    If(nactive.gt.0) xstart = noa+1 - nactive
    xstart1 = max(2,xstart)

    !: only beta electrons will be ionized by default

    if ( IP_alpha_beta ) then
       singlesA : do x=xstart, noa
          istate = istate + 1 
          hole_index(istate,1) = -x
          hole_index(istate,2) = 0
          part_index(istate,1) = 0
       end do singlesA
    end if

    singlesB : do x=xstart, nob
       istate = istate + 1 
       hole_index(istate,1) = x
       hole_index(istate,2) = 0
       part_index(istate,1) = 0
    end do singlesB
    
    if ( IP_alpha_beta ) then
 
!:  alpha,alpha -> alpha
        do x=xstart1, noa
          do i=1, x-1
             do a=1, nva
                istate = istate + 1 
                hole_index(istate,1) = -x
                hole_index(istate,2) = -i
                part_index(istate,1) = -a 
             end do
          end do
       end do 
 
!:  alpha,beta -> beta
        do x=1, noa
          do i=1, nob
            if(x.gt.natrim .or. i.gt.nbtrim) then
             do a=1, nvb
                istate = istate + 1 
                hole_index(istate,1) = -x
                hole_index(istate,2) = i
                part_index(istate,1) = a 
             end do
           end if
         end do
       end do 
    end if
    
!:  beta,alpha -> alpha
      do x=1, nob
        do i=1, noa
          if(x.gt.nbtrim .or. i.gt.natrim) then
           do a=1, nva
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = -i
             part_index(istate,1) = -a
          end do
         end if
        end do 
      end do 
    
!:  beta,beta -> beta
      do x=xstart1, nob
        do i=1, x-1
          do a=1, nvb
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = i
             part_index(istate,1) = a
          end do
        end do
      end do 

      if(QsocA2B) then

!:    alpha,alpha -> beta
      if ( IP_alpha_beta ) then
        do x=xstart1, noa
          do i=1, x-1
             do a=1, nvb
                istate = istate + 1 
                hole_index(istate,1) = -x
                hole_index(istate,2) = -i
                part_index(istate,1) = a 
             end do
          end do
        end do 
      end if
 
!:    beta,beta -> alpha
      do x=xstart1, nob
        do i=1, x-1
          do a=1, nva
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = i
             part_index(istate,1) = -a
          end do
        end do
      end do 

      end if
 
!:      do i=1, istate 
!:        write(42,123) i,hole_index(i,1),hole_index(i,2),part_index(i,1) 
!:      end do 
!: 123  format(4i5)
    
    if ( istate.ne.nstates) write(iout,'(A)')  ' ERROR ERROR ERROR IN INDEX ASSIGNMENT IN SOCIP' 
    open(unit=100,file='INDEX.OUT')
    do i=1, nstates
       write(100,"(i7,i7,i7)") hole_index(i,1), hole_index(i,2), part_index(i,1)
    end do
    close(100)    

    call writeme_ham0( 'socip_index','imdone' )
    
  end subroutine get_socip_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DIAGONALIZE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine diagonalize

    !: call lapack dysev to get eigenstates and eigenvalues    
    
    use util
    implicit none
    
    call write_header( 'diagonalize', 'get_ham0', 'enter' )    

    call cpu_time(start)
    
    call writeme_ham0( 'diagonalize' , 'calling' )

    call eigen_mat(iout,nstates,cis_vec,nstates,cis_eig,QeigenDC)

    call cpu_time(finish)

    write(iout,"(' LAPACK diagonalization time:',f12.4,' seconds')") finish-start
    flush(iout)
 
    call write_header( 'diagonalize','get_ham0','leave' )
    
  end subroutine diagonalize
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ZDIAGONALIZE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine Zdiagonalize

    !: call lapack dysev to get eigenstates and eigenvalues    
    
    use util
    implicit none
    
    call write_header( 'diagonalize', 'get_ham0', 'enter' )    
    
    call cpu_time(start)

    call testherm(nstates,Zcis_vec,'Ham')
      
    call writeme_ham0( 'diagonalize' , 'calling' )

    call Zeigen_mat(iout,nstates,Zcis_vec,nstates,cis_eig,QeigenDC)
      
    call cpu_time(finish)
    
    write(iout,"(' LAPACK diagonalization time:',f12.4,' seconds')") finish-start
    flush(iout)
 
    call write_header( 'diagonalize','get_ham0','leave' )
    
  end subroutine Zdiagonalize
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_HAM0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_ham0(myroutine,option)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option
    

    select case( trim(myroutine) )
    case( 'cis' ) 
       select case( trim(option) )
       case('form')    ; write(iout,"(' forming ',i0,' x ',i0,' CIS Hamiltonian')") nstates, nstates
       end select
    case( 'cis_index' )
       select case( trim(option) ) 
       case( 'imhere' ) ; write(iout,"(' assigning hole_index(nstates,1) and part_index(nstates,1) for CIS')")
       case( 'imdone' ) ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'ip' )
       select case( trim(option) )
       case('form')    ; write(iout,"(' forming ',i0,' x ',i0,' IP-CISD Hamiltonian')") nstates, nstates
       end select
    case( 'ip_index' )
       select case( trim(option) )
       case( 'imhere' ) ; write(iout,"(' assigning hole_index(nstates,2) and part_index(nstates,1) for IP-CISD')")
       case( 'imdone' ) ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'cisd_index' )
       select case( trim(option) )
       case('imhere') ; write(iout,"(' assigning hole_index(nstates,2) and part_index(nstates,2) for CISD')")
       case('imdone') ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'diagonalize' )
       select case( trim(option) )
       case('calling')  ; write(iout,'(A)') ' calling LAPACK diagonalization'
       end select
    end select

    
    flush(iout)


  end subroutine writeme_ham0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ERRORS_GETHAM0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_getham0(myroutine,option,mystring)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option
    character(*), optional,intent(in) :: mystring

    
    select case( trim(myroutine) )
    case( 'diagonalize' )
       select case( trim(option) )
       case( 'illegal' ) 
          write(iout,'(A)') ' ERROR:  illegal value in '//trim(adjustl(mystring))//' ham0 element'
          go to 100
       case('no_converge')
          write(iout,'(A)') ' ERROR:  dysev failed to converge!  INFO = '//trim(adjustl(mystring))
          go to 100
       case('success')
          write(iout,'(A)') ' INFO = '//trim(adjustl(mystring))
          go to 200 
       end select
    end select
    
    
100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " SOMEDAY YOU'LL SEE THE REASON WHY THERE'S GOOD IN GOODBYE - unknown"
    stop


200 continue

    
  end subroutine errors_getham0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module get_ham0
