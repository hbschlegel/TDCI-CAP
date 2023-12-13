module get_ham0_cisd
  
  use global_variables
  implicit none


  integer(8) :: double_startAA, double_startBB, double_startAB, double_startBA


contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_0_IJAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_0(sigma)

    use read_integrals
    implicit none
    
    real(8), intent(inout) :: sigma(nstates)

    integer(8) :: i, j, a, b, ia
    

    sigma = 0.d0

    
    !: [H]0 = iJ --> aB  <iJ||aB> 
    do i=1, noa
       do J=1, NOB
          ia = state_index(-i,j) - 1 
          do a=1, nva
             do B=1, NVB
                ia = ia + 1
                sigma(ia) = get_dijabAB(i,J,a,B)
             end do
          end do
       end do
    end do

    
    !: [H]0 = i<j --> a<b   <ij||ab> = <ib||aj> = -<ib||ja>    
    ia = double_startAA - 1
    do i=1, noa
       do j=(i+1), noa
          do a=1, nva
             do b=(a+1), nva
                ia = ia + 1
                sigma(ia) = get_dijabAA(i,j,a,b)
             end do
          end do
       end do
    end do

    
    !: [H]0 = I<J --> A<B  <IJ||AB> = <IB||AJ> = -<IB||JA>
    do I=1, NOB
       do J=(I+1), NOB
          do A=1, NVB
             do B=(A+1), NVB
                ia = ia + 1
                sigma(ia) = get_dijabBB(I,J,A,B)
             end do
          end do
       end do
    end do

    
    close(100)

  end subroutine get_cisd_0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CIS_ia_alpha
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_ia_AA(i_in,a_in,sigma)
    

    use read_integrals
    implicit none

    integer(8), intent(in) :: i_in, a_in
    real(8), intent(inout) :: sigma(nstates)
    
    real(8) :: rdum
    integer(8) :: ia, jb, jb2, ij
    integer(8) :: i, j, k, l, a, b, c, d

    
    sigma = 0.d0
    
    i = i_in
    a = a_in
    
    
    !: singles
    jb = 1 
    
    
    !: k --> c  -<ka||ic>   ;  k=i will take care of later
    do k=1, noa
       do c=1, nva 
          jb = jb + 1
          sigma(jb) = - get_diajbAA(k,a,i,c)
       end do
    end do
    
    !:  K  --> C  <iK||aC>
    do K=1, NOB
       do C=1, NVB
          jb = jb + 1 
          sigma(jb) = get_dijabAB(i,K,a,C)
       end do
    end do
    
    
    !: doubles 

    
    !: k < l  -->  c < a  ;  k < l  -->  a < c 
    do k=1, noa
       do l=(k+1), noa
          ij = state_index(-k,-l) - 1
          !: k < l  -->  c < a  ; k=i, l=i will take care of later
          do c=1, (a-1)
             jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
             sigma(jb) = get_dijkaAA(k,l,i,c)
          end do
          !: k < l  -->  a < c
          jb = state_index(-k,-l) - 1 + (a-1)*nva - a*(a-1)/2
          do c=(a+1), nva
             jb = jb + 1
             sigma(jb) = - get_dijkaAA(k,l,i,c)
          end do
       end do
    end do
    
                           
    
    !: k < i  -->  c < a  ;  k < i  -->  a < c 
    do k=1, (i-1)
       ij = state_index(-k,-i) - 1 
       !: k < i  -->  c < a 
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = get_dijkaAA(k,i,i,c) + get_diabcAA(k,a,c,a)
       end do
       !: k < i  -->  a < c 
       jb = state_index(-k,-i) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva  
          jb = jb + 1 
          sigma(jb) = - get_dijkaAA(k,i,i,c) - get_diabcAA(k,a,c,a)
       end do
    end do


    !: i < k  -->  c < a   ;  i < k  -->  a < c 
    do k=(i+1), noa
       ij = state_index(-i,-k) - 1
       !: i < k  -->  c < a
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = - get_dijkaAA(k,i,i,c) - get_diabcAA(k,a,c,a)
       end do
       !:  i < k  -->  a < c
       jb = state_index(-i,-k) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva 
          jb = jb + 1 
          sigma(jb) = get_dijkaAA(k,i,i,c) + get_diabcAA(k,a,c,a)
       end do
    end do    


    !: k < i  -->  c < d 
    do k=1, (i-1) 
       do c=1, nva
          if ( c.ne.a ) then             
             jb = state_index(-k,-i) - 1 + (c-1)*nva - c*(c-1)/2 
             do d=(c+1), nva
                jb = jb + 1 
                if ( d.ne.a ) sigma(jb) = get_diabcAA(k,a,c,d)
             end do
          end if
       end do
    end do

    !: i < k  -->  c < d 
    do k=(i+1), noa
       do c=1, nva 
          if ( c.ne.a ) then
             jb  = state_index(-i,-k) - 1 + (c-1)*nva - c*(c-1)/2
             do d=(c+1), nva 
                jb = jb + 1 
                if ( d.ne.a ) sigma(jb) = - get_diabcAA(k,a,c,d)
             end do
          end if
       end do
    end do
    
    
    !: alpha beta doubles


    !: kL  -->  aD  ; k=i will take care of later
    do k=1, noa
       do L=1, NOB
          jb = state_index(-k,L) + (a-1)*nvb - 1 
          do D=1, NVB 
             jb = jb + 1 
             sigma(jb)  = - get_dijkaAB(k,L,i,D)
          end do
       end do
    end do
    
    
    !: iL  -->  cD  ;  c=a will take care of later
    do L=1, NOB
       do c=1, nva 
          jb = state_index(-i,L) - 1 + (c-1)*NVB
          do D=1, NVB 
             jb = jb + 1 
             sigma(jb)  = get_diabcBA(L,a,D,c)
          end do
       end do
    end do
    
    
    !: iL  -->  aD 
    do L=1, nob
       jb = state_index(-i,L) - 1 + (a-1)*NVB
       do D=1, NVB 
          jb  = jb + 1 
          sigma(jb)  = - get_dijkaAB(i,L,i,D) + get_diabcBA(L,a,D,a)
       end do
    end do
    
    flush(iout)

    !: DIAGONAL  


    ia = (i-1)*nva + a + 1 
    rdum = - orben(i) + orben(noa+a)
    rdum = rdum - get_diajbAA(i,a,i,a)
    sigma(ia) = rdum

    
  end subroutine get_cisd_ia_AA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CIS_ia_beta
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_ia_BB(i_in,a_in,sigma)
    

    use read_integrals
    implicit none

    integer(8), intent(in) :: i_in, a_in
    real(8), intent(inout) :: sigma(nstates)
    
    real(8) :: rdum
    integer(8) :: ia, jb, jb2, ij
    integer(8) :: i, j, k, l, a, b, c, d
    

    sigma = 0.d0
    
    i = i_in
    a = a_in
    
    
    !: singles
    jb = 1 
    
    
    !: k  --> c  < kI || cA >
    do k=1, noa
       do c=1, nva
          jb = jb + 1 
          sigma(jb) = get_dijabAB(k,I,c,A)
       end do
    end do
    
    !:  K  --> C  < iK || aC >
    do K=1, NOB
       do C=1, NVB
          jb = jb + 1 
          sigma(jb) = - get_diajbBB(K,A,I,C)
       end do
    end do
    
    
    !: doubles 
        
    
    !: K < L  -->  C < A  ;  K < L  -->  A < C
    do K=1, NOB
       do L=(K+1), NOB
          ij = state_index(K,L) - 1 
          !: K < L  -->  C < A  ;  K=I, L=I will take care of later
          do C=1, (A-1)
             jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
             sigma(jb) = get_dijkaBB(K,L,I,C)
          end do
          !: K < L  -->  A < C
          jb = state_index(K,L) - 1 + (A-1)*NVB - A*(A-1)/2
          do C=(A+1), NVB
             jb = jb + 1 
             sigma(jb) = - get_dijkaBB(K,L,I,C)
          end do
       end do
    end do
    
    
    !: K < I  -->  C < A  ;  K < I  -->  A < C 
    do K=1, (I-1)
       ij = state_index(K,I) - 1 
       !: K < I  -->  C < A 
       do C=1, (A-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = get_dijkaAA(K,I,I,C) + get_diabcBB(K,A,C,A)
       end do
       !: K < I  -->  A < C 
       jb = state_index(K,I) - 1 + (A-1)*NVB - A*(A-1)/2
       do C=(A+1), NVB
          jb = jb + 1 
          sigma(jb) = - get_dijkaBB(K,I,I,C) - get_diabcBB(K,A,C,A)
       end do
    end do
    

    !: I < K  -->  C < A   ;  I < K  -->  A < C 
    do K=(I+1), NOB
       ij = state_index(I,K) - 1  
       !: I < K  -->  C < A
       do C=1, (A-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = - get_dijkaBB(K,I,I,C) - get_diabcBB(K,A,C,A)
       end do
       !:  I < K  -->  A < C
       jb = state_index(I,K) - 1 + (A-1)*NVB - A*(A-1)/2
       do C=(A+1), NVB
          jb = jb + 1 
          sigma(jb) = get_dijkaBB(K,I,I,C) + get_diabcBB(K,A,C,A)
       end do
    end do

    !: K < I  -->  C < D  
    do K=1, (I-1) 
       do C=1, NVB
          if ( C.ne.A ) then
             jb = state_index(K,I) - 1 + (C-1)*NVB - C*(C-1)/2
             do D=(C+1), NVB
                jb = jb + 1 
                if ( D.ne.A ) sigma(jb) = get_diabcBB(K,A,C,D)
             end do
          end if
       end do
    end do
    
    !: I < K  -->  C < D 
    do K=(I+1), NOB
       do C=1, NVB
          if ( C.ne.A ) then
             jb = state_index(I,K) - 1 + (C-1)*NVB - C*(C-1)/2
             do D=(C+1), NVB
                jb = jb + 1 
                if ( D.ne.A ) sigma(jb) = - get_diabcBB(K,A,C,D)
             end do
          end if
       end do
    end do   
    

    !: alpha beta doubles
    
    
    !: lK -->  dA  ;  K=I will take care of later
    do l=1, noa 
       do K=1, NOB
          do d=1, nva
             jb = state_index(-l,K) + (d-1)*NVB + A - 1 
             sigma(jb) = - get_dijkaBA(K,l,I,d)
          end do
       end do
    end do

    !:  lI  -->  dC  ;  C=A will take care of later
    do l=1, noa
       do d=1, nva
          jb = state_index(-l,I) - 1 + (d-1)*NVB
          do C=1, NVB
             jb = jb + 1 
             !jb = state_index(-l,I) + (d-1)*NVB + C - 1 
             sigma(jb) = get_diabcAB(l,A,d,C)
          end do
       end do
    end do
    
    
    !: lI  -->  dA
    do l=1, noa 
       do d=1, nva
          jb = state_index(-l,I) + (d-1)*NVB + A - 1 
          sigma(jb) = - get_dijkaBA(I,l,I,d) + get_diabcAB(l,A,d,A)
       end do
    end do
    

    !: DIAGONAL  
    
    
    ia = (I-1)*NVB + A + 1 + noanva
    rdum = - orben(nrorb+I) + orben(nrorb+NOB+A)
    rdum = rdum - get_diajbBB(I,A,I,A)
    sigma(ia) = rdum
    
    
    
  end subroutine get_cisd_ia_BB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_ijab_alpha_alpha
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_ijab_AAAA(i_in, j_in, a_in, b_in, sigma)

    use read_integrals
    implicit none

    integer(8), intent(in) :: i_in, j_in, a_in, b_in
    real(8), intent(inout) :: sigma(nstates)

    integer(8) :: ia, jb, jb2, ij
    integer(8) :: i, j, k, l, a, b, c, d
    real(8) :: rdum
    
    
    sigma = 0.d0
    i = min( i_in, j_in )
    j = max( i_in, j_in )
    a = min( a_in, b_in )
    b = max( a_in, b_in )
    
    
    !: with ground state
    sigma(1) = get_dijabAA(i,j,a,b)
    
    !: single
    

    !: i->c
    jb = state_index(-i,0) - 1 
    do c=1, nva
       jb = jb + 1
       sigma(jb) = - get_diabcAA(j,c,a,b)
    end do
    
    !: j->c
    jb = state_index(-j,0) - 1
    do c=1, nva
       jb = jb + 1
       sigma(jb) = get_diabcAA(i,c,a,b)
    end do
    
    !: k->a  ;  k=i k=j will take care of later
    do k=1, noa
       jb = (k-1)*nva + a + 1 
       sigma(jb) = - get_dijkaAA(i,j,k,b) 
    end do
    
    !: k->b  ;  k=i k=j will take care of later
    do k=1, noa
       jb = (k-1)*nva + b + 1
       sigma(jb) = get_dijkaAA(i,j,k,a)
    end do
    

    !: i->a 
    jb = (i-1)*nva + a + 1 
    sigma(jb) = - get_dijkaAA(i,j,i,b) + get_diabcAA(j,a,b,a)

    !: i->b
    jb = (i-1)*nva + b + 1 
    sigma(jb) = get_dijkaAA(i,j,i,a) - get_diabcAA(j,b,a,b)

    !: j->a 
    jb = (j-1)*nva + a + 1 
    sigma(jb) = - get_dijkaAA(i,j,j,b) - get_diabcAA(i,a,b,a)
    
    !: j->b
    jb = (j-1)*nva + b + 1 
    sigma(jb) = get_dijkaAA(i,j,j,a) + get_diabcAA(i,b,a,b)


    
    !: double i < j  -->  c < d   ;  c,d = a,b dont worry will take care of later
    jb = state_index(-i,-j) - 1 
    do c=1, nva
       do d=(c+1), nva
          jb = jb + 1 
          sigma(jb) = get_dabcdAA(a,b,c,d)
       end do
    end do          
    

    !: double k < l  -->  a < b  ;  k,l = i,j don't worry will take care of later
    do k=1, noa
       do l=(k+1), noa
          jb = state_index(-k,-l) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a) 
          sigma(jb) = get_dijklAA(k,l,i,j)
       end do
    end do
    
    
    !: double i, a    
    
    
    do k=1, (i-1)
       !: k < i  -->  c < a 
       ij = state_index(-k,-i) - 1 
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = - get_diajbAA(k,b,j,c)
       end do
       !: k < i  -->  a < c  ; c=b will take care of later
       jb = state_index(-k,-i) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva
          jb = jb + 1
          sigma(jb) = get_diajbAA(k,b,j,c)
       end do
    end do



    do k=(i+1), noa
       !: i < k  -->  c < a   ;  k=j will take care of later
       ij = state_index(-i,-k) - 1
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = get_diajbAA(k,b,j,c)
       end do
       !: i < k  -->  a < c   ;  c=b will take care of later
       jb = state_index(-k,-i) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva
          jb = jb + 1 
          sigma(jb) = - get_diajbAA(k,b,j,c)
       end do
    end do
    
    
    !: double i, b
    

    do k=1, (i-1) 
       !: k < i  -->  c < b  ;  c=a will take care of later
       ij = state_index(-k,-i) - 1
       do c=1, (b-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (b-c)
          sigma(jb) = get_diajbAA(k,a,j,c)
       end do
       !: k < i  -->  b < c
       jb = state_index(-k,-i) - 1 + (b-1)*nva - b*(b-1)/2
       do c=(b+1), nva
          jb = jb + 1 
          sigma(jb) = - get_diajbAA(k,a,j,c)
       end do
    end do
    
    
    do k=(i+1), noa
       !: i < k  -->  c < b   ; k=j will take care of later  ;  c=a will take care of later
       ij = state_index(-i,-k) - 1 
       do c=1, (b-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (b-c)
          sigma(jb) = - get_diajbAA(k,a,j,c)
       end do
       !: i < k  ->  b < c
       jb = state_index(-i,-k) - 1 + (b-1)*nva - b*(b-1)/2
       do c=(b+1), nva
          jb = jb + 1 
          sigma(jb) = get_diajbAA(k,a,j,c)
       end do
    end do

    
    !: double j, a 
    

    do k=1, (j-1)
       !: k < j  -->  c < a   ;  k=i will take care of later
       ij = state_index(-k,-j) - 1 
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = get_diajbAA(k,b,i,c)
       end do
       !: k < j  -->  a < c   ;  c=b will take care of later
       jb = state_index(-k,-j) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva
          jb = jb + 1 
          sigma(jb) = - get_diajbAA(k,b,i,c)
       end do
    end do
    

    do k=(j+1), noa
       !: j < k  -->  c < a 
       ij = state_index(-j,-k) - 1 
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = - get_diajbAA(k,b,i,c)
       end do
       !: j < k  -->  a < c   ;  c=b will take care of later
       jb = state_index(-j,-k) - 1 + (a-1)*nva - a*(a-1)/2 
       do c=(a+1), nva
          jb = jb + 1 
          sigma(jb) = get_diajbAA(k,b,i,c)
       end do
    end do
    

    !: double j, b
    

    do k=1, (j-1)
       !: k < j  -->  c < b   ;  k=i, c=a will take care of later 
       ij = state_index(-k,-j) - 1 
       do c=1, (b-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (b-c)
          sigma(jb) = - get_diajbAA(k,a,i,c)
       end do
       !: k < j  -->  b < c 
       jb = state_index(-k,-j) - 1 + (b-1)*nva - b*(b-1)/2
       do c=(b+1), nva
          jb = jb + 1
          sigma(jb) = get_diajbAA(k,a,i,c)
       end do
    end do

    
    do k=(j+1), noa
       !: j < k  -->  c < b  ;  c=a will take care of later
       ij = state_index(-j,-k) - 1 
       do c=1, (b-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (b-c)
          sigma(jb) = get_diajbAA(k,a,i,c)
       end do
       !: j < k  -->  b < c 
       jb = state_index(-j,-k) - 1 + (b-1)*nva - b*(b-1)/2
       do c=(b+1), nva
          jb = jb + 1 
          sigma(jb) = - get_diajbAA(k,a,i,c)
       end do
    end do
    
    
    !: doubles alpha beta

    
    !: iK  -->  aC , Ki  -->  Ca
    do K=1, NOB
       jb = state_index(-i,K) - 1 + (a-1)*NVB
       do C=1, NVB
          jb = jb + 1 
          sigma(jb) = get_dijabAB(j,K,b,C)
       end do
    end do
    
    
    !: jK  -->  bC  ,  Kj  -->  Cb
    do K=1, NOB
       jb = state_index(-j,K) - 1 + (b-1)*NVB
       do C=1, NVB
          jb = jb + 1
          sigma(jb) = get_dijabAB(i,K,a,C)
       end do
    end do
    
    
    !: iK  -->  bC  ,  Ki  -->  Cb
    do K=1, NOB
       jb = state_index(-i,K) - 1 + (b-1)*NVB
       do C=1, NVB
          jb = jb + 1 
          sigma(jb) = - get_dijabAB(j,K,a,C)
       end do
    end do
    
    
    !: jK  -->  aC  ,  Kj  -->  Ca
    do K=1, NOB
       jb = state_index(-j,K) - 1 + (a-1)*NVB
       do C=1, NVB
          jb = jb + 1 
          sigma(jb) = - get_dijabAB(i,K,b,C)
       end do
    end do
    
    
    !: doubles k<i  -->  a<b  ,  k<i  -->  a<b

    
    !: k < i  --> a < b
    do k=1, (i-1)
       jb = state_index(-k,-i) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a)
       rdum = - get_dijklAA(k,i,j,i)
       rdum = rdum + get_diajbAA(k,a,j,a)
       rdum = rdum + get_diajbAA(k,b,j,b)
       sigma(jb) = rdum
    end do

    
    !: i < k  -->  a < b  ;  k=j will take care of later
    do k=(i+1), noa
       jb = state_index(-i,-k) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a)
       rdum = get_dijklAA(k,i,j,i)
       rdum = rdum - get_diajbAA(k,a,j,a)
       rdum = rdum - get_diajbAA(k,b,j,b)
       sigma(jb) = rdum
    end do

    
    !: doubles k<j  --> a<b  ,  j<k  --> a<b
    

    !: k < j  -->  a < b  ;  k=i will take care of later
    do k=1, (j-1)
       jb = state_index(-k,-j) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a)
       rdum = get_dijklAA(k,j,i,j)
       rdum = rdum - get_diajbAA(k,a,i,a)
       rdum = rdum - get_diajbAA(k,b,i,b)
       sigma(jb) = rdum
    end do

    
    !:  j < k  -->  a < b
    do k=(j+1), noa
       jb = state_index(-j,-k) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a)
       rdum = - get_dijklAA(k,j,i,j)
       rdum = rdum + get_diajbAA(k,a,i,a)
       rdum = rdum + get_diajbAA(k,b,i,b)
       sigma(jb) = rdum
    end do

    
    !: doubles i<j -->  c<a  ,  i<j  -->  a<c
    
    
    !:  i < j  -->  c < a 
    ij = state_index(-i,-j) - 1 
    do c=1, (a-1)       
       jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
       rdum = get_diajbAA(i,b,i,c)
       rdum = rdum + get_diajbAA(j,b,j,c)
       rdum = rdum - get_dabcdAA(b,a,c,a)
       sigma(jb) = rdum
    end do

    
    !:  i < j  -->  a < c   ;  c=b will take care of later
    ij = state_index(-i,-j) - 1 
    do c=(a+1), nva
       jb = ij + (a-1)*nva - a*(a-1)/2 + (c-a)
       rdum = - get_diajbAA(i,b,i,c)
       rdum = rdum - get_diajbAA(j,b,j,c)
       rdum = rdum + get_dabcdAA(b,a,c,a)
       sigma(jb) = rdum
    end do

    
    !: doubles i<j  --> c<b  ,  i<j  -->  b<c
    
    !: i < j  --> c < b    ;  c=a will take care of later
    ij = state_index(-i,-j) - 1 
    do c=1, (b-1)
       jb = ij + (c-1)*nva - c*(c-1)/2 + (b-c)
       rdum = - get_diajbAA(i,a,i,c)
       rdum = rdum - get_diajbAA(j,a,j,c)
       rdum = rdum + get_dabcdAA(a,b,c,b)
       sigma(jb) = rdum
    end do

    
    !: i < j  -->  b < c
    ij = state_index(-i,-j) - 1 
    do c=(b+1), nva
       jb = ij + (b-1)*nva - b*(b-1)/2 + (c-b)
       rdum = get_diajbAA(i,a,i,c)
       rdum = rdum + get_diajbAA(j,a,j,c)
       rdum = rdum - get_dabcdAA(a,b,c,b)
       sigma(jb) = rdum
    end do
    
    
    !: *** DIAGONAL  *** :!
    

    rdum = - orben(i) - orben(j) + orben(noa+a) + orben(noa+b)
    rdum = rdum + get_dijklAA(i,j,i,j) + get_dabcdAA(a,b,a,b)
    rdum = rdum - get_diajbAA(i,a,i,a) - get_diajbAA(i,b,i,b)
    rdum = rdum - get_diajbAA(j,a,j,a) - get_diajbAA(j,b,j,b)    
    
    ia = state_index(-i,-j) - 1 + (a-1)*nva - a*(a-1)/2 + (b-a)
    sigma(ia) = rdum
           
    

  end subroutine get_cisd_ijab_AAAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_ijab_beta_beta
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_ijab_BBBB(i_in, j_in, a_in, b_in, sigma)

    use read_integrals
    implicit none

    integer(8), intent(in) :: i_in, j_in, a_in, b_in
    real(8), intent(inout) :: sigma(nstates)

    integer(8) :: ia, jb, jb2, ij
    integer(8) :: i, j, k, l, a, b, c, d
    real(8) :: rdum
    
    
    sigma = 0.d0
    i = min( i_in, j_in )
    j = max( i_in, j_in )
    a = min( a_in, b_in )
    b = max( a_in, b_in )

    
    !: with ground state
    sigma(1) = get_dijabBB(i,j,a,b)   
    
    
    !: single
    

    !: I->C
    jb = state_index(I,0) - 1 
    do C=1, NVB
       jb = jb + 1
       sigma(jb) = - get_diabcBB(J,C,A,B)
    end do
    
    !: J->C
    jb = state_index(J,0) - 1 
    do C=1, NVB
       jb = jb + 1 
       sigma(jb) = get_diabcBB(I,C,A,B)
    end do
    
    !: K->A  ;  K!=J will take care of later
    do K=1, (I-1)
       jb = (K-1)*NVB + A + 1 + noanva
       sigma(jb) = - get_dijkaBB(I,J,K,B) 
    end do
    do K=(I+1), NOB
       jb = (K-1)*NVB + A + 1 + noanva
       sigma(jb) = - get_dijkaBB(I,J,K,B)
    end do

    !: K->B  ;  K!=J will take care of later
    do K=1, (I-1)
       jb = (K-1)*NVB + B + 1 + noanva
       sigma(jb) = get_dijkaBB(I,J,K,A)
    end do
    do K=(I+1), NOB
       jb = (K-1)*NVB + B + 1 + noanva
       sigma(JB) = get_dijkaBB(I,J,K,A)
    end do
    
    !: I->A 
    jb = (I-1)*NVB + A + 1 + noanva
    sigma(jb) = - get_dijkaBB(I,J,I,B) + get_diabcBB(J,A,B,A)

    !: I->B
    jb = (I-1)*NVB + B + 1 + noanva
    sigma(jb) = get_dijkaBB(I,J,I,A) - get_diabcBB(J,B,A,B)

    !: J->A 
    jb = (J-1)*NVB + A + 1 + noanva
    sigma(jb) = - get_dijkaBB(I,J,J,B) - get_diabcBB(I,A,B,A)

    !: J->B
    jb = (J-1)*NVB + B + 1 + noanva
    sigma(jb) = get_dijkaBB(I,J,J,A) + get_diabcBB(I,B,A,B)
    
    

    !: double I < J  -->  C < D  ;  C=A, D=B will take care of later
    jb = state_index(I,J) - 1 
    do C=1, NVB
       do D=(C+1), NVB
          jb = jb + 1 
          sigma(jb) = get_dabcdBB(A,B,C,D)
       end do
    end do
    

    !: double K < L  -->  A < B  ;  K=I, L=I will take care of later
    do K=1, NOB
       do L=(K+1), NOB
          jb = state_index(K,L) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
          sigma(jb) = get_dijklBB(K,L,I,J)
       end do
    end do
    
    
    
    !: double i, a    
    
    
    do K=1, (I-1)
       !: K < I  -->  C < A 
       ij = state_index(K,I) - 1 
       do c=1, (a-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = - get_diajbBB(K,B,J,C)
       end do
       !: K < I  -->  A < C  ;  C=B will take care of later
       jb = state_index(K,I) - 1 + (A-1)*NVB - A*(A-1)/2 
       do C=(A+1), NVB
          jb = jb + 1
          sigma(jb) = get_diajbBB(K,B,J,C)
       end do
    end do

    do K=(I+1), NOB
       !: I < K  -->  C < A  ;  K=J will take care of later
       ij = state_index(I,K) - 1
       do C=1, (A-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = get_diajbBB(K,B,J,C)
       end do
       !: I < K  -->  A < C  ;  C=B will take care of later
       jb = state_index(I,K) - 1 + (A-1)*NVB - A*(A-1)/2
       do C=(A+1), NVB
          jb = jb + 1 
          sigma(jb) = - get_diajbBB(K,B,J,C)
       end do
    end do
    
    
    !: double i, b
    

    do K=1, (I-1) 
       !: K < I  -->  C < B  ;  C=A will take care of later
       ij = state_index(K,I) - 1
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = get_diajbBB(K,A,J,C)
       end do
       !: K < I  -->  B < C
       jb = state_index(K,I) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), NVB
          jb = jb + 1 
          sigma(jb) = - get_diajbBB(K,A,J,C)
       end do
    end do
    
    do K=(I+1), NOB
       !: I < K  -->  C < B  ;  K=J, C=A will take care of later
       ij = state_index(I,K) - 1 
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = - get_diajbBB(K,A,J,C)
       end do
       !: I < K  ->  B < C
       jb = state_index(I,K) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), NVB
          jb = jb + 1 
          sigma(jb) = get_diajbBB(K,A,J,C)
       end do
    end do
    

    !: double j, a 


    do K=1, (J-1)
       !: K < J  -->  C < A  ;  K=I will take care of later
       ij = state_index(K,J) - 1 
       do C=1, (A-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = get_diajbBB(K,B,I,C)
       end do
       !: K < J  -->  A < C  ;  C=B will take care of later
       jb = state_index(K,J) - 1 + (A-1)*NVB - A*(A-1)/2
       do C=(A+1), NVB
          jb = jb + 1 
          sigma(jb) = - get_diajbBB(K,B,I,C)
       end do
    end do
    
    do K=(J+1), NOB
       !: J < K  -->  C < A 
       ij = state_index(J,K) - 1 
       do C=1, (A-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C)
          sigma(jb) = - get_diajbBB(K,B,I,C)
       end do       
       !: J < K  -->  A < C  ;  C=B will take care of later
       jb = state_index(J,K) - 1 + (A-1)*NVB - A*(A-1)/2
       do C=(A+1), NVB
          jb = jb + 1 
          sigma(jb) = get_diajbBB(K,B,I,C)
       end do
    end do
    

    !: double j, b
    

    do K=1, (J-1)
       !: K < J  -->  C < B  ;  K=I, C=A will take care of later
       ij = state_index(K,J) - 1 
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = - get_diajbBB(K,A,I,C)
       end do
       !: K < J  -->  B < C 
       jb = state_index(K,J) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), NVB
          jb = jb + 1 
          sigma(jb) = get_diajbBB(K,A,I,C)
       end do
    end do
    
    do K=(J+1), NOB
       !: J < K  -->  C < B  ;  C=A will take care of later
       ij = state_index(J,K) - 1 
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = get_diajbBB(K,A,I,C)
       end do
       !: J < K  -->  B < C 
       jb = state_index(J,K) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), NVB
          jb = jb + 1 
          sigma(jb) = - get_diajbBB(K,A,I,C)
       end do
    end do
    
    
    !: doubles alpha beta

    
    !: kI  -->  cA
    do k=1, noa
       ij  = state_index(-k,I) - 1 
       do c=1, nva
          jb = ij + (c-1)*NVB + A 
          sigma(jb) = get_dijabAB(k,J,c,B)
       end do
    end do
    
    !: kJ  -->  cB
    do k=1, noa
       ij = state_index(-k,J) - 1 
       do c=1, nva
          jb = ij + (c-1)*NVB + B
          sigma(jb) = get_dijabAB(k,I,c,A)
       end do
    end do

    !: kI  -->  cB
    do k=1, noa
       ij = state_index(-k,I) - 1
       do c=1, nva
          jb = ij + (c-1)*NVB + B
          sigma(jb) = - get_dijabAB(k,J,c,A)
       end do
    end do
    
    !: kJ  -->  cA
    do k=1, noa
       ij = state_index(-k,J) - 1 
       do c=1, nva
          jb = ij + (c-1)*NVB + A
          sigma(jb) = - get_dijabAB(k,I,c,B)
       end do
    end do
    
    
    !: doubles k<i  -->  a<b  ,  k<i  -->  a<b

    !: K < I  --> A < B
    do K=1, (I-1)
       jb = state_index(K,I) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
       rdum = - get_dijklBB(K,I,J,I)
       rdum = rdum + get_diajbBB(K,A,J,A)
       rdum = rdum + get_diajbBB(K,B,J,B)
       sigma(jb) = rdum
    end do

    !: I < K  -->  A < B  ;  K=J will take care of later
    do K=(I+1), NOB
       jb = state_index(I,K) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
       rdum = get_dijklBB(K,I,J,I)
       rdum = rdum - get_diajbBB(K,A,J,A)
       rdum = rdum - get_diajbBB(K,B,J,B)
       sigma(jb) = rdum
    end do


    !: doubles k<j  --> a<b  ,  j<k  --> a<b
    

    !: K < J  -->  A < B  ;  K=I will take care of later
    do K=1, (J-1)
       jb = state_index(K,J) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
       rdum = get_dijklBB(K,J,I,J)
       rdum = rdum - get_diajbBB(K,A,I,A)
       rdum = rdum - get_diajbBB(K,B,I,B)
       sigma(jb) = rdum
    end do
    
    !:  J < K  -->  A < B
    do K=(J+1), NOB
       jb = state_index(J,K) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
       rdum = - get_dijklBB(K,J,I,J)
       rdum = rdum + get_diajbBB(K,A,I,A)
       rdum = rdum + get_diajbBB(K,B,I,B)
       sigma(jb) = rdum
    end do


    !: doubles i<j -->  c<a  ,  i<j  -->  a<c
    
    
    !:  I < J  -->  C < A 
    ij = state_index(I,J)
    do C=1, (A-1)       
       jb = ij + (C-1)*NVB - C*(C-1)/2 + (A-C) - 1 
       rdum = get_diajbBB(I,B,I,C)
       rdum = rdum + get_diajbBB(J,B,J,C)
       rdum = rdum - get_dabcdBB(B,A,C,A)
       sigma(jb) = rdum
    end do

    !:  I < J  -->  A < C  ;  C=B will take care of later
    do C=(A+1), NVB
       jb = ij + (A-1)*NVB - A*(A-1)/2 + (C-A) - 1 
       rdum = - get_diajbBB(I,B,I,C)
       rdum = rdum - get_diajbBB(J,B,J,C)
       rdum = rdum + get_dabcdBB(B,A,C,A)
       sigma(jb) = rdum
    end do
    

    !: doubles i<j  --> c<b  ,  i<j  -->  b<c

    
    !: I < J  --> C < B  ;  C=A will take care of later
    ij = state_index(I,J)
    do C=1, (B-1)
       jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C) - 1 
       rdum = - get_diajbBB(I,A,I,C)
       rdum = rdum - get_diajbBB(J,A,J,C)
       rdum = rdum + get_dabcdBB(A,B,C,B)
       sigma(jb) = rdum
    end do
    
    !: I < J  -->  B < C
    do C=(B+1), NVB
       jb = ij + (B-1)*NVB - B*(B-1)/2 + (C-B) - 1 
       rdum = get_diajbBB(I,A,I,C)
       rdum = rdum + get_diajbBB(J,A,J,C)
       rdum = rdum - get_dabcdBB(A,B,C,B)
       sigma(jb) = rdum
    end do
    

    !: *** DIAGONAL  *** :!
    

    rdum = - orben(nrorb+I) - orben(nrorb+J) + orben(nrorb+NOB+A) + orben(nrorb+NOB+B)
    rdum = rdum + get_dijklBB(I,J,I,J) + get_dabcdBB(A,B,A,B)
    rdum = rdum - get_diajbBB(I,A,I,A) - get_diajbBB(I,B,I,B)
    rdum = rdum - get_diajbBB(J,A,J,A) - get_diajbBB(J,B,J,B)    
    
    ia = state_index(I,J) + (A-1)*NVB - A*(A-1)/2 + (B-A) - 1 
    sigma(ia) = rdum
    
    

  end subroutine get_cisd_ijab_BBBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_IJAB_ALPHA_BETA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_ijab_ABAB( i_in, j_in, a_in, b_in, sigma )
    
    use read_integrals
    implicit none

    integer(8), intent(in) :: i_in, j_in, a_in, b_in
    real(8), intent(inout) :: sigma(nstates)

    integer(8) :: jb, ij
    integer(8) :: i, j, k, l, a, b, c, d
    real(8) :: rdum
    

    sigma = 0.d0
    i = i_in !: alpha
    j = j_in !: beta
    a = a_in !: alpha
    b = b_in !: beta

    
    sigma(1) = get_dijabAB(i,J,a,B)


    !: singles
    
    
    !:  k --> a  ;  k=i will take care of later
    do k=1, noa
       jb = (k-1)*nva + a + 1 
       sigma(jb) = - get_dijkaAB(i,J,k,B)
    end do

    !: i --> c  ;  c=a will take care of later
    jb = (i-1)*nva + 1 
    do c=1, nva
       jb = jb + 1 
       sigma(jb) = get_diabcBA(J,c,B,a)
    end do

    !: K --> B  ;  K=J will take care of later
    do K=1, NOB
       jb = 1 + noanva + (K-1)*NVB + B
       sigma(jb) = - get_dijkaBA(J,i,K,a)
    end do
    
    !: J --> C  ;  C=B will take care of later
    jb = 1 + noanva + (J-1)*NVB
    do C=1, NVB
       jb = jb + 1 
       sigma(jb) = get_diabcAB(i,C,a,B)
    end do
    
    !: i --> a 
    jb = (i-1)*nva + a + 1 
    rdum = - get_dijkaAB(i,J,i,B) + get_diabcBA(J,a,B,a)
    sigma(jb) = rdum

    !: J --> B
    jb = noanva + (J-1)*NVB + B + 1 
    rdum = - get_dijkaBA(J,i,J,a) + get_diabcAB(i,B,a,B)
    sigma(jb) = rdum
    

    
    !: doubles
    

    !: k<i  -->  a<c  ;  k<i  -->  c<a 
    do k=1, (i-1)
       !: k<i --> c<a 
       ij = state_index(-k,-i) - 1
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = get_dijabAB(k,J,c,B)
       end do
       !:  k<i -->  a<c 
       jb = state_index(-k,-i) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva
          jb = jb + 1
          sigma(jb) = - get_dijabAB(k,J,c,B)
       end do
    end do
          
    !: i<k  -->  a<c  ;  i<k  -->  c<a 
    do k=(i+1), noa
       ij = state_index(-i,-k) - 1 
       !: i<k  -->  c<a
       do c=1, (a-1)
          jb = ij + (c-1)*nva - c*(c-1)/2 + (a-c)
          sigma(jb) = - get_dijabAB(k,J,c,B)
       end do
       !: i<k  -->  a<c
       jb = state_index(-i,-k) - 1 + (a-1)*nva - a*(a-1)/2
       do c=(a+1), nva
          jb = jb + 1 
          sigma(jb) = get_dijabAB(k,J,c,B)
       end do
    end do

    
    !: K<J  -->  C<B  ;  K<J  -->  B<C
    do K=1, (J-1)
       ij = state_index(K,J) - 1 
       !: K<J  -->  C<B
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = get_dijabAB(i,K,a,C)
       end do
       !: K<J  -->  B<C
       jb = state_index(K,J) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), nvb
          jb = jb + 1 
          sigma(jb) = - get_dijabAB(i,K,a,C)
       end do
    end do
       
    !: J<K  -->  C<B  ; J<K --> B<C
    do K=(J+1), NOB
       ij = state_index(J,K) - 1 
       !: J<K  -->  C<B
       do C=1, (B-1)
          jb = ij + (C-1)*NVB - C*(C-1)/2 + (B-C)
          sigma(jb) = - get_dijabAB(i,K,a,C)
       end do
       !: J<K --> B<C
       jb = state_index(J,K) - 1 + (B-1)*NVB - B*(B-1)/2
       do C=(B+1), NVB
          jb = jb + 1 
          sigma(jb) = get_dijabAB(i,K,a,C)
       end do
    end do


    !: iJ -->  cD 
    jb = state_index(-i,J) - 1 
    do c=1, nva
       do D=1, NVB
          jb = jb + 1 
          sigma(jb) = get_dabcdAB(a,B,c,D)
       end do
    end do


    !: kL -->  aB
    do k=1, noa
       do L=1, NOB
          jb = state_index(-k,L) - 1 + (a-1)*NVB + B
          sigma(jb) = get_dijklAB(k,L,i,J)
       end do
    end do


    !: iK -->  aC  ;  K=J, C=B will take care of later
    do K=1, NOB
       jb = state_index(-i,K) + (a-1)*NVB - 1 
       do C=1, NVB
          jb = jb + 1          
          sigma(jb) = - get_diajbBB(K,B,J,C)
       end do
    end do

    !: iK --> cB  
    do K=1, NOB
       ij = state_index(-i,K) - 1 
       do c=1, nva
          jb = ij + (c-1)*NVB + B
          sigma(jb) = - get_diajbBA(K,a,J,c)
       end do
    end do


    !: kJ -->  cB  ;  k=i, c=a will take care of later
    do k=1, noa
       ij = state_index(-k,J) - 1 
       do c=1, nva
          jb = ij + (c-1)*NVB + B 
          sigma(jb) = - get_diajbAA(k,a,i,c)
       end do
    end do
    

    !: iK  -->  aC  ;  K=J, C=B will take care of later
    do K=1, NOB
       jb = state_index(-i,K) + (a-1)*NVB - 1 
       do C=1, NVB
          jb = jb + 1
          rdum = - get_diajbBB(K,B,J,C)
          sigma(jb) = rdum
       end do
    end do
       

    !: kJ  -->  cB  ;  k=i, c=a will take care of later
    do k=1, noa
       ij = state_index(-k,J) - 1
       do c=1, nva
          jb = ij + (c-1)*NVB + B 
          sigma(jb) = - get_diajbAA(k,a,i,c)
       end do
    end do


    !: kJ  -->  aC
    do k=1, noa
       ij = state_index(-k,J) - 1
       do C=1, NVB
          jb = ij + (a-1)*NVB + C
          sigma(jb) = - get_diajbAB(k,B,i,C)
       end do
    end do

    !: iJ  -->  aC  ;  C=B will take care of later
    do C=1, NVB
       jb  = state_index(-i,J) + (a-1)*NVB + C - 1 
       rdum = - get_diajbAB(i,B,i,C)
       rdum = rdum - get_diajbBB(J,B,J,C)
       rdum = rdum + get_dabcdAB(a,B,a,C)
       sigma(jb)  = rdum
    end do


    !: iK  -->  aB  ;  K=J will take care of later
    do K=1, nob
       jb  = state_index(-i,K) + (a-1)*NVB + B - 1 
       rdum = get_dijklAB(i,K,i,J) 
       rdum = rdum - get_diajbBA(K,a,J,a)
       rdum = rdum - get_diajbBB(K,B,J,B)
       sigma(jb)  = rdum
    end do


    !: iJ  -->  cB  ;  c=a will take care of later
    do c=1, nva
       jb  = state_index(-i,J) + (c-1)*NVB + B - 1 
       rdum = - get_diajbAA(i,a,i,c)
       rdum = rdum - get_diajbBA(J,a,J,c)
       rdum = rdum + get_dabcdAB(a,B,c,B)
       sigma(jb)  = rdum
    end do
    
    
    !: kJ  -->  aB  ;  k=i will take care of later
    do k=1, noa
       jb  = state_index(-k,J) + (a-1)*NVB + B - 1 
       rdum = get_dijklAB(k,J,i,J)
       rdum = rdum - get_diajbAA(k,a,i,a)
       rdum = rdum - get_diajbAB(k,B,i,B)
       sigma(jb)  = rdum
    end do
              

    !: diagonal

    jb = state_index(-i,J) + (a-1)*NVB + B - 1 
    rdum = - orben(i) - orben(nrorb+J) + orben(noa+a) + orben(nrorb+NOB+B)
    rdum = rdum + get_dijklAB(i,J,i,J) + get_dabcdAB(a,B,a,B)
    rdum = rdum - get_diajbBA(J,a,J,a) - get_diajbAB(i,B,i,B)
    rdum = rdum - get_diajbAA(i,a,i,a) - get_diajbBB(J,B,J,B)
    sigma(jb) = rdum

    
  end subroutine get_cisd_ijab_ABAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_DIAGONAL
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_diagonal( sigma )

    use read_integrals
    implicit none

    real(8), parameter :: scf = -75.9835736457
    real(8), intent(inout) :: sigma(nstates)
    
    real(8)    :: rdum
    integer(8) :: ia, i, a, j, b
    

    sigma = 0.d0
    ia = 1 
    
    !: alpha singles 
    open( unit=100,file='diagonals.out' )    

    write( 100, '(i5,f20.10)') 1, scf
    
    do i=1, noa
       do a=1, nva
          ia = ia + 1 
          sigma(ia) = - orben(i) + orben(noa+a) - get_diajbAA(i,a,i,a)
          write(100,'(5(i5,1x),f20.10)') ia, -i, 0, -a, 0, sigma(ia)
       end do
    end do
    
    !: beta singles 
    do i=1, nob
       do a=1, nvb
          ia = ia + 1 
          sigma(ia) = - orben(nrorb+i) + orben(nrorb+nob+a) - get_diajbBB(i,a,i,a)
          write(100,'(5(i5,1x),f20.10)') ia, i, 0, a, 0, sigma(ia)
       end do
    end do


    !: alpha beta doubles
    do i=1, noa
       do J=1, NOB
          do a=1, nva
             do B=1, NVB                
                ia = ia + 1 
                rdum = - orben(i) - orben(nrorb+J) + orben(noa+a) + orben(nrorb+NOB+B)
                rdum = rdum + get_dijklAB(i,J,i,J) + get_dabcdAB(a,B,a,B)
                rdum = rdum - get_diajbBA(J,a,J,a) - get_diajbAB(i,B,i,B)
                rdum = rdum - get_diajbAA(i,a,i,a) - get_diajbBB(J,B,J,B)
                sigma(ia) = rdum
                write(100,'(5(i5,1x),f20.10)') ia, -i, j, -a, b, sigma(ia)
             end do
          end do
       end do
    end do

    
    !: alpha doubles
    do i=1, noa
       do j=(i+1), noa
          do a=1, nva
             do b=(a+1), nva
                ia = ia + 1 
                rdum = - orben(i) - orben(j) + orben(noa+a) + orben(noa+b)
                rdum = rdum + get_dijklAA(i,j,i,j) + get_dabcdAA(a,b,a,b)
                rdum = rdum - get_diajbAA(i,a,i,a) - get_diajbAA(j,a,j,a)
                rdum = rdum - get_diajbAA(i,b,i,b) - get_diajbAA(j,b,j,b)
                sigma(ia) = rdum
                write(100,'(5(i5,1x),f20.10)') ia, -i, -j, -a, -b, sigma(ia)
             end do
          end do
       end do
    end do


    !: beta doubles
    do I=1, NOB
       do J=(I+1), NOB
          do A=1, NVB
             do B=(A+1), NVB
                ia = ia + 1 
                rdum = - orben(nrorb+I) - orben(nrorb+J) + orben(nrorb+NOB+A) + orben(nrorb+NOB+B)
                rdum = rdum + get_dijklBB(I,J,I,J) + get_dabcdBB(A,B,A,B)
                rdum = rdum - get_diajbBB(I,A,I,A) - get_diajbBB(I,B,I,B)
                rdum = rdum - get_diajbBB(J,A,J,A) - get_diajbBB(J,B,J,B)    
                sigma(ia) = rdum
                write(100,'(5(i5,1x),f20.10)') ia, i, j, a, b, sigma(ia)
             end do
          end do
       end do
    end do


    close(100)
    

  end subroutine get_cisd_diagonal
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CISD_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cisd_index

    implicit none

    integer(8) :: i, j, a, b, ia


    
    ia = 1 !: ground state
    hole_index = 0
    part_index = 0


    allocate( state_index(-noa:nob,-noa:nob) ) 

    !: temp for debugging
    allocate( cisd_indices(-nva:nvb, -nva:nvb, -noa:nob, -noa:nob) ) 

    
    open(unit=100,file='OUT_CISD_INDICES')

    
    !: singles

    !: alpha
    do i=1, noa
       state_index(-i,0) = ia + 1 
       state_index(0,-i) = ia + 1 
       do a=1, nva
          ia = ia + 1 
          write(100,"(5(i5))") ia, -i, 0, -a, 0
          cisd_indices(-a,0,-i,0) = ia 
          cisd_indices(0,-a,0,-i) = ia 
       end do
    end do
    
    !: beta
    do i=1, nob
       state_index(i,0) = ia + 1 
       state_index(0,i) = ia + 1 
       do a=1, nvb
          ia = ia + 1 
          write(100,"(5(i5))") ia, i, 0, a, 0
          cisd_indices(a,0,i,0) = ia 
          cisd_indices(0,a,0,i) = ia 
       end do
    end do
    
    !: doubles


    !: alpha beta
    double_startAB = ia + 1
    do i=1, noa
       do j=1, nob
          state_index(-i,j) = ia + 1 
          state_index(j,-i) = ia + 1 
          do a=1, nva
             do b=1, nvb
                ia = ia + 1 
                write(100,"(5(i5))") ia, -i, j, -a, b
                cisd_indices(-a,b,-i,j) = ia
                cisd_indices(b,-a,j,-i) = ia 
             end do
          end do
       end do
    end do
    

    !: alpha alpha
    double_startAA = ia + 1 
    do i=1, noa
       do j=i+1, noa
          state_index(-i,-j) = ia + 1 
          state_index(-j,-i) = ia + 1 
          do a=1, nva
             do b=a+1, nva
                ia = ia + 1 
                write(100,"(5(i5))") ia, -i, -j, -a, -b
                cisd_indices(-b,-a,-j,-i) = ia 
                cisd_indices(-a,-b,-i,-j) = ia 
             end do
          end do
       end do
    end do
    
    !: beta beta
    double_startBB = ia + 1 
    do i=1, nob
       do j=i+1, nob
          state_index(i,j) = ia + 1 
          state_index(j,i) = ia + 1 
          do a=1, nvb
             do b=a+1, nvb
                ia = ia + 1 
                write(100,"(5(i5))") ia, i, j, a, b
                cisd_indices(b,a,j,i) = ia 
                cisd_indices(a,b,i,j) = ia                 
             end do
          end do
       end do
    end do

    close(100)
   
    write(iout,'(A)') ' '     
    write(iout,'(4(a5,1x),a15)') 'i', 'j', 'a', 'b', 'index'
    
    !: singles    
    do i=1, noa
       write(iout,'(4(i5,1x),i15)') -i, 0, -1, 0, state_index(-i,0)
    end do
    do i=1, nob
       write(iout,'(4(i5,1x),i15)') i, 0, 1, 0, state_index(i,0)
    end do

    !: doubles alpha alpha
    do i=1, noa
       do j=(i+1), noa
          write(iout,'(4(i5,1x),i15)') -i, -j, -1, -2, state_index(-i,-j)
       end do
    end do
    !: beta beta
    do i=1, nob
       do j=(i+1), nob
          write(iout,'(4(i5,1x),i15)') i, j, 1, 2, state_index(i,j)
       end do
    end do
    !: alpha beta
    do i=1, noa
       do j=1, nob
          write(iout,'(4(i5,1x),i15)') -i, j, -1, 1, state_index(-i,j)
       end do
    end do
    write(iout,'(A)') ''    
    

    
  end subroutine get_cisd_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine write_cisd_hamiltonian( sigma, iout2 )

    implicit none

    integer(8), intent(in) :: iout2
    real(8), intent(in) :: sigma(nstates)

    integer(8) :: i, j, a, b, ia


    !write(iout2,'(i5,f20.10)') 1, sigma(1)
    write(iout2,'(f20.7)') sigma(1)
    
    !: alpha singles
    do i=1, noa
       do a=1, nva
          ia = cisd_indices(-a,0,-i,0)
          !write(iout2,"(3(i5),f20.10)") ia, -i, -a, sigma(ia)
          write(iout2,"(f20.7)") sigma(ia)
       end do
    end do

    !: beta singles
    do i=1, nob
       do a=1, nvb
          ia = cisd_indices(a,0,i,0)
          !write(iout2,"(3(i5),f20.10)") ia, i, a, sigma(ia)
          write(iout2,"(f20.7)") sigma(ia)
       end do
    end do

    !: alpha beta double
    do i=1, noa 
       do j=1, nob
          do a=1, nva
             do b=1, nvb
                ia = cisd_indices(-a,b,-i,j)
                !write(iout2,"(5(i5),f20.10)") ia, -i,j,-a,b, sigma(ia)
                write(iout2,"(f20.7)") sigma(ia)
             end do
          end do
       end do
    end do
    
    !: alpha alpha double
    do i=1, noa 
       do j=(i+1), noa
          do a=1, nva
             do b=(a+1), nva
                ia = cisd_indices(-a,-b,-i,-j)
                !write(iout2,"(5(i5),f20.10)") ia, -i,-j,-a,-b, sigma(ia)
                write(iout2,"(f20.7)") sigma(ia)
             end do
          end do
       end do
    end do
    
    !: beta beta double
    do i=1, nob 
       do j=(i+1), nob
          do a=1, nvb
             do b=(a+1), nvb
                ia = cisd_indices(a,b,i,j)
                !write(iout2,"(5(i5),f20.10)") ia,i,j,a,b, sigma(ia)
                write(iout2,"(f20.7)") sigma(ia)
             end do
          end do
       end do
    end do
    


  end subroutine write_cisd_hamiltonian
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module get_ham0_cisd

