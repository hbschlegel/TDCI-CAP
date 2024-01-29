module util

  
  implicit none
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DET00
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine det00(noa,nva,mo_matAA,myval,nob,nvb,mo_matBB)

    ! <C> get < 0 | 1e_operator | 0 > matrix element

    implicit none    
    
    real(8), intent(inout) :: myval
    integer(8), intent(in) :: noa, nva
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva)
    ! <C> if unrestricted 
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)

    
    integer(8) :: i 
    
    
    myval = 0.d0
    do i=1, noa 
       myval = myval + mo_matAA(i,i)
    end do
    
    if ( present(nob) ) then
       do i=1, nob
          myval = myval + mo_matBB(i,i)
       end do
    end if
       
        
  end subroutine det00
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE AO2MO
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine ao2mo(nbasis,nrorb,ao_mat,mo_mat,rotmat)  !nrorbB,mo_matBB,rotmatBB)
    

    ! <C> Transform one electron matrix from AO (lower triangle) to MO (full) 
    ! <C> mo_mat = [rotmat] [ao_mat] [rotmat]^T

    implicit none
    
    integer(8), intent(in) :: nbasis, nrorb
    real(8),    intent(in) :: ao_mat(nbasis*nbasis)
    real(8),    intent(in) :: rotmat(nbasis,nrorb)
    real(8), intent(inout) :: mo_mat(nrorb,nrorb)
    !integer(8), optional, intent(in) :: nrorbB
    !real(8), optional, intent(in)    :: rotmatBB(nbasis,nrorbB)
    !real(8), optional, intent(inout) :: mo_matBB(nrorbB,nrorbB)
    
    
    integer(8) :: iao, jao, kmo, lmo, ij, i
    real(8)    :: aodum, rotmatdum, rdum
    real(8)    :: tmp(nrorb,nbasis)
    
    
    mo_mat = 0.d0  
    tmp    = 0.d0
    
    ! <C> ao is stored as a 1D array, consisting of the lower half of the square matrix
    do iao=1, nbasis    
       do jao=1, nbasis 
          
          ij  = iao*(iao-1)/2 + jao
          if (jao.gt.iao) ij = jao*(jao-1)/2 + iao

          aodum = ao_mat(ij)
          tmp(:,jao) = tmp(:,jao) + rotmat(iao,:) * aodum
          
       end do
    end do
    
    mo_mat = matmul( tmp,rotmat )
    
    
  end subroutine ao2mo
  !: ---------------------------!
  !: SUBROUTINE AO2MO_COMPLEX   !
  !: ---------------------------!
  subroutine ao2mo_complex(nbasis,nrorb,ao_mat,mo_mat,rotmatAA,rotmatBB)  
    

    !: used for SOC ao-->mo transformation
    
    implicit none
    
    integer(8), intent(in) :: nbasis, nrorb
    complex(8), intent(in) :: ao_mat(nbasis,nbasis)
    real(8),    intent(in) :: rotmatAA(nbasis,nrorb), rotmatBB(nbasis,nrorb)
    complex(8), intent(inout) :: mo_mat(nrorb,nrorb)    
    
    integer(8) :: iao, jao, kmo, lmo, ij, i
    complex(8) :: aodum, cdum
    complex(8) :: tmp(nrorb,nbasis)
    
    
    mo_mat = dcmplx( 0.d0, 0.d0 )
    tmp    = dcmplx( 0.d0, 0.d0 )
    
    ! <C> ao is stored as a 1D array, consisting of the lower half of the square matrix
    do iao=1, nbasis    
       do jao=1, nbasis 
          aodum = ao_mat(iao,jao)
          do kmo=1, nrorb
             tmp(kmo,jao) = tmp(kmo,jao) + rotmatAA(iao,kmo)*aodum
          end do
       end do
    end do
    
    do jao=1, nbasis
       do lmo=1, nrorb
          do kmo=1, nrorb
             mo_mat(kmo,lmo) = mo_mat(kmo,lmo) + tmp(kmo,jao)*rotmatBB(jao,lmo)
          end do
       end do
    end do

    
  end subroutine ao2mo_complex
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FORM1H
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine form1det(noa,nva,nstates,i2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CIS states          !
    ! <C> Xmat(I,A,0,0) = X(I,A)Sqrt(2)  FOR UNRESTRICTED NO SQRT(2)!
    ! <C> Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)         !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
    
    real(8), intent(inout) :: i2a_mat(nstates*nstates)

    integer(8), optional, intent(in) :: nob,nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)


    real(8) :: rdum, sqrt2
    integer(8) :: ii, jj, aa, bb, i, j, a, b
    integer(8) :: k, ia2, ia, jb, kk
    

    sqrt2  = dsqrt(2.d0)
    if( present(nob) ) sqrt2 = 1.d0    
    

    i2a_mat    = 0.d0
    i2a_mat(1) = mo00    

    ! <C> < ia | h(1) | 0 > == 1/sqrt(2) { <ia|h(1)|0> + <IA|h(1)|0> } = 2/sqrt(2) <ia|h(1)|0>
    ! <C> < ia | h(1) | 0 > == < a | h(1) | i >
    ia0 : do ia = 2, nstates

       ia2 = (ia-1)*nstates + 1
       ii  = hole(ia)  ;  i = abs(ii)
       aa  = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       
       if ( ii.lt.0 .and. aa.lt.0 ) then
          i2a_mat(ia)  = sqrt2 * mo_matAA(i,a)
          i2a_mat(ia2) = sqrt2 * mo_matAA(i,a)
       end if
       if ( ii.gt.0 .and. aa.gt.0 ) then
          i2a_mat(ia)  = sqrt2 * mo_matBB(i,a)
          i2a_mat(ia2) = sqrt2 * mo_matBB(i,a)
       end if

    end do ia0

    
    ia_jb : do ia = 2, nstates

       ii = hole(ia)  ;  i = abs(ii)
       aa = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       

       do jb = 2 , nstates

          jj = hole(jb)  ;  j = abs(jj)
          bb = part(jb)  ;  b = abs(bb) + noa
          if ( bb.gt.0 ) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          rdum = 0.d0

          if ( ii.eq.jj ) then
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = mo_matAA(a,b)
             if ( bb.gt.0 .and. bb.gt.0 ) rdum = mo_matBB(a,b)
          end if

          if ( aa.eq.bb ) then
             if ( ii.lt.0 .and. jj.lt.0 ) rdum = rdum - mo_matAA(i,j)
             if ( ii.gt.0 .and. jj.gt.0 ) rdum = rdum - mo_matBB(i,j)
          end if

          if ( ia.eq.jb ) rdum = rdum + mo00
          i2a_mat(k) = rdum
          
       end do
    end do ia_jb



  end subroutine form1det
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ZFORM1H
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine Zform1det(noa,nva,nstates,Zi2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CIS states          !
    ! <C> Xmat(I,A,0,0) = X(I,A)Sqrt(2)  FOR UNRESTRICTED NO SQRT(2)!
    ! <C> Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)         !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00  

    integer(8), intent(in) :: nob,nvb
    real(8),    intent(in) :: mo_matBB(nob+nvb,nob+nvb)

    complex(8), intent(inout) :: Zi2a_mat(nstates*nstates)
    

    complex(8) :: cdum
    integer(8) :: ii, jj, aa, bb, i, j, a, b
    integer(8) :: k, ia2, ia, jb, kk

    
    Zi2a_mat    = dcmplx( 0.d0,0.d0 )
    Zi2a_mat(1) = dcmplx( mo00, 0.d0 )


    !: < ia | h(1) | 0 > == < a | h(1) | i >
    ia0 : do ia = 2, nstates

       ia2 = (ia-1)*nstates + 1
       ii  = hole(ia)  ;  i = abs(ii)
       aa  = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       
        if ( ii.lt.0 .and. aa.lt.0 ) then
           Zi2a_mat(ia)  = dcmplx( mo_matAA(i,a), 0.d0 )
           Zi2a_mat(ia2) = dcmplx( mo_matAA(i,a), 0.d0 )
        end if
        if ( ii.gt.0 .and. aa.gt.0 ) then
           Zi2a_mat(ia)  = dcmplx( mo_matBB(i,a), 0.d0 )
           Zi2a_mat(ia2) = dcmplx( mo_matBB(i,a), 0.d0 )
        end if
       
    end do ia0

    
    ia_jb : do ia = 2, nstates

       ii = hole(ia)  ;  i = abs(ii)
       aa = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob

       do jb = 2 , nstates
          
          jj = hole(jb)  ;  j = abs(jj)
          bb = part(jb)  ;  b = abs(bb) + noa
          if ( bb.gt.0 ) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          cdum = dcmplx( 0.d0,0.d0 )
          
          if ( ii.eq.jj ) then
              if ( aa.lt.0 .and. bb.lt.0 ) cdum = dcmplx( mo_matAA(a,b),0.d0 )
              if ( aa.gt.0 .and. bb.gt.0 ) cdum = dcmplx( mo_matBB(a,b),0.d0 )
          end if
          
          if ( aa.eq.bb ) then
              if ( ii.lt.0 .and. jj.lt.0 ) cdum = cdum - dcmplx( mo_matAA(i,j),0.d0 )
              if ( ii.gt.0 .and. jj.gt.0 ) cdum = cdum - dcmplx( mo_matBB(i,j),0.d0 )
          end if
          
          if ( ia.eq.jb ) cdum = cdum + dcmplx( mo00,0.d0 )
          zi2a_mat(k) = cdum

       end do
    end do ia_jb
    


  end subroutine Zform1det
  !==================================================================!
  !==================================================================!
  subroutine form1det_ip(noa,nva,nstates,i2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CISD-IP states      !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
    
    real(8), intent(inout) :: i2a_mat(nstates*nstates)

    integer(8), optional, intent(in) :: nob,nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)


    real(8) :: rdum
    integer(8) :: ii, jj, aa, bb, xx, yy, i, j, a, b, x, y
    integer(8) :: k, ia, jb


    i2a_mat = 0.d0
    
    do ia=1, nstates
       
       xx = hole(ia,1)   ;  x = abs(xx)
       ii = hole(ia,2)   ;  i = abs(ii)
       aa = part(ia)     ;  a = abs(aa) + noa 
       if (aa.gt.0) a = abs(aa) + nob


       do jb=1, nstates
          
          yy = hole(jb,1)  ;  y = abs(yy)
          jj = hole(jb,2)  ;  j = abs(jj)
          bb = part(jb)    ;  b = abs(bb) + noa 
          if (bb.gt.0) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          rdum = 0.d0


          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
              !: < x |A| y >
              if ( xx.lt.0 .and. yy.lt.0 ) rdum = - mo_matAA(y,x)
              !: < X |A| Y >
              if ( xx.gt.0 .and. yy.gt.0 ) rdum = - mo_matBB(Y,X)
             go to 78
          end if SS


          !: singles doubles
          
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 ) then             
             !: < i->a , x |A| x >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) rdum = mo_matAA(a,i) 
             !: < I->A, x |A| x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) rdum = mo_matBB(A,I) 
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 ) then
             !: < x |A| j->b, x >  =  < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) rdum = mo_matAA(b,j)
             !: < x |A| J->B, x >  =  < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) rdum = mo_matBB(b,j)
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
              !: < y->a, x |A| y >  =  - < a |A| x >
              if ( xx.lt.0 .and. aa.lt.0 ) rdum = - mo_matAA(a,x)
              !: < Y->A, X |A| Y >  =  - < A |A| X >
              if ( xx.gt.0 .and. aa.gt.0 ) rdum = - mo_matBB(A,X)
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
              !: < x |A| x->b, y >  =  - < y | A | b >
              if ( yy.lt.0 .and. bb.lt.0 ) rdum = - mo_matAA(b,y)
              !: < X |A| X->B, Y >  =  - < Y |A| B >
              if ( yy.gt.0 .and. bb.gt.0 ) rdum = - mo_matBB(B,Y)
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) rdum = - mo_matAA(i,j)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( ii.gt.0 .and. jj.gt.0 ) rdum = - mo_matBB(I,J)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum + mo_matAA(a,b)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum + mo_matBB(A,B)
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             !: < i->a, x |A| i->a, y >  =  - < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) rdum = rdum - mo_matAA(x,y)
             !: < I->A, X |A| I->A, Y >  =  - < Y |A| X >
             if ( xx.gt.0 .and. yy.gt.0 ) rdum = rdum - mo_matBB(X,Y)
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             !: < i->a, x |A| x->a, y >  =  < y |A| i >
             if ( ii.lt.0 .and. yy.lt.0 ) rdum = rdum + mo_matAA(i,y)
             !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
             if ( ii.gt.0 .and. yy.gt.0 ) rdum = rdum + mo_matBB(I,Y) 
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             !: < y->a, x |A| j->a, y >  =  < j |A| x >
             if ( jj.lt.0 .and. xx.lt.0 ) rdum = rdum + mo_matAA(j,x)
             !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
             if ( jj.gt.0 .and. xx.gt.0 ) rdum = rdum + mo_matBB(J,X)
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             !: < y->a, x |A| x->b, y >  =  - < a |A| b >
             if ( aa.lt.0 .and. bb.lt.0 ) rdum = rdum - mo_matAA(a,b)
             !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
             if ( aa.gt.0 .and. bb.gt.0 ) rdum = rdum - mo_matBB(A,B)
          end if kdelta_iy_xj
          
78        continue
          
          if( ia.eq.jb ) rdum = rdum + mo00
          
          i2a_mat(k) = rdum
          
       end do
    end do

    
  end subroutine form1det_ip
  !==================================================================!
  !==================================================================!
  subroutine Zform1det_ip(noa,nva,nstates,Zi2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CISD-IP SOC states  !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
    
    complex(8), intent(inout) :: Zi2a_mat(nstates*nstates)

    integer(8), optional, intent(in) :: nob,nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)


    complex(8) :: cdum
    integer(8) :: ii, jj, aa, bb, xx, yy, i, j, a, b, x, y
    integer(8) :: k, ia, jb


    Zi2a_mat    = dcmplx( 0.d0,0.d0 )
    
    do ia=1, nstates
       
       xx = hole(ia,1)   ;  x = abs(xx)
       ii = hole(ia,2)   ;  i = abs(ii)
       aa = part(ia)     ;  a = abs(aa) + noa 
       if (aa.gt.0) a = abs(aa) + nob


       do jb=1, nstates
          
          yy = hole(jb,1)  ;  y = abs(yy)
          jj = hole(jb,2)  ;  j = abs(jj)
          bb = part(jb)    ;  b = abs(bb) + noa 
          if (bb.gt.0) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          cdum = dcmplx( 0.d0,0.d0 )


          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
             !: < x |A| y >
             if ( xx.lt.0 .and. yy.lt.0 ) cdum = - dcmplx( mo_matAA(y,x), 0.d0)
             !: < X |A| Y >
             if ( xx.gt.0 .and. yy.gt.0 ) cdum = - dcmplx( mo_matBB(Y,X), 0.d0)
             go to 78
          end if SS


          !: singles doubles
          
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 ) then             
             !: < i->a , x |A| x >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) cdum = dcmplx( mo_matAA(a,i), 0.d0)
             !: < I->A, x |A| x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) cdum = dcmplx( mo_matBB(A,I), 0.d0)
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 ) then
             !: < x |A| j->b, x >  =  < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) cdum = dcmplx( mo_matAA(j,b), 0.d0)
             !: < x |A| J->B, x >  =  < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) cdum = dcmplx( mo_matBB(J,B), 0.d0)
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
             !: < y->a, x |A| y >  =  - < a |A| x >
             if ( xx.lt.0 .and. aa.lt.0 ) cdum = - dcmplx( mo_matAA(a,x), 0.d0)
             !: < Y->A, X |A| Y >  =  - < A |A| X >
             if ( xx.gt.0 .and. aa.gt.0 ) cdum = - dcmplx( mo_matBB(A,X), 0.d0)
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
             !: < x |A| x->b, y >  =  - < y | A | b >
             if ( yy.lt.0 .and. bb.lt.0 ) cdum = - dcmplx( mo_matAA(y,b), 0.d0)
             !: < X |A| X->B, Y >  =  - < Y |A| B >
             if ( yy.gt.0 .and. bb.gt.0 ) cdum = - dcmplx( mo_matBB(Y,B), 0.d0)
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) cdum = - dcmplx( mo_matAA(j,i), 0.d0)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( II.gt.0 .and. JJ.gt.0 ) cdum = - dcmplx( mo_matBB(J,I), 0.d0)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) cdum = cdum + dcmplx( mo_matAA(a,b), 0.d0)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( AA.gt.0 .and. BB.gt.0 ) cdum = cdum + dcmplx( mo_matBB(A,B), 0.d0)
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             !: < i->a, x |A| i->a, y >  =  - < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) cdum = cdum - dcmplx( mo_matAA(y,x), 0.d0)
             !: < I->A, X |A| I->A, Y >  =  - < Y |A| X >
             if ( XX.gt.0 .and. YY.gt.0 ) cdum = cdum - dcmplx( mo_matBB(Y,X), 0.d0)
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             !: < i->a, x |A| x->a, y >  =  < y |A| i >
             if ( ii.lt.0 .and. yy.lt.0 ) cdum = cdum + dcmplx( mo_matAA(y,i), 0.d0)
             !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
             if ( II.gt.0 .and. YY.gt.0 ) cdum = cdum + dcmplx( mo_matBB(Y,I), 0.d0)
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             !: < y->a, x |A| j->a, y >  =  < j |A| x >
             if ( jj.lt.0 .and. xx.lt.0 ) cdum = cdum + dcmplx( mo_matAA(j,x), 0.d0)
             !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
             if ( JJ.gt.0 .and. XX.gt.0 ) cdum = cdum + dcmplx( mo_matBB(J,X), 0.d0)
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             !: < y->a, x |A| x->b, y >  =  - < a |A| b >
             if ( aa.lt.0 .and. bb.lt.0 ) cdum = cdum - dcmplx( mo_matAA(a,b), 0.d0)
             !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
             if ( AA.gt.0 .and. BB.gt.0 ) cdum = cdum - dcmplx( mo_matBB(A,B), 0.d0)
          end if kdelta_iy_xj
          
78        continue
          
          if( ia.eq.jb ) cdum = cdum + dcmplx( mo00, 0.d0)
          
          Zi2a_mat(k) = cdum
          
       end do
    end do
    call testherm(nstates,Zi2a_mat,'TD_IP')

    
  end subroutine Zform1det_ip
  !==================================================================!
  !==================================================================!
  subroutine Zform1det_soc(noa,nva,nstates,Zi2a_mat,Zmo_matAA,hole,part,nob,nvb,Zmo_matBB,Zmo_matAB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CIS SOC states      !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nob, nvb, nstates
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    complex(8), intent(in) :: Zmo_matAA(noa+nva,noa+nva)
    complex(8), intent(in) :: Zmo_matBB(nob+nvb,nob+nvb)
    complex(8), intent(in) :: Zmo_matAB(noa+nva,nob+nvb)

    complex(8), intent(inout) :: Zi2a_mat(nstates*nstates)

    complex(8) :: cdum
    integer(8) :: ii, jj, aa, bb, i, j, k, a, b, ia, jb

    Zi2a_mat    = dcmplx( 0.d0,0.d0 )
    
    do ia=1, nstates
       
       ii = hole(ia,1)   ;  i = abs(ii)
       aa = part(ia)     ;  a = abs(aa) + noa 
       if (aa.gt.0) a = abs(aa) + nob

       do jb=1, nstates
          
          jj = hole(jb,1)  ;  j = abs(jj)
          bb = part(jb)    ;  b = abs(bb) + noa 
          if (bb.gt.0) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          cdum = dcmplx( 0.d0,0.d0 )

          !: ground state | singles 
          
          GS1: if ( jj.eq.0 .and. aa.ne.0 ) then             
             !:write(42,124) "GS1    ",ii,aa,jj,bb
             !: < i->a |A| 0 >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) cdum = Zmo_matAA(a,i)
             !: < I->A |A| 0 >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) cdum = Zmo_matBB(A,I)
             !: < i->A |A| 0 >  =  < A |A| I > 
             if ( ii.lt.0 .and. aa.gt.0 ) cdum = Dconjg(Zmo_matAB(i,A))
             !: < I->a |A| 0 >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.lt.0 ) cdum = Zmo_matAB(a,I)
             go to 78
          end if GS1
          
          GS2: if ( ii.eq.0 .and. bb.ne.0 ) then
             !:write(42,124) "GS2    ",ii,aa,jj,bb
             !: < 0 |A| j->b > = < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) cdum = Zmo_matAA(j,b)
             !: < 0 |A| J->B > = < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) cdum = Zmo_matBB(J,B)
             !: < 0 |A| j->B > = < j |A| B >
             if ( jj.lt.0 .and. bb.gt.0 ) cdum = Zmo_matAB(j,B)
             !: < 0 |A| J->b > = < J |A| B >
             if ( jj.gt.0 .and. bb.lt.0 ) cdum = Dconjg(Zmo_matAB(b,J))
             go to 78
          end if GS2

          !: singles | singles 

          SS1 : if ( aa.eq.bb ) then
             !:write(42,124) "SS1",ii,aa,jj,bb
             !: < i->a |A| j->a >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) cdum = - Zmo_matAA(j,i)
             !: < I->A |A| J->A >  =  - < J |A| I > 
             if ( ii.gt.0 .and. jj.gt.0 ) cdum = - Zmo_matBB(J,I)
             !: < i->A |A| J->A >  =  - < J |A| i > 
             if ( ii.lt.0 .and. jj.gt.0 ) cdum = - Dconjg(Zmo_matAB(i,J))
             !: < I->A |A| j->A >  =  - < j |A| I > 
             if ( ii.gt.0 .and. jj.lt.0 ) cdum = - Zmo_matAB(j,I)
          end if SS1
          
          SS2 : if ( ii.eq.jj  ) then
             !:write(42,124) "SS2",ii,aa,jj,bb
             !: < i->a |A| i->b >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) cdum = cdum + Zmo_matAA(a,b)
             !: < I->A |A| I->B >  =  < A |A| B > 
             if ( aa.gt.0 .and. bb.gt.0 ) cdum = cdum + Zmo_matBB(A,B)
             !: < I->a |A| I->B >  =  < a |A| B > 
             if ( aa.lt.0 .and. bb.gt.0 ) cdum = cdum + Zmo_matAB(a,B)
             !: < I->A |A| I->b >  =  < A |A| b > 
             if ( aa.gt.0 .and. bb.lt.0 ) cdum = cdum + Dconjg(Zmo_matAB(b,A))
          end if SS2
          
78        continue

          k  = (ia-1)*nstates + jb
!:          write(42,"(6i5,6f16.10)") ii,aa,jj,bb,ia,jb,k,cdum
          Zi2a_mat(k) = cdum
          
       end do
    end do

!: confirm hermitian
    call testherm(nstates,Zi2a_mat,'SOC')
  124    format(A12,6i3)        

  end subroutine Zform1det_soc
  !==================================================================!
  !==================================================================!
  subroutine Zform1det_socip(noa,nva,nstates,Zi2a_mat,Zmo_matAA,hole,part,nob,nvb,Zmo_matBB,Zmo_matAB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CISD-IP SOC states  !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nob, nvb, nstates
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    complex(8), intent(in) :: Zmo_matAA(noa+nva,noa+nva)
    complex(8), intent(in) :: Zmo_matBB(nob+nvb,nob+nvb)
    complex(8), intent(in) :: Zmo_matAB(noa+nva,nob+nvb)

    complex(8), intent(inout) :: Zi2a_mat(nstates*nstates)

    complex(8) :: cdum
    integer(8) :: ii, jj, aa, bb, xx, yy, i, j, a, b, x, y
    integer(8) :: k,k1, ia, jb


    Zi2a_mat    = dcmplx( 0.d0,0.d0 )
    
    do ia=1, nstates
       
       xx = hole(ia,1)   ;  x = abs(xx)
       ii = hole(ia,2)   ;  i = abs(ii)
       aa = part(ia)     ;  a = abs(aa) + noa 
       if (aa.gt.0) a = abs(aa) + nob


       do jb=1, nstates
          
          yy = hole(jb,1)  ;  y = abs(yy)
          jj = hole(jb,2)  ;  j = abs(jj)
          bb = part(jb)    ;  b = abs(bb) + noa 
          if (bb.gt.0) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          cdum = dcmplx( 0.d0,0.d0 )

          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
             !:write(42,124) "SS     ",ii,aa,xx,jj,bb,yy
             !: < 0,x |A| 0,y > = < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) cdum = - Zmo_matAA(y,x)
             !: < 0,X |A| 0,Y > = < Y |A| X >
             if ( xx.gt.0 .and. yy.gt.0 ) cdum = - Zmo_matBB(Y,X)
             !: < 0,x |A| 0,Y > = < Y |A| x >
             if ( xx.lt.0 .and. yy.gt.0 ) cdum = - Dconjg(Zmo_matAB(x,Y))
             !: < 0,X |A| 0,y > = < y |A| X >
             if ( xx.gt.0 .and. yy.lt.0 ) cdum = - Zmo_matAB(y,X)
             go to 78
          end if SS

          !: singles doubles
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 .and. aa.ne.0 ) then             
             !:write(42,124) "SD1    ",ii,aa,xx,jj,bb,yy
             !: < i->a, x |A| 0,x >  =  < a |A| i >
             if ( ii.lt.0 .and. aa.lt.0 ) cdum = Zmo_matAA(a,i)
             !: < I->A, x |A| 0,x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.gt.0 ) cdum = Zmo_matBB(A,I)
             !: < i->A, x |A| 0,x >  =  < A |A| I > 
             if ( ii.lt.0 .and. aa.gt.0 ) cdum = Dconjg(Zmo_matAB(i,A))
             !: < I->a, x |A| 0,x >  =  < A |A| I > 
             if ( ii.gt.0 .and. aa.lt.0 ) cdum = Zmo_matAB(a,I)
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 .and. bb.ne.0 ) then
             !:write(42,124) "SD2    ",ii,aa,xx,jj,bb,yy
             !: < 0,x |A| j->b, x > = < 0,X |A| j->b, X > = < j |A| b >
             if ( jj.lt.0 .and. bb.lt.0 ) cdum = Zmo_matAA(j,b)
             !: < 0,x |A| J->B, x > = < 0,X |A| J->B, X > = < J |A| B >
             if ( jj.gt.0 .and. bb.gt.0 ) cdum = Zmo_matBB(J,B)
             !: < 0,x |A| j->B, x > = < 0,X |A| j->B, X > = < j |A| B >
             if ( jj.lt.0 .and. bb.gt.0 ) cdum = Zmo_matAB(j,B)
             !: < 0,x |A| J->b, x > = < 0,X |A| J->b, X > = < J |A| B >
             if ( jj.gt.0 .and. bb.lt.0 ) cdum = Dconjg(Zmo_matAB(b,J))
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 .and. aa.ne.0 ) then
             !:write(42,124) "SD3    ",ii,aa,xx,jj,bb,yy
             !: < y->a, x |A| 0,y >  =  - < a |A| x >
             if ( xx.lt.0 .and. aa.lt.0 ) cdum = - Zmo_matAA(a,x)
             !: < Y->A, X |A| 0,Y >  =  - < A |A| X >
             if ( xx.gt.0 .and. aa.gt.0 ) cdum = - Zmo_matBB(A,X)
             !: < y->a, X |A| 0,y >  =  - < a |A| X >
             if ( xx.gt.0 .and. aa.lt.0 ) cdum = - Zmo_matAB(a,X)
             !: < Y->A, x |A| 0,Y >  =  - < A |A| x >
             if ( xx.lt.0 .and. aa.gt.0 ) cdum = - Dconjg(Zmo_matAB(x,A))
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 .and. bb.ne.0 ) then
             !:write(42,124) "SD4    ",ii,aa,xx,jj,bb,yy
             !: < 0,x |A| x->b, y >  =  - < y | A | b >
             if ( yy.lt.0 .and. bb.lt.0 ) cdum = - Zmo_matAA(y,b)
             !: < 0,X |A| X->B, Y >  =  - < Y |A| B >
             if ( yy.gt.0 .and. bb.gt.0 ) cdum = - Zmo_matBB(Y,B)
             !: < 0,x |A| x->b, Y >  =  - < Y | A | b >
             if ( yy.gt.0 .and. bb.lt.0 ) cdum = - Dconjg(Zmo_matAB(b,Y))
             !: < 0,X |A| X->B, y >  =  - < y | A | B >
             if ( yy.lt.0 .and. bb.gt.0 ) cdum = - Zmo_matAB(y,B)
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             !:write(42,124) "delta_xy_ab",ii,aa,xx,jj,bb,yy
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 .and. jj.lt.0 ) cdum = - Zmo_matAA(j,i)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( ii.gt.0 .and. jj.gt.0 ) cdum = - Zmo_matBB(J,I)
             !: < i->A, x |A| J->A, x >  =  - < J |A| i > 
             if ( ii.lt.0 .and. jj.gt.0 ) cdum = - Dconjg(Zmo_matAB(i,J))
             !: < I->A, x |A| j->A, x >  =  - < j |A| I > 
             if ( ii.gt.0 .and. jj.lt.0 ) cdum = - Zmo_matAB(j,I)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj .and. aa.ne.0 .and. bb.ne.0  ) then
             !:write(42,124) "delta_xy_ij",ii,aa,xx,jj,bb,yy
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 .and. bb.lt.0 ) cdum = cdum + Zmo_matAA(a,b)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( aa.gt.0 .and. bb.gt.0 ) cdum = cdum + Zmo_matBB(A,B)
             !: < I->a, x |A| I->B, x >  =  < a |A| B > 
             if ( aa.lt.0 .and. bb.gt.0 ) cdum = cdum + Zmo_matAB(a,B)
             !: < I->A, x |A| I->b, x >  =  < A |A| b > 
             if ( aa.gt.0 .and. bb.lt.0 ) cdum = cdum + Dconjg(Zmo_matAB(b,A))
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             !:write(42,124) "delta_ij_ab",ii,aa,xx,jj,bb,yy
             !: < i->a, x |A| i->a, y > = < I->A, x |A| I->A, y >  =  - < y |A| x >
             if ( xx.lt.0 .and. yy.lt.0 ) cdum = cdum - Zmo_matAA(y,x)
             !: < i->a, X |A| i->a, Y > = < I->A, X |A| I->A, Y >  =  - < Y |A| X >
             if ( xx.gt.0 .and. yy.gt.0 ) cdum = cdum - Zmo_matBB(Y,X)
             !: < i->a, x |A| i->a, Y > = < I->A, x |A| I->A, Y >  =  - < Y |A| x >
             if ( xx.lt.0 .and. yy.gt.0 ) cdum = cdum - Dconjg(Zmo_matAB(x,Y))
             !: < i->a, X |A| i->a, y > = < I->A, X |A| I->A, y >  =  - < y |A| X >
             if ( xx.gt.0 .and. yy.lt.0 ) cdum = cdum - Zmo_matAB(y,X)
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             !:write(42,124) "delta_xj_ab",ii,aa,xx,jj,bb,yy
             !: < i->a, x |A| x->a, y >  =  < y |A| i >
             if ( ii.lt.0 .and. yy.lt.0 ) cdum = cdum + Zmo_matAA(y,i)
             !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
             if ( ii.gt.0 .and. yy.gt.0 ) cdum = cdum + Zmo_matBB(Y,I)
             !: < i->a, x |A| x->a, Y >  =  < Y |A| i >
             if ( ii.lt.0 .and. yy.gt.0 ) cdum = cdum + Dconjg(Zmo_matAB(i,Y))
             !: < I->A, X |A| X->A, y >  =  < y |A| I >
             if ( ii.gt.0 .and. yy.lt.0 ) cdum = cdum + Zmo_matAB(y,I)
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             !:write(42,124) "delta_iy_ab",ii,aa,xx,jj,bb,yy
             !: < y->a, x |A| j->a, y >  =  < j |A| x >
             if ( xx.lt.0 .and. jj.lt.0 ) cdum = cdum + Zmo_matAA(j,x)
             !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
             if ( xx.gt.0 .and. jj.gt.0 ) cdum = cdum + Zmo_matBB(J,X)
             !: < y->a, X |A| j->a, y >  =  < j |A| X >
             if ( xx.gt.0 .and. jj.lt.0 ) cdum = cdum + Zmo_matAB(j,X)
             !: < Y->A, x |A| J->A, Y >  =  < J |A| x >
             if ( xx.lt.0 .and. jj.gt.0 ) cdum = cdum + Dconjg(Zmo_matAB(x,J))
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj .and. aa.ne.0 .and. bb.ne.0 ) then
             !:write(42,124) "delta_iy_xj",ii,aa,xx,jj,bb,yy
             !: < y->a, x |A| x->b, y >  =  - < a |A| b >
             if ( aa.lt.0 .and. bb.lt.0 ) cdum = cdum - Zmo_matAA(a,b)
             !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
             if ( aa.gt.0 .and. bb.gt.0 ) cdum = cdum - Zmo_matBB(A,B)
             !: < y->a, X |A| X->B, y >  =  - < a |A| B >
             if ( aa.lt.0 .and. bb.gt.0 ) cdum = cdum - Zmo_matAB(a,B)
             !: < Y->A, x |A| x->b, Y >  =  - < A |A| b >
             if ( aa.gt.0 .and. bb.lt.0 ) cdum = cdum - Dconjg(Zmo_matAB(b,A))
          end if kdelta_iy_xj
          
78        continue

          k  = (ia-1)*nstates + jb
!:          write(42,123) ii,aa,xx,jj,bb,yy,ia,jb,k,cdum
          Zi2a_mat(k) = cdum
          
       end do
    end do

!: confirm hermitian
    call testherm(nstates,Zi2a_mat,'SOC_IP')
  123    format(9i5,6f20.12)
  124    format(A12,6i3)        

  end subroutine Zform1det_socip
  !==================================================================!
  !==================================================================!
  subroutine form1cis0(n,m,A,B,C)
    
    ! <C> ==========================================================!
    ! <C>  Transformation B^T A B, using work array C               !
    ! <C>  A(M*M), B(M*N), C(M*N), results written back to A(N*N)   !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in) :: m, n
    real(8), intent(inout) :: A(m*m)
    real(8), intent(in)    :: B(m*n), C(m*n)
    real(8)                :: zero, one

    zero = 0.d0
    one  = 1.d0

    call dgemm('n','n,',m,n,m,one,A,m,B,m,zero,C,m) 
    call dgemm('t','n,',n,n,m,one,B,m,C,m,zero,A,n) 

  end subroutine form1cis0
  !==================================================================!
  !==================================================================!
  subroutine Zform1cis0(n,m,A,B,C)
    
    ! <C> ==========================================================!
    ! <C>  Complex transformation B^T A B, using work array C       !
    ! <C>  A(M*M), B(M*N), C(M*N), results written back to A(N*N)   !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in)    :: m, n
    complex(8), intent(inout) :: A(m*m)
    complex(8), intent(in)    :: B(m*n), C(m*n)
    complex(8)                :: zero, one

    zero = dcmplx(0.d0,0.d0)
    one  = dcmplx(1.d0,0.d0)

    call zgemm('n','n,',m,n,m,one,A,m,B,m,zero,C,m) 
    call zgemm('c','n,',n,n,m,one,B,m,C,m,zero,A,n) 

  end subroutine Zform1cis0
  !==================================================================!
  !==================================================================!
  subroutine form1cis(nstates,nstuse,mat,cis_vec)
    
    ! <C> ==========================================================!
    ! <C>  Transformation B^T A B, with result written back to A    !
    ! <C>  A(N*N), B(N*M), V(N), results written back to A(M*M)     !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in) :: nstates, nstuse
    real(8), intent(in)    :: cis_vec(nstates*nstates)
    real(8), intent(inout) :: mat(nstates*nstates)
    
    integer(8) :: istate, jstate, ii, jj, iuse, kuse, kk, m, n, juse, jjj
    real(8)    :: dum
    real(8)    :: tmpmat(nstates)

    ! <C>  A*B and store in A
    do istate=1, nstates       

       ! <C> store a row from A to tmpmat
       do jstate=1,nstates
          jj = (jstate-1)*nstates
          tmpmat(jstate) = mat(jj+istate)
       end do
       
       ! <C> mat*cis_vec
       do kuse=1, nstuse
          kk = (kuse-1)*nstates
          mat(istate+kk) = dot_product( tmpmat(:), cis_vec(kk+1:kk+nstates) )
       end do
       
    end do

    
    ! <C>  B^T*(A*B) and store in A
    do kuse=1, nstuse       

       kk = (kuse-1)*nstates
       tmpmat(:) = mat(kk+1:kk+nstates)

       do iuse=1, nstuse
          ii = (iuse-1)*nstates
          mat(iuse+kk) = dot_product( tmpmat(:), cis_vec(ii+1:ii+nstates) )
       end do

    end do
    

    ! <C> Pack down if nstuse < nstates
    if(nstuse.lt.nstates) then
       do juse=2, nstuse
          jj  = (juse-1)*nstuse
          jjj = (juse-1)*nstates          
          do iuse=1, nstuse
             mat(iuse+jj)=mat(iuse+jjj)
          end do
       end do
    end if


  end subroutine form1cis
  !==================================================================!
  !==================================================================!
  subroutine Zform1cis(nstates,nstuse,Zmat,Zcis_vec)
    
    ! <C> ==========================================================!
    ! <C>  Transformation B^T A B, with result written back to A    !
    ! <C>  A(N*N), B(N*M), V(N), results written back to A(M*M)     !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in) :: nstates, nstuse
    complex(8), intent(in) :: Zcis_vec(nstates*nstates)
    complex(8), intent(inout) :: Zmat(nstates*nstates)
    
    integer(8) :: istate, jstate, ii, jj, iuse, kuse, kk, m, n, juse, jjj
    real(8)    :: rdum
    complex(8) :: cdum, Ztmpmat(nstates)
    
    
    ! <C>  A*B and store in A
    do istate=1, nstates       
       
       ! <C> store a row from A to tmpmat
       !: Zmat initially all real
       do jstate=1, nstates
          jj = (jstate-1)*nstates
          Ztmpmat(jstate) = Zmat(jj+istate)
       end do
       
       ! <C> mat*cis_vec
       do kuse=1, nstuse
          kk = (kuse-1)*nstates
          cdum = dcmplx( 0.d0,0.d0 )
          do jstate=1, nstates
             cdum = cdum + Ztmpmat(jstate)*Zcis_vec(kk+jstate)
          end do
          Zmat(kk+istate) = cdum
          !Zmat(istate+kk) = dot_product( Zcis_vec(kk+1:kk+nstates) , Ztmpmat(:) ) 
       end do
       
    end do

    
    ! <C>  B^T*(A*B) and store in A
    do kuse=1, nstuse       
       
       kk = (kuse-1)*nstates
       do jstate=1, nstates
          Ztmpmat(jstate) = Zmat(kk+jstate)
       end do
       
       !Ztmpmat(:) = Zmat(kk+1:kk+nstates)
       
       do iuse=1, nstuse
          ii = (iuse-1)*nstates
          cdum = dcmplx(0.d0,0.d0)
          do jstate=1, nstates
             cdum = cdum + dconjg( Zcis_vec(ii+jstate) ) * Ztmpmat(jstate)
          end do
          !Zmat(kk+iuse) = dot_product( Ztmpmat(:), Zcis_vec(ii+1:ii+nstates) ) 
          Zmat(kk+iuse) = cdum
       end do
       
    end do
    
    
    ! <C> Pack down if nstuse < nstates
    if(nstuse.lt.nstates) then
       do juse=2, nstuse
          jj  = (juse-1)*nstuse
          jjj = (juse-1)*nstates          
          do iuse=1, nstuse
             Zmat(iuse+jj) = Zmat(iuse+jjj)
          end do
       end do
    end if


  end subroutine Zform1cis
  !==================================================================!
  !==================================================================!
  subroutine eigen_mat(iout,n,a,lda,w,QeigenDC)

    !: call lapack dysev or dysevd to get eigenstates and eigenvalues    
    
    implicit none
    
    integer(8) :: iout, n, lda, lwork, liwork, info, iwork_tmp(2)
    real(8)    :: a(:), w(:), work_tmp(2)
    logical    :: QeigenDC
    integer(8), allocatable :: iwork(:)
    real(8),    allocatable :: work(:)
    
!:    If(QeigenDC) then

      ! Divide and Conquer Lapack eigensolver
      lwork = -1
      liwork = -1
      call dsyevd('vectors','upper',n,a,lda,w,work_tmp,lwork,iwork_tmp,liwork,info)
      lwork = int(work_tmp(1))
      liwork = int(iwork_tmp(1))

      info  = 10
      write(iout,"(' WORK and IWORK space allocated for dsyevd:  ',2i10)") lwork, liwork
      allocate( work(lwork) )
      allocate( iwork(liwork) )
      call dsyevd('vectors','upper',n,a,lda,w,work,lwork,iwork,liwork,info)
      call fix_phase(n,a)
      deallocate( iwork )
      deallocate( work )

!:    else
!:
!:      ! Standard Lapack eigensolver
!:      lwork = -1
!:      call dsyev('vectors','upper',n,a,lda,w,work_tmp,lwork,info)
!:      lwork = int(work_tmp(1))
!:      lwork = n*n
!:
!:      info  = 10
!:      write(iout,"(' WORK space allocated for dsyev:  ',i10)") lwork
!:      allocate( work(lwork) )
!:      call dsyev('vectors','upper',n,a,lda,w,work,lwork,info)
!:      call fix_phase(n,a)
!:      deallocate( work )
!:
!:    end if
    
    if( info.ne.0 ) write( iout,"('dsyev digonalization error')") info

  end subroutine eigen_mat
  !==================================================================!
  !==================================================================!
  subroutine Zeigen_mat(iout,n,Za,lda,w,QeigenDC)

    !: call lapack zheev or zheevd to get eigenstates and eigenvalues    
    
    implicit none
    
    integer(8) :: iout, n, lda, lwork, lrwork, liwork, info, iwork_tmp(2)
    real(8)    :: w(:), rwork_tmp(2)
    complex(8) :: Za(:), work_tmp(2)
    logical    :: QeigenDC
    integer(8), allocatable :: iwork(:)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)

!:    If(QeigenDC) then

      ! Divide and Conquer Lapack complex eigensolver
      lwork = -1
      lrwork = -1
      liwork = -1
      call zheevd('vectors','upper',n,Za,lda,w,work_tmp,lwork,rwork_tmp,lrwork,iwork_tmp,liwork,info)
      lwork =  int(work_tmp(1))
      lrwork = int(rwork_tmp(1))
      liwork = int(iwork_tmp(1))

      info  = 10
      write(iout,"(' WORK, RWORK and IWORK space allocated for zheevd:  ',3i10)") lwork, lrwork, liwork
      allocate( work(lwork) )
      allocate( rwork(lrwork) ) ;  rwork = 0.d0
      allocate( iwork(liwork) )
      call zheevd('vectors','upper',n,Za,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
      call Zfix_phase(n,Za)
      deallocate( work )
      deallocate( rwork )
      deallocate( iwork )

!:    else
!:
!:      ! standard Lapack complex eigensolver
!:      lwork = -1
!:      lrwork = -1
!:      call zheev('vectors','upper',n,Za,lda,w,work_tmp,lwork,rwork_tmp,info)
!:      lwork  = int(work_tmp(1))
!:      lrwork = 3*n-2
!:      lwork  = n*n
!:
!:      info  = 10
!:      write(iout,"(' WORK and RWORK space allocated for zheev:  ',2i10)") lwork, lrwork
!:      allocate( work(lwork) )
!:      allocate( rwork(lrwork) ) ;  rwork = 0.d0
!:      call zheev('vectors','upper',n,Za,lda,w,work,lwork,rwork,info)
!:      call Zfix_phase(n,Za)
!:      deallocate( work )
!:      deallocate( rwork )
!:
!:    end if 
      
    if( info.ne.0 ) write( iout,"('zheev digonalization error')") info 

  end subroutine Zeigen_mat
  !==================================================================!
  !==================================================================!
  subroutine exp_mat(iout,nstuse,mat,expmat,hdt,QeigenDC)

    ! <C> get exponential of Hermetian matrix
    ! <C> diag_mat = [R]^T [mat] [R]
    ! <C> Exp_mat  = [R] exp( [diag_mat] ) [R]^T 
    ! <C> Exp_mat  = [R] [exp(eig)] [R]^T
    
    implicit none
    
    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: mat(nstuse*nstuse), hdt
    real(8), intent(inout) :: expmat(nstuse*nstuse)
    logical,    intent(in) :: QeigenDC
    
    integer(8) :: iout,i,k,kk,j,jj
    real(8)    :: mat_k, exp_eig
    real(8)    :: eig(nstuse), zero, one
    real(8),  allocatable :: mat_cp(:)

    allocate( mat_cp(nstuse*nstuse) )
    mat_cp = mat

    ! <C> diagonalize matrix
    call eigen_mat(iout,nstuse,mat_cp,nstuse,eig,QeigenDC)
    
    ! <C> Form the exponential of the matrix    
!:    expmat = 0.d0    
!:    do j=1, nstuse
!:       jj   = (j-1)*nstuse
!:       exp_eig = dexp( hdt*eig(j) )
!:       do k=1, nstuse
!:          kk = (k-1)*nstuse
!:          mat_k = mat_cp(jj+k) * exp_eig
!:          expmat(kk+1:kk+nstuse) = expmat(kk+1:kk+nstuse) + mat_k * mat_cp(jj+1:jj+nstuse)
!:       end do
!:    end do
 
    zero = 0.d0
    one  =1.d0
    do j=1, nstuse
       jj   = (j-1)*nstuse
       exp_eig = dexp( 0.5d0*hdt*eig(j) )
       call dscal(nstuse,exp_eig,mat_cp(jj+1:jj+nstuse),1)
    end do
    call dgemm('n','t',nstuse,nstuse,nstuse,one,mat_cp,nstuse,mat_cp,nstuse,zero,expmat,nstuse)
 
    deallocate( mat_cp )

  end subroutine exp_mat
  !==================================================================!
  !==================================================================!
  subroutine Zexp_mat(iout,nstuse,Zmat,Zexpmat,hdt,QeigenDC)

    ! <C> get exponential of Hermetian matrix
    ! <C> diag_mat = [R]^T [mat] [R]
    ! <C> Exp_mat  = [R] exp( [diag_mat] ) [R]^T 
    ! <C> Exp_mat  = [R] [exp(eig)] [R]^T
    
    implicit none
    
    integer(8),    intent(in) :: nstuse
    real(8),       intent(in) :: hdt
    complex(8),    intent(in) :: Zmat(nstuse*nstuse)
    complex(8), intent(inout) :: Zexpmat(nstuse*nstuse)
    logical,       intent(in) :: QeigenDC
    
    integer(8) :: iout,i,k,kk,j,jj
    real(8)    :: exp_eig
    complex(8) :: mat_k, zero, one
    
    real(8)    :: eig(nstuse)
    complex(8), allocatable :: mat_cp(:)
    
    allocate( mat_cp(nstuse*nstuse) )
    mat_cp = Zmat

    ! <C> diagonalize matrix
    call Zeigen_mat(iout,nstuse,mat_cp,nstuse,eig,QeigenDC)
    
    ! <C> Form the exponential of the matrix    
!:    Zexpmat = dcmplx( 0.d0,0.d0 )
!:    do j=1, nstuse
!:       jj   = (j-1)*nstuse
!:       exp_eig = dexp( hdt*eig(j) )
!:       do k=1, nstuse
!:          kk = (k-1)*nstuse
!:          mat_k = mat_cp(jj+k) * exp_eig
!:          Zexpmat(kk+1:kk+nstuse) = Zexpmat(kk+1:kk+nstuse) + mat_k * dconjg(mat_cp(jj+1:jj+nstuse))
!:       end do
!:    end do
 
    zero = dcmplx(0.d0,0.d0)
    one  = dcmplx(1.d0,0.d0)
    do j=1, nstuse
       jj   = (j-1)*nstuse
       exp_eig = dexp( 0.5d0*hdt*eig(j) )
       do k=1, nstuse
          mat_cp(jj+k) = mat_cp(jj+k) * exp_eig
       end do
!:       call zscal(nstuse,dcmplx(exp_eig,0.d0),mat_cp(jj+1:jj+nstuse),1)
    end do
    call zgemm('n','c',nstuse,nstuse,nstuse,one,mat_cp,nstuse,mat_cp,nstuse,zero,Zexpmat,nstuse)
!:    Zexpmat = dconjg(Zexpmat) 
    deallocate( mat_cp )

  end subroutine Zexp_mat
  !==================================================================!
  !==================================================================!
  subroutine matmult(A,B,nsize)
   
    
    ! <C> commented section does
    ! <C> Matrix Multiplication A*B=C
    ! <C> C is written back to A       

    
    integer(8), intent(in) :: nsize
    real(8),    intent(in) :: B(nsize*nsize) !B(nsize,nsize)
    real(8), intent(inout) :: A(nsize*nsize) !A(nsize,nsize)

    integer(8) :: i,ii,j,jj,k,kk
    real(8) :: rdum, scratch(nsize)    !scratch(nsize,nsize)

    
    !scratch = matmul( A,B )
    !A = scratch

    do i=1, nsize
       ii=(i-1)*nsize
       do j=1, nsize
          scratch(j) = A(ii+j)
       end do
       do k=1, nsize
          kk = (k-1)*nsize
          rdum = 0.d0
          do j=1, nsize
             rdum = rdum+scratch(j)*B(j+kk)
          end do
          A(ii+k) = rdum
       end do
    end do


    do i=2, nsize
       ii=(i-1)*nsize
       do j=1, i
          jj = (j-1)*nsize
          rdum = A(i+jj)
          A(i+jj) = A(ii+j)
          A(ii+j) = rdum
       end do
    end do

    
  end subroutine matmult
  !==================================================================!
  !==================================================================!
  subroutine Zmatmult(A,B,nsize)
   
    
    ! <C> commented section does
    ! <C> Matrix Multiplication A*B=C
    ! <C> C is written back to A       

    
    integer(8), intent(in) :: nsize
    complex(8), intent(in) :: B(nsize*nsize) !B(nsize,nsize)
    complex(8), intent(inout) :: A(nsize*nsize) !A(nsize,nsize)

    integer(8) :: i,ii,j,jj,k,kk
    complex(8) :: cdum, scratch(nsize)    !scratch(nsize,nsize)

    
    !scratch = matmul( A,B )
    !A = scratch
    
    do i=1, nsize
       ii=(i-1)*nsize
       do j=1, nsize
          scratch(j) = A(ii+j)
       end do
       do k=1, nsize
          kk = (k-1)*nsize
          cdum = dcmplx( 0.d0,0.d0 )
          do j=1, nsize
!:             cdum = cdum+scratch(j)*B(j+kk)
             cdum = cdum+scratch(j)*dconjg(B(j+kk))
          end do
          A(ii+k) = cdum
       end do
    end do


    do i=2, nsize
       ii=(i-1)*nsize
       do j=1, i
          jj = (j-1)*nsize
          cdum = A(i+jj)
!:          A(i+jj) = A(ii+j)
!:          A(ii+j) = cdum
          A(i+jj) = dconjg(A(ii+j))
          A(ii+j) = dconjg(cdum)
       end do
    end do

    
  end subroutine Zmatmult
  !==================================================================!
  subroutine testherm(n,a,label)

    implicit none

    integer(8), intent(in)    :: n
    complex(8), intent(in)    :: A(n,n)
    character(*), intent(in)  :: label
    integer(4) :: i,j

 5  Format(' TestHerm: Matrix ',A10)
10  Format(A10,' not Hermitian: i,A(i,i) ',I6,2F15.10)
20  Format(A10,' not Hermitian: i,j,A(i,j) ',2I6,4F15.10)

    write(42,5) trim(label)
    do i=1,n
      if(abs(aimag(A(i,i))).gt.1.d-7) write(42,10) trim(label),i,a(i,i)
      do j=1,i
        if(abs(A(i,j)-DConjg(A(j,i))).gt.1.d-7) &
          write(42,20) trim(label),i,j,a(i,j),a(j,i)
      end do
    end do

  end subroutine testherm
  !==================================================================!
  subroutine math_ham

    use readintegrals

    implicit none
    integer(8) :: istate, xx, ii, aa, yy, jj, bb, ia, jb

     do ia=1, nstates

       xx = hole_index(ia,1) 
       ii = hole_index(ia,2) 
       aa = part_index(ia,1) 

       do jb=1, nstates

          yy = hole_index(jb,1) 
          jj = hole_index(jb,2) 
          bb = part_index(jb,1) 

          istate = (ia-1)*nstates + jb
          write(42,123) ii,aa,xx,jj,bb,yy,ia,jb,Zcis_vec(istate)

       end do
    end do
 123      format(8i5,3f20.12)

  end subroutine math_ham
  !=================================================================!
    subroutine get_proj_ion(iout,ip_states,ip_vec,Zion_coeff,Zproj_ion)

    implicit none
 
    integer(8), intent(in)    :: iout,ip_states
    real(8), intent(in)       :: ip_vec(ip_states*ip_states)
    complex(8), intent(in)    :: Zion_coeff(ip_states*ip_states)
    complex(8), intent(inout) :: Zproj_ion(ip_states*ip_states)

    integer(8) :: i,j,k
    
!:    write(iout,*) " enter get_proj_ion",nstates

    Zproj_ion = dcmplx(0.d0,0.d0)
    do k = 1,ip_states
      do i = 1,ip_states
        do j = 1,ip_states
          Zproj_ion(i+(k-1)*ip_states) = Zproj_ion(i+(k-1)*ip_states) + &
            ip_vec(j+(i-1)*ip_states)*Zion_coeff(j+(k-1)*ip_states)
        end do
      end do
    end do
 
!:    write(iout,*) " leave get_proj_ion"
 
    end subroutine get_proj_ion
  !=================================================================!
    subroutine get_Zproj_ion(iout,ip_states,Zip_vec,Zion_coeff,Zproj_ion)

    implicit none
 
    integer(8), intent(in)    :: iout,ip_states
    complex(8), intent(in)    :: Zion_coeff(ip_states*ip_states),Zip_vec(ip_states*ip_states)
    complex(8), intent(inout) :: Zproj_ion(ip_states*ip_states)

    integer(8) :: i,j,k
    
!:    write(iout,*) " enter get_Zproj_ion",nstates

    Zproj_ion = dcmplx(0.d0,0.d0)
    do k = 1,ip_states
      do i = 1,ip_states
        do j = 1,ip_states
          Zproj_ion(i+(k-1)*ip_states) = Zproj_ion(i+(k-1)*ip_states) + &
            dconjg(Zip_vec(j+(i-1)*ip_states))*Zion_coeff(j+(k-1)*ip_states)
        end do
      end do
    end do
 
!:    write(iout,*) " leave get_Zproj_ion"
 
    end subroutine get_Zproj_ion
  !=================================================================!
      subroutine get_polar(istates,nstates,tdx,tdy,tdz,cis_eig,polar)
 
      ! calculate polarizability lowest istates using sum over nstates states
 
      implicit none

      integer(8), intent(in)    ::  istates,nstates
      real(8),    intent(in)    ::  tdx(nstates*nstates),tdy(nstates*nstates),tdz(nstates*nstates)
      real(8),    intent(in)    ::  cis_eig(nstates) 
      real(8),    intent(inout) ::  polar(nstates,6)
 
      integer(8) i,ii,j
      real(8) deltaE
 
      polar = 0.d0
 
      do i = 1, istates
        ii = (i-1)*nstates
        do j = 1, nstates
          If( i.ne.j ) then
            deltaE = cis_eig(i) - cis_eig(j)
            If( abs(deltaE).gt.1.d-5 ) then
              polar(i,1) = polar(i,1) - & 
                2.d0*tdx(ii+j)*tdx(ii+j)/deltaE
              polar(i,2) = polar(i,2) - & 
                2.d0*tdy(ii+j)*tdy(ii+j)/deltaE
              polar(i,3) = polar(i,3) - & 
                2.d0*tdz(ii+j)*tdz(ii+j)/deltaE
              polar(i,4) = polar(i,4) - & 
                2.d0*tdx(ii+j)*tdy(ii+j)/deltaE
              polar(i,5) = polar(i,5) - & 
                2.d0*tdx(ii+j)*tdz(ii+j)/deltaE
              polar(i,6) = polar(i,6) - & 
                2.d0*tdy(ii+j)*tdz(ii+j)/deltaE
            end If
          end If
        end Do
      end Do
 
      end subroutine get_polar
  !=================================================================!
    subroutine get_Zpolar(istates,nstates,Ztdx,Ztdy,Ztdz,cis_eig,polar)
 
      ! calculate polarizability lowest istates using sum over nstates states
 
      implicit none

      integer(8), intent(in)    ::  istates,nstates
      complex(8), intent(in)    ::  Ztdx(nstates*nstates),Ztdy(nstates*nstates),Ztdz(nstates*nstates)
      real(8),    intent(in)    ::  cis_eig(nstates) 
      real(8),    intent(inout) ::  polar(nstates,6)
 
      integer(8) i,ii,j
      real(8) deltaE
 
      polar = 0.d0
 
      do i = 1, istates
        ii = (i-1)*nstates
        do j = 1, nstates
          If( i.ne.j ) then
            deltaE = cis_eig(i) - cis_eig(j)
            If( abs(deltaE).gt.1.d-5 ) then
              polar(i,1) = polar(i,1) - & 
                2.d0*real(dconjg(Ztdx(ii+j))*Ztdx(ii+j))/deltaE
              polar(i,2) = polar(i,2) - & 
                2.d0*real(dconjg(Ztdy(ii+j))*Ztdy(ii+j))/deltaE
              polar(i,3) = polar(i,3) - & 
                2.d0*real(dconjg(Ztdz(ii+j))*Ztdz(ii+j))/deltaE
              polar(i,4) = polar(i,4) - & 
                real(dconjg(Ztdx(ii+j))*Ztdy(ii+j) + &
                     dconjg(Ztdy(ii+j))*Ztdx(ii+j))/deltaE
              polar(i,5) = polar(i,5) - & 
                real(dconjg(Ztdx(ii+j))*Ztdz(ii+j) + &
                     dconjg(Ztdz(ii+j))*Ztdx(ii+j))/deltaE
              polar(i,6) = polar(i,6) - & 
                real(dconjg(Ztdy(ii+j))*Ztdz(ii+j) + &
                     dconjg(Ztdz(ii+j))*Ztdy(ii+j))/deltaE
            end If
          end If
        end Do
      end Do
 
    end subroutine get_Zpolar
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FIX_PHASE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine fix_phase(n,c)
 
    !: multiply eigenvectors so that the largest coefficient is real and positive  
    
    implicit none
    
    integer(8), intent(in) :: n
    real(8), intent(inout) :: c(n*n)
 
    integer(8) :: i, imax, j, jj
    real(8) :: cmax, cmax1
    real(8) :: phase,factor

    factor=1.d0
!:    write(42,*) " phase factor =",factor

    do j = 1, n
      jj = (j-1)*n
      cmax = 0.d0
      imax = 0
      do i = 1, n
        cmax1 = c(i+jj)*c(i+jj)
        if(cmax1 .gt. cmax+1.d-7) then
          imax = i
          cmax = cmax1
        end if
      end do
      if(imax.gt.0) then
        phase = factor*c(imax+jj)/dsqrt(cmax)
      else
        phase = factor
      end if
      do i = 1, n
        c(i+jj) = phase*c(i+jj)
      end do
    end do

  end subroutine fix_phase 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ZFIX_PHASE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine Zfix_phase(n,c)
 
    !: multiply eigenvectors so that the largest coefficient is real and positive  
    
    implicit none
    
    integer(8), intent(in) :: n
    complex(8), intent(inout) :: c(n*n)
 
    integer(8) :: i, imax, j, jj
    real(8) :: cmax, cmax1
    complex(8) :: phase,factor

    factor=dcmplx(1.d0,0.d0)
    factor=factor/dcmplx(dsqrt(dble(dconjg(factor)*factor))) 
!:    write(42,*) " phase factor =",factor

    do j = 1, n
      jj = (j-1)*n
      cmax = 0.d0
      imax = 0
      do i = 1, n
        cmax1 = dconjg(c(i+jj))*c(i+jj)
        if(cmax1 .gt. cmax+1.d-7) then
          imax = i
          cmax = cmax1
        end if
      end do
      if(imax.gt.0) then
        phase = factor*dconjg(c(imax+jj))/dsqrt(cmax)
      else
        phase = factor
      end if
      do i = 1, n
        c(i+jj) = phase*c(i+jj)
      end do
    end do

  end subroutine Zfix_phase 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FIX_ORDER
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine fix_order(na,nb,ntot,c1,c2)
 
    !: for a pair of complex eigenvectors with degenerate eigenvalues, 
    !: fix order so alpha spin comes before beta spin  
    !: first na coefficient are taken as alpha, last nb coefficients are taken as beta
    
    implicit none
    
    integer(8), intent(in) :: na, nb, ntot
    complex(8), intent(inout) :: c1(ntot), c2(ntot)
 
    integer(8) :: i
    real(8)    :: spin1, spin2
    complex(8) :: cdum

    spin1 = 0.d0
    spin2 = 0.d0
    do i = 1, na
      spin1 = spin1 + dble(dconjg(c1(i))*c1(i))
      spin2 = spin2 + dble(dconjg(c2(i))*c2(i))
      end do
    do i = ntot-nb+1, ntot
      spin1 = spin1 - dble(dconjg(c1(i))*c1(i))
      spin2 = spin2 - dble(dconjg(c2(i))*c2(i))
      end do

      If(spin2.gt.spin1) then
        do i = 1, ntot
          cdum  = c1(i)
          c1(i) = c2(i)
          c2(i) = cdum
          end do
        end if

  end subroutine fix_order 
end module util
    
