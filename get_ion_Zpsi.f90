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
       cdum = cdum + Zcis_vec(i+kk)*psi_det(i)
     end do
     psi(k) = psi(k) + cdum
!:     norm3 = norm3 + dble(dconjg(psi(k))*psi(k))
   end do
!:       write(iout,"('norm',8f16.12)") norm0,norm1,norm2,norm3,norm2-norm1,norm3-norm0,rate
!:     write(iout,"(9f12.7)") (psi(i),i=1,nstuse)

  end subroutine get_Zion_psi1
end module analyze_psi
