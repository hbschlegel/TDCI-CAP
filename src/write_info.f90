module write_info 


  use variables_global
  use util
  implicit none


  !: contains subroutine write_specifics1
  !: contains subroutine write_specifics2
  !: contains subroutine write_summary
  !: contains subroutine write_mo_energies
  !: contains subroutine write_field_shape
  !: contains subroutine write_ham0
  !: contains subroutine write_mo_elements
  !: contains subroutine save_restart_bin

  
contains
  !: --------------------------- :!
  !: subroutine write_specifics1 :!
  !: --------------------------- :!
  subroutine write_specifics1
    

    implicit none

    real(8)    :: fstrnth
    integer(8) :: i, ifield, idir, ij
    

    call write_header( 'write_specifics1','write_info','enter' )
    

    !: job title taken from TDCI.dat
    write(iout,'(A)') ' JOB TITLE: '//trim(job_title)

    !: molecule info
    write(iout,"(5x,'charge = ',i0,'  multiplicity = ',i0,'  natoms = ',i0)") ICharg, Multip, natoms
    write(iout,"(5x,'coordinates in Angstroms')")

    do iatom=1, natoms
       write(iout,"(5x,a4,2x,3(f17.11,2x))") &
            myatom(iatom),          &
            xcoord(iatom)*bohr2ang, &
            ycoord(iatom)*bohr2ang, &
            zcoord(iatom)*bohr2ang
    end do
    
    write(iout,'(A)') '     ground state dipole in au:'
    write(iout,"(5x,'xdip = ',f10.4,'  ydip = ',f10.4,'  zdip = ',f10.4)") dipx00, dipy00, dipz00
    write(iout,'(A)') '     ground state dipole in Debye:'
    write(iout,"(5x,'xdip = ',f10.4,'  ydip = ',f10.4,'  zdip = ',f10.4)") &
         dipx00*audip2debye, dipy00*audip2debye, dipz00*audip2debye

    !: jobtype info
    if( unrestricted )      write(iout,"(A)") ' unrestricted'
    if( .not.unrestricted ) write(iout,"(A)") ' restricted'

    select case( trim(jobtype) )
    case( flag_cis)   ; write(iout,'(A)') ' CIS'
    case( flag_tda)   ; write(iout,'(A)') ' TDA'
    case( flag_ip)    ; write(iout,'(A)') ' CISD-IP'
    case( flag_cisd ) ; write(iout,'(A)') ' CISD'
    end select
    
    !: system size info   
    write(iout,"(5x,'nbasis   = ',i0)") nbasis
    write(iout,"(5x,'noa      = ',i0)") noa
    write(iout,"(5x,'nva      = ',i0)") nva
    write(iout,"(5x,'nob      = ',i0)") nob
    write(iout,"(5x,'nvb      = ',i0)") nvb
    write(iout,"(5x,'nactive  = ',i0)") nactive
    write(iout,"(5x,'nvirtual = ',i0)") nvirtual
    write(iout,"(5x,'norb     = ',i0)") norb
    write(iout,"(5x,'nrorb    = ',i0)") nrorb
    write(iout,"(5x,'noanva   = ',i0)") noanva
    write(iout,"(5x,'nobnvb   = ',i0)") nobnvb
    write(iout,"(5x,'nstates  = ',i0)") nstates
    
    !: field info
    write( iout, '(A)' )  ' FIELD VARIABLES'
    do i=1, nemax
       fstrnth = tdciresults((i-1)*ndir+1)%fstrength0
       write(iout,100) i, fstrnth, fstrnth**2 * auflux2si
    end do
100 format( 5x,'max field strength',i3,' = ',f9.5,'  au',1x,'laser intensity ',es12.4,' W/cm2' )

    if ( envelope.eq.'cirl' .or. envelope.eq.'cirr' ) then
       write(iout,"(5x, 'nonlinear light')") 
       write(iout,"(5x, 'EULER angle euler    = ',f10.4)") euler
       write(iout,"(5x, 'ellipticity, E1/E2 = ',f10.4)")  ellipt
    end if
    
    if ( envelope.ne.'none' .and. envelope.ne.'stat' ) then
       write(iout,"(5x,'laser frequency       = ',f10.4,' au',es12.4,' Hz')") omega, omega*s2autime
       write(iout,"(5x,'laser wavelength      = ',f10.4,' nm')")              light *2.d0*pi / omega * autime2s * m2nm
       write(iout,"(5x,'pulse duration        = ',f10.4,' au',f12.4,' fs')")  field_duration * dble(dt), field_duration * dble(dt) *autime2s * s2fs

       select case ( trim(envelope) )
       case( 'none' ) ; write(iout,"(5x,'envelope none ',7x,' =   null pulse')")      
       case( 'cos2' ) ; write(iout,"(5x,'envelope cos2 ',7x,' =   cosine pulse')")    
       case( 'gaus' ) ; write(iout,"(5x,'envelope gaus ',7x,' =   gaussian pulse')")  
       case( 'stat' ) ; write(iout,"(5x,'envelope stat ',7x,' =   static pulse')")    
       case( 'band' ) ; write(iout,"(5x,'envelope band ',7x,' =   Bandrauk pulse')")  
       case( 'sin2' ) ; write(iout,"(5x,'envelope sin2 ',7x,' =   squared pulse')")   
       case( 'trap' ) ; write(iout,"(5x,'envelope trap ',7x,' =   trapezoidal pulse')")  
       case( 'cirl' ) ; write(iout,"(5x,'envelope cirl ',7x,' =   left circular sine squared pulse')") 
       case( 'cirr' ) ; write(iout,"(5x,'envelope cirr ',7x,' =   right circular sine squared pulse')")  
       end select

       write(iout,"(5x,'# of laser cycles     = ',i10,' cycles')")               ncyc
       write(iout,"(5x,'# timesteps/cycle     = ',i10)")                         int(period)
       write(iout,"(5x,'# timesteps for pulse = ',i10)")                         int(field_duration)
       write(iout,"(5x,'write results every     ',i10,  ' step   (total ',i0,' every ',f6.4,' fs)')") outstep, ndata, dble(outstep)*dt*autime2s* s2fs
       write(iout,"(5x,'timestep size         = ',f10.3,' au',f12.4,' fs')")     dt, dt * autime2s* s2fs
       write(iout,"(5x,'TOTAL SIMULATION TIME = ',f10.3,' au',f12.4,' fs (nstep=',i0,')')") &
            dt*dble(nstep), dt*dble(nstep)*autime2s*s2fs, nstep
    else
       write(iout,"(5x,'pulse duration        = ',f10.4,' au',f12.4,' fs')") field_duration * dble(dt), field_duration * dble(dt) *autime2s * s2fs
       write(iout,"(5x,'timestep size         = ',f10.3,' au',f12.4,' fs')")    dt, dt * autime2s* s2fs
       write(iout,"(5x,'TOTAL SIMULATION TIME = ',f10.3,' au',f12.4,' fs (nstep=',i0,')')") &
            dt*dble(nstep), dt*dble(nstep)*autime2s*s2fs, nstep
    end if

    !: field direction
    write( iout,'(A)' ) ' FIELD DIRECTIONS'
    if ( envelope.ne.'cirl' .and. envelope.ne.'cirr' ) then

       write(iout,50) '#' , 'E_0(au)' , 'x0' , 'y0' , 'z0' , 'theta0' , 'phi0'
       do i=1, nemax*ndir
          write(iout,60) i, tdciresults(i)%fstrength0, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%theta0, tdciresults(i)%phi0
       end do
       
    else
       write(iout,70) '#' , 'E_0(au)' , 'E_1(au)' , 'E_2(au)' ,           &
            'x0' , 'y0' , 'z0' , 'x1' , 'y1' , 'z1' , 'x2' , 'y2' , 'z2' , &
            'theta0' , 'phi0' , 'theta1' , 'phi1' , 'theta2' , 'phi2'
       do i=1, nemax*ndir
          write(iout,80) i, &
               tdciresults(i)%fstrength0, tdciresults(i)%fstrength1, tdciresults(i)%fstrength2, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%x1, tdciresults(i)%y1, tdciresults(i)%z1, &
               tdciresults(i)%x2, tdciresults(i)%y2, tdciresults(i)%z2, &
               tdciresults(i)%theta0, tdciresults(i)%phi0, &
               tdciresults(i)%theta1, tdciresults(i)%phi1, &
               tdciresults(i)%theta2, tdciresults(i)%phi2
       end do

    end if
       
50  format( 1x, a4, a8,   2x, ,a9  ,1x,a9,  1x,a9,  1x,2x,a9  ,1x,a9,  1x )
60  format( 1x, i4, f8.4, 2x, ,f9.4,1x,f9.4,1x,f9.4,1x,2x,f9.4,1x,f9.4,1x )
70  format( 1x,a4,a8,  2x,a8,  2x,a8,  2x,3(a9  ,1x,a9  ,1x,a9  ,1x,2x), 3(a9,  1x,a9,  1x,2x) )
80  format( 1x,i4,f8.4,2x,f8.4,2x,f8.4,2x,3(f9.4,1x,f9.4,1x,f9.4,1x,2x), 3(f9.4,1x,f9.4,1x,2x) )
       

    call write_header( 'write_specifics1','write_info','leave' )
    

  end subroutine write_specifics1
  !: ----------------------------- :!
  !: subroutine write_specific2    :!
  !: ----------------------------- :!
  subroutine write_specifics2


    implicit none
    integer(8) :: i, j, ii, ia, jb, istates
    real(8) rdum
    
    
    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip )  go to 78
    

    call write_header( 'write_specifics2','write_info','enter' )

    istates = 20
    write(iout,'(A)') ' ENERGIES (eV)'
    do i=1, istates
       write(iout,50) i, cis_eig(i)*au2eV
    end do

    write(iout,'(A)') ' NSTUSE' 
    write(iout,"(5x,'NstUse          = ',i12)") nstuse
    write(iout,"(5x,'Nstates         = ',i12)") nstates
    write(iout,"(5x,'cis_eig(nstuse) = ',f12.6,' au',f12.6,' eV')") cis_eig(nstuse), cis_eig(nstuse)*au2eV

    write(iout,'(A)') ' ENERGIES and EXPECTATION VALUES'
    write(iout,40) '#', '<i|H0|i> eV  ', '<i|mu_x|i> au', '<i|mu_y|i> au', '<i|mu_z|i> au', '<i|Vabs|i> au', &
     '<i|a_xx|i> au', '<i|a_yy|i> au', '<i|a_zz|i> au', '<i|a_xy|i> au', '<i|a_xz|i> au', '<i|a_yz|i> au'
    call get_polar(istates,nstates,tdx,tdy,tdz,cis_eig,polar)
    do i=1, istates
       ii = (i-1)*nstuse + i
       write(iout,50) i, cis_eig(i)*au2eV, tdx(ii), tdy(ii), tdz(ii), abp(ii), &
            (polar(i,j),j=1,6)
          do ia=1, nstates
             if ( trim(jobtype).eq.flag_cis .and. abs(cis_vec(ia+nstates*(i-1))).ge.0.05D0 ) &
               write(iout,61) hole_index(ia,1), part_index(ia,1), cis_vec(ia+nstates*(i-1))
             if ( trim(jobtype).eq.flag_ip .and. abs(cis_vec(ia+nstates*(i-1))).ge.0.05D0 ) &
               write(iout,60) hole_index(ia,1), hole_index(ia,2), part_index(ia,1), cis_vec(ia+nstates*(i-1))
          end do
60  format( 7x,i3,', ',i3,' -> ',i6,2x,f15.10,2x,f15.10 )
61  format( 7x,i3,' -> ',i6,2x,f15.10,2x,f15.10 )
    end do

    write(iout,'(A)') ' ENERGIES AND TRANSITION MATRIX ELEMENTS'
    write(iout,41) '#', '<i|H0|i> eV', 'i   j','<i|mu_x|j> au', '<i|mu_y|j> au', '<i|mu_z|j> au', '<i|Vabs|j> au'
    do i=1, istates
       ii = i
       write(iout,51) i, cis_eig(i)*au2eV, 1, i, tdx(i), tdy(i), tdz(i), abp(i)
       ii = ii + nstates
       if(i.gt.1) write(iout,52) 2,i,tdx(ii), tdy(ii), tdz(ii), abp(ii)
       ii = ii + nstates
       if(i.gt.2) write(iout,52) 3,i,tdx(ii), tdy(ii), tdz(ii), abp(ii)
       ii = ii + nstates
       if(i.gt.3) write(iout,52) 4,i,tdx(ii), tdy(ii), tdz(ii), abp(ii)
       ii = ii + nstates
       if(i.gt.4) write(iout,52) 5,i,tdx(ii), tdy(ii), tdz(ii), abp(ii)
       ii = ii + nstates
       if(i.gt.5) write(iout,52) 6,i,tdx(ii), tdy(ii), tdz(ii), abp(ii)
    end do

!:    do i=1,nstates
!:      rdum = 0.d0
!:      do ia=1,8
!:        rdum=rdum+cis_vec(ia+(i-1)*nstates)**2
!:      end do
!:      if(rdum.gt.1.d-3) write(iout,"(i5,f12.4,9f10.6)") i,cis_eig(i)*au2ev,rdum,(cis_vec(ia+(i-1)*nstates),ia=1,8)
!:    end do

40  format( a5,1x,5(1x,a15)/6x,6(1x,a15) )
41  format( a5,1x,a15,a8,4(1x,a15) )
50  format( i5,1x,5(1x,f15.10)/6x,6(1x,f15.8) )
51  format( i5,1x,f15.10,2i4,4(1x,f15.10) )
52  format( 21x,2i4,4(1x,f15.10) )

    call write_header( 'write_specifics2','write_info','leave' )
    return

78  continue
    call Zwrite_specifics2
    
    
  end subroutine write_specifics2
  !: ----------------------------- :!
  !: subroutine Zwrite_specific2    :!
  !: ----------------------------- :!
  subroutine Zwrite_specifics2


    implicit none
    integer(8) :: i, j, ii, iii, ia, jb, istates
    real(8) :: szestm, s2estm, dip2
    
    call write_header( 'Zwrite_specifics2','write_info','enter' )


    write(iout,'(A)') ' NSTUSE' 
    write(iout,"(5x,'NstUse          = ',i12)") nstuse
    write(iout,"(5x,'Nstates         = ',i12)") nstates
    write(iout,"(5x,'cis_eig(nstuse) = ',f12.6,' au',f12.6,' eV')") cis_eig(nstuse), cis_eig(nstuse)*au2eV
    
    istates = 20
    write(iout,'(A)') ' ENERGIES (eV)'
    do i=1, istates
       write(iout,50) i, cis_eig(i)*au2eV
    end do

    write(iout,'(A)') ' ENERGIES and EXPECTATION VALUES'
    write(iout,40) '#  ', '<i|H0|i> eV', ' Sz   ', &
         'R<i|mu_x|i> au', 'I<i|mu_x|i> au', &
         'R<i|mu_y|i> au', 'I<i|mu_y|i> au', &
         'R<i|mu_z|i> au', 'I<i|mu_z|i> au', &
         'R<i|Vabs|i> au', 'I<i|mu_z|i> au', &
         '<i|a_xx|i> eV', '<i|a_yy|i> au', '<i|a_zz|i> au', &
         '<i|a_xy|i> au', '<i|a_xz|i> au', '<i|a_yz|i> au'
    call get_Zpolar(istates,nstates,Ztdx,Ztdy,Ztdz,cis_eig,polar)
    do i=1, istates
       ii = (i-1)*nstates
       iii = (i-1)*nstuse + i
        if( trim(jobtype).eq.flag_soc ) then
          szestm = 0.d0
          s2estm = 0.d0
          do ia = 1, nstates
            if(hole_index(ia,1).lt.0.and.part_index(ia,1).gt.0) &
              szestm=szestm-abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0.and.part_index(ia,1).lt.0) &
              szestm=szestm+abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1)*part_index(ia,1).lt.0) &
              s2estm=s2estm+2.d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0.and.part_index(ia,1).gt.0) &
              s2estm=s2estm+abs(Zcis_vec(ia+ii)-Zcis_vec(ia+ii-noanva))**2 
          end do
        end if
        if( trim(jobtype).eq.flag_socip ) then
          szestm = 0.d0
          s2estm = 0.75d0
          do ia = 1, nstates
            if(hole_index(ia,1).lt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,2).lt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,2).gt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(part_index(ia,1).lt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(part_index(ia,1).gt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).lt.0.and.hole_index(ia,2).lt.0.and.part_index(ia,1).gt.0) &
              s2estm=s2estm+3.d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0.and.hole_index(ia,2).gt.0.and.part_index(ia,1).lt.0) &
              s2estm=s2estm+3.d0*abs(Zcis_vec(ia+ii))**2 
!:            if(hole_index(ia,1).gt.0.and.hole_index(ia,2).lt.0.and.part_index(ia,1).lt.0) &
!:              s2estm=s2estm+1.5d0*abs(Zcis_vec(ia+ii)-Zcis_vec(ia+ii-noa*nobnvb))**2 
          end do
        end if
        write(iout,49) i, cis_eig(i)*au2eV, szestm, Ztdx(iii), Ztdy(iii), Ztdz(iii), Zabp(iii), &
            (polar(i,j),j=1,6)
          do ia=1, nstates
             if ( abs(Zcis_vec(ia+nstates*(i-1))).ge.0.01D0 ) then
               if ( trim(jobtype).eq.flag_soc ) &
                 write(iout,60) hole_index(ia,1),  part_index(ia,1), Zcis_vec(ia+nstates*(i-1))
               if ( trim(jobtype).eq.flag_socip ) &
                 write(iout,61) hole_index(ia,1), hole_index(ia,2), part_index(ia,1), Zcis_vec(ia+nstates*(i-1))
            end if
          end do
    end do

    write(iout,'(A)') ' ENERGIES AND TRANSITION MATRIX ELEMENTS'
    write(iout,41) '#', '<i|H0|i> eV', ' Sz     i,j', 'TransDipole**2 ', &
         'R<i|mu_x|j> au', 'I<i|mu_x|j> au', &
         'R<i|mu_y|j> au', 'I<i|mu_y|j> au', &
         'R<i|mu_z|j> au', 'I<i|mu_z|j> au', &
         'R<i|Vabs|j> au', 'I<i|Vabs|j> au'
    do i=1, istates
        ii = (i-1)*nstates
        if( trim(jobtype).eq.flag_soc ) then
          szestm = 0.d0
          do ia = 1, nstates
            if(hole_index(ia,1).lt.0.and.part_index(ia,1).gt.0) &
              szestm=szestm-abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0.and.part_index(ia,1).lt.0) &
              szestm=szestm+abs(Zcis_vec(ia+ii))**2 
          end do
        end if
        if( trim(jobtype).eq.flag_socip ) then
          szestm = 0.d0
          do ia = 1, nstates
            if(hole_index(ia,1).lt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,1).gt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,2).lt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(hole_index(ia,2).gt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(part_index(ia,1).lt.0) szestm=szestm+0.5d0*abs(Zcis_vec(ia+ii))**2 
            if(part_index(ia,1).gt.0) szestm=szestm-0.5d0*abs(Zcis_vec(ia+ii))**2 
          end do
        end if
       dip2 = Dconjg(Ztdx(i))*Ztdx(i) &
            + Dconjg(Ztdy(i))*Ztdy(i) &
            + Dconjg(Ztdz(i))*Ztdz(i)
       write(iout,53) i, cis_eig(i)*au2eV, szestm, i, 1, dip2, Ztdx(i), Ztdy(i), Ztdz(i), Zabp(i)
       If(i.gt.1) then
         ii = i + nstuse
         dip2 = Dconjg(Ztdx(ii))*Ztdx(ii) &
              + Dconjg(Ztdy(ii))*Ztdy(ii) &
              + Dconjg(Ztdz(ii))*Ztdz(ii)
         write(iout,51) i,2,dip2,Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
         end if
       If(i.gt.2) then
         ii = i + 2*nstuse
         dip2 = Dconjg(Ztdx(ii))*Ztdx(ii) &
              + Dconjg(Ztdy(ii))*Ztdy(ii) &
              + Dconjg(Ztdz(ii))*Ztdz(ii)
         write(iout,51) i,3,dip2,Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
         end if
       If(i.gt.3) then
         ii = i + 3*nstuse
         dip2 = Dconjg(Ztdx(ii))*Ztdx(ii) &
              + Dconjg(Ztdy(ii))*Ztdy(ii) &
              + Dconjg(Ztdz(ii))*Ztdz(ii)
         write(iout,51) i,4,dip2,Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
         end if
       If(i.gt.4) then
         ii = i + 3*nstuse
         dip2 = Dconjg(Ztdx(ii))*Ztdx(ii) &
              + Dconjg(Ztdy(ii))*Ztdy(ii) &
              + Dconjg(Ztdz(ii))*Ztdz(ii)
         write(iout,51) i,5,dip2,Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
         end if
       If(i.gt.5) then
         ii = i + 3*nstuse
         dip2 = Dconjg(Ztdx(ii))*Ztdx(ii) &
              + Dconjg(Ztdy(ii))*Ztdy(ii) &
              + Dconjg(Ztdz(ii))*Ztdz(ii)
         write(iout,51) i,6,dip2,Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
         end if
        end do

40  format( a5,1x,a15,3x,a8,8(1x,a15)/32x,6(1x,a15) )
41  format( a5,1x,a15,a11,2x,9(1x,a15) )
49  format( i5,1x,f15.10,f10.2,8(1x,f15.10)/31x,6(1x,f15.10) )
50  format( i5,1x,f15.10,f5.2,i4,i2,8(1x,f15.10)/31x,6(1x,f15.10) )
51  format( 26x,i4,i2,9(1x,f15.10) )
52  format( 9(1x,f15.10) )
53  format( i5,1x,f15.10,f5.2,i4,i2,9(1x,f15.10) )
60  format( 7x,i3,' -> ',i6,2x,f15.10,2x,f15.10 )
61  format( 7x,i3,', ',i3,' -> ',i6,2x,f15.10,2x,f15.10 )
70  format( 3x,' S**2 = ',f5.2,4x,' Sz = ',f5.2 )
80  format( 3x,' Sz = ',f5.2 )

    call write_header( 'Zwrite_specifics2','write_info','leave' )
    
    
  end subroutine Zwrite_specifics2
  !: ----------------------------- :!
  !: SUBROUTINE WRITE_SUMMARY      :!
  !: ----------------------------- :!
  subroutine write_summary

    implicit none

    integer(8) :: ifield
    character(1000) :: cformat

    call write_header( 'write_summary', 'write_info', 'enter' )
    
    write(iout,'(A)') ' SUMMARY'
    cformat = '    #'//'   e_max'//'         x'//'         y'//'         z' //'      norm'//'      dipx'//'      dipy'//'      dipz'

    write(iout,'(A)') trim(cformat)
    do ifield=1, nemax*ndir
       write(iout,"(i5,1x,f7.5,100(f10.4))") ifield, &
            tdciresults(ifield)%fstrength0,       &
            tdciresults(ifield)%x0, &
            tdciresults(ifield)%y0, &
            tdciresults(ifield)%z0, &
            tdciresults(ifield)%norm0, &
            tdciresults(ifield)%dipx,  &
            tdciresults(ifield)%dipy,  &
            tdciresults(ifield)%dipz
    end do
    

    call write_header( 'write_summary', 'write_info', 'leave' )


  end subroutine write_summary
  !: ----------------------------- :!
  !: subroutine write_mo_energies  :! 
  !: ----------------------------- :!
  subroutine write_mo_energies


    implicit none


    character(100) :: myout 
    integer(8) :: i

    
    call write_header( 'write_mo_energies','write_info','enter' )

    
    myout = trim(outputfile)//'_MO_ENERGIES'
    open( unit=100,file=trim(myout) )


    !: restricted MOs
    if ( .not.unrestricted )  then
       write(100,"(' HOMO =', i0, 5x, ' LUMO =', i0)") noa, noa+1
       write(100,"(a10,a20,1x,a20)") ' # MO' , 'restricted (au)' , 'restricted (eV)'
       do i=1, nrorb
          write(100,"(i10,2(f20.10,1x))") i, orben(i), orben(i)*au2eV
       end do
       
    !: unrestricted MOs
    else if ( unrestricted ) then
       write(100,"(' alpha HOMO = ', i0, 5x, ' alpha LUMO = ', i0)") noa, noa+1
       write(100,"(' beta  HOMO = ', i0, 5x, ' beta  LUMO = ', i0)") nob, nob+1
       write(100,"(a10,4(a20,1x))") ' # MO','alpha (au)','alpha (eV)','beta (au)','beta (eV)'
       do i=1, nrorb
          write(100,"(i10,4(f20.10,1x))") i, orben(i), orben(i)*au2eV, orben(i+nrorb), orben(i+nrorb)*au2eV
       end do
    end if

    close(100)
    write(iout,'(A)') " MO orbital energies written out to file "//"'"//trim(myout)//"'"


    call write_header( 'write_mo_energies','write_info','leave' )

    
  end subroutine write_mo_energies
  !: ----------------------------- :!
  !: subroutine write_field_shape  :!
  !: ----------------------------- :!
  subroutine write_field_shape
  
    
    implicit none
    
    
    character(100) :: myout_shape 
    character(100) :: myout_dir   
    integer(8) :: i, istp
    real(8)    :: rdum, one

    call write_header( 'write_field_shape','write_info','enter' )

    myout_shape = trim(outputfile)//'_FIELD_SHAPE'
    myout_dir   = trim(outputfile)//'_FIELD_DIR'
    
    field_shape : if ( .not.linear ) then

       open( unit=100,file=trim(myout_shape) )

!: something weird with the compilation for some machines in calculating 1.d0/dsqrt(1.d0+ellipt**2)
!: seems to work if split into several steps with writes and calls in between
    one = 1.d0
    rdum = ellipt
    write(42,*)nstep,dt,autime2s,s2fs,ellipt,rdum,fvect1(1),fvect1(nstep),fvect2(1),fvect2(nstep)
    flush(42)
    rdum=rdum**2
    write(42,*) " rdum1",rdum
    flush(42)
    rdum=rdum+one
    rdum=sqrt(rdum)
       write(100,"( a10,a10,a10,a20,a20,a20 )") 'step#','time(au)','time(fs)','envelope','E1(t)','E2(t)'
       do istp=1, nstep
          write(100,100) istp, dble(istp)*dt, dble(istp)*dt*autime2s*s2fs, env(istp), &
               fvect2(istp)/rdum, fvect1(istp)/rdum
       end do

       close(100)
       
    else
       
       open( unit=100,file=trim(myout_shape) )

       write(100,"( a10,a10,a10,a20,a20 )") 'step#','time(au)','time(fs)','envelope','E1(t)'       
       do istp=1, nstep
          write(100,200) istp, dble(istp)*dt, dble(istp)*dt*autime2s*s2fs, env(istp), fvect1(istp)
       end do

       close(100)
              
    end if field_shape

100 format( i10,f10.2,f10.4,e20.4,e20.4,e20.4 )
200 format( i10,f10.2,f10.4,e20.4,e20.4 )
    
    
    write(iout, '(A)') " pulse shape written to file "//"'"//trim(myout_shape)//"'"
    flush(iout)
    

    !: field direction
    open( unit=101, file=trim(myout_dir) )

    field_dir : if ( envelope.ne.'cirl' .and. envelope.ne.'cirr' ) then

       write(101,50) '#' , 'E_0(au)' , 'x0' , 'y0' , 'z0' , 'theta0' , 'phi0'
       flush(101)
       do i=1, nemax*ndir
          write(101,60) i, &
               tdciresults(i)%fstrength0, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%theta0, tdciresults(i)%phi0
          flush(101)
       end do
       
    else

       write(101,70) '#' , 'E_0(au)' , 'E_1(au)' , 'E_2(au)' ,           &
            'x0' , 'y0' , 'z0' , 'x1' , 'y1' , 'z1' , 'x2' , 'y2' , 'z2' , &
            'theta0' , 'phi0' , 'theta1' , 'phi1' , 'theta2' , 'phi2'
       flush(101)
       do i=1, nemax*ndir
          write(101,80) i, &
               tdciresults(i)%fstrength0, tdciresults(i)%fstrength1, tdciresults(i)%fstrength2, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%x1, tdciresults(i)%y1, tdciresults(i)%z1, &
               tdciresults(i)%x2, tdciresults(i)%y2, tdciresults(i)%z2, &
               tdciresults(i)%theta0, tdciresults(i)%phi0, &
               tdciresults(i)%theta1, tdciresults(i)%phi1, &
               tdciresults(i)%theta2, tdciresults(i)%phi2
          flush(101)
       end do
       
    end if field_dir

    
50  format( 1x, a4, a8,   2x, '(',a9  ,1x,a9,  1x,a9,  1x,')',2x,'(',a9  ,1x,a9,  1x,')' )
60  format( 1x, i4, f8.4, 2x, '(',f9.4,1x,f9.4,1x,f9.4,1x,')',2x,'(',f9.4,1x,f9.4,1x,')' )
70  format( 1x,a4,a8,  2x,a8,  2x,a8,  2x,3('(',a9  ,1x,a9  ,1x,a9  ,1x,')',2x), 3(' (',a9,  1x,a9,  1x,')',2x) )
80  format( 1x,i4,f8.4,2x,f8.4,2x,f8.4,2x,3('(',f9.4,1x,f9.4,1x,f9.4,1x,')',2x), 3(' (',f9.4,1x,f9.4,1x,')',2x) )

    
    write(iout, '(A)') " field direction written to file "//"'"//trim(myout_dir)//"'"
    

    call write_header( 'write_field_shape','write_info','leave' )


  end subroutine write_field_shape
  !: ----------------------- :!
  !: SUBROUTINE WRITE_HAM0   :! 
  !: ----------------------- :!
  subroutine write_ham0

    implicit none

    real(8), parameter :: thres  = 0.00d0
    real(8), parameter :: thres1 = thres * dsqrt(2.d0)
    
    integer(8)     :: ia, ii, aa, xx, jj, bb, istate, i, j, a, b
    real(8), allocatable    :: coeffs(:,:)
    complex(8), allocatable :: Zcoeffs(:,:)
    character(5)   :: cii, caa, cxx, aorb, xorb, cjj, cbb
    character(100) :: myout
    
    
    call write_header( 'write_ham0','write_info','enter' )

    myout = trim(outputfile)//'_HAM0'

    open( unit=100,file=trim(myout) )

    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
       allocate( Zcoeffs(nstates,nstates) )
       Zcoeffs = reshape( Zcis_vec, (/ nstates, nstates/) )
    else
       allocate( coeffs(nstates,nstates) )
       coeffs = reshape( cis_vec,(/ nstates, nstates /) )
    end if

    select case ( trim(jobtype) )
    case( flag_cis , flag_tda )
       
       if ( .not.unrestricted ) then
          restricted : do i=1, nstates
             write(100,100) i, cis_eig(i)*au2eV
             do ia=1, nstates
                if( abs(coeffs(ia,i)).ge.thres1 ) &
                     write(100,200) abs(hole_index(ia,1)), abs(part_index(ia,1)), coeffs(ia,i)/sqrt(2.d0)
             end do
          end do restricted
       else
          unrestricted : do i=1, nstates
             write(100,100) i, cis_eig(i)*au2eV
             do ia=1, nstates
                if ( abs(coeffs(ia,i)).ge.thres ) &
                     write(100,200) hole_index(ia,1), part_index(ia,1), coeffs(ia,i)
             end do
          end do unrestricted
       end if
       
    case( flag_soc )
       do i=1, nstates
          write(100,100) i, cis_eig(i)*au2eV
          do ia=1, nstates
             if ( abs(Zcoeffs(ia,i)).ge.thres ) &
                  write(100,200) hole_index(ia,1), part_index(ia,1), Zcoeffs(ia,i)
          end do
       end do
       
    case( flag_ip )

       ip : do i=1, nstates
          write(100,100) i, cis_eig(i)*au2eV
          do ia=1, nstates
             if ( abs(coeffs(ia,i)).ge.thres ) then
                xx = hole_index(ia,1) 
                ii = hole_index(ia,2) 
                aa = part_index(ia,1) 
                write(100,310) xx, ii, aa, coeffs(ia,i)
             end if
          end do
       end do ip

    case( flag_socip )
       do i=1, nstates
          write(100,100) i, cis_eig(i)*au2eV
          do ia=1, nstates
             if ( abs(Zcoeffs(ia,i)).ge.thres ) &
                  write(100,320) hole_index(ia,1), hole_index(ia,2), part_index(ia,1), Zcoeffs(ia,i)
          end do
       end do
       
    case( flag_cisd ) 
       
       do ia=1, nstates


          write(100,100) ia, cis_eig(ia)*au2eV

          if( abs(coeffs(1,ia).ge.thres) ) write(100,350) ii, aa, jj, bb, coeffs(1,ia)
          
          istate = 1 
          
          alpha_singles : do i=1, noa
             do a=1, nva
                istate = istate + 1 
                if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) -i, -a, 0, 0, coeffs(istate,ia)
             end do
          end do alpha_singles

          beta_singles : do i=1, nob
             do a=1, nvb
                istate = istate + 1 
                if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) i, a, 0, 0, coeffs(istate,ia)
             end do
          end do beta_singles

          alpha_beta_doubles : do i=1, noa
             do j=1, nob
                do a=1, nva
                   do b=1, nvb
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) j, b, -i, -a, coeffs(istate,ia)
                   end do
                end do
             end do
          end do alpha_beta_doubles

          alpha_alpha_doubles : do i=1, noa
             do j=(i+1), noa
                do a=1, nva
                   do b=(a+1), nva
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) -i, -a, -j, -b, coeffs(istate,ia)
                   end do
                end do
             end do
          end do alpha_alpha_doubles
          
          beta_beta_doubles : do i=1, nob
             do j=(i+1), nob
                do a=1, nvb
                   do b=(a+1), nvb
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) i, a, -j, -b, coeffs(istate,ia)
                   end do
                end do
             end do
          end do beta_beta_doubles
          
       end do
       
    end select

    close(100)
    
100 format(' Excited state ',i10,' : ',f20.10,' eV ')
200 format( 7x,i3,' -> ',i3,2x,f15.10,2x,f15.10 )
300 format( 7x,i5,' -> ',i5,2x,f10.7 )
310 format( 7x,i5,' -> ','inf',2x,i5,' ->',i5,2x,f10.7 )
320 format( 7x,i5,' -> ','inf',2x,i5,' ->',i5,2x,f15.10,2x,f15.10 )
350 format( 7x,i5,' -> ',i5,2x,i5,' ->',i5,2x,f10.7 )
    
    write(iout,'(A)') ' finished writing eigen-info to file '//"'"//trim(myout)//"'"

    call write_header( 'write_ham0','write_info','leave' )
    

  end subroutine write_ham0
  !: ---------------------------- :!
  !: SUBROUTINE WRITE_MO_ELEMENTS :!
  !: ---------------------------- :!
  subroutine write_mo_elements
  
    implicit none

    integer(8)     :: i, j
    character(100) :: myout, myout2
    character(1)   :: ov1, ov2 !: occupied or virtual

    
    call write_header( 'write_mo_elements','write_info','enter' )


    if( .not. unrestricted ) then
       myout = trim(outputfile)//'_MO_ELEMENTS'
       open( unit=100,file=trim(myout) )
    else
       myout  = trim(outputfile)//'_MO_ELEMENTS_ALPHA'
       myout2 = trim(outputfile)//'_MO_ELEMENTS_BETA'
       open( unit=100,file=trim(myout) )
       open( unit=200,file=trim(myout2) )
    end if

    
    write(100,50) 'i','j', '<i|mu_x|j>(au)', '<i|mu_y|j>(au)', '<i|mu_z|j>(au)', '<i|Vabs|j>(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=i, (noa+nva)
          ov2 = 'v' ; if ( j.le.noa ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, dipxmoa(j,i), dipymoa(j,i), dipzmoa(j,i), vabsmoa(j,i) 
       end do
    end do
    close(100)
    
    write(iout,'(A)') " MO matrix elements written out to file "//"'"//trim(myout)//"'"

    if ( unrestricted ) then

       write(200,50) 'i','j', '<i|mu_x|j>(au)', '<i|mu_y|j>(au)', '<i|mu_z|j>(au)', '<i|Vabs|j>(au)'
       do i=1, (nob+nvb)
          ov1 = 'v' ; if ( i.le.nob ) ov1 = 'o'
          do j=i, (nob+nvb) 
             ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
             write(200,100) i, ov1, j, ov2, dipxmob(j,i), dipymob(j,i), dipzmob(j,i), vabsmob(j,i)
          end do
       end do
       close(200)
       write(iout,'(A)') " MO matrix elements written out to file "//"'"//trim(myout2)//"'"
       
    end if


    if( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) call write_soc_elements

             
50  format( 2(a7),4(a20) )
100 format( 2(i5,a2,1x),4(f20.10) )
    
    call write_header( 'write_mo_elements','write_info','leave' )
    

  end subroutine write_mo_elements
  !:-----------------------------:!
  !: SUBROUTINE WRITE_SOC_ELEMENTS
  !:-----------------------------:!
  subroutine write_soc_elements

    implicit none

    integer(8)     :: i, j
    character(100) :: myout
    character(1)   :: ov1, ov2 !: occupied or virtual

    
    call write_header( 'write_soc_elements','write_info','enter' )
    

    myout = trim(outputfile)//'_SOC_ELEMENTS_AA'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_alpha', 'j_alpha', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (noa+nva)
          ov2 = 'v' ; if ( j.le.noa ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoAA(j,i)
       end do
    end do
    close(100)

    write(iout,'(A)') " alpha alpha SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"


    myout = trim(outputfile)//'_SOC_ELEMENTS_BB'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_beta', 'j_beta', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (nob+nvb)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (nob+nvb)
          ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoBB(j,i)
       end do
    end do
    close(100)
    
    write(iout,'(A)') " beta  beta  SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"


    myout = trim(outputfile)//'_SOC_ELEMENTS_AB'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_alpha', 'j_beta', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (nob+nvb)
          ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoAB(j,i)
       end do
    end do
    close(100)

    write(iout,'(A)') " alpha beta  SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"    


50  format( 2(a7,1x),2(a20) )
100 format( 2(i5,a2,1x),2(f20.10) )
    

  end subroutine write_soc_elements
  !: --------------------------- :!
  !: SUBROUTINE SAVE_RESTART_BIN
  !: --------------------------- :!
  subroutine save_restart_bin


    implicit none

    integer(8), parameter :: zero = 0

    character(100) :: myout
    integer(8) :: i, j, ij, iflag


    call write_header( 'save_restart_bin','write_info','enter' )


    myout = trim(outputfile)//'_RESTART.bin'
    open( unit=50,file=trim(myout),form='unformatted' )


   !:  restart not available for soc_cis yet
   !:  if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
   !:    write(iout,'(A)') ' WARNING:  RESTART NOT AVAILABLE FOR SOC_CIS or SOCIP_CISD'
   !:    call write_header( 'save_restart_bin','write_info','leave' )
   !:    return
   !:  end if
    
    
    write(50) nstates, nstuse
    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
      write(50) Zcis_vec
      write(50) cis_eig
      write(50) Ztdx
      write(50) Ztdy
      write(50) Ztdz
      write(50) Zabp
      write(50) Zexp_abp
      write(50) Zip_vec
    else
      write(50) cis_vec
      write(50) cis_eig
      write(50) tdx
      write(50) tdy
      write(50) tdz
      write(50) abp
      write(50) exp_abp
      write(50) ip_vec
    end if
    
    close(50)
    write( iout,'(A)' ) ' restart binary file written out to file '//"'"//trim(myout)//"'"
    write( iout,'(A)' ) ' to read binary file '
    write( iout,"(5x,'nstates, nstuse')" ) 
    write( iout,"(5x,'cis_vec, cis_eig, tdx, tdy, tdz, abp, exp_abp')") 
    

    myout = trim(outputfile)//'_RESTART_MO.bin'
    open( unit=50, file=trim(myout), form='unformatted' )

    if( unrestricted ) write(50) noa, nva, nob, nvb
    if( .not.unrestricted ) write(50) noa, nva, zero, zero
    write(50) ( vabsmoa(:,i), i=1, nrorb )
    write(50) ( dipxmoa(:,i), i=1, nrorb )
    write(50) ( dipymoa(:,i), i=1, nrorb )
    write(50) ( dipzmoa(:,i), i=1, nrorb )
    if ( unrestricted ) then
       write(50) ( vabsmob(:,i), i=1, nrorb )
       write(50) ( dipxmob(:,i), i=1, nrorb )
       write(50) ( dipymob(:,i), i=1, nrorb )
       write(50) ( dipzmob(:,i), i=1, nrorb )
    end if
    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) then
       write(50) ( socmoAA(:,i), i=1, nrorb )
       write(50) ( socmoBB(:,i), i=1, nrorb )
       write(50) ( socmoAB(:,i), i=1, nrorb )
    end if
    
    close(50)
    write( iout,'(A)' ) ' restart binary file written out to file '//"'"//trim(myout)//"'"
    write( iout,'(A)' ) ' to read binary file '
    write( iout,"(5x,'noa, nva, nob, nvb (nob=nvb=0 for restricted)' )" ) 
    write( iout,"(5x,'vabsmoa, dipxmoa, dipymoa, dipzmoa')")
    if ( unrestricted ) &
      write( iout,"(5x,'vabsmob, dipxmob, dipymob, dipzmob')")
    if ( trim(jobtype).eq.flag_soc .or. trim(jobtype).eq.flag_socip ) &
      write( iout,"(5x,'socmoAA, socmoBB, socmoAB')")
    
    call write_header( 'save_restart_bin','write_info','leave' )


  end subroutine save_restart_bin
  !: ------------------------ :!
  !: END MODULE WRITE_INFO    :!
  !: ------------------------ :!
end module write_info
