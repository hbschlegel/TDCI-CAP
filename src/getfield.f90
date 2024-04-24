module getfield

  use variables_global
  implicit none  

  !: contains subroutine  check_emax
  !: contains subroutine  get_lindirection
  !: contains subroutine  get_circdirection
  !: contains subroutine  shape_field
  !: contains subroutines null_env, cos_env, trap_env, gau_env, cw_env, cs_env, sin2_env, circ_env
  !: contains subroutine  writeme_getfield
  !: contains subroutine  errors_getfield
  
  
contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine check_emax


    !: check if enput is here
    implicit none
    
    real(8)    :: moenergy
    real(8)    :: guessfield, dphi, dtheta
    integer(8) :: nphi, ntheta, iphi, itheta, i

    
    call write_header( 'check_emax','getfield','enter' )
    
    
    if( read_emax(1).lt.-100.0 ) then
       call errors_getfield( 'check_emax','no_enput' )
       nemax = 3 
       moenergy = dabs( Mol%orben(noa) )
       if( unrestricted ) moenergy = dabs( Mol%orben(nrorb+nob) )

       !: guessed field strength on Keldysh parameter gamma = 1.0
       guessfield = omega * dsqrt( 2.d0 * moenergy )
       open(unit=100,file=trim(fieldfile))
       if ( envelope.ne.'cirl' .and. envelope.ne.'cirr' ) then
          read_emax(1) = guessfield - 0.01
          read_emax(2) = guessfield
          read_emax(3) = guessfield + 0.01
          call errors_getfield( 'check_emax','no_lin' )
       else
          guessfield = guessfield * dsqrt(2.d0)
          read_emax(1) = guessfield - 0.01
          read_emax(2) = guessfield
          read_emax(3) = guessfield + 0.01
          call errors_getfield( 'check_emax','no_cir' )
       end if
       
    else
       
       write(iout,"(' emax #',i3,' = ', f7.4, ' au')" ) (i, read_emax(i), i=1, nemax )

    end if
    
    
    call write_header( 'check_emax','getfield','leave' )

    
  end subroutine check_emax
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_LINDIRECTION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_lindirection
    

    !: read E-field propgation in from 'enput' ; record in polar coordinates

    
    implicit none    

    
    integer(8)  :: ifield, idir, ij, total, mystat
    real(8)     :: tmpx, tmpy, tmpz, dnorm, theta, phi
    
    
    call write_header( 'get_lindirection','getfield','enter' )

    
    total = nemax * ndir

    
    if ( trim(dirform).eq.'cart' ) then

       !: convert Cartesian to Polar
       cart : do idir=1, ndir       
          tmpx = read_x(idir) 
          tmpy = read_y(idir)
          tmpz = read_z(idir)
          dnorm = dsqrt( tmpx**2 + tmpy**2 + tmpz**2 )
          theta = dacos(tmpz/dnorm) * rad2deg
          phi   = datan2(tmpy,tmpx) * rad2deg
          do ifield=1, nemax
             ij = (ifield-1)*ndir + idir
             tdciresults(ij)%fstrength0 = read_emax(ifield)
             tdciresults(ij)%theta0     = theta
             tdciresults(ij)%phi0       = phi
             tdciresults(ij)%x0         = tmpx / dnorm
             tdciresults(ij)%y0         = tmpy / dnorm
             tdciresults(ij)%z0         = tmpz / dnorm
          end do
       end do cart
       
    else if ( trim(dirform).eq.'polar' ) then
       
       !: read in direction of field-propagation in polar coords
       polar : do idir=1, ndir 
          theta = read_theta(idir)
          phi   = read_phi(idir)
          theta = theta * deg2rad
          phi   = phi   * deg2rad
          tmpx  = dsin(theta) * dcos(phi)
          tmpy  = dsin(theta) * dsin(phi)
          tmpz  = dcos(theta)
          dnorm = sqrt( tmpx**2 + tmpy**2 + tmpz**2 )
          do ifield=1, nemax
             ij = (ifield-1)*ndir + idir
             tdciresults(ij)%fstrength0 = read_emax(ifield)
             tdciresults(ij)%theta0     = theta * rad2deg
             tdciresults(ij)%phi0       = phi   * rad2deg
             tdciresults(ij)%x0         = tmpx / dnorm
             tdciresults(ij)%y0         = tmpy / dnorm
             tdciresults(ij)%z0         = tmpz / dnorm
          end do
       end do polar
       
    end if
    
    write(iout,'(A)') ' finished assigning datatype tdciresults'
    call write_header( 'get_lindirection','getfield','leave' )
    

  end subroutine get_lindirection
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_CIRCDIRECTION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_circdirection
 
   
    !: for circularly polarized light. 
    !: theta,phi correspond to direction of light propagation
    !: two field vectors are perpendicular to each other and to direction of light prop.
    !: https://www.youtube.com/watch?v=tmtGEHTBSdQ
    !: according to Euler's rotation theorem
    implicit none
    

    integer(8) :: ifield, idir, ij
    real(8)    :: tmpx, tmpy, tmpz, tmp, theta, phi
    real(8)    :: alltheta(ndir), allphi(ndir)


    call write_header( 'get_circdirection','getfield','enter' )
    write(iout,'(A)') ' generating E1 and E2 for each propagation direcion '


    !: retrieve x, y, z vectors from get_linvector to light prop direction 
    !: define as E_x in field frame
    alltheta(:) = tdciresults(1:ndir)%theta0 * deg2rad
    allphi(:)   = tdciresults(1:ndir)%phi0   * deg2rad    
    euler       = euler * deg2rad
    
    edir : do idir=1, ndir 
       estrength : do ifield=1, nemax
          ij = (ifield-1)*ndir + idir 
          theta = alltheta(idir)
          phi   = allphi(idir)
          !: first field vector
          tdciresults(ij)%fstrength1 = tdciresults(ij)%fstrength0 / dsqrt( 1.d0 + ellipt**2 )
          tdciresults(ij)%x1         = dcos(phi)*dcos(theta)*dcos(euler) - dsin(phi)*dsin(euler)
          tdciresults(ij)%y1         = dsin(phi)*dcos(theta)*dcos(euler) + dcos(phi)*dsin(euler)
          tdciresults(ij)%z1         = -dsin(theta)*dcos(euler)
          tdciresults(ij)%theta1     = dacos(tdciresults(ij)%z1) * rad2deg
          tdciresults(ij)%phi1       = datan2(tdciresults(ij)%y1,tdciresults(ij)%x1) * rad2deg
          if ( tdciresults(ij)%theta1.eq.0.d0 )   tdciresults(ij)%phi1 = 0.d0 
          if ( tdciresults(ij)%theta1.eq.180.d0 ) tdciresults(ij)%phi1 = 0.d0 
          !: second field vector from cross product 
          tdciresults(ij)%fstrength2 = tdciresults(ij)%fstrength0 * ellipt / dsqrt( 1.d0 + ellipt**2 )
          tdciresults(ij)%x2         = -dcos(phi)*dcos(theta)*dsin(euler) - dsin(phi)*dcos(euler)
          tdciresults(ij)%y2         = -dsin(phi)*dcos(theta)*dsin(euler) + dcos(phi)*dcos(euler) 
          tdciresults(ij)%z2         = dsin(theta)*dsin(euler)
          tdciresults(ij)%theta2     = dacos(tdciresults(ij)%z2) * rad2deg
          tdciresults(ij)%phi2       = datan2(tdciresults(ij)%y2,tdciresults(ij)%x2) * rad2deg
          if ( tdciresults(ij)%theta2.eq.0.d0 )   tdciresults(ij)%phi2 = 0.d0
          if ( tdciresults(ij)%theta2.eq.180.d0 ) tdciresults(ij)%phi2 = 0.d0
       end do estrength
    end do edir
    
    !: store euler in degrees
    euler = euler * rad2deg
    
    write(iout,'(A)') ' finished setting up E_1 and E_2 '
    call write_header( 'get_circdirection','getfield','leave' )
    
    
  end subroutine get_circdirection
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SHAPE_FIELD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine shape_field


    !: period and field_duration set in SUBROUTINE set_variables in module initialize 
    implicit none
    

    call write_header( 'shape_field','getfield','enter' )
    
    select case ( trim(envelope) ) 
    case( 'none' ) ; call null_env
    case( 'cos2' ) ; call cos_env
    case( 'trap' ) ; call trap_env
    case( 'gaus' ) ; call gau_env
    case( 'stat' ) ; call cw_env
    case( 'band' ) ; call cs_env
    case( 'ramp' ) ; call ramp_env
    case( 'sin2' ) ; call sin2_env
    case( 'cirl' ) ; call circ_env(-1.d0, 1.d0 )
    case( 'cirr' ) ; call circ_env( 1.d0, 1.d0 )
    end select
    
    
        
    call write_header( 'shape_field','getfield','leave' )
    

  end subroutine shape_field
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE NULL_ENV 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine null_env
    
    !: envelope=0 : no pulse
    implicit none
    
    write(iout,'(A)') ' in subroutine null_env for null pulse'
    fvect1 = 0.d0
    write(iout, '(A)') ' leaving subroutine null_env'
    
  end subroutine null_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE COS_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine cos_env
    
    !: envelope=1 : generates *Cosine Envelope* pulse
    implicit none
    
    integer(8) :: istp
    real(8)    :: tau
    
    write(iout,'(A)') ' in subroutine cos_env for cosine envelope'           
    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0
    tau = 2.d0*pi / field_duration

    do istp=1, nstep
       if ( dble(istp) < field_duration ) then
          Mol%field_env(istp) = 0.5d0 - 0.5d0*dcos( dble(istp) * tau )
          fvect1(istp) = Mol%field_env(istp)*dcos( omega*dt*dble(istp) - phase )
       end if
    end do
    
    write(iout, '(A)') ' leaving subroutine cos_env'

    
  end subroutine cos_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE TRAP_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine trap_env
    
    !: envelope=2 : generates *Trapezoidal Evelope* pulse

    implicit none
    
    
    integer(8) :: istp
    

    write(iout,'(A)') ' in subroutine trap_env for trapezoidal envelope'    
    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0
    
    do istp=1, nstep
       if ( dble(istp) < field_duration ) then

          if ( dble(istp) < (ncyc-1)*period ) then
             !: slope up in first cycle
             if ( dble(istp) < int(period) ) then
                Mol%field_env(istp) = dble(istp)/period
                Mol%field_env(istp) = 1.d0
                fvect1(istp) = Mol%field_env(istp)*dsin( omega*dt*dble(istp) - phase )
             else
                !: flat/constant
                Mol%field_env(istp) = 1.d0
                fvect1(istp) = 1.d0*dsin( omega*dt*dble(istp) - phase )
             end if
          else 
             ! < C> slope down in last cycle
             Mol%field_env(istp) = ( field_duration - dble(istp) ) / period !: = (ncyc*period-istp)/period
             fvect1(istp) = Mol%field_env(istp)*dsin( omega*dt*dble(istp) - phase)
          end if
          
       end if
    end do
    

    write(iout, '(A)') ' leaving subroutine trap_env'
    
    
  end subroutine trap_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GAU_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine gau_env
    

    !: envelope=3 : generates *Gaussian Evelope* pulse

    
    implicit none
    

    real(8), parameter :: wdth = 16.0d0*log(2.0d0)
    real(8), parameter :: norm = 16.d0 / 15.0d0 
        
    integer(8) :: istp
    real(8)    :: frac
    

    write(iout,'(A)') ' in subroutine gau_env for Gaussian envelope'
    
    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0
    
    do istp=1, nstep
       if ( dble(istp) < field_duration ) then
          frac = dble(istp)/field_duration
          Mol%field_env(istp) = norm * ( dexp(-wdth*(frac-0.5d0 )**2) - 1.d0/16.d0 )
          fvect1(istp) = Mol%field_env(istp) * dsin( omega*dt*dble(istp) - phase )
       end if
    end do
    
    
    write(iout, '(A)') ' leaving subroutine gau_env'
    
    
  end subroutine gau_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE CW_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine cw_env

    !: envelope=4 : generates *Static* pulse
    
    implicit none    
    
    real(8), parameter :: frac = 2.d0/3.d0
    integer(8), parameter :: ramp_step = 16000
    integer(8) :: istp
    real(8) :: temp
    
    write(iout,'(A)') ' in subroutine cw_env for static pulse' 
    Mol%field_env = 0.d0
    fvect1 = 0.d0 

    do istp = 1, nstep
      if ( istp .le. int(frac*dble(ramp_step)*0.05d0/dt) ) then
        !: Suddently the nvidia compiler hates this line???
        !: NVFORTRAN-F-0000-Internal compiler error. Could not locate uplevel
        !:   instance for stblock 1214
        !fvect1(istp) = 1.d0 - ( 1.d0 - dble(istp)/(frac*dble(ramp_step)*0.05d0/dt) )**4
        temp = frac*dble(ramp_step)*0.05d0/dt
        fvect1(istp) = 1.d0 - ( 1.d0 - dble(istp)/temp)**4
      else
        fvect1(istp) = 1.d0
      end if
    end do

    Mol%field_env = fvect1
    
    write(iout, '(A)') ' leaving subroutine cw_env'

  end subroutine cw_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE CS_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine cs_env
    
    !: envelope=5 : generates *Bandrauk* pulse (linear rise then constant CW)


    implicit none

    
    integer(8) :: istp
    

    write(iout,'(A)') ' in subroutine cs_env for Bandrauk pulse (linear rise to constant)'
    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0
    
    do istp = 1, nstep
       if ( dble(istp) < field_duration ) then
          Mol%field_env(istp) = dble(istp) / field_duration
          fvect1(istp) = Mol%field_env(istp)*dcos( omega*dt*dble(istp) - phase )
       else
          Mol%field_env(istp) = 1.d0
          fvect1(istp) = 1.d0*dcos( omega*dt*dble(istp) - phase )
       end if
    end do

    
    write(iout, '(A)') ' leaving subroutine cs_env'
    
    
  end subroutine cs_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE RAMP_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine ramp_env

    !: generates a linear ramp over the entire duration of the pulse

    implicit none

    integer(8) :: istp
    real(8) :: denom

    write(iout, '(A)') ' in subroutine ramp_env for linear ramping field'
    write(iout, *) nstep
    Mol%field_env = 0.d0
    fvect1 = 0.d0
    denom = nstep-50
    write(iout, *) denom

    do istp = 1, nstep
      if ( istp < denom) then
        fvect1(istp) = dble(istp)/denom
        Mol%field_env(istp) = dble(istp)/denom
      else
        fvect1(istp) = 1.d0
        Mol%field_env(istp) = 1.d0
      end if
    end do

    write(iout, '(A)') ' leaving subroutine ramp_env'
  end subroutine ramp_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SIN2_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine sin2_env

    
    !: envelope=6 : generates *sin^2 squared* pulse ; to generate cos^2, set phase = pi/2

    
    implicit none

    
    integer(8) :: istp
    

    write(iout,'(A)') ' in subroutine sin2_env for sin squared pulse'

    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0

    do istp = 1,nstep
       if ( dble(istp) < field_duration ) then
          Mol%field_env(istp) = dsin( pi*dble(istp) / field_duration )**2
          fvect1(istp) = Mol%field_env(istp) * dsin( omega*dt*dble(istp)-phase )
       end if
    end do

    
    write(iout, '(A)') ' leaving subroutine sin2_env'
    
    
  end subroutine sin2_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE CIRC_ENV
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine circ_env(xpol,ypol)

    
    implicit none
    
    
    real(8), intent(in) :: xpol, ypol
    integer(8) :: istp
    

    if ( xpol.lt.0.d0 ) write(iout,'(A)') ' in subroutine circ_env for LEFT circularly/ellptclly polarized'
    if ( xpol.gt.0.d0 ) write(iout,'(A)') ' in subroutine circ_env for RIGHT circularly/ellptclly polarized'
    
    Mol%field_env = 0.d0 
    fvect1 = 0.d0 
    fvect2 = 0.d0

    do istp=1, nstep       
       if ( dble(istp) < field_duration ) then
          Mol%field_env(istp) = dsin( pi*dble(istp) / field_duration )**2
          fvect1(istp) = xpol * Mol%field_env(istp) * dsin( omega*dt*dble(istp) - phase )
          fvect2(istp) = ypol * Mol%field_env(istp) * dcos( omega*dt*dble(istp) - phase )
       end if
    end do
    
    write(iout, '(A)') ' leaving subroutine circ_env'
           
    
  end subroutine circ_env
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_GETFIELD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_getfield(myroutine,option)
    
    implicit none


    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option

    integer(8) :: i


    select case( trim(myroutine) )
    end select

    
  end subroutine writeme_getfield
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ERRORS_GETFIELD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_getfield(myroutine,option)


    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option

    
    select case( trim(myroutine) )
    !: subroutine check_enput
    case( 'check_emax' ) 
       select case(option)
       case('no_enput')  
          write(iout,'(A)') " WARNING: Setting default emax:  "
          write(iout,'(A)') " WARNING: *** NEMAX =  3 *** "
          go to 200
       case('no_lin') 
          write(iout,'(A)') ' set default field strengths for linearly polarized light'
          go to 200
       case('no_cir') 
          write(iout,'(A)') ' set default field strengths for circularly polarized light'
          go to 200
       end select

    end select
    
100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " IF WE WILL BE QUIET AND READY ENOUGH,"
    write(iout,'(A)') " WE WILL FIND COMPENSATION IN EVERY DISAPPOINTMENT"
    write(iout,'(A)') " - Henry David Thoreau"
    stop

200 continue

    
  end subroutine errors_getfield
end module getfield
