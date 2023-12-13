module units 

  !: module units
  !: contains units conversion parameters and functions


  implicit none  


  complex(8), parameter :: eye = dcmplx(0.d0,1.d0)


  real(8), parameter :: pi      = 2.d0*dacos(0.d0) !: pi
  real(8), parameter :: deg2rad = pi / 180.d0     !: multiply for degrees --> radians
  real(8), parameter :: rad2deg = 1.d0/deg2rad    !: multiply for radians --> degrees
  

  !: all constants and conversion factors from NIST
  
  
  !: fundamental physical constants
  real(8), parameter :: planck   = 6.626070040d-34  !: Planck's constant in Js
  real(8), parameter :: hbar     = planck/2.d0/pi   !: Planck's constant divided by 2pi
  real(8), parameter :: light    = 299792458d0      !: speed of light in m/s
  real(8), parameter :: avogrado = 6.022140857d23   !: Avogradro's number 1/mole
  real(8), parameter :: kb       = 1.38064852d-23   !: Boltzmann's constant in J/K
  
  
  !: energy conversion
  real(8), parameter :: wvnbr2joules = 1.d0/planck/light !: multiply for cm-1 --> J
  real(8), parameter :: eV2joules    = 1.6021766208d-19  !: multiply for eV   --> J
  real(8), parameter :: au2joules    = 4.359744650d-18   !: multiply for E_h  --> J
  real(8), parameter :: au2eV        = 27.21138602d0     !: multiply for E_h  --> eV

  real(8), parameter :: joules2eV    = 1.d0/eV2joules    !: multiply for J --> eV
  real(8), parameter :: joules2wvnbr = 1.d0/wvnbr2joules !: multiply for J --> cm-1
  real(8), parameter :: joules2au    = 1.d0/au2joules    !: multiply for J --> E_h  
  real(8), parameter :: eV2au        = 1.d0/au2eV        !: multiply for E_h --> eV

  
  !: length conversion
  real(8), parameter :: ang2m  = 1.d-10            !: multiply for Angstrom  --> m
  real(8), parameter :: nm2m   = 1.d-9             !: multiply for nanometer --> m
  real(8), parameter :: bohr2m = 0.52917721067d-10 !: multiply for bohr      --> m
  
  real(8), parameter :: m2ang  = 1.d0/ang2m        !: multiply for m --> Angstrom
  real(8), parameter :: m2nm   = 1.d0/nm2m         !: multiply for m --> nanometer
  real(8), parameter :: m2bohr = 1.d0/bohr2m       !: multiply for m --> bohr

  real(8), parameter :: bohr2ang = bohr2m * m2ang  !: multiply for bohr     --> Angstrom
  real(8), parameter :: ang2bohr = 1.d0/bohr2ang   !: multiply for Angstrom --> bohr 


  !: time conversion
  real(8), parameter :: autime2s = 2.418884326509d-17  !: multiply for au_time --> s (hbar/E_h)
  real(8), parameter :: ps2s     = 1.d-12              !: multiply for picosecond --> s
  real(8), parameter :: fs2s     = 1.d-15              !: multiply for femtosecond --> s
  real(8), parameter :: as2s     = 1.d-18              !: multiply for attosecond --> s

  real(8), parameter :: s2autime = 1.d0/autime2s !: multiply for s --> au_time
  real(8), parameter :: s2ps     = 1.d0/ps2s     !: multiply for s --> picosecond
  real(8), parameter :: s2fs     = 1.d0/fs2s     !: multiply for s --> femtosecond 
  real(8), parameter :: s2as     = 1.d0/as2s     !: multiply for s --> attosecond 
  real(8), parameter :: au2fs    = autime2s*s2fs !: multiply for au_time to fs


  !: dipole conversion
  real(8), parameter :: audip2sidip = 8.478353552d-30 !: multiply for au_dipole --> Cm (e * bohr)
  real(8), parameter :: debye2sidip = 3.335641d-30    !: multiply for debye     --> Cm

  real(8), parameter :: sidip2audip = 1.d0/audip2sidip !: multiply for Cm --> au_dipole
  real(8), parameter :: sidip2debye = 1.d0/debye2sidip !: multiply for Cm --> debye
  
  real(8), parameter :: audip2debye = audip2sidip * sidip2debye !: multiply for au_dipole --> debye
  real(8), parameter :: debye2audip = 1.d0/audip2debye          !: multiply for debye     --> au_dipole


  !: field strenght conversion
  real(8), parameter :: aufield2si1 = 5.142206707d11 !: multiply for au_field --> V/m ( E_h/ (e*bohr) )
  real(8), parameter :: aufield2si2 = 514.2206707d0  !: multiply for au_field --> V/nm
  
  real(8), parameter :: sifield_au1 = 1.d0/aufield2si1 !: multiply for V/m  --> au_field
  real(8), parameter :: sifield_au2 = 1.d0/aufield2si2 !: multiply for V/nm --> au_field
  
  
  !: velocity conversion
  real(8), parameter :: auvel2si = 2.18769126277d6 !: multiply for au_velocity --> m/s ( bohr * E_h / hbar )
  real(8), parameter :: sivel2au = 1.d0/auvel2si   !: multiply for m/s         --> au_velocity
  
  
  !: momentum conversion
  real(8), parameter :: aumom2si = 1.992851882d-24 !: multiply for au_momentum --> kg m /s ( hbar / bohr )
  real(8), parameter :: simom2au = 1.d0/aumom2si   !: multiply for SI momentum --> au_momentum
  

  !: energy flux conversion
  real(8), parameter :: auflux2si = 3.50944758d16    !: multiply for au    --> W/cm2 
  real(8), parameter :: siflux2au = 1.d0/auflux2si   !: multiply for W/cm2 --> au


  
  !: computer stuff
  integer(8), parameter :: &
       byte2bit  = 8 ,    & !: 1 byte = 8 bits
       mega2byte = 1.d6 , & !: 1 megabyte = 10^6 byte
       giga2byte = 1.d9     !: 1 gigabyte = 10^9 byte
  

contains
  !: ------------------------------- :!
  !: FUNCTION CONVERT_FREQ
  !: ------------------------------- :!
  real function convert_freq(myfreq) 

    
    !: converst freq in nm to au time


    implicit none
    real(8), intent(inout) :: myfreq
    
    
    myfreq = 2.d0 * pi * light / ( myfreq * nm2m ) * autime2s    
    convert_freq = myfreq

        
  end function convert_freq
  !: ------------------------------- :!
  !: FUNCTION CONVERT_TIME
  !: ------------------------------- :!
  real function convert_time(mytime,inunits,outunits)

    
    !: converts time in some seconds to au time


    implicit none
    real(8), intent(inout) :: mytime
    character(*), intent(in) :: inunits, outunits    

    
    if( trim(outunits).eq.'au' ) then
       select case( trim(inunits) )
       case( 'ps' ) 
          mytime = mytime * ps2s * s2autime
       case( 'fs' )
          mytime = mytime * fs2s * s2autime
       case( 'as' )
          mytime = mytime * as2s * s2autime
       case default
          mytime = -1.d20
       end select
    end if


    convert_time = mytime

    
  end function convert_time  
  !: ------------------------------- :!
  !: FUNCTION CONVERT_ANGLES
  !: ------------------------------- :!
  real function convert_angles(myangle) 


    !: converts angle from deg to radians


    implicit none
    real(8), intent(inout)  :: myangle
    

    myangle = myangle * deg2rad    
    convert_angles = myangle


  end function convert_angles
  ! ****************************************************************! 
  real function convert_energy(myenergy,inunits,outunits)

    
    !: converts energy in ev or cm-1 to hartrees


    implicit none
    real(8), intent(inout) :: myenergy
    character(*), intent(in) :: inunits, outunits
  

    if( trim(outunits).eq.'au' ) then
       select case( trim(inunits) )
       case( 'ev' )
          myenergy = myenergy * ev2au
       case( 'cm-1' )
          myenergy = myenergy * wvnbr2joules * joules2au
       case default
          myenergy = 0.d0          
       end select
    end if

    
    convert_energy = myenergy

    
  end function convert_energy
  ! ****************************************************************! 
end module units
