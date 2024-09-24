module variables_setup

  !: module variables_setup
  !: these are not required for the propagation step
  
  implicit none


  type MolInfo

    !: read in from TDCI.dat
    integer(4) :: unrstedflag, vabsflag
    integer(4) :: ICharg, Multip, natoms
    integer(8) :: NAE
    real(8),      allocatable :: xcoord(:), ycoord(:), zcoord(:)
    character(3), allocatable :: myatom(:)
    character(4000) :: job_title

    !: envelope for field
    real(8), allocatable :: field_env(:) 

    !: <0|A|0> elements
    real(8) :: dipx00, dipy00, dipz00, vabs00

    real(8), allocatable :: &
         orben(:),  & !: MO orbital energies
         vabsao(:), & !: Vabs matrix elements in AO basis
         dipxao(:), & !: x dipole matrix elements in AO basis
         dipyao(:), & !: y dipole matrix elements in AO basis
         dipzao(:), & !: z dipole matrix elements in AO basis
         cmo_a(:),  & !: MO-LCAO coefficients for alpha 
         cmo_b(:),  & !: MO-LCAO coefficients for beta       
         overlap(:)   !: AO overlap matrix
 

    complex(8), allocatable :: &
         socxao(:,:), socyao(:,:), soczao(:,:) !: move to variables_setup module later

    
    real(8), allocatable :: &
         dijabAA(:,:), & !: <ij||ab> AAAA BUCKET 1
         dijabAB(:,:), & !: <ij||ab> ABAB BUCKET 2
         dijabBB(:,:), & !: <ij||ab> BBBB BUCKET 3 
         dijklAA(:,:), & !: <ij||kl> AAAA BUCKET 4
         diajbAA(:,:), & !: <ia||jb> AAAA BUCKET 5
         dijklAB(:,:), & !: <ij||kl> ABAB BUCKET 6
         diajbAB(:,:), & !: <ia||jb> ABAB BUCKET 7
         diajbBA(:,:), & !: <ia||jb> BABA BUCKET 8
         dijklBB(:,:), & !: <ij||kl> BBBB BUCKET 9
         diajbBB(:,:), & !: <ia||jb> BBBB BUCKET 10
         dijkaAA(:,:), & !: <ij||ka> AAAA BUCKET 11       
         dijkaAB(:,:), & !: <ij||ka> ABAB BUCKET 12
         dijkaBA(:,:), & !: <ij||ka> BABA BUCKET 13
         dijkaBB(:,:), & !: <ij||ka> BBBB BUCKET 14
         diabcAA(:,:), & !: <ia||bc> AAAA BUCKET 15 
         diabcAB(:,:), & !: <ia||bc> ABAB BUCKET 16
         diabcBA(:,:), & !: <ia||bc> BABA BUCKET 17
         diabcBB(:,:), & !: <ia||bc> BBBB BUCKET 18
         dabcdAA(:,:), & !: <ab||cd> AAAA BUCKET 19
         dabcdAB(:,:), & !: <ab||cd> ABAB BUCKET 20
         dabcdBB(:,:)    !: <ab||cd> BBBB BUCKET 21


    real(8), allocatable :: &
         dipxmoa(:,:) , & !: x dipole matrix elements in MO basis for alpha
         dipymoa(:,:) , & !: y dipole matrix elements in MO basis for alpha
         dipzmoa(:,:) , & !: z dipole matrix elements in MO basis for alpha
         vabsmoa(:,:) , & !: Vabs matrix elements in MO basis for alpha  
         dipxmob(:,:) , & !: x dipole matrix elements in MO basis for beta
         dipymob(:,:) , & !: y dipole matrix elements in MO basis for beta
         dipzmob(:,:) , & !: z dipole matrix elements in MO basis for beta
         vabsmob(:,:)     !: Vabs matrix elements in MO basis for alpha
    

    complex(8), allocatable :: & 
         socmoAA(:,:), socmoBB(:,:), socmoAB(:,:) 


  end type MolInfo


         

    
  
end module variables_setup
  
