module readintegrals


  use variables_global
  implicit none

  
  integer(8) :: noa2, nva2, nob2, nvb2, noa3, nob3, nva3, nvb3, ntt
  real(8)    :: rskip
  character(1000) :: cskip
  character(1000), parameter :: form1e = "(1x,'allocated rank 1 array ',a8,' of size ',i11)"
  character(1000), parameter :: form2e = "(1x,'allocated rank 2 array ',a8,' of size ',i11,i11)"


contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_size_variables
    

    implicit none

    
    noa2 = noa*(noa+1)/2
    nva2 = nva*(nva+1)/2
    nob2 = nob*(nob+1)/2
    nvb2 = nvb*(nvb+1)/2
    noa3 = noa*(noa-1)/2
    nob3 = nob*(nob-1)/2
    nva3 = nva*(nva-1)/2
    nvb3 = nvb*(nvb-1)/2
    ntt  = nbasis*(nbasis+1)/2    


  end subroutine get_size_variables
!: paste below
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE Read_MatrixElements
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_matrixelements
    

    !: adapted from files in gauopen 
    implicit none
    
      Parameter (MaxAt=100,MaxBf=1000,MaxArr=100000000)
      logical        :: iamhere
      logical        :: EOF, AOInts, ASym, DAOInts
      integer(4)     :: LFN,IU,IVers,NLab,NBsuse,NE,Len12L,Len4L,IOpCl,ICGU
      integer(4)     :: NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot
      integer(4)     :: NAE,NBE,NI,NR,NRI,NTOT,LenBuf,n1,n2,n3,n4,n5,i,j,ij
      character(100) :: cskip
      character(80)  :: filename
      character(64)  :: LabFil, GVers, Title, CBuf
      integer(4)     :: maxat,maxbf,maxarr,ian,iattyp,ibfatm,ibftyp
      real(8)        :: atmchg,atmwgt,c
      real(8), allocatable :: RArr(:), density(:)
      Dimension IAn(MaxAt),IAtTyp(MaxAt),AtmChg(MaxAt),C(3,MaxAt), &
       IBfAtm(MaxBf),IBfTyp(MaxBf),AtmWgt(MaxAt)


 1000 Format(' File ',A,' IU=',I6)
 1010 Format(' Label ',A,' IVers=',I2,' NLab=',I2,' Version=',A,/,     &
        ' Title ',A60,/,' Mol%NAtoms=',I6,' NBasis=',I6,' NBsUse=',I6,       &
        ' ICharg=',I6,' Multip=',I6,' NE=',I6,' Len12L=',I1,' Len4L=', &
        I1,' IOpCl=',I6,' ICGU=',I3)
 1110 Format(' Label ',A48,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',I8, &
        ' LenBuf=',I8,' N=',5I6,' ASym=',L1)

      call write_header( 'read_matrixelements','initialize','enter' )
!   
!
!     Open the file.  IU is returned with the Fortran unit number
!     on which the file is open.
!
      !: read MatrixElements.dat
      Call Open_Read(trim(tdcidatfile),IU,LabFil,IVers,NLab,GVers,Mol%job_title, &
        Mol%natoms,nbasis,NBsUse,Mol%ICharg,Mol%Multip,NE,Len12L,Len4L,IOpCl,ICGU)
!      Write(IOut,1000) trim(filename), IU
      If(IU.le.0) Stop
      Write(IOut,1010) trim(LabFil),IVers,NLab,trim(GVers),           &
        trim(Mol%job_title),Mol%natoms,nbasis,NBsUse,Mol%ICharg,Mol%Multip,NE,Len12L, &
        Len4L,IOpCl,ICGU
!
!     Read the header records, which contain basic information about
!     the molecule such as number of atoms, atomic numbers, Cartesian
!     coordinates, and information about the basis functions.
!
      Call Rd_Head(IU,NLab,Mol%natoms,nbasis,IAn,IAtTyp,AtmChg,C,IBfAtm,  &
        IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,&
        NBTot)
 
      ntt = nbasis*(nbasis+1)/2 
      if(Mol%Multip.ge.0) then
        NAE = (NE+Mol%Multip-1)/2
        NBE = (NE-Mol%Multip+1)/2
      else
        NAE = (NE+Mol%Multip+1)/2
        NBE = (NE-Mol%Multip-1)/2
      end if
 
      Mol%NAE = NAE
      Mol%unrstedflag = IOpCl
      Mol%vabsflag = 1
 
      NOA = NAE - NFC
      NOB = NBE - NFC
      NROrb = NBsUse - NFC - NFV
      NVA = NRORB - NOA
      NVB = NRORB - NOB

!:      write(42,*) nae,nbe,noa,nob,nva,nvb,nrorb

    !: allocate coordinate arrays
    allocate( Mol%xcoord(Mol%natoms), Mol%ycoord(Mol%natoms), Mol%zcoord(Mol%natoms), Mol%myatom(Mol%natoms) )    
    do i = 1, Mol%natoms
      Mol%xcoord(i) = c(1,i)
      Mol%ycoord(i) = c(2,i)
      Mol%zcoord(i) = c(3,i)
      write(Mol%myatom(i),'(I3)') ian(i)
      !:write(42,*) ian(i),atmchg(i),xcoord(i),ycoord(i),zcoord(i)
    end do
 
      call get_size_variables
      allocate( RArr(MaxArr) )
   10 Call Rd_Labl(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,NRI,EOF)
      Write(IOut,1110) CBuf,NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym
!:      write(iout,*) '***',trim(CBUF),'***'
!:      flush(iout)
      If(NTot.gt.MaxArr) Write(IOut,*) " *** file too large *** increase MaxArr"
      If(.not.EOF .and. NTot.le.MaxArr) then
        select case ( trim(CBuf) )
        case ( trim('ALPHA ORBITAL ENERGIES') )
          allocate( Mol%orben(2*nrorb) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          do i = 1,NRorb
            Mol%orben(i) = RArr(NFC+i)
!:            write(42,*) i,Mol%orben(i)
            end do
        case ( trim('BETA ORBITAL ENERGIES') )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          do i = 1,NRorb
            Mol%orben(NRorb+i) = RArr(NFC+i)
!:            write(42,*) NRorb+i,Mol%orben(NRorb+i)
            end do
        case ( trim('TRANS MO COEFFICIENTS') )
          allocate( Mol%cmo_a(nrorb*nbasis) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          Mol%cmo_a = RArr(1:nrorb*nbasis)
          If(IOpCl.eq.1) then
            allocate( Mol%cmo_b(nrorb*nbasis) )
            Mol%cmo_b = RArr(nrorb*nbasis+1:2*nrorb*nbasis)
            end if
        case ( trim('ALPHA SCF DENSITY MATRIX') )
          allocate( density(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,density)
        case ( trim('BETA SCF DENSITY MATRIX') )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          If(IOpCl.eq.1) density = density + RArr(1:ntt)
        case ( trim('DIPOLE INTEGRALS') )
          allocate( Mol%dipxao(ntt), Mol%dipyao(ntt), Mol%dipzao(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          Mol%dipxao = RArr(1:ntt)
          Mol%dipyao = RArr(ntt+1:2*ntt)
          Mol%dipzao = RArr(2*ntt+1:3*ntt)
        case ( trim(' File   617'), trim('FILE 617 REALS'), trim('FILE 751 REALS') )
           allocate(Mol%socxao(nbasis,nbasis),Mol%socyao(nbasis,nbasis),Mol%soczao(nbasis,nbasis))
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          ij = 0
          do i = 1,nbasis
            do j = 1,i
              ij = ij + 1
              Mol%socxao(i,j) = dcmplx( 0.d0, RArr(ij) )
              Mol%socxao(j,i) = dconjg(Mol%socxao(i,j))
              Mol%socyao(i,j) = dcmplx( 0.d0, RArr(ij+ntt) )
              Mol%socyao(j,i) = dconjg(Mol%socyao(i,j))
              Mol%soczao(i,j) = dcmplx( 0.d0, RArr(ij+2*ntt) )
              Mol%soczao(j,i) = dconjg(Mol%soczao(i,j))
              end do
            end do
        case ( trim(' File   828'), trim('FILE 828 REALS'), trim('FILE 870 REALS') )
          allocate( Mol%vabsao(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,Mol%vabsao)
!:          write(iout,*) 'file 828', (vabsao(i),i=1,ntt)
        case ( trim(' File     1'), trim('FILE 1 REALS') )
          allocate( Mol%dijabAA(nva3,noa3) )
          Call Rd_RBf2(IU,NTot,nva3,noa3,LenBuf,Mol%dijabAA)
        case ( trim(' File     2'), trim('FILE 2 REALS') )
          allocate( Mol%dijabAB(nob*nvb,noa*nva) )
          Call Rd_RBf2(IU,NTot,nob*nvb,noa*nva,LenBuf,Mol%dijabAB)
!:          write(42,*) 'Bucket 2 head',(dijabAB(i,1),i=1,100)
!:          write(42,*) 'Bucket 2 tail',(dijabAB(nob*nvb-100+i,noa*nva),i=1,100)
        case ( trim(' File     3'), trim('FILE 3 REALS') )
          allocate( Mol%dijabBB(nvb3,nob3) )
          Call Rd_RBf2(IU,NTot,nvb3,nob3,LenBuf,Mol%dijabBB)
        case ( trim(' File     4'), trim('FILE 4 REALS') )
          allocate( Mol%dijklAA(noa3,noa3) )
          Call Rd_RBf3(IU,NTot,noa3,noa3,LenBuf,Mol%dijklAA)
        case ( trim(' File     5'), trim('FILE 5 REALS') )
          allocate( Mol%diajbAA(noa*nva,noa*nva) )
          Call Rd_RBf2(IU,NTot,noa*nva,noa*nva,LenBuf,Mol%diajbAA)
!:          write(42,*) 'Bucket 5 head',(diajbAA(i,1),i=1,100)
!:          write(42,*) 'Bucket 5 tail',(diajbAA(nob*nvb-100+i,noa*nva),i=1,100)
        case ( trim(' File     6'), trim('FILE 6 REALS') )
          allocate( Mol%dijklAB(nob2,noa2) )
          Call Rd_RBf2(IU,NTot,nob2,noa2,LenBuf,Mol%dijklAB)
        case ( trim(' File     7'), trim('FILE 7 REALS') )
          allocate( Mol%diajbAB(nvb2,noa2) )
          Call Rd_RBf2(IU,NTot,nvb2,noa2,LenBuf,Mol%diajbAB)
        case ( trim(' File     8'), trim('FILE 8 REALS') )
          allocate( Mol%diajbBA(nva2,nob2) )
          Call Rd_RBf2(IU,NTot,nva2,nob2,LenBuf,Mol%diajbBA)
        case ( trim(' File     9'), trim('FILE 9 REALS') )
          allocate( Mol%dijklBB(nob3,nob3) )
          Call Rd_RBf3(IU,NTot,nob3,nob3,LenBuf,Mol%dijklBB)
        case ( trim(' File    10'), trim('FILE 10 REALS') )
          allocate( Mol%diajbBB(nob*nvb,nob*nvb) )
          Call Rd_RBf2(IU,NTot,nob*nvb,nob*nvb,LenBuf,Mol%diajbBB)
!:          write(42,*) 'Bucket 10 head',(Mol%diajbBB(i,1),i=1,100)
!:          write(42,*) 'Bucket 10 tail',(Mol%diajbBB(nob*nvb-100+i,noa*nva),i=1,100)
        case ( trim(' File    11'), trim('FILE 11 REALS') )
          allocate( Mol%dijkaAA(noa*nva,noa3))
          Call Rd_RBf2(IU,NTot,noa*nva,noa3,LenBuf,Mol%dijkaAA)
        case ( trim(' File    12'), trim('FILE 12 REALS') )
          allocate( Mol%dijkaAB(nob*nvb,noa2))
          Call Rd_RBf2(IU,NTot,nob*nvb,noa2,LenBuf,Mol%dijkaAB)
        case ( trim(' File    13'), trim('FILE 13 REALS') )
          allocate( Mol%dijkaBA(noa*nva,nob2))
          Call Rd_RBf2(IU,NTot,noa*nva,nob2,LenBuf,Mol%dijkaBA)
        case ( trim(' File    14'), trim('FILE 14 REALS') )
          allocate( Mol%dijkaBB(nob*nvb,nob3))
          Call Rd_RBf2(IU,NTot,nob*nvb,nob3,LenBuf,Mol%dijkaBB)
        case ( trim(' File    15'), trim('FILE 15 REALS') )
          allocate( Mol%diabcAA(nva3,noa*nva))
          Call Rd_RBf2(IU,NTot,nva3,noa*nva,LenBuf,Mol%diabcAA)
        case ( trim(' File    16'), trim('FILE 16 REALS') )
          allocate( Mol%diabcAB(nvb2,noa*nva))
          Call Rd_RBf2(IU,NTot,nvb2,noa*nva,LenBuf,Mol%diabcAB)
        case ( trim(' File    17'), trim('FILE 17 REALS') )
          allocate( Mol%diabcBA(nva2,nob*nvb))
          Call Rd_RBf2(IU,NTot,nva2,nob*nvb,LenBuf,Mol%diabcBA)
        case ( trim(' File    18'), trim('FILE 18 REALS') )
          allocate( Mol%diabcBB(nvb3,nob*nvb))
          Call Rd_RBf2(IU,NTot,nvb3,nob*nvb,LenBuf,Mol%diabcBB)
        case ( trim(' File    19'), trim('FILE 19 REALS') )
          allocate( Mol%dabcdAA(nva3,nva3) )
          Call Rd_RBf3(IU,NTot,nva3,nva3,LenBuf,Mol%dabcdAA)
        case ( trim(' File    20'), trim('FILE 20 REALS') )
          allocate( Mol%dabcdAB(nvb2,nva2) )
          Call Rd_RBf2(IU,NTot,nvb2,nva2,LenBuf,Mol%dabcdAB)
        case ( trim(' File    21'), trim('FILE 21 REALS') )
          allocate( Mol%dabcdBB(nvb3,nvb3) )
          Call Rd_RBf3(IU,NTot,nvb3,nvb3,LenBuf,Mol%dabcdBB)
        case default
          write(iout,*) 'skipping ',CBUF
          Call Rd_Skip(IU,NTot,LenBuf)
        end select
        Goto 10
        endIf
        deallocate( RArr )
!:        deallocate(orben,cmo_a,socxao,socyao,soczao,vabsao)
!:        If(IOpCl.eq.1) deallocate(cmo_b)
!:        deallocate(dijabAB,diajbAA,diajbBB)
      Call Close_MatF(IU)
!
!    Calculate dipole moment
!
    Mol%dipx00 = 0.d0
    Mol%dipy00 = 0.d0
    Mol%dipz00 = 0.d0
    do i = 1,ntt
      Mol%dipx00 = Mol%dipx00 - 2*density(i)*Mol%dipxao(i)
      Mol%dipy00 = Mol%dipy00 - 2*density(i)*Mol%dipyao(i)
      Mol%dipz00 = Mol%dipz00 - 2*density(i)*Mol%dipzao(i)
    end do
    do i = 1,nbasis
      ij = i*(i+1)/2
      Mol%dipx00 = Mol%dipx00 + density(ij)*Mol%dipxao(ij)
      Mol%dipy00 = Mol%dipy00 + density(ij)*Mol%dipyao(ij)
      Mol%dipz00 = Mol%dipz00 + density(ij)*Mol%dipzao(ij)
    end do
    do i = 1,Mol%natoms
      Mol%dipx00 = Mol%dipx00 + AtmChg(i)*Mol%xcoord(i)
      Mol%dipy00 = Mol%dipy00 + AtmChg(i)*Mol%ycoord(i)
      Mol%dipz00 = Mol%dipz00 + AtmChg(i)*Mol%zcoord(i)
    end do
    write(42,*) " dipole ",Mol%dipx00,Mol%dipy00,Mol%dipz00
    deallocate( density )
!:    deallocate(dipxao,dipyao,dipzao)
!:    deallocate(xcoord,ycoord,zcoord,myatom)
            
    call write_header( 'read_matrixelements','initialize','leave' )
    
    
  end subroutine read_matrixelements

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_vabs_ao( Qstore )
    

    implicit none    
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables    
    

    if ( Qstore ) then

       allocate( Mol%vabsao(ntt) )
       write(iout,form1e) 'Mol%vabsao', ntt  ;   call track_mem( ntt )
       
       read(10, '(A)') cskip
       read(10,*) ( Mol%vabsao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    

  end subroutine read_vabs_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipx_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables

    
    if ( Qstore ) then

       allocate( Mol%dipxao(ntt) )
       write(iout,form1e) 'Mol%dipxao', ntt  ;  call track_mem( ntt )       
       
       read(10,'(A)') cskip
       read(10,*) ( Mol%dipxao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    

  end subroutine read_dipx_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipy_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables


    if ( Qstore ) then

       allocate( Mol%dipyao(ntt) )
       write(iout,form1e) 'Mol%dipyao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       read(10,*) ( Mol%dipyao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else

       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if

    
  end subroutine read_dipy_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipz_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables

    
    if ( Qstore ) then

       allocate( Mol%dipzao(ntt) )
       write(iout,form1e) 'Mol%dipzao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       read(10,*) ( Mol%dipzao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else

       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    

  end subroutine read_dipz_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socx_ao( Qstore )
    

    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables


    if ( Qstore ) then

       !allocate( Mol%socxao(ntt), r_array(ntt) )
       allocate( Mol%socxao(nbasis,nbasis) ) 
       write(iout,form1e) 'Mol%socxao', ntt  ;  call track_mem( ntt )       

       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             Mol%socxao(i,j) = dcmplx( 0.d0, rdum )
             Mol%socxao(j,i) = dconjg( Mol%socxao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !Mol%socxao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
       !deallocate( r_array )

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    
    
  end subroutine read_socx_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socy_ao( Qstore ) 


    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables
    

    if ( Qstore ) then

       !allocate( Mol%socyao(ntt), r_array(ntt) )
       allocate( Mol%socyao(nbasis,nbasis) )
       write(iout,form1e) 'Mol%socyao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             Mol%socyao(i,j) = dcmplx( 0.d0, rdum )
             Mol%socyao(j,i) = dconjg( Mol%socyao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !Mol%socyao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
       !deallocate( r_array )

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    

  end subroutine read_socy_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socz_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables


    if ( Qstore ) then

       !allocate( Mol%soczao(ntt), r_array(ntt) )
       allocate( Mol%soczao(nbasis,nbasis) )
       write(iout,form1e) 'Mol%soczao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             Mol%soczao(i,j) = dcmplx( 0.d0, rdum )
             Mol%soczao(j,i) = dconjg( Mol%soczao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !Mol%soczao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

       !deallocate( r_array ) 

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if

    
  end subroutine read_socz_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_orben( Qstore )
    
    
    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i

    
    call get_size_variables

    
    if ( Qstore ) then

       allocate( Mol%orben(norb) )
       write(iout,form1e) 'Mol%orben', norb  ;  call track_mem( norb )
       
       read(10,'(A)') cskip
       read(10,*) ( Mol%orben(i) , i=1, norb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else
       
       read(10,'(A)') cskip
       read(10,*) ( rskip , i=1, norb )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    
    
  end subroutine read_orben
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_cmo( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables


    if ( Qstore ) then

       allocate( Mol%cmo_a(nrorb*nbasis) )
       write(iout,form1e) 'Mol%cmo_a', nrorb*nbasis  ;  call track_mem( nrorb*nbasis )
       
       if ( unrestricted ) then
          allocate( Mol%cmo_b(nrorb*nbasis) ) 
          write(iout,form1e) 'Mol%cmo_b', nrorb*nbasis  ;  call track_mem( nrorb*nbasis )
       end if
       
       read(10,'(A)') cskip
       read(10,*) ( Mol%cmo_a(i) , i=1, nbasis*nrorb )       
       if (unrestricted) read(10,*) ( Mol%cmo_b(i) , i=1, nbasis*nrorb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else
       
       read(10,'(A)') cskip
       read(10,*) ( rskip , i=1, nbasis*nrorb )       
       if (unrestricted) read(10,*) ( rskip, i=1, nbasis*nrorb )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    
  end subroutine read_cmo
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_custom_ham

    implicit none
    integer(8) :: i, j


    !: write out warning
    write(iout,'(A)') ' directly reading in the Hamiltonian from TDCI.dat '
    write(iout,'(A)') ' WARNING:  change subroutine read_custom_ham if necessary'
    write(iout,'(A)') ' '
    flush(iout)

    
    !: change if not TDA elements read from Gaussian output
    cis_vec(1) = 0.d0 


    !: read from TDCI.dat
    read(10,'(A)') cskip
    do i=2, nstates
       read(10,*) ( cis_vec( (i-1)*nstates+j ), j=2,nstates )
    end do


    !: make sure Ham is Hermetian
    do i=1, nstates
       do j=(i+1), nstates
          cis_vec( (i-1)*nstates+j ) = 0.5d0 * ( cis_vec((i-1)*nstates+j) + cis_vec((j-1)*nstates+i) )
          cis_vec( (j-1)*nstates+i ) = cis_vec( (i-1)*nstates+j )
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_custom_ham
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_1
    !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%dijabAA(nva3,noa3) )            ; Mol%dijabAA = 0.d0 ! Bucket 2    
    write(iout,form2e) 'Mol%dijabAA', nva3, noa3  ; call track_mem( nva3*noa3 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijabAA(j,i), j=1, nva3 ), i=1, noa3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_1
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_2
    !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%dijabAB(nobnvb,noanva) )            ;  Mol%dijabAB = 0.d0 ! Bucket 2    
    write(iout,form2e) 'Mol%dijabAB', nobnvb, noanva  ;  call track_mem( nobnvb*noanva )
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijabAB(j,i), j=1, nobnvb ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
!:          write(42,*) 'Bucket 2 head',(Mol%dijabAB(i,1),i=1,100)
!:          write(42,*) 'Bucket 2 tail',(Mol%dijabAB(nob*nvb-100+i,noa*nva),i=1,100)
    
    
  end subroutine read_bucket_2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_3
    !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB
    
    implicit none
    integer(8) :: i, j 
    

    call get_size_variables

    
    allocate( Mol%dijabBB(nvb3,nob3) )            ;  Mol%dijabBB = 0.d0 ! Bucket 2    
    write(iout,form2e) 'Mol%dijabBB', nvb3, nob3  ;  call track_mem( nvb3*nob3 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijabBB(j,i), j=1, nvb3 ), i=1, nob3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_3
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_4
    !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L,IJ.LE.KL  IJKL
    

    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%dijklAA(noa3,noa3) )            ;  Mol%dijklAA = 0.d0 
    write(iout,form2e) 'Mol%dijklAA', noa3, noa3  ;  call track_mem( noa3*noa3 )
    

    read(10,'(A)') cskip
    do i=1, noa3
       do j=i, noa3
          read(10,*) Mol%dijklAA(j,i)
          Mol%dijklAA(i,j) = Mol%dijklAA(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_4
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_5
    !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%diajbAA(noanva,noanva) )            ;  Mol%diajbAA = 0.d0 ! Bucket 5 
    write(iout,form2e) 'Mol%diajbAA', noanva, noanva  ;  call track_mem( noanva*noanva )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diajbAA(j,i), j=1, noanva ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
!:          write(42,*) 'Bucket 5 head',(Mol%diajbAA(i,1),i=1,100)
!:          write(42,*) 'Bucket 5 tail',(Mol%diajbAA(nob*nvb-100+i,noa*nva),i=1,100)

    
  end subroutine read_bucket_5
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_6
    !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%dijklAB(nob2,noa2) )            ;  Mol%dijklAB = 0.d0 
    write(iout,form2e) 'Mol%dijklAB', nob2, noa2  ;  call track_mem( nob2*noa2 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijklAB(j,i), j=1, nob2 ), i=1, noa2 )
    write(iout,'(A)') ' finishd reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_6
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_7
    !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( Mol%diajbAB(nvb2,noa2) )            ;  Mol%diajbAB = 0.d0 
    write(iout,form2e) 'Mol%diajbAB', nvb2, noa2  ;  call track_mem( nvb*noa2 )


    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diajbAB(j,i), j=1, nvb2 ), i=1, noa2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_7
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_8
    !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( Mol%diajbBA(nva2,nob2) )            ;  Mol%diajbBA = 0.d0 
    write(iout,form2e) 'Mol%diajbBA', nva2, nob2  ;  call track_mem( nva2*nob2 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diajbBA(j,i), j=1, nva2 ), i=1, nob2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_8
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_9
    !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( Mol%dijklBB(nob3,nob3) )            ;  Mol%dijklBB = 0.d0 
    write(iout,form2e) 'Mol%dijklBB', nob3, nob3  ;  call track_mem( nob3*nob3 )


    read(10,'(A)') cskip
    do i=1, nob3
       do j=i, nob3
          read(10,*) Mol%dijklBB(j,i)
          Mol%dijklBB(i,j) = Mol%dijklBB(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_9
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!  
  subroutine read_bucket_10
    !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( Mol%diajbBB(nobnvb,nobnvb) )            ;  Mol%diajbBB = 0.d0 ! Bucket 10
    write(iout,form2e) 'Mol%diajbBB', nobnvb, nobnvb  ;  call track_mem( nobnvb*nobnvb ) 

    
    if(unrestricted) then
       read(10,'(A)') cskip
       read(10,*) ( ( Mol%diajbBB(j,i), j=1, nobnvb ), i=1, nobnvb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    end if
!:          write(42,*) 'Bucket 10 head',(Mol%diajbBB(i,1),i=1,100)
!:          write(42,*) 'Bucket 10 tail',(Mol%diajbBB(nob*nvb-100+i,noa*nva),i=1,100)

    
  end subroutine read_bucket_10
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_11
    !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( Mol%dijkaAA(noanva,noa3))             ;  Mol%dijkaAA = 0.d0 
    write(iout,form2e) 'Mol%dijkaAA', noa3, noanva  ;  call track_mem( noa3*noanva )


    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijkaAA(j,i), j=1, noanva ), i=1, noa3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_11
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_12
    !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( Mol%dijkaAB(nobnvb,noa2))             ;  Mol%dijkaAB = 0.d0 
    write(iout,form2e) 'Mol%dijkaAB', nobnvb, noa2  ;  call track_mem( noa2*nobnvb )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijkaAB(j,i), j=1, nobnvb ), i=1, noa2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_12
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_13
    !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%dijkaBA(noanva,nob2))             ;  Mol%dijkaBA = 0.d0 
    write(iout,form2e) 'Mol%dijkaBA', noanva, nob2  ;  call track_mem( nob2*noanva )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijkaBA(j,i), j=1, noanva ), i=1, nob2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_13
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_14
    !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( Mol%dijkaBB(nobnvb,nob3))             ;  Mol%dijkaBB = 0.d0 
    write(iout,form2e) 'Mol%dijkaBB', nobnvb, nob3  ;  call track_mem( nob3*nobnvb )


    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dijkaBB(j,i), j=1, nobnvb ), i=1, nob3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_14
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_15
    !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( Mol%diabcAA(nva3,noanva))             ;  Mol%diabcAA = 0.d0 
    write(iout,form2e) 'Mol%diabcAA', noanva, nva3  ;  call track_mem( noanva*nva3 )
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diabcAA(j,i), j=1, nva3 ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_15
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_16
    !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( Mol%diabcAB(nvb2,noanva))             ;  Mol%diabcAB = 0.d0
    write(iout,form2e) 'Mol%diabcAB', nvb2, noanva  ;  call track_mem( noanva*nvb2 ) 
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diabcAB(j,i), j=1, nvb2 ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_16
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_17
    !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables    

    
    allocate( Mol%diabcBA(nva2,nobnvb))             ;  Mol%diabcBA = 0.d0 
    write(iout,form2e) 'Mol%diabcBA', nva2, nobnvb  ;  call track_mem( nobnvb*nva2 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diabcBA(j,i), j=1, nva2 ), i=1, nobnvb )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_17
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_18
    !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables
    

    allocate( Mol%diabcBB(nvb3,nobnvb))              ;  Mol%diabcBB = 0.d0 
    write(iout,form2e)  'Mol%diabcBB', nvb3, nobnvb  ;  call track_mem( nobnvb*nvb3 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%diabcBB(j,i), j=1, nvb3 ), i=1, nobnvb )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_18
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_19
    !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables


    allocate( Mol%dabcdAA(nva3,nva3) )            ;  Mol%dabcdAA = 0.d0 
    write(iout,form2e) 'Mol%dabcdAA', nva3, nva3  ;  call track_mem( nva3*nva3 )
    

    read(10,'(A)') cskip
    do i=1, nva3
       do j=i, nva3
          read(10,*) Mol%dabcdAA(j,i)
          Mol%dabcdAA(i,j) = Mol%dabcdAA(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_19
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_20
    !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD

    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( Mol%dabcdAB(nvb2,nva2) )            ;  Mol%dabcdAB = 0.d0 
    write(iout,form2e) 'Mol%dabcdAB', nvb2, nva2  ;  call track_mem( nva2*nvb2 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( Mol%dabcdAB(j,i), j=1, nvb2 ), i=1, nva2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_20
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_21
    !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( Mol%dabcdBB(nvb3,nvb3) )            ;  Mol%dabcdBB = 0.d0
    write(iout,form2e) 'Mol%dabcdBB', nvb3, nvb3  ;  call track_mem( nvb3*nvb3 )
    
    
    read(10,'(A)') cskip
    do i=1, nvb3
       do j=i, nvb3
          read(10,*) Mol%dabcdBB(j,i)
          Mol%dabcdBB(i,j) = Mol%dabcdBB(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_21
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  !: GET INTEGRALS FOR EACH BUCKET.  LOOK get last function INTEGER(8) GET_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabAA( i, j, a, b )
    !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB

    integer(8), intent(in) :: i, j, a, b
    integer(8) :: ij, ab

    !: <ij||ab> DijabAA( ab, ij )

    ij = get_index( i,j, noa, 'lt' )
    ab = get_index( a,b, nva, 'lt' )

    get_dijabAA = Mol%dijabAA(ab,ij)

  end function get_dijabAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabBB( i, j, a, b )
    !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB

    integer(8), intent(in) :: i, j, a, b
    integer(8) :: ij, ab

    !: <IJ||AB> DijabBB( AB, IJ )

    IJ = get_index( I, J, nob, 'lt' )
    AB = get_index( A, B, nvb, 'lt' )
    
    get_dijabBB = Mol%dijabBB(AB,IJ)
    
  end function get_dijabBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabAB( i, j, a, b )
    !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <iJ||aB> DijabAB( JB, ia )

    ia = get_index( i,a, nva, 'all' )
    jb = get_index( j,b, nvb, 'all' )

    get_dijabAB = Mol%dijabAB(jb,ia)

  end function get_dijabAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklAA( i, j, k, l )
    !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L  IJKL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: kl, ij
    real(8) :: sign

    !: <ij||kl> DijklAA( k<l, i<j )

    sign = 1.d0
    if ( l.lt.k ) sign = - 1.d0 * sign
    if ( j.lt.i ) sign = - 1.d0 * sign

    kl = get_index( k,l, noa, 'lt' )
    ij = get_index( i,j, noa, 'lt' )

    get_dijklAA = sign * Mol%dijklAA(kl,ij)

  end function get_dijklAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbAA( i, a, j, b )
    !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <ia||jb> DiajbAA( jb, ia )

    ia = get_index( i,a, nva, 'all' )
    jb = get_index( j,b, nva, 'all' )

    get_diajbAA = Mol%diajbAA(jb,ia)

  end function get_diajbAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklAB( i, j, k, l )
    !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: jl, ik

    if (.not.allocated(Mol%dijklAB)) then
      write(iout,*) "Mol%dijklAB is NOT allocated" ; flush(iout)
    end if

    !: <iJ||kL> DijklAB( J<=L,i<=k )

    jl = get_index( j, l, nob, 'le' )
    ik = get_index( i, k, noa, 'le' )

    get_dijklAB = Mol%dijklAB(jl,ik)

  end function get_dijklAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbAB( i, a, j, b )
    !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ab, ij

    !: <iA||jB> DiajbAB( A<=B, i<=j )

    ab = get_index( a,b, nvb, 'le' )
    ij = get_index( i,j, noa, 'le' )

    get_diajbAB = Mol%diajbAB(ab,ij)

  end function get_diajbAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbBA( i, a, j, b )
    !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ab, ij

    !: <Ia||Jb> DiajbBA( a<=b, I<=J )

    ab = get_index( a,b, nva, 'le' )
    ij = get_index( i,j, nob, 'le' )

    get_diajbBA = Mol%diajbBA(ab,ij)

  end function get_diajbBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklBB( i, j, k, l )
    !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: kl, ij
    real(8) :: sign

    if (.not.allocated(Mol%dijklBB)) then
      write(iout,*) "Mol%dijklBB is NOT allocated" ; flush(iout)
    end if

    !: <IJ||KL> DijklBB( K<L,I<J )

    sign = 1.d0
    if ( l.lt.k ) sign = -1.d0 * sign
    if ( j.lt.i ) sign = -1.d0 * sign

    kl = get_index( k,l, nob, 'lt' )
    ij = get_index( i,j, nob, 'lt' )

    get_dijklBB = sign * Mol%dijklBB(kl,ij)

  end function get_dijklBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbBB( i, a, j, b )
    !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <IA||JB> DiajbBB( JB, IA )

    ia = get_index( i,a, nvb, 'all' )
    jb = get_index( j,b, nvb, 'all' )

    get_diajbBB = Mol%diajbBB(jb,ia)

  end function get_diajbBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaAA( i, j, k, a )
    !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ka, ij
    real(8) :: sign

    !: <ij||ka> DijkaAA( ka, i<j )

    sign = 1.d0
    if ( j.lt.i ) sign = -1.d0 * sign

    ka = get_index( k, a, nva, 'all' )
    ij = get_index( i, j, noa, 'lt' )

    get_dijkaAA = sign * Mol%dijkaAA(ka,ij)

  end function get_dijkaAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaAB( i, j, k, a )
    !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ja, ik

    !: <iJ||kA> DijkaAB( JA, i<=k )

    ja = get_index( j,a, nvb, 'all' )
    ik = get_index( i,k, noa, 'le' )

    get_dijkaAB = Mol%dijkaAB(ja,ik)

  end function get_dijkaAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaBA( i, j, k, a )
    !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ja, ik

    !: <Ij||Ka> DijkaBA( ja, I<=K )

    ja = get_index( j, a, nva, 'all' )
    ik = get_index( i, k, nob, 'le' )

    get_dijkaBA = Mol%dijkaBA(ja,ik)

  end function get_dijkaBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaBB( i, j, k, a )
    !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ka, ij
    real(8) :: sign

    !: <IJ||KA> DijkaBB( KA, I<J )

    sign = 1.d0
    if ( j.lt.i ) sign = -1.d0  * sign

    ka = get_index( k, a, nvb, 'all' )
    ij = get_index( i, j, nob, 'lt' )

    get_dijkaBB = sign * Mol%dijkaBB(ka,ij)

  end function get_dijkaBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcAA( i, a, b, c )
    !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: bc, ia
    real(8) :: sign

    !: <ia||bc> DiabcAA( b<c, ia )

    sign = 1.d0
    if( c.lt.b ) sign = -1.d0 * sign

    bc = get_index( b, c, nva, 'lt'  )
    ia = get_index( i, a, nva, 'all' )

    get_diabcAA = sign * Mol%diabcAA(bc,ia)

  end function get_diabcAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcAB( i, a, b, c )
    !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: ib, ac

    !: <iA||bC> DiabcAB( A<=C, ib)

    ac = get_index( a, c, nvb, 'le'  )
    ib = get_index( i, b, nva, 'all' )

    get_diabcAB = Mol%diabcAB(ac,ib)

  end function get_diabcAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcBA( i, a, b, c )
    !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: ac, ib

    !: <Ia||Bc> DiabcBA( a<=c, IB )

    ac = get_index( a, c, nva, 'le' )
    ib = get_index( i, b, nvb, 'all' )

    get_diabcBA = Mol%diabcBA(ac,ib)

  end function get_diabcBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcBB( i, a, b, c )
    !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: bc, ia
    real(8) :: sign

    !: <IA||BC> DiabcBB( B<C, IA )

    sign = 1.d0
    if ( c.lt.b ) sign = -1.d0 * sign

    bc = get_index( b, c, nvb, 'lt'  )
    ia = get_index( i, a, nvb, 'all' )

    get_diabcBB = sign * Mol%diabcBB(bc,ia)

  end function get_diabcBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdAA( a, b, c, d )
    !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: ab, cd
    real(8) :: sign

    !: <ab||cd> DabcdAA( c<d, a<b )

    sign = 1.d0
    if ( d.lt.c ) sign = -1.d0 * sign
    if ( b.lt.a ) sign = -1.d0 * sign

    ab = get_index( a, b, nva, 'lt' )
    cd = get_index( c, d, nva, 'lt' )

    get_dabcdAA = sign * Mol%dabcdAA(cd,ab)

  end function get_dabcdAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdAB( a, b, c, d )
    !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: ac, bd

    !: <aB||cD> DabcdAB( B<=D, a<=c )

    bd = get_index( b, d, nvb, 'le' )
    ac = get_index( a, c, nva, 'le' )

    get_dabcdAB = Mol%dabcdAB(bd,ac)

  end function get_dabcdAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdBB( a, b, c, d )
    !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: cd, ab
    real(8) :: sign

    !: <AB||CD> DabcdBB( C<D, A<B )

    sign = 1.d0
    if ( d.lt.c ) sign = -1.d0 * sign
    if ( b.lt.a ) sign = -1.d0 * sign

    cd = get_index( c, d, nvb, 'lt' )
    ab = get_index( a, b, nvb, 'lt' )

    get_dabcdBB = sign * Mol%dabcdBB(cd,ab)

  end function get_dabcdBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  integer(8) function get_index( i, j, nn, opt )

    integer(8), intent(in) :: i,j, nn
    character(*), intent(in) :: opt

    integer(8) :: big, small

    if ( trim(opt).eq.'lt' ) then
       big   = max( i,j )
       small = min( i,j )
       get_index = (small-1) * nn - small*(small-1)/2 + (big-small)
       return
    else if ( trim(opt).eq.'le' ) then
       big   = max( i,j )
       small = min( i,j )
       get_index = (small-1) * nn - (small-1)*(small-2)/2 + (big-small) + 1
       return
    else if ( trim(opt).eq.'all' ) then
       get_index = (i-1) * nn + j
       return
    end if

  end function get_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module readintegrals
