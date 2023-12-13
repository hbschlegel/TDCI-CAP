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
      integer(4)     :: NAE,NBE,i
      character(100) :: cskip
      character(80)  :: filename
      character(64)  :: LabFil, GVers, Title, CBuf
      integer(4)     :: maxat,maxbf,maxarr,ian,iattyp,ibfatm,ibftyp
      real(8)        :: atmchg,atmwgt,c
      real(8), allocatable :: RArr(:), density(:)
      real(8), allocatable :: xcoord(:),ycoord(:),zcoord(:),myatom(:),orben(:)
      real(8), allocatable :: cmo_a(:),cmo_b(:),dipxao(:),dipyao(:),dipzao(:)
      real(8), allocatable :: socxao(:),socyao(:),soczao(:)
      real(8), allocatable :: vabsao(:),dijabaa(:,:),dijabab(:,:),dijabbb(:,:)
      real(8), allocatable :: diajbaa(:,:),diajbab(:,:),diajbbb(:,:),diajbba(:,:)
      real(8), allocatable :: dijklaa(:,:),dijklab(:,:),dijklbb(:,:)
      real(8), allocatable :: dijkaaa(:,:),dijkaab(:,:),dijkaba(:,:),dijkabb(:,:)
      real(8), allocatable :: diabcaa(:,:),diabcab(:,:),diabcba(:,:),diabcbb(:,:)
      real(8), allocatable :: dabcdaa(:,:),dabcdab(:,:),dabcdbb(:,:)
      Dimension IAn(MaxAt),IAtTyp(MaxAt),AtmChg(MaxAt),C(3,MaxAt), &
       IBfAtm(MaxBf),IBfTyp(MaxBf),AtmWgt(MaxAt)

 1000 Format(' File ',A,' IU=',I6)
 1010 Format(' Label ',A,' IVers=',I2,' NLab=',I2,' Version=',A,/,     &
        ' Title ',A,/,' NAtoms=',I6,' NBasis=',I6,' NBsUse=',I6,       &
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
      filename = 'MatrixElements.dat' 
      Call Open_Read(trim(filename),IU,LabFil,IVers,NLab,GVers,job_title, &
        natoms,nbasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU)
      Write(IOut,1000) trim(filename), IU
      If(IU.le.0) Stop
      Write(IOut,1010) trim(LabFil),IVers,NLab,trim(GVers),       &
        trim(Title),natoms,nbasis,NBsUse,ICharg,Multip,NE,Len12L, &
        Len4L,IOpCl,ICGU
!
!     Read the header records, which contain basic information about
!     the molecule such as number of atoms, atomic numbers, Cartesian
!     coordinates, and information about the basis functions.
!
      Call Rd_Head(IU,NLab,natoms,nbasis,IAn,IAtTyp,AtmChg,C,IBfAtm,  &
        IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,&
        NBTot)

      ntt = nbasis*(nbasis+1)/2 
      If(Multip.ge.0) then
        NAE = (NE+Multip-1)/2
        NBE = (NE-Multip+1)/2
      else
        NAE = (NE+Multip+1)/2
        NBE = (NE-Multip-1)/2
        endIf
 
      unrstedflag = IOpCl
      vabsflag = 1
 
      NOA = NAE - NFC
      NOB = NBE - NFC
      NROrb = NBsUse - NFC - NFV
      NVA = NRORB - NOA
      NVB = NRORB - NOB

    !: allocate coordinate arrays
    allocate( xcoord(natoms), ycoord(natoms), zcoord(natoms), myatom(natoms) )    
    do i = 1, natoms
         xcoord(i) = c(1,i)
         ycoord(i) = c(2,i)
         zcoord(i) = c(3,i)
         myatom(i) = ian(i)
         end do
 
      call get_size_variables
      allocate( RArr(MaxArr) )
   10 Call Rd_Labl(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,NRI,EOF)
      Write(IOut,1110) CBuf,NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym
      If(NTot.gt.MaxArr) Write(IOut,*) " *** file too large *** increase MaxArr"
      If(.not.EOF .and. NTot.le.MaxArr) then
        select case ( trim(CBuf) )
        case ( trim('ALPHA ORBITAL ENERGIES') )
          allocate( orben(nrorb) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          do i = 1,NRorb
            orben(i) = RArr(NFC+i)
            end do
        case ( trim('BETA ORBITAL ENERGIES') )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          do i = 1,NRorb
            orben(NRorb+i) = RArr(NFC+i)
            end do
        case ( trim('TRANS MO COEFFICIENTS') )
          allocate( cmo_a(nrorb*nbasis) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          cmo_a = RArr(1:nrorb*nbasis)
          If(IOpCl.eq.1) then
            allocate( cmo_b(nrorb*nbasis) )
            cmo_b = RArr(nrorb*nbasis+1:2*nrorb*nbasis)
            end if
        case ( trim('ALPHA SCF DENSITY MATRIX') )
          allocate( density(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,density)
        case ( trim('BETA SCF DENSITY MATRIX') )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          density = density + RArr
        case ( trim('DIPOLE INTEGRALS') )
          allocate( dipxao(ntt), dipyao(ntt), dipzao(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          dipxao = RArr(1:ntt)
          dipyao = RArr(ntt+1:2*ntt)
          dipzao = RArr(2*ntt+1:3*ntt)
        case ( trim('617') )
          allocate( socxao(ntt), socyao(ntt), soczao(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          socxao = RArr(1:ntt)
          socyao = RArr(ntt+1:2*ntt)
          soczao = RArr(2*ntt+1:3*ntt)
        case ( trim('828') )
          allocate( vabsao(ntt) )
          Call Rd_RBuf(IU,NTot,LenBuf,vabsao)
        case ( trim('1') )
          allocate( dijabAA(nva3,noa3) )
          Call Rd_RBf2(IU,NTot,nva3,noa3,LenBuf,dijabAA)
        case ( trim('2') )
          allocate( dijabAB(nobnvb,noanva) )
          Call Rd_RBf2(IU,NTot,nobnvb,noanva,LenBuf,dijabAB)
        case ( trim('3') )
          allocate( dijabBB(nvb3,nob3) )
          Call Rd_RBf2(IU,NTot,nvb3,nob3,LenBuf,dijabBB)
        case ( trim('4') )
          allocate( dijklAA(noa3,noa3) )
          Call Rd_RBf3(IU,NTot,noa3,noa3,LenBuf,dijklAA)
        case ( trim('5') )
          allocate( diajbAA(noanva,noanva) )
          Call Rd_RBf2(IU,NTot,noanva,noanva,LenBuf,diajbAA)
        case ( trim('6') )
          allocate( dijklAB(nob2,noa2) )
          Call Rd_RBf2(IU,NTot,nob2,noa2,LenBuf,dijklAB)
        case ( trim('7') )
          allocate( diajbAB(nvb2,noa2) )
          Call Rd_RBf2(IU,NTot,nvb2,noa2,LenBuf,diajbAB)
        case ( trim('8') )
          allocate( diajbBA(nva2,nob2) )
          Call Rd_RBf2(IU,NTot,nva2,nob2,LenBuf,diajbBA)
        case ( trim('9') )
          allocate( dijklBB(nob3,nob3) )
          Call Rd_RBf3(IU,NTot,nob3,nob3,LenBuf,dijklBB)
        case ( trim('10') )
          allocate( diajbBB(nobnvb,nobnvb) )
          Call Rd_RBf2(IU,NTot,nobnvb,nobnvb,LenBuf,diajbBB)
        case ( trim('11') )
          allocate( dijkaAA(noanva,noa3))
          Call Rd_RBf2(IU,NTot,noanva,noa3,LenBuf,dijkaAA)
        case ( trim('12') )
          allocate( dijkaAB(nobnvb,noa2))
          Call Rd_RBf2(IU,NTot,nobnvb,noa2,LenBuf,dijkaAB)
        case ( trim('13') )
          allocate( dijkaBA(noanva,nob2))
          Call Rd_RBf2(IU,NTot,noanva,nob2,LenBuf,dijkaBA)
        case ( trim('14') )
          allocate( dijkaBB(nobnvb,nob3))
          Call Rd_RBf2(IU,NTot,nobnvb,nob3,LenBuf,dijkaBB)
        case ( trim('15') )
          allocate( diabcAA(nva3,noanva))
          Call Rd_RBf2(IU,NTot,nva3,noanva,LenBuf,diabcAA)
        case ( trim('16') )
          allocate( diabcAB(nvb2,noanva))
          Call Rd_RBf2(IU,NTot,nvb2,noanva,LenBuf,diabcAB)
        case ( trim('17') )
          allocate( diabcBA(nva2,nobnvb))
          Call Rd_RBf2(IU,NTot,nva2,nobnvb,LenBuf,diabcBA)
        case ( trim('18') )
          allocate( diabcBB(nvb3,nobnvb))
          Call Rd_RBf2(IU,NTot,nvb3,nobnvb,LenBuf,diabcBB)
        case ( trim('19') )
          allocate( dabcdAA(nva3,nva3) )
          Call Rd_RBf3(IU,NTot,nva3,nva3,LenBuf,dabcdAA)
        case ( trim('20') )
          allocate( dabcdAB(nvb2,nva2) )
          Call Rd_RBf3(IU,NTot,nvb2,nva2,LenBuf,dabcdAB)
        case ( trim('21') )
          allocate( dabcdBB(nvb3,nvb3) )
          Call Rd_RBf3(IU,NTot,nvb3,nvb3,LenBuf,dabcdBB)
        case default
          Call Rd_Skip(IU,NTot,LenBuf)
        end select
        Goto 10
        endIf
        deallocate( RArr )
      Call Close_MatF(IU)
!
!    Calculate dipole moment
!
     dipx00 = 0.d0
     dipy00 = 0.d0
     dipz00 = 0.d0
     do i = 1,ntt
       dipx00 = dipx00 - 2*density(i)*dipxao(i)
       dipy00 = dipy00 - 2*density(i)*dipyao(i)
       dipz00 = dipz00 - 2*density(i)*dipzao(i)
       end do
     do i = 1,nbasis
       ij = i*(i+1)/2
       dipx00 = dipx00 + density(i)*dipxao(i)
       dipy00 = dipy00 + density(i)*dipyao(i)
       dipz00 = dipz00 + density(i)*dipzao(i)
       end do
     do i = 1,natoms
       dipx00 = dipx00 + AtmChg(i)*xcoord(i)
       dipy00 = dipy00 + AtmChg(i)*ycoord(i)
       dipz00 = dipz00 + AtmChg(i)*zcoord(i)
       end do
    deallocate( density )
            
    call write_header( 'read_matrixelements','initialize','leave' )
    
    
  end subroutine read_matrixelements
