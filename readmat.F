*Deck ReadMat
      Program ReadMat
      Implicit Real*8(A-H,O-Z)
C
C     This example program reads a matrix element file.  It illustrates
C     use of the Fortran interface.  The name of the file is the first
C     argument on the command line.
C
      Parameter (LStr=64,IOut=6,MaxAt=100,MaxBf=1000,MaxArr=100000000,
     $  MaxAr5=(MaxArr/5))
      Character*(LStr) LabFil, GVers, Title, CBuf, FileName*80
      Logical EOF, AOInts, ASym, DAOInts
      Dimension IAn(MaxAt),IAtTyp(MaxAt),AtmChg(MaxAt),C(3,MaxAt),
     $  IBfAtm(MaxBf),IBfTyp(MaxBf),AtmWgt(MaxAt),IArr(MaxArr),
     $  RArr(MaxArr),RArr3(3,MaxAr5),IArr5(5,MaxAr5)
      Complex*16 CArr(MaxArr)
      Equivalence (IArr(1),RArr(1),CArr(1),RArr3(1,1))
 1000 Format(' File ',A,' IU=',I6)
 1010 Format(' Label ',A,' IVers=',I2,' NLab=',I2,' Version=',A,/,
     $  ' Title ',A,/,' NAtoms=',I6,' NBasis=',I6,' NBsUse=',I6,
     $  ' ICharg=',I6,' Multip=',I6,' NE=',I6,' Len12L=',I1,' Len4L=',
     $  I1,' IOpCl=',I6,' ICGU=',I3)
 1030 Format('    IAn=',20I4)
 1040 Format(' IAtTyp=',5I12)
 1050 Format(' AtmChg=',5F12.6)
 1060 Format('      C=',3F12.6)
 1070 Format(' IBfAtm=',10I8)
 1080 Format(' IBfTyp=',10I8)
 1090 Format(' AtmWgt=',5F12.6)
 1100 Format(' NFC=',I6,' NFV=',I6,' ITran=',I6,' IDum9=',I6,/,
     $  ' NShlAO=',I6,' NPrmAO=',I7,' NShlDB=',I6,' NPrmDB=',I7,
     $  ' NBTot=',I7)
 1110 Format(' Label ',A48,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',I8,
     $  ' LenBuf=',I8,' N=',5I6,' ASym=',L1)
 1120 Format(' Length=',I12,' exceeds MaxArr=',I12,'.')
 1140 Format(' LR=',I12,' but NTot=',I12,'.')
 1150 Format('    IArr=',10I12)
 1160 Format('    RArr=',10F12.6)
 1170 Format('    CArr=',10F12.6)
 1180 Format(5I9,3D15.6)
C
      FileName = ' '
      Call Get_Command_Argument(1,FileName,LFN,IStat)
      If(IStat.ne.0) then
        Write(6,*) 'Error fetching argument, status=',IStat
        Stop
        endIf
C
C     Open the file.  IU is returned with the Fortran unit number
C     on which the file is open.
C
      Call Open_Read(FileName(1:LFN),IU,LabFil,IVers,NLab,GVers,Title,
     $  NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU)
      Write(IOut,1000) FileName(1:LFN), IU
      If(IU.le.0) Stop
      LLab = Len_Trim(LabFil)
      LVers = Len_Trim(GVers)
      LTit = Len_Trim(Title)
      Write(IOut,1010) LabFil(1:LLab),IVers,NLab,GVers(1:LVers),
     $  Title(1:LTit),NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,
     $  Len4L,IOpCl,ICGU
C
C     Read the header records, which contain basic information about
C     the molecule such as number of atoms, atomic numbers, Cartesian
C     coordinates, and information about the basis functions.
C
      Call Rd_Head(IU,NLab,NAtoms,NBasis,IAn,IAtTyp,AtmChg,C,IBfAtm,
     $  IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,
     $  NBTot)
      Write(IOut,1030) (IAn(I),I=1,NAtoms)
      Write(IOut,1040) (IAtTyp(I),I=1,NAtoms)
      Write(IOut,1050) (AtmChg(I),I=1,NAtoms)
      Write(IOut,1060) ((C(J,I),J=1,3),I=1,NAtoms)
      Write(IOut,1070) (IBfAtm(I),I=1,NBasis)
      Write(IOut,1080) (IBfTyp(I),I=1,NBasis)
      Write(IOut,1090) (AtmWgt(I),I=1,NAtoms)
      Write(IOut,1100) NFC, NFV, ITran, IDum9, NShlAO, NPrmAO, NShlDB,
     $  NPrmDB, NBTot
C
C     Loop over objects in the file, printing each.  Any object we
C     do not recognize is skipped using Rd_Skip.
C
   10 Call Rd_Labl(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,
     $  NRI,EOF)
      Write(IOut,1110) CBuf,NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym
      If(.not.EOF) then
        LR = LenArr(N1,N2,N3,N4,N5)
        If((NI*NR).eq.0.and.LR.ne.NTot) then
          Write(IOut,1140) LR, NTot
          Call Rd_Skip(IU,NTot,LenBuf)
        else if(LR.gt.MaxArr) then
          Write(IOut,1120) LR, MaxArr
          Call Rd_Skip(IU,NTot,LenBuf)
        else if(NI.eq.1.and.NR.eq.0) then
          Call Rd_IBuf(IU,NTot,LenBuf,IArr)
          Write(IOut,1150) (IArr(I),I=1,LR)
        else if(NI.eq.0.and.NR.eq.1.and.NRI.eq.1) then
          Call Rd_RBuf(IU,NTot,LenBuf,RArr)
          Write(IOut,1160) (RArr(I),I=1,LR)
        else if(NI.eq.0.and.NR.eq.1.and.NRI.eq.2) then
          Call Rd_CBuf(IU,NTot,LenBuf,CArr)
          Write(IOut,1170) (CArr(I),I=1,LR)
        else if(AOInts(CBuf)) then
          Call Rd_2EN(IU,NR,LR,NR*LR,NTot,LenBuf,RArr)
          Write(IOut,1160) (RArr(I),I=1,NR*LR)
        else if(DAOInts(CBuf)) then
          Call Rd_SpA(IU,NI,NR,NTot,LenBuf,IArr5,RArr3)
          Write(IOut,1180)
     $      ((IArr5(J,I),J=1,5),(RArr3(J,I),J=1,3),I=1,NTot)
        else if(NI.eq.1) then
          If(NRI.eq.1) then
            Call Rd_RInd(IU,NR,LR,NTot,LenBuf,LNZ,RArr)
            Write(IOut,1160) (RArr(I),I=1,NR*LNZ)
          else
            Call Rd_RInd(IU,NRI*NR,LR,NTot,LenBuf,LNZ,CArr)
            Write(IOut,1170) (CArr(I),I=1,NR*LNZ)
            endIf
        else
          Call Rd_Skip(IU,NTot,LenBuf)
          endIf
        Goto 10
        endIf
      Call Close_MatF(IU)
      End
