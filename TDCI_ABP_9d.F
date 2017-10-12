C
C FORTRAN code to propagate a molecule in an electric field
C using Absorbing Boundary Potential to calculate ionization rates
C
C Initial version written by J.Sonk, Sept 2012
C Extensions/modifications by P.Krause, H.B.Schlegel
C Restricted and unrestricted CIS, parallel implementation by H.B.Schlegel
C--------------------------------------------------------------------------------
C Program requirements:
C A CIS calultation of the molecule with GDV using link l914 
C Latest version of l914 at grid.wayne.edu:/home/groups/hbsq/share/subst*
C which creates the file TDCI.dat (contains orbital energies, MO
C coefficients, absorbing potential matix elements (AO), dipole matrix
C elements (AO), (ij||ab) (MO) and (ia||jb) (MO)
C Also needed is the set-up file named 'input'
C--------------------------------------------------------------------------------
C Compile:
C gfortran TDCI_ABP.F -L. -llapack -lrefblas -ltmglib -o tdci_abp
C
C or
C
C pgf95 TDCI_ABP-dipoles-new8c.F  -Minfo -Mneginfo -time -fast -O3 -Miomutex -Mconcur=allcores 
C -mp=allcores -Munroll -Mvect  -llapack -lblas -o tdcis-test8c>compile.log
C
C Execute: ./tdci_abp 
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Bugs:
C
C********************************************************************************
      PROGRAM MainIon
C********************************************************************************
      CHARACTER*80 PATH,FILENAME,EFile,TDCIdata
      CHARACTER*160 FULLNAME
      INTEGER::NStates,MaxStep,Ncyc,Ecnt,N,ct,envelope,IOut,iwrite,
     $ NStUse,outstep,Lscratch,Lscratch1,NEmax,Ndir,DirForm,INFO,I,J,IJ
      INTEGER::NBasis,NRorb,NOA,NVA,NOB,NVB,NOANVA,NOBNVB,NTT,len,Ndata
      COMPLEX*16,ALLOCATABLE,DIMENSION(:)::Psi,Psi1,Psi2,CIdata
      REAL*8,ALLOCATABLE,DIMENSION(:)::OrbEn,CMOA,CMOB,
     $ VabsAO,DipXAO,DipYAO,DipZAO,
     $ scratch,scratch1,scratch2,scratch3,scratch4,
     $ VabsMOA,DipXMOA,DipYMOA,DipZMOA,VabsMOB,DipXMOB,DipYMOB,DipZMOB,
     $ Dijab,DiajbA,DiajbB,Hmat,CISvec,CISeig,EmaxVec,
     $ DirxVec,DiryVec,DirzVec,ABPmat,TDmat,TDX,TDY,TDZ,
     $ ExpABP,EigABP,evect,fvect1,fvect2,Pop,PopInt,TDCIresults
      REAL*8 Pi,emax,delta,omega,period,phase,tt,tf,EVec(18),ScanDir,
     $ aufs,dipx,dipy,dipz,dirx,diry,dirz,norm,ABPScal,
     $ ExpectA,Vabs00,Fac,Theta,Phi,NormD,tmp,tmpx,tmpy,tmpz
     $ Heuristic,EigMax,Ionization,Zero,One,Two
      REAL start,finish,start1,finish1
      LOGICAL OutExists,efExists,UnRestricted,Linear
      INTEGER OMP_get_num_threads,OMP_get_num_procs,OMP_get_thread_num,
     $ NThread1,NThread2,NThread3,NThread4,NThreads
      COMMON IOut
      Zero = 0.D0
      One = 1.D0
      Two = 2.D0
      Pi = Two*ACOS(Zero)
      aufs = 0.024188843
      IOut=42
      Linear=.True.
 1000 FORMAT(F25.18)
 1001 FORMAT(2F25.18)
 2000 FORMAT(83F25.18)
 3000 FORMAT(3F25.18)
 4000 FORMAT(3(F20.16,'+',F20.16,'i') )
 5000 FORMAT('Total Elapse time: ', F12.3,' seconds.')
 6000 FORMAT(7I10)
 6001 FORMAT(D22.12)
 1233 FORMAT(2X,A21,F21.11)
 1234 FORMAT(2X,A21,F12.4)
 1235 FORMAT(2X,A21,F12.4,1X,A2)
 1236 FORMAT(2X,A21,2X,'(',F5.3,',',F5.3,',',F5.3,')')
 1237 FORMAT(2X,A23,2X,'(',F7.3,',',F7.3,',',F7.3,')')
 1238 FORMAT(" Unrestricted CIS, NOA=",I5," NVA=",I5,
     $    " NOB=",I5," NVB=",I5)
 1239 FORMAT("Restricted CIS, NOA=",I5," NVA=",I5)
 1240 Format(" NStates = ",I8," NStUse = ",I8," Max(CISeig) = ",F10.4)
 1241 Format(" Heuristic lifetime = ",F8.4," Sqrt(CISeig-Ionization)")
 2345 FORMAT(2X,A21,I7)
 
      NAMELIST /Field/ omega, ncyc, envelope, phase, NEmax, NDir,
     $ DirForm
      NAMELISt /System/ EigMax, Ionization, maxstep, delta, outstep,
     $ iwrite, dipx, dipy, dipz, Heuristic
C
C     !==========================================================!
C     ! Creating log file and debug directory                    !
C     !==========================================================!
C
      OPEN(IOut,FILE='tdci-abp.log')
c      IOut = 6
c     INQUIRE( FILE="./debug/.", exist=OutExists )
c     IF (OutExists .eqv. .FALSE. ) then
c         call system('mkdir debug')
c     END IF
c     PATH='debug'
      CALL DNT(start)
C
C     !==========================================================!
C     ! Checking for the setup file 'input' and reading in info   !
C     !==========================================================!
C
C     Check to see if 'input' is present 
      INQUIRE(FILE='input', EXIST=efExists)
      IF ( efExists .eqv. .FALSE.) THEN
         STOP 'Setup file is not present'
      END IF
C
C     Open and read 'input'      
      OPEN(10,FILE='input')
      READ(10,NML=Field)
      READ(10,NML=System)
      
      Ndata = NEmax*Ndir
      If(EigMax.eq.Zero) EigMax = 1.D40
      If(envelope.gt.6) Linear=.False.
      ALLOCATE(EmaxVec(NEmax))
      ALLOCATE(DirxVec(NDir))
      ALLOCATE(DiryVec(NDir))
      ALLOCATE(DirzVec(NDir))
      ALLOCATE(TDCIresults(27*Ndata))
      Do I = 1,NEmax
        Read(10,*) EmaxVec(I)
        end Do
C
C     Read polarization directions
C     if NDir > 0, read (x,y,z) vector and normalize
C     if NDir < 0, read polar angles theta and phi
C
      Fac = Pi/180.D0
      If(DirForm.gt.Zero) then
        Do I = 1,NDir
          Read(10,*) DirxVec(I),DiryVec(I),DirzVec(I)
          NormD = Sqrt(DirxVec(I)**2+DiryVec(I)**2+DirzVec(I)**2)
          DirxVec(I) = DirxVec(I)/NormD
          DiryVec(I) = DiryVec(I)/NormD
          DirzVec(I) = DirzVec(I)/NormD
          Theta = ACos(DirzVec(I))/Fac
          Phi = ATan2(DiryVec(I),DirxVec(I))/Fac
          Do J = 1,NEmax
            IJ = I+(J-1)*Ndir
            TDCIresults(IJ) = EmaxVec(J)
            TDCIresults(IJ+Ndata) = Theta
            TDCIresults(IJ+2*Ndata) = Phi
            TDCIresults(IJ+3*Ndata) = DirxVec(I)
            TDCIresults(IJ+4*Ndata) = DiryVec(I)
            TDCIresults(IJ+5*Ndata) = DirzVec(I)
            end Do
          end Do
      else
        Do I = 1,NDir
          Read(10,*) Theta,Phi
C          Write(6,*) I,Theta,Phi
          DirxVec(I) = Sin(Fac*Theta)*Cos(Fac*Phi)
          DiryVec(I) = Sin(Fac*Theta)*Sin(Fac*Phi)
          DirzVec(I) = Cos(Fac*Theta)
C          Write(6,*) DirxVec(I),DiryVec(I),DirzVec(I)
          Do J = 1,NEmax
            IJ = I+(J-1)*Ndir
            TDCIresults(IJ) = EmaxVec(J)
            TDCIresults(IJ+Ndata) = Theta
            TDCIresults(IJ+2*Ndata) = Phi
            TDCIresults(IJ+3*Ndata) = DirxVec(I)
            TDCIresults(IJ+4*Ndata) = DiryVec(I)
            TDCIresults(IJ+5*Ndata) = DirzVec(I)
            end Do
          end Do
      endIf
      write(iout,*) " finsihed reading input"
C
C     for circularly polarized pulses, set up two vectors
C     perpendicular to the propagation direction secified by the input
C
      If(.not.Linear) then
       Do I = 1,NDir
        tmp=Sqrt(DirxVec(I)**2+DiryVec(I)**2)
        If(tmp.gt.1.D-8) then
         tmpx = -DiryVec(I)/tmp
         tmpy = DirxVec(I)/tmp
         tmpz = Zero
        else
         tmpx = One
         tmpy = Zero
         tmpz = Zero
         endIf
        Do J = 1,NEmax
         IJ = I+(J-1)*Ndir
         TDCIresults(IJ+6*Ndata) = TDCIresults(IJ)
         TDCIresults(IJ+7*Ndata) = ACos(tmpz)/Fac
         TDCIresults(IJ+8*Ndata) = ATan2(tmpy,tmpx)/Fac
         TDCIresults(IJ+9*Ndata) = tmpx
         TDCIresults(IJ+10*Ndata) = tmpy
         TDCIresults(IJ+11*Ndata) = tmpz
         TDCIresults(IJ+12*Ndata) = TDCIresults(IJ)
         TDCIresults(IJ+15*Ndata) = DiryVec(I)*tmpz-DirzVec(I)*tmpy
         TDCIresults(IJ+16*Ndata) = DirzVec(I)*tmpx-DirxVec(I)*tmpz
         TDCIresults(IJ+17*Ndata) = DirxVec(I)*tmpy-DiryVec(I)*tmpx
         TDCIresults(IJ+13*Ndata) = ACos(TDCIresults(IJ+17*Ndata))/Fac
         TDCIresults(IJ+14*Ndata) =
     $    ATan2(TDCIresults(IJ+16*Ndata),TDCIresults(IJ+15*Ndata))/Fac
         TDCIresults(IJ+18*Ndata) = Zero
         end Do
        end Do
        endIf
      ScanDir = 0.0D0
C
C     Calculate period, phase, and pulse duration, and final propa time     
C
      ABPScal = 1.0D0
      period  = 2.0D0*Pi/(omega*delta)
C      phase   = 0.0D0*Pi/2.0D0
      tt      =(ncyc*period)
      tf      =(delta*maxstep)
C
C     If Maxstep is less than 1, calculate 1.5* the pulse duration needed
C     for the requested number of cycles with a given period padded by 2      
C
      IF (maxstep .le. 1) THEN
         maxstep = Int(1.5D0 * ncyc * period + 2)
         WRITE(IOut,*)'Setting maxstep to: ',maxstep
      END IF
C
C Write simulation information to log file
C
      Write(IOut,*) "Log file for tdci-abp"
      Write(IOut,*) "TDCIS-ABP: all expecation values. PERMANENT DIPOLE"
      Write(IOut,*) "_______________________________________________"
      Write(IOut,*) "Propagation setup:"
      Write(IOut,2345) " MaxStep =",maxstep
      Write(IOut,1234) " delta t =",delta
      Write(IOut,1235) " Final time =",(delta*maxstep),"au"
      Write(IOut,1235) " Final time =",(delta*maxstep*aufs),"fs"
      WRITE(IOut,*) " Wavefunction analyzed every ",outstep," steps"
      WRITE(IOut,*) " Write flag for CIcoeff,dipole,ion,pop,norm.dat ",
     $ iwrite
      Write(IOut,*) "_______________________________________________"
      Write(IOut,*) "Pulse setup:"
      Write(IOut,1234) " Emax =",emax
      Write(IOut,1236) " Polarized (x,y,z) =",dirx,diry,dirz
      Write(IOut,1234) " Omega =",omega
      Write(IOut,2345) " Nr of Cycles =",ncyc
      Write(IOut,2345) " Envelope =",envelope
      Write(IOut,1234) " Phase =",phase
      Write(IOut,*) "_______________________________________________"
      Write(IOut,1237) "Molecular Dipole moment:", dipx,dipy,dipz
C
      ALLOCATE(evect(maxstep))
      ALLOCATE(fvect1(maxstep))
      ALLOCATE(fvect2(maxstep))
C
C     Calling Field generating Subroutines
C
      IF(envelope.eq.1) THEN
        WRITE(IOut,*)'                    --> cosine pulse'      
        CALL CosEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.2) THEN
        WRITE(IOut,*)'                    --> trapezoidal pulse'
        CALL TrapEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.3) THEN
        WRITE(IOut,*)'                    --> gaussian pulse'
        CALL GauEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.4) THEN
        WRITE(IOut,*)'                    --> static pulse'
        CALL CwEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.0) THEN
        WRITE(IOut,*)'                    --> "Null Pulse"'
        CALL NullEnGen(maxstep,evect,fvect1)
      ELSE IF (envelope.eq.5) THEN
        WRITE(IOut,*)'                    --> "Bandrauk" pulse'
        CALL CsEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.6) THEN
        WRITE(IOut,*)'                    --> sine squared pulse'
        CALL Cos2EnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $ evect,fvect1)
      ELSE IF (envelope.eq.7) THEN
        WRITE(IOut,*)'                  --> right circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect1,fvect2,One,One)
      ELSE IF (envelope.eq.8) THEN
        WRITE(IOut,*)'                  --> left circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect1,fvect2,-One,One)
      ELSE IF (envelope.eq.9) THEN
        WRITE(IOut,*)'                  --> testx circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect1,fvect2,One,Zero)
        Do I = 1,NDir
          Do J = 1,NEmax
            IJ = I+(J-1)*Ndir
            TDCIresults(IJ+9*Ndata)=DirxVec(I)
            TDCIresults(IJ+10*Ndata)=DiryVec(I)
            TDCIresults(IJ+11*Ndata)=DirzVec(I)
            End Do
          End Do
      ELSE IF (envelope.eq.10) THEN
        WRITE(IOut,*)'                  --> testy circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect2,fvect1,One,Zero)
        Do I = 1,NDir
          Do J = 1,NEmax
            IJ = I+(J-1)*Ndir
            TDCIresults(IJ+15*Ndata)=DirxVec(I)
            TDCIresults(IJ+16*Ndata)=DiryVec(I)
            TDCIresults(IJ+17*Ndata)=DirzVec(I)
            End Do
          End Do
      ELSE IF (envelope.eq.11) THEN
        WRITE(IOut,*)'                --> testxyR circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect2,fvect1,One,One)
        Do I = 1,NDir
          Do J = 1,NEmax
            IJ = I+(J-1)*Ndir
            tmpx=TDCIresults(IJ+9*Ndata)
            tmpy=TDCIresults(IJ+10*Ndata)
            tmpz=TDCIresults(IJ+11*Ndata)
            TDCIresults(IJ+9*Ndata)=TDCIresults(IJ+15*Ndata)
            TDCIresults(IJ+10*Ndata)=TDCIresults(IJ+16*Ndata)
            TDCIresults(IJ+11*Ndata)=TDCIresults(IJ+17*Ndata)
            TDCIresults(IJ+15*Ndata)=tmpx
            TDCIresults(IJ+16*Ndata)=tmpy
            TDCIresults(IJ+17*Ndata)=tmpz
            End Do
          End Do
      ELSE IF (envelope.eq.12) THEN
        WRITE(IOut,*)'                --> testxyL Circular sine squared'
        CALL CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,phase,tt,
     $   evect,fvect2,fvect1,One,-One)
            tmpx=TDCIresults(IJ+9*Ndata)
            tmpy=TDCIresults(IJ+10*Ndata)
            tmpz=TDCIresults(IJ+11*Ndata)
            TDCIresults(IJ+9*Ndata)=TDCIresults(IJ+15*Ndata)
            TDCIresults(IJ+10*Ndata)=TDCIresults(IJ+16*Ndata)
            TDCIresults(IJ+11*Ndata)=TDCIresults(IJ+17*Ndata)
            TDCIresults(IJ+15*Ndata)=tmpx
            TDCIresults(IJ+16*Ndata)=tmpy
            TDCIresults(IJ+17*Ndata)=tmpz
      ELSE 
        WRITE(IOut,*)'                    !! undefined envelope !!'
        STOP
      END IF
C
      OPEN(11,file='efield.dat')
      If(Linear) then
        WRITE(11,1000)fvect1
      else
        WRITE(11,1001)(fvect1(i),fvect2(i),i=1,maxstep)
        endIf
      CLOSE(11)
      WRITE(IOut,*)'Pulse written to: efield.dat (Format:',
     $ ' field strength)'
      Write(IOut,*) "_______________________________________________"
C
C     Read Orbital energies and coefficients, 
C     Vabsorb and dipole moment matrix in the AO basis
C     (ij||ab) and (ia||jb) in the MO basis
C
      TDCIdata = 'TDCI.dat'
      OPEN(10,FILE=TDCIdata)
      Read(10,6000) NBasis,NRorb,NOA,NVA,NOB,NVB,I
      UnRestricted = .False.
      If(I.ne.0) UnRestricted=.True.
      NTT = NBasis*(NBasis+1)/2
      NOANVA = NOA*NVA
      NOBNVB = NOB*NVB
      If(UnRestricted) then
        NStates = NOANVA + NOBNVB + 1
        Write(IOut,1238) NOA,NVA,NOB,NVB
      else
        NStates = NOANVA + 1
        Write(IOut,1239) NOA,NVA
        endif
      NstUse = NStates
C
      Call CPU_TIME(start1)
      len = NRorb
      If(UnRestricted) len = 2*NRorb
      ALLOCATE(OrbEn(len))
      ALLOCATE(CMOA(NRorb*NBasis))
      If(UnRestricted) ALLOCATE(CMOB(NRorb*NBasis))
      ALLOCATE(VabsAO(NBasis*NBasis))
      ALLOCATE(DipXAO(NBasis*NBasis))
      ALLOCATE(DipYAO(NBasis*NBasis))
      ALLOCATE(DipZAO(NBasis*NBasis))
      ALLOCATE(Dijab(NOANVA*NOBNVB))
      Do i=1,NOANVA*NOBNVB
        Dijab(i) = 0.D0
        endDo
      ALLOCATE(DiajbA(NOANVA*NOANVA))
      Do i=1,NOANVA*NOANVA
        DiajbA(i) = 0.D0
        endDo
      If(UnRestricted) then
        ALLOCATE(DiajbB(NOBNVB*NOBNVB))
        Do i=1,NOBNVB*NOBNVB
          DiajbB(i) = 0.D0
          endDo
        endIf
C
      Read(10,6001) (OrbEn(I),I=1,len)
      Write(iout,*) (OrbEn(I),I=1,len)
      Read(10,6001) (CMOA(I),I=1,NRorb*NBasis)
C      write(6,*) "CMOA ", (CMOA(I),I=1,NRorb*NBasis)
      If(UnRestricted) Read(10,6001) (CMOB(I),I=1,NRorb*NBasis)
C      If(UnRestricted) write(6,*) "CMOB ", (CMOB(I),I=1,NRorb*NBasis)
      Write(iout,*) " finished reading orbitals"
      Read(10,6001) (VabsAO(I),I=1,NTT)
      Read(10,6001) (DipXAO(I),I=1,NTT)
      Read(10,6001) (DipYAO(I),I=1,NTT)
      Read(10,6001) (DipZAO(I),I=1,NTT)
      Write(iout,*) " finished one electron integrals"
      If(NOB.ne.0) Read(10,6001) (Dijab(I),I=1,NOANVA*NOBNVB)
      If(NOA.ne.0) Read(10,6001) (DiajbA(I),I=1,NOANVA*NOANVA)
      If(UnRestricted.and.NOB.ne.0) 
     $ Read(10,6001) (DiajbB(I),I=1,NOBNVB*NOBNVB)
      write(iout,*) " finished reading TDCI.dat",
     $ noa,nva,nob,nvb,noanva,nobnvb
      CLOSE(10)
      Call CPU_TIME(finish1)
      Write(IOut,*) " >>> TDCI.data read =             ",finish1-start1
C
C     Form CIS Hamiltonian matrix and diagonalize 
C
      Call CPU_TIME(start1)
      ALLOCATE(CISvec(NStates*NStates))
      ALLOCATE(CISeig(NStates))
      If(UnRestricted) then
        Call FormHU(NOA,NVA,NOANVA,NOB,NVB,NOBNVB,CISvec,OrbEn,
     $    Dijab,DiajbA,DiajbB)
      else
        Call FormH(NOA,NVA,NOANVA,CISvec,OrbEn,Dijab,DiajbA)
        endIf
      write(iout,*) " finished building CI Hamiltonian"
      Lscratch = -1
      CALL DSYEV('Vectors','Upper',NStates,CISvec,NStates,CISeig,
     $ CISeig,Lscratch,INFO) 
      Lscratch = MAX(3*NStates,INT(CISeig(1)))
      Write(IOut,*) " Lscratch= ",Lscratch
      ALLOCATE(scratch(Lscratch))
      CALL DSYEV('Vectors','Upper',NStates,CISvec,NStates,CISeig,
     $ scratch,Lscratch,INFO) 
C undo simple fix that ensure the ground state is the first eigenvector
        CISeig(1) = 0.D0
      Call CPU_TIME(finish1)
      write(iout,*) " finished CI calculation"
      NStUse = NStates
      If(EigMax.gt.0) then
        Do I = 1,NStates
          If(CISeig(i).lt.Ionization+EigMax) NStUse = I
          endDo
      else
        If(EigMax.ne.Zero) NStUse = IAbs(EigMax)
        endIf
      Write(Iout,1240) NStates,NStUse,CISeig(NStUse)
      Write(IOut,*) "_______________________________________________"
      Write(IOut,*) " >>> CIS calculation =            ",finish1-start1
C
C     Transform one electron integrals from AO to MO basis
C     (VabsAO,DipXAO,DipYAO,DipZAO to VabsMO,DipXMO,DipYMO,DipZMO)
C
C     Form absorbing potential and transition dipole matrices
C     for singly excited determinants (ABPmat,TDX,TDY,TDZ)
C
C     Transform matrices to CIS states (ABPmat,TDX,TDY,TDZ) 
C
      Call CPU_TIME(start1)
      Lscratch1 = Max(NBasis*NBasis,NStates)
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
      NThread1 = OMP_get_thread_num()
      Write(IOut,*) " Vabs: Thread ",NThread1
      ALLOCATE(scratch1(Lscratch1))
      ALLOCATE(ABPmat(NStates*NStates))
      ALLOCATE(VabsMOA(NRorb*NRorb))
      Vabs00 = 0.D00
      If(UnRestricted) then
        ALLOCATE(VabsMOB(NRorb*NRorb))
        Call TransAOMO(NBasis,NROrb,VabsAO,VabsMOA,CMOA,scratch1)
        Call TransAOMO(NBasis,NROrb,VabsAO,VabsMOB,CMOB,scratch1)
        Do i=1,NOA
          Vabs00 = Vabs00 + VabsMOA(NRorb*(I-1)+I)
          End Do
        Do i=1,NOB
          Vabs00 = Vabs00 + VabsMOB(NRorb*(I-1)+I)
          End Do
        Call FormMatU(NOA,NVA,NOB,NVB,NRorb,ABPmat,VabsMOA,VabsMOB,
     $    Vabs00)
      else
        Call TransAOMO(NBasis,NROrb,VabsAO,VabsMOA,CMOA,scratch1)
        Do i=1,NOA
          Vabs00 = Vabs00 + 2.D0*VabsMOA(NRorb*(I-1)+I)
          End Do
        Call FormMat(NOA,NVA,NRorb,ABPmat,VabsMOA,Vabs00)
        endIf
c      Call TransCIS(NStates,NStUse,ABPmat,CISvec,scratch1)
      Call Trans1(NStates,NStUse,ABPmat,CISvec,scratch1)
!$OMP SECTION
      NThread2 = OMP_get_thread_num()
      Write(IOut,*) " DipX: Thread ",NThread2
      ALLOCATE(scratch2(Lscratch1))
      ALLOCATE(TDX(NStates*NStates))
      ALLOCATE(DipXMOA(NRorb*NRorb))
      If(UnRestricted) then
        ALLOCATE(DipXMOB(NRorb*NRorb))
        Call TransAOMO(NBasis,NROrb,DipXAO,DipXMOA,CMOA,scratch2)
        Call TransAOMO(NBasis,NROrb,DipXAO,DipXMOB,CMOB,scratch2)
        Call FormMatU(NOA,NVA,NOB,NVB,NRorb,TDX,DipXMOA,DipXMOB,dipx)
      else
        Call TransAOMO(NBasis,NROrb,DipXAO,DipXMOA,CMOA,scratch2)
        Call FormMat(NOA,NVA,NRorb,TDX,DipXMOA,dipx)
        endIf
c      Call TransCIS(NStates,NStUse,TDX,CISvec,scratch2)
      Call Trans1(NStates,NStUse,TDX,CISvec,scratch2)
!$OMP SECTION
      NThread3 = OMP_get_thread_num()
      Write(IOut,*) " DipY: Thread ",NThread3
      ALLOCATE(scratch3(Lscratch1))
      ALLOCATE(TDY(NStates*NStates))
      ALLOCATE(DipYMOA(NRorb*NRorb))
      If(UnRestricted) then
        ALLOCATE(DipYMOB(NRorb*NRorb))
        Call TransAOMO(NBasis,NROrb,DipYAO,DipYMOA,CMOA,scratch3)
        Call TransAOMO(NBasis,NROrb,DipYAO,DipYMOB,CMOB,scratch3)
        Call FormMatU(NOA,NVA,NOB,NVB,NRorb,TDY,DipYMOA,DipYMOB,dipy)
      else
        Call TransAOMO(NBasis,NROrb,DipYAO,DipYMOA,CMOA,scratch3)
        Call FormMat(NOA,NVA,NRorb,TDY,DipYMOA,dipy)
        endIf
c      Call TransCIS(NStates,NStUse,TDY,CISvec,scratch3)
      Call Trans1(NStates,NStUse,TDY,CISvec,scratch3)
!$OMP SECTION
      NThread4 = OMP_get_thread_num()
      Write(IOut,*) " DipZ: Thread ",NThread4
      ALLOCATE(scratch4(Lscratch1))
      ALLOCATE(TDZ(NStates*NStates))
      ALLOCATE(DipZMOA(NRorb*NRorb))
      If(UnRestricted) then
        ALLOCATE(DipZMOB(NRorb*NRorb))
        Call TransAOMO(NBasis,NROrb,DipZAO,DipZMOA,CMOA,scratch4)
        Call TransAOMO(NBasis,NROrb,DipZAO,DipZMOB,CMOB,scratch4)
        Call FormMatU(NOA,NVA,NOB,NVB,NRorb,TDZ,DipZMOA,DipZMOB,dipz)
      else
        Call TransAOMO(NBasis,NROrb,DipZAO,DipZMOA,CMOA,scratch4)
        Call FormMat(NOA,NVA,NRorb,TDZ,DipZMOA,dipz)
        endIf
c      Call TransCIS(NStates,NStUse,TDZ,CISvec,scratch4)
      Call Trans1(NStates,NStUse,TDZ,CISvec,scratch4)
!$OMP END SECTIONS
!$OMP END PARALLEL
      Call CPU_TIME(finish1)
      Write(IOut,*) " >>> AO->MO->CIS transformation = ",
     $ finish1-start1
      write(6,*) " finished transformation"
C
 6003 Format("  State   Energy(au)      gs -> es Transition Dipole")
 6004 Format(I7,F15.10,3F10.4)
 6005 Format("  State  State    es -> es Transition Dipole")
 6006 Format(2I7,3F10.4)
      Write(IOut,6003)
      Do I = 1,11
        Write(IOut,6004) I-1,CISeig(I),TDX(I),TDY(I),TDZ(I)
        end Do
      Write(IOut,6005)
      Do I = 3,11
        Do J = 2,I
          Write(IOut,6006) I-1,J-1,TDX(NStUse*(I-1)+J),
     $    TDY(NStUse*(I-1)+J),TDZ(NStUse*(I-1)+J)
          end Do
        end Do
C
C     If Heuristic.gt.0, add Heuristic*Sqrt(CISig(I)) 
C     to diagonal of ABPmat for unbound states
C
      If(Heuristic.gt.Zero) then
        Write(IOut,1241) Heuristic
        Do I = 1,NStUse
          If(CISeig(I).gt.Ionization) ABPmat(NStUse*(I-1)+I) = 
     $      ABPmat(NStUse*(I-1)+I)+Heuristic*Sqrt(CISeig(I)-Ionization)
          endDo
        endIf
C
C     Form matrix exponential of ABPmat (ExpABP)
C
      Call CPU_TIME(start1)
      ALLOCATE(ExpABP(NStUse*NStUse))
      ALLOCATE(EigABP(NStUse))
C      CALL MATEXP(ABPmat*(-delta/2.0D0),ExpABP,NStUse,Lscratch,scratch1)
      FAC = -delta/2.D0
      CALL MatExp1(ABPmat,EigABP,ExpABP,Fac,NStUse,Lscratch,scratch)
      Call CPU_TIME(finish1)
      Write(IOut,*) " >>> ExpABP time =                ",finish1-start1
C 
C     TDCIS-CAP propagation
C
      len1 = maxstep/outstep
      If(iwrite/10000.gt.0) then
        ALLOCATE(CIdata(Ndata*len1*NStUse))
      else
        ALLOCATE(CIdata(2))
        endIf
      If(Linear) then
        CALL Prop1(CISeig,CISvec,ExpABP,TDX,TDY,TDZ,NStUse,len,len1,
     $   fvect1,maxstep,delta,outstep,NOA,NVA,NOB,NVB,NEmax,Ndir,
     $   TDCIresults,CIdata,Ndata,NStates,iwrite,Lscratch,UnRestricted)
      else
        CALL Prop2(CISeig,CISvec,ExpABP,TDX,TDY,TDZ,NStUse,len,len1,
     $   fvect1,fvect2,maxstep,delta,outstep,NOA,NVA,NOB,NVB,NEmax,
     $   NDir,TDCIresults,CIdata,Ndata,NStates,iwrite,Lscratch,
     $   UnRestricted)
         endIf
C
      Write(Iout,*) " after Prop"
      OPEN(29,File='TDCIresults.dat')
      If(Linear) then
        Write(29,*) (TDCIresults(j),j=1,15*Ndata)
      else
        Write(29,*) (TDCIresults(j),j=1,27*Ndata)
        endIf
      CLOSE(29)
      If(iwrite/10000.gt.0) then
        OPEN(30,File='CIcoeff.dat')
        Write(30,*) (CIdata(j),j=1,Ndata*len1*NStUse)
        CLOSE(30)
        endIf
      CALL DNT(finish)
      WRITE(IOut,5000)finish-start
      WRITE(IOut,*)' '
      CLOSE(IOut)
      END PROGRAM MainIon
C********************************************************************
C
C      !==========================================================!
C      !  Propagation of the TD-CI wfn (linear polarization)      !
C      !==========================================================!
      SUBROUTINE Prop1(H,CISvec,ExpABP,TDX,TDY,TDZ,N,len,len1,fvect,
     $ maxstep,delta,outstep,NOA,NVA,NOB,NVB,NEmax,Ndir,
     $ TDCIresults,CIdata,Ndata,NStates,iwrite,Lscratch,UnRestricted)
      INTEGER N,maxstep,idata,i,j,k,outstep,NEmax,Ndir,NStates,
     $ NOA,NVA,NOB,NVB,NRorb,NTT1,len,len1,Lscratch,INFO,Ndata,iwrite
      COMPLEX*16 ExpH(N), A
      COMPLEX*16 Psi(N),Psi1(NStates),Psi2(NStates),CIdata(*) 
      REAL*8 Pop(len),PopI(len),PopInt(len),NormData(maxstep/100),
     $  PopData(len,len1),IonData(len+3,len1),MuData(3,len1),
     $  HP(N*N),TDvals(N),scratch(Lscratch)
      REAL*8 H(N),CISvec(N*N),ExpABP(N*N),TDX(N*N),TDY(N*N),TDZ(N*N),
     $  fvect(maxstep),emax,efield,delta,Norm,temp,Temp1,ExpectA,
     $  Zero,Half,One,TDCIresults(Ndata,*)
      REAL start,finish,start1,finish1
      Character Char
      Character*4 DirStr,EmaxStr
      Character*60 NormFile,PopFile,IonFile,DipoleFile
      INTEGER OMP_get_num_threads,OMP_get_num_procs,OMP_get_thread_num,
     $ NThread
      LOGICAL UnRestricted
C      COMMON IOut
C
 9991 Format(" Thread",I3,"  propagation for direction",I4,
     $ " and intensity",I4)
 9992 Format(" Dir= (",F8.5,",",F8.5,",",F8.5,") Emax=",F8.5)
 9993 FORMAT(10(' ',F7.5))
 9994 FORMAT("    t                   Psi(1)",
     $ "                     Norm",
     $ "                       mu_z")
 9995 FORMAT(24(' ',F7.5))
 9996 FORMAT(F9.2, ' ',4D20.10)
 9997 FORMAT(F9.2, ' ',24(' ',F7.5))
 9998 FORMAT(F9.2, ' ',F26.10, ' ',F26.10, ' ',F26.10)
 9999 FORMAT(F9.2, ' ',F26.10) 
C
      IOut=42
      write(IOut,*) "Start Prop1"
      Zero = 0.0D0
      Half = 0.5D0
      One  = 1.0D0
      NRorb = NOA+NVA
      NTT1 = (NOA+3)*(NOA+4)/2
      Lscratch = N*N
      len = NRorb
      If(UnRestricted) len = 2*NRorb
      len1 = maxstep/outstep
C
C     ExpH = Exp(-iH*delta/2)
C
      Do i=1,N
        temp = -Half*H(i)*delta
        ExpH(i)=DCMPLX(DCOS(temp),DSIN(temp))
        End Do
C
C     Start loop over directions
C
!$OMP PARALLEL DEFAULT(NONE),
!$OMP& PRIVATE(i,j,jj,k,kk,NThread,Psi,Psi1,Psi2,HP,efield,temp,temp1,
!$OMP& A,NORM,normfile,popfile,ionfile,dipolefile,TDvals,
!$OMP& scratch,IDir,DirStr,dirx,diry,dirz,IEmax,EmaxStr,emax,
!$OMP& Pop,PopI,PopInt,idata,NormData,PopData,IonData,MuData,
!$OMP& INFO,start,finish),
!$OMP& SHARED(ExpH,CISvec,ExpABP,TDX,TDY,TDZ,fvect,
!$OMP& NDir,NEmax,delta,Half,Zero,One,
!$OMP& NOA,NVA,NOB,NVB,NRorb,N,maxstep,len,len1,Lscratch,outstep,
!$OMP& iwrite,IOut,UnRestricted,TDCIresults)
!$OMP DO
      Do IDir=1,NDir
      DirStr=Char(48+IDir/100)//Char(48+Mod(IDir/10,10))//
     $ Char(48+Mod(IDir,10))
      dirx = TDCIresults(IDir,4)
      diry = TDCIresults(IDir,5)
      dirz = TDCIresults(IDir,6)
C
C     Form the transition dipole in the (dirx,diry,dirz) direction,
C     diagonalize and multiply eigenvectors by ExpABP
C 
      Call CPU_TIME(start)
      Do I = 1,N*N
         HP(I) = dirx*TDX(I)+diry*TDY(I)+dirz*TDZ(I)
         end Do
      CALL DSYEV('Vectors','Upper',N,HP,N,TDvals,scratch,Lscratch,
     $ INFO) 
      CALL MatMult1(HP,ExpABP,scratch,N)
      Call CPU_TIME(finish)
C
C     Start loop over intensities
C
      Do IEmax=1,NEmax
      EmaxStr=Char(48+IEmax/100)//Char(48+Mod(IEmax/10,10))//
     $ Char(48+Mod(IEmax,10))
      emax = TDCIresults(1+(IEmax-1)*NDir,1)
C
!$OMP CRITICAL 
      NThread = omp_get_thread_num()
      If(IEmax.eq.1) then
        Write(IOut,*) " >>> TD diag and TDvec*ExpABP  =  ",finish-start
        endif
      WRITE(IOut,9991) NThread,IDir,IEmax
      Write(Iout,9992) dirx,diry,dirz,emax
!$OMP END CRITICAL
C
C     Initialize Ion population and CI wave function
      Do i=1,len
        PopInt(i)=0.D0
        End Do
      Do i=1,N
        Psi(i)=DCMPLX(0.D0,0.D0)
        End Do
      Psi(1)=DCMPLX(1.D0,0.D0)
C
C     Begin Looping over the Time Steps
      Call CPU_TIME(start)
C      i=1
C      If(len.lt.0) then
      DO i=1, maxstep
C
        efield = emax*(fvect(i)+fvect(i+1))*Half
C       write(6,*) i,Psi(i),norm,efield
        temp1 = delta*efield
        Do j = 1, N
          Psi(j) = ExpH(j)*Psi(j)
          End Do
        Do j = 1, N
          Psi1(j) = DCMPLX(0.D0,0.D0)
          End Do
c      write(6,*) " N,Psi,Psi1,HP,Tvals,temp",N,Psi(1),Psi1(1),HP(1),
c     $ TDvals(1),temp
c      If(len.lt.0) then
        Do j = 1, N
          jj = (j-1)*N
          Do k = 1, N
            Psi1(k) = Psi1(k) + HP(jj+k)*Psi(j)
            End Do
          End Do
C      If(IEmax.eq.0) then
        Do j = 1, N
          temp = temp1 * TDvals(j)
          Psi1(j) = DCMPLX(DCOS(temp),DSIN(temp))*Psi1(j)
          End Do
        Do k = 1, N
          kk = N*(k-1)
          A = DCMPLX(Zero,Zero)
          Do j = 1, N
            A = A + HP(j+kk)*Psi1(j)
            End Do
            Psi(k) = ExpH(k)*A
          End Do
C
        CALL NORMsr(PSI, NORM, N, 'I')
        If(Mod(i,100).eq.0) NormData(i/100) = NORM
C       If(Mod(i,100).eq.0) Write(IOut,*) i,Psi(1),NORM,temp1
        If(Mod(i,outstep).eq.0) then
          idata=i/outstep
          Call Analysis(NOA,NVA,NOB,NVB,N,len,
     $    Psi,Psi1,Psi2,CISvec,ExpABP,TDX,TDY,TDZ,Pop,PopI,PopInt,
     $    idata,NORM,delta,efield,PopData,IonData,MuData,
     $    UnRestricted)
C         Write(6,9998)i*delta,ABS(Psi(1)),NORM,MuData(3,idata) 
          If(iwrite/10000.gt.0) then
            j = IDir + (IEmax-1)*NDir
            Do k = 1, N
              CIdata(N*len1*(j-1)+N*(idata-1)+k) = Psi(k)
              endDo
            endIf
          endif
c      endif
        END DO
C      endif
C
C     End of loop over time steps
C
      temp = Zero
      Do i=1,len
        temp = temp + PopInt(i)
        END DO
      temp = (NOA+NOB)/temp
      Do i=1,len
        PopInt(i)= PopInt(i)*temp
        END DO
C
C     Write Norm, Pop, Ion and Dipole data to files
C
!$OMP CRITICAL 
      i = IDir + (IEmax-1)*NDir
      TDCIresults(i,7) = NORM
      TDCIresults(i,8) = PopInt(NOA)
      TDCIresults(i,9) = PopInt(NOA-1)
      TDCIresults(i,10) = PopInt(NOA-2)
      TDCIresults(i,11) = PopInt(NOA-3)
      If(Unrestricted) then
        TDCIresults(i,12) = PopInt(NRorb+NOB)
        TDCIresults(i,13) = PopInt(NRorb+NOB-1)
        TDCIresults(i,14) = PopInt(NRorb+NOB-2)
        TDCIresults(i,15) = PopInt(NRorb+NOB-3)
      else
        TDCIresults(i,12) = PopInt(NOA-4)
        TDCIresults(i,13) = PopInt(NOA-5)
        TDCIresults(i,14) = PopInt(NOA-6)
        TDCIresults(i,15) = PopInt(NOA-7)
      endIf
C     Write(IOut,*) (TDCIresults(i,j),j=1,15)
      If(Mod(iwrite,10).ne.0) then
       NormFile='norm-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(18,File=Trim(NormFile))
        WRITE(18,9999) Zero, One
        Do idata = 1,maxstep/100
          WRITE(18,9999)100*idata*delta,NormData(idata)
          END DO
        WRITE(IOut,*)'Norm written to: ',Trim(NormFile),
     $    ' (Format: t, norm)'
       CLOSE(18)
       endIf
      If(Mod(iwrite/10,10).ne.0) then
       PopFile='pop-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(19,File=Trim(PopFile))
        Do idata = 1,len1
          WRITE(19,9997)outstep*idata*delta,(PopData(j,idata),j=1,len)
          END DO
        WRITE(IOut,*)'MO populations written to:',Trim(PopFile),
     $    ' (Format: t,pop(1),pop(2),pop(3),etc)'
       CLOSE(19)
       endIf
      If(Mod(iwrite/100,10).ne.0) then
       IonFile='ion-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(20,File=Trim(IonFile))
        Do idata = 1,len1
          WRITE(20,9996)outstep*idata*delta,IonData(1,idata),
     $      IonData(2,idata),IonData(3,idata)
          WRITE(20,9995)(IonData(j+3,idata),j=1,len)
          END DO
        WRITE(IOut,*)'Instantaneous rate & ion populations written to:',
     $   Trim(IonFile),' (Format: t, efield, rate, norm',
     $   ' pop(1), pop(2), etc)'
        Write(IOut,*) "  Integrated MO populations for the ion"
c        Write(IOut,9993) (PopInt(i),i=1,len)
       CLOSE(20)
       endIf
      If(Mod(iwrite/1000,10).ne.0) then
       DipoleFile='mu-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(28,File=Trim(DipoleFile))
        WRITE(28,9998) Zero,TDX(1),TDY(1),TDZ(1)
        Do idata = 1,len1
          WRITE(28,9998)outstep*idata*delta,MuData(1,idata),
     $      MuData(2,idata),MuData(3,idata)
          END DO
        WRITE(IOut,*)'Dipole moments written to: ',Trim(DipoleFile),
     $    ' (Format: t, mu_x, mu_y, mu_z)'
       CLOSE(28)
       endIf
      Write(Iout,*) "  Final Norm = ",NORM
      Call CPU_TIME(finish)
      Write(IOut,*) " >>> Propagation time =           ",finish-start
!$OMP END CRITICAL 
C
C     End of loop over intensities
C
      End Do
C
C     End of loop over directions
C
      End Do
!$OMP END DO
!$OMP END PARALLEL 
C      OPEN(29,File='TDCIresults.dat')
C      Do i = 1,NData
C        Write(29,*) (TDCIresults(i,j),j=1,15)
C        End Do
C      CLOSE(29)
      WRITE(IOut,*)'Done propagation1'
      WRITE(IOut,*)' '
      END SUBROUTINE Prop1

C********************************************************************
C
C      !==========================================================!
C      !  Propagation of the TD-CI wfn (circular polarization)    !
C      !==========================================================!
      SUBROUTINE Prop2(H,CISvec,ExpABP,TDX,TDY,TDZ,N,len,len1,
     $ fvect1,fvect2,maxstep,delta,outstep,NOA,NVA,NOB,NVB,NEmax,
     $ Ndir,TDCIresults,CIdata,Ndata,NStates,iwrite,Lscratch,
     $ UnRestricted)
      INTEGER N,maxstep,idata,i,j,k,outstep,NEmax,Ndir,NStates,
     $ NOA,NVA,NOB,NVB,NRorb,NTT1,len,len1,Lscratch,INFO,Ndata,iwrite
      COMPLEX*16 ExpH(N), A
      COMPLEX*16 Psi(N),Psi1(NStates),Psi2(NStates),CIdata(*) 
      REAL*8 Pop(len),PopI(len),PopInt(len),NormData(maxstep/100),
     $ PopData(len,len1),IonData(len+3,len1),MuData(3,len1),
     $ HP1(N*N),HP2(N*N),TDvals1(N),TDvals2(N),scratch(Lscratch)
      REAL*8 H(N),CISvec(N*N),ExpABP(N*N),TDX(N*N),TDY(N*N),TDZ(N*N),
     $ fvect1(maxstep),fvect2(maxstep),delta,Norm,temp,Temp1,Temp2,
     $ ExpectA,Zero,Half,One,TDCIresults(Ndata,*),
     $ dirx1,diry1,dirz1,emax1,efield1,dirx2,diry2,dirz2,emax2,efield2
      REAL start,finish,start1,finish1
      Character Char
      Character*4 DirStr,EmaxStr
      Character*60 NormFile,PopFile,IonFile,DipoleFile
      INTEGER OMP_get_num_threads,OMP_get_num_procs,OMP_get_thread_num,
     $ NThread
      LOGICAL UnRestricted
C      COMMON IOut
C
 9991 Format(" Thread",I3,"  propagation for direction",I4,
     $ " and intensity",I4)
 9992 Format(" Dir= (",F8.5,",",F8.5,",",F8.5,") Emax=",F8.5)
 9993 FORMAT(10(' ',F7.5))
 9994 FORMAT("    t                   Psi(1)",
     $ "                     Norm",
     $ "                       mu_z")
 9995 FORMAT(24(' ',F7.5))
 9996 FORMAT(F9.2, ' ',4D20.10)
 9997 FORMAT(F9.2, ' ',24(' ',F7.5))
 9998 FORMAT(F9.2, ' ',F26.10, ' ',F26.10, ' ',F26.10)
 9999 FORMAT(F9.2, ' ',F26.10) 
C
      IOut=42
      write(IOut,*) "Start Prop2"
      Zero = 0.0D0
      Half = 0.5D0
      One  = 1.0D0
      NRorb = NOA+NVA
      NTT1 = (NOA+3)*(NOA+4)/2
      Lscratch = N*N
      len = NRorb
      If(UnRestricted) len = 2*NRorb
      len1 = maxstep/outstep
C
C     ExpH = Exp(-iH*delta/2)
C
      Do i=1,N
        temp = -Half*H(i)*delta
        ExpH(i)=DCMPLX(DCOS(temp),DSIN(temp))
        End Do
C
C     Start loop over directions
C
!$OMP PARALLEL DEFAULT(NONE),
!$OMP& PRIVATE(i,j,jj,k,kk,NThread,Psi,Psi1,Psi2,HP1,HP2,
!$OMP& efield1,efield2,temp,temp1,temp2,A,NORM,normfile,
!$OMP& popfile,ionfile,dipolefile,TDvals1,TDvals2,scratch,
!$OMP& IDir,DirStr,dirx1,diry1,dirz1,dirx2,diry2,dirz2,
!$OMP& IEmax,EmaxStr,emax1,emax2,
!$OMP& Pop,PopI,PopInt,idata,NormData,PopData,IonData,MuData,
!$OMP& INFO,start,finish),
!$OMP& SHARED(ExpH,CISvec,ExpABP,TDX,TDY,TDZ,fvect1,fvect2,
!$OMP& NDir,NEmax,delta,Half,Zero,One,
!$OMP& NOA,NVA,NOB,NVB,NRorb,N,maxstep,len,len1,Lscratch,outstep,
!$OMP& iwrite,IOut,UnRestricted,TDCIresults,CIdata)
!$OMP DO
C
C     Circular polarization defined by two directions
C
      Do IDir=1,NDir
      DirStr=Char(48+IDir/100)//Char(48+Mod(IDir/10,10))//
     $ Char(48+Mod(IDir,10))
      dirx1 = TDCIresults(IDir,10)
      diry1 = TDCIresults(IDir,11)
      dirz1 = TDCIresults(IDir,12)
      dirx2 = TDCIresults(IDir,16)
      diry2 = TDCIresults(IDir,17)
      dirz2 = TDCIresults(IDir,18)
C
C     Form the transition dipole in the 
C     (dirx1,diry1,dirz1) and (dirx2,diry2,dirz2) directions,
C     diagonalize and form products
C 
      Call CPU_TIME(start)
      Do I = 1,N*N
         HP1(I) = dirx1*TDX(I)+diry1*TDY(I)+dirz1*TDZ(I)
         HP2(I) = dirx2*TDX(I)+diry2*TDY(I)+dirz2*TDZ(I)
         end Do
      CALL DSYEV('Vectors','Upper',N,HP1,N,TDvals1,scratch,Lscratch,
     $ INFO) 
      CALL DSYEV('Vectors','Upper',N,HP2,N,TDvals2,scratch,Lscratch,
     $ INFO) 
      CALL MatMult1(HP2,HP1,scratch,N)
      CALL MatMult1(HP1,ExpABP,scratch,N)
c      Do j=1,10
c        Write(IOut,9993) (HP2((j-1)*N+i),i=1,10)
c        End Do
      Call CPU_TIME(finish)
C
C     Start loop over intensities
C
      Do IEmax=1,NEmax
      EmaxStr=Char(48+IEmax/100)//Char(48+Mod(IEmax/10,10))//
     $ Char(48+Mod(IEmax,10))
      emax1 = TDCIresults(1+(IEmax-1)*NDir,7)
      emax2 = TDCIresults(1+(IEmax-1)*NDir,13)
C
!$OMP CRITICAL 
      NThread = omp_get_thread_num()
      If(IEmax.eq.1) then
        Write(IOut,*) " >>> TD diag and TDvec*ExpABP  =  ",finish-start
        endif
      WRITE(IOut,9991) NThread,IDir,IEmax
      Write(Iout,9992) dirx1,diry1,dirz1,emax1
      Write(Iout,9992) dirx2,diry2,dirz2,emax2
!$OMP END CRITICAL
C
C     Initialize Ion population and CI wave function
      Do i=1,len
        PopInt(i)=0.D0
        End Do
      Do i=1,N
        Psi(i)=DCMPLX(0.D0,0.D0)
        End Do
      Psi(1)=DCMPLX(1.D0,0.D0)
C
C     Begin Looping over the Time Steps
      Call CPU_TIME(start)
C      i=1
C      If(len.lt.0) then
      DO i=1, maxstep
C
        efield1 = emax1*(fvect1(i)+fvect1(i+1))*Half
        efield2 = emax2*(fvect2(i)+fvect2(i+1))*Half
        temp1 = delta*efield1*half
        temp2 = delta*efield2
        Do j = 1, N
          Psi(j) = ExpH(j)*Psi(j)
          End Do
        Do j = 1, N
          Psi1(j) = DCMPLX(0.D0,0.D0)
          End Do
c      write(6,*) " N,Psi,Psi1,HP,Tvals,temp",N,Psi(1),Psi1(1),HP1(1),
c     $ TDvals1(1),temp1,HP2(1),TDvals2(1),temp2
        Do j = 1, N
          jj = (j-1)*N
          Do k = 1, N
            Psi1(k) = Psi1(k) + HP1(jj+k)*Psi(j)
            End Do
          End Do
        Do j = 1, N
          temp = temp1 * TDvals1(j)
          Psi1(j) = DCMPLX(DCOS(temp),DSIN(temp))*Psi1(j)
          End Do
        Do j = 1, N
          Psi(j) = DCMPLX(0.D0,0.D0)
          End Do
        Do j = 1, N
          jj = (j-1)*N
          Do k = 1, N
            Psi(k) = Psi(k) + HP2(jj+k)*Psi1(j)
            End Do
          End Do
        Do j = 1, N
          temp = temp2 * TDvals2(j)
          Psi(j) = DCMPLX(DCOS(temp),DSIN(temp))*Psi(j)
          End Do
        Do k = 1, N
          kk = N*(k-1)
          A = DCMPLX(Zero,Zero)
          Do j = 1, N
            A = A + HP2(j+kk)*Psi(j)
            End Do
            Psi1(k) = A
          End Do
          Do j = 1, N
          temp = temp1 * TDvals1(j)
          Psi1(j) = DCMPLX(DCOS(temp),DSIN(temp))*Psi1(j)
          End Do
      Do k = 1, N
          kk = N*(k-1)
          A = DCMPLX(Zero,Zero)
          Do j = 1, N
            A = A + HP1(j+kk)*Psi1(j)
            End Do
            Psi(k) = ExpH(k)*A
          End Do
C
        CALL NORMsr(PSI, NORM, N, 'I')
        If(Mod(i,100).eq.0) NormData(i/100) = NORM
c        If(Mod(i,100).eq.0) Write(IOut,*) i,Psi(1),NORM,temp1,temp2
        If(Mod(i,outstep).eq.0) then
          idata=i/outstep
          Call Analysis(NOA,NVA,NOB,NVB,N,len,
     $    Psi,Psi1,Psi2,CISvec,ExpABP,TDX,TDY,TDZ,Pop,PopI,PopInt,
     $    idata,NORM,delta,efield1,PopData,IonData,MuData,
     $    UnRestricted)
C          Write(IOut,9998)i*delta,ABS(Psi(1)),NORM,MuData(3,idata) 
          If(iwrite/10000.gt.0) then
            j = IDir + (IEmax-1)*NDir
            Do k = 1, N
              CIdata(N*len1*(j-1)+N*(idata-1)+k) = Psi(k)
              endDo
            endIf
          endif
c      endif
        END DO
C      endif
C
C     End of loop over time steps
C
      temp = Zero
      Do i=1,len
        temp = temp + PopInt(i)
        END DO
      temp = (NOA+NOB)/temp
      Do i=1,len
        PopInt(i)= PopInt(i)*temp
        END DO
C
C     Write Norm, Pop, Ion and Dipole data to files
C
!$OMP CRITICAL 
      i = IDir + (IEmax-1)*NDir
      TDCIresults(i,19) = NORM
      TDCIresults(i,20) = PopInt(NOA)
      TDCIresults(i,21) = PopInt(NOA-1)
      TDCIresults(i,22) = PopInt(NOA-2)
      TDCIresults(i,23) = PopInt(NOA-3)
      If(Unrestricted) then
        TDCIresults(i,24) = PopInt(NRorb+NOB)
        TDCIresults(i,25) = PopInt(NRorb+NOB-1)
        TDCIresults(i,26) = PopInt(NRorb+NOB-2)
        TDCIresults(i,27) = PopInt(NRorb+NOB-3)
      else
        TDCIresults(i,24) = PopInt(NOA-4)
        TDCIresults(i,25) = PopInt(NOA-5)
        TDCIresults(i,26) = PopInt(NOA-6)
        TDCIresults(i,27) = PopInt(NOA-7)
      endIf
C     Write(IOut,*) (TDCIresults(i,j),j=1,15)
      If(Mod(iwrite,10).ne.0) then
       NormFile='norm-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(18,File=Trim(NormFile))
        WRITE(18,9999) Zero, One
        Do idata = 1,maxstep/100
          WRITE(18,9999)100*idata*delta,NormData(idata)
          END DO
        WRITE(IOut,*)'Norm written to: ',Trim(NormFile),
     $    ' (Format: t, norm)'
       CLOSE(18)
       endIf
      If(Mod(iwrite/10,10).ne.0) then
       PopFile='pop-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(19,File=Trim(PopFile))
        Do idata = 1,len1
          WRITE(19,9997)outstep*idata*delta,(PopData(j,idata),j=1,len)
          END DO
        WRITE(IOut,*)'MO populations written to:',Trim(PopFile),
     $    ' (Format: t,pop(1),pop(2),pop(3),etc)'
       CLOSE(19)
       endIf
      If(Mod(iwrite/100,10).ne.0) then
       IonFile='ion-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(20,File=Trim(IonFile))
        Do idata = 1,len1
          WRITE(20,9996)outstep*idata*delta,IonData(1,idata),
     $      IonData(2,idata),IonData(3,idata)
          WRITE(20,9995)(IonData(j+3,idata),j=1,len)
          END DO
        WRITE(IOut,*)'Instantaneous rate & ion populations written to:',
     $   Trim(IonFile),' (Format: t, efield1, rate, norm',
     $   ' pop(1), pop(2), etc)'
        Write(IOut,*) "  Integrated MO populations for the ion"
        Write(IOut,9993) (PopInt(i),i=1,len)
       CLOSE(20)
       endIf
      If(Mod(iwrite/1000,10).ne.0) then
       DipoleFile='mu-e'//Trim(EmaxStr)//'-d'//Trim(DirStr)//'.dat'
       OPEN(28,File=Trim(DipoleFile))
        WRITE(28,9998) Zero,TDX(1),TDY(1),TDZ(1)
        Do idata = 1,len1
          WRITE(28,9998)outstep*idata*delta,MuData(1,idata),
     $      MuData(2,idata),MuData(3,idata)
          END DO
        WRITE(IOut,*)'Dipole moments written to: ',Trim(DipoleFile),
     $    ' (Format: t, mu_x, mu_y, mu_z)'
       CLOSE(28)
       endIf
      Write(Iout,*) "  Final Norm = ",NORM
      Call CPU_TIME(finish)
      Write(IOut,*) " >>> Propagation time =           ",finish-start
!$OMP END CRITICAL 
C
C     End of loop over intensities
C
      End Do
C
C     End of loop over directions
C
      End Do
!$OMP END DO
!$OMP END PARALLEL 
C      OPEN(29,File='TDCIresults.dat')
C      Do i = 1,NData
C        Write(29,*) (TDCIresults(i,j),j=1,27)
C        End Do
C      CLOSE(29)
      WRITE(IOut,*)'Done propagation2'
      WRITE(IOut,*)' '
      END SUBROUTINE Prop2
 
C      !==========================================================!
C      ! Analysis of the wavefunction                             !
C      !==========================================================!
      Subroutine Analysis(NOA,NVA,NOB,NVB,N,len,Psi,Psi1,Psi2,
     $ CISvec,ExpABP,TDX,TDy,TDZ,Pop,PopI,PopInt,idata,norm,
     $ delta,efield,PopData,IonData,MuData,UnRestricted)
      COMPLEX*16 Psi(*),Psi1(*),Psi2(*)
      REAL*8 CISvec(*),ExpABP(*),TDX(*),TDY(*),TDZ(*),
     $ Pop(len),PopI(len),PopInt(*),
     $ PopData(len,*),IonData(len+3,*),MuData(3,*),
     $ delta,Rate,ExpectA,efield,norm,norm1
      INTEGER idata,j,k,NOA,NVA,NOB,NVB,N,len
      LOGICAL UnRestricted
C
      If(UnRestricted) then
        Call PopUCIS(NOA,NVA,NOB,NVB,N,Psi,Psi1,CISvec,Pop)
        Call PopIonU(NOA,NVA,NOB,NVB,N,Psi,Psi1,Psi2,
     $    ExpABP,CISvec,Rate,norm1,PopI,PopInt)
      else
        Call PopCIS(NOA,NVA,N,Psi,Psi1,CISvec,Pop)
        Call PopIon(NOA,NVA,N,Psi,Psi1,Psi2,ExpABP,CISvec,Rate,
     $    norm1,PopI,PopInt)
        endIf
      IonData(1,idata) = efield
      IonData(2,idata) = Rate/delta
      IonData(3,idata) = norm1
      Do j=1,len
        PopData(j,idata) = Pop(j)
        IonData(j+3,idata) = PopI(j)
        End Do
      MuData(1,idata) = -ExpectA(Psi,TDX,N)/(norm**2)
      MuData(2,idata) = -ExpectA(Psi,TDY,N)/(norm**2)
      MuData(3,idata) = -ExpectA(Psi,TDZ,N)/(norm**2)
      END SUBROUTINE Analysis

C      !==========================================================!
C      ! Calculate Expectation value <A> = Conjugate(PSI)*A*PSI   !
C      !==========================================================!
      FUNCTION ExpectA(PSI, A, N)
      INTEGER i, ii, j, ij, N
      REAL*8 :: ExpectA
      REAL*8 A(N*N), B
      COMPLEX*16 PSI(N), C
C
      B = 0.D0
      DO i = 2,N
        ii = (i-1)*N
        C = DCONJG(PSI(i))
        DO j = 1,(i-1)
          ij = ii + j
          B = B + A(ij)*REAL(C*PSI(j))
        END DO
      END DO
      B = 2*B
      DO i = 1,N
        ii = (i-1)*N + i
        B = B + A(ii)*REAL(PSI(i)*DCONJG(PSI(i)))
      END DO
      ExpectA = B
      END FUNCTION ExpectA

C      !==========================================================!
C      ! Form CIS Hamiltonian                                     !
C      ! Hmat(I,A,J,B) = delta(i,j)delta(a,b)(en(a)-en(i))        !
C      !                 +(ij||ab)-(ib||ja)                       !
C      ! storage order Dijab -> iajb, Diajb -> iajb               !
C      !==========================================================!
      SUBROUTINE FormH(NOA,NVA,NOANVA,Hmat,En,Dijab,Diajb)
      INTEGER I,J,K,KK,L,A,B,IA,IB,JA,JB,NOA,NVA,NOANVA,NStates
      REAL*8 Hmat(*),En(*),
     $ Dijab(NOANVA,*),Diajb(NOANVA,*)
C
      NStates = NOANVA + 1
      DO K = 1, NStates
        Hmat((K-1)*NStates+1) = 0.D0
        Hmat(K) = 0.D0
        END DO
C simple fix to make sure that the ground state is the first eigenvector
C (need to undo this fix after diagonalization)
        Hmat(1) = -100.D0
C
      DO K = 2, NStates
        KK = (K-1)*NStates
        IA = K-1
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)
        DO L = 2, NStates
          JB = L-1
          J = (JB-1)/NVA+1
          B = JB-NVA*(J-1)
          IB = NVA*(I-1)+B
          JA = NVA*(J-1)+A
          Hmat(KK+L) = Dijab(IA,JB)-Diajb(IB,JA)
          END DO
        Hmat(KK+K) = Hmat(KK+K)+En(A+NOA)-En(I)
        END DO
      END SUBROUTINE FormH

C      !==========================================================!
C      ! Form unrestricted CIS Hamiltonian                        !
C      !                                                          !
C      ! Hmat(I,A,J,B)=delta(i,j)delta(a,b)(en(a)-en(i))-(ib||ja) !
C      !                                                          !
C      !              Gaussian symmetry:  sequence:               !
C      !                                                          !
C      ! Bucket 2  Dijab   abab  all        iajb                  !
C      ! Bucket 5  DiajbA  aaaa  all        iajb                  !
C      ! Bucket 10 DiajbB  bbbb  all        iajb                  !
C      !==========================================================!
      SUBROUTINE FormHU(NOA,NVA,NOANVA,NOB,NVB,NOBNVB,Hmat,En,
     $ Dijab,DiajbA,DiajbB)
      INTEGER I,J,K,KK,L,LL,A,B,IA,IB,JA,JB,
     $ NOA,NVA,NOANVA,NOB,NVB,NOBNVB,NRorb,NStates
      REAL*8 Hmat(*),En(*),
     $ Dijab(NOBNVB,*),DiajbA(NOANVA,*),DiajbB(NOBNVB,*)
C
      NRorb = NOA + NVA
      NStates = NOANVA + NOBNVB + 1
C
      DO K = 1, NStates
        Hmat((K-1)*NStates+1) = 0.D0
        Hmat(K) = 0.D0
        END DO
C simple fix to make sure that the ground state is the first eigenvector
C (need to undo this fix after diagonalization)
        Hmat(1) = -100.D0
C
C     alpha,alpha block
C
      DO IA = 1, NOANVA
        K = IA+1
        KK = (K-1)*NStates
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)
        DO JB = 1, NOANVA
          L = JB+1
          J = (JB-1)/NVA+1
          B = JB-NVA*(J-1)
          IB = NVA*(I-1)+B
          JA = NVA*(J-1)+A
          Hmat(KK+L) = -DiajbA(IB,JA)
          END DO
        Hmat(KK+K) = Hmat(KK+K)+En(A+NOA)-En(I)
        END DO
      if(NOBNVB.ne.0) then
C
C     beta,beta block
C
      DO IA = 1, NOBNVB
        K = IA+NOANVA+1
        KK = (K-1)*NStates
        I = (IA-1)/NVB+1
        A = IA-NVB*(I-1)
        DO JB = 1, NOBNVB
          L = JB+NOANVA+1
          J = (JB-1)/NVB+1
          B = JB-NVB*(J-1)
          IB = NVB*(I-1)+B
          JA = NVB*(J-1)+A
          Hmat(KK+L) = -DiajbB(IB,JA)
          END DO
        Hmat(KK+K) = Hmat(KK+K)+En(A+NOB+NRorb)-En(I+NRorb)
        END DO
C
C     alpha,beta and beta,alpha blocks
C
      DO IA = 1, NOANVA
        K = IA+1
        KK = (K-1)*NStates
        DO JB = 1, NOBNVB
          L = JB+NOANVA+1
          LL = (L-1)*NStates
          Hmat(KK+L) = Dijab(JB,IA)
          Hmat(LL+K) = Hmat(KK+L)
          END DO
        END DO
        endIf
      END SUBROUTINE FormHU

C      !==========================================================!
C      ! Form one electron matrix elements for CIS states         !
C      ! Xmat(I,A,0,0) = X(I,A)Sqrt(2)                            !
C      ! Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)        !
C      !==========================================================!
      SUBROUTINE FormMat(NOA,NVA,NRorb,Xmat,X,X00)
      INTEGER I,J,K,KK,L,A,B,IA,JB,NOA,NVA,NStates
      REAL*8 Xmat(*),X(NRorb,*),X00
C
      NStates = NOA*NVA + 1
      XMat(1) = X00
      DO K = 2, NStates
        KK = NStates*(K-1)
        IA = K-1
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        Xmat(KK+1) = X(I,A)*Sqrt(2.D0)
        Xmat(K) = Xmat(KK+1)
        END DO
C
      DO K = 2, NStates
        KK = NStates*(K-1)
        IA = K-1
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        DO L = 2, NStates
          JB = L-1
          J = (JB-1)/NVA+1
          B = JB-NVA*(J-1)+NOA
          Xmat(KK+L) = 0.D0
          If(I.eq.J) Xmat(KK+L) = Xmat(KK+L)+X(A,B)
          If(A.eq.B) Xmat(KK+L) = Xmat(KK+L)-X(I,J)
          END DO
        Xmat(KK+K) = Xmat(KK+K) + X00
        END DO
      END SUBROUTINE FormMat

C      !==========================================================!
C      ! Form one electron matrix elements for UCIS states        !
C      ! Xmat(I,A,0,0) = X(I,A)                                   !
C      ! Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)        !
C      !==========================================================!
      SUBROUTINE FormMatU(NOA,NVA,NOB,NVB,NRorb,Xmat,XA,XB,X00)
      INTEGER I,J,K,KK,L,A,B,IA,JB,
     $ NOA,NVA,NOB,NVB,NOANVA,NOBNVB,NRorb,NStates
      REAL*8 Xmat(*),XA(NRorb,*),XB(NRorb,*),X00
C
      NOANVA = NOA*NVA
      NOBNVB = NOB*NVB
      NStates = NOANVA + NOBNVB + 1
C
      DO K = 1, NStates*NStates
        XMat(K) = 0.D0
        END DO
      XMat(1) = X00
C
C     Ground state and alpha excitations
C
      DO IA = 1, NOANVA
        K = IA+1
        KK = NStates*(K-1)
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        Xmat(KK+1) = XA(I,A)
        Xmat(K) = Xmat(KK+1)
        END DO
C
C     Ground state and beta excitations
C
      DO IA = 1, NOBNVB
        K = IA+NOANVA+1
        KK = NStates*(K-1)
        I = (IA-1)/NVB+1
        A = IA-NVB*(I-1)+NOB
        Xmat(KK+1) = XB(I,A)
        Xmat(K) = Xmat(KK+1)
        END DO
C
C     alpha to alpha excitations
C
      DO IA = 1, NOANVA
        K = IA+1
        KK = NStates*(K-1)
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        DO JB = 1, NOANVA
          L = JB+1
          J = (JB-1)/NVA+1
          B = JB-NVA*(J-1)+NOA
          If(I.eq.J) Xmat(KK+L) = Xmat(KK+L)+XA(A,B)
          If(A.eq.B) Xmat(KK+L) = Xmat(KK+L)-XA(I,J)
          END DO
        Xmat(KK+K) = Xmat(KK+K) + X00
        END DO
C
C     beta to beta excitations
C
      DO IA = 1, NOBNVB
        K = IA+NOANVA+1
        KK = NStates*(K-1)
        I = (IA-1)/NVB+1
        A = IA-NVB*(I-1)+NOB
        DO JB = 1, NOBNVB
          L = JB+NOANVA+1
          J = (JB-1)/NVB+1
          B = JB-NVB*(J-1)+NOB
          If(I.eq.J) Xmat(KK+L) = Xmat(KK+L)+XB(A,B)
          If(A.eq.B) Xmat(KK+L) = Xmat(KK+L)-XB(I,J)
          END DO
        Xmat(KK+K) = Xmat(KK+K) + X00
        END DO
      END SUBROUTINE FormMatU

C      !==========================================================!
C      ! Calculate MO occupancies for a CIS vector                !
C      !==========================================================!
      SUBROUTINE PopCIS(NOA,NVA,N,Psi,Psi1,CISvec,Pop)
      INTEGER I,II,A,AA,K,L,KK,N,NOA,NVA,NRorb,NStates
      COMPLEX*16 Psi(*),Psi1(*)
      REAL*8 CISvec(*),Pop(*),C,Norm,Zero,One,Two
      COMMON IOut
C
      Zero = 0.D0
      One = 1.D0
      Two = 2.D0
      NStates = NOA*NVA + 1
      NRorb = NOA+NVA
C
C     Transform Psi from CIS basis to determinental basis
C
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        end Do
      Do K = 1, N
        KK = NStates*(K-1)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*Psi(K)
          end Do
      end Do
C
C     Normalize
C
      Norm = Zero
      Do K = 1, NStates
        Norm = Norm + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Norm =One/Sqrt(Norm)
      Do K = 1, NStates
        Psi1(K) = Psi1(K)*Norm
        end Do
C
C     Calculate populations in the MO basis
C
      Do I = 1, NOA
        Pop(I) = Two
        end Do
      Do I = NOA+1, NRorb
        Pop(I) = Zero
        end Do
      DO K = 2, NStates
        I = (K-2)/NVA+1
        A = K-1-NVA*(I-1)+NOA
        C = REAL(Psi1(K)*DCONJG(Psi1(K)))
        Pop(A) = Pop(A)+C
        Pop(I) = Pop(I)-C
      end Do
      END SUBROUTINE PopCIS

C      !==========================================================!
C      ! Calculate MO occupancies for a UCIS vector                !
C      !==========================================================!
      SUBROUTINE PopUCIS(NOA,NVA,NOB,NVB,N,Psi,Psi1,CISvec,Pop)
      INTEGER I,II,IA,A,AA,K,L,KK,N,
     $ NOA,NVA,NOB,NVB,NOANVA,NOBNVB,NRorb,NStates
      COMPLEX*16 Psi(*),Psi1(*)
      REAL*8 CISvec(*),Pop(*),C,Norm,Zero,One
C
      Zero = 0.D0
      One = 1.D0
      NRorb = NOA+NVA
      NOANVA = NOA*NVA
      NOBNVB = NOB*NVB
      NStates = NOANVA + NOBNVB + 1
C
C     Transform Psi from CIS basis to determinental basis
C
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        end Do
      Do K = 1, N
        KK = NStates*(K-1)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*Psi(K)
          end Do
      end Do
C
C     Normalize
C
      Norm = Zero
      Do K = 1, NStates
        Norm = Norm + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Norm =One/Sqrt(Norm)
      Do K = 1, NStates
        Psi1(K) = Psi1(K)*Norm
        end Do
C
C     Calculate populations in the MO basis
C
      Do I = 1, 2*NRorb
        Pop(I) = Zero
        end Do
      Do I = 1, NOA
        Pop(I) = One
        end Do
      Do I = 1, NOB
        Pop(I+NRorb) = One
        end Do
C
C     Alpha MOs
C
      DO IA = 1,NOANVA
        K = IA+1
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        C = REAL(Psi1(K)*DCONJG(Psi1(K)))
        Pop(A) = Pop(A)+C
        Pop(I) = Pop(I)-C
      end Do
C
C     Beta MOs
C
      DO IA = 1,NOBNVB
        K = IA+NOANVA+1
        I = (IA-1)/NVB+1
        A = IA-NVB*(I-1)+NOB
        C = REAL(Psi1(K)*DCONJG(Psi1(K)))
        Pop(A+NRorb) = Pop(A+NRorb)+C
        Pop(I+NRorb) = Pop(I+NRorb)-C
      end Do
      END SUBROUTINE PopUCIS

C      !==========================================================!
C      ! Calculate the density for a CIS vector                   !
C      !==========================================================!
      SUBROUTINE DensCIS(NOA,NVA,N,Psi,Psi1,CISvec,Dens)
      INTEGER I,II,J,IJ,A,AA,B,K,L,KK,N,NOA,NVA,NRorb,NStates
      COMPLEX*16 Psi(*),Psi1(*)
      REAL*8 CISvec(*),Dens(*),C,Norm,Zero,One,Two
      COMMON IOut
C
      Zero = 0.D0
      One = 1.D0
      Two = 2.D0
      NStates = NOA*NVA + 1
      NRorb = NOA+NVA
C
C     Transform Psi from CIS basis to determinental basis
C
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        end Do
      Do K = 1, N
        KK = NStates*(K-1)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*Psi(K)
          end Do
      end Do
C
C     Normalize
C
      Norm = Zero
      Do K = 1, NStates
        Norm = Norm + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Norm =One/Sqrt(Norm)
      Do K = 1, NStates
        Psi1(K) = Psi1(K)*Norm
        end Do
C
C     Form density matrix in the MO basis
C
      Do IJ = 1, NRorb*NRorb
        Dens(IJ) = Zero
        end Do
      Do I = 1, NOA
        II = NRorb*(I-1)+I
        Dens(II) = Two
        end Do
C
C     GS to ES contributions
C
      DO K = 2, NStates
        I = (K-2)/NVA+1
        II = NRorb*(I-1)
        A = K-1-NVA*(I-1)+NOA
        AA = NRorb*(A-1)
        C = REAL(Psi1(K)*DCONJG(Psi1(1)))*Sqrt(Two)
        Dens(II+A) = Dens(II+A)+C
        Dens(AA+I) = Dens(AA+I)+C
        end Do
C
C     ES to ES contributions
C
      DO K = 2, NStates
        I = (K-2)/NVA+1
        II = NRorb*(I-1)
        A = K-1-NVA*(I-1)+NOA
        AA = NRorb*(A-1)
        DO L = 2, NStates
          J = (L-2)/NVA+1
          B = L-1-NVA*(J-1)+NOA
          C = REAL(Psi1(K)*DCONJG(Psi1(L)))
          If(I.eq.J) Dens(AA+B) = Dens(AA+B)+C
          If(A.eq.B) Dens(II+J) = Dens(II+J)-C
          END DO
        END DO
C
C     Convert from Square to Lower Triangle
C
      IJ = 0
      Do I = 1,NRorb
        II = NRorb*(I-1)
        Do J = 1,I
        IJ = IJ+1
        Dens(IJ) = Dens(II+J)
        end Do
      end Do
      END SUBROUTINE DensCIS

C      !==========================================================!
C      ! Calculate the density for a UCIS vector                  !
C      !==========================================================!
      SUBROUTINE DensUCIS(NOA,NVA,NOB,NVB,N,Psi,Psi1,CISvec,
     $ DensA,DensB)
      INTEGER I,II,IA,J,IJ,JB,A,AA,B,K,L,KK,
     $ N,NOA,NVA,NOB,NVB,NOANVA,NOBNVB,NRorb,NStates
      COMPLEX*16 Psi(*),Psi1(*)
      REAL*8 CISvec(*),DensA(*),DensB(*),C,Norm,Zero,One
C
      Zero = 0.D0
      One = 1.D0
      NRorb = NOA+NVA
      NOANVA = NOA*NVA
      NOBNVB = NOB*NVB
      NStates = NOANVA + NOBNVB + 1
C
C     Transform Psi from CIS basis to determinental basis
C
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        end Do
      Do K = 1, N
        KK = NStates*(K-1)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*Psi(K)
          end Do
      end Do
C
C     Normalize
C
      Norm = Zero
      Do K = 1, NStates
        Norm = Norm + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Norm =One/Sqrt(Norm)
      Do K = 1, NStates
        Psi1(K) = Psi1(K)*Norm
        end Do
C
C     Form density matrix in the MO basis
C
      Do IJ = 1, NRorb*NRorb
        DensA(IJ) = Zero
        DensB(IJ) = Zero
        end Do
      Do I = 1, NOA
        II = NRorb*(I-1)+I
        DensA(II) = One
        end Do
      Do I = 1, NOB
        II = NRorb*(I-1)+I
        DensB(II) = One
        end Do
C
C     GS to ES contributions for alpha
C
      DO IA = 1, NOANVA
        K = IA+1
        I = (IA-1)/NVA+1
        II = NRorb*(I-1)
        A = IA-NVA*(I-1)+NOA
        AA = NRorb*(A-1)
        C = REAL(Psi1(K)*DCONJG(Psi1(1)))
        DensA(II+A) = DensA(II+A)+C
        DensA(AA+I) = DensA(AA+I)+C
        end Do
C
C     GS to ES contributions for beta
C
      DO IA = 1, NOBNVB
        K = IA+NOANVA+1
        I = (IA-1)/NVB+1
        II = NRorb*(I-1)
        A = IA-NVB*(I-1)+NOB
        AA = NRorb*(A-1)
        C = REAL(Psi1(K)*DCONJG(Psi1(1)))
        DensB(II+A) = DensB(II+A)+C
        DensB(AA+I) = DensB(AA+I)+C
        end Do
C
C     ES to ES contributions for alpha
C
      DO IA = 1, NOANVA
        K = IA+1
        I = (IA-1)/NVA+1
        II = NRorb*(I-1)
        A = IA-NVA*(I-1)+NOA
        AA = NRorb*(A-1)
        DO JB = 1, NOANVA
          L = JB+1
          J = (JB-1)/NVA+1
          B = JB-NVA*(J-1)+NOA
          C = REAL(Psi1(K)*DCONJG(Psi1(L)))
          If(I.eq.J) DensA(AA+B) = DensA(AA+B)+C
          If(A.eq.B) DensA(II+J) = DensA(II+J)-C
          END DO
        END DO
C
C     ES to ES contributions for beta
C
      DO IA = 1, NOBNVB
        K = IA+NOANVA+1
        I = (IA-1)/NVB+1
        II = NRorb*(I-1)
        A = IA-NVB*(I-1)+NOB
        AA = NRorb*(A-1)
        DO JB = 1, NOBNVB
          L = JB+NOANVA+1
          J = (JB-1)/NVB+1
          B = JB-NVB*(J-1)+NOB
          C = REAL(Psi1(K)*DCONJG(Psi1(L)))
          If(I.eq.J) DensB(AA+B) = DensB(AA+B)+C
          If(A.eq.B) DensB(II+J) = DensB(II+J)-C
          END DO
        END DO
C
C     Convert from Square to Lower Triangle
C
      IJ = 0
      Do I = 1,NRorb
        II = NRorb*(I-1)
        Do J = 1,I
        IJ = IJ+1
        DensA(IJ) = DensA(II+J)
        DensB(IJ) = DensB(II+J)
        end Do
      end Do
      END SUBROUTINE DensUCIS

C      !==========================================================!
C      ! Calculate MO occupancies for the absorbed wavefunction   !
C      !==========================================================!
      SUBROUTINE PopIon(NOA,NVA,N,Psi,Psi1,Psi2,ExpABP,CISvec,
     $ Rate,Norm,Pop,PopInt)
      INTEGER I,A,K,L,KK,N,NOA,NVA,NStates
      COMPLEX*16 Psi(*),Psi1(*),Psi2(*),C
      REAL*8 ExpABP(*),CISvec(*),Pop(*),PopInt(*),D,Norm,Norm1,Norm2,
     $ Rate,Zero,One,Two
      COMMON IOut
C
      Zero = 0.D0
      One = 1.D0
      Two = 2.D0
      NStates = NOA*NVA + 1
      NRorb = NOA+NVA
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        Psi2(K) = DCMPLX(Zero,Zero)
        end Do
C
C     Normalize Psi
C
      Norm = Zero
      Do K = 1, N
        Norm = Norm + REAL(Psi(K)*DCONJG(Psi(K)))
        end Do
      Norm = Sqrt(Norm)
C
C     Multiply by the exponential of the absorbing potential
C     Psi2 = Exp(V Dt) Psi
C
      Do K = 1, N
        KK = N*(K-1)
        C = Psi(K)
        Do L = 1, N
          Psi2(L) = Psi2(L) + ExpABP(KK+L)*C
          end Do
        end Do
C
C     Decrease in Norm is related to the instantaneous rate
C
      Norm1 = Zero
      Do K = 1, N
        Norm1 = Norm1 + REAL(Psi2(K)*DCONJG(Psi2(K)))
        end Do
      Norm1 = Sqrt(Norm1)
      Rate = Norm - Norm1 
C
C     Transform Psi-Exp(V Dt)Psi from CIS basis to determinental basis
C
      Do K = 1, N
        KK = NStates*(K-1)
        C = Psi(K) - Psi2(K)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*C
          end Do
      end Do
C
C     Normalize and calculate populations in the MO basis
C
      Norm2 = Zero
      Do K = 1, NStates
        Norm2 = Norm2 + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Do I = 1, NOA
        Pop(I) = Two
        end Do
      Do I = NOA+1, NRorb
        Pop(I) = Zero
        end Do
      DO K = 2, NStates
        I = (K-2)/NVA+1
        A = K-1-NVA*(I-1)+NOA
        D = REAL(Psi1(K)*DCONJG(Psi1(K)))/Norm2
        Pop(A) = Pop(A)+D
        Pop(I) = Pop(I)-D
        end Do
      Do I = 1, NRorb
        PopInt(I) = PopInt(I) + Rate*Pop(I)
        end Do
      END SUBROUTINE PopIon

C      !==========================================================!
C      ! Calculate MO occupancies for the absorbed UCIS wfn       !
C      !==========================================================!
      SUBROUTINE PopIonU(NOA,NVA,NOB,NVB,N,Psi,Psi1,Psi2,
     $ ExpABP,CISvec,Rate,Norm,Pop,PopInt)
      INTEGER I,A,K,L,KK,N,NOA,NVA,NOB,NVB,NOANVA,NOBNVB,
     $ NStates
      COMPLEX*16 Psi(*),Psi1(*),Psi2(*),C
      REAL*8 ExpABP(*),CISvec(*),Pop(*),PopInt(*),
     $ D,Norm,Norm1,Norm2,Rate,Zero,One
C
      Zero = 0.D0
      One = 1.D0
      NRorb = NOA+NVA
      NOANVA = NOA*NVA
      NOBNVB = NOB*NVB
      NStates = NOANVA + NOBNVB + 1
      Do K = 1, NStates
        Psi1(K) = DCMPLX(Zero,Zero)
        Psi2(K) = DCMPLX(Zero,Zero)
        end Do
C
C     Normalize Psi
C
      Norm = Zero
      Do K = 1, N
        Norm = Norm + REAL(Psi(K)*DCONJG(Psi(K)))
        end Do
      Norm = Sqrt(Norm)
C
C     Multiply by the exponential of the absorbing potential
C     Psi2 = Exp(V Dt) Psi
C
      Do K = 1, N
        KK = N*(K-1)
        C = Psi(K)
        Do L = 1, N
          Psi2(L) = Psi2(L) + ExpABP(KK+L)*C
          end Do
        end Do
C
C     Decrease in Norm is related to the instantaneous rate
C
      Norm1 = Zero
      Do K = 1, N
        Norm1 = Norm1 + REAL(Psi2(K)*DCONJG(Psi2(K)))
        end Do
      Norm1 = Sqrt(Norm1)
      Rate = Norm - Norm1 
C
C     Transform Psi-Exp(V Dt)Psi from CIS basis to determinental basis
C
      Do K = 1, N
        KK = NStates*(K-1)
        C = Psi(K) - Psi2(K)
        Do L = 1, NStates
          Psi1(L) = Psi1(L) + CISvec(KK+L)*C
          end Do
      end Do
C
C     Normalize and calculate populations in the MO basis
C
      Norm2 = Zero
      Do K = 1, NStates
        Norm2 = Norm2 + REAL(Psi1(K)*DCONJG(Psi1(K)))
        end Do
      Do I = 1, 2*NRorb
        Pop(I) = Zero
        end Do
      Do I = 1, NOA
        Pop(I) = One
        end Do
      Do I = 1, NOB
        Pop(I+NRorb) = One
        end Do
C
C     Alpha MOs
C
      DO IA = 1,NOANVA
        K = IA+1
        I = (IA-1)/NVA+1
        A = IA-NVA*(I-1)+NOA
        D = REAL(Psi1(K)*DCONJG(Psi1(K)))/Norm2
        Pop(A) = Pop(A)+D
        Pop(I) = Pop(I)-D
      end Do
C
C     Beta MOs
C
      DO IA = 1,NOBNVB
        K = IA+NOANVA+1
        I = (IA-1)/NVB+1
        A = IA-NVB*(I-1)+NOB
        D = REAL(Psi1(K)*DCONJG(Psi1(K)))/Norm2
        Pop(A+NRorb) = Pop(A+NRorb)+D
        Pop(I+NRorb) = Pop(I+NRorb)-D
      end Do
      Do I = 1, 2*NRorb
        PopInt(I) = PopInt(I) + Rate*Pop(I)
        end Do
      END SUBROUTINE PopIonU
 
C      !==========================================================!
C      ! Transform one electron matrix from AO (lower triangle)   !
C      ! to MO (full)                                             !
C      !==========================================================!
      SUBROUTINE TransAOMO(NBasis,NRorb,XAO,XMO,CMO,scratch)
      INTEGER I,J,K,L,NBasis,NRorb
      REAL*8 XAO(*),XMO(NRorb,*),CMO(Nbasis,*),scratch(NRorb,*)
C
      Do L = 1,NRorb
        Do K = 1,NRorb
          XMO(K,L) = 0.D0
          End Do
        End Do
      Do J = 1,NBasis
        Do K = 1,NRorb
          scratch(K,J) = 0.D0
          End Do
        End Do
      Do I = 1,NBasis
        Do J = 1,NBasis
          IJ = I*(I-1)/2+J
          If(J.gt.I) IJ = J*(J-1)/2+I
          Do K = 1,NRorb
            scratch(K,J)=scratch(K,J)+CMO(I,K)*XAO(IJ)
            End Do
          End Do
        End Do
      Do J = 1,NBasis
        Do L = 1,NRorb
          Do K = 1,NRorb
            XMO(K,L)=XMO(K,L)+scratch(K,J)*CMO(J,L)
            End Do
          End Do
        End Do
      END SUBROUTINE TransAOMO

C      !==========================================================!
C      ! Transform matrix X(NStates,NStates) from single det      !
C      ! to CIS states (NstUse*NstUse)                            !
C      !==========================================================!
      SUBROUTINE TransCIS(NStates,NStUse,X,CISvec,scratch)
      INTEGER NStates,NStUse
      REAL*8 X(*),CISvec(*),scratch(*)
      Real*8 Zero,One
C
      Zero = 0.D0
      One = 1.D0
      Call DGEMM('n','n',NStates,NStUse,NStates,One,X,NStates,
     $ CISvec,NStates,Zero,scratch,NStates)
      Call DGEMM('t','n',NstUse,NstUse,NStates,One,CISvec,NStates,
     $ scratch,NStates,Zero,X,NstUse)
C
      END SUBROUTINE TransCIS

C      !==========================================================!
C      ! Transformation B^T A B, with result written back to A    !
C      | A(N*N), B(N*M), V(N), results written back to A(M*M)     !
C      !==========================================================!
      SUBROUTINE Trans1(N,M,A,B,V)
      INTEGER i,ii,j,jj,jjj,k,kk,m,n
      REAL*8 A(*),B(*),V(*),C
C
C  A*B and store in A
C
      Do i=1,N
        Do j=1,N
          jj=(j-1)*N
          V(j)=A(i+jj)
          End Do
        Do k=1,M
          kk=(k-1)*N
          C=0.D0
          Do j=1,N
            C=C+V(j)*B(j+kk)
            End Do
          A(i+kk)=C
          End Do
        End Do
C
C  B^T*(A*B) and store in A
C
      Do k=1,M
        kk=(k-1)*N
        Do j=1,N
          V(j)=A(j+kk)
          End Do
        Do i=1,M
          ii=(i-1)*N
          C=0.D0
          Do j=1,N
            C=C+B(ii+j)*V(j)
            End Do
          A(i+kk)=C
          End Do
        End Do
C
C  Pack down if M < N
C
      If(M.lt.N) then
        Do j=2,M
          jj=(j-1)*M
          jjj=(j-1)*N
          Do i=1,M
            A(i+jj)=A(i+jjj)
            End Do
          End Do
        EndIf
C
      END SUBROUTINE Trans1

C      !==========================================================!
C      ! Norm of a Complex Vector                                 !
C      !==========================================================!
      SUBROUTINE NORMsr(X, NORM, N, RorI)
      INTEGER N, i
      REAL*8 X(2*N), NORM, total
      CHARACTER*1 RorI
      total = 0.0
      IF ( RorI.eq.'R') THEN
           DO i =1, N
             total = total + ( X(i)*X(i) ) 
           END DO
      ELSE IF ( RorI.eq.'I') THEN
           DO i=1, N
             total=total+ ( ( X(2*i)*X(2*i) ) + ( X(2*i-1)* X(2*i-1) ) )
           END DO
      ELSE
          WRITE(*,*) 'Error in Norm'
      END IF
      NORM = DSQRT(total)
      END SUBROUTINE NORMsr

C      !==========================================================!
C      ! Complex Matrix Vector Multiplication Wrapper             !
C      !==========================================================!
      SUBROUTINE MatVec(InMat, InVec, OutVec, N, IT)
      INTEGER N, ONE
      CHARACTER*(*) IT
      COMPLEX*16 InMat(N*N), InVec(N), OutVec(N), iONE, iZERO
C
       ONE = 1
      iONE = DCMPLX(1.0,0.0)
      iZERO= DCMPLX(0.0,0.0)
C      CALL ZGEMM(IT,'n',N,ONE,N,iONE,InMat,N,InVec,N,iZERO,OutVec,N)
      CALL ZGEMV(IT,N,N,iONE,InMat,N,InVec,1,iZERO,OutVec,1)
      END SUBROUTINE MatVec

C      !==========================================================!
C      ! Square Real Matrix times Complex Vector                  !
C      !==========================================================!
      SUBROUTINE MatVec1(InMat, InVec, OutVec, N, IT)
      INTEGER N, I, J, IJ
      CHARACTER*(*) IT
      Real*8 InMat(N*N)
      COMPLEX*16 InVec(N), OutVec(N), A
C
      If(IT.eq.'t'.or.IT.eq.'T') then
        Do J = 1, N
          IJ = N*(J-1)
          A = DCMPLX(0.0D0,0.0D0)
          Do I = 1, N
            IJ = IJ + 1
            A = A + InMat(IJ)*InVec(I)
            End Do
            OutVec(J) = A
          End Do
      Else
        IJ = 0
        Do I = 1, N
          OutVec(I) = DCMPLX(0.0D0,0.0D0)
          Do J = 1, N
            IJ = IJ + 1
            OutVec(J) = OutVec(J) + InMat(IJ)*InVec(I)
            End Do
          End Do
        End If
      END SUBROUTINE MatVec1
C
C      !==========================================================!
C      ! Complex Matrix Multiplication Wrapper                    !
C      !==========================================================!
      SUBROUTINE iMatMult( MatA, AT, MatB, BT, MatC, N, RealorImg)
C
C     MatMult Multiplies the matrices MatA and MatB and stores the
C     product in MatC.  AT and BT are flags which indicate if MatA and
C     MatB repsectively should be transposed or not, indicated by 't' or
C     'n' repsectively
C
      CHARACTER*(*) AT, BT, RealorImg
      INTEGER N
      COMPLEX*16 MatA(N*N), MatB(N*N ), MatC(N*N),D1,D2
C
      D1= 1.0D0
      D2= 0.0D0
C
      CALL ZGEMM(AT, BT, N, N, N,D1,MatA,N,MatB,N,D2,MatC,N)
C
      END SUBROUTINE iMatMult      
C
C      !==========================================================!
C      ! Matrix Multiplication Wrapper                            !
C      !==========================================================!
      SUBROUTINE MatMult( MatA, AT, MatB, BT, MatC, N, RealorImg)
C
C     MatMult Multiplies the matrices MatA and MatB and stores the
C     product in MatC.  AT and BT are flags which indicate if MatA and
C     MatB repsectively should be transposed or not, indicated by 't' or
C     'n' repsectively
C
      CHARACTER*(*) AT, BT, RealorImg
      INTEGER N
      REAL*8 MatA( N, N ), MatB( N, N ), MatC( N, N ),D1,D2
C
      D1= 1.0
      D2= 0.0
C
      CALL DGEMM(AT, BT, N, N, N,D1,MatA,N,MatB,N,D2,MatC,N)
C
      END SUBROUTINE MatMult      
C
C      !==========================================================!
C      ! Matrix Multiplication A(j,i)B(j,k)=C(i,k)                ! 
C      ! C is written back to A                                   !
C      !==========================================================!
      SUBROUTINE MatMult1(A,B,V,N)
C
      INTEGER I,II,J,JJ,K,KK,N
      REAL*8 A(N*N),B(N*N),V(N),C
C
      Do I=1,N
        II=(I-1)*N
        Do J=1,N
          V(J)=A(II+J)
          End Do
        Do K=1,N
          KK=(K-1)*N
          C=0.D0
          Do J=1,N
            C=C+V(J)*B(J+KK)
            End Do
          A(II+K)=C
          End Do
        End Do
      Do I=2,N
        II=(I-1)*N
        Do J=1,I
          JJ=(J-1)*N
          C=A(I+JJ)
          A(I+JJ)=A(II+J)
          A(II+J)=C
          End Do
      End Do
C
      END SUBROUTINE MatMult1
C
C      !==========================================================!
C      ! Matrix Diagonalization Wrapper                           !
C      !==========================================================!
      SUBROUTINE MATDIAG( INMat, VECS, VALS, N, LWMAX, WORK )
C
C     MatDiag Diagonalizes the N by N symmetric matrix INMat and 
C     places the eigenvectors and eigenvalues in VECS and VALS
C     respectively.
C
      INTEGER          LWORK, N, INFO, I, J, LWMAX
      REAL*8 INMat( N, N ), VECS( N, N ), VALS( N ), WORK( * ), 
     $ TEMP( N , N )
C
C     Make backup of original matrix
      do I = 1, N
         do J = 1, N
            TEMP(I,J) = INMat(I,J)
         end do
      end do
      LWORK = -1
      CALL DSYEV( 'Vectors', 'Upper', N, INMat, N, VALS, WORK, LWORK,
     $ INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) ) 
C
C     Solve eigenproblem.
      CALL DSYEV( 'Vectors', 'Upper', N, INMat, N, VALS, WORK, LWORK,
     $ INFO )
C
C     Check for convergence.
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
C
C     Move Eigenvectors to Vecs and Original Matrix Back to INMat
      do I = 1 , N
         do J = 1, N
            VECS(I,J) = INMat(I,J)
            INMat(I,J) = TEMP(I,J)
         end do
      end do
C
      END SUBROUTINE MATDIAG      
C
C      !=============================================================!
C      ! B(i) = Exp(I c*A(i)) * B(i) (c and A are Real, B is complex)!
C      !=============================================================!
      SUBROUTINE ExpVec(c,A,B,N)
      INTEGER N,i
      REAL*8 c,A(*)
      COMPLEX*16 B(*)
C      WRITE(*,*)c,A(1),DCOS(c*A(1)),DSIN(c*A(1)),B(1)
      DO i=1,N
        B(i)=DCMPLX(DCOS(c*A(i)),DSIN(c*A(i)))*B(i)
      END DO
      END SUBROUTINE ExpVec
C
C      !==========================================================!
C      ! Generating the pulse envelope *Cosine Envelope*          !
C      !==========================================================!
      Subroutine CosEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect(*)
      do cnt = 1,maxstep
        if (cnt < tt) then
          evect(cnt)= 0.5D0 - DCOS(2.0D0*Pi*cnt/(ncyc*period))/2.0D0
        else 
          evect(cnt)=0.0D0
        end if
      end do
      do cnt= 1,maxstep
        fvect(cnt)=evect(cnt)*DCOS((omega*delta*cnt)-phase)
      end do  
      End Subroutine CosEnGen

C      !==========================================================!
C      ! Generating the pulse envelope *Gaussian Envelope*        !
C      !==========================================================!
      Subroutine GauEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect(*)
      do cnt = 1,maxstep
        if (cnt <tt) then
          evect(cnt)=(DEXP(-16.0D0*DLOG(2.0D0)*(cnt/(ncyc*period)-0.5D0)
     $**2.0D0)-1.0D0/16.0D0)/(15.0D0/16.0D0)
        else
          evect(cnt)=0.0
        end if
      end do
      do cnt= 1,maxstep
        fvect(cnt)=evect(cnt)*DSIN((omega*delta*cnt)-phase)
      end do  
      End Subroutine GauEnGen


C      !==========================================================!
C      ! Generating the pulse envelope *Trapezoidal Envelope*     !
C      !==========================================================!
      Subroutine TrapEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect(*)
      do cnt = 1,maxstep
        if (cnt <tt) then
          if (cnt < (ncyc-1)*period) then
            if (cnt < period) then
              evect(cnt)=cnt/period
            else
              evect(cnt)=1.0
            end if
          else
            evect(cnt)=(ncyc*period-cnt)/period
          end if
        else
          evect(cnt)=0.0
        end if
      end do
      do cnt= 1,maxstep
        fvect(cnt)=evect(cnt)*DSIN((omega*delta*cnt)-phase)
      end do  
      End Subroutine TrapEnGen
      
C      !==========================================================!
C      ! Generating the pulse envelope *Static Field*             !
C      !==========================================================!
      SUBROUTINE CwEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      INTEGER maxstep,ncyc,cnt,ifield,mxstep10
      REAL*8 Pi,delta,emax,omega,period,phase,tt,evect(*),fvect(*)
      DO cnt = 1, maxstep
         IF ( cnt .le. (2D0*maxstep/3D0) ) THEN
          fvect(cnt)=(1.0D0-(1.0D0-3.0D0*cnt/(2.0D0*maxstep))**4.0D0)
         ELSE
          fvect(cnt)=1.0D0
         END IF
      END DO
c      DO cnt = 1, maxstep
c         fvect(cnt) = emax*evect(cnt)
c      END DO
      END SUBROUTINE CwEnGen
C
C      !==========================================================!
C      ! Generating the Null Pulse                                !
C      !==========================================================!
      SUBROUTINE NullEnGen(maxstep,evect,fvect)
      INTEGER maxstep,cnt
      REAL*8 evect(*),fvect(*)
      DO cnt = 1, maxstep
         evect(cnt) = 0.0D0
         fvect(cnt) = 0.0D0
      END DO
      END SUBROUTINE NullEnGen
C
C      !==========================================================!
C      ! Generating the "Bandrauk pulse"                          !
C      !    *Envelope: linear rise then constant CW*              !
C      !==========================================================!
      Subroutine CsEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect(*)
      do cnt = 1,maxstep
        if (cnt < tt) then
          evect(cnt)= cnt/(ncyc*period)
        else 
          evect(cnt)=1.0D0
        end if
      end do
      do cnt= 1,maxstep
      fvect(cnt)=evect(cnt)*DCOS((omega*delta*cnt)-phase)
      end do  
      End Subroutine CsEnGen
C      !==========================================================!
C      ! Generating the pulse envelope *Sine squared Envelope*    !
C      !==========================================================!
      Subroutine Cos2EnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect(*)
      do cnt = 1,maxstep
        if (cnt < tt) then
          evect(cnt)= DSIN(Pi*cnt/(ncyc*period))**2.0D0
        else 
          evect(cnt)=0.0D0
        end if
      end do
      do cnt= 1,maxstep
        fvect(cnt)=evect(cnt)*DSIN((omega*delta*cnt)-phase)
      end do  
      End Subroutine Cos2EnGen
 
C      !==========================================================!
C      ! Generating the pulse envelope *Sine squared Envelope*    !
C      ! circularly polarized (XPol=1. for right, -1. for left)   !
C      !==========================================================!
      Subroutine CircEnGen(maxstep,ncyc,Pi,emax,delta,omega,period,
     $phase,tt,evect,fvect1,fvect2,XPol,YPol)
      implicit none
      integer :: maxstep,ncyc,cnt
      real*8 :: Pi,delta,emax,omega,period,phase,tt
      real*8 :: evect(*),fvect1(*),fvect2(*),XPol,YPol
      do cnt = 1,maxstep
        if (cnt < tt) then
          evect(cnt)= DSIN(Pi*cnt/(ncyc*period))**2.0D0
        else 
          evect(cnt)=0.0D0
        end if
      end do
      do cnt= 1,maxstep
        fvect1(cnt)=XPol*evect(cnt)*DSIN((omega*delta*cnt)-phase)
        fvect2(cnt)=YPol*evect(cnt)*DCOS((omega*delta*cnt)-phase)
      end do  
      End Subroutine CircEnGen

C
C      !==========================================================!
C      !Matrix Exp                                                !
C      !==========================================================!
      SUBROUTINE MatExp(InMat, EXPMAT, N, Lscratch, scratch)
      INTEGER I, N, Lscratch
      REAL*8 InMat(*), EXPMAT(*), VECS(N*N), VALS(N), VALMAT(N*N), D,
     $ TEMP1(N*N), scratch(*)
  100 Format(4f16.8)
      D = 1.0
C
C     Diagonalize the Matrix, find eigenVECtorS and eigenVALueS
      CALL MATDIAG(InMat, VECS, VALS, N, Lscratch, scratch)
C
C     Exponentiate & Convert the eigenVALueS vector into a diagonal eigenVALue MATrix
      CALL VectoExp(VALS, VALMAT, N, D)
C
C     Multiply the Exp eigenVALueMATrix by the eigenVECtorS      
      CALL MatMult(VECS, 'n', VALMAT, 'n', TEMP1, N, 'R')
      CALL MatMult(TEMP1,'n', VECS, 't', EXPMAT, N, 'R')
C
      END SUBROUTINE MATEXP
C
C      !==========================================================!
C      !Matrix Exp                                                !
C      !==========================================================!
      SUBROUTINE MatExp1(Mat, Eig, ExpMat, Fac, N, LWORK, WORK)
      INTEGER I, J, K, N, Lscratch
      REAL*8 Mat(*), Eig(*), ExpMat(*), WORK(*), Fac, D, Zero
C
C     Diagonalize the Matrix
C
      CALL DSYEV('Vectors','Upper',N,Mat,N,Eig,WORK,LWORK,INFO)
C
C     Form the exponential of the matrix
C
      Zero = 0.D0
      Do I=1,N*N
        ExpMat(I) = Zero
        End Do
      Do J=1,N
        JJ=(J-1)*N
        D=DExp(Fac*Eig(J))
        Do K=1,N
          KK=(K-1)*N
          Do I=1,N
            ExpMat(I+KK)=ExpMat(I+KK)+Mat(I+JJ)*D*Mat(K+JJ)
            End Do
          End Do
        End Do
      END SUBROUTINE MATEXP1
C
C      !==========================================================!
C      ! Convert vector into Exp Diagonal Matrix                  !
C      !==========================================================!
      SUBROUTINE VectoExp(Vec, ExpMat, NStates, delta)
      INTEGER NStates, i
      REAL*8 Vec(*), delta, ExpMat(*)
      DO i=1,NStates*NStates
          ExpMat(i)=0.0
      END DO
      DO i=1, NStates
C        WRITE(IOut,*)exp(Vec(i))
        ExpMat(i+(i-1)*NStates)= DEXP(Vec(i)*delta)
      END DO
      END SUBROUTINE VectoExp
C
C      !==========================================================!
C      ! Retrieve Date and Time                                   !
C      !==========================================================!
      SUBROUTINE DNT(QTIME)
      CHARACTER date*8, time*10, zone*5, month*3, day*3, ampm*2
      INTEGER values(8)
      REAL QTIME
      COMMON IOut
 100  FORMAT(A,X,I2,X,I4)
 200  FORMAT(I2,':',I2,':',I2,':',I3,X,A)
 300  FORMAT(A,X,I2,X,I4,3X,I2,':',I2,':',I2,':',I3,X,A)
      CALL DATE_AND_TIME(date,time,zone,values)
      SELECT CASE (values(2))
         CASE (1)
                 month = 'Jan'
         CASE (2)
                 month = 'Feb'
         CASE (3)
                 month = 'Mar'
         CASE (4)
                 month = 'Apr'
         CASE (5)
                 month = 'May'
         CASE (6)
                 month = 'Jun'
         CASE (7)
                 month = 'Jul'
         CASE (8)
                 month = 'Aug'
         CASE (9)
                 month = 'Sep'
         CASE (10)
                 month = 'Oct'
         CASE (11)
                 month = 'Nov'
         CASE (12)
                 month = 'Dec'
      END SELECT
      IF (values(5) >= 12) THEN
              ampm = "PM"
              values(5) = MOD(values(5),12)
      ELSE IF (values(5) == 0) THEN
              ampm = "AM"
              values(5) = 12
      ELSE
              ampm = "AM"
      END IF
      WRITE(IOut,300)month,values(3),values(1),values(5),values(6),
     $ values(7),values(8),ampm
      CALL CPU_TIME(QTIME)
      END SUBROUTINE DNT
