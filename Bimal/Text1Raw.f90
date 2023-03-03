program textRaw
IMPLICIT REAL*8 (A-H,O-Z)                                                 
      PARAMETER (NMIE=10000, NPL=2*NMIE)
      REAL*8 AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BET1(NPL),BET2(NPL),MRR,MRI,LAM
      AA=0.6D0
      BB=0.2D0
      AA1=0.044D0
      AA2=0.36D0
      BB1=DLOG(1.96D0)*DLOG(0.432D0)
      BB2=DLOG(2.37D0)*DLOG(0.741D0)
      GAM=1D0                                                               
      LAM=0.532D0                                            
      MRR=1.141D0                                               
      MRI=0.0036D0                                            
      NDISTR=6
      NK=100
      N=1000
      NP=40
      R1=0.01D0
      R2=15D0       
      DDELT=1D-7
      WRITE (6,1000) DDELT
 1000 FORMAT ('NUMERICAL ACCURACY =',D8.2)
      CALL SPHER (AA,BB,GAM,LAM,MRR,MRI,R1,R2,N,NP,NDISTR,NK,L1,AL1,AL2,AL3,AL4,BET1,BET2,CEXT,CSCA,AREA,VOL,RVW,RMEAN,AA1,BB1,AA2,BB2,DDELT)
      QE=CEXT/AREA
      LMAX=L1-1
      WRITE (6,1001) AREA,VOL,RVW,RMEAN
      CALL MATR (AL1,AL2,AL3,AL4,BET1,BET2,LMAX)
 1001 FORMAT ('<G> = ',D12.6,'  <V> = ',D12.6,'  Rvw = ',D12.6,'  <R> = ',D12.6)
      STOP                                                         
      END 
 
      SUBROUTINE MATR (A1,A2,A3,A4,B1,B2,LMAX)
      PARAMETER (NMIE=10000, NPL=2*NMIE, NPNA=19)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)
      L1MAX=LMAX+1
      D6=DSQRT(6D0)*0.25D0
      PRINT 1003
 1003 FORMAT(' ',5X,'<',11X,'F11',6X,'     F33',6X,'     F12',6X,'     F34')
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      TB=-DB
      TAA=-DA
      DO 500 I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1)*DFLOAT(L*L-4)
            PL4=1D0/(DFLOAT(L)*DFLOAT(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F33=(F2-F3)*0.5D0
10002    FORMAT (4F15.6)
         PRINT 1004,TB,F11,F33,F12,F34
        WRITE (10,1004) TB,F11,F33,F12,F34
  500 CONTINUE
 1004 FORMAT(' ',F6.2,6F14.5)
      RETURN
      END
             
      SUBROUTINE SPHER (AA,BB,GAM,LAM,MRR,MRII,R1,R2,N,NP,NDISTR, NK,L1,AL1,AL2,AL3,AL4,BET1,BET2,CEXT,CSCA,AREA, VOLUME,RVW,RMEAN,AA1,BB1,AA2,BB2,DDELT)              
      PARAMETER (NGRAD=100000, NMIE=10000, NPL=2*NMIE,  NDRDI=3*NMIE)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 LAM,MRR,MRI,MRII,  PSI(NMIE),HI(NDRDI),RPSI(NDRDI),AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BET1(NPL),&        
BET2(NPL),P1(NPL),P2(NPL),P3(NPL),P4(NPL),DR(NDRDI),DI(NDRDI),AR(NMIE),AI(NMIE),BR(NMIE),BI(NMIE),X(NPL),W(NPL),XX(NPL),YY(NGRAD),WY(NGRAD),&
COEFF1(NMIE),COEFF2(NMIE),COEFF3(NMIE),COEF1(NPL),COEF2(NPL),COEF3(NPL),COEF4(NPL),COEF5(NPL),COEF6(NPL),COEF7(NPL),COEF8(NPL),&              
F11(NPL),F33(NPL),F12(NPL),F34(NPL),PIN(NMIE),TAUN(NMIE)                                              
      COMMON /COEFF/ COEFF1, COEFF2, COEFF3                                  
      COMMON /COEF/ COEF1,COEF2,COEF3,COEF4,COEF5,COEF6,COEF7,COEF8,D6        
      COMMON /P/ P1,P2,P3,P4                                  
                        
      MRI=-MRII
      IF (NDISTR.EQ.3) CALL POWER (AA,BB,R1,R2)
                                                                         
      IF(NK.GT.NPL) PRINT 3336                                              
      IF(NK.GT.NPL) STOP                                                      
 3336 FORMAT ('NK IS GREATER THAN NPL. EXECUTION TERMINATED')              
      CALL GAUSS (NK,0,0,X,W)                                                   
      NNK=N*NK                                                                  
      NNPK=NP*NK
      IF (NDISTR.NE.5) NNPK=0
      NNK=NNK+NNPK
      IF (NNK.GT.NGRAD) PRINT 3334                                             
      IF (NNK.GT.NGRAD) STOP                                                   
 3334 FORMAT ('NNK IS GREATER THAN NGRAD, EXECUTION TERMINATED')               
      PI=3.1415927D0                                               
      WN=2D0*PI/LAM                                                             
   25 RX=R2*WN                                                                  
      M=IDINT(RX+4.05D0*RX**0.33333D0+8D0)                                     
      IF (M.GT.NMIE) PRINT 3335                                                 
      IF (M.GT.NMIE) STOP                                                       
 3335 FORMAT ('TOO MANY MIE COEFFICIENTS. INCREASE NMIE.')                
      DO 27 I=1,M                                                               
         DD=DFLOAT(I)                                                           
         DD1=DD+1D0                                                             
         COEFF1(I)=DD1/DD                                                       
         DD2=(2D0*DD+1D0)                                                       
         COEFF2(I)=DD2                                                          
         COEFF3(I)=0.5D0*DD2/(DD*DD1)                                           
   27 CONTINUE                                                                  
      NG=2*M-1                                                                  
      L1MAX=2*M                                                                 
      CM=1D0/(MRR*MRR+MRI*MRI)                                                  
      RMR=MRR*CM                                                                
      RMI=-MRI*CM                                                               

      IF (NDISTR.NE.5) GO TO 30
      Z1=R1/DFLOAT(NP)                                                      
      Z2=Z1*0.5D0                                                               
      DO 28 I=1,NK                                                              
         XX(I)=Z2*X(I)+Z2                     
   28 CONTINUE                                                                  
      DO 29 J=1,NP                                                              
         J1=J-1                                                                 
         ZJ=Z1*DFLOAT(J1)                                                    
         IJ=J1*NK                                                               
         DO 29 I=1,NK                                                           
            I1=I+IJ                                                             
            YY(I1)=XX(I)+ZJ
            WY(I1)=W(I)*Z2                                                      
   29 CONTINUE                                                                  
   30 Z1=(R2-R1)/DFLOAT(N)                                                      
      Z2=Z1*0.5D0                                                               
      DO 32 I=1,NK                                                              
         XX(I)=Z2*X(I)+Z2                        
   32 CONTINUE                                                                  
      DO 34 J=1,N                                                               
         J1=J-1                                                                 
         ZJ=Z1*DFLOAT(J1)+R1                                                    
         IJ=J1*NK+NNPK                                                         
         DO 34 I=1,NK                                                           
            I1=I+IJ                                                             
            YY(I1)=XX(I)+ZJ  
            WY(I1)=W(I)*Z2                                                      
   34 CONTINUE                                                                  

      CALL DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,REFF,VEFF,AREA,VOLUME,RVW,RMEAN,PI,AA1,BB1,AA2,BB2)
      CEXT=0D0                                                                  
      CSCA=0D0                                                                  
      CALL GAUSS (NG,0,0,X,W)                                                   
      DO 40 I=1,NG                                                              
         F11(I)=0D0                                                             
         F33(I)=0D0                                                             
         F12(I)=0D0                                                             
         F34(I)=0D0                                                             
   40 CONTINUE                                 
                                                                                
      DO 100 I=1,NNK                                                            
          Y=YY(I)                                                               
          ZW=WY(I)                                                              
          RX= Y*WN                                                              
          RXR=MRR*RX                                                            
          RXI=MRI*RX                                                            
          CDD=RXR*RXR+RXI*RXI                                                   
          CD=DSQRT(CDD)                                                         
          DC=DCOS(RX)                                                           
          DS=DSIN(RX)                                                           
          CX=1D0/RX                                                             
          CXR=RXR/CDD                                                           
          CXI=-RXI/CDD                                                          
                                                                                   
          M1=IDINT(RX+4.05D0*RX**0.33333D0+8D0)                                
          M2=M1+2+IDINT(1.2D0*DSQRT(RX))+5                                      
          M3=M2-1                                                               
          IF (M2.GT.NDRDI) PRINT 3338                                           
          IF (M2.GT.NDRDI) STOP                                                 
 3338     FORMAT ('M2.GT.NDRDI. EXECUTION TERMINATED')         
          QMAX=DMAX1(DFLOAT(M1),CD)                                             
          M4=IDINT(6.4D0*QMAX**0.33333D0+QMAX)+8                                
          IF (M4.GT.NDRDI) PRINT 3337                                           
          IF (M4.GT.NDRDI) STOP                                                 
 3337     FORMAT ('M4.GT.NDRDI. EXECUTION TERMINATED')         
          D4=DFLOAT(M4+1)                                                       
          DR(M4)=D4*CXR                                                         
          DI(M4)=D4*CXI                                                         
          HI(1)=DS+DC*CX                                                        
          HI(2)=3.0D0*HI(1)*CX-DC                                               
          PSI(1)=CX*DS-DC                                                       
          RPSI(M2)=RX/DFLOAT(2*M2+1)                                            
          DO 45 J=2,M3                                                          
              J1=M2-J+1                                                         
   45         RPSI(J1)=1D0/(DFLOAT(2*J1+1)*CX-RPSI(J1+1))                       
          DO 46 J=2,M4                                                          
              J1=M4-J+2                                                         
              J2=J1-1                                                           
              DJ=DFLOAT(J1)                                                     
              FR=DJ*CXR                                                         
              FI=DJ*CXI                                                         
              OR=DR(J1)+FR                                                      
              OI=DI(J1)+FI                                                      
              ORI=1D0/(OR*OR+OI*OI)                                             
              DR(J2)=FR-OR*ORI                                                  
   46         DI(J2)=FI+OI*ORI                                                  
          M2=M1-1                                                               
          DO 47 J=2,M2                                                          
              J1=J-1                                                            
              J2=J+1                                                            
   47         HI(J2)=DFLOAT(2*J+1)*HI(J)*CX-HI(J1)                              
          DO 48 J=2,M1                                                          
   48         PSI(J)=RPSI(J)*PSI(J-1)                                           
          PSI1=PSI(1)                                                           
          DR1=DR(1)                                                             
          DI1=DI(1)                                                             
          HI1=HI(1)                                                             
          OR=DR1*RMR-DI1*RMI+CX                                                 
          OR1=OR*PSI1-DS                                                        
          OI=DR1*RMI+DI1*RMR                                                    
          OI1=OI*PSI1                                                           
          OR2=OR*PSI1-OI*HI1-DS                                                 
          OI2=OR*HI1+OI*PSI1-DC                                                 
          OAB=1D0/(OR2*OR2+OI2*OI2)                                             
          AR(1)=(OR1*OR2+OI1*OI2)*OAB                                           
          AI(1)=(OR2*OI1-OR1*OI2)*OAB                                           
          OR=MRR*DR1-MRI*DI1+CX                                                 
          OI=MRR*DI1+MRI*DR1                                                    
          OR1=OR*PSI1-DS                                                        
          OR2=OR*PSI1-OI*HI1-DS                                                 
          OI1=OI*PSI1                                                           
          OI2=OR*HI1+OI*PSI1-DC                                                 
          OAB=1D0/(OR2*OR2+OI2*OI2)                                             
          BR(1)=(OR1*OR2+OI1*OI2)*OAB                                           
          BI(1)=(OR2*OI1-OR1*OI2)*OAB                                           
          DO 50 J=2,M1                                                          
              J1=J-1                                                            
              DJ=DFLOAT(J)*CX                                                   
              PSI1=PSI(J)                                                       
              PSI2=PSI(J1)                                                      
              HI1=HI(J)                                                         
              HI2=HI(J1)                                                        
              DR1=DR(J)                                                         
              DI1=DI(J)                                                         
              OR=DR1*RMR-DI1*RMI+DJ                                             
              OI=DR1*RMI+DI1*RMR                                                
              OR1=OR*PSI1-PSI2                                                  
              OI1=OI*PSI1                                                       
              OR2=OR*PSI1-OI*HI1-PSI2                                           
              OI2=OR*HI1+OI*PSI1-HI2                                            
              OAB=1D0/(OR2*OR2+OI2*OI2)                                         
              YAR=(OR1*OR2+OI1*OI2)*OAB                                         
              YAI=(OR2*OI1-OR1*OI2)*OAB                                         
              AR(J)=YAR                                                         
              AI(J)=YAI                                                         
              OR=MRR*DR1-MRI*DI1+DJ                                             
              OI=MRR*DI1+MRI*DR1                                                
              OR1=OR*PSI1-PSI2                                                  
              OR2=OR*PSI1-OI*HI1-PSI2                                           
              OI1=OI*PSI1                                                       
              OI2=OR*HI1+OI*PSI1-HI2                                            
              OAB=1D0/(OR2*OR2+OI2*OI2)                                         
              YBR=(OR1*OR2+OI1*OI2)*OAB                                         
              YBI=(OR2*OI1-OR1*OI2)*OAB                                         
              BR(J)=YBR                                                         
              BI(J)=YBI                                                         
              YAR=YAR*YAR+YAI*YAI+YBR*YBR+YBI*YBI                               
   50     CONTINUE                                                              
                                                                                        
   55     CE=0D0                                                                
          CS=0D0                                                                
          DO 60 J=1,M1                                                          
              CJ=COEFF2(J)                                                      
              ARJ=AR(J)                                                         
              AIJ=AI(J)                                                         
              BRJ=BR(J)                                                         
              BIJ=BI(J)                                                         
              CDA=ARJ*ARJ+AIJ*AIJ                                               
              CDB=BRJ*BRJ+BIJ*BIJ                                               
              CE=CE+CJ*(ARJ+BRJ)                                                
              CS=CS+CJ*(CDA+CDB)                                                
              CJ=COEFF3(J)                                                      
              AR(J)=CJ*(ARJ+BRJ)                                                
              AI(J)=CJ*(AIJ+BIJ)                                                
              BR(J)=CJ*(ARJ-BRJ)                                                
              BI(J)=CJ*(AIJ-BIJ)                                                
   60     CONTINUE                                                              
          CEXT=CEXT+ZW*CE                                                       
          CSCA=CSCA+ZW*CS                                                       
          DO 70 K=1,NG                                                          
              CALL ANGL (M1,X(K),PIN,TAUN)                                      
              SPR=0D0                                                           
              SPI=0D0                                                           
              SMR=0D0                                                           
              SMI=0D0                                                           
              DO 65 J=1,M1                                                      
                  PJ=PIN(J)                                                     
                  TJ=TAUN(J)                                                    
                  PP=TJ+PJ                                                      
                  PM=TJ-PJ                                                      
                  SPR=SPR+AR(J)*PP                                              
                  SPI=SPI+AI(J)*PP                                              
                  SMR=SMR+BR(J)*PM                                              
                  SMI=SMI+BI(J)*PM                                              
   65         CONTINUE                                                          
              D1=(SPR*SPR+SPI*SPI)*ZW                                           
              D2=(SMR*SMR+SMI*SMI)*ZW                                           
              F11(K)=F11(K)+D1+D2                                               
              F33(K)=F33(K)+D1-D2                                               
              F12(K)=F12(K)+(SPR*SMR+SPI*SMI)*ZW*2D0                            
              F34(K)=F34(K)+(SPR*SMI-SPI*SMR)*ZW*2D0                            
   70     CONTINUE                                                              
  100 CONTINUE                                                                  
                                                                                     
      DD=2D0/CSCA                                                               
      DO 120 I=1,NG                                                             
          F11(I)=F11(I)*DD                                                      
          F12(I)=F12(I)*DD                                                      
          F33(I)=F33(I)*DD                                                      
          F34(I)=F34(I)*DD                                                      
  120 CONTINUE                                                                            
                                                                                
      VOL=2D0*PI/(WN*WN)                                                        
      CEXT=CEXT*VOL                                                             
      CSCA=CSCA*VOL                                                             
      ALB=CSCA/CEXT                                                                        
                                                                                
      DO 150 L1=3,L1MAX                                                         
          L=L1-1                                                                
          COEF1(L1)=1D0/DFLOAT(L+1)                                             
          COEF2(L1)=DFLOAT(2*L+1)                                               
          COEF3(L1)=1D0/DSQRT(DFLOAT((L+1)*(L+1)-4))                            
          COEF4(L1)=DSQRT(DFLOAT(L*L-4))                                        
          COEF5(L1)=1D0/(DFLOAT(L)*DFLOAT((L+1)*(L+1)-4))                       
          COEF6(L1)=DFLOAT(2*L+1)*DFLOAT(L*(L+1))                               
          COEF7(L1)=DFLOAT((2*L+1)*4)                                           
          COEF8(L1)=DFLOAT(L+1)*DFLOAT(L*L-4)                                   
  150 CONTINUE
      DO 160 L1=1,L1MAX
          AL1(L1)=0D0                                                           
          AL2(L1)=0D0                                                           
          AL3(L1)=0D0                                                           
          AL4(L1)=0D0                                                           
          BET1(L1)=0D0                                                          
          BET2(L1)=0D0                                                          
  160 CONTINUE                                                                  
      D6=0.25D0*DSQRT(6D0)                                                      
      DO 300 I=1,NG                                                             
          CALL GENER (X(I),L1MAX)                                               
          WI=W(I)                                                               
          FF11=F11(I)*WI                                                        
          FF33=F33(I)*WI                                                        
          FF12=F12(I)*WI                                                        
          FF34=F34(I)*WI                                                        
          FP=FF11+FF33                                                          
          FM=FF11-FF33                                                          
          DO 260 L1=1,L1MAX                                                     
              P1L1=P1(L1)                                                       
              P4L1=P4(L1)                                                       
              AL1(L1)=AL1(L1)+FF11*P1L1                                         
              AL4(L1)=AL4(L1)+FF33*P1L1                                         
              AL2(L1)=AL2(L1)+FP*P2(L1)                                         
              AL3(L1)=AL3(L1)+FM*P3(L1)                                         
              BET1(L1)=BET1(L1)+FF12*P4L1                                       
              BET2(L1)=BET2(L1)+FF34*P4L1                                       
  260     CONTINUE                                                              
  300 CONTINUE                                                                  
      DO 350 L1=1,L1MAX                                                         
          CL=DFLOAT(L1-1)+0.5D0                                                 
          L=L1                                                                  
          AL1(L1)=AL1(L1)*CL                                                    
          A2=AL2(L1)*CL*0.5D0                                                   
          A3=AL3(L1)*CL*0.5D0                                                   
          AL2(L1)=A2+A3                                                         
          AL3(L1)=A2-A3                                                         
          AL4(L1)=AL4(L1)*CL                                                    
          BET1(L1)=BET1(L1)*CL                                                  
          BET2(L1)=-BET2(L1)*CL                                                  
          IF (DABS(AL1(L1)).LE.DDELT) GO TO 400                     
  350 CONTINUE                                                                  
  400 L1MAX=L                                                                   
                                                                                     
  450 PRINT 504,R1,R2                                                           
  504 FORMAT ('R1 =',d12.6,'   R2=',d12.6)                                    
      PRINT 515, REFF,VEFF                                              
  515 FORMAT ('REFF=',d12.6,'   VEFF=',d12.6)      
      PRINT 500,LAM,MRR,MRII                                     
  500 FORMAT('LAM =',F 8.4,'  MRR=',F7.3,'  MRI=',F8.5)            
  505 FORMAT (' X1 =',d12.6,'   X2=',d12.6)                                    
      PRINT 506,NK,N,NP   
  506 FORMAT('NK=',I4,'   N=',I4,'   NP=',I4)  
      Q1=AL1(2)/3D0                                                             
      PRINT 511,Q1                                                         
  511 FORMAT('<COS> =',d12.6)                                
      PRINT 508,CEXT,CSCA,ALB                                                  
  508 FORMAT('CEXT=',d12.6,'   CSCA=',d12.6,'   ALBEDO =',d12.6)    
      PRINT 555,M                                                               
      PRINT 550                                                                 
      PRINT 510                                                                 
  510 FORMAT('*********   EXPANSION COEFFICIENTS   *********')                  
  570 FORMAT('   S     ALPHA 1    ALPHA 2    ALPHA 3',  '    ALPHA 4     BETA 1     BETA 2')                       
      PRINT 570                                                                 
      DO 520 L=1,L1MAX                                                          
         LL=L-1                                                                 
         PRINT 728,LL,AL1(L),AL2(L),AL3(L),AL4(L),BET1(L),BET2(L)               
  520 CONTINUE                                                                  
  728 FORMAT(I4,1X,6F11.5)
  550 FORMAT(' ')                                                               
  555 FORMAT('MAXIMAL ORDER OF MIE COEFFICIENTS = ',I4)                       
      WRITE (10,580) ALB,L1MAX
      DO L=1,L1MAX
         WRITE (10,575) AL1(L),AL2(L),AL3(L),AL4(L),BET1(L),BET2(L)           
      ENDDO   
  575 FORMAT(6D15.7)
  580 FORMAT(D14.8,I8)
      RETURN                                                                    
      END                                                                       
                                                                                                                                                         
      SUBROUTINE ANGL (NMAX,U,PIN,TAUN)                                         
      PARAMETER (NMIE=10000)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 PIN(NMAX),TAUN(NMAX),COEFF1(NMIE)                                  
      COMMON /COEFF/ COEFF1                                                     
      P1=0D0                                                                    
      P2=1D0                                                                    
      N=1                                                                       
    5 S=U*P2                                                                    
      T=S-P1                                                                    
      TAUN(N)=DFLOAT(N)*T-P1                                                    
      PIN(N)=P2                                                                 
      P1=P2                                                                     
      P2=S+COEFF1(N)*T                                                          
      N=N+1                                                                     
      IF(N.LE.NMAX) GO TO 5                                                     
      RETURN                                                                    
      END                                                                                
                                                                                
      SUBROUTINE GENER (U,L1MAX)                                                
      PARAMETER (NMIE=10000, NPL=2*NMIE)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 P1(NPL),P2(NPL),P3(NPL),P4(NPL),  COEF1(NPL),COEF2(NPL),COEF3(NPL),COEF4(NPL), COEF5(NPL),COEF6(NPL),COEF7(NPL),COEF8(NPL)                
      COMMON /COEF/ COEF1,COEF2,COEF3,COEF4,COEF5,COEF6,COEF7,COEF8,D6          
      COMMON /P/ P1,P2,P3,P4                                                    
      DUP=1D0+U                                                                 
      DUM=1D0-U                                                                 
      DU=U*U                                                                    
      P1(1)=1D0                                                                 
      P1(2)=U                                                                   
      P1(3)=0.5D0*(3D0*DU-1D0)                                                  
      P2(1)=0D0                                                                 
      P2(2)=0D0                                                                 
      P2(3)=0.25D0*DUP*DUP                                                      
      P3(1)=0D0                                                                 
      P3(2)=0D0                                                                 
      P3(3)=0.25D0*DUM*DUM                                                      
      P4(1)=0D0                                                                 
      P4(2)=0D0                                                                 
      P4(3)=D6*(DU-1D0)                                                         
      LMAX=L1MAX-1                                                              
      DO 100 L1=3,LMAX                                                          
         C1=COEF1(L1)                                                           
         C2=COEF2(L1)                                                           
         C3=COEF3(L1)                                                           
         C4=COEF4(L1)                                                           
         C5=COEF5(L1)                                                           
         C6=COEF6(L1)                                                           
         C7=COEF7(L1)                                                           
         C8=COEF8(L1)                                                           
         CU1=C2*U                                                               
         CU2=C6*U                                                               
         L2=L1+1                                                                
         L3=L1-1                                                                
         DL=DFLOAT(L3)                                                          
         P1(L2)=C1*(CU1*P1(L1)-DL*P1(L3))                                       
         P2(L2)=C5*((CU2-C7)*P2(L1)-C8*P2(L3))                                  
         P3(L2)=C5*((CU2+C7)*P3(L1)-C8*P3(L3))                                  
         P4(L2)=C3*(CU1*P4(L1)-C4*P4(L3))                                       
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
                                                                                
      SUBROUTINE DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,REFF,VEFF,AREA,VOLUME,RVW,RMEAN,PI,AA1,BB1,AA2,BB2)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 YY(NNK),WY(NNK)                                                    
      IF (NDISTR.EQ.2) GO TO 100                                                
      IF (NDISTR.EQ.3) GO TO 200                                                
      IF (NDISTR.EQ.4) GO TO 300                                                
      IF (NDISTR.EQ.5) GO TO 360                                                
      IF (NDISTR.EQ.6) GO TO 380                                                                                                     
 1001 FORMAT('MODIFIED GAMMA DISTRIBUTION, ALPHA=',F6.4,'  r_c=', F6.4,'  GAMMA=',F6.4)                                                   
      A2=AA/GAM                                                                 
      DB=1D0/BB
      DO 50 I=1,NNK                                                             
         X=YY(I)                                                             
         Y=X**AA                                                                
         X=X*DB
         Y=Y*DEXP(-A2*(X**GAM))                                                 
         WY(I)=WY(I)*Y                                                       
   50 CONTINUE                                                                  
      GO TO 400                                                                 
  100 PRINT 1002,AA,BB                                                          
 1002 FORMAT('LOG NORMAL DISTRIBUTION, r_g=',F8.4,'  [ln(sigma_g)]**2=', F6.4)                                                              
      DA=1D0/AA                                                                 
      DO 150 I=1,NNK                                                            
         X=YY(I)                                                                
         Y=DLOG(X*DA)                                                          
         Y=DEXP(-Y*Y*0.5D0/BB)/X                                             
         WY(I)=WY(I)*Y                                                          
  150 CONTINUE                                                                  
      GO TO 400                                                                 
  200 PRINT 1003                                                                
 1003 FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')                 
      DO 250 I=1,NNK                                                            
         X=YY(I)                                                                
         WY(I)=WY(I)/(X*X*X)                                                 
  250 CONTINUE                                                                  
      GO TO 400                                                                 
  300 PRINT 1004,AA,BB     
 1004 FORMAT ('GAMMA DISTRIBUTION,  a=',F8.4,'  b=',F6.4)
      B2=(1D0-3D0*BB)/BB                                                        
      DAB=1D0/(AA*BB)                                                          
      DO 350 I=1,NNK                                                            
         X=YY(I)                                                                
         X=(X**B2)*DEXP(-X*DAB)                                                 
         WY(I)=WY(I)*X                                                       
  350 CONTINUE                                                                  
      GO TO 400                                                                 
  360 PRINT 1005,BB                                                             
 1005 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  ALPHA=',D10.4)
      DO 370 I=1,NNK                                                            
         X=YY(I)                                                                
         IF (X.LE.R1) WY(I)=WY(I)
         IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE                                                                  
      GO TO 400                                                                 
  380 PRINT 1006                                    
      PRINT 1007,AA1,BB1                                    
      PRINT 1008,AA2,BB2,GAM                                    
 1006 FORMAT('BIMODAL VOLUME LOG NORMAL DISTRIBUTION')
 1007 FORMAT('r_g1=',F8.4,'  [ln(sigma_g1)]**2=',F6.4)
 1008 FORMAT('r_g2=',F8.4,'  [ln(sigma_g2)]**2=',F6.4,' gamma=',F7.3)
      DA1=1D0/AA1       
      DA2=1D0/AA2                  
      DO I=1,NNK                                                            
         X=YY(I)                                                                
         Y1=DLOG(X*DA1)                                                   
         Y2=DLOG(X*DA2)                                           
         Y1=DEXP(-Y1*Y1*0.5D0/BB1)                                             
         Y2=DEXP(-Y2*Y2*0.5D0/BB2)                                             
         WY(I)=WY(I)*(Y1+GAM*Y2)/(X*X*X*X)                        
      ENDDO                                                                     
  400 CONTINUE                                                                  
      SUM=0D0
      DO 450 I=1,NNK
         SUM=SUM+WY(I)
  450 CONTINUE
      SUM=1D0/SUM
      DO 500 I=1,NNK
         WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0D0
      DO 550 I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
  550 CONTINUE
      REFF=0D0
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      REFF=REFF/G
      VEFF=0D0
      VOLUME=0D0
      RVW=0D0
      RMEAN=0D0
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
         VOLUME=VOLUME+X*X*X*WY(I)
         RVW=RVW+X*X*X*X*WY(I)
         RMEAN=RMEAN+X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)
      AREA=G*PI  
      RVW=RVW/VOLUME
      VOLUME=VOLUME*4D0*PI/3D0
      RETURN                                                                    
      END                                                  
                                                                                
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.check*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA', ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
    PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
 
      SUBROUTINE POWER (A,B,R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      COMMON AA,BB
      AA=A
      BB=B
      AX=1D-12
      BX=A  
      R1=ZEROIN (AX,BX,F,0D0)
      R2=(1D0+B)*2D0*A-R1
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION F(R1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON A,B
      R2=(1D0+B)*2D0*A-R1
      F=(R2-R1)/DLOG(R2/R1)-A
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION ZEROIN (AX,BX,F,TOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
   15 A=AX
      B=BX
      FA=F(A)
      FB=F(B)
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
   35 A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (DABS(XM).LE.TOL1) GO TO 90
   44 IF (FB.EQ.0D0) GO TO 90
   45 IF (DABS(E).LT.TOL1) GO TO 70
   46 IF (DABS(FA).LE.DABS(FB)) GO TO 70
   47 IF (A.NE.C) GO TO 50
   48 S=FB/FA
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=DABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70
   64 IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70
   65 E=D
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (DABS(D).GT.TOL1) B=B+D
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB=F(B)
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
   85 GO TO 30
   90 ZEROIN=B
      RETURN
      END