SUBROUTINE SPACE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO SOLVE FOR THE SMALLEST EIGENVALUES-- ASSUMED .GT. 0 --  . 
! .        AND CORRESPONDING EIGENVECTORS IN THE GENERALIZED          . 
! .        EIGENPROBLEM USING THE SUBSPACE ITERATION METHOD           . 
! .                                                                   . 
! .  - - INPUT VARIABLES - -                                          . 
! .        A(NWK)    = STIFFNESS MATRIX IN COMPACTED FORM (ASSUMED    . 
! .                    POSITIVE DEFINITE)                             . 
! .        B(NWM)    = MASS MATRIX IN COMPACTED FORM                  . 
! .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        . 
! .                    ELEMENTS OF STIFFNESS MATRIX A                 . 
! .        R(NN,NC)  = STORAGE FOR EIGENVECTORS                       . 
! .        EIGV(NC)  = STORAGE FOR EIGENVALUES                        . 
! .        TT(NN)    = WORKING VECTOR                                 . 
! .        W(NN)     = WORKING VECTOR                                 . 
! .        AR(NNC)   = WORKING MATRIX STORING PROJECTION OF K         . 
! .        BR(NNC)   = WORKING MATRIX STORING PROJECTION OF M         . 
! .        VEC(NC,NC)= WORKING MATRIX                                 . 
! .        D(NC)     = WORKING VECTOR                                 . 
! .        RTOLV(NC) = WORKING VECTOR                                 . 
! .        BUP(NC)   = WORKING VECTOR                                 . 
! .        BLO(NC)   = WORKING VECTOR                                 . 
! .        BUPC(NC)  = WORKING VECTOR                                 . 
! .        NN        = ORDER OF STIFFNESS AND MASS MATRICES           . 
! .        NNM       = NN + 1                                         . 
! .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF            . 
! .                    STIFFNESS MATRIX                               . 
! .        NWM       = NUMBER OF ELEMENTS BELOW SKYLINE OF            . 
! .                    MASS MATRIX                                    . 
! .                      I. E. NWM=NWK FOR CONSISTENT MASS MATRIX     . 
! .                            NWM=NN  FOR LUMPED MASS MATRIX         . 
! .        NROOT     = NUMBER OF REQUIRED EIGENVALUES AND EIGENVECTORS. 
! .        RTOL      = CONVERGENCE TOLERANCE ON EIGENVALUES           . 
! .                    ( 1.E-06 OR SMALLER )                          . 
! .        NC        = NUMBER OF ITERATION VECTORS USED               . 
! .                    (USUALLY SET TO MIN(2*NROOT, NROOT+8), BUT NC  . 
! .                    CANNOT BE LARGER THAN THE NUMBER OF MASS       . 
! .                    DEGREES OF FREEDOM)                            . 
! .        NNC       = NC*(NC+1)/2 DIMENSION OF STORAGE VECTORS AR,BR . 
! .        NITEM     = MAXIMUM NUMBER OF SUBSPACE ITERATIONS PERMITTED. 
! .                    (USUALLY SET TO 16)                            . 
! .                    THE PARAMETERS NC AND/OR NITEM MUST BE         . 
! .                    INCREASED IF A SOLUTION HAS NOT CONVERGED      . 
! .        IFSS      = FLAG FOR STURM SEQUENCE CHECK                  . 
! .                      EQ.0  NO CHECK                               . 
! .                      EQ.1  CHECK                                  . 
! .        IFPR      = FLAG FOR PRINTING DURING ITERATION             . 
! .                      EQ.0  NO PRINTING                            . 
! .                      EQ.1  PRINT                                  . 
! .        NSTIF     = SCRATCH FILE                                   . 
! .        IOUT      = UNIT USED FOR OUTPUT                           . 
! .                                                                   . 
! .  - - OUTPUT - -                                                   . 
! .        EIGV(NROOT) = EIGENVALUES                                  . 
! .        R(NN,NROOT) = EIGENVECTORS                                 . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . 
! .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . 
! .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . 
! .   SINGLE PRECISION ARITHMETIC.                                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      INTEGER MAXA(NNM)   
      DIMENSION A(NWK),B(NWM),R(NN,NC),TT(NN),W(NN),EIGV(NC),D(NC)
      DIMENSION VEC(NC,NC),AR(NNC),BR(NNC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC)  
      integer nprint

      nprint = 1
!                                                                       
!     SET TOLERANCE FOR JACOBI ITERATION                               
      TOLJ=1.0D-12    
      !write(*,*) "Passei #1" 
      !write(IOUT,*) "Passei #1"  
      write(IOUT,*) "NWK=",NWK 
      write(IOUT,*) "NWM=",NWM
      write(IOUT,*) "NNM=",NNM
      write(IOUT,*) "NN=",NN
      write(IOUT,*) "NC=",NC 
      write(IOUT,*) "NROOT=",NROOT 
      write(IOUT,*) "NITEM=",NITEM         
!
!
      R=0.0D0
      EIGV=0.0D0
      TT=0.0D0
      W=0.0D0 
      AR=0.0D0
      BR=0.0D0
      VEC=0.0D0
      D=0.0D0
      RTOLV=0.0D0
      BUP=0.0D0
      BLO=0.0D0
      BUPC=0.0D0      
!                                                                       
!     INITIALIZATION                                                    
!                                                                     
      ICONV=0                                                           
      NSCH=0                                                            
      NSMAX=12                                                          
      N1=NC + 1                                                         
      NC1=NC - 1                                                        
      REWIND NSTIF                                                      
      WRITE (NSTIF) A                                                   
      do I=1,NC                                                       
      D(I)=0.0D0
      enddo                                                           
!                                                                       
!     ESTABLISH STARTING ITERATION VECTORS                              
!                                                                       
      ND=NN/NC                                                          
      IF (NWM.GT.NN) GO TO 4                                            
      J=0                                                               
      do I=1,NN                                                       
      II=MAXA(I)                                                        
      R(I,1)=B(I)                                                       
      IF (B(I).GT.0) J=J + 1                                            
      W(I)=B(I)/A(II)
      enddo                                                   
      IF (NC.LE.J) GO TO 16                                             
      WRITE (IOUT,1007)                                                 
      GO TO 800                                                         
    4 do I=1,NN                                                      
      II=MAXA(I)                                                        
      R(I,1)=B(II)                                                      
      W(I)=B(II)/A(II)
      enddo                                                  
   16 do J=2,NC                                                      
      do I=1,NN                                                      
      R(I,J)=0.
      enddo
      enddo                                                         
!                                                                       
      L=NN - ND                                                         
      do J=2,NC                                                      
      RT=0.                                                             
      do I=1,L                                                       
      IF (W(I).LT.RT) GO TO 40                                          
      RT=W(I)                                                           
      IJ=I                                                              
   40 CONTINUE
      enddo                                                          
      do I=L,NN                                                      
      IF (W(I).LE.RT) GO TO 50                                          
      RT=W(I)                                                           
      IJ=I                                                              
   50 CONTINUE
      enddo                                                          
      TT(J)=DFLOAT(IJ)                                                   
      W(IJ)=0.0D0                                                          
      L=L - ND                                                          
      R(IJ,J)=1.0D0
      enddo                                                        
!                                                                       
      WRITE (IOUT,1008)                                                 
      WRITE (IOUT,1002) (TT(J),J=2,NC)                                  
!                                                                       
!     A RANDOM VECTOR IS ADDED TO THE LAST VECTOR                       
!                                                                       
      PI=3.141592654D0                                                  
      XX=0.5D0                                                          
      do K=1,NN                                                      
      XX=(PI + XX)**5.0D0                                                   
      IX=INT(XX)                                                        
      XX=XX - DFLOAT(IX)
!     eliminando XX                                                 
      R(K,NC)=R(K,NC) + XX*0.0D0
      enddo    
!
!                                            
!                                                                       
!     FACTORIZE MATRIX A INTO (L)*(D)*(L(T))                            
!                                                                       

        if (nprint == 1) then
        write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
            "75.00%"
        endif

      ISH=0                                                             
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)                                  

        if (nprint == 1) then
        write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
            "85.00%"
        endif

!                                                                       
! - - - S T A R T   O F   I T E R A T I O N   L O O P                   
!                                                                       
      NITE=0                                                            
      TOLJ2=1.0D-24                                                     
  100 NITE=NITE + 1
      !write(*,*) "Passei #1"   
      !write(IOUT,*) "Passei #1"                                                    
      IF (IFPR.EQ.0) GO TO 90                                           
      WRITE (IOUT,1010) NITE                                            
!                                                                       
!     CALCULATE THE PROJECTIONS OF A AND B                              
!                                                                       
   90 IJ=0                                                              
      do J=1,NC                                                     
      do K=1,NN                                                     
      TT(K)=R(K,J)
      enddo                                                      
      CALL REDBAK (A,TT,MAXA,NN)                                        
      do I=J,NC                                                     
      ART=0.                                                            
      do K=1,NN                                                     
      ART=ART + R(K,I)*TT(K)
      enddo                                            
      IJ=IJ + 1                                                         
      AR(IJ)=ART
      enddo                                                        
      do K=1,NN                                                     
      R(K,J)=TT(K)
      enddo                                                      
!  110 CONTINUE
      enddo                                                          
      IJ=0                                                              
      do J=1,NC                                                     
      CALL MULT (TT,B,R(1,J),MAXA,NN,NWM)                               
      do I=J,NC                                                     
      BRT=0.0D0                                                            
      do K=1,NN                                                     
      BRT=BRT + R(K,I)*TT(K)
      enddo                                            
      IJ=IJ + 1                                                         
      BR(IJ)=BRT
      enddo                                                        
      IF (ICONV.GT.0) GO TO 160                                         
      do K=1,NN                                                     
      R(K,J)=TT(K)
      enddo                                                      
  160 CONTINUE
      enddo                                                          
!                                                                       
!     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS                       
!                                                                       
      IF (IFPR.EQ.0) GO TO 320                                          
      IND=1                                                             
  210 WRITE (IOUT,1020)                                                 
      II=1                                                              
      do I=1,NC                                                     
      ITEMP=II + NC - I                                                 
      WRITE (IOUT,1005) (AR(J),J=II,ITEMP)                              
      II=II + N1 - I
      enddo                                                    
      WRITE (IOUT,1030)                                                 
      II=1                                                              
      do I=1,NC                                                     
      ITEMP=II + NC - I                                                 
      WRITE (IOUT,1005) (BR(J),J=II,ITEMP)                              
      II=II + N1 - I 
      enddo                                                   
      IF (IND.EQ.2) GO TO 350                                           
!                                                                       
  320 CALL JACOBI (AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,IOUT)
       
!                                                                       
      IF (IFPR.EQ.0) GO TO 350                                          
      WRITE (IOUT,1040)                                                 
      IND=2                                                             
      GO TO 210                                                         
!                                                                       
!     ARRANGE EIGENVALUES IN ASCENDING ORDER                            
!                                                                       
  350 IS=0                                                              

        if (nprint == 1) then
        write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
            "90.00%"
        endif

      II=1                                                              
      do I=1,NC1                                                    
      ITEMP=II + N1 - I                                                 
      IF (EIGV(I+1).GE.EIGV(I)) GO TO 360                               
      IS=IS + 1                                                         
      EIGVT=EIGV(I+1)                                                   
      EIGV(I+1)=EIGV(I)                                                 
      EIGV(I)=EIGVT                                                     
      BT=BR(ITEMP)                                                      
      BR(ITEMP)=BR(II)                                                  
      BR(II)=BT                                                         
      do K=1,NC                                                     
      RT=VEC(K,I+1)                                                     
      VEC(K,I+1)=VEC(K,I)                                               
  370 VEC(K,I)=RT
      enddo                                                       
  360 II=ITEMP
      enddo                                                          
      IF (IS.GT.0) GO TO 350                                            
      IF (IFPR.EQ.0) GO TO 375                                          
      WRITE (IOUT,1035)                                                 
      WRITE (IOUT,1006) (EIGV(I),I=1,NC)                                
!                                                                       
!     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)           
!        OR     FINAL EIGENVECTOR APPROXIMATIONS (ICONV.GT.0)           
!                                                                       
  375 do I=1,NN                                                     
      do J=1,NC                                                     
      TT(J)=R(I,J)
      enddo                                                      
      do K=1,NC                                                     
      RT=0.                                                             
      do L=1,NC                                                     
      RT=RT + TT(L)*VEC(L,K)
      enddo                                            
  424 R(I,K)=RT
      enddo 
      enddo                                                        
!  420 CONTINUE
!                                                                       
!     CALCULATE ERROR BOUNDS AND CHECK FOR CONVERGENCE OF EIGENVALUES   
!                                                                       
      do I=1,NC                                                     
      VDOT=0.0D0                                                           
      do J=1,NC                                                     
      VDOT=VDOT + VEC(I,J)*VEC(I,J)
      enddo                                     
      EIGV2=EIGV(I)*EIGV(I)                                             
      DIF=VDOT - EIGV2                                                  
      RDIF=MAX(DIF,TOLJ2*EIGV2)/EIGV2                                   
      RDIF=DSQRT(RDIF)                                                   
      RTOLV(I)=RDIF                                                     
!  380 CONTINUE 
      enddo                                                         
      IF (IFPR.EQ.0 .AND. ICONV.EQ.0) GO TO 385                         
      WRITE (IOUT,1050)                                                 
      WRITE (IOUT,1005) (RTOLV(I),I=1,NC)                               
  385 IF (ICONV.GT.0) GO TO 500                                         
!                                                                       
      do I=1,NROOT                                                  
      IF (RTOLV(I).GT.RTOL) GO TO 400                                   
 ! 390 CONTINUE
      enddo                                                          
      WRITE (IOUT,1060) RTOL                                            
      ICONV=1                                                           
      GO TO 100                                                         
  400 IF (NITE.LT.NITEM) GO TO 100                                      
      WRITE (IOUT,1070)                                                 
      ICONV=2                                                           
      IFSS=0                                                            
      GO TO 100                                                         
!                                                                       
! - - - E N D   O F   I T E R A T I O N   L O O P                       
!                                                                       
  500 WRITE (IOUT,1100)                                                 
      WRITE (IOUT,1006) (EIGV(I),I=1,NROOT)                             
      WRITE (IOUT,1110)                                                 
      do J=1,NROOT                                                  
      WRITE (IOUT,1005) (R(K,J),K=1,NN)
      enddo                                 

        if (nprint == 1) then
        write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
            "95.00%"
        endif
!                                                                       
!     CALCULATE AND PRINT ERROR MEASURES                                
!                                                                       
      REWIND NSTIF                                                      
      READ (NSTIF) A                                                    
!                                                                       
      do L=1,NROOT                                                  
      RT=EIGV(L)                                                        
      CALL MULT(TT,A,R(1,L),MAXA,NN,NWK)                                
      VNORM=0.0D0                                                          
      do I=1,NN                                                     
      VNORM=VNORM + TT(I)*TT(I)
      enddo                                         
      CALL MULT(W,B,R(1,L),MAXA,NN,NWM)                                 
      WNORM=0.0D0                                                          
      do I=1,NN                                                     
      TT(I)=TT(I) - RT*W(I)                                             
      WNORM=WNORM + TT(I)*TT(I)
      enddo                                         
      VNORM=SQRT(VNORM)                                                 
      WNORM=SQRT(WNORM)                                                 
      D(L)=WNORM/VNORM                                                  
!  580 CONTINUE 
      enddo                                                         
      WRITE (IOUT,1115)                                                 
      WRITE (IOUT,1005) (D(I),I=1,NROOT)                                
!                                                                       
!     APPLY STURM SEQUENCE CHECK                                        
!                                                                       
      IF (IFSS.EQ.0) GO TO 900                                          
      CALL SCHECK (EIGV,RTOLV,BUP,BLO,BUPC,NC,NEI,RTOL,SHIFT,IOUT)    
!                                                                       
      WRITE (IOUT,1120) SHIFT                                           
!                                                                       
!     SHIFT MATRIX A                                                    
!                                                                       
      REWIND NSTIF                                                      
      READ (NSTIF) A                                                    
      IF (NWM.GT.NN) GO TO 645                                          
      do I=1,NN                                                     
      II=MAXA(I)                                                        
      A(II)=A(II) - B(I)*SHIFT
      enddo                                          
      GO TO 660                                                         
  645 do I=1,NWK                                                    
      A(I)=A(I) - B(I)*SHIFT
      enddo                                            
!                                                                       
!     FACTORIZE SHIFTED MATRIX                                          
!                                                                       
  660 ISH=1                                                             
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)                                  
!                                                                       
!     COUNT NUMBER OF NEGATIVE DIAGONAL ELEMENTS                        
!                                                                       
      NSCH=0                                                            
      do I=1,NN                                                     
      II=MAXA(I)                                                        
      IF (A(II).LT.0.) NSCH=NSCH + 1                                    
!  664 CONTINUE 
      enddo                                                         
      IF (NSCH.EQ.NEI) GO TO 670                                        
      NMIS=NSCH - NEI                                                   
      WRITE (IOUT,1130) NMIS                                            
      GO TO 900                                                         
  670 WRITE (IOUT,1140) NSCH                                            
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 1002 FORMAT (' ',10F10.0)                                              
 1005 FORMAT (' ',12E11.4)                                              
 1006 FORMAT (' ',6E22.14)                                              
 1007 FORMAT (///,' STOP, NC IS LARGER THAN THE NUMBER OF MASS DEGREES OF FREEDOM')                                     
 1008 FORMAT (///,' DEGREES OF FREEDOM EXCITED BY UNIT STARTING ITERATION VECTORS')                                      
 1010 FORMAT (//,' I T E R A T I O N   N U M B E R ',I8)                
 1020 FORMAT (/,' PROJECTION OF A (MATRIX AR)')                         
 1030 FORMAT (/,' PROJECTION OF B (MATRIX BR)')                         
 1035 FORMAT (/,' EIGENVALUES OF AR-LAMBDA*BR')                         
 1040 FORMAT (//,' AR AND BR AFTER JACOBI DIAGONALIZATION')             
 1050 FORMAT (/,' ERROR BOUNDS REACHED ON EIGENVALUES')                 
 1060 FORMAT (///,' CONVERGENCE REACHED FOR RTOL ',E10.4)               
 1070 FORMAT (' *** NO CONVERGENCE IN MAXIMUM NUMBER OF ITERATIONS PERMITTED')             
 1100 FORMAT (///,' THE CALCULATED EIGENVALUES ARE')                    
 1115 FORMAT (//,' ERROR MEASURES ON THE EIGENVALUES')                  
 1110 FORMAT (//,' THE CALCULATED EIGENVECTORS ARE',/)                  
 1120 FORMAT (///,' CHECK APPLIED AT SHIFT ',E22.14)                    
 1130 FORMAT (//,' THERE ARE ',I8,' EIGENVALUES MISSING')               
 1140 FORMAT (//,' WE FOUND THE LOWEST ',I8,' EIGENVALUES')             
!                                                                       
END SUBROUTINE SPACE
