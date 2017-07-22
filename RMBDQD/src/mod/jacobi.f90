SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)       
! ..................................................................... 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE            . 
! .        GENERALIZED JACOBI ITERATION                               . 
! ..................................................................... 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)                       
!                                                                       
!     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES                    
!                                                                       
      N1=N + 1                                                          
      II=1                                                              
      do I=1,N                                                       
      IF (A(II).GT.0. .AND. B(II).GT.0.) GO TO 4                        
      WRITE (IOUT,2020) II,A(II),B(II)                                  
      GO TO 800                                                         
    4 D(I)=A(II)/B(II)                                                  
      EIGV(I)=D(I)                                                      
      II=II + N1 - I
      enddo                                                    
      do I=1,N                                                       
      do J=1,N                                                       
      X(I,J)=0.
      enddo                                                         
      X(I,I)=1.
      enddo                                                         
      IF (N.EQ.1) GO TO 900                                             
!                                                                       
!     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION                      
!                                                                       
      NSWEEP=0                                                          
      NR=N - 1                                                          
   40 NSWEEP=NSWEEP + 1                                                 
      IF (IFPR.EQ.1) WRITE (IOUT,2000) NSWEEP                           
!                                                                       
!     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE  
!     ZEROING                                                           
!                                                                       
      EPS=(0.01D0)**(NSWEEP*2)                                             
      do J=1,NR                                                     
      JP1=J + 1                                                         
      JM1=J - 1                                                         
      LJK=JM1*N - JM1*J/2                                               
      JJ=LJK + J                                                        
      do K=JP1,N                                                    
      KP1=K + 1                                                         
      KM1=K - 1                                                         
      JK=LJK + K                                                        
      KK=KM1*N - KM1*K/2 + K                                            
      EPTOLA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                
      EPTOLB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                
      IF (EPTOLA.LT.EPS .AND. EPTOLB.LT.EPS) GO TO 210                  
!                                                                       
!     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA 
!     AND CG                                                            
!                                                                       
      AKK=A(KK)*B(JK) - B(KK)*A(JK)                                     
      AJJ=A(JJ)*B(JK) - B(JJ)*A(JK)                                     
      AB=A(JJ)*B(KK) - A(KK)*B(JJ)                                      
      SCALE=A(KK)*B(KK)                                                 
      ABCH=AB/SCALE                                                     
      AKKCH=AKK/SCALE                                                   
      AJJCH=AJJ/SCALE                                                   
      CHECK=(ABCH*ABCH+4.0D0*AKKCH*AJJCH)/4.0D0                             
      !IF (CHECK) 50,60,60 
      if (CHECK.LT.0) then
        go to 50
      else
        go to 60
      endif                                              
   50 WRITE (IOUT,2020) JJ,A(JJ),B(JJ)                                  
      GO TO 800                                                         
   60 SQCH=SCALE*DSQRT(CHECK)                                            
      D1=AB/2.0D0 + SQCH                                                   
      D2=AB/2.0D0 - SQCH                                                   
      DEN=D1                                                            
      IF (DABS(D2).GT.DABS(D1)) DEN=D2                                    
      !IF (DEN) 80,70,80
      if (DEN.LT.0) then
        go to 80
      elseif (DEN.GT.0) then
        go to 80
      else
        go to 70
      endif                                                 
   70 CA=0.0D0                                                             
      CG=-A(JK)/A(KK)                                                   
      GO TO 90                                                          
   80 CA=AKK/DEN                                                        
      CG=-AJJ/DEN                                                       
!                                                                       
!     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL 
!     ELEMENT                                                           
!                                                                       
!   90 IF (N-2) 100,190,100 
   90 if ((N-2).LT.0) then
        go to 100
      elseif ((N-2).GT.0) then
        go to 100
      else
        go to 190
      endif                                             
!  100 IF (JM1-1) 130,110,110
  100 if ((JM1-1).LT.0) then
        go to 130
      else
        go to 110
      endif                                            
  110 do I=1,JM1                                                    
      IM1=I - 1                                                         
      IJ=IM1*N - IM1*I/2 + J                                            
      IK=IM1*N - IM1*I/2 + K                                            
      AJ=A(IJ)                                                          
      BJ=B(IJ)                                                          
      AK=A(IK)                                                          
      BK=B(IK)                                                          
      A(IJ)=AJ + CG*AK                                                  
      B(IJ)=BJ + CG*BK                                                  
      A(IK)=AK + CA*AJ                                                  
      B(IK)=BK + CA*BJ
      enddo                                                  
!  130 IF (KP1-N) 140,140,160 
  130 if ((KP1-N).LE.0) then
        go to 140
      else
        go to 160
      endif                                           
  140 LJI=JM1*N - JM1*J/2                                               
      LKI=KM1*N - KM1*K/2                                               
      do I=KP1,N                                                    
      JI=LJI + I                                                        
      KI=LKI + I                                                        
      AJ=A(JI)                                                          
      BJ=B(JI)                                                          
      AK=A(KI)                                                          
      BK=B(KI)                                                          
      A(JI)=AJ + CG*AK                                                  
      B(JI)=BJ + CG*BK                                                  
      A(KI)=AK + CA*AJ                                                  
      B(KI)=BK + CA*BJ
      enddo                                                  
!  160 IF (JP1-KM1) 170,170,190
  160 if ((JP1-KM1).LE.0) then
        go to 170
      else
        go to 190
      endif                                          
  170 LJI=JM1*N - JM1*J/2                                               
      do I=JP1,KM1                                                  
      JI=LJI + I                                                        
      IM1=I - 1                                                         
      IK=IM1*N - IM1*I/2 + K                                            
      AJ=A(JI)                                                          
      BJ=B(JI)                                                          
      AK=A(IK)                                                          
      BK=B(IK)                                                          
      A(JI)=AJ + CG*AK                                                  
      B(JI)=BJ + CG*BK                                                  
      A(IK)=AK + CA*AJ                                                  
      B(IK)=BK + CA*BJ
      enddo                                                  
  190 AK=A(KK)                                                          
      BK=B(KK)                                                          
      A(KK)=AK + 2.0D0*CA*A(JK) + CA*CA*A(JJ)                              
      B(KK)=BK + 2.0D0*CA*B(JK) + CA*CA*B(JJ)                              
      A(JJ)=A(JJ) + 2.0D0*CG*A(JK) + CG*CG*AK                              
      B(JJ)=B(JJ) + 2.0D0*CG*B(JK) + CG*CG*BK                              
      A(JK)=0.0D0                                                          
      B(JK)=0.0D0                                                          
!                                                                       
!     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION                 
!                                                                       
      do I=1,N                                                      
      XJ=X(I,J)                                                         
      XK=X(I,K)                                                         
      X(I,J)=XJ + CG*XK                                                 
      X(I,K)=XK + CA*XJ
      enddo                                                 
  210 CONTINUE
      enddo
      enddo                                                          
!                                                                       
!     UPDATE THE EIGENVALUES AFTER EACH SWEEP                           
!                                                                       
      II=1                                                              
      do I=1,N                                                      
      IF (A(II).GT.0.0D0 .AND. B(II).GT.0.0D0) GO TO 215                      
      WRITE (IOUT,2020) II,A(II),B(II)                                  
      GO TO 800                                                         
  215 EIGV(I)=A(II)/B(II)                                               
      II=II + N1 - I
      enddo                                                    
      IF (IFPR.EQ.0) GO TO 230                                          
      WRITE (IOUT,2030)                                                 
      WRITE (IOUT,2010) (EIGV(I),I=1,N)                                 
!                                                                       
!     CHECK FOR CONVERGENCE                                             
!                                                                       
  230 do I=1,N                                                      
      TOL=RTOL*D(I)                                                     
      DIF=DABS(EIGV(I)-D(I))                                             
      IF (DIF.GT.TOL) GO TO 280                                         
 ! 240 CONTINUE
      enddo                                                          
!                                                                       
!     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS        
!     REQUIRED                                                          
!                                                                       
      EPS=RTOL**2.0D0                                                       
      do J=1,NR                                                     
      JM1=J - 1                                                         
      JP1=J + 1                                                         
      LJK=JM1*N - JM1*J/2                                               
      JJ=LJK + J                                                        
      do K=JP1,N                                                    
      KM1=K - 1                                                         
      JK=LJK + K                                                        
      KK=KM1*N - KM1*K/2 + K                                            
      EPSA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                  
      EPSB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                  
      IF (EPSA.LT.EPS .AND. EPSB.LT.EPS) GO TO 250                      
      GO TO 280                                                         
  250 CONTINUE 
      enddo
      enddo                                                         
!                                                                       
!     SCALE EIGENVECTORS                                                
!                                                                       
  255 II=1                                                              
      do I=1,N                                                      
      BB=DSQRT(B(II))                                                    
      do K=1,N                                                      
      X(K,I)=X(K,I)/BB
      enddo                                                  
      II=II + N1 - I
      enddo                                                    
      GO TO 900                                                         
!                                                                       
!     UPDATE  D  MATRIX AND START NEW SWEEP, IF ALLOWED                 
!                                                                       
  280 do I=1,N                                                      
      D(I)=EIGV(I)
      enddo                                                      
      IF (NSWEEP.LT.NSMAX) GO TO 40                                     
      GO TO 255                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 2000 FORMAT (//,' SWEEP NUMBER IN *JACOBI* = ',I8)                     
 2010 FORMAT (' ',6E20.12)                                              
 2020 FORMAT (' *** ERROR *** SOLUTION STOP',' II = ',I8,' A(II) = ',E20.12,' B(II) = ',E20.12)        
 2030 FORMAT (/,' CURRENT EIGENVALUES IN *JACOBI* ARE',/)               
END SUBROUTINE JACOBI
