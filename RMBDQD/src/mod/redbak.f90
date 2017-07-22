SUBROUTINE REDBAK (A,V,MAXA,NN)                                   
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS            . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION A(1),V(1),MAXA(1)                                       
!                                                                       
      do N=1,NN                                                     
      KL=MAXA(N) + 1                                                    
      KU=MAXA(N+1) - 1                                                  
      !IF (KU-KL) 400,410,410
      if ((KU-KL).LT.0) then
        go to 400
      else
        go to 410
      endif                                            
  410 K=N                                                               
      C=0.                                                              
      do KK=KL,KU                                                   
      K=K - 1                                                           
      C=C + A(KK)*V(K)
      enddo                                                  
      V(N)=V(N) - C                                                     
  400 CONTINUE
      enddo                                                          
!                                                                       
      do N=1,NN                                                     
      K=MAXA(N)                                                         
      V(N)=V(N)/A(K)
      enddo                                                    
      IF (NN.EQ.1) GO TO 900                                            
      N=NN                                                              
      do L=2,NN                                                     
      KL=MAXA(N) + 1                                                    
      KU=MAXA(N+1) - 1                                                  
      !IF (KU-KL) 500,510,510
      if ((KU-KL).LT.0) then
        go to 500
      else
        go to 510
      endif                                            
  510 K=N                                                               
      do KK=KL,KU                                                   
      K=K - 1                                                           
      V(K)=V(K) - A(KK)*V(N)
      enddo                                            
  500 N=N - 1
      enddo                                                           
!                                                                       
  900 RETURN                                                            
END SUBROUTINE REDBAK
