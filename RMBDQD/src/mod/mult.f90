SUBROUTINE MULT (TT,B,RR,MAXA,NN,NWM)                             
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT   . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION TT(1),B(1),RR(1),MAXA(1)                                
!                                                                       
      IF (NWM.GT.NN) GO TO 20                                           
      do I=1,NN                                                      
      TT(I)=B(I)*RR(I)
      enddo                                                  
      GO TO 900                                                         
!                                                                       
   20 do I=1,NN                                                      
      TT(I)=0.0D0
      enddo                                                          
      do I=1,NN                                                     
      KL=MAXA(I)                                                        
      KU=MAXA(I+1) - 1                                                  
      II=I + 1                                                          
      CC=RR(I)                                                          
      do KK=KL,KU                                                   
      II=II - 1                                                         
      TT(II)=TT(II) + B(KK)*CC
      enddo 
      enddo                                         
      IF (NN.EQ.1) GO TO 900                                            
      do I=2,NN                                                     
      KL=MAXA(I) + 1                                                    
      KU=MAXA(I+1) - 1                                                  
      !IF (KU-KL) 200,210,210
      if ((KU-KL).LT.0) then
        go to 200
      else
        go to 210
      endif                                            
  210 II=I                                                              
      AA=0.0D0                                                             
      do KK=KL,KU                                                   
      II=II - 1                                                         
      AA=AA + B(KK)*RR(II)
      enddo                                              
      TT(I)=TT(I) + AA                                                  
  200 CONTINUE
      enddo                                                          
!                                                                       
  900 RETURN                                                            
END SUBROUTINE MULT
