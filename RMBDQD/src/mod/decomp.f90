SUBROUTINE DECOMP (A,MAXA,NN,ISH,IOUT)                            
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF               . 
! .        STIFFNESS MATRIX                                           . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION A(1),MAXA(1)                                            
      IF (NN.EQ.1) GO TO 900                                            
!                                                                       
      do N=1,NN                                                     
      KN=MAXA(N)                                                        
      KL=KN + 1                                                         
      KU=MAXA(N+1) - 1                                                  
      KH=KU - KL                                                        
      !IF (KH) 304,240,210 
      if (KH.LT.0) then
          go to 304
      elseif (KH.EQ.0) then
          go to 240
      else
      go to 210
      endif                                              
  210 K=N - KH                                                          
      IC=0                                                              
      KLT=KU                                                            
      do J=1,KH                                                     
      IC=IC + 1                                                         
      KLT=KLT - 1                                                       
      KI=MAXA(K)                                                        
      ND=MAXA(K+1) - KI - 1                                             
    !  IF (ND) 260,260,270
      if (ND.LE.0) then
        go to 260
      else
        go to 270
      endif                                               
  270 KK=MIN0(IC,ND)                                                    
      C=0.                                                              
      do L=1,KK                                                     
      C=C + A(KI+L)*A(KLT+L)
      enddo                                            
      A(KLT)=A(KLT) - C                                                 
  260 K=K + 1
      enddo                                                           
  240 K=N                                                               
      B=0.                                                              
      do KK=KL,KU                                                   
      K=K - 1                                                           
      KI=MAXA(K)                                                        
      C=A(KK)/A(KI)                                                     
      IF (DABS(C).LT.1.E07) GO TO 290                                    
      WRITE (IOUT,2010) N,C                                             
      GO TO 800                                                         
  290 B=B + C*A(KK)                                                     
      A(KK)=C 
      enddo                                                          
      A(KN)=A(KN) - B                                                   
!  304 IF (A(KN)) 310,310,200
  304 if (A(KN).LE.0) then
        go to 310
      else
        go to 200
      endif                                            
  310 IF (ISH.EQ.0) GO TO 320                                           
      IF (A(KN).EQ.0.) A(KN)=-1.E-16                                    
      GO TO 200                                                         
  320 WRITE (IOUT,2000) N,A(KN)                                         
      GO TO 800                                                         
  200 CONTINUE
      enddo                                                          
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 2000 FORMAT (//' STOP - STIFFNESS MATRIX',//,' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,' PIVOT = ',E20.12) 
 2010 FORMAT (//' STOP - STURM SEQUENCE CHECK FAILED',I8,//,' MULTIPLIER = ',E20.8)
 END SUBROUTINE DECOMP
