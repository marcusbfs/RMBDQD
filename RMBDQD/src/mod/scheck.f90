SUBROUTINE SCHECK (EIGV,RTOLV,BUP,BLO,BUPC,NC,NEI,RTOL,SHIFT,IOUT)                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO EVALUATE SHIFT FOR STURM SEQUENCE CHECK                 . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION EIGV(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),NEIV(NC)    
!                                                                       
      FTOL=0.001D0                                                         
!                                                                       
      do I=1,NC                                                     
      BUP(I)=EIGV(I)*(1.0D0+FTOL)                                          
      BLO(I)=EIGV(I)*(1.0D0-FTOL)
      enddo                                          
      NROOT=0                                                           
      do I=1,NC                                                     
      IF (RTOLV(I).LT.RTOL) NROOT=NROOT + 1 
      enddo                            
      IF (NROOT.GE.1) GO TO 200                                         
      WRITE (IOUT,1010)                                                 
      GO TO 800                                                         
!                                                                       
!      FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS                         
!                                                                       
  200 do I=1,NROOT                                                  
      NEIV(I)=1 
      enddo                                                        
      IF (NROOT.NE.1) GO TO 260                                         
      BUPC(1)=BUP(1)                                                    
      LM=1                                                              
      L=1                                                               
      I=2                                                               
      GO TO 295                                                         
  260 L=1                                                               
      I=2                                                               
  270 IF (BUP(I-1).LE.BLO(I)) GO TO 280                                 
      NEIV(L)=NEIV(L) + 1                                               
      I=I + 1                                                           
      IF (I.LE.NROOT) GO TO 270                                         
  280 BUPC(L)=BUP(I-1)                                                  
      IF (I.GT.NROOT) GO TO 290                                         
      L=L + 1                                                           
      I=I + 1                                                           
      IF (I.LE.NROOT) GO TO 270                                         
      BUPC(L)=BUP(I-1)                                                  
  290 LM=L                                                              
      IF (NROOT.EQ.NC) GO TO 300                                        
  295 IF (BUP(I-1).LE.BLO(I)) GO TO 300                                 
      IF (RTOLV(I).GT.RTOL) GO TO 300                                   
      BUPC(L)=BUP(I)                                                    
      NEIV(L)=NEIV(L) + 1                                               
      NROOT=NROOT + 1                                                   
      IF (NROOT.EQ.NC) GO TO 300                                        
      I=I + 1                                                           
      GO TO 295                                                         
!                                                                       
!      FIND SHIFT                                                       
!                                                                       
  300 WRITE (IOUT,1020)                                                 
      WRITE (IOUT,1005) (BUPC(I),I=1,LM)                                
      WRITE (IOUT,1030)                                                 
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                
      LL=LM - 1                                                         
      IF (LM.EQ.1) GO TO 310                                            
  330 do I=1,LL                                                     
      NEIV(L)=NEIV(L) + NEIV(I)
      enddo                                         
      L=L - 1                                                           
      LL=LL - 1                                                         
      IF (L.NE.1) GO TO 330                                             
  310 WRITE (IOUT,1040)                                                 
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                
      L=0                                                               
      do I=1,LM                                                     
      L=L + 1                                                           
      IF (NEIV(I).GE.NROOT) GO TO 350                                   
!  340 CONTINUE
      enddo                                                          
  350 SHIFT=BUPC(L)                                                     
      NEI=NEIV(L)                                                       
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 1005 FORMAT (' ',6E22.14)                                              
 1006 FORMAT (' ',6I22)                                                 
 1010 FORMAT (' *** ERROR ***  SOLUTION STOP IN *SCHECK*',/,' NO EIGENVALUES FOUND',/)                                
 1020 FORMAT (///,' UPPER BOUNDS ON EIGENVALUE CLUSTERS')               
 1030 FORMAT (//,' NO. OF EIGENVALUES IN EACH CLUSTER')                 
 1040 FORMAT (' NO. OF EIGENVALUES LESS THAN UPPER BOUNDS')             
END SUBROUTINE SCHECK 
