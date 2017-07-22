subroutine formnf(nn,nodof,nf)
   ! reform nf
      implicit none
      integer :: nodof,nn
      integer,intent(in out)::nf(nodof,nn)
      integer:: i,j,m
      m=0
      do j=1,ubound(nf,2)
         do i=1,ubound(nf,1)
            if(nf(i,j)/=0) then
               m=m+1; nf(i,j)=m
            end if
         end do
      end do
      return
end subroutine formnf
!
!
subroutine num_to_g(iel,nels,num,nn,nod,nodof,ndof,fnf,g)
      !finds the g vector from num and nf
      implicit none
      integer :: fnod,fnodof,fab(3),ii
      integer::i,k,nod,nodof,nn,ndof,nels,iel
      integer::num(nod),fnf(nodof,nn)  ; integer,intent(out):: g(ndof)
      fnod=ubound(num,1) ; fnodof=ubound(fnf,1)
      !write(*,*) "tentativa",fnf(:,iel)
      do i = 1 , nod
         k = i*nodof
         fab=fnf(:,num(i))
!
         !write(*,*) "tentativafab",fnf
!
         do ii=1,3
!
         g(k-nodof+ii) = fab(ii)
!
         end do
      end do
      !write(*,*) "tentativafabio",g
      return
end subroutine num_to_g
!
!
subroutine fkdiag(a,b,kdiag,g)
   ! finds the maximum bandwidth for each freedom
     !implicit none
     integer::a,b
     integer::g(a)
     integer::kdiag(b)
     integer::idof,i,iwp1,j,im,k
    
     idof=size(g)
      
     do i = 1,idof
       iwp1=1
       if(g(i)/=0) then
         do j=1,idof
           if(g(j)/=0) then
              im=g(i)-g(j)+1
              if(im>iwp1) iwp1=im
           end if
         end do
         k=g(i)
         if(iwp1>kdiag(k))kdiag(k)=iwp1
       end if
     end do
     write(11,*) "tentativafabio", kdiag
     return
end subroutine fkdiag
!
!
subroutine dpsample(element,ndim,nip,s,wt)
   ! returns the local coordinates of the integrating points
     !implicit none 
     double precision,intent(out)::s(nip,ndim),wt(nip)  
     character(*),intent(in):: element
     integer::nip
     double precision:: root3, r15 , w(3),v(9),b,c   
     root3 = 1.0/dsqrt(3.0d0)
     r15 = 0.2*dsqrt(15.0d0) 
     !nip = ubound( s , 1 ) 
     w = (/5./9.,8./9.,5./9./)
     v=(/5./9.*w,8./9.*w,5./9.*w/)
     select case (element)
          case('line')
            select case(nip)
                 case(1)
                   s(1,1)=0.  ;  wt(1)=2.   
                 case(2)
                   s(1,1)=root3  ; s(2,1)=-s(1,1)  ;  wt(1)=1.  ; wt(2)=1.
                 case(3)
                   s(1,1)=r15 ; s(2,1)=.0     ; s(3,1)=-s(1,1)
                   wt = w
                 case(4)
                   s(1,1)=.861136311594053  ; s(2,1)=.339981043584856  
                   s(3,1)=-s(2,1)  ; s(4,1)=-s(1,1)
                   wt(1)=.347854845137454 ; wt(2)=.652145154862546 
                   wt(3)=wt(2) ; wt(4)=wt(1)
                 case(5)
                   s(1,1)=.906179845938664 ; s(2,1)=.538469310105683  
                   s(3,1)=.0 ; s(4,1)=-s(2,1) ; s(5,1)=-s(1,1)
                   wt(1)=.236926885056189 ; wt(2)=.478628670499366
                   wt(3)=.568888888888889 ; wt(4)=wt(2) ; wt(5)=wt(1)
                 case(6)
                   s(1,1)=.932469514203152 ; s(2,1)=.661209386466265 
                   s(3,1)=.238619186083197
                   s(4,1)=-s(3,1) ; s(5,1)=-s(2,1) ; s(6,1)=-s(1,1)
                   wt(1)=.171324492379170 ; wt(2)=.360761573048139 
                   wt(3)=.467913934572691
                   wt(4)=wt(3); wt(5)=wt(2) ; wt(6)=wt(1)
                 case default
                     print*,"wrong number of integrating points for a line"
            end select
          case('triangle') 
            select case(nip)
                case(1)   ! for triangles weights multiplied by .5
                   s(1,1)=1./3.  ; s(1,2)=1./3.  ;  wt(1)= .5
                case(3)   
                   s(1,1)=.5 ;  s(1,2)=.5 ;  s(2,1)=.5  
                   s(2,2)=0.;  s(3,1)=0.  ;  s(3,2)=.5
                   wt(1)=1./3.  ;  wt(2)=wt(1) ; wt(3)=wt(1)   ; wt = .5*wt
                case(6)
                  s(1,1)=.816847572980459  ; s(1,2)=.091576213509771
                  s(2,1)=s(1,2);  s(2,2)=s(1,1) ;  s(3,1)=s(1,2); s(3,2)=s(1,2)
                  s(4,1)=.108103018168070 ;  s(4,2)=.445948490915965
                  s(5,1)=s(4,2) ;   s(5,2)=s(4,1) ;  s(6,1)=s(4,2)  ; s(6,2)=s(4,2)
                  wt(1)=.109951743655322 ;   wt(2)=wt(1)  ;   wt(3)=wt(1)
                  wt(4)=.223381589678011 ;   wt(5)=wt(4)  ;   wt(6)=wt(4)    ; wt = .5*wt
                case(7)
                  s(1,1)=1./3. ; s(1,2)=1./3.
                  s(2,1)=.797426985353087 ;s(2,2)=.101286507323456
                  s(3,1)=s(2,2) ;  s(3,2)=s(2,1) 
                  s(4,1)=s(2,2) ;  s(4,2)=s(2,2)
                  s(5,1)=.470142064105115 ;   s(5,2)=.059715871789770
                  s(6,1)=s(5,2) ; s(6,2)=s(5,1)
                  s(7,1)=s(5,1);  s(7,2)=s(5,1)
                  wt(1)=.225 ; wt(2)=.125939180544827 ;  wt(3)=wt(2);  wt(4)=wt(2)
                  wt(5)=.132394152788506;  wt(6)=wt(5)      ;  wt(7)=wt(5)     ;wt = .5*wt
              case(12)
                  s(1,1)=.873821971016996 ; s(1,2)=.063089014491502
                  s(2,1)=s(1,2) ;  s(2,2)=s(1,1);  s(3,1)=s(1,2) ;  s(3,2)=s(1,2)
                  s(4,1)=.501426509658179 ;  s(4,2)=.249286745170910
                  s(5,1)=s(4,2); s(5,2)=s(4,1)   ;  s(6,1)=s(4,2) ;  s(6,2)=s(4,2)
                  s(7,1)=.636502499121399 ;      s(7,2)=.310352451033785
                  s(8,1)=s(7,1) ;  s(8,2)=.053145049844816 ;  s(9,1)=s(7,2) ; s(9,2)=s(7,1)
                  s(10,1)=s(7,2) ; s(10,2)=s(8,2) ; s(11,1)=s(8,2);   s(11,2)=s(7,1)
                  s(12,1)=s(8,2) ;  s(12,2)=s(7,2)
                  wt(1)=.050844906370207 ; wt(2)=wt(1); wt(3)=wt(1)
                  wt(4)=.116786275726379 ; wt(5)=wt(4); wt(6)=wt(4)
                  wt(7)=.082851075618374 ; wt(8:12)=wt(7)           ; wt = .5*wt
              case(16)
                  s(1,1)=1./3. ;  s(1,2)=1./3.  ;  s(2,1)=.658861384496478
                  s(2,2)=.170569307751761 ; s(3,1)=s(2,2)   ;  s(3,2)=s(2,1)
                  s(4,1)=s(2,2)  ; s(4,2)=s(2,2)
                  s(5,1)=.898905543365938 ; s(5,2)=.050547228317031
                  s(6,1)=s(5,2);  s(6,2)=s(5,1) ; s(7,1)=s(5,2)  ;  s(7,2)=s(5,2)
                  s(8,1)=.081414823414554; s(8,2)=.459292588292723
                  s(9,1)=s(8,2)  ;  s(9,2)=s(8,1);  s(10,1)=s(8,2) ;  s(10,2)=s(8,2)
                  s(11,1)=.008394777409958; s(11,2)=.263112829634638
                  s(12,1)=s(11,1)    ;  s(12,2)=.728492392955404
                  s(13,1)=s(11,2) ;   s(13,2)=s(11,1)  ;  s(14,1)=s(11,2); s(14,2)=s(12,2)
                  s(15,1)=s(12,2) ;  s(15,2)=s(11,1) ;  s(16,1)=s(12,2) ;  s(16,2)=s(11,2)
                  wt(1)=.144315607677787 ; wt(2)=.103217370534718 ; wt(3)=wt(2); wt(4)=wt(2)
                  wt(5)=.032458497623198 ; wt(6)=wt(5)   ;  wt(7)=wt(5)
                  wt(8)=.095091634267284 ; wt(9)=wt(8)   ;  wt(10)=wt(8)
                  wt(11)=.027230314174435 ; wt(12:16) = wt(11)  ;     wt = .5*wt
              case default
                  print*,"wrong number of integrating points for a triangle"
            end select
          case ('quadrilateral')
            select case (nip)
              case(1)
                  s(1,1) = .0 ; wt(1) = 4.
              case(4)
                  s(1,1)=-root3; s(1,2)= root3
                  s(2,1)= root3; s(2,2)= root3
                  s(3,1)=-root3; s(3,2)=-root3
                  s(4,1)= root3; s(4,2)=-root3
                  wt = 1.0
              case(9)
                  s(1:7:3,1) = -r15; s(2:8:3,1) = .0
                  s(3:9:3,1) =  r15; s(1:3,2)   = r15
                  s(4:6,2)   =  .0 ; s(7:9,2)   =-r15
                  wt= v
              case default
                  print*,"wrong number of integrating points for a quadrilateral"
            end select
          case('tetrahedron')    
             s=0.0d0
             wt=0.0d0
             !write (*,*) "pontos de integracao 222", nip
              select case(nip)
              case(1)          ! for tetrahedra weights multiplied by 1/6
                  s(1,1)=.25    ; s(1,2)=.25  ;  s(1,3)=.25   ; wt(1)=1./6.
              case(4)
                  s(1,1)=.58541020 ; s(1,2)=.13819660  ;  s(1,3)=s(1,2)
                  s(2,2)=s(1,1) ; s(2,3)=s(1,2)  ;  s(2,1)=s(1,2)
                  s(3,3)=s(1,1) ; s(3,1)=s(1,2)  ;  s(3,2)=s(1,2)
                  s(4,1)=s(1,2) ; s(4,2)=s(1,2)  ;  s(4,3)=s(1,2) 
                  wt(1)=.25/6. ;  wt(2)=.25/6. ;  wt(3)=.25/6. ;  wt(4)=.25/6.
              case(5)
                  s(1,1)=.25  ;  s(1,2)=.25   ; s(1,3)=.25 ;  s(2,1)=.5
                  s(2,2)=1./6. ;  s(2,3)=s(2,2);  s(3,2)=.5
                  s(3,3)=1./6.  ;   s(3,1)=s(3,3)   ;   s(4,3)=.5
                  s(4,1)=1./6. ;    s(4,2)=s(4,1);    s(5,1)=1./6.
                  s(5,2)=s(5,1) ;  s(5,3)=s(5,1) 
                  wt(1)=-.8  ;  wt(2)=9./20. ;   wt(3:5)=wt(2)   ; wt =wt/6.  
              case(6)
                  wt = 4./3.        ;  s(6,3) = 1.
                  s(1,1)=-1. ;s(2,1)=1. ; s(3,2)=-1. ; s(4,2)=1. ;  s(5,3)=-1. 
              case default
                  print*,"wrong number of integrating points for a tetrahedron 11"
            end select
          case('hexahedron')
            select case ( nip )
              case(1)
                  s(1,1) = .0 ; wt(1) = 8.    
              case(8)   
                  s(1,1)= root3;s(1,2)= root3;s(1,3)= root3
                  s(2,1)= root3;s(2,2)= root3;s(2,3)=-root3
                  s(3,1)= root3;s(3,2)=-root3;s(3,3)= root3
                  s(4,1)= root3;s(4,2)=-root3;s(4,3)=-root3
                  s(5,1)=-root3;s(5,2)= root3;s(5,3)= root3
                  s(6,1)=-root3;s(6,2)=-root3;s(6,3)= root3
                  s(7,1)=-root3;s(7,2)= root3;s(7,3)=-root3
                  s(8,1)=-root3;s(8,2)=-root3;s(8,3)=-root3
                  wt = 1.0                                                
              case(14)
                  s=0.0;wt=0.0
                  b=0.795822425754222     ;      c=0.758786910639328
                  wt(1:6)=0.886426592797784   ; wt(7:) =  0.335180055401662
                  s(1,1)=-b ; s(2,1)=b  ;  s(3,2)=-b ;   s(4,2)=b
                  s(5,3)=-b   ;     s(6,3)=b
                  s(7:,:) = c
                  s(7,1)=-c  ;  s(7,2)=-c  ; s(7,3)=-c ; s(8,2)=-c ;   s(8,3)=-c
                  s(9,1)=-c  ;  s(9,3)=-c  ; s(10,3)=-c; s(11,1)=-c
                  s(11,2)=-c ;  s(12,2)=-c ; s(13,1)=-c
              case(15)
                  b=1.     ;      c=0.674199862
                  wt(1)=1.564444444 ;  wt(2:7)=0.355555556  ; wt(8:15)=0.537777778
                  s(2,1)=-b  ;    s(3,1)=b  ;    s(4,2)=-b  ;    s(5,2)=b
                  s(6,3)=-b  ;    s(7,3)=b  ;    s(8:,:)=c  ;    s(8,1)=-c
                  s(8,2)=-c  ;    s(8,3)=-c ;    s(9,2)=-c  ;    s(9,3)=-c
                  s(10,1)=-c ;    s(10,3)=-c  ;  s(11,3)=-c ;    s(12,1)=-c
                  s(12,2)=-c ;    s(13,2)=-c  ;  s(14,1)=-c                          
              case(27)
                  s=0.0;wt=0.0
                  wt = (/5./9.*v,8./9.*v,5./9.*v/)
                  s(1:7:3,1) = -r15; s(2:8:3,1) = .0
                  s(3:9:3,1) =  r15; s(1:3,3)   = r15
                  s(4:6,3)   =  .0 ; s(7:9,3)   =-r15
                  s(1:9,2)   = -r15
                  s(10:16:3,1) = -r15; s(11:17:3,1) = .0
                  s(12:18:3,1) =  r15; s(10:12,3)   = r15
                  s(13:15,3)   =  .0 ; s(16:18,3)   =-r15
                  s(10:18,2)   = .0   
                  s(19:25:3,1) = -r15; s(20:26:3,1) = .0
                  s(21:27:3,1) =  r15; s(19:21,3)   = r15
                  s(22:24,3)   =  .0 ; s(25:27,3)   =-r15
                  s(19:27,2)   =  r15
               case default
                  print*,"wrong number of integrating points for a hexahedron" 
             end select
           case default
             print*,"not a valid element type" 
     end select
     return
end subroutine dpsample
!
!
subroutine dpdeemat(dee,nst,e,v)
   ! returns the elastic dee matrix for given ih
   ! ih=3,plane strain; =4,axisymmetry or plane strain elastoplasticity
   ! =6 , three dimensional
     implicit none
     integer :: nst
     double precision,intent(in)::e,v
     double precision,intent(out)::dee(nst,nst)
   ! local variables
     double precision::v1,v2,c,vv; integer :: i,ih;  dee=0.0  ; ih = ubound(dee,1)
              v1 = 1. - v; c = e/((1.+v)*(1.-2.*v))
     select case (ih)
          case(3)
             dee(1,1)=v1*c; dee(2,2)=v1*c; dee(1,2)=v*c; dee(2,1)=v*c
             dee(3,3)=.5*c*(1.-2.*v)
          case(4)
             dee(1,1)=v1*c; dee(2,2)=v1*c; dee(4,4)=v1*c
             dee(3,3)=.5*c*(1.-2.*v) ; dee(1,2)=v*c; dee(2,1)=v*c
             dee(1,4)=v*c; dee(4,1)=v*c; dee(2,4)=v*c; dee(4,2)=v*c
          case(6)
             v2=v/(1.-v); vv=(1.-2.*v)/(1.-v)*.5
             do i=1,3; dee(i,i)=1.;end do; do i=4,6; dee(i,i)=vv; end do
             dee(1,2)=v2; dee(2,1)=v2; dee(1,3)=v2; dee(3,1)=v2
             dee(2,3)=v2; dee(3,2)=v2
             dee = dee*e/(2.*(1.+v)*vv)
          case default
             print*,'wrong size for dee matrix'
     end select
     return
end subroutine dpdeemat
!
!
subroutine dpshape_fun(fun,ndim,nip,nod,points,i)
     !implicit none  
     integer,intent(in):: i; double precision,intent(in)::points(nip,ndim)
     double precision :: eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3     !local variables
     double precision :: t1,t2,t3,t4,t5,t6,t7,t8,t9,x,y,z
     double precision :: zeta,xi0,eta0,zeta0; integer::xii(20),etai(20),zetai(20),l,ndim,nod
     double precision,intent(out)::fun(nod)
     ndim = ubound(points , 2); nod = ubound(fun , 1 )  
     !write(*,*) "aqui",ndim
     select case (ndim)
       case(1) ! one dimensional cases
         xi=points(i,1)
         select case(nod)
           case(2)
             t1=-1.-xi ; t2=1.-xi
             fun(1)=t2/2. ; fun(2)=-t1/2.
           case(3)
             t1=-1.-xi ; t2=-xi ; t3=1.-xi
             fun(1)=t2*t3/2. ; fun(2)=-t1*t3 ; fun(3)=t1*t2/2.
           case(4)
             t1=-1.-xi ; t2=-1./3.-xi ; t3=1./3.-xi ; t4=1.-xi
             fun(1)=t2*t3*t4*9./16.  ; fun(2)=-t1*t3*t4*27./16.
             fun(3)=t1*t2*t4*27./16. ; fun(4)=-t1*t2*t3*9./16.
           case(5)
             t1=-1.-xi ; t2=-0.5-xi ; t3=-xi ; t4=0.5-xi ; t5=1.-xi
             fun(1)=t2*t3*t4*t5*2./3. ; fun(2)=-t1*t3*t4*t5*8./3.
             fun(3)=t1*t2*t4*t5*4. ; fun(4)=-t1*t2*t3*t5*8./3.
             fun(5)=t1*t2*t3*t4*2./3.
           case default
             print*,"wrong number of nodes in shape_fun"
         end select
       case(2) ! two dimensional cases
         c1=points(i,1); c2=points(i,2); c3=1.-c1-c2 
         xi=points(i,1);  eta=points(i,2)
         etam=.25*(1.-eta); etap=.25*(1.+eta)
         xim=.25*(1.-xi); xip=.25*(1.+xi)
         select case(nod)
           case(3)
             fun = (/c1,c3,c2/)  
           case(6)
             fun(1)=(2.*c1-1.)*c1 ;  fun(6)=4.*c1*c2 ;  fun(5)=(2.*c2-1.)*c2
             fun(4)=4.*c2*c3      ;  fun(3)=(2.*c3-1.)*c3 ; fun(2)=4.*c3*c1
           case(15)
             t1=c1-.25  ;  t2=c1-.5 ;  t3=c1-.75   ;   t4=c2-.25
             t5=c2-.5   ;  t6=c2-.75 ;  t7=c3-.25  ;   t8=c3-.5 ;  t9=c3-.75
             fun(1)=32./3.*c1*t1*t2*t3   ;  fun(12)=128./3.*c1*c2*t1*t2
             fun(11)=64.*c1*c2*t1*t4     ;  fun(10)=128./3.*c1*c2*t4*t5
             fun(9)=32./3.*c2*t4*t5*t6   ;  fun(8)=128./3.*c2*c3*t4*t5
             fun(7)=64.*c2*c3*t4*t7      ;  fun(6)=128./3.*c2*c3*t7*t8
             fun(5)=32./3.*c3*t7*t8*t9   ;  fun(4)=128./3.*c3*c1*t7*t8
             fun(3)=64.*c3*c1*t1*t7      ;  fun(2)=128./3.*c3*c1*t1*t2
             fun(13)=128.*c1*c2*t1*c3    ;  fun(15)=128.*c1*c2*c3*t4
             fun(14)=128.*c1*c2*c3*t7      
           case(4)
             fun=(/4.*xim*etam,4.*xim*etap,4.*xip*etap,4.*xip*etam/)
           case(8)
             fun=(/4.*etam*xim*(-xi-eta-1.),32.*etam*xim*etap,&
                  4.*etap*xim*(-xi+eta-1.),32.*xim*xip*etap, &
                  4.*etap*xip*(xi+eta-1.), 32.*etap*xip*etam,&
                  4.*xip*etam*(xi-eta-1.), 32.*xim*xip*etam/)
           case(9)
             etam = eta - 1.; etap= eta + 1.; xim = xi - 1.; xip = xi + 1.
             fun=(/.25*xi*xim*eta*etam,-.5*xi*xim*etap*etam,&
                   .25*xi*xim*eta*etap,-.5*xip*xim*eta*etap,&
                   .25*xi*xip*eta*etap,-.5*xi*xip*etap*etam,&
                   .25*xi*xip*eta*etam,-.5*xip*xim*eta*etam,xip*xim*etap*etam/)
           case default
             print*,"wrong number of nodes in shape_fun"
         end select
       case(3) ! three dimensional cases == tetraedro
         xi=points(i,1); eta=points(i,2); zeta=points(i,3)
         etam=1.-eta ;  xim=1.-xi  ;  zetam=1.-zeta
         etap=eta+1. ;  xip=xi+1.   ;  zetap=zeta+1.
         select case(nod)
           case(4)
             fun(1)=xi   ;   fun(2)= eta ;  fun(3)=zeta 
             fun(4)=1.-fun(1)-fun(2)-fun(3)
           case(8)
             fun=(/.125*xim*etam*zetam,.125*xim*etam*zetap,.125*xip*etam*zetap,&
                   .125*xip*etam*zetam,.125*xim*etap*zetam,.125*xim*etap*zetap,&
                   .125*xip*etap*zetap,.125*xip*etap*zetam/)
           case(14) !type 6 element
             x = points(i,1);  y = points(i,2);  z = points(i,3)
             fun(1)=((x*y+x*z+2.*x+y*z+2.*y+2.*z+2.)*(x-1.)*(y-1.)*(z-1.))/8.
             fun(2)=((x*y-x*z-2.*x+y*z+2.*y-2.*z-2.)*(x-1.)*(y+1.)*(z-1.))/8.
             fun(3)=((x*y+x*z+2.*x-y*z-2.*y-2.*z-2.)*(x+1.)*(y-1.)*(z-1.))/8.
             fun(4)=((x*y-x*z-2.*x-y*z-2.*y+2.*z+2.)*(x+1.)*(y+1.)*(z-1.))/8.
             fun(5)=-((x*y-x*z+2.*x-y*z+2.*y-2.*z+2.)*(x-1.)*(y-1.)*(z+1.))/8.
             fun(6)=-((x*y+x*z-2.*x-y*z+2.*y+2.*z-2.)*(x-1.)*(y+1.)*(z+1.))/8.
             fun(7)=-((x*y-x*z+2.*x+y*z-2.*y+2.*z-2.)*(x+1.)*(y-1.)*(z+1.))/8.
             fun(8)=-((x*y+x*z-2.*x+y*z-2.*y-2.*z+2.)*(x+1.)*(y+1.)*(z+1.))/8.
             fun(9)=-((x+1.)*(x-1.)*(y+1.)*(y-1.)*(z-1.))/2.
             fun(10)=((x+1.)*(x-1.)*(y+1.)*(y-1.)*(z+1.))/2.
             fun(11)=-((x+1.)*(x-1.)*(y-1.)*(z+1.)*(z-1.))/2.
             fun(12)=((x+1.)*(x-1.)*(y+1.)*(z+1.)*(z-1.))/2.
             fun(13)=-((x-1.)*(y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
             fun(14)=((x+1.)*(y+1.)*(y-1.)*(z+1.)*(z-1.))/2.      
           case(20)
             xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
             etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
             zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
             do l=1,20 
                xi0=xi*xii(l); eta0=eta*etai(l); zeta0=zeta*zetai(l)
                if(l==4.or.l==8.or.l==16.or.l==20) then
                   fun(l)=.25*(1.-xi*xi)*(1.+eta0)*(1.+zeta0)
                else if(l>=9.and.l<=12)then
                   fun(l)=.25*(1.+xi0)*(1.-eta*eta)*(1.+zeta0)
                else if(l==2.or.l==6.or.l==14.or.l==18) then
                   fun(l)=.25*(1.+xi0)*(1.+eta0)*(1.-zeta*zeta)
                else
                   fun(l)=.125*(1.+xi0)*(1.+eta0)*(1.+zeta0)*(xi0+eta0+zeta0-2)
                end if
             end do
           case default
             print*,"wrong number of nodes in shape_fun"
         end select
       case default
         print*,"wrong number of dimensions in shape_fun"  
     end select
     return
end subroutine dpshape_fun
!
!
subroutine dpshape_der(der,ndim,nod,nip,points,i)
     implicit none
     integer :: xii(20), etai(20), zetai(20) ,l,ndim , nod, nip   ! local variables
     integer,intent(in):: i; double precision,intent(in)::points(nip,ndim)
     double precision,intent(out)::der(ndim,nod)
     double precision::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3 ! local variables
     double precision:: t1,t2,t3,t4,t5,t6,t7,t8,t9 ,x2p1,x2m1,e2p1,e2m1,zetam,zetap,x,y,z
     integer:: iii(20),lll
     ndim = ubound(der , 1); nod = ubound(der , 2)
     select case (ndim)
       case(1) ! one dimensional case
         xi=points(i,1)
         select case (nod)
           case(2)
             der(1,1)=-0.5 ; der(1,2)=0.5
           case(3)
             t1=-1.-xi ; t2=-xi  ; t3=1.-xi
             der(1,1)=-(t3+t2)/2.  ; der(1,2)=(t3+t1)    
             der(1,3)=-(t2+t1)/2.   
           case(4)
             t1=-1.-xi ; t2=-1./3.-xi ; t3=1./3.-xi ; t4=1.-xi
             der(1,1)=-(t3*t4+t2*t4+t2*t3)*9./16.     
             der(1,2)=(t3*t4+t1*t4+t1*t3)*27./16. 
             der(1,3)=-(t2*t4+t1*t4+t1*t2)*27./16. 
             der(1,4)=(t2*t3+t1*t3+t1*t2)*9./16.   
           case(5)
             t1=-1.-xi ; t2=-0.5-xi ; t3=-xi ; t4=0.5-xi ; t5=1.-xi
             der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*2./3.   
             der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*8./3.
             der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*4. 
             der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*8./3.
             der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*2./3.
           case default
             print*,"wrong number of nodes in shape_der"        
         end select
       case(2)      ! two dimensional elements
         xi=points(i,1); eta=points(i,2) ; c1=xi ; c2=eta ; c3=1.-c1-c2
         etam=.25*(1.-eta); etap=.25*(1.+eta); xim=.25*(1.-xi); xip=.25*(1.+xi)
         x2p1=2.*xi+1. ;   x2m1=2.*xi-1. ;  e2p1=2.*eta+1. ;   e2m1=2.*eta-1.    
         select case (nod)
           case(3)
             der(1,1)=1.;der(1,3)=0.;der(1,2)=-1.
             der(2,1)=0.;der(2,3)=1.;der(2,2)=-1.
           case(6) 
             der(1,1)=4.*c1-1. ;  der(1,6)=4.*c2;  der(1,5)=0.  ; der(1,4)=-4.*c2
             der(1,3)=-(4.*c3-1.);  der(1,2)=4.*(c3-c1);   der(2,1)=0.
             der(2,6)=4.*c1 ; der(2,5)=4.*c2-1.; der(2,4)=4.*(c3-c2)
             der(2,3)=-(4.*c3-1.)  ; der(2,2)=-4.*c1
           case(15)                          
             t1=c1-.25  ;  t2=c1-.5 ;  t3=c1-.75   ;   t4=c2-.25
             t5=c2-.5   ;  t6=c2-.75 ;  t7=c3-.25  ;   t8=c3-.5 ;  t9=c3-.75
             der(1,1)=32./3.*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
             der(1,12)=128./3.*c2*(t2*(t1+c1)+c1*t1) ;  der(1,11)=64.*c2*t4*(t1+c1)
             der(1,10)=128./3.*c2*t4*t5  ; der(1,9)=0. ; der(1,8)=-128./3.*c2*t4*t5
             der(1,7)=-64.*c2*t4*(t7+c3) ; der(1,6)=-128./3.*c2*(t8*(t7+c3)+c3*t7)
             der(1,5)=-32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
             der(1,4)=128./3.*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
             der(1,3)=64.*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
             der(1,2)=128./3.*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
             der(1,13)=128.*c2*(c3*(t1+c1)-c1*t1) ;  der(1,15)=128.*c2*t4*(c3-c1)
             der(1,14)=128.*c2*(c3*t7-c1*(t7+c3))
             der(2,1)=0.0 ;  der(2,12)=128./3.*c1*t1*t2;  der(2,11)=64.*c1*t1*(t4+c2)
             der(2,10)=128./3.*c1*(t5*(t4+c2)+c2*t4)
             der(2,9)=32./3.*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
             der(2,8)=128./3.*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
             der(2,7)=64.*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
             der(2,6)=128./3.*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
             der(2,5)=-32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
             der(2,4)=-128./3.*c1*(t8*(t7+c3)+c3*t7)
             der(2,3)=-64.*c1*t1*(t7+c3)  ;  der(2,2)=-128./3.*c1*t1*t2
             der(2,13)=128.*c1*t1*(c3-c2)
             der(2,15)=128.*c1*(c3*(t4+c2)-c2*t4)
             der(2,14)=128.*c1*(c3*t7-c2*(c3+t7))        
           case (4)                                                              
             der(1,1)=-etam; der(1,2)=-etap; der(1,3)=etap; der(1,4)=etam
             der(2,1)=-xim; der(2,2)=xim; der(2,3)=xip; der(2,4)=-xip
           case(8)
             der(1,1)=etam*(2.*xi+eta); der(1,2)=-8.*etam*etap
             der(1,3)=etap*(2.*xi-eta); der(1,4)=-4.*etap*xi
             der(1,5)=etap*(2.*xi+eta); der(1,6)=8.*etap*etam
             der(1,7)=etam*(2.*xi-eta); der(1,8)=-4.*etam*xi
             der(2,1)=xim*(xi+2.*eta); der(2,2)=-4.*xim*eta
             der(2,3)=xim*(2.*eta-xi); der(2,4)=8.*xim*xip
             der(2,5)=xip*(xi+2.*eta); der(2,6)=-4.*xip*eta
             der(2,7)=xip*(2.*eta-xi); der(2,8)=-8.*xim*xip   
           case(9)
             etam = eta - 1.; etap = eta + 1.; xim = xi - 1.; xip = xi + 1.
             der(1,1)=.25*x2m1*eta*etam  ;   der(1,2)=-.5*x2m1*etap*etam
             der(1,3)=.25*x2m1*eta*etap  ;   der(1,4)=-xi*eta*etap
             der(1,5)=.25*x2p1*eta*etap  ;   der(1,6)=-.5*x2p1*etap*etam
             der(1,7)=.25*x2p1*eta*etam  ;   der(1,8)=-xi*eta*etam
             der(1,9)=2.*xi*etap*etam    ;   der(2,1)=.25*xi*xim*e2m1
             der(2,2)=-xi*xim*eta        ;   der(2,3)=.25*xi*xim*e2p1
             der(2,4)=-.5*xip*xim*e2p1   ;   der(2,5)=.25*xi*xip*e2p1
             der(2,6)=-xi*xip*eta        ;   der(2,7)=.25*xi*xip*e2m1
             der(2,8)=-.5*xip*xim*e2m1   ;   der(2,9)=2.*xip*xim*eta
           case default
             print*,"wrong number of nodes in shape_der"        
         end select
       case(3)  ! three dimensional elements
         xi=points(i,1); eta=points(i,2); zeta=points(i,3)
         etam=1.-eta ; xim=1.-xi;  zetam=1.-zeta
         etap=eta+1. ; xip=xi+1. ;  zetap=zeta+1.
         select case (nod)
           case(4)
             der(1:3,1:4) = .0
             der(1,1)=1.;  der(2,2)=1.  ;  der(3,3)=1.
             der(1,4)=-1. ;  der(2,4)=-1. ;  der(3,4)=-1. 
           case(8)
             der(1,1)=-.125*etam*zetam    ;   der(1,2)=-.125*etam*zetap
             der(1,3)=.125*etam*zetap     ;   der(1,4)=.125*etam*zetam
             der(1,5)=-.125*etap*zetam    ;   der(1,6)=-.125*etap*zetap
             der(1,7)=.125*etap*zetap     ;   der(1,8)=.125*etap*zetam
             der(2,1)=-.125*xim*zetam     ;   der(2,2)=-.125*xim*zetap
             der(2,3)=-.125*xip*zetap     ;   der(2,4)=-.125*xip*zetam
             der(2,5)=.125*xim*zetam      ;   der(2,6)=.125*xim*zetap
             der(2,7)=.125*xip*zetap      ;   der(2,8)=.125*xip*zetam
             der(3,1)=-.125*xim*etam      ;   der(3,2)=.125*xim*etam
             der(3,3)=.125*xip*etam       ;   der(3,4)=-.125*xip*etam
             der(3,5)=-.125*xim*etap      ;   der(3,6)=.125*xim*etap
             der(3,7)=.125*xip*etap       ;   der(3,8)=-.125*xip*etap  
           case(14) ! type 6 element
             x= points(i,1)    ;   y= points(i,2)  ;    z= points(i,3) 
             der(1,1)=((2.*x*y+2.*x*z+4.*x+y*z+y+z)*(y-1.)*(z-1.))/8.
             der(1,2)=((2.*x*y-2.*x*z-4.*x+y*z+y-z)*(y+1.)*(z-1.))/8.
             der(1,3)=((2.*x*y+2.*x*z+4.*x-y*z-y-z)*(y-1.)*(z-1.))/8.
             der(1,4)=((2.*x*y-2.*x*z-4.*x-y*z-y+z)*(y+1.)*(z-1.))/8.
             der(1,5)=-((2.*x*y-2.*x*z+4.*x-y*z+y-z)*(y-1.)*(z+1.))/8.
             der(1,6)=-((2.*x*y+2.*x*z-4.*x-y*z+y+z)*(y+1.)*(z+1.))/8.
             der(1,7)=-((2.*x*y-2.*x*z+4.*x+y*z-y+z)*(y-1.)*(z+1.))/8.
             der(1,8)=-((2.*x*y+2.*x*z-4.*x+y*z-y-z)*(y+1.)*(z+1.))/8.
             der(1,9)=-(y+1.)*(y-1.)*(z-1.)*x  ;   der(1,10)=(y+1.)*(y-1.)*(z+1.)*x
             der(1,11)=-(y-1.)*(z+1.)*(z-1.)*x ;   der(1,12)=(y+1.)*(z+1.)*(z-1.)*x
             der(1,13)=-((y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
             der(1,14)=((y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
             der(2,1)=((2.*x*y+x*z+x+2.*y*z+4.*y+z)*(x-1.)*(z-1.))/8.
             der(2,2)=((2.*x*y-x*z-x+2.*y*z+4.*y-z)*(x-1.)*(z-1.))/8.
             der(2,3)=((2.*x*y+x*z+x-2.*y*z-4.*y-z)*(x+1.)*(z-1.))/8.
             der(2,4)=((2.*x*y-x*z-x-2.*y*z-4.*y+z)*(x+1.)*(z-1.))/8.
             der(2,5)=-((2.*x*y-x*z+x-2.*y*z+4.*y-z)*(x-1.)*(z+1.))/8.
             der(2,6)=-((2.*x*y+x*z-x-2.*y*z+4.*y+z)*(x-1.)*(z+1.))/8.
             der(2,7)=-((2.*x*y-x*z+x+2.*y*z-4.*y+z)*(x+1.)*(z+1.))/8.
             der(2,8)=-((2.*x*y+x*z-x+2.*y*z-4.*y-z)*(x+1.)*(z+1.))/8.
             der(2,9)=-(x+1.)*(x-1.)*(z-1.)*y
             der(2,10)=(x+1.)*(x-1.)*(z+1.)*y
             der(2,11)=-((x+1.)*(x-1.)*(z+1.)*(z-1.))/2.
             der(2,12)=((x+1.)*(x-1.)*(z+1.)*(z-1.))/2.
             der(2,13)=-(x-1.)*(z+1.)*(z-1.)*y
             der(2,14)=(x+1.)*(z+1.)*(z-1.)*y
             der(3,1)=((x*y+2.*x*z+x+2.*y*z+y+4.*z)*(x-1.)*(y-1.))/8.
             der(3,2)=((x*y-2.*x*z-x+2.*y*z+y-4.*z)*(x-1.)*(y+1.))/8.
             der(3,3)=((x*y+2.*x*z+x-2.*y*z-y-4.*z)*(x+1.)*(y-1.))/8.
             der(3,4)=((x*y-2.*x*z-x-2.*y*z-y+4.*z)*(x+1.)*(y+1.))/8.
             der(3,5)=-((x*y-2.*x*z+x-2.*y*z+y-4.*z)*(x-1.)*(y-1.))/8.
             der(3,6)=-((x*y+2.*x*z-x-2.*y*z+y+4.*z)*(x-1.)*(y+1.))/8.
             der(3,7)=-((x*y-2.*x*z+x+2.*y*z-y+4.*z)*(x+1.)*(y-1.))/8.
             der(3,8)=-((x*y+2.*x*z-x+2.*y*z-y-4.*z)*(x+1.)*(y+1.))/8.
             der(3,9)=-((x+1.)*(x-1.)*(y+1.)*(y-1.))/2.
             der(3,10)=((x+1.)*(x-1.)*(y+1.)*(y-1.))/2.
             der(3,11)=-(x+1.)*(x-1.)*(y-1.)*z  ; der(3,12)=(x+1.)*(x-1.)*(y+1.)*z
             der(3,13)=-(x-1.)*(y+1.)*(y-1.)*z  ; der(3,14)=(x+1.)*(y+1.)*(y-1.)*z
           case(20)
             xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
             etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
             zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
             do l=1,20
               xi0=xi*xii(l); eta0=eta*etai(l); zeta0=zeta*zetai(l)
               if(l==4.or.l==8.or.l==16.or.l==20) then
                   der(1,l)=-.5*xi*(1.+eta0)*(1.+zeta0)
                   der(2,l)=.25*etai(l)*(1.-xi*xi)*(1.+zeta0)
                   der(3,l)=.25*zetai(l)*(1.-xi*xi)*(1.+eta0)
               else if(l>=9.and.l<=12)then
                   der(1,l)=.25*xii(l)*(1.-eta*eta)*(1.+zeta0)
                   der(2,l)=-.5*eta*(1.+xi0)*(1.+zeta0)
                   der(3,l)=.25*zetai(l)*(1.+xi0)*(1.-eta*eta)
               else if(l==2.or.l==6.or.l==14.or.l==18) then
                   der(1,l)=.25*xii(l)*(1.+eta0)*(1.-zeta*zeta)
                   der(2,l)=.25*etai(l)*(1.+xi0)*(1.-zeta*zeta)
                   der(3,l)=-.5*zeta*(1.+xi0)*(1.+eta0)
               else
                   der(1,l)=.125*xii(l)*(1.+eta0)*(1.+zeta0)*(2.*xi0+eta0+zeta0-1.)
                   der(2,l)=.125*etai(l)*(1.+xi0)*(1.+zeta0)*(xi0+2.*eta0+zeta0-1.)
                   der(3,l)=.125*zetai(l)*(1.+xi0)*(1.+eta0)*(xi0+eta0+2.*zeta0-1.)
               end if
             end do 
           case default
             print*,"wrong number of nodes in shape_der"        
         end select
       case default
         print*,"wrong number of dimensions in shape_der"
     end select
     return
end subroutine dpshape_der
!
!
subroutine dpinvert(a,ndim,d)
     implicit none
     integer :: ndim
     double precision,intent(in out)::a(ndim,ndim)
     double precision,allocatable:: b(:,:)
     double precision::d
     integer::i,j,k,n,np
     n= ubound(a,1)
     allocate(b(n,n))
     b=a
     if(n.eq.1) then
       d=1.0
       if(d.eq.0.0) then 
          write(11,*) "determinante matriz jacobiano igual a zero"
          stop
       endif
       a(1,1)=1/b(1,1)
     endif
     if(n.eq.2) then
       d=b(1,1)*b(2,2)-b(1,2)*b(2,1)
       if(d.eq.0.0) then 
          write(11,*) "determinante matriz jacobiano igual a zero"
          stop
       endif   
       a(1,1)=b(2,2)/d
       a(1,2)=-b(1,2)/d
       a(2,1)=-b(2,1)/d
       a(2,2)=b(1,1)/d
     endif
     if(n.eq.3) then
       d=0
       d=d+b(1,1)*(b(3,3)*b(2,2)-b(3,2)*b(2,3))
       d=d-b(2,1)*(b(3,3)*b(1,2)-b(3,2)*b(1,3))
       d=d+b(3,1)*(b(2,3)*b(1,2)-b(2,2)*b(1,3))
       if(d.eq.0.0) then 
          write(11,*) "determinante matriz jacobiano igual a zero"
          stop
       endif   
       a(1,1)=+(b(3,3)*b(2,2)-b(3,2)*b(2,3))/d
       a(1,2)=-(b(3,3)*b(1,2)-b(3,2)*b(1,3))/d
       a(1,3)=+(b(2,3)*b(1,2)-b(2,2)*b(1,3))/d
       a(2,1)=-(b(3,3)*b(2,1)-b(3,1)*b(2,3))/d
       a(2,2)=+(b(3,3)*b(1,1)-b(3,1)*b(1,3))/d
       a(2,3)=-(b(2,3)*b(1,1)-b(2,1)*b(1,3))/d
       a(3,1)=+(b(3,2)*b(2,1)-b(3,1)*b(2,2))/d
       a(3,2)=-(b(3,2)*b(1,1)-b(3,1)*b(1,2))/d
       a(3,3)=+(b(2,2)*b(1,1)-b(2,1)*b(1,2))/d
     endif   
     if(n.ge.4) then
       write(11,*) "Ordem da matriz é superior a 3, dpinvert precisar ser reprogramada"
       stop
     endif
     return         
end subroutine dpinvert
!
!
subroutine dpbeemat(bee,nod,nst,ndof,ndim,deriv)
   ! bee matrix for 2-d elasticity or elastoplasticity (ih=3 or 4 respectively)
   ! or for 3-d (ih = 6)
     implicit none
     integer :: nod,nst,ndof,ndim
     double precision,intent(in)::deriv(ndim,nod)
     double precision,intent(out)::bee(nst,ndof)
   ! local variables
     integer::k,l,m,n , ih
     double precision::x,y,z
     bee=0.
     ih = ubound(bee,1)
     nod = ubound(deriv,2)
     select case (ih)
       case(3,4)
        do m=1,nod
           k=2*m; l=k-1; x=deriv(1,m); y=deriv(2,m)
           bee(1,l)=x; bee(3,k)=x; bee(2,k)=y; bee(3,l)=y
        end do
       case(6)
        do m=1,nod
           n=3*m;  k=n-1; l=k-1
           x=deriv(1,m); y=deriv(2,m); z=deriv(3,m)
           bee(1,l)=x; bee(4,k)=x; bee(6,n)=x
           bee(2,k)=y; bee(4,l)=y; bee(5,n)=y
           bee(3,n)=z; bee(5,k)=z; bee(6,l)=z
        end do
       case default
        print*,'wrong dimension for nst in bee matrix'        
     end select   
     return
end subroutine dpbeemat
!
!
subroutine dpecmat(ecm,fun,nod,ndof,nodof)
     implicit none
     integer,intent(in):: nodof,ndof,nod
     double precision,intent(in)::fun(nod)
     double precision,intent(out)::ecm(ndof,ndof)
     integer:: i,j
     double precision::nt(ndof,nodof),tn(nodof,ndof)
       nt = .0; tn = .0
       do i = 1 , nod ; do j = 1 , nodof
          nt((i-1)*nodof+j,j) = fun(i); tn(j,(i-1)*nodof+j) = fun(i)
       end do; end do
       ecm = matmul( nt , tn )
       return
end subroutine dpecmat
!
!
subroutine dpfsparv(a,b,c,bk,km,g,kdiag)
   ! assembly of element matrices into skyline global matrix
     !implicit none
     integer :: a,b,c
     double precision::km(b,b)
     integer::g(b),kdiag(c)
     double precision::bk(a)
     integer::i,idof,k,j,iw,ival
     idof=ubound(g,1)
     do i=1,idof
        k=g(i)
        if(k/=0) then
          do j=1,idof
            if(g(j)/=0) then
               iw=k-g(j)
               if(iw>=0) then
                   ival=kdiag(k)-iw
                   bk(ival)=bk(ival)+km(i,j) 
                end if
            end if
          end do
        end if
     end do
     return
end subroutine dpfsparv
!
!
