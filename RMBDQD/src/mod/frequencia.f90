! Subrotina que, dado um E, encontra uma frequencia correspondente
subroutine firstfreqE(nfile, e, frequency)

    implicit none
    character(len=50) :: nfile
    character(len=12) :: mode
    integer::nels,neq,nn,nr,nip,nodof=3,nod=4,nst=6,ndof,i,k,iel,ndim=3,ijk
    integer::nband,fab,band
    integer::j, aa
    double precision:: e,v,rho,det,xsinal, frequency
    character(len=15):: element='tetrahedron'

    real :: rdescarte
    integer :: idescarte
    integer :: nex
    integer nprint, OS
    !logical :: fexist

    !----------------------------- dynamic arrays--------------------------------

    integer,allocatable:: nf(:,:),g_num(:,:),kdiag(:),num(:),g(:),g_g(:,:),kdiagB(:)
    !
    double precision,allocatable :: g_coord(:,:),kbv(:),mbv(:),points(:,:),weights(:),dee(:,:), &
        fun(:),der(:,:),jac(:,:),bee(:,:),deriv(:,:),km(:,:),emm(:,:), &
        ecm(:,:),coord(:,:),XHELP(:)
    double precision, allocatable:: R(:,:),EIGV(:)
    integer :: NC,NROOT
    integer:: NNC,NWK,NWM,NNM
    double precision:: RTOL
    integer:: NITEM,IFSS,IFPR,NSTIF,IOUT  
    double precision, allocatable:: TT(:),W(:),AR(:),BR(:),VEC(:,:),D(:),&
        RTOLV(:),BUP(:),BLO(:),BUPC(:)

    !--------------------------input and initialisation---------------------------

    nprint = 1
    call get_OS(OS)

    if (nprint == 1) then
    write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
        " 0.00%"
    endif

    open (unit=10,file=trim(nfile),status='old',action='read')
    open (unit=11,file='Modal.res',status='replace',action='write')
    open (unit=16,file='jacobi.res',status='replace',action='write')
    read(10,*) mode
    read(10,*) idescarte, rdescarte, rdescarte, rdescarte, NROOT, rdescarte, idescarte
    read (10,*) nels,nn,nip,v,rho 

    ndof=nod*nodof !numero de nos x graus de liberdade  = numero de equacoes por elemento   

    allocate (nf(nodof,nn),g_coord(ndim,nn),g_num(nod,nels),num(nod),g(ndof),g_g(ndof,nels))
    allocate (points(nip,ndim),weights(nip),dee(nst,nst),fun(nod),der(ndim,nod),jac(ndim,ndim))
    allocate (bee(nst,ndof),deriv(ndim,nod),km(ndof,ndof),emm(ndof,ndof),&
        ecm(ndof,ndof),coord(nod,ndim))
    read(10,*) g_coord
    read(10,*) g_num
    nf=1
    read(10,*) nr
    if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
    !
    call formnf(nn,nodof,nf)
    neq=maxval(nf)
    allocate(kdiag(neq))

    do i=1, neq
        kdiag(i) = 0
    enddo

    close(10)

    if (nprint == 1) then
    write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
        " 5.00%"
    endif

    !-------- loop the elements to find nband and set up global arrays ------------
    !
    nband=0

    do iel =1,nels
    !
    fab=ubound(nf,1)
    num=g_num(:,iel)
    call num_to_g(iel,nels,num,nn,nod,nodof,ndof,nf,g)
    g_g(:,iel)=g
    g_num(:,iel)=num
    !
    call fkdiag(ndof,neq,kdiag,g)

    !
    band= maxval(g,1,g>0)-minval(g,1,g>0)
    if(nband<band) then
        nband=band
    end if
    !

    if (mod(iel,500) == 0) then
        if (OS == 1) call execute_command_line("echo ' ' > Modal.res")
    endif

    if (nprint == 1) then
    write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
        5.0 + 35.0*(real(iel)/real(nels)), "%" ! 40%
    endif

    enddo

    if (OS == 1) then
        call execute_command_line("echo ' ' > Modal.res")
    endif

    kdiag(1)=1
    do i=2,neq
        kdiag(i)=kdiag(i)+kdiag(i-1)

    if (nprint == 1) then
        write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
            40.0 + 5.0*(real(i)/real(neq)), "%" ! 45%
    endif

    end do
    !
    !
    aa=kdiag(neq)

    allocate(kbv(aa),mbv(aa))
    allocate (kdiagB(neq+1),XHELP(aa)) ! alterado XHELP(neq)=>XHELP(aa)

    !
    kbv=0.0
    mbv=0.0
    !
    call dpsample(element,ndim,nip,points,weights)
    call dpdeemat(dee,nst,e,v)

    !----------------- element stiffness integration and assembly------------------

    do iel=1,nels
        num = g_num(:,iel)
        g = g_g(:12,iel ) 
        coord =transpose(g_coord(:,num))
        km=0.0 
        emm=0.0   
        do i=1,nip
            call dpshape_fun(fun,ndim,nip,nod,points,i)
            call dpshape_der(der,ndim,nod,nip,points,i)
            jac=matmul(der,coord) 
            call dpinvert(jac,ndim,det) 
            deriv = matmul(jac,der)
            call dpbeemat(bee,nod,nst,ndof,ndim,deriv)
            km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
            call dpecmat(ecm,fun,nod,ndof,nodof)
            emm=emm+ecm*det*weights(i)*rho
        end do
        xsinal=0.0
        do i=1,ndof
            xsinal=xsinal+emm(i,i)/dabs(emm(i,i))
        enddo  
        if(xsinal.lt.0.0) emm=-emm
        xsinal=0.0
        do i=1,ndof
            xsinal=xsinal+km(i,i)/dabs(km(i,i))
        enddo  
        if(xsinal.lt.0.0) km=-km    
        !write(11,*) "aa",  aa
        call dpfsparv (aa,ndof,neq,kbv,km,g,kdiag)  
        call dpfsparv (aa,ndof,neq,mbv,emm,g,kdiag)
        !

    if (nprint == 1) then
        write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
            45.0 + 22.0*(.02+ real(iel)/real(nels)), "%"! 68%
    endif

    enddo 
    !
    ! forming the kdiag, kbv and mbv  consistent with SSPACE Subroutine
    !   
    kdiagB(1)=1
    kdiagB(2)=2
    do i=3,neq+1
        kdiagB(i)=kdiag(i-1)+1
    enddo   
    !
    !
    do i=2,neq
        ijk=0
        do j=kdiag(i),kdiag(i-1)+1,-1
            ijk=ijk+1
            XHELP(ijk)=kbv(j) 
            if (ijk.gt.1) then
                if(dabs(XHELP(ijk)/XHELP(1)).LT.1D-10) XHELP(ijk)=0.0D0
            endif   
        enddo
        ijk=0
        do j=kdiagB(i),kdiagB(i+1)-1
            ijk=ijk+1
            kbv(j)=XHELP(ijk)
        enddo
    enddo    
    !
    !  
    do i=2,neq
        ijk=0
        do j=kdiag(i),kdiag(i-1)+1,-1
            ijk=ijk+1
            XHELP(ijk)=mbv(j)
            if (ijk.gt.1) then
                if(dabs(XHELP(ijk)/XHELP(1)).LT.1D-10) XHELP(ijk)=0.0D0
            endif
        enddo
        ijk=0
        do j=kdiagB(i),kdiagB(i+1)-1
            ijk=ijk+1
            mbv(j)=XHELP(ijk)
        enddo
    enddo

    if (nprint == 1) then
    write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
        70.0, "%"! 70%
    endif

    deallocate(kdiag,num,g_num,weights,g_coord,km,nf, points, &
        dee,coord,jac,der,deriv,g,g_g,fun,ecm,emm,XHELP)

    !------------------------------find eigenvalues--------------------------------


    !NROOT=7   -> agora é um dado de entrada
    nex = NROOT
    NC=MIN(2*NROOT,NROOT+8)
    NNC=NC*(NC+1)/2
    NWK=kdiagB(NEQ+1)-1
    NWM=kdiagB(NEQ+1)-1
    !! write(*,*) "NWK=",NWK
    NNM=NEQ+1
    RTOL=1.0D-6
    NITEM=16
    IFSS=0
    IFPR=0
    NSTIF=3
    IOUT=11

    allocate(EIGV(NC),R(NEQ,NC))
    allocate(TT(NEQ),W(NEQ),AR(NNC),BR(NNC),VEC(NC,NC),D(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC))

    call SPACE (kbv,mbv,kdiagB,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ,NNM,NWK, &
        NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)

    !print*, "after space"

    frequency = DSQRT(EIGV(nex))/(2.0D0*3.141592654D0)
    !print*, "Frequency: ", frequency

    deallocate(kdiagB, kbv, mbv,bee)
    deallocate(TT, W, AR, BR, VEC, D, RTOLV, BUP, BLO, BUPC)
    deallocate(R, EIGV)

    close(11, status='delete')
    close(16, status='delete')

    if (nprint == 1) then
    write(*,fmt="(a1)", advance="no") achar(13)
    endif

endsubroutine firstfreqE

subroutine firstfreqPoisson(nfile, v, frequency)
    implicit none

    character(len=50) :: nfile
    character(len=12) :: mode
    integer::nels,neq,nn,nr,nip,nodof=3,nod=4,nst=6,ndof,i,k,iel,ndim=3,ijk
    integer::nband,fab,band
    integer::j,aa
    double precision:: e,v,rho,det,xsinal, frequency
    character(len=15):: element='tetrahedron'
    logical :: fexist

    real :: rdescarte
    integer :: idescarte, nex, OS
    integer :: nprint

    !----------------------------- dynamic arrays--------------------------------

    integer,allocatable:: nf(:,:),g_num(:,:),kdiag(:),num(:),g(:),g_g(:,:),kdiagB(:)
    !
    double precision,allocatable :: g_coord(:,:),kbv(:),mbv(:),points(:,:),weights(:),dee(:,:), &
        fun(:),der(:,:),jac(:,:),bee(:,:),deriv(:,:),km(:,:),emm(:,:), &
        ecm(:,:),coord(:,:),XHELP(:)
    double precision, allocatable:: R(:,:),EIGV(:)
    integer :: NC,NROOT
    integer:: NNC,NWK,NWM,NNM
    double precision:: RTOL
    integer:: NITEM,IFSS,IFPR,NSTIF,IOUT  
    double precision, allocatable:: TT(:),W(:),AR(:),BR(:),VEC(:,:),D(:),RTOLV(:),&
        BUP(:),BLO(:),BUPC(:)

    !--------------------------input and initialisation---------------------------                         

    nprint = 1
    call get_OS(OS)

    if (nprint == 1) then
    write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
        " 0.00%"
    endif

    open (unit=10,file=trim(nfile),status='old',action='read')
    open (unit=11,file='Modal.res',status='replace',action='write')
    open (unit=16,file='jacobi.res',status='replace',action='write')
    read(10,*) mode
    read(10,*) idescarte, rdescarte, rdescarte, rdescarte, NROOT, rdescarte, idescarte
    read (10,*) nels,nn,nip,e,rho 

    ndof=nod*nodof !numero de nos x graus de liberdade  = numero de equacoes por elemento   

    allocate (nf(nodof,nn),g_coord(ndim,nn),g_num(nod,nels),num(nod),g(ndof),g_g(ndof,nels))
    allocate (points(nip,ndim),weights(nip),dee(nst,nst),fun(nod),der(ndim,nod),jac(ndim,ndim))
    allocate (bee(nst,ndof),deriv(ndim,nod),km(ndof,ndof),emm(ndof,ndof),ecm(ndof,ndof),&
        coord(nod,ndim))
    read(10,*) g_coord
    read(10,*) g_num
    nf=1
    read(10,*) nr
    if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
    !
    call formnf(nn,nodof,nf)
    neq=maxval(nf)
    allocate(kdiag(neq))

    do i=1, neq
        kdiag(i) = 0
    enddo

    close(10)

    if (nprint == 1) then
    write(*,fmt="(a1,t6,a)", advance="no") achar(13), &
        " 5.00%"
    endif

    !-------- loop the elements to find nband and set up global arrays ------------
    !
    nband=0

    do iel =1,nels
    !
    fab=ubound(nf,1)
    num=g_num(:,iel)
    call num_to_g(iel,nels,num,nn,nod,nodof,ndof,nf,g)
    g_g(:,iel)=g
    g_num(:,iel)=num
    !
    call fkdiag(ndof,neq,kdiag,g)

    !
    band= maxval(g,1,g>0)-minval(g,1,g>0)
    if(nband<band) then
        nband=band
    end if
    !

    if (mod(iel,500) == 0) then
        if (OS == 1 ) call execute_command_line("echo ' ' > Modal.res")
    endif

    if (nprint == 1) then
    write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
        5.0 + 35.0*(real(iel)/real(nels)), "%" ! 40%
    endif

    enddo

    if (OS == 1) then
        call execute_command_line("echo ' ' > Modal.res")
    endif

    kdiag(1)=1
    do i=2,neq
        kdiag(i)=kdiag(i)+kdiag(i-1)

        if (nprint == 1) then
            write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
                40.0 + 5.0*(real(i)/real(neq)), "%" ! 45%
        endif
    end do
    !
    !
    aa=kdiag(neq)

    allocate(kbv(aa),mbv(aa))
    allocate (kdiagB(neq+1),XHELP(aa)) ! alterado XHELP(neq)=>XHELP(aa)

    !
    kbv=0.0
    mbv=0.0
    !
    call dpsample(element,ndim,nip,points,weights)
    call dpdeemat(dee,nst,e,v)

    !----------------- element stiffness integration and assembly------------------

    do iel=1,nels
        num = g_num(:,iel)
        g = g_g(:12,iel ) 
        coord =transpose(g_coord(:,num))
        km=0.0 
        emm=0.0   
        do i=1,nip
            call dpshape_fun(fun,ndim,nip,nod,points,i)
            call dpshape_der(der,ndim,nod,nip,points,i)
            jac=matmul(der,coord) 
            call dpinvert(jac,ndim,det) 
            deriv = matmul(jac,der)
            call dpbeemat(bee,nod,nst,ndof,ndim,deriv)
            km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
            call dpecmat(ecm,fun,nod,ndof,nodof)
            emm=emm+ecm*det*weights(i)*rho
        end do
        xsinal=0.0
        do i=1,ndof
            xsinal=xsinal+emm(i,i)/dabs(emm(i,i))
        enddo  
        if(xsinal.lt.0.0) emm=-emm
        xsinal=0.0
        do i=1,ndof
            xsinal=xsinal+km(i,i)/dabs(km(i,i))
        enddo  
        if(xsinal.lt.0.0) km=-km    
        !write(11,*) "aa",  aa
        call dpfsparv (aa,ndof,neq,kbv,km,g,kdiag)  
        call dpfsparv (aa,ndof,neq,mbv,emm,g,kdiag)
        !

    if (nprint == 1) then
        write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
            45.0 + 22.0*(.02+ real(iel)/real(nels)), "%"! 68%
    endif

    enddo 
    !
    ! forming the kdiag, kbv and mbv  consistent with SSPACE Subroutine
    !   
    kdiagB(1)=1
    kdiagB(2)=2
    do i=3,neq+1
        kdiagB(i)=kdiag(i-1)+1
    enddo   
    !
    !
    do i=2,neq
        ijk=0
        do j=kdiag(i),kdiag(i-1)+1,-1
            ijk=ijk+1
            XHELP(ijk)=kbv(j) 
            if (ijk.gt.1) then
                if(dabs(XHELP(ijk)/XHELP(1)).LT.1D-10) XHELP(ijk)=0.0D0
            endif   
        enddo
        ijk=0
        do j=kdiagB(i),kdiagB(i+1)-1
            ijk=ijk+1
            kbv(j)=XHELP(ijk)
        enddo
    enddo    
    !
    !  
    do i=2,neq
        ijk=0
        do j=kdiag(i),kdiag(i-1)+1,-1
            ijk=ijk+1
            XHELP(ijk)=mbv(j)
            if (ijk.gt.1) then
                if(dabs(XHELP(ijk)/XHELP(1)).LT.1D-10) XHELP(ijk)=0.0D0
            endif
        enddo
        ijk=0
        do j=kdiagB(i),kdiagB(i+1)-1
            ijk=ijk+1
            mbv(j)=XHELP(ijk)
        enddo
    enddo

    if (nprint == 1) then
    write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
        70.0, "%"! 70%
    endif

    deallocate(kdiag,num,g_num,weights,g_coord,km,nf, points, &
        dee,coord,jac,der,deriv,g,g_g,fun,ecm,emm,XHELP)

    !------------------------------find eigenvalues--------------------------------

    !NROOT=7
    nex = NROOT
    NC=MIN(2*NROOT,NROOT+8)
    NNC=NC*(NC+1)/2
    NWK=kdiagB(NEQ+1)-1
    NWM=kdiagB(NEQ+1)-1
    !! write(*,*) "NWK=",NWK
    NNM=NEQ+1
    RTOL=1.0D-6
    NITEM=16
    IFSS=0
    IFPR=0
    NSTIF=3
    IOUT=11

    allocate(EIGV(NC),R(NEQ,NC))
    allocate(TT(NEQ),W(NEQ),AR(NNC),BR(NNC),VEC(NC,NC),D(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC))

    call SPACE (kbv,mbv,kdiagB,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ,NNM,NWK, &
        NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)

    frequency = DSQRT(EIGV(nex))/(2.0D0*3.141592654D0)

    deallocate(kdiagB, kbv, mbv,bee)
    deallocate(TT, W, AR, BR, VEC, D, RTOLV, BUP, BLO, BUPC)
    deallocate(R, EIGV)

    close(11, status='delete')
    close(16, status='delete')

    if (nprint == 1) then
    write(*,fmt="(a1,t6,f0.2,a)", advance="no") achar(13), &
        100.0, "%"! 100%
    write(*,fmt="(a1)", advance="no") achar(13)
    endif

    return
    
end subroutine firstfreqPoisson
