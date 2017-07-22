module commonUtils
    contains

        subroutine delete_file(nfile)
            implicit none
            character(len=*) :: nfile
            logical :: fexist
            integer :: u

            inquire(file=nfile, exist=fexist)
            if (fexist .eqv. .true.) then
                open(newunit=u, file=nfile)
                close(u,status='delete')
            endif
            
        end subroutine delete_file

        function checkPoisson(v) result (r)
            implicit none
            double precision :: v
            logical :: r

            if (v < 0.5d0 .and. v > -1.0d0 ) then
                r = .true.
            else
                r = .false.
            endif

        endfunction checkPoisson

        function lower( string ) result (new) ! função para colocar a palavra em letras minúsculas.
                                              ! Retirada de um forum
            character(len=*)           :: string

            character(len=len(string)) :: new
            integer                    :: i
            integer                    :: k, length

            length = len(string)
            new    = string
            do i = 1,len(string)
            k = iachar(string(i:i))
            if ( k >= iachar('A') .and. k <= iachar('Z') ) then
                k = k + iachar('a') - iachar('A')
                new(i:i) = achar(k)
            endif
            enddo
        end function lower 

end module
