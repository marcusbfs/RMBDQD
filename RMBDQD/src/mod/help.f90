module help

contains

subroutine help_options()

    implicit none

    character(len=*), parameter :: prog = 'Nome do programa (a definir)'
    character(len=*), parameter :: t = '    '
    
    write(*,'(a)') 'SYNOPSIS'
    write(*,'(3a)') t, prog, ' [ModelFile.dat <txt, ASCII>]'

    print*,
    write(*,'(a)') 'DESCRIPTION'
    write(*,'(2a)') t, "This program takes an input file that contains physical, geometrical and"
    write(*,'(2a)') t, "mesh information about a body and produces through a Finite Element Method"
    write(*,'(2a)') t, "simulation the Young's modulus or the Poisson's ratio related to those"
    write(*,'(2a)') t, "properties, focusing on the fundamental frequency."

    print*,
    write(*,'(a)') 'INPUT MODEL FILE'
    write(*,'(2a)') t, "=========================.txt ====================================="
    print*,
    write(*,'(2a)') t, "find_elastic/find_poisson"
    print*,
    write(*,'(2a)') t, "1/2/3    R1    R2    fundamental_frequency(Hz)     modal_frequency_integer"
    print*,
    write(*,'(2a)') t, "tolerance    max_iterations_value"
    print*,
    write(*,'(2a)') t, "n_elements    n_nodes    n_GaussianQuadrature_points"
    print*,
    write(*,'(2a)') t, "poisson_ratio/young_modulus(Pa)    volumetric_mass_density(kg/m^3)"
    print*,
    print*,
    print*,
    write(*,'(2a)') t, "x1            y1            z1"
    write(*,'(2a)') t, "x2            y2            z2"
    write(*,'(2a)') t, ".             .             . "
    write(*,'(2a)') t, ".             .             . "
    write(*,'(2a)') t, ".             .             . "
    write(*,'(2a)') t, "xn_nodes-1    yn_nodes-1    zn_nodes-1"
    write(*,'(2a)') t, "xn_nodes      yn_nodes      zn_nodes"
    print*,
    print*,
    write(*,'(2a)') t, "conectivity_information (n_elements)"
    print*,
    print*,
    write(*,'(2a)') t, "n_boundary_conditions"
    write(*,'(2a)') t, "node_number    x_val    y_val    z_val  (n_boundary_conditions)"
    print*,
    write(*,'(2a)') t, "======================= end .txt ==================================="
    
    stop

end subroutine help_options

endmodule help
