!  Afina_fort.f90 
!
!  FUNCTIONS:
!  Afina_fort - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Afina_fort
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
include "Storage.f90"
include "For_splayn.f90"
include "Geometry.f90"


program Afina_fort
    use STORAGE
    use GEOMETRY
    implicit none

    integer :: N, i
    N = 700 ! 400

    par_N = N
    call Init_Setka(gl_S_in, 2 * N)
    call Init_Setka(gl_S_out, N)
    gl_S_in%area = 1
    gl_S_out%area = 2
	call ALL_konstruct(N)

    call Read_surface(7)
    call Print_gran_setka(gl_S_in)

    do i = 1, 300
        print*, "step = ", i
        call Calc_grans(gl_S_in)
        call Calc_grans(gl_S_out)

        call Set_Matrix(gl_S_in)
        call Set_Matrix(gl_S_out)
        ! call Print_matrix_real(gl_S_in%MM)
        call Culc_equ(gl_S_in)
        call Culc_equ(gl_S_out)
        call Move_surface()
    end do


    call Calc_grans(gl_S_in)
    call Calc_grans(gl_S_out)
    call Set_Matrix(gl_S_in)
    call Set_Matrix(gl_S_out)
    call Culc_equ(gl_S_in)
    call Culc_equ(gl_S_out)
    call Save_surface(8)
    call Print_Pressure()


    ! call Test_Splayn()


    ! call Print_test()

    call Print_Solution()

    call Print_yzel_all()
	call Print_gran_setka(gl_S_in)
	call Print_gran_setka(gl_S_out)
	

    ! Variables

    ! Body of Afina_fort

    



    print *, 'Hello World'
    ! print*, gl_S_in%N
    ! print*, gl_S_in%u


    call Dell_Setka(gl_S_in)
    call Dell_Setka(gl_S_out)
	
	pause

end program Afina_fort

