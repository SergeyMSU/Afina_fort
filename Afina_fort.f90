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
include "Geometry.f90"


program Afina_fort
    use STORAGE
    use GEOMETRY
    implicit none

    integer :: N
    N = 400


    call Init_Setka(gl_S_in, 2 * N)
    call Init_Setka(gl_S_out, N)
    gl_S_in%area = 1
    gl_S_out%area = 2
	call ALL_konstruct(N)

    ! print*, "u = ", gl_S_out%u
    ! print*, "un = ", gl_S_out%un

    call Set_Matrix(gl_S_in)
    call Set_Matrix(gl_S_out)
    ! call Print_matrix_real(gl_S_in%MM)
    call Culc_equ(gl_S_in)
    call Culc_equ(gl_S_out)
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

