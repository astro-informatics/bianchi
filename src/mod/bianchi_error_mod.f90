!------------------------------------------------------------------------------
! bianchi_error_mod  -- BIANCHI library error class
!
!! Functionality to handle errors that may occur in the bianchi library. 
!! Public bianchi error codes are defined, with corresponding private error 
!! comments and default halt execution status.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 June 2005
!
! Revisions:
!   June 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

module bianchi_error_mod

  use s2_types_mod, only: S2_STRING_LEN

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: bianchi_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: BIANCHI_ERROR_NUM = 12

  integer, public, parameter :: &
    BIANCHI_ERROR_NONE = 0, &
    BIANCHI_ERROR_INIT = 1, &
    BIANCHI_ERROR_NOT_INIT = 2, &
    BIANCHI_ERROR_INIT_FAIL = 3, &
    BIANCHI_ERROR_MEM_ALLOC_FAIL = 4, &
    BIANCHI_ERROR_SKY_QUAD_INVALID = 5, &
    BIANCHI_ERROR_SKY_QUAD_STEP_EXCEED = 6, &
    BIANCHI_ERROR_SKY_N_MISSING = 7, &
    BIANCHI_ERROR_SIM_PARAM_INVALID = 8, &
    BIANCHI_ERROR_SIM_NARG = 9, &
    BIANCHI_ERROR_PLM1TABLE_THETA_INVALID = 10, &
    BIANCHI_ERROR_PLM1TABLE_L_INVALID = 11

  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.
  !! Comment associated with each error type.
  character(len=S2_STRING_LEN), parameter :: &
    error_comment(BIANCHI_ERROR_NUM) = &
      (/ & 
      'No error                                                                 ', &
      'Attempt to initialise object that has already been initialised           ', &
      'Object not initialised                                                   ', &
      'Object initialisation failed                                             ', &
      'Memory allocation failed                                                 ', &
      'Quadrature type invalid                                                  ', &
      'Exceeded number of steps limit whem computing quadrature                 ', &
      'Must specify N if using direct quadrature                                ', &
      'Invalid simulation parameter                                             ', &
      'Invalid number of command line parameters                                ', &
      'Invalid theta for Plm1 lookup table (not on regular grid)                ', &
      'Invalid l for Plm1 lookup table (out of range)                           ' &
      /) 
  
  !! Default program halt status of each error type.
  logical, parameter :: &
    halt_default(BIANCHI_ERROR_NUM) = &
      (/ &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true., & 
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.  /)
  
  
  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! bianchi_error
    !
    !! Display error message corresponding to error_code and halt program 
    !! execution if required.
    !!
    !! Variables:
    !!   - error_code: Integer error code.
    !!   - [procedure]: Procedure name where bianchi_error called from.  Displayed 
    !!     when error message printed to screen.
    !!   - [comment_add]: If present, additional comment to append to default 
    !!     error comment.
    !!   - [comment_out]: If present the error comment is copied to comment_out
    !!     on output.
    !!   - [halt_in]: If present overrides default halt value.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_error(error_code, procedure, comment_add, &
      comment_out, halt_in)

      integer, intent(in) :: error_code
      character(len=*), intent(in), optional :: procedure, comment_add
      character(len=*), intent(inout), optional :: comment_out
      logical, intent(in), optional :: halt_in

      logical :: halt
      character(len=*), parameter :: comment_prefix = 'BIANCHI_ERROR: '

      !---------------------------------------
      ! Display error message
      !---------------------------------------

      if(present(procedure)) then

        if(present(comment_add)) then
	  write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            '''', &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            ''''
        end if
 
     else

        if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
            trim(error_comment(error_code+1)), &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
        end if

      end if

      ! Copy error comment if comment_out present.
      if(present(comment_out)) comment_out = error_comment(error_code+1)

      !---------------------------------------
      ! Halt program execution if required
      !---------------------------------------
      
      if( present(halt_in) ) then
        halt = halt_in
      else
        halt = halt_default(error_code+1)
      end if

      if( halt ) then
        write(*,'(a,a,a,a,a)') comment_prefix, &
          '  Halting program execution ', &
          'due to error ''', trim(error_comment(error_code+1)), ''''
        stop
      end if

    end subroutine bianchi_error


end module bianchi_error_mod
