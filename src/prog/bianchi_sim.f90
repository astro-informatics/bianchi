!------------------------------------------------------------------------------
! bianchi_sim -- BIANCHI simulation program
!
!! Simulate a Bianchi VII_h model of the CMB.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 June 2005
!
! Revisions:
!   June 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

program bianchi_sim

  use s2_types_mod
  use s2_sky_mod, only: S2_SKY_FILE_TYPE_MAP, S2_SKY_FILE_TYPE_SKY
  use bianchi_sky_mod
  use bianchi_error_mod
  use pix_tools, only: nside2npix

  use extension, only: getArgument, nArguments
  use paramfile_io, only: paramfile_handle, parse_init, parse_int, &
    parse_real, parse_double, parse_lgt, parse_string, concatnl

  implicit none

  character(len=S2_STRING_LEN) :: filename_param
  character(len=S2_STRING_LEN) :: description
  character(len=S2_STRING_LEN) :: line
  type(paramfile_handle) :: handle

  real(s2_dp), parameter :: OMEGA0_LOWER = 0d0
  real(s2_dp), parameter :: OMEGA0_UPPER = 1d0
  real(s2_dp), parameter :: OMEGA0_DEFAULT = 0.5d0
  real(s2_dp), parameter :: X_LOWER = 0.06d0
  real(s2_dp), parameter :: X_UPPER = 10d0
  real(s2_dp), parameter :: X_DEFAULT = 0.55d0
  real(s2_dp), parameter :: ZE_LOWER = 1d2
  real(s2_dp), parameter :: ZE_UPPER = 1d4
  real(s2_dp), parameter :: ZE_DEFAULT = 1d3
  real(s2_dp), parameter :: S12H_DEFAULT = 1d0
  real(s2_dp), parameter :: S13H_DEFAULT = 1d0
  real(s2_sp), parameter :: FWHM_DEFAULT = 330d0
  logical, parameter :: RHAND_DEFAULT = .true.
  character(len=*), parameter :: BIANCHI_SKY_QUAD_DIRECT_STR = 'direct'
  character(len=*), parameter :: BIANCHI_SKY_QUAD_QTRAP_STR = 'qtrap'
  character(len=*), parameter :: BIANCHI_SKY_QUAD_QSIMP_STR = 'qsimp'
  character(len=*), parameter :: INIT_TYPE_REAL = 'real'
  character(len=*), parameter :: INIT_TYPE_ALM = 'alm'
  character(len=*), parameter :: FILE_TYPE_MAP_STR = 'map'
  character(len=*), parameter :: FILE_TYPE_SKY_STR = 'sky'

  character(len=S2_STRING_LEN) :: filename_out
  character(len=S2_STRING_LEN) :: filetype_str = FILE_TYPE_MAP_STR
  integer :: filetype = S2_SKY_FILE_TYPE_MAP

  type(bianchi_sky) :: b
  real(s2_dp) :: omega0, x, zE, s12H, s13H
  integer :: nside, N, quad_AB, quad_IAB, lmax
  ! Default angles in radians for user input, but converted to radians later.
  real(s2_sp) :: alpha=0e0, beta=-90e0, gamma=0e0 
  logical :: rhand = .true.
  logical :: apply_beam = .false.
  real(s2_sp) :: fwhm = FWHM_DEFAULT
  character(len=S2_STRING_LEN) :: quad_str = BIANCHI_SKY_QUAD_QSIMP_STR
  character(len=S2_STRING_LEN) :: init_type = INIT_TYPE_REAL

  ! Set default parameter values.
  filename_out = 'sky.fits'
  omega0 = 0.5d0
  x = 0.55d0
  zE = 1d3
  s12H = 1d0
  s13H = 1d0
  nside = 64
  lmax = 64
  N = 100
  rhand = .true.


  write(*,'(a)') '**********************************************'
  write(*,'(a)') 'BIANCHI VII_h rotating universe CMB simulation'
  write(*,'(a)') 'Jason McEwen                    September 2005'
  write(*,'(a)') 'mcewen@mrao.cam.ac.uk                         '
  write(*,'(a)') '**********************************************'


  !---------------------------------------
  ! Parse parameters
  !---------------------------------------
 
  ! Initialise file parser.
  if(nArguments() == 0) then
     filename_param = ''
  else
     if(nArguments() /= 1) then
        call bianchi_error(BIANCHI_ERROR_SIM_NARG, 'bianchi_sim', &
             comment_add='Usage: bianchi_sim [input parameter filename]')
     end if
     call getArgument(1, filename_param)
  end if
  handle = parse_init(trim(filename_param))

  ! Get omega0.
  write(line,'(a,f4.1,a,f4.1,a)') '(In the range ', &
       OMEGA0_LOWER, ' < omega0 < ', OMEGA0_UPPER, ')'
  description = concatnl('', &
       'Enter omega0: ', &
       line)
1 continue
  omega0 = parse_double(handle, 'omega0', &
       default=OMEGA0_DEFAULT, descr=description)
  if(omega0 <=  OMEGA0_LOWER .or. omega0 >=  OMEGA0_UPPER) then
     if(handle%interactive) goto 1
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='omega0 invalid')

  end if
      
  ! Get x.
  write(line,'(a,f4.1,a,f4.1,a)') '(In the range ', &
       X_LOWER, ' <= x <= ', X_UPPER, ')'
  description = concatnl('', &
       'Enter x: ', &
       line)
2 continue
  x = parse_double(handle, 'x', &
       default=X_DEFAULT, descr=description)
  if(x <  X_LOWER .or. x >  X_UPPER) then
     if(handle%interactive) goto 2
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='x invalid')

  end if

 ! Get zE.
  write(line,'(a,e8.2,a,e8.2,a)') '(In the range ', &
       ZE_LOWER, ' <= zE <= ', ZE_UPPER, ')'
  description = concatnl('', &
       'Enter zE: ', &
       line)
3 continue
  zE = parse_double(handle, 'zE', &
       default=ZE_DEFAULT, descr=description)
  if(zE <  ZE_LOWER .or. zE >  ZE_UPPER) then
     if(handle%interactive) goto 3
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='zE invalid')

  end if

  ! Get s12H.
  description = concatnl('', &
       'Enter s12H: ')
4 continue
  s12H = parse_double(handle, 's12H', &
       default=S12H_DEFAULT, descr=description)
  if(s12H < 0d0) then
     if(handle%interactive) goto 4
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='s12H invalid')
  end if

  ! Get s13H.
  description = concatnl('', &
       'Enter s13H: ')
5 continue
  s13H = parse_double(handle, 's13H', &
       default=S13H_DEFAULT, descr=description)
  if(s13H < 0d0) then
     if(handle%interactive) goto 5
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='s13H invalid')
  end if

  ! Get right-handedness status.
  description = concatnl('', &
       'Enter right-handedness status (logical): ')
  rhand = parse_lgt(handle, 'rhand', &
       default=RHAND_DEFAULT, descr=description)

  ! Get alpha.
  description = concatnl('', &
       'Enter alpha (degrees): ')
6 continue
  alpha = parse_real(handle, 'alpha', &
       default=alpha, descr=description)
  if(alpha < -360d0 .or. alpha > 360) then
     if(handle%interactive) goto 6
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='alpha invalid')
  end if
  alpha = alpha / 180e0 * pi

  ! Get beta.
  description = concatnl('', &
       'Enter beta (degrees): ')
7 continue
  beta = parse_real(handle, 'beta', &
       default=beta, descr=description)
  if(beta < -180d0 .or. beta > 180) then
     if(handle%interactive) goto 7
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='beta invalid')
  end if
  beta = beta / 180e0 * pi

  ! Get gamma.
  description = concatnl('', &
       'Enter gamma (degrees): ')
8 continue
  gamma = parse_real(handle, 'gamma', &
       default=gamma, descr=description)
  if(gamma < -360d0 .or. gamma > 360) then
     if(handle%interactive) goto 8
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='gamma invalid')
  end if
  gamma = gamma / 180e0 * pi

  ! Get init_type, i.e. harmonic or real space.
  description = concatnl('', &
       'Enter initialisation type (init_type={real; alm}): ')
9  continue
  init_type = parse_string(handle, 'init_type', &
       default=trim(init_type), descr=description)
! write(*,*) 'trim(init_type): ', trim(init_type)
! write(*,*) 'trim(INIT_TYPE_REAL): ', trim(INIT_TYPE_REAL)
! write(*,*) 'trim(INIT_TYPE_ALM): ', trim(INIT_TYPE_ALM)

  if(trim(init_type) /= INIT_TYPE_REAL .and. trim(init_type) /= INIT_TYPE_ALM) then
     if(handle%interactive) goto 9
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='Init type string invalid')
  end if

  ! Get apply_beam.
  description = concatnl('', &
    'Enter apply_beam status (logical): ')
  apply_beam = parse_lgt(handle, 'apply_beam', &
    default=apply_beam, descr=description)

  if(apply_beam) then
    ! Get beam fwhm.
    description = concatnl('', &
       'Enter beam fwhm (arcmin): ')
15  continue
    fwhm = parse_real(handle, 'fwhm', &
        default=FWHM_DEFAULT, descr=description)
    if(fwhm <= 0d0) then
      if(handle%interactive) goto 15
      call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
            comment_add='fwhm invalid')
    end if
  end if

  ! Get lmax if compute bianchi model in alm space or applying beam.
  if(trim(init_type) == INIT_TYPE_ALM .or. apply_beam) then
     description = concatnl("", &
          "Enter the maximum harmonic l (lmax) for the simulated sky: ")
10   continue
     lmax = parse_int(handle, 'lmax', &
          default=lmax, descr=description)
     if(lmax < 0) then
        if(handle%interactive) goto 10
        call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
             comment_add='lmax invalid')
     endif
  end if

  ! Get nside.
  description = concatnl("", &
       "Enter the resolution parameter (nside) for the simulated sky: ", &
       "(npix = 12*nside**2, where nside must be a power of 2)")
11 continue
  nside = parse_int(handle, 'nside', &
       default=nside, descr=description)
  if(nside2npix(nside) < 0) then
     if(handle%interactive) goto 11
     call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='nside invalid')
  endif
  
  ! Get quadrature AB type.
  description = concatnl('', &
       'Enter quadrature type (quad_AB={direct; qtrap; qsimp}): ')
12  continue
  quad_str = parse_string(handle, 'quad_AB', &
       default=trim(quad_str), descr=description)

  ! Set quad integer status.
  select case (quad_str)
     
    case (BIANCHI_SKY_QUAD_DIRECT_STR)
       quad_AB = BIANCHI_SKY_QUAD_DIRECT

    case (BIANCHI_SKY_QUAD_QTRAP_STR)
      quad_AB = BIANCHI_SKY_QUAD_QTRAP

    case (BIANCHI_SKY_QUAD_QSIMP_STR)
      quad_AB = BIANCHI_SKY_QUAD_QSIMP

    case default
       if(handle%interactive) goto 12
       call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, 'bianchi_sim')

  end select

  ! Get quadrature IAB type.
  if(trim(init_type) == INIT_TYPE_ALM) then
    quad_str = BIANCHI_SKY_QUAD_DIRECT_STR   ! Reset default.
    description = concatnl('', &
      'Enter quadrature type (quad_IAB={direct; qtrap; qsimp}): ')
13  continue
    quad_str = parse_string(handle, 'quad_IAB', &
        default=trim(quad_str), descr=description)

    ! Set quad integer status.
    select case (quad_str)
     
      case (BIANCHI_SKY_QUAD_DIRECT_STR)
        quad_IAB = BIANCHI_SKY_QUAD_DIRECT

      case (BIANCHI_SKY_QUAD_QTRAP_STR)
        quad_IAB = BIANCHI_SKY_QUAD_QTRAP

      case (BIANCHI_SKY_QUAD_QSIMP_STR)
        quad_IAB = BIANCHI_SKY_QUAD_QSIMP

      case default
        if(handle%interactive) goto 13
        call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, 'bianchi_sim')

    end select
  end if

    
  if(quad_AB == BIANCHI_SKY_QUAD_DIRECT .or. quad_IAB == BIANCHI_SKY_QUAD_DIRECT) then

     ! Get N for direct quadrature.
     description = concatnl("", &
          'Enter the numerical integration resolution N: ', &
          '(N = number of terms in approximation)')
14   continue
     N = parse_int(handle, 'N', &
          default=N, descr=description)
     if(N < 2) then
        if(handle%interactive) goto 14
        call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
             comment_add='N too small')
     endif

  end if
  
  ! Get filename_out.
  description = concatnl('', &
       'Enter filename_out: ')
  filename_out = parse_string(handle, 'filename_out', &
       default=trim(filename_out), descr=description)

  ! Get output file type: map or sky.
  description = concatnl('', &
       'Enter output file type (filetype={map; sky}): ')
16  continue
  filetype_str = parse_string(handle, 'filetype', &
       default=trim(filetype_str), descr=description)

  ! Set filetype integer status.
  select case (filetype_str)
     
    case (FILE_TYPE_MAP_STR)
      filetype = S2_SKY_FILE_TYPE_MAP

    case (FILE_TYPE_SKY_STR)
      filetype = S2_SKY_FILE_TYPE_SKY

    case default
       if(handle%interactive) goto 16
       call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
         comment_add='Invalid output file type')

  end select


  !---------------------------------------
  ! Run simulation and save sky
  !---------------------------------------

  ! Simulated bianchi sky.
  write(*,'(a)')

  ! Initialise with specified method.
  select case (trim(init_type))
     
    case (INIT_TYPE_REAL)
       write(*,'(a)') 'Computing BIANCHI simulation in real space...'
       b = bianchi_sky_init(omega0, x, zE, s12H, s13H, rhand, nside, quad_AB, N)

    case (INIT_TYPE_ALM)
       write(*,'(a)') 'Computing BIANCHI simulation in harmonic space...'
       b = bianchi_sky_init_alm(omega0, x, zE, s12H, s13H, rhand, lmax, &
         quad_AB, quad_IAB, N, alpha, beta, gamma)

    case default
       call bianchi_error(BIANCHI_ERROR_SIM_PARAM_INVALID, 'bianchi_sim', &
          comment_add='Init type string invalid')

  end select

  write(*,'(a)') 'Simulation complete'
  write(*,'(a)')

  ! Perform rotation if computing simulation in real space (rotation performed
  ! in init routine if compute simulaiton in harmonic space).
  if(trim(init_type) == INIT_TYPE_REAL) then
    call bianchi_sky_rotate(b, alpha, beta, gamma)
  end if 
  
  ! Apply beam if required.
  if(apply_beam) then
      
    ! Recompute alms if necessary (if not necessary does nothing).
    ! Note, if compute alms directly but then perform
    ! bianchi_sky_rotate, alms are removed and must be recomputed.
    call bianchi_sky_compute_alm(b, lmax, lmax)

    ! If map already present then apply beam will recompute the map.
    call bianchi_sky_apply_beam(b, fwhm, lmax)
    
  end if

  ! Compute map if computed bianchi simulation in alm space.
  if(trim(init_type) == INIT_TYPE_ALM) then
    call bianchi_sky_compute_map(b, nside)
  end if
  
  ! Save sky.
  call bianchi_sky_write(b, filename_out, filetype)
  write(*,'(a,a)') 'Simulated map written to ', trim(filename_out)
  write(*,'(a)')

  ! Write bianchi simulation variables to standard output.
  write(*,'(a)') 'Simulation parameters:'
  call bianchi_sky_param_write(b)

  ! Free bianchi object.
  call bianchi_sky_free(b)

end program bianchi_sim
