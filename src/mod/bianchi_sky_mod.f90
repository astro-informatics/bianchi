!------------------------------------------------------------------------------
! bianchi_sky_mod -- BIANCHI library sky class
!
!! Provides functionality to simulate a Bianchi VII_h model of the CMB.
!! Uses the s2_sky module to create healpix sky maps.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 June 2005
!
! Revisions:
!   June 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module bianchi_sky_mod

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod
  use bianchi_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    bianchi_sky_init, &
    bianchi_sky_init_alm, &
    bianchi_sky_free, &
    bianchi_sky_write, &
    bianchi_sky_rotate, &
    bianchi_sky_param_write, &
    bianchi_sky_compute_map, &
    bianchi_sky_compute_alm, &
    bianchi_sky_get_alm, &
    bianchi_sky_apply_beam


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Quadrature type: Direct (rectangular).
  integer, public, parameter :: BIANCHI_SKY_QUAD_DIRECT = 1

  !! Quadrature type: Trapezium rule.
  integer, public, parameter :: BIANCHI_SKY_QUAD_QTRAP = 2

  !! Quadrature type: Simpson's rule.
  integer, public, parameter :: BIANCHI_SKY_QUAD_QSIMP = 3

  ! Logical to disable the lookup table
  logical, parameter :: BIANCHI_DISABLE_PLM1TABLE = .false.

  ! CMB T to convert Delta_T / T map to just Delta_T map.
#ifdef MILLIK
  ! produce maps in units of mK.
  real(s2_dp), parameter :: BIANCHI_CMB_T = 2.725d3
#else
  ! Produce Delta_T / T maps.
  real(s2_dp), parameter :: BIANCHI_CMB_T = 1d0 
#endif


  !---------------------------------------
  ! Data types
  !---------------------------------------

  !! - init: Initialisation status.
  !! - omega0: Density (input parameter).
  !! - x: Related to characteristic wavelength over which basis vectors 
  !!   change orientation (input parameter).
  !! - zE: Red shift (input parameter).
  !! - s12H: Normalised shear 12 (normalised to Hubble constant)
  !!   (input parameter).
  !! - s13H: Normalised shear 13 (normalised to Hubble constant) 
  !!   (input parameter).
  !! - rhand: Logical specifying handedness of map (true=right).
  !! - wH: Normalised vorticity (normalised to Hubble constant).
  !! - h: Related to characteristic wavelength over which basis vectors 
  !!   change orientation.
  !! - tau0: Conformal time of photon reception.
  !! - tauE: Conformal time of photon emission.
  !! - sky: Simulated map.

  type, public :: bianchi_sky
     private
     logical :: init = .false.
     real(s2_dp) :: omega0 = 0.0d0
     real(s2_dp) :: x = 0.0d0
     real(s2_dp) :: zE = 0.0d0
     real(s2_dp) :: s12H = 0.0d0
     real(s2_dp) :: s13H = 0.0d0
     logical :: rhand = .true.
     real(s2_dp) :: wH = 0.0d0
     real(s2_dp) :: h = 0.0d0
     real(s2_dp) :: tau0 = 0.0d0
     real(s2_dp) :: tauE = 0.0d0
     type(s2_sky) :: sky
  end type bianchi_sky


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! bianchi_sky_init
    !
    !! Initialise bianchi object by performing a bianchi simulation, computed
    !! in real space, with the specified parameters.
    !!
    !! Variables:
    !!  - omega0: Input omega0 parameter (see bianchi data type for
    !!    explanation).
    !!  - x: Input x parameter (see bianchi data type for explanation).
    !!  - zE: Input zE parameter (see bianchi data type for explanation).
    !!  - s12H: Input s12H parameter (see bianchi data type for explanation).
    !!  - s13H: Input s13H parameter (see bianchi data type for explanation).
    !!  - rhand: Logical to specify handedness of map.
    !!  - nside: Nside of Healpix map to generate.
    !!  - quad: Quadrature rule to use.
    !!  - [N]: Resolution of numerical integration (number of terms used to
    !!    evaluate integral).  If not specified default of 10 is used.
    !!  - b: Initialised bianchi object with all parameters and simulated
    !!    map calculated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_init(omega0, x, zE, s12H, s13H, rhand, nside, &
      quad, N) result(b)

      use pix_tools, only: pix2ang_ring, nside2npix, in_ring 

      real(s2_dp), intent(in) :: omega0, x, zE, s12H, s13H
      integer, intent(in) :: nside
      logical, intent(in) :: rhand
      integer, intent(in) :: quad
      integer, intent(in), optional :: N
      type(bianchi_sky) :: b

      real(s2_dp) :: c1, c3, Atheta0, Btheta0
      real(s2_dp) :: theta0, phi0, thetaOB, phiOB
      integer :: npix, ipix, fail, Nuse = 10
      real(s2_dp), allocatable :: map(:)
      integer, parameter :: PERCENT_LOOP = 1000
      real(s2_sp) :: handedness_sign = +1d0

      real(s2_dp), parameter :: PHI_CENTRE = 0.0d0  ! Extract extire rings.
      real(s2_dp), parameter :: PHI_DIFF = PI
      integer :: iring, nring, max_pix_per_ring, npixring
      integer, allocatable :: ipixring(:)

      ! Check object not already initialised.
      if(b%init) then
        call bianchi_error(BIANCHI_ERROR_INIT, 'bianchi_sky_init')
        return
      end if

      ! Check N specified if using direct quadrature.
      if(quad == BIANCHI_SKY_QUAD_DIRECT &
          .and. .not. present(N)) then
        call bianchi_error(BIANCHI_ERROR_SKY_N_MISSING, 'bianchi_sky_init')
      end if
      
      ! Initialise parameters passed as arguments.
      b%omega0 = omega0
      b%x = x
      b%zE = zE
      b%s12H = s12H
      b%s13H = s13H
      b%rhand = rhand

      ! Compute other variables based on these parameters.
      b%h = b%x**2.0d0 * (1.0d0 - b%omega0)
      b%wH = sqrt(1d0+b%h) * sqrt(1d0+9d0*b%h) &
           * sqrt(b%s12H**2d0 + b%s13H**2d0) &
           / ( 6d0 * b%x**2d0 * b%omega0   ) 
      b%tau0 = 2.0d0 / sqrt(b%h) * asinh( sqrt(1/omega0 - 1.0d0) ) 
      b%tauE = 2.0d0 / sqrt(b%h) &
           * asinh(sqrt( (1/omega0 - 1.0d0)/(1.0d0 + b%zE) ) ) 
      c1 = bianchi_sky_comp_c1(b%omega0, b%x)
      c3 = bianchi_sky_comp_c3(b%omega0, b%h)

      ! Initialise healpix map settings.
      npix = nside2npix(nside)
      allocate(map(0:npix-1), stat=fail)
      map = 0d0
      if(fail /= 0) then
         call bianchi_error(BIANCHI_ERROR_MEM_ALLOC_FAIL, 'bianchi_sky_init')
      end if

      ! Set handedness sign.
      if(b%rhand) then
         handedness_sign = 1d0
      else
         handedness_sign = -1d0
      end if

      ! Set ring parameter values.
      nring = 4*nside - 1
      max_pix_per_ring = 4*nside

      ! Allocate space for ring pixel indices.
      allocate(ipixring(0:max_pix_per_ring-1), stat=fail)
      if(fail /= 0) then
         call bianchi_error(BIANCHI_ERROR_MEM_ALLOC_FAIL, 'bianchi_sky_init')
      end if
      ipixring = 0

      ! Compute map value for each pixel in each ring.
      ! Do ring at a time since A and B functions of theta so once need 
      ! to compute once for each ring.
      do iring = 1, nring    ! Checked this and ring is indeed indexed from 1.

         if(mod(iring,nring/4) == 0) then
            write(*,'(a,f5.1,a)') ' Percent complete: ', &
                 iring/real(nring,s2_sp)*100e0, '%'
         end if

         ! Get ring_pix in order to determine theta of ring.
         call in_ring(nside, iring, PHI_CENTRE, PHI_DIFF, ipixring, &
              & npixring, nest=0) ! Ring scheme!!!!!!!!!

         ! Calculate thetaOB for current ring.
         call pix2ang_ring(nside, ipixring(0), thetaOB, phiOB)

         ! Set photon theta angle from observation angle.
         theta0 = PI - thetaOB

         ! Compute A(theta) and B(theta)
         Atheta0 = bianchi_sky_comp_A(theta0, b%tauE, b%tau0, b%h, &
            b%zE, c1, c3, quad, N)
         Btheta0 = bianchi_sky_comp_B(theta0, b%tauE, b%tau0, b%h, &
            b%zE, c1, c3, quad, N)

         ! Compute map valus for all pixels in current ring.
         do ipix = 0,npixring-1

            ! Calculate thetaOB for current ring.
            call pix2ang_ring(nside, ipixring(ipix), thetaOB, phiOB)

            ! Set photon phi angle from observation angle.
            phi0 = phiOB - PI

            ! Compute map pixel value.
            map(ipixring(ipix)) = &
              ( b%s12H * Atheta0 + b%s13H * Btheta0 ) * sin(phi0) &
              + handedness_sign * ( b%s12H * Btheta0 - b%s13H * Atheta0 ) &
                * cos(phi0)

         end do

      end do

      write(*,'(a)') ' Percent complete: 100.0%'

      ! Convert Delta_T/T map computed to Delta_T map.
      map = BIANCHI_CMB_T * map

      ! Initialise sky object with map.
      b%sky = s2_sky_init(real(map,s2_sp), nside, S2_SKY_RING)

      ! Set initialised status.
      b%init = .true.

      ! Free memory.
      deallocate(map)
      deallocate(ipixring)

    end function bianchi_sky_init


    !--------------------------------------------------------------------------
    ! bianchi_sky_init_alm
    !
    !! Initialise bianchi object by performing a bianchi simulation, computed
    !! in harmonic space, with the specified parameters.
    !!
    !! Notes:
    !!   - Uses lookup table for Plms for N=100, quad_IA=direct.
    !     If try case quad_AB = direct, quad_IAB /= direct, N=100 will try
    !     to use lookup table when invalid.  Simply do not use this case
    !     (since it will not give accurate estimates for integrals either)
    !     or set the global variable to disable the loopup table and
    !     recompile.  Adding functionality to account for this case would
    !     involve passing quad_IAB to low levels of the program where it is
    !     otherwise not required and would introduce unnecessary complexity. 
    !        -> Actually fine now since use A_grid/B_grid to indicate whether
    !           using direct quadrature for IA/IB.
    !!
    !! Variables:
    !!  - omega0: Input omega0 parameter (see bianchi data type for
    !!    explanation).
    !!  - x: Input x parameter (see bianchi data type for explanation).
    !!  - zE: Input zE parameter (see bianchi data type for explanation).
    !!  - s12H: Input s12H parameter (see bianchi data type for explanation).
    !!  - s13H: Input s13H parameter (see bianchi data type for explanation).
    !!  - rhand: Logical to specify handedness of map.
    !!  - nside: Nside of Healpix map to generate.
    !!  - lmax: Maximum harmonic l to compute alms up to.
    !!  - quad_AB: Quadrature rule to use for AB integration.
    !!  - quad_IAB: Quadrature rule to use for IAB integration.
    !!  - [N]: Resolution of numerical integration (number of terms used to
    !!    evaluate integral).  If not specified default of 10 is used.
    !!  - [beam_fwhm]: Full-width-half-maximum of Gaussian beam to apply.
    !!    If not present then no beam is applied.
    !!    (Must be specified in arcmin.)
    !!  - b: Initialised bianchi object with all parameters and simulated
    !!    map calculated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_init_alm(omega0, x, zE, s12H, s13H, rhand, lmax, &
      quad_AB, quad_IAB, N, alpha, beta, gamma) result(b)

      use s2_dl_mod, only: s2_dl_beta_operator

      real(s2_dp), intent(in) :: omega0, x, zE, s12H, s13H
      integer, intent(in) :: lmax
      logical, intent(in) :: rhand
      integer, intent(in) :: quad_AB
      integer, intent(in) :: quad_IAB
      integer, intent(in), optional :: N
      real(s2_sp), intent(in), optional :: alpha, beta, gamma
      type(bianchi_sky) :: b

      complex(s2_spc), allocatable :: alm(:,:)
      
      real(s2_dp), pointer :: dl(:, :) => null() 
      integer :: fail, l, Nuse = 10
      real(s2_dp) :: IA, IB, c1, c3
      real(s2_dp) :: handedness_sign = 1d0, lsign = 1d0
      real(s2_dp), allocatable :: A_grid(:), B_grid(:)
      real(s2_dp) :: theta, theta_inc
      integer :: itheta

      real(s2_sp), parameter :: ZERO_TOL = 1e-4
      integer :: m
      complex(s2_spc), allocatable :: alm_rotate(:,:)
      complex(s2_dpc) :: Dm_p1, Dm_m1
      complex(s2_dpc) :: icmpx

      ! Set I=sqrt(-1) since can't seem to do when declare.
      icmpx = cmplx(0d0,1d0)

!      ** Temporary code for generating lookup table.
!               real(s2_dp) :: Plm1_table(1:64,0:99)
!               do l = 1,64
!                 do itheta = 0,99
!                   theta = itheta / 100.0d0 * pi
!                   Plm1_table(l,itheta) = plgndr(l,1,cos(theta))
!                   write(*,'(a,i2,a,i2,a,f20.16,a)') 'data PLM1_TABLE(', &
!                     l, ',',itheta, ') / ', Plm1_table(l,itheta), '/'
!                 end do
!               end do
!               stop
!      ** End of temporary code.

      ! Check N specified if using direct quadrature.
      if((quad_AB == BIANCHI_SKY_QUAD_DIRECT &
          .or. quad_IAB == BIANCHI_SKY_QUAD_DIRECT) &
          .and. .not. present(N)) then
        call bianchi_error(BIANCHI_ERROR_SKY_N_MISSING, 'bianchi_sky_init_alm')
      end if
      if(present(N)) Nuse = N

      ! Initialise parameters passed as arguments.
      b%omega0 = omega0
      b%x = x
      b%zE = zE
      b%s12H = s12H
      b%s13H = s13H
      b%rhand = rhand

      ! Compute other variables based on these parameters.
      b%h = b%x**2.0d0 * (1.0d0 - b%omega0)
      b%wH = sqrt(1d0+b%h) * sqrt(1d0+9d0*b%h) &
           * sqrt(b%s12H**2d0 + b%s13H**2d0) &
           / ( 6d0 * b%x**2d0 * b%omega0   ) 
      b%tau0 = 2.0d0 / sqrt(b%h) * asinh( sqrt(1/omega0 - 1.0d0) ) 
      b%tauE = 2.0d0 / sqrt(b%h) &
           * asinh(sqrt( (1/omega0 - 1.0d0)/(1.0d0 + b%zE) ) ) 
      c1 = bianchi_sky_comp_c1(b%omega0, b%x)
      c3 = bianchi_sky_comp_c3(b%omega0, b%h)

      ! Set handedness sign.
      if(b%rhand) then
         handedness_sign = 1d0
      else
         handedness_sign = -1d0
      end if

      ! Initialise healpix alm settings.
      allocate(alm(0:lmax,0:lmax), stat=fail)
      alm = cmplx(0e0, 0e0)
      if(fail /= 0) then
        call bianchi_error(BIANCHI_ERROR_MEM_ALLOC_FAIL, &
          'bianchi_sky_init_alm')
      end if

      ! If IAB integrals computed over regular grid then can precompute
      ! A(theta) and B(theta) terms.
      ! Shoudl be considerably fasted than using other type of
      ! quadrature since for the other case one must recompute
      ! A(theta) and B(theta) for each of the (unknown) thetas,
      if(quad_IAB == BIANCHI_SKY_QUAD_DIRECT) then

        allocate(A_grid(0:Nuse-1), stat=fail)
        allocate(B_grid(0:Nuse-1), stat=fail)
        if(fail /= 0) then
          call bianchi_error(BIANCHI_ERROR_MEM_ALLOC_FAIL, &
            'bianchi_sky_init_alm')
        end if
        A_grid = 0d0
        B_grid = 0d0
        theta = 0d0
        theta_inc = (pi - 0d0) / real(Nuse, s2_dp)

        do itheta = 0,Nuse-1
          A_grid(itheta) = bianchi_sky_comp_A(theta, b%tauE, b%tau0, b%h, &
            b%zE, c1, c3, quad_AB, N)
          B_grid(itheta) = bianchi_sky_comp_B(theta, b%tauE, b%tau0, b%h, &
            b%zE, c1, c3, quad_AB, N) 
          theta = theta + theta_inc
        end do

      end if
      
      ! Compute alms.
      ! Compute alms from 1, since zero for l=0 (no dc component).
      lsign  = +1d0        ! lsign has value (-1)**l for each l loop.
      do l = 1, lmax

         !write(*,'(a,i3,a,i3)') ' Computing l = ', l, ' of ', lmax

         ! Invert lsign.
         lsign = - lsign

         ! Compute integrals.
         if(quad_IAB == BIANCHI_SKY_QUAD_DIRECT) then
            ! Use precomputed A(theta) and B(theta).
            IA = bianchi_sky_comp_IA(b%tauE, b%tau0, b%h, &
              b%zE, c1, c3, l, quad_AB, quad_IAB, Nuse, A_grid)
            IB = bianchi_sky_comp_IB(b%tauE, b%tau0, b%h, &
              b%zE, c1, c3, l, quad_AB, quad_IAB, Nuse, B_grid) 
         else
            IA = bianchi_sky_comp_IA(b%tauE, b%tau0, b%h, &
              b%zE, c1, c3, l, quad_AB, quad_IAB, N)
            IB = bianchi_sky_comp_IB(b%tauE, b%tau0, b%h, &
              b%zE, c1, c3, l, quad_AB, quad_IAB, N)
         end if
         
         ! Compute alms for given l.  Only m=1 is non-zero.
         alm(l,1) = - lsign * pi &
              * cmplx(- handedness_sign * (b%s12H * IB - b%s13H * IA), &
                      (b%s12H * IA + b%s13H * IB))

      end do

      ! Convert Delta_T/T map computed to Delta_T map.
      alm = BIANCHI_CMB_T * alm

      ! Rotate alms if Euler angles present and at least one angle non-zero.
      if(present(alpha) .and. present(beta) .and. present(gamma) &
          .and. (abs(alpha) > ZERO_TOL .or. &
          abs(beta) > ZERO_TOL .or. abs(gamma) > ZERO_TOL) ) then
      
        allocate(dl(-lmax:lmax, -lmax:lmax), stat=fail)
        allocate(alm_rotate(0:lmax, 0:lmax), stat=fail)
        if(fail /= 0) then
          call bianchi_error(BIANCHI_ERROR_MEM_ALLOC_FAIL, &
            'bianchi_sky_init_alm')
        end if
        alm_rotate = 0e0
        
        ! Perform rotation in harmonic space, noting Bianchi alms only
        ! non-zero for m=+/-m.
        do l = 1,lmax   
  
          call s2_dl_beta_operator(dl, real(beta,s2_dp), l)

          do m = 0,l

            Dm_p1 = exp(-icmpx*m*alpha) * dl(m, 1) * exp(-icmpx*gamma)
            Dm_m1 = exp(-icmpx*m*alpha) * dl(m,-1) * exp( icmpx*gamma)
  
            alm_rotate(l,m) = - Dm_m1 * conjg(alm(l,1)) + Dm_p1 * alm(l,1)
  
          end do
  
        end do

        ! Initialise sky object with rotated alms.
        b%sky = s2_sky_init(alm_rotate, lmax, lmax)

        ! Deallocate memory used for rotating alms.
        deallocate(alm_rotate)
        deallocate(dl)

      else

        ! Initialise sky object with alms.
        b%sky = s2_sky_init(alm, lmax, lmax)
        
      end if

      ! Set initialised status.
      b%init = .true.

      ! Free memory.
      deallocate(alm)
      
      if(allocated(A_grid)) deallocate(A_grid)
      if(allocated(B_grid)) deallocate(B_grid)
      
    end function bianchi_sky_init_alm


    !--------------------------------------------------------------------------
    ! bianchi_sky_free
    !
    !! Free all memory associated with a bianchi object.
    !!
    !! Variables:
    !!  - b: Bianchi object to free.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_free(b)

      type(bianchi_sky), intent(inout) :: b

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_free')
      end if 

      ! Free sky.
      call s2_sky_free(b%sky)

      ! Reset attributes.
      b%omega0 = 0.0d0
      b%x = 0.0d0
      b%zE = 0.0d0
      b%s12H = 0.0d0
      b%s13H = 0.0d0
      b%rhand = .true.
      b%wH = 0.0d0
      b%h = 0.0d0
      b%tau0 = 0.0d0
      b%tauE = 0.0d0
      b%init = .false.

    end subroutine bianchi_sky_free


    !--------------------------------------------------------------------------
    ! bianchi_sky_param_write
    !
    !! Write parameters of the Bianchi simulation to standard output.
    !!
    !! Variables:
    !!  - b: Bianchi object containing parameter atrtributes to write.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_param_write(b)

      type(bianchi_sky), intent(in) :: b

      write(*,'(a,e12.5)') ' b%omega0: ', b%omega0
      write(*,'(a,e17.5)') ' b%x: ', b%x
      write(*,'(a,e16.5)') ' b%zE: ', b%zE
      write(*,'(a,e14.5)') ' b%s12H: ', b%s12H
      write(*,'(a,e14.5)') ' b%s13H: ', b%s13H
      write(*,'(a,e16.5)') ' b%wH: ', b%wH
      write(*,'(a,l13)') ' b%rhand: ', b%rhand
      write(*,'(a,e17.5)') ' b%h: ', b%h
      write(*,'(a,e14.5)') ' b%tau0: ', b%tau0
      write(*,'(a,e14.5)') ' b%tauE: ', b%tauE

    end subroutine bianchi_sky_param_write


    !--------------------------------------------------------------------------
    ! bianchi_sky_compute_map
    !
    !! Compute the map of the bianchi sky, assuming the alms are already
    !! defined.
    !!
    !! Variables:
    !!  - b: Bianchi object containing alms to compute map of.
    !!  - nside: Healpix nside to compute map at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_compute_map(b, nside)

      type(bianchi_sky), intent(inout) :: b
      integer, intent(in) :: nside

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_compute_map')
      end if

      ! Compute bianchi sky alms.
      call s2_sky_compute_map(b%sky, nside)

    end subroutine bianchi_sky_compute_map


    !--------------------------------------------------------------------------
    ! bianchi_sky_compute_alm
    !
    !! Compute the alm of the bianchi sky, assuming the map is already
    !! defined.
    !!
    !! Variables:
    !!  - b: Bianchi object containing sky to compute alms of.
    !!  - lmax: Maximum harmonic l to consider when computing alms.
    !!  - mmax: Maximum harmonic m to consider when computing alms.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_compute_alm(b, lmax, mmax)

      type(bianchi_sky), intent(inout) :: b
      integer, intent(in) :: lmax, mmax

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_compute_alm')
      end if

      ! Compute bianchi sky alms.
      call s2_sky_compute_alm(b%sky, lmax, mmax)

    end subroutine bianchi_sky_compute_alm


    !--------------------------------------------------------------------------
    ! bianchi_sky_get_alm
    !
    !! Get alms contained in a bianchi sky object.  Alms must already be
    !! computed.
    !!
    !! Notes:
    !!   - Error occurs if alms are not already computed.
    !!
    !! Variables:
    !!   - b: Bianchi object containing sky that in turn contained the alms
    !!     to get.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_get_alm(b, alm)

      type(bianchi_sky), intent(in) :: b
      complex(s2_spc), intent(out) :: alm(:,:)

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_get_alm')
      end if 

      ! Get alms from sky.
      call s2_sky_get_alm(b%sky, alm)

    end subroutine bianchi_sky_get_alm


    !--------------------------------------------------------------------------
    ! bianchi_sky_apply_beam
    !
    !! Apply Gaussian beam with specified FWHM.  (FWHM must be passed in
    !! arcmin.)
    !!
    !! Notes:
    !!   - Error occurs if alms are not already computed.
    !!
    !! Variables:
    !!   - b: Bianchi object containing sky that is to be comvolved with
    !!     the beam.
    !!   - fwhm: Gaussian beam FWHM to use (specified in arcmin).
    !!   - lmax: Maximum harmonic l to consider.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_apply_beam(b, fwhm, lmax)

      use s2_pl_mod
      use s2_vect_mod, only: s2_vect_arcmin_to_rad

      type(bianchi_sky), intent(inout) :: b
      real(s2_sp), intent(in) :: fwhm
      integer, intent(in) :: lmax

      real(s2_sp) :: fwhm_rad
      type(s2_pl) :: beam

      ! Convert beam_fwhm to radians.
      fwhm_rad = s2_vect_arcmin_to_rad(fwhm)

      ! Create beam.
      beam = s2_pl_init_guassian(fwhm_rad, lmax)

      ! Apply beam.
      call s2_sky_conv(b%sky, beam)

      ! Free beam.
      call s2_pl_free(beam)

    end subroutine bianchi_sky_apply_beam


    !--------------------------------------------------------------------------
    ! bianchi_sky_write
    !
    !! Write the bianchi simulated sky to a file. 
    !!
    !! Variables:
    !!  - b: Bianchi object to save sky of.
    !!  - filename: Name of output file.
    !!  - file_type: Type of output file, either fits map or sky file 
    !!    (see s2_sky_mod for more details).  Integer flag that may be 
    !!    either S2_SKY_FILE_TYPE_MAP or S2_SKY_FILE_TYPE_SKY.
    !!  - [comment]: Optional comment to append to file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_write(b, filename, file_type, comment)

      type(bianchi_sky), intent(inout) :: b
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_type
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_write')
      end if 

      select case(file_type)

         case(S2_SKY_FILE_TYPE_MAP)
            call s2_sky_write_map_file(b%sky, filename, comment)
!            call s2_sky_write_alm_file(b%sky, filename, comment)

         case(S2_SKY_FILE_TYPE_SKY)
            call s2_sky_io_fits_write(filename, b%sky, comment)

         case default
            call s2_error(S2_ERROR_SKY_FILE_INVALID, 'bianchi_sky_write', &
              comment_add='Invalid file type specifier')

      end select

    end subroutine bianchi_sky_write


    !--------------------------------------------------------------------------
    ! bianchi_sky_rotate
    !
    !! Rotate the simulated sky of the bianchi object.  
    !!
    !! Variables:
    !!  - b: Bianchi object containing sky to be rotated.
    !!  - alpha: Alpha Euler angle of the rotation.
    !!  - beta: Beta Euler angle of the rotation.
    !!  - gamma: Gamma Euler angle of the rotation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi_sky_rotate(b, alpha, beta, gamma)

      type(bianchi_sky), intent(inout) :: b
      real(s2_sp) :: alpha, beta, gamma

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi_error(BIANCHI_ERROR_NOT_INIT, 'bianchi_sky_rotate')
      end if 

      call s2_sky_rotate(b%sky, alpha, beta, gamma)

    end subroutine bianchi_sky_rotate


    !--------------------------------------------------------------------------
    ! Bianchi parameter functions
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_s
    !
    !! Compute s variable in bianchi sky simulation. 
    !!
    !! Variables:
    !!  - tau: Conformal time to evaluate s for.
    !!  - tau0: tau0 parameter (see bianchi data type for explanation).
    !!  - theta0: theta0 position to evaluate s for.
    !!  - h: h parameter (see bianchi data type for explanation).
    !!  - s: s function value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_s(tau, tau0, theta0, h) result(s)

      real(s2_dp), intent(in) :: tau, tau0, theta0, h
      real(s2_dp) :: s

      s = tan(theta0/2.0d0) * exp( -sqrt(h) * (tau-tau0) )

    end function bianchi_sky_comp_s


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_psi
    !
    !! Compute psi variable in bianchi sky simulation. 
    !!
    !! Variables:
    !!  - tau: Conformal time to evaluate psi for.
    !!  - tau0: tau0 parameter (see bianchi data type for explanation).
    !!  - theta0: theta0 position to evaluate psi for.
    !!  - h: h parameter (see bianchi data type for explanation).
    !!  - psi: psi function value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_psi(tau, tau0, theta0, h) result(psi)

      real(s2_dp), intent(in) :: tau, tau0, theta0, h
      real(s2_dp) :: psi

      psi = (tau-tau0) &
           - 1/sqrt(h) * log( sin(theta0/2.0d0)**2.0d0 + &
           exp(2*sqrt(h)*(tau-tau0))*cos(theta0/2.0d0)**2.0d0 )

    end function bianchi_sky_comp_psi


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_c1
    !
    !! Compute c1 variable in bianchi sky simulation.
    !!
    !! Variables:
    !!  - omega0: omega0 parameter (see bianchi data type for explanation).
    !!  - x: x parameter (see bianchi data type for explanation).
    !!  - c1: Value of constant c1 computed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_c1(omega0, x) result(c1)

      real(s2_dp), intent(in) :: omega0, x
      real(s2_dp) :: c1

      c1 = 1/(3.0d0 * omega0 * x)

    end function bianchi_sky_comp_c1


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_c2
    !
    !! Compute c2 variable in bianchi sky simulation.
    !!
    !! Variables:
    !!  - sE: sE parameter (value of s function computed at emission).
    !!  - zE: zE parameter (see bianchi data type for explanation).
    !!  - c2: Value of constant c2 computed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_c2(sE, zE) result(c2)

      real(s2_dp), intent(in) :: sE, zE
      real(s2_dp) :: c2

      c2 = 2.0d0 * sE * (1d0 + zE) / (1d0 + sE**2.0d0)

    end function bianchi_sky_comp_c2


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_c3
    !
    !! Compute c3 variable in bianchi sky simulation.
    !!
    !! Variables:
    !!  - omega0: omega0 parameter (see bianchi data type for
    !!    explanation).
    !!  - h: h parameter (see bianchi data type for
    !!    explanation).
    !!  - c3: Value of constant c3 computed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_c3(omega0, h) result(c3)

      real(s2_dp), intent(in) :: omega0, h
      real(s2_dp) :: c3

      c3 = 4.0d0 * sqrt(h) * (1.0d0-omega0)**(3.0d0/2.0d0) / (omega0**2.0d0)

    end function bianchi_sky_comp_c3


    !--------------------------------------------------------------------------
    ! Bianchi integrals and integrands
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_A
    !
    !! Compute A(theta) for Bianchi simulation.
    !!
    !! Variables:
    !!  - theta: theta value to eavluate B(theta) for.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals.
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - Atheta: Value of Atheta term evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_A(theta, tauE, tau0, h, zE, &
         c1, c3, quad, N) result(Atheta)

      real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: quad
      integer, optional, intent(in) :: N
      real(s2_dp) :: Atheta

      integer :: Nuse = 100, i
      real(s2_dp) :: IA = 0d0
      real(s2_dp) :: tau, tau_inc, dtau
      real(s2_dp) :: psiE, sE, c2

      ! Set integral resolution if present, else stick with default.
      if(present(N)) Nuse = N

      ! Compute integral using specified quadrature rule.
      select case (quad)

        case (BIANCHI_SKY_QUAD_DIRECT)

           ! Set discretisation width.
           tau_inc = (tau0 - tauE) / real(Nuse, s2_dp)
           dtau = tau_inc    ! Simple rectangular quadrature rule

           tau = tauE
           IA = 0.0d0

           do i = 1, Nuse

              IA = IA + bianchi_sky_integrand_A(tau, tau0, theta, h) * dtau

              tau = tau + tau_inc

           end do

        case (BIANCHI_SKY_QUAD_QTRAP)

           IA = qtrap_a(bianchi_sky_integrand_A, tau0, theta, &
                h, tauE, tau0)

        case (BIANCHI_SKY_QUAD_QSIMP)

           IA = qsimp_a(bianchi_sky_integrand_A, tau0, theta, &
                h, tauE, tau0)

        case default

           call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, &
                'bianchi_sky_comp_A')

      end select

      ! Compute variables that depend on theta.
      psiE = bianchi_sky_comp_psi(tauE, tau0, theta, h)
      sE   = bianchi_sky_comp_s(tauE, tau0, theta, h)
      c2   = bianchi_sky_comp_c2(sE, zE)

      ! Compute A(theta) term.
      Atheta = &
           c1 * ( sin(theta) &
                  - c2 * (cos(psiE) - 3d0*sqrt(h)*sin(psiE)) ) &
           + c3 * IA

    end function bianchi_sky_comp_A


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_B
    !
    !! Compute B(theta) for Bianchi simulation.
    !!
    !! Variables:
    !!  - theta: theta value to eavluate B(theta) for.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals.
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - Btheta: Value of Btheta term evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_B(theta, tauE, tau0, h, zE, &
         c1, c3, quad, N) result(Btheta)

      real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: quad
      integer, optional, intent(in) :: N
      real(s2_dp) :: Btheta

      integer :: Nuse = 100, i
      real(s2_dp) :: IB = 0d0
      real(s2_dp) :: tau, tau_inc, dtau
      real(s2_dp) :: psiE, sE, c2

      ! Set integral resolution if present, else stick with default.
      if(present(N)) Nuse = N

      ! Compute integral using specified quadrature rule.
      select case (quad)

        case (BIANCHI_SKY_QUAD_DIRECT)

           ! Set discretisation width.
           tau_inc = (tau0 - tauE) / real(Nuse, s2_dp)
           dtau = tau_inc    ! Simple rectangular quadrature rule

           tau = tauE
           IB = 0.0d0

           do i = 1, Nuse

              IB = IB + bianchi_sky_integrand_B(tau, tau0, theta, h) * dtau

              tau = tau + tau_inc

           end do

        case (BIANCHI_SKY_QUAD_QTRAP)

           IB = qtrap_a(bianchi_sky_integrand_B, tau0, theta, &
                h, tauE, tau0)

        case (BIANCHI_SKY_QUAD_QSIMP)

           IB = qsimp_a(bianchi_sky_integrand_B, tau0, theta, &
                h, tauE, tau0)

        case default

           call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, &
                'bianchi_sky_comp_B')

      end select

      ! Compute variables that depend on theta.
      psiE = bianchi_sky_comp_psi(tauE, tau0, theta, h)
      sE   = bianchi_sky_comp_s(tauE, tau0, theta, h)
      c2   = bianchi_sky_comp_c2(sE, zE)

      ! Compute B(theta) term.
      Btheta = &
           c1 * ( 3d0*sqrt(h)*sin(theta) &
                  - c2 * (sin(psiE) + 3d0*sqrt(h)*cos(psiE)) ) &
           - c3 * IB

    end function bianchi_sky_comp_B


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_IA
    !
    !! Compute IA_l integral required in Bianchi simulation.
    !!
    !! Variables:
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute IA_l for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals (required to give to routine used to compute
    !!     A(theta)).
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - IA: Value of IA_l integral evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_IA(tauE, tau0, h, zE, c1, c3, l, &
        quad_AB, quad_IAB, N, A_grid) result(IA)

      real(s2_dp), intent(in) :: tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB, quad_IAB
      integer, intent(in), optional :: N
      real(s2_dp), intent(in), optional :: A_grid(0:)
      real(s2_dp) :: IA

      integer :: Nuse = 10, i
      real(s2_dp) :: theta, theta_inc, dtheta

      IA = 0d0

       ! Set integral resolution if present, else stick with default.
      if(present(N)) Nuse = N

      ! Compute integral using specified quadrature rule.
      select case (quad_IAB)

        case (BIANCHI_SKY_QUAD_DIRECT)

          ! Set discretisation width.
          theta = 0d0
          theta_inc = (pi - 0d0) / real(Nuse, s2_dp)
          dtheta = theta_inc

          do i = 1, Nuse

            IA = IA + bianchi_sky_integrand_IA(theta, tauE, tau0, h, &
              zE, c1, c3, l, quad_AB, Nuse, A_grid) * dtheta

            theta = theta + theta_inc

           end do

        case (BIANCHI_SKY_QUAD_QTRAP)

           IA = qtrap_ia(bianchi_sky_integrand_IA, tauE, tau0, h, &
                zE, c1, c3, l, quad_AB, 0d0, real(pi,s2_dp), Nuse)

        case (BIANCHI_SKY_QUAD_QSIMP)

           IA = qsimp_ia(bianchi_sky_integrand_IA,tauE, tau0, h, &
                zE, c1, c3, l, quad_AB, 0d0, real(pi,s2_dp), Nuse)

        case default

           call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, &
                'bianchi_sky_comp_IA')

      end select

      IA = IA * sqrt(  ( 2d0*l+1d0 ) / real( 4d0*pi*l*(l+1d0), s2_dp ) ) 

    end function bianchi_sky_comp_IA


    !--------------------------------------------------------------------------
    ! bianchi_sky_comp_IB
    !
    !! Compute IB_l integral required in Bianchi simulation.
    !!
    !! Variables:
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute IB_l for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals (required to give to routine used to compute
    !!     B(theta)).
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - IB: Value of IB_l integral evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_comp_IB(tauE, tau0, h, zE, c1, c3, l, &
        quad_AB, quad_IAB, N, B_grid) result(IB)

      real(s2_dp), intent(in) :: tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB, quad_IAB
      integer, intent(in), optional :: N
      real(s2_dp), intent(in), optional :: B_grid(0:)
      real(s2_dp) :: IB

      integer :: Nuse = 10, i
      real(s2_dp) :: theta, theta_inc, dtheta

      IB = 0d0

       ! Set integral resolution if present, else stick with default.
      if(present(N)) Nuse = N

      ! Compute integral using specified quadrature rule.
      select case (quad_IAB)

        case (BIANCHI_SKY_QUAD_DIRECT)

          ! Set discretisation width.
          theta = 0d0
          theta_inc = (pi - 0d0) / real(Nuse, s2_dp)
          dtheta = theta_inc

          do i = 1, Nuse

            IB = IB + bianchi_sky_integrand_IB(theta, tauE, tau0, h, &
              zE, c1, c3, l, quad_AB, Nuse, B_grid) * dtheta

            theta = theta + theta_inc

           end do

        case (BIANCHI_SKY_QUAD_QTRAP)

           IB = qtrap_ia(bianchi_sky_integrand_IB, tauE, tau0, h, &
                zE, c1, c3, l, quad_AB, 0d0, real(pi,s2_dp), Nuse)

        case (BIANCHI_SKY_QUAD_QSIMP)

           IB = qsimp_ia(bianchi_sky_integrand_IB,tauE, tau0, h, &
                zE, c1, c3, l, quad_AB, 0d0, real(pi,s2_dp), Nuse)

        case default

           call bianchi_error(BIANCHI_ERROR_SKY_QUAD_INVALID, &
                'bianchi_sky_comp_IB')

      end select

      IB = IB * sqrt(  ( 2d0*l+1d0 ) / real( 4d0*pi*l*(l+1d0), s2_dp ) ) 

    end function bianchi_sky_comp_IB


    !--------------------------------------------------------------------------
    ! bianchi_sky_integrand_A
    !
    !! Compute integrand of A integral.
    !!
    !! Variables:
    !!  - tau: Conformal time to evaluate integrand for.
    !!  - tau0: tau0 parameter (see bianchi data type for explanation).
    !!  - theta0: theta0 position to evaluate integrand for.
    !!  - h: h parameter (see bianchi data type for explanation).
    !!  - integrand: Integrand value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_integrand_A(tau, tau0, theta0, h) result(integrand)

      real(s2_dp), intent(in) :: tau, tau0, theta0, h
      real(s2_dp) :: integrand

      real(s2_dp) :: s, psi

      s = bianchi_sky_comp_s(tau, tau0, theta0, h)
      psi = bianchi_sky_comp_psi(tau, tau0, theta0, h)

      integrand = s * (1d0-s**2d0) * sin(psi) &
                  / ( (1d0+s**2d0)**2d0 * sinh(sqrt(h)*tau/2d0)**4d0 )

    end function bianchi_sky_integrand_A


    !--------------------------------------------------------------------------
    ! bianchi_sky_integrand_B
    !
    !! Compute integrand of B integral.
    !!
    !! Variables:
    !!  - tau: Conformal time to evaluate integrand for .
    !!  - tau0: tau0 parameter (see bianchi data type for explanation).
    !!  - theta0: theta0 position to evaluate integrand for.
    !!  - h: h parameter (see bianchi data type for explanation).
    !!  - integrand: Integrand value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_integrand_B(tau, tau0, theta0, h) result(integrand)

      real(s2_dp), intent(in) :: tau, tau0, theta0, h
      real(s2_dp) :: integrand

      real(s2_dp) :: s, psi

      s = bianchi_sky_comp_s(tau, tau0, theta0, h)
      psi = bianchi_sky_comp_psi(tau, tau0, theta0, h)

      integrand = s * (1d0-s**2d0) * cos(psi) &
                  / ( (1d0+s**2d0)**2d0 * sinh(sqrt(h)*tau/2d0)**4d0 )

    end function bianchi_sky_integrand_B


    !--------------------------------------------------------------------------
    ! bianchi_sky_integrand_IA
    !
    !! Compute integrand of IA_l integral.
    !!
    !! Notes:
    !!   - Computation of this integrand required the computation of
    !!     A(theta) and the associated integral requried to compute A(theta).
    !!
    !! Variables:
    !!   - theta: theta angle to evaluate integrand for.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute IA_l for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals (required to give to routine used to compute
    !!     A(theta)).
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - integrand: Integrand value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_integrand_IA(theta, tauE, tau0, h, zE, &
      c1, c3, l, quad_AB, N, A_grid) result(integrand)
      
      use bianchi_plm1table_mod

      real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB
      integer, optional, intent(in) :: N
      real(s2_dp), optional, intent(in) :: A_grid(0:)
      real(s2_dp) :: integrand

      integer, parameter :: m = 1
      real(s2_dp) :: Atheta
      integer :: Nlocal = 0
      real(s2_dp) :: itheta_dp
      integer :: itheta
      real(s2_dp), parameter :: ZERO_TOL = 1d-4

      if(present(N)) Nlocal = N

      if(present(A_grid)) then

        ! Use precomputed Atheta.

        itheta_dp = theta / pi * real(Nlocal,s2_dp)
        itheta = nint(itheta_dp)

        ! Check itheta is an integer (within error limits).
        if(abs(itheta_dp - itheta) > ZERO_TOL) then
          call bianchi_error(BIANCHI_ERROR_PLM1TABLE_THETA_INVALID, &
            'bianchi_sky_integrand_IA')
        end if
      
        Atheta = A_grid(itheta)

      else
      
        ! Compute Atheta.
        Atheta = bianchi_sky_comp_A(theta, tauE, tau0, h, zE, &
          c1, c3, quad_AB, N)
          
      end if 

      if(l>0 .and. l<=BIANCHI_PLM1TABLE_LMAX &
          .and. Nlocal == BIANCHI_PLM1TABLE_NTHETA &
          .and. .not. BIANCHI_DISABLE_PLM1TABLE .and. present(A_grid)) then
        ! Ahha, lookup table case.
        integrand = Atheta * bianchi_plm1table_getval(l,theta) * sin(theta)
      else  
        integrand = Atheta * plgndr(l,m,cos(theta)) * sin(theta)
!        integrand = Atheta * hm_plm_rec(l,m,cos(theta)) * sin(theta)

      end if
    
    end function bianchi_sky_integrand_IA


    !--------------------------------------------------------------------------
    ! bianchi_sky_integrand_IB
    !
    !! Compute integrand of IB_l integral.
    !!
    !! Notes:
    !!   - Computation of this integrand required the computation of
    !!     B(theta) and the associated integral requried to compute B(theta).
    !!
    !! Variables:
    !!   - theta: theta angle to evaluate integrand for.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute IB_l for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integrals (required to give to routine used to compute
    !!     B(theta)).
    !!   - [N]:  Number of terms in use in direct quadrature.
    !!   - integrand: Integrand value evaluated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi_sky_integrand_IB(theta, tauE, tau0, h, zE, &
      c1, c3, l, quad_AB, N, B_grid) result(integrand)

      use bianchi_plm1table_mod

      real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB
      integer, optional, intent(in) :: N
      real(s2_dp), optional, intent(in) :: B_grid(0:)
      real(s2_dp) :: integrand
      
      integer, parameter :: m = 1
      real(s2_dp) :: Btheta
      integer :: Nlocal = 0
      real(s2_dp) :: itheta_dp
      integer :: itheta
      real(s2_dp), parameter :: ZERO_TOL = 1d-4

      if(present(N)) Nlocal = N

      if(present(B_grid)) then

        ! Use precomputed Btheta.

        itheta_dp = theta / pi * real(Nlocal,s2_dp)
        itheta = nint(itheta_dp)

        ! Check itheta is an integer (within error limits).
        if(abs(itheta_dp - itheta) > ZERO_TOL) then
          call bianchi_error(BIANCHI_ERROR_PLM1TABLE_THETA_INVALID, &
            'bianchi_sky_integrand_IB')
        end if
      
        Btheta = B_grid(itheta)

      else

        Btheta = bianchi_sky_comp_B(theta, tauE, tau0, h, zE, &
              c1, c3, quad_AB, N)
      
      end if
              
      if(l>0 .and. l<=BIANCHI_PLM1TABLE_LMAX .and. &
          Nlocal == BIANCHI_PLM1TABLE_NTHETA .and. &
          .not. BIANCHI_DISABLE_PLM1TABLE .and. present(B_grid)) then
        ! Ahha, lookup table case.
        integrand = Btheta * bianchi_plm1table_getval(l,theta) * sin(theta)
      else      
        integrand = Btheta * plgndr(l,m,cos(theta)) * sin(theta)
!        integrand = Btheta * hm_plm_rec(l,m,cos(theta)) * sin(theta)

      end if

    end function bianchi_sky_integrand_IB


    !--------------------------------------------------------------------------
    ! Numerical integration routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! trapzd_a
    !
    !! Computes nth stage of refinement of extended trapezoidal rule.
    !! Adapted from numerical recipes for the application at hand.
    !! For use in computing integrals required to compute A(theta) and
    !! B(theta).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     This routine computes the nth stage of refinement of an extended 
    !!     trapezoidal rule. func is input as the name of the function to be 
    !!     integrated between limits a and b, also input. When called with
    !!     n=1, the routine returns as s the crudest estimate of 
    !!     int_b^a f(x)dx. Subsequent calls with n=2,3,... (in that sequential
    !!     order) will improve the accuracy of s by adding 2n-2 additional
    !!     interior points. s should not be modified between sequential calls.
    !!   - Adapted for use in evaluating integrals required to compute
    !!     A *or* B, (i.e. also takes parameters required to compute these
    !!    integrals).
    !!   - For our application tau0 and b will be the same, but clear and more
    !!     general to program as two separate variables.
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - theta0: theta0 position to evaluate integral refinement for.
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - s: Value of integral computed to current order.
    !!   - n: Refinement order.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine trapzd_a(func, tau0, theta0, h, a, b, s, n)

      real(s2_dp), intent(in) :: tau0, theta0, h
      real(s2_dp), intent(IN) :: a, b
      real(s2_dp), intent(INOUT) :: s
      integer, intent(IN) :: n
      interface
         function func(tau, tau0, theta0, h) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: tau, tau0, theta0, h
           real(s2_dp) :: integrand
         end function func
      end interface

      real(s2_dp) :: del, fsum, x
      integer :: it, j

      if (n == 1) then
         s = 0.5d0 * (b-a)*(func(a, tau0, theta0, h) + func(b, tau0, theta0, h))
      else
         it = 2**(n-2)
         del = (b-a) / real(it, s2_dp)  
                             ! This is the spacing of the points to be added.
         x = a + 0.5d0*del
         fsum = 0d0
         do j = 1,it
            fsum = fsum + func(x, tau0, theta0, h)
            x = x + del
         end do
         s = 0.5d0 * (s + (b-a)*fsum/real(it,s2_dp)) 
                            ! This replaces s by its refined value.
      end if

    end subroutine trapzd_a


    !--------------------------------------------------------------------------
    ! qtrap_a
    !
    !! Computes the integral of the function func from a to b using the 
    !! trapezoid method.  Adapted from numerical recipes for the application 
    !! at hand.
    !! For use in computing integrals required to compute A(theta) and
    !! B(theta).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Returns the integral of the function func from a to b. 
    !!     The parameter EPS should be set to the desired fractional accuracy
    !!     and JMAX so that 2 to the power JMAX-1 is the maximum allowed
    !!     number of steps. Integration is performed by the trapezoidal rule.
    !!   - Adapted for use in evaluating integrals required to compute
    !!     A *or* B, (i.e. also takes parameters required to compute these
    !!    integrals).
    !!   - For our application tau0 and b will be the same, but clear and more
    !!     general to program as two separate variables.!!
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - theta0: theta0 position to evaluate integral refinement for.
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - integral: Value of the evaluated integral.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function qtrap_a(func, tau0, theta0, h, a, b) result(integral)

      real(s2_dp), intent(in) :: tau0, theta0, h
      real(s2_dp), intent(IN) :: a, b
      interface
         function func(tau, tau0, theta0, h) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: tau, tau0, theta0, h
           real(s2_dp) :: integrand
         end function func
      end interface
      real(s2_dp) :: integral

      integer, parameter :: JMAX=20
      real(s2_dp), parameter :: EPS=1d-3
      real(s2_dp) :: olds
      integer :: j

      olds = 0.0    !Initial value of olds is arbitrary.

      do j = 1,JMAX
         call trapzd_a(func, tau0, theta0, h, a,b,integral,j)
         if (j > 5) then     ! Avoid spurious early convergence.
            if (abs(integral-olds) < EPS*abs(olds) .or. &
                 (integral == 0.0 .and. olds == 0.0)) RETURN
         end if
         olds = integral
      end do

      ! If reach JMAX without finishing then call error.
      call bianchi_error(BIANCHI_ERROR_SKY_QUAD_STEP_EXCEED, 'qtrap_a')

    end function qtrap_a


    !--------------------------------------------------------------------------
    ! qsimp_a
    !
    !! Computes the integral of the function func from a to b using
    !! Simpson's rule.  Adapted from numerical recipes for the application 
    !! at hand.
    !! For use in computing integrals required to compute A(theta) and
    !! B(theta).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Returns the integral of the function func from a to b.
    !!     The parameter EPS should be set to the desired fractional 
    !!     accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
    !!     allowed number of steps. Integration is performed by Simpson's rule.
    !!   - Adapted for use in evaluating integrals required to compute
    !!     A *or* B, (i.e. also takes parameters required to compute these
    !!    integrals).
    !!   - For our application tau0 and b will be the same, but clear and more
    !!     general to program as two separate variables.
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - theta0: theta0 position to evaluate integral refinement for.
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - integral: Value of the evaluated integral.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function qsimp_a(func, tau0, theta0, h, a, b) result(integral)

      real(s2_dp), intent(in) :: tau0, theta0, h
      real(s2_dp), intent(IN) :: a, b
      interface
         function func(tau, tau0, theta0, h) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: tau, tau0, theta0, h
           real(s2_dp) :: integrand
         end function func
      end interface
      real(s2_dp) :: integral

      integer, parameter :: JMAX=20
      real(s2_dp), parameter :: EPS=1d-3
      integer :: j
      real(s2_dp) :: os,ost,st

      ost=0.0
      os= 0.0
      do j=1,JMAX
         call trapzd_a(func, tau0, theta0, h, a, b, st, j)
         integral = (4d0*st-ost)/3d0     !Compare equation (4.2.4).
         if (j > 5) then               !Avoid spurious early convergence.
            if (abs(integral-os) < EPS*abs(os) .or. &
                 (integral == 0.0 .and. os == 0.0)) RETURN
         end if
         os=integral
         ost=st
      end do

      ! If reach JMAX without finishing then call error.
      call bianchi_error(BIANCHI_ERROR_SKY_QUAD_STEP_EXCEED, 'qsimp_a')

    end function qsimp_a


    !--------------------------------------------------------------------------
    ! trapzd_ia
    !
    !! Computes nth stage of refinement of extended trapezoidal rule.
    !! Adapted from numerical recipes for the application at hand.
    !! For use in computing integrals IA_l *and* IB_l (although function
    !! name is only give the a subscript).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     This routine computes the nth stage of refinement of an extended 
    !!     trapezoidal rule. func is input as the name of the function to be 
    !!     integrated between limits a and b, also input. When called with
    !!     n=1, the routine returns as s the crudest estimate of 
    !!     int_b^a f(x)dx. Subsequent calls with n=2,3,... (in that sequential
    !!     order) will improve the accuracy of s by adding 2n-2 additional
    !!     interior points. s should not be modified between sequential calls.
    !!   - Adapted for use in evaluating integrals IA_l *or* IB_l, (i.e. also
    !!     takes parameters required to compute these integrals).
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute integral for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integral.
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - s: Value of integral computed to current order.
    !!   - n: Refinement order.
    !!   - [N_quad]:  Number of terms in use in direct quadrature.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine trapzd_ia(func, tauE, tau0, h, zE, c1, c3, l, &
        quad_AB, a, b, s, n, N_quad)

      real(s2_dp), intent(in) :: tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB
      real(s2_dp), intent(IN) :: a, b
      real(s2_dp), intent(INOUT) :: s
      integer, intent(IN) :: n
      integer, intent(in), optional :: N_quad
      interface
         function func(theta, tauE, tau0, h, zE, c1, c3, l, &
              quad, N, AB_grid) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta, tauE, tau0, h
           real(s2_dp), intent(in) :: zE, c1, c3
           integer, intent(in) :: l, quad
           integer, optional, intent(in) :: N
           real(s2_dp), optional, intent(in) :: AB_grid(0:)
           real(s2_dp) :: integrand
         end function func
      end interface

      real(s2_dp) :: del, fsum, x
      integer :: it, j

      if (n == 1) then
         s = 0.5d0 * (b-a)*(func(a, tauE, tau0, h, zE, c1, c3, l, &
             quad_AB, N_quad) &
           + func(b, tauE, tau0, h, zE, c1, c3, l, quad_AB, N_quad))
      else
         it = 2**(n-2)
         del = (b-a) / real(it, s2_dp)  
                             ! This is the spacing of the points to be added.
         x = a + 0.5d0*del
         fsum = 0d0
         do j = 1,it
            fsum = fsum + func(x, tauE, tau0, h, zE, c1, c3, l, &
              quad_AB, N_quad)
            x = x + del
         end do
         s = 0.5d0 * (s + (b-a)*fsum/real(it,s2_dp)) 
                            ! This replaces s by its refined value.
      end if

    end subroutine trapzd_ia


    !--------------------------------------------------------------------------
    ! qtrap_ia
    !
    !! Computes the integral of the function func from a to b using the 
    !! trapezoid method.  Adapted from numerical recipes for the application 
    !! at hand.
    !! For use in computing integrals IA_l *and* IB_l (although function
    !! name is only give the a subscript).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Returns the integral of the function func from a to b. 
    !!     The parameter EPS should be set to the desired fractional accuracy
    !!     and JMAX so that 2 to the power JMAX-1 is the maximum allowed
    !!     number of steps. Integration is performed by the trapezoidal rule.
    !!   - Adapted for use in evaluating integrals IA_l *or* IB_l, (i.e. also
    !!     takes parameters required to compute these integrals).
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute integral for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integral.
    !!   - [N_quad]:  Number of terms in use in direct quadrature.
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - integral: Value of the evaluated integral.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function qtrap_ia(func, tauE, tau0, h, zE, c1, c3, l, &
        quad_AB, a, b, N_quad) result(integral)

      real(s2_dp), intent(in) :: tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB
      real(s2_dp), intent(IN) :: a, b
      integer, intent(in), optional :: N_quad
      interface
         function func(theta, tauE, tau0, h, zE, c1, c3, l, &
            quad, N, AB_grid) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
           integer, intent(in) :: l, quad
           integer, optional, intent(in) :: N
           real(s2_dp), optional, intent(in) :: AB_grid(0:)
           real(s2_dp) :: integrand
         end function func
      end interface
      real(s2_dp) :: integral

      integer, parameter :: JMAX=20
      real(s2_dp), parameter :: EPS=1d-3
      real(s2_dp) :: olds
      integer :: j

      olds = 0.0    !Initial value of olds is arbitrary.

      do j = 1,JMAX
         call trapzd_ia(func, tauE, tau0, h, zE, c1, c3, l, &
           quad_AB, a,b,integral,j, N_quad)
         if (j > 5) then     ! Avoid spurious early convergence.
            if (abs(integral-olds) < EPS*abs(olds) .or. &
                 (integral == 0.0 .and. olds == 0.0)) RETURN
         end if
         olds = integral
      end do

      ! If reach JMAX without finishing then call error.
      call bianchi_error(BIANCHI_ERROR_SKY_QUAD_STEP_EXCEED, 'qtrap_ia')

    end function qtrap_ia


    !--------------------------------------------------------------------------
    ! qsimp_ia
    !
    !! Computes the integral of the function func from a to b using 
    !! Simpson's rule.  Adapted from numerical recipes for the application 
    !! at hand.
    !! For use in computing integrals IA_l *and* IB_l (although function
    !! name is only give the a subscript).
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Returns the integral of the function func from a to b.
    !!     The parameter EPS should be set to the desired fractional 
    !!     accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
    !!     allowed number of steps. Integration is performed by Simpson's rule.
    !!   - Adapted for use in evaluating integrals required to compute
    !!     A *or* B, (i.e. also takes parameters required to compute these
    !!    integrals).
    !!
    !! Variables:
    !!   - func: "Pointer" to integrand function.
    !!   - tauE: tauE parameter (see bianchi data type for explanation).
    !!   - tau0: tau0 parameter (see bianchi data type for explanation).
    !!   - h: h parameter (see bianchi data type for explanation).
    !!   - zE: zE parameter (see bianchi data type for explanation).
    !!   - c1: Bianchi simulation c1 parameter (constant parameter).
    !!   - c3: Bianchi simulation c3 parameter (constant parameter).
    !!   - l: Harmonic l to compute integral for.
    !!   - quad: Integer specifying quadrature rule to use when
    !!     evaluating integral.
    !!   - [N_quad]:  Number of terms in use in direct quadrature.
    !!   - a: Lower limit to evalue definite integral for.
    !!   - b: Upper limit to evalue definite integral for.
    !!   - integral: Value of the evaluated integral.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function qsimp_ia(func, tauE, tau0, h, zE, c1, c3, l, &
        quad_AB, a, b, N_quad) result(integral)

      real(s2_dp), intent(in) :: tauE, tau0, h, zE, c1, c3
      integer, intent(in) :: l, quad_AB
      real(s2_dp), intent(IN) :: a, b
      integer, intent(in), optional :: N_quad
      interface
         function func(theta, tauE, tau0, h, zE, c1, c3, l, &
            quad, N, AB_grid) result(integrand)
           use s2_types_mod
           real(s2_dp), intent(in) :: theta, tauE, tau0, h, zE, c1, c3
           integer, intent(in) :: l, quad
           integer, optional, intent(in) :: N
           real(s2_dp), optional, intent(in) :: AB_grid(0:)
           real(s2_dp) :: integrand
         end function func
      end interface
      real(s2_dp) :: integral

      integer, parameter :: JMAX=20
      real(s2_dp), parameter :: EPS=1d-3
      integer :: j
      real(s2_dp) :: os,ost,st

      ost=0.0
      os= 0.0
      do j=1,JMAX
         call trapzd_ia(func,tauE, tau0, h, zE, c1, c3, l, &
          quad_AB, a,b,st,j,N_quad)
         integral=(4d0*st-ost)/3d0     !Compare equation (4.2.4).
         if (j > 5) then               !Avoid spurious early convergence.
            if (abs(integral-os) < EPS*abs(os) .or. &
                 (integral == 0.0 .and. os == 0.0)) RETURN
         end if
         os=integral
         ost=st
      end do

      ! If reach JMAX without finishing then call error.
      call bianchi_error(BIANCHI_ERROR_SKY_QUAD_STEP_EXCEED, 'qsimp_ia')

    end function qsimp_ia


    !--------------------------------------------------------------------------
    ! Special function routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! plgndr
    !
    !! Computes the associated Legendre function for l and m.
    !!  Adapted from numerical recipes 
    !!
    !! Notes:
    !!   - Numerical recipies comment:
    !!     Computes the associated Legendre polynomial P_m^l(x).
    !!
    !! Variables:
    !!   - l: Legendre function l parameter.
    !!   - m: Legendre function m parameter.
    !!   - x: Point to evaluate specified Legendre funtion at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function plgndr(l,m,x)

      integer :: l, m
      real(s2_dp) :: plgndr, x

      integer :: i, ll
      real(s2_dp) ::  fact, pll, pmm, pmmp1, somx2

      if(m<0 .or. m>l .or. abs(x)>1d0) stop 'bad arguments in plgndr'

      pmm=1.                     ! Compute Pmm.
      if(m.gt.0) then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
         enddo
      endif
      if(l.eq.m) then
         plgndr=pmm
      else
         pmmp1=x*(2*m+1)*pmm      ! Compute Pm m+1.
         if(l.eq.m+1) then
            plgndr=pmmp1
         else                    ! Compute Pm  l , l > m+ 1.
            do ll=m+2,l
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            enddo
            plgndr=pll
         endif
      endif
      return
    end function plgndr


    !--------------------------------------------------------------------------
    ! asinh
    !
    !! Implementation of asinh, since no intrinsic Fortran function.
    !!
    !! Variables:
    !!  - x: Asinh argument.
    !!  - y: Asinh result.
    !
    !! @author J. D. McEwen
    !! @version 0.1 June 2005
    !
    ! Revisions:
    !   June 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function asinh(x) result(y)

      real(s2_dp), intent(in) :: x
      real(s2_dp) :: y

      y = log(x + sqrt(x**2.0d0 + 1.0d0))

    end function asinh


end module bianchi_sky_mod
