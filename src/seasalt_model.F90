!===============================================================================
! Seasalt for Bulk Aerosol Model
!===============================================================================
module seasalt_model

! TODO: clean up and move stuff to the object
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,       only: pcols, pver
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use aerosol_properties_mod, only: aerosol_properties
  use fanci_aerosol_properties_mod, only: fanci_aerosol_properties

  implicit none
  private

  ! public variables
  public :: seasalt_specprop_ndx

  ! public procedures
  public :: seasalt_init
  public :: seasalt_emis

  ! initialize public

  ! module variables
  integer :: seasalt_specprop_ndx = 0

  contains

   !=============================================================================
   !=============================================================================
   subroutine seasalt_init(aero_props)
     use cam_history,   only: addfld, fieldname_len
     use constituents,  only: cnst_get_ind
     use string_utils,  only: int2str

     type(fanci_aerosol_properties), intent(in) :: aero_props

     ! local variables
     integer              :: ispecprop, ibin, irange
     integer              :: specprop_ndx
     character(len=10)    :: type

     character(len=fieldname_len) :: dummy
     character(len=*), parameter :: subname = 'seasalt_init'

    ! Initalize seasalt vars

     do ispecprop = 1, aero_props%nspecies_tot()
        type = aero_props%spectype(ispecprop)

        if (trim(type) == 'seasalt') then
            seasalt_specprop_ndx = ispecprop
        end if
    end do

   end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed, u10in, srf_temp, ocnfrc, ncol, cflx, aero_props )
! TODO: https://acp.copernicus.org/articles/11/4587/2011/acp-11-4587-2011.pdf

    ! dummy arguments
    real(r8), intent(in)    :: u10cubed(:)
    real(r8), intent(in)    :: u10in(:)
    real(r8), intent(in)    :: srf_temp(:)
    real(r8), intent(in)    :: ocnfrc(:)
    integer,  intent(in)    :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    type(fanci_aerosol_properties), intent(in) :: aero_props

    !--------CMS (Clarke, Monahan, and Smith source function) from CARMA: carma_model_mod.F90 /src/physics/carma/models/sea_salt/carma_model_mod.F90

    ! ------------------------------------------------------------
    ! ----  Clarke Source Function. Coefficients for Ai    -------
    ! ------------------------------------------------------------
    real(r8), parameter :: beta01 =-5.001e3_r8
    real(r8), parameter :: beta11 = 0.808e6_r8
    real(r8), parameter :: beta21 =-1.980e7_r8
    real(r8), parameter :: beta31 = 2.188e8_r8
    real(r8), parameter :: beta41 =-1.144e9_r8
    real(r8), parameter :: beta51 = 2.290e9_r8
    real(r8), parameter :: beta02 = 3.854e3_r8
    real(r8), parameter :: beta12 = 1.168e4_r8
    real(r8), parameter :: beta22 =-6.572e4_r8
    real(r8), parameter :: beta32 = 1.003e5_r8
    real(r8), parameter :: beta42 =-6.407e4_r8
    real(r8), parameter :: beta52 = 1.493e4_r8
    real(r8), parameter :: beta03 = 4.498e2_r8
    real(r8), parameter :: beta13 = 0.839e3_r8
    real(r8), parameter :: beta23 =-5.394e2_r8
    real(r8), parameter :: beta33 = 1.218e2_r8
    real(r8), parameter :: beta43 =-1.213e1_r8
    real(r8), parameter :: beta53 = 4.514e-1_r8
    real(r8)            :: A1                                 ! Coefficient Ak in Clarkes's source function
    real(r8)            :: A2
    real(r8)            :: A3

    ! --------------------------------------------------------------------
    ! ---- constants in calculating the particle wet radius [Gerber, 1985]
    ! --------------------------------------------------------------------
    real(r8), parameter :: c1   = 0.7674_r8        ! .
    real(r8), parameter :: c2   = 3.079_r8         ! .
    real(r8), parameter :: c3   = 2.573e-11_r8     ! .
    real(r8), parameter :: c4   = -1.424_r8        ! constants in calculating the particle wet radius

    ! ---------------------------------------------
    ! coefficient A1, A2 in Andreas's Source funcion
    ! ---------------------------------------------
    real(r8)            ::A1A92
    real(r8)            ::A2A92                         ! coefficient as desired by Andreas source function

    ! ---------------------------------------------
    ! Parameter in Monahan
    ! ---------------------------------------------
    real(r8)            :: B_mona

    ! ---------------------------------------------
    ! coefficient in Smith's Source funcion
    ! ---------------------------------------------
    real(r8), parameter ::  f1 = 3.1_r8
    real(r8), parameter ::  f2 = 3.3_r8
    real(r8), parameter ::  r1 = 2.1_r8
    real(r8), parameter ::  r2 = 9.2_r8

    ! Fluxes for Clarke, Monahan, Smith
    real(r8)            :: Monahan, Clarke, Smith             ! dF/dr [#/m2/s/um]
    real(r8)            :: ncflx                              ! dF/dr [#/m2/s/um]

    ! Local
    real(r8)            :: wcap                               ! whitecap coverage
    real(r8)            :: rpdry, rpdry_cm                    ! dry radius
    real(r8)            :: r80, r80cm
    real(r8), allocatable :: bin_centers(:)
    real(r8)            :: bin_width_um
    real(r8), parameter :: xkar = 0.4_r8                      ! Von Karman constant
    real(r8)            :: u14, cd_smith         ! 14m wind velocity, friction velocity and drag
    real(r8), allocatable :: cflx_tmp(:,:)
    integer             :: sslt_bin
    ! running variables
    integer             :: ibin, icol, irange

    character(len=*), parameter :: subname = 'seasalt_emis'
    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------


    if (.not. aero_props%is_active('seasalt')) return

! -> get q-array indices for emissions from aero_props%sped_bin_q_ndx(ispec, spec_bin_ndx) and spec_mmr_q_ndx(ispec, spec_range_ndx)

    allocate(bin_centers(aero_props%nbins()))
    allocate(cflx_tmp(pcols, aero_props%spec_nbin(spectype='seasalt')))
    cflx_tmp    = 0._r8

    bin_centers = aero_props%bin_centers(aero_props%nbins()) * 1.e6_r8 ! to um

    sslt_bin = 0
    do ibin = aero_props%spec_bin_ndx(seasalt_specprop_ndx, 1), aero_props%spec_bin_ndx(seasalt_specprop_ndx, aero_props%spec_nbin(spectype='seasalt'))
        sslt_bin = sslt_bin + 1
        rpdry = bin_centers(ibin)               ! [um]
        rpdry_cm = rpdry * 1.e-4_r8             ! [cm]
! todo: get 80% RH radius differently?
        r80cm   = (c1 *  (rpdry_cm) ** c2 / (c3 * rpdry_cm ** c4 - log10(0.8_r8)) + (rpdry_cm)**3._r8) ** (1._r8/3._r8) ! [cm]
        r80     = r80cm *1.e4_r8                ! [um]

        bin_width_um = (aero_props%bin_bounds(ibin, 2) - aero_props%bin_bounds(ibin, 1)) * 1.e6_r8 ! bin width in [um]

        A1 = beta01 + beta11*(2._r8*rpdry) + beta21*(2._r8*rpdry)**2 + &
                beta31*(2._r8*rpdry)**3 + beta41*(2._r8*rpdry)**4 + beta51*(2._r8*rpdry)**5
        A2 = beta02 + beta12*(2._r8*rpdry) + beta22*(2._r8*rpdry)**2 + &
                 beta32*(2._r8*rpdry)**3 + beta42*(2._r8*rpdry)**4 + beta52*(2._r8*rpdry)**5
        A3 = beta03 + beta13*(2._r8*rpdry) + beta23*(2._r8*rpdry)**2 + &
                 beta33*(2._r8*rpdry)**3 + beta43*(2._r8*rpdry)**4 + beta53*(2._r8*rpdry)**5


        do icol = 1, ncol

        ! Initialize surface fluxes
            ncflx       = 0.0_r8
            Monahan     = 0.0_r8
            Clarke      = 0.0_r8
            Smith       = 0.0_r8

            ! white cap fraction
            wcap = 3.84e-6_r8 * u10cubed(icol)      ! in percent, ie., 75%, wcap = 0.75

            !****************************************
            !        Hoppel correction factor
            !        Smith drag coefficients and etc
            !****************************************
            if (u10in(icol) .le. 10._r8) then
                cd_smith = 1.14e-3_r8
            else
                cd_smith = (0.49_r8 + 0.065_r8 * u10in(icol)) * 1.e-3_r8
            end if

! TODO: replaced rdry with rpdry here
            !Clarke
            if (rpdry .lt. 0.066_r8) then
              Clarke = A1 * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rpdry .ge. 0.066_r8) .and. (rpdry .lt. 0.6_r8)) then
              Clarke = A2 * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rpdry .ge. 0.6_r8) .and. (rpdry .lt. 4.0_r8)) then
              Clarke = A3 * 1.e4_r8 * wcap                      ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)                 ! dF/dr    [#/s/m2/um]
            end if

            !Monahan
            B_Mona = (0.38_r8 - log10(r80)) / 0.65_r8
            Monahan = 1.373_r8 * u10cubed(icol) * r80 ** (-3._r8) * &
                 (1._r8 + 0.057_r8 * r80**1.05_r8) * 10._r8 ** (1.19_r8 * exp(- B_Mona **2))

            !Smith
            u14 = u10in(icol) * (1._r8 + cd_smith**0.5_r8 / xkar*log(14._r8 / 10._r8))  ! 14 meter wind
            A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
            A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
            Smith = A1A92*exp(-f1 *(log(r80 / r1))**2) + A2A92*exp(-f2 * (log(r80 / r2))**2)     ! dF/dr   [#/m2/s/um]

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            !     CMS1 or CMS2
            !%%%%%%%%%%%%%%%%%%%%%%%%%
  !          if (rdry .lt. 0.1_r8) then   ! originally cut at 0.1 um
            ! ***CMS1*****
            if (rpdry .lt. 1._r8) then    ! cut at 1.0 um
            ! ***CMS2*****
  !          if (rdry .lt. 2._r8) then    ! cut at 2.0 um
                ncflx = Clarke
            else
                if (u10in(icol) .lt. 9._r8) then
                    ncflx = Monahan
                else
                    if (Monahan .gt. Smith) then
                        ncflx = Monahan
                    else
                        ncflx = Smith
                    end if
                end if
            end if

            ! Apply Hoppel correction
!            ncflx = ncflx * fref ! fref is = 1._r8 in CARMA

          ! convert ncflx [#/m^2/s/um] to surfaceFlx [kg/m^2/s]

            ! Number flux
            cflx_tmp(icol, sslt_bin) = ncflx * bin_width_um * ocnfrc(icol)  ! [#/m2/s]

            cflx(icol, aero_props%spec_bin_q_ndx(seasalt_specprop_ndx, sslt_bin)) = cflx(icol, aero_props%spec_bin_q_ndx(seasalt_specprop_ndx, sslt_bin)) &
                                                                                + cflx_tmp(icol, sslt_bin)
            ! mass flux

            do irange = 1, aero_props%spec_nrange(seasalt_specprop_ndx)

            if (aero_props%bins2ranges(aero_props%spec_bin_ndx(seasalt_specprop_ndx, sslt_bin)) == aero_props%spec_range_ndx(seasalt_specprop_ndx, irange)) then

                cflx(icol, aero_props%spec_mmr_q_ndx(seasalt_specprop_ndx, irange)) = &
                                                                                    cflx(icol, aero_props%spec_mmr_q_ndx(seasalt_specprop_ndx, irange)) &
                                                                                  + cflx_tmp(icol,sslt_bin) &
                                                                                  * aero_props%density(seasalt_specprop_ndx) &
                                                                                  * aero_props%particle_volume(aero_props%spec_bin_ndx(seasalt_specprop_ndx, sslt_bin))
                end if
            end do
        end do

    end do


  end subroutine seasalt_emis

end module seasalt_model
