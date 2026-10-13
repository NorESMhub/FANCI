!===============================================================================
! FANCI Sectional Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,      only: r8 => shr_kind_r8
  use constituents,      only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,            only: pcols, pver, pverp
  use cam_abortutils,    only: endrun
  use cam_logfile,       only: iulog
  use perf_mod,          only: t_startf, t_stopf
  use camsrfexch,        only: cam_in_t, cam_out_t
  use aerodep_flx,       only: aerodep_flx_prescribed
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,    only: physics_buffer_desc
  use physconst,         only: gravit, rair
  use dust_model,        only: dust_specprop_ndx
  use seasalt_model,     only: seasalt_specprop_ndx
  use spmd_utils,        only: masterproc
  use physics_buffer,    only: pbuf_get_field, pbuf_get_index, pbuf_get_chunk
  use cam_history,       only: outfld
  use infnan,            only: nan, assignment(=)
  use fanci_aerosol_properties_mod, only: fanci_aerosol_properties
  use fanci_aerosol_state_mod, only: fanci_aerosol_state, aero_state_ptr
  use string_utils,      only: int2str

  implicit none
  private

  public :: aero_model_get_state
  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea    ! tropospheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea   ! stub

  ! name of the aerosol scheme
  public :: aero_modelname
  character(len=*), parameter :: aero_modelname = 'fanci'

 ! Misc private data

  !integer :: so4_ndx, cb2_ndx, oc2_ndx, nit_ndx
  !integer :: soa_ndx, soai_ndx, soam_ndx, soab_ndx, soat_ndx, soax_ndx

  ! aerosol_nl Namelist variables
!  character(len=16), allocatable :: wetdep_list(:)
!  character(len=16), allocatable :: drydep_list(:)

  !integer :: ndrydep = 0
  integer :: nwetdep = 0
  logical :: drydep_lq(pcnst) = .false.
  logical :: wetdep_lq(pcnst) = .false.

  real(r8) :: aer_sol_facti(pcnst) ! in-cloud solubility factor
  real(r8) :: aer_sol_factb(pcnst) ! below-cloud solubility factor
  real(r8) :: aer_scav_coef(pcnst)

  integer :: fracis_idx = 0
  integer :: prain_idx  = 0

  type(fanci_aerosol_properties), pointer :: aero_props=>null()
  !type(fanci_aerosol_state), pointer :: aero_state=>null()

  type(aero_state_ptr), allocatable :: master_aero_state(:)

contains

  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)
    use dust_model,      only: dust_readnl

    use oslo_aero_control, only: oslo_aero_ctl_readnl ! TODO: use the file in oslo_aero directly? currently this is a copy in the local chemistry.F90
    use aero_wetdep_cam, only: aero_wetdep_readnl
    ! filepath for file containing namelist input
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer                     :: unitn, ierr, ind, pos
    character(len=50)           :: tmp

    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! modal_aero: read aerosol_nl: aer_drydep_list, modal_strat_sulfate, modal_accum_coarse_exch, seasalt_emis_scale
    ! oslo aero: read aerosol_nl: sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    call oslo_aero_ctl_readnl(nlfile)

    call aero_wetdep_readnl(nlfile)

    call dust_readnl(nlfile)

    if (.not. aerodep_flx_prescribed()) then
        aero_props => fanci_aerosol_properties(nlfile) ! calls constructor function in fanci_aerosol_properties
    end if

    ! initialize props
    ! nlfile input optional
    ! TODO: initialize gasaerexch

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()
    character(len=*), parameter :: subname = 'aero_model_register'

    ! modal_aero: modal_aero_data_reg: allocate all kind of stuff, etc with number modes
    ! oslo aero: aero_register: lots of cnst_get_ind calls for all tracers, set aerosol types

    if (masterproc) then
        write(iulog,*) subname//' nothing to do here..'
    end if

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d, phys_state)

    use mo_chem_utls,   only: get_inv_ndx, get_spc_ndx
    use chem_mods,   only : gas_pcnst
    use mo_tracname,              only: solsym
    use cam_history,    only: addfld, add_default, horiz_only
    use phys_control,   only: phys_getopts
    use dust_model,     only: dust_init
    use seasalt_model,  only: seasalt_init
    use string_utils,   only: int2str
    !use mo_setsox,      only : setsox, has_sox
  use ppgrid,               only: begchunk, endchunk, pcols, pver
    use aero_deposition_cam, only: aero_deposition_cam_init
    use aer_drydep_mod,  only: inidrydep
    use aero_wetdep_cam, only: aero_wetdep_init

    use oslo_aero_ocean, only: oslo_aero_ocean_init ! TODO: DMS, add to build-namelist and chemistry.F90 and as well

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state),    intent(in)    :: phys_state(begchunk:endchunk)     ! Physics state variables

    ! local vars
    integer           :: m, id, ierr, ibin, ispec, ichunk, lchnk, icnst
    integer           :: specprop_ndx, q_ndx, spec_bin_ndx
    integer           :: ind, irange
    logical           :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical           :: history_dust    ! Output dust
    logical           :: history_chemistry ! Output Chemistry
    character(len=2)  :: unit_basename ! Units 'kg' or '1'
    character(len=10) :: aerosol_names(500), aerosol_cw_names(500)
    character(len=20) :: dummy
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    character(len=12), parameter :: subname = 'aero_model_init'

    fracis_idx      = pbuf_get_index('FRACIS')
    prain_idx       = pbuf_get_index('PRAIN')

    call oslo_aero_ocean_init() ! DMS

    ! TODO: fix this :)
    call phys_getopts( history_aerosol_out   = history_aerosol, &
                       history_dust_out      = history_dust,    &
                       history_chemistry_out = history_chemistry   )

    if (.not. aerodep_flx_prescribed()) then
        aero_props => fanci_aerosol_properties() ! calls constructor function in fanci_aerosol_properties
        allocate(master_aero_state(begchunk:endchunk))
        do lchnk = begchunk, endchunk
            phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
            master_aero_state(lchnk)%ptr => fanci_aerosol_state(phys_state(lchnk), phys_buffer_chunk)
        end do
  !     aero_state => fanci_aerosol_state(phys_state, pbuf)
      !  end do
        call aero_deposition_cam_init(aero_props)
    end if

    call dust_init(aero_props)
    call seasalt_init(aero_props)

! TODO: if drydep active:
    call inidrydep(rair, gravit)
    do specprop_ndx = 1, aero_props%nspecies_tot()
        do spec_bin_ndx = 1,aero_props%spec_nbin(specprop_ndx)
            q_ndx = aero_props%spec_bin_q_ndx(specprop_ndx, spec_bin_ndx)
!            drydep_lq(q_ndx) = .true.
!            wetdep_lq(q_ndx) = .true.
        end do
    end do

    dummy = 'RAM1'
    call addfld (dummy,horiz_only, 'A','frac','RAM1')
    if ( history_aerosol ) then
        call add_default (dummy, 1, ' ')
    endif
    dummy = 'airFV'
    call addfld (dummy,horiz_only, 'A','frac','FV')
    if ( history_aerosol ) then
        call add_default (dummy, 1, ' ')
    endif

    if (aero_props%is_active('dust')) then

        do irange = 1, aero_props%spec_nrange(dust_specprop_ndx)
            dummy = trim(aero_props%spec_tracernames(dust_specprop_ndx, irange)) // 'SF'
            call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dummy)//' dust surface emission')
            if (history_aerosol) then
                call add_default (dummy, 1, ' ')
            endif
        end do

       dummy = 'DSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol .or. history_dust) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,horiz_only, 'A','frac','Soil erodibility factor')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif
    endif

    if (aero_props%is_active('seasalt')) then

        do irange = 1, aero_props%spec_nrange(seasalt_specprop_ndx)
            dummy = trim(aero_props%spec_tracernames(seasalt_specprop_ndx, irange)) // 'SF'
            call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dummy)//' seasalt surface emission')
            if (history_aerosol) then
                call add_default (dummy, 1, ' ')
            endif
        end do

        dummy = 'SSTSFMBL'
        call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
        if (history_aerosol) then
           call add_default (dummy, 1, ' ')
        endif
    endif

    aerosol_names = ''
    aerosol_cw_names = ''
    ind = 0
    do ispec = 1, aero_props%nspecies_tot()
        do irange = 1, aero_props%spec_nrange(ispec)
            ind = ind+1
            aerosol_names(ind) = aero_props%spec_tracernames(ispec, irange)
            aerosol_cw_names(ind) = trim(aero_props%spec_tracernames(ispec, irange))//'_cw'
        end do
    end do

    do ibin = 1, aero_props%nbins()
        aerosol_names(ind + ibin) = "num_"//int2str(ibin)
        aerosol_cw_names(ind + ibin) = "num_"//trim(int2str(ibin))//'_cw'
    end do

    do m = 1,ibin+ind-1

       ! units
       if (aerosol_names(m)(1:3) == 'num') then
          unit_basename = '#'
       else
          unit_basename = 'kg'
       endif

       call addfld (trim(aerosol_names(m))//'DDF', horiz_only,  'A',unit_basename//'/m2/s ', &
            trim(aerosol_names(m))//' dry deposition flux at bottom (grav + turb)')
       call addfld (trim(aerosol_names(m))//'TBF', horiz_only,  'A',unit_basename//'/m2/s',  &
            trim(aerosol_names(m))//' turbulent dry deposition flux')
       call addfld (trim(aerosol_names(m))//'GVF', horiz_only,  'A',unit_basename//'/m2/s ', &
            trim(aerosol_names(m))//' gravitational dry deposition flux')
       call addfld (trim(aerosol_names(m))//'DTQ', (/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(aerosol_names(m))//' dry deposition')
       call addfld (trim(aerosol_names(m))//'DDV', (/ 'lev' /), 'A','m/s',                   &
            trim(aerosol_names(m))//' deposition velocity')

       call addfld (trim(aerosol_cw_names(m)), (/ 'lev' /), 'A', unit_basename//'/kg ', &
           trim(aerosol_cw_names(m))//' in cloud water')
       call addfld (trim(aerosol_cw_names(m))//'DDF', horiz_only,  'A',unit_basename//'/m2/s ', &
           trim(aerosol_cw_names(m))//' dry deposition flux at bottom (grav + turb)')
       call addfld (trim(aerosol_cw_names(m))//'TBF', horiz_only,  'A',unit_basename//'/m2/s',  &
           trim(aerosol_cw_names(m))//' turbulent dry deposition flux')
       call addfld (trim(aerosol_cw_names(m))//'GVF', horiz_only,  'A',unit_basename//'/m2/s ', &
           trim(aerosol_cw_names(m))//' gravitational dry deposition flux')

       if ( history_aerosol.or.history_chemistry ) then
          call add_default (trim(aerosol_names(m))//'DDF', 1, ' ')
          call add_default (trim(aerosol_cw_names(m)), 1, ' ')
       endif
       if ( history_aerosol ) then
          call add_default (trim(aerosol_names(m))//'TBF', 1, ' ')
          call add_default (trim(aerosol_names(m))//'GVF', 1, ' ')
          call add_default (trim(aerosol_cw_names(m))//'TBF', 1, ' ')
          call add_default (trim(aerosol_cw_names(m))//'GVF', 1, ' ')
          call add_default (trim(aerosol_cw_names(m))//'DDF', 1, ' ')
       endif
    enddo

    call aero_wetdep_init()

    if(aero_props%is_active('sulfate')) then
        do icnst = 1, gas_pcnst
            call addfld('GS_'//trim(solsym(icnst)),horiz_only, 'A','kg/m2/s', &
                trim(solsym(icnst))//' gas chemistry tendency')
            if (history_aerosol .or. history_chemistry) then
                call add_default ('GS_'//trim(solsym(icnst)), 1, ' ')
            endif
        end do

        ! TODO: if active(sulfate or nitrate) call gasaerexch_init()
        return
    end if
    end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
function aero_model_get_state(lchnk) result(state_ptr)
    use aerosol_state_mod, only: aerosol_state
    integer, intent(in) :: lchnk
    class(aerosol_state), pointer :: state_ptr
    state_ptr => master_aero_state(lchnk)%ptr
end function aero_model_get_state

  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use dust_sediment_mod, only: dust_sediment_tend
    use aer_drydep_mod,    only: d3ddflux, calcram
! TODO: add cw stuff
! TODO: use ndrydep/nwetdep?

!    use modal_aero_data,   only: nspec_amode -> species in each mode
!    use modal_aero_data,   only: numptr_amode -> r-array index for nr mixing ratio for aerosol
!    use modal_aero_data,   only: lmassptr_amode -> r-array index for the mixing ratio
!       (moles-x/mole-air) for chemical species l in aerosol mode m
!       that is in clear air or interstitial air (but not in cloud water)
    use aero_deposition_cam,only: aero_deposition_cam_setdry

! TODO: move out state, cam_in, cam_out, ptend, pbuf -> e.g. move to different subroutine for now
    ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

     ! local decarations

    integer :: lchnk                   ! chunk identifier
    integer :: ncol                    ! number of atmospheric columns
    integer :: jvlc                    ! index for last dimension of vlc_xxx arrays
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: lspec                   ! index for aerosol number / chem-mass / water-mass
    integer :: m                       ! aerosol mode index
    integer :: mm                      ! tracer index
    integer :: i
    integer :: ibin, nbins, icol, ilev, irange, ispec, ierr, bin_ndx

    real(r8) :: sflx(pcols)
    real(r8) :: sflx_num(pcols)
    real(r8) :: sflx_range_species(pcols)
    real(r8),allocatable :: sflx_range(:,:)

    real(r8) :: tvs(pcols,pver)
    real(r8) :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp) ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity ! TODO: get rid of last dimension?
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_trb(pcols,4)          ! dep velocity
    real(r8) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8) :: massfrac(pcols, pver)
    real(r8), allocatable :: bin_centers(:)

    real(r8) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
!    real(r8), pointer :: fldcw(:,:)
!    real(r8), pointer :: dgncur_awet(:,:,:)
!    real(r8), pointer :: wetdens(:,:,:)
!    real(r8), pointer :: qaerwat(:,:,:)

    real(r8) :: bin_mmr_tend(pcols, pver)
    real(r8) :: bin_mmr_tot(pcols, pver)
    real(r8), pointer :: bin_num(:,:)                ! number concentration in #/kg
    real(r8), pointer :: range_mmr(:,:)              ! mass mixing ratio in kg/kg
    real(r8), allocatable :: bin_num_tend(:,:,:)     ! pcols, pver, nbins
    real(r8), allocatable :: range_mmr_tend(:,:,:) ! pcols, pver, nranges
    character(len=15) :: species_tracername

    character(len=*), parameter :: subname = 'aero_model_drydep'

    nbins = aero_props%nbins()

    allocate(bin_num_tend(pcols, pver, aero_props%nbins()), stat=ierr)
    if( ierr /= 0 ) then
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating bin_num_tend")
    end if
    allocate(range_mmr_tend(pcols, pver, aero_props%nranges()), stat=ierr)
    if( ierr /= 0 ) then
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating range_mmr_tend")
    end if
    allocate(bin_centers(nbins), stat=ierr)
    if( ierr /= 0 ) then
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating bin_centers")
    end if
    allocate(sflx_range(pcols, aero_props%nranges()), stat=ierr)
    if( ierr /= 0 ) then
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating sflx_range")
    end if

    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    aerdepdryis = 0._r8
    aerdepdrycw = 0._r8
    bin_num_tend = 0._r8
    range_mmr_tend = 0._r8
    sflx_range = 0._r8
    sflx = 0._r8
    sflx_num = 0._r8
    sflx_range_species = 0._r8

! TODO MAKE AERDEPDRYIS and AERDEPDRYCW

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
                  state%pdel(:,pver),fvin,fv)

    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )

    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', lq=drydep_lq)

    tvs(:ncol,:) = state%t(:ncol,:)
    rho(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

! TODO: get deposition velocities for cloud stuff
!    rad_drop(:,:) = 5.0e-6_r8
!    dens_drop(:,:) = rhoh2o
!    sg_drop(:,:) = 1.46_r8
!
    dens_aer(:,:) = 0._r8
    bin_centers = aero_props%bin_centers(nbins)

    do irange = 1, aero_props%nranges()
        call master_aero_state(lchnk)%ptr%update_range(irange, ncol)
    end do

    irange = 1
    do ibin = 1, nbins  ! main loop over aerosol size bins aero
        irange = aero_props%bins2ranges(ibin)

        do lphase = 1, 2 ! interstitial/cloud borne forms
            if (lphase == 1) then ! interstitial
                ! reset tmp arrays
                bin_mmr_tend = 0._r8
                bin_mmr_tot = 0._r8
                bin_num => null()

! TODO: use WET radius and density in future!!
                rad_aer(1:ncol,:) = bin_centers(ibin)
                dens_aer(1:ncol,:) = master_aero_state(lchnk)%ptr%bin_dry_density(ibin, ncol)
                jvlc = 1 ! TODO: remove since we don't use masses ?

                ! calculate deposition velocities of single particles
                call aero_depvel_part(ncol,state%t(:,:), state%pmid(:,:), ram1, fv, &
                             vlc_dry(:,:,jvlc), vlc_trb(:, jvlc), vlc_grv(:,:,jvlc), &
                             rad_aer(:,:), dens_aer(:,:), lchnk)
            ! if lphase == 2 then cloud-borne

    ! loop through species_in_bin
    ! do some weird jvlc stuff -> find "mm", tracer index => use ncnst_tot
    ! jvlc = 1 number dry
    ! jvlc = 2
    ! jvlc = 3
    ! jvlc = 4
                ! get total mass mixing ratio in a bin with #/kg and total density
                do icol = 1, ncol
                    do ilev = 1, pver
                        bin_mmr_tot(icol, ilev) = master_aero_state(lchnk)%ptr%ambient_total_bin_mmr(aero_props, ibin, icol, ilev)
                    end do
                end do
                ! convert velocity to Pa/s
                pvmzaer(:ncol,1)=0._r8
                pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

                ! calculate deposition fluxes NOTE: "dust_sediment_tend" is valid for all aerosol, not just dust
                ! state%q has been changed to bin_mmr_tot (intent(in))
                ! ptend%q has been changed to bin_mmr_tend(pcols, pver)
                call master_aero_state(lchnk)%ptr%get_ambient_num(ibin, bin_num)

                call dust_sediment_tend(ncol, dt, state%pint(:,:), state%pmid, state%pdel, state%t, bin_num(:,:), pvmzaer, bin_num_tend(:,:, ibin), sflx_num )
                call dust_sediment_tend(ncol, dt, state%pint(:,:), state%pmid, state%pdel, state%t, bin_mmr_tot(:,:), pvmzaer, bin_mmr_tend(:,:), sflx )
                ! calculate #/kg tendency and put tendency to state
                bin_num = bin_num + bin_num_tend(:ncol,:, ibin)

                dep_trb = 0._r8
                dep_grv = 0._r8

                do i=1, ncol
                    if ( vlc_dry(i,pver,jvlc) /= 0._r8 ) then
                        dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                        dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                    end if
                end do

                call outfld( 'num_'//trim(int2str(ibin))//'DDF', sflx_num, pcols, lchnk)
                call outfld( 'num_'//trim(int2str(ibin))//'TBF', dep_trb, pcols, lchnk)
                call outfld( 'num_'//trim(int2str(ibin))//'GVF', dep_grv, pcols, lchnk)
                call outfld( 'num_'//trim(int2str(ibin))//'DTQ', bin_num_tend(:ncol,:,ibin), pcols, lchnk)
                mm = aero_props%indexer(ibin, 0)
                aerdepdryis(:ncol, mm) = sflx(:ncol)
          !      write(6,*)"DEBUG: maxval for number: ", maxval(aerdepdryis(:ncol,mm))
            end if
        end do
        ! add up mass in a range
        range_mmr_tend(:ncol,:,irange) = range_mmr_tend(:ncol,:,irange) + bin_mmr_tend(:ncol,:)
        sflx_range(:,irange) = sflx_range(:,irange) + sflx
    end do
        ! calculate the tendency for each species/range
    do irange = 1, aero_props%nranges()
        do ispec = 1, aero_props%range_nspecies(irange)
            sflx_range_species = 0._r8
            species_tracername = ''
            range_mmr => null()

            ! mass fraction of each species
            massfrac(:ncol,:) = master_aero_state(lchnk)%ptr%aero_range_state(irange)%massfrac(:,:,ispec)

            ! move tendency into aero range state
            bin_ndx = aero_props%range_bounds(irange, 1)
            call master_aero_state(lchnk)%ptr%get_ambient_mmr_0list(ispec, bin_ndx, range_mmr)
            range_mmr(:ncol,:) = range_mmr(:ncol,:) + range_mmr_tend(:ncol,:,irange)*massfrac(:ncol, :)

            ! use mass fractions at lowest level to get surface flux TODO: sedimentation out of higher layers?
            do icol = 1, ncol
!                do ilev = 1, pver
                    if (massfrac(icol, pver) > 0._r8) then
                        sflx_range_species(icol) = sflx_range(icol,irange) * massfrac(icol, pver)
                    end if
!                end do
            end do

            species_tracername = master_aero_state(lchnk)%ptr%aero_range_state(irange)%range_name(ispec)

            call outfld( trim(species_tracername)//'DDF', sflx_range_species, pcols, lchnk)
            call outfld( trim(species_tracername)//'DTQ', range_mmr_tend(:ncol,:,irange)*massfrac(:ncol, :), pcols, lchnk)

            ibin = aero_props%range_bounds(irange, 1)
            mm = aero_props%indexer(ibin, ispec)
            ! skip the "filler spots" in the index array
            if ( mm > 0 ) then
                aerdepdryis(:ncol, mm) = sflx_range_species(:ncol)
            end if
        end do ! species loop
    end do ! range loop

! rebin bulk fluxes for 'dust'
    ! rebin_bulk_fluxes prep
    ! Mass of species in bin

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call aero_deposition_cam_setdry(aerdepdryis, aerdepdrycw, cam_out)
    endif

  end subroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)
      use aero_wetdep_cam, only: aero_wetdep_tend

      ! args

      type(physics_state), intent(in)    :: state       ! Physics state variables
      real(r8),            intent(in)    :: dt          ! time step
      real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
      type(cam_out_t),     intent(inout) :: cam_out     ! export state
      type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
      type(physics_buffer_desc), pointer :: pbuf(:)
                 integer  :: lchnk                    ! chunk identifier
      lchnk = state%lchnk

      ! smb: TODO REMOVE THIS.
      !if (nwetdep<1) return
      call aero_wetdep_tend(state, dt, dlf, cam_out, ptend, pbuf)

  end subroutine aero_model_wetdep
!           use wetdep,        only : wetdepa_v1, wetdep_inputs_set, wetdep_inputs_t
!
!           ! args
!
!           type(physics_state), intent(in)    :: state       ! Physics state variables
!           real(r8),            intent(in)    :: dt          ! time step
!           real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
!           type(cam_out_t),     intent(inout) :: cam_out     ! export state
!           type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
!           type(physics_buffer_desc), pointer :: pbuf(:)
!
!           ! local vars
!
!           integer  :: ncol                     ! number of atmospheric columns
!           integer  :: lchnk                    ! chunk identifier
!           integer  :: m,mm, i,k
!
!           real(r8) :: sflx_tot_dst(pcols)
!           real(r8) :: sflx_tot_slt(pcols)
!
!           real(r8) :: iscavt(pcols, pver)
!           real(r8) :: scavt(pcols, pver)
!           real(r8) :: scavcoef(pcols,pver)     ! Dana and Hales coefficient (/mm) (0.1)
!           real(r8) :: sflx(pcols)              ! deposition flux
!
!           real(r8) :: icscavt(pcols, pver)
!           real(r8) :: isscavt(pcols, pver)
!           real(r8) :: bcscavt(pcols, pver)
!           real(r8) :: bsscavt(pcols, pver)
!
!           real(r8) :: sol_factb, sol_facti
!
!           real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
!           real(r8) :: cldv(pcols,pver)         ! cloudy volume undergoing scavenging
!           real(r8) :: cldvcu(pcols,pver)       ! Convective precipitation area at the top interface of current layer
!           real(r8) :: cldvst(pcols,pver)       ! Stratiform precipitation area at the top interface of current layer
!
!           real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble
!
!           type(wetdep_inputs_t) :: dep_inputs  ! obj that contains inputs to wetdepa routine
!
!           character(len=*), parameter :: subname = 'aero_model_wetdep'
!
!           call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
!
!           call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)
!
!       if (nwetdep<1) return
!
!       call endrun(subname//":: is not yet implemented")
!
!         endsubroutine aero_model_wetdep
!
  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  state, mmr, radmean, relhum, pmid, temp, strato_sad, sulfate,  m, ltrop, &
                  dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total, reff_trop )

    use mo_constants, only : pi, avo => avogadro


    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: m(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)

    ! local vars

    integer  :: i,k
    real(r8) :: rho_air
    real(r8) :: v, n, n_exp, r_rd, r_sd
    real(r8) :: dm_sulf, dm_sulf_wet, log_sd_sulf, sfc_sulf, sfc_nit
    real(r8) :: dm_orgc, dm_orgc_wet, log_sd_orgc, sfc_oc, sfc_soa
    real(r8) :: sfc_soai, sfc_soam, sfc_soab, sfc_soat, sfc_soax
    real(r8) :: dm_bc, dm_bc_wet, log_sd_bc, sfc_bc
    real(r8) :: rxt_sulf, rxt_nit, rxt_oc, rxt_soa
    real(r8) :: c_n2o5, c_ho2, c_no2, c_no3
    real(r8) :: s_exp

    !-----------------------------------------------------------------
    ! 	... parameters for log-normal distribution by number
    ! references:
    !   Chin et al., JAS, 59, 461, 2003
    !   Liao et al., JGR, 108(D1), 4001, 2003
    !   Martin et al., JGR, 108(D3), 4097, 2003
    !-----------------------------------------------------------------
    real(r8), parameter :: rm_sulf  = 6.95e-6_r8        ! mean radius of sulfate particles (cm) (Chin)
    real(r8), parameter :: sd_sulf  = 2.03_r8           ! standard deviation of radius for sulfate (Chin)
    real(r8), parameter :: rho_sulf = 1.7e3_r8          ! density of sulfate aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_orgc  = 2.12e-6_r8        ! mean radius of organic carbon particles (cm) (Chin)
    real(r8), parameter :: sd_orgc  = 2.20_r8           ! standard deviation of radius for OC (Chin)
    real(r8), parameter :: rho_orgc = 1.8e3_r8          ! density of OC aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_bc    = 1.18e-6_r8        ! mean radius of soot/BC particles (cm) (Chin)
    real(r8), parameter :: sd_bc    = 2.00_r8           ! standard deviation of radius for BC (Chin)
    real(r8), parameter :: rho_bc   = 1.0e3_r8          ! density of BC aerosols (kg/m3) (Chin)

    real(r8), parameter :: mw_so4 = 98.e-3_r8     ! so4 molecular wt (kg/mole)

    integer  ::  irh, rh_l, rh_u
    real(r8) ::  factor, rfac_sulf, rfac_oc, rfac_bc, rfac_ss
    logical :: zero_aerosols

    !-----------------------------------------------------------------
    ! 	... table for hygroscopic growth effect on radius (Chin et al)
    !           (no growth effect for mineral dust)
    !-----------------------------------------------------------------
    real(r8), dimension(7) :: table_rh, table_rfac_sulf, table_rfac_bc, table_rfac_oc, table_rfac_ss

    character(len=*), parameter :: subname = 'aero_model_surfarea'

    call endrun(subname//":: is not yet implemented")

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! stub
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( state, ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    !-------------------------------------------------------------------------
    ! provides WET stratospheric aerosol surface area info for modal aerosols
    ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
    ! Copied from bulk_aero
    !-------------------------------------------------------------------------

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

    character(len=*), parameter :: subname = 'aero_model_strat_surfarea'

    strato_sad(:,:) = 0._r8
    reff_strat(:,:) = 0._r8

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev,       &
      delt, reaction_rates, tfld, pmid, pdel, mbar, relhum, zm, qh2o, cwat,     &
      cldfr, cldnum, airdens, invariants, del_h2so4_gasprod, vmr0, vmr, pbuf )

    use chem_mods,    only : gas_pcnst, adv_mass
    use mo_tracname,  only: solsym
    use time_manager, only: get_nstep

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state), intent(in) :: state       ! Physics state variables
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(:)
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)              ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)            ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)            ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)          ! mixing ratios ( vmr )

    type(physics_buffer_desc), pointer :: pbuf(:)


    ! local vars
!    integer, parameter :: nmodes_aq_chem = 1
    integer  :: icol,ilev
    integer  :: l_aero
    integer  :: imode,icnst,itrac
    integer  :: nstep
    real(r8) :: wrk(ncol)
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)   ! cloud-borne aerosol (vmr)
    real(r8) :: del_h2so4_aeruptk(ncol,pver)
    real(r8) :: del_h2so4_aqchem(ncol,pver)
   ! real(r8) :: mmr_cond_vap_start_of_timestep(pcols,pver,N_COND_VAP)
   ! real(r8) :: mmr_cond_vap_gasprod(pcols,pver,N_COND_VAP)
    real(r8) :: del_soa_lv_gasprod(ncol,pver)
    real(r8) :: del_soa_sv_gasprod(ncol,pver)
    real(r8) :: dvmrdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: dvmrcwdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: mmr_tend_ncols(ncol, pver, gas_pcnst)
    real(r8) :: mmr_tend_pcols(pcols, pver, gas_pcnst)
    integer  :: cond_vap_idx
!    real(r8) :: aqso4(ncol,nmodes_aq_chem)   ! aqueous phase chemistry
!    real(r8) :: aqh2so4(ncol,nmodes_aq_chem) ! aqueous phase chemistry
 !   real(r8) :: aqso4_h2o2(ncol)             ! SO4 aqueous phase chemistry due to H2O2
 !   real(r8) :: aqso4_o3(ncol)               ! SO4 aqueous phase chemistry due to O3
    real(r8) :: xphlwc(ncol,pver)            ! pH value multiplied by lwc
    real(r8) :: delt_inverse                 ! 1 / timestep
    real(r8), pointer :: pblh(:)
    character(len=32) :: name

    integer :: l_h2so4, l_dms, l_so2
    character(len=*), parameter :: subname = 'aero_model_gasaerexch'

    nstep = get_nstep()

    delt_inverse = 1.0_r8 / delt

    ! indices of gases to condense
    call cnst_get_ind('H2SO4'  ,l_h2so4, abort=.true.)
    ! gas phase species
    call cnst_get_ind('SO2'    ,l_so2,   abort=.true.) !sulfur dioxide
    call cnst_get_ind('DMS'    ,l_dms,   abort=.true.) !dimethyl sulfide

    ! Get height of boundary layer (needed for boundary layer nucleation)
   ! call pbuf_get_field(pbuf, pblh_idx, pblh)


  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
     use oslo_aero_control, only: dms_from_ocn ! DMS
     use constituents,      only: cnst_get_ind, sflxnam
     use oslo_aero_ocean,   only: oslo_aero_dms_emis ! DMS
     use dust_model,        only: dust_emis
     use seasalt_model,     only: seasalt_emis

     ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer  :: lchnk, ncol
    integer  :: m, irange
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    integer  :: pndx_fdms  ! DMS surface flux physics index
    real(r8) :: u10cubed(pcols), u10(pcols)
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model

    character(len=*), parameter :: subname = 'aero_model_emissions'

    lchnk = state%lchnk
    ncol  = state%ncol

    if (aero_props%is_active('dust')) then

        call dust_emis( lchnk, ncol, cam_in%dstflx, cam_in%cflx, aero_props )

       ! some dust emis diagnostics ...
        sflx(:)=0._r8
        do irange = 1, aero_props%spec_nrange(dust_specprop_ndx)
            sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,aero_props%spec_mmr_q_ndx(dust_specprop_ndx,irange))
            call outfld(trim(aero_props%spec_tracernames(dust_specprop_ndx, irange))//'SF',cam_in%cflx(:,aero_props%spec_mmr_q_ndx(dust_specprop_ndx,irange)),pcols, lchnk)
        end do
    ! TODO: make this work for bins? but that should include all species..
            !        do ibin = 1, dust_nbin
    !            m = (ibin)
    !            sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,m)
    !            call outfld(trim(aero_props%bin//'SF',cam_in%cflx(:,m),pcols, lchnk)
        call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
        call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )

    endif

    if (aero_props%is_active('seasalt')) then
        u10(:ncol) = sqrt(state%u(:ncol,pver)**2+state%v(:ncol,pver)**2)
        ! move the winds to 10m high from the midpoint of the gridbox:
        ! follows Tie and Seinfeld and Pandis, p.859 with math.

        u10(:ncol)=u10(:ncol)*log(10._r8/z0)/log(state%zm(:ncol,pver)/z0)

        ! we need them to the 3.41 power, according to Gong et al., 1997:
        u10cubed(:ncol)=u10(:ncol)**3.41_r8

        sflx(:)=0._r8
        call seasalt_emis(u10cubed, u10, cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx, aero_props)

        do irange = 1, aero_props%spec_nrange(seasalt_specprop_ndx)

            sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,aero_props%spec_mmr_q_ndx(seasalt_specprop_ndx,irange))
            call outfld( trim( aero_props%spec_tracernames(seasalt_specprop_ndx, irange) )//'SF',cam_in%cflx( :,aero_props%spec_mmr_q_ndx(seasalt_specprop_ndx,irange) ),pcols, lchnk )
        end do
    end if


  end subroutine aero_model_emissions

  subroutine aero_depvel_part( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
                                     radius_part, density_part, lchnk )

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    !
    ! !USES
    !
    use physconst,     only: pi,boltz, gravit, rair
    use mo_drydep,     only: n_land_type, fraction_landuse
    use ref_pres,       only: top_lev => clim_modal_aero_top_lev

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8), intent(in) :: t(pcols,pver)       !atm temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8), intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume/number) particle radius (m)
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material (kg/m3)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk

    real(r8), intent(out) :: vlc_trb(pcols)       !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(pcols,pver)       !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(pcols,pver)       !dry deposn velocity (m/s)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k,ix                !indices
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: slp_crc(pcols,pver) ![frc] Slip correction factor
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface


    integer  :: lt
    real(r8) :: lnd_frc
    real(r8) :: wrk1, wrk2, wrk3

    ! constants
    real(r8) gamma(11)      ! exponent of schmidt number
!   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
!              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
!              0.50d+00/
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
!   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
!               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
!             100.00d+00/
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
!   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
!                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
!                         -1.00d+00/
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
!   data iwet/1,   -1,   -1,   -1,   -1,  &
!            -1,   -1,   -1,    1,   -1,  &
!             1/
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet


    vlc_trb = 0._r8
    vlc_grv = 0._r8
    vlc_dry = 0._r8

    !------------------------------------------------------------------------
    do k=top_lev,pver ! radius_part is not defined above top_lev
       do i=1,ncol

! use a maximum radius of 50 microns when calculating deposition velocity
! TODO: limit max radius??
! no dispersion for sectional

          rho=pmid(i,k)/rair/t(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(i,k) = 1.0_r8 + mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*radius_part(i,k)/(mfp_atm(i,k)))) / &
                  radius_part(i,k)   ![frc] Slip correction factor SeP97 p. 464

          vlc_grv(i,k) = (4.0_r8/18.0_r8) * radius_part(i,k)*radius_part(i,k)*density_part(i,k)* &
                  gravit*slp_crc(i,k) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466

          vlc_dry(i,k)=vlc_grv(i,k)
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       dff_aer = boltz * t(i,k) * slp_crc(i,k) / &    ![m2 s-1]
                 (6.0_r8*pi*vsc_dyn_atm(i,k)*radius_part(i,k)) !SeP97 p.474
       shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972

       wrk2 = 0._r8
       wrk3 = 0._r8
       do lt = 1,n_land_type
          lnd_frc = fraction_landuse(i,lt,lchnk)
          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(lt))
             if (radius_collector(lt) > 0.0_r8) then
!       vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) / (gravit*radius_collector(lt))
                interception = 2.0_r8*(radius_part(i,k)/radius_collector(lt))**2.0_r8
             else
!       non-vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**2.0_r8

             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(i) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i,k)

             wrk1 = 1.0_r8 / rss_trb
             wrk2 = wrk2 + lnd_frc*( wrk1 )
             wrk3 = wrk3 + lnd_frc*( wrk1 + vlc_grv(i,k) )
          endif
       enddo  ! n_land_type
       vlc_trb(i) = wrk2
       vlc_dry(i,k) = wrk3
    enddo !ncol

    return
  end subroutine aero_depvel_part

end module aero_model
