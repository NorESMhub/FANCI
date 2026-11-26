module aero_wetdep_cam
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,    only: cam_out_t
  use physics_buffer,only: physics_buffer_desc, pbuf_get_index, pbuf_set_field, pbuf_get_field
  use constituents,  only: pcnst, cnst_name, cnst_get_ind
  use phys_control,  only: phys_getopts
  use ppgrid,        only: pcols, pver
  use physconst,     only: gravit

  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use infnan,        only: nan, assignment(=)

  use cam_history,   only: addfld, add_default, horiz_only, outfld
  use wetdep,        only: wetdep_init


  use aerosol_properties_mod, only: aero_name_len
  use aerosol_properties_mod, only: aerosol_properties
  use fanci_aerosol_properties_mod, only: fanci_aerosol_properties

  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use fanci_aerosol_state_mod, only: fanci_aerosol_state

  ! Unified convective processing (aero_convproc, convproc_do_aer) is not implemented for FANCI.
  ! It could be added later; the full wet deposition side of it is in CAM's generic
  ! chemistry/aerosol/aero_wetdep_cam.F90 (and in this file's git history before option 3).
  ! aero_convproc assumes per-bin mass tracers, FANCI carries mass per range.

  implicit none                                 
  private

  public :: aero_wetdep_readnl
  public :: aero_wetdep_init
  public :: aero_wetdep_tend


  real(r8), parameter :: NOTSET = -huge(1._r8)
  real(r8) :: sol_facti_cloud_borne   = NOTSET
  real(r8) :: sol_factb_interstitial  = NOTSET
  real(r8) :: sol_factic_interstitial = NOTSET

  integer :: fracis_idx = -1


  logical :: wetdep_active = .false.
  integer :: nwetdep = 0
  logical :: convproc_do_aer = .false.   ! only read to stop the run if it is set (not supported)

  logical,allocatable :: aero_cnst_lq(:,:)
  integer,allocatable :: aero_cnst_id(:,:)
  logical, public, protected :: wetdep_lq(pcnst) ! set flags true for constituents with non-zero tendencies

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8),allocatable :: scavimptblnum(:,:)
  real(r8),allocatable :: scavimptblvol(:,:)

  integer :: nele_tot            ! total number of aerosol elements
  class(aerosol_properties), pointer :: aero_props=>null()

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_character, mpi_real8, mpi_integer, mpi_success
    use spmd_utils,      only: mpi_logical

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'aero_wetdep_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /aero_wetdep_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aero_wetdep_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aero_wetdep_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)

       ! ============================
       ! Log namelist options
       ! ============================
       write(iulog,*) subname,' namelist settings: '
       write(iulog,*) '   sol_facti_cloud_borne  : ',sol_facti_cloud_borne
       write(iulog,*) '   sol_factb_interstitial : ',sol_factb_interstitial
       write(iulog,*) '   sol_factic_interstitial: ',sol_factic_interstitial
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(sol_facti_cloud_borne, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_facti_cloud_borne')
    end if
    call mpi_bcast(sol_factb_interstitial, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_factb_interstitial')
    end if
    call mpi_bcast(sol_factic_interstitial, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_factic_interstitial')
    end if

    call mpi_bcast(nwetdep, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: nwetdep')
    end if

    wetdep_active = .true. !nwetdep>0

    if (masterproc) then
       write(iulog,*) subname,' wetdep_active = ',wetdep_active,' nwetdep = ',nwetdep
    endif

    !call aero_convproc_readnl(nlfile)

  end subroutine aero_wetdep_readnl

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_init( )

    character(len=*), parameter :: subrname = 'aero_wetdep_init'

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    character(len=aero_name_len) :: tmpname
    character(len=aero_name_len) :: tmpname_cw

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry

    integer :: l,m, id, astat
    character(len=2) :: binstr

    fracis_idx = pbuf_get_index('FRACIS')

    if (.not.wetdep_active) return

    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      convproc_do_aer_out = convproc_do_aer)


    aero_props => fanci_aerosol_properties()
    if (.not.associated(aero_props)) then
        call endrun(subrname//' : construction of aero_props fanci_aerosol_properties object failed')
    end if
    if (convproc_do_aer) then
        call endrun(subrname// &
                ' : convproc_do_aer=.true. is not supported for the fanci scheme: '// &
                   'unified convective processing (aero_convproc) is not implemented in FANCI. '// &
                'Set convproc_do_aer=.false.')
    end if

    nele_tot = aero_props%ncnst_tot()

    allocate(aero_cnst_lq(aero_props%nbins(),0:maxval(aero_props%nspecies())), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate aero_cnst_lq array')
    end if
    aero_cnst_lq(:,:) = .false.

    allocate(aero_cnst_id(aero_props%nbins(),0:maxval(aero_props%nspecies())), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate aero_cnst_id array')
    end if
    aero_cnst_id(:,:) = -1

    wetdep_lq = .false.

    do m = 1, aero_props%nbins()
       write(binstr,'(i2.2)') m
       call addfld('SOLFACTB'//binstr,  (/ 'lev' /), 'A', '1', 'below cld sol fact')

       do l = 0, aero_props%nspecies(m)

          if (l == 0) then   ! number
             call aero_props%num_names( m, tmpname, tmpname_cw)
          else
             call aero_props%mmr_names( m,l, tmpname, tmpname_cw)
          end if

          call cnst_get_ind(tmpname, id, abort=.false.)
          aero_cnst_id(m,l) = id
          aero_cnst_lq(m,l) = id > 0
          if (id > 0) then
             wetdep_lq(id) = .true.
          end if

          ! units --
          if (l==0) then
             unit_basename = ' 1' ! for num
          else
             unit_basename = 'kg'
          endif

          if (.not. tmpname == '') then
              call add_hist_fields(tmpname, unit_basename)
              call add_hist_fields(tmpname_cw, unit_basename)

          end if
       end do
    end do

    allocate(scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, aero_props%nbins()), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate scavimptblnum array')
    end if
    allocate(scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, aero_props%nbins()), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate scavimptblvol array')
    end if
    scavimptblnum = nan
    scavimptblvol = nan

    call wetdep_init()


    call init_bcscavcoef()


  contains

    subroutine add_hist_fields(name,baseunits)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: baseunits

if (masterproc) then
         write(iulog,*) subrname//' : adding history fields for ',trim(name)
      end if

      call addfld (trim(name)//'SFWET', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux at surface')
      call addfld (trim(name)//'SFWETC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (convective) at surface')
      call addfld (trim(name)//'SFSIC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (incloud, convective) at surface')
      call addfld (trim(name)//'SFSIS', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (incloud, stratiform) at surface')
      call addfld (trim(name)//'SFSBC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (belowcloud, convective) at surface')
      call addfld (trim(name)//'SFSBS', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (belowcloud, stratiform) at surface')

      call addfld (trim(name)//'SFSEC', &
              horiz_only,  'A',unit_basename//'/m2/s','Wet deposition flux (precip evap, convective) at surface')


      call addfld (trim(name)//'WET',(/ 'lev' /), 'A',baseunits//'/kg/s ','wet deposition tendency')
      call addfld (trim(name)//'INS',(/ 'lev' /), 'A',baseunits//'/kg/s ','insol frac')

      call addfld (trim(name)//'SIC',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' ic wet deposition')
      call addfld (trim(name)//'SIS',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' is wet deposition')
      call addfld (trim(name)//'SBC',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' bc wet deposition')
      call addfld (trim(name)//'SBS',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' bs wet deposition')

      if ( history_aerosol .or. history_chemistry ) then
         call add_default (trim(name)//'SFWET', 1, ' ')
      endif
      if ( history_aerosol ) then
         call add_default (trim(name)//'SFSEC', 1, ' ')
         call add_default (trim(name)//'SFSIC', 1, ' ')
         call add_default (trim(name)//'SFSIS', 1, ' ')
         call add_default (trim(name)//'SFSBC', 1, ' ')
         call add_default (trim(name)//'SFSBS', 1, ' ')

      endif

    end subroutine add_hist_fields

  end subroutine aero_wetdep_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_tend( state, dt, dlf, cam_out, ptend, pbuf)
    use wetdep, only: wetdepa_v2, wetdep_inputs_set, wetdep_inputs_t
    use aerodep_flx, only: aerodep_flx_prescribed
    use aero_deposition_cam, only: aero_deposition_cam_setwet

    type(physics_state), target, intent(in) :: state  ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state, gets the surface wet deposition fluxes
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer of this chunk

    character(len=*), parameter :: subrname = 'aero_wetdep_tend'
    type(wetdep_inputs_t) :: dep_inputs                !cloud and precipitation fields gathered by wetdep_inputs_set: cloud fractions, rain production, evaporation, precipitation areas
    real(r8), pointer :: fracis(:,:,:)                 ! fraction of transported species that are not scavenged (pcols, pver, pcnst)
    real(r8), target :: fracis_nadv(pcols,pver)        ! fraction of not-transported aerosols (cloud-borne)

    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)


    real(r8) :: sol_factb(pcols, pver)    ! below-cloud solubility factor (fraction that can be captured)


    real(r8) :: dqdt_tmp(pcols,pver)      ! temporary array to hold tendency of one tracer from one call

    real(r8) :: iscavt(pcols, pver)       ! total in-cloud scavenging; an output that is never used afterwards
    real(r8) :: icscavt(pcols, pver)      ! in-cloud scavenging: convective(diagnostics)
    real(r8) :: isscavt(pcols, pver)      ! in-cloud scavenging: stratiform (diagnostics)
    real(r8) :: bcscavt(pcols, pver)      ! below-cloud scavenging: convective(diagnostics)
    real(r8) :: bsscavt(pcols, pver)      ! below-cloud scavenging: stratiform (diagnostics)

    real(r8) :: diam_wet(state%ncol, pver)! wet diameter of the bin, which picks the impaction coefficient
    logical  :: isprx(pcols,pver)         ! true where precipitation falls into the layer from above
    real(r8) :: prec(pcols)               ! precipitation flux summed from the top down, used to set isprx

    !real(r8) :: rtscavt(pcols, pver, 0:nspec_max) ! saved cloud-borne resuspension

    integer :: ncol, lchnk, ndx, mm
    integer :: i,k

    real(r8) :: q_tmp(pcols,pver)         ! temporary array to hold "most current" mixing ratio for 1 species

    type(ptr2d_t) :: raer(nele_tot)       !pointers to interstitial tracers (into state%q)
    type(ptr2d_t) :: qqcw(nele_tot)       !pointers to cloud-borne tracers (in the state object)

    character(len=aero_name_len) :: aname, cname   ! interstitial / cloud-borne tracer name

    real(r8) :: qqcw_in(pcols,pver)       ! the cloud-borne amount that matches the interstitial tracer, passed as qqcw
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud

    character(len=2) :: binstr            ! bin number as 2 characters, for SOLFACTB01…
    real(r8) :: aerdepwetcw(pcols,pcnst)  ! aerosol surface wet deposition (cloud water) sent to cam_out
    real(r8) :: aerdepwetis(pcols,pcnst)  ! aerosol surface wet deposition (interstitial) sent to cam_out


    class(aerosol_state), pointer :: aero_state ! the FANCI master state for this chunk (via the abstract class)
    type(fanci_aerosol_state), pointer :: fanci_aero_state ! same object, FANCI type (for FANCI-only methods)
    ! FANCI interstitial mass (per-bin wet deposition summed onto the range)
    type(fanci_aerosol_properties), pointer :: fanci_aero_props ! FANCI properties with the FANCI type, so nranges, range_bounds and so on are available
    integer  :: ibin, irange, ispec, b1, b2    ! bin, range, species and first and last bin in range (b1,b2)
    real(r8) :: bin_mmr_tot(pcols,pver)        ! total dry mass of one bin (kg/kg)
    real(r8), allocatable :: bin_cw_mmr(:,:,:) ! (pcols,pver,nbins) cloud-borne mass of each bin, before any cloud-borne removal
    real(r8) :: bin_fracis(pcols,pver)         ! fraction not scavenged, from the bin's mass call
    real(r8) :: range_mmr(pcols,pver) !interstitial mass of the range, summed over species
    real(r8) :: range_mmr_tend(pcols,pver)   ! mass tendency of the range, summed over bins
    real(r8) :: range_fracis(pcols,pver), range_wsum(pcols,pver) !range_fracis: mass-weighted fraction not scavenged; range_wsum: sum of bin mass (its weight)
    real(r8) :: massfrac(pcols,pver)           ! mass fraction of one species in the range
    real(r8), pointer :: mmr_s(:,:)
    real(r8) :: range_icscavt(pcols,pver), range_isscavt(pcols,pver)   ! diagnostics summed over bins
    real(r8) :: range_bcscavt(pcols,pver), range_bsscavt(pcols,pver)

    nullify(fanci_aero_props)
    nullify(aero_state)
    nullify(fanci_aero_state)

    if (.not.wetdep_active) return


    aero_props => fanci_aerosol_properties()
    fanci_aero_state => fanci_aerosol_state(state, pbuf)
    aero_state       => fanci_aero_state

    lchnk = state%lchnk
    ncol = state%ncol
   ! TODO: remove cam dependency, initialises ptend
    call physics_ptend_init(ptend, state%psetcols, subrname, lq=wetdep_lq)
   ! TODO: remove cam dependency. This code gets input for wet deposition from host model
    call wetdep_inputs_set( state, pbuf, dep_inputs )
   ! TODO: remove cam dependecy. This is fraction not scavenged, per constituent; read later by e.g. deep convection
    call pbuf_get_field(pbuf, fracis_idx, fracis)
   ! Gets the interstitial (raer) and cloud borne (qqcw) aerosols
    call aero_state%get_states( aero_props, raer, qqcw )
   

    ! initialises aerosol surface wet deposition (interstitial and cloud borne)
    aerdepwetis(:,:) = 0.0_r8
    aerdepwetcw(:,:) = 0.0_r8

   ! Calculating summed precipitation from top down.
    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true. ! there is precip
       elsewhere
          isprx(:ncol,k) = .false.  ! there is not precip
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
                    *state%pdel(:ncol,k)/gravit
    end do

    f_act_conv = 0._r8 ! Initializing prescribed aerosol activation fraction for convective cloud 
    scavcoefnv = nan    ! initializing Dana and Hales coefficient (/mm) for cloud-borne num & vol (0), interstitial num (1), interstitial vol (2)


    ! Unified convective processing (convproc_do_aer) is not implemented for FANCI;
    ! aero_wetdep_init stops if it is set. See the note at the top of this module.



    ! FANCI wet deposition: a loop over ranges, and over the bins of each range.
    !   bin_loop (each bin, with the bin's own size, coefficients and factors):
    !     1. interstitial number                     -> ptend
    !     2. interstitial total dry mass of the bin  -> summed onto the range
    !     3. cloud-borne number                      -> in place (state object)
    !   species loop (each species of the range):
    !     4. interstitial mass: range tendency split by mass fraction -> ptend
    !     5. cloud-borne mass (does not depend on size)               -> in place
    ! The interstitial calls see cloud-borne amounts from before any cloud-borne removal:
    ! bin_cw_mmr is computed for all bins first, and a bin's cloud-borne number is removed
    ! only after that bin's interstitial calls.

    ! The solubility factors for each phase are set in scavenge (see contains below).
    fanci_aero_props => fanci_aerosol_properties()


    ! cloud-borne mass of each bin (its share of the range's cloud-borne mass).
    ! Computed for all bins BEFORE any cloud-borne removal: cldbrne_total_bin_mmr
    ! sums the cloud-borne number of all bins in the range.
    allocate(bin_cw_mmr(pcols, pver, aero_props%nbins()))
    bin_cw_mmr = 0._r8
    do ibin = 1, aero_props%nbins()
        do k = 1, pver
            do i = 1, ncol
                bin_cw_mmr(i,k,ibin) = fanci_aero_state%cldbrne_total_bin_mmr(aero_props, ibin, i, k)
            end do
        end do
    end do

    range_loop: do irange = 1, fanci_aero_props%nranges()

        b1 = fanci_aero_props%range_bounds(irange,1)   ! leader bin: owns the mass tracers
        b2 = fanci_aero_props%range_bounds(irange,2)

        ! range totals; the mass tracers live on the leader bin
        range_mmr    = 0._r8
        do ispec = 1, fanci_aero_props%range_nspecies(irange)
            call aero_state%get_ambient_mmr(ispec, b1, mmr_s)
            range_mmr(:ncol,:) = range_mmr(:ncol,:) + mmr_s(:ncol,:)
        end do


        range_mmr_tend = 0._r8
        range_fracis   = 0._r8
        range_wsum     = 0._r8
        range_icscavt  = 0._r8
        range_isscavt  = 0._r8
        range_bcscavt  = 0._r8
        range_bsscavt  = 0._r8

        bin_loop: do ibin = b1, b2
            !this bin's own total dry mass
            do k = 1, pver
                do i = 1, ncol
                    bin_mmr_tot(i,k) = aero_state%ambient_total_bin_mmr(aero_props, ibin, i, k)
                end do
            end do
            ! this bin's own size, coefficients and factors (used for both number and mass)
            diam_wet = aero_state%wet_diameter(ibin, ncol, pver)
            scavcoefnv = 0.0_r8
            ! below cloud impaction scavenging coefs:
            ! get scav coef for interstitial number and volume, index 0 (cloud-borne) stays 0 (cloud borne aerosols are not below cloud scavenged)
            call get_bcscavcoefs( ibin, ncol, isprx, diam_wet, scavcoefnv(:,:,1), scavcoefnv(:,:,2) )

            ! below-cloud solubility factor (fraction that can be captured)
            if ( sol_factb_interstitial /= NOTSET ) then
                sol_factb(:ncol,:) = sol_factb_interstitial
            else
                sol_factb(:ncol,:) = aero_state%sol_factb_interstitial( ibin, ncol, pver, aero_props )
            end if
            write(binstr,'(i2.2)') ibin
            call outfld('SOLFACTB'//binstr, sol_factb, pcols, lchnk)
            !fraction activated in convective cloud
            f_act_conv(:ncol,:) = aero_state%convcld_actfrac( ibin, 1, ncol, pver )

            ! interstitial NUMBER of this bin: same size, coefficients and factors as the
            ! mass call below, but with the number impaction coefficient (scavcoefnv(:,:,1))
            ndx = aero_cnst_id(ibin, 0)
            if (ndx > 0) then
                mm = aero_props%indexer(ibin, 0)  ! index 0 is number
                call aero_props%num_names( ibin, aname, cname )
                q_tmp(:ncol,:)   = raer(mm)%fld(:ncol,:) + ptend%q(:ncol,:,ndx)*dt  ! most current number
                qqcw_in(:ncol,:) = qqcw(mm)%fld(:ncol,:)                            ! matching cloud-borne number
               ! Also fills icscavt, isscavt, bcscavt and bsscavt (declared in aero_wetdep_tend):
                call scavenge( .false., &                   ! .true. for stratiform cloud-borne aerosol
                              q_tmp, qqcw_in, &             ! number (#/kg) or mass (kg/kg) + equivalent for cloud-borne
                              scavcoefnv(:,:,1), &          ! scavenging coefficient for below-cloud impaction, here for number
                              sol_factb, f_act_conv, &      ! solubility factors below cloud, activation factor convective
                              dqdt_tmp, &                   ! tendency of tracer (out)
                              fracis(:,:,ndx) )             ! fraction not scavenged (out)
                ptend%q(:ncol,:,ndx) = ptend%q(:ncol,:,ndx) + dqdt_tmp(:ncol,:)
                call wetdep_output( aname, ndx, .false., dqdt_tmp, icscavt, isscavt, bcscavt, bsscavt, &
                                    fracis(:,:,ndx) )
            end if

            ! interstitial MASS of this bin (total dry mass, volume impaction coefficient).
            ! wetdepa_v2 goes through the column from the top and follows the convective and
            ! stratiform precipitation and the tracer it carries; where rain evaporates, the
            ! same fraction of the carried tracer is released. Four processes: in-cloud and
            ! below-cloud, convective and stratiform.
            ! Also fills icscavt, isscavt, bcscavt and bsscavt (declared in aero_wetdep_tend):
            call scavenge( .false., &             ! not cloud-borne
                    bin_mmr_tot, &                ! bin mass mmr
                    bin_cw_mmr(:,:,ibin), &       ! its cloud-borne mass
                    scavcoefnv(:,:,2), &          ! volume (mass) impaction coefficient
                    sol_factb, f_act_conv, &
                    dqdt_tmp, bin_fracis )        ! outputs

            ! add up the bins of the range
            range_mmr_tend(:ncol,:) = range_mmr_tend(:ncol,:) + dqdt_tmp(:ncol,:)
            range_fracis(:ncol,:)   = range_fracis(:ncol,:) + bin_fracis(:ncol,:)*bin_mmr_tot(:ncol,:)
            range_wsum(:ncol,:)     = range_wsum(:ncol,:)   + bin_mmr_tot(:ncol,:)
            range_icscavt(:ncol,:)  = range_icscavt(:ncol,:) + icscavt(:ncol,:)
            range_isscavt(:ncol,:)  = range_isscavt(:ncol,:) + isscavt(:ncol,:)
            range_bcscavt(:ncol,:)  = range_bcscavt(:ncol,:) + bcscavt(:ncol,:)
            range_bsscavt(:ncol,:)  = range_bsscavt(:ncol,:) + bsscavt(:ncol,:)

            ! cloud-borne NUMBER of this bin (stratiform in-cloud only; does not depend on size).
            ! Done last for this bin: the interstitial calls above must see the cloud-borne (to determine clear air aerosol)
            ! number from before removal. Applied in place to the state object.
            if (ndx > 0) then
                q_tmp(:ncol,:) = qqcw(mm)%fld(:ncol,:)          ! copy: fld is (ncol,pver)
                ! Also fills icscavt, isscavt, bcscavt and bsscavt (declared in aero_wetdep_tend):
                call scavenge( .true., q_tmp, q_tmp, &          ! tracer_cw not used for cloud-borne
                               scavcoefnv(:,:,0), &             ! 0: no below-cloud for cloud-borne
                               sol_factb, f_act_conv, &         ! not used for cloud-borne
                               dqdt_tmp, fracis_nadv )
                call update_cldbrne( qqcw(mm)%fld, dqdt_tmp )   ! clip at zero and apply
                call wetdep_output( cname, ndx, .true., dqdt_tmp, icscavt, isscavt, bcscavt, bsscavt, &
                                    fracis_nadv )
            end if
        end do bin_loop

        ! mass weighted fraction not scavenged (used by convective transport)
        where (range_wsum(:ncol,:) > 0._r8)
            range_fracis(:ncol,:) = range_fracis(:ncol,:)/range_wsum(:ncol,:)
        elsewhere
            range_fracis(:ncol,:) = 1._r8
        end where

        ! split the range tendency over the species and apply it
        do ispec = 1, fanci_aero_props%range_nspecies(irange)

            ndx = aero_cnst_id(b1, ispec)
            if (ndx < 1) cycle

            call aero_state%get_ambient_mmr(ispec, b1, mmr_s)
            massfrac = 0._r8
            where (range_mmr(:ncol,:) > 0._r8) &
                    massfrac(:ncol,:) = mmr_s(:ncol,:)/range_mmr(:ncol,:)

            dqdt_tmp(:ncol,:)    = range_mmr_tend(:ncol,:)*massfrac(:ncol,:)
            ptend%q(:ncol,:,ndx) = ptend%q(:ncol,:,ndx) + dqdt_tmp(:ncol,:)
            fracis(:ncol,:,ndx)  = range_fracis(:ncol,:)

            call aero_props%mmr_names( b1, ispec, aname, cname)

            call wetdep_output( aname, ndx, .false., dqdt_tmp,       &
                    range_icscavt*massfrac, range_isscavt*massfrac,  &
                    range_bcscavt*massfrac, range_bsscavt*massfrac, range_fracis )
            ! cloud-borne MASS of this species (stratiform in-cloud only). Removal does not
            ! depend on size, so it is done once per species on the range, not per bin.
            ! Applied in place to the state object.
            mm = aero_props%indexer(b1, ispec)              ! position in raer/qqcw of this species
            q_tmp(:ncol,:) = qqcw(mm)%fld(:ncol,:)          ! copy: fld is (ncol,pver)
            ! Also fills icscavt, isscavt, bcscavt and bsscavt (declared in aero_wetdep_tend):
            call scavenge( .true.,     &           ! cloud borne true
                           q_tmp,      &           ! tracer mass
                           q_tmp,      &           ! tracer_cw not used for cloud-borne
                           scavcoefnv(:,:,0), &             ! 0: no below-cloud for cloud-borne
                           sol_factb, f_act_conv, &         ! not used for cloud-borne
                           dqdt_tmp, fracis_nadv )          ! outputs
            call update_cldbrne( qqcw(mm)%fld, dqdt_tmp )   ! clip at zero and apply
            call wetdep_output( cname, ndx, .true., dqdt_tmp, icscavt, isscavt, bcscavt, bsscavt, &
                                fracis_nadv )

        end do ! species loop
    end do range_loop
    deallocate(bin_cw_mmr)


    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not. aerodep_flx_prescribed()) then
       call aero_deposition_cam_setwet(aerdepwetis, aerdepwetcw, cam_out)
    endif

  contains
        ! one wetdepa_v2 call; solubility factors set from the phase
        ! Also fills icscavt, isscavt, bcscavt and bsscavt (declared in aero_wetdep_tend).
        subroutine scavenge( cldbrn, tracer, tracer_cw, scavcoef, sfb_bin, factconv_bin, dqdt, insol )

            logical,  intent(in)  :: cldbrn                   ! .true. for stratiform cloud-borne aerosol
            real(r8), intent(in)  :: tracer(pcols,pver)       ! number (#/kg) or mass (kg/kg)
            real(r8), intent(in)  :: tracer_cw(pcols,pver)    ! matching cloud-borne amount (used only for interstitial calls)
            real(r8), intent(in)  :: scavcoef(pcols,pver)     ! below-cloud impaction coefficient
            real(r8), intent(in)  :: sfb_bin(pcols,pver)      ! below-cloud solubility factor of the bin
            real(r8), intent(in)  :: factconv_bin(pcols,pver) ! convective activation fraction of the bin
            real(r8), intent(out) :: dqdt(pcols,pver)         ! tendency of tracer
            real(r8), intent(out) :: insol(pcols,pver)        ! fraction not scavenged

            real(r8) :: sfb(pcols,pver), sfi(pcols,pver), sfic(pcols,pver), fac(pcols,pver)

            if ( cldbrn ) then
                sfb  = 0._r8                     ! no below-cloud: cloud-borne is in cloud
                sfi  = sol_facti_cloud_borne     ! stratiform in-cloud
                sfic = 0._r8                     ! no convective in-cloud
                fac  = 0._r8
            else
                sfb  = sfb_bin
                sfi  = 0._r8                     ! no stratiform in-cloud: activated aerosol is cloud-borne
                sfic = sol_factic_interstitial
                fac  = factconv_bin
            end if

            dqdt(:ncol,:) = 0._r8

            call wetdepa_v2(state%pmid, state%q(:,:,1), state%pdel, &
                    dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                    dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                    dep_inputs%evapr, dep_inputs%totcond, tracer, dt, &
                    dqdt, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                    dlf, insol, sfb, ncol, scavcoef, &
                    is_strat_cloudborne=cldbrn, qqcw=tracer_cw, f_act_conv=fac, &
                    icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                    sol_facti_in=sfi, sol_factic_in=sfic )
            ! Not passed (convproc only, not implemented for FANCI): convproc_do_aer, rcscavt,
            ! rsscavt, convproc_do_evaprain_atonce_in, bergso_in. Without them wetdepa_v2 uses
            ! partial-evaporation resuspension and no Bergeron correction (as FANCI did before).

        end subroutine scavenge

        ! apply a cloud-borne tendency in place, without going below zero
        subroutine update_cldbrne( fld, dqdt )

            real(r8), intent(inout) :: fld(:,:)            ! cloud-borne field (in the state object)
            real(r8), intent(inout) :: dqdt(pcols,pver)    ! tendency; reduced where it would go negative

            where ( fld(:ncol,:) + dqdt(:ncol,:)*dt < 0._r8 ) dqdt(:ncol,:) = -fld(:ncol,:)/dt
            fld(:ncol,:) = fld(:ncol,:) + dqdt(:ncol,:)*dt

        end subroutine update_cldbrne

        ! history output and surface fluxes for one tracer
        subroutine wetdep_output( name, ndx, cldbrn, dqdt, ic, is, bc, bs, insol )

            character(len=*), intent(in) :: name   ! tracer name (interstitial or cloud-borne)
            integer,  intent(in) :: ndx            ! constituent index of the interstitial tracer
            logical,  intent(in) :: cldbrn         ! .true. for cloud-borne
            real(r8), intent(in) :: dqdt(pcols,pver)
            real(r8), intent(in) :: ic(pcols,pver), is(pcols,pver), bc(pcols,pver), bs(pcols,pver)
            real(r8), intent(in) :: insol(pcols,pver)

            real(r8) :: flx(pcols)

            call outfld( trim(name)//'WET', dqdt,  pcols, lchnk )
            call outfld( trim(name)//'SIC', ic,    pcols, lchnk )
            call outfld( trim(name)//'SIS', is,    pcols, lchnk )
            call outfld( trim(name)//'SBC', bc,    pcols, lchnk )
            call outfld( trim(name)//'SBS', bs,    pcols, lchnk )
            call outfld( trim(name)//'INS', insol, pcols, lchnk )

            call column_integral( dqdt, flx )
            if ( cldbrn ) then
                call outfld( trim(name)//'SFWET', flx, pcols, lchnk )
                aerdepwetcw(:ncol,ndx) = flx(:ncol)
            else
                aerdepwetis(:ncol,ndx) = aerdepwetis(:ncol,ndx) + flx(:ncol)
                call outfld( trim(name)//'SFWET', aerdepwetis(:ncol,ndx), ncol, lchnk )
            end if

            call column_integral( ic, flx )
            call outfld( trim(name)//'SFSIC', flx, pcols, lchnk )
            call column_integral( is, flx )
            call outfld( trim(name)//'SFSIS', flx, pcols, lchnk )
            call column_integral( bc, flx )
            call outfld( trim(name)//'SFSBC', flx, pcols, lchnk )
            call column_integral( bs, flx )
            call outfld( trim(name)//'SFSBS', flx, pcols, lchnk )

        end subroutine wetdep_output

        ! column integral of a tendency (kg/kg/s) -> flux (kg/m2/s)
        subroutine column_integral( tend, flux )

            real(r8), intent(in)  :: tend(:,:)
            real(r8), intent(out) :: flux(:)

            integer :: kk

            flux = 0._r8
            do kk = 1, pver
                flux(:ncol) = flux(:ncol) + tend(:ncol,kk)*state%pdel(:ncol,kk)/gravit
            end do

        end subroutine column_integral

        ! below cloud impaction scavenging coefs
        subroutine get_bcscavcoefs( m, ncol, isprx, diam_wet, scavcoefnum, scavcoefvol )

            integer,intent(in) :: m, ncol
            logical,intent(in):: isprx(:,:)
            real(r8), intent(in) :: diam_wet(:,:)
            real(r8), intent(out) :: scavcoefnum(:,:), scavcoefvol(:,:)

            integer i, k, jgrow
            real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum

            do k = 1, pver
               do i = 1, ncol

                  ! do only if no precip
                  if ( isprx(i,k) .and. diam_wet(i,k)>0.0_r8) then
                     !
                     ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

                     dumdgratio = diam_wet(i,k)/aero_props%scav_diam(m)
                     if ((dumdgratio >= 0.99_r8) .and. (dumdgratio <= 1.01_r8)) then
                        scavimpvol = scavimptblvol(0,m)
                        scavimpnum = scavimptblnum(0,m)
                     else
                        xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                        jgrow = int( xgrow )
                        if (xgrow < 0._r8) jgrow = jgrow - 1
                        if (jgrow < nimptblgrow_mind) then
                           jgrow = nimptblgrow_mind
                           xgrow = jgrow
                        else
                           jgrow = min( jgrow, nimptblgrow_maxd-1 )
                        end if

                        dumfhi = xgrow - jgrow
                        dumflo = 1._r8 - dumfhi

                        scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                                     dumfhi*scavimptblvol(jgrow+1,m)
                        scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                                     dumfhi*scavimptblnum(jgrow+1,m)

                     end if

                     ! impaction scavenging removal amount for volume
                     scavcoefvol(i,k) = exp( scavimpvol )
                     ! impaction scavenging removal amount to number
                     scavcoefnum(i,k) = exp( scavimpnum )

                  else
                     scavcoefvol(i,k) = 0._r8
                     scavcoefnum(i,k) = 0._r8
                  end if

               end do
            end do

        end subroutine get_bcscavcoefs

  end subroutine aero_wetdep_tend

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine init_bcscavcoef( )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    ! Simone Tilmes Nov 2021
    ! added modifications for bin model, assuming sigma = 1.
    !
    !-----------------------------------------------------------------------

    use mo_constants,   only: pi

    !   local variables
    integer nnfit_maxd
    parameter (nnfit_maxd=27)

    integer m, jgrow, nnfit
    integer lunerr

    real(r8) dg0, dg0_cgs, press, dg0_base, &
         rhodryaero, rhowetaero, rhowetaero_cgs, &
         scavratenum, scavratevol, logsig,                &
         temp, wetdiaratio, wetvolratio

    real(r8) :: xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
    real(r8) :: xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)

    character(len=*), parameter :: subname = 'aero_wetdep_cam::init_bcscavcoef'

    lunerr = iulog
    dlndg_nimptblgrow = log( 1.25_r8 )

    ! bin model: main loop over aerosol bins

    modeloop: do m = 1, aero_props%nbins()

       ! for setting up the lookup table, use the dry density of the first species
       ! -- assume the first species of the mode/bin is the dominate species
       call aero_props%get(m,1,density=rhodryaero)

       dg0_base = aero_props%scav_diam(m)

       logsig = aero_props%alogsig(m)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          dg0 = dg0_base*wetdiaratio

          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

          !
          !   compute impaction scavenging rates at 1 temp-press pair and save
          !
          nnfit = 0

          temp = 273.16_r8
          press = 0.75e6_r8   ! dynes/cm2
          rhowetaero = rhodryaero

          dg0_cgs = dg0*1.0e2_r8   ! m to cm

          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3

          call calc_1_impact_rate( &
               dg0_cgs, logsig, rhowetaero_cgs, temp, press, &
               scavratenum, scavratevol, lunerr )

          nnfit = nnfit + 1
          if (nnfit > nnfit_maxd) then
             write(lunerr,9110)
             call endrun(subname//' : nnfit > nnfit_maxd')
          end if
9110      format( '*** subr. init_bcscavcoef -- nnfit too big' )

          xxfitnum(1,nnfit) = 1._r8
          yyfitnum(nnfit) = log( scavratenum )

          xxfitvol(1,nnfit) = 1._r8
          yyfitvol(nnfit) = log( scavratevol )

         !depends on both bins and different species
          scavimptblnum(jgrow,m) = yyfitnum(1)
          scavimptblvol(jgrow,m) = yyfitvol(1)

       enddo growloop
    enddo modeloop

  contains

    !===============================================================================
    subroutine calc_1_impact_rate(          &
         dg0, logsig, rhoaero, temp, press, &
         scavratenum, scavratevol, lunerr )
      !
      !   routine computes a single impaction scavenging rate
      !	for precipitation rate of 1 mm/h
      !
      !   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
      !   sigmag = geometric standard deviation of size distrib.
      !   rhoaero = density of aerosol particles (g/cm^3)
      !   temp = temperature (K)
      !   press = pressure (dyne/cm^2)
      !   scavratenum = number scavenging rate (1/h)
      !   scavratevol = volume or mass scavenging rate (1/h)
      !   lunerr = logical unit for error message
      !
      use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, rgas => rgas_cgs

      implicit none

      !   subr. parameters
      integer, intent(in) :: lunerr
      real(r8), intent(in) :: dg0, logsig, rhoaero, temp, press
      real(r8), intent(out) :: scavratenum, scavratevol

      !   local variables
      integer nrainsvmax
      parameter (nrainsvmax=50)
      real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
           vfallrainsv(nrainsvmax)

      integer naerosvmax
      parameter (naerosvmax=51)
      real(r8) aaerosv(naerosvmax), &
           ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

      integer i, ja, jr, na, nr
      real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
      real(r8) anumsum, avolsum, cair, chi
      real(r8) d, dr, dum, dumfuchs, dx
      real(r8) ebrown, eimpact, eintercept, etotal, freepath
      real(r8) precip, precipmmhr, precipsum
      real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
      real(r8) scavsumnum, scavsumnumbb
      real(r8) scavsumvol, scavsumvolbb
      real(r8) schmidt, sqrtreynolds, sstar, stokes, sx
      real(r8) taurelax, vfall, vfallstp
      real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair

      rlo = .005_r8
      rhi = .250_r8
      dr = 0.005_r8
      nr = 1 + nint( (rhi-rlo)/dr )
      if (nr > nrainsvmax) then
         write(lunerr,9110)
         call endrun(subname//' : nr > nrainsvmax')
      end if

9110  format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

      precipmmhr = 1.0_r8
      precip = precipmmhr/36000._r8

      ag0 = dg0/2._r8
      sx = logsig
      xg0 = log( ag0 )
      if (aero_props%model_is('fanci')) then
          na=1
      else

        xg3 = xg0 + 3._r8*sx*sx

        xlo = xg3 - 4._r8*sx
        xhi = xg3 + 4._r8*sx
        dx = 0.2_r8*sx

        dx = max( 0.2_r8*sx, 0.01_r8 )
        xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
        xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

        na = 1 + nint( (xhi-xlo)/dx )
      end if

      if (na > naerosvmax) then
         write(lunerr,9120)
         call endrun(subname//' : na > naerosvmax')
      end if

9120  format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

      !   air molar density
      cair = press/(rgas*temp)
      !   air mass density
      rhoair = 28.966_r8*cair
      !   molecular freepath
      freepath = 2.8052e-10_r8/cair
      !   air dynamic viscosity
      airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
           ((temp/296.16_r8)**1.5_r8)
      !   air kinemaic viscosity
      airkinvisc = airdynvisc/rhoair
      !   ratio of water viscosity to air viscosity (from Slinn)
      xmuwaterair = 60.0_r8

      !
      !   compute rain drop number concentrations
      !	rrainsv = raindrop radius (cm)
      !	xnumrainsv = raindrop number concentration (#/cm^3)
      !		(number in the bin, not number density)
      !	vfallrainsv = fall velocity (cm/s)
      !
      precipsum = 0._r8
      do i = 1, nr
         r = rlo + (i-1)*dr
         rrainsv(i) = r
         xnumrainsv(i) = exp( -r/2.7e-2_r8 )

         d = 2._r8*r
         if (d <= 0.007_r8) then
            vfallstp = 2.88e5_r8 * d**2._r8
         else if (d <= 0.025_r8) then
            vfallstp = 2.8008e4_r8 * d**1.528_r8
         else if (d <= 0.1_r8) then
            vfallstp = 4104.9_r8 * d**1.008_r8
         else if (d <= 0.25_r8) then
            vfallstp = 1812.1_r8 * d**0.638_r8
         else
            vfallstp = 1069.8_r8 * d**0.235_r8
         end if

         vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
         vfallrainsv(i) = vfall
         precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
      end do
      precipsum = precipsum*pi*1.333333_r8

      rnumsum = 0._r8
      do i = 1, nr
         xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
         rnumsum = rnumsum + xnumrainsv(i)
      end do

      !
      !   compute aerosol concentrations
      !	aaerosv = particle radius (cm)
      !	fnumaerosv = fraction of total number in the bin (--)
      !	fvolaerosv = fraction of total volume in the bin (--)
      !
      anumsum = 0._r8
      avolsum = 0._r8
      do i = 1, na
          if (aero_props%model_is('fanci')) then

              ! Assume just one mean value and one bin:
              a = ag0 ! radius of bin center.
              ynumaerosv(i) = 1._r8 !all in one
          else
            x = xlo + (i-1)*dx
            a = exp( x )
            dum = (x - xg0)/sx
            ynumaerosv(i) = exp( -0.5_r8*dum*dum )
          end if
         aaerosv(i) = a

         yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
         anumsum = anumsum + ynumaerosv(i)
         avolsum = avolsum + yvolaerosv(i)
      end do

      do i = 1, na
         ynumaerosv(i) = ynumaerosv(i)/anumsum
         yvolaerosv(i) = yvolaerosv(i)/avolsum
      end do
      !
      !   compute scavenging
      !
      scavsumnum = 0._r8
      scavsumvol = 0._r8
      !
      !   outer loop for rain drop radius
      !
      jr_loop: do jr = 1, nr

         r = rrainsv(jr)
         vfall = vfallrainsv(jr)

         reynolds = r * vfall / airkinvisc
         sqrtreynolds = sqrt( reynolds )

         !
         !   inner loop for aerosol particle radius
         !
         scavsumnumbb = 0._r8
         scavsumvolbb = 0._r8

         ja_loop: do ja = 1, na

            a = aaerosv(ja)

            chi = a/r

            dum = freepath/a
            dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
            taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)

            aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
            aerodiffus = boltz_cgs*temp*taurelax/aeromass

            schmidt = airkinvisc/aerodiffus
            stokes = vfall*taurelax/r

            ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
                 (reynolds*schmidt)

            dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
                 (1._r8 + xmuwaterair/sqrtreynolds)
            eintercept = 4._r8*chi*(chi + dum)

            dum = log( 1._r8 + reynolds )
            sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
            eimpact = 0._r8
            if (stokes > sstar) then
               dum = stokes - sstar
               eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
            end if

            etotal = ebrown + eintercept + eimpact
            etotal = min( etotal, 1.0_r8 )

            rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

            scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
            scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

         enddo ja_loop

         scavsumnum = scavsumnum + scavsumnumbb
         scavsumvol = scavsumvol + scavsumvolbb

      enddo jr_loop

      scavratenum = scavsumnum*3600._r8
      scavratevol = scavsumvol*3600._r8

    end subroutine calc_1_impact_rate

  end subroutine init_bcscavcoef

end module aero_wetdep_cam
