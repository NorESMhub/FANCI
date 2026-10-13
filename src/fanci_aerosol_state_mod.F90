module fanci_aerosol_state_mod
! TODO: make cloud borne tracers
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_spfn_mod, only: erf => shr_spfn_erf
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use fanci_aerosol_properties_mod, only: fanci_aerosol_properties
  use physconst,  only: rhoh2o, mwh2o

  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use ppgrid,         only: pver, pcols
  use string_utils,   only: int2str
  use ppgrid,         only: begchunk, endchunk
  use constituents, only: cnst_get_ind, qmin
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

  implicit none

  private

  public :: fanci_aerosol_state
  public :: aero_state_ptr


  integer, allocatable :: nr_copies(:) ! internal count of state object copies
  real(r8), parameter   :: dryvol_min = 1e-26_r8   ! m3_aer/kg_air; range_dry_volume at/below this => no aerosol

  type aero_state_ptr
    type(fanci_aerosol_state), pointer :: ptr => null()
  end type aero_state_ptr

  type aerosol_range_state ! one instance per range
     character(len=16), allocatable :: range_name(:)
     integer, allocatable  :: transport_ndx(:)
     integer, allocatable  :: spec_ndx(:)             ! same length as transport_ndx or third dimension of mmr(:,:,range_nspecies), indices corresponding to species properties object
     real(r8), allocatable :: dry_density(:,:)        ! density of the species mixture in a range without water, ncol, pver
     real(r8), allocatable :: hygroscopicity(:,:)     ! hygroscopicity of the species mixture
     real(r8), pointer     :: mmr_cw(:, :, :) => null()         ! (kg/kg) (ncol, pver, range_nspecies) cloud borne stuff, not transported
     real(r8), allocatable :: massfrac(:,:,:)         ! mass fraction of each species
     ! ...

  end type aerosol_range_state

  type, extends(aerosol_state) :: fanci_aerosol_state
     !private

     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     type(fanci_aerosol_properties), pointer :: sec_aero_props => null()
     real(r8), pointer :: bin_numconc_cw(:,:,:) => null()
     real(r8), pointer :: wet_radius(:,:,:) => null()           ! wet radius at bin center -> wet_radius*2 = dgnumwet
     real(r8), pointer :: qaerwat(:,:,:) => null()              ! aerosol water concentration (g/g)
     integer, allocatable  :: num_transport_ndx(:)
     type(aerosol_range_state), allocatable :: aero_range_state(:)
     integer :: ncol
     integer :: lchnk

   contains

     procedure :: get_transported               ! done
     procedure :: set_transported               ! done
     procedure :: ambient_total_bin_mmr         ! done
     procedure :: get_ambient_mmr_0list         ! done
     procedure :: get_ambient_mmr_rlist   ! TODO
     procedure :: get_cldbrne_mmr               ! done
     procedure :: get_ambient_num               ! done
     procedure :: get_cldbrne_num               ! done
     procedure :: get_states                    ! done
     procedure :: icenuc_size_wght_arr    ! TODO
     procedure :: icenuc_size_wght_val    ! TODO
     procedure :: icenuc_type_wght        ! TODO
     procedure :: update_bin              ! TODO
     procedure :: hetfrz_size_wght        ! TODO
     procedure :: hygroscopicity                ! mostly done, implement list_ndx
     procedure :: water_uptake            ! TODO
     procedure :: dry_volume                    ! done
     procedure :: wet_volume                    ! done (copied from modal)
     procedure :: water_volume                  ! done (copied from modal)
     procedure :: wet_diameter                  ! done (wet_radius stored internally)
     procedure :: convcld_actfrac         ! TODO
     procedure :: wgtpct                  ! TODO
     procedure :: bin_dry_density               ! done
     procedure :: update_range                  ! done
     procedure :: bin_species_mmr               !done, duplicate of get_ambient_mmr?
     final :: destructor

  end type fanci_aerosol_state

  interface fanci_aerosol_state
     procedure :: constructor
  end interface fanci_aerosol_state

  real(r8), parameter :: rho_aer_fallback = 1000._r8   ! kg/m3; only used where a range holds no aerosol
  real(r8), parameter :: rh2odens = 1._r8/rhoh2o
  type(aero_state_ptr), allocatable :: master_aero_state(:)

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state, pbuf, copy) result(newobj)
    use mo_tracname, only: solsym
    use constituents, only: cnst_get_ind
    use string_utils, only: int2str

    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(fanci_aerosol_state), pointer :: newobj
    type(fanci_aerosol_properties), target :: aero_props
    logical, intent(in), optional :: copy
    logical :: make_copy
    integer :: ierr, irange, solsym_ndx, ispec, ubar_ndx, ibin, ispecprops, lchnk
    character(len=:), allocatable :: num_name
    character(len=10)             :: speciesname, speciesname_props
    logical :: solsym_found
    integer :: icol, iver
    character(len=*), parameter :: subname = 'constructor'

    make_copy = .false.

    if ( present(copy) ) then
        make_copy = copy
    end if

    if ( .not. allocated(master_aero_state) ) then
        allocate(master_aero_state(begchunk:endchunk))
    end if

    if ( .not. allocated(nr_copies) ) then
        allocate(nr_copies(begchunk:endchunk))
        nr_copies = 0
    end if

    lchnk = state%lchnk
    if ( associated(master_aero_state(lchnk)%ptr) .and. .not. make_copy ) then
        newobj => master_aero_state(lchnk)%ptr
    else

        allocate(newobj,stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        newobj%state => state
        newobj%pbuf => pbuf
        newobj%sec_aero_props => fanci_aerosol_properties()

        newobj%ncol = state%ncol
        newobj%lchnk = state%lchnk

        allocate(newobj%bin_numconc_cw(newobj%ncol, pver, newobj%sec_aero_props%nbins()), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        allocate(newobj%wet_radius(newobj%ncol, pver, newobj%sec_aero_props%nbins()), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        allocate(newobj%qaerwat(newobj%ncol, pver, newobj%sec_aero_props%nbins()), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if
        ! temporary fix: no hygroscopic growth yet -> wet size = dry size.
        ! use bin_centers, alread in m
        ! TODO: use real wet radius!
        do icol = 1, newobj%ncol
            do iver = 1, pver
                newobj%wet_radius(icol, iver,:) = newobj%sec_aero_props%bin_centers(newobj%sec_aero_props%nbins())
            end do
        end do
        newobj%qaerwat(:, :,:) = 0._r8

        allocate(newobj%aero_range_state(newobj%sec_aero_props%nranges()), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        allocate(newobj%num_transport_ndx(newobj%sec_aero_props%nbins()), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        newobj%bin_numconc_cw = 0._r8
        newobj%num_transport_ndx = 0

        do irange = 1, newobj%sec_aero_props%nranges()
            allocate(newobj%aero_range_state(irange)%mmr_cw(newobj%ncol, pver, newobj%sec_aero_props%range_nspecies(irange) ), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%massfrac(newobj%ncol, pver, newobj%sec_aero_props%range_nspecies(irange) ), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%dry_density(newobj%ncol, pver), stat=ierr)
            if( ierr /= 0 )then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%hygroscopicity(newobj%ncol, pver), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%range_name(newobj%sec_aero_props%range_nspecies(irange)), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%transport_ndx(newobj%sec_aero_props%range_nspecies(irange)), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if
            allocate(newobj%aero_range_state(irange)%spec_ndx(newobj%sec_aero_props%range_nspecies(irange)), stat=ierr)
            if( ierr /= 0 ) then
                nullify(newobj)
                return
            end if

            newobj%aero_range_state(irange)%dry_density = 0._r8
            newobj%aero_range_state(irange)%hygroscopicity = 0._r8
            newobj%aero_range_state(irange)%mmr_cw = 0._r8
            newobj%aero_range_state(irange)%massfrac = 0._r8
            newobj%aero_range_state(irange)%range_name = ''
            newobj%aero_range_state(irange)%transport_ndx = 0
            newobj%aero_range_state(irange)%spec_ndx = 0

            ispec=0
            do solsym_ndx = 1, size(solsym)

                ubar_ndx = index(solsym(solsym_ndx), '_R'//int2str(irange))

                if ( ubar_ndx > 0 ) then
                    ispec = ispec + 1

                    if ( ispec > newobj%sec_aero_props%range_nspecies(irange) ) then
                        call endrun(subname//':: ERROR : number of species is larger than number of species in range')
                    end if

                    newobj%aero_range_state(irange)%range_name(ispec) = solsym(solsym_ndx)
                    call cnst_get_ind(solsym(solsym_ndx), newobj%aero_range_state(irange)%transport_ndx(ispec), abort=.false.)
                    if ( newobj%aero_range_state(irange)%transport_ndx(ispec) < 0 ) then
                        call endrun(subname//":: ERROR: transport array index for"//trim(solsym(solsym_ndx))//" not found")
                    end if

                    ! add index to connect to the species objects in aero_props
                    speciesname = solsym(solsym_ndx)

                    do ispecprops = 1, newobj%sec_aero_props%range_nspecies(irange)
                        ! translate irange to ibin
                        ibin = newobj%sec_aero_props%range_bounds(irange, 1)
                        call newobj%sec_aero_props%get(bin_ndx=ibin,species_ndx=ispecprops, specname=speciesname_props)
                        if ( speciesname(1:ubar_ndx-1) == trim( speciesname_props ) ) then ! before ubar_ndx -> species name
                            newobj%aero_range_state(irange)%spec_ndx(ispec) = ispecprops
                        end if
                    end do

                end if
            end do

        end do

        do ibin = 1, newobj%sec_aero_props%nbins()
            num_name = 'num_'//int2str(ibin)
            solsym_found = .false.
            do solsym_ndx = 1, size(solsym)
                if ( trim(solsym(solsym_ndx)) == trim(num_name) ) then
                    solsym_found = .true.
                    exit
                end if
            end do
            if ( .not. solsym_found ) then
                call endrun(subname//':: ERROR: bin '//trim(num_name)//' not found')
            end if
            call cnst_get_ind(num_name, newobj%num_transport_ndx(ibin), abort = .false.)
            if ( newobj%num_transport_ndx(ibin) < 0 ) then
                call endrun(subname//" :: ERROR: transport array index for "//trim(num_name)//' not found')
            end if
        end do


        if (.not. associated(master_aero_state(lchnk)%ptr)) then
            master_aero_state(lchnk)%ptr => newobj
            if (nr_copies(lchnk) > 2) then
 !               write(iulog,*) 'constructor: fanci_aerosol_state_mod: number of copies = ', nr_copies(lchnk)
                call endrun(subname//':: ERROR: master_aero_state already nr_copies > 2')
            end if
        end if

        nr_copies(lchnk) = nr_copies(lchnk) + 1

    end if


  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(fanci_aerosol_state), intent(inout) :: self

    character(len=*), parameter :: subname = 'destructor'

    nullify(self%state)
    nullify(self%pbuf)
    if (allocated(self%aero_range_state)) then
        deallocate(self%aero_range_state)
    end if
    if (associated(self%bin_numconc_cw)) then
        deallocate(self%bin_numconc_cw)
        nullify(self%bin_numconc_cw)
    end if
    if (associated(self%wet_radius)) then
        deallocate(self%wet_radius)
        nullify(self%wet_radius)
    end if
    if (associated(self%qaerwat)) then
        deallocate(self%qaerwat)
        nullify(self%qaerwat)
    end if

    if (masterproc) then
        write(iulog,*) 'destructor: fanci_aerosol_state_mod: number of copies = ', nr_copies(self%lchnk)
    end if

    nr_copies(self%lchnk) = nr_copies(self%lchnk) - 1

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol
  ! This updates the transported aerosol constituent array to match the aerosol model state.
  ! constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(fanci_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    integer              :: irange, ispec, ibin

    character(len=*), parameter :: subname = 'set_transported'

    ! update the range properties
    do irange = 1, self%sec_aero_props%nranges()
        call self%update_range(irange=irange, ncol=self%ncol)
    end do
    ! TODO add something like:
    !do ibin=1, self%sec_aero_props%nbins()
    !   call self%water_uptake(self%sec_aero_props, 0, ibin, self%ncol, pver, gnumwet, qaerwat)
    !   self%wet_radius(:,:,ibin) = 0.5_r8 * dgnumwet
    !   self%qaerwat(:,:,ibin)    = qaerwat
    ! end do
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(fanci_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    integer               :: irange, ispec, ibin

    character(len=*), parameter :: subname = 'get_transported'

    return

  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    integer :: ncol, irange

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

    character(len=*), parameter :: subname = 'ambient_total_bin_mmr'

    ncol = self%state%ncol
    irange = self%sec_aero_props%bins2ranges(bin_ndx)

    ! (kg_tot_aerosol_in_bin / kg_air)
    mmr_tot = self%state%q(col_ndx,lyr_ndx,self%num_transport_ndx(bin_ndx)) &
                * self%sec_aero_props%particle_volume(bin_ndx) &
                * self%aero_range_state(irange)%dry_density(col_ndx, lyr_ndx)

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_0list(self, species_ndx, bin_ndx, mmr)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer   :: mmr(:,:)     ! mass mixing ratios (ncol,nlev)
    integer             :: range

    character(len=*), parameter :: subname = 'get_ambient_mmr_0list'

    range = self%sec_aero_props%bins2ranges(bin_ndx)
    mmr => self%state%q(:self%ncol,:,self%aero_range_state(range)%transport_ndx(species_ndx))

  end subroutine get_ambient_mmr_0list

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given radiation diagnostics
  ! list index, species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_rlist(self, list_ndx, species_ndx, bin_ndx, mmr)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx     ! rad climate list index
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    character(len=*), parameter :: subname = 'get_ambient_mmr_rlist'

    call endrun(subname//' is not yet implemented')

  end subroutine get_ambient_mmr_rlist

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer   :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)
    integer             :: range

    character(len=*), parameter :: subname = 'get_cldbrne_mmr'

    range = self%sec_aero_props%bins2ranges(bin_ndx)
    mmr => self%aero_range_state(range)%mmr_cw(:,:,species_ndx)

  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities

    character(len=*), parameter :: subname = 'get_ambient_num'

    num => self%state%q(:self%ncol,:,self%num_transport_ndx(bin_ndx))

  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    character(len=*), parameter :: subname = 'get_cldbrne_num'

    num => self%bin_numconc_cw(:,:,bin_ndx)

  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: ibin,ispec, indx

    character(len=*), parameter :: subname = 'get_states'

    do ibin = 1, aero_props%nbins()
        indx = aero_props%indexer(ibin, 0)
        call self%get_ambient_num(ibin, raer(indx)%fld)
        call self%get_cldbrne_num(ibin, qqcw(indx)%fld)

        if ( self%sec_aero_props%nmasses(ibin) /= 0 ) then
            do ispec = 1, self%sec_aero_props%nmasses(ibin)
                indx = aero_props%indexer(ibin, ispec)
                call self%get_ambient_mmr(ispec,ibin, raer(indx)%fld)
                call self%get_cldbrne_mmr(ispec,ibin, qqcw(indx)%fld)
            end do
        end if
    end do

  end subroutine get_states

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_arr(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=*), parameter :: subname = 'icenuc_size_wght_arr'

    call endrun(subname//' is not yet implemented')
! ??
  end subroutine icenuc_size_wght_arr

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and vertical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    character(len=*), parameter :: subname = 'icenuc_size_wght_val'

    call endrun(subname//' is not yet implemented')

  end subroutine icenuc_size_wght_val

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)

    use aerosol_properties_mod, only: aerosol_properties

    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    character(len=*), parameter :: subname = 'icenuc_type_wght'

    call endrun(subname//' is not yet implemented')

  end subroutine icenuc_type_wght

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
    real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
    integer, intent(in) :: tnd_ndx                ! tendency index
    real(r8),intent(in) :: dtime                  ! time step size (sec)
    real(r8),intent(inout) :: tend(:,:,:)         ! tendency

    character(len=*), parameter :: subname = 'update_bin'

    call endrun(subname//' is not yet implemented')

    ! bin_num(bin_ndx) = bin_num + tendency
  end subroutine update_bin

  !------------------------------------------------------------------------------
  ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
  ! as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels

    real(r8) :: wght(ncol,nlev)                 !

    character(len=*), parameter :: subname = 'hetfrz_size_wght'

    call endrun(subname//' is not yet implemented')

  end function hetfrz_size_wght

  !------------------------------------------------------------------------------
  ! returns hygroscopicity for a given radiation diagnostic list number and
  ! bin number
  !------------------------------------------------------------------------------
  subroutine hygroscopicity(self, list_ndx, bin_ndx, kappa)
    class(fanci_aerosol_state), intent(in) :: self
! TODO: what is list_ndx?
    integer, intent(in) :: list_ndx        ! rad climate list number
    integer, intent(in) :: bin_ndx         ! bin number
    real(r8), intent(out) :: kappa(:,:)                 !
    integer :: irange

    character(len=*), parameter :: subname = 'hygroscopicity'

    irange = self%sec_aero_props%bins2ranges(bin_ndx)
    kappa = self%aero_range_state(irange)%hygroscopicity

  end subroutine hygroscopicity

  !------------------------------------------------------------------------------
  ! returns aerosol wet diameter and aerosol water concentration for a given
  ! radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  subroutine water_uptake(self, aero_props, list_ndx, bin_ndx, ncol, nlev, dgnumwet, qaerwat)

    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: list_ndx             ! rad climate/diags list number
    integer, intent(in) :: bin_ndx              ! bin number
    integer, intent(in) :: ncol                 ! number of columns
    integer, intent(in) :: nlev                 ! number of levels
    real(r8),intent(out) :: dgnumwet(ncol,nlev) ! aerosol wet diameter (m)
    real(r8),intent(out) :: qaerwat(ncol,nlev)  ! aerosol water concentration (g/g)

    character(len=*), parameter :: subname = 'water_uptake'
! TODO: get sulfate to work first!!
    call endrun(subname//' is not yet implemented')

  end subroutine water_uptake

  !------------------------------------------------------------------------------
  ! aerosol dry volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function dry_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    character(len=*), parameter :: subname = 'dry_volume'

    vol = self%state%q(:ncol,:nlev,self%num_transport_ndx(bin_ndx)) &
        * self%sec_aero_props%particle_volume(bin_ndx)

  end function dry_volume

  !------------------------------------------------------------------------------
  ! aerosol wet volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function wet_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dryvol(ncol,nlev)
    real(r8) :: watervol(ncol,nlev)

    character(len=*), parameter :: subname = 'wet_volume'

    dryvol = self%dry_volume(aero_props, list_ndx, bin_ndx, ncol, nlev)
    watervol = self%water_volume(aero_props, list_ndx, bin_ndx, ncol, nlev)

    vol = watervol + dryvol

  end function wet_volume

  !------------------------------------------------------------------------------
  ! aerosol water volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function water_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(fanci_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dgnumwet(ncol,nlev)
    real(r8) :: qaerwat(ncol,nlev)

    character(len=*), parameter :: subname = 'water_volume'

    call self%water_uptake(aero_props, list_ndx, bin_ndx, ncol, nlev, dgnumwet, qaerwat)

    vol(:ncol,:nlev) = qaerwat(:ncol,:nlev)*rh2odens
    where (vol<0._r8)
       vol = 0._r8
    end where

  end function water_volume

  !------------------------------------------------------------------------------
  ! aerosol wet diameter
  !------------------------------------------------------------------------------
  function wet_diameter(self, bin_ndx, ncol, nlev) result(diam)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: diam(ncol,nlev)

    character(len=*), parameter :: subname = 'wet_diameter'

    diam(:ncol,:nlev) = 2._r8 * self%wet_radius(:ncol, :nlev, bin_ndx)
  end function wet_diameter

  !------------------------------------------------------------------------------
  ! prescribed aerosol activation fraction for convective cloud
  !------------------------------------------------------------------------------
  function convcld_actfrac(self, ibin, ispc, ncol, nlev) result(frac)

    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: ispc   ! species index
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of vertical levels

    real(r8) :: frac(ncol,nlev)

    character(len=*), parameter :: subname = 'convcld_actfrac'

    !TODO make these into namelist options and make it more nuanced...
    ! Ideally we would make a fraction per range maybe? Or I don't see why it could not call the activation code for a
    ! fixed (and high) supersaturation? Maybe expensive?

    frac=0.8

  end function convcld_actfrac

  !------------------------------------------------------------------------------
  ! aerosol weight precent of H2SO4/H2O solution
  !------------------------------------------------------------------------------
  function wgtpct(self, ncol, nlev) result(wtp)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) ::  ncol, nlev
    real(r8) :: wtp(ncol,nlev)  ! weight precent of H2SO4/H2O solution for given icol, ilev
    character(len=*), parameter :: subname = 'wgtpct'

    call endrun(subname//' is not yet implemented')

  end function wgtpct

  function bin_dry_density(self, bin_ndx, ncol) result(ddens)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in)   :: bin_ndx
    integer, intent(in)   :: ncol                 ! number of columns
    integer               :: irange
    real(r8)              :: ddens(ncol,pver)

    character(len=*), parameter :: subname = 'dry_density'

    irange = self%sec_aero_props%bins2ranges(bin_ndx)
    ddens = self%aero_range_state(irange)%dry_density

  end function bin_dry_density

  subroutine update_range(self, irange, ncol, mmr_tend)
    class(fanci_aerosol_state), intent(inout) :: self
    real(r8), optional, intent(in) :: mmr_tend(:,:,:) ! shape aero_props%range_nspecies (ncol, pver, range_nspecies) -> one array for one range
    integer, intent(in)  :: irange
    integer, intent(in)  :: ncol                 ! number of columns
    integer              :: ispec, ibin, ispecprop
    real(r8)             :: range_dry_volume(ncol, pver)
    real(r8)             :: range_total_mmr(ncol, pver)

    range_dry_volume = 0._r8
    range_total_mmr = 0._r8

    ! update mmr of each component if mmr tendency has been passed as an argument
    ! else: update other state variables with mmr from before
    if ( present(mmr_tend) ) then
        self%state%q(:ncol,:pver, self%aero_range_state(irange)%transport_ndx(:)) = &
        self%state%q(:ncol,:pver, self%aero_range_state(irange)%transport_ndx(:)) + mmr_tend(:,:,:)
! TODO: fix this, use ptend!!?
    end if
    ! reset density and hygroscopicity
    self%aero_range_state(irange)%dry_density = 0._r8
    self%aero_range_state(irange)%hygroscopicity = 0._r8
    self%aero_range_state(irange)%massfrac = 0._r8

    range_total_mmr = sum(self%state%q(:ncol,:pver, self%aero_range_state(irange)%transport_ndx(:)), dim=3) ! sum over species (kg/kg)

    ! calculate mass fractions
    do ispec = 1, self%sec_aero_props%range_nspecies(irange)
        where ( range_total_mmr /= 0._r8 )
            self%aero_range_state(irange)%massfrac(:,:,ispec) = self% state%q(:ncol,:pver, self%aero_range_state(irange)%transport_ndx(ispec)) &
                                                                 / range_total_mmr
        end where
    end do

    ! volume of all aerosol /kg_air in a range (m3_aer/kg_air)
    do ibin = self%sec_aero_props%range_bounds(irange, 1), self%sec_aero_props%range_bounds(irange, 2)
        range_dry_volume = range_dry_volume + self%dry_volume(self%sec_aero_props, 1, ibin, self%state%ncol, pver) !TODO: change input parameters
    end do

    !where (range_dry_volume /= 0._r8)
    where (range_dry_volume > dryvol_min)
        self%aero_range_state(irange)%dry_density = range_total_mmr / range_dry_volume
    elsewhere
        self%aero_range_state(irange)%dry_density = rho_aer_fallback  ! fall back value, set at top.
    end where
    !self%aero_range_state(irange)%dry_density = range_total_mmr/range_dry_volume
    if ( self%sec_aero_props%range_nspecies(irange) /= 0 ) then
        do ispec = 1,self%sec_aero_props%range_nspecies(irange)
            ispecprop = self%aero_range_state(irange)%spec_ndx(ispec)
! TODO: source, total hygroscopicity parameter kappa_tot = SUM_OVER_ALL_SPECIES(volume_i/volume_tot * kappa_i)
            where(range_dry_volume /= dryvol_min)
            self%aero_range_state(irange)%hygroscopicity = self%aero_range_state(irange)%hygroscopicity &
                + self%state%q(:ncol,:pver, self%aero_range_state(irange)%transport_ndx(ispec)) &
                / self%sec_aero_props%density(ispecprop) &
                / range_dry_volume * self%sec_aero_props%kappa(ispecprop) ! TODO: probably not the least ugly way to do this
            end where
        end do
    end if
  end subroutine update_range

  function bin_species_mmr(self, bin_ndx, species_ndx, col_ndx, lyr_ndx) result(spec_mmr_in_bin)
    class(fanci_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx
    integer, intent(in) :: species_ndx
    integer, intent(in) :: col_ndx, lyr_ndx
    integer             :: irange
    real(r8)            :: spec_mmr_in_bin

    irange = self%sec_aero_props%bins2ranges(bin_ndx)
    spec_mmr_in_bin = self%aero_range_state(irange)%massfrac(col_ndx, lyr_ndx, species_ndx) &
                        * self%ambient_total_bin_mmr(self%sec_aero_props, bin_ndx, col_ndx, lyr_ndx)

  end function bin_species_mmr

end module fanci_aerosol_state_mod
