module fanci_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len

  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
   use string_utils, only: int2str

  implicit none

  private

  public :: fanci_aerosol_properties

  type aerosol_species_properties
     character(len=:), allocatable :: specname        ! e.g. DU
     character(len=:), allocatable :: spectype        ! e.g. 'dust'
     integer                 :: nbin            ! number of bins containing species
     integer                 :: nrange          ! nr of ranges containing species
     integer, allocatable    :: range_ndx(:)    ! indices of ranges containing species
     integer, allocatable    :: bin_ndx(:)      ! bin indices containing species
     real(r8)                :: density
     real(r8)                :: molecular_weight ! kg/kmol
     real(r8)                :: kappa
     logical                 :: mixed           ! true if internally mixed
     !TODO: add hydrophilic/phobic
     character(len=10), allocatable :: tracernames(:) ! e.g. DU_R3
     end type aerosol_species_properties

  type, extends(aerosol_properties) :: fanci_aerosol_properties
     private
     integer               :: nranges_ = 0
     integer               :: nspecies_tot_ = 0
     integer, allocatable  :: range_nspecies_(:)
     real(r8), allocatable :: bin_centers_(:) ! radii at bin center (m)
     real(r8), allocatable :: bin_bounds_(:,:)! radii at bin bounds (m)
     integer, allocatable  :: range_bounds_(:,:) ! index of bins at range bounds
     integer, allocatable  :: bins2ranges_(:) ! range index for each bin
     real(r8), allocatable :: particle_volume_(:) ! volume of a single particle in a bin (center radius) in m3
     type(aerosol_species_properties), allocatable :: aer_spec_prop(:)

   contains
     procedure :: number_transported
     procedure :: get                       ! partially implemented
     procedure :: amcube                    ! copied from CARMA
     procedure :: density                   ! done
     procedure :: actfracs                  ! from CARMA -> if we use effective critical supersat, this needs to be changed
     procedure :: num_names                 ! done
     procedure :: mmr_names                 ! done
     procedure :: amb_num_name          ! TODO
     procedure :: amb_mmr_name          ! TODO
     procedure :: species_type          ! TODO
     procedure :: icenuc_updates_num    ! TODO
     procedure :: icenuc_updates_mmr    ! TODO
     procedure :: apply_number_limits   ! ??TODO, maybe
     procedure :: hetfrz_species        ! TODO
     procedure :: optics_params         ! TODO
     procedure :: nbins_rlist           ! TODO
     procedure :: nspecies_tot              ! done
     procedure :: range_nspecies            ! done
     procedure :: nranges                   ! done
     procedure :: bin_centers               ! done
     procedure :: bin_bounds                ! done
     procedure :: particle_volume           ! done
     procedure :: spec_range_ndx            ! done
     procedure :: spec_bin_ndx              ! done
     procedure :: spec_nrange               ! done
     procedure :: spec_nbin                 ! done
     procedure :: spec_tracernames          ! done
     procedure :: bins2ranges               ! done
     procedure :: nspecies_per_bin_rlist ! TODO
     procedure :: alogsig_rlist          ! TODO
     procedure :: soluble               ! TODO
     procedure :: min_mass_mean_rad     ! TODO
     procedure :: bin_name              ! TODO
     procedure :: scav_diam             ! TODO
     procedure :: resuspension_resize   ! TODO
     procedure :: rebin_bulk_fluxes         ! done
     procedure :: hydrophilic           ! TODO: currently simple approximation
     procedure :: model_is                  ! done
     procedure :: is_active                  ! done
     procedure :: range_bounds              ! done
     procedure :: kappa                     ! done
     procedure :: molecular_weight          ! done
     procedure :: specname                  ! done
     procedure :: spectype                  ! done
     procedure :: spec_bin_q_ndx  ! done
     procedure :: spec_mmr_q_ndx  ! done
     final :: destructor
  end type fanci_aerosol_properties

  interface fanci_aerosol_properties
     procedure :: constructor
  end interface fanci_aerosol_properties

  logical, parameter :: debug = .false.
  type(fanci_aerosol_properties), pointer :: prop_obj => null()


contains
  !------------------------------------------------------------------------------
  function constructor(nlfile) result(newobj)

    use mpi,               only: mpi_integer, mpi_real8, mpi_character, mpi_logical, MPI_SUCCESS
    use spmd_utils,        only: mstrid=>masterprocid, mpicom
    use string_utils,      only: int2str
    use namelist_utils,    only: find_group_name
    use infnan,            only: nan, assignment(=)

    type(fanci_aerosol_properties), pointer :: newobj

    character(len=*),optional, intent(in) :: nlfile
    integer                      :: ncnst_tot=0
    integer,allocatable          :: nspecies(:) ! nspecies per bin (given by base_object)
    integer,allocatable          :: nmasses(:)  ! nspecies for first bin in a range, 0 elsewhere
    real(r8),allocatable         :: alogsig(:) ! given by base obj
    real(r8),allocatable         :: f1(:) ! given by base obj: Abdul-Razzak 1998 eq 28
    real(r8),allocatable         :: f2(:) ! given by base obj: Abdul-Razzak 1998 eq 29
    real(r8),allocatable         :: bin_centers(:)
    real(r8),allocatable         :: bin_bounds(:,:)
    integer,allocatable          :: range_bounds(:,:)
    integer,allocatable          :: bins2ranges(:)
    integer                      :: bin_ndx(100)
    integer                      :: ierr, unitn, pos, lower, upper
    integer                      :: ind, ibin, irange, ispec

    ! namelist variables
    integer               :: fanci_nbins
    integer               :: fanci_nranges
    integer               :: fanci_nspecies_tot
    integer               :: fanci_nspecies(500) ! range_nspecies TODO: make allocatable!!
    integer, parameter    :: strlen=50

    character(len=strlen) :: fanci_bin_centers(500)
    character(len=strlen) :: fanci_bin_bounds(500)
    character(len=strlen) :: fanci_range_bounds(500)

    ! namelist aerosol species variables
    type(aerosol_species_properties), allocatable :: fanci_species_properties(:)
    character(len=10)                     :: fanci_aerosol_name
    character(len=10)                     :: fanci_aerosol_type
    character(len=10)                     :: fanci_aerosol_range
    real(r8)                              :: fanci_aerosol_density
    real(r8)                              :: fanci_aerosol_weight
    real(r8)                              :: fanci_aerosol_kappa
    logical                               :: fanci_aerosol_mixed

    character(len=aero_name_len) :: spectype

    character(len=*), parameter :: subname = 'constructor'

!==================================================================================================
! Read Namelists
!==================================================================================================

      ! Namelists (constructed in bin_config.py)
    namelist /fanci_properties_nl/ fanci_nspecies_tot, &
                                            fanci_nspecies, &
                                            fanci_nbins, &
                                            fanci_nranges, &
                                            fanci_bin_bounds, &
                                            fanci_bin_centers, &
                                            fanci_range_bounds

    namelist /fanci_properties_aerosol_nl/ fanci_aerosol_name, &
                                            fanci_aerosol_type, &
                                            fanci_aerosol_range, &
                                            fanci_aerosol_density, &
                                            fanci_aerosol_weight, &
                                            fanci_aerosol_mixed, &
                                            fanci_aerosol_kappa
    if ( present(nlfile) ) then

        if ( associated(prop_obj) ) then
            call endrun(subname//':: ERROR fanci_aerosol_properties has already been initialized')
        end if
        ! initialize variables
        fanci_nspecies_tot = 0
        fanci_nspecies = 0
        fanci_nbins = 0
        fanci_nranges = 0
        fanci_bin_centers = ''
        fanci_bin_bounds = ''
        fanci_range_bounds = ''

        ! read aerosol properties namelist
        if (masterproc) then
            open(newunit=unitn, file=trim(nlfile), status='old')
            call find_group_name(unitn, 'fanci_properties_nl', ierr)
            if ( ierr /= 0 ) then
                close(unitn)
                call endrun(subname//":: ERROR could not find group 'fanci_properties_nl'")
            end if

            read(unitn, fanci_properties_nl, iostat=ierr)
            if ( ierr /= 0 ) then
                close(unitn)
                call endrun(subname // ':: ERROR reading fanci_properties_nl namelist')
            end if
        end if

!==================================================================================================
! Broadcast fanci_properties
!==================================================================================================
        call MPI_Bcast(fanci_nspecies_tot, 1, mpi_integer, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_nspecies_tot'")
        end if

        call MPI_Bcast(fanci_nspecies, size(fanci_nspecies), mpi_integer, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_nspecies'")
        end if

        call MPI_Bcast(fanci_nbins, 1, mpi_integer, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_nbins'")
        end if

        call MPI_Bcast(fanci_nranges, 1, mpi_integer, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_nranges'")
        end if

        call MPI_Bcast(fanci_bin_bounds, strlen*size(fanci_bin_bounds), mpi_character, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_bin_bounds'")
        end if

        call MPI_Bcast(fanci_bin_centers, strlen*size(fanci_bin_centers), mpi_character, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_bin_centers")
        end if

        call MPI_Bcast(fanci_range_bounds, strlen*size(fanci_range_bounds), mpi_character, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_range_bounds")
        end if

!==================================================================================================
! Read and broadcast aerosol species properties namelist in sequence
!==================================================================================================
 !    allocate array for species objects
        allocate(fanci_species_properties(fanci_nspecies_tot), stat=ierr)
        if(ierr/=0) then
            if (masterproc) close(unitn)
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating fanci_species_properties")
        end if

        do ispec=1,fanci_nspecies_tot

            ! namelist variables
            fanci_aerosol_name = ''
            fanci_aerosol_type = ''
            fanci_aerosol_range = ''
            fanci_aerosol_density = 0.0_r8
            fanci_aerosol_weight = 0.0_r8
            fanci_aerosol_kappa = 0.0_r8
            fanci_aerosol_mixed = .false.
            lower = 0
            upper = 0

            ! read namelist
            if (masterproc) then
                call find_group_name(unitn, 'fanci_properties_aerosol_nl', ierr)
                if (ierr /= 0) then
                    close(unitn)
                    call endrun(subname//":: ERROR could not find group 'fanci_properties_aerosol_nl' for species number " &
                    //int2str(ispec)//'of '//int2str(fanci_nspecies_tot) )
                end if

                read(unitn, fanci_properties_aerosol_nl, iostat=ierr)
                if (ierr /= 0) then
                    close(unitn)
                    call endrun(subname // ":: ERROR reading 'fanci_properties_aerosol_nl' for species number " &
                    //int2str(ispec)//'of '//int2str(fanci_nspecies_tot) )
                end if

            end if

            ! broadcast
            call MPI_Bcast(fanci_aerosol_name, len(fanci_aerosol_name), mpi_character, mstrid, mpicom, ierr)
            if ( ierr /= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_name'")
            end if

            call MPI_Bcast(fanci_aerosol_type, len(fanci_aerosol_type), mpi_character, mstrid, mpicom, ierr)
            if ( ierr /= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_type'")
            end if

            call MPI_Bcast(fanci_aerosol_range, len(fanci_aerosol_range), mpi_character, mstrid, mpicom, ierr)
            if ( ierr /= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_range'")
            end if

            call MPI_Bcast(fanci_aerosol_density, 1, mpi_real8, mstrid, mpicom, ierr)
            if ( ierr/= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_density'")
            end if

            call MPI_Bcast(fanci_aerosol_weight, 1, mpi_real8, mstrid, mpicom, ierr)
            if ( ierr/= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_weight'")
            end if

            call MPI_Bcast(fanci_aerosol_kappa, 1, mpi_real8, mstrid, mpicom, ierr)
            if ( ierr/= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_kappa'")
            end if

            call MPI_Bcast(fanci_aerosol_mixed, 1, mpi_logical, mstrid, mpicom, ierr)
            if ( ierr/= MPI_SUCCESS ) then
                call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'fanci_aerosol_mixed'")
            end if

            ! initialize properties object variables
            fanci_species_properties(ispec)%specname = ''
            fanci_species_properties(ispec)%density = 0.0_r8
            fanci_species_properties(ispec)%molecular_weight = 0.0_r8
            fanci_species_properties(ispec)%kappa = 0.0_r8
            fanci_species_properties(ispec)%mixed = .false.
            fanci_species_properties(ispec)%nbin = 0
            fanci_species_properties(ispec)%nrange = 0

            pos = index(fanci_aerosol_range, ':')

            if ( pos == 0 ) then
                call endrun(subname//":: ERROR invalid format for 'fanci_aerosol_range'")
            end if

            read(fanci_aerosol_range(1:pos-1), *) lower
            read(fanci_aerosol_range(pos+1:), *) upper

            fanci_species_properties(ispec)%nrange = upper - lower + 1
            fanci_species_properties(ispec)%specname = fanci_aerosol_name
            fanci_species_properties(ispec)%spectype = fanci_aerosol_type
            fanci_species_properties(ispec)%density = fanci_aerosol_density
            fanci_species_properties(ispec)%molecular_weight = fanci_aerosol_weight
            fanci_species_properties(ispec)%kappa = fanci_aerosol_kappa
            fanci_species_properties(ispec)%mixed = fanci_aerosol_mixed

            ! allocate and initialize range_ndx, and tracernames
            allocate(fanci_species_properties(ispec)%range_ndx(fanci_species_properties(ispec)%nrange), stat=ierr)
            if(ierr/=0) then
                call endrun(subname// ": ERROR "//int2str(ierr)//" allocating fanci_species_properties range_ndx")
            end if
            fanci_species_properties(ispec)%range_ndx = 0

            do ind = 1, (upper - lower + 1)
                fanci_species_properties(ispec)%range_ndx(ind) = lower + ind - 1
            end do

        end do

        if (masterproc) then
            close(unitn)
        end if

!==================================================================================================
! Allocate local arrays
!==================================================================================================

        ! allocate with size nbins
        allocate( nspecies(fanci_nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating nspecies")
            return
        end if

        allocate( nmasses(fanci_nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating nmasses")
            return
        end if

        allocate( bins2ranges(fanci_nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating bins2ranges'")
            return
        end if

        allocate( bin_centers(fanci_nbins), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_centers'")
           return
        end if

        allocate( bin_bounds(fanci_nbins, 2), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_bounds'")
           return
        end if

        allocate( range_bounds(fanci_nranges, 2), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating range_bounds'")
           return
        end if

! TODO: what to do with these, since we are not using them currently?
        allocate( alogsig(fanci_nbins), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating alogsig'")
           return
        end if

        allocate( f1(fanci_nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f1'")
           return
        end if

        allocate( f2(fanci_nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f2'")
           return
        end if
!==================================================================================================
! parse bin info
!==================================================================================================
        ! TODO: figure out what to do with these
        ! alogsig -> alogsig(m) = log(sigmag(m))
        ! f1 -> f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m)) f28 -> if sigma = 1, then this is 0.5
        ! f2 -> f2(m) = 1._r8 + 0.25_r8*alogsig(m)               f29 -> if sigma = 1, then this is 1
        alogsig = 0._r8
        f1 = 0.5_r8
        f2 = 1._r8 !nan

        ! TODO: check indexer_ var in aerosol_properties mod. same as the indices from chemical pp?
        do ibin=1,fanci_nbins
            read(fanci_bin_centers(ibin),'(F10.5)') bin_centers(ibin)
            pos = index(fanci_bin_bounds(ibin), ':')
            read(fanci_bin_bounds(ibin)(1:pos-1), '(F10.5)') bin_bounds(ibin,1)
            read(fanci_bin_bounds(ibin)(pos+1:), '(F10.5)') bin_bounds(ibin,2)
        end do

        ! parse range bounds
        do irange=1,fanci_nranges
            pos = index(fanci_range_bounds(irange), ':')
            read(fanci_range_bounds(irange)(1:pos-1), '(I3)') range_bounds(irange,1)
            read(fanci_range_bounds(irange)(pos+1:), '(I3)') range_bounds(irange,2)
        end do

        ncnst_tot = fanci_nbins + sum(fanci_nspecies(:fanci_nranges)) ! nbins + sum(nspecies)

        nmasses = 0
        ! initialize bins2ranges array
        do irange=1,fanci_nranges
            do ibin = range_bounds(irange, 1), range_bounds(irange, 2)
                bins2ranges(ibin) = irange
                nspecies(ibin) = fanci_nspecies(irange)
            end do
            nmasses(range_bounds(irange,1)) = nspecies(range_bounds(irange, 1))
        end do

        ! fill species objects with bin info
        do ispec = 1, fanci_nspecies_tot
            ind = 0
            fanci_species_properties(ispec)%nbin = 0
            do ibin = 1, fanci_nbins
                if (bins2ranges(ibin) >= fanci_species_properties(ispec)%range_ndx(1) .and. &
                    bins2ranges(ibin) <= maxval(fanci_species_properties(ispec)%range_ndx)) then
                    ! nr of bins containing each species (e.g. dust_nbin in dust_model.F90)
                    fanci_species_properties(ispec)%nbin = fanci_species_properties(ispec)%nbin + 1
                    ind = ind+1
                    ! indices of bins containing each species
                    bin_ndx(ind) = ibin
                end if
            end do
! TODO: I moved the 6 lines below here, after species%nbin is initialized to make the bin_ndx allocatable but feels a bit messy..
            allocate(fanci_species_properties(ispec)%bin_ndx(fanci_species_properties(ispec)%nbin), stat=ierr)
            if(ierr/=0) then
                call endrun(subname// ": ERROR "//int2str(ierr)//" allocating fanci_species_properties bin_ndx")
            end if
            fanci_species_properties(ispec)%bin_ndx = bin_ndx(:fanci_species_properties(ispec)%nbin)

            allocate(fanci_species_properties(ispec)%tracernames(fanci_species_properties(ispec)%nrange), stat=ierr)
            if(ierr/=0) then
                call endrun(subname// ": ERROR "//int2str(ierr)//" allocating fanci_species_properties tracernames")
            end if
            fanci_species_properties(ispec)%tracernames = ''

            do irange = 1, fanci_species_properties(ispec)%nrange
                ! tracernames for each species - potentially check if consistent with cnst_get_ind?
                fanci_species_properties(ispec)%tracernames(irange) = &
                trim(fanci_species_properties(ispec)%specname)//'_R'//&
                trim(int2str(fanci_species_properties(ispec)%range_ndx(irange)))
            end do

        end do

!==================================================================================================
! Allocate and initialize newobj
!==================================================================================================

        allocate(newobj,stat=ierr)
        if ( ierr/=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj")
        end if

        allocate(newobj%bins2ranges_(fanci_nbins), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%bins2ranges")
        end if

        allocate(newobj%range_nspecies_(fanci_nranges), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%range_nspecies")
        end if

        allocate(newobj%bin_centers_(fanci_nbins), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%bin_centers")
        end if

        allocate(newobj%bin_bounds_(fanci_nbins, 2), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%bin_bounds")
        end if

        allocate(newobj%range_bounds_(fanci_nranges, 2), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%range_bounds")
        end if

        allocate(newobj%aer_spec_prop(fanci_nspecies_tot), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%aer_spec_prop")
        end if

        allocate(newobj%particle_volume_(fanci_nbins), stat=ierr)
        if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%particle_volume")
        end if

        do ispec = 1, fanci_nspecies_tot
            allocate(newobj%aer_spec_prop(ispec)%bin_ndx(fanci_species_properties(ispec)%nbin), stat=ierr)
            if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%aer_spec_prop%bin_ndx")
            end if
            allocate(newobj%aer_spec_prop(ispec)%range_ndx(fanci_species_properties(ispec)%nrange), stat=ierr)
            if( ierr /=0 ) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating newobj%aer_spec_prop%range_ndx")
            end if
        end do

        newobj%bins2ranges_ = bins2ranges
        newobj%nranges_ = fanci_nranges
        newobj%nspecies_tot_ = fanci_nspecies_tot
        newobj%range_nspecies_ = fanci_nspecies(:fanci_nranges)
        newobj%bin_centers_ = bin_centers(:fanci_nbins) * 1.e-9_r8             ! nm to m
        newobj%bin_bounds_ = bin_bounds(:fanci_nbins, :) * 1.e-9_r8            ! nm to m
        newobj%range_bounds_ = range_bounds(:fanci_nranges, :)
        newobj%aer_spec_prop = fanci_species_properties(:fanci_nspecies_tot)
        newobj%particle_volume_ = 4._r8/3._r8*pi*(newobj%bin_centers_**3)

        ! deallocate local variables
        if (allocated(bin_centers)) deallocate(bin_centers)
        if (allocated(bin_bounds)) deallocate(bin_bounds)
        if (allocated(range_bounds)) deallocate(range_bounds)
        if (allocated(bins2ranges)) deallocate(bins2ranges)
        if (allocated(fanci_species_properties)) deallocate(fanci_species_properties)

        call newobj%initialize(fanci_nbins, ncnst_tot, nspecies, nmasses, alogsig, f1, f2, ierr)
!==================================================================================================
! Report
!==================================================================================================
!    call MPI_Barrier(mpicom, ierr)
        if (masterproc) then
            write(iulog,*) 'fanci aerosol properties: '
            write(iulog,*) 'nbins = ', newobj%nbins()
            write(iulog,*) 'nranges = ', newobj%nranges_
            write(iulog,*) 'ncnst_tot = ', newobj%ncnst_tot()
            write(iulog,*) 'nspecies_tot = ', newobj%nspecies_tot_
            write(iulog,*) 'nspecies = ', newobj%nspecies()

            do irange=1,fanci_nranges,5
                write(iulog,*) 'range_nspecies = ', newobj%range_nspecies_(irange:min(irange+4, fanci_nranges))
            end do

            do ibin=1,fanci_nbins,5
                write(iulog,*) 'bin_centers in nm = ', newobj%bin_centers_(ibin:min(ibin+4, fanci_nbins))*1e9
            end do

            do ibin=1,fanci_nbins,5
                write(iulog,*) 'bins2ranges = ',newobj%bins2ranges_(ibin:min(ibin+4, fanci_nbins))
            end do

            do ibin=1,fanci_nbins,5
                write(iulog,*) 'particle_volume in m3 = ',newobj%particle_volume_(ibin:min(ibin+4, fanci_nbins))
            end do

            write(iulog,*) 'Parameters for activation calculation are uniform across all bins for the fanci model: '
            write(iulog,*) 'f1 for abdul-razzak & gahn = ', f1(1) ! they are the same everywhere for fanci model
            write(iulog,*) 'f2 for abdul-razzak & gahn = ', f2(1)
            write(iulog,*) 'alogsig is = ', alogsig(1)

            do ibin=1,fanci_nbins
                write(iulog,*) 'bin_bounds in nm = ', newobj%bin_bounds_(ibin,1)*1e9, &
                                        ' : ', newobj%bin_bounds_(ibin,2)*1e9
            end do
            do irange=1,fanci_nranges !TODO FIX format
                write(iulog,*) 'range_bounds = ', newobj%range_bounds_(irange,1), &
                                        ' : ', newobj%range_bounds_(irange,2)
            end do

            do ind = 1,fanci_nspecies_tot
                ! TODO (low priority): fix the format :D
                write(iulog ,*) 'fanci aerosol species properties: '
                write(iulog ,*) 'species name = ', newobj%aer_spec_prop(ind)%specname
                write(iulog ,*) 'species type = ', newobj%aer_spec_prop(ind)%spectype
                write(iulog ,*) 'range indices = ', newobj%aer_spec_prop(ind)%range_ndx(1), ' : ', newobj%aer_spec_prop(ind)%range_ndx(newobj%aer_spec_prop(ind)%nrange)
                write(iulog ,*) 'density = ', newobj%aer_spec_prop(ind)%density
                write(iulog ,*) 'molecular_weight = ', newobj%aer_spec_prop(ind)%molecular_weight
                write(iulog ,*) 'kappa = ', newobj%aer_spec_prop(ind)%kappa
                write(iulog ,*) 'mixed = ', newobj%aer_spec_prop(ind)%mixed
                write(iulog ,*) 'nrange = ', newobj%aer_spec_prop(ind)%nrange
                write(iulog ,*) 'nbin = ', newobj%aer_spec_prop(ind)%nbin
                write(iulog ,*) 'bin_indices = ', newobj%aer_spec_prop(ind)%bin_ndx(1), ' : ', maxval(newobj%aer_spec_prop(ind)%bin_ndx)
                write(iulog ,*) 'tracer names = ', newobj%aer_spec_prop(ind)%tracernames(:newobj%aer_spec_prop(ind)%nrange)
            end do
        end if
        prop_obj => newobj
    else
        if ( .not. associated(prop_obj) ) then
            call endrun("Internal Error: fanci_aerosol_properties has not been initialized")
        end if
        newobj => prop_obj

    end if

    ! deallocate local variables
    if (allocated(nspecies)) deallocate(nspecies)
    if (allocated(alogsig)) deallocate(alogsig)
    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(fanci_aerosol_properties), intent(inout) :: self

    character(len=*), parameter :: subname = 'destructor'

    if (allocated(self%range_nspecies_)) then
        deallocate(self%range_nspecies_)
    end if

    if (allocated(self%bin_centers_)) then
        deallocate(self%bin_centers_)
    end if

    if (allocated(self%bin_bounds_)) then
        deallocate(self%bin_bounds_)
    end if

    if (allocated(self%range_bounds_)) then
        deallocate(self%range_bounds_)
    end if

    if (allocated(self%bins2ranges_)) then
        deallocate(self%bins2ranges_)
    end if

    if (allocated(self%aer_spec_prop)) then
        deallocate(self%aer_spec_prop)
    end if

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(fanci_aerosol_properties), intent(in) :: self
    character(len=*), parameter :: subname = 'number_transported'

    number_transported = self%nbins() + sum(self%range_nspecies_)

  end function number_transported

  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !  species type
  !  species name
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, list_ndx, density, hygro, &
                 spectype, specname, specmorph, refindex_sw, refindex_lw)

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    integer, optional, intent(in) :: list_ndx  ! climate or a diagnostic list number
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity
    character(len=*), optional, intent(out) :: spectype  ! species type
    character(len=*), optional, intent(out) :: specname  ! species name
    character(len=*), optional, intent(out) :: specmorph ! species morphology
    complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
    complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices

    integer :: ilist
    integer :: range_ndx_in, specprop_ndx, spec_counter, ispec
    character(len=*), parameter :: subname = 'get'

    range_ndx_in = self%bins2ranges(bin_ndx)

    ! find species properties index
    ! the input species_ndx is the index in the aero_props%indexer
    ! to find the index of the correct entry in the species_properties object array,
    ! this needs to be converted to specprop_ndx, since not all species live in all bins/ranges
    ! not every species_ndx corresponds to the same species
    ! every specprop_ndx DOES correspond to the same species

    spec_counter = 0                                                            ! make a counter for # active species in a range
    specprop_ndx = 0                                                            ! initialize index for species properties object array

    do ispec = 1, self%nspecies_tot()                                           ! loop through the species properties
        if ( any(self%aer_spec_prop(ispec)%range_ndx == range_ndx_in) ) then    ! find out if the species is in the requested range
            spec_counter = spec_counter + 1                                     ! if so, add to the counter
            if ( spec_counter == species_ndx ) then                             ! if counter is equal to species index
                specprop_ndx = ispec                                            ! then we know what species is meant!
            end if
        end if
    end do

    if (present(list_ndx)) then
        call endrun(subname//' list_ndx in fanci_aerosol_properties is not yet implemented')
    end if

    if (present(density)) then
        if ( specprop_ndx /= 0) then
            density = self%aer_spec_prop(specprop_ndx)%density
        else
            density = 0._r8
        end if
    end if

    if (present(hygro)) then
!        call endrun(subname//' hygro in fanci_aerosol_properties is not yet implemented')
! TODO: this is the hygroscopicity for each species.. should be mixed?
        if (specprop_ndx /= 0) then
            hygro = self%aer_spec_prop(specprop_ndx)%kappa
        else
            hygro = 0._r8
        end if

    end if

    if (present(spectype)) then
        if ( specprop_ndx /= 0) then
            spectype = self%aer_spec_prop(specprop_ndx)%spectype
        else
            spectype = ''
        end if
    end if

    if (present(specname)) then
        if ( specprop_ndx /= 0) then
            specname = self%aer_spec_prop(specprop_ndx)%specname
        else
            specname = ''
        end if
    end if

    if (present(specmorph)) then
        call endrun(subname//' specmorph in fanci_aerosol_properties is not yet implemented')
    end if

    if (present(refindex_sw)) then
        call endrun(subname//' refindex_sw in fanci_aerosol_properties is not yet implemented')
    end if

    if (present(refindex_lw)) then
        call endrun(subname//' refindex_lw in fanci_aerosol_properties is not yet implemented')
    end if

  end subroutine get

  !------------------------------------------------------------------------
  ! returns optics type and table parameters
  !------------------------------------------------------------------------
  subroutine optics_params(self, list_ndx, bin_ndx, opticstype, extpsw, abspsw, asmpsw, absplw, &
       refrtabsw, refitabsw, refrtablw, refitablw, ncoef, prefr, prefi, sw_hygro_ext_wtp, &
       sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, wgtpct, nwtp, &
       sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, &
       corefrac, bcdust, kap, relh, nfrac, nbcdust, nkap, nrelh )

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: list_ndx            ! rad climate/diags list

    character(len=*), optional, intent(out) :: opticstype

    ! refactive index table parameters
    real(r8),  optional, pointer     :: extpsw(:,:,:,:) ! short wave specific extinction
    real(r8),  optional, pointer     :: abspsw(:,:,:,:) ! short wave specific absorption
    real(r8),  optional, pointer     :: asmpsw(:,:,:,:) ! short wave asymmetry factor
    real(r8),  optional, pointer     :: absplw(:,:,:,:) ! long wave specific absorption
    real(r8),  optional, pointer     :: refrtabsw(:,:)  ! table of short wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitabsw(:,:)  ! table of short wave imaginary refractive indices for aerosols
    real(r8),  optional, pointer     :: refrtablw(:,:)  ! table of long wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitablw(:,:)  ! table of long wave imaginary refractive indices for aerosols
    integer,   optional, intent(out) :: ncoef  ! number of chebychev polynomials
    integer,   optional, intent(out) :: prefr  ! number of real refractive indices in table
    integer,   optional, intent(out) :: prefi  ! number of imaginary refractive indices in table

    ! hygrowghtpct table parameters
    real(r8),  optional, pointer     :: sw_hygro_ext_wtp(:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_ssa_wtp(:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_asm_wtp(:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_ext_wtp(:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: wgtpct(:)   ! weight precent of H2SO4/H2O solution
    integer,   optional, intent(out) :: nwtp        ! number of weight precent values

    ! hygrocoreshell table parameters
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: corefrac(:) ! core fraction dimension values
    real(r8),  optional, pointer     :: bcdust(:)   ! bc/(bc + dust) fraction dimension values
    real(r8),  optional, pointer     :: kap(:)      ! hygroscopicity dimension values
    real(r8),  optional, pointer     :: relh(:)     ! relative humidity dimension values
    integer,   optional, intent(out) :: nfrac       ! core fraction dimension size
    integer,   optional, intent(out) :: nbcdust     ! bc/(bc + dust) fraction dimension size
    integer,   optional, intent(out) :: nkap        ! hygroscopicity dimension size
    integer,   optional, intent(out) :: nrelh       ! relative humidity dimension size

    character(len=*), parameter :: subname = 'optics_params'

    call endrun(subname//' is not yet implemented')

  end subroutine optics_params

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    character(len=*), parameter :: subname = 'amcube'
! copied from CARMA
    amcube = 3._r8/(4._r8*pi)*volconc/numconc

  end function amcube

  !------------------------------------------------------------------------------
  ! returns density for a species
  !------------------------------------------------------------------------------
  real(r8) function density(self, species_ndx)

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx                          ! element in the sec_aero_props array

    if (species_ndx > 0) then
        density = self%aer_spec_prop(species_ndx)%density
    else
        density = 0._r8
    end if
  end function density

  !------------------------------------------------------------------------------
  ! returns kappa for a species (hygroscopicity)
  !------------------------------------------------------------------------------
  real(r8) function kappa(self, species_ndx)

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx

    if (species_ndx > 0) then
       kappa = self%aer_spec_prop(species_ndx)%kappa
    else
       kappa = 0._r8
    end if

  end function kappa

  !------------------------------------------------------------------------------
  ! returns molecular_weight for a species
  !------------------------------------------------------------------------------
  real(r8) function molecular_weight(self, species_ndx)

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx

    molecular_weight = self%aer_spec_prop(species_ndx)%molecular_weight

  end function molecular_weight

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    character(len=*), parameter :: subname = 'actfracs'
! TODO: check activated fraction

    fn = 0._r8
    fm = 0._r8

    if (smc < smax) then
       fn = 1._r8
       fm = 1._r8
    end if
!    call endrun(subname//' is not yet implemented')

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    character(len=*), parameter :: subname = 'num_names'

    name_a = 'num_'//trim(int2str(bin_ndx))
    name_c = 'num_'//trim(int2str(bin_ndx))//'_cw'

  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    integer :: range_ndx, ispec, specprop_ndx, spec_counter
    character(len=10) :: specname
    character(len=*), parameter :: subname = 'mmr_names'


    ! TODO: only for for first bin in range
    range_ndx = self%bins2ranges(bin_ndx)
    if (self%range_bounds(range_ndx,1) == bin_ndx) then

! TODO: make this into a function!
    spec_counter = 0
    do ispec = 1, self%nspecies_tot()                                           ! loop through the species properties
        if ( any(self%aer_spec_prop(ispec)%range_ndx == range_ndx) ) then    ! find out if the species is in the requested range
            spec_counter = spec_counter + 1                                     ! if so, add to the counter
                if ( spec_counter == species_ndx ) then                             ! if counter is equal to species index
                specprop_ndx = ispec                                            ! then we know what species is meant!
                specname = self%aer_spec_prop(specprop_ndx)%specname
                name_a = trim(specname)//'_R'//trim(int2str(range_ndx))
                name_c = trim(name_a)//'_cw'
            end if
        end if
    end do

else
    name_a = ''
    name_c = ''

end if
  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens
    character(len=*), parameter :: subname = 'amb_num_name'

    call endrun(subname//' is not yet implemented')

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR
    character(len=*), parameter :: subname = 'amb_mmr_name'

    call endrun(subname//' is not yet implemented')

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: spectype ! species type
    character(len=*), parameter :: subname = 'species_type'

    call self%get(bin_ndx, species_ndx, spectype=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=aero_name_len) :: spectype
    character(len=aero_name_len) :: modetype
    integer :: spc_ndx
    character(len=*), parameter :: subname = 'icenuc_updates_num'

    res = .false.

    call endrun(subname//' is not yet implemented')

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=32) :: spectype
    character(len=32) :: modetype
    character(len=*), parameter :: subname = 'icenuc_updates_mmr'

    res = .false.

    call endrun(subname//' is not yet implemented')

  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(fanci_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index
    character(len=*), parameter :: subname = 'apply_number_limits'
! TODO: is this a modal thing?

    return

    !call endrun(subname//' is not yet implemented')

  end subroutine apply_number_limits

  !------------------------------------------------------------------------------
  ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
  ! the particles' ability to act as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_species(self, bin_ndx, spc_ndx) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    integer, intent(in) :: spc_ndx  ! species number

    logical :: res
    character(len=*), parameter :: subname = 'hetfrz_species'

    call endrun(subname//' is not yet implemented')

  end function hetfrz_species

  !------------------------------------------------------------------------------
  ! returns TRUE if soluble
  !------------------------------------------------------------------------------
  logical function soluble(self,bin_ndx)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), parameter :: subname = 'soluble'

    call endrun(subname//' is not yet implemented')

  end function soluble

  !------------------------------------------------------------------------------
  ! returns minimum mass mean radius (meters)
  !------------------------------------------------------------------------------
  function min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    real(r8) :: minrad  ! meters

    integer :: nbins
    character(len=*), parameter :: subname = 'min_mass_mean_rad'

    call endrun(subname//' is not yet implemented')

  end function min_mass_mean_rad

  !------------------------------------------------------------------------------
  ! returns the total number of bins for a given radiation list index
  !------------------------------------------------------------------------------
  function nbins_rlist(self, list_ndx)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx  ! radiation list number

    integer :: res
    character(len=*), parameter :: subname = 'nbins_rlist'

    call endrun(subname//' is not yet implemented')

  end function nbins_rlist

  !------------------------------------------------------------------------------
  ! returns the total number of species objects
  !------------------------------------------------------------------------------
  function nspecies_tot(self)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer :: res
    character(len=*), parameter :: subname = 'nspecies_tot'

    res = self%nspecies_tot_

  end function nspecies_tot

  !------------------------------------------------------------------------------
  ! returns the total number of species in each range
  !------------------------------------------------------------------------------
  function range_nspecies(self, irange)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: irange
    integer :: res
    character(len=*), parameter :: subname = 'range_nspecies'

    res = self%range_nspecies_(irange)

  end function range_nspecies

  !------------------------------------------------------------------------------
  ! returns the total number of ranges
  !------------------------------------------------------------------------------
  function nranges(self)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer :: res
    character(len=*), parameter :: subname = 'nranges'

    res = self%nranges_

  end function nranges

  !------------------------------------------------------------------------------
  ! returns bin centers
  !------------------------------------------------------------------------------
  function bin_centers(self, nbins) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: nbins
    real(r8) :: res(nbins)
    character(len=*), parameter :: subname = 'bin_centers'

    res = self%bin_centers_(:min(nbins,size(self%bin_centers_)))

  end function bin_centers

  !------------------------------------------------------------------------------
  ! returns bin bounds
  !------------------------------------------------------------------------------
  function bin_bounds(self, bin_ndx, bound) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx
    integer, intent(in) :: bound     ! 1 for lower, 2 for upper
    real(r8) :: res
    character(len=*), parameter :: subname = 'bin_bounds'

    res = self%bin_bounds_(bin_ndx, bound)

  end function bin_bounds

  !------------------------------------------------------------------------------
  ! returns range bounds
  !------------------------------------------------------------------------------
  integer function range_bounds(self, irange, bound) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: irange
    integer, intent(in) :: bound ! 1 for lower bound, 2 for upper bound
    character(len=*), parameter :: subname = 'range_bounds'

    res = self%range_bounds_(irange, bound)

  end function range_bounds

  !------------------------------------------------------------------------------
  ! returns volume of a particle in a bin
  !------------------------------------------------------------------------------
  function particle_volume(self, ibin) result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: ibin
    real(r8) :: res
    character(len=*), parameter :: subname = 'particle_volume'

    res = self%particle_volume_(ibin)

  end function particle_volume
  !------------------------------------------------------------------------------
  ! returns the index of a range
  !------------------------------------------------------------------------------
  function spec_range_ndx(self, specprop_ndx, ind)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: specprop_ndx ! location in the aer_spec_props object
    integer, intent(in) :: ind          ! location in the aer_spec_props(species)%range_ndx array
    integer :: res
    character(len=*), parameter :: subname = 'spec_range_ndx'

    res = self%aer_spec_prop(specprop_ndx)%range_ndx(ind)

  end function spec_range_ndx

  !------------------------------------------------------------------------------
  ! returns the bin indices for a species
  !------------------------------------------------------------------------------
  function spec_bin_ndx(self, specprop_ndx, ind)  result(res)
    ! TODO: needed?
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: ind, specprop_ndx
    integer :: ibin
    integer :: res
    character(len=*), parameter :: subname = 'spec_bin_ndx'

    res = self%aer_spec_prop(specprop_ndx)%bin_ndx(ind)

  end function spec_bin_ndx

  !------------------------------------------------------------------------------
  ! returns the upper or lower range idx TODO: change to array of idices?
  !------------------------------------------------------------------------------
  function spec_nrange(self, species_ndx)  result(res)
    ! TODO: needed?
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx
    integer :: res
    character(len=*), parameter :: subname = 'spec_nrange'

    res = self%aer_spec_prop(species_ndx)%nrange

  end function spec_nrange

  !------------------------------------------------------------------------------
  ! returns number of bins containing species
  !------------------------------------------------------------------------------
  function spec_nbin(self, specprop_ndx, spectype)  result(res)
    ! TODO: needed?
    class(fanci_aerosol_properties), intent(in) :: self
    integer, optional, intent(in)   :: specprop_ndx
    character(len=*), optional, intent(in) :: spectype
    integer :: ispecprop
    integer :: res
    character(len=*), parameter :: subname = 'spec_nbin'


    if (present(spectype)) then
        do ispecprop = 1, self%nspecies_tot()
            if (trim(spectype) == trim(self%aer_spec_prop(ispecprop)%spectype)) then
                res = self%aer_spec_prop(ispecprop)%nbin
            end if
        end do
    else if (present(specprop_ndx)) then
        res = self%aer_spec_prop(specprop_ndx)%nbin
    else
        call endrun("spec_nbin could not be defined")
    end if

  end function spec_nbin

  !------------------------------------------------------------------------------
  ! returns number of bins containing species
  !------------------------------------------------------------------------------
  function spec_tracernames(self, species_ndx, spec_range_ndx)  result(res)
    ! TODO: needed?
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx, spec_range_ndx
    character(len=10) :: res
    character(len=*), parameter :: subname = 'spec_tracernames'

    res = self%aer_spec_prop(species_ndx)%tracernames(spec_range_ndx)

  end function spec_tracernames

  !------------------------------------------------------------------------------
  ! returns the total number of species objects
  !------------------------------------------------------------------------------

  integer function bins2ranges(self, ibin)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in)         :: ibin
    character(len=*), parameter :: subname = 'bins2ranges'

    res = self%bins2ranges_(ibin)

  end function bins2ranges
  !------------------------------------------------------------------------------
  ! returns number of species in a bin for a given radiation list index
  !------------------------------------------------------------------------------
  function nspecies_per_bin_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    integer :: res
    character(len=*), parameter :: subname = 'nspecies_per_bin_rlist'

    call endrun(subname//' is not yet implemented')

  end function nspecies_per_bin_rlist

  !------------------------------------------------------------------------------
  ! returns the natural log of geometric standard deviation of the number
  ! distribution for radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function alogsig_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: res
    character(len=*), parameter :: subname = 'alogsig_rlist'

    res = 0._r8                     ! sigma is 1, so logsig is 0
    call endrun(subname//' is not yet implemented')

  end function alogsig_rlist

  !------------------------------------------------------------------------------
  ! returns name for a given radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, list_ndx,  bin_ndx) result(name)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=32) name
    character(len=*), parameter :: subname = 'bin_name'

    call endrun(subname//' is not yet implemented')

  end function bin_name

  !------------------------------------------------------------------------------
  ! returns scavenging diameter (cm) for a given aerosol bin number
  !------------------------------------------------------------------------------
  function scav_diam(self, bin_ndx) result(diam)
    use modal_aero_data, only: dgnum_amode

    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: diam

    character(len=*), parameter :: subname = 'scav_diam'

    diam = self%bin_centers_(bin_ndx) * 100._r8 * 2._r8  ! convert to cm

  end function scav_diam

  !------------------------------------------------------------------------------
  ! adjust aerosol concentration tendencies to create larger sizes of aerosols
  ! during resuspension
  !------------------------------------------------------------------------------
  subroutine resuspension_resize(self, dcondt)

    use modal_aero_data, only:  mode_size_order

    class(fanci_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: dcondt(:)
    character(len=*), parameter :: subname = 'resuspension_resize'

    call endrun(subname//' is not yet implemented')

  end subroutine resuspension_resize

  !------------------------------------------------------------------------------
  ! returns bulk deposition fluxes of the specified species type
  ! rebinned to specified diameter limits
  !------------------------------------------------------------------------------
  subroutine rebin_bulk_fluxes(self, bulk_type, dep_fluxes, diam_edges, bulk_fluxes, &
                               error_code, error_string)
    use infnan, only: nan, assignment(=)

    class(fanci_aerosol_properties), intent(in) :: self
    character(len=*),intent(in) :: bulk_type       ! aerosol type to rebin
    real(r8), intent(in) :: dep_fluxes(:)          ! kg/m2 -> for each bin!
    real(r8), intent(in) :: diam_edges(:)          ! meters
    real(r8), intent(out) :: bulk_fluxes(:)        ! kg/m2
    integer,  intent(out) :: error_code            ! error code (0 if no error)
    character(len=*), intent(out) :: error_string  ! error string

    integer :: irange, ispecprop, nbulk, imas, ibin, mm, ibulk
    logical :: type_not_found
    real(r8) :: spec_density, density
    character(len=10) :: spectype
    real(r8), allocatable :: bin_dep_flux(:), species_bin_dep_flux(:), adjusted_diameter(:)
    real(r8), allocatable :: range_volume(:), range_dep_flux(:), spec_dep_flux_in(:), species_massfrac(:), spec_dep_flux(:), range_density(:)

    character(len=*), parameter :: subname = 'rebin_bulk_fluxes'

    allocate(bin_dep_flux(self%nbins()), species_bin_dep_flux(self%nbins()), adjusted_diameter(self%nbins()))
    allocate(range_volume(self%nranges()), range_dep_flux(self%nranges()), spec_dep_flux_in(self%nranges()))
    allocate(range_density(self%nranges()), species_massfrac(self%nranges()), spec_dep_flux(self%nranges()))

    bin_dep_flux = 0._r8
    range_dep_flux = 0._r8
    range_volume = 0._r8
    spec_dep_flux = 0._r8
    spec_dep_flux_in = 0._r8
    species_bin_dep_flux = 0._r8
    species_massfrac = 0._r8
    adjusted_diameter = 0._r8
    range_density = 0._r8

    ! initialize
    error_code = 0
    error_string = ' '
    bulk_fluxes = 0._r8
    type_not_found = .true.

    nbulk = size(bulk_fluxes)

    ! Find out if the species exists
    do ispecprop = 1, self%nspecies_tot()
        if ( self%aer_spec_prop(ispecprop)%spectype == bulk_type ) then
            type_not_found = .false.
            exit
        end if
    end do

    ! if species doesn't exist, make an error
    if (type_not_found) then
        bulk_fluxes(:) = nan
        error_code = 1
        write(error_string,*) 'ERROR:: aerosol_properties:: ',subname,' bulk_type: ', bulk_type, ' not found'
    else
        do ibin = 1, self%nbins()
            irange = self%bins2ranges(ibin)                                                                 ! find range corresponding to bin
            do imas = 0, self%nmasses(ibin)                                                                      ! loop through the indexer
                mm = self%indexer(ibin,imas)
                if (imas == 0) then                                                                         ! when imas = 0, index is for bin number conc
                    bin_dep_flux(ibin) = dep_fluxes(mm)
                    range_volume(irange) = range_volume(irange) + dep_fluxes(mm) * self%particle_volume(ibin) ! sum up the volume for density calculation
                else
                    range_dep_flux(irange) = range_dep_flux(irange) + dep_fluxes(mm)                        ! sum up the deposition flux in kg/m2 for the range
                    call self%get(ibin, mm, density=density, spectype=spectype)                             ! get the species density and species type
                    if (trim(spectype) == trim(bulk_type)) then                                                         ! get dep flux for species to be re-binned
                        spec_dep_flux_in(irange) = spec_dep_flux(irange) + dep_fluxes(mm)
                        spec_density = density
                    end if
                end if
            end do
        end do

        where (range_volume /= 0._r8 .and. range_dep_flux /= 0._r8)
            range_density = range_dep_flux/range_volume
            species_massfrac = spec_dep_flux_in/range_dep_flux
        end where

        ! get the mass deposited for one species in each bin
        do ibin = 1, self%nbins()
            irange = self%bins2ranges(ibin)
            species_bin_dep_flux(ibin) = species_massfrac(irange) * bin_dep_flux(ibin) * self%particle_volume(ibin)*range_density(irange) ! mass fraction * mass in a bin

            ! since all aerosol is mixed and only one species is re-binned here, we need to find the new smaller radius
            if (spec_density /= 0._r8 .and. bin_dep_flux(ibin) /= 0._r8) then
                adjusted_diameter(ibin) = ( ( ( species_bin_dep_flux(ibin) * 3._r8 ) / ( spec_density * bin_dep_flux(ibin) * 4._r8 * pi ) )**(1._r8/3._r8) ) * 2 ! convert to diameter
                if ( adjusted_diameter(ibin) < diam_edges(1) ) then
                    bulk_fluxes(1) = bulk_fluxes(1) + species_bin_dep_flux(ibin)
                else if ( adjusted_diameter(ibin) > diam_edges(size(diam_edges))) then
                    bulk_fluxes(size(bulk_fluxes)) = bulk_fluxes(size(bulk_fluxes)) + species_bin_dep_flux(ibin)
                else
                    do ibulk = 1, size(bulk_fluxes)-1
                        if ( adjusted_diameter(ibin) > diam_edges(ibulk) .and. adjusted_diameter(ibin) < diam_edges(ibulk+1) ) then
                            bulk_fluxes(ibulk) = bulk_fluxes(ibulk) + species_bin_dep_flux(ibin)
                        end if
                    end do
                end if
            end if
        end do

    end if

  end subroutine rebin_bulk_fluxes

  !------------------------------------------------------------------------------
  ! Returns TRUE if bin is hydrophilic, otherwise FALSE
  !------------------------------------------------------------------------------
  logical function hydrophilic(self, bin_ndx)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx ! bin number

    character(len=aero_name_len) :: modetype
    character(len=*), parameter :: subname = 'hydrophilic'

! TODO: check sizing
! depends on size -> modal: "accum"
    hydrophilic = ( self%bin_centers_(bin_ndx)*1e9 > 50._r8 .and. self%bin_centers_(bin_ndx)*1e9 < 500._r8 )

  end function hydrophilic

  !------------------------------------------------------------------------------
  ! returns TRUE if fanci aerosol representation
  !------------------------------------------------------------------------------
  pure logical function model_is(self, query)
    class(fanci_aerosol_properties), intent(in) :: self
    character(len=*),               intent(in) :: query

    if (trim(query) == 'fanci') then
       model_is = .true.
    else
       model_is = .false.
    end if

  end function model_is
  !------------------------------------------------------------------------------
  ! returns TRUE if given species is active
  !------------------------------------------------------------------------------
  pure logical function is_active(self, spectype)
    class(fanci_aerosol_properties), intent(in) :: self
    character(len=*),                    intent(in) :: spectype
    integer :: ispec

    is_active = .false.

    do ispec = 1, self%nspecies_tot_
        if (trim(self%aer_spec_prop(ispec)%spectype) == trim(spectype)) then
            is_active = .true.
        end if
    end do

  end function is_active

  !------------------------------------------------------------------------------
  ! Returns specname with aer_spec_props index
  !------------------------------------------------------------------------------
  function specname(self, specprop_ndx) result(species_name)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in)           :: specprop_ndx
    character(len=:), allocatable :: species_name

    species_name = self%aer_spec_prop(specprop_ndx)%specname

  end function specname
  !------------------------------------------------------------------------------
  ! Returns spectype with aer_spec_props index
  !------------------------------------------------------------------------------
  function spectype(self, specprop_ndx) result(species_type)
    class(fanci_aerosol_properties), intent(in) :: self
    integer, intent(in)           :: specprop_ndx
    character(len=:), allocatable :: species_type

    species_type = self%aer_spec_prop(specprop_ndx)%spectype

  end function spectype

  !------------------------------------------------------------------------------
  ! Returns CAM q array index for species bin number concentrations
  !------------------------------------------------------------------------------
  function spec_bin_q_ndx(self, specprop_ndx, spec_bin_ndx) result(bin_q_ndx)
     use cam_history, only: fieldname_len
     use constituents, only: cnst_get_ind

     class(fanci_aerosol_properties), intent(in) :: self
     integer, intent(in) :: specprop_ndx, spec_bin_ndx
     character(len=fieldname_len) :: dummy
     integer             :: bin_q_ndx

     bin_q_ndx = 0
     dummy = 'num_'//int2str(self%aer_spec_prop(specprop_ndx)%bin_ndx(spec_bin_ndx))
     call cnst_get_ind(dummy, bin_q_ndx)

  end function spec_bin_q_ndx
  !------------------------------------------------------------------------------
  ! Returns CAM q array index for species range mmr
  !------------------------------------------------------------------------------
  function spec_mmr_q_ndx(self, specprop_ndx, spec_range_ndx) result(mmr_q_ndx)
     use cam_history, only: fieldname_len
     use constituents, only: cnst_get_ind

     class(fanci_aerosol_properties), intent(in) :: self
     integer, intent(in) :: specprop_ndx, spec_range_ndx
     character(len=fieldname_len) :: dummy
     integer             :: mmr_q_ndx

     mmr_q_ndx = 0

     dummy = self%aer_spec_prop(specprop_ndx)%tracernames(spec_range_ndx)
     call cnst_get_ind(dummy, mmr_q_ndx)

  end function spec_mmr_q_ndx


end module fanci_aerosol_properties_mod
