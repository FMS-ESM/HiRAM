
module ocean_model_mod

!-----------------------------------------------------------------------

use time_manager_mod, only: time_type, operator(+), operator(>), &
                            get_date, set_time

use          fms_mod, only: file_exist, open_restart_file, &
                            close_file, mpp_pe, mpp_root_pe, mpp_npes,         &
                            write_version_number, stdlog, error_mesg, WARNING, FATAL, &
                            check_nml_error, write_data, set_domain, NOTE, &
                            field_exist, get_mosaic_tile_grid, read_data,  &
                            field_size

#ifdef INTERNAL_FILE_NML
use          mpp_mod, only: input_nml_file
#else
use          fms_mod, only: open_namelist_file
#endif

use  fms_io_mod,      only: get_restart_io_mode

use  amip_interp_mod, only: amip_interp_type, amip_interp_new,  &
                            amip_interp_del, get_amip_sst

use  mpp_domains_mod, only: domain1d, domain2d, mpp_define_domains,  &
                            mpp_get_compute_domain, mpp_get_compute_domains,  &
                            mpp_get_domain_components, mpp_get_pelist,  &
                            CYCLIC_GLOBAL_DOMAIN, mpp_define_layout
use          mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, stdout,mpp_error, mpp_chksum

use    constants_mod, only: PI, RADIUS

use mosaic_mod,       only: get_mosaic_ntiles, get_mosaic_grid_sizes, get_mosaic_xgrid
use mosaic_mod,       only: get_mosaic_xgrid_size, calc_mosaic_grid_area

use coupler_types_mod,only: coupler_2d_bc_type

#ifdef SCM
use     scm_forc_mod, only: do_specified_tskin, TSKIN
#endif

implicit none
private

!-----------------------------------------------------------------------
!----------------- public interfaces -----------------------------------

public :: ocean_model_init, ocean_model_end, update_ocean_model, ocean_public_type, &
          ice_ocean_boundary_type, ocean_grids_type, &
          ocean_model_flux_init, ocean_model_init_sfc, ocean_stock_pe, ocean_model_restart
public :: ice_ocn_bnd_type_chksum, ocean_public_type_chksum

public    ocean_model_data_get
interface ocean_model_data_get
   module procedure ocean_model_data1D_get 
   module procedure ocean_model_data2D_get 
end interface

!
! the following type is for data exchange with the new coupler
! it is defined here but declared in coupler_main and allocated in flux_init
!
type ice_ocean_boundary_type
  real, dimension(:,:), pointer :: u_flux =>NULL(), &
                                   v_flux =>NULL(), &
                                   t_flux =>NULL(), &
                                   q_flux =>NULL(), &
                                   salt_flux =>NULL(), &
                                   lw_flux =>NULL(), &
                                   sw_flux_vis_dir =>NULL(), &
                                   sw_flux_vis_dif =>NULL(), &
                                   sw_flux_nir_dir =>NULL(), &
                                   sw_flux_nir_dif =>NULL(), &
                                   lprec =>NULL(), &
                                   fprec  =>NULL()
  real, dimension(:,:), pointer :: runoff =>NULL(), &
                                   calving  =>NULL(), &
                                   runoff_hflx  =>NULL(), &
                                   calving_hflx  =>NULL()
  real, dimension(:,:), pointer :: p  =>NULL()
  ! "data" is collective field for "named" fields above
  real, dimension(:,:,:), pointer :: data  =>NULL()
  integer :: xtype             !REGRID, REDIST or DIRECT used by coupler
  type(coupler_2d_bc_type)      :: fluxes  ! array of fields used for additional tracers
end type ice_ocean_boundary_type

!-----------------------------------------------------------------------

 type ocean_grids_type
    real,    pointer, dimension(:)   :: lon_bnd =>NULL(), lat_bnd =>NULL()
    real,    pointer, dimension(:,:)   :: lon =>NULL(), lat =>NULL()
    logical, pointer, dimension(:,:) :: mask  =>NULL()
 end type

!   lon_bnd    = longitude boundaries for grid boxes
!   lat_bnd    = latitude  boundaries for grid boxes
!   mask       = land-sea mask for grid boxes
!
!    note: longitude/latitude is in radians
!          mask is true for ocean points
!
!-----------------------------------------------------------------------

type ocean_public_type
   type (domain2d)               :: Domain
   type (ocean_grids_type)       :: Global, Data
   real, pointer, dimension(:,:) :: t_surf =>NULL() , &
                                    frazil =>NULL() ,  &
                                    u_surf =>NULL() , &
                                    v_surf =>NULL() , &
                                    s_surf =>NULL() , &
                                    area   =>NULL() , &
                                    sea_lev =>NULL()
   logical, pointer, dimension(:,:) :: maskmap =>NULL()! A pointer to an array indicating which
                                                       ! logical processors are actually used for
                                                       ! the ocean code. The other logical
                                                       ! processors would be all land points and
                                                       ! are not assigned to actual processors.
                                                       ! This need not be assigned if all logical
                                                       ! processors are used. This variable is dummy and need 
                                                       ! not to be set, but it is needed to pass compilation.
   type (time_type)              :: Time, &
                                    Time_step
   logical :: is_ocean_pe
   integer, pointer :: pelist(:) =>NULL()
   integer, dimension(3)            :: axes    
   type(coupler_2d_bc_type)         :: fields  ! array of fields used for additional tracers
end type ocean_public_type

!  Global = grid information for the global data grid
!  Data   = grid information for the local data grid
!
!  t_surf      = surface temperature on the local ocean model grid
!  frazil      = frazil on the local ocean model grid
!  u_surf      = zonal ocean current on the local ocean model grid
!  v_surf      = meridional ocean current on the local ocean model grid
!
!  Time      = current ocean model time
!  Time_step = ocean model time step
!
!    Notes: Global grid information will be the same on all processors
!           The data grid is global, but only the portion local to
!             the current processor is store in Data.
!-----------------------------------------------------------------------

  type, public ::  ocean_state_type; private
     ! This type is private, and can therefore vary between different ocean models.
     ! All information entire ocean state may be contained here, although it is not
     ! necessary that this is implemented with all models.
     logical       :: is_ocean_pe = .false.       ! .true. on processors that run the ocean model.
  end type ocean_state_type

!------- namelist ---------
   logical :: do_netcdf_restart = .true.
   logical :: use_climo_sst  = .false.
   logical :: use_annual_sst = .false.
   integer, dimension(2) :: layout = (/ 0, 0 /)
   character(len=64) :: interp_method  = "conservative" ! default, conservative scheme

!  layout =  domain decomposition (# X-PEs by # Y-PEs)
!             layout = (/0,0/) will use default rules
!
   namelist /ocean_model_nml/ do_netcdf_restart, use_climo_sst, use_annual_sst, layout, interp_method

!-----------------------------------------------------------------------
!--------------------- private below here ------------------------------

!  ---- version number -----

   character(len=128) :: version = '$Id: ocean_model.F90,v 17.0.2.1.4.1 2010/11/15 18:32:26 bw Exp $'
   character(len=128) :: tagname = '$Name: hiram_20101115_bw $'

!-----------------------------------------------------------------------
!------ model resolution parameters ------

   type (amip_interp_type), save :: Amip

   logical :: module_is_initialized=.false.
   logical :: stock_warning_issued =.false.

contains

!#######################################################################

 !#######################################################################
! <SUBROUTINE NAME="update_ocean_model">
!
! <DESCRIPTION>
! Update in time the ocean model fields. 
!   This subroutine uses the forcing in Ice_ocean_boundary to advance the
! ocean model's state from the input value of Ocean_state (which must be for
! time time_start_update) for a time interval of Ocean_coupling_time_step,
! returning the publicly visible ocean surface properties in Ocean_sfc and
! storing the new ocean properties in Ocean_state.
!
! Arguments: 
!  Ice_ocean_boundary - A structure containing the various forcing
!                                 fields coming from the ice. It is intent in.
!  Ocean_state - A structure containing the internal ocean state.
!  Ocean_sfc - A structure containing all the publicly visible ocean
!                        surface fields after a coupling time step.
!  time_start_update - The time at the beginning of the update step.
!  Ocean_coupling_time_step - The amount of time over which to advance
!                                       the ocean.

! Note: although several types are declared intent(inout), this is to allow for
!   the possibility of halo updates and to keep previously allocated memory.
!   In practice, Ice_ocean_boundary is intent in, Ocean_state is private to
!   this module and intent inout, and Ocean_sfc is intent out.
! </DESCRIPTION>
!
  subroutine update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, &
       time_start_update, Ocean_coupling_time_step)
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
    type(ocean_state_type),        pointer       :: Ocean_state
    type(ocean_public_type),       intent(inout) :: Ocean_sfc
    type(time_type), intent(in)                  :: time_start_update
    type(time_type), intent(in)                  :: Ocean_coupling_time_step
    
    integer :: num_ocean_calls, no

 !check if required boundary fields have been initialized
      if ( .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_vis_dir) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_vis_dif) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_nir_dir) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_nir_dif) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%lw_flux) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%fprec)   .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%calving) .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%t_flux)  .OR. &
        .NOT.ASSOCIATED(Ice_ocean_boundary%q_flux) ) &
          call error_mesg( 'Update_ocean_model', &
           'Ice_ocean_boundary not correctly initialized.', FATAL )

!-----------------------------------------------------------------------
! ---- update time ----

      Ocean_sfc%Time = Ocean_sfc%Time + Ocean_coupling_time_step
!-----------------------------------------------------------------------
!----- update sst, set currents -----

      call set_ocean_model_state (Ocean_sfc)

!-----------------------------------------------------------------------

 end subroutine update_ocean_model

!#######################################################################

 subroutine set_ocean_model_state (Ocean_sfc)

   type (ocean_public_type), intent(inout) :: Ocean_sfc

!-----------------------------------------------------------------------
!----- get new sea surface temperature ------

       call get_amip_sst ( Ocean_sfc%Time, Amip, Ocean_sfc%t_surf )

#ifdef SCM
!--- for single column model -------------------------------------!
!--- initialize surface temperature to observed value ------------!  
       if (do_specified_tskin) then     
          Ocean_sfc%t_surf = TSKIN
       end if
!-----------------------------------------------------------------!
#endif

!-----------------------------------------------------------------------
!----- currents -----

   Ocean_sfc%u_surf = 0.0
   Ocean_sfc%v_surf = 0.0

!----- dummy out frazil ??? -----

   Ocean_sfc%frazil = 0.0
   Ocean_sfc%area   = 1.0

   Ocean_sfc%s_surf  = 0.0
   Ocean_sfc%sea_lev = 0.0
!-----------------------------------------------------------------------

 end subroutine set_ocean_model_state

!#######################################################################

subroutine ocean_model_init (Ocean, Ocean_state, Time_init, Time)
!
! <DESCRIPTION>
! Initialize the ocean model. 
! Arguments: 
!  Ocean (inout)  - A structure containing various publicly visible ocean
!                    surface properties after initialization.
!  Ocean_state (pointer)- A structure whose internal contents are private
!                    to ocean_model_mod that may be used to contain all
!                    information about the ocean's interior state.
!  Time_init (in) - The start time for the coupled model's calendar.
!  Time_in   (in) - The time at which to initialize the ocean model.
! </DESCRIPTION>

  type (ocean_public_type), intent(inout) :: Ocean
  type (ocean_state_type),  pointer       :: Ocean_state
  type (time_type),         intent(in)    :: Time_init, Time

  integer                              :: siz(4)
  integer                              :: unit, ierr, io, nlon, nlat
  integer                              :: isd, ied, jsd, jed
  integer                              :: i, j, npes
  logical, allocatable, dimension(:,:) :: global_mask
  integer, allocatable, dimension(:)   :: xextent, yextent
  real,    allocatable, dimension(:,:) :: geo_lonv, geo_latv, rmask, netr_mask
  real,    allocatable, dimension(:,:) :: geo_lont, geo_latt
  real,  allocatable, dimension(:,:,:) :: x_vert_T, y_vert_T
  real,    allocatable, dimension(:,:) :: tmpx, tmpy, garea
  real, allocatable, dimension(:)      :: xgrid_area(:)
  integer, allocatable, dimension(:)   :: i1, j1, i2, j2
  character(len=256)                   :: err_mesg
  character(len=80)                    :: domainname
  character(len=256)                   :: grid_file = "INPUT/grid_spec.nc"
  character(len=256)                   :: ocean_mosaic, tile_file
  character(len=256)                   :: axo_file      ! atmosXocean exchange grid file
  integer                              :: nx(1), ny(1)
  integer                              :: ntiles, nfile_axo, nxgrid, n, m
  integer                              :: grid_version
  integer, parameter                   :: VERSION_0 = 0  ! grid file with field geolon_t
  integer, parameter                   :: VERSION_1 = 1  ! grid file with field x_T
  integer, parameter                   :: VERSION_2 = 2  ! mosaic file


   if(module_is_initialized) return

     Ocean%Time      = Time

!   ----- read namelist -----

     if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=ocean_model_nml, iostat=io)
        ierr = check_nml_error(io, 'ocean_model_nml')
#else
        unit = open_namelist_file ( )
        ierr=1
         do while (ierr /= 0)
           read  (unit, nml=ocean_model_nml, iostat=io, end=10)
           ierr = check_nml_error(io,'ocean_model_nml')
        enddo
 10     call close_file (unit)
#endif
     endif
     call get_restart_io_mode(do_netcdf_restart)

!   ----- write version number and namelist -----
     call write_version_number(version, tagname)

     if ( mpp_pe() == mpp_root_pe() ) then
          write (stdlog(),nml=ocean_model_nml)
     endif

    !--- get the grid size 
    if(field_exist(grid_file, 'geolon_t')) then
       grid_version = VERSION_0 
       call field_size( grid_file, 'geolon_t', siz)  
       nlon = siz(1)
       nlat = siz(2)        
    else if(field_exist(grid_file, 'x_T')) then
       grid_version = VERSION_1
       call field_size( grid_file, 'x_T', siz)
       nlon = siz(1)
       nlat = siz(2) 
    else if(field_exist(grid_file, 'ocn_mosaic_file') ) then ! read from mosaic file
       grid_version = VERSION_2
       call read_data(grid_file, "ocn_mosaic_file", ocean_mosaic)
       ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
       ntiles = get_mosaic_ntiles(ocean_mosaic)
       if(ntiles .NE. 1) call error_mesg('ocean_model_init', ' ntiles should be 1 for ocean mosaic, contact developer', FATAL)
       call get_mosaic_grid_sizes( ocean_mosaic, nx, ny)
       nlon = nx(1)
       nlat = ny(1)
    else
       call error_mesg('ocean_model_init','x_T, geolon_t, ocn_mosaic_file does not exist in file '//trim(grid_file), FATAL )
    end if

    if( nlon .LE. 0 .or. nlat .LE. 0 ) call error_mesg('ocean_model_init', 'nlon and nlat should be a positive integer.', FATAL)

!-----------------------------------------------------------------------
!----- set up global storage and local storage -----

   allocate ( Ocean%Global%lon_bnd (nlon+1)  ,  &
              Ocean%Global%lat_bnd (nlat+1)  ,  &
              Ocean%Global%lon     (nlon, nlat),  &
              Ocean%Global%lat     (nlon, nlat),  &    
              Ocean%Global%mask    (nlon, nlat), &
              netr_mask            (nlon, nlat))

   allocate (rmask(nlon,nlat), geo_lont(nlon,nlat), geo_latt(nlon,nlat), &
             geo_lonv(1:nlon+1,1:nlat+1), geo_latv(1:nlon+1,1:nlat+1) )

!--- domain decompsition -----------------------------------------------
   npes = mpp_npes()

!---- domain decomposition ----

    if( layout(1).EQ.0 .AND. layout(2).EQ.0 ) &
         call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
    if( layout(1).NE.0 .AND. layout(2).EQ.0 )layout(2) = mpp_npes()/layout(1)
    if( layout(1).EQ.0 .AND. layout(2).NE.0 )layout(1) = mpp_npes()/layout(2)

    domainname = 'AMIP Ocean'
    call mpp_define_domains ( (/1,nlon,1,nlat/), layout, Ocean%Domain, name=domainname )
    call set_domain(Ocean%domain)

!----- grid information -----
    select case (grid_version)
    case(VERSION_0)
       call read_data(grid_file, "geolon_t",      geo_lont, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_t",      geo_latt, no_domain=.TRUE. )
       call read_data(grid_file, "geolon_vert_t", geo_lonv, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_vert_t", geo_latv, no_domain=.TRUE. )
       call read_data(grid_file, "wet",      rmask,     no_domain=.TRUE.)
    case(VERSION_1)
       allocate (x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4) ) 
       call read_data(grid_file, "x_T", geo_lont, no_domain=.TRUE. )
       call read_data(grid_file, "y_T", geo_latt, no_domain=.TRUE. )
       call read_data(grid_file, "x_vert_T", x_vert_t, no_domain=.TRUE.)
       call read_data(grid_file, "y_vert_T", y_vert_t, no_domain=.TRUE. )
       geo_lonv(1:nlon,1:nlat) = x_vert_t(1:nlon,1:nlat,1)
       geo_lonv(nlon+1,1:nlat) = x_vert_t(nlon,1:nlat,2)
       geo_lonv(1:nlon,nlat+1) = x_vert_t(1:nlon,nlat,4)
       geo_lonv(nlon+1,nlat+1) = x_vert_t(nlon,nlat,3)
       geo_latv(1:nlon,1:nlat) = y_vert_t(1:nlon,1:nlat,1)
       geo_latv(nlon+1,1:nlat) = y_vert_t(nlon,1:nlat,2)
       geo_latv(1:nlon,nlat+1) = y_vert_t(1:nlon,nlat,4)
       geo_latv(nlon+1,nlat+1) = y_vert_t(nlon,nlat,3)
       deallocate(x_vert_t, y_vert_t)
       call read_data(grid_file, "wet",      rmask,     no_domain=.TRUE.)
    case(VERSION_2)
       call get_mosaic_tile_grid(tile_file, ocean_mosaic, Ocean%Domain )
       allocate(tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
       allocate(garea(nlon, nlat))
       call read_data(tile_file, "x", tmpx, no_domain=.TRUE.)
       call read_data(tile_file, "y", tmpy, no_domain=.TRUE.)
       do j = 1, nlat
          do i = 1, nlon
             geo_lont(i,j) = tmpx(i*2,j*2)
             geo_latt(i,j) = tmpy(i*2,j*2)
          end do
       end do
       do j = 1, nlat+1
          do i = 1, nlon+1
             geo_lonv(i,j) = tmpx(i*2-1,j*2-1)
             geo_latv(i,j) = tmpy(i*2-1,j*2-1)
          end do
       end do

       call calc_mosaic_grid_area(geo_lonv*pi/180., geo_latv*pi/180., garea )
       garea = garea/(4*PI*RADIUS*RADIUS)  ! scale the earth are to be 1
       call field_size(grid_file, "aXo_file", siz)
       nfile_axo = siz(2)
       rmask = 0.0
       do n = 1, nfile_axo
          call read_data(grid_file, "aXo_file", axo_file, level=n)
          axo_file = 'INPUT/'//trim(axo_file)
          nxgrid = get_mosaic_xgrid_size(axo_file)
          allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
          call get_mosaic_xgrid(aXo_file, i1, j1, i2, j2, xgrid_area)
          do m = 1, nxgrid
             i = i2(m); j = j2(m)
             rmask(i,j) = rmask(i,j) + xgrid_area(m)
          end do
          deallocate(i1, j1, i2, j2, xgrid_area)
       end do
       rmask = rmask/garea    

       deallocate(tmpx, tmpy, garea)
    end select

!--- for conservation interpolation, the grid should be rectangular ----
   if(trim(interp_method) == "conservative" ) then
      err_mesg = 'Bilinear interpolation must be used for a tripolar grid'
      do i=1,nlon+1
        if(any(geo_lonv(i,1) /= geo_lonv(i,:)))  &
            call error_mesg ('ocean_model_init',err_mesg,FATAL)
      enddo
      do j=1,nlat+1
         if(any(geo_latv(1,j) /= geo_latv(:,j)))  &
            call error_mesg ('ocean_model_init',err_mesg,FATAL)
      enddo
   endif

!--- define the ice data -----------------------------------------------
   Ocean%Global%mask = .false.
   where (rmask > 0) Ocean%Global%mask = .true.
   Ocean%Global%lon_bnd = geo_lonv(:,1)*pi/180.
   Ocean%Global%lat_bnd = geo_latv(1,:)*pi/180.
   Ocean%Global%lon = geo_lont*pi/180.
   Ocean%Global%lat = geo_latt*pi/180.

!--- release the memory ------------------------------------------------
deallocate(geo_lonv, geo_latv, geo_lont, geo_latt, rmask )

!----- write (to standard output?) domain decomposition -----

     if (allocated(xextent))  deallocate ( xextent, yextent )

     if ( mpp_pe() == mpp_root_pe() ) then
          allocate ( xextent(layout(1)), yextent(layout(2)) )
          call compute_extent ( Ocean%Domain, layout, xextent, yextent )
          write (stdout(),100)
          write (stdout(),110) xextent
          write (stdout(),120) yextent
      100 format ('OCEAN DATA DOMAIN DECOMPOSITION')
      110 format ('  X-AXIS = ',24i4,/,(11x,24i4))
      120 format ('  Y-AXIS = ',24i4,/,(11x,24i4))
          deallocate ( xextent, yextent )
     endif

!----- allocate for local (compute) domain ------

   call mpp_get_compute_domain ( Ocean%Domain, isd, ied, jsd, jed )

   allocate ( Ocean%Data%lon_bnd (isd:ied+1),      &
              Ocean%Data%lat_bnd (jsd:jed+1),      &
              Ocean%Data%lon (isd:ied, jsd:jed),   &
              Ocean%Data%lat (isd:ied, jsd:jed),   &
              Ocean%Data%mask    (isd:ied,jsd:jed) )

!  ---- set up local grid -----
   Ocean%Data%lon_bnd = Ocean%Global%lon_bnd(isd:ied+1)
   Ocean%Data%lat_bnd = Ocean%Global%lat_bnd(jsd:jed+1)
   Ocean%Data%lon = Ocean%Global%lon(isd:ied, jsd:jed)
   Ocean%Data%lat = Ocean%Global%lat(isd:ied, jsd:jed)

!------------ done domain decomposition --------------------------------
!=======================================================================
!-----------------------------------------------------------------------

   allocate ( Ocean%t_surf (isd:ied,jsd:jed), &
              Ocean%u_surf (isd:ied,jsd:jed), &
              Ocean%v_surf (isd:ied,jsd:jed), &
              Ocean%frazil (isd:ied,jsd:jed), &
              Ocean%area   (isd:ied,jsd:jed), &
              Ocean%s_surf (isd:ied,jsd:jed), &
              Ocean%sea_lev(isd:ied,jsd:jed))

     !Ocean%s_surf  = 0.0
     !Ocean%sea_lev = 0.0
!-----------------------------------------------------------------------

!---- maybe this should be on restart? -----
     ! Ocean%frazil = 0.0
     ! Ocean%area   = 1.0
! ---- initialize ocean model ocean/land mask -----


!  ---- read ocean mask from the restart file (must be present) -----
!      if ( file_exist('INPUT/ocean_model.res.nc')) then
!         if (mpp_pe() == mpp_root_pe()) call error_mesg ('ocean_model_mod', &
!           'Reading NetCDF formatted restart file: INPUT/ocean_model.res.nc', NOTE)
!        netr_mask = 0.0
!        call read_data('INPUT/ocean_model.res.nc', 'mask', netr_mask, no_domain=.TRUE.)
!        where(netr_mask .GT. 0.0)
!           Ocean%Global%mask = .true.
!        endwhere
!        deallocate(netr_mask)
!     else
!        if ( file_exist('INPUT/ocean_model.res')) then
!           if (mpp_pe() == mpp_root_pe()) call error_mesg ('ocean_model_mod', &
!                'Reading native formatted restart file.', NOTE)
!
!           unit = open_restart_file ('INPUT/ocean_model.res', 'read')
!
!        ---- read global field ----
!           read ( unit ) Ocean%Global%mask
!
!20         call close_file (unit)
!
!         endif
!      endif

!   ------ define Data masks ------

   !allocate (global_mask (nlon, nlat) )

   !global_mask = Ocean%Global%mask

   !Ocean%Data %mask    = global_mask    (isd:ied,jsd:jed)
    Ocean%Data %mask    = Ocean%Global%mask(isd:ied,jsd:jed)
    
   !deallocate ( global_mask )
!  ---- initialize other modules ----

   if(trim(interp_method) == "conservative") then
      Amip = amip_interp_new ( Ocean%Data%lon_bnd,            &
                               Ocean%Data%lat_bnd,            &
                               Ocean%Data%mask,               &
                               interp_method = interp_method, &
                               use_climo=use_climo_sst,       &
                               use_annual=use_annual_sst )
   else if(trim(interp_method) == "bilinear") then
      Amip = amip_interp_new ( Ocean%Data%lon,                &
                               Ocean%Data%lat,                &
                               Ocean%Data%mask,               &
                               interp_method = interp_method, &
                               use_climo=use_climo_sst,       &
                               use_annual=use_annual_sst )
   else
      call error_mesg('ice_model_init', 'interp_method should be conservative or bilinear', &
                      FATAL)
   endif

   call get_amip_sst ( Ocean%Time, Amip, Ocean%t_surf )

#ifdef SCM
!--- for single column model -------------------------------------!
!--- initialize surface temperature to observed value ------------!  
       if (do_specified_tskin) then     
          Ocean%t_surf = TSKIN
       end if
!-----------------------------------------------------------------!
#endif

!-----------------------------------------------------------------------
!----- set the initial state -------------

   call set_ocean_model_state (Ocean)

!-----------------------------------------------------------------------
   allocate(Ocean_state)
!------------------------------------------

   module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine ocean_model_init

!#######################################################################

 subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type),           intent(inout) :: Ocean_sfc
  type(ocean_state_type),            pointer       :: Ocean_state
  type(time_type),                   intent(in)    :: Time
!   This subroutine terminates the model run, saving the ocean state in a
! restart file and deallocating any data associated with the ocean.

! Arguments: Ocean_sfc - An ocean_public_type structure that is to be
!                        deallocated upon termination.
!  (inout)   Ocean_state - A pointer to the structure containing the internal
!                          ocean state to be deallocated upon termination.
!  (in)      Time - The model time, used for writing restarts.

  !real, dimension(:, :), allocatable :: netr_mask
   integer :: unit

   if(.not.module_is_initialized) return

  !allocate(netr_mask(size(Ocean_sfc%Global%mask, 1), size(Ocean_sfc%Global%mask, 2)))
  !netr_mask = 0.0
  !where ( Ocean_sfc%Global%mask)
  !   netr_mask = 1.0
  !endwhere

  !if( do_netcdf_restart ) then
  !   if (mpp_pe() == mpp_root_pe()) call error_mesg ('ocean_model_mod', &
  !        'Writing NetCDF formatted restart file: RESTART/ocean_model.res.nc', NOTE)
  !   call write_data('RESTART/ocean_model.res.nc', 'mask', netr_mask, no_domain=.true.)
  !else
  !   call set_domain(Ocean_sfc%Domain)
  !   if (mpp_pe() == mpp_root_pe()) call error_mesg ('ocean_model_mod', &
  !        'Writing native formatted restart file.', NOTE)
  !  unit = open_restart_file ('RESTART/ocean_model.res', 'write')
  !  call write_data ( unit, Ocean_sfc%Data%mask )

!    ---- write out model fields ----
!    ---- set domain for global i/o ----

  !  call close_file (unit)
  !endif
  !deallocate(netr_mask)
  call amip_interp_del(Amip)
  module_is_initialized = .false.
     
 end subroutine ocean_model_end

!#######################################################################
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
  subroutine ocean_model_restart(Ocean_state, timestamp)
     type(ocean_state_type),    pointer     :: Ocean_state
     character(len=*), intent(in), optional :: timestamp

    call error_mesg('ocean_model_restart(ocean_model_mod)', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

!#######################################################################
! dummy interface for ESM coupler
subroutine ocean_model_init_sfc(Ocean_state, Ocean)

type(ocean_state_type), pointer          :: Ocean_state
type(ocean_public_type), intent(in)      :: Ocean

return
end subroutine ocean_model_init_sfc

!#######################################################################
subroutine ocean_model_flux_init(Ocean_state)
type(ocean_state_type), pointer       :: Ocean_state

return
end subroutine ocean_model_flux_init
!#######################################################################

 subroutine compute_extent (Domain, layout, xsizelist, ysizelist) 
 type (domain2D), intent(in) :: Domain
 integer, intent(in) :: layout(2)
 integer, intent(out), optional :: xsizelist(:), ysizelist(:)
 integer, dimension(0:layout(1)*layout(2)-1) :: xsize, ysize
 integer :: i, j, xlist(layout(1)), ylist(layout(2))
 type (domain1D) :: Xdom, Ydom

   call mpp_get_compute_domains   ( Domain, xsize=xsize, ysize=ysize )
   call mpp_get_domain_components ( Domain, Xdom, Ydom )
   call mpp_get_pelist ( Xdom, xlist ) 
   call mpp_get_pelist ( Ydom, ylist ) 

     do i = 1, layout(1)
       xsizelist(i) = xsize(xlist(i))
     enddo

     do j = 1, layout(2)
       ysizelist(j) = ysize(ylist(j))
     enddo

 end subroutine compute_extent

!#######################################################################
! dummy routine

!
subroutine ocean_stock_pe(Ocean_state, index, value, time_index)
  type(ocean_state_type),pointer     :: Ocean_state
  integer,               intent(in)  :: index
  real,                  intent(out) :: value
  integer, optional,     intent(in)  :: time_index

  value = 0.0
  if (.not.associated(Ocean_state)) return
  if (.not.Ocean_state%is_ocean_pe) return
  if(.not.stock_warning_issued) then
     call error_mesg('ocean_stock_pe','Stocks not yet implemented. Returning zero.',NOTE)
     stock_warning_issued = .true.
  endif

  end subroutine ocean_stock_pe


subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc
  
  array2D(isc:,jsc:) = 0.0
  
end subroutine ocean_model_data2D_get

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value

  value = 0.0

end subroutine ocean_model_data1D_get


!#######################################################################

subroutine ice_ocn_bnd_type_chksum(id, timestep, iobt)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_ocean_boundary_type), intent(in) :: iobt
 integer ::   n,m, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_land_boundary_type):: ', id, timestep
    write(outunit,100) 'iobt%u_flux         ', mpp_chksum( iobt%u_flux         )
    write(outunit,100) 'iobt%v_flux         ', mpp_chksum( iobt%v_flux         )
    write(outunit,100) 'iobt%t_flux         ', mpp_chksum( iobt%t_flux         )
    write(outunit,100) 'iobt%q_flux         ', mpp_chksum( iobt%q_flux         )
    write(outunit,100) 'iobt%salt_flux      ', mpp_chksum( iobt%salt_flux      )
    write(outunit,100) 'iobt%lw_flux        ', mpp_chksum( iobt%lw_flux        )
    write(outunit,100) 'iobt%sw_flux_vis_dir', mpp_chksum( iobt%sw_flux_vis_dir)
    write(outunit,100) 'iobt%sw_flux_vis_dif', mpp_chksum( iobt%sw_flux_vis_dif)
    write(outunit,100) 'iobt%sw_flux_nir_dir', mpp_chksum( iobt%sw_flux_nir_dir)
    write(outunit,100) 'iobt%sw_flux_nir_dif', mpp_chksum( iobt%sw_flux_nir_dif)
    write(outunit,100) 'iobt%lprec          ', mpp_chksum( iobt%lprec          )
    write(outunit,100) 'iobt%fprec          ', mpp_chksum( iobt%fprec          )
    write(outunit,100) 'iobt%runoff         ', mpp_chksum( iobt%runoff         )
    write(outunit,100) 'iobt%calving        ', mpp_chksum( iobt%calving        )
    write(outunit,100) 'iobt%p              ', mpp_chksum( iobt%p              )

100 FORMAT("CHECKSUM::",A32," = ",Z20)
    do n = 1, iobt%fluxes%num_bcs  !{
       do m = 1, iobt%fluxes%bc(n)%num_fields  !{
          write(outunit,101) 'iobt%',trim(iobt%fluxes%bc(n)%name), &
               trim(iobt%fluxes%bc(n)%field(m)%name), &
               mpp_chksum(iobt%fluxes%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_ocn_bnd_type_chksum

subroutine ocean_public_type_chksum(id, timestep, ocn)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_public_type), intent(in) :: ocn
 integer ::   n,m, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(ocean_type):: ', id, timestep
    write(outunit,100) 'ocean%t_surf   ',mpp_chksum(ocn%t_surf )
    write(outunit,100) 'ocean%s_surf   ',mpp_chksum(ocn%s_surf )
    write(outunit,100) 'ocean%u_surf   ',mpp_chksum(ocn%u_surf )
    write(outunit,100) 'ocean%v_surf   ',mpp_chksum(ocn%v_surf )
    write(outunit,100) 'ocean%sea_lev  ',mpp_chksum(ocn%sea_lev)
    write(outunit,100) 'ocean%frazil   ',mpp_chksum(ocn%frazil )

    do n = 1, ocn%fields%num_bcs  !{
       do m = 1, ocn%fields%bc(n)%num_fields  !{
          write(outunit,101) 'ocean%',trim(ocn%fields%bc(n)%name), &
               trim(ocn%fields%bc(n)%field(m)%name), &
               mpp_chksum(ocn%fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)


100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine ocean_public_type_chksum

end module ocean_model_mod
