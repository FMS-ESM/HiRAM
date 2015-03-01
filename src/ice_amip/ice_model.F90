module ice_model_mod

use   ice_albedo_mod, only:  ice_albedo_init, ice_albedo

use ocean_albedo_mod, only:  compute_ocean_albedo_new

use  ocean_rough_mod, only:  compute_ocean_roughness, fixed_ocean_roughness

use          fms_mod, only: file_exist, open_restart_file, close_file, &
                            mpp_pe, mpp_root_pe, mpp_npes, write_version_number, stdlog,   &
                            error_mesg, FATAL, check_nml_error, read_data, write_data,     &
                            NOTE, WARNING, field_exist, field_size, get_mosaic_tile_grid, stdout

use fms_io_mod,       only: save_restart, register_restart_field, restart_file_type, &
                            restore_state, set_domain, nullify_domain, query_initialized, &
                            get_restart_io_mode

use mpp_mod,          only: mpp_chksum

#ifdef INTERNAL_FILE_NML
use          mpp_mod, only: input_nml_file
#else
use          fms_mod, only: open_namelist_file
#endif

use    constants_mod, only: hlv, hlf, tfreeze, pi, radius

use  mpp_domains_mod, only: domain1d, domain2d, mpp_define_domains, mpp_get_compute_domain, &
                            mpp_get_compute_domains, mpp_get_domain_components, mpp_get_pelist, &
                            mpp_define_layout, mpp_define_io_domain

use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data

use time_manager_mod, only: time_type, operator(+)

use mosaic_mod,       only: get_mosaic_ntiles, get_mosaic_grid_sizes, calc_mosaic_grid_area
use mosaic_mod,       only: get_mosaic_xgrid_size, get_mosaic_xgrid

use  amip_interp_mod, only: amip_interp_type, amip_interp_new, amip_interp_del, get_amip_ice
use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
implicit none
private

public :: ice_data_type, ocean_ice_boundary_type,               &
          atmos_ice_boundary_type, land_ice_boundary_type,      &
          ice_model_init, ice_model_end, update_ice_model_fast, &
          update_ice_model_slow_up, update_ice_model_slow_dn,   &
          ice_stock_pe, cell_area, ice_model_restart,           &
          ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum, &
          lnd_ice_bnd_type_chksum, ice_data_type_chksum

type ice_data_type
  type(domain2d)                        :: Domain

   logical                              :: pe

   real,    pointer, dimension(:)       :: glon_bnd =>NULL(), &
                                           glat_bnd =>NULL(), &
                                           lon_bnd =>NULL() , &
                                           lat_bnd =>NULL()

   real,    pointer, dimension(:,:)     :: glon =>NULL(), &
                                           glat =>NULL(), &
                                           lon =>NULL(), &
                                           lat =>NULL()

   logical, pointer, dimension(:,:)     :: gmask =>NULL(), &
                                           mask =>NULL()

   logical, pointer, dimension(:,:,:)   :: ice_mask =>NULL()

   real,    pointer, dimension(:,:,:,:) :: temp =>NULL()

   real,    pointer, dimension(:,:,:)   :: part_size =>NULL(), &
                                           t_surf =>NULL(), &
                                           albedo =>NULL(), &
                                           albedo_vis_dir =>NULL(), &
                                           albedo_nir_dir =>NULL(), &
                                           albedo_vis_dif =>NULL(), &
                                           albedo_nir_dif =>NULL(), &
                                           rough_mom =>NULL(),&
                                           rough_heat =>NULL(), &
                                           rough_moist =>NULL(),  &
                                           frazil =>NULL(),  &
                                           u_surf =>NULL(),  &
                                           v_surf =>NULL()

   real,    pointer, dimension(:,:,:)   :: flux_u_bot =>NULL(), &
                                           flux_v_bot =>NULL(), &
                                           flux_t_bot =>NULL(),   &
                                           flux_q_bot =>NULL(), &
                                           flux_lh_bot =>NULL(), &
                                           flux_sw_bot =>NULL(), &
                                           flux_sw_vis_bot =>NULL(), &
                                           flux_sw_dir_bot =>NULL(), &
                                           flux_sw_dif_bot =>NULL(), &
                                           flux_sw_vis_dir_bot =>NULL(), &
                                           flux_sw_vis_dif_bot =>NULL(), &
                                           flux_sw_nir_dir_bot =>NULL(), &
                                           flux_sw_nir_dif_bot =>NULL(), &
                                           flux_lw_bot =>NULL(), &
                                           lprec_bot =>NULL(), &
                                           fprec_bot =>NULL(), &
                                           runoff_bot =>NULL()

   real,    pointer, dimension(:,:  )   :: flux_u =>NULL(), &
                                           flux_v =>NULL(), &
                                           flux_t =>NULL(), &
                                           flux_q =>NULL(), &
                                           flux_lh =>NULL(), &
                                           flux_sw =>NULL(), &
                                           flux_sw_vis =>NULL(), &
                                           flux_sw_dir =>NULL(), &
                                           flux_sw_dif =>NULL(), &
                                           flux_sw_vis_dir =>NULL(), &
                                           flux_sw_vis_dif =>NULL(), &
                                           flux_sw_nir_dir =>NULL(), &
                                           flux_sw_nir_dif =>NULL(), &
                                           flux_lw =>NULL(), &
                                           lprec =>NULL(), &
                                           fprec =>NULL(), &
                                           p_surf =>NULL(), &
                                           runoff =>NULL(), &
                                           calving =>NULL(), &
                                           runoff_hflx =>NULL(), &
                                           calving_hflx =>NULL(), &
                                           area =>NULL(), &
                                           flux_salt =>NULL()
  logical, pointer, dimension(:,:) :: maskmap =>NULL()   ! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used. This variable is dummy and need 
                                                         ! not to be set, but it is needed to pass compilation.

   integer                              :: avg_kount

   real,    pointer, dimension(:,:,:)   :: thickness =>NULL()

   type (time_type)                     :: Time_Init, Time,  &
                                           Time_step_fast,   &
                                           Time_step_slow
   integer, dimension(3)              :: axes
   type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
   type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
   type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging

end type ice_data_type

type :: ocean_ice_boundary_type
  real, dimension(:,:),   pointer :: u =>NULL(), &
                                     v =>NULL(), &
                                     t =>NULL(), &
                                     s =>NULL(), &
                                     frazil =>NULL(), &
                                     sea_level =>NULL()
  real, dimension(:,:,:), pointer :: data =>NULL()
  integer                         :: xtype
  type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
end type ocean_ice_boundary_type

type :: atmos_ice_boundary_type
  real, dimension(:,:,:), pointer :: u_flux =>NULL(), &
                                     v_flux =>NULL(), &
                                     u_star =>NULL(), &
                                     t_flux =>NULL(), &
                                     q_flux =>NULL(), &
                                     lw_flux =>NULL(), &
                                     sw_flux_vis_dir =>NULL(), &
                                     sw_flux_vis_dif =>NULL(), &
                                     sw_flux_nir_dir =>NULL(), &
                                     sw_flux_nir_dif =>NULL(), &
                                     lprec =>NULL(), &
                                     fprec =>NULL()
  real, dimension(:,:,:), pointer :: dhdt =>NULL(), &
                                     dedt =>NULL(), &
                                     drdt =>NULL(), &
                                     coszen =>NULL(), &
                                     p =>NULL(), &
                                     data =>NULL()
  integer                         :: xtype
  type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
end type atmos_ice_boundary_type

type :: land_ice_boundary_type
  real, dimension(:,:),   pointer :: runoff =>NULL(), &
                                     calving =>NULL(), &
                                     runoff_hflx =>NULL(), &
                                     calving_hflx =>NULL()
  real, dimension(:,:,:), pointer :: data =>NULL()
  integer :: xtype
end type land_ice_boundary_type

character(len=128) :: version = '$Id: ice_model.F90,v 17.0.2.1.2.1.4.1 2010/11/15 18:34:53 bw Exp $'
character(len=128) :: tagname = '$Name: hiram_20101115_bw $'

character(len=80) :: restart_format = 'amip ice model restart format 02'
logical :: module_is_initialized = .false.
logical :: stock_warning_issued  = .false.

! id's for diagnostics
integer :: id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, &
           id_runoff, id_calving, id_evap, id_fax, id_fay, &
           id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif, &
           id_sw_nir_dir, id_sw_nir_dif
logical :: sent

!-----------------------------------------------------------------------
!
!  use_climo_ice           = use monthly climatological amip ice mask
!  use_annual_ice          = use annual climatology amip ice mask
!  temp_ice_freeze         = temperature at which sea ice melts
!  no_ice                  = run with no ice (only open water)
!  use_leads               = use fraction ice coverage (i.e., leads) if it exists
!  roughness_ice           = constant roughness for all ice
!  specified_ice_thickness = constant thickness for specified ice

integer, parameter :: num_lev  = 1
integer, parameter :: num_part = 2

real    :: diff                     = 2.092  
real    :: thickness_min            = 0.10 
real    :: specified_ice_thickness  = 2.0
real    :: temp_ice_freeze          = -1.66    ! was 271.5
real    :: roughness_ice            = 1.e-4
logical :: no_ice                   = .false.
logical :: use_leads                = .false.
logical :: use_climo_ice            = .false.
logical :: use_annual_ice           = .false.
integer, dimension(2) :: layout     = (/ 0, 0 /)
integer, dimension(2) :: io_layout  = (/ 0, 0 /)
character(len=64) :: interp_method  = "conservative" ! default conservative scheme
logical :: do_netcdf_restart        = .true.
character(len=128) :: axisname_x    = 'xt'  ! x-axis name of temperature grid
character(len=128) :: axisname_y    = 'yt'  ! y-axis name of temperature grid
character(len=128) :: axisname_xb   = 'xb'  ! x-axis bounds name of temperature grid
character(len=128) :: axisname_yb   = 'yb'  ! y-axis bounds name of temperature grid

namelist /ice_model_nml/ do_netcdf_restart, diff, thickness_min, &
                         specified_ice_thickness,                &
                         temp_ice_freeze, roughness_ice,         &
                         use_climo_ice, use_annual_ice,          &
                         no_ice, use_leads, layout, interp_method,  &
                         axisname_x, axisname_y, axisname_xb, axisname_yb, &
                         io_layout

real, parameter :: latent = hlv + hlf
type(amip_interp_type), save :: Amip
real, allocatable, dimension(:,:) ::  cell_area  ! grid cell area; sphere frac.

integer :: id_restart_albedo
integer :: mlon, mlat, npart ! global grid size
type(restart_file_type), save :: Ice_restart

contains

!=============================================================================================
  subroutine ice_model_init( Ice, Time_Init, Time, Time_step_fast, Time_step_slow )
    type(ice_data_type), intent(inout) :: Ice
    type(time_type)    , intent(in)    :: Time_Init, Time, Time_step_fast, Time_step_slow

    real, allocatable, dimension(:,:)   :: geo_lonv, geo_latv, rmask
    real, allocatable, dimension(:,:)   :: geo_lont, geo_latt
    real, allocatable, dimension(:,:,:) :: x_vert_T, y_vert_T
    real, allocatable, dimension(:)     :: glonb, glatb
    real, allocatable, dimension(:)     :: xb, yb ! 1d global grid for diag_mgr
    real, allocatable, dimension(:,:)   :: tmpx, tmpy, garea
    real, allocatable, dimension(:)     :: xgrid_area(:)
    integer, allocatable, dimension(:)  :: i1, j1, i2, j2
    integer                             :: io, ierr, unit, siz(4)
    integer                             :: nlon, nlat, is, ie, js, je, i, j, k
    character(len=80)                   :: control
    character(len=80)                   :: domainname
    character(len=256)                  :: err_mesg
    character(len=64)                   :: fname = 'INPUT/ice_model.res.nc'
    character(len=64)                   :: lvltag
    character(len=256)                  :: grid_file='INPUT/grid_spec.nc'
    character(len=256)                  :: ocean_mosaic, tile_file
    character(len=256)                  :: axo_file      ! atmosXocean exchange grid file
    integer                             :: nx(1), ny(1)
    integer                             :: ntiles, nfile_axo, nxgrid, n, m
    integer                             :: grid_version
    integer, parameter                  :: VERSION_0 = 0  ! grid file with field geolon_t
    integer, parameter                  :: VERSION_1 = 1  ! grid file with field x_T
    integer, parameter                  :: VERSION_2 = 2  ! mosaic file

    if(module_is_initialized) return

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ice_model_nml, iostat=io)
    ierr = check_nml_error(io, 'ice_model_nml')
#else
    if ( file_exist( 'input.nml' ) ) then
       unit = open_namelist_file ( )
       ierr = 1
       do while ( ierr /= 0 )
          read ( unit,  nml = ice_model_nml, iostat = io, end = 10 )
          ierr = check_nml_error ( io, 'ice_model_nml' )
       enddo
10     continue
       call close_file (unit)
    endif
#endif

    call get_restart_io_mode(do_netcdf_restart)

    call write_version_number (version, tagname)
    if ( mpp_pe() == mpp_root_pe() ) then
       write (stdlog(), nml=ice_model_nml)
    endif

   !if (num_part /= 2) call error_mesg ('ice_model_init','this version must have num_part = 2', FATAL)
   !if (num_lev  /= 1) call error_mesg ('ice_model_init','this version must have num_lev = 1', FATAL)

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
       if(ntiles .NE. 1) call error_mesg('ice_model_init', ' ntiles should be 1 for ocean mosaic, contact developer', FATAL)
       call get_mosaic_grid_sizes( ocean_mosaic, nx, ny)
       nlon = nx(1)
       nlat = ny(1)
    else
       call error_mesg('ice_model_init','x_T, geolon_t, ocn_mosaic_file does not exist in file '//trim(grid_file), FATAL )
    end if

    if( nlon .LE. 0 .or. nlat .LE. 0 ) call error_mesg('ice_model_init', 'nlon and nlat should be a positive integer.', FATAL)

    !----- set up global storage and local storage -----
    allocate(Ice%gmask(nlon,nlat),    rmask(nlon,nlat) )
    allocate(geo_lonv(nlon+1,nlat+1), geo_latv(nlon+1,nlat+1) )
    allocate(geo_lont(nlon,  nlat),   geo_latt(nlon,  nlat)   )
    allocate(Ice%glon_bnd(nlon+1),    Ice%glat_bnd(nlat+1)    )
    allocate(Ice%glon(nlon,nlat),     Ice%glat(nlon,nlat)     )
    allocate(glonb(nlon+1), glatb(nlat+1), xb(nlon+1), yb(nlat+1) )

    !-------------------- domain decomposition -----------------------------
    if( layout(1).EQ.0 .AND. layout(2).EQ.0 ) &
         call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
    if( layout(1).NE.0 .AND. layout(2).EQ.0 )layout(2) = mpp_npes()/layout(1)
    if( layout(1).EQ.0 .AND. layout(2).NE.0 )layout(1) = mpp_npes()/layout(2)
    domainname = 'AMIP Ice'
    call mpp_define_domains( (/1,nlon,1,nlat/), layout, Ice%Domain, name=domainname )
    call mpp_define_io_domain (Ice%Domain, io_layout)
    call set_domain (Ice%Domain)
    call mpp_get_compute_domain( Ice%Domain, is, ie, js, je )

    !---------------- read ice cell areas from grid_spec.nc or ----------------
    !---------------- calculate the area for mosaic grid file  ----------------
    allocate (cell_area(is:ie, js:je))
    cell_area = 0.0
    select case (grid_version)
    case(VERSION_0)
       call read_data(grid_file, "geolon_t",      geo_lont, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_t",      geo_latt, no_domain=.TRUE. )
       call read_data(grid_file, "geolon_vert_t", geo_lonv, no_domain=.TRUE. )
       call read_data(grid_file, "geolat_vert_t", geo_latv, no_domain=.TRUE. )
       call read_data(grid_file, "wet",      rmask,     no_domain=.TRUE.)
       call read_data(grid_file, 'AREA_OCN', cell_area, Ice%Domain)
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
       call read_data(grid_file, 'AREA_OCN', cell_area, Ice%Domain)
    case(VERSION_2)
       call get_mosaic_tile_grid(tile_file, ocean_mosaic, Ice%Domain )
       allocate(tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
       allocate(garea(nlon,nlat))
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
       cell_area(is:ie,js:je) = garea(is:ie,js:je)  ! zero out the area of the cell on land, maybe not needed  ! kerr
       do j = js, je
          do i = is, ie
             if(rmask(i,j) == 0.0) cell_area(i,j) = 0.0
          end do
       end do
       deallocate(tmpx, tmpy, garea)
    end select

    !--- xb and yb is for diagnostics --------------------------------------
    xb = sum(geo_lonv,2)/(nlat+1)
    yb = sum(geo_latv,1)/(nlon+1)

    !--- for conservation interpolation, the grid should be rectangular ----
    if(trim(interp_method) == "conservative" ) then
       err_mesg = 'Bilinear interpolation must be used for a tripolar grid'
       do i=1,nlon+1
          if(any(geo_lonv(i,1) /= geo_lonv(i,:)))  &
               call error_mesg ('ice_model_init',err_mesg,FATAL)
       enddo
       do j=1,nlat+1
          if(any(geo_latv(1,j) /= geo_latv(:,j)))  &
               call error_mesg ('ice_model_init',err_mesg,FATAL)
       enddo
    endif
    !--- define the ice data -----------------------------------------------
    Ice%gmask = .false.
    where ( rmask > 0 ) Ice%gmask = .true.

    glonb = geo_lonv(:,1)*pi/180.
    glatb = geo_latv(1,:)*pi/180.
    Ice%glon_bnd = glonb
    Ice%glat_bnd = glatb

    Ice % glon = geo_lont*pi/180.
    Ice % glat = geo_latt*pi/180.

    !-----------------------------------------------------------------------

    allocate ( Ice%lon_bnd     (is:ie+1)                        , &
               Ice%lat_bnd     (js:je+1)                        , &
               Ice%lon         (is:ie, js:je)                   , &
               Ice%lat         (is:ie, js:je)                   , &
               Ice%ice_mask    (is:ie, js:je, num_part)         , &
               Ice%temp        (is:ie, js:je, num_part, num_lev), &
               Ice%part_size   (is:ie, js:je, num_part)         , &
               Ice%albedo      (is:ie, js:je, num_part)         , &
            Ice%albedo_vis_dir (is:ie, js:je, num_part)         , &
            Ice%albedo_nir_dir (is:ie, js:je, num_part)         , &
            Ice%albedo_vis_dif (is:ie, js:je, num_part)         , &
            Ice%albedo_nir_dif (is:ie, js:je, num_part)         , &
               Ice%rough_mom   (is:ie, js:je, num_part)         , &
               Ice%rough_heat  (is:ie, js:je, num_part)         , &
               Ice%rough_moist (is:ie, js:je, num_part)         , &
               Ice%u_surf      (is:ie, js:je, num_part)         , &
               Ice%v_surf      (is:ie, js:je, num_part)         , &
               Ice%thickness   (is:ie, js:je, num_part)         , &
               Ice%mask        (is:ie, js:je) )

    Ice%t_surf => Ice%temp (:,:,:,1)

    Ice%lon_bnd    = Ice%glon_bnd(is:ie+1)
    Ice%lat_bnd    = Ice%glat_bnd(js:je+1)
    Ice%lon        = Ice%glon(is:ie, js:je)
    Ice%lat        = Ice%glat(is:ie, js:je)
    Ice%mask       = Ice%gmask(is:ie, js:je)
    Ice%Time           = Time
    Ice%Time_init      = Time_init
    Ice%Time_step_fast = Time_step_fast
    Ice%Time_step_slow = Time_step_slow
    Ice%avg_kount = 0

    allocate ( Ice%flux_u_bot  (is:ie, js:je, num_part) , &
               Ice%flux_v_bot  (is:ie, js:je, num_part) , &
               Ice%flux_t_bot  (is:ie, js:je, num_part) , &
               Ice%flux_q_bot  (is:ie, js:je, num_part) , &
               Ice%flux_lh_bot (is:ie, js:je, num_part) , &
               Ice%flux_sw_bot (is:ie, js:je, num_part) , &
        Ice%flux_sw_vis_bot    (is:ie, js:je, num_part) , &
        Ice%flux_sw_dir_bot    (is:ie, js:je, num_part) , &
        Ice%flux_sw_dif_bot    (is:ie, js:je, num_part) , &
        Ice%flux_sw_vis_dir_bot(is:ie, js:je, num_part) , &
        Ice%flux_sw_vis_dif_bot(is:ie, js:je, num_part) , &
        Ice%flux_sw_nir_dir_bot(is:ie, js:je, num_part) , &
        Ice%flux_sw_nir_dif_bot(is:ie, js:je, num_part) , &
               Ice%flux_lw_bot (is:ie, js:je, num_part) , &
               Ice%lprec_bot   (is:ie, js:je, num_part) , &
               Ice%fprec_bot   (is:ie, js:je, num_part) , &
               Ice%runoff_bot  (is:ie, js:je, num_part) , &
               Ice%frazil      (is:ie, js:je, num_part)   )

    allocate ( Ice%flux_u    (is:ie, js:je) , &
               Ice%flux_v    (is:ie, js:je) , &
               Ice%flux_t    (is:ie, js:je) , &
               Ice%flux_q    (is:ie, js:je) , &
               Ice%flux_lh   (is:ie, js:je) , &
               Ice%flux_sw   (is:ie, js:je) , &
         Ice%flux_sw_vis     (is:ie, js:je) , &
         Ice%flux_sw_dir     (is:ie, js:je) , &
         Ice%flux_sw_dif     (is:ie, js:je) , &
         Ice%flux_sw_vis_dir (is:ie, js:je) , &
         Ice%flux_sw_vis_dif (is:ie, js:je) , &
         Ice%flux_sw_nir_dir (is:ie, js:je) , &
         Ice%flux_sw_nir_dif (is:ie, js:je) , &
               Ice%flux_lw   (is:ie, js:je) , &
               Ice%lprec     (is:ie, js:je) , &
               Ice%fprec     (is:ie, js:je) , &
               Ice%p_surf    (is:ie, js:je) , &
               Ice%runoff    (is:ie, js:je) , &
               Ice%calving   (is:ie, js:je) , &
             Ice%runoff_hflx (is:ie, js:je) , &
             Ice%calving_hflx(is:ie, js:je) , &
             Ice%area        (is:ie, js:je) , &
               Ice%flux_salt (is:ie, js:je)   )

Ice%area = cell_area  * 4*PI*RADIUS*RADIUS

    !-----------------------------------------------------------------------
    !  -------- read restart --------
if (do_netcdf_restart) call ice_register_restart(Ice, 'ice_model.res.nc')

if (file_exist('INPUT/ice_model.res.nc') ) then
  !if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
  !         'Reading NetCDF formatted restart file: INPUT/ice_model.res.nc', NOTE)
   call restore_state(Ice_restart)

   if (.not. query_initialized(Ice_restart, id_restart_albedo)) then
      if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
                'Initializing diffuse and direct streams to albedo', NOTE)
    ! initialize ocean values - ice values initialized below
      Ice%albedo_vis_dir(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_nir_dir(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_vis_dif(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_nir_dif(:,:,1) = Ice%albedo(:,:,1)
   endif

else
    if (file_exist('INPUT/ice_model.res')) then
       if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
            'Reading native formatted restart file.', NOTE)

       unit = open_restart_file ('INPUT/ice_model.res', 'read')

       read  (unit) control

       ! must use correct restart version with native format
       if (trim(control) /= trim(restart_format)) call error_mesg &
            ('ice_model_init', 'invalid restart format', FATAL)

       read  (unit) mlon, mlat, npart

       !     ---- restart resolution must be consistent with input args ---
       if (mlon /= nlon .or. mlat /= nlat .or. npart /= 2)  &
            call error_mesg ('ice_model_init',           &
            'incorrect resolution on restart', FATAL)

       call read_data ( unit, Ice%part_size  )
       call read_data ( unit, Ice%temp       )
       call read_data ( unit, Ice%thickness  )
       call read_data ( unit, Ice%albedo     )

       call read_data ( unit, Ice%albedo_vis_dir )
       call read_data ( unit, Ice%albedo_nir_dir )
       call read_data ( unit, Ice%albedo_vis_dif )
       call read_data ( unit, Ice%albedo_nir_dif )

       call read_data ( unit, Ice%rough_mom  )
       call read_data ( unit, Ice%rough_heat )
       call read_data ( unit, Ice%rough_moist)
       call read_data ( unit, Ice%u_surf     )
       call read_data ( unit, Ice%v_surf     )
       call read_data ( unit, Ice%frazil     )
       call read_data ( unit, Ice%flux_u_bot )
       call read_data ( unit, Ice%flux_v_bot )

       call close_file (unit)

    else

       !     ----- no restart - no ice -----

       Ice%temp      = tfreeze + temp_ice_freeze
       Ice%thickness = 0.0
       Ice%part_size         = 0.0
       Ice%part_size (:,:,1) = 1.0
       Ice%albedo     = 0.14
     ! initialize ocean values - ice values initialized below
       Ice%albedo_vis_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_vis_dif(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dif(:,:,1) = Ice%albedo(:,:,1)
       Ice%rough_mom  = 0.0004
       Ice%rough_heat = 0.0004
       Ice%rough_moist= 0.0004
       Ice%u_surf     = 0.0
       Ice%v_surf     = 0.0
       Ice%frazil     = 0.0

       !     --- open water roughness (initially fixed) ---

       call fixed_ocean_roughness ( Ice%mask, Ice%rough_mom  (:,:,1), &
            Ice%rough_heat (:,:,1), &
            Ice%rough_moist(:,:,1)  )

    endif
endif

!! set ice partiton values to that of Ice%albedo.
   Ice%albedo_vis_dir(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_nir_dir(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_vis_dif(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_nir_dif(:,:,2) = Ice%albedo(:,:,2)

    !---- initialization of ice mask (actually where ice exists) -----
    Ice%ice_mask = .false.
    where (Ice%mask     (:,:)           .and.     &
         Ice%part_size(:,:,2) > 0.0   .and.     &
         Ice%thickness(:,:,2) >= thickness_min) &
         Ice%ice_mask(:,:,2) = .true.


    call ice_albedo_init (tfreeze)

    if(trim(interp_method) == "conservative") then
       Amip = amip_interp_new ( Ice%lon_bnd,     Ice%lat_bnd,  &
            Ice%mask(:,:),                 &
            interp_method = interp_method, &
            use_climo=use_climo_ice,       &
            use_annual=use_annual_ice     )
    else if(trim(interp_method) == "bilinear") then
       Amip = amip_interp_new ( Ice%lon,     Ice%lat,          &
            Ice%mask(:,:),                 &
            interp_method = interp_method, &
            use_climo=use_climo_ice,       &
            use_annual=use_annual_ice     )
    else
       call error_mesg('ice_model_init', 'interp_method should be conservative or bilinear', &
            FATAL)
    endif

    ! --- diagnostics ---

    call ice_diag_init (Ice, xb, yb)

    !--- release the memory ------------------------------------------------
    deallocate(geo_lonv, geo_latv, geo_lont, geo_latt, glonb, glatb, rmask, xb, yb )
    call nullify_domain()

    module_is_initialized = .true.

  end subroutine ice_model_init
!=============================================================================================
subroutine print_layout (npes, layout, Domain)
integer, intent(in) :: npes, layout(2)
type(domain2d), intent(in) :: Domain
integer, dimension(0:npes-1) :: xsize, ysize
integer :: i, j, xlist(layout(1)), ylist(layout(2))
type (domain1D) :: Xdom, Ydom

call mpp_get_compute_domains   ( Domain, xsize=xsize, ysize=ysize )
call mpp_get_domain_components ( Domain, Xdom, Ydom )
call mpp_get_pelist ( Xdom, xlist )
call mpp_get_pelist ( Ydom, ylist ) 

write (*,100)             
write (*,110) (xsize(xlist(i)),i=1,layout(1))
write (*,120) (ysize(ylist(j)),j=1,layout(2))
                               
100 format ('ICE MODEL DOMAIN DECOMPOSITION')
110 format ('  X-AXIS = ',24i4,/,(11x,24i4))
120 format ('  Y-AXIS = ',24i4,/,(11x,24i4))

end subroutine print_layout
!=============================================================================================
subroutine update_ice_model_fast( Atmos_boundary, Ice )
type(atmos_ice_boundary_type), intent(in) :: Atmos_boundary
type (ice_data_type), intent(inout) :: Ice

real, dimension(size(Ice%flux_u_bot,1),size(Ice%flux_u_bot,2),size(Ice%flux_u_bot,3)) ::  &
      ts_new, gamma, flux_i, t_dt_surf, flux_t_new, flux_q_new, flux_lw_new, flux_sw_new, &
      flux_sw_vis_new, flux_sw_dir_new, flux_sw_dif_new, &
      flux_sw_vis_dir_new, flux_sw_vis_dif_new, &
      flux_u_new, flux_v_new, lprec_new, fprec_new, flux_lh_new

!-----------------------------------------------------------------------
!
!   updates ice model on the atmospheric (fast) time step 
!   averages input quantities to be seen by the ocean
!
!    flux_u  = zonal wind stress
!    flux_v  = meridional wind stress
!    flux_sw = net shortwave radiation (down-up) 
!    flux_sw_vis = net visible shortwave radiation (down-up)
!    flux_sw_dir = net direct shortwave radiation (down-up)
!    flux_sw_dif = net diffuse shortwave radiation (down-up)
!    flux_sw_vis_dir = net visible direct shortwave radiation (down-up)
!    flux_sw_vis_dif = net visible diffuse shortwave radiation (down-up)
!    flux_sw_nir_dir = net near IR direct shortwave radiation (down-up)
!    flux_sw_nir_dif = net near IR diffuse shortwave radiation (down-up)
!    flux_lw = net longwave radiation (down-up) 
!    flux_t  = sensible heat flux
!    flux_q  = specific humidity flux
!    flux_lh = latent heat flux
!    lprec   = liquid precipitiation rate (kg/m2/s)
!    fprec   = frozen precipitiation rate (kg/m2/s)
!    coszen  = cosine of the zenith angle
!
!-----------------------------------------------------------------------
!----- set up local copies of fluxes for modification -----

flux_u_new  = Atmos_boundary%u_flux
flux_v_new  = Atmos_boundary%v_flux
flux_t_new  = Atmos_boundary%t_flux
flux_q_new  = Atmos_boundary%q_flux
flux_lh_new = Atmos_boundary%q_flux*hlv
flux_lw_new = Atmos_boundary%lw_flux
flux_sw_new = Atmos_boundary%sw_flux_vis_dir + Atmos_boundary%sw_flux_vis_dif + &
              Atmos_boundary%sw_flux_nir_dir + Atmos_boundary%sw_flux_nir_dif
flux_sw_vis_new = Atmos_boundary%sw_flux_vis_dir + Atmos_boundary%sw_flux_vis_dif
flux_sw_dir_new = Atmos_boundary%sw_flux_vis_dir + Atmos_boundary%sw_flux_nir_dir
flux_sw_dif_new = Atmos_boundary%sw_flux_vis_dif + Atmos_boundary%sw_flux_nir_dif
flux_sw_vis_dir_new = Atmos_boundary%sw_flux_vis_dir
flux_sw_vis_dif_new = Atmos_boundary%sw_flux_vis_dif
lprec_new   = Atmos_boundary%lprec
fprec_new   = Atmos_boundary%fprec

!----- implicit update of ice surface temperature -----

ts_new = 0.0

where (Ice%ice_mask)
  gamma = diff / max(Ice%thickness,thickness_min)
  flux_i = gamma * (tfreeze + temp_ice_freeze - Ice%t_surf)

  t_dt_surf = (flux_i + Atmos_boundary%lw_flux + &
               Atmos_boundary%sw_flux_vis_dir + &
               Atmos_boundary%sw_flux_vis_dif + &
               Atmos_boundary%sw_flux_nir_dir + &
               Atmos_boundary%sw_flux_nir_dif - &
               Atmos_boundary%t_flux - Atmos_boundary%q_flux*latent)           &
           / (Atmos_boundary%dhdt + Atmos_boundary%dedt*latent + Atmos_boundary%drdt + gamma)

  ts_new = Ice%t_surf + t_dt_surf
  flux_lh_new = flux_lh_new + Atmos_boundary%q_flux*hlf
endwhere

!   ----- compute new fluxes (adjusted for temp change) -----
!              (longwave up has negative sign)

where (Ice%ice_mask .and. ts_new > tfreeze )
  t_dt_surf   = tfreeze - Ice%t_surf
  Ice%t_surf  = tfreeze
  flux_t_new  = flux_t_new  + t_dt_surf * Atmos_boundary%dhdt
  flux_q_new  = flux_q_new  + t_dt_surf * Atmos_boundary%dedt
  flux_lh_new = flux_lh_new + t_dt_surf * Atmos_boundary%dedt*latent
  flux_lw_new = flux_lw_new - t_dt_surf * Atmos_boundary%drdt
endwhere

where (Ice%ice_mask .and. ts_new <= tfreeze)
  Ice%t_surf  = Ice%t_surf  + t_dt_surf
  flux_t_new  = flux_t_new  + t_dt_surf * Atmos_boundary%dhdt
  flux_q_new  = flux_q_new  + t_dt_surf * Atmos_boundary%dedt
  flux_lh_new = flux_lh_new + t_dt_surf * Atmos_boundary%dedt*latent
  flux_lw_new = flux_lw_new - t_dt_surf * Atmos_boundary%drdt
endwhere

!-----------------------------------------------------------------------
!------ update ocean/ice surface parameters --------

!  ---- over ice -----

where (Ice%ice_mask)
  Ice%rough_mom   = roughness_ice
  Ice%rough_heat  = roughness_ice
  Ice%rough_moist = roughness_ice
endwhere

call ice_albedo (Ice%ice_mask, Ice%thickness, Ice%t_surf, Ice%albedo)

!!! EIther define some or all of the additional albedoes in ice_albedo, 
!!! or define them  upon return, based on the values returned.
!call ice_albedo (Ice%ice_mask, Ice%thickness, Ice%t_surf, Ice%albedo, &
!                 Ice%albedo_vis_dir, Ice%albedo_nir_dir,   &
!               Ice%albedo_vis_dif, Ice%albedo_nir_dif)

!! FOR now, simply set all to be the same:
     Ice%albedo_vis_dir = Ice%albedo
     Ice%albedo_nir_dir = Ice%albedo
     Ice%albedo_vis_dif = Ice%albedo
     Ice%albedo_nir_dif = Ice%albedo

!  ---- over ocean -----
!  store values into ice-free partition (n=1)

call compute_ocean_roughness ( Ice%mask,   &       
            Atmos_boundary%u_star(:,:,1),  &
            Ice%rough_mom(:,:,1), Ice%rough_heat(:,:,1), Ice%rough_moist(:,:,1) )

!!! EIther define some or all of the additional albedoes in 
!!   compute_ocean_albedo, or define them  upon return, based on the 
!!   values returned.

!call compute_ocean_albedo ( Ice%mask, Atmos_boundary%coszen(:,:,1), Ice%albedo(:,:,1) )
call compute_ocean_albedo_new ( Ice%mask, Atmos_boundary%coszen(:,:,1), &
                                Ice%albedo_vis_dir(:,:,1), Ice%albedo_vis_dif(:,:,1),      &
                                Ice%albedo_nir_dir(:,:,1), Ice%albedo_nir_dif(:,:,1), Ice%lat )

!-----------------------------------------------------------------------
!----- average fluxes to be used by the ocean model -----
!-----------------------------------------------------------------------

    call sum_bottom_quantities ( Ice, flux_u_new,  flux_v_new,  &
                                      flux_t_new,  flux_q_new,  flux_lh_new,  &
                                      flux_sw_new, flux_lw_new, &
                                      lprec_new,   fprec_new,   &
                                      flux_sw_vis_new, flux_sw_dir_new,&
                                      flux_sw_dif_new,   &
                                      flux_sw_vis_dir_new,&
                                      flux_sw_vis_dif_new )

!-----------------------------------------------------------------------
!--------- advance time -----------------

Ice%Time = Ice%Time + Ice%Time_step_fast

!--------- do diagnostics here ----------

end subroutine update_ice_model_fast
!=============================================================================================

subroutine sum_bottom_quantities ( Ice, flux_u,  flux_v,  &
                                        flux_t,  flux_q,  flux_lh, &
                                        flux_sw, flux_lw, &
                                        lprec,   fprec,   &
                                        flux_sw_vis, flux_sw_dir,&
                                        flux_sw_dif,   &
                                        flux_sw_vis_dir,&
                                        flux_sw_vis_dif )

type (ice_data_type), intent(inout)  :: Ice
real, intent(in), dimension(:,:,:)   :: flux_u,  flux_v,  &
                                        flux_t,  flux_q,  flux_lh, &
                                        flux_sw, flux_lw, &
                                        lprec,   fprec,   &
                                        flux_sw_vis, flux_sw_dir,&
                                        flux_sw_dif,   &
                                        flux_sw_vis_dir,&
                                        flux_sw_vis_dif

if (Ice%avg_kount == 0) call zero_bottom_quantities (Ice)

Ice%flux_u_bot  = Ice%flux_u_bot  + flux_u
Ice%flux_v_bot  = Ice%flux_v_bot  + flux_v
Ice%flux_t_bot  = Ice%flux_t_bot  + flux_t
Ice%flux_q_bot  = Ice%flux_q_bot  + flux_q
Ice%flux_lh_bot = Ice%flux_lh_bot + flux_lh
Ice%flux_sw_bot = Ice%flux_sw_bot + flux_sw
Ice%flux_sw_vis_bot = Ice%flux_sw_vis_bot + flux_sw_vis
Ice%flux_sw_dir_bot = Ice%flux_sw_dir_bot + flux_sw_dir
Ice%flux_sw_dif_bot = Ice%flux_sw_dif_bot + flux_sw_dif
Ice%flux_sw_vis_dir_bot = Ice%flux_sw_vis_dir_bot + flux_sw_vis_dir
Ice%flux_sw_vis_dif_bot = Ice%flux_sw_vis_dif_bot + flux_sw_vis_dif
Ice%flux_sw_nir_dir_bot = Ice%flux_sw_nir_dir_bot + flux_sw_vis_dir
Ice%flux_sw_nir_dif_bot = Ice%flux_sw_nir_dif_bot + flux_sw_vis_dif
Ice%flux_lw_bot = Ice%flux_lw_bot + flux_lw
Ice%lprec_bot   = Ice%lprec_bot   + lprec
Ice%fprec_bot   = Ice%fprec_bot   + fprec

Ice%avg_kount = Ice%avg_kount + 1

end subroutine sum_bottom_quantities
!=============================================================================================
subroutine zero_bottom_quantities ( Ice )
type (ice_data_type), intent(inout) :: Ice

Ice%avg_kount = 0
Ice%flux_u_bot  = 0.0
Ice%flux_v_bot  = 0.0
Ice%flux_t_bot  = 0.0
Ice%flux_q_bot  = 0.0
Ice%flux_lh_bot = 0.0
Ice%flux_sw_bot = 0.0
Ice%flux_sw_vis_bot = 0.0
Ice%flux_sw_dir_bot = 0.0
Ice%flux_sw_dif_bot = 0.0
Ice%flux_sw_vis_dir_bot = 0.0
Ice%flux_sw_vis_dif_bot = 0.0
Ice%flux_sw_nir_dir_bot = 0.0
Ice%flux_sw_nir_dif_bot = 0.0
Ice%flux_lw_bot = 0.0
Ice%lprec_bot   = 0.0
Ice%fprec_bot   = 0.0

end subroutine zero_bottom_quantities
!=============================================================================================
subroutine update_ice_model_slow_up( Ocean_boundary, Ice )
type(ocean_ice_boundary_type), intent(in) :: Ocean_boundary
type(ice_data_type),           intent(inout) :: Ice

call ice_bottom_to_ice_top ( Ice, &
                             Ocean_boundary%t,      &
                             Ocean_boundary%frazil, &
                             Ocean_boundary%u,      &
                             Ocean_boundary%v )

end subroutine update_ice_model_slow_up
!=============================================================================================
subroutine ice_bottom_to_ice_top ( Ice, t_surf_ice_bot, frazil_ice_bot, &
                                   u_surf_ice_bot, v_surf_ice_bot  )

type (ice_data_type), intent(inout) :: Ice
real, dimension(:,:), intent(in)    :: t_surf_ice_bot,  frazil_ice_bot, &
                                       u_surf_ice_bot,  v_surf_ice_bot
!-----------------------------------------------------------------------
!                 pass ocean state through ice
!            store values into ice-free partition (n=1)

where (Ice%part_size(:,:,1) > .00001)
  Ice%temp  (:,:,1,1) = t_surf_ice_bot
  Ice%frazil(:,:,1)   = frazil_ice_bot
  Ice%u_surf(:,:,1)   = u_surf_ice_bot
  Ice%v_surf(:,:,1)   = v_surf_ice_bot
endwhere

 end subroutine ice_bottom_to_ice_top
!=============================================================================================
subroutine update_ice_model_slow_dn( Atmos_boundary, Land_boundary, Ice )
type(atmos_ice_boundary_type), intent(in   ) :: Atmos_boundary
type(land_ice_boundary_type),  intent(in   ) :: Land_boundary
type(ice_data_type),           intent(inout) :: Ice

real, dimension(size(Ice%mask,1),size(Ice%mask,2)) :: frac
real :: frac_cutoff, frac_floor

!-----------------------------------------------------------------------

!----- compute average fluxes -----

   call avg_bottom_quantities ( Ice )

!
! Flux diagnostics
!
  if (id_sh   >0) sent = send_data(id_sh,    Ice%flux_t,  Ice%Time, mask=Ice%mask)
  if (id_lh   >0) sent = send_data(id_lh,    Ice%flux_lh, Ice%Time, mask=Ice%mask)
  if (id_evap >0) sent = send_data(id_evap,  Ice%flux_q,  Ice%Time, mask=Ice%mask)
  if (id_sw   >0) sent = send_data(id_sw,    Ice%flux_sw, Ice%Time, mask=Ice%mask)
  if (id_sw_vis   >0) sent = send_data(id_sw_vis,    Ice%flux_sw_vis, Ice%Time, mask=Ice%mask)
  if (id_sw_dir   >0) sent = send_data(id_sw_dir,    Ice%flux_sw_dir, Ice%Time, mask=Ice%mask)
  if (id_sw_dif   >0) sent = send_data(id_sw_dif,    Ice%flux_sw_dif, Ice%Time, mask=Ice%mask)
  if (id_sw_vis_dir   >0) sent = send_data(id_sw_vis_dir,    Ice%flux_sw_vis_dir, Ice%Time, mask=Ice%mask)
  if (id_sw_vis_dif   >0) sent = send_data(id_sw_vis_dif,    Ice%flux_sw_vis_dif, Ice%Time, mask=Ice%mask)
  if (id_sw_vis_dir   >0) sent = send_data(id_sw_nir_dir,    Ice%flux_sw_nir_dir, Ice%Time, mask=Ice%mask)
  if (id_sw_vis_dif   >0) sent = send_data(id_sw_nir_dif,    Ice%flux_sw_nir_dif, Ice%Time, mask=Ice%mask)
  if (id_lw   >0) sent = send_data(id_lw,    Ice%flux_lw, Ice%Time, mask=Ice%mask)
  if (id_snofl>0) sent = send_data(id_snofl, Ice%fprec,   Ice%Time, mask=Ice%mask)
  if (id_rain >0) sent = send_data(id_rain,  Ice%lprec,   Ice%Time, mask=Ice%mask)

  if (id_fax  >0) sent = send_data(id_fax,   Ice%flux_u,  Ice%Time, mask=Ice%mask)
  if (id_fay  >0) sent = send_data(id_fay,   Ice%flux_v,  Ice%Time, mask=Ice%mask)

!----- quantities from land model ----

  Ice%runoff  = Land_boundary%runoff
  Ice%calving = Land_boundary%calving

  if (id_runoff >0) sent = send_data(id_runoff,  Ice%runoff,  Ice%Time, mask=Ice%mask)
  if (id_calving>0) sent = send_data(id_calving, Ice%calving, Ice%Time, mask=Ice%mask)

!-----------------------------------------------------------------------
!----- modify fluxes to be used by the ocean model -----
!-----------------------------------------------------------------------

!----- where there is ice, ocean will not feel fluxes ? -----

where (Ice%ice_mask(:,:,2))
  Ice%flux_u = 0.0
  Ice%flux_v = 0.0
endwhere

!-----------------------------------------------------------------------
!---- get the specified ice field -----

      call get_amip_ice (Ice%Time, Amip, frac)


!  --- turn off sea-ice ??? ---
   if (no_ice) frac = 0.0

!  --- set constants for determining ice fraction
   if (use_leads) then
       frac_cutoff = 1.e-6 ! machine dependent value
       frac_floor = 0.0
   else
      !--- discretize (0. or 1.) ----
       frac_cutoff = 0.5
       frac_floor = 1.0
   endif

!  --- determine which grid boxes have ice coverage ---
   where ( Ice%mask(:,:) .and. frac > frac_cutoff )
!     --- ice ---
      Ice%part_size(:,:,2) = min(max(frac_floor,frac),1.0)
      Ice%part_size(:,:,1) = 1.0 - Ice%part_size(:,:,2)
      Ice%thickness(:,:,2) = specified_ice_thickness
      Ice%ice_mask (:,:,2) = .true.
   elsewhere
!     --- no ice ---
      Ice%part_size(:,:,1) = 1.0
      Ice%part_size(:,:,2) = 0.0
      Ice%thickness(:,:,2) = 0.0
      Ice%ice_mask (:,:,2) = .false.
   endwhere

end subroutine update_ice_model_slow_dn
!=============================================================================================
subroutine avg_bottom_quantities ( Ice )
type(ice_data_type), intent(inout) :: Ice
real :: divid

!----- compute average fluxes -----

if (Ice%avg_kount == 0) call error_mesg ('avg_bottom_quantities', &
                      'no ocean model fluxes have been averaged', FATAL)

divid = 1./float(Ice%avg_kount)

Ice%flux_u_bot  = Ice%flux_u_bot  * divid
Ice%flux_v_bot  = Ice%flux_v_bot  * divid
Ice%flux_t_bot  = Ice%flux_t_bot  * divid
Ice%flux_q_bot  = Ice%flux_q_bot  * divid
Ice%flux_lh_bot = Ice%flux_lh_bot * divid
Ice%flux_sw_bot = Ice%flux_sw_bot * divid
Ice%flux_sw_vis_bot = Ice%flux_sw_vis_bot * divid
Ice%flux_sw_dir_bot = Ice%flux_sw_dir_bot * divid
Ice%flux_sw_dif_bot = Ice%flux_sw_dif_bot * divid
Ice%flux_sw_vis_dir_bot = Ice%flux_sw_vis_dir_bot * divid
Ice%flux_sw_vis_dif_bot = Ice%flux_sw_vis_dif_bot * divid
Ice%flux_sw_nir_dir_bot = Ice%flux_sw_nir_dir_bot * divid
Ice%flux_sw_nir_dif_bot = Ice%flux_sw_nir_dif_bot * divid
Ice%flux_lw_bot = Ice%flux_lw_bot * divid
Ice%lprec_bot   = Ice%lprec_bot   * divid
Ice%fprec_bot   = Ice%fprec_bot   * divid

Ice%flux_t  = all_avg( Ice%flux_t_bot , Ice%part_size )
Ice%flux_q  = all_avg( Ice%flux_q_bot , Ice%part_size )
Ice%flux_lh = all_avg( Ice%flux_lh_bot, Ice%part_size )
Ice%flux_sw = all_avg( Ice%flux_sw_bot, Ice%part_size )
Ice%flux_sw_vis = all_avg( Ice%flux_sw_vis_bot, Ice%part_size )
Ice%flux_sw_dir = all_avg( Ice%flux_sw_dir_bot, Ice%part_size )
Ice%flux_sw_dif = all_avg( Ice%flux_sw_dif_bot, Ice%part_size )
Ice%flux_sw_vis_dir = all_avg( Ice%flux_sw_vis_dir_bot, Ice%part_size )
Ice%flux_sw_vis_dif = all_avg( Ice%flux_sw_vis_dif_bot, Ice%part_size )
Ice%flux_sw_nir_dir = all_avg( Ice%flux_sw_nir_dir_bot, Ice%part_size )
Ice%flux_sw_nir_dif = all_avg( Ice%flux_sw_nir_dif_bot, Ice%part_size )
Ice%flux_lw = all_avg( Ice%flux_lw_bot, Ice%part_size )
Ice%fprec   = all_avg( Ice%fprec_bot  , Ice%part_size )
Ice%lprec   = all_avg( Ice%lprec_bot  , Ice%part_size )

Ice%flux_u  = all_avg( Ice%flux_u_bot , Ice%part_size )
Ice%flux_v  = all_avg( Ice%flux_v_bot , Ice%part_size )

!--- set count to zero and fluxes will be zeroed before the next sum

Ice%avg_kount = 0

end subroutine avg_bottom_quantities
!=============================================================================================
function all_avg(x,part)
real, dimension(:,:,:) :: x, part
real, dimension(size(x,1), size(x,2)) :: all_avg
integer :: k

if(any(shape(x) /= shape(part))) then
  call error_mesg('all_avg','input arguments "x" and "part" are not dimensioned the same',FATAL)
endif

all_avg = 0.
do k=1,size(x,3)
  all_avg = all_avg + part(:,:,k)*x(:,:,k)
enddo

return
end function all_avg
!=============================================================================================
subroutine ice_model_end(Ice)
type(ice_data_type), intent(inout) :: Ice
character(len=64) :: fname='RESTART/ice_model.res.nc'
integer :: unit, k
character(len=64) :: lvltag

if(.not.module_is_initialized) return
if( do_netcdf_restart) then
   if(mpp_pe() == mpp_root_pe() ) then
      call error_mesg ('ice_model_mod', 'Writing NetCDF formatted restart file: RESTART/ice_model.res.nc', NOTE)
   endif
   call save_restart(Ice_restart)
else

   if(mpp_pe() == mpp_root_pe() ) then
      call error_mesg ('ice_model_mod', 'Writing native formatted restart file.', NOTE)
   endif
  unit = open_restart_file ('RESTART/ice_model.res', 'write')
  if ( mpp_pe() == mpp_root_pe() ) then
    write (unit) restart_format
    write (unit) size(Ice%gmask,1), size(Ice%gmask,2), num_part
  endif

  call set_domain (Ice%Domain)
  call write_data ( unit, Ice%part_size  )
  call write_data ( unit, Ice%temp       )
  call write_data ( unit, Ice%thickness  )
  call write_data ( unit, Ice%albedo     )
! code to output the new albedos
  call write_data ( unit, Ice%albedo_vis_dir )
  call write_data ( unit, Ice%albedo_nir_dir )
  call write_data ( unit, Ice%albedo_vis_dif )
  call write_data ( unit, Ice%albedo_nir_dif )

  call write_data ( unit, Ice%rough_mom  )
  call write_data ( unit, Ice%rough_heat )
  call write_data ( unit, Ice%rough_moist)
  call write_data ( unit, Ice%u_surf     )
  call write_data ( unit, Ice%v_surf     )
  call write_data ( unit, Ice%frazil     )
  call write_data ( unit, Ice%flux_u_bot )
  call write_data ( unit, Ice%flux_v_bot )
  call close_file ( unit )
endif

deallocate(cell_area)
call amip_interp_del(Amip)

module_is_initialized = .false.

end subroutine ice_model_end

 !#######################################################################
  ! <SUBROUTINE NAME="ice_model_restart">
  ! <DESCRIPTION>
  !  dummy routine
  ! </DESCRIPTION>
  subroutine ice_model_restart(time_stamp)
    character(len=*),         intent(in), optional :: time_stamp

    call error_mesg ('ice_model_restart in ice_model_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)    


  end subroutine ice_model_restart
  ! </SUBROUTINE>

!=============================================================================================

subroutine ice_diag_init (Ice, xb, yb)
type(ice_data_type), intent(in) :: Ice
real   , intent(in) :: xb(:), yb(:)

integer :: nlon, nlat
integer :: id_xv, id_yv, id_xb, id_xt, id_yb, id_yt, axv(2), axt(2)
real, parameter :: missing = -1e34

  nlon = size(xb(:))-1
  nlat = size(yb(:))-1

! define axes

 !id_xv = diag_axis_init('xv', xb(2:nlon+1), 'degrees_E', 'X','longitude', &
 !                                 set_name='ice', Domain2=Ice%Domain )
 !id_yv = diag_axis_init('yv', yb(2:nlat+1), 'degrees_N', 'Y','latitude',  &
 !                                 set_name='ice', Domain2=Ice%Domain )
 !axv = (/ id_xv, id_yv /)

  id_xb = diag_axis_init(trim(axisname_xb), xb, 'degrees_E', 'X', 'longitude', &
                                    set_name='ice', Domain2=Ice%Domain )
  id_xt = diag_axis_init(trim(axisname_x), (xb(1:nlon)+xb(2:nlon+1))/2, 'degrees_E', 'X', &
                     'longitude',set_name='ice',edges=id_xb,Domain2=Ice%Domain)
  id_yb = diag_axis_init(trim(axisname_yb), yb, 'degrees_N', 'Y', 'latitude', &
                                   set_name='ice', Domain2=Ice%Domain )
  id_yt = diag_axis_init(trim(axisname_y), (yb(1:nlat)+yb(2:nlat+1))/2, 'degrees_N', 'Y', &
                     'latitude',set_name='ice', edges=id_yb,Domain2=Ice%Domain)
  axt  = (/ id_xt, id_yt /)

! register fields

  id_sh = register_diag_field('ice_model','SH' ,axt, Ice%Time, &
                              'sensible heat flux', 'W/m^2', &
                              missing_value=missing)
  id_lh = register_diag_field('ice_model','LH' ,axt, Ice%Time, &
                             'latent heat flux', 'W/m^2', missing_value=missing)
  id_sw = register_diag_field('ice_model','SW' ,axt, Ice%Time, &
                              'short wave heat flux', 'W/m^2',  &
                              missing_value=missing)
  id_sw_vis = register_diag_field('ice_model','SWvis' ,axt, Ice%Time, &
                              'visible short wave heat flux', 'W/m^2', &
                              missing_value=missing)
  id_sw_dir = register_diag_field('ice_model','SWdir' ,axt, Ice%Time, &
                              'direct short wave heat flux', 'W/m^2',  &
                              missing_value=missing)
  id_sw_dif = register_diag_field('ice_model','SWdif' ,axt, Ice%Time, &
                              'diffuse short wave heat flux', 'W/m^2', &
                              missing_value=missing)
  id_sw_vis_dir = register_diag_field('ice_model','SWvisdir' ,axt, Ice%Time, &
                              'vis dir short wave heat flux', 'W/m^2', &
                              missing_value=missing)
  id_sw_vis_dif = register_diag_field('ice_model','SWvisdif' ,axt, Ice%Time, &
                              'vis diff short wave heat flux', 'W/m^2',  &
                              missing_value=missing)
  id_sw_nir_dir = register_diag_field('ice_model','SWnirdir' ,axt, Ice%Time, &
                              'NIR dir short wave heat flux', 'W/m^2', &
                              missing_value=missing)
  id_sw_nir_dif = register_diag_field('ice_model','SWnirdif' ,axt, Ice%Time, &
                              'NIR diff short wave heat flux', 'W/m^2',  &
                              missing_value=missing)
  id_lw = register_diag_field('ice_model','LW' ,axt, Ice%Time, &
                              'long wave heat flux over ice', 'W/m^2', &
                              missing_value=missing)
  id_snofl = register_diag_field('ice_model','SNOWFL' ,axt, Ice%Time, &
                                 'rate of snow fall', 'kg/(m^2*s)', &
                                 missing_value=missing)
  id_rain  = register_diag_field('ice_model','RAIN' ,axt, Ice%Time, &
                                 'rate of rain fall', 'kg/(m^2*s)', &
                                 missing_value=missing)
  id_runoff= register_diag_field('ice_model','RUNOFF' ,axt, Ice%Time, &
                                 'liquid runoff', 'kg/(m^2*s)', &
                                 missing_value=missing)
  id_calving = register_diag_field('ice_model','CALVING',axt, Ice%Time, &
                                 'frozen runoff', 'kg/(m^2*s)', &
                                 missing_value=missing)
  id_evap = register_diag_field('ice_model','EVAP',axt, Ice%Time, &
                                 'evaporation', 'kg/(m^2*s)', &
                                 missing_value=missing)
! wind stress at t points
  id_fax = register_diag_field('ice_model', 'FA_X', axt, Ice%Time, &
                               'air stress on ice - x component', 'Pa', &
                               missing_value=missing)
  id_fay = register_diag_field('ice_model', 'FA_Y', axt, Ice%Time, &
                               'air stress on ice - y component', 'Pa', &
                               missing_value=missing)

end subroutine ice_diag_init

!=============================================================================================

 subroutine ice_register_restart (Ice, restart_file)
 type(ice_data_type), intent(inout) :: Ice
 character(len=*), intent(in) :: restart_file
 integer                      :: id_restart

  !id_restart = register_restart_field (Ice_restart, restart_file, 'mlon',      mlon)
  !id_restart = register_restart_field (Ice_restart, restart_file, 'mlon',      mlat)
  !id_restart = register_restart_field (Ice_restart, restart_file, 'num_part',  num_part)
   id_restart = register_restart_field (Ice_restart, restart_file, 'part_size', Ice%part_size, domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'thickness', Ice%thickness, domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'albedo',    Ice%albedo,    domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'temp_1',    Ice%temp(:,:,:,1), domain=Ice%domain)

   ! albedo streams
   id_restart_albedo = register_restart_field (Ice_restart, restart_file, 'albedo_vis_dir', Ice%albedo_vis_dir, &
                                                 domain=Ice%domain, mandatory=.false.)
   id_restart        = register_restart_field (Ice_restart, restart_file, 'albedo_nir_dir', Ice%albedo_nir_dir, &
                                                 domain=Ice%domain, mandatory=.false.)
   id_restart        = register_restart_field (Ice_restart, restart_file, 'albedo_vis_dif', Ice%albedo_vis_dif, &
                                                 domain=Ice%domain, mandatory=.false.)
   id_restart        = register_restart_field (Ice_restart, restart_file, 'albedo_nir_dif', Ice%albedo_nir_dif, &
                                                 domain=Ice%domain, mandatory=.false.)

   id_restart = register_restart_field (Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,   domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,  domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'rough_moist', Ice%rough_moist, domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'u_surf',      Ice%u_surf,      domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'v_surf',      Ice%v_surf,      domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'frazil',      Ice%frazil,      domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'flux_u_bot',  Ice%flux_u_bot,  domain=Ice%domain)
   id_restart = register_restart_field (Ice_restart, restart_file, 'flux_v_bot',  Ice%flux_v_bot,  domain=Ice%domain)

 end subroutine ice_register_restart

!=============================================================================================

! dummy routine
 subroutine ice_stock_pe(Ice, index, value)
 type(ice_data_type), intent(in) :: Ice
 integer, intent(in) :: index
 real, intent(out)   :: value

 value = 0.0
 if(.not.Ice%pe) return

 if(.not.stock_warning_issued) then
   call error_mesg('ice_stock_pe','Stocks not yet implemented. Returning zero.',NOTE)
   stock_warning_issued = .true.
 endif

 end subroutine ice_stock_pe
!=============================================================================================

subroutine ice_data_type_chksum(id, timestep, data_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_data_type), intent(in) :: data_type
    integer ::   n, m, outunit
    
    outunit = stdout()

100 FORMAT("CHECKSUM::",A32," = ",Z20)
    write(outunit,*) 'BEGIN CHECKSUM(ice_data_type):: ', id, timestep
    write(outunit,100) 'ice_data_type%part_size          ',mpp_chksum(data_type%part_size          )
    write(outunit,100) 'ice_data_type%t_surf             ',mpp_chksum(data_type%t_surf             )
    write(outunit,100) 'ice_data_type%albedo             ',mpp_chksum(data_type%albedo             )
    write(outunit,100) 'ice_data_type%albedo_vis_dir     ',mpp_chksum(data_type%albedo_vis_dir     )
    write(outunit,100) 'ice_data_type%albedo_nir_dir     ',mpp_chksum(data_type%albedo_nir_dir     )
    write(outunit,100) 'ice_data_type%albedo_vis_dif     ',mpp_chksum(data_type%albedo_vis_dif     )
    write(outunit,100) 'ice_data_type%albedo_nir_dif     ',mpp_chksum(data_type%albedo_nir_dif     )
    write(outunit,100) 'ice_data_type%rough_mom          ',mpp_chksum(data_type%rough_mom          )
    write(outunit,100) 'ice_data_type%rough_heat         ',mpp_chksum(data_type%rough_heat         )
    write(outunit,100) 'ice_data_type%rough_moist        ',mpp_chksum(data_type%rough_moist        )
    write(outunit,100) 'ice_data_type%frazil             ',mpp_chksum(data_type%frazil             )
    write(outunit,100) 'ice_data_type%u_surf             ',mpp_chksum(data_type%u_surf             )
    write(outunit,100) 'ice_data_type%v_surf             ',mpp_chksum(data_type%v_surf             )
    write(outunit,100) 'ice_data_type%flux_u_bot         ',mpp_chksum(data_type%flux_u_bot         )
    write(outunit,100) 'ice_data_type%flux_v_bot         ',mpp_chksum(data_type%flux_v_bot         )
    write(outunit,100) 'ice_data_type%flux_t_bot         ',mpp_chksum(data_type%flux_t_bot         )
    write(outunit,100) 'ice_data_type%flux_q_bot         ',mpp_chksum(data_type%flux_q_bot         )
    write(outunit,100) 'ice_data_type%flux_lh_bot        ',mpp_chksum(data_type%flux_lh_bot        )
    write(outunit,100) 'ice_data_type%flux_sw_bot        ',mpp_chksum(data_type%flux_sw_bot        )
    write(outunit,100) 'ice_data_type%flux_sw_vis_bot    ',mpp_chksum(data_type%flux_sw_vis_bot    )
    write(outunit,100) 'ice_data_type%flux_sw_dir_bot    ',mpp_chksum(data_type%flux_sw_dir_bot    )
    write(outunit,100) 'ice_data_type%flux_sw_dif_bot    ',mpp_chksum(data_type%flux_sw_dif_bot    )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dir_bot',mpp_chksum(data_type%flux_sw_vis_dir_bot)
    write(outunit,100) 'ice_data_type%flux_sw_vis_dif_bot',mpp_chksum(data_type%flux_sw_vis_dif_bot)
    write(outunit,100) 'ice_data_type%flux_sw_nir_dir_bot',mpp_chksum(data_type%flux_sw_nir_dir_bot)
    write(outunit,100) 'ice_data_type%flux_sw_nir_dif_bot',mpp_chksum(data_type%flux_sw_nir_dif_bot)
    write(outunit,100) 'ice_data_type%flux_lw_bot        ',mpp_chksum(data_type%flux_lw_bot        )
    write(outunit,100) 'ice_data_type%lprec_bot          ',mpp_chksum(data_type%lprec_bot          )
    write(outunit,100) 'ice_data_type%fprec_bot          ',mpp_chksum(data_type%fprec_bot          )
    write(outunit,100) 'ice_data_type%runoff_bot         ',mpp_chksum(data_type%runoff_bot         )

    write(outunit,100) 'ice_data_type%flux_u             ',mpp_chksum(data_type%flux_u             )
    write(outunit,100) 'ice_data_type%flux_v             ',mpp_chksum(data_type%flux_v             )
    write(outunit,100) 'ice_data_type%flux_t             ',mpp_chksum(data_type%flux_t             )
    write(outunit,100) 'ice_data_type%flux_q             ',mpp_chksum(data_type%flux_q             )
    write(outunit,100) 'ice_data_type%flux_lh            ',mpp_chksum(data_type%flux_lh            )
    write(outunit,100) 'ice_data_type%flux_sw            ',mpp_chksum(data_type%flux_sw            )
    write(outunit,100) 'ice_data_type%flux_sw_vis        ',mpp_chksum(data_type%flux_sw_vis        )
    write(outunit,100) 'ice_data_type%flux_sw_dir        ',mpp_chksum(data_type%flux_sw_dir        )
    write(outunit,100) 'ice_data_type%flux_sw_dif        ',mpp_chksum(data_type%flux_sw_dif        )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dir    ',mpp_chksum(data_type%flux_sw_vis_dir    )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dif    ',mpp_chksum(data_type%flux_sw_vis_dif    )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dir    ',mpp_chksum(data_type%flux_sw_nir_dir    )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dif    ',mpp_chksum(data_type%flux_sw_nir_dif    )
    write(outunit,100) 'ice_data_type%flux_lw            ',mpp_chksum(data_type%flux_lw            )
    write(outunit,100) 'ice_data_type%lprec              ',mpp_chksum(data_type%lprec              )
    write(outunit,100) 'ice_data_type%fprec              ',mpp_chksum(data_type%fprec              )
    write(outunit,100) 'ice_data_type%p_surf             ',mpp_chksum(data_type%p_surf             )
    write(outunit,100) 'ice_data_type%runoff             ',mpp_chksum(data_type%runoff             )
    write(outunit,100) 'ice_data_type%calving            ',mpp_chksum(data_type%calving            )
    write(outunit,100) 'ice_data_type%flux_salt          ',mpp_chksum(data_type%flux_salt          )

    do n = 1, data_type%ocean_fields%num_bcs  !{
       do m = 1, data_type%ocean_fields%bc(n)%num_fields  !{
          !write(outunit,101) 'ice%', m, n, mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
          write(outunit,101) 'ice%',trim(data_type%ocean_fields%bc(n)%name), &
               trim(data_type%ocean_fields%bc(n)%field(m)%name), &
               mpp_chksum(data_type%ocean_fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum


subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_ice_boundary_type), intent(in) :: bnd_type
    integer ::   n, m, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(ocean_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'ocn_ice_bnd_type%u        ',mpp_chksum(bnd_type%u        )
    write(outunit,100) 'ocn_ice_bnd_type%v        ',mpp_chksum(bnd_type%v        )
    write(outunit,100) 'ocn_ice_bnd_type%t        ',mpp_chksum(bnd_type%t        )
    write(outunit,100) 'ocn_ice_bnd_type%s        ',mpp_chksum(bnd_type%s        )
    write(outunit,100) 'ocn_ice_bnd_type%frazil   ',mpp_chksum(bnd_type%frazil   )
    write(outunit,100) 'ocn_ice_bnd_type%sea_level',mpp_chksum(bnd_type%sea_level)
!    write(outunit,100) 'ocn_ice_bnd_type%data     ',mpp_chksum(bnd_type%data     )
100 FORMAT("CHECKSUM::",A32," = ",Z20)

    do n = 1, bnd_type%fields%num_bcs  !{
       do m = 1, bnd_type%fields%bc(n)%num_fields  !{
          write(outunit,101) 'oibt%',trim(bnd_type%fields%bc(n)%name), &
               trim(bnd_type%fields%bc(n)%field(m)%name), &
               mpp_chksum(bnd_type%fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ocn_ice_bnd_type_chksum

subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_ice_boundary_type), intent(in) :: bnd_type
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'atm_ice_bnd_type%u_flux          ',mpp_chksum(bnd_type%u_flux)          
    write(outunit,100) 'atm_ice_bnd_type%v_flux          ',mpp_chksum(bnd_type%v_flux)
    write(outunit,100) 'atm_ice_bnd_type%u_star          ',mpp_chksum(bnd_type%u_star)
    write(outunit,100) 'atm_ice_bnd_type%t_flux          ',mpp_chksum(bnd_type%t_flux)
    write(outunit,100) 'atm_ice_bnd_type%q_flux          ',mpp_chksum(bnd_type%q_flux)
    write(outunit,100) 'atm_ice_bnd_type%lw_flux         ',mpp_chksum(bnd_type%lw_flux)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dir ',mpp_chksum(bnd_type%sw_flux_vis_dir)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dif ',mpp_chksum(bnd_type%sw_flux_vis_dif)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dir ',mpp_chksum(bnd_type%sw_flux_nir_dir)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dif ',mpp_chksum(bnd_type%sw_flux_nir_dif)
    write(outunit,100) 'atm_ice_bnd_type%lprec           ',mpp_chksum(bnd_type%lprec)
    write(outunit,100) 'atm_ice_bnd_type%fprec           ',mpp_chksum(bnd_type%fprec)
    write(outunit,100) 'atm_ice_bnd_type%dhdt            ',mpp_chksum(bnd_type%dhdt)
    write(outunit,100) 'atm_ice_bnd_type%dedt            ',mpp_chksum(bnd_type%dedt)
    write(outunit,100) 'atm_ice_bnd_type%drdt            ',mpp_chksum(bnd_type%drdt)
    write(outunit,100) 'atm_ice_bnd_type%coszen          ',mpp_chksum(bnd_type%coszen)
    write(outunit,100) 'atm_ice_bnd_type%p               ',mpp_chksum(bnd_type%p)
!    write(outunit,100) 'atm_ice_bnd_type%data            ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine atm_ice_bnd_type_chksum

subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_boundary_type), intent(in) :: bnd_type
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(land_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'lnd_ice_bnd_type%runoff  ',mpp_chksum(bnd_type%runoff)
    write(outunit,100) 'lnd_ice_bnd_type%calving ',mpp_chksum(bnd_type%calving)
!    write(outunit,100) 'lnd_ice_bnd_type%data    ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine lnd_ice_bnd_type_chksum


end module ice_model_mod
