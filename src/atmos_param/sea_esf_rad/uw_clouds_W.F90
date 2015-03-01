!FDOC_TAG_GFDL

                 module uw_clouds_W_mod
! <CONTACT EMAIL="fei.liu@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!          uw shallow convection cloud radiative properties module
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use  time_manager_mod, only: time_type
use           mpp_mod, only: input_nml_file
use           fms_mod, only: open_namelist_file, file_exist,   &
                             check_nml_error, error_mesg,   &
                             close_file, FATAL, NOTE, &
                             WARNING, mpp_pe, mpp_root_pe, &
                             write_version_number, stdlog
use     constants_mod, only: DENS_H2O, RDGAS, TFREEZE
use rad_utilities_mod, only: longwave_control_type, Lw_control, &
                             shortwave_control_type, Sw_control,&
                             microphysics_type,  &
                             microrad_properties_type, &
                             cld_specification_type, &
                             cloudrad_control_type, Cldrad_control

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          uw shallow convection cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: uw_clouds_W.F90,v 17.0.6.2 2010/09/07 16:17:19 wfc Exp $'
   character(len=128)  :: tagname =  '$Name: hiram_20101115_bw $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          uw_clouds_W_init,   &
          uw_clouds_W_end , uw_clouds_amt

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: dummy = .true.


namelist /uw_clouds_W_nml /     &
                                     dummy                          


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

real, parameter :: K_LAND  = 1.143
real, parameter :: K_OCEAN = 1.077
real, parameter :: N_LAND  = 250.E+06
real, parameter :: N_OCEAN = 100.E+06
real, parameter :: N_MIN   = 1.0e06
real, parameter :: QMIN    = 1.0e-10
real, parameter :: QAMIN   = 1.0e-2

  logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





! <SUBROUTINE NAME="uw_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_W_init  (pref, lonb, latb, axes, Time)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
!  <IN NAME="axes" TYPE="integer">
! 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine uw_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:),  intent(in) :: pref
real, dimension(:,:),  intent(in) :: lonb, latb
integer, dimension(4), intent(in) :: axes
type(time_type),       intent(in) :: Time

      integer            :: unit, ierr, io, logunit

     if (module_is_initialized) return
!---------------------------------------------------------------------
!-----  read namelist  ------
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=uw_clouds_W_nml, iostat=io)
      ierr = check_nml_error(io,"uw_clouds_W_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=uw_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'uw_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif
#endif

      if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         logunit = stdlog()
         write (logunit,nml=uw_clouds_W_nml)
      endif

!---------------------------------------------------------------------

       module_is_initialized = .true.


end subroutine uw_clouds_W_init

! <SUBROUTINE NAME="uw_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine uw_clouds_W_end
       
!----------------------------------------------------------------------
!    uw_clouds_W_end is the destructor for uw_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------


end subroutine uw_clouds_W_end


!#################################################################


!---------------------------------------------------------------------

! <SUBROUTINE NAME="uw_clouds_amt">
!  <OVERVIEW>
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_amt (is, ie, js, je, Shallow_microphys)   
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine uw_clouds_amt (is, ie, js, je,   &
                   shallow_cloud_area, shallow_liquid, shallow_ice, &
                   shallow_droplet_number, land, pfull, tkel, &
                   Shallow_microphys)

!---------------------------------------------------------------------
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!----------------------------------------------------------------------

integer,                 intent(in)    :: is,ie,js,je
real, dimension(:,:,:),  intent(in)    :: shallow_cloud_area,  &
                                          shallow_liquid, shallow_ice, &
                                          shallow_droplet_number
real, dimension(:,:),    intent(in)    :: land
real, dimension(:,:,:),  intent(in)    :: pfull, tkel
type(microphysics_type), intent(inout) :: Shallow_microphys

!--------------------------------------------------------------------------
!   local variables:

      real, dimension (size(shallow_liquid,1),                             &
                       size(shallow_liquid,2),                             &
                       size(shallow_liquid,3)) ::  liq_local, area_local,  &
                                                   ice_local, rho
      real, dimension (size(shallow_liquid,1),                             &
                       size(shallow_liquid,2)) ::  k_ratio, ndrops
      integer     :: ix, jx, kx
      integer     :: i, j, k

!----------------------------------------------------------------------
!    define vertical index.
!----------------------------------------------------------------------
      ix = size(Shallow_microphys%size_drop,1)
      jx = size(Shallow_microphys%size_drop,2)
      kx = size(Shallow_microphys%size_drop,3)

!---------------------------------------------------------------------
!    define k_ratio as appropriate mix of land and ocean values.
!----------------------------------------------------------------------
      k_ratio(:,:) = K_LAND*land(:,:) + K_OCEAN*(1.0 - land(:,:))
      ndrops(:,:)  = N_LAND*land(:,:) + N_OCEAN*(1.0 - land(:,:))
      
!---------------------------------------------------------------------
!    limit shallow cloud area and cloud amount so that non-zero clouds 
!    exist only in boxes with non-zero cloud areas.
!---------------------------------------------------------------------
      do k=1,kx                                      
        do j=1,jx
          do i=1,ix
            area_local(i,j,k) = 0.
            if ( (shallow_cloud_area(i,j,k) > QAMIN) .and.  &
                 (shallow_liquid(i,j,k) > QMIN) ) then
              liq_local(i,j,k) = shallow_liquid(i,j,k)
              area_local(i,j,k) = shallow_cloud_area(i,j,k)
            else
              liq_local(i,j,k) = 0.0                     
            endif

            if ( (shallow_cloud_area(i,j,k) > QAMIN) .and.  &
                 (shallow_ice(i,j,k) > QMIN)) then
              ice_local(i,j,k) = shallow_ice(i,j,k)
              area_local(i,j,k) = shallow_cloud_area(i,j,k)
            else
              ice_local(i,j,k) = 0.0                     
            endif
          end do
        end do
      end do
          
!----------------------------------------------------------------------
!    define droplet diameter based on cloud amount and droplet number. 
!    use formula as in cloud_rad_mod.
!----------------------------------------------------------------------
      if (Cldrad_control%do_liq_num) then
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (liq_local(i,j,k) > QMIN) then
              Shallow_microphys%size_drop(i,j,k)  =  2.0* &
                       (k_ratio(i,j)*620350.49*(liq_local(i,j,k)/DENS_H2O/ &
                       MAX(shallow_droplet_number(i,j,k), N_MIN*   &
                               max(area_local(i,j,k), QAMIN)/ &
                       (pfull(i,j,k)/RDGAS/tkel(i,j,k))))**(1./3.) )
            else
              Shallow_microphys%size_drop(i,j,k)  =  0.0 
            endif
          end do
        end do
      end do
      Shallow_microphys%droplet_number = shallow_droplet_number
      else
!----------------------------------------------------------------------
!  case of non-prognostic droplet number
!----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (liq_local(i,j,k) > QMIN) then
              Shallow_microphys%size_drop(i,j,k)  =    &
                 2.0* k_ratio(i,j)*620350.49*(pfull(i,j,k)*  &
                 shallow_liquid(i,j,k)/area_local(i,j,k)/RDGAS/  &
                 tkel(i,j,k)/DENS_H2O/ndrops(i,j))**(1./3.)
            else
              Shallow_microphys%size_drop(i,j,k)  =  0.0 
            endif
              Shallow_microphys%droplet_number(i,j,k) =   &
                                        ndrops(i,j)/(pfull(i,j,k)/  &
                                                    (RDGAS*tkel(i,j,k)))
          end do
        end do
      end do
      endif

!--------------------------------------------------------------------
!    if ice crystals are present, define their effective size, which
!    is a function of temperature. for ice clouds the effective radius
!    is taken from the formulation in Donner (1997, J. Geophys. Res., 
!    102, pp. 21745-21768) which is based on Heymsfield and Platt (1984)
!    with enhancement for particles smaller than 20 microns.  
!
!              T Range (K)               Reff (microns) 
!     -------------------------------    --------------
!
!     tfreeze-25. < T                      100.6           
!     tfreeze-30. < T <= Tfreeze-25.        80.8         
!     tfreeze-35. < T <= Tfreeze-30.        93.5             
!     tfreeze-40. < T <= Tfreeze-35.        63.9            
!     tfreeze-45. < T <= Tfreeze-40.        42.5           
!     tfreeze-50. < T <= Tfreeze-45.        39.9         
!     Tfreeze-55  < T <= Tfreeze-50         21.6          
!                   T <= Tfreeze-55.        20.2        
!
!--------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
      if (Cldrad_control%using_fu2007) then
!+yim Fu's parameterization of dge
        Shallow_microphys%size_ice(i,j,k) =   &
                            47.05 + 0.6624*(tkel(i,j,k) - TFREEZE) +  &
                            0.001741*(tkel(i,j,k)-TFREEZE)**2 
      else ! (using_fu2007)
            if (ice_local(i,j,k) > QMIN) then
              if (tkel(i,j,k) > TFREEZE - 25. ) then
                Shallow_microphys%size_ice(i,j,k) = 100.6    
              else if (tkel(i,j,k) >  TFREEZE - 30. .and. &
                       tkel(i,j,k) <= TFREEZE - 25.) then
                Shallow_microphys%size_ice(i,j,k) = 80.8       
              else if (tkel(i,j,k) >  TFREEZE - 35. .and. &
                       tkel(i,j,k) <= TFREEZE - 30.) then
                Shallow_microphys%size_ice(i,j,k)  = 93.5     
              else if (tkel(i,j,k) >  TFREEZE - 40. .and. &
                       tkel(i,j,k) <= TFREEZE - 35.) then
                Shallow_microphys%size_ice(i,j,k) = 63.9      
              else if (tkel(i,j,k) >  TFREEZE - 45. .and. &
                       tkel(i,j,k) <= TFREEZE - 40.) then
                Shallow_microphys%size_ice(i,j,k) = 42.5    
              else if (tkel(i,j,k) >  TFREEZE - 50. .and. &
                       tkel(i,j,k) <= TFREEZE - 45.) then
                Shallow_microphys%size_ice(i,j,k) = 39.9       
              else if (tkel(i,j,k) >  TFREEZE - 55. .and. &
                       tkel(i,j,k) <= TFREEZE - 50.) then
                Shallow_microphys%size_ice(i,j,k) = 21.6        
              else
                Shallow_microphys%size_ice(i,j,k) = 20.2        
              end if
            else
              Shallow_microphys%size_ice(i,j,k) = 0.0
            endif
      endif  ! (using_fu2007)
          end do
        end do
      end do

!---------------------------------------------------------------------
!    convert the cloud and ice amounts from kg(h2o) / kg(air) to 
!    g(h2o) / m**3, as required for use in the microphys_rad routines
!    which compute cloud radiative properties.
!---------------------------------------------------------------------
      rho(:,:,:) = pfull(:,:,1:kx)/(RDGAS*tkel(:,:,1:kx))
      Shallow_microphys%conc_drop = 1.0e03*rho*liq_local  
      Shallow_microphys%conc_ice  = 1.0e03*rho* ice_local  

!---------------------------------------------------------------------
!    define the cloud area to be that after adjustment for trivial cloud
!    amounts.
!---------------------------------------------------------------------
      Shallow_microphys%cldamt = area_local


!---------------------------------------------------------------------



end subroutine uw_clouds_amt  



!####################################################################


                     end module uw_clouds_W_mod

