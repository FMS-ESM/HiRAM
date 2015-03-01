module ocean_tpm_mod  !{
! 
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Null version of Ocean tracer package module
!</OVERVIEW>
!
!<DESCRIPTION>
! Null version of Ocean tracer package module
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!

!
!       Place tracer modules here
!

!
!       force all variables to be "typed"
!

implicit none

!
!       Set all variables to be private by default

private

!
!       Private routines
!

!
!       Public routines
!

!public ocean_tpm_bbc
!public ocean_tpm_end
!public ocean_tpm_init
public ocean_tpm_flux_init
!public ocean_tpm_sbc
!public ocean_tpm_source
!public ocean_tpm_start
!public ocean_tpm_tracer
public ocean_tpm_init_sfc
!public ocean_tpm_sum_sfc
!public ocean_tpm_avg_sfc
!public ocean_tpm_zero_sfc
!public ocean_tpm_sfc_end

!
!       private parameters
!

character(len=48), parameter                    :: mod_name = 'ocean_tpm_mod'

!
!       Public variables
!
!
!       Private variables
!

character(len=128) :: version = '$Id: ocean_tpm.F90,v 13.0 2006/03/28 21:23:57 fms Exp $'
character(len=128) :: tagname = '$Name: hiram_20101115_bw $'

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_init_sfc">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
!
!       Note: this subroutine should be merged into ocean_tpm_start
! </DESCRIPTION>
!

subroutine ocean_tpm_init_sfc  !{

return

end subroutine ocean_tpm_init_sfc  !}
! </SUBROUTINE> NAME="ocean_tpm_init_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>
!

subroutine ocean_tpm_flux_init  !{

return

end subroutine ocean_tpm_flux_init  !}
! </SUBROUTINE> NAME="ocean_tpm_flux_init"

end module ocean_tpm_mod  !}
