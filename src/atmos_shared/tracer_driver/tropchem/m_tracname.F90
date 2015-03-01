
      module m_tracname_mod
!-----------------------------------------------------------
!       ... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use mo_grid_mod,   only : pcnst
      use chem_mods_mod, only : grpcnt

      implicit none

character(len=128), parameter :: version     = '$Id: m_tracname.F90,v 13.0.14.1 2010/03/25 00:31:42 pjp Exp $'
character(len=128), parameter :: tagname     = '$Name: hiram_20101115_bw $'
logical                       :: module_is_initialized = .false.

      save

      character(len=8) :: tracnam(pcnst)          ! species names
      character(len=8) :: natsnam(max(1,grpcnt))  ! names of non-advected trace species

      end module m_tracname_mod
