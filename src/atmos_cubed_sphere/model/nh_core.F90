module nh_core_mod

! Notes:
! Using k_top=2 to treat the top layer hydrostatically so that delz will
! be computed using hydrostatic balance (instead of the update by
! advection of height using extrapolated winds at the model top)
!
! To do list:
! include moisture effect in pt
!------------------------------

   use fms_mod, only: error_mesg, FATAL

   implicit none
   private

   public Riem_Solver, Riem_Solver_C, update_dz_c, update_dz_d
   real, parameter:: dz_max = -0.5               ! (meters)

CONTAINS 

  subroutine update_dz_c(is, ie, js, je, km, ng, area, zh, ut, vt, dz_in, dz_out, wk)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt, zh
  real, intent(in ):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in ):: dz_in (is:ie,js:je,km) 
  real, intent(out):: dz_out(is:ie,js:je,km) 
  real, intent(out):: wk(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array
! Local Work array:
  real, dimension(is:ie+1,js:je  ):: xfx, fx
  real, dimension(is:ie  ,js:je+1):: yfx, fy
  integer  i, j, k

  call error_mesg('update_dz_c','The null version of update_dz_c should not be called.',FATAL)

  end subroutine update_dz_c



  subroutine update_dz_d(hord, is, ie, js, je, km, ng, npx, npy, area, zh, crx, cry, xfx, yfx, delz, wk, delp, n_sponge)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: hord, n_sponge
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout) ::delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(  out) ::   wk(is:ie,js:je,km)  ! work array
!-----------------------------------------------------
! Local array:
  real, dimension(is:   ie+1, js-ng:je+ng):: crx_adv, xfx_adv
  real, dimension(is-ng:ie+ng,js:   je+1 ):: cry_adv, yfx_adv
  real, dimension(is:ie+1,js:je  ):: fx
  real, dimension(is:ie  ,js:je+1):: fy
  real  delx(is:ie+1,km), dely(is-ng:ie+ng,km)
  real :: ra_x(is:ie,js-ng:je+ng)
  real :: ra_y(is-ng:ie+ng,js:je)
!--------------------------------------------------------------------
  integer  i, j, k, iord, isd, ied, jsd, jed, lm

  call error_mesg('update_dz_d','The null version of update_dz_d should not be called.',FATAL)

  end subroutine update_dz_d


  subroutine Riem_Solver_C(dt, is, ie, js, je, km, ng, akap, cp, ptop, hs, w, delz, pt, delp, gz, pk, ip)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ip       ! ip==1 pk is full pressure
   real, intent(in):: dt,  akap, cp, ptop
   real, intent(in):: delz(is:ie,js:je,km)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::       hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS 
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pk
! Local:
  real, dimension(is:ie,km  ):: pm, dm, dz2
  real, dimension(is:ie,km+1):: pem, pk2
  real gama, rgrav, ptk
  integer i, j, k
  integer m_split_c

  call error_mesg('Riem_Solver_C','The null version of Riem_Solver_C should not be called.',FATAL)

  end subroutine Riem_Solver_C


  subroutine Riem_Solver(dt, is, ie, js, je, km, ng, akap, cp, ptop, hs, peln, w, delz, pt, delp, gz, pkc, pk, pe, last_call, ip)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       pkc: full non-hydrostatic pressure
!--------------------------------------------
   integer, intent(in):: is, ie, js, je, km, ng
   integer, intent(in):: ip      ! ip==0 pkc is perturbation pressure
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: last_call
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pkc
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: peln(is:ie,km+1,js:je)           ! ln(pe)
! Local:
  real, dimension(is:ie,km):: pm, dm, dz2
  real :: pem(is:ie,km+1)
  real gama, rgrav, ptk
  integer i, j, k

  call error_mesg('Riem_Solver','The null version of Riem_Solver should not be called.',FATAL)

  end subroutine Riem_Solver


end module nh_core_mod

