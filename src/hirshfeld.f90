! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Hirshfeld integration
module hirshfeld
  implicit none

  private

  public :: hirsh_grid
  public :: hirsh_nogrid
  public :: hirsh_weights
  public :: voronoi_grid
  ! iterative Hirshfeld (Hirshfeld-I, Bultinck 2007)
  public :: hirsh_i_driver
  public :: hirsh_i_eval
  public :: hirsh_i_cleanup
  ! mesh / charge-aware reference helpers (used by molecular HI-XDM)
  public :: hirsh_i_prepare
  public :: hirsh_i_refrho
  public :: hirsh_i_qfloor
  public :: hirsh_i_cache_clean

  interface
     module subroutine hirsh_grid(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine hirsh_grid
     module subroutine hirsh_weights(s,bas,idb,w)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(in) :: bas
       integer, intent(in) :: idb
       real*8, intent(out) :: w(:,:,:)
     end subroutine hirsh_weights
     module subroutine voronoi_grid(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine voronoi_grid
     module subroutine hirsh_nogrid()
     end subroutine hirsh_nogrid
     module subroutine hirsh_i_driver(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine hirsh_i_driver
     module subroutine hirsh_i_eval(bas,idcel,dist,rho)
       use types, only: basindat
       type(basindat), intent(in) :: bas
       integer, intent(in) :: idcel
       real*8, intent(in) :: dist
       real*8, intent(out) :: rho
     end subroutine hirsh_i_eval
     module subroutine hirsh_i_cleanup(bas)
       use types, only: basindat
       type(basindat), intent(inout) :: bas
     end subroutine hirsh_i_cleanup
     module subroutine hirsh_i_prepare(s,qcel,wfcdir)
       use systemmod, only: system
       type(system), intent(in) :: s
       real*8, intent(in) :: qcel(:)
       character(len=*), intent(in) :: wfcdir
     end subroutine hirsh_i_prepare
     module function hirsh_i_refrho(iz,qreal,dist) result(rho)
       integer, intent(in) :: iz
       real*8, intent(in) :: qreal, dist
       real*8 :: rho
     end function hirsh_i_refrho
     module function hirsh_i_qfloor(iz,wfcdir) result(qf)
       integer, intent(in) :: iz
       character(len=*), intent(in) :: wfcdir
       real*8 :: qf
     end function hirsh_i_qfloor
     module subroutine hirsh_i_cache_clean()
     end subroutine hirsh_i_cache_clean
  end interface

end module hirshfeld
