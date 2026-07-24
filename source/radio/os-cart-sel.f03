!-----------------------------------------------------------------------------------------
subroutine calc_dflt(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! Do nothing

end subroutine calc_dflt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine comp_select_cart_det( this )

  implicit none

  class(cart_Detector), intent(inout) :: this

  integer :: sele

  sele = 0
  if (this%cmpts(6)) sele = sele + 32
  if (this%cmpts(5)) sele = sele + 16
  if (this%cmpts(4)) sele = sele + 8
  if (this%cmpts(3)) sele = sele + 4
  if (this%cmpts(2)) sele = sele + 2
  if (this%cmpts(1)) sele = sele + 1

  select case(sele)
  case(1)
    this % fill_emf_loc => calc_comp_E1
  case(2)
    this % fill_emf_loc => calc_comp_E2
  case(3)
    this % fill_emf_loc => calc_comp_E1E2
  case(4)
    this % fill_emf_loc => calc_comp_E3
  case(5)
    this % fill_emf_loc => calc_comp_E1E3
  case(6)
    this % fill_emf_loc => calc_comp_E2E3
  case(7)
    this % fill_emf_loc => calc_comp_E1E2E3
#ifdef __HAS_RAD_BFLD__
  case(8)
    this % fill_emf_loc => calc_comp_B1
  case(9)
    this % fill_emf_loc => calc_comp_E1B1
  case(10)
    this % fill_emf_loc => calc_comp_E2B1
  case(11)
    this % fill_emf_loc => calc_comp_E1E2B1
  case(12)
    this % fill_emf_loc => calc_comp_E3B1
  case(13)
    this % fill_emf_loc => calc_comp_E1E3B1
  case(14)
    this % fill_emf_loc => calc_comp_E2E3B1
  case(15)
    this % fill_emf_loc => calc_comp_E1E2E3B1
  case(16)
    this % fill_emf_loc => calc_comp_B2
  case(17)
    this % fill_emf_loc => calc_comp_E1B2
  case(18)
    this % fill_emf_loc => calc_comp_E2B2
  case(19)
    this % fill_emf_loc => calc_comp_E1E2B2
  case(20)
    this % fill_emf_loc => calc_comp_E3B2
  case(21)
    this % fill_emf_loc => calc_comp_E1E3B2
  case(22)
    this % fill_emf_loc => calc_comp_E2E3B2
  case(23)
    this % fill_emf_loc => calc_comp_E1E2E3B2
  case(24)
    this % fill_emf_loc => calc_comp_B1B2
  case(25)
    this % fill_emf_loc => calc_comp_E1B1B2
  case(26)
    this % fill_emf_loc => calc_comp_E2B1B2
  case(27)
    this % fill_emf_loc => calc_comp_E1E2B1B2
  case(28)
    this % fill_emf_loc => calc_comp_E3B1B2
  case(29)
    this % fill_emf_loc => calc_comp_E1E3B1B2
  case(30)
    this % fill_emf_loc => calc_comp_E2E3B1B2
  case(31)
    this % fill_emf_loc => calc_comp_E1E2E3B1B2
  case(32)
    this % fill_emf_loc => calc_comp_B3
  case(33)
    this % fill_emf_loc => calc_comp_E1B3
  case(34)
    this % fill_emf_loc => calc_comp_E2B3
  case(35)
    this % fill_emf_loc => calc_comp_E1E2B3
  case(36)
    this % fill_emf_loc => calc_comp_E3B3
  case(37)
    this % fill_emf_loc => calc_comp_E1E3B3
  case(38)
    this % fill_emf_loc => calc_comp_E2E3B3
  case(39)
    this % fill_emf_loc => calc_comp_E1E2E3B3
  case(40)
    this % fill_emf_loc => calc_comp_B1B3
  case(41)
    this % fill_emf_loc => calc_comp_E1B1B3
  case(42)
    this % fill_emf_loc => calc_comp_E2B1B3
  case(43)
    this % fill_emf_loc => calc_comp_E1E2B1B3
  case(44)
    this % fill_emf_loc => calc_comp_E3B1B3
  case(45)
    this % fill_emf_loc => calc_comp_E1E3B1B3
  case(46)
    this % fill_emf_loc => calc_comp_E2E3B1B3
  case(47)
    this % fill_emf_loc => calc_comp_E1E2E3B1B3
  case(48)
    this % fill_emf_loc => calc_comp_B2B3
  case(49)
    this % fill_emf_loc => calc_comp_E1B2B3
  case(50)
    this % fill_emf_loc => calc_comp_E2B2B3
  case(51)
    this % fill_emf_loc => calc_comp_E1E2B2B3
  case(52)
    this % fill_emf_loc => calc_comp_E3B2B3
  case(53)
    this % fill_emf_loc => calc_comp_E1E3B2B3
  case(54)
    this % fill_emf_loc => calc_comp_E2E3B2B3
  case(55)
    this % fill_emf_loc => calc_comp_E1E2E3B2B3
  case(56)
    this % fill_emf_loc => calc_comp_B1B2B3
  case(57)
    this % fill_emf_loc => calc_comp_E1B1B2B3
  case(58)
    this % fill_emf_loc => calc_comp_E2B1B2B3
  case(59)
    this % fill_emf_loc => calc_comp_E1E2B1B2B3
  case(60)
    this % fill_emf_loc => calc_comp_E3B1B2B3
  case(61)
    this % fill_emf_loc => calc_comp_E1E3B1B2B3
  case(62)
    this % fill_emf_loc => calc_comp_E2E3B1B2B3
  case(63)
    this % fill_emf_loc => calc_comp_E1E2E3B1B2B3
#endif
  case default
    this % fill_emf_loc => calc_dflt
  end select

end subroutine comp_select_cart_det
!-----------------------------------------------------------------------------------------