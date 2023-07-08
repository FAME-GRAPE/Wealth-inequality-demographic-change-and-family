! WHAT  : declaration of global (ALL subroutines and functions MAY USE them) parameters and variables
! TAKE  : none 
! DO    : definition of a switch in comments 
! RETURN: nothing
MODULE global_vars
IMPLICIT NONE
   save
    integer, parameter ::  n_iter_ss =  50
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: bigJ = 16
	integer, parameter :: bigT = 2+bigJ+1
    real(dp), parameter :: zbar = real(80/bigJ) ! 4 periods model scale parameter
    integer :: iter, i, j, s
    character(5) :: version
    
    integer :: experiment_no
    character(4) :: experiment
    character(4) :: no_steady
    
    character(128) :: cwd, cwd_r, cwd_w
       
    integer :: switch_ss_write                              ! 0 - do not save big csv files with steady state, 1 save
    
    integer :: switch_run_1, switch_run_2                   ! 0 = don't run the first/second steady state; 1 = run the first/second steady state/
    integer :: switch_param_1, switch_param_2               ! 0 = with old parameters; 1 = with new parameters (pi,gam,N,jbar)                
    
    integer :: switch_epsilon_corr                          ! 0 do not worry that variances shift means, 1 - correct to keep mean constant

    integer :: switch_unequal_bequest
    
    integer :: switch_labor_choice                         ! 0 = no labor choice (phi = 1) , 1 =  labor choice determined by 0<phi<1
    integer :: switch_persistent_delta                     ! 0 = AR1 shocks to patience, 1 = permanent types assigned at birth
    integer :: switch_income_risk
    integer :: switch_discount_risk
    integer :: switch_return_risk
    
    real*8  :: switch_fix_labor                             ! 0 = endogenous labor, other number (=0.33 for US) fix labor force participation
                                                           

    ! Deklaracja zmiennych wczytywanych
    real(dp), dimension(bigJ) :: omega_ss 


! Deklaracje zmiennych, ktore nam zostaja po steady state'ach
    real(dp) :: k_ss_1, r_ss_1, r_bar_ss_1, w_bar_ss_1, t1_ss_1
    real(dp) :: k_ss_2, r_ss_2, r_bar_ss_2, w_bar_ss_2, t1_ss_2
    real(dp), dimension(bigJ) :: l_ss_j_1, s_ss_j_1, c_ss_j_1, b_ss_j_1, l_ss_pen_j_1
    real(dp), dimension(bigJ) :: l_ss_j_2, s_ss_j_2, c_ss_j_2, b_ss_j_2, l_ss_pen_j_2
    
! Parametry
    real(dp) ::  delta, depr, theta, phi, up_ss, err_ss_tol
    real(dp) :: tl_ss, t1_ss_old, t1_ss_new, alpha_ss_old, alpha_ss_new
    real(dp) :: jbar_ss_old, jbar_ss_new, gam_ss_old, gam_ss_new, nu_ss_old, nu_ss_new, epsilon_correction_ss_old, epsilon_correction_ss_new
    real(dp), dimension(bigJ) :: pi_ss_old, pi_ss_new, pi_weight_ss_old, pi_weight_ss_new
    real(dp) :: superstar_factor_1, superstar_factor_2

! transition variables
    real(dp), dimension(bigJ) :: omega
    real(dp), dimension(bigJ,2) :: pi, pi_weight
  
 ! pfi 
    real*8, parameter  :: fi = (5d0**(1d0/2d0)-1d0)/2d0
    integer, parameter :: n_a = 110, n_sp = 5, n_sd = 5, n_sr = 5, n_beq = 2
    real*8, parameter  ::  zipf = 1.5d0  
    real*8 :: const_zipf
    real*8 :: zeta_p, a_l, a_u, a_grow,  poss_ass_sum_ss(bigJ),  sigma_nu_p, n_sp_initial, sigma_nu_r, n_sr_initial,&
                zeta_r, r_ss_, zeta_d, n_sd_initial, sigma_nu_d, &
               pi_ir(n_sr,n_sr), n_sr_value(n_sr), pi_id(n_sd,n_sd), n_sd_value(n_sd), prob_norm(n_sr),  prob_norm_d(n_sd), pi_id_init(n_sd), pi_ir_init(n_sr), pi_ip_init(n_sp)
    real(dp) :: sigma2_epsilon_ss_old, sigma2_epsilon_ss_new ! steady state variances
    real*8   :: pi_ip_ss_old(n_sp,n_sp), n_sp_value_ss_old(n_sp), pi_ip_ss_new(n_sp,n_sp), n_sp_value_ss_new(n_sp), pi_ip_init_ss_old(n_sp), pi_ip_init_ss_new(n_sp) ! steady state shock realizations and transition probabilities
    real*8   :: pi_ip(n_sp,n_sp), n_sp_value(n_sp) ! holder to make this code compatibile with older subroutines
    
end module global_vars