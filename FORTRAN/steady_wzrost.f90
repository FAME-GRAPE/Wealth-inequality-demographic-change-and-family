! WHAT   : stady state routine for PAYG 
! TAKE   : data on mortality [[pi]], retirement age [[jbar]], change in technological progress [[gamma = z_t/z_(t-1)]], the size of the cohort [[N_ss_j]], pension system parameters [[t1 - ZUS, [[t2 - OFE]]  <- constant during iteration
!          and guess for capital, labor force participation [[l_]], consumption, bequest , nu_ss population growth rate in steady state, N_ss_j population age structure
! DO     : during iterations update guess values to new ones in order to find optimal allocation in steady state; iterations end when difference between old and new value is really small (err_ss < 1e-8)
! RETURN : k_ss, c_ss, r_ss etc. that fullfil the condition of optimal allocation. It is one of the edge of transition path

MODULE steady_state
use global_vars
!use individual_vf
use pfi_trans


IMPLICIT NONE

CONTAINS

subroutine steady(switch_unequal_bequest, param_ss, k_ss, r_ss, r_bar_ss,  w_bar_ss, l_ss_j, s_ss_j, c_ss_j, b_ss_j, t1_ss)
    real(dp) :: k_ss_new, err_ss, &
                jbar_ss, alpha, gam_ss, N_ss, nu_ss, bigl_ss, y_ss, consumption_ss_gross, &
                savings_ss, bequest_ss, income_ss, &
                Tax_ss, sum_b_ss, sum_priv_sv_ss, valor_mult_ss,  replacement_ss
    real(dp), dimension(bigj) :: pi_ss, life_exp, pi_weight_ss
	real(dp), dimension(bigj) :: N_ss_j, bequest_left_ss_j, bequest_ss_j, bequest_ss_j_old

    integer, intent(in)   :: param_ss
    integer, intent(in)   :: switch_unequal_bequest			
    real(dp), intent(out) :: k_ss, r_ss, r_bar_ss, w_bar_ss, t1_ss
    real(dp), dimension(bigj), intent(out) :: l_ss_j, s_ss_j, c_ss_j, b_ss_j
    
    ! pension system 
     real(dp) ::  avg_wl, tot_contrib

     
     real(dp), dimension(bigj, n_a) :: V_ss_j

    if (param_ss == 0) then ! 0 = with old parameters;  1 = with new parameters
        alpha = alpha_ss_old
        gam_ss = gam_ss_old
        pi_ss = pi_ss_old 
        pi_weight_ss = pi_weight_ss_old
        jbar_ss = jbar_ss_old
        nu_ss =  nu_ss_old
        t1_ss = t1_ss_old
        pi_ip = pi_ip_ss_old
        n_sp_value = n_sp_value_ss_old
        pi_ip_init = pi_ip_init_ss_old
   
    else 
        alpha = alpha_ss_new
        gam_ss = gam_ss_new
        pi_ss = pi_ss_new
        pi_weight_ss = pi_weight_ss_new
        jbar_ss = jbar_ss_new
        nu_ss = nu_ss_new
        t1_ss = t1_ss_new
        pi_ip = pi_ip_ss_new
        n_sp_value = n_sp_value_ss_new
        pi_ip_init = pi_ip_init_ss_new

    endif
    
!normalized structure of population such as N_ss_j(1) = 1 (number of 20 years old) 
    N_ss_j(1) = 1.0_dp 
    do j = 2, bigj
        N_ss_j(j) = nu_ss**(-j+1)*pi_weight_ss(j)/pi_weight_ss(1)
    enddo


    valor_mult_ss = (1 + (nu_ss*gam_ss - 1))/gam_ss 
    N_ss = sum(N_ss_j(1:bigJ))       

    ! guess 
    r_bar_ss = (1 + 0.05_dp)**(zbar) - 1 !(1 + 0.078_dp)**(zbar) - 1
    k_ss = ((r_bar_ss + depr)/(alpha*zbar))**(1/(alpha - 1))
    w_bar_ss = zbar*(1 - alpha)*k_ss**alpha
    LabIncAVG_ss_vfi = 0.33*w_bar_ss 

    bequest_ss = 0.0_dp
    bequest_ss_j = 0.0_dp
    bequest_left_ss_j = 0.0_dp
    bequest_ss_j_old = 0.0_dp
    
    
!!! ITERATIONS STARTS     
do iter = 1,n_iter_ss,1
                      
    r_bar_ss = zbar*alpha*k_ss**(alpha - 1) - depr
    w_bar_ss = zbar*(1 - alpha)*k_ss**alpha
    y_ss =zbar* k_ss**(alpha)
    
    if (r_bar_ss < 0) then
        r_bar_ss = 0
    endif
    
    r_ss = 1 + r_bar_ss  

    sum_priv_sv_ss = k_ss*gam_ss*nu_ss 

! no interest is added when switch_unequal_bequest == 1

    if (switch_unequal_bequest==0) then
        do j = 2,bigJ,1
            bequest_left_ss_j(j-1) = (pi_weight_ss(j-1) - pi_weight_ss(j))*(r_ss*s_ss_j(j-1))/gam_ss
        enddo
        bequest_left_ss_j(bigJ) = (pi_weight_ss(bigJ))*(r_ss*s_ss_j(bigJ))/gam_ss
        bequest_ss = sum(bequest_left_ss_j(1:bigJ))
        
        bequest_ss_j_old = bequest_ss_j
        bequest_ss_j(1) = 0d0
    
        do j = 2,bigJ,1
            bequest_ss_j(j) = up_ss*bequest_ss_j_old(j) + (1 - up_ss)*bequest_left_ss_j(j-1)/pi_weight_ss(j) 
        enddo  

    elseif (switch_unequal_bequest==1) then
        bequest_ss_j(1) = 0d0
        do j = 2,bigJ,1
            bequest_ss_j(j) = 0d0
            bequest_left_ss_j(j-1) = (pi_weight_ss(j-1) -   pi_weight_ss(j))*s_ss_j(j-1) / nu_ss**(j-1)
        enddo
        bequest_left_ss_j(bigJ) = pi_weight_ss(bigJ)*s_ss_j(bigj) * nu_ss**(-bigj+1)
        bequest_ss = sum(bequest_left_ss_j(1:bigJ))
    endif  
      

    ! household problem 
        w_pom_ss_vfi = (1.0_dp - t1_ss)*w_bar_ss 
        r_ss_vfi = r_bar_ss  
        bequest_ss_vfi =  bequest_ss
        gam_ss_vfi = gam_ss
        pi_ss_vfi = pi_ss
        b_ss_j_vfi = b_ss_j
        bequest_ss_j_vfi =  bequest_ss_j
        jbar_ss_vf = ceiling(jbar_ss)
        N_ss_j_vfi =  N_ss_j
        call agent_vf()
        c_ss_j = c_ss_j_vfi
        l_ss_j = l_ss_j_vfi
        s_ss_j(:) =  s_pom_ss_j_vfi(:)
        
    ! Futher aggregation        
    LabIncAVG_ss_vfi =  sum(N_ss_j(1:jbar_ss-1)*l_ss_j_vfi(1:jbar_ss-1)*w_pom_ss_vfi(1:jbar_ss-1))/sum(N_ss_j(1:jbar_ss-1))
    bigl_ss = sum(N_ss_j*l_ss_j(1:jbar_ss-1))        
    consumption_ss_gross  = sum(c_ss_j*N_ss_j(1:bigJ))/bigl_ss
    savings_ss = sum(N_ss_j*s_ss_j(1:bigJ))/bigl_ss !to samo co sum_priv_sv_ss w FF
    ! pension system 
    b_ss_j(1:jbar_ss-1) = 0
    tot_contrib = sum(N_ss_j*t1_ss*w_bar_ss*l_ss_pen_j(1:jbar_ss-1))
    b_ss_j(jbar_ss:bigJ)  = tot_contrib/sum(N_ss_j(jbar_ss:BigJ))
    sum_b_ss = sum(b_ss_j*N_ss_j(1:bigJ))/bigl_ss
 
         
    k_ss_new = savings_ss/(gam_ss*nu_ss)
    err_ss = abs(k_ss_new - k_ss)
    k_ss = up_ss*k_ss + (1 - up_ss)*k_ss_new

        if (mod(iter,1) == 0) then
            print*, iter, 'err_ss:', err_ss, 'feas_ss:', abs((y_ss - consumption_ss_gross)/y_ss - ((nu_ss*gam_ss+depr-1)*k_ss)/y_ss)
        endif
        if (err_ss < err_ss_tol ) then
            exit
        endif

enddo 

    replacement_ss = b_ss_j(jbar_ss)/((1 - t1_ss)*w_bar_ss*l_ss_pen_j(jbar_ss-1))   
 
    include 'Print_steady_db.f90'
  
end subroutine steady

END MODULE steady_state