!Discount rate
    zeta_d = 0.95d0
    sigma_nu_d = 0.001d0 

    do id = 1 , n_sd, 1
        pi_id_init(id) = 1.0d0 / n_sd
    enddo 
    
    pi_id = 1d0
    n_sd_value = 0d0

     if (n_sd> 1) then 
        call discretize_AR(zeta_d, 0d0, sigma_nu_d, n_sd_value, pi_id)
       ! approximately stationary dist
       do t = 1,10, 1
            pi_id = matmul(pi_id,pi_id)
       enddo
       pi_id_init = pi_id(1,:)/sum(pi_id(1,:))
       call discretize_AR(zeta_d, 0d0, sigma_nu_d, n_sd_value, pi_id)
     endif
 
    if (switch_persistent_delta == 1) then
        call normal_discrete_1(n_sd_value, prob_norm_d, 1d0, sigma_nu_d)
         pi_id(:,:) = 0.0d0
        do  s = 1, n_sd, 1
            pi_id_init(s) = prob_norm_d(s)
            pi_id(s,s) = 1.0d0     
        enddo
    endif
    
    
! ----------------------------------------------------------------------------------------   
!Labor productivity
    zeta_p = 0.95d0**zbar   ! autoregression
    
    pi_ip = 0d0
    n_sp_value = 0d0

    ! definie initial distributions
        n_sp_initial = int(n_sp/2)+1
        n_sr_initial = int(n_sr/2)+1
        n_sd_initial = int(n_sd/2)+1
    
        if (n_sp>5) then
            n_sp_initial = int((n_sp-2)/2)+1  ! do not allow for people to be born as superstars
    endif
    
    sigma2_epsilon_ss_old = 0.013304d0
    sigma2_epsilon_ss_new = 0.025936d0

    ! implement correction of epsilons
        if (switch_epsilon_corr == 1) then
            epsilon_correction_ss_old = - (sigma2_epsilon_ss_old/(1.0d0 - zeta_p ** 2d0)) / 2.0d0
            epsilon_correction_ss_new = - (sigma2_epsilon_ss_new/(1.0d0 - zeta_p ** 2d0)) / 2.0d0
        else
            epsilon_correction_ss_old = 0.0d0
            epsilon_correction_ss_new = 0.0d0
    endif
      ! schock definition
        if (n_sp>5) then 
            ! get steady state shock realizations and transition matrices
            call discretize_AR(zeta_p, epsilon_correction_ss_old, sigma2_epsilon_ss_old, n_sp_value_ss_old(1:n_sp-2), pi_ip_ss_old(1:n_sp-2,1:n_sp-2))
            call discretize_AR(zeta_p, epsilon_correction_ss_new, sigma2_epsilon_ss_new, n_sp_value_ss_new(1:n_sp-2), pi_ip_ss_new(1:n_sp-2,1:n_sp-2))
        
            n_sp_value_ss_old = exp(n_sp_value_ss_old) 
            n_sp_value_ss_new = exp(n_sp_value_ss_new)     

            n_sp_value = exp(n_sp_value)  

            pi_i_6 = 5e-3
            pi_6_6 = 0.975d0
            pi_6_7 = 0.008d0
            pi_7_7 = 0.4d0

            n_sp_value_ss_old(n_sp-1) = superstar_factor_1*n_sp_value_ss_old(n_sp-2)
            n_sp_value_ss_old(n_sp) = superstar_factor_2*n_sp_value_ss_old(n_sp-1)
            n_sp_value_ss_new(n_sp-1) = superstar_factor_1*n_sp_value_ss_new(n_sp-2)
            n_sp_value_ss_new(n_sp) = superstar_factor_2*n_sp_value_ss_new(n_sp-1)
        
            pi_ip_ss_old = (1d0-pi_i_6)*pi_ip_ss_old
            pi_ip_ss_new = (1d0-pi_i_6)*pi_ip_ss_new
        
        do s=1, n_sp-2,1
            pi_ip_ss_old(s,n_sp-1) = pi_i_6
            pi_ip_ss_new(s,n_sp-1) = pi_i_6
        enddo
        
        pi_ip_ss_old(n_sp-1,n_sp-1) = pi_6_6
        pi_ip_ss_old(n_sp-1,n_sp)   = pi_6_7
        pi_ip_ss_old(n_sp-1,3)      = 1d0 - pi_6_7 -  pi_6_6 ! note it goes back to point = 3!
        pi_ip_ss_old(n_sp,n_sp)     = pi_7_7  
        pi_ip_ss_old(n_sp,n_sp-1)   = 1d0 - pi_7_7
    
        pi_ip_ss_new(n_sp-1,n_sp-1) = pi_6_6
        pi_ip_ss_new(n_sp-1,n_sp)   = pi_6_7
        pi_ip_ss_new(n_sp-1,3)      = 1d0 - pi_6_7 -  pi_6_6 ! note it goes back to point = 3!
        pi_ip_ss_new(n_sp,n_sp)      = pi_7_7  
        pi_ip_ss_new(n_sp,n_sp-1)    = 1d0 - pi_7_7
    

        do ip = 1 , n_sp, 1
            pi_ip_init_ss_old(ip) = pi_ip_ss_old(n_sp_initial,ip)
            pi_ip_init_ss_new(ip) = pi_ip_ss_new(n_sp_initial,ip)
        enddo

    
    elseif (n_sp>1)  then
        ! get steady state shock realizations and transition matrices
        call discretize_AR(zeta_p, epsilon_correction_ss_old, sigma2_epsilon_ss_old, n_sp_value_ss_old(1:n_sp), pi_ip_ss_old(1:n_sp,1:n_sp))
        call discretize_AR(zeta_p, epsilon_correction_ss_new, sigma2_epsilon_ss_new, n_sp_value_ss_new(1:n_sp), pi_ip_ss_new(1:n_sp,1:n_sp))
        n_sp_value_ss_old = exp(n_sp_value_ss_old) 
        n_sp_value_ss_new = exp(n_sp_value_ss_new)     
    else      
        pi_ip = 1d0
        n_sp_value = 1d0      
    endif
    ! now do initial things
    pi_ip_init_ss_old = pi_ip_ss_old(n_sp_initial,:)
    pi_ip_init_ss_new = pi_ip_ss_new(n_sp_initial,:)
! ----------------------------------------------------------------------------------------   
!Interest rate
    sigma_nu_r = 0.008d0 ** 2.0d0
        
    do ir = 1 , n_sr, 1
        pi_ir_init(ir) = 1.0d0 / n_sr
    enddo 
    pi_ir = 1d0
    n_sr_value = 0d0

    if (n_sr >1) then 
         call normal_discrete_1(n_sr_value, prob_norm, 0d0, sigma_nu_r)
         pi_ir(:,:) = 0.0d0
            do  s = 1, n_sr, 1
                pi_ir_init(s) = prob_norm(s)
                pi_ir(s,:) = prob_norm   
            enddo
    endif   

! CASE FOR DETERMINISTIC MODEL 

if (switch_income_risk == 0) then
    n_sp_value_ss_old(:)  = 1.0d0
    n_sp_value_ss_new(:)  = 1.0d0
endif
    
if (switch_discount_risk == 0) then
    n_sd_value(:) = 0.0d0
endif
    
if (switch_return_risk == 0) then
    n_sr_value(:) = 0.0d0
endif