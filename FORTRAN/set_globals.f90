! WHAT   : declare values for switches
! DO     : read initial values of switches and base (from base scenario run) values of basic variables
! RETURN : clean values for next run 

MODULE set_global_vars
USE global_vars
USE get_data
use pfi_trans

IMPLICIT NONE
CONTAINS

subroutine globals 
    real*8 :: pi_i_6, pi_6_6, pi_6_7, pi_7_7 ! super stars 
    integer :: ibeq
    
call chdir(cwd_r)
 
    experiment_no             = 0
    
!!! DEBUG_SWITCH
       switch_labor_choice      = 1         ! 0 = no labor choice (phi = 1) , 1 =  labor choice determined by 0<phi<1
       switch_fix_labor         = 0         ! if labor is fixed it is fixed to this number
       switch_unequal_bequest   = 0         ! 0 - bequests given by people of age j to people with age j-1, distributed equally; 1 - bequests given by all people to j=1, unequal distribution
       switch_persistent_delta  = 0
       switch_epsilon_corr      = 1
     
 ! loads some preset experiments based on expriment_no
      select case (experiment_no)  
            case(0)
            experiment = 'all_'
            switch_income_risk       = 1
            switch_discount_risk     = 1
            switch_return_risk       = 1
      end select


 ! switches related to transition experiments 
    switch_ss_write = 1        ! 0 - do not save big csv files with steady state, 1 save
    switch_run_1 = 1            ! 0 = don't run old steady state; 1 = run old steady state
    switch_run_2 = 1            ! 0 = don't run new steady state; 1 = run new steady state

! note: transition path is run only if the second steady state is run
    switch_param_1 = 0          ! 0 = with old parameters; 1 = with new parameters  
    switch_param_2 = 1          ! 0 = with old parameters; 1 = with new parameters  

!!!!!!!!!!!!!!!!!!!
    err_ss_tol = 1e-11
    
    up_ss = 0.7d0 
 
    superstar_factor_1 = 6.5d0
    superstar_factor_2 = 25.0d0

    depr = (1.0_dp + 0.050_dp)**zbar - 1.0_dp 
        
    theta = 2.0_dp
    delta =  (0.999_dp)**zbar 
    phi = 0.30d0 

    if (switch_labor_choice == 0) then
        phi  = 1.00_dp 
    endif

    ! grid definition 
    a_l    = 0.0d0   !dla bigJ = 80, a_l = -2d0, inaczej -8d0
    a_u    = 50d0   !dla bigJ = 80, a_u = 10d0, inaczej 30d0
    a_grow = 0.04d0 !dla bigJ = 80, a_grow = 0.05d0, inaczej 0.04d0        

    call read_data(omega_ss, pi, pi_weight)
    
    include 'shocks_parameters.f90'
    
    pi_ss_old = pi(:,1)
    pi_weight_ss_old = pi_weight(:,1)
  
    pi_ss_new = pi(:,2)
    pi_weight_ss_new = pi_weight(:,2)
  
    nu_ss_old =1.01d0
    nu_ss_new =1d0
    
    jbar_ss_old = 10 
    jbar_ss_new = 10 
    
    t1_ss_old =  0.0d0
    t1_ss_new =  0.0d0
    
    gam_ss_old = 1.12d0
    gam_ss_new = 1.03d0
    
    alpha_ss_old = 0.36d0
    alpha_ss_new = 0.43d0
    
    const_zipf = 0.0d0
    
    do ibeq=1,n_beq,1
            const_zipf = const_zipf + 1 / ibeq**(zipf)
    enddo
    
call chdir(cwd_w)    
end subroutine globals

end module set_global_vars