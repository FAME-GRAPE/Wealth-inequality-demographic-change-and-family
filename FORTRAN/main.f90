!##############################################################################
! Program is a one of the product of the 
! "Wealth inequality, demographic change and family" grand
!  The support of National Center for Science (grant #2016/22/E/HS4/00129) 
!   is greatefully acknowledged. 
!##############################################################################
    
    program olg2
    use global_vars
    use set_global_vars   
    use steady_state
    use get_data
    use clock 
 
    implicit none

    ! set paths for inputs and outputs
    call getcwd(cwd)
    cwd_r = trim(cwd)//"/Data"
    cwd_w = trim(cwd)//"/Results"    
    
call globals         ! globals is a subroutine in set_global_vars module                                                                

call tic()

    version = 'base_'
    write (*,*) 'We are performing base simulations path' 
    call globals         

    if (switch_run_1 == 1) then
        call steady(switch_unequal_bequest, switch_param_1, k_ss_1, r_ss_1, r_bar_ss_1, w_bar_ss_1, l_ss_j_1, s_ss_j_1, c_ss_j_1, b_ss_j_1, t1_ss_1)       
    endif ! run_1
    write(*,*) V_ss_j_vfi(1)
        
    switch_run_1 = 0 ! to be sure that we are running 2nd ss, in pfi procedure we are fullfiling 2nd part of transition
    if (switch_run_2 == 1) then
       call steady(switch_unequal_bequest, switch_param_2, k_ss_2, r_ss_2, r_bar_ss_2, w_bar_ss_2, l_ss_j_2, s_ss_j_2, c_ss_j_2, b_ss_j_2, t1_ss_2)
    endif ! run_2
    
    write (*,*) 'computations completed' 

call toc()

read*
endprogram olg2