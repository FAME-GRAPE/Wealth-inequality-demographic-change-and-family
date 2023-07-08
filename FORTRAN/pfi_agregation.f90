!***************************************************************************************
! find aggegate variables for steady state

    subroutine aggregation_ss()
    
        implicit none
        
        integer :: ij, ial, iar, ia_last, tt
        real*8 :: check_e(bigJ), check_euler_cum, euler_max(bigJ), w_sum(0:bigJ), bc,  ERHS_ss(bigJ, 0:n_a,  n_sp, n_sr,n_sd),  sum_y, dist 
        real*8 :: sum_help, p_1_5(bigJ)
    
       ! write(*,*) c_ss(1,0,3)
        ! calculate cohort aggregates
        savings_cohort_ten = 0d0
        c_ss_j_vfi(:) = 0d0
        V_ss_j_vfi(:) = 0d0
        l_ss_j_vfi(:) = 0d0
        s_pom_ss_j_vfi(:) = 0d0
        lab_ss_j_vfi(:) = 0d0
        asset_pom_ss_j(:) = 0d0
        l_ss_pen_j(:) = 0d0
        w_sum(0) = 0d0
        ERHS_ss = 0d0
        top_ten(:) = 0d0 
        top_ten_coh(:) =0d0
        top_100 =0d0
        savings_top_100 = 0d0
        lab_high = 0d0
        gini_weight_sv = 0d0
        gini_weight_consumption = 0d0

        share_neg = 0d0
        share_nonpos = 0d0

        lw_ss_j_vfi(:) = 0d0
        t = 10
        tt =100
        if (n_sp > 5) then
            do j=1,bigJ,1
                p_1_5(j) = sum(prob_ss(j,:,1:5,:,:))
            enddo
        else 
            do j=1,bigJ,1
                p_1_5(j) = 1d0
            enddo
        endif
        savings_top_ten =0d0
        cons_proc_top_ten_coh(:)= 0d0
        consumption_top_ten = 0d0
        call linear_int(0d0, ial, iar, dist, sv, n_a, a_grow)
        ial = min(ial, n_a)
        iar = min(iar, n_a)
        dist = min(dist, 1d0)
                                        
        ia_last = n_a
            do ia = n_a, 0 , -1
                    do j = 1, bigJ 
                        do ip = 1, n_sp, 1  
                            do ir=1, n_sr, 1
                                do id = 1, n_sd, 1
                                    gini_weight_sv(j,ia) = gini_weight_sv(j,ia) + prob_ss(j, ia, ip, ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                    gini_weight_consumption(j, ia, ip,ir, id) = gini_weight_consumption(j, ia, ip,ir, id) + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                    gini_income(j, ia, ip,ir, id) = omega_ss(j)*n_sp_value(ip)*l_ss(j, ia, ip, ir, id)*w_pom_ss_vfi(j)
                                    if (top_ten(t) >= 0.1d0) then 
                                        t = t-1
                                        top_ten(t) =  top_ten(t) + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                        savings_top_ten(t) = savings_top_ten(t) + sv(ia)*prob_ss(j, ia, ip, ir, id)*N_ss_j_vfi(j)
                                    else
                                        top_ten(t) =  top_ten(t) + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                        savings_top_ten(t) = savings_top_ten(t) + sv(ia)*prob_ss(j, ia, ip, ir, id)*N_ss_j_vfi(j)
                                    endif
                                    
                                        
                                    if (ia < ial) then
                                        share_neg = share_neg +  prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                    endif
                                
                                    if (ia <= iar) then
                                        share_nonpos =  share_nonpos + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                    endif
                                    
                                     if (top_100(tt) >= 0.01d0) then 
                                        tt = tt-1
                                        top_100(tt) =  top_100(tt) + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                        savings_top_100(tt) = savings_top_100(tt) + sv(ia)*prob_ss(j, ia, ip, ir, id)*N_ss_j_vfi(j)
                                    else
                                        top_100(tt) =  top_100(tt) + prob_ss(j, ia, ip,ir, id)*N_ss_j_vfi(j)/sum(N_ss_j_vfi)
                                        savings_top_100(tt) = savings_top_100(tt) + sv(ia)*prob_ss(j, ia, ip, ir, id)*N_ss_j_vfi(j)
                                    endif
                                    if (top_ten_coh(j) <= 0.1d0 ) then 
                                
                                        top_ten_coh(j) =  top_ten_coh(j) + prob_ss(j, ia, ip,ir, id)
                                        savings_cohort_ten(1,j) = savings_cohort_ten(1,j) + sv(ia)*prob_ss(j, ia, ip, ir, id)*10d0
                                        savings_cohort_ten(3,j) = savings_cohort_ten(3,j) + sv(ia)*prob_ss(j, ia, ip, ir, id)
                                        consumption_top_ten(1,j) =  consumption_top_ten(1,j) + c_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)*10d0
                                        consumption_top_ten(3,j) =  consumption_top_ten(3,j) + c_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)
                                    elseif(top_ten_coh(j) >= 0.9d0)then
                                        top_ten_coh(j) =  top_ten_coh(j) + prob_ss(j, ia, ip,ir, id)
                                        savings_cohort_ten(2,j) = savings_cohort_ten(2,j) + sv(ia)*prob_ss(j, ia, ip, ir, id)*10
                                        savings_cohort_ten(3,j) = savings_cohort_ten(3,j) + sv(ia)*prob_ss(j, ia, ip, ir, id)
                                        consumption_top_ten(3,j) =  consumption_top_ten(3,j) + c_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)
                                        consumption_top_ten(2,j) =  consumption_top_ten(2,j) + c_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)*10
                                    else
                                        top_ten_coh(j) =  top_ten_coh(j) + prob_ss(j, ia, ip,ir, id)
                                        savings_cohort_ten(3,j) = savings_cohort_ten(3,j) + sv(ia)*prob_ss(j, ia, ip, ir, id)
                                        consumption_top_ten(3,j) =  consumption_top_ten(3,j) + c_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)
                                    endif
                            
                                    c_ss_j_vfi(j) = c_ss_j_vfi(j) + c_ss(j, ia, ip, ir, id)*prob_ss(j, ia, ip,ir, id)
                            
                                    if(ip<6)then
                                        l_ss_pen_j(j) = l_ss_pen_j(j) + omega_ss(j)*n_sp_value(ip)*l_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)/p_1_5(j)
                                    endif 
                                    if(ip>=6)then
                                        lab_high(j) = lab_high(j) + l_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)/(1d0-p_1_5(j))
                                    endif
                                    l_ss_j_vfi(j) = l_ss_j_vfi(j) + omega_ss(j)*n_sp_value(ip)*l_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip,ir, id)
                                    lab_ss_j_vfi(j) = lab_ss_j_vfi(j) + l_ss(j, ia, ip, ir, id)*prob_ss(j, ia, ip, ir, id) 
                                    lw_ss_j_vfi(j) = lw_ss_j_vfi(j) + omega_ss(j)*n_sp_value(ip)*w_pom_ss_vfi(j)*l_ss(j, ia, ip,ir, id)*prob_ss(j, ia, ip, ir, id) 
                                    s_pom_ss_j_vfi(j) = s_pom_ss_j_vfi(j) + svplus_ss(j, ia, ip, ir, id)*prob_ss(j, ia, ip, ir, id)
                                    asset_pom_ss_j(j) = asset_pom_ss_j(j) + sv(ia)*prob_ss(j, ia, ip, ir, id) 
                                    V_ss_j_vfi(j)  =  V_ss_j_vfi(j) + V_ss(j, ia, ip, ir, id)*prob_ss(j, ia, ip, ir, id)
                                    sum_y = sum_y + prob_ss(j, ia, ip, ir, id)  
                                enddo
                            enddo
                        enddo
                        if(j>1)then ! todo18 
                            w_sum(j) = w_sum(j-1) + b_ss_j_vfi(j) + l_ss_j_vfi(j)*w_pom_ss_vfi(j) + bequest_ss_j_vfi(j) 
                        else
                            w_sum(j) = 0d0 + b_ss_j_vfi(j) + l_ss_j_vfi(j)*w_pom_ss_vfi(j) + bequest_ss_j_vfi(j) 
                        endif
                    ia_last = ia
                enddo
            enddo
            
            do j = 1,bigJ-1,1
                check_e(j)=sum(ERHS_ss(j,:,:, :,:)) ! check euler (only of theta = 1) with other theta marginal utility is not linear 
            enddo
            check_euler_cum = sum(check_e)
       
            bc = sum(c_ss_j_vfi(1:bigJ))- (r_ss_vfi/gam_ss_vfi-1d0)*sum(s_pom_ss_j_vfi(1:bigJ)) - w_sum(bigJ) ! budget constraint 
            euler_max = abs(check_e)

        open(unit = 104, file= "euler_ss.csv")
            do j = 2,bigJ,1
                write(104, '(F20.10)') check_e(j)
            enddo
        close(104)
              
    end subroutine

!!***************************************************************************************