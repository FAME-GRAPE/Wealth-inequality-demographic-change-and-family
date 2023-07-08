!***************************************************************************************

! get distribution for every gridpoint and state for every age
! Steady state
    subroutine get_distribution_ss()

    implicit none

    integer :: j, ial, iar, ind
    real*8 :: dist, ia_initial(n_beq), p_initial(n_beq), const

    ! set distribution to zero
    prob_ss = 0d0
    
    if((switch_unequal_bequest==1))then
        
        do ind=n_beq,1,-1
            p_initial(ind) = 1d0/ind**(zipf)/const_zipf !zipf law p.d.f. 
            ia_initial(ind) = 1d0/(2d0**(n_beq-ind))* bequest_ss_vfi/(p_initial(ind)*N_ss_j_vfi(1)) ! 1d0/(2d0**(n_beq-ind+1)) - number of people in one sub-cohort,  bequest_ss_vfi - sum of bequest, 
            !(p_initial(ind)*edu_rate*N_ss_j_vfi(1)) -  part of bequest which agent from ind- subcohort inherit  
        enddo
            ia_initial(2) = 1d0/(2d0**(n_beq-2))* bequest_ss_vfi/(p_initial(2)*N_ss_j_vfi(1))
            ia_initial(1) = 0d0
            
            const = sum(ia_initial)
            const = sum(ia_initial * p_initial)
            const = sum(p_initial)
            
        do ind=1,n_beq,1
            call linear_int(ia_initial(ind), ial, iar, dist, sv, n_a, a_grow)
            ial = min(ial, n_a)
            iar = min(iar, n_a)
            dist = min(dist, 1d0)

            prob_ss(1, ial, :, :, :)  = prob_ss(1, ial, :, :, :)  + p_initial(ind)*dist ! y-ss rename f_dens_ss ! poczatkowy rozk³ad
            prob_ss(1, iar, :,  :, :) = prob_ss(1, iar, :, :, :) +  p_initial(ind)*(1d0 - dist)
        enddo
            ! need to amend it to get initial dispersion
            const = 0.0d0
            do ia = 0, n_a, 1
                do ip = 1 , n_sp, 1
                    do ir=1, n_sr, 1
                        do id=1,n_sd,1 
                        prob_ss(1, ia, ip, ir, id) = prob_ss(1, ia, ip, ir, id) * pi_ip_init(ip) * pi_id_init(id) * pi_ir_init(ir)
                        const = const + prob_ss(1, ia, ip, ir, id)
                        enddo
                    enddo
                enddo
            enddo
                
    else  
            ! get initial distribution in age 1
            call linear_int(0d0, ial, iar, dist, sv, n_a, a_grow)
            ial = min(ial, n_a)
            iar = min(iar, n_a)
            dist = min(dist, 1d0)
            
            ! need to amend it to get initial dispersion
            prob_ss(1, ial, :, :, :) = dist ! y-ss rename f_dens_ss ! poczatkowy rozk³ad
            prob_ss(1, iar, :,  :, :) = 1d0 - dist
            
            do ia = 0, n_a, 1
                do ip = 1 , n_sp, 1
                    do ir=1, n_sr, 1
                        do id=1,n_sd,1 
                            prob_ss(1, ia, ip, ir, id) = prob_ss(1, ia, ip, ir, id) * pi_ip_init(ip) * pi_id_init(id) * pi_ir_init(ir)
                        enddo
                    enddo
                enddo
            enddo
                
    endif
            

              
    ! successively compute distribution over ages
    do j = 2, bigJ  
    ! iterate over yesterdays gridpoints
        do ia = 0, n_a, 1
            do ip = 1 , n_sp, 1
                do ir=1, n_sr, 1
                    do id=1,n_sd,1 
                        ! interpolate yesterday's savings decision
                        call linear_int(svplus_ss(j-1, ia, ip, ir, id), ial, iar, dist, sv, n_a, a_grow)
                        ! restrict values to grid just in case               
                        dist = min(abs(dist), 1d0)
                        ! redistribute households
                        do ip_p = 1, n_sp,1
                            do ir_r=1, n_sr, 1
                                do id_d =1, n_sd, 1
                                    prob_ss(j, ial, ip_p, ir_r, id_d) = prob_ss(j, ial, ip_p, ir_r, id_d) &
                                                                                + pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id, id_d)*dist      *prob_ss(j-1, ia, ip, ir, id)
                                    prob_ss(j, iar, ip_p, ir_r, id_d) = prob_ss(j, iar, ip_p, ir_r, id_d) &
                                                                                + pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id, id_d)*(1d0-dist)*prob_ss(j-1, ia, ip, ir, id) 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo 
    
 end subroutine