!*******************************************************************************************
! find futur assets for every age, assets grid point, state  
! steady state
subroutine household_endo()

implicit none
real*8 :: available, EV_prim, c_opt, lab_income
real*8 :: av, wage, foc(2), optimal_choice(2)
real*8 :: dist, c_help, l_help

       
do ia = 0, n_a, 1 
    do ip=1, n_sp, 1 ! labor shock realisation 
        do ir =1, n_sr,1 ! intrest rate realisation
            do id = 1, n_sd,1 ! discount rate realisation 
                c_ss(bigj, ia, ip, ir, id) = max(((1d0+n_sr_value(ir)+r_ss_vfi)*sv(ia)/gam_ss_vfi + b_ss_j_vfi(bigJ) + bequest_ss_j_vfi(bigJ)), 1d-10)
                l_ss(bigj, ia, ip, ir, id) = 0d0
                svplus_ss(bigj, ia, ip, ir, id)=0d0
                V_ss(bigj, ia, ip, ir, id) = valuefunc(0d0, c_ss(bigj, ia, ip, ir, id), l_ss(bigJ,ia, ip, ir, id), bigJ, ip, ir, id)
            enddo
        enddo
    enddo
enddo
do ia=0, n_a, 1 
    do ip = 1, n_sp, 1
            do ir = 1, n_sr, 1
            do id = 1, n_sd, 1   
                EV_prim = 0d0
                EV_ss(bigj, ia, ip, ir, id)  = 0d0
                do ir_r= 1, n_sr, 1
                    do id_d= 1, n_sd, 1
                        do ip_p = 1, n_sp,1
                            if(theta == 1_dp)then
                                EV_prim = EV_prim + (1d0+r_ss_vfi+n_sr_value(ir_r))/gam_ss_vfi*pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)/c_ss(bigj, ia, ip_p, ir_r, id_d)    
                            else
                                EV_prim = EV_prim + (1d0+r_ss_vfi+n_sr_value(ir_r))/gam_ss_vfi*pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)*c_ss(bigj, ia, ip_p, ir_r, id_d)**(phi -theta*phi -1)
                            endif 
                            EV_ss(bigj, ia, ip, ir, id)  = EV_ss(bigj, ia, ip, ir, id) + pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)*V_ss(bigj, ia, ip_p, ir_r, id_d)
                        enddo
                    enddo
                enddo
                if(theta == 1)then
                    RHS_ss(bigj, ia, ip, ir, id) = 1d0/((delta+n_sd_value(id))*pi_ss_vfi_cond(bigJ)*EV_prim)
                else
                    RHS_ss(bigj, ia, ip, ir, id) = (delta+n_sd_value(id))*pi_ss_vfi_cond(bigJ)*EV_prim
                    EV_ss(bigj, ia, ip, ir, id)  = ((1d0-theta)*EV_ss(bigj, ia, ip, ir, id))**(1d0/(1d0-theta)) 
                endif
                    
            enddo
        enddo
    enddo
enddo

do j = bigJ-1, 1, -1
    poss_ass_sum_ss(j) = 0d0
        do i= j, bigJ, 1
            if(i < jbar_ss_vf)then
                poss_ass_sum_ss(j) = poss_ass_sum_ss(j) + ((w_pom_ss_vfi(i)*omega_ss(j)*n_sp_value(1)) + bequest_ss_j_vfi(i))/((1d0+n_sr_value(1)+r_ss_vfi)/gam_ss_vfi)**(i-j) 
            else
                poss_ass_sum_ss(j) = poss_ass_sum_ss(j) + (b_ss_j_vfi(i)                              + bequest_ss_j_vfi(i))/((1d0+n_sr_value(1)+r_ss_vfi)/gam_ss_vfi)**(i-j)              
            endif         
        enddo   
    do ia=0, n_a, 1
        do ip=1, n_sp,1
            do ir = 1, n_sr,1
                do id =1, n_sd,1
                    if((sv(ia)*(1d0+n_sr_value(ir)+r_ss_vfi)/gam_ss_vfi + poss_ass_sum_ss(j))  <a_l)then    
                        c_ss(j, ia, ip, ir, id) = 1d-10 
                        if(j < jbar_ss_vf)then
                            l_ss(j, ia, ip, ir, id) = 1d0 
                            lab_income = omega_ss(j)*n_sp_value(ip)*w_pom_ss_vfi(j)
                        else
                            l_ss(j,ia, ip, ir, id) = 0d0
                            lab_income = 0d0
                        endif
                        sv_tempo(j, ia, ip, ir, id) = (c_ss(j, ia, ip, ir, id)+sv(ia)-lab_income&
                                                                -b_ss_j_vfi(j)- bequest_ss_j_vfi(j))/((1d0+n_sr_value(ir)+r_ss_vfi)/gam_ss_vfi)
                    else 
                        if(j>=jbar_ss_vf) then ! retired thus labor choice is trivial 
                                l_ss(j, ia, ip, ir, id) = 0d0
                                lab_income = 0d0 
                                if(theta == 1)then ! consumption can be calculated diractly from RHS
                                    c_ss(j, ia, ip, ir, id) = max(RHS_ss(j+1, ia, ip, ir, id),1d-15)
                                else
                                    c_ss(j, ia, ip, ir, id) = max(RHS_ss(j+1, ia, ip, ir, id)**(1d0/(phi -theta*phi -1)),1d-15)
                                endif
                        else
                                wage            = omega_ss(j)*n_sp_value(ip)*w_pom_ss_vfi(j)
                                if(theta == 1)then
                                        c_ss(j, ia, ip, ir, id) = max(RHS_ss(j+1, ia, ip, ir, id),1d-15)
                                        c_opt = c_ss(j, ia, ip, ir, id) 
                                        l_ss(j, ia, ip, ir, id) = optimal_labor(c_opt, wage, phi)
                                else
                                        c_ss(j, ia, ip, ir, id) = max(((1-phi)/(phi*wage))**((1-theta)*(1-phi)/theta)*RHS_ss(j+1, ia, ip, ir, id)**(1d0/(-theta)),1d-15) !optimal_choice(1)
                                        c_opt = c_ss(j, ia, ip, ir, id) 
                                        l_ss(j, ia, ip, ir, id) = optimal_labor(c_opt, wage, phi)
                                endif
                                lab_income = wage*l_ss(j, ia, ip, ir, id)

                                    
                        endif   
                        sv_tempo(j, ia, ip, ir, id) = (c_ss(j, ia, ip, ir, id)+sv(ia)&
                                                                - lab_income-b_ss_j_vfi(j)&
                                                                - bequest_ss_j_vfi(j))/((1d0+n_sr_value(ir)+r_ss_vfi)/gam_ss_vfi )                
                    endif
                enddo
            enddo
        enddo
    enddo

    do ip=1, n_sp, 1
        do ir=1, n_sr, 1
            do id=1, n_sd, 1
                call change_grid_piecewise_lin_spline(sv_tempo(j,:, ip, ir, id), sv,   sv, svplus_ss(j,:, ip,ir, id))
            enddo
        enddo  
    enddo
        
     do ia=0, n_a, 1       
            do ip=1, n_sp, 1
                 do ir=1, n_sr, 1
                    do id= 1, n_sd, 1 
                        if(svplus_ss(j, ia, ip, ir, id)<a_l)then
                            svplus_ss(j, ia, ip, ir, id) = a_l
                        endif
                           available = (1d0+n_sr_value(ir)+r_ss_vfi)*sv(ia)/gam_ss_vfi + b_ss_j_vfi(j)+ bequest_ss_j_vfi(j) &
                                       - svplus_ss(j, ia, ip, ir, id)
                        if(j>=jbar_ss_vf) then
                            c_ss(j, ia, ip, ir, id) = max(available, 1e-10)
                            l_ss(j, ia, ip, ir, id)=0d0
                            lab_income = 0d0
                        else                     
                            wage =  w_pom_ss_vfi(j)*omega_ss(j)*n_sp_value(ip)
                            if (switch_fix_labor == 0) then 
                                foc = foc_intratemp(available, wage, 0.001d0)
                            else
                                foc = foc_intratemp(available, wage, switch_fix_labor)
                            endif
                            c_ss(j, ia, ip, ir, id) = foc(1)! max(phi*(available+wage), 1d-8) 
                            l_ss(j, ia, ip, ir, id) = foc(2) !optimal_labor(c_opt, wage, phi)
                            lab_income = wage*l_ss(j, ia, ip, ir, id)
                        endif
                        pi_com = pi_ss_vfi_cond(j)
                        V_ss(j, ia, ip, ir, id) = valuefunc(svplus_ss(j, ia, ip, ir, id), c_ss(j, ia, ip, ir, id), l_ss(j, ia, ip, ir, id), j,  ip, ir, id) 
                    enddo
                enddo
            enddo
   
            do ip=1, n_sp, 1
                do ir=1, n_sr, 1
                    do id =1, n_sd,1
                            EV_prim = 0d0
                            EV_ss(j, ia, ip, ir, id)  = 0d0
                            do ip_p=1, n_sp, 1
                                do ir_r=1, n_sr, 1
                                    do id_d=1, n_sd, 1  
                                        c_help = c_ss(j, ia, ip_p, ir_r, id_d) 
                                        l_help = l_ss(j, ia, ip_p, ir_r, id_d) 
                                        if(theta == 1_dp)then
                                            EV_prim = EV_prim + (1d0+r_ss_vfi+ n_sr_value(ir_r))/gam_ss_vfi*pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)*1/c_help
                                        else
                                            if(j<jbar_ss_vf)then !base  on D:\Dropbox (UW)\NCN EMERYT\__model\egm\CRRA
                                                EV_prim =  EV_prim + (1d0+r_ss_vfi+n_sr_value(ir_r))/gam_ss_vfi*pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)&
                                                                    *((1-l_help)/c_help)**((1d0-theta)*(1d0-phi))*c_help**(-theta)
                                            else
                                                EV_prim = EV_prim + (1d0+r_ss_vfi+n_sr_value(ir_r))/gam_ss_vfi*pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)&
                                                          *c_help**(phi -theta*phi -1)
                                            endif
                                        endif                                  
                                        EV_ss(j, ia, ip, ir, id)  = EV_ss(j, ia, ip, ir, id) + pi_ip(ip, ip_p)*pi_ir(ir, ir_r)*pi_id(id,id_d)&
                                                                           *V_ss(j, ia, ip_p, ir_r, id_d)
                                    enddo
                                enddo
                            enddo
                        if(theta==1_dp)then
                            RHS_ss(j, ia, ip, ir, id)=1d0/((delta+n_sd_value(id)) *pi_ss_vfi_cond(j)*EV_prim)  
                        else
                            RHS_ss(j, ia, ip, ir, id)= (delta+n_sd_value(id)) *pi_ss_vfi_cond(j)*EV_prim
                        endif 
                    
                        if (theta == 1) then 
                            EV_ss(j, ia, ip, ir, id) = EV_ss(j, ia, ip, ir, id)
                        else 
                            EV_ss(j, ia, ip, ir, id) = ((1d0-theta)*EV_ss(j, ia, ip, ir, id))**(1d0/(1d0-theta))
                        endif
                    enddo
               enddo  
            enddo
    enddo
enddo
end subroutine