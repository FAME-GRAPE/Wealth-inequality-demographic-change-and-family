    
    if ( switch_run_1 == 1) then 
    no_steady = 'ss1_'    
    else
    no_steady = 'ss2_'
    endif
    
    
    
    
    ! write on the screen
        write(*,*) '******* STEADY STATE PAYG *******'
        write(*,*)
        write(*,*) '*********************************'
        write(*,*) 'Calibration:'
        write(*,*)
        
        ! most of it has to be removed, we do not care about these... 
        if  (bigJ == 16)  then ! US
            write(*,'(A25,F10.7,A)') ' 100*sum_b/y =  ', 100*sum_b_ss/y_ss, '  |  Should be 5.2%'
        else
            write(*,'(A25,F10.7,A)') ' 100*sum_b/y = ', 100*sum_b_ss/y_ss, '  |  Should be 5%'
        endif
        if (bigJ == 16)  then
            write(*,'(A25,F10.7,A)') ' average hours =  ', 100*sum(N_ss_j(1:jbar_ss-1)*lab_ss_j_vfi(1:jbar_ss-1))/sum(N_ss_j(1:jbar_ss-1)), '  |  Should be 33%' 
        else
            write(*,'(A25,F10.7,A)') ' average hours = ', bigl_ss/sum(N_ss_j(1:jbar_ss-1)), '  |  Should be 56.8%' 
        endif
            if (bigJ == 4) then
                write(*,'(A20,F10.7,A)') ' 1 + r_bar_ss = ', (1 + r_bar_ss)**0.05_dp, '  |  Should be 7.8%'
                write(*,'(A20,F10.7,A)') ' r_ss = ', r_ss**0.05_dp
            elseif (bigJ == 20)  then   
                write(*,'(A20,F10.7,A)') ' 1 + r_bar_ss = ', (1 + r_bar_ss)**0.25_dp, '  |  Should be 7.8%'
                write(*,'(A20,F10.7,A)') ' r_ss = ', r_ss**0.25_dp
           elseif (bigJ == 16)  then   ! US
                write(*,'(A25,F10.7,A)') ' 1 + r_bar_ss = ', 100*((1 + r_bar_ss)**0.2_dp -1d0), '  |  Should be 5.5%'
                write(*,'(A25,F10.7,A)') ' r_ss = ', r_ss**0.2_dp
            else
                write(*,'(A20,F10.7,A)') ' 1 + r_bar_ss = ', 1 + r_bar_ss, '  |  Should be 7.8%'
                write(*,'(A20,F10.7,A)') ' r_ss = ', r_ss
                write(*,'(A20,F10.7,A)') ' investment rate = ', (y_ss - consumption_ss_gross)/y_ss, '  |  Should be 21%'
                write(*,'(A20,F10.7,A)') ' investment rate = ', ((gam_ss+depr-1)*k_ss)/y_ss
  
            endif       
        write(*,'(A20,F10.7,A)') 'gam_ss*nu_ss = ',   gam_ss*nu_ss
        write(*,'(A20,F10.7,A)') 'r= ',   1 + r_bar_ss
        
        write(*,*) '*********************************'
        write(*,*)
        write(*,'(A30,F10.7,A)') ' err_ss = ', err_ss
        write(*,'(A30,F10.7,A)') ' k_ss = ', k_ss
        write(*,'(A30,F10.7,A)') ' y_ss = ', y_ss
        write(*,'(A30,F10.7,A)') ' K_ss = ', k_ss*bigl_ss/N_ss_j(1)
        write(*,'(A30,F10.7,A)') ' L_ss = ', bigl_ss/N_ss_j(1)
        write(*,'(A30,F10.7,A)') ' w_bar_ss = ', w_bar_ss
        write(*,'(A30,F10.7,A)') ' capital output ratio = ', k_ss/(y_ss/zbar)


        write(*,'(A30,F10.7,A)') ' capital labour ratio = ', k_ss/bigl_ss
        write(*,'(A30,F10.7,A)') ' labour share = ', w_bar_ss/y_ss
        write(*,'(A30,F10.7,A)') ' capital share = ', ((r_bar_ss + depr)*k_ss/y_ss)
        write(*,'(A30,F10.7,A)') ' l_ss_pen_j(jbar-1) = ', l_ss_pen_j(jbar_ss-1)
        write(*,'(A30,F10.7,A)') ' l_ss_j(jbar-1) = ', l_ss_j(jbar_ss-1)
        write(*,*) ' u_ss = ', V_ss_j_vfi(1)
        write(*,'(A30,F16.7,A)') ' N_ss = ', N_ss
        write(*,'(A30,F16.7,A)') ' bigl_ss = ', bigl_ss
        write(*,'(A30,F16.7,A)') ' bequest_ss = ', bequest_ss
        write(*,'(A30,F16.7,A)') ' savings_top_ten = ',  savings_top_ten(10)/sum(N_ss_j*asset_pom_ss_j)
        write(*,'(A30,F16.7,A)') ' savings_top_100 = ',  savings_top_100(100)/sum(N_ss_j*asset_pom_ss_j)
        write(*,'(A30,F16.7,A)') 'private_wealth/y_ss ratio =', sum(N_ss_j*asset_pom_ss_j)/y_ss
        write(*,'(A30,F16.7,A)') 'share negative assets =', share_neg
        write(*,'(A30,F16.7,A)') 'share nonpositive assets =', share_nonpos
        write(*,'(A30,F16.7,A)') 'borrowing limit to LabIncAVG_ss_vfi', sv(0) / LabIncAVG_ss_vfi
        write(*,'(A30,F16.7,A)') ' top_ten = ',  top_ten(10)
        write(*,'(A30,F16.7,A)') ' phi = ',  phi
        write(*,'(A30,F10.7,A)') ' delta = ', delta
        write(*,'(A30,F16.7,A)') ' t1 = ',  t1_ss
        write(*,'(A30,F16.7,A)') ' LabIncAVG_ss_vfi = ',  LabIncAVG_ss_vfi
        write(*,'(A30,F16.7,A)') ' replacement = ',  replacement_ss
        write(*,'(A30,F16.7,A)') ' feasibility = ', abs((y_ss - consumption_ss_gross)/y_ss - ((nu_ss*gam_ss+depr-1)*k_ss)/y_ss) 
        write(*,*) '********************************************'

    OPEN (unit=666, FILE = version//experiment//no_steady//"aggregates.csv")

 
    write(666, '(A)') "Outcomes"
    write(666, '(A)') "y;k/y;c/y;i/y;bigl;r;beq/y;gam_ss;average hours;r-g;replacement;"
    write(666, '(F20.10,A)', advance='no') y_ss, ";"
    write(666, '(F20.10,A)', advance='no') k_ss/y_ss, ";"
    write(666, '(F20.10,A)', advance='no') consumption_ss_gross/(y_ss), ";"
    write(666, '(F20.10,A)', advance='no') ((gam_ss+depr-1)*k_ss)/y_ss, ";"
    write(666, '(F20.10,A)', advance='no') bigl_ss, ";"
    write(666, '(F20.10,A)', advance='no') r_ss, ";"
    write(666, '(F20.10,A)', advance='no') bequest_ss/y_ss, ";"
    write(666, '(F20.10,A)', advance='no') gam_ss, ";"
    write(666, '(F20.10,A)', advance='no') bigl_ss/sum(N_ss_j(1:jbar_ss-1)), ";"
    write(666, '(F20.10,A)', advance='no') r_ss-gam_ss, ";"
    write(666, '(F20.10,A)', advance='no') replacement_ss, ";"

    write(666, '(A)') ""
    write(666, '(A)') "Lifecycle"
    write(666, '(A)') "yr;c;l;s;V;disc"
    do j = 1, bigJ
        write(666, '(I2,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10)') j, ";", c_ss_j(j), ";", l_ss_j(j), ";", s_ss_j(j), ";", V_ss_j_vfi(j), ";", delta**(j-1)*(pi_ss(j)/pi_ss(1))
    enddo
CLOSE(666)

    
if (switch_ss_write == 1) then
    open(unit = 106, FILE = version//experiment//no_steady//"gini_weight.csv")
    write(106, '(A)') "weight;age;assets"
            do j = 1, bigJ, 1
                do ia = 0, n_a, 1
        write(106, '(F20.10,A,I5,A,F20.10)') &
                    gini_weight_sv(j,ia), ";", & ! weight
                    j , ";",  & !age
                    sv(ia)  !assets
                enddo        
            enddo
    close(106)

    open(unit = 107, FILE = version//experiment//no_steady//"prob.csv")
    write(107, '(A)') "prob;age;asset;inc_shock;ret_shock;disc_shock"
   
            do j = 1, bigJ, 1
                do ia = 0, n_a, 1
                    do ip = 1, n_sp, 1
                        do ir = 1, n_sr, 1
                            do id = 1, n_sd, 1
                            write(107, '(F20.10,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5)') &
                            prob_ss(j, ia, ip, ir, id), ";", & ! probability
                            j , ";",  & !age
                            ia , ";",  & !asset
                            ip , ";",  & !income
                            ir , ";",  & !return
                            id  !discount
                            enddo        
                        enddo
                    enddo
                enddo
            enddo


    close(107)

    open(unit = 108, FILE = version//experiment//no_steady//"mass.csv")
    write(108, '(A)') "mass;cons;hours;labinc;labinc_pretax;age;asset;inc_shock;ret_shock;disc_shock"
        do i = 1, bigT, 1 !todo18
            do j = 1, bigJ, 1
                do ia = 0, n_a, 1
                    do ip = 1, n_sp, 1
                        do ir = 1, n_sr, 1
                            do id = 1, n_sd, 1
                            write(108, '(F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,I5,A,I5,A,I5,A,I5,A,I5)') &
                            prob_ss(j, ia, ip, ir, id)*n_ss_j(j)/sum(n_ss_j(:)), ";", & ! mass
                            c_ss(j, ia, ip, ir, id), ";", & !consumption
                            l_ss(j, ia, ip, ir, id), ";", & !hours
                            lab_income_ss(j, ia, ip, ir, id), ";", & !lab income
                            lab_income_pretax_ss(j, ia, ip, ir, id), ";", & !lab income
                            j , ";",  & !age
                            ia , ";",  & !asset
                            ip , ";",  & !income
                            ir , ";",  & !return
                            id  !discount
                            enddo        
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endif
    