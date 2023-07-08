!|==============================================================================================================================================|
!| module solove consumer problem (given rbar, wbar, taxes and other parameters) using policy iterations function, Right Hand Side (RHS) method |                         
!|functions: euler_trans, euler_transsec, euler_ss, eulersec_ss - finding optimal futur assets for given grid point                             |
!|function margu - marginal utility for consumption                                                                                             |
!|function year - find year in which consumer will have ijj years, having in year it, ij years                                                  |    
!|subroutines: linear_int - linear evaluation assets on the asset grid point                                                  |
!|subroutines: household_ss, household_trans - find futur assets for every age, assets grid point                                               |
!|subroutines: interpolate_ss, interpolate_trans - interpolate RHS                                                                              |
!|subroutines: get_distribution_ss, get_distribution_trans - get distribution of consumers on age, gridpoints                                   |    
!|subroutines: aggregation_ss, aggregation_trans - agreggate consumption, labour, savings                                                       |
!|subroutine output - print to CSV                                                                                                              |
!|==============================================================================================================================================|
module pfi_trans
    ! modules
    use assertions
    use errwarn
    use gaussian_int
    use linint
    use matrixtools
    use minimization
    use normalProb
    use polynomial
    use simplex
    use clock
    use rootfinding
    use splines
    use AR_discrete
    use global_vars
    
    
implicit none

!definition of variables
integer ::  ia, ip, ir, id, ir_r, ip_p, id_d
real*8  :: sv(0:n_a)
integer ::  n_a_1, n_a_2

!steady state variables
real*8, dimension(bigJ, 0:n_a, n_sp, n_sr, n_sd) :: V_ss, EV_ss, RHS_ss,  svplus_ss, l_ss, c_ss, lab_income_ss, lab_income_pretax_ss, sv_tempo,  prob_ss, gini_weight_consumption

real*8, dimension(bigJ) :: V_ss_j_vfi, c_ss_j_vfi, s_pom_ss_j_vfi, l_ss_j_vfi, lab_ss_j_vfi, l_ss_pen_j, b_ss_j_vfi, &
                           bequest_ss_j_vfi, pi_ss_vfi, pi_ss_vfi_cond, &
                           lw_ss_j_vfi, w_pom_ss_vfi, lab_high 

real*8 ::   r_ss_vfi,  LabIncAVG_ss_vfi, gam_ss_vfi, bequest_ss_vfi 
real*8 ::   gini_weight_sv(bigJ, 0:n_a)
integer ::  jbar_ss_vf

!comunication variables
integer :: j_com, ia_com, ip_com, i_com, ir_com
real*8 :: c_com, l_com, pi_com, c_ss_com, l_ss_com

! top_ten
real*8, dimension(bigJ) :: N_ss_j_vfi, asset_pom_ss_j, top_ten_coh, cons_proc_top_ten_coh 
real*8, dimension(bigJ, 0:n_a, n_sp, n_sr,n_sd) :: gini_income
real*8 :: savings_top_ten(10), top_ten(10), savings_cohort_ten(3,bigJ), &
          consumption_top_ten(3, bigJ), top_100(100),  savings_top_100(100), &
          top_ten_trans(10,bigT), savings_top_ten_trans(10, bigT), wspl(bigT), l_pen_j(bigJ,bigT), asset_trans(bigJ,bigT) , &
          gini_weight_trans(bigJ,0:n_a, bigT), share_neg, share_nonpos
          
integer :: t
contains 


!*******************************************************************************************
! marginal utility of consumption
    function margu(cons, labor, tc)
        real*8, intent(in) :: cons, labor, tc !cosumption , labor and consumption tax
        real*8 :: margu, leis
        leis=1d0-labor
        margu = 1/tc*phi*(cons**phi*leis**(1d0-phi))**(1d0-theta)/cons
    end function 

!*******************************************************************************************
! value function - steady state!
 
    function valuefunc(sv_plus, cons, lab, j, ip, ir, id)
        implicit none

        integer, intent(in) :: j, ip, ir, id
        real*8, intent(in) :: sv_plus, cons, lab
        real*8 :: valuefunc, dist, c_help, l_help, pr
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)  
        l_help = min(max(lab,1d-10), 1d0-1d-10)
        ! get tomorrows utility
        call linear_int(sv_plus, ial, iar, dist, sv, n_a, a_grow)

        
        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(j < bigJ)then
            if (theta == 1) then
                !valuefunc = log(max(dist*EV_ss(j+1, ial, ip,ir, id) + &
                              !(1d0-dist)*EV_ss(j+1, iar, ip,ir, id), 1d-10))
                valuefunc =   dist*EV_ss(j+1, ial, ip, ir, id) + &
                              (1d0-dist)*EV_ss(j+1, iar, ip, ir, id) 
            else
                valuefunc =     max(dist*EV_ss(j+1, ial, ip, ir, id) + &
                              (1d0-dist)*EV_ss(j+1, iar, ip, ir, id), 1d-10)**(1d0-theta)/(1d0-theta)
            endif
        endif    
        ! add todays part and discount
        if (theta == 1) then
            valuefunc = (phi*log(c_help) + (1d0-phi)*log(1d0-l_help)) +  (delta+n_sd_value(id))*pi_com*valuefunc 
        else
            valuefunc = (c_help**phi*(1d0-l_help)**(1d0-phi))**(1d0-theta)/(1d0-theta) + (delta+n_sd_value(id))*pi_com*valuefunc 
        endif
    end function

        
  !*******************************************************************************************  
    ! find optimal labor for given consumption, wage, preference for leisure and labor tax parameters
     function optimal_labor(c, w, phi)
     
        real *8 :: c, w, phi
        real *8 :: l
        real *8 :: optimal_labor
        real *8 :: c_mult

        c_mult = (1-phi)/phi
        l = 1 - c*c_mult/w
  
        optimal_labor = min(1d0,max(l,0d0))
     end function
    
!*******************************************************************************************

    subroutine agent_vf()

        implicit none

        integer :: iter
        
            call initialize_trans
                
            ! solve the household problem
            call household_endo()
    
            ! calculate the distribution of households over state space
            call get_distribution_ss()
    
            ! aggregate individual decisions
            call aggregation_ss()
          
    end subroutine

!   
!*******************************************************************************************
! adjust policy function for endogenous grid to exogenous one by linear interpolation
subroutine change_grid_piecewise_lin_spline(endo_grid, egzo_grid, f, g)
    real*8 :: endo_grid(0:n_a), egzo_grid(0:n_a), f(0:n_a), g(0:n_a), coef1, coef2
    integer :: first, last, ik, it
        
do ia = 0, n_a, 1 
    first=0
    last=n_a
    do while(abs(last-first)>2)
        ik=floor(fi*first+(1d0-fi)*last)
        it=floor((1d0-fi)*first+fi*last)
        if(endo_grid(ik)<endo_grid(it)) then
            if(endo_grid(it)<=egzo_grid(ia)) then
                first=it
            else if(endo_grid(ik)<=egzo_grid(ia) .and. egzo_grid(ia)<endo_grid(it)) then
                first=ik
                last=it
            else if(egzo_grid(ia)<endo_grid(ik)) then
                last=ik
            endif
        else
             if(endo_grid(ik)<=egzo_grid(ia)) then
                first=ik
            else if(endo_grid(it)<=egzo_grid(ia) .and. egzo_grid(ia)<endo_grid(ik)) then
                first=it
                last=ik
            else if(egzo_grid(ia)<endo_grid(it)) then
                last=it
            endif
        endif
    enddo
    ik=first
    if(last-first == 2 .and. endo_grid(first+1)<egzo_grid(ia)) then
        ik=ik+1
    endif
    if(ik==0 .and. endo_grid(ik)>egzo_grid(ia)) then
        g(ia)=f(0)
    else
        coef1=abs(endo_grid(ik+1)-egzo_grid(ia))/abs(endo_grid(ik+1)-endo_grid(ik))
        coef2=abs(endo_grid(ik)-egzo_grid(ia))/abs(endo_grid(ik+1)-endo_grid(ik))
        g(ia)=coef1*f(ik)+coef2*f(ik+1)
        if(g(ia)>egzo_grid(n_a)) then 
            g(ia) = f(n_a)
        endif
    endif
enddo    
end subroutine

!******************************************************************************************
    subroutine linear_int(assets, ial, iar, dist, grid, grid_size, grid_grow)

        implicit none
        integer :: grid_size
        real*8 :: assets, grid_grow
        real*8 :: grid(0:grid_size)
        integer :: ial, iar
        real*8 :: dist, real_ia

        real_ia = grid_Val_Inv(assets, grid(0), grid(grid_size), grid_size, grid_grow)
        
        ial = floor(real_ia)
        ial = max(ial, 0)
        ial = min(ial, grid_size-1)
        iar = ial+1
        dist = 1d0 - (assets-grid(ial))/(grid(iar)-grid(ial))

    end subroutine
!******************************************************************************************
    
   subroutine initialize_trans() 
    implicit none 
        
    real*8 :: sv_pom(0:n_a), Brackets(2)
   

    pi_ss_vfi_cond = 1d0 !todo 
    
    do i = 1 , bigT, 1    
        do j=2, bigJ        

            pi_ss_vfi_cond(j) =  pi_ss_vfi(j)/pi_ss_vfi(j-1) 
        enddo
    enddo 
    ! size of the asset grid
    sv(0:n_a) = grid_cons(a_l, a_u, n_a, a_grow)
    
   end subroutine
   
!*******************************************************************************************
function foc_intratemp(av, w_tax, l_guess)
implicit none
    real*8 :: av, w_tax, l_guess
    real*8 ::  c, l , foc_intratemp(2)
    real*8:: l0, del, taxinc, res, lp, tax, mrate, dres
    integer :: maxit, i

    !besed on Heathcote, Jonathan, Kjetil Storesletten, and Giovanni L. Violante. "Optimal tax progressivity: An analytical framework." The Quarterly Journal of Economics 132.4 (2017): 1693-1754
    ! https://www.nber.org/papers/w19899.pdf
    ! https://www.sas.upenn.edu/~dkrueger/research/LafferCurves.pdf

    ! TAKE a look at ncn emeryt\model\prog_income_tax.lyx
      maxit  = 10
      l0  = l_guess
      del = 1d-8
      if (switch_fix_labor == 0d0)  then 
          do i= 1, maxit
            taxinc = w_tax*l0
            res    = (1-phi)/phi*(av+taxinc)-w_tax*(1.0d0-l0)
            lp     = l0+del
            taxinc = w_tax*lp
            dres   = ((1-phi)/phi*(av+taxinc)-w_tax*(1.0d0-lp)-res)/del
            l0     = min(max(l0-res/dres, del), 1d0-del)
          enddo  
    endif  
      l = l0
      taxinc = w_tax*l0
      c = max(av+taxinc, del)
      if ((c .NE. c) .or. (l .NE. l) ) then
         write(*,*) 'problem'  
      endif
      foc_intratemp(1) = c
      foc_intratemp(2) = min(max(l, 0d0), 1d0-del)

endfunction 

!*******************************************************************************************

include "pfi_household_problem.f90"

include "pfi_distribution.f90"
 
include "pfi_agregation.f90"

include "pfi_print.f90"

end module