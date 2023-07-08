! WHAT   :  read initial data and inital values from world without the reform
! TAKE   :  data files and output files from base scenario (without reform) 
! DO     :  read data from files to variables and parameters 
! RETURN :  base variable CRUCIAL to the next run on the path 

MODULE get_data
use global_vars
IMPLICIT NONE
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_data(omega_ss_d, pi_d, pi_weight_d)
      real(dp), dimension(bigJ), intent(out) :: omega_ss_d
      real(dp), dimension(bigJ,2), intent(out) ::  pi_d, pi_weight_d

      call chdir(cwd_r)

! -------------------------------- age specyfic productivity -------------------------------
     OPEN (unit=3, FILE = "_data_omega_ok.txt")    
       do j = 1, bigJ, 1
        read(3,*) omega_ss_d(j)
      end do
    close(3)
    
! -------------------------------- population -------------------------------
    ! UNCONDITIONAL SURVIVAL PROBABILITIE 
    OPEN (unit=3, FILE = "_data_pi.txt")   
    do i = 1, 2, 1
       do j = 1, bigJ, 1
        read(3,*) pi_d(j,i)
       end do
    enddo
    close(3)

    ! pi_weight  is needed to calculate relative masses in steady states
    pi_weight_d = pi_d

end subroutine read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module  get_data