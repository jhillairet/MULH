module constants

use precision_def
implicit none
real (long) , parameter :: c = 299792458.0_long			! units m s-1
real (long) , parameter :: mu_0 = 1.25663706143591729e-6_long	! units m kg s^-2 A^-2
real (long) , parameter :: eps_0 = 8.85418781762038985e-12_long	! units m^-3 kg^-1 s^4 A^2
real (long) , parameter :: pi = 3.14159265358979323_long
real (long) , parameter :: e = -1.60217646e-19_long		! units C
real (long) , parameter :: me = 9.10938188e-31_long		! units kg
real (long) , parameter :: eme = e/me				! units C kg^-1

end module constants
