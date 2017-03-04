module def_types

use precision_def

! Info about SEE model
type sey_t
   integer Emax, seec, ReRr, E1, E_0, E_0p
   real(long) deltamax, kse, ks, Eom, p_n, a_lara, delta_b, z_lara
end type

! Info about geometry of waveguide and mesh
type geo_t
   real(long) a, b, lw, R_max, d
   integer NOC_PML, NOC(3), m_PML
end type

! Info about seed particles
type seed_t
   integer Np, Nt, vth, fmax, nsamplev
   real(long) vra
end type 

! Info about particle simulated
type particle_t
   integer pcol(4)
   real(long) px1(3), pv1(3), pv2(3), pB0(3), pB1(3), pB2(3), pE1(3), pE2(3)
end type

! Saving outputs
type output_s
  integer :: psave
end type
   

end module def_types
