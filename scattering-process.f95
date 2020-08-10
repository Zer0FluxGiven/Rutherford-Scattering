!------------------------------------------------------------
! INTRODUCTION
!------------------------------------------------------------

! This program simulates Rutherford Scattering in 3 Dimensions.

!Written by Andrew Murphy, 2020


!-----------------------------------------------------------
!	FUNCTIONS
!-----------------------------------------------------------
!This function produces the Hamiltonian for a repulsive Kepler-Force (H = p**/2m + k/r) 
!we use this function to calculate the energy of our particle
function H(pc, k, m)
implicit none
real :: k, m, H
real, dimension(6) :: pc
H = (0.5/m)*(pc(4)*pc(4) + pc(5)*pc(5) + pc(6)*pc(6)) + k/sqrt(pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))
end function H

!The following section is a work in progress -- my naivite in FORTRAN has reuslted in several compiler errors wtih these array-valued functions
!Hopefully I figure out whats going wrong here and  I can incorperate these functions into the main code--Hopefully...


!This function returns the cross-product of two 3D vectors
!real, dimension(3) function cross(a,b) 
!real, dimension(3) :: a, b
!product = (/ (a(2)*b(3) - a(3)*b(2)), (a(3)*b(1) - a(1)*b(2)) , (a(1)*b(2) - a(2)*b(1)) /)
!end function cross

!This function returns the 3D Angular Momentum vector of a particular phase-space coordinate (L = r x p)
!real, dimension(3) function L(pc) 
!real, dimension(3) :: L, cross
!real, dimension(6) :: pc
!ang_mom = cross((/ pc(1), pc(2) , pc(3)/) , (/ pc(4), pc(5), pc(6) /))
!end function L

!This function returns the 3D Laplace-Runge-Lenz vector of a particular phase-space coordinate
!real, dimension(3) function A(pc, k, m) 
!real :: k, m, r
!real, dimension(3) :: L, LRL
!real, dimension(6) :: pc
!r = sqrt(pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))
!LRL = cross( (/ pc(4), pc(5), pc(6) /) , L(pc)) - (m*k/r)*(/ pc(1) , pc(2) , pc(3) /)
!end function A

!This function returns the magnitude of the difference between two 3D vectors.
!real function mag_diff(a, b)
!real :: mag_diff
!real, dimension(3) :: a, b
!mag_diff = sqrt((b(1) - a(1))**2 + (b(2) - a(2))**2 + (b(3) - a(3))**2)
!end function mag_diff



!----------------------------------------------------------
!	SUBROUTINES
!----------------------------------------------------------
!This subroutine represents Hamilton's Equations for a repulsive Kepler-Force.
!'pc' is the 6-dimesional phase-space coordinate (x, y, z, px, py, pz) of the particle
!'pq-dot' contains  the time-derivatives of the phase-space coordinates, given by Hamilton's Equations (dH/dx = -dpx/dt; dH/dpx = dx/dt)
subroutine Hamiltons_eqs(pc, pq_dot, k, m)
implicit none
real :: k, m
real, dimension(6):: pc, pq_dot
pq_dot(1) = pc(4)/m
pq_dot(2) = pc(5)/m
pq_dot(3) = pc(6)/m
pq_dot(4) = k*pc(1)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5) 
pq_dot(5) = k*pc(2)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5) 
pq_dot(6) = k*pc(3)/((pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3))**1.5) 
end subroutine Hamiltons_eqs

!---------------------------------------------------
!	MAIN CODE
!---------------------------------------------------

program rscattering
implicit none
real, parameter :: k=100.0, m=1.0, s_max =1.0, z_0 = -100.0, p_0=100.0, pi=3.1416 
integer, parameter:: N_particles=10000 ,N_steps=10000
!our parameters include: 'k' = force-constant (ZZ'e^2 for nuclear scattering), 'm' = mass
!'s_max' = beam-radius (maximum impact parameter), 'z_0' = initial z-coordinate
!'p_0' = initial z-momentum, and 'pi' . Note: since pi is only used to generate random initial coordinate, a precise value is not too neccessary
!'N_particles' = number of particles to be simulated, 'N_steps' = number of time-steps in our Euler-Method Algortihm

real:: s_0, phi_0, H, dt
real, dimension(6):: pc_0, pc, pq_dot
real, dimension(N_particles,5) :: data_vals
integer :: i, j

dt = abs((4*m*p_0/(z_0))/N_steps) !our maximum time is abs(4m*p_0/z_0), which is ""sufficeiently large"" to capture the scattering dynamics 

do i=1,N_particles
	s_0 = s_max * rand()											!Generate random impact parameter
	phi_0 = 2*pi * rand()											!Generate random azimuthal angle 
	pc_0 = (/ s_0*cos(phi_0), s_0*sin(phi_0), z_0, 0.0, 0.0, p_0 /)	!create initial phase-space coordinate
	pc = pc_0											!record initial impact-parameter
	
	do j=1,N_steps								!Begin Euler-Method Algortihm for approximating the trajectory 
		call Hamiltons_eqs(pc, pq_dot, k, m)	!Call subroutine to generate phase-space differentials
		pc = pc + pq_dot*dt						!Increment the phase-space coordinate of the particle
	end do
	
	data_vals(i,1) = acos(pc(3)/sqrt(pc(1)*pc(1) + pc(2)*pc(2) + pc(3)*pc(3)))
	data_vals(i,2) = s_0
	data_vals(i,3) = H(pc, k, m) - H(pc_0, k, m)
!	data_vals(i,4) = mag_diff(L(pc_0),L(pc))
!	data_vals(i,5) = mag_diff(A(pc_0,k,m),A(pc,k,m))
end do	

open(10, file="data\theta_vals.txt", status='replace')	!list of scattering angles
open(11, file="data\s_vals.txt", status='replace')		!list of initial impact parameters
open(12, file="data\Ediff_vals.txt", status='replace')	!list of energy-drifts
!open(13, file="data\Ldiff_vals.txt", status='replace')	!list of angular momentum-drifts
!open(14, file="data\Adiff_vals.txt", status='replace')	!list of LRL-drifts

do i=1, N_particles				!Write data to files
	write(10,*) data_vals(i,1)
	write(11,*) data_vals(i,2)
	write(12,*) data_vals(i,3)
!	write(13,*) data_vals(i,4)
!	write(14,*) data_vals(i,5)
end do
end program rscattering