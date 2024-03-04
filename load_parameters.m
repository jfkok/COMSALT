function load_parameters(d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function initializes various global parameters       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global u_fr_thr_soil P_reb alpha_ej_fract z_s_frac; %these are global parameters
global rho_p pi rho mu g P tau beta_time mu_angle mu_alpha_R sigma_alpha_R surf_angle a_splash min_lift;
global gamma init_mu_angle spin_sigma spin_mu sigma_W_prop sigma_U_prop elevation d_mean beta cohesion;

%Physical parameters
rho_p = 2650; %density of mineral sand/dust in kg/m3
P_reb = 0.96; %the probability that a high-speed particle will rebound, after Anderson and Haff (1991)
gamma = 1; %sets the constant in s/m with which the no-rebound probability decreases with increasing particle speed, after Andersons and Haff (1991)
alpha_ej_fract = 2.5; %the fraction of the impacting momentum not spent on the rebounding grain which is spent on ejecting surface particles
a_splash = 0.020; %this is the proportionality constant relating the number of ejected surface particles and the impact speed
tau = 600; %the electric relaxation time of the atmosphere in seconds, from Roble and Tzur (1986)
pi = 3.1415926536; %the ratio of the circumference to the diameter of a circle
sigma_W_prop = 1.4; %the proportionality between the standard deviation of the vertical wind speed and the friction speed.
sigma_U_prop = 2.5; %the proportionality between the standard deviation of the horizontal wind speed and the friction speed.
beta_time = 1; %the ratio of the Lagrangian and Eulerian time scales.  Values in the literature range from ~0.3 - 3 (Anfossi et al., 2006).
z_s_frac = 0.95; %the fraction of mass flux that occurs below the height z_s, which is used to compute the Lagrangian time scale based on the analogous flow of saltation layers and canopies [Raupach, 1991; Raupach et al., 1996; Leuning et al., 2000]
mu_angle = 40; %the average angle in degrees with which particles rebound from the surface 
init_mu_angle = 50; %the average angle in degrees with which particles are lifted/ejected from the surface
mu_alpha_R = sqrt(0.45); %the average restitution coefficient alpha_R
sigma_alpha_R = sqrt(0.22); %the standard deviation of the restitution coefficient alpha_R (alpha_R itself is calculated following the Word document "Theory of spherical particles colliding with bed of similar particles"
spin_sigma = 2*pi*500; %the standard deviation of the average spin (perpendicular to the xz plane) upon leaving the bed, after experiments
spin_mu = 2*pi*400; %the average spin (perpendicular to the xz plane) upon leaving the bed, after experiments
surf_angle = 0; %the angle in radians between the surface plane and the plane perpendicular to gravity

%physical constants of nature:
g = 9.8; %gravitational acceleration in m/s2
T = 300;  %the temperature in Kelvin
mu = 0.00751*(T/291.15)^(1.5)/(T+120); %the viscosity in kg/(m s) from the Sutherland relation
m_molar = 0.0289; %molar weight of a mol of air
P0 = 101325;  %Sea-level pressure in Pa
scale_height = 8.31*T/(m_molar*g); %the scale height of the atmosphere in m
P = P0*exp(-elevation/scale_height);  %the pressure in Pascal assuming hydrostatic balance
rho = P/(8.31447*T)*m_molar; %air density in kg/m3 from the ideal gas law

%the following sets beta, which scales the cohesive forces on soil particles, which affects the splashing of soil particles
if (cohesion) %checking if the cohesion switch is on
    if (Mars) %checking if the Mars switch is on
        beta = 6*10^(-5); %this scales the interparticle forces for Martian soils after Kok (2010, GRL)
    else %if the Mars switch is off, the simulation is for Earth
        beta = 1.5*10^(-4); %this scales the interparticle forces for dry Earth soils after Kok and Renno (2006, GRL)
    end %if, checking if the cohesion switch is on
else %in this case the cohesion switch is off, so beta is set to zero
    beta = 0; %setting beta to zero
end %if, checking if the cohesion switch is on

%the threshold friction velocity of the soil, following Shao and Li (2000)
u_fr_thr_soil = 0.10*sqrt(((rho_p-rho)/rho)*g*d_mean); %the fluid threshold for loose sand, after Bagnold (1941)
u_fr_thr_d = sqrt(0.0123*rho_p/rho)*sqrt(g*d+6*beta./(pi*rho_p*d)); %the fluid threshold for each particle bin
[u_fr_thr_min,min_lift] = min(u_fr_thr_d); %this determines the particle bin that is most easily lifted, which is used for initializing the model
