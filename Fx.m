function Fx = Fx(Cd,d,V_rel_x,V_rel_z,V_rel_mag,spin,dU)

%Fx calculates the horizontal force (i.e., parallel to the surface)

%these are global parameters
global rho rho_p pi mu Saffman g surf_angle;

Fd = -(pi/8)*(d^2)*V_rel_mag*rho*Cd*V_rel_x;  %the aerodynamic drag force

omega_p = spin*d/V_rel_mag; %a parameter used to calculate the aerodynamic spin force
Re_p = V_rel_mag*d*rho/mu; %the particle reynolds number
Cl_omega = 1-(0.675+0.15*(1+tanh(0.28*(omega_p-2))))*tanh(0.18*sqrt(Re_p)); %the normalized spin coefficient, after Loth (2008), eq. 16
Fmag = -(pi/8)*d^3*rho*Cl_omega*spin*V_rel_z; %the Magnus lift force due to particle rotation
F_grav = sin(surf_angle)*g*((rho_p-rho)/rho_p); %the projection of the gravitational force along the surface (= 0 for a flat surface)

if (Saffman==false) %checking whether the Saffman force should be calculated, which is computationally intensive and has little effect on trajectories
    Fx = Fd + Fmag + F_grav; %the total force. In this case the Saffman force, which is generally very small, is neglected
else
    %the following calculates the function J, which is necessary to compute the Saffman lift force that arises from 
    %shear in the flow.  After McLaughlin(1991), table 1 and equations 3.26 and 3.27
    epsilon = sqrt(mu*abs(dU)/(rho*V_rel_mag^2)); %a parameter used in calculating the Saffman force
    if (epsilon == 0) %calculating the Saffman force based on the range in which epsilon lies
        J = 0;
    elseif (epsilon < 0.025)
        J = -32*pi^2*(epsilon^5)*log(epsilon^(-2));
    elseif (epsilon < 0.05)
        J = -0.0000133+(epsilon - 0.025)*(-0.0002845+0.0000133)/0.025;
    elseif (epsilon < 0.1)
        J = -0.0002845+(epsilon - 0.05)*(-0.004658+0.0002845)/0.05;
    elseif (epsilon < 0.15)
        J = -0.004658+(epsilon - 0.1)*(-0.01458+0.004658)/0.05;
    elseif (epsilon < 0.2)
        J = -0.01458+(epsilon - 0.15)*(-0.01247+0.01458)/0.05;
    elseif (epsilon < 0.25)
        J = -0.01247+(epsilon - 0.20)*(0.02782+0.01247)/0.05;
    elseif (epsilon < 0.30)
        J = 0.02782+(epsilon - 0.25)*(0.1179-0.02782)/0.05;
    elseif (epsilon < 0.4)
        J = 0.1179+(epsilon - 0.30)*(0.4076-0.1179)/0.1;
    elseif (epsilon < 0.5)
        J = 0.4076+(epsilon - 0.40)*(0.7350-0.4076)/0.1;
    elseif (epsilon < 0.6)
        J = 0.7350+(epsilon - 0.50)*(1.0236-0.7350)/0.1;
    elseif (epsilon < 0.7)
        J = 1.0236+(epsilon - 0.60)*(1.2554-1.0236)/0.1;
    elseif (epsilon < 0.8)
        J = 1.2554+(epsilon - 0.70)*(1.436-1.2554)/0.1;
    elseif (epsilon < 0.9)
        J = 1.436+(epsilon - 0.80)*(1.576-1.436)/0.1;   
    elseif (epsilon < 1)
        J = 1.576+(epsilon - 0.90)*(1.686-1.576)/0.1;    
    elseif (epsilon < 1.5)
        J = 1.686+(epsilon - 1.0)*(1.979-1.686)/0.5;
    elseif (epsilon < 2)
        J = 1.979+(epsilon - 1.5)*(2.094-1.979)/0.5;
    elseif (epsilon < 5)
        J = 2.094+(epsilon - 2)*(2.227-2.094)/3;
    else        
        J = 2.255-0.6463*epsilon^(-2);
    end %if, calculating the Saffman force based on the range in which epsilon lies
    Fsaff = -1.615*d^2*J*sqrt(rho*mu*abs(dU))*V_rel_z; %the Saffman lift force due to the shearing flow
    Fx = Fd + Fsaff + Fmag + F_grav; %the total force including the Saffman force
end %if, checking whether the Saffman force should be calculated, which is computationally intensive and has little effect on trajectories