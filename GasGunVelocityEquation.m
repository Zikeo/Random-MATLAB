function [Pressure] = GasGunVelocityEquation(V,M,L,d,a0,y,E,varagin) %#ok<INUSD>
% Zak Wilde, ASU, 11/14/18, Version 0.1
% This computes the pressure required to achieve the desired velocity from
% a Gas Gun. The types of parameters are related to the sabot, the gun
% barrel  and the Gas used.  M, is the mass of the projectile (grams).
% L and d are the barrel length and inner diameter (ft and in). a0 and y
% are the sound speed and specific heat of the gas used. V is velocity is
% in (m/s). E is efficency of gun (1 being perfect)

if nargin == 0
    V =700; % Desired Velocity  (m/s)
    % p0=2e7;
    a0 = 972; % Speed of sound (m/s) (He) https://www.engineeringtoolbox.com/speed-sound-gases-d_1160.html
    y = 1.667; %Ratio of specific heats (He)https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html

    M=2;
    L=5;% GG Barrel Length ft
    d = 1/2; % GG Inner Diameter (in)

    disp('No inputs: showing example')
end

M=M*1e-3; % Target Mass (kg)
L=L*0.3048;% GG Barrel Length m
S=pi/4*(d*0.0254)^2; %GG Barrel cross section m^2)
GGeq=@(p0,S,L,M,a0,y,v)-p0*S*L/M + 2*a0^2/(y+1)*(1+(((y+1)/(2*a0))*v/E-1)/((1-(y-1)*v/E/(2*a0))^((y+1)/(y-1))));

if length(V) == 1
    Pressure=fzero(@(p0)GGeq(p0,S,L,M,a0,y,V),100)/6894.757; % Divison converts to PSI
    fprintf('Pressure = %4.2f (psi), V = %4.2f (m/s), M = %4.2f (g), L = %4.2f (m), d = %4.2f (in), S = %4.6f (m^2), a0 = %4.2f (m/s), y = %4.2f \n'  ,[Pressure,V,M*1000,L,d,S,a0,y])
elseif length(V) > 1
    Pressure=ones(1,length(V))*NaN;
    for i=1:length(V)
        [Pressure(i)] = fzero(@(p0)GGeq(p0,S,L,M,a0,y,V(i)),100)/6894.757; % Divison converts to PSI
    end
end