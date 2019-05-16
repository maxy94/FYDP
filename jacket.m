function [sol,type,ReL,ReT] = jacket(di,As,V_m,V_r,wallT)
format long
%% Params
%{
di  - Inner diameter of reactor, m
As  - Aspect Ratio (L/D)
V_m - Volumetric flow rate of coolant solution, m^3/s
V_r - Volumetric flow rate of reactor bulk liquid, m^3/s
U   - Desired overall heat transfer coefficient, W/m^2.K
wallT - Vessel Wall thickness, m
%}

%% Reactant Properties to calculate (Solution,56,pg227)
U = 35;
k_reactant = 0.1474;               %W/(m.K) Thermal Conductivity of Reactant (50 C)
vis_reactant =  0.7109*(10^-3);    %Pa.s Viscosity of Reactant
den_reactant = 902.3;              %kg/m^3 Reactant Liq Density
cp_reactant = 1803;                %J/kg.K Specific Heat Capacity of Liquid

Re_r = di*(V_r/(pi()/4*di^2))*den_reactant/vis_reactant;
Pr_r = vis_reactant*cp_reactant/k_reactant;
h_inf = (k_reactant/di)*0.023*Re_r^0.8*Pr_r^0.3;
h_L = (1+(di/(di*As))^0.7)*h_inf;

if (1/U-1/h_L) < 0
    error('Overall Heat Transfer Coefficient (U) specified is bigger than calculated reactant side specific heat transfer coefficient.Specified:%g Calculated k of Reactant:%g\n',U,h_L)
end
%% Solver
% Solve with Laminar Conditions first
[solL,ReL] = jacketsolver(di,As,V_m,h_L,U,wallT,true);
% Solve for Turbulent Conditions
[solT,ReT] = jacketsolver(di,As,V_m,h_L,U,wallT,false);

% Filter Solutions based on flow region
if ReL < 2300
    sol = solL;
    type = 'Laminar';
elseif ReL>4000 && ReT > 4000
    sol = solT;
    type = 'Turbulent';
else
    sol = (solL+solT)/2;
    type = 'Transitional';
    warning('Transitional or unknown Regime. Average solution taken.\nLaminar Solution:%g (Re:%g)\nTurbulent Solution:%g (Re:%g)',solL,ReL,solT,ReT);
end
end

function [solution, Rey] = jacketsolver(di,As,V,h_L,U,wallT,isLam) %di = Diameter of vessel, V = Vol flow rate of coolant, As = Aspect
%% Input Properties 
cp = 4186;                         %J/(kg.K) Specific Heat Capacity of water                  
wallK = 50.2;                      %W/(m.K) Thermal Conductivity of Wall, Steel
den_cool = 1000;                   %kg/m^3 Density of Coolant
vis_cool = 0.001;                  %Pa.s Dyn viscosity of Coolant
k_cool = 0.6145;                   %W/(m.K) Thermal Conductivity of Coolant (30 C)
ff = 0.0002;                       %m^2 K/W https://checalc.com/solved/jacketHeat.html
%% Calculations of Dimensionless Numbers
syms d
assume(d,'real')
%d= 0.3;
do = d+di+d;                        %m Outer Diameter 
L = di*As;                          %m Length of Reactor
M = V*den_cool;                     %m^3/s Mass Flow Rate
kw = wallK/wallT;                   %W/(m^2.K) Heat transfer number
if isLam
    eqD = do-di;                    %m, Eqivalent Diameter (Laminar)
else
    eqD = (do.^2-di.^2)/di;        	%m, Eqivalent Diameter (Turbulent)
end
A = (pi()/4)*(do^2-di^2);           %m^2 Cross sectional Area of jacket

Re = eqD*(V/A)*den_cool/vis_cool;   %Reynolds Number
Pr = cp*vis_cool/k_cool;            %Prantl Number
Gr = M*cp/(k_cool*L);               %Grashof Number

%% Equations to Solve
h_v = (k_cool/eqD)*1.02*Re^0.45*Pr^0.33*(eqD/L)^0.4*(do/di)^0.8*Gr^0.05;
eqn = 1/(1/h_v+1/h_L+1/kw+ff) == U;
%Us = 1/(1/h_v+1/kw)
solution = double(solve(eqn,d)); %Convert symbolic vars back to number
clear d

%% Reynolds Number Calculation
do = solution*2+di;
if isLam
    eqD = do-di;                 %m, Eqivalent Diameter (Laminar)
else
    eqD = (do^2-di^2)/di;        %m, Eqivalent Diameter (Turbulent)
end
A = (pi()/4)*(do^2-di^2);
Rey = eqD*(V/A)*den_cool/vis_cool;
end