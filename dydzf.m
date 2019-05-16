function dydz=dydzf(z,y,Dia,m_water,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,OverallU,Area,Oxygen_frac)
%% Chemical Properties and Constants
g = 9.81;                          %m/s^2 grav acc
R = 8.314;                         %J/mol.K Gas Constant
dir=1;                             %Direction: 1 for counter current, -1 for cocurrent
uL =V_LiqIn/Area;                  %m/s, Bulk liquid Flow
LiqDen = 902.3;                    %kg/m^3 Bulk Liq Density
cp_water = 4186;                   %J/kg.K Specific Heat Capacity of water
cp_Liq = 1803;                     %J/kg.K Specific Heat Capacity of Liquid
H_Rxn=-85000;                       %J/mol Heat of Reaction,Positive for Exothermic reaction
visL= 0.7109*(10^-3);              %Pa.s Viscosity of liquid phase
visG= 0.02109*(10^-3);             %Pa.s Viscosity of Gas Phase
SurfT= 24.3*0.001;                 %N/m Surface Tension
Ea = 54428.4;                      %J/mol Activaton Energy
Dab = (4.8*(10^-9));               %m^2/s Diffusvity coeff of oxygen
aw = pi()*Dia;                     %m^3/m^2 Specific Heat transfer area
kw = OverallU;%41/WallThickness;       %W/(m.K) Heat transfer number
dk = Dia;                          %m Diameter of heat transfer pipe
kVisL = visL/LiqDen;               %m^2/s Kinematic Viscosity Liquid phase

%% Correlations and Equations 
H = 11044.4*exp((1/y(2))-(1/323.15));         %Pa.m^3/mol Henry Constant
k = 3.83*10^(-3)*exp((Ea/R)*(1/323-1/y(2)));  %m^3/mol.s Kinetic rate constant
EL=2.7/10000*((Dia*100)^1.4)*((y(6)*100)^0.3); %m^2/s
MW = y(7)*32+(1-y(7))*28 ; %MW of gas (g/mol)
GasDen = y(5)*MW/1000;       %kg/m^3
Gas_Holdup =  y(6) /(0.24 + 1.35 * ((y(6) + uL)^0.93)); %uG/(0.24+1.35*(uG+uL)^0.93))
%0.672.*((y(1).*visL./SurfT).^0.578).*(((g.*visL.^4)./(LiqDen.*SurfT.^3)).^-0.131).*((GasDen./LiqDen).^0.062).*((visG./visL)^0.107);
EG = 9.36*(10^-5)*((Dia*100)^1.33)*((y(6)*100/Gas_Holdup)^3.56)/10000; %m^2/s
Kla =(0.6*((kVisL/Dab)^0.5) * (((g*(0.6^2)*LiqDen)/(SurfT))^0.62)*(((g*(0.6^3))/(kVisL^2))^0.31)*(Gas_Holdup^1.1)*Dab)/(0.6^2); %AAA
%0.6.*(Dab.^0.5).*(kVisL.^-0.12).*((SurfT./LiqDen).^-0.62).*(Dia.^0.17).*(g.^0.93).*(Gas_Holdup.^1.1); %Liq mass transfer
Liq_Holdup = 1-Gas_Holdup;     %Liq Holdup (eL)
%k= A.*exp(-Ea./(R.*y(6)));    %Rate constant
%if Liq_Holdup > 1 || Gas_Holdup > 1
%    warning('Liquid Holdup or Gas Holdup >1\nGas Holdup: %s Liquid Holdup:%s',Liq_Holdup,Gas_Holdup);
%end
if y(9) <= 0 || y(11) <= 0
    r = 0;
else
    r = k*y(9)*y(11);            %Rate Law
end
eqO2 = y(7)*y(1)/H;            %Eqilibrium Conc of oxygen in liquid phase yi*P*H
LamdaEff = EL*cp_Liq*LiqDen;   %Effective thermal conductivity or heat dispersion coefficient



%% Odes
%{ 
Notes:
y(1)= P || y(2)= T || y(3)= dT/dz || y(4)= Tw
y(5) = cG || y(6) = Ug || y(7) = ya || y(8) = dy/dz
y(9) = cAL || y(10)= dcAL/dz || y(11)= cBL || y(12)= dcBL/dz 
y(13) = cCL || y(14) = dcCL/dz, A=O2, B=THEAQH2, C=H2O2
%}
%Pressure balance
dydz(1) = (Liq_Holdup*LiqDen+Gas_Holdup*GasDen)*-g; %y(1) = P
%Energy Balance
dydz(2,1) = y(3); %y(2) = T
dydz(3,1)= (1/(Liq_Holdup*LamdaEff))*(kw*aw*(y(2)-y(4))-dir*LiqDen*cp_Liq*uL*y(3)-Liq_Holdup*-H_Rxn*r); %y(3) = dT/dz
dydz(4,1)= (kw*aw*(y(2)-y(4)))/(m_water*cp_water); %y(4)=Tw Jung et al Paper
%Conc. of Gas balance (CG)
dydz(5,1) = ((y(2)*dydz(1))-(y(1)*y(3)))/((y(2)^2)*R); %y(5) = cG
%Overall balance
dydz(6,1)= (1/y(5))*(-y(6)*dydz(5)-Kla*(eqO2-y(9))); %y(6) = ug
dydz(7,1) = y(8) ; %y(7) = ya, fraction of A in gas phase
dydz(8,1) = ((1/(Gas_Holdup*EG))*(y(6)*(y(7)*dydz(5)+y(5)*y(8))+y(7)*y(5)*dydz(6)+Kla*(eqO2-y(9)))-(y(8)*dydz(5)))/y(5); % y(8) = ya',dydz(8) = ya"
%O2 Balance
dydz(9,1) = y(10); %y(9) = cAL
dydz(10,1) = (1/(Liq_Holdup*EL))*(Liq_Holdup*r-dir*uL*y(10)-Kla*(eqO2-y(9))); %y(10) = dcAL/dz
%Theaqh2 Balance  
dydz(11,1)= y(12); %y(11) = cBL
dydz(12,1)= (1/(Liq_Holdup*EL))*(Liq_Holdup*r-dir*uL*y(12)); %y(12) = dcBL/dz
%H2O2 Balance
dydz(13,1)= y(14); %y(13) = h202 conc (cCL)
dydz(14,1)= (1/(Liq_Holdup*EL))*(-Liq_Holdup*r-dir*uL*y(14)); %y(14) = dcCL/dz
end