function res=bcf(ya,yb,Tin,Pin,Dia,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,Area,Oxygen_frac) %ya=initial, z=0, yb = limits at z=z , z=length
%{ 
Notes:
y(1)= P || y(2)= T || y(3)= dT/dz || y(4)= Tw
y(5) = cG || y(6) = Ug || y(7) = ya || y(8) = dy/dz
y(9) = cAL || y(10)= dcAL/dz || y(11)= cBL || y(12)= dcBL/dz 
y(13) = cCL || y(14) = dcCL/dz, A=O2, B=THEAQH2, C=H2O2
%}
uL = V_LiqIn/Area;
eGb = yb(6)/(0.24+1.35*((yb(6)+uL)^0.93));
eLb = 1 - eGb;
ELb = 2.7*((Dia*100)^1.4)*((yb(6)*100)^0.3)/10000;


res(1) = yb(1)-Pin;
res(2,1) = yb(2)-Tin;
res(3,1) = ya(3)-0;
res(4,1) = ya(4)-303;
res(5,1) = ya(5)-F_AirIn/V_O2In;
res(6,1) = ya(6)-V_O2In/Area;
res(7,1) = ya(7)-Oxygen_frac;
res(8,1) = yb(8)-0;  
res(9,1) = F_a0/Area-(yb(9)*uL+eLb*ELb*yb(10));
res(10,1) = ya(10)-0;
res(11,1) = F_b0/Area-(yb(11)*uL+eLb*ELb*yb(12));
res(12,1) = ya(12)-0;
res(13,1) = 0 -(yb(13)*uL+eLb*ELb*yb(14));
res(14,1) = ya(14)-0;
end
