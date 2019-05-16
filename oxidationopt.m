function oxidationopt(excel,v1,v2,v3,v4,v5,v6,v7)
close all
format short

%% Reactor Inputs and Parameters
%Fixed Flow rates
F_b0 = 616.83*1000/3600;        %mol/s Molar Flow rate of THEAQH2 Inlet, Liq Phase
F_a0 = 2.8314*1000/3600;        %mol/s Molar Flow rate of Oxygen Inlet, Liq Phase
F_aOut = 9.332*1000/3600;       %mol/s Molar Flow rate of Oxygen in Reactor Product
F_bOut = 33.3*1000/3600;        %mol/s Molar Flow rate of Oxygen in Reactor Product
V_LiqIn = 968.4/3600;           %m^3/s Total Liquid Volumetric flow rate
V_LiqOut= 984.4/3600;           %m^3/s Total Liquid Volumetric flow rate
Oxygen_frac = 0.936;              %Mole Fraction of Oxygen in inlet to reactor

%Reactor Variables (Input with commas)
Var{1} = [50]+273.15;           %Tin, K Temp of inlet (Enter Values in Celsius)
Var{2} = [4]*100000;            %Pin, N/m^2(Pa) Pressure of inlet (Enter Values in bar)
Var{3} = [50];                	%W/m^2.K Overall Heat transfer coeff thru wall
Var{4} = [3];                    %Diameter,m
Var{5} = [4];                     %Aspect Ratio, L/D
Var{6} = [1.7];                   %m^3/s Volumetric flow rate of oxygen feed, V_O2In
Var{7} = [0.1];           %m^3/s V_water,Total Cooling Volumetric Flow

% Excel Exports
filename='oxidation.xlsx';
headers = {'Inlet Tempreture/K', 'Inlet Pressure/Pa','Wall Thickness/m','Diameter/m',... 
    'Aspect Ratio', 'O2 Inlet Vol. Flow Rate/ m^-3/s', 'Cooling Water Vol. Flow Rate/ m^3/s','Conversion/%','maxT/K','Gas Flow Vel/ ms^-1'};
  
%% Wrapper Code
[row col] = size(Var);

%Check if the user specifies any input parameters for variables and use
%it if it exists
varCount = 0;
for v=1:col
    if exist(strcat('v',num2str(v)),'var')
        varCount = varCount+1;
    end
end

%Excel Exports check
excelexport = false;
if exist('excel','var') && (ischar(excel) || excel)
    excelexport = true;
    sheetCount = 1;
    if islogical(excel)
        excel = 'Data';
    end
    fprintf('Excel Data Export Active\n')
end

%Check if optional input params are used
if varCount > 0 && varCount ~= col
    error('Not Enough Input Parameters Entered! 8 Parameters required. First Parameter(Excel Export) must be false/true, followed by 7 params.')
elseif varCount==col
    fprintf('Single Variable Problem Input\n')
    Var{1} = [v1]+273.15;
    Var{2} = [v2]*100000;
    Var{3} = [v3];
    Var{4} = [v4];
    Var{5} = [v5];
    Var{6} = [v6];
    Var{7} = [v7];
    [row col] = size(Var);
end
       
%Output names and units for Optimizer Report
PropName = {'Inlet Tempreture:', 'Inlet Pressure:','Overall U:','Diameter:',... 
    'Aspect Ratio:', 'O2 Inlet Vol. Flow Rate:', 'Cooling Water Vol. Flow Rate:','Max. Conversion:', 'Reactor Volume:'};

PropName(3,:) = {'K','Pa','m','m','','m^3','m^3','%','m^3'}; %Units for corresponding Values

PlotTitle = {'Pressure', 'Tempreture', 'dT/dz', 'Wall Tempreture Tw', 'Gas Concentration cG', 'Gas Superficial Velocity uG',...
            'Mole Fraction of O2 ya', 'dya/dz',  'Conc. of O2 cAL','dcAL/dz', 'Conc. of THEAQH2 cBL', 'dcBL/dz','Conc of H2O2 cCL', 'ccCL/dz'};
PlotUnits = {'Pa', 'K', 'K/m','K','mol/m^3','m/s','val','m^-1','mol/m^3','val','mol/m^3','val','mol/m^3','val'};

%Variables Initialization
j = 1;
k = 0;
i_conv = 1;
isSingleVar = true;
jacobianErrors=0;
VarNo = 5;
vol = [];

%Multi Variable Problem Checker
for m = 1:col
    if numel(Var{m}) > 1
        k = k+1;
        if k == 1
            VarNo=m;
            list = Var{m};
        elseif k > 1
            isSingleVar = false;
        end
    end 
end

%While there are 7 for loops, add more element to only those that are most
%sensitive!!
for a = 1:numel(Var{1})
    Tin=Var{1}(a); 
    for b = 1:numel(Var{2})
        Pin=Var{2}(b);
        for c = 1:numel(Var{3})
            OverallU = Var{3}(c);
            for d = 1:numel(Var{4})
                Dia=Var{4}(d);
                for e = 1:numel(Var{5})
                    Aspect = Var{5}(e);
                    for f = 1:numel(Var{6})
                         V_O2In = Var{6}(f);
                        for g = 1:numel(Var{7})
                            V_water = Var{7}(g);
                            
                            Length = Dia*Aspect;       %L/D 
                            Area = pi()*Dia/4;         %m^2 Area of reactor
                            m_water = V_water*1000;    %kg/s Mass flow rate of cooling water
                            Area = pi()*Dia/4;         %m^2 Area of reactor
                            F_AirIn = ((Pin*V_O2In)/(8.314*Tin))/Oxygen_frac; %mol/s
                            try
                                solve %Nested Solver function
                            
                            catch err %We catch any jacobian errors and warn instead
                                if(strcmp(err.identifier,'MATLAB:bvp4c:SingJac'))
                                    conversion(i_conv) = 0;
                                    i_conv = i_conv+1;
                                    jacobianErrors = jacobianErrors+1;
                                    alpha = [a,b,c,d,e,f,g];
                                    for h = 1:col
                                        errorList{h,jacobianErrors} = Var{h}(alpha(h));
                                    end
                                    fprintf('Jacobian Error Encountered: Skipping Iteration\nVariables used\n')
                                    fprintf('%g ',errorList{:,jacobianErrors})
                                    fprintf('\n')
                                else
                                    rethrow(err) %Throw error for those of non-jacobian nature
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Solver per case
    function solve
        sol = solver(Dia,m_water,Length,Tin,Pin,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,OverallU,Area,Oxygen_frac); %Loop over this wrapper function to optimize script later on
        z = sol.x; y = sol.y; 

        %molA(j) = F_AirIn*0.936+F_a0-y(1,end)*y(6,end)*Area/(8.314*y(2,end))*0.936-y(9,1)*V_LiqOut;
        %molB(j) = F_b0-y(11,1)*V_LiqOut;

        if isSingleVar %Plot only if a single variable problem is concerned
            for i = 1:14
                figure(1)
                hold on
                subplot(5,3,i)
                plot(z,y(i,:),'linewidth', 1.25)
                title(PlotTitle{i}),
                xlabel('z/m'),
                ylabel(PlotUnits{i})
            end
            if exist('list','var')
                legend(sprintfc('%g',list));
            end
        end
        
        
        
        %Set Conversion to 0 if it becomes >100% (infeasible)
        currentConv = (F_b0/V_LiqIn-y(11,1))*100./(F_b0/V_LiqIn);
        if not(isreal(currentConv)) || (currentConv > 100)
            currentConv = 0;
        end
        
        maximumtemp = max(y(2,:));
        outputUG = y(6,end);
        pressuredrop = y(1,1)- y(1,end);
        outletLiqTemp = y(2,1);
        
        if isSingleVar
            %molo2In = F_AirIn*V_O2In*0.932+F_a0*V_LiqIn
            %molo2Out = y(9,1)*V_LiqOut + y(7,end)*(y(1,end)*
            conversion(i_conv) = currentConv;
            vol(i_conv) = (Var{4}(d)*pi()/4)*Var{4}(d)*Var{5}(e);
            i_conv = i_conv+1;
            maxT(i_conv) = maximumtemp;
            pDrop(i_conv) = pressuredrop;
            outTemp(i_conv) = outletLiqTemp;
        else
            conversion(a,b,c,d,e,f,g) = currentConv;
            maxT(a,b,c,d,e,f,g) = maximumtemp;
        end
        
        % Sanitize Data
        if excelexport && currentConv >= 94.6  && currentConv < 100 && maximumtemp<333.15 %60degC
            if exist('excelMatrix','var')
                excelMatrix(end+1,:)=[Tin,Pin,OverallU,Dia,Aspect,V_O2In,V_water,currentConv,maximumtemp,outputUG];
            else
                excelMatrix(1,:)=[Tin,Pin,OverallU,Dia,Aspect,V_O2In,V_water,currentConv,maximumtemp,outputUG];
            end
        end
        
        
        
   %{
        cap = xxxxxxxxx; %Insert cost formulae here!!!!
        com= xxxxxx;
        land = xxxx;
        r = 2.5/100;
        
        cost = land;
        for z=1:12
            if i == 1 || i == 2
                cost = cost + cap/((1+r)^z);
            else
                cost = cost + com/((1+r)^z);
            end
        end     
        %We increase the cost for any less than desirable conversion to
        %filter them out from min search function
        
        if currentConv < 94.6 
            cost = cost*1e+20;
        end
        cost(a,b,c,d,e,f,g) = cost
        %}
        
        j=j+1;
    end
%% Data Outputs handling
    %Only Plot conversion graph if it is Single Var Problem
if isSingleVar
    vol
    conversion
    maxT
    pDrop
    outTemp
    hold off
    figure(2)
    plot(Var{VarNo},conversion,'linewidth', 1.25)
    title('Conversion Vs Single Variable'),
    xlabel('Variable'),
    ylabel('Conversion')
%Use maxmat/minmat (MATLAB repo code) to get optimized cost/conv.
else
    [maximum,index] = maxmat(conversion);
    %[max,index] = minmat(cost);
    %Lowest_Cost_Conversion = conversion(index)
    PropI = ones(1,col);
    PropI(1:numel(index)) = index;
    %{
    1. Inlet Temp
    2. Inlet Pressure
    3. Wall thickness
    4. Diameter
    5. Aspect Ratio
    6. Inlet Vol. Flow rate
    7. Cooling Water Flow rate
    %}
    for o=1:numel(PropI)
        Prop(o) = Var{o}(PropI(o));
    end
    Opt_Volume = Prop(4)*Prop(5)*(pi()*Prop(4)^2)/4;
    Prop = [Prop,maximum,Opt_Volume]; %Add more Calculated Values here as needed
    for o=1:numel(PropI)
        Prop(o) = Var{o}(PropI(o));
    end
    %Can't touch this
    PropName(2,:) = num2cell(Prop);
    
% Print Data Outputs
    fprintf('Optimal Parameters:\n\n')
    fprintf('%s %g %s\n',PropName{:})
end
    fprintf('Jacobian Errors Encountered: %g\nJacobian Error List:\n T    P   Dw  D AR V_O2 V_C\n',jacobianErrors)
    if jacobianErrors > 0
        fprintf('%g %g %g %g %g %g %g\n',errorList{:})
    end
    
end



%% BVP4C Solver
function sol = solver(Dia,m_water,Length,Tin,Pin,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,OverallU,Area,Oxygen_frac)
%{ 
Notes:
y(1)= P || y(2)= T || y(3)= dT/dz || y(4)= Tw
y(5) = cG || y(6) = Ug || y(7) = ya || y(8) = dy/dz
y(9) = cAL || y(10)= dcAL/dz || y(11)= cBL || y(12)= dcBL/dz
y(13) = cCL || y(14) = dcCL/dz
%}
%% Code
fh1 = @(z,y) dydzf(z,y,Dia,m_water,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,OverallU,Area,Oxygen_frac);
fh2 = @(ya, yb) bcf(ya,yb,Tin,Pin,Dia,F_b0,V_LiqIn,V_O2In,F_a0,F_AirIn,F_aOut,F_bOut,V_LiqOut,Area,Oxygen_frac);
z = 0:0.01:Length;

%Initial Guess at z=0
initguess=[500000 323 0 303 F_AirIn/V_O2In V_O2In/Area Oxygen_frac 0 F_aOut/V_LiqOut 0.5 F_bOut/V_LiqOut 10 0 0];
%initguess=[0 0 0 0 0 0 0 0 0 0 2 0 0 0];

solinit = bvpinit(z, initguess);
opts = bvpset('NMax',10000,'RelTol',1e-5); %Set Iterative step to 10000 to prevent tolerance issues
sol = bvp4c(fh1, fh2, solinit,opts);
end