function  [ N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, CO2Vector,CO2Vector_over, C_EXTVector] = ...
calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector,C_MNVector, N_MNVector, CO2Vector,CO2Vector_over)

    previousConcentration_C =sum(C_MNVector+ C_SVector + C_BVector + CO2Vector);
    previousConcentration_N =sum(N_MNVector+ N_SVector + N_BVector);
    %previousConcentration_N =sum(C_MNVector+ C_SVector + C_BVector);
    % K_Cliquid: Michaelis-Menten-Konstante/ half saturation constant
    % v_Cliquid: chemical composition of the substrate/ Maximal uptake rate
    tau = parameters.tau_ode;
    %parameters.eta_PAR=1 available organic N is directly assimilated prior to mineralization
   
    
    value2vector = ones(g.numT, 1);
    
%     C_N_S= C_SVector./N_SVector;
%     C_N_B= C_BVector./N_BVector;
%     
% 
%     uptake_rate = bacteriaVector * parameters.v_Cliquid;
%     U = uptake_rate .* ((C_SVector)./(C_SVector + value2vector .* parameters.K_Cliquid)) .* C_BVector;
%     R_GE = parameters.Resp_GE * U;
%     R_M = parameters.Resp_Maint * C_BVector; 
%     R_O = 0; %overflow resp -> strong imbalance of resources
%     BD = 0; 
%     
%     G = U - R_GE - R_M - R_O;
%     CUE = G./U;
%     
%     U_Const = (uptake_rate .* C_SVector)./(C_SVector + value2vector .* parameters.K_Cliquid) ;
%     exponent = ((1-parameters.Resp_GE) * U_Const - value2vector .* parameters.Resp_Maint);
%     changeFract = exp(exponent.* tau);
%     C_BVector_ = changeFract .* C_BVector;
%     
%     a = 1.4 * exp(-0.25);
    C_BVector_old = C_BVector;
    N_BVector_old = N_BVector;

    C_SVector_old = C_SVector;
    N_SVector_old = N_SVector;
    C_EXTVector =  zeros(g.numT, 1);
    bact_idx = find(C_BVector >= parameters.minConC_B);
   % bact_idx = find(C_BVector > 0);
    %check ob bact_idx ist nachbar mit solid, Dann wird EPS produziert

    for i = 1:numel(bact_idx)
        bact_i = bact_idx(i);
         % if(C_SVector(bact_i)< 0)
         %    C_SVector(bact_i)
         % end
        C_EXTVector(bact_i) = parameters.Resp_GE * (parameters.v_Cliquid *C_SVector(bact_i))/(C_SVector(bact_i) + parameters.K_Cliquid) *  C_BVector(bact_i);
        % if( C_SVector(bact_i) < 0 )
        %     fprintf("ahh %f %f \n", C_SVector(bact_i), bact_i)
        % end
        tspan = [0:tau];
        %if(C_SVector(bact_i)< 0)
            % C_SVector(bact_i)
       % end
        y0 = [C_SVector(bact_i), C_BVector(bact_i), N_SVector(bact_i), N_BVector(bact_i), C_MNVector(bact_i), N_MNVector(bact_i), CO2Vector(bact_i), CO2Vector_over(bact_i)]';
        %if(any(y0 < 0))
        %    fprintf("y0 < 0")
        %end
        %options = odeset(RelTol=1e-8,AbsTol=1e-10);
        options = odeset(RelTol=1e-8,AbsTol=1e-9);
        opts = odeset( options) ;%, 'NonNegative',1:7
        [t,y] = ode45(@(t,Y) MMKfunction(t,Y, parameters), tspan, y0, opts);%ode89
        
        %y_2_4 = [y(:,2),y(:,4)];
        
        %plot(t,y,'-o')
        %legend('C_S', 'C_B', 'N_S', 'N_B')
        %legend('C_B','N_B')
        % if(y(end,1) < 0)
        % fprintf("verlgeih C_S %f %f \n",C_SVector(bact_i),y(end,1))
        % end
        % if(y(end,3) < 0)
        % fprintf("verlgeih N_S %f %f \n",N_SVector(bact_i),y(end,3))
        % end
        %fprintf("verlgeih N_S %f %f \n",N_SVector(bact_i),y(end,3))
        
        C_S= max(y(end,1),0);
        C_SVector(bact_i)= C_S;
        C_B = max(y(end,2),0);
        C_BVector(bact_i) = C_B;
        N_S = max(y(end,3),0);
        N_SVector(bact_i) = N_S;
        N_B = max(y(end,4),0);
        N_BVector(bact_i) = N_B;

        C_MN =max(y(end,5),0);
        C_MNVector(bact_i) = C_MN;
        N_MN = max(y(end,6),0);
        N_MNVector(bact_i) = N_MN;
        
        CO2 = y(end,7);
        CO2Vector(bact_i) = CO2;

        CO2_over = y(end,8);
        CO2Vector_over(bact_i) = CO2_over;
    end


    
    



    currentConcentration_C = sum(C_MNVector+ C_SVector + C_BVector + CO2Vector);
    %currentConcentration_C = sum(C_MNVector+ C_SVector + C_BVector);
   %%if(abs(previousConcentration_N - currentConcentration_C) < 0)
   if(abs(previousConcentration_C - currentConcentration_C) > 10^-10)
             %abs(previousConcentration_C - currentConcentration_C)
              fprintf('Falsch MB C %f', abs(previousConcentration_C - currentConcentration_C))
   end

    currentConcentration_N = sum(N_MNVector+ N_SVector + N_BVector);
   if(abs(previousConcentration_N - currentConcentration_N) > 10^-10)
             abs(previousConcentration_N - currentConcentration_N)
            fprintf('Falsch MB N %f', abs(previousConcentration_N - currentConcentration_N))
    end


%difference_MB = abs(previousConcentration_N - currentConcentration_C)
    %C_N_B = C_B/N_B;
    %C_N_S = C_S/N_S;
    
    
   % C_BVector = spreadConcentration(g, C_BVector,  bulkVector, parameters.maxConcC_B);
   % N_BVector = C_BVector ./ parameters.C_N_B;
    
end
function dYdt = MMKfunction(t,Y, parameters)
     persistent last_t
    dt = 1;
    if  ~isempty(last_t)
        dt = t - last_t;  % Compute ΔT
        %fprintf('Current timestep: %.6f\n', dt);
    end
    last_t = t;

    % C_S = abs(Y(1));
    % C_B = abs(Y(2));
    % N_S = abs(Y(3));
    % N_B = abs(Y(4));
    % C_MN = abs(Y(5));
    % N_MN = abs(Y(6));
    % CO2 = abs(Y(7));
%t
    C_S = Y(1);
    C_B = Y(2);
    N_S = Y(3);
    N_B = Y(4);
    C_MN = Y(5);
    N_MN = Y(6);
    CO2 = Y(7);
    CO2_over = Y(8);



    if(N_B == 0 && C_B > 0)
        C_N_B = Inf;
    elseif(N_B == 0 && C_B == 0)
            C_N_B = 1;
    else 
        C_N_B= C_B/N_B;
    end
  
    assert((C_N_B - 10.0) < 0.000000001, "ERROR C_N_B ", C_N_B)
    %C_S = max(C_S, 0);
    %N_S = max(N_S, 0);

    %t
    %C_S
    if(N_S < 0)

       % N_B = N_B + abs(N_S);
       % C_B = C_B + abs(N_S) * 10;
       % C_S = C_S - abs(N_S) * 10;
       N_S =0;
      % fprintf("N_S kleiner 0")
    end
    if(C_S < 0)

       % C_B = C_B + abs(C_S);
       % N_B = N_B + abs(C_S)/10;
       % N_S = N_S - abs(C_S)/ 10;
       C_S =0;
      % fprintf("C_S kleiner 0")
    end
    if(C_B < 0)
     %   fprintf("C_B kleiner 0")
        C_B = 0;
    end

 if(C_S < 0|| N_S <0 )
        C_S_dt = 0;
        C_B_dt = 0;
        N_S_dt = 0;
        N_B_dt = 0;
        C_MN_dt = 0;
        N_MN_dt = 0;
        CO2_dt = 0;
        CO2_over_dt =0;
%  if(C_S < 0 )       
% fprintf("C_S negativ")
%  else
%      fprintf("N_S negativ")
% end
        dYdt = [C_S_dt;
             C_B_dt;
            N_S_dt; 
             N_B_dt;
             C_MN_dt;
             N_MN_dt;
             CO2_dt;
             CO2_over_dt];
        return
   end

    
   

    U = parameters.v_Cliquid *(C_S/(C_S + parameters.K_Cliquid)) * C_B;
    if(U < 0)
        aaa
    end
    if(C_S/(C_S + parameters.K_Cliquid)  < 10^-20)
        U = 0;
    end
    % if(U < 0)
       % U
    %end
    %U = 0.1 * C_S;
    %parameters.Resp_GE = 0;
    %parameters.Resp_Maint = 0;
    
    
    R_GE = parameters.Resp_GE * U;
    R_M = parameters.Resp_Maint * C_B; 
    R_O = 0; %overflow resp -> strong imbalance of resources
    BD = parameters.BD * C_B;
    EXT = 0;
     
 
    if(N_S == 0 && C_S > 0)
        C_N_S = Inf;
    elseif(C_S == 0)
            C_N_S = 1;
    else 
        C_N_S= C_S/N_S;
    end
       

    

    if(N_MN == 0 && C_MN > 0)
        C_N_MN = Inf;
    elseif(N_MN == 0 && C_MN == 0)
            C_N_MN = 1;
    else 
        C_N_MN= C_MN/N_MN;
    end
     


   C_B_dt = U - R_GE- R_M - R_O;
    %% C overflow hypothesis (CO)
  if(C_B_dt/C_N_B < N_S)
      R_O = 0;
  else
      R_O = C_B_dt - N_S *C_N_B;
      %fprintf("R_O")
  end
    C_eq = (1-parameters.Resp_GE) * U - R_O - R_M ;
    N_eq = ((1-parameters.Resp_GE) * U - R_O - R_M)/C_N_B;
    
    
    %abs(C_eq - N_eq)
  
    % assert(abs(C_eq - N_eq) < 0.000001, N_S =0;"ERROR C_N_B balance ", abs(C_eq - N_eq))
    
    %U
    %U - R_GE- R_M -R_O
    %BD 
    %R_O
    %R_GE
    %R_M

    C_B_dt = U - R_GE- R_M - R_O - BD- EXT;
    N_B_dt = C_B_dt/C_N_B;
   
    
    
    N_S_dt = -(U - R_GE- R_M - R_O)/C_N_B;
   
    C_S_dt = -U;

    % if(N_S_dt > 0)
    %  fprintf("asdf")
    % end
    % if(C_S + C_S_dt < 0)
    %    fprintf("check C")
    %    C_S_dt =-C_S;
    % end
    %  if(N_S + N_S_dt < 0)
    %   fprintf("check N")
    %   N_S_dt =-N_S;
    % end

    C_MN_dt = BD;%*parameters.C_N_NM/C_N_B;
    N_MN_dt = BD/C_N_B;


    CO2_dt = R_GE+ R_M + R_O;
    CO2_over_dt = R_O;
    equalC = C_S_dt + C_B_dt + C_MN_dt + CO2_dt;
    if(abs(equalC)> 10^-20)
    error( "ERROR C_N_B C balance equal %f \n", equalC);
    end
    equalN = N_S_dt + N_B_dt + N_MN_dt;
    if(abs(equalN)> 10^-20)
    error( "ERROR C_N_B N balance equal %f \n", equalN);
    end

   %  if(abs(C_eq - N_eq) > 10^-10)
   %    fprintf("--------------------------------")
   %      error("ERROR C_N_B balance %f \n", abs(C_eq - N_eq))
   %      C_S_dt = 0;
   %      C_B_dt = 0;
   %      N_S_dt = 0;
   %      N_B_dt = 0;
   %      C_MN_dt = 0;
   %      N_MN_dt = 0;
   %      CO2_dt = 0;
   %      CO2_over_dt = 0;
   % end

    dYdt = [C_S_dt;
             C_B_dt;
            N_S_dt; 
             N_B_dt;
             C_MN_dt;
             N_MN_dt;
             CO2_dt;
             CO2_over_dt];
         
end