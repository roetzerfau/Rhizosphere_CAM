function  [ N_SVector, C_SVector, N_BVector, C_BVector, MNsolidVector] = ...
calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector, bulkVector, MNsolidVector)


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
    
    bact_idx = find(C_BVector >= parameters.minConC_B);
    for i = 1:numel(bact_idx)
        bact_i = bact_idx(i);
        tspan = [0:tau];
        y0 = [C_SVector(bact_i), C_BVector(bact_i), N_SVector(bact_i), N_BVector(bact_i) ]';
        [t,y] = ode45(@(t,Y) MMKfunction(t,Y,parameters), tspan, y0);
        
        %y_2_4 = [y(:,2),y(:,4)];
        
        %plot(t,y,'-o')
        %legend('C_S', 'C_B', 'N_S', 'N_B')
        %legend('C_B','N_B')

        C_S= y(end,1);
        C_SVector(bact_i)= C_S;
        C_B = y(end,2);
        C_BVector(bact_i) = C_B;
        N_S = y(end,3);
        N_SVector(bact_i) = N_S;
        N_B = y(end,4);
        N_BVector(bact_i) = N_B;
        
    end
    %C_N_B = C_B/N_B;
    %C_N_S = C_S/N_S;
    
    
   % C_BVector = spreadConcentration(g, C_BVector,  bulkVector, parameters.maxConcC_B);
   % N_BVector = C_BVector ./ parameters.C_N_B;
    
end
function dYdt = MMKfunction(t,Y, parameters)
    C_S = Y(1);
    C_B = Y(2);
    N_S = Y(3);
    N_B = Y(4);

    U = (parameters.v_Cliquid *C_S)/(C_S + parameters.K_Cliquid) * C_B;
    %U = 0.1 * C_S;
    %parameters.Resp_GE = 0;
    %parameters.Resp_Maint = 0;
    
    
    R_GE = parameters.Resp_GE * U;
    R_M = parameters.Resp_Maint * C_B; 
    R_O = 0; %overflow resp -> strong imbalance of resources
    BD = 0;
    EX = 0;
     
    
    G = U - R_GE - R_M - R_O - EX;
    CUE = G/U;
    
    
    if(N_S == 0 && C_S > 0)
        C_N_S = Inf;
    elseif(N_S == 0 && C_S == 0)
            C_N_S = 1;
    else 
        C_N_S= C_S/N_S;
    end
       
    if(N_B == 0 && C_B > 0)
        C_N_B = Inf;
    elseif(N_B == 0 && C_B == 0)
            C_N_B = 1;
    else 
        C_N_B= C_B/N_B;
    end
    
    C_N_CR = C_N_B/(1-parameters.Resp_GE);
     
    %Phi = U * (parameters.eta_PAR/C_N_S - 1/C_N_CR);
    if(U == 0) % TODO this case
        Phi = R_M / C_N_B;
    end
    Phi = U * (parameters.eta_PAR/C_N_S - CUE/C_N_B);
    
    % Organic N in excess (Phi > 0) -- Organic C in excess (Phi < 0)
    C_N_IMM = parameters.eta_PAR * C_N_CR; %-> Phi = 0   
      
   % N-limitation occurs:
   C_N_LIM = 0;
   
    %% C overflow hypothesis (CO)
     IMM_max = 0;
%     if (abs(Phi) <= IMM_max)
%         R_O = 0;
%     else
%         R_O = U * (1-parameters.Resp_GE) - C_N_B *(parameters.eta_PAR * U /C_N_S + IMM_max);
%     end
%     
    if (Phi >= IMM_max)
        R_O = 0;
    else
         %R_O = C_N_B * abs(Phi);
         R_O = U * (1-parameters.Resp_GE) - C_N_B *(parameters.eta_PAR * U /C_N_S + IMM_max);
         %Phi = - IMM_max;
         Phi = R_M / C_N_B ;
    end
   
    %R_O = -C_N_B * (parameters.eta_PAR * U/C_N_S - Phi) - R_M +  (1-parameters.Resp_GE) * U;
    C_eq = (1-parameters.Resp_GE) * U - R_O - R_M ;
    N_eq = C_N_B * (parameters.eta_PAR * U/C_N_S - Phi);
    
    %abs(C_eq - N_eq)
    assert(abs(C_eq - N_eq) < 0.0001, "ERROR C_N_B balance ", abs(C_eq - N_eq))
    
    
   
    N_B_dt = parameters.eta_PAR * U/C_N_S - Phi - BD/C_N_B;
    C_B_dt = U - R_GE- R_M - R_O - BD;
    C_S_dt = -U;
    %C_S_dt = 0;
    N_S_dt = -(parameters.eta_PAR * U/C_N_S - Phi - BD/C_N_B);
   % N_S_dt = 0;
    
    %a = [C_N_B, Phi, R_O]
    dYdt = [C_S_dt; % C_S -U *0.01
             C_B_dt; % C_B
            N_S_dt; % N_S  %-U/C_N_B
             N_B_dt]; % N_B         % U/C_N_S   U * CUE/C_N_B
         
end