%%%%%%%%%%%%%%%% Rabbit AVN Model Structure & Parameters %%%%%%%%%%%%%%%%%
% 1D multi-functional model of rabbit AV node with dual pathways
% v 1.06 (c) 2022
%     Maxim Ryzhii, University of Aizu, Japan
%     Elena Ryzhii, Fukushima Medical University, Japan
%
% Code for the paper "A compact multi-functional model of the rabbit 
% atrioventricular node with dual pathways", 
% Frontiers in Physiology, 14 (2023). DOI: 10.3389/fphys.2023.1126648
%
% Tested with MATLAB R2022b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = avn_fun(t,y,Mode,S,ap,tmp_distrib) 
global tN     % pulse counter

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbitAVN''\n');
   return
end

y1 = reshape(y,[S.iHB2,4]);
Cx = zeros(S.iHB2,2); 
Cy = zeros(S.iHB2,1);
Cyp1 = zeros(S.iHB2,1);
Cyp2 = zeros(S.iHB2,1);
dydt = zeros(S.iHB2,4);

tdur = ap.pulse_width;    % pulse duration (s)
tau = 0.00002;            % pulse time constant
pulse = 0;

% Set S1-S2 stimulus
if Mode >= 10
    r = tmp_distrib(tN);
    if t < r+tdur*10 
            tmod = mod(t-r/2,r)-r/2;       % modulo time
            yup = 1./(1+exp(-tmod/tau));   % this part goes up, from 0 to 1, as tmod crosses 0
            ydown = 1./(1+exp((tmod-tdur)/tau)); 
            pulse = ap.pulse_amp*yup*ydown;
    elseif tN < length(tmp_distrib)
        tN = tN+1;
    end
end

%%%% "Horizontal" coupling matrix update (with asymmetry)
% First row
Cx(    S.iSN,1) = ap.dx(S.iSN,1).*(y1(S.iSN+1,1) -y1(S.iSN,1).*ap.adx(S.iSN,1)); 

Cx(S.iSN+1:S.iHB2-1,1) = ap.dx(S.iSN:S.iHB2-2,1).*( y1( S.iSN:S.iHB2-2,1) - ...    % Left
                         y1(S.iSN+1:S.iHB2-1,1).*ap.adx(S.iSN:S.iHB2-2,1)) + ...   % Left
                         ap.dx(S.iSN+1:S.iHB2-1,1).*(-y1(S.iSN+1:S.iHB2-1,1) + ... % Right
                         y1(S.iSN+2:S.iHB2,  1).*ap.adx(S.iSN+1:S.iHB2-1,1));      % Right
Cx(    S.iHB2,1) = 1.0.*ap.dx(S.iHB2-1,1).*(y1(S.iHB2-1,1) - ...
                         y1(S.iHB2  ,1).*ap.adx(S.iHB2-1,1));  

% Second row
Cx(S.iSN+1:S.iHB2-1,2) = ap.dx(S.iSN:S.iHB2-2,2).*( y1(S.iSN:S.iHB2-2,3) - ...     % Left
                         y1(S.iSN+1:S.iHB2-1,3).*ap.adx(S.iSN:S.iHB2-2,2)) + ...   % Left
                         ap.dx(S.iSN+1:S.iHB2-1,2).*(-y1(S.iSN+1:S.iHB2-1,3) + ... % Right
                         y1(S.iSN+2:S.iHB2,  3).*ap.adx(S.iSN+1:S.iHB2-1,2));      % Right

%%%%  "Vertical" coupling, AM3<->SP1 and SP10<->PB (with asymmetry)
Cy(S.iAM3) = ap.dy(S.iAM3).*(y1(S.iAM3,1) - y1(S.iSlow1,3).*ap.ady(S.iSlow1));  
Cy(S.iPB)  = ap.dy(S.iPB ).*(-y1(S.iSlow2,3) + y1(S.iPB, 1).*ap.ady(S.iPB) );  

%%%%  
switch Mode   % Pulses 
     case {20,30} 
        Cyp2(S.iAM3-1) = pulse;        % Applied to the bottom row, AM*
        Cy(S.iAM3-1)  = ap.dy(S.iAM3-1).*(y1(S.iAM3-1, 1)- y1(S.iAM3-1,3));  % "Vertical" coupling AM*<->AM2
     case {40,41,42,43,44,45}
        Cyp1(S.imp_Antero) = pulse;   % Applied to the top row
     case {50,51,52,53,54,55}
        Cyp1(S.imp_Retro)  = pulse;   % Applied to the top row
end

% Calculation of derivatives
% First raw (SN-PS-AM-FP-PB-HB) 
dydt(:,1) = ap.Tc(:,1).*(ap.k(:,1).*y1(:,1).*(y1(:,1)-ap.a1(:,1)).*(1.0-y1(:,1))  -y1(:,1).*y1(:,2)) +Cx(:,1) -Cy(:) +Cyp1(:);
dydt(:,2) = ap.Tc(:,1).*(ap.eps0(:,1)+ap.mu1(:,1).*y1(:,2)./(ap.mu2(:,1)+y1(:,1))).*(-y1(:,2)-ap.k(:,1).*y1(:,1).*(y1(:,1)-ap.a2(:,1)-1.0));
% Second raw (AM*, SP1-SP10)
dydt(:,3) = ap.Tc(:,2).*(ap.k(:,2).*y1(:,3).*(y1(:,3)-ap.a1(:,2)).*(1.0-y1(:,3))  -y1(:,3).*y1(:,4)) +Cx(:,2) +Cy(:) +Cyp2(:);
dydt(:,4) = ap.Tc(:,2).*(ap.eps0(:,2)+ap.mu1(:,2).*y1(:,4)./(ap.mu2(:,2)+y1(:,3))).*(-y1(:,4)-ap.k(:,2).*y1(:,3).*(y1(:,3)-ap.a2(:,2)-1.0));

dydt = reshape(dydt,[S.iHB2*4,1]);

end % function end





