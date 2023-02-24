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
function [Str,ap]= avn_data1()

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbitAVN''\n');
   return
end

% Structure: Indexes of model cells
Str.iSN        = 1;  
Str.iAM3       = 7;          % Last AM cell
Str.iSlow1     = 7;          % First cell of SP
Str.ipSlow     = 13;         % First SP pacemaker
Str.iSlow2     = 16;         % Last cell of SP /last SP pacemaker
Str.iPB        = 16;         % PB
Str.iHB2       = 22;         % Last cell of HB
Str.pointA     = Str.iAM3;   % Measurement point A
Str.pointB     = Str.iPB+1;  % Measurement point B
Str.imp_Antero = Str.iAM3-1; % Impulse application point for Anterograde
Str.imp_Retro  = Str.iPB+1;  % Impulse application point for Retrograde
Str.iAbl       = 11;         % Ablation place (cell)

TC0 = 1000.0/4.0;            % Time constant for physical time
TC  = zeros(Str.iHB2,2);
TC(:,1) = TC0; 
TC(Str.iAM3-2:Str.iPB,2) = TC0;

a_base = 0.13;
a1 = a_base.*ones(Str.iHB2,2);

Str.a1_normal = -0.0694;   % Normal sinus ~166 bpm
Str.a1_max    = -0.1556;   % Max sinus bpm ~350 bpm

a1(Str.iSN,1)      = Str.a1_normal; 
a1(Str.iSN+1,1)    = 0.08;
a1(Str.iPB-1,1)    = -0.0368; % Pacemaker FP
a1(Str.iPB,1)      = -0.0368; % Pacemaker PB
a1(Str.ipSlow,2)   = -0.0369; % Pacemakers SP
a1(Str.ipSlow+1,2) = -0.0369; % 
a1(Str.ipSlow+2,2) = -0.0368; %
a1(Str.ipSlow+3,2) = -0.0368; %

a2 = a_base.*ones(Str.iHB2,2);
a2(Str.iSN,1)      = 0.20; % SN
a2(Str.iSN+1,1)    = 0.16; % PSN1
a2(Str.iPB-1,1)    = 0.28; % Pacemaker FP
a2(Str.iPB,1)      = 0.28; % Pacemaker PB
a2(Str.ipSlow,2)   = 0.28; % Pacemakers SP
a2(Str.ipSlow+1,2) = 0.28; 
a2(Str.ipSlow+2,2) = 0.28; 
a2(Str.ipSlow+3,2) = 0.28; 

mu1 = [0.160 0.180 0.200 0.200 0.220 0.220 0.220 ...        % SN--PS--AM3
       0.150 0.150 0.150 0.150 0.150 0.140 0.140 0.065 ...  % FP
       0.060 ...                                            % PB
       0.020 0.020 0.020 0.020 0.020 0.020;                 % HB
       0     0     0     0     0     0.220  ...             % AM* 
       0.230 0.220 0.215 0.205 0.190 0.175 ...              % SP
       0.070 0.070 0.060 0.060 ...                          % SP pacemakers
       0     0     0     0     0     0];
       
mu2 = [4.20 4.20 3.30 2.50 1.60 1.60 1.60 ...        % SN--PS--AM3
       1.70 1.80 1.85 1.95 2.05 2.15 2.25 4.60 ...   % FP
       4.60 ...                                      % PB
       0.35 0.10 0.10 0.10 0.10 0.10;                % HB
       0    0    0    0    0    1.60 ...             % AM*
       1.80 1.80 1.80 1.80 1.80 1.80 ...             % SP
       2.50 3.30 3.70 3.90 ...                       % SP pacemakers
       0    0    0    0    0    0 ];
    
k = [ 9  9 10 14 18 18 18 ...          % SN--PS1--AM3
      18 18 18 18 18 15 15  8  ...      % FP
       8 ...                            % PB
      10 10 10 10 10 10;  ...           % HB
       0  0  0  0  0 18 ...             % AM*
      11 11 11 11 11 11 ...             % SP
       8  8  8  8 ...                   % SP pacemakers
       0  0  0  0  0  0];

eps0 = [ 0.040 0.036 0.025 0.025 0.018 0.018 0.018 ...       % SN--PS--AM3
         0.011 0.011 0.011 0.011 0.011 0.011 0.011 0.042 ... % FP
         0.042 ...                                           % PB
         0.042 0.050 0.050 0.050 0.050 0.050;                % HB
         0     0     0     0     0     0.018 ...             % AM*
         0.060 0.060 0.060 0.060 0.060 0.060 ...             % SP 
         0.044 0.044 0.042 0.042 ...                         % SP pacemakrs
         0     0     0     0     0     0 ];

% Coupling x 
dx = [ 41 45 47 56 68 95 ...             % SN--PS--AM3 
       95 84 70 64 62 60 59 59 55 ...    % FP--PB
       80 50 50 48 48 48;                % HB
        0  0  0  0  0  0 ...             % AM*
       59 57 57 56 56 53 32 32 32 ...    % SP
        0  0  0  0  0  0];

% Coupling y 
dy = [ 0  0  0  0  0  0 69 ...         % AM2<->AM*, AM3<->SP1
       0  0  0  0  0  0  0  0  42 ...  % SP10<->PB
       0  0  0  0  0];

% Coupling asymmetry x
adx = [ 1    1    1    1    0.75 0.53 ...                  % SN--PS--AM3 
        0.62 0.75 0.80 0.85 0.95 1.00 1.25 1.40 1.50 ...   % FP--PB
        0.65 0.70 1    1    1    1;                        % HB
        0    0    0    0    0    1 ...                     % AM*
        1    1    1    1    1    1    1    1    1    1 ... % SP
        0    0    0    0    0 ];

% Coupling asymmetry y
ady = [ 0   0   0   0   0   1   1 ...           % AM2<->AM*, AM3<->SP1
        0   0   0   0   0   0   0   0   1 ...   % SP10<->PB
        0   0   0   0   0   0];

drnd = [1;  0.075; 0.125];  % Parameters for rnd: Pattern number, range

ap.pulse_width = 0.001;  % [s]
ap.pulse_amp =  280.0; 
ap.Tc  = TC;
ap.a2   = a2; 
ap.a1   = a1;
ap.mu1 = mu1';
ap.mu2 = mu2';
ap.k  = k';
ap.eps0 = eps0';
ap.dx  = dx';
ap.dy  = dy'; 
ap.adx = adx'; 
ap.ady = ady';
ap.drnd = drnd;

