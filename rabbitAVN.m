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
clear all

global tN; 
tN = 1;  % pulse counter

Tictoc = 0; 
model_name = 'avn_data1';

%%%%%%  Select Mode from the list below: %%%%
Mode = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S2_interval = 0.175; % Test S1-S2 interval (for a single point calculation) [s]

% Mode:  
% 1 -  normal sinus rhythm 166 bpm
% 2 -  tachycardia (with highest achievable sinus rhythm)
% 3 -  automaticity in AVN 100 bpm
% 20 - periodic pacing, Wenckebach 5:4
% 30 - uniform random pacing 75-150 ms, atrial fibrillation
% 40 - S1S2 protocol (Control),     single point, anterograde conduction
% 41 - S1S2 protocol (FP ablation), single point, anterograde conduction
% 42 - S1S2 protocol (SP ablation), single point, anterograde conduction
% 43 - S1S2 protocol (Control),     full conduction curve, anterograde conduction
% 44 - S1S2 protocol (FP ablation), full conduction curve, anterograde conduction
% 45 - S1S2 protocol (SP ablation), full conduction curve, anterograde conduction
%
% 50 - S1S2 protocol (Control),     single point, retrograde conduction
% 51 - S1S2 protocol (FP ablation), single pont, retrograde conduction
% 52 - S1S2 protocol (SP ablation), single pont, retrograde conduction
% 53 - S1S2 protocol (Control),     full conduction curve, retrograde conduction
% 54 - S1S2 protocol (FP ablation), full conduction curve, retrograde conduction
% 55 - S1S2 protocol (SP ablation), full conduction curve, retrograde conduction
%

f = str2func(model_name);     
[Str,ap]=f();             
 
S1_interval = 0.360;   % Basic S1-S1 interval [s]
nS1 = 10;              % Number of S1 pulses 
nS2 = 1;
nS = nS1+nS2;
NFP_cells = Str.iPB-Str.iAM3-1;      % Number of FP cells 
NSP_cells = Str.iSlow2-Str.iSlow1+1; % Number of SP cells 
DN_cells = NSP_cells-NFP_cells;

V_cut = 0.25;               % Upstroke voltage level [dimensionless]
min_interval(1:3) = 0.030;  % Minimal interval between laddergram points [s]
max_interval(1:2) = 0.120;  % Maximal interval between laddergram points [s]

tmp_distrib = 1;
dx1_12 = ap.dx(Str.iAbl,1);  % Store FP coupling at ablation point
dx2_12 = ap.dx(Str.iAbl,2);  % Store SP coupling at ablation point

Asymmetry = 1; % 0 = OFF
%%%%%%%%%%%%%%%%%%%%%%%%
% ap.dx(Str.iAbl,1)=0;   % Uncomment for Fast PW ablation
% ap.dx(Str.iAbl,2)=0;   % Uncomment for Slow PW ablation
%%%%%%%%%%%%%%%%%%%%%%%%

if Asymmetry == 0
    ap.adx(:,1:2) = 1.0;  % Sets all asymmetry coefs to 1
    ap.ady(:) = 1.0;
    model_name = [model_name 'na'];
end

fprintf("= Dataset %s, Mode=%d: ",model_name(1:end),Mode);
switch Mode
case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normal sinus rhythm
    fprintf("Sinus rhythm ==\n");
    End_time = 4.0; % Total simulation time [s]
    plot_time.start = 3.26;
    plot_time.end = 3.82; 

 case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Highest sinus rhythm
    fprintf("Maximal sinus rhythm ==\n");
    End_time = 5.2; % Total simulation time [s]
    plot_time.start = 3.6;
    plot_time.end = 5.1;
    ap.dx(Str.iAbl,1) = dx1_12;  % No FP Ablation
    ap.dx(Str.iAbl,2) = dx2_12;  % No SP Ablation
    ap.a1(Str.iSN,1)  = Str.a1_max;

case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AVN automaticity
    fprintf("AVN automaticity ==\n");
    End_time = 5.1; % Total simulation time [s]
    plot_time.start = 3.6;
    plot_time.end = 5.1;
    ap.dx(Str.iAbl,1) = dx1_12;  % No FP Ablation
    ap.dx(Str.iAbl,2) = dx2_12;  % No SP Ablation
    ap.a1(Str.iSN,1)  = 0.001; 

case 20  %%%%%%%%%% Atrial Flutter/Fast regular pacing/Wenckebach
    fprintf("Atrial flutter ==\n");
    ap.dx(Str.iAbl,1) = dx1_12; 
    ap.dx(Str.iAbl,2) = dx2_12; 
    ap.dy(Str.imp_Antero) = 95; 
    period = 0.1250;      % Stimulation period [s]
    End_time = 7.0;       % Total simulation time [s]
    plot_time.start = 3.0;
    plot_time.end = 6.5; 
    tmp_distrib = NaN(1,[]);
    t_shift = 0; % Time shift of first pulse 
    i = 0; 
    tt = t_shift;
    while tt < End_time
        i = i+1;
        tmp_distrib(i) = tt;
        tt = tt + period;
        if(tt > End_time), break; end
    end
    min_interval(1:3) = 0.0370;  % Minimal interval between points [s]
    % 1 = SP, 2 = SN-AM, 3 = FP
    if ap.dx(Str.iAbl,2)==0, min_interval(1) = 0.0250; end 
    max_interval(1:2) = 0.120;   % Maximal interval between points [s]
    
case 30 %%%%%%%%%%%%%%%%%%%%%%%%  Atrial fibrillation
    fprintf("Atrial fibrillation ==\n");
    ap.dy(Str.imp_Antero) = 95;  % Sets coupling for AM* cell
    ap.a1(Str.iSN,1)  = -0.0800; % Sets higher sinus rate

    rng(ap.drnd(1),'v5uniform'); % Reset random sequence for reproducibility
    End_time = 11.0; %  Total simulation time [s]
    plot_time.start = 5.5;
    plot_time.end = 10.5;
    t_shift = 0; % Time shift of first pulse 
    tmp_distrib = NaN(1,[]);
    i = 0; 
    tt = t_shift;
    while tt < End_time
        tt = tt + random('Uniform',ap.drnd(2),ap.drnd(3));
        i = i+1;
        tmp_distrib(i) = tt;
    end
    min_interval(1) = 0.0370; % Minimal interval between points in SP [s]
    if ap.dx(Str.iAbl,2) == 0, min_interval(1) = 0.0250; end;
    min_interval(2) = 0.0220;  % in SN-AM
    min_interval(3) = 0.0250;  % in FP
    max_interval(1:2) = 0.120;   % Maximal interval between points [s]
 

case {40,41,42,43,44,45}    %%%%%%%%%%%  Single/Multiple Atrial S1S2 pacing
    name = ' Control ';
    if Mode == 41 || Mode == 44      % FP ablation
        ap.dx(Str.iAbl,1)=0; 
        name = ' FP ablation ';
    elseif Mode == 42 || Mode == 45  % SP ablation
        ap.dx(Str.iAbl,2)=0;  
        name = ' SP ablation ';
    end
    if Mode < 43 
        name = ['Single point, Atrial pacing,' name];
    else
        name = ['Conduction curve, Atrial pacing,' name];
        if ap.dx(Str.iAbl,2) ~= 0
           S2_interval = [0.36 0.34 0.32 0.30 0.28 0.26 0.24 0.22 0.20 ...
               0.18 0.17 0.16 0.15 0.145 0.140 0.135 0.130 0.125 0.120 ...
               0.115 0.110 0.105 0.100  0.098  0.096 0.094 0.092];         
        else
           S2_interval = [0.36 0.34 0.32 0.30 0.28 0.26 0.24 0.22 0.20 ...
            0.18 0.17 0.16 0.15 0.145 0.140 0.135 0.132 0.130 0.128 0.126];         
        end
        Condtime_array = NaN(1,length(S2_interval));
    end
    fprintf("%s \n",name);
    tmp_distrib = NaN(1,nS);
    t_shift = 0.0;
    tt = 0;
    for i = 1:nS
         if i <= nS1
            tt = t_shift + (i-1)*S1_interval;
         else 
            tt = tt + S2_interval(1);
         end
        tmp_distrib(i) = tt;
    end
    End_time        = tmp_distrib(end) +1.1; % Total simulation time [s]
    plot_time.start = tmp_distrib(end-1) -0.04;
    plot_time.end   = t_shift + 3.8; 

case {50,51,52,53,54,55}  %%%%%%%%%%%%%%%%%% Single/Multiple HB S1S2 pacing
    name = ' Control ';
    if Mode==51 || Mode==54       % FP ablation 
        ap.dx(Str.iAbl,1) = 0; 
        name = ' FP ablation ';
    elseif Mode==52 || Mode==55   % SP ablation 
        ap.dx(Str.iAbl,2) = 0;  
        name = ' SP ablation ';
    end

    if Mode < 53 
        name = ['Single point, HB pacing,' name];
    else
        name = ['Conduction curve, HB pacing,' name];
        % S1-S2 intervals 
        S2_interval = [0.36 0.34 0.32 0.30 0.28 0.26 0.24 0.22 ...
                       0.20 0.19 0.18 0.17 0.165 0.16 0.155 0.15];
        Condtime_array = NaN(1,length(S2_interval));
    end
    fprintf("%s \n",name);
    tmp_distrib = NaN(1,nS);
    t_shift = 0.0;
    tt = 0;
    for i = 1:nS
         if i <= nS1
            tt = t_shift + (i-1)*S1_interval;
         else 
            tt = tt + S2_interval(1);
         end
        tmp_distrib(i) = tt; 
    end

    End_time        = tmp_distrib(end) + 1.4; % Total simulation time [s]
    plot_time.start = tmp_distrib(end-1) - 0.04; 
    plot_time.end   = t_shift + 3.8;

otherwise  %%%%%%%%%%%%%%% Error
    fprintf('Wrong Mode number\n');
    return
end


t_step = 0.1*1.0e-3;  % For post-analysis 
opts = odeset('Reltol',1e-7,'AbsTol',1e-10, 'Stats','off');
%opts = odeset('Reltol',1e-6,'AbsTol',1e-8, 'Stats','off', ...
%    'NormControl','on','InitialStep',1e-6);  % test
initcond = 0.0001.*ones(Str.iHB2*4,1);
tspan = [0 t_step End_time];
t1 = 0:t_step:End_time; 

NN = 1;
switch Mode % Conduction curves
    case {43,44,45,53,54,55}
       NN = length(S2_interval);
end

for nn = 1:NN  % Loop for multiple calculations
  switch Mode  % If conduction curves - change S1S2 interval
    case {43,44,45,53,54,55}
       tmp_distrib(end) = tmp_distrib(end-1) + S2_interval(nn);
       tN = 1;
  end

  if Tictoc, tic; end
  fprintf('= Solving ODEs ...\n');
%%%%%%%%%%%%%%%%% Solving ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  sol = ode23(@(t,y) avn_fun(t,y,Mode,Str,ap,tmp_distrib), ...
                     tspan, initcond, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  y11 = sol.y';
  y11 = reshape(y11,[length(y11),Str.iHB2,4]);
  if Tictoc, toc; end

  y1 = 0.001.*ones(Str.iHB2,length(t1));
  y2 = 0.001.*ones(Str.iHB2,length(t1));
  for i = 1:Str.iHB2
    y1(i,:) = interp1(sol.x,y11(:,i,1),t1); % u in 1st row  (FP)
    y2(i,:) = interp1(sol.x,y11(:,i,3),t1); % u in 2nd row  (SP)
  end

  y1c = zeros(Str.iHB2,length(t1)); % "Cut" arrays
  y2c = zeros(Str.iHB2,length(t1));

% Only for images of S1S2 stimulation pulses 
  if     Mode == 20 || Mode == 30,  ii = Str.imp_Antero-1;
  elseif Mode >= 40 && Mode <= 45,  ii = Str.imp_Antero;
  elseif Mode >= 50,                ii = Str.imp_Retro;   
  end 
  if Mode >= 20 
    for i = 1:length(t1)
      for j = 1:length(tmp_distrib)  
          if(t1(i) >= tmp_distrib(j) && t1(i) < tmp_distrib(j)+ ap.pulse_width*2)
                y2(ii,i) = 0.8;
          end
       end
    end    
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of excitation latencies (Fast/Slow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('= Creating laddergram and calculating latencies...\n');
  locs_F = NaN(Str.iHB2,[]);
  locs2  = NaN(Str.iHB2,[]);
  locs_S = NaN(Str.iHB2+DN_cells,[]);
  alen = NaN(Str.iHB2,1);
  blen = NaN(Str.iHB2,1);

  locs_indexF = 1:Str.iHB2;
  locs_indexS = 1:Str.iHB2+DN_cells;

  y1c = min(V_cut,y1);    % Cut arrays, to get max V point
  y2c = min(V_cut,y2);

  for i = 1:Str.iHB2
    MPH = V_cut-0.01; 
    MPP = V_cut-0.01; 
    if i == Str.iSN,  MPP = V_cut*2/3; end

    [pks,lcs] = findpeaks(y1c(i,1:end),t1(1:end),'MinPeakHeight',MPH, ...
            'MinPeakProminence',MPP,'MinPeakDistance',0.03);
    if ~isnan(lcs)
        alen(i) = length(lcs);
        locs_F(i,1:length(lcs)) = lcs();
    else
        fprintf('No peaks (Fast) found\n');
        %return
    end
  end
 for i = Str.iSlow1:Str.iSlow2
    [pks2,lcs2] = findpeaks(y2c(i,1:end),t1(1:end),'MinPeakHeight',MPH, ...
            'MinPeakProminence',MPP,'MinPeakDistance',0.03);
    if ~(lcs2)
        fprintf('No peaks (Slow) found\n');
    %    return
    end
    blen(i) = length(lcs2);
    locs2(i,1:length(lcs2)) = lcs2();
  end

  al = max(alen);
  bl = max(blen); 
  locs_S(1:Str.iAM3,1:al)                = locs_F(1:Str.iAM3,1:al);
  locs_S(Str.iSlow1+1:Str.iSlow2+1,1:bl) = locs2(Str.iAM3:Str.iPB,1:bl);
  locs_S(Str.iSlow2+2:Str.iHB2+2,1:al)   = locs_F(Str.iPB:Str.iHB2,1:al);

  for i = Str.iAM3+1:Str.iPB+2
    locs_indexS(i) = Str.iAM3 +(i-Str.iAM3)/(Str.iPB-Str.iAM3+2)*(Str.iPB-Str.iAM3);
  end
  for i = Str.iPB+3:Str.iHB2+2
    locs_indexS(i) = locs_indexF(i-2);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Mode <= 30 
    start_SN = floor(length (nonzeros(locs_F(Str.iSN,:)))/2.0);
    start_PB = floor(length (nonzeros(locs_F(Str.iPB,:)))/2.0);
    start_PS = floor(length (nonzeros(locs_F(Str.ipSlow,:)))/2.0);
    start_HB = floor(length (nonzeros(locs_F(Str.iHB2,:)))/2.0);
    start_AM = floor(length (nonzeros(locs_F(Str.iAM3,:)))/2.0);

    Period_SN = mean(diff(nonzeros(locs_F(Str.iSN,start_SN:end-1))));
    Bpm_SN = 60.0/Period_SN; 
    Period_PB = mean(diff(nonzeros(locs_F(Str.iPB,start_PB:end-1))));
    Bpm_PB = 60.0/Period_PB;
    Period_HB2 = mean(diff(nonzeros(locs_F(Str.iHB2,start_HB:end-1))));
    Bpm_HB2 = 60./Period_HB2;
    Period_pSlow = mean(diff(nonzeros(locs2(Str.ipSlow,start_PS:end-1))));
    Bpm_pSlow = 60./Period_pSlow;
    Period_AM3 = mean(diff(nonzeros(locs2(Str.iAM3,start_AM:end-1))));
    Bpm_AM3 = 60./Period_AM3;

    fprintf('Rates at: SN=%0.1f AM3=%0.1f SP7=%0.1f PB=%0.1f HB6=%0.1f [bpm]\n',...
             Bpm_SN,Bpm_AM3,Bpm_pSlow,Bpm_PB,Bpm_HB2);
    fprintf('Ratios: SN/HB=%0.2f Atria/HB=%0.2f\n', Bpm_SN/Bpm_HB2,Bpm_AM3/Bpm_HB2);
  end

% Get conduction latencies 
  iTime.pointA  = 0;
  iTime.pointB  = 0;
  S1_condtime = 0;
  S2_condtime = 0;

  if Mode ~= 30
    iend = min(size(locs_F,2),size(locs2,2))-1;
    if Mode >= 40 && Mode < 50 
      for i = 2:size(locs_F,2)
          if abs(locs_F(Str.imp_Antero,i) - tmp_distrib(end)) < min_interval
              iend = i;
              break;
          end
      end
    elseif Mode >= 50 
      for i = 2:min(size(locs_F,2),size(locs2,2))
          if abs(locs_F(Str.imp_Retro,i) - tmp_distrib(end)) < min_interval
              iend = i;
              break;
          end
      end
    end

    if Mode ==1 %<= 2

      e1 = nnz(locs_F(Str.iSN,:))- nnz(locs_F(Str.pointA,:));
      e2 = nnz(locs_F(Str.iSN,:))- nnz(locs_F(Str.iPB,:));
      e3 = nnz(locs_F(Str.iSN,:))- nnz(locs_F(Str.pointB,:));

      iTime.pointA = 1.e3*(locs_F(Str.pointA,iend-1-e1) - locs_F(Str.iSN,iend-1) );
      iTime.PB     = 1.e3*(locs_F(Str.iPB,iend-e2)    - locs_F(Str.iSN,iend) );
      iTime.pointB = 1.e3*(locs_F(Str.pointB,iend-e3) - locs_F(Str.iSN,iend) );
      iTime.HB2    = 1.e3*(locs_F(Str.iHB2,iend-e3)   - locs_F(Str.iSN,iend) );
      iTime.Slow_max = 1.e3*(max(locs2(Str.iSlow1:Str.iSlow2,iend-e3)) - locs_F(Str.iSN,iend) );
    end

    S1_condtime = 1.e3*(locs_F(Str.pointB,iend-1) - locs_F(Str.pointA,iend-1));
    S2_condtime = 1.e3*(locs_F(Str.pointB,iend)   - locs_F(Str.pointA,iend));
    Delta = S2_condtime-S1_condtime;
    if Mode >= 50
      S1_condtime = -S1_condtime; S2_condtime = -S2_condtime; 
      Delta = -Delta;
      if Delta < -0.9, S2_condtime = NaN; Delta = NaN; end
    end
    if (S2_condtime > 200 || S2_condtime < 40), S2_condtime = NaN; Delta = NaN; end
  
    if Mode >= 40
      fprintf(' S1/S2=%0.0f/%0.0f  S1_condtime=%0.2f S2_condtime=%0.2f [ms]\n', ...
          1.e3*( tmp_distrib(end-1)- tmp_distrib(end-2)), ...    
          1.e3*( tmp_distrib(end)- tmp_distrib(end-1)), ...
          S1_condtime,S2_condtime);
    elseif  Mode ==1 % <= 2
      fprintf('Latency: pointA=%0.1f max_SP=%0.1f PB=%0.1f pointB=%0.1f Total=%0.1f S2_conduction=%0.1f [ms]\n',...
              iTime.pointA,iTime.Slow_max,iTime.PB,iTime.pointB,iTime.HB2,S2_condtime);
    end
 end


%%%%%%%%%%%%%%%%%% Plot of Laddergrams and Action Potentials %%%%%%%%%%%%%

  avn_plot(Mode,Str,ap,al,bl,t1,y1,y2,locs_F,locs_S,locs_indexS, ...
            min_interval,max_interval,tmp_distrib,iTime,nS,plot_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  switch Mode
    case {43,44,45,53,54,55} 
    Condtime_array(nn) = S2_condtime;
  end

end %%%%%%%%%%% END of Loop %%%%%%%%%%%%%%%

% Conduction time plot
switch Mode 
    case {43,44,45,53,54,55}
    figure(5);
    plot(S2_interval*1.0e3,Condtime_array,'--o');
    xlabel('S1S2 interval (ms)');
    ylabel('Conduction time (ms)');
    set(gca,'FontSize',10,'FontName','Arial');
    title(name);
    fname = ['Mode' num2str(Mode,'%d') '_data.txt'];
    fileID = fopen(fname,"w");
    fprintf(fileID,"S1S2  Cond_time [ms]\n");
    for i = 1:length(Condtime_array)
        if ~isnan(Condtime_array(i))
           fprintf(fileID,"%0.1f  %0.1f\n",S2_interval(i)*1.0e3,Condtime_array(i));
        end
    end
    fclose(fileID);
end
