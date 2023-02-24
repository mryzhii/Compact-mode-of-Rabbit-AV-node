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
function avn_plot(Mode,S,ap,al,bl,t1,y1,y2,locs_F,locs_S,locs_indexS,...
    min_interval,max_interval,tmp_distrib,iTime,...
    nS,plot_time)

stack = dbstack('-completenames');
if numel(stack) <2
   fprintf('Wrong call. Start with ''rabbitAVN''\n');
   return
end

fprintf('= Analysing laddergram and plotting...  ');

fs1 = 10; fs2 = 10; 
fn1 = 'Arial';
main_col = '#888888';
pulse_col = '#000000';
osc_color = '#000000';
point_color = '#F0B0B0';
point_color1 = '#900000';

plot_height = 420; 
plot_scale = 1.0; 

Title0 = {'Sinus rhythm','Tachycardia','ANV automaticity','Atrial flutter / Wenckebach',...
          'Atrial fibrillation','Atrial pacing','HB pacing'};
Title1 = {', FP ablation',', SP ablation'};
lw1 = 0.5; 
ms = 1.5; % marker size
ms2 = 4;  % stimulation marker
msp = 3;  % pacemaker marker size 
mcol = '#555555'; % SN pacemaker marker color
lw2 = 0.3;
lw3 = 0.3;
NFP_cells = S.iPB-S.iAM3-1;
NSP_cells = S.iSlow2-S.iSlow1+1;
DN_cells = NSP_cells - NFP_cells;

%%%%%%% FIGURES %%%%%%%%%%%%%%%
f1 = figure(1);
clf;
axis tight
f1.Position = [30 277 450 plot_height]; 
%set(groot, 'defaulttextinterpreter',         'none');  
set(groot, 'defaultAxesTickLabelInterpreter','tex');  

nskip0 = 5;
xlim1 = plot_time.start;
xlim2 = plot_time.end; 
switch Mode
    case 1
        x_pos = -0.065;
        nskip  = 25; % skip points in AP plots
        xtk = 0.1;
        Title = Title0{1};
    case 2
        x_pos = -0.065;
        nskip  = 25;
        xtk = 0.2;
        Title = Title0{2};
    case 3
        x_pos = -0.065;
        nskip  = 25; % skip points in AP plots
        xtk = 0.2;
        Title = Title0{3};

    case {20,30}
        x_pos = -0.035;
        nskip  = 50;
        xtk = 0.5;
        f1.Position =      [30 277 900 plot_height];   
        Title = Title0{4};
        if Mode==30, Title = Title0{5}; end

    case {40,41,42,43,44,45}
        x_pos = -0.065;
        nskip = 25;
        xtk = 0.1;
        if( tmp_distrib(end)-tmp_distrib(end-1) >0.0995 && ...
               tmp_distrib(end)-tmp_distrib(end-1) <0.123  && Mode<43 )   % AVNRT
            x_pos = -0.035;
            xlim2 = plot_time.end +0.2; 
            f1.Position =      [30 277 900 plot_height]; 
        end

        Title = Title0{6};
        if (Mode==41 || Mode==42), Title = [Title Title1{Mode-40}]; end
        if (Mode==44 || Mode==45), Title = [Title Title1{Mode-43}]; end

    case {50,51,52,53,54,55}
        x_pos = -0.065;
        nskip = 25;
        xtk = 0.1;
        Title = Title0{7};
        if (Mode==51 || Mode==52), Title = [Title Title1{Mode-50}]; end
        if (Mode==54 || Mode==55), Title = [Title Title1{Mode-53}]; end
    otherwise
        x_pos = -0.065;
        nskip  = 25; 
        xtk = 0.2;
        Title = ' ';
end


startt = min(find(t1 > xlim1));
endt =   max(find(t1 < xlim2));

ysub1 = fix((S.iHB2+1)/1.0); 
ysub2 = fix((S.iHB2+1)*plot_scale); 
ysub3 = fix((S.iHB2+1+DN_cells)*plot_scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Top panel, ladder diagrams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = subplot(ysub1+ysub2+ysub3,1,1:ysub1); 
hold on
set(gcf,'color','w');
ylim([-(S.iHB2+1) 0]);
xlim([xlim1 xlim2]);
x_1 = [0 xlim2];
y_1 = [-S.pointA -S.pointA];
line(x_1,y_1,'Color',point_color,'LineWidth',0.5); 
y_1 = [-S.pointB -S.pointB];
line(x_1,y_1,'Color',point_color,'LineWidth',0.5); 
plot(xlim2,-S.pointA,'<','MarkerFaceColor',point_color1, ...
                   'MarkerEdgeColor',point_color1,'MarkerSize',ms2);
plot(xlim2,-S.pointB,'<','MarkerFaceColor',point_color1, ...
                   'MarkerEdgeColor',point_color1,'MarkerSize',ms2);

title(Title,'FontSize',fs1+1)
hh = gca;
set(gca,'YTickMode','manual');
hh.YAxis.TickLength=[0 0];
set(gca,'YTick',[-20 -14 -11 -6 -1]) 
set(gca,'YTickLabel',{'HB','\color{red}FP','\color{blue}SP','AM','SN'})
set(gca,'XTick',[0:xtk:xlim2]) 
set(gca,'XTickLabel',[]) ;
set(gca,'FontName',fn1);
set(gca,'FontSize',fs1);

% ---- Fast PW
pos_F = zeros([],1);
index = 0;
for  n = 1:S.iHB2-1
    for j = 1:al  
        for i = 1:al 
            if(n<S.iAM3)
                min_int = min_interval(3);  % min interval between AP ponts
            else
                min_int = min_interval(2); 
            end
            if abs(locs_F(n,j)-locs_F(n+1,i)) < min_int % select a pair of points
                x = [locs_F(n,j) locs_F(n+1,i)];
                y = [-n -n-1];
                if n < S.iAM3   % draw line in SN-AM part
                    line(x,y,'Color',main_col,'LineWidth',lw1,'Marker','o','MarkerSize',ms);
                    if (n==S.iSN && ap.a1(S.iSN,1)<0)  % if SN is active pacemeker cell
                        plot(locs_F(n,j),-S.iSN,'o','MarkerFaceColor',mcol, ...
                            'MarkerEdgeColor',mcol,'MarkerSize',msp); 
                    end
                elseif n >= S.iAM3 && n<=S.iPB   %  draw line in FP part
                    line(x,y,'Color','r','LineWidth',lw1,'Marker','o','MarkerSize',ms);
                    if (n==S.iPB && ap.a1(S.iPB,1)<0) ||  (n==S.iPB-1 && ap.a1(S.iPB-1,1)<0)
                      plot(locs_F(n,j),-n,'o','MarkerFaceColor','r', ...
                          'MarkerEdgeColor','r','MarkerSize',msp); 
                    end
                else            %  draw line in HB part 
                    if(Mode>=50) col = main_col;
                    else         col = 'r'; 
                    end
                    line(x,y,'Color',col,'LineWidth',lw1,'Marker','o','MarkerSize',ms);
                end
                if n == S.iPB-1  % if last FP earlier than PB
                    if x(1)~=0 && x(2)-x(1)>0
                        index= index +1;
                        pos_F(index)=x(1);
                    end
                end
            end
        end
    end
end


% ---- Slow PW
for i = 1:S.iHB2+1,  blen2(i) = size(locs_S,2);  end
pos_S = zeros([],1);
index = 0; 
for n = 1:S.iHB2+DN_cells-1
    for j = 1:al
        for i = 1:bl 
            if (abs(locs_S(n,j)-locs_S(n+1,i)) < min_interval(1)) % select a pair of points 
                x = [locs_S(n,j) locs_S(n+1,i)];
                y = [-locs_indexS(n) -locs_indexS(n+1)];
                if n>=S.iAM3 && n<S.iPB+2                      %  draw line in SP part
                      line(x,y,'Color','b','LineWidth',lw1,'Marker','o','MarkerSize',ms); 
                      if ( (n==S.ipSlow+1 && ap.a1(S.ipSlow,  2)<0) || ... % if pacemeker cell 
                           (n==S.ipSlow+2 && ap.a1(S.ipSlow+1,2)<0) || ... % in SP is active,
                           (n==S.ipSlow+3 && ap.a1(S.ipSlow+2,2)<0) || ... % draw big circle
                           (n==S.ipSlow+4 && ap.a1(S.ipSlow+3,2)<0) ) 
                          plot(locs_S(n,j),-locs_indexS(n),'o', ...
                              'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',msp); 
                      end
                end
                if n == S.iPB+1 % check 
                    if x(1) < x(2) 
                       index1 = 0;
                       for k = 1:length(pos_F)
                           a = abs(x(1)-pos_F(k));
                           if abs(x(1)-pos_F(k)) < max_interval(1)
                               if x(1)<pos_F(k)
                                  index = index +1;
                                  pos_S(index) = x(1);
                                  index1 = k+1;
                                  break;
                               else
                                   index1 = 1;
                                   break
                               end
                           end
                       end
                       if index1 == length(pos_F) || index1==0
                          index = index +1;
                          pos_S(index) = x(1);
                          % break;
                       end
                    end
                end

                if n >= S.iPB+2 && x(1)~=0  %  redraw line for SP in HB part if SP leads 
                    if(Mode>=50) col = main_col;
                    else col = 'b';
                    end
                    for k = 1:length(nonzeros(pos_S))
                       if x(1)-pos_S(k)>0 && x(1)-pos_S(k)  < max_interval(1)
                           line(x,y,'Color',col,'LineWidth',lw1,'Marker','o','MarkerSize',ms); 
                      end
                    end
                end
            end
        end
    end
end

% Mark external signal input place 
if  Mode >= 20 && Mode < 50
    m = ms2;
    if Mode <= 30, m = ms2-1; end
    plot(tmp_distrib-0.010,-S.imp_Antero,'>',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',m);
elseif Mode >= 50
    plot(tmp_distrib-0.010,-S.imp_Retro,'>',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',ms2);
end

% Numbers and labels on the ladder plot
switch Mode
case 1
    dx = locs_F(1,10)-locs_F(1,9);
    xpos = locs_F(1,10)-0.014-fs1*0.001; 
    ypos = -S.iSN;
    txt = [num2str(0,'%.0f')];
    text('Units','data','Position',[xpos ypos],'String',txt, ...
        'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

    txt = [num2str( 1.0/dx*60.0,'%.0f bpm')];
    ypos = -4;
    xpos = locs_F(1,10) +dx/2.0;
    text('Units','data','Position',[xpos ypos],'String',txt, ...
        'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

    txt = [num2str(iTime.pointA,'%.0f')];
    xpos = locs_F(S.pointA,10)-0.014-2*fs1*0.001;  
    ypos = -S.pointA;
    text('Units','data','Position',[xpos ypos],'String',txt, ...
        'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

    if(ap.dx(S.iAbl,1)~=0 && ap.dx(S.iAbl,2)~=0)
     txt = [num2str(iTime.Slow_max,'%.0f')];
     xpos = locs_F(S.pointB,10)+0.012;
     ypos = -S.pointB+4;
     text('Units','data','Position',[xpos ypos],'String',txt, ...
         'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);
    end

    txt = [num2str(iTime.PB,'%.0f')]; 
    xpos = locs_F(S.iPB,10)-0.016-length(txt)*fs1*0.001; 
    ypos = -S.iPB;
    text('Units','data','Position',[xpos ypos],'String',txt, ...
        'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

    xpos = locs_F(S.iHB2,10)-0.005+1*fs1*0.001; 
    ypos = -S.iHB2+1;
    txt = [num2str(iTime.HB2,'%.0f')]; %  ' ms'];
    text('Units','data','Position',[xpos ypos],'String',txt, ...
        'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

case  {2,3}
      ypos = -(S.pointB)+4.0;
      d = 7;
      if Mode == 2
          d = 24; 
          ypos = -(S.pointA)+4.0;
      end
      tt1 = locs_F(S.pointB,d);
      tt2 = locs_F(S.pointB,d+1);
      txt = [num2str( 1.0/(tt2-tt1)*60.0,'%.0f bpm')];
      xpos = tt1 + (tt2-tt1)/2 - length(txt)*0.01;
      text('Units','data','Position',[xpos ypos],'String',txt, ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

case {40,41,42,43,44,45}    
      txt = [num2str(round(1000.0*(tmp_distrib(end)-tmp_distrib(end-1)+0.0001)),'%.0f')]; 
      xpos = tmp_distrib(end-1) + (tmp_distrib(end)-tmp_distrib(end-1))/2 - length(txt)*0.01;
      ypos = -S.imp_Antero+2.5;
      text('Units','data','Position',[xpos ypos],'String',txt, ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);
      dx = 0;
      if(tmp_distrib(end)-tmp_distrib(end-1) == 0.110), dx = 0.005; end 
      xpos = tmp_distrib(end-1)-0.020-dx;
      text('Units','data','Position',[xpos ypos],'String','S1', ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);
      xpos = tmp_distrib(end)-0.020-dx;
      text('Units','data','Position',[xpos ypos],'String','S2', ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

      dx = 0.003;
      tt1 = locs_F(S.pointB,nS-1);
      nnS = size(locs_F,2);
      ypos = -S.imp_Retro+2.0;
      if nnS >= nS 
         tt2 = locs_F(S.pointB,nS);
         if(tt2>0 && tt2 < tmp_distrib(end)+0.200)
           txt = [num2str(round(1000.0*(tt2-tt1)+0.0001),'%.0f')]; 
           xpos = tt1 + (tt2-tt1)/2- length(txt)*0.01;
           text('Units','data','Position',[xpos ypos],'String',txt, ...
               'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);
           xpos = tt2-0.050-dx;
           text('Units','data','Position',[xpos ypos],'String','H2', ...
               'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);
        end
      end
      xpos = tt1-0.042-dx;
      text('Units','data','Position',[xpos ypos],'String','H1', ...
          'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);

case {50,51,52,53,54,55}
      dx = 0;
      txt = [num2str(round(1000.0*(tmp_distrib(end)-tmp_distrib(end-1)+0.0001)),'%.0f')]; 
      xpos = tmp_distrib(end-1) + (tmp_distrib(end)-tmp_distrib(end-1))/2 - length(txt)*0.005;
      ypos = -S.imp_Retro+2.5;
      text('Units','data','Position',[xpos ypos],'String',txt, ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);
      xpos = tmp_distrib(end-1)-0.025;
      text('Units','data','Position',[xpos ypos],'String','S1', ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);
      xpos = tmp_distrib(end)-0.025;
      text('Units','data','Position',[xpos ypos],'String','S2', ...
          'Color',[0.0,0.0,0.0],'FontSize',fs1,'FontName',fn1);

      tt1 = locs_F(S.pointA,nS-1);
      nnS = size(locs_F,2);
      ypos = -S.pointA+2;      
      
      if nnS > nS-1 
         tt2 = locs_F(S.pointA,nS);
         if(tt2>0 && tt2 < tmp_distrib(end)+0.200)
           txt = [num2str(round(1000.0*(tt2-tt1)+0.0001),'%.0f')]; 
           xpos = tt1 + (tt2-tt1)/2- length(txt)*0.01;
           text('Units','data','Position',[xpos ypos],'String',txt, ...
               'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);
           xpos = tt2-0.025-dx;
           text('Units','data','Position',[xpos ypos],'String','A2', ...
               'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);
        end
      end
      xpos = tt1-0.025-dx;
      text('Units','data','Position',[xpos ypos],'String','A1', ...
          'Color',[0.0,0.0,0.0],'fontsize',fs2,'fontname',fn1);
end

del = 0.3; % Distance between AP lines
lw = lw1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Middle panel, action potentials via Fast %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2 = subplot(ysub1+ysub2+ysub3,1,ysub1+1:ysub1+ysub2);
cla(h2);
hold on

switch Mode
case {1,2,3}
    ylim([-(S.iHB2+1)*del 1+1*del]);
case {20,30}  % Atrial flutter and Atrial fibrillation -> Stimuli and AM* on the top
    plot(t1(startt:nskip:endt), y2(S.imp_Antero,startt:nskip:endt)+del, ...
        'LineWidth',lw,'Color',main_col);
    plot(t1(startt:nskip0:endt),y2(S.imp_Antero-1,startt:nskip0:endt)+2*del, ...
        'LineWidth',lw,'Color',pulse_col);
    ylim([-(S.iHB2+1)*del 1+2*del]);
case {40,41,42,43,44,45}           % Anterograde -> Stimuli on the top
    plot(t1(startt:nskip0:endt),y2(S.imp_Antero,startt:nskip0:endt)+del, ...
        'LineWidth',lw,'Color',pulse_col);
    ylim([-(S.iHB2+1)*del 1+1*del]);    
case {50,51,52,53,54,55}        % Retrogarde -> Stimuli on the bottom
    plot(t1(startt:nskip0:endt),y2(S.imp_Retro,startt:nskip0:endt)-S.iHB2*del, ...
        'LineWidth',lw,'Color',pulse_col);
    ylim([-(S.iHB2+1)*del 1+del]);
end

for i = S.iSN : S.iHB2  %-1
                                          lw = lw1; col = main_col;
    if (i<=S.iAM3),                       lw = lw1; col = main_col; end
    if (i==S.pointA || i==S.pointB),      lw = lw1; col = point_color1; end
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(i-1), ...
        'LineWidth',lw,'Color',col);
end

xlim([xlim1 xlim2]);
hh = gca;
hh.YAxis.TickLength=[0 0];
set(gca,'YTickMode','manual');
set(gca,'YTick',[-5.7 -3.3  0.0]) 
set(gca,'YTickLabel',{'HB','\color{red}FP','SN'})
set(gca,'XTick',[0:xtk:xlim2]) 
set(gca,'XTickLabel',[]) ;
set(gca,'FontName',fn1);
set(gca,'FontSize',fs1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Bottom panel, action potentials via Slow %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = subplot(ysub1+ysub2+ysub3,1,ysub1+ysub2+1:ysub1+ysub2+ysub3);
cla(h3);
hold on

for i = S.iSN : S.iAM3
                                       lw = lw1; col = main_col;
    if (i<=S.iAM3),                    lw = lw1; col = main_col; end
    if (i==S.pointA),                  lw = lw1; col = point_color1;  end 
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(i-1), ...
        'LineWidth',lw,'Color',col);
    y21(i,startt:nskip:endt) = y1(i,startt:nskip:endt);
end
j1 = i;
for i = S.iSlow1:S.iSlow2
    lw = lw1; 
    col = main_col;
    plot(t1(startt:nskip:endt),y2(i,startt:nskip:endt)-del*(j1+i-S.iSlow1), ...
        'LineWidth',lw,'Color',col);
    y21(S.iAM3+1-S.iSlow1+i,startt:nskip:endt) = y2(i,startt:nskip:endt);
end
j2 = i;
for i = S.iPB : S.iHB2
    lw = lw1; 
    col = main_col;
    if  (i==S.pointB)  lw = lw1; col = point_color1;  end
    plot(t1(startt:nskip:endt),y1(i,startt:nskip:endt)-del*(j2+i-S.iPB+1), ...
        'LineWidth',lw,'Color',col);
end

ylim([-(S.iHB2+2)*del 1]);
xlim([xlim1 xlim2]);
hh = gca;
hh.YAxis.TickLength=[0 0];
set(gca,'YTickMode','manual');
set(gca,'YTick',[-6.3 -3.3  0.0]) 
set(gca,'YTickLabel',{'HB','\color{blue}SP','SN'})
set(gca,'XTick',[0:xtk:xlim2]) 
xlabel('Time (s)');
set(gca,'FontName',fn1);
set(gca,'FontSize',fs1);
drawnow

fprintf('Done\n');

