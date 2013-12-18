pos% plot histograms of classification output and save the figures


%% plot histogram of results
%clear all
%load(mfilename,'class','logp','mc','rc','nc1','nc2')

% construct histogram
r = -50:50; 
nc = length(mc);
m1 = errorV; % radii'-mc(rc,1); 
h1 = histc(m1,r); h1 = h1/sum(h1); r1 = 26*r/nc1;
m2 = errorR; % speed'-mc(rc,2); 
h2 = histc(m2,r); h2 = h2/sum(h2); r2 = 101*r/nc2;

%m1 = mc(classR,1)-mc(rc,1); h1 = histc(m1,r); h1 = h1/sum(h1); r1 = 26*r/nc1;
%m2 = mc(classR,2)-mc(rc,2); h2 = histc(m2,r); h2 = h2/sum(h2); r2 = 101*r/nc2;

mu1 = mean(m1)*7.2 
si1 = std(m1)*7.2 
mu2 = mean(m2) 
si2 = std(m2)

% plot results
%figure(5); clf;
let = {'\bf A','\bf B','\bf C','\bf D','\bf E','\bf F','\bf G','\bf H'};

figure(1); grid on; box on; hold on
bar(r2,100*h2,'r'); 
axis([-101 101 0 15]); set(gca,'Fontsize',16);axis square
xlabel('Error (mm)','Fontsize',17); ylabel('Frequency (%)','Fontsize',17); grid on; box on;
title('Distance errors','Fontsize',18)
text(0.6,0.9,{['\mu=',num2str(mu2,'%3.1f'),'mm'],['\sigma=',num2str(si2,'%3.1f'),'mm']},'Units','Normalized','Fontsize',16)
%xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
%set(gca,'position',[xpos(1) ypos(1) 0.38 0.36],'units','normalized')    
%text(-0.275,1.10,let{1},'Units','Normalized','Fontsize',10)
print('-depsc','featureHistR')

Xlabels = 7.2*(-20:20:20);
figure(2); grid on; box on; hold on
bar(r1,100*h1,'r'); %drawnow;
axis([-26 26 0 30]); set(gca,'Fontsize',16);axis square
set(gca,'XTickLabel',{'-140','-70','0','70','140'});
xlabel('error (mm/sec)','Fontsize',17); ylabel('Frequency (%)','Fontsize',17); 
title('Speed errors','Fontsize',18)
text(0.6,0.9,{['\mu=',num2str(mu1,'%3.1f'),'mm/s'],['\sigma=',num2str(si1,'%3.1f'),'mm/s']},'Units','Normalized','Fontsize',16)
%xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
%set(gca,'position',[xpos(2) ypos(1) 0.38 0.36],'units','normalized')    
%text(-0.275,1.10,let{2},'Units','Normalized','Fontsize',10)
Ylabels = 7.2*(-25:5:25);
print('-depsc','featureHistS')
%subplot(1,3,3); 
  figure(3);
  grid on; box on; hold on
  plot(101/nc2*m2+randn(1,nc),26/nc1*m1+randn(1,nc),'LineStyle','none','Marker', '.','color','r'); %,'MarkerEdgeColor' ,'r'  , ...
                                                                                                   %  'MarkerFaceColor' , 'w' ,'Markersize',4); %drawnow;
  axis([-101 101 -26 26]); set(gca,'Fontsize',16); axis square
  set(gca,'YTickLabel',Ylabels)
  set(gca,'XTick',linspace(-100,100,9)); 
  title('Scatterplot of speed-distance classification errors','Fontsize',18)
  ylabel('speed error (mm/s)','Fontsize',17); xlabel('distance error (mm)','Fontsize',17)
print('-depsc','featureHistScat')


while(0);

  %subplot(1,3,3); 
  figure(3);
  grid on; box on; hold on
  plot(101/nc2*m2+randn(1,nc),26/nc1*m1+randn(1,nc),'LineStyle','none','Marker', '.','color','r'); %,'MarkerEdgeColor' ,'r'  , ...
                                                                                                   %  'MarkerFaceColor' , 'w' ,'Markersize',4); %drawnow;
  axis([-101 101 -26 26]); set(gca,'Fontsize',16); axis square
  set(gca,'YTickLabel',Ylabels)
  set(gca,'XTick',linspace(-100,100,9)); 
  title('Scatterplot of speed-distance classification errors','Fontsize',18)
  ylabel('speed error (mm/s)','Fontsize',17); xlabel('distance error (mm)','Fontsize',17)





  %xpos = [0.10 0.57]; ypos = [0.59,0.08]; 
  %set(gca,'position',[xpos(1) ypos(2) 0.88 0.36],'units','normalized')    
  %text(-0.12,1.10,let{3},'Units','Normalized','Fontsize',10)

  %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 4]);
  %print('-r300','-djpeg',mfilename) 



  % % +++++++++++++++++++++++++++++
  % % Plotting error for different speeds (static method)
  piss = reshape(errorR,101,26);
  gay = mean(piss);
  meh = std(piss);
  figure(2);
  errorbar(gay,meh,'color','k', ...
           'Marker','o','MarkerEdgeColor','k', ...
           'MarkerFaceColor' , 'w' ,'Markersize',5);
  set(gca,'Fontsize',16);
  Ylabels = 7.2*(-25:5:20);
  set(gca,'YTick',(-25:5:20));
  set(gca,'YTickLabel',Ylabels);
  set(gca,'XTick',linspace(0,30,7));
  set(gca,'XTickLabel',{'0',' ','80',' ','160',' ','240'});
  axis square; axis([0 30 -25 20]);
  grid on;
  box on;
  xlabel('Speed')
  ylabel('Error')
  title('Distance errors','Fontsize',18)
  % +++++++++++++++++++++++++++++





  % plot histogram
  %c = 5*101+20;
  %figure(3); clf; hold on; box on; grid off; colormap(hot); 
  %f = 1/4; clims = [f*clims(1)+(1-f)*clims(2) clims(2)]; 
  %imagesc(reshape(logp{c},[nc2,nc1])',clims)
  %ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
  %plot(mc(c,2),mc(c,1),'xg','Markersize',20); 
  %plot(mc(class(c),2),mc(class(c),1),'ob','Markersize',20); 
  %set(gca,'YTick',0:5:30); 
  %set(gca,'XTick',0:10:100); 
  %axis([0.5 101.5 0.5 26.5]); set(gca,'Fontsize',8)
  %text(-0.0175,-0.0275,'0','Units','Normalized','Fontsize',8)
  %title('Compound eye of speed-distance classification','Fontsize',10)
  %set(gca,'position',[0.06 0.11 0.92 0.82],'units','normalized')    

  %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 3.5]);
  %print('-r300','-djpeg',[mfilename,'_1']) 

  



  % STATM PLOTS
  load statFreq2;
  clear statM; 
  for i = 1:4; 
    for j = 1:6; 
      statM(i,j) = mean(stat(i,j,:));
    end;
  end

  statM(:,7) = [-12.1 7.5 -5.43 36.22]; %Feat [ 0.4 3.6 4.6 11.9];
                                     %Temp [-0.17 2.18 0.4 4.92]; 
                                     %Freq[-12.1 7.5 -5.43 36.22];
                                     %NB[-0.47 1.78 0.5 4.9]

  figure(1); errorbar(statM(3,:),statM(4,:),'color','k', ...
                      'Marker','o','MarkerEdgeColor','k', ...
                      'MarkerFaceColor' , 'w' ,'Markersize',8); ...
      axis square; grid on; title('Distance error','Fontsize',18);
  xlabel('Training set size','Fontsize',18); 
  ylabel('Error','Fontsize',18);
  set(gca,'XTick',0:1:8,'XTickLabel',{' ','5%','10%','20%','40%', ...
                      '60%','80%','100%',' '});
  set(gca,'YTick',-50:10:50,'Fontsize',16);
  axis([0 8 -50 50]);

  figure(2); errorbar(7.2*statM(1,:),7.2*statM(2,:),'color','k', ...
                      'Marker','o','MarkerEdgeColor','k', ...
                      'MarkerFaceColor' , 'w' ,'Markersize',8); 
  axis square; grid on; title('Speed error','Fontsize',18);
  xlabel('Training set size','Fontsize',18); 
  ylabel('Error','Fontsize',18); 
  set(gca,'XTick',0:1:8,'XTickLabel',{' ','5%','10%','20%','40%', ...
                      '60%','80%','100%',' '},'Fontsize',16);
  set(gca,'YTick',-180:20:20,'Fontsize',16);
  axis([0 8 -180 20]);
end


while(0);
  % Comparative plots
  s = [ -89.98 -1.22 0.09  -3.38;53.99 15.70    25.80  12.79];

  r = [0.25   -1.78  0.4  1.65  0.47 ; 8.8  34.33  4.92    7.84  4.85]; 
  
  
  figure(1);errorbar(r(1,:),r(2,:),'LineStyle','none','color','k', ...
                     'Marker','o','MarkerEdgeColor','k', ...
                     'MarkerFaceColor' , 'b' ,'Markersize',8); ...
      axis square; grid on; title('Distance error','Fontsize',18);xlabel('Method','Fontsize',18); ylabel('Error','Fontsize',18);
  set(gca,'XTick',0:1:6,'XTickLabel', {' ','Static','Freq','Temp','Feat', ...
                      'NB',' '});
  set(gca,'YTick',-50:10:50);
  axis([0 6 -50 50]);set(gca,'Fontsize',16);
  
  figure(2);errorbar(s(1,:),s(2,:),'LineStyle','none','color','k', ...
                     'Marker','o','MarkerEdgeColor','k', ...
                     'MarkerFaceColor' , 'b' ,'Markersize',8); ...
      axis square; grid on; title('Speed error','Fontsize',18); ...
      xlabel('Method','Fontsize',18); ylabel('Error','Fontsize',18);
  set(gca,'YTick',-160:20:40);
  set(gca,'XTick',0:1:7,'XTickLabel', {' ','Freq','Temp','Feat', ...
                      'NB',' '});
  axis([0 5 -160 40]);set(gca,'Fontsize',16);

  % PLOTTING DATA FOR INSPECTION/PAPER
  indS = [0 44 88 132 176 220 264];        % Angle index
  sindS = [0 4 8 12 16 20 24 28 32 36 40]; % Speed Index
  ord = [1 3 5 7];
  figure(1);clf;
  for i = 1:4; %[1 3 5 7];
               % for j = [1 6 11];

    subplot(4,4,(4*i - 3)); plot(squeeze(trainSub{indS(ord(i))+sindS(1)+1})); hold all;

    plot(squeeze(trainSub{indS(ord(i))+sindS(6)+1}));
    plot(squeeze(trainSub{indS(ord(i))+sindS(11)+1}));
    axis([0 3000 -0.2 2]);%  axis square;
    xlabel('time(ms)','Fontsize',6); ylabel('Hall output','Fontsize',6); 
    title(['Texture: P80, Angle:',num2str(ord(i)*10,'%d'),' \circ'],'Fontsize',8); 
    %   legend ('S = 36mm/s','S = 72mm/s','S = 108mm/s'); 
    set(gca,'Fontsize',6); 
    %  legend('boxoff')
    

    subplot(4,4,(4*i - 2)); plot(squeeze(trainSub{indS(ord(i))+sindS(1)+2})); hold all;

    plot(squeeze(trainSub{indS(ord(i))+sindS(6)+2}));
    plot(squeeze(trainSub{indS(ord(i))+sindS(11)+2}));
    axis([0 3000 -0.2 2]);%  axis square;
    xlabel('time(ms)','Fontsize',6); ylabel('Hall output','Fontsize',6); 
    title(['Texture: P180, Angle:',num2str(ord(i)*10,'%d'),' \circ'],'Fontsize',8); 
    %legend ('S = 36mm/s','S = 72mm/s','S = 108mm/s'); 
    set(gca,'Fontsize',6); 
    %   legend('boxoff')

    subplot(4,4,(4*i - 1)); plot(squeeze(trainSub{indS(ord(i))+sindS(1)+3})); hold all;

    plot(squeeze(trainSub{indS(ord(i))+sindS(6)+3}));
    plot(squeeze(trainSub{indS(ord(i))+sindS(11)+3}));
    axis([0 3000 -0.2 2]);%  axis square;
    xlabel('time(ms)','Fontsize',6); ylabel('Hall output','Fontsize',6); 
    title(['Texture: P600, Angle:',num2str(ord(i)*10,'%d'),' \circ'],'Fontsize',8); 
    %legend ('S = 36mm/s','S = 72mm/s','S = 108mm/s'); 
    set(gca,'Fontsize',6); 
    %   legend('boxoff')

    subplot(4,4,(4*i)); plot(squeeze(trainSub{indS(ord(i))+sindS(1)+4})); hold all;

    plot(squeeze(trainSub{indS(ord(i))+sindS(6)+4}));
    plot(squeeze(trainSub{indS(ord(i))+sindS(11)+4}));
    axis([0 3000 -0.2 2]);%  axis square;
    xlabel('time(ms)','Fontsize',6); ylabel('Hall output','Fontsize',6); 
    title(['Texture: Smooth, Angle:',num2str(ord(i)*10,'%d'),' \circ'],'Fontsize',8);
    
    %legend ('S = 36mm/s','S = 72mm/s','S = 108mm/s'); 
    set(gca,'Fontsize',6);
    %legend('boxoff')
    
  end; print('-depsc','angletexsigsNoL')


  
  
  % Plotting errors with respect to speed/radius
  V = reshape(errorV,101,26);
  R = reshape(errorR,101,26);
  figure(1);
  %bar(mean(V),'b')
  errorbar(mean(V),std(V)/sqrt(101),'ob','MarkerEdgeColor','b','MarkerFaceColor','b');
  grid on
  axis([0 27 -13 13]);
  %axis([0 27 -4 2]); 
  set(gca,'Fontsize',16); axis square;
  %set(gca,'XTick',[0 5 10 15 20 25],'XTickLabel',[0,36,72,108,144,180,216]);
  %set(gca,'YTick',[-4.67 -2.78 -1.389 0 1.389 2.78],'YTickLabel',[-30,-20,-10,0,10,20]);
  set(gca,'XTick',[0 5 10 15 20 25],'XTickLabel',[0,36,72,108,144,180,216]);
  set(gca,'YTick',[-10 -5 0 5 10],'YTickLabel',[-72,-36,0,36,72]);
 
  
  xlabel('Contact speed (mm/sec)','Fontsize',17); ylabel('Error (mm/s)','Fontsize',17); 
title('Speed errors, with respect to true speed','Fontsize',18)
print('-depsc','speederrorwrtspeedEB')

figure(2);
  %bar(fliplr(mean(R')),'b')
  errorbar(fliplr(mean(R')),fliplr(std(R'))/sqrt(26),'ob','MarkerEdgeColor','b','MarkerFaceColor','b');
  grid on
  %axis([0 101 -10 8]); 
  axis([0 101 -50 50]); 
  set(gca,'Fontsize',16); axis square;
  set(gca,'XTick',[0 20 40 60 80 100],'XTickLabel',[90 110 130 150 170 190]);
  %set(gca,'YTick',[-10 -8 -6 -4 -2 0 2 4 6 8],'YTickLabel',[-10 -8 -6 -4 -2 0 2 4 6 8]);
  set(gca,'YTick',[-40 -20 0 20 40],'YTickLabel',[-40 -20 0 20 40]);
  xlabel('Radial distance to contact','Fontsize',17); ylabel('Error (mm)','Fontsize',17); 
title('Radius errors, with respect to true radius','Fontsize',18)
print('-depsc','radiuserrorwrtradiusEB')
end;