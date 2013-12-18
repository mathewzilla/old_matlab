%% load roomba data and save in standard format
%       cell array: {trials, classes}
%       elements array: [time, stream]

%% example analysis

clear all; clc; %path(path,'../functions'); 

load('../data_XY','data'); 
[ni,nc1,nc2] = size(data); [nt,nn] = size(data{1,1});

ri = 1:ni; rc1 = 1:1:nc1; rc2 = 1:3:nc2; % every 4th trial
nc1 = length(rc1); nc2 = length(rc2); nc = nc1*nc2; rc = 1:nc;

% onsets
for c1 = rc1
    for c2 = rc2
        for i = ri
            t(i,c1,c2) = min(find(data{i,c1,c2}(:,1)>0.1))-100;
        end
    end
end

% smoothed velocity
g = @(n,si) exp(-(-n:n).^2/(2*si^2))/sum( exp((-n:n).^2/(2*si^2)) );   
for c1 = rc1
    for c2 = rc2 
        for i = ri 
            vdata{i,c1,c2} = diff(data{i,c1,c2}); 
            vdata{i,c1,c2} = [vdata{i,c1,c2}(1,:);vdata{i,c1,c2}];
            vdata{i,c1,c2} = filter(g(50,15)/sum(g(50,15)),1,vdata{i,c1,c2});
        end
    end
end

% training data: truncate and concatenate (3 trials); just x data
c = 0;
for c1 = rc1
    for c2 = rc2
        c = c + 1; data_train{c} = [];
        for i = 2:4
            rt = t(i,c1,c2)+(1:750);
            data_train{c} = [data_train{c};vdata{i,c1,c2}(rt,1)];
        end
    end
end
range = [min(min(cell2mat(data_train))) max(max(cell2mat(data_train)))]

% test data: truncate 
c = 0; 
for c1 = rc1
    for c2 = rc2
        c = c + 1; data_test{c} = [];
        for i = 1 
            rt = t(i,c1,c2)+(1:750);
            data_test{c} = [data_test{c};vdata{i,c1,c2}(rt,1)];
            test(c) = c; 
        end
    end
end

% training: naive X 
d = linspace(-0.1,0.1,501); d = d(:);
ns = 10; %ns = []; % smooth over 10 samples
for c = rc
    logl{c} = train_SNB(data_train{c},d,ns); 
end

% useful label for plotting
mc1 = ceil(reshape(1:nc,[nc1,nc2])/nc2); 
mc2 = 1+rem(reshape(0:nc-1,[nc1,nc2]),nc2);
mc = [mc1(:),mc2(:)];

% testing - plot results
for c = rc
    for c1 = rc; logp{c}(c1) = test_SNB(data_test{c},d,logl{c1}); end
    [ig,class(c)] = max(squeeze(logp{c}),[],2); 
    
    % plot results
    figure(3); clf; hold on; colormap(hot); 
    clims = [min(min(logp{c})) max(max(logp{c}))];
    f = 1/4; clims = [f*clims(1)+(1-f)*clims(2) clims(2)]; 
    imagesc(reshape(logp{c},[nc2,nc1])',clims)
    ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
    plot(mc(c,2),mc(c,1),'xg','Markersize',20); 
    plot(mc(class(c),2),mc(class(c),1),'ob','Markersize',20);
    set(gca,'YTick',0:5:30,'YTickLabel',{'0','5','10','15','20','25','30'}); 
    set(gca,'XTick',linspace(0,nc2,11),'XTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'}); 
    axis([0.5 nc2+0.5 0.5 nc1+0.5])
end

% save results
save(mfilename,'class','logp','mc','rc','nc1','nc2')

%% plot histogram of results
clear all
load(mfilename,'class','logp','mc','rc','nc1','nc2')

% construct histogram
r = -30:30; 
nc = length(mc);
m1 = mc(class,1)-mc(rc,1); h1 = histc(m1,r); h1 = h1/sum(h1); r1 = 26*r/nc1;
m2 = mc(class,2)-mc(rc,2); h2 = histc(m2,r); h2 = h2/sum(h2); r2 = 101*r/nc2;
mu1 = mean(m1); si1 = std(m1); mu2 = mean(m2); si2 = std(m2);

% plot results
figure(1); clf 
let = {'\bf A','\bf B','\bf C','\bf D','\bf E','\bf F','\bf G','\bf H'};

subplot(2,2,1); grid on; box on; hold on
bar(r2,100*h2,'r'); 
axis([-101 101 0 15]); set(gca,'Fontsize',6)
xlabel('error (mm)','Fontsize',7); ylabel('frequency (%)','Fontsize',7); grid on; box on;
title('Distance errors','Fontsize',8)
text(0.6,0.9,{['\mu=',num2str(mu2,'%3.1f'),'mm'],['\sigma=',num2str(si2,'%3.1f'),'mm']},'Units','Normalized','Fontsize',6)
xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
set(gca,'position',[xpos(1) ypos(1) 0.38 0.36],'units','normalized')    
text(-0.275,1.10,let{1},'Units','Normalized','Fontsize',10)

subplot(2,2,2); grid on; box on; hold on
bar(r1,100*h1,'r'); 
axis([-26 26 0 30]); set(gca,'Fontsize',6)
xlabel('error (mm/sec)','Fontsize',7); ylabel('frequency (%)','Fontsize',7); 
title('Speed errors','Fontsize',8)
text(0.6,0.9,{['\mu=',num2str(mu1,'%3.1f'),'mm/s'],['\sigma=',num2str(si1,'%3.1f'),'mm/s']},'Units','Normalized','Fontsize',6)
xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
set(gca,'position',[xpos(2) ypos(1) 0.38 0.36],'units','normalized')    
text(-0.275,1.10,let{2},'Units','Normalized','Fontsize',10)

subplot(2,2,[3,4]); grid on; box on; hold on
plot(101/nc2*m2+randn(nc,1),26/nc1*m1+randn(nc,1),'.r','Markersize',5); 
axis([-101 101 -26 26]); set(gca,'Fontsize',6)
set(gca,'XTick',linspace(-100,100,9)); 
title('Scatterplot of speed-distance classification errors','Fontsize',8)
ylabel('speed error (mm/s)','Fontsize',7); xlabel('distance error (mm)','Fontsize',7)
xpos = [0.10 0.57]; ypos = [0.59,0.08]; 
set(gca,'position',[xpos(1) ypos(2) 0.88 0.36],'units','normalized')    
text(-0.12,1.10,let{3},'Units','Normalized','Fontsize',10)

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 4]);
print('-r300','-djpeg',mfilename) 

% plot histogram
c = 5*101+20;
figure(3); clf; hold on; box on; grid off; colormap(hot); 
clims = [min(min(logp{c})) max(max(logp{c}))];
f = 1/4; clims = [f*clims(1)+(1-f)*clims(2) clims(2)]; 
imagesc(reshape(logp{c},[nc2,nc1])',clims)
ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
plot(mc(c,2),mc(c,1),'xg','Markersize',20); 
plot(mc(class(c),2),mc(class(c),1),'ob','Markersize',20); 
set(gca,'YTick',0:5:30); 
set(gca,'XTick',0:10:100); 
axis([0.5 101.5 0.5 26.5]); set(gca,'Fontsize',8)
text(-0.0175,-0.0275,'0','Units','Normalized','Fontsize',8)
title('Compound eye of speed-distance classification','Fontsize',10)
set(gca,'position',[0.06 0.11 0.92 0.82],'units','normalized')    

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 3.5]);
print('-r300','-djpeg',[mfilename,'_1']) 
