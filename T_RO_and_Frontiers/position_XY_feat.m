%% load roomba data and save in standard format
%       cell array: {trials, classes}
%       elements array: [time, stream]

%% example analysis
%while(0);
  clear all; clc; %path(path,'../functions'); 

  load('data_XY','data'); 
  [ni,nc1,nc2] = size(data); [nt,nn] = size(data{1,1}); % 4 by 26 by
                                                        % 101. 750 by 2
  tic
    ri = 1:ni; rc1 = 1:1:nc1; rc2 = 1:1:nc2; % was every 4th trial (radius)
    nc1 = length(rc1); nc2 = length(rc2); nc = nc1*nc2; rc = 1:nc;

    % onsets
    for c1 = rc1
      for c2 = rc2
        for i = ri
          t(i,c1,c2) = min(find(data{i,c1,c2}(:,1)>0.1))-100;
        end
      end
    end



    % MAT need it to be separate files, not concatenated. ------

    for i = 1:3;
      c= 0;
      for c1 = rc1
        for c2 = rc2
          c = c+1;
          rt = t(i+1,c1,c2)+(1:750);
          data_train{i,c} = data{i+1,c1,c2}(rt,1);
          indRadius(i,c) = c2;
          indVelocity(i,c) = c1;
        end
      end
    end
    % ----------------------------------------------------------

    range = [min(min(cell2mat(data_train))) ...
             max(max(cell2mat(data_train)))];


    % test data: truncate 
    c = 0; 
    for c1 = rc1
      for c2 = rc2
        c = c + 1; data_test{c} = [];
        for i = 1 
          rt = t(i,c1,c2)+(1:750);
          data_test{c} = [data_test{c};data{i,c1,c2}(rt,1)];
          test(c) = c; 
        end
      end
    end
    % =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    % for looping
    clear t
    ratios = [0.5 0.1 0.2 0.4 0.6 0.8 1]; % Different training set sizes
    dims = [nc1, nc2];
    for m = 7;  %1:7; 
      Zipp = ratios(m);
      for n = 1; %:5;
        m
        n
        
        % =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+


        % =============================================
        % For running on reduced training data
        clear sub_data_train;
        clear sub_indRadius;
        clear sub_indVelocity;



        tot = nc1*nc2;
        sub = ceil(Zipp*tot); % subset of total;
        totR = randperm(tot);
        subN = totR(1:sub);
        sortSubN = sort(subN);
        for i = 1:sub;
          sub_data_train{1,i} = data_train{1,sortSubN(i)};
          sub_data_train{2,i} = data_train{2,sortSubN(i)};
          sub_data_train{3,i} = data_train{3,sortSubN(i)};
          sub_indRadius(:,i) = indRadius(:,sortSubN(i));
          sub_indVelocity(:,i) = indVelocity(:,sortSubN(i));
        end
        out = train_feature(sub_data_train,sub_indRadius,sub_indVelocity,dims);
        % =============================================



        % MAT training: feature classifier ------------------------
        %out = train_feature(data_train,indRadius,indVelocity,dims);
        % ---------------------------------------------------------

        coeffs = out{1,1};
        coeffsR = coeffs(1,:);
        coeffsS = coeffs(2,:);

        %features = out{2,1};
        % useful label for plotting 
        % MAT: like meshgrid(1:26,1:26)
        mc1 = ceil(reshape(1:nc,[nc1,nc2])/nc2); 
        mc2 = 1+rem(reshape(0:nc-1,[nc1,nc2]),nc2);
        mc = [mc1(:),mc2(:)];

        % testing: feature classifier  ----------------------------
        clear class
        for c1 = rc; 
          class{c1} = test_feature(data_test{c1},coeffs); 
        end
        class = cell2mat(class);

        %end;

        radius = class(1,:); 
        speed = class(2,:);

        % bound max and min of reports to reasonable values
        for i = rc;
          if radius(i)<1;
            radius(i) = 1;
          end;
          if radius(i)>101;
            radius(i) = 101;
          end;
        end;

        for i = rc;
          if speed(i)<1;
            speed(i) = 1;
          end;
          if speed(i)>26;
            speed(i) = 26;
          end;
        end;
        % --------------------------------------------------------
      
      end;
    end;
    
 
%end;


%      for i = 1:2626;
%        radius(i) = indRadius(class(i));
%        speed(i) = indVelocity(class(i));
%      end

      
      
      errorR = radius-indRadius(1,:);
      errorV = speed-indVelocity(1,:);

%while(0);
        %==========
        stat(1,m,n) = mean(errorV);
        stat(2,m,n) = std(errorV);
        stat(3,m,n) = mean(errorR);
        stat(4,m,n) = std(errorR);
      toc;
      t(m*5-5+n) = toc;
%end;

%load statFeat
%clear statM; 
%for i = 1:4; 
%  for j = 1:6; 
%    statM(i,j) = mean(stat(i,j,:));
%  end;
%end
%statM
%==========

%radii = (radius+3)/4; 

%imagesc(reshape(radius,26,26));


% correct = zeros(26,26);
% for c1 = rc1
%     for c2 = rc2
%       correct(c1,(c2+3)/4) = c2;
%     end
% end

%for c1 = rc; 
%  error{c1} = test_error(data_test{c1},features); 
%end

%R = zeros(1,rc);
%S = zeros(1,rc);
%for c = rc;  
%      err = error{c};
%      [ig,classR(c)] = min(err(:,1)); 
%      [ig,classS(c)] = min(err(:,2)); 
%    % plot results
%    figure(3); clf; hold on; colormap(hot); 
%    clims1 = [min(err(:,1)) max(err(:,1))];
%    clims2 = [min(err(:,2)) max(err(:,2))];
%    f = 1/4; clims1 = [f*clims1(1)+(1-f)*clims1(2) clims1(2)]; 
%    f = 1/4; clims2 = [f*clims2(1)+(1-f)*clims2(2) clims2(2)]; 


%  For imaging regression output

%  imagesc(reshape(radius,[nc2,nc1]));
%  ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
%  plot((mc(c,2))*4-3,mc(c,1),'xg','Markersize',20); 
%  plot(radius(c),speed(c),'ob','Markersize',20);
%    err = reshape(err,26,26,2);
%      [meh R(c)] = min(min(err(:,:,1)));  
%      [meh S(c)] = min(min(err(:,:,1)'));

%    imagesc(-err(:,:,1));
%    ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
%    plot(mc(c,2),mc(c,1),'xg','Markersize',20); 
%    plot(mc(classR(c),2),mc(classR(c),1),'og','Markersize',20); % Min error between input and stored features (in R only)

%    plot(R(c),S(c),'ob','Markersize',20); % Min error between input
% and stored features (in
% R and S)
%    plot(radii(c),speed(c),'or','Markersize',20);  % Feature class
% output. R
% corrected for
% reduced
% dataset

%    set(gca,'YTick',0:5:30,'YTickLabel',{'0','5','10','15','20','25','30'}); 
%    set(gca,'XTick',linspace(0,nc2,11),'XTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'}); 
%    axis([0.5 nc2+0.5 0.5 nc1+0.5])
%  pause(0.1);

%end

% save results
%save(mfilename,'class','logp','mc','rc','nc1','nc2')

%% plot histogram of results
%clear all
%load(mfilename,'class','logp','mc','rc','nc1','nc2')

% ======================================

%% construct histogram
%r = -30:30; 
%nc = length(mc);
%m1 = radii'-mc(rc,2); h1 = histc(m1,r); h1 = h1/sum(h1); r1 = 26*r/nc1;
%m2 = speed'-mc(rc,1); h2 = histc(m2,r); h2 = h2/sum(h2); r2 = 101*r/nc2;

%%m1 = mc(classR,1)-mc(rc,1); h1 = histc(m1,r); h1 = h1/sum(h1); r1 = 26*r/nc1;
%%m2 = mc(classR,2)-mc(rc,2); h2 = histc(m2,r); h2 = h2/sum(h2); r2 = 101*r/nc2;

%mu1 = mean(m1); si1 = std(m1); mu2 = mean(m2); si2 = std(m2);

%% plot results
%figure(3); clf 
%let = {'\bf A','\bf B','\bf C','\bf D','\bf E','\bf F','\bf G','\bf H'};

%subplot(2,2,1); grid on; box on; hold on
%bar(r2,100*h2,'r'); 
%axis([-101 101 0 15]); set(gca,'Fontsize',6)
%xlabel('error (mm)','Fontsize',7); ylabel('frequency (%)','Fontsize',7); grid on; box on;
%title('Distance errors','Fontsize',8)
%text(0.6,0.9,{['\mu=',num2str(mu2,'%3.1f'),'mm'],['\sigma=',num2str(si2,'%3.1f'),'mm']},'Units','Normalized','Fontsize',6)
%xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
%set(gca,'position',[xpos(1) ypos(1) 0.38 0.36],'units','normalized')    
%text(-0.275,1.10,let{1},'Units','Normalized','Fontsize',10)

%subplot(2,2,2); grid on; box on; hold on
%bar(r1,100*h1,'r'); 
%axis([-26 26 0 30]); set(gca,'Fontsize',6)
%xlabel('error (mm/sec)','Fontsize',7); ylabel('frequency (%)','Fontsize',7); 
%title('Speed errors','Fontsize',8)
%text(0.6,0.9,{['\mu=',num2str(mu1,'%3.1f'),'mm/s'],['\sigma=',num2str(si1,'%3.1f'),'mm/s']},'Units','Normalized','Fontsize',6)
%xpos = [0.10 0.60]; ypos = [0.59,0.08]; 
%set(gca,'position',[xpos(2) ypos(1) 0.38 0.36],'units','normalized')    
%text(-0.275,1.10,let{2},'Units','Normalized','Fontsize',10)

%subplot(2,2,[3,4]); grid on; box on; hold on
%plot(101/nc2*m2+randn(nc,1),26/nc1*m1+randn(nc,1),'.r','Markersize',5); 
%axis([-101 101 -26 26]); set(gca,'Fontsize',6)
%set(gca,'XTick',linspace(-100,100,9)); 
%title('Scatterplot of speed-distance classification errors','Fontsize',8)
%ylabel('speed error (mm/s)','Fontsize',7); xlabel('distance error (mm)','Fontsize',7)
%xpos = [0.10 0.57]; ypos = [0.59,0.08]; 
%set(gca,'position',[xpos(1) ypos(2) 0.88 0.36],'units','normalized')    
%text(-0.12,1.10,let{3},'Units','Normalized','Fontsize',10)

% ======================================

%%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 4]);
%%print('-r300','-djpeg',mfilename) 

%% plot histogram
%%c = 5*101+20;
%%figure(3); clf; hold on; box on; grid off; colormap(hot); 
%%f = 1/4; clims = [f*clims(1)+(1-f)*clims(2) clims(2)]; 
%%imagesc(reshape(logp{c},[nc2,nc1])',clims)
%%ylabel('speed (mm/sec)','Fontsize',10); xlabel('distance (mm)')
%%plot(mc(c,2),mc(c,1),'xg','Markersize',20); 
%%plot(mc(class(c),2),mc(class(c),1),'ob','Markersize',20); 
%%set(gca,'YTick',0:5:30); 
%%set(gca,'XTick',0:10:100); 
%%axis([0.5 101.5 0.5 26.5]); set(gca,'Fontsize',8)
%%text(-0.0175,-0.0275,'0','Units','Normalized','Fontsize',8)
%%title('Compound eye of speed-distance classification','Fontsize',10)
%%set(gca,'position',[0.06 0.11 0.92 0.82],'units','normalized')    

%%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 3.5]);
%%print('-r300','-djpeg',[mfilename,'_1']) 
