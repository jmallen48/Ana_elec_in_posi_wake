%% Margaux's Header, etc.

% clear all;
% 
% header = '/Users/margauxschmeltz/Desktop/SLAC/Data Analysis';
% prefix = '/Volumes/PWFA_4big';
% expt = 'E200';
% data_path = ['/nas/nas-li20-pm00/' expt '/2014/'];
% day = '20140629';
% data_set = 'E200_13543';

%% Margaux's Prefix, etc.
% user = 'e157user';
% 
% %%
% user = 'corde';
% 
% %% Load data
% 
% if ~exist(prefix)
%     system(['mkdir ' prefix]);
%     system(['/usr/local/bin/sshfs ' user '@quickpicmac3.slac.stanford.edu:/Volumes/PWFA_4big/' ' ' prefix]);
%     pause(2);
% end
% load([prefix data_path day '/' data_set '/' data_set '.mat']);
% cmap = custom_cmap()
%%
clear all;

%%
cmap = custom_cmap;

%% Data number (last two digits)
k = 43;

%% Load data
prefix = '/Users/jmallen/FACET/Electrons_in_Positron_Wake/Data_Analysis';
nas = '/nas/nas-li20-pm00';
experiment = '/E200';
year = '/2014';
date = '/20140629';
data_set = ['E200_135' num2str(k)];
path = [prefix nas experiment year date '/' data_set '.mat'];
load(path);

%% BetaL zoom

% Create ROI using pixel values
roi.top = 1;
roi.bottom = 170;
roi.left = 515;
roi.right = 685;

zoom_height = roi.bottom - roi.top + 1;
zoom_width = roi.right - roi.left + 1;

% Not really used...
% roi_process.top = 1;
% roi_process.bottom = 200;
roi_process.left = 90;
roi_process.right = 140;
process_width = roi_process.right - roi_process.left + 1;

pix_inf = 384.5;

%% BetaL images

% Choose what to plot. To plot set == 1, to skip set == 0.
plot_zoom = 0;
plot_process = 1;
plot_fullBETAL = 0;
shots_plot = 50; %number of shots plotted per step

POS_x = zeros(1,6*shots_plot);
POS_y = zeros(1,6*shots_plot);

BETAL = data.raw.images.BETAL;
[~,epics_index,BETAL_index]=intersect(data.raw.scalars.PATT_SYS1_1_PULSEID.UID,BETAL.UID);
USTORO = data.raw.scalars.GADC0_LI20_EX01_CALC_CH2_.dat(epics_index);
DSTORO = data.raw.scalars.GADC0_LI20_EX01_CALC_CH3_.dat(epics_index);
PYRO = data.raw.scalars.BLEN_LI20_3014_BRAW.dat(epics_index);
DIFFTORO = USTORO - DSTORO;
USTORO_goodshots = zeros(6,shots_plot);
DIFFTORO_goodshots = zeros(6,shots_plot);
PYRO_goodshots = zeros(6,shots_plot);
BETAL_charge = zeros(6,shots_plot);


list_shots = {};

list_shots{1} = 2:2:100;
% list_shots{1} = [02, 04, 06, 08, 10, 12, 14, 16, 18, 20];
% list_shots{1} = [22, 24, 26, 28, 30, 32, 34, 36, 38, 40];
% list_shots{1} = [42, 44, 46, 48, 50, 52, 54, 56, 58, 60];
% list_shots{1} = [62, 64, 66, 68, 70, 72, 74, 76, 78, 80];
% list_shots{1} = [82, 84, 86, 88, 90, 92, 94, 96, 98, 100];
%list_shots{1} = [10, 22, 24, 26, 50, 56, 62, 76, 88, 92];

list_shots{2} = 102:2:200;
% list_shots{2} = [102, 104, 106, 108, 110, 112, 114, 116, 118, 120];
% list_shots{2} = [122, 124, 126, 128, 130, 132, 134, 136, 138, 140];
% list_shots{2} = [142, 144, 146, 148, 150, 152, 154, 156, 158, 160];
% list_shots{2} = [162, 164, 166, 168, 170, 172, 174, 176, 178, 180];
% list_shots{2} = [182, 184, 186, 188, 190, 192, 194, 196, 198, 200];
%list_shots{2} = [108, 110, 122, 124, 128, 130, 134, 138, 156, 160];

list_shots{3} = 202:2:300;
%list_shots{3} = [202, 204, 206, 208, 210, 212, 214, 216, 218, 220];
% list_shots{3} = [222, 224, 226, 228, 230, 232, 234, 236, 238, 240];
% list_shots{3} = [242, 244, 246, 248, 250, 252, 254, 256, 258, 260];
% list_shots{3} = [262, 264, 266, 268, 270, 272, 274, 276, 278, 280];
% list_shots{3} = [282, 284, 286, 288, 290, 292, 294, 296, 298, 300];

list_shots{4} = 302:2:400;
% list_shots{4} = [302, 304, 306, 308, 310, 312, 314, 316, 318, 320];
% list_shots{4} = [322, 324, 326, 328, 330, 332, 334, 336, 338, 340];
% list_shots{4} = [342, 344, 346, 348, 350, 352, 354, 356, 358, 360];
% list_shots{4} = [362, 364, 366, 368, 370, 372, 374, 376, 378, 380];
% list_shots{4} = [382, 384, 386, 388, 390, 392, 394, 396, 398, 400];
%list_shots{4} = [306, 310, 320, 322, 324, 330, 338, 340, 342, 346];

list_shots{5} = 402:2:500;
% list_shots{5} = [402, 404, 406, 408, 410, 412, 414, 416, 418, 420];
% list_shots{5} = [422, 424, 426, 428, 430, 432, 434, 436, 438, 440];
% list_shots{5} = [442, 444, 446, 448, 450, 452, 454, 456, 458, 460];
% list_shots{5} = [462, 464, 466, 468, 470, 472, 474, 476, 478, 480];
% list_shots{5} = [482, 484, 486, 488, 490, 492, 494, 496, 498, 500];
%list_shots{5} = [402, 404, 406, 408, 410, 412, 420, 422, 424, 440];
% Bad Shots: 434, 446, 464, 472, 480, 484, 494, 500

list_shots{6} = 502:2:600;
% list_shots{6} = [502, 504, 506, 508, 510, 512, 514, 516, 518, 520];
% list_shots{6} = [522, 524, 526, 528, 530, 532, 534, 536, 538, 540];
% list_shots{6} = [542, 544, 546, 548, 550, 552, 554, 556, 558, 560];
% list_shots{6} = [562, 564, 566, 568, 570, 572, 574, 576, 578, 580];
% list_shots{6} = [582, 584, 586, 588, 590, 592, 594, 596, 598, 600];
%list_shots{6} = [504, 510, 512, 514, 516, 522, 530, 532, 538, 542];
% Very Bad Shots: 502, 508, 526, 536, 540, 544, 548, 554, 556, 562, 570,
% 574, 580, 582, 584, 592
% Not very intense: 506, 518, 520, 524, 550, 558, 560, 568, 572, 578, 586

% BETAL = data.raw.images.BETAL;
% [~,epics_index,BETAL_index]=intersect(data.raw.scalars.PATT_SYS1_1_PULSEID.UID,BETAL.UID);

W1=zeros(zoom_height,zoom_width*shots_plot);
W2=zeros(zoom_height,zoom_width*shots_plot);
W3=zeros(zoom_height,zoom_width*shots_plot);
W4=zeros(zoom_height,zoom_width*shots_plot);
W5=zeros(zoom_height,zoom_width*shots_plot);
W6=zeros(zoom_height,zoom_width*shots_plot);
% c4=0; %column counter
% c3=0;
% c2=0;
% c1=0;
% c5=0;
% c6=0;

for i = 1:6 %numel(BETAL.dat)
    
    %if (mod(i,100)==0) %100 shots per step
    c=0; %column counter
    for j=1:shots_plot
        
        B5D36_BDES = data.raw.scalars.step_value.dat(list_shots{i}(j));
        image = double(imread([prefix BETAL.dat{list_shots{i}(j)}]));
        BETAL_zoom = image(roi.top:roi.bottom,roi.left:roi.right);
        BETAL_process = image(roi.top:roi.bottom+10,roi.left-10:roi.right+10);
        BETAL_process = medfilt2(BETAL_process,[3 3]);
%         BETAL_process = conv2(BETAL_process,ones(5)/5^2,'same');
        
        %Process
        for k=1:zoom_height
            
            mean_left = 0;
            mean_right = 0;
            mean_width = 5;
            mean_height = 1;
            for m=1:mean_width
                for l=1:mean_height
                    mean_left = mean_left + BETAL_process(k+l-1,m);
                    mean_right = mean_right + BETAL_process(k+l-1,zoom_width-m+1); %pretty sure it's supposed to be zoom_width+20-m+1
                end
            end
            bg = linspace(mean_left/mean_width/mean_height, mean_right/mean_width/mean_height, zoom_width+20);
            BETAL_process(k,:) = BETAL_process(k,:) - bg;
            
        end
        
        BETAL_process = conv2(BETAL_process,ones(3)/3^2,'same');
        BETAL_process = BETAL_process(1:end-10,11:end-10);
        
        %BETAL_process2 = BETAL_process(:,roi_process.left:roi_process.right);
        
        E_AXIS = get_BETAL_axis(B5D36_BDES,pix_inf);
        esub = E_AXIS(roi.top:roi.bottom);
        
        if (B5D36_BDES==1)
            W1(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub1 = esub;
            c=c+1;
        elseif (B5D36_BDES==2)
            W2(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub2 = esub;
            c=c+1;
        elseif (B5D36_BDES==3)
            W3(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub3 = esub;
            c=c+1;
        elseif (B5D36_BDES==4)
            W4(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub4 = esub;
            c=c+1;
        elseif (B5D36_BDES==5)
            W5(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub5 = esub;
            c=c+1;
        elseif (B5D36_BDES==6)
            W6(:,((c*zoom_width+1):((c+1)*zoom_width)))=BETAL_process(:,:);
            esub6 = esub;
            c=c+1;
        end
        
        %Positron beam position
        
        %Before QS
        pos_bef_x = data.raw.scalars.BPMS_LI20_3265_X.dat(list_shots{i}(j));
        pos_bef_y = data.raw.scalars.BPMS_LI20_3265_Y.dat(list_shots{i}(j));
        
        %After QS
        pos_aft_x = data.raw.scalars.BPMS_LI20_3315_X.dat(list_shots{i}(j));
        pos_aft_y = data.raw.scalars.BPMS_LI20_3315_Y.dat(list_shots{i}(j));
        
        disp_x = pos_aft_x - pos_bef_x;
        disp_y = pos_aft_y - pos_bef_y;
        POS_x((i-1)*shots_plot+j) = disp_x;
        POS_y((i-1)*shots_plot+j) = disp_y;
        
        %Plotting
        if plot_process
            
            subplot(1,2,1);
            pcolor(1:(zoom_width),esub,BETAL_process);shading interp;
            colorbar; caxis([0 30]);
            %daspect([1 1 1]); axis ij ;
            xlim([1 zoom_width]);
            ylim([esub(1) esub(end)]);
            %pause(0.1);
            subplot(1,2,2);
            pcolor(BETAL_zoom); shading flat;
            colorbar; caxis([0 300]);
            pause(0.3);
            
        elseif plot_zoom
            
            pcolor(BETAL_zoom); shading flat;
            colorbar; caxis([0 300]);
            daspect([1 1 1]); axis ij ;
            %xlim([1 zoom_width]);
            %ylim([1 zoom_height]);
            pause(0.3);
            
        elseif plot_fullBETAL
            
            figure();
            imagesc(image);
            colorbar; caxis([0 300]);
            pause(0.1);
            
        end %if plot
        
        BETAL_charge(i,j) = sum(sum(BETAL_process));
        USTORO_goodshots(i,j) = USTORO(list_shots{i}(j));
        DIFFTORO_goodshots(i,j) = DIFFTORO(list_shots{i}(j));
        PYRO_goodshots(i,j) = PYRO(list_shots{i}(j));
                
    end %for j
    
%     figure;
%     [valsdiff,indexdiff]=sort(DIFFTORO_goodshots);
%     
%     plot(BETAL_charge,valsdiff,'s');
%     plot(valsdiff,BETAL_charge(indexdiff));
%     xlabel('USTORO - DSTORO'); ylabel('BETAL SUM');
%     title (['B5D36 = ' num2str(i)]);
    
    %     end %if mod i
    
end %for i

%% Comparison BETAL charges and diff TORO
% 
% BETAL_charge=BETAL_charge(DIFFTORO_goodshots>0);
% DIFFTORO_goodshots=DIFFTORO_goodshots(DIFFTORO_goodshots>0);
%plot (abs(BETAL_charge-DIFFTORO_goodshots));
[valsdiff,indexdiff]=sort(DIFFTORO_goodshots);
plot(valsdiff,BETAL_charge(indexdiff));
xlabel('USTORO - DSTORO'); ylabel('BETAL SUM');
% plot(DIFFTORO_goodshots);

%% Average charge / number of counts in the injected electron beam

charge1 = sum(sum(W1))/shots_plot;
charge2 = sum(sum(W2))/shots_plot;
charge3 = sum(sum(W3(esub3>4.8,:)))/shots_plot;
charge4 = sum(sum(W4(esub4>7.1,:)))/shots_plot;
charge6 = sum(sum(W6(esub6>9.6,:)))/shots_plot;

total_charge = charge1+charge2+charge3+charge4+charge5+charge6;

e_per_count = 1; % Need to calibratre the # of electrons per count on BETAL
total_charge = e_per_count * total_charge;

% figure
% set(gcf,'color','white');
% set(gcf, 'PaperPosition', [0.25, 2.5, 10, 4]);
% colormap(cmap.wbgyr);
% 
% pcolor(1:zoom_width*shots_plot,esub1,W1); shading flat; colorbar; caxis([0 30]);
% 
% set(gca, 'Box', 'Off')
% set(gca, 'fontsize', 18)
% set(gca, 'XTickLabel', [])
% ylabel('E (GeV)', 'fontsize', 20)
% title('Dipole magnetic field = 20 mT', 'fontsize', 20)

%% Waterfall

figure(1);
cmap = custom_cmap();

colormap(cmap.wbgyr);
subplot(3,2,1);
pcolor(1:zoom_width*shots_plot,esub1,W1); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 1 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,2);
pcolor(1:zoom_width*shots_plot,esub2,W2); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 2 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,3);
pcolor(1:zoom_width*shots_plot,esub3,W3); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 3 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,4); 
pcolor(1:zoom_width*shots_plot,esub4,W4); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 4 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,5);
pcolor(1:zoom_width*shots_plot,esub5,W5); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 5 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,6);
pcolor(1:zoom_width*shots_plot,esub6,W6); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 6 GeV'; ylabel 'E (GeV)' ;

subtitle('10-GeV-class electron acceleration in a positron wake');

%% Plot with y axis in mm

%Get the axis in mm - 35mm = 309px in BetaL 
y_top = (pix_inf-roi.top+1)*35/309;
y_bottom = (pix_inf-roi.bottom+1)*35/309;
y_axis=linspace(y_top,y_bottom,zoom_height);
esub = y_axis;

figure(2);
cmap = custom_cmap();

%subtitle('10-GeV-class electron acceleration in a positron wake');
colormap(cmap.wbgyr);
% subplot(3,2,1);
% pcolor(1:zoom_width*shots_plot,esub,W1); shading flat; colorbar; caxis([0 50]);
% title 'B5D36 = 1 GeV'; ylabel 'E (GeV)' ;
% subplot(3,2,2);
% pcolor(1:zoom_width*shots_plot,esub,W2); shading flat; colorbar; caxis([0 50]);
% title 'B5D36 = 2 GeV'; ylabel 'E (GeV)' ;
% subplot(3,2,3);
% pcolor(1:zoom_width*shots_plot,esub,W3); shading flat; colorbar; caxis([0 50]);
% title 'B5D36 = 3 GeV'; ylabel 'E (GeV)' ;
% subplot(3,2,4); 
% pcolor(1:zoom_width*shots_plot,esub,W4); shading flat; colorbar; caxis([0 50]);
% title 'B5D36 = 4 GeV'; ylabel 'E (GeV)' ;
subplot(1,2,1);
pcolor(1:zoom_width*shots_plot,esub,W5); shading flat; colorbar; caxis([0 25]);
title 'B5D36 = 5 GeV'; ylabel 'E (GeV)' ;
subplot(1,2,2);
pcolor(1:zoom_width*shots_plot,esub,W6); shading flat; colorbar; caxis([0 25]);
title 'B5D36 = 6 GeV'; ylabel 'E (GeV)' ;


%% Plot Positron beam position

figure(3);
plot(POS_x, '*');
title('X Dispersion by QS','fontsize',16);
ylabel ('X (after OS) - X (before QS)','fontsize',13);
xlabel ('QS from 2.35 GeV to 20.35 GeV','fontsize',13)
Xtick_position = zeros(1,7);
Xtick_position(1) = 1;
Xtick_position(7) = 301;
for i=1:5
    Xtick_position(i+1) = 50*i +1;
    line([Xtick_position(i+1) Xtick_position(i+1)],[-1.5 1.5],'color','k','linestyle','--');
end
set(gca, 'XTick', (Xtick_position(2)+Xtick_position(1))/2);
set(gca,'XTicklabel',data.raw.metadata.param.dat{1}.PV_scan_list(1));

figure(4);
plot(POS_y, '*');
title('Y Dispersion by QS','fontsize',16);
ylabel ('Y (after OS) - Y (before QS)','fontsize',13);
xlabel ('QS from 2.35 GeV to 20.35 GeV','fontsize',13)
Xtick_position = zeros(1,8);
Xtick_position(1) = 1;
Xtick_position(7) = 71;
for i=1:5
    Xtick_position(i+1) = 50*i +1;
    line([Xtick_position(i+1) Xtick_position(i+1)],[-1.5 1.5],'color','k','linestyle','--');
end
set(gca, 'XTick', (Xtick_position(2)+Xtick_position(1))/2);
set(gca,'XTicklabel',data.raw.metadata.param.dat{1}.PV_scan_list(1));

%%

figure(1);
cmap = custom_cmap();

figure
set(gcf,'color','white');
set(gcf, 'PaperPosition', [0.25, 2.5, 10*2.54, 4*2.54]);
colormap(cmap.bg);

pcolor(1:zoom_width*shots_plot,esub6,W6); shading flat; colorbar; caxis([0 30]);

set(gca, 'Box', 'Off')
set(gca, 'fontsize', 18)
set(gca, 'XTickLabel', [])
ylabel('E (GeV)', 'fontsize', 20)
title('Dipole magnetic field = 120 mT', 'fontsize', 20)

%%
saveas(gcf, 'Dipole_at_6_GeV_matrix.png');











