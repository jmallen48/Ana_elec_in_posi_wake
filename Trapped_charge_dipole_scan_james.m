clear variables;

cmap = custom_cmap; % Custom color map for plots

data_digits = 543; % Data number (last three digits)

% Load data
% prefix = '/Users/jmallen/Work/FACET/PWFA_5big'; % If using mount
prefix = '/Users/jmallen/FACET/Electrons_in_Positron_Wake/Data'; % if using data on this machine
nas = '/nas/nas-li20-pm00';
experiment = '/E200';
year = '/2014';
date = '/20140629';
data_set = ['E200_13' num2str(data_digits)];
data_set_str = ['E200\_13' num2str(data_digits)];
path = [prefix nas experiment year date '/' data_set '.mat'];
load(path);

%% Common UIDs
epics_ind = [];
cam_ind = [];
cams = [];
scalars = [];
[epics_ind, cams, cam_ind] = select_cams(data);
scalars = getScalars(data,epics_ind);

%% Find Laser On/Off
cam_field = 'E224_Probe';
cam_struct = data.raw.images.(cam_field);
sum_px = zeros(1,numel(cam_struct.dat)); % create empty vector
% laser_roi = select_ROI(cam_struct);
laser_roi.top = 1;
laser_roi.bottom = 734;
laser_roi.left = 1;
laser_roi.right = 1292;
h = waitbar(0,'Please wait...','Name','Finding Laser On/Off',... % make waitbar
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
for w = 1:numel(cam_struct.dat)
    im = imread([prefix cam_struct.dat{w}]);
    if cam_struct.X_ORIENT(1)  % Check orientation of images
        im = fliplr(im);
    end
    if cam_struct.Y_ORIENT(1)
        im = flipud(im);
    end
    bg_struct = load([prefix cam_struct.background_dat{1}]);
    bg = bg_struct.img;
    im = im-bg;
    sum_px(w) = sum(sum(im(laser_roi.top:laser_roi.bottom,... % Sum total pixels
        laser_roi.left:laser_roi.right)));
    if getappdata(h,'canceling') % Check if cancel button clicked on
        delete(h);
        return
    end
    waitbar(w/numel(cam_struct.dat));
end
delete(h);

laser_status_raw = nlfilter(sum_px,[1 3],@check_avg); % Create vector with 1's where above local avg and 0's where below
cam_num = strcmp(cam_field,cams);
laser_status = laser_status_raw(cam_ind(:,cam_num));

%% BetaL ROI

% Create ROI using pixel values
roi.top = 1;
roi.bottom = 170;
roi.left = 515;
roi.right = 685;

% Create ROI using function
% roi = select_ROI(BETAL);

zoom_height = roi.bottom - roi.top + 1;
zoom_width = roi.right - roi.left + 1;

% Not really used...
% roi_process.top = 1;
% roi_process.bottom = 200;
roi_process.left = 90;
roi_process.right = 140;
process_width = roi_process.right - roi_process.left + 1;

pix_inf = 384.5; % Infinite energy axis pixel

%% BetaL images

% Choose what to plot. To plot set == 1, to skip set == 0.
plot_zoom = 0;
plot_process = 0;
plot_fullBETAL = 0;

% Choose laser on or off data
laser_on_off = 1    ;

if laser_on_off == 1
    laser = find(laser_status);  % define laser on indices. laser status already uses common UIDs
else
    laser = find(laser_status == 0); % define laser off indices. laser status already uses common UIDs
end

scan_list = data.raw.metadata.param.dat{1}.PV_scan_list; % list of scan values
step_value = data.raw.scalars.step_value.dat(epics_ind);
step_value = step_value(laser); % vector of scan value for each shot, only laser on
nstep = numel(scan_list);
step_ind = cell(1,numel(scan_list)); % create empty cell
im_array = cell(1,numel(scan_list)); % create empty cell

% Create empty cells for position difference before/after plasma
POS_x = cell(1,nstep);
POS_y = cell(1,nstep);

% Get BetaL images for either laser on or off
cam_field = 'BETAL';
cam_struct = data.raw.images.(cam_field);
cam_num = strcmp(cam_field,cams);
BETAL_on_off = cam_struct.dat(cam_ind(:,cam_num));
BETAL = BETAL_on_off(laser);
bg_struct = load([prefix cam_struct.background_dat{1}]);
bg_image = im2double(bg_struct.img);
x_orient = 0;
y_orient = 0;
if cam_struct.X_ORIENT(1)  % Check orientation of images
        x_orient = 1;
end
if cam_struct.Y_ORIENT(1)
        y_orient = 1;
end
     
% Get E224_Probe images for either laser on or off
cam_field = 'E224_Probe';
cam_struct = data.raw.images.(cam_field);
cam_num = strcmp(cam_field,cams);
E224_on_off = cam_struct.dat(cam_ind(:,cam_num));
E224 = E224_on_off(laser);

% Get scalars for either laser on or off
USTORO_vec = scalars.GADC0_LI20_EX01_CALC_CH2_(laser);
DSTORO_vec = scalars.GADC0_LI20_EX01_CALC_CH3_(laser);
PYRO_vec = scalars.BLEN_LI20_3014_BRAW(laser);
x_3265 = scalars.BPMS_LI20_3265_X(laser);
y_3265 = scalars.BPMS_LI20_3265_Y(laser);
x_3315 = scalars.BPMS_LI20_3315_X(laser);
y_3315 = scalars.BPMS_LI20_3315_Y(laser);
DIFFTORO_vec = USTORO_vec - DSTORO_vec;

% Create empty to be filled
USTORO = cell(1,nstep);
DSTORO = cell(1,nstep);
DIFFTORO = cell(1,nstep);
PYRO = cell(1,nstep);
BETAL_charge = cell(1,nstep);
wf = cell(1,nstep); % waterfall empty cell
wf_sum = cell(1,nstep); % sum image horizontally before putting into waterfall, empty cell
spectrum_median = cell(1,nstep); % mean of each wf_sum cell
spectrum_mean = cell(1,nstep); % mean of each wf_sum cell
esub = cell(1,nstep); % energy axis empty cell

h = waitbar(0,'Please wait...','Name','Processing BetaL Images',... % make waitbar
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

for i = 1:nstep    
    step_ind{i} = find(step_value==scan_list(i)); % fill cell with separate vectors of indices for each step
    
    % Make place holder matrices
    im_array{i} = zeros(zoom_height,zoom_width,numel(step_ind{i}));
    wf{i} = zeros(zoom_height,zoom_width*numel(step_ind{i}));
    wf_sum{i} = zeros(zoom_height,numel(step_ind{i}));
    POS_x{i} = zeros(numel(step_ind{i}));
    POS_y{i} = zeros(numel(step_ind{i}));
    USTORO{i} = zeros(numel(step_ind{i}));
    DSTORO{i} = zeros(numel(step_ind{i}));
    DIFFTORO{i} = zeros(numel(step_ind{i}));
    PYRO{i} = zeros(numel(step_ind{i}));
    BETAL_charge{i} = zeros(numel(step_ind{i}));
    % B5D36_BDES = step_value(step_ind{i}(1));
    B5D36_BDES = 1;
    spectrum_median{i} = zeros(zoom_height);
    spectrum_mean{i} = zeros(zoom_height);
    
    % BetaL energy axis
    E_AXIS = get_BETAL_axis(B5D36_BDES,pix_inf);
    esub{i} = E_AXIS(roi.top:roi.bottom); % energy axis per step for roi
    c=0; %column counter for filling waterfall
    width = 10; % Choose width of outer border and height to process for each loop
    process_height = 1; % Choose width of outer border and height to process for each loop
    height = 10; % Choose additional height for BETAL_process
    
    for j=1:numel(step_ind{i})
        image = double(imread([prefix BETAL{step_ind{i}(j)}]));
        if x_orient  % Check orientation of images
            image = fliplr(image);
        end
        if y_orient
            image = flipud(image);
        end
        image = image-bg_image;
        if x_orient  % Check orientation of images
            image = fliplr(image);
        end
        if y_orient
            image = flipud(image);
        end
        BETAL_zoom = image(roi.top:roi.bottom,roi.left:roi.right);
        BETAL_process = image(roi.top:roi.bottom+height,roi.left-width:roi.right+width); % Outer border of roi used to process
        BETAL_process = medfilt2(BETAL_process,[3 3]); % Median filter
        
        % Process BetaL images
        for k=1:zoom_height+10-(process_height-1)
            
            mean_left = 0;
            mean_right = 0;
            for m=1:width
                for l=1:process_height
                    mean_left = mean_left + BETAL_process(k+l-1,m);
                    mean_right = mean_right + BETAL_process(k+l-1,zoom_width+20-m+1);
                end
            end
            % background row using outer border of roi:
            bg = linspace(mean_left/width/process_height, mean_right/width/process_height, zoom_width+width*2); 
            BETAL_process(k,:) = BETAL_process(k,:) - bg;
            
            if getappdata(h,'canceling') % Check if cancel button clicked on
                delete(h);
                return
            end
            
        end
        
        % BETAL_process = conv2(BETAL_process,ones(3)/3^2,'same'); % What's this for? 2D convolution
        BETAL_process = BETAL_process(1:end-height,width+1:end-width);
        % BETAL_process(BETAL_process<=20) = 0;
        
        %BETAL_process2 = BETAL_process(:,roi_process.left:roi_process.right);
        
        % Fill cells and increase column counter for waterfall cell
        im_array{i}(:,:,j) = BETAL_zoom;
        wf{i}(:,((c*zoom_width+1):((c+1)*zoom_width))) = BETAL_process;
        wf_sum{i}(:,j) = sum(BETAL_process,2);
        c=c+1;
        
        %Positron beam position
        
        %Before QS
        pos_bef_x = x_3265(step_ind{i}(j));
        pos_bef_y = y_3265(step_ind{i}(j));
        
        %After QS
        pos_aft_x = x_3315(step_ind{i}(j));
        pos_aft_y = y_3315(step_ind{i}(j));
        
        disp_x = pos_aft_x - pos_bef_x;
        disp_y = pos_aft_y - pos_bef_y;
        POS_x{i}(j) = disp_x;
        POS_y{i}(j) = disp_y;
        
        BETAL_charge{i}(j) = sum(sum(BETAL_process));
        USTORO{i}(j) = USTORO_vec(step_ind{i}(j));
        DIFFTORO{i}(j) = DIFFTORO_vec(step_ind{i}(j));
        PYRO{i}(j) = PYRO_vec(step_ind{i}(j));
        
        E224_im = imread([prefix E224{step_ind{i}(j)}]);
        
        %Plotting
        if plot_process
            
            subplot(1,3,1);
            pcolor(BETAL_zoom); shading flat;
            colorbar; caxis([0 300]);
            title('Raw Image');
            
            subplot(1,3,2);
            pcolor(1:(zoom_width),esub{i},BETAL_process);shading interp;
            colorbar; caxis([0 30]);
            %daspect([1 1 1]); axis ij ;
            xlim([1 zoom_width]);
            ylim([esub{i}(1) esub{i}(end)]);
            title('BG subtracted');
            ylabel('Energy (GeV)');

            subplot(1,3,3);
            imagesc(E224_im);
            if laser_on_off == 1
                title('Laser Visibility: On');
            else
                title('Laser Visibility: Off');
            end
            set(findobj(gcf,'type','axes'),'FontSize',14);
            
            pause(0.1);
            
            if getappdata(h,'canceling') % Check if cancel button clicked on
                delete(h);
                return
            end
            
            
        elseif plot_zoom
            
            pcolor(BETAL_zoom); shading flat;
            colorbar; caxis([0 300]);
            daspect([1 1 1]); axis ij ;
            %xlim([1 zoom_width]);
            %ylim([1 zoom_height]);
            pause(0.3);
            
        elseif plot_fullBETAL
            
            imagesc(image);
            colorbar; caxis([0 300]);
            pause(0.1);
            
        end %if plot
    end %for j
    spectrum_median{i} = median(wf_sum{i},2);
    spectrum_mean{i} = mean(wf_sum{i},2);
%     figure;
%     [valsdiff,indexdiff]=sort(DIFFTORO_goodshots);
%     
%     plot(BETAL_charge,valsdiff,'s');
%     plot(valsdiff,BETAL_charge(indexdiff));
%     xlabel('USTORO_vec - DSTORO_vec'); ylabel('BETAL SUM');
%     title (['B5D36 = ' num2str(i)]);
    
    %     end %if mod i
    waitbar(i/nstep);
end %for i
delete(h);
%% Plot Spectrum
for i = 1:6
%     colormap(cmap.wbgyr);
    figure(data_digits*1);
    
    subplot(1,2,1);
    hold on;
    plot(esub{i},spectrum_mean{i});
    xlabel('Energy (GeV)');
    ylabel('Average Charge (arb. units)');
    title([data_set ' Electron Energy Spectrum Mean; First 6 Steps'],'Interpreter','none');
    set(gca,'FontSize',12);
    hold off;
    
    subplot(1,2,2);
    hold on;
    plot(esub{i},spectrum_median{i});
    xlabel('Energy (GeV)');
    ylabel('Average Charge (arb. units)');
    title([data_set ' Electron Energy Spectrum Median; First 6 Steps'],'Interpreter','none');
    set(gca,'FontSize',12);
    hold off;
    
    figure(data_digits*2)
    subplot(2,3,i);
    pcolor(1:numel(step_ind{i}),esub{i},wf_sum{i}); shading flat; colorbar;
end

%% Plot Charge
figure(data_digits*3);
hold on;
for i = 1:numel(step_ind{2})
    plot(esub{2},wf_sum{2}(:,i));
end
hold off;

%% Comparison BETAL charges and diff TORO
% 
% BETAL_charge=BETAL_charge(DIFFTORO_goodshots>0);
% DIFFTORO_goodshots=DIFFTORO_goodshots(DIFFTORO_goodshots>0);
%plot (abs(BETAL_charge-DIFFTORO_goodshots));
[valsdiff,indexdiff]=sort(DIFFTORO);
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
pcolor(1:zoom_width*numel(step_ind{1}),esub{1},wf{1}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 1 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,2);
pcolor(1:zoom_width*numel(step_ind{2}),esub{2},wf{2}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 2 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,3);
pcolor(1:zoom_width*numel(step_ind{3}),esub{3},wf{3}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 3 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,4); 
pcolor(1:zoom_width*numel(step_ind{4}),esub{4},wf{4}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 4 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,5);
pcolor(1:zoom_width*numel(step_ind{5}),esub{5},wf{5}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 5 GeV'; ylabel 'E (GeV)' ;
subplot(3,2,6);
pcolor(1:zoom_width*numel(step_ind{6}),esub{6},wf{6}); shading flat; colorbar; caxis([0 50]);
title 'B5D36 = 6 GeV'; ylabel 'E (GeV)' ;

if laser_on_off == 1
    subtitle(['Laser On ' data_set_str ' Waterfall plots of first 6 steps']);
else
    subtitle(['Laser Off ' data_set_str ' Waterfall plots of first 6 steps']);
end

set(findobj(gcf,'type','axes'),'FontSize',14);

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
% saveas(gcf, 'Dipole_at_6_GeV_matrix.png');

%% Compare laser visibility plot to other PVs

% Laser power
PV = data.raw.scalars.PMTR_LA20_10_PWR.dat(epics_ind);

% Charge after IP
% PV = data.raw.scalars.GADC0_LI20_EX01_CALC_CH2_.dat(epics_ind)-...
%     data.raw.scalars.GADC0_LI20_EX01_CALC_CH3_.dat(epics_ind);
hold on;
plot(PV(1:25));
plot(laser_status(1:25)*2*mean(PV));
legend('Laser Power','Laser Visibility');
title('Laser On/Off vs. Shot Number');
set(gca,'FontSize',18);
hold off;

%% Compare USOTR mean position to BPM 3156
cam_field = 'USOTR';
cam_struct = data.raw.images.(cam_field);
cam_num = find(strcmp(cam_field,cams));
roi = select_ROI(cam_struct);
x_axis = linspace(roi.left,roi.right,roi.right-roi.left+1);
y_axis = linspace(roi.top,roi.bottom,roi.bottom-roi.top+1);
nshots = numel(cam_struct.dat(cam_ind(:,cam_num)));
centroid = zeros(2,nshots); % X coordinate is columns so it's centroid(2,:) and Y is centroid(1,:)

bg_struct = load([prefix cam_struct.background_dat{1}]);
bg = bg_struct.img;

h = waitbar(0,'Please wait...','CreateCancelBtn',... % Make waitbar
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
for i = 1:nshots
    im = imread([prefix cam_struct.dat{cam_ind(i,cam_num)}]);
    im = im - bg;
    im = im(roi.top:roi.bottom,roi.left:roi.right);
%     im = im/sum(im(:));
%     im = im2double(im);
%     [m,n]=size(im);
%     [I,J]=ndgrid(1:m,1:n);
%     centroid(:,i)=[dot(I(:),im(:)),dot(J(:),im(:))];
    proj_x = mean(im,1);
    norm = sum(proj_x);
    centroid(2,i) = sum(proj_x.*x_axis)/norm;
    proj_y = mean(im,2);
    norm = sum(proj_y);
    centroid(1,i) = sum(proj_y.*y_axis)/norm;
    if getappdata(h,'canceling') % Check if cancel button clicked on
        delete(h);
        return
    end
    waitbar(i/nshots);
end
delete(h);

bpm_3156 = scalars.BPMS_LI20_3156_X;

hold on;
scatter(centroid(1,:),bpm_3156);
hold off;
%% Display images
cam_field = 'USOTR';
cam_struct = data.raw.images.(cam_field);
for i = 1:numel(cam_struct.dat)
    im = imread([prefix cam_struct.dat{i}]);
    if cam_struct.X_ORIENT(1)  % Check orientation of images
        im = fliplr(im);
    end
    if cam_struct.Y_ORIENT(1)
        im = flipud(im);
    end
    bg_struct = load([prefix cam_struct.background_dat{1}]);
    bg = bg_struct.img;
    im = im-bg;
    imagesc(im);
%     caxis([0 250]);
    colormap(cmap.wbgyr);
    colorbar;
    title(i);
    pause();
end

%%
function f = check_avg(x)
    if x(2) < mean(x)
        f = 0;
    elseif x(2) > mean(x)
        f = 1;
    end
end










