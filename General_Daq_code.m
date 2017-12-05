%% Script for DAQ 1D Scan

clear variables;
cmap = custom_cmap;

k = 543; % Data number (last three digits)
j = 'USOTR'; % Choose camera (or leave empty '')
%% Load data
% header = '/Users/jmallen/FACET/PWFA_5big'; % If using mount
header = '/Users/jmallen/FACET/Electrons_in_Positron_Wake/Data'; % if using data on this machine
nas = '/nas/nas-li20-pm00';
experiment = '/E200';
year = '/2014';
date = '/20140629';
data_set = ['E200_13' num2str(k)];
path = [header nas experiment year date '/' data_set '.mat'];
load(path);

%% Select cameras, get common UIDs, and get scalars
epics_ind = [];
cam_ind = [];
cams = [];
scalars = [];
[epics_ind, cams, cam_ind] = select_cams(data);
scalars = getScalars(data,epics_ind);

%% Store images in cell array
if isempty(j)
    disp(cams);
    prompt = 'Display images for which camera?';
    j = input(prompt,'s'); % Asks user which camera to use, requires text input matching one of the displayed cameras
end
i = find(strcmp(j,cams)); % Finds matching index to user input in 'cams'
fldnm = cams{i};
cam_struct = data.raw.images.(fldnm); % Stores all images for desired camera
% roi = select_ROI(im); % Select ROI, click top left, then bottom right of desired ROI

if cam_struct.X_ORIENT(cam_ind(1,i)) % Check orientation of x and y axes
    x_orient = 1;
else
    x_orient = 0;
end

if cam_struct.Y_ORIENT(cam_ind(1,i))
    y_orient = 1;
else
    y_orient = 0;
end

scan_list = data.raw.metadata.param.dat{1}.PV_scan_list; % list of scan values
step_value = data.raw.scalars.step_value.dat(epics_ind); % vector of scan value for each shot
xnp = cam_struct.ROI_XNP(1);
ynp = cam_struct.ROI_YNP(1);
step_ind = cell(1,numel(scan_list)); % create empty cell
im_array = cell(1,numel(scan_list)); % create empty cell

h = waitbar(0,'Please wait...');
for w = 1:numel(scan_list)
    step_ind{w} = find(step_value==scan_list(w)); % fill cell with separate vectors of indices for each step
    im_array{w} = zeros(numel(step_ind{w}),ynp,xnp);
    for u = 1:numel(step_ind{w})
        im = imread([header cam_struct.dat{step_ind{w}(u)}]);
        if y_orient
            im = flipud(im);
        end
        if x_orient
            im = fliplr(im);
        end
        bg_struct = load([header cam_struct.background_dat{1}]);
        bg = bg_struct.img;
        im = im-bg;
        im_array{w}(u,:,:) = im; % fill cell with image matrix
        %E224_2 = data(data.raw.images.E224_2;
        %cam_loc = strcmp('E224_2',cams{i});
        %do = [];
        %E224_2_ANA(i) = BAISG_2016(E224_2,1,E224_2_roi,header,do,cam_ind{i}(:,cam_loc));
        %step_val_axis = [step_val_axis scalars{i}.step_value];
    end
    waitbar(w/numel(scan_list));
end
close(h);

%% Calculation test
avg_px = cell(1,numel(scan_list)); % create empty cell
for m = 1:numel(scan_list)
    avg_px{m} = zeros(numel(step_ind{m}));
    for n = 1:numel(step_ind{m})
        avg_px{m}(n) = mean(mean(im_array{m}(n,:,:))); % fill cell with calculated value
    end
end

%% Plot calculation test
figure(102);
hold on;
for b = 1:numel(scan_list) % Plot each vector in the same figure
    plot(avg_px{b});
end
hold off;

%% Nfilter test single image
% roi = select_ROI(cam_struct);
roi.top = 1;
roi.bottom = 250;
roi.left = 640;
roi.right = 720;
test_im = squeeze(im_array{1}(30,roi.top:roi.bottom,roi.left:roi.right));

post_im = nlfilter(test_im,[21 21],@mean_dif);

figure(1);
imagesc(test_im);
title('Before');
colormap(cmap.wbgyr);
caxis([0 750]);
colorbar;

figure(2);
imagesc(post_im);
title('After');
colormap(cmap.wbgyr);
caxis([0 750]);
colorbar;

%% Gradient test single image
roi.top = 1;
roi.bottom = 250;
roi.left = 640;
roi.right = 720;
test_im = squeeze(im_array{1}(30,roi.top:roi.bottom,roi.left:roi.right));

%test_im = squeeze(im_array{1}(30,:,:));

[gmag,gdir] = imgradientxy(test_im);

figure(1);
imagesc(test_im);
title('Before');
colormap(cmap.wbgyr);
caxis([0 750]);
colorbar;

figure(2);
imagesc(gdir);
title('After');
colormap(cmap.wbgyr);
caxis([-200 200]);
colorbar;

%% Gradient test loop
for t = 5:numel(im_array)
    for r = 1:numel(im_array{t}(:,1,1))
%         roi.top = 1;
%         roi.bottom = 250;
%         roi.left = 640;
%         roi.right = 720;
%         test_im = squeeze(im_array{t}(r,roi.top:roi.bottom,roi.left:roi.right));
        
        test_im = squeeze(im_array{t}(r,:,:));
        
        [gx,gy] = imgradientxy(test_im);
        
        subplot(1,2,1);
        imagesc(test_im);
        title('Image');
        colormap(cmap.wbgyr);
        caxis([0 750]);
        colorbar;
        
        subplot(1,2,2);
        imagesc(gx);
        title('GX');
        colormap(cmap.wbgyr);
        caxis([0 500]);
        colorbar;
        pause(0.2);
    end
end

%% Edge test loop
for t = 5:numel(im_array)
    for r = 1:numel(im_array{t}(:,1,1))
%         roi.top = 1;
%         roi.bottom = 250;
%         roi.left = 640;
%         roi.right = 720;
%         test_im = squeeze(im_array{t}(r,roi.top:roi.bottom,roi.left:roi.right));
        
        test_im = squeeze(im_array{t}(r,:,:));
        
        BW = edge(test_im);
        
        subplot(1,2,1);
        imagesc(test_im);
        title('Image');
        colormap(cmap.wbgyr);
        caxis([0 750]);
        colorbar;
        
        subplot(1,2,2);
        imagesc(BW);
        title('Edge');
        colormap(cmap.wbgyr);
        colorbar;
        pause(0.5);
    end
end
%% Display images
for t = 1:numel(im_array)
    for r = 1:numel(im_array{t}(:,1,1))
        imagesc(squeeze(im_array{t}(r,:,:)));
        colormap(cmap.wbgyr);
        caxis([0 1500]);
        colorbar;
        title(sprintf([func2str(data.raw.metadata.param.dat{1}.fcnHandle) ' Step value: ' ...
            num2str(data.raw.metadata.param.dat{1}.PV_scan_list(t)) ...
            '; Shot number: ' num2str(r)]),'Interpreter','none');
%         title(sprintf([func2str(data.raw.metadata.param.dat{1}.fcnHandle2) ' Step value: ' ...
%             num2str(data.raw.metadata.param.dat{1}.PV_scan_list2(ceil(t/6))) '; ' ...
%             func2str(data.raw.metadata.param.dat{1}.fcnHandle) ' Step value: ' ...
%             num2str(data.raw.metadata.param.dat{1}.PV_scan_list1(ceil(t-6*floor((t-1)/6)))) ...
%             '; Shot number: ' num2str(r)]),'Interpreter','none');
        pause(0.01);
    end
end

%% Plot charge at each toroid
e_charge = 1.60217662E-19;
ch0 = data.raw.scalars.GADC0_LI20_EX01_CALC_CH0_.dat(epics_ind)*e_charge;
ch2 = data.raw.scalars.GADC0_LI20_EX01_CALC_CH2_.dat(epics_ind)*e_charge;
ch3 = data.raw.scalars.GADC0_LI20_EX01_CALC_CH3_.dat(epics_ind)*e_charge;
hold on;
plot(ch0,'k','linewidth',1);
plot(ch2,'r','linewidth',1);
plot(ch3,'g','linewidth',1);
legend({'ch0','ch2','ch3'},'FontSize',20);
xlabel('Shot Number','FontSize',20);
ylabel('Charge [C]','FontSize',20);
set(gca,'FontSize',18);
title([data_set ' Charge per Shot for Each Toroid'],'Interpreter','none');
hold off;
%% In file functions

function f = mean_dif(x)
    mid_ind = floor((size(x)+1)/2);
    z = mean(x(:))-x(mid_ind(1),mid_ind(2));
    if z <= 0
        f = 0;
    else
        f = z;
    end
end

    