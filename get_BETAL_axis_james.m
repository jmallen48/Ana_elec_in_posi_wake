%% %Visualisation of 20.35GeV e- on CmosFAR
   %y(20.35) depends on the QS value and on the B5D36 value
   %This function plots y(E,QS,B5D36).
   %QS Scan : plots y(E=20.35GeV,QS,B5D36=3GeV)
   %Dipole Scan : plots y(E=20.35GeV,QS=-10GeV,B5D36)

clear all;
QS_Scan = 1; % 0 : Dipole Scan

%%

header = '/Users/margauxschmeltz/Desktop/SLAC/Data Analysis';
prefix = '/Volumes/PWFA_4big';
expt = 'E200';
data_path = ['/nas/nas-li20-pm00/' expt '/2014/'];
day = '20140629';

%E200_13540 QS Scan
%E200_13543 Dipole Scan

if QS_Scan
    data_set = 'E200_13540';
else %Dipole_Scan
    data_set = 'E200_13543';
end

%%
user = 'e157user';

%% Load data

if ~exist(prefix)
    system(['mkdir ' prefix]);
    system(['/usr/local/bin/sshfs ' user '@quickpicmac3.slac.stanford.edu:/Volumes/PWFA_4big/' ' ' prefix]);
    pause(2);
end
load([prefix data_path day '/' data_set '/' data_set '.mat']);
cmap = custom_cmap();


%% Select UID

CMOS_FAR=data.raw.images.CMOS_FAR;
[~,epics_index,CMOS_FAR_index]=intersect(data.raw.scalars.PATT_SYS1_1_PULSEID.UID,CMOS_FAR.UID);

%% CMOS FAR Zoom (in pixels)

if QS_Scan
    roi.top = 2250;
    roi.bottom = 2559;
    roi.left = 170;
    roi.right = 580;
    
else %Dipole_Scan
    roi.top = 1750;
    roi.bottom = 2559;
    roi.left = 200;
    roi.right = 500;
end

zoom_height = roi.bottom - roi.top + 1;
zoom_width = roi.right - roi.left + 1;

%% Find the position of E0 in CMOS FAR

calib=data.raw.images.CMOS_FAR.RESOLUTION(1); %um/pix
list_shots = {};

if QS_Scan
    
    shots_plot = 50;
    nsteps=7;
    
    %13540 Step 1 : laser on for odd shots (i)
    list_shots{1} = 2:2:100;
    %list_shots{1} = [8, 10, 18, 30, 32, 48, 56, 80, 88, 92];
    
    %13540 Step 2 : laser on for odd shots (i)
    list_shots{2} = 102:2:200;
    %list_shots{2} = [112, 114, 116, 122, 140, 148, 150, 156, 168, 170];
    
    %13540 Step 3 : laser on for odd shots (i)
    list_shots{3} = 202:2:300;
    %list_shots{3} = [202, 210, 218, 248, 254, 284, 292, 298, 222, 258];
    
    %13540 Step 4 : laser on for even shots (p)
    list_shots{4} = 301:2:399;
    %list_shots{4} = [303, 317, 357, 373, 393, 309, 395, 397, 301, 323];
    
    %13540 Step 5 : laser on for even shots (p)
    list_shots{5} = 401:2:499;
    %list_shots{5} = [401, 403, 405, 407, 409, 411, 413, 415, 417, 419];
    
    %13540 Step 6 : laser on for odd shots (i)
    list_shots{6} = 502:2:600;
    %list_shots{6} = [502, 504, 506, 508, 510, 512, 514, 516, 518, 520];
    
    %13540 Step 7 : laser on for even shots (p)
    list_shots{7} = 601:2:699;
    %list_shots{7} = [601, 603, 605, 613, 615, 619, 621, 623, 625, 627];

else %Dipole_Scan
    
    shots_plot = 50;
    nsteps=6;
    
    %laser on for odd shots (i)
    list_shots{1} = 2:2:100;
    list_shots{2} = 102:2:200;
    list_shots{3} = 202:2:300;
    list_shots{4} = 302:2:400;
    list_shots{5} = 402:2:500;
    list_shots{6} = 502:2:600;
    
end

E0_position=ones(1,nsteps);
B5D36=ones(1,nsteps);
QS=ones(1,nsteps);

for i=1:nsteps
    
    if QS_Scan
        QS(i) = data.raw.scalars.step_value.dat(list_shots{i}(1));
    else
        B5D36(i) = data.raw.scalars.step_value.dat(list_shots{i}(1));
    end
    max_indices=ones(1,shots_plot);
    
    for j=1:shots_plot
    
        image = double(imread([prefix CMOS_FAR.dat{list_shots{i}(j)}]));
        image = rot90(image);
        img_zoom = image(roi.top:roi.bottom,roi.left:roi.right);
%         figure(1);
%         imagesc(img_zoom);
%         pause(0.2);
        y_prof = mean(img_zoom,2);
        [~,max_index] = max(y_prof);
        max_indices(j)=max_index;
        
    end
    
    %%From index in pixels in the zoom to index in um
    pix_zoom = mean(max_indices);
    pix_total = pix_zoom+roi.top-1;
    E0_position(i)=pix_total;%(390-pix_total)*calib;

end

%%
figure(2)

if QS_Scan
    plot(QS,E0_position,'s');xlabel ('QS value');ylabel ('y (pix)');
else
    plot(B5D36,E0_position,'s');xlabel ('B5D36 value');ylabel ('y (pix)');
end

%%
F=fit(QS',E0_position','poly2');

plot(QS,E0_position,'s');
hold on;
plot(F);
xlabel ('QS value');ylabel ('y (pix)');


%% Matlab solver

syms a b Emax

eta=60.9;

%Equation 1
QS=5.35;
B5D36=3;
eq1 = (a*(QS/Emax)^2+(eta*B5D36+b*QS)/Emax == 37);

%Equation 2
QS=8.35;
B5D36=3;
eq2 = (a*(QS/Emax)^2+(eta*B5D36+b*QS)/Emax == 20);

%Equation 3
QS=10.35;
B5D36=5;
eq3 = (a*(QS/Emax)^2+(eta*B5D36+b*QS)/Emax == 26);

%Equation 4
QS=10.35;
B5D36=6;
eq4 = (a*(QS/Emax)^2+(eta*B5D36+b*QS)/Emax == 36);

Sa=zeros(1,4);
Sb=zeros(1,4);
SEmax=zeros(1,4);

[Sa(1),Sb(1),SEmax(1)]=solve(eq1,eq2,eq3);
[Sa(2),Sb(2),SEmax(2)]=solve(eq1,eq2,eq4);
[Sa(3),Sb(3),SEmax(3)]=solve(eq1,eq3,eq4);
[Sa(4),Sb(4),SEmax(4)]=solve(eq2,eq3,eq4);
%[Sa,Sb,Sc,SEmax]=solve(eq1,eq2,eq3,eq4);
%eval([Sa,Sb,Sc,SEmax])

figure(3); hold on;
%plot(Sa,'s','color','k'); hold on; plot(Sb,'*','color','b') ; plot(SEmax,'.','color','r');

plot(mean(Sa)*ones(1,4),'color','k');
plot(mean(Sb)*ones(1,4),'color','b');
plot(mean(SEmax)*ones(1,4),'color','r');
legend(num2str(mean(Sa)),num2str(mean(Sb)),num2str(mean(SEmax)));

%% Approx : delta1=0

syms b Emax

eta=60.9;

%Equation 1
QS=5.35;
B5D36=3;
eq1 = ((eta*B5D36+b*QS)/Emax == 37);

%Equation 2
QS=8.35;
B5D36=3;
eq2 = ((eta*B5D36+b*QS)/Emax == 20);

%Equation 3
QS=10.35;
B5D36=5;
eq3 = ((eta*B5D36+b*QS)/Emax == 26);

%Equation 4
QS=10.35;
B5D36=6;
eq4 = ((eta*B5D36+b*QS)/Emax == 36);

Sb=zeros(1,6);
SEmax=zeros(1,6);

[Sb(1),SEmax(1)]=solve(eq1,eq2);
[Sb(2),SEmax(2)]=solve(eq1,eq3);
[Sb(3),SEmax(3)]=solve(eq1,eq4);
[Sb(4),SEmax(4)]=solve(eq2,eq3);
[Sb(5),SEmax(5)]=solve(eq2,eq4);
[Sb(6),SEmax(6)]=solve(eq3,eq4);

figure(3); hold on;
plot(Sb,'*','color','b') ; plot(SEmax,'.','color','r');
plot(mean(Sb)*ones(1,6),'color','b');
plot(mean(SEmax)*ones(1,6),'color','r');

%%
D2=5;
D3=12.12;
f10=2.5936; %3.19;
f20=-4.0416; %-4.33;
eta_dipole_0 = -55e-3;

    
data=1E-3*[37 20 26 36];
data=1E-3*[37 20 26 36];
weight = [1 1 1 1];
x=[5.35 8.53 10.35 10.35];
y=[3 3 5 6];
coeff0(1)=1;
coeff0(2)=1;
coeff0(3)=1;
   
F = quadratic_fit(x,y,data,weight,coeff0);
plot(data,'s'); hold on;
plot(F(1)*x.^2+F(2)*x-F(3)*y,'*');

Emax=eta_dipole_0/F(3)
a=F(1)*Emax*Emax;
b=F(2)*Emax;

delta1 = -a*f10*f20/(D2*D3)
delta2 = f20/D3*(b-delta1*(D2+D3)/f10)

QS=10.35;
B5D36=6;
y = delta1*(D2+D3)/f10*QS/Emax - delta1*D2*D3/(f10*f20)*(QS/Emax)^2 + (delta2*D3/f20)*QS/Emax - eta_dipole_0*B5D36/Emax
hold off


