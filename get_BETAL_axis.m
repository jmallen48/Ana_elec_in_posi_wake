function E_AXIS = get_BETAL_axis(B5D36_BDES, pix_inf)

% pix_inf = 384.5;
RES = 1E-6 * 35E3/309; %m/pixel (hauteur de la coupure (height of the cut?) 35mm = 309 pix )
BETAL_Z = 2015.32;
THETA0 = 6e-3;
E0 = B5D36_BDES;
z_B5D36 = 2005.65085; % middle of magnet
L = BETAL_Z - z_B5D36;
D0 = THETA0 * L;

pix0 = pix_inf - D0/RES;
E_AXIS = E0 ./ (1-((1:734)-pix0)*RES/D0);
