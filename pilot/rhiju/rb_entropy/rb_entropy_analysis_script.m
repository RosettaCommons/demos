set(gcf, 'PaperPositionMode','auto','color','white','Position',[50 50 800 800]);
clf;

xyz  = load( 'xyz.txt' );
Ixyz = eig( xyz' * xyz ); % moments of inertia
delx = 0.1; N = size( xyz, 1 );
x = [0:delx:10];

PLOT_ALL = 1;
if PLOT_ALL
% Results of 2D simulation
subplot(2,2,1);
rmsd = load( '2D_rotations/rmsd.txt' ); Nvals = length( rmsd );
h = hist( rmsd, x ) / Nvals / delx; % dp/dr.
loglog( x, h, 'o' ); hold on; 
% low RMSD approximation (constant)
pred_curve1 = ones(1,length(x)) * sqrt(N) /pi / sqrt( sum( Ixyz ) );
% full expression
pred_curve2 = (N*x)/(pi*sum(Ixyz) ) ./ sqrt( 1 - (1 - N*x.^2/(2*sum(Ixyz))).^2);
plot( x, pred_curve1, 'k' ); 
plot( x, pred_curve2, 'r' ); hold off
ylim([0.05 2]);
set(gca,'fontweight','bold','xgrid','on','ygrid','on');
xlabel( 'RMSD (Angstroms)');
ylabel( 'Probability distribution (dP/dRMSD)');
title( '2D rotations' ); legend( 'simulated', 'low RMSD analytical', 'analytical',2 )


% Results of 3D simulation
subplot(2,2,3);
rmsd = load( '3D_rotations/rmsd.txt' ); Nvals = length( rmsd );
h = hist( rmsd, x ) / Nvals / delx; % dp/dr.
loglog( x, h, 'o' ); hold on; 
% 2D example (for testing)
% low RMSD approximation (constant)
pred_curve1 = N^(3/2) * x.^2 /(2*pi) / sqrt(Ixyz(1)+Ixyz(2)) / sqrt(Ixyz(1)+Ixyz(3)) / sqrt(Ixyz(2)+Ixyz(3));
plot( x, pred_curve1, 'k' ); hold off
ylim([1e-6 2]);
set(gca,'fontweight','bold','xgrid','on','ygrid','on');
xlabel( 'RMSD (Angstroms)');
ylabel( 'Probability distribution (dP/dRMSD)');
title( '3D rotations' ); legend( 'simulated', 'low RMSD analytical',2 );
end

% Results of 6D simulation
subplot(1,2,2);
molar = (1e27/6.022e23);

rmsd1 = load( '6D_1AngstromMax/rmsd.txt' ); Nvals = length( rmsd1 );
box_size = 1.0; C = molar*1.0/((2.0*box_size)^3);
h1 = hist( rmsd1, x ) / Nvals / delx/ C; % dp/dr.
loglog( x, h1, 'o' ); hold on; 

rmsd2 = load( '6D_2AngstromMax/rmsd.txt' ); Nvals = length( rmsd2 );
box_size = 2.0; C = molar*1.0/((2.0*box_size)^3);
h2 = hist( rmsd2, x ) / Nvals / delx / C; % dp/dr.
loglog( x, h2, 'o','color',[0 0.7 0] ); hold on; 

rmsd3 = load( '6D_3AngstromMax/rmsd.txt' ); Nvals = length( rmsd3 );
box_size = 3.0; C = molar*1.0/((2.0*box_size)^3);
h3 = hist( rmsd3, x ) / Nvals / delx / C; % dp/dr.
loglog( x, h3, 'o','color',[0.5 0 0] ); hold on; 


% 2D example (for testing)
% low RMSD approximation (constant)
pred_curve1 = (1.0/molar) * N^(3/2) * x.^5 * pi /8 / sqrt(Ixyz(1)+Ixyz(2)) / sqrt(Ixyz(1)+Ixyz(3)) / sqrt(Ixyz(2)+Ixyz(3));
plot( x, pred_curve1, 'k' ); hold off
ylim([1e-8 2]);
set(gca,'fontweight','bold','xgrid','on','ygrid','on');
xlabel( 'RMSD (Angstroms)');
ylabel( 'Probability distribution per Molar (dP/dRMSD/C)');
title( '6D rotations' ); legend( 'simulated 1.0A', 'simulated 2.0A', 'low RMSD analytical',2 );