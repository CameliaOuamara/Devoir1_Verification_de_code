% MATLAB script to launch a fiber structure generation and the corresponding LBM simulation
%
%INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place). 
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
% 
% NX: domain lateral size in grid cell

seed= 101;
deltaP= 0.1 ; % pressure drop in Pa

% Generation d'une valeur de porosite suivant une distribution normale 
mean_poro = 0.9;
std_poro = 7.5e-3;
%poro= normrnd(mean_poro,std_poro) ; A VOIR PLUS TARD POUR LA DISTRIBUTION
%DE LA POROSITE
poro = 0.9 ;

mean_fiber_d= 12.5 ; % in microns
std_d= 2.85 ; % in microns
filename= 'fiber_mat.tiff' ;

% --------------------- Vérification de solution ------------------------
% Generation de trois solutions [f1, f2, f3] avec f1 la solution sur le
% maillage le plus fin

% Vecteur nombre de noeuds
r = 2;
maillage_fin = 200;
NX_vect = [maillage_fin, maillage_fin/r, maillage_fin/(2*r)];

% Initialisation vecteur solution
f = zeros(1, numel(NX_vect));

for i = 1:numel(NX_vect) ;
    nx = NX_vect(i)
    dx= 2e-4/nx ;
    [d_equivalent]=Generate_sample(seed,filename,mean_fiber_d,std_d,poro,nx,dx);
    f(i) = LBM(filename,nx,deltaP,dx,d_equivalent);
end
figure(3)
plot(NX_vect, f)
% Calcul ordre de convergence observé
p_hat = log((f(3)-f(2))/(f(2)-f(1)))/log(r)

% Ordre de convergence formel
p_f = 2;

% Calcul du Grid Convergence Index
ratio = abs((p_hat-p_f)/p_f);

if ratio <= 0.1
    GCI = (1.25/(r^p_f-1)) * abs(f(2)-f(1));
else 
    p = min(max(0.5, p_hat), p_f);
    GCI = (3/(r^p-1)) * abs(f(2)-f(1));

end
% Q A) Incertitude numerique
u_num = GCI/2;

%% --------------------- Propagation d'incertitude ------------------------
seed= 101;
deltaP= 0.1 ; % pressure drop in Pa
mean_poro = 0.9;
std_poro = 7.5e-3;
filename= 'fiber_mat.tiff' ;
mean_fiber_d= 12.5 ; % in microns
std_d= 2.85 ; % in microns

% Propagation des incertitudes sur les données d'entrées de porosité a
% l'aide de la methode de Monte-Carlo
% Etape 1 : Definir la PDF pour les variables d'entrée
mu = 0.9;
sigma =7.5e-3;

% Etape 2 : Generer un echantillon aleatoire de type lhs
sample = 100;
lhs = lhsdesign(sample,1);

porosity_sample = zeros(1, numel(lhs));
for i = 1:numel(lhs) 
    porosity_sample (i) = Inverse_CDF(lhs(i), mu, sigma);
end

% Etape 3 : Evaluer la SRQ pour chaque valeur de l'echantillon
NX = 100;
dx= 2e-4/NX ;
SRQ = zeros(1, numel(lhs));
for i = 1:numel(lhs)
    [d_equivalent]=Generate_sample(seed,filename,mean_fiber_d,std_d,porosity_sample(i),NX,dx);
    SRQ(i) = LBM(filename,NX,deltaP,dx,d_equivalent);
end

% Etape 4 : Tracer la CDF de la SRQ
figure(4);
cdfplot(SRQ);
title('CDF Permeabilite pour propagation de la porosité');

% Evaluation de u_input
S_bar = sum(SRQ)/numel(lhs);

% Q B) Incertitude sur les variables d'entree
u_input_2 = sum((SRQ-S_bar).^2)/(numel(lhs)-1); % Attention ici cest u_input^2

% ----------------- Incertitudes données experimentales--------------------
D = 80.6; % Mediane des valeurs experimentales
exp_sigma = 14.7; % exponentiel de l'ecart type

sr = (exp(D+2*exp_sigma)-exp(D-2*exp_sigma))/2; % Pas sure
br = 10;

% Q C) Incertitude sur les données experimentales
u_D = sqrt(sr^2+br^2);


%----------------- Erreur de simulation--------------------
% Q D)
E = S_bar - D;

% ----------------- Erreur de modelisation--------------------
u_val = sqrt(u_num^2+u_input_2+u_D^2);
k=2;
Upper_bound = E+k*u_val
Lower_bound = E-k*u_val


%--------------------------- Fonctions -----------------------------
function inv_cdf = Inverse_CDF(p, mu, sigma)
inv_cdf = norminv(p,mu,sigma)
end