% MATLAB script to launch a convergence study  in space using the provided
% LBM function as a black box
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


%% definition of convergence study arrays
nx_array = [50,100,200,400]; % number of points
dx_array = 2.0e-4./nx_array; % grid size in m

%% other inputs
seed=54321;
deltaP= 0.1 ; % pressure drop in Pa
poro= 0.9 ;
mean_fiber_d= 12.5 ; % in microns
std_d= 2.85 ; % in microns
filename= 'fiber_mat.tiff' ;

%% tracking output variable arrays
poro_out_array = zeros(1,length(nx_array));
Re_out_array = zeros(1,length(nx_array));
k_out_array = zeros(1,length(nx_array));



%% generation of the fiber structure
for i_convstud = 1:length(nx_array)
    tic
    disp(i_convstud)
    NX = nx_array(i_convstud);
    dx = dx_array(i_convstud);
    filename_temp = strcat('fiber_mat',int2str(i_convstud),'.tiff')
    [d_equivalent]=Generate_sample(seed,filename_temp,mean_fiber_d,std_d,poro,NX,dx);

    % calculation of the flow field and the permeability from Darcy Law
    outtemp = LBM(filename_temp,NX,deltaP,dx,d_equivalent);
    poro_out_array(i_convstud) = outtemp(1);
    Re_out_array(i_convstud) = outtemp(2);
    k_out_array(i_convstud) = outtemp(3);
    toc
    
end

%% writing output to text files
header = {'delta_x', 'porosity', 'k'}
filename = 'tempoutputfile.dat';
writeDataToTextFile(filename, header, dx_array, poro_out_array, k_out_array);

