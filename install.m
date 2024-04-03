% Install ABIP LP and QCP
% Intel MKL is required
clear; clc; close all;

addpath(genpath('.'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in MKL path here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mkl_path = '';
mkl_path = '/opt/intel/oneapi/mkl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abip_home = pwd;
build_qcp = false;

% Install ABIP-LP
fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
fprintf("Building and installing the ABIP-LP part. \n");
cd(fullfile(abip_home, 'src', 'abip-lp'));
make_abip;

% Install ABIP-QCP

fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
fprintf("Building and installing the ABIP-QCP part. \n");
cd(fullfile(abip_home, 'src', 'abip-qcp'));
if exist(mkl_path, 'dir')
    make_abip_qcp;
    build_qcp = true;
else
    fprintf("QCP is not installed due to missing MKL library. \n");
end 

% Go back
cd(abip_home);

% Run a toy example
fprintf("Running a toy example \n");

try 
    test_abip_install(1);
catch
    fprintf("ABIP installation failed \n");
    return;
end % End try

cd(abip_home);

fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
fprintf("Successfully installed ABIP \n");
