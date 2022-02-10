%% GEOMETRY EXECUTION PROCESS
% Input data check
% Output aircraft and components figures

%% aircraft name here
ac.name = Aircraft.Certification.Aircraft_Name.value;
str = ['### aircraft name is ' ac.name];
disp ('-----------------------------------------------------------------')
disp ('-----------------------------------------------------------------')
disp(' ');
disp(str);

%% Checking results directory
main_dir = pwd;                                 % current directory
res_dir  = [main_dir '\' ac.name '_results'];   % results directory
if exist(res_dir,'dir')
    rmdir(res_dir,'s')                          % removing directory
end
mkdir(res_dir)                                  % making directory
Aircraft.res_dir = res_dir;

%% Aicraft Components in-out
FuselageOut(Aircraft);
WingOut(Aircraft);
HorizontalOut(Aircraft);
VerticalOut(Aircraft);
EngineOut(Aircraft);
LandingGearOut(Aircraft);
%temporary whole aircraft
AircraftOut(Aircraft)

clear Aircraft.res_dir 
