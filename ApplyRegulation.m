function Aircraft = ApplyRegulation(Aircraft)
%% DESCRIPTION OF THE FUNCTION
% Aircraft = ApplyRegulation(Aircraft)
%   This functions allows the user to automatically switch through various
%   regulations, basing the switch activation on a char value stored inside
%   a struct variable called 'Aircraft'. The field inside which the char
%   value is stored is 
%                    Aircraft.Certification.Regulation.value
%   Possible switch case are:
%   ---> CSVLA
%   ---> CS22 (currently not available) 
%   ---> CS23 (currently not available)
%   Inside each case a properly defined script is called, to apply
%   calculations methods stored inside class files. To manage in a more
%   clear manner all the operations every case will change the working
%   directory to the one which contains all the files needed to carry out
%   the calculations. Outputs are stored inside the struct variable
%   'Aircraft' and in a sub-folder called Output in .pdf files format.
%% SWITCH CASE TO SELECT APPLICABLE REGULATION 
switch (Aircraft.Certification.Regulation.value)
    % CASE 1: Very Light Aircraft
    case 'CSVLA'
        % Change working directory
        cd csvla
        % The 'dir' variable contains working directory path saved as a
        % char value
        dir = pwd;
        % Store working directory inside the log file
        fprintf('-----------------');
        fprintf('\n');
        fprintf('### Current directory ###');
        fprintf('\n');
        fprintf('%s\n', dir);
        % Call the class csvla.m
        fprintf('-----------------');
        fprintf('\n');
        fprintf('### Flight Envelope - per CS - VLA ###');
        fprintf('\n');
        obj = csvla;
        % Apply all the methods required
        CalcFlightEnv
        fprintf('-----------------');
        fprintf('\n');
        fprintf('### Balancing loads - per CS - VLA ###');
        fprintf('\n');
        CalcBalancLoads
        % CHANGE DIRECTORY TO CALCULATES C_l = C_l(y)
        cd .. 
        cd utilities\OpenVSP
        % cd utilities\v2
        % The 'dir' variable contains working directory path saved as a
        % char value
        dir = pwd;
        % Store working directory inside the log file
        fprintf('-----------------');
        fprintf('\n');
        fprintf('### Current directory ###');
        fprintf('\n');
        fprintf('%s\n', dir);
        
        %% STARTING OPEN VSP CALCULATION
        % Main_UAS
        % +++ ISTRUZIONE IMPORTANTE +++
        disp(" ++++ STARTING OPEN VSP CALCULATIONS ++++ ");
        diary off
        Main_PROSIB_WT 
        % +++ ISTRUZIONE IMPORTANTE +++
        % LOADING DATA INSIDE VARIABLE STRUCT
        % CHANGE DIRECTORY AND LOAD AERODATA
        FromTableToStructAircraft
        % CHANGE DIRECTORY TO CALCULATES SHEAR AND BENDING MOMENT
        cd .. 
        cd ..
        cd ..
        diary CMFlightLoads.txt
        disp(" ++++ FIGURE 9 - OPEN VSP RESULTS ++++ ");
        cd csvla
        % The 'dir' variable contains working directory path saved as a
        % char value
        dir = pwd;
        % Store working directory inside the log file
        fprintf('-----------------');
        fprintf('\n');
        fprintf('### Current directory ###');
        fprintf('\n');
        fprintf('%s\n', dir);
        % CALCULATE SHEAR, BENDING, AND TORSION MOMENT
        CalcShearBendTorsMom
        % CALCULATE UNSYMMETRICAL LOADS
        CalcUnsymmLoads
        % CALCULATE HORIZONTAL TAIL LOADS 
        CalcHorizTailLoads
        % CALCULATE VERTICAL TAIL LOADS
        CalcVertTailLoads
        % CALCULATE SUPPLEMENTARY CONDITIONS 
        CalcSuppLoads
        % FLIGHT ENVELOPE WITH FLAPS DEPLOYED
        CalcFlapsEnvelope
%     case 'CS23'
%     case 'CS22'
end
%% DISPLAYING RESULTS
end

