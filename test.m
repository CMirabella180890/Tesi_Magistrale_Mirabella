clear all; close all; clc;
%% COMMENTO INIZIALE
%  Prova del 18/07/2021 -- Implementazione della corretta variabile
%  struttura 'Aircraft', come velivolo di prova è stato selezionato un
%  sistema unmanned da certificare sotto la norma CS-VLA.

% Log command window output. Avoid appending data to existing file.
if exist('CMFlightLoads.txt','file')
    diary OFF
    delete CMFlightLoads.txt
end
diary CMFlightLoads.txt
%% STORING THE WORKING DIRECTORY 
dir = pwd;
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Working directory ###"); 
fprintf('\n');
fprintf('%s\n', dir);
fprintf("--------------------------------------");
fprintf('\n\n');
%% DEFINING A GLOBAL VARIABLE AIRCRAFT
global Aircraft
% Printing on screen the correct output
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Aircraft Flight Loads ###");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf('\n');
%% Input specification - TLAR
InputSource = "From File"; 
% "From File"
% "Custom"
%do not edit
% filename = 'Aircraft_test'; %do not edit
% filename = 'drone_vla';     %do not edit
filename = 'TecnamP92_input';  %do not edit

%% INIZIALIZATION OF AircraftSTRUCT VARIABLE (Just the minimal values)
cd initialization 
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);
Aircraft = FlightLoadsInitialize();
fprintf('\n');
fprintf('Certification used from file.');
% Aircraft = FromFileCertification_fun(Aircraft,filename);
Aircraft = FromFileCertification_funUNDERTEST(Aircraft,filename);
fprintf('\n');
pause(5/1000);
%% Pause
cd .. 
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);
pause(5/1000);
%% Aerodynamic data 
% In this section a convenient series of option are provided to insert all
% the necessary aerodynamic data required to evaluate the horizontal tail
% balancing loads which are related to the aerodynamic force that the
% horizontal empennage must produce to achieve a trimmed flight condition.
% These balancing loads ARE NOT RELATED to aerodynamic design of the
% aircraft empennage; airworthiness prescriptions of Subpart C are 
% concerned to structural design and sizing of the aircraft.
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% % Call the class aero_model.m
% % obj1 = aero_model; 
% % cd .. 
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
%% CALLING FUNCTIONS TO PLOT THE V - N DIAGRAMS
% Printing on screen the correct output
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Airworthiness regulations applied ###");
fprintf('\n');
fprintf('%s\n', Aircraft.Certification.Regulation.value);
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
%% CALLING FUNCTION TO APPLY ALL THE CALCULATIONS METHODS
Aircraft = ApplyRegulation(Aircraft);

%% CALCULATING SHEAR AND BENDING MOMENT ON THE WING

%% GEMOETRY SUBROUTINE 
cd .. 
cd utilities\Geometry
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% CALLING RepGenMain FUNCTION
Main_Geometry;

close all;

% %% REPORT GENERATOR 
% cd .. 
% cd RepGen
% % cd utilities\RepGen
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % CALLING RepGenMain FUNCTION
% RepGenMain(Aircraft);