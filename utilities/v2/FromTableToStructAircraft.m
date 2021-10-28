%% LIFT AND DRAG ALONG THE SPAN - CL = 1 
% In this script we use the Open VSP output to obtain the lift and drag
% coefficient distribution along the main wing span, in correspondence to a
% resultant lift coefficient equal to one.

res_dir = [Aircraft.Certification.Aircraft_Name.value '_results'];
cd(res_dir)
polar   = importdata([Aircraft.Certification.Aircraft_Name.value '_DegenGeom.polar']);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value    = polar.data(:,5);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDo.value   = polar.data(:,6);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDi.value   = polar.data(:,7);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDtot.value = polar.data(:,8);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CS.value    = polar.data(:,9);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.E.value     = polar.data(:,10);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.e.value     = polar.data(:,11);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFx.value   = polar.data(:,12);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFy.value   = polar.data(:,13);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFz.value   = polar.data(:,14);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMx.value   = polar.data(:,15);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMy.value   = polar.data(:,16);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMz.value   = polar.data(:,17);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMl.value   = polar.data(:,18);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMm.value   = polar.data(:,19);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMn.value   = polar.data(:,20);

%% LOCAL COEFFICIENTS FROM LOD FILE
% extract wing only aerodynamic loads for each angle of attack

% number of rows to read for each angle of attack
if isfield(ac,'hor') % if horizontal tail exist
    nRows = sum(ac.wing.utess)-length(ac.wing.utess) + ...
        sum(ac.hor.utess)-length(ac.hor.utess);
else % otherwise consider only wing
    nRows = sum(ac.wing.utess)-length(ac.wing.utess);
end

% number of lines to skip depending on the number of components
% (for some reason PodGeom are duplicates in addition to symmetry - at least in v3.19)
nomi = fieldnames(ac);
nWing = 0;
nFuse = 0;
nPod = 0;
for ii = 1:length(nomi)
    if isstruct(ac.(nomi{ii}))
        vspID = ac.(nomi{ii}).id;
        if (~isempty(strfind(ac.config,'W')) && strcmp(vspID,'wing')) || ...
                (~isempty(strfind(ac.config,'H')) && strcmp(vspID,'horizontal'))
            nWing = nWing + 1;
        end
        if ~isempty(strfind(ac.config,'B')) && strcmp(vspID,'fuselage')
            nFuse = nFuse + 1;
        end
        if ~isempty(strfind(ac.config,'P')) && strcmp(vspID,'pod')
            nPod = nPod + length(ac.(nomi{ii}).yloc) * 4; % yes they are added 4 times, don't know why
        end
    end
end

nLines = 18 +  + nRows+ nWing + nFuse + nPod;

if length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value) > 1
    startRow = [18, 18+(1:(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value)-1))*nLines];
    %    startRow = [18, 125, 232];
else
    startRow = 18; % it should start from 18
end
endRow = startRow + nRows -1;

cd .. 
copyfile ReadLODfile.m DroneVLA_results
cd(res_dir)

lod = ReadLODfile([Aircraft.Certification.Aircraft_Name.value '_DegenGeom.lod'],startRow,endRow);

for ii = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value)
    loads = lod(1+nRows*(ii-1):nRows*ii,:);
    for jj = 1:sum(ac.wing.utess)-length(ac.wing.utess) % only wing 1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Yavg.value(ii,jj) = loads(jj,3);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(ii,jj)   = loads(jj,6);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(ii,jj)   = loads(jj,7);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(ii,jj)   = loads(jj,13);
    end
end