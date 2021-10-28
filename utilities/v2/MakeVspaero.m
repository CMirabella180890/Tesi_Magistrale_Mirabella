function MakeVspaero(name,set)

% Makes vspaero file.
%
% name: aircraft name
% set : setup parameters

fileid = fopen([name '_DegenGeom.vspaero'],'w+');     % opening vspaero file

% filling in vspaero file with simulation parameters
fprintf(fileid,'%s\r\n',['Sref = ' num2str(set.Sref)]);
fprintf(fileid,'%s\r\n',['Cref = ' num2str(set.Cref)]);
fprintf(fileid,'%s\r\n',['Bref = ' num2str(set.Bref)]);
fprintf(fileid,'%s\r\n',['X_cg = ' num2str(set.X_cg)]);
fprintf(fileid,'%s\r\n',['Y_cg = ' num2str(set.Y_cg)]);
fprintf(fileid,'%s\r\n',['Z_cg = ' num2str(set.Z_cg)]);
fprintf(fileid,'%s\r\n',['Mach = ' num2str(set.Mach)]);
fprintf(fileid,'%s','AoA = ');
for i = 1:numel(set.AoA)-1
    fprintf(fileid,'%s',[num2str(set.AoA(i)) ',']);
end
fprintf(fileid,'%s\r\n',num2str(set.AoA(end)));
fprintf(fileid,'%s\r\n',['Beta = ' num2str(set.Beta)]);
fprintf(fileid,'%s\r\n',['Vinf = ' num2str(set.Vinf)]);
fprintf(fileid,'%s\r\n',['Rho = ' num2str(set.rho)]);
fprintf(fileid,'%s\r\n',['ReCref = ' num2str(set.ReCref)]);
fprintf(fileid,'%s\r\n',['ClMax = ' num2str(set.ClMax)]);
fprintf(fileid,'%s\r\n',['MaxTurningAngle = ' num2str(set.MaxTurningAngle)]);
fprintf(fileid,'%s\r\n',['Symmetry = ' set.Symmetry]);
fprintf(fileid,'%s\r\n',['FarDist = ' num2str(set.FarDist)]);
fprintf(fileid,'%s\r\n',['NumWakeNodes = ' num2str(set.NumWakeNodes)]);
fprintf(fileid,'%s\r\n',['WakeIters = ' num2str(set.WakeIters)]);
% Write rotor data to vspaero file only if RPM > 0
if isfield(set,'rpm')
    enabledRot = find(set.rpm); % IDs of enabled rotors
    nRot = numel(enabledRot);
else
    nRot = 0;
end
fprintf(fileid,'%s\r\n',['NumberOfRotors = ' num2str(nRot)]);
if nRot > 0
    for ii = 1:nRot
        fprintf(fileid,'%s\r\n',['PropElement_' num2str(ii)]);
        fprintf(fileid,'%s\r\n',num2str(ii));
        fprintf(fileid,'%s\r\n',[num2str(set.xloc(enabledRot(ii)),'%.4f'),' ',...
            num2str(set.yloc(enabledRot(ii)),'%.4f'),' ',...
            num2str(set.zloc(enabledRot(ii)),'%.4f')]);
        fprintf(fileid,'%s\r\n',[num2str(set.xdir(enabledRot(ii)),'%.4f'),' ',...
            num2str(set.ydir(enabledRot(ii)),'%.4f'),' ',...
            num2str(set.zdir(enabledRot(ii)),'%.4f')]);
        fprintf(fileid,'%s\r\n',num2str(set.dia(enabledRot(ii))/2)); % radius
        fprintf(fileid,'%s\r\n',num2str(set.hub(enabledRot(ii))/2)); % radius
        fprintf(fileid,'%s\r\n',num2str(set.rpm(enabledRot(ii))));
        fprintf(fileid,'%s\r\n',num2str(set.ct(enabledRot(ii))));
        fprintf(fileid,'%s\r\n',num2str(set.cp(enabledRot(ii))));
    end
end
if isfield(set,'movables') % check if the field has been defined from main file
    if ~isempty(set.movables)
        fprintf(fileid,'%s\r\n',['NumberOfControlGroups = ',num2str(length(set.movables))]);
        for ii = 1:length(set.movables)
            switch set.movables(ii)
                case 'f'
                    fprintf(fileid,'%s\r\n','Flap'); % CS Group custom name
                    fprintf(fileid,'%s\r\n',[set.wingID,'_Surf0_flap,',set.wingID,'_Surf1_flap']);
                    fprintf(fileid,'%s\r\n','1, -1'); % symmetric deflection, positive te down
                    fprintf(fileid,'%s\r\n',num2str(set.flapRot));
                case 'a'
                    % TODO
                case 'e'
                    fprintf(fileid,'%s\r\n','Elevator');
                    fprintf(fileid,'%s\r\n',[set.horID,'_Surf0_elevator,',set.horID,'_Surf1_elevator']);
                    fprintf(fileid,'%s\r\n','-1, 1'); % symmetric deflection, positive te up
                    fprintf(fileid,'%s\r\n',num2str(set.elevRot));
                case 'r'
                    % TODO
            end
        end
    end
end
fprintf(fileid,'%s\r\n','Preconditioner = Matrix');
fprintf(fileid,'%s\r\n',['Vortex Lift = ',set.VRTXLift]);
fprintf(fileid,'%s\r\n',['LE Suction = ',set.LESuction]);
fprintf(fileid,'%s\r\n',['Karman-Tsien Correction = ',set.KT]);

fclose(fileid);                                           % closing vspaero file

end
