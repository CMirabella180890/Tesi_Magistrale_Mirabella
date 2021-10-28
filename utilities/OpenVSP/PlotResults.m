function PlotResults(set,res,name,config)

% Plots results.
%
% set: setup parameters
% res: simulation results

plotDegenPlate([name '_DegenGeom'])
saveas(gcf,[name '3D_Plate.png']);

plotDegenStick([name '_DegenGeom'])
saveas(gcf,[name '3D_Stick.png']);

plotDegenSurf([name '_DegenGeom'])
saveas(gcf,[name '3D_Surf.png']);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
plot(set.AoA,res.CL,'o-')
axis auto,xlabel('\alpha'), ylabel('C_L'), grid on, grid minor

subplot(2,2,2)
plot(res.CDtot,res.CL,'o-')
axis auto,xlabel('C_D'), ylabel('C_L'), grid on, grid minor

subplot(2,2,3)
plot(set.AoA,res.CMy,'o-')
axis auto,xlabel('\alpha'), ylabel('CM_y'), grid on, grid minor

subplot(2,2,4)
plot(res.Yavg',res.Cl') % plot wing loads for each AoA
axis auto, xlabel('Y'), ylabel('C_l'), grid on, grid minor

for ii = 1:length(res.CL)
    legendNames{ii} = ['AoA = ',num2str(set.AoA(ii)),' deg'];
end
legend(legendNames,'location','northeastoutside')

% Write figure title
if isfield(set,'flapRot') && isfield(set,'elevRot')
    sgtitle([name '  ' config, ...
        '  M = ' num2str(set.Mach) '  Re = ' num2str(set.ReCref/1e6,'%.2f') 'E6', ...
        '  flap = ', num2str(set.flapRot),' deg', ...
        '  elev = ', num2str(set.elevRot),' deg', ...
        ],'interpreter','none')
elseif isfield(set,'flapRot')
        sgtitle([name '  ' config, ...
        '  M = ' num2str(set.Mach) '  Re = ' num2str(set.ReCref/1e6,'%.2f') 'E6', ...
        '  flap = ', num2str(set.flapRot),' deg', ...
        ],'interpreter','none')
elseif isfield(set,'elevRot')    
        sgtitle([name '  ' config, ...
        '  M = ' num2str(set.Mach) '  Re = ' num2str(set.ReCref/1e6,'%.2f') 'E6', ...
        '  elev = ', num2str(set.elevRot),' deg', ...
        ],'interpreter','none')
else
    sgtitle([name '  ' config, ...
        '  M = ' num2str(set.Mach) '  Re = ' num2str(set.ReCref/1e6,'%.2f') 'E6', ...
        ],'interpreter','none')
end

saveas(gcf,[name '_plots.png']);

end