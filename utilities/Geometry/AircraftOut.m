function AircraftOut(Aircraft)
%Moving to Results directory
cd(Aircraft.res_dir)

comp = 'Fuselage';
if isfield(Aircraft.Geometry, comp) == 1
    fus_top     = openfig([comp '-Top-View.fig']);
    %edit figure title
    title('Aircraft Top-View')
    fus_side    = openfig([comp '-Side-View.fig']);
    %edit figure title
    title('Aircraft Side-View')
    fus_front   = openfig([comp '-Front-View.fig']);
    %edit figure title
    title('Aircraft Front-View')
    fus_3D      = openfig([comp '3D.fig']);
    %edit figure title
    title('Aircraft 3D-View')
else
    disp([comp 'does not exist'])
    disp([comp 'must be assigned'])
    return %fuselage must be always assigned in geometry
end

comp = 'Wing';
if isfield(Aircraft.Geometry, comp) == 1
    win_top = openfig([comp '-Top-View.fig'],'invisible');
    win_side = openfig([comp '-Side-View.fig'],'invisible');
    win_front = openfig([comp '-Front-View.fig'],'invisible');
    win_3D = openfig([comp '3D.fig'],'invisible');
    copyobj(get(get(win_top,'children'),'children'), get(fus_top,'children'));
    copyobj(get(get(win_side,'children'),'children'), get(fus_side,'children'));
    copyobj(get(get(win_front,'children'),'children'), get(fus_front,'children'));
    copyobj(get(get(win_3D,'children'),'children'), get(fus_3D,'children'));
    close(win_top); close(win_side); close(win_front); close(win_3D);
else
    disp([comp 'does not exist'])
end


comp = 'Horizontal';
if isfield(Aircraft.Geometry, comp) == 1
    hor_top = openfig([comp '-Top-View.fig'],'invisible');
    hor_side = openfig([comp '-Side-View.fig'],'invisible');
    hor_front = openfig([comp '-Front-View.fig'],'invisible');
    hor_3D = openfig([comp '3D.fig'],'invisible');
    copyobj(get(get(hor_top,'children'),'children'), get(fus_top,'children'));
    copyobj(get(get(hor_side,'children'),'children'), get(fus_side,'children'));
    copyobj(get(get(hor_front,'children'),'children'), get(fus_front,'children'));
    copyobj(get(get(hor_3D,'children'),'children'), get(fus_3D,'children'));
    close(hor_top); close(hor_side); close(hor_front); close(hor_3D);
else
    disp([comp 'does not exist'])
end


comp = 'Vertical';
if isfield(Aircraft.Geometry, comp) == 1
    ver_top = openfig([comp '-Top-View.fig'],'invisible');
    ver_side = openfig([comp '-Side-View.fig'],'invisible');
    ver_front = openfig([comp '-Front-View.fig'],'invisible');
    ver_3D = openfig([comp '3D.fig'],'invisible');
    copyobj(get(get(ver_top,'children'),'children'), get(fus_top,'children'));
    copyobj(get(get(ver_side,'children'),'children'), get(fus_side,'children'));
    copyobj(get(get(ver_front,'children'),'children'), get(fus_front,'children'));
    copyobj(get(get(ver_3D,'children'),'children'), get(fus_3D,'children'));
    close(ver_top); close(ver_side); close(ver_front);  close(ver_3D);
else
    disp([comp 'does not exist'])
end

comp = 'Engine';
if isfield(Aircraft.Geometry, comp) == 1
    eng_top = openfig([comp '-Top-View.fig'],'invisible');
    eng_side = openfig([comp '-Side-View.fig'],'invisible');
    eng_front = openfig([comp '-Front-View.fig'],'invisible');
    eng_3D = openfig([comp '3D.fig'],'invisible');
    copyobj(get(get(eng_top,'children'),'children'), get(fus_top,'children'));
    copyobj(get(get(eng_side,'children'),'children'), get(fus_side,'children'));
    copyobj(get(get(eng_front,'children'),'children'), get(fus_front,'children'));
    copyobj(get(get(eng_3D,'children'),'children'), get(fus_3D,'children'));
    close(eng_top); close(eng_side); close(eng_front); close(eng_3D);
else
    disp([comp 'does not exist'])
end

comp = 'LandingGear';
%if isfield(Aircraft.Geometry, comp) == 1
eng_top = openfig([comp '-Top-View.fig'],'invisible');
eng_side = openfig([comp '-Side-View.fig'],'invisible');
eng_front = openfig([comp '-Front-View.fig'],'invisible');
eng_3D = openfig([comp '3D.fig'],'invisible');
copyobj(get(get(eng_top,'children'),'children'), get(fus_top,'children'));
copyobj(get(get(eng_side,'children'),'children'), get(fus_side,'children'));
copyobj(get(get(eng_front,'children'),'children'), get(fus_front,'children'));
copyobj(get(get(eng_3D,'children'),'children'), get(fus_3D,'children'));
close(eng_top); close(eng_side); close(eng_front); close(eng_3D);
% else
%     disp([comp 'does not exist'])
% end
disp ('-----------------------------------------------------------------')
disp ('-----------------------------------------------------------------')
disp('Saving complete aircaft figures...')

%2d
saveas(fus_top,'Top-View.fig');
saveas(fus_top,'Top-View.png');
saveas(fus_side,'Side-View.fig');
saveas(fus_side,'Side-View.png');
saveas(fus_front,'Front-View.fig');
saveas(fus_front,'front-View.png');

%3D
saveas(fus_3D,'Aircraft3D.fig');
saveas(fus_3D,'Aircraft3D.png');
view(0,0.5)
title('Aircraft 3D Side-View')
saveas(gcf,'SideView3D.fig');
saveas(gcf,'SideView3D.png');
title('Aircraft 3D Front-View')
view(-90,0.0)
saveas(gcf,'FrontView3D.fig');
saveas(gcf,'FrontView3D.png');
title('Aircraft 3D Top-View')
view(0,90.0)
saveas(gcf,'TopView3D.fig');
saveas(gcf,'TopView3D.png');

cd ..

if exist('engine.mat') == 2
    movefile('engine.mat', Aircraft.res_dir);
end
if exist('fuselage.mat') == 2
    movefile('fuselage.mat', Aircraft.res_dir);
end
if exist('wing.mat') == 2
    movefile('wing.mat', Aircraft.res_dir);
end
if exist('vertical.mat') == 2
    movefile('vertical.mat', Aircraft.res_dir);
end
if exist('horizontal.mat') == 2
    movefile('horizontal.mat', Aircraft.res_dir);
end
if exist('landing_gear.mat') == 2
    movefile('landing_gear.mat', Aircraft.res_dir);
end

disp ('-----------------------------------------------------------------')
disp ('-----------------------------------------------------------------')
disp ('End Geometry utilities run')

end