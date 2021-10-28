function MoveResults(name,directory)
% Moves results from main directory to results directory.
%
% name     : aircraft name
% directory: directory you're moving files to

fclose('all');
movefile([name '.vspscript'],directory);
movefile([name '_MassProps.txt'],directory);
movefile([name '_DegenGeom.*'],directory);
movefile([name '.vsp3'],directory);

try
    movefile('*.png',directory)
catch
    warning('No PNG image files to move')
end

try
    movefile([name '_ParasiteDrag.vspscript'],directory);
    movefile('Unnamed_*',directory)
    try % rename parasite drag files
        cd(directory)
        movefile('Unnamed_ParasiteBuildUp.csv',[name,'_ParasiteBuildUp.csv'])
        movefile('Unnamed_DragBuild.tsv',[name,'_DragBuild.tsv'])
        movefile('Unnamed_CompGeom.csv',[name,'_CompGeom.csv'])
        movefile('Unnamed_CompGeom.txt',[name,'_CompGeom.txt'])
        cd ..
    catch
        warning('Could not rename parasite drag file')
    end
catch
    warning('No parasite drag files to move')
end

end

