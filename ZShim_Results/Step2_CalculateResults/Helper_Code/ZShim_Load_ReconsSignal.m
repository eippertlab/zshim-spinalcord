function signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode)

% load the extracted signal (.txt file) from different reconstructed volumes
% average for each slice, put the results in a matrix
% output: 3D matrix, (subjects x number of reconstructed volumes x slices)

% ----------
% Inputs:
% ----------

% processeddatapath: fullpath, string, processed data location
% subjects:          structure,subjects directory
% reconsName:        string, name of the reconstrcuted volume
% reconsMode:        cell of strings, reconstructed acc to which selection


% Merve Kaptan, mkaptan@cbs.mpg.de

% cut where all subjects have slices
ZMin = 710;
ZMax = 935;

for sub = 1:size(subjects,1)
    
    outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
    cd(outdir)
    
    for zmode = 1:numel(reconsMode)
        
        if ~contains(reconsMode{zmode}, '_normalized_GM.txt')
            
            tmp = load([ 'normalizedVols' filesep reconsName '_' reconsMode{zmode} '_normalized_GM.txt']);
            
        else
            
            tmp = load([ 'normalizedVols' filesep reconsMode{zmode}]);
            
        end
        
        tmp = tmp(3:4,:)';
        minz = find(tmp(:,1)==ZMin);
        maxz = find(tmp(:,1) == ZMax);
        tmp = tmp(minz(1):maxz(end),:);
        tmp(:,1) = tmp(:,1) + 1;
        
        slices = unique(tmp(:,1));
        
        for s = 1:size(slices,1)
            
            signalsli = [];
            
            signalsli =  nanmean(tmp(find(tmp(:,1) == slices(s)),2));
            
            results(sub,zmode,s) = signalsli;
            
        end
        clear tmp;
        
    end
    
end

signal = results;
