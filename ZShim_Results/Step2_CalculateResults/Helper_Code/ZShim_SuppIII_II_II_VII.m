function results = ZShim_SuppIII_II_II_VII(scttemplatepath,rawdatapath,processdatapath,recalculateResults, fslDir)

% Use the histogram-based z-shim selections (see the files in the following folder: ZShim_SuppIII_II_II_VII_HistogramEvaluation) to
% calculate new zshims, create artificial volumes based on these modified
% zshims and extract signal.
% Compare modified implementation signal to previous one and return the result.
% If recalculateResults is set to False, no new data will be created, but
% rather the pre-saved results will be loaded.

% ----------
% Inputs:
% ----------

% scttemplatepath:   fullpath, string, PAM50 template location
% rawdatapath:       fullpath, string, raw data location
% processeddatapath: fullpath, string, processed data location
% recalculateResults: True or False

% Merve Kaptan, mkaptan@cbs.mpg.de

if recalculateResults
    
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    for s = 1:size(subjects,1)
        

        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
    
        % save the output of IDL script (ZShim_SuppIII_II_II_VII_HistogramEvaluation/EvalzShimsFromB0.pro) as a .mat file
        
         tmp = load('zShimsMean_vs_Histogram.dat');
         hist_picks = tmp(:,3);
         hist_picks(hist_picks<=0) = 1;
         save('AutoSelection_gre_histogram.mat', 'hist_picks')
        
        ZrefName = 'ZRef1.nii.gz';
        
        % reconstruct the artificial volumes based on eroded and dilated masks
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_gre_histogram.mat' , 'ZShim1_Recons_gre_histogram')
     
        % to normalize these to template space
        % 1. register to moco image (as moco mask will be used)
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_gre_histogram.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o gre_histogram_MocoReg.nii.gz']);
        
        % 2. apply the warping field (moco volume to PAM50 --> to normalize to the template space)
        system(['sct_apply_transfo -i gre_histogram_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline ' ...
            ' -o normalizedVols' filesep 'gre_histogram_normalized.nii.gz'])
        
        % 3. extract the signal in the gray matter
        system(['fslmeants -i normalizedVols' filesep 'gre_histogram_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_gre_histogram_normalized_GM.txt ' ...
            ' --showall']);
       
    end
    
    % extract signal
    reconsName  = 'ZShim1_Recons';
    reconsMode    = {'gre', 'gre_histogram',  'no'};
    recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);
    
    for r = 1:numel(reconsMode)
        
        eval([reconsMode{r}  '= recons_signal(' num2str(r) ');'])
        
    end
    
else
    
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'FM_histogram', 'results.mat'))
end


results_histogram1  =  ZShim_CalculateResults(no,gre,gre_histogram,0,1);
 
results{1,1} = 'results_histogram- original vs histogram';
results{1,2} = results_histogram1;

results_histogram2  =  ZShim_CalculateResults(no,gre_histogram,manual,0,1);

results{2,1} = 'results_histogram- histogram vs manual';
results{2,2} = results_histogram2;

end
