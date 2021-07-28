function results = ZShim_SuppIII_II_II_III(scttemplatepath,rawdatapath,processdatapath,recalculateResults)

if recalculateResults
    
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    teShift = {'plus', 'minus'};
    
    zshimSteps = 21;
    fitSli = 8;
    beta = 0;
    res = [2.2 1 1]*1e-3;
    basis = [1:4];
    totalSli = 8;
    smoothfactor = '1';
    %%
    
    for t = 1:numel(teShift)
        
        for sub = 1:size(subjects,1)
            
            outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
            cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
            
            fmFilename = ['phase' smoothfactor '.nii.gz'];
            maskFilename = 'T2_mask.nii.gz';
            
            fm = read_avw(fmFilename);
            mask = read_avw(maskFilename);
            mask = logical(mask);
            
            dTE = 2.46e-3; %difference in echo time
            HzMax = 1/(2*dTE);
            fmMax = max(abs(fm(:)));
            fm = fm*HzMax/fmMax;
            
            basis = [1:4];
            res = [2.2 1 1]*1e-3; % the resolution of the phase map in meter
            
            epiSli = 24;          % number of axial slices in the EPI acquisition
            fmSli = 180;          % number of axial slices in the field map acquisition (after swapping the dimensions)
            
            % Cut the field map 5 slices above and 4 slices below of the coverage of the EPI slice
            % In our acquisition the slice stack was placed in a way that
            % the center of the EPI volume (24 slices) corresponds to the center of the
            % fm acquisition (180 sagittal slices) --> therefore,the middle point of the
            % acquisitions are identical
            fm = fm(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);
            mask = mask(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);
            
            % Compute the shims - main code
            [shims,~,~,~,shimSli] = compute_shims(fm,mask,basis,fitSli*res(3),res,beta);
            Hz_scale = [1 5e-3 5e-3 5e-3]; % this is the thickness of our EPI slice
            shims = diag(Hz_scale(basis))*shims;
            
            shimY = shims(find(basis==3),:)/5e-3;   % [Hz/m]
            Zshim = shims(find(basis==4),:);        % [Hz]
            
            dk = 1/64e-3; % [1/m] FOV 128mm, GRAPPA factor 2, previoulsy 2 pi so we had it in raidans, but everywhere else we had HZ
            dt = 0.93e-3; % [s] Echo spacing 0.93 ms, should be secionds!!
            G_PE = dk/dt; % [Hz/m]
            
            if isequal(teShift{t}, 'plus')
                Q = 1 + shimY/G_PE;
            elseif isequal(teShift{t}, 'minus')
                Q = 1 - shimY/G_PE;
            end
            
            % Check range of Q to force effective TE to lie within readout
            TE = 40e-3; % in seconds
            N0 = 24; % Number of readout lines before center of k-space (assuming 128 matrix, GRAPPA 2, Partial Fourier 7/8)
            N1 = 32; % Number of readout lines after center of k-space (assuming 128 matrix, GRAPPA 2, Partial Fourier 7/8)
            T0 = TE - N0*dt; % Start of readout
            T1 = TE + N1*dt; % End of readout
            
            if Q>TE/T0 % Case when refocusing would happen before start of readout - i.e. center of k-space is never reached
                Q = TE/T0;
            elseif Q<TE/T1 % Case when refocusing would happen after end of readout - again center of k-space is never reached
                Q = TE/T1;
            end
            
            Zshim = Zshim .* Q;           % [Hz]
            mP = 0.21e-3;                 % [T/m] maximum Z-shim gradient - corresponding to gradient moment yielding refocusing at TE
            mP = mP/floor(zshimSteps/2);  % [T/m] single Z-shim step gradient moment
            gamma = 42.6e6;               % [Hz/T] gyromagnetic ratio
            slice = 5e-3;                 % [m] this is the thickness of our EPI slice
            mHz = mP*gamma*slice;         % [Hz] Z-shim frequency range within slice
            pHz = mHz;                    % [Hz]
            Zshim_picks = Zshim/pHz + ceil(zshimSteps/2);
            Zshim_picks =round(Zshim_picks);
            
            %% Create results matrix fitArray
            % rescaled, rounded Zshim picks for each slice in the saggital fieldmap
            % (will be used as an input to sequence)
            fitArray(:,1)= Zshim_picks;
            % which slices were used for fitting - output from main code
            fitArray(:,2) = shimSli;
            
            % which 5 slices of the field map correspond to the one EPI slice
            if totalSli == 4
                vector = [NaN(1,3) repmat([1 2 3 4 5], 1, 24) NaN(1,3)]';
            elseif totalSli == 8
                vector = [NaN(1,1) repmat([1 2 3 4 5], 1, 24) NaN(1,1)]';
            elseif totalSli == 12
                vector = repmat([1 2 3 4 5], 1, 24)';
                vector([1 end]) = [];
            end
            
            fitArray(:,3) = vector;
            
            auto_gre = fitArray((fitArray(:,3)==3),1);
            save(['AutoSelection_FM_gre_TE_' teShift{t}  '.mat'], 'auto_gre')
            
            clear fitArray Zshim_picks fm mask Zshim_picks auto_gre
            
            
        end
    end
    
    %% Reconstruct artificial volumes and process as before
    
    for sub = 1:size(subjects,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
        
        % copy the files to the reconstruction of z-shim values folder
        copyfile('*TE_plus.mat', outdir)
        copyfile('*TE_minus.mat', outdir)
        
        cd(outdir)
        
        ZrefName = 'ZRef1.nii.gz';
        
        % reconstruct the artificial volumes based on plus and minus TE shift
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_gre_TE_plus.mat' , 'ZShim1_Recons_gre_TE_plus')
        
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_gre_TE_minus.mat' , 'ZShim1_Recons_gre_TE_minus')
        
        % to normalize these to template space
        % 1. register to moco image (as moco mask will be used)
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_gre_TE_plus.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o gre_TE_plus_MocoReg.nii.gz']);
        
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_gre_TE_minus.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o gre_TE_minus_MocoReg.nii.gz']);
        
        % 2. apply the warping field (moco volume to PAM50 --> to normalize to the template space)
        system(['sct_apply_transfo -i gre_TE_plus_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline ' ...
            ' -o normalizedVols' filesep 'gre_TE_plus_normalized.nii.gz'])
        
        system(['sct_apply_transfo -i gre_TE_minus_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline' ...
            ' -o normalizedVols' filesep 'gre_TE_minus_normalized.nii.gz'])
        
        % 3. extract the signal in the gray matter
        system(['fslmeants -i normalizedVols' filesep 'gre_TE_plus_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_gre_TE_plus_normalized_GM.txt ' ...
            ' --showall']);
        
        system(['fslmeants -i normalizedVols' filesep 'gre_TE_minus_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_gre_TE_minus_normalized_GM.txt ' ...
            ' --showall']);
        
    end
    
    reconsName  = 'ZShim1_Recons';
    reconsMode    = {'gre', 'gre_TE_plus', 'gre_TE_minus'};
    recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);
    
    for r = 1:numel(reconsMode)
        
        eval([reconsMode{r}  '= recons_signal(' num2str(r) ');'])
        
    end
    
else
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'FM_effectiveecho', 'results.mat'))
end

results_effectiveTE =  ZShim_CalculateResults(gre,gre_TE_plus,gre_TE_minus,0,1);

results{1,1} = 'results_gradientsAP';
results{1,2} = results_effectiveTE;

end

