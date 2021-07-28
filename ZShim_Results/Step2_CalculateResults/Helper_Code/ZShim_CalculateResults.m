function results = ZShim_CalculateResults(baseline,manipulated1,manipulated2,calculateVar,numberofcomp)
% take signal matrix for each z-shim condition as an input and
% 1. do the statistical tests (t-test) for mean signal & coefficient of
%   variation
% 2. calculate percent increase in the mean signal & decrease in the
%   coefficient of variation
% 3. provide confidence intervals for percent differences in mean signal
%   and coefficient of variation
%
%----------------------
% inputs
% ---------------------
% baseline:         baseline data
% manipulated1:     data obtained under first manipulated z-shim condition
% manipulated2:     data obtained under second manipulated z-shim condition 
%                   (EPI or FM-based)
%                   N subjects x N slices or 
%                   N subjects x 1 (average of slices)
% calculateVar:     1 or 0 (whether the same parameters should be calculated for 
%                   the variance as well, not possible if manipulated2 is entered)
% numberofcomp:     number of comparisons (for Bonferroni-correction of p value)
%
% ---------------------
% output
% ---------------------
% results:          a cell matrix containing:
%                   (corrected) p-value, t-value, degrees of freedom,
%                   percent increase/decrease, bootstrapped 95% CI of
%                   1) mean signal
%                   2) coefficient of variation
%                   3) variance (if selected as input option)

% Merve Kaptan, mkaptan@cbs.mpg.de; mervekaptan5@gmail.com
% 04.07.2021


if isempty(manipulated2)
    
    baselineMean = mean(baseline,2);
    manipulated1Mean = mean(manipulated1,2);
    
    zshimSubjectMeans = [baselineMean manipulated1Mean];
    zshimSubjectStds = [std(baseline,[],2)  std(manipulated1,[],2)];
    
    baselineCoeffVar = zshimSubjectStds(:,1)./zshimSubjectMeans(:,1);
    manipulated1CoeffVar = zshimSubjectStds(:,2)./zshimSubjectMeans(:,2);
    zshimSubjectCoeffVar = [baselineCoeffVar manipulated1CoeffVar];
    
    % t test for mean signal
    [~,p,~,stats] = ttest(manipulated1Mean, baselineMean ,'Tail','right');
    
    results{1,1} = 'mean signal';
    results{1,2} = 'p';
    results{1,3} = 't-value';
    results{1,4} = 'dof';
    
    if numberofcomp > 1
        p = p * numberofcomp;
    end
    
    results{2,2} = p;
    results{2,3} = stats.tstat;
    results{2,4} = stats.df;
    
    % percent increase for mean signal
    zshimGlobalMeans = nanmean(zshimSubjectMeans);
    meanG = zshimGlobalMeans(1);
    meanW = zshimGlobalMeans(2);
    meanP = 100 * (meanW/meanG);
    meanP_Increase = meanP - 100;
    
    results{1,5} = 'percent increase';
    results{2,5} = round(meanP_Increase,1);
    
    % Confidence intervals for percent differences in mean signal
    meanGall = zshimSubjectMeans(:,1);
    meanWall = zshimSubjectMeans(:,2);
    
    if or(sum(isnan(meanGall)),sum(isnan(meanWall)))
        meanGall(isnan(meanGall)) = [];
        meanWall(isnan(meanWall)) = [];
    end
    
    CIs = bootci(10000, {@pD_Group_inc, meanGall, meanWall},'type', 'per');
    
    results{1,6} = 'bootstrapped CIs (95%) of percent differences';
    results{2,6} =  round(CIs,1);
    
    % t test for coefficient of variation
    [~,p,~,stats] = ttest(baselineCoeffVar, manipulated1CoeffVar ,'Tail','right');
    
    results{3,1} = 'coefficient of var';
    results{3,2} = 'p';
    results{3,3} = 't-value';
    results{3,4} = 'dof';
    
    results{4,2} = p;
    results{4,3} = stats.tstat;
    results{4,4} = stats.df;
    
    % percent decrease for coefficient of variation
    zshimGlobalCoeffVar = nanmean(zshimSubjectCoeffVar);
    CoeffvarG = zshimGlobalCoeffVar(1);
    CoeffvarW = zshimGlobalCoeffVar(2:end);
    CoeffvarP = 100*CoeffvarW./CoeffvarG;
    CoeffvarP_Decrease = 100 - CoeffvarP;
    
    results{3,5} = 'percent decrease';
    results{4,5} = round(CoeffvarP_Decrease,1);
    
    % Confidence intervals for percent differences in coefficient of
    % variation
    coeffGall = zshimSubjectCoeffVar(:,1);
    coeffWall = zshimSubjectCoeffVar(:,2:end);
    
    if or(sum(isnan(coeffGall)),sum(isnan(coeffWall)))
        coeffGall(isnan(coeffGall)) = [];
        coeffWall(isnan(coeffWall)) = [];
    end
    
    if sum(coeffGall)~=0 || sum(coeffWall)~=0 
        CIs = bootci(10000, {@pD_Group_dec, coeffGall, coeffWall}, 'type', 'per');
        
        results{3,6} = 'bootstrapped CIs (95%) of percent differences';
        results{4,6} =  round(CIs,1);
    end
    
    if  calculateVar
        baselineVar  = var(baseline,[],2);
        manipulated1Var = var(manipulated1,[],2);
        
        zshimSubjectVariances = [baselineVar manipulated1Var];
        
        % t test for variance
        [~,p,~,stats] = ttest(baselineVar, manipulated1Var ,'Tail','right');
        
        results{5,1} = 'var';
        results{5,2} = 'p';
        results{5,3} = 't-value';
        results{5,4} = 'dof';
        
        results{6,2} = p;
        results{6,3} = stats.tstat;
        results{6,4} = stats.df;
        
        % percent decrease for variance
        zshimGlobalVariances = mean(zshimSubjectVariances);
        varG = zshimGlobalVariances(1);
        varW = zshimGlobalVariances(2:end);
        varP = 100*varW./varG;
        varP_Decrease = 100 - varP;
        
        results{5,5} = 'percent decrease';
        results{6,5} = round(varP_Decrease,1);
        
        % Confidence intervals for percent differences in variance
        varGall = zshimSubjectVariances(:,1);
        varWall = zshimSubjectVariances(:,2:end);
        
        if or(sum(isnan(varGall)),sum(isnan(varWall)))
            varGall(isnan(varGall)) = [];
            varWall(isnan(varWall)) = [];
        end
        
        if sum(varGall)~=0 || sum(varWall)~=0 
            CIs = bootci(10000, {@pD_Group_dec, varGall, varWall}, 'type', 'per');
            
            results{5,6} = 'bootstrapped CIs (95%) of percent differences';
            results{6,6} =  round(CIs,1);
        end
        
    end
    
elseif ~isempty(manipulated2)
    
    baselineMean = mean(baseline,2);
    manipulated1Mean = mean(manipulated1,2);
    manipulated2Mean = mean(manipulated2,2);
    
    zshimSubjectMeans = [baselineMean manipulated1Mean,manipulated2Mean];
    zshimSubjectStds = [std(baseline,[],2)  std(manipulated1,[],2) std(manipulated2,[],2)];
    
    baselineCoeffVar = zshimSubjectStds(:,1)./zshimSubjectMeans(:,1);
    manipulated1CoeffVar = zshimSubjectStds(:,2)./zshimSubjectMeans(:,2);
    manipulated2CoeffVar = zshimSubjectStds(:,3)./zshimSubjectMeans(:,3);
    zshimSubjectCoeffVar = [baselineCoeffVar manipulated1CoeffVar manipulated2CoeffVar];
    
    % t tests for mean signals
    [~,p1,~,stats1] = ttest(manipulated1Mean, baselineMean ,'Tail','right');
    [~,p2,~,stats2] = ttest(manipulated2Mean, baselineMean ,'Tail','right');
    [~,p3,~,stats3] = ttest(manipulated1Mean, manipulated2Mean );
    
    results{1,1} = 'mean signal-manipulated1vsbaseline';
    results{1,2} = 'p';
    results{1,3} = 't-value';
    results{1,4} = 'dof';
    
    results{3,1} = 'mean signal-manipulated2vsbaseline';
    results{3,2} = 'p';
    results{3,3} = 't-value';
    results{3,4} = 'dof';
    
    results{5,1} = 'mean signal-manipulated1vsmanipulated2';
    results{5,2} = 'p';
    results{5,3} = 't-value';
    results{5,4} = 'dof';
    
    if numberofcomp > 1
        p1 = p1 * numberofcomp;
        p2 = p2 * numberofcomp;
        p3 = p3 * numberofcomp;
        
    end
    
    results{2,2} = p1;
    results{2,3} = stats1.tstat;
    results{2,4} = stats1.df;
    
    results{4,2} = p2;
    results{4,3} = stats2.tstat;
    results{4,4} = stats2.df;
    
    results{6,2} = p3;
    results{6,3} = stats3.tstat;
    results{6,4} = stats3.df;
    
    % percent increases for mean signals
    zshimGlobalMeans = mean(zshimSubjectMeans);
    meanG = zshimGlobalMeans(1);
    meanW = zshimGlobalMeans(2:end);
    meanP = 100 * (meanW/meanG);
    meanP_Increase = meanP - 100;
    
    results{1,5} = 'percent increase';
    results{2,5} = round(meanP_Increase(1),1);
    results{3,5} = 'percent increase';
    results{4,5} = round(meanP_Increase(2),1);
    
    % Confidence intervals for percent differences in mean signals
    meanGall = zshimSubjectMeans(:,1);
    meanWall = zshimSubjectMeans(:,2:end);
    
    CIs = bootci(10000, {@pD_Group_inc, meanGall, meanWall}, 'type', 'per');
    
    results{1,6} = 'bootstrapped CIs (95%) of percent differences';
    results{2,6} =  round(CIs(:,1),1);
    results{3,6} = 'bootstrapped CIs (95%) of percent differences';
    results{4,6} =  round(CIs(:,2),1);
    
    % t tests for coefficients of variation
    [~,p1,~,stats1] = ttest(baselineCoeffVar, manipulated1CoeffVar ,'Tail','right');
    [~,p2,~,stats2] = ttest(baselineCoeffVar, manipulated2CoeffVar ,'Tail','right');
    [~,p3,~,stats3] = ttest( manipulated1CoeffVar, manipulated2CoeffVar );
    
    results{7,1} = 'coefficient of var- baselinevsmanipulated1';
    results{7,2} = 'p';
    results{7,3} = 't-value';
    results{7,4} = 'dof';
    
    results{9,1} = 'coefficient of var- baselinevsmanipulated2';
    results{9,2} = 'p';
    results{9,3} = 't-value';
    results{9,4} = 'dof';
    
    results{11,1} = 'coefficient of var- manipulated1vsmanipulated2';
    results{11,2} = 'p';
    results{11,3} = 't-value';
    results{11,4} = 'dof';
    
    if numberofcomp > 1
        p1 = p1 * numberofcomp;
        p2 = p2 * numberofcomp;
        p3 = p3 * numberofcomp;
        
    end
    
    results{8,2} = p1;
    results{8,3} = stats1.tstat;
    results{8,4} = stats1.df;
    
    results{10,2} = p2;
    results{10,3} = stats2.tstat;
    results{10,4} = stats2.df;
    
    results{12,2} = p3;
    results{12,3} = stats3.tstat;
    results{12,4} = stats3.df;
    
    % percent decreases for coefficients of variation
    zshimGlobalCoeffVar = mean(zshimSubjectCoeffVar);
    CoeffvarG = zshimGlobalCoeffVar(1);
    CoeffvarW = zshimGlobalCoeffVar(2:end);
    CoeffvarP = 100*CoeffvarW./CoeffvarG;
    CoeffvarP_Decrease = 100 - CoeffvarP;
    
    results{7,5} = 'percent decrease';
    results{8,5} = round(CoeffvarP_Decrease(1),1);
    results{9,5} = 'percent decrease';
    results{10,5} = round(CoeffvarP_Decrease(2),1);
    
    % Confidence intervals for percent differences in coefficients of
    % variation
    coeffGall = zshimSubjectCoeffVar(:,1);
    coeffWall = zshimSubjectCoeffVar(:,2:end);

    if or(sum(coeffGall)~=0,sum(coeffWall)~=0)
        CIs = bootci(10000, {@pD_Group_dec, coeffGall, coeffWall}, 'type', 'per');
        
        results{7,6} = 'bootstrapped CIs (95%) of percent differences';
        results{8,6} = round(CIs(:,1),1);
        results{9,6} = 'bootstrapped CIs (95%) of percent differences';
        results{10,6} = round(CIs(:,2),1);

    end
    
end

end