% One degree of freedom:
% Examine how [metric] varies as a function of differential drug selective pressure (alpha1/alpha2), and how this is
% affected by varying overall drug selective pressure (alpha_s)
close all;clear all;clc;

%settings
saveFigs = false;
plotsDir = './plots/';
plotRatios = true;
plotKineticsInd = false;
plotKineticsOverall = false;
odefunc = @Models.ODEexp;
dosefunc = @Models.dose_exp;
useAnalytical = false;
solvefunc = @Models.solved_ODEexp; %applicable if using analytical solution
evalfunc = @Models.evaluate_ODEexp; %applicable if using analytical solution

%initial parameter values
P1_0 = 0.5;
P2_0 = 1-P1_0;
k1 = 0.2;
k2 = 0.2;
k_s = 0.2;
alpha_s = [0.25,0.5,0.75];

%derived parameter values
[~,t_starts,t_ends] = dosefunc(1);
t_d1 = t_starts(1); % time for first dose
t_d1span = t_ends(1) - t_starts(1); % duration of first dose
P1_d1 = P1_0*exp((k1)*t_d1);
P2_d1 = P2_0*exp((k2)*t_d1);
Ps1_d1 = P1_0*exp(k_s*t_d1);
Ps2_d1 = P2_0*exp(k_s*t_d1);
plotSettings = struct('save', saveFigs, ...
                      'plotsDir', plotsDir, ...
                      'kineticsInd', plotKineticsInd, ...
                      'kineticsOverall', plotKineticsOverall);
y0 = [P1_0 P2_0]; %initial conditions

%initializations
ratio_alpha2alpha1Set = {};
slopesLists = {};
finalTumorSizeLists = {};
amtTumorReductionLists = {};
solvedEqn = '';
if(useAnalytical)
    solvedEqn = solvefunc(odefunc, dosefunc, y0)
end

%default settings/initializations
markercolors = [61/255 110/255 217/255; 34/255 73/255 156/255; 14/255 42/255 102/255];

for mm = 1:length(alpha_s)

    %finding starting beta values that will give real alpha values
    min_alpha1 = k1 - (1/t_d1span)*log( (1/P1_d1)*(Ps1_d1+Ps2_d1)*exp((k_s - alpha_s(mm))*t_d1span) );
    if min_alpha1 <= 0, min_alpha1 = 0; end

    %lists of alpha1 and alpha2 values for asymmetric treatment
    starts = min_alpha1 + 0.001;
    ends = alpha_s(mm);
    alpha1List = 2.^(log2(starts):log2(ends/starts)/30:log2(ends));
    alpha2List = k2 - (1/t_d1span)*log( (1/P2_d1)*((Ps1_d1+Ps2_d1)*exp((k_s - alpha_s(mm))*t_d1span) - ...
                                                    P1_d1*exp((k1 - alpha1List)*t_d1span)) );
    ratio_alpha2alpha1 = alpha2List./alpha1List;
    ratio_alpha2alpha1Set{mm} = ratio_alpha2alpha1;
    
    %collect all metrics for asymmetric treatment
    slopesList = [];
    finalTumorSizeList = [];
    amtTumorReductionList = [];

    for nn = 1:length(alpha1List)
        alpha1 = alpha1List(nn);
        alpha2 = alpha2List(nn);
        
        params = [k1 k2 alpha1 alpha2];
        results = Sim.simulate(odefunc, dosefunc, y0, params, solvedEqn, evalfunc, plotSettings);

        slopesList(nn) = results.slope;
        finalTumorSizeList(nn) = results.finalTumorSize;
        amtTumorReductionList(nn) = results.amtTumorReduction;

    end

    slopesLists{mm} = slopesList;
    finalTumorSizeLists{mm} = finalTumorSizeList;
    amtTumorReductionLists{mm} = amtTumorReductionList;
end

%generate plots
if plotRatios
    legend_txt = cellstr([repmat('overall kill = ',length(alpha_s),1), ...
                        num2str(alpha_s')]);
                    
    figure;
    hold on;
    for ii = 1:length(ratio_alpha2alpha1Set)
        plot(log2(ratio_alpha2alpha1Set{ii}),slopesLists{ii}, 'o', ...
            'markersize', 12, 'markerfacecolor', markercolors(ii, :), 'color','k')
    end
    set(gca,'fontsize',15)
    xlabel('log_{2}(\alpha/\beta)','fontsize',15);
    ylabel('Rate of change in tumor reduction','fontsize',15);
    title('Effects of alpha_s on [metrics] vs. alpha2/alpha1');
    legend(legend_txt, 'location', 'best');
    hold off;
    if(saveFigs), saveas(gcf, ['./plots_ratio_99/' num2str(S0) '_slopes.jpg']); end

    figure;
    hold on
    for ii = 1:length(ratio_alpha2alpha1Set)
        plot(log2(ratio_alpha2alpha1Set{ii}),log10(finalTumorSizeLists{ii}), 'o', ...
            'markersize', 12, 'markerfacecolor', markercolors(ii, :), 'color','k');
    end
    set(gca,'fontsize',15);
    xlabel('log_{2}(\alpha/\beta)','fontsize',15);
    ylabel('log_{10}(Final Tumor Size)','fontsize',15)
    legend(legend_txt, 'location', 'best');
    hold off
    if(saveFigs), saveas(gcf, ['./plots_ratio_99/' num2str(S0) '_finaltumorsize.jpg']); end

    figure;
    hold on
    for ii = 1:length(ratio_alpha2alpha1Set)
        plot(log2(ratio_alpha2alpha1Set{ii}),amtTumorReductionLists{ii}, 'o', ...
            'markersize', 12, 'markerfacecolor', markercolors(ii, :), 'color','k');
    end
    set(gca,'fontsize',15);
    xlabel('log_{2}(\alpha/\beta)','fontsize',15);
    ylabel('Final % Tumor Reduction','fontsize',15);
    legend(legend_txt, 'location', 'best');
    hold off
    if(saveFigs), saveas(gcf, ['./plots_ratio_99/' num2str(S0) '_amttumorreduct.jpg']); end
end