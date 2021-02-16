% One degree of freedom:
% Examine how [metric] varies as a function of differential drug selective pressure (alpha1/alpha2), and how this is
% affected by varying differential growth rate (k1/k2 ratio)
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
k_s = [0.14];
alpha_s = [0.5];

%derived parameter values
[~,t_starts,t_ends] = dosefunc(1);
t_d1 = t_starts(1); % time for first dose
t_d1span = t_ends(1) - t_starts(1); % duration of first dose
plotSettings = struct('save', saveFigs, ...
                      'plotsDir', plotsDir, ...
                      'kineticsInd', plotKineticsInd, ...
                      'kineticsOverall', plotKineticsOverall);
y0 = [P1_0 P2_0]; %initial conditions

%initializations
solvedEqn = '';
if(useAnalytical)
    solvedEqn = solvefunc(odefunc, dosefunc, y0)
end

%default settings/initializations
markercolors = [61/255 110/255 217/255];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulations
for ii = 1:length(k_s) %loop through each k_s
    %finding maximum allowable k1 value given constraint
    max_k1 = (1/t_d1) * log( (1/P1_0)*(P1_0+P2_0)*exp(k_s(ii)*t_d1) );
    
    %initialize data structures, for storing metric values from treatments
    alpha1Set = {};
    alpha2Set = {};
    k1Set = {};
    k2Set = {};
    ratio_k2k1Set = {};
    ratio_alpha2alpha1Set = {};
    slopesSet = {};
    finalTumorSizeSet = {};
    amtTumorReductionSet = {};
    
    %lists of k1 and k2 values for asymmetric treatment
    starts = 0.05;
    ends = k_s(ii);
    k1List = 2.^(log2(starts):log2(ends/starts)/5:log2(ends));
    %k1List = starts:(ends-starts)/3:ends;
    k2List = (1/t_d1) * log( (1/P2_0)*((P1_0+P2_0)*exp(k_s(ii)*t_d1) - P1_0*exp(k1List*t_d1)) );
    ratio_k2k1 = k2List./k1List;
    
    for jj = 1:length(k1List) %loop through each k2/k1 ratio
        %derived parameter values
        k1 = k1List(jj);
        k2 = k2List(jj);
        P1_d1 = P1_0*exp(k1*t_d1); %P1 under asymmetric growth k1 at time point d1
        P2_d1 = P2_0*exp(k2*t_d1); %P2 under asymmetric growth k2 at time point d1
        Ps1_d1 = P1_0*exp(k_s(ii)*t_d1); %P1 under symmetric growth ks at time point d1
        Ps2_d1 = P2_0*exp(k_s(ii)*t_d1); %P2 under symmetric growth ks at time point d1
        
        for mm = 1:length(alpha_s) %loop through alpha_s
            %finding minimum allowable alpha1 values given constraints (so alpha1 values are defined)
            min_alpha1 = k1 - (1/t_d1span)*log( (1/P1_d1)*(Ps1_d1+Ps2_d1)*exp((k_s(ii) - alpha_s(mm))*t_d1span) );
            if min_alpha1 <= 0, min_alpha1 = 0; end
            
            %lists of alpha1 and alpha2 values for asymmetric treatment
            starts = min_alpha1 + 0.001;
            ends = alpha_s(mm);
            alpha1List = 2.^(log2(starts):log2(ends/starts)/15:log2(ends));
            alpha2List = k2 - (1/t_d1span)*log( (1/P2_d1)*((Ps1_d1+Ps2_d1)*exp((k_s(ii) - alpha_s(mm))*t_d1span) - ...
                                                            P1_d1*exp((k1 - alpha1List)*t_d1span)) );
            ratio_alpha2alpha1 = alpha2List./alpha1List;
            
            for nn = 1:length(alpha1List) %loop through each alpha2/alpha1 ratio
                alpha1 = alpha1List(nn);
                alpha2 = alpha2List(nn);
                
                %simulation
                params = [k1 k2 alpha1 alpha2];
                results = Sim.simulate(odefunc, dosefunc, y0, params, solvedEqn, evalfunc, plotSettings);
                
                %collect results
                alpha1Set{jj,nn} = alpha1;
                alpha2Set{jj,nn} = alpha2;
                k1Set{jj,nn} = k1;
                k2Set{jj,nn} = k2;
                ratio_k2k1Set{jj,nn} = ratio_k2k1(jj);
                ratio_alpha2alpha1Set{jj,nn} = ratio_alpha2alpha1(nn);
                slopesSet{jj,nn} = results.slope;
                finalTumorSizeSet{jj,nn} = results.finalTumorSize;
                amtTumorReductionSet{jj,nn} = results.amtTumorReduction;
            end
            
            metricsSet = struct('alpha1', cell2mat(alpha1Set), ...
                                'alpha2', cell2mat(alpha2Set), ...
                                'k1', cell2mat(k1Set), ...
                                'k2', cell2mat(k2Set), ...
                                'ratio_k2k1', cell2mat(ratio_k2k1Set), ...
                                'ratio_alpha2alpha1', cell2mat(ratio_alpha2alpha1Set), ...
                                'slopes', cell2mat(slopesSet), ...
                                'fTumorSize', cell2mat(finalTumorSizeSet), ...
                                'amtTumorReduction', cell2mat(amtTumorReductionSet));
            metricsResults{ii,mm} = metricsSet;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
if plotRatios
    for ii = 1:length(k_s)
        for mm = 1:length(alpha_s)
            figure;
            hold on;
            k1k2ratiosN = length(metricsResults{ii,mm}.ratio_k2k1(:,1));
            for jj = 1:k1k2ratiosN
                %flip the values so the scatter marker of higher alpha2/alpha1 overlays the lower ones
                scatter(flip(log2(metricsResults{ii,mm}.ratio_alpha2alpha1(jj,:))), ...
                        flip(metricsResults{ii,mm}.slopes(jj,:)), ...
                        125, ...
                        'markerfacecolor', markercolors.*jj/k1k2ratiosN, ...
                        'markeredgecolor','k', ...
                        'linewidth', 1.25);
            end
            set(gca,'fontsize',15);
            xlabel('log_{2}(\alpha2/\alpha1)','fontsize',15);
            ylabel({'Drug sensitivity','(Rate of change in tumor reduction)'},'fontsize',15);
            title('Effects of k2/k1 on [metrics] vs. alpha2/alpha1');
            xmax = get(gca,'xlim');
            xmax = xmax(2);
            ymin = get(gca,'ylim');
            ymin = ymin(1);
            axis([0 xmax ymin 0.05]); %scale axis so minimum x is 0 and maximum y is 0.05
            hl = legend(cellstr([repmat('k2/k1 = ',k1k2ratiosN,1) ...
                              num2str(round(metricsResults{1,1}.ratio_k2k1(:,1)))]), 'location', 'best');
            set(findobj(hl,'type','patch'),'MarkerSize', sqrt(125));
            hold off;
            if(saveFigs), saveas(gcf, [plotsDir 'ks_' num2str(k_s(ii)) '-' 'alphas_' num2str(alpha_s(mm)) '.jpg']); end
            
            
            multiplier = (1:k1k2ratiosN)/k1k2ratiosN;
            [k1k2N n] = size(metricsResults{ii,mm}.slopes);
            markercolorsArray = markercolors'*[reshape([multiplier'*ones(1,n)]',1,n*length(multiplier))];
            markercolorsArray = markercolorsArray';
            
            %{
            figure;
            hold on;
            x = reshape(metricsResults{ii,mm}.ratio_k2k1', k1k2N*n, 1);
            y = reshape(log2(metricsResults{ii,mm}.ratio_alpha2alpha1)', k1k2N*n, 1);
            z = reshape(metricsResults{ii,mm}.slopes', k1k2N*n, 1);
            for i = 1:length(x) %had to be a loop because markerfacecolor cannot take a matrix
                scatter3(x(i), ...
                         y(i), ...
                         z(i), ...
                         125, ...
                         'markerfacecolor', markercolorsArray(i,:), ...
                         'markeredgecolor','k');
            end
            set(gca,'fontsize',10);
            ylabel('log_{2}(\alpha2/\alpha1)','fontsize',10);
            zlabel({'Drug sensitivity','(Rate of change in tumor reduction)'},'fontsize',10);
            xlabel('k2/k1','fontsize',10);
            title('Effects of k2/k1 on [metrics] vs. alpha2/alpha1','fontsize',10);
            grid on;
            view(3);
            hold off;
            if(saveFigs), saveas(gcf, [plotsDir 'ks_' num2str(k_s(ii)) '-' 'alphas_' num2str(alpha_s(mm)) '.jpg']); end
            %}
            
            figure;
            hold on;
            y = metricsResults{ii,mm}.ratio_k2k1';
            x = log2(metricsResults{ii,mm}.ratio_alpha2alpha1)';
            z = metricsResults{ii,mm}.slopes';
            surf(x,y,z);
            set(gca,'fontsize',10);
            xlabel('log_{2}(\alpha2/\alpha1)','fontsize',10);
            zlabel({'Drug sensitivity','(Rate of change in tumor reduction)'},'fontsize',10);
            ylabel('k2/k1','fontsize',10);
            title('Effects of k2/k1 on [metrics] vs. alpha2/alpha1','fontsize',10);
            grid on;
            view(2);
            shading faceted;
            %camlight right;
            colormap hot;
            view([37 32]);
            hold off;
            if(saveFigs), saveas(gcf, [plotsDir 'ks_' num2str(k_s(ii)) '-' 'alphas_' num2str(alpha_s(mm)) '.jpg']); end
            
        end
    end
end
 

