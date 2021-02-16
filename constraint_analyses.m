% Constraint analyses
close all;clear all;clc;

%settings
saveFigs = false;
plotsDir = './plots/';
plotFigs = true;
odefunc = @Models.ODEexp;
dosefunc = @Models.dose_exp;
useAnalytical = false;
solvefunc = @Models.solved_ODEexp; %applicable if using analytical solution
evalfunc = @Models.evaluate_ODEexp; %applicable if using analytical solution

%initial parameter values
P1_0_val = 0.5;
P2_0_val = 1-P1_0_val;
ks_val = 0.14;
alphas_val = 0.5;

%derived parameter values
syms P1_0 P2_0 P1_d1 P2_d1 Ps1_d1 Ps2_d1 k1 k2 ks alpha1 alpha2 alphas t;
[~,t_starts,t_ends] = dosefunc(1);
t_d1 = t_starts(1); % time for first dose
t_d1span = t_ends(1) - t_starts(1); % duration of first dose
y0 = [P1_0 P2_0]; %initial conditions

%equations
symGrowthEqn = P1_0*exp(k1*t) + P2_0*exp(k2*t) == (P1_0+P2_0)*exp(ks*t);
treatSysGrowthEqn = P1_d1*exp(k1*t-alpha1*t) + P2_d1*exp(k2*t-alpha2*t) == (Ps1_d1+Ps2_d1)*exp(ks*t-alphas*t);

%solved equations for specific variables
symGrowthEqn_k2 = solve(symGrowthEqn, k2);
symGrowthEqn_ks = solve(symGrowthEqn, ks);
treatSysGrowthEqn_alpha1 = solve(treatSysGrowthEqn, alpha1);
treatSysGrowthEqn_alpha2 = solve(treatSysGrowthEqn, alpha2);
treatSysGrowthEqn_ks = solve(treatSysGrowthEqn, ks);
treatSysGrowthEqn_alphas = solve(treatSysGrowthEqn, alphas);

%solved equations, for max/min by setting expression inside ln equal to zero
max_k1Eqn = log((exp(ks*t)*(P1_0 + P2_0))/P1_0)/t;
min_alpha1Eqn = -(log((exp(ks*t - alphas*t)*(Ps1_d1 + Ps2_d1))/P1_d1) - k1*t)/t;
min_alpha2Eqn = -(log((exp(ks*t - alphas*t)*(Ps1_d1 + Ps2_d1))/P2_d1) - k2*t)/t;

%additional derived parameter values, based on k1=k2=ks
k1_val = ks_val;
k2_val = ks_val;
params = Models.getParams_ODEexp(dosefunc, P1_0_val, P2_0_val, ks_val, alphas_val, k1_val, k2_val, 0.001, true);

P1_d1_val = params.P1_d1;
P2_d1_val = params.P2_d1;
Ps1_d1_val = params.Ps1_d1;
Ps2_d1_val = params.Ps2_d1;
min_alpha1 = params.min_alpha1;
min_alpha2 = params.min_alpha2;
min_alpha1_theo = params.min_alpha1_theo;
min_alpha2_theo = params.min_alpha2_theo;
max_alpha1 = params.max_alpha1;
max_alpha2 = params.max_alpha2;
min_alpha2alpha1 = params.min_alpha2alpha1;
max_alpha2alpha1 = params.max_alpha2alpha1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examines relationship between k1 and k2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% k2 versus k1 plots with single ks value; given P1_0, P2_0, ks, and dose function
figure;
k1List = 0:0.01:ks_val;
k2List = eval(subs(symGrowthEqn_k2, {t, P1_0, P2_0, ks, k1}, {t_d1, P1_0_val, P2_0_val, ks_val, k1List}));
plot(k1List, k2List, 'k', 'linewidth', 2);
xlabel('k1');
ylabel('k2');
title('k2 versus k1; given P_1^0, P_2^0, k_s, and dose function');

%% k2 versus k1 plots with varying ks; given P1_0, P2_0, and dose function
ksEqn_subs = eval(subs(symGrowthEqn_ks, {t, P1_0, P2_0}, {t_d1, P1_0_val, P2_0_val}));
[k1List, k2List] = meshgrid(0:0.01:ks_val, 0:0.01:ks_val);
ksList = eval(subs(ksEqn_subs,{k1,k2},{k1List,k2List}));

%{
figure;
%ezsurf(ksEqn_subs,[0 ks_val],6); %resolution for ez-based plots are low
surf(k1List, k2List, ksList);
shading interp;
view(2);
xlabel('k1');
ylabel('k2');
zlabel('ks');
h = colorbar;
ylabel(h,'ks');
%}

figure;
%ezcontour(ksEqn_subs,[0 ks_val],6); %resolution for ez-based plots are low
contour(k1List, k2List, ksList, 6, 'ShowText', 'on', 'linewidth', 2);
xlabel('k1');
ylabel('k2');
zlabel('ks');
h = colorbar;
ylabel(h,'ks');
title('k2 versus k1 with varying k_s; given P_1^0, P_2^0, and dose function');

%% maximum allowable k1
figure;
hold on;
ksList = 0:0.01:0.14;
max_k1List = eval(subs(max_k1Eqn, {t, P1_0, P2_0, ks}, {t_d1, P1_0_val, P2_0_val, ksList}));
plot(ksList, max_k1List, 'k', 'linewidth', 2);
h = refline([1 0]);
set(h,'color',[0.7 0.7 0.7]);
set(h,'linestyle',':');
xlabel('ks');
ylabel('maximum k1 allowable');
title('maximum k1 versus varying k_s');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examines relationship between alpha1 and alpha2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alpha2 versus alpha1 plots with single alphas value; given P1_0, P2_0, ks, alphas, and dose function
% examined at symmetric growth with k1=k2=ks
figure;
hold on;
alpha1List = min_alpha1:0.005:alphas_val;
alpha2List = eval(subs(treatSysGrowthEqn_alpha2, ...
                       {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, alpha1, alphas, k1, k2, ks}, ...
                       {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alpha1List, alphas_val, k1_val, k2_val, ks_val}));
plot(alpha1List, alpha2List, 'k', 'linewidth', 2);
h = line([min_alpha1_theo min_alpha1_theo], [0 max(get(gca,'ylim'))],'color',[0.7 0.7 0.7],'linestyle',':');
h = line(get(gca,'xlim'), [min_alpha2_theo min_alpha2_theo],'color',[0.7 0.7 0.7],'linestyle',':');
xlabel('\alpha_1');
ylabel('\alpha_2');
title('\alpha_1 versus \alpha_2; given P_1^0, P_2^0, k_s, \alpha_s, and dose function');
hold off;

%% alpha2 versus alpha1 plots with varying ks; given P1_0, P2_0, alphas, and dose function
ksEqn_subs = eval(subs(treatSysGrowthEqn_ks, ...
                       {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, alphas, k1, k2}, ...
                       {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alphas_val, k1_val, k2_val}));
[alpha1List, alpha2List] = meshgrid(0:0.01:alphas_val, 0:0.01:alphas_val);
ksList = eval(subs(ksEqn_subs,{alpha1,alpha2},{alpha1List,alpha2List}));

%{
figure;
surf(alpha1List, alpha2List, ksList);
shading interp;
view(2);
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('ks');
h = colorbar;
ylabel(h,'ks');
title('\alpha_1 versus \alpha_2 with varying k_s; given P_1^0, P_2^0, \alpha_s, and dose function');
%}

figure;
h = contour(alpha1List, alpha2List, ksList, 6, 'ShowText', 'on', 'linewidth', 2);
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('ks');
h = colorbar;
ylabel(h,'ks');
title('\alpha_1 versus \alpha_2 with varying k_s; given P_1^0, P_2^0, \alpha_s, and dose function');

%% alpha2 versus alpha1 plots with varying alphas; given P1_0, P2_0, k2, and dose function
alphasEqn_subs = eval(subs(treatSysGrowthEqn_alphas, ...
                       {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, ks, k1, k2}, ...
                       {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, ks_val, k1_val, k2_val}));
[alpha1List, alpha2List] = meshgrid(0:0.01:alphas_val, 0:0.01:alphas_val);
alphasList = eval(subs(alphasEqn_subs,{alpha1,alpha2},{alpha1List,alpha2List}));

%{
figure;
surf(alpha1List, alpha2List, alphasList);
shading interp;
view(2);
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('\alpha_s');
h = colorbar;
ylabel(h,'\alpha_s');
title('\alpha_1 versus \alpha_2 with varying \alpha_s; given P_1^0, P_2^0, k_s, and dose function');
%}

figure;
contour(alpha1List, alpha2List, alphasList, 6, 'ShowText', 'on', 'linewidth', 2);
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('\alpha_s');
h = colorbar;
ylabel(h,'\alpha_s');
title('\alpha_1 versus \alpha_2 with varying \alpha_s; given P_1^0, P_2^0, k_s, and dose function');

%% minimum allowable alpha1
figure;
hold on;
alphasList = 0:0.01:0.9;
min_alpha1List = eval(subs(min_alpha1Eqn, ...
                           {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, alphas, k1, ks}, ...
                           {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alphasList, k1_val, ks_val}));
plot(alphasList, min_alpha1List, 'k', 'linewidth', 2);
h = refline([0 0]);
set(h,'color',[0.7 0.7 0.7]);
set(h,'linestyle',':');
xlabel('\alpha_s');
ylabel('minimum \alpha_1 allowable');
title('minimum \alpha_1 versus varying \alpha_s');
hold off;

%% maxmimum log2(alpha2/alpha2)
starts = 0.001;
ends = ks_val;
k1List = 2.^(log2(starts):log2(ends/starts)/10:log2(ends));
k2List = eval(subs(symGrowthEqn_k2, {t, P1_0, P2_0, ks, k1}, {t_d1, P1_0_val, P2_0_val, ks_val, k1List}));
k2k1RatioList = k2List./k1List;
k1_val = k1List;
k2_val = k2List;

markercolors = [61/255 110/255 217/255];
tolList = logspace(-1,-6,6); %tolerance, the value to be added to the minimum alpha1
tolListN = length(tolList);

figure;
hold on;
for idx = 1:tolListN
    params = Models.getParams_ODEexp(dosefunc, P1_0_val, P2_0_val, ks_val, alphas_val, k1_val, k2_val, tolList(idx), false);
    %{
    scatter(params.min_alpha1, log2(params.max_alpha2alpha1), 100, ...
            'markeredgecolor','k', ...
            'markerfacecolor', markercolors.*idx/tolListN);
    %}
    plot(params.min_alpha1, log2(params.max_alpha2alpha1), '--o', 'color', markercolors.*idx/tolListN, 'linewidth', 2, 'markersize', 10);
end
legend(cellstr([repmat('tol = ',length(tolList),1) num2str(tolList')]), 'location', 'best');
xlabel('Minimum \alpha1');
ylabel('Maximum log2(\alpha2/\alpha1)');
hold off;

figure;
hold on;
for idx = 1:length(tolList)
    params = Models.getParams_ODEexp(dosefunc, P1_0_val, P2_0_val, ks_val, alphas_val, k1_val, k2_val, tolList(idx), false);
    %{
    scatter(log2(params.max_alpha2alpha1), k2k1RatioList, 100, ...
            'markeredgecolor','k', ...
            'markerfacecolor', markercolors.*idx/tolListN);
    %}
    plot(log2(params.max_alpha2alpha1), k2k1RatioList, '--o', 'color', markercolors.*idx/tolListN, 'linewidth', 2, 'markersize', 10);
end
legend(cellstr([repmat('tol = ',length(tolList),1) num2str(tolList')]), 'location', 'best');
xlabel('Maximum log2(\alpha2/\alpha1)');
ylabel('k2/k1');
hold off;

figure;
hold on;
for idx = 1:length(tolList)
    params = Models.getParams_ODEexp(dosefunc, P1_0_val, P2_0_val, ks_val, alphas_val, k1_val, k2_val, tolList(idx), false);
    %{
    scatter(log2(params.max_alpha2alpha1), k2k1RatioList, 100, ...
            'markeredgecolor','k', ...
            'markerfacecolor', markercolors.*idx/tolListN);
    %}
    plot(log2(params.min_alpha2alpha1), k2k1RatioList, '--o', 'color', markercolors.*idx/tolListN, 'linewidth', 2, 'markersize', 10);
end
legend(cellstr([repmat('tol = ',length(tolList),1) num2str(tolList')]), 'location', 'best');
xlabel('Minimum log2(\alpha2/\alpha1)');
ylabel('k2/k1');
hold off;

