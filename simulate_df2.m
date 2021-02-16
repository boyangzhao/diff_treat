% Two degrees of freedom: with floating alpha and beta
close all;clear all;clc;

%settings
saveFigs = false;
plotsDir = './plots/';
plotKineticsInd = false;
plotKineticsOverall = false;
odefunc = @Models.ODEexp;
dosefunc = @Models.dose_exp;
useAnalytical = true;
solvefunc = @Models.solved_ODEexp; %applicable if using analytical solution
evalfunc = @Models.evaluate_ODEexp; %applicable if using analytical solution
symAlphaBeta = 0.25;

%initial parameter values
S0 = 0.5;
R0 = 0.5;
k1a = 0.14;
k2a = 0.14;

%derived parameters
alphabetaf = 0:0.025:symAlphaBeta;
[alphaf, betaf] = meshgrid(alphabetaf);
plotSettings = struct('save', saveFigs, ...
                      'plotsDir', plotsDir, ...
                      'kineticsInd', plotKineticsInd, ...
                      'kineticsOverall', plotKineticsOverall);
y0 = [S0 R0]; %initial conditions

%initializations
slopesMatrix = zeros(length(alphabetaf),length(alphabetaf));
finalTumorSizeMatrix = zeros(length(alphabetaf),length(alphabetaf));
amtTumorReductionMatrix = zeros(length(alphabetaf),length(alphabetaf));
solvedEqn = '';
if(useAnalytical)
    solvedEqn = solvefunc(odefunc, dosefunc, y0);
end

for ii = 1:length(alphabetaf)
    for jj = 1:length(alphabetaf)
   
        beta = betaf(ii,jj);
        alpha = alphaf(ii,jj);
        
        params = [k1a k2a alpha beta];
        results = Sim.simulate(odefunc, dosefunc, y0, params, solvedEqn, evalfunc, plotSettings);
    
        slopesMatrix(ii,jj) = results.slope;
        finalTumorSizeMatrix(ii,jj) = results.finalTumorSize;
        amtTumorReductionMatrix(ii,jj) = results.amtTumorReduction;
        
    end
end
    
%for plotting alpha dependent on beta
%initial parameter values
d = 1;
S8a = 1.5324;
S8s = 1.5324;
R8a = 1.5324;
R8s = 1.5324;

%symmetric parameters
k1a = 0.14;
k2a = 0.14;

%asymmetric parameters
k1s = 0.14;
k2s = 0.14;

t = 4;

%generate lists of alpha and beta values for asymmetric treatment
minBeta = (k2a - log((R8s*exp(t*(k2s - symAlphaBeta*d)) + S8s*exp(t*(k1s - symAlphaBeta*d)))/R8a)/t)/d;
betaList = (minBeta+0.001):0.01:symAlphaBeta; 
alphaList = (k1a - log((R8s.*exp(t.*(k2s - symAlphaBeta.*d)) - R8a.*exp(t.*(k2a - betaList.*d)) + S8s.*exp(t.*(k1s - symAlphaBeta.*d)))./S8a)./t)./d;
ratioAlphaToBeta = alphaList./betaList;

%figures: contour and surface plots for 2 degrees of freedom 
%plots for final tumor size
figure;
surf(alphaf,betaf,finalTumorSizeMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Final Tumor Size', 'fontsize', 15)
colorbar('location','EastOutside')
view(2)
shading interp
if(saveFigs), saveas(gcf,'float_finaltumorsize_surf.jpg'); end

figure;
plot(alphaList,betaList,'k--','linewidth',2)
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Final Tumor Size', 'fontsize', 15)
if(saveFigs), saveas(gcf,'float_finaltumorsize_alphabeta.jpg'); end

figure;
contour(alphaf,betaf,finalTumorSizeMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Final Tumor Size', 'fontsize', 15)
colorbar('location','EastOutside')
if(saveFigs), saveas(gcf,'float_finaltumorsize_contour.jpg'); end

%plots for final percent tumor reduction
figure;
surf(alphaf,betaf,amtTumorReductionMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Final Percent Tumor Reduction', 'fontsize', 15)
colorbar('location','EastOutside')
view(2)
shading interp
if(saveFigs), saveas(gcf,'float_finalpercenttumorred_surf.jpg'); end

figure; 
plot(alphaList,betaList,'k','linewidth', 3)
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Final Percent Tumor Reduction', 'fontsize', 15)
if(saveFigs), saveas(gcf,'float_finalpercenttumorred_alphabeta.jpg'); end

figure;
contour(alphaf,betaf,amtTumorReductionMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha')
ylabel('\beta')
title('Final Percent Tumor Reduction', 'fontsize', 15)
colorbar('location','EastOutside')
if(saveFigs), saveas(gcf,'float_finalpercenttumorred_contour.jpg'); end

%plot rate of change of tumor size
figure;
surf(alphaf,betaf,slopesMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha')
ylabel('\beta')
title('Rate of Change of Tumor Size', 'fontsize', 15)
colorbar('location','EastOutside')
view(2)
shading interp
if(saveFigs), saveas(gcf,'float_ratechangetumorsize_surf.jpg'); end

figure; 
plot(alphaList,betaList,'k','linewidth', 3)
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Rate of Change of Tumor Size', 'fontsize', 15)
if(saveFigs), saveas(gcf,'float_ratechangetumorsize_alphabeta.jpg'); end

figure;
contour(alphaf,betaf,slopesMatrix)
set(gca,'xlim',[alphabetaf(1) alphabetaf(end)])
set(gca,'fontsize',15)
xlabel('\alpha', 'fontsize', 15)
ylabel('\beta', 'fontsize', 15)
title('Rate of Change of Tumor Size', 'fontsize', 15)
colorbar('location','EastOutside')
if(saveFigs), saveas(gcf,'float_ratechangetumorsize_contour.jpg'); end
