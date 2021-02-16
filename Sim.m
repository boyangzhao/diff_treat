classdef Sim < handle
    
    methods(Static)
        function m = getMax(x)
            m = max(max(x));
            x = round(x*10)/10.0;
        end
        
        function [starts,ends] = getTreatEnds(x)
            %retrieve treatment cycle start and end indices given a course of dose
            %treatment; only use this if the dosefunc does not supply this information
            %input: treatment cycle (binary vector)
            %output: vector of last time point for each treatment cycle (values refer to indices in x)
            starts = [];
            ends = [];
            t = find(x == 1);
            starts = t(1);
            if(length(t) > 1)
                for i = 2:length(t)             
                    diff = t(i)-t(i-1);
                    if(diff > 1)
                        ends = [ends;t(i-1)];
                        starts = [starts;t(i)];
                    end
                end
            end
            if(ends(end) ~= t(end))
                ends = [ends;t(end)];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulation methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function results = simulate(odefunc, dosefunc, y0, params, solvedEqn, evaluatefunc, plotSettings)
            if nargin < 7
                plotSettings = struct('save', false, ...
                                      'plotsDir', './', ...
                                      'kineticsInd', false, ...
                                      'kineticsOverall', false);
            end
            
            tspan = 0:0.001:35;

            %indices of treatment cycles
            %[~,starts,ends] = dosefunc(tspan); 
            [starts,ends] = Sim.getTreatEnds(dosefunc(tspan));
            
            if(~isempty(solvedEqn))
                %solve analytically
                [t,y] = evaluatefunc(solvedEqn, tspan, params);
            else
                %solve numerically
                options = odeset('RelTol',1e-8,'AbsTol',[1e-8 1e-8]);
                [t,y] = ode45(odefunc, tspan, y0, options, params, dosefunc);
            end
            if(plotSettings.kineticsInd), Sim.plotKineticsInd(t, y, dosefunc, params, plotSettings); end

            %calculate overall tumor metrics
            y_sum = sum(y,2);
            finalTumorSize = y_sum(end);
            amtreduction = (y_sum(starts)-y_sum(ends))./y_sum(starts);
            amtTumorReduction = amtreduction(end);
            fitfunc = @(k,x) k(1).*x+k(2);
            kfit = nlinfit(1:length(amtreduction),amtreduction',fitfunc,[1 0]);
            slope = kfit(1);
            if(plotSettings.kineticsOverall), Sim.plotKineticsOverall(t, y, y_sum, amtreduction, kfit, params, plotSettings); end
            
            results = struct('slope', slope, ...
                             'finalTumorSize', finalTumorSize, ...
                             'amtTumorReduction', amtTumorReduction);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotKineticsOverall(t, y, y_sum, amtreduction, kfit, params, plotSettings)
            figure(1);
            plot(t, y_sum, 'k--','linewidth', 2);
            set(gca,'fontsize', 16, 'YLim',[0 Sim.getMax([sum(y,2)])]);
            xlabel('Time', 'fontsize', 16);
            ylabel('Overall Population Size (a.u.)', 'fontsize', 16);
            title(['\alpha_1=' num2str(params(1,3)) '; \alpha_2=' num2str(params(1,4)) ...
                 '; \alpha_2/\alpha_1 =' num2str(params(1,4)/params(1,3)) ';k1=' num2str(params(1,1)) '; k2=' num2str(params(1,2))...
                   'k2/k2= ' num2str(params(1,2)/params(1,2))], 'fontsize', 10)
            if(plotSettings.save), saveas(gcf, [plotSettings.plotsDir 'alpharatio_' num2str(params(1,4)/params(1,3))...
                    '_kratio_' num2str(params(1,2)/params(1,1)) '_kineticsoverlay.jpg']); end

            %Endpoint tumor size
            figure;
            title('Endpoint tumor size comparison');
            bar(y_sum(end),0.8,'k');
            set(gca,'XTick', 1:2);
            set(gca,'XLim',[0.4 2.6]);
            ylabel('Final tumor size');
            if(plotSettings.save), saveas(gcf, [plotSettings.plotsDir num2str(uid) '_endpointturmorsize.jpg']); end

            %Amount of tumor reduction
            figure(2);
            bar(amtreduction','grouped','k');
            ylabel('% Tumor reduction', 'fontsize', 16);
            title(['%Tumor Reduction:' '\alpha_1=' num2str(params(1,3)) '; \alpha_2=' num2str(params(1,4)) ...
                 '; \alpha_2/\alpha_1 =' num2str(params(1,4)/params(1,3)) ';k1=' num2str(params(1,1)) '; k2=' num2str(params(1,2))...
                   'k2/k2= ' num2str(params(1,2)/params(1,2))], 'fontsize', 10)
            if(saveFigs), saveas(gcf, ['./plots_kinetics_tumorcomp_0.9/' 'alpharatio_' num2str(params(1,4)/params(1,3))...
                    '_kratio_' num2str(params(1,2)/params(1,1)) '_tumorreduction.jpg']); end

            Rate of adaptation
            figure;
            bar(kfit(1),0.8,'k');
            ylabel('Rate of change in tumor reduction', 'fontsize', 16);
            title(['Rate of Change in Tumor Reduction:' '\alpha=' num2str(params(1,3)) '; \beta=' num2str(params(1,4)) ...
            '; \alpha/\beta =' num2str(params(1,3)/params(1,4))],'fontsize', 14)
            set(gca, 'fontsize', 16)
            if(plotSettings.save), saveas(gcf, [plotSettings.plotsDir num2str(params(1,3)) '_' num2str(params(1,4)) '_rateofadaptation.jpg']); end
        end
        
        function plotKineticsInd(t, y, dosefunc, params, plotSettings)
            %Drug dose
            figure(3);
            subplot(2,1,1);
            plot(t,dosefunc(t),'k','linewidth',2);
            title('Drug', 'fontsize', 14);
            set(gca,'YLim',[0 1.2], 'fontsize', 16);
            ylabel('Concentration', 'fontsize', 16);
            xlabel('Time', 'fontsize', 16);

            %Subpopulation plot
            subplot(2,1,2);
            y_max = Sim.getMax(y);
            plot(t, y(:,1), 'b', t, y(:,2), 'r', 'linewidth', 2);
            legend('sensitive','resistant');
            set(legend,'fontsize',14)
            set(gca,'YLim',[0 y_max], 'fontsize', 16);
            title(['Subpopulations:' '\alpha_1=' num2str(params(1,3)) '; \alpha_2=' num2str(params(1,4)) ...
                 '; \alpha_2/\alpha_1 =' num2str(params(1,4)/params(1,3)) ';k1=' num2str(params(1,1)) '; k2=' num2str(params(1,2))...
                   'k2/k2= ' num2str(params(1,2)/params(1,2))], 'fontsize', 10)
            ylabel('Subpopulation size', 'fontsize', 16);
            xlabel('Time', 'fontsize', 16);
            if(plotSettings.save), saveas(gcf, [plotSettings.plotsDir 'alpharatio_' num2str(params(1,4)/params(1,3))...
                    '_kratio_' num2str(params(1,2)/params(1,1)) '_dosesubpopulationkinetics.jpg']); end
            
            %Tumor composition
            figure(4);
            subplot(1,2,1);
            pie(y(1,:)./sum(y(1,:)));
            title('Initial tumor composition', 'fontsize', 12);
            colormap([50/255 108/255 191/255;232/255 44/255 12/255]);

            subplot(1,2,2);
            pie(y(end,:)./sum(y(end,:)));
            title(['Final tumor composition' '\alpha_1=' num2str(params(1,3)) '; \alpha_2=' num2str(params(1,4)) ...
                 '; \alpha_2/\alpha_1 =' num2str(params(1,4)/params(1,3)) ';k1=' num2str(params(1,1)) '; k2=' num2str(params(1,2))...
                   'k2/k2= ' num2str(params(1,2)/params(1,2))], 'fontsize', 10);
            colormap([50/255 108/255 191/255;232/255 44/255 12/255]);
            if(plotSettings.save), saveas(gcf, [plotSettings.plotsDir 'alpharatio_' num2str(params(1,4)/params(1,3))...
                    '_kratio_' num2str(params(1,2)/params(1,1)) '_tumorcomposition.jpg']); end
        end
    end
end