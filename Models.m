classdef Models < handle
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dose models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [d, starts, ends] = dose_exp(t)
            d = heaviside(t-8) - heaviside(t-12) + ...
                 heaviside(t-15) - heaviside(t-19) + ...
                 heaviside(t-22) - heaviside(t-26) + ...
                 heaviside(t-29) - heaviside(t-33);
            
             starts = [8,15,22,29];
             ends = [12,19,26,33];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kinetics models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Exponential growth model
        function dydt = ODEexp(t,y,param,dose)
            pcell = num2cell(param);
            [k1, k2, alpha1, alpha2] = pcell{:};

            d = dose(t);

            ycell = num2cell(y);
            [P1, P2]=ycell{:};

            dP1 = k1*P1-d.*alpha1.*P1;
            dP2 = k2*P2-d.*alpha2.*P2;

            dydt = [dP1;dP2];
        end
        
        function solvedEqn = solved_ODEexp(odefunc, dosefunc, y0)
            ycell = num2cell(y0);
            [P1_0, P2_0]=ycell{:};
            
            %note 'alpha' and 'beta' are reserved
            syms P1(t) P2(t) k1 k2 alpha1 alpha2;
            solvedEqn = dsolve(diff(P1)==k1*P1-alpha1*P1*(dosefunc(t)), ...
                               diff(P2)==k2*P2-alpha2*P2*(dosefunc(t)),...
                               P1(0)==P1_0, P2(0)==P2_0);
        end
        
        function [tspan, y] = evaluate_ODEexp(solvedEqn, tspan, param)
            pcell = num2cell(param);
            [k1_val, k2_val, alpha1_val, alpha2_val] = pcell{:};
            
            %assumes that tspan given is a horizontal vector (1xn), need to
            %change this to vertical vector (nx1)
            tspan = tspan';
            
            syms P1(t) P2(t) k1 k2 alpha1 alpha2 t;
            P1solved = eval(subs(solvedEqn.P1,{t,k1,k2,alpha1,alpha2},{tspan,k1_val,k2_val,alpha1_val,alpha2_val}));
            P2solved = eval(subs(solvedEqn.P2,{t,k1,k2,alpha1,alpha2},{tspan,k1_val,k2_val,alpha1_val,alpha2_val}));
            
            y = [P1solved, P2solved];
        end
        
        function params = getParams_ODEexp(dosefunc, ...
                                           P1_0_val, ...
                                           P2_0_val, ...
                                           ks_val, ...
                                           alphas_val, ...
                                           k1_val, ...
                                           k2_val, ...
                                           tolerance, ...
                                           verboseOutput)
            %Given P1_0, P2_0, alphas, ks, k1, and k2
            %calculate bounds for alpha1 and alpha2
            
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

            %derive parameters
            P1_d1_val = P1_0_val.*exp(k1_val.*t_d1);
            P2_d1_val = P2_0_val.*exp(k2_val.*t_d1);
            Ps1_d1_val = P1_0_val.*exp(ks_val.*t_d1);
            Ps2_d1_val = P2_0_val.*exp(ks_val.*t_d1);
            
            min_alpha1_theo = eval(subs(min_alpha1Eqn, ...
                                   {t, P1_d1, Ps1_d1, Ps2_d1, alphas, k1, ks}, ...
                                   {t_d1span, P1_d1_val, Ps1_d1_val, Ps2_d1_val, alphas_val, k1_val, ks_val}));
            min_alpha1 = min_alpha1_theo;
            min_alpha1(min_alpha1<0) = 0;
            min_alpha1 = min_alpha1 + tolerance;
            
            min_alpha2_theo = eval(subs(min_alpha2Eqn, ...
                                   {t, P2_d1, Ps1_d1, Ps2_d1, alphas, k2, ks}, ...
                                   {t_d1span, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alphas_val, k2_val, ks_val}));
            min_alpha2 = min_alpha2_theo;
            min_alpha2(min_alpha2<0) = 0;
            min_alpha2 = min_alpha2 + tolerance;
            
            max_alpha1 = eval(subs(treatSysGrowthEqn_alpha1, ...
                                   {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, alphas, k1, k2, ks, alpha2}, ...
                                   {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alphas_val, k1_val, k2_val, ks_val, min_alpha2}));
            max_alpha2 = eval(subs(treatSysGrowthEqn_alpha2, ...
                                   {t, P1_d1, P2_d1, Ps1_d1, Ps2_d1, alphas, k1, k2, ks, alpha1}, ...
                                   {t_d1span, P1_d1_val, P2_d1_val, Ps1_d1_val, Ps2_d1_val, alphas_val, k1_val, k2_val, ks_val, min_alpha1}));
            min_alpha2alpha1 = min_alpha2./max_alpha1;
            max_alpha2alpha1 = max_alpha2./min_alpha1;
            
            params = struct('P1_0', P1_0_val, ...
                            'P2_0', P2_0_val, ...
                            'alphas', alphas_val, ...
                            'ks', ks_val, ...
                            'k1', k1_val, ...
                            'k2', k2_val, ...
                            'P1_d1', P1_d1_val, ...
                            'P2_d1', P2_d1_val, ...
                            'Ps1_d1', Ps1_d1_val, ...
                            'Ps2_d1', Ps2_d1_val, ...
                            'min_alpha1', min_alpha1, ...
                            'min_alpha1_theo', min_alpha1_theo, ...
                            'max_alpha1', max_alpha1, ...
                            'min_alpha2', min_alpha2, ...
                            'min_alpha2_theo', min_alpha2_theo, ...
                            'max_alpha2', max_alpha2, ...
                            'min_alpha2alpha1', min_alpha2alpha1, ...
                            'max_alpha2alpha1', max_alpha2alpha1);
                        
            if(verboseOutput)
                %summary: parameter values
                disp('**** initial parameter values ****');
                disp(sprintf('P1_0: %0.2f', P1_0_val));
                disp(sprintf('P2_0: %0.2f', P2_0_val));
                disp(sprintf('alphas: %0.2f', alphas_val));
                disp(sprintf('ks: %0.2f', ks_val));
                disp(sprintf('k1: %0.2f', k1_val));
                disp(sprintf('k2: %0.2f', k2_val));
                disp(sprintf('Tolerance: %0.6f', tolerance));
                
                disp(' ');
                disp('**** given the provided k1 and k2 ****');
                disp(sprintf('P1_d1: %0.2f', P2_d1_val));
                disp(sprintf('P2_d1: %0.2f', P2_d1_val));
                disp(sprintf('Ps1_d1: %0.2f', Ps1_d1_val));

                %summary: parameter bounds
                disp(sprintf('Theoretical minimum alpha1: %0.2f', min_alpha1_theo));
                disp(sprintf('Theoretical minimum alpha2: %0.2f', min_alpha2_theo));
                disp(sprintf('Minimum alpha1: %0.2f', min_alpha1));
                disp(sprintf('Minimum alpha2: %0.2f', min_alpha2));
                disp(sprintf('Maximum alpha1: %0.2f', max_alpha1));
                disp(sprintf('Maximum alpha2: %0.2f', max_alpha2));
                disp(sprintf('Minimum alpha2/alpha1: %0.2f', min_alpha2alpha1));
                disp(sprintf('Maximum alpha2/alpha1: %0.2f', max_alpha2alpha1));
            end
        end
        
        function solutions_ODEexp() %for reference use only
            disp('***** symmetric growth, solve for k2 ***** ')
            syms P1_0 P2_0 k1 k2 ks t
            solve(P1_0*exp(k1*t) + P2_0*exp(k2*t) == (P1_0+P2_0)*exp(ks*t), k2)
            %solution: log((exp(ks*t)*(P1_0 + P2_0) - P1_0*exp(k1*t))/P2_0)/t

            disp('***** maximum k1 allowed ***** ')
            syms P1_0 P2_0 k1 ks t
            solve((exp(ks*t)*(P1_0 + P2_0) - P1_0*exp(k1*t))/P2_0 == 0, k1)
            %solution: log((exp(ks*t)*(P1_0 + P2_0))/P1_0)/t

            disp('***** symmetric treatment (with symmetric growth), solve for alpha2 ***** ')
            syms P1_d1 P2_d1 Ps1_d1 Ps2_d1 k1 alpha1 k2 alpha2 ks alphas t
            solve(P1_d1*exp(k1*t-alpha1*t) + P2_d1*exp(k2*t-alpha2*t) == (Ps1_d1+Ps2_d1)*exp(ks*t-alphas*t), alpha2)
            %solution: -(log((exp(ks*t - alphas*t)*(Ps1_d1 + Ps2_d1) - P1_d1*exp(k1*t - alpha1*t))/P2_d1) - k2*t)/t
            
            disp('***** minimum alpha1 (with symmetric growth) allowed ***** ')
            syms P1_d1 P2_d1 Ps1_d1 Ps2_d1 k1 alpha1 alpha2 ks alphas t
            solve((exp(ks*t - alphas*t)*(Ps1_d1 + Ps2_d1) - P1_d1*exp(k1*t - alpha1*t))/P2_d1 == 0, alpha1)
            %solution: -(log((exp(ks*t - alphas*t)*(Ps1_d1 + Ps2_d1))/P1_d1) - k1*t)/t
            
            disp('***** symmetric treatment, solve for alpha2 ***** ')
            syms P1_d1 P2_d1 k1 alpha1 k2 alpha2 alphas t
            solve(P1_d1*exp(k1*t-alpha1*t) + P2_d1*exp(k2*t-alpha2*t) == P1_d1*exp(k1*t-alphas*t)+P2_d1*exp(k2*t-alphas*t), alpha2)
            %solution: -(log((P1_d1*exp(k1*t - alphas*t) - P1_d1*exp(k1*t - alpha1*t) + P2_d1*exp(k2*t - alphas*t))/P2_d1) - k2*t)/t
            
            disp('***** minimum alpha1 allowed ***** ')
            syms P1_d1 P2_d1 k1 alpha1 alpha2 k1 k2 alphas t
            solve((P1_d1*exp(k1*t - alphas*t) - P1_d1*exp(k1*t - alpha1*t) + P2_d1*exp(k2*t - alphas*t))/P2_d1 == 0, alpha1)
            %solution: -(log((P1_d1*exp(k1*t - alphas*t) + P2_d1*exp(k2*t - alphas*t))/P1_d1) - k1*t)/t
        end

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Gompertz growth model
        function dydt = ODE_Gomp(t,y,param,dose)
            pcell = num2cell(param);
            [alpha, beta] = pcell{:};

            d = dose(t);

            ycell = num2cell(y);
            [S, R]=ycell{:};

            dS = .15.*log(1./S).*S-d.*alpha.*S;
            dR = .15.*log(1./R).*R-d.*beta.*R;

            dydt = [dS;dR];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Logistic growth model
        function dydt = ODE_logistic(t,y,param,dose)
            pcell = num2cell(param);
            [alpha, beta] = pcell{:};

            d = dose(t);

            ycell = num2cell(y);
            [S, R]=ycell{:};

            dS = S.*(1-S)-d.*alpha.*S;
            dR = R.*(1-R)-d.*beta.*R;

            dydt = [dS;dR];
        end
    end
end
