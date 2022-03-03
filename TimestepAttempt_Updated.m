% Variables
% q_r, q_d, q_ev
q_r_d = [40];
for i = 1:length(q_r_d)
    experiment_number = i;
    q_d_max = q_r_d(i);
    q_r_max = q_r_d(i);

    tic
    % maximum energy
    q_cap = 350; % 350kW

    % maximum discharge power (grid tie)
    % q_d_max = 40;

    % maximum recharge power (grid tie)
    % q_r_max = 40;

    % discharge cost $/kWh
    C_d = 0;

    % recharge cost $/kwh
    C_r = 0;

    % discount rate
    discount = 0;

    % charging efficiency
    gamma_c = 0.87;

    % storage efficiency
    gamma_s = 0.87;

    % Regulation Factors: 
    mu_rd = 0.5;
    mu_ru = 0.5;

    alpha_ru = 0.5;
    alpha_rd = 0.5;

    gamma_ru = mu_ru * alpha_ru;

    gamma_rd = mu_rd * alpha_rd;

    % Resource Adequacy Penalty
    Penalty = 5000;


    % EV charging unit price
    EV_unitPrice = 0.55; % $/kWh


    %% Initialize
    disp('Reading Input Data')

    % Input data
    LMP_CA = readtable('2020 MTVIEW N001 LMP DATA.xlsx', 'Sheet', 'Data', 'ReadVariableNames', true, 'PreserveVariableNames', true);
    EV_demand = readtable('EV demand input.xlsx', 'Sheet', 'Demand', 'ReadVariableNames', true,'PreserveVariableNames', true);
    months = ["January", "February", "March", "April", "May", "June", "July",...
        "August", "September", "October", "November", "December"];


    ra_month = zeros(12,1);
    experiment_folder = sprintf("experiment_noev_number_%i_%i", experiment_number, q_r_d(i));
    status = mkdir(experiment_folder);
    filename = sprintf(experiment_folder + '/Optimal Solution_%i.xlsx',experiment_number);
    years_CA = max(LMP_CA{:, 'Year'});
    start = 1;
    disp('Beginning Variable Initialization')

    y = 1;
    % Call LMP f
    lastRow = find(LMP_CA{:, 'Year'} == y, 1, 'last');
    YearLMP = LMP_CA{start:lastRow, 7:end};
    YearLMP_vec = 1e-3* reshape(YearLMP.',1,[]);
    YearLMP_vec(1,1611) = (mean(YearLMP_vec(1,1609:1610))+ mean(YearLMP_vec(1,1612:1613)))/2; %filling in missing nan data
    %EV Discharge Minimum
%     EV_disch_min = repmat(transpose(EV_demand{:, 2}), 1, size(YearLMP, 1));
    % NO EV CHARGING
    EV_disch_min = repmat(zeros(length(transpose(EV_demand{:, 2}))), 1, size(YearLMP, 1));
    P_RA =[ 3.07 3.06 3.02 3.13 3.22 2.91 4.32 4.45 5.94 3.71 3.24 3.27];
    prevInd = 0;
    A_prev = 0;
    monthlyMaxRevenue = zeros(12,1);


    area_graph = zeros(length(YearLMP_vec), 6);
    area_graph_monthly = zeros(12,6);
    for m = 1:12
        disp('Month = ')
        disp(m)
        startInd = find(LMP_CA{:, 'Month'} == months(m), 1, 'first');
        endInd = find(LMP_CA{:, 'Month'} == months(m), 1, 'last');
        numDays = endInd - startInd + 1;
        vecLength = numDays * 24;
        f_start = prevInd + 1;
        f_end = f_start + vecLength - 1;
        prevInd = f_end;    
        reshapeLength = 4*vecLength; % for 4 (EV, ab, recharge, RA)

        A_start = A_prev + 1;
        A_end = A_start + reshapeLength - 1;
        A_prev = A_end;
        j = [];
        for i = 1:numDays
            j = [j, 16+24*(i-1), 17+24*(i-1), 18+24*(i-1), 19+24*(i-1), 20+24*(i-1)];
        end
        bin = zeros(vecLength,1);
        for i =1:vecLength
            if sum(j==i) > 0
                bin(i) = randi([0 1], 1,1);
            end
        end
        %Random Vector generates numDays x 1 array adding up to 1 for Resource Adequacy
        randVec = rand(sum(bin),1);
        ind = 1;
        gamma_ra = zeros(vecLength,1);
        for i = 1:vecLength
            if bin(i)==1
                gamma_ra(i) = randVec(ind);
                ind = ind + 1;
            end
        end
        disp("Beginning optimization")
    %     cvx_solver mosek
        cvx_begin
            variable q_d(vecLength);    %Discharge Energy (kWh)
            variable q_r(vecLength);    %Charging Energy (kWh)
            variable q_ev(vecLength);   %EV Energy (kWH)
            variable S_t(vecLength+1);  %State of Charge (kWh)
            variable q_ra;              %Resource Adequacy (kWh)
    %         variable cha(vecLength) binary;
    %         variable q_ru(vecLength);    %Regulation Up (kWh)
    %         variable q_rd(vecLength);    %Regulation Down (kWh)
    %         variable gamma_ra_sc(vecLength);
            minimize(-sum(EV_unitPrice .* q_ev) - YearLMP_vec(f_start:f_end) * q_d + YearLMP_vec(f_start:f_end) * q_r - P_RA(m) * q_ra - YearLMP_vec(f_start:f_end)*(gamma_ra .*q_ra));
            % The above is minimize (-revenue from EV Charging - revenue from
            % Arbitrage + cost of recharging - Revenue from RA + Penalty from
            % RA. 
            % We need to add RU/RD Terms.
            subject to
                for i =1:vecLength
                   % State of charge constraint for each timestep
                   S_t(i+1) == gamma_s * S_t(i) + gamma_c * q_r(i) - q_ev(i) - q_d(i) - gamma_ra(i)*q_ra; % + gamma_ru * q_ru + gamma_rd * q_rd;
                end
                %Initial SOC Constraint
                S_t(1) == 0;
                S_t(vecLength+1) == 0;
                % Charging Power Constraint
                0 <= q_r <= q_r_max; 
                %0 <= q_r + alpha_rd*gamma_rd*q_rd <= q_r_max (alternate)
                % General positivity constraints
                q_d >= 0;
                q_ra >= 0;
                S_t >= 0;
                S_t <= q_cap;
                %q_ru >= 0;
                %q_rd >= 0;
                % gamma_ra_sc should be summing to 1 (as it is a percentage of
                % RA called upon at time t in a month m. 

                for i = 1:vecLength
                    % Discharging Power Constraint
                    q_ev(i) + q_d(i) + gamma_ra(i) * q_ra <= q_d_max;
                    q_r(i) + q_d(i) + gamma_ra(i)*q_ra <= q_d_max;
                    %Alternate: q_ev(i) + q_d(i) + gamma_ra_sc(i)*q_ra + alpha_ru*gamma_ru*q_ru<= q_d_max;
                    % EV Charging Minimum Constraint
                    q_ev(i) == EV_disch_min(i);     
                    % If gamma_ra(i) >= gamma_ra_sc(i), then there is a lack of
                    % power and we have a penalty => RHS >= 0, t >= 0.
                end
         cvx_end

        ra_month(m) = q_ra;
        powerOpt.discharging = q_d;
        powerOpt.recharging = q_r;
        powerOpt.EVPower = q_ev;
        monthlyMaxRevenue(m) = cvx_optval;
        if months(m) == "July" 
            figure;
            subplot(3,1,1)
            plot(1:7*24,q_ev(1:7*24));
            hold on
            plot(1:7*24,q_d(1:7*24));
            plot(1:7*24, gamma_ra(1:7*24)*q_ra);
            hold off
            ylabel("Energy (kWh)")
            legend("EV Charging", "Arbitrage", "Resource Adequacy")

            subplot(3,1,2)
            a = YearLMP_vec(f_start:f_end);
            plot(1:7*24,a(1:7*24))
            ylabel("LMP Price ($/kWh)")

            subplot(3,1,3)
            plot(1:7*24,S_t(1:7*24))
            hold on
            plot(1:7*24,q_r(1:7*24))
            plot(1:7*24, -1.*(q_ev(1:7*24)+q_d(1:7*24)+gamma_ra(1:7*24).*q_ra))
            xlabel("Time (in hours)")
            ylabel("Energy (in kWh)")
            legend("SOC","Recharge Curve","Discharge Curves")
            saveas(gca,experiment_folder + "/July_stats.png")
        end
        full_powerOpt = zeros(vecLength,6);
        full_powerOpt(1:vecLength,1) = q_d;
        full_powerOpt(1:vecLength,2) = q_r;
        full_powerOpt(1:vecLength,3) = q_ev;
        full_powerOpt(1:vecLength,4) = S_t(2:end);
        full_powerOpt(1:vecLength,5) = YearLMP_vec(f_start:f_end);
        area_graph(f_start:f_end,1) = -1* transpose(YearLMP_vec(f_start:f_end)) .* q_d + transpose(YearLMP_vec(f_start:f_end)) .* q_r;
%         area_graph(f_start:f_end,2) = transpose(YearLMP_vec(f_start:f_end)) .* q_r ;
        area_graph(f_start:f_end,2) = -1*q_ev .* EV_unitPrice - 0.001;
        area_graph(f_start:f_end,3) =  -1 * P_RA(m) .* (gamma_ra .* q_ra);
        area_graph(f_start:f_end,4) = -1 * transpose(YearLMP_vec(f_start:f_end)) .* (gamma_ra .*q_ra);
%         area_graph(f_start:f_end,5) = -1* transpose(YearLMP_vec(f_start:f_end)) .* q_d + transpose(YearLMP_vec(f_start:f_end)) .* q_r  + q_ev .* EV_unitPrice + -1 * P_RA(m) .* (gamma_ra .* q_ra) +  -1 * transpose(YearLMP_vec(f_start:f_end)) .* (gamma_ra .*q_ra);
        area_graph_monthly(m,1) = sum(-1* transpose(YearLMP_vec(f_start:f_end)) .* q_d) - sum(transpose(YearLMP_vec(f_start:f_end)) .* q_r );
%         area_graph_monthly(m,2) = sum(transpose(YearLMP_vec(f_start:f_end)) .* q_r );
        area_graph_monthly(m,2) = sum(-q_ev .* EV_unitPrice )- 0.001;
        area_graph_monthly(m,3) = sum( -1 * P_RA(m) .* q_ra) ;
        area_graph_monthly(m,4) = sum(-1 * transpose(YearLMP_vec(f_start:f_end)) .* (gamma_ra .*q_ra));
        
%         area_graph_monthly(m,5) = sum(-1* transpose(YearLMP_vec(f_start:f_end)) .* q_d + transpose(YearLMP_vec(f_start:f_end)) .* q_r  + q_ev .* EV_unitPrice + -1 * P_RA(m) .* (gamma_ra .* q_ra) +  -1 * transpose(YearLMP_vec(f_start:f_end)) .* (gamma_ra .*q_ra));
%         area(area_graph)
        full_powerOpt(1,6) = sum(EV_unitPrice .* q_ev); 
        full_powerOpt(2,6) = YearLMP_vec(f_start:f_end) * q_d;
        full_powerOpt(3,6) = YearLMP_vec(f_start:f_end) * q_r; 
        full_powerOpt(4,6) = P_RA(m) * q_ra;
        full_powerOpt(5,6) = YearLMP_vec(f_start:f_end)*(gamma_ra .*q_ra);
        writetable(array2table(full_powerOpt),filename,'Sheet',months(m),'Range','A1');
    end
    area_graph_yearly(1) = sum(area_graph_monthly(1:end,1));
    area_graph_yearly(2) = sum(area_graph_monthly(1:end,2));
    area_graph_yearly(3) = sum(area_graph_monthly(1:end,3));
    area_graph_yearly(4) = sum(area_graph_monthly(1:end,4));
    
    figure;
    area(-1*area_graph_monthly);
    xlabel('Month of Year');
    ylabel('Revenue');
    xlim([1,12]);
    legend({'Arbitrage', 'EV Charging', 'RA Revenue (from Utility)', 'RA Revenue (from discharge)'});
    title(sprintf('Monthly Revenue Streams for Interconnection Capacity %i kW', q_r_max));
    %     figure;
%     area(-1*area_graph);
    figure;
    pie(-1*area_graph_yearly)
    legend({'Arbitrage', 'EV Charging', 'RA Revenue (from Utility)', 'RA Revenue (from discharge)'});
    title(sprintf('Pie Chart of Revenue Streams Percentages for Interconnection Capacity %i kW', q_r_max));
    save(experiment_folder + "/ra_month_total_revenue.mat", 'ra_month', 'monthlyMaxRevenue')
    total_profit = sum(-1 .* monthlyMaxRevenue)
    figure, plot(-monthlyMaxRevenue)
    set(gca,'xtick',1:12,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
    legend("Cumulative Profit", "Monthly Profit");
    xlabel("Month of Year")
    ylabel("Profit (in $)")
    title("Profit for a Single Unit per Month")
    saveas(gca,experiment_folder+"/profits.jpg")
    figure, plot(ra_month)
    set(gca,'xtick',1:12,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
    xlabel("Month of Year")
    ylabel("Power Capacity (kW)")
    title("Resource Adequacy Committed Capacity vs Month of Year")
    saveas(gca,experiment_folder+"/ra_month.jpg")
    figure;
    subplot(3,1,1)
    plot(1:7*24,q_ev(1:7*24));
    hold on
    plot(1:7*24,q_d(1:7*24));
    plot(1:7*24, gamma_ra(1:7*24)*q_ra);
    hold off
    ylabel("Energy (kWh)")
    legend("EV Charging", "Arbitrage", "Resource Adequacy")

    subplot(3,1,2)
    a = YearLMP_vec(f_start:f_end);
    plot(1:7*24,a(1:7*24))
    ylabel("LMP Price ($/kWh)")

    subplot(3,1,3)
    plot(1:7*24,S_t(1:7*24))
    hold on
    plot(1:7*24,q_r(1:7*24))
    plot(1:7*24, -1.*(q_ev(1:7*24)+q_d(1:7*24)+gamma_ra(1:7*24).*q_ra))
    xlabel("Time (in hours)")
    ylabel("Energy (in kWh)")
    legend("SOC","Recharge Curve","Discharge Curves")
    saveas(gca,experiment_folder+"/December_stats.png")
    toc
end