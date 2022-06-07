clear; clc;
figure(1); clf(1); 
figure(2); clf(2);
figure(3); clf(3);

colors(1:6,1:3) = [[1 0.47 0]; [0.4 0 1]; [0.07 1 0]; [0.9 0 1]; [0.96 1 0]; [1 0 0.2]];
fill_colors = colors.*0.8;

% this code will run three different scenarios
% scenario 1, three levels of density dependence
% scenario 2, vary the steepness of the trade-off function
% scenario 3, same as scenario 1 but start the trait above the ESS

% this code requires the following additional files to run
% GEMv2_logistic_growthalt_4 - the core GEM function
% LGalt_model - the bd_logistic function file
% QG_model_LG_alt - the bd_logisitc quantitative genetics function file
% jbfill - an area plotting function
% pick_individuals - individual selection function called from GEMv2_logistic_growthalt_3
% medians_and_cis - gets medians and cis called from GEMv2_logistic_growthalt_3
% NEEA_through_time - calculates the running NEEA from LRS

% pick your scenario here
scenario = 1;

% find bslope to achieve target Ks
% given also bmax = 1.8, dmin = 0.3, and to_slope = dmin/bmax^2
    target_K = [10 20 40];
    bslope = (1.8 - 0.3)./(2.*target_K); % density dependence

% parameters and replicates
    dslope = bslope; % density dependence of deaths = that of births
    to_slope = 0.3/1.8^2; % specify slope of trade-off
    cv = [0.3 0.3 0.3]; % vector of cv's, one for each version
    h_2 = [0.75 0.75 0.75]; % vector of h2's, one for each version
    
    num_replicates = 5; % number of GEM simulations
    t_max = 300; % time span to run simulations 
    tspan = [0 t_max]; % start and end times
    
    if scenario == 1
        b = 1.8; % specify initial birth trait
        d = to_slope*b^2; % initial death trait
        y0 = [5 5 5]; % set starting abundances
        topyaxis = [80 80 80]; % set the top of the y axis
        Rcull = [0 0 0]; % set culling level
        cr = 0; % set culling rate
        titles = {'K_{init}=10','K_{init}=20','K_{init}=40'};
    elseif scenario == 2
        bslope = [0.06, 0.06, 0.06]; % density dependence
        dslope = bslope; % density dependence
        b = 1.8;
        to_slopes = [0.06 0.09 0.12]; % specify slope of trade-off
        ds = to_slopes.^2;
        y0 = [5 5 5]; % set starting abundances
        Rcull = [0 0 0]; % set culling level
        cr = 0; % set culling rate
        topyaxis = [50 50 50]; % set the top of the y axis
        titles = {'s=0.06','s=0.09','s=0.12'};
    elseif scenario == 3
        b_vector = [1.8 6]; % specify initial birth trait
        cv = [0.3 (0.3*1.8)/6]; % vector of cv's, one for each version
        bslope = [bslope(1) bslope(1)];
        dslope = bslope; % density dependence
        d_vector = to_slope.*b_vector.^2; % initial death trait
        y0 = [5 5]; % set starting abundances
        topyaxis = [30 30]; % set the top of the y axis
        Rcull = [0 0]; % set culling level
        cr = 0; % set culling rate
        titles = {'b_{max,init}=1.8','b_{max,init}=6'};
    end
    
for j = 1:length(bslope)
    
    if scenario == 2 
        to_slope = to_slopes(j);
        d = ds(j);
    elseif scenario == 3
        b = b_vector(j);
        d = d_vector(j);
    end

    K = (b - d)/(dslope(j) + bslope(j)); % calculate carrying capacity
    b_ESS = 1/(2*to_slope); % calculate ESS birth max trait
    d_ESS = to_slope*b_ESS^2; % calculate ESS death min trait
    K_ESS = (b_ESS - d_ESS)/(dslope(j) + bslope(j)); % calculate K at ESS
    
    % actually call the GEM function
    [x_dist,stand_times, R_data_out, x_data_out, x_var_data_out] = GEMv2_logistic_growthalt_4(j, b, d, to_slope, bslope(j), dslope(j), cv(j), h_2(j), num_replicates, y0(j), t_max, Rcull(j), cr);

    realized_initial_b = x_data_out(1,1); % get the average of the realized starting epsilons 
    realized_initial_b_var = x_var_data_out(1,1);
    realized_initial_d = to_slope*realized_initial_b^2;
    realized_initial_K = (realized_initial_b - realized_initial_d)/(bslope(j) + dslope(j));

    ode = @(t,y) LGalt_model(t,y,realized_initial_b,realized_initial_d,bslope(j),dslope(j),cr,Rcull(j)); % compile function and call
        [t1,y1] = ode45(ode, tspan, y0(j)); % return time and population density vectors

    ode = @(t,y) QG_model_LG_alt(t,y,realized_initial_d,bslope(j),dslope(j),h_2(j),realized_initial_b_var,to_slope,cr,Rcull(j));
        y0_QG = [y0(j) realized_initial_b];
        [t2,y2] = ode45(ode, tspan, y0_QG); % return time and population density vectors
    
    % call function to calculate NEEA through time
    LRSmax = NEEA_through_time(b,d,to_slope,stand_times,R_data_out)

        
        
    figure(1);
    hold on;
    subplot(5,3,j);
        box on;
        jbfill(stand_times,R_data_out(2,:),R_data_out(3,:),fill_colors(2,:),'w',1,0.2); hold on;
        h1 = plot(stand_times,R_data_out(1,1:length(stand_times)),'-','Color',colors(2,:),'LineWidth',2);
        if scenario == 1 || 4
            h2 = plot(t2,y2(:,1),'-','Color',colors(1,:)','LineWidth',2);
        end
        %plot([0 t_max],[realized_initial_K realized_initial_K],'--k');
        %plot([0 t_max],[K_ESS K_ESS],'--','Color','k');
        ylim([0 topyaxis(j)]);
        xlim([0 t_max]);
        title(titles(j));

    subplot(5,3,j+3);
        box on;
        pop_sizes = R_data_out(1,:);
        TEA_at_N = pop_sizes.*bslope(j) + sqrt((bslope(j).*pop_sizes).^2 + pop_sizes.*dslope(j)./(to_slope));        

        jbfill(stand_times,x_data_out(2,:),x_data_out(3,:),fill_colors(2,:),'w',1,0.2); hold on;
        plot(stand_times,x_data_out(1,:),'-','Color',colors(2,:),'LineWidth',2);
        if scenario == 1 || 4
            plot(t2,y2(:,size(y2,2)),'-','Color',colors(1,:),'LineWidth',2);
        end
        h4 = plot(stand_times,TEA_at_N,'--','Color',colors(4,:),'LineWidth',1);
        h5 = plot([0 t_max],[b_ESS b_ESS],'--','Color',colors(1,:));
        ylim([0 7]);
        xlim([0 t_max]);

        if scenario == 1 || 2
            tea_pop = mean(R_data_out(1,50:end));
            TEA = tea_pop*bslope(j) + sqrt((bslope(j)*tea_pop)^2 + tea_pop*dslope(j)/(to_slope));        
        elseif scenario == 2 || 3
            TEA = cull(j)*bslope(j) + sqrt((bslope(j)*cull(j))^2 + cull(j)*dslope(j)/(to_slope));  
        end
        %h5 = plot([0 t_max],[TEA TEA],'--','Color',colors(4,:));
        
    subplot(5,3,j+6); 
        box on;
        jbfill(stand_times,x_var_data_out(2,:),x_var_data_out(3,:),fill_colors(2,:),'w',1,0.2); hold on;
        plot(stand_times,x_var_data_out(1,:),'-','Color',colors(2,:),'LineWidth',2);
        plot([0 t_max],[realized_initial_b_var realized_initial_b_var],'--k');
        if scenario < 3
            ylim([0 0.8]);
        else
            ylim([0 3]);
        end
        xlim([0 t_max]);

    subplot(5,3,j+9);
        box on; hold on;
        % plot points
        plot(x_dist(find(x_dist(:,2)==0 & x_dist(:,4)>350),1),x_dist(find(x_dist(:,2)==0 & x_dist(:,4)>350),3),'.k');
        %plot(x_dist(find(x_dist(:,2)==0 & x_dist(:,4)<1),1),x_dist(find(x_dist(:,2)==0 & x_dist(:,4)<1),3),'.','Color',[0.5 0.5 0.5]);
        % plot vertical bars for attractors
        plot([b_ESS b_ESS],[0 20],'-','Color',colors(1,:),'LineWidth',2);
        plot([TEA TEA],[0 20],'--','Color',colors(4,:),'LineWidth',2);

        ylim([0 20]);
        xlim([0 10]);

        realized_b_dist = x_dist(find(x_dist(:,2)==0),1);
        realized_d_dist = to_slope.*realized_b_dist.^2;
        realized_lrs_dist = x_dist(find(x_dist(:,2)==0),3);
        time_of_birth_dist = x_dist(find(x_dist(:,2)==0),4);
        n_at_birth = x_dist(find(x_dist(:,2)==0),5);
        expected_lrs = max((realized_b_dist - bslope(j).*n_at_birth),0)./(realized_d_dist + dslope(j).*n_at_birth);
        lrs_diff = realized_lrs_dist - expected_lrs;
        expected_at_ess = (b_ESS - bslope(j).*K_ESS)./(d_ESS + dslope(j).*K_ESS);

    subplot(5,3,j+12);
        box on; hold on;
        xx = 1:100;
        TEA_at_N = xx.*bslope(j) + sqrt((bslope(j).*xx).^2 + xx.*dslope(j)./(to_slope));        
        plot(xx,TEA_at_N,'-k');
        plot(R_data_out(1,1:length(stand_times)),x_data_out(1,:),'-','Color',colors(2,:),'LineWidth',2);
        plot(K_ESS,b_ESS,'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:));

    if scenario == 3
        pop_sizes = R_data_out(1,:);
        TEA_at_N = pop_sizes.*bslope(j) + sqrt((bslope(j).*pop_sizes).^2 + pop_sizes.*dslope(j)./(to_slope));        

        figure(2);
        subplot(2,2,j);
            box on; hold on;
            fit_gradient_QG = 1-2.*to_slope.*x_data_out(1,:);
            fit_gradient_LRS = 1./(d + dslope(j).*R_data_out(1,1:length(stand_times)));
            z1 = plot(stand_times,fit_gradient_QG,':','Color',colors(1,:),'LineWidth',2);
            %plot(stand_times,fit_gradient_LRS,'-','Color',colors(2,:),'LineWidth',2);
            plot([0 t_max],[0 0],'--','Color','k');
            title(titles(j));
        subplot(2,2,j+2);
            box on;
            jbfill(stand_times,x_data_out(2,:),x_data_out(3,:),fill_colors(2,:),'w',1,0.2); hold on;
            z2 = plot(stand_times,x_data_out(1,:),'-','Color',colors(2,:),'LineWidth',2);
            z3 = plot(stand_times,TEA_at_N,'--','Color',colors(4,:),'LineWidth',1);
            z4 = plot(t2,y2(:,size(y2,2)),'-','Color',colors(1,:),'LineWidth',2);
            z5 = plot([0 t_max],[b_ESS b_ESS],'--','Color',colors(1,:));
            ylim([0 7]);
            xlim([0 t_max]);
    end
end

figure(1);subplot(5,3,8);
    xlabel('Time','FontSize',12);
figure(1);subplot(5,3,1);
    ylabel('Abundance','FontSize',12);
figure(1);subplot(5,3,4);
    ylabel('Trait mean','FontSize',12);
figure(1);subplot(5,3,7);
    ylabel('Trait variance','FontSize',12);
figure(1);subplot(5,3,10);
    ylabel({['Lifetime reproductive'],['success']},'FontSize',12);
figure(1);subplot(5,3,13);
    ylabel('b_{max}','FontSize',12);
figure(1);subplot(5,3,11);
    xlabel('Birth trait','FontSize',12);
figure(1);subplot(5,3,14);
    xlabel('Abundance','FontSize',12);
figure(1);subplot(5,3,13);
    ylabel('Birth trait','FontSize',12);
figure(1);subplot(5,3,3);
if scenario == 1
    legend([h1 h2 h4 h5],'GEM','QG','ESS','TEA','Location','Best');
elseif scenario == 2
    legend([h1 h4 h5],'GEM','ESS','TEA','Location','Best');
end

figure(2);
subplot(2,2,3);
    xlabel('Time','FontSize',12);
    ylabel('b_{max}','FontSize',12);
subplot(2,2,4);
    xlabel('Time','FontSize',12);
    legend([z1 z2 z3 z4 z5],'QG fitness gradient','GEM','TEA','QG','ESS','Location','Best');
subplot(2,2,1);
    ylabel('Fitness gradient','FontSize',12);

%     gtext('K_{ESS}');
%     gtext('K_{init}');
%     gtext('b_{init}');
%     gtext('ESS');
%     gtext('TEA');

shg;

gtext('TEAs');