clear variables
close all

addpath ('./util')
%% Read evolved parameters
rootDir = '../evolved_parameters/';
BVs = [16, 36, 64, 100, 144];
displayPlot = 1;
   
paramsAll = cell(length(BVs),1);
nParams = zeros(length(BVs),1);
nTrials = zeros(length(BVs),1);
for b = 1:length(BVs)
    bv = BVs(b);
    name = [num2str(bv) 'BV'];
    params = csvread([rootDir name '/params.csv']);
    paramsAll{b} = params;
    nParams(b) = size(params,2);
    nTrials(b) = size(params,1);
end
%% Plot STDP curves
modes = ["EE" "EI" "IE"];

% y-axis
W_m = cell(length(BVs), length(modes));
W_p = cell(length(BVs), length(modes));
% x-axis
t_neg = 0:1:100;
t_pos = 0:1:100;

% variables to store AUC values
AUC_m = cell(length(BVs), length(modes));
AUC_p = cell(length(BVs), length(modes));

for m = 1:length(modes)
    mode = modes(m);
    
    switch mode
        case 'EE'
            cs = [1 4 7 8];
        case 'EI'
            cs = [2 5 9 10];
        case 'IE'
            cs = [6 3 11 12];
    end

    for b = 1:length(BVs)
        W_p{b,m} = zeros(nTrials(b),length(t_pos)); 
        W_m{b,m} = zeros(nTrials(b),length(t_neg));
        
        AUC_m{b,m} = zeros(nTrials(b),1); 
        AUC_p{b,m} = zeros(nTrials(b),1);
        
        for i = 1:nTrials(b)
            p = paramsAll{b}(i,cs);
            A_p = p(1);
            A_m = p(2);
            tao_p = p(3);
            tao_m = p (4);
            
            % compute weight decay from STDP parameters
            w_p = A_p*exp(-t_pos/tao_p);
            w_m = A_m*exp((-t_neg)/tao_m);

            W_m{b,m}(i,:) = w_m; 
            W_p{b,m}(i,:) = w_p;
            
            AUC_p{b,m}(i,:) = abs(trapz(t_pos,w_p));
            AUC_m{b,m}(i,:) = abs(trapz(t_neg,w_m));
        end
    
        if displayPlot
            figure(1)
            subplot(length(BVs), length(modes), 3*(b-1)+m)
            indi = 1:nTrials(b);

            [aline_p, aFill_p] = stdshade([W_p{b,m}(indi,:)],...
                0.2,c1,[t_pos]);
            hold on
            [aline_m, aFill_m] = stdshade([flip(W_m{b,m}(indi,:),2)],...
                0.2,c2,[flip(-t_neg)]);
            
            xlimit = get(gca, 'xlim');
            

            ymax = 4*10^(-3);
            ylimit = [-ymax ymax];
            ylim(ylimit)
            
            plot(xlimit, [0 0], 'LineWidth',1, 'Color', [105,105,105]/255)
            plot([0 0], ylimit, 'LineWidth',1,'Color', [105,105,105]/255)
            
            hold off

            if b == 1
                if mode == "EE"
                    title('MT-MSTd (E-STDP)')
                elseif mode == "EI"
                    title('MSTd-Inh (E-STDP)')
                elseif mode == "IE"
                    title('Inh-MSTd (I-STDP)')
                end 
            end
            if b == length(BVs)
                xlabel('t_{post} - t_{pre} (ms)')
            end
            if m == 1
                ylabel('\Delta w')
            end
        end
    end
end
%% Calculate mean and std of AUCs
meanAUC_m = cellfun(@mean,AUC_m);
stdAUC_m = cellfun(@std,AUC_m);
meanAUC_p = cellfun(@mean,AUC_p);
stdAUC_p = cellfun(@std,AUC_p);
%% Visualize AUCs
meanAUC_m_f = flipud(meanAUC_m);
meanAUC_p_f = flipud(meanAUC_p);
stdAUC_m_f = flipud(stdAUC_m);
stdAUC_p_f = flipud(stdAUC_p);
BVs_f = fliplr(BVs);

h2 = figure(3);
clf
for m = 1:length(modes)
    subplot(1,3,m)
    mode = modes(m);
    y = [meanAUC_m_f(:,m), meanAUC_p_f(:,m)];
    err = [stdAUC_m_f(:,m), stdAUC_p_f(:,m)];
    
    if mode == "IE"
        y = flip(y,2);
        err = flip(err,2);
    end
    
    ngroups = size(y, 1);
    nbars = size(y, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    
    hold on
    w = 0.4;
    x1 = (1:5)-w/2;
    x2 = (1:5)+w/2;
    b1 = barh(x1,y(:,1), 'FaceColor','flat', 'BarWidth', w, 'FaceColor',c2);
    b2 = barh(x2,y(:,2), 'FaceColor','flat', 'BarWidth', w, 'FaceColor',c1);
    errorbar(y(:,1), x1, err(:,1), 'horizontal','k', 'linestyle', 'none');
    errorbar(y(:,2), x2, err(:,2), 'horizontal','k', 'linestyle', 'none');
    
    xlim([0 0.25])
    ylim([0.5 5.5])
    if m == 1
        ylabel('Number of MSTd Neurons')
    end
    xlabel('Area over LTD/ LTP')
    
    set(gca,'YTick', 1:5);
    set(gca,'YTickLabel', BVs_f);
    
    if m == 1
        f=get(gca,'Children');
        legend([f(3),f(4)],'LTP','LTD')
    end
        
    if mode == "EE"
        title('MT-MSTd (E-STDP)')
    elseif mode == "EI"
        title('MSTd-Inh (E-STDP)')
    elseif mode == "IE"
        title('Inh-MSTd (I-STDP)')
    end 
    hold off
end
h2.Position = [1 1 900 450];