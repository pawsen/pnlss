function hfig = fPlotSubSpace(data,models,G,covG,freq,fs,na,max_r)
% freq: excited frequencies

hfig = {}; % handle for storing plots

KSSs = data.KSSs;
KLMs = data.KLMs;
stabilizedSS = data.stabilizedSS;
unstablesLM = data.unstablesLM;

[m,p,F] = size(G);
min_na = min(na); % Minimal model order that is scanned
max_na = max(na); % Maximal model order that is scanned

figure
semilogy(KSSs','.') % Cost function of all subspace models
hold on;
semilogy(stabilizedSS','o','Color',[0.5 0.5 0.5]) % Encircle stabilized models in gray
% if any KLMs are not NaN, then we did optimize
optimize = sum(any(~isnan(KLMs)));
if optimize > 0
    try
        set(gca,'ColorOrderIndex',1) % Restart from the same color as in the KSSs plot
    catch % Older Matlab versions restart by default from the same color, and 'ColorOrderIndex' is not available
    end
    semilogy(KLMs','*') % Cost function of all LM-optimized models
    try
        set(gca,'ColorOrderIndex',1) % Restart from the same color as in the KSSs plot
    catch % Older Matlab versions restart by default from the same color, and 'ColorOrderIndex' is not available
    end
    semilogy(unstablesLM','ko') % Encircle unstable models
end
ylabel('V_{WLS}')
xlabel('r')
legend(cellstr([repmat('n = ',max_na - min_na + 1,1) num2str((min_na:max_na)')]));
title({'Cost functions of subspace models (dots, stablilized models encircled in gray)', ...
    'and of LM-optimized models (stars, unstable models encircled in black)'})
xlim([min_na (max_r + 1)])
set(gcf,'Name','Summary')
set(gcf,'NumberTitle','off')
hfig(end+1) = {{gcf, "sub_cost"}};

for n = na % Plot FRFs of the best model for each model order
    figure;
    model = models{n};
    try % Maybe no stable models were estimated
        A = model{1};
        B = model{2};
        C = model{3};
        D = model{4};
        GModel = fss2frf(A,B,C,D,freq/fs);
        fPlotFrfMIMO(GModel,freq,'b-','LineWidth',4);
        fPlotFrfMIMO(G,freq,'r.');
        fPlotFrfMIMO(G-GModel,freq,'.','Color',[0.5,0.5,0.5]); % grey color
        temp = zeros(m*p,F);
        for i = 1:m*p
            temp(i,:) = covG(i,i,:);
        end
        temp2 = zeros(size(G));
        for f = 1:F
            temp2(:,:,f) = reshape(temp(:,f),p,m);
        end
        stdG = sqrt(temp2);
        fPlotFrfMIMO(stdG,freq,'k.');
        set(gcf,'Name',['n = ' num2str(n)])
        set(gcf,'NumberTitle','off')
        
        legend('BLA (parametric)','BLA (nonpar)','residual','standard deviation')
    catch
        gca;
        title('No stable models');
    end
    hfig(end+1) = {{gcf, sprintf('subspace_n%d',n)}};
end

end