% Aqui esta para graficar lo que auierop
% to plot the weight
figure%,stem(patterns.all.aversive(:,1),'filled'),hold on
stem(patterns.all.aversive(:,6),'filled'),ylim([-0.8 0.8]),hold on
ylabel('Weight'),xlim([0 76])
yline(1/sqrt(75))
xline(numberD)
cond = abs(patterns.all.aversive(:,6)) > (1/sqrt(75));
cond = flip(cond);

%                     period = [aversiveTS_run(1)./1000+200  aversiveTS_run(1)./1000+205];


% for REM
% [i ii] = max(REM.aversive(:,2) - REM.aversive(:,1));
% i = ((REM.aversive(ii,2) - REM.aversive(ii,1))/2) + REM.aversive(ii,1);
% period = [i-2.5 , i+2.5];
period = [REM.all(5,1)+23 , REM.all(5,1)+28];
% %for NREM
% period = [ripplesD(5000,2)-2.5 ripplesD(5000,2)+2.5];
% period = [6726 6731]

L1 = Restrict(dHPC,period);
L2 = Restrict(vHPC1,period);

figure,
subplot(211),plot(L1(:,1),L1(:,2)),hold on
L1 = Restrict(ripplesD(:,2),period);
for i = 1 :length(L1)
    xline(L1(i))
end
subplot(212),plot(L2(:,1),L2(:,2))
L2 = Restrict(ripplesV(:,2),period);
for i = 1 :length(L2)
    xline(L2(i))
end
xlim(period)

group = [group_vHPC;group_dHPC];
s = [spks_dHPC ; spks_vHPC];
C = size(group_dHPC,1);
figure ()
ii = 0;

for i = 1 : length(group)
    cluster = group(i,1);
%     celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%     if celltype
        ii = ii+1;
        spks = s(s(:,1)==cluster,2);
        spks = Restrict(spks,period);
        spks = spks';
        xspikes = repmat(spks,3,1);
        yspikes = nan(size(xspikes));
        
        if not(isempty(yspikes))
            yspikes(1,:) = ii-1;
            yspikes(2,:) = ii;
        end
        
        %                             if i >= C
        if i>length(group_vHPC)
            plot(xspikes,yspikes,'Color','r','LineWidth',1),hold on
        else
            plot(xspikes,yspikes,'Color','k','LineWidth',1),hold on
        end
        %                             yline(ii)
%     end
end
xlim(period)
% ylim([0 75])
% ylabel('Neuron#')
% figure

in = InIntervals(bins,period);
plot(bins(in),(a(6,in)))%,ylim([-15 30])
ylabel('Assembly activity')


a = assembly_activity(patterns.all.aversive, Spikes');
