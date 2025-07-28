[sortedPValues, Idx] = sort(p_values_most_conservative_comps_scaled(1:200));
[sortedMostSig,sortedMostSigIdx] = sort(min(reshape(p_values_most_conservative_comps_vec,5,[]),[],1));
sortedMostSig(p_values_most_conservative_comps_vec_sorted_ind) = p_values_most_conservative_comps_vec_sorted;
[sortedMostSigAbridged, sortedMostSigAbridgedIdx] = sort(min(reshape(sortedMostSig,5,[]),[],1));
figure
% plot(1:200, sortedPValues)
hold on
plot(1:3,sortedMostSigAbridged(1:3),'ko','MarkerFaceColor','k');
plot(4:40,sortedMostSigAbridged(4:40),'ko');
xlim([0,41])
ylim([10^(-13),10^2])
yticks(10.^[-10,-5,0])
for i = 1:40
    AbridgedSortedNames{i} = erase(AllCrossStatsNames{sortedMostSigAbridgedIdx(i),2}(1:(strfind(AllCrossStatsNames{sortedMostSigAbridgedIdx(i),2}, '>')-2)),'split ');
end
set(gca,'XTick',1:40,'XTickLabel',AbridgedSortedNames,'YScale','Log')
ylabel('Holm-Bonferroni Corrected p Value')
hold on
ylim([10^(-14), 1])
% yticklabels([10.^([-3:1]*5)])
plot([0:40],[0.05./(205-5*[0:40])],'k--')
% yline(0.05,'k--')
set(gca, 'color', 'none');
% set(gca,'visible','off')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
set(gca,'box','off')