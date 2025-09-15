% RecordNew = rmoutliers(C);
% 
% nbins = 30;
% figure
% h = histogram(RecordNew(:,1),nbins,'Normalization','pdf','FaceColor','r','FaceAlpha', 0.3);
% legend(['$\hat{\pi}_+$'],'Interpreter','latex','FontSize',30)
% xlim([0.875, 0.925])
% figure
% h = histogram(RecordNew(:,2),nbins,'Normalization','pdf','FaceColor','b','FaceAlpha', 0.3);
% legend(['$\hat{\pi}_-$'],'Interpreter','latex','FontSize',30)
% xlim([0.375, 0.425])
% figure
% h = histogram(RecordNew(:,3),nbins,'Normalization','pdf','FaceColor','y','FaceAlpha', 0.8);
% legend(['$\hat{\alpha}$'],'Interpreter','latex','FontSize',30)
% xlim([0.27, 0.33])
% mean(RecordNew)
% std(RecordNew)

% %%%%%%%
RecordNew1 = rmoutliers(Record1);
RecordNew2 = rmoutliers(Record2);
nbins = 30;
figure
h = histogram(RecordNew2(:,3),nbins,'Normalization','pdf','FaceColor','b','FaceAlpha', 0.3);
hold on
h = histogram(RecordNew1(:,3),nbins,'Normalization','pdf','FaceColor','r','FaceAlpha', 0.3);
legend('$M_{3,K}$ when $\alpha = 0.3$', '$M_{3,K}$ when $\alpha = 0.6$','Interpreter','latex','FontSize',26)