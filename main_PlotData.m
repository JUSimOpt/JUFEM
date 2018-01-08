%clear; clc; close all


% nuu = [0.45,0.475,linspace(0.49,0.499,10)];
% intRules = {'Full','Reduced','IsoStab','IsoVolStab'};
% 
% nx = length(nuu);
% displacements = zeros(nx,4);
% 
% for ix = 1:nx
%     nu = nuu(ix);
%     tipDisps = [];
%     
%     for i = 1:4
%         filename = ['results//HansboBeam_nu_',num2str(nu),'_',intRules{i},'.mat'];
%         load(filename)
%         u = OUT.step.inc(end).nodhistvar.u;
%         U = [u(1:3:end),u(2:3:end),u(3:3:end)];
%         P = model.mesh.P;
%         ind = find(P(:,1)==max(P(:,1)) & P(:,3)==min(P(:,3)));
%         tipDisp = min(U(ind,3));
%         displacements(ix,i) = tipDisp;
%     end
%     1;
% end
% 
% %%
% xfigure; hold on;view(2)
% msize = 6;lwidth = 1;
% plot(nuu,displacements(:,1),'-o','MarkerSize',msize,'LineWidth',lwidth)
% plot(nuu,displacements(:,2),'-*','MarkerSize',msize,'LineWidth',lwidth)
% plot(nuu,displacements(:,3),'-s','MarkerSize',msize,'LineWidth',lwidth)
% plot(nuu,displacements(:,4),'-+','MarkerSize',msize,'LineWidth',lwidth)
% legend('Full integration',...
%        'One point integration of the volumetric term',...
%        'One point integration with stabilized isochoric term',...
%         'One point integration with stabilization')
% set(gca,'FontName','Times New Roman','FontSize',14)
% ylabel('Tip displacement');
% xlabel('\nu');
% grid on

%% h
% intRules = {'Full','Reduced','IsoStab','IsoVolStab'};
% hh = zeros(4,1);
% nx = 4;
% displacements = zeros(nx,4);
% for ix = 1:nx
%     tipDisps = [];
%     for i = 1:4
%         filename = ['results//HansboBeam_R',num2str(i-1),'_',intRules{ix},'.mat'];
%         load(filename)
%         u = OUT.step.inc(end).nodhistvar.u;
%         U = [u(1:3:end),u(2:3:end),u(3:3:end)];
%         P = model.mesh.P;
%         nodes = model.mesh.nodes;
%         XC = P(nodes(1,:),:);
%         dx = norm(XC(2,:)-XC(1,:));
%         dy = norm(XC(4,:)-XC(1,:));
%         dz = norm(XC(5,:)-XC(1,:));
%         h = (dx*dy*dz)^(1/3);
%         hh(i) = h;
%         ind = find(equalToPrecision(P(:,1),max(P(:,1))) & equalToPrecision(P(:,3),min(P(:,3))) );
%         tipDisp = min(U(ind,3));
%         displacements(ix,i) = tipDisp;
%     end
%     1;
% end
% 
% %
% xfigure; hold on;view(2)
% msize = 6;lwidth = 1;
% plot(hh,displacements(1,:),'-o','MarkerSize',msize,'LineWidth',lwidth)
% plot(hh,displacements(2,:),'-*','MarkerSize',msize,'LineWidth',lwidth)
% plot(hh,displacements(3,:),'-s','MarkerSize',msize,'LineWidth',lwidth)
% plot(hh,displacements(4,:),'-+','MarkerSize',msize,'LineWidth',lwidth)
% legend('Full integration',...
%        'One point integration of the volumetric term',...
%        'One point integration with stabilized isochoric term',...
%         'One point integration with stabilization')
% set(gca,'FontName','Times New Roman','FontSize',14)
% ylabel('Tip displacement');
% xlabel('$h$','Interpreter','latex')
% grid on

%% Loads
% intRules = {'Full','Reduced','IsoStab','IsoVolStab'};
% displacements = zeros(4,10);
% f = linspace(1,10,10); %10; %GN/mm^3
% for i = 1:4
%     filename = ['results//HansboBeam_F_',intRules{i},'.mat'];
%     load(filename)
%     for j = 1:10
%         u = OUT.step.inc(j).nodhistvar.u;
%         U = [u(1:3:end),u(2:3:end),u(3:3:end)];
%         P = model.mesh.P;
%         ind = find(equalToPrecision(P(:,1),max(P(:,1))) & equalToPrecision(P(:,3),min(P(:,3))) );
%         tipDisp = min(U(ind,3));
%         displacements(i,j) = tipDisp;
%     end
%     
% end
% xfigure; hold on;view(2)
% msize = 6;lwidth = 1;
% plot(displacements(1,:),f,'-o','MarkerSize',msize,'LineWidth',lwidth)
% plot(displacements(2,:),f,'-*','MarkerSize',msize,'LineWidth',lwidth)
% plot(displacements(3,:),f,'-s','MarkerSize',msize,'LineWidth',lwidth)
% plot(displacements(4,:),f,'-+','MarkerSize',msize,'LineWidth',lwidth)
% legend('Full integration',...
%        'One point integration of the volumetric term',...
%        'One point integration with stabilized isochoric term',...
%         'One point integration with stabilization')
% set(gca,'FontName','Times New Roman','FontSize',14)
% ylabel('Volume force ($GN/mm^3$)','Interpreter','latex');
% xlabel('Tip displacement','Interpreter','latex')
% grid on
