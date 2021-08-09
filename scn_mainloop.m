%% Interaction among SCN, HPA, HPT, HNS
% Model parameters and plots



%% First get datafiles for all figures


clear all;

% these are the studied conditions
conditions_NPY = {'NPYlo','NPYhi'};		% low and high values of ARC -> TRH efficacy
conditions_DA = {'DAlo','DAmed'};		% fast and slow HNS oscillation (i.e. low and intermediate values of striatal dopamine)
conditions_ZG = {'SP','LP','SH'};		% zeitgeber period and photoperiod combinations
values_NPY = [0.1 25];
values_DA = [28.5714 114.2857];
values_ZG_period = [1440 1440 1100];
values_ZG_photoperiod = [0.2 0.8 0.8];


scn = [];
hpa = [];
hpt = [];
hns = [];
x_ini = [];
SAVE_FLAG = 1;

% numerics
n_steps = 525600;
dt=1;						% time step (min)
t=dt*[1:n_steps];			% time array

% Initial conditions for the fast variables
x_ini(1) = 1;	% SCN activity
x_ini(2) = 0;				% SCN velocity
x_ini(3) = 0.0001;		% CRH
x_ini(4) = 0.03;		% ACTH
x_ini(5) = 400;			% Cortisol
x_ini(9) = 0.1;			% EYA3
x_ini(10) = 0.1;		% w1
x_ini(11) = 0.1;		% Melatonin
x_ini(12) = 0.1;		% w2
x_ini(13) = 1;			% ARC
x_ini(14) = 0;			% VTA
x_ini(15) = 0;			% NAc

% HPA fixed parameters
hpa.a1 = 0.063;
hpa.a2 = 21;
hpa.a3 = 210;
hpa.k1deg = 0.03639;
hpa.k2deg = 0.09704;
hpa.k3deg = 0.02079;
hpa.betaA = 0.1;
hpa.Cort0 = 150;

% HPT fixed parameters
hpt.a1 = 0.0000117;
hpt.a2 = 0.000039;
hpt.a3 = 0.039;
hpt.k1deg = 0.000006758;
hpt.k2deg = 0.000018022;
hpt.k3deg = 0.000003862;
hpt.betaT = 0.1;
hpt.T30 = 150;
hpt.epsilon = 0.08;
hpt.bE = 0.75;
hpt.c = 0.7;
hpt.d = 0.8;
hpt.bM = 1;
hpt.FhN_tauE = 33.3333;
hpt.FhN_tauM = 66.6666;


% HNS fixed parameters
hns.a1 = 2;
hns.a2 = 10;
hns.a3 = 2;
hns.b1 = 8;
hns.b2 = 10;
hns.b3 = 8;
hns.c3 = 2;
hns.rho1 = -6;
hns.rho2 = 0;
hns.rho3 = -6;

% some parameters change value according to the condition
for NPY_counter = 1:length(conditions_NPY)
	hns.betaN = values_NPY(NPY_counter);
	for DA_counter = 1:length(conditions_DA)
		hns.tauN = values_DA(DA_counter);
		for ZG_counter = 1:length(conditions_ZG)
			disp(['Running: ' conditions_NPY{NPY_counter} ', ' conditions_DA{DA_counter} ', ' conditions_ZG{ZG_counter}]);
			scn.period = values_ZG_period(ZG_counter);
			hpt.photoperiod = values_ZG_photoperiod(ZG_counter);
			data_filename = ['timeseries_' conditions_NPY{NPY_counter} '_' conditions_DA{DA_counter} '_' conditions_ZG{ZG_counter} '.mat'];

			% Initial conditions for the slow variables
			if strcmp(conditions_NPY(NPY_counter),'NPYlo')
				if strcmp(conditions_DA(DA_counter),'DAlo')
					if strcmp(conditions_ZG(ZG_counter),'SP')
						% NPYlo,DAlo,SP initial data
						x_ini(6) = 0.015;		% TRH
						x_ini(7) = 0.025;		% TSH
						x_ini(8) = 265;			% T3
					elseif strcmp(conditions_ZG(ZG_counter),'LP')
						% NPYlo,DAlo,LP initial data
						x_ini(6) = 0.00035;	% TRH
						x_ini(7) = 0.039;		% TSH
						x_ini(8) = 362;		% T3
					else
						% NPYlo,DAlo,SH initial data
						x_ini(6) = 0.016;		% TRH
						x_ini(7) = 0.025;		% TSH
						x_ini(8) = 265;			% T3
					end
				else
					if strcmp(conditions_ZG(ZG_counter),'SP')
						% NPYlo,DAmed,SP initial data
						x_ini(6) = 0.029;		% TRH
						x_ini(7) = 0.038;		% TSH
						x_ini(8) = 265;			% T3
					elseif strcmp(conditions_ZG(ZG_counter),'LP')
						% NPYlow initial data
						x_ini(6) = 0.00035;	% TRH
						x_ini(7) = 0.039;		% TSH
						x_ini(8) = 362;		% T3
					else
						% NPYlo,DAlo,SH initial data
						x_ini(6) = 0.016;		% TRH
						x_ini(7) = 0.025;		% TSH
						x_ini(8) = 265;			% T3
					end
				end
			else
				if strcmp(conditions_DA(DA_counter),'DAlo')
					if strcmp(conditions_ZG(ZG_counter),'SP')
						% NPYhi,DAlo,SP initial data
						x_ini(6) = 0.0159;		% TRH
						x_ini(7) = 0.0159;		% TSH
						x_ini(8) = 153;		% T3
					elseif strcmp(conditions_ZG(ZG_counter),'LP')
						% NPYhi,DAlo,LP initial data
						x_ini(6) = 0.00037;	% TRH
						x_ini(7) = 0.037;		% TSH
						x_ini(8) = 276;		% T3
					else
						% NPYhi,DAlo,SH initial data
						x_ini(6) = 0.02;	% TRH
						x_ini(7) = 0.024;		% TSH
						x_ini(8) = 161;		% T3
					end
				else
					if strcmp(conditions_ZG(ZG_counter),'SP')
						% NPYhi,DAmed,SP initial data
						x_ini(6) = 0.017;		% TRH
						x_ini(7) = 0.022;		% TSH
						x_ini(8) = 189;		% T3
					elseif strcmp(conditions_ZG(ZG_counter),'LP')
						% NPYhi,DAmed,LP initial data
						x_ini(6) = 0.00055;	% TRH
						x_ini(7) = 0.05;		% TSH
						x_ini(8) = 277;		% T3
					else
						% NPYhigh,DAmed initial data
						x_ini(6) = 0.022;	% TRH
						x_ini(7) = 0.025;		% TSH
						x_ini(8) = 159;		% T3
					end
				end
			end

			options = odeset('reltol',1e-6);
			x = ode45(@(t,x)scn_modelkernel(t,x,scn,hpa,hpt,hns),t,x_ini,options);

			% Save time evolution of selected variables
			scn.t_axis = x.x;
			scn.x = x.y(1,:);
			hpa.Cort = x.y(5,:);
			hpt.T3 = x.y(8,:);
			hpt.EYA3 = x.y(9,:);
			hpt.Mel = x.y(11,:);
			hns.ARC = x.y(13,:);
			if SAVE_FLAG==1
				save(data_filename,'scn','hpa','hpt','hns');
			end


		end
	end
end



%% plot SCN-HPT, model features
% Figure 2A
% comparison between SP,LP,NPYhi (DAlo)


clear all;

sp.timeshift_t = 145;			% in days
sp.duration_t = 90;				% in days
lp.timeshift_t = 195;			% in days
lp.duration_t = 90;				% in days
NPYhiDAlo.timeshift_t = 205;			% in days
NPYhiDAlo.duration_t = 90;				% in days

% short-photoperiod data
load_filename = 'timeseries_NPYlo_DAlo_SP.mat';
load(load_filename);
sp.timeshift_n = find(scn.t_axis/60/24>sp.timeshift_t,1);
sp.duration_n = find(scn.t_axis/60/24>sp.duration_t,1);
sp.scn.t_axis = scn.t_axis(sp.timeshift_n:sp.timeshift_n+sp.duration_n)/60/24;		% in days
sp.hpt.T3 = hpt.T3(sp.timeshift_n:sp.timeshift_n+sp.duration_n);
sp.hpt.photoperiod = hpt.photoperiod;

% long-photoperiod data
load_filename = 'timeseries_NPYlo_DAlo_LP.mat';
load(load_filename);
lp.timeshift_n = find(scn.t_axis/60/24>lp.timeshift_t,1);
lp.duration_n = find(scn.t_axis/60/24>lp.duration_t,1);
lp.scn.t_axis = scn.t_axis(lp.timeshift_n:lp.timeshift_n+lp.duration_n)/60/24;
lp.hpt.T3 = hpt.T3(lp.timeshift_n:lp.timeshift_n+lp.duration_n);
lp.hpt.photoperiod = hpt.photoperiod;

% NPYhi,DAlo,LP data
load_filename = 'timeseries_NPYhi_DAlo_LP.mat';
load(load_filename);
NPYhiDAlo.timeshift_n = find(scn.t_axis/60/24>NPYhiDAlo.timeshift_t,1);
NPYhiDAlo.duration_n = find(scn.t_axis/60/24>NPYhiDAlo.duration_t,1);
NPYhiDAlo.scn.t_axis = scn.t_axis(NPYhiDAlo.timeshift_n:NPYhiDAlo.timeshift_n+NPYhiDAlo.duration_n)/60/24;
NPYhiDAlo.hpt.T3 = hpt.T3(NPYhiDAlo.timeshift_n:NPYhiDAlo.timeshift_n+NPYhiDAlo.duration_n);
NPYhiDAlo.hpt.photoperiod = hpt.photoperiod;




fsize = 12;
fig_factor = 1;
lwidth1 = 1;
lwidth2 = 3;

figure(1);
clf(1);
set(gcf,'Units','centimeters');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[10 10],'paperposition',[0 0 10 10]);
set(gca,'fontsize',fsize);

y_scale = max([sp.hpt.T3 lp.hpt.T3 NPYhiDAlo.hpt.T3]);
plot(sp.scn.t_axis-sp.scn.t_axis(1),sp.hpt.T3/y_scale,'k--','linewidth',lwidth2);
hold on;
plot(lp.scn.t_axis-lp.scn.t_axis(1),lp.hpt.T3/y_scale,'k-','linewidth',lwidth2);
plot(NPYhiDAlo.scn.t_axis-NPYhiDAlo.scn.t_axis(1),NPYhiDAlo.hpt.T3/y_scale,'k:','linewidth',lwidth2);
set(gca,'xtick',[0:30:sp.duration_t],'ytick',[0:0.25:1],'fontsize',fsize);
xlim([0 sp.duration_t]);
ylim([0 1]);
xlabel('Time (days)');
ylabel('T3 concentration (a.u.)');
legend('SP','LP','LP, fasting','location','southeast');
title('SCN-HPT-HNS interaction','fontweight','bold');


print -dpdf 'figure_2A.pdf';



%% plot SCN-HPT, model predictions
% Figure 2B
% comparison between SP,LP,SH (NPYlo,DAlo)


clear all;

sp.timeshift_t = 145;			% in days
sp.duration_t = 90;				% in days
lp.timeshift_t = 195;			% in days
lp.duration_t = 90;				% in days
sh.timeshift_t = 195;			% in days
sh.duration_t = 90;				% in days

% short-photoperiod data
load_filename = 'timeseries_NPYlo_DAlo_SP.mat';
load(load_filename);
sp.timeshift_n = find(scn.t_axis/60/24>sp.timeshift_t,1);
sp.duration_n = find(scn.t_axis/60/24>sp.duration_t,1);
sp.scn.t_axis = scn.t_axis(sp.timeshift_n:sp.timeshift_n+sp.duration_n)/60/24;		% in days
sp.hpt.T3 = hpt.T3(sp.timeshift_n:sp.timeshift_n+sp.duration_n);
sp.hpt.photoperiod = hpt.photoperiod;

% long-photoperiod data
load_filename = 'timeseries_NPYlo_DAlo_LP.mat';
load(load_filename);
lp.timeshift_n = find(scn.t_axis/60/24>lp.timeshift_t,1);
lp.duration_n = find(scn.t_axis/60/24>lp.duration_t,1);
lp.scn.t_axis = scn.t_axis(lp.timeshift_n:lp.timeshift_n+lp.duration_n)/60/24;
lp.hpt.T3 = hpt.T3(lp.timeshift_n:lp.timeshift_n+lp.duration_n);
lp.hpt.photoperiod = hpt.photoperiod;

% short-period data
load_filename = 'timeseries_NPYlo_DAlo_SH.mat';
load(load_filename);
sh.timeshift_n = find(scn.t_axis/60/24>sh.timeshift_t,1);
sh.duration_n = find(scn.t_axis/60/24>sh.duration_t,1);
sh.scn.t_axis = scn.t_axis(sh.timeshift_n:sh.timeshift_n+sh.duration_n)/60/24;
sh.hpt.T3 = hpt.T3(sh.timeshift_n:sh.timeshift_n+sh.duration_n);
sh.hpt.photoperiod = hpt.photoperiod;



fsize = 12;
fig_factor = 1;
lwidth1 = 1;
lwidth2 = 3;

figure(2);
clf(2);
set(gcf,'Units','centimeters');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[10 10],'paperposition',[0 0 10 10]);
set(gca,'fontsize',fsize);

y_scale = max([sp.hpt.T3 lp.hpt.T3 sh.hpt.T3]);
plot(sp.scn.t_axis-sp.scn.t_axis(1),sp.hpt.T3/y_scale,'k--','linewidth',lwidth2);
hold on;
plot(lp.scn.t_axis-lp.scn.t_axis(1),lp.hpt.T3/y_scale,'k-','linewidth',lwidth2);
plot(sh.scn.t_axis-sh.scn.t_axis(1),sh.hpt.T3/y_scale,'k:','linewidth',lwidth2);
set(gca,'xtick',[0:30:sp.duration_t],'ytick',[0:0.25:1],'fontsize',fsize);
xlim([0 sp.duration_t]);
ylim([0 1]);
xlabel('Time (days)');
ylabel('T3 concentration (a.u.)');
legend('SP','LP','SH','location','southeast');
title('SCN-HPT interaction','fontweight','bold');


print -dpdf 'figure_2B.pdf';


%% plot SCN-HPT-HNS, model predictions
% Figure 2C
% comparison between NPYlo/NPYhi,DAlo/DAmed (LP)


clear all;

NPYloDAlo.timeshift_t = 195;			% in days
NPYloDAlo.duration_t = 90;				% in days
NPYloDAmed.timeshift_t = 195;			% in days
NPYloDAmed.duration_t = 90;				% in days
NPYhiDAlo.timeshift_t = 205;			% in days
NPYhiDAlo.duration_t = 90;				% in days
NPYhiDAmed.timeshift_t = 250;			% in days
NPYhiDAmed.duration_t = 90;				% in days

% NPYlo,DAlo,LP data
load_filename = 'timeseries_NPYlo_DAlo_LP.mat';
load(load_filename);
NPYloDAlo.timeshift_n = find(scn.t_axis/60/24>NPYloDAlo.timeshift_t,1);
NPYloDAlo.duration_n = find(scn.t_axis/60/24>NPYloDAlo.duration_t,1);
NPYloDAlo.scn.t_axis = scn.t_axis(NPYloDAlo.timeshift_n:NPYloDAlo.timeshift_n+NPYloDAlo.duration_n)/60/24;		% in days
NPYloDAlo.hpt.T3 = hpt.T3(NPYloDAlo.timeshift_n:NPYloDAlo.timeshift_n+NPYloDAlo.duration_n);
NPYloDAlo.hpt.photoperiod = hpt.photoperiod;

% NPYlo,DAmed,LP data
load_filename = 'timeseries_NPYlo_DAmed_LP.mat';
load(load_filename);
NPYloDAmed.timeshift_n = find(scn.t_axis/60/24>NPYloDAmed.timeshift_t,1);
NPYloDAmed.duration_n = find(scn.t_axis/60/24>NPYloDAmed.duration_t,1);
NPYloDAmed.scn.t_axis = scn.t_axis(NPYloDAmed.timeshift_n:NPYloDAmed.timeshift_n+NPYloDAmed.duration_n)/60/24;
NPYloDAmed.hpt.T3 = hpt.T3(NPYloDAmed.timeshift_n:NPYloDAmed.timeshift_n+NPYloDAmed.duration_n);
NPYloDAmed.hpt.photoperiod = hpt.photoperiod;

% NPYhi,DAlo,LP data
load_filename = 'timeseries_NPYhi_DAlo_LP.mat';
load(load_filename);
NPYhiDAlo.timeshift_n = find(scn.t_axis/60/24>NPYhiDAlo.timeshift_t,1);
NPYhiDAlo.duration_n = find(scn.t_axis/60/24>NPYhiDAlo.duration_t,1);
NPYhiDAlo.scn.t_axis = scn.t_axis(NPYhiDAlo.timeshift_n:NPYhiDAlo.timeshift_n+NPYhiDAlo.duration_n)/60/24;
NPYhiDAlo.hpt.T3 = hpt.T3(NPYhiDAlo.timeshift_n:NPYhiDAlo.timeshift_n+NPYhiDAlo.duration_n);
NPYhiDAlo.hpt.photoperiod = hpt.photoperiod;

% NPYhi,DAmed,LP data
load_filename = 'timeseries_NPYhi_DAmed_LP.mat';
load(load_filename);
NPYhiDAmed.timeshift_n = find(scn.t_axis/60/24>NPYhiDAmed.timeshift_t,1);
NPYhiDAmed.duration_n = find(scn.t_axis/60/24>NPYhiDAmed.duration_t,1);
NPYhiDAmed.scn.t_axis = scn.t_axis(NPYhiDAmed.timeshift_n:NPYhiDAmed.timeshift_n+NPYhiDAmed.duration_n)/60/24;
NPYhiDAmed.hpt.T3 = hpt.T3(NPYhiDAmed.timeshift_n:NPYhiDAmed.timeshift_n+NPYhiDAmed.duration_n);
NPYhiDAmed.hpt.photoperiod = hpt.photoperiod;



fsize = 12;
fig_factor = 1;
lwidth1 = 1;
lwidth2 = 3;

figure(3);
clf(3);
set(gcf,'Units','centimeters');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[10 10],'paperposition',[0 0 10 10]);
set(gca,'fontsize',fsize);

y_scale = max([NPYloDAlo.hpt.T3 NPYloDAmed.hpt.T3 NPYhiDAlo.hpt.T3 NPYhiDAmed.hpt.T3]);
plot(NPYloDAlo.scn.t_axis-NPYloDAlo.scn.t_axis(1),NPYloDAlo.hpt.T3/y_scale,'k-','linewidth',lwidth2);
hold on;
% plot(NPYloDAmed.scn.t_axis-NPYloDAmed.scn.t_axis(1),NPYloDAmed.hpt.T3/y_scale,'k-','linewidth',lwidth2);
plot(NPYhiDAlo.scn.t_axis-NPYhiDAlo.scn.t_axis(1),NPYhiDAlo.hpt.T3/y_scale,'k:','linewidth',lwidth2);
plot(NPYhiDAmed.scn.t_axis-NPYhiDAmed.scn.t_axis(1),NPYhiDAmed.hpt.T3/y_scale,'k--','linewidth',lwidth2);
set(gca,'xtick',[0:30:NPYloDAlo.duration_t],'ytick',[0:0.25:1],'fontsize',fsize);
xlim([0 NPYloDAlo.duration_t]);
ylim([0 1]);
xlabel('Time (days)');
ylabel('T3 concentration (a.u.)');
% legend('NPYlow,DAlow','NPYlow,DAmed','NPYhigh,DAlow','NPYhigh,DAmed','location','southeast');
% legend('Fed, normal DA','Fed, intermediate DA','Fasting, normal DA','Fasting, intermediate DA','location','southeast');
legend('Fed, normal DA','Fasting, normal DA','Fasting, intermediate DA','location','southeast');
title('SCN-HPT-HNS interaction','fontweight','bold');


print -dpdf 'figure_2C.pdf';



%% plot SCN-HPA-HPT-HNS, short timescales
% Supplementary Figure 1


clear all;

load('timeseries_NPYlo_DAlo_SP.mat');

transient_t = 10*24;		% in hours
transient_n = find(scn.t_axis/60>transient_t,1);
end_t = 13*24;				% in hours
end_n = find(scn.t_axis/60>end_t,1);
t_axis_short = (scn.t_axis(transient_n:end_n)-scn.t_axis(transient_n))/60;		% in hours
scn.x_short = scn.x(transient_n:end_n);
hpa.Cort_short = hpa.Cort(transient_n:end_n);
hpt.EYA3_short = hpt.EYA3(transient_n:end_n);
hpt.Mel_short = hpt.Mel(transient_n:end_n);
hns.ARC_short = hns.ARC(transient_n:end_n);


fsize = 12;
fig_factor = 1;
lwidth = 1;

figure(4);
clf(4);
set(gcf,'Units','centimeters');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[10 20],'paperposition',[0 0 10 20]);
set(gca,'fontsize',fsize);

ptt = suptitle({'SCN-HPA-HPT-HNS interaction','Short timescales (NPYlow,DAlow,SP)'});
set(ptt,'fontweight','bold');

subplot(3,1,1);
plot(t_axis_short,scn.x_short,'k-.','linewidth',lwidth);
hold on;
plot(t_axis_short,hpa.Cort_short/max(hpa.Cort_short),'k-','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
xlabel('Time (hours)');
ylabel('CORT concentration (a.u.)');
legend('SCN','CORT','location','southwest');


subplot(3,1,2);
plot(t_axis_short,scn.x_short,'k-.','linewidth',lwidth);
hold on;
plot(t_axis_short,hns.ARC_short/max(hns.ARC_short),'k-','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
xlabel('Time (hours)');
ylabel('ARC activity (a.u.)');
legend('SCN','ARC','location','southwest');


subplot(3,1,3);
plot(t_axis_short,scn.x_short,'k-.','linewidth',lwidth);
hold on;
plot(t_axis_short,hpt.Mel_short/max(hpt.Mel_short),'k-','linewidth',lwidth);
plot(t_axis_short,hpt.EYA3_short/max(hpt.EYA3_short),'k--','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
xlabel('Time (hours)');
ylabel('Concentration (a.u.)');
legend('SCN','Mel','EYA3','location','southwest');





print -dpdf 'suppfig_1.pdf';


%% plot SCN-HPT prediction, short timescales
% Supplementary Figure 2


clear all;

transient_t = 10*24;		% in hours
end_t = 13*24;				% in hours


% LP data
load('timeseries_NPYlo_DAlo_LP.mat');
transient_n = find(scn.t_axis/60>transient_t,1);
end_n = find(scn.t_axis/60>end_t,1);
lp.t_axis_short = (scn.t_axis(transient_n:end_n)-scn.t_axis(transient_n))/60;		% in hours
lp.scn.x_short = scn.x(transient_n:end_n);
lp.hpt.EYA3_short = hpt.EYA3(transient_n:end_n);
lp.hpt.Mel_short = hpt.Mel(transient_n:end_n);
lp.yscale = max([lp.hpt.EYA3_short lp.hpt.Mel_short]);
lp.photoperiod = hpt.photoperiod;

% SP data
load('timeseries_NPYlo_DAlo_SP.mat');
transient_n = find(scn.t_axis/60>transient_t,1);
end_n = find(scn.t_axis/60>end_t,1);
sp.t_axis_short = (scn.t_axis(transient_n:end_n)-scn.t_axis(transient_n))/60;		% in hours
sp.scn.x_short = scn.x(transient_n:end_n);
sp.hpt.EYA3_short = hpt.EYA3(transient_n:end_n);
sp.hpt.Mel_short = hpt.Mel(transient_n:end_n);
sp.yscale = max([sp.hpt.EYA3_short sp.hpt.Mel_short]);
sp.photoperiod = hpt.photoperiod;

% SH data
load('timeseries_NPYlo_DAlo_SH.mat');
transient_n = find(scn.t_axis/60>transient_t,1);
end_n = find(scn.t_axis/60>end_t,1);
sh.t_axis_short = (scn.t_axis(transient_n:end_n)-scn.t_axis(transient_n))/60;		% in hours
sh.scn.x_short = scn.x(transient_n:end_n);
sh.hpt.EYA3_short = hpt.EYA3(transient_n:end_n);
sh.hpt.Mel_short = hpt.Mel(transient_n:end_n);
sh.yscale = max([sh.hpt.EYA3_short sh.hpt.Mel_short]);
sh.photoperiod = hpt.photoperiod;

y_scale = max([lp.yscale sp.yscale sh.yscale]);


fsize = 12;
fig_factor = 1;
lwidth = 1;

figure(5);
clf(5);
set(gcf,'Units','centimeters');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[10 20],'paperposition',[0 0 10 20]);
set(gca,'fontsize',fsize);

% ptt = suptitle({'SCN-HPT interaction','Short timescales (NPYlow,DAlow)'});
% set(ptt,'fontweight','bold');


subplot(3,1,1);
plot(lp.t_axis_short,((lp.scn.x_short-(1-lp.photoperiod))>0),'k:','linewidth',2*lwidth);
hold on;
plot(lp.t_axis_short,lp.scn.x_short,'k-.','linewidth',lwidth);
plot(lp.t_axis_short,lp.hpt.Mel_short/y_scale,'k-','linewidth',lwidth);
plot(lp.t_axis_short,lp.hpt.EYA3_short/y_scale,'k--','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
ylim([-0.1 1.1]);
xlabel('Time (hours)');
ylabel('Concentration (a.u.)');
legend('Photoperiod','SCN','Mel','EYA3','location','southwest');
title({'SCN-HPT interaction','Short timescales: LP'},'fontweight','bold');


subplot(3,1,2);
plot(sp.t_axis_short,((sp.scn.x_short-(1-sp.photoperiod))>0),'k:','linewidth',2*lwidth);
hold on;
plot(sp.t_axis_short,sp.scn.x_short,'k-.','linewidth',lwidth);
plot(sp.t_axis_short,sp.hpt.Mel_short/y_scale,'k-','linewidth',lwidth);
plot(sp.t_axis_short,sp.hpt.EYA3_short/y_scale,'k--','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
ylim([-0.1 1.1]);
xlabel('Time (hours)');
ylabel('Concentration (a.u.)');
title('SP','fontweight','bold');


subplot(3,1,3);
plot(sh.t_axis_short,((sh.scn.x_short-(1-sh.photoperiod))>0),'k:','linewidth',2*lwidth);
hold on;
plot(sh.t_axis_short,sh.scn.x_short,'k-.','linewidth',lwidth);
plot(sh.t_axis_short,sh.hpt.Mel_short/y_scale,'k-','linewidth',lwidth);
plot(sh.t_axis_short,sh.hpt.EYA3_short/y_scale,'k--','linewidth',lwidth);
set(gca,'xtick',[transient_t:12:end_t]-transient_t,'xticklabel',[0:12:end_t-transient_t]);
xlim([0 end_t-transient_t]);
ylim([-0.1 1.1]);
xlabel('Time (hours)');
ylabel('Concentration (a.u.)');
title('SH','fontweight','bold');


print -dpdf 'suppfig_2.pdf';


%%
