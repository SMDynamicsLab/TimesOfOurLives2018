%% Interaction among SCN, HPA, HPT, HNS
% Model parameters and plots


clear all;

% initial values
x_ini = [];

% SCN variables and parameters
scn.gamma = 1;
scn.period = 1440;%1100;%1440;%1100;		% LP and SP
% scn.period = 1100;%1440;%1100;			% SH
scn.amplitude = 1;
x_ini(1) = scn.amplitude;
x_ini(2) = 0;

% HPA variables and parameters
hpa.gamma = 0.21;%0.2;		% time scale
hpa.a1 = 0.3;%0.03;
hpa.a2 = 100;%10;
hpa.a3 = 1000;%200;
hpa.k1deg = log(2)/4;%0.173;
hpa.k2deg = 10*log(2)/15;%log(2)/15;%0.046;
hpa.k3deg = 10*log(2)/70;%log(2)/70;%0.0099;
% hpa.K = 100;%150;%100;
% hpa.Hillcoeff = 12;%16;
hpa.beta = 0.1;%1;
hpa.Cort0 = 150;
x_ini(3) = 0.0001;%0.001;		% CRH initial concentration
x_ini(4) = 0.03;%0.05;			% ACTH initial concentration
x_ini(5) = 400;%300;			% Cortisol initial concentration
% hpa.circadian_constant = 1;

% HPT variables and parameters
hpt.gamma = 0.000039;%0.0005/2;%0.0005/4;%0.0005/12;%0.0005/14;%0.2;		% time scale
hpt.a1 = 0.3;%0.03;
hpt.a2 = 100;%10;
hpt.a3 = 1000;%200;
hpt.k1deg = log(2)/4;%0.173;
hpt.k2deg = 10*log(2)/15;%log(2)/15;%0.046;
hpt.k3deg = 10*log(2)/70;%log(2)/70;%0.0099;
% hpt.K = 100;%150;%100;
% hpt.Hillcoeff = 12;%16;
hpt.beta = 0.1;%1;
hpt.T30 = 150;
% hpt.a4 = log(2)/300;
% hpt.photoperiod = 0.8;%0.2;%0.8;%0.2;		% LP and SH
hpt.photoperiod = 0.2;%0.8;%0.2;			% SP
% hpt.a5 = log(2)/300;
hpt.tau = 12.5;
hpt.b = 0.75;%5;%0.75;
hpt.c = 0.7;
hpt.d = 0.8;
% hpt.I = 0.5;
hpt.pineal = 1;%0.55;
hpt.FhN1 = 0.03;%0.065;%0.08;%0.065;
hpt.FhN2 = 0.015;%0.022;%0.035;
% % LP initial data
% x_ini(6) = 0.144;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.015;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 184;%200;%250;%200;				% T3 initial concentration
% % SP initial data
% x_ini(6) = 0.865;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.00085;%0.00084;%0.0085;%0.007;     		% TSH initial concentration
% x_ini(8) = 9;%12;%200;				% T3 initial concentration
% % SH initial data
% x_ini(6) = 0.866;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.001;%0.007;     		% TSH initial concentration
% x_ini(8) = 11;%200;				% T3 initial concentration

% % LP initial data
% x_ini(6) = 0.000004;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.045;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 500;%200;%250;%200;				% T3 initial concentration
% % SP initial data
% x_ini(6) = 0.0005;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.06;%0.00084;%0.0085;%0.007;     		% TSH initial concentration
% x_ini(8) = 350;%12;%200;				% T3 initial concentration
% % SH initial data
% x_ini(6) = 0.00015;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.01;%0.007;     		% TSH initial concentration
% x_ini(8) = 215;%200;				% T3 initial concentration

% % LP initial data
% x_ini(6) = 0.0005;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.055;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 355;%200;%250;%200;				% T3 initial concentration
% % SP initial data
% x_ini(6) = 0.015;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.025;%0.00084;%0.0085;%0.007;     		% TSH initial concentration
% x_ini(8) = 265;%12;%200;				% T3 initial concentration
% % SH initial data
% x_ini(6) = 0.016;%0.001;         	% TRH initial concentration
% x_ini(7) = 0.025;%0.007;     		% TSH initial concentration
% x_ini(8) = 265;%200;				% T3 initial concentration

% % NPYlow initial data
% x_ini(6) = 0.00035;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.039;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 362;%200;%250;%200;				% T3 initial concentration
% % NPYhigh initial data
% x_ini(6) = 0.00037;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.037;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 276;%200;%250;%200;				% T3 initial concentration
% % NPYlow,SP initial data
% x_ini(6) = 0.023;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.033;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 328;%200;%250;%200;				% T3 initial concentration
% % NPYhigh,SP initial data
% x_ini(6) = 0.019;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.018;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 152;%200;%250;%200;				% T3 initial concentration
% % NPYhigh,DAmed initial data
% x_ini(6) = 0.00048;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.047;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 309;%200;%250;%200;				% T3 initial concentration
% % NPYhigh,DAmed,SP initial data
% x_ini(6) = 0.022;%0.2;%0.04;%0.001;           % TRH initial concentration
% x_ini(7) = 0.027;%0.02;%0.04;%0.007;			% TSH initial concentration
% x_ini(8) = 202;%200;%250;%200;				% T3 initial concentration
% NPYlow,DAmedSP initial data
x_ini(6) = 0.029;%0.2;%0.04;%0.001;           % TRH initial concentration
x_ini(7) = 0.038;%0.02;%0.04;%0.007;			% TSH initial concentration
x_ini(8) = 265;%200;%250;%200;				% T3 initial concentration
x_ini(9) = 0.1;						% EYA3 initial concentration
x_ini(10) = 0.1;					% w1 initial concentration
x_ini(11) = 0.1;					% Melatonin initial concentration
x_ini(12) = 0.1;					% w2 initial concentration
% hpt.circadian_constant = 1;

% HNS variables and parameters
% hns.gamma = 0.035;%0.035;%0.0005/12;%0.0005/14;%0.2;		% time scale, DA-low (normal)
hns.gamma = 0.035/4;%0.035;%0.0005/12;%0.0005/14;%0.2;		% time scale, DA-med
hns.a1 = 2;%10;
hns.a2 = 10;%2;%10;
hns.a3 = 2;%2;%10;
hns.b1 = 8;%10;
hns.b2 = 10;%8;%10;
hns.b3 = 8;%8;
hns.c3 = 2;%2;
hns.rho1 = -6;%0;%-3;%0;
hns.rho2 = 0;%-6;
hns.rho3 = -6;%-6;
hns.beta = 0.1;%25;%0.1;%25;%0.1;		% NPY-low
% hns.beta = 25;%0.1;%25;%0.1;			% NPY-high
x_ini(13) = 1;
x_ini(14) = 0;
x_ini(15) = 0;


% numerics
n_steps = 525600;%200000;%525600;%5000000;%50000000;%35000000;%30000000;%800000;
dt=1;					% time step (min)
t=dt*[1:n_steps];			% time array



%% run


% options = odeset('reltol',1e-12);
options = odeset('reltol',1e-6);
x = ode45(@(t,x)scn_kernel(t,x,scn,hpa,hpt,hns),t,x_ini,options);



%% plot SCN-HPA

t_axis = x.x/60;
scn.x = x.y(1,:);
scn.y = x.y(2,:);
hpa.CRH = x.y(3,:);
hpa.ACTH = x.y(4,:);
hpa.Cort = x.y(5,:);


fsize = 12;
fig_factor = 1;
lwidth = 2;

figure(10);
clf(10);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HPA');
set(ptt,'fontsize',14,'fontweight','bold');


subplot(3,2,1);
plot(t_axis,scn.x/max(scn.x)*max(hpa.CRH),'k-','linewidth',lwidth);
hold on;
plot(t_axis,hpa.CRH,'c.-','linewidth',lwidth);
xlabel('Time (hours)');
ylabel('CRH concentration (nmol/L)');


subplot(3,2,3);
plot(t_axis,scn.x/max(scn.x)*max(hpa.ACTH),'k-','linewidth',lwidth);
hold on;
plot(t_axis,hpa.ACTH,'m.-','linewidth',lwidth);
xlabel('Time (hours)');
ylabel('ACTH concentration (nmol/L)');


subplot(3,2,5);
plot(t_axis,scn.x/max(scn.x)*max(hpa.Cort),'k-','linewidth',lwidth);
hold on;
plot(t_axis,hpa.Cort,'g.-','linewidth',lwidth);
xlabel('Time (hours)');
ylabel('Cortisol concentration (nmol/L)');



subplot(3,2,[2 4 6]);

hpa.CRHnull1_xaxis = [0:500];
hpa.CRHnull1_yaxis = hpa.CRHnull1(hpa.CRHnull1_xaxis);
hpa.ACTHnull_xaxis = hpa.CRHnull1_xaxis;
hpa.ACTHnull_yaxis = hpa.ACTHnull(hpa.ACTHnull_xaxis);

plot(hpa.CRHnull1_xaxis,hpa.CRHnull1_yaxis,'k--');
hold on;
plot(hpa.ACTHnull_xaxis,hpa.ACTHnull_yaxis,'k--');
plot(hpa.Cort,hpa.CRH,'m.-','linewidth',lwidth);
xlim(hpa.CRHnull1_xaxis([1 end]));
ylim([0 0.0005]);
xlabel('Cortisol concentration (nmol/L)');
ylabel('CRH concentration (nmol/L)');
% legend('hpaxnull','hpaynull');



% print -dpdf 'infeccion_suave.pdf';
% print -dpng -r150 'infeccion_suave.png';



%% plot SCN-HPT

SAVE_FLAG = 0;
% save_filename = 'timeseries_SP.mat';
% save_filename = 'timeseries_LP.mat';
% save_filename = 'timeseries_SH.mat';

% save_filename = 'timeseries_NPYlo.mat';
% save_filename = 'timeseries_NPYhi.mat';
% save_filename = 'timeseries_NPYlo_DAmed.mat';
% save_filename = 'timeseries_NPYhi_DAmed.mat';
% save_filename = 'timeseries_NPYlo_SP.mat';
% save_filename = 'timeseries_NPYhi_SP.mat';
% save_filename = 'timeseries_NPYlo_DAmed_SP.mat';
% save_filename = 'timeseries_NPYhi_DAmed_SP.mat';


transient_t = 0;%24*60;				% in minutes
% transient_t = round(n_steps/3)*dt;				% in minutes
transient_n = 1;%find(x.x>transient_t,1);
% "short" versions to show detail (first days)
short_t_end_t = transient_t + 24*60*10;						% in minutes
short_t_end_n = find(x.x>short_t_end_t,1);
t_axis_short = (x.x(transient_n:short_t_end_n)-x.x(transient_n))/60;				% in hours
scn.x_short = x.y(1,transient_n:short_t_end_n);
scn.y_short = x.y(2,transient_n:short_t_end_n);
hpt.EYA3_short = x.y(9,transient_n:short_t_end_n);
hpt.w1_short = x.y(10,transient_n:short_t_end_n);
hpt.Mel_short = x.y(11,transient_n:short_t_end_n);
hpt.w2_short = x.y(12,transient_n:short_t_end_n);
hns.ARC_short = x.y(13,transient_n:short_t_end_n);
% "long" versions to show evolution across weeks
scn.t_axis = (x.x(transient_n:end)-x.x(transient_n))/60/24;				% in days
scn.x = x.y(1,transient_n:end);
% scn.y = x.y(2,transient_n:end);
hpt.TRH = x.y(6,transient_n:end);
hpt.TSH = x.y(7,transient_n:end);
hpt.T3 = x.y(8,transient_n:end);
hpt.EYA3 = x.y(9,transient_n:end);
hpt.Mel = x.y(11,transient_n:end);
hpt.w2 = x.y(12,transient_n:end);


fsize = 12;
fig_factor = 1;
lwidth = 2;

figure(11);
clf(11);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HPT');
set(ptt,'fontsize',14,'fontweight','bold');


% short-term evolution
subplot(3,2,1);
plot(t_axis_short,((scn.x_short-(1-hpt.photoperiod))>0)*max(hpt.Mel_short),'k:','linewidth',lwidth);
hold on;
plot(t_axis_short,scn.x_short,'k-.','linewidth',lwidth);
plot(t_axis_short,hpt.Mel_short,'r-','linewidth',lwidth);
set(gca,'xtick',[0:6]*12);
xlabel('Time (hours)');
ylabel('Melatonin concentration (nmol/L ??)');
legend('Light','SCN activity','Melatonin');

subplot(3,2,3);
ptp1 = plot(t_axis_short,((scn.x_short-(1-hpt.photoperiod))>0)*max(hpt.EYA3_short),'k:','linewidth',lwidth);
hold on;
ptp2 = plot(t_axis_short,scn.x_short,'k-.','linewidth',lwidth);
ptp3 = plot(t_axis_short,hpt.w2_short,'b-','linewidth',lwidth);
ptp4 = plot(t_axis_short,hpt.EYA3_short,'g-','linewidth',lwidth);
ptp5 = plot(t_axis_short,hns.ARC_short,'r-','linewidth',lwidth);
set(gca,'xtick',[0:6]*12);
xlabel('Time (hours)');
ylabel('EYA3 concentration (nmol/L ??)');
legend([ptp3 ptp4],'w2','EYA3');


% long-term evolution
subplot(3,2,2);
% plot(scn.t_axis,scn.x/max(scn.x)*max(hpt.TRH),'k-','linewidth',lwidth);
% hold on;
plot(scn.t_axis,hpt.TRH,'b.-','linewidth',lwidth);
xlabel('Time (days)');
ylabel('TRH concentration (nmol/L)');


subplot(3,2,4);
% plot(scn.t_axis,hpt.EYA3/max(hpt.EYA3)*max(hpt.TSH),'k-','linewidth',lwidth);
% hold on;
% plot(scn.t_axis,hpt.EYA3.^8,'k:','linewidth',lwidth);
plot(scn.t_axis,hpt.TSH,'r.-','linewidth',lwidth);
xlabel('Time (days)');
ylabel('TSH concentration (nmol/L)');


subplot(3,2,6);
% plot(scn.t_axis,scn.x/max(scn.x)*max(hpt.T3),'k-','linewidth',lwidth);
% hold on;
plot(scn.t_axis,hpt.T3,'g.-','linewidth',lwidth);
xlabel('Time (days)');
ylabel('T3 concentration (nmol/L)');

% subplot(3,2,5);
% % plot(scn.t_axis,scn.x/max(scn.x)*max(hpt.T3),'k-','linewidth',lwidth);
% % hold on;
% plot3(hpt.TRH,hpt.TSH,hpt.T3,'g.-','linewidth',lwidth);
% xlabel('Time (days)');
% ylabel('T3 concentration (nmol/L)');


if SAVE_FLAG==1
	save(save_filename,'scn','hpa','hpt','hns');
end



%% plot SCN-HNS

t_axis = x.x/60;
transient_t = 200;		% in hours
transient_n = find(t_axis>transient_t);
t_axis_stationary = t_axis(transient_n:end);
scn.x = x.y(1,transient_n:end);
scn.y = x.y(2,transient_n:end);
hns.ARC = x.y(13,transient_n:end);
hns.VTA = x.y(14,transient_n:end);
hns.NAc = x.y(15,transient_n:end);


fsize = 12;
fig_factor = 1;
lwidth = 2;

figure(12);
clf(12);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HNS');
set(ptt,'fontsize',14,'fontweight','bold');


% subplot(3,2,1);
% plot(t_axis_stationary,scn.x/max(scn.x)*max(hns.ARC),'k-','linewidth',lwidth);
plot(t_axis_stationary,scn.x,'k-','linewidth',lwidth);
hold on;
plot(t_axis_stationary,hns.ARC,'c.-','linewidth',lwidth);
plot(t_axis_stationary,hns.VTA,'m.-','linewidth',lwidth);
plot(t_axis_stationary,hns.NAc,'b.-','linewidth',lwidth);
xlabel('Time (hours)');
ylabel('firing rate');
legend('SCN','ARC','VTA','NAc');



%% plot SCN-HPT
% comparison between SP,LP,SH


clear all;

% short-photoperiod data
load_filename = 'timeseries_SP.mat';
load(load_filename);
transient_t = 48;%60*12;				% in hours
short_term_duration = 3;				% in days
sp.transient_n = find(scn.t_axis*24>transient_t,1);
sp.short_t_end_t = transient_t + short_term_duration*24;						% in hours
sp.short_t_end_n = find(scn.t_axis*24>sp.short_t_end_t,1);
sp.t_axis_short = (scn.t_axis(sp.transient_n:sp.short_t_end_n)-scn.t_axis(sp.transient_n))*24;				% in hours
sp.scn.x_short = scn.x(sp.transient_n:sp.short_t_end_n);
sp.hpt.EYA3_short = hpt.EYA3(sp.transient_n:sp.short_t_end_n);
sp.hpt.Mel_short = hpt.Mel(sp.transient_n:sp.short_t_end_n);
sp.hpt.w2_short = hpt.w2(sp.transient_n:sp.short_t_end_n);
sp.scn.t_axis = scn.t_axis(sp.transient_n:end);
sp.hpt.T3 = hpt.T3(sp.transient_n:end);
sp.hpt.TSH = hpt.TSH(sp.transient_n:end);
sp.hpt.photoperiod = hpt.photoperiod;

% long-photoperiod data
load_filename = 'timeseries_LP.mat';
load(load_filename);
lp.transient_n = find(scn.t_axis*24>transient_t,1);
lp.short_t_end_t = transient_t + short_term_duration*24;						% in hours
lp.short_t_end_n = find(scn.t_axis*24>lp.short_t_end_t,1);
lp.t_axis_short = (scn.t_axis(lp.transient_n:lp.short_t_end_n)-scn.t_axis(lp.transient_n))*24;				% in hours
lp.scn.x_short = scn.x(lp.transient_n:lp.short_t_end_n);
lp.hpt.EYA3_short = hpt.EYA3(lp.transient_n:lp.short_t_end_n);
lp.hpt.Mel_short = hpt.Mel(lp.transient_n:lp.short_t_end_n);
lp.hpt.w2_short = hpt.w2(lp.transient_n:lp.short_t_end_n);
lp.scn.t_axis = scn.t_axis(lp.transient_n:end);
lp.hpt.T3 = hpt.T3(lp.transient_n:end);
lp.hpt.TSH = hpt.TSH(lp.transient_n:end);
lp.hpt.photoperiod = hpt.photoperiod;

% short-period data
load_filename = 'timeseries_SH.mat';
load(load_filename);
% transient_t = 12;%60*12;				% in hours
sh.transient_n = find(scn.t_axis*24>transient_t,1);
sh.short_t_end_t = transient_t + short_term_duration*24;						% in hours
sh.short_t_end_n = find(scn.t_axis*24>sh.short_t_end_t,1);
sh.t_axis_short = (scn.t_axis(sh.transient_n:sh.short_t_end_n)-scn.t_axis(sh.transient_n))*24;				% in hours
sh.scn.x_short = scn.x(sh.transient_n:sh.short_t_end_n);
sh.hpt.EYA3_short = hpt.EYA3(sh.transient_n:sh.short_t_end_n);
sh.hpt.Mel_short = hpt.Mel(sh.transient_n:sh.short_t_end_n);
sh.hpt.w2_short = hpt.w2(sh.transient_n:sh.short_t_end_n);
sh.scn.t_axis = scn.t_axis(sh.transient_n:end);
sh.hpt.T3 = hpt.T3(sh.transient_n:end);
sh.hpt.TSH = hpt.TSH(sh.transient_n:end);
sh.hpt.photoperiod = hpt.photoperiod;




fsize = 12;
fig_factor = 1;
lwidth1 = 1;
lwidth2 = 3;

figure(2);
clf(2);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HPT');
set(ptt,'fontsize',14,'fontweight','bold');


% short-term evolution
subplot(3,2,1);
plot(lp.t_axis_short,((lp.scn.x_short-(1-lp.hpt.photoperiod))>0)*max(lp.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(lp.t_axis_short,lp.scn.x_short,'k-.','linewidth',lwidth1);
plot(lp.t_axis_short,lp.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(lp.t_axis_short,lp.hpt.EYA3_short,'m-','linewidth',lwidth2);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('Long photoperiod (LP)','fontweight','bold');


subplot(3,2,3);
plot(sp.t_axis_short,((sp.scn.x_short-(1-sp.hpt.photoperiod))>0)*max(sp.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(sp.t_axis_short,sp.scn.x_short,'k-.','linewidth',lwidth1);
plot(sp.t_axis_short,sp.hpt.Mel_short,'c--','linewidth',lwidth2);
plot(sp.t_axis_short,sp.hpt.EYA3_short,'m--','linewidth',lwidth2);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
legend('Light','SCN activity','Melatonin','EYA3');
title('Short photoperiod (SP)','fontweight','bold');


subplot(3,2,5);
plot(sh.t_axis_short,((sh.scn.x_short-(1-sh.hpt.photoperiod))>0)*max(sh.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(sh.t_axis_short,sh.scn.x_short,'k-.','linewidth',lwidth1);
plot(sh.t_axis_short,sh.hpt.Mel_short,'c-.','linewidth',lwidth2);
plot(sh.t_axis_short,sh.hpt.EYA3_short,'m-.','linewidth',lwidth2);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('Prediction: Short period (SH)','fontweight','bold');



% long-term evolution
subplot(2,2,[2 4]);
% plot(lp.scn.t_axis,lp.hpt.T3,'b-','linewidth',lwidth2);
% hold on;
% plot(sp.scn.t_axis,sp.hpt.T3,'b--','linewidth',lwidth2);
% plot(sh.scn.t_axis,sh.hpt.T3,'b-.','linewidth',lwidth2);
plot(lp.scn.t_axis,lp.hpt.TSH,'b-','linewidth',lwidth2);
hold on;
plot(sp.scn.t_axis,sp.hpt.TSH,'b--','linewidth',lwidth2);
plot(sh.scn.t_axis,sh.hpt.TSH,'b-.','linewidth',lwidth2);
% plot(lp.scn.t_axis,lp.hpt.TSH,'r-','linewidth',lwidth2);
% hold on;
% plot(sp.scn.t_axis,sp.hpt.TSH,'r--','linewidth',lwidth2);
% plot(sh.scn.t_axis,sh.hpt.TSH,'r-.','linewidth',lwidth2);
set(gca,'fontsize',fsize);
xlabel('Time (days)');
ylabel('T3 concentration (nmol/L ??)');
legend('LP','SP','SH','location','northeast');
title('Seasonal evolution','fontweight','bold');


% print -dpdf 'scn-hpt.pdf';



%% plot SCN-HPT-HNS
% comparison between NPY-high/low, DA-med/low, SP/LP


clear all;

% NPY-low data
load_filename = 'timeseries_NPYlo.mat';
load(load_filename);
transient_t = 48;%60*12;				% in hours
short_term_duration = 3;				% in days
NPYlo.transient_n = find(scn.t_axis*24>transient_t,1);
NPYlo.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYlo.short_t_end_n = find(scn.t_axis*24>NPYlo.short_t_end_t,1);
NPYlo.t_axis_short = (scn.t_axis(NPYlo.transient_n:NPYlo.short_t_end_n)-scn.t_axis(NPYlo.transient_n))*24;				% in hours
NPYlo.scn.x_short = scn.x(NPYlo.transient_n:NPYlo.short_t_end_n);
NPYlo.hpt.EYA3_short = hpt.EYA3(NPYlo.transient_n:NPYlo.short_t_end_n);
NPYlo.hpt.Mel_short = hpt.Mel(NPYlo.transient_n:NPYlo.short_t_end_n);
NPYlo.hpt.w2_short = hpt.w2(NPYlo.transient_n:NPYlo.short_t_end_n);
NPYlo.hns.ARC_short = hns.ARC_short(NPYlo.transient_n:NPYlo.short_t_end_n);
NPYlo.scn.t_axis = scn.t_axis(NPYlo.transient_n:end);
NPYlo.hpt.T3 = hpt.T3(NPYlo.transient_n:end);
NPYlo.hpt.TSH = hpt.TSH(NPYlo.transient_n:end);
NPYlo.hpt.photoperiod = hpt.photoperiod;
NPYlo.hns.beta = hns.beta;


% NPY-high data
load_filename = 'timeseries_NPYhi.mat';
load(load_filename);
NPYhi.transient_n = find(scn.t_axis*24>transient_t,1);
NPYhi.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYhi.short_t_end_n = find(scn.t_axis*24>NPYhi.short_t_end_t,1);
NPYhi.t_axis_short = (scn.t_axis(NPYhi.transient_n:NPYhi.short_t_end_n)-scn.t_axis(NPYhi.transient_n))*24;				% in hours
NPYhi.scn.x_short = scn.x(NPYhi.transient_n:NPYhi.short_t_end_n);
NPYhi.hpt.EYA3_short = hpt.EYA3(NPYhi.transient_n:NPYhi.short_t_end_n);
NPYhi.hpt.Mel_short = hpt.Mel(NPYhi.transient_n:NPYhi.short_t_end_n);
NPYhi.hpt.w2_short = hpt.w2(NPYhi.transient_n:NPYhi.short_t_end_n);
NPYhi.hns.ARC_short = hns.ARC_short(NPYhi.transient_n:NPYhi.short_t_end_n);
NPYhi.scn.t_axis = scn.t_axis(NPYhi.transient_n:end);
NPYhi.hpt.T3 = hpt.T3(NPYhi.transient_n:end);
NPYhi.hpt.TSH = hpt.TSH(NPYhi.transient_n:end);
NPYhi.hpt.photoperiod = hpt.photoperiod;
NPYhi.hns.beta = hns.beta;

% NPY-low, DA-med data
load_filename = 'timeseries_NPYlo_DAmed.mat';
load(load_filename);
NPYloDAmed.transient_n = find(scn.t_axis*24>transient_t,1);
NPYloDAmed.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYloDAmed.short_t_end_n = find(scn.t_axis*24>NPYloDAmed.short_t_end_t,1);
NPYloDAmed.t_axis_short = (scn.t_axis(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n)-scn.t_axis(NPYloDAmed.transient_n))*24;				% in hours
NPYloDAmed.scn.x_short = scn.x(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n);
NPYloDAmed.hpt.EYA3_short = hpt.EYA3(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n);
NPYloDAmed.hpt.Mel_short = hpt.Mel(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n);
NPYloDAmed.hpt.w2_short = hpt.w2(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n);
NPYloDAmed.hns.ARC_short = hns.ARC_short(NPYloDAmed.transient_n:NPYloDAmed.short_t_end_n);
NPYloDAmed.scn.t_axis = scn.t_axis(NPYloDAmed.transient_n:end);
NPYloDAmed.hpt.T3 = hpt.T3(NPYloDAmed.transient_n:end);
NPYloDAmed.hpt.TSH = hpt.TSH(NPYloDAmed.transient_n:end);
NPYloDAmed.hpt.photoperiod = hpt.photoperiod;
NPYloDAmed.hns.beta = hns.beta;


% NPY-high, DA-med data
load_filename = 'timeseries_NPYhi_DAmed.mat';
load(load_filename);
NPYhiDAmed.transient_n = find(scn.t_axis*24>transient_t,1);
NPYhiDAmed.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYhiDAmed.short_t_end_n = find(scn.t_axis*24>NPYhiDAmed.short_t_end_t,1);
NPYhiDAmed.t_axis_short = (scn.t_axis(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n)-scn.t_axis(NPYhiDAmed.transient_n))*24;				% in hours
NPYhiDAmed.scn.x_short = scn.x(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n);
NPYhiDAmed.hpt.EYA3_short = hpt.EYA3(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n);
NPYhiDAmed.hpt.Mel_short = hpt.Mel(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n);
NPYhiDAmed.hpt.w2_short = hpt.w2(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n);
NPYhiDAmed.hns.ARC_short = hns.ARC_short(NPYhiDAmed.transient_n:NPYhiDAmed.short_t_end_n);
NPYhiDAmed.scn.t_axis = scn.t_axis(NPYhiDAmed.transient_n:end);
NPYhiDAmed.hpt.T3 = hpt.T3(NPYhiDAmed.transient_n:end);
NPYhiDAmed.hpt.TSH = hpt.TSH(NPYhiDAmed.transient_n:end);
NPYhiDAmed.hpt.photoperiod = hpt.photoperiod;
NPYhiDAmed.hns.beta = hns.beta;


% NPY-low,SP data
load_filename = 'timeseries_NPYlo_SP.mat';
load(load_filename);
NPYlo_SP.transient_n = find(scn.t_axis*24>transient_t,1);
NPYlo_SP.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYlo_SP.short_t_end_n = find(scn.t_axis*24>NPYlo_SP.short_t_end_t,1);
NPYlo_SP.t_axis_short = (scn.t_axis(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n)-scn.t_axis(NPYlo_SP.transient_n))*24;				% in hours
NPYlo_SP.scn.x_short = scn.x(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n);
NPYlo_SP.hpt.EYA3_short = hpt.EYA3(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n);
NPYlo_SP.hpt.Mel_short = hpt.Mel(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n);
NPYlo_SP.hpt.w2_short = hpt.w2(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n);
NPYlo_SP.hns.ARC_short = hns.ARC_short(NPYlo_SP.transient_n:NPYlo_SP.short_t_end_n);
NPYlo_SP.scn.t_axis = scn.t_axis(NPYlo_SP.transient_n:end);
NPYlo_SP.hpt.T3 = hpt.T3(NPYlo_SP.transient_n:end);
NPYlo_SP.hpt.TSH = hpt.TSH(NPYlo_SP.transient_n:end);
NPYlo_SP.hpt.photoperiod = hpt.photoperiod;
NPYlo_SP.hns.beta = hns.beta;


% NPY-high,SP data
load_filename = 'timeseries_NPYhi_SP.mat';
load(load_filename);
NPYhi_SP.transient_n = find(scn.t_axis*24>transient_t,1);
NPYhi_SP.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYhi_SP.short_t_end_n = find(scn.t_axis*24>NPYhi_SP.short_t_end_t,1);
NPYhi_SP.t_axis_short = (scn.t_axis(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n)-scn.t_axis(NPYhi_SP.transient_n))*24;				% in hours
NPYhi_SP.scn.x_short = scn.x(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n);
NPYhi_SP.hpt.EYA3_short = hpt.EYA3(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n);
NPYhi_SP.hpt.Mel_short = hpt.Mel(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n);
NPYhi_SP.hpt.w2_short = hpt.w2(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n);
NPYhi_SP.hns.ARC_short = hns.ARC_short(NPYhi_SP.transient_n:NPYhi_SP.short_t_end_n);
NPYhi_SP.scn.t_axis = scn.t_axis(NPYhi_SP.transient_n:end);
NPYhi_SP.hpt.T3 = hpt.T3(NPYhi_SP.transient_n:end);
NPYhi_SP.hpt.TSH = hpt.TSH(NPYhi_SP.transient_n:end);
NPYhi_SP.hpt.photoperiod = hpt.photoperiod;
NPYhi_SP.hns.beta = hns.beta;


% NPY-low,DAmed,SP data
load_filename = 'timeseries_NPYlo_DAmed_SP.mat';
load(load_filename);
NPYloDAmed_SP.transient_n = find(scn.t_axis*24>transient_t,1);
NPYloDAmed_SP.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYloDAmed_SP.short_t_end_n = find(scn.t_axis*24>NPYloDAmed_SP.short_t_end_t,1);
NPYloDAmed_SP.t_axis_short = (scn.t_axis(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n)-scn.t_axis(NPYloDAmed_SP.transient_n))*24;				% in hours
NPYloDAmed_SP.scn.x_short = scn.x(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n);
NPYloDAmed_SP.hpt.EYA3_short = hpt.EYA3(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n);
NPYloDAmed_SP.hpt.Mel_short = hpt.Mel(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n);
NPYloDAmed_SP.hpt.w2_short = hpt.w2(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n);
NPYloDAmed_SP.hns.ARC_short = hns.ARC_short(NPYloDAmed_SP.transient_n:NPYloDAmed_SP.short_t_end_n);
NPYloDAmed_SP.scn.t_axis = scn.t_axis(NPYloDAmed_SP.transient_n:end);
NPYloDAmed_SP.hpt.T3 = hpt.T3(NPYloDAmed_SP.transient_n:end);
NPYloDAmed_SP.hpt.TSH = hpt.TSH(NPYloDAmed_SP.transient_n:end);
NPYloDAmed_SP.hpt.photoperiod = hpt.photoperiod;
NPYloDAmed_SP.hns.beta = hns.beta;


% NPY-high,DAmed,SP data
load_filename = 'timeseries_NPYhi_DAmed_SP.mat';
load(load_filename);
NPYhiDAmed_SP.transient_n = find(scn.t_axis*24>transient_t,1);
NPYhiDAmed_SP.short_t_end_t = transient_t + short_term_duration*24;						% in hours
NPYhiDAmed_SP.short_t_end_n = find(scn.t_axis*24>NPYhiDAmed_SP.short_t_end_t,1);
NPYhiDAmed_SP.t_axis_short = (scn.t_axis(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n)-scn.t_axis(NPYhiDAmed_SP.transient_n))*24;				% in hours
NPYhiDAmed_SP.scn.x_short = scn.x(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n);
NPYhiDAmed_SP.hpt.EYA3_short = hpt.EYA3(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n);
NPYhiDAmed_SP.hpt.Mel_short = hpt.Mel(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n);
NPYhiDAmed_SP.hpt.w2_short = hpt.w2(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n);
NPYhiDAmed_SP.hns.ARC_short = hns.ARC_short(NPYhiDAmed_SP.transient_n:NPYhiDAmed_SP.short_t_end_n);
NPYhiDAmed_SP.scn.t_axis = scn.t_axis(NPYhiDAmed_SP.transient_n:end);
NPYhiDAmed_SP.hpt.T3 = hpt.T3(NPYhiDAmed_SP.transient_n:end);
NPYhiDAmed_SP.hpt.TSH = hpt.TSH(NPYhiDAmed_SP.transient_n:end);
NPYhiDAmed_SP.hpt.photoperiod = hpt.photoperiod;
NPYhiDAmed_SP.hns.beta = hns.beta;




% long photoperiod
fsize = 12;
fig_factor = 1;
lwidth1 = 1;
lwidth2 = 3;

figure(2);
clf(2);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HPT - HNS');
set(ptt,'fontsize',14,'fontweight','bold');


% short-term evolution
subplot(2,3,1);
plot(NPYlo.t_axis_short,((NPYlo.scn.x_short-(1-NPYlo.hpt.photoperiod))>0)*max(NPYlo.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYlo.t_axis_short,NPYlo.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYlo.t_axis_short,NPYlo.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYlo.t_axis_short,NPYlo.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYlo.t_axis_short,NPYlo.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYlo.t_axis_short,1*2./(1+exp(NPYlo.hns.beta*(NPYlo.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-low','fontweight','bold');


subplot(2,3,2);
plot(NPYhi.t_axis_short,((NPYhi.scn.x_short-(1-NPYhi.hpt.photoperiod))>0)*max(NPYhi.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYhi.t_axis_short,NPYhi.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYhi.t_axis_short,NPYhi.hpt.Mel_short,'c--','linewidth',lwidth2);
plot(NPYhi.t_axis_short,NPYhi.hpt.EYA3_short,'m--','linewidth',lwidth2);
plot(NPYhi.t_axis_short,NPYhi.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYhi.t_axis_short,1*2./(1+exp(NPYhi.hns.beta*(NPYhi.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
legend('Light','SCN activity','Melatonin','EYA3','ARC','ARC input to TRH');
title('NPY-high','fontweight','bold');


subplot(2,3,4);
plot(NPYloDAmed.t_axis_short,((NPYloDAmed.scn.x_short-(1-NPYloDAmed.hpt.photoperiod))>0)*max(NPYloDAmed.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYloDAmed.t_axis_short,NPYloDAmed.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYloDAmed.t_axis_short,NPYloDAmed.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYloDAmed.t_axis_short,NPYloDAmed.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYloDAmed.t_axis_short,NPYloDAmed.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYloDAmed.t_axis_short,1*2./(1+exp(NPYloDAmed.hns.beta*(NPYloDAmed.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-low,DA-med','fontweight','bold');


subplot(2,3,5);
plot(NPYhiDAmed.t_axis_short,((NPYhiDAmed.scn.x_short-(1-NPYhiDAmed.hpt.photoperiod))>0)*max(NPYhiDAmed.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYhiDAmed.t_axis_short,NPYhiDAmed.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYhiDAmed.t_axis_short,NPYhiDAmed.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYhiDAmed.t_axis_short,NPYhiDAmed.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYhiDAmed.t_axis_short,NPYhiDAmed.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYhiDAmed.t_axis_short,1*2./(1+exp(NPYhiDAmed.hns.beta*(NPYhiDAmed.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-high,DA-med','fontweight','bold');



% long-term evolution
subplot(2,3,[3 6]);
plot(NPYlo.scn.t_axis,NPYlo.hpt.T3,'b-','linewidth',lwidth2);
hold on;
plot(NPYhi.scn.t_axis,NPYhi.hpt.T3,'b--','linewidth',lwidth2);
plot(NPYloDAmed.scn.t_axis,NPYloDAmed.hpt.T3+10,'b:','linewidth',lwidth2);
plot(NPYhiDAmed.scn.t_axis,NPYhiDAmed.hpt.T3,'b-.','linewidth',lwidth2);
% plot(sh.scn.t_axis,sh.hpt.T3,'b-.','linewidth',lwidth2);
% plot(NPYlo.scn.t_axis,NPYlo.hpt.TSH,'r-','linewidth',lwidth2);
% hold on;
% plot(NPYhi.scn.t_axis,NPYhi.hpt.TSH,'r--','linewidth',lwidth2);
% plot(sh.scn.t_axis,sh.hpt.TSH,'r-.','linewidth',lwidth2);
set(gca,'fontsize',fsize);
xlabel('Time (days)');
ylabel('T3 concentration (nmol/L ??)');
% legend('LP','SP','SH','location','east');
% legend('NPY-low','NPY-high','NPY-low, DA-high','location','east');
legend('NPY-low','NPY-high','NPY-low,DA-med + 10','NPY-high,DA-med','location','southwest');
title('Long-term evolution','fontweight','bold');


% print -dpdf 'scn-hpt-hns.pdf';



% short photoperiod

figure(3);
clf(3);
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',fig_factor*[pos(3), pos(4)],'paperposition',fig_factor*[0 0 pos(3), pos(4)]);
set(gca,'fontsize',fsize);

ptt = suptitle('SCN - HPT - HNS');
set(ptt,'fontsize',14,'fontweight','bold');


% short-term evolution
subplot(2,3,1);
plot(NPYlo_SP.t_axis_short,((NPYlo_SP.scn.x_short-(1-NPYlo_SP.hpt.photoperiod))>0)*max(NPYlo_SP.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYlo_SP.t_axis_short,NPYlo_SP.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYlo_SP.t_axis_short,NPYlo_SP.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYlo_SP.t_axis_short,NPYlo_SP.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYlo_SP.t_axis_short,NPYlo_SP.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYlo_SP.t_axis_short,1*2./(1+exp(NPYlo_SP.hns.beta*(NPYlo_SP.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-low,SP','fontweight','bold');


subplot(2,3,2);
plot(NPYhi_SP.t_axis_short,((NPYhi_SP.scn.x_short-(1-NPYhi_SP.hpt.photoperiod))>0)*max(NPYhi_SP.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYhi_SP.t_axis_short,NPYhi_SP.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYhi_SP.t_axis_short,NPYhi_SP.hpt.Mel_short,'c--','linewidth',lwidth2);
plot(NPYhi_SP.t_axis_short,NPYhi_SP.hpt.EYA3_short,'m--','linewidth',lwidth2);
plot(NPYhi_SP.t_axis_short,NPYhi_SP.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYhi_SP.t_axis_short,1*2./(1+exp(NPYhi_SP.hns.beta*(NPYhi_SP.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
legend('Light','SCN activity','Melatonin','EYA3','ARC','ARC input to TRH');
title('NPY-high,SP','fontweight','bold');


subplot(2,3,4);
plot(NPYloDAmed_SP.t_axis_short,((NPYloDAmed_SP.scn.x_short-(1-NPYloDAmed_SP.hpt.photoperiod))>0)*max(NPYloDAmed_SP.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYloDAmed_SP.t_axis_short,NPYloDAmed_SP.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYloDAmed_SP.t_axis_short,NPYloDAmed_SP.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYloDAmed_SP.t_axis_short,NPYloDAmed_SP.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYloDAmed_SP.t_axis_short,NPYloDAmed_SP.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYloDAmed_SP.t_axis_short,1*2./(1+exp(NPYloDAmed_SP.hns.beta*(NPYloDAmed_SP.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-low,DA-med','fontweight','bold');


subplot(2,3,5);
plot(NPYhiDAmed_SP.t_axis_short,((NPYhiDAmed_SP.scn.x_short-(1-NPYhiDAmed_SP.hpt.photoperiod))>0)*max(NPYhiDAmed_SP.hpt.Mel_short),'k:','linewidth',lwidth1);
hold on;
plot(NPYhiDAmed_SP.t_axis_short,NPYhiDAmed_SP.scn.x_short,'k-.','linewidth',lwidth1);
plot(NPYhiDAmed_SP.t_axis_short,NPYhiDAmed_SP.hpt.Mel_short,'c-','linewidth',lwidth2);
plot(NPYhiDAmed_SP.t_axis_short,NPYhiDAmed_SP.hpt.EYA3_short,'m-','linewidth',lwidth2);
plot(NPYhiDAmed_SP.t_axis_short,NPYhiDAmed_SP.hns.ARC_short,'g-','linewidth',lwidth1);
plot(NPYhiDAmed_SP.t_axis_short,1*2./(1+exp(NPYhiDAmed_SP.hns.beta*(NPYhiDAmed_SP.hns.ARC_short-0))),'k-','linewidth',lwidth1);
set(gca,'xtick',[0:6]*12,'fontsize',fsize);
ylim([0 2]);
xlabel('Time (hours)');
ylabel('Melatonin, EYA3 (nmol/L ??)');
% legend('Light','SCN activity','Melatonin','EYA3');
title('NPY-high,DA-med,SP','fontweight','bold');



% long-term evolution
subplot(2,3,[3 6]);
plot(NPYlo_SP.scn.t_axis,NPYlo_SP.hpt.T3,'b-','linewidth',lwidth2);
hold on;
plot(NPYhi_SP.scn.t_axis,NPYhi_SP.hpt.T3,'b--','linewidth',lwidth2);
plot(NPYloDAmed_SP.scn.t_axis,NPYloDAmed_SP.hpt.T3+10,'b:','linewidth',lwidth2);
plot(NPYhiDAmed_SP.scn.t_axis,NPYhiDAmed_SP.hpt.T3,'b-.','linewidth',lwidth2);
set(gca,'fontsize',fsize);
xlabel('Time (days)');
ylabel('T3 concentration (nmol/L ??)');
% legend('LP','SP','SH','location','east');
% legend('NPY-low','NPY-high','NPY-low, DA-high','location','east');
legend('NPY-low','NPY-high,SP','NPY-low,DA-med + 10','NPY-high,DA-med','location','northeast');
title('Long-term evolution, SP','fontweight','bold');


% print -dpdf 'scn-hpt-hns.pdf';



%%

