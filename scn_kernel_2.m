%% SCN, HPA, HPT, HPG, HNS
% Model equations


function dxdt = scn_kernel_2(t,z,scn,hpa,hpt,hns)

	% SCN variables and parameters
	scn.x = z(1);		% SCN activity
	scn.y = z(2);		% SCN velocity

	% HPA variables and parameters
	hpa.CRH = z(3);		% CRH concentration
	hpa.ACTH = z(4);	% ACTH concentration
	hpa.Cort = z(5);	% Cort concentration

	% HPT variables and parameters
	hpt.TRH = z(6);		% TRH concentration
	hpt.TSH = z(7);		% TSH concentration
	hpt.T3 = z(8);		% T3 concentration
	hpt.EYA3 = z(9);	% EYA3 concentration
	hpt.w1 = z(10);		% W1 concentration
	hpt.Mel = z(11);	% Melatonin concentration
	hpt.w2 = z(12);		% W2 concentration

	% HNS variables and parameters
	hns.ARC = z(13);	% activity of ARC neurons
	hns.VTA = z(14);	% activity of VTA neurons
	hns.NAc = z(15);	% activity of NAc neurons

	dxdt = [...
			% SCN system
			scn.y;...
			-(2*pi/scn.period)^2*(scn.x-0.5);...

			% HPA system
			hpa.a1*scn.x^4/(1 + exp(hpa.betaA*(hpa.Cort-hpa.Cort0))) - hpa.k1deg*hpa.CRH;...
			hpa.a2*hpa.CRH - hpa.k2deg*hpa.ACTH;...
			hpa.a3*hpa.ACTH - hpa.k3deg*hpa.Cort;...

			% HPT system
			(hpt.a1*scn.x/(1 + exp(hpt.betaT*(hpt.T3-hpt.T30))))*2/(1 + exp(hns.betaN*(hns.ARC-0))) - hpt.k1deg*hpt.TRH;...
			hpt.a2*hpt.TRH*(hpt.EYA3^4-0) - hpt.k2deg*hpt.TSH;...
			hpt.a3*hpt.TSH - hpt.k3deg*hpt.T3;...
			(hpt.EYA3-2 - (hpt.EYA3-2)^3/3 - hpt.w1 + hpt.bE*scn.x*0.1^2/(0.1^2 + hpt.w2^2))/hpt.FhN_tauE;...
			((hpt.EYA3-2 + hpt.c - hpt.d*hpt.w1)*hpt.epsilon)/hpt.FhN_tauE;...
			(hpt.Mel-2.5 - (hpt.Mel-2.5)^3/3 - 1.5*hpt.w2 + hpt.bM*((scn.x-(1-hpt.photoperiod))<0))/hpt.FhN_tauM;...
			((hpt.Mel-2.5 + hpt.c - hpt.d*hpt.w2)*hpt.epsilon)/hpt.FhN_tauM;...

			% HNS system
			(-hns.ARC + 1/(1 + exp(-(hns.rho1*(0.75+(1-scn.x)) + hns.a1*hns.ARC + hns.b1*hns.NAc))))/hns.tauN;...
			(-hns.VTA + 1/(1 + exp(-(hns.rho2 + hns.a2*hns.VTA - hns.b2*hns.ARC))))/hns.tauN;...
			(-hns.NAc + 1/(1 + exp(-(hns.rho3 + hns.a3*hns.NAc + hns.b3*hns.VTA - hns.c3*hns.ARC))))/hns.tauN];

end


%%
