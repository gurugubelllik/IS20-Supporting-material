function [measures] = jitter_shimmer_feat(A0)

mean_Ampl = mean(A0);

% Mean absolute difference of successive cycles
measures(1) = mean(abs(diff(A0)));

% Mean absolute difference of successive cycles - expressed in percent (%)
measures(2) = 100*mean(abs(diff(A0)))/mean_Ampl;

% Perturbation quotient
[Ampl_PQ3] = perq1(A0,3);
measures(3) = Ampl_PQ3.classical_Schoentgen;
measures(4) = Ampl_PQ3.classical_Baken;
measures(5) = Ampl_PQ3.generalised_Schoentgen;

[Ampl_PQ5]=perq1(A0,5); % Use 5 cycle samples (Schoentgen)
measures(6)=Ampl_PQ5.classical_Schoentgen;
measures(7)=Ampl_PQ5.classical_Baken;
measures(8)=Ampl_PQ5.generalised_Schoentgen;

[Ampl_PQ11]=perq1(A0,11); % Use 11 cycle samples (Schoentgen)
measures(9)=Ampl_PQ11.classical_Schoentgen;
measures(10)=Ampl_PQ11.classical_Baken;
measures(11)=Ampl_PQ11.generalised_Schoentgen;

% zeroth order perturbation
measures(12) = mean(abs(A0-mean_Ampl));

% Shimmer(dB)
measures(13) = mean(20*(abs(log10((A0(1:end-1))./(A0(2:end))))));

% CV
measures(14)=mean((diff(A0)).^2)/(mean_Ampl)^2;

% TKEO
measures(15) = mean(abs(TKEO(A0)));
measures(16) = std(TKEO(A0));
Ampl_TKEO_prc = prctile(TKEO(A0),[5 25 50 75 95]);
measures(17)=Ampl_TKEO_prc(1);
measures(18)=Ampl_TKEO_prc(2);
measures(19)=Ampl_TKEO_prc(3);
measures(20)=Ampl_TKEO_prc(4);

% AM
measures(21) = (max(A0)-min(A0))/(max(A0)+min(A0));

measures(22) = Ampl_TKEO_prc(4) - Ampl_TKEO_prc(1);
end
