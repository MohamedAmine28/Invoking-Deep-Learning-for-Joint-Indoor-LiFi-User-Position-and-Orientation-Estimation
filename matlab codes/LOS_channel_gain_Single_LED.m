function [H_LOS] = LOS_channel_gain_Single_LED(UE_position,UE_orientation,PD_position,Psi_Fov,Phi_Fov,H0,m)
%-------------------------------------------------------------------------%
UE = UE_position; 
alpha = UE_orientation(1); 
beta = UE_orientation(2);
gamma = UE_orientation(3);
PD = PD_position;
%-------------------------------------------------------------------------%
Dmtr = 7*10^-3;                               % Cellphone thickness
d1 = 0; 
d2 = 3.5*10^-2;
d3 = 5*10^-2;
%-------------------------------------------------------------------------%
n_T = [0;0;1];
LED = UE + [d3;0;Dmtr];
%---------------------------------------------------------------------%
R_alpha = [cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0; 0 0 1;];
R_beta = [1 0 0;0 cos(beta) -sin(beta);0 sin(beta) cos(beta)];
R_gamma = [cos(gamma) 0 sin(gamma);0 1 0;-sin(gamma) 0 cos(gamma)];
%---------------------------------------------------------------------%
n_T = R_alpha*R_beta*R_gamma*n_T;
LEDnew = R_alpha*R_beta*R_gamma*(LED-UE)+UE;
n_R = [0;0;-1];
%---------------------------------------------------------------------%
%                               LoS
%---------------------------------------------------------------------%
TxRxvec = (PD-LEDnew);
d = sqrt(sum(TxRxvec.^2));
cos_psi = -dot(n_R,TxRxvec)/d;
cos_psi(cos_psi<cos(Psi_Fov))=0;
cos_phi = dot(n_T,TxRxvec)/d;
cos_phi(cos_phi<cos(Phi_Fov))=0;
H_LOS = H0/(d^2).*(cos_phi^m)*cos_psi;
end