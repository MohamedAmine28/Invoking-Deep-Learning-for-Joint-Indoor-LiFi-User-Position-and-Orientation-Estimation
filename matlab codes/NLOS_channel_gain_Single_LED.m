function [H_NLOS] = NLOS_channel_gain_Single_LED(UE_position,UE_orientation,PD_position,dimension,Psi_Fov,Phi_Fov,rho,H0,m)
%-------------------------------------------------------------------------%
L = dimension(1); 
W = dimension(2); 
H = dimension(3);
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
n_1=[1;0;0];
n_2=[-1;0;0];
n_3=[0;1;0];
n_4=[0;-1;0];
%---------------------------------------------------------------------%
%                              NLoS
%---------------------------------------------------------------------%
N=50; 
%---------------------------------------------------------------------%
%                              Wall 1
%---------------------------------------------------------------------%
XY = [-L*ones(1,N+1);linspace(-W,W,N+1)];
Z = linspace(0,H,N+1);
A = (2*W/N)*(H/N);
Wall_pos = combvec(XY,Z); 
TxWall = (Wall_pos-LEDnew);
d = sqrt(sum(TxWall.^2,1));
cos_psi1 = -(n_1'*TxWall)./d;
cos_phi1 = (n_T'*TxWall)./d;
cos_phi1(cos_phi1<cos(Phi_Fov))=0;
H_TX_Wall = ((m+1)/(2*pi)*A)*1./(d.^2).*(cos_phi1.^m).*cos_psi1;
WallRx = (PD-Wall_pos);
d = sqrt(sum(WallRx.^2,1));
cos_psi2 = -(n_R'*WallRx)./d;
cos_psi2(cos_psi2<cos(Psi_Fov))=0;
cos_phi2 = (n_1'*WallRx)./d;
H_Wall_RX = H0*1./(d.^2).*(cos_phi2.^m).*cos_psi2;
H1 = rho*(H_TX_Wall.*H_Wall_RX);
H1 = sum(H1);
%---------------------------------------------------------------------%
%                              Wall 2
%---------------------------------------------------------------------%
XY = [L*ones(1,N+1);linspace(-W,W,N+1)];
Z = linspace(0,H,N+1);
A = (2*W/N)*(H/N);
Wall_pos = combvec(XY,Z); 
TxWall = (Wall_pos-LEDnew);
d = sqrt(sum(TxWall.^2,1));
cos_psi1 = -(n_2'*TxWall)./d;
cos_phi1 = (n_T'*TxWall)./d;
cos_phi1(cos_phi1<cos(Phi_Fov))=0;
H_TX_Wall = ((m+1)/(2*pi)*A)*1./(d.^2).*(cos_phi1.^m).*cos_psi1;
WallRx = (PD-Wall_pos);
d = sqrt(sum(WallRx.^2,1));
cos_psi2 = -(n_R'*WallRx)./d;
cos_psi2(cos_psi2<cos(Psi_Fov))=0;
cos_phi2 = (n_2'*WallRx)./d;
H_Wall_RX = H0*1./(d.^2).*(cos_phi2.^m).*cos_psi2;
H2 = rho*(H_TX_Wall.*H_Wall_RX);
H2 = sum(H2);
%---------------------------------------------------------------------%
%                              Wall 3
%---------------------------------------------------------------------%
XY = [linspace(-L,L,N+1);-W*ones(1,N+1)];
Z = linspace(0,H,N+1);
A = (2*L/N)*(H/N);
Wall_pos = combvec(XY,Z); 
TxWall = (Wall_pos-LEDnew);
d = sqrt(sum(TxWall.^2,1));
cos_psi1 = -(n_3'*TxWall)./d;
cos_phi1 = (n_T'*TxWall)./d;
cos_phi1(cos_phi1<cos(Phi_Fov))=0;
H_TX_Wall = ((m+1)/(2*pi)*A)*1./(d.^2).*(cos_phi1.^m).*cos_psi1;
WallRx = (PD-Wall_pos);
d = sqrt(sum(WallRx.^2,1));
cos_psi2 = -(n_R'*WallRx)./d;
cos_psi2(cos_psi2<cos(Psi_Fov))=0;
cos_phi2 = (n_3'*WallRx)./d;
H_Wall_RX = H0*1./(d.^2).*(cos_phi2.^m).*cos_psi2;
H3 = rho*(H_TX_Wall.*H_Wall_RX);
H3 = sum(H3);
%---------------------------------------------------------------------%
%                              Wall 4
%---------------------------------------------------------------------%
XY = [linspace(-L,L,N+1);W*ones(1,N+1)];
Z = linspace(0,H,N+1);
A = (2*L/N)*(H/N);
Wall_pos = combvec(XY,Z); 
TxWall = (Wall_pos-LEDnew);
d = sqrt(sum(TxWall.^2,1));
cos_psi1 = -(n_4'*TxWall)./d;
cos_phi1 = (n_T'*TxWall)./d;
cos_phi1(cos_phi1<cos(Phi_Fov))=0;
H_TX_Wall = ((m+1)/(2*pi)*A)*1./(d.^2).*(cos_phi1.^m).*cos_psi1;
WallRx = (PD-Wall_pos);
d = sqrt(sum(WallRx.^2,1));
cos_psi2 = -(n_R'*WallRx)./d;
cos_psi2(cos_psi2<cos(Psi_Fov))=0;
cos_phi2 = (n_4'*WallRx)./d;
H_Wall_RX = H0*1./(d.^2).*(cos_phi2.^m).*cos_psi2;
H4 = rho*(H_TX_Wall.*H_Wall_RX);
H4 = sum(H4);
%---------------------------------------------------------------------%
H_NLOS = H1+H2+H3+H4;
end