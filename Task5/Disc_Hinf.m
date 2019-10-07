%Task 5
%Part 2
%Function for calculating H-inf for discrete-time linear dynamic system

function [Hd_inf,w0_dinf] = Disc_Hinf(Ad,Bd,Cd,Dd,low_boundd,upp_boundd,told,ts)
Ac = -(eye(size(Ad))+Ad)^(-1)*(eye(size(Ad))-Ad);
Bc = sqrt(2)*(eye(size(Ad))+Ad)^(-1)*Bd;
Cc = sqrt(2)*Cd*(eye(size(Ad))+Ad)^(-1);
Dc = Dd - Cd*(eye(size(Ad))+Ad)^(-1)*Bd;

[H_inf,w_0] = Cont_Hinf(Ac,Bc,Cc,Dc,low_boundd,upp_boundd,told);

Hd_inf = H_inf;
for i=1:length(w_0)
bb(i) = log((1+1j*w_0(i))/(1-1j*w_0(i)))/(1j*ts); 
end
w0_dinf = max(bb);
end
