%Task 5
%Part 1
%Function for calculating H-inf for continous-time linear dynamic system


function [H_inf,w_0] = Cont_Hinf(A,B,C,D,low_bound,upp_bound,tol)

size_A = size(A,1);
size_D = size(D,2);
cl = 0;

K = (upp_bound-low_bound)/tol;
for i=1:1:K
   A_clp11(i) = {zeros(size(A,1))};
   A_clp12(i) = {zeros(size(B,1))};
   A_clp21(i) = {zeros(size(C,2))};
   A_clp22(i) = {zeros(size(A,1))};
   A_clp(i) = {zeros(size(A,1)+size(C,2),size(A,1)+size(B,1))};
   Aclp_eval(i) = {zeros(size(A,1)+size(C,2),1)};
   real_eval(i) = 0;
end

low_bound = max(svd(D));
for gamma=low_bound+0.01:0.005:upp_bound
    
    gamma;
    cl = cl+1;
    D_gamma = (gamma^2)*eye(size_D)-(D'*D);

A_clp11(cl) = {A + B*D_gamma^(-1)*D'*C};
A_clp12(cl) = {-B*D_gamma^(-1)*B'};
A_clp21(cl) = {C'*C + (C'*D*D_gamma^(-1)*D'*C)};
A_clp22(cl) = {-A'-(C'*D*D_gamma^(-1)*B')};

A_clp(cl) = {[cell2mat(A_clp11(cl)) cell2mat(A_clp12(cl));cell2mat(A_clp21(cl)) cell2mat(A_clp22(cl))]};
cell2mat(A_clp(cl));
Aclp_eval(cl) = {eig(cell2mat(A_clp(cl)))};

real_eval(cl) = prod(real(cell2mat(Aclp_eval(cl))));
real_eval1(cl) = {real(cell2mat(Aclp_eval(cl)))};
real(cell2mat(Aclp_eval(cl)));

if cl>1
    if (min(abs(real(cell2mat(Aclp_eval(cl-1)))))<tol) && (min(abs(real(cell2mat(Aclp_eval(cl)))))>tol)
        H_inf = gamma;
        a = min(abs(real(cell2mat(Aclp_eval(cl)))));
        b = cell2mat(Aclp_eval(cl));
        for i = 1:length(b)
            if abs(real(b(i)))-a<=0.001
               w_0(i) = abs(imag(b(i)));
            end
        end
    end
end


end



end


% function [r,v] = Disc_Hinf(Ad,Bd,Cd,Dd,low_boundd,upp_boundd,told)
% Ac = -(eye(size(Ad))+Ad)^(-1)*(eye(size(Ad))-Ad);
% Bc = sqrt(2)*(eye(size(Ad))+Ad)^(-1)*Bd;
% Cc = sqrt(2)*Cd*(eye(size(Ad))+Ad)^(-1);
% Dc = Dd - Cd*(eye(size(Ad))+Ad)^(-1)*Bd;
% 
%     function [H_inf,w_0] = Cont_Hinf(Ac,Bc,Cc,Dc,low_boundd,upp_boundd,told)
%     end
% 
% 
% end
