%task 4
%part 1, 2, 3

clc
close all
clear all

A = [0.4716   -0.6076   -0.4275    0.1885    0.0153   -0.0821    0.0162;
     0.6293   -0.2447    0.6089   -0.2044   -0.0029    0.0332   -0.0214;
    -0.4472   -0.6227    0.2729    0.2165   -0.0660    0.1804    0.0711;
     0.1570    0.1962    0.0765    0.4630   -0.2127    0.3362    0.3385;
     0.0113    0.0282    0.0643    0.2262    0.9333    0.2383   -0.0712;
    -0.0141   -0.0458   -0.1085   -0.4733    0.2502   -0.1849    0.7381;
    -0.0090    0.0242    0.1388    0.2258    0.0863   -0.7779   -0.0508];

B = [0.1793    0.2278;
    -0.1221   -0.1967;
    -0.1211    0.0072;
    -0.1914    0.1808;
     0.0117    0.0031;
    -0.1003    0.0581;
    -0.0239    0.0382];

C = [-0.2684   -0.2054    0.1079   -0.0107    0.0156   -0.0758    0.0058;
     -0.0610   -0.1264   -0.1394   -0.2236   -0.0116    0.1183   -0.0449];

load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
U1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(U1>0); %%% find index where pulse occurs

load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
U2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(U1);
y12 = y12/max(U2);
y21 = y21/max(U1);
y22 = y22/max(U2);
U1 = U1/max(U1);
U2 = U2/max(U2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets

t1 = [0:N-1]*ts - 1;
h = {[0 0;0 0]};

H2_ActualP=0;

for i=1:1:N
    h(i) = {[y11(i) y12(i);y21(i) y22(i)]};
    H2_ActualP = norm(cell2mat(h(i)),'fro')^2+H2_ActualP;
end

H2_ActualP = sqrt(H2_ActualP);

load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;

y1_scaled = 0.5*y1;
y2_scaled = 0.5*y2;
u1_scaled = 0.5*u1;
u2_scaled = 0.5*u2;

ts=1/40;
N=length(y1);
t=[0:N-1]*ts-1;
y_rms = 0;
Ryy0 = [0 0;0 0];
Go = zeros(7);
Gc= zeros(7);

p = (N-1)/2;
for q=-p:1:p
    y_rms = y_rms + trace([y1_scaled(q+p+1);y2_scaled(p+q+1)]*[y1_scaled(q+p+1);y2_scaled(p+q+1)]');
    Ryy0 = Ryy0 + [y1_scaled(q+p+1);y2_scaled(p+q+1)]*[y1_scaled(q+p+1);y2_scaled(p+q+1)]';
end

%calculating Gramians (other methods could be used
for k=1:1:30000
    Go = Go + (A')^(k-1)*C'*C*A^(k-1);
    Gc = Gc + (A)^(k-1)*B*B'*A'^(k-1);
end

y_rms = sqrt(y_rms/(2*p));
trac_Ryy0 = sqrt(trace(Ryy0)/(2*p));


fprintf('*******************\nTask 4\n*******************\n\n\n\n\n');


fprintf('*************\nPart 1\n*************\n');
fprintf('The RMS value for y is: %0.4f\n\n\n',trac_Ryy0);


NormH2_Po = sqrt(trace(B'*Go*B));
NormH2_Pc = sqrt(trace(C*Gc*C'));
fprintf('*************\nPart 2\n*************\n');
fprintf('The H2 norm caluclated from observability Gramian is: %0.4f\n',NormH2_Po);
fprintf('The H2 norm caluclated from Controlability Gramian is: %0.4f\n\n\n',NormH2_Pc);


fprintf('*************\nPart 3\n*************\n');
fprintf('The H2 norm caluclated from experimental pulse response data is: %0.4f\n',H2_ActualP);


