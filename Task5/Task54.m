%Task 5
%Part 4

clear all
clc

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
D = [0 0;0 0];
low_boundd = 0.05;
upp_boundd = 20;
told = 0.0001;

load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs

load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets

t = [0:N-1]*ts - 1;



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
D = [0 0;0 0];

%ts=1/40;
%N=401;
wnyq = (2*pi)/(2*ts);
w = [0:(N-1)]*(wnyq/N);

y11f = fft(y11)./fft(u1);
y11f2 = y11f(1:(N-1)/2 + 1);
y12f = fft(y12)./fft(u2);
y12f2 = y12f(1:(N-1)/2 + 1);
y21f = fft(y21)./fft(u1);
y21f2 = y21f(1:(N-1)/2 + 1);
y22f = fft(y22)./fft(u2);
y22f2 = y22f(1:(N-1)/2 + 1);
N = length(y11f);
om = (2*pi)*[0:(N-1)/2]/(ts*N); %%%% frequency vector in hertz

for i=1:1:length(om)
tt_emp = [y11f2(i) y12f2(i);y21f2(i) y22f2(i)];
a_emp = svd(tt_emp);
sing1_emp(i) = a_emp(1);
sing2_emp(i) = a_emp(2);
end


for i=1:1:length(w)
    


h(i) = {C*(exp(1i*(w(i)*ts))*eye(7)-A)^(-1)*B};
tt_model = C*(exp(1i*(w(i)*ts))*eye(7)-A)^(-1)*B;

a_model = svd(tt_model);


sing1_model(i) = a_model(1);
sing2_model(i) = a_model(2);



h11(i) = tt_model(1,1);
h12(i) = tt_model(1,2);
h21(i) = tt_model(2,1);
h22(i) = tt_model(2,2);
end


%emperical
[Hd_inf,w0_dinf] = Disc_Hinf(A,B,C,D,low_boundd,upp_boundd,told,1/40);


G = abs(h11);
%loglog(w,G)
loglog(w,sing1_model)
hold on
loglog(w,sing2_model)
hold on
loglog(om,sing1_emp)
hold on
loglog(om,sing2_emp)
hold on 
plot(w0_dinf,Hd_inf,'r.','markersize',15)
xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('Singular values','FontSize',14,'Interpreter','Latex');
legend('7-state Model response','7-state Model response','Emperical frequency response data','Emperical frequency response data')
grid on


