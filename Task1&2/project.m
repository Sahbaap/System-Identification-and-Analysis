%Sahba Aghajani Pedram
%Started 11/13/2016
%Project ME270A
%This file is for tasks 1 & 2

clear
close all

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

%plotting the actual system impulse response
figure(1);
subplot(311)
plot(t,u1,'b*','LineWidth',2)
ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])


subplot(312)
plot(t,y11,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])


subplot(313)
plot(t,y21,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])


figure(2);
subplot(311)
plot(t,u2,'b*','LineWidth',2)
ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])


subplot(312)
plot(t,y12,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])


subplot(313)
plot(t,y22,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])

M = cell(20);
Njj = [20,40,80,100];
Col = [[0.0,0.0,1.0];[1.0,0.0,0.0];[0.0,1.0,0.0];[0.5,0.5,0.5]];


%task 1
%part 1
%calculating Hankel matrices
for i_n=1:1:length(Njj)

    Nj = Njj(i_n);
    R = Col(i_n,:);
              for i=1:1:(2*Nj)
              M(i) = {[y11(i+mi) y12(i+mi);y21(i+mi) y22(i+mi)]};
              end


for j=1:1:Nj
    for k=1:1:Nj
       hh(j,k)=M(k+j-1); 
       hhh(j,k)=M(k+j); 
    end
end


%Calculating Hankel matrices: H20,H40,H80,and H100
if i_n == 1 
        HH20 = cell2mat(hh);
        HHH20 = cell2mat(hhh);
        s20 = svd(HH20);
        s = svd(HH20);
               else if i_n == 2
                     HH40 = cell2mat(hh);
                     HHH40 = cell2mat(hhh);
                     s40 = svd(HH40);
                     s = svd(HH40);
                              else if i_n == 3
                                        HH80 = cell2mat(hh);
                                        HHH80 = cell2mat(hhh);
                                        s80 = svd(HH80);
                                        s = svd(HH80);
                                                 else if i_n == 4 
                                                          HH100 = cell2mat(hh);
                                                          HHH100 = cell2mat(hhh);
                                                          s = svd(HH100);
                                                          s100 = svd(HH100);
                                                     end
                                  end
                 end
end




%plotting Hankel singular values in semilog format
figure(11)
semilogy(s,'.','markersize',20,'color',R)
axis([0 40 10^(-3) 1])
grid on
hold on
xlabel('singular value index');
ylabel('Hankel singular value');
end
legend('H_{20}','H_{40}','H_{80}','H_{100}')


O100 = cell(5);
O100_l = cell(5);
C100 = cell(5);
C100_r = cell(5);


%dimension of system model %we considered values from 1 to 20
N_order = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];


%Calculating system's matrices and max of its abs. value of the e-val of A
for i_ns=1:1:length(N_order)
           ns = N_order(i_ns);
           [U,S,V] = svd(HH100);
           O100(i_ns) = {U(:,1:ns)*S(1:ns,1:ns)^(0.5)};
           O100_l(i_ns) = {S(1:ns,1:ns)^(-0.5)*U(:,1:ns)'};
           C100(i_ns) = {S(1:ns,1:ns)^(0.5)*V(:,1:ns)'};
           C100_r(i_ns) = {V(:,1:ns)*S(1:ns,1:ns)^(-0.5)};
           A(i_ns) = {cell2mat(O100_l(i_ns))*HHH100*cell2mat(C100_r(i_ns))};
           BDum = cell2mat(C100(i_ns));
           B(i_ns) = {BDum(:,1:2)};
           CDum = cell2mat(O100(i_ns));
           C(i_ns) = {CDum(1:2,:)};
           evalmax(i_ns) = max(abs(eig(cell2mat(A(i_ns)))));
end


% task 1
% part 2
%calculating response for model dimensions of 6,7,10,20
for j=[6,7,10,20]

%N_order = 20;
N_order = j;
X1(1) = {zeros(N_order,1)};
X2(1) = {zeros(N_order,1)};

%Calculating system's response for given order of the system
for i_x=1:1:N
    k = i_x;
    U1in(k) = {[u1(k);0]};
    U2in(k) = {[0;u2(k)]};
    X1(k+1) = {(cell2mat(A(N_order))*cell2mat(X1(k)) + cell2mat(B(N_order))*cell2mat(U1in(k)))};
    X2(k+1) = {(cell2mat(A(N_order))*cell2mat(X2(k)) + cell2mat(B(N_order))*cell2mat(U2in(k)))};
    Y1(k) = {cell2mat(C(N_order))*cell2mat(X1(k))};
    Y2(k) = {cell2mat(C(N_order))*cell2mat(X2(k))};
end



YY1 = cell2mat(Y1);
UU1 = cell2mat(U1in);

YY2 = cell2mat(Y2);
UU2 = cell2mat(U2in);

Y11 = YY1(1,:);
Y21 = YY1(2,:);
Y12 = YY2(1,:);
Y22 = YY2(2,:);

%plotting actual vs. modeled system's response

if j==6
    figure(121)
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', '\bf{$n_{s}=6$}', ...
    'FontSize',16,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Interpreter','Latex')
end

if j==7
    figure(122)
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', '\bf{$n_{s}=7$}', ...
    'FontSize',16,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Interpreter','Latex')
end

if j==10
    figure(123)
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', '\bf{$n_{s}=10$}', ...
    'FontSize',16,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Interpreter','Latex')
end

if j==20
    figure(124)
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', '\bf{$n_{s}=20$}', ...
    'FontSize',16,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Interpreter','Latex')
end

%2,1 element 
subplot(223)
plot(t,Y21,'b*','LineWidth',2);
hold on
plot(t,y21,'ro','LineWidth',2)
axis([-0.2 2 -0.1 0.1])
legend('Simulated impulse response','Actual impulse response')
grid on

%1,1 element
subplot(221)
plot(t,Y11,'b*','LineWidth',2);
hold on
plot(t,y11,'ro','LineWidth',2)
axis([-0.2 2 -0.1 0.1])
legend('Simulated impulse response','Actual impulse response')
grid on

%1,2 element
subplot(222)
plot(t,Y12,'b*','LineWidth',2);
hold on
plot(t,y12,'ro','LineWidth',2)
axis([-0.2 2 -0.1 0.1])
legend('Simulated impulse response','Actual impulse response')
grid on

%2,2 element
subplot(224)
plot(t,Y22,'b*','LineWidth',2);
hold on
plot(t,y22,'ro','LineWidth',2)
axis([-0.2 2 -0.1 0.1])
legend('Simulated impulse response','Actual impulse response')
grid on

end

%task 1
%part 3 and 4
%calculating the frequency response of the model for various dimensions
wnyq = (2*pi)/(2*ts);
w = [0:(N-1)]*(wnyq/N);

for i=1:1:length(w)
hhhh1(i) = {cell2mat(C(6))*(exp(1i*(w(i)*ts))*eye(6)-cell2mat(A(6)))^(-1)*cell2mat(B(6))};
hhhh2(i) = {cell2mat(C(7))*(exp(1i*(w(i)*ts))*eye(7)-cell2mat(A(7)))^(-1)*cell2mat(B(7))};
hhhh3(i) = {cell2mat(C(10))*(exp(1i*(w(i)*ts))*eye(10)-cell2mat(A(10)))^(-1)*cell2mat(B(10))};
hhhh4(i) = {cell2mat(C(20))*(exp(1i*(w(i)*ts))*eye(20)-cell2mat(A(20)))^(-1)*cell2mat(B(20))};
end



HHHH(1) = {cell2mat(hhhh1)};
HHHH(2) = {cell2mat(hhhh2)};
HHHH(3) = {cell2mat(hhhh3)};
HHHH(4) = {cell2mat(hhhh4)};



%calculating phase and magnitude for Bode plot
for i=1:1:4
HHHH1 = cell2mat(HHHH(i));
a = HHHH1(2,2:2:end);
Angle11(i) = {phase(HHHH1(1,1:2:end))*(180/pi)};
Angle12(i) = {phase(HHHH1(1,2:2:end))*(180/pi)};
Angle21(i) = {phase(HHHH1(2,1:2:end))*(180/pi)};
Angle22(i) = {phase(HHHH1(2,2:2:end))*(180/pi)};
Mag11(i) = {abs(HHHH1(1,1:2:end))};
Mag12(i) = {abs(HHHH1(1,2:2:end))};
Mag21(i) = {abs(HHHH1(2,1:2:end))};
Mag22(i) = {abs(HHHH1(2,2:2:end))};
end


%calculating the frequency response of the actual system based on
%input-output data

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


%1,1 element
figure(131)
subplot(212)
semilogx(w,cell2mat(Angle11(1)),w,cell2mat(Angle11(2)),w,cell2mat(Angle11(3)),w,cell2mat(Angle11(4)),om,phase(y11f2)*(180/pi))
xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Phase$ (degree)','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on
subplot(211)
loglog(w,cell2mat(Mag11(1)),w,cell2mat(Mag11(2)),w,cell2mat(Mag11(3)),w,cell2mat(Mag11(4)),om,abs(y11f2))
%xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Magnitude$','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on

%2,1 element
figure(132)
subplot(211)
loglog(w,cell2mat(Mag12(1)),w,cell2mat(Mag12(2)),w,cell2mat(Mag12(3)),w,cell2mat(Mag12(4)),om,abs(y12f2))
%xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Magnitude$','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on
subplot(212)
semilogx(w,cell2mat(Angle12(1)),w,cell2mat(Angle12(2)),w,cell2mat(Angle12(3)),w,cell2mat(Angle12(4)),om,phase(y12f2)*(180/pi))
xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Phase$ (degree)','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on

%1,2 element
figure(133)
subplot(211)
loglog(w,cell2mat(Mag21(1)),w,cell2mat(Mag21(2)),w,cell2mat(Mag21(3)),w,cell2mat(Mag21(4)),om,abs(y21f2))
%xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Magnitude$','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on
subplot(212)
semilogx(w,cell2mat(Angle21(1)),w,cell2mat(Angle21(2)),w,cell2mat(Angle21(3)),w,cell2mat(Angle21(4)),om,phase(y21f2)*(180/pi))
xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Phase$ (degree)','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
 grid on
 
 %2,2 element
figure(134)
subplot(211)
loglog(w,cell2mat(Mag22(1)),w,cell2mat(Mag22(2)),w,cell2mat(Mag22(3)),w,cell2mat(Mag22(4)),om,abs(y22f2))
%xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Magnitude$','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on
subplot(212)
semilogx(w,cell2mat(Angle22(1)),w,cell2mat(Angle22(2)),w,cell2mat(Angle22(3)),w,cell2mat(Angle22(4)),om,(phase(y22f2))*(180/pi))
xlabel('$Frequency$ (radian/second)','FontSize',14,'Interpreter','Latex');
ylabel('$Phase$ (degree)','FontSize',14,'Interpreter','Latex');
legend('Simulated (n_{s}=6)','Simulated (n_{s}=7)','Simulated (n_{s}=10)','Simulated (n_{s}=20)','Actual')
grid on



%Task 2
%Part 1

%calculating matrices
N_order_task2 = 7;
mat_zero = zeros(N_order_task2+2);
NegativeC = -cell2mat(C(N_order_task2));
cell_zero = cell(2);
cell_zero(1) = A(N_order_task2);
cell_zero(3) = B(N_order_task2);
cell_zero(2) = {NegativeC};
cell_zero(4) = {[0 0;0 0]};
mat_zero = cell2mat(cell_zero);

L = [eye(N_order_task2) zeros(N_order_task2,2);zeros(2,N_order_task2) zeros(2,2)];
% calculating zeros of the system
ZerosOfSys = eig(mat_zero,L);
ZerosOfSys_sorted = sort(ZerosOfSys(~isinf(ZerosOfSys)),'descend'); %sorting from largest to smallest and eliminating inf zeros
fprintf('*******************\nTask 2\n*******************\n\n\n');
fprintf('*************\nPart 1\n*************\n');
fprintf('The zeros of the 7-state model are: %0.4f\n');
ZerosOfSys_sorted 



%task 2
%Part 2
%plotting the zeros and e-vals of the system
figure(21)
plot(ZerosOfSys_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
hold on
%plotting unit circle
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');

%task 2
%part 3
Evalsd = eig(cell2mat(A(N_order_task2)));
Evalsc = (1/ts)*log(Evalsd);
SSSSS = abs(Evalsc);
fprintf('*************\nPart 3\n*************\n');
fprintf('equivalent continuous-time eigenvalues are: \n');
Evalsc


%task 2
%4
b7 = cell2mat(B(N_order_task2));
c7 = -cell2mat(C(N_order_task2));

mat_zero11 = zeros(N_order_task2+1);
cell_zero11 = cell(2);

cell_zero11(1) = A(N_order_task2);
cell_zero11(3) = {b7(:,1)};
cell_zero11(2) = {c7(1,:)};
cell_zero11(4) = {0};
mat_zero11 = cell2mat(cell_zero11);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys11 = eig(mat_zero11,L);
ZerosOfSys11_sorted = sort(ZerosOfSys11(~isinf(ZerosOfSys11)),'descend'); %sorting from largest to smallest and eliminating inf zeros



%creating Hankel matrix for channel 11,12,21,22
Nj = 100;
for i=1:1:Nj
    for j=1:1:Nj
         H11(i,j) = y11(i+j+mi-1);
         H12(i,j) = y12(i+j+mi-1);
         H21(i,j) = y21(i+j+mi-1);
         H22(i,j) = y22(i+j+mi-1);
    end
end



% plotting poles/zeros for different channels
%1,1 channel
figure(241)
subplot(1,2,1)
plot(ZerosOfSys11_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')
subplot(1,2,2)
s11 = svd(H11);
semilogy(s11,'.','markersize',20);
grid on
xlabel('Singular value index');
ylabel('Hankel singular value');


%1,2 channel
mat_zero12 = zeros(N_order_task2+1);
cell_zero12 = cell(2);

cell_zero12(1) = A(N_order_task2);
cell_zero12(3) = {b7(:,2)};
cell_zero12(2) = {c7(1,:)};
cell_zero12(4) = {0};
mat_zero12 = cell2mat(cell_zero12);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys12 = eig(mat_zero12,L);
ZerosOfSys12_sorted = sort(ZerosOfSys12(~isinf(ZerosOfSys12)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(242)
subplot(1,2,1)
plot(ZerosOfSys12_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')
subplot(1,2,2)
s12 = svd(H12);
semilogy(s12,'.','markersize',20);
grid on
xlabel('Singular value index');
ylabel('Hankel singular value');




%2,1 channel
mat_zero21 = zeros(N_order_task2+1);
cell_zero21 = cell(2);

cell_zero21(1) = A(N_order_task2);
cell_zero21(3) = {b7(:,1)};
cell_zero21(2) = {c7(2,:)};
cell_zero21(4) = {0};
mat_zero21 = cell2mat(cell_zero21);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys21 = eig(mat_zero21,L);
ZerosOfSys21_sorted = sort(ZerosOfSys21(~isinf(ZerosOfSys21)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(243)
subplot(1,2,1)
plot(ZerosOfSys21_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')
subplot(1,2,2)
s21 = svd(H21);
semilogy(s21,'.','markersize',20);
grid on
xlabel('Singular value index');
ylabel('Hankel singular value');


%2,2 channel
mat_zero22 = zeros(N_order_task2+1);
cell_zero22 = cell(2);

cell_zero22(1) = A(N_order_task2);
cell_zero22(3) = {b7(:,2)};
cell_zero22(2) = {c7(2,:)};
cell_zero22(4) = {0};
mat_zero22 = cell2mat(cell_zero22);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys22 = eig(mat_zero22,L);
ZerosOfSys22_sorted = sort(ZerosOfSys22(~isinf(ZerosOfSys22)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(244)
subplot(1,2,1)
plot(ZerosOfSys22_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')

subplot(1,2,2)
s22 = svd(H22);
semilogy(s22,'.','markersize',20);
grid on
xlabel('Singular value index');
ylabel('Hankel singular value');







%task 2
%part 5
% Plotting zeros/poles for state dimension 8
N_order_task2 = 8;
b8 = cell2mat(B(N_order_task2));
c8 = -cell2mat(C(N_order_task2));

%1,1 element
mat_zero811 = zeros(N_order_task2+1);
cell_zero811 = cell(2);

cell_zero811(1) = A(N_order_task2);
cell_zero811(3) = {b8(:,1)};
cell_zero811(2) = {c8(1,:)};
cell_zero811(4) = {0};
mat_zero811 = cell2mat(cell_zero811);

L8 = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys811 = eig(mat_zero811,L8);
ZerosOfSys811_sorted = sort(ZerosOfSys811(~isinf(ZerosOfSys811)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(251)
plot(ZerosOfSys811_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')

%1,2 element
mat_zero812 = zeros(N_order_task2+1);
cell_zero812 = cell(2);

cell_zero812(1) = A(N_order_task2);
cell_zero812(3) = {b8(:,2)};
cell_zero812(2) = {c8(1,:)};
cell_zero812(4) = {0};
mat_zero812 = cell2mat(cell_zero812);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys812 = eig(mat_zero812,L8);
ZerosOfSys812_sorted = sort(ZerosOfSys812(~isinf(ZerosOfSys812)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(252)
plot(ZerosOfSys812_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')

%2,1 element
mat_zero821 = zeros(N_order_task2+1);
cell_zero821 = cell(2);

cell_zero821(1) = A(N_order_task2);
cell_zero821(3) = {b8(:,1)};
cell_zero821(2) = {c8(2,:)};
cell_zero821(4) = {0};
mat_zero821 = cell2mat(cell_zero821);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys821 = eig(mat_zero821,L);
ZerosOfSys821_sorted = sort(ZerosOfSys821(~isinf(ZerosOfSys821)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(253)
plot(ZerosOfSys821_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')


%2,2 element
mat_zero822 = zeros(N_order_task2+1);
cell_zero822 = cell(2);

cell_zero822(1) = A(N_order_task2);
cell_zero822(3) = {b8(:,2)};
cell_zero822(2) = {c8(2,:)};
cell_zero822(4) = {0};
mat_zero822 = cell2mat(cell_zero822);

L = [eye(N_order_task2) zeros(N_order_task2,1);zeros(1,N_order_task2) zeros(1,1)];
% calculating zeros of the system
ZerosOfSys822 = eig(mat_zero822,L);
ZerosOfSys822_sorted = sort(ZerosOfSys822(~isinf(ZerosOfSys822)),'descend'); %sorting from largest to smallest and eliminating inf zeros
figure(254)
plot(ZerosOfSys822_sorted,'bo')
grid on
hold on
plot(eig(cell2mat(A(N_order_task2))),'r*')
xlabel('$Real$','FontSize',14,'Interpreter','Latex');
ylabel('$Imag$','FontSize',14,'Interpreter','Latex');
hold on
theta = linspace(0,2*pi,100);
x = cos(theta);
y = sin(theta);
plot(x,y,'Color',[0,0.7,0])
axis square
axis equal
axis([-3 1.5 -2.25 2.25])
legend('Zeros','Poles','Unit circle')


