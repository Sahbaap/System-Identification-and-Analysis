
%Task 3
%Part 4
clc
close all
clear all


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

load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;

ts=1/40;
N=length(y1);
t=[0:N-1]*ts-1;
Neg_secs=-0.1;
Pos_secs=0.1;
Ave_time = ((Pos_secs/ts)-(Neg_secs/ts))/2;
for k=(Neg_secs/ts):1:(Pos_secs/ts)
Ruu_cell(k+Ave_time+1)={[0 0;0 0]};
Ryu_cell(k+Ave_time+1)={[0 0;0 0]};
T(k+Ave_time+1)=0;
Ruu11(k+Ave_time+1)=0;
Ruu12(k+Ave_time+1)=0;
Ruu21(k+Ave_time+1)=0;
Ruu22(k+Ave_time+1)=0;
end


for k=(Neg_secs/ts):1:(Pos_secs/ts)
    T(k+Ave_time+1) = k/40;
    if k<0
    p = (Neg_secs/(2*ts))+(N-1+k)/2;
    end
    if k>=0
        p = -(Pos_secs/(2*ts))+(N-1)/2;
    end
    
    
for q=-p:1:p-k
Ruu_cell(k+Ave_time+1) = {[u1(q+p+k+Ave_time+1);u2(q+p+k+Ave_time+1)]*[u1(q+p+Ave_time+1);u2(q+p+Ave_time+1)]' + cell2mat(Ruu_cell(k+Ave_time+1))};
Ryu_cell(k+Ave_time+1) = {[y1(q+p+k+Ave_time+1);y2(q+p+k+Ave_time+1)]*[u1(q+p+Ave_time+1);u2(q+p+Ave_time+1)]' + cell2mat(Ryu_cell(k+Ave_time+1))};
end
A = cell2mat(Ruu_cell(k+Ave_time+1));
B = cell2mat(Ryu_cell(k+Ave_time+1));
Ruu11(k+Ave_time+1) = A(1,1)/(2*p);
Ruu12(k+Ave_time+1) = A(1,2)/(2*p);
Ruu21(k+Ave_time+1) = A(2,1)/(2*p);
Ruu22(k+Ave_time+1) = A(2,2)/(2*p);

Ryu11(k+Ave_time+1) = B(1,1)/(2*p);
Ryu12(k+Ave_time+1) = B(1,2)/(2*p);
Ryu21(k+Ave_time+1) = B(2,1)/(2*p);
Ryu22(k+Ave_time+1) = B(2,2)/(2*p);
end

%task 3
%Part 4

%1,1 element
figure(34)
subplot(2,2,1)
plot(T,Ryu11/max(Ruu11),'b.','Markersize',15)
hold on 
plot(t1,y11,'r.','Markersize',15)
axis([-0.2 2 -0.1 0.1])
grid on
xlabel('$time (s)$','FontSize',14,'Interpreter','Latex');
ylabel('$y_1 (u_1)$','FontSize',14,'Interpreter','Latex');
legend('System''s Pulse response from Cross-Correlation','Actual System''s pulse response')

%1,2 element
subplot(2,2,2)
plot(T,Ryu12/max(Ruu22),'b.','Markersize',15)
hold on 
plot(t1,y12,'r.','Markersize',15)
axis([-0.2 2 -0.1 0.1])
grid on
xlabel('$time (s)$','FontSize',14,'Interpreter','Latex');
ylabel('$y_1 (u_2)$','FontSize',14,'Interpreter','Latex');
legend('System''s Pulse response from Cross-Correlation','Actual System''s pulse response')

%2,1 element
subplot(2,2,3)
plot(T,Ryu21/max(Ruu11),'b.','Markersize',15)
hold on 
plot(t1,y21,'r.','Markersize',15)
axis([-0.2 2 -0.1 0.1])
grid on
xlabel('$time (s)$','FontSize',14,'Interpreter','Latex');
ylabel('$y_2 (u_1)$','FontSize',14,'Interpreter','Latex');
legend('System''s Pulse response from Cross-Correlation','Actual System''s pulse response')

%2,2 element
subplot(2,2,4)
plot(T,Ryu22/max(Ruu22),'b.','Markersize',15)
hold on 
plot(t1,y22,'r.','Markersize',15)
axis([-0.2 2 -0.1 0.1])
grid on
xlabel('$time (s)$','FontSize',14,'Interpreter','Latex');
ylabel('$y_2 (u_2)$','FontSize',14,'Interpreter','Latex');
legend('System''s Pulse response from Cross-Correlation','Actual System''s pulse response')

