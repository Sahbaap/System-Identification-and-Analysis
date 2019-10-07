% Sahba Aghajani Pedram
% Task 3
% Part 1,2,3


clc
close all
clear all

load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts=1/40;
N=length(y1);
t=[0:N-1]*ts-1;

for k=-200:1:200
Ruu_cell(k+201)={[0 0;0 0]};
Ryu_cell(k+201)={[0 0;0 0]};
T(k+201)=0;
Ruu11(k+201)=0;
Ruu12(k+201)=0;
Ruu21(k+201)=0;
Ruu22(k+201)=0;
end

for k=-200:1:200
    T(k+201) = k/40;
    if k<0
    p = -100+(N-1+k)/2;
    end
    if k>=0
        p = -100+(N-1)/2;
    end
    
    
for q=-p:1:p-k
Ruu_cell(k+201) = {[u1(q+p+k+201);u2(q+p+k+201)]*[u1(q+p+201);u2(q+p+201)]' + cell2mat(Ruu_cell(k+201))};
Ryu_cell(k+201) = {[y1(q+p+k+201);y2(q+p+k+201)]*[u1(q+p+201);u2(q+p+201)]' + cell2mat(Ruu_cell(k+201))};
end
A = cell2mat(Ruu_cell(k+201));
Ruu11(k+201) = A(1,1)/(2*p);
Ruu12(k+201) = A(1,2)/(2*p);
Ruu21(k+201) = A(2,1)/(2*p);
Ruu22(k+201) = A(2,2)/(2*p);
end


fprintf('*******************\nTask 3\n*******************\n\n\n');

%part 1
fprintf('*************\nPart 1\n*************\n');
fprintf('The mean for first input signal is %f\n',mean(u1));
fprintf('The mean for second input signal is %f\n\n',mean(u2));



%part 2
figure(32)
subplot(2,2,1)
plot(T,Ruu11,'.')
xlabel('$second$','FontSize',14,'Interpreter','Latex');
ylabel('$R_{u,u}(1,1)$','FontSize',14,'Interpreter','Latex');
axis([-5 5 -1 5])
grid on
subplot(2,2,2)
plot(T,Ruu12,'.')
xlabel('$second$','FontSize',14,'Interpreter','Latex');
ylabel('$R_{u,u}(1,2)$','FontSize',14,'Interpreter','Latex');
grid on
subplot(2,2,3)
plot(T,Ruu21,'.')
xlabel('$second$','FontSize',14,'Interpreter','Latex');
ylabel('$R_{u,u}(2,1)$','FontSize',14,'Interpreter','Latex');
grid on
subplot(2,2,4)
plot(T,Ruu22,'.')
xlabel('$second$','FontSize',14,'Interpreter','Latex');
ylabel('$R_{u,u}(2,2)$','FontSize',14,'Interpreter','Latex');
axis([-5 5 -1 5])
grid on


%part 3
fprintf('*************\nPart 3\n*************\n')
fprintf('The auto-corelation matrix for input u with k=0 is:\n')
cell2mat(Ruu_cell(201))/(2*p)