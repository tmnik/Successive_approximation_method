% Ознакомится с решением задачи синтеза оптимального управления методом последовательных приближений.
% Используя алгоритмы М2-3 или М4 метода последовательных приближений, выполнить синтез оптимального управления для задачи

clc;
close all;
clear all;
%all variables
syms t T p1 p2 p3 p4 x1 x2 x3 x4 u1 u2;
 
% Исходная система в пространстве состояний
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
 
% Проверка управляемости
Co   = ctrb(A,B);
unco = length(A)-rank(Co)


X0 = [0; 3; 0.5; -0.5];   % Начальные значения
% U0 = [0; 0];            % Начальное приближение
U0 = [-1; -1];

x = [x1; x2; x3; x4];
u = [u1; u2];
p = [p1; p2; p3; p4];
 
Jk = x(1)*x(1);  %G(tk)
Fk = x(2);       %F(tk)
 
f = A*x + B*u;
n = size(f);     % n=4
 
H = 0;           % Гамильтониан
for i = 1:1:n
    H = H + p(i)*f(i);
end
fprintf('H = %s\n',char(H)); 

% Частная производная по времени
dJk_dt_par = diff(Jk,t);
dFk_dt_par = diff(Fk,t);
 
dJk_dt = dJk_dt_par;
dFk_dt = dFk_dt_par;
 
% Частная производная по х
for i=1:1:n
    diff_p(i) = - diff(H,x(i));
    dJk_dx_par(i) = diff(Jk,x(i));  
    dFk_dx_par(i) = diff(Fk,x(i));
    
    % Полная производная по времени
    dJk_dt = dJk_dt + f(i)*dJk_dx_par(i); 
    dFk_dt = dFk_dt + f(i)*dFk_dx_par(i); 
end
 fprintf('diff_P=\n'); disp(diff_p);

fprintf('Частная производная Jk по времени=\n'); disp(dJk_dt_par);
fprintf('Частная производная Fk по времени =\n'); disp(dFk_dt_par);
fprintf('Частная производная Jk по x=\n'); disp(dJk_dx_par);
fprintf('Частная производная Fk по x=\n'); disp(dFk_dx_par);
fprintf('Полная производная Jk по времени=\n'); disp(dJk_dt);
fprintf('Полная производная Fk по времени=\n'); disp(dFk_dt);

 

for i = 1:1:n
    p_T(i) = dJk_dt/dFk_dt * dFk_dx_par(i) - dJk_dx_par(i);
end
 
fprintf('p(T) = %s\n');disp(p_T);
 
syms t x1(t) x2(t) x3(t) x4(t);
x = [x1; x2; x3; x4];
diff(x) == A*x+B*U0;
[x11, x12, x13, x14] = dsolve(diff(x) == A*x+B*U0,x(0) == X0);
x1_n=[x11; x12; x13; x14]; 
 
fprintf('Решение системы = %s\n');disp(x1_n);

Fk1   = x1_n(2)==0;
s_1   = solve(Fk1,t);
sol_1 = double(s_1);
t     = sol_1(2);
fprintf('T1 = \n');disp(t);

Jk1   = t^2/4;
fprintf('J1 = \n');disp(Jk1);

p1_n =subs(-2*x1_n(1));
p2_n = subs((2*x1_n(1)*x1_n(3))/x1_n(4));
p3_n = 0;
p4_n = 0;
p_n  = [p1_n;p2_n;p3_n;p4_n];
fprintf('Значение сопряженных уравнений для первого приближения p_n = \n');
disp(p_n);


