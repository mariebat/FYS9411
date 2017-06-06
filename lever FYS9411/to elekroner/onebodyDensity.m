clear
close all 

filenameIN = 'positions_2el_w1_noJ'; 
suffixIN = '.txt'; 
filename = strcat(filenameIN,suffixIN); 

ID = fopen(filename); 
A = textscan(ID,'%f %f %f %f'); 

x1 = A{1}; 
y1 = A{2};
x2 = A{3};
y2 = A{4};

n = length(x1); 
rho = linspace(0,0,n); 

xmin = min(x1); 
xmax = max(x1); 
ymin = min(y1); 
ymax = max(y1); 

alpha = 0.689664;
beta = 0.399237;
omega = 1.0; 
a = 0.0; 

for i = 1:n 
    if(mod(n,i) == 0) 
        disp(i); 
    end
    %psi2 = @(x,y) (exp(-0.5*alpha*omega.*(x1(i)^2 + y1(i)^2 + x.^2 + y.^2)).*exp(a*sqrt((x1(i)-x).^2 + (y1(i) - y).^2)./(1+beta*sqrt((x1(i)-x).^2 + (y1(i) - y).^2)))).^2;
    psi2 = @(x,y) (exp(-0.5*alpha*omega.*(x1(i)^2 + y1(i)^2 + x.^2 + y.^2)).*exp(a*sqrt((x1(i)-x).^2 + (y1(i) - y).^2)./(1+beta*sqrt((x1(i)-x).^2 + (y1(i) - y).^2)))).^2;   
    rho(i) = integral2(psi2,xmin,xmax,ymin,ymax); 
end

figure(01) 
plot(x1,rho)
axis([-4 4 0 14])
title('One-body density')
xlabel('x_1') 
ylabel('\rho(x_1,y_1)')

figure(02) 
plot(y1,rho)
axis([-4 4 0 14])
title('One-body density')
xlabel('y_1') 
ylabel(' \rho(x_1,y_1)')




