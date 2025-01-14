close all; clear; clc;

% Trochoidal waves

lambda = 3;
g = 9.81;

k = 2*pi/lambda;
c = sqrt(g/k);

X = @(a,b,t) (a + (exp(k*b))/k*sin(k*(a+c*t)));
Y = @(a,b,t) (b - (exp(k*b))/k*cos(k*(a+c*t)));

x_mesh = 0:10;
y_mesh = 0:10;
t = 0;

map = zeros(size(x_mesh), size(y_mesh));

figure

for i=x_mesh
    for j=y_mesh
        x = X(i, j, t);
        y = Y(i, j, t);
        plot(x, y)
        hold on
    end
end

figure
plot(X(x_mesh, y_mesh, t), Y(x_mesh, y_mesh, t));