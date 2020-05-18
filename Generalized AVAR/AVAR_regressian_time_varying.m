clear all
clc
close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);

data = readtable('C:\Users\hossein_haeri\Desktop\Forgetful Databases\ring_road_data_generator/vehicles_data.csv','ReadVariableNames',true);
r = 200;                  % ring radius (m)


% x represents query point (test point)
p = (0:0.05:2*pi)';
x_lat = r*cos(p)/800000 + 40.861785;
x_lon = r*sin(p)/800000 - 77.836118;


ep = 0.1;

x = abs(table2array(data(:,{'ang_pos'})));
t = abs(table2array(data(:,{'time'})));
y = table2array(data(:,{'friction'}));

figure(1)
scatter3(t,x,y,'.')
xlabel('$t[s]$')
ylabel('$x[rad]$')
zlabel('friction')

n = numel(x);
t_max = max(t);
t_min = min(t);
x_max = max(x);
x_min = min(x);
t_range = t_max - t_min;
x_range = x_max - x_min;

gamma = 2;

m_list_t = t_range/gamma;
m_list_x = x_range/gamma;
M = 10;

for i= 1:M-1
        m_list_t = [m_list_t(1)/gamma m_list_t];
        m_list_x = [m_list_x(1)/gamma m_list_x];
end


avar = zeros(M);

%%
for i_a= 1:numel(m_list_x)
for i_b= 1:numel(m_list_t) 
    a = m_list_x(i_a);
    b = m_list_t(i_b);
    y_mean = zeros(numel(x),1);
    d = zeros(numel(x));

    for i= 1:numel(x)
        counter = 0;
        for j= 1:numel(x)
            d(i,j) = ((x(i)-x(j))/a)^2 + ((t(i)-t(j))/b)^2;
            if d(i,j) < 1
                y_mean(i) = y_mean(i) + y(j);
                counter = counter + 1;
            end
        end
        y_mean(i) = y_mean(i)/counter;
    end

    E = NaN;
    counter = 0;
    for i= 1:numel(x)
        for j= 1:numel(x)
            if d(i,j) < 1+ep && d(i,j) > 1-ep && i~=j
                if isnan(E)
                    E = (y_mean(i) - y_mean(j))^2;
                else
                E = E + (y_mean(i) - y_mean(j))^2;
                end
                counter = counter + 1;
            end
        end
    end
    avar(i_a,i_b) = E/(2*counter);

end
end


figure(2)
[X, T] = meshgrid(m_list_x,m_list_t);
surf(X, T, avar)
xlabel('$m_x$')
ylabel('$m_t$')
zlabel('variance')
% scatter3(avar)
set(gca,'xscale','log')
set(gca,'yscale','log')
% grid on



