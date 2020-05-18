clear all
clc
close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

data = readtable('vehicles_data.csv','ReadVariableNames',true);
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
xlabel('$t$')
ylabel('$x$')
zlabel('y')


n = numel(x);
t_max = max(t);
t_min = min(t);
x_max = max(x);
x_min = min(x);
t_range = t_max - t_min;
x_range = x_max - x_min;

gamma = 2;

m_list_t = t_range/gamma^2;
m_list_x = x_range/gamma^2;

M = 10;



for i= 1:M-1
        m_list_t = [m_list_t(1)/gamma m_list_t];
        m_list_x = [m_list_x(1)/gamma m_list_x];
end


avar = zeros(numel(m_list_t),numel(m_list_x));


%%




for m_t_counter= 1:numel(m_list_t)
for m_x_counter= 1:numel(m_list_x)
        
    m_t = m_list_t(m_t_counter);
    m_x = m_list_x(m_x_counter);

    grid_length_t = round(t_range/m_t);
    grid_length_x = round(x_range/m_x);
    cumsum = zeros(grid_length_t,grid_length_x);
    counter = zeros(grid_length_t,grid_length_x);
    for i= 1:n
        
        basket_t = floor((t(i)-t_min)/m_t) + 1;
        basket_x = floor((x(i)-x_min)/m_x) + 1;
        
        basket_t = min(basket_t,grid_length_t);
        basket_x = min(basket_x,grid_length_x);
        
        cumsum(basket_t,basket_x) = cumsum(basket_t,basket_x) + y(i);
        counter(basket_t,basket_x) = counter(basket_t,basket_x) + 1;
    end
    

    y_bar = cumsum ./ counter;


    
    
    
    
    c = 0;
    E = 0;
    for i_t= 1:grid_length_t-1
    for i_x= 1:grid_length_x-1
        
        if ~isnan(y_bar(i_t,i_x))
            if ~isnan(y_bar(i_t+1,i_x))
                E = E + (y_bar(i_t,i_x) -  y_bar(i_t+1,i_x))^2;
                c  = c + 1;
            end
            if ~isnan(y_bar(i_t,i_x+1))
                E = E + (y_bar(i_t,i_x) -  y_bar(i_t,i_x+1))^2;
                c  = c + 1;
            end
        end
    end
    end

    E = E/c;
    
    avar(m_t_counter,m_x_counter) = E;
    
    
%     if m_t_counter == 10 && m_x_counter == 10
%         figure(3)
%         heatmap(y_bar)
%     end
    
end
end


figure(2)
[T, X] = meshgrid(m_list_t,m_list_x);
surf(T, X, avar)
xlabel('$m_x$')
ylabel('$m_t$')
zlabel('variance')
set(gca,'xscale','log')
set(gca,'yscale','log')

grid on


figure(1)
draw_block([pi, 5, 10], m_list_t(8),m_list_x(8),20)




