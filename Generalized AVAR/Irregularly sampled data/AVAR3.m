function avar = AVAR3(data_t,data_x,data_y,tau_list_t,tau_list_x)

%% FUNCTION INPUTS %%
% data_t: data time stamps
% data_x: data values
% tau_list: list of window lengths which AVAR needs to be evaluated with 
%% FUNCTION OUTPUTS %%
% avar: Allan variance of the data evaluated with each window length in tau_list 


%%
avar = zeros(numel(tau_list_t),numel(tau_list_x));
t_min = data_t(1);
t_max = data_t(end);
x_min = data_x(1);
x_max = data_x(end);


% resulation of sliding integration
dt = 5*mean(diff(data_t));
dx = 6*mean(diff(data_x));

% for each window length (tau) in tau_list
for tau_indx_x= 1:numel(tau_list_x)
for tau_indx_t= 1:numel(tau_list_t)
    tau_x = tau_list_x(tau_indx_x);
    tau_t = tau_list_t(tau_indx_t);
    total_weights = 0;
    E = 0;

    % for each sliding time t
    for x= x_min+tau_x:dx:x_max-tau_x
    for t= t_min+tau_t:dt:t_max-tau_t

        
        % extract data points which fall into the two adjacent windows 1 and 2 
        y_1 = data_y(t>data_t & data_t>t-tau_t & x>data_x & data_x>x-tau_x);
        y_2 = data_y(t+tau_t>data_t & data_t>t & x>data_x & data_x>x-tau_x);
        y_3 = data_y(t>data_t & data_t>t-tau_t & x+tau_x>data_x & data_x>x);
        
        % count how many points each windows contains
        c_1 = numel(y_1);
        c_2 = numel(y_2);
        c_3 = numel(y_3);
        
        % calculate the weight
        weight_12 = c_1 * c_2;
        weight_13 = c_1 * c_3;
        
        % compute average values of each window
        y_bar_1 = mean(y_1,'all');
        y_bar_2 = mean(y_2,'all');
        y_bar_3 = mean(y_3,'all');
        
        % if weight is nonzero then
        if ~weight_12==0
%                 weight = 1;
                %  add the weighted squared difference of averages to E
                E = E + weight_12*(y_bar_1 -  y_bar_2)^2;
                % keep track of total weights
                total_weights  = total_weights + weight_12;
        end
        if ~weight_13==0
%                 weight = 1;
                %  add the weighted squared difference of averages to E
                E = E + weight_13*(y_bar_1 -  y_bar_3)^2;
                % keep track of total weights
                total_weights  = total_weights + weight_13;
        end
    end
    end
    % normalize expected value with respect of the weights
    E = 0.5*E/total_weights;
    avar(tau_indx_t,tau_indx_x) = E;
end
end
end




