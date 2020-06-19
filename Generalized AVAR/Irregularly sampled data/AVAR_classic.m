function avar = AVAR_classic(tau,x,m_list)

%% FUNCTION INPUTS %%
% data_t: data time stamps
% data_x: data values
% tau_list: list of window lengths which AVAR needs to be evaluated with 
%% FUNCTION OUTPUTS %%
% avar: Allan variance of the data evaluated with each window length in tau_list 


%%
n = numel(x);
avar = zeros(numel(m_list),1);
% t_min = data_t(1);
% t_max = data_t(end);

% resulation of sliding integration
% dt = mean(diff(data_t));

% for each window length (m) in m_list
for m_indx= 1:numel(m_list)        
    m = m_list(m_indx);
    
    E = 0;
    for k= 1:n-2*m
        E = E + (x(k)-2*(k+m)+x(k+2*m)).^2;
    end
    avar(m_indx) = E ./ (tau(m_indx)^2*(n-2*m));
end

end




