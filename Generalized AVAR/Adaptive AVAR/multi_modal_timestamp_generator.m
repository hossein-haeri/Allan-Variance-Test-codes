function timestamps = multi_modal_timestamp_generator(num_simul_timesteps, num_timestamps,num_modes)

mu = rand(1,num_modes);
std = rand(1,num_modes)/num_modes*0.5;
t = (0:1/num_simul_timesteps:1);
y = ones(1,numel(t))/8;
for i= 1:num_modes
   y = y + exp(-((t-mu(i))./std(i)).^2);
end
y = y./sum(y);

s = RandStream('mlfg6331_64');
timestamps = sort(datasample(s,t,num_timestamps,'Weights',y));
end

