function calcAxonAlongTracts(in_6000, in_30000, grads_6000, grads_30000, out_name)

% Read powder averaged signal along tracts
data_6000 = importdata(in_6000);
data_30000 = importdata(in_30000);

% read gradients along tracts
g_6000 = importdata(grads_6000);
g_30000 = importdata(grads_30000);

% sequence parameters
delta = [15, 15];
Delta = [29.25, 29.25];

ar = zeros(size(data_6000),'like',data_6000);
beta = zeros(size(data_6000),'like',data_6000);

disp('Starting axon diameter calculation');
tic
for i = 1:numel(data_6000)
    data_vx = [data_6000(i), data_30000(i)];
    g_vx = [g_6000(i), g_30000(i)];
    if all(data_vx)
        [ar(i), beta(i)] = getAxonRadius(delta,Delta,g_vx,data_vx,'Neumann');
    end
end
toc

[path, name] = fileparts(in_6000);
savename = fullfile(path,out_name);
writematrix(ar, savename);

end