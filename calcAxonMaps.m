function calcAxonMaps(in_6000, in_30000, in_bval, in_bvec, in_grad_dev)

% Read powder averaged signal
data_6000 = niftiread(in_6000);
data_6000 = double(data_6000(:,:,:,1));
data_6000(isnan(data_6000)) = 0;
data_30000 = niftiread(in_30000);
data_30000 = double(data_30000(:,:,:,1));
data_30000(isnan(data_30000)) = 0;

% change nifti header
nifti_info = niftiinfo(in_6000);
nifti_info.ImageSize = nifti_info.ImageSize(1:3);
nifti_info.PixelDimensions = nifti_info.PixelDimensions(1:3);

% sequence parameters
bval1 = 6;
bval2 = 30.45;
gmax = 280;
delta = [15, 15];
Delta = [29.25, 29.25];
g = [sqrt(bval1/bval2)*gmax, gmax];

% correct gradient with factors from nonlinearity correction
bval = importdata(in_bval);
bvec = importdata(in_bvec);
if size(bvec,2) ~= 3
    bvec = bvec.';
end    
grad_dev = niftiread(in_grad_dev);
corr_fac = correctgradient(bval, bvec, grad_dev);
g_corr =  {corr_fac{1} * g(1), corr_fac{2} * g(2)};
g_corr_1 = g_corr{1};
g_corr_2 = g_corr{2};

% save gradient maps for profiling
[path, name] = fileparts(in_6000);
savename = fullfile(path,'grads_6000');
nifti_info.Datatype = class(g_corr_1);
niftiwrite(g_corr_1, savename, nifti_info);

savename = fullfile(path,'grads_30000');
nifti_info.Datatype = class(g_corr_2);
niftiwrite(g_corr_2, savename, nifti_info);

ar = zeros(size(data_6000),'like',data_6000);
beta = zeros(size(data_6000),'like',data_6000);

% skip zeros in volumes
disp('Starting axon diameter calculation');
tic
parfor i = 1:numel(data_6000)
    data_vx = [data_6000(i), data_30000(i)];
    g_vx = [g_corr_1(i), g_corr_2(i)];
    if all(data_vx)
        [ar(i), beta(i)] = getAxonRadius(delta,Delta,g_vx,data_vx,'Neumann');
    end
end
toc

savename = fullfile(path,'AxonRadiusMap');
nifti_info.Datatype = class(ar);
niftiwrite(ar, savename, nifti_info);

end