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

% correct gradient for slightly incorrect b-value calculation in getAxonRadius
gyroMagnRatio =  267.513*10^(-6);
q = g*gyroMagnRatio;
b = (q.*delta).^2.*(Delta - delta/3); % from getAxonRadius
g(1) = g(1) * sqrt(bval1/b(1));
g(2) = g(2) * sqrt(bval2/b(2));

% correct gradient with factors from nonlinearity correction
bval = importdata(in_bval);
bvec = importdata(in_bvec);
grad_dev = niftiread(in_grad_dev);
corr_fac = correctgradient(bval, bvec, grad_dev);
g_corr =  {corr_fac{1} * g(1), corr_fac{2} * g(2)};

ar = zeros(size(data_6000),'like',data_6000);
beta = zeros(size(data_6000),'like',data_6000);

% skip zeros in volumes
ix = find(data_6000);
disp('Starting axon diameter calculation');
for i = 1:size(ix,1)
    ix_vx = ix(i);
    data_vx = [data_6000(ix_vx), data_30000(ix_vx)];
    g_vx = [g_corr{1}(ix_vx), g_corr{2}(ix_vx)];
    [ar(ix_vx), beta(ix_vx)] = getAxonRadius(delta,Delta,g_vx,data_vx,'Neumann');
    if mod(i,floor(size(ix,1)/10))<1e-2
       fprintf('Progress: %d %%\n', floor(i/size(ix,1)*100));
    end
end

[path, name] = fileparts(in_6000);
name = extractBefore(name, ".");
savename = strcat(path,'/',name,'_axonMap');

nifti_info.Datatype = class(ar);
niftiwrite(ar, savename, nifti_info);

end