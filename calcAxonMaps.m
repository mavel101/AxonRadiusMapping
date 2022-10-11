function calcAxonMaps(input_6000, input_30000)

data_6000 = niftiread(input_6000);
data_6000 = double(data_6000(:,:,:,1));
data_6000(isnan(data_6000)) = 0;
data_30000 = niftiread(input_30000);
data_30000 = double(data_30000(:,:,:,1));
data_30000(isnan(data_30000)) = 0;

nifti_info = niftiinfo(input_6000);
nifti_info.ImageSize = nifti_info.ImageSize(1:3);
nifti_info.PixelDimensions = nifti_info.PixelDimensions(1:3);

delta = [15, 15];
Delta = [29.25, 29.25];
g = [sqrt(6000/30450)*280, 280];

ar = zeros(size(data_6000),'like',data_6000);
beta = zeros(size(data_6000),'like',data_6000);

% skip zeros in volumes
ix = find(data_6000);
for i = 1:size(ix,1)
    ix_vx = ix(i);
    data_vx = [data_6000(ix_vx), data_30000(ix_vx)];
    [ar(ix_vx), beta(ix_vx)] = getAxonRadius(delta,Delta,g,data_vx,'Neumann');
end

[path, name] = fileparts(input_6000);
name = extractBefore(name, ".");
savename = strcat(path,'/',name,'_axonMap');

niftiwrite(ar, savename, nifti_info);

end