function fitSH(data_file, mask_file, noise_file, bval, bvec, fname_out)

% load data
nii = niftiread(mask_file);
mask = nii>0;
nii = niftiread(data_file);
y = vectorize(single(nii), mask);
nii = niftiread(noise_file);
sigma = vectorize(single(nii), mask);
bvec = importdata(bvec).'; 
bval = importdata(bval);

[~, nvox] = size(y);

% MLE fitting per shell, extract zeroth order SH coefficient only and
% normalize by b0

bvalues = [0, 6000, 30450];
N = 0;
disp('Starting spherical harmonic fit');
for i = 1:numel(bvalues)
    list = bval == bvalues(i);
    dirs = bvec(list, :);
    
    if bvalues(i)==0
        b0 = fitSphericalMean_mle(dirs, y(list, :), sigma, 0);
    else
        N = N+1;
        K0(N,:) = fitSphericalMean_mle(dirs, y(list, :), sigma, 6);
        b(N, 1) = bvalues(i);
    end
end

K0_1 = vectorize(K0(1,:),mask);
K0_2 = vectorize(K0(2,:),mask);
K0_1(isnan(K0_1))=0;
K0_2(isnan(K0_2))=0;

nifti_info = niftiinfo(mask_file);
nifti_info.Datatype = class(K0_1);

niftiwrite(K0_1, strcat(fname_out, '_sh_b6000_powderavg'), nifti_info);
niftiwrite(K0_2, strcat(fname_out, '_sh_b30000_powderavg'), nifti_info);

end
