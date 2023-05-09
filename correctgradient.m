function [corr_fac] = correctgradient(bval,bvec,gd)

bval1 = 6000;
bval2 = 30450;

dims=size(gd);
dims=dims(1:3);
gd=reshape(gd,[dims,3,3]);
gd=permute(gd,[5,4,1,2,3]);

% allow some tolerance
bidxs{1}=(bval<bval1+100 & bval>bval1-100);
bidxs{2}=(bval<bval2+100 & bval>bval2-100);

corr_fac=cell(2,1);
for bidx=1:length(bidxs)
    corr_fac{bidx}=zeros(dims);
    for n=1:numel(corr_fac{bidx})
        [i,j,k]=ind2sub(dims,n);
        L=eye(3)+squeeze(gd(:,:,i,j,k));

        bvecL=bvec(bidxs{bidx},:)*L;

        corr_fac{bidx}(n)=sqrt(mean(sum(bvecL.^2,2)));
    end
end