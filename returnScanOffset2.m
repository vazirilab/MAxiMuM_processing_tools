function correction = returnScanOffset2(Iin,dim)

if numel(size(Iin)) == 3
    Iin = mean(Iin,3);
elseif numel(size(Iin)) == 4
    Iin = mean(mean(Iin,4),3);
end

n = 8;

switch dim
    case 1
        Iv1 = Iin(1:2:end,:);
        Iv2 = Iin(2:2:end,:);
        
        Iv1 = Iv1(1:min([size(Iv1,1) size(Iv2,1)]),:); 
        Iv2 = Iv2(1:min([size(Iv1,1) size(Iv2,1)]),:); 
        
        buffers = zeros(size(Iv1,1),n);

        Iv1 = cat(2,buffers,Iv1,buffers);
        Iv2 = cat(2,buffers,Iv2,buffers);

        Iv1 = reshape(Iv1',[],1);
        Iv2 = reshape(Iv2',[],1);
        
    case 2
        Iv1 = Iin(:,1:2:end);
        Iv2 = Iin(:,2:2:end);

        Iv1 = Iv1(:,1:min([size(Iv1,2) size(Iv2,2)])); 
        Iv2 = Iv2(:,1:min([size(Iv1,2) size(Iv2,2)]),:); 
        
        buffers = zeros(n,size(Iv1,2));

        Iv1 = cat(1,buffers,Iv1,buffers);
        Iv2 = cat(1,buffers,Iv2,buffers);
        
        Iv1 = reshape(Iv1,[],1);
        Iv2 = reshape(Iv2,[],1);
end



Iv1 = Iv1-mean(Iv1); Iv1(Iv1<0) = 0;
Iv2 = Iv2-mean(Iv2); Iv2(Iv2<0) = 0;

[r,lag] = xcorr(Iv1,Iv2,n,'unbiased');

% figure; plot(lag,r)

[~,ind] = max(r);
correction = lag(ind);
