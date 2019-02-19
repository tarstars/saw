function [c] = makeTetragonalMaterialTensor(c11, c12, c13, c33, c66, c44)
c2d=[[c11 c12 c13   0   0   0];
     [c12 c11 c13   0   0   0];
     [c13 c13 c33   0   0   0];
     [  0   0   0 c44   0   0];
     [  0   0   0   0 c44   0];
     [  0   0   0   0   0 c66]];
 
c=zeros(3,3,3,3);
 
for p = 1:3
    for q = 1:3
        for r = 1:3
            for s = 1:3
                c(p, q, r, s) = c2d(fromPairToIndex(p,q), fromPairToIndex(r,s));
            end
        end
    end
end
 
                 
 
 
    