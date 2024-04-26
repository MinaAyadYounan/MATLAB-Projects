function output = blur(img,w)
mat = double(img);
k = 2*w+1;
[row,col] = size(img);
A = [];
result = [];
for i = 1:row
    for ii = 1:col
             a = i - w;
             c = ii - w;
             Rend = a+k -1;
             Cend = c+k -1 ;
               if a < 1
                   a = 1;
               end          
               if Rend > row
                   Rend = row;
               end
                 if c <1
                   c= 1;
               end          
               if Cend > col
                   Cend = col;
               end
                  A =mat([a:Rend],[c:Cend]);
                  result(i,ii) = mean(A(:));                  
    end
end
output = uint8(result);
end