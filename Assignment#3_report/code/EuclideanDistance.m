function Distance = EuclideanDistance(x1, x2)

sum = 0;

for counterDimension = 1:size(x1, 2),
    sum = sum + (x1(1,counterDimension) - x2(1,counterDimension))^2;
end

Distance = sqrt(sum);
       
end
