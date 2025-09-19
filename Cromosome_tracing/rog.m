function yy = rog(Trace)
Mean = mean(Trace,1); %finds the average coordinate in x, y, z
Sum = 0;
N = 0;
for i = 1:size(Trace,1) %1 to number of coordinates in Trace
    Sum = Sum + sum((Trace(i,:)-Mean).^2);
    N = N+1;
end
yy = (Sum/N)^0.5;