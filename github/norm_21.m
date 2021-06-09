function result = norm_21(data)

B = data.*data;
c = sum(B,2);
D = sqrt(c);
result = sum(D);
end