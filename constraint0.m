function constraint=constraint0(k)
constraint = (sum(k(4:5)) >=0) & (k(1)*k(4)+k(2)*k(5) +  k(3)*(1-k(4)-k(5)) <= 0);
