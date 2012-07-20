C1(t,b,a) = (b**4) - (a**4) * (sin(2*t)**2);
C(t,b,a)  = (a**2) * cos(2*t) + sqrt(C1(t,b,a));
X(t,b,a)  = cos(t) * sqrt(C(t,b,a));
Y(t,b,a)  = sin(t) * sqrt(C(t,b,a));
