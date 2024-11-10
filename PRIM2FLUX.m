function F=PRIM2FLUX(Q,c)

% UNPACK VARIABLES
w = Q(:,1);
v = Q(:,2);

% BUILD LINEAR WAVE EQUATION FLUX VECTOR
f1 = -v;
f2 = -c^2*w;

F=[f1 f2];