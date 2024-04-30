function fgp = EL_gravityvector(params,Lo,Ln,c)

% Function to calculate the gravity force vector for each element.

p = params.p;
L = Lo;
Lca = Ln * c;

fgp = [0, -L*p/2, -(L*p/12)*Lca, 0, -L*p/2, (L*p/12)*Lca]';

end
