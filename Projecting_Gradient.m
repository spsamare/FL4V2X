function projectedParam....
    = Projecting_Gradient( ...
    localParam, ...
    this_sample ...
    )

a = localParam(1);
b = localParam(2);
c = localParam(3);
q = this_sample;

x = 0;
y = 0;
z = 0; znext = 1;

t = eps; %offset
opts1=  optimset('display','off');

while abs(z-znext) > 1e-3
    %disp(abs(z-znext));
    znext = z;
    %% fix z, find x,y
    H = 2*eye(2);
    f = -2*[a; b];
    lb = [0; 0.1];
    Acons = [z t-1];
    bcons = z*q;
    
    xy = quadprog(H,f,Acons,bcons,[],[],lb,[],[],opts1);
    x = xy(1);
    y = xy(2);
    
    %% fix x,y, find z
    z0 = (t-1)*y/(q-x);
    if q-x == 0
        z = c;
    elseif q-x > 0 % z >= z0
        z = max(c, z0);
    else % z <= z0
        z = min(c, z0);
    end    
end

projectedParam = [x y z];
% disp('done');
end