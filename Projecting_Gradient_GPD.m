function projectedParam....
    = Projecting_Gradient_GPD( ...
    localParam, ...
    this_sample ...
    )

b = localParam(1);
c = localParam(2);
q = this_sample;

LB_scale = 10;
UB_scale = 70;

UB_shape = 0.1;

opts1=  optimset('display','off');
e = 0;%1e-3;

if (b + c*q < 0)|| (b<=LB_scale) || (c>=UB_shape)
    H = 2*eye(2);
    f = -2*[b; c];
    lb = []; %moved to a constraint
    ub = []; %moved to a constraint
    Acons = [ ...
        1 0; ...
        -1 0; ...
        0 1; ....
        -1*ones(size(q)) -q ...
        ];
    bcons = [...
        UB_scale-e; ...
        -LB_scale-e; ...
        UB_shape-e; ...
        -e*ones(size(q)) ...
        ];
    
    yz = quadprog(H,f,Acons,bcons,[],[],lb,ub,[],opts1);
    
    yz_rounded = round(yz*10000)/10000;
    
    projectedParam = [yz_rounded(1) yz_rounded(2)];
    % disp('Projected');
    
else
    projectedParam = [b c];
end
% disp('done');
end