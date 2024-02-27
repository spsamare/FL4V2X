function GradientVal ...
    = Gradient_MaxLikelihood_GPD( ...
    evtParam, ...
    sample ....
    )

GradientVal = zeros( 1, 2);

scale = evtParam(1);
shape = evtParam(2);

exponent_term = ( 1 + shape*sample/scale )^(-1/shape);
if (exponent_term < 1e-3)
    exponent_term = 0;
end

GradientVal(1) = ( (1 + 1/shape) * exponent_term^shape - 1/shape ) ...
    / scale;

GradientVal(2) = (shape + 1)*(1 - exponent_term^shape) / shape^2 + ...
    log(exponent_term) / shape;

if (sum(isinf(GradientVal)+isnan(GradientVal))>0) ...
        || (abs(shape)<1e-2)
    
    GradientVal(1) = (scale - sample) / scale^2;
    
    GradientVal(2) = -shape; %undefined
    
end



end