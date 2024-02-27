function GradientVal ...
    = Gradient_MaxLikelihood_GEVD( ...
    evtParam, ...
    sample ....
    )

GradientVal = zeros( 1, 3);

location = evtParam(1);
scale = evtParam(2);
shape = evtParam(3);

test_term = ( 1 + shape*(sample-location)/scale )^(-1/shape);
threshold = 1e-1;
if (test_term<threshold)
    exponent_term = exp( -(sample-location)/scale );
    common_term = (exponent_term - 1)/scale;
    
    GradientVal(1) = common_term;
    
    GradientVal(2) = ( common_term*(sample-location) + 1 )/scale;
    
    GradientVal(3) = 0; % does not exist
else
    exponent_term = ( 1 + shape*(sample-location)/scale )^(-1/shape);
    common_term = (exponent_term - shape - 1);
    
    GradientVal(1) = common_term * exponent_term^shape / scale;
    
    GradientVal(2) = (sample-location) * common_term * exponent_term^shape /...
        scale^2 + 1/scale;
    
    GradientVal(3) = common_term * ( exponent_term^shape - 1 ) / shape^2 ...
        + (common_term + shape) * log(exponent_term) / shape; %???????????
end



end