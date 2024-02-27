function Vals = Test_MLE( scale, shape, samples)

Vals = zeros(1, 2);
n = length(samples);

Vals(1) = abs( sum(( samples / scale )./( 1 + scale * samples / shape )) - n /( 1 + shape ) );

Vals(2) = abs( sum( log( 1 + shape * samples / scale )) - n * shape );

end