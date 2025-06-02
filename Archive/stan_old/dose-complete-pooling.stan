//
// Complete pooling dose model
// 2023-02-17
//

data{
    int<lower=0> N;
    array[N] int id;
    vector[N] y;
    vector[N] x;
    vector[N] dose;
}
parameters{
    array[2] real a;
    array[2] real b;
    real<lower=0> s;
}
model{
    vector[N] mu;
    s ~ exponential( 1 );
    b ~ normal( 0 , 5 );
    a ~ normal( 5 , 5 );
    for ( i in 1:N ) {
        mu[i] = a[dose[i]] + b[dose[i]] * x[i];
    }
    y ~ normal( mu , s );
}

