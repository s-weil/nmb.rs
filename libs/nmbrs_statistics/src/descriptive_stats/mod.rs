mod array_stats;
mod sorted_array_stats;

pub use array_stats::{
    covariance, mean, sum, variance, Covariance, Mean, Sum, Variance, VarianceBias,
};

/* TODOs:
- splt into descriptive and inferential stats and ordered and unordered stats
- math formulas
- combined stats (returning mean and std) for perf reasons. e.g. provide a class which then sorts,
    holds the sorted array, the mean, etc, etc and provides the stats
- furhter metrics */
