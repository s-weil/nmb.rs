use std::{
    ops::{Add, Mul},
    usize,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Vector<const D: usize> {
    v: [f64; D],
}

impl<const D: usize> Vector<D> {
    pub fn new(v: [f64; D]) -> Self {
        Self { v }
    }
}

impl<const D: usize> From<[f64; D]> for Vector<D> {
    fn from(v: [f64; D]) -> Self {
        Self { v }
    }
}

type Scalar = f64;

impl<const D: usize> Add for Vector<D> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut v = rhs.v.clone();
        for i in 0..D {
            v[i] += self.v[i];
        }
        Self { v }
    }
}

impl<const D: usize> Mul<Scalar> for Vector<D> {
    type Output = Self;

    fn mul(self, rhs: Scalar) -> Self::Output {
        let mut v = self.v.clone();
        for i in 0..D {
            v[i] *= rhs * self.v[i];
        }
        Self { v }
    }
}

/// Convenicence syntax.
/// Write `V![3; 1.1, 2.2, 3.3]` for the $3$-dimensional vector `[1.1, 2.2, 3.3]`.
macro_rules! V {
    ( $ d : expr; $ ( $ x : expr), +  ) => {
        Vector::<$d>::new([ $ ( $ x ) , + ])
    };
}

#[cfg(test)]
mod tests {
    use super::Vector;

    #[test]
    fn add() {
        assert_eq!(
            Vector::<2>::new([1.0, 1.0]) + Vector::<2>::new([2.0, 2.0]),
            Vector::<2>::new([3.0, 3.0])
        );

        assert_eq!(
            Vector::<2>::new([2.0, 2.0]) + [1.0, 1.0].into(),
            [3.0, 3.0].into()
        );

        assert_eq!(
            V![2; 1.0, 1.0] + V![2; 2.0, 2.0],
            Vector::<2>::new([3.0, 3.0])
        );

        assert_eq!(V![2; 1.0, 1.0] + V![2; 2.0, 2.0], V![2; 3.0, 3.0]);
    }
}
