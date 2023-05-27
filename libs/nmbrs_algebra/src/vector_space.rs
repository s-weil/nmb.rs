use crate::{
    algebraic_extensions::{AddIdentity, Inverse, NumericField},
    NumericGroup, NumericRing,
};
use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub},
    usize,
};

// TODO: improve on trait bounds below

/// Representation of a [`Vector Space`](https://en.wikipedia.org/wiki/Vector_space).
///
pub trait VectorSpace:
    AddIdentity + Add<Output = Self> + Inverse + Mul<Self::Field, Output = Self> + Sized
{
    type Field: NumericField;

    // 0 for add, 1 for mul
    // associativity
    // distributivity for add and mul
}

impl<F: NumericField + Copy + MulAssign + AddAssign, const D: usize> VectorSpace for Vector<D, F> {
    type Field = F;
}

impl<F> VectorSpace for F
where
    F: NumericField,
{
    type Field = F;
}

pub trait VectorSpaceF32: VectorSpace<Field = f32> {}
pub trait VectorSpaceF64: VectorSpace<Field = f64> {}

// impl<F, V> Mul<V> for F
// where
//     F: NumericField,
//     V: VectorSpace<F>,
// {
//     fn mul(self, rhs: V) -> Self::Output {
//         rhs * self
//     }
// }

// TODO: impl for ndarray and nalgebra vectors

#[derive(Debug, Clone, PartialEq)]
pub struct Vector<const D: usize, F> {
    v: [F; D],
}

impl<const D: usize, F> Vector<D, F> {
    pub fn new(v: [F; D]) -> Self {
        Self { v }
    }
}

impl<const D: usize, F> From<[F; D]> for Vector<D, F> {
    fn from(v: [F; D]) -> Self {
        Self { v }
    }
}

impl<const D: usize, F> Copy for Vector<D, F> where F: Copy {}
// impl<const D: usize, F> Clone for Vector<D, F> where F: Clone {}

impl<const D: usize, F: NumericGroup + Copy> AddIdentity for Vector<D, F> {
    fn zero() -> Self {
        [F::zero(); D].into()
    }
}

impl<const D: usize, F: NumericGroup + AddAssign + Copy> Add for Vector<D, F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut v = rhs.v.clone();
        for i in 0..D {
            v[i] += self.v[i];
        }
        Self { v }
    }
}

impl<const D: usize, F: NumericGroup + Copy> Neg for Vector<D, F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut v = self.v.clone();
        for i in 0..D {
            v[i] = -self.v[i];
        }
        Self { v }
    }
}

impl<const D: usize, F: NumericGroup + Copy> Sub for Vector<D, F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut v = rhs.v.clone();
        for i in 0..D {
            v[i] = self.v[i] - rhs.v[i];
        }
        Self { v }
    }
}

impl<const D: usize, F: NumericRing + MulAssign + Copy> Mul<F> for Vector<D, F> {
    type Output = Self;

    fn mul(self, rhs: F) -> Self::Output {
        let mut v = self.v.clone();
        for i in 0..D {
            v[i] *= rhs;
        }
        Self { v }
    }
}

/// Convenicence syntax.
///
/// Write `V![3; 1.1, 2.2, 3.3]` for the $3$-dimensional vector `[1.1, 2.2, 3.3]`.
macro_rules! V {
    ( $d:expr; $( $x:expr ), +  ) => {
        Vector::<$d, f64>::new([ $ ( $x ) , + ])
    };
}

#[cfg(test)]
mod tests {
    use super::Vector;

    #[test]
    fn add() {
        assert_eq!(
            Vector::<2, f64>::new([1.0, 1.0]) + Vector::<2, f64>::new([2.0, 2.0]),
            Vector::<2, f64>::new([3.0, 3.0])
        );

        assert_eq!(
            Vector::<2, f64>::new([2.0, 2.0]) + [1.0, 1.0].into(),
            [3.0, 3.0].into()
        );

        assert_eq!(
            V![2; 1.0, 1.0] + V![2; 2.0, 2.0],
            Vector::<2, f64>::new([3.0, 3.0])
        );

        assert_eq!(V![2; 1.0, 1.0] + V![2; 2.0, 2.0], V![2; 3.0, 3.0]);
    }

    #[test]
    fn scalar() {
        assert_eq!(
            Vector::<2, f64>::new([2.0, 3.0]) * 2.0,
            Vector::<2, f64>::new([4.0, 6.0])
        );

        assert_eq!(Vector::<2, f64>::new([2.0, 3.0]) * 2.0, [4.0, 6.0].into());

        assert_eq!(V![2; 2.0, 3.0] * 2.0, V![2; 4.0, 6.0]);
    }
}
