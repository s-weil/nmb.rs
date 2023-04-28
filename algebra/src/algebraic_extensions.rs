use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait AddIdentity: Sized {
    fn zero() -> Self;
}

/// Mimic features of a [(mathematical) semigroup, or monoid](https://en.wikipedia.org/wiki/Monoid).
pub trait NumericSemiGroup: Sized + AddIdentity + Add<Output = Self> + PartialEq {
    /// Assicoativity can only be asserted and must be assumed.
    fn is_associative(a: Self, b: Self, c: Self) -> bool
    where
        Self: Copy,
    {
        (a + b) + c == a + (b + c)
    }
}

impl NumericSemiGroup for usize {}
impl NumericSemiGroup for i8 {}
impl NumericSemiGroup for i16 {}
impl NumericSemiGroup for i32 {}
impl NumericSemiGroup for i64 {}
impl NumericSemiGroup for f32 {}
impl NumericSemiGroup for f64 {}

/// Mimic features of a [(mathematical) group](https://en.wikipedia.org/wiki/Group_(mathematics)#Definition).
///
/// Technically the bound `Sub` is not required, but it is implied (from `Add` + `Neg`) and added for convenience.
pub trait NumericGroup: NumericSemiGroup + Neg<Output = Self> + Sub<Output = Self> {}

impl<T> NumericGroup for T where T: NumericSemiGroup + Neg<Output = Self> + Sub<Output = Self> {}
// impl NumericGroup for i8 {}
// impl NumericGroup for i16 {}
// impl NumericGroup for i32 {}
// impl NumericGroup for i64 {}
// impl NumericGroup for f32 {}
// impl NumericGroup for f64 {}

pub trait MulIdentity: Sized {
    fn one() -> Self;
}

/// Mimic features of a [(mathematical) ring](https://en.wikipedia.org/wiki/Ring_(mathematics)#Definition),
/// without ring axioms (Abelian group, Associativity, Distributivity)
pub trait NumericRing: NumericGroup + MulIdentity + Mul<Output = Self> {
    fn is_commutative(a: Self, b: Self) -> bool
    where
        Self: Copy,
    {
        a * b == b * a
    }

    fn is_mul_associative(a: Self, b: Self, c: Self) -> bool
    where
        Self: Copy,
    {
        (a * b) * c == a * (b * c)
    }

    fn is_distributive(a: Self, b: Self, c: Self) -> bool
    where
        Self: Copy,
    {
        a * (b + c) == a * b + a * c
    }
}

impl<T> NumericRing for T where T: NumericGroup + MulIdentity + Mul<Output = Self> {}

/// Mimic features of a [(mathematical) field](https://en.wikipedia.org/wiki/Field_(mathematics)).
pub trait NumericField: NumericRing + Div<Output = Self> {
    fn inverse(a: Self) -> Self
    where
        Self: Copy,
    {
        if a == Self::zero() {
            panic!("Cannot divide by zero")
        }
        Self::one() / a
    }
}

impl<T> NumericField for T where T: NumericRing + Div<Output = Self> {}

#[macro_export]
macro_rules! impl_add_identity {
    ($impl_type:ty) => {
        impl AddIdentity for $impl_type {
            fn zero() -> Self {
                0 as $impl_type
            }
        }
    };
}

#[macro_export]
macro_rules! impl_mul_identity {
    ($impl_type:ty) => {
        impl MulIdentity for $impl_type {
            fn one() -> Self {
                1 as $impl_type
            }
        }
    };
}

// implement MulIdentity
impl_mul_identity! { usize }
impl_mul_identity! { i8 }
impl_mul_identity! { i16 }
impl_mul_identity! { i32 }
impl_mul_identity! { i64 }
impl_mul_identity! { f32 }
impl_mul_identity! { f64 }

// implement AddIdentity
impl_add_identity! { usize}
impl_add_identity! { i8 }
impl_add_identity! { i16 }
impl_add_identity! { i32 }
impl_add_identity! { i64 }
impl_add_identity! { f32 }
impl_add_identity! { f64 }

// impl NumericRing for f32 {}
// impl NumericRing for f64 {}
// impl NumericField for f32 {}
// impl NumericField for f64 {}
