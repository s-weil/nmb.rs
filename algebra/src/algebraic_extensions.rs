use std::ops::{Add, Div, Mul, Neg};

/// Mimic some features of a ring, without ring axioms (Abelian group, Associativity, Distributivity)
pub trait NumericRing:
    Sized + AddIdentity + Add<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
}

pub trait MulIdentity: Sized {
    fn one() -> Self;
}

pub trait NumericField: NumericRing + MulIdentity + Div<Output = Self> //Option<Self>>
{
}

pub trait AddIdentity: Sized {
    fn zero() -> Self;
}

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
impl_mul_identity! { i8 }
impl_mul_identity! { i16 }
impl_mul_identity! { i32 }
impl_mul_identity! { i64 }
impl_mul_identity! { f32 }
impl_mul_identity! { f64 }

// implement AddIdentity
impl_add_identity! { i8 }
impl_add_identity! { i16 }
impl_add_identity! { i32 }
impl_add_identity! { i64 }
impl_add_identity! { f32 }
impl_add_identity! { f64 }

impl NumericRing for f32 {}
impl NumericRing for f64 {}
impl NumericField for f32 {}
impl NumericField for f64 {}
