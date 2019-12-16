//! Identity
use super::operator::{Operator, Addition, Multiplication};

/// A type that is equipped with identity.
pub trait Identity<O: Operator>
{
    /// The identity element.
    fn id() -> Self;
}

macro_rules! impl_identity
{
    ($t:ty; $v:expr; $($s:ty),*) =>
    {
        $(
        impl Identity<$t> for $s
        {
            fn id() -> Self
            {
                return $v;
            }
        }
        )*
    };
}

impl_identity!(Addition; 0; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128);
impl_identity!(Addition; 0.; f32, f64);
impl_identity!(Multiplication; 1; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128);
impl_identity!(Multiplication; 1.; f32, f64);