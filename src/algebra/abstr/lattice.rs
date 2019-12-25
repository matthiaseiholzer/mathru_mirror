


pub trait Lattice : PartialOrd
{

}

macro_rules! impl_lattice
{
    ($($t:ty),*) =>
    {
        $(
        impl Lattice for $t
        {

        }
        )*
    };
}

impl_lattice!(u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);