use crate::{
    algebra::abstr::{
        cast,
        cast::{FromPrimitive, NumCast, ToPrimitive},
        AbelianGroup, AbelianGroupAdd, AbelianGroupMul, Addition, CommutativeRing,
        Complex as ComplexT, Field, Group, GroupAdd, GroupMul, Identity, Loop, Magma, MagmaAdd,
        MagmaMul, Monoid, MonoidAdd, MonoidMul, Multiplication, One, Quasigroup, Real, Ring,
        Scalar, Semigroup, SemigroupAdd, SemigroupMul, Sign, Zero,
    },
    elementary::{Exponential, Hyperbolic, Power, Trigonometry},
    algebra::abstr::AbsDiffEq,
};

use std::{
    cmp::Ordering,
    fmt,
    fmt::Display,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[cfg(feature = "blaslapack")]
use crate::algebra::abstr::{Blas, Lapack};

#[doc(hidden)]
#[macro_export]
macro_rules! Complex {
    ($r:expr, $i:expr) => {
        Complex::new($r, $i);
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! Complex32 {
    ($r:expr, $i:expr) => {
        (Complex::new($r, $i));
    };
}
#[doc(hidden)]
#[macro_export]
macro_rules! Complex64 {
    ($r:expr, $i:expr) => {
        (Complex::new($r, $i));
    };
}

pub type Complex32 = Complex<f32>;
pub type Complex64 = Complex<f64>;

/// Complex number in cartesian form
#[derive(Debug, Clone, Copy)]
pub struct Complex<T>
    where T: Real
{
    /// Real portion of the complex number
    re: T,
    /// Imaginary  portion of the complex number
    im: T,
}

impl<T> Complex<T> where T: Real
{
    /// Create a new complex
    pub fn new(re: T, im: T) -> Complex<T>
    {
        Complex { re, im }
    }
}

/// Returns -self = -Re(self) - i Im(self)
impl<T> Neg for Complex<T> where T: Real /* Neg<Output = T> */
{
    type Output = Complex<T>;

    fn neg(self: Self) -> Complex<T>
    {
        Complex { re: -self.re,
                  im: -self.im }
    }
}

/// Complex numbers can not be compared
///
/// # Panics
///
/// # FIXME
impl<T> PartialOrd for Complex<T> where T: Real /* PartialOrd */
{
    fn partial_cmp(self: &Self, other: &Self) -> Option<Ordering>
    {
        self.re.partial_cmp(&other.re)
    }
}

impl<T> ComplexT for Complex<T> where T: Real
{
    /// Returns the complex conjuagte
    /// conj(self) = Re(self) - i Im(self)
    fn conj(self: Self) -> Complex<T>
    {
        Complex { re: self.re,
                  im: -self.im }
    }

    /// Returns the argument of the complex number
    fn arg(self: Self) -> Self
    {
        Complex { re: self.im.arctan2(self.re),
                  im: T::zero() }
    }
}

//impl<T> Real for Complex<T>
//	where T: Real + Copy + Add<T, Output = T> + Mul<T, Output = T>
//{
//	/// Returns the smallest integer greater than or equal to a number.
//	fn ceil(self: &Self) -> Self
//	{
//		unimplemented!();
//	}
//
//	/// Returns the largest integer less than or equal to a number.
//	fn floor(self: &Self) -> Self
//	{
//		unimplemented!();
//	}
//}

#[cfg(feature = "blaslapack")]
impl<T> Lapack for Complex<T> where T: Real
{
    fn xgehrd(_n: i32,
              _ilo: i32,
              _ihi: i32,
              _a: &mut [Self],
              _lda: i32,
              _tau: &mut [Self],
              _work: &mut [Self],
              _lwork: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgehrd_work_size(_n: i32,
                        _ilo: i32,
                        _ihi: i32,
                        _a: &mut [Self],
                        _lda: i32,
                        _tau: &mut [Self],
                        _info: &mut i32)
                        -> i32
    {
        unimplemented!();
    }

    fn xorghr(_n: i32,
              _ilo: i32,
              _ihi: i32,
              _a: &mut [Self],
              _lda: i32,
              _tau: &[Self],
              _work: &mut [Self],
              _lwork: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xorghr_work_size(_n: i32,
                        _ilo: i32,
                        _ihi: i32,
                        _a: &mut [Self],
                        _lda: i32,
                        _tau: &[Self],
                        _info: &mut i32)
                        -> i32
    {
        unimplemented!();
    }

    fn xgeev(_jobvl: u8,
             _jobvr: u8,
             _n: i32,
             _a: &mut [Self],
             _lda: i32,
             _wr: &mut [Self],
             _wi: &mut [Self],
             _vl: &mut [Self],
             _ldvl: i32,
             _vr: &mut [Self],
             _ldvr: i32,
             _work: &mut [Self],
             _lwork: i32,
             _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgeev_work_size(_jobvl: u8,
                       _jobvr: u8,
                       _n: i32,
                       _a: &mut [Self],
                       _lda: i32,
                       _wr: &mut [Self],
                       _wi: &mut [Self],
                       _vl: &mut [Self],
                       _ldvl: i32,
                       _vr: &mut [Self],
                       _ldvr: i32,
                       _info: &mut i32)
                       -> i32
    {
        unimplemented!();
    }

    fn xgetrf(_m: i32, _n: i32, _a: &mut [Self], _lda: i32, _ipiv: &mut [i32], _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgeqrf(_m: i32,
              _n: i32,
              _a: &mut [Self],
              _lda: i32,
              _tau: &mut [Self],
              _work: &mut [Self],
              _lwork: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgeqrf_work_size(_m: i32,
                        _n: i32,
                        _a: &mut [Self],
                        _lda: i32,
                        _tau: &mut [Self],
                        _info: &mut i32)
                        -> i32
    {
        unimplemented!();
    }

    fn xorgqr(_m: i32,
              _n: i32,
              _k: i32,
              _a: &mut [Self],
              _lda: i32,
              _tau: &mut [Self],
              _work: &mut [Self],
              _lwork: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xorgqr_work_size(_m: i32,
                        _n: i32,
                        _k: i32,
                        _a: &mut [Self],
                        _lda: i32,
                        _tau: &mut [Self],
                        _info: &mut i32)
                        -> i32
    {
        unimplemented!();
    }

    fn xgetri(_n: i32,
              _a: &mut [Self],
              _lda: i32,
              _ipiv: &mut [i32],
              _work: &mut [Self],
              _lwork: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgetri_work_size(_n: i32,
                        _a: &mut [Self],
                        _lda: i32,
                        _ipiv: &mut [i32],
                        _info: &mut i32)
                        -> i32
    {
        unimplemented!();
    }

    fn xpotrf(_uplo: char, _n: i32, _a: &mut [Self], _lda: i32, _info: &mut i32)
    {
        unimplemented!();
    }

    fn xgetrs(_n: i32,
              _nrhs: i32,
              _a: &mut [Self],
              _lda: i32,
              _ipiv: &mut [i32],
              _b: &mut [Self],
              _ldb: i32,
              _info: &mut i32)
    {
        unimplemented!();
    }

    fn xggev(
        n: i32,
        a: &mut [Self],
        lda: i32,
        b: &mut [Self],
        ldb: i32,
        alphar: &mut [Self],
        alphai: &mut [Self],
        beta: &mut [Self],
        vl: &mut [Self],
        ldvl: i32,
        vr: &mut [Self],
        ldvr: i32,
        work: &mut [Self],
        lwork: i32,
        info: &mut i32)
    {
        unimplemented!();
    }

    fn xggev_work_size(n: i32,
                       a: &mut [Self],
                       lda: i32,
                       b: &mut [Self],
                       ldb: i32,
                       alphar: &mut [Self],
                       alphai: &mut [Self],
                       beta: &mut [Self],
                       vl: &mut [Self],
                       ldvl: i32,
                       vr: &mut [Self],
                       ldvr: i32,
                       info: &mut i32 ) -> i32
    {
        unimplemented!();
    }
}

#[cfg(feature = "blaslapack")]
impl<T> Blas for Complex<T> where T: Real
{
    fn xgemm(_transa: u8,
             _transb: u8,
             _m: i32,
             _n: i32,
             _k: i32,
             _alpha: Self,
             _a: &[Self],
             _lda: i32,
             _b: &[Self],
             _ldb: i32,
             _beta: Self,
             _c: &mut [Self],
             _ldc: i32)
    {
        unimplemented!();
    }

    fn xtrsm(_side: char,
             _uplo: char,
             _transa: char,
             _diag: char,
             _m: i32,
             _n: i32,
             _alpha: Self,
             _a: &[Self],
             _lda: i32,
             _b: &mut [Self],
             _ldb: i32)
    {
        unimplemented!();
    }

    fn xscal(_n: i32, _a: Self, _x: &mut [Self], _inc: i32)
    {
       unimplemented!();
    }

    fn xaxpy(_n: i32, _a: Self, _x: &[Self], _incx: i32, _y: &mut [Self], _incy: i32)
    {
        unimplemented!();
    }
}

//impl<T> Field for Complex<T>
//	where T: Real
//{
//	fn epsilon() -> Complex<T>
//	{
//		Complex
//		{
//			re: T::epsilon(),
//			im: T::epsilon()
//		}
//	}
//}

impl<T> Sign for Complex<T> where T: Real
{
    fn sign(self: &Self) -> Self
    {
        unimplemented!()
    }

    /// Absolute value of the complex number
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     algebra::abstr::{cast::ToPrimitive, Sign},
    ///     num::Complex,
    /// };
    ///
    /// let a: Complex<f64> = Complex::new(1.0, 2.0);
    /// let refer: f64 = (5.0_f64).sqrt();
    /// assert_eq!(refer, a.abs().to_f64());
    /// ```
    fn abs(self: &Self) -> Self
    {
        let root: T = T::from_f64(2.0);
        Complex { re: (self.re * self.re + self.im * self.im).root(root),
                  im: T::zero() }
    }

    fn is_positive(self: &Self) -> bool
    {
        unimplemented!();
    }

    fn is_negative(&self) -> bool
    {
        unimplemented!();
    }
}

//impl<T> Semiring for Complex<T>
//	where T: Field
//{
//
//}

impl<T> Scalar for Complex<T> where T: Real
{
}

/// Compares to complex numbers
impl<T> PartialEq for Complex<T> where T: Real
{
    fn eq<'a, 'b>(self: &'a Self, rhs: &Self) -> bool
    {
        if self.re == rhs.re && self.im == rhs.im
        {
            return true;
        }
        false
    }
}

impl<T> Display for Complex<T> where T: Real /* Display */
{
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result
    {
        write!(f, "{} + {}i", self.re, self.im)
    }
}

/// Returns 0 + i0
impl<T> Zero for Complex<T> where T: Real /* Zero */
{
    fn zero() -> Self
    {
        return Complex::new(T::zero(), T::zero());
    }
}

/// Adds two complex numbers
///
/// (a + ib) + (c + id) = (a + c) + i(b + d)
impl<T> Add for Complex<T> where T: Real /* Add<T, Output = T> */
{
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Complex<T>
    {
        Complex { re: self.re + rhs.re,
                  im: self.im + rhs.im }
    }
}

impl<'a, 'b, T> Add<&'b Complex<T>> for &'a Complex<T>
    where T: Real /* Add<T, Output = T> + Clone + Copy */
{
    type Output = Complex<T>;

    fn add(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re + rhs.re,
                  im: self.im + rhs.im }
    }
}

impl<T> AddAssign for Complex<T> where T: Real /* AddAssign */
{
    fn add_assign<'a>(self: &'a mut Self, other: Self)
    {
        self.re += other.re;
        self.im += other.im;
    }
}

/// Returns 1 + i0
impl<T> One for Complex<T> where T: Real
{
    fn one() -> Self
    {
        return Complex::new(T::one(), T::zero());
    }
}

/// Multiplies two complex numbers
///
/// (a + ib)(c + id) = (ac - bd) + i(bc + ad)
impl<T> Mul for Complex<T>
    where T: Real /* Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone + Copy */
{
    type Output = Complex<T>;

    fn mul(self: Self, other: Self) -> Self::Output
    {
        Complex { re: self.re * other.re - self.im * other.im,
                  im: self.im * other.re + self.re * other.im }
    }
}

impl<'a, 'b, T> Mul<&'b Complex<T>> for &'a Complex<T>
    where T: Real /* Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone + Copy */
{
    type Output = Complex<T>;

    fn mul(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re * rhs.re - self.im * rhs.im,
                  im: self.im * rhs.re + self.re * rhs.im }
    }
}

impl<T> MulAssign for Complex<T>
    where T: Real /* Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone + Copy */
{
    fn mul_assign(self: &mut Self, other: Self)
    {
        let temp: Complex<T> = *self;
        self.re = temp.re * other.re - temp.im * other.im;
        self.im = temp.im * other.re + temp.re * other.im;
    }
}

///Subtracts two complex numbers
///
/// (a + ib) - (c + id) = (a - c) + i(b - d)
impl<T> Sub for Complex<T> where T: Real /* Sub<T, Output = T> + Clone + Copy */
{
    type Output = Complex<T>;

    fn sub(self: Self, rhs: Self) -> Self::Output
    {
        &self - &rhs
    }
}

impl<'a, 'b, T> Sub<&'b Complex<T>> for &'a Complex<T>
    where T: Real /* Sub<T, Output = T> + Clone + Copy */
{
    type Output = Complex<T>;

    fn sub(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re - rhs.re,
                  im: self.im - rhs.im }
    }
}

impl<T> SubAssign for Complex<T> where T: Real /* SubAssign + Clone + Copy */
{
    fn sub_assign<'a>(self: &'a mut Self, other: Self)
    {
        self.re -= other.re;
        self.im -= other.im;
    }
}

/// Divides two complex numbers
impl<T> Div for Complex<T> where T: Real
{
    type Output = Complex<T>;

    fn div(self: Self, rhs: Self) -> Self::Output
    {
        &self / &rhs
    }
}

impl<'a, 'b, T> Div<&'b Complex<T>> for &'a Complex<T> where T: Real
{
    type Output = Complex<T>;

    fn div(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        let divisor: T = rhs.re * rhs.re + rhs.im * rhs.im;

        if divisor == T::zero()
        {
            panic!("rhs * rhs.conj() == 0")
        }
        let quot: Complex<T> = self * &Complex::new(rhs.re, -rhs.im);

        return Complex::new(quot.re / divisor, quot.im / divisor);
    }
}

impl<T> DivAssign for Complex<T> where T: Real
{
    fn div_assign<'a>(self: &'a mut Self, other: Self)
    {
        *self = *self / other;
    }
}

macro_rules! impl_to_primitive {
    ($ty:ty, $to:ident) => {
        fn $to(&self) -> $ty
        {
            return self.re.$to();
            //            if self.im == T::zero()
            //            {
            //            	self.re.$to()
            //            }
            //            else
            //            {
            //            	None
            //            }
        }
    };
}

/// Returns None if Complex part is non-zero
impl<T: ToPrimitive> ToPrimitive for Complex<T> where T: Real
{
    impl_to_primitive!(u8, to_u8);

    impl_to_primitive!(u16, to_u16);

    impl_to_primitive!(u32, to_u32);

    impl_to_primitive!(u64, to_u64);

    impl_to_primitive!(u128, to_u128);

    impl_to_primitive!(i8, to_i8);

    impl_to_primitive!(i16, to_i16);

    impl_to_primitive!(i32, to_i32);

    impl_to_primitive!(i64, to_i64);

    impl_to_primitive!(i128, to_i128);

    impl_to_primitive!(f32, to_f32);

    impl_to_primitive!(f64, to_f64);
}

impl<T> Exponential for Complex<T> where T: Real
{
    /// Returns the euler number represented as a complex number
    fn e() -> Self
    {
        Complex { re: T::e(),
                  im: T::zero() }
    }

    ///Exponential function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential, num::Complex};
    ///
    /// let z: Complex<f64> = Complex::new(1.0, 2.0);
    /// let a: Complex<f64> = z.exp();
    /// ```
    fn exp(self: Self) -> Self
    {
        let k: T = self.re.exp();
        Complex { re: k * self.im.cos(),
                  im: k * self.im.sin() }
    }

    ///Logiarithm function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(5.0_f64.powf(0.5_f64).ln(), 2.0_f64.atan());
    ///
    /// assert_eq!(refer, a.ln());
    /// ```
    fn ln(self: Self) -> Self
    {
        Complex { re: self.abs().re.ln(),
                  im: self.arg().re }
    }
}

impl<T> Trigonometry for Complex<T> where T: Real
{
    /// Returns the mathematic constant PI, represented as a complex number
    fn pi() -> Self
    {
        Complex { re: T::pi(),
                  im: T::zero() }
    }

    /// Sinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = (-(-b).exp() * a.sin() - b.exp() * a.sin()) / -2.0;
    /// let im: f64 = ((-b).exp() * a.cos() - b.exp() * a.cos()) / -2.0;
    ///
    /// let uut: Complex<f64> = z.sin();
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn sin(self: Self) -> Self
    {
        let a: Complex<T> = Complex!(-self.im, self.re);
        let b: Complex<T> = Complex!(self.im, -self.re);
        let c: Complex<T> = Complex!(T::zero(), T::one() + T::one());

        (a.exp() - b.exp()) / c
    }

    /// Cosinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = ((-b).exp() * a.cos() + b.exp() * (-a).cos()) / 2.0;
    /// let im: f64 = ((-b).exp() * a.sin() + b.exp() * (-a).sin()) / 2.0;
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// let uut: Complex<f64> = z.cos();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn cos(self: Self) -> Self
    {
        let a: Complex<T> = Complex!(-self.im, self.re);
        let b: Complex<T> = Complex!(self.im, -self.re);
        let c: Complex<T> = Complex!(T::one() + T::one(), T::zero());

        (a.exp() + b.exp()) / c
    }

    /// Tangens function
    ///
    /// # Arguments
    ///
    /// self \in \mathbb{C} \setminus \{ k\pi + \frac{\pi}{2} | k \in \mathbb{Z}
    /// \}
    ///
    /// # Panics
    ///
    /// if the argument bounds are not fulfilled
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let refer: Complex<f64> = z.sin() / z.cos();
    ///
    /// let uut: Complex<f64> = z.tan();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn tan(self: Self) -> Self
    {
        return self.sin() / self.cos();
    }

    /// Cotangens function
    ///
    /// # Arguments
    ///
    /// self: \mathbb{C} \setminus \{ \frac{k * \pi}{2} | k \in \mathbb{Z} \}
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.tan();
    ///
    /// assert_eq!(refer, a.cot());
    /// ```
    fn cot(self: Self) -> Self
    {
        Complex::one() / self.tan()
    }

    /// Secant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.cos();
    ///
    /// assert_eq!(refer, a.sec());
    /// ```
    fn sec(self: Self) -> Self
    {
        Complex::one() / self.cos()
    }

    /// Cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.sin();
    ///
    /// assert_eq!(refer, a.csc());
    /// ```
    fn csc(self: Self) -> Self
    {
        Complex::one() / self.sin()
    }

    /// Inverse sinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(0.4270785863924768, 1.5285709194809995);
    ///
    /// assert_eq!(refer, a.arcsin());
    /// ```
    fn arcsin(self: Self) -> Self
    {
        let mi: Complex<T> = Complex!(T::zero(), -T::one());
        let iz: Complex<T> = Complex!(-self.im, self.re);
        let exp: Complex<T> = Complex!(T::one() / (T::one() + T::one()), T::zero());

        mi * (iz + (Complex::one() - self * self).pow(exp)).ln()
    }

    /// Inverse cosinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0_f64, 0.0_f64) - a.arcsin();
    ///
    /// assert_eq!(refer, a.arccos());
    /// ```
    fn arccos(self: Self) -> Self
    {
        Complex!(T::pi() / (T::one() + T::one()), T::zero()) - self.arcsin()
    }

    /// Inverse tangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(0.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0,
    ///                                        (4.0_f64 / 5.0_f64).atanh() / 2.0_f64);
    ///
    /// assert_eq!(refer, a.arctan());
    /// ```
    fn arctan(self: Self) -> Self
    {
        //		let iz: Complex<T> = Complex!(-self.im, self.re);
        //		let f: Complex<T> = Complex!(T::zero(), T::one() / (T::one() + T::one()));
        //		((Complex::one() - iz).ln() - (Complex::one() + iz).ln()) * f
        //		let k: Complex<T> = Complex!(T::zero(), T::one() + T::one());
        //		((Complex::one() + iz).ln() - (Complex::one() - iz).ln()) / k

        let two: T = T::one() + T::one();
        let re: T;

        if self.re == T::zero()
        {
            if self.im.abs() <= T::one()
            {
                re = T::zero()
            }
            else
            {
                if self.im > T::zero()
                {
                    re = T::pi() / two;
                }
                else
                {
                    re = -T::pi() / two;
                }
            }
        }
        else
        {
            if self.re > T::zero()
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                     + T::pi() / two)
                    / two
            }
            else
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                     - T::pi() / two)
                    / two
            }
        }

        let im: T =
            ((two * self.im) / (self.re * self.re + self.im * self.im + T::one())).artanh() / two;

        Complex!(re, im)
    }

    fn arctan2(self: Self, _other: Self) -> Self
    {
        unimplemented!()
    }

    /// Inverse cotangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arctan();
    ///
    /// assert_eq!(refer, a.arccot());
    /// ```
    fn arccot(self: Self) -> Self
    {
        if self.re == T::zero()
        {
            if self.im == T::one() || self.im == -T::one()
            {
                panic!()
            }
        }
        (Complex::one() / self).arctan()
    }

    /// Inverse secant function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-1, 0, 1}
    ///
    /// # Panics
    ///
    /// iff self = -1 or self = 0 or self = 1
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arccos();
    ///
    /// assert_eq!(refer, a.arcsec());
    /// ```
    fn arcsec(self: Self) -> Self
    {
        if self.im == T::zero()
        {
            if self.re == -T::one() || self.re == T::zero() || self.re == T::one()
            {
                panic!()
            }
        }

        (Complex::one() / self).arccos()
    }

    /// Inverse cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arcsin();
    ///
    /// assert_eq!(refer, a.arccsc());
    /// ```
    fn arccsc(self: Self) -> Self
    {
        (Complex::one() / self).arcsin()
    }
}

impl<T> Power for Complex<T> where T: Real
{
    /// Power
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Power, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let b: Complex<f64> = Complex::new(-2.0_f64, -1.0_f64);
    /// let refer: Complex<f64> = Complex::new(-0.6006033457684014, -0.07399065302898929);
    ///
    /// assert_eq!(refer, a.pow(b));
    /// ```
    fn pow(self: Self, exp: Self) -> Self
    {
        let r: T = self.abs().re;
        let phi: T = self.arg().re;
        let k: T = r.pow(exp.re) * (-exp.im * phi).exp();
        let theta: T = r.ln() * exp.im + exp.re * phi;
        let re: T = k * theta.cos();
        let im: T = k * theta.sin();
        Complex!(re, im)
    }

    fn root(self: Self, _root: Self) -> Self
    {
        unimplemented!();
    }

    fn sqrt(self: Self) -> Self
    {
        unimplemented!();
    }
}

impl<T> Hyperbolic for Complex<T> where T: Real
{
    /// Hyperbolic sine
    fn sinh(self: Self) -> Self
    {
        Complex!(T::zero(), -T::one()) * Complex!(-self.im, self.re).sin()
    }

    /// Hyperbolic cosine
    fn cosh(self: Self) -> Self
    {
        Complex!(-self.im, self.re).cos()
    }

    /// Hyperbolic tangens
    fn tanh(self: Self) -> Self
    {
        self.sinh() / self.cosh()
    }

    /// Hyperbolic cotangens
    fn coth(self: Self) -> Self
    {
        self.cosh() / self.sinh()
    }

    /// Hyperbolic secant
    fn sech(self: Self) -> Self
    {
        Complex!(-self.im, self.re).sec()
    }

    /// Hyperbolic cosecant
    fn csch(self: Self) -> Self
    {
        Complex!(T::zero(), -T::one()) * Complex!(-self.im, self.re).csc()
    }

    /// Hyperbolic inverse sine
    ///
    /// # Arguments
    ///
    /// # Panics
    fn arsinh(self: Self) -> Self
    {
        let p: Complex<T> = Complex!(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self + Complex::one()).pow(p)).ln()
    }

    /// Hyperbolic inverse cosine
    ///
    /// # Argument
    ///
    /// # Panics
    fn arcosh(self: Self) -> Self
    {
        let p: Complex<T> = Complex!(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self - Complex::one()).pow(p)).ln()
    }

    /// Inverse hyperbolic tangent
    ///
    /// # Arguments
    ///
    /// # Panics
    fn artanh(self: Self) -> Self
    {
        let f: Complex<T> = Complex!(T::one() / (T::one() + T::one()), T::zero());
        ((Complex::one() + self) / (Complex::one() - self)).ln() * f
    }

    /// Inverse hyperbolic cosecant
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    fn arcoth(self: Self) -> Self
    {
        let f: Complex<T> = Complex!(T::one() / (T::one() + T::one()), T::zero());
        ((self + Complex::one()) / (self - Complex::one())).ln() * f
    }

    /// Hyperbolic inverse secant
    fn arsech(self: Self) -> Self
    {
        (Complex::one() / self).arcosh()
    }

    // Hyperbolic inverse cosecant
    fn arcsch(self: Self) -> Self
    {
        (Complex::one() / self).arsinh()
    }
}

/// A generic trait for converting a number to a value.
impl<T> FromPrimitive for Complex<T> where T: Real
{
    /// Convert an `i64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i64(_n: i64) -> Self
    {
        unimplemented!();
    }

    /// Convert an `i128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i128(_n: i128) -> Self
    {
        unimplemented!();
    }

    /// Convert an `u64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u64(n: u64) -> Self
    {
        Complex { re: cast::cast(n),
                  im: T::zero() }
    }

    /// Convert an `u128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u128(n: u128) -> Self
    {
        Complex { re: cast::cast(n),
                  im: T::zero() }
    }

    /// Convert a `f64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_f64(n: f64) -> Self
    {
        Complex { re: cast::cast(n),
                  im: T::zero() }
    }
}

/// An interface for casting between machine scalars.
impl<T> NumCast for Complex<T> where T: Real
{
    /// Creates a number from another value that can be converted into
    /// a primitive via the `ToPrimitive` trait.
    fn from<K: ToPrimitive>(n: K) -> Self
    {
        Complex { re: cast::cast(n.to_f64()),
                  im: T::zero() }
    }
}

impl<T> Identity<Addition> for Complex<T> where T: Identity<Addition> + Real
{
    fn id() -> Self
    {
        return Complex::new(Identity::<Addition>::id(), Identity::<Addition>::id());
    }
}

impl<T> Identity<Multiplication> for Complex<T>
    where T: Identity<Multiplication> + Identity<Addition> + Real
{
    fn id() -> Self
    {
        return Complex::new(Identity::<Multiplication>::id(), Identity::<Addition>::id());
    }
}

impl<T> Magma<Addition> for Complex<T> where T: Real
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re + rhs.re, self.im + rhs.im);
    }
}

impl<T> AbsDiffEq for Complex<T>
    where T: Real
{
    type Epsilon = Self;

    fn default_epsilon() -> Self
    {
        unimplemented!()
    }

    fn abs_diff_eq(&self, _other: &Self, _epsilon: Self) -> bool
    {
        unimplemented!();
    }
}

impl<T> MagmaAdd for Complex<T> where T: Real
{
}

impl<T> Magma<Multiplication> for Complex<T> where T: Real
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re * rhs.re - self.im * rhs.im,
                            self.re * rhs.im + self.im * rhs.re);
    }
}

impl<T> MagmaMul for Complex<T> where T: Real
{
}

impl<T> Quasigroup<Addition> for Complex<T> where T: Real
{
}

impl<T> Quasigroup<Multiplication> for Complex<T> where T: Real
{
}

impl<T> AbelianGroup<Addition> for Complex<T> where T: Real
{
}

impl<T> AbelianGroupAdd for Complex<T> where T: Real
{
}

impl<T> AbelianGroup<Multiplication> for Complex<T> where T: Real
{
}
impl<T> AbelianGroupMul for Complex<T> where T: Real
{
}

impl<T> Loop<Addition> for Complex<T> where T: Real
{
}

impl<T> Loop<Multiplication> for Complex<T> where T: Real
{
}

impl<T> CommutativeRing for Complex<T> where T: Real
{
}

impl<T> Ring for Complex<T> where T: Real
{
}

impl<T> Monoid<Addition> for Complex<T> where T: Real
{
}

impl<T> MonoidAdd for Complex<T> where T: Real
{
}

impl<T> Monoid<Multiplication> for Complex<T> where T: Real
{
}
impl<T> MonoidMul for Complex<T> where T: Real
{
}

impl<T> Semigroup<Addition> for Complex<T> where T: Real
{
}

impl<T> SemigroupAdd for Complex<T> where T: Real
{
}

impl<T> Semigroup<Multiplication> for Complex<T> where T: Real
{
}

impl<T> SemigroupMul for Complex<T> where T: Real
{
}

impl<T> Field for Complex<T> where T: Real
{
}

impl<T> Group<Addition> for Complex<T> where T: Real
{
}

impl<T> GroupAdd for Complex<T> where T: Real
{
}

impl<T> Group<Multiplication> for Complex<T> where T: Real
{
}

impl<T> GroupMul for Complex<T> where T: Real
{
}
