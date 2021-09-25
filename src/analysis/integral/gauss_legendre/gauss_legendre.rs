use crate::algebra::abstr::Real;
use crate::analysis::integral::gauss_legendre::root_weight::RootWeight;

///
/// ```math
/// \int_{a}^{b}f(x)\,dx = \frac{b - a}{2}\sum_{i=1}^{n}f(\frac{b - a}{2}x_i + \frac{a + b}{2})a_i
/// ```
///
/// <https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss-Legendre_quadrature>
pub struct GaussLegendre<T>
{
    root_weight: RootWeight<T>
}


impl<T> GaussLegendre<T>
    where T: Real
{

    /// # Arguments
    ///
    ///
    /// # Panics
    ///
    /// Panics if n < 1 || n > 9
    ///
    /// # Examples
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::analysis::integral::gauss_legendre::GaussLegendre;
    ///
    /// let gl: GaussLegendre<f64> = GaussLegendre::new(1);
    /// let f: fn(f64) -> f64 = | x | {x};
    ///
    /// let integral: f64 = gl.integrate(f, 2.0, 4.0);
    ///
    /// assert_relative_eq!(integral, 6.0)
    /// # }
    /// ```
    pub fn new(n: u8) -> GaussLegendre<T>
    {
        if n < 1 || n > 9
        {
            panic!("n is not within the limits");
        }


        return GaussLegendre{
            root_weight: RootWeight::new(n)
        }
    }

    ///
    /// ```math
    /// \int_{a}^{b}f(x)\,dx = \frac{b - a}{2}\sum_{i=1}^{n}f(\frac{b - a}{2}x_i + \frac{a + b}{2})a_i
    /// ```
    ///
    /// # Arguments
    /// * a
    /// * b
    /// #
    pub fn integrate<F>(self: &Self, f: F, a: T, b: T) -> T
        where
            F: Fn(T) -> T,
    {

        let sum = self.root_weight.clone().fold(
        T::zero(),
        |s,  (x_i, a_i)  | {
            let x = (b - a) / T::from_f64(2.0) * x_i + (a + b) / T::from_f64(2.0);
            s + f(x) * a_i
        }
        );

        (b - a) / T::from_f64(2.0) * sum
    }
}
