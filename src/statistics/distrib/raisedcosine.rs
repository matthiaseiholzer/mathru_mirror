use crate::algebra::abstr::Real;
use crate::statistics::distrib::Continuous;

/// Raised Cosine distribution
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Raised_cosine_distribution">https://en.wikipedia.org/wiki/Raised_cosine_distribution</a>
pub struct RaisedCosine<T>
{
    mu: T,
    s: T,
}

impl<T> RaisedCosine<T> where T: Real
{
    /// Creates a probability distribution
    ///
    /// # Arguments
    ///
    /// * `mu`
    /// * `s` > 0.0
    ///
    /// # Panics
    ///
    /// if s < 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::RaisedCosine;
    /// use std::f64::consts::PI;
    ///
    /// let mu: f64 = PI;
    /// let s: f64 = 0.5 * PI;
    /// let distrib: RaisedCosine<f64> = RaisedCosine::new(mu, s);
    /// ```
    pub fn new(mu: T, s: T) -> RaisedCosine<T>
    {
        if s < T::zero()
        {
            panic!();
        }
        RaisedCosine { mu, s }
    }
}

impl<T> Continuous<T> for RaisedCosine<T> where T: Real
{
    /// Probability density function
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x
    ///
    /// # Panics
    ///
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, RaisedCosine};
    ///
    /// let distrib: RaisedCosine<f64> = RaisedCosine::new(-1.2, 1.5);
    /// let x: f64 = 5.0;
    /// let p: f64 = distrib.pdf(x);
    /// ```
    fn pdf(self: &Self, x: T) -> T
    {
        if (self.mu - self.s) <= x && x < (self.mu + self.s)
        {
            return (T::one() + (T::pi() * (x - self.mu) / self.s).cos())
                   / (T::from_f64(2.0) * self.s);
        }

        return T::zero();
    }

    /// Cumulative distribution function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, RaisedCosine};
    /// use std::f64::consts::PI;
    ///
    /// let distrib: RaisedCosine<f64> = RaisedCosine::new(1.0, PI);
    /// let x: f64 = PI / 2.0;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf(self: &Self, x: T) -> T
    {
        if (self.mu - self.s) <= x && x <= (self.mu + self.s)
        {
            let k: T = (x - self.mu) / self.s;
            return (T::one() + k + T::one() / T::pi() * (k * T::pi()).sin()) / T::from_f64(2.0);
        }
        else
        {
            if x < (self.mu - self.s)
            {
                return T::zero();
            }
            else
            {
                return T::one();
            }
        }
    }

    /// Quantile function of inverse cdf
    fn quantile<'a, 'b>(self: &'a Self, _p: T) -> T
    {
        unimplemented!();
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, RaisedCosine};
    ///
    /// let distrib: RaisedCosine<f64> = RaisedCosine::new(-2.0, 0.5);
    /// let mean: f64 = distrib.mean();
    /// ```
    fn mean<'a>(self: &'a Self) -> T
    {
        return self.mu;
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, RaisedCosine};
    /// use std::f64::consts::PI;
    ///
    /// let distrib: RaisedCosine<f64> = RaisedCosine::new(2.0, PI);
    /// let var: f64 = distrib.variance();
    /// ```
    fn variance(self: &Self) -> T
    {
        return self.s * self.s * (T::from_f64(1.0 / 3.0) - T::from_f64(2.0) / (T::pi() * T::pi()));
    }

    ///
    fn skewness(self: &Self) -> T
    {
        return T::zero();
    }

    /// Median is the value separating the higher half from the lower half of a
    /// probability distribution.
    fn median(self: &Self) -> T
    {
        return self.mu;
    }

    ///
    fn entropy(self: &Self) -> T
    {
        unimplemented!();
    }
}
