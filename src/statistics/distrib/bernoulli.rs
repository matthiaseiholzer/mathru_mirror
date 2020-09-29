//! Bernoulli Distribution
use crate::{algebra::abstr::Real, statistics::distrib::Discrete};
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Bernoulli distribution
///
/// Fore more information:
///
/// <a href="https://en.wikipedia.org/wiki/Bernoulli_distribution">https://en.wikipedia.org/wiki/Bernoulli_distribution</a>
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Bernoulli<T>
{
    p: T,
}

impl<T> Bernoulli<T> where T: Real
{
    /// Create a probability distribution with p(X=1) = p
    ///
    /// # Arguments
    ///
    /// * `p` Probability that random varibale X=1, 0.0 <= p <= 1.0
    ///
    /// # Panics
    ///
    /// if p < 0 || p > 1.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::Bernoulli;
    ///
    /// let distrib: Bernoulli<f64> = Bernoulli::new(0.2);
    /// ```
    pub fn new(p: T) -> Bernoulli<T>
    {
        if p < T::zero() || p > T::one()
        {
            panic!()
        }

        return Bernoulli { p };
    }
}

impl<T> Discrete<T, u8, T> for Bernoulli<T> where T: Real
{
    /// Probability mass function of the Bernoulli distribution
    ///
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x &isin; {0, 1}
    ///
    /// # Panics
    ///
    /// if x &notin; {0, 1}
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Bernoulli, Discrete};
    ///
    /// let distrib: Bernoulli<f64> = Bernoulli::new(0.2);
    /// let x: u8 = 0;
    /// let p: f64 = distrib.pmf(x);
    /// ```
    fn pmf<'a>(self: &'a Self, x: u8) -> T
    {
        if (x == 1) || (x == 0)
        {
            if x == 0
            {
                return T::one() - self.p;
            }
            else
            {
                return self.p;
            }
        }
        else
        {
            panic!()
        }
    }

    /// Cumulative distribution function of the Bernoulli distribution
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x &isin; {0, 1}
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Bernoulli, Discrete};
    ///
    /// let distrib: Bernoulli<f64> = Bernoulli::new(0.2);
    /// let x: f64 = 0.4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf<'a>(self: &'a Self, x: T) -> T
    {
        if x >= T::one()
        {
            return T::one();
        }

        if x <= T::zero()
        {
            return T::zero();
        }
        else
        {
            return T::one() - self.p;
        }
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Bernoulli, Discrete};
    ///
    /// let distrib: Bernoulli<f64> = Bernoulli::new(0.2);
    /// let mean: f64 = distrib.mean();
    /// ```
    fn mean<'a>(self: &'a Self) -> T
    {
        self.p
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Bernoulli, Discrete};
    ///
    /// let distrib: Bernoulli<f64> = Bernoulli::new(0.2);
    /// let var: f64 = distrib.variance();
    /// ```
    fn variance<'a>(self: &'a Self) -> T
    {
        self.p * (T::one() - self.p)
    }
}
