use crate::statistics::distrib::Continuous;
use rand;
use crate::algebra::abstr::Real;

/// Exponential distribution
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Exponential_distribution">https://en.wikipedia.org/wiki/Exponential_distribution</a>
///
pub struct Exponential<T>
{
	lambda: T
}

impl<T> Exponential<T>
	where T: Real
{
	/// Creates a probability distribution
    ///
    /// # Arguments
    ///
    /// * `lambda`  > 0.0
    ///
    /// # Panics
    ///
    /// if lambda <= 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::Exponential;
    ///
    /// let distrib: Exponential<f64> = Exponential::new(0.3);
    /// ```
	pub fn new(lambda: T) -> Exponential<T>
	{
		if lambda <= T::zero()
		{
			panic!()
		}

		Exponential
		{
			lambda: lambda
		}
	}

	pub fn from_data<'a>(data: &'a Vec<T>) -> Self
    {
        let lambda: T = T::one() / Exponential::calc_mean(data);

        return Exponential::new(lambda)
    }

    fn calc_mean<'a>(data: &'a Vec<T>) -> T
    {
        let mut sum: T = T::zero();

        for x in data.iter()
        {
            sum += *x;
        }

        return sum / T::from_u64(data.len() as u64).unwrap();
    }
}

impl<T> Continuous<T> for Exponential<T>
	where T: Real
{
    /// Probability density function
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x &isin; &#x2115 | x > 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Exponential};
    ///
    /// let distrib: Exponential<f64> = Exponential::new(0.3);
    /// let x: f64 = 5.0;
    /// let p: f64 = distrib.pdf(x);
    /// ```
    fn pdf(self: &Self, x: T) -> T
	{
		if x < T::zero()
		{
			return T::zero()
		}

		return self.lambda * (-self.lambda * x).exp();
	}

    /// Cumulative distribution function
    ///
    /// # Arguments
    ///
    /// * `x` Random variable
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Exponential};
    ///
    /// let distrib: Exponential<f64> = Exponential::new(0.3);
    /// let x: f64 = 0.4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf(self: &Self, x: T) -> T
	{
		if x < T::zero()
		{
			return T::zero()
		}

		return T::one() - (-x * self.lambda).exp();
	}

	/// Quantile function of inverse cdf
    fn quantile(self: &Self, p: T) -> T
    {
    	return -(T::one() - p).ln() / self.lambda;
    }

  	/// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Exponential};
    ///
    /// let distrib: Exponential<f64> = Exponential::new(0.2);
    /// let mean: f64 = distrib.mean();
    /// ```
	fn mean(self: &Self) -> T
	{
		return T::one() / self.lambda
	}

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Exponential};
    ///
    /// let distrib: Exponential<f64> = Exponential::new(0.2);
    /// let var: f64 = distrib.variance();
    /// ```
	fn variance(self: &Self) -> T
	{
		return T::one() / self.lambda.pow(&T::from_u8(2).unwrap())
	}
}


impl<T> Exponential<T>
	where T: Real
{
 	pub fn random(self: &Self) -> T
    {
  		let y: T = T::from_f64(rand::random::<f64>()).unwrap();
   		let p: T = self.quantile(y);

   		return p;
	}
}