use crate::stats::distrib::{Continuous, Normal};
use crate::algebra::linear::{Vector};
use crate::stats::distrib::T as TD;


/// T-Test
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Student%27s_t-test">https://en.wikipedia.org/wiki/Student%27s_t-test</a>
///
/// # Example
/// ```
/// use mathru;
/// use mathru::stats::distrib::{Distribution, Normal};
/// use mathru::stats::test::T;
///
/// let rv1 = Normal::new(1.0, 0.5).random_vector(100);
/// let rv2 = Normal::new(1.0, 0.5).random_vector(100);
///
/// //Test with sample with identical means
/// let mut measure: T = T::test_independence_unequal_variance(&rv1, &rv2);
/// println!("{}", measure.t());
/// measure = T::test_independence_equal_variance(&rv1, &rv2);
/// println!("{}", measure.t());
///
/// // TEst with different equal mean, but unequal variances
/// let rv3 = Normal::new(1.0, 1.5).random_vector(100);
/// measure = T::test_independence_unequal_variance(&rv1, &rv3);
/// println!("{}", measure.t());
/// measure = T::test_independence_equal_variance(&rv1, &rv3);
/// println!("{}", measure.t());
///
/// // When the sample size is not equal anymore
/// //the equal variance t-statistic is no longer equal to the unequal variance t-statistic:
///	let rv4 = Normal::new(2.0, 0.5).random_vector(300);
/// measure = T::test_independence_unequal_variance(&rv1, &rv4);
/// println!("{}", measure.t());
/// measure = T::test_independence_equal_variance(&rv1, &rv4);
/// println!("{}", measure.t());
///
/// //t-Test with different mean, variance and sample size
///	let rv5 = Normal::new(2.0, 1.0).random_vector(300);
/// measure = T::test_independence_unequal_variance(&rv1, &rv5);
/// println!("{}", measure.t());
/// measure = T::test_independence_equal_variance(&rv1, &rv5);
/// println!("{}", measure.t());
/// ```
pub struct T
{
	p: f64,
	t: f64
}

impl T
{
	pub fn t(self: &Self) -> f64
	{
		self.t
	}

	pub fn p_value(self: &Self) -> f64
	{
		return self.p
	}


	/// Calculates the T-test for the means of two independent samples of scores
	///
	/// This is a two-sided test for the null hypothesis that two independent samples have identical expected values.
	/// It is assumed, that the populations have identical variances.
	pub fn test_independence_equal_variance(x: &Vector<f64>, y: &Vector<f64>) -> T
	{
		let (_, n_x) : (usize, usize) = x.dim();
		let (_, n_y) : (usize, usize) = y.dim();

		let x_dist: Normal = Normal::from_data(&x);
		let y_dist: Normal = Normal::from_data(&y);

		let mean_x: f64 = x_dist.mean();
		let mean_y: f64 = y_dist.mean();

		let df: usize = n_x + n_y - 2;

		let s_x_squared: f64 = x_dist.variance();
		let s_y_squared: f64 = y_dist.variance();

		let nomin: f64 = ((n_x - 1) as f64) * s_x_squared + ((n_y - 1) as f64) * s_y_squared;
		let denom: f64 = (df) as f64;

		let s_p: f64 = (nomin / denom).sqrt();

		let t: f64 = (mean_x - mean_y) / (s_p * (1.0 / (n_x as f64) + 1.0 / (n_y as f64)).sqrt());
		T
		{
			p: 0.0,
			t: t
		}
	}

	/// Calculates the T-test for the means of two independent samples of scores
	///
	/// This is a two-sided test for the null hypothesis that two independent samples have identical expected values.
	/// It is assumed, that the populations have NOT identical variances. It performs the Welch’s t-test
	pub fn test_independence_unequal_variance(x: &Vector<f64>, y: &Vector<f64>) -> T
	{
		let (_, n_x) : (usize, usize) = x.dim();
		let (_, n_y) : (usize, usize) = y.dim();

		let x_dist: Normal = Normal::from_data(&x);
		let y_dist: Normal = Normal::from_data(&y);

		let mean_x: f64 = x_dist.mean();
		let mean_y: f64 = y_dist.mean();

		let s_x_squared: f64 = x_dist.variance();
		let s_y_squared: f64 = y_dist.variance();

		let term1: f64 = s_x_squared / (n_x as f64) + s_y_squared / (n_y as f64);

		let df: f64 =  term1 * term1 / (s_x_squared * s_x_squared /
		 ((n_x * n_x * (n_x - 1)) as f64) +	s_y_squared * s_y_squared / ((n_y * n_y * (n_y - 1)) as f64));

		let s_p: f64 = term1.sqrt();

		let t: f64 = (mean_x - mean_y) / s_p ;

		let p: f64 = 2.0 * TD::new(df).cdf(-t.abs());
		T
		{
			p: p,
			t: t
		}
	}

}
