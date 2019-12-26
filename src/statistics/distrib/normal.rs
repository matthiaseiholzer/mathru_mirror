use crate::statistics::distrib::{Distribution, Continuous};
use crate::algebra::linear::{Vector};
use crate::special;
use rand;
use std::f64;

/// Normal distribution
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Normal_distribution">https://en.wikipedia.org/wiki/Normal_distribution</a>
///
pub struct Normal
{
    mean: f64,
    variance: f64,
}

impl Normal
{
     /// Creates a probability distribution
    ///
    /// # Arguments
    ///
    /// * `mean`: Expected value
    /// * `variance`: variance > 0.0
    ///
    /// # Panics
    ///
    /// if variance <= 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::Normal;
    ///
    /// let distrib: Normal = Normal::new(0.3, 0.2);
    /// ```
    pub fn new(mean: f64, variance: f64) -> Self
    {
        if variance <= 0.0
        {
            panic!();
        }

        Normal
        {
            mean,
            variance
        }
    }

    ///
    /// It is assumed that data are normal distributed.
    ///
    /// data.len() >= 2
    ///
    pub fn from_data<'a>(data: &'a Vector<f64>) -> Self
    {
        let (_m, n): (usize, usize) = data.dim();
        if n < 2
        {
                panic!()
        }

        let mean: f64 = Normal::calc_mean(data);
        let variance: f64 = Normal::calc_variance(data, mean);

        return Normal::new(mean, variance)
    }

    fn calc_mean<'a>(data: &'a Vector<f64>) -> f64
    {
        let (_m, n) = data.dim();
        let mut sum: f64 = 0.0;

        for x in data.iter()
        {
            sum = sum + *x;
        }

        return sum / (n as f64)
    }

    fn calc_variance<'a>(data: &'a Vector<f64>, mean: f64) -> f64
    {
        let (_m, n): (usize, usize) = data.dim();
        let mut sum: f64 = 0.0;

        for x in data.iter()
        {
            sum += (x - mean).powi(2);
        }

        return sum / ((n - 1) as f64)
    }
}

impl Continuous<f64, f64> for Normal
{

    /// Probability density function
    ///
    /// # Arguments
    ///
    /// * `x`: Random variable x &isin; &#x2115
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal = Normal::new(0.3, 0.2);
    /// let x: f64 = 5.0;
    /// let p: f64 = distrib.pdf(x);
    /// ```
    fn pdf<'a, 'b>(self: &'a Self, x: f64) -> f64
    {
        let z: f64 = -0.5 * ((x - self.mean) / self.variance).powf(2.0);
        let f: f64 = 1.0 / (self.variance * (2.0 * f64::consts::PI).sqrt());
        f * z.exp()
    }

    /// Cumulative distribution function
    ///
    /// # Arguments
    ///
    /// * `x`: Random variable
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal = Normal::new(0.3, 0.2);
    /// let x: f64 = 0.4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf<'a, 'b>(self: &'a Self, x: f64) -> f64
    {
        let k: f64 = (x - self.mean) / ((2.0 * self.variance).sqrt());
        let prob: f64 =  0.5 * (1.0 + special::erf(k));
        prob
    }

    /// Quantile: function of inverse cdf
    ///
    /// The Percentage Points of the Normal Distribution
    ///  Author(s): Michael J. Wichura
    /// Year 1988
    /// Journal of the Royal Statistical Society
    /// 0.0 < p < 1.0
    ///
    /// # Panics
    ///
    /// if  p <= 0.0 || p >= 1.0
    fn quantile<'a>(self: &'a Self, p: f64) -> f64
    {
        if p <= 0.0 || p >= 1.0
        {
            panic!();
        }

        let mut ppnd16: f64;
        let mut r: f64;
        let q: f64 = p - 0.5;
        if q.abs() <= 0.425
        {
            let r: f64 = 0.180625 - q * q;
            ppnd16 = q * Normal::r1(r);
        }
        else
        {
            if q < 0.0
            {
                r = p
            }
            else
            {
                r = 1.0 - p
            }
            if r <= 0.0
            {
                ppnd16 = 0.0;
                return ppnd16
            }
            r = (-r.ln()).sqrt();
            if r <= 5.
            {
                r -= 1.6;
                ppnd16 = Normal::r2(r);
            }
            else
            {
                r -= 5.0;
                ppnd16 = Normal::r3(r);
            }
            if q <= 0.0
            {
                ppnd16 = -ppnd16;

            }
        }

        return self.mean + self.variance.sqrt() * ppnd16;
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru;
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal = Normal::new(0.0, 0.2);
    /// let mean: f64 = distrib.mean();
    /// ```
	fn mean<'a>(self: &'a Self) -> f64
    {
        return self.mean
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru;
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal = Normal::new(0.0, 0.2);
    /// let var: f64 = distrib.variance();
    ///
    /// ```
	fn variance<'a>(self: &'a Self) -> f64
    {
        return self.variance
    }
}


impl Distribution for  Normal
{
    ///
    ///  See Knuth The Art of Computer Programming Vol 2 3.4.1 C Algorithm P
    ///
    fn random(self: &Self) -> f64
    {
        let mut s : f64 = 1.0;
        let mut v1 : f64 = 1.0;
        let mut v2 : f64;
        while s >= 1.0
        {
            let u1 : f64 = rand::random::<f64>();
            let u2 : f64 = rand::random::<f64>();
            v1 = 2.0 * u1 - 1.0;
            v2 = 2.0 * u2 - 1.0;
            s =  v1 * v1 + v2 * v2
        }
        let x1: f64 = v1 * (-2.0 * s.ln() / s).sqrt();
        x1 * self.variance.sqrt() + self.mean
    }
}


impl Normal
{

    fn r1(r: f64) -> f64
    {
        let value: f64 = (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
        value
    }

    fn r2(r: f64) -> f64
    {
        let value: f64 = (((((((r * 7.7454501427834140764e-4 +
                       0.0227238449892691845833) * r + 0.24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + 0.0151986665636164571966) * r +
                       0.14810397642748007459) * r + 0.68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);

        value
    }

    fn r3(r: f64) -> f64
    {
        let value: f64 = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      0.0012426609473880784386) * r + 0.026532189526576123093) *
                    r + 0.29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + 0.0148753612908506148525)
                     * r + 0.13692988092273580531) * r +
                    0.59983220655588793769) * r + 1.);

        return value
    }
}