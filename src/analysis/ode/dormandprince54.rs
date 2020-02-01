use crate::algebra::linear::{Vector};
use crate::algebra::abstr::Sign;
use crate::algebra::abstr::Real;
use super::explicit_method::{ExplicitAdaptiveMethod};
use super::ExplicitODE;

/// Solves an ordinary differential equation using the 4th order Runge-Kutta-Dormand-Prince algorithm.
///
///<a href="https://en.wikipedia.org/wiki/Dormand-Prince_method">https://en.wikipedia.org/wiki/Dormand-Prince_method</a>
pub struct DormandPrince54<T>
{
    abs_tol: T,
    h_0: T,
    n_max: u32,
}

impl<T> DormandPrince54<T>
    where T: Real
{
    /// Creates a DormandPrince54 instance, also known as explicit Runge-Kutta method of order 5(4) with step-size control
    pub fn new(abs_tol: T, begin_step_size: T, n_max: u32) -> DormandPrince54<T>
    {
        return DormandPrince54
        {
            abs_tol: abs_tol,
            h_0: begin_step_size,
            n_max: n_max,
        }
    }

    pub fn get_begin_step_size(self: &Self) -> &T
    {
        return &self.h_0;
    }

    pub fn set_begin_step_size(self: &mut Self, step_size: T)
    {
        self.h_0 = step_size
    }

    pub fn solve<F>(self: &Self, prob: &F) -> Result<(Vec<T>, Vec<Vector<T>>), &'static str>
        where F: ExplicitODE<T>,
    {
        let t_span: (T, T) = prob.time_span();
        let t_start: T = t_span.0;
        let t_stop: T = t_span.1;
        if t_start > t_stop
        {
            panic!();
        }

        let mut x_n: Vector<T> = prob.init_cond();
        let mut t_n: T = t_start;
        let mut h_n: T = self.h_0;

        let mut t_vec: Vec<T> = Vec::new();
        t_vec.push(t_n);

        let mut res_vec: Vec<Vector<T>> = Vec::new();
        res_vec.push(x_n.clone());

        let mut n: u32 = 0;

        while n < self.n_max && t_n < t_stop
        {
            h_n = h_n.min(t_stop - t_n);
            //
            let (x_n_new, x_ne): (Vector<T>, Vector<T>) = self.do_step(prob, &t_n, &x_n, &h_n);
            let error_n: T = self.calc_error(&x_n_new, &x_ne);

            if error_n <= self.abs_tol
            {
                t_n = t_n + h_n;
                x_n = x_n_new;
                t_vec.push(t_n);
                res_vec.push(x_n.clone());
				n = n + 1;
            }

            //Update step size
            h_n = self.calc_step_size(h_n, error_n);
        }

        if t_n < t_stop
        {
            return Err("Maxmimum number of iterations reached");
        }
        return Ok((t_vec, res_vec));
    }

    /// calculate new step size
    fn calc_step_size(self: &Self, h_n: T, error_n: T) -> T
    {
        return T::from_f64(0.9).unwrap() * h_n * (self.abs_tol / error_n).pow(&(T::one() / T::from_f64(5.0).unwrap()))
    }

    ///
    /// ```math
    /// \lvert \lvert y - \hat{y} \rvert \rvert_{\infty}
    /// ```
    pub fn calc_error(self: &Self, y: &Vector<T>, y_hat: &Vector<T>) -> T
    {
        let diff: Vector<T> = (y - y_hat).abs();

        let idx_max: usize = diff.argmax();
        return *diff.get(idx_max);
    }
}

impl<T> ExplicitAdaptiveMethod<T> for DormandPrince54<T>
    where T: Real
{

    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>
    {
        // k_1 = hf(t_n, x_n)
        let k_1: Vector<T> = &prob.func(t_n, x_n) * h;

        // k_2 = h f(t_n + h/5, x_n + k1/5);
        let k_2: Vector<T> = &prob.func(&(*t_n + *h / T::from_f64(1.0 / 5.0).unwrap()), &(x_n + &(&k_1 / &T::from_f64(1.0
        / 5.0).unwrap()))) * h;

        //k_3 = h f(t_n + h3/10, x_n + k_1*3/40 + k2*9/40)
        let k_31: Vector<T> = x_n + &(&k_1 * &T::from_f64(3.0 / 40.0).unwrap());
        let k_32: Vector<T> = k_31 + (&k_2 * &T::from_f64(9.0 / 40.0).unwrap());
        let k_3: Vector<T> = &prob.func(&(*t_n + *h * T::from_f64(3.0 / 10.0).unwrap()),  &k_32) * h;

        // k_4 = h f(t_n + 4/5h, x_n + 44/45k_1 - 56/15k_2 + 32/9k_3)
        let k_41: Vector<T> = x_n + &(&k_1 * &T::from_f64(44.0 / 45.0).unwrap());
        let k_42: Vector<T> = k_41 - (&k_2 * &T::from_f64(56.0 / 15.0).unwrap());
        let k_43: Vector<T> = k_42 + (&k_3 * &T::from_f64(32.0 / 9.0).unwrap());
        let k_4: Vector<T> = &prob.func(&(*t_n + *h * T::from_f64(4.0 / 5.0).unwrap()), &k_43) * h;

        // k_5 = h f(t_n + h, x_n + 19372/6561k_1 - 25360/2187k_2 + 64448/6561k_3 - 212/729k_4)
        let k_51: Vector<T> = x_n + &(&k_1 * &T::from_f64(19372.0 / 6561.0).unwrap());
        let k_52: Vector<T> = k_51 - (&k_2 * &T::from_f64( 25360.0 / 2187.0).unwrap());
        let k_53: Vector<T> = k_52 + (&k_3 * &T::from_f64(64448.0 / 6561.0).unwrap());
        let k_54: Vector<T> = k_53 - (&k_4 * &T::from_f64(212.0 / 729.0).unwrap());
        let k_5: Vector<T> = &prob.func(&(*t_n + *h * T::from_f64(8.0 / 9.0).unwrap()) , &k_54) * h;

        // k_6 = h f(t_n + h, x_n + 9017/3168k_1 - 355/33k_2 + 46732/5247k_3 + 49/176k_4 - 5103/18656k_5)
        let k_61: Vector<T> = x_n + &(&k_1 * &T::from_f64(9017.0 / 3168.0).unwrap());
        let k_62: Vector<T> = k_61 - (&k_2 * &T::from_f64(355.0 / 33.0).unwrap());
        let k_63: Vector<T> = k_62 + (&k_3 * &T::from_f64(46732.0 / 5247.0).unwrap());
        let k_64: Vector<T> = k_63 + (&k_4 * &T::from_f64(49.0 / 176.0).unwrap());
        let k_65: Vector<T> = k_64 - (&k_5 * &T::from_f64(5103.0 / 18656.0).unwrap());
        let k_6: Vector<T> = &prob.func(&(*t_n + *h), &k_65) * h;

        // k_7 = h f(t_n + h, x_n + 35/384k_1 + 500/1113k_3 + 125/192k_4 - 2187/6784k_5 + 11/84k_6)
        let k_71: Vector<T> = x_n + &(&k_1 * &T::from_f64(35.0 / 384.0).unwrap());
        let k_72: Vector<T> = k_71 + (&k_3 * &T::from_f64(500.0 / 1113.0).unwrap());
        let k_73: Vector<T> = k_72 + (&k_4 * &T::from_f64(125.0 / 192.0).unwrap());
        let k_74: Vector<T> = k_73 - (&k_5 * &T::from_f64(2187.0 / 6784.0).unwrap());
        let k_75: Vector<T> = k_74 + (&k_6 * &T::from_f64(11.0 / 84.0).unwrap());
        let k_7: Vector<T> = &prob.func(&(*t_n + *h), &k_75) * h;

        let rkdp4: Vector<T> = k_75;

        // y_(n +1) = x_n + 5179/57600k_1 + 7571/16695k_3 + 393/640k_4 - 92097/339200k_5 + 187/2100k_6 + 1/40k_7
        let rkdp5_1: Vector<T> = x_n + &(&k_1 * &T::from_f64(5179.0 / 57600.0).unwrap());
        let rkdp5_2: Vector<T> = rkdp5_1 + (&k_3 * &T::from_f64(7571.0 / 16695.0).unwrap());
        let rkdp5_3: Vector<T> = rkdp5_2 + (&k_4 * &T::from_f64(393.0 / 640.0).unwrap());
        let rkdp5_4: Vector<T> = rkdp5_3 - (&k_5 * &T::from_f64(92097.0 / 339200.0).unwrap());
        let rkdp5_5: Vector<T> = rkdp5_4 + (&k_6 * &T::from_f64(187.0 / 2100.0).unwrap());
        let rkdp5: Vector<T> = rkdp5_5 + (&k_7 * &T::from_f64(1.0 / 40.0).unwrap());
        return (rkdp4, rkdp5);
    }


    fn order(self: &Self) -> (u8, u8)
    {
        return (4, 5);
    }

}