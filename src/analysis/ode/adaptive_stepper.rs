//! Adaptive step size stepper

use std::default::Default;
use crate::algebra::linear::{Vector};
use crate::algebra::abstr::Real;
use super::ExplicitODE;
use super::explicit_method::ExplicitAdaptiveMethod;

/// Adaptive step size stepper
pub struct AdaptiveStepper<T>
{
     /// Step size
    n_max: u32,
    h_0: T,
    fac: T,
    fac_min: T,
    fac_max: T,
    /// abs_tol: Absolute tolerance. This is the tolerance on local error estimates, not necessarily the global error.
    /// Defaults to 1e-6.
    abs_tol: T,
    ///reltol: Relative tolerance. This is the tolerance on local error estimates, not necessarily the global error.
    /// Defaults to 1e-3.
    rel_tol: T,
}

impl<T> Default for AdaptiveStepper<T>
    where T: Real
{
    fn default() -> AdaptiveStepper<T>
    {
        return AdaptiveStepper::new(1000, T::from_f64(0.02).unwrap(), T::from_f64(0.8).unwrap(), T::from_f64(0.001).unwrap(), T::from_f64
        (3.0).unwrap(), T::from_f64(10e-6).unwrap(), T::from_f64(10e-3).unwrap());
    }
}

impl<T> AdaptiveStepper<T>
    where T: Real
{
    /// Creates an instance with the given step siz
    ///a
    ///
    /// # Param
    ///
    /// * 'fac_min':
    /// *'fac_max': 1.5 <= fac_max <= 5.0
    pub fn new(n_max: u32, h_0: T, fac: T, fac_min: T, fac_max: T, abs_tol: T, rel_tol: T) -> AdaptiveStepper<T>
    {
        return AdaptiveStepper
        {
            n_max: n_max,
            h_0: h_0,
            fac: fac,
            fac_min: fac_min,
            fac_max: fac_max,
            abs_tol: abs_tol,
            rel_tol: rel_tol,
        }
    }

    /// Returns the aboslute tolerance
    pub fn get_abs_tol(self: &Self) -> &T
    {
        return &self.abs_tol;
    }

    /// Returns the relative tolerance
    pub fn get_rel_tol(self: &Self) -> &T
    {
        return &self.rel_tol;
    }

    /// Sets the aboslute tolerance
    ///
    /// # Parameters
    ///
    /// * 'abs_tol': abs_tol >= 0.0
    /// # Panics
    ///
    /// if 'abs_tol' < 0.0
    pub fn set_abs_tol(self: &mut Self, abs_tol: T)
    {
        if abs_tol < T::zero()
        {
            panic!();
        }
        self.abs_tol = abs_tol;
    }

    /// Sets the relative tolerance
    ///
    /// # Parameters
    ///
    /// * 'rel_tol': rel_tol >= 0.0
    /// # Panics
    ///
    /// if 'rel_tol' < 0.0
    pub fn set_rel_tol(self: &mut Self, rel_tol: T)
    {
        if rel_tol < T::zero()
        {
            panic!();
        }
        self.rel_tol = rel_tol;
    }


    /// Solves `func` using the 4th order Runge-Kutta-Fehlberg algorithm.
    ///
    /// # Arguments
    ///
    /// * 'func' is an explict oridnary diffential equation
    /// * 'init' is the initial value at the time 't_start'
    /// * 't_span' Time span t_span.0 = t_start, t_span.1 = t_stop
    ///
    /// # Return
    ///
    /// The solver returns a vector and a matrix, containing the times used in each step of the
    /// algorithm and the respectful values for that time.
    ///
    /// # Panic
    ///
    /// if t_span.0 > t_span.1
    pub fn solve<F, M>(self: &Self, prob: &F, method: &M) -> Result<(Vec<T>, Vec<Vector<T>>), &'static str>
        where F: ExplicitODE<T>,
               M: ExplicitAdaptiveMethod<T>
    {
        let t_span: (T, T) = prob.time_span();
        let t_start: T = t_span.0;
        let t_stop: T = t_span.1;
        if t_start > t_stop
        {
            panic!();
        }

        let order: (u8, u8) = method.order();
        let q: T = T::from_u8(order.0.max(order.1)).unwrap();
        let l: T = T::one() / (q + T::one());
        let mut x_n: Vector<T> = prob.init_cond();
        let mut t_n: T = t_start;
        let mut h: T = self.h_0;

        let mut t_vec: Vec<T> = Vec::new();
        t_vec.push(t_n);

        let mut res_vec: Vec<Vector<T>> = Vec::new();
        res_vec.push(x_n.clone());

        let mut n: u32 = 0;

        while n < self.n_max && t_n < t_stop
        {
            h = h.min(t_stop - t_n);
            //
            let (x_n_new, x_ne): (Vector<T>, Vector<T>) = method.do_step(prob, &t_n, &x_n, &h);
            let err: T = self.calc_error(&x_n_new, &x_ne);


            if err <= T::one()
            {
                t_n = t_n + h;
                x_n = x_n_new;
                t_vec.push(t_n);
                res_vec.push(x_n.clone());
				n = n + 1;
            }

            if err != T::zero()
            {
                let mut s: T = self.fac * (T::one() / err).pow(&l);

                if s < self.fac_min
                {
                    s = self.fac_min;
                }

                if s > self.fac_max
                {
                    s = self.fac_max;
                }

                h = s * h;
            }

        }
        return Ok((t_vec, res_vec));
    }

    fn calc_error(self: &Self, y: &Vector<T>, y_hat: &Vector<T>) -> T
    {
        let (_m, n) = y.dim();

        let mut sum: T = T::zero();

        for i in 0..n
        {
            let y_i: T = *y.get(i);
            let y_hat_i: T = *y_hat.get(i);
            let sc: T = self.abs_tol + y_i.abs().max(y_hat_i.abs()) * self.rel_tol;

            let k: T = (y_i - y_hat_i) / sc;
            sum += k * k;

        }

        let p = (sum / T::from_f64(n as f64).unwrap()).pow(&T::from_f64(0.5).unwrap());
        return p;
    }


}