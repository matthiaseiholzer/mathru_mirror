//! Solves an implicit ODE equation using backward Euler.
use super::{implicit_method::ImplicitFixedStepSizeMethod, ImplicitODE};
use crate::{
    algebra::{
        abstr::Real,
        linear::{Matrix, Vector},
    },
    analysis::{
        differential_equation::ordinary::fixed_stepper::ImplicitFixedStepper, Function, Jacobian,
        NewtonRaphson,
    },
};
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Solves an ODE using backward Euler
///
/// <a href="https://en.wikipedia.org/wiki/Backward_Euler_method">https://en.wikipedia.org/wiki/Backward_Euler_method</a>
/// # Example
///
/// For this example, we want to solve the following stiff ordinary
/// differiential equation:
/// ```math
/// 0 = -4(y(t) -2) - y(t)^{'} = f(t, y, y^{'})
/// ```
/// The inial condition is $`y(0) = 1.0`$ and we solve it in the interval
/// $`\lbrack 0, 2\rbrack`$.\ The following equation is the closed solution for
/// this ODE:
/// ```math
/// y(t) = 2 - e^{-t}
/// ```
///
/// ```
/// # #[macro_use]
/// # extern crate mathru;
/// # fn main()
/// # {
/// use mathru::{
///     algebra::linear::{Matrix, Vector},
///     analysis::differential_equation::ordinary::{ImplicitEuler, ImplicitODE},
/// };
///
/// pub struct ImplicitODEExample
/// {
///     time_span: (f64, f64),
///     init_cond: Vector<f64>,
/// }
///
/// impl Default for ImplicitODEExample
/// {
///     fn default() -> ImplicitODEExample
///     {
///         return ImplicitODEExample { time_span: (0.0, 2.0),
///                                     init_cond: vector![1.0] };
///     }
/// }
///
/// impl ImplicitODE<f64> for ImplicitODEExample
/// {
///     fn func(self: &Self, _t: f64, x: &Vector<f64>) -> Vector<f64>
///     {
///         let result = (x * &-4.0) + 8.0;
///         return result;
///     }
///
///     fn time_span(self: &Self) -> (f64, f64)
///     {
///         return self.time_span;
///     }
///
///     fn init_cond(self: &Self) -> Vector<f64>
///     {
///         return self.init_cond.clone();
///     }
///
///     fn jacobian(self: &Self, _t: f64, _input: &Vector<f64>) -> Matrix<f64>
///     {
///         let jacobian = matrix![-4.0];
///         return jacobian;
///     }
/// }
///
/// // We instanciate Euler's backward algorithm with a stepsize of 0.001
/// let step_size: f64 = 0.0001;
/// let solver: ImplicitEuler<f64> = ImplicitEuler::new(step_size);
///
/// let problem: ImplicitODEExample = ImplicitODEExample::default();
///
/// // Solve the ODE
/// let (t, x): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem).unwrap();
///
///
/// # }
/// ```
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ImplicitEuler<T>
{
    stepper: ImplicitFixedStepper<T>,

    root_finder: NewtonRaphson<T>,
}

impl<T> ImplicitEuler<T> where T: Real
{
    /// Creates a backward Euler instance
    pub fn new(step_size: T) -> ImplicitEuler<T>
    {
        return ImplicitEuler { stepper: ImplicitFixedStepper::new(step_size),
                               root_finder: NewtonRaphson::new(100, T::from_f64(0.00000001)) };
    }

    pub fn solve<F>(self: &Self, prob: &F) -> Result<(Vec<T>, Vec<Vector<T>>), ()>
        where F: ImplicitODE<T>
    {
        return self.stepper.solve(prob, self);
    }

    pub fn get_step_size(self: &Self) -> &T
    {
        return self.stepper.get_step_size();
    }

    pub fn set_step_size(self: &mut Self, step_size: T)
    {
        self.stepper.set_step_size(step_size)
    }
}

impl<T> ImplicitFixedStepSizeMethod<T> for ImplicitEuler<T> where T: Real
{
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ImplicitODE<T>
    {
        let t: T = *t_n + *h;
        let ie_helper = ImplicitEulerHelper::new(prob, &t, x_n, h);
        let x_n = self.root_finder.find_root(&ie_helper, x_n).unwrap();

        return x_n;
    }

    /// Euler's method is a first order method
    fn order(self: &Self) -> u8
    {
        return 1;
    }
}

/// this structure is implemented, therewith it is possible to implement the
/// traits needed by NewtonRaphson without exposing this traits.
struct ImplicitEulerHelper<'a, T, F>
    where T: Real,
          F: ImplicitODE<T>
{
    function: &'a F,
    t: &'a T,
    x: &'a Vector<T>,
    h: &'a T,
}

impl<'a, T, F> ImplicitEulerHelper<'a, T, F>
    where T: Real,
          F: ImplicitODE<T>
{
    pub fn new(function: &'a F,
               t: &'a T,
               x: &'a Vector<T>,
               h: &'a T)
               -> ImplicitEulerHelper<'a, T, F>
    {
        return ImplicitEulerHelper { function, t, x, h };
    }
}

impl<'a, T, F> Function<Vector<T>> for ImplicitEulerHelper<'a, T, F>
    where T: Real,
          F: ImplicitODE<T>
{
    type Codomain = Vector<T>;

    ///$`x(t_{n+1}) = y(t_n) + hf(t_{n+1}, x(t_{n+1})`$
    ///
    /// g(z) = y(t_n) + hf(t_{n+1}, z) - z)`$
    fn eval(self: &Self, z: &Vector<T>) -> Vector<T>
    {
        let result = &(self.x + &(&self.function.func(*self.t, z) * self.h)) - z;
        return result;
    }
}

impl<'a, T, F> Jacobian<T> for ImplicitEulerHelper<'a, T, F>
    where T: Real,
          F: ImplicitODE<T>
{
    /// $` \frac{\partial g(z)}{\partial z} = h \frac{\partial f(t_{n+1},
    /// z)}{\partial z} - I`$
    fn jacobian(self: &Self, z: &Vector<T>) -> Matrix<T>
    {
        let (m, _n): (usize, usize) = z.dim();
        let jacobian: Matrix<T> = self.function.jacobian(*self.t, z) * *self.h - Matrix::one(m);

        return jacobian;
    }
}
