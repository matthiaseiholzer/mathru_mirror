use crate::algebra::linear::{Vector, Matrix};
use crate::algebra::linear::matrix::Solve;
use crate::optimization::{OptimResult, Jacobian, Hessian};
use std::marker::PhantomData;
use crate::algebra::abstr::Real;


/// Newton's method
///
/// ```math
/// f \colon \mathbb{R}^n \to \mathbb{R}
/// ```
/// is a twice differentiable function.
///
/// Newton's method solves the minimization problem
///
/// ```math
///  f(x) \to min
/// ```
///
/// input: $` f \colon \mathbb{R}^{n} \to \mathbb{R} `$ with initial approximation $` x_{0} \in \mathbb{R}^{n} `$
///
/// output: $` x_{k} `$
///
/// 1. Initialization: Choose $` \sigma \in (0, 1) `$
///
///     $` \rho > 0, k := 0 `$
/// 2. Solve de equation system $` \nabla^{2} f(x_{k})d_{k} = -\nabla f(x_{k}) `$
/// 3. If the euqation is not solvable, or the condition
/// $` \nabla f(x_{k})^{T}d_{k} \leq -\rho \lvert \lvert \nabla f(x_k) \rvert \rvert_{2}^{2} `$  is not fullfilled
///
///     Than $` d_{k} := \nabla f(x_{k}) `$
/// 4. $` \alpha_{k} := 1 `$
/// 5. while  $` f(x_{k} + \alpha_{k}d_{k}) > f(x_{k}) + \sigma \alpha_{k} \nabla f(x_{k})^{T}d_{k} `$
///
///     set $` \alpha_{k} `$
/// 6. $` x_{k + 1} := x_{k} + d_{k} `$
/// 7. $` k := k + 1 `$  go to 2.
///
/// ```
pub struct Newton<T>
{
    iters: u64,
    sigma: T,
    rho: T,
    __phantom: PhantomData<T>
}

impl<T> Newton<T>
{
    /// Creates an instance of newtons method
    ///
    /// # Arguments
    ///
    /// * 'iters': Number of iterations
    pub fn new(iters: u64, sigma: T, rho: T) -> Newton<T>
    {
        Newton
        {
            iters: iters,
            sigma: sigma,
            rho: rho,
            __phantom: PhantomData
        }
    }

}

impl<T> Newton<T>
    where T: Real
{
    /// Minimize function func
    ///
    /// # Arguments
    ///
    /// * 'func0': Function to be minimized
    /// * 'x_0': Initial guess for the minimum
    ///
    /// # Return
    ///
    /// local minimum
    pub fn minimize<F: Jacobian<T> + Hessian<T>>(self: &Self, func: &F, x_0: &Vector<T>) -> OptimResult<Vector<T>>
    {
        let mut x_n: Vector<T> = x_0.clone();

        for _i in 0..self.iters
        {
            let hessian_x_n: Matrix<T> = func.hessian(&x_n);
            let grad_x_n: Vector<T> = func.jacobian(&x_n).get_row(0).transpose();
            let res_solve: Option<Vector<T>> = hessian_x_n.solve(&-grad_x_n.clone());
            let d_k: Vector<T>;

            match res_solve
            {
                Some(d_k_temp) =>
                {
                    let grad_x_n_abs: T = grad_x_n.dotp(&grad_x_n);
                    let grad_d_k_temp: T = grad_x_n.dotp(&d_k_temp);
                    if grad_d_k_temp <= -self.rho * grad_x_n_abs
                    {
                        println!("IN");
                        d_k = d_k_temp;
                    }
                    else
                    {
                        d_k = -grad_x_n.clone();
                    }
                }
                None =>
                {
                    d_k = -grad_x_n.clone();
                }
            }
            let mut alpha: T = T::one();
            //Backtracking line search
            //Armijo–Goldstein condition
            let mut r: Vector<T> = &x_n + &(&d_k * &alpha);
            println!("d_k: {} x_n: {}, r: {}", d_k, x_n, r);
            let mut f_r: T = *func.eval(&r).get(0);
            let f_x_n: T = *func.eval(&x_n).get(0);
            let temp: T = grad_x_n.dotp(&d_k);
            println!("{}", d_k);
            while f_r > f_x_n + temp * (self.sigma * alpha)
            {
                alpha = alpha / T::from_f64(2.0).unwrap();
                r = &x_n + &(&d_k * &alpha);
                f_r = *func.eval(&r).get(0);
            }

            //Make step
            println!("{}, {}", d_k, alpha);
            x_n = x_n + d_k * alpha;

        }

        return OptimResult::new(x_n);
    }

}