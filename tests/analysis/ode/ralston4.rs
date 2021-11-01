use mathru::{
    algebra::linear::Vector,
    analysis::differential_equation::ordinary::{ExplicitODE, FixedStepper, Ralston4},
};

use super::problem::{ExplicitODE1, ExplicitODE2};

#[test]
fn fn1()
{
    let problem: ExplicitODE1 = ExplicitODE1::default();
    let solver: FixedStepper<f64> = FixedStepper::new(0.01);
    let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Ralston4::default()).unwrap();
    let len: usize = y.len();

    let time_span: (f64, f64) = problem.time_span();
    let init_cond: Vector<f64> = problem.init_cond();

    assert_relative_eq!(time_span.1, t[len - 1], epsilon=0.000000001);
    assert_relative_eq!(init_cond[0] * (2.0 * time_span.1).exp(), y[len - 1][0], epsilon=0.2);
}

#[test]
fn fn2()
{
    let problem: ExplicitODE2 = ExplicitODE2::default();
    let solver: FixedStepper<f64> = FixedStepper::new(0.1);
    let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Ralston4::default()).unwrap();

    let len: usize = y.len();

    let time_span: (f64, f64) = problem.time_span();

    assert_relative_eq!(time_span.1, t[len - 1], epsilon=0.00000001);
    assert_relative_eq!(time_span.1.tan(), y[len - 1][0], epsilon=0.1);
}
