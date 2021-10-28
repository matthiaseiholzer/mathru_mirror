//! Some experiments why mathru::algebra::linear::Vector is so much slower than std::Vec
//!

use mathru::{algebra::linear::Vector, vector};
use criterion::Criterion;
use mathru::analysis::differential_equation::ordinary::{ExplicitEuler, ExplicitODE, FixedStepper};



criterion_group!(euler, bench_vector_euler, bench_vec_euler);

fn vec_dgl(state: &[f64]) -> Vec<f64> {
    vec![state[0]]
}

pub fn vec_euler() -> f64 {
    let mut state = vec![1.];
    let dt: f64 = 0.01;
    for _ in 1..10000 {
        let diff = vec_dgl(&state);
        state[0] += diff[0] * dt;
    }
    state[0]
}

fn vector_dgl(state: &Vector<f64>) -> Vector<f64> {
    vector![state[0]]
}

pub fn vector_euler() -> f64 {
    let mut state = vector![1.];
    let dt: f64 = 0.01;
    for _ in 1..10000 {
        let diff = vector_dgl(&state);
        state[0] += diff[0] * dt;
    }
    state[0]
}

fn bench_vec_euler(bench: &mut Criterion) {
    bench.bench_function("bench_vec_euler", move |bh| {
        bh.iter(vec_euler);
    });
}

fn bench_vector_euler(bench: &mut Criterion) {
    bench.bench_function("bench_vector_euler", move |bh| {
        bh.iter(vector_euler);
    });
}

#[inline]
pub fn dgl_x(_x: f64, v: f64) -> f64 {
    v
}
#[inline]
pub fn dgl_v(x: f64, _v: f64) -> f64 {
    -x
}

pub struct LibDgl {
    time_span: (f64, f64),
    init_cond: Vector<f64>,
}

impl ExplicitODE<f64> for LibDgl {
    #[inline]
    fn func(&self, _t: &f64, state: &Vector<f64>) -> Vector<f64> {
        let x = state[0];
        let v = state[1];
        vector![dgl_x(x,v); dgl_v(x,v)]
    }
    #[inline]
    fn time_span(&self) -> (f64, f64) {
        self.time_span
    }
    #[inline]
    fn init_cond(&self) -> Vector<f64> {
        self.init_cond.clone()
    }
}

pub fn lib_euler() -> (f64, f64) {
    let solver: FixedStepper<f64> = FixedStepper::new(0.01);
    let problem = LibDgl {
        time_span: (0.0, 100.),
        init_cond: vector![1.;0.],
    };
    let (_t, y): (Vec<f64>, Vec<Vector<f64>>) =
        solver.solve(&problem, &ExplicitEuler::default()).unwrap();
    (y[y.len() - 1][0], y[y.len() - 1][1])
}

fn bench_lib(bench: &mut Criterion) {
    bench.bench_function("bench_lib", move |bh| {
        bh.iter(lib_euler);
    });
}


