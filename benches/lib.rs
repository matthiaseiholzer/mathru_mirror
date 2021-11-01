//#![feature(test)]
#![allow(unused_macros)]
extern crate mathru;

#[macro_use]
extern crate criterion;


mod algebra;
mod analysis;


criterion_main!(
    analysis::vector_bench::euler,
    algebra::linear::vector::add::vector_add,
    algebra::linear::vector::sub::vector_sub ,
    algebra::linear::vector::mul::vector_mul*/);
// criterion_main!(algebra::linear::matrix::matrix, analysis::ode::ode,);
