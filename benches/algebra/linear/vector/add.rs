use mathru::{algebra::linear::Vector};
use std::ops::Add;
use criterion::Criterion;

criterion_group!(vector_add, bench_vector_add_assign_scalar, /*bench_vec_addAssign_scalar, */bench_vector_add_scalar, bench_vec_add_scalar, bench_vector_add_vector, bench_vec_add_vec);

fn bench_vector_add_assign_scalar(bench: &mut Criterion)
{
    bench.bench_function("bench_vector_add_assign_scalar", move |bh| {
        bh.iter(vector_add_assign_scalar);
    });
}

fn vector_add_assign_scalar()
{
    let mut vec: Vector<f64> = Vector::new_column(vec![3.0; 100000]);

    vec += 3.0f64;
}

// fn bench_vec_add_assign_scalar(bench: &mut Criterion)
// {
//     bench.bench_function("bench_vec_add_assign_scalar", move |bh| {
//         bh.iter(vec_add_assign_scalar);
//     });
// }
//
// fn vec_add_assign_scalar()
// {
//     let mut vec: Vec<f64> = vec![3.0; 100000];
//
//     vec += 3.0f64;
// }

fn bench_vector_add_scalar(bench: &mut Criterion)
{
    bench.bench_function("bench_vector_add_scalar", move |bh| {
        bh.iter(vector_add_scalar);
    });
}

fn vector_add_scalar()
{
    let mut vec: Vector<f64> = Vector::new_column(vec![3.0; 100000]);

    let _ = (&mut vec) + &3.0f64;
}

fn bench_vec_add_scalar(bench: &mut Criterion)
{
    bench.bench_function("bench_vec_add_scalar", move |bh| {
        bh.iter(vec_add_scalar);
    });
}

fn vec_add_scalar()
{
    let mut vec: Vec<f64> = vec![3.0; 100000];

    vec.iter_mut().for_each(|a: &mut f64| *a += 3.0);
}

fn bench_vector_add_vector(bench: &mut Criterion)
{
    bench.bench_function("bench_vector_add_vector", move |bh| {
        bh.iter(vector_add_vector);
    });
}

fn vector_add_vector()
{
    let mut vec1: Vector<f64> = Vector::new_column(vec![3.0; 100000]);
    let vec2: Vector<f64> = Vector::new_column(vec![3.0; 100000]);

    let _: &mut Vector<f64> = (&mut vec1).add(&vec2);
}


fn bench_vec_add_vec(bench: &mut Criterion)
{
    bench.bench_function("bench_vec_add_vec", move |bh| {
        bh.iter(vec_add_vec);
    });
}

// fn vec_add_vec()
// {
//     let mut vec1: Vec<f64> = vec![3.0; 100000];
//     let vec2: Vec<f64> = vec![3.0; 100000];
//
//     let _: Vec<f64>= vec1.iter_mut().zip(&vec2).map(|(a, b)| *a + *b).collect();
// }

fn vec_add_vec()
{
    let mut vec1: Vec<f64> = vec![3.0; 100000];
    let vec2: Vec<f64> = vec![3.0; 100000];

    // for i in 0..vec1.len() {
    //     vec1[i] += vec2[i];
    // }
    let _: Vec<f64>= vec1.iter_mut().zip(&vec2).map(|(a, b)| *a += *b).collect();
}
