use mathru::algebra::linear::Vector;
use std::ops::Add;

#[test]
fn add_owner()
{
    let a: Vector<f32> = Vector::new_column(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let b: Vector<f32> = Vector::new_column(vec![1.0, 4.0, -1.0, 0.0, -7.0]);
    let res_ref: Vector<f32> = Vector::new_column(vec![2.0, 6.0, 2.0, 4.0, -2.0]);

    let res: Vector<f32> = a + b;

    assert_relative_eq!(res_ref, res);
}

#[test]
fn add_borrow()
{
    let a: Vector<f32> = Vector::new_column(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let b: Vector<f32> = Vector::new_column(vec![1.0, 4.0, -1.0, 0.0, -7.0]);
    let res_ref: Vector<f32> = Vector::new_column(vec![2.0, 6.0, 2.0, 4.0, -2.0]);

    let res: Vector<f32> = &a + &b;

    assert_relative_eq!(res, res_ref);
}

#[test]
fn scalar_add_owner()
{
    let a: Vector<f32> = Vector::new_column(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let res_ref: Vector<f32> = Vector::new_column(vec![6.0, 7.0, 8.0, 9.0, 10.0]);

    let res: Vector<f32> = a.add(5.0);

    assert_relative_eq!(res, res_ref);
}

// #[test]
// fn scalar_add_borrow()
// {
//     let mut a: Vector<f32> = Vector::new_column(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
//     let res_ref: Vector<f32> = Vector::new_column(vec![6.0, 7.0, 8.0, 9.0, 10.0]);
//
//     &a + &5.0f32;
//
//     // assert_relative_eq!(*res, res_ref);
// }

#[test]
fn add_assign_matrix()
{
    let mut a: Vector<f32> = vector![1.0, -2.0, -3.0];

    let b: Vector<f32> = vector![6.0, 3.0, 2.0];

    let sum_ref = vector![7.0, 1.0, -1.0];
    a += b;

    assert_relative_eq!(sum_ref, a);
}

#[test]
fn add_assign_scalar()
{
    let mut a: Vector<f32> = vector![ 1.0, -2.0, -3.0];

    a += 5.0f32;

    let sum: Vector<f32> = vector![6.0, 3.0, 2.0];

    assert_relative_eq!(sum, a);
}