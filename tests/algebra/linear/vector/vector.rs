
#[cfg(test)]
mod vector_test
{
    use mathru::algebra::linear::{Vector, Matrix};
    use mathru::elementary::Power;


    #[test]
    fn macro_vector_column()
    {
        let vec: Vector<f32> = vector![1.0; 2.0; 3.0];

        let vec_ref: Vector<f32> = Vector::new_column(3, vec![1.0, 2.0, 3.0]);

        assert_eq!(vec_ref, vec);
    }

    #[test]
    fn macro_vector_row()
    {
        let vec: Vector<f32> = vector![1.0, 2.0, 3.0];

        let vec_ref: Vector<f32> = Vector::new_row(3, vec![1.0, 2.0, 3.0]);

        assert_eq!(vec_ref, vec);
    }

    #[test]
    fn serde_0()
    {
        let mat:Vector<f64> = vector![1.0; 2.0; 3.0];
        let serialized = serde_json::to_string(&mat).unwrap();

        let mat_s: Vector<f64> = serde_json::from_str(&serialized).unwrap();

        assert_eq!(mat_s, mat);
    }

    #[test]
    fn zero()
    {
        let rows: usize = 5;

        let m_zero : Vector<f32> = Vector::zero(rows);
        let (m, n) = m_zero.dim();
        assert_eq!(m, rows);
        assert_eq!(n, 1);

        for i in 0..rows
        {
            assert_eq!(*(m_zero.get(&i)), 0.0);
        }
    }

    #[test]
    fn partial_eq0()
    {
        let lhs: Vector<f32> = Vector::new_column(2, vec![1.0, 2.0]);
        let rhs: Vector<f32> = Vector::new_column(2, vec![1.0, 2.0]);
        assert_eq!(lhs, rhs);
    }

     #[test]
    fn partial_eq1()
    {
        let lhs: Vector<u32> = Vector::new_column(2, vec![1, 2]);
        let rhs: Vector<u32> = Vector::new_column(2, vec![1, 2]);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn add_owner()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let b : Vector<f32> = Vector::new_column(dim, vec![1.0, 4.0, -1.0, 0.0, -7.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![2.0, 6.0, 2.0, 4.0, -2.0]);

        let res : Vector<f32> = a + b;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn add_borrow()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let b : Vector<f32> = Vector::new_column(dim, vec![1.0, 4.0, -1.0, 0.0, -7.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![2.0, 6.0, 2.0, 4.0, -2.0]);

        let res : Vector<f32> = &a + &b;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn scalar_add_owner()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![6.0, 7.0, 8.0, 9.0, 10.0]);

        let res : Vector<f32> = a + 5.0;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn scalar_add_borrow()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![6.0, 7.0, 8.0, 9.0, 10.0]);

        let res : Vector<f32> = &a + &5.0;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn sub()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let b : Vector<f32> = Vector::new_column(dim, vec![1.0, 4.0, -1.0, 0.0, -7.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![0.0, -2.0, 4.0, 4.0, 12.0]);

        let res : Vector<f32> = a - b;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn sub_ref()
    {
        let dim: usize = 5;
        let a : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let b : Vector<f32> = Vector::new_column(dim, vec![1.0, 4.0, -1.0, 0.0, -7.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![0.0, -2.0, 4.0, 4.0, 12.0]);

        let res : Vector<f32> = &a - &b;

        for i in 0..dim
        {
            assert_eq!(*(res.get(&i)), *(res_ref.get(&i)));
        }
    }

    #[test]
    fn get_0()
    {
        let dim: usize = 5;
        let res : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        for i in 0..dim
        {
            assert_eq!(*(res_ref.get(&i)), *(res.get(&i)));
        }
    }

    #[test]
    fn get_1()
    {
        let dim: usize = 5;
        let res : Vector<f32>  = Vector::new_row(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let res_ref : Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        for i in 0..dim
        {
            assert_eq!(*(res_ref.get(&i)), *(res.get(&i)));
        }
    }

    #[test]
    fn get_slice_0()
    {
        let dim: usize = 5;
        let res: Vector<f32>  = Vector::new_column(dim, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let res_ref: Vector<f32>  = Vector::new_column(3, vec![3.0, 4.0, 5.0]);

        let slice: Vector<f32> = res.get_slice(2, 4);

        assert_eq!(res_ref, slice);
    }

    #[test]
    fn set_slice_0()
    {
        let mut a: Vector<f32>  = Vector::new_column(5, vec![1.0, 2.0, 3.0, 4.0, 5.0]);

        let b: Vector<f32>  = Vector::new_column(3, vec![ -3.0, -4.0, -5.0]);
        let res_ref: Vector<f32>  = Vector::new_column(5, vec![1.0, 2.0, -3.0, -4.0, -5.0]);

        a.set_slice(&b, 2);

        assert_eq!(res_ref, a);
    }

    #[test]
    fn transpose()
    {
        let (m_ref, n_ref) : (usize, usize) = (4, 1);
        let  vec : Vector<f32> = Vector::new_column(m_ref, vec![2.0, 6.0, -2.5, 0.0]);
        let vec_trans : Vector<f32>= vec.transpose();

        let (m, n) : (usize, usize) =  vec_trans.dim();

        assert_eq!((n_ref, m_ref), (m,n));
    }


    #[test]
    fn dotp_0()
    {
        let a: Vector<f32> = Vector::new_column(4, vec![-1.0, -3.0, 6.0, -1.0]);
        let b: Vector<f32> = Vector::new_column(4, vec![-2.0, -5.0, -3.0, 2.0]);
        let dotp_ref : f32 = -3.0;
        let dotp : f32 = a.dotp(&b);
        assert_eq!(dotp_ref, dotp);
    }

    #[test]
    fn dotp_1()
    {
        let a: Vector<f32> = Vector::new_column(4, vec![-1.0, -3.0, 6.0, -1.0]);
        let dotp_ref : f32 = 47.0;
        let dotp : f32 = a.dotp(&a);
        assert_eq!(dotp_ref, dotp);
    }

//    #[test]
//    fn crossp()
//    {
//        let a: Vector<f32> = Vector::new_column(&3, vec![-1.0), -3.0),
//                                                                              6.0)]);
//        let b: Vector<f32> = Vector::new_column(&3, vec![-2.0), -5.0),
//                                                                              Real::new
//            (-3.0)]);
//        let crossp_ref : f32 = -3.0);
//        let cross : Vector<f32> = a.crossp(&b);
//        assert_eq!(crossp_ref, croossp);
//    }

    #[test]
    fn dyadp()
    {
        let x: Vector<f32> = Vector::new_column(3, vec![1.0, 3.0, 2.0]);
        let y: Vector<f32> = Vector::new_column(4, vec![2.0, 1.0, 0.0, 3.0]);
        let dyadp_ref : Matrix<f32> = Matrix::new(3, 4, vec![2.0, 1.0, 0.0, 3.0, 6.0, 3.0, 0.0, 9.0, 4.0, 2.0, 0.0, 6.0]);

        let p: Matrix<f32> =  x.dyadp(&y);

        assert_eq!(dyadp_ref, p);
    }

    #[test]
    fn p_norm()
    {
        let p : f32 = 2.0;
        let v : Vector<f32> = Vector::new_column(4, vec![-2.0, -5.0, -3.0, 2.0]);
        let p_norm_ref : f32 = 42.0.pow(&0.5);
        let p_norm : f32 =  v.p_norm(&p);
        assert_eq!(p_norm_ref, p_norm);
    }

    #[test]
    fn scalar_mul_owner()
    {
        let v = vector![1.0; 2.0; 3.0; 4.0];
        let prod_ref = vector![-0.5; -1.0; -1.5; -2.0];

        let res = v * -0.5;

        assert_eq!(prod_ref, res);
    }

    #[test]
    fn scalar_mul_borrow()
    {
        let v = vector![1.0; 2.0; 3.0; 4.0];
        let prod_ref = vector![-0.5; -1.0; -1.5; -2.0];

        let res = &v * &-0.5;

        assert_eq!(prod_ref, res);
    }

    #[test]
    fn matrix_mul_owner()
    {
        let m = matrix![1.0, 2.0; 3.0, 4.0];
        let v = vector![1.0, 2.0];
        let prod_ref = vector![7.0, 10.0];

        let res = v * m;

        assert_eq!(prod_ref, res);
    }

    #[test]
    fn matrix_mul_borrow()
    {
        let m = matrix![1.0, 2.0; 3.0, 4.0];
        let v = vector![1.0, 2.0];
        let prod_ref = vector![7.0, 10.0];

        let res = &v * &m;

        assert_eq!(prod_ref, res);
    }

    #[test]
    fn argmax()
    {
        let m = vector![1.0, -2.0, 3.0, -4.0];

        assert_eq!(2, m.argmax());
    }

    #[test]
    fn argmin()
    {
        let m = vector![1.0, -2.0, 3.0, -4.0];

        assert_eq!(3, m.argmin());
    }
}