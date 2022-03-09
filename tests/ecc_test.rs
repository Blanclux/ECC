//!
//! Elliptic Curve Calculation Test Program
//!
//  Wrtten by blanclux
//  This software is distributed on an "AS IS" basis WITHOUT WARRANTY OF ANY KIND.
extern crate ecc;

use ibig::{ibig, IBig, UBig};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::{env, str};
use std::time::SystemTime;

use ecc::ec_param::{EcParam, ParamOp};
use ecc::ecc_a::EcpA;
use ecc::ecc_j::EcpJ;
use ecc::{EcAxis, EcOp, Point};

#[test]
fn ecc_test() {
    println!("< EC Test >");
    let args: Vec<String> = env::args().collect();
    let ecname = if args.len() == 2 {
        &args[1]
    } else {
        "secp160k1"
    };

    let param = EcParam::new(ecname);
    println!(" EC Name: {}", param.get_ecid());

    println!("\n< ibig version >");
    test_ibig(ecname, &param);

    println!("\n< bigint version >");
    test_bigint(ecname, &param);
}

fn test_ibig(ecname: &str, param: &EcParam) {
    let eca = EcpA::new(ecname);

    let a: IBig = param.get_a();
    println!(" a: {:X}", &a);
    let b: IBig = param.get_b();
    println!(" b: {:X}", &b);
    let p: IBig = param.get_prime();
    println!(" p: {:X}", &p);
    let n: IBig = param.get_order();
    println!(" n: {:X}", &n);
    let h: IBig = param.get_cofactor();
    println!(" h: {:02X}", &h);
    let gx: IBig = param.get_gx();
    println!(" gx: {:X}", &gx);
    let gy: IBig = param.get_gy();
    println!(" gy: {:X}", &gy);

    println!("\n< 2D (Affine) >");

    let p_0 = Point {
        axis: EcAxis::Affine,
        x: ibig!(0),
        y: ibig!(0),
        z: ibig!(1),
    };
    assert!(p_0.is_zero());
    assert!(eca.on_curve(&p_0));

    let p_g = Point {
        axis: EcAxis::Affine,
        x: gx,
        y: gy,
        z: ibig!(1),
    };
    assert!(eca.on_curve(&eca.p_g));
    assert!(eca.on_curve(&p_g));
    println!(" G: {}", &p_g);

    let g = get_g_ib(param);
    let _g_p = eca.decode_point(&g);
    println!("> decode_point");
    println!(" G: {}", &_g_p);

    do_test(&eca);
    do_perform(&eca, &n);

    println!("\n< 3D (Jacobian) >");
    let gx = param.get_gx();
    let gy = param.get_gy();

    let ecj = EcpJ::new(ecname);

    let p_0 = Point {
        axis: EcAxis::Proj,
        x: ibig!(0),
        y: ibig!(0),
        z: ibig!(0),
    };
    assert!(p_0.is_zero());
    assert!(ecj.on_curve(&p_0));

    let p_g = Point {
        axis: EcAxis::Proj,
        x: gx,
        y: gy,
        z: ibig!(1),
    };
    assert!(ecj.on_curve(&ecj.p_g));
    assert!(ecj.on_curve(&p_g));
    println!(" G: {}", &p_g);

    let g = get_g_ib(param);
    let _g_p = ecj.decode_point(&g);
    println!("> decode_point");
    println!(" G: {}", &_g_p);

    println!("> normalize");
    let p_r = ecj.gen_point();
    assert!(ecj.on_curve(&p_r));
    let mut p_a = ecj.add(&p_r, &p_r);
    println!(" R: {}", &p_a);
    ecj.normalize(&mut p_a);
    println!(" R: {}", &p_a);

    println!("> to_affine");
    let mut p_r = ecj.gen_point();
    p_r = ecj.add(&p_r, &p_r);
    println!(" R: {}", &p_r);
    let p_a = ecj.to_affine(&p_r);
    println!(" R: {}", &p_a);

    do_test(&ecj);
    do_perform(&ecj, &n);
}

fn get_g_ib(param: &EcParam) -> Vec<u8> {
    let s = str::from_utf8(param.g).ok().unwrap();
    let g = UBig::from_str_radix(s, 16).unwrap();

    g.to_be_bytes()
}

fn test_bigint(ecname: &str, param: &EcParam) {
    let eca = EcpA::new(ecname);

    let a: BigInt = param.get_a();
    println!(" a: {:X}", &a);
    let b: BigInt = param.get_b();
    println!(" b: {:X}", &b);
    let p: BigInt = param.get_prime();
    println!(" p: {:X}", &p);
    let n: BigInt = param.get_order();
    println!(" n: {:X}", &n);
    let h: BigInt = param.get_cofactor();
    println!(" h: {:02X}", &h);
    let gx: BigInt = param.get_gx();
    println!(" gx: {:X}", &gx);
    let gy: BigInt = param.get_gy();
    println!(" gy: {:X}", &gy);

    println!("\n< 2D (Affine) >");

    let p_0 = Point {
        axis: EcAxis::Affine,
        x: BigInt::zero(),
        y: BigInt::zero(),
        z: BigInt::zero(),
    };
    assert!(p_0.is_zero());
    assert!(eca.on_curve(&p_0));

    let p_g = Point {
        axis: EcAxis::Affine,
        x: gx,
        y: gy,
        z: BigInt::one(),
    };
    assert!(eca.on_curve(&eca.p_g));
    assert!(eca.on_curve(&p_g));
    println!(" G: {}", &p_g);

    let g = get_g_bi(param);
    let _g_p = eca.decode_point(&g);
    println!("> decode_point");
    println!(" G: {}", &_g_p);

    do_test_bi(&eca);
    do_perform_bi(&eca, &n);

    println!("\n< 3D (Jacobian) >");
    let gx = param.get_gx();
    let gy = param.get_gy();

    let ecj = EcpJ::new(ecname);

    let p_0 = Point {
        axis: EcAxis::Proj,
        x: BigInt::zero(),
        y: BigInt::zero(),
        z: BigInt::zero(),
    };
    assert!(p_0.is_zero());
    assert!(ecj.on_curve(&p_0));

    let p_g = Point {
        axis: EcAxis::Proj,
        x: gx,
        y: gy,
        z: BigInt::one(),
    };
    assert!(ecj.on_curve(&ecj.p_g));
    assert!(ecj.on_curve(&p_g));
    println!(" G: {}", &p_g);

    let g = get_g_bi(param);
    let _g_p = ecj.decode_point(&g);
    println!("> decode_point");
    println!(" G: {}", &_g_p);

    println!("> normalize");
    let p_r = ecj.gen_point();
    assert!(ecj.on_curve(&p_r));
    let mut p_a = ecj.add(&p_r, &p_r);
    println!(" R: {}", &p_a);
    ecj.normalize(&mut p_a);
    println!(" R: {}", &p_a);

    println!("> to_affine");
    let mut p_r = ecj.gen_point();
    p_r = ecj.add(&p_r, &p_r);
    println!(" R: {}", &p_r);
    let p_a = ecj.to_affine(&p_r);
    println!(" R: {}", &p_a);

    do_test_bi(&ecj);
    do_perform_bi(&ecj, &n);
}

fn do_test(ecp: &impl EcOp<IBig>) {
    let p_0 = ecp.get_zero();
    assert!(ecp.is_zero(&p_0));
    assert!(ecp.on_curve(&p_0));

    let p_g = ecp.get_gp();
    assert!(ecp.on_curve(&p_g));

    // P - P = O
    println!("> P - P = O");
    let mut p_q = ecp.negate(&p_g);
    let mut p_r = ecp.add(&p_g, &p_q);
    assert!(ecp.is_zero(&p_r));
    // Q = G + G
    println!("> Add: G + G");
    p_q = ecp.add(&p_g, &p_g);
    assert!(ecp.on_curve(&p_q));
    // Q = 2 G (double)
    println!("> Double: 2 * G");
    p_q = ecp.double(&p_g);
    assert!(ecp.on_curve(&p_q));
    // Q = 2 G (mul)
    println!("> Mul: 2 * G");
    p_q = ecp.mul_bin(&p_g, &ibig!(2));
    assert!(ecp.on_curve(&p_q));
    // Q = 2G + G
    p_q = ecp.add(&p_q, &p_g);
    assert!(ecp.on_curve(&p_q));
    // R = 3 G
    p_r = ecp.mul(&p_g, &ibig!(3));
    assert!(ecp.on_curve(&p_r));
    println!("> 2 * G + G = 3 * G");
    assert!(ecp.equals(&p_q, &p_r));

    println!("> n * G = O (mul_bin)");
    p_q = ecp.mul_bin(&p_g, &ecp.get_order());
    //ecp.print(" n * G", &p_q);
    assert!(ecp.on_curve(&p_q));
    assert!(ecp.is_zero(&p_q));

    println!("> n * G = O (mul)");
    p_q = ecp.mul(&p_g, &ecp.get_order());
    //ecp.print(" n * G", &p_q);
    assert!(ecp.on_curve(&p_q));
    assert!(ecp.is_zero(&p_q));
}

fn get_g_bi(param: &EcParam) -> Vec<u8> {
    let g = BigInt::parse_bytes(param.g, 16).unwrap();
    g.to_bytes_be().1
}

/// Perormance test
fn do_perform(ecp: &impl EcOp<IBig>, k: &IBig) {
    println!("> EC mul performance");
    println!(" k = {}", k);
    let count = 25;

    let mut p_q = ecp.get_gp();
    let p_g = ecp.get_gp();

    // Binary method
    let s_time = SystemTime::now();
    for _ in 1..count {
        p_q = ecp.mul_bin(&p_g, k);
    }
    let e_time = s_time.elapsed().expect("Clock may have gone backwards");
    assert!(ecp.on_curve(&p_q));
    println!(" time: {:?}  (bin)", e_time / count);

    // Windows method
    let s_time = SystemTime::now();
    for _ in 1..count {
        p_q = ecp.mul(&p_g, k);
    }
    let e_time = s_time.elapsed().expect("Clock may have gone backwards");
    assert!(ecp.on_curve(&p_q));
    println!(" time: {:?}  (win)", e_time / count);
}

fn do_test_bi(ecp: &impl EcOp<BigInt>) {
    let p_0 = ecp.get_zero();
    assert!(ecp.is_zero(&p_0));
    assert!(ecp.on_curve(&p_0));

    let p_g = ecp.get_gp();
    assert!(ecp.on_curve(&p_g));

    // P - P = O
    println!("> P - P = O");
    let mut p_q = ecp.negate(&p_g);
    let mut p_r = ecp.add(&p_g, &p_q);
    assert!(ecp.is_zero(&p_r));
    // Q = G + G
    println!("> Add: G + G");
    p_q = ecp.add(&p_g, &p_g);
    assert!(ecp.on_curve(&p_q));
    // Q = 2 G (double)
    println!("> Double: 2 * G");
    p_q = ecp.double(&p_g);
    assert!(ecp.on_curve(&p_q));
    // Q = 2 G (mul)
    println!("> Mul: 2 * G");
    p_q = ecp.mul_bin(&p_g, &BigInt::from(2));
    assert!(ecp.on_curve(&p_q));
    // Q = 2G + G
    p_q = ecp.add(&p_q, &p_g);
    assert!(ecp.on_curve(&p_q));
    // R = 3 G
    p_r = ecp.mul(&p_g, &BigInt::from(3));
    assert!(ecp.on_curve(&p_r));
    println!("> 2 * G + G = 3 * G");
    assert!(ecp.equals(&p_q, &p_r));

    println!("> n * G = O (mul_bin)");
    p_q = ecp.mul_bin(&p_g, &ecp.get_order());
    //ecp.print(" n * G", &p_q);
    assert!(ecp.on_curve(&p_q));
    assert!(ecp.is_zero(&p_q));

    println!("> n * G = O (mul)");
    p_q = ecp.mul(&p_g, &ecp.get_order());
    //ecp.print(" n * G", &p_q);
    assert!(ecp.on_curve(&p_q));
    assert!(ecp.is_zero(&p_q));
}

/// Perormance test
fn do_perform_bi(ecp: &impl EcOp<BigInt>, k: &BigInt) {
    println!("> EC mul performance");
    println!(" k = {}", k);
    let count = 25;

    let mut p_q = ecp.get_gp();
    let p_g = ecp.get_gp();

    // Binary method
    let s_time = SystemTime::now();
    for _ in 1..count {
        p_q = ecp.mul_bin(&p_g, k);
    }
    let e_time = s_time.elapsed().expect("Clock may have gone backwards");
    assert!(ecp.on_curve(&p_q));
    println!(" time: {:?}  (bin)", e_time / count);

    // Windows method
    let s_time = SystemTime::now();
    for _ in 1..count {
        p_q = ecp.mul(&p_g, k);
    }
    let e_time = s_time.elapsed().expect("Clock may have gone backwards");
    assert!(ecp.on_curve(&p_q));
    println!(" time: {:?}  (win)", e_time / count);
}
