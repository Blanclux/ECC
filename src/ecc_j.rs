//!
//! Elliptic Curve Calculation (Projective - Jacobian)
//!
//  Written by blanclux
//  This software is distributed on an "AS IS" basis WITHOUT WARRANTY OF ANY KIND.
use super::ec_param::{EcParam, ParamOp};
use super::number::Number;
use super::{EcAxis, EcOp, Point};

/// Elliptic curve  (Projective - Jacobian)
#[derive(Clone)]
pub struct EcpJ<T: Number> {
    pub a: T,
    pub b: T,
    pub p: T,
    pub n: T,
    pub h: T,
    pub p_zero: Point<T>,
    pub p_g: Point<T>,
}

impl<T: Number> EcOp<T> for EcpJ<T> {
    /// Get zero point
    fn get_zero(&self) -> Point<T> {
        Point {
            axis: EcAxis::Proj,
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }
    /// Get zero point
    #[inline]
    fn is_zero(&self, p_p: &Point<T>) -> bool {
        p_p.is_zero()
    }

    /// Copy point
    fn set(&self, p_p1: &mut Point<T>, p_p2: &Point<T>) {
        *p_p1 = p_p2.clone();
    }

    /// Point P is on curve ?
    #[allow(clippy::many_single_char_names)]
    fn on_curve(&self, p_p: &Point<T>) -> bool {
        let (l, r);
        let (z2, z4, z6);

        if p_p.z.is_one() {
            return true;
        }
        // Y^2 = X^3 + aXZ^2 + bZ^6
        l = p_p.y.mul_ref(&p_p.y) % &self.p;
        z2 = p_p.z.mul_ref(&p_p.z) % &self.p;
        z4 = z2.mul_ref(&z2) % &self.p;
        z6 = z2.mul_ref(&z4) % &self.p;
        r = ((p_p.x.mul_ref(&p_p.x) * &p_p.x) + self.a.mul_ref(&p_p.x) * &z4 + self.b.mul_ref(&z6))
            % &self.p;

        l == r
    }

    /// Convert 3D to 2D
    fn to_affine(&self, p_p: &Point<T>) -> Point<T> {
        if p_p.z.is_one() {
            return p_p.clone();
        }

        // x = x * z^(-2)
        let x = (Number::mod_inv(&(p_p.z.mul_ref(&p_p.z)), &self.p) * &p_p.x) % &self.p;
        // y = y * z^(-3)
        let y = (Number::mod_inv(&Number::mod_pow(&p_p.z, &T::from(3u32), &self.p), &self.p)
            * &p_p.y)
            % &self.p;

        Point {
            axis: EcAxis::Proj,
            x,
            y,
            z: T::one(),
        }
    }

    fn normalize(&self, p_p: &mut Point<T>) {
        if p_p.z.is_one() || p_p.z.is_zero() {
            return;
        }

        // x = x * z^(-2)
        p_p.x = (Number::mod_inv(&(p_p.z.mul_ref(&p_p.z)), &self.p) * &p_p.x) % &self.p;
        // y = y * z^(-3)
        p_p.y = (Number::mod_inv(&Number::mod_pow(&p_p.z, &T::from(3u32), &self.p), &self.p)
            * &p_p.y)
            % &self.p;
        p_p.z = T::one();
    }

    /// Check equality
    fn equals(&self, p1: &Point<T>, p2: &Point<T>) -> bool {
        let mut result = true;
        let (mut z1, mut z2);

        if p1.z.is_one() && p2.z.is_one() {
            return p1.x == p2.x && p1.y == p2.y;
        }
        if p1.z.is_one() {
            // z2^2
            z2 = p2.z.mul_ref(&p2.z) % &self.p;
            // x1 * z2^2 = x2 ?
            result = result && (p2.x == p1.x.mul_ref(&z2) % &self.p);
            z2 = z2.mul_ref(&p1.z) % &self.p;
            // y1 * z2^3 = y1 ?
            result = result && (p1.y == p1.y.mul_ref(&z2) % &self.p);
            return result;
        }
        z1 = p1.z.mul_ref(&p1.z) % &self.p; // z1^2
        if p2.z.is_one() {
            // P1.x * z^2 = x ?
            result = result && (p1.x == p2.x.mul_ref(&z1) % &self.p);
            z1 = z1.mul_ref(&p1.z);
            // P1.y * z^3 = y ?
            result = result && (p1.y == p2.y.mul_ref(&z1) % &self.p);
        } else {
            // z2^2
            z2 = p2.z.mul_ref(&p2.z) % &self.p;
            // x1 * z2^2 = x2 * z1^2 ?
            result = result && (p1.x.mul_ref(&z2) % &self.p == p2.x.mul_ref(&z1) % &self.p);
            // z1^3
            z1 = p1.z.mul_ref(&z1) % &self.p;
            // z2^3
            z2 = p2.z.mul_ref(&z2) % &self.p;
            // y1 * z2^3 = y2 * z1^3 ?
            result = result && (p1.y.mul_ref(&z2) % &self.p == p2.y.mul_ref(&z1) % &self.p);
        }
        result
    }

    /// Negate P
    fn negate(&self, p_p: &Point<T>) -> Point<T> {
        Point {
            axis: EcAxis::Proj,
            x: p_p.x.clone(),
            y: self.p.clone() - p_p.y.clone(),
            z: p_p.z.clone(),
        }
    }

    /// EC point double : Q = 2 * P
    fn double(&self, p_p: &Point<T>) -> Point<T> {
        let (mut t1, mut t2, mut t3) = (p_p.x.clone(), p_p.y.clone(), p_p.z.clone());
        let (mut t4, mut t5);

        if t2.is_zero() || t3.is_zero() {
            return Point {
                axis: EcAxis::Proj,
                x: T::one(),
                y: T::one(),
                z: T::zero(),
            };
        }
        if self.a == self.p.sub_ref(&T::from(3u32)) % &self.p {
            t4 = t3.mul_ref(&t3) % &self.p;
            t5 = t1.sub_ref(&t4) % &self.p;
            t4 = t1.add_ref(&t4) % &self.p;
            t5 = t4.mul_ref(&t5) % &self.p;
            t4 = T::from(3u32) * &t5 % &self.p;
        } else {
            t4 = self.a.clone();
            t5 = t3.mul_ref(&t3) % &self.p;
            t5 = t5.mul_ref(&t5) % &self.p;
            t5 = t4.mul_ref(&t5) % &self.p;
            t4 = t1.mul_ref(&t1) % &self.p;
            t4 = T::from(3u32) * &t4 % &self.p;
            t4 = t4.add_ref(&t5) % &self.p;
        }
        t3 = t2.mul_ref(&t3) % &self.p;
        t3 = Number::mod_cal(&(T::from(2u32) * &t3), &self.p);
        t2 = t2.mul_ref(&t2) % &self.p;
        t5 = t1.mul_ref(&t2) % &self.p;
        t5 = T::from(4u32) * &t5 % &self.p;
        t1 = t4.mul_ref(&t4) % &self.p;
        t1 = Number::mod_cal(&(t1.sub_ref(&(T::from(2u32) * &t5))), &self.p);
        t2 = t2.mul_ref(&t2) % &self.p;
        t2 = T::from(8u32) * &t2 % &self.p;
        t5 = t5.sub_ref(&t1) % &self.p;
        t5 = t4.mul_ref(&t5) % &self.p;
        t2 = Number::mod_cal(&(t5.sub_ref(&t2)), &self.p);

        Point {
            axis: EcAxis::Proj,
            x: t1,
            y: t2,
            z: t3,
        }
    }

    /// EC Point<T> add : Q = P1 + P2
    fn add(&self, p_p1: &Point<T>, p_p2: &Point<T>) -> Point<T> {
        if p_p1.z.is_zero() {
            return p_p2.clone();
        }
        if p_p2.z.is_zero() {
            return p_p1.clone();
        }
        if self.equals(p_p1, p_p2) {
            return self.double(p_p1);
        }

        let (mut t1, mut t2, mut t3) = (p_p1.x.clone(), p_p1.y.clone(), p_p1.z.clone());
        let (mut t4, mut t5) = (p_p2.x.clone(), p_p2.y.clone());
        let (mut t6, mut t7): (T, T);

        t6 = T::one();
        if !p_p2.z.is_one() {
            t6 = p_p2.z.clone() % &self.p;
            t7 = t6.mul_ref(&t6) % &self.p;
            t1 = t1.mul_ref(&t7) % &self.p; // U0 = P1.x * P2.z^2
            t7 = t6.mul_ref(&t7) % &self.p;
            t2 = t2.mul_ref(&t7) % &self.p; // S0 = P1.y * P2.z^3
        }

        t7 = t3.mul_ref(&t3) % &self.p;
        t4 = t4.mul_ref(&t7) % &self.p; // U1 = P2.x * P1.z^2
        t7 = t3.mul_ref(&t7) % &self.p;
        t5 = t5.mul_ref(&t7) % &self.p; // S1 = P2.y * P1.z^3
        t4 = t1.sub_ref(&t4) % &self.p; // W = U0 - U1
        t5 = t2.sub_ref(&t5) % &self.p; // R = S0 - S1

        if t4.is_zero() {
            if t5.is_zero() {
                return Point {
                    axis: EcAxis::Proj,
                    x: T::zero(),
                    y: T::zero(),
                    z: T::zero(),
                };
            } else {
                return Point {
                    axis: EcAxis::Proj,
                    x: T::one(),
                    y: T::one(),
                    z: T::zero(),
                };
            }
        }
        t1 = (T::from(2u32) * &t1 - &t4) % &self.p; // T = U0 + U1
        t2 = (T::from(2u32) * &t2 - &t5) % &self.p; // M = S0 + S1
        if !p_p2.z.is_one() {
            t3 = t3.mul_ref(&t6) % &self.p;
        }

        t3 = Number::mod_cal(&(t3.mul_ref(&t4)), &self.p); // Z3 = z1*z2*W
        t7 = t4.mul_ref(&t4) % &self.p;
        t4 = t4.mul_ref(&t7) % &self.p;
        t7 = t1.mul_ref(&t7) % &self.p;
        t1 = t5.mul_ref(&t5) % &self.p;
        t1 = Number::mod_cal(&(t1.sub_ref(&t7)), &self.p); // X3 = R^2 - T*W^2

        t7 = (t7.sub_ref(&(T::from(2u32) * &t1))) % &self.p; // V = T*W^2 - 2*x3
        t5 = t5.mul_ref(&t7) % &self.p;
        t4 = t2.mul_ref(&t4) % &self.p;
        t2 = t5.sub_ref(&t4) % &self.p;
        t2 = Number::mod_inv(&T::from(2u32), &self.p) * &t2 % &self.p; // Y3 = (V*R - M*W^3)/2
        t2 = Number::mod_cal(&t2, &self.p);

        Point {
            axis: EcAxis::Proj,
            x: t1,
            y: t2,
            z: t3,
        }
    }

    /// Generate a random point
    fn gen_point(&self) -> Point<T> {
        let (mut x, mut y, mut y2): (T, T, T);
        let p = self.p.clone();

        loop {
            x = Number::gen_rand(&T::from(1u32), &(p.sub_ref(&T::from(1u32))));

            // y^2 = x^3 + ax + b = (x^2 + a) * x + b
            y2 = Number::mod_cal(&((x.mul_ref(&x) + &self.a) * &x + &self.b), &self.p);
            y = Number::mod_sqrt(&y2, &self.p);

            let pt = Point {
                axis: EcAxis::Affine,
                x,
                y,
                z: T::one(),
            };
            if self.on_curve(&pt) {
                return pt;
            }
        }
    }

    /// Find points on curve at x
    /// - returns: ((x, y, 1), (x, -y, 1)) or not found exception
    fn point_from_x(&self, x: &T, yt: u32) -> Point<T> {
        assert!(x < &self.p);

        // y^2 = x^3 + ax + b = (x^2 + a) * x + b
        let y2 = Number::mod_cal(&((x.mul_ref(x) + &self.a) * x + &self.b), &self.p);
        let mut y = Number::mod_sqrt(&y2, &self.p);

        let bit0 = y.clone() % T::from(2u32);
        if bit0 != T::from(yt) {
            y = self.p.sub_ref(&y);
        }

        Point {
            axis: EcAxis::Proj,
            x: x.clone(),
            y,
            z: T::one(),
        }
    }

    /// Find points on curve at (x, y)
    fn point_from_xy(&self, x: &T, y: &T) -> Point<T> {
        Point {
            axis: EcAxis::Proj,
            x: x.clone(),
            y: y.clone(),
            z: T::one(),
        }
    }

    /// OS2ECPP: Decode a point on this curve which has been encoded using point
    /// compression(X9.62) returning the point.
    #[allow(unused_assignments)]
    fn decode_point(&self, enc: &[u8]) -> Point<T> {
        if enc.len() == 1 {
            if enc[0] == 0 {
                return self.get_zero();
            } else {
                panic!("Invalid point compression");
            }
        }

        let mut x = T::zero();
        let mut y = T::zero();

        match enc[0] {
            // compressed
            0x02 | 0x03 => {
                let ln = enc.len() - 1;
                let xb = &enc[1..ln + 1];
                let x = T::from_bytes_be(xb);
                let yt = enc[0] & 1;
                return self.point_from_x(&x, yt.into());
            }
            // uncompressed
            0x04 => {
                let ln = (enc.len() - 1) / 2;
                let xb = &enc[1..(ln + 1)];
                let yb = &enc[(ln + 1)..enc.len()];
                x = T::from_bytes_be(xb);
                y = T::from_bytes_be(yb);
            }
            _ => panic!("Invalid point compression"),
        }
        Point {
            axis: EcAxis::Proj,
            x,
            y,
            z: T::one(),
        }
    }

    /// Returns the field element encoded with point compression.
    fn get_encoded(&self, p: &Point<T>) -> Vec<u8> {
        let q = self.to_affine(p);

        let y = q.y;
        let x = q.x;
        // y~ = y mod 2
        let yt = y % T::from(2u32);
        let hd = if yt == T::one() { 0x03 } else { 0x02 };

        let xb = x.to_bytes_be();
        let mut po = vec![0u8; xb.len() + 1];
        po[0] = hd;
        po[1..].copy_from_slice(&xb[..xb.len()]);

        po
    }

    /// Order of point P
    fn calc_order(&self, p_p: &Point<T>) -> T {
        let mut m = T::one();

        while m < self.p.add_ref(&T::one()) {
            if self.mul(p_p, &m).is_zero() {
                return m;
            }
            m = m + T::one();
        }
        panic!("Invalid order");
    }

    /// Get generator Point<T>
    #[inline]
    fn get_gp(&self) -> Point<T> {
        self.p_g.clone()
    }

    /// Get order
    #[inline]
    fn get_order(&self) -> T {
        self.n.clone()
    }

    /// Get cofactor
    #[inline]
    fn get_cofac(&self) -> T {
        self.h.clone()
    }

    /// Display point data
    fn print(&self, s: &str, p: &Point<T>) {
        println!("{}: [{}, {}, {}]", s, p.x, p.y, p.z);
    }
}

impl<T: Number> EcpJ<T> {
    pub fn new(ec_name: &str) -> EcpJ<T> {
        let ec = EcParam::new(ec_name);

        EcpJ {
            a: ec.get_a(),
            b: ec.get_b(),
            p: ec.get_prime(),
            n: ec.get_order(),
            h: ec.get_cofactor(),
            p_zero: Point {
                axis: EcAxis::Proj,
                x: T::zero(),
                y: T::zero(),
                z: T::zero(),
            },
            p_g: Point {
                axis: EcAxis::Proj,
                x: ec.get_gx(),
                y: ec.get_gy(),
                z: T::one(),
            },
        }
    }
}
