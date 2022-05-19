//!
//! Elliptic Curve Point Parameter Interface
//!
//  Written by blanclux
//  This software is distributed on an "AS IS" basis WITHOUT WARRANTY OF ANY KIND.

use std::str;

use super::number::Number;

///   EC curve parameters
///   id,  a  b,  p,  g (04||X||Y),  order,   h
pub struct EcParam<'a> {
    pub id: &'a str,
    pub a: &'a [u8],
    pub b: &'a [u8],
    pub p: &'a [u8],
    pub g: &'a [u8],
    pub n: &'a [u8],
    pub h: &'a [u8],
}

static CURVE_PARAMS: [EcParam; 8] = [
    EcParam {
        id: "secp160k1",
        a: b"0000000000000000000000000000000000000000",
        b: b"0000000000000000000000000000000000000007",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFAC73",
        g: b"043B4C382CE37AA192A4019E763036F4F5DD4D7EBB938CF935318FDCED6BC28286531733C3F03C4FEE",
        n: b"0100000000000000000001B8FA16DFAB9ACA16B6B3",
        h: b"01",
    },
    EcParam {
        id: "secp160r1",
        a: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFC",
        b: b"1C97BEFC54BD7A8B65ACF89F81D4D4ADC565FA45",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFF",
        g: b"044A96B5688EF573284664698968C38BB913CBFC8223A628553168947D59DCC912042351377AC5FB32",
        n: b"0100000000000000000001F4C8F927AED3CA752257",
        h: b"01"
    },
    EcParam {
        id: "secp192k1",
        a: b"000000000000000000000000000000000000000000000000",
        b: b"000000000000000000000000000000000000000000000003",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFEE37",
        g: b"04DB4FF10EC057E9AE26B07D0280B7F4341DA5D1B1EAE06C7D9B2F2F6D9C5628A7844163D015BE86344082AA88D95E2F9D",
        n: b"FFFFFFFFFFFFFFFFFFFFFFFE26F2FC170F69466A74DEFD8D",
        h: b"01"
    },
    EcParam {
        id: "secp224k1",
        a: b"00000000000000000000000000000000000000000000000000000000",
        b: b"00000000000000000000000000000000000000000000000000000005",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFE56D",
        g: b"04A1455B334DF099DF30FC28A169A467E9E47075A90F7E650EB6B7A45C7E089FED7FBA344282CAFBD6F7E319F7C0B0BD59E2CA4BDB556D61A5",
        n: b"010000000000000000000000000001DCE8D2EC6184CAF0A971769FB1F7",
        h: b"01"
    },
    EcParam {
    	id: "secp224r1",
        a: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE",
        b: b"B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001",
        g: b"04B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34",
        n: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D",
        h: b"01"
    },
    EcParam {
        id: "secp256k1",
        a: b"0000000000000000000000000000000000000000000000000000000000000000",
        b: b"0000000000000000000000000000000000000000000000000000000000000007",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
        g: b"0479BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
        n: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
        h: b"01"
    },
    EcParam {
        id: "secp384r1",
        a: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC",
        b: b"B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF",
        p: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF",
        g: b"04AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB73617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F",
        n: b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973",
        h: b"01"
    },
    EcParam {
        id: "secp512r1",
        a: b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC",
        b: b"0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00",
        p: b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
        g: b"0400C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650",
        n: b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409",
        h: b"01"
    },
];

impl EcParam<'_> {
    pub fn new(ecid: &str) -> EcParam {
        let mut id: &str = CURVE_PARAMS[0].id;
        let mut ai: &[u8] = CURVE_PARAMS[0].a;
        let mut bi: &[u8] = CURVE_PARAMS[0].b;
        let mut pi: &[u8] = CURVE_PARAMS[0].p;
        let mut gi: &[u8] = CURVE_PARAMS[0].g;
        let mut ni: &[u8] = CURVE_PARAMS[0].n;
        let mut hi: &[u8] = CURVE_PARAMS[0].h;

        for prm in CURVE_PARAMS.iter() {
            if prm.id == ecid {
                id = prm.id;
                ai = prm.a;
                bi = prm.b;
                pi = prm.p;
                gi = prm.g;
                ni = prm.n;
                hi = prm.h;
            }
        }
        EcParam {
            id,
            a: ai,
            b: bi,
            p: pi,
            g: gi,
            n: ni,
            h: hi,
        }
    }

    pub fn get_ecid(&self) -> &str {
        self.id
    }

    pub fn set_parameters(&mut self, id: &str) {
        for prm in CURVE_PARAMS.iter() {
            if prm.id == id {
                self.id = prm.id;
                self.a = prm.a;
                self.b = prm.b;
                self.p = prm.p;
                self.g = prm.g;
                self.n = prm.n;
                self.h = prm.h;
            }
        }
    }

    pub fn check_parameters(&self, id: &str) -> bool {
        for prm in CURVE_PARAMS.iter() {
            if prm.id == id {
                return true;
            }
        }
        false
    }
}

pub trait ParamOp<T: Number> {
    fn get_prime(&self) -> T;
    fn get_a(&self) -> T;
    fn get_b(&self) -> T;
    fn get_order(&self) -> T;
    fn get_cofactor(&self) -> T;
    fn get_gx(&self) -> T;
    fn get_gy(&self) -> T;
}

impl<T: Number> ParamOp<T> for EcParam<'_> {
    fn get_prime(&self) -> T {
        T::from_bytes_radix(self.p, 16)
    }

    fn get_a(&self) -> T {
        T::from_bytes_radix(self.a, 16)
    }

    fn get_b(&self) -> T {
        T::from_bytes_radix(self.b, 16)
    }

    fn get_order(&self) -> T {
        T::from_bytes_radix(self.n, 16)
    }

    fn get_cofactor(&self) -> T {
        T::from_bytes_radix(self.h, 16)
    }

    fn get_gx(&self) -> T {
        let ln = (self.g.len() - 2) / 2;
        let gb = &self.g[2..(ln + 2)];
        T::from_bytes_radix(gb, 16)
    }

    fn get_gy(&self) -> T {
        let ln = (self.g.len() - 2) / 2;
        let gb = &self.g[(ln + 2)..self.g.len()];
        T::from_bytes_radix(gb, 16)
    }
}
