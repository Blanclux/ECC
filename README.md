# ECC  

Elliptic Curve Calculation Library over Finite Fields

## Features

- The library implements elliptic curve operations in pure Rust.
- Support point add, double, mul and other operations for Affine and Projective coordinate
- Using num-bigint crate and ibig crate for multiprecition integer

## modules

- lib  
  Point structure and tarit  
- number  
  Number trait (the general interface for number operations)
  for BigInt and IBig
- ec_param  
  Elliptic curve parameters
- ecc_a  
  Elliptic curve operations (Affine)  
- ecc_j  
  Elliptic curve operations (Projective - Jacobian)  
