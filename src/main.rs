use ecgfp5::field::{GFp, GFp5};
use ecgfp5::curve::Point;
use rand::Rng;

fn sgn(x: GFp5) -> bool {
    x.0.iter().map(|f| f.to_u64() % 2u64).sum::<u64>() % 2 == 1
}

fn is_square(x: GFp5) -> bool {
    let leg = x.legendre();
    (leg.square() - leg).iszero() == u64::MAX
}

fn sswu(u: GFp5) -> Point {

    let two_thirds: GFp5 = GFp5::from_u64(6148914689804861441, 0, 0, 0, 0).0;

    // coefficients for double-odd form y^2 = x(x^2 + Ax + B)
    let A: GFp5 = GFp5::from_u64(2, 0, 0, 0, 0).0;
    let B: GFp5 = GFp5::from_u64(0, 263, 0, 0, 0).0;

    // coefficients for Short Weierstrass form Y^2 = X^3 + A_sw*x + B_sw
    // A_sw = (3B - A^2)/3
    // B_sw = A(2A^2 -9B)/27 
    let A_sw: GFp5 = GFp5::from_u64(6148914689804861439, 263, 0, 0, 0).0;
    let B_sw: GFp5 = GFp5::from_u64(15713893096167979237, 6148914689804861265, 0, 0, 0).0;

    // Z computed using SageMath
    // Z_sw = -4 - z = 18446744069414584317 + 18446744069414584320*z
    let Z_sw: GFp5 = GFp5::from_u64(GFp::MOD-4, GFp::MOD-1, 0, 0, 0).0;

    let denom_part = Z_sw*u.square();
    let denom = denom_part.square() + denom_part;
    let tv1 = denom.invert();
    
    let x1 = match tv1.iszero() == u64::MAX {
        true => B_sw / (Z_sw * A_sw),
        false => (-B_sw/A_sw)*(GFp5::ONE + tv1)
    };
    let x2 = denom_part * x1;

    // g(x) = X^3 + A_sw*X + B_sw
    let gx1 = x1*x1.square() + A_sw*x1 + B_sw;
    let gx2 = x2*x2.square() + A_sw*x2 + B_sw;

    let (x_sw, y_pos) = match is_square(gx1) {
        true => (x1, gx1.sqrt().0),
        false => (x2, gx2.sqrt().0),
        _ => unreachable!(),
    };

    let x_cand = x_sw - two_thirds;
    let y_cand = match sgn(y_pos) == sgn(u) {
        true => y_pos,
        false => -y_pos,
    };

    Point::decode(y_cand/x_cand).0
}


fn main() {
    let mut rng = rand::thread_rng();
    let x: Vec<u64> = (0..5).map(|_| rng.gen::<u64>()).collect();
    let u = GFp5::from_u64_reduce(x[0], x[1], x[2], x[3], x[4]);
    println!("u: {}", u);
    println!("point: {}", sswu(u));

}


// SageMath code to generate Z
// from RFC 9380 Appendix H.2

/*
    p = 18446744069414584321
    P.<x> = PolynomialRing(GF(p))
    m = x^5 - 3 
    F.<z> = GF(p^5, modulus=m)

    one = F.vector_space(map=False).gen(0); one
    z = F.vector_space(map=False).gen(1); z

    a = 2*F(one)
    b = 263*F(z)

    print(a)
    print(b)

    A = (3*b - a*a) / 3
    B = a*(2*a*a - 9*b) / 27

    print(A)
    print(B)

    # Arguments:
    # - F, a field object, e.g., F = GF(2^521 - 1)
    # - A and B, the coefficients of the curve y^2 = x^3 + A * x + B
    def find_z_sswu(F, A, B):
    R.<xx> = F[] # Polynomial ring over F
    g = xx^3 + F(A) * xx + F(B) # y^2 = g(x) = x^3 + A * x + B
    ctr = F.gen()
    while True:
    for Z_cand in (F(ctr), F(-ctr)):
        # Criterion 1: Z is non-square in F.
        if is_square(Z_cand):
        continue
        # Criterion 2: Z != -1 in F.
        if Z_cand == F(-1):
        continue
        # Criterion 3: g(x) - Z is irreducible over F.
        if not (g - Z_cand).is_irreducible():
        continue
        # Criterion 4: g(B / (Z * A)) is square in F.
        if is_square(g(B / (Z_cand * A))):
        return Z_cand
    ctr += 1
    
    Z = find_z_sswu(F, A, B)
    print(Z)
    print(-Z)


    -------------------------------------------------------------------
     a: 2
     b: 263*z
     A: 263*z + 6148914689804861439
     B: 6148914689804861265*z + 15713893096167979237
     Z: 18446744069414584320*z + 18446744069414584317
    -Z: z + 4
*/