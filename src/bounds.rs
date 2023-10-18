use rust_decimal_macros::dec;
use rust_decimal::Decimal;
use rust_decimal::MathematicalOps;

const PI: Decimal = dec!(3.1415926535897932384626433);
// All of the following are rounded up

/* 2 * (sqrt(2) + 4 / sqrt(5) - 2) */
const A53: Decimal = dec!(2.4061358887458537125780821953);
/* 27 - 8*sqrt(2)/3/sqrt(5)-10*sqrt(5)/3*/
const B53: Decimal = dec!(17.859891989577565851732288138);
/* (2*sqrt(2)+20*sqrt(5))/3 */
const C53: Decimal = dec!(15.849928891580661982629862905);
// I'm done writing these by hand god damn
const D53: Decimal = dec!(-0.3346474215415660147243670508);
const E53: Decimal = dec!(5.4920159537189991993955118232);

const ONEPLUS: Decimal = dec!(1.0000000000000002220446049250);
const ONEMINUS: Decimal = dec!(0.9999999999999988897769753748);

const EPS: Decimal = dec!(0);
const DEPTH: Decimal = dec!(10);

// Caching these so the macro doesn't have to run repeatedly
const DEC_2: Decimal = dec!(2);
const DEC_1: Decimal = dec!(1);

// This is the main function. 
// I just want it called manually by a different main() method.
pub fn init(u128 argc, char)
{
    
}

fn upper1 (k: Decimal) -> Decimal
{
    let val = k / PI;
    return val.sqrt().unwrap_or(DEC_1);
}

fn lower1 (k: Decimal) -> Decimal
{
    let val = ((k / PI) + dec!(4)).sqrt().unwrap_or(DEC_1); /* x < sqrt(Pi*(k-1) + 4) */
    return (val - DEC_2) / PI;
}

// This is better than lower1() for k >= 76
fn lower2 (k: Decimal) -> Decimal
{
    assert!(k > dec!(6));
    let mut val = (((Decimal::from(k) - dec!(6)) * PI) + DEC_2).sqrt().unwrap();
    val -= DEC_2.sqrt().unwrap();
    return val / PI;
}

// This is better for k >= 2798
// (c is the same as k, I just wanted to convert it from u128 to decimal)
fn lower3 (k: Decimal) -> Decimal
{
    let twopi = DEC_2 * PI;
    let mut r = lower2(k);
    let mut l = k;
    let mut f = k;
    let mut lk = k;
    
    assert!(dec!(2798) <= k && k < dec!(13647034));
    // l2 is better (i.e. larger) for k >= 76
    // The maximal number of loops for 2798 <= k < 13647034 is 5 for several values up to k=15656
    while r > l
    {
        l = r; /* Save old lower bound */
        f = B53 + C53 / l + D53 / (l * l) + E53 / (l * l * l);
        // Solve the 2nd degree equation: Pi*r^2 + a*r + f-k = 0
        r = A53 * A53 + DEC_2 * twopi * (lk - f);
        r = (r.sqrt().unwrap_or(DEC_1) - A53) / twopi;
    }
    let small = Decimal::from_scientific("8.88178419700125e-16").unwrap();
    r -= small * r;

    // Missing some assertions from the OG but it should be fine, right?
    assert!(r >= dec!(21));
    return r;
}

/* This lower bound is better than lower3() for k >= 13647034,
and better than lower3a() for k >= 93598905 */
// k must be between 1e7 and 1e8 for this to guaranteed work
fn lower4 (k: Decimal) -> Decimal
{
    let mut r = dec!(0);
    let num = dec!(30.84274723);
    let two_thirds = dec!(0.666666666666666666667);
    let mut oldr = r;
    // Introducing this wacky boolean because otherwise it would automatically terminate
    // In the original C code this was a do while loop
    assert!(dec!(10000000) <= k && k < dec!(100000000));
    let mut test = false;
    while !test
    {
        oldr = r;
        /* We perform two iterations, since if r is smaller than the root r0
         * of k = Pi*r^2 + 30.84274723*r^(2/3), then the next iteration is larger
         * than r0, and vice-versa. */
        r = (k - num * r.powd(two_thirds) / PI).sqrt().unwrap();
        r = (k - num * r.powd(two_thirds) / PI).sqrt().unwrap();

        if r != oldr
        {
            test = true
        }
    }
    let small_2 = Decimal::from_scientific("4.44089209850063e-16").unwrap();
    r -= small_2 * r;

    return r;
}

fn lower3a (k: Decimal) -> Decimal
{
    let mut r_0 = DEC_1;
    let mut i = DEC_1;
    assert!(DEC_2 <= k);
    
    r_0 = if k < dec!(76) {lower1(k)} else {lower2(k)};

    let m = if r_0 <= dec!(18.94427191) {DEC_2} 
            else {((r_0 - DEC_2) / DEC_2 / (r_0 - DEC_1).sqrt().unwrap()).floor()};

    assert!(m >= DEC_2);
    // T >= 4 * r * (m - the sum from i=1..m of (i / sqrt(i^2 + 1)))
    // We want a lower bound for T, so we need to round i/sqrt(i^2 + 1) upwards
    r_0 = r_0.floor();
    let mut sum = dec!(0);
    let mut temp = dec!(0);
    let mut i = DEC_1;
    // I have to use a while loop here because iteration is not implemented for Decimal
    while i < m
    {
        temp = ((i * i) + DEC_1).sqrt().unwrap_or(DEC_1);
        temp = DEC_1 / temp;
        sum += temp;
        i += DEC_1;
    }
    // I probably didn't need to introduce this x term but the original code was weirdly
    // written (Because C) and I'm trying to follow it as best as I can
    let mut x = sum;
    x = m - sum; // m - sum(...)
    x -= DEC_1; // Accounts for 4r term
    x *= DEC_2; // 4 * (m - add(...))
    /* Now we have to solve k = Pi*r^2 - a*r - b with a = m and b = 8m - 1
    */
    let mut y = PI * (k - (dec!(8) * m) - DEC_1);
    y *= dec!(4);
    let mut z = x * x;
    z += y;
    z = z.sqrt().unwrap_or(DEC_1); // a^2 + 4*Pi*(b + k)^(1/2)
    z += x;
    z /= PI;
    return dec!(0.5) * z;
}

fn lower3b (k: Decimal) -> Decimal
{
    let four = dec!(4);
    assert!(DEC_2 <= k);
    let r_0 = if k < dec!(76) {lower1(k)} else {lower2(k)};

    /* m is 2m from the paper */
    let m = if r_0 <= dec!(18.94427191) {four} 
            else {((r_0 - DEC_2) / (r_0 - DEC_1).sqrt().unwrap()).floor()};

    let mut z = r_0;
    let mut sum = dec!(0); // Will contain the sum (This is x53 in the original code)
    let mut i = dec!(3);
    let mut temp = dec!(0);
    // I have to use a while loop here because iteration is not implemented for Decimal
    while i <= m
    {
        // i/2/sqrt((i/2)^2 + 1) = i/sqrt(i^2 + 4)
        // temp is y53 in the original code
        temp = ((i * i) + four).sqrt().unwrap();
        temp = i / temp;
        sum += temp;
        i += DEC_1;
    }
    let root_2 = DEC_2.sqrt().unwrap();
    sum += root_2; // sqrt(2) + add(...)
    sum = ((m - sum) * DEC_2) - four; // -4 takes into account the 4r term

    /* Now we have to solve k = Pi*r^2 - a*r - b with a = sum and b = 4m-1
    *  Which has root 1/2/Pi*(a + (a^2 + 4*Pi*(b + k))^(1/2)) */
    temp = PI * (k - four * m - DEC_1);
    temp *= four;
    let mut z = (sum * sum) + temp; // a^2+4*Pi*(b+k)
    z = z.sqrt().unwrap() + sum;
    z /= PI;
    return dec!(0.5) * z;
}

/** Return an upper bound for the number of points in the circle of center
    xc+i*yc and radius r, with xcmin <= xc <= xcmax, ycmin <= yc <= ycmax */
fn upper_bound_k (xcmin: Decimal, 
    xcmax: Decimal, 
    ycmin: Decimal, 
    ycmax: Decimal, 
    r: Decimal) -> Decimal
{
    let xmin = (xcmin - r).ceil();
    let xmax = xcmax - r; // No need to floor here
    let r2 = (r * r) * ONEPLUS; /* r2 >= r^2 */

    let mut x = xmin;
    let mut s = dec!(0);
    while x <= xmax
    {
        /* t = x - xc */
        let mut tmin;
        // There was also some stuff about a tmax variable
        // but it was all commented out lol
        if x <= xcmin
        {
            // tmax = tmin * tmin; /* not used below */
            let tmax = x - xcmin;
            tmin = (tmax * tmax) * ONEMINUS;
        }
        else if xcmax <= x /* 0 <= tmin */
        {
            tmin = x - xcmax;
            tmin = (tmin * tmin) * ONEMINUS;
        }
        else
        {
            tmin = dec!(0);
        }
        /* r3 = upper(r2 - (x-xc)^2) */
        let r3 = r2 - tmin;
        let sr3 = r3.sqrt().unwrap() * ONEPLUS;
        let ymin = (ycmin - sr3).ceil();
        let ymax = (ycmax + sr3).floor();
        s += ymax - ymin + DEC_1;

        x += DEC_1;
    }
    return s;
}

/* return a lower bound for the number of points in the circle of center 
   xcmin + i*ycmin and radius r */
fn lower_bound_k (xcmin: Decimal,
    ycmin: Decimal,
    r: Decimal) -> Decimal
{
    let xmin = (xcmin - r).ceil();
    let xmax = xcmin + r; // No need to floor here
    let r2 = (r * r) * ONEMINUS;

    let mut x = xmin;
    let mut s = dec!(0);
    while x <= xmax
    {
        let tmin = x - xcmin;
        let tmax = (tmin * tmin) * ONEPLUS;
        let r3 = r2 - tmax;
        let sr3 = r3.sqrt().unwrap() * ONEMINUS;
        let ymin = (ycmin - sr3).ceil();
        let ymax = (ycmin + sr3).floor();
        s += ymax - ymin + DEC_1;
        x += DEC_1;
    }
    return s;
}

/** Return an upper bound for the number of points in the circle of center xc + i*yc
    and radius r, with xcmin <= xc <= xcmax, ycmin <= yc <= ycmax, by doing at most
    2^depth recursive cuttings.
    At the end, we want to decide if the upper bound k is >= threshold or not 
*/
fn upper_bound_k_rec (xcmin: Decimal,
    xcmax: Decimal,
    ycmin: Decimal,
    ycmax: Decimal,
    r: Decimal,
    depth: Decimal,
    threshold: Decimal) -> Decimal
{
    let mut k = upper_bound_k(xcmin, xcmax, ycmin, ycmax, r);
    if k < threshold
    {
        return k; // Subdividing will not improve the upper bound
    }
    if depth == dec!(0) // k >= threshold, but reached end of recursion
    {
        return k;
    }
    if ycmax <= xcmin /* ensure xc <= yc (symmetry) */
    {
        return threshold - DEC_1;
    }
    if xcmax - xcmin >= ycmax - ycmin // Cut x in two
    {
        let xcmid = (xcmin + xcmax) * dec!(0.5);
        k = upper_bound_k_rec(xcmin, xcmid, ycmin, ycmax, r, depth - DEC_1, threshold);
        if k >= threshold
        {
            return k;
        }
        else
        {
            return upper_bound_k_rec(xcmid, xcmax, ycmin, ycmax, r, depth - DEC_1, threshold);
        }
    }
    else // Cut y in two
    {
        let ycmid = (ycmin + ycmax) * dec!(0.5);
        k = upper_bound_k_rec(xcmin, xcmax, ycmin, ycmid, r, depth - DEC_1, threshold);
        if k >= threshold /* Same as above */
        {
            return k;
        }
        else
        {
            return upper_bound_k_rec(xcmin, xcmax, ycmid, ycmax, r, depth - DEC_1, threshold);
        }
    }
}

/** Return a lower bound for the number of points in the circle of center xc+i*yc
*   and radius r, with xcmin <= xc <= xcmax, ycmin <= yc <= ycmax, by doing at most
*   2^depth recursive cuttings. 
*/
fn lower_bound_k_rec (xcmin: Decimal,
    xcmax: Decimal,
    ycmin: Decimal,
    ycmax: Decimal,
    r: Decimal,
    depth: Decimal,
    threshold: Decimal) -> Decimal
{
    let mut k = lower_bound_k(xcmin, ycmin, r);
    if k >= threshold
    {
        return k; // Subdividing will not improve the upper bound
    }
    if depth == dec!(0) // k < threshold, but reached end of recursion
    {
        return k;
    }
    if ycmax <= xcmin /* ensure xc <= yc (symmetry) */
    {
        return k;
    }
    if xcmax - xcmin >= ycmax - ycmin // Cut x in two
    {
        let xcmid = (xcmin + xcmax) * dec!(0.5);
        k = lower_bound_k_rec(xcmin, xcmid, ycmin, ycmax, r, depth - DEC_1, threshold);
        if k >= threshold
        {
            return k;
        }
        else
        {
            return lower_bound_k_rec(xcmid, xcmax, ycmin, ycmax, r, depth - DEC_1, threshold);
        }
    }
    else // Cut y in two
    {
        let ycmid = (ycmin + ycmax) * dec!(0.5);
        k = lower_bound_k_rec(xcmin, xcmax, ycmin, ycmid, r, depth - DEC_1, threshold);
        if k >= threshold /* Same as above */
        {
            return k;
        }
        else
        {
            return upper_bound_k_rec(xcmin, xcmax, ycmid, ycmax, r, depth - DEC_1, threshold);
        }
    }
}

/** This code is useless but I'm keeping it in case it's needed suddenly */
fn lower5 (k: Decimal, depth: Decimal) -> Decimal
{
    return lower2(k);
}
/** Same with this code lol */
fn upper5 (k: Decimal, depth: Decimal) -> Decimal
{
    return upper1(k);
}

/** return an upper bound of r_k and if rk2 <> NULL, put in rk2
*   an upper bound of r_k^2.
*   We use r_k < sqrt((k-1)/Pi)
*/
fn upper (k: Decimal) -> Decimal
{
    // I couldn't find a good equivalent to this code in rust.
    return upper1(k);
}

fn lower (k: Decimal) -> Decimal
{
    if k < dec!(76)
    {
        return lower1(k);
    }
    if k < dec!(2798)
    {
        return lower2(k);
    }
    if k < dec!(13647034)
    {
        return lower3(k);
    }
    else
    {
        return lower4(k);
    }
}