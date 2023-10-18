use rust_decimal_macros::dec;
use rust_decimal::Decimal;
use rust_decimal::MathematicalOps;

const THREE: Decimal = dec!(3);

pub fn sqrt(arg: Decimal) -> Decimal
{
    let mut val = arg;
    for i in 1..25
    {
        val = ((val * val * val) + (THREE * arg * val)) / (THREE * (val * val) + arg);
    }
    return val;
}