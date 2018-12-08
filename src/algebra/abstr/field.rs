use algebra::abstr::Ring;
use std::ops::{Sub, SubAssign, Div, DivAssign};

/// Field
///
/// <a href="https://en.wikipedia.org/wiki/Field_(mathematics)">https://en.wikipedia.org/wiki/Field_(mathematics)</a>
pub trait Field : Ring + Sub<Self, Output = Self> + SubAssign<Self> + Div<Self, Output = Self> + DivAssign<Self>
{

}
