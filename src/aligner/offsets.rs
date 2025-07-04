use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, Not, Shr};

use num::{FromPrimitive, Unsigned, One, Bounded};
use num::traits::{SaturatingAdd, SaturatingSub};

pub trait OffsetType: FromPrimitive + Unsigned + PartialEq + Eq
    + PartialOrd + Ord + Default + Copy + Hash + Debug + Bounded + SaturatingSub + SaturatingAdd
    + Not<Output=Self> + BitAnd<Output=Self> + Shr<Output=Self>
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn as_isize(&self) -> isize;
    fn increase_one(&self) -> Self;
}

impl OffsetType for u8 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u16 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u32 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }
    
    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u64 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }
    
    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

#[cfg(test)]
mod tests {
    use super::OffsetType;

    #[test]
    fn test_offset_u8() {
        let o = <u8 as OffsetType>::new(5);
        assert_eq!(o.as_usize(), 5);
        assert_eq!(o.as_isize(), 5);
        assert_eq!(o.increase_one(), <u8 as OffsetType>::new(6));
    }

    #[test]
    fn test_offset_u16() {
        let o = <u16 as OffsetType>::new(10);
        assert_eq!(o.as_usize(), 10);
        assert_eq!(o.increase_one(), <u16 as OffsetType>::new(11));
    }

    #[test]
    fn test_offset_u32() {
        let o = <u32 as OffsetType>::new(42);
        assert_eq!(o.as_usize(), 42);
        assert_eq!(o.increase_one(), <u32 as OffsetType>::new(43));
    }

    #[test]
    fn test_offset_u64() {
        let o = <u64 as OffsetType>::new(7);
        assert_eq!(o.as_usize(), 7);
        assert_eq!(o.increase_one(), <u64 as OffsetType>::new(8));
    }
}
