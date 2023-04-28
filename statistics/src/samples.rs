// TODO: use ExactSizeIterator

/// This trait is similiar to `ExactSizeIterator` but we adjust it for out own needs.
// pub trait Samples<T>: std::iter::Iterator<Item = T> {
//     fn is_empty(&self) -> bool {
//         self.len() == 0
//     }
//     fn len(&self) -> usize;
// }

/// Should be the same as AsRef<&[T]> essentially
pub trait AsSlice<T> {
    fn as_slice(&self) -> &[T];
}

impl<'a, T> AsSlice<T> for &'a [T] {
    fn as_slice(&self) -> &[T] {
        self
    }
}

impl<T> AsSlice<T> for Vec<T> {
    fn as_slice(&self) -> &[T] {
        self
    }
}

// TODO. rename to samples
// pub trait ExactSizeIteratorExt<T>: AsRef<T>
// where
//     T: ExactSizeIterator,
// {
//     fn is_empty(&self) -> bool {
//         self.as_ref().len() == 0
//     }
// }
// impl<'a, T> Samples<'a, T> for Vec<T> {
//     fn len(&self) -> usize {
//         (self as &Vec<T>).len()
//     }
// }

// impl<'a, T> Samples<'a, T> for &'a Vec<T>
// where
//     T: Copy,
// {
//     fn len(&self) -> usize {
//         (self as &'a Vec<T>).len()
//     }
// }

// impl<'a, T> Samples<'a, T> for &'a [T] {
//     fn len(&self) -> usize {
//         (self as &[T]).len()
//     }
// }
