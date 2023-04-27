// TODO: use ExactSizeIterator

pub trait Samples<'a, T: 'a>: std::iter::IntoIterator<Item = &'a T> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn len(&self) -> usize;
}

// impl<'a, T> Samples<'a, T> for Vec<T> {
//     fn len(&self) -> usize {
//         (self as &Vec<T>).len()
//     }
// }

impl<'a, T> Samples<'a, T> for &'a Vec<T>
where
    T: Copy,
{
    fn len(&self) -> usize {
        (self as &'a Vec<T>).len()
    }
}

impl<'a, T> Samples<'a, T> for &'a [T] {
    fn len(&self) -> usize {
        (self as &[T]).len()
    }
}
