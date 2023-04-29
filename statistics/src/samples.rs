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
