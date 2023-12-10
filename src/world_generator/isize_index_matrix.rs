pub trait IsizeIndexMatrix<T> {
    fn at(&self, i: (isize, isize)) -> &T;
    fn at_mut(&mut self, i: (isize, isize)) -> &mut T;
    fn at_checked(&self, i: (isize, isize)) -> Option<&T>;
    fn at_mut_checked(&mut self, i: (isize, isize)) -> Option<&mut T>;
}
impl<T> IsizeIndexMatrix<T> for Vec<Vec<T>> {
    fn at(&self, i: (isize, isize)) -> &T {
        &self[i.0 as usize][i.1 as usize]
    }
    fn at_mut(&mut self, i: (isize, isize)) -> &mut T {
        &mut self[i.0 as usize][i.1 as usize]
    }

    fn at_checked(&self, i: (isize, isize)) -> Option<&T> {
        if i.0 < 0 || i.1 < 0 || i.0 >= self.len() as isize || i.1 >= self[i.0 as usize].len() as isize {
            None
        } else {
            Some(self.at(i))
        }
    }

    fn at_mut_checked(&mut self, i: (isize, isize)) -> Option<&mut T> {
        if i.0 < 0 || i.1 < 0 || i.0 >= self.len() as isize || i.1 >= self[i.0 as usize].len() as isize {
            None
        } else {
            Some(self.at_mut(i))
        }
    }
}
