
// pub trait IsizeIndex<T> {
//     fn get_el(&self, i : isize) -> &T;
//     fn get_mut_el(&mut self, i : isize) -> &mut T;
// }
// impl<T> IsizeIndex<T> for Vec<T> where Vec{
//     fn get_el(&self, i : isize) -> &T {
//         &self[i as usize]
//     }
//     fn get_mut_el(&mut self, i: isize) -> &mut T {
//         &mut self[i as usize]
//     }
// }

pub trait IsizeIndexMatrix<T> {
    fn at(&self, i : (isize,isize)) -> &T;
    fn at_mut(&mut self, i : (isize,isize)) -> &mut T;
}
impl<T> IsizeIndexMatrix<T> for Vec<Vec<T>> {
    fn at(&self, i: (isize, isize)) -> &T {
        &self[i.0 as usize][i.1 as usize]
    }
    fn at_mut(&mut self, i: (isize, isize)) -> &mut T {
        &mut self[i.0 as usize][i.1 as usize]
    }
}