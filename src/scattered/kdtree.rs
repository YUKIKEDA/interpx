//! A k-d tree with a search for the nearest `k` points that does not
//! allocate.
//!
//! `kiddo`, which finds the one nearest point for
//! [`NearestNDInterpolator`](crate::NearestNDInterpolator), returns several
//! points in a `Vec` of its own, new at every call. The `neighbors` mode of the
//! RBF interpolant asks for `k` points at every evaluation, and an evaluation
//! is not to allocate: this tree writes them to a buffer of the caller.
//!
//! It is slower than `kiddo` (1.2 to 1.5 times for `k` points, 1.6 to 2.1
//! times for one), which is why the nearest-neighbor interpolant stays with
//! `kiddo`.

/// Ranges of this many points or fewer are scanned, and not split.
const LEAF: usize = 16;

/// The deepest a tree gets: every split halves a range, and no range holds
/// more than `usize::MAX` points.
const DEPTH: usize = usize::BITS as usize;

/// A point found by [`Tree::nearest`]: the square of its distance from the
/// query, and the place it was given at.
///
/// The order is by distance, and among points at one distance by the place
/// given: the nearest `k` points are then one set, whatever the layout of the
/// tree.
pub(crate) type Found = (f64, u32);

/// Points in `D` dimensions, in an order that lets the nearest ones to a query
/// be found without looking at them all.
///
/// The node of a range is its middle point. It splits the range along one
/// axis, the axes taking turns from one level to the next: the points before
/// it are not above it along that axis, and those after it are not below.
#[derive(Debug, Clone)]
pub(crate) struct Tree<const D: usize> {
    /// The points, in the order of the tree.
    points: Vec<[f64; D]>,
    /// The place each of them was given at.
    given: Vec<u32>,
}

impl<const D: usize> Tree<D> {
    /// Builds the tree. `points` has no more than `u32::MAX` entries, all
    /// finite.
    pub(crate) fn new(points: &[[f64; D]]) -> Self {
        let mut given: Vec<u32> = (0..points.len() as u32).collect();
        // Ranges still to split, with the axis of their node.
        let mut ranges = vec![(0, points.len(), 0)];
        while let Some((lo, hi, axis)) = ranges.pop() {
            if hi - lo <= LEAF || D == 0 {
                continue;
            }
            let middle = (hi - lo) / 2;
            given[lo..hi].select_nth_unstable_by(middle, |a, b| {
                points[*a as usize][axis].total_cmp(&points[*b as usize][axis])
            });
            let next = (axis + 1) % D;
            ranges.push((lo, lo + middle, next));
            ranges.push((lo + middle + 1, hi, next));
        }
        Self {
            points: given.iter().map(|place| points[*place as usize]).collect(),
            given,
        }
    }

    /// Finds the `k` points nearest to `query` (all of them when there are
    /// fewer), and leaves them in `found`, in no order.
    ///
    /// `found` is cleared first. It grows to `k` entries and is not grown
    /// after that: with room for `k` this does not allocate.
    pub(crate) fn nearest(&self, query: &[f64; D], k: usize, found: &mut Vec<Found>) {
        found.clear();
        if k == 0 {
            return;
        }
        // Ranges still to look at: the range, the axis of its node, the least
        // squared distance from the query to its box, and how far the query
        // is outside that box along each axis.
        let mut ranges = [(0, 0, 0, 0.0, [0.0; D]); DEPTH];
        ranges[0] = (0, self.points.len(), 0, 0.0, [0.0; D]);
        let mut waiting = 1;
        while waiting > 0 {
            waiting -= 1;
            let (mut lo, mut hi, mut axis, least, outside) = ranges[waiting];
            // `found` is a heap with its farthest point first. A range as
            // far as that point is still looked at: it may hold a point at
            // the same distance that was given earlier.
            if found.len() == k && least > found[0].0 {
                continue;
            }
            // Down the side of the query, leaving the other side of every
            // node for later.
            while hi - lo > LEAF && D > 0 {
                let middle = lo + (hi - lo) / 2;
                self.offer(query, middle, k, found);
                let apart = query[axis] - self.points[middle][axis];
                let next = (axis + 1) % D;
                let (near, far) = if apart < 0.0 {
                    ((lo, middle), (middle + 1, hi))
                } else {
                    ((middle + 1, hi), (lo, middle))
                };
                // Along this axis the query is `apart` outside the far side,
                // in place of what it was outside the whole range.
                let mut beyond = outside;
                beyond[axis] = apart;
                let farther = least - outside[axis] * outside[axis] + apart * apart;
                ranges[waiting] = (far.0, far.1, next, farther, beyond);
                waiting += 1;
                (lo, hi, axis) = (near.0, near.1, next);
            }
            for place in lo..hi {
                self.offer(query, place, k, found);
            }
        }
    }

    /// Puts the point at `place` of the tree among the `k` found so far, when
    /// it is nearer than the farthest of them.
    fn offer(&self, query: &[f64; D], place: usize, k: usize, found: &mut Vec<Found>) {
        // From 0.0: the sum of no terms would be -0.0.
        let squared = self.points[place]
            .iter()
            .zip(query)
            .fold(0.0, |sum, (p, q)| sum + (p - q) * (p - q));
        let entry = (squared, self.given[place]);
        if found.len() < k {
            // A new leaf of the heap, moved up while it is farther than its
            // parent.
            found.push(entry);
            let mut child = found.len() - 1;
            while child > 0 {
                let parent = (child - 1) / 2;
                if found[parent] >= found[child] {
                    break;
                }
                found.swap(parent, child);
                child = parent;
            }
        } else if entry < found[0] {
            // In place of the farthest, moved down while a child is farther.
            found[0] = entry;
            let mut parent = 0;
            loop {
                let mut farthest = parent;
                for child in [2 * parent + 1, 2 * parent + 2] {
                    if child < k && found[child] > found[farthest] {
                        farthest = child;
                    }
                }
                if farthest == parent {
                    break;
                }
                found.swap(parent, farthest);
                parent = farthest;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Found, Tree};

    /// A small generator, so that the points are the same on every run.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }
    }

    /// The reference: every point, sorted by distance and then by place.
    fn scan<const D: usize>(points: &[[f64; D]], query: &[f64; D], k: usize) -> Vec<Found> {
        let mut all: Vec<Found> = points
            .iter()
            .enumerate()
            .map(|(place, p)| {
                let squared = p.iter().zip(query).map(|(p, q)| (p - q) * (p - q)).sum();
                (squared, place as u32)
            })
            .collect();
        all.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
        all.truncate(k);
        all
    }

    fn sorted(mut found: Vec<Found>) -> Vec<(u64, u32)> {
        found.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
        found
            .iter()
            .map(|(d, place)| (d.to_bits(), *place))
            .collect()
    }

    fn check<const D: usize>(points: &[[f64; D]], queries: &[[f64; D]], what: &str) {
        let tree = Tree::new(points);
        let mut found = Vec::new();
        for k in [0, 1, 2, 5, 16, 17, 40, points.len(), points.len() + 3] {
            for query in queries {
                tree.nearest(query, k, &mut found);
                let want = scan(points, query, k);
                assert_eq!(found.len(), k.min(points.len()), "{what}, k = {k}");
                assert_eq!(sorted(found.clone()), sorted(want), "{what}, k = {k}");
            }
        }
    }

    #[test]
    fn the_nearest_points_are_those_a_scan_of_all_points_finds() {
        let mut rng = Lcg(11);
        for n in [1, 2, 16, 17, 33, 100, 400] {
            let flat: Vec<[f64; 1]> = (0..n).map(|_| [rng.next()]).collect();
            let plane: Vec<[f64; 2]> = (0..n).map(|_| [rng.next(), 10.0 * rng.next()]).collect();
            let space: Vec<[f64; 4]> = (0..n)
                .map(|_| core::array::from_fn(|_| rng.next() - 0.5))
                .collect();
            // Queries among the points, at a point, and far outside.
            let mut inside = |d: usize| -> Vec<f64> { (0..d).map(|_| rng.next()).collect() };
            let q1 = [[inside(1)[0]], flat[0], [-5.0], [7.0]];
            let p = inside(2);
            let q2 = [[p[0], 10.0 * p[1]], plane[n / 2], [-3.0, 40.0]];
            let p = inside(4);
            let q4 = [
                [p[0] - 0.5, p[1] - 0.5, p[2] - 0.5, p[3] - 0.5],
                space[n - 1],
                [9.0; 4],
            ];
            check(&flat, &q1, &format!("{n} points on a line"));
            check(&plane, &q2, &format!("{n} points in a plane"));
            check(&space, &q4, &format!("{n} points in four dimensions"));
        }
    }

    #[test]
    fn points_at_one_distance_are_taken_in_the_order_given() {
        // A lattice: a query at the center of a cell has four nearest points,
        // and one at a point of the lattice has four at the next distance.
        let mut lattice = Vec::new();
        for i in 0..12 {
            for j in 0..12 {
                lattice.push([f64::from(i), f64::from(j)]);
            }
        }
        let queries = [[5.5, 5.5], [3.0, 4.0], [0.0, 0.0], [11.5, 0.5], [6.5, 6.0]];
        check(&lattice, &queries, "a lattice");
        // The same point many times: the first ones given.
        let same = vec![[1.0, 2.0, 3.0]; 50];
        let tree = Tree::new(&same);
        let mut found = Vec::new();
        tree.nearest(&[0.0, 0.0, 0.0], 7, &mut found);
        let mut places: Vec<u32> = found.iter().map(|(_, place)| *place).collect();
        places.sort_unstable();
        assert_eq!(places, [0, 1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn without_dimensions_every_point_is_at_the_distance_zero() {
        let points = vec![[0.0_f64; 0]; 40];
        let tree = Tree::new(&points);
        let mut found = Vec::new();
        tree.nearest(&[], 3, &mut found);
        let mut places: Vec<u32> = found.iter().map(|(_, place)| *place).collect();
        places.sort_unstable();
        assert_eq!(places, [0, 1, 2]);
        assert!(found.iter().all(|(d, _)| d.to_bits() == 0.0_f64.to_bits()));
    }

    #[test]
    fn a_buffer_with_room_is_not_grown() {
        let mut rng = Lcg(5);
        let points: Vec<[f64; 3]> = (0..500)
            .map(|_| core::array::from_fn(|_| rng.next()))
            .collect();
        let tree = Tree::new(&points);
        let mut found = Vec::with_capacity(30);
        let capacity = found.capacity();
        for _ in 0..50 {
            let query = core::array::from_fn(|_| rng.next());
            tree.nearest(&query, 30, &mut found);
            assert_eq!((found.len(), found.capacity()), (30, capacity));
        }
    }
}
