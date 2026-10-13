//! What the stored forms of the `serde` feature are made of, and how they are
//! read back (docs/design.ja.md section 6.5).
//!
//! Every interpolator is stored through a private form of its own: a version,
//! and the state it was built into. Reading one back runs the checks of the
//! type's constructors, so that what is loaded is what could have been built.
//!
//! The scalars, the values, and the arrays of them are written by [`Stored`]
//! and not by `Serialize`: `serde` has no array of every length, and a bound on
//! `Serialize` would change the public traits with the feature.

use serde::de::{self, Deserializer, SeqAccess, Visitor};
use serde::ser::{SerializeSeq, SerializeTuple, Serializer};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::fmt;
use std::marker::PhantomData;

/// The most elements a sequence is given room for before they are read: the
/// length a stored form claims is not believed.
const ROOM: usize = 4096;

/// A part of a stored form: a scalar, or an array or a sequence of parts.
pub trait Stored: Sized {
    /// Writes the part.
    fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error>;

    /// Reads the part.
    fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error>;

    /// A part to be written over while an array is read.
    fn blank() -> Self;
}

macro_rules! scalar {
    ($t:ty) => {
        impl Stored for $t {
            fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
                self.serialize(to)
            }

            fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
                Self::deserialize(from)
            }

            fn blank() -> Self {
                0.0
            }
        }
    };
}

scalar!(f32);
scalar!(f64);

/// A part as `serde` writes it.
struct Written<'a, X>(&'a X);

impl<X: Stored> Serialize for Written<'_, X> {
    fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
        self.0.store(to)
    }
}

/// A part as `serde` reads it.
struct Read<X>(X);

impl<'de, X: Stored> Deserialize<'de> for Read<X> {
    fn deserialize<D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
        X::load(from).map(Self)
    }
}

/// An array is a tuple: its length is the type's, and is not written.
impl<X: Stored, const N: usize> Stored for [X; N] {
    fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
        let mut tuple = to.serialize_tuple(N)?;
        for part in self {
            tuple.serialize_element(&Written(part))?;
        }
        tuple.end()
    }

    fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
        struct Parts<X, const N: usize>(PhantomData<X>);

        impl<'de, X: Stored, const N: usize> Visitor<'de> for Parts<X, N> {
            type Value = [X; N];

            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "an array of {N} elements")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                let mut parts: [X; N] = core::array::from_fn(|_| X::blank());
                for (i, part) in parts.iter_mut().enumerate() {
                    let Some(Read(read)) = seq.next_element()? else {
                        return Err(de::Error::invalid_length(i, &self));
                    };
                    *part = read;
                }
                Ok(parts)
            }
        }

        from.deserialize_tuple(N, Parts(PhantomData))
    }

    fn blank() -> Self {
        core::array::from_fn(|_| X::blank())
    }
}

/// Writes a sequence of parts.
fn store_all<X: Stored, S: Serializer>(parts: &[X], to: S) -> Result<S::Ok, S::Error> {
    let mut seq = to.serialize_seq(Some(parts.len()))?;
    for part in parts {
        seq.serialize_element(&Written(part))?;
    }
    seq.end()
}

impl<X: Stored> Stored for Vec<X> {
    fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
        store_all(self, to)
    }

    fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
        struct Parts<X>(PhantomData<X>);

        impl<'de, X: Stored> Visitor<'de> for Parts<X> {
            type Value = Vec<X>;

            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a sequence")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                let mut parts = Vec::with_capacity(seq.size_hint().unwrap_or(0).min(ROOM));
                while let Some(Read(part)) = seq.next_element()? {
                    parts.push(part);
                }
                Ok(parts)
            }
        }

        from.deserialize_seq(Parts(PhantomData))
    }

    fn blank() -> Self {
        Self::new()
    }
}

/// A sequence that is borrowed to be written and owned once read.
impl<X: Stored + Clone> Stored for Cow<'_, [X]> {
    fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
        store_all(self, to)
    }

    fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
        Vec::load(from).map(Cow::Owned)
    }

    fn blank() -> Self {
        Cow::Owned(Vec::new())
    }
}

/// An array that is borrowed to be written and owned once read.
impl<X: Stored + Clone, const N: usize> Stored for Cow<'_, [X; N]> {
    fn store<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
        (**self).store(to)
    }

    fn load<'de, D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
        <[X; N]>::load(from).map(Cow::Owned)
    }

    fn blank() -> Self {
        Cow::Owned(<[X; N]>::blank())
    }
}

/// For `#[serde(with = "crate::stored")]` on a field that is a part.
pub(crate) fn serialize<X: Stored, S: Serializer>(part: &X, to: S) -> Result<S::Ok, S::Error> {
    part.store(to)
}

/// For `#[serde(with = "crate::stored")]` on a field that is a part.
pub(crate) fn deserialize<'de, X: Stored, D: Deserializer<'de>>(from: D) -> Result<X, D::Error> {
    X::load(from)
}

/// Checks the version a stored form carries against the one this build writes.
///
/// # Errors
///
/// An error of the deserializer when they differ: a form of another version
/// is not guessed at.
pub(crate) fn version<E: de::Error>(found: u32, known: u32) -> Result<(), E> {
    if found == known {
        Ok(())
    } else {
        Err(E::custom(format_args!(
            "unknown format version {found}: this build reads version {known}"
        )))
    }
}

/// Checks the degree of the piecewise polynomial or the spline inside a stored
/// interpolator against the one its type evaluates.
///
/// # Errors
///
/// An error of the deserializer when they differ.
pub(crate) fn degree<E: de::Error>(found: usize, of_the_type: usize) -> Result<(), E> {
    if found == of_the_type {
        Ok(())
    } else {
        Err(E::custom(format_args!(
            "degree {found} where this type holds degree {of_the_type}"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::{Read, Stored, Written};

    /// Through JSON and back.
    fn round_trip<X: Stored>(part: &X) -> Result<X, serde_json::Error> {
        let text = serde_json::to_string(&Written(part))?;
        serde_json::from_str::<Read<X>>(&text).map(|read| read.0)
    }

    #[test]
    fn arrays_of_every_length_and_sequences_of_them_come_back() -> Result<(), serde_json::Error> {
        let none: [f64; 0] = [];
        assert_eq!(round_trip(&none)?, none);
        let long: [f32; 40] = core::array::from_fn(|i| i as f32 + 0.5);
        assert_eq!(round_trip(&long)?, long);
        let nested = vec![[1.5_f64, -2.0], [0.1, 4.0]];
        assert_eq!(round_trip(&nested)?, nested);
        let axes = [vec![0.0_f64, 1.0], vec![], vec![3.0]];
        assert_eq!(round_trip(&axes)?, axes);
        assert_eq!(
            serde_json::to_string(&Written(&nested))?,
            "[[1.5,-2.0],[0.1,4.0]]"
        );
        Ok(())
    }

    #[test]
    fn an_array_of_another_length_is_an_error() {
        assert!(serde_json::from_str::<Read<[f64; 3]>>("[1.0, 2.0]").is_err());
        assert!(serde_json::from_str::<Read<[f64; 1]>>("[1.0, 2.0]").is_err());
        assert!(serde_json::from_str::<Read<Vec<[f64; 2]>>>("[[1.0, 2.0], [3.0]]").is_err());
    }

    #[test]
    fn a_version_or_a_degree_that_is_not_the_known_one_is_an_error() {
        type E = serde::de::value::Error;
        assert!(super::version::<E>(1, 1).is_ok());
        assert!(super::version::<E>(2, 1).is_err());
        assert!(super::degree::<E>(3, 3).is_ok());
        assert!(super::degree::<E>(1, 3).is_err());
    }
}
