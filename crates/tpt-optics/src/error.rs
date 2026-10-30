//! One error type for the whole suite.
//!
//! Each crate keeps its own precise error type. [`Error`] wraps them, so
//! code that calls several crates can return a single type. Each variant
//! has a coarse [`ErrorKind`] for dispatch, and [`Error::hint`] gives a
//! suggested fix when the crate can offer one.

use alloc::string::String;
use core::fmt;

/// Coarse category of an error.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    /// Text could not be parsed or decoded (file formats).
    Parse,
    /// An argument or a prescription is invalid.
    Input,
    /// A glass or material name could not be resolved.
    Glass,
    /// A numerical method refused its input (size, conditioning, evanescence).
    Numerics,
}

/// Any error from the suite.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// A dispersion table is malformed.
    #[cfg(feature = "core")]
    Table(tpt_optics_core::dispersion::TableError),
    /// A prescription could not be built.
    #[cfg(feature = "geometric")]
    Build(tpt_optics_geometric::builder::BuildError),
    /// A Zemax `.zmx` file is malformed or unsupported.
    #[cfg(feature = "io")]
    Zmx(tpt_optics_io::zmx::ZmxError),
    /// A CODE V `.seq` file is malformed or unsupported.
    #[cfg(feature = "io")]
    Seq(tpt_optics_io::seq::SeqError),
    /// An AGF glass catalogue is malformed.
    #[cfg(feature = "io")]
    Agf(tpt_optics_io::glass::AgfError),
    /// A text file is not valid UTF-8 or UTF-16.
    #[cfg(feature = "io")]
    Decode(tpt_optics_io::text::DecodeError),
    /// A supercell photonic-crystal problem is invalid.
    #[cfg(feature = "waveguide")]
    Supercell(tpt_optics_waveguide::supercell::SupercellError),
    /// A coupled-mode helper refused its input.
    #[cfg(feature = "waveguide")]
    Coupled(tpt_optics_waveguide::coupled::CoupledError),
    /// The vector mode solver refused or failed.
    #[cfg(feature = "waveguide")]
    VectorFem(tpt_optics_waveguide::vector_fem::VectorFemError),
    /// A parameter sweep is malformed.
    #[cfg(feature = "design")]
    Sweep(tpt_optics_design::sweep::SweepError),
    /// A named-parameter operation failed.
    #[cfg(feature = "design")]
    Param(tpt_optics_design::params::ParamError),
}

impl Error {
    /// Coarse category of this error.
    #[must_use]
    pub fn kind(&self) -> ErrorKind {
        match self {
            #[cfg(feature = "core")]
            Error::Table(_) => ErrorKind::Input,
            #[cfg(feature = "geometric")]
            Error::Build(tpt_optics_geometric::builder::BuildError::UnknownGlass { .. }) => {
                ErrorKind::Glass
            }
            #[cfg(feature = "geometric")]
            Error::Build(_) => ErrorKind::Input,
            #[cfg(feature = "io")]
            Error::Zmx(_) | Error::Seq(_) | Error::Agf(_) | Error::Decode(_) => ErrorKind::Parse,
            #[cfg(feature = "waveguide")]
            Error::Supercell(_) | Error::Coupled(_) | Error::VectorFem(_) => ErrorKind::Numerics,
            #[cfg(feature = "design")]
            Error::Sweep(_) | Error::Param(_) => ErrorKind::Input,
            #[allow(unreachable_patterns)]
            _ => ErrorKind::Input,
        }
    }

    /// A suggested fix, when the underlying error offers one. Currently this
    /// is the closest glass name for an unknown glass.
    #[must_use]
    pub fn hint(&self) -> Option<String> {
        match self {
            #[cfg(feature = "geometric")]
            Error::Build(tpt_optics_geometric::builder::BuildError::UnknownGlass {
                suggestion: Some(near),
                ..
            }) => Some(alloc::format!("did you mean '{near}'?")),
            _ => None,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            #[cfg(feature = "core")]
            Error::Table(e) => write!(f, "{e}"),
            #[cfg(feature = "geometric")]
            Error::Build(e) => write!(f, "{e}"),
            #[cfg(feature = "io")]
            Error::Zmx(e) => write!(f, "{e}"),
            #[cfg(feature = "io")]
            Error::Seq(e) => write!(f, "{e}"),
            #[cfg(feature = "io")]
            Error::Agf(e) => write!(f, "{e}"),
            #[cfg(feature = "io")]
            Error::Decode(e) => write!(f, "{e}"),
            #[cfg(feature = "waveguide")]
            Error::Supercell(e) => write!(f, "{e}"),
            #[cfg(feature = "waveguide")]
            Error::Coupled(e) => write!(f, "{e}"),
            #[cfg(feature = "waveguide")]
            Error::VectorFem(e) => write!(f, "{e}"),
            #[cfg(feature = "design")]
            Error::Sweep(e) => write!(f, "{e}"),
            #[cfg(feature = "design")]
            Error::Param(e) => write!(f, "{e}"),
            #[allow(unreachable_patterns)]
            _ => write!(f, "tpt-optics error"),
        }
    }
}

impl core::error::Error for Error {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            #[cfg(feature = "core")]
            Error::Table(e) => Some(e),
            #[cfg(feature = "geometric")]
            Error::Build(e) => Some(e),
            #[cfg(feature = "io")]
            Error::Zmx(e) => Some(e),
            #[cfg(feature = "io")]
            Error::Seq(e) => Some(e),
            #[cfg(feature = "io")]
            Error::Agf(e) => Some(e),
            #[cfg(feature = "io")]
            Error::Decode(e) => Some(e),
            #[cfg(feature = "waveguide")]
            Error::Supercell(e) => Some(e),
            #[cfg(feature = "waveguide")]
            Error::Coupled(e) => Some(e),
            #[cfg(feature = "waveguide")]
            Error::VectorFem(e) => Some(e),
            #[cfg(feature = "design")]
            Error::Sweep(e) => Some(e),
            #[cfg(feature = "design")]
            Error::Param(e) => Some(e),
            #[allow(unreachable_patterns)]
            _ => None,
        }
    }
}

/// Result alias for the suite.
pub type Result<T> = core::result::Result<T, Error>;

#[cfg(feature = "core")]
impl From<tpt_optics_core::dispersion::TableError> for Error {
    fn from(e: tpt_optics_core::dispersion::TableError) -> Self {
        Error::Table(e)
    }
}

#[cfg(feature = "geometric")]
impl From<tpt_optics_geometric::builder::BuildError> for Error {
    fn from(e: tpt_optics_geometric::builder::BuildError) -> Self {
        Error::Build(e)
    }
}

#[cfg(feature = "io")]
impl From<tpt_optics_io::zmx::ZmxError> for Error {
    fn from(e: tpt_optics_io::zmx::ZmxError) -> Self {
        Error::Zmx(e)
    }
}

#[cfg(feature = "io")]
impl From<tpt_optics_io::seq::SeqError> for Error {
    fn from(e: tpt_optics_io::seq::SeqError) -> Self {
        Error::Seq(e)
    }
}

#[cfg(feature = "io")]
impl From<tpt_optics_io::glass::AgfError> for Error {
    fn from(e: tpt_optics_io::glass::AgfError) -> Self {
        Error::Agf(e)
    }
}

#[cfg(feature = "io")]
impl From<tpt_optics_io::text::DecodeError> for Error {
    fn from(e: tpt_optics_io::text::DecodeError) -> Self {
        Error::Decode(e)
    }
}

#[cfg(feature = "waveguide")]
impl From<tpt_optics_waveguide::supercell::SupercellError> for Error {
    fn from(e: tpt_optics_waveguide::supercell::SupercellError) -> Self {
        Error::Supercell(e)
    }
}

#[cfg(feature = "waveguide")]
impl From<tpt_optics_waveguide::coupled::CoupledError> for Error {
    fn from(e: tpt_optics_waveguide::coupled::CoupledError) -> Self {
        Error::Coupled(e)
    }
}

#[cfg(feature = "waveguide")]
impl From<tpt_optics_waveguide::vector_fem::VectorFemError> for Error {
    fn from(e: tpt_optics_waveguide::vector_fem::VectorFemError) -> Self {
        Error::VectorFem(e)
    }
}

#[cfg(feature = "design")]
impl From<tpt_optics_design::sweep::SweepError> for Error {
    fn from(e: tpt_optics_design::sweep::SweepError) -> Self {
        Error::Sweep(e)
    }
}

#[cfg(feature = "design")]
impl From<tpt_optics_design::params::ParamError> for Error {
    fn from(e: tpt_optics_design::params::ParamError) -> Self {
        Error::Param(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "geometric")]
    #[test]
    fn unknown_glass_has_kind_glass_and_a_hint() {
        use tpt_optics_geometric::builder::BuildError;
        let e: Error = BuildError::UnknownGlass {
            name: String::from("bk8"),
            suggestion: Some(String::from("bk7")),
        }
        .into();
        assert_eq!(e.kind(), ErrorKind::Glass);
        assert_eq!(e.hint().as_deref(), Some("did you mean 'bk7'?"));
        assert!(alloc::format!("{e}").contains("bk8"));
    }

    #[cfg(feature = "waveguide")]
    #[test]
    fn numerical_refusals_have_kind_numerics() {
        use tpt_optics_waveguide::vector_fem::VectorFemError;
        let e: Error = VectorFemError::Singular.into();
        assert_eq!(e.kind(), ErrorKind::Numerics);
        assert_eq!(e.hint(), None);
    }

    #[cfg(feature = "io")]
    #[test]
    fn parse_errors_have_kind_parse_and_expose_their_source() {
        use core::error::Error as _;
        let e: Error = tpt_optics_io::text::DecodeError::InvalidUtf8.into();
        assert_eq!(e.kind(), ErrorKind::Parse);
        assert!(e.source().is_some());
    }
}
