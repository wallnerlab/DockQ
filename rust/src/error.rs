use std::fmt;

/// Errors are returned by library operations; only the executable prints/exits.
#[derive(Debug, Clone, PartialEq)]
pub struct DockError(pub String);
impl fmt::Display for DockError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
impl std::error::Error for DockError {}
impl From<std::io::Error> for DockError {
    fn from(e: std::io::Error) -> Self {
        Self(e.to_string())
    }
}
pub type Result<T> = std::result::Result<T, DockError>;
