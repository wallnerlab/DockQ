//! File loading and compression dispatch.

use std::fs;
use std::io::Read;
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::{DockError, Result};
use crate::io::{mmcif, pdb};
use crate::types::{LoadOptions, Structure};

/// Load a PDB or mmCIF structure from a plain-text or gzip-compressed file.
///
/// Parsing is attempted as PDB first and mmCIF second, matching the baseline
/// loader's format fallback. Input and decompression errors are reported before
/// parser dispatch so they retain their original cause and path context.
pub fn load_structure(path: impl AsRef<Path>, options: &LoadOptions) -> Result<Structure> {
    let path = path.as_ref();
    let text = read_structure_text(path)?;

    let mut structure = match pdb::parse_pdb(&text, options) {
        Ok(structure) => structure,
        Err(pdb_error) => match mmcif::parse_mmcif(&text, options) {
            Ok(structure) => structure,
            Err(mmcif_error) => {
                return Err(DockError(format!(
                    "could not parse structure {} as PDB ({}) or mmCIF ({})",
                    path.display(),
                    pdb_error,
                    mmcif_error
                )));
            }
        },
    };

    structure.id = path.to_string_lossy().into_owned();
    Ok(structure)
}

fn read_structure_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|error| {
        DockError(format!(
            "failed to read structure {}: {error}",
            path.display()
        ))
    })?;
    decode_structure_text(path, &bytes)
}

fn decode_structure_text(path: &Path, bytes: &[u8]) -> Result<String> {
    let is_gzip = path.as_os_str().to_string_lossy().ends_with(".gz");
    let decoded = if is_gzip {
        let mut decoder = GzDecoder::new(bytes);
        let mut text = String::new();
        decoder.read_to_string(&mut text).map_err(|error| {
            DockError(format!(
                "failed to decompress structure {}: {error}",
                path.display()
            ))
        })?;
        text
    } else {
        String::from_utf8(bytes.to_vec()).map_err(|error| {
            DockError(format!(
                "structure {} is not valid UTF-8: {error}",
                path.display()
            ))
        })?
    };
    Ok(decoded)
}

#[cfg(test)]
mod tests {
    use super::decode_structure_text;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use std::path::Path;

    #[test]
    fn decodes_plain_text() {
        assert_eq!(
            decode_structure_text(Path::new("model.pdb"), b"ATOM\n").unwrap(),
            "ATOM\n"
        );
    }

    #[test]
    fn decodes_gzip_text() {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"data_example\n").unwrap();
        let compressed = encoder.finish().unwrap();
        assert_eq!(
            decode_structure_text(Path::new("model.cif.gz"), &compressed).unwrap(),
            "data_example\n"
        );
    }

    #[test]
    fn reports_corrupt_gzip_with_path() {
        let error = decode_structure_text(Path::new("broken.cif.gz"), b"not gzip").unwrap_err();
        assert!(error.to_string().contains("broken.cif.gz"));
        assert!(error.to_string().contains("decompress"));
    }

    #[test]
    fn reports_invalid_plain_utf8_with_path() {
        let error = decode_structure_text(Path::new("broken.pdb"), &[0xff]).unwrap_err();
        assert!(error.to_string().contains("broken.pdb"));
        assert!(error.to_string().contains("UTF-8"));
    }
}
