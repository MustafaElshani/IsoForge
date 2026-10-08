// Core error types for IsoForge

use thiserror::Error;

#[derive(Error, Debug)]
pub enum IsoForgeError {
    #[error("BAM/SAM I/O error: {0}")]
    BamError(String),

    #[error("Parquet I/O error: {0}")]
    ParquetError(String),

    #[error("Invalid strand: {0}")]
    InvalidStrand(String),

    #[error("Invalid position type: {0}. Expected 'start' or 'end'")]
    InvalidPosType(String),

    #[error("Database integrity error: {0}")]
    IntegrityError(String),

    #[error("Configuration error: {0}")]
    ConfigError(String),

    #[error("Read QC failure: {0}")]
    QcError(String),

    #[error("No matching gene found for transcript")]
    NoGeneMatch,

    #[error("Arrow error: {0}")]
    ArrowError(#[from] arrow::error::ArrowError),

    #[error("Parquet IO error: {0}")]
    ParquetIoError(#[from] parquet::errors::ParquetError),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("HTS error: {0}")]
    HtsError(String),

    #[error("{0}")]
    Other(String),
}

pub type IsoForgeResult<T> = Result<T, IsoForgeError>;
