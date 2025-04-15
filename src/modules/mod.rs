pub mod match_module;
pub mod distance;
pub mod coverage;
pub mod miniset;
pub mod collection;

// 重命名match模块，因为match是Rust关键字
pub use self::match_module as match_mod; 