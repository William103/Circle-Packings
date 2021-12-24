use structopt::StructOpt;

/// General tool for working with polyhedral circle packings. Can be used to create high detail svg
/// images of packings, compute fractal dimensions using the circle counting method, or do
/// computations on gram matrices, generators, etc.
#[derive(StructOpt, Debug)]
#[structopt(name = "Polyhedral Circle Packings")]
pub enum Polyhedral {
    /// Compute fractal dimension of a polyhedral circle packing
    Dimension {
        /// File containing the gram matrix as well as information about the faces
        #[structopt(name = "data file")]
        data_file: String,

        /// Activate debug mode
        #[structopt(short, long)]
        debug: bool,

        /// Number of sample points in linear regression
        #[structopt(short, long, default_value = "50")]
        n: usize,

        /// Maximum curvature of circles
        #[structopt(short, long, default_value = "1000000")]
        max: f64,

        /// Whether or not to time execution
        #[structopt(short, long)]
        time: bool,

        /// Recursion cap (0 means no cap)
        #[structopt(short, long, default_value = "0")]
        recursion_depth: usize,
    },

    /// Create an image of a polyhedral circle packing
    Picture {
        /// File containing the gram matrix as well as information about the faces
        #[structopt(name = "data file")]
        data_file: String,

        /// Activate debug mode
        #[structopt(short, long)]
        debug: bool,

        /// Maximum curvature of circles
        #[structopt(short, long, default_value = "1000")]
        max: f64,

        /// Whether or not to time execution (automatically enabled by --debug)
        #[structopt(short, long)]
        time: bool,

        /// Recursion cap (0 means no cap)
        #[structopt(short, long, default_value = "0")]
        recursion_depth: usize,

        /// Name of the output file
        #[structopt(short, long, default_value = "output.svg")]
        output_file: String,

        /// Width in pixels
        #[structopt(short, long, default_value = "1000")]
        width: usize,

        /// Height in pixels
        #[structopt(short, long, default_value = "1000")]
        height: usize,
    },
}
