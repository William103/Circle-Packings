# Instructions for running code
# Compiling
You need to compile the code before you can run it. You can compile it locally, but to do that you would need to install Rust, which is probably not worth it. I would recommend `ssh`-ing into the lab computers and using them to compile it. In either case, to compile the code, simply run the command `cargo build --release` from this directory. If you want to run it locally, or can't `ssh` into a lab computer, let me know what platform you're on (Windows or Mac) and I'll see if I can compile it for that platform and send the binary to you.

# Data
Next you will need to give the code the right data. Inside this directory there is a folder called `data` which has all the data. For each polyhedron, there is a corresponding file called `polyhedron.txt` (e.g. `tetrahedron.txt`) with the corresponding data. First is the Gram matrix in the format Mathematica spits out. Next is the face data. It should correspond with the Gram matrix and be in Mathematica's list format, but *it is indexed starting from 0* unlike Mathematica. The file `fractal_dimension.nb` has some functions that make this convenient. You can simply type up the gram matrix, then run `FindFace[graphFromG[G]]-1` and that will spit out exactly what you want. Then put this somewhere in the data folder where it makes sense.

# Actually running
## Fractal dimension
I wrote a little script called `driver.sh` that runs this code on everything in the data folder that isn't already computed and stores the result in a folder called `output`. You should just be able to run `./driver.sh <max_curvature>` (I'd recommend `./driver.sh 1e3` to see if it works, then `./driver.sh 1e6` and/or `./driver.sh 1e7`). If you want to run it manually, you can run `./target/release/circle_counting dimension --help` and it will tell you what the options and their defaults are. You can also tweak the options in `driver.sh`, but that can have some unexpected effects.

## Computing Generators/Root Tuple 
This program can also compute generators (both geometric and algebraic, but mathematica can probably get more accurate results for the algebraic generators) and a bounded root tuple. The command for that is `./target/release/circle_counting generators <data_file>`. If you want the output formatted like python (e.g. for use with Emmi's McMullen code), add the `--format=python` flag. The formats `mathematica`, `c`, and `human-readable` are also handled. I would also recommend adding `> output.txt` or something like that to the end to save the results in a text file.

## Picture
If you want to generate a picture, you can run `./target/release/circle_counting picture --help` and it will tell you how to run that. I would recommend just running `./target/release/circle_counting picture <data_file>`, as the defaults are pretty good.
