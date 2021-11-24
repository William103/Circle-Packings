# Instructions for running code
# Compiling
You need to compile the code before you can run it. You can compile it locally, but to do that you would need to install Rust, which is probably not worth it. I would recommend `ssh`-ing into the lab computers and using them to compile it. In either case, to compile the code, simply run the command `cargo build --release` from this directory.

# Data
Next you will need to give the code the right data. Inside this directory there is a folder called `data` which has all the data. For each polyhedron, there is a corresponding file called `polyhedron.txt` (e.g. `tetrahedron.txt`) with the corresponding data. First is the Gram matrix in the format Mathematica spits out. Next is the face data. It should correspond with the Gram matrix and be in Mathematica's list format, but *it is indexed starting from 0* unlike Mathematica. If I remember correctly there is some function in some Mathematica file somewhere in this repo that has a function that finds this face list from the Gram matrix. The output of that *minus 1* is what goes here.

*Warning: most of the data files are outdated, cube is a good, trustworthy, updated example to follow.*

# Actually running
I would recommend running the command `./target/release/circle_counting --max=<max> --n=<n> --debug ./data/<polyhedron data file>` where `max` is the maximum curvature and `n` is the number of data points to collect. I'd recommend starting with `max` around `1e3` or `1e4` to see if things are working, then slowly bumping it up to around `1e6`. `n` isn't as important, but somewhere around `100` is more than enough.

# Output
The output is fairly self-explanatory. At first it will spit out a bunch of matrices. These are stuff like generators and the inital tuple it generated, which can be safely ignored unless something really crazy is going on. While computing it won't actually print anything (this might change in the future), but don't worry it should be fine. After it finishes it will print out the sample points it checked and the counts for each of them as well as some regression info, the time, and the final fractal dimension. If you want a lot less output, you can replace the `--debug` flag for a `--time` flag and it will only print out the final time and fractal dimension.
