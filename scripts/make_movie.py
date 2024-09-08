'''
Script to create movie from frames in path_to_frames and save it in output_path.
Assumes frames are saved as .png, and creates an mp4 file at 20fps
Usage:
    make_movie.py <path_to_frames>
    make_movie.py <path_to_frames> [--name=<name>]
    make_movie.py <path_to_frames> [--output_path=<dir>] [--name=<name>]

Options:
    --output_path=<dir>  Output Directory [default: ./]
    --name=<name>  Movie Name [default: movie.mp4]

'''
import imageio
from docopt import docopt
import pathlib
import glob

args = docopt(__doc__)

path_to_frames = pathlib.Path(args['<path_to_frames>']).absolute()
output_path = pathlib.Path(args['--output_path']).absolute()


# create output directory if it doesn't yet exist
if not output_path.exists():
    output_path.mkdir()

# create writer to make movie
fps = 20
output = str(output_path.joinpath(args['--name']))
w = imageio.get_writer(output, fps= fps)

# iterate over frames in directory, in alphanumeric order
frames = glob.glob(str(path_to_frames.joinpath('./*.png')))
for frame in sorted(frames):
    w.append_data(imageio.imread(str(frame)))
w.close()
